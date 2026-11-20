# This is a function to binning peak count, and clustering.
suppressPackageStartupMessages({
  library(Matrix)
  library(Signac)
  library(Seurat)
  library(tidyr)
  library(ggplot2)
  library(futile.logger)
})
# source("funs_find_subclones.R") 
# source("mydataProcess.R") 
# source("figurePlot.R")


#' @title binClust()
#' @description Process the input matrix into a bin count matrix if required, then calculate 
#' the ratio relative to the reference, and perform PCA dimensional reduction and clustering 
#' based on the ratio.
#' @param cellMeta data.frame of cell annotation, group information on the first column, cell IDs on rownames
#' @param ref_group_names the group name as reference
#' @param chrs analysis for specific chromosomes
#' @param cytoBandFile the reference annotation of chromosome location, BED format with columns of c("chrom","chromStart","chromEnd","name","gieStain")
#' @param seu_resol A non-negative numeric resolution, or a vector of candidate resolutions.
#' @param silhouette_npcs Maximum number of PCs used to evaluate candidate resolutions.
#' @param silhouette_max_cells Maximum number of cells for silhouette evaluation (Inf for all cells).
#' @param silhouette_seed Random seed for sampling and candidate clustering.
#' @export
binClust <- function(inputMat,cellMeta,ref_group_names=NULL,
                     chrs = NULL,
                     delimit = "-",
                     nFeature.min = 50,
                     #cellProp_cutoff = 0.5,
                     Center = TRUE,
                     onRatio=TRUE,
                     outdir=".",
                     window = 5,
                     seu_resol = 0.8,
                     outputFigures = TRUE,
                     plt4peak = FALSE,
                     OutMat = TRUE,
                     Figure_name = NULL,
                     noRefShow = FALSE,
                     count_lim = NULL,
                     Figure.width=10,
                     Figure.height=6,
                     cytoBand = NULL,
                     global_norm = FALSE,
                     silhouette_npcs = 30L,
                     silhouette_max_cells = 5000L,
                     silhouette_seed = 123L
                     ) {
  suppressMessages({
    library(ggplot2)
    library(dplyr)
    library(stringr)
  })
  appender_func <- flog.logger()$appender
  log_path <- environment(appender_func)$file
  if(!is.null(log_path)){
    if(file.exists(log_path)){flog.appender(appender.file(log_path))}else{flog.appender(appender.console())}
  }else{
    flog.appender(appender.console())
  }
  flog.info("\n")
  if(outdir != "." && !file.exists(outdir)){
    flog.info(paste0("Creating output path ", outdir))
    dir.create(outdir,recursive=T)
  }
  if(!is.null(cytoBand)){
    if (Reduce("|", is(cytoBand) == "character")) {
      if(substr(cytoBand, nchar(cytoBand)-2, nchar(cytoBand)) %in% c("txt","tsv","csv")){
        cytoBand <- read.table(cytoBand, sep="\t", header=F, check.names=FALSE)
        # cytoBand <- read.table(connection <- gzfile(cytoBand, 'rt'), sep="\t", header=F, check.names=FALSE)
        # close(connection)
        cytoBand <- as.data.frame(cytoBand)
      }
    }else if(Reduce("|", is(cytoBand) %in% c("data.frame"))){
    }else{
      message("CytoBand isn't recognized as data.frame, or filename")
    }
    #check the columns of cytoBand
    if(ncol(cytoBand)==5){
      colnames(cytoBand)=c("chrom","chromStart","chromEnd","name","gieStain")
    }else{
      stop("Error, the columns of cytoBand do NOT match the columns of the cytoBand table in the database underlying the UCSC Genome Browser")
    }
  }
  
  inputMat <- as.matrix(inputMat)
  # extract tumor cells
  if(!is.null(ref_group_names)){
    cells_t <- rownames(cellMeta)[!cellMeta[,1] %in% ref_group_names]
    cells_ref <- rownames(cellMeta)[cellMeta[,1] %in% ref_group_names]
  }else{
    cells_t <- colnames(inputMat)
    cells_ref <- cells_t
  }
  
  
  if(onRatio){
    mt_ana <- inputMat
  }else{
    mt_ana <- inputMat[,cells_t,drop=F]
  }

  if(window>1){
    flog.info(sprintf("Binning count matrix by averaging the count of %s peaks.", window))
    #run for each arm region
    bins_input <- rownames(mt_ana)
    bins_input <- gsub("_|-|:","_",bins_input)
    bins_chr <- sapply(strsplit(bins_input, "_|-|:"), "[", 1)
    bins_start <- sapply(strsplit(bins_input, "_|-|:"), "[", 2)
    bins_end <- sapply(strsplit(bins_input, "_|-|:"), "[", 3)
    
    bins_bed <- data.frame(chr=bins_chr,bins_start,bins_end)
    bins_anno <- align_Grange2bin(bed.query=bins_bed,bed.subject=cytoBand)
    bins_anno$arm <- paste(bins_anno$chr,str_sub(bins_anno$name,1,1),sep="_")
    bins_anno <- bins_anno[match(bins_input,bins_anno$binID),]
    mt_ana <- as.data.frame(mt_ana)
    mt_ana$region <- bins_anno$arm
    mtls <- split(mt_ana,mt_ana$region)
    mtls <- mtls[unique(mt_ana$region)]
    mtx_bin <- do.call(rbind,lapply(mtls,function(x,window){
      x <- x[-ncol(x)]
      x <- as.matrix(x)
      mtx_bin_x <- suppressMessages(smooth_by_chromosome(x,wl=window,method="winMean",win_step=window,adjust_zero =TRUE))
      return(mtx_bin_x)
    },window))
    #mtx_bin <- suppressMessages(smooth_by_chromosome(mt_ana,wl=window,method="winMean",win_step=window,adjust_zero =TRUE))
  }else{
    mtx_bin <- mt_ana
  }

  flog.info(sprintf("Number of bins: %s", dim(mtx_bin)[1]))
  flog.info(sprintf("Number of cells: %s", dim(mtx_bin)[2]))

  #normalization on tumor cells
  flog.info("Normalize the bin count matrix if perform PCA on ratio (onRatio=TRUE).")
  
  if(global_norm){
    countsMnormal = t(t(mtx_bin)*mean(colSums(mtx_bin))/colSums(mtx_bin))
  }else if(!is.null(ref_group_names)){
    mtx_bin.t <- mtx_bin[,cells_t,drop=F]
    countsMnormal.t = t(t(mtx_bin.t)*mean(colSums(mtx_bin.t))/colSums(mtx_bin.t))
    mtx_bin.n <- mtx_bin[,cells_ref,drop=F]
    countsMnormal.n = t(t(mtx_bin.n)*mean(colSums(mtx_bin.n))/colSums(mtx_bin.n))
    countsMnormal <- cbind(countsMnormal.t,countsMnormal.n)
    
  }else{
    countsMnormal = t(t(mtx_bin)*mean(colSums(mtx_bin))/colSums(mtx_bin))
  }
  
  if(!is.null(chrs)){
    flog.info(paste0("Extract peaks from ",paste(chrs,collapse = ","),"."))
    countsMnormal <- countsMnormal[grepl(paste0("^",chrs,delimit,collapse = "|"),rownames(countsMnormal)),,drop=F]
  
    }
  
  if(onRatio){
    ##cells from ref_group_names are combined into bulk, and the rest are kept at the single cell level to calculate the ratio

    cellMeta_new <- data.frame(row.names = rownames(cellMeta),group=rownames(cellMeta))
    cellMeta_new[cells_ref,1] <- "reference"
    cell_info2 <- data.frame(rownames(cellMeta_new),cellMeta_new)
    if(!is.null(ref_group_names)){
      cell_info2_r <- cell_info2
    }else{
      cell_info2_r <- data.frame(rownames(cellMeta_new),rownames(cellMeta))
    }
    flog.info("Calculate single-cell level bin-ratio relative to reference of 'ref_group_names'.") 
    mt_b <- pseudo_bulk_v2(countsMnormal,group_ano = cell_info2,method ="mean",adjust_zero =TRUE)
    ref.value <- mt_b[,"reference",drop=F]
    
    ratio_res <- CNratio_peak(countsMnormal,cell_ano=cell_info2_r,ref.value=ref.value,colname_ref="reference",Norm_by_ref=FALSE) 
    is.na(ratio_res) <- sapply(ratio_res, is.infinite)
    is.na(ratio_res) <- sapply(ratio_res, is.nan)
    ratio_res <- na.omit(ratio_res)
    clt_matrix <- ratio_res
  }else{
    clt_matrix <- countsMnormal
  }
  
  
  # #cell scoring on Zero count in chr
  # cellScore <- cell_score(data = clt_matrix,prop_cutoff=cellProp_cutoff )
  # cells_used <- cellScore$cell[cellScore$label %in% "Low"]
  # cells_filered <- cellScore$cell[!cellScore$label %in% "Low"]
  # clt_matrix_filt <- clt_matrix[,cells_used,drop=F]

  # Initialize selection outputs, also for the low-feature fallback.
  if (!is.numeric(seu_resol) || length(seu_resol) < 1L ||
      any(!is.finite(seu_resol)) || any(seu_resol < 0)) {
    stop("seu_resol must be a non-empty, finite, non-negative numeric vector.")
  }
  cell_info_ls <- NULL
  best_resolution <- NA_real_
  resolution_eval <- NULL

  if(nrow(clt_matrix)>=nFeature.min){
    if(nrow(clt_matrix)>1e4){
       #run for each arm region
      bins <- rownames(clt_matrix)
      bins <- gsub("_|-|:","_",bins)
      bins_chr <- sapply(strsplit(bins, "_|-|:"), "[", 1)
      bins_start <- sapply(strsplit(bins, "_|-|:"), "[", 2)
      bins_end <- sapply(strsplit(bins, "_|-|:"), "[", 3)
      
      bins_bed <- data.frame(chr=bins_chr,bins_start,bins_end)
      bins_anno <- align_Grange2bin(bed.query=bins_bed,bed.subject=cytoBand)
      bins_anno$arm <- paste(bins_anno$chr,str_sub(bins_anno$name,1,1),sep="_")
      bins_anno <- bins_anno[match(bins,bins_anno$binID),]
      mt_ana <- as.data.frame(clt_matrix)
      mt_ana$region <- bins_anno$arm
      mtls <- split(mt_ana,mt_ana$region)
      mtls <- mtls[unique(mt_ana$region)]
      clt_matrix <- do.call(rbind,lapply(mtls,function(x){
        x <- x[-ncol(x)]
        x <- as.matrix(x)
        mtx_bin_x <- suppressMessages(smooth_by_chromosome(x,wl=5,method="winMean",win_step=5,adjust_zero =TRUE))
        return(mtx_bin_x)
      }))
    }

    if (length(seu_resol) == 1L) {
      cell_info_ls <- suppressMessages(find_subclust(clt_matrix,resolution = seu_resol,sep_by = c("-", "-")))      
      best_resolution <- seu_resol
    }else{
      # Compare all candidate resolutions in the same PCA space.
      if (length(silhouette_npcs) != 1L ||
          !is.finite(silhouette_npcs) || silhouette_npcs < 1 ||
          silhouette_npcs != floor(silhouette_npcs)) {
        stop("silhouette_npcs must be a positive integer.")
      }
      if (length(silhouette_max_cells) != 1L ||
          is.na(silhouette_max_cells) || silhouette_max_cells < 3 ||
          (is.finite(silhouette_max_cells) &&
           silhouette_max_cells != floor(silhouette_max_cells))) {
        stop("silhouette_max_cells must be an integer >= 3 or Inf.")
      }
      resolutions <- sort(unique(seu_resol))
      n_cells <- ncol(clt_matrix)
      cells <- colnames(clt_matrix)

      set.seed(silhouette_seed)
      idx <- if (n_cells > silhouette_max_cells) {
        sort(sample.int(n_cells, silhouette_max_cells))
      } else {
        seq_len(n_cells)
      }
      d <- NULL
      if (length(idx) >= 3L && nrow(clt_matrix) >= 1L) {
        x <- t(clt_matrix)
        n_pcs <- min(silhouette_npcs, nrow(x) - 1L, ncol(x))
        # Use truncated SVD if available; otherwise use base R PCA.
        if (requireNamespace("irlba", quietly = TRUE) &&
            n_pcs < min(dim(x)) && n_pcs < 0.5 * min(dim(x))) {
          pcs <- tryCatch(
            irlba::prcomp_irlba(x, n = n_pcs,
                               center = TRUE, scale. = FALSE)$x,
            error = function(e) stats::prcomp(x, rank. = n_pcs,
                                             center = TRUE, scale. = FALSE)$x
          )
        } else {
          pcs <- stats::prcomp(x, rank. = n_pcs,
                                center = TRUE, scale. = FALSE)$x
        }
        d <- stats::dist(pcs[idx, , drop = FALSE])
        rm(x, pcs)
        if (!any(d > 0)) d <- NULL
      }

      resolution_eval <- data.frame(
        resolution = resolutions,
        n_clusters = NA_integer_,
        n_eval = length(idx),
        n_eval_clusters = NA_integer_,
        silhouette = NA_real_,
        status = rep("not_evaluated", length(resolutions)),
        stringsAsFactors = FALSE
      )
      best_score <- -Inf
      fallback_result <- NULL
      fallback_resolution <- NA_real_

      for (i in seq_along(resolutions)) {
        res <- resolutions[i]
        set.seed(silhouette_seed)
        tmp <- tryCatch(
          suppressMessages(find_subclust(clt_matrix, resolution = res, sep_by = c("-", "-"))),
          error = function(e) e
        )
        if (inherits(tmp, "error")) {
          resolution_eval$status[i] <- paste("clustering_error:", conditionMessage(tmp))
          flog.info(sprintf("Resolution %g failed: %s", res, conditionMessage(tmp)))
          next
        }

        meta <- tmp$meta_cluster
        if (is.null(meta) || !"subCluster" %in% colnames(meta) ||
            is.null(rownames(meta)) || anyDuplicated(rownames(meta))) {
          resolution_eval$status[i] <- "invalid_meta_cluster"
          next
        }
        # Match by cell ID, never by metadata row order.
        labs_all <- as.character(meta$subCluster[match(cells, rownames(meta))])
        valid <- !is.na(labs_all) & nzchar(labs_all) & labs_all != "reference"
        k <- length(unique(labs_all[valid]))
        resolution_eval$n_clusters[i] <- k
        if (k == 0L) {
          resolution_eval$status[i] <- "no_cluster_labels"
          next
        }
        if (is.null(fallback_result)) {
          fallback_result <- tmp
          fallback_resolution <- res
        }

        labs_eval <- labs_all[idx]
        if (!all(valid[idx])) {
          resolution_eval$status[i] <- "missing_or_reference_labels_in_sample"
          next
        }
        labs_eval <- as.integer(factor(labs_eval))
        k_eval <- length(unique(labs_eval))
        resolution_eval$n_eval_clusters[i] <- k_eval
        if (is.null(d) || k_eval < 2L || k_eval >= length(idx)) {
          resolution_eval$status[i] <- "silhouette_undefined"
          next
        }
        if (k_eval != k) {
          resolution_eval$status[i] <- "sample_missing_clusters"
          next
        }
        score <- mean(cluster::silhouette(labs_eval, d)[, "sil_width"])
        if (!is.finite(score)) {
          resolution_eval$status[i] <- "nonfinite_silhouette"
          next
        }
        resolution_eval$silhouette[i] <- score
        resolution_eval$status[i] <- "ok"
        flog.info(sprintf("Resolution = %g; clusters = %d; silhouette = %.6f",
                          res, k, score))

        # Keep only the best full find_subclust result, not just its labels.
        # Because resolutions are sorted, equal scores favor lower resolution.
        if (score > best_score) {
          best_score <- score
          best_resolution <- res
          cell_info_ls <- tmp
        }
      }

      if (!is.finite(best_score)) {
        if (is.null(fallback_result)) {
          stop("find_subclust produced no usable clustering for any resolution.")
        }
        # A one-cluster result has no defined silhouette, for example.
        warning("No valid silhouette score; using the smallest successful resolution.")
        cell_info_ls <- fallback_result
        best_resolution <- fallback_resolution
      }
      resolution_eval$selected <- resolution_eval$resolution == best_resolution
      flog.info(sprintf("Selected resolution: %g; number of clusters: %d; silhouette: %s",
                        best_resolution,
                        resolution_eval$n_clusters[resolution_eval$selected],
                        format(resolution_eval$silhouette[resolution_eval$selected], digits = 6)))

    }
    cell_info <- cell_info_ls$meta_cluster 



  }else{
    cell_info <- cellMeta
    cell_info$subCluster <- NA
    cell_info[colnames(clt_matrix),"subCluster"] <- "0"
    cell_info_ls <- list(meta_cluster = cell_info)
  }

  flog.info("Number of cells in clusters:")
  flog.info(paste(table(cell_info$subCluster),collapse = ","))
  
  cellMeta$subCluster <- NA
  row_idx <- match(rownames(cellMeta),rownames(cell_info))
  #cellMeta$subCluster <- as.character(cellMeta$subCluster)
  cellMeta$subCluster[!is.na(row_idx)] <- as.character(cell_info$subCluster[na.omit(row_idx)])
 # cellMeta$subCluster[rownames(cellMeta)%in%cells_filered] <- "removed"
  cellMeta$subCluster[is.na(cellMeta$subCluster)] <- "reference"
  cellMeta <- suppressWarnings(cellMeta[order(as.numeric(cellMeta$subCluster)),])
  
  if(onRatio & noRefShow){
   # mt_plot <- mtx_bin #with normal cells
    mt_plot <- mtx_bin[,cells_t,drop=F]
    
  }else{
    # mtx_bin_ref <- suppressMessages(smooth_by_chromosome(inputMat[,cells_ref],wl=window,method="winMean",win_step=window,adjust_zero =TRUE))
    # mt_plot <- cbind(mtx_bin,mtx_bin_ref)
    mt_plot <- mtx_bin
  }
  if(!is.null(chrs)){
    mt_plot <- mt_plot[grepl(paste0("^",chrs,delimit,collapse = "|"),rownames(mt_plot)),,drop=F]
  }
  
  
  if(outputFigures){
    left_anno_cols <- list()
    color_r <- suppressWarnings(get_group_color_palette("Set3")(length(unique(cellMeta$subCluster))))
    names(color_r) <- sort(unique(cellMeta$subCluster))
    left_anno_cols[["subCluster"]] <- color_r
    if(is.null(Figure_name)){
      fileout_nameC = paste0("heatmap_subClust_binCount")
      fileout_nameR = paste0("heatmap_subClust_binRatio")
    }else{
      fileout_nameC = paste0(Figure_name,"_binCount")
      fileout_nameR = paste0(Figure_name,"_binRatio")
    }
    if(plt4peak & is.null(Figure_name)){
      fileout_name2C = paste0("heatmap_subClust_peakCount")
    }else{
      fileout_name2C = paste0(Figure_name,"_peakCount")
    }

    # saveRDS(mt_plot,paste0(outdir,"/binmtx.rds"))
    # write.table(cellMeta,paste0(outdir,"/cell_anno_subClust.txt"),col.names = T,row.names = T,quote=F,sep="\t")
    #outplot_path <- paste0(outdir,"/Figures");if(!file.exists(outplot_path)){dir.create(outplot_path,recursive=T)}
    if(noRefShow){
      ano <- cellMeta[!cellMeta$subCluster %in%"reference",]
      ano$subCluster <- factor(ano$subCluster,levels = sort(as.numeric(unique(ano$subCluster))))
    }else{
      ano <- cellMeta
      num_label <- suppressWarnings(sort(as.numeric(unique(ano$subCluster))))
      ano$subCluster <- factor(ano$subCluster,levels = c(num_label,"reference"))
    }
    
     pt <- heatmap4peakMt(mt_plot[,rownames(ano)],meta_info=ano[,c("subCluster"),drop=F],
                       max_lim=NULL,
                       value.type = "count",max.value=count_lim,
                       legend_titles="bin count", width=Figure.width,height=Figure.height,
                       col_list=left_anno_cols,
                       clust_rows=F,clustering_method_rows="ward.D2",sep_by="_|-|:",
                       outdir=outdir,
                       fileout_name = fileout_nameC)
    if(is.null(chrs)){plt4peak=FALSE}
    if(plt4peak){
      if(!is.null(chrs)){
       inputMat_chr <- inputMat[grepl(paste0("^",chrs,delimit,collapse = "|"),rownames(inputMat)),,drop=F]
      }else{
        inputMat_chr = inputMat
      }
      pt <- heatmap4peakMt(inputMat_chr[,rownames(ano)],meta_info=ano[,c("subCluster"),drop=F],
                           max_lim=NULL,
                           value.type = "any",legend_titles="peak count", width=Figure.width,height=Figure.height,
                           col_list=left_anno_cols,
                           clust_rows=F,clustering_method_rows="ward.D2",sep_by="_|-|:",
                           outdir=outdir,
                           fileout_name = fileout_name2C)
    }
    if(onRatio){
    ratio_res <- as.matrix(ratio_res)
    ano2 <- cellMeta[!cellMeta$subCluster %in%"reference",]
    ano2$subCluster <- factor(ano2$subCluster,levels = sort(as.numeric(unique(ano2$subCluster))))
    
    mt_plt <- ratio_res[,rownames(ano2)]
    sm_win <- ifelse(nrow(mt_plt)>1e4,ceiling(20*nrow(mt_plt)/1e4),10)
    
    mt_plt <- smooth_and_denoise(mt_plt,window=sm_win)
    
    pt <- heatmap4peakMt(mt_plt,meta_info=ano2[,c("subCluster"),drop=F],
                         max_lim=NULL,
                         value.type = "ratio",legend_titles="bin ratio", width=Figure.width,height=Figure.height,
                         col_list=left_anno_cols,
                         clust_rows=F,clustering_method_rows="ward.D2",sep_by="_|-|:",
                         outdir=outdir,
                         fileout_name = fileout_nameR)

    }
    
  }
  if(OutMat){
    binCount_OutName <- paste0(outdir,"/binCount_Mtx.rds")
    saveRDS(mtx_bin,binCount_OutName)
    if(onRatio){
      ratio_res <- as.matrix(ratio_res)
      ano2 <- cellMeta[!cellMeta$subCluster %in%"reference",]
      binRatio_OutName <- paste0(outdir,"/binRatio_Mtx.rds")
      saveRDS(ratio_res[,rownames(ano2)],binRatio_OutName)
    }
  }

  binRatio <- if (onRatio) {
    ratio_cells <- rownames(cellMeta)[!cellMeta$subCluster %in%"reference"]
    ratio_res[, ratio_cells, drop = FALSE]
  } else {
    NULL
  }



  return(list(cellMeta=cellMeta,binCount=mtx_bin,binCount.nom=countsMnormal,binRatio=binRatio,
    best_resolution = best_resolution,resolution_eval = resolution_eval
    ))
}


  
    
    
    
  
