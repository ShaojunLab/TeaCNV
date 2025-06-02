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
                     global_norm = FALSE
                     ) {
  suppressMessages({
    library(ggplot2)
    library(dplyr)
    library(stringr)
  })
  appender_func <- flog.logger()$appender
  log_path <- environment(appender_func)$file
  if((!is.null(log_path)) & file.exists(log_path)){flog.appender(appender.file(log_path))}else{
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

    cell_info_ls <- suppressMessages(find_subclust(clt_matrix,resolution = seu_resol,sep_by = c("-", "-")))
    cell_info <- cell_info_ls$meta_cluster 
  }else{
    cell_info <- cellMeta
    cell_info$subCluster <- NA
    cell_info[colnames(clt_matrix),"subCluster"] <- "0"
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
                       clust_rows=F,clustering_method_rows="ward.D2",sep_by="-",
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
                           clust_rows=F,clustering_method_rows="ward.D2",sep_by="-",
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
                         clust_rows=F,clustering_method_rows="ward.D2",sep_by="-",
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
  if(onRatio){
    return(list(cellMeta=cellMeta,binCount=mtx_bin,binCount.nom=countsMnormal,binRatio=ratio_res[,rownames(ano2)]))
  }else{
    return(list(cellMeta=cellMeta,binCount=mtx_bin,binCount.nom=countsMnormal,binRatio=NULL))
  }
}


  
    
    
    
  
