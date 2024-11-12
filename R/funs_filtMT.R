#!/usr/bin/env Rscript

# This is a program to filter peaks from scATAC count matrix.
suppressMessages({
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(ggsci)
  library(plyranges)
  library(GenomicRanges)
  library(futile.logger)
})


#' @title blacklist.rm()
#' @description remove peaks in blacklist
#' @param mat count matrix with peaks on rows, cells on columns. 
#' Peaks with "chrx-xxx-xxx" format were default separated with "-"
#' @param blacklist bed format file with path
#' @export
#' 
blacklist.rm <- function(mat,blacklist,peak.sep="-"){
  #mat <- mat[grepl(paste(paste0("chr",c(seq(1:22),"X")),collapse="|"),rownames(mat)),]
  if (Reduce("|", is(blacklist) == "character")) {
    if(substr(blacklist, nchar(blacklist)-3, nchar(blacklist)) == ".bed"){
      blacklist <- read_bed(blacklist)
    }else{
      stop("Cannot find a blacklist file with the suffix '.bed'.")
    }
  }
  bins <- data.frame(peak=rownames(mat))
  #remove bins within blacklist
  chromInfo <- separate(bins, peak, into = c("seqnames", "start","end"), sep = peak.sep)
  bins <- data.frame(bins,chromInfo)
  #make a Grange object
  peak_gr <- GRanges(seqnames = bins$seqnames,
                     ranges = IRanges(start = as.numeric(bins$start), end = as.numeric(bins$end)),
                     strand = "*")
  # Find the intersection
  message("Subtracting Blacklist...")
  inblack <- findOverlaps(peak_gr, blacklist)
  idx <- setdiff(1:length(peak_gr), S4Vectors::queryHits(inblack))
  # Remove the blacklist regions from peak_gr
  peak_gr_use <- peak_gr[idx]
  keep_seqnames <- seqnames(peak_gr_use)
  keep_ranges <- ranges(peak_gr_use)
  bins_kp <- paste(keep_seqnames,start(keep_ranges),end(keep_ranges),sep=peak.sep)
  
  mtx_filt <- mat[rownames(mat)%in%bins_kp,]
 # mtx_filt <- as.matrix(mtx_filt)
  return(mtx_filt)
}


#' @title filt_peak_perChr()
#' @description Remove peaks with abnormally high total count in all cells per chromosome
#' @param mat peak-cell count matrix
#' @param zscore.lim the Zscore cutoff value for distribution of peak total count
#' @param split_by the delimiter of peak name, e.g. "chrx-xxx-xxx"
#' @export
filt_peak_perChr <- function(mat,zscore.lim=2,split_by="-"){
  suppressMessages({ require(data.table)})
  raw_peaks <- rownames(mat)
  
  mat_nom <- t(t(mat)*mean(colSums(mat))/colSums(mat)) 
  
  bed_data <- tstrsplit(raw_peaks, split_by)
  bed_data <- as.data.frame(bed_data)
  colnames(bed_data) <- c("chromosome", "start", "end")
  chroms <- unique(bed_data$chromosome)
  peaks_keep <- c()
  for(chrj in chroms){
    mtx_chr <- mat_nom[grepl(paste0("^",chrj,split_by),rownames(mat_nom)),,drop=F]
    if(ncol(mtx_chr)>1){
      totalcount <- rowSums(mtx_chr)
    }else{
      totalcount <- mtx_chr[,1]
      names(totalcount)<- rownames(mtx_chr)
    }
    zscore <- scale(totalcount)
    peaks <- names(totalcount)[abs(zscore) <= zscore.lim]
    peaks_keep <- c(peaks_keep,peaks)
  }
  #order peaks
  rowInfo <- tstrsplit(peaks_keep, split_by)
  row <- as.data.frame(rowInfo)
  colnames(row) <- c("chromosome", "start", "end")
  row <- row[grepl("chr",row$chromosome),]
  row$chromnum <- sapply(strsplit(row$chromosome,'hr'),'[',2)
  row <- row[order(as.integer(gsub("[^0-9]", "", row$chromnum)),as.numeric(row$start)),]
  row$peak <- paste(row$chromosome,row$start,row$end,sep=split_by)
  mat_filt <- mat[row$peak,,drop=F]
  return(mat_filt)
  
}


#' @title filt_by_prop()
#' @export

filt_by_prop <- function(cellprop.lim,cell_anno=NULL,ref_group_names=NULL,mtx_qc_filt0){
  if(!is.null(cell_anno)){
    min_cells_prop_t <- cellprop.lim
    min_cells_prop_n <- cellprop.lim
    if(!is.null(ref_group_names)){
      cells_t <- rownames(cell_anno)[!cell_anno[,1]%in%ref_group_names]
      cells_n <- rownames(cell_anno)[(cell_anno[,1]%in%ref_group_names)]
      
      # chromnum=sapply(strsplit(rownames(mtx_qc_filt0),'-|_|:'),'[',1)
         
      peaks_keep <- c()
      if(length(cells_t)>0){
        # mtx_t <- mtx_qc_filt0[,cells_t]

        peak_t_keep <- rownames(mtx_qc_filt0)[rowSums(mtx_qc_filt0[,cells_t]>0)>min_cells_prop_t*length(cells_t)]
        peaks_keep <- peak_t_keep
      }
      if(length(cells_n)>0){
        #remove outlier peaks 
        # nCount_peak <- rowSums(mtx_qc_filt0[,cells_n],na.rm=TRUE)
        # nCount_cutoff <- quantile(nCount_peak,c(0.05,0.95))
        # peaks_keep_n <- which(nCount_peak>=nCount_cutoff[1] & nCount_peak<= nCount_cutoff[2] )
        # mtx_n <- mtx_qc_filt0[,cells_n][peaks_keep_n,,drop=F]
        mtx_n <- mtx_qc_filt0[,cells_n]
        peak_n_keep <- rownames(mtx_n)[rowSums(mtx_n>0)>min_cells_prop_n*length(cells_n)]
        peaks_keep <- intersect(peaks_keep,peak_n_keep)
      }
      
    }else{
      peaks_keep <- rownames(mtx_qc_filt0)[rowSums(mtx_qc_filt0>0)>cellprop*ncol(mtx_qc_filt0)]
    }
    
  }else{
    peaks_keep <- rownames(mtx_qc_filt0)[rowSums(mtx_qc_filt0>0)>cellprop*ncol(mtx_qc_filt0)]
  }
  return(peaks_keep)
}



#' @title FiltPeak()

#' @param blacklist_file .bed file of blacklist
#' @param cell_anno annotation of cells with cell groups information on the first column.
#' @export
FiltPeak <- function(
    mtx,
    cell_anno=NULL,
    ref_group_names=NULL,
    cellprop = 0.05,
    outheatmap = FALSE,
    blacklist_file = NULL,
    ChrRemove = c('chrX', 'chrY', 'chrM'),
    prop_reset=TRUE
    ) {
  ###1 remove peaks from blacklist
  if(!is.null(blacklist_file)){
    mtx_qc_filt0 <- blacklist.rm(mtx,blacklist_file)
    }else{
     mtx_qc_filt0 <- mtx
    }
  peaks.raw <- data.frame(peak=rownames(mtx_qc_filt0))
  row.sep <- separate(peaks.raw, peak, into = c("seqnames", "start","end"), sep = "-")
  peaks.raw <- cbind(peaks.raw,row.sep)
  #remove peaks from chrY,chrM
  if(!is.null(ChrRemove) & any(peaks.raw$seqnames %in% ChrRemove)){
    peaks.raw <- peaks.raw[-which(peaks.raw$seqnames %in% ChrRemove),,drop=F]
    mtx_qc_filt0 <- mtx_qc_filt0[peaks.raw$peak,]
  }
  ###2. Filtering peaks by the number (or proportion) of cells that the peak showing in the dataset
  print(paste0("Number of raw cells: ",ncol(mtx)))
  

  peaks_keep <- filt_by_prop(cellprop.lim=cellprop,cell_anno,ref_group_names,mtx_qc_filt0)
  
  if(length(peaks_keep)==0){
    stop("No peak retained after filtering, please check data and parameter 'cellprop'.")
  }
  if(prop_reset){
    if(nrow(mtx_qc_filt0)>1e4 & length(peaks_keep)<1e4 & cellprop>=0.02){
      message("The number of filtered peaks is too small. Reset the 'cellprop' to a smaller value. ")
      while(length(peaks_keep)<1e4){
        cellprop <- cellprop- 0.005
        peaks_keep <- filt_by_prop(cellprop.lim=cellprop,cell_anno,ref_group_names,mtx_qc_filt0)
        if(length(peaks_keep)>=1e4 | cellprop<= 0.01){break}
      }
      print(paste0("Reset 'cellprop'= ",cellprop))
    }
    
  }
    #order peaks by chromosome
  row <- data.frame(peak=peaks_keep)
  rowInfo <- separate(row, peak, into = c("seqnames", "start","end"), sep = "-")
  row <- cbind(row,rowInfo)
  row$chromnum <- sapply(strsplit(row$seqnames,'hr'),'[',2)
  row <- row[order(as.integer(gsub("[^0-9]", "", row$chromnum)),as.numeric(row$start)),]
  

  mtx_qc_filt <- mtx_qc_filt0[row$peak,]
  print("Number of peaks and cells after filtering: ")
  print(dim(mtx_qc_filt))
  print(paste0(nrow(mtx_qc_filt)," (",round(100*nrow(mtx_qc_filt)/nrow(mtx_qc_filt0),2)," %) peaks retained by keeping ",cellprop*100,"% in provided cells"))
  
  mtx_qc_filt <- as.matrix(mtx_qc_filt)
  
  #saveRDS(mtx_qc_filt,"mtx_qc_filt.rds")
   
  # #heatmap for peaks
  if(outheatmap){
    fileout_name = paste0("heatmap_PeakCount.pdf")
    heatmap4peakMt(mat=mtx_qc_filt, max_lim= FALSE,sep_by="-",
                   outdir= outdir,value.type="count",fileout_name=fileout_name,
                   width=10,height=6)
  }
  print("Fitering peaks done!")
  
  return(mtx_qc_filt)
}

#' @title FiltCell.mt()
#' @export
FiltCell.mt <- function(mtx,
                     nCount_quantile_lim=c(0.05,0.95),
                     nFeature_quantile_lim=c(0.05,0.95)){
  #if(!is.matrix(mtx)){mtx <- as.matrix(mtx)}
  nCount_per_cell = colSums(mtx)
  nFeature_per_cell = colSums(mtx>0)
  
  low_ncount <- quantile(nCount_per_cell, probs = nCount_quantile_lim[1],na.rm=T)
  hig_ncount <- quantile(nCount_per_cell, probs = nCount_quantile_lim[2],na.rm=T)
  cells_keep1 <- which(nCount_per_cell>=low_ncount & nCount_per_cell<= hig_ncount)
  
  low_nFeature <- quantile(nFeature_per_cell, probs = nFeature_quantile_lim[1],na.rm=T)
  hig_nFeature <- quantile(nFeature_per_cell, probs = nFeature_quantile_lim[2],na.rm=T)
  cells_keep2 <- which(nFeature_per_cell>=low_nFeature & nFeature_per_cell<= hig_nFeature)
  cells_keep <- intersect(cells_keep1,cells_keep2)
  mtx_filt <- mtx[,cells_keep]
  
  n_raw_cells <- ncol(mtx)
  n_remove <- n_raw_cells - length(cells_keep)
  
  #flog.info(sprintf("Removing %g (%g %%) of cells", n_remove, round(n_remove/n_raw_cells * 100,2)))
  
  return(mtx_filt)
}


#' @title FiltCell.obj()
#' @return list of filtered cells and plots
#' @export
FiltCell.obj <- function(object,
                        assay="peaks",
                        cell_meta =NULL,
                        filt_by=paste0(c("nCount_","nFeature_"),assay),
                        filt_per_group = NULL,
                        by_density = TRUE,
                        by_quantile = FALSE,
                        by_zscore.x = TRUE,
                        zscore.min = -1,
                        nCount_quantile_lim=c(0.05,0.95),
                        nFeature_quantile_lim=c(0.05,0.95),
                        log_=TRUE){
  suppressMessages({
    library(MASS)
    library(ggplot2)
    library(viridis)
    library(progress)
    library(cowplot)
    library(ggExtra)
    })
  theme_set(theme_bw(base_size = 16))
  # Get density of points in 2 dimensions.
  # @param x A numeric vector.
  # @param y A numeric vector.
  # @param n Create a square n by n grid to compute density.
  # @return The density within each square.
  get_density <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
  }
  if(Reduce("|", is(object) %in% c("Seurat"))){
    mat <- object[[assay]]@counts
    if(is.null(cell_meta))cell_meta <- object@meta.data
  }else if(Reduce("|", is(object) %in% c("dgCMatrix", "matrix"))){
    mat <- object
    if(is.null(cell_meta)){
      cell_meta <- data.frame(row.names = colnames(mat),filtGroup=rep("total",ncol(mat))) 
    }
     
  }
  if(!is.null(filt_per_group)){
    cell_meta$filtGroup <- cell_meta[,filt_per_group]
  }else{
    cell_meta$filtGroup <- "total"
  }
  if(!any(grepl(paste(filt_by,collapse = "|"),colnames(cell_meta)))){
    message(paste0("Use nCount and nFeature calculated from data instead."))
    nCount_per_cell = colSums(mat)
    nFeature_per_cell = colSums(mat>0)
    filt_by <- c("nCount","nFeature")
    cell_meta$nCount <- nCount_per_cell
    cell_meta$nFeature <- nFeature_per_cell
    
  }
Group_size <- table(cell_meta$filtGroup)
Group_size <- Group_size[Group_size>10]

  Groups <- names(Group_size)
  plt_ls <- list()
  cells_keep <- c()
  pb <- progress_bar$new(total = length(Groups), format = "[:bar] :percent :elapsed")
  for(i in 1:length(Groups)){
    grp <- Groups[i]
    cells_select <- rownames(cell_meta)[cell_meta$filtGroup == grp]
    mtx_sub <- mat[,cells_select]
    cell_anno_grp <- cell_meta[cells_select,,drop=F]
    

    #check data QC
    # binary_mat <- mtx_sub >0
    # #nFeature_peaks
    # nFeature_per_cell = colSums(binary_mat)
    # nCells_per_peak = rowSums(binary_mat)
    # binary_mat_overlapPeaks <- binary_mat[nCells_per_peak>1,]
    # overlapPeaks.N <- colSums(binary_mat_overlapPeaks)
    # FractionOverlapPeaks_per_cell <- overlapPeaks.N/nrow(mtx_sub)
    # cell_anno_grp$FractionOverlapPeaks_per_cell <- FractionOverlapPeaks_per_cell
    # cell_anno_grp$nFeature_per_cell <- nFeature_per_cell
    # p2 <- ggplot(cell_anno_grp,aes(FractionOverlapPeaks_per_cell)) + geom_density(alpha=.4, fill="#999999")+xlab("Fraction overlaping peaks per cell")
    # p1 <- ggplot(cell_anno_grp,aes(log10(nFeature_per_cell))) + geom_density(alpha=.4, fill="#999999")+xlab("Log nFeature per cell")
    # plot_den <- cowplot::plot_grid(p1, p2,ncol = 1)
    # title <- cowplot::ggdraw() +
    #   cowplot::draw_label(grp,fontface = 'bold',hjust = 0.5) +
    #   theme(plot.margin = margin(0, 0, 0, 0))
    # plot_den_com <- cowplot::plot_grid(title, plot_den,ncol = 1,rel_heights = c(0.1, 1),align = "hv")

    cells_keep.grp <- rownames(cell_anno_grp)
    
    # if(by_fraction){
    #   if(is.null(fraction.min)){
    #     fraction.min <- find_1st_Trough(cell_anno_grp$FractionOverlapPeaks_per_cell,log_scale = F)
    #   }
    #   cells_keep.grp <- rownames(cell_anno_grp)[cell_anno_grp$FractionOverlapPeaks_per_cell >= fraction.min]
    #   mtx_sub <- mtx_sub[,cells_keep.grp]
    # }
    
    if(by_zscore.x){
      x = cell_anno_grp[,filt_by[1]]
      x_scale <- scale(x)
      row_idx <- x_scale>=zscore.min
      if(all(row_idx)){
        #cutoff <- find_1st_Trough(x,log_scale = log_)
        row_idx <- x_scale >= 0
      }

      cells_keep.z <- rownames(cell_anno_grp)[row_idx]
      cells_keep.grp <- intersect(cells_keep.grp,cells_keep.z)
      mtx_sub <- mtx_sub[,cells_keep.grp]
    }
    
    if(by_quantile){
      mtx_sub.q <- suppressWarnings({FiltCell.mt(mat[,cells_select],nCount_quantile_lim,nFeature_quantile_lim)})
      cells_keep.q <- colnames(mtx_sub.q)
      cells_keep.grp <- intersect(cells_keep.grp,cells_keep.q)
      mtx_sub <- mtx_sub[,cells_keep.grp]
    }

    
    cell_anno_sub <- cell_anno_grp[cells_keep.grp,,drop=F]
    cell_ano_rm <- cell_anno_grp[!rownames(cell_anno_grp)%in%cells_keep.grp,,drop=F ]
    dat_rm <- data.frame(x=as.numeric(cell_ano_rm[,filt_by[1]]),y=as.numeric(cell_ano_rm[,filt_by[2]]))
    dat <- data.frame(x=as.numeric(cell_anno_sub[,filt_by[1]]),y=as.numeric(cell_anno_sub[,filt_by[2]]))

    if(log_){
      dat$x <- log10(dat$x)
      dat_rm$x <- log10(dat_rm$x)
      xlabel <- paste0("Log ",filt_by[1])
      dat$y <- log10(dat$y)
      dat_rm$y <- log10(dat_rm$y)
      ylabel <- paste0("Log ",filt_by[2])
    }else{xlabel = filt_by[1];ylabel = filt_by[2]}
    
    dat$density <- get_density(dat$x, dat$y, n = 100)
    if(nrow(dat_rm)>0){
      dat_rm$density <- NA
      dat2 <- rbind(dat,dat_rm)
    }else{
      dat2 <- dat 
    }

    
    gg <- ggplot(dat2,aes(x, y)) + geom_point(aes(color = density)) + scale_color_viridis()+
      labs(x=xlabel,y=ylabel,title = grp)
    
    if(by_density){
      level.min <- quantile(dat2$density,0.18,na.rm=T)
      level.max <- max(dat2$density)

      gg <- gg+
        geom_density_2d(data=dat,aes(x, y),breaks = c(level.min, level.max)) 
      
      gb <- ggplot_build(gg)

      level_cutoff <- min(gb$data[[2]]$level)
      row_ind <- which(dat$density>=level_cutoff)
      cells_keep.grp <- rownames(cell_anno_sub)[row_ind]
    }else{
      cells_keep.grp <- rownames(cell_anno_sub)
    }
    

    plot_com <- ggMarginal(gg, type = "density",margins="x")
    
    plt_ls[[as.character(grp)]] <- plot_com
    cells_keep <- c(cells_keep,cells_keep.grp)
    pb$tick()
  }
  
  n_raw_cells <- ncol(mat)
  n_remove <- n_raw_cells - length(cells_keep)
  flog.info(sprintf("Removing %g (%g %%) of cells", n_remove, round(n_remove/n_raw_cells * 100,2)))
  
  res <- list(cell_keep=cells_keep,plots = plt_ls)
  return(res)
}
  



#' @title cell_score()
#' @export

cell_score <- function(data,ZeroProp.max=0.5,byBin=TRUE,binWin=20,outFigure=TRUE,outdir="./",cell_meta=NULL,max.value=NULL){
  suppressMessages({
    library(dplyr)
    library(tidyr)
    library(MatrixGenerics)
  })
  if(outdir != "." && !file.exists(outdir)){
    flog.info(paste0("Creating output path ", outdir))
    dir.create(outdir,recursive=T)
  }
  data <- as.matrix(data)
  if(byBin){
    mtx_bin <- suppressMessages(smooth_by_chromosome(data,wl=binWin,method="winMean",win_step=binWin,adjust_zero =TRUE))
    data <- mtx_bin
  }
  chromInfo_sub <- data.frame(peak=rownames(data))
  chromInfo_sub <- separate(chromInfo_sub, peak, into = c("chrom", "start","end"), sep = "-|:|_")
  chrom <- unique(chromInfo_sub$chrom)
  cell_scores <- data.frame(cell = colnames(data), score = numeric(ncol(data)))
  chrNz <- c()
  for (i in 1:length(chrom)) {
    chrj <- chrom[i]
    mtx_sub_chrj <- data[grepl(paste0(paste0("^",chrj),c("-",":","_"),collapse = "|"),rownames(data)),]
    Nz <- colSums(mtx_sub_chrj==0)
    ProZ <- Nz/nrow(mtx_sub_chrj)
    chrNz <- rbind(chrNz,ProZ)
  }
  maxNz <- apply(chrNz, 2, max)
  #maxNz <- MatrixGenerics::colMaxs(chrNz)
  cell_scores$score <- maxNz
  cell_scores$label <- ifelse(maxNz<=ZeroProp.max,"Low","Hig")
  rownames(cell_scores) <- cell_scores$cell
  
  if(outFigure){
    # p <- ggplot(cell_scores,aes(x=score, fill=label)) +
    #   geom_density(alpha=0.6)+
    #   scale_fill_manual(values=c("#8338ec","#3a86ff")) +
    #   theme_minimal() +
    #   labs(fill="") 
    # ggsave(paste0(outdir,"/density_CellScore.pdf"),p,width = 5,height = 4)
    
    left_anno_cols <- list()
    color_r <- suppressWarnings(get_group_color_palette("Set3")(length(unique(cell_scores$label))))
    names(color_r) <- unique(cell_scores$label)
    left_anno_cols[["label"]] <- color_r
    fileout_nameC = paste0("heatmap_binCount_",ZeroProp.max)
    
    if(!is.null(cell_meta)){
      matched_indices <- match(rownames(cell_meta),cell_scores$cell)
      cell_meta$label <- cell_scores$label[matched_indices]
      cell_scores_ord <- cell_meta[order(cell_meta$label,cell_meta[,1]),]
    }else{
      cell_scores_ord <- cell_scores[order(cell_scores$label),c("label"),drop=F]
    }
    
    pt <- heatmap4peakMt(data[,rownames(cell_scores_ord)],meta_info=cell_scores_ord,
                         max_lim=NULL,
                         value.type = "count",max.value=max.value,
                         legend_titles="bin count", width=10,height=6,
                         col_list=left_anno_cols,
                         clust_rows=F,clustering_method_rows="ward.D2",sep_by="-",
                         outdir=outdir,
                         fileout_name = fileout_nameC)
  }
  
  return(cell_scores)
}


#' @title cell_VarScore()
#' @description score cell based on the arm-levle weighted standard deviation
#' @param data peak- cell matrix (row:peak,column:cell)
#' @param genomic_ref Bed file of with genomic range for arms of chromosome (e.g. cytoBand_hg38.tsv)
#' @export
cell_VarScore <- function(data,genomic_ref=NULL,Zscore.Label=1){
  suppressMessages({
    library(dplyr)
    library(tidyr)
  })
  if(!is.null(genomic_ref)){
    if (Reduce("|", is(genomic_ref) == "character")) {
      if(substr(genomic_ref, nchar(genomic_ref)-2, nchar(genomic_ref)) %in% c(".gz","txt","tsv","csv")){
        genomic_ref <- read.table(connection <- gzfile(genomic_ref, 'rt'), sep="\t", header=F, check.names=FALSE)
        close(connection)
        genomic_ref <- as.data.frame(genomic_ref)
      }
    }else if(Reduce("|", is(genomic_ref) %in% c("data.frame"))){
    }else{
      message("CytoBand isn't recognized as data.frame, or filename")
    }
  }
    colnames(genomic_ref)=c("chrom","chromStart","chromEnd","name","gieStain")
    genomic_pq <- genomic_ref[grep("acen",genomic_ref$gieStain),]
    genomic_pq <- genomic_pq %>%
      group_by(chrom)%>%
      top_n(1,wt=chromStart)%>%as.data.frame()
    
  calc_arm_variation <- function(arm_data) {
    apply(arm_data, 2, sd) 
  }
  
  chromInfo_sub <- data.frame(peak=rownames(data))
  chromInfo_sub <- separate(chromInfo_sub, peak, into = c("chrom", "start","end"), sep = "-|:|_")
  chrom <- unique(chromInfo_sub$chrom)
  
  #Initialize an empty DataFrame to store the total score of each cell
  cell_scores <- data.frame(cell = colnames(data), score = numeric(ncol(data)))
  armSD <- c()
  metadata <- c()
  for (i in 1:length(chrom)) {
    chrj <- chrom[i]
    genomic_pos <- genomic_pq[genomic_pq$chrom==chrj,"chromStart"]
    mtx_sub_chrj <- data[grepl(paste0(paste0("^",chrj),c("-",":","_"),collapse = "|"),rownames(data)),]
    peak_end <- sapply(strsplit(rownames(mtx_sub_chrj),"-|:|_"),"[",3)
    for(arm in c("p","q")){
      if(arm == "p"){
           peak_p_idx <- as.numeric(peak_end)<=genomic_pos
         }else{
           peak_p_idx <- as.numeric(peak_end)>genomic_pos
         }
      mtx_sub_chrj_arm <- mtx_sub_chrj[peak_p_idx,,drop=F]
      arm_variation <- calc_arm_variation(mtx_sub_chrj_arm)
      armSD <- rbind(armSD,arm_variation)
      genomic_info <- c(chrj,arm,nrow(mtx_sub_chrj_arm))
      metadata <- rbind(metadata,genomic_info)
    }
  }
  metadata <- as.data.frame(metadata)
  colnames(metadata) <- c("chr","arm","NFeature")
  weights <- as.numeric(metadata$NFeature)
  weighted_sds <- apply(armSD, 2, weighted.mean, weights)
  cell_scores$score <- weighted_sds
  cell_scores$Zscore <- scale(cell_scores$score)
  cell_scores$label <- ifelse(cell_scores$Zscore <= Zscore.Label,"Low","Hig")
  return(cell_scores)
}



#' @description find the intersect of two genomic regions
#' @title bedtools_intersect()
#' @export

bedtools_intersect <- function(stringsA,stringsB,wa=TRUE,sepA = c("-", "-"),sepB=c("-", "-")){
  suppressMessages({
    library(Signac)
    library(dplyr)})
  grA <- StringToGRanges(stringsA,sep=sepA)
  grB <- StringToGRanges(stringsB,sep=sepB)
  grB$length <- ranges(grB)@width
  
  Ovr <- findOverlaps(grA, grB,type="any",select="all",ignore.strand=T)
  idxA <- S4Vectors::queryHits(Ovr)
  idxB <- S4Vectors::subjectHits(Ovr)
  #If multiple subjectHits for same  queryHits, select the longer one
  DF <- data.frame(A=stringsA[idxA],B=stringsB[idxB],widthB=grB$length[idxB])
  DF <- DF %>%
    group_by(A) %>%
    filter(widthB == max(widthB))%>%
    as.data.frame()
  
  dfA <- data.frame(A=stringsA)
  dfA <- left_join(dfA,DF[,1:2],by="A",keep=F)
  if(!wa){
    dfA <- na.omit(dfA)
  }
  return(dfA)
}

