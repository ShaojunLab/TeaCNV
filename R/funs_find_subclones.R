
#' @title CNratio_peak()
#' @description calculate peak-cluster level CN ratio
#' @param mt combined single cell matrix
#' @param cell_ano a data.frame has same number of cells on row, cell ID (first column) and sub-cluster annotation (second column) on columns
#' This function will convert single-cell matrix to cluster-level matrix first based on sub-cluster annotation 
#' @param colname_ref the name of group taken as reference.
#' @param ref.value numeric data.frame with the length of nrow(mat)
#' @export
CNratio_peak <- function(mt,cell_ano,ref.value=NULL,colname_ref="adjacent",Norm_by_ref=FALSE){
  mt_b <- pseudo_bulk_v2(mt,group_ano = cell_ano,method ="mean",adjust_zero =TRUE)
  # ratio matrix
  sub_cl <- colnames(mt_b)
  ratio_res <- data.frame(matrix(nrow = nrow(mt_b), ncol = 0))
  if(!is.null(ref.value)){
    if(nrow(ref.value)!=nrow(mt_b)){
      message("The nrow of 'ref.value' is not equal to nrow of 'mt'!")
      stop()
    }
    
  }else{
    ref.value <- mt_b[,colname_ref,drop=F]
  }
  if(Norm_by_ref){
    mt_b <- t(t(mt_b)*sum(ref.value)/colSums(mt_b)) 
  }
  
  mt_b.obs <- as.matrix(mt_b[,!colnames(mt_b)%in%colname_ref])
  ratio_res <- sweep(mt_b.obs,1,ref.value[,1],FUN="/")
  
  
  # for(j in 1:length(sub_cl)){
  #   c2 <- as.character(sub_cl[j])
  #   if(!c2 %in%  colname_ref){
  #     r <- data.frame(mt_b[,which(colnames(mt_b)%in%c2),drop=F]/ref.value)
  #     colnames(r) <- c2
  #     ratio_res <- cbind(ratio_res,r)
  #   }
  # }
  colnames(ratio_res) <- colnames(mt_b)[!colnames(mt_b)%in%colname_ref]
  rownames(ratio_res) <- rownames(mt_b)
  return(ratio_res)
}

#' @title count2ratio()
#' @param mt_obs and mt_ref: peak-cell count matrix of tumor and normal reference, which rows should be the same
#' @param win The number of features to be combined into one bin.
#' @param colname_ref the name of group taken as reference.
#' @param chr_select This is only needed when performed on a specific chromosome.
#' @return cell_info a data.frame with cell ID and sub-cluster annotation on columns
#' @export
#source("mydataProcess.R",chdir = T) 
count2ratio <- function(mt_obs,mt_ref,colname_ref="adjacent",win=10,chr_select=NULL){
  suppressMessages({
    library(Signac)
    })
  mtx_sub_lim <- as.matrix(mt_obs)
  bound_up <- round(quantile(mtx_sub_lim,0.99));bound_up
  mtx_sub_lim[mtx_sub_lim>bound_up] <- bound_up
  
  bound_up <- round(quantile(mt_ref,0.99));bound_up
  mt_ref[mt_ref>bound_up] <- bound_up
  mt_ana <- cbind(mtx_sub_lim,mt_ref) 
  
  mtx_bin <- smooth_by_chromosome(mt_ana,wl=win,method="winMean",win_step=win,adjust_zero =TRUE)
  #normalization
  countsMnormal = t(t(mtx_bin)*mean(colSums(mtx_bin))/colSums(mtx_bin))
  #process cell annotatio dataframe
  cell_info_new <- data.frame(row.names = colnames(mtx_sub_lim),group=colnames(mtx_sub_lim))
  cell_info_ref <- data.frame(row.names = colnames(mt_ref),group=rep(colname_ref,ncol(mt_ref)))
  cell_info_new <- rbind(cell_info_new,cell_info_ref)
  cell_info2 <- data.frame(rownames(cell_info_new),cell_info_new)
  
  #generate pseudo bulk level of reference cells, generate one single cell matrix with one "cell" of bulk-reference mean count
  ratio_res <- CNratio_peak(mt = countsMnormal,cell_ano = cell_info2,colname_ref=colname_ref)
    
  if(!is.null(chr_select)){
    ratio_res_chr <- ratio_res[grepl(paste0("^",chr_select,"_|^",chr_select,"-|^",chr_select,":"),rownames(ratio_res)),]
  }else{
    ratio_res_chr <-ratio_res
  }
  bound_up <- round(quantile(as.matrix(ratio_res_chr),0.99));bound_up #6
  ratio_res_chr[ratio_res_chr>bound_up] <- bound_up
  
  return(ratio_res_chr)
}

#' 
#' @description  Find sub-cluster based on bin-level ratio using PCA method of seurat
#' @title find_subclust()
#' @export
find_subclust <- function(mt,add_metadata=NULL,resolution = 0.8,sep_by = c("-", "-"),assay="binCount"){
  suppressMessages({
    library(Signac)
  })
  set.seed(123)
  mt <- as.matrix(mt)
  ### 对 mt 基于PCA聚类 
  atac_assay <- CreateChromatinAssay(counts = mt, 
                                     sep = sep_by,
                                     genome = "GRCh38",
                                     min.cells = 1)
  
  obj_sub2 <- CreateSeuratObject(
    counts = atac_assay,
    assay = assay,
    meta.data = add_metadata) #obj_sub@meta.data[colnames(mtx_sub),,drop=F])
  
  obj_sub2 <- FindTopFeatures(obj_sub2, min.cutoff = 'q5')
  obj_sub2 <- ScaleData(obj_sub2,features = rownames(obj_sub2) )
  Ncells <- ncol(obj_sub2)
  if(Ncells<50 & Ncells>1){
    npcs.pca=floor(Ncells/2)
    dims.max=min(npcs.pca,15)
    n.neighbors = min(30,Ncells)
  }else{npcs.pca=50;dims.max=15;n.neighbors=30}
  obj_sub2 <- RunPCA(obj_sub2, features = VariableFeatures(object = obj_sub2),npcs=npcs.pca)
  #ElbowPlot(obj_sub2)
  obj_sub2 <- FindNeighbors(obj_sub2, dims = 1:dims.max)
  obj_sub2 <- FindClusters(obj_sub2, resolution = resolution,algorithm =1 )
  obj_sub2 <- RunUMAP(obj_sub2, reduction = 'pca',dims = 1:dims.max,n.neighbors=n.neighbors)
  
  cell_info <- data.frame(row.names = rownames(obj_sub2[[]]),
                          subCluster=obj_sub2$seurat_clusters)
  cell_info$subCluster <- factor(cell_info$subCluster,levels = sort(unique(as.numeric(as.character(cell_info$subCluster)))))
  cell_info <-  cell_info[order(cell_info$subCluster),,drop=F]
  return(list(object=obj_sub2, meta_cluster = cell_info))
}

#' reclustering for subset of matrix by specific cells, and update cell_info in the first column.
#' cell_info:
#'                    subCluster
#' AAACTCGGTCAGTGTT-1          0
#' AAACTCGTCTCGACAA-1          0
#' AAACTGCAGTTATGAG-1          0
#' @title reClustByGroup()
#' @export
reClustByGroup <- function(mat,cell_info,cells=NULL,resolution = 0.8,sep_by = c("-", "-"),assay="binCount"){
  cell_info[,1] <- as.character(cell_info[,1])
  if(is.null(cells)){cells=colnames(mat)}
   mat_new <- mat[,cells]
  newClust <- find_subclust(mat_new,resolution=resolution,sep_by=sep_by,assay=assay)
  info_new <- newClust$meta_cluster
  row_idx <- match(rownames(info_new),rownames(cell_info))
  raw_clust <-cell_info[cells,1]
  
  new_label <- paste(raw_clust,info_new[,1],sep=".")
  cell_info[row_idx,1] <- new_label
  
  # cell_info$subCluster <- factor(cell_info$subCluster,levels = sort(unique(as.numeric(as.character(cell_info$subCluster)))))
  cell_info <-  cell_info[order(as.numeric(as.character(cell_info$subCluster))),,drop=F]
  # cell_info <- droplevels(cell_info)
  return(cell_info)
}


#' @title huber_loss()
#' @export
# Huber Loss Function in R
huber_loss <- function(y_true, y_pred, delta = 1.35) {
  # Calculate the absolute differences
  diff <- abs(y_true - y_pred)
  
  # Initialize a vector to store the Huber loss values
  loss <- numeric(length(diff))
  
  # Calculate loss for each point
  for (i in seq_along(diff)) {
    if (diff[i] <= delta) {
      # Quadratic loss for small differences
      loss[i] <- 0.5 * diff[i]^2
    } else {
      # Linear loss for large differences
      loss[i] <- delta * (diff[i] - 0.5 * delta)
    }
  }
  
  # Return the mean of the loss values
  mean(loss)
}

#' @title segment4df_v1()
#' @description segmentation for each column (cluster) of data.frame
#' @param df data.frame with cluster on column, peak or bin on row
#' @param reference file path for segmentation if acen_remove=TRUE
#' @param segValue.method the method of take the value of segment based on peaks on it.
#' segValue.method is "mean", "median" or "density". For "density", we will take the density$x from the peak of density as the segMean.
#' @param seg.count.lim the maximum number of segments on whole genome. Default is 80.
#' cytoBand <- "/Volumes/B/reference/cytoBand_hg38.tsv"
#' @export
segment4df_v1 <- function(df,cytoBand=NULL,win=0,outpath="./",
                       doFilt=TRUE, 
                       segLen_min = 1e6,
                       segSize_min = 5,
                       outPlot=TRUE,
                       genome="hg38",
                       acen_remove= TRUE,
                       plot_ylim = NULL,#c(0,5),
                       color_dot = TRUE,
                       filename="",
                       seg.method="robseg", #"robseg",#"gaussian"
                       segValue.method="density",
                       seg.count.lim =80,
                       penalty = 1,
                       add_yline=NULL){
  require(ggpubr)
  library(dplyr)
  SegScore_total <- list(seg_score_binLevel=list(),
                         seg_score_segLevel=list(),
                         segPlots=list(),
                         mse=c(),
                         mse.score=c())
  for(k in 1:(ncol(df))){
    colnm <- colnames(df)[k]
    ce_count <- df[,k,drop=F]
    if(win > 0 ){
      ce_count <- smooth_by_chromosome(ce_count,wl=win,method="winMean",win_step=win)
    }

    is.na(ce_count) <- sapply(ce_count, is.infinite)
    is.na(ce_count) <- sapply(ce_count, is.nan)
    ce_count <- na.omit(ce_count)

    if(length(penalty)>1){
      penalty <- round(seq(min(penalty),max(penalty),length.out=5),2)
    }  
    
    mse_res <- c()
    segment_counts <- c()
    segScore_ls <- list()
    for(i in 1:length(penalty)){
      p <- penalty[i]
      seg_bulk<- segment_by_chr_location(ce_count,cytoBand=cytoBand,
                                         acen_remove=acen_remove,
                                         method=seg.method,#"robseg",#"gaussian", #
                                         penalty=p)
      segScore = seg_score(seg_bulk,method=segValue.method)
      #calculte the weighted mean standard error for segmentation
      seg_df <- segScore$seg_score_binLevel[[1]]
      mse_errors <- mse_loss(seg_df,y_true.term = "binRatio",
                              y_pred.term = "SegMean",weight.term ="length_bin")
      #residual sum of squares
      segScore$mse <- mse_errors
      segScore$segment_counts <- dim(segScore$seg_score_segLevel[[1]])[1]
      mse_res <- c(mse_res,mse_errors)
      segment_counts <- c(segment_counts,segScore$segment_counts)
      segScore_ls[[i]] <- segScore

    }
   # Loss Functions: mse
    names(mse_res) <- penalty
    err.score <- mse_res[segment_counts<= seg.count.lim]
    if(length(err.score)>0){
      opt_penalty <- names(err.score)[which.min(err.score)]
      opt_index <- which(names(mse_res) == opt_penalty) 
    }else{
      opt_index <- length(mse_res)
      err.score <- 0
    }


    score <- ifelse(segment_counts > seg.count.lim,mse_res+min(err.score),mse_res)
    segScore <- segScore_ls[[opt_index]]
    segScore$mse.score <- score[opt_index]
      
    segScore$seg_score_segLevel[[1]]$SegSD[is.na(segScore$seg_score_segLevel[[1]]$SegSD)] <- 0
    segScore$seg_score_segLevel[[1]]$SegSE[is.na(segScore$seg_score_segLevel[[1]]$SegSE)] <- 0
    
    df.seg<- segScore$seg_score_segLevel[[1]]
    df.seg<-df.seg[,c("segName","Num_bins")]
    
    seg_df <- segScore$seg_score_binLevel[[1]]
    seg_df <-dplyr::left_join(seg_df,df.seg,by="segName")
    rownames(seg_df) <- seg_df$binID
    
    is.na(seg_df) <- sapply(seg_df, is.infinite)
    is.na(seg_df) <- sapply(seg_df, is.nan)
    #seg_df <- seg_df[,1:(ncol(seg_df)-5)]
    
    if(doFilt){
      seg_df <- seg_df[seg_df$length_seg>=segLen_min,,drop=F]
      seg_df <- seg_df[seg_df$Num_bins>=segSize_min,,drop=F]
      
      if(nrow(seg_df)==0){
        stop(paste0("No segments left on the cutoff of length ",segLen_min))
      }
      

      seg_score_segLevel_filt <- segScore$seg_score_segLevel[[1]][,1:7]
      segnames_keep <- seg_score_segLevel_filt$segName[seg_score_segLevel_filt$Num_bins>=segSize_min &seg_score_segLevel_filt$length_seg >=segLen_min  ]
      seg_score_segLevel_filt <- seg_score_segLevel_filt[seg_score_segLevel_filt$segName%in%segnames_keep,,drop=F]
      SegScore_total$seg_score_segLevel_filt[[colnm]] <- na.omit(seg_score_segLevel_filt)
      
      seg_score_binLevel_filt <- segScore$seg_score_binLevel[[1]]
      seg_score_binLevel_filt <-seg_score_binLevel_filt[seg_score_binLevel_filt$segName%in%segnames_keep,,drop=F]
      SegScore_total$seg_score_binLevel_filt[[colnm]] <- na.omit(seg_score_binLevel_filt)
    }
    
    SegScore_total$mse <- c(SegScore_total$mse,segScore$mse)
    SegScore_total$mse.score <- c(SegScore_total$mse.score,segScore$mse.score)
    seg_df <- na.omit(seg_df)
    if(is.null(plot_ylim)){
      plot_ylim <- round(quantile(seg_df$SegMean,c(0,1)),2)
      plot_ylim <- c(max(0,(plot_ylim[1]-0.5)),plot_ylim[2]+0.5)
    }
    mse_score <- signif(segScore$mse.score,3)
    p1 <- seg_plot(seg_df,name.data=paste0(filename,".",colnm),
                   genome=genome, add_yline=add_yline,
                 color_dot =color_dot,plotDir=outpath,ylab = "Ratio",outPlot=outPlot,ylim = plot_ylim)
    #ggsave(paste0(outpath,"/",filename,".png"),p1,height = 3,width=10,dpi = 300)
    SegScore_total$segPlots[[colnm]] <- p1
    # name.data=paste0("seg_subClone_ratio_",colnm,"(",nrow(seg_df)," peaks)")
    # pp=ggarrange(p1$p1,hist_bulk,p1$p2,  align ="h",ncol = 3, nrow = 1,widths = c(15,3.6,4),heights=3 )
    # ggsave(paste0(outplot,"/",name.data,".png"),pp,height = 3,width=14,dpi = 300)
    SegScore_total$seg_score_binLevel[[colnm]] <- segScore$seg_score_binLevel[[colnm]]
    SegScore_total$seg_score_segLevel[[colnm]] <- segScore$seg_score_segLevel[[colnm]]

    
    rm(seg_df,ce_count)
  }
  return(SegScore_total)
}

#' @title segment4df()
#' @export
#' segmentation per chromosome
segment4df <- function(df,cytoBand=NULL,win=0,outpath="./",
                       doFilt=TRUE, 
                       segLen_min = 1e6,
                       segSize_min = 5,
                       outPlot=TRUE,
                       genome="hg38",
                       acen_remove= TRUE,
                       plot_ylim = NULL,#c(0,5),
                       color_dot = TRUE,
                       filename="",
                       seg.method='PELT',#"strucchange", #"robseg",#"gaussian"
                       segValue.method="median",
                       seg.count.lim =80,
                       penalty = 1,
                       add_yline=NULL,
                       color_seg_gradient = FALSE,
                       rmNA=TRUE){
  require(ggpubr)
  library(dplyr)
  SegScore_total <- list(seg_score_binLevel=list(),
                         seg_score_segLevel=list(),
                         segPlots=list(),
                         mse=c(),
                         mse.score=c())
  chroms <- unique(sapply(strsplit(rownames(df),":|_|-"), "[",1))
  
  for(k in 1:(ncol(df))){
    colnm <- colnames(df)[k]
    ce_count <- df[,k,drop=F]
    if(win > 0 ){
      ce_count <- smooth_by_chromosome(ce_count,wl=win,method="winMean",win_step=win)
    }
    
    is.na(ce_count) <- sapply(ce_count, is.infinite)
    is.na(ce_count) <- sapply(ce_count, is.nan)
    if(rmNA){
      ce_count <- na.omit(ce_count)
    }
   
    
    if(length(penalty)>1){
      if(seg.method == "gaussian"){
        penalty[which.min(penalty)] <- 0.5
        penalty[which.max(penalty)] <- 1
        penalty <- sort(unique(c(1,round(seq(min(penalty),max(penalty),length.out=6),2))))
      }else if(seg.method=="robseg"){
        penalty <-  sort(unique(c(1,round(seq(min(penalty),max(penalty),length.out=6),2))))
      }else{
        penalty=1
      }
    }  
    
    mse_res <- c()
    segment_counts <- c()
    segScore_ls <- list()
    for(p_i in 1:length(penalty)){
      p <- penalty[p_i]
      seg_bulk<- segment_by_chr_location(ce_count,cytoBand=cytoBand,
                                         acen_remove=acen_remove,
                                         method=seg.method,#"robseg",#"gaussian", #
                                         segSize_min=segSize_min,
                                         penalty=p)
      segScore = seg_score(seg_bulk,method=segValue.method)
      #calculte the weighted mean standard error for segmentation
      seg_df_p <- segScore$seg_score_binLevel[[1]]
      mse_errors <- mse_loss(seg_df_p,y_true.term = "binRatio",
                             y_pred.term = "SegMean",weight.term ="length_bin")
      #residual sum of squares
      segScore$mse <- mse_errors
      segScore$segment_counts <- dim(segScore$seg_score_segLevel[[1]])[1]
      mse_res <- c(mse_res,mse_errors)
      segment_counts <- c(segment_counts,segScore$segment_counts)
      segScore_ls[[p_i]] <- segScore
      
    }
   # seg_plot(segScore_ls[[1]]$seg_score_binLevel[[1]],ylab = "Ratio",outPlot=F)$p1
    
    # Loss Functions: mse
    names(mse_res) <- penalty
    err.score <- mse_res[segment_counts<= seg.count.lim]
    if(length(err.score)>0){
      opt_penalty <- names(err.score)[which.min(err.score)]
      opt_index <- which(names(mse_res) == opt_penalty) 
    }else{
      opt_index <- length(mse_res)
      err.score <- 0
    }
    
    score <- ifelse(segment_counts > seg.count.lim,mse_res+min(err.score),mse_res)
    segScore <- segScore_ls[[opt_index]]
    segScore$mse.score <- score[opt_index]
    segScore$seg_score_segLevel[[1]]$SegSD[is.na(segScore$seg_score_segLevel[[1]]$SegSD)] <- 0
    segScore$seg_score_segLevel[[1]]$SegSE[is.na(segScore$seg_score_segLevel[[1]]$SegSE)] <- 0

    df.seg<- segScore$seg_score_segLevel[[1]]
    df.seg<-df.seg[,c("segName","Num_bins")]
    
    seg_df <- segScore$seg_score_binLevel[[1]]
    seg_df <-dplyr::left_join(seg_df,df.seg,by="segName")
    rownames(seg_df) <- seg_df$binID
    
    is.na(seg_df) <- sapply(seg_df, is.infinite)
    is.na(seg_df) <- sapply(seg_df, is.nan)
    #seg_df <- seg_df[,1:(ncol(seg_df)-5)]
    
    if(doFilt){
      seg_df <- seg_df[seg_df$length_seg>=segLen_min,,drop=F]
      seg_df <- seg_df[seg_df$Num_bins>=segSize_min,,drop=F]
      
      if(nrow(seg_df)==0){
        stop(paste0("No segments left on the cutoff of length ",segLen_min))
      }
      
      
      seg_score_segLevel_filt <- segScore$seg_score_segLevel[[1]][,1:7]
      segnames_keep <- seg_score_segLevel_filt$segName[seg_score_segLevel_filt$Num_bins>=segSize_min &seg_score_segLevel_filt$length_seg >=segLen_min  ]
      seg_score_segLevel_filt <- seg_score_segLevel_filt[seg_score_segLevel_filt$segName%in%segnames_keep,,drop=F]
      
      seg_score_binLevel_filt <- segScore$seg_score_binLevel[[1]]
      seg_score_binLevel_filt <-seg_score_binLevel_filt[seg_score_binLevel_filt$segName%in%segnames_keep,,drop=F]
      
      if(rmNA){
        SegScore_total$seg_score_segLevel_filt[[colnm]] <- na.omit(seg_score_segLevel_filt)
        SegScore_total$seg_score_binLevel_filt[[colnm]] <- na.omit(seg_score_binLevel_filt)
      }
      
      
         
    }
    
    SegScore_total$mse <- c(SegScore_total$mse,segScore$mse)
    SegScore_total$mse.score <- c(SegScore_total$mse.score,segScore$mse.score)
    seg_df <- seg_df[,!grepl("seg_SD|seg_score",colnames(seg_df))]
    if(rmNA){
     seg_df <- na.omit(seg_df)
    }
    if(is.null(plot_ylim)){
      plot_ylim <- round(quantile(seg_df$SegMean,c(0,1),na.rm=T),2)
      plot_ylim <- c(max(0,(plot_ylim[1]-0.5)),plot_ylim[2]+0.5)
    }
    mse_score <- signif(segScore$mse.score,3)
    # if(outPlot){
      p1 <- seg_plot(na.omit(seg_df),name.data=paste0(filename,".",colnm),
                     genome=genome, add_yline=add_yline,
                     color_dot =color_dot,plotDir=outpath,ylab = "Ratio",outPlot=outPlot,ylim = plot_ylim,color_seg_gradient=color_seg_gradient)
      #ggsave(paste0(outpath,"/",filename,".png"),p1,height = 3,width=10,dpi = 300)
    # }else{
    #   p1 <- NULL
    # }

    SegScore_total$segPlots[[colnm]] <- p1
    # name.data=paste0("seg_subClone_ratio_",colnm,"(",nrow(seg_df)," peaks)")
    # pp=ggarrange(p1$p1,hist_bulk,p1$p2,  align ="h",ncol = 3, nrow = 1,widths = c(15,3.6,4),heights=3 )
    # ggsave(paste0(outplot,"/",name.data,".png"),pp,height = 3,width=14,dpi = 300)
    SegScore_total$seg_score_binLevel[[colnm]] <- segScore$seg_score_binLevel[[colnm]]
    SegScore_total$seg_score_segLevel[[colnm]] <- segScore$seg_score_segLevel[[colnm]]
    
    
    rm(seg_df,ce_count)
  }
  return(SegScore_total)
}

#' @title seg4clusters()
#' @description segmentation for clusters
#' @param mat single-cell matrix with cells on columns for generate cluster-level pseudo-bulk matrix
#' @param cell_info data.frame with cells as rownames, and clusters on the first column
#' @param relative_method  ratio or Subtraction for relative value
#' @export
seg4clusters <- function(mat,cell_info,
                         cytoBand=NULL,
                         col_ref=NULL,
                         method ="mean",
                         adjust_zero=TRUE,
                         seg.penalty = 1,
                         seg.method="PELT",#"robseg", #"gaussian",
                         segSize.min = 10,
                         doFilt=TRUE,
                         outplot_path="./",
                         outname=NULL,
                         outFigure=TRUE,
                         color_seg_gradient=TRUE,
                         relative_method="ratio",
                         plot.ylim=NULL){
  if(!file.exists(outplot_path) & outFigure){dir.create(outplot_path,recursive=T)}
  if(!is.null(cytoBand)){
    if (Reduce("|", is(cytoBand) == "character")) {
      if(substr(cytoBand, nchar(cytoBand)-2, nchar(cytoBand)) %in% c(".gz","txt","tsv","csv")){
        cytoBand <- read.table(connection <- gzfile(cytoBand, 'rt'), sep="\t", header=F, check.names=FALSE)
        close(connection)
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
  cell_info_new <- data.frame(cellid=rownames(cell_info),group=as.character(cell_info[,1]))
  mtx_qcs_bulk <- pseudo_bulk_v2(mat[,cell_info_new$cellid],group_ano = cell_info_new,method =method,adjust_zero=adjust_zero)
  if(is.null(col_ref)){
    clt.sd <- colSds(as.matrix(mtx_qcs_bulk),na.rm =T)
    col_ref <- colnames(mtx_qcs_bulk)[which.min(clt.sd)]
    print(col_ref)
  }
  if(is.null(outname)){
    outname <- paste0("segRatio_subClone_vs",col_ref,".pdf")
  }
  
  col_obs <- colnames(mtx_qcs_bulk)[!colnames(mtx_qcs_bulk) %in% col_ref]
  if(relative_method=="ratio"){
    mtx_qcR_bulk0 <- as.matrix(mtx_qcs_bulk[,col_obs]/mtx_qcs_bulk[,col_ref]) 
    
  }else{
    mtx_qcR_bulk0 <- as.matrix(mtx_qcs_bulk[,col_obs]-mtx_qcs_bulk[,col_ref]) 
  }
   
  colnames(mtx_qcR_bulk0) <- col_obs
  rownames(mtx_qcR_bulk0)<- rownames(mtx_qcs_bulk)
  is.na(mtx_qcR_bulk0) <- sapply(mtx_qcR_bulk0, is.infinite)
  is.na(mtx_qcR_bulk0) <- sapply(mtx_qcR_bulk0, is.nan)
  mtx_qcR_bulk0 <- na.omit(mtx_qcR_bulk0)
  
  SegScore_total <- suppressMessages(segment4df(mtx_qcR_bulk0,   #mtx_qcs_bulk,
                               cytoBand = cytoBand,
                               outpath = outplot_path,
                               outPlot = F,
                               filename = paste0("segRatio_Clust_vs",col_ref),
                               seg.method=seg.method,
                               segSize_min = segSize.min,
                               penalty = seg.penalty,
                               color_dot = F,
                               doFilt=doFilt,
                               color_seg_gradient = color_seg_gradient,
                               plot_ylim = plot.ylim))
  if(outFigure){
    p_ls <- lapply(SegScore_total$segPlots,function(x){y=x$ggarranged_p;return(y)})
    pcom <- cowplot::plot_grid(plotlist=p_ls,ncol = 1,align="v")
    ggsave(paste0(outplot_path, "/",outname),pcom, width=14, height=2*ncol(mtx_qcR_bulk0),device = pdf,bg="white")
    
  }
  return(SegScore_total)
}

#' @title seg_assess()
#' @param seg_ls list output from seg4clusters()
#' @export
seg_assess <- function(seg_ls,sd.cutoff=NULL,outname=paste0("./hist.seg_segRatio.pdf"),outplot=TRUE){

  #filterout bad subCluster
  ls_new= seg_ls$seg_score_segLevel
  segSD <- do.call(rbind,lapply(1:length(seg_ls$seg_score_segLevel),function(x,ls_new){
    df <- seg_ls$seg_score_segLevel[[x]]
    df$subCluster=names(seg_ls$seg_score_segLevel)[x]
    return(df)
  },ls_new))
  
  if(is.null(sd.cutoff)){
    sd.cutoff=0.1
  }else if(is.data.frame(sd.cutoff)){
    segSD <- segSD %>%
      left_join(sd.cutoff, by = "subCluster") %>%
      as.data.frame()
  }
  
  #filter out subCluster where mean SegSD >0.1
  segSD <- segSD %>%
    dplyr::group_by(subCluster)%>%
    dplyr:: mutate(meanSD=round(mean(SegSD),2),SegMean_r=round(SegMean,1)) %>%
    as.data.frame()
  
  segSD_filt <- segSD %>%
    dplyr::group_by(subCluster)%>%
    dplyr::filter(meanSD<=sd.cutoff) %>%
    as.data.frame()
  
  clt_rm<- which(!names(seg_ls$seg_score_binLevel) %in% segSD_filt$subCluster)
  cltName_rm <- names(seg_ls$seg_score_binLevel)[clt_rm]
  if(length(clt_rm)>0){
    SegScore_total_new <- list(seg_score_binLevel=seg_ls$seg_score_binLevel[-clt_rm],
                               seg_score_segLevel=seg_ls$seg_score_segLevel[-clt_rm],
                               segPlots=seg_ls$segPlots[-clt_rm],
                               subCluster=unique(segSD_filt$subCluster))
  }else{
    SegScore_total_new <- list(seg_score_binLevel=seg_ls$seg_score_binLevel,
                               seg_score_segLevel=seg_ls$seg_score_segLevel,
                               segPlots=seg_ls$segPlots,
                               subCluster=unique(segSD_filt$subCluster))
  }
  if(length(SegScore_total_new$seg_score_binLevel)>1){
    df_seg_combin <- do.call(rbind,lapply(SegScore_total_new$seg_score_segLevel,function(x)x))
   # ylim = c(0.5,1.5)
    ylim = c(min(df_seg_combin$SegMean)*0.9,max(df_seg_combin$SegMean)*1.1)
    hist_bulk <- ggplot(df_seg_combin,aes(x=SegMean))+
      geom_histogram(aes(y=..density..,weight =length_seg),bins=50,alpha=0.8, position = 'identity',colour = "#e9ecef", fill = "#404080")+ #binwidth = 0.1,
      theme_minimal()+
      labs(fill="",x="",y="",title = paste0("i vs ",col_ref))+
      xlim(ylim)+
      coord_flip(xlim=ylim)+
      geom_density()+
      theme(
        plot.title = element_text(size=10,angle=0,vjust = 0,margin = margin(0,0,0,0))
      )
    histdata <- ggplot_build(hist_bulk)$data[[1]]
    highest_bin.x <- signif(histdata$x[which.max(histdata$count)],3)
    highest_bin.y <- signif(histdata$y[which.max(histdata$count)],3)
    
    # highest_bin = densMode(df_seg_combin$SegMean)
    # highest_bin.x = signif(highest_bin$x,3)
    # highest_bin.y = signif(highest_bin$y,3)
    
    hist_bulk <- hist_bulk +
      geom_vline(xintercept = highest_bin.x, linetype = "dashed", color = "darkgrey")+
      annotate("text", x = highest_bin.x, y = highest_bin.y*0.9, label = round(highest_bin.x,2),hjust = 0, vjust = -0.5)
    if(outplot){
      ggsave(outname,hist_bulk,width = 3,height = 3.3)
    }
   
    #find clusters where ratio distance between segments smaller than primary segments
    primary_seg <- segSD_filt[segSD_filt$SegMean_r %in% round(highest_bin.x,1),]
    primary_seg <- primary_seg %>% 
      dplyr::group_by(SegMean_r) %>%
      dplyr::summarise(length_seg_mean=mean(length_seg)) %>%
      as.data.frame()
    highest_bin.x_primary <- primary_seg$SegMean_r[which.max(primary_seg$length_seg_mean)]
    
    dist_seg <-  abs(round(highest_bin.x[1]-highest_bin.x[2],1))
    segSD_filt2 <- segSD_filt %>%
      dplyr::group_by(subCluster)%>%
      dplyr::filter(any(SegMean_r %in% highest_bin.x_primary )) %>%
      as.data.frame()
    
    clt_filt <-  segSD_filt %>%
      dplyr::group_by(subCluster)%>%
      dplyr::filter(any(SegMean_r %in% highest_bin.x_primary )) %>%
      dplyr::summarize(min_diff = abs(min(diff(SegMean_r),na.rm = TRUE)))%>%
      dplyr::filter(min_diff>=dist_seg) %>%
      as.data.frame()
    
    clts_reClust <- unique(segSD_filt$subCluster)[!unique(segSD_filt$subCluster) %in% c(clt_filt$subCluster,col_ref)]       
    
    return(list(seg.df=segSD,
                segMean.primary=highest_bin.x_primary,
                cltName_rm= cltName_rm,
                clusters.reclust=clts_reClust,
                dist=dist_seg))
  
  }else{return(list(seg.df=segSD,
                    segMean.primary=NULL,
                    cltName_rm= cltName_rm,
                    clusters.reclust=NULL,
                    dist=NULL)
               )}
  
}




#‘ @description add columns of chrosome absolute start, end and segLen for "segName" of a data.frame
#' @title add_segLen()
#' @export
add_segLen <- function(df,genome="hg38",true_loc=TRUE){
  suppressMessages({
    library(BSgenome)
    library(dplyr)})
  g <- getBSgenome(genome, masked=FALSE)
  chrLine <- data.frame(chrom=1:24, length=seqlengths(g)[1:24])
  chromnum=sapply(strsplit(rownames(df),'-|_|:'),'[',1)
  
  if(grepl('chr',chromnum[1])){
    chromnum=sapply(strsplit(chromnum,'hr'),'[',2) ## if chrommnum with chr
  }else{
    chromnum=chromnum
  }
  chromnum <- gsub("X","23",chromnum)
  chromnum <- gsub("Y","24",chromnum)
  df$chrom <- chromnum
  chrLine <- chrLine[chrLine$chrom %in% as.character(unique(chromnum)),,drop=F] 
  genomeLength <- sum(chrLine)
  if(nrow(chrLine)>1){
    chrPosStart = c(0,chrLine[1:(nrow(chrLine)-1),2])
  }else(
    chrPosStart = 0
  )
  chrLine$chrPosStart_abs <- cumsum(chrPosStart)
  df <- merge(df,chrLine,by="chrom",all.x=T)
  
 
  df$segStrat <- as.numeric(sapply(strsplit(df$segName,'_'),'[',2))+as.numeric(df$chrPosStart_abs)
  df$segEnd <- as.numeric(sapply(strsplit(df$segName,'_'),'[',3))+as.numeric(df$chrPosStart_abs)
  df <- df[order(as.integer(df$chrom),as.numeric(df$chrPosStart_abs),as.numeric(df$segStrat)),]
  true_loc = TRUE
  if(true_loc){
    df$segLen <- as.numeric(df$segEnd)-as.numeric(df$segStrat)
  }else{
    df$segLen <- as.numeric(df$length_seg)
  }   
  segLen_sum <- sum(unique(df[,c("segID","segLen")]$segLen))
  df$wei <- df$segLen/segLen_sum
  
  return(df)
  
}

#' @title .get_ref_mean()
#' @keywords internal function
#' @param matrix input matrix, group in columns (cells) 
#' @param ref_indice numeric indices for group
#' @export
.get_ref_mean <- function(matrix,ref_indice,zero.rm=FALSE){
  library(matrixStats)
  RN <- rownames(matrix)
  matrix_data <- matrix[,ref_indice,drop=FALSE]
  if(zero.rm){
    peak_ref_grp_mean <- apply(matrix_data, 1, function(x){
      #the condition mean of Zero-Inflated distribution
      x=x[!is.na(x)]
      if(all(x==0)){
        adjusted_mean=0
      }else{
        non_zero_mean <- mean(x[x != 0],na.rm=TRUE)
        zero_proportion <- sum(as.numeric(x) == 0) / length(x)
        adjusted_mean <- non_zero_mean *(1 - zero_proportion)
      }
      return(adjusted_mean)
    } )
  }else{
    peak_ref_grp_mean <- rowMeans(matrix_data)
  }
  
  peak_ref_grp_sd <- apply(matrix_data, 1, sd)
  peak_ref_grp_sum <- rowSums(matrix_data)
  peak_ref_grp_means_sd <- data.frame(mean=peak_ref_grp_mean,SD=peak_ref_grp_sd,sum=peak_ref_grp_sum)
  rownames(peak_ref_grp_means_sd) <- RN
  return(peak_ref_grp_means_sd)
}

#' @title find_mean()
#' @description Sort the cells according to the expression level from low to high, 
#' increase the number of cells from the first one until the median value of the 
#' row mean of the selected cells reaches the threshold, and return the selected cells
#' @param mt matrix with cells in column and features in row
#' @param order_col order columns from low to high (l2h) or high to low (h2l)
#' @export
find_mean <- function(mt,cutoff,order_col="l2h",zero.rm=TRUE){
  if(order_col %in% "l2h"){
  mt <- mt[,order(colSums(mt))]
  }else{
    mt <- mt[,order(colSums(mt),decreasing = TRUE)]
  }
  median_count <- c()
  for(i in 1:ncol(mt)){
    mti <- mt[,1:i,drop=F]
    mean_mti <- .get_ref_mean(mti,1:i,zero.rm=zero.rm)
    mean_mti_median <- median(mean_mti$mean)
    median_count <- c(median_count,mean_mti_median)
  ##
  if(i >= 10){
    if(order_col%in%"l2h"){
      if(all(median_count[(i-2):i] >= cutoff)){
        break
      }
    }else{
      if(all(median_count[(i-2):i] <= cutoff)){
        break
      }
    }
  }
  
  }
  if(i==ncol(mt)){
    return(list(meanCount=mean_mti_median,cols_within=colnames(mt)))
  }else{
    return(list(meanCount=median_count[(i-2)],cols_within=colnames(mt)[1:(i-2)]))
  }
}

#' @title filt_peaks_bySum()
#' @param mat matrix with peaks on row, cells on clumn
#' @export
filt_peaks_bySum <- function(mat,rows,zscore.lim=2){
  totalcount <- rowSums(mat)
  zscore <- scale(totalcount)
  peaks <- names(totalcount)[abs(zscore) < zscore.lim]
  return(peaks)
}


#' @title findRefClone()
#' @description find Reference Clone based On Clonal Segmentation for Ratio
#' @param mat input matrix with cell barcodes on columns
#' @param cell_metadata data.frame with cell barcodes as rownames and cell annotation on the fist column
#' @export
#' 
findRefClone <- function(mat,cell_metadata,seg.method="gaussian",seg.penalty=0.4,
                         outplot_path="./"
                         ){
  col_ref.sd <- c()
  col_ref.global <- c()
  primary_segRatio.global <- c()
  SegScore_total.ls <- list()
  groups <- unique(as.character(cell_metadata[,1]))
  for(col_ref in groups){
    SegScore_total <- seg4clusters(mtx_chr[,rownames(cell_metadata)],cell_metadata,col_ref,
                                   doFilt=T,
                                   seg.penalty = seg.penalty,
                                   seg.method=seg.method,
                                   outplot_path=outplot_path,
                                   outname=paste0("segRatio_groupi_vs",col_ref,".pdf"))
    SegScore_total.ls[[as.character(col_ref)]] <- SegScore_total
    seg_df <- do.call(rbind,lapply(SegScore_total$seg_score_segLevel,function(x)x))
    p = plot_hist(seg_df,x="SegMean",wei="length_seg",title=paste0("segRatio_vs",col_ref),
                  ylim = c(min(seg_df$SegMean),max(seg_df$SegMean)),
                  label_primary=T)
    histdata <- ggplot_build(p)$data[[1]]
    histdata <- na.omit(histdata)
    count_max <- max(histdata$count)
    highest_bin.x <- signif(median(histdata$x[histdata$count == count_max]),2)
    primary_segRatio <- highest_bin.x[1]
    
    primary_segRatio.global <- c(primary_segRatio.global,primary_segRatio)
    
    if((!is.na(primary_segRatio)) & primary_segRatio==1){
      col_ref.global <- cbind(col_ref.global,col_ref)
      ls_new= SegScore_total$seg_score_segLevel
      segSD <- do.call(rbind,lapply(1:length(SegScore_total$seg_score_segLevel),function(x,ls_new){
        df <- SegScore_total$seg_score_segLevel[[x]]
        df$subCluster=names(SegScore_total$seg_score_segLevel)[x]
        return(df)
      },ls_new))
      meanSD <- median(segSD$SegSD,na.rm=TRUE)
      col_ref.sd <- cbind(col_ref.sd,meanSD)
    }
    rm(SegScore_total)
  }
    names(primary_segRatio.global) <- groups
    col_ref<- col_ref.global[which.min(col_ref.sd)]
    if(is.null(col_ref)){
      dist2one <- abs(1-primary_segRatio.global)
      col_ref <- names(dist2one)[which.min(dist2one)]
    }
    SegScore_total <- SegScore_total.ls[[col_ref]]
    
    return(list(refGroup=col_ref,
                SegScore_total=SegScore_total,
                SegScore_total.ls=SegScore_total.ls))
  
}


#' @title findClonePerChr()
#' @description sub-clustering cells per chromosome count matrix if the method is "perChr"
#' @param mat peak-cell count matrix of whole genome. peaks as rownames with "chrx-xxx-xxx", "chrx:xxx-xxx" or "chrx_xxx_xxx" format
#' @param ref.value 
#' @param ChrExclude the list of chromosome should be excluded from analysis. Default = c('chrX', 'chrY', 'chrM')
#' @export
findClonePerChr <- function(mat,
                            ref.value,
                            outdir=".",
                            ChrExclude=c('chrX', 'chrY', 'chrM'),
                            clust_resolution=1
                            ){
  suppressMessages({
    library(futile.logger)
  })
  if(outdir != "." && !file.exists(outdir)){dir.create(outdir,recursive=T)}
  if (Reduce("|", is(ref.value) == "character")) {
    flog.info(sprintf("Loading segment ratio file: %s", ref.value))
    if(substr(ref.value, nchar(ref.value)-3, nchar(ref.value)) == ".rds") {
      ref.value <- readRDS(ref.value)
    }
  }else if(Reduce("|", is(ref.value) == "data.frame")){
    ref.value <- ref.value
  }else{
    stop("Can not find ref.value as a data.frame or filename of '.rds'.")
  }
  
  
  chrs <- unique(sapply(strsplit(rownames(mat),":|-|_"),"[",1))
  chrs <- chrs[grepl("chr",chrs) &(!grepl(paste(ChrExclude,collapse = "|"),chrs))]
  flog.info(paste0("Starting inferring sub-clonal CNV ratio for chromosome matrix."))
  for(chrj in chrs){
    flog.info(paste0("Analysing ",chrj))
    outdir_j <- paste0(outdir,"/",chrj);if(!file.exists(outdir_j)){dir.create(outdir_j,recursive=T)}
    mtx_chr <- mat[grepl(paste0("^",chrj,"-"),rownames(mat)),]
    ref.value_chr <- ref.value[rownames(mtx_chr),]
    #clustering
    cell_info_ls <- suppressMessages(find_subclust(mtx_chr,resolution = clust_resolution,sep_by = c("-", "-")))
    cell_info <- cell_info_ls$meta_cluster
    write.csv(cell_info,paste0(outdir_j,"/cell_info_by",chrj,".csv"),row.names = T)
    color_r <- suppressWarnings(get_group_color_palette("Set3")(length(unique(cell_info$subCluster))))
    names(color_r) <- levels(cell_info$subCluster)
    left_anno_cols[["subCluster"]] <- color_r
    fileout_name = paste0("heat_nomCount_subClust_",chrj)
    ht <- heatmap4peakMt(mtx_chr[,rownames(cell_info)],meta_info=cell_info[,c("subCluster"),drop=F],
                         value.type = "count",legend_titles="peak count", width=3,height=5,
                         outdir=outdir_j,
                         fileout_name=fileout_name,
                         col_list=left_anno_cols,
                         clust_rows=F,clustering_method_rows="ward.D2",sep_by="-")
    ###4. segment for sub-cluster using   segment4df() on Ratio
    cell_ano2_all <- data.frame(cellid=rownames(cell_info),group=as.character(cell_info[,"subCluster"]))
    mtx_chr_bulk <- pseudo_bulk_v2(mtx_chr[,cell_ano2_all$cellid],group_ano = cell_ano2_all,method ="mean",adjust_zero=TRUE)
    #segmentation on ratio
    col_obs <- colnames(mtx_chr_bulk)
    mtx_chr_bulk0 <- as.matrix(mtx_chr_bulk[,col_obs]/ref.value_chr) 
    is.na(mtx_chr_bulk0) <- sapply(mtx_chr_bulk0, is.infinite)
    is.na(mtx_chr_bulk0) <- sapply(mtx_chr_bulk0, is.nan)
    mtx_chr_bulk0 <- na.omit(mtx_chr_bulk0)
   # any(is.na(mtx_chr_bulk0)) 
    ylim = c(min(mtx_chr_bulk0),max(mtx_chr_bulk0))
    SegScore_total <- segment4df(mtx_chr_bulk0,  
                                 cytoBand = cytoBand,
                                 outpath = outdir_j,
                                 outPlot = F,
                                 filename = "segRatio_subClust")
    p_ls <- lapply(SegScore_total$segPlots,function(x){y=x$ggarranged_p;return(y)})
  }
}

#' @title findClone.SortCells()
#' @export
findClone.SortCells <- function(mtx_sub,
                                df_seg_baseline,
                                ref.value,
                                cytoBand=NULL,
                                Ncs = 20,
                                er = 0.005,
                                outdir="."){
  if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
  df_seg <- df_seg_baseline[!is.na(df_seg_baseline$ratio_map),]
  initialCN <- as.data.frame(unique(cbind(df_seg$ratio_map,df_seg$integerCNV)))
  colnames(initialCN) <- c("ratio","CN")
  initialCN <- initialCN[order(initialCN$CN),]
  initialCN$base=0
  initialCN1 <- peakIndex(df_seg,initialCN)
  
  CNest <- data.frame(initialCN1[,c(1,3,2)])
  cloneDelt <- (CNest$ratio[2]-CNest$ratio[1])/(CNest[2,3]-CNest[1,3])
  ratio_interval <- CNest$ratio
  min_level <- min(ratio_interval)-cloneDelt
  ratio.min <- ifelse(min_level>=0,min_level,ifelse(min_level>-0.1 &min_level<0.1, 0,min(ratio_interval)))
  ratio.max <- max(ratio_interval)+cloneDelt
  ratio_interval <- rev(c(ratio.min,ratio_interval,ratio.max)) # from high to low
  ratio_interval <- unique(round(ratio_interval,2))
  cloneBase <- round(CNest$ratio[CNest$base==1],2)
  ###find sub-clone by sorting cells to ratio
  df_seg_clean <- unique(df_seg[,c("segName","chrom","Chromosome","SegMean","segStrat","segEnd","ratio_map","integerCNV","segLen")])
  df_seg_clean <- df_seg_clean[order(as.numeric(df_seg_clean$chrom),df_seg_clean$segStrat),]
  
  for(k in 1:nrow(df_seg_clean)){
    print(k)
    segk <-df_seg_clean$segName[k]
    
    peaks <- df_seg$binID[df_seg$segName == segk]
    peaks <- intersect(peaks, rownames(mtx_sub))
    
    if(length(peaks)<10){
      next
    }else{
      mtx_seg <- mtx_sub[peaks,,drop=F]
      mtx_seg <- mtx_seg[,order(colSums(mtx_seg),decreasing = T),drop=F]
      cell_info <- data.frame(row.names =colnames(mtx_seg),cellID = colnames(mtx_seg),subCluster=NA,score=NA)
      outdir_j <- paste0(outdir,"/",segk);if(!file.exists(outdir_j)){dir.create(outdir_j,recursive=T)}
      setwd(outdir_j)
      
      #group cells from high to low
      col.start <- 1
      col.end <- col.start + Ncs - 1
      
      while (col.start < ncol(mtx_seg)) {
        #print(col.start)
        ratio_ls <- mt2ratio_by_col(mtx_seg, col_idx = c(col.start,col.end),ref.value = ref.value[peaks,])
        ratio <- ratio_ls$segRatio
        #where is ratio in ratio_interval
        ratio_idx <- which((ratio_interval-ratio)<=0)[1]-1
        dist2level <- abs(ratio_interval-ratio)
        nearest_level <- ratio_interval[which.min(dist2level)]
        cell_score <- dist2level[which.min(dist2level)]
        if(!is.na(ratio_idx)){
          if(ratio_idx != 0 ){
            level_up <- ratio_interval[ratio_idx]
          }else{level_up <- ratio_interval[1]}
          level_dn <- ratio_interval[ratio_idx + 1]
        }else{level_up = level_dn <- min(ratio_interval) }
        
        if((col.end-col.start)==( Ncs - 1)){
          cell_info$score[col.start:col.end] <- cell_score
        }else{
          cell_info$score[col.end] <- cell_score
        }
        #If the ratio is down to cloneBase, then switch the finding direction from low to high.
        if(ratio > (cloneBase - er) & ratio < (cloneBase + er) ){
          break
        }else{
          cell_used <- cell_info$cellID[col.start:col.end]
          if(ratio > (level_up - er) & ratio < (level_up + er)){
            cell_info[cell_used,"subCluster"] <- level_up
            
            col.start <- col.end + 1
            if(col.start>= ncol(mtx_seg)){break}
            col.end <- col.start + Ncs - 1
            if(col.end >= ncol(mtx_seg)){col.end <- ncol(mtx_seg)}
            
          }else if(ratio > (level_dn - er) & ratio < (level_dn + er)){
            cell_info[cell_used,"subCluster"] <- level_dn
            
            col.start <- col.end + 1
            if(col.start>= ncol(mtx_seg)){break}
            col.end <- col.start + Ncs - 1
            if(col.end >= ncol(mtx_seg)){col.end <- ncol(mtx_seg)}
            
          }else{
            col.end <- col.end +1 
            if(col.end>=ncol(mtx_seg)){break}
          }
        }
      }
      #Find ratio from low to high
      if(ratio > (cloneBase - er) & ratio < (cloneBase + er)){
        cells_keep <- rownames(cell_info)[is.na(cell_info$subCluster)]
        mtx_seg <- mtx_seg[,cells_keep,drop=F]
        mtx_seg <- mtx_seg[,order(colSums(mtx_seg),decreasing = F),drop=F]
        if(ncol(mtx_seg)<=Ncs){
          ratio_ls <- mt2ratio_by_col(mtx_seg,ref.value = ref.value[peaks,])
          ratio <- ratio_ls$segRatio
          ratio_idx <- which((ratio_interval-ratio)<=0)[1]-1
          dist2level <- abs(ratio_interval-ratio)
          nearest_level <- ratio_interval[which.min(dist2level)]
          cell_score <- dist2level[which.min(dist2level)]
          cell_info$score[is.na(cell_info$subCluster)] <- cell_score
          cell_info$subCluster[is.na(cell_info$subCluster)] <- nearest_level
        }else{
          col.start <- 1
          col.end <- col.start + Ncs - 1
          while (col.start < ncol(mtx_seg)) {
            #print(col.start)
            ratio_ls <- mt2ratio_by_col(mtx_seg, col_idx = c(col.start,col.end),ref.value = ref.value[peaks,])
            ratio <- ratio_ls$segRatio
            
            #where is ratio in ratio_interval
            ratio_idx <- which((ratio_interval-ratio)<=0)[1]-1
            dist2level <- abs(ratio_interval-ratio)
            nearest_level <- ratio_interval[which.min(dist2level)]
            cell_score <- dist2level[which.min(dist2level)]
            
            if(!is.na(ratio_idx)){
              if(ratio_idx != 0 ){
                level_up <- ratio_interval[ratio_idx]
              }else{level_up <- ratio_interval[1]}
              level_dn <- ratio_interval[ratio_idx + 1]
            }else{level_up = level_dn <- min(ratio_interval) }
            
            if((col.end-col.start)==( Ncs - 1)){
              
              cell_info$score[cell_info$cellID %in% colnames(mtx_seg)[col.start:col.end]] <- cell_score
            }else{
              cell_info$score[cell_info$cellID %in% colnames(mtx_seg)[col.end]] <- cell_score
            }
            #If the ratio is down to cloneBase, then switch the finding direction from low to high.
            if(ratio > (cloneBase - er) & ratio < (cloneBase + er) ){
              break
            }else{
              cell_used <- colnames(mtx_seg)[col.start:col.end]
              if(ratio > (level_up - er) & ratio < (level_up + er)){
                cell_info[cell_used,"subCluster"] <- level_up
                
                col.start <- col.end + 1
                if(col.start>= ncol(mtx_seg)){break}
                col.end <- col.start + Ncs - 1
                if(col.end >= ncol(mtx_seg)){col.end <- ncol(mtx_seg)}
                
              }else if(ratio > (level_dn - er) & ratio < (level_dn + er)){
                cell_info[cell_used,"subCluster"] <- level_dn
                
                col.start <- col.end + 1
                if(col.start>= ncol(mtx_seg)){break}
                col.end <- col.start + Ncs - 1
                if(col.end >= ncol(mtx_seg)){col.end <- ncol(mtx_seg)}
                
              }else{
                col.end <- col.end +1 
                if(col.end>=ncol(mtx_seg)){break}
              }
            }
          }
        }
        
        if(any(is.na(cell_info$subCluster))){
          cells_left <- cell_info$cellID[is.na(cell_info$subCluster)]
          if(length(cells_left)>10){
            ratio_ls <- mt2ratio_by_col(mtx_seg[,cells_left,drop=F],ref.value = ref.value[peaks,])
            ratio <- ratio_ls$segRatio
            dist2level <- abs(ratio_interval-ratio)
            nearest_level <- ratio_interval[which.min(dist2level)]
            cell_score <- dist2level[which.min(dist2level)]
            cell_info$subCluster[is.na(cell_info$subCluster)] <- nearest_level
            
          }else{
            cell_info$subCluster[is.na(cell_info$subCluster)] <-cloneBase
          }
        }
        any(is.na(cell_info$score))
        
      }else if(all(is.na(cell_info$subCluster))){
        # cell_info$subCluster[is.na(cell_info$subCluster)] <- ratio_interval[which.min(abs(ratio_interval-ratio))]
      }else{
        cell_info$score[is.na(cell_info$subCluster)] <- 0
        cell_info$subCluster[is.na(cell_info$subCluster)] <- nearest_level
      }
      
      table(cell_info$subCluster,useNA="ifany")
      write.csv(cell_info,paste0(outdir_j,"/cell_info_by",segk,".csv"),row.names = T)
      
      ###4. segment for sub-cluster using   segment4df()
      if(length(unique(cell_info$subCluster))>1){
        cell_info2 <- cell_info[,"subCluster",drop=F]

        mt_ana <- mtx_sub[,rownames(cell_info2)] # calculating sub-clone ratio using the same matrix with bulk ratio
        
        cell_info2b <- data.frame(rownames(cell_info2),cell_info2)
        ratio_res_pl <- CNratio_peak(mt_ana,cell_info2b,ref.value=ref.value[rownames(mt_ana),,drop=F])
        
        outplot_path <- paste0(outdir_j,"/subCluster_ratio");if(!file.exists(outplot_path)){dir.create(outplot_path,recursive=T)}
        SegScore_total <- segment4df(ratio_res_pl,
                                     cytoBand = cytoBand,
                                     doFilt=F,
                                     outpath = outplot_path,
                                     filename = "segRatio_finalsubClone")
        #Process result
        res_ls <- list()
        for(cl in names(SegScore_total$seg_score_binLevel)){
          df_seg_C1 <- SegScore_total$seg_score_binLevel[[cl]]
          df_seg_C1_add <- add_segLen(df_seg_C1)
          df_seg_C1_add <- df_seg_C1_add[order(as.integer(df_seg_C1_add$chrom),as.numeric(df_seg_C1_add$chrPosStart_abs),as.numeric(df_seg_C1_add$Start)),]
          res_ls[[cl]] <- df_seg_C1_add
        }
        saveRDS(res_ls,paste0(outdir_j,"/segRatio_",segk,".rds"))
        rm(mtx_seg)
      
    }
    }
  }
}

#' @title compare2RefSeg()
#' @description the proportion of peaks (or bins) with different integerCNV  between df_query and df_ref
#' @param df_ref data.frame with c("binID","ratio_map","integerCNV") columns
#' @param df_query data.frame with "binID", "SegMean" columns
#' @export

compare2RefSeg <- function(df_ref,df_query,df_binRatio=df_query,delta_base=NULL,diff_prop_lim=0.01,reCalcuCNA=FALSE){
  CNbase_bin <- df_ref[,c("binID","ratio_map","integerCNV")]
  CNbase <- unique(df_ref[,c("ratio_map","integerCNV")])
  CNbase <- peakIndex(df_ref,CNbase)
  CNbase <- CNbaseline.fill(CNbase,delt=delta_base)
  # if(is.null(delta_base)){
  #   delta_base <- (CNbase$ratio[2]-CNbase$ratio[1])/(CNbase$CN[2]-CNbase$CN[1])
  # }
  
  bin_indx <- match(df_query$binID,CNbase_bin$binID)
  if(all(is.na(bin_indx))){
    stop("'binID' does NOT match between df_ref and df_query.")
  }
  df_query <- cbind(df_query,CNbase_bin[bin_indx,2:3])
  df_query$integerCNV_ref <- df_query$integerCNV 
  df_query$ratio_map_ref <- df_query$ratio_map
  #fill NA with neighbor
  df_seg_vsRef.ls <- split(df_query,df_query$Chromosome)
  df_seg_vsRef_new <- do.call(rbind,lapply(df_seg_vsRef.ls,function(x){
    x$integerCNV_ref <- fill_na_with_next(x$integerCNV_ref)
    x$integerCNV_ref <-fill_na_with_previous(x$integerCNV_ref)
    x$ratio_map_ref <- fill_na_with_next(x$ratio_map_ref)
    x$ratio_map_ref <-fill_na_with_previous(x$ratio_map_ref)
    
    return(x)
    }))
  if(reCalcuCNA | all(!grepl("ratio_map",colnames(df_seg_vsRef_new)))){
    df_seg_vsRef_new$SegMean_2ref <- df_seg_vsRef_new$SegMean*df_seg_vsRef_new$ratio_map_ref
    
    df_seg_vsRef_new <- CNV.align(df_seg_vsRef_new,CNbase,align_by = "SegMean_2ref")
    #df_seg_vsRef_new$integerCNV <- round(df_seg_vsRef_new$SegMean*df_seg_vsRef_new$integerCNV_ref)
    # if(is.na(delta_base)|delta_base==0){
    #   message("No delt provide, use ratio_map of df_ref for output.")
    #   
    # }else{
    #   df_seg_vsRef_new$ratio_map2 <- delta_base*df_seg_vsRef_new$integerCNV  ##TBC
    #   ratio_map_adj <- global_offset_correct(df_seg_vsRef_new$ratio_map,dat.ref=df_ref$ratio_map)
    #   df_seg_vsRef_new$ratio_map2 <- df_seg_vsRef_new$ratio_map -ratio_map_adj
    # }
    
  }
  
  df_seg_vsRef_new$diff <- abs(df_seg_vsRef_new$integerCNV_ref-df_seg_vsRef_new$integerCNV)
  df_seg_vsRef_new$diff[df_seg_vsRef_new$diff>0] <- 1
  diff_prop <- sum(df_seg_vsRef_new$diff)/nrow(df_seg_vsRef_new)

  if(diff_prop<= diff_prop_lim){
    return(list(df_seg_update = df_ref,prop_diff = diff_prop))
  }else{
    #update segmentation by combine those from clt_opt (baseline-cluster) and those from SegScore_OptRef
    segs_ref <- unique(df_ref$segName)
    segs_vsRef <- unique(df_query$segName)
    binSet <- df_binRatio[,c("binID","Chromosome","Start","End")]
    seg_process <- seg_split(unique(c(segs_ref,segs_vsRef)),binSet)
    df_binC <- seg_process$seg_binLevel
    
    #add binRatio and update SegMean
    binC_indx <- match(df_binC$binID,df_binRatio$binID)
    df_binC$binRatio <- df_binRatio[binC_indx,"binRatio"]
    # df_binC <- df_binC %>%
    #   dplyr::group_by(segName) %>%
    #   mutate(SegMean=median(binRatio,na.rm = T))%>% as.data.frame()
    
    bin_indx2 <- match(df_binC$binID,df_seg_vsRef_new$binID)
    df_binC <- cbind(df_binC,df_seg_vsRef_new[bin_indx2,c("SegMean_2ref","ratio_map","integerCNV")])
    colnames(df_binC)[grepl("SegMean_2ref",colnames(df_binC))] <- "SegMean"
    
    df_binC.ls <- split(df_binC,df_binC$Chromosome)
    df_seg_C_new <- do.call(rbind,lapply(df_binC.ls,function(x){
      x$integerCNV <- fill_na_with_next(x$integerCNV)
      x$integerCNV <-fill_na_with_previous(x$integerCNV)
      x$ratio_map <- fill_na_with_next(x$ratio_map)
      x$ratio_map <- fill_na_with_previous(x$ratio_map)
      return(x)
    }))
    df_seg_C_new$chrom <- gsub("chr","",df_seg_C_new$Chromosome)
    df_seg_C_new$chrom <- gsub("X","23",df_seg_C_new$chrom)
    df_seg_C_new$chrom <- gsub("Y","24",df_seg_C_new$chrom)
    df_seg_C_new$chrom <- as.numeric(df_seg_C_new$chrom)
    df_seg_C_new <- df_seg_C_new[order(df_seg_C_new$chrom,df_seg_C_new$Start),,drop=F]
    rownames(df_seg_C_new) <- df_seg_C_new$binID
    return(list(df_seg_update = df_seg_C_new,prop_diff = diff_prop))
  }
}



