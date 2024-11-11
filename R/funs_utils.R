##Functions utilized for HCC project.

#' @title global_offset_correct()
#' @description Correct the global CN ratio offset by the distance between 2 times the maximum peak of the histogram and the nearest integer CNV.
#' if provided a reference correction, then use the maximum peak of the reference data as target to calculate the distance.
#' @param dat a numeric vector of data to be correction
#' @param breaks a integer number for hist()
#' @return the global CN ratio offset
global_offset_correct <- function(dat,dat.ref=NULL){
  cen <- median(dat,na.rm = T)
  dist1 <- abs(cen-max(dat))
  dist2 <- abs(cen-min(dat))
  d <- max(dist1,dist2)*1.2
  ylim_l <- ifelse((cen-d) >=0,cen-d,0)
  ylim_h <- cen+d
  ylim=c(ylim_l,ylim_h)
  breakshist <- seq(min(ylim),max(ylim),length.out=80)
  hist_obj <- hist(dat,breaks=breakshist,plot=FALSE)
  highest_bin <- hist_obj$mids[which.max(hist_obj$counts)]
  if(is.null(dat.ref)){
    mv_dist <- (2*highest_bin-floor(2*highest_bin))/2
    if(mv_dist>0.25){mv_dist <- (2*highest_bin-ceiling(2*highest_bin))/2}
    if(mv_dist==0.25){mv_dist <- 0}
  }else{
    hist_obj.ref <- hist(dat.ref,breaks=breaks,plot=FALSE)
    highest_bin.ref <- hist_obj.ref$mids[which.max(hist_obj.ref$counts)]
    mv_dist <- highest_bin-highest_bin.ref
  }

  return(mv_dist)
}


#' @title infer_intCNV()
#' @param CNVdata input data.frame, bins on row, The columns should be 'chrom, start, end,ratio,breakpoint'
#'                    segName chrom                    binID Chromosome     Start       End binRatio
#' 1 chr1_144551589_169368392     1 chr1-155261887-155262788       chr1 155261887 155262788 2.647409
#' 2 chr1_144551589_169368392     1 chr1-155250380-155251334       chr1 155250380 155251334 4.546832
# length_bin segID  SegMean length_seg seg_geno_fraction seg_SD seg_score breakpoint    length
# 1        901     3 3.403457     179850        0.01924034  1.147  166.5213          0 248956422
# 2        954     3 3.403457     179850        0.01924034  1.147  166.5213          0 248956422
#' @param SegLen_cutoff quantile of segment length to be filtered out
#' 
# source("dataProcess.R")
# source("results.check.R")
# source("findmode.R")
# source("optEstimation.R")
infer_intCNV <- function(CNVdata,SegLen_cutoff=NULL,shift_correc = TRUE,delta_cutoff=0.4,delta_lim= c(0.4,0.7)){
  Dsets <- seq(min(delta_lim),max(delta_lim),by=0.005)
  if(!is.null(SegLen_cutoff)){
    ###filter segments
    CNVdata <- CNVdata %>%
      dplyr::group_by(segID) %>%
      do(data.frame(., w = length(.$segID)))%>%
      dplyr::filter(length_seg>=SegLen_cutoff)%>%
      as.data.frame()
  }
  CNVres=dataProcess(CNVdata)
  seg.dat=CNVres$seg
  seg.dat=data.frame(chr=seg.dat$chr,start=seg.dat$start,end=seg.dat$end,ratio=as.numeric(seg.dat$ratio),sd=as.numeric(seg.dat$SD),w=as.numeric(seg.dat$bin),abspos=as.numeric(seg.dat$abspos),absend=as.numeric(seg.dat$absend))
  seg.dat$sd[is.na(seg.dat$sd)]=0
  
  seg.dat_new <- seg.dat
  if(shift_correc){
    mv_dist <- global_offset_correct(seg.dat$ratio)
    
    #update the ratio
    if(mv_dist!=0){
      seg.dat_new$ratio <- as.numeric(seg.dat_new$ratio) - (mv_dist)
      seg.dat_new <- seg.dat_new[seg.dat_new$ratio>0,,drop=F]
    }
  }else{
    mv_dist = 0
  }
  output<-tryCatch(doCNV1_v2(seg.dat_new),error=function(e) NA )
  if(any(is.na(output))){stop("Faild to Estimate!")}
  
  if(output$CNVestimate$delta < delta_cutoff){
    outer_res=c()
    for (d in Dsets ){
      newInteCNV <- output$seg.dat$ratio/d
      outer = sum((newInteCNV - output$seg.dat$integerCNV)^ 2)
      outer_res <- c(outer_res,outer)
    }
    min_index=which.min(outer_res)
    delta_new <- Dsets[min_index]
    integerCNV <- round(output$seg.dat$ratio/delta_new)
    output$seg.dat$integerCNV <- integerCNV
    if(mv_dist!=0){
      output$seg.dat$CNV <- delta_new*integerCNV + mv_dist
    }else{ output$seg.dat$CNV <- delta_new*integerCNV}
    
    output$CNVestimate$delta <- delta_new
  }
  return(output)
}

#' @title mt2ratio_by_col()
#' @description calculate bulk-level ratio from single-cell count matrix based on global mean
#' @param mat count maxtrix with cell on column, and feature on row
#' @param col_idx integer vector of column start index and end index
#' @param ref.value numeric data.frame with the length of nrow(mat)
#source("fun_pseudoBulk_mat.R")
mt2ratio_by_col <- function(mat,col_idx=NULL,ref.value=NULL){
  if(is.null(col_idx)){
    col.start <- 1
    col.end <- ncol(mat)
  }else{
    col.start <- col_idx[1]
    col.end <- col_idx[2]
  }
  mat <- mat[,col.start:col.end,drop=F]
  mtx_b <- pseudo_bulk_v2(mat,group_ano = data.frame(colnames(mat),group=1),method ="mean",adjust_zero=TRUE)
  count.med <-  median(mtx_b[,1],na.rm=TRUE)
  if(is.null(ref.value)){
    ref.value <- count.med
  }
  binRatio <-  mtx_b/ref.value
  binRatio_finit <- binRatio$`1`
  is.na(binRatio_finit) <- sapply(binRatio_finit, is.infinite)
  is.na(binRatio_finit) <- sapply(binRatio_finit, is.nan)
  segRatio <- median(binRatio_finit,na.rm = T)
  return(list(binRatio=binRatio$`1`,segRatio=segRatio,median.count=count.med))
}


#' @title genGIF()
#' @description generate gif figures from a image path
#' 
genGIF <- function(path,outfile="./figure.gif"){
  require(magick)
  ## list file names and read in
  imgs <- list.files(path, full.names = TRUE)
  img_list <- lapply(imgs, image_read)
  img_joined <- image_join(img_list)
  ## animate at 2 frames per second
  img_animated <- image_animate(img_joined, fps = 2)
  ## view animated image
  #img_animated
  ## save to disk
  image_write(image = img_animated,
              path = outfile)
}

#' @title Correct_by_length()
#' @description correct count by peak length
#' @param mat peak-cell count matrix, rowname is in "chrx-xxx-xxx" format
#' @param unit per base, kb or Mb, default is 1000b (kb).
Correct_by_length <- function(mat,unit=1000){
  rows <- rownames(mat)
  start <- as.numeric(sapply(strsplit(rows,"_|:|-"),"[",2))
  end <- as.numeric(sapply(strsplit(rows,"_|:|-"),"[",3))
  len <- end-start
  
  mat_new<- mat*unit/len
  return(mat_new)
}

#' @title genProliferationScore()
genProliferationScore <- function(object,gene.list=NULL,assay="RNA"){
  library(Seurat)
  DefaultAssay(object) <- assay
  if(is.null(gene.list)){
    genes <- c("ZWINT", "E2F1", "FEN1", "FOXM1", "H2AFZ", "HMGB2", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MKI67", "MYBL2", "PCNA", "PLK1", "CCND1", "AURKA", "BUB1", "TOP2A", "TYMS", "DEK", "CCNB1","CCNE1")
    genes = intersect(genes,rownames(object))
    if(length(genes)==0){stop("No genes found in features.")}
    gene.list <- list(genes)
  }
  
  object <- AddModuleScore(
    object = object,
    features = gene.list,
    name = 'ProliferationScore')
  gsub('ProliferationScore1','ProliferationScore',colnames(object@meta.data))->colnames(object@meta.data)
  return(object)
}


#' @title cnv_score()
#' @description Score cells based on their CNV profile, then impute malignant cells.
#' @param mat CNV matrix (cells in columns, genes or peaks in rows)
#' @param diploid.base diploid value
cnv_score <- function(mat,diploid.base=NULL,diploid.base.sd=NULL,Var.cutoff=NULL,ref_cells=NULL,plot_out=TRUE,outdir="./",scaleFactor=1){
  if(is.null(diploid.base)){
    diploid.base <- mean(colMeans(mat,na.rm = TRUE))
  }
  #degree of variation
  culm_diff_use <- apply(mat,2,function(x){ 
    y=(as.numeric(x)-diploid.base)^2
    sigma <- sum(y)*scaleFactor/length(y)
    #sigma <- sum(y)
    return(sigma)
  })
  #Degree of intracellular variation
  cnv_sd <- apply(mat,2,function(x){ 
    return(sd(x))
  })
  culm_diff <- data.frame(CNV_score= culm_diff_use,CNV_SD=cnv_sd)
  value <- "CNV_score"
  cutoff <- Var.cutoff
  if(!is.null(cutoff)){
    ref_cells <- rownames(culm_diff)[culm_diff$CNV_score<cutoff]
    culm_diff$cnv_group <- ifelse(rownames(culm_diff)%in% ref_cells,"reference","observation")
    
  }else{
    if(!is.null(ref_cells)){
      culm_diff$cnv_group <- ifelse(rownames(culm_diff)%in% ref_cells,"reference","observation")
      Vref <- culm_diff[culm_diff$cnv_group=="reference",value]
      Vobs <- culm_diff[culm_diff$cnv_group=="observation",value]
      #find the cross point
      xlim <- c(quantile(Vref,0.01),quantile(Vref,0.99))
      df <- merge(
        as.data.frame(density(Vref, from = xlim[1], to = xlim[2])[c("x", "y")]),
        as.data.frame(density(Vobs, from = xlim[1], to = xlim[2])[c("x", "y")]),
        by = "x", suffixes = c(".a", ".b")
      )
      df$comp <- as.numeric(df$y.a > df$y.b)
      df$cross <- c(NA, diff(df$comp))
      cutoff <- df[which(df$cross != 0), "x"][1] #the cutoff to determine tumor cells from normal
      
      
    }else if(!is.null(diploid.base.sd)){
      cutoff <- diploid.base.sd^2
      ref_cells <- rownames(culm_diff)[culm_diff$CNV_score<cutoff]
      culm_diff$cnv_group <- ifelse(rownames(culm_diff)%in% ref_cells,"reference","observation")
    }else{
      cutoff <- find_1st_Trough(culm_diff_use)
      ref_cells <- rownames(culm_diff)[culm_diff$CNV_score<cutoff]
      culm_diff$cnv_group <- ifelse(rownames(culm_diff)%in% ref_cells,"reference","observation")
      
    }
  }

  culm_diff$cnv_impute <- ifelse((culm_diff$CNV_score>=cutoff) & culm_diff$cnv_group%in%"observation","Malignant","Non-malignant")
  culm_diff <- culm_diff[,!grepl("cnv_group",colnames(culm_diff)),drop=F]
  if(plot_out){
    
    p <- ggplot(culm_diff) +
      geom_density(aes(x = get(value), fill = cnv_impute),alpha=0.5)+
      theme_bw()+
      theme(legend.position="right",plot.title = element_text(hjust = 0.5),)+
      geom_vline(xintercept=cutoff, linetype="dashed", color = "red")+
      ggtitle("aneuploidy score")+
      xlab(value)+
      xlim(0, quantile(culm_diff[[value]],0.99))
    ggsave(paste0(outdir,'/CNVscore_distribution.pdf'),p,height = 3, width = 5)
  }

  return(culm_diff)
}



