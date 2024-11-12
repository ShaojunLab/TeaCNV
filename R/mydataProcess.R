#data pre-process
#input single cell data to matrix
#' @title load_scdata()
#' @export
load_scdata <- function(mtx_path,feature_path,barcode_path){
  features <- readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature)
  barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)
  
  mtx_ans <- Matrix::readMM(mtx_path) %>%
    magrittr::set_rownames(features$feature) %>%
    magrittr::set_colnames(barcodes$barcode)
  print(paste0("Raw counts: ", dim(mtx_ans)[1]," peaks and ",dim(mtx_ans)[2]," cells")) #574952 peaks and 3759
  return(mtx_ans)
}




#smoothing
#' @title smooth_by_chromosome()
#' @description uses the smooth_by_harmonic_mean() function to perform smoothing operations for sparse, large scATAC data.
#' @param mat_obj counts matrix with peaks(rows) and cells (columns).The row names should be with "chr1_247368351_247368851" format.
#' @param wl window length, which is the number of neighbor features to use for smoothing.
#' @param chr_annotation the chr for each feature
#' @param Ncores the number of cores if doParallel. Default is 1.
#' @param method smooth by "sum", "harmonic_mean", "maxDensity", or "loess" in a sliding window (these methods don't change the nrow of raw matrix)
#' or using "sum2bin" to sum the number of "wl" peaks(bin) and convert raw peak-cell matrix to bin-cell matrix ("sum2bin" decrease the nrow of raw matrix),
#' or using "winMean" to smooth by the mean in a sliding window, convert the peak-cell matrix to bin-cell matrix if win_step >1.
#' or using "weighted.mean" to smooth by the weighted mean in a sliding window, take length of peak as weight
#' default "sum2bin".
#' @param log_scale if method is "maxDensity",log_scale=TRUE, means density distribution based on log transformed value
#' @param win_step if method is "winMean", win_step is moving step of a sliding window
#' @keywords internal
#' @noRd
#' @export

smooth_by_chromosome <-function(mat_obj,wl=200,chr_annotation=NULL,Ncores=1,
                                  method="sum2bin",log_scale=TRUE,win_step=1,
                                  span=0.33,
                                  adjust_zero=TRUE
                                  ){
  if(!is.matrix(mat_obj)){ mat_obj <- as.matrix(mat_obj)}
  if(is.null(chr_annotation)){
    chr_annotation=sapply(strsplit(rownames(mat_obj),'_|-|:'),'[',1)
  }else if(length(chr_annotation)!=dim(mat_obj)[1]){
    stop("Error, the number of chromosome annotation for peaks does NOT match to the size of data")
  }
  #order chr
  chrs_uniq = unlist(unique(chr_annotation)) 
  if(grepl('chr',chrs_uniq[1])){
    chromnum=sapply(strsplit(chrs_uniq,'hr'),'[',2) ## if chrommnum with chr
    chromnum=chromnum[order(as.integer(gsub("[^0-9]", "", chromnum)))]
    chrs <- paste0("chr",chromnum)
  }else{
    chrs=chrs_uniq[order(as.integer(gsub("[^0-9]", "", chrs_uniq)))]
    }
  new_Mat <- c()
  for (chr in chrs) {
    chr_peaks_indices = which(chr_annotation %in% chr)
    #flog.info(paste0("smooth_by_chromosome: ",chr))
    chr_data=mat_obj[chr_peaks_indices, , drop=FALSE]
    if (nrow(chr_data) > 1) {
      if(method=="sum"){
        smoothed_chr_data <- .smooth_by_sum(chr_data,wl=wl,cores=Ncores)
      }
      if(method=="maxDensity"){
        smoothed_chr_data <- .smooth_by_maxDensity(matrix=chr_data,wl=wl,cores=Ncores,log_scale=log_scale)
      }
      if(method=="sum2bin"){
        smoothed_chr_data <- .smooth_by_sum2bin(Mat=chr_data,wl=wl)
      }
      if(method=="winMean"){
        smoothed_chr_data <- .smooth_by_windowMean(Mat=chr_data,wl=wl,step=win_step,zero.rm=adjust_zero)
      }
      if(method=="weighted.mean"){
        smoothed_chr_data <- .smooth_by_weighted_mean(Mat=chr_data,wl=wl,step=win_step)
      }
      if(method=="median"){
        library(matrixStats)
        smoothed_chr_data <- .smooth_by_windowMedian(Mat=chr_data,wl=wl,step=win_step)
      }
      if(method=="harmonic_mean"){
        smoothed_chr_data <- .smooth_by_harmonic_mean(Mat=chr_data,wl=wl,step=win_step)
      }
      if(method=="loess"){
        if(nrow(chr_data)<=10){span=1}
        smoothed_chr_data <- .smooth_by_loess(Mat=chr_data,span)
      }
    }else{
      smoothed_chr_data <- chr_data
    }
    new_Mat <- rbind(new_Mat,smoothed_chr_data)
  }
  return(new_Mat)

}

#' @title .smooth_by_weighted_mean()
#' @description calculate mean with weighted by the length of peaks within the window
#' wl is the window length, which is the number of peaks(default wl=20)
#' @keywords internal
#' @export

.smooth_by_weighted_mean<-function(Mat,wl=NULL,step=1,sep_by="-"){
  if(is.null(wl)){wl=20}
  row_start <- seq(1,nrow(Mat),by=step)
  mat_new=do.call(
    rbind,lapply(seq_len(length(row_start)), function(i){
      starti <- row_start[i]
      if(wl%%2==0){
        j1=starti-wl/2+1
      }else{j1=starti-wl%/%2}
      j2=starti+wl%/%2
      if(j1<=0)j1=1
      if(j2>=nrow(Mat))j2=nrow(Mat)
      if(starti+step-1<=nrow(Mat)){
        lastbin <- starti+step-1
      }else{lastbin <- nrow(Mat) }
      bini_end<-as.numeric(sapply(strsplit(rownames(Mat)[lastbin],":|_|-"),'[',3))
      newM=Mat[j1:j2,,drop=F]
      rownames(newM)<-rownames(Mat)[j1:j2]
      colnames(newM)<-colnames(Mat)
      chri <- unique(sapply(strsplit(rownames(newM),":|_|-"),'[',1))
      bini_start <- as.numeric(sapply(strsplit(rownames(Mat)[starti],":|_|-"),'[',2))
      binname <- paste(chri,bini_start,bini_end,sep=sep_by)
      #weighted_mean
      peak_s <- as.numeric(sapply(strsplit(rownames(newM),":|_|-"),'[',2))
      peak_e <- as.numeric(sapply(strsplit(rownames(newM),":|_|-"),'[',3))
      len_weight <- peak_e-peak_s
      
      weighted_sums <- colSums(newM * len_weight)
      total_weight <- sum(len_weight)
      weighted_means <- weighted_sums / total_weight
      
      
      bin_Hmean <- t(as.matrix(weighted_means))
      rownames(bin_Hmean) <- binname
      return(bin_Hmean)
    }))
  return(mat_new)
}
  
  
#'
#' @title .smooth_by_loess()
#' @description  loess smoothing
#' The size of the neighborhood can be controlled using the span argument, 
#' which ranges between 0 to 1. It controls the degree of smoothing. 
#' So, the greater the value of span, more smooth is the fitted curve. (default=1/3)
#' @param Mat counts matrix
#' @export
.smooth_by_loess <- function(Mat,span=0.33){
  mat_new=do.call(
    cbind,lapply(seq_len(ncol(Mat)), function(i){
      data_i <- data.frame(yvalue=Mat[,i],xindex=1:length(Mat[,i]))
      loessMod <- loess(yvalue ~ xindex, data=data_i, span=span) # 10% smoothing span
      #loessMod2 <- lowess(data_i$xindex,data_i$yvalue,f=span) # 10% smoothing span
      yvalue.new <- predict(loessMod) 
      # plot(data_i$xindex,data_i$yvalue,cex=0.5)
      # lines(x=data_i$xindex,y=yvalue.new,lty=1,col="blue")
      # lines(loessMod2, col=2, lty=2)
      # # Add a legend
      # legend(60, 0.19, legend=c("LOESS", "LOWESS"),
      #        col=c("blue", "red"), lty=1:2, cex=0.5)
       return(yvalue.new)
    }))
  colnames(mat_new)<-colnames(Mat)
  rownames(mat_new)<-rownames(Mat)
  return(mat_new)
}
#'
#' @title .smooth_by_sum()
#' @description calculate sum within the window
#' wl is the window length, which is the number of peaks(default wl=200)
#' @keywords internal
#' @export
#' 
.smooth_by_sum <- function(matrix,wl=NULL,cores=1){
  suppressMessages(library(doParallel))
  suppressMessages(library(doMC))
  mc <- getOption("mc.cores", cores)
  matrix0 <- do.call(
    cbind,mclapply(1:dim(matrix)[2], function(x){
      if(x%%100==0){cat(sprintf("Processing %1.f th cell, %1.f%% is done",x, x*100/dim(matrix)[2]),"\n")}
      if(is.null(wl)){wl=200} 
      cellx=as.vector(matrix[,x])
      cellx_new=c()
      for(i in 1:length(cellx)){
        if(wl%%2==0){
          j1=i-wl/2+1
        }else{j1=i-wl%/%2}
        j2=i+wl%/%2
        if(j1<=0)j1=1
        if(j2>=length(cellx))j2=length(cellx)
        wx=cellx[j1:j2]
        cellx_new[i]=sum(wx)
      }
      return(cellx_new)
    }, mc.cores = mc))
  colnames(matrix0)<-colnames(matrix)
  rownames(matrix0)<-rownames(matrix)
  return(matrix0)
}

#'
#' @title .smooth_by_sum2bin()
#' @description calculate sum within the window,
#' and shrink the raw matrix (rows) by combining the peaks to bins decided by the start and end site of peaks within the window
#' wl is the window length, which is the number of peaks(default wl=200)
#' @keywords internal
#' @export
#' 
.smooth_by_sum2bin <- function(Mat,wl=NULL,sep_by="-"){
  row_start <- seq(1,nrow(Mat),by=wl)
  mat_new=do.call(
    rbind,lapply(seq_len(length(row_start)), function(i){
    starti <- row_start[i]
    if(starti+wl<=nrow(Mat)){
      endi=starti+wl-1
    }else{endi=nrow(Mat)}
    if((nrow(Mat)-starti)>2){
      newM <- Mat[starti:endi, ]
      rownames(newM)<-rownames(Mat)[starti:endi]
      colnames(newM)<-colnames(Mat)
      chri <- unique(sapply(strsplit(rownames(newM),":|_|-"),'[',1))
      bini_start <- as.numeric(sapply(strsplit(rownames(newM)[1],":|_|-"),'[',2))
      bini_end<-as.numeric(sapply(strsplit(rownames(newM)[nrow(newM)],":|_|-"),'[',3))
      binname <- paste(chri,bini_start,bini_end,sep=sep_by)
      bin_sum <- t(as.matrix(colSums(newM)))
      rownames(bin_sum) <- binname
      return(bin_sum)
    }
    }))
  return(mat_new)
}

#'
#' @title .smooth_by_windowMean()
#' @description calculate MEAN within the window,
#' and convert the raw matrix (rows) by combining the peaks to bins decided by the start and end site of peaks within a window step
#' @param wl window length, which is the number of neighbor features to use for smoothing.
#' @param step the moving step of a sliding window
#' @keywords internal
#' Example:
#' Mat=matrix(seq(1:150),nrow=30)
#' rownames(Mat)=paste(rep("chr1",30),seq(1:30),seq(201,230),sep="_")
#' @export

.smooth_by_windowMean <- function(Mat,wl=NULL,step=1,zero.rm=FALSE,sep_by="-"){
  row_start <- seq(1,nrow(Mat),by=step)
  mat_new=do.call(
    rbind,lapply(seq_len(length(row_start)), function(i){
      starti <- row_start[i]
      if(wl%%2==0){
        j1=starti-wl/2+1
      }else{j1=starti-wl%/%2}
      j2=starti+wl%/%2
      if(j1<=0)j1=1
      if(j2>=nrow(Mat))j2=nrow(Mat)
      if(starti+step-1<=nrow(Mat)){
        lastbin <- starti+step-1
      }else{lastbin <-nrow(Mat) }
      bini_end<-as.numeric(sapply(strsplit(rownames(Mat)[lastbin],":|_|-"),'[',3))
      newM=Mat[j1:j2,,drop=F]
      rownames(newM)<-rownames(Mat)[j1:j2]
      colnames(newM)<-colnames(Mat)
      chri <- unique(sapply(strsplit(rownames(newM),":|_|-"),'[',1))
      bini_start <- as.numeric(sapply(strsplit(rownames(Mat)[starti],":|_|-"),'[',2))
       binname <- paste(chri,bini_start,bini_end,sep=sep_by)
       if(zero.rm){
         bin_mean <- t(as.matrix(apply(newM, 2, function(x){
           #the condition mean of Zero-Inflated distribution
           if(all(x==0)){
             adjusted_mean=0
           }else{
             non_zero_mean <- mean(x[x != 0])
             zero_proportion <- sum(x == 0) / length(x)
             adjusted_mean <- non_zero_mean *(1 - zero_proportion)
           }
           return(adjusted_mean)
         } )))
       }else{
      bin_mean <- t(as.matrix(colMeans(newM,na.rm = TRUE)))
       }
      rownames(bin_mean) <- binname
      return(bin_mean)
  }))
  return(mat_new)
}


#' @title .smooth_by_windowMedian()
#' @export
.smooth_by_windowMedian <- function(Mat,wl=NULL,step=1,sep_by="-"){
  row_start <- seq(1,nrow(Mat),by=step)
  mat_new=do.call(
    rbind,lapply(seq_len(length(row_start)), function(i){
      starti <- row_start[i]
      if(wl%%2==0){
        j1=starti-wl/2+1
      }else{j1=starti-wl%/%2}
      j2=starti+wl%/%2
      if(j1<=0)j1=1
      if(j2>=nrow(Mat))j2=nrow(Mat)
      if(starti+step-1<=nrow(Mat)){
        lastbin <- starti+step-1
      }else{lastbin <-nrow(Mat) }
      bini_end<-as.numeric(sapply(strsplit(rownames(Mat)[lastbin],":|_|-"),'[',3))
      newM=Mat[j1:j2,,drop=F]
      rownames(newM)<-rownames(Mat)[j1:j2]
      colnames(newM)<-colnames(Mat)
      chri <- unique(sapply(strsplit(rownames(newM),":|_|-"),'[',1))
      bini_start <- as.numeric(sapply(strsplit(rownames(Mat)[starti],":|_|-"),'[',2))
      binname <- paste(chri,bini_start,bini_end,sep=sep_by)
      bin_mean <- t(as.matrix(colMedians(newM,na.rm = TRUE)))
      rownames(bin_mean) <- binname
      return(bin_mean)
    }))
  colnames(mat_new) <- colnames(Mat)
  return(mat_new)
}

#' @title .smooth_by_weighted_mean()
#' @description calculate mean with weighted by the length of peaks within the window
#' wl is the window length, which is the number of peaks(default wl=20)
#' @keywords internal
#' @export

.smooth_by_weighted_mean<-function(Mat,wl=NULL,step=1,sep_by="-"){
  if(is.null(wl)){wl=20}
  row_start <- seq(1,nrow(Mat),by=step)
  mat_new=do.call(
    rbind,lapply(seq_len(length(row_start)), function(i){
      starti <- row_start[i]
      if(wl%%2==0){
        j1=starti-wl/2+1
      }else{j1=starti-wl%/%2}
      j2=starti+wl%/%2
      if(j1<=0)j1=1
      if(j2>=nrow(Mat))j2=nrow(Mat)
      if(starti+step-1<=nrow(Mat)){
        lastbin <- starti+step-1
      }else{lastbin <- nrow(Mat) }
      bini_end<-as.numeric(sapply(strsplit(rownames(Mat)[lastbin],":|_|-"),'[',3))
      newM=Mat[j1:j2,,drop=F]
      rownames(newM)<-rownames(Mat)[j1:j2]
      colnames(newM)<-colnames(Mat)
      chri <- unique(sapply(strsplit(rownames(newM),":|_|-"),'[',1))
      bini_start <- as.numeric(sapply(strsplit(rownames(Mat)[starti],":|_|-"),'[',2))
      binname <- paste(chri,bini_start,bini_end,sep=sep_by)
      #weighted_mean
      peak_s <- as.numeric(sapply(strsplit(rownames(newM),":|_|-"),'[',2))
      peak_e <- as.numeric(sapply(strsplit(rownames(newM),":|_|-"),'[',3))
      len_weight <- peak_e-peak_s
      
      weighted_sums <- colSums(newM * len_weight)
      total_weight <- sum(len_weight)
      weighted_means <- weighted_sums / total_weight
      
    
      bin_Hmean <- t(as.matrix(weighted_means))
      rownames(bin_Hmean) <- binname
      return(bin_Hmean)
    }))
  return(mat_new)
}

#' @title .smooth_by_harmonic_mean()
#' @description calculate harmonic mean with weighted by the number of non-zero within the window
#' wl is the window length, which is the number of peaks(default wl=20)
#' @keywords internal
#' @export

.smooth_by_harmonic_mean<-function(Mat,wl=NULL,step=1,sep_by="-"){
  row_start <- seq(1,nrow(Mat),by=step)
  mat_new=do.call(
    rbind,lapply(seq_len(length(row_start)), function(i){
      starti <- row_start[i]
      if(wl%%2==0){
        j1=starti-wl/2+1
      }else{j1=starti-wl%/%2}
      j2=starti+wl%/%2
      if(j1<=0)j1=1
      if(j2>=nrow(Mat))j2=nrow(Mat)
      if(starti+step-1<=nrow(Mat)){
        lastbin <- starti+step-1
      }else{lastbin <- nrow(Mat) }
      bini_end<-as.numeric(sapply(strsplit(rownames(Mat)[lastbin],":|_|-"),'[',3))
      newM=Mat[j1:j2,,drop=F]
      rownames(newM)<-rownames(Mat)[j1:j2]
      colnames(newM)<-colnames(Mat)
      chri <- unique(sapply(strsplit(rownames(newM),":|_|-"),'[',1))
      bini_start <- as.numeric(sapply(strsplit(rownames(Mat)[starti],":|_|-"),'[',2))
      binname <- paste(chri,bini_start,bini_end,sep=sep_by)
      #harmonic_mean
      nonZero <- colSums(newM != 0)
      bin_Hmean <- apply(newM,2,function(wx){
        if(any(wx!=0)){
          Hmean <- length(wx[wx!=0])/sum(1/wx[wx!=0])
          Hmean_weighted <-Hmean*length(wx[wx!=0])/wl   #weighted by the proportion
          x_new=Hmean_weighted
        }else{x_new=0}
      })
      bin_Hmean <- t(as.matrix(bin_Hmean))
      rownames(bin_Hmean) <- binname
      return(bin_Hmean)
  }))
  return(mat_new)
}
#


#'
#' @title .smooth_by_maxDensity()
#' @description  To smooth for sparse, large scATAC data per chromosome, identify the max density distributions within a window,
#'  and take the value under the max density to replace all the values within this window region.
#' @param matrix counts matrix with peaks(rows) and cells (columns).The row names should be with "chr1_247368351_247368851" format.
#' @param wl window length, which is the number of neighbor features to use for smoothing.
#' @param cores the number of cores if doParallel. Default is 1.
#' @keywords internal
#' @export
#' 
.smooth_by_maxDensity <- function(matrix,wl=NULL,cores=1,log_scale){
  suppressMessages(library(doParallel))
  suppressMessages(library(doMC))
  mc <- getOption("mc.cores", cores)
  matrix0 <- do.call(
    cbind,mclapply(1:dim(matrix)[2], function(x){
      if(x%%100==0){cat(sprintf("Processing %1.f th cell, %1.f%% is done",x, x*100/dim(matrix)[2]),"\n")}
      if(is.null(wl)){wl=200} 
      cellx=as.vector(matrix[,x])
      cellx_new=c()
      for(i in 1:length(cellx)){
        if(wl%%2==0){
          j1=i-wl/2+1
        }else{j1=i-wl%/%2}
        j2=i+wl%/%2
        if(j1<=0)j1=1
        if(j2>=length(cellx))j2=length(cellx)
        wx=cellx[j1:j2]
        if(length(wx<=1)){
          cellx_new[i]=wx
        }else{
          cellx_new[i]=.find_max_density(wx,log_scale=log_scale)
          }
      }
      return(cellx_new)
    }, mc.cores = mc))
  colnames(matrix0)<-colnames(matrix)
  rownames(matrix0)<-rownames(matrix)
  return(matrix0)
}

#' @title .find_max_density()
#' @param counts_analysis vector of data
#' @param log_scale if TRUE, density distribution based on log transformed value
#' @return the value at the max density 
#' @export
#'
.find_max_density <- function(counts_analysis,log_scale=T){
  if(log_scale){counts_analysis=log10(counts_analysis[!is.na(counts_analysis)])}
  Dy <- density(counts_analysis[!is.na(counts_analysis)])
  #plot(Dy)
  dy=Dy$y
  ind_max <-which(dy==max(dy))
  if(is.na(ind_max)){
    return(ind_max)
  }else{
    if(log_scale){
      return(10^(Dy$x[ind_max]))
    }else{return(Dy$x[ind_max])}
  }
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
      if(all(x==0)){
        adjusted_mean=0
      }else{
        non_zero_mean <- mean(x[x != 0])
        zero_proportion <- sum(x == 0) / length(x)
        adjusted_mean <- non_zero_mean *(1 - zero_proportion)
      }
      return(adjusted_mean)
    } )
  }else{
    peak_ref_grp_mean <- rowMeans(matrix_data,na.rm=TRUE)
  }
  
  peak_ref_grp_sd <- rowSds(matrix_data,useNames = TRUE)
  peak_ref_grp_sum <- rowSums(matrix_data,na.rm=TRUE)
  peak_ref_grp_means_sd <- data.frame(mean=peak_ref_grp_mean,SD=peak_ref_grp_sd,sum=peak_ref_grp_sum)
  rownames(peak_ref_grp_means_sd) <- RN
  return(peak_ref_grp_means_sd)
}


#' @title gen_bin_cell_mtx()
#' @description Generate bin by cell matrix for scATAC-seq data from peak by cell matrix
#' The matrix also can be generated by https://github.com/seasoncloud/Basic_CNV_SNV/blob/main/scripts/Gen_bin_cell_atac.R using "fragments.tsv.gz" file
#' @param mat_obj A count matrix of peak (or fragment) by cell
#' @param bin_size each bin size, i.e., sliding window size
#' @param window_step the step of the sliding window
#' @param chrom_size the data.frame of chromosome size with two columns: chr and size
#' @description if do fixed window bin, the bin_bed can be generate by "bedtools makewindows -g chrom.sizes -w 1000000 > genome.1Mb.bed"
#' or GenomicRanges::tileGenome() function with tilewidth=bin_size for chrosome size
#' @export
#' 
gen_bin_cell_mtx <- function(mat_obj,bin_size=NULL,window_step=NULL,
                             chrom_size=NULL,doFixBin=FALSE,sep_by="-"){
  require(GenomicRanges)
  if(!is.matrix(mat_obj)){ mat_obj <- as.matrix(mat_obj)}
  if(is.null(bin_size)){bin_size=1e6} #1Mb
  # extract peaks bed
  subject=GenomicRanges::GRanges(sapply(strsplit(rownames(mat_obj),":|_|-"),'[',1),
                  IRanges(as.numeric(sapply(strsplit(rownames(mat_obj),":|_|-"),'[',2))+1,
                          as.numeric(sapply(strsplit(rownames(mat_obj),":|_|-"),'[',3))))
  
  if(doFixBin){
    if(is.null(chrom_size)){
      stop(print("chrom_size must be specified"))
    }else{
      ## check bin bed chr column
      #order chr
      chrom_size=chrom_size[order(as.integer(gsub("[^0-9]", "", chrom_size[,1]))),]
      chrom_size=chrom_size[1:22,]
      seqlength=chrom_size[,2]
      names(seqlength)<- chrom_size[,1]
      subject_bins <- tileGenome(seqlength, tilewidth=bin_size, cut.last.tile.in.chrom=T)
      }
  }else{
    #convert bin matrix: for each chromosome, split the region from first peak start site to the end of last peak site into 1M bins
    rowidx=cumsum(seqnames(subject)@lengths)
    subject_size=c()
    for(chr in as.character(seqnames(subject)@values)){
      chr_idx <- which(seqnames(subject)==chr)
      start_idx=chr_idx[1]
      end_idx=chr_idx[length(chr_idx)]
      subject_size_chr <- c(chr,as.numeric(start(subject)[start_idx]),as.numeric(end(subject)[end_idx]))
      subject_size <-rbind(subject_size,subject_size_chr)
    }
    bins <- GenomicRanges::GRanges(as.character(subject_size[,1]),
                    IRanges(as.numeric(subject_size[,2]),
                            as.numeric(subject_size[,3])))
    sw=slidingWindows(x = bins, width = bin_size, step = window_step)
    subject_bins <- sw@unlistData
  }
  ## extract overlapping bins
  ov.bins <- subsetByOverlaps(subject_bins,subject)
  ov=findOverlaps(ov.bins, subject )
  ov=as.matrix(ov)
  smooth_mat=matrix(ncol=ncol(mat_obj), nrow=length(ov.bins))
  for(ii in 1:length(ov.bins)){
    ri=colSums(mat_obj[ov[which(ov[,1]==ii),2],, drop=FALSE]) #sum the total count within each bin
    smooth_mat[ii,]=ri
  }
  colnames(smooth_mat)=colnames(mat_obj)
  rownames(smooth_mat)=paste(as.character(seqnames(ov.bins)),start(ov.bins), end(ov.bins),sep=sep_by)
  invisible(gc())
  return(smooth_mat)
}



if (!require(strucchange)) install.packages("strucchange")
library(strucchange)
#find and return breakpoint indices
find_breakpoints.strucchange <- function(data,segSize.lim=NULL) {
  bp_results <- breakpoints(data ~ 1, h=segSize.lim)
  bp_indices <- bp_results$breakpoints
  bp_indices <- unique(c(bp_indices,length(data)))
  return(bp_indices)
}

find_breakpoints.changepoint <- function(x,method="PELT",penalty = "Manual",pen.value = "1.5 * log(n)",minseglen=5){
  if (!require(changepoint)) install.packages("changepoint")
  suppressMessages(library(changepoint))
  if(minseglen <= floor(length(x)/2)){
   # if(minseglen>(length(x)/2)){minseglen=floor(length(x)/2)}
    mv_pelt <- cpt.meanvar(x,method=method,penalty = penalty, pen.value = pen.value,minseglen=minseglen) 
    seg_end <- mv_pelt@cpts
  }else{
    seg_end <-length(x) 
  }
  return(seg_end)
}

find_breakpoints.robseg <- function(x,robseg_loss="Outlier",
                                    penalty=1){
  if (!require(robseg)) install.packages("robseg")
  suppressMessages(library(robseg))
  est.sd <- mad(diff(x)/sqrt(2))
  if(est.sd==0){est.sd <- 1e-6}
  ### Robust Fpop with the Biweight loss
  res.ou <- Rob_seg.std(x = x/est.sd,  
                        loss = robseg_loss, 
                        lambda = 1.5*penalty*log(length(x)), 
                        lthreshold=3)
  seg_end <- res.ou$t.est
  return(seg_end)
}


#' @title align_Grange2bin()
#' @param bed.query A bed file or a data frame of query regions
#' @param bed.subject A bed file or a data frame of subject regions
#' @return A data frame of the same order as the reference bed file, with the corresponding bin ID for each query region.
#' @export

align_Grange2bin = function(bed.query,bed.subject){
  options(scipen = 999)
  suppressMessages({
    library(GenomicRanges)
  })
  if (Reduce("|", is(bed.query) == "character")) {
    if(substr(bed.query, nchar(bed.query)-3, nchar(bed.query)) == ".bed") {
      bed.query <- read.table(bed.query,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
    }else if (Reduce("|", is(bed.query ) %in% c("data.frame"))) {
      bed.query  <- as.data.frame(bed.query )
    }
  }
  ## check bin bed chr column
  if(!grepl("chr",bed.query[1,1])){
    bed.query[,1]=paste0("chr",bed.query[,1])
  }
  if(!grepl("chr",bed.subject[1,1])){
    bed.subject[,1]=paste0("chr",bed.subject[,1])
  }
  ## convert to GRanges
  query=GRanges(bed.query[,1], IRanges(as.numeric(bed.query[,2]),as.numeric(bed.query[,3])))
  subject=GRanges(bed.subject[,1], IRanges(as.numeric(bed.subject[,2]),as.numeric(bed.subject[,3])))
  ## align
  align=findOverlaps(query,subject)
  idx_query<- queryHits(align)
  idx_subject<- subjectHits(align)
  
  colnames(bed.subject)[1:3] <-c("chr.subject","start.subject","end.subject")
  bed.query2 <- cbind(bed.query[idx_query,],bed.subject[idx_subject,])
  colnames(bed.query2)[1:3] = colnames(bed.query)[1:3] <-c("chr","start","end")
  bed.query2$binID <- paste(bed.query2$chr,bed.query2$start,bed.query2$end,sep="_")
  bed.query2 <- bed.query2[!duplicated(bed.query2$binID), ]
  
  bed.query$binID <- paste(bed.query[,1],bed.query[,2],bed.query[,3],sep="_")
  bed_final <- left_join(bed.query,bed.query2,by=intersect(colnames(bed.query),colnames(bed.query2)))
  
  return(as.data.frame(bed_final))
}
                                         
                                         
                                         

#' @title segment_by_chr_location()
#'
#' @description segment by computing the changepoint for an ordered vector
#'  (e.g., peaks or bins across chromosome p arm and q arm location)
#'  
#' @param mat_obj matrix with ordered by chromosome 'chr?-startCoordinate-endCoordinate' as rownames, cell Id as colnames. NA value was not allowed in the matrix.
#' @param cytoBand the reference annotation of chromosome location, BED format
#' @param acen_remove whether the features which located in Centromere region were removed (default is T)
#' @param clust_by_cell_group if TRUE, take the mean of group as the feature value
#' @param cell_group if clust_by_cell_group=T, a data.frame of cell_group needs to be provided, the first column is cell ID, the second column is cell groups
#' @param method use 'opart' package to segment which provides 2 functions:(i) opart_gaussian，(ii) opart_poisson
#' we use opart_gaussian(data, penalty) as default
#' @param  penalty A non-negative real number indicating penalty parameter, given the square loss (to
#' minimize) / gaussian likelihood (to maximize).较高的penalty通常会导致检测到的分段较少
#' @param min.Nseg the minimum number of segments
#' @param robseg_loss c("L1","L2","Huber","Outlier")
#' @return res a list of the input analysis matrix, and the changepoint matrix:a matrix of the changepoints (ends) of segment for each cell, value 1 means the changepoints
#'
#' @keywords internal
#' @noRd
#' @export
#'
segment_by_chr_location <- function(mat_obj,cytoBand=NULL,
                                    acen_remove=TRUE,
                                    clust_by_cell_group=FALSE,
                                    cell_group=NA,
                                    method=c("gaussian","poisson","robseg","strucchange","PELT"),
                                    robseg_loss="Outlier",
                                    segSize_min = 5,
                                    penalty=1
                                    ){
  suppressMessages({
    library(dplyr)
  })
  if(!is.matrix(mat_obj)){ mat_obj <- as.matrix(mat_obj)}
  myrowname <- rownames(mat_obj)
  raw_mat <- mat_obj
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
    cyto_gr <- cytoBand[which(cytoBand$gieStain %in% "acen"),]
    colnames(cyto_gr)[1:3]<- c("chr","start","end")
    
    
    # remove Centromere region 
    if(acen_remove){
      cytoBand = cytoBand[which(!cytoBand$gieStain %in% "acen"),]
    }
  }

  
 
  
  #segment for clusters
  if(clust_by_cell_group){
    if(any(is.na(cell_group))){
      stop("Error, clust_by_cell_group=T, but data.frame 'cell_group' with two columns of cell ID and cell group isn't provided")
    }else{
      unique_groups <- unique(cell_group[,2])
      groups=table(cell_group[,2])[unique_groups]
      mat_groups <- c()
      for(i in 1:length(unique_groups)){
        grp <- unique_groups[i]
        group_cells <- cell_group[which(cell_group[,2]  %in% grp),1]
        group_cells_idx <- which(colnames(mat_obj) %in% group_cells)
        group_feature_means <- .get_ref_mean(mat_obj,group_cells_idx)
        mat_groups <- cbind(mat_groups,group_feature_means$mean)
      }
      mat_groups <- as.data.frame(mat_groups)
      colnames(mat_groups) <- unique_groups
      mat_obj <- mat_groups
      rownames(mat_obj) <- myrowname
      
    }
  }
  #divide the features based on chromosome region
  chr_region_mat <- map2chrLocatio(mat_obj,cytoBand,merge_cyto=TRUE)
  changepoint_mat=matrix(ncol=ncol(mat_obj), nrow=nrow(mat_obj))
  colnames(changepoint_mat)=colnames(mat_obj)
  rownames(changepoint_mat)=rownames(mat_obj)

  for(ii in 1:length(chr_region_mat$subject_overlap)){ #for each region (p and q arm in each chromosome)
    mat_obj_sub=mat_obj[chr_region_mat$overlap[which(chr_region_mat$overlap[,1]==ii),2],, drop=F] 
    if(nrow(mat_obj_sub)>1){
      changepoint_mat_sub = apply(mat_obj_sub,2,function(x){
        if(method=="gaussian"){
          seg_end <- opart::opart_gaussian(x, penalty=penalty)$end.vec
        }else if(method=="poisson"){
          seg_end <- opart::opart_poisson(x, penalty=penalty)$end.vec
        }else if(method=="robseg"){
          seg_end <- find_breakpoints.robseg(x,robseg_loss=robseg_loss,penalty = penalty)
        }else if(method=="strucchange"){
          if(length(x)<2){seg_end <- length(x)
          }else{
            if(segSize_min>(length(x)/2)){segSize_min=max(floor(length(x)/2),2)}
            seg_end <- find_breakpoints.strucchange(x,segSize.lim=segSize_min)
            }
        }else if(method=="PELT"){
          seg_end <- find_breakpoints.changepoint(x,minseglen=segSize_min)
        }else{
          stop("method should be one of 'gaussian', 'poisson', 'strucchange', 'robseg' or 'PELT'")
        }
        y=rep(0,length(x))
        y[seg_end] <- 1
        return(y)
      })
      rownames(changepoint_mat_sub)=rownames(mat_obj_sub)
    }else{
      changepoint_mat_sub <- as.matrix(mat_obj_sub)
    }
    rowInd_changepoint_mat <- which(rownames(changepoint_mat) %in% rownames(changepoint_mat_sub))
    changepoint_mat[rowInd_changepoint_mat,] <- changepoint_mat_sub
    changepoint_mat[is.na(changepoint_mat)] <- 0
    
  }

  #add break by acen
  bins_input <- rownames(changepoint_mat)
  bins_chr <- sapply(strsplit(bins_input, "_|-|:"), "[", 1)
  bins_start <- sapply(strsplit(bins_input, "_|-|:"), "[", 2)
  bins_end <- sapply(strsplit(bins_input, "_|-|:"), "[", 3)
  
  bins_bed <- data.frame(chr=bins_chr,bins_start,bins_end)
  bins_acen <- align_Grange2bin(bed.query=cyto_gr,bed.subject=bins_bed)
  bins_acen <- na.omit(bins_acen)
  
  bins_input2 <- gsub("_|-|:","_",bins_input)
  binID_acen <- paste(bins_acen$chr.subject,bins_acen$start.subject,bins_acen$end.subject,sep ="_")
  bins_indx <- which(bins_input2 %in% binID_acen)
  changepoint_mat[bins_indx,] <- 1
  changepoint_mat[(bins_indx-1),] <- 1
  
  invisible(gc())
  res=list(data_raw=raw_mat,data_matrix=mat_obj,changepoint_matrix=changepoint_mat)
  if(method %in%c("gaussian","poisson")){
    res$penalty=penalty
  }
   if(clust_by_cell_group){
    res$Ncells_groups=groups
    res$metadata = cell_group
    }else{res$Ncells_groups=NA;res$metadata =NA}
  return(res)
}

#' @title map2chrLocatio()
#' @description map to chromosome location,return a GRanges subject
#' @param mat_obj matrix with 'chr?-startCoordinate-endCoordinate' as rownames,cell Id as colnames
#' @return res, list of GRanges subject_syto, subject_observe, subject_overlap, and overlap row index
#' @export
map2chrLocatio <- function(mat_obj,cytoBand=NULL,merge_cyto=TRUE){
  suppressMessages(require(GenomicRanges))
  subject=GenomicRanges::GRanges(sapply(strsplit(rownames(mat_obj),":|_|-"),'[',1),
          IRanges::IRanges(as.numeric(sapply(strsplit(rownames(mat_obj),":|_|-"),'[',2))+1,
                          as.numeric(sapply(strsplit(rownames(mat_obj),":|_|-"),'[',3))))
  subject_syto=GenomicRanges::GRanges(cytoBand$chrom,
               IRanges::IRanges(as.numeric(cytoBand$chromStart)+1,
                                as.numeric(cytoBand$chromEnd)))
  subject_syto$name <- cytoBand$name
  
  if(merge_cyto){
    subject_syto <- GenomicRanges::reduce(subject_syto) # merge ranges for p arm and q arm, separately
  }
  
  ## extract overlapping bins
  ov.bins <- subsetByOverlaps(subject_syto,subject)
  ov=findOverlaps(ov.bins, subject )
  ov=as.matrix(ov)
  colnames(ov) = c("ref_chr_cytoband","observed_matrix_row")
  res=list(subject_syto=subject_syto,subject_observe=subject,subject_overlap=ov.bins,overlap=ov)
  return(res)
}

#' @title get_seg_data()
#' @description get segment dataframe from the bin-cell matrix after segmentation with changepoint information
#' @param mtx_seg a list with 'data_matrix' bin-cell matrix; 'changepoint_matrix' with 1 means the segment end position; 'Ncells_groups',cell group if analysis based on groups (NA default)
#' @param outdir path for output file
#' @param Project  Project for the file name output
#' @param level the output file is "bin" level or "segment" (default) level
#' @export
get_seg_data = function(mtx_seg,outdir,Project,level="segment",plot=TRUE){
  suppressMessages(require(ggplot2))
  suppressMessages(require(ggpubr))
  suppressMessages(require(plyr))
  segfile_dir <- paste0("./outdir")
  if(!file.exists(segfile_dir)){dir.create(segfile_dir,recursive=T)}
  if(! level %in% "segment"){level="bin"}

  changepoint_mat <- mtx_seg$changepoint_matrix
  mtx_qcCR_sub <- mtx_seg$data_matrix
  Ncell_group <- mtx_seg$Ncells_groups
  bin_df <-data.frame(chr=sapply(strsplit(rownames(mtx_qcCR_sub),'-|_|:'),'[',1),
                      start=as.numeric(sapply(strsplit(rownames(mtx_qcCR_sub),'-|_|:'),'[',2)),
                      end=as.numeric(sapply(strsplit(rownames(mtx_qcCR_sub),'-|_|:'),'[',3)))
  bin_df$length <- bin_df$end-bin_df$start
  rownames(bin_df) <- rownames(mtx_qcCR_sub)
  bin_len <- bin_df$length
  genom_len <- sum(bin_len)

  for(j in 1:ncol(mtx_qcCR_sub)){
    if(!is.na(Ncell_group)){
      sub_cj <- as.numeric(colnames(mtx_qcCR_sub))[j]
      Ncell_cj <- Ncell_group[which(names(Ncell_group)==sub_cj)]
      plottitle=paste("Segmentation: ",Project,", sub_",sub_cj," (",Ncell_cj," cells)")
    }else{
      sub_cj <- colnames(mtx_qcCR_sub)[j]   
      plottitle=paste("Segmentation: ",Project,", sub_",sub_cj)
    }
   
    # png(paste0(plotDir,"/",de_noise_lb,"segment_",file_seg,"_sub_",sub_cj,".png"),width = 10, height = 3,units = 'in',res= 300)
    y_value = mtx_qcCR_sub[,j]
    
    cgp = changepoint_mat[,j]
    seg_end = which(cgp==1)
    seg_n <- length(seg_end)
    
    seg_value4bin=seg_id4bin=seg_len4bin= seg_name4bin<- rep(NA,nrow(mtx_qcCR_sub))
    segSet <- c()
    seg_start <- 1
    for(k in 1:seg_n){
      end_k <- seg_end[k]
      seg_k <- mean(y_value[seg_start:end_k])

      
      seg_mean <- seg_k
      seg_SE <- sd(y_value[seg_start:end_k])/sqrt(length(y_value[seg_start:end_k]))
      seg_nBin <- end_k-seg_start+1
      
      chr_k <- bin_df$chr[seg_start] 
      seg_start_loc <- bin_df$start[seg_start]
      seg_end_loc <- bin_df$end[end_k]
      seg_len_k <- sum(bin_df$length[seg_start:end_k])
      seg_name_k <- paste(chr_k,seg_start_loc,seg_end_loc,sep="_")
      
      segSet_k <- c(seg_name_k,chr_k,seg_start_loc,seg_end_loc,seg_len_k,seg_nBin,seg_mean,seg_SE)
      segSet <-rbind(segSet,segSet_k)

      
      seg_value4bin[seg_start:end_k] <- seg_k
      seg_id4bin[seg_start:end_k] <- k
      seg_len4bin[seg_start:end_k] <- seg_len_k
      seg_name4bin[seg_start:end_k] <- seg_name_k

      seg_start <- end_k+1
    }
    colnames(segSet)=c("segName","Chromosome","Start","End","cumul_SegLength","Num_bins","SegMean","SegSE")
    segOut <- paste0(segfile_dir,"/",Project,"_",level,"_Level.txt")

    if(!level %in% "segment"){
        geno_fraction4bin <- seg_len4bin/genom_len
        seg_bin_df <- data.frame(binID=rownames(mtx_qcCR_sub),binRatio=y_value,
                             length_bin=bin_df$length,
                             segID=seg_id4bin,segRatio=seg_value4bin,
                             segName=seg_name4bin,
                             length_seg=seg_len4bin,
                             seg_geno_fraction=geno_fraction4bin)
        write.table(seg_bin_df,segOut,col.names=T,row.names=F,quote=F,sep="\t")
      }else{
        write.table(segSet,segOut,col.names=T,row.names=F,quote=F,sep="\t")
      }
   
    #  ## scatter plot
    if(plot){
      outfile <- paste0(segfile_dir,"/segment_",Project,"_",sub_cj,".png")
      chromnum=sapply(strsplit(rownames(mtx_qcCR_sub),'-|_|:'),'[',1)
      chr_name=as.character(unique(chromnum))
      if(grepl('chr',chromnum[1])){
        chromnum=sapply(strsplit(chromnum,'hr'),'[',2) ## if chrommnum with chr
      }else{
        chromnum=chromnum}
      chromnum=chromnum[order(as.integer(gsub("[^0-9]", "", chromnum)))] #order by chr
      chromnum_xaxis = table(chromnum)[order(as.integer(gsub("[^0-9]", "", names(table(chromnum)))))]#order by chr
      chrline=cumsum(chromnum_xaxis)
      maploc=1:nrow(mtx_qcCR_sub)
      seg_geno_fraction <- as.numeric(segSet[,"cumul_SegLength"])/genom_len
      segments_df <- data.frame(SegMean=as.numeric(segSet[,"SegMean"]),seg_geno_fraction=seg_geno_fraction)
      segments_df_summary <- aggregate(segments_df[,-1,drop=F],segments_df["SegMean"],sum)

      breakss<-data.frame(c(0,chrline[1:(length(chrline))-1])+chromnum_xaxis*0.5)[,2]
      a<-data.frame(cbind(maploc,y_value,seg_value4bin))
      ylim=c(min(a$y_value,segments_df$segRatio),max(a$y_value,segments_df$segRatio))
      p1 <- ggplot(a) +
        geom_point(aes(x = factor(maploc), y = y_value), alpha = 1, shape = 16,size=1) +
        theme_bw() +
        labs(x = "", y = "Relative coverage")+
        ylim(ylim)+
        ggtitle(plottitle)+
        geom_line(aes(x = maploc,y=seg_value4bin),col="red")+
        geom_vline(xintercept =chrline,col="blue")+
        geom_hline(yintercept =1,col="grey",linetype = 2)+
        theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
              legend.position = "bottom",plot.title = element_text(size = 15),
              axis.title=element_text(size=15),
              axis.text = element_text(size = 15))+
        scale_x_discrete(breaks=breakss,labels=sapply(strsplit(chr_name,'hr'),'[',2))
      
      breakshist <- seq(min(ylim),max(ylim),length.out=20)
      p3=ggplot(segSet,aes(x=SegMean))+
        geom_histogram(aes(y=..density..),color="#E69F00",fill="wheat",breaks=breakshist)+
        labs(x = "", y = "Density")+
        ggtitle("Summary histogram")+
        coord_cartesian(xlim=ylim)+
        theme_bw() +
        theme_classic()+
        theme(legend.position = "bottom",plot.title = element_text(size = 15),
              axis.title=element_text(size=15),
              axis.text = element_text(size = 12))+
        coord_flip()
      pp=ggarrange(p1, p3, align ="h",ncol = 2, nrow = 1,widths = c(15,3),heights=3 )
      ggsave(outfile,pp,height = 3,width=14,dpi = 300)
    }
  }
  if(!level %in% "segment"){
    return(seg_bin_df)
    }else{
      return(segSet)
    }
}




#' @title map2chrLocatio()
#' @description convert gene symbol to gnomic Coordinate, vecter format
#' @param refFile Bed format in first three columns "chr" "start" "end"
#' @export
symbol2Coordinate = function(vect,refFile){
  if (Reduce("|", is(refFile) == "character")) {
    if(substr(refFile, nchar(refFile)-2, nchar(refFile)) %in% c(".gz","txt")){
      ref <- read.table(connection <- gzfile(refFile, 'rt'), sep="\t", header=F, check.names=FALSE)
      close(connection)
      ref <- as.data.frame(ref)
    }
  }else if(Reduce("|", is(refFile) %in% c("data.frame"))){
  }else{
    stop("Error, refFile isn't recognized as data.frame, or filename")
  }
  if(ncol(ref)==4){colnames(ref)=c("chr","start","end","symbol")}
  chrOrder<-paste0("chr",c(1:22,"X","Y"))
  ref=ref[ref[,1] %in% chrOrder,]
  chromnum <- ref[,1]
  chromnum=chromnum[order(as.integer(gsub("[^0-9]", "", sapply(strsplit(chromnum,"chr"),'[',2))))] 
  ref$chr<-factor(ref$chr, levels=chrOrder)
  ref=ref[order(ref$chr,ref$start),] #order by chr  # get 23118 genes   5
  ref$location = paste(ref$chr,ref$start,ref$end,sep="_")
  ref=ref[,c("symbol","location")]
  indx <- match(vect,ref$symbol)
  vect_ano=data.frame(symbol=vect,location=ref$location[indx])
  location = vect_ano$location

  return(location)
}

#' @title merge_atac()
#' @description merge ATAC objects
#' @export
#' 
merge_atac <- function(obj1,obj2,add_cell_id = NULL){
  suppressMessages({
    library(Signac)
    library(Seurat)
  })
  if(Reduce("|", is(obj2) %in% c("list"))){
    obj2_ls <- obj2
  }else if(Reduce("|", is(obj2) %in% c("Seurat"))){
    obj2_ls <-list(obj2)
  }
    
    obj <- merge(
    x = obj1,
    y = obj2_ls,
    add.cell.ids = add_cell_id
  )
  obj$dataset <- obj$orig.ident
  return(obj)
}


#' @title atac_integration()
#' @description merge ATAC objects
#' @export
#' 
atac_integration <- function(obj1,obj2,peakwidth.lim =c(200,1e4),reduction="rlsi",dims=2:30){
  #1.1 Creating a common peak set
  DefaultAssay(obj1) <- "ATAC"
  peaks_d1 <- rownames(obj1)
  gr_d1 <-  StringToGRanges(peaks_d1, sep = c("-", "-"))
  
  DefaultAssay(obj2) <- "ATAC"
  peaks_d2 <- rownames(obj2)
  gr_d2 <-  StringToGRanges(peaks_d2, sep = c("-", "-"))
  combined.peaks <- reduce(x = c(gr_d1, gr_d2))
  # Filter out bad peaks based on length
  peakwidths <- width(combined.peaks)
  combined.peaks <- combined.peaks[peakwidths  < peakwidth.lim[2] & peakwidths > peakwidth.lim[1]]
  # merge all datasets, adding a cell ID to make sure cell names are unique
  # if(is.null(add_cell_id)){
  #   add_cell_id <- c("obj1","obj2")
  # }
  peaks <- keepStandardChromosomes(combined.peaks, pruning.mode = "coarse")
  peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
  
  frags1 <- Fragments(obj1)
  counts1 <- FeatureMatrix(fragments=frags1,features=peaks, cells=colnames(obj1) )
  obj1_assay <- CreateChromatinAssay(counts1, fragments = frags1)
  obj1_new <- CreateSeuratObject(obj1_assay, assay = "ATAC", meta.data=obj1@meta.data,project =Project(obj1))
  obj1_new <- FindTopFeatures(obj1_new, min.cutoff = 10)
  obj1_new <- RunTFIDF(obj1_new)
  obj1_new <- RunSVD(obj1_new)
  
  frags2<- Fragments(obj2)
  counts2 <- FeatureMatrix(fragments=frags2,features=peaks, cells=colnames(obj2) )
  obj2_assay <- CreateChromatinAssay(counts2, fragments = frags2)
  obj2_new <- CreateSeuratObject(obj2_assay, assay = "ATAC", meta.data=obj2@meta.data,project =Project(obj2))
  obj2_new <- FindTopFeatures(obj2_new, min.cutoff = 10)
  obj2_new <- RunTFIDF(obj2_new)
  obj2_new <- RunSVD(obj2_new)
  
  obj.combined <- merge(obj1_new, obj2_new)
  obj.combined <- FindTopFeatures(obj.combined, min.cutoff = 10)
  obj.combined <- RunTFIDF(obj.combined)
  obj.combined <- RunSVD(obj.combined)
  # find integration anchors
  integration.anchors <- FindIntegrationAnchors(
    object.list = list(obj1_new, obj2_new),
    anchor.features = rownames(obj1_new),
    reduction = reduction,
    dims = dims
  )
  integrated <- IntegrateEmbeddings(
    anchorset = integration.anchors,
    reductions = obj.combined[["lsi"]],
    new.reduction.name = "lsi",
    dims.to.integrate = 1:30
  )
  
  integrated$dataset <- integrated$orig.ident
  
  integrated <- RunUMAP(integrated, reduction = "lsi", dims = 2:30)
  # obj <- merge(
  #   x = obj1,
  #   y = list(obj2),
  #   add.cell.ids = add_cell_id,
  #   merge.data = FALSE
  # )
  
  return(integrated)
}


#' @title dataProcess()
##input relative copy number ratio and breakpoint data
##chrom, start, end,ratio,breakpoint index
#' @param method "mean", "median" or "density" for calculate the value of segment.
#' @export
dataProcess <- function(CNVdata,method="median"){
  suppressMessages(library(bit64))
  coln <- colnames(CNVdata)  
  if(any(grepl('chrom|Chromosome',coln,ignore.case = T))){
    coln[grepl('chrom|Chromosome',coln,ignore.case = T)] <- 'chrom'
  }
  if(any(grepl('start',coln,ignore.case = T))){
    coln[grepl('start',coln,ignore.case = T)] <- 'start'
  }
  if(any(grepl('end',coln,ignore.case = T))){
    coln[grepl('end',coln,ignore.case = T)] <- 'end'
  }
  if(any(grepl('binRatio|ratio',coln,ignore.case = T))){
    coln[grepl('binRatio|ratio',coln,ignore.case = T)] <- 'ratio'
  }
  if(any(grepl('breakpoint',coln,ignore.case = T))){
    coln[grepl('breakpoint',coln,ignore.case = T)] <- 'breakpoint'
  }
  if(sum(grepl('chrom|start|ratio|breakpoint',coln))<4){
    stop("The columns should be 'chrom, start, end,ratio,breakpoint'")
  }
  colnames(CNVdata) <- coln
  #check breakpoint
  CNVdata <- CNVdata %>%
    group_by(segName) %>%
    mutate(breakpoint2 = if_else(row_number() == n(), 1, 0)) %>%
    ungroup()%>%as.data.frame()
  if(any(CNVdata$breakpoint != CNVdata$breakpoint2)){
    CNVdata$breakpoint <- CNVdata$breakpoint2
  }
  chromoPos=data.frame(chr=CNVdata$chrom,start=CNVdata$start,end=CNVdata$end,breakpoint=CNVdata$breakpoint)
  subseg=c()
  breakpoints=which(CNVdata$breakpoint==1)
  for (i in 1:length(breakpoints)){
    chr=chromoPos$chr[breakpoints[i]]
    if (i ==1){
      start=chromoPos$start[1]
      startindex=1
    }else{
      if(chromoPos$chr[breakpoints[i-1]+1] %in% chr){
        start=chromoPos$start[breakpoints[i-1]+1]
        startindex=breakpoints[i-1]+1 
      }else{
        chr_first_ind <- which(chromoPos$chr %in% chr)[1]
        start=chromoPos$start[chr_first_ind]
        startindex=chr_first_ind
        #add the last segment for the last chr with "0" breakpoint
        chr_mis <- chromoPos$chr[chr_first_ind-1]
        start_mis <- chromoPos$start[breakpoints[i-1]+1]
        end_mis <- chromoPos$end[chr_first_ind-1]
        startindex_mis <- breakpoints[i-1]+1
        endindex_mis <- chr_first_ind-1
        ratio_mis <- median(CNVdata$ratio[startindex_mis:endindex_mis])
        SD_mis <- sd(CNVdata$ratio[startindex_mis:endindex_mis])
        subseg=rbind(subseg,c(chr_mis,start_mis,end_mis,ratio_mis,SD_mis,endindex_mis-startindex_mis+1))
      }
    }
    end=chromoPos$end[breakpoints[i]]
    endindex=breakpoints[i]
    
    if(method=="median"){
     ratio=median(CNVdata$ratio[startindex:endindex],na.rm = T)
    }else if(method=="mean"){
      ratio <- mean(CNVdata$ratio[startindex:endindex],na.rm =TRUE)
    }else{
      if(length(CNVdata$ratio[startindex:endindex])>2){
        dens_seg <- density(CNVdata$ratio[startindex:endindex])
        max_index <- which.max(dens_seg$y)
        ratio <- dens_seg$x[max_index]
      }else{
        ratio <- mean(CNVdata$ratio[startindex:endindex],na.rm =TRUE)
      }
      
    }
    
    SD=sd(CNVdata$ratio[startindex:endindex])
    subseg=rbind(subseg,c(chr,start,end,ratio,SD,endindex-startindex+1))
  }
  subseg=subseg[subseg[,1]!=24,]
  subseg=data.frame(chr=subseg[,1],start=subseg[,2],end=subseg[,3],ratio=subseg[,4],SD=subseg[,5],bin=subseg[,6])
  chromosome=paste("chr",c(1:22,"X"),sep="")
  chromosome=chromosome[chromosome %in%unique(chromoPos$chr)]
  endlength <- c()
  for(i in chromosome){
    subdata=max(as.numeric(chromoPos$end[chromoPos$chr == i]))
    endlength <- cbind(endlength,subdata)
  }
  
  # endlength=as.numeric(apply(chromosome, function(i,chromoPos){
  #   subdata=max(as.numeric(chromoPos$end[chromoPos$chr == i]))
  #   return(subdata)
  # },chromoPos))
  
  endSum=c()
  endSum[1]=endlength[1]
  for(i in 2:length(endlength)){
    endSum[i]=sum(endlength[1:i])
  }
  chromoPos$abspos=as.numeric(chromoPos$start)
  for (i in 2:length(chromosome)){
    chromoPos$abspos[chromoPos$chr==chromosome[i]]=as.numeric(chromoPos$start[chromoPos$chr==chromosome[i]])+endSum[i-1]
  }
  subseg$abspos=as.numeric(subseg$start)
  subseg$absend=as.numeric(subseg$end)
  for (i in 2:length(chromosome)){
    subseg$abspos[subseg$chr==chromosome[i]]=as.numeric(subseg$start[subseg$chr==chromosome[i]])+endSum[i-1]
    subseg$absend[subseg$chr==chromosome[i]]=as.numeric(subseg$end[subseg$chr==chromosome[i]])+endSum[i-1]
  }
  return(list(seg=subseg,CNV=CNVdata,chromoPos=chromoPos))
}


#' @title seg_filter()
#' @description filter out the segment with low number of bins (min_bins_per_seg >1)
#' @param CNVdata a data.frame with segID or breakpoint in column
#' @export
seg_filter <- function(CNVdata,min_bins_per_seg=1){
  coln <- colnames(CNVdata) 
  if(any(grepl('segID',coln))){
    nbs <- table(CNVdata$segID)
    segID_removed <- names(nbs)[nbs<=min_bins_per_seg]
    CNVdata <- CNVdata[!CNVdata$segID %in% segID_removed,,drop=F]
  }else if(any(grepl('breakpoint',coln,ignore.case = T))){
    breakpoints<- which(CNVdata$breakpoint==1)
    breakpoints <- c(0,breakpoints)
    nbs <- diff(breakpoints)
    rowInd <- cumsum(nbs)
    names(nbs) <- rowInd
    segID_removed <- as.numeric(names(nbs)[nbs<=min_bins_per_seg])
    CNVdata <- CNVdata[-segID_removed,,drop=F]
  }else{message("No criteria for filtering segments. Return raw CNVdata!")}
  if(exists("segID_removed")&length(segID_removed)>0){
    print(paste0(length(segID_removed)," segmengs were removed ( ",length(nbs),"segments in total)."))
  }
  return(CNVdata)
}


#' @title split_grange()
#' @param bins_input strings of genomic regions,"chrx-xxx-xxx" format
#' @param cytoBand the reference annotation of chromosome location, BED format
#' @return bed, with the fourth column is bins_input, the second column is new bin ID, which
#' @export
#' breaks by acen
split_grange <- function(bins_input,cytoBand=NULL){
  options(scipen = 999)
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
  cyto_gr <- cytoBand[which(cytoBand$gieStain %in% "acen"),]
  colnames(cyto_gr)[1:3]<- c("chr","start","end")

  bins_chr <- sapply(strsplit(bins_input, "_|-|:"), "[", 1)
  bins_start <- sapply(strsplit(bins_input, "_|-|:"), "[", 2)
  bins_end <- sapply(strsplit(bins_input, "_|-|:"), "[", 3)
  bins_bed <- data.frame(chr=bins_chr,bins_start,bins_end)
  bins_acen2 <- align_Grange2bin(bed.query=bins_bed,bed.subject=cyto_gr)
  acen_rowInd <- which(!is.na(bins_acen2$gieStain))
  bins_acen2_filt <- bins_acen2[acen_rowInd,,drop=F]
  
  bins_acen2_new <- do.call(rbind,apply(bins_acen2_filt,1,function(x){
    start_bin <- as.numeric(x[2])
    end_bin <- as.numeric(x[3])
    start_acen <- as.numeric(x[6])
    end_acen <- as.numeric(x[7])

    if(start_bin>=start_acen){
      if(end_bin>end_acen){
        start <- end_acen+1; end <- end_bin
        }else{
         start <- NA; end <- NA 
        }
    }else{
      if(end_bin<=end_acen){
        start <- start_bin ; end <- start_acen-1
      }else{
        start <- c(start_bin,end_acen+1) ;
        end <- c(start_acen-1,end_bin)
      }
    }
    df <- data.frame(chr=x[1],start.new=start,end.new= end,binID=x[4])
    return(df)
  }))
  bins_df <- bins_acen2[,1:4]
  bins_df <- left_join(bins_df,bins_acen2_new,by=intersect(colnames(bins_df),colnames(bins_acen2_new)))
  rowInd <- which(!is.na(bins_df$start.new))
  bins_df$start[rowInd] <- bins_df$start.new[rowInd]
  bins_df$end[rowInd] <- bins_df$end.new[rowInd]
  bins_df$changed <- 0
  bins_df$changed[rowInd] <- 1
  bins_df$binID_new <- paste(bins_df[,1],bins_df[,2],bins_df[,3],sep="-")
  bins_df <- bins_df[,c(1:4,7:8)]
  colnames(bins_df)[4] <- "bins_input"
  return(bins_df)
}

#' breaks table of segments by acen
break_seg_tb <- function(seg_dat,cytoBandFile){
  segs <- seg_dat$segName
  segs_ano <- suppressWarnings(split_grange(segs,cytoBand=cytoBandFile))
  colnames(segs_ano)[2:4] <- c("start_new","end_new","segName")
  segs_ano$binID_new <- gsub("-","_",segs_ano$binID_new)
  seg_dat_new <- merge(seg_dat,segs_ano,by=intersect(colnames(seg_dat),colnames(segs_ano)),all=T)
  # seg_dat_new <- seg_dat_new %>%
  # dplyr::group_by(segName) %>%
  # fill(ratio,sd,relativeCN,integerCN, .direction = "downup") %>%
  # ungroup()
  row_changed <- which(seg_dat_new$changed==1)
  diff_start <- as.numeric(seg_dat_new$start_new) - as.numeric(seg_dat_new$start)
  diff_end <- as.numeric(seg_dat_new$end_new) - as.numeric(seg_dat_new$end)
  seg_dat_new$abspos[row_changed] <- as.numeric(seg_dat_new$abspos[row_changed])+diff_start[row_changed]
  seg_dat_new$absend[row_changed]<- as.numeric(seg_dat_new$absend[row_changed])+diff_end[row_changed]
  seg_dat_new$start[row_changed] <- seg_dat_new$start_new[row_changed]
  seg_dat_new$end[row_changed] <- seg_dat_new$end_new[row_changed]
  seg_dat_new$segName_raw <- seg_dat_new$segName
  seg_dat_new$segName <- seg_dat_new$binID_new
  seg_dat_new$chrom <- gsub("chr","",seg_dat_new$chr)
  seg_dat_new$chrom <- gsub("X|x","23",seg_dat_new$chrom)
  seg_dat_new$chrom <- gsub("Y|y","24",seg_dat_new$chrom)
  seg_dat_new <- seg_dat_new[order(as.numeric(seg_dat_new$chrom),as.numeric(seg_dat_new$start)),]
  seg_dat_new2 <- seg_dat_new[,c(colnames(seg_dat),"segName_raw","changed","chrom")]

  return(seg_dat_new2)
}


#' @title find_acen_bins()
#' @export
find_acen_bins <- function(bins,cytoBandFile){
  suppressPackageStartupMessages({
    library(dplyr)
  })

  bins_ano <- split_grange(bins,cytoBand=cytoBandFile)
  colnames(bins_ano)[2:4] <- c("start_new","end_new","binID")
  bins_ano$binID <- gsub("_","-",bins_ano$binID)
  bins_ano <- bins_ano[bins_ano$changed == 1,,drop=FALSE]
  binID_changed <- unique(bins_ano$binID)
  return(binID_changed)
}

#' @title find_1st_Trough()
#' @export
find_1st_Trough <- function(counts_analysis,log_scale=T){
  if(log_scale){counts_analysis=log10(counts_analysis)}
  Dy <- density(counts_analysis)
  # hist_obj <- hist(counts_analysis,breaks=500,plot=T)
  # highest_bin <- hist_obj$mids[which.max(hist_obj$counts)]
  # 
  #plot(Dy)
  dy=Dy$y
  dy_left <- dy[1:which(dy==max(dy))]
  ind_down <-which(diff(dy_left)<0)
  
  N_Trough = length(which(diff(ind_down)>1))
  if(N_Trough>0){
    dy_firstTrough <- ind_down[which(diff(ind_down)>1)[1]]+1
  }else{
    dy_firstTrough <- ind_down[length(ind_down)]
  }
  
  if(length(dy_firstTrough)==0){
    dy_right <- dy[which(dy==max(dy)):length(dy)]
    ind_up <- which(diff(dy_right)>0)[1]+which(dy==max(dy))-1
    dy_firstTrough <- ind_up
  }
  
  if(length(dy_firstTrough)==0){
    return(NA)
  }else{
    if(log_scale){
      return(10^(Dy$x[dy_firstTrough]))
    }else{return(Dy$x[dy_firstTrough])}
  }
}

#' @title denoise_via_ref_mean_sd_perFeature()
#'
#' @description Define noise based on the standard deviation of the reference cell expression data.
#' The range to remove noise would be mean +- sdev * sd_amplifier
#' where sd_amplifier expands the range around the mean to be removed as noise.
#' Data points defined as noise are set to zero.
#'
#' @param mat_obj matrix
#' @param ref_cells the reference cell id, same format with col names of mat_obj
#'
#' @param sd_amplifier multiplicative factor applied to the standard deviation to alter the noise
#'                     range (default: 1)
#'
#' @param noise_logistic uses a logistic (sigmoidal) function to noise removal.
#'
#' @return mat_obj
#'
#' @keywords internal
#' @noRd
#' @export
#'

denoise_via_ref_mean_sd_perFeature <- function(mat_obj,ref_cells=NULL, sd_amplifier=1.5) {
  if(!is.matrix(mat_obj)){ mat_obj <- as.matrix(mat_obj)}
  if (!is.null(ref_cells)) {
    ref_idx = which(colnames(mat_obj) %in% ref_cells)
    flog.info("denoising using mean(normal) +- sd_amplifier * sd(normal) per feature per cell across all data")
  }
  else {
    ref_idx = seq(1,ncol(mat_obj))
    flog.info("-no reference cells specified... using mean and sd of all cells as proxy for denoising")
  }
  vals = mat_obj[,ref_idx]
  
  mean_ref_vals = apply(vals, 1, function(x) mean(x, na.rm=TRUE))
  mean_ref_sd <- apply(vals, 1, function(x) sd(x, na.rm=TRUE)) * sd_amplifier
  
  upper_bound = mean_ref_vals + mean_ref_sd
  lower_bound = mean_ref_vals - mean_ref_sd
  
  flog.info(paste0(":: **** clear_noise_via_ref_quantiles **** : removing noise between bounds: mean_ref +- mean_ref_sd"))
  
  deno_mat_obj <- t(sapply(seq_len(nrow(mat_obj)),function(x){
    mat_obj_x <- mat_obj[x,]
    #mat_obj_x[mat_obj_x > lower_bound[x] & mat_obj_x < upper_bound[x]] = mean_ref_vals[x]
    mat_obj_x[mat_obj_x > lower_bound[x] & mat_obj_x < upper_bound[x]] = 1
    return(mat_obj_x)
  }))
  deno_mat_obj=as.matrix(deno_mat_obj)
  rownames(deno_mat_obj)=rownames(mat_obj)
  colnames(deno_mat_obj)=colnames(mat_obj)
  return(deno_mat_obj)
}


#' @title denoise_via_global_mean()
#' @param mat_obj cells in columns
#' @export
denoise_via_global_mean <- function(mat_obj,sd_amplifier=1){
  if(!is.matrix(mat_obj)){ mat_obj <- as.matrix(mat_obj)}
  mean_ref_vals = apply(mat_obj, 1, function(x) mean(x, na.rm=TRUE))
  mean_ref_sd = sd(mean_ref_vals)*sd_amplifier
  upper_bound = mean(mean_ref_vals) + mean_ref_sd
  lower_bound = mean(mean_ref_vals) - mean_ref_sd
  deno_mat_obj <- t(sapply(seq_len(nrow(mat_obj)),function(x){
    mat_obj_x <- mat_obj[x,]
    mat_obj_x[mat_obj_x > lower_bound[x] & mat_obj_x < upper_bound[x]] = mean(mean_ref_vals)
    return(mat_obj_x)
  }))
  deno_mat_obj=as.matrix(deno_mat_obj)
  rownames(deno_mat_obj)=rownames(mat_obj)
  colnames(deno_mat_obj)=colnames(mat_obj)
  #center
  d <- mean(mean_ref_vals)-1
  deno_mat_obj <- deno_mat_obj-d
  deno_mat_obj[deno_mat_obj<0] <- 0
  return(deno_mat_obj)
}


#' @title smooth_and_denoise()
#' @param mat_obj cells in columns
#' @export
smooth_and_denoise <- function(mat_obj,window=10,sd_amplifier=1){
  mt_de <- smooth_by_chromosome(mat_obj,wl=window,method="winMean")
  mt_de <- denoise_via_global_mean(mt_de,sd_amplifier=sd_amplifier)
  return(mt_de)
}




