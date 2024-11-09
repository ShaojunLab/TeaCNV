#'
#' @title seg_score()
#' @description score each segment based on the segment length (bin size within the segment),
#' and the standard deviation of bin values in this segment
#' @param data an list output from function 'segment_by_chr_location()', which contains 'changepoint_matrix', 'data_matrix',
#'  and 'Ncells_groups'. The rownames of matrix are 'chr1_237513_1530360' format.
#'  @param method is "mean", "median" or "density". For "density", we will take the density$x from the peak of density as the segMean.
#'  @return res list of two data.frame: segment score at bin-level and segment-level
seg_score = function(data,outname=NA,outdir="./",method="density"){
  #filter chr
  changepoint_mat <- data$changepoint_matrix
  changepoint_mat <-changepoint_mat[grepl("chr",rownames(changepoint_mat)),,drop=F]
  data.mtx <- data$data_matrix
  data.mtx <-data.mtx[grepl("chr",rownames(data.mtx)),,drop=F]
  #Ncell_group <- data$Ncells_groups
  
  bin_df <-data.frame(chr=sapply(strsplit(rownames(data.mtx),'-|_|:'),'[',1),
                      start=as.numeric(sapply(strsplit(rownames(data.mtx),'-|_|:'),'[',2)),
                      end=as.numeric(sapply(strsplit(rownames(data.mtx),'-|_|:'),'[',3)))
  bin_df$length <- bin_df$end-bin_df$start
  rownames(bin_df) <- rownames(data.mtx)
  bin_len <- bin_df$length
  genom_len <- sum(bin_len)
  
  res = list()
  for(j in 1:ncol(data.mtx)){
    sub_cj <- colnames(data.mtx)[j]
    changp_cj <- changepoint_mat[,j]
    #Ncell_cj <- Ncell_group[which(names(Ncell_group) %in% sub_cj)]
    y_value = data.mtx[,j]
    
    cgp = changepoint_mat[,j]
    seg_end = which(cgp==1)
    #the number of bins of each segment
    seg_size = c(seg_end[1],diff(seg_end))
    seg_n <- length(seg_end)
    
    seg_value4bin=seg_id4bin=seg_len4bin= seg_name4bin= seg_SD4bin=seg_score4bin<- rep(NA,nrow(data.mtx))
    segSet <- c()
    seg_start <- 1
    for(k in 1:seg_n){
      end_k <- seg_end[k]
      if(method=="mean"){
        seg_mean <- mean(y_value[seg_start:end_k],na.rm =TRUE)
      }else if(method=="median"){
        seg_mean <- median(y_value[seg_start:end_k],na.rm =TRUE)
      }else{
        if(length(y_value[seg_start:end_k])>2){
          dens_seg <- density(y_value[seg_start:end_k])
          max_index <- which.max(dens_seg$y)
          seg_mean <- dens_seg$x[max_index]
        }else{
          seg_mean <- mean(y_value[seg_start:end_k],na.rm =TRUE)
        }

      }
      
      seg_value4bin[seg_start:end_k] <- seg_mean
      seg_id4bin[seg_start:end_k] <- k
      
      seg_SD <- sd(y_value[seg_start:end_k])
      seg_SE <- seg_SD/sqrt(length(y_value[seg_start:end_k]))
      seg_nBin <- end_k-seg_start+1
      seg_score <- seg_nBin/seg_SD
      
      chr_k <- bin_df$chr[seg_start] 
      seg_start_loc <- bin_df$start[seg_start]
      seg_end_loc <- bin_df$end[end_k]
      seg_len_k <- seg_end_loc-seg_start_loc#sum(bin_df$length[seg_start:end_k])
      seg_name_k <- paste(chr_k,seg_start_loc,seg_end_loc,sep="_")
      
      segSet_k <- cbind(seg_name_k,chr_k,seg_start_loc,seg_end_loc,seg_len_k,seg_nBin,seg_mean,seg_SE,seg_SD,seg_score)
      segSet <-rbind(segSet,segSet_k)
      seg_len4bin[seg_start:end_k] <- seg_len_k
      seg_name4bin[seg_start:end_k] <- seg_name_k
      seg_SD4bin[seg_start:end_k] <- seg_SD
      seg_score4bin[seg_start:end_k] <- seg_score
      
      seg_start <- end_k+1
      rm(chr_k,seg_start_loc,seg_end_loc)
    }
    geno_fraction4bin <- seg_len4bin/genom_len
    seg_bin_df <- data.frame(binID=rownames(data.mtx),
                            Chromosome=bin_df$chr,
                            Start=bin_df$start,
                            End=bin_df$end,
                            binRatio=y_value,
                            length_bin=bin_df$length,
                             segID=seg_id4bin,SegMean=seg_value4bin,
                             segName=seg_name4bin,
                             length_seg=seg_len4bin,
                             seg_geno_fraction=geno_fraction4bin,
                             seg_SD = seg_SD4bin,
                             seg_score = seg_score4bin,
                             breakpoint = changp_cj)
    
    segSet = as.data.frame(segSet)
    colnames(segSet)=c("segName","Chromosome","Start","End","length_seg","Num_bins","SegMean","SegSE","SegSD","Score")
    columns <-c("Start","End","length_seg","Num_bins","SegMean","SegSE","SegSD","Score")
    segSet[,columns] <- sapply(segSet[,columns], as.numeric)
    #if(!is.na(outname)){
    #  segOut_seglevel <- paste0(outdir,"/",outname,"_",sub_cj,"_segLevel.txt")
    #  segOut_binlevel <- paste0(outdir,"/",outname,"_",sub_cj,"_binLevel.txt")

    #  #write.table(seg_bin_df,segOut_binlevel,col.names=T,row.names=F,quote=F,sep="\t")
    #  write.table(segSet,segOut_seglevel,col.names=T,row.names=F,quote=F,sep="\t")
    #}
    res$seg_score_binLevel[[sub_cj]] <- seg_bin_df
    res$seg_score_segLevel[[sub_cj]]<- segSet
  }
  #res=list(seg_score_binLevel,seg_score_segLevel)
  return(res)
}
