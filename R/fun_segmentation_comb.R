suppressMessages({require(data.table)})

#' @title chromosome_order()
#' @export
chromosome_order <- function(chrom) {
  sapply(chrom, function(x) {
    if (x == "chrX") {
      return(23)
    } else if (x == "chrY") {
      return(24)
    } else {
      return(as.numeric(sub("chr", "", x)))
    }
  })
}

#' @title strings2bed()
#' @description convert strings (e.g., "chrx-xxx-xxx") to bed format
#' @export

strings2bed <- function(strings,split_by = "_|-|:",check.sort = TRUE){
    bin_bed <- strsplit(strings, "_|-|:")
    bin_bed <- do.call(rbind,bin_bed)
    bin_bed <- as.data.frame(bin_bed)
    colnames(bin_bed) <- c("chromosome", "start", "end")
    bin_bed$start <- as.numeric(bin_bed$start)
    bin_bed$end <- as.numeric(bin_bed$end)
    rownames(bin_bed) <- strings
    if(check.sort){
        bin_bed <- bin_bed %>%
        mutate(Chromosome_numeric = chromosome_order(chromosome)) %>%
        arrange(Chromosome_numeric, start, end) %>%
        select(-Chromosome_numeric)  %>%
        as.data.frame()

    }
    return(bin_bed)
}


#' @title DEsegs()
#' @description Test the difference in ratio values of all bins on each segment between the two groups.
#' @param data.frame with columns: bin name, segment name, ratio of group1, ratio of group2
#' @return mode=2 in the output data.frame means the segment with differential mean ratio between groups.
#' @export
DEsegs <- function(df,group1="refC",group2="obsC",
    wide2long_by="binID",
    p.adj_cutoff = 0.05,
    ratio_diff_cutoff=NULL,
    na_rm=TRUE){
    df.long <- suppressMessages(reshape2::melt(df, variable_name = "binID"))
    colnames(df.long) <- c("binID","segID","group","ratio")
    if(na_rm){
      df.long <- na.omit(df.long)
    }
    groupMean <- df.long %>%
                dplyr::group_by(segID,group)%>%
                dplyr::summarise(groupMean = median(ratio,na.rm=TRUE))%>% as.data.frame()
    diffMean <-  groupMean%>%
        group_by(segID)%>%
        summarise(diffMean =abs(groupMean[group%in%group2]- groupMean[group%in%group1]))%>% as.data.frame()

    ###differential test
    stat.test <- df.long %>%
        dplyr::add_count(segID, group) %>% 
        dplyr::filter(n >= 2) %>%
        dplyr::select(-n) %>%
        dplyr::ungroup()%>%
        dplyr::group_by(segID) %>%
        t_test(ratio ~ group) %>%
        adjust_pvalue(method = "BH") %>%
        add_significance()
    stat.test <- as.data.frame(stat.test)
    stat.test <- left_join(stat.test,diffMean,by="segID")       
    stat.test$mode <- 1

    if(!is.null(ratio_diff_cutoff)){
        stat.test$mode[stat.test$p.adj < p.adj_cutoff & stat.test$diffMean>=ratio_diff_cutoff] <- 2 #differential
     
    }else{
        stat.test$mode[stat.test$p.adj < p.adj_cutoff ] <- 2 #differential
    }

    return(stat.test)
} 


#' @title segCNV_refine()
#' @export

segCNV_refine <- function(clonalRes,CNest.ref,cytoBand,
    seg_method="PELT", # "robseg",#"gaussian", #
    seg.count.lim = 80,
    penalty = 1,
    segValue_method="median",
    FiltSeg= TRUE,
    SegLen_min=1e6,
    SegSize_min = 5,
    p.adj_cutoff = 0.05,
    ratio_diff_cutoff=NULL){
  
    clones <- names(clonalRes)
    diploidy_clones <- do.call(cbind,lapply(clonalRes,function(x){x$diploidy}))
    aneuploidy_clones <- clones[!diploidy_clones]
    cur_ref <- which(!diploidy_clones)[1]
    if(!is.na(cur_ref)){
      ref_clone <- clones[cur_ref]
      clones_ordered <- c(aneuploidy_clones,clones[diploidy_clones])
      res_ref <- clonalRes[[ref_clone]]
      
      BinSegRatio_ref <- res_ref$input_BinSegRatio
      BinSegRatio_ref <- BinSegRatio_ref[,!grepl("relativeCN|integerCN",colnames(BinSegRatio_ref))]
      
      segDat_ref <- res_ref$seg.dat
      segDat_ref<-segDat_ref[,c("segName","relativeCN","integerCN")]
      BinSegRatio_ref <-dplyr::left_join(BinSegRatio_ref,segDat_ref,by="segName")
      rownames(BinSegRatio_ref)<- BinSegRatio_ref$binID
      CNest_refclone <- unique(data.frame(ratio=BinSegRatio_ref$relativeCN,CN=BinSegRatio_ref$integerCN))
      CNest_refclone <- na.omit(CNest_refclone[order(CNest_refclone$CN),])
      if(nrow(CNest_refclone)>1){CNest_refclone <-CNbaseline.fill(CNest_refclone)}
        
      delt.ref <- (CNest.ref$ratio[2]-CNest.ref$ratio[1])/(CNest.ref$CN[2]-CNest.ref$CN[1])
      if(is.null(ratio_diff_cutoff)){ratio_diff_cutoff <- round(0.7*delt.ref,2)}
      
      if(length(clones_ordered)>1){ 
          for(cl in 2:length(clones_ordered)){
              cur_clone <- clones_ordered[cl]
              res_obs <- clonalRes[[cur_clone]]
              BinSegRatio_obs <- res_obs$input_BinSegRatio
              BinSegRatio_obs <- BinSegRatio_obs[,!grepl("relativeCN|integerCN|^w$",colnames(BinSegRatio_obs))]
              segDat_obs <- res_obs$seg.dat
              segDat_obs<-segDat_obs[,c("segName","relativeCN","integerCN","w")]
              BinSegRatio_obs <-dplyr::left_join(BinSegRatio_obs,segDat_obs,by="segName")
              rownames(BinSegRatio_obs)<- BinSegRatio_obs$binID
  
              CNest_obs <- unique(data.frame(ratio=segDat_obs$relativeCN,CN=segDat_obs$integerCN))
              CNest_obs <- CNest_obs[order(CNest_obs$CN),]
              if(nrow(CNest_obs)>1){
                CNest_obs <- CNbaseline.fill(CNest_obs) 
              }
  
              if(cur_clone %in% clones[diploidy_clones]){                
                  clonalRes[[cur_clone]]$seg.dat$integerCN <- 2
                  clonalRes[[cur_clone]]$seg.dat$relativeCN <- CNest_obs$ratio[CNest_obs$CN == 2]
              }else{               
                  ###
                  ### correct segments' CNV for ref_clone
                  if(cl == 2){
                      seg_df_ref <- BinSegRatio_ref
                      index <- match(seg_df_ref$binID,BinSegRatio_obs$binID)
                      seg_df_ref$obsC <- BinSegRatio_obs$binRatio[index]
                      seg_df_ref <- seg_df_ref[,c("binID","segName","binRatio","obsC")]
                      DEsegs_ref<- DEsegs(seg_df_ref,group1="binRatio",group2="obsC",
                          p.adj_cutoff = p.adj_cutoff,
                          ratio_diff_cutoff=ratio_diff_cutoff)
                      segs_diff_ref<- DEsegs_ref$segID[DEsegs_ref$mode==2]
                      #DEsegs_ref[grepl("chr16",DEsegs_ref$segID),]
                      # p0 <- seg_plot(na.omit(BinSegRatio_ref),name.data=paste0(ref_clone),
                      #             show_specific_seg=segs_diff_ref,
                      #              genome="hg38", ylim=NULL,
                      #              color_dot =F,plotDir="./",ylab = "Ratio",outPlot=F,color_seg_gradient=F)
                      # p0$ggarranged_p
  
                      ###combine common segs for ref_clone
                      seg_df_ref$segName_raw <- seg_df_ref$segName
                      seg_df_ref$segName[!seg_df_ref$segName %in% segs_diff_ref] <- NA
                      seg_df_ref$chr <- sapply(strsplit(seg_df_ref$binID,"_|-|:"),"[",1)
                      seg_df_ref <- seg_df_ref %>%
                          dplyr::group_by(chr) %>%
                          dplyr::mutate(
                              group = cumsum(!is.na(segName)),
                              group_label = ifelse(is.na(segName), paste0(chr,"_", group), as.character(segName))
                          ) %>%
                          dplyr::group_by(chr, group) %>%
                          dplyr::mutate(group_label = ifelse(is.na(segName), paste0(chr,"_", cur_group_id()), group_label)) %>%  
                          ungroup() %>%
                          dplyr::group_by(chr, group_label) %>% 
                          dplyr::mutate(Nbins = dplyr::n())%>%
                          ungroup() %>%
                          dplyr::select(-group)%>%as.data.frame()
                      ###Merge tiny fragments
                      row_index <- which(seg_df_ref$Nbins<5)
                      seg_df_ref <- seg_df_ref %>%
                          dplyr::group_by(chr) %>%
                          dplyr::mutate(segName = ifelse(Nbins < 5,  na.locf(na.locf(segName, na.rm = FALSE), fromLast = TRUE), segName))%>% 
                          ungroup() %>%
                          as.data.frame()
                      rownames(seg_df_ref)<- seg_df_ref$binID
  
                      regions_com <- unique(seg_df_ref[,c("segName","group_label")])
                      regions_com <- regions_com$group_label[is.na(regions_com$segName)]
                      for(x in 1:length(regions_com)){
                          region_x <- regions_com[x]
                          bins_x <- seg_df_ref$binID[seg_df_ref$group_label %in% region_x]
                          segs_x <- unique(seg_df_ref$segName_raw[seg_df_ref$group_label %in% region_x])
                          seg_datRef <- na.omit(unique(BinSegRatio_ref[bins_x,c("segName","integerCN","SegMean")]))
                          xCNV_ref <- seg_datRef$integerCN
                          row_index_ref <- which(BinSegRatio_ref$segName %in% segs_x)
                          row_index_refSeg <- which(res_ref$seg.dat$segName %in% segs_x)
                          
                          bin_dat_obs <- BinSegRatio_obs[bins_x,,drop=F ]
                          xCNV_obs_prop <- prop.table(table(bin_dat_obs$integerCN)) * 100
  
                          if(all(xCNV_ref==xCNV_ref[1])){
                            if(length(xCNV_ref)==1){
                                #First,compare with neighbor segment 
                                chr <- strsplit(segs_x,":|-|_")[[1]][1]
                                seg_dat_ref_chr <- res_ref$seg.dat[res_ref$seg.dat$chr %in% chr,c("chr","segName", "w","relativeCN", "integerCN"),drop=F]
                                segs_x_index <- which(seg_dat_ref_chr$segName ==seg_datRef$segName)
                                Neighbors_index <- c(segs_x_index-1,segs_x_index+1)
                                Neighbors_index <- Neighbors_index[Neighbors_index!=0&Neighbors_index<=nrow(seg_dat_ref_chr)]
                                Neighbors_index <- Neighbors_index[which.max(seg_dat_ref_chr$w[Neighbors_index])]
                                if(length(Neighbors_index)>0){
                                  Nei_xCNV_ref <- seg_dat_ref_chr$integerCN[Neighbors_index]
                                  Nei_seg_ref <- seg_dat_ref_chr$segName[Neighbors_index]
                                }else{
                                  Nei_xCNV_ref <- -1 #label NA
                                }
                                if(xCNV_ref != Nei_xCNV_ref){
                                
                                  #Are the groups consistent？
                                  if(length(xCNV_obs_prop)>0){
                                      dominant_xCNV_obs <- as.numeric(names(xCNV_obs_prop)[which.max(xCNV_obs_prop)])
                                  }else{dominant_xCNV_obs<- NA}
                                  #correcting CNV based on integrated groups, if the CNV of the same fragment in the two groups is inconsistent
                                  if((!is.na(dominant_xCNV_obs))&dominant_xCNV_obs!=xCNV_ref[1]){
                                    mat_sub <- seg_df_ref[bins_x,c("binRatio","obsC"),drop=F]
                                    mat_sub <- data.frame(meanratio=rowMeans(mat_sub,na.rm=TRUE))
                                    regionsMean_integ <- median(mat_sub[,1],na.rm=TRUE)
                                    xCNV_new <- CNest.ref$CN[which.min(abs(CNest.ref$ratio-regionsMean_integ))]
                                    BinSegRatio_ref$integerCN[row_index_ref] <- xCNV_new
                                    BinSegRatio_ref$relativeCN[row_index_ref] <- CNest_refclone$ratio[CNest_refclone$CN == xCNV_new]
    
                                    res_ref$seg.dat$integerCN[row_index_refSeg] <- xCNV_new
                                    res_ref$seg.dat$relativeCN[row_index_refSeg] <- CNest_refclone$ratio[CNest_refclone$CN == xCNV_new]
                                  }
                                }
                            }   
                              
                          }else{
                              #Check the ratio difference between fragments in region_x
                              bin_dat_ref <- BinSegRatio_ref[bins_x,c("binID","binRatio","segName","SegMean","integerCN")]
                              xCNV_ref_prop <- bin_dat_ref %>%
                                dplyr::count(segName,integerCN) %>% 
                                dplyr::mutate(CNprop = n / sum(n))%>% 
                                na.omit()
                              dominant_xCNV_ref <- xCNV_ref_prop[which.max(xCNV_ref_prop$CNprop),,drop=F]
                              SegMean_domin<- unique(BinSegRatio_ref$SegMean[BinSegRatio_ref$segName %in%dominant_xCNV_ref$segName])
  
                              for(seg_y in segs_x[!segs_x%in%dominant_xCNV_ref$segName]){
                                  bins_y <- bin_dat_ref$binID[bin_dat_ref$segName %in% seg_y]
                                  SegMean_y <- unique(BinSegRatio_ref$SegMean[BinSegRatio_ref$segName %in%seg_y])
                                  integerCN_y <- unique(BinSegRatio_ref$integerCN[BinSegRatio_ref$segName %in%seg_y])
                                  integerCN_y[is.na(integerCN_y)] <- -1 
                                  
                                  SegMean_diff <- abs(SegMean_domin-SegMean_y)
                                  if((integerCN_y != dominant_xCNV_ref$integerCN)|SegMean_diff>=ratio_diff_cutoff){
                                    
                                    row_index_ref_y <- which(BinSegRatio_ref$segName %in% seg_y)
                                    row_index_refSeg_y <- which(res_ref$seg.dat$segName %in% seg_y)
  
                                    if(integerCN_y>0 & length(row_index_ref_y)>=SegSize_min){
                                        segs_x.test <- bin_dat_ref[bin_dat_ref$segName %in% c(seg_y,dominant_xCNV_ref$segName),,drop=F] %>%
                                            t_test(binRatio ~ segName) %>%
                                            adjust_pvalue(method = "BH") %>%
                                            add_significance()%>%as.data.frame()
                                        if(segs_x.test$p.adj >= p.adj_cutoff ){ #| SegMean_diff<ratio_diff_cutoff
                                            BinSegRatio_ref$integerCN[row_index_ref_y] <- dominant_xCNV_ref$integerCN
                                            BinSegRatio_ref$relativeCN[row_index_ref_y] <- CNest_refclone$ratio[CNest_refclone$CN == dominant_xCNV_ref$integerCN]
        
                                            res_ref$seg.dat$integerCN[row_index_refSeg_y] <- dominant_xCNV_ref$integerCN
                                            res_ref$seg.dat$relativeCN[row_index_refSeg_y] <- CNest_refclone$ratio[CNest_refclone$CN == dominant_xCNV_ref$integerCN]
                 
                                        }else{
                                          #check if integerCN_y equals with integerCN_y_obs
                                          bin_dat_obs_y <- BinSegRatio_obs[bins_y,,drop=F ]                               
                                          xCNV_obs_prop_y <- prop.table(table(bin_dat_obs_y$integerCN)) * 100
                                          if(length(xCNV_obs_prop_y)>0){
                                            dominant_xCNV_obs_y <- as.numeric(names(xCNV_obs_prop_y)[which.max(xCNV_obs_prop_y)])
                                          }else{dominant_xCNV_obs_y<- NA}
                                          mat_sub_y <- seg_df_ref[bins_y,c("binID","segName","binRatio","obsC"),drop=F]
                                          DEsegs_y<- DEsegs(mat_sub_y,group1="binRatio",group2="obsC",
                                                              p.adj_cutoff = p.adj_cutoff,na_rm=FALSE)
  
                                          if(DEsegs_y$p.adj >= p.adj_cutoff &(!is.na(dominant_xCNV_obs_y))&dominant_xCNV_obs_y!=integerCN_y){
                                            
                                            mat_sub_y <- data.frame(meanratio=rowMeans(mat_sub_y[,c("binRatio","obsC"),drop=F],na.rm=TRUE))
                                            regionsMean_integ_y <- median(mat_sub_y[,1],na.rm=TRUE)
                                            xCNV_new_y <- CNest_refclone$CN[which.min(abs(CNest_refclone$ratio-regionsMean_integ_y))]
                                            
                                            row_index_ref_y <- which(BinSegRatio_ref$segName %in% seg_y)
                                            row_index_refSeg_y <- which(res_ref$seg.dat$segName %in% seg_y)
                                            BinSegRatio_ref$integerCN[row_index_ref_y] <- xCNV_new_y
                                            BinSegRatio_ref$relativeCN[row_index_ref_y] <- CNest_refclone$ratio[CNest_refclone$CN == xCNV_new_y]
                                            
                                            res_ref$seg.dat$integerCN[row_index_refSeg_y] <- xCNV_new_y
                                            res_ref$seg.dat$relativeCN[row_index_refSeg_y] <- CNest_refclone$ratio[CNest_refclone$CN == xCNV_new_y]
                                            
                                            
                                          }
                                            
                                        }
        
                                    }#else{
                                    #     BinSegRatio_ref$integerCN[row_index_ref_y] <- dominant_xCNV_ref$integerCN
                                    #     BinSegRatio_ref$relativeCN[row_index_ref_y] <- CNest_refclone$ratio[CNest_refclone$CN == dominant_xCNV_ref$integerCN]
                                    # 
                                    #     res_ref$seg.dat$integerCN[row_index_refSeg_y] <- dominant_xCNV_ref$integerCN
                                    #     res_ref$seg.dat$relativeCN[row_index_refSeg_y] <- CNest_refclone$ratio[CNest_refclone$CN == dominant_xCNV_ref$integerCN]
                                    # 
                                    # }
                                  }
                                  
                              }    
                          }
                      }
                      BinSegRatio_ref <- BinSegRatio_ref[,!grepl("relativeCN|integerCN",colnames(BinSegRatio_ref))]
                      res_ref$input_BinSegRatio <- BinSegRatio_ref
                      clonalRes[[ref_clone]] <- res_ref
                  
                      #correct segs_diff_ref regions for ref_clone
  
  
                  }
  
                  ###
                  ### correct segments' CNV for obs_clone
                  seg_df <- BinSegRatio_obs
                  index1 <- match(seg_df$binID,BinSegRatio_ref$binID)
                  seg_df$refC <- BinSegRatio_ref$binRatio[index1]
                  ###add obs ratio
                  index2 <- match(seg_df$binID,BinSegRatio_obs$binID)
                  seg_df$obsC <- BinSegRatio_obs$binRatio[index2]
                  seg_df<- seg_df[seg_df$w>=SegSize_min,,drop=F] #filt segs
                  # seg_df<- seg_df[seg_df$Num_bins>=round(SegSize_min/5),,drop=F] #filt segs
  
                  subdata <- seg_df[,c("binID","segName","refC","obsC")]
                  stat.test <- DEsegs(subdata,ratio_diff_cutoff=ratio_diff_cutoff)
                  segs_diff <- stat.test$segID[stat.test$mode==2]
                  # stat.test[grepl("chr1_",stat.test$segID),]
  
                  # p1 <- seg_plot(na.omit(seg_df),name.data=paste0(cur_clone),
                  #             show_specific_seg=segs_diff,
                  #              genome="hg38", ylim=NULL,
                  #              color_dot =F,plotDir="./",ylab = "Ratio",outPlot=F,color_seg_gradient=F)
                  # p1$ggarranged_p
  
                  ###replace the relativeCN and integerCN for the bins in segs_com
                  segs_com <- unique(BinSegRatio_obs$segName[!BinSegRatio_obs$segName %in% segs_diff])
              
                  for(si in 1:length(segs_com)){
                      segi <- segs_com[si]
                      row_index_obs <- which(BinSegRatio_obs$segName %in% segi)
                      row_index_obsSeg <- which(res_obs$seg.dat$segName %in% seg_y)
  
                      # cnvObs_s <- unique(BinSegRatio_obs$integerCN[row_index_obs])
  
                      bins_s <- BinSegRatio_obs$binID[row_index_obs]
                      ref_s <- BinSegRatio_ref[BinSegRatio_ref$binID %in%bins_s,,drop=F ]
                      refCN <- prop.table(table(ref_s$integerCN)) * 100
                      if(length(refCN)>0){
                          dominant_CN <- as.numeric(names(refCN)[which.max(refCN)])
                          #replace relativeCN and integerCN for s
                          BinSegRatio_obs$integerCN[row_index_obs] <- dominant_CN
                          if(dominant_CN %in% CNest_obs$CN){
                              BinSegRatio_obs$relativeCN[row_index_obs] <- CNest_obs$ratio[CNest_obs$CN == dominant_CN]
                          }else{
                              delta_obs <- (CNest_obs$ratio[2]-CNest_obs$ratio[1])/(CNest_obs$CN[2]-CNest_obs$CN[1])
                              CNest_obs_ratio <- CNest_obs$ratio[1]+(dominant_CN-CNest_obs$CN[1])*delta_obs
                              BinSegRatio_obs$relativeCN[row_index_obs] <- CNest_obs_ratio
                          }
                                          
                          res_obs$seg.dat$integerCN[row_index_obsSeg] <- dominant_CN
                          res_obs$seg.dat$relativeCN[row_index_obsSeg] <- CNest_obs$ratio[CNest_obs$CN == dominant_CN]
                      }   
                  }
                  BinSegRatio_obs <- BinSegRatio_obs[,!grepl("relativeCN|integerCN",colnames(BinSegRatio_obs))]
                  res_obs$input_BinSegRatio <- BinSegRatio_obs
                  clonalRes[[cur_clone]] <- res_obs
              }
          }
  
      }
  
  
      clusterout <-  lapply(names(clonalRes),function(cluster,clonalRes){
        res_ini <- clonalRes[[cluster]]
        seg_dat <- res_ini$seg.dat
        CNest <- unique(data.frame(ratio=seg_dat$relativeCN,CN=seg_dat$integerCN))
        CNest <- peakIndex(res_ini$input_BinSegRatio,CNest)
        CNest <- CNest[order(CNest$CN),]
        ploidy <- sum(seg_dat$integerCN*seg_dat$w,na.rm=TRUE)/sum(seg_dat$w,na.rm=TRUE)
          
        res_ini$CNest <- CNest
        res_ini$ploidy <- ploidy
    
        return(res_ini)
        },clonalRes)
    }
    return(clonalRes)

}

