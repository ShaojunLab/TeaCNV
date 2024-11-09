suppressPackageStartupMessages({
  library(reshape)
  library(fpc)
  library(LaplacesDemon)
  library(diptest)
  library(mixtools)
  library(dplyr)
  library(tidyverse)
  library(rstatix)
  library(ggpubr)
})

#' @title regionSPlit()
#' @description check if there is breakpoint for genomic region 
#' @param regiondata: the ratio of bin
#' @param mu: the expected ratio of modes
#' @return update mu for the new split region
#'          data: the ratio of new split region with new region cluster ID
#'          clusterindex: 0 no further split
#' 

regionSPlit <- function(regiondata,mu){
  if (length(mu) == 1){
    cluster.info <- list(mu=mu,data=regiondata,clusterindex=0)
  }else{
    ##mu has three cases if length(mu) > 1:
    ##(1) the region has breakpoint so bin distribution differ;
    ##(2) the region has subclone so bin distribution is uniform;
    ##(3) both (1) and (2)
    if (max(regiondata$pos)-min(regiondata$pos) < 5){# no check if the region less than 5MB
      mu <- mean(mu)
      cluster.info <- list(mu=mu,data=regiondata,clusterindex=0)
    }else{
      if (max(mu)-min(mu) < 0.25){ # no check if the difference between modes less than 0.25
        mu <- mean(mu)
        cluster.info <- list(mu=mu,data=regiondata,clusterindex=0)
      }else{
        mixres <- normalmixEM(regiondata$ratio, lambda = 1/length(mu), mu = mu, sigma = 0.1)
        if (min(mixres$lambda) <= 0.2){
          mu <- mean(mu)
          cluster.info <- list(mu=mu,data=regiondata,clusterindex=0)
        }else{
          if (max(mixres$mu)-min(mixres$mu) < 0.25){
            mu <- mean(regiondata$ratio)
            cluster.info <- list(mu=mu,data=regiondata,clusterindex=0)
          }else{
            breakpoint <- NULL
            if (length(mu) == 2){##case (1) or (2)
              comres <- apply(mixres$posterior,1,which.max)
              pos1 <- regiondata$pos[comres==1]
              pos2 <- regiondata$pos[comres==2]
              p.value <- wilcox.test(pos1,pos2)$p.value
              if (p.value < 0.01){#case(1)
                break1 <- quantile(pos1,pro=0.05)
                break2 <- quantile(pos2,pro=0.05)
                breakpoint <- max(break1,break2)
              }##else case (2)
            }else{##case (3)
              comres <- apply(mixres$posterior,1,which.max)
              pos1 <- regiondata$pos[comres==1]
              pos2 <- regiondata$pos[comres==2]
              pos3 <- regiondata$pos[comres==3]
              pos.include <- c()
              p.value <- c()
              if (length(pos1) > 10 & length(pos2) > 10){
                p.value1 <- wilcox.test(pos1,pos2)$p.value 
                pos.include <- c(pos.include,1)
                p.value <- c(p.value,p.value1)
              }
              if (length(pos2) > 10 & length(pos3) > 10){
                p.value2 <- wilcox.test(pos2,pos3)$p.value
                pos.include <- c(pos.include,2)
                p.value <- c(p.value,p.value2)
              }
              if (length(pos3) > 10 & length(pos1) > 10){
                p.value3 <- wilcox.test(pos1,pos3)$p.value
                pos.include <- c(pos.include,3)
                p.value <- c(p.value,p.value3)
              }
              p.index <- which(p.value < 0.01)
              if (length(p.index) == 3){
                break1 <- quantile(pos1,pro=0.05)
                break2 <- quantile(pos2,pro=0.05)
                break3 <- quantile(pos3,pro=0.05)
                breakpoint <- c(break1,break2,break3)
                breakpoint <- breakpoint[order(breakpoint)]
                breakpoint <- breakpoint[2:3]
              }else if (length(p.index) == 2){
                k <- setdiff(c(1,2,3),p.index)
                if (k == 1){
                  pos11 <- c(pos1,pos2)
                  pos21 <- pos3
                }else if (k == 2){
                  pos11 <- c(pos2,pos3)
                  pos21 <- pos1
                }else{
                  pos11 <- c(pos1,pos3)
                  pos21 <- pos2
                }
                break1 <- quantile(pos11,pro=0.05)
                break2 <- quantile(pos21,pro=0.05)
                breakpoint <- max(break1,break2)
              }else{
                if (pos.include == 1){
                  pos11 <- pos1
                  pos21 <- pos2
                }else if (pos.include == 2){
                  pos11 <- pos2
                  pos21 <- pos3
                }else{
                  pos11 <- pos1
                  pos21 <- pos3
                }
                break1 <- quantile(pos11,pro=0.05)
                break2 <- quantile(pos21,pro=0.05)
                breakpoint <- max(break1,break2)
              }
            }
            if (!is.null(breakpoint)){
              cluster.info <- list()
              if (length(breakpoint) == 1){
                cluster1 <- regiondata[regiondata$pos < breakpoint,]
                cluster2 <- regiondata[regiondata$pos >= breakpoint,]
                mu1 <- Modes(cluster1$ratio)$modes
                mu2 <- Modes(cluster2$ratio)$modes
                cluster1$cluster <- paste0(cluster1$cluster,".1")
                cluster2$cluster <- paste0(cluster2$cluster,".2")
                clusterindex1 <- length(mu1)-1
                clusterindex2 <- length(mu2)-1
                cluster.info[[length(cluster.info)+1]] <- list(mu=mu1,data=cluster1,clusterindex=clusterindex1)
                cluster.info[[length(cluster.info)+1]] <- list(mu=mu2,data=cluster2,clusterindex=clusterindex2)
              }else{
                cluster1 <- regiondata[regiondata$pos < breakpoint[1],]
                cluster2 <- regiondata[regiondata$pos >= breakpoint[1]&regiondata$pos < breakpoint[2],]
                cluster3 <- regiondata[regiondata$pos >= breakpoint[2],]
                mu1 <- Modes(cluster1$ratio)$modes
                mu2 <- Modes(cluster2$ratio)$modes
                mu3 <- Modes(cluster3$ratio)$modes
                cluster1$cluster <- paste0(cluster1$cluster,".1")
                cluster2$cluster <- paste0(cluster2$cluster,".2")
                cluster3$cluster <- paste0(cluster3$cluster,".3")
                clusterindex1 <- length(mu1)-1
                clusterindex2 <- length(mu2)-1
                clusterindex3 <- length(mu3)-1
                cluster.info[[length(cluster.info)+1]] <- list(mu=mu1,data=cluster1,clusterindex=clusterindex1)
                cluster.info[[length(cluster.info)+1]] <- list(mu=mu2,data=cluster2,clusterindex=clusterindex2)
                cluster.info[[length(cluster.info)+1]] <- list(mu=mu3,data=cluster3,clusterindex=clusterindex3)
              }
            }else{
              cluster.info <- list(mu=mu,data=regiondata,clusterindex=0)
            }
          }
        }
      }
    }
  }
  return(cluster.info)
}

#' @title regionMerge()
#' @description merg region split result
#' @param cluster.res: output of regionSPlit from all initial bin region


regionMerge <- function(cluster.res){
  chrratio <- c()
  mu <- list()
  clusterindex <- list()
  for (i in 1:length(cluster.res)){
    res <- cluster.res[[i]]
    if (is.null(names(res))){
      for (j in 1:length(res)){
        chrratio <- rbind(chrratio,res[[j]]$data)
        mu[[length(mu)+1]] <- res[[j]]$mu
        clusterindex[[length(clusterindex)+1]] <- res[[j]]$clusterindex
      }
    }else{
      chrratio <- rbind(chrratio,res$data)
      mu[[length(mu)+1]] <- res$mu
      clusterindex[[length(clusterindex)+1]] <- res$clusterindex
    }
  }
  clusterID <- unique(chrratio$cluster)
  index <- match(chrratio$cluster,clusterID)
  chrratio$cluster <- index
  merge.res <- list(data=chrratio,mu=mu,clusterindex=clusterindex)
  return(merge.res)
}


#' @title regionSubclone()
#' @description check the subclonality for genome region
#' 
regionSubclone <- function(chrratio,cluster,overall.mu){
  subratio<- chrratio[chrratio$cluster==cluster,]
  if (dim(subratio)[1] > 50){
    cloneID <- unique(subratio$clone)
    mu.ref <- overall.mu[[cluster]]
    seg.subclone <- do.call(rbind,lapply(cloneID, function(clone,subratio,mu.ref){
      sub <- subratio[subratio$clone==clone,]
      mu<- Modes(sub$ratio)$modes
      if (length(mu) > 1){
        if (max(mu)-min(mu) < 0.25){
          mu=mean(sub$ratio)
        }else{
          if (length(sub$ratio) < 30){
            mu=mean(sub$ratio)
          }else{
            mixres <- normalmixEM(sub$ratio, lambda = 1/length(mu), mu = mu[order(mu)], sigma = 0.1)
            if (min(mixres$lambda)>0.4){
              if (max(mixres$mu)-min(mixres$mu)<0.25){
                mu=mean(sub$ratio)
              }
            }else{
              mu=mean(sub$ratio)
            }
          }
        }
      }
      if (length(mu) == 1){
        if (is.na(mu)){
          mu = mean(sub$ratio)
        }
        diff <- abs(mu - mu.ref)
        diff <- min(diff)
        if (diff > 0.25){
          cloneindex <- 3
        }else{
          cloneindex <- 1
        }
        overall.index <- mu.ref[which.min(diff)]
      }else{
        if (length(mu.ref)>1){
          cloneindex <- NA
        }else{
          cloneindex <- 2
        }
        overall.index <- NA
      }
      return(c(cloneindex,overall.index))
    },subratio,mu.ref))
    seg.subclone <- seg.subclone[!is.na(seg.subclone[,1])&!is.na(seg.subclone[,2]),]
    #subratio$ratio <- abs(subratio$ratio)
    if (length(unique(seg.subclone[,1]))>1){
      p.value <- bartlett.test(ratio ~ clone, data = subratio)$p.value
      if (p.value < 0.05){
        return(cluster) 
      }else{
        return(NULL)
      }
    }else{
      if (unique(seg.subclone[,1])==2){
        p.value <- bartlett.test(ratio ~ clone, data = subratio)$p.value
        if (p.value < 0.05){
          return(cluster) 
        }else{
          return(NULL)
        }
      }
      if (unique(seg.subclone[,1])==3){
        p.value <- bartlett.test(ratio ~ clone, data = subratio)$p.value
        if (p.value < 0.05){
          return(cluster) 
        }else{
          return(NULL)
        }
      }
      if (unique(seg.subclone[,1])==1){
        if (length(unique(seg.subclone[!is.na(seg.subclone[,2]),2]))>1){
          p.value <- bartlett.test(ratio ~ clone, data = subratio)$p.value
          if (p.value < 0.05){
            return(cluster) 
          }else{
            return(NULL)
          }
        }else{
          return(NULL)
        }
      }
    } 
  }else{
    return(NULL)
  }
}


#' @title genomicRegion()
#' @description collect all genome region  
#' @param ratiodata: initial clone list
#' @return all genomic binID
#' 
genomicRegion <- function(ratiodata,chr){
  bindata <- c()
  for (i in 1:length(ratiodata)){
    subdata <- ratiodata[[i]]$input_BinSegRatio
    subdata <- subdata[subdata$Chromosome==chr,]
    bindata <- rbind(bindata,subdata[,c(1,3,4)])
  }
  bindata <- unique(bindata)
  bindata <- bindata[order(bindata$Start),]
  return(bindata)
}


#' @title cloneRatio()
#' @description collect overall ratio at clone level  
#' @param ratiodata: initial clone list
#' @return all ratio value frm all clone and genomic binID
#' 

cloneRatio <- function(ratiodata,bindata,chr){
  chrratio <- c()
  for (i in 1:length(ratiodata)){
    subdata <- ratiodata[[i]]$input_BinSegRatio
    subdata <- subdata[subdata$Chromosome==chr,]
    index <- match(subdata$binID,bindata$binID)
    chrratio <- rbind(chrratio,data.frame(posindex = index,ratio =subdata$binRatio,clone=i))
  }
  chrratio$pos <- bindata$Start[chrratio$posindex]/1000000
  #chrratio$ratio[chrratio$ratio>4]=4
  return(chrratio)
}


#' @title relativeRatio.Ref()
#' @description relative rato data 
#' @param ratiodata: all ratio value frm all clone and genomic binID
#' 

relativeRatio.Ref <- function(chrratio,ref=1){
  posindex <- unique(chrratio$posindex)
  cloneID <- unique(chrratio$clone)
  chrratiowide <- reshape(chrratio[,1:3], idvar = "clone", timevar = "posindex", direction = "wide")
  refclone <- chrratiowide[ref,2:dim(chrratiowide)[2]]
  ratio.ref <- do.call(rbind,apply(chrratiowide[,2:dim(chrratiowide)[2]],1,function(x,refclone){
    abs(x-refclone)
  },refclone))
  ratio.ref <- as.data.frame(ratio.ref)
  row.names(ratio.ref) <- cloneID
  colnames(ratio.ref) <- posindex
  ratio.ref$clone <- cloneID
  ratio.reflong <- melt(ratio.ref, id = c("clone"), variable_name = "posindex")
  ratio.reflong$pos <- chrratio$pos[match(ratio.reflong$posindex,chrratio$posindex)]
  return(ratio.reflong)
}


#' @title genomicPartition()
#' @description genome bin partition based on overall ratio  
#' @param chrratio: initial clone list
#' @return genomic region partitioning
#' 

genomicPartition <- function(chrratio){
  db <- fpc::dbscan(chrratio[,c(4,2)], eps = 8, MinPts = 20)
  chrratio$cluster = db$cluster
  clusterID <- unique(db$cluster)
  clusterID <- clusterID[clusterID!=0]
  cluster.res <- lapply(clusterID, function(cluster,chrratio){
    sub <- chrratio[chrratio$cluster==cluster,]
    mu <- Modes(sub$ratio)$modes
    cluster.info <- regionSPlit(sub,mu)
    return(cluster.info)
  },chrratio)
  
  merge.res <- regionMerge(cluster.res)
  return(merge.res)
}

#' @title renewGroup()
renewGroup <- function(list_data){
  is_equal <- function(x, y) {
    all(sort(x) == sort(y))
  }
  
  # 
  new_index <- numeric(length(list_data))
  current_index <- 1
  
  
  for (i in seq_along(list_data)) {
    if (i == 1) {
      new_index[i] <- current_index
      next
    }
    # 
    is_assigned <- FALSE
    for (j in 1:(i-1)) {
      if (is_equal(list_data[[i]], list_data[[j]])) {
        new_index[i] <- new_index[j]
        is_assigned <- TRUE
        break
      }
    }
    # 
    if (!is_assigned) {
      current_index <- current_index + 1
      new_index[i] <- current_index
    }
  }
  
  # 
  results_df <- data.frame(
    OldIndex = names(list_data),
    NewIndex = new_index
  )
  return(results_df)
}

#' @title mergeClones_0()
mergeClones_0 <- function(ratiodata,chrom = c(1:22),doPlot=FALSE,outdir="./"){
  suppressPackageStartupMessages({
    library(reshape)
    library(fpc)
    library(LaplacesDemon)
    library(diptest)
    library(mixtools)
    library(dplyr)
    library(tidyverse)
    library(rstatix)
    library(ggpubr)
  })
  chrom <- paste0("chr",chrom)
  set.seed(123)
  clonematrix <- matrix(0,nr=length(ratiodata),nc=length(ratiodata))
  for (chr in chrom){
    bindata <- genomicRegion(ratiodata,chr)
    chrratio <- cloneRatio(ratiodata,bindata,chr)
    merge.res <- genomicPartition(chrratio)
    while(sum(unlist(merge.res$clusterindex))!=0){
      clusterID <- unique(merge.res$data$cluster)
      cluster.res <- lapply(clusterID, function(cluster,merge.res){
        chrratio <- merge.res$data
        sub <- chrratio[chrratio$cluster==cluster,]
        if (merge.res$clusterindex[[cluster]] == 0){
          cluster.info <- list(mu=merge.res$mu[[cluster]],data=sub,clusterindex=0)
        }else{
          mu <- Modes(sub$ratio)$modes
          cluster.info <- regionSPlit(sub,mu)
        }
        return(cluster.info)
      },merge.res)
      merge.res <- regionMerge(cluster.res)
    }
    chrratio <- merge.res$data
    regioncount <- chrratio %>% dplyr::count(cluster,clone)
    regionID <- unique(regioncount$cluster[regioncount$n > 30])
    chrratio <- chrratio[chrratio$cluster %in% regionID,]
    if (dim(chrratio)[1] > 30){
      stat.test <- chrratio %>%
        group_by(cluster) %>%
        t_test(ratio ~ clone) %>%
        adjust_pvalue(method = "BH") %>%
        add_significance()
      row.index <- as.numeric(stat.test$group1[stat.test$p.adj < 0.01])
      col.index <- as.numeric(stat.test$group2[stat.test$p.adj < 0.01])
      cloneindex <- cbind(row.index,col.index)
      for (i in 1:length(row.index)){
        clonematrix[row.index[i],col.index[i]] <- 1
        clonematrix[col.index[i],row.index[i]] <- 1
      }
    }
  }
  clonemerge <- lapply(1:dim(clonematrix)[1], function(j,clonematrix){
    return(which(clonematrix[j,]==0))
  },clonematrix)
  cloneID <- 1:length(ratiodata)
  for (i in cloneID){
    if (length(clonemerge[[i]])>1){
      for (j in setdiff(clonemerge[[i]],i)){
        if (!(i %in% clonemerge[[j]])){
          clonemerge[[i]] <- setdiff(clonemerge[[i]],j)
        }
      }
    }
  }
  for (i in cloneID){
    if (length(clonemerge[[i]])>1){
      newclone <- clonemerge[[i]]
      for (j in setdiff(clonemerge[[i]],i)){
        newclone <- union(newclone,clonemerge[[j]])
      }
      for (j in setdiff(clonemerge[[i]],i)){
        clonemerge[[j]] <- newclone
      }
    }
  }
  clonemerge <- lapply(clonemerge, function(ele){
    return(ele-1)
  })
  names(clonemerge) <- names(ratiodata)
  if(doPlot){
    clonematrix <- matrix(0,nr=length(ratiodata),nc=length(ratiodata))
    row.names(clonematrix) <- names(ratiodata)
    colnames(clonematrix) <- names(ratiodata)
    for (i in 1:length(clonemerge)){
      clonematrix[i,clonemerge[[i]]+1] <- 1
    }
    breaks <- seq(0,1,by=0.2)
    colorss <- colorRampPalette(c("blue","red"))(length(breaks))
    pdf(file = paste0(outdir,"mergeclone.heatmap.pdf"),width = 4,height = 4)
    pheatmap::pheatmap(clonematrix,breaks = breaks,color = colorss)
    dev.off()
  }
  df_clone <- renewGroup(clonemerge)
  return(df_clone)
}

#' @title mergeClones()
#' @param ratiodata list of "step03.CNV_res_clonal.rds"
#' @param segdata list of clonal_RelativeRatio_segment
mergeClones <- function(ratiodata,segdata,doPlot=FALSE,outdir="./",Zscore.cutoff=1.28,
  p.adj.cutoff=0.1,p.cutoff=NULL){
  suppressPackageStartupMessages({
    library(reshape)
    library(fpc)
    library(LaplacesDemon)
    library(diptest)
    library(mixtools)
    library(dplyr)
    library(tidyverse)
    library(rstatix)
    library(ggpubr)
  })
  cloneID <- names(ratiodata)
  clonematrix <- matrix(0,nr=length(ratiodata),nc=length(ratiodata))
  colnames(clonematrix) <- cloneID
  rownames(clonematrix) <- cloneID
  for (clone in cloneID){
    subseg <- segdata[[clone]]
    otherclone <- names(subseg)
    otherclone <- otherclone[otherclone%in%cloneID]
    for (C2 in otherclone){
      subdata <- subseg[[C2]]
      subdata <- subdata[,c("binID","segName")]
      C1data <- ratiodata[[clone]]$input_BinSegRatio
      index1 <- match(subdata$binID,C1data$binID)
      subdata$C1 <- C1data$binRatio[index1]
      C2data <- ratiodata[[C2]]$input_BinSegRatio
      index1 <- match(subdata$binID,C2data$binID)
      subdata$C2 <- C2data$binRatio[index1]
      subdata.long <- suppressMessages(reshape2::melt(subdata, variable_name = "binID"))
      colnames(subdata.long) <- c("binID","segID","group","ratio")
      segID <- unique(subdata.long$segID)
      moderes <- do.call(rbind,lapply(segID, function(seg,subdata.long){
        subseg <- subdata.long[subdata.long$segID==seg,]
       
        if (dim(subseg)[1]>30){
          mu <- Modes(subseg$ratio)$modes
          return(c(seg,length(mu)))
        }
        #
       # dip.test(subseg$ratio)
        # subseg %>% t_test(ratio ~ group) %>%
        # adjust_pvalue(method = "BH") 
        # ggplot(subseg, aes(x = ratio, fill = group)) +
        #   geom_histogram(position = "dodge", bins = 20) 
        
      },subdata.long))
      moderes <- as.data.frame(moderes)
      colnames(moderes) <- c("segID","mode")
      index <- match(subdata.long$segID,moderes$segID)
      subdata.long <- subdata.long[!is.na(index),]

      ###
      groupMean <- subdata.long %>%
        group_by(segID,group)%>%
        summarise(groupMean = mean(ratio,na.rm=TRUE))%>% as.data.frame()
      diffMean <-  groupMean%>%
        group_by(segID)%>%
        summarise(diffMean =abs(groupMean[group=="C1"]- groupMean[group=="C2"]))%>% as.data.frame()
      ###
      group_counts <- subdata.long %>%
        group_by(segID, group) %>%
        dplyr::count(group) %>%as.data.frame()
      
      stat.test <- subdata.long %>%
        group_by(segID) %>%
        t_test(ratio ~ group) %>%
        adjust_pvalue(method = "BH") %>%
        add_significance()
      stat.test <- as.data.frame(stat.test)
      index <- match(moderes$segID,stat.test$segID)
      stat.test <- stat.test[index,]
      stat.test$mode <- moderes$mode
      ###
      stat.test <- left_join(stat.test,diffMean)
      z_scores <- (stat.test$diffMean - mean(stat.test$diffMean)) / sd(stat.test$diffMean)
      stat.test$diffMean.Zscore  <- z_scores
      if(!is.null(p.cutoff)){
        stat.test$mode[stat.test$p < p.cutoff & stat.test$diffMean.Zscore>Zscore.cutoff] <- 2
        if (dim(stat.test[stat.test$p < p.cutoff & stat.test$mode > 1,])[1] > 0){
          clonematrix[clone,C2] <- 1
        }
      }else{
          stat.test$mode[stat.test$p.adj < p.adj.cutoff & stat.test$diffMean.Zscore>Zscore.cutoff] <- 2 ### 90% quantile (Zscore=1.28)
          if (dim(stat.test[stat.test$p.adj < p.adj.cutoff & stat.test$mode > 1,])[1] > 0){
            clonematrix[clone,C2] <- 1
          }
      }
      
      ###
      
    }
  }
  clonematrix <- 1-clonematrix
  mergeRes <- lapply(cloneID, function(clone,clonematrix){
    colnames(clonematrix)[clonematrix[clone,]==1]
  },clonematrix)
  names(mergeRes) <- cloneID
  
  
  for (i in cloneID){
    if (length(mergeRes[[i]])>1){
      for (j in setdiff(mergeRes[[i]],i)){
        if (!(i %in% mergeRes[[j]])){
          mergeRes[[i]] <- setdiff(mergeRes[[i]],j)
        }
      }
    }
  }
  
  for (i in cloneID) {
    otherclone <- mergeRes[[i]]
    otherclone <- setdiff(otherclone,i)
    if (length(otherclone) > 1){
      ele <- unlist(lapply(otherclone, function(C2,mergeRes){
        if (length(intersect(setdiff(mergeRes[[C2]],C2),setdiff(otherclone,C2)))==0){
          return(C2)
        }
      },mergeRes))
      if (!is.null(ele)){
        mergeRes[[i]] <- setdiff(mergeRes[[i]],ele)
      }
    }
  }
  
  for (i in cloneID){
    if (length(mergeRes[[i]])>1){
      newclone <- mergeRes[[i]]
      for (j in setdiff(mergeRes[[i]],i)){
        newclone <- union(newclone,mergeRes[[j]])
      }
      for (j in newclone){
        mergeRes[[j]] <- newclone
      }
    }
  }
  
   

   if(doPlot){
    clonematrix <- matrix(0,nr=length(cloneID),nc=length(cloneID))
    row.names(clonematrix) <- cloneID
    colnames(clonematrix) <- cloneID
    for (i in cloneID){
      clonematrix[i,mergeRes[[i]]] <- 1
    }
    breaks <- seq(0,1,by=0.2)
    colorss <- colorRampPalette(c("blue","red"))(length(breaks))
    pdf(file = paste0(outdir,"/mergeclone.heatmap.newsegment.pdf"),width = 4,height = 4)
    pheatmap::pheatmap(clonematrix,breaks = breaks,color = colorss)
    dev.off()
  }
  df_clone <- renewGroup(mergeRes)
  return(df_clone)

}

#' @title refine_clones()
#' @description This function refines the clonal CNA estimates based on the paired-comparison of segment ratios.
#' @param ratiodata A list of data.frames containing the bin-level ratio estimates for each clone.
#' @param seg_compare_ls A list of data.frames containing the segment ratios for each clone compared to the best clone.
#' @param best_clone_ratiodata The name of the best clone for the segment ratio comparison in ratiodata.
refine_clones <- function(ratiodata,seg_compare_ls,best_clone_ratiodata,CNest.ref){
  #subseg= seg_ls[[1]]
  delta.ref <- (CNest.ref$ratio[2]-CNest.ref$ratio[1])/(CNest.ref$CN[2]-CNest.ref$CN[1])
  C1data <- ratiodata[[best_clone_ratiodata]]$input_BinSegRatio
  C1data <- C1data[,!grepl("relativeCN$|integerCN",colnames(C1data))]
  
  C1seg_dat <- ratiodata[[best_clone_ratiodata]]$seg.dat
  C1seg_dat <- C1seg_dat[,c("segName","relativeCN","integerCN")]
  C1data <- left_join(C1data,C1seg_dat,by="segName")
  rownames(C1data) <- C1data$binID
  
  #mode
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  otherclone <- names(seg_compare_ls)
  for (C2 in otherclone){
    C2data <- ratiodata[[C2]]$input_BinSegRatio
    C2data <- C2data[,!grepl("integerCN|relativeCN",colnames(C2data))]
    
    C2seg_dat<- ratiodata[[C2]]$seg.dat
    C2seg_dat2<- C2seg_dat[,c("segName","relativeCN","integerCN")]
    
    CNest <- unique(data.frame(ratio=C2seg_dat2$relativeCN,CN=C2seg_dat2$integerCN))
    CNest <- CNest[order(CNest$CN),]
    delta <- (CNest$ratio[2]-CNest$ratio[1])/(CNest$CN[2]-CNest$CN[1])
    
    C2data <- left_join(C2data,C2seg_dat2,by="segName")
    rownames(C2data) <- C2data$binID
    
    index2C1 <- match(C1data$binID,C2data$binID)
    C2data <- C2data[index2C1,] #same order with C1
    
    subdata <- seg_compare_ls[[C2]]
    subdata <- subdata[,c("binID","segName","SegMean")]
    index1 <- match(subdata$binID,C1data$binID)
    subdata$C1 <- C1data$binRatio[index1]
    index1 <- match(subdata$binID,C2data$binID)
    subdata$C2 <- C2data$binRatio[index1]
    subdata.long <- suppressMessages(melt(subdata,measure.vars=c("C1","C2"), variable_name = "binID"))
    colnames(subdata.long) <- c("binID","segID","SegMean","group","ratio")
    segID <- unique(subdata.long$segID)
    moderes <- do.call(rbind,lapply(segID, function(seg,subdata.long){
      subseg <- subdata.long[subdata.long$segID==seg,]
      
      if (dim(subseg)[1]>30){
        mu <- Modes(subseg$ratio)$modes
        return(c(seg,length(mu)))
      }
    },subdata.long))
    moderes <- as.data.frame(moderes)
    colnames(moderes) <- c("segID","mode")
    index <- match(subdata.long$segID,moderes$segID)
    subdata.long <- subdata.long[!is.na(index),]
    
    groupMean <- subdata.long %>%
      group_by(segID,group)%>%
      summarise(groupMean = mean(ratio,na.rm=TRUE))%>% as.data.frame()
    diffMean <-  groupMean%>%
      group_by(segID)%>%
      summarise(diffMean =abs(groupMean[group=="C1"]- groupMean[group=="C2"]))%>% as.data.frame()
    ###
    group_counts <- subdata.long %>%
      group_by(segID, group) %>%
      dplyr::count(group) %>%as.data.frame()
    
    stat.test <- subdata.long %>%
      group_by(segID) %>%
      t_test(ratio ~ group) %>%
      adjust_pvalue(method = "BH") %>%
      add_significance()
    stat.test <- as.data.frame(stat.test)
    index <- match(moderes$segID,stat.test$segID)
    stat.test <- stat.test[index,]
    stat.test$mode <- moderes$mode
    stat.test <- left_join(stat.test,diffMean)
    
    stat.test_res <- unique(stat.test[,c("segID","p.adj","mode","diffMean")])
    subdata <- left_join(subdata,stat.test_res,by=join_by(segName==segID))
    
    subdata$label <- 0
    lab_indx <- which(abs(subdata$SegMean)>(0.5*delta.ref)|(subdata$mode>1&subdata$p.adj<0.1))
    subdata$label[lab_indx] <- 1

    # subdata$label <- ifelse(abs(subdata$SegMean)<=(delta.ref-0.1),0,1)
    subdata_diff <- subdata[subdata$label==1,,drop=FALSE]
    
    subdata_diff <- subdata_diff %>%
      add_count(segName) %>%
      as.data.frame()
    colnames(subdata_diff)[colnames(subdata_diff)=="n"] <- "w"
    subdata_diff <- subdata_diff %>%
      dplyr::filter(w>2) %>%
      as.data.frame()
    subdata_diff_segs <- unique(subdata_diff$segName)
    
    if(length(subdata_diff_segs)>0){
      C2indx <- c()
      for(i in 1:length(subdata_diff_segs)){
        seg_i <- subdata_diff_segs[i]
        dif_bins_inC2 <- subdata$binID[subdata$segName %in% seg_i]
        C2data_diffSeg <- C2data[C2data$binID %in% dif_bins_inC2,]
        diffSegs <- unique(C2data_diffSeg$segName)
        
        row_indx <- which(C2data$segName %in% diffSegs)
        # relativeCN_new <- getmode(C2data_diffSeg$relativeCN)
        # integerCN_new <- getmode(C2data_diffSeg$integerCN)
        # SegMean_new <- getmode(C2data_diffSeg$SegMean)
        # 
        # 
        # 
        # C2data$relativeCN[row_indx] <- relativeCN_new
        # C2data$integerCN[row_indx] <- integerCN_new
        # C2data$SegMean[row_indx] <- SegMean_new
        
        C2indx <- c(C2indx,row_indx)
      }
    }
   
    C2indx_part2 <- setdiff(1:nrow(C2data),C2indx)

    integerCN_new2 <- C1data$integerCN[C2indx_part2]
    C2data$integerCN[C2indx_part2] <- integerCN_new2
    
    
    ratio_index <- match(integerCN_new2,CNest$CN)
    new_relativeCN <- CNest$ratio[ratio_index]
    
    if(any(is.na(ratio_index))){
      newCN_NA <- integerCN_new2[is.na(ratio_index)]
      newCN_NA <- newCN_NA[!is.na(newCN_NA)]
      if(length(newCN_NA)>0){
        newRatio <- (newCN_NA-CNest$CN[1])*delta + CNest$ratio[1]
        new_relativeCN[is.na(ratio_index)] <- newRatio
      }
    }
    
    C2data$relativeCN[C2indx_part2] <- new_relativeCN

    segs <- unique(C2data$segName)
    for(segj in segs){
      C2data_sub <- C2data[C2data$segName==segj,]
      row_indx2 <- which(C2data$segName %in% segj)
      relativeCN_new <- getmode(C2data_sub$relativeCN)
      integerCN_new <- getmode(C2data_sub$integerCN)
      SegMean_new <- getmode(C2data_sub$SegMean)
      C2data$relativeCN[row_indx2] <- relativeCN_new
      C2data$integerCN[row_indx2] <- integerCN_new
      C2data$SegMean[row_indx2] <- SegMean_new
    }
    C2data_seg_new <- na.omit(unique(C2data[,c("segName","relativeCN","integerCN")]))
    
    C2seg_dat <- C2seg_dat[,!grepl("relativeCN$|integerCN|seg_SD|seg_score",colnames(C2seg_dat))]
    C2seg_dat <- left_join(C2seg_dat,C2data_seg_new,by="segName")
    
    CNest <- unique(data.frame(ratio=C2seg_dat$relativeCN,CN=C2seg_dat$integerCN))
    CNest <- CNest[order(CNest$CN),]
    CNest <- peakIndex(C2data,CNest)
    ploidy <- sum(C2seg_dat$integerCN*C2seg_dat$w,na.rm=TRUE)/sum(C2seg_dat$w,na.rm=TRUE)
    
    ratiodata[[C2]]$seg.dat <- C2seg_dat
    ratiodata[[C2]]$input_BinSegRatio <- C2data
    ratiodata[[C2]]$ploidy <- ploidy
    ratiodata[[C2]]$CNest <- CNest
      
  }
  return(ratiodata)
    
}




