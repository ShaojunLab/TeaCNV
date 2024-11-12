
library(dplyr)

#' @title findccloneNum()
#' @description the optimal cluster number
#' @param CNdata: cell-bin integerCN matrix
#' @param n: the maximumal number of clusters
#' @export
#' 
findccloneNum <- function(CNdata,n=10){
  set.seed(123)
  bincount <- apply(CNdata, 2, function(x){length(x[is.na(x)])})
  CNdata <- CNdata[,bincount==0]
  wss <- numeric(n)
  for (i in 1:n) {
    # Fit the model: km.out
    km.out <- tryCatch(kmeans(CNdata, centers = i, nstart = 20),error=function(e) NA )
    if(all(is.na(km.out))){
      wss[i] <-NA
    }else{
      # Save the within cluster sum of squares
      wss[i] <- km.out$tot.withinss
    }
  }
  wss <- wss[!is.na(wss)]
  diff <- wss[2:length(wss)]-wss[1:(length(wss)-1)]
  d <- diff[2:length(diff)]-diff[1:(length(diff)-1)]
  knum <- NULL
  for (i in 2:(length(d)-1)){
    if(i>=2){
      if (abs(d[i]) < abs(d[i-1])/10){
        knum <- c(knum,i)
      }
    }
  }
  if (!is.null(knum)){
    if (length(knum)==1){
      return(knum+1)
    }else{
      knum <- knum[d[knum]>0]
      if (length(knum) > 0){
        return(knum[which.min(d[knum])]+1)
      }else{
        kindex <- which(d > 0)
        return(which.min(d[kindex])+1)
      }
    }
  }else{
    kindex <- which(d > 0)
    return(which.min(d[kindex])+1)
  }
}



#' @title cellRatio_corrected()
#' @description cell level corrected ratio
#' @param cellratioID: initial cell ratio
#' @param finalres: clonal estimation
#' @return cell-bin matrix: corrected Ratio
#' @export

cellRatio_corrected <- function(cellratioID,finalres,refClone=NULL){
  cellres <- lapply(names(cellratioID), function(cluster,cellratioID,finalres,refClone){
    celldata <- cellratioID[[cluster]]
    if(!is.null(refClone)){
      clonalest.Ref <- finalres[[refClone]]
    }else{
      clonalest.Ref <- NULL
    }
    
    if(cluster %in% names(finalres)){
      clonalest <- finalres[[cluster]]
    }else if(!is.null(refClone)){
      clonalest <- finalres[[refClone]]
    }else{
      clonalest <- finalres[[1]]
    }
    segratio.correct <- cellPloidyCorrection(celldata,clonalest,clonalest.Ref)
    return(segratio.correct)
  },cellratioID,finalres,refClone)
  names(cellres)  <- names(cellratioID)
  
  binID <- do.call(rbind,lapply(names(cellratioID), function(cluster,cellratioID){
    celldata <- cellratioID[[cluster]]
    binID <- do.call(rbind,(lapply(1:length(celldata), function(j,celldata){
      return(cbind(celldata[[j]]$binID,celldata[[j]]$segName))
    },celldata)))
    binID <- unique(binID)
  },cellratioID))
  
  binID <- as.data.frame(binID)
  colnames(binID) <- c("binID","segID")
  binID <- unique(binID)
  binpos <- do.call(rbind,strsplit(binID$binID,split="-"))
  binpos <- cbind(binpos,binID$segID)
  chromo <- unique(binpos[,1])
  bin <- do.call(rbind,lapply(chromo, function(chr,binpos){
    subdata <- binpos[binpos[,1]==chr,]
    bindata <- data.frame(chr=subdata[,1],start=as.numeric(subdata[,2]),end=as.numeric(subdata[,3]),segID=subdata[,4])
    bindata <- bindata[order(bindata[,2]),]
    return(bindata)
  },binpos))
  bin$binID <- paste(bin$chr,bin$start,bin$end,sep="-")
  binID <- unique(bin$binID)
  
  cellbinRatio <- do.call(rbind,lapply(1:length(cellres), function(j,cellres,bin,binID){
    res <- cellres[[j]]
    matrix(NA,nc=length(binID),nr=dim(res)[1])
    
    segID <- colnames(res)
    binRatio <- matrix(NA,nc=length(binID),nr=dim(res)[1])
    index1 <- match(bin$segID,segID)
    bindata <- res[,index1[!is.na(index1)]]
    colnames(bindata) <- bin$binID[!is.na(index1)]
    index <- match(colnames(bindata),binID)
    binRatio[,index]=bindata
    row.names(binRatio) <- row.names(bindata)
    colnames(binRatio) <- binID
    return(binRatio)
  },cellres,bin,binID))
  return(cellbinRatio)
}



#' @title cellintegerCN()
#' @description cell level integerCN estimation
#' @param cellratioID: cell ratio
#' @param finalres: clonal estimation
#' @return cell-bin matrix: integerCN
#' @export

cellintegerCN <- function(cellratioID,finalres,refClone=NULL){
  cellres <- lapply(names(cellratioID), function(cluster,cellratioID,finalres,refClone){
    celldata <- cellratioID[[cluster]]
    if(!is.null(refClone)){
      clonalest.Ref <- finalres[[refClone]]
    }else{
      clonalest.Ref <- NULL
    }
    
    if(cluster %in% names(finalres)){
      clonalest <- finalres[[cluster]]
    }else if(!is.null(refClone)){
      clonalest <- finalres[[refClone]]
    }else{
      clonalest <- finalres[[1]]
    }
    segratio.correct <- cellPloidyCorrection(celldata,clonalest,clonalest.Ref)
    cellCN <- cellEst(segratio.correct,clonalest,clonalest.Ref)
    return(cellCN)
  },cellratioID,finalres,refClone)
  names(cellres)  <- names(cellratioID)
  
  binID <- do.call(rbind,lapply(names(cellratioID), function(cluster,cellratioID){
    celldata <- cellratioID[[cluster]]
    binID <- do.call(rbind,(lapply(1:length(celldata), function(j,celldata){
      return(cbind(celldata[[j]]$binID,celldata[[j]]$segName))
    },celldata)))
    binID <- unique(binID)
  },cellratioID))
  
  binID <- as.data.frame(binID)
  colnames(binID) <- c("binID","segID")
  binID <- unique(binID)
  binpos <- do.call(rbind,strsplit(binID$binID,split="-"))
  binpos <- cbind(binpos,binID$segID)
  chromo <- unique(binpos[,1])
  bin <- do.call(rbind,lapply(chromo, function(chr,binpos){
    subdata <- binpos[binpos[,1]==chr,]
    bindata <- data.frame(chr=subdata[,1],start=as.numeric(subdata[,2]),end=as.numeric(subdata[,3]),segID=subdata[,4])
    bindata <- bindata[order(bindata[,2]),]
    return(bindata)
  },binpos))
  bin$binID <- paste(bin$chr,bin$start,bin$end,sep="-")
  binID <- unique(bin$binID)
  
  cellbinCN <- do.call(rbind,lapply(1:length(cellres), function(j,cellres,bin,binID){
    res <- cellres[[j]]
    matrix(NA,nc=length(binID),nr=dim(res)[1])
   
    segID <- colnames(res)
    binCN <- matrix(NA,nc=length(binID),nr=dim(res)[1])
    index1 <- match(bin$segID,segID)
    bindata <- res[,index1[!is.na(index1)]]
    colnames(bindata) <- bin$binID[!is.na(index1)]
    index <- match(colnames(bindata),binID)
    binCN[,index]=bindata
    row.names(binCN) <- row.names(bindata)
    colnames(binCN) <- binID
    return(binCN)
  },cellres,bin,binID))
  return(cellbinCN)
}



#' @title cellPloidyCorrection()
#' @description centering cell ratio both segment and cell level based on the clonal expectation
#' @param celldata cell level ratio
#' @param clonalest: clonal ploidy estimation
#' @param clonalest.Ref reference clone if the user specifies
#' @return corrected ratio matrix
#' @export


cellPloidyCorrection <- function(celldata,clonalest,clonalest.Ref=NULL){
  if(!is.null(clonalest.Ref)){
    clonal.seg <- clonalest.Ref$seg.dat
  }else{
    clonal.seg <- clonalest$seg.dat
  }
  
  segID <- clonalest$seg.dat$segName
  
  CNest <- as.data.frame(clonal.seg %>% group_by(integerCN) %>% summarise_at(vars(sd), list(name = mean)))
  CNest$ratio <- clonal.seg$relativeCN[match(CNest[,1],clonal.seg$integerCN)]
  names(CNest) <- c("integerCN","sd","ratio")
  overallratio <- sum(clonal.seg$relativeCN*clonal.seg$w,na.rm = TRUE)/sum(clonal.seg$w,na.rm = TRUE)
  ###cell ploidy centering
  segploidy.correct <- do.call(rbind,lapply(1:length(celldata), function(j,celldata,clonal.seg,overallratio){
    subcell <- celldata[[j]]
    cellseg <- unique(subcell[,c("SegMean","segName")])
    
    index <- match(cellseg$segName,clonal.seg$segName)
    #if segname do not match
    if((!is.null(clonalest.Ref)) & any(is.na(index))){
      Nseg <- as.data.frame(table(celldata[[j]]$segName))
      colnames(Nseg) <- c("segName","w")
      cellseg <- suppressMessages(left_join(cellseg,Nseg))
    }else{
      cellseg$w <- clonal.seg$w[index]
    }
    
    ratio <- sum(cellseg$SegMean*cellseg$w,na.rm = TRUE)/sum(cellseg$w,na.rm = TRUE)
    bias <- overallratio - ratio
    cellseg$newratio <- cellseg$SegMean+bias
    
    index <- match(segID,cellseg$segName)
    return(cellseg$newratio[index])
  },celldata,clonal.seg,overallratio))
  row.names(segploidy.correct) <- names(celldata)
  colnames(segploidy.correct) <- segID
  ###segment centering
  cellmean <- apply(segploidy.correct, 2, function(x){return(mean(x[!is.na(x)]))})
  delt <- clonalest$seg.dat$relativeCN-cellmean
  segratio.correct <- apply(segploidy.correct, 1, function(x,delt){return(x+delt)},delt)
  return(segploidy.correct)
}



#' @title cellEst()
#' @description integer CN estimation at cell level
#' @param segratio.correct: corrected ratio matrix
#' @param clonalest: clonal ploidy estimation
#' @return cell-integerCN matrix
#' @export

cellEst <- function(segratio.correct,clonalest,clonalest.Ref=NULL){
  if(!is.null(clonalest.Ref)){
    clonal.seg <- clonalest.Ref$seg.dat
  }else{
    clonal.seg <- clonalest$seg.dat
  }
  segID <- clonalest$seg.dat$segName
  CNest <- as.data.frame(clonal.seg %>% group_by(integerCN) %>% summarise_at(vars(sd), list(name = mean)))
  CNest$ratio <- clonal.seg$relativeCN[match(CNest[,1],clonal.seg$integerCN)]
  names(CNest) <- c("integerCN","sd","ratio")
  cellCN <- apply(segratio.correct,2,function(x,CNest){
    pp <- apply(CNest, 1, function(y,x){
      p <- dnorm(x,y[3],y[2])
      return(p)
    },x)
    colnames(pp) <- CNest$integerCN
    pprob <- apply(pp, 1, function(p){
      return(p/sum(p))
    })
    pprob <- t(pprob)
    integerCNindex <- apply(pprob, 1, function(x){
      if (sum(is.na(x))==length(x)){
        return(NA)
      }else{
        return(which.max(x))
      }
    })
    cellCN <- CNest$integerCN[integerCNindex]
    return(cellCN)
  },CNest)
  colnames(cellCN) <- segID
  row.names(cellCN) <- row.names(segratio.correct)
  return(cellCN)
} 


#' @title minPloidysummary()
#' @export
minPloidysummary <- function(clusterout){
  delt_set <- do.call(cbind,lapply(1:length(clusterout), function(j,clusterout){
    if (!is.null(clusterout[[j]])){
      res <- clusterout[[j]]
      CNest <- res$newCN
      if(nrow(CNest)>1){
        delt <- (CNest$ratio[2]-CNest$ratio[1])/(CNest$CN[2]-CNest$CN[1])
        return(delt)
      }
    }
    
  },clusterout))
  delt_mean <- median(delt_set)
  
  clustersta <- do.call(rbind,lapply(1:length(clusterout), function(j,clusterout,delt_mean){
    if (!is.null(clusterout[[j]])){
      res <- clusterout[[j]]
      #input_BinSegRatio <- res$input_BinSegRatio
      seg.dat <- res$seg.dat
      CNest <- res$newCN
      baseratio <- CNest$ratio[1]
      baseCN <- CNest$CN[1]
      if(nrow(CNest)>1){
        delt <- (CNest$ratio[2]-CNest$ratio[1])/(CNest$CN[2]-CNest$CN[1])
      }else{
        baseCN <- 2
        delt <- delt_mean
      }
      
      #seg.dat1 <- seg.dat[,-grep("CNV",colnames(seg.dat))]
      integerCN <- round((seg.dat$ratio-baseratio)/delt)+baseCN
      seg.dat$CN <- integerCN
      CNcount <- as.data.frame(aggregate(seg.dat$w, by=list(CN=seg.dat$CN), FUN=sum))
      CNcount$frac <- CNcount$x/sum(CNcount$x)
      CNstate <- unique(integerCN)
      CNstate <- CNstate[order(CNstate)]
      ratio <- baseratio+(CNstate-baseCN)*delt
      CNest <- data.frame(ratio=ratio,CN=CNstate)
      CN_dominant <- CNcount$CN[which.max(CNcount$frac)]
      ratio_dominant  <- CNest$ratio[CNest$CN==CN_dominant]
      return(c(j,baseratio,delt,min(CNcount$CN[CNcount$frac>0.001]),CN_dominant,ratio_dominant))
    }
  },clusterout,delt_mean))
  return(clustersta)
}


#' @title ploidyRefine()
#' @description the final estimation : select the optimal ploidy estimation for each cluster and correct the negative integer CN
#' @param sampleres initial estimation list including all clusters:
#'                                                        $CNVestimate: initial CNV estimation parameters
#'                                                        $seg.dat: segmentatio data
#'                                                        $input_BinSegRatio: bin data
#' @return for the input cluster output: $input_BinSegRatio
#'                                        $seg.dat
#'                                        $CNest: ratio corresponding to integerCN
#'                                        $ploidy: the overall ploidy
#' @export

ploidyRefine <- function(sampelres,delt.lim = 0.3,minCN.frac=0.01){
  #sampelres <- initialRes[[ID]]
  clusterout <- lapply(names(sampelres), function(cluster,sampelres){
    if ("seg.dat" %in% names(sampelres[[cluster]])){
      seg.dat <- sampelres[[cluster]]$seg.dat
      df_seg_C1 <- sampelres[[cluster]]$input_BinSegRatio
      integerCNV <- sampelres[[cluster]]
      clusterEst <- PloidyCorrect(integerCNV,delt.lim = delt.lim)
      
      if (length(clusterEst) > 1){
        likEst <- optimalPloidy(clusterEst)
        if(all(!is.na(likEst$AIC))){
          kindex <- likEst$est[which.min(likEst$AIC)]
        }else{
          kindex <- 1
        }
        output <- clusterEst[[kindex]]
        output$integerCN$base=0
        output$integerCN$base[output$integerCN$ratio==output$base]=1
      }else{
        output <- clusterEst[[1]]
        output$integerCN$base=0
        output$integerCN$base[output$integerCN$ratio==output$base]=1
      }
      return(output)
    }
  },sampelres)
  
  baseratio <- unlist(lapply(1:length(clusterout), function(j,clusterout){
    a <- clusterout[[j]]$integerCN$ratio[clusterout[[j]]$integerCN$base==1]
    return(a)
  },clusterout))
  
  baseCN <- unique(round(baseratio*2))
  if (length(baseCN) > 1){
    getmode <- function(v) {
      uniqv <- unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    baseCN <- getmode(round(baseratio*2))
  }
  CNchange <- NULL
  while (1) {
    if (!is.null(CNchange)){
      break
    }else{
      baseindex <- 0
      for (i in 1:length(clusterout)){
        #print(i)
        if ("integerCN" %in% names(clusterout[[i]])){
          CNest <- clusterout[[i]]$integerCN
          k <- which(CNest$base==1)
          newCN <- data.frame(ratio=CNest$ratio,CN=c(1:dim(CNest)[1])-k+baseCN)
          clusterout[[i]]$newCN <- newCN
          seg.dat <- clusterout[[i]]$seg.dat
          if(nrow(newCN)>1){
            delt <- newCN$ratio[2]-newCN$ratio[1]
            ratio0 <- newCN$ratio[k]-baseCN*delt
            dis <- round((seg.dat$ratio-ratio0)/delt)
          }else{
            dis=0
          }

          if (min(dis) < 0){
            posindex <- which(dis < 0)
            seg.len <- seg.dat$w[posindex]
            if (sum(seg.len > 10) == length(seg.len)){
              baseindex <- 1
            }
          }
        }
      }
      if (baseindex == 1){
        baseCN <- baseCN + 1
      }else{
        CNchange = 1
      }
    }
  }
  
  outcorrect <- ploidyUpdate(clusterout,baseCN,delt.lim = delt.lim)
  clusterout <- outcorrect$clusterout
  clusterout_summary <- minPloidysummary(clusterout)
  baseCN <- outcorrect$baseCN
  finalres <- lapply(1:length(clusterout), function(j,clusterout,baseCN){
    if (!is.null(clusterout[[j]])){
      res <- clusterout[[j]]
      input_BinSegRatio <- res$input_BinSegRatio
      seg.dat <- res$seg.dat
      CNest <- res$newCN
      CNest <- peakIndex(input_BinSegRatio,initialCN = CNest)
      baseCN_j <- CNest$CN[CNest$base==1]
      
      if (baseCN == baseCN_j){
        baseratio <- CNest$ratio[CNest$CN==baseCN]
      }else if(sum(CNest$base)==1){
        baseratio <- CNest$ratio[CNest$CN==baseCN_j]
       # baseCN_j <- baseCN
      }else{
        baseratio <- CNest$ratio[1]
        baseCN_j <- CNest$CN[1]
      }
      if(nrow(CNest)>1){
        delt <- (CNest$ratio[2]-CNest$ratio[1])/(CNest$CN[2]-CNest$CN[1])
      }else{
        delt <- clusterout_summary[j,3]
      }
      
      if(delt < delt.lim){
        delt_original<- delt
        while(delt < delt.lim){
          delt <- delt + delt_original
        }
      }else if(delt > 2*delt.lim){
        delt <- delt/2
        baseratio <- CNest$ratio[1]
        baseCN_j <- CNest$CN[1]
      }
      seg.dat1 <- seg.dat[,-grep("CNV",colnames(seg.dat))]
      integerCN <- round((seg.dat$ratio-baseratio)/delt)+baseCN_j
      integerCN[integerCN < 0] = 0
      
      CNstate <- unique(integerCN)
      CNstate <- CNstate[order(CNstate)]
      ratio <- baseratio+(CNstate-baseCN_j)*delt
      CNestnew <- data.frame(ratio=ratio,CN=CNstate)
      CNestnew <- CNestnew[CNestnew$CN>=0,]
      CNestnew <- peakIndex(input_BinSegRatio,initialCN = CNestnew)
      
      relativeCN <- baseratio+(integerCN-baseCN_j)*delt
      seg.dat1$relativeCN <- relativeCN
      seg.dat1$integerCN <- integerCN
      
      #edited 0514
      CNcount <- as.data.frame(aggregate(seg.dat1$w, by=list(ratio=seg.dat1$relativeCN,CN=seg.dat1$integerCN), FUN=sum))
      CNcount$frac <- CNcount$x/sum(CNcount$x)
      CNcount <- CNcount[CNcount$frac > minCN.frac,]
      m <- min(CNcount$CN)
      if( nrow(CNcount)==1){
        CNestnew$CN <-2
        CNcount$CN <- 2
      }else if (m <=0 & nrow(CNcount)>1 ){
        CNcount$CN <- CNcount$CN+(1-m)
        CNestnew$CN <- CNestnew$CN+(1-m)
      }else{
        if (m >1 & CNcount$CN[which.max(CNcount$frac)]>2){
          CNestnew$CN <- CNestnew$CN-(m-1)
        }else if (CNcount$CN[which.max(CNcount$frac)]==1){
          CNestnew$CN <- CNestnew$CN+1
          CNcount$CN <- CNcount$CN+1
        }
      }
      
      baseCN_j_new <- CNestnew$CN[CNestnew$base==1]
      baseratio_new <- CNestnew$ratio[CNestnew$base==1]
      integerCN <- round((seg.dat$ratio-baseratio_new)/delt)+baseCN_j_new
      integerCN[integerCN < 0] = 0
      relativeCN <- baseratio_new+(integerCN-baseCN_j_new)*delt
      seg.dat1$relativeCN <- relativeCN
      seg.dat1$integerCN <- integerCN
      #unique(seg.dat1[,c("relativeCN","integerCN")])
      
      ploidy <- sum(seg.dat1$integerCN*seg.dat1$w,na.rm=TRUE)/sum(seg.dat1$w,na.rm=TRUE)
      a <- list(input_BinSegRatio = input_BinSegRatio, seg.dat = seg.dat1, CNest = CNestnew, ploidy = ploidy)
    }else{
      a <- NULL
    }
    return(a)
  },clusterout,baseCN)
  names(finalres) <- names(sampelres)
  
  listnames <- c()
  res <- list()
  for (i in 1:length(finalres)){
    if (!is.null(finalres[[i]])){
      if (!is.null(finalres[[i]]$input_BinSegRatio)){
        res[[length(res)+1]]=finalres[[i]]
        listnames <- c(listnames,names(finalres)[i])
      }
    }
  }
  names(res) <- listnames
  
  return(res)
  
}




#' @title refine_segCN.Bayes()
#' @description This function refines the segCN results using a Bayesian classifier.
#' @param clonalres A list of clonal segCN results.
#' @param seg_dat_ref A data frame of segmented data with columns "ratio" and "integerCN".
#' @return A list of refined segCN results.
#install.packages("e1071")
#' @export

#library(e1071)
refine_segCN.Bayes <- function(clonalres,seg_dat_ref) {
  library(e1071)
  #naiveBayes
  Nsample <- 1:nrow(seg_dat_ref)
  diploidy_index <- which(seg_dat_ref$integerCN==2)
  aneuploidy_indices <- setdiff(Nsample, diploidy_index)
  set.seed(123)
  train_indices <- sample(diploidy_index, ceiling(length(diploidy_index)*0.8))
  
  test_indices <- setdiff(Nsample, train_indices)
  train_indices <- sort(c(train_indices,aneuploidy_indices))
  
  
  train_data <- seg_dat_ref[train_indices,"ratio",drop=FALSE]
  train_labels <- factor(seg_dat_ref[train_indices,"integerCN"])
  test_data <- seg_dat_ref[test_indices,"ratio",drop=FALSE]
  test_labels <- factor(seg_dat_ref[test_indices,"integerCN"])
  model <- naiveBayes(train_data, train_labels)
  # predictions <- predict(model, test_data)
  #print(table(predictions, test_labels))
  # Visualize the results
  # plot(seg_dat_ref$ratio, col=factor(seg_dat_ref$integerCN), pch=19, main="Bayesian Classifier")
  # points(test_data, col=predictions, pch=8, cex=2)
  # legend("topleft", legend=c("Class 1", "Class 2"), col=1:2, pch=19)
  # legend("topright", legend=c("Predicted Class 1", "Predicted Class 2"), col=1:2, pch=8)
  
  
  clusterout <- lapply(1:length(clonalres), function(j,clonalres){
    if (!is.null(clonalres[[j]])){
      res <- clonalres[[j]]
      seg_dat <- res$seg.dat
      colnames(seg_dat) <- gsub("^CNV$","relativeCN",colnames(seg_dat))
      colnames(seg_dat) <- gsub("^integerCNV$","integerCN",colnames(seg_dat))
      
      CNest <- unique(data.frame(ratio=seg_dat$relativeCN,CN=seg_dat$integerCN))
      CNest <- CNest[order(CNest$CN),]
      
      delta <- (CNest$ratio[2]-CNest$ratio[1])/(CNest$CN[2]-CNest$CN[1])
      
      Prediction_set <- seg_dat[,"ratio",drop=FALSE]
      predictions <- predict(model, Prediction_set)
      seg_dat$integerCN_predict <- predictions
      diff_row <- which(seg_dat$integerCN!=seg_dat$integerCN_predict)
      newCN <- as.numeric(seg_dat$integerCN_predict[diff_row])
      ratio_index <- match(newCN,CNest$CN)
      new_relativeCN <- CNest$ratio[ratio_index]
      
      if(any(is.na(ratio_index))){
        newCN_NA <- newCN[is.na(ratio_index)]
        newRatio <- (newCN_NA-CNest$CN[1])*delta + CNest$ratio[1]
        new_relativeCN[is.na(ratio_index)] <- newRatio
      }
      
      seg_dat$relativeCN[diff_row] <- new_relativeCN
      seg_dat$integerCN[diff_row] <- newCN
      seg_dat <- seg_dat[,!grepl("integerCN_predict",colnames(seg_dat)),drop=FALSE]
      CNest <- unique(data.frame(ratio=seg_dat$relativeCN,CN=seg_dat$integerCN))
      CNest <- CNest[order(CNest$CN),]
      CNest <- peakIndex(res$input_BinSegRatio,CNest)

      ploidy <- sum(seg_dat$integerCN*seg_dat$w,na.rm=TRUE)/sum(seg_dat$w,na.rm=TRUE)
      res$ploidy <- ploidy
      res$CNest <- CNest
      res$seg.dat <- seg_dat
      
      clonalres[[j]] <- res
    }else{
      clonalres[[j]] <- NULL
    }
    return(clonalres[[j]])
  },clonalres)
  names(clusterout) <- names(clonalres)
  return(clusterout)
}


#' @title refine_segCN.dist()

#' @description This function refines segCN based on distance to neighbour segments.
#' @param clonalres A list of clonal segCN results.
#' @return A list of refined segCN results.
#' @export
refine_segCN.dist <- function(clonalres) {
  clusterout <- lapply(1:length(clonalres), function(j,clonalres){
    if (!is.null(clonalres[[j]])){
      res <- clonalres[[j]]
      seg_dat <- res$seg.dat
      colnames(seg_dat) <- gsub("^CNV$","relativeCN",colnames(seg_dat))
      colnames(seg_dat) <- gsub("^integerCNV$","integerCN",colnames(seg_dat))
      
      CNest <- unique(data.frame(ratio=seg_dat$relativeCN,CN=seg_dat$integerCN))
      CNest <- CNest[order(CNest$CN),]
      
      delta <- (CNest$ratio[2]-CNest$ratio[1])/(CNest$CN[2]-CNest$CN[1])
      
      for(i in 1:nrow(seg_dat)){
        segi <- seg_dat[i,"segName"]
        segi_Len <- as.numeric(seg_dat[i,"end"])-as.numeric(seg_dat[i,"start"])
        # chri <- seg_dat[i,"chr"]
        ratio_i <- seg_dat[i,"ratio"]
        CN_i <- seg_dat[i,"integerCN"]
        k <- seq(i-2,i+2)
        k <- k[k>=1 & k<=nrow(seg_dat)&k!=i]
        # chrj <- seg_dat[j,"chr"]
        # j <- j[chrj==chri]
        neighbour_segs <- seg_dat[k,"segName"]
        neighbour_segsLen <- as.numeric(seg_dat[k,"end"])-as.numeric(seg_dat[k,"start"])
        ratio_neighbours <- seg_dat[k,"ratio"]
        CN_neighbours <- seg_dat[k,"integerCN"]
        relativeCN_neighbours <- seg_dat[k,"relativeCN"]

        dist_neighbours <- (ratio_i-ratio_neighbours)
        dist_CNratio <- (ratio_i-seg_dat[i,"relativeCN"])
        
        if(any(abs(dist_neighbours)<abs(dist_CNratio))|abs(dist_CNratio)>=delta){
          min_index <- which.min(abs(dist_neighbours))
          if(CN_i != CN_neighbours[min_index] ){
           # if((dist_CNratio < 0 & dist_neighbours[min_index]>0) |(dist_CNratio > 0 & dist_neighbours[min_index]>0)
           if(length(unique(CN_neighbours))==1 & any(segi_Len <= neighbour_segsLen) ){
             seg_dat[i,"integerCN"] <- CN_neighbours[min_index]
             seg_dat[i,"relativeCN"] <- relativeCN_neighbours[min_index]
           }else if( CN_i !=2 & segi_Len <= neighbour_segsLen[min_index] ){
             seg_dat[i,"integerCN"] <- CN_neighbours[min_index]
             seg_dat[i,"relativeCN"] <- relativeCN_neighbours[min_index]
           }
          }
        }
      }
      CNest <- unique(data.frame(ratio=seg_dat$relativeCN,CN=seg_dat$integerCN))
      CNest <- CNest[order(CNest$CN),]
      CNest <- peakIndex(res$input_BinSegRatio,CNest)
      
      ploidy <- sum(seg_dat$integerCN*seg_dat$w,na.rm=TRUE)/sum(seg_dat$w,na.rm=TRUE)
      res$ploidy <- ploidy
      res$CNest <- CNest
      res$seg.dat <- seg_dat
      
      clonalres[[j]] <- res
    }else{
      clonalres[[j]] <- NULL
    }
    return(clonalres[[j]])
  },clonalres)
  names(clusterout) <- names(clonalres)
  return(clusterout)
}


#' @title refine_segCN.dist_v2()
#' @export
refine_segCN.dist_v2 <- function(clonalres,Nadjacent=2) {
  clusterout <- lapply(1:length(clonalres), function(j,clonalres){
    if (!is.null(clonalres[[j]])){
      #print(paste0("clone",j))
      res <- clonalres[[j]]
      seg_dat <- res$seg.dat
      colnames(seg_dat) <- gsub("^CNV$","relativeCN",colnames(seg_dat))
      colnames(seg_dat) <- gsub("^integerCNV$","integerCN",colnames(seg_dat))
      seg_dat <- na.omit(seg_dat)
      
      CNest <- unique(data.frame(ratio=seg_dat$relativeCN,CN=seg_dat$integerCN))
      CNest <- CNest[order(CNest$CN),]
      SegRatioSD <- seg_dat %>%
        group_by(integerCN)%>%
        summarise(Rsd=sd(ratio))
      SegRatioSD <- na.omit(SegRatioSD)
      CNest <- left_join(CNest,SegRatioSD,by=c("CN"="integerCN"))%>%as.data.frame()
      CNest <- na.omit(CNest)
      # diploidy_lo <- round(min(seg_dat$ratio[seg_dat$integerCN==2]),2)
      # diploidy_hi <- round(max(seg_dat$ratio[seg_dat$integerCN==2]),2)
      diploidy_lo <- round(median(seg_dat$ratio[seg_dat$integerCN==2],na.rm=T)-1.5*SegRatioSD$Rsd[SegRatioSD$integerCN==2],2)
      diploidy_hi <- round(median(seg_dat$ratio[seg_dat$integerCN==2],na.rm=T)+1.5*SegRatioSD$Rsd[SegRatioSD$integerCN==2],2)
      
      delta <- (CNest$ratio[2]-CNest$ratio[1])/(CNest$CN[2]-CNest$CN[1])
      bin_dat <- res$input_BinSegRatio
      bin_dat <- bin_dat[,!grepl("integerCNV|integerCN|relativeCN",colnames(bin_dat)),drop=F]
      bin_dat <-left_join(bin_dat,seg_dat[,c("segName","relativeCN","integerCN")],by="segName")
      binseg_dat <- unique(bin_dat[,c("segName","relativeCN","integerCN","SegMean")])
      
      #Correct False aneuploidy
      for(i in 1:nrow(seg_dat)){
        #print(i)
        segi <- seg_dat[i,"segName"]
       
        
        segi_Len <- as.numeric(seg_dat[i,"end"])-as.numeric(seg_dat[i,"start"])
        w_i <- seg_dat[i,"w"]
        #chri <- seg_dat[i,"chr"]
        #ratio_i <- seg_dat[i,"ratio"]
        ratio_i <- unique(bin_dat$SegMean[bin_dat$segName %in% segi])
        CN_i <- seg_dat[i,"integerCN"]
        k <- seq(i-Nadjacent,i+Nadjacent)
        k <- k[k>=1 & k<=nrow(seg_dat)&k!=i]
        #chrj <- seg_dat[k,"chr"]
        #k <- k[chrj==chri]
        neighbour_segs <- seg_dat[k,"segName"]
        neighbour_w <- seg_dat$w[k]
        neighbour_segsLen <- as.numeric(seg_dat[k,"end"])-as.numeric(seg_dat[k,"start"])
        ratio_neighbours <- binseg_dat$SegMean[match(neighbour_segs,binseg_dat$segName)]
        CN_neighbours <- seg_dat[k,"integerCN"]
        relativeCN_neighbours <- seg_dat[k,"relativeCN"]
        
        dist_neighbours <- (ratio_i-ratio_neighbours)
        dist_CNratio <- (ratio_i-seg_dat[i,"relativeCN"])
        min_index <- which.min(abs(dist_neighbours))
        
        # bins_i <- bin_dat$binID[bin_dat$segName %in% c(segi,neighbour_segs[min_index])]
        # #Between-group
        # #otherGroups <- setdiff(1:length(clonalres),j)
        # otherGroups <- 1:length(clonalres)
        # 
        # neighbour_segs_groups <- do.call(rbind,lapply(otherGroups,function(x,clonalres,bins_i){
        #   #print(x)
        #   datk <- clonalres[[x]]
        #   seg_gr <- unique(datk$input_BinSegRatio$segName[datk$input_BinSegRatio$binID %in% bins_i])
        #   seg_gr <- datk$seg.dat$segName[datk$seg.dat$segName %in% seg_gr]
        #   seg_gr_CN<- datk$seg.dat$integerCN[datk$seg.dat$segName %in% seg_gr]
        #   seg_gr_w <- datk$seg.dat$w[datk$seg.dat$segName %in% seg_gr]
        #   relativeCN <-  datk$seg.dat$relativeCN[datk$seg.dat$segName %in% seg_gr]
        #   length <- as.numeric(datk$seg.dat$end[datk$seg.dat$segName %in% seg_gr])-as.numeric(datk$seg.dat$start[datk$seg.dat$segName %in% seg_gr])
        #   df <- data.frame(segName=seg_gr,integerCN=seg_gr_CN,w=seg_gr_w,relativeCN=relativeCN,length=length)
        #   df$clone <- names(clonalres)[[x]]
        #   return(df)
        # },clonalres,bins_i))
        # 
        # longSegIndx <- which(neighbour_segs_groups$w>=w_i)
        # if(length(unique(neighbour_segs_groups$integerCN[longSegIndx]))>2){
        #   mainCN <- round(mean(neighbour_segs_groups$integerCN[longSegIndx],na.rm=T))
        # }else if(length(longSegIndx)==1){
        #   w_cutoff <- quantile(neighbour_segs_groups$w,0.6)
        #   mainCN <- round(mean(neighbour_segs_groups$integerCN[neighbour_segs_groups$w>=w_cutoff],na.rm=T))
        # }else if(length(unique(neighbour_segs_groups$integerCN[longSegIndx]))==2 & CN_i!=2){
        #   neighbour_segs_groups_filt <- neighbour_segs_groups[longSegIndx,,drop=F]
        #   neighbour_segs_groups_filt <- neighbour_segs_groups_filt[neighbour_segs_groups_filt$clone != names(clonalres)[[j]],,drop=F]
        #            
        #   mainCN <- round(mean(neighbour_segs_groups_filt$integerCN,na.rm=T))
        # }


        # # mainCN <- neighbour_segs_groups$integerCN[which.max(neighbour_segs_groups$w)]
        # w_max <- neighbour_segs_groups$w[which.max(neighbour_segs_groups$w)]
        # 
        # 
        if(CN_i !=2 &(any(abs(dist_neighbours)<abs(dist_CNratio))|abs(dist_CNratio)>=delta)){
         # CN_mainSeg <- neighbour_segs_groups$integerCN[which.max(neighbour_segs_groups$w)]
        
            seg_dat[i,"integerCN"] <- CN_neighbours[min_index] 
            seg_dat[i,"relativeCN"] <- relativeCN_neighbours[min_index]

        }
        if(CN_i >2 & (ratio_i>=diploidy_lo & ratio_i <= diploidy_hi)){
          # if((dist_CNratio < 0 & dist_neighbours[min_index]>0) |(dist_CNratio > 0 & dist_neighbours[min_index]>0)
          seg_dat[i,"integerCN"] <- 2 
          seg_dat[i,"relativeCN"] <- CNest$ratio[CNest$CN==2]
        }
      }


      CNest <- unique(data.frame(ratio=seg_dat$relativeCN,CN=seg_dat$integerCN))
      CNest <- CNest[order(CNest$CN),]
      CNest <- peakIndex(res$input_BinSegRatio,CNest)
      
      ploidy <- sum(seg_dat$integerCN*seg_dat$w,na.rm=TRUE)/sum(seg_dat$w,na.rm=TRUE)
      res$ploidy <- ploidy
      res$CNest <- CNest
      res$seg.dat <- seg_dat
      
      clonalres[[j]] <- res
    }else{
      clonalres[[j]] <- NULL
    }
    return(clonalres[[j]])
  },clonalres)
  names(clusterout) <- names(clonalres)
  return(clusterout)
}




#' @title ploidyRefine.ref()
#' @description the final estimation : estimate clonal CNV based on reference ploidy
#' @param sampleres initial estimation list including all clusters:
#'                                                        $CNVestimate: initial CNV estimation parameters
#'                                                        $seg.dat: segmentatio data
#'                                                        $input_BinSegRatio: bin data
#' @return for the input cluster output: $input_BinSegRatio
#'                                        $seg.dat
#'                                        $CNest: ratio corresponding to integerCN
#'                                        $ploidy: the overall ploidy
#' @export

ploidyRefine.ref <- function(sampelres,CNest.ref,minCN.frac=0.01,seg_dat_ref=NULL){
  baseCN<- CNest.ref$CN[CNest.ref$base==1]
  baseRatio <- CNest.ref$ratio[CNest.ref$base==1]
  delt.ref <- (CNest.ref$ratio[2]-CNest.ref$ratio[1])/(CNest.ref$CN[2]-CNest.ref$CN[1])
  if(!is.null(seg_dat_ref)){
    overallratio <- sum(seg_dat_ref$relativeCN*seg_dat_ref$w,na.rm = TRUE)/sum(seg_dat_ref$w,na.rm = TRUE)
    
  }
    
  
  clusterout <-  lapply(names(sampelres),function(cluster,sampelres){
    res_ini <- sampelres[[cluster]]
    seg_dat <- res_ini$seg.dat
    seg_dat$segName <- paste(seg_dat$chr,seg_dat$start,seg_dat$end,sep="_")
    colnames(seg_dat) <- gsub("^CNV$","relativeCN",colnames(seg_dat))
    colnames(seg_dat) <- gsub("^integerCNV$","integerCN",colnames(seg_dat))
    CNest <- unique(data.frame(ratio=seg_dat$relativeCN,CN=seg_dat$integerCN))
    CNest <- CNest[order(CNest$CN),]
    res_ini$CNest <- CNest
    res_ini$seg.dat <- seg_dat
    return(res_ini)
    },sampelres)
  # clusterout <- lapply(sort(names(sampelres)), function(cluster,sampelres){
  #   if ("seg.dat" %in% names(sampelres[[cluster]])){
  #     #print(cluster)
  #     seg.dat <- sampelres[[cluster]]$seg.dat
  #     df_seg_C1 <- sampelres[[cluster]]$input_BinSegRatio
  #     integerCNV <- sampelres[[cluster]]
  #     clusterEst <- PloidyCorrect(integerCNV,delt.lim = delt.ref)
  #     if (length(clusterEst) > 1){
  #       likEst <- optimalPloidy(clusterEst)
  #       kindex <- likEst$est[which.min(likEst$AIC)]
  #       output <- clusterEst[[kindex]]
  #       output$integerCN$base=0
  #       output$integerCN$base[output$integerCN$ratio==output$base]=1
  #     }else{
  #       output <- clusterEst[[1]]
  #       output$integerCN$base=0
  #       output$integerCN$base[output$integerCN$ratio==output$base]=1
  #     }
  #     return(output)
  #   }
  # },sampelres)
  names(clusterout) <- names(sampelres)
  # 
  #clusterout2 <- clusterout
  clusterout <- lapply(1:length(clusterout), function(j,clusterout,baseCN,delt.ref){
    if (!is.null(clusterout[[j]])){
      res <- clusterout[[j]]
      seg.dat <- res$seg.dat
      CNestraw <- res$CNest
      CNestraw <- peakIndex(res$input_BinSegRatio,CNestraw)
      if(!(identical(CNest.ref[,1], CNestraw[,1]) & identical(CNest.ref[,2], CNestraw[,2]))){
      #correct dominant ratio
      breakshist <- seq(0,max(res$input_BinSegRatio$SegMean[is.finite(res$input_BinSegRatio$SegMean)]),length.out=100)
      hist_data <- hist(res$input_BinSegRatio$SegMean[is.finite(res$input_BinSegRatio$SegMean)], breaks = breakshist, plot = FALSE)
      max_density_index <- which.max(hist_data$density)
      peak_x <- hist_data$mids[max_density_index]
      # peak_x <- mean(res$input_BinSegRatio$SegMean,na.rm=TRUE)
      bias <- peak_x - CNestraw$ratio[CNestraw$base==1]
      CNestraw$ratio <- CNestraw$ratio + bias
      
      #correct clonal average ratio based on ref.
      ratio_clonalMean <- sum(seg.dat$ratio*seg.dat$w,na.rm = TRUE)/sum(seg.dat$w,na.rm = TRUE)
      if(!is.null(seg_dat_ref)){
        bias2 <- overallratio - ratio_clonalMean
      }else{
        bias2 <- 0
      }
      baseRatio <- baseRatio - bias2
      dis <- abs(CNestraw$ratio-baseRatio)
    
      baseratio <- CNestraw$ratio[which.min(dis)]
      delt <- delt.ref
      integerCN <- round((seg.dat$ratio-baseratio)/delt)+baseCN
      integerCN[integerCN < 0] = 0
      
      CNstate <- unique(integerCN)
      CNstate <- CNstate[order(CNstate)]
      ratio <- baseratio+(CNstate-baseCN)*delt
      CNest <- data.frame(ratio=ratio,CN=CNstate)
      CNest <- peakIndex(res$input_BinSegRatio,CNest)
      df_seg_C1 <- res$input_BinSegRatio
      df_seg_C1$segLen <- df_seg_C1$End - df_seg_C1$Start
      #fnew<- SegFrac(df_seg_C1,CNest,by_term = "SegMean")
      

      relativeCN <- baseratio +(integerCN-baseCN)*delt
      seg.dat$relativeCN <- relativeCN
      seg.dat$integerCN <- integerCN
      
      CNest <- unique(data.frame(ratio=seg.dat$relativeCN,CN=seg.dat$integerCN))
      CNest <- CNest[order(CNest$CN),]
      CNest <- peakIndex(df_seg_C1,CNest)
      
      
      #edited 0514
      CNcount <- as.data.frame(aggregate(seg.dat$w, by=list(ratio=seg.dat$relativeCN,CN=seg.dat$integerCN), FUN=sum))
      CNcount$frac <- CNcount$x/sum(CNcount$x)
      CNcount <- CNcount[CNcount$frac > minCN.frac,]
      m <- min(CNcount$CN)
      if( nrow(CNcount)==1){
        CNest$CN <-2
        CNcount$CN <- 2
      }else if (m <=0 & nrow(CNcount)>1 ){
        CNcount$CN <- CNcount$CN+(1-m)
        CNest$CN <- CNest$CN+(1-m)
      }else{
        if (m >1 & CNcount$CN[which.max(CNcount$frac)]>2){
          CNest$CN <- CNest$CN-(m-1)
        }else if (CNcount$CN[which.max(CNcount$frac)]==1){
          CNest$CN <- CNest$CN+1
          CNcount$CN <- CNcount$CN+1
        }
      }
      baseCN_j_new <- CNcount$CN[which.max(CNcount$frac)]
      baseratio_new <- CNcount$ratio[which.max(CNcount$frac)]
      integerCN <- round((seg.dat$ratio-baseratio_new)/delt)+baseCN_j_new
      integerCN[integerCN < 0] = 0
      relativeCN <- baseratio_new+(integerCN-baseCN_j_new)*delt
      seg.dat$relativeCN <- relativeCN
      seg.dat$integerCN <- integerCN
      
      ploidy <- sum(seg.dat$integerCN*seg.dat$w,na.rm=TRUE)/sum(seg.dat$w,na.rm=TRUE)
      
      clusterout[[j]]$seg.dat <- seg.dat
      clusterout[[j]]$CNest <- CNest
      clusterout[[j]]$ploidy <- ploidy
      }
      
      }else{
        clusterout[[j]] <- NULL
      }
      return(clusterout[[j]])
    },clusterout,baseCN,delt.ref)
  names(clusterout) <- names(sampelres)
  
  # clustersta <- minPloidysummary(clusterout)

  # finalres <- lapply(1:length(clusterout), function(j,clusterout,baseCN){
  #   if (!is.null(clusterout[[j]])){
  #     res <- clusterout[[j]]
  #     input_BinSegRatio <- res$input_BinSegRatio
  #     seg.dat <- res$seg.dat
  #     CNest <- res$newCN
  #     if (baseCN %in% CNest$CN){
  #       baseratio <- CNest$ratio[CNest$CN==baseCN]
  #     }else{
  #       baseratio <- CNest$ratio[1]
  #       baseCN <- CNest$CN[1]
  #     }
  #     if(nrow(CNest)>1){
  #       delt <- (CNest$ratio[2]-CNest$ratio[1])/(CNest$CN[2]-CNest$CN[1])
  #     }else{delt <- delt.ref}
  #    
  #     seg.dat1 <- seg.dat[,-grep("CNV",colnames(seg.dat))]
  #     integerCN <- round((seg.dat$ratio-baseratio)/delt)+baseCN
  #     integerCN[integerCN < 0] = 0
  #     relativeCN <- baseratio+(integerCN-baseCN)*delt
  #     seg.dat1$relativeCN <- relativeCN
  #     seg.dat1$integerCN <- integerCN
  #     ploidy <- sum(seg.dat1$integerCN*seg.dat1$w,na.rm=TRUE)/sum(seg.dat1$w,na.rm=TRUE)
  #     a <- list(input_BinSegRatio = input_BinSegRatio, seg.dat = seg.dat1, CNest = CNest, ploidy = ploidy)
  #   }else{
  #     a <- NULL
  #   }
  #   return(a)
  # },clusterout,baseCN)
  # names(finalres) <- names(clusterout)
  # 
  
  finalres <- clusterout
  listnames <- c()
  res <- list()
  for (i in 1:length(finalres)){
    if (!is.null(finalres[[i]])){
      if (!is.null(finalres[[i]]$input_BinSegRatio)){
        res[[length(res)+1]]=finalres[[i]]
        listnames <- c(listnames,names(finalres)[i])
      }
    }
  }
  names(res) <- listnames
  
  return(res)
  
  
}



#' @title ploidyUpdate()
#' @export

ploidyUpdate <- function(clusterout,baseCN,delt.lim = NULL){
  baseCN_raw <- baseCN
  clonesummary <- function(clusterout){
    clustersta <- do.call(rbind,lapply(1:length(clusterout), function(j,clusterout){
      if (!is.null(clusterout[[j]])){
        res <- clusterout[[j]]
        #input_BinSegRatio <- res$input_BinSegRatio
        seg.dat <- res$seg.dat
        CNest <- res$newCN
        baseratio <- CNest$ratio[1]
        baseCN <- CNest$CN[1]
        delt <- (CNest$ratio[2]-CNest$ratio[1])/(CNest$CN[2]-CNest$CN[1])
        seg.dat1 <- seg.dat[,-grep("CNV",colnames(seg.dat))]
        integerCN <- round((seg.dat$ratio-baseratio)/delt)+baseCN
        seg.dat1$CN <- integerCN
        CNstate <- unique(integerCN)
        CNstate <- CNstate[order(CNstate)]
        ratio <- baseratio+(CNstate-baseCN)*delt
        CNest <- data.frame(ratio=ratio,CN=CNstate)
        CN0index <- 0
        if (0%in%CNest$CN){
          CN0index <- 1
        }
        CNneg <- 0
        if (-1 %in% CNest$CN){
          CNneg <- 1
        }
        c0frac <- sum(seg.dat1$w[seg.dat1$CN==0])/sum(seg.dat1$w)
        c1frac <- sum(seg.dat1$w[seg.dat1$CN<0])/sum(seg.dat1$w)
        ploidy <- sum(seg.dat1$CN*seg.dat1$w,na.rm = TRUE)/sum(seg.dat1$w,na.rm = TRUE)
        return(c(j,baseratio,delt,CN0index,CNneg,c0frac,c1frac,ploidy))
      }
    },clusterout))
    return(clustersta)
  }
  clustersta <- clonesummary(clusterout)
  
  ###correction for CN<=0
  referenceinedx <- which(clustersta[,4]==0&clustersta[,5]==0)
 
  if (length(referenceinedx)==0){
    deltmin <- mean(clustersta[,3],na.rm = TRUE)
  }else{
    deltset <- clustersta[referenceinedx,3]
    deltset <- deltset[!is.na(deltset)]
    if(length(deltset)>0){
      deltmin <- min(deltset)
    }else{
      deltmin <- mean(clustersta[,3],na.rm = TRUE)
    }
  }
  candiateindex <- clustersta[clustersta[,6]>0.01|clustersta[,7] > 0.01,1]
  candiateindex <- candiateindex[!is.na(candiateindex)]
  if (length(candiateindex)>0){
    for (i in candiateindex){
      res <- clusterout[[i]]
      CNest <- res$newCN
      seg.dat <- res$seg.dat
      baseratio <- CNest$ratio[CNest$CN==baseCN]
      delt <- (CNest$ratio[2]-CNest$ratio[1])/(CNest$CN[2]-CNest$CN[1])
      #seg.dat1 <- seg.dat[,-grep("CNV",colnames(seg.dat))]
      integerCN <- round((seg.dat$ratio-baseratio)/delt)+baseCN
      CNstate <- unique(integerCN)
      CNstate <- CNstate[order(CNstate)]
      ratio <- baseratio+(CNstate-baseCN)*delt
      CNest <- data.frame(ratio=ratio,CN=CNstate)
      ########################################################
      if (baseCN > 1){
        deltnew <- (baseratio-min(CNest$ratio))/(baseCN-1)
      }else{
        deltnew <- (baseratio-min(CNest$ratio))
      }
      if (abs(deltnew-deltmin) < deltmin/2){
          Mindex <- floor((max(CNest$ratio)-baseratio)/deltnew)
          mindex <- (min(CNest$ratio)-baseratio)/deltnew
          CNstate <- c((baseCN+mindex):(baseCN+Mindex))
          ratio <- min(CNest$ratio)+(CNstate-1)*deltnew
          CNest <- data.frame(ratio=ratio,CN=CNstate)
          res$newCN <- CNest
          clusterout[[i]] <- res
      }
    }
    clustersta <- clonesummary(clusterout)
    candiateindex <- clustersta[clustersta[,6]>0.01|clustersta[,7]>0.01,1]
    candiateindex <- candiateindex[!is.na(candiateindex)]
    if (length(candiateindex) > 0){
      referenceinedx <- which(clustersta[,4]==0&clustersta[,5]==0)
      if (length(referenceinedx)==0){
        ratioreference <- mean(clustersta[,2],na.rm = TRUE)
      }else{
        ratioreference <- mean(clustersta[referenceinedx,2],na.rm = TRUE)
      }
      for (i in candiateindex){
        res <- clusterout[[i]]
        CNest <- res$newCN
        df_seg_C1 <- res$input_BinSegRatio
        df_seg_C1$segLen <- df_seg_C1$End - df_seg_C1$Start
        CNest$base=0
        delt <- (CNest$ratio[2]-CNest$ratio[1])/(CNest$CN[2]-CNest$CN[1])
        fraw <- SegFrac(df_seg_C1,CNest,by_term = "SegMean")
        Qraw <- round((max(seg.dat$ratio)-min(seg.dat$ratio))/delt)
        baseratio <- ratioreference
        seg.dat <- res$seg.dat
        seg.dat1 <- seg.dat[,-grep("CNV",colnames(seg.dat))]
        integerCN <- round((seg.dat$ratio-baseratio)/delt)+baseCN
        CNstate <- unique(integerCN)
        CNstate <- CNstate[order(CNstate)]
        ratio <- baseratio+(CNstate-baseCN)*delt
        CNest <- data.frame(ratio=ratio,CN=CNstate)
        if (baseCN > 1){
          delt <- (baseratio-min(CNest$ratio))/(baseCN-1)
        }else{
          delt <- (baseratio-min(CNest$ratio))
        }
        #delt <- (baseratio-min(CNest$ratio))/(baseCN-1)
        if (abs(delt-deltmin) < deltmin/2){
          Mindex <- floor((max(CNest$ratio)-baseratio)/delt)
          mindex <- (min(CNest$ratio)-baseratio)/delt
          CNstate <- c((baseCN+mindex):(baseCN+Mindex))
          ratio <- min(CNest$ratio)+(CNstate-1)*delt
          CNest <- data.frame(ratio=ratio,CN=CNstate)
          CNest$base=0
          fnew<- SegFrac(df_seg_C1,CNest,by_term = "SegMean")
          Qnew <- round((max(seg.dat$ratio)-min(seg.dat$ratio))/delt)
          if (fnew/Qnew > fraw/Qraw){
            res$newCN <- CNest
            clusterout[[i]] <- res
          }
        }else{
          kk <- which.min(abs(CNest$ratio-baseratio))
          CNest$CN <- CNest$CN +baseCN - CNest$CN[kk]
          CNest$base=0
          fnew<- SegFrac(df_seg_C1,CNest,by_term = "SegMean")
          Qnew <- round((max(seg.dat$ratio)-min(seg.dat$ratio))/delt)
          if (fnew/Qnew > fraw/Qraw){
            res$newCN <- CNest
            clusterout[[i]] <- res
          }
        }
      }
    }
  }
  

  #clusterout_step1 = clusterout ###
  #clusterout <- clusterout_step1 ###

  ###Correct the minimum CN to 1 for each cluster
  clustersta <- minPloidysummary(clusterout)
  #clustersta <- cbind(clustersta,rep(baseCN,length=dim(clustersta)[1]))
  if (sum(clustersta[,4] > 1)>0){
    #if (sum(clustersta[,4]>1)/dim(clustersta)[1]>0.4){
    #  CNdiff <- (unique(clustersta[clustersta[,4]>1,4])-1)
    #  for (i in 1:length(clusterout)){
    #    res <- clusterout[[i]]
    #    CNest <- res$newCN
    #    CNest$CN <- CNest$CN-CNdiff
    #    clusterout[[i]]$newCN <- CNest
    #  }
    #  clustersta[,5] <- clustersta[,5] - CNdiff
    #}else{
      kk <- clustersta[clustersta[,4]>1,1]
      for (j in kk){
        CNdiff <- clustersta[clustersta[,1]==j,4]-1
        res <- clusterout[[j]]
        #hist_seg1(clusterout[[j]]$input_BinSegRatio,clusterout[[j]],ylim=NULL,color_limit = c(0,3),cluster = j,length.out = 100)
        CNest <- res$newCN
        CNest$CN <- CNest$CN-CNdiff
        clusterout[[j]]$newCN <- CNest
        clustersta[clustersta[,1]==j,5] = clustersta[clustersta[,1]==j,5]-CNdiff
      }
      clustersta[clustersta[,5]<2,5]=2
    #}
  
    ratioreference <- median(clustersta[,2])
    candiateindex <- clustersta[abs(clustersta[,2]-ratioreference) > 0.15,1]
    #clusterout_step2 <- clusterout

    for (i in candiateindex){
      res <- clusterout[[i]]
      CNestraw <- res$newCN
      baseratio <- ratioreference
      delt <- (CNestraw$ratio[2]-CNestraw$ratio[1])/(CNestraw$CN[2]-CNestraw$CN[1])
      seg.dat <- res$seg.dat
      #fill CN if if CN distance larger than 2
      if((CNestraw$CN[2]-CNestraw$CN[1])>1){
        CNestraw <- CNbaseline.fill(CNestraw,delt)
      }
      #seg.dat1 <- seg.dat[,-grep("CNV",colnames(seg.dat))]
      integerCN <- round((seg.dat$ratio-CNestraw$ratio[CNestraw$CN==clustersta[clustersta[,1]==i,5]])/delt)+clustersta[clustersta[,1]==i,5]
      #Q <- max(integerCN)
      CNstate <- unique(integerCN)
      CNstate <- CNstate[order(CNstate)]
      ratio <- baseratio+(CNstate-clustersta[clustersta[,1]==i,5])*delt
      #baseratio+(CNstate-CNestraw$CN[1])*delt
      CNest <- data.frame(ratio=ratio,CN=CNstate)
      delt <- (baseratio-min(CNest$ratio))/(clustersta[clustersta[,1]==i,5]-1)
      Mindex <- floor((max(CNest$ratio)-baseratio)/delt)
      mindex <- (min(CNest$ratio)-baseratio)/delt
      CNstate <- c((clustersta[clustersta[,1]==i,5]+mindex):(clustersta[clustersta[,1]==i,5]+Mindex))
      ratio <- min(CNest$ratio)+(CNstate-1)*delt
      CNestnew <- data.frame(ratio=ratio,CN=CNstate)
      #integerCN <- round((seg.dat$ratio-baseratio)/delt)+baseCN
      #Q <- max(integerCN)
      CNestnew$base <- 0
      CNestraw$base <- 0
      df_seg_C1 <- res$input_BinSegRatio
      df_seg_C1$segLen <- df_seg_C1$End - df_seg_C1$Start
      fnew <- SegFrac(df_seg_C1,CNestnew,by_term = "SegMean")
      fraw <- SegFrac(df_seg_C1,CNestraw,by_term = "SegMean")
      if (fnew > fraw){
        res$newCN <- CNestnew[,1:2]
      }
      clusterout[[i]] <- res
    }
  }
  # 
  clustersta <- minPloidysummary(clusterout)
  baseCN <- round(median(clustersta[,5]))
  clustersta1 <- clonesummary(clusterout)
  candiateindex <- clustersta[clustersta1[,6] > 0.01|clustersta1[,7]>0.01,1]
  candiateindex <- candiateindex[!is.na(candiateindex)]
  if (length(candiateindex) > 0){
    referenceinedx <- which(clustersta1[,4]==0&clustersta1[,5]==0)
    if (length(referenceinedx)==0){
      ratioreference <- mean(clustersta1[,2],na.rm=TRUE)
      deltreference <- mean(clustersta1[,3],na.rm=TRUE)
    }else{
      ratioreference <- mean(clustersta1[referenceinedx,2],na.rm=TRUE)
      deltreference <- mean(clustersta1[referenceinedx,3],na.rm=TRUE)
    }
    for (i in candiateindex){
      res <- clusterout[[i]]
      seg.dat <- res$seg.dat
      CNestraw <- res$newCN
      CNestraw$base <- 0
      df_seg_C1 <- res$input_BinSegRatio
      df_seg_C1$segLen <- df_seg_C1$End - df_seg_C1$Start
      fraw <- SegFrac(df_seg_C1,CNestraw,by_term = "SegMean")
      baseratio <- ratioreference
      if(nrow(CNestraw)>1){
        delt <- (CNestraw$ratio[2]-CNestraw$ratio[1])/(CNestraw$CN[2]-CNestraw$CN[1])
      }else{
        delt <- deltreference
      }
      
      seg.dat <- res$seg.dat
      #seg.dat1 <- seg.dat[,-grep("CNV",colnames(seg.dat))]
      integerCN <- round((seg.dat$ratio-baseratio)/delt)+baseCN
      CNstate <- unique(integerCN)
      CNstate <- CNstate[order(CNstate)]
      ratio <- baseratio+(CNstate-baseCN)*delt
      CNest <- data.frame(ratio=ratio,CN=CNstate)
      CNest <- CNest[CNest$CN>=0,]
      CNest$base <- 0
      if (abs(delt - deltreference) < 0.1){
        fnew <- SegFrac(df_seg_C1,CNest,by_term = "SegMean")
        #fnew <- fnew/((max(seg.dat$ratio)-min(seg.dat$ratio))/delt)
        if (fraw > fnew){
          CNest <- CNestraw
        }
      }else{
        #Mindex <- floor((max(CNest$ratio)-CNest$ratio[CNest$CN==clustersta[i,5]])/deltreference)
        deltcount <- floor((max(CNest$ratio)-min(CNest$ratio))/deltreference)
        if(length(deltcount)<=1){
          deltcount <- deltcount+1
        }
        CNstate <- c(1:deltcount)
        ratio <- min(CNest$ratio)+(CNstate-min(CNest$ratio))*deltreference
        CNest <- data.frame(ratio=ratio,CN=CNstate)
        integerCN <- round((seg.dat$ratio-CNest$ratio[1])/deltreference)+CNest$CN[1]
        CNstate <- unique(integerCN)
        CNstate <- CNstate[order(CNstate)]
        CNest <- CNest[CNest$CN%in%CNstate,]
        CNest$base <- 0
        fnew <- SegFrac(df_seg_C1,CNest,by_term = "SegMean")
        #fnew <- fnew*deltreference
        if (fraw > fnew){
          CNest <- CNestraw
        }
      }
      if(nrow(CNest)>1){
        delt <- (CNest$ratio[2]-CNest$ratio[1])/(CNest$CN[2]-CNest$CN[1])
        integerCN <- round((seg.dat$ratio-CNest$ratio[CNest$CN==baseCN])/delt)+baseCN
        if (min(integerCN) > 1){
          CNest$CN <- CNest$CN+(1-min(CNest$CN))
        }
      }else{
        CNest$CN <- 2
      }

      clusterout[[i]]$newCN <- CNest[,1:2]
      
    }
  }
  
  ####clusterout_step4 = clusterout
  
  {
    #Round1:similar clonal ploidy
    clustersta <- minPloidysummary(clusterout)
    clustersta2 <- clonesummary(clusterout)
    
    delt_mean <- median(clustersta[,3],na.rm=T)
    ploidy <- round(median(clustersta2[,8],na.rm=T))
    ploidy_mean <- mean(ploidy)
    candiateindex_row <- which(clustersta2[,3]<delt_mean | round(clustersta2[,3]/delt_mean)!=1|round(abs(clustersta2[,8]-ploidy_mean))>=1)
    candiateindex <-clustersta2[candiateindex_row,1]
    referenceinedx <- which(!clustersta2[,1] %in%candiateindex)
    if(length(referenceinedx)>0){
      delt_ref<- median(clustersta2[referenceinedx,3],na.rm=T)
      ratioreference <- median(clustersta2[referenceinedx,2],na.rm=T)
      baseCN_ref <- round(median(clustersta[referenceinedx,5],na.rm=T))
    }else{
      delt_ref<- delt_mean
      ratioreference <- mean(clustersta2[,2],na.rm=T)
      baseCN_ref <- round(mean(clustersta[,5],na.rm=T))
    }
    
    # candiateindex <- which(round(clustersta2[,3]/delt_mean)>=2)
    for(i in candiateindex){
      res <- clusterout[[i]]
      CNestraw <- res$newCN
      delt <- delt_ref
      seg.dat <- res$seg.dat
      if(nrow(CNestraw)>1 & (CNestraw$CN[2]-CNestraw$CN[1])>1){
        CNestraw <- CNbaseline.fill(CNestraw,delt)
      }
      baseratio <- ratioreference
      integerCN <- round((seg.dat$ratio-baseratio)/delt)+baseCN_ref
      CNstate <- unique(integerCN)
      CNstate <- CNstate[order(CNstate)]
      ratio <- baseratio+(CNstate-baseCN_ref)*delt
      CNest <- data.frame(ratio=ratio,CN=CNstate)
      CNest$base <- 0
      CNest <- CNest[CNest$CN>=0,]
      
      clusterout[[i]]$newCN <- CNest
    }
    
   # clusterout_R1 <- clusterout
    ###Round2: correct clone-delt based on CNbase_raw
    baseRatio <- unlist(lapply(clusterout,function(x)x$base))
    clustersta <- minPloidysummary(clusterout)
    clustersta2 <- clonesummary(clusterout)
    candiateindex <- clustersta2[which(abs(clustersta[,5]-round(2*baseRatio))>1),1] 
    referenceinedx <- which(!clustersta2[,1] %in%candiateindex)
    if(length(referenceinedx)>0){
      delt_ref<- median(clustersta2[referenceinedx,3],na.rm=T)
      ratioreference <- median(clustersta2[clustersta2[,1]%in%referenceinedx,2],na.rm=T)
      baseCN_ref <- round(median(clustersta[clustersta2[,1]%in%referenceinedx,5],na.rm=T))
    }else{
      delt_ref<- delt_mean
      ratioreference <- median(clustersta2[,2],na.rm=T)
      baseCN_ref <- round(median(clustersta[,5],na.rm=T))
    }
    
  for(i in candiateindex){
    res <- clusterout[[i]]
    seg.dat <- res$seg.dat
    CNestraw <- res$newCN
    CNestraw$base <- 0
    df_seg_C1 <- res$input_BinSegRatio
    df_seg_C1$segLen <- df_seg_C1$End - df_seg_C1$Start
    fraw <- SegFrac(df_seg_C1,CNestraw,by_term = "SegMean")
    
    mainRatio <- clustersta[clustersta2[,1]%in%i,6]
    mainCN_ref <- round(2*mainRatio)
    if(mainRatio>CNestraw$ratio[1] &mainCN_ref>CNestraw$CN[1]){
      delt <- (mainRatio-CNestraw$ratio[1])/(mainCN_ref-CNestraw$CN[1])

    }else{
      delt <- delt_ref
      mainCN_ref <- baseCN_ref
      mainRatio <- ratioreference
     
    }
    integerCN <- round((seg.dat$ratio-mainRatio)/delt)+mainCN_ref
    CNstate <- unique(integerCN)
    CNstate <- CNstate[order(CNstate)]
    ratio <- mainRatio+(CNstate-mainCN_ref)*delt
    
    CNest <- data.frame(ratio=ratio,CN=CNstate)
    CNest <- CNest[CNest$CN>=0,]
    CNest$base <- 0
    fnew <- SegFrac(df_seg_C1,CNest,by_term = "SegMean")
    if (fraw > fnew){
      CNest <- CNestraw
    }
    if(nrow(CNest)>1){
      delt <- (CNest$ratio[2]-CNest$ratio[1])/(CNest$CN[2]-CNest$CN[1])
      integerCN <- round((seg.dat$ratio-CNest$ratio[1])/delt)+CNest$CN[1]
      if (min(integerCN) > 1){
        CNest$CN <- CNest$CN+(1-min(CNest$CN))
      }
    }else{
      CNest$CN <- 2
    }
    
    
    clusterout[[i]]$newCN <- CNest[,1:2]
  }
    

  } 
  
  #clusterout_step3 <- clusterout
  clustersta <- minPloidysummary(clusterout)
  delt_mean <- median(clustersta[,3],na.rm=T)
  baseCN <- c()
  for (j in 1:length(clusterout)){
    if (!is.null(clusterout[[j]])){
      res <- clusterout[[j]]
      seg.dat <- res$seg.dat
      CNest <- res$newCN
      if(nrow(CNest)>1){
        delt <- (CNest$ratio[2]-CNest$ratio[1])/(CNest$CN[2]-CNest$CN[1])
      }else{
        delt <- delt_mean
      }
      
      integerCN <- round((seg.dat$ratio-CNest$ratio[1])/delt)+CNest$CN[1]
      seg.dat$CN <- integerCN
      CNcount <- as.data.frame(aggregate(seg.dat$w, by=list(CN=seg.dat$CN,CNV=seg.dat$CNV), FUN=sum))
      CNcount$frac <- CNcount$x/sum(CNcount$x)
      CNcount <- CNcount[CNcount$frac > 0.01,]
      m <- min(CNcount$CN)
      if( nrow(CNcount)==1){
        CNest$CN <-2
        res$newCN <- CNest
        clusterout[[j]] <- res
        CNcount$CN <- 2
      }else if (m <=0 & nrow(CNcount)>1 ){
        CNcount$CN <- CNcount$CN+(1-m)
        CNest$CN <- CNest$CN+(1-m)
        res$newCN <- CNest
        clusterout[[j]] <- res
      }else{
        if (m >1 & CNcount$CN[which.max(CNcount$frac)]>2){
          CNest$CN <- CNest$CN-(m-1)
          if(nrow(CNest)>1){
            delt <- (CNest$ratio[2]-CNest$ratio[1])/(CNest$CN[2]-CNest$CN[1])
          }else{
            delt <- delt_mean
          }
          
          integerCN <- round((seg.dat$ratio-CNest$ratio[1])/delt)+CNest$CN[1]
          seg.dat$CN <- integerCN
          CNcount1 <- as.data.frame(aggregate(seg.dat$w, by=list(CN=seg.dat$CN), FUN=sum))
          CNcount1$frac <- CNcount1$x/sum(CNcount1$x)
          m <- min(CNcount1$CN)
          if( nrow(CNcount1)==1){
            CNest$CN <-2
            res$newCN <- CNest
            clusterout[[j]] <- res
            CNcount$CN <- 2
          }else if (m ==1 & nrow(CNcount1)>1){
            res$newCN <- CNest
            clusterout[[j]] <- res
            CNcount$CN <- CNcount$CN-1
          }else {
            if (0 %in% CNcount1$CN){
              if (CNcount1$frac[CNcount1$CN==0] < 0.005){
                res$newCN <- CNest
                clusterout[[j]] <- res
                CNcount$CN <- CNcount$CN-1
              }
            }
          }
        }else if (CNcount$CN[which.max(CNcount$frac)]==1){
          CNest$CN <- CNest$CN+1
          res$newCN <- CNest
          clusterout[[j]] <- res
          CNcount$CN <- CNcount$CN+1
        }else{
          if(nrow(CNest)>1){
            delt <- (CNest$ratio[2]-CNest$ratio[1])/(CNest$CN[2]-CNest$CN[1])
          }else{
            delt <- delt_mean
          }
          
          integerCN <- round((seg.dat$ratio-CNest$ratio[1])/delt)+CNest$CN[1]
          seg.dat$CN <- integerCN
          CNcount1 <- as.data.frame(aggregate(seg.dat$w, by=list(CN=seg.dat$CN), FUN=sum))
          CNcount1$frac <- CNcount1$x/sum(CNcount1$x)
          m <- min(CNcount1$CN)
          if (0 %in% CNcount1$CN){
            if (CNcount1$frac[CNcount1$CN==0] > 0.01){
              CNest$CN <- CNest$CN+1
              res$newCN <- CNest
              clusterout[[j]] <- res
              CNcount$CN <- CNcount$CN+1
            }
          }
        }
      }
      
      baseCN <- c(baseCN,sum(CNcount$CN*CNcount$x)/sum(CNcount$x))
      
      
      #if (1%in% CNcount$CN){
      #  if (CNcount$frac[CNcount$CN==1]>0.30){
      #    CNest$CN <- CNest$CN+1
      #    res$newCN <- CNest
      #    clusterout[[j]] <- res
      #    CNcount$CN <- CNcount$CN+1
      #  }else if (CNcount$frac[CNcount$CN==1]<0.01 & CNcount$CN[which.max(CNcount$frac)]>2){
      #    CNest$CN <- CNest$CN-1
      #    res$newCN <- CNest
      #    clusterout[[j]] <- res
      #    CNcount$CN <- CNcount$CN-1
      #  }else{
      #    CNcount <- CNcount[CNcount$frac > 0.01,]
      #    if (min(CNcount$CN) < 1){
      #      CNest$CN <- CNest$CN+(1-min(CNcount$CN))
      #      res$newCN <- CNest
      #      clusterout[[j]] <- res
      #      CNcount$CN <- CNcount$CN+(1-min(CNcount$CN))
      #    }
      #  }
      #}
      #baseCN <- c(baseCN,sum(CNcount$CN*CNcount$x)/sum(CNcount$x)) 
    }
  }
  baseCN <- round(min(baseCN))
  return(list(clusterout=clusterout,baseCN=baseCN))
}


#' @title hist_seg2()
#' @description plot the barplot of ratio as well as corresponding integer CN
#' @param clusterEst:the colonal estimation
#' @export
hist_seg2 <- function(clusterEst,ylim=NULL,color_limit = c(0.5,4),cluster=NULL,length.out=100,
                      color_hist_gradient=FALSE,hist_color=NULL,gtitle=NULL){
  suppressMessages({
    library(BSgenome)
    library(dplyr)})
  #colorss <- c("black","blue","grey50","green", "yellow", "orange","brown")

  df_seg_C1 <- clusterEst$input_BinSegRatio
  CNseg.dat <- clusterEst$seg.dat
  ploidy <- sum(CNseg.dat$integerCN*CNseg.dat$w)/sum(CNseg.dat$w)
  #CNseg.dat$segName <- paste(CNseg.dat$chr,CNseg.dat$start,CNseg.dat$end,sep="_")
  colnames(CNseg.dat) <- gsub("^CNV$","relativeCN", colnames(CNseg.dat))
  colnames(CNseg.dat) <- gsub("^integerCNV","integerCN", colnames(CNseg.dat))
  
  CNseg.names <- colnames(CNseg.dat)
  a <- length(grep("integerCN",CNseg.names))
  if (a==1){
  CNseg.dat_slim <- CNseg.dat[,c("segName","relativeCN","integerCN"),drop=F]
  }else{
  CNseg.dat_slim <- CNseg.dat[,c("segName","relativeCN","integerCN",paste0("integerCN",1:(a-1))),drop=F]
  }
  
  g <- getBSgenome("hg38", masked=FALSE)
  chrLine <- data.frame(chrom=1:24, length=seqlengths(g)[1:24])
  chromnum=sapply(strsplit(rownames(df_seg_C1),'-|_|:'),'[',1)
  if(grepl('chr',chromnum[1])){
    chromnum=sapply(strsplit(chromnum,'hr'),'[',2) ## if chrommnum with chr
  }else{
    chromnum=chromnum
  }
  chromnum <- gsub("X","23",chromnum)
  chromnum <- gsub("Y","24",chromnum)
  df_seg_C1$chrom <- chromnum
  chrLine <- chrLine[chrLine$chrom %in% as.character(unique(chromnum)),,drop=F] 
  genomeLength <- sum(chrLine)
  if(nrow(chrLine)>1){
    chrPosStart = c(0,chrLine[1:(nrow(chrLine)-1),2])
  }else(
    chrPosStart = 0
  )
  chrLine$chrPosStart_abs <- cumsum(chrPosStart)
  df_seg_C1 <- merge(df_seg_C1,chrLine,by="chrom",all.x=T)
  df_seg_C1 <- df_seg_C1[order(as.integer(df_seg_C1$chrom),df_seg_C1$chrPosStart_abs),]
  df_seg_C1$segStrat <- as.numeric(sapply(strsplit(df_seg_C1$segName,'_'),'[',2))+as.numeric(df_seg_C1$chrPosStart_abs)
  df_seg_C1$segEnd <- as.numeric(sapply(strsplit(df_seg_C1$segName,'_'),'[',3))+as.numeric(df_seg_C1$chrPosStart_abs)
  
  df_seg_C1 <- merge(df_seg_C1,CNseg.dat_slim,by="segName",all.x=T)
  colnames(df_seg_C1)[grepl("relativeCN",colnames(df_seg_C1))] <- "ratio_map"
  
  if(is.null(ylim)){ylim=c(0,max(df_seg_C1$SegMean[is.finite(df_seg_C1$SegMean)])+0.5)}
  breakshist <- seq(0,max(df_seg_C1$SegMean[is.finite(df_seg_C1$SegMean)]),length.out=length.out)
  
  if (a == 1){
  segSet_C1 = unique(df_seg_C1[,c("binID","segID","SegMean","segName","ratio_map","integerCN"),drop=F]) ###
  }else{
   segSet_C1 = unique(df_seg_C1[,c("binID","segID","SegMean","segName","ratio_map","integerCN",paste0("integerCN",1:(a-1))),drop=F]) ###
  }
  
  
  true_loc = TRUE
  if(true_loc){
    df_seg_C1$segLen <- as.numeric(df_seg_C1$segEnd)-as.numeric(df_seg_C1$segStrat)
  }else{
    df_seg_C1$segLen <- as.numeric(df_seg_C1$length_seg)
  }   
  
  segLen_sum <- sum(unique(df_seg_C1[,c("segID","segLen")]$segLen))
  #df_seg_C1$wei <- df_seg_C1$segLen/segLen_sum
  
  df_seg_C1$wei <- df_seg_C1$length_bin/segLen_sum
  
  integerCNVcolumn <- grep("integerCN",colnames(df_seg_C1))
  labs <- unique(df_seg_C1[integerCNVcolumn])
  labs <- labs[!is.na(labs[,1]),]
  if (a >1){
   labs <- apply(labs, 1, function(x){
     paste(x,collapse=" - ")
   })
  }
  if(is.null(cluster)){cluster=""}

  if(color_hist_gradient){
    if(is.null(hist_color)){
      colorss <- c("#023e8a","grey50", "#cb793a", "#9a031e","#6a040f","#370617")
      #colorss <- c("black","#023e8a","grey50","#fcdc4d", "#cb793a", "#9a031e","#5f0f40")
      
    }else{
      colorss <- hist_color
    }
  }else{
    if(is.null(hist_color)){
      colorss <- "#472d30"
    }else{
      colorss <- hist_color[1]
    }
    
  }
  
  hist_plt <- ggplot(df_seg_C1,aes(x=SegMean,fill = ..x..))+
    geom_vline(xintercept = unique(df_seg_C1$ratio_map),col="grey",linetype = "dotted")
  
  if(color_hist_gradient){
    hist_plt <- hist_plt+
      geom_histogram(aes(y=..density.., weight = wei),breaks=breakshist)+
      scale_fill_gradientn(name = "",colours = colorss,#c("#003049", "#669bbc","#7e7f83","#c1121f","#780000"),
                         limits  = color_limit,oob = scales::squish)
  }else{
    hist_plt <- hist_plt+
      geom_histogram(aes(y=..density.., weight = wei),color=colorss,fill=colorss,breaks=breakshist)
  }
  if(is.null(gtitle)){
    gtitle <- paste0("Ploidy of Cluster ",cluster,": ",round(ploidy,2))
  }
  hist_plt <- hist_plt+
    ggtitle(gtitle)+
    coord_cartesian(xlim=ylim)+
    theme_bw() +
    theme_classic()+
    coord_flip(xlim=ylim)+
    scale_x_continuous(
      sec.axis = sec_axis( trans = ~.,breaks=unique(df_seg_C1$ratio_map)[!is.na(unique(df_seg_C1$ratio_map))], labels =labs, name="")
    )+
    theme(plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
          plot.title = element_text(size = 15), #legend.position = "bottom",
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.title=element_text(size=15),
          axis.text = element_text(size = 12))+
    labs(x = "", y = "",fill="")
  hist_plt
  return(hist_plt)
}



#' @title PloidyCorrect()
#' @description the relative ratio baseline and distance of integerCN states correction
#' @param integerCNV: the initial estimation from ABSOLUTE model
#' @return corrected estimation 
#' @export

PloidyCorrect <- function(integerCNV,delt.lim =0.3){
  df_seg_C1 <- integerCNV$input_BinSegRatio
  CNseg.dat <- integerCNV$seg.dat
  CNseg.dat$segName <- paste(CNseg.dat$chr,CNseg.dat$start,CNseg.dat$end,sep="_")
  colnames(CNseg.dat) <- gsub("^relativeCN$","CNV",colnames(CNseg.dat))
  colnames(CNseg.dat) <- gsub("^integerCN$","integerCNV",colnames(CNseg.dat))
  
  CNseg.dat_slim <- CNseg.dat[,c("segName","CNV","integerCNV"),drop=F]
  df_seg_C2 <- merge(df_seg_C1,CNseg.dat_slim,by="segName",all.x=T)
  colnames(df_seg_C2)[grepl("^CNV$",colnames(df_seg_C2))] <- "ratio_map"
  # true_loc = TRUE
  # if(true_loc){
  df_seg_C2$segLen <- as.numeric(df_seg_C2$End)-as.numeric(df_seg_C2$Start)
  # }else{
  #   df_seg_C2$segLen <- as.numeric(df_seg_C2$length_seg)
  # }   
  segLen_sum <- sum(df_seg_C2$segLen)
  df_seg_C2$wei <- df_seg_C2$segLen/segLen_sum
  df_seg_C2 <- df_seg_C2[!is.na(df_seg_C2$ratio_map),]
  initialCN <- as.data.frame(unique(cbind(df_seg_C2$ratio_map,df_seg_C2$integerCNV)))
  colnames(initialCN) <- c("ratio","CN")
  initialCN <- initialCN[order(initialCN$CN),]
  initialCN$base=0
  
  if (dim(initialCN)[1]!=1){
    initialCN1 <- CNupdate(df_seg_C2,initialCN)
    disupdate <- DeltSeek(df_seg_C2,initialCN1)
    disupdate=disupdate[disupdate$delt>delt.lim,,drop=FALSE]
    
    CNest <- FinalCN(initialCN1,disupdate,df_seg_C2)
  }else{
    initialCN1 <- initialCN
    initialCN1$base = 1
    CNest <- list(initialCN[,1:2])
  }
  clusterEst <- lapply(1:length(CNest), function(i,CNest,CNseg.dat){
    delt <- (CNest[[i]]$ratio[2]-CNest[[i]]$ratio[1])/(CNest[[i]][2,2]-CNest[[i]][1,2])
    CNseg.dat$integerCNV <- CNest[[i]][1,2]+round((CNseg.dat$ratio-CNest[[i]]$ratio[1])/delt)
    CN0=CNest[[i]]$ratio[1]-(CNest[[i]][1,2])*delt
    CNseg.dat$CNV <- CNseg.dat$integerCNV*delt+CN0
    CNseg.names <- colnames(CNseg.dat)
    if (dim(CNest[[i]])[2]>2){
      CN <- do.call(cbind,lapply(3:dim(CNest[[i]])[2], function(j,i,CNest,delt,CNseg.dat){
        CN0=CNest[[i]][1,j]
        CN <- CN0+round((CNseg.dat$ratio-CNest[[i]]$ratio[1])/delt)
        CN[CN<0] <- 0
        return(CN)
      },i,CNest,delt,CNseg.dat))
      CNseg.dat <- cbind(CNseg.dat,CN)
      colnames(CNseg.dat) <- c(CNseg.names,paste0("integerCNV",c(1:(dim(CNest[[i]])[2]-2))))
    }
    a <- list(input_BinSegRatio = integerCNV$input_BinSegRatio,seg.dat = CNseg.dat,integerCN = CNest[[i]],base=initialCN1$ratio[initialCN1$base==1])
    return(a)
  },CNest,CNseg.dat)
  
  return(clusterEst)
}



#' @title optimalPloidy()
#' @description ploidy estimation selecting 
#' @param clusterEst: the corrected CN estimation
#' @return the likelihood for each ploidy estimation
#' @export

optimalPloidy <- function(clusterEst){
  likout <- do.call(rbind,lapply(1:length(clusterEst), function(j,clusterEst){
    seg.dat <- clusterEst[[j]]$seg.dat
    seg.dat$segCount <- round((as.numeric(seg.dat$end)-as.numeric(seg.dat$start))/500000)
    CN <- seg.dat %>% group_by(CNV) %>% summarise(sum_count=sum(segCount),sd_CNV=sd(ratio),
                                                  .groups = 'drop')
    CN <- as.data.frame(CN)
    CN$sd_CNV[is.na(CN$sd_CNV)] <- 0
    CN$fre <- CN$sum_count/sum(CN$sum_count)
    ratio <- rep(seg.dat$ratio,times=seg.dat$segCount)
    CNratio <- as.numeric(CN[,1])
    pro <- as.numeric(CN[,4])
    sd0 <- as.numeric(CN$sd_CNV)
    a <- do.call(cbind,lapply(1:length(pro), function(i,pro,CNratio,ratio){
      a <- pro[i]*dnorm(ratio,mean=CNratio[i],sd=sd0[i])
      return(a)
    },pro,CNratio,ratio))
    loglik <- sum(log(apply(a,1,sum)))
    RSS <- sum((seg.dat$ratio-seg.dat$CNV)^2*seg.dat$segCount)/sum(seg.dat$segCount)
    paraN <- length(CNratio)
    AIC <- 2*paraN-2*loglik
    deltAIC <- 2*paraN+sum(seg.dat$segCount)*log(RSS)
    delt <- clusterEst[[j]]$integerCN$ratio[2]-clusterEst[[j]]$integerCN$ratio[1]
    output <- data.frame(loglik = loglik, paraN = paraN,RSS=RSS,AIC=AIC,deltAIC=deltAIC,delt=delt)
    return(output)
  },clusterEst))
  likout$est <- c(1:length(clusterEst))
  likout <- likout[order(likout$AIC),]
  dAIC <- likout$AIC-min(likout$AIC)
  likout$weight <- exp(-dAIC)/sum(exp(-dAIC))
  return(likout)
}



#' @title CNupdate()
#' @description estimating the corrected ratio baseline and distance used in PloidyCorrect
#' @param df_seg_C2: segmentation data
#' @param initialCN: the initial estimation from ABSOLUTE model 
#' @return the likelihood for each ploidy estimation
#' @export

CNupdate <- function(df_seg_C2,initialCN){
  initialCN1 <- peakIndex(df_seg_C2,initialCN)
  ratioindex <- df_seg_C2[df_seg_C2$integerCNV==initialCN1$CN[initialCN1$base==1],]
  ratioindex <- as.data.frame(table(ratioindex$SegMean))
  ratioindex[,3] <- sapply(ratioindex[,1], function(value,df_seg_C2){
    subdata <- df_seg_C2[df_seg_C2$SegMean[is.finite(df_seg_C2$SegMean)]==value,]
    a <- sum(subdata$segLen*subdata$wei)/sum(df_seg_C2$segLen)
    return(a)
  },df_seg_C2)
  ratioindex[,1] <- as.numeric(as.character(ratioindex[,1]))
  ratioindex$score <- ratioindex[,2]*ratioindex[,3]
  ratioindex <- ratioindex[order(ratioindex[,2],decreasing=T),]
  ratioindex <- ratioindex[1,1]
  diff1 <- abs(initialCN$ratio[initialCN$CN==initialCN1$CN[initialCN1$base==1]]-ratioindex)
  diff2 <- abs(initialCN1$ratio[initialCN1$CN==initialCN1$CN[initialCN1$base==1]]-ratioindex)
  if (min(diff1,diff2) > 0.06){
    initialCN2 <- initialCN1
    d <- ratioindex-initialCN1$ratio[initialCN1$base==1]
    initialCN2$ratio <- initialCN1$ratio+d
    f1 <- SegFrac(df_seg_C2,initialCN1) 
    f2 <- SegFrac(df_seg_C2,initialCN2) 
    if (f2 > f1){
      initialCN1 <- initialCN2
    }
  }else{
    if (diff1 <= diff2){
      initialCN$base <- initialCN1$base
      f1 <- SegFrac(df_seg_C2,initialCN) 
      f2 <- SegFrac(df_seg_C2,initialCN1) 
      if (f2 < f1){
        initialCN1 <- initialCN
      }
    }else{
      initialCN2 <- initialCN1
      d <- ratioindex-initialCN1$ratio[initialCN1$base==1]
      initialCN2$ratio <- initialCN1$ratio+d
      f2 <- SegFrac(df_seg_C2,initialCN2)
      f1 <- SegFrac(df_seg_C2,initialCN1)
      if (f2 > f1){
        initialCN1 <- initialCN2
      }
    }
  }
  return(initialCN1)
}



#' @title SegFrac()
#' @description estimating explained genome fraction based on the integerCN estimation used in PloidyCorrect
#' @param df_seg_C1: segmentation data
#' @param initialCN: the initial estimation from ABSOLUTE model 
#' @return the fraction of genome 
#' @export

SegFrac <- function(df_seg_C1,initialCN,by_term = "SegMean"){
  ratioFre <- as.data.frame(table(df_seg_C1[,by_term]))
  ratioFre$Var1 <- as.numeric(as.character(ratioFre$Var1))
  if(nrow(initialCN)>1){
    if (max(df_seg_C1$SegMean) > max(initialCN$ratio)){
      delt=abs(initialCN$ratio[1]-initialCN$ratio[2])
      k <- ceiling((max(df_seg_C1$SegMean)-max(initialCN$ratio))/delt)
      initialCN <- rbind(initialCN,data.frame(ratio=max(initialCN$ratio)+c(1:k)*delt,CN=max(initialCN$CN)+c(1:k),base=rep(0,length=k)))
    }
    delt <- initialCN$ratio[2]-initialCN$ratio[1]
  }else{delt=0.5}

  if (delt/2 < 0.1){
    cutoff <- delt/2-0.05
  }else{
    cutoff <- 0.1
  }
  a <- sum(sapply(initialCN$ratio, function(j,df_seg_C1,cutoff){
    dis <- abs(j-as.numeric(df_seg_C1$SegMean))
    if (min(dis) < cutoff){
      a <- sum(df_seg_C1$segLen[dis < cutoff])/sum(df_seg_C1$segLen)
    }else{
      a <- 0
    }
    return(a)
  },df_seg_C1,cutoff))
  return(a)
}



#' @title DeltSeek()
#' @description finding the possible distance between integerCN states used in PloidyCorrect
#' @param df_seg_C1: segmentation data
#' @param initialCN1: the estimation of ratio baseline 
#' @return all distance value
#' @export

DeltSeek <- function(df_seg_C1,initialCN1,by_term = "SegMean"){
  genoFrac0 <- SegFrac(df_seg_C1,initialCN1)
  base0 <- initialCN1$ratio[initialCN1$base==1]
  delt0 <- data.frame(delt=(initialCN1$ratio[2]-initialCN1$ratio[1]),genoFrac = genoFrac0)
  pp <- density(df_seg_C1[,by_term],width=0.05)
  dens <- pp$y
  ratioX <- pp$x
  breapRatio <- c()
  for (ii in 2:(length(dens)-1)){
    if (dens[ii-1] <= dens[ii] & dens[ii+1] <= dens[ii]){
      breapRatio <- rbind(breapRatio,c(ratioX[ii],dens[ii]))
    }
  }
  breapRatio <- breapRatio[breapRatio[,2]>0.1,,drop=F]
  breapRatio1 <- breapRatio[order(breapRatio[,2],decreasing=T),]
  j=2
  deltobs <- NULL
  while(j <= dim(breapRatio1)[1]){
    deltobs <- abs(breapRatio1[1,1]-breapRatio1[j,1])
    if (deltobs >= 0.25){
      break
    }else{
      j=j+1
    }
  }
  if (is.null(deltobs)){
    deltobs <- max(abs(breapRatio1[1,1]-breapRatio1[,1]))
  }
  
  if (dim(breapRatio)[1]==0){
    genoFrac <- SegFrac(df_seg_C1,initialCN1)
    delt <- data.frame(delt=abs(initialCN1$ratio[1]-initialCN1$ratio[2]),genoFrac=genoFrac)
  }else{
    k <- which.min(abs(breapRatio[,1]-base0))
    dis <- abs(breapRatio[1,1]-breapRatio[k,1])
    delt <- abs(breapRatio[1,1]-initialCN1$ratio)
    if (min(delt) < 0.3){
      kk = which(delt < 0.3)
    }else{
      kk = which.min(delt)
    }
    breaks1 <- abs(initialCN1$CN[initialCN1$base==1]-initialCN1$CN[kk])
    breaks1 <- breaks1[breaks1!=0]
    if(any(breaks1 !=0)){
      if (dis == 0){
        delt=deltobs
      }else{
        delt <- dis/breaks1
      }
    }
    baseCN = initialCN1$CN[initialCN1$base==1]
    if (length(delt)==1){
      if (abs(delt - deltobs) < 0.01){
        delt <- data.frame(delt=delt,genoFrac=genoFrac0)
      }else{
        deltnew <- c()
        if (delt != deltobs){
          if (delt < deltobs){
            signindex <- 1
          }else{
            signindex <- -1
          }
          i=0
          while(abs(deltobs-(delt+signindex*0.01*i))>0.01){
            d <- delt+signindex*0.01*i
            initialCN1$ratio <- base0 + (initialCN1$CN-baseCN)*d
            genoFrac <- SegFrac(df_seg_C1,initialCN1)
            deltnew <- rbind(deltnew,data.frame(delt=d,genoFrac=genoFrac))
            i=i+1
          }
        }
        delt <- data.frame(delt=deltnew[which.max(deltnew[,2]),1],genoFrac=deltnew[which.max(deltnew[,2]),2]) 
      }
    }else{
      delt <- do.call(rbind,lapply(delt[delt>0], function(d,df_seg_C1,initialCN1){
        deltnew <- c()
        if (d != deltobs){
          if (d < deltobs){
            signindex <- 1
          }else{
            signindex <- -1
          }
          i=0
          while(abs(deltobs-(d+signindex*0.01*i))>0.01){
            dd <- d+signindex*0.01*i
            initialCN1$ratio <- base0 + (initialCN1$CN-baseCN)*dd
            genoFrac <- SegFrac(df_seg_C1,initialCN1)
            deltnew <- rbind(deltnew,data.frame(delt=dd,genoFrac=genoFrac))
            i=i+1
          }
        }
        return(c(deltnew[which.max(deltnew[,2]),1],deltnew[which.max(deltnew[,2]),2]))
      },df_seg_C1,initialCN1))
    }
  }
  
  
  delt <- data.frame(delt =delt[,1],genoFrac=delt[,2])
  delt <- rbind(delt,delt0)
  delt <- unique(delt)
  
  return(delt)
}





#' @title peakIndex()
#' @description finding peak of ration distribution
#' @param df_seg_C1: segmentation data
#' @param initialCN: the estimation of ratio baseline 
#' @return possible ratio vale corresponding to integerCN
#' @export

peakIndex <- function(df_seg_C1,initialCN,index_col="SegMean"){
  colnames(initialCN)[grepl("ratio",colnames(initialCN),ignore.case = T)] <- "ratio"
  colnames(initialCN)[grepl("CN",colnames(initialCN),ignore.case = T)] <- "CN"
  initialCN <- initialCN[order(initialCN$ratio),,drop=F]
  initialCN$base = 0
  pp <- density(df_seg_C1[,index_col],width=0.05)
  dens <- pp$y
  ratioX <- pp$x
  ratio.pos <- unlist(lapply(2:(length(ratioX)), function(i,ratioX,dens){
    if (dens[i-1]<dens[i] & dens[i+1]<dens[i]){
      return(i)
    }
  },ratioX,dens))
  ratioPeak <- data.frame(ratio=ratioX[ratio.pos],dens=dens[ratio.pos])
  ratioPeak <- ratioPeak[order(ratioPeak$dens,decreasing=T),]
  delt <- abs(ratioPeak$dens[2:dim(ratioPeak)[1]]-ratioPeak$dens[1:(dim(ratioPeak)[1]-1)])
  deltd <- abs(delt[2:length(delt)]-delt[1:(length(delt)-1)])
  peakpos <- which.max(deltd)
  
  if(length(deltd)>2){
    for (i in 1:(length(deltd)-2)){
      if ( deltd[i] < 0.1 & deltd[i+1] < 0.1){
        peakpos <- i-1
        break
      }
    }
  }


  ratioPeak <- ratioPeak[1:peakpos,]
  ratiomax <-ratioPeak$ratio[which.max(ratioPeak$dens)]
  ratioPeak <- ratioPeak[ratioPeak$ratio >= ratiomax-0.25 & ratioPeak$ratio <= ratiomax+0.25,]
  
  a <- mean(ratioPeak$ratio)
  dis <- na.omit(abs(a-initialCN$ratio))
  if (min(dis) < 0.02){
    initialCN$base[which.min(dis)]=1
    return(initialCN)
  }else{
    base0 <- a
    baseCN <- initialCN$CN[which.min(dis)]
    initialCN1 <- initialCN
    if(nrow(initialCN)>1){
      delt= abs(initialCN1$ratio[1]-initialCN1$ratio[2])/(initialCN1$CN[2]-initialCN1$CN[1])
      initialCN1$ratio <- base0 + (initialCN1$CN-baseCN)*delt
      initialCN1$base[which(initialCN1$ratio==base0)]=1
    }else{
      initialCN1$base=1
    }
    return(initialCN1) 
  }
}


#' @title FinalCN()
#' @description determine the final integerCN
#' @param initialCN1: the estimation of ratio baseline 
#' @param disupdate: all possible distance
#' @param df_seg_C1: segmentation data
#' @return CNest: the final integerCN baseline
#' @export
#' 
FinalCN <- function(initialCN1,disupdate,df_seg_C1){
  cnm <- colnames(disupdate)
  delt <- (initialCN1$ratio[2]-initialCN1$ratio[1])/(initialCN1$CN[2]-initialCN1$CN[1])
  genoFrac0 <- SegFrac(df_seg_C1,initialCN1)
  disupdate <- rbind(disupdate,c(delt,genoFrac0))
  colnames(disupdate) <- cnm
  
  disupdate <- unique(disupdate)
  disupdate <- disupdate[order(disupdate[,2],decreasing=T),]
  if (length(disupdate$genoFrac[disupdate$genoFrac>0.2])==0){
    disupdate <- disupdate[which.max(disupdate$genoFrac),]
  }else{
    disupdate <- disupdate[disupdate[,2]>0.2,]
  }
  
  CNres <- lapply(1:dim(disupdate)[1],function(x,initialCN1,df_seg_C1){
    delt <- as.numeric(disupdate[x,1])
    base0 <- initialCN1$ratio[initialCN1$base==1]
    baseCN <- initialCN1$CN[initialCN1$base==1]
    initialCN1$ratio <- base0 + (initialCN1$CN-baseCN)*delt
    k <- ceiling((max(df_seg_C1$SegMean)-max(initialCN1$ratio))/delt)
    if (k > 0){
      initialCN1 <- rbind(initialCN1,data.frame(ratio=max(initialCN1$ratio)+c(1:k)*delt,CN=max(initialCN1$CN)+c(1:k),base=rep(0,length=k)))
    }
    initialCN2 <-ploidyEst(initialCN1,df_seg_C1)
    minratio <- min(df_seg_C1$SegMean)
    k <- which.min(abs(initialCN2$ratio-minratio))
    initialCN2 <- initialCN2[k:dim(initialCN2)[1],]
    CNindex <- do.call(cbind,lapply(2:dim(initialCN2)[2], function(j,initialCN2,df_seg_C1,delt){
      CN0 <- initialCN2[1,1]-(initialCN2[1,j]*delt)
      CN <- round((df_seg_C1$SegMean-CN0)/delt)
      if (min(CN)<0){
        if (min(initialCN2[,j])<0){
          kk <- abs(min(initialCN2[,j]))
          initialCN2[,j] <- initialCN2[,j]+kk
        }
      }
      return(initialCN2[,j])
    },initialCN2,df_seg_C1,delt))
    columIndex <- NULL
    if (dim(CNindex)[2]>1){
      for (j in 2:dim(CNindex)[2]){
        if (sum(round(CNindex[,j])!=round(CNindex[,j-1]))==0){
          columIndex <- c(columIndex,j)
        }
      }
    }
    if (!is.null(columIndex)){
      initialCN2 <- cbind(initialCN2[,1],CNindex[,-columIndex])
    }
    initialCN2 <- as.data.frame(initialCN2)
    colnames(initialCN2) <- c("ratio",paste0("CN",c(1:(dim(initialCN2)[2]-1))))
    
    
    maxratio <- as.numeric(quantile(df_seg_C1$SegMean,c(0.99)))
    k <- which.min(abs(initialCN2$ratio-maxratio))
    initialCN2 <- initialCN2[1:k,]
    maxratio <- as.numeric(quantile(df_seg_C1$SegMean,c(0.05)))
    k <- which.min(abs(initialCN2$ratio-maxratio))
    initialCN2 <- initialCN2[k:dim(initialCN2)[1],]
    
    return(list(initialCN=initialCN2))
  },initialCN1,df_seg_C1)
  
  disupdate1 <- do.call(rbind,apply(disupdate, 1, function(x,disupdate){
    dis <- abs(disupdate$delt-as.numeric(x[1]))
    subres <- disupdate[dis<0.02,]
    if (sum(dis<0.02)>1){
      subres <- subres[subres[,1]>0.1,]
      subres <- subres[which.max(subres[,2]),]
    }
    return(subres)
  },disupdate))
  disupdate1 <- unique(disupdate1)
  
  disIndex <- match(disupdate1$delt,disupdate$delt)
  
  CNres <- lapply(disIndex, function(i,CNres){
    return(CNres[[i]]$initialCN)
  },CNres)
  return(CNres)
}



#' @title ploidyEst()
#' @description the overall ploidy of tumor clone
#' @param CNest: the final integerCN baseline
#' @param df_seg_C1: segmentation data
#' @return whole genome ploidy
#' @export
#' 

ploidyEst <- function(CNest,df_seg_C1){
  #CNest <- lapply(1:length(CNest), function(i,df_seg_C1){
  CNest1 <- CNest
  ploidy <- round(CNest1$ratio[CNest1$base==1]*2)
  if (ploidy != CNest1$CN[CNest1$base==1]){
    d= ploidy-CNest1$CN[CNest1$base==1]
    CNest1 <- cbind(CNest1,CNest1$CN+d)
  }
  base0 <- CNest1$ratio[CNest1$base==1]-(CNest1$CN[CNest1$base==1]-1)*abs(CNest1$ratio[1]-CNest1$ratio[2])
  base1 <- base0
  mratio <- df_seg_C1$SegMean[df_seg_C1$SegMean<base1]
  dis <- abs(base1-mratio)
  prob <- length(dis[dis>0.2])/dim(df_seg_C1)[1]
  while(1){
    if (prob < 0.01){
      break
    }else{
      base1 <-  base1-abs(CNest1$ratio[1]-CNest1$ratio[2])
      mratio <- df_seg_C1$SegMean[df_seg_C1$SegMean<base1]
      dis <- abs(base1-mratio)
      prob <- length(dis[dis>0.2])/dim(df_seg_C1)[1]
    }
  }
  if (base1 != base0){
    a <- round((CNest1$ratio-base1)/abs(CNest1$ratio[1]-CNest1$ratio[2])+1)
    if (sum(a!=CNest1[,dim(CNest1)[2]])!=0){
      CNest1 <- cbind(CNest1,a)
    }
  }
  CNest1 <- CNest1[,-3]
  colnames(CNest1) <- c("ratio",paste("CN",c(1:(dim(CNest1)[2]-1)),sep=""))
  CNest1.names <- colnames(CNest1)
  a <- CNest1$ratio*2-1
  aindex <- which(a<0.5&a>0)
  if (length(aindex)>1){
    CN1index <- aindex[which.min(a[aindex])]
  }else{
    CN1index <- aindex
  }
  dis <- CNest1$CN1[CN1index]-1
  CNadd <- NULL
  CNadd <-unlist(lapply(2:dim(CNest1)[2], function(j,CNest1,dis){
    b <- sum(CNest1[,j]-(CNest1$CN1-dis))
    if (b == 0){
      CNadd <- 1
      return(CNadd)
    }
  },CNest1,dis))
  if (is.null(CNadd)){
    CNest1 <- cbind(CNest1,CNest1$CN1-dis)
    colnames(CNest1) <- c(CNest1.names,paste0("CN",length(CNest1.names):((dim(CNest1)[2])-1)))
  }
  return(CNest1)
  #},df_seg_C1)
  return(CNest1)
}







#' @title hist_seg()
#' @export

#df_seg_C1 is the ratio(segScore$seg_score_binLevel[[i]]); 
#integCNV_Ci is the CNV estimation result from doCNV1_v2()

hist_seg <- function(df_seg_C1,integerCNV,ylim=NULL,color_limit = c(0,3),cluster,length.out=80){
  suppressMessages({
    library(BSgenome)
    library(dplyr)})
  colorss <- c("black","blue","grey50","green", "yellow", "orange","brown")
  CNseg.dat <- integerCNV$seg.dat
  CNseg.dat$segName <- paste(CNseg.dat$chr,CNseg.dat$start,CNseg.dat$end,sep="_")
  CNseg.names <- colnames(CNseg.dat)
  a <- length(grep("integerCNV",CNseg.names))
  if (a==1){
    CNseg.dat_slim <- CNseg.dat[,c("segName","CNV","integerCNV"),drop=F]
  }else{
    CNseg.dat_slim <- CNseg.dat[,c("segName","CNV","integerCNV",paste0("integerCNV",1:(a-1))),drop=F]
  }
  
  g <- getBSgenome("hg38", masked=FALSE)
  chrLine <- data.frame(chrom=1:24, length=seqlengths(g)[1:24])
  chromnum=sapply(strsplit(rownames(df_seg_C1),'-|_|:'),'[',1)
  if(grepl('chr',chromnum[1])){
    chromnum=sapply(strsplit(chromnum,'hr'),'[',2) ## if chrommnum with chr
  }else{
    chromnum=chromnum
  }
  chromnum <- gsub("X","23",chromnum)
  chromnum <- gsub("Y","24",chromnum)
  df_seg_C1$chrom <- chromnum
  chrLine <- chrLine[chrLine$chrom %in% as.character(unique(chromnum)),,drop=F] 
  genomeLength <- sum(chrLine)
  if(nrow(chrLine)>1){
    chrPosStart = c(0,chrLine[1:(nrow(chrLine)-1),2])
  }else(
    chrPosStart = 0
  )
  chrLine$chrPosStart_abs <- cumsum(chrPosStart)
  df_seg_C1 <- merge(df_seg_C1,chrLine,by="chrom",all.x=T)
  df_seg_C1 <- df_seg_C1[order(as.integer(df_seg_C1$chrom),df_seg_C1$chrPosStart_abs),]
  df_seg_C1$segStrat <- as.numeric(sapply(strsplit(df_seg_C1$segName,'_'),'[',2))+as.numeric(df_seg_C1$chrPosStart_abs)
  df_seg_C1$segEnd <- as.numeric(sapply(strsplit(df_seg_C1$segName,'_'),'[',3))+as.numeric(df_seg_C1$chrPosStart_abs)
  
  df_seg_C1 <- merge(df_seg_C1,CNseg.dat_slim,by="segName",all.x=T)
  colnames(df_seg_C1)[grepl("^CNV$",colnames(df_seg_C1))] <- "ratio_map"
  
  #qq <- quantile(df_seg_C1$SegMean[is.finite(df_seg_C1$SegMean)],prob=c(0.01,0.99))
  
  #ylim=c(0,as.numeric(qq[2]))
  
  if(is.null(ylim)){ylim=c(0,max(df_seg_C1$SegMean[is.finite(df_seg_C1$SegMean)])+0.5)}
  breakshist <- seq(min(ylim),max(ylim),length.out=80)
  
  if (a == 1){
    segSet_C1 = unique(df_seg_C1[,c("binID","segID","SegMean","segName","ratio_map","integerCNV"),drop=F]) ###
  }else{
    segSet_C1 = unique(df_seg_C1[,c("binID","segID","SegMean","segName","ratio_map","integerCNV",paste0("integerCNV",1:(a-1))),drop=F]) ###
  }
  
  
  true_loc = TRUE
  if(true_loc){
    df_seg_C1$segLen <- as.numeric(df_seg_C1$segEnd)-as.numeric(df_seg_C1$segStrat)
  }else{
    df_seg_C1$segLen <- as.numeric(df_seg_C1$length_seg)
  }   
  
  segLen_sum <- sum(unique(df_seg_C1[,c("segID","segLen")]$segLen))
  df_seg_C1$wei <- df_seg_C1$segLen/segLen_sum
  
  integerCNVcolumn <- grep("integerCNV",colnames(df_seg_C1))
  labs <- unique(df_seg_C1[integerCNVcolumn])
  labs <- labs[!is.na(labs[,1]),]
  if (a >1){
    labs <- apply(labs, 1, function(x){
      paste(x,collapse=" - ")
    })
  }
  hist_plt <- ggplot(df_seg_C1,aes(x=SegMean,fill = ..x..))+
    geom_vline(xintercept = unique(df_seg_C1$ratio_map),col="grey",linetype = "dotted")+ 
    geom_histogram(aes(y=..density.., weight = wei),breaks=breakshist)+
    ggtitle(paste0("Summary histogram of cluster ",cluster))+
    coord_cartesian(xlim=ylim)+
    theme_bw() +
    theme_classic()+
    coord_flip(xlim=ylim)+
    #scale_fill_manual(values = colorss) +
    scale_fill_gradientn(name = "",colours = colorss,#c("#003049", "#669bbc","#7e7f83","#c1121f","#780000"),
                         limits  = color_limit,oob = scales::squish)+
    # scale_fill_gradient2(name = "", low = "#004e89",mid = "#e9d8a6", high = "#c1121f",midpoint = 1,
    #                      breaks = c(0, 1, ceiling(max(a[,p_x]))),labels= c(0,1,ceiling(max(a[,p_x]))))+
    scale_x_continuous(
      sec.axis = sec_axis( trans = ~.,breaks=unique(df_seg_C1$ratio_map)[!is.na(unique(df_seg_C1$ratio_map))], labels =labs, name="")
    )+
    theme(plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
          plot.title = element_text(size = 15), #legend.position = "bottom",
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.title=element_text(size=15),
          axis.text = element_text(size = 12))+
    labs(x = "", y = "",fill="")
  hist_plt
  return(hist_plt)
}




#' @title hist_seg1()
#' @export

hist_seg1 <- function(df_seg_C1,integerCNV,ylim=NULL,color_limit = c(0,3),cluster,length.out=100){
  suppressMessages({
    library(BSgenome)
    library(dplyr)})
  colorss <- c("black","blue","grey50","green", "yellow", "orange","brown")
  CNseg.dat <- integerCNV$seg.dat
  CNseg.dat$segName <- paste(CNseg.dat$chr,CNseg.dat$start,CNseg.dat$end,sep="_")
  CNseg.names <- colnames(CNseg.dat)
  a <- length(grep("integerCNV",CNseg.names))
  if (a==1){
    CNseg.dat_slim <- CNseg.dat[,c("segName","CNV","integerCNV"),drop=F]
  }else{
    CNseg.dat_slim <- CNseg.dat[,c("segName","CNV","integerCNV",paste0("integerCNV",1:(a-1))),drop=F]
  }
  
  g <- getBSgenome("hg38", masked=FALSE)
  chrLine <- data.frame(chrom=1:24, length=seqlengths(g)[1:24])
  chromnum=sapply(strsplit(rownames(df_seg_C1),'-|_|:'),'[',1)
  if(grepl('chr',chromnum[1])){
    chromnum=sapply(strsplit(chromnum,'hr'),'[',2) ## if chrommnum with chr
  }else{
    chromnum=chromnum
  }
  chromnum <- gsub("X","23",chromnum)
  chromnum <- gsub("Y","24",chromnum)
  df_seg_C1$chrom <- chromnum
  chrLine <- chrLine[chrLine$chrom %in% as.character(unique(chromnum)),,drop=F] 
  genomeLength <- sum(chrLine)
  if(nrow(chrLine)>1){
    chrPosStart = c(0,chrLine[1:(nrow(chrLine)-1),2])
  }else(
    chrPosStart = 0
  )
  chrLine$chrPosStart_abs <- cumsum(chrPosStart)
  df_seg_C1 <- merge(df_seg_C1,chrLine,by="chrom",all.x=T)
  df_seg_C1 <- df_seg_C1[order(as.integer(df_seg_C1$chrom),df_seg_C1$chrPosStart_abs),]
  df_seg_C1$segStrat <- as.numeric(sapply(strsplit(df_seg_C1$segName,'_'),'[',2))+as.numeric(df_seg_C1$chrPosStart_abs)
  df_seg_C1$segEnd <- as.numeric(sapply(strsplit(df_seg_C1$segName,'_'),'[',3))+as.numeric(df_seg_C1$chrPosStart_abs)
  
  df_seg_C1 <- merge(df_seg_C1,CNseg.dat_slim,by="segName",all.x=T)
  colnames(df_seg_C1)[grepl("^CNV$",colnames(df_seg_C1))] <- "ratio_map"
  
  if(is.null(ylim)){ylim=c(0,max(df_seg_C1$SegMean[is.finite(df_seg_C1$SegMean)])+0.5)}
  breakshist <- seq(0,max(df_seg_C1$SegMean[is.finite(df_seg_C1$SegMean)]),length.out=length.out)
  
  if (a == 1){
    segSet_C1 = unique(df_seg_C1[,c("binID","segID","SegMean","segName","ratio_map","integerCNV"),drop=F]) ###
  }else{
    segSet_C1 = unique(df_seg_C1[,c("binID","segID","SegMean","segName","ratio_map","integerCNV",paste0("integerCNV",1:(a-1))),drop=F]) ###
  }
  
  
  true_loc = TRUE
  if(true_loc){
    df_seg_C1$segLen <- as.numeric(df_seg_C1$segEnd)-as.numeric(df_seg_C1$segStrat)
  }else{
    df_seg_C1$segLen <- as.numeric(df_seg_C1$length_seg)
  }   
  
  segLen_sum <- sum(unique(df_seg_C1[,c("segID","segLen")]$segLen))
  #df_seg_C1$wei <- df_seg_C1$segLen/segLen_sum
  
  df_seg_C1$wei <- df_seg_C1$length_bin/segLen_sum
  
  integerCNVcolumn <- grep("integerCNV",colnames(df_seg_C1))
  labs <- unique(df_seg_C1[integerCNVcolumn])
  labs <- labs[!is.na(labs[,1]),]
  if (a >1){
    labs <- apply(labs, 1, function(x){
      paste(x,collapse=" - ")
    })
  }
  
  hist_plt <- ggplot(df_seg_C1,aes(x=SegMean,fill = ..x..))+
    geom_vline(xintercept = unique(df_seg_C1$ratio_map),col="grey",linetype = "dotted")+ 
    geom_histogram(aes(y=..density.., weight = wei),breaks=breakshist)+
    ggtitle(paste0("Summary histogram of cluster ",cluster))+
    coord_cartesian(xlim=ylim)+
    theme_bw() +
    theme_classic()+
    coord_flip(xlim=ylim)+
    #scale_fill_manual(values = colorss) +
    scale_fill_gradientn(name = "",colours = colorss,#c("#003049", "#669bbc","#7e7f83","#c1121f","#780000"),
                         limits  = color_limit,oob = scales::squish)+
    # scale_fill_gradient2(name = "", low = "#004e89",mid = "#e9d8a6", high = "#c1121f",midpoint = 1,
    #                      breaks = c(0, 1, ceiling(max(a[,p_x]))),labels= c(0,1,ceiling(max(a[,p_x]))))+
    scale_x_continuous(
      sec.axis = sec_axis( trans = ~.,breaks=unique(df_seg_C1$ratio_map)[!is.na(unique(df_seg_C1$ratio_map))], labels =labs, name="")
    )+
    theme(plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
          plot.title = element_text(size = 15), #legend.position = "bottom",
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.title=element_text(size=15),
          axis.text = element_text(size = 12))+
    labs(x = "", y = "",fill="")
  hist_plt
  return(hist_plt)
}



