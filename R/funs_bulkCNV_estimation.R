
# This is a program to estimate clonal integer CNVs from ratio.
suppressMessages({
  library(dplyr)
  library(numDeriv)
  library(ggpubr)
  library(BSgenome)
  library(futile.logger)
  library(ggplot2)
  library(gridExtra)
  })

# source("dataProcess.R",chdir = T)
# source("results.check.R",chdir = T)
# source("figurePlot.R",chdir = T)
# source("findmode.R",chdir = T)
# source("optEstimation.R",chdir = T)
# source("funs_utils.R",chdir = T)
# source("integerCNV.correction.R",chdir = T)


#' @title CNV_esti_ini()
#' @description estimate initial clonal-level CNVs
#' @export
CNV_esti_ini <- function(outdir=".",
                         segScore_file="segScore_median.rds",
                         filt_seg = TRUE,
                         #quantile_seg_cutoff=0.04,
                         size_seg_cutoff = 5,
                         length_seg_cutoff =1e6,
                         segValue_method="median",
                         genome="hg38",
                         true_loc = TRUE,
                         outplot_name = NULL,
                         outputFigure = TRUE){
  if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
  if (Reduce("|", is(segScore_file) == "character")) {
    flog.info(sprintf("Loading segment ratio file: %s", segScore_file))
    if(substr(segScore_file, nchar(segScore_file)-3, nchar(segScore_file)) == ".rds") {
      segScore <- readRDS(segScore_file)
    }
  }else if(Reduce("|", is(segScore_file) == "list")){
    segScore <- segScore_file
  }else{
    stop("Can not find segScore_file as a list or filename of '.rds'.")
  }
  CNVres_total <- list()
  if(any(grepl("seg_score_binLevel",names(segScore)))){
    ls_binseg <- segScore$seg_score_binLevel
  }else{
    stop("Can not find seg_score_binLevel in the list of segScore_file.")
  }
  p_ls <- list()
  for(i in 1:length(ls_binseg)){
    Cname <- names(ls_binseg)[i]
    CNVdata <- ls_binseg[[i]]
    #nbs <- sort(table(CNVdata$segID))
    #hist(CNVdata$length_seg/1000,breaks=100,xlim = c(0,100))
    if(filt_seg){
      ###filter segments
      # length_seg_cutoff <- quantile(CNVdata$length_seg,quantile_seg_cutoff)
      
      CNVdata <- CNVdata %>%
        dplyr::group_by(segID) %>%
        dplyr::do(data.frame(., w = length(.$segID)))%>%
        dplyr::filter(w>=size_seg_cutoff,length_seg>=length_seg_cutoff)%>%
        as.data.frame()
    }
    
    flog.info(paste0("Estimating ",Cname,"..."))
    flog.info(paste0("Median length_seg is ",median(CNVdata$length_seg/1e3)," kb"))
    #CNVdata <- seg_filter(CNVdata,min_bins_per_seg=1)
    CNVres=dataProcess(CNVdata,method=segValue_method)
    
    seg.dat=CNVres$seg
    seg.dat=data.frame(chr=seg.dat$chr,start=seg.dat$start,end=seg.dat$end,ratio=as.numeric(seg.dat$ratio),sd=as.numeric(seg.dat$SD),w=as.numeric(seg.dat$bin),abspos=as.numeric(seg.dat$abspos),absend=as.numeric(seg.dat$absend))
    seg.dat$sd[is.na(seg.dat$sd)]=0
    
    seg.dat_new <- seg.dat
    output<-tryCatch(doCNV1_v2(seg.dat_new),error=function(e) NA )
    if(any(is.na(output))){
      message(paste0("Faild to Estimate for cluster ",Cname,"!"))
      CNVres_total[[Cname]] <- output 
      #return(NA)
    }else{
      seg_df <- output$seg.dat
      colnames(seg_df)[4] <- "SegMean"
      if(is.null(outplot_name)){
        outplot_name_i= paste0("bulkCNV_",Cname)
      }else{
        outplot_name_i= paste0(outplot_name,"_",Cname)
      }
      if(outputFigure){
        plt <- cellLineSeg_v3(CNVres=CNVres,integerCNV=output,outdir=outdir,fname=outplot_name_i, ggarrange_width=c(10,5),
                              height = 3,width=10,outPlot=FALSE)
        p_ls[[Cname]]<- plt$ggarranged_p
      }
      output[["input_BinSegRatio"]] <- ls_binseg[[i]]
      CNVres_total[[Cname]] <- output 
      
    }
    rm(CNVres,output,seg.dat_new)
  }
  if(outputFigure){
    pcom <- cowplot::plot_grid(plotlist=p_ls,ncol = 1,align="v")
    ggsave(paste0(outdir, "/",outplot_name,".png"),pcom, width=14, height=3*length(ls_binseg),device = png,bg="white")
    #saveRDS(CNVres_total,paste0(outdir,"/CNV_res_initial.rds"))
  }
  return(CNVres_total)
}

  ###------------------------------------###
  #merge segment ratio and integer CNV results
  ###------------------------------------###

#' @title addCNV2BinDF()
#' @export

addCNV2BinDF <- function(output.cnv,ratio.df,genome="hg38",true_loc=TRUE){
  ratio.df$row.name <- rownames(ratio.df)
  g <- getBSgenome(genome, masked=FALSE)
  CNseg.dat <- output.cnv$seg.dat
  CNseg.dat$segName <- paste(CNseg.dat$chr,CNseg.dat$start,CNseg.dat$end,sep="_")
  CNseg.dat_slim <- CNseg.dat[,c("segName","CNV","integerCNV"),drop=F]

  chrLine <- data.frame(chrom=1:24, length=seqlengths(g)[1:24])
  chromnum=sapply(strsplit(rownames(ratio.df),'-|_|:'),'[',1)
  if(grepl('chr',chromnum[1])){
      chromnum=sapply(strsplit(chromnum,'hr'),'[',2) ## if chrommnum with chr
  }else{
    chromnum=chromnum
  }
  chromnum <- gsub("X","23",chromnum)
  chromnum <- gsub("Y","24",chromnum)
  ratio.df$chrom <- chromnum
  chrLine <- chrLine[chrLine$chrom %in% as.character(unique(chromnum)),,drop=F] 
  genomeLength <- sum(chrLine)
  if(nrow(chrLine)>1){
    chrPosStart = c(0,chrLine[1:(nrow(chrLine)-1),2])
  }else(
    chrPosStart = 0
  )
  chrLine$chrPosStart_abs <- cumsum(chrPosStart)
  ratio.df <- merge(ratio.df,chrLine,by="chrom",all.x=T)
  ratio.df <- ratio.df[order(as.integer(ratio.df$chrom),ratio.df$chrPosStart_abs),]
  ratio.df$segStrat <- as.numeric(sapply(strsplit(ratio.df$segName,'_'),'[',2))+as.numeric(ratio.df$chrPosStart_abs)
  ratio.df$segEnd <- as.numeric(sapply(strsplit(ratio.df$segName,'_'),'[',3))+as.numeric(ratio.df$chrPosStart_abs)
  
  ratio.df <- merge(ratio.df,CNseg.dat_slim,by="segName",all.x=T)
  colnames(ratio.df)[grepl("^CNV$",colnames(ratio.df))] <- "ratio_map"
  if(true_loc){
    ratio.df$segLen <- as.numeric(ratio.df$segEnd)-as.numeric(ratio.df$segStrat)
  }else{
    ratio.df$segLen <- as.numeric(ratio.df$length_seg)
  }  
  segLen_sum <- sum(unique(ratio.df[,c("segID","segLen")]$segLen))
  ratio.df$wei <- ratio.df$segLen/segLen_sum
  rownames(ratio.df) <- ratio.df$row.name
  ratio.df <- ratio.df[,!grepl("row.name",colnames(ratio.df)),drop=F]
  ratio.df <- ratio.df[order(as.numeric(ratio.df$chrom),as.numeric(ratio.df$Start)),,drop=F]
  
  return(ratio.df)

}


#' @title CNV_esti()
#' @param segScore_file the file with segment ratio of each bin for each cell
#' @param length_seg_cutoff if filt_seg=TRUE, set the quantile threshold length of segments
#' @param genome the reference genome for the genomic location annotation
#' @param offset.adjust Adjust the ratio to nearest integer.
#' @export
CNV_esti <- function(outdir=".",
                     segScore_file="segScore_median.rds",
                     filt_seg = TRUE,
                     offset.adjust = TRUE,
                     #quantile_seg_cutoff=0.04,
                     length_seg_cutoff =1e6,
                     genome="hg38",
                     true_loc = TRUE,
                     outplot_name = NULL,
                     delta.min = 0.4,
                     delta.max = 0.7,
                     outputFigure = TRUE
                     
){
  if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
  # outplot_path <- paste0(outdir,"/Figures");if(!file.exists(outplot_path)){dir.create(outplot_path,recursive=T)}
  
  
  if (Reduce("|", is(segScore_file) == "character")) {
    flog.info(sprintf("Loading segment ratio file: %s", segScore_file))
    if(substr(segScore_file, nchar(segScore_file)-3, nchar(segScore_file)) == ".rds") {
      segScore <- readRDS(segScore_file)
    }
  }else if(Reduce("|", is(segScore_file) == "list")){
    segScore <- segScore_file
  }else{
    stop("Can not find segScore_file as a list or filename of '.rds'.")
  }
  
  CNVres_total <- list()
  
  if(any(grepl("seg_score_binLevel",names(segScore)))){
    ls_binseg <- segScore$seg_score_binLevel
  }else{
    stop("Can not find seg_score_binLevel in the list of segScore_file.")
  }
  
  #mv_dist <- 0
  p_ls <- list()
  for(i in 1:length(ls_binseg)){
    Cname <- names(ls_binseg)[i]
    CNVdata <- ls_binseg[[i]]
    #nbs <- sort(table(CNVdata$segID))
    #hist(CNVdata$length_seg/1000,breaks=100,xlim = c(0,100))
    if(filt_seg){
      ###filter segments
      # length_seg_cutoff <- quantile(CNVdata$length_seg,quantile_seg_cutoff)
      
      CNVdata <- CNVdata %>%
        dplyr::group_by(segID) %>%
        do(data.frame(., w = length(.$segID)))%>%
        dplyr::filter(length_seg>=length_seg_cutoff)%>%
        as.data.frame()
    }
    
    flog.info(paste0("Estimating ",Cname,"..."))
    flog.info(paste0("Median length_seg is ",median(CNVdata$length_seg/1e3)," kb"))
    #CNVdata <- seg_filter(CNVdata,min_bins_per_seg=1)
    CNVres=dataProcess(CNVdata)
    seg.dat=CNVres$seg
    seg.dat=data.frame(chr=seg.dat$chr,start=seg.dat$start,end=seg.dat$end,ratio=as.numeric(seg.dat$ratio),sd=as.numeric(seg.dat$SD),w=as.numeric(seg.dat$bin),abspos=as.numeric(seg.dat$abspos),absend=as.numeric(seg.dat$absend))
    seg.dat$sd[is.na(seg.dat$sd)]=0
    
    if(offset.adjust){
      #Correct the overall deviation of the ratio according to the distance between the maximum peak and the nearest integer copy number
      mv_dist <- global_offset_correct(CNVdata$SegMean)
    }else{
      mv_dist <- 0
    }
    
    #update the ratio
    seg.dat_new <- seg.dat
    if(mv_dist!=0){
      seg.dat_new$ratio <- as.numeric(seg.dat_new$ratio) - (mv_dist)
      seg.dat_new <- seg.dat_new[seg.dat_new$ratio>=0,,drop=F]
      #seg.dat$ratio[seg.dat$ratio<0] <- 0
      
    }
    
    # seg_plot(CNVdata,name.data=paste0(Cname,".step0-Ratio"),plotDir=outdir,value.bin="binRatio",value.segment="SegMean")
    # binDFnew <- CNVdata
    # binDFnew$binRatio <- as.numeric(binDFnew$binRatio) - (mv_dist)
    # binDFnew$SegMean <- as.numeric(binDFnew$SegMean) - (mv_dist)
    # seg_plot(binDFnew,name.data=paste0(Cname,".step1-RatioCorrect"),plotDir=outdir,value.bin="binRatio",value.segment="SegMean")
        
    Q.max <- 2*round(max(seg.dat_new$ratio))
    if(Q.max>8){
      Q.max <- 8
    }else if(Q.max==0){
      Q.max <- 2
    }
    output<-tryCatch(doCNV1_v2(seg.dat_new,Q=Q.max),error=function(e) NA )
    if(any(is.na(output))){
      message(paste0("Faild to Estimate for cluster ",Cname,"!"))
      CNVres_total[[Cname]] <- output 
      #return(NA)
      }else{
    
    # 
    delta_initial <- output$CNVestimate$delta
    #mse_initial <- mse_loss(seg_df,y_true.term="SegMean",y_pred.term="CNV",weight.term="length")
    if(delta_initial>=delta.min & delta_initial<=delta.max){
      delta.min2 <- max(round(delta_initial,2)-0.05,delta.min)
      delta.max2 <- min(round(delta_initial,2) + 0.1,delta.max)
    }else{
      delta.min2 <- delta.min
      delta.max2 <- delta.max
    }

    Dsets <- seq(delta.min2,delta.max2,by=0.002)

    
    seg_df <- output$seg.dat
    seg_df$length_seg <- as.numeric(seg_df$end)-as.numeric(seg_df$start)
    colnames(seg_df)[4] <- "SegMean"
    mse_res <- c()
    for (d in Dsets ){
      integerCNV <- round(output$seg.dat$ratio/d)
      ratio_map<- d*integerCNV
      seg_df$ratio_map <- ratio_map
      
      mse_d <- mse_loss(seg_df,y_true.term="SegMean",y_pred.term="ratio_map",weight.term="length_seg")
      mse_res <- c(mse_res,mse_d)
    }
    delta_new <- Dsets[which.min(mse_res)]
    integerCNV <- round(output$seg.dat$ratio/delta_new)
    output$seg.dat$integerCNV <- integerCNV
    
    if(mv_dist!=0){
      output$seg.dat$CNV <- delta_new*integerCNV + mv_dist
    }else{ output$seg.dat$CNV <- delta_new*integerCNV}
    
    output$CNVestimate$delta <- delta_new
    # }
    
    ## update ploidy choice
    seg_df <- output$seg.dat
    colnames(seg_df)[4] <- "SegMean"
    initialCN <- unique(seg_df[,c("CNV","integerCNV")])
    colnames(initialCN) <- c("ratio","CN")
    initialCN <- peakIndex(seg_df,initialCN)
    CN_base <- initialCN$CN[initialCN$base==1]
    CN_min <- initialCN$CN[1]
    CNest <- ploidyEst1(initialCN,seg_df)
    
   
      if(is.null(outplot_name)){
        outplot_name_i= paste0("bulkCNV_",Cname)
      }else{
        outplot_name_i= paste0(outplot_name,"_",Cname)
      }
      
      plt <- cellLineSeg_v3(CNVres,output,outdir,outplot_name_i, ggarrange_width=c(10,5),
                            height = 3,width=10,outPlot=FALSE)
      
      p_ls[[Cname]]<- plt$ggarranged_p
      CNVres_total[[Cname]] <- output 
      
      rm(CNVres,output,seg.dat_new,mv_dist)
  }
  } 
  if(outputFigure){
    pcom <- cowplot::plot_grid(plotlist=p_ls,ncol = 1,align="v")
    ggsave(paste0(outdir, "/",outplot_name,".pdf"),pcom, width=14, height=3*length(ls_binseg),device = pdf,bg="white")

}
    
  
  ###------------------------------------###
  #merge segment ratio and integer CNV results
  ###------------------------------------###
  g <- getBSgenome(genome, masked=FALSE)
  res_ls <- list()
  for(groupi in names(CNVres_total)){
    integCNV_Ci <- CNVres_total[[groupi]]
    if(all(!is.na(integCNV_Ci))){
      
    CNseg.dat <- integCNV_Ci$seg.dat
    CNseg.dat$segName <- paste(CNseg.dat$chr,CNseg.dat$start,CNseg.dat$end,sep="_")
    CNseg.dat_slim <- CNseg.dat[,c("segName","CNV","integerCNV"),drop=F]
    
    #groupi <- paste0(clust_i,"_adjacent")
    df_seg_C1 <- ls_binseg[[groupi]]
    
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
    
    
    if(true_loc){
      df_seg_C1$segLen <- as.numeric(df_seg_C1$segEnd)-as.numeric(df_seg_C1$segStrat)
    }else{
      df_seg_C1$segLen <- as.numeric(df_seg_C1$length_seg)
    }   
    
    segLen_sum <- sum(unique(df_seg_C1[,c("segID","segLen")]$segLen))
    df_seg_C1$wei <- df_seg_C1$segLen/segLen_sum
    
    res_ls[[groupi]] <- df_seg_C1
    }
  }

  return(res_ls)
  # saveRDS(res_ls,paste0(outdir,"/CNV_res.rds"))
}



#' 

#' @title mse_loss()
#' @description the Mean Squared Error (MSE) of integer CNV estimation to the mapping ratio 
#' @param seg.cnv.df data.frame with SegName and SegMean in columns
#'  segName  SegMean ratio_map
  # chr1_827077_120942799 1.426278      1.47
  # chr1_148951389_248907553 1.421164      1.47
  # chr2_19449_91766890 1.451020      1.47
  # chr2_95076467_108240173 1.368492      1.47
#' @export

mse_loss<- function(seg.cnv.df,
                     y_true.term = NULL,
                     y_pred.term = NULL,
                     weighted=TRUE,
                     weight.term = NULL){
  if(is.null(y_true.term)){y_true.term="SegMean"}
  if(is.null(y_pred.term)){y_pred.term="ratio_map"}
  seg.cnv.df <- seg.cnv.df[(!is.na(seg.cnv.df[,y_true.term])) & (!is.na(seg.cnv.df[,y_pred.term])),,drop=F]
  if(nrow(seg.cnv.df)>0){
    y_true <- as.numeric(seg.cnv.df[,y_true.term])
    y_pred <- as.numeric(seg.cnv.df[,y_pred.term])
    if(weighted){
      if(is.null(weight.term)){weight.term="length_seg"}
      if(!any(grepl(weight.term,colnames(seg.cnv.df)))){
        seg.cnv.df[,weight.term] <- as.numeric(seg.cnv.df$end) -as.numeric(seg.cnv.df$start) 
      }
      len <- as.numeric(seg.cnv.df[,weight.term])
      mse <- sum(len*(y_true - y_pred)^2)/sum(len)
    }else{
      mse <- mean((y_true - y_pred)^2,na.rm = TRUE)
    } 
  }else{
    stop(paste0("No data in colums ",y_true.term, " or ",y_pred.term,"!"))
  }
  return(mse)
}

#' @title peakIndex()
#' @export
peakIndex <- function(df_seg_C1,initialCN,index_col="SegMean"){
  colnames(initialCN)[grepl("ratio",colnames(initialCN),ignore.case = T)] <- "ratio"
  colnames(initialCN)[grepl("CN",colnames(initialCN),ignore.case = T)] <- "CN"
  initialCN <- initialCN[order(initialCN$ratio),,drop=F]
  initialCN$base = 0
  pp <- density(df_seg_C1[,index_col],width=0.1)
  dens <- pp$y
  ratioX <- pp$x
  a <- ratioX[which.max(dens)]
  dis <- na.omit(abs(a-initialCN$ratio))
  if (min(dis) < 0.05){
    initialCN$base[which.min(dis)]=1
    return(initialCN)
  }else{
    base0 <- a
    baseCN <- initialCN$CN[which.min(dis)]
    if(nrow(initialCN)>1){
      delt= abs(initialCN$ratio[1]-initialCN$ratio[2])/(initialCN$CN[2]-initialCN$CN[1])
      initialCN$ratio <- base0 + (initialCN$CN-baseCN)*delt
      initialCN$base[which(initialCN$ratio==base0)]=1
    }else{
      initialCN$base=1
    }

    
    return(initialCN)
  }
}

#' @title ploidyEst1()
#' @export
ploidyEst1 <- function(CNest,df_seg_C1){
  ##average ploidy
  ploidy <- round(CNest$ratio[CNest$base==1]*2)
  if (ploidy != CNest$CN[CNest$base==1]){
    d= ploidy-CNest$CN[CNest$base==1]
    CNest <- cbind(CNest,CNest$CN+d)
  }
  ##minimal ploidy
  base0 <- CNest$ratio[CNest$base==1]-(CNest$CN[CNest$base==1]-1)*abs(CNest$ratio[1]-CNest$ratio[2])
  base1 <- base0
  mratio <- df_seg_C1$SegMean[df_seg_C1$SegMean<base1]
  dis <- abs(base1-mratio)
  prob <- length(dis[dis>0.2])/dim(df_seg_C1)[1]
  while(1){
    if (prob < 0.01){
      break
    }else{
      base1 <-  base1-abs(CNest$ratio[1]-CNest$ratio[2])
      mratio <- df_seg_C1$SegMean[df_seg_C1$SegMean<base1]
      dis <- abs(base1-mratio)
      prob <- length(dis[dis>0.2])/dim(df_seg_C1)[1]
    }
  }
  if (base1 != base0){
    a <- (CNest$ratio-base1)/abs(CNest$ratio[1]-CNest$ratio[2])+1
    CNest <- cbind(CNest,a)
  }
  ##overall ploidy
  intePloidy <- round(CNest$ratio*2)
  k <- which(CNest[CNest$base==1,]==ploidy)
  CNindex <- unlist(sapply(k, function(j,CNest,intePloidy){
    fre <- sum(intePloidy==CNest[,j])/length(intePloidy)
    if (fre < 0.75){
      return(j)
    }
  },CNest,intePloidy))
  if (!is.null(CNindex)){
    newCN <- do.call(cbind,lapply(CNindex, function(j,CNest,intePloidy){
      indexorder <- order(abs(CNest$ratio*2-intePloidy))
      ploidyIndex = which(CNest$base==1)
      indexorder = indexorder[indexorder!=ploidyIndex]
      delt=abs(CNest$ratio[1]-CNest$ratio[2])
      newCN <- c()
      finish <- c()
      for (ele in indexorder){
        if (!(ele %in% finish)){
          deltG <- abs(CNest[ele,j]-ploidy)
          if (deltG > 1){
            a <- rep(NA,length=dim(CNest)[1])
            a[ploidyIndex]=ploidy
            base0 = intePloidy[ele]
            base1Ratio = 0
            if (base0 < ploidy & (ploidyIndex+deltG) <= dim(CNest)[1]){
              base1Ratio <- CNest$ratio[ploidyIndex+deltG]
            }else if (base0 > ploidy & (ploidyIndex-deltG) >=1){
              base1Ratio <- CNest$ratio[ploidyIndex-deltG]
            }
            if (base1Ratio != 0){
              if (base1Ratio %in% CNest$ratio & abs(ploidy-intePloidy[which(CNest$ratio == base1Ratio)])==abs(ploidy-base0)){
                a[ele]=base0
                a[which(CNest$ratio == base1Ratio)]=intePloidy[which(CNest$ratio == base1Ratio)]
              }
              if (length(a[!is.na(a)])>1){
                newCN=cbind(newCN,a)
              } 
            }
          } 
        }
        finish <- c(finish,ele,ploidyIndex+ploidyIndex-ele)
      }
      return(newCN)
    },CNest,intePloidy))
    if (!is.null(newCN)){
      CNest <- cbind(CNest,newCN)
    }
  }
  if (dim(CNest)[2]>3){
    CNest <- CNest[,c(1,3,2,4:dim(CNest)[2])]
  }else{
    CNest <- CNest[,c(1,3,2)] 
  }
  colnames(CNest) <- c("ratio","base",paste("CN",c(1:(dim(CNest)[2]-2)),sep=""))
  return(CNest)
}


#' @title CNplot()
#' @export
CNplot <- function(df_seg_C1,CNest,clust_i,ylim=NULL){
  require(ggplot2)
  if(is.null(ylim)){
    ylim <- round(quantile(df_seg_C1$SegMean,c(0,1)),2)
    ylim <- c(max(0,(ylim[1]-0.5)),ylim[2]+0.5)
  }
  breakshist <- seq(min(ylim),max(ylim),length.out=80)
  integers_map <- CNest$ratio
  if(dim(CNest)[2]>3){
    integers_label <- apply(CNest[,3:dim(CNest)[2]],1,function(x){paste(x,collapse = "     ")})
  }else{
    integers_label <- CNest$CN
  }
  ###histogram for bulk-level segment ratio
  while (length(dev.list()) > 0) {
    dev.off()
  }
  hist_bulk <- ggplot(df_seg_C1,aes(x=SegMean,fill = ..x..))+
    geom_histogram(aes(y=..density.., weight = wei),breaks=breakshist,color="darkgrey",fill="grey")+
    geom_vline(xintercept = integers_map,col="grey",linetype = "dotted")+
    labs(x = "", y = "",fill="")+
    ggtitle(paste0("C",clust_i," ratio"))+
    coord_cartesian(xlim=ylim)+
    theme_bw() +
    theme_classic()+
    coord_flip(xlim=ylim)+
    scale_x_continuous(
      sec.axis = sec_axis( trans = ~.,breaks=integers_map, labels =integers_label, name=""))+  
    theme(plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
          # legend.position = "bottom",
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          plot.title = element_text(size = 15),
          axis.title=element_text(size=15),
          axis.text = element_text(size = 10))+
    labs(x = "", y = "",fill="")
  
  return(hist_bulk)
  
}


#' @title bestCNVres.select()
#' @description selection the best integer CNV result from  multiple clusters based on the the Mean Squared Error (MSE) of integer CNV estimation to the mapping ratio 
#' @param inteCNV_res list of integer CNV result with length of N clusters
#' @export
bestCNVres.select <- function(inteCNV_res){
  mse <- do.call(c,lapply(inteCNV_res, function(x){
    #cnv <- unique(x[,c("segName","SegMean","ratio_map")])
    cnv <- x
    mse <- mse_loss(cnv)
  }))
  clust_i <- names(inteCNV_res)[which.min(mse)]
  df_seg_C1 <- inteCNV_res[[clust_i]]
  df_seg_C1 <- df_seg_C1[!is.na(df_seg_C1$ratio_map),]
 # update ploidy choice
  initialCN <- as.data.frame(unique(cbind(df_seg_C1$ratio_map,df_seg_C1$integerCNV)))
  colnames(initialCN) <- c("ratio","CN")
  initialCN <- initialCN[order(initialCN$CN),]
  initialCN$base=0
  initialCN1 <- peakIndex(df_seg_C1,initialCN)
  
  if (dim(initialCN)[1]!=1){
    disupdate <- DeltSeek(df_seg_C1,initialCN1) ##edit: 1216
    CNest <- FinalCN(initialCN1,disupdate,df_seg_C1) ##edit: 1216
    CNest <-ploidyEst1(CNest,df_seg_C1)
   # CNest <-ploidyEst1(initialCN1,df_seg_C1)
  }else{
    CNest <- initialCN[,1:2]
  }
  while (length(dev.list()) > 0) {
    dev.off()
  }
   hist_bulk2 <- CNplot(df_seg_C1,initialCN[,1:2],paste0(clust_i,"(best)"))
  if(ncol(CNest)>3){
    hist_bulk1 <- CNplot(df_seg_C1,CNest,paste0(clust_i,"(ploidy choice)"))
    p <- grid.arrange(hist_bulk2, hist_bulk1, nrow = 1)
  }else{
    p <- grid.arrange(hist_bulk2, grid::nullGrob(), nrow = 1)
  }
  while (length(dev.list()) > 0) {
    dev.off()
  }
  df_seg_C1 <- df_seg_C1[order(as.numeric(df_seg_C1$chrom),as.numeric(df_seg_C1$Start)),]
  return(list(cluster.select=clust_i,
              inteCNV_res = df_seg_C1,
              ploidy = CNest,
              CNV_histogram = p
              ))
}




seg_plot.CNV <- function(inteCNV_res,
                         outdir=".",
                         color_hist = TRUE,
                         plot_colors = NULL,
                         color_limit = c(0,3),
                         ylim = NULL,
                         add_score = TRUE,
                         outplot_name="CNV_ratio_mapping"){
  suppressPackageStartupMessages({library(cowplot)})
  if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
  if(is.null(plot_colors)){
    #colorss <- c("black","blue","grey40","grey60","green3","green","greenyellow", "yellow", "yellow1", "yellow2", "orange","orange2","orangered","red","red3","brown","#780000")
    colorss <- c("black","blue","grey50","green", "yellow", "orange","brown")
  }else{
    colorss <- plot_colors
  }
  while (length(dev.list()) > 0) {
    dev.off()
  }

  p_ls <- list()
  clusterID <- names(inteCNV_res)
  # mse.all <- do.call(c,lapply(inteCNV_res, function(x){
  #   mse <- mse_loss(x)
  # }))
  # bar_df <- data.frame(cluster=clusterID,mse=mse.all)
  # bar_df$score <-  round(-log(bar_df$mse),2)
  # score_col_lim <- c(min(bar_df$score),max(bar_df$score))
  for (clust_i in clusterID){
    df_seg_C1 <- inteCNV_res[[clust_i]]
    df_seg_C1 <- df_seg_C1[,!grepl("seg_SD|seg_score",colnames(df_seg_C1)),drop=F]
    df_seg_C1 <- na.omit(df_seg_C1)
    #df_seg_C1.s <- unique(df_seg_C1[,c("segName","SegMean","ratio_map","integerCNV")])
    #df_seg_C1 <- df_seg_C1[!is.na(df_seg_C1$ratio_map),]
    if(!any(grepl("chrom",colnames(df_seg_C1)))){
      chromnum=sapply(strsplit(df_seg_C1$binID,'-|_|:'),'[',1)
      if(grepl('chr',chromnum[1])){
        chromnum=sapply(strsplit(chromnum,'hr'),'[',2) ## if chrommnum with chr
      }else{
        chromnum=chromnum
      }
      chromnum <- gsub("X","23",chromnum)
      chromnum <- gsub("Y","24",chromnum)
      df_seg_C1$chrom <- chromnum
    }
    df_seg_C1 <- df_seg_C1[order(as.numeric(df_seg_C1$chrom),as.numeric(df_seg_C1$Start)),]
    # if(add_score){
    #   bar <- bar_df[bar_df$cluster==clust_i,,drop=F]
    #   p_bar <- ggplot(data=bar,aes(x=cluster, y=score,fill=score)) +
    #     geom_bar(stat="identity", width=0.3)+
    #     ylim(c(0,ceiling(max(bar_df$score))))+
    #     coord_flip()+
    #     theme_classic()+ 
    #     scale_fill_gradientn(name = "",colours = c( "#d4ea62","#ed711e","#e9295c"),
    #                          limits  = score_col_lim,
    #                          oob = scales::squish)
    #   # legend <- get_legend( p_bar + guides(color = guide_legend(nrow = 1)) +
    #   #     theme(legend.position = "bottom"))   
    # }
    
    if(is.null(ylim)){
      ylim <- round(quantile(df_seg_C1$SegMean,c(0,1)),2)
      ylim <- c(max(0,(ylim[1]-0.5)),ylim[2]+0.5)
    }
    if(add_score){
      mse = bar_df[bar_df$cluster==clust_i,"mse",drop=F]
      mse <- signif(mse,3)
      add_text <- paste0("C",clust_i,": MSE=",mse)
    }else{
      add_text <-  paste0("C",clust_i)
    }
    pl_scRatio <- seg_plot(df_seg_C1,name.data=add_text,
                           value.bin="binRatio",
                           value.segment="SegMean",
                           ylab="Relative CN",
                           plot_dots = F,
                           ylim=ylim,
                           plot_seg = T,
                           seg_col = "#DF2935",
                           color_hist = TRUE,
                           plotDir=outdir,
                           outPlot=F,
                           plot_hist=T)
    data <-pl_scRatio$p1$data
    data$integerCNV <-  df_seg_C1$integerCNV
    data$ratio_map <- df_seg_C1$ratio_map
    p1 <-  pl_scRatio$p1 +
      geom_segment(
        data = data,
        mapping = aes(x=segStart, y=ratio_map, xend=segEnd, yend=ratio_map),
        col="#080708",alpha = 0.8,na.rm =T,size = 1,
        inherit.aes = FALSE)
    breakshist <- seq(min(ylim),max(ylim),length.out=80)
    label <- unique(data[,c("integerCNV","ratio_map")])
    my_theme <- theme(plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
                      plot.title = element_text(size = 15), #legend.position = "bottom",
                      axis.line.y.right = element_blank(),
                      axis.ticks.y.right = element_blank(),
                      axis.title=element_text(size=15),
                      axis.text = element_text(size = 12),
                      axis.text.y.right = element_text(color = "#080708"))
    if(color_hist){
    p3=ggplot(df_seg_C1,aes(x=SegMean,fill = ..x..))+
      geom_vline(xintercept = label$ratio_map,col="grey",linetype = "dotted")+  #a$value.segment.map
      geom_histogram(aes(y=..density.., weight = length_seg),breaks=breakshist)+ #color="#E69F00",fill="#DF2935",
      ggtitle("")+
      coord_cartesian(xlim=ylim)+
      theme_bw() +
      theme_classic()+
      coord_flip(xlim=ylim)+
      scale_fill_gradientn(name = "",colours = colorss,#c("#003049", "#669bbc","#7e7f83","#c1121f","#780000"),
                           limits  = color_limit,
                           oob = scales::squish)+
      scale_x_continuous(
        sec.axis = sec_axis( trans = ~.,breaks=label$ratio_map, labels =label$integerCNV, name="")
      )+
      my_theme+
      labs(x = "", y = "",fill="")
    
    # legend2 <- get_legend( p3 + guides(color = guide_legend(nrow = 1)) +
    #     theme(legend.position = "bottom"))   
    
    }else{
      p3=ggplot(df_seg_C1,aes(x=SegMean))+
        geom_vline(xintercept = label$ratio_map,col="grey",linetype = "dotted")+  #a$value.segment.map
        geom_histogram(aes(y=..density.., weight = length_seg),breaks=breakshist)+ #color="#E69F00",fill="#DF2935",
        ggtitle("")+
        coord_cartesian(xlim=ylim)+
        theme_bw() +
        theme_classic()+
        coord_flip(xlim=ylim)+
        scale_x_continuous(
          sec.axis = sec_axis( trans = ~.,breaks=label$ratio_map, labels =label$integerCNV, name="")
        )+
        my_theme+
        labs(x = "", y = "",fill="")
    }
    
    # pp=ggarrange(p1, p3, align ="h",ncol = 2, nrow = 1,widths = c(15,5),heights=3 )+
    #   theme(plot.margin=unit(c(0,0,0,0),"cm"))
    
    pp <- cowplot::plot_grid(p1+theme(plot.margin=unit(c(0,0,0,0.5),"cm")),
                            p3+theme(plot.margin=unit(c(0,0,0,0),"cm")),
                            axis="tblr",align="h",
                            rel_widths=c(15,5),rel_heights=c(3,3))
    
    # if(add_score){
    #   null <- ggplot()+theme_void()+
    #     theme(plot.margin=unit(c(0,0,-1,0),"cm"))
    #   #ppb <- ggarrange(p_bar.com,pp,nrow = 2,heights=c(1,3),align = "v",axis = "l")
    #   
    #  # legend.all <- cowplot::plot_grid(legend, legend2, ncol = 2, align = "h")+theme(plot.margin=unit(c(0,0,-1,0),"cm"))
    #   ppb <- plot_grid(null,p_bar+ theme(legend.position = "none")+
    #                      theme(plot.margin=unit(c(0.5,2,-1,0),"cm")),
    #                    p1+theme(plot.margin=unit(c(0,0,0,0.5),"cm")),
    #                     p3+theme(plot.margin=unit(c(0,0.5,0,0),"cm")),
    #                     ncol = 2, align = 'hv',axis="bl",rel_heights = c(1.5,3),rel_widths=c(15,5))
    #   pp <- ppb
    # }
    
    p_ls[[clust_i]]<- pp
  }
  pcom <- cowplot::plot_grid(plotlist=p_ls,ncol = 1,align="v")
  ggsave(paste0(outdir, "/",outplot_name,".png"),pcom, width=14, height=4*length(clusterID),device = png,bg="white",dpi = 300)
  while (length(dev.list()) > 0) {
    dev.off()
  }
  return(bar_df)
}  


#' @title CNbaseline.fill()
#' @description fill the blank CN level based on known delt and CNs 
#' @param CNbase data.frame with columns of 'ratio' and 'CN'
#' @export
CNbaseline.fill <- function(CNbase,delt=NULL){
  CNbase <- as.data.frame(CNbase)
  colnames(CNbase)[1:2] <- c("ratio","CN")
  CNbase <- CNbase[order(CNbase$CN),]
  cn.min <- min(CNbase$CN)
  cn.max <- max(CNbase$CN)
  
  if(nrow(CNbase)>1){
  dis <- diff(CNbase$CN)
  if(dis[1]==1){
    delt <- (CNbase$ratio[2]-CNbase$ratio[1])
  }else{
    fold = dis[1]
    delt <- (CNbase$ratio[2]-CNbase$ratio[1])/fold
  }
  
  }else if(is.null(delt)){
    stop("Can not fill CNs with only one CN level provided!")
  }
  # if(cn.min>1){
  #   cn.min <- cn.min-1
  #   ratio.min <- CNbase$ratio[1]-delt
  # }else{
    ratio.min <- CNbase$ratio[1]
  # }
  ratio.max <- CNbase$ratio[nrow(CNbase)]
  ratio <- seq(ratio.min,ratio.max,by=delt)
  CN <- seq(cn.min,cn.max,by=1)
  
  CNbase_new <- data.frame(ratio=ratio,CN=CN)

  if(any(grepl("base",colnames(CNbase)))){
    ratiobase = CNbase$ratio[CNbase$base==1]
    CNbase_new$base = 0
    CNbase_new$base[CNbase_new$ratio==ratiobase] <- 1
  }
    
  return(CNbase_new)
}


