#Figures
#CNVres is cell level process data
#seg.dat1 == cellIntegerCNV
#library("ComplexHeatmap")
library('RColorBrewer')
suppressMessages(library(dplyr))

#' @title segPlot()
#' @export

segPlot<-function(CNVres,cellIntegerCNV,delta,outdir,fname,ylim=NULL){
  plotDir<-outdir
  if(!dir.exists(plotDir)){dir.create(plotDir,recursive=T)}
  png(paste0(plotDir,"/",fname,".png"),width = 16, height = 4,units = 'in',res= 300)
  #cell_level
  ##cellres=CNVres
  seg.dat=CNVres$seg
  chromoPos=CNVres$chromoPos
  endSum <- seg.dat %>% group_by(chr) %>%
    dplyr::summarise(chrline=max(absend))%>% as.data.frame
  endSum <- endSum$chrline
  endSum <- endSum[order(endSum)]
  CNV=CNVres$CNV
  #box()
  par(fig=c(0,0.8,0,1))
  if(is.null(ylim)){ylim=c(0,max(CNV$ratio)*1.3)}
  plot(chromoPos$abspos,CNV$ratio,xlim=c(min(chromoPos$abspos),max(chromoPos$abspos)),ylim=ylim,
    col=rgb(0.5,0.5,1,alpha = 0.5),pch=20,cex=0.4,ylab="Ratio",axes = F,xlab="chrom",main=fname)
  segments(seg.dat$abspos,as.numeric(seg.dat$ratio),seg.dat$absend,as.numeric(seg.dat$ratio),col = rgb(0.6,0,0.2,alpha = 0.7),lwd=2)
  segments(seg.dat$absend[1:(length(seg.dat$absend)-1)],as.numeric(seg.dat$ratio[1:(length(seg.dat$absend)-1)]),seg.dat$abspos[2:(length(seg.dat$absend))],as.numeric(seg.dat$ratio[2:(length(seg.dat$absend))]),col = rgb(0.7,0,0.2,alpha = 0.7),lwd=2)
  abline(v=endSum)
  axis(side=2)
  #box()
  labelss<-c(1:22,"X")
  pos<-c(0,endSum[1:(length(endSum)-1)])+(endSum-c(0,endSum[1:(length(endSum)-1)]))/2
  axis(side=1,labels = labelss,at=pos)
  
  #seg.dat1=output$seg.dat
  seg.dat1=cellIntegerCNV
  #cellIntegerCNV
  #print("aa")
  seg.dat1$CNV=seg.dat1$integerCNV*delta
  CNVfrac=aggregate(seg.dat1$w, by=list(Category=seg.dat1$CNV), FUN=sum)
  CNVfrac$frac=CNVfrac$x/sum(CNVfrac$x)
  #print(dim(CNVfrac))
  #print(dim(seg.dat1))
  CNVfrac$integerCNV=as.numeric(levels(factor(seg.dat1$integerCNV)))
  #print("here")
  par(fig=c(0.75,1,0,1),new=TRUE)
  plot(0,0,xlim=c(0,max(CNVfrac$frac)),ylim=c(0,max(CNV$ratio)*1.3),col="white",pch=20,ylab="",axes = F,xlab="")
  segments(0,CNVfrac$Category,CNVfrac$frac,CNVfrac$Category,col = "purple",lwd=4)
  axis(side=4,at=CNVfrac$Category,labels = CNVfrac$integerCNV)
  axis(side=2)
  box()  
  par(mfrow=c(1,1))
  box() 
  dev.off()
}


#' @title cellLineSeg()
#' @export

cellLineSeg<-function(CNVres,integerCNV,outdir,fname,ylim=NULL){
  plotDir<-outdir
  if(!dir.exists(plotDir)){dir.create(plotDir,recursive=T)}
  png(paste0(plotDir,"/",fname,".png"),width = 16, height = 4,units = 'in',res= 300)
  seg.dat <- integerCNV$seg.dat
  chromoPos=CNVres$chromoPos
  endSum <- seg.dat %>% group_by(chr) %>%
    dplyr::summarise(chrline=max(absend))%>% as.data.frame
  endSum <- endSum$chrline
  endSum <- endSum[order(endSum)]
  CNV=CNVres$CNV
  #box()
  par(fig=c(0,0.8,0,1))
  if(is.null(ylim)){ylim=c(0,max(CNV$ratio)*1.3)}
  plot(chromoPos$abspos,CNV$ratio,xlim=c(min(chromoPos$abspos),max(chromoPos$abspos)),ylim=ylim,
       col=rgb(0.5,0.5,1,alpha = 0.5),pch=20,cex=0.4,ylab="Ratio",axes = F,xlab="chrom",main=fname)
  segments(seg.dat$abspos,as.numeric(seg.dat$ratio),seg.dat$absend,as.numeric(seg.dat$ratio),col = rgb(0.6,0,0.2,alpha = 0.7),lwd=2)
  segments(seg.dat$absend[1:(length(seg.dat$absend)-1)],as.numeric(seg.dat$ratio[1:(length(seg.dat$absend)-1)]),seg.dat$abspos[2:(length(seg.dat$absend))],as.numeric(seg.dat$ratio[2:(length(seg.dat$absend))]),col = rgb(0.7,0,0.2,alpha = 0.7),lwd=2)
  abline(v=endSum)
  axis(side=2)
  #box()
  labelss<-c(1:22,"X")
  pos<-c(0,endSum[1:(length(endSum)-1)])+(endSum-c(0,endSum[1:(length(endSum)-1)]))/2
  axis(side=1,labels = labelss,at=pos)
  
  seg.dat1=output$seg.dat
  CNVfrac=aggregate(seg.dat1$w, by=list(Category=seg.dat1$CNV), FUN=sum)
  CNVfrac$frac=CNVfrac$x/sum(CNVfrac$x)
  CNVfrac$integerCNV=as.numeric(levels(factor(seg.dat1$integerCNV)))
  
  par(fig=c(0.75,1,0,1),new=TRUE)
  plot(0,0,xlim=c(0,max(CNVfrac$frac)),ylim=c(0,max(CNV$ratio)*1.3),col="white",pch=20,
       ylab="",axes = F,xlab="")
  segments(0,CNVfrac$Category,CNVfrac$frac,CNVfrac$Category,col = "purple",lwd=4)
  axis(side=4,at=CNVfrac$Category,labels = CNVfrac$integerCNV)
  axis(side=2,at=0:4,labels=0:4)
  box()
  par(mfrow=c(1,1))
  box() 
  
  dev.off()

}


#' @title cellLineSeg_v3()
#' @description plot dots and histogram of segment ratio at true chromosome location, 
#' remove bad segment from histogram and hist based on the true genomic length of segments,set a global color
#' @param CNVres list out put from dataProcess()function
#' @param integerCNV list output from doCNV1_v2() function
#' @export

cellLineSeg_v3<-function(binRatio.df = NULL,
                         integerCNV.df=NULL,
                         CNVres=NULL,integerCNV=NULL,
                         outdir,fname,ylim=NULL,
                         genome="hg38",
                         ylab="Relative CN",
                         plot_seg = TRUE,
                         color_seg_gradient = TRUE,
                         seg.color = "#0d3b66",
                         ggarrange_width=c(10,5),#c(10,5),
                         height = 2,width=15,
                         plot_dots = TRUE,
                         color_dot =FALSE,
                         plot_hist = TRUE,
                         color_hist_gradient = TRUE,
                         hist_color = NULL,
                         color_limit = c(1,6),
                         label_gene=NULL,
                         label_gene_color="blue",
                         annotations=NULL,
                         outPlot=TRUE,length.out=100,
                         value.bin="ratio",
                         value.segment="SegMean",
                         label_CN=TRUE
){
  suppressMessages({
    library(BSgenome)
    library(ggplot2)
    library(dplyr)
    library(ggpubr)
  })
  g <- getBSgenome(genome, masked=FALSE)
  chrLine <- data.frame(chrom=1:24, length=seqlengths(g)[1:24])
  if(is.null(integerCNV.df)){
    seg.dat <- integerCNV$seg.dat
  }else{
    seg.dat <- integerCNV.df
  } 
  colnames(seg.dat)[grepl("^CNV$",colnames(seg.dat))] <- "relativeCN"
  colnames(seg.dat)[grepl("^integerCNV$",colnames(seg.dat))] <-"integerCN"
  
  if(is.null(binRatio.df)){
    data.bin=CNVres$CNV
  }else{
    data.bin <- binRatio.df
  }
  colnames(data.bin)[grepl("chrom|Chromosome",colnames(data.bin))] <- "chrom"
  colnames(data.bin)[grepl("^Start$",colnames(data.bin),ignore.case = TRUE)] <- "start"
  colnames(data.bin)[grepl("^binRatio$",colnames(data.bin),ignore.case = TRUE)] <- "ratio"
  ploidy <- sum(seg.dat$integerCN*seg.dat$w)/sum(seg.dat$w)
  
  chr_name=as.character(unique(data.bin$chrom))
  #chromnum=sapply(strsplit(rownames(CNV),'-|_|:'),'[',1)
  chromnum <- data.bin$chrom
  if(grepl('chr',chromnum[1])){
    chromnum=sapply(strsplit(chromnum,'hr'),'[',2) ## if chrommnum with chr
  }else{
    chromnum=chromnum
  }
  chromnum <- gsub("X","23",chromnum)
  chromnum <- gsub("Y","24",chromnum)
  data.bin$chrom <- chromnum
  chrLine <- chrLine[chrLine$chrom %in% as.character(unique(chromnum)),,drop=F] 
  genomeLength <- sum(chrLine)
  if(nrow(chrLine)>1){
    chrPosStart = c(0,chrLine[1:(nrow(chrLine)-1),2])
  }else(
    chrPosStart = 0
  )
  chrLine$chrPosStart_abs <- cumsum(chrPosStart)
  data.bin <- merge(data.bin,chrLine,by=intersect(colnames(data.bin),colnames(chrLine)),all.x=T)
  data.bin <- data.bin[order(as.integer(data.bin$chrom),data.bin$start),]
  data.bin$seg_start_abs <- as.numeric(sapply(strsplit(data.bin$segName,'_|-|:'),'[',2))+as.numeric(data.bin$chrPosStart_abs)
  data.bin$seg_end_abs <- as.numeric(sapply(strsplit(data.bin$segName,'_|-|:'),'[',3))+as.numeric(data.bin$chrPosStart_abs)
  
  seg.dat$segName <- paste(seg.dat$chr,seg.dat$start,seg.dat$end,sep="_")
  if(any(grepl("SegMean|ratio",colnames(seg.dat)))){
    colnames(seg.dat)[grepl("SegMean|ratio",colnames(seg.dat))] <- "SegMean"
    seg.dat_slim <- seg.dat[,c("segName","SegMean","relativeCN","integerCN"),drop=F]
    data.bin <- data.bin[,!grepl("relativeCN|integerCN|SegMean",colnames(data.bin))]
  }else{
    seg.dat_slim <- seg.dat[,c("segName","relativeCN","integerCN"),drop=F]
    data.bin <- data.bin[,!grepl("relativeCN|integerCN",colnames(data.bin))]
  }
  
  data.bin <- merge(data.bin,seg.dat_slim,by="segName",all.x=T)
  data.bin <- data.bin[!is.na(data.bin$integerCN),,drop=F]
  
  data.bin <- data.bin[,!grepl("segID",colnames(data.bin)),drop=F]
  data.bin <- data.bin %>%
    group_by(segName) %>%
    dplyr::mutate(segID = cur_group_id()) %>%
    ungroup()%>% as.data.frame()
  data.bin <- data.bin[order(as.integer(data.bin$chrom),data.bin$chrPosStart_abs),]
    
  maploc <- data.bin$chrPosStart_abs+data.bin$start
  if(grepl('chr',chr_name[1])){
    chr_name <- sapply(strsplit(chr_name,'hr'),'[',2)
  }
  chr_name_nm <- gsub("X","23",chr_name)
  chr_name_nm <- gsub("Y","24",chr_name_nm)
  which_chrom <- which(as.character(chrLine$chrom) %in% chr_name_nm)
  
  breakss<-cumsum(chrPosStart)+chrLine[1:(nrow(chrLine)),2]*0.5
  breakss<-breakss[which_chrom]

  a<-data.frame(maploc=maploc,data.bin[,value.bin],data.bin[,value.segment],data.bin[,c("segID","seg_start_abs","seg_end_abs","relativeCN","integerCN","length_bin")])
  colnames(a)=c('axis.loc','value.bin','value.segment','segID',"segStrat","segEnd","ratio_map","integerCN","length_bin")
  # a$colour<- ifelse(a$integerCNV==2,"neutral",ifelse(a$integerCNV>2,"gain","loss"))
  # a$colour <- factor(a$colour,levels=c("loss","neutral","gain"))
  
  if(is.null(ylim)){ylim=c(0,max(data.bin$SegMean[is.finite(data.bin$SegMean)])+0.5)}
  breakshist <- seq(0,max(data.bin$SegMean[is.finite(data.bin$SegMean)]),length.out=100)
  a <- a %>%
    group_by(integerCN) %>%
    mutate(value.segment.map = median(value.segment))%>%
    as.data.frame()
  if(color_hist_gradient){
    if(is.null(hist_color)){
      colorss <- c("#023e8a","grey50", "#cb793a", "#9a031e","#6a040f","#370617")

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
  
  
  chr_name_label <- chr_name
  chr_name_label[as.character(chr_name_label)%in%c("19","21")] <- ""
  
  if(plot_dots){
    p1 <- ggplot2::ggplot(a) +
      geom_point(aes(x = axis.loc, y = value.bin), colour="darkgrey", alpha = 0.6, shape = 16,size=0.4) 
  }else{
    p1 <- ggplot2::ggplot(a,aes(x = axis.loc, y = value.bin)) 
  }
  
  p1 <- p1+
    theme_bw() +
    labs(x = "", y = ylab)+
    coord_cartesian(ylim = ylim)+
    ggtitle(fname)+
    geom_hline(yintercept =c(1,2,3),col="grey",linetype = "dashed")+
    #geom_vline(xintercept =chrLine$chrPosStart_abs,col="grey")+
    theme(plot.margin = unit(c(0.5,0,0,0.5), 'lines'), #c(top, right, bottom, left)
          panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.position = "none",plot.title = element_text(size = 15),
          axis.title=element_text(size=15),
          axis.text = element_text(size = 15))+
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE),breaks=breakss,labels=chr_name_label)
  odd_positions <- chrLine$chrPosStart_abs[seq(1, length(chrLine$chrPosStart_abs), by = 2)]
  even_positions <- chrLine$chrPosStart_abs[seq(2, length(chrLine$chrPosStart_abs), by = 2)]
  ## highlight region data
  rects <- data.frame(start=odd_positions, end=even_positions)
  p1 <- p1+
    geom_rect(data=rects, inherit.aes=FALSE, 
              aes(xmin=start, xmax=end, ymin=-Inf,ymax=Inf),
              color="transparent", fill="grey50", alpha=0.3)

  if(plot_seg){

    if(color_seg_gradient){

     if(max(a$integerCN,na.rm=TRUE)<=6){
       names(colorss) <- seq(1:6)
     }else{
       add_Ncolor <- as.numeric(max(a$integerCN,na.rm=TRUE))-6
       colorss_add <- colorRampPalette(c(colorss[6],"black"))(add_Ncolor)
       colorss <- c(colorss,colorss_add)
       names(colorss) <- seq(1,max(a$integerCN,na.rm=TRUE))
       color_limit <- c(1,max(a$integerCN,na.rm=TRUE))
     }
     
      p1 <- p1+
        #geom_line(aes(x = maploc,y=value.segment),col="red")+
        geom_hline(yintercept =1,col="grey",linetype = 2)+
        geom_segment(
          data = a, 
          mapping = aes(x=segStrat, y=value.segment, xend=segEnd, yend=value.segment,colour=integerCN), 
          na.rm =T,size = 1.5,
          inherit.aes = FALSE)+
        scale_color_gradientn(name = "",colours = colorss,#c("#003049", "#669bbc","#7e7f83","#c1121f","#780000"),
                              limits  = color_limit,
                              oob = scales::squish)
    }else{
      p1 <- p1+
        #geom_line(aes(x = maploc,y=value.segment),col="red")+
        geom_hline(yintercept =1,col="grey",linetype = 2)+
        geom_segment(
          data = a, 
          mapping = aes(x=segStrat, y=value.segment, xend=segEnd, yend=value.segment), 
          colour = seg.color,
          na.rm =T,size = 1,
          inherit.aes = FALSE)
    }
    
    
    
  }
  if(!is.null(label_gene)){
    gene_loc <- annotate_gene(label_gene,genome=genome,annotations=annotations)
    gene_chr <- as.character(seqnames(gene_loc))
    label_x <- chrLine[chrLine$chrom==gene_chr,"chrPosStart_abs"]+start(gene_loc)
    label_y <- unique(data.bin[data.bin$chrom==gene_chr&data.bin$seg_start_abs<=label_x&data.bin$seg_end_abs>=label_x,"SegMean"])[1]
    
    p1=p1+
      geom_point(aes(x = label_x, y = label_y), alpha = 1, shape = 17,size=1,colour=label_gene_color)+
      geom_text(aes(label_x,label_y,label = label_gene),color=label_gene_color,vjust = -1)+
      geom_segment(aes(x = label_x, y = label_y, xend = label_x, yend = label_y -1),
                   color = label_gene_color,  # 指定指引线颜色
                   linetype = "dashed")  
    
  }
  if(plot_hist){
    a$segLen <- as.numeric(a$segEnd)-as.numeric(a$segStrat)
    
    segLen_sum <- sum(unique(a[,c("segID","segLen")]$segLen))
    a$wei <- a$length_bin/segLen_sum
    integerCNVcolumn <- grep("integerCN",colnames(a))
    labs <- unique(a[integerCNVcolumn])
    labs <- labs[!is.na(labs[,1]),]
    
    if(plot_seg){
      p_x <- "value.segment"
    }else{
      p_x <- "value.bin"
    }
    
    clone_name <- strsplit(fname,":|_|-")[[1]][1]
    # a<- a[!is.na(a$integerCN),,drop=F]
    # a$integerCN <- factor(a$integerCN,levels = sort(unique(a$integerCN )))
    p3=ggplot(a,aes(x=get(p_x),weight = wei))+
        geom_histogram(aes(y=after_stat(density), ),breaks=breakshist)

    
    plot_build <- ggplot_build(p3)
    hist_data <- plot_build$data[[1]]
    a$bin <- cut(a$value.segment, breaks = breakshist, include.lowest = TRUE)
    a$bin_center <- with(a, (as.numeric(bin) - 0.5) * diff(breakshist)[1] + min(breakshist))
    
    CN_ano <- na.omit(unique(a[,c("bin_center","integerCN")]))
    find_nearest <- function(value, reference_vector) {
      reference_vector[which.min(abs(reference_vector - value))]
    }
    nearest_values <- sapply(hist_data$x, find_nearest, reference_vector = CN_ano$bin_center)
    CNindex <- match(nearest_values,CN_ano$bin_center)
    CNgroup <- CN_ano$integerCN[CNindex]
    hist_data$group <- factor(CNgroup)
    p4<- ggplot(hist_data, aes(x = x, y = density, fill = group))

    if(color_hist_gradient){
      if(max(CN_ano$integerCN)<=6){
        names(colorss) <- seq(1:6)
      }else{
        add_Ncolor <- as.numeric(max(CN_ano$integerCN))-6
        colorss_add <- colorRampPalette(c(colorss[6],"black"))(add_Ncolor)
        colorss <- c(colorss,colorss_add)
        names(colorss) <- seq(1,max(CN_ano$integerCN))
      }
     
      if(label_CN){
        p4=p4+
          geom_vline(xintercept = na.omit(unique(a$ratio_map)),col="grey",linetype = "dotted")
      }
      p3 <- p4 +
        geom_bar(stat = "identity", color = NA) +
        scale_fill_manual(values = colorss, breaks = breakshist)
      
      # p3 <- p3+  #a$value.segment.map
      #   geom_histogram(aes(y=after_stat(density), weight = wei),breaks=breakshist)+ #color="#E69F00",fill="wheat",
      #   #scale_fill_manual(values = colorss) 
      #   scale_fill_gradientn(name = "",colours = colorss,#c("#003049", "#669bbc","#7e7f83","#c1121f","#780000"),
      #                        limits  = color_limit,oob = scales::squish)
    }else{
      if(label_CN){
        p3=p3+
          geom_vline(xintercept = na.omit(unique(a$ratio_map)),col="grey",linetype = "dotted")
      }
      p3 <- p3+  #a$value.segment.map
        geom_histogram(aes(y=after_stat(density), weight = wei),color=colorss,fill=colorss,breaks=breakshist)
    }
    
    if(label_CN){
       p3=p3+ scale_x_continuous(
          sec.axis = sec_axis( trans = ~.,breaks=unique(a$ratio_map)[!is.na(unique(a$ratio_map))], labels =labs, name="")
        )
    }
    p3 <- p3 +
      ggtitle(paste0("Ploidy of ",clone_name,": ",round(ploidy,2)))+
      coord_cartesian(xlim=ylim)+
      theme_bw() +
      theme_classic()+
      coord_flip(xlim=ylim)+
      #guides(fill = guide_legend( reverse = TRUE)) +
      theme(plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
          plot.title = element_text(size = 15), #legend.position = "none",
          axis.line.y.right = element_blank(),
          axis.ticks.y.right = element_blank(),
          axis.title=element_text(size=15),
          axis.text = element_text(size = 12))+
    labs(x = "", y = "",fill="")
    pp=ggarrange(p1, p3, align ="h",ncol = 2, nrow = 1,widths =ggarrange_width,heights=3 )
    p_ls <- list(ggarranged_p=pp,p1=p1,p2=p3)
  }else{
    pp <- p1
    p_ls <- list(p1=p1)
  }
  if(outPlot){
    ggsave(paste0(outdir,"/",fname,".png"),pp,height = height,width=width,dpi = 300)
  }
  return(p_ls)
}


#' @title hist_seg()
#' @export
#df_seg_C1 is the ratio(segScore$seg_score_binLevel[[i]]); 
#integCNV_Ci is the CNV estimation result from doCNV1_v2()
hist_seg <- function(df_seg_C1,integerCNV,ylim=NULL,color_limit = c(0,3)){
  suppressMessages({
  library(BSgenome)
  library(dplyr)})
colorss <- c("black","blue","grey50","green", "yellow", "orange","brown")
CNseg.dat <- integerCNV$seg.dat
CNseg.dat$segName <- paste(CNseg.dat$chr,CNseg.dat$start,CNseg.dat$end,sep="_")
CNseg.dat_slim <- CNseg.dat[,c("segName","CNV","integerCNV"),drop=F]

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


breakshist <- seq(0,max(df_seg_C1$SegMean[is.finite(df_seg_C1$SegMean)]),length.out=80)

segSet_C1 = unique(df_seg_C1[,c("binID","segID","SegMean","segName","ratio_map","integerCNV"),drop=F]) ###

true_loc = TRUE
if(true_loc){
  df_seg_C1$segLen <- as.numeric(df_seg_C1$segEnd)-as.numeric(df_seg_C1$segStrat)
}else{
  df_seg_C1$segLen <- as.numeric(df_seg_C1$length_seg)
}   

segLen_sum <- sum(unique(df_seg_C1[,c("segID","segLen")]$segLen))
df_seg_C1$wei <- df_seg_C1$segLen/segLen_sum

hist_plt <- ggplot(df_seg_C1,aes(x=SegMean,fill = ..x..))+
geom_vline(xintercept = unique(df_seg_C1$ratio_map),col="grey",linetype = "dotted")+ 
  geom_histogram(aes(y=after_stat(density), weight = wei),breaks=breakshist)+
      ggtitle("Summary histogram")+
      coord_cartesian(xlim=ylim)+
      theme_bw() +
      theme_classic()+
      coord_flip(xlim=ylim)+
      scale_fill_gradientn(name = "",colours = colorss,
                           limits  = color_limit,oob = scales::squish)+
      scale_x_continuous(
        sec.axis = sec_axis( trans = ~.,breaks=unique(df_seg_C1$ratio_map), labels =unique(df_seg_C1$integerCNV), name="")
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



#' @title hist_custom()
#' @export

hist_custom <- function(data,x,
  color_hist=FALSE, color_col=NULL,row_filt_by=NULL,
  add_density=TRUE,scale_factor=1,adjust=1,
  log_y=FALSE,
  plotbreak=FALSE,y_cut=NULL,max_lim_y=NULL,stepby=NULL,
  cols_fill= NA,
  title_x=NULL,title_y="density",binwidth=0.1,xlim=NA,
  legend.position="none", legend.title="",
  title="",
  color_hist_by_value=TRUE){
  suppressMessages({
    library(ggplot2)
    library(dplyr)
    library(ggbreak) 
  })

  if(color_hist){
    if(color_hist_by_value){
      colorss <- c("#023e8a","grey50", "#cb793a", "#9a031e","#6a040f","#370617")
      if(!is.null(color_col)){colname_color=color_col}else{stop("Provide a colname as 'color_col'.")}
      if(max(data[,color_col])<=6){
        names(colorss) <- seq(1:6)
      }else{
        add_Ncolor <- as.numeric(max(data[,color_col]))-6
        colorss_add <- colorRampPalette(c(colorss[6],"black"))(add_Ncolor)
        colorss <- c(colorss,colorss_add)
        names(colorss) <- seq(1,max(data[,color_col]))
      }
    }else{
      colorss <- cols_fill
    }


    if(!log_y){

      p1 <- ggplot(data,aes(x=get(x)))+
        geom_histogram(mapping = aes(y = ..count.. / (sum(..count..) * ..width..),fill=factor(get(color_col))),colour = "white",  binwidth = binwidth)+
        scale_fill_manual(values = colorss)
      
    }else{
      p1 <- ggplot(data,aes(x=get(x)))+
            geom_histogram(aes(y=log10(..count.. / (sum(..count..) * ..width..)+1),fill=factor(get(color_col))),colour = "white",  binwidth = binwidth)+
        scale_fill_manual(values = colorss)
      
    }
  }
  if(plotbreak){
    if(!is.null(y_cut) &(!is.null(max_lim_y)) &(!is.null(stepby))){
      p1 <- p1 +
        scale_y_break(y_cut)+  ###
        scale_y_continuous(limits = c(0, max_lim_y), breaks = seq(0, max_lim_y, by = stepby))
  
    }else{
      print("Provide 'y_cut','max_lim_y', and 'stepby' for breaking y-axis.")
    }
  }
  if(add_density){
    if(!log_y){
      p1 <-  p1 +
            geom_density(data=data,aes(x=get(x),y=scale_factor*after_stat(density)),lwd = 0.5,linetype ="dashed",colour = "grey",adjust = adjust)
    }else{
       p1 <-  p1 +
            geom_density(data=data,aes(x=get(x),y=scale_factor*log10(after_stat(density)+1)),lwd = 0.5,linetype ="dashed",colour = "grey",adjust = adjust)
 
    }
  }

  p1 <- p1 +
      xlim(xlim)+
      labs(title=title,x=title_x,y=title_y)+
      theme_classic()+ 
      theme(plot.background=element_blank(),
        legend.position=legend.position,
        plot.title = element_text(hjust = 0.5),
              axis.title = element_text(color="black",size=7),
              axis.text=element_text(color="black",size=7),
              axis.line.y.right = element_blank(),
              axis.text.y.right=element_blank(),
              axis.ticks.y.right=element_blank(),
        axis.ticks=element_line(color="black",linewidth=0.5)
            )+
  guides(fill=guide_legend(title=legend.title))
  return(p1)

}

plot.format=theme(plot.background=element_blank(),panel.grid=element_blank(),panel.background=element_blank(),panel.border=element_rect(color="black",linewidth=0.5,fill=NA),axis.line=element_blank(),axis.ticks=element_line(color="black",linewidth=0.5),axis.text=element_text(color="black",size=7),axis.title=element_text(color="black",size=7),plot.title=element_text(color="black",size=7),legend.background=element_blank(),legend.key=element_blank(),legend.text=element_text(color="black",size=7),legend.title=element_text(color="black",size=7))


#' @title segPlot4CNV()
#'plot segment with estimated CNV
#' @param integerCNV a vector of integerCNV with the same length of CNVres.
#' @export
segPlot4CNV<-function(CNVres,integerCNV,outdir,fname){
plotDir<-outdir
if(!dir.exists(plotDir)){dir.create(plotDir,recursive=T)}
png(paste0(plotDir,"/",fname,".png"),width = 16, height = 4,units = 'in',res= 300)
seg.dat <- CNVres$seg
chromoPos=CNVres$chromoPos
endSum <- seg.dat %>% group_by(chr) %>%
  dplyr::summarise(chrline=max(absend))%>% as.data.frame
endSum <- endSum$chrline
endSum <- endSum[order(endSum)]
CNV=CNVres$CNV
CNV$integerCNV <- integerCNV
#box()
par(fig=c(0,1,0,1))
plot(chromoPos$abspos,CNV$integerCNV,xlim=c(min(chromoPos$abspos),max(chromoPos$abspos)),ylim=c(0,max(CNV$integerCNV)),
     col=rgb(0.5,0.5,1,alpha = 0.5),pch=20,cex=0.4,ylab="CNV",axes = F,xlab="chrom",main=fname)
abline(v=endSum)
axis(side=2)
#box()
labelss<-c(1:22,"X")
pos<-c(0,endSum[1:(length(endSum)-1)])+(endSum-c(0,endSum[1:(length(endSum)-1)]))/2
axis(side=1,labels = labelss,at=pos)
box()

dev.off()

}

#' @title segPlot1()
#' @export
#maybe not use
segPlot1<-function(celltmp,cellIntegerCNV,delta,outdir,fname){
  plotDir<-outdir
  if(!dir.exists(plotDir)){dir.create(plotDir,recursive=T)}
 png(paste0(plotDir,"/",fname,".png"),width = 15, height = 4,units = 'in',res= 300)
  
  y_value=celltmp$seg.mean.LOWESS
  cpg=celltmp$seg_matrix
  seg_end=which(cpg==1)
  seg_value <- rep(NA,nrow(celltmp))
  seg_start <- 1
  for(k in seg_end){
    seg_k <- mean(y_value[seg_start:k])
    seg_value[seg_start:k] <- seg_k
    seg_start <- k+1
  }
  
  maploc=1:nrow(celltmp)
  chromnum=sapply(strsplit(rownames(celltmp),'-|_|:'),'[',1)
  chr_name=as.character(unique(chromnum))
  if(grepl('chr',chromnum[1])){
    chromnum=sapply(strsplit(chromnum,'hr'),'[',2) ## if chrommnum with chr
  }else{
    chromnum=chromnum}
  chromnum=chromnum[order(as.integer(gsub("[^0-9]", "", chromnum)))] #order by chr
  chromnum_xaxis = table(chromnum)[order(as.integer(gsub("[^0-9]", "", names(table(chromnum)))))]#order by chr
  chrline=cumsum(chromnum_xaxis)
  

  ###
  plot(x=maploc, y=y_value, ylab="Ratio",
       main=paste("Segmentation"),
       axes = F,xlab="chrom", pch=20, cex=0.5, ylim=c(0,3),col=rgb(0.5,0.5,1,alpha = 0.5),)
  points(x=maploc, y=seg_value, type='l',col = rgb(0.7,0,0.2,alpha = 0.7), lwd=1)
  
  
  if(length(unique(chromnum))==1){
    axis(side=1, at=chromnum_xaxis*0.5, labels = sapply(strsplit(chr_name,'hr'),'[',2), cex.axis=1)
  }else{
    axis(side=1, at=c(0,chrline[1:(length(chrline))-1])+chromnum_xaxis*0.5, labels = sapply(strsplit(chr_name,'hr'),'[',2), cex.axis=1)}
  abline(v=chrline, col='blue', lty=1, lwd=1)
  abline(v=1, col='blue', lty=1, lwd=1)
  axis(side=2)
  

}


#' @title get_group_color_palette()
#' @export
get_group_color_palette <- function (col_nm="Set3") {
  return(colorRampPalette(RColorBrewer::brewer.pal(12, col_nm)))}

#' @title add_ends()
#' @export
add_ends <- function(dat,bin_size=50000){
  dat<-dat[order(dat[,1],dat[,2]),]
  for(j in 1:(nrow(dat)-1)){
    if (dat$chrom[j]==dat$chrom[j+1]){
      dat$end[j]=dat$chrompos[j+1]-1
    }else{
      dat$end[j]<-dat$chrompos[j]+bin_size
    }
  }
  dat$end[nrow(dat)]<-dat$chrompos[nrow(dat)]+bin_size
  return(dat)
}
#heatmap plot for matrix 
#' @title HeatmapPlot()
#' @param dat  ratio or CNV matrix with columns are cells and rows are features,first three columns are c("chr","start","end"), other columns are cells
#' @param row_anno data.frame with group annotation in columns and cellID as rownames
#' @param ref_group_names group name in the first column of row_anno
#' @param custom_colors list of paramters for color values
#' @param label_genes vector of genes
#' @param annotations GRanges object
#' @export
HeatmapPlot<-function(dat,plotDir,type="any",fname=NULL,
                      clust=FALSE,
                      clustering_method_rows = "complete",
                      row_anno=NULL,
                      width = 20, height = 9,
                      show_row_names=FALSE,
                      left_anno_col=NULL,
                      color_bars= "discrete",#"continuous"
                      n_breaks = 6,
                      column.title = "Genomic Region",
                      legend.direction = "vertical",  #"horizontal"
                      legend_side = "right",
                      legend_titles = NULL,
                      show_legend_row = TRUE,
                      draw_normal = FALSE,
                      ref_group_names=NULL,
                      max.legend.value = NULL,
                      custom_colors=NULL,
                      device="png",
                      label_genes=NULL,
                      annotations=NULL,
                      label_color="black"
                      ){
  nd <- list("ComplexHeatmap","circlize","ggplot2","ggsci","futile.logger")
  #Plus.library(nd)
  lapply(nd, function(pkg) invisible(require(pkg, character.only = TRUE)))
  
  if(!dir.exists(plotDir)){dir.create(plotDir,recursive=T)}
  
  dat <- as.data.frame(dat)
  chr_idx <- which(tolower(colnames(dat))%in% tolower(c("chr","chrom")))
  if(length(chr_idx)!=0)colnames(dat)[chr_idx] <- "chrom"
  start_idx <- which(tolower(colnames(dat))%in% tolower(c("loc.start","start","chrompos")))
  if(length(start_idx)!=0)colnames(dat)[start_idx] <- "chrompos"
  end_idx <- which(tolower(colnames(dat)) %in% tolower(c("loc.end","end","stop")))
  if(length(end_idx)!=0){
    colnames(dat)[end_idx] <- "end"
  }else{
    dat<-add_ends(dat)
  }
  if(!is.numeric(dat$chrompos)){
    dat$chrompos <- as.numeric(as.character(dat$chrompos))
  }
  if(!is.numeric(dat$end)){
    dat$end <- as.numeric(as.character(dat$end))
  }
  
  #sapply(strsplit(cellfiles[,1],'\\.'),'[',1)
  if(any(grepl("chr",dat$chrom))){dat$chrom=sapply(strsplit(as.character(dat$chrom),'chr'),'[',2)}
  dat<-dat[dat$chrom %in%c(1:23,"X",paste0("chr",c(seq(1,22),"X"))),]
  get_group_color_palette <- function (col_nm="Set3") {
    return(colorRampPalette(RColorBrewer::brewer.pal(12, col_nm)))}
  
  color <- rep(c("grey90", "grey50"), length.out = length(unique(dat$chrom)))
  text_colors <- rep(c("black", "white"), length.out = length(unique(dat$chrom)))
  
  dat$chrom[dat$chrom %in% c(23,"X","chrX")] <- "X"
  dat$chrom <- factor(dat$chrom,levels=c(1:22,"X"))
  dat<- dat[order(dat$chrom,dat$chrompos),]
  chr_cluster<-droplevels(dat$chrom)
  #chr_cluster[chr_cluster %in% 23] <-"X"
  #levels(chr_cluster)<- 1:23
  top_color <- HeatmapAnnotation(
    cluster = anno_block(gp = gpar(fill = color,col = "NA"),
                         labels = levels(droplevels(chr_cluster)), 
                         #labels_gp = gpar(cex = 0.9, col = "black"),
                         labels_gp = gpar(col = text_colors),
                         height = unit(0.5, "cm")
                         )
    )
  if(!is.null(row_anno)){
    if(draw_normal){
      if(!is.null(ref_group_names) ){
        cells_ref <- rownames(row_anno)[row_anno[,1]%in% ref_group_names]
      }else{
        cells_ref <- c()
        flog.warn("ref_group_names is NULL, no reference cells were plotted separately.")
      }
      dat_ref <- dat[,cells_ref]
      
    }
    if(is.null(left_anno_col)){
      col_map_func <- function(dat){
        res = list()
        for(n in colnames(dat)){
          tmp_name = sort(unique(dat[, n]))
          tmp_col = pal_ucscgb(alpha = 0.8)(length(tmp_name))
          names(tmp_col) = tmp_name
          res[[n]] = tmp_col
        }
        return(res)
      }
    left_anno_col <-  col_map_func(row_anno)
  }
  
    row_color <- rowAnnotation(df=row_anno,
                              col=left_anno_col,
                              gp = gpar(col = "NA"),
                              annotation_legend_param = list(nrow = 30),
                              show_legend = show_legend_row
    )
  }else{row_color <- NULL}
  
  dat$chrom<-NULL
  dat$chrompos<-NULL
  dat$abspos<-NULL
  dat$start <- NULL
  dat$end <- NULL
  if(is.null(max.legend.value)){
    mt <- as.matrix(t(dat))
    mt[is.null(mt)] <- NA
    max.legend.value <- ceiling(quantile(mt,0.98,na.rm=TRUE))
  }
  
  if(is.null(custom_colors)){
    if(type%in% "ratio"){
      if(max.legend.value>=4){max.legend.value=4}

        colorss <- colorRampPalette(c("#3d8bff","#ccf7f4", "white","#F5F0D4","#ffba08","#FF6600","#d00000","#9d0208","#6a040f"))(9)
       

      colors <-colorRamp2(c(0,0.5,1,1.5,2,2.5,3,3.5,max.legend.value),colorss)
      if(max.legend.value<4 &max.legend.value>2){

        colorss <- colorss[1:(2+(2*max.legend.value-1))]
        colors <-colorRamp2(c(0,0.5,seq(1,max.legend.value,by=0.5)),colorss)
      }else if(max.legend.value<=2 &max.legend.value>1){
        x.center <- 1
        quantiles = quantile(mt[mt != x.center], c(0.01,0.25,0.75, 0.99),na.rm =TRUE)
        # determine max distance from the center.
        delta = max( abs( c(x.center - quantiles[1],  quantiles[4] - x.center) ) )
        
        at_brk <- c(quantiles[1],quantiles[2],x.center,quantiles[3],quantiles[4])
        colorss <- colorRampPalette(colors = c("#3d8bff","#ccf7f4","white","#F5F0D4","#FF6600"))(length(at_brk))
        colors <- circlize::colorRamp2(at_brk,colorss)
      }else if(max.legend.value==1){
        colorss <- colorss[1:3]
        colors <-colorRamp2(c(0,0.5,1),colorss)
      }
      at_brk <- c(0:max.legend.value)
      label_brk <- c(as.character(c(0:(max.legend.value-1))),paste0(max.legend.value,"+"))
      color_bars <- "continuous"
      titles<-"Ratio"
    }else if(type%in% "CNV"){
      CN_max <- quantile(as.numeric(as.matrix(dat[!is.null(dat)])),1,na.rm=T)
      colorss <- colorRampPalette(c("#3d8bff","#C4D8F5", "grey95","#f9dcc4","#f79d65","#f27059","#85182a"))(7)
      if(CN_max>7){
        colorss_plus <- colorRampPalette(c("#85182a","black"))(3)
        colorss <- c(colorss,colorss_plus[2:(length(colorss_plus)-1)])
        colors <-colorRamp2(c(seq(1,7),CN_max),colorss)
        at_brk <- c(1:7,CN_max)
        label_brk <- c(as.character(c(1:7)),paste0(CN_max))
      }else{
        colors <-colorRamp2(c(1,2,3,4,5,6,7),colorss)
        at_brk <- c(1:7)
        label_brk <- c(as.character(c(1:6)),"7+")
      }

      #colorss <- adjustcolor(colorss, alpha.f = 0.8)
      CN_mean <- quantile(as.numeric(as.matrix(dat[!is.null(dat)])),0.5,na.rm=T)
      if(CN_mean<=2){
        colorss <- colorss[-1]
        colors <-colorRamp2(c(1,2,3,4,5,6),colorss)
        at_brk <- c(1:6)
        label_brk <- c(as.character(c(1:5)),"6+")
      }

      
      color_bars <- "discrete"
      titles<-"IntegerCNV"
    }else if(type %in% "count"){
      dat_p <- as.matrix(dat[,!grepl("chrom|chrompos|abspos|start|end",colnames(dat))])
      if(max.legend.value>=5){max.legend.value=5}
      #at_brk=unique(round(sort(seq(1,max.value,length=n_breaks))));
      #colorss <- c(colorRampPalette(colors = c("#FFDC00","#FF6600","#d00000","#9d0208","#6a040f"))(5)) #"#FFE800",#FFE5E0"
      colorss <- c(colorRampPalette(colors = c("#FF6600","#d00000","#9d0208","#6a040f"))(5)) #"#FFE800",#FFE5E0"
      colorss <- c("white",colorss) #5E95BE
      
      if(max.legend.value<5 & max.legend.value>1){
        colorss <- colorss[1:(max.legend.value+1)]
      }else if(max.legend.value==1){
        colorss <- colorss[1:2]
      }
      at_brk <- c(0:max.legend.value)
      colors <- colorRamp2(at_brk,colorss)
      
      label_brk <- c(as.character(c(0:(max.legend.value-1))),paste0(max.legend.value,"+"))
      color_bars <- "discrete"
      titles<-"counts"
    }else{
      nb_breaks <- n_breaks
      dat_p <- as.matrix(dat[,!grepl("chrom|chrompos|abspos|start|end",colnames(dat))])
      at_brk=round(unique(sort(seq(min(dat_p),max(dat_p),length=nb_breaks))),2);
      colors <- c(colorRampPalette(colors = c("#FFDC00","#FF6600","red4"))(length(at_brk))) #"#FFE800",#FFE5E0"
      colors <- c("white",colors) #5E95BE
      at_brk=sort(unique(c(0,at_brk)))
      colors <- colors[1:length(at_brk)]
      label_brk <- c(as.character(at_brk[1:(length(at_brk)-1)]),paste0(at_brk[length(at_brk)],"+"))
      
      color_bars <- color_bars
      titles<-"value"
    }
  }else{
    colors <- custom_colors$colors
    color_bars <- custom_colors$color_bars
    titles<- custom_colors$titles
    at_brk<-custom_colors$at_brk
    label_brk <- custom_colors$label_brk
  }
  if(!is.null(legend_titles)){titles <- legend_titles}

  if(!is.null(label_genes)){
    if(is.null(annotations)){
      require(Signac)
      require(EnsDb.Hsapiens.v86)
      annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
      seqlevelsStyle(annotations) <- "UCSC"
    }
    peaks <- data.frame(Region=rownames(dat))
    peaks <- peaks %>%
    tidyr::separate(
      col = Region,
      into = c("chr", "start", "end"),
      sep = ":|-|_",   
      remove = FALSE
    ) %>%
    dplyr::mutate(
      start = as.integer(start),
      end = as.integer(end)
    ) %>%
    dplyr::select(chr, start, end, everything())
    peak_ranges <- GRanges(
      seqnames = peaks$chr,
      ranges = IRanges(
        start = peaks$start,
        end = peaks$end
      )
    )
    gene_miss <- label_genes[!label_genes %in% annotations$gene_name]
    if(length(gene_miss)>0){
      message(paste0(paste(gene_miss,collapse=",")," not found in annotations."))
    }
    Genes_gr <- annotations[annotations$gene_name %in% label_genes,,drop=F]

    hits <- findOverlaps(Genes_gr, peak_ranges)
    gene_annotation <- data.frame(
      gene = Genes_gr$gene_name[queryHits(hits)],
      peak_index = subjectHits(hits)  # 对应矩阵的列索引
    )%>% 
      distinct(gene, .keep_all = TRUE)

    gene_annotation_res <- gene_annotation %>%
    group_by(peak_index) %>%
    summarise(gene = paste(gene, collapse = ", "))



    bottom_col <-  HeatmapAnnotation(
      gene = anno_mark(
        at = gene_annotation_res$peak_index, 
        labels = gene_annotation_res$gene,
        side = "bottom",
        labels_gp = gpar(fontsize = 10,
                        col = label_color
                        ))
      )

  }else{
    bottom_col <- NULL
  }

  ht_plot = Heatmap(as.matrix(t(dat)),
                    cluster_rows = clust,
                    clustering_method_rows = clustering_method_rows,
                    cluster_columns = F,
                    show_row_names = show_row_names,
                    row_title = NULL,  ##
                    column_split = chr_cluster,
                    heatmap_legend_param = list(
                      direction = legend.direction,
                      at = at_brk,
                      labels = label_brk,
                      title = titles,
                      title_position = "leftcenter-rot", # 图例标题位置
                      legend_height = unit(3, "cm"), #图例长度
                      color_bar = color_bars),
                    top_annotation = top_color,
                    bottom_annotation=bottom_col,
                    left_annotation = row_color,
                    column_title = column.title,
                    show_column_names = FALSE,
                    column_title_side = c("bottom"),
                    col=colors)
  
  if(draw_normal){
    if(ncol(dat_ref)>0.5*ncol(dat) &ncol(dat_ref)>3000){
      n_cells <- min(round(0.25*ncol(dat)),3000)
      sample_col <- sample(1:ncol(dat_ref),n_cells)
      dat_ref <- dat_ref[,sample_col]
    }
    ht_normal = Heatmap(t(dat_ref), 
                        col=colors,
                        cluster_rows = F,cluster_columns = F,
                        show_column_names = F,show_row_names = F,
                        column_split = chr_cluster,
                        show_heatmap_legend=FALSE,
                        row_title = "Ref.(Cells)",
                        row_title_side = c("right"),
                        row_title_rot = 90,
                        column_title = NULL)
  }
  # annotation_to_remove <- colnames(row_anno)[sapply(row_anno,function(x) length(unique(x)))>50]
  # annotation_index <- which(names(ht_plot@left_annotation) == annotation_to_remove)
  # 
  if(!is.null(fname)){
    if(device=="png"){
      png(paste0(plotDir,"/",fname,".png"),width = width,height = height,units = 'in',res= 300)
    }else{
      pdf(paste0(plotDir,"/",fname,".pdf"),width = width,height = height)
    }
    if(draw_normal){
      ComplexHeatmap::draw(ht_normal%v%ht_plot, padding = unit(c(20, 10, 10, 10), "mm"),heatmap_legend_side = legend_side)
      dev.off()
    }else{
      ComplexHeatmap::draw(ht_plot, heatmap_legend_side = legend_side,annotation_legend_side = legend_side) # 图例位置
      dev.off()
    }

  }
  return(ht_plot)
}

#' @title abspos()
#' @export
abspos <- function(df){
  suppressMessages(library(bit64))
  coln <- colnames(df)  
  if(any(grepl('chr|chrom|Chromosome',coln,ignore.case = T))){
    coln[grepl('chr|chrom|Chromosome',coln,ignore.case = T)] <- 'chrom'
  }
  if(any(grepl('start',coln,ignore.case = T))){
    coln[grepl('start',coln,ignore.case = T)] <- 'start'
  }
  if(any(grepl('end',coln,ignore.case = T))){
    coln[grepl('end',coln,ignore.case = T)] <- 'end'
  }
  colnames(df) <- coln
  
  chromoPos=data.frame(chr=df$chrom,start=df$start,end=df$end)
  chromosome=unique(df$chrom)
  endlength=as.numeric(sapply(chromosome, function(i,chromoPos){
    subdata=max(as.numeric(chromoPos$end[chromoPos$chr==i]))
    return(subdata)
  },chromoPos))
  
  endSum=c()
  endSum[1]=endlength[1]
  if(length(endlength)>1){
    for(i in 2:length(endlength)){
      endSum[i]=sum(endlength[1:i])
    }
  }
  
  chromoPos$abspos=as.numeric(chromoPos$start)
  if(length(chromosome)>1){
    for (i in 2:length(chromosome)){
      chromoPos$abspos[chromoPos$chr==chromosome[i]]=as.numeric(chromoPos$start[chromoPos$chr==chromosome[i]])+endSum[i-1]
    }
  }
  
  df$abspos=as.numeric(df$start)
  df$absend=as.numeric(df$end)
  if(length(chromosome)>1){
    for (i in 2:length(chromosome)){
      df$abspos[df$chr==chromosome[i]]=as.numeric(df$start[df$chr==chromosome[i]])+endSum[i-1]
      df$absend[df$chr==chromosome[i]]=as.numeric(df$end[df$chr==chromosome[i]])+endSum[i-1]
    }
  }
  return(df)
}

#' @title plot_peak()
#' @export
plot_peak<-function(y0,aaa,main0="",breaks0=30,Zcut=-Inf)
{
  y<-y0
  if(sum(y<Zcut)>0)
  {
    zz<-min(y[which(y>Zcut)])-2
    y[which(y<=Zcut)]<-zz
  }
  mm<-min(y)
  MM<-max(y)
  diff_c<-MM-mm
  MM<-MM+diff_c*0.2
  mm<-mm-diff_c*0.2
  x<-seq(mm,MM,by=(MM-mm)/1000)
  h<-hist(y,breaks=breaks0,plot=F)
  plot(h,col="lightblue",main=main0)
  n<-length(y)*(h$breaks[2]-h$breaks[1])
  for(i in 1:nrow(aaa))
  {
    z<-dnorm(x,aaa[i,2],aaa[i,3])*aaa[i,1]*n
    abline(v=aaa[i,2],col=c(i+1),lwd=2)
    points(z~x,type="l",col=c(i+1),lwd=2)
  }
}


#' @title plot_combine_seg()
#' @export 
plot_combine_seg <- function(clonal.res,ylim=NULL,outplot_name=NULL,show_dots=FALSE,outdir="./",ggarrange_width=c(2.5,1)){
  if(is.null(outplot_name)){
    outplot_name <- "clonalCN_seg.pdf"
  }
  plt_clones <- lapply(as.character(sort(as.numeric(names(clonal.res)))),function(x,clonal.res,outdir){
    outplot_name_i <- x
    ratiores <- clonal.res[[outplot_name_i]]$input_BinSegRatio
    integerCNVres <- clonal.res[[outplot_name_i]]$seg.dat
    score <-  round(clonal.res[[outplot_name_i]]$score$score,2)
    plt <- cellLineSeg_v3(binRatio.df=ratiores,integerCNV.df=integerCNVres,outdir=outdir,
                          ylim = ylim,
                          fname=paste0("C",outplot_name_i,": score=",score),
                          plot_dots = show_dots,
                          ggarrange_width=ggarrange_width,
                          height = 2,width=10,outPlot=FALSE)
    return(plt$ggarranged_p)
  },clonal.res,outdir)
  pcom <- cowplot::plot_grid(plotlist=plt_clones,ncol = 1,align="v")
  ggsave(paste0(outdir, "/",outplot_name),pcom, width=14, height=2*length(plt_clones),device = pdf,bg="white")
}

#' @title replotInfercnv()
#' @param expr.file data.frame or matrix of tumor cell cnv (row:gene, column:cell)
#' @param gene.order.file path of gene.order.file
#' @param metadata metsdata od cell annotation (cell order)
#' @param ref.mat data.frame or matrix of normal cell cnv (if draw_normal=TRUE)
#' @param names output name
#' @export 
replotInfercnv <- function(expr.mat,
                           gene.order.file,
                           metadata=NULL,
                           left_anno_by=NULL,
                           row_split_by=NULL,
                           names="Heatmap",
                           left_anno_col=NULL,
                           clust_rows=F,
                           draw_normal=FALSE,
                           ref.mat=NULL,
                           width=10,height=5,
                           show_legend_row = TRUE,
                           max.legend.value = NULL,
                           clustering_method_rows = "complete",
                           show_row_names=FALSE,
                           direction = "vertical",
                           titles = "inferCNV score",
                           column.title = "Genomic Region",
                           legend_side = "right",
                           color.breaks=NULL,
                           plotDir="./",
                           device="png",
                           type="ratio",
                           color_bars= "continuous" #"discrete"
                           ){
  nd <- list("data.table","ComplexHeatmap","dplyr","circlize","ggplot2","ggsci")
  lapply(nd, require, character.only = TRUE)
  if(!dir.exists(plotDir)){dir.create(plotDir,recursive=T)}
  gene_order <- fread(gene.order.file)
  expr.mat <- expr.mat %>% as.data.frame()
  inter_genes <- intersect(rownames(expr.mat),gene_order$V1)
  gene_order <- gene_order[gene_order$V1 %in% inter_genes,]
  expr.mat <- expr.mat[gene_order$V1,]
  pos <- grep(paste0("chr",1:22,collapse = "|"),gene_order$V2)
  expr.mat <- expr.mat[pos,]
  gene_order <- gene_order[pos,]
  
  # top annotation
  color <- rep(c("grey90", "grey50"), length.out =22)
  text_colors <- rep(c("black", "white"), length.out = 22)
  top_color <- HeatmapAnnotation(
    cluster = anno_block(gp = gpar(fill = color,col = "NA"),
                         labels = 1:22,
                         labels_gp = gpar(col = text_colors),
                         height = unit(0.5, "cm")
  ))
  #expr.mat order by metadata
  if(!is.null(metadata)){
    expr.mat <- expr.mat[,match(rownames(metadata),colnames(expr.mat))]
    if(is.null(left_anno_by)){left_anno_by=colnames(metadata)}
    left_anno <- metadata[,left_anno_by,drop=F]
    
    if(is.null(left_anno_col)){
      col_map_func <- function(dat){
        res = list()
        for(n in colnames(dat)){
          tmp_name = sort(unique(dat[, n]))
          tmp_col = pal_igv("default",alpha = 0.8)(length(tmp_name))
          names(tmp_col) = tmp_name
          res[[n]] = tmp_col
        }
        return(res)
      }
      left_anno_col <-  col_map_func(left_anno)
    }
    
    row_color <- rowAnnotation(df=metadata,
                               col=left_anno_col,
                               gp = gpar(col = "NA"),
                               annotation_legend_param = list(nrow = 30),
                               show_legend = show_legend_row
    )
    
    if(!is.null(row_split_by)){row_split_by=metadata[,row_split_by]}
  }else{
    row_color=NULL
  }
  
  if(type%in% "ratio"){
    if(is.null(color.breaks)){
      x.center <- 1
      quantiles = quantile(expr.mat[expr.mat != x.center], c(0.01, 0.99),na.rm =TRUE)
      # determine max distance from the center.
      delta = max( abs( c(x.center - quantiles[1],  quantiles[2] - x.center) ) )
      low_threshold = x.center - delta
      high_threshold = x.center + delta
      # bk <-seq(0.5,2,length.out=16)
      bk <- c(low_threshold,x.center,high_threshold)
    }else{
      bk <- color.breaks
    }
    cols <- colorRampPalette(colors = c("#2A72B2", "grey95", "#6a040f"))(length(bk))
    colors <- circlize::colorRamp2(bk,cols)
    label_brk <- as.character(round(bk,2))
    color_bars <- "continuous"
  }else if(type%in% "CNV"){
    CN_max <- quantile(as.numeric(as.matrix(expr.mat[!is.null(expr.mat)])),1,na.rm=T)
      # colorss <-  c("#1A4E7A","#2A72B2", "#F2F2F2","#B97B77", "#6a040f","#3F0209")
      colorss <-  c("#1a4fb3", "#6699FF", "grey95","#F6DAC3","#EE9863","#EB6363")
      colors <-colorRamp2(c(0,1,2,3,4,5),colorss)
      bk <- c(0,1,2,3,4,5)
      label_brk <- c(as.character(c(0,1,2,3,4,5)))
      
      color_bars <- "discrete"
      titles <- "CNV"
  }
  
  ht_plot <- Heatmap(as.matrix(t(expr.mat)),
                     cluster_rows = clust_rows,
                     clustering_method_rows = clustering_method_rows,
                     cluster_columns = F,
                     show_row_names = show_row_names,
                     row_split = row_split_by,
                     row_title = NULL,
                     column_split = factor(gene_order$V2, paste("chr",1:22,sep = "")),
                     heatmap_legend_param = list(
                       direction = direction,
                       at = bk,
                       labels = label_brk,
                       title = titles,
                       title_position = "leftcenter-rot", # 
                       legend_height = unit(3, "cm"), #
                       color_bar = color_bars),
                     top_annotation = top_color,
                     left_annotation = row_color,
                     column_title = column.title,
                     show_column_names = FALSE,
                     column_title_side = c("bottom"),
                     col=colors)
  if(draw_normal){
    if(ncol(ref.mat)>0.5*ncol(dat) &ncol(ref.mat)>3000){
      n_cells <- min(round(0.25*ncol(dat)),3000)
      sample_col <- sample(1:ncol(ref.mat),n_cells)
      ref.mat <- ref.mat[,sample_col]
    }
    ref.mat <- ref.mat %>% as.data.frame()
    ref.mat <- ref.mat[gene_order$V1,]
    ht_normal = Heatmap(t(ref.mat), 
                        col=colors,
                        cluster_rows = F,cluster_columns = F,
                        show_column_names = F,show_row_names = F,
                        column_split = factor(gene_order$V2, paste("chr",1:22,sep = "")),
                        show_heatmap_legend=FALSE,
                        row_title = "Ref.(Cells)",
                        row_title_side = c("right"),
                        row_title_rot = 90,
                        column_title = NULL)
    
  }
  if(!is.null(names)){
    if(device=="png"){
      png(paste0(plotDir,"/",names,".png"),width = width,height = height,units = 'in',res= 300)
    }else{
      pdf(paste0(plotDir,"/",names,".pdf"),width = width,height = height)
    }

    if(draw_normal){
      ComplexHeatmap::draw(ht_normal%v%ht_plot, padding = unit(c(20, 10, 10, 10), "mm"),heatmap_legend_side = legend_side)
      dev.off()
    }else{
      ComplexHeatmap::draw(ht_plot, heatmap_legend_side = legend_side,annotation_legend_side = legend_side) # 图例位置
      dev.off()
    }    
  }
  return(ht_plot)
}




#' @title plot_genome()
#' @description plot summed count of CNAs (amplifications or deletions)
#' @export 
plot_genome <- function(
    chr_bg_color=NULL,
    show_x_axis=TRUE,
    special_region=NULL,
    linewidth=0.5,
    vline=T,
    exclude=c("chrX","chrY","chrM")){
  chrinfo = DNAcopy::cytoBand
  chrinfo$chromNum <- gsub("chr0","chr",chrinfo$chromNum)
  chrinfo$chr_len = chrinfo$chromEnd - chrinfo$chromStart
  chrinfo <- chrinfo %>%
    filter(!chromNum %in% exclude ) 
  
  chrinfo$absEnd = cumsum(as.numeric(chrinfo$chr_len))
  chrinfo$absStart = chrinfo$absEnd-chrinfo$chr_len
  tmp_fun = function(x,y){
    (min(x)+max(y))/2
  }
  chr_center_pos = chrinfo %>%
    group_by(chromNum) %>%
    summarise(center = tmp_fun(absStart, absEnd))
  
  chrinfo_other = list(
    chr_center_pos = as.numeric(chr_center_pos$center), # 
    chr_name = as.vector(chr_center_pos$chromNum), # 
    chr_band = chrinfo[chrinfo$chromStart==0, 'absStart'] # 
  )
  chrinfo_other$chr_name <- gsub("chr","",chrinfo_other$chr_name)
  chrinfo_other$chr_name[as.character(chrinfo_other$chr_name)%in%c("19","21")] <- ""
    
  chrinfo_chr_absstart = chrinfo[, c('chromNum', 'absStart', 'absEnd')] %>%
    group_by(chromNum) %>%
    summarise(absStart=min(absStart), absEnd=max(absEnd)) %>% as.data.frame()
  rownames(chrinfo_chr_absstart) = chrinfo_chr_absstart$chromNum
  
  chrstart_min = 0 # 
  chrstart_max = max(chrinfo$absEnd) # 
  
  ##### 
  plot_panel = ggplot()+
    theme_bw()+
    theme(panel.grid = element_blank())
  
  ##### 
  plot_chr=plot_panel
  if(vline){
    plot_chr = plot_chr + geom_vline(xintercept = chrinfo_other$chr_band, color='black', linewidth=0.5, linetype='dashed')
  }
  plot_chr = plot_chr +
    scale_x_continuous(breaks=chr_center_pos$center, labels = chrinfo_other$chr_name,
                       limits=c(chrstart_min, chrstart_max),
                       expand = expansion(mult = c(0,0)))+
    labs(x='chromosome')
  
  if(!show_x_axis){
    #####
    plot_chr = plot_chr+
      theme(axis.text.x = element_blank(),
            axis.ticks.x  = element_blank(),
            axis.title.x = element_blank())
    
    
  }
  #####
  if(is.null(chr_bg_color)){
    bg_color = rep(c('gray', 'white'), 12)
    for(i in 1:nrow(chrinfo_chr_absstart)){
      tmp_grop = rasterGrob(bg_color[i],
                            width=unit(1,"npc"), height = unit(1,"npc"),
                            interpolate = TRUE)
      plot_chr = plot_chr + annotation_custom(tmp_grop,
                                              xmin=chrinfo_chr_absstart$absStart[i],
                                              xmax=chrinfo_chr_absstart$absEnd[i],ymin=-Inf,ymax=Inf)
    }
  }else if(is.character(chr_bg_color)){
    bg_color = rep(chr_bg_color, 24)
    for(i in 1:nrow(chrinfo_chr_absstart)){
      tmp_grop = rasterGrob(bg_color[i],
                            width=unit(1,"npc"), height = unit(1,"npc"),
                            interpolate = TRUE)
      plot_chr = plot_chr + annotation_custom(tmp_grop,
                                              xmin=chrinfo_chr_absstart$absStart[i],
                                              xmax=chrinfo_chr_absstart$absEnd[i],ymin=-Inf,ymax=Inf)
    }
  }else if(is.list(chr_bg_color)){
    bk = chr_bg_color$bk
    col = chr_bg_color$col
    alpha = chr_bg_color$alpha
    col_num = chr_bg_color$value
    col_fun = circlize::colorRamp2(bk, col,transparency = alpha)
    
    for(i in 1:nrow(chrinfo_chr_absstart)){
      tmp_grop = rasterGrob(col_fun(col_num[i]),
                            width=unit(1,"npc"), height = unit(1,"npc"),
                            interpolate = TRUE)
      plot_chr = plot_chr + annotation_custom(tmp_grop,
                                              xmin=chrinfo_chr_absstart$absStart[i],
                                              xmax=chrinfo_chr_absstart$absEnd[i],ymin=-Inf,ymax=Inf)
    }
  }else{
    bg_color = chr_bg_color
    for(i in 1:nrow(chrinfo_chr_absstart)){
      tmp_grop = rasterGrob(bg_color[i],
                            width=unit(1,"npc"), height = unit(1,"npc"),
                            interpolate = TRUE)
      plot_chr = plot_chr + annotation_custom(tmp_grop,
                                              xmin=chrinfo_chr_absstart$absStart[i],
                                              xmax=chrinfo_chr_absstart$absEnd[i],ymin=-Inf,ymax=Inf)
    }
    
  }
  ##### (chr:start-end)
  if(!is.null(special_region)){
    tmp_chr = strsplit(special_region, ':')[[1]][1]
    tmp_start = strsplit(special_region,':|-')[[1]][2]
    tmp_end = strsplit(special_region,':|-')[[1]][3]
    if(tmp_start=='none'){
      tmp_start=0
    }else{
      tmp_start = as.numeric(tmp_start)
    }
    if(tmp_end=='none'){
      tmp_end = chrinfo_chr_absstart[tmp_chr, 'absEnd'] - chrinfo_chr_absstart[tmp_chr, 'absStart']
    }else{
      tmp_end = as.numeric(tmp_end)
    }
    chr_start = chrinfo_chr_absstart[tmp_chr, 'absStart'] + tmp_start
    chr_end = chrinfo_chr_absstart[tmp_chr, 'absStart'] + tmp_end
    #print(c(chr_start,chr_end))
    plot_chr = plot_chr+scale_x_continuous(limits = c(chr_start, chr_end))
  }
  return(plot_chr)
}

plot_theme <- theme(plot.background=element_blank(),
    panel.background=element_blank(),
    panel.grid = element_blank(), 
    axis.title.x = element_blank(), 
    axis.title.y = element_text(size = 7),
    axis.ticks=element_line(color="black",linewidth=0.5),
    axis.text=element_text(color="black",size=7),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    legend.background=element_blank(),
    legend.key=element_blank(),
    legend.text=element_text(color="black",size=7),
    legend.title=element_blank())

#plots
#plot for count matrix 
#' @title plot4Count()
#' @description Data visualization for the # of reads(or peaks) per cell, and # of reads (or cells) per peak (column : cell)
#' @export 
plot4Count <- function(matrix,outdir="./",outfile="rawdataCheck.pdf",
                       addline=FALSE,
                       format = "pdf",
                       min_Nreads_per_cell=NULL,
                       min_Nfeature_per_cell=NULL,
                       min_Nreads_per_feature=NULL,
                       min_Ncell_per_feature=NULL){
  if(!is.matrix(matrix)){matrix <- as.matrix(matrix)}
  suppressMessages(library(pheatmap))
  suppressMessages(library(RColorBrewer))
  plotDir <- outdir
  if(!file.exists(plotDir)){dir.create(plotDir,recursive=T)}
  reads_per_cell = colSums(matrix)
  reads_per_peak = rowSums(matrix)
  peaks_per_cell = colSums(matrix>0)
  cells_per_peak = rowSums(matrix>0)
  if(format=="pdf"){
    pdf(file=paste0(plotDir,"/",outfile),width = 8,height = 8)
  }else{
    png(file=paste0(plotDir,"/",outfile),width = 8,height =8,units = 'in',res= 300)
  }
  par(mfrow=c(3,2))
  #layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  hist(log10(reads_per_cell),main='reads per cell',col='wheat',breaks = 20,prob = TRUE)
  lines(density(log10(reads_per_cell)), # density plot
        lwd = 2, # thickness of line
        col = "chocolate3")
  abline(v = log10(mean(reads_per_cell)), col = "blue",lty=2)
  text(x=log10(mean(reads_per_cell))+0.2,y=max(density(log10(reads_per_cell))$y-0.1),paste0("Mean: ",round(mean(reads_per_cell))),cex = .75,col="blue")
  
  if(addline){
    abline(v = log10(min_Nreads_per_cell), col = "red",lty=2)
    }
  hist(log10(peaks_per_cell), main='peaks per cell', col='wheat',breaks = 20,prob = TRUE)
  lines(density(log10(peaks_per_cell)), # density plot
        lwd = 2, # thickness of line
        col = "chocolate3")
  abline(v = log10(mean(peaks_per_cell)), col = "blue",lty=2)
  text(x=log10(mean(peaks_per_cell))+0.2,y=max(density(log10(peaks_per_cell))$y-0.1),paste0("Mean: ",round(mean(peaks_per_cell))),cex = .75,col="blue")
  if(addline){abline(v = log10(min_Nfeature_per_cell), col = "red",lty=2)}
  
  hist(log10(reads_per_peak),main='reads per peak',col='wheat',breaks=20,prob = TRUE)
  lines(density(log10(reads_per_peak)), # density plot
        lwd = 2, # thickness of line
        col = "chocolate3")
  abline(v = log10(mean(reads_per_peak)), col = "blue",lty=2)
  text(x=log10(mean(reads_per_peak)) + 0.2, y=max(density(log10(reads_per_peak))$y*0.6),paste0("Mean: ",round(mean(reads_per_peak))),cex = .75,col="blue")
  
  if(addline){abline(v = log10(min_Nreads_per_feature), col = "red",lty=2)}
  hist(log10(cells_per_peak),main='cells per peak',col='wheat',breaks=20,prob = TRUE)
  lines(density(log10(cells_per_peak)), # density plot
        lwd = 2, # thickness of line
        col = "chocolate3")
  abline(v = log10(mean(cells_per_peak)), col = "blue",lty=2)
  text(x=log10(mean(cells_per_peak)) + 0.2, y=max(density(log10(cells_per_peak))$y*0.6),paste0("Mean: ",round(mean(cells_per_peak))),cex = .75,col="blue")
  if(addline){abline(v = log10(min_Ncell_per_feature), col = "red",lty=2)}
  plot(reads_per_cell, peaks_per_cell, log='xy', col='wheat')
  
  dev.off()
}

#' @title Vlnplot4Count()
#' @export 
Vlnplot4Count <- function(matrix,project=NULL){
  require(ggplot2)
  require(cowplot)
  plotDir <- paste0("./plots")
  if(!file.exists(plotDir)){dir.create(plotDir,recursive=T)}
  if (is.null(project)) {project="Project"}
  if(!is.matrix(matrix)){matrix <- as.matrix(matrix)}
  ncells=dim(matrix)[2]
  
  plotmain="pct_cells_per_peak"
  cells_per_peak = rowSums(matrix>0)
  pct_cells_per_peak <- cells_per_peak*100/ncells
  df <- data.frame(pct_cells_per_peak=pct_cells_per_peak)  
  df$Project=project
  # p <- ggplot(df, aes(x=Project, y=pct_cells_per_peak,fill=Project)) + 
  #   geom_violin()+labs(x = "", y = "", fill = NULL)  +
  #   theme_cowplot()+ 
  #   geom_jitter(shape=16,size=0.5, position=position_jitter(0.2))+
  #   theme(legend.position="none",plot.title = element_text(hjust = 0.5))+
  #   ggtitle(plotmain)
  # ggsave(paste0(plotDir,"/",plotmain,".pdf"),p,width =8,height = 6)
  pdf(file=paste0(plotDir,"/hist_",plotmain,".pdf"),width = 8,height = 6)
  hist(pct_cells_per_peak,main=plotmain,col='wheat',breaks=50,prob = TRUE)
  lines(density(pct_cells_per_peak), # density plot
        lwd = 2, # thickness of line
        col = "chocolate3")
  dev.off()
}

#'
#' @title plot_segment()
#' @description plot segmentation across genome
#' @param main  main=paste("Segmentation: ",file_seg,", sub_",sub_cj)
#' @param y_value for point plot
#' @param seg_value for line
#' @param chromnum chromosome annotation with the same length of y_value
#' @param filename output paste0(plotDir,"/",de_noise_lb,"segment_",file_seg,"_sub_",sub_cj,".png")
#' @export 

plot_segment <- function(y_value,seg_value, ylab,main,ylim,chromnum,filename){
  png(filename,width = 10, height = 3,units = 'in',res= 300)
  
  chromnum_xaxis = table(chromnum)[order(as.integer(gsub("[^0-9]", "", names(table(chromnum)))))]#order by chr
  chrline=cumsum(chromnum_xaxis)
  chr_name=as.character(unique(chromnum))
  maploc=seq_len(y_value)
  plot(x=maploc, y=y_value, ylab="Relative coverage",
     main=main,
     ylim= range(y_value),
     xaxt="n", pch=20, cex=0.3, xlab = "" ,cex.main=1)
  points(x=maploc, y=seg_value, type='l', col='red', lwd=1)
  if(length(unique(chromnum))==1){
    axis(side=1, at=chromnum_xaxis*0.5, labels = sapply(strsplit(chr_name,'hr'),'[',2), cex.axis=1)
  }else{
    axis(side=1, at=c(0,chrline[1:(length(chrline))-1])+chromnum_xaxis*0.5, labels = sapply(strsplit(chr_name,'hr'),'[',2), cex.axis=1)}
  abline(v=chrline, col='blue', lty=1, lwd=1)
  abline(v=1, col='blue', lty=1, lwd=1)
  dev.off()
}


library(ComplexHeatmap)
library(GenomicRanges)
library(dplyr)
library(tidyr)

#' @title heatmap4peakMt()
#' @description heatmap for the 'mat' ordered by row of meta_info
#' @param mat matrix with peak in row and cells in columns, rownames pattern "chr-xxx-xxx"
#' @param sep_by a symbol for the row names separation
#' @param label_genes vector of genes
#' @param annotations GRanges object
#' @export 
heatmap4peakMt <- function(mat,meta_info=NULL,max_lim=NULL,sep_by="-",outdir="./",value.type="count",max.value = NULL,color_bars= "discrete",
                           n_breaks = 5,
                           width=10,height=6,col_list=NULL,fileout_name=NULL,clust_rows=FALSE,clustering_method_rows = "complete",
                           legend_direction = "vertical", legend_titles=NULL,heatmap_legend_side="right",
                           show_legend_row = TRUE,
                           column.title = "Genomic Region",draw_normal = FALSE,
                           ref_group_names=NULL,custom_colors=NULL,device="png",
                           label_genes=NULL,annotations=NULL,label_color="black"
                           ){
  
  if(!is.null(max_lim)){
    mat[mat>max_lim] <- max_lim
    if(all[mat]==0){
      stop("All the values are zero. Please check 'max_lim' or input 'mat'.")
    }
  }
  
  if(!is.null(meta_info)){
    if(is.data.frame(meta_info)){
      mat <- mat[,rownames(meta_info),drop=F]
    }else{
      print("'meta_info' is not a data.frame. Please check!")
    }
    
  }
  # 
  chromInfo <- data.frame(seg=rownames(mat))
  chromInfo_new <- tidyr::separate(chromInfo, seg, into = c("chrom", "start","end"), sep =sep_by)
  if(all(is.na(chromInfo_new$end))|all(is.na(chromInfo_new$start))){
    message("Please check the separate symbol of row names of 'mat', and re-set 'sep_by'.")
    stop()
  }
  chromInfo_new2 <- abspos(chromInfo_new)
  mtx_add <- data.frame(chromInfo_new2[,c("chrom","start","abspos","end")],mat)
  p <- HeatmapPlot(mtx_add,plotDir=outdir,type=value.type,max.legend.value=max.value,fname=fileout_name,color_bars= color_bars,
              n_breaks=n_breaks,
              row_anno = meta_info,width=width,height=height,left_anno_col=col_list,
              clust = clust_rows,
              clustering_method_rows =clustering_method_rows,
              column.title = column.title,
              legend.direction=legend_direction,legend_side = heatmap_legend_side,legend_titles=legend_titles,show_legend_row=show_legend_row,
              draw_normal = draw_normal,
              ref_group_names=ref_group_names,custom_colors=custom_colors,device=device,
              label_genes=label_genes,annotations=annotations,label_color=label_color)
  return(p)
}

