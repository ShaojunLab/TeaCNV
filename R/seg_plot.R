#'
#' @title seg_plot()
#' @description segment Plot (on absolute genomic position)
#' @param data.bin data.frame of bin level segment with chromosome location as rownames (rowname is important!), such as 'chr1_237513_1530360'
#' @param true_loc plot segments and dots at true genomic location. If not true, plot by order.
#' 
seg_plot = function(data.bin,name.data="",num_cells=NULL,plotDir=NULL,
                    show_specific_seg=NULL,
                    specific_seg_color="green",
                    value.bin="binRatio",
                    value.segment="SegMean",
                    ylab="Relative CN",
                    true_loc=TRUE,
                    genome="hg38",
                    plot_seg = TRUE,
                    seg_col = "#a52a2a",
                    plot_hist = TRUE,
                    plot_dots = TRUE,
                    log_ratio = FALSE,
                    color_dot =FALSE,
                    color_hist = TRUE,
                    color_hist_gradient = TRUE,
                    plot_colors = NULL,
                    color_limit = c(0,3),
                    ylim=NULL,
                    add_yline = NULL,
                    label_gene=NULL,
                    label_gene_color="blue",
                    annotations=NULL,
                    height = 3,width=14,
                    outPlot=TRUE,
                    device = "png",
                    color_seg_gradient=FALSE,
                    disconnected=TRUE
                    ){
  suppressMessages({
    library(BSgenome)
    library(dplyr)})
  g <- getBSgenome(genome, masked=FALSE)
  chrLine <- data.frame(chrom=1:24, length=seqlengths(g)[1:24])
  #rownames(data.bin)<- data.bin$

  if(is.null(plotDir)){plotDir="./"}
  data.bin <- as.data.frame(data.bin)
  data.bin <- data.bin[,!grepl("seg_SD|seg_score",colnames(data.bin)),drop=F]
  if(!any(grepl("segID",colnames(data.bin))) & any(grepl("segName",colnames(data.bin))) ){
    data.bin$segID <- data.bin$segName
  }
  if((!any(grepl("chr",rownames(data.bin)))) & any(grepl("binID",colnames(data.bin)))){
    rownames(data.bin) <- data.bin$binID
  }
  chromnum=sapply(strsplit(rownames(data.bin),'-|_|:'),'[',1)
  chr_name=as.character(na.omit(unique(chromnum)))
  if(grepl('chr',chromnum[1])){
    chromnum=sapply(strsplit(chromnum,'hr'),'[',2) ## if chrommnum with chr
  }else{
    chromnum=chromnum
    }
  chromnum <- gsub("X","23",chromnum)
  chromnum <- gsub("Y","24",chromnum)
  data.bin$chrom <- chromnum
  chr_name <- sapply(strsplit(chr_name,'hr'),'[',2)
  # chr_name <- na.omit(chr_name)
  chr_name_nm <- gsub("X","23",chr_name)
  chr_name_nm <- gsub("Y","24",chr_name_nm)
  chromnum <- na.omit(chromnum)
 
  #
  if(all(!grepl("^Start$",colnames(data.bin)))){
    data.bin$Start <- as.numeric(sapply(strsplit(rownames(data.bin),'-|_|:'),'[',2))
  }
  if(true_loc){
    chrLine <- chrLine[chrLine$chrom %in% as.character(unique(chromnum)),,drop=F] 
    genomeLength <- sum(chrLine)
    if(nrow(chrLine)>1){
      chrPosStart = c(0,chrLine[1:(nrow(chrLine)-1),2])
    }else{chrPosStart = 0}

    chrLine$chrPosStart_abs <- cumsum(chrPosStart)

    chrLine_mpOrder <- match(data.bin$chrom,chrLine$chrom)
    
    data.bin <- cbind(data.bin,chrLine[chrLine_mpOrder,,drop=F])
    
    #data.bin <- merge(data.bin,chrLine,by=merge_by,all.x=T)
    data.bin <- data.bin[order(as.integer(data.bin$chrom),data.bin$chrPosStart_abs,data.bin$Start),]
    if(plot_seg){
       data.bin$seg_start_abs <- as.numeric(sapply(strsplit(data.bin$segName,'-|_|:'),'[',2))+as.numeric(data.bin$chrPosStart_abs)
       data.bin$seg_end_abs <- as.numeric(sapply(strsplit(data.bin$segName,'-|_|:'),'[',3))+as.numeric(data.bin$chrPosStart_abs)
       
     }
    # data.bin <- na.omit(data.bin)
    
    maploc <- data.bin$chrPosStart_abs+data.bin$Start
    intercept_x =chrLine$chrPosStart_abs
    breakss<-cumsum(chrPosStart)+chrLine[1:(nrow(chrLine)),2]*0.5
    which_chrom <- which(as.character(chrLine$chrom) %in% chr_name_nm)
    breakss<-breakss[which_chrom]
    if(plot_seg){
    a<-data.frame(maploc=maploc,data.bin[,value.bin],data.bin[,value.segment],data.bin$segID,data.bin$seg_start_abs,data.bin$seg_end_abs,data.bin$segName)
    colnames(a)=c('axis.loc','value.bin','value.segment','segID',"segStart","segEnd",'segName')
    }else{
      if(all(!grepl("^End$",colnames(data.bin)))){
        data.bin$End <- as.numeric(sapply(strsplit(rownames(data.bin),'-|_|:'),'[',3))
      }
      a<-data.frame(maploc=maploc,data.bin[,c(value.bin,"Start","End")])
      colnames(a)=c('axis.loc','value.bin','Start','End')
      a$length_seg = a$End - a$Start
      
    }
    
    
  }else{
   
    chromnum=chromnum[order(as.integer(gsub("[^0-9]", "", chromnum)))] #order by chr
    chromnum_xaxis = table(chromnum)[order(as.integer(gsub("[^0-9]", "", names(table(chromnum)))))]#order by chr
    chrline=cumsum(chromnum_xaxis)
    maploc=1:nrow(data.bin)
    breakss<-as.numeric(c(0,chrline[-length(chrline)])+chromnum_xaxis*0.5)
    intercept_x =chrline
    data.bin <- data.bin[order(as.integer(data.bin$chrom),data.bin$Start),]
    a<-data.frame(maploc=maploc,data.bin[,value.bin],data.bin[,value.segment],data.bin$segID,data.bin$length_seg,data.bin$segName)
    colnames(a)=c('axis.loc','value.bin','value.segment','segID','length_seg','segName')
    a$rank <- maploc
    if(!is.factor( a$segID)){
      a$segID <- factor(a$segID,levels = sort(unique(as.numeric(as.character(a$segID)))))
    }
    a <- a %>%
      dplyr::group_by(segID) %>%
      dplyr::mutate(segStart=min(rank),segEnd=max(rank))%>%
      as.data.frame()
  }
  if(!is.null(show_specific_seg)){
    a$axis_seg_specific = NA
    whichRow = which(a$segName %in% show_specific_seg)
    a$axis_seg_specific[whichRow] <- maploc[whichRow] 
     line_df.seg <- a[complete.cases(a), ]
  }
  
  if(is.null(add_yline)){
    if(log_ratio){
      a$value.bin <- log1p(a$value.bin)
      a$value.segment <- log1p(a$value.segment)
      add_yline = log1p(1)
    }else{add_yline = 1} 
  }

  
  if(is.null(ylim)){
    if(plot_dots & (!plot_seg)){
      ylim=c(quantile(a$value.bin,0.02),quantile(a$value.bin,0.95))
    }else{
      cen <- median(a$value.segment,na.rm = T)
      dist1 <- abs(cen-max(a$value.segment))
      dist2 <- abs(cen-min(a$value.segment))
      d <- max(dist1,dist2)*1.2
      ylim_l <- ifelse((cen-d) >=0,cen-d,0)
      ylim_h <- cen+d
      ylim=c(ylim_l,ylim_h)
    }
   
  }
 
  if(!is.null(num_cells)){
    ptitle = paste(name.data," (",num_cells," cells)")
  }else{
    ptitle = paste(name.data)
  }
  

  chromnum=sapply(strsplit(rownames(data.bin),'-|_|:'),'[',1)
  chr_name_label=as.character(na.omit(unique(chromnum)))
  chr_name_label <- sapply(strsplit(chr_name_label,'hr'),'[',2)
  chr_name_label <- na.omit(chr_name_label)
  if(length(chr_name_label)>1){
    chr_name_label[as.character(chr_name_label)%in%c("19","21")] <- ""
  }

    p1 <- ggplot(a) +
    #geom_point(aes(x = axis.loc, y = value.bin), alpha = 1, shape = 16,size=1) +
    geom_hline(yintercept =add_yline,col="grey",linetype = 2)+
    theme_bw() +
    labs(x = "", y = ylab)+
    coord_cartesian(ylim = ylim)+
    ggtitle(ptitle)+
    #geom_vline(xintercept =intercept_x,col="grey")+
    theme(plot.margin = unit(c(0.5,0,0,0.5), 'lines'), #c(top, right, bottom, left)
      panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.position = "none",plot.title = element_text(size = 15),
          axis.title=element_text(size=15),
          axis.text = element_text(size = 15))+
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE),breaks=breakss,labels=chr_name_label)
   
    if(length(intercept_x)>1){

      odd_positions <- intercept_x[seq(1, length(intercept_x), by = 2)]
     
      even_positions <- intercept_x[seq(2, length(intercept_x), by = 2)]

      if(length(odd_positions)!=length(even_positions)){
        nrow <- min(length(odd_positions),length(even_positions))
      }else{nrow=length(odd_positions)}
      ## highlight region data
      rects <- data.frame(start=odd_positions[1:nrow], end=even_positions[1:nrow])
      p1 <- p1+
        geom_rect(data=rects, inherit.aes=FALSE, 
                  aes(xmin=start, xmax=end, ymin=-Inf,ymax=Inf),
                  color="transparent", fill="grey50", alpha=0.3)
      }
 # brks <- unique(c(seq(0.1,0.6,length.out =2),seq(0.8,1.2,length.out =2),seq(1.5,2,length.out =3),2.5,3))
    if(color_hist_gradient){
      if(is.null(plot_colors)){
        colorss <- c("#023e8a","grey50", "#cb793a", "#9a031e","#6a040f","#370617")
        #colorss <- c("black","#023e8a","grey50","#fcdc4d", "#cb793a", "#9a031e","#5f0f40")
        
      }else{
        colorss <- plot_colors
      }
    }else{
      if(is.null(plot_colors)){
        colorss <- "#472d30"
      }else{
        colorss <- plot_colors[1]
      }
      
    }

  
  if(plot_dots){

    if(color_dot){

      p1 <- p1+
        geom_point(aes(x = axis.loc, y = value.bin,colour=value.segment), alpha = 1, shape = 16,size=0.5)+
        scale_color_gradientn(name = "",colours = colorss,#c("#003049", "#669bbc","#7e7f83","#c1121f","#780000"),
                            limits  = color_limit,
                            oob = scales::squish)
    }else{
      
      p1 <- p1+
        geom_point(aes(x = axis.loc, y = value.bin),colour="darkgrey", alpha = 1, shape = 16,size=0.5)
    
    }
  }

  
  if(plot_seg){
    if(color_seg_gradient){
      p1 <- p1+
        geom_segment(
          data = a, 
          mapping = aes(x=segStart, y=value.segment, xend=segEnd, yend=value.segment,colour=value.segment), 
          na.rm =T,size = 1,
          inherit.aes = FALSE)+
        scale_color_gradientn(name = "",colours = colorss,#c("#003049", "#669bbc","#7e7f83","#c1121f","#780000"),
                              limits  = color_limit,
                              oob = scales::squish)
    }else{ 
    p1 <- p1+
      #geom_line(aes(x = maploc,y=value.segment),col="red")+
      
      geom_segment(
        data = a, 
        mapping = aes(x=segStart, y=value.segment, xend=segEnd, yend=value.segment), 
        colour=seg_col,na.rm =T,size = 1,
        inherit.aes = FALSE)
    }
    
    if(!disconnected){
      dat2<- data.frame(x=a$segEnd[1:(length(a$segEnd)-1)], y=a$value.segment[1:(length(a$segEnd)-1)],
                              xend=a$segStart[2:length(a$segEnd)], yend=a$value.segment[2:length(a$segEnd)])
      p1 <- p1 +
        geom_segment(
          data = dat2,
          mapping = aes(x=x, y=y,
                        xend=xend, yend=yend), 
          colour=seg_col,na.rm =T,size = 1,
          inherit.aes = FALSE)
    }
  }


  if(!is.null(show_specific_seg)){ 
    p1=p1+geom_segment(
      data = line_df.seg, 
      mapping = aes(x=segStart, y=value.segment, xend=segEnd, yend=value.segment), 
      col=specific_seg_color,na.rm =T,size = 1.5,
      inherit.aes = FALSE)
  }
  if(!is.null(label_gene)){
    gene_loc <- annotate_gene(label_gene,genome=genome,annotations=annotations)
    gene_chr <- as.character(seqnames(gene_loc))
    gene_chr <- sapply(strsplit(gene_chr,'hr'),'[',2)
    if(true_loc){
    label_x <- chrLine[chrLine$chrom==gene_chr,"chrPosStart_abs"]+start(gene_loc)
    label_y <- unique(data.bin[data.bin$chrom==gene_chr&data.bin$seg_start_abs<=label_x&data.bin$seg_end_abs>=label_x,value.segment])[1]
    }else{
      g_start <- start(ranges(gene_loc))
      g_end <- end(ranges(gene_loc))
      label_x <- which(data.bin$chrom==gene_chr&data.bin$Start>=g_start)[1]
      label_y <- data.bin[data.bin$chrom==gene_chr&data.bin$Start>=g_start,value.segment][1]
      
    }

    p1=p1+
      geom_point(aes(x = label_x, y = label_y), alpha = 1, shape = 17,size=1,colour=label_gene_color)+
      geom_text(aes(label_x,label_y,label = label_gene),color=label_gene_color,vjust = -1)+
      geom_segment(aes(x = label_x, y = label_y, xend = label_x, yend = label_y -0.5),
                   color = label_gene_color,  # 指定指引线颜色
                   linetype = "dashed")  
    
  }
  
  if(plot_hist){
    
  breakshist <- seq(min(ylim),max(ylim),length.out=80)
  #segSet = unique(data.bin[,c("binID","segID",value.segment,"segName"),drop=F]) ###
  if(true_loc){
    if(plot_seg){
      a$segLen <- as.numeric(a$segEnd)-as.numeric(a$segStart)
    }else{
     a$segLen <- as.numeric(a$length_seg) #Here, length_seg is the length of peak
    }
  }else{
    a$segLen <- as.numeric(a$length_seg) 
  }  
  
  if(!plot_seg){
    a$segID <- seq(1:nrow(a))
  }
  segLen_sum <- sum(unique(a[,c("segID","segLen")]$segLen))
  a$wei <- a$segLen/segLen_sum
  if(plot_seg){
    p_x <- "value.segment"
  }else{
    p_x <- "value.bin"
  }


  if(color_hist){
    if(color_hist_gradient){
      p3 <- ggplot(a,aes(x=get(p_x),fill = ..x..))+
        geom_histogram(aes(y=..density.., weight = wei),breaks=breakshist)+ #color="#E69F00",fill="wheat",
        labs(x = "", y = "",fill="")+
        ggtitle("Summary histogram")+
        coord_cartesian(xlim=ylim)+
        theme_bw() +
        theme_classic()+
        theme(plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
              # legend.position = "bottom",
              plot.title = element_text(size = 15),
              axis.title=element_text(size=15),
              axis.text = element_text(size = 12))+
        coord_flip(xlim=ylim)+
      scale_fill_gradientn(name = "",colours = colorss,#c("#003049", "#669bbc","#7e7f83","#c1121f","#780000"),
                           limits  = color_limit,
                           oob = scales::squish)
    }else{
      p3=ggplot(a,aes(x=get(p_x)))+
        geom_histogram(aes(y=..density.., weight = wei),color=colorss,fill=colorss,breaks=breakshist)+ #color="#E69F00",fill="wheat",
        labs(x = "", y = "",fill="")+
        ggtitle("Summary histogram")+
        coord_cartesian(xlim=ylim)+
        theme_bw() +
        theme_classic()+
        theme(plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
              # legend.position = "bottom",
              plot.title = element_text(size = 15),
              axis.title=element_text(size=15),
              axis.text = element_text(size = 12))+
        coord_flip(xlim=ylim)
  }
  }else{
    p3=ggplot(a,aes(x=get(p_x)))+
      geom_histogram(aes(y=..density.., weight = wei),breaks=breakshist)+ #color="#E69F00",fill="wheat",
      labs(x = "", y = "",fill="")+
      ggtitle("Summary histogram")+
      coord_cartesian(xlim=ylim)+
      theme_bw() +
      theme_classic()+
      theme(plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
            # legend.position = "bottom",
            plot.title = element_text(size = 15),
            axis.title=element_text(size=15),
            axis.text = element_text(size = 12))+
      coord_flip(xlim=ylim)
  }
  pp=ggarrange(p1, p3, align ="h",ncol = 2, nrow = 1,widths = c(15,4),heights=3 )
  p_ls <- list(ggarranged_p=pp,p1=p1,p2=p3)
  }else{
    pp <- p1
  p_ls <- list(p1=p1)
  }
  if(outPlot){
    ggsave(paste0(plotDir,"/",name.data,".",device),pp,height = height,width=width,dpi = 300,device = device)
  }
  return(p_ls)
}


#' @title seg_plot_dot()
seg_plot_dot <- function(data.bin,name.data="",plotDir=NULL,ylab="Count",genome="hg38",outplot=TRUE){
  suppressMessages(require(BSgenome))
  g <- getBSgenome(genome, masked=FALSE)
  chrLine <- data.frame(chrom=1:24, length=seqlengths(g)[1:24])

  
  if(is.null(plotDir)){plotDir="./"}
  #order data.bin by chromosome
  row <- data.frame(rn=rownames(data.bin))
  
  chromnum=sapply(strsplit(rownames(data.bin),'-|_|:'),'[',1)
  row$seqnames <- chromnum
  row$start <- sapply(strsplit(rownames(data.bin),'-|_|:'),'[',2)
  row$end <-  sapply(strsplit(rownames(data.bin),'-|_|:'),'[',3)

    
  if(grepl('chr',chromnum[1])){
    chromnum=sapply(strsplit(chromnum,'hr'),'[',2) ## if chrommnum with chr
  }else{
    chromnum=chromnum
  }
  chr_name=as.character(unique(chromnum))
  chr_name= na.omit(chr_name)
  
  chromnum <- gsub("X","23",chromnum)
  chromnum <- gsub("Y","24",chromnum)
  row$chrom <- chromnum
  data.bin$chrom <- chromnum
  
  chrLine <- chrLine[chrLine$chrom %in% as.character(unique(chromnum)),,drop=F] 
  genomeLength <- sum(chrLine)
  if(nrow(chrLine)>1){
    chrPosStart = c(0,chrLine[1:(nrow(chrLine)-1),2])
  }else(
    chrPosStart = 0
  )
  
  chrLine$chrPosStart_abs <- cumsum(chrPosStart)
  
  
  row <- merge(row,chrLine,by="chrom",all.x=T)
  row <- row[order(as.integer(row$chrom),as.numeric(row$chrPosStart_abs)),]
  data.bin <- as.data.frame(data.bin)[row$rn,,drop=F]
  
  chrLine_mpOrder <- match(data.bin$chrom,chrLine$chrom)
  data.bin <- cbind(data.bin,chrLine[chrLine_mpOrder,,drop=F])

  ## scatter plot
  chr_name_nm <- gsub("X","23",chr_name)
  chr_name_nm <- gsub("Y","24",chr_name_nm)
  which_chrom <- which(as.character(chrLine$chrom) %in% chr_name_nm)
  data.bin <- na.omit(data.bin)
  if(all(!grepl("^Start$",colnames(data.bin)))){
    data.bin$Start <- as.numeric(sapply(strsplit(rownames(data.bin),'-|_|:'),'[',2))
  }
  maploc <- data.bin$chrPosStart_abs+data.bin$Start
  intercept_x =chrLine$chrPosStart_abs
  breakss<-cumsum(chrPosStart)+chrLine[1:(nrow(chrLine)),2]*0.5
  which_chrom <- which(as.character(chrLine$chrom) %in% chr_name_nm)
  breakss<-breakss[which_chrom]
  
  chr_name_label <- chr_name
  chr_name_label[as.character(chr_name_label)%in%c("19","21")] <- ""
  
  
  a<-data.frame(maploc=maploc,data.bin)
  colnames(a)=c('axis.loc','value.bin')
  ylim=c(min(a$value.bin),max(a$value.bin))
  p1 <- ggplot(a) +
    geom_point(aes(x = axis.loc, y = value.bin), alpha = 1, shape = 16,size=0.6) +
    theme_bw() +
    labs(x = "", y = ylab)+
    ylim(ylim)+
    ggtitle(paste(name.data))+
    #geom_line(aes(x = maploc,y=value.segment),col="red")+
    geom_vline(xintercept =intercept_x,col="grey")+
    #geom_hline(yintercept =1,col="grey",linetype = 2)+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.position = "bottom",plot.title = element_text(size = 15),
          axis.title=element_text(size=15),
          axis.text = element_text(size = 15))+
    scale_x_continuous(breaks=breakss,labels=chr_name_label)
 
  breakshist <- seq(min(ylim),max(ylim),length.out=100)
  a$segLen <- as.numeric(a$segEnd)-as.numeric(a$segStart)
  segLen_sum <- sum(unique(a[,c("segID","segLen")]$segLen))
  a$wei <- a$segLen/segLen_sum

  p3=ggplot(a,aes(x=value.bin))+
    geom_histogram(aes(y=..density.., weight = wei),color="#E69F00",fill="wheat",breaks=breakshist)+
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
  if(outplot){
    ggsave(paste0(plotDir,"/",name.data,".png"),pp,height = 3,width=14,dpi = 300)
  }else{
    return(pp)
  }
 
}

#' @title annotate_gene()
#' @description genomic annotation for gene symbol
#' @param gene_name gene symbol
#' @param genome string of reference genome
#' @param annotations GRanges object
#' @return a GRanges object of gene annotation
#' 
annotate_gene <- function(gene_name,genome="hg38",annotations=NULL){
  suppressMessages({
    library(BSgenome)
    library(EnsDb.Hsapiens.v86)})
  if(is.null(annotations)&genome%in%c("hg38","GRCh38")){
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  }
  anno <- annotations[annotations$gene_name%in%gene_name,,drop=F]
  gene_loc <- CollapseToLongestTranscript(anno)
  return(gene_loc)
}

#' @title CollapseToLongestTranscript()
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @import data.table
CollapseToLongestTranscript <- function(ranges) {
  library(data.table)
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

#' @title plot_point()
plot_point <- function(df,x_aes,y_aes,colour_item=NULL,size_item=NULL,size_lab=size_item,
                       shape_item=NULL,shape_label = NULL,
                       color_value,color_gradient=FALSE,
                       color_gradientn_name="CN",
                       color_gradientn_limits= c(1,6),
                       x_lab,y_lab,
                       label_term=NULL,gtitle,log_x=FALSE,log_y=FALSE,
                       legend_position="rignt",
                       legend_color=NULL,
                       lm_test=TRUE,
                       outfile=NULL,
                       out_width = 6,out_height = 5,
                       xlim.set=TRUE,
                       ylim.set = TRUE){
  suppressPackageStartupMessages({
    library(ggpubr)
    library(ggpmisc)
    library(ggplot2)
    library("colorspace")
  })
  if(xlim.set){
    if(log_x){
      x_min = log1p(min(df[,x_aes],na.rm = T))
      #if(is.infinite(x_min))x_min=0
      x_max = log1p(max(df[,x_aes],na.rm = T))
    }else{
      x_min = min(df[,x_aes],na.rm = T)
      x_max = max(df[,x_aes],na.rm = T)
    }
    xlim_countrna <-c(x_min,x_max)
  }else{
    xlim_countrna=NULL
  }
if(ylim.set){
  if(log_y){
    y_min = log1p(min(df[,y_aes],na.rm = T))
    #if(is.infinite(y_min))y_min=0
    y_max = log1p(max(df[,y_aes],na.rm = T))
  }else{
    y_min = min(df[,y_aes],na.rm = T)
    y_max = max(df[,y_aes],na.rm = T)
  }
  ylim_countrna <- c(y_min,y_max)
}else{
  ylim_countrna=NULL
}

  
  
  
  if(log_x&log_y){
    p <- ggplot(df,aes(log1p(get(x_aes)),log1p(get(y_aes))))
  }else if(log_x &(!log_y)){
    p <- ggplot(df,aes(log1p(get(x_aes)),get(y_aes)))
  }else if(!(log_x) &log_y){
    p <- ggplot(df,aes(get(x_aes),log1p(get(y_aes))))
  }else{p <- ggplot(df,aes(get(x_aes),get(y_aes)))}
  
  aes_params <- ifelse(is.null(size_item),list(colour = colour_item,shape_item),list(colour = colour_item,shape=shape_item,size=size_item))
  p <- p+
    geom_point(aes_string(colour = colour_item,size=size_item,shape=shape_item))+
    labs(x=x_lab,y=y_lab,colour=colour_item,size=size_lab,shape=shape_label)
  
  # if(!is.null(colour_item)){
  #   if(!is.null(size_item)){
  #     p <- p+
  #       geom_point(aes(colour = get(colour_item),size=get(size_item)))+
  #       labs(x=x_lab,y=y_lab,colour=colour_item,size=size_lab)
  #   }else{
  #     p <- p+
  #       geom_point(aes(colour = get(colour_item)))+
  #       labs(x=x_lab,y=y_lab,colour=colour_item)
  #   }
  # }else{
  #   p <- p+
  #     geom_point()+
  #     labs(x=x_lab,y=y_lab)
  # }
  if(xlim.set){
    p <- p+
      xlim(xlim_countrna)
  }
  if(ylim.set){
    p <- p+
      ylim(ylim_countrna)
  }
  p <- p+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))+
    ggtitle(gtitle)
  if(color_gradient){
    p <- p+
      scale_color_gradientn(name = color_gradientn_name,colours = color_value,limits  = color_gradientn_limits,oob = scales::squish)
  }else if(!is.null(colour_item)){
    p <- p+
      scale_color_manual(values=color_value)
  }
  if(!is.null(label_term)){
    p <- p+
      ggrepel::geom_text_repel(
        aes(label=get(label_term)),size = 2.5, 
        box.padding = 0.5, 
        point.padding = 0.5, 
       # min.segment.length = 0.5, #短线段可以省略
        segment.color = "black", #segment.colour = NA, 不显示线段
        show.legend = F
      )
      # geom_label_repel(aes(label = get(label_term), size = NULL,),
      #                       max.overlaps = getOption("ggrepel.max.overlaps", default = 20),show_guide=F)
  }
  if(!is.null(legend_color)){
    
    p2 <- p+
      theme(legend.position=legend_position)+
      guides(colour = legend_color )
  }else{p2<- p}
  if(lm_test){
    p2 <- p2+
      geom_smooth(method="lm", se=F, fullrange=F, colour="grey70", size=0.5,linetype=0)+
      stat_poly_eq(aes(label =  paste(..adj.rr.label..,..p.value.label.., sep = "~~~~")))
  }
  if(!is.null(outfile)){
    if(!grepl(".pdf",outfile)){outfile=paste0(outfile,".pdf")}
    ggsave(outfile,p2,width = out_width,height = out_height)
  }
  return(p2)
  
}


#' @title abs_pos()
#' @description add columns of chrosome absolute start and end for "segName" of a data.frame
#' segName can be chr1-100037672-100038546 or chr1_100037672_100038546 format
abs_pos <- function(data.bin,genome="hg38",true_loc=TRUE){
    suppressMessages({
    library(BSgenome)
    library(dplyr)})
  g <- getBSgenome(genome, masked=FALSE)
  chrLine <- data.frame(chrom=1:24, length=seqlengths(g)[1:24])
  chromnum=sapply(strsplit(data.bin$segName,'-|_|:'),'[',1)
  chr_name=as.character(unique(chromnum))
  if(grepl('chr',chromnum[1])){
    chromnum=sapply(strsplit(chromnum,'hr'),'[',2) ## if chrommnum with chr
  }else{
    chromnum=chromnum
    }
  chromnum <- gsub("X","23",chromnum)
  chromnum <- gsub("Y","24",chromnum)
  data.bin$chrom <- chromnum
  chr_name <- sapply(strsplit(chr_name,'hr'),'[',2)
  chr_name_nm <- gsub("X","23",chr_name)
  chr_name_nm <- gsub("Y","24",chr_name_nm)

  if(true_loc){
    chrLine <- chrLine[chrLine$chrom %in% as.character(unique(chromnum)),,drop=F] 
    genomeLength <- sum(chrLine)
    if(nrow(chrLine)>1){
      chrPosStart = c(0,chrLine[1:(nrow(chrLine)-1),2])
    }else(
      chrPosStart = 0
    )
    chrLine$chrPosStart_abs <- cumsum(chrPosStart)
    data.bin <- merge(data.bin,chrLine,by="chrom",all.x=T)
    data.bin <- data.bin[order(as.integer(data.bin$chrom),data.bin$chrPosStart_abs),]
    data.bin$segStart <- as.numeric(sapply(strsplit(data.bin$segName,'-|_|:'),'[',2))+as.numeric(data.bin$chrPosStart_abs)
    data.bin$segEnd <- as.numeric(sapply(strsplit(data.bin$segName,'-|_|:'),'[',3))+as.numeric(data.bin$chrPosStart_abs)
    
    #data.bin$Start <- sapply(strsplit(data.bin$segName,'-|_|:'),'[',2) #for segment start( not bin)
    
    #maploc <- data.bin$chrPosStart_abs+data.bin$Start
    intercept_x =chrLine$chrPosStart_abs
    breakss<-cumsum(chrPosStart)+chrLine[1:(nrow(chrLine)-1),2]*0.5
    which_chrom <- which(as.character(chrLine$chrom) %in% chr_name_nm)
    breakss<-breakss[which_chrom]

    #data.bin$maploc <- maploc

  }else{
   
    chromnum=chromnum[order(as.integer(gsub("[^0-9]", "", chromnum)))] #order by chr
    chromnum_xaxis = table(chromnum)[order(as.integer(gsub("[^0-9]", "", names(table(chromnum)))))]#order by chr
    chrline=cumsum(chromnum_xaxis)
    maploc=1:nrow(data.bin)
    breakss<-as.numeric(c(0,chrline[-length(chrline)])+chromnum_xaxis*0.5)
    intercept_x =chrline
    data.bin$maploc <- maploc
    data.bin$rank <- maploc
    data.bin <- data.bin %>%
      dplyr::group_by(segName) %>%
      dplyr::mutate(segStart=min(rank),segEnd=max(rank))%>%
      as.data.frame()
  }
return(data.bin)
}

#' @title seg_compare()
#‘ @df data.frame input with c("segID","SegMean","segName","Start") on columns
#                         segID  SegMean                segName      Start   
#chr1-100037672-100038546     3 2.578573 chr1_9823588_121519947    100037672 
#chr1-100132556-100133506     3 2.578573 chr1_9823588_121519947    100132556 
seg_compare <- function(df1,df2,genome="hg38",ylim=NULL,ylab="Relative CN",ggtitle="",
  hline=c(1,2),
  seg_col=c("darkgrey","red"),
  alpha=0.5){
  suppressMessages({
    library(BSgenome)
    library(dplyr)
    library(ggplot2)
    })
  g <- getBSgenome(genome, masked=FALSE)
  chrLine <- data.frame(chrom=1:24, length=seqlengths(g)[1:24])
  if(nrow(chrLine)>1){
      chrPosStart = c(0,chrLine[1:(nrow(chrLine)-1),2])
    }else{
      chrPosStart = 0
    }

    chrLine$chrPosStart_abs <- cumsum(chrPosStart)

  df1 <- abs_pos(df1)
  df2 <- abs_pos(df2)

  if(is.null(ylim)){
    ylim=c(min(a$SegMean),max(a$SegMean)+1)
  }
  chromnum <- sapply(strsplit(df1$segName,'-|_|:'),'[',1)
  chr_name <- as.character(unique(chromnum))
  chr_name <- sapply(strsplit(chr_name,'hr'),'[',2)
  chr_name_nm <- gsub("X","23",chr_name)
  chr_name_nm <- gsub("Y","24",chr_name_nm)
  chr_name_label <- chr_name
  chr_name_label[as.character(chr_name_label)%in%c("19","21")] <- ""
  breakss<-cumsum(chrPosStart)+chrLine[1:(nrow(chrLine)-1),2]*0.5
  which_chrom <- which(as.character(chrLine$chrom) %in% chr_name_nm)
  breakss<-breakss[which_chrom]

  intercept_x =chrLine$chrPosStart_abs
  p1 <- ggplot() +
    #geom_point(aes(x = axis.loc, y = value.bin), alpha = 1, shape = 16,size=1) +
    geom_hline(yintercept =hline,col="grey",linetype = 2)+
    theme_bw() +
    labs(x = "", y = ylab)+
    ylim(ylim)+
    ggtitle(ggtitle)+
    geom_vline(xintercept =intercept_x,col="grey")+
    theme(plot.margin = unit(c(0.5,0,0,0.5), 'lines'), #c(top, right, bottom, left)
      panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
          legend.position = "none",plot.title = element_text(size = 15),
          axis.title=element_text(size=15),
          axis.text = element_text(size = 15))+
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE),breaks=breakss,labels=chr_name_label)

  p1 <- p1 +
    #geom_line(aes(x = maploc,y=value.segment),col="red")+ 
    geom_segment(
      data = df1, 
      mapping = aes(x=segStart, y=SegMean, xend=segEnd, yend=SegMean), 
      col=seg_col[1],na.rm =T,size = 1, alpha = alpha,
      inherit.aes = FALSE) +
    geom_segment(
      data = df2, 
      mapping = aes(x=segStart, y=SegMean, xend=segEnd, yend=SegMean), 
      col=seg_col[2],na.rm =T,size = 1,alpha = alpha,
      inherit.aes = FALSE)
return(p1)
}

#' @title seg_split()
#' @description split intersecting segments
#' @param binSet data.frame with bins in rowname
#'                                     binID Chromosome   Start     End
# chr1-826972-1013895   chr1-826972-1013895       chr1  826972 1013895
# chr1-1019104-1116833 chr1-1019104-1116833       chr1 1019104 1116833
# chr1-1121514-1350035 chr1-1121514-1350035       chr1 1121514 1350035
#' 
seg_split <- function(segs_input,binSet=NULL,split_by="_|-|:"){
  suppressMessages({
    require(data.table)
    require(Matrix)
  })
  bed_data <- tstrsplit(segs_input, split_by)
  bed_data <- as.data.frame(bed_data)
  colnames(bed_data) <- c("chromosome", "start", "end")
  bed_data$start <- as.numeric(bed_data$start)
  bed_data$end <- as.numeric(bed_data$end)
  bed_data <- bed_data[order(bed_data$chromosome, bed_data$start, bed_data$end), ]
if(! any(is(binSet) %in% "data.frame")){
  binSet <- data.frame(row.names = binSet,binID = binSet)
  binSet$Chromosome <- sapply(strsplit(rownames(binSet),split_by),"[",1)
  binSet$Start <- as.numeric(sapply(strsplit(rownames(binSet),split_by),"[",2))
  binSet$End <- as.numeric(sapply(strsplit(rownames(binSet),split_by),"[",3))
}
  mt_tag <- Matrix(0,nrow = nrow(binSet),ncol = nrow(bed_data),sparse = TRUE)
  rownames(mt_tag) <- rownames(binSet)
  for(i in 1:ncol(mt_tag)){
    chr <- bed_data$chromosome[i]
    seg_start <- bed_data$start[i]
    seg_end <- bed_data$end[i]
    strat_idx <- which(binSet$Start==seg_start & binSet$Chromosome==chr)
    if(length(strat_idx)==0){
      strat_idx <- which(binSet$Chromosome==chr & binSet$Start<= seg_start)
      strat_idx <- strat_idx[length(strat_idx)]
      if(length(strat_idx)==0){
        strat_idx <- which(binSet$Chromosome==chr)[1]
      }
    }
    end_idx <- which(binSet$Chromosome==chr & binSet$End==seg_end)
    if(length(end_idx)==0){
      end_idx <- which(binSet$Chromosome==chr & binSet$End>= seg_end)[1]
      if(length(end_idx)==0|is.na(end_idx)){
        end_idx <- which(binSet$Chromosome==chr)
        end_idx <- end_idx[length(end_idx)]
      }
    }
    mt_tag[strat_idx:end_idx,i] <- 1
    
  }

  colnames(mt_tag) <- paste(bed_data$chromosome, bed_data$start, bed_data$end,sep ="_")
  mt_tag <- mt_tag[rowSums(mt_tag)>0,,drop=F]
  row_dist <- c()
  for(j in 1:(nrow(mt_tag)-1)){
    dist_row <- sum(abs(mt_tag[j+1,]-mt_tag[j,]))
    row_dist <- c(row_dist,dist_row)
  }
  row_dist <- c(row_dist,1)
  brek_point <- ifelse(row_dist>0,1,0)
  seg_end = which(brek_point==1)
  seg_n <- length(seg_end)
  segSet <- c()
  seg_start <- 1
  binSet$segName <-NA
  binSet$length_seg <- NA
  for(k in 1:seg_n){
    end_k <- seg_end[k]
    chr_k <- binSet$Chromosome[seg_start] 
    seg_start_loc <- binSet$Start[seg_start]
    seg_end_loc <- binSet$End[end_k]
    seg_len_k <- seg_end_loc-seg_start_loc#
    seg_name_k <- paste(chr_k,seg_start_loc,seg_end_loc,sep="_")
    segSet_k <- cbind(seg_name_k,chr_k,seg_start_loc,seg_end_loc,seg_len_k)
    segSet <-rbind(segSet,segSet_k)
   
    binSet$segName[seg_start:end_k] <- seg_name_k
    binSet$length_seg[seg_start:end_k] <- as.numeric(seg_len_k)
    seg_start <- end_k+1
    
    rm(chr_k,seg_start_loc,seg_end_loc)
  }
  segSet = as.data.frame(segSet)
  colnames(segSet)=c("segName","Chromosome","Start","End","length_seg")
  columns <-c("Start","End","length_seg")
  segSet[,columns] <- sapply(segSet[,columns], as.numeric)
  
  res <- list(seg_binLevel = binSet,seg_segLevel= segSet)
  return(res)
}

#' @title seg_overlap_prop()
#' @description get the percentage of overlap in each region
#' @param region1 and region2, genomic region with format "chrx-xxx-xxx"
seg_overlap_prop <- function(region1,region2,split_by="_|-|:"){
  suppressMessages({require(data.table)})
  segs <- c(region1,region2)
  bed_data <- tstrsplit(segs, split_by)
  bed_data <- as.data.frame(bed_data)
  colnames(bed_data) <- c("chromosome", "start", "end")
  bed_data$start <- as.numeric(bed_data$start)
  bed_data$end <- as.numeric(bed_data$end)
  bed_data <- bed_data[order(bed_data$chromosome, bed_data$start, bed_data$end), ]
  bed_data$length <- bed_data$end-bed_data$start
  if(bed_data$chromosome[1]!=bed_data$chromosome[2]){prop=0}else{
    if(bed_data$start[2]>bed_data$end[1]){prop=0}else{
      if(bed_data$end[2]<=bed_data$end[1]){
        length_ovelap <- bed_data$length[2]
      }else{
        length_ovelap <- bed_data$end[1]-bed_data$start[2]
      }
      prop <- length_ovelap/bed_data$length
    }
  }
  
  return(prop)
}

#' @title seg_slimByAlign()
#' @description Align and unify highly similar fragments
#' @param Start_site_error the start site distance for same segment (bp)
#' @param End_site_error 
seg_slimByAlign <- function(segs, 
                            split_by = "_|-|:",
                            Start_site_error=1e6,
                            End_site_error=1e6){
  suppressMessages({require(data.table)})
  bed_data <- tstrsplit(segs, split_by)
  bed_data <- as.data.frame(bed_data)
  colnames(bed_data) <- c("chromosome", "start", "end")
  bed_data$start <- as.numeric(bed_data$start)
  bed_data$end <- as.numeric(bed_data$end)
  bed_data$raw <- segs
  bed_data$new <- bed_data$raw
  bed_data$length <- bed_data$end -  bed_data$start
  bed_data <- bed_data[order(bed_data$chromosome, bed_data$start, bed_data$end), ]
  
  replace_idx <- c()
  for(k in 1:(length(segs)-1)){
    if(!k %in% replace_idx){
      segk <- bed_data$raw[k]
      start_k <- bed_data$start[k]
      end_k <- bed_data$end[k]
      chr_k <- bed_data$chromosome[k]
    }else{
      segk <- bed_data$new[k]
      segk_split <- unlist(tstrsplit(segk, split_by))
      start_k <- as.numeric(segk_split[2])
      end_k <-as.numeric(segk_split[3])
      chr_k <- segk_split[1]
    }
    len_k <- end_k- start_k
    
    for(l in (k+1):length(segs)){
      if(!l %in% replace_idx){
        segl <- bed_data$raw[l]
        start_l <- bed_data$start[l]
        end_l <- bed_data$end[l]
        chr_l <- bed_data$chromosome[l]
        
      }else{
        segl <-  bed_data$new[l]
        segl_split <- unlist(tstrsplit(segl, split_by))
        start_l <- as.numeric(segl_split[2])
        end_l <-as.numeric(segl_split[3])
        chr_l <- segl_split[1]
        }
      len_l <- end_l- start_l
      
      dist_start <- abs(start_l - start_k)
      dist_end <-  abs(end_l - end_k)
      if(chr_k==chr_l) {
        which_keep <- which.min(c(len_k,len_l))
        
      idx_keep <- ifelse(which_keep==1,k,l) 
      idx_change <- ifelse(which_keep==1,l,k) #replace the long seg with short seg
      #replace the longer segment with the shorter one if over 95% region are identical
      if(dist_start<=Start_site_error & dist_end <= End_site_error){
        bed_data$new[idx_change] <- bed_data$raw[idx_keep]
        replace_idx <- cbind(replace_idx,idx_change)
      }
    }
    }
  }
  bed_data$replaced <- bed_data$raw!=bed_data$new
  return(bed_data)
}

#' @title correlat_plot()
#' @param seu_obj Seurat object
#' @param mt_atac (Optional)scATAC count matrix (normalized), peaks on row and cells on column
#' @param mt_rna (Optional)scRNA count matrix (normalized), genes on row and cells on column
#' @param cellMeta data.frame with cell group annotation on the second column, cellID on the rowname and first column
#' @param CN_res final CNV result output by runEiCNVs()
correlat_plot <- function(seu_obj=NULL,
                          mt_atac=NULL,
                          mt_rna=NULL,
                          cnv_rna=NULL,                          
                          assay.atac="peaks",
                          assay.rna = "RNA",
                          cellMeta,
                          CN_res,
                          FiltSegbyLen =TRUE,
                          segLen_lim = 1e6,
                          outdir = "./",
                          dotPlot4Seg=TRUE,
                          label_term= NULL#"Cluster"#"Chromosome"
                          ){
  suppressMessages({
    library(ggplot2)
    library(ggpubr)
    library(tidyr)
    library(ggsci)
    library(Seurat)
    library(Signac)
    library(dplyr)
    library(plyranges)
    library(ggrepel)
    library(GenomicRanges)
  })
  if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
  colnames(cellMeta) <- c("cellID","Cluster")
  rownames(cellMeta) <- cellMeta$cellID
  if(is.null(mt_atac) & !is.null(seu_obj)){
    mt_atac <- seu_obj[[assay.atac]]@counts
  }
   mt_atac <- mt_atac[rowSums(mt_atac>0)>1,,drop=F]
   mt_atac<- mt_atac[,rownames(cellMeta)]

  nCount_peaks <- colSums(mt_atac)
  nFeature_peaks <- colSums(mt_atac>0)
  cellMeta$nCount_ATAC <- nCount_peaks
  cellMeta$nFeature_ATAC <- nFeature_peaks

  if(is.null(mt_rna) & !is.null(seu_obj)){
    mt_rna <- seu_obj[[assay.rna]]@counts   
  }
  if(!is.null(mt_rna)){
    mt_rna <- mt_rna[rowSums(mt_rna>0)>1,,drop=F]
    mt_rna <- mt_rna[,rownames(cellMeta)]
    nCount_RNAc <- colSums(mt_rna)
    nFeature_RNAc <- colSums(mt_rna>0)
    cellMeta$nCount_RNA <- nCount_RNAc
    cellMeta$nFeature_RNA <- nFeature_RNAc
  }
  
  mt_atac2_bulk <- pseudo_bulk_v2(mt_atac,group_ano = cellMeta,
                                  method ="mean",adjust_zero=TRUE)
  Ncol <- ncol(mt_atac2_bulk)
  mt_atac2_bulk$peakID <- rownames(mt_atac2_bulk)
  mt_atac2_bulk_long <- tidyr::gather(mt_atac2_bulk,Cluster,Clmean_nCountATAC_peak,0:Ncol, factor_key=TRUE)
  mt_atac2_bulk_long$sampleID <- sapply(strsplit(as.character(mt_atac2_bulk_long$Cluster),'_|:|-'),'[',1)
  
  ##(2) process to segment matrix 
  #combine segment CNV of all clone
  # clonal ATAC cnv res
  segMT_long <- do.call(rbind, lapply(seq_along(cnvres), function(i) {
    data_i <- cnvres[[i]]$clonalest
    name_i <- names(cnvres)[i]
    binDF <- do.call(rbind,lapply(seq_along(data_i), function(j){
      clonej <- names(data_i)[j]
      dfj <- data_i[[j]]$input_BinSegRatio
      dfj$ploidy <- data_i[[j]]$ploidy
      dfj$clone <- clonej
      dfj_cn <-  data_i[[j]]$seg.dat
      dfj_cn <- dfj_cn[,c("segName","relativeCN","integerCN")]
      dfj<- left_join(dfj,dfj_cn,by="segName")
      dfj <- dfj %>%
        dplyr::add_count(clone,segName) %>%
        as.data.frame()
      colnames(dfj)[colnames(dfj)=="n"] <- "nSeg"
      return(dfj)
    }))
    binDF$sampleID <- name_i
    binDF$clone_final <-  paste(binDF$sampleID,binDF$clone,sep="_")
    binDF <- binDF[,!grepl("segName_raw$",colnames(binDF))]
    return(binDF)
  }))

  colnames(segMT_long) <- gsub("clone_final","Cluster",colnames(segMT_long))
  #add binID annottaion to peakID
  require(GenomicRanges)
  peak_gr=GenomicRanges::GRanges(sapply(strsplit(rownames(mt_atac2_bulk),":|_|-"),'[',1),
                    IRanges(as.numeric(sapply(strsplit(rownames(mt_atac2_bulk),":|_|-"),'[',2)),
                            as.numeric(sapply(strsplit(rownames(mt_atac2_bulk),":|_|-"),'[',3))))
  
  df_final <- c()
  for(i in 1:length(cnvres)){
    clonei <- names(cnvres)[i]
    segMT_longi <- segMT_long[segMT_long$sampleID==clonei,]
    segs <- as.character(segMT_longi$segName) #binID not the same accross SampleID
    seg_gr <- GenomicRanges::GRanges(sapply(strsplit(segs,":|_|-"),'[',1),
                                     IRanges(as.numeric(sapply(strsplit(segs,":|_|-"),'[',2)),
                                             as.numeric(sapply(strsplit(segs,":|_|-"),'[',3))))
    mcols(seg_gr) <- segMT_longi
    
    hits <- findOverlaps(peak_gr, seg_gr)
    idx_seg <- subjectHits(hits)
    idx_peak <- queryHits(hits)
    peak_gr2 <- peak_gr[idx_peak,]
    peak_gr2$peakID <- rownames(mt_atac2_bulk)[idx_peak]
    mcols(peak_gr2) <- c(mcols(peak_gr2), mcols(seg_gr[idx_seg,]))
    df_peak <- data.frame(mcols(peak_gr2))
    
    rm(peak_gr2,seg_gr,segMT_longi)
    invisible(gc())
    
    mt_atac2_bulk_long_i <- mt_atac2_bulk_long %>%
      filter(sampleID==clonei)%>%
      left_join(df_peak[,c("Chromosome","peakID","segName","ploidy","relativeCN","integerCN","sampleID","Cluster")]) 
    mt_atac2_bulk_long_i <- unique(na.omit(mt_atac2_bulk_long_i))

    mt_atac2_bulk_long_i <- mt_atac2_bulk_long_i %>%
      dplyr::group_by(Cluster,segName) %>%
      dplyr::mutate(Clmean_nCountATAC_seg = mean(Clmean_nCountATAC_peak,na.rm = T))%>%
      as.data.frame()
    segDF <- unique(mt_atac2_bulk_long_i[,c("Cluster","sampleID","Chromosome","segName","ploidy","relativeCN","integerCN","Clmean_nCountATAC_seg")])
    segDF$segStart <- as.numeric(sapply(strsplit(segDF$segName,"_"),"[",2))
    segDF$segEnd <- as.numeric(sapply(strsplit(segDF$segName,"_"),"[",3))
    segDF$segLen <- segDF$segEnd-segDF$segStart
    rm(mt_atac2_bulk_long_i)
    df_final <- rbind(df_final,segDF)
    }

  ##ATAC dataframe
  df_final <- unique(df_final)
  df_final <- df_final[!is.na(df_final$segName),,drop=F]
  
   
  color_c <- pal_ucscgb(alpha = 0.8)(length(sort(unique(df_final$Cluster))))
  names(color_c) <- sort(unique(df_final$Cluster))
  
  if(FiltSegbyLen){
    df_final_filt <- df_final[df_final$segLen>1e6,,drop=F]
  }else{df_final_filt <- df_final}

  df_final_filt$segLen_MB <- df_final_filt$segLen/1e6
  df_final_filt$Clmean_nCountATAC_segnom <- df_final_filt$Clmean_nCountATAC_seg/df_final_filt$segLen_MB
  
  df_final_filt$label <- df_final_filt[,label_term]
  df_final_filt$label[duplicated(df_final_filt$label)] <- NA
  #df_final_filt$label <- as.numeric(sapply(strsplit(df_final_filt$label,"_"),"[",2))
  if(!is.null(label_term)){
    label_term <- "label"
  }

  #Plot for nCount_RNA vs nCount_ATAC
  if(dotPlot4Seg){    
    #integerCN vs nCountATAC
    plot_point(df_final_filt,x_aes = "integerCN",y_aes="Clmean_nCountATAC_segnom",gtitle= paste0("Segment level"),
               x_lab = "integer CN",y_lab = "nCount_ATAC",colour_item = "Cluster",size_item = "segLen_MB",size_lab="length(Mb)",
               color_value = color_c,outfile = paste0(outdir,"/atacIntCNV_nCountATAC_segment.pdf"),
               log_x=F,log_y=F,label_term=label_term)
      
  }
 
  return(df_final_filt)
}





