# Genome-independent plotting helpers for TeaCNV.
.teacnv_plot_genome <- function(genome, auto_install = TRUE) {
  if (inherits(genome, "BSgenome")) return(genome)
  if (!is.character(genome) || length(genome) != 1L ||
      is.na(genome) || !nzchar(trimws(genome)))
    stop("genome must be a genome assembly or BSgenome package name.")
  genome <- trimws(genome)
  packages <- c(hg19 = "BSgenome.Hsapiens.UCSC.hg19",
                hg38 = "BSgenome.Hsapiens.UCSC.hg38",
                mm10 = "BSgenome.Mmusculus.UCSC.mm10",
                mm39 = "BSgenome.Mmusculus.UCSC.mm39")
  key <- tolower(genome)
  pkg <- if (key %in% names(packages)) unname(packages[key]) else
    if (startsWith(genome, "BSgenome.")) genome else NULL
  required <- c("BSgenome", "GenomeInfoDb", pkg)
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    if (!auto_install) stop("Missing packages: ", paste(missing, collapse = ", "))
    if (!requireNamespace("BiocManager", quietly = TRUE))
      utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
    message("Installing genome packages: ", paste(missing, collapse = ", "))
    BiocManager::install(missing, ask = FALSE, update = FALSE)
    missing <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) stop("Unable to load: ", paste(missing, collapse = ", "))
  }
  # For other assemblies a previously installed short name is accepted.
  BSgenome::getBSgenome(if (is.null(pkg)) genome else pkg,
                        masked = FALSE, load.only = TRUE)
}

.teacnv_plot_chr <- function(x) {
  x <- sub("^chr", "", trimws(as.character(x)), ignore.case = TRUE)
  k <- toupper(x) %in% c("X", "Y", "M", "MT")
  x[k] <- toupper(x[k]); x[x %in% "MT"] <- "M"
  x
}

.teacnv_plot_parse <- function(ids) {
  ids <- as.character(ids)
  m <- regmatches(ids, regexec("^(.+)[_:-]([0-9]+)[_:-]([0-9]+)$", ids))
  if (!length(m) || any(lengths(m) != 4L))
    stop("Invalid genomic IDs; expected chr1-0-1000, chr1_0_1000 or chr1:0-1000.")
  z <- data.frame(chrom = .teacnv_plot_chr(vapply(m, `[`, "", 2L)),
                  Start = as.numeric(vapply(m, `[`, "", 3L)),
                  End = as.numeric(vapply(m, `[`, "", 4L)))
  if (anyNA(z$chrom) || any(!nzchar(z$chrom)) ||
      any(!is.finite(z$Start) | !is.finite(z$End) | z$End <= z$Start))
    stop("Invalid genomic coordinates.")
  z
}

.teacnv_plot_num <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  suppressWarnings(as.numeric(x))
}

.teacnv_plot_coordinates <- function(data.bin, genome, true_loc, auto_install) {
  if (is.atomic(data.bin) && is.null(dim(data.bin))) {
    data.bin <- data.frame(value = unname(data.bin), row.names = names(data.bin))
  }
  d <- as.data.frame(data.bin)
  if (!nrow(d)) stop("data.bin has no rows.")
  # Prefer an explicit binID, then genomic row names, then coordinate columns.
  id <- if ("binID" %in% names(d)) d$binID else rownames(d)
  valid_id <- length(id) == nrow(d) && !anyNA(id) &&
    all(grepl("^(.+)[_:-]([0-9]+)[_:-]([0-9]+)$", id))
  if (valid_id) {
    z <- .teacnv_plot_parse(id)
  } else {
    find_col <- function(keys) {
      k <- match(keys, tolower(names(d))); k <- k[!is.na(k)]
      if (length(k)) names(d)[k[1L]] else NA_character_
    }
    cols <- c(find_col(c("chrom", "chr", "chromosome", "seqnames")),
              find_col("start"), find_col("end"))
    if (anyNA(cols)) stop("Provide genomic row names, binID or chromosome/start/end columns.")
    z <- data.frame(chrom = .teacnv_plot_chr(d[[cols[1L]]]),
                    Start = .teacnv_plot_num(d[[cols[2L]]]),
                    End = .teacnv_plot_num(d[[cols[3L]]]))
  }
  if (anyNA(z$chrom) || any(!nzchar(z$chrom)) ||
      any(!is.finite(z$Start) | !is.finite(z$End) | z$Start < 0 | z$End <= z$Start))
    stop("Invalid bin coordinates.")
  used <- unique(z$chrom)
  autos <- used[grepl("^[0-9]+$", used)]
  autos <- autos[order(as.numeric(autos))]
  special <- intersect(c("X", "Y", "M"), used)
  other <- sort(setdiff(used, c(autos, special)))
  chr_order <- c(autos, special, other)
  if (true_loc) {
    g <- .teacnv_plot_genome(genome, auto_install)
    lens <- GenomeInfoDb::seqlengths(g)
    ref <- .teacnv_plot_chr(names(lens))
    k <- match(chr_order, ref)
    if (anyNA(k)) stop("Chromosomes not found in genome: ",
                       paste(chr_order[is.na(k)], collapse = ", "))
    if (any(vapply(chr_order, function(x) sum(ref == x) != 1L, logical(1))))
      stop("Ambiguous chromosome names in the selected genome.")
    len <- as.numeric(lens[k])
    if (any(!is.finite(len) | len <= 0)) stop("Invalid reference chromosome lengths.")
    if (any(z$End > len[match(z$chrom, chr_order)]))
      stop("Bin coordinates exceed reference chromosome lengths; check genome assembly.")
  } else {
    # Equal bin spacing: genome sequence lengths are not needed.
    len <- as.numeric(table(factor(z$chrom, levels = chr_order)))
  }
  chrLine <- data.frame(chrom = chr_order, length = len,
                        chrPosStart_abs = c(0, head(cumsum(len), -1L)))
  chrLine$chrPosEnd_abs <- chrLine$chrPosStart_abs + len
  d$chrom <- z$chrom; d$Start <- z$Start; d$End <- z$End
  d <- d[order(match(d$chrom, chr_order), d$Start, d$End), , drop = FALSE]
  d$chrPosStart_abs <- chrLine$chrPosStart_abs[match(d$chrom, chr_order)]
  d$axis.loc <- if (true_loc) d$chrPosStart_abs + d$Start else seq_len(nrow(d)) - 0.5
  list(data = d, chrLine = chrLine,
       breaks = chrLine$chrPosStart_abs + len/2, labels = chr_order)
}

.teacnv_plot_limits <- function(x, ylim = NULL) {
  if (is.null(ylim)) ylim <- range(x[is.finite(x)])
  if (length(ylim) != 2L || any(!is.finite(ylim)) || ylim[1] > ylim[2])
    stop("Invalid y limits or no finite plotting values.")
  if (ylim[1] == ylim[2]) {
    pad <- max(abs(ylim[1])*0.05, 0.05)
    ylim <- ylim + c(-pad, pad)
  }
  ylim
}

.teacnv_plot_save <- function(p, plotDir, name.data, height, width, device) {
  if (is.null(plotDir)) plotDir <- "."
  if (!dir.exists(plotDir)) dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(file.path(plotDir, paste0(name.data, ".", device)),
                  plot = p, height = height, width = width, dpi = 300, device = device)
}

#'
#' @title seg_plot()
#' @description segment Plot (on absolute genomic position)
#' @param data.bin data.frame of bin level segment with chromosome location as rownames (rowname is important!), such as 'chr1_237513_1530360'
#' @param true_loc plot segments and dots at true genomic location. If not true, plot by order.
#' @export
#' 
# seg_plot = function(data.bin,name.data="",num_cells=NULL,plotDir=NULL,
#                     show_specific_seg=NULL,
#                     specific_seg_color="green",
#                     value.bin="binRatio",
#                     value.segment="SegMean",
#                     ylab="Relative CN",
#                     true_loc=TRUE,
#                     genome="hg38",
#                     plot_seg = TRUE,
#                     seg_col = "#a52a2a",
#                     plot_hist = TRUE,
#                     plot_dots = TRUE,
#                     log_ratio = FALSE,
#                     color_dot =FALSE,
#                     color_hist = TRUE,
#                     color_hist_gradient = TRUE,
#                     plot_colors = NULL,
#                     color_limit = c(0,3),
#                     ylim=NULL,
#                     add_yline = NULL,
#                     label_gene=NULL,
#                     label_gene_color="blue",
#                     annotations=NULL,
#                     height = 3,width=14,
#                     outPlot=TRUE,
#                     device = "png",
#                     color_seg_gradient=FALSE,
#                     disconnected=TRUE
#                     ){
#   suppressMessages({
#     library(BSgenome)
#     library(dplyr)})
#   g <- getBSgenome(genome, masked=FALSE)
#   chrLine <- data.frame(chrom=1:24, length=seqlengths(g)[1:24])
#   #rownames(data.bin)<- data.bin$

#   if(is.null(plotDir)){plotDir="./"}
#   data.bin <- as.data.frame(data.bin)
#   data.bin <- data.bin[,!grepl("seg_SD|seg_score",colnames(data.bin)),drop=F]
#   if(!any(grepl("segID",colnames(data.bin))) & any(grepl("segName",colnames(data.bin))) ){
#     data.bin$segID <- data.bin$segName
#   }
#   if((!any(grepl("chr",rownames(data.bin)))) & any(grepl("binID",colnames(data.bin)))){
#     rownames(data.bin) <- data.bin$binID
#   }
#   chromnum=sapply(strsplit(rownames(data.bin),'-|_|:'),'[',1)
#   chr_name=as.character(na.omit(unique(chromnum)))
#   if(grepl('chr',chromnum[1])){
#     chromnum=sapply(strsplit(chromnum,'hr'),'[',2) ## if chrommnum with chr
#   }else{
#     chromnum=chromnum
#     }
#   chromnum <- gsub("X","23",chromnum)
#   chromnum <- gsub("Y","24",chromnum)
#   data.bin$chrom <- chromnum
#   chr_name <- sapply(strsplit(chr_name,'hr'),'[',2)
#   # chr_name <- na.omit(chr_name)
#   chr_name_nm <- gsub("X","23",chr_name)
#   chr_name_nm <- gsub("Y","24",chr_name_nm)
#   chromnum <- na.omit(chromnum)
 
#   #
#   if(all(!grepl("^Start$",colnames(data.bin)))){
#     data.bin$Start <- as.numeric(sapply(strsplit(rownames(data.bin),'-|_|:'),'[',2))
#   }
#   if(true_loc){
#     chrLine <- chrLine[chrLine$chrom %in% as.character(unique(chromnum)),,drop=F] 
#     genomeLength <- sum(chrLine)
#     if(nrow(chrLine)>1){
#       chrPosStart = c(0,chrLine[1:(nrow(chrLine)-1),2])
#     }else{chrPosStart = 0}

#     chrLine$chrPosStart_abs <- cumsum(chrPosStart)

#     chrLine_mpOrder <- match(data.bin$chrom,chrLine$chrom)
    
#     data.bin <- cbind(data.bin,chrLine[chrLine_mpOrder,,drop=F])
    
#     #data.bin <- merge(data.bin,chrLine,by=merge_by,all.x=T)
#     data.bin <- data.bin[order(as.integer(data.bin$chrom),data.bin$chrPosStart_abs,data.bin$Start),]
#     if(plot_seg){
#        data.bin$seg_start_abs <- as.numeric(sapply(strsplit(data.bin$segName,'-|_|:'),'[',2))+as.numeric(data.bin$chrPosStart_abs)
#        data.bin$seg_end_abs <- as.numeric(sapply(strsplit(data.bin$segName,'-|_|:'),'[',3))+as.numeric(data.bin$chrPosStart_abs)
       
#      }
#     # data.bin <- na.omit(data.bin)
    
#     maploc <- data.bin$chrPosStart_abs+data.bin$Start
#     intercept_x =chrLine$chrPosStart_abs
#     breakss<-cumsum(chrPosStart)+chrLine[1:(nrow(chrLine)),2]*0.5
#     which_chrom <- which(as.character(chrLine$chrom) %in% chr_name_nm)
#     breakss<-breakss[which_chrom]
#     if(plot_seg){
#     a<-data.frame(maploc=maploc,data.bin[,value.bin],data.bin[,value.segment],data.bin$segID,data.bin$seg_start_abs,data.bin$seg_end_abs,data.bin$segName)
#     colnames(a)=c('axis.loc','value.bin','value.segment','segID',"segStart","segEnd",'segName')
#     }else{
#       if(all(!grepl("^End$",colnames(data.bin)))){
#         data.bin$End <- as.numeric(sapply(strsplit(rownames(data.bin),'-|_|:'),'[',3))
#       }
#       a<-data.frame(maploc=maploc,data.bin[,c(value.bin,"Start","End")])
#       colnames(a)=c('axis.loc','value.bin','Start','End')
#       a$length_seg = a$End - a$Start
      
#     }
    
    
#   }else{
   
#     chromnum=chromnum[order(as.integer(gsub("[^0-9]", "", chromnum)))] #order by chr
#     chromnum_xaxis = table(chromnum)[order(as.integer(gsub("[^0-9]", "", names(table(chromnum)))))]#order by chr
#     chrline=cumsum(chromnum_xaxis)
#     maploc=1:nrow(data.bin)
#     breakss<-as.numeric(c(0,chrline[-length(chrline)])+chromnum_xaxis*0.5)
#     intercept_x =chrline
#     data.bin <- data.bin[order(as.integer(data.bin$chrom),data.bin$Start),]
#     a<-data.frame(maploc=maploc,data.bin[,value.bin],data.bin[,value.segment],data.bin$segID,data.bin$length_seg,data.bin$segName)
#     colnames(a)=c('axis.loc','value.bin','value.segment','segID','length_seg','segName')
#     a$rank <- maploc
#     if(!is.factor( a$segID)){
#       a$segID <- factor(a$segID,levels = sort(unique(as.numeric(as.character(a$segID)))))
#     }
#     a <- a %>%
#       dplyr::group_by(segID) %>%
#       dplyr::mutate(segStart=min(rank),segEnd=max(rank))%>%
#       as.data.frame()
#   }

  
#   if(is.null(add_yline)){
#     if(log_ratio){
#       a$value.bin <- log1p(a$value.bin)
#       a$value.segment <- log1p(a$value.segment)
#       add_yline = log1p(1)
#     }else{add_yline = 1} 
#   }

  
#   if(is.null(ylim)){
#     if(plot_dots & (!plot_seg)){
#       ylim=c(quantile(a$value.bin,0.02),quantile(a$value.bin,0.95))
#     }else{
#       cen <- median(a$value.segment,na.rm = T)
#       dist1 <- abs(cen-max(a$value.segment))
#       dist2 <- abs(cen-min(a$value.segment))
#       d <- max(dist1,dist2)*1.2
#       ylim_l <- ifelse((cen-d) >=0,cen-d,0)
#       ylim_h <- cen+d
#       ylim=c(ylim_l,ylim_h)
#     }
   
#   }
 
#   if(!is.null(num_cells)){
#     ptitle = paste(name.data," (",num_cells," cells)")
#   }else{
#     ptitle = paste(name.data)
#   }
  

#   chromnum=sapply(strsplit(rownames(data.bin),'-|_|:'),'[',1)
#   chr_name_label=as.character(na.omit(unique(chromnum)))
#   chr_name_label <- sapply(strsplit(chr_name_label,'hr'),'[',2)
#   chr_name_label <- na.omit(chr_name_label)
#   if(length(chr_name_label)>1){
#     chr_name_label[as.character(chr_name_label)%in%c("19","21")] <- ""
#   }

#     p1 <- ggplot(a) +
#     #geom_point(aes(x = axis.loc, y = value.bin), alpha = 1, shape = 16,size=1) +
#     geom_hline(yintercept =add_yline,col="grey",linetype = 2)+
#     theme_bw() +
#     labs(x = "", y = ylab)+
#     coord_cartesian(ylim = ylim)+
#     ggtitle(ptitle)+
#     #geom_vline(xintercept =intercept_x,col="grey")+
#     theme(plot.margin = unit(c(0.5,0,0,0.5), 'lines'), #c(top, right, bottom, left)
#       panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
#           legend.position = "none",plot.title = element_text(size = 15),
#           axis.title=element_text(size=15),
#           axis.text = element_text(size = 15))+
#     scale_x_continuous(guide = guide_axis(check.overlap = TRUE),breaks=breakss,labels=chr_name_label)
   
#     if(length(intercept_x)>1){

#       odd_positions <- intercept_x[seq(1, length(intercept_x), by = 2)]
     
#       even_positions <- intercept_x[seq(2, length(intercept_x), by = 2)]

#       if(length(odd_positions)!=length(even_positions)){
#         nrow <- min(length(odd_positions),length(even_positions))
#       }else{nrow=length(odd_positions)}
#       ## highlight region data
#       rects <- data.frame(start=odd_positions[1:nrow], end=even_positions[1:nrow])
#       p1 <- p1+
#         geom_rect(data=rects, inherit.aes=FALSE, 
#                   aes(xmin=start, xmax=end, ymin=-Inf,ymax=Inf),
#                   color="transparent", fill="grey50", alpha=0.3)
#       }
#  # brks <- unique(c(seq(0.1,0.6,length.out =2),seq(0.8,1.2,length.out =2),seq(1.5,2,length.out =3),2.5,3))
#     if(color_hist_gradient){
#       if(is.null(plot_colors)){
#         colorss <- c("#023e8a","grey50", "#cb793a", "#9a031e","#6a040f","#370617")
#         #colorss <- c("black","#023e8a","grey50","#fcdc4d", "#cb793a", "#9a031e","#5f0f40")
        
#       }else{
#         colorss <- plot_colors
#       }
#     }else{
#       if(is.null(plot_colors)){
#         colorss <- "#472d30"
#       }else{
#         colorss <- plot_colors[1]
#       }
      
#     }

  
#   if(plot_dots){

#     if(color_dot){

#       p1 <- p1+
#         geom_point(aes(x = axis.loc, y = value.bin,colour=value.segment), alpha = 1, shape = 16,size=0.5)+
#         scale_color_gradientn(name = "",colours = colorss,#c("#003049", "#669bbc","#7e7f83","#c1121f","#780000"),
#                             limits  = color_limit,
#                             oob = scales::squish)
#     }else{
      
#       p1 <- p1+
#         geom_point(aes(x = axis.loc, y = value.bin),colour="darkgrey", alpha = 1, shape = 16,size=0.5)
    
#     }
#   }

  
#   if(plot_seg){
#     if(color_seg_gradient){
#       p1 <- p1+
#         geom_segment(
#           data = a, 
#           mapping = aes(x=segStart, y=value.segment, xend=segEnd, yend=value.segment,colour=value.segment), 
#           na.rm =T,size = 1,
#           inherit.aes = FALSE)+
#         scale_color_gradientn(name = "",colours = colorss,#c("#003049", "#669bbc","#7e7f83","#c1121f","#780000"),
#                               limits  = color_limit,
#                               oob = scales::squish)
#     }else{ 
#     p1 <- p1+
#       #geom_line(aes(x = maploc,y=value.segment),col="red")+
      
#       geom_segment(
#         data = a, 
#         mapping = aes(x=segStart, y=value.segment, xend=segEnd, yend=value.segment), 
#         colour=seg_col,na.rm =T,size = 1,
#         inherit.aes = FALSE)
#     }
    
#     if(!disconnected){
#       dat2<- data.frame(x=a$segEnd[1:(length(a$segEnd)-1)], y=a$value.segment[1:(length(a$segEnd)-1)],
#                               xend=a$segStart[2:length(a$segEnd)], yend=a$value.segment[2:length(a$segEnd)])
#       p1 <- p1 +
#         geom_segment(
#           data = dat2,
#           mapping = aes(x=x, y=y,
#                         xend=xend, yend=yend), 
#           colour=seg_col,na.rm =T,size = 1,
#           inherit.aes = FALSE)
#     }
#   }


#   if(!is.null(show_specific_seg)){ 
#     # a$axis_seg_specific = NA
#     # whichRow = which(a$segName %in% show_specific_seg)
#     # a$axis_seg_specific[whichRow] <- maploc[whichRow] 
#     # line_df.seg <- a[complete.cases(a), ]

#     show_seg_chr <- sapply(strsplit(show_specific_seg,'-|_|:'),'[',1)
#     show_seg_chr <- gsub("chr","",show_seg_chr)
#     show_seg_chr <- gsub("X","23",show_seg_chr)
#     show_seg_chr <- gsub("Y","24",show_seg_chr)
#     show_seg_start <- as.numeric(sapply(strsplit(show_specific_seg,'-|_|:'),'[',2))
#     show_seg_end <- as.numeric(sapply(strsplit(show_specific_seg,'-|_|:'),'[',3))
#     df_show <- data.frame(Chromosome=show_seg_chr,cnv_start=show_seg_start,cnv_end=show_seg_end)
#     line_df.seg <- cbind(data.bin[,c("chrom","Start","End","chrPosStart_abs")],a) %>%
#       dplyr::inner_join(df_show%>% dplyr::select(Chromosome, cnv_start, cnv_end), 
#         join_by(chrom == Chromosome,
#           Start <= cnv_end,
#           End >= cnv_start)) %>%
#       dplyr::mutate(
#         seg_show = paste(chrom, cnv_start, cnv_end, sep = "_"),
#         seg_show_start = cnv_start+chrPosStart_abs,
#         seg_show_end = cnv_end +chrPosStart_abs
#       ) %>%
#       dplyr::select(-cnv_start, -cnv_end)
        
#     p1=p1+geom_segment(
#       data = line_df.seg, 
#       mapping = aes(x=seg_show_start, y=value.segment, xend=seg_show_end, yend=value.segment), 
#       col=specific_seg_color,na.rm =T,size = 1.5,
#       inherit.aes = FALSE)
#   }
#   if(!is.null(label_gene)){
#     gene_loc <- annotate_gene(label_gene,genome=genome,annotations=annotations)
#     gene_chr <- as.character(seqnames(gene_loc))
#     gene_chr <- sapply(strsplit(gene_chr,'hr'),'[',2)
#     if(true_loc){
#     label_x <- chrLine[chrLine$chrom==gene_chr,"chrPosStart_abs"]+start(gene_loc)
#     label_y <- unique(data.bin[data.bin$chrom==gene_chr&data.bin$seg_start_abs<=label_x&data.bin$seg_end_abs>=label_x,value.segment])[1]
#     }else{
#       g_start <- start(ranges(gene_loc))
#       g_end <- end(ranges(gene_loc))
#       label_x <- which(data.bin$chrom==gene_chr&data.bin$Start>=g_start)[1]
#       label_y <- data.bin[data.bin$chrom==gene_chr&data.bin$Start>=g_start,value.segment][1]
      
#     }

#     p1=p1+
#       geom_point(aes(x = label_x, y = label_y), alpha = 1, shape = 17,size=1,colour=label_gene_color)+
#       geom_text(aes(label_x,label_y,label = label_gene),color=label_gene_color,vjust = -1)+
#       geom_segment(aes(x = label_x, y = label_y, xend = label_x, yend = label_y -0.5),
#                    color = label_gene_color,  # 指定指引线颜色
#                    linetype = "dashed")  
    
#   }
  
#   if(plot_hist){
    
#   breakshist <- seq(min(ylim),max(ylim),length.out=80)
#   #segSet = unique(data.bin[,c("binID","segID",value.segment,"segName"),drop=F]) ###
#   if(true_loc){
#     if(plot_seg){
#       a$segLen <- as.numeric(a$segEnd)-as.numeric(a$segStart)
#     }else{
#      a$segLen <- as.numeric(a$length_seg) #Here, length_seg is the length of peak
#     }
#   }else{
#     a$segLen <- as.numeric(a$length_seg) 
#   }  
  
#   if(!plot_seg){
#     a$segID <- seq(1:nrow(a))
#   }
#   segLen_sum <- sum(unique(a[,c("segID","segLen")]$segLen))
#   a$wei <- a$segLen/segLen_sum
#   if(plot_seg){
#     p_x <- "value.segment"
#   }else{
#     p_x <- "value.bin"
#   }


#   if(color_hist){
#     if(color_hist_gradient){
#       p3 <- ggplot(a,aes(x=get(p_x),fill = ..x..))+
#         geom_histogram(aes(y=..density.., weight = wei),breaks=breakshist)+ #color="#E69F00",fill="wheat",
#         labs(x = "", y = "",fill="")+
#         ggtitle("Summary histogram")+
#         coord_cartesian(xlim=ylim)+
#         theme_bw() +
#         theme_classic()+
#         theme(plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
#               # legend.position = "bottom",
#               plot.title = element_text(size = 15),
#               axis.title=element_text(size=15),
#               axis.text = element_text(size = 12))+
#         coord_flip(xlim=ylim)+
#       scale_fill_gradientn(name = "",colours = colorss,#c("#003049", "#669bbc","#7e7f83","#c1121f","#780000"),
#                            limits  = color_limit,
#                            oob = scales::squish)
#     }else{
#       p3=ggplot(a,aes(x=get(p_x)))+
#         geom_histogram(aes(y=..density.., weight = wei),color=colorss,fill=colorss,breaks=breakshist)+ #color="#E69F00",fill="wheat",
#         labs(x = "", y = "",fill="")+
#         ggtitle("Summary histogram")+
#         coord_cartesian(xlim=ylim)+
#         theme_bw() +
#         theme_classic()+
#         theme(plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
#               # legend.position = "bottom",
#               plot.title = element_text(size = 15),
#               axis.title=element_text(size=15),
#               axis.text = element_text(size = 12))+
#         coord_flip(xlim=ylim)
#   }
#   }else{
#     p3=ggplot(a,aes(x=get(p_x)))+
#       geom_histogram(aes(y=..density.., weight = wei),breaks=breakshist)+ #color="#E69F00",fill="wheat",
#       labs(x = "", y = "",fill="")+
#       ggtitle("Summary histogram")+
#       coord_cartesian(xlim=ylim)+
#       theme_bw() +
#       theme_classic()+
#       theme(plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
#             # legend.position = "bottom",
#             plot.title = element_text(size = 15),
#             axis.title=element_text(size=15),
#             axis.text = element_text(size = 12))+
#       coord_flip(xlim=ylim)
#   }
#   pp=ggpubr::ggarrange(p1, p3, align ="h",ncol = 2, nrow = 1,widths = c(15,4),heights=3 )
#   p_ls <- list(ggarranged_p=pp,p1=p1,p2=p3)
#   }else{
#     pp <- p1
#   p_ls <- list(p1=p1)
#   }
#   if(outPlot){
#     ggsave(paste0(plotDir,"/",name.data,".",device),pp,height = height,width=width,dpi = 300,device = device)
#   }
#   return(p_ls)
# }
seg_plot <- function(data.bin, name.data = "", num_cells = NULL, plotDir = NULL,
                     show_specific_seg = NULL, specific_seg_color = "green",
                     value.bin = "binRatio", value.segment = "SegMean",
                     ylab = "Relative CN", true_loc = TRUE, genome = "hg38",
                     plot_seg = TRUE, seg_col = "#a52a2a", plot_hist = TRUE,
                     plot_dots = TRUE, log_ratio = FALSE, color_dot = FALSE,
                     color_hist = TRUE, color_hist_gradient = TRUE,
                     plot_colors = NULL, color_limit = c(0,3), ylim = NULL,
                     add_yline = NULL, label_gene = NULL, label_gene_color = "blue",
                     annotations = NULL, height = 3, width = 14, outPlot = TRUE,
                     device = "png", color_seg_gradient = FALSE,
                     disconnected = TRUE, auto_install = TRUE) {
  pdat <- .teacnv_plot_coordinates(data.bin, genome, true_loc, auto_install)
  d <- pdat$data; chrLine <- pdat$chrLine
  if (!value.bin %in% names(d)) stop("Missing bin value column: ", value.bin)
  a <- data.frame(axis.loc = d$axis.loc, chrom = d$chrom,
                  Start = d$Start, End = d$End,
                  value.bin = .teacnv_plot_num(d[[value.bin]]))
  if (plot_seg) {
    if (!all(c(value.segment, "segName") %in% names(d)))
      stop("plot_seg=TRUE requires segName and ", value.segment, ".")
    s <- .teacnv_plot_parse(d$segName)
    if (any(s$chrom != d$chrom)) stop("segName chromosome does not match its bin.")
    if (true_loc && any(s$End > chrLine$length[match(d$chrom, chrLine$chrom)]))
      stop("Segment exceeds reference chromosome length.")
    a$value.segment <- .teacnv_plot_num(d[[value.segment]])
    a$segName <- paste(s$chrom, s$Start, s$End, sep = "_")
    a$segID <- match(a$segName, unique(a$segName))
    a$segStart <- if (true_loc) s$Start + d$chrPosStart_abs else
      ave(a$axis.loc, a$segID, FUN = min) - 0.5
    a$segEnd <- if (true_loc) s$End + d$chrPosStart_abs else
      ave(a$axis.loc, a$segID, FUN = max) + 0.5
    a$segStart_local <- s$Start; a$segEnd_local <- s$End
  } else {
    a$value.segment <- a$value.bin
    a$segID <- seq_len(nrow(a))
    a$segName <- paste(d$chrom, d$Start, d$End, sep = "_")
    a$segStart_local <- a$Start; a$segEnd_local <- a$End
    a$segStart <- if (true_loc) d$chrPosStart_abs + a$Start else a$axis.loc - 0.5
    a$segEnd <- if (true_loc) d$chrPosStart_abs + a$End else a$axis.loc + 0.5
  }
  if (log_ratio) {
    a$value.bin <- suppressWarnings(log1p(a$value.bin))
    a$value.segment <- suppressWarnings(log1p(a$value.segment))
  }
  if (is.null(add_yline)) add_yline <- if (log_ratio) log1p(1) else 1
  v <- if (plot_seg) a$value.segment else a$value.bin
  v <- v[is.finite(v)]
  if (!length(v)) stop("No finite values available for plotting.")
  if (is.null(ylim)) {
    if (plot_dots && !plot_seg) ylim <- as.numeric(stats::quantile(v, c(.02,.95))) else {
      cen <- stats::median(v); delta <- max(abs(v-cen))*1.2
      ylim <- c(if (all(v >= 0)) max(0, cen-delta) else cen-delta, cen+delta)
    }
  }
  ylim <- .teacnv_plot_limits(v, ylim)
  ptitle <- if (is.null(num_cells)) name.data else paste0(name.data, " (", num_cells, " cells)")
  colorss <- if (!is.null(plot_colors)) plot_colors else if (color_hist_gradient)
    c("#023e8a", "grey50", "#cb793a", "#9a031e", "#6a040f", "#370617") else "#472d30"
  gradient <- if (length(colorss) >= 2L) colorss else rep(colorss, 2L)
  odd <- seq.int(1L, nrow(chrLine), by = 2L)
  rects <- data.frame(start = chrLine$chrPosStart_abs[odd], end = chrLine$chrPosEnd_abs[odd])
  p1 <- ggplot2::ggplot(a) + ggplot2::theme_bw() +
    ggplot2::geom_rect(data = rects, inherit.aes = FALSE,
                       ggplot2::aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                       fill = "grey50", alpha = .3, colour = NA) +
    ggplot2::geom_hline(yintercept = add_yline, colour = "grey", linetype = 2) +
    ggplot2::labs(x = "", y = ylab, title = ptitle) +
    ggplot2::coord_cartesian(ylim = ylim, xlim = c(0, tail(chrLine$chrPosEnd_abs,1))) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "none",
                   plot.title = ggplot2::element_text(size = 15),
                   axis.title = ggplot2::element_text(size = 15),
                   axis.text = ggplot2::element_text(size = 15)) +
    ggplot2::scale_x_continuous(breaks = pdat$breaks, labels = pdat$labels,
                                guide = ggplot2::guide_axis(check.overlap = TRUE))
  if (plot_dots) {
    if (color_dot) p1 <- p1 + ggplot2::geom_point(ggplot2::aes(axis.loc, value.bin, colour = value.segment), size = .5, na.rm = TRUE)
    else p1 <- p1 + ggplot2::geom_point(ggplot2::aes(axis.loc, value.bin), colour = "darkgrey", size = .5, na.rm = TRUE)
  }
  segs <- unique(a[,c("chrom", "segID", "segStart", "segEnd", "segStart_local", "segEnd_local", "value.segment")])
  if (plot_seg) {
    if (color_seg_gradient) p1 <- p1 + ggplot2::geom_segment(data = segs, inherit.aes = FALSE,
      ggplot2::aes(x = segStart, xend = segEnd, y = value.segment, yend = value.segment, colour = value.segment), linewidth = 1, na.rm = TRUE)
    else p1 <- p1 + ggplot2::geom_segment(data = segs, inherit.aes = FALSE,
      ggplot2::aes(x = segStart, xend = segEnd, y = value.segment, yend = value.segment), colour = seg_col, linewidth = 1, na.rm = TRUE)
    if (!disconnected && nrow(segs) > 1L) {
      segs <- segs[order(match(segs$chrom, chrLine$chrom), segs$segStart), , drop = FALSE]
      i <- seq_len(nrow(segs)-1L); i <- i[segs$chrom[i] == segs$chrom[i+1L]]
      if (length(i)) {
        links <- data.frame(x=segs$segEnd[i], xend=segs$segStart[i+1L],
                            y=segs$value.segment[i], yend=segs$value.segment[i+1L])
        p1 <- p1 + ggplot2::geom_segment(data=links, inherit.aes=FALSE,
          ggplot2::aes(x=x, xend=xend, y=y, yend=yend), colour=seg_col, linewidth=1, na.rm=TRUE)
      }
    }
  }
  if ((plot_dots && color_dot) || (plot_seg && color_seg_gradient))
    p1 <- p1 + ggplot2::scale_color_gradientn(name="", colours=gradient, limits=color_limit, oob=scales::squish)
  if (!is.null(show_specific_seg)) {
    q <- .teacnv_plot_parse(show_specific_seg)
    highlights <- lapply(seq_len(nrow(q)), function(i) {
      if (true_loc) {
        j <- which(segs$chrom == q$chrom[i] & segs$segStart_local < q$End[i] & segs$segEnd_local > q$Start[i])
        if (!length(j)) return(NULL)
        offset <- chrLine$chrPosStart_abs[match(q$chrom[i], chrLine$chrom)]
        data.frame(x=pmax(segs$segStart_local[j], q$Start[i])+offset,
                   xend=pmin(segs$segEnd_local[j], q$End[i])+offset, y=segs$value.segment[j])
      } else {
        j <- which(a$chrom == q$chrom[i] & a$Start < q$End[i] & a$End > q$Start[i])
        if (!length(j)) return(NULL)
        do.call(rbind, lapply(split(j, a$segID[j]), function(k)
          data.frame(x=min(a$axis.loc[k])-.5, xend=max(a$axis.loc[k])+.5, y=a$value.segment[k[1L]])))
      }
    })
    highlights <- do.call(rbind, highlights)
    if (NROW(highlights)) p1 <- p1 + ggplot2::geom_segment(data=highlights, inherit.aes=FALSE,
      ggplot2::aes(x=x,xend=xend,y=y,yend=y), colour=specific_seg_color, linewidth=1.5, na.rm=TRUE)
  }
  if (!is.null(label_gene)) {
    annotate <- get0("annotate_gene", mode="function", inherits=TRUE)
    if (is.null(annotate) && requireNamespace("TeaCNV", quietly=TRUE))
      annotate <- get0("annotate_gene", envir=asNamespace("TeaCNV"), inherits=FALSE)
    for (gene in label_gene) {
      pos <- tryCatch({
        if (is.null(annotate)) stop("annotate_gene() is not available")
        gr <- annotate(gene, genome=genome, annotations=annotations)
        if (!length(gr)) stop("gene not found")
        gc <- .teacnv_plot_chr(as.character(GenomeInfoDb::seqnames(gr))[1L])
        gs <- BiocGenerics::start(gr)[1L]
        j <- which(a$chrom == gc & a$segStart_local <= gs & a$segEnd_local >= gs & is.finite(a$value.segment))
        if (!length(j)) stop("gene does not overlap plotted segments")
        if (true_loc) x <- chrLine$chrPosStart_abs[match(gc, chrLine$chrom)] + gs else {
          k <- j[which.min(abs(a$Start[j]-gs))]; x <- a$axis.loc[k]
        }
        data.frame(x=x,y=a$value.segment[j[1L]],label=gene)
      }, error=function(e) {warning("Skipping gene ", gene, ": ", conditionMessage(e), call.=FALSE); NULL})
      if (!is.null(pos)) p1 <- p1 + ggplot2::geom_point(data=pos, inherit.aes=FALSE,
        ggplot2::aes(x=x,y=y), shape=17, colour=label_gene_color) +
        ggplot2::geom_text(data=pos, inherit.aes=FALSE, ggplot2::aes(x=x,y=y,label=label), vjust=-1, colour=label_gene_color) +
        ggplot2::geom_segment(data=pos, inherit.aes=FALSE, ggplot2::aes(x=x,xend=x,y=y,yend=y-.5), linetype="dashed", colour=label_gene_color)
    }
  }
  if (plot_hist) {
    # Each observed bin contributes its own genomic width, not the entire segment length.
    h <- data.frame(value=if (plot_seg) a$value.segment else a$value.bin, weight=a$End-a$Start)
    h <- h[is.finite(h$value) & is.finite(h$weight) & h$weight>0, , drop=FALSE]
    # Filter to the chosen display limits so histogram breaks always cover observations.
    h <- h[h$value>=ylim[1L] & h$value<=ylim[2L], , drop=FALSE]
    breaks <- seq(ylim[1L], ylim[2L], length.out=80L)
    p3 <- ggplot2::ggplot(h, ggplot2::aes(x=value, weight=weight))
    if (color_hist && color_hist_gradient) p3 <- p3 + ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(density), fill=ggplot2::after_stat(x)), breaks=breaks, na.rm=TRUE) +
      ggplot2::scale_fill_gradientn(name="", colours=gradient, limits=color_limit, oob=scales::squish)
    else if (color_hist) p3 <- p3 + ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(density)), breaks=breaks, fill=colorss[1L], colour=colorss[1L], na.rm=TRUE)
    else p3 <- p3 + ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(density)), breaks=breaks, na.rm=TRUE)
    p3 <- p3 + ggplot2::labs(x="",y="",title="Summary histogram",fill="") +
      ggplot2::theme_classic() + ggplot2::coord_flip(xlim=ylim) +
      ggplot2::theme(plot.title=ggplot2::element_text(size=15),axis.text=ggplot2::element_text(size=12))
    pp <- ggpubr::ggarrange(p1,p3,align="h",ncol=2,nrow=1,widths=c(15,4),heights=3)
    p_ls <- list(ggarranged_p=pp,p1=p1,p2=p3)
  } else { pp <- p1; p_ls <- list(p1=p1) }
  if (outPlot) .teacnv_plot_save(pp,plotDir,name.data,height,width,device)
  return(p_ls)
}



#' @title seg_plot_dot()
#' @export

# seg_plot_dot <- function(data.bin,name.data="",plotDir=NULL,ylab="Count",genome="hg38",outplot=TRUE){
#   suppressMessages(require(BSgenome))
#   g <- getBSgenome(genome, masked=FALSE)
#   chrLine <- data.frame(chrom=1:24, length=seqlengths(g)[1:24])

  
#   if(is.null(plotDir)){plotDir="./"}
#   #order data.bin by chromosome
#   row <- data.frame(rn=rownames(data.bin))
  
#   chromnum=sapply(strsplit(rownames(data.bin),'-|_|:'),'[',1)
#   row$seqnames <- chromnum
#   row$start <- sapply(strsplit(rownames(data.bin),'-|_|:'),'[',2)
#   row$end <-  sapply(strsplit(rownames(data.bin),'-|_|:'),'[',3)

    
#   if(grepl('chr',chromnum[1])){
#     chromnum=sapply(strsplit(chromnum,'hr'),'[',2) ## if chrommnum with chr
#   }else{
#     chromnum=chromnum
#   }
#   chr_name=as.character(unique(chromnum))
#   chr_name= na.omit(chr_name)
  
#   chromnum <- gsub("X","23",chromnum)
#   chromnum <- gsub("Y","24",chromnum)
#   row$chrom <- chromnum
#   data.bin$chrom <- chromnum
  
#   chrLine <- chrLine[chrLine$chrom %in% as.character(unique(chromnum)),,drop=F] 
#   genomeLength <- sum(chrLine)
#   if(nrow(chrLine)>1){
#     chrPosStart = c(0,chrLine[1:(nrow(chrLine)-1),2])
#   }else(
#     chrPosStart = 0
#   )
  
#   chrLine$chrPosStart_abs <- cumsum(chrPosStart)
  
  
#   row <- merge(row,chrLine,by="chrom",all.x=T)
#   row <- row[order(as.integer(row$chrom),as.numeric(row$chrPosStart_abs)),]
#   data.bin <- as.data.frame(data.bin)[row$rn,,drop=F]
  
#   chrLine_mpOrder <- match(data.bin$chrom,chrLine$chrom)
#   data.bin <- cbind(data.bin,chrLine[chrLine_mpOrder,,drop=F])

#   ## scatter plot
#   chr_name_nm <- gsub("X","23",chr_name)
#   chr_name_nm <- gsub("Y","24",chr_name_nm)
#   which_chrom <- which(as.character(chrLine$chrom) %in% chr_name_nm)
#   data.bin <- na.omit(data.bin)
#   if(all(!grepl("^Start$",colnames(data.bin)))){
#     data.bin$Start <- as.numeric(sapply(strsplit(rownames(data.bin),'-|_|:'),'[',2))
#   }
#   maploc <- data.bin$chrPosStart_abs+data.bin$Start
#   intercept_x =chrLine$chrPosStart_abs
#   breakss<-cumsum(chrPosStart)+chrLine[1:(nrow(chrLine)),2]*0.5
#   which_chrom <- which(as.character(chrLine$chrom) %in% chr_name_nm)
#   breakss<-breakss[which_chrom]
  
#   chr_name_label <- chr_name
#   chr_name_label[as.character(chr_name_label)%in%c("19","21")] <- ""
  
  
#   a<-data.frame(maploc=maploc,data.bin)
#   colnames(a)=c('axis.loc','value.bin')
#   ylim=c(min(a$value.bin),max(a$value.bin))
#   p1 <- ggplot(a) +
#     geom_point(aes(x = axis.loc, y = value.bin), alpha = 1, shape = 16,size=0.6) +
#     theme_bw() +
#     labs(x = "", y = ylab)+
#     ylim(ylim)+
#     ggtitle(paste(name.data))+
#     #geom_line(aes(x = maploc,y=value.segment),col="red")+
#     geom_vline(xintercept =intercept_x,col="grey")+
#     #geom_hline(yintercept =1,col="grey",linetype = 2)+
#     theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
#           legend.position = "bottom",plot.title = element_text(size = 15),
#           axis.title=element_text(size=15),
#           axis.text = element_text(size = 15))+
#     scale_x_continuous(breaks=breakss,labels=chr_name_label)
 
#   breakshist <- seq(min(ylim),max(ylim),length.out=100)
#   a$segLen <- as.numeric(a$segEnd)-as.numeric(a$segStart)
#   segLen_sum <- sum(unique(a[,c("segID","segLen")]$segLen))
#   a$wei <- a$segLen/segLen_sum

#   p3=ggplot(a,aes(x=value.bin))+
#     geom_histogram(aes(y=..density.., weight = wei),color="#E69F00",fill="wheat",breaks=breakshist)+
#     labs(x = "", y = "Density")+
#     ggtitle("Summary histogram")+
#     coord_cartesian(xlim=ylim)+
#     theme_bw() +
#     theme_classic()+
#     theme(legend.position = "bottom",plot.title = element_text(size = 15),
#           axis.title=element_text(size=15),
#           axis.text = element_text(size = 12))+
#     coord_flip()
#   pp=ggarrange(p1, p3, align ="h",ncol = 2, nrow = 1,widths = c(15,3),heights=3 )
#   if(outplot){
#     ggsave(paste0(plotDir,"/",name.data,".png"),pp,height = 3,width=14,dpi = 300)
#   }else{
#     return(pp)
#   }
 
# }
seg_plot_dot <- function(data.bin, name.data = "", plotDir = NULL,
                         ylab = "Count", genome = "hg38", outplot = TRUE,
                         auto_install = TRUE, value.bin = NULL) {
  # Preserve original signal column names before adding coordinate metadata.
  raw <- if (is.atomic(data.bin) && is.null(dim(data.bin)))
    data.frame(value=unname(data.bin),row.names=names(data.bin)) else as.data.frame(data.bin)
  if (is.null(value.bin)) {
    candidate <- setdiff(names(raw), names(raw)[tolower(names(raw)) %in%
      c("chrom","chr","chromosome","seqnames","start","end","binid","segid","segname","length","length_bin","length_seg","chrposstart_abs","chrposend_abs","axis.loc")])
    if (length(candidate) == 1L) value.bin <- candidate else {
      numeric_cols <- candidate[vapply(raw[candidate],is.numeric,logical(1))]
      if (length(numeric_cols) != 1L) stop("Specify value.bin to identify the signal column.")
      value.bin <- numeric_cols[1L]
    }
  }
  if (!is.character(value.bin) || length(value.bin)!=1L || !value.bin %in% names(raw))
    stop("Missing or invalid value.bin column.")
  pdat <- .teacnv_plot_coordinates(raw, genome, TRUE, auto_install)
  d <- pdat$data; chrLine <- pdat$chrLine
  a <- data.frame(axis.loc=d$axis.loc,value.bin=.teacnv_plot_num(d[[value.bin]]),weight=d$End-d$Start)
  a <- a[is.finite(a$value.bin) & is.finite(a$axis.loc), , drop=FALSE]
  if (!nrow(a)) stop("No finite bin values available.")
  ylim <- .teacnv_plot_limits(a$value.bin)
  p1 <- ggplot2::ggplot(a,ggplot2::aes(x=axis.loc,y=value.bin)) +
    ggplot2::geom_point(alpha=1,shape=16,size=.6) + ggplot2::theme_bw() +
    ggplot2::labs(x="",y=ylab,title=name.data) +
    ggplot2::geom_vline(xintercept=chrLine$chrPosStart_abs,colour="grey") +
    ggplot2::coord_cartesian(ylim=ylim,xlim=c(0,tail(chrLine$chrPosEnd_abs,1L))) +
    ggplot2::scale_x_continuous(breaks=pdat$breaks,labels=pdat$labels,
                               guide=ggplot2::guide_axis(check.overlap=TRUE)) +
    ggplot2::theme(panel.grid=ggplot2::element_blank(),legend.position="bottom",
                   plot.title=ggplot2::element_text(size=15),axis.title=ggplot2::element_text(size=15),
                   axis.text=ggplot2::element_text(size=15))
  p3 <- ggplot2::ggplot(a,ggplot2::aes(x=value.bin,weight=weight)) +
    ggplot2::geom_histogram(ggplot2::aes(y=ggplot2::after_stat(density)),colour="#E69F00",fill="wheat",
                           breaks=seq(ylim[1L],ylim[2L],length.out=100L)) +
    ggplot2::labs(x="",y="Density",title="Summary histogram") +
    ggplot2::theme_classic() + ggplot2::coord_flip(xlim=ylim) +
    ggplot2::theme(plot.title=ggplot2::element_text(size=15),axis.text=ggplot2::element_text(size=12))
  pp <- ggpubr::ggarrange(p1,p3,align="h",ncol=2,nrow=1,widths=c(15,3),heights=3)
  if (outplot) { .teacnv_plot_save(pp,plotDir,name.data,3,14,"png"); return(invisible(pp)) }
  return(pp)
}


# Genome-aware gene annotation for TeaCNV.
# No BSgenome sequence or Signac installation is needed for annotation alone.

.teacnv_annotation_assembly <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- sub('^bsgenome\\.[^.]+\\.[^.]+\\.', '', x)
  x <- sub('\\.masked$', '', x)
  x[grepl('^grch37(\\.p[0-9]+)?$', x)] <- 'hg19'
  x[grepl('^grch38(\\.p[0-9]+)?$', x)] <- 'hg38'
  x[grepl('^grcm38(\\.p[0-9]+)?$', x)] <- 'mm10'
  x[grepl('^grcm39(\\.p[0-9]+)?$', x)] <- 'mm39'
  x
}

.teacnv_annotation_require <- function(pkgs, auto_install) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) return(invisible(NULL))
  if (!isTRUE(auto_install))
    stop('Missing required packages: ', paste(missing, collapse = ', '), call. = FALSE)
  if (!requireNamespace('BiocManager', quietly = TRUE))
    utils::install.packages('BiocManager', repos = 'https://cloud.r-project.org')
  message('Installing required packages: ', paste(missing, collapse = ', '))
  BiocManager::install(missing, ask = FALSE, update = FALSE)
  missing <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop('Unable to load: ', paste(missing, collapse = ', '), call. = FALSE)
  invisible(NULL)
}

# In-memory cache of database handles; AnnotationHub also caches downloads on disk.
.teacnv_annotation_db <- local({
  cache <- new.env(parent = emptyenv())
  function(key, auto_install, ensembl_release) {
    pkgs <- c(hg19 = 'EnsDb.Hsapiens.v75', hg38 = 'EnsDb.Hsapiens.v86',
              mm10 = 'EnsDb.Mmusculus.v79')
    id <- if (key == 'mm39') paste0(key, '_', ensembl_release) else key
    if (exists(id, envir = cache, inherits = FALSE))
      return(get(id, envir = cache, inherits = FALSE))
    if (key %in% names(pkgs)) {
      pkg <- unname(pkgs[key])
      .teacnv_annotation_require(c('ensembldb', 'AnnotationFilter', pkg), auto_install)
      db <- getExportedValue(pkg, pkg)
    } else if (key == 'mm39') {
      .teacnv_annotation_require(c('ensembldb', 'AnnotationFilter', 'AnnotationHub'), auto_install)
      ah <- AnnotationHub::AnnotationHub()
      hits <- AnnotationHub::query(ah, c('EnsDb', 'Mus musculus'))
      md <- as.data.frame(S4Vectors::mcols(hits))
      if (!all(c('genome', 'title', 'rdataclass', 'species') %in% names(md)))
        stop('AnnotationHub metadata lacks assembly/title/class fields.', call. = FALSE)
      release_pattern <- paste0('(^|[^0-9])', ensembl_release, '([^0-9]|$)')
      idx <- which(.teacnv_annotation_assembly(md$genome) == 'mm39' &
                   md$rdataclass == 'EnsDb' & tolower(md$species) == 'mus musculus' &
                   grepl(release_pattern, md$title))
      if (length(idx) != 1L)
        stop('Expected one mouse GRCm39 EnsDb for Ensembl ', ensembl_release,
             '; found ', length(idx), '. Check AnnotationHub or supply annotations.', call. = FALSE)
      message('Using AnnotationHub: ', names(hits)[idx], ' (', md$title[idx], ')')
      db <- hits[[names(hits)[idx]]]
    } else {
      stop('No built-in annotation for ', key, '; provide annotations.', call. = FALSE)
    }
    actual <- unique(na.omit(as.character(GenomeInfoDb::genome(GenomeInfoDb::seqinfo(db)))))
    if (!length(actual) || !all(.teacnv_annotation_assembly(actual) == key))
      stop('Annotation database assembly mismatch: expected ', key, '; found ',
           paste(actual, collapse = ', '), call. = FALSE)
    assign(id, db, envir = cache)
    db
  }
})


# #' @title annotate_gene()
# #' @description genomic annotation for gene symbol
# #' @param gene_name gene symbol
# #' @param genome string of reference genome
# #' @param annotations GRanges object
# #' @return a GRanges object of gene annotation
# #' @export
# #' 
# annotate_gene <- function(gene_name,genome="hg38",annotations=NULL){
#   suppressMessages({
#     library(BSgenome)
#     library(EnsDb.Hsapiens.v86)})
#   if(is.null(annotations)&genome%in%c("hg38","GRCh38")){
#     annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#   }
#   anno <- annotations[annotations$gene_name%in%gene_name,,drop=F]
#   gene_loc <- CollapseToLongestTranscript(anno)
#   return(gene_loc)
# }

#' Locate genes using assembly-matched annotations
#' @param gene_name Character vector of exact gene symbols.
#' @param genome Assembly name, e.g. hg19, hg38, mm10, mm39.
#' @param annotations Optional GRanges with gene_name (or gene_symbol) metadata,
#'   or an EnsDb object. Supplied annotations take precedence.
#' @param auto_install Install missing dependencies using BiocManager.
#' @param ensembl_release Ensembl release for mm39 via AnnotationHub (default 109).
#' @return GRanges, one longest transcript span per gene locus where available.
#' @export
annotate_gene <- function(gene_name, genome = 'hg38', annotations = NULL,
                          auto_install = TRUE, ensembl_release = 109L) {
  if (!is.logical(auto_install) || length(auto_install) != 1L || is.na(auto_install))
    stop('auto_install must be TRUE or FALSE.', call. = FALSE)
  if (inherits(genome, 'BSgenome')) {
    ids <- unique(na.omit(as.character(GenomeInfoDb::genome(genome))))
    if (length(ids) != 1L) stop('Cannot determine one assembly from BSgenome.')
    genome <- ids
  }
  if (!is.character(genome) || length(genome) != 1L || is.na(genome) || !nzchar(trimws(genome)))
    stop('genome must be a reference assembly name.', call. = FALSE)
  key <- .teacnv_annotation_assembly(genome)
  if (!is.numeric(ensembl_release) || length(ensembl_release) != 1L ||
      !is.finite(ensembl_release) || ensembl_release <= 0 || ensembl_release != floor(ensembl_release))
    stop('ensembl_release must be a positive integer.', call. = FALSE)
  gene_name <- unique(trimws(as.character(gene_name)))
  gene_name <- gene_name[!is.na(gene_name) & nzchar(gene_name)]
  .teacnv_annotation_require(c('GenomicRanges', 'GenomeInfoDb', 'IRanges', 'S4Vectors', 'BiocGenerics'), auto_install)
  if (!length(gene_name)) return(GenomicRanges::GRanges())
  if (is.null(annotations)) annotations <- .teacnv_annotation_db(key, auto_install, ensembl_release)
  # Validate supplied assembly metadata; never silently relabel another assembly.
  if (inherits(annotations, 'EnsDb')) {
    .teacnv_annotation_require(c('ensembldb', 'AnnotationFilter'), auto_install)
    actual <- unique(na.omit(as.character(GenomeInfoDb::genome(GenomeInfoDb::seqinfo(annotations)))))
    if (!length(actual) || !all(.teacnv_annotation_assembly(actual) == key))
      stop('EnsDb assembly does not match ', genome, '.', call. = FALSE)
    gr <- ensembldb::transcripts(annotations,
      filter = AnnotationFilter::GeneNameFilter(gene_name),
      columns = c('tx_id', 'gene_id', 'gene_name', 'gene_biotype'), return.type = 'GRanges')
  } else if (inherits(annotations, 'GRanges')) {
    actual <- unique(na.omit(as.character(GenomeInfoDb::genome(annotations))))
    actual <- actual[nzchar(actual)]
    if (length(actual) && !all(.teacnv_annotation_assembly(actual) == key))
      stop('Supplied GRanges assembly does not match ', genome, '.', call. = FALSE)
    if (!length(actual)) warning('Supplied GRanges has no assembly metadata; verify coordinates match ', genome, '.', call. = FALSE)
    gr <- annotations
    if (!'gene_name' %in% names(S4Vectors::mcols(gr))) {
      if ('gene_symbol' %in% names(S4Vectors::mcols(gr))) gr$gene_name <- gr$gene_symbol
      else stop('annotations must contain gene_name or gene_symbol metadata.', call. = FALSE)
    }
    gr <- gr[as.character(gr$gene_name) %in% gene_name]
  } else stop('annotations must be a GRanges or EnsDb object.', call. = FALSE)
  if (!length(gr)) {
    warning('No requested genes found in ', genome, ': ', paste(gene_name, collapse = ', '), call. = FALSE)
    return(gr)
  }
  # Reassemble transcript spans from exon-level annotations when needed.
  # If explicit transcript records exist, use their coordinates directly.
  nm <- as.character(gr$gene_name)
  chrom <- as.character(GenomeInfoDb::seqnames(gr))
  st <- as.character(GenomicRanges::strand(gr))
  gid <- if ('gene_id' %in% names(S4Vectors::mcols(gr))) as.character(gr$gene_id) else nm
  gid[is.na(gid) | !nzchar(gid)] <- nm[is.na(gid) | !nzchar(gid)]
  locus <- paste(nm, gid, chrom, st, sep = '\r')
  tx_col <- intersect(c('tx_id', 'transcript_id', 'tx_name'), names(S4Vectors::mcols(gr)))
  selected <- lapply(split(seq_along(gr), locus), function(ii) {
    r <- gr[ii]
    if (length(tx_col)) {
      tx <- as.character(S4Vectors::mcols(r)[[tx_col[1L]]])
      valid <- !is.na(tx) & nzchar(tx)
      if (any(valid)) {
        groups <- split(which(valid), tx[valid])
        spans <- vapply(groups, function(i) max(BiocGenerics::end(r)[i]) - min(BiocGenerics::start(r)[i]) + 1, numeric(1))
        k <- groups[[which.max(spans)]]
        ans <- r[k[1L]]
        ans <- GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(ans),
          ranges = IRanges::IRanges(start = min(BiocGenerics::start(r)[k]), end = max(BiocGenerics::end(r)[k])),
          strand = GenomicRanges::strand(ans), seqinfo = GenomeInfoDb::seqinfo(ans))
        S4Vectors::mcols(ans) <- S4Vectors::mcols(r[k[1L]])
        return(ans)
      }
    }
    # Gene-level or transcript-level GRanges without transcript IDs.
    if ('type' %in% names(S4Vectors::mcols(r))) {
      type <- tolower(as.character(r$type))
      keep <- which(type %in% c('transcript', 'mrna'))
      if (!length(keep)) keep <- which(type == 'gene')
      if (!length(keep)) {
        warning('No transcript IDs or gene/transcript records for ', nm[ii[1L]], '; skipping.', call. = FALSE)
        return(r[0])
      }
      r <- r[keep]
    }
    r[which.max(BiocGenerics::width(r))]
  })
  ans <- do.call(c, unname(selected))
  # Prefer primary chromosome loci where available; never choose an alternate
  # contig over a known primary locus merely because its transcript is longer.
  if (key %in% c('hg19', 'hg38', 'mm10', 'mm39') && length(ans)) {
    ch <- sub('^chr', '', as.character(GenomeInfoDb::seqnames(ans)), ignore.case = TRUE)
    primary <- ch %in% c(as.character(seq_len(if (key %in% c('hg19','hg38')) 22L else 19L)), 'X','Y','M','MT')
    genes_primary <- as.character(ans$gene_name)[primary]
    ans <- ans[primary | !as.character(ans$gene_name) %in% genes_primary]
  }
  ans <- ans[order(match(as.character(ans$gene_name), gene_name), -BiocGenerics::width(ans))]
  missing <- setdiff(gene_name, as.character(ans$gene_name))
  if (length(missing)) warning('Genes not found: ', paste(missing, collapse = ', '), call. = FALSE)
  ans
}



#' @title CollapseToLongestTranscript()
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @import data.table
#' @export
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
#' @export

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


.teacnv_abs_layout <- function(data_list, genome, auto_install) {
  # One union axis for all data sets, using the plotting helpers from this file.
  seg_ids <- unique(unlist(lapply(data_list, function(d) d$.seg_key), use.names = FALSE))
  .teacnv_plot_coordinates(data.frame(binID = seg_ids), genome,
                           true_loc = TRUE, auto_install = auto_install)$chrLine
}

.teacnv_abs_apply <- function(d, chrLine = NULL, true_loc = TRUE) {
  if (true_loc) {
    k <- match(d$chrom, chrLine$chrom)
    if (anyNA(k)) stop("Shared chromosome layout is missing input chromosomes.")
    if (any(d$.seg_local_end > chrLine$length[k]))
      stop("Segment coordinates exceed reference lengths; check genome assembly.")
    d$length <- chrLine$length[k]
    d$chrPosStart_abs <- chrLine$chrPosStart_abs[k]
    d$chrPosEnd_abs <- chrLine$chrPosEnd_abs[k]
    d <- d[order(k, d$.bin_order, d$.seg_local_start, d$.seg_local_end), , drop = FALSE]
    d$segStart <- d$.seg_local_start + d$chrPosStart_abs
    d$segEnd <- d$.seg_local_end + d$chrPosStart_abs
    d$maploc <- d$.bin_order + d$chrPosStart_abs
  } else {
    # Equal spacing follows the same chromosome ordering as seg_plot without BSgenome.
    u <- unique(d$chrom)
    autos <- u[grepl("^[0-9]+$", u)]
    autos <- autos[order(as.numeric(autos))]
    special <- intersect(c("X", "Y", "M"), u)
    chr_order <- c(autos, special, sort(setdiff(u, c(autos, special))))
    d <- d[order(match(d$chrom, chr_order), d$.bin_order,
                 d$.seg_local_start, d$.seg_local_end), , drop = FALSE]
    d$maploc <- d$rank <- seq_len(nrow(d))
    d$segStart <- ave(d$rank, d$.seg_key, FUN = min)
    d$segEnd <- ave(d$rank, d$.seg_key, FUN = max)
  }
  d[c(".seg_local_start", ".seg_local_end", ".seg_key", ".bin_order")] <- NULL
  d
}


# #' @title abs_pos()
# #' @description add columns of chrosome absolute start and end for "segName" of a data.frame
# #' segName can be chr1-100037672-100038546 or chr1_100037672_100038546 format
# #' @export

# abs_pos <- function(data.bin,genome="hg38",true_loc=TRUE){
#     suppressMessages({
#     library(BSgenome)
#     library(dplyr)})
#   g <- getBSgenome(genome, masked=FALSE)
#   chrLine <- data.frame(chrom=1:24, length=seqlengths(g)[1:24])
#   chromnum=sapply(strsplit(data.bin$segName,'-|_|:'),'[',1)
#   chr_name=as.character(unique(chromnum))
#   if(grepl('chr',chromnum[1])){
#     chromnum=sapply(strsplit(chromnum,'hr'),'[',2) ## if chrommnum with chr
#   }else{
#     chromnum=chromnum
#     }
#   chromnum <- gsub("X","23",chromnum)
#   chromnum <- gsub("Y","24",chromnum)
#   data.bin$chrom <- chromnum
#   chr_name <- sapply(strsplit(chr_name,'hr'),'[',2)
#   chr_name_nm <- gsub("X","23",chr_name)
#   chr_name_nm <- gsub("Y","24",chr_name_nm)

#   if(true_loc){
#     chrLine <- chrLine[chrLine$chrom %in% as.character(unique(chromnum)),,drop=F] 
#     genomeLength <- sum(chrLine)
#     if(nrow(chrLine)>1){
#       chrPosStart = c(0,chrLine[1:(nrow(chrLine)-1),2])
#     }else(
#       chrPosStart = 0
#     )
#     chrLine$chrPosStart_abs <- cumsum(chrPosStart)
#     data.bin <- merge(data.bin,chrLine,by="chrom",all.x=T)
#     data.bin <- data.bin[order(as.integer(data.bin$chrom),data.bin$chrPosStart_abs),]
#     data.bin$segStart <- as.numeric(sapply(strsplit(data.bin$segName,'-|_|:'),'[',2))+as.numeric(data.bin$chrPosStart_abs)
#     data.bin$segEnd <- as.numeric(sapply(strsplit(data.bin$segName,'-|_|:'),'[',3))+as.numeric(data.bin$chrPosStart_abs)
    
#     #data.bin$Start <- sapply(strsplit(data.bin$segName,'-|_|:'),'[',2) #for segment start( not bin)
    
#     #maploc <- data.bin$chrPosStart_abs+data.bin$Start
#     intercept_x =chrLine$chrPosStart_abs
#     breakss<-cumsum(chrPosStart)+chrLine[1:(nrow(chrLine)-1),2]*0.5
#     which_chrom <- which(as.character(chrLine$chrom) %in% chr_name_nm)
#     breakss<-breakss[which_chrom]

#     #data.bin$maploc <- maploc

#   }else{
   
#     chromnum=chromnum[order(as.integer(gsub("[^0-9]", "", chromnum)))] #order by chr
#     chromnum_xaxis = table(chromnum)[order(as.integer(gsub("[^0-9]", "", names(table(chromnum)))))]#order by chr
#     chrline=cumsum(chromnum_xaxis)
#     maploc=1:nrow(data.bin)
#     breakss<-as.numeric(c(0,chrline[-length(chrline)])+chromnum_xaxis*0.5)
#     intercept_x =chrline
#     data.bin$maploc <- maploc
#     data.bin$rank <- maploc
#     data.bin <- data.bin %>%
#       dplyr::group_by(segName) %>%
#       dplyr::mutate(segStart=min(rank),segEnd=max(rank))%>%
#       as.data.frame()
#   }
# return(data.bin)
# }

#' Calculate absolute or rank-based genomic positions
#' @param data.bin Data frame with segName and optional bin coordinates.
#' @param genome Genome assembly or BSgenome object.
#' @param true_loc Use physical coordinates rather than equal bin ranks.
#' @param auto_install Install missing genome packages on demand.
#' @return Data frame with segment start/end and plotting coordinates.
#' @export
abs_pos <- function(data.bin, genome = "hg38", true_loc = TRUE,
                    auto_install = TRUE) {
  d <- .teacnv_abs_input(data.bin)
  layout <- if (true_loc) .teacnv_abs_layout(list(d), genome, auto_install) else NULL
  .teacnv_abs_apply(d, layout, true_loc)
}


# #' @title seg_compare()
# #' @param df1 data.frame input with c("segID","SegMean","segName","Start") on columns
# #                         segID  SegMean                segName      Start   
# #chr1-100037672-100038546     3 2.578573 chr1_9823588_121519947    100037672 
# #chr1-100132556-100133506     3 2.578573 chr1_9823588_121519947    100132556 
# #' @export

# seg_compare <- function(df1,df2,genome="hg38",ylim=NULL,ylab="Relative CN",ggtitle="",
#   hline=c(1,2),
#   seg_col=c("darkgrey","red"),
#   alpha=0.5){
#   suppressMessages({
#     library(BSgenome)
#     library(dplyr)
#     library(ggplot2)
#     })
#   g <- getBSgenome(genome, masked=FALSE)
#   chrLine <- data.frame(chrom=1:24, length=seqlengths(g)[1:24])
#   if(nrow(chrLine)>1){
#       chrPosStart = c(0,chrLine[1:(nrow(chrLine)-1),2])
#     }else{
#       chrPosStart = 0
#     }

#     chrLine$chrPosStart_abs <- cumsum(chrPosStart)

#   df1 <- abs_pos(df1)
#   df2 <- abs_pos(df2)

#   if(is.null(ylim)){
#     ylim=c(min(a$SegMean),max(a$SegMean)+1)
#   }
#   chromnum <- sapply(strsplit(df1$segName,'-|_|:'),'[',1)
#   chr_name <- as.character(unique(chromnum))
#   chr_name <- sapply(strsplit(chr_name,'hr'),'[',2)
#   chr_name_nm <- gsub("X","23",chr_name)
#   chr_name_nm <- gsub("Y","24",chr_name_nm)
#   chr_name_label <- chr_name
#   chr_name_label[as.character(chr_name_label)%in%c("19","21")] <- ""
#   breakss<-cumsum(chrPosStart)+chrLine[1:(nrow(chrLine)-1),2]*0.5
#   which_chrom <- which(as.character(chrLine$chrom) %in% chr_name_nm)
#   breakss<-breakss[which_chrom]

#   intercept_x =chrLine$chrPosStart_abs
#   p1 <- ggplot() +
#     #geom_point(aes(x = axis.loc, y = value.bin), alpha = 1, shape = 16,size=1) +
#     geom_hline(yintercept =hline,col="grey",linetype = 2)+
#     theme_bw() +
#     labs(x = "", y = ylab)+
#     ylim(ylim)+
#     ggtitle(ggtitle)+
#     geom_vline(xintercept =intercept_x,col="grey")+
#     theme(plot.margin = unit(c(0.5,0,0,0.5), 'lines'), #c(top, right, bottom, left)
#       panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
#           legend.position = "none",plot.title = element_text(size = 15),
#           axis.title=element_text(size=15),
#           axis.text = element_text(size = 15))+
#     scale_x_continuous(guide = guide_axis(check.overlap = TRUE),breaks=breakss,labels=chr_name_label)

#   p1 <- p1 +
#     #geom_line(aes(x = maploc,y=value.segment),col="red")+ 
#     geom_segment(
#       data = df1, 
#       mapping = aes(x=segStart, y=SegMean, xend=segEnd, yend=SegMean), 
#       col=seg_col[1],na.rm =T,size = 1, alpha = alpha,
#       inherit.aes = FALSE) +
#     geom_segment(
#       data = df2, 
#       mapping = aes(x=segStart, y=SegMean, xend=segEnd, yend=SegMean), 
#       col=seg_col[2],na.rm =T,size = 1,alpha = alpha,
#       inherit.aes = FALSE)
# return(p1)
# }
#' Compare two CNV segment profiles on a shared reference axis
#' @param df1,df2 Data frames containing segName and SegMean.
#' @param genome Genome assembly or BSgenome object.
#' @param ylim Optional y-axis limits.
#' @param ylab Y-axis label.
#' @param ggtitle Plot title.
#' @param hline Horizontal reference lines.
#' @param seg_col Two colours for the profiles.
#' @param alpha Segment opacity.
#' @param auto_install Install missing genome packages on demand.
#' @return ggplot object.
#' @export
seg_compare <- function(df1, df2, genome = "hg38", ylim = NULL,
                        ylab = "Relative CN", ggtitle = "", hline = c(1, 2),
                        seg_col = c("darkgrey", "red"), alpha = 0.5,
                        auto_install = TRUE) {
  d1 <- .teacnv_abs_input(df1)
  d2 <- .teacnv_abs_input(df2)
  if (!"SegMean" %in% names(d1) || !"SegMean" %in% names(d2))
    stop("Both input data frames must contain SegMean.", call. = FALSE)
  if (length(seg_col) != 2L) stop("seg_col must contain two colours.")
  # Critical: never calculate independent offsets for the two samples.
  chrLine <- .teacnv_abs_layout(list(d1, d2), genome, auto_install)
  d1 <- .teacnv_abs_apply(d1, chrLine)
  d2 <- .teacnv_abs_apply(d2, chrLine)
  d1$SegMean <- .teacnv_plot_num(d1$SegMean)
  d2$SegMean <- .teacnv_plot_num(d2$SegMean)
  # Remove duplicate bin-level repetitions of the same plotted segment.
  cols <- c("chrom", "segStart", "segEnd", "SegMean")
  d1 <- unique(d1[is.finite(d1$SegMean), cols, drop = FALSE])
  d2 <- unique(d2[is.finite(d2$SegMean), cols, drop = FALSE])
  if (!nrow(d1) || !nrow(d2)) stop("Each profile must have finite SegMean values.")
  vals <- c(d1$SegMean, d2$SegMean)
  if (is.null(ylim)) ylim <- c(min(vals), max(vals) + 1)
  ylim <- .teacnv_plot_limits(vals, ylim)
  ggplot2::ggplot() + ggplot2::theme_bw() +
    ggplot2::geom_hline(yintercept = hline, colour = "grey", linetype = 2) +
    ggplot2::geom_vline(xintercept = chrLine$chrPosStart_abs, colour = "grey") +
    ggplot2::geom_segment(data = d1, inherit.aes = FALSE,
      ggplot2::aes(x = segStart, xend = segEnd, y = SegMean, yend = SegMean),
      colour = seg_col[1L], linewidth = 1, alpha = alpha, na.rm = TRUE) +
    ggplot2::geom_segment(data = d2, inherit.aes = FALSE,
      ggplot2::aes(x = segStart, xend = segEnd, y = SegMean, yend = SegMean),
      colour = seg_col[2L], linewidth = 1, alpha = alpha, na.rm = TRUE) +
    ggplot2::labs(x = "", y = ylab, title = ggtitle) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::scale_x_continuous(
      limits = c(0, tail(chrLine$chrPosEnd_abs, 1L)),
      breaks = chrLine$chrPosStart_abs + chrLine$length/2,
      labels = chrLine$chrom, guide = ggplot2::guide_axis(check.overlap = TRUE)) +
    ggplot2::theme(plot.margin = grid::unit(c(.5, 0, 0, .5), "lines"),
      panel.grid = ggplot2::element_blank(), legend.position = "none",
      plot.title = ggplot2::element_text(size = 15),
      axis.title = ggplot2::element_text(size = 15), axis.text = ggplot2::element_text(size = 15))
}



#' @title seg_split()
#' @description split intersecting segments
#' @param binSet data.frame with bins in rowname
#'                                     binID Chromosome   Start     End
# chr1-826972-1013895   chr1-826972-1013895       chr1  826972 1013895
# chr1-1019104-1116833 chr1-1019104-1116833       chr1 1019104 1116833
# chr1-1121514-1350035 chr1-1121514-1350035       chr1 1121514 1350035
#' @export
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
#' @export
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
#' @export
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
#' @export
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




#' @param mt_atac scATAC count matrix (normalized), peaks on row and cells on column
#' @param mt_rna scRNA count matrix (normalized), genes on row and cells on column
#' @param cellMeta data.frame with cell annotation on the second column, cellID on the rowname and first column
#' @param segRatio_res segmentation result for the group-level which annotated by cellMeta
#' @param gene_ls vector of gene names for plot
#' @export
correlat_plot_v0 <- function(mt_atac,
                          mt_rna,
                          cnv_rna=NULL,
                          cellMeta,
                          segRatio_res,
                          annotations=NULL,
                          genome="hg38",
                          segRatio_Correct=TRUE,
                          df_seg4Correct=NULL,
                          expand_distance = 0,
                          FiltSegbyLen =TRUE,
                          segLen_lim = 1e6,
                          outdir = "./",
                          dotPlot4Seg=TRUE,
                          dotPlot4Gene = TRUE,
                          gene_ls=NULL,
                          label_term="Cluster"#"Chromosome"
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
  #covert to bulk
  mt_atac2_bulk <- pseudo_bulk_v2(mt_atac[,rownames(cellMeta)],group_ano = cellMeta,
                                  method ="mean",adjust_zero=TRUE)
  Ncol <- ncol(mt_atac2_bulk)
  mt_atac2_bulk$binID <- rownames(mt_atac2_bulk)
  mt_atac2_bulk_long <- tidyr::gather(mt_atac2_bulk,Cluster,Clmean_nCountATAC_peak,0:Ncol, factor_key=TRUE)
  mt_atac2_bulk_long$Chromosome <- sapply(strsplit(mt_atac2_bulk$binID,"-"),"[",1)
  mt_atac2_bulk_long$Start <- sapply(strsplit(mt_atac2_bulk$binID,"-"),"[",2)
  mt_atac2_bulk_long$End <- sapply(strsplit(mt_atac2_bulk$binID,"-"),"[",3)
  
  ##(2) process to segment matrix 
  segMT_long <- c()
  for(i in 1:length(segRatio_res)){
    Cname <- names(segRatio_res)[i]
    CNVdata <- segRatio_res[[i]]
    CNVdata <-CNVdata[CNVdata$binID %in% rownames(mt_atac),,drop=F]
    CNVdata$Cluster <- Cname
    CNVdata <- CNVdata[,c("binID","Chromosome","Start","End","binRatio","segName","SegMean","Cluster")]
    count_i <- mt_atac2_bulk[CNVdata$binID,Cname]
    CNVdata$Clmean_nCountATAC_peak <- count_i
    if(segRatio_Correct){
      CNVdata$SegMean <- as.numeric(CNVdata$Cluster)
      if(!is.null(df_seg4Correct)){
        df_seg4Correct <- df_seg4Correct[match(CNVdata$binID,df_seg4Correct$binID),,drop=F]
        CNVdata$segName <- df_seg4Correct$segName
      }
    }
    
    CNVdata <- CNVdata %>%
      dplyr::add_count(Cluster,segName) %>%
      as.data.frame()
    colnames(CNVdata)[colnames(CNVdata)=="n"] <- "nFeature_atac_seg"
    segMT_long <- rbind(segMT_long,CNVdata)
  }

  segMT_long <- segMT_long %>%
    dplyr::group_by(Cluster,segName) %>%
    dplyr::mutate(med_nCountATAC_seg = mean(Clmean_nCountATAC_peak,na.rm = T))%>%
    as.data.frame()
  by_cols <- intersect(colnames(segMT_long),colnames(mt_atac2_bulk_long))
  mt_atac2_bulk_long <- mt_atac2_bulk_long[mt_atac2_bulk_long$Chromosome%in%unique(segMT_long$Chromosome),,drop=F]
  
  ### to be used for plot
  mt_atac2_bulk_long <- merge(mt_atac2_bulk_long,segMT_long,by=by_cols,all=T) 
  mt_atac2_bulk_long$segStart <- as.numeric(sapply(strsplit(mt_atac2_bulk_long$segName,"_"),"[",2))
  mt_atac2_bulk_long$segEnd <- as.numeric(sapply(strsplit(mt_atac2_bulk_long$segName,"_"),"[",3))
  mt_atac2_bulk_long$segLen <- mt_atac2_bulk_long$segEnd-mt_atac2_bulk_long$segStart
  
  ##ATAC dataframe
  segMT_slim <- unique(mt_atac2_bulk_long[,c("Chromosome","segStart","segEnd","segName","SegMean","Cluster","segLen","med_nCountATAC_seg","nFeature_atac_seg")])
  segMT_slim <- segMT_slim[!is.na(segMT_slim$segName),,drop=F]
  
  #Find segment and gene intesect regions
  seg_gr <- GRanges(seqnames = segMT_slim$Chromosome,
                    ranges = IRanges(start = as.numeric(segMT_slim$segStart), end = as.numeric(segMT_slim$segEnd)),
                    strand = "*")
  seg_gr$segName <- segMT_slim$segName
  seg_gr$Cluster <- segMT_slim$Cluster
  
  peak_gr <- GRanges(seqnames = mt_atac2_bulk_long$Chromosome,
                     ranges = IRanges(start = as.numeric(mt_atac2_bulk_long$Start), end = as.numeric(mt_atac2_bulk_long$End)),
                     strand = "*")
  mcols(peak_gr) <- mt_atac2_bulk_long[,c("segName","SegMean","Cluster","binID","binRatio")]
  peak_gr$Clmean_nCountATAC_peak <- mt_atac2_bulk_long$Clmean_nCountATAC_peak
  
  if (is.null(x = annotations)) {
    #annotations <- Annotation(object = object[[peak.assay]])
    if(genome=="hg38"){
      require(EnsDb.Hsapiens.v86)
      annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    }else if (is.null(x = annotations)) {
      stop("Gene annotations not found")
    }
  }
  
  genes <- rownames(mt_rna)
  gene_loc <- annotate_gene(genes,genome="hg38",annotations=annotations)
  gene_loc_expan <- Extend(x = gene_loc,upstream = expand_distance,downstream = expand_distance)
  gr_gene_region <- paste(seqnames(gene_loc_expan),ranges(gene_loc_expan),sep="-")
  gene_loc_expan$gene_region <- gr_gene_region
  
  hits <- findOverlaps(gene_loc_expan, peak_gr)
  idx_seg <- subjectHits(hits)
  idx_gene <- queryHits(hits)
  ranges_ov <- gene_loc_expan[idx_gene,]
  mcols(ranges_ov) <- c(mcols(ranges_ov), mcols(peak_gr[idx_seg,]))
  
  #data.frame of gene and segment annotation
  gene_seg <- data.frame(segName=mcols(ranges_ov)$segName,gene_name=mcols(ranges_ov)$gene_name,Cluster=mcols(ranges_ov)$Cluster,
                         Clmean_nCountATAC_peak=mcols(ranges_ov)$Clmean_nCountATAC_peak,binID=mcols(ranges_ov)$binID, SegMean=mcols(ranges_ov)$SegMean)
  #gene expression: bulk 
  mt_rna_bulk <- pseudo_bulk_v2(mt_rna[,rownames(cellMeta)],group_ano = cellMeta,method ="mean",adjust_zero=TRUE)
  mt_rna_bulk$gene_name <- rownames(mt_rna_bulk)
  mt_rna_bulk_long <- tidyr::gather(mt_rna_bulk,Cluster,Clmean_nCountRNA_gene,0:Ncol, factor_key=TRUE)
  mt_rna_bulk_long <- merge(mt_rna_bulk_long,gene_seg,by=c("gene_name","Cluster"),all=T) #one gene may have multiple peaks in the same cluster
  
  gene_seg_dat <-  unique(mt_rna_bulk_long[,c("gene_name","Cluster","segName","Clmean_nCountRNA_gene")])
  gene_seg_dat <- gene_seg_dat[!is.na(gene_seg_dat$segName),]
  gene_seg_dat <- gene_seg_dat %>%
    dplyr::group_by(Cluster,segName) %>%
    dplyr::mutate(Cluster_nCountRNA_seg=mean(Clmean_nCountRNA_gene,na.rm = T))%>%
    ungroup()%>%
    add_count(Cluster,segName) %>%
    as.data.frame()
  colnames(gene_seg_dat)[colnames(gene_seg_dat)=="n"] <- "Cluster_nFeatureRNA_seg"
  
  mt_rna_bulk_long <- merge(mt_rna_bulk_long,gene_seg_dat,by=intersect(colnames(mt_rna_bulk_long),colnames(gene_seg_dat)),all=T) 
  
  ###RNA df
  mt_rna_bulk_long <- mt_rna_bulk_long %>%
    dplyr::group_by(Cluster,gene_name) %>%
    dplyr::mutate(Cluster_nCountATAC_gene=mean(Clmean_nCountATAC_peak,na.rm = T),  ####Mean
                  Cluster_segName_gene = median(SegMean,na.rm = T))%>%
    as.data.frame()
  
  mt_rna_bulk_l_slim <- mt_rna_bulk_long[!is.na(mt_rna_bulk_long$segName),,drop=F]
  mt_rna_bulk_l_slim <- unique(mt_rna_bulk_l_slim[,c("Cluster","segName","Cluster_nCountRNA_seg","Cluster_nFeatureRNA_seg")])
  #add RNA to ATAC dataframe
  segMT_slim <- merge(segMT_slim,mt_rna_bulk_l_slim,by=intersect(colnames(mt_rna_bulk_l_slim),colnames(segMT_slim)),all=T)
  
  if(!is.null(cnv_rna)){
    cnv_rna_bulk <- pseudo_bulk_v2(cnv_rna[,rownames(cellMeta)],group_ano = cellMeta,method ="median",adjust_zero=TRUE)
    cnv_rna_bulk$gene_name <- rownames(cnv_rna_bulk)
    cnv_rna_bulk_long <- tidyr::gather(cnv_rna_bulk,Cluster,med_inferCNV,0:Ncol, factor_key=TRUE)
    cnv_rna_bulk_long <- merge(cnv_rna_bulk_long,gene_seg,by=c("gene_name","Cluster"),all=T)
    
    ###CNV df
    cnv_seg <- unique(cnv_rna_bulk_long[,c("gene_name","Cluster","segName","med_inferCNV")])
    cnv_seg <- cnv_seg %>%
      dplyr::group_by(Cluster,segName) %>%
      dplyr::mutate(med_inferCNV_seg=median(med_inferCNV,na.rm = T))%>%
      as.data.frame()
    
    cnv_ge <- unique(cnv_seg[,c("gene_name","Cluster","med_inferCNV")])%>%
      dplyr::group_by(Cluster,gene_name) %>%
      dplyr::mutate(Cluster_inferCNV_gene=median(med_inferCNV,na.rm = T))%>%
      as.data.frame()
    cnv_seg <- merge(cnv_seg,cnv_ge,by=intersect(colnames(cnv_seg),colnames(cnv_ge)),all=T)
    cnv_seg <- cnv_seg[!is.na(cnv_seg$segName),]
    
    cnv_rna_bulk_long <- merge(cnv_rna_bulk_long,cnv_seg,by=intersect(colnames(cnv_rna_bulk_long),colnames(cnv_seg)),all=T)
    
    cnv_rna_bulk_l_slim <- cnv_rna_bulk_long[!is.na(cnv_rna_bulk_long$segName),,drop=F]
    cnv_rna_bulk_l_slim <- unique(cnv_rna_bulk_l_slim[,c("Cluster","segName","med_inferCNV_seg")])
    #add inferCNV to ATAC dataframe
    segMT_slim <- merge(segMT_slim,cnv_rna_bulk_l_slim,by=intersect(colnames(cnv_rna_bulk_l_slim),colnames(segMT_slim)),all=T)
    
  }
  
  color_c <- pal_ucscgb(alpha = 0.8)(length(sort(unique(segMT_slim$Cluster))))
  names(color_c) <- sort(unique(segMT_slim$Cluster))
  
  if(FiltSegbyLen){
    segMT_slim_filt <- segMT_slim[segMT_slim$segLen>1e6,,drop=F]
  }else{segMT_slim_filt <- segMT_slim}

  segMT_slim_filt$segLen_MB <- segMT_slim_filt$segLen/1e6
  segMT_slim_filt$Cluster_nFeatureRNA_segnom <- segMT_slim_filt$Cluster_nFeatureRNA_seg/segMT_slim_filt$segLen_MB
  segMT_slim_filt$nFeature_atac_segnom <- segMT_slim_filt$nFeature_atac_seg/segMT_slim_filt$segLen_MB
  
  #Plot for nCount_RNA vs nCount_ATAC
  if(dotPlot4Seg){
    outdir_seg <- paste0(outdir,"/seg-level");if(!file.exists(outdir_seg)){dir.create(outdir_seg,recursive=T)}
    
    plot_point(segMT_slim_filt,x_aes = "Cluster_nCountRNA_seg",y_aes="med_nCountATAC_seg",gtitle= paste0("Segment"),
               x_lab = "lognnCount_RNA",y_lab = "lognCount_ATAC",colour_item = "Cluster",size_item = "segLen_MB",size_lab="length(Mb)",
               color_value = color_c,outfile = paste0(outdir_seg,"/nCountATAC_nCountRNA_segment.pdf"),
               log_x=T,log_y=T,label_term=label_term)
    # plot_point(segMT_slim_filt,x_aes = "Cluster_nFeatureRNA_segnom",y_aes="nFeature_atac_segnom",gtitle= paste0("Segment"),
    #            x_lab = "lognFeature_RNA",y_lab = "lognFeature_ATAC",colour_item = "Cluster",size_item = "segLen_MB",size_lab="length(Mb)",
    #            color_value = color_c,outfile = paste0(outdir_seg,"/nFeatureATAC_nFeatureRNA_segment.pdf"),
    #            log_x=T,log_y=T,label_term="Chromosome")
    
    #Figure5: dot:segment (atacCNV vs nCountRNA/nFeatureRNA)
    plot_point(segMT_slim_filt,x_aes = "Cluster_nCountRNA_seg",y_aes="SegMean",gtitle= paste0("Segment"),
               x_lab = "lognCount_RNA",y_lab = "Relative CN (ATAC)",colour_item = "Cluster",size_item = "segLen_MB",size_lab="length(Mb)",
               color_value = color_c,outfile = paste0(outdir_seg,"/atacCNV_vs_nCountRNA_segment.pdf"),
               log_x=T,log_y=F,label_term=label_term)
    # plot_point(segMT_slim_filt,x_aes = "Cluster_nFeatureRNA_seg",y_aes="SegMean",gtitle= paste0("Segment"),
    #            x_lab = "lognFeature_RNA",y_lab = "Relative CN (ATAC)",colour_item = "Cluster",size_item = "segLen_MB",size_lab="length(Mb)",
    #            color_value = color_c,outfile = paste0(outdir_seg,"/atacCNV_vs_nFeatureRNA_segment.pdf"),
    #            log_x=T,log_y=F,label_term="Chromosome")
    
    #Figure4: dot:segment (CNV_atac vs CNV_rna)
    if(!is.null(cnv_rna)){
      plot_point(segMT_slim_filt,x_aes = "med_inferCNV_seg",y_aes="SegMean",gtitle= paste0("Segment"),
                 x_lab = "Relative CN (RNA)",y_lab = "Relative CN (ATAC)",colour_item = "Cluster",size_item = "segLen_MB",size_lab="length(Mb)",
                 color_value = color_c,outfile = paste0(outdir_seg,"/atacCNV_vs_rnaCNV_segment.pdf"),
                 log_x=F,log_y=F,label_term=label_term)
    }
      
  }

    #Figure2: dot:gene
    ###add inferCNV, ATAC and Ncells to RNA
    Ncell <- cellMeta%>%
      dplyr::count(Cluster)%>%
      as.data.frame()
    colnames(Ncell)[2]<- "Ncells"
    
    if(!is.null(cnv_rna)){
    mt_rna_bulk_long2 <- merge(mt_rna_bulk_long,cnv_rna_bulk_long,by=intersect(colnames(mt_rna_bulk_long),colnames(cnv_rna_bulk_long)),all.x=T)
    df_gene0<- unique(mt_rna_bulk_long2[,c("gene_name","Cluster","segName","binID","Cluster_segName_gene","Clmean_nCountRNA_gene","Cluster_nCountRNA_seg",
                                           "Cluster_nFeatureRNA_seg","Cluster_nCountATAC_gene","med_inferCNV")])
    df_gene0 <- merge(df_gene0,Ncell,by="Cluster",all=T)
    df_gene<- unique(df_gene0[,c("Cluster","gene_name","Clmean_nCountRNA_gene","Cluster_nCountATAC_gene",
                                 "med_inferCNV","Cluster_segName_gene","Ncells")]) 
    }else{
      mt_rna_bulk_long2 <- mt_rna_bulk_long
      df_gene0<- unique(mt_rna_bulk_long2[,c("gene_name","Cluster","segName","binID","Cluster_segName_gene","Clmean_nCountRNA_gene","Cluster_nCountRNA_seg",
                                             "Cluster_nFeatureRNA_seg","Cluster_nCountATAC_gene")])
      df_gene0 <- merge(df_gene0,Ncell,by="Cluster",all=T)
      df_gene<- unique(df_gene0[,c("Cluster","gene_name","Clmean_nCountRNA_gene","Cluster_nCountATAC_gene",
                                   "Cluster_segName_gene","Ncells")])  
    }
    if(dotPlot4Gene){
      outdir_g <- paste0(outdir,"/gene-level");if(!file.exists(outdir_g)){dir.create(outdir_g,recursive=T)}
      
    #dotplot Gene-level for Clusters
    if(is.null(gene_ls)){
      gene_ls <- unique(df_gene$gene_name)
    }
    for(ge in gene_ls){
      df_gene_gi <- unique(df_gene[df_gene$gene_name==ge,c("Cluster","gene_name","Cluster_segName_gene","Clmean_nCountRNA_gene","Cluster_nCountATAC_gene","Ncells"),drop=F])
      plot_point(df_gene_gi,x_aes = "Clmean_nCountRNA_gene",y_aes="Cluster_segName_gene",gtitle=ge,size_item="Ncells",size_lab="N cells",
                 x_lab = "lognCount_RNA",y_lab = "Relative CN (ATAC)",colour_item = "Cluster",color_gradient=F,
                 color_value = color_c,outfile = paste0(outdir_g,"/atacCNV_vs_nCountRNA_gene_",ge,".pdf"),
                 log_x=T,log_y=F,label_term="Cluster",legend_position="right",legend_color="none")

      
      plot_point(df_gene_gi,x_aes = "Clmean_nCountRNA_gene",y_aes="Cluster_nCountATAC_gene",gtitle= ge,
                 size_item="Ncells",size_lab="N cells",
                 x_lab = "lognCount_RNA",y_lab = "lognCount_ATAC",colour_item = "Cluster",color_gradient=F,
                 color_value = color_c,outfile = paste0(outdir_g,"/nCountATAC_nCountRNA_gene_",ge,".pdf"),
                 log_x=T,log_y=T,label_term="Cluster",legend_position="right",legend_color="none")
      # if(!is.null(cnv_rna)){
      #   df_gene_gi <- unique(df_gene[df_gene$gene_name==ge,c("Cluster","gene_name","Cluster_segName_gene","Clmean_nCountRNA_gene","Cluster_nCountATAC_gene","med_inferCNV","Ncells"),drop=F])
      #   plot_point(df_gene_gi,x_aes = "med_inferCNV",y_aes="Cluster_segName_gene",gtitle= ge,size_item="Ncells",size_lab="N cells",
      #              x_lab = "Relative CN (RNA)",y_lab = "Relative CN (ATAC)",colour_item = "Cluster",color_gradient=F,
      #              color_value = color_c,outfile = paste0(outdir_g,"/atacCNV_vs_rnaCNV_gene_",ge,".pdf"),
      #              log_x=F,log_y=F,label_term="Cluster",legend_position="right",legend_color="none")
      # }
      # 
    }
                      
  }
  
  return(list(segPlot.data=segMT_slim_filt,
              genePlot.data=mt_rna_bulk_long2
              ))
}


