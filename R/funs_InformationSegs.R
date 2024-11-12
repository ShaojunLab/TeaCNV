#Plotting density of CNVs across all clones for information segments
suppressMessages({
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(ggsci)
  library(LaplacesDemon)
  library(Seurat)
  library(Signac)
  library(futile.logger)
  library(Matrix) 
})
cols_Palette <- c("#B0D9A5","#A6DAEF","#D9BDD8","#E58579","#8AB1D2","#F9E9A4","#F1AEA7","#9D9ECD","#C9C780")
# source("fun_Grange.R")

custom_theme <-
    theme_classic()+ 
    theme(plot.background=element_blank(),
          legend.position='right',
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(color="black",size=15),
          axis.text=element_text(color="black",size=12),
          axis.line.y.right = element_blank(),
          axis.text.y.right=element_blank(),
          axis.ticks.y.right=element_blank(),
          axis.ticks=element_line(color="black",linewidth=0.5),
          plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
          axis.line.x.top = element_blank(),
          axis.ticks.x.top = element_blank()
        )#+guides(fill=guide_legend(title=legend.title))


#' @title densityPlot_segs()
#' @param clonalest_ls list of clonal-CNVs results
#' @return list of plots
#' @export
densityPlot_segs <- function(clonalest_ls,
  seg.bed,
  cols_Palette=NULL,
  xlim=NULL){
  res <- c()
  for(i in 1:length(clonalest_ls)){
    df_bin <- clonalest_ls[[i]]$input_BinSegRatio
    df_bin <- df_bin[,c("binID","Chromosome","Start","End","binRatio","length_bin","SegMean","segName")]
    df_seg <- clonalest_ls[[i]]$seg.dat
    df_seg <- df_seg[,c("segName","relativeCN","integerCN")]
    df_bin <- dplyr::left_join(df_bin,df_seg,by=c("segName"))
    df_bin$clone <- names(clonalest_ls)[i]
    df_bin$wei <-  df_bin$length_bin/sum(df_bin$length_bin)
    res <- rbind(res,df_bin)
  }

  if(is.null(cols_Palette)){
      cols_Palette <- c("#B0D9A5","#A6DAEF","#D9BDD8","#E58579","#8AB1D2","#F9E9A4","#F1AEA7","#9D9ECD","#C9C780")
    }
    if(length(cols_Palette)<length(clonalest_ls)){
      generate_random_colors <- function(N) {
        rgb(runif(N), runif(N), runif(N))
      }
      N <- length(clonalest_ls)-length(cols_Palette)
      colors_add <- generate_random_colors(N)
      cols_Palette <- c(cols_Palette,colors_add)
    }
    color_r <- cols_Palette[1:length(clonalest_ls)]
    names(color_r) <- sort(unique(res$clone))
    left_anno_cols <- list()
    left_anno_cols[["clone"]] <- color_r

  pl <- list()
  for(segi in 1:nrow(seg.bed)){
      chr_plt <- seg.bed[segi,1]
      segstart <- seg.bed[segi,2]
      segend <- seg.bed[segi,3]
      seg_plt <- paste(chr_plt,segstart,segend,sep="-")

      res_sub <- res[res$Chromosome %in% chr_plt & res$Start >= segstart & res$End<=segend,]
      res_sub <- res_sub %>%
        dplyr::group_by(integerCN) %>%
        dplyr::mutate(relativeCN_mean =mean(relativeCN))%>%
        as.data.frame()
    integerCNVcolumn <- grep("integerCN",colnames(res_sub))
    labs <- unique(res_sub[integerCNVcolumn])
    labs <- labs[!is.na(labs[,1]),]
    if(is.null(xlim)){
      x_min <- floor(quantile(as.numeric(res_sub$binRatio),0.01,na.rm=TRUE))
      x_max <- ceiling(quantile(as.numeric(res_sub$binRatio),0.99,na.rm=TRUE))
      xlim <- c(x_min,x_max)
    }
 

  custom_theme <-
    theme_classic()+ 
    theme(plot.background=element_blank(),
          legend.position='right',
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(color="black",size=15),
          axis.text=element_text(color="black",size=12),
          axis.line.y.right = element_blank(),
          axis.text.y.right=element_blank(),
          axis.ticks.y.right=element_blank(),
          axis.ticks=element_line(color="black",linewidth=0.5),
          plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
          axis.line.x.top = element_blank(),
          axis.ticks.x.top = element_blank()
        )#+guides(fill=guide_legend(title=legend.title))

    ht2= ggplot() + 
          geom_density(data = subset(res_sub), aes(x = binRatio, color = clone), 
                       lwd = 0.5, adjust =1.5) +
          geom_vline(xintercept = unique(res_sub$relativeCN_mean),col="grey",linetype = "dotted")+
          scale_x_continuous(
            limits = xlim,
              sec.axis = sec_axis( trans = ~.,breaks=unique(res_sub$relativeCN_mean)[!is.na(unique(res_sub$relativeCN_mean))], labels =labs, name="")
            )+
          scale_color_manual(values = left_anno_cols[["clone"]])+
          labs(title=seg_plt,x="Ratio(Bin)")+custom_theme
        
        pl[[seg_plt]] <- ht2
  }
  return(pl)
}



#' @title find_infoSeg()
#' @param clonalest_ls list of clonal-CNVs results
#' @export
find_infoSeg <- function(clonalest_ls,
  seg_len_min=2e6,
  doPlot=TRUE,
  outdir="./",
  cols_Palette=NULL,
  p.adj_cutoff = 0.05,
  ratio_diff_cutoff=NULL){
  res <- c()
  for(i in 1:length(clonalest_ls)){
    df_bin <- clonalest_ls[[i]]$input_BinSegRatio
    df_bin <- df_bin[,c("binID","Chromosome","Start","End","binRatio","length_bin","SegMean","segName")]
    df_seg <- clonalest_ls[[i]]$seg.dat
    df_seg <- df_seg[,c("segName","relativeCN","integerCN")]
    df_bin <- left_join(df_bin,df_seg,by=c("segName"))
    df_bin$clone <- names(clonalest_ls)[i]
    df_bin$wei <-  df_bin$length_bin/sum(df_bin$length_bin)
    res <- rbind(res,df_bin)
  }
  #select differential CNV region accross clones
  res_cnv <- unique(res[,c("Chromosome","segName","integerCN","clone")])
  res_cnv <- na.omit(res_cnv[res_cnv$integerCN != "2",])
  res_cnv$segstart <- as.numeric(sapply(strsplit(res_cnv$segName,"_|-|:"),"[",2))
  res_cnv$segend <- as.numeric(sapply(strsplit(res_cnv$segName,"_|-|:"),"[",3))
  filtered_res_cnv <- res_cnv %>%
    dplyr::group_by(Chromosome, integerCN) %>%
    dplyr::filter(n_distinct(clone) != length(unique(res_cnv$clone))) %>%
    ungroup()%>%
    as.data.frame()

  # table(filtered_res_cnv$segName,filtered_res_cnv$clone)
  merged_df <- filtered_res_cnv %>%
    # Create a grouping variable to indicate whether the integerCN of adjacent rows is the same
    dplyr::group_by(Chromosome,clone)%>%
    dplyr::mutate(newSeg = cumsum(integerCN != lag(integerCN, default = first(integerCN)) | segstart >= lag(segend + seg_len_min, default = first(segstart)))) %>%
    ungroup()%>%
    dplyr::group_by(Chromosome,newSeg, integerCN) %>%
    summarise(
      start = min(segstart),    # 
      end = max(segend)          # 
    ) %>%
    ungroup()%>%
    dplyr::distinct(Chromosome, start, end, .keep_all = TRUE)   %>%
    #filter(Chromosome=="chr22")%>%
    as.data.frame()


  chr_mum <- gsub("chr","",merged_df$Chromosome)
  chr_mum <- gsub("X|x","23",chr_mum)
  chr_mum <- gsub("Y|y","24",chr_mum)
  merged_df$chr_mum  <- as.numeric(chr_mum)
  merged_df <- merged_df[order(merged_df$chr_mum,merged_df$start),]

  if(doPlot){
    if(is.null(cols_Palette)){
      cols_Palette <- c("#B0D9A5","#A6DAEF","#D9BDD8","#E58579","#8AB1D2","#F9E9A4","#F1AEA7","#9D9ECD","#C9C780")
    }
    if(length(cols_Palette)<length(clonalest_ls)){
      generate_random_colors <- function(N) {
        rgb(runif(N), runif(N), runif(N))
      }
      N <- length(clonalest_ls)-length(cols_Palette)
      colors_add <- generate_random_colors(N)
      cols_Palette <- c(cols_Palette,colors_add)
    }
    color_r <- cols_Palette[1:length(clonalest_ls)]
    names(color_r) <- sort(unique(res$clone))
    left_anno_cols <- list()
    left_anno_cols[["clone"]] <- color_r
  }
  pl <- list()
  merged_df$InfoSeg <- 0
  merged_df$diffClones <- NA
  for(segi in 1:nrow(merged_df)){
    chr_plt <- merged_df$Chromosome[segi]
    segstart <- merged_df$start[segi]
    segend <- merged_df$end[segi]
    seg_plt <- paste(chr_plt,segstart,segend,sep="_")

    res_sub <- res[res$Chromosome %in% chr_plt & res$Start >= segstart & res$End<=segend,]
    # table(res_sub$segName,res_sub$clone)
    # table(res_sub$segName,res_sub$integerCN)

    #test differential ratio of Infomation segmens
    clones <- unique(res_sub$clone)
    if(length(clones)>1){
      diffClones <- c()

      for(cli_indx in 1:(length(clones)-1)){
        cli <- clones[cli_indx]
        for(clj_index in (cli_indx+1):(length(clones))){
          clj <- clones[clj_index]
          df.long <- res_sub[res_sub$clone %in% c(cli, clj),c("binID","segName","clone","binRatio"),drop=F]
          df.long <- na.omit(df.long)
          colnames(df.long) <- c("binID","segID","group","ratio")
          groupMean <- df.long %>%
              dplyr::group_by(group)%>%
              dplyr::summarise(groupMean = median(ratio,na.rm=TRUE))%>% as.data.frame()
          diffMean <-  groupMean%>%
              summarise(diffMean =abs(groupMean[group%in%cli]- groupMean[group%in%clj]))%>% as.data.frame()

          ###differential test
          stat.test <- df.long %>%
              dplyr::group_by(group) %>%
              dplyr::filter(dplyr::n() >= 2) %>%
              dplyr::ungroup()%>%
              t_test(ratio ~ group) %>%
              adjust_pvalue(method = "BH") %>%
              add_significance()
          stat.test <- as.data.frame(stat.test)
          stat.test$diffMean <- diffMean    
          stat.test$mode <- 1
          if(!is.null(ratio_diff_cutoff)){
            stat.test$mode[stat.test$p.adj < p.adj_cutoff & stat.test$diffMean>=ratio_diff_cutoff] <- 2 #differential
         
          }else{
              stat.test$mode[stat.test$p.adj < p.adj_cutoff ] <- 2 #differential
          }
          if(stat.test$mode==2){
            merged_df$InfoSeg[segi] <- 1
            diffClones <- c(diffClones,paste(cli,clj,sep="_"))
          }
        }
      }
      merged_df$diffClones[segi] <- paste(diffClones,collapse = ",")
      rm(diffClones)
    }

    res_sub <- res_sub %>%
    dplyr::group_by(integerCN) %>%
    dplyr::mutate(relativeCN_mean =mean(relativeCN))%>%
    as.data.frame()

    if(stat.test$mode==2 & doPlot ){
      pltdata <- res_sub
      integerCNVcolumn <- grep("integerCN",colnames(pltdata))
      labs <- unique(pltdata[integerCNVcolumn])
      labs <- labs[!is.na(labs[,1]),]

      ##Custermrized

      custom_theme <-
        theme_classic()+ 
        theme(plot.background=element_blank(),
              legend.position='right',
              plot.title = element_text(hjust = 0.5),
              axis.title = element_text(color="black",size=15),
              axis.text=element_text(color="black",size=12),
              axis.line.y.right = element_blank(),
              axis.text.y.right=element_blank(),
              axis.ticks.y.right=element_blank(),
              axis.ticks=element_line(color="black",linewidth=0.5),
              plot.margin = unit(c(0.5,0.5,0,0), 'lines'),
              axis.line.x.top = element_blank(),
              axis.ticks.x.top = element_blank()
            )#+guides(fill=guide_legend(title=legend.title))
      ht2= ggplot() + 
        geom_density(data = subset(pltdata), aes(x = binRatio, color = clone), 
                     lwd = 0.5, adjust =1.5) +
        # geom_density(data = subset(pltdata, clone %in%c("1","2") ), aes(x = binRatio, color = clone), 
        #              lwd = 0.5, adjust =1) +
        geom_vline(xintercept = unique(pltdata$relativeCN_mean),col="grey",linetype = "dotted")+
        scale_x_continuous(
          limits = c(0,3),
            sec.axis = sec_axis( trans = ~.,breaks=unique(pltdata$relativeCN_mean)[!is.na(unique(pltdata$relativeCN_mean))], labels =labs, name="")
          )+
        scale_color_manual(values = left_anno_cols[["clone"]])+
        labs(title=seg_plt,x="Ratio(Bin)")+custom_theme
      
      pl[[seg_plt]] <- ht2
    }
  }
  merged_df <- merged_df[merged_df$InfoSeg>0,,drop=F]
  merged_df <- merged_df[,c("Chromosome","start","end","diffClones","integerCN")]

  if(doPlot){
    ncol <- ifelse(length(pl)>4,4,length(pl))
    width <- ifelse(length(pl)>4,14,3.5*length(pl))
    pcom <- cowplot::plot_grid(plotlist=pl,ncol = ncol,align="v")
    ggsave(paste0(outdir, "/DensityPlot_ColnalBinRatio_InfoSegs.pdf"),pcom, width=width, height=2.5*ceiling(length(pl)/4),device = pdf,bg="white")  

  }

  #generate clonal CNV matrix for InfoSegs
  res2 <- res[,c("Chromosome","Start","End","integerCN","clone")]
  res_wide <- res2 %>%
   tidyr::pivot_wider(names_from = clone, values_from = integerCN)%>%as.data.frame()
  res2 <- align_Grange2bin(res_wide,merged_df[,1:4])
  res2 <- res2[!is.na(res2$diffClones),,drop=FALSE]
  rownames(res2) <- res2$binID
  InfoSegs <- paste(res2$chr.subject,res2$start.subject,res2$end.subject,sep="_")
  CN_InfoSegs <- res2[,4:(ncol(res2) - 5),drop=F]
  CN_InfoSegs$InfoSegs <- InfoSegs
  CN_InfoSegs <- na.omit(CN_InfoSegs)

  if(doPlot){
  clone_info <- data.frame(row.names=colnames(CN_InfoSegs)[1:(ncol(CN_InfoSegs)-1)],clone=colnames(CN_InfoSegs)[1:(ncol(CN_InfoSegs)-1)])
    cols_Palette <- c("#B0D9A5","#A6DAEF","#D9BDD8","#E58579","#8AB1D2","#F9E9A4","#F1AEA7","#9D9ECD","#C9C780")
    clone_info <- data.frame(row.names=colnames(CN_InfoSegs)[1:(ncol(CN_InfoSegs)-1)],clone=colnames(CN_InfoSegs)[1:(ncol(CN_InfoSegs)-1)])
    color_r <- cols_Palette[1:length(unique(clone_info$clone))]
    names(color_r) <- sort(unique(clone_info$clone))
    
    left_anno_cols <- list()
    left_anno_cols[["clone"]] <- color_r
    height <- ifelse(ncol(CN_InfoSegs)>2,0.35*(ncol(CN_InfoSegs))+0.5,ifelse(ncol(CN_InfoSegs)==1,1.2,1.5))
    p_cloneCN <- heatmap4peakMt(mat=CN_InfoSegs[,1:(ncol(CN_InfoSegs)-1)],
                          meta_info=clone_info,
                          sep_by="_",
                          outdir= outdir,
                          value.type="CNV",
                          clust_rows=F,
                          show_legend_row = T,
                          legend_titles="integer CN",
                          fileout_name=paste0("heatmap_cloneCNA_InfoSegs"),
                          col_list=left_anno_cols,
                          column.title = NULL,
                          width=10,height=height,device="pdf") 
  }
  return(list(merged_df=merged_df,CN_InfoSegs_mt=CN_InfoSegs))
}





