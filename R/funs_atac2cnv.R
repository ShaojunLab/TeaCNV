# This is a function to inferring clonal CNVs from scATAC-seq.
suppressMessages({
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(ggsci)
  library(Signac)
  library(futile.logger)
})


#' @title CNratioInfer()
#' @description infer CN ratio for group and segmentation
#' @export

CNratioInfer <- function(input_matrix,inputcellano,
                         cytoBand,
                         outdir=NULL,
                         group2bulkby="subCluster",
                         min_cells_in_group = 1,
                         normallabel="reference",
                         cells.ref = NULL,
                         seg_method="PELT", # "robseg",#"gaussian", #
                         seg.count.lim = 80,
                         penalty = 1,
                         segValue_method="density",
                         FiltSeg= TRUE,
                         SegLen_min=1e6,
                         SegSize_min = 5,
                         plot_ylim= NULL,
                         outFigure = TRUE,
                         outFigureName="CNRatio_subCluster.pdf",
                         color_dot=FALSE,
                         output_ref.value = TRUE,
                         outRefValueName="meanCount_Ref.rds",
                         outRDS=TRUE,
                         RDSname = "segScore_median.rds") {
  if (is.null(outdir)) {
    flog.error("Error, outdir is NULL, please provide a path.")
    stop("outdir is NULL, please provide a path.")
  }
  if(outdir != "." && !file.exists(outdir)){dir.create(outdir,recursive=T)}
  flog.info(paste0("Starting inferring clonal CNV ratio for count matrix."))
  
  
    ###-----------------###
    #4. generate bulk-level mat group by cluster
    ###-----------------###
    cell_ano2_all <- data.frame(cellid=rownames(inputcellano),group=as.character(inputcellano[,group2bulkby]))
    if(!is.null(cells.ref)){
      cell_ano_b <- data.frame(cellid=rownames(inputcellano),group="obs")
      cell_ano_b$group[cell_ano_b$cellid %in%cells.ref ] <- "ref"
      mtx_ref_bulk <- pseudo_bulk_v2(input_matrix[,cell_ano_b$cellid],group_ano = cell_ano_b,method ="mean",adjust_zero=TRUE)
      ref.value <- mtx_ref_bulk[,"ref",drop=F]
    }else if((!is.na(normallabel)) & any(cell_ano2_all$group %in% normallabel)){
      cell_ano2_all$group[cell_ano2_all$group %in% normallabel] <- normallabel
      col_ref <- normallabel
    }else{
      col_ref <- NA
    }
    element_counts <- table(cell_ano2_all$group)
    cell_ano2_all_filt <- subset(cell_ano2_all, cell_ano2_all$group %in% names(element_counts[element_counts >= min_cells_in_group]))
    
    mtx_qcs_bulk <- pseudo_bulk_v2(input_matrix[,cell_ano2_all_filt$cellid],group_ano = cell_ano2_all_filt,method ="mean",adjust_zero=TRUE)
    ###-----------------###
    #5. cluster-level CN Ratio
    ###-----------------### 
    if(!is.null(cells.ref)){
      col_obs <- colnames(mtx_qcs_bulk)
    }else if(!is.na(col_ref)){
      ref.value <- mtx_qcs_bulk[,col_ref,drop=F]
      col_obs <- colnames(mtx_qcs_bulk)[!grepl(col_ref,colnames(mtx_qcs_bulk))]
    }else{
      ref.value <-data.frame(rowMeans(mtx_qcs_bulk,na.rm = TRUE))
      col_obs <- colnames(mtx_qcs_bulk)
    }
    if(output_ref.value){
      saveRDS(ref.value,paste0(outdir,"/",outRefValueName))
    }
    
    mtx_qcR_bulk0 <- as.matrix(mtx_qcs_bulk[,col_obs]/ref.value[,1]) 
    
    colnames(mtx_qcR_bulk0) <- col_obs
    rownames(mtx_qcR_bulk0)<- rownames(mtx_qcs_bulk)
    is.na(mtx_qcR_bulk0) <- sapply(mtx_qcR_bulk0, is.infinite)
    is.na(mtx_qcR_bulk0) <- sapply(mtx_qcR_bulk0, is.nan)
    mtx_qcR_bulk0 <- na.omit(mtx_qcR_bulk0)


  # mtx_qcR_bulk_clean <- apply(mtx_qcR_bulk0,2,function(x,max_threshold,min_threshold){
  #   bins <- x
  #   is.na(bins) <- sapply(bins, is.infinite)
  #   is.na(bins) <- sapply(bins, is.nan)
  #   bin_value <- na.omit(bins)
  #   #max_threshold <- round(quantile(bins,0.9))
  #   bins[(!is.na(bins)) &(bins>=max_threshold)] <- NA
  #   bins[(!is.na(bins)) &(bins<=min_threshold)] <- NA
  #   return(bins)
  # },max_threshold,min_threshold)
  # 
  #outplot_path <- paste0(outdir,"/Figures");if(!file.exists(outplot_path)){dir.create(outplot_path,recursive=T)}
  # if(is.null(plot_ylim)){
  #   plot_ylim <- round(quantile(mtx_qcR_bulk0,c(0.01,0.99)),2)
  # }
  SegScore_total <- suppressMessages(segment4df(mtx_qcR_bulk0,  
                                 cytoBand = cytoBand,
                                 plot_ylim = plot_ylim,
                                 doFilt= FiltSeg,
                                 segLen_min= SegLen_min,
                                 segSize_min = SegSize_min,
                                 color_dot = color_dot,
                                 outpath = outdir,
                                 outPlot = FALSE,
                                seg.method=seg_method,
                                segValue.method=segValue_method,
                                 penalty = penalty,
                                seg.count.lim =seg.count.lim,
                                 filename = "subClust"))
  if(outFigure){
    p_ls <- lapply(SegScore_total$segPlots,function(x){y=x$ggarranged_p;return(y)})
    pcom <- cowplot::plot_grid(plotlist=p_ls,ncol = 1,align="v")
    ggsave(paste0(outdir, "/",outFigureName),pcom, width=14, height=3*ncol(mtx_qcR_bulk0),device = pdf,bg="white")
    
  }


  SegScore_total$ref.value <- ref.value

 
  if(outRDS){
    saveRDS(SegScore_total,paste0(outdir,"/",RDSname))
  }
  return(SegScore_total)
  flog.info(print("Ratio matrix done!"))
}

#' @title CNratioInfer.sc()
#' @description infer CN ratio for individual single cell
#' @param prop.zero numeric vector with the same length of rows of input_matrix.
#' @export
CNratioInfer.sc <- function(input_matrix,
                            cell.name,
                            inputcellano,
                            cytoBandFile,
                            group2bulkby="subCluster",
                            normallabel="reference",
                            seg.ref=NULL,
                            seg_method="robseg",
                            segValue.method="density",
                            outdir=".",
                            scale_value = TRUE,
                            prop_zero = NULL,
                            outPlot = TRUE,
                            outFigureName.prefix="segRatio_sc",
                            add_yline = NULL){
  if(outdir != "." && !file.exists(outdir)){dir.create(outdir,recursive=T)}
  cell_ano2_all <- data.frame(row.names = rownames(inputcellano),cellid=rownames(inputcellano),group=as.character(inputcellano[,group2bulkby]))
  if((!is.na(normallabel)) & any(cell_ano2_all$group %in% normallabel)){
    cell_ano2_all$group[cell_ano2_all$group %in% normallabel] <-"reference"
    col_ref <- "reference"
  }else{
    stop("No normal reference provided!")
  }

  cells_ref <-  rownames(cell_ano2_all)[cell_ano2_all$group%in%"reference"]
  cells_obs <- rownames(cell_ano2_all)[!cell_ano2_all$group %in% "reference"]
  col_indx <- match(cell.name,colnames(input_matrix))
  mtx_i <- input_matrix[,col_indx,drop=F]
  if(scale_value){
    #scale peak count by non-zero percentage
    if(is.null(prop_zero)){
      Nzero <- rowSums(input_matrix[,cells_obs]==0)
      prop_zero <- Nzero/ncol(input_matrix[,cells_obs])
    }
    mtx_ic <- mtx_i*(1-prop_zero)
  }else{
    mtx_ic <- mtx_i
  }
  mtx_ic <- mtx_ic[mtx_ic!=0,,drop=F]
  peak_used <- rownames(mtx_ic)
  mtx_ref <- input_matrix[peak_used,c(cells_ref),drop=F]
  mtx_sc <- cbind(mtx_ic,mtx_ref)
  
  cellMeta <- inputcellano[c(cell.name,cells_ref),,drop=F]
  cellMeta_new <- data.frame(row.names = rownames(cellMeta),group=rownames(cellMeta))
  cellMeta_new[cells_ref,1] <- "reference"
  cell_info2 <- data.frame(rownames(cellMeta_new),cellMeta_new)
  #flog.info("Calculate single-cell level ratio relative to reference.") 
  ratio_res <- CNratio_peak(mtx_sc,cell_ano=cell_info2,colname_ref="reference",Norm_by_ref=F) 
  is.na(ratio_res) <- sapply(ratio_res, is.infinite)
  is.na(ratio_res) <- sapply(ratio_res, is.nan)
  ratio_res <- na.omit(ratio_res)
  
  
  #remove outliers
  # ratio_res_clean <- filt_peak_perChr(ratio_res,zscore.lim = 2,split_by = "-")
  #seg_plot(data.bin=ratio_res_clean,value.bin=colnames(ratio_res),plot_seg=F,outPlot=F)
  if(is.null(seg.ref)){
    SegScore_sc <- suppressMessages(segment4df(ratio_res,  
                                               cytoBand = cytoBandFile,
                                               doFilt= F,
                                               segLen_min= 1,
                                               color_dot = F,
                                               outpath = outdir,
                                               outPlot = outPlot,
                                               seg.method=seg_method,
                                               add_yline = add_yline,
                                               filename = outFigureName.prefix))
  }else{
   
    seg.df <- seg.ref[rownames(ratio_res),"segName",drop=F]
    Nseg <- seg.df %>%
      dplyr::group_by(segName) %>%
      dplyr::add_count() %>%
      dplyr::mutate(LastRow = row_number() == n)
    changepoint_matrix <- data.frame(row.names =rownames(ratio_res),breakpoint= rep(0,nrow(ratio_res)) )
    changepoint_matrix$breakpoint[Nseg$LastRow] <- 1
    
    seg_bulk <- list(data_matrix=ratio_res,changepoint_matrix=changepoint_matrix)
    SegScore_sc = seg_score(seg_bulk,method=segValue.method)
    #seg_plot(data.bin=segScore$seg_score_binLevel[[1]],value.bin="binRatio",plot_seg=T,outPlot=F,plotDir=outdir,name.data=names(ratio_res))
    
    #replace the segName by seg.ref
    bin_df <- SegScore_sc$seg_score_binLevel[[1]]
    bin_df_ref <- seg.ref[rownames(ratio_res),c("binID","segName"),drop=F]
    colnames(bin_df_ref)[2] <- "segName.ref"
    bin_df <- dplyr::left_join(bin_df,bin_df_ref,by="binID")
    segNames <- unique(bin_df[,c("segName","segName.ref")])
    
    bin_df <-  bin_df[,!grepl("^segName$",colnames(bin_df)),drop=F]
    colnames(bin_df)[grepl("segName.ref",colnames(bin_df))] <- "segName"
    rownames(bin_df) <- bin_df$binID
      
    seg_df <- SegScore_sc$seg_score_segLevel[[1]]
    seg_df <- dplyr::left_join(seg_df,segNames,by="segName")
    seg_df <- seg_df[,!grepl("^segName$",colnames(seg_df)),drop=F]
    colnames(seg_df)[grepl("segName.ref",colnames(seg_df))] <- "segName"
    
    SegScore_sc$seg_score_binLevel[[1]] <- bin_df
    SegScore_sc$seg_score_segLevel[[1]] <- seg_df
  }
  
  return(SegScore_sc)
}



#' @title pcaClustCN()
#' @description Perform (1) convert peak to bin matrix if customized required;
#' (2) Calculate single-cell level ratio relative to bulk normal reference; 
#' (3) PCA clustering based on whole genome bin-cell ratio matrix;
#' (4) Calculate sub-cluster level CN ratio
#' @return list of sub-cluster level segmentation, ratio, and CNVs
#' @param cnv.ref.df data.frame of reference CNV and ratio.map if 'optimalClust' is TRUE
#' @export
pcaClustCN <- function(mtx_qcn_lim,
                       cell_anno,
                       ref_group_names,
                       chrs=NULL,
                       ref_cells=NULL,
                       delimit = "-",
                       cytoBandFile=NULL,
                       binSize = NULL, 
                       ClustOnRatio=TRUE,
                       seu_resolution=1,
                       min_cells_in_group=10,
                       OutBinMat =TRUE,
                       outputFigures=TRUE,
                       heatmap_name = NULL,
                       count_lim = NULL,
                       Figure.width=10,
                       Figure.height=6,
                       FigureName_prefix=NULL,
                       FiltSeg= TRUE,
                       segRatio_on_bin = TRUE,
                       SegLen_min=1e6,
                       SegSize_min = 5,
                       seg_method = "PELT",# "robseg",#"gaussian",#
                       segValue_method="density",
                       seg.count.lim=80,
                       penalty = 1,
                       CorlorDot4Ratio = TRUE,
                       outdir=".",
                       optimalClust=FALSE,
                       mse.min = 0.01,
                       cnv.ref.df =NULL,
                       global.norm=FALSE
                       ){
  message("Start bin-ratio-PCA-clustering...")
  if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
  outdir_sub <- outdir
  if(!is.null(cytoBandFile)){
    if (Reduce("|", is(cytoBandFile) == "character")) {
      if(substr(cytoBandFile, nchar(cytoBandFile)-2, nchar(cytoBandFile)) %in% c("txt","tsv","csv")){
        cytoBandFile <- read.table(cytoBandFile, sep="\t", header=F, check.names=FALSE)
        # cytoBand <- read.table(connection <- gzfile(cytoBand, 'rt'), sep="\t", header=F, check.names=FALSE)
        # close(connection)
        cytoBandFile <- as.data.frame(cytoBandFile)
      }
    }else if(Reduce("|", is(cytoBandFile) %in% c("data.frame"))){
    }else{
      message("cytoBandFile isn't recognized as data.frame, or filename")
    }
    #check the columns of cytoBand
    if(ncol(cytoBandFile)==5){
      colnames(cytoBandFile)=c("chrom","chromStart","chromEnd","name","gieStain")
    }else{
      stop("Error, the columns of cytoBandFile do NOT match the columns of the cytoBandFile table in the database underlying the UCSC Genome Browser")
    }
  } 
  #outdir_sub <- paste0(outdir,"/subCluster");if(!file.exists(outdir_sub)){dir.create(outdir_sub,recursive=T)}
  if(is.null(FigureName_prefix)){
   
    CNVresFile_subC <- paste0(outdir_sub,"/CNV_res_subClust.rds")
    outFigureName <- paste0("subClust_CNRatio.pdf")
    SegScore_total_subCName <- paste0("segScore_median_subClust.rds")
  }else{
    
    CNVresFile_subC <- paste0(outdir_sub,"/",FigureName_prefix,".CNV_res_subClust.rds")
    outFigureName <- paste0(FigureName_prefix,".subClust_CNRatio.pdf")
    SegScore_total_subCName <- paste0(FigureName_prefix,".segScore_median_subClust.rds")
   # outCNVnm<- paste0(FigureName_prefix,".CNV_subClust")
  }

  message("(1-3) Process matrix, normalization, ratio, PCA-clustering...")
 if(is.null(binSize)){
   binSize <- ceiling(dim(mtx_qcn_lim)[1]/2000)
 }
  
  clust_res <- binClust(inputMat=mtx_qcn_lim,cellMeta=cell_anno,
                            ref_group_names=ref_group_names,
                            chrs=chrs,
                            delimit=delimit,
                            onRatio = ClustOnRatio,
                            outdir=outdir_sub,
                            window = binSize,
                            seu_resol = seu_resolution,
                            OutMat =OutBinMat,
                            outputFigures = outputFigures,
                            plt4peak = TRUE, ###
                            count_lim = count_lim,
                            Figure_name = heatmap_name,
                            Figure.width=Figure.width,
                            Figure.height=Figure.height,
                            cytoBand=cytoBandFile,
                            global_norm=global.norm) #Figure_name = paste0(FigureName_prefix,".heatmap_subClust")
  cell_anno_new <- clust_res$cellMeta
  cell_anno$subCluster <- cell_anno_new[rownames(cell_anno),"subCluster"]
  subClusters <- sort(unique(cell_anno$subCluster))
  
  if(segRatio_on_bin){
    mtx4ratio <- clust_res$binCount.nom
    mtx4ratio <- mtx4ratio[,rownames(cell_anno_new)]
  }else{
    mtx4ratio <- mtx_qcn_lim[,rownames(cell_anno_new)]
    mtx4ratio <-  t(t(mtx4ratio)*mean(colSums(mtx4ratio))/colSums(mtx4ratio))
  }
  
  if(!is.null(chrs)){
    mtx4ratio <- mtx4ratio[grepl(paste0("^",chrs,delimit,collapse = "|"),rownames(mtx4ratio)),,drop=F]
  }
  
  
  message(" (4) Calculate sub-cluster level CN ratio.")
  SegScore_total_subC <- CNratioInfer(input_matrix=mtx4ratio,
                                      inputcellano=cell_anno_new,
                                      cytoBand=cytoBandFile,
                                      outdir=outdir_sub,
                                      group2bulkby="subCluster",
                                      min_cells_in_group = min_cells_in_group,
                                      normallabel="reference",  #"reference",
                                      cells.ref = ref_cells,
                                      FiltSeg= FiltSeg,
                                      SegLen_min=SegLen_min,
                                      SegSize_min=SegSize_min,
                                      seg_method=seg_method,
                                      segValue_method=segValue_method,
                                      seg.count.lim = seg.count.lim,
                                      penalty = penalty,
                                      outFigure = outputFigures,
                                      outFigureName=outFigureName,
                                      color_dot=CorlorDot4Ratio,
                                      output_ref.value = FALSE,
                                      outRDS=FALSE,
                                      RDSname = SegScore_total_subCName)

  if(optimalClust){
    SegRatio <- SegScore_total_subC
    cells_ref <- rownames(cell_anno)[cell_anno[,1] %in% ref_group_names]
    
    
    if(is.null(cnv.ref.df)){stop("Please provide cnv.ref.df  with columns with c('ratio_map','integerCNV') for optimal clustering.")}
   
    Rn=1
    while(Rn<=10){

      CNV_res_subCj <- SegRatio$seg_score_binLevel_filt
  
      plts <- unlist(lapply(1:length(CNV_res_subCj),function(x,CNV_res_subCj,cnv.ref.df,chrs){
        clusti <- names(CNV_res_subCj)[x]
        p <- plot_ratio2CNV.chr(CNV_res_subCj[[x]],cnv.ref.df,chrs=chrs,name.data=clusti,doPlot=FALSE)
        return(p)
      },CNV_res_subCj,cnv.ref.df,chrs))
      
      clust_select <- names(CNV_res_subCj)[which(plts >mse.min)]
      clust.tb <- table(cell_anno_new$subCluster)
      clust_select_smallclt <- names(clust.tb)[clust.tb<min_cells_in_group &names(clust.tb)!="reference"]
      
      cells_select <- rownames(cell_anno_new)[cell_anno_new$subCluster %in% c(clust_select,clust_select_smallclt)]

      #re-name the subCluster in cell_anno
      fixRow_idx <- which(!cell_anno$subCluster %in% c(clust_select,"reference"))
      
      cell_anno$subCluster[fixRow_idx] <- paste0(Rn,".",cell_anno$subCluster)[fixRow_idx]
      
      if(length(cells_select)>=min_cells_in_group){
        mtx_filt <- mtx_qcn_lim[,c(cells_select,cells_ref)]
        
        cell_anno_new <- binClust(inputMat=mtx_filt,cellMeta=cell_anno[c(cells_select,cells_ref),,drop=F],
                              ref_group_names=ref_group_names,
                              chrs=chrs,
                              delimit=delimit,
                              onRatio = ClustOnRatio,
                              outdir=outdir_sub,
                              window = binSize,
                              seu_resol = seu_resolution,
                              OutMat =OutBinMat,
                              outputFigures = F)
        cell_anno[rownames(cell_anno_new),"subCluster"] <- cell_anno_new$subCluster
        
        SegRatio <- CNratioInfer(input_matrix=mtx_filt[,rownames(cell_anno_new)],
                                            inputcellano=cell_anno_new,
                                            cytoBand=cytoBandFile,
                                            outdir=outdir_sub,
                                            group2bulkby="subCluster",
                                            min_cells_in_group = min_cells_in_group,
                                            normallabel="reference",
                                            FiltSeg= FiltSeg,
                                            SegLen_min=SegLen_min,
                                            seg_method=seg_method,
                                            penalty = penalty,
                                            outFigureName="Round1.pdf",
                                            color_dot=CorlorDot4Ratio,
                                            output_ref.value = FALSE,
                                            outRDS=F)
      
        Rn=Rn+1
      }else{break}

  }
  } 
  
  return(list(SegScore_total_subC=SegScore_total_subC,
              cell_anno_new=cell_anno_new,
              binCount=clust_res$binCount,
              binCount.norm = clust_res$binCount.nom,
              binRatio=clust_res$binRatio))
}


#' @title find_nearst()
#' @export

find_nearst <- function(df_seg_C1,initialCN,align_by="SegMean"){
  CNs <- initialCN$CN
  Ratio_tar <- initialCN$ratio
  Ratio_obs <- df_seg_C1[,align_by]
  dis_mt <- do.call(cbind,lapply(1:length(Ratio_tar),function(x,Ratio_obs,Ratio_tar){
    d= abs(Ratio_obs-Ratio_tar[x])
    return(d)
  },Ratio_obs,Ratio_tar))
  rownames(dis_mt) <- rownames(df_seg_C1)
  indice <- which(rowSums(!is.na(dis_mt))>0)
  min_indices <- unlist(apply(dis_mt, 1, which.min))
  df_seg_C1$ratio_map[indice] <- Ratio_tar[min_indices]
  df_seg_C1$integerCNV[indice] <- CNs[min_indices]
  return(df_seg_C1)
}


#' @title CNV.align()
#' @description align CNVs based on the baseline
#' @param rato.correct subclone ratio correction
#' @export
CNV.align <- function(df_seg_C1,initialCN,subclone.rawBase=NULL,rato.correct=TRUE,align_by = "SegMean"){
  colnames(initialCN)[grepl("ratio",colnames(initialCN),ignore.case = T)] <- "ratio"
  colnames(initialCN)[grepl("CN",colnames(initialCN),ignore.case = T)] <- "CN"
  initialCN <- initialCN[order(initialCN$ratio),,drop=F]

   if(any(diff(initialCN$CN)>1)){
    initialCN <- CNbaseline.fill(initialCN)
   }
  
  df_seg_C1 <- find_nearst(df_seg_C1,initialCN,align_by)
  mse <- mse_loss(df_seg_C1)
  
  if(rato.correct){
    if(all(!grepl("base",colnames(initialCN)))){
      stop("Run 'peakIndex' for initialCN first!")
    }
    cloneBase <- initialCN$ratio[initialCN$base==1]
    cloneBase.cn <- initialCN$CN[initialCN$base==1]
    #find peak of subclone
    pp <- density(df_seg_C1[,align_by],width=0.1)
    dens <- pp$y
    ratioX <- pp$x
    subcloneBase <- ratioX[which.max(dens)]
    if(!is.null(subclone.rawBase)){
      subcloneBase.cn <- subclone.rawBase$CN[which(subclone.rawBase$base==1)]
      diff.cn <- subcloneBase.cn - cloneBase.cn
    }else{diff.cn=0}

    diff = subcloneBase-cloneBase

    df_seg_C2 <- df_seg_C1
    if(diff.cn <=1 & (abs(diff) > 0.01)){
      df_seg_C2[,align_by] <- df_seg_C2[,align_by]-diff
      df_seg_C2 <- find_nearst(df_seg_C2,initialCN,align_by)
      mse.corrected <- mse_loss(df_seg_C2)
      if(mse.corrected< mse){
        df_seg_C1=df_seg_C2
      }
    }
  }
  
  return(df_seg_C1)
}



#' @title combine_res_chr()
#' @description combine the segment ratio result from each chromosome
#' @return matrix of peak-cell seg.ratio combined different cluster of cells from each chromosome
#' @param data.dir contains chr folders
#' @param mtx_seg NA matrix with colname with cell and rowname with peak
#' @export
combine_res_chr <- function(data.dir, mtx_seg,chrs,df_meta_slim,
                            groups_col="subclone",
                            outdir=".",outFile=TRUE,estract_by = "SegMean"){
  if(!file.exists(outdir)){dir.create(outdir,recursive=T)}
  subclones <- unique(df_meta_slim$subclone)
  ###(1)create original peak-level segment ratio list of subCluster
  seg_res_ls <- lapply(subclones,function(x){
    seg <- c()
  })
  names(seg_res_ls) <-subclones 
  mtx_seg_total <- c()
  
  for(chrj in chrs){
    outdir_j <- paste0(data.dir,"/",chrj)
    SegScore_total <- readRDS(paste0(outdir_j,"/segRatio_",chrj,".rds"))
    binIDs <- do.call(rbind,lapply(SegScore_total,function(x)return(x[,c("binID","Chromosome","Start","End")])))
    binIDs <- unique(binIDs)
    binIDs <- binIDs[binIDs$Chromosome %in% chrj,,drop=F]
    
    binIDs <- binIDs[order(binIDs$Start),,drop=F]
    
    mtx_seg_chr <- matrix(nrow=nrow(binIDs),ncol = ncol(mtx_seg))
    rownames(mtx_seg_chr) <- binIDs$binID
    colnames(mtx_seg_chr) <- colnames(mtx_seg)
    
    csvfile <- paste0(outdir_j,"/cell_info_by",chrj,".csv")
    df_subclone <- read.csv(csvfile,row.names = 1)
    df_subclone <- df_subclone[!df_subclone$subCluster%in%"reference",,drop=F]
    clusterset <- as.character(unique(df_subclone[,groups_col]))

    for(i in 1:length(clusterset)){
      groupi <- as.character(clusterset[i])
      cells_i <- rownames(df_subclone)[df_subclone[,groups_col] %in% groupi]
      # cnv_i <- CNV_res[[groupi]]
      # cnv_i_chr <- cnv_i[grepl(paste0(chrj,"-"),cnv_i$binID),,drop=F]
      # 
      # intCNV <- cnv_i_chr$integerCNV
      # mt_cnvi <- matrix(rep(intCNV,times=length(cells_i)),ncol=length(cells_i))
      # mtx_cnv[cnv_i_chr$binID,cells_i] <- mt_cnvi
      
      seg_df_bin <- SegScore_total[[groupi]]
      seg_df_bin <- seg_df_bin[grepl(paste0(chrj,"-"),seg_df_bin$binID),,drop=F]
      cols_keep <- c("binID","Chromosome","Start","End","segName","length_seg","binRatio","SegMean","ratio_map","integerCNV") #,"chrom"
      seg_df_bin <- seg_df_bin[,cols_keep]
      
      segmean<- seg_df_bin[,estract_by]
      mt_segmean <- matrix(rep(segmean,times=length(cells_i)),ncol=length(cells_i))
      mtx_seg_chr[seg_df_bin$binID,cells_i] <- mt_segmean

      #update segment ratio
      sub_label <- paste(chrj,groupi,sep=".")
      subclone_idx_2fill <- as.character(df_meta_slim$subclone[df_meta_slim[,chrj] %in% sub_label])
      seg_res_ls[subclone_idx_2fill] <- lapply(subclone_idx_2fill,function(x,seg_res_ls,seg_df_bin){
        seg_res_ls[[x]] <- rbind(seg_res_ls[[x]],seg_df_bin)
      },seg_res_ls,seg_df_bin)
    }
    
    mtx_seg_total <- rbind(mtx_seg_total,mtx_seg_chr)
  }
  if(outFile){
    saveRDS(seg_res_ls,paste0(outdir,"/SegRatio_binLevel_combinChr.rds"))
    saveRDS(mtx_seg_total,paste0(outdir,"/segMtx_combinChr.",estract_by,".rds"))
  }

  return(list(seg_res_ls=seg_res_ls,mtx_seg=mtx_seg_total))
}

#' @param res.ls list of resuts of multiple clusters
#' @description the names of res.ls should in first column of metadata
ls.colummn2mtx <- function(res.ls,metadata,extract_by="integerCNV"){
  clusters <- names(res.ls)
  metadata <- metadata[metadata[,1] %in% clusters,,drop=F]
  bins <- unique(unlist(lapply(res.ls, function(x)as.character(x$binID))))
  bins.df <- data.frame(peak=bins)
  chromInfo <- separate(bins.df, peak, into = c("seqnames", "start","end"), sep = "-")
  bins.df$seqnames <- gsub("chr","",chromInfo$seqnames)
  bins.df$seqnames[bins.df$seqnames=="X"] <- 23
  bins.df$seqnames[bins.df$seqnames=="Y"] <- 24
  bins.df$seqnames <- as.numeric(bins.df$seqnames)
  bins.df$start <- as.numeric(chromInfo$start)
  bins.df$end <- as.numeric(chromInfo$end)
  bins.df <- bins.df[order(bins.df$seqnames,bins.df$start),]
  mtx <- matrix(nrow=nrow(bins.df),ncol = nrow(metadata))
  rownames(mtx) <- bins.df$peak
  colnames(mtx) <- rownames(metadata)
  for(clt in clusters){
    cells_select <- rownames(metadata)[metadata[,1] %in% clt]
    res_df <- res.ls[[clt]]
    rows_select <- res_df$binID
    mtx[rows_select,cells_select] <- res_df[,extract_by]
  }
  #Fill NA with neighbors
  mtx2 <- apply(mtx,2,function(x){
    x <- fill_na_with_previous(x)
    x <-  fill_na_with_next(x)
  })
  
 # mtx2 <- na.omit(mtx)
  return(mtx2)
}


module_scCNesti <- function(mtx_bin,cellMeta,cells_obs,cells_ref,clonal_res,
                            cloneName_ref= NULL,
                            group_by="subCluster",
                            cytoBandFile,
                            numCores = 4,
                            outdir="./"){
  if(!is.null(cloneName_ref)){
    flog.info("Calculate single-cell level ratio relative to reference.") 
  }
 if(length(clonal_res)==0){
   print("No clonal CNV successful estimated in 'clonal_res'.")
   return(NA)
   }
  #calculate correction factor
  #scale peak count by non-zero percentage
  Nzero <- rowSums(mtx_bin[,cells_obs]==0)
  prop.zero <- Nzero/ncol(mtx_bin[,cells_obs])
  subClusts <- sort(unique(cellMeta[cells_obs,group_by]))
  if(length(subClusts)<numCores){
    numCores <- length(subClusts)
  }
  cl <- makeCluster(numCores)
  registerDoParallel(cl)  
  segRatio.ls <- foreach(k=1:length(subClusts),.combine='c',
                         .packages = c("futile.logger", "dplyr"),
                         .export = c("CNratioInfer.sc","CNratio_peak","pseudo_bulk_v2",".get_ref_mean","seg_score")) %dopar% {

    subc <- as.character(subClusts[k])
    #print(subc)
    if(subc %in% names(clonal_res)){
      CNVres_df <- clonal_res[[subc]]$input_BinSegRatio
    }else if(!is.null(cloneName_ref)){
      CNVres_df <- clonal_res[[cloneName_ref]]$input_BinSegRatio
    }else{
      CNVres_df <- clonal_res[[1]]$input_BinSegRatio
    }
    
    cells_obs_sub <- rownames(cellMeta)[cellMeta[,group_by] %in% subc]
    
    mtx_bin_in <- mtx_bin[rownames(CNVres_df),,drop=F]
    prop.zero_in <- prop.zero[rownames(CNVres_df)]
    
    SegScore_ls <- list(seg_score_binLevel=list())
    for(j in 1:length(cells_obs_sub)){
      SegScore_sc <- CNratioInfer.sc(input_matrix= mtx_bin_in,
                                     cell.name =cells_obs_sub[j],
                                     inputcellano = cellMeta,
                                     cytoBandFile = cytoBandFile,
                                     group2bulkby=group_by,
                                     normallabel="reference",
                                     prop_zero = prop.zero_in,
                                     seg.ref = CNVres_df,
                                     seg_method="robseg",
                                     outPlot = F,
                                     outFigureName.prefix="segRatio_sc",
                                     outdir = outdir)
      df_seg_C1 <- SegScore_sc$seg_score_binLevel[[1]]
      SegScore_ls$seg_score_binLevel[[cells_obs_sub[j]]] <- df_seg_C1
    }
    
    list(SegScore_ls$seg_score_binLevel)
  }
  stopCluster(cl)
  names(segRatio.ls) <- subClusts
  
  cellbinRatio <- cellRatio_corrected(segRatio.ls,clonal_res,refClone = cloneName_ref)
  cellbinCN <- cellintegerCN(segRatio.ls,clonal_res,refClone = cloneName_ref)  
  return(list(cellbinCN=cellbinCN,cellbinRatio=cellbinRatio))
}



#' @title clone_similarity()
#' @export

clone_similarity <- function(cellbinCN,cellMeta,group_by="subCluster",Csize.min=20){
  library(MixGHD)
  bincount <- apply(cellbinCN, 2, function(x){length(x[is.na(x)])})
  cellbinCN2 <- cellbinCN[,bincount==0]
  #pheatmap::pheatmap(cellbinCN2,cluster_cols = F,cluster_rows = T,show_rownames = F,show_colnames = F,annotation_row = cell_anno_new)
  set.seed(123)
  optK <- findccloneNum(cellbinCN2)
  while(1){
    set.seed(123)
    kcluster <- kmeans(cellbinCN2, optK, iter.max = 10, nstart = 1)
    #Check the size of each cluster and make adjustments
    cluster_sizes <- table(kcluster$cluster)
    small_clusters <- names(cluster_sizes[cluster_sizes < Csize.min])
    keep_clusters <-  as.numeric(names(cluster_sizes[cluster_sizes >= Csize.min]))
    if(optK==1 | length(keep_clusters)>1){
      break
    }else{
      optK <- optK-1
    }
  }

  if (length(small_clusters) > 0 &length(keep_clusters)>0) {
    for( k in small_clusters){
      
      cluster_indices <- which(as.character(kcluster$cluster) == k)
      distances <- sapply(keep_clusters, function(center) sum((kcluster$centers[center,] - kcluster$centers[as.numeric(k),])^2))
      nearest_cluster <- keep_clusters[which.min(distances)]
      kcluster$cluster[cluster_indices] <- nearest_cluster
      
      cluster_indices2 <- which(kcluster$cluster == nearest_cluster)
      kcluster$centers[nearest_cluster,] <- colMeans(cellbinCN2[cluster_indices2,])

    }
    
    #re-name cluster
    unique_sorted_clusters <- sort(unique(kcluster$cluster))
    new_clusters <- match(kcluster$cluster, unique_sorted_clusters)
    kcluster$cluster <- new_clusters
  }
  
  
  
  #table(kcluster$cluster)
  #pheatmap::pheatmap(cellbinCN2[row.names(kclusterres),],cluster_cols = F,cluster_rows = F,show_rownames = F,show_colnames = F,annotation_row = kclusterres)
  
  cluster1 <- cellMeta[row.names(cellbinCN2),group_by]
  cluster2 <- as.character(kcluster$cluster)
  score <- ARI(cluster1, cluster2)
  cellMeta_new <-cellMeta[row.names(cellbinCN2),group_by,drop=F]
  cellMeta_new[,group_by] <-  cluster2
  return(list(score=score,new_cluster=cellMeta_new))
}


#' @title clonal_ploidy()
#' @description add clonal ploidy and diploidy annotation to clonal results
#' @export
clonal_ploidy <- function(clone_res_update, baseCN_frac.min = 0.95){
  clone_ploidy <- unlist(lapply(clone_res_update,function(x)x$ploidy))
  if(length(clone_res_update)>2){
    clone_ploidy_mean <- median(clone_ploidy)
  }else{
    clone_ploidy_mean <- NULL
  }
  
  clone_res_update <-  lapply(clone_res_update,function(x,clone_ploidy_mean){
    cl_score <- clone_scoring(x,ploidy.ref=clone_ploidy_mean)
    x$score <- cl_score
    seg_dat <- x$seg.dat
    CNest <- x$CNest
    #CNest <- peakIndex(x$input_BinSegRatio,CNest)
    baseCN <- CNest$CN[CNest$base==1]
    baseCN_frac <- sum(seg_dat$w[seg_dat$integerCN==baseCN])/sum(seg_dat$w)
    if(baseCN_frac>baseCN_frac.min){
      x$diploidy <- TRUE
    }else{x$diploidy <- FALSE}
    x$score$baseCN_frac <- baseCN_frac
    return(x)
  },clone_ploidy_mean)
  return(clone_res_update)
}


#' @title celloutput()
#' @param cells.obs the input cells
#' @param cluster1 group annotation for cells.obs
#' @param cluster2 updated group annotation for cells.obs
#' @export
celloutput <- function(cells.obs,cluster1,cluster2,clone_res,mtx_bin,
                       CNest.ref,
                       min_cells_in_group=20,
                       penalty.lim=c(0.5,1.5),
                       seg_method="PELT",
                       FiltSeg=TRUE,
                       SegLen_min=2e6,
                       SegSize_min = 10,
                       seg.count.lim = 80,
                       cytoBandFile = NULL,
                       outdir="./",
                       minCN.frac=0.01,
                       seg_dat_ref,
                       correct_by_dist=FALSE,
                       delt_lim = 0.1,
                       segValue_method="median",
                       p.adj_cutoff = 0.05,
                       ratio_diff_cutoff=NULL,
                       diploidy_cutoff = 0.95,
                       Diploid.domin=TRUE,
                       ...){
  clone_res_raw <- clone_res
  cellinfo <- data.frame(cellname =cells.obs)
  row.names(cellinfo) <- cells.obs
  # if(length(unique(cluster1))==length(unique(cluster2))){
  #   cellinfo$clone <- cluster1
  # }else{
    cellinfo$clone <- cluster2
    ##Repeat step 1, step 2, step 3 according to the new clones
    cells_ref <- colnames(mtx_bin)[!colnames(mtx_bin) %in% cells.obs]
    cellinfo_ref <- data.frame(row.names=cells_ref,cellname =cells_ref)
    cellinfo_ref$clone <- "reference"
    cell_meta <- rbind(cellinfo,cellinfo_ref)
    clone_ratio <- CNratioInfer(input_matrix=mtx_bin[,rownames(cell_meta)],
                                        inputcellano=cell_meta,
                                        cytoBand=cytoBandFile,
                                        outdir=outdir,
                                        group2bulkby="clone",
                                        min_cells_in_group = min_cells_in_group,
                                        normallabel="reference",  #"reference",
                                        FiltSeg= FiltSeg,
                                        SegLen_min=SegLen_min,
                                        SegSize_min=SegSize_min,
                                        seg_method=seg_method,
                                        segValue_method=segValue_method, #default is "density"
                                        seg.count.lim = seg.count.lim,
                                        penalty = penalty.lim,
                                        outFigure=FALSE,
                                        output_ref.value = FALSE,
                                        outRDS=FALSE)



    CNV_clone.ini <- CNV_esti_ini(outdir=outdir,
                                     segScore_file=clone_ratio,
                                     filt_seg = FiltSeg,
                                     length_seg_cutoff =SegLen_min,
                                     segValue_method=segValue_method,
                                    outputFigure=F)
    if(!all(unlist(lapply(CNV_clone.ini,is.na)))){
      
      clone_res <- ploidyRefine(CNV_clone.ini,delt.lim =delt_lim,minCN.frac=minCN.frac)
      clone_res <- ploidyRefine.ref(clone_res,CNest.ref,minCN.frac,seg_dat_ref,Diploid.domin=Diploid.domin)
      
      #clone_res <- refine_segCN.Bayes(clone_res,seg_dat_ref)
      if(correct_by_dist){
        clone_res <- refine_segCN.dist_v2(clone_res)
      }

    
    }else{
      clone_res <- clone_res_raw
      cellinfo$clone <- cluster1
    }
   
  # }
  
  clone_res <-  clonal_ploidy(clone_res,baseCN_frac.min=diploidy_cutoff)
  clone_res_update <- segCNV_refine(
        clonalRes=clone_res,
        CNest.ref=CNest.ref,
        seg_method=seg_method,
        seg.count.lim = seg.count.lim,
        penalty = penalty.lim,
        segValue_method=segValue_method,
        FiltSeg= FiltSeg,
        SegLen_min=SegLen_min,
        SegSize_min = SegSize_min,
        p.adj_cutoff = p.adj_cutoff,
        ratio_diff_cutoff=ratio_diff_cutoff)

    clone_res_update <-  clonal_ploidy(clone_res_update,baseCN_frac.min=diploidy_cutoff)
  
  
  binID <- rownames(mtx_bin)
  binID <- do.call(rbind,strsplit(binID,split="-|:|_"))
  binsize <- as.numeric(binID[,3])-as.numeric(binID[,2])
  # cellploidy <- apply(cellbinCN2, 1, function(x,binsize){
  #   sum(x*binsize,na.rm=TRUE)/sum(binsize,na.rm=TRUE)
  # },binsize)
  # index <- match(cellinfo$cellname,row.names(cellbinCN2))
  # cellinfo$ploidy <- cellploidy[index]
  # cellinfo <- cellinfo[order(cellinfo$clone,cellinfo$ploidy),]
  cellinfo <- cellinfo[order(cellinfo$clone),]
  
  # cellbinCN2 <- cellbinCN2[cellinfo$cellname,]

  binCount_norm <- mtx_bin[,cellinfo$cellname]
  
  outputres <- list(#cellCN = cellbinCN2,
                    cellinfo = cellinfo,
                    clonalest = clone_res_update,
                    cellbinCount.norm=binCount_norm)
}


#' @title clone_scoring()
#' @param cloneRes the output list of ploidyRefine(), i.e. clonal_res
#' @export
clone_scoring <- function(cloneRes,scaleFactor=1,y_true.term="ratio",y_pred.term="relativeCN",
                          CNbase =cloneRes$CNest,ploidy.ref=NULL){
  seg.dat <- cloneRes$seg.dat
  score_mse <- mse_loss(seg.dat,y_true.term=y_true.term,y_pred.term=y_pred.term)
  
  df_seg_C1 <- cloneRes$input_BinSegRatio
  df_seg_C1$segLen <- df_seg_C1$End - df_seg_C1$Start
  
  CNbase <- CNbase[order(CNbase$CN),]
  initialCN1 <- peakIndex(df_seg_C1,CNbase)
  initialCN1 <- na.omit(initialCN1)
  Seg_Frac <- SegFrac(df_seg_C1,initialCN1)
  delt <- (initialCN1$ratio[2]-initialCN1$ratio[1])/(initialCN1$CN[2]-initialCN1$CN[1])
  param_delt <- ifelse(delt<=0.3,delt,1)
  
  Nseg <- nrow(cloneRes$seg.dat)
  if(!is.null(ploidy.ref)){
    ploidy <- as.numeric(cloneRes$ploidy)
    ploidy_dis <- abs(ploidy-ploidy.ref)
  }else{ploidy_dis=0}
  

  score <- scaleFactor*Seg_Frac*(-log(score_mse))*param_delt/(log1p(Nseg)*exp(ploidy_dis))
  score_df <- data.frame(score=score,mse=score_mse,Frac=Seg_Frac,Nseg=Nseg,ploidy_diff=ploidy_dis,delt=delt)
  return(score_df)
}
