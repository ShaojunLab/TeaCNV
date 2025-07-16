options(expressions=10000)
workpath=system.file(package = "TeaCNV")
setwd(workpath)
suppressPackageStartupMessages({
  library(Matrix)
  library(Signac)
  library(Seurat)
  library(tidyr)
  library(ggplot2)
  library(futile.logger)
  library(gridExtra)
  library(MatrixGenerics)
  library(MixGHD)
  library(foreach)
  library(doParallel)
  library(plyranges)
  library(ggExtra)
  library(stringr)
  library(dplyr)
  library(R.utils)
})
# source("./R/funs_filtMT.R",chdir = T)
# source("./R/mydataProcess.R",chdir = T) 
# source("./R/funs_atac2cnv.R",chdir = T)
# source("./R/funs_binClust.R",chdir = T)
# source("./R/seg_plot.R",chdir = T)
# source("./R/figurePlot.R",chdir = T)
# source("./R/seg_score.R",chdir = T)
# source("./R/fun_pseudoBulk_mat.R",chdir = T)
# source("./R/funs_find_subclones.R",chdir = T)
# source("./R/funs_utils.R",chdir = T)


# #CNV esti
# source("./R/findmode.R",chdir = T)
# source("./R/optEstimation.R",chdir = T)
# source("./R/results.check.R",chdir = T)
# source("./R/funs_bulkCNV_estimation.R",chdir = T)
# source("./R/ploidy.correction.v2.R",chdir = T)
# source("./R/genome.subclonality.function.R",chdir = T)
# source("./R/fun_segmentation_comb.R",chdir = T)
# source("./R/funs_InformationSegs.r")

TeaCNV <- methods::setClass("TeaCNV",
                               slots = c(data.rawpeak="ANY",
                                         data.reference = "ANY",
                                         data.filt="ANY",
                                         data.lim="ANY",
                                         data.binRatio = "ANY",
                                         data.binCount = "ANY",
                                         data.binCount.norm ="ANY",
                                         cell_anno.filt="data.frame",
                                         ref.value="ANY",
                                         plots.filt = "ANY",
                                         options = "list",
                                         logVar="ANY"))

#' @title CreateTeaCNVObject()
#' @param FiltCell_by Filtering method if FiltCell=TRUE, c("Zscore","density","quantile"),default is "Zscore".
#' @param Correct_by_length Set to FALSE if the input matrix feature bins are of equal length
#' @export
CreateTeaCNVObject <- function(input,
                                  annotationFile,
                                  ref_group_names,
                                  assay = "ATAC",
                                  delim = "\t",
                                  genome = "hg38",
                                  blacklistFile = NULL,
                                  cytoBandFile = NULL,
                                  ChrRemove = c('chrX', 'chrY', 'chrM'),
                                  FiltCell = TRUE,
                                  FiltCell_by = "quantile",
                                  CellnCount_quant=c(0.05,0.95),
                                  CellnFeature_quant=c(0.05,0.95),
                                  cellproplim = 0.05,
                                  count_lim = 4,
                                  prop_reset=TRUE,
                                  Correct_by_length=TRUE
){
  packages <- c("stringr", "dplyr", "tidyr", "Matrix", "irlba","plyranges","futile.logger")
  invisible(lapply(packages, require, character.only = TRUE))
  
  logVar <- ""
  log_appender <- function(level, ...) {
    logVar <<- paste0(logVar, level)
  }
  flog.appender(log_appender)
  
  ## input data
  if (Reduce("|", is(input) == "character")) {
    flog.info(sprintf("Loading count matrix: %s", input))
    if (substr(input, nchar(input)-2, nchar(input)) == ".gz") {
      raw.data <- read.table(connection <- gzfile(input, 'rt'), sep=delim, header=TRUE, row.names=1, check.names=FALSE)
      close(connection)
      raw.data <- as.matrix(raw.data)
    }else if(substr(input, nchar(input)-3, nchar(input)) == ".rds") {
      raw.data <- readRDS(input)
    }else {
      raw.data <- read.table(input, sep=delim, header=TRUE, row.names=1, check.names=FALSE)
      raw.data <- as.matrix(raw.data)
    }
  }else if(Reduce("|", is(input) %in% c("dgCMatrix", "matrix"))){
    raw.data <- as.matrix(input)
  }else if (Reduce("|", is(input) %in% c("data.frame"))) {
    raw.data <- as.matrix(input)
  }else if(Reduce("|", is(input) %in% c("Seurat"))){
    raw.data <- as.matrix(input[[assay]]@counts)
  } else {
    stop("Can not find input as a matrix, data.frame, Seurat object, or filename")
  }
  
  ##annotations file
  if (Reduce("|", is(annotationFile) == "character")) {
    flog.info(sprintf("Loading cell annotations file: %s", annotationFile))
    cell_anno <- read.table(annotationFile, header=FALSE, row.names=1, sep=delim, stringsAsFactors=FALSE, colClasses = c('character', 'character'))
  }else if (Reduce("|", is(annotationFile) %in% c("dgCMatrix", "matrix", "data.frame"))) {
    cell_anno <- annotationFile
  }else {
    stop("Can not find annotationFile as a matrix, data.frame, or filename")
  }
  
  if (rownames(cell_anno)[1] == "V1") {
    cell_anno = cell_anno[-1, , drop=FALSE]
  }
  
  if(is.null(blacklistFile)){
    if(genome %in% "hg38"){
      blacklistFile <- system.file("data", "hg38-blacklist.v2.bed", package = "TeaCNV")
    }else if(genome %in% "hg19"){
      blacklistFile <- system.file("data", "hg19-blacklist.v2.bed", package = "TeaCNV")
    }
  }
  
  if(is.null(cytoBandFile)){
    if(genome %in% "hg38"){
      cytoBandFile <- system.file("data", "cytoBand_hg38.tsv", package = "TeaCNV")
    }else if(genome %in% "hg19"){
      cytoBandFile <- system.file( "data", "cytoBand_hg19.txt", package = "TeaCNV")
    }
    
  }
  
  flog.info("\n")
  flog.info("START: Filtering cells by nCount and nFeature.")
  plts <- list()
  cells_t <- rownames(cell_anno)[!cell_anno[,1]%in%ref_group_names]
  cells_n <- rownames(cell_anno)[(cell_anno[,1]%in%ref_group_names)]
  data.t <- raw.data[,cells_t]
  data.n <- raw.data[,cells_n]
  
  if(FiltCell){
    by_zscore.x = by_density = by_quantile = FALSE
    if(FiltCell_by=="Zscore"){
      by_zscore.x=TRUE
    }else if(FiltCell_by=="density"){
      by_density = TRUE
    }else if(FiltCell_by=="quantile"){
      by_quantile = TRUE
    }
    filt.t <- FiltCell.obj(object=data.t,
                           assay=assay,
                           by_density = by_density,
                           by_quantile = by_quantile,
                           by_zscore.x = by_zscore.x,
                           nCount_quantile_lim = CellnCount_quant,
                           nFeature_quantile_lim=CellnFeature_quant,
                           log_ = TRUE)
    mtx_filt.t <- data.t[,filt.t$cell_keep]
    plts[["obs"]] <- filt.t$plots$total
    
    filt.n <- FiltCell.obj(object=data.n,
                           assay=assay,
                           by_density = by_density,
                           by_quantile = by_quantile,
                           by_zscore.x = by_zscore.x,
                           nCount_quantile_lim = CellnCount_quant,
                           nFeature_quantile_lim=CellnFeature_quant,
                           log_ = TRUE)
    mtx_filt.n <- data.n[,filt.n$cell_keep]
    plts[["ref"]] <- filt.n$plots$total
    mtx_filt <- cbind(mtx_filt.t,mtx_filt.n)
  }else{
    mtx_filt <- raw.data
  }


  mtx_filt <- blacklist.rm(mtx_filt,blacklistFile)
  
  cell_anno_filt <- cell_anno[colnames(mtx_filt),,drop=F]
  if(ncol(mtx_filt)>10){
    mtx_filt <- mtx_filt[rowSums(mtx_filt)>10,]
  }else{
    mtx_filt <- mtx_filt[rowSums(mtx_filt)>1,]
  }
 
  mtx_qc_filt <- FiltPeak(mtx=mtx_filt,
                          cell_anno=cell_anno_filt,
                          ref_group_names=ref_group_names,
                          cellprop = cellproplim,
                          outheatmap = FALSE,
                          blacklist_file = blacklistFile,
                          ChrRemove = ChrRemove,
                          prop_reset=prop_reset)
  
  n_raw_peaks <- nrow(raw.data)
  n_remove_pk <- n_raw_peaks - nrow(mtx_qc_filt)
  flog.info(sprintf("Keep %g peaks, removing %g %% of peaks",nrow(mtx_qc_filt), n_remove_pk/n_raw_peaks * 100))
  flog.info("Correct count by the length of peak to count per kb.")
  if(Correct_by_length){
    mtx_qc_filt <- Correct_by_length(mtx_qc_filt)
  }
  
  flog.info(paste0("The maximum count is ",round(max(mtx_qc_filt)),"."))
  flog.info(paste0("Set maximum boundary for count to ",count_lim,"."))
  #count_lim <- round(quantile(mtx_qc_filt,0.99));
  mtx_qc_filt[mtx_qc_filt>count_lim] <- count_lim
  
  object <- new(
    Class = "TeaCNV",
    data.rawpeak= mtx_filt,
    data.reference = NULL,
    data.filt = mtx_qc_filt,
    data.lim = NULL,
    data.binRatio = NULL,
    data.binCount = NULL,
    data.binCount.norm = NULL,
    cell_anno.filt=cell_anno_filt,
    ref.value = NULL,
    plots.filt = plts,
    options = list("ref_group_names"=ref_group_names,
                   "genome"=genome,
                   "cytoBandFile"=cytoBandFile,
                   "count_lim"=count_lim,
                   "CurrentStep"=0),
    logVar = logVar
  )
  flog.info("\n")
  
  return(object)
}

#' @title EffectCorrectByRef()
#' @param cell_anno data.frame of cell annotation, containing the cell name(row.name), 
#' the cell group classification(first column) and the sample ID (second column), tab-delimited
#' @param NormalTypeList cell groups in the first column of cell_anno
#' @param sampleID_normal Normal samples in the second column of cell_anno
#' @export
EffectCorrectByRef <- function(mat,cell_anno,sampleID_normal,NormalTypeList,sampleIDs_to_correct,Ncells.min=10){
  library(Matrix) 
  library(progress)
  invisible(gc())
  #mat <- as.matrix(mat)
  cell_anno <- cell_anno[colnames(mat),,drop=F] #match cell order

  get_factor <- function(data,cell_meta,sample.in,NormalTypeList){
    scaleFactor <- c()
    for(i in 1:length(NormalTypeList)){
      type_i <- NormalTypeList[i]
      cells_ref4cor_i <- rownames(cell_meta)[cell_meta[,2] %in% sample.in & cell_meta[,1] %in% type_i]
     if(length(cells_ref4cor_i)>1){
       mt_ref4cor_i <- data[,cells_ref4cor_i,drop=F]
       sf <- mean(colSums(mt_ref4cor_i,na.rm = TRUE),na.rm = TRUE)
     }else{sf=NA}
      scaleFactor <- c(scaleFactor,sf)
    }
    names(scaleFactor) <- NormalTypeList
    return(scaleFactor)
  }
  
  cell_scaleFactor <- get_factor(mat,cell_anno,sampleID_normal,NormalTypeList)

  pb <- progress_bar$new(total = length(sampleIDs_to_correct), format = "[:bar] :percent :elapsed")
  
  # mat_cor <- matrix(NA,nrow = nrow(mat),ncol = ncol(mat))
  mat_cor <-sparseMatrix(i = integer(0),j = integer(0), dims = c(nrow(mat), ncol(mat)),x = numeric(0))
  rownames(mat_cor) <- rownames(mat)
  colnames(mat_cor) <- colnames(mat) 
  tag <- c()
  cell_used <- c()
  
  for(j in sampleIDs_to_correct){
    cells_j  <- rownames(cell_anno)[cell_anno[,2]%in%j]
    cells_n <- rownames(cell_anno)[cell_anno[,1] %in% NormalTypeList & cell_anno[,2]%in%j]
    N_cells_in_type <- table(cell_anno[cells_n,1])
    N_cells_in_type <- N_cells_in_type[N_cells_in_type>=Ncells.min]
    #flog.info(paste0(length(cells_n), " normal cells in sample ",j,"."))
    if(length(cells_n)>=Ncells.min){
      cell_scaleFactor_j <- get_factor(mat,cell_anno,j,names(N_cells_in_type))
      cell_scaleFactor_j <-cell_scaleFactor_j[!is.na(cell_scaleFactor_j)]
      if(length(cell_scaleFactor_j)>0){
        total_factors <- cell_scaleFactor[names(cell_scaleFactor_j)]/cell_scaleFactor_j
        total_factors_mean <- mean(total_factors,na.rm = TRUE)
        #cell_scaleFactor_j <- mean(colSums(mat[,cells_n,drop=F],na.rm = TRUE),na.rm = TRUE)
        mat_cor[,cells_j] <- sweep(mat[, cells_j, drop = FALSE], 2, total_factors_mean, `*`)
        tag <- cbind(tag,"corrected")
      }

    }else{
      flog.warn(paste0("Number of normal cells found in sample ",j," is less than 10, the original matrix returned!"))
      tag <- cbind(tag,"raw")
      mat_cor[,cells_j] <- mat[,cells_j]
    }
    cell_used <- c(cell_used,cells_j)
    pb$tick()
  }
  # cell_left <- rownames(cell_anno)[!rownames(cell_anno) %in% cell_used]
  # cell_left <- unique(c(cell_left,cells_ref4cor))
  na_index <- which(
    vapply(
      seq_len(ncol(mat_cor)), 
      function(j) all(is.na(mat_cor[, j])), 
      logical(1)
    )
  )
  mat_cor[,na_index] <- mat[,na_index]
  mat_cor <- as(mat_cor, "sparseMatrix") 
  return(list(matrix=mat_cor,tag=tag))
  invisible(gc())
}



#' @title runTeaCNV()
#' 
#' @param input_obj  the TeaCNV object or matrix of peaks (rows) vs. cells (columns) containing the raw counts
#'                    It can be a filename or the data directly.
#' @param annotationFile annotation of cells, tab-delimited.   
#' @param seg_method method for segmentation, including "robseg" and "gaussian".
#' @param penalty paramter if seg_method is "gaussian"
#' @scFactor the scale factor for scFactor*SegSize_min, which is the final minimum size of segmentation. (Default is 1)
#' @export
runTeaCNV <- function(
    input_obj,
    outdir=NULL,
    output_obj_name = "TeaCNV.obj",
    #cell filtering settings
    cellFilt_ZeroProp =0.6,
    #batch effect correction settings
    correctTn5Batch = FALSE,
    corrected_group_by=NULL,
    sampleID_normal=NULL,
    batchTag = "raw",
    global.norm=TRUE,
    # segmentation setting
    binSize = 1,
    cytoBandFile = NULL,
    seg_method = "PELT",
    segValue.method="median",
    penalty.lim = c(0.5,1.5),
    seg.count.lim = 80,
    CorlorDot4Ratio = TRUE,
    FiltSeg= TRUE,
    SegLen_min = 2e6,
    SegSize_min = 20,
    #find sub-clone settings
    Zscore_cutoff = 1.28,
    ClustOnRatio = TRUE,
    OutBinMat =TRUE,
    seu_resolution = 1,
    outputFigures = TRUE,
    min_cells_in_group = 20,
    ChrRemove = c('chrX', 'chrY', 'chrM'),
    numCores = 6,
    #Clone optimization setting
    minCN_frac = 0.01,
    delt_lim =0.25,
    StopStep = 4,
    choice = "",
    scFactor=1,
    diploidy_pct.max=0.95,
    p.adj_CloneMerge=0.1
){
  packages <- c("Signac", "Seurat", "tidyr", "Matrix", "irlba","plyranges","futile.logger","ggplot2","MatrixGenerics","MixGHD","stringr","data.table","ComplexHeatmap","RColorBrewer","changepoint")
  invisible(lapply(packages, require, character.only = TRUE))
  
  if (is.null(outdir)) {
    flog.error("Error, outdir is NULL, please provide a path.")
    stop("outdir is NULL, please provide a path.")
  }
  if(outdir != "." && !file.exists(outdir)){
    flog.info(paste0("Creating output path ", outdir))
    dir.create(outdir,recursive=T)
  }
  flog.appender(appender.file(paste0(outdir,"/LogInfo.log")))
  flog.info(input_obj@logVar)
  
  startTime <- Sys.time()
  CurrentStep <- as.numeric(input_obj@options$CurrentStep)
  
  input_obj@options$cellFilt_ZeroProp = cellFilt_ZeroProp
  input_obj@options$seg_method = seg_method
  input_obj@options$penalty.lim = penalty.lim
  input_obj@options$seg.count.lim = seg.count.lim
  input_obj@options$CorlorDot4Ratio = CorlorDot4Ratio
  input_obj@options$FiltSeg= FiltSeg
  input_obj@options$SegLen_min = SegLen_min
  input_obj@options$ClustOnRatio = ClustOnRatio
  input_obj@options$OutBinMat =OutBinMat
  input_obj@options$seu_resolution = seu_resolution
  input_obj@options$outputFigures = outputFigures
  input_obj@options$min_cells_in_group = min_cells_in_group
  input_obj@options$numCores = numCores
  input_obj@options$delt_lim = delt_lim
  input_obj@options$segValue.method = segValue.method
  input_obj@options$minCN.frac = minCN_frac
  
  ref_group_names <- input_obj@options$ref_group_names
  count_lim <- input_obj@options$count_lim
  cytoBandFile <- input_obj@options$cytoBandFile
  # plots_filtCells <- cowplot::plot_grid(plotlist = input_obj@plots.filt,
  #                                       labels = paste0(input_obj@plots.filt,"."),
  #                                       ncol = 2,align="hv")
  # ggsave(paste0(outdir, "/inputCells_filtered.pdf"),plots_filtCells, width=12, height=4,device = pdf,bg="white")
  # 

  mtx_qc_filt <- input_obj@data.filt
  cell_anno <- input_obj@cell_anno.filt

  runstep = 1
  mtx_filt_File <- paste0(outdir,"/step01.mtx_filt.rds")
  cell_anno_updateFile <- paste0(outdir,"/step01.cellAnno.csv")
  if(CurrentStep < runstep){
    flog.info("\n\n")
    flog.info("Step 1- Data process: continue filtering cells.")
    flog.info("Filtering cells by cell score.")
   
    cell_anno$group <- ifelse(cell_anno[,1]%in%ref_group_names,"reference","observed")
    cells_ref <-  rownames(cell_anno)[cell_anno[,1]%in%ref_group_names]
    cells_obs <- rownames(cell_anno)[cell_anno$group %in% "observed"]
    mtx_obs <- mtx_qc_filt[,cells_obs]
    cellScore <- cell_score(data = mtx_obs,ZeroProp.max=cellFilt_ZeroProp,outdir = outdir,outFigure = TRUE,max.value=count_lim) #,cell_meta=cell_anno[cells_obs,,drop=F]
    cells_used <- cellScore$cell[cellScore$label %in% "Low"]
    cells_filered <- cellScore$cell[!cellScore$label %in% "Low"]
    matrix_filt <- mtx_qc_filt[,c(cells_used,cells_ref),drop=F]
    cell_anno$group[rownames(cell_anno) %in%cells_filered ] <- "removed"
    write.csv(cell_anno,cell_anno_updateFile,row.names = T)
    
    cell_anno_filt <- cell_anno[!cell_anno$group %in% "removed",,drop=F]
    
    input_obj@cell_anno.filt <- cell_anno_filt
    input_obj@data.filt <- matrix_filt
    #Step 2:
    if(correctTn5Batch){
      flog.info("Varies Tn5 enzyme effect correction for filtered matrix.")
      if(!is.null(corrected_group_by)){
        samples_to_correct <- unique(cell_anno_filt[,corrected_group_by])
        if(is.null(NormalTypeList)){NormalTypeList=ref_group_names}
        #NormalTypeList is reference groups in the first column of cell_anno_filt
        #samples_to_correct is sampleIDs in the second column of cell_anno_filt
        correc_res<- EffectCorrectByRef(matrix_filt,cell_anno_filt,sampleID_normal=sampleID_normal,NormalTypeList=NormalTypeList,sampleIDs_to_correct=samples_to_correct)
        mtx_qc_filt2 <- correc_res$matrix
        batchTag <- correc_res$tag
      }else{
        stop(paste0("No 'corrected_group_by' provided for batch effect correction."))
      }
      mtx_qcn_lim <- as.matrix(mtx_qc_filt2)
 
    }else{mtx_qcn_lim <- as.matrix(matrix_filt) }
    
    ##3. set Max boundary
    flog.info(sprintf("Setting max of peek count boundary: %g",count_lim))
    #bound_up <- round(quantile(mtx_qcn_lim,0.99));bound_up
    mtx_qcn_lim[mtx_qcn_lim>count_lim] <- count_lim
    
    saveRDS(mtx_qcn_lim,mtx_filt_File)
  
    input_obj@data.lim <- mtx_qcn_lim
    input_obj@options$CurrentStep <- runstep
    CurrentStep <- runstep
    saveRDS(input_obj,paste0(outdir,"/",output_obj_name))
    
  }else{ 
      if(!is.null(input_obj@data.lim)){
      mtx_qcn_lim <- input_obj@data.lim
    }
    cell_anno_filt <- input_obj@cell_anno.filt
  }
  
  if(!exists('mtx_qcn_lim') & file.exists(mtx_filt_File)){
    mtx_qcn_lim <- readRDS(mtx_filt_File)
  }else if(!exists('mtx_qcn_lim')){
    stop(paste0("No ",mtx_filt_File," or ",output_obj_name," found in ",outdir))
  } 
  if((!exists("cell_anno_filt")) & file.exists(cell_anno_updateFile)){
    cell_anno <- read.csv(cell_anno_updateFile,row.names = 1)
    cell_anno_filt <- cell_anno[!cell_anno$group %in% "removed",,drop=F]
    flog.info("Load cell_anno_filt.")
  }else if(!exists("cell_anno_filt")){
    stop(paste0(cell_anno_updateFile," not found in ",outdir))
  }

  flog.info("\n\n")

  N_peaks <- dim(mtx_qcn_lim)[1]
  if(is.null(binSize)){
    binSize <- ceiling(N_peaks/2000)
  }
  if(N_peaks>1e4 &N_peaks<1.5e4){
    SegSize_min= ifelse(binSize<6,round(100/binSize),20)
  }else if(N_peaks>=1.5e4){
    SegSize_min= ifelse(binSize<6,round(120/binSize),20)
  }
   
  input_obj@options$binSize <- binSize
  input_obj@options$SegSize_min = SegSize_min

  ####step2:Initial subclone determination, CN ratio calculation, segmentation
  runstep=runstep+1

  flog.info(sprintf("Number of peaks is: %g, binSize is: %g",N_peaks,binSize))

  if(StopStep>=runstep){
    outdir_sub <- paste0(outdir,"/subCluster");if(!file.exists(outdir_sub)){dir.create(outdir_sub,recursive=T)}
    cell_anno_new_file <- paste0(outdir_sub,"/cell_info_subClust.csv")
    CNV_res_ls_file <- paste0(outdir_sub,"/step02.Initilal_subClust.rds")
    
    if(CurrentStep < runstep){
      flog.info("\n\n")
      flog.info("Step 2: Initial subclone determination, CN ratio calculation and segmentation.")
      CNV_res_ls <-  pcaClustCN(mtx_qcn_lim,
                 cell_anno_filt,
                 ref_group_names,
                 cytoBandFile=cytoBandFile,
                 ClustOnRatio=ClustOnRatio,
                 seu_resolution=seu_resolution,
                 min_cells_in_group = min_cells_in_group,
                 binSize = binSize,
                 outputFigures=TRUE,
                 OutBinMat =OutBinMat,
                 heatmap_name = paste0("heatmap_subClust"),
                 count_lim = count_lim,
                 segRatio_on_bin=TRUE,
                 seg_method = seg_method,
                 segValue_method=segValue.method,
                 seg.count.lim = seg.count.lim,
                 SegSize_min = SegSize_min,
                 penalty = penalty.lim,
                 CorlorDot4Ratio = TRUE,
                 FiltSeg= FiltSeg,
                 SegLen_min = SegLen_min,
                 outdir=outdir_sub,
                 global.norm=global.norm)
      SegScore_total_subC <- CNV_res_ls$SegScore_total_subC
      
      saveRDS(CNV_res_ls,CNV_res_ls_file)
      
      input_obj@data.binRatio <- CNV_res_ls$binRatio
      input_obj@data.binCount <- CNV_res_ls$binCount
      input_obj@data.binCount.norm <- CNV_res_ls$binCount.norm
      
      cell_anno_new <- CNV_res_ls$cell_anno_new
      write.csv(cell_anno_new,cell_anno_new_file,row.names = T)
      while (length(dev.list()) > 0) {
        dev.off()
      }
      
      input_obj@cell_anno.filt <- cell_anno_new
      input_obj@ref.value <- SegScore_total_subC$ref.value
      input_obj@options$CurrentStep <- runstep
      CurrentStep <- runstep
      saveRDS(input_obj,paste0(outdir,"/",output_obj_name))

    }else if(file.exists(CNV_res_ls_file)){
      CNV_res_ls <- readRDS(CNV_res_ls_file)
      SegScore_total_subC <- CNV_res_ls$SegScore_total_subC
    }else{
      stop(paste0("No ",CNV_res_ls_file," found."))
    }
    

    bins <- rownames(input_obj@data.binCount.norm)
    start <- sapply(strsplit(bins, "_|-|:"), "[", 2)
    end <- sapply(strsplit(bins, "_|-|:"), "[", 3)
    length <- as.numeric(end)-as.numeric(start)
    meanLn <- median(length) 
    flog.info("\n")
    flog.info(sprintf("The median length of each bin is: %g",meanLn))
    
    

    if(!exists("cell_anno_new")){
      if("cell_anno.filt" %in% slotNames(input_obj)){
        cell_anno_new <- input_obj@cell_anno.filt
      }else if(file.exists(cell_anno_new_file)){
        cell_anno_new <- read.csv(cell_anno_new_file,row.names=1)
      }else{
        stop(paste0("No ",cell_anno_new_file," found."))
      }
    }
    left_anno_cols <- list()
    color_r <- suppressWarnings(get_group_color_palette("Set3")(length(unique(cell_anno_new$subCluster))))
    names(color_r) <- sort(unique(cell_anno_new$subCluster))
    left_anno_cols[["subCluster"]] <- color_r
  }

  flog.info("\n")
  
  #Step 3: clone-level integer copy numbers estimation
  runstep=runstep+1
  if(StopStep>=runstep){
    CNVresFile_subC <- paste0(outdir_sub,"/step03.CNV_res_clonal.rds")
    CNVresFile_subC_update <- paste0(outdir_sub,"/step03.CNV_res_clonal_update.rds")
    if(CurrentStep < runstep){
      flog.info("\n\n")
      flog.info("Step 3: Clone-level integer copy numbers estimation.")
      
      ###3.1 initial estimate
      CNV_res_subC.ini <- CNV_esti_ini(outdir=outdir_sub,
                               segScore_file=SegScore_total_subC,
                               filt_seg = TRUE,
                               length_seg_cutoff =1e6,
                               segValue_method=segValue.method,
                               genome="hg38",
                               true_loc = TRUE,
                               outputFigure = FALSE,
                               outplot_name=paste0("CNV_subClust_initial"))
      saveRDS(CNV_res_subC.ini,paste0(outdir_sub,"/CNV_res_subC.ini.rds"))
      if(all(unlist(lapply(CNV_res_subC.ini,is.na)))){
        message("Failed to estimate CNV for any clone, try filtering the data using a stricter threshold!")
        return(0)
      }
      ###3.2 ploidyRefine
      clonal_res <- ploidyRefine(CNV_res_subC.ini,delt.lim =delt_lim )

      clone_ploidy <- unlist(lapply(clonal_res,function(x)x$ploidy))
      if(length(clonal_res)>2){
        clone_ploidy_mean <- median(clone_ploidy)
      }else{
        clone_ploidy_mean <- NULL
      }
     
      ###3.3 scoring clones
      clonal_res <-  lapply(clonal_res,function(x,clone_ploidy_mean){
        cl_score <- clone_scoring(x,ploidy.ref=clone_ploidy_mean)
        x$score <- cl_score
        seg_dat <- x$seg.dat
        seg_dat <- seg_dat[!is.na(seg_dat$integerCN),,drop=F]
        CNest <- x$CNest
        baseCN <- CNest$CN[CNest$base==1]
        baseCN_frac <- sum(seg_dat$w[seg_dat$integerCN==baseCN])/sum(seg_dat$w)
        if(baseCN_frac>diploidy_pct.max){
          x$diploidy <- TRUE
        }else{x$diploidy <- FALSE}
        x$score$baseCN_frac <- baseCN_frac
        return(x)
      },clone_ploidy_mean)
      
      #unlist(lapply(clonal_res,function(x)x$score$score))
      saveRDS(clonal_res,CNVresFile_subC)
      
      #3.4 clones paired-comparison
      cell_anno_cl <- cell_anno_new[,c("subCluster"),drop=FALSE]
      cell_anno_cl <- cell_anno_cl[!cell_anno_cl$subCluster %in% "reference",,drop=F]
      clones <- sort(unique(cell_anno_cl$subCluster))
      Nclone <- length(clones)
      outdir_clt2 <- paste0(outdir_sub,"/segRatio_compare")
      ifelse(!dir.exists(file.path(outdir_clt2)), dir.create(file.path(outdir_clt2)), FALSE)
      mtx_bin <- input_obj@data.binCount.norm
      
      if(Nclone>1){
        seg_ls <- list()
        for(col_ref in clones){
          SegScore_total <- seg4clusters(mtx_bin[,rownames(cell_anno_cl)],cell_anno_cl,
                                         cytoBand =cytoBandFile,
                                         col_ref,
                                         doFilt=T,
                                         seg.penalty =c(0.5,1),
                                         seg.method=seg_method,
                                         segSize.min = SegSize_min,
                                         outplot_path=outdir_clt2,
                                         outFigure=FALSE,
                                         outname=paste0("segRatio_clone_vsC",col_ref,".pdf"))
          
          seg_ls[[col_ref]] <- SegScore_total$seg_score_binLevel
          
        }
        saveRDS(seg_ls,paste0(outdir_sub,"/clonal_RelativeRatio_segment.rds"))
        cluster_new <- suppressMessages(mergeClones(ratiodata=clonal_res,segdata =seg_ls,doPlot=FALSE,outdir=outdir,Zscore.cutoff=Zscore_cutoff,p.adj.cutoff=p.adj_CloneMerge))

      }else{
         cluster_new <- data.frame(subCluster=cell_anno_new$subCluster,clone_merged="1")
      }
      
      #3.5 merge clones
      cluster_rm <- clones[!clones %in% names(clonal_res)]
      cells_rm <- rownames(cell_anno_cl)[cell_anno_cl$subCluster %in% cluster_rm]
      
      colnames(cluster_new) <- c("subCluster","clone_merged")
      cell_anno_new$row <- rownames(cell_anno_new)
      cell_anno_new <- cell_anno_new[,!grepl("clone_merged",colnames(cell_anno_new)),drop=F]
      cell_anno_new <- left_join(cell_anno_new,cluster_new,by="subCluster")
      rownames(cell_anno_new) <- cell_anno_new$row
      cell_anno_new$clone_merged[rownames(cell_anno_new) %in% cells_rm] <- "removed"
      cell_anno_new$clone_merged[is.na(cell_anno_new$clone_merged)] <- "reference"
      
      input_obj@cell_anno.filt <- cell_anno_new
      write.csv(cell_anno_new,cell_anno_new_file,row.names = T)

           
      clone_scores <- unlist(lapply(clonal_res,function(x)x$score$score))
      diploidy_cl <- unlist(lapply(clonal_res,function(x)x$diploidy))
      diploidy_indx <- which(diploidy_cl)
      if(length(diploidy_indx)>0){
        clone_scores <- clone_scores[-diploidy_indx]
      }
      
      best_clone <- names(clone_scores)[which.max(clone_scores)]
      cat("\nThe subgroup scores are:\n")
      cat(clone_scores)
      cat("\n")
      cat(paste0("\nThe best subgroup is C",as.character(best_clone)))
      input_obj@options$bestClone <- best_clone 
      input_obj@options$CurrentStep <- runstep
      CurrentStep <- runstep
      saveRDS(input_obj,paste0(outdir,"/",output_obj_name))
     
      ###3.5 visualization
      plotdir <- paste0(outdir,"/Figures")
      ifelse(!dir.exists(file.path(plotdir)), dir.create(file.path(plotdir)), FALSE)

      flog.info("\n")
      flog.info(paste0("Output figures of initial subgroup-level CN estimation to ", outdir))
      #plot segment (initial)
      plot_combine_seg(clonal_res,outplot_name="subgroup_CN_initial.pdf",show_dots=TRUE,
                       outdir=plotdir)
      

    }else if(file.exists(CNVresFile_subC)){
      clonal_res <- readRDS(CNVresFile_subC)
    }else{
      stop(paste0("No ",CNVresFile_subC," found."))
    }

    cat("\nSet up one subgroup as the reference to estimate CNVs based on figure 'subgroup_CN_initial.pdf'.\nNOTE:Do NOT selcet a diploidy subgroup!")
    # cat("\n\nWhether to manually select the optimal subgroup as the reference ? (input [the number after 'C' as the clone name] or [Enter] for automatic selection:)\n")
    
    choice <- ""
    # tryCatch({
    #   choice <- withTimeout({
    #     as.character(readline(prompt = "best_clone: "))
    #   }, timeout = 5, onTimeout = "silent")
    # }, TimeoutException = function(ex) {
    #   message("No input detected within 5 seconds. Continuing...")
    # })
    
    ###
    if(choice!=""){
    	#choice <- as.character(input_obj@options$bestClone)  
    	
    }else{
      choice <- as.character(input_obj@options$bestClone) 
      
    }
    if(choice!=best_clone){
      best_clone <- choice
      input_obj@options$bestClone <- best_clone
      saveRDS(input_obj,paste0(outdir,"/",output_obj_name))
    }
    
    ###
    flog.info(paste0("\nThe reference subgroup is: ",choice))
    #
    flog.info(paste0("\nUpdate all subgroup-level CNVs based on the CNVs of reference subgroup..."))
    clone.ref <- choice
    CNest.ref <- clonal_res[[clone.ref]]$CNest
    CNest.ref <- peakIndex(clonal_res[[clone.ref]]$input_BinSegRatio,CNest.ref,index_col="SegMean")
    ploidy.ref <- clonal_res[[clone.ref]]$ploidy
    seg_dat_ref <- clonal_res[[clone.ref]]$seg.dat
    
    
    clonal_res_updated <- ploidyRefine.ref(clonal_res,CNest.ref,minCN.frac=minCN_frac)

    clone_ploidy <- unlist(lapply(clonal_res_updated,function(x)x$ploidy))
    if(length(clonal_res_updated)>2){
      clone_ploidy_mean <- median(clone_ploidy)
    }else{
      clone_ploidy_mean <- NULL
    }
    
    ###scoring clones
    clonal_res_updated <-  lapply(clonal_res_updated,function(x,clone_ploidy_mean){
      cl_score <- clone_scoring(x,ploidy.ref=clone_ploidy_mean)
      x$score <- cl_score
      seg_dat <- x$seg.dat
      seg_dat <- seg_dat[!is.na(seg_dat$integerCN),,drop=F]
      CNest <- x$CNest
      baseCN <- CNest$CN[CNest$base==1]
      baseCN_frac <- sum(seg_dat$w[seg_dat$integerCN==baseCN])/sum(seg_dat$w)
      if(baseCN_frac>0.95){
        x$diploidy <- TRUE
      }else{x$diploidy <- FALSE}
      x$score$baseCN_frac <- baseCN_frac
      return(x)
    },clone_ploidy_mean)
    
    saveRDS(clonal_res_updated,CNVresFile_subC_update)
    
    
    
    flog.info("\n")
    # flog.info(paste0("Output figures of initial clone-level CN estimation to ", outdir))
    # #plot segment (initial-updated)
    # plot_combine_seg(clonal_res_updated,
    #                  outplot_name="clonalCN_initial_update.pdf",show_dots=TRUE,
    #                  outdir=outdir)
    # plot_combine_seg(clonal_res_updated,
    #                  outplot_name="clonalCN_initial_update_noDots.pdf",show_dots=FALSE,
    #                  outdir=outdir)
    clonal_res <- clonal_res_updated
  }
  
  
  #Step 4: Single-cell level integer copy number estimation
  runstep=runstep+1
  if(StopStep>=runstep){
    CNVresFile_sc <- paste0(outdir,"/final.CNVres.rds")
    #4.1 initial scCNV
    cells_obs <- rownames(cell_anno_new)[cell_anno_new$group %in% "observed"]
    cells_ref <-  rownames(cell_anno_new)[cell_anno_new$group%in%"reference"]
    mtx_bin <- input_obj@data.binCount.norm

    if(CurrentStep < runstep){
      flog.info("\n\n")
      flog.info("Step 4: Final clonal absolute copy number estimation. Iterative optimization of clone structures.")
      
      ###re-estimate clonal CNVs for merged clones
      cluster1 <- cell_anno_new[cells_obs,"subCluster"]
      cluster2 <- cell_anno_new[cells_obs,"clone_merged"]
      outres <- celloutput(cells_obs,cluster1,cluster2,clone_res=clonal_res,mtx_bin=mtx_bin,
                           CNest.ref= CNest.ref,
                           min_cells_in_group,
                           penalty.lim,
                           seg_method,
                           FiltSeg,
                           SegLen_min,
                           SegSize_min*scFactor,
                           seg.count.lim,
                           cytoBandFile,
                           outdir=outdir,
                           minCN.frac=minCN_frac,
                           seg_dat_ref=seg_dat_ref,
                           segValue_method=segValue.method)
   
      cellbinCount <- input_obj@data.binCount[,outres$cellinfo$cellname]
      outres$cellbinCount <- cellbinCount
      outres$cellbinRatio_raw <- input_obj@data.binRatio
      outres$cellinfo$clone <- factor(outres$cellinfo$clone,levels=sort(unique(as.numeric(outres$cellinfo$clone))))
      outres$CNest.ref <- CNest.ref


      # #plot segment (final)
      # plot_combine_seg(outres$clonalest,ylim=NULL,
      #                  outplot_name="clonalCNV_final_refine0.pdf",show_dots=TRUE,
      #                  outdir=outdir)
      
      
      #
      clonal_res_final <- outres$clonalest
      clonal_res_final <- lapply(clonal_res_final,function(x){
        clonal_res <- x$input_BinSegRatio
        clonal_res <- clonal_res[,!grepl("relativeCN|integerCN",colnames(clonal_res))]
        x$input_BinSegRatio<- clonal_res
        return(x)
      })

      clonal_res_final <-  lapply(clonal_res_final,function(x){
        seg_dat <- x$seg.dat
        seg_dat <- seg_dat[!is.na(seg_dat$integerCN),,drop=F]
        CNest <- x$CNest
        baseCN <- 2
        baseCN_frac <- sum(seg_dat$w[seg_dat$integerCN==baseCN])/sum(seg_dat$w)
        aneuploidy_score <- 1- baseCN_frac
        x$score$aneuploidy_score <- aneuploidy_score
        return(x)
      })
      outres$clonalest <- clonal_res_final
      aneuploidy_score <- unlist(lapply(outres$clonalest,function(x){
        x$score$aneuploidy_score
      }))
      clone_anno <- data.frame(clone_merged=names(outres$clonalest),aneuploidy_score=aneuploidy_score)
      cell_anno_new <- left_join(cell_anno_new,clone_anno,by=intersect(colnames(cell_anno_new),colnames(clone_anno)))
      rownames(cell_anno_new) <- cell_anno_new$row
      input_obj@cell_anno.filt <- cell_anno_new
      write.csv(cell_anno_new,cell_anno_new_file,row.names = T)
      

      sm_win <- ifelse(binSize<6,ceiling(60/binSize),10)
      plt_mt2 <- smooth_and_denoise(input_obj@data.binRatio,window=sm_win)
      outres$cellbinRatio_deNoise <- plt_mt2
      
      saveRDS(outres,CNVresFile_sc)
      saveRDS(input_obj,paste0(outdir,"/",output_obj_name))
      
      
      input_obj@options$CurrentStep <- runstep
      CurrentStep <- runstep
     
    }else if(file.exists(CNVresFile_sc)){
      outres <- readRDS(CNVresFile_sc)
      plt_mt2 <- outres$cellbinRatio_final
      

    }
  }

   
  #Step 5: Visualization
  #segment clonal level Ratio snd CNVs
  plot_combine_seg(outres$clonalest,
                   outplot_name="clonalCN_final.pdf",show_dots=TRUE,
                   outdir=plotdir)  
  plot_combine_seg(outres$clonalest,
                   outplot_name="clonalCN_final_noDots.pdf",show_dots=FALSE,
                   outdir=plotdir) 

  #segment clonal CNVs
    p_ls <- lapply(names(outres$clonalest),function(cluster,outdir){
      clonal_res <- outres$clonalest[[cluster]]$input_BinSegRatio
      integerCNV <- outres$clonalest[[cluster]]$seg.dat
      score <-  round(outres$clonalest[[cluster]]$score$score,2)
      ymax <- quantile(integerCNV$integerCN,0.99)
      plt <- cellLineSeg_v3(binRatio.df=clonal_res,integerCNV.df=integerCNV,outdir=plotdir,
                            ylab="CNVs",ylim=c(0,ymax),
                            plot_dots = FALSE,
                            plot_hist = FALSE,
                            color_seg_gradient=T, seg.color = "#0d3b66",
                            color_hist_gradient=T,hist_color= NULL,
                            fname=paste0("C",cluster,": score=",score), 
                            height = 2,width=10,outPlot=FALSE,
                            color_limit = c(1,6),
                            value.bin="integerCN",
                            value.segment="integerCN",
                            label_CN=FALSE)
      return(plt$p1+ theme(plot.margin = unit(c(0.5, 0.5, 0, 0.5),'lines')))
    },outdir)
    names(p_ls) <- names(outres$clonalest)
    pcom <- cowplot::plot_grid(plotlist=p_ls,ncol = 1,align="v")
    ggsave(paste0(plotdir, "/integerCN.pdf"),pcom, width=10, height=2*length(p_ls),device = pdf,bg="white")

    flog.info("\n")
    flog.info(paste0("Output figures of raw single-level relative CN (Ratio) to ", outdir))

    cellMeta <- cell_anno_new[!cell_anno_new$subCluster %in%"reference",,drop=F]
    cellMeta  <- cellMeta[!cellMeta$clone_merged %in%"removed",,drop=F]
    cellMeta$subCluster <- factor(cellMeta$subCluster,levels=sort(unique(as.numeric(cellMeta$subCluster))))
    cellMeta$clone_merged <- factor(cellMeta$clone_merged,levels=sort(unique(as.numeric(cellMeta$clone_merged))))
    cellMeta <- cellMeta[order(cellMeta$clone_merged,cellMeta$subCluster),,drop=F]
    colnames(cellMeta) <- gsub("clone_merged","clone",colnames(cellMeta))
    cellMeta <- cellMeta %>%
      add_count(clone) %>%
      as.data.frame()
    colnames(cellMeta)[colnames(cellMeta)=="n"] <- "Ncells"
    rownames(cellMeta)<- cellMeta$row
    cellMeta$CellProportion <- cellMeta$Ncells/nrow(cellMeta)
    cloneInfo_table <- unique(cellMeta[,c("clone","aneuploidy_score","Ncells","CellProportion")])
    write.csv(cloneInfo_table,paste0(outdir,"/cloneInfo_table.csv"),row.names=F)
    
    left_anno_cols <- list()
    cellmeta2 <- outres$cellinfo[,c("clone"),drop=F]
    cellmeta2 <- na.omit(cellmeta2)
    color_r2 <- suppressWarnings(get_group_color_palette("Set3")(length(unique(cellmeta2$clone))))
    names(color_r2) <- sort(unique(as.numeric(cellmeta2$clone)))
    left_anno_cols[['clone']] <- color_r2
    

    p_scRatio <- heatmap4peakMt(mat=plt_mt2[,rownames(cellmeta2)],
                                meta_info=cellmeta2,
                                sep_by="-",
                                outdir= plotdir,value.type="ratio",
                                clust_rows=F,clustering_method_rows = "ward.D2",
                                show_legend_row = T,
                                fileout_name=paste0("heatmap_CNratio"),
                                col_list=left_anno_cols,
                                width=10,height=5,device="pdf")

    ###heatmap: clonal CNVs
    CNmt <- do.call(cbind,lapply(names(outres$clonalest),function(cluster){
      clonal_res <- outres$clonalest[[cluster]]$input_BinSegRatio
      integerCNV <- outres$clonalest[[cluster]]$seg.dat
      integerCNV <- integerCNV[,c("segName","relativeCN","integerCN")]
      clonal_res <- left_join(clonal_res,integerCNV,by="segName")
      rownames(clonal_res) <- clonal_res$binID
      return(clonal_res$integerCN)
      }))
    rownames(CNmt) <- rownames(outres$clonalest[[1]]$input_BinSegRatio)
    colnames(CNmt) <- names(outres$clonalest)
    
    # Ncol <- round(cloneInfo_table$CellProportion*100)
    # new_CNmt <- do.call(cbind, lapply(1:ncol(CNmt), function(x) {
    #   matrix(rep(CNmt[, i], Ncol[x]), ncol = Ncol[x])
    # }))
    # rownames(new_CNmt) <- rownames(CNmt)
    # colnames(new_CNmt) <- seq(1:ncol(new_CNmt))
    
    # clone_info <- data.frame(row.names=colnames(new_CNmt),clone=rep(colnames(CNmt), Ncol))
    # color_r <- suppressWarnings(get_group_color_palette("Set3")(length(unique(clone_info$clone))))
    # names(color_r) <- sort(unique(as.numeric(clone_info$clone)))
    # left_anno_cols[["clone"]] <- color_r
    # p_cloneCN <- heatmap4peakMt(mat=new_CNmt,
    #                       meta_info=clone_info,
    #                       sep_by="-",
    #                       outdir= plotdir,value.type="CNV",
    #                       clust_rows=F,clustering_method_rows = "ward.D2",
    #                       show_legend_row = T,
    #                       fileout_name=paste0("heatmap_cloneCNV-1_SizeScale"),
    #                       col_list=left_anno_cols,
    #                       width=10,height=3,device="pdf")

    clone_info <- data.frame(row.names=colnames(CNmt),clone=colnames(CNmt))
    color_r <- suppressWarnings(get_group_color_palette("Set3")(length(unique(clone_info$clone))))
    names(color_r) <- sort(unique(clone_info$clone))
    left_anno_cols <- list()
    left_anno_cols[["clone"]] <- color_r
    height <- ifelse(ncol(CNmt)>2,0.35*(ncol(CNmt))+0.5,ifelse(ncol(CNmt)==1,1.2,1.5))
    p_cloneCN <- heatmap4peakMt(mat=CNmt,
                                meta_info=clone_info,
                                sep_by="-",
                                outdir= plotdir,
                                value.type="CNV",
                                clust_rows=F,
                                show_legend_row = T,
                                legend_titles="integer CN",
                                fileout_name=paste0("heatmap_cloneCNV"),
                                col_list=left_anno_cols,
                                column.title = NULL,
                                width=10,height=height,device="pdf") #height=5

    #plot combined segmentation and CNV
    pl <- c() 
    pl_data <- c()
    ymax <- quantile(unlist(lapply(outres$clonalest,function(x){x$seg.dat$integerCN})),0.99)
    for(cluster in sort(names(outres$clonalest))){
        clonal_res <- outres$clonalest[[cluster]]$input_BinSegRatio
        clonal_res <- clonal_res[,!grepl("relativeCN|integerCN",colnames(clonal_res))]
        integerCNV <- outres$clonalest[[cluster]]$seg.dat
        clonal_res <- merge(clonal_res,integerCNV[,c("relativeCN", "integerCN","segName")],by="segName",all.x=T)
         pl_scRatio <- seg_plot(clonal_res,name.data="",
                        value.bin="integerCN",
                        value.segment="integerCN",
                        ylab="CNVs",
                        plot_dots = F,
                        ylim=c(0,ymax),
                        add_yline = 2,
                        plot_seg = T,
                        seg_col = left_anno_cols[["clone"]][cluster],
                        plotDir=plotdir,outPlot=F,
                        color_hist=T,plot_hist=T,
                        color_hist_gradient=F,
                        plot_colors=left_anno_cols[["clone"]][cluster],
                        device = "pdf",disconnected=F)
       p1 <- pl_scRatio$p1
       pl[[cluster]] <- p1
       pl_data[[cluster]] <- pl_scRatio$p1$data
    }
    p1 <- pl[[1]]
    if(length(pl)>1){
        for(pl_i in 2:length(pl_data)){
            cluster <- names(pl_data)[pl_i]
            a <- pl_data[[cluster]]
            dat2<- data.frame(x=a$segEnd[1:(length(a$segEnd)-1)], y=a$value.segment[1:(length(a$segEnd)-1)],
                              xend=a$segStart[2:length(a$segEnd)], yend=a$value.segment[2:length(a$segEnd)])
            p1 <- p1 +
                geom_segment(
                data = dat2,
                mapping = aes(x=x, y=y,
                                xend=xend, yend=yend), 
                colour=left_anno_cols[["clone"]][cluster],na.rm =T,size = 1,
                inherit.aes = FALSE)
        }
    }
    ggsave(paste0(plotdir,"/integerCN_combined.pdf"),p1,height = 2,width=10)

  ###clonal paired-comparison
    cell_anno_cl <- outres$cellinfo[,c("clone"),drop=FALSE]
    cell_anno_cl <- cell_anno_cl[!is.na(cell_anno_cl$clone),,drop=F]
    clones <- as.character(sort(as.numeric(as.character(unique(cell_anno_cl$clone)))))
    Nclone <- length(clones)
    outdir_clt4 <- paste0(outdir,"/clone_segRatio_compare")
    ifelse(!dir.exists(file.path(outdir_clt4)), dir.create(file.path(outdir_clt4)), FALSE)
    mtx_bin <- input_obj@data.binCount.norm
    
    if(Nclone>1){
      seg_ls <- list()
      for(col_ref in clones){
        SegScore_total <- seg4clusters(mtx_bin[,rownames(cell_anno_cl)],cell_anno_cl,
                                       cytoBand =cytoBandFile,
                                       col_ref,
                                       doFilt=T,
                                       seg.penalty = c(0.5,1),
                                       seg.method = seg_method,
                                       segSize.min = SegSize_min,
                                       outplot_path=outdir_clt4,
                                       outname=paste0("segRatio_cloneX_vsC",col_ref,".pdf"))
        
        seg_ls[[col_ref]] <- SegScore_total$seg_score_binLevel
        
      }
      # saveRDS(seg_ls,paste0(outdir,"/clonal_RelativeRatio_segment.rds"))
    }

  
  endTime <- Sys.time()
  duration <- endTime - startTime
  flog.info("\n\n")
  flog.info(paste0("Program running time: ",duration))
  flog.info("\n\n")
  flog.info("All analysis done!")
}



