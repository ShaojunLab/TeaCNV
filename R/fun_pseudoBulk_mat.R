#' @title gen_pseudo_mat()
#' @description generate pseudobulk expression for seurat object based on N neighbors RNA expression
#' @param obj seurat object, need to run FindNeighbors previously.
#' @param expMatrix sparse Matrix of class "dgCMatrix" for expression
gen_pseudo_mat <- function(obj,randomseeds,Neighbors=100,assay="RNA",expMatrix=obj[[assay]]@data,rowValue="rowSums"){
  seedNeighbors <- lapply(randomseeds,function(i,obj){
    return(TopNeighbors(obj@neighbors[[paste0(assay,".nn")]],i,n=Neighbors))
  },obj)
  if(rowValue!="rowSums"){
    pseudoExp <- expMatrix[,randomseeds,drop=F]
    pseudoExp <- as.matrix(pseudoExp)
    pseudoExp <- do.call(cbind,lapply(1:length(randomseeds), function(i,randomseeds,seedNeighbors,expMatrix){
      subcell=c(randomseeds[i],seedNeighbors[[i]])
      subexp <- expMatrix[,subcell,drop=F]
      subexp <- as.matrix(subexp)
      return(rowMedians(subexp,na.rm=TRUE))
    },randomseeds,seedNeighbors,expMatrix))
  }else{
    pseudoExp <- do.call(cbind,lapply(1:length(randomseeds), function(i,randomseeds,seedNeighbors,expMatrix){
      subcell=c(randomseeds[i],seedNeighbors[[i]])
      subexp <- expMatrix[,subcell,drop=F]
      return(rowSums(subexp,na.rm=TRUE))
    },randomseeds,seedNeighbors,expMatrix))
  }
 
  rownames(pseudoExp) <- rownames(expMatrix)
  colnames(pseudoExp)<- randomseeds
  return(pseudoExp)
}

#' @title pseudo_bulk_v2()
#'
#' @description segment by computing the changepoint for an ordered vector
#'  (e.g., peaks or bins across chromosome p arm and q arm location)
#'  
#' @param mat_obj matrix with features as rownames, cell Id as colnames
#' @param group_ano a data.frame of cell_group needs to be provided, the first column is cell ID, the second column is cell groups
#' @param type sum or median expression of cells in a group c("mean","sum","median","density")
pseudo_bulk_v2 <- function(mat_obj,group_ano,method="sum",adjust_zero=TRUE){
  if(!is.matrix(mat_obj)){ mat_obj <- as.matrix(mat_obj)}
  myrowname <- rownames(mat_obj)
  if(is.factor(group_ano[,2])){
    unique_groups <- levels(group_ano[,2])
  }else{
    unique_groups <- unique(as.character(group_ano[,2]))
  }
  groups=table(group_ano[,2])[unique_groups]
  mat_groups <- c()
  for(i in 1:length(unique_groups)){
    grp <- unique_groups[i]
    group_cells <- group_ano[which(as.character(group_ano[,2]) %in% grp),1]
    group_cells_idx <- which(colnames(mat_obj) %in% group_cells)
    if(!any(method%in%c("mean","density","sum","median"))){
      message('method must be one of c("mean","density","sum","median"). Default is sum.')
      stop()
    }
    if(method=="mean"){
      group_feature_means <- .get_ref_mean(mat_obj,group_cells_idx,zero.rm=adjust_zero)
      mat_groups <- cbind(mat_groups,group_feature_means$mean)
    }else if(method=="density"){
      group_feature_means <- .get_ref_density(mat_obj,group_cells_idx,zero.rm=adjust_zero)
      mat_groups <- cbind(mat_groups,group_feature_means$densi_peak)
    }else{
      group_feature_means <- .get_ref_median(mat_obj,group_cells_idx,zero.rm=adjust_zero)
      if(method=="sum"){
        mat_groups <- cbind(mat_groups,group_feature_means$sum)
      }else if(method=="median"){
        mat_groups <- cbind(mat_groups,group_feature_means$mean)
      }
    } 
  }
  mat_groups <- as.data.frame(mat_groups)
  colnames(mat_groups) <- unique_groups
  rownames(mat_groups) <- myrowname
  return(mat_groups)
}

#' @title .get_ref_median()
#'@keywords internal function
#'@param matrix input matrix, group in columns (cells) 
#'@param ref_indice numeric indices for group
.get_ref_median <- function(matrix,ref_indice,zero.rm=FALSE){
  library(matrixStats)
  RN <- rownames(matrix)
  matrix_data <- matrix[,ref_indice,drop=FALSE]
  if(zero.rm){
    peak_ref_grp_median <- apply(matrix_data, 1, function(x){
      if(all(x==0)){
        adjusted_mean=0
      }else{
        non_zero_mean <- median(x[x != 0],na.rm=TRUE)
        zero_proportion <- sum(x == 0,na.rm=TRUE) / length(x)
        adjusted_mean <- non_zero_mean * (1 - zero_proportion)
      }
      return(adjusted_mean)
      })
  }else{
    peak_ref_grp_median <- rowMedians(matrix_data,na.rm=TRUE)
  }
  
  peak_ref_grp_sd <- rowSds(matrix_data,na.rm=TRUE)
  peak_ref_grp_sum <- rowSums(matrix_data,na.rm=TRUE)
  peak_ref_grp_means_sd <- data.frame(mean=peak_ref_grp_median,SD=peak_ref_grp_sd,sum=peak_ref_grp_sum)
  rownames(peak_ref_grp_means_sd) <- RN
  return(peak_ref_grp_means_sd)
}

#'@keywords internal function
#'@param matrix input matrix, group in columns (cells) 
#'@param ref_indice numeric indices for group
.get_ref_mean <- function(matrix,ref_indice,zero.rm=FALSE){
  library(matrixStats)
  RN <- rownames(matrix)
  matrix_data <- matrix[,ref_indice,drop=FALSE]
  if(zero.rm){
    peak_ref_grp_mean <- apply(matrix_data, 1, function(x){
      #the condition mean of Zero-Inflated distribution
      if(all(x==0)){
        adjusted_mean=0
      }else{
        non_zero_mean <- mean(x[x != 0],na.rm=TRUE)
        zero_proportion <- sum(x == 0,na.rm=TRUE) / length(x)
        adjusted_mean <- non_zero_mean *(1 - zero_proportion)
      }
      return(adjusted_mean)
    } )
  }else{
    peak_ref_grp_mean <- rowMeans(matrix_data,na.rm=TRUE)
  }
  
  peak_ref_grp_sd <- rowSds(matrix_data,na.rm=TRUE)
  peak_ref_grp_sum <- rowSums(matrix_data,na.rm=TRUE)
  peak_ref_grp_means_sd <- data.frame(mean=peak_ref_grp_mean,SD=peak_ref_grp_sd,sum=peak_ref_grp_sum)
  rownames(peak_ref_grp_means_sd) <- RN
  return(peak_ref_grp_means_sd)
}

#' @title .get_ref_density()
#'@keywords internal function
#'@param matrix input matrix, group in columns (cells) 
#'@param ref_indice numeric indices for group
.get_ref_density <- function(matrix,ref_indice,zero.rm=FALSE){
  library(matrixStats)
  RN <- rownames(matrix)
  matrix_data <- matrix[,ref_indice,drop=FALSE]
  if(zero.rm){
    peak_ref_grp <- apply(matrix_data, 1, function(x){
      x <- x[x != 0]
      if(length(x)>1){
        density <- density(x)
        max_index <- which.max(density$y)
        peak_x <- density$x[max_index]
      }else{
        peak_x <- x
      }
    return(peak_x)
    })
  }else{
    peak_ref_grp <- apply(matrix_data, 1, function(x){
      if(length(x)>1){
        density <- density(x)
        max_index <- which.max(density$y)
        peak_x <- density$x[max_index]
      }else{
        peak_x <- x
      }
      return(peak_x)
    })
  }
  
  peak_ref_grp_sd <- rowSds(matrix_data,na.rm=TRUE)
  peak_ref_grp_sum <- rowSums(matrix_data,na.rm=TRUE)
  peak_ref_grp_means_sd <- data.frame(densi_peak=peak_ref_grp,SD=peak_ref_grp_sd,sum=peak_ref_grp_sum)
  rownames(peak_ref_grp_means_sd) <- RN
  return(peak_ref_grp_means_sd)
}




