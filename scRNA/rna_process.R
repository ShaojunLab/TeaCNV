#' @Description Convert gene symbols in RNA matrix to genomic coordinates
#' @param mat matrix with gene name as rowname
#' @param reference genes.gtf.gz file with path
rnaMTprocess <- function(mat,reference,chr_remove=c("chrM","chrX","chrY")){	
	annotations_all <- rtracklayer::import(reference)
	annotations <- annotations_all[annotations_all$type=="gene",]
	annotations <- annotations[!grepl(paste0("^Mt|^MT|^X|^Y|^GL|^KI|",paste0(chr_remove,collapse ="|")),as.character(seqnames(annotations))), ]
	annotations <- annotations[annotations$gene_name %in% rownames(mat),]
	mat <- mat[rownames(mat)%in%annotations$gene_name, ]
	return(mat)
}