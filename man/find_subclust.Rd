% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs_find_subclones.R
\name{find_subclust}
\alias{find_subclust}
\title{find_subclust()}
\usage{
find_subclust(
  mt,
  add_metadata = NULL,
  resolution = 0.8,
  sep_by = c("-", "-"),
  assay = "binCount"
)
}
\description{
Find sub-cluster based on bin-level ratio using PCA method of seurat
}
