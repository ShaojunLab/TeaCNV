% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs_filtMT.R
\name{FiltPeak}
\alias{FiltPeak}
\title{FiltPeak()}
\usage{
FiltPeak(
  mtx,
  cell_anno = NULL,
  ref_group_names = NULL,
  cellprop = 0.05,
  outheatmap = FALSE,
  blacklist_file = NULL,
  ChrRemove = c("chrX", "chrY", "chrM"),
  prop_reset = TRUE
)
}
\arguments{
\item{cell_anno}{annotation of cells with cell groups information on the first column.}

\item{blacklist_file}{.bed file of blacklist}
}
\description{
FiltPeak()
}
