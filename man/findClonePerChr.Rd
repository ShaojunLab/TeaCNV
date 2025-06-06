% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funs_find_subclones.R
\name{findClonePerChr}
\alias{findClonePerChr}
\title{findClonePerChr()}
\usage{
findClonePerChr(
  mat,
  ref.value,
  outdir = ".",
  ChrExclude = c("chrX", "chrY", "chrM"),
  clust_resolution = 1
)
}
\arguments{
\item{mat}{peak-cell count matrix of whole genome. peaks as rownames with "chrx-xxx-xxx", "chrx:xxx-xxx" or "chrx_xxx_xxx" format}

\item{ref.value}{}

\item{ChrExclude}{the list of chromosome should be excluded from analysis. Default = c('chrX', 'chrY', 'chrM')}
}
\description{
sub-clustering cells per chromosome count matrix if the method is "perChr"
}
