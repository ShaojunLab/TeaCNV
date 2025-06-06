% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_segmentation_comb.R
\name{segCNV_refine}
\alias{segCNV_refine}
\title{segCNV_refine()}
\usage{
segCNV_refine(
  clonalRes,
  CNest.ref,
  cytoBand,
  seg_method = "PELT",
  seg.count.lim = 80,
  penalty = 1,
  segValue_method = "median",
  FiltSeg = TRUE,
  SegLen_min = 1e+06,
  SegSize_min = 5,
  p.adj_cutoff = 0.05,
  ratio_diff_cutoff = NULL
)
}
\description{
segCNV_refine()
}
