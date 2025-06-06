% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genome.subclonality.function.R
\name{refine_clones}
\alias{refine_clones}
\title{refine_clones()}
\usage{
refine_clones(ratiodata, seg_compare_ls, best_clone_ratiodata, CNest.ref)
}
\arguments{
\item{ratiodata}{A list of data.frames containing the bin-level ratio estimates for each clone.}

\item{seg_compare_ls}{A list of data.frames containing the segment ratios for each clone compared to the best clone.}

\item{best_clone_ratiodata}{The name of the best clone for the segment ratio comparison in ratiodata.}
}
\description{
This function refines the clonal CNA estimates based on the paired-comparison of segment ratios.
}
