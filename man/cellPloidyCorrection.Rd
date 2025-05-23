% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ploidy.correction.v2.R
\name{cellPloidyCorrection}
\alias{cellPloidyCorrection}
\title{cellPloidyCorrection()}
\usage{
cellPloidyCorrection(celldata, clonalest, clonalest.Ref = NULL)
}
\arguments{
\item{celldata}{cell level ratio}

\item{clonalest.Ref}{reference clone if the user specifies}

\item{clonalest:}{clonal ploidy estimation}
}
\value{
corrected ratio matrix
}
\description{
centering cell ratio both segment and cell level based on the clonal expectation
}
