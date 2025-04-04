% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mydataProcess.R
\name{.smooth_by_loess}
\alias{.smooth_by_loess}
\title{.smooth_by_loess()}
\usage{
.smooth_by_loess(Mat, span = 0.33)
}
\arguments{
\item{Mat}{counts matrix}
}
\description{
loess smoothing
The size of the neighborhood can be controlled using the span argument,
which ranges between 0 to 1. It controls the degree of smoothing.
So, the greater the value of span, more smooth is the fitted curve. (default=1/3)
}
