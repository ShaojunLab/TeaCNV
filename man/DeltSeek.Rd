% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ploidy.correction.v2.R
\name{DeltSeek}
\alias{DeltSeek}
\title{DeltSeek()}
\usage{
DeltSeek(df_seg_C1, initialCN1, by_term = "SegMean")
}
\arguments{
\item{df_seg_C1:}{segmentation data}

\item{initialCN1:}{the estimation of ratio baseline}
}
\value{
all distance value
}
\description{
finding the possible distance between integerCN states used in PloidyCorrect
}
