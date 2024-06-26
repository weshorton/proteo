% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classifyProteo.R
\name{classifyProteo}
\alias{classifyProteo}
\title{Classify Proteomics Data}
\usage{
classifyProteo(
  data_dt,
  lfcCol_v = "logFC",
  pvalCol_v = "FDR",
  lfc_v = 0.5,
  pval_v = c(0.01, 0.05),
  newName_v = "diffExp"
)
}
\arguments{
\item{data_dt}{proteomics data.table with differential expression results}

\item{lfcCol_v}{column in data_dt indicating the fold change results (logFC for TMT; Z-score for Silac)}

\item{pvalCol_v}{column in data_dt indicating the significance value (FDR for TMT; BH_adjusted for Silac)}

\item{lfc_v}{value for cut-off for fold change. Will be included in plot as vertical dotted line}

\item{pval_v}{vector of length 2! first value indicates the threshold for 'high' significance, second for 'low'}

\item{newName_v}{name of new column}
}
\value{
data_dt with one new column
}
\description{
Group data_dt into 3 levels each for up/down expression:
low sig: 0.1 > p >= 0.05
medium sig: 0.05 > p >= 0.01
high sig: 0.01 > p
}
\details{
Add a new column in data_dt summarizing the direction and level of each entry.
}
