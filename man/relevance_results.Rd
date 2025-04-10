% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stability_analysis.r
\name{relevance_results}
\alias{relevance_results}
\title{calculate the relevance of biclusters}
\usage{
relevance_results(row_c, col_c, true_r, true_c)
}
\arguments{
\item{row_c}{row clusters}

\item{col_c}{column clusters}

\item{true_r}{true row clusters}

\item{true_c}{true column clusters}
}
\value{
vector of relevance values for each bicluster
}
\description{
calculate the relevannce of biclusters obtained from
sub-sampled data, using the initial biclustering
as the ground truth
}
