% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/obtain_bicl.r
\name{check_biclusters}
\alias{check_biclusters}
\title{Check biclusters}
\usage{
check_biclusters(data, output_f, num_repeats)
}
\arguments{
\item{data}{list of data matrices}

\item{output_f}{list of F matrices}

\item{num_repeats}{number of repeats in removing spurious biclusters}
}
\value{
list containing the JSD scores,
average and max thresholds
}
\description{
Determine scores and thresholds of the factorisation
to be used to remove biclusters
}
