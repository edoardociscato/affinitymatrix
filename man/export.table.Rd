% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Descriptive.R
\name{export.table}
\alias{export.table}
\title{Export an affinitymatrix table}
\usage{
export.table(tabular, name = "table", path = getwd())
}
\arguments{
\item{tabular}{A long string corresponding to the output of
\code{\link{show.affinity.matrix}}, \code{\link{show.diagonal}} or
\code{\link{show.test}}, or one of the two elements of
\code{\link{show.saliency}} (\code{U.table} or \code{V.table}).}

\item{name}{A string indicating the name of the txt file. Defaults to
\code{"affinity_matrix"}.}

\item{path}{A string indicating the path where to save the txt file. Defaults
to current path.}
}
\value{
The function stores a long string in LaTeX style that can be
processed in the standard LaTeX tabular environment in a txt file in
located in \code{path}.
}
\description{
The function stores a LaTeX style table in a txt file.
}
