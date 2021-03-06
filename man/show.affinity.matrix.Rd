% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Descriptive.R
\name{show.affinity.matrix}
\alias{show.affinity.matrix}
\title{Print affinity matrix}
\usage{
show.affinity.matrix(
  res,
  labels_x = paste0("Trait ", 1:Kx),
  labels_y = paste0("Trait ", 1:Ky),
  pr = 0.05
)
}
\arguments{
\item{res}{A list corresponding to the output of
\code{\link{estimate.affinity.matrix}},
\code{\link{estimate.affinity.matrix.lowrank}} or
\code{\link{estimate.affinity.matrix.unipartite}}.}

\item{labels_x}{A vector of strings indicating the names of men's matching
variables. Defaults to \code{"Trait k"} for every \code{k} matching
variable.}

\item{labels_y}{A vector of strings indicating the names of women's matching
variables. Defaults to \code{"Trait k"} for every \code{k} matching
variable.}

\item{pr}{A probability indicating the two-tailed significance level required
for an estimated parameter to be printed in boldface. Defaults to 0.05
and can be set to 0 to avoid printing any estimate in boldface.}
}
\value{
The function returns a long string in LaTeX style that can be
processed in the standard LaTeX tabular environment in order to display
the estimates of the affinity matrix \code{Aopt}.
}
\description{
This function prints the estimated affinity matrix in LaTeX style. Standard
errors are printed below the elements of the affinity matrix. Estimates that
are significant at the \code{pr} level are printed in boldface: this format
feature can be avoided by setting \code{pr} to 0.
}
