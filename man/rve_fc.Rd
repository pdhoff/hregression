% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/hreg_bayes.R
\name{rve_fc}
\alias{rve_fc}
\title{Bayes Hierarchical Linear Regression Posterior Simulation}
\usage{
rve_fc(Y, X, B, ss0 = 1, nu0 = 1)
}
\arguments{
\item{Y}{An \code{r} by \code{n} matrix.}

\item{X}{An \code{r} by \code{p} by \code{n} array.}

\item{B}{An \code{r} by \code{p} matrix.}
}
\value{
A scalar simulated
from the posterior distribution.
}
\description{
This functions returns a simulation from the posterior distribution
of the shared error variance.
}
\author{
Peter Hoff
}

