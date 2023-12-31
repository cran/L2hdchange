% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nbd_functions.R
\name{get_breaks.nbd}
\alias{get_breaks.nbd}
\title{Obtain the time-stamps and spatial locations with breaks}
\usage{
\method{get_breaks}{nbd}(estobj)
}
\arguments{
\item{estobj}{An S3 object of class 'no_nbd' or 'nbd' generated
by \code{\link[=est_hdchange]{est_hdchange()}}.}
}
\value{
A list containing the total number of breaks \eqn{\widehat{R}} and the
spatial-temporal location of the break \eqn{(\hat{\tau}_{r},\hat{s}_{r})}. See Algorithm
2 of Li et al. (2023).
}
\description{
Obtain the time-stamps and spatial locations with breaks
}
\references{
Li, J., Chen, L., Wang, W. and Wu, W.B., 2022. \eqn{\ell^2} Inference for Change Points in High-Dimensional Time Series via a Two-Way MOSUM.
\emph{arXiv preprint arXiv:2208.13074}.
}
