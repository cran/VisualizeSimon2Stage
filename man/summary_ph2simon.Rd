% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary_ph2simon}
\alias{summary_ph2simon}
\alias{summary.ph2simon4}
\alias{summary.ph2simon}
\title{Summarize a Simon's Two-Stage Design}
\usage{
\method{summary}{ph2simon4}(object, ...)

\method{summary}{ph2simon}(object, ...)
}
\arguments{
\item{object}{a \link[clinfun]{ph2simon} or \linkS4class{ph2simon4} object}

\item{...}{parameters of function \code{\link[=ph2simon4]{ph2simon4()}}, most importantly \code{type}}
}
\value{
Function \code{\link[=summary.ph2simon]{summary.ph2simon()}} returns a \link[base]{list} with three (3) elements
\describe{
\item{\code{'design'}}{\link[base]{integer} \link[base]{matrix}}
\item{\code{'EN'}}{\link[base]{double} \link[base]{matrix}}
\item{\code{'p'}}{\link[base]{double} \link[base]{matrix}}
}
}
\description{
Summarize a Simon's two-stage design
}
\keyword{internal}
