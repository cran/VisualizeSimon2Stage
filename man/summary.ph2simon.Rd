% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ph2simon_etc.R
\name{summary.ph2simon}
\alias{summary.ph2simon}
\title{Summarize a Simon's Two-Stage Design}
\usage{
\method{summary}{ph2simon}(object, ...)
}
\arguments{
\item{object}{\link[clinfun]{ph2simon} object}

\item{...}{potential parameters, currently not in use}
}
\value{
Function \link{summary.ph2simon} returns a \link[base]{list} with three (3) elements
\describe{
\item{\code{'design'}}{\link[base]{integer} \link[base]{matrix}}
\item{\code{'EN'}}{\link[base]{double} \link[base]{matrix}}
\item{\code{'p'}}{\link[base]{double} \link[base]{matrix}}
}
}
\description{
Summarize a Simon's two-stage design
}
\examples{
library(clinfun)
(x = ph2simon(pu = .2, pa = .4, ep1 = .05, ep2 = .1)) 
summary(x)

}
