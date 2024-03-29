% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ppc_fit.R
\name{df_ppc}
\alias{df_ppc}
\alias{df_ppc.fitTK}
\alias{ppc}
\alias{ppc.fitTK}
\title{PPC data.frame}
\usage{
df_ppc(fit, ...)

\method{df_ppc}{fitTK}(fit, ...)

ppc(fit, ...)

\method{ppc}{fitTK}(fit, ...)
}
\arguments{
\item{fit}{And object returned by fitTK}

\item{\dots}{Additional arguments}
}
\value{
A data frame with median and 95\\% credible interval

a plot of class \code{ggplot}
}
\description{
This is the generic \code{ppc} S3 method for plots of the predicted
values along with 95\\% credible intervals
versus the observed values for \code{fitTK} objects.
}
\details{
The black points show the observed number of survivors (pooled
replicates, on \eqn{X}-axis) against the corresponding predicted
number (\eqn{Y}-axis). Predictions come along with 95\\% prediction
intervals, which are depicted in green when they contain the
observed value and in red otherwise. Samples with equal observed
value are shifted on the \eqn{X}-axis. For that reason, the
bisecting line (y = x), is represented by steps when observed
values are low. That way we ensure green intervals do intersect the
bisecting line.
}
