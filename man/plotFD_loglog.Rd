% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFUN.R
\name{plotFD_loglog}
\alias{plotFD_loglog}
\title{Plot output from fluctuation analyses based on log-log regression}
\usage{
plotFD_loglog(
  fd.OUT,
  title = "",
  subtitle = "",
  xlabel = "Bin size",
  ylabel = "Fluctuation",
  logBase = NA,
  doPlot = TRUE,
  returnPlot = FALSE
)
}
\arguments{
\item{fd.OUT}{Output from one of the \code{fd_} functions that use log-log regression to get scaling exponents.}

\item{title}{Plot title}

\item{subtitle}{Plot subtitle}

\item{xlabel}{x label}

\item{ylabel}{y label}

\item{logBase}{base of the log used}

\item{doPlot}{Display the plot (A plot object is always returned invisibly)}

\item{returnPlot}{return a \link{ggplot2} object}
}
\value{
A ggplot object
}
\description{
Plot output from fluctuation analyses based on log-log regression
}
