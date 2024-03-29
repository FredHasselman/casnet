% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rn.R
\name{rn_recSpec}
\alias{rn_recSpec}
\title{Recurrence Time Spectrum}
\usage{
rn_recSpec(
  RN,
  fitRange = NULL,
  fs = 1,
  doPlot = TRUE,
  returnPlot = FALSE,
  returnPLAW = FALSE,
  returnInfo = FALSE,
  silent = TRUE,
  noTitle = FALSE,
  tsName = "y"
)
}
\arguments{
\item{RN}{A thresholded recurrence matrix generated by function \code{rn()}}

\item{fitRange}{If \code{NULL} the entire range will be used for log-log slope. If a 2-element vector of integers, this will represent the range of recurrence times to use for fitting the log=log slope (e.g. \code{c(1,50)} would fit the first 50 recurrence times).}

\item{fs}{Sample rate (default = \code{1})}

\item{doPlot}{Should a plot of the recurrence time spectrum be produced?}

\item{returnPlot}{Return ggplot2 object (default = \code{FALSE})}

\item{returnPLAW}{Return the power law data (default = \code{FALSE})}

\item{returnInfo}{Return all the data used in SDA (default = \code{FALSE})}

\item{silent}{Silent-ish mode}

\item{noTitle}{Do not generate a title (only the subtitle)}

\item{tsName}{Name of y added as a subtitle to the plot}
}
\value{
A vector of frequencies of recurrence times and a plot (if requested)
}
\description{
Get the recurrence time distribution from a recurrence network.
}
\seealso{
Other Distance matrix operations (recurrence network): 
\code{\link{mat_di2bi}()},
\code{\link{mat_di2ch}()},
\code{\link{mat_di2we}()},
\code{\link{rn}()},
\code{\link{rn_phases}()},
\code{\link{rn_plot}()}
}
\concept{Distance matrix operations (recurrence network)}
