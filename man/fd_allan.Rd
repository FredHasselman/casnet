% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fd.R
\name{fd_allan}
\alias{fd_allan}
\title{Allan Variance Analysis}
\usage{
fd_allan(
  y,
  fs = stats::tsp(stats::hasTsp(y))[3],
  useSD = FALSE,
  doPlot = FALSE,
  returnPlot = FALSE,
  returnPLAW = FALSE,
  returnInfo = FALSE,
  silent = FALSE,
  noTitle = FALSE,
  tsName = "y"
)
}
\arguments{
\item{y}{A numeric vector or time series object}

\item{fs}{Sample frequency in Hz}

\item{useSD}{Use the standarddeviation instead of variance?}

\item{doPlot}{Return the log-log scale versus fluctuation plot with linear fit (default = \code{TRUE}).}

\item{returnPlot}{Return ggplot2 object (default = \code{FALSE})}

\item{returnPLAW}{Return the power law data (default = \code{FALSE})}

\item{returnInfo}{Return all the data used in DFA (default = \code{FALSE})}

\item{silent}{Silent-ish mode}

\item{noTitle}{Do not generate a title (only the subtitle)}

\item{tsName}{Name of y added as a subtitle to the plot}
}
\value{
A dataframe with the Allan Factor (variance), Alan standard deviation and error due to bin size
}
\description{
Allan Variance Analysis
}
\seealso{
Other Fluctuation Analyses: 
\code{\link{fd_RR}()},
\code{\link{fd_dfa}()},
\code{\link{fd_mfdfa}()},
\code{\link{fd_psd}()},
\code{\link{fd_sda}()},
\code{\link{fd_sev}()}
}
\concept{Fluctuation Analyses}
