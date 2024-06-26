% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts.R
\name{ts_slice}
\alias{ts_slice}
\title{Slice a Matrix}
\usage{
ts_slice(y, epochSz = 4, overlap = NA, removeUnequal = FALSE)
}
\arguments{
\item{y}{A matrix with timeseries as columns}

\item{epochSz}{Epoch size}

\item{removeUnequal}{Do not return bins whose length is not equal to \code{epochSz} (default = \code{FALSE})}
}
\value{
A list with epochs
}
\description{
Slices rows of a matrix into a list of matrices representing epochs of length \code{epochSz}.
}
\seealso{
Other Time series operations: 
\code{\link{ts_center}()},
\code{\link{ts_changeindex}()},
\code{\link{ts_checkfix}()},
\code{\link{ts_detrend}()},
\code{\link{ts_diff}()},
\code{\link{ts_discrete}()},
\code{\link{ts_duration}()},
\code{\link{ts_embed}()},
\code{\link{ts_integrate}()},
\code{\link{ts_levels}()},
\code{\link{ts_peaks}()},
\code{\link{ts_permtest_block}()},
\code{\link{ts_permtest_transmat}()},
\code{\link{ts_rasterize}()},
\code{\link{ts_sd}()},
\code{\link{ts_slopes}()},
\code{\link{ts_standardise}()},
\code{\link{ts_sumorder}()},
\code{\link{ts_symbolic}()},
\code{\link{ts_trimfill}()},
\code{\link{ts_windower}()}
}
\author{
Fred Hasselman
}
\concept{Time series operations}
