% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts.R
\name{ts_center}
\alias{ts_center}
\title{Center a vector}
\usage{
ts_center(numvec, na.rm = TRUE, type = c("mean", "median")[1])
}
\arguments{
\item{numvec}{A numeric vector}

\item{na.rm}{Set the \code{na.rm} field}

\item{type}{Center on the \code{"mean"} (default) or the \code{"median"} of the vector.}
}
\value{
A mean or median centered vector
}
\description{
Center a vector
}
\seealso{
Other Time series operations: 
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
\code{\link{ts_slice}()},
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
