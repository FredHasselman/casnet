% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ts.R
\name{ts_checkfix}
\alias{ts_checkfix}
\title{Check and/or Fix a vector}
\usage{
ts_checkfix(
  y,
  checkNumericVector = TRUE,
  checkWholeNumbers = FALSE,
  checkTimeVector = FALSE,
  checkPow2 = FALSE,
  checkScale = FALSE,
  checkSummationOrder = FALSE,
  checkNonStationarity = FALSE,
  checkNonHomogeneity = FALSE,
  fixNumericVector = FALSE,
  fixWholeNumbers = FALSE,
  fixTimeVector = FALSE,
  fixPow2 = FALSE,
  fixNA = TRUE,
  fixScale = FALSE,
  fixSummationOrder = FALSE,
  fixNonStationarity = FALSE,
  fixNonHomogeneity = FALSE
)
}
\arguments{
\item{y}{A time series object or numeric vector}

\item{checkNumericVector}{is 1D numeric vector?}

\item{checkWholeNumbers}{contains only wholenumbers?}

\item{checkTimeVector}{has time vector?}

\item{checkPow2}{length is power of 2?}

\item{checkScale}{checkScale}

\item{checkSummationOrder}{checkSummationOrder}

\item{checkNonStationarity}{checkNonStationarity}

\item{checkNonHomogeneity}{checkNonHomogeneity}

\item{fixNumericVector}{return a 1D numeric vector (WARNING: Data frames and Matrices with NCOL > 1 wil be converted to long form)}

\item{fixWholeNumbers}{fixWholeNumber}

\item{fixTimeVector}{fixTimeVector}

\item{fixPow2}{foxPow2}

\item{fixNA}{fixNA}

\item{fixScale}{fixScale}

\item{fixSummationOrder}{fixSummationOrder}

\item{fixNonStationarity}{fixNonStationarity}

\item{fixNonHomogeneity}{fixNonHomogeneity}
}
\value{
A 'check' report and/or a 'fixed' vector y.
}
\description{
Check and/or Fix a vector
}
\seealso{
Other Time series operations: 
\code{\link{ts_center}()},
\code{\link{ts_changeindex}()},
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
\concept{Time series operations}
