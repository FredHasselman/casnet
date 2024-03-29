% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fd.R
\name{inf_SampEn}
\alias{inf_SampEn}
\title{Sample Entropy}
\usage{
inf_SampEn(
  y,
  m = 2,
  r = 0.2,
  D = NA,
  fs = NULL,
  standardise = c("none", "mean.sd", "median.mad")[1],
  transformBefore = TRUE,
  removeTrend = c("no", "poly", "adaptive", "bridge")[1],
  polyOrder = 1,
  relativeEntropy = FALSE,
  returnInfo = FALSE,
  silent = FALSE
)
}
\arguments{
\item{y}{A numeric vector or time series object.}

\item{m}{The size of the window in which tho evaluate whether a pattern repeats (default = \code{2})}

\item{r}{A factor that will determine the threshold for similarity of values, calculated as r x D (default = \code{0.2})}

\item{D}{Commonly the standard deviation of the time series, the similarity threshold will be calculated as r x D. Note that if the series is detrended and/or standardised and \code{D = NA} the standard deviation will be calculated after the transformations (default = \code{NA})}

\item{fs}{Sample rate}

\item{standardise}{Standardise the series using \code{\link[=ts_standardise]{ts_standardise()}} with \code{adjustN = FALSE} (default = "mean.sd")}

\item{transformBefore}{Detrend/standardise before coarse graining. If set to \code{FALSE}, each coarsegrained series will be detrended/standardised separately (default = \code{TRUE})}

\item{removeTrend}{Method to use for global detrending (default = \code{"poly"})}

\item{polyOrder}{Order of global polynomial trend to remove if \code{removeTrend = "poly"}. If \code{removeTrend = "adaptive"} polynomials \code{1} to \code{polyOrder} will be evaluated and the best fitting curve (R squared) will be removed (default = \code{1})}

\item{relativeEntropy}{The relative entropy, SampEn / (-1 * log(1/length(y))) will be returned (default = \code{FALSE})}

\item{returnInfo}{Return all the data used in SDA (default = \code{FALSE})}

\item{silent}{Silent-ish mode (default = \code{FALSE})}
}
\value{
The sample entropy (SampEn) of the time series y.
}
\description{
Sample Entropy
}
\examples{


y <- rnorm(100)

# Similarity threshold is r x D = 0.2 * sd(y)
inf_SampEn(y)

# Similarity threshold is r = 0.2
inf_SampEn(y, D = 1)

}
\seealso{
info_MSE

Other Information based complexity measures: 
\code{\link{inf_MSE}()}
}
\concept{Information based complexity measures}
