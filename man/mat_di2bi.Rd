% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrixOps.R
\name{mat_di2bi}
\alias{mat_di2bi}
\title{Distance to binary matrix}
\usage{
mat_di2bi(distmat, emRad = NA, convMat = FALSE)
}
\arguments{
\item{distmat}{Distance matrix}

\item{emRad}{The radius or threshold value}

\item{convMat}{Should the matrix be converted from a \code{distmat} object of class \code{\link[Matrix:Matrix]{Matrix::Matrix()}} to \code{\link[base:matrix]{base::matrix()}} (or vice versa)}
}
\value{
A (sparse) matrix with only 0s and 1s
}
\description{
Distance matrix to binary matrix based on threshold value
}
\seealso{
Other Distance matrix operations (recurrence plot): 
\code{\link{bandReplace}()},
\code{\link{createCorridor}()},
\code{\link{mat_di2ch}()},
\code{\link{mat_di2we}()},
\code{\link{mat_hamming}()},
\code{\link{rp}()},
\code{\link{rp_lineDist}()},
\code{\link{rp_nzdiags}()},
\code{\link{rp_plot}()},
\code{\link{rp_size}()}

Other Distance matrix operations (recurrence network): 
\code{\link{mat_di2ch}()},
\code{\link{mat_di2we}()},
\code{\link{rn}()},
\code{\link{rn_phases}()},
\code{\link{rn_plot}()},
\code{\link{rn_recSpec}()}
}
\concept{Distance matrix operations (recurrence network)}
\concept{Distance matrix operations (recurrence plot)}
