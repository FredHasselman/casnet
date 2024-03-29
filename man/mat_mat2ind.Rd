% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrixOps.R
\name{mat_mat2ind}
\alias{mat_mat2ind}
\title{Matrix to indexed data frame}
\usage{
mat_mat2ind(mat)
}
\arguments{
\item{mat}{A matrix}
}
\value{
A data.frame with two index columns named \code{"Var1"} and \code{"Var2"}
}
\description{
Mimics the default behaviour of \code{reshape2::melt()}
}
\examples{
mat_mat2ind(as.matrix(1:100,ncol=10))

}
