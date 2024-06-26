% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFUN.R
\name{plotRN_phaseProjection}
\alias{plotRN_phaseProjection}
\title{Plot Phase Space Projection}
\usage{
plotRN_phaseProjection(
  RNdist,
  phaseOutput,
  epochColours = NULL,
  showEpochLegend = TRUE,
  epochLabel = "Phase",
  excludeOther = TRUE,
  excludeNorec = TRUE
)
}
\arguments{
\item{RNdist}{A distance matrix (unthresholded) created with \link{rn} or \link{rp}}

\item{phaseOutput}{Output from function \link{rn_phases}}

\item{excludeOther}{Exclude the default Phase "Other"}

\item{excludeNorec}{Exclude the default Phase "No recurrence"}

\item{PhaseOut}{Output from function \link{rn_phases} based on the thresholded matrix in \code{RNdist}}
}
\value{
Dataframe with coordinates
}
\description{
2D umap projection of multidimensional Phase Space. Categories "No recurrence" and "Other" are removed.
}
