#' A Toolbox for Studying Complex Adaptive Systems and NETworks
#'
#' @author Fred Hasselman
#'
#' @name casnet
#' @docType package
#' @import DescTools
#' @importFrom magrittr %>%
#' @importFrom plyr .
#' @import ggplot2
#' @importFrom igraph %--%
#'
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1"){utils::globalVariables(c("."))}
