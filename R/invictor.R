#' Rose tinted infix
#'
#' @param x If (an element of) \code{x} is any of \code{Inf,-Inf,NA,NaN,NULL,length(x)==0}, it will return/replace the value of \code{y}; otherwise \code{x}.
#' @param y The value to return/replace for \code{x} in case of catastrophy \code{>00<}
#'
#' @export
#' @author Fred Hasselman
#' @description When your functions wear these rose tinted glasses, the world will appear to be a nicer, fluffier place.
#'
#' @seealso purrrr::%||%
#' @keywords internal
#' @examples
#'
#' Inf %00% NA
#'
#' numeric(0) %00% ''
#'
#' NA %00% 0
#'
#' NaN %00% NA
#'
#' NULL %00% NA
#'
`%00%` <- function(x, y) {
  if (length(x) == 0) {
    x <- y
  } else {
    for (i in seq_along(x)) {
      l0 <- isna <- isnan <- isinf <- isnll <- isTryError <- FALSE
      if (length(x[i]) == 0) {
        l0 <- TRUE
      } else {
        if (all(is.na(x[i]))) {isna <- TRUE}
        if (all(is.nan(x[i]))) {isnan <- TRUE}
        if (all(is.infinite(x[i]))) {isinf <- TRUE}
        if (all(is.null(x[i]))) {isnll <- TRUE}
        if (all(class(x[i]) %in% "try-error")) {isTryError <- TRUE}
      }
      if (any(l0, isna, isnan, isinf, isnll, isTryError)) {
        x[i] <- y
      }
    }
  }
  return(x)
}

#' Inside interval
#'
#' Decide if a value \code{x} falls inside an interval \code{j[1],j[2]} that can be open or closed on the left and/or the right. Either a logical vector equal to \code{x}, or the actual values are extracted, when the `.`-versions are used.
#'
#' @param x A vector
#' @param j A 2-element numeric vector indicating a range
#'
#' @note Package `DescTools` provides similar functions
#' @keywords internal
#'
#' @name insiders
#'
#' @return Logical vector of length \code{x}, or, values in the range \code{j}
#'
#' @examples
#'
#' # Closed interval
#' 0:5 %[]% c(1,5)  # logical vector
#' 0:5 %[.]% c(1,5) # extract values
#'
#' # Open interval
#' 0:5 %()% c(1,5)
#' 0:5 %(.)% c(1,5)
#'
#' # Closed interval left
#' 0:5 %[)% c(1,5)
#' 0:5 %[.)% c(1,5)
#'
#' # Closed interval right
#' 0:5 %(]% c(1,5)
#' 0:5 %(.]% c(1,5)
#'
#'
NULL
# >NULL


# Insiders ----

#'  In closed interval
#'
#' @rdname insiders
#' @export
#' @keywords internal
#'
`%[]%` <- function(x, j) {
  if(all(length(j) == 2, is.numeric(x), is.numeric(j))){
    rng <- sort(j)
    x >= rng[1] & x <= rng[2]
  }
}

#'  In open interval
#'
#' @rdname insiders
#' @export
#' @keywords internal
#'
#'
`%()%` <- function(x, j) {
  if (all(length(j) == 2, is.numeric(x), is.numeric(j))) {
    rng <- sort(j)
    x > rng[1] & x < rng[2]
  }
}



# Outsiders -----

#' Outside interval
#'
#' Decide if a value \code{x} falls outside an interval \code{j[1],j[2]} that can be open or closed on the left and/or the right. Either a logical vector equal to \code{x}, or the actual values are extracted,
#'
#' @param x A vector
#' @param j A range
#'
#' @note Package `DescTools` provides similar functions
#'
#' @name outsiders
#'
#' @return logical vector of length x, or, values of x outside the range j
#' @keywords internal
#'
#' @examples
#'
#' # Closed interval
#' 5%][%c(1,5)
#' 5%].[%c(1,5)
#'
#' # Open interval
#' 5%)(%c(1,5)
#' 5%).(%c(1,5)
#'
#' # Half-losed interval left
#' 5%](%c(1,5)
#' 5%].(%c(1,5)
#'
#' # Half-losed interval right
#' 5%)[%c(1,5)
#' 5%).[%c(1,5)
#'
#'
NULL
# >NULL


#'  Not in closed interval
#'
#' @rdname outsiders
#' @export
#' @keywords internal
#'
`%][%` <- function(x, j) {
  return(!x%()%j)
}

#'  Not in open interval
#'
#' @rdname outsiders
#' @export
#' @keywords internal
#'
`%)(%` <- function(x, j) {
  return(!x%[]%j)
}
