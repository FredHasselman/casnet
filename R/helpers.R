# HELPERS ----

# Misc. ----

# Check package
checkPkg <- function(pkg){
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package '", pkg,"' is needed for this function to work. Please install it."), call. = FALSE)
  } else {
    requireNamespace(pkg)
  }
}

# Convert decimal point
c2p <- function(text,N=1){
  if(!is.character(text)){text<-as.character(text)}
  if(sum(grepl("[,]",text))>=N){text <- gsub(",",".",text)}
  return(text)
}

#' Get all combinations
#'
#' @param x A vector with values.
#' @param y An optional vector with values
#'
#' @return a data frame with pairs of values
#' @export
#'
#' @examples
#'
#' getPairs(x=1:6)
#'
#' getPairs(x=1:3, y=c("A","B","C"))
#'
getPairs <- function (x, y = NULL){
  if (is.null(y)) {
    return(data.frame(t(utils::combn(x, 2L)), stringsAsFactors = FALSE))
  } else {
    return(expand.grid(x, y, stringsAsFactors = FALSE))
  }
}


# Count missing values in x
nmissing <- function(x){
  sum(is.na(x))
}


#' Numeric factor to numeric vector
#'
#' Converts a factor with numeric levels to a numeric vector, using the values of the levels.
#'
#' @param x A factor based on numeric values.
#' @param sortUnique Should the unique character values be sorted? (default = `FALSE`)
#' @param keepNA Keep NA values (`TRUE`), or remove them (default = `FALSE`)
#'
#' @return A numeric vector with factor levels as names.
#' @export
#'
#' @examples
#'
#' f <- factor(round(runif(10,0,9)))
#' as.numeric_factor(f)
#'
#' # Add NAs
#' f <- factor(c(round(runif(9,0,9)),NA))
#' as.numeric_factor(f)
#' as.numeric_factor(f, keepNA = TRUE)
#'
#'
as.numeric_factor <- function(x, keepNA = FALSE, sortUnique = FALSE){
  idNA <- is.na(x)
  if(!is.factor(x)){stop("Not a factor, use: `as.numeric_character()`")}
  out <- try(as.numeric(levels(x))[x])
  if(sum(is.na(out))>sum(idNA)){
    stop("Factor levels are not numeric!")
  }
  if(!keepNA){
    out <- out[!is.na(out)]
  }
  names(out) <- as.character(out)
  return(out)
}



#' Character vector to named numeric vector
#'
#' Converts a character vector to a named numeric vector, with the character elements as names.
#'
#' @param x A character vector
#' @param sortUnique Should the unique character values be sorted? (default = `FALSE`)
#' @param keepNA Keep NA values (`TRUE`), or remove them (default = `FALSE`)
#'
#' @return A named numeric vector
#' @export
#'
#' @examples
#'
#' f <- letters
#' as.numeric_character(f)
#'
as.numeric_character <- function(x, sortUnique = FALSE, keepNA = FALSE){
  IDna <- is.na(x)
  xx <- as.character(x[!IDna])
  labels.char <- unique(as.character(xx))
  if(sortUnique){labels.char <- sort(labels.char)}
  labels.num <- seq_along(labels.char)
  names(x) <- x
  xx <- x
  plyr::laply(labels.num, function(num){
    xx[xx%in%labels.char[num]]<<-num}
  )
  xx<-as.numeric(xx)
  names(xx) <- x
  return(xx)
}


#' Discrete (factor or character) to numeric vector
#'
#' Converts a factor with numeric levels, or, a character vector with numeric values to a numeric vector using [as.numeric_factor], or, [as.numeric_character] respectively. If an unnamed numeric vector is passed, it will be returned as a named numeric vector, if this vector is continuous, it will be returned discretised (by calling [ts_discrete]), the labels will be rounded to `signif(x, digits = 4).
#'
#' @param x A factor with levels that are numeric, or, a character vector representing numbers.
#' @param keepNA Keep NA values (`TRUE`), or remove them (default = `FALSE`)
#' @param sortUnique Should the unique character/factor level values be sorted? (default = `FALSE`)
#'
#' @return A named numeric vector with original factor levels / character values / numeric values as names.
#' @export
#'
#' @examples
#'
#' # Continuous
#' i <- runif(10,0,9)
#' as.numeric_discrete(i)
#'
#' # Integer
#' as.numeric_discrete(round(i))
#'
#' # Factor with NAs
#' f <- factor(c(round(runif(9,0,9)),NA))
#' as.numeric_discrete(f)
#' as.numeric_discrete(f, keepNA = FALSE)
#'
#' # Character vector
#' c <- c("Thank","you", "for", "the flowers")
#' as.numeric_discrete(c)
#' as.numeric_discrete(c, sortUnique = TRUE)
#'
#' c <- c("Thank","you", "for", "the", "flowers")
#' as.numeric_discrete(c)
#' as.numeric_discrete(c, sortUnique = TRUE)
#'
#'
as.numeric_discrete <- function(x, keepNA = FALSE, sortUnique = FALSE){

  idNA <- is.na(x)
  x <- x[!idNA]

  if(is.date(x)){
    x <- as.character(x)
  }

  if(plyr::is.discrete(x)){
    if(is.factor(x)){
      if(suppressWarnings(all(is.na(as.numeric(levels(x)))))){
        x <- as.character(x)
      } else {
        y <- as.numeric_factor(x, keepNA = keepNA, sortUnique = sortUnique)
      }
    }
    if(is.character(x)){
      y <- as.numeric_character(x, keepNA = keepNA, sortUnique = sortUnique)
    }
  } else {
    if(is.numeric(x)){
      if(!all(is.wholenumber(x))){
        y <- ts_discrete(x)
        if(is.null(names(y))){
          names(y) <- paste(signif(x,4))
        }
      } else { # wholenumber
        y <- x
        if(is.null(names(y))){
          names(y) <- paste(x)
        }
      }
    } else { # numeric
      stop("Variable is not a factor, character vector, or, unnamed numeric vector.")
    }
  } # discrete

  if(keepNA){
    tmp <- rep(NA,length(idNA))
    tmp[!idNA] <- y
    names(tmp)[!idNA] <- names(y)
    names(tmp)[idNA] <- NA
    return(tmp)
  } else {
    return(y)
  }
}

#' It's a Date!
#'
#'  Check if a vector is of the Date format.
#'
#' @param x A vector
#'
#' @return TRUE if it's a Date.
#'
#' @export
#'
is.date <- function(x){
  inherits(x, 'Date')
}



# doChecks <- function(y,standardise=FALSE,center=FALSE,){
#   nextPow2 <- nextn(y,factors = 2)
#   prevPow2 <- nextn(y,factors = 2)-1
# }


#' Wrapper for filtfilt
#'
#' @param TS A time series
#' @param f A filter
#'
#' @return A filtered signal
#' @export
#' @keywords internal
fltrIT <- function(TS,f){
  return(signal::filtfilt(f=f,x=TS))
}


#' Elastic Scaler - A Flexible Rescale Function
#'
#' @description The 'elastic scaler'will rescale numeric vectors (1D, or columns in a matrix or data.frame) to a user defined minimum and maximum, either based on the extrema in the data, or, a minimum and maximum defined by the user.
#'
#' @param x   Input vector or data frame.
#' @param mn  Minimum value of original, defaults to `min(x, na.rm = TRUE)` if set to `NA`.
#' @param mx  Maximum value of original, defaults to `max(x, na.rm = TRUE)` if set to `NA`.
#' @param lo  Minimum value to rescale to, defaults to `0`.
#' @param hi  Maximum value to rescale to, defaults to `1`.
#' @param groupwise If `x` is a data frame with `2+` columns, `mn = NA` and/or `mx = NA` and `groupwise = TRUE`, scaling will occur
#' @param keepNA Keep `NA` values?
#' @param boundaryPrecision If set to `NA` the precision of the input will be the same as the precision of the output. This can cause problems when detecting values that lie just outside of, or, exactly on boundaries given by `lo` and `hi`, e.g. after saving the data using a default precision. Setting `boundaryPrecision` to an integer value will ensure that the boundaries of the new scale are given by `round(..., digits = boundaryPrecision)`. Alternatively one could just round all the output after rescaling to a desired precision (default = `NA`)
#' @param tol The tolerance for deciding wether a value is on the boundary `lo` or `hi` (default = `.Machine$double.eps^0.5)`)
#'
#' @details Three uses:
#' \enumerate{
#' \item elascer(x)             - Scale x to data range: min(x.out)==0;      max(x.out)==1
#' \item elascer(x,mn,mx)       - Scale x to arg. range: min(x.out)==mn==0;  max(x.out)==mx==1
#' \item elascer(x,mn,mx,lo,hi) - Scale x to arg. range: min(x.out)==mn==lo; max(x.out)==mx==hi
#' }
#'
#' @return scaled inout
#' @export
#'
#' @examples
#' # Works on numeric objects
#' somenumbers <- cbind(c(-5,100,sqrt(2)),c(exp(1),0,-pi))
#'
#' # Using the defaults:
#' # 1. mn and mx are derived globally (groupWise = FALSE)
#' # 2. values rescaled to hi and lo are integers, 0 and 1 respectively
#' elascer(somenumbers)
#'
#' # If the data contain values < mn they will return as < lo
#' elascer(somenumbers,mn=-100)
#' # If the data contain values > mx they will return > hi
#' elascer(somenumbers,mx=99)
#'
#' # Effect of setting groupWise
#' elascer(somenumbers,lo=-1,hi=1)
#' elascer(somenumbers,lo=-1,hi=1, groupwise = TRUE)
#'
#' elascer(somenumbers,mn=-10,mx=100,lo=-1,hi=4)
#' elascer(somenumbers,mn= NA,mx=100,lo=-1,hi=4, groupwise = TRUE)
#'
#' # Effect of setting boundaryPrecision
#' x <- rbind(1/3, 1/7)
#'
#' re1 <- elascer(x, lo = 0, hi = 1/13, boundaryPrecision = NA)
#' max(re1)==0.07692308 # FALSE
#' max(re1)==1/13       # TRUE
#'
#' re2 <- elascer(x, lo = 0, hi = 1/13, boundaryPrecision = 8)
#' max(re2)==0.07692308 # TRUE
#' max(re2)==1/13       # FALSE
#'
elascer <- function(x,
                    mn=NA,
                    mx=NA,
                    lo=0,
                    hi=1,
                    groupwise = FALSE,
                    keepNA = TRUE,
                    boundaryPrecision = NA,
                    tol= .Machine$double.eps^0.5){

  doGroupwise <- FALSE
  mnNA <- FALSE
  mxNA <- FALSE
  UNLIST      <- FALSE
  if(length(dim(x))<2){UNLIST <- TRUE}
  if(any(is.na(mn),is.na(mx))){
    if(groupwise&NCOL(x)>1){doGroupwise <- TRUE}
    if(is.na(mn)){
      mnNA <- TRUE
      mn <- min(x,na.rm=TRUE)%00%NA
    }
    if(is.na(mx)){
      mxNA <- TRUE
      mx <- max(x,na.rm=TRUE)%00%NA
    }
    if(all(is.na(mn),is.na(mx))){
      warning(paste("Only missing values in vector. Returning x."))
      return(x)
    }
  }
  x <- as.data.frame(x)
  isNA <- plyr::colwise(is.na)(x)
  u <- x
  for(i in 1:NCOL(x)){
    if(doGroupwise){
      if(mnNA){mn<-min(x[,i],na.rm=TRUE)%00%NA}
      if(mxNA){mx<-max(x[,i],na.rm=TRUE)%00%NA}
      if(all(is.na(mn),is.na(mx))){warning(paste("Only missing values in column: ",i))}
    }
    if(mn>mx){warning("Minimum (mn) >= maximum (mx).")}
    if(lo>hi){warning("Lowest scale value (lo) >= highest scale value (hi).")}
    if(mn==mx){
      u[,i]<-rep(hi,length(x[,i]))
    } else {
      u[,i]<-(((x[,i]-mn)*(hi-lo))/((mx-mn))+lo)
      if(!keepNA){
        id<-stats::complete.cases(u[,i])
        u[!id,i]<-mn
      }
    }
    if(!is.na(boundaryPrecision)){
      # Make sure the end points of the scale are the same as passed in the argument
      idLo <- u[,i]%[]%c(lo-tol,lo+tol)
      u[idLo,i] <- round(u[idLo,i], digits = boundaryPrecision)
      idHi <- u[,i]%[]%c(hi-tol,hi+tol)
      u[idHi,i] <- round(u[idHi,i], digits = boundaryPrecision)
    }
  }
  if(UNLIST){
    u <- as.numeric(u[,1])
  }
  return(u)
}


# Rmd2htmlWP <- function(infile, outfile, sup = T) {
#   require(markdown)
#   require(knitr)
#   mdOpt <- markdownHTMLOptions(default = T)
#   mdOpt <- mdOpt[mdOpt != "mathjax"]
#   mdExt <- markdownExtensions()
#   mdExt <- mdExt[mdExt != "latex_math"]
#   if (sup == T) {
#     mdExt <- mdExt[mdExt != "superscript"]
#   }
#   knit2html(input = infile, output = outfile, options = c(mdOpt), extensions = c(mdExt))
# }

#' Heaviside step function
#'
#' @param value value
#'
#' @return heaviside change
#' @export
#' @keywords internal
#'
Heaviside <- function(value){

  if (value>0){
    h=1
  }
  else {h=0}
  return(h)
}


#' @title Return Error
#'
#' @param expr The expression
#' @return A list object with 2 fields: 'return' and 'message'.
#' @export
#' @keywords internal
return_error <- function(expr){

  Error  <- NULL

  catch_it <- function(error){
    Error <<- error
    invokeRestart("muffleWarning")
  }

  list(value = withCallingHandlers(tryCatch(expr, Error = function(e) e), warning = catch_it), warning = Error)
}


#' Repeat Copies of a Matrix
#'
#' @param X A matrix
#' @param m Multiply `dim(X)[1]` m times
#' @param n Multiply `dim(X)[2]` n times
#'
#' @return A repeated matrix
#' @export
#'
repmat <- function(X,m,n){
  #  REPMAT R equivalent of repmat (matlab)
  #  FORMAT
  #  DESC
  #  description not available.

  if (is.matrix(X))
  {
    mx = dim(X)[1]
    nx = dim(X)[2]
    out <- matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=TRUE)
  }
  else if (is.vector(X)) {
    mx = 1
    nx = length(X)
    out <- matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=TRUE)
  }
  else if (length(X) == 1)
  {
    out <- matrix(X,m, n)
  }
  return (out)
}



#' HurvichDeo
#'
#' Adaptation from discontinued package `fractal` by William Constantine and Donald Percival.
#'
#' @param nr NFT
#' @param spec sdf
#' @param A A
#' @param delta delta
#'
#' @return Hurvich-Deo estimate
#'
#' @export
#' @keywords internal
#'
HurvichDeo <- function(nr, spec, A=0.3, delta=6/7){

  # define frequencies
  omega <- ((2 * pi)/nr) * (1:nr)

  # notation as in reference: construct X matrix
  L  <- round(A * nr^delta)
  c1 <- rep(1, L)
  c2 <- log(abs(2 * sin(omega[1:L]/2)))
  c3 <- omega[1:L]^2/2
  X  <- cbind(c1, c2, c3)

  # construct least squares product matrix and find 3rd row:
  TX    <- t(X)
  ProdX <- solve(TX %*% X) %*% TX
  b     <- ProdX[3,]

  # evaluate K.hat, C.hat, and optimum m, as in Reference;
  # use first L non-zero frequencies:
  K.hat <- as.numeric(b %*% log(spec[1:L]))
  if (K.hat == 0)
    stop("Method fails: K.hat=0")
  C.hat <- ((27/(128 * pi^2))^0.2) * (K.hat^2)^-0.2

  round(C.hat * nr^0.8)
}


#' Standardised Dispersion Analysis
#'
#' Adaptation from discontinued package `fractal` by William Constantine and Donald Percival.
#'
#' @param x x
#' @param front front
#'
#' @return output
#' @export
#'
#' @keywords internal
#'
SDA <- function(x, front=FALSE){

  n.sample <- length(x)
  if (n.sample < 32)
    stop("Time series must contain at least 32 points")

  scale <- as.integer(scale_log(scale.min=1, scale.max=n.sample, scale.ratio=2))
  sd <- unlist(lapply(scale, function(m, x, n.sample, front){
    nscale   <- n.sample %/% m
    dyad     <- nscale * m
    n.offset <- splus2R::ifelse1(front, n.sample - dyad, 0)
    plyr::colwise(ts_sd)(as.data.frame(aggregate_bins(x[seq(1 + n.offset, length=dyad)], by=m, FUN=mean)))
  },
  x=x,
  n.sample=n.sample,
  front=front))

  list(sd=as.numeric(sd), scale=scale)
}


#' scale_log
#'
#' Adaptation from discontinued package `ifultools` by William Constantine.
#'
#' @param scale.min min
#' @param scale.max max
#' @param scale.ratio ratio
#' @param scale.res reolution
#' @param coerce coerce?
#'
#' @return output
#' @export
#' @keywords internal
#'
scale_log <- function(scale.min, scale.max, scale.ratio=2, scale.res=NULL, coerce=NULL){
  # check input arguments
  if (!is.null(scale.res)){
    scale.ratio <- 1 / scale.res + 1
  }

  # check input arguments
  if (scale.ratio == 1){stop("scale.ratio cannot be unity")}
  if (scale.ratio <= 0){stop("scale.ratio must be positive")}

  n.scale <- splus2R::ifelse1(scale.ratio > 1, ceiling(log(scale.max/scale.min)/log(scale.ratio)),
                              ceiling(log(scale.min/scale.max)/log(scale.ratio)))

  scale    <- vector("numeric", length=n.scale)
  scale[1] <- splus2R::ifelse1(scale.ratio > 1, scale.min, scale.max)

  for (i in seq(2, n.scale)){
    scale[i] <- scale[i - 1] * scale.ratio
  }

  if (!is.null(coerce) && is.function(coerce))
    scale <- coerce(scale)

  unique(scale)
}


#' aggregate_bins
#'
#' Adaptation from discontinued package `ifultools` by William Constantine.
#'
#' @param x x
#' @param by by
#' @param FUN FUNction
#' @param moving moving
#' @param ... additional
#'
#' @return output
#' @export
#' @keywords internal
#'
aggregate_bins <- function(x, by, FUN, moving=FALSE, ...){

  x <- as.numeric(x)
  if(!by%[]%c(1,length(x))){
    stop("Range of x not on specified interval")
  }
  if (!is.function(FUN)){
    stop("FUN must be a function")
  }

  n.sample <- length(x)

  if (!moving){
    n.usable <- (n.sample %/% by) * by

    z <- as.vector(apply(matrix(x[1:n.usable], nrow=by), MARGIN=2, FUN=FUN, ...))
    if (n.usable != n.sample){z <- c(z, as.vector(FUN(x[(n.usable+1):n.sample], ...)))}
  }else{

    moving  <- as.integer(moving)
    if(!moving%[]%c(1, n.sample)){
      stop("Range of x not on specified interval")
    }
    index   <- seq(moving)
    n.group <- as.integer(floor((n.sample - moving)/ by + 1))
    last    <- (n.group - 1)*by + moving

    if (!n.group){
      stop("Moving window not possible with input parameters")
    }

    z <- vector("numeric", length=n.group)

    for (i in seq(n.group)){
      z[i] <- FUN(x[index + by*(i-1)], ...)
    }

    if (last < n.sample){
      z <- c(z, as.vector(FUN(x[(last+1):n.sample], ...)))
    }
  }

  return(z)
}

