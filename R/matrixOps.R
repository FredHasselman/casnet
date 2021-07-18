
# Matrix operations ----

#' @title Distance to binary matrix
#'
#' @description Distance matrix to binary matrix based on threshold value
#'
#' @param distmat Distance matrix
#' @param emRad The radius or threshold value
#' @param convMat Should the matrix be converted from a `distmat` object of class [Matrix::Matrix()] to [base::matrix()] (or vice versa)
#'
#' @return A (sparse) matrix with only 0s and 1s
#'
#' @export
#'
#' @family Distance matrix operations (recurrence plot)
#' @family Distance matrix operations (recurrence network)
#'
di2bi <- function(distmat, emRad = NA, convMat = FALSE){

  matPack <- FALSE
  # if already Matrix do not convert to matrix
  if(any(grepl("Matrix",class(distmat)))){
    matPack <- TRUE
    convMat <- TRUE
  }

  distmat <- rp_checkfix(distmat, checkAUTO = TRUE, fixAUTO = TRUE)
  #attributes(distmat)

  if(any(is.na(distmat))){
    NAij  <- Matrix::which(is.na(distmat), arr.ind=TRUE)
  } else {
    NAij <- NA
  }

  # RP <- matrix(0,dim(distmat)[1],dim(distmat)[2])
  # RP[as.matrix(distmat <= emRad)] <- 1
  if(is.na(emRad)){emRad <- est_radius(distmat)$Radius}
  if(emRad==0){emRad <- .Machine$double.eps}
  # Always use sparse representation for conversion to save memory load
  ij  <- Matrix::which(distmat <= emRad, arr.ind=TRUE)

  if(NROW(ij)>0){

    xij <- data.frame(y =  sapply(seq_along(ij[,1]),function(r){distmat[ij[[r,1]],ij[[r,2]]]}), ij)
    suppressMessages(RP <- Matrix::sparseMatrix(x=rep(1,length(xij$y)),i=xij$row,j=xij$col, dims = dim(distmat)))

    # Simple check
    if(!all(as.vector(RP)==0|as.vector(RP)==1)){warning("Matrix did not convert to a binary (0,1) matrix!!")}

  } else {

    RP <- matrix(0,dim(distmat)[1],dim(distmat)[2])
  }

  if(convMat&matPack){RP <- Matrix::as.matrix(RP)}

  suppressMessages(RP <- rp_copy_attributes(source = distmat,  target = RP))
  attributes(RP)$emRad <- emRad
  attributes(RP)$NAij  <- NAij

  return(RP)
}


#' Distance 2 weighted matrix
#'
#' Distance matrix to weighted matrix based on threshold value
#'
#' @param distmat Distance matrix
#' @param emRad The radius or threshold value
#' @param theiler Use a theiler window around the line of identity / synchronisation to remove high auto-correlation at short time-lags (default = `0`)
#' @param convMat convMat Should the matrix be converted from a `distmat` object of class [Matrix::Matrix()] to [base::matrix()] (or vice versa)
#'
#' @return A matrix with 0s and values < threshold distance value
#'
#' @export
#'
#' @family Distance matrix operations (recurrence plot)
#' @family Distance matrix operations (recurrence network)
#'
di2we <- function(distmat, emRad, theiler = 0, convMat = FALSE){

  matPack <- FALSE
  if(any(grepl("Matrix",class(distmat)))){
    matPack <- TRUE
    convMat <- TRUE
  }

  # RP <- NetComp::matrix_threshold(distmat,threshold = emRad, minval = 1, maxval = 0)
  if(emRad==0) emRad <- .Machine$double.eps
  # RP <- distmat #matrix(0,dim(distmat)[1],dim(distmat)[2])
  # RP[distmat <= emRad] <- 0

  if(any(is.na(distmat))){
    NAij  <- Matrix::which(is.na(distmat), arr.ind=TRUE)
  } else {
    NAij <- NA
  }

  ij  <- Matrix::which(distmat <= emRad, arr.ind=TRUE)

  if(NROW(ij)>0){

    # Always use sparse representation for conversion to save memory load
    #xij <- data.frame(y =  sapply(which(distmat > emRad, arr.ind=TRUE)[,1],function(r){distmat[ij[[r,1]],ij[[r,2]]]}), which(distmat > emRad, arr.ind=TRUE))
    xij <- data.frame(y =  sapply(seq_along(ij[,1]),function(r){distmat[ij[[r,1]],ij[[r,2]]]}), ij)

    suppressWarnings(RP <- Matrix::sparseMatrix(x=xij$y,i=xij$row,j=xij$col, dims = dim(distmat)))

    #  if(!all(as.vector(RP)==0|as.vector(RP)==1)){warning("Matrix did not convert to a binary (0,1) matrix!!")}
  } else {

    RP <- matrix(0,dim(distmat)[1],dim(distmat)[2])
  }

  if(convMat&matPack){RP <- Matrix::as.matrix(RP)}

  RP <- rp_copy_attributes(source = distmat,  target = RP)
  attributes(RP)$emRad <- emRad
  attributes(RP)$NAij <- NAij

  return(RP)
}




#' Distance 2 chromatic matrix
#'
#' Distance matrix to chromatic matrix based on unordered categorical series
#'
#' @param distmat Distance matrix
#' @param y One of the dimensions (as a data frame or matrix) of the RP which must contain unique unordered categorical values
#' @param emRad The radius or threshold value
#' @param theiler Use a theiler window around the line of identity / synchronisation to remove high auto-correlation at short time-lags (default = `0`)
#' @param convMat convMat Should the matrix be converted from a `distmat` object of class [Matrix::Matrix()] to [base::matrix()] (or vice versa)
#'
#' @return A matrix with 0s and the unordered categorical values that are recurring
#'
#' @export
#'
#' @family Distance matrix operations (recurrence plot)
#' @family Distance matrix operations (recurrence network)
#'
di2ch <- function(distmat, y, emRad, theiler = 0, convMat = FALSE){

  matPack <- FALSE
  if(any(grepl("Matrix",class(distmat)))){
    matPack <- TRUE
    convMat <- TRUE
  }

  if(any(is.na(distmat))){
    NAij  <- Matrix::which(is.na(distmat), arr.ind=TRUE)
  } else {
    NAij <- NA
  }

  if(emRad==0) emRad <- .Machine$double.eps

  ij  <- Matrix::which(distmat <= emRad, arr.ind=TRUE)

  if(NROW(ij)>0){

    # Always use sparse representation for conversion to save memory load
    #xij <- data.frame(y =  sapply(which(distmat > emRad, arr.ind=TRUE)[,1],function(r){distmat[ij[[r,1]],ij[[r,2]]]}), which(distmat > emRad, arr.ind=TRUE))
    xij <- data.frame(y =  sapply(seq_along(ij[,1]),function(r){y[ij[r,1]]}), ij)

    suppressWarnings(RP <- Matrix::sparseMatrix(x=xij$y,i=xij$row,j=xij$col, dims = dim(distmat)))


    #  if(!all(as.vector(RP)==0|as.vector(RP)==1)){warning("Matrix did not convert to a binary (0,1) matrix!!")}


  } else {

    RP <- matrix(0,dim(distmat)[1],dim(distmat)[2])
  }

  if(convMat&matPack){RP <- Matrix::as.matrix(RP)}

  RP <- rp_copy_attributes(source = distmat,  target = RP)
  attributes(RP)$emRad <- emRad
  attributes(RP)$NAij <- NAij

  return(RP)
}




#' Matrix to indexed data frame
#'
#' Mimics the default behaviour of `reshape2::melt()`
#'
#' @param mat A matrix
#'
#' @return A data.frame with two index columns named `"Var1"` and `"Var2"`
#'
#' @export
#'
#' @examples
#'
#' mat2ind(as.matrix(1:100,ncol=10))
#'
mat2ind <- function(mat){
  if(all(!is.matrix(mat),!attr(class(mat),"package")%in%"Matrix")){stop("Input has to be a matrix.")}
  ind <- which(!is.na(mat), arr.ind = TRUE)
  out <- cbind.data.frame(ind, mat[!is.na(mat)])
  colnames(out) <- c("Var1","Var2","value")
  return(out)
}



#' Calculate Hamming distance
#'
#' @param X A matrix (of coordinates)
#' @param Y A matrix (of coordinates)
#' @param embedded Do X and/or Y represent surrogate dimensions of an embedded time series?
#'
#' @return A hamming-distance matrix of X, or X and Y. Useful for ordered and unordered categorical data.
#' @export
#'
#' @family Distance matrix operations (recurrence plot)
#' @author Fred Hasselman
#'
dist_hamming <- function(X, Y=NULL, embedded=TRUE) {
  if ( missing(Y) ) {
    if(embedded){
      # X and Y represent delay-embeddings of a timeseries
      if(which.max(dim(X))==1){X <- t(X)}
    }
    uniqs <- unique(as.vector(X))
    U <- X == uniqs[1]
    H <- t(U) %*% U
    for ( uniq in uniqs[-1] ) {
      U <- X == uniq
      H <- H + t(U) %*% U
    }
  } else {
    if(embedded){
      # X and Y represent delay-embeddings of a timeseries
      if(which.max(dim(X))==1){X <- t(X)}
      if(which.max(dim(Y))==1){Y <- t(Y)}
    }
    uniqs <- union(X, Y)
    H <- t(X == uniqs[1]) %*% (Y == uniqs[1])
    for ( uniq in uniqs[-1] ) {
      H <- H + t(X == uniq) %*% (Y == uniq)
    }
  }
  NROW(X) - H
}


#' Replace matrix diagonals
#'
#' Sets a band of matrix diagonals to any given value
#'
#' @param mat A Matrix
#' @param lower Lower diagonal to be included in the band (should be \eqn{\le 0})
#' @param upper Upper diagonal to be included in the band (should be \eqn{\ge 0})
#' @param value A single value to replace all values in the selected band (default = `NA`)
#' @param silent Operate in silence, only (some) warnings will be shown (default = `TRUE`)
#'
#' @return A matrix in which the values in the selected diagonals have been replaced
#'
#' @export
#'
#' @family Distance matrix operations (recurrence plot)
#'
#' @author Fred Hasselman
#'
#'
#' @examples
#' # Create a 10 by 10 matrix
#' library(Matrix)
#' m <- Matrix(rnorm(10),10,10)
#'
#' bandReplace(m,-1,1,0)   # Replace diagonal and adjacent bands with 0 (Theiler window of 1)
bandReplace <- function(mat, lower, upper, value = NA, silent=TRUE){

  # if(lower>0){lower=-1*lower
  # warning("lower > 0 ...\n using: -1*lower")
  # }
  # if(upper<0){upper=abs(upper)
  # warning("upper > 0 ...\n using: abs(upper)")
  # }

  # if(all(lower==0,upper==0)){
  #   #diag(mat) <- value
  #   if(!silent){message(paste0("lower and upper are both 0 (no band, just diagonal)\n using: diag(mat) <- ",round(value,4),"..."))}
  # }

  tmp <- mat
  delta <- col(mat)-row(mat)
  indc  <- delta >= lower & delta <= upper
  suppressMessages(mat[indc] <- value)
  mat <- methods::as(mat, "dgCMatrix")
  mat <- rp_copy_attributes(source = tmp, target = mat, source_remove = c("names", "row.names", "class","dim", "dimnames","x","i","p"))
  rm(tmp)

  return(mat)
}



#' Set theiler window on a distance matrix or recurrence matrix.
#'
#' @inheritParams rp
#' @inheritParams rp_measures
#'
#' @return The matrix with the diagonals indicated in the `theiler` argument set to either `max(RM)+1` (if `RM` is a distance matrix) or `0` (if `RM` is a recurrence matrix).
#'
#' @export
#'
setTheiler <- function(RM, theiler = NA, silent = FALSE){

  checkPkg("Matrix")
  # Check auto-recurrence
  RM <- rp_checkfix(RM, checkAUTO = TRUE, fixAUTO = TRUE)

  skip <- FALSE

  # Theiler
  if(!is.na(attributes(RM)$theiler%00%NA)){
    if(is.numeric(attributes(RM)$theiler)){
      if(!silent){message(paste0("Value found in attribute 'theiler'... assuming a theiler window of size: ",attributes(RM)$theiler," was already removed."))}
      skip <- TRUE
    } else {
      if(!silent){message(paste0("Value found in attribute 'theiler' is not numeric (",attributes(RM)$theiler,"), setting: 'theiler <- NA'"))}
      theiler <- NA
    }
  } else {
    if(!is.numeric(theiler)){
      if(!is.na(theiler)){
        if(!silent){message(paste0("Value passed to 'theiler' is not numeric (",attributes(RM)$theiler,"), setting: 'theiler <- NA'"))}
      theiler <- NA
    }
   }
  }

  if(is.na(theiler)){
    if(attributes(RM)$AUTO){
      theiler <- 1
    } else {
      theiler <- 0
    }
  }

  if(length(theiler==1)){
    if(theiler < 0){theiler <- 0}
    if(length(Matrix::diag(RM))<length(-theiler:theiler)){
      if(!silent){message("Ignoring theiler window, it is larger than the matrix...")}
      skip <- TRUE
    }
  }

  if(length(theiler==2)){
    if(length(Matrix::diag(RM))<length(min(theiler):max(theiler))){
      if(!silent){message("Ignoring theiler window, it is larger than the matrix...")}
      skip <- TRUE
    }
  }

  if(length(theiler>2)){
    if(any(length(Matrix::diag(RM))<abs(theiler))){
      if(!silent){message("Ignoring theiler window, it contains diagonals that are larger than the matrix...")}
      skip <- TRUE
    }
  }

  if(!skip){

    if(all(as.vector(RM)==0|as.vector(RM)==1)){
      value <- 0
    } else {
      value <- max(Matrix::as.matrix(RM), na.rm = TRUE)+1
    }

  #LOS <- diag(RM)

    if(length(theiler)==1){
      if(theiler > 1){
        theiler <- (theiler-1)
        RM <- bandReplace(mat = RM, lower = -theiler, upper = theiler, value = value)
      } else {
        if(theiler != 0){
          RM <- bandReplace(mat = RM, lower = 0, upper = 0, value = value)
        }
      }
    }

  if(length(theiler)==2){
      RM <- bandReplace(mat = RM, lower = min(theiler), upper = max(theiler), value = value)
    }

    if(length(theiler)>2){
      theiler <- sort(theiler)
      for(d in seq_along(theiler)){
      RM <- bandReplace(mat = RM, lower = theiler[d], upper = theiler[d], value = value)
      }
    }

    attr(RM,"theiler") <- theiler

  } # skip

  return(RM)
}
