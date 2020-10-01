
# Matrix operations ----

#' @title Distance to binary matrix
#'
#' @description Distance matrix to binary matrix based on threshold value
#'
#' @param distmat Distance matrix
#' @param emRad The radius or threshold value
#' @param theiler Use a theiler window around the line of identity / synchronisation to remove high auto-correlation at short time-lags (default = `0`)
#' @param convMat Should the matrix be converted from a `distmat` object of class [Matrix::Matrix()] to [base::matrix()] (or vice versa)
#'
#' @return A (sparse) matrix with only 0s and 1s
#'
#' @export
#'
#' @family Distance matrix operations (recurrence plot)
#' @family Distance matrix operations (recurrence network)
#'
di2bi <- function(distmat, emRad, theiler = 0, convMat = FALSE){

  matPack <- FALSE
  # if already Matrix do not convert to matrix
  if(any(grepl("Matrix",class(distmat)))){
    matPack <- TRUE
    convMat <- TRUE
  }

  if(is.na(theiler%00%NA)){
    if(!is.null(attributes(distmat)$theiler)){
      message(paste0("Value found in attribute 'theiler'... assuming a theiler window of size:",attributes(distmat)$theiler,"was already removed."))
      #theiler <- attr(RM,"theiler")
    }
    theiler <- 0
  } else {
    if(theiler < 0){
      theiler <- 0
    } else {
      if(length(diag(distmat))<length(-theiler:theiler)){
        message("Ignoring theiler window larger than matrix...")
        theiler <- 0
      } else {
        attr(distmat,"theiler") <- theiler
      }
    }
  }

  distmat <- bandReplace(distmat,-theiler,theiler,0)


  # RP <- matrix(0,dim(distmat)[1],dim(distmat)[2])
  # RP[as.matrix(distmat <= emRad)] <- 1
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

  return(RP)
}


#' Distance 2 weighted matrix
#'
#' Distance matrix to weighted matrix based on threshold value
#'
#' @param distmat Distance matrix
#' @param emRad The radius or threshold value
#' @param theiler Use a theiler window around the line of identity / synchronisation to remove high auto-correlation at short time-lags (default = `0`)
#' @param convMat convMat Should the matrix be converted from a `distmat` obkect of class [Matrix::Matrix()] to [base::matrix()] (or vice versa)
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

  if(is.na(theiler%00%NA)){
    if(!is.null(attributes(distmat)$theiler)){
      message(paste0("Value found in attribute 'theiler'... assuming a theiler window of size:",attributes(distmat)$theiler,"was already removed."))
      #theiler <- attr(RM,"theiler")
    }
    theiler <- 0
  } else {
    if(theiler <= 0){
      theiler <- 0
    } else {
      if(length(diag(distmat))<length(-theiler:theiler)){
        message("Ignoring theiler window larger than matrix...")
        theiler <- 0
      } else {
        attr(distmat,"theiler") <- theiler
      }
    }
  }

  distmat <- bandReplace(distmat,-theiler,theiler,0)

  # RP <- NetComp::matrix_threshold(distmat,threshold = emRad, minval = 1, maxval = 0)
  if(emRad==0) emRad <- .Machine$double.eps
  # RP <- distmat #matrix(0,dim(distmat)[1],dim(distmat)[2])
  # RP[distmat <= emRad] <- 0

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

  return(RP)
}


#' Matrix to indexed data frame
#'
#' Mimics the default behaviour of [reshape2::melt]
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
  if(lower>0){lower=-1*lower
  warning("lower > 0 ...\n using: -1*lower")
  }
  if(upper<0){upper=abs(upper)
  warning("upper > 0 ...\n using: abs(upper)")
  }

  if(all(lower==0,upper==0)){
    #diag(mat) <- value
    if(!silent){message(paste0("lower and upper are both 0 (no band, just diagonal)\n using: diag(mat) <- ",round(value,4),"..."))}
  }

  delta <- col(mat)-row(mat)
  mat[delta >= lower & delta <= upper] <- value

  return(mat)
}
