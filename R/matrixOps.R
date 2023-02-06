
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
mat_di2bi <- function(distmat, emRad = NA, convMat = FALSE){

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
    if(!all(stats::na.exclude(as.vector(RP))%in%c(0,1))){warning("Matrix did not convert to a binary (0,1) matrix!!")}

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
mat_di2we <- function(distmat, emRad, theiler = 0, convMat = FALSE){

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
mat_di2ch <- function(distmat, y, emRad, theiler = 0, convMat = FALSE){

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
#' mat_mat2ind(as.matrix(1:100,ncol=10))
#'
mat_mat2ind <- function(mat){
  if(all(!is.matrix(mat),!attr(class(mat),"package")%in%"Matrix")){stop("Input has to be a matrix.")}
  ind <- which(!is.na(mat), arr.ind = TRUE)
  out <- cbind.data.frame(ind, mat[!is.na(mat)])
  colnames(out) <- c("Var1","Var2","value")
  return(out)
}



#' Get indices of matrix diagonals, rows, or columns
#'
#' @param Xlength X dim
#' @param Ylength Y dim
#' @param index index of diagonal, row or column
#' @param diagonal diagonal
#' @param horizontal horizontal
#'
#' @return list
#'
#' @export
#'
#' @examples
#'
mat_ind <- function(Xlength, Ylength, index, diagonal = FALSE){

  if(diagonal){
    if(index>=0){
      Xx <- seq(1,(Xlength-index))
      Yy <- seq((index+1),Ylength)
    } else {
      index <- abs(index)
      # Yy <- seq(1,(Xlength-index))
      # Xx <- seq((index+1),Ylength)
      Yy <- seq(1,(Xlength-index))
      Xx <- seq((index+1),Ylength)
    }
    return(list(r = Xx, c = Yy))
  } else {
   # if(symmetrical){
      Yy <- seq(1,Ylength)
      Xx <- index
  #  } else {
  #     Yy <- seq(1,Ylength)
  #     Xx <- index
  #   }

    # if(horizontal){
    #   return(list(r = Yy, c = Xx))
    # } else {
      return(list(r = Xx, c = Yy))
   # }
  }
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
mat_hamming <- function(X, Y=NULL, embedded=TRUE) {
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
setTheiler <- function(RM, theiler = NA, silent = FALSE, chromatic = FALSE){

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

  if(length(theiler)==1){
    if(theiler < 0){theiler <- 0}
    if(length(Matrix::diag(RM))<length(-theiler:theiler)){
      if(!silent){message("Ignoring theiler window, it is larger than the matrix...")}
      skip <- TRUE
    }
    if(theiler == 0){skip <- TRUE}
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

    if(all(as.vector(RM)==0|as.vector(RM)==1)|chromatic){
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




#' Course grain a matrix for plotting
#'
#'
#' @param RM A (recurrence) matrix
#' @param target_height How many rows? (default = `NROW(RM)/2`)
#' @param target_width How many columns? (default = `NCOL(RM)/2`)
#' @param summary_func How to summarise the values in subset `X` of `RM`. If set to `NA`, the function will try to pick a summary function based on the cell values: If `RM` is a distance matrix,  `mean(X, na.rm = TRUE)` will be used; If it is a binary matrix `ifelse(mean(X, na.rm = TRUE)>recurrence_threshold,1,0)`, a categorical matrix (`categorical = TRUE`, or, matrix attribute `chromatic = TRUE`) will pick the most frequent category in the subset `attributes(ftable(X))$col.vars$x[[which.max(ftable(X))]]`. (default = `NA`)
#' @param recurrence_threshold For a binary matrix the mean of the cells to be summarised will vary between `0` and `1`, which essentially represents the recurrence rate for that subset of the matrix. If `NA` the threshold will be set to a value that in most cases should return a plot with a similar `RR` as the original plot. (default = `NA`)
#' @param categorical If set to `TRUE`, will force `summary_func` to select the most frequent value. If `NA` the matrix attribute `chromatic` will be used. If `chromatic` is not present, all values in the matrix have to be whole numbers as determined by `plyr::is.discrete()`. (default = `NA`)
#' @param output_type The output format for `plyr::vapply()`. (default = `0.0`)
#' @param n_core Number of cores for parallel processing. Set to `NA` to automatically choose cores. (default = `1`)
#' @param silent Silt-ish mode (default = `FALSE`)
#'
#' @note This code was inspired by code published in a blog post by Guillaume Devailly on 29-04-2020 (https://gdevailly.netlify.app/post/plotting-big-matrices-in-r/)
#'
#' @return A coursegrained matrix of size `target_width` by `target_height`.
#'
#' @export
#'
#' @examples
#'
#' # Continuous
#' RMc1 <- rp(cumsum(rnorm(200)))
#' rp_plot(RMc1)
#' RMc2 <- mat_coursegrain(RMc1)
#' rp_plot(RMc2)
#'
#' # Binary
#' RMb1 <- rp(cumsum(rnorm(200)), emRad = NA)
#' rp_plot(RMb1, plotMeasures = TRUE)
#' # Reported RQA measures in rp_plot will be based on the full matrix
#' rp_plot(RMb1, maxSize = 100^2, plotMeasures = TRUE)
#' # Plotting the coursegrained matrix itself will yield different values
#' RMb2 <- mat_coursegrain(RMb1)
#' rp_plot(RMb2, plotMeasures = TRUE)
#'
#' # Categorical
#' RMl1 <- rp(y1 = round(runif(100, min = 1, max = 3)), chromatic = TRUE)
#' rp_plot(RMl1)
#' RMl2 <- mat_coursegrain(RMl1, categorical = TRUE)
#' rp_plot(RMl2)
#'
mat_coursegrain <- function(RM,
                               target_height = round(NROW(RM)/2),
                               target_width = round(NCOL(RM)/2),
                               summary_func = NA,
                               recurrence_threshold = NA,
                               categorical = NA,
                               output_type = 0.0, #vapply style
                               n_core = 1, # parallel processing
                               silent = FALSE){

  # square?
  if(NROW(RM)==NCOL(RM)){
    if(target_height!=target_width){
      stop("Original matrix is square, but target dimensions are not!")
    }
  }

  if(target_height > NROW(RM) | target_width > NCOL(RM)){
    stop("Input matrix must be bigger than target width and height.")
  }

  if(is.na(recurrence_threshold)){
    if(all(stats::na.exclude(as.vector(RM))%in%c(0,1))){
      RR <- rp_measures(RM)
      recurrence_threshold <- mean(c(RR$RR,RR$SING_rate,RR$DET,RR$LAM_hv), na.rm = TRUE)
      rm(RR)
      #((NROW(RM)*NCOL(RM))/(target_height*target_width))
    }
  } else {
    if(recurrence_threshold%)(%c(0,1)){
      recurrence_threshold <- mean(RM, na.rm = TRUE)
    }
  }

  if(is.na(categorical)){
    if(attributes(RM)$chromatic|all(plyr::is.discrete(RM))){
      categorical <- TRUE
    } else {
      categorical <- FALSE
    }
  }

  if(is.na(n_core)){
    n_core <- parallel::detectCores()-1
  }
  if(n_core>1){
    available_core <- parallel::detectCores()
    if(!n_core%[]%c(2,available_core)){
      n_core <- (available_core-1)
    }
  }

  if(is.na(summary_func)){
    if(all(stats::na.exclude(as.vector(RM))%in%c(0,1))){
      summary_func <- function(x){ifelse(mean(x, na.rm = TRUE)>recurrence_threshold,1,0)}
      if(!silent){message("Binary matrix... using summary function 'ifelse(mean(x, na.rm = TRUE)>recurrence_threshold,1,0)' for coursegraining.")}
    } else {
      if(categorical){
        summary_func <- function(x){as.numeric(attributes(stats::ftable(x))$col.vars$x[[which.max(stats::ftable(x))]])}
        if(!silent){message("Categorical matrix... using summary function 'attributes(ftable(x))$col.vars[[which.max(ftable(x))]]' for coursegraining.")}
      } else {
        if(!all(plyr::is.discrete(RM))){
          summary_func <- function(x){mean(x, na.rm = TRUE)}
          if(!silent){message("Continuous matrix... using summary function 'mean(x, na.rm = TRUE)' for coursegraining.")}
        }
      }
    }
  }
  tmpMat <- as.matrix(RM)

  seq_height <- round(seq(1, NROW(RM), length.out = target_height + 1))
  seq_width  <- round(seq(1, NCOL(RM), length.out = target_width  + 1))

  # subMat <- list()
  # m <- 0
  # for(i in seq_len(target_height)){
  #   for(j in seq_len(target_width)){
  #     m <- m + 1
  #     subMat[[m]] <- list(x = seq(seq_height[i], seq_height[i + 1]),
  #                         y = seq(seq_width[j], seq_width[j + 1]))
  #   }
  #  }
  #
  # tmpMat <- matrix(plyr::laply(subMat, function(f) summary_func(as.numeric(RM[f$x,f$y]))),nrow = target_height, ncol = target_width)

#  rp_plot(tmpMat, courseGrain = FALSE)

  # complicated way to write a double for loop
  tmpMat <- do.call(rbind, parallel::mclapply(seq_len(target_height), function(i) { # i is row
    vapply(seq_len(target_width), function(j) { # j is column
      summary_func(as.numeric(RM[seq(seq_height[i], seq_height[i + 1]), seq(seq_width[j] , seq_width[j + 1])]))
    }, output_type)
  }, mc.cores = n_core))


  # tmpMat <- purrr::map(seq_len(target_height), ~ purrr::map(seq_len(target_width), ~ summary_func(as.numeric(RM[seq(seq_height[.x], seq_height[.x + 1]), seq(seq_width[.y] , seq_width[.y + 1])])), .y=.x))
  #

  tmpMat <- Matrix::as.matrix(tmpMat)
  tmpMat <- rp_copy_attributes(source = RM, target = tmpMat)
  tmpMat <- rp_checkfix(tmpMat, checkS4 = TRUE, fixS4 = TRUE, checkAUTO = TRUE, fixAUTO = TRUE)
  rm(RM)

  return(tmpMat)
}


