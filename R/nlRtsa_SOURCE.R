#' @title Distance to binary matrix
#'
#' @description Distance matrix to binary matrix based on threshold value
#'
#' @param distmat Distance matrix
#' @param emRad The radius or threshold value
#' @param theiler = Use a theiler window around the line of identity / synchronisation to remove high auto-correlation at short time-lags (default = `0`)
#' @param convMat Should the matrix be converted from a `distmat` obkect of class [Matrix::Matrix()] to [base::matrix()] (or vice versa)
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

  # RP <- matrix(0,dim(distmat)[1],dim(distmat)[2])
  # RP[as.matrix(distmat <= emRad)] <- 1
  if(emRad==0){emRad <- .Machine$double.eps}
  # Always use sparse representation for conversion to save memory load
  ij  <- Matrix::which(distmat <= emRad, arr.ind=TRUE)

  if(NROW(ij)>0){

    xij <- data.frame(y =  sapply(seq_along(ij[,1]),function(r){distmat[ij[[r,1]],ij[[r,2]]]}), ij)
    suppressWarnings(RP <- Matrix::sparseMatrix(x=rep(1,length(xij$y)),i=xij$row,j=xij$col, dims = dim(distmat)))

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
#' @param convMat convMat Should the matrix be converted from a `distmat` obkect of class [Matrix::Matrix()] to [base::matrix()] (or vice versa)
#'
#' @return A matrix with 0s and values < threshold distance value
#'
#' @export
#'
#' @family Distance matrix operations (recurrence plot)
#' @family Distance matrix operations (recurrence network)
#'
di2we <- function(distmat, emRad, convMat = FALSE){

  matPack <- FALSE
  if(any(grepl("Matrix",class(distmat)))){
    matPack <- TRUE
    convMat <- TRUE
  }

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



# Complex Networks ---


SWtest0 <- function(g){
  Nreps <- 10;
  histr  <- vector("integer",Nreps)
  target<- round(mean(igraph::degree(g)))
  now   <- target/2
  for(i in 1:Nreps){
    gt      <- igraph::watts.strogatz.game(dim=1, size=length(igraph::degree(g)), nei=now, 0)
    histr[i] <- round(mean(igraph::degree(gt)))
    ifelse(histr[i] %in% histr,break,{
      ifelse(histr[i]>target,{now<-now-1},{
        ifelse(histr[i]<target,{now<-now+1},{
          break})
      })
    })
  }
  return(gt)
}


# SWtestV <- function(g,N){
#  return(list(cp=igraph::transitivity(g,type="global"),cpR=igraph::transitivity(igraph::rewire(g,mode=c("simple"),niter=N),type="global"),lp=igraph::average.path.length(g), lpR=igraph::average.path.length(igraph::rewire(g,mode=c("simple"),niter=N))))
# }

#' Small World test
#'
#' @param g An igraph object
#' @param p p
#' @param N N
#'
#' @export
#'
SWtestE <- function(g,p=1,N=20){
  values <- matrix(nrow=N,ncol=6,dimnames=list(c(1:N),c("cp","cpR","cp0","lp","lpR","lp0")))

  for(n in 1:N) {
    gt<-SWtest0(g)
    values[n,] <- c(igraph::transitivity(g,type="localaverage"),igraph::transitivity(igraph::rewire(g,igraph::each_edge(p=p)),type="localaverage"),igraph::transitivity(gt,type="localaverage"),igraph::average.path.length(g),igraph::average.path.length(igraph::rewire(g,igraph::each_edge(p=p))),igraph::average.path.length(gt))}
  values[n,values[n,]==0] <- NA #values[n,values[n,]==0]+1e-8}

  values   <- cbind(values,(values[,1]/values[,2])/(values[,4]/values[,5]),(values[,1]/values[,3]),(values[,4]/values[,6]),((values[,1]/values[,3])/values[,2])/((values[,4]/values[,6])/values[,5]))
  valuesSD <- data.frame(matrix(apply(values[,1:10],2,FUN = stats::sd,na.rm=TRUE),nrow=1,ncol=10,dimnames=list(c(1),c("cp","cpR","cp0","lp","lpR","lp0","SWI","cp:cp0","lp:lp0","SWIn"))))
  valuesAV <- data.frame(matrix(colMeans(values[,1:10],na.rm=T),nrow=1,ncol=10,dimnames=list(c(1),c("cp","cpR","cp0","lp","lpR","lp0","SWI","cp:cp0","lp:lp0","SWIn"))))
  return(list(valuesAV=valuesAV,valuesSD=valuesSD,valuesSE=valuesSD/sqrt(N)))
}



# PLFsmall <- function(g){
#
#   if(length(igraph::V(g))>100){stop("Vertices > 100, no need to use PLFsmall, use a binning procedure")}
#
#   d <- igraph::degree(g)
#
#   y <- graphics::hist(d,breaks=0.5:(max(d)+0.5),plot=FALSE)$counts
#   if(length(y)<2){
#     warning("Less than 2 points in Log-Log regression... alpha=0")
#     alpha <- 0
#   } else {
#     if(length(y)==2){
#       warning("Caution... Log-Log slope is a bridge (2 points)")
#       chop <- 0
#     } else {
#       chop <- 1
#     }
#     alpha <- stats::coef(stats::lm(rev(log1p(y)[1:(length(y)-chop)]) ~ log1p(1:(length(y)-chop))))[2]
#   }
#
#   return(alpha)
# }


#' Import GridWare files
#'
#' @param gwf_name Name of the GridWare project file. A directory named `../gwf_name_trjs` must be present at the location of the project file.
#' @param delta_t Time between two samples or sampling frequency
#' @param returnOnlyData Just return the data, do not return a list object with data, variable info and preferences.
#' @param saveLongFormat Save the long format trajectory data as a `.csv` file in the same location as `gwf_name`
#'
#' @return A data frame containing State Space Grid trajectories, or a list object with additional info.
#' @export
#'
#' @family State Space Grid functions
#'
ssg_gwf2long <- function(gwf_name, delta_t = 0.01,returnOnlyData = TRUE, saveLongFormat = FALSE){
  gwf_lines <- readr::read_lines(gwf_name)

  var_list_b <- which(grepl("<Config>", gwf_lines, fixed = TRUE))+1
  var_list_e <- which(grepl("MinReturns", gwf_lines, fixed = TRUE))-1

  Nvars   <- length(gwf_lines[var_list_b:var_list_e])
  MaxCols <-  max(plyr::laply(as.list(gwf_lines[var_list_b:var_list_e]), function(li){length(strsplit(li,"\t")[[1]])}))

  var_list <- as.list(gwf_lines[var_list_b:var_list_e])
  varinfo <-   plyr::ldply(var_list, function(li){
    tmp <- strsplit(li,"\t")[[1]]
    dat<- tibble::as.tibble(data.frame(var.name = tmp[3],
                                       var.role = tmp[1],
                                       var.type = tmp[2]))
    if(dat$var.type%in%"integer"|dat$var.role%in%"state"){
      dat$min <- tmp[4]
      dat$max <- tmp[5]
    } else {
      dat$min <- NA
      dat$max <- NA
    }
    if(dat$var.type%in%"categorical"){
      Bvals <- 4
    } else {
      Bvals <- 6
    }
    Evals <- length(tmp)-Bvals
    Nvals <- length(tmp[Bvals:(Bvals+Evals)])
    dat2 <- tibble::as.tibble(matrix(NA,ncol=MaxCols-NCOL(dat), dimnames = list(NULL,paste0("value",1:(MaxCols-NCOL(dat))))))
    if(length(tmp)>3){
      dat2[1,1:Nvals] <- tmp[Bvals:(Bvals+Evals)]
    }
    return(cbind.data.frame(dat,dat2))
  })
  varinfo <- tibble::as.tibble(varinfo)

  conf_list_b <- which(grepl("MinReturns", gwf_lines, fixed = TRUE))
  conf_list_e <- which(grepl("</Config>", gwf_lines, fixed = TRUE))-1
  conf_list <-  readr::read_tsv(paste0(c("setting\tvalue",gwf_lines[conf_list_b:conf_list_e]),collapse = "\n"))

  traj_list_b <- which(grepl("<Trajectories>", gwf_lines, fixed = TRUE))+1
  traj_list_e <- which(grepl("</Trajectories>", gwf_lines, fixed = TRUE))-1
  traj_list <- readr::read_tsv(paste0(gwf_lines[traj_list_b:traj_list_e],collapse = "\n"))

  state_vars <- sum(varinfo$var.role=="state")+1
  dirname <- paste0(gsub(".gwf","",gwf_name,fixed = TRUE),"_trjs")
  if(dir.exists(dirname)){
    trj_data <- dir(dirname, pattern = "([.]trj)$", full.names =TRUE, include.dirs = FALSE)
    names(trj_data) <-  dir(dirname, pattern = "([.]trj)$", full.names =FALSE, include.dirs = FALSE)

    suppressMessages(suppressWarnings(data_long <- plyr::ldply(trj_data, function(p){
      tmp <- readr::read_tsv(p)
      for(l in 1:(NROW(tmp)-1)){
        return(data.frame(time = seq(from = tmp$Onset[l],to = (tmp$Onset[l+1]-delta_t), by = delta_t), tmp[l,-1]))
      }
    }, .id = "Filename"))
    )
    suppressMessages(suppressWarnings(data_long <- dplyr::left_join(x = data_long,traj_list,by="Filename")))


    attr(data_long,"Preferences") <-  conf_list
    attr(data_long,"Variable definitions") <-  varinfo
  }

  if(returnOnlyData){
    return(data_long)
  } else {
    return(list(variable_config = varinfo,
                settings_config = conf_list,
                trajectories_list  = traj_list,
                data_long          = data_long))
  }
}


#' Winnowing procedure for SSG
#'
#' Find attractor states in a State Space Grid using a winnowing procedure.
#'
#' @param durations A data frame obtained by function [ts_duration()]
#' @param screeCut Cutoff based on a scree plot.
#'
#' @return Attractor frame
#' @export
#'
#' @family State Space Grid functions
#'
ssg_winnowing <- function(durations, screeCut){

  durations$duration.time[is.na(durations$duration.time)] <- 0
  winnowing <- durations %>% dplyr::filter(.data$duration.time>0)
  Ncells <- NROW(winnowing)
  winnowingList <- scree <- heterogeneity <- list()
  removed <- 0
  run <- 1
  exp <- sum(winnowing$duration.time, na.rm = TRUE)/NROW(winnowing)
  #obs <- winnowing$duration.time-exp
  heterogeneity[[run]] <- sum((winnowing$duration.time-exp)^2 / exp, na.rm = TRUE) / NROW(winnowing)
  scree[[run]] <- 1
  removed <- Ncells - NROW(winnowing)
  winnowingList[[run]] <- winnowing
  while(removed!=(Ncells-1)){
    run <- run +1

    winnowing <- winnowing %>% dplyr::filter(.data$duration.time>min(winnowing$duration.time, na.rm = TRUE))
    removed <- Ncells - NROW(winnowing)
    exp <- sum(winnowing$duration.time, na.rm = TRUE)/NROW(winnowing)

    heterogeneity[[run]] <- sum((winnowing$duration.time-exp)^2 / exp, na.rm = TRUE) / NROW(winnowing)
    scree[[run]]         <-  (heterogeneity[[run-1]]-heterogeneity[[run]])/heterogeneity[[run-1]]
    winnowingList[[run]] <- winnowing
  }

  useRun <- which(unlist(scree)<.5)[1]

  return(list(
    useRun = which(unlist(scree)<screeCut)[1],
    attractors = winnowingList[[useRun]],
    scree = scree,
    heterogeneity = heterogeneity
  )
  )
}


#' Add expected factor labels to observed values
#'
#'
#'
#' @param observed_Ncat obsN
#' @param observed_labels obsL
#' @param expected_Ncat expN
#' @param expected_labels expL
#' @param varname varname
#'
#' @family State Space Grid functions
#'
#' @return character vector
#' @export
#'
factor_obs_exp <- function(observed_Ncat, observed_labels, expected_Ncat=0, expected_labels="", varname = ""){

  LABS <- ""
  warningMessage <- ""
  if(expected_Ncat>0){

    if(expected_Ncat==observed_Ncat){
      if(!identical(sort(observed_labels),sort(expected_labels))){
        #  warningMessage <- paste0("User defined and observed state labels are different. User defined labels will be used.")
        newN <- sum(expected_labels%in%observed_labels, na.rm = TRUE)
        if(newN < observed_Ncat){
          warningMessage <- paste("Different user defined state labels from those observed. New user defined labels will be added.")
          LABS <- expected_labels[!expected_labels%in%observed_labels]
          LABS0 <- observed_labels[observed_labels%in%expected_labels]
          origin <- "expected"
          code <-  paste0("c('",paste0(LABS0, collapse="','"),"','",paste0(LABS, collapse="','"),"')")
        }
      } else {
        LABS <- observed_labels
        origin <- "observed"
        code <- paste0("factor(",varname,", levels = c('",paste0(LABS, collapse="','"),"'))")
      }
    }

    if(expected_Ncat<observed_Ncat){
      warningMessage <- paste("Fewer user defined state labels than observed. Observed labels will be used.")
      LABS <- observed_labels
      origin <- "observed"
      code <- paste0("factor(",varname,", levels = c('",paste0(LABS, collapse="','"),"'))")
    }
    if(expected_Ncat>observed_Ncat){
      newN <- sum(expected_labels%in%observed_labels, na.rm = TRUE)
      if(newN <observed_Ncat){
        warningMessage <- paste("More user defined state labels than observed. New user defined labels will be added.")
        LABS <- expected_labels[!expected_labels%in%observed_labels]
        LABS0 <- observed_labels[observed_labels%in%expected_labels]
        origin <- "expected"
        code <-  paste0("c('",paste0(LABS0, collapse="','"),"','",paste0(LABS, collapse="','"),"')")
      } else {
        warningMessage <- paste("More user defined state labels than observed. User defined labels will be used.")
        LABS <- expected_labels
        origin <- "expected"
        code <- paste0("factor(",varname,", levels = c('",paste0(LABS, collapse="','"),"'))")
      }
    }
  } else {
    LABS <- observed_labels
    origin <- "observed"
    code <- paste0("factor(",varname,", levels = c('",paste0(LABS, collapse="','"),"'))")
  }
  attr(LABS,"varname") <- varname
  attr(LABS,"warningMessage") <- warningMessage
  attr(LABS,"origin") <- origin
  attr(LABS,"code") <- code
  return(LABS)
}



# Toy models ----------

#' Examples of dynamical growth models (maps)
#'
#' Autocatlytic Growth: Iterating differential equations (maps)
#'
#' @param Y0    Initial value.
#' @param r    Growth rate parameter.
#' @param k    Carrying capacity.
#' @param N    Length of the time series.
#' @param type    One of: "driving" (default), "damping", "logistic", "vanGeert1991".
#'
#' @return A timeseries object of length N.
#' @export
#'
#' @author Fred Hasselman
#'
#' @family autocatalytic growth functions
#'
#' @examples
#' # The logistic map in the chaotic regime
#' growth_ac(Y0 = 0.01, r = 4, type = "logistic")
growth_ac <- function(Y0 = 0.01, r = 1, k = 1, N = 100, type = c("driving", "damping", "logistic", "vanGeert")[1]){
  # Create a vector Y of length N, which has value Y0 at Y[1]
  if(N>1){
    Y <- as.numeric(c(Y0, rep(NA,N-2)))
    # Conditional on the value of type ...
    switch(type,
           # Iterate N steps of the difference function with values passed for Y0, k and r.
           driving  = sapply(seq_along(Y), function(t) Y[[t+1]] <<- r * Y[t] ),
           damping  = k + sapply(seq_along(Y), function(t) Y[[t+1]] <<- - r * Y[t]^2 / k),
           logistic = sapply(seq_along(Y), function(t) Y[[t+1]] <<- r * Y[t] * ((k - Y[t]) / k)),
           vanGeert = sapply(seq_along(Y), function(t) Y[[t+1]] <<- Y[t] * (1 + r - r * Y[t] / k))
    )}
  return(stats::ts(Y))
}

#' Examples of conditional dynamical growth models (maps)
#'
#' Conditional Autocatlytic Growth: Iterating differential equations (maps)
#'
#' @param Y0 Initial value
#' @param r Growth rate parameter
#' @param k Carrying capacity
#' @param cond Conditional rules passed as a data.frame of the form: cbind.data.frame(Y = ..., par = ..., val = ...)
#' @param N Length of the time series
#'
#' @export
#'
#' @author Fred Hasselman
#'
#' @family autocatalytic growth functions
#'
#' @examples
#' # Plot with the default settings
#' library(lattice)
#' xyplot(growth_ac_cond())
#'
#' # The function can take a set of conditional rules
#' # and apply them sequentially during the iterations.
#' # The conditional rules are passed as a `data.frame`
#'
#' (cond <- cbind.data.frame(Y = c(0.2, 0.6), par = c("r", "r"), val = c(0.5, 0.1)))
#' xyplot(growth_ac_cond(cond=cond))
#'
#' # Combine a change of `r` and a change of `k`
#'
#' (cond <- cbind.data.frame(Y = c(0.2, 1.99), par = c("r", "k"), val = c(0.5, 3)))
#' xyplot(growth_ac_cond(cond=cond))
#'
#' # A fantasy growth process
#'
#' cond <- cbind.data.frame(Y = c(0.1, 1.99, 1.999, 2.5, 2.9),
#' par = c("r", "k", "r", "r","k"),
#' val = c(0.3, 3, 0.9, 0.1, 1.3))
#'
#' xyplot(growth_ac_cond(cond=cond))
growth_ac_cond <- function(Y0 = 0.01, r = 0.1, k = 2, cond = cbind.data.frame(Y = 0.2, par = "r", val = 2), N = 100){
  # Create a vector Y of length N, which has value Y0 at Y[1]
  Y <- c(Y0, rep(NA, N-1))
  # Iterate N steps of the difference equation with values passed for Y0, k and r.
  cnt <- 1
  for(t in seq_along(Y)){
    # Check if the current value of Y is greater than the threshold for the current conditional rule in cond
    if(Y[t] > cond$Y[cnt]){
      # If the threshold is surpassed, change the parameter settings by evaluating: cond$par = cond$val
      eval(parse(text = paste(cond$par[cnt], "=", cond$val[cnt])))
      # Update the counter if there is another conditional rule in cond
      if(cnt < nrow(cond)){cnt <- cnt + 1}
    }
    # Van Geert growth model
    Y[[t+1]] <- Y[t] * (1 + r - r * Y[t] / k)
  }
  return(stats::ts(Y))
}



# HELPERS ----


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


#' Get some nice colours
#'
#' @param Ncols Number of colours
#'
#' @return A list of colours
#' @export
#'
#' @examples
#'
#' # Plot all available colours
#' df <- expand.grid(x=1:8,y=1:8)[1:58,]
#' plot(df$x,df$y,pch=15,col=getColours(Ncols=58),cex=5, axes = FALSE, ann = FALSE)
#' text(df$x,df$y,paste(1:58),col="white")
#'
getColours <- function(Ncols = 20){
  pl <- c("#A65141","#ABB2C4","#C89284","#7896BC","#E1CAC2","#536489","#ECE3CD","#575E7A","#BEBEC4","#C4AD75","#30261E","#BC964D","#3D434D","#D29E55","#B4B6A7","#9B7342","#E8D69D","#48211A","#EBDAB4","#3D4B68","#E8DCCF","#3E688C","#9D3D2D","#507074","#9A756C","#546959","#93725F","#858563","#E1D3B1","#A6B4BB","#78828C","#AA9C6A","#91A193","#CDB16E","#1B528D","#B8A98F","#6E432D","#4F5B71","#D9B196","#20304F","#827561","#98A39E","#8B4F31","#7E7B65","#C1B6A1","#775E45","#C0B3A2","#5A524A","#BBA27F","#3A3935","#C9BDA3","#626961","#8A4F42","#8D8A77","#947B5E","#5D3C30","#AA8470","#493A2C")

  if(Ncols<=58){
    cols <- pl[1:Ncols]
  } else {
    stop("Currently the max. number of colours is 58.")
  }
  return(cols)
}


#' Generate noise series with power law scaling exponent
#'
#' @param y Time series to use as a 'model'. If specified, `N` will be `N = length(y)`, and the series will be constructed based on `stats::fft(y)`.
#' @param alpha The log-log spectral slope, the scaling exponent. Use `0` for white noise, negative numbers for anti-persistant noises: `-1` for \eqn{\frac{1}{f}} noise, positive numbers for persistent noises, e.g. `1` for blue noise.
#' @param N Length of the time series
#' @param standardise Forces scaling of the output to the range `[-1, 1]`, consequently the power law will not necessarily extend right down to `0Hz`.
#' @param randomPower If `TRUE` phases will be deterministic, uniformly distributed in `[-pi,pi]`. If `FALSE`, the spectrum will be stochastic with a Chi-square distribution. If `y` is not `NULL` this argument will be ignored.
#' @param seed Provide an integer number to set the seed for the random number generator in order to get reproducible results. If `NA` (default) no user defined seed will be set,
#'
#' @return Time series with a power law of alpha.
#' @export
#'
#' @note Adapted from a Matlab script called `powernoise.m` by Max Little. The script contained the following commented text:
#'
#' With no option strings specified, the power spectrum is
#  deterministic, and the phases are uniformly distributed in the range
#  -pi to +pi. The power law extends all the way down to 0Hz (DC)
#  component. By specifying the 'randpower' option string however, the
#  power spectrum will be stochastic with Chi-square distribution. The
#  'normalize' option string forces scaling of the output to the range
#  [-1, 1], consequently the power law will not necessarily extend
#  right down to 0Hz.
#
#  (cc) Max Little, 2008. This software is licensed under the
#
#  Attribution-Share Alike 2.5 Generic Creative Commons license:
#  http://creativecommons.org/licenses/by-sa/2.5/
#  If you use this work, please cite:
#
#  Little MA et al. (2007), "Exploiting nonlinear recurrence and fractal
#  scaling properties for voice disorder detection", Biomed Eng Online, 6:23
#
#  As of 20080323 markup
#  If you use this work, consider saying hi on comp.dsp
#  Dale B. Dalrymple
#'
noise_powerlaw <- function(y = NULL, alpha=-1, N=512, standardise = FALSE, randomPower = FALSE, seed = NA){

  if(alpha%)(%c(-2,2)){
    warning("If alpha is outside [-2,2], result are likely not accurate.")
  }

  if(!is.null(y)){
    N <- length(y)
  }

  if(!is.na(seed)){
    if(is.wholenumber(seed)){
      set.seed(seed)
    }
  }

  alpha <- -1 * alpha

  # Nyquist
  N2 <- floor(N/2)-1
  f  <- 2:(N2+1)
  A2 <- 1/(f^(alpha/2))

  if(!is.null(y)){

    p2 <- stats::fft(y)[f]
    d2 <- A2 * p2

  } else {

    if(!randomPower){
      p2 <- (stats::runif(n = N2)-0.5)*2*pi
      d2 <- A2*exp(1i*p2)
    } else {
      # 20080323 update
      p2 <- complex(real= stats::rnorm(n = N2), imaginary = stats::rnorm(n = N2))
      d2 <- A2*p2
    }
  }
  d <- c(1, d2, 1/((N2+2)^alpha), rev(Conj(d2)))
  x <- Re(stats::fft(d, inverse = TRUE)/length(d))

  if (standardise){
    x <- ((x - min(x))/(max(x) - min(x)) - 0.5) * 2
  }

  return(x)
}

#' Generate fractional Gaussian noise
#'
#' @param H Hurst exponent
#' @param N Length of noise series
#' @param mu Mean
#' @param sigma SD
#'
#' @return fGn
#' @export
#'
noise_fGn <- function(H=0.5, N = 512, mu = NULL, sigma = NULL){

  # Determine whether fGn or fBn should be produced.
  if(H%)[%c(0,1)){stop("H must be in (0,1] for fGn!")}
  fBn = 0

  # Calculate the fGn.
  if (H == 0.5){
    y <- stats::rnorm(N)  # If H=0.5, then fGn is equivalent to white Gaussian noise.
    } else {
    Nfft <- 2^ceiling(log2(2*(N-1)))
    NfftHalf <- round(Nfft/2);
    k <- c(0:NfftHalf, (NfftHalf-1):-1:1)
    Zmag <- 0.5 * ( (k+1)^(2*H) - 2.*k^(2*H) + (abs(k-1))^(2*H) )
    rm(k)
    Zmag <- Re(stats::fft(Zmag))
    if ( any(Zmag < 0) ){
      stop('The fast Fourier transform of the circulant covariance had negative values.')
    }
    Zmag <- sqrt(Zmag);
    # Store N and H values in persistent variables for use during subsequent calls to this function.
    Nlast <- N
    Hlast <- H

    #Z = Zmag*(rnorm(1,Nfft) + i*randn(1,Nfft))
    Z <- Zmag*complex(real= stats::rnorm(n = Nfft), imaginary = stats::rnorm(n = Nfft))
    y <- Re(stats::fft(Z,inverse = TRUE)/length(Z)) * sqrt(Nfft)
    rm(Z)
    y[1:N] <- y

    }

    # Change the standard deviation.
    if(!is.null(sigma)){
      y = y * sigma
    }
    # Change the mean.
  if(!is.null(mu)){
      y = y + mu
  }

  return(y)

}


#' Generate fractional Brownian motion
#'
#' @param H Hurst exponent
#' @param N Length of noise series
#' @param mu Mean
#' @param sigma SD
#'
#' @return fBm
#' @export
#'
noise_fBm <- function(H=1.5, N = 512, mu = NULL, sigma = NULL){

  # Determine whether fGn or fBn should be produced.
  if(H%)[%c(1,2)){stop("H must be in (1,2] for fBm!")}

 y <- noise_fGn(H = H, mu = NULL,sigma = NULL)

 # Change the standard deviation.
 if(!is.null(sigma)){
   y = y * sigma
 }
 # Change the mean.
 if(!is.null(mu)){
   y = y + mu
 }

 return(cumsum(y))
}

noise_multifractal <- function(type = c("fGn","fBn","PSD"), N = 512){

  # numb1=10;
  # numb2=1000;
  # alpha=2;
  # N1=numb1*(numb2+N);
  #
  # mGnSum=zeros(numb1,1);
  # mGnSumm=zeros(numb1*(numb2-1),1);
  # mGn=zeros(N,1);
  #
  # R=randn(N1,1);
  # for t=1:N,
  # for n=1:numb1;
  # mGnSum(n)=(n^(Ht(t)-(1/alpha)))*R(1+(numb1*(numb2+t))-n);
  # end;
  # mGnSum1=sum(mGnSum);
  # for nn=1:(numb1*(numb2-1));
  # mGnSumm(nn)=(((numb1+nn)^(Ht(t)-(1/alpha)))-nn^(Ht(t)-(1/alpha)))*R(1+(numb1*(numb2-1+t))-nn);
  # end;
  # mGnSum2=sum(mGnSumm);
  # mGn(t)=((numb1^(-Ht(t)))/gamma(Ht(t)-(1/alpha)+1))*(mGnSum1+mGnSum2);
  # end;
  #
  # mBm=cumsum(mGn);

  # switch(type,
  #        fGn =
  #
  #        )


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
  return(y)
}


#' Create Cauchy Flight
#'
#' Creates a Cauchy flight by taking increments from the Cauchy distributions
#' implemented as the stable distribution ([stabledist::rstable()]) with index paramter `alpha = 1`
#' and skewness parameter `beta = 0`.
#'
#' @param N Length of time series (default = `1000`)
#' @param ndims Number of dimensions (default = `2`)
#' @param alpha Index of stability parameter in `(0,2]`
#' @param beta  Skewness parameter in `[-1,1]`
#' @param scale Scale parameterin `(0,Inf)`
#' @param location Location (shift) parameter in `[-Inf,Inf]`
#'
#' @return A data frame with `ndims` columns and `N` rows.
#'
#' @export
#'
#' @examples
#'
#' df <- flight_Cauchy()
#' plot(density(diff(df$dim1)))
#' plot(df$dim1, df$dim2, type = "l")
#'
flight_Cauchy <- function(N=1000, ndims = 2, alpha = 1, beta = 0, scale = 1, location = 0){
  if(!all(alpha%(]%c(0,2),beta%[]%c(-1,1),scale%()%c(0,Inf),location%()%c(-Inf,Inf))){
    stop(paste("\nValid parameter values:\n alpha in (0,2]\n beta in [-1,1]\n scale in (0,Inf)\n location in (-Inf,Inf)"))}
  if(alpha>1){warning("Not a Cauchy distribution if alpha > 1")}

  df <- data.frame(matrix(rep(NA,N*ndims),ncol = ndims, dimnames = list(NULL,paste0("dim",1:ndims))))
  for(c in 1:NCOL(df)){df[,c] <- cumsum(stabledist::rstable(n=N, alpha = alpha, beta = beta, gamma = scale, delta = location, pm = 2))}
  return(df)
}

#' Create Rayleigh Flight (Brownian Motion)
#'
#' Creates a Rayleigh flight by taking increments from the Normal distributions
#' implemented as the stable distribution ([stabledist::rstable()]) with index paramter `alpha = 2`
#' and skewness parameter `beta = 0`.
#'
#' @inheritParams flight_Cauchy
#'
#' @return A data frame with `ndims` columns and `N` rows.
#'
#' @export
#'
#' @examples
#'
#' df <- flight_Rayleigh()
#' plot(density(diff(df$dim1)))
#' plot(df$dim1, df$dim2, type = "l")
#'
flight_Rayleigh <- function(N=1000, ndims = 2, alpha = 2, beta = 0,scale = 1, location = 0){
  if(!all(alpha%(]%c(0,2),beta%[]%c(-1,1),scale%()%c(0,Inf),location%()%c(-Inf,Inf))){
    stop(paste("\nValid parameter values:\n alpha in (0,2]\n beta  in [-1,1]\n scale in (0,Inf)\n location in (-Inf,Inf)"))}
  if(alpha!=2){warning("Not a Rayleigh (Normal) distribution if alpha != 2")}

  df <- data.frame(matrix(rep(NA,N*ndims),ncol = ndims, dimnames = list(NULL,paste0("dim",1:ndims))))
  for(c in 1:NCOL(df)){df[,c] <- cumsum(stabledist::rstable(n=N, alpha = alpha, beta = beta, gamma = scale, delta = location, pm = 2))}
  return(df)
}

#' Create a Levy-Pareto flight
#'
#' Creates a Rayleigh flight by taking increments from the Normal distributions
#' implemented as the stable distribution ([stabledist::rstable()]) with index paramter `alpha = 1.5`
#' and skewness parameter `beta = 0`.
#'
#' Note that the increments are not strictly from the distribution called **the** Levy distribution, but rather **a**
#' a Levy-with-Pareto-tail-type distribution (i.e. when `1 < alpha < 2`). Use `alpha = 1/2` and `beta = 1` if **the** Levy distribution
#' is required.
#'
#' @inheritParams flight_Cauchy
#'
#' @return A data frame with `ndims` columns and `N` rows.
#'
#' @export
#'
#' @examples
#'
#' # Levy-Pareto
#' df <- flight_LevyPareto()
#' plot(density(diff(df$dim1)))
#' plot(df$dim1, df$dim2, type = "l")
#'
#' # "The" Levy distribution
#' df <- flight_LevyPareto(alpha = 1/2, beta = 1)
#' plot(density(diff(df$dim1)))
#' plot(df$dim1, df$dim2, type = "l")
#'
flight_LevyPareto <- function(N=1000, ndims = 2, alpha = 1.5, beta = 0, scale = 1, location = 0){
  if(!all(alpha%(]%c(0,2),beta%[]%c(-1,1),gamma%()%c(0,Inf),location%()%c(-Inf,Inf))){
    stop(paste("\nValid parameter values:\n alpha in (0,2]\n beta  in [-1,1]\n scale in (0,Inf)\n location in (-Inf,Inf)"))}
  if(alpha%)(%c(1,2)){warning("Not a Levy-Pareto distribution if alpha <= 1 or alpha >= 2")}

  df <- data.frame(matrix(rep(NA,N*ndims),ncol = ndims, dimnames = list(NULL,paste0("dim",1:ndims))))
  for(c in 1:NCOL(df)){df[,c] <- cumsum(stabledist::rstable(n=N, alpha = alpha, beta = beta, gamma = scale, delta = location, pm = 2))}
  return(df)
}


# Convert decimal point
c2p <- function(text,N=1){
  if(!is.character(text)){text<-as.character(text)}
  if(sum(grepl("[,]",text))>=N){text <- gsub(",",".",text)}
  return(text)
}

# Count missing values in x
nmissing <- function(x){
  sum(is.na(x))
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
elascer <- function(x,mn=NA,mx=NA,lo=0,hi=1,groupwise = FALSE, keepNA = TRUE, boundaryPrecision = NA, tol= .Machine$double.eps^0.5){

  doGroupwise <- FALSE
  mnNA <- FALSE
  mxNA <- FALSE
  UNLIST      <- FALSE
  if(length(dim(x))<2){UNLIST <- TRUE}
  if(any(is.na(mn),is.na(mx))){
    if(groupwise&NCOL(x)>1){doGroupwise <- TRUE}
    if(is.na(mn)){
      mnNA <- TRUE
      mn <- min(x,na.rm=TRUE)
      }
    if(is.na(mx)){
      mxNA <- TRUE
      mx <- max(x,na.rm=TRUE)
      }
  }
  x <- as.data.frame(x)
  isNA <- plyr::colwise(is.na)(x)
  u <- x
  for(i in 1:NCOL(x)){
    if(doGroupwise){
      if(mnNA){mn<-min(x[,i],na.rm=TRUE)}
      if(mxNA){mx<-max(x[,i],na.rm=TRUE)}
    }
    if(mn>mx){warning("Minimum (mn) >= maximum (mx).")}
    if(lo>hi){warning("Lowest scale value (lo) >= highest scale value (hi).")}
    if(mn==mx){
      u[,i]<-rep(mx,length(x[,i]))
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


# [copied from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/ ]
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multi_PLOT <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    grid::grid.newpage()
    grid::grid.draw(plots[[1]])
    #print(plots[[1]])

  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      grid::grid.draw(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
    }
  }
}

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


#' TRY ... CATCH
#'
#' @details
#'  In longer simulations, aka computer experiments,
#'  you may want to
#'  1) catch all errors and warnings (and continue)
#'  2) store the error or warning messages
#'
#'  Here's a solution  (see \R-help mailing list, Dec 9, 2010):
#'
#' Catch *and* save both errors and warnings, and in the case of
#' a warning, also keep the computed result.
#'
#' @title tryCatch both warnings (with value) and errors
#' @param expr an \R expression to evaluate
#' @return a list with 'value' and 'warning', where value' may be an error caught.
#' @author Martin Maechler; Copyright (C) 2010-2012  The R Core Team
#' @export
#' @keywords internal
try_CATCH <- function(expr){
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
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

