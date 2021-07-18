# package casnet ----
#
# Time Series HELPERS
#
# ts functions


#' Permutation Test: Block Randomisation
#'
#' Use block randomistion to get a permutation test evaluation of the deviation of an observed value at each time point from a target value. To do block permutation without any tests, pass `NULL` for argument `targetValue`.
#'
#' @param y1 Time series 1. The goal of the permutation test will be to decide whether the difference `y1-targetValue != 0` for each time point, given `alpha`.
#' @param y2 An optional second time series. If this timeseries is provided then the goal of the permutation test will be the to decide wether the difference `y2-y1 != targetValue` for each time point, given `alpha`.
#' @param targetValue The target value for the permutation test. If `NULL`, the function will return a data frame with the block randomised surrogates  columns (default = `0`)
#' @param Nperms Number of permutations (default = `19`)
#' @param sim Value passed to the `sim` argument of [boot::tsboot()] valid options are: `"model","fixed","geom","scramble"` (default = `"geom"`)
#' @param l Block sizes to use, see [boot::tsboot()] for details (default = `3`)
#' @param alpha Alpha level for deciding significance  (default = `0.05`)
#' @param returnBootObject Return the `boot` object (default = `FALSE`)
#' @param ... Other arguments passed to function [boot::tsboot()]
#'
#' @return A data frame with the difference time series and variables indicating N and significance.
#'
#' @family Time series operations
#'
#' @export
#'
#' @examples
#'
#' \donttest{set.seed(4321)
#' y1 <- rnorm(5000)
#' y2 <- y1-(mean(y1)+rnorm(1))
#'
#' ts_permtest_block(y1 = y1, y2 = y2)}
#'
ts_permtest_block <- function(y1, y2 = NULL, targetValue = 0, Nperms = 19, sim = "geom", l = 3, alpha = .05, returnBootObject = FALSE,...){

  checkPkg("boot")

  tsSame <- function(y){
    y
  }

  #  dotArgs <- list(...)
  #  Args <- methods::formalArgs(boot::tsboot)
  #   Args <- Args[!Args%in%c("tseries","statistic","R","l","sim")]
  #   nameOK  <- names(dotArgs)%in%Args
  #   if(!all(nameOK)){
  #     dotArgs    <- formals(boot::tsboot)
  #     nameOk <- rep(TRUE,length(dotArgs))
  #   }
  #
  #   dotArgs$RM <- dmat
  #   do.call(rp_plot, dotArgs[nameOk])
  # }

  #dat <- ts(data.frame(y1 = y1, y2=y2))

  returnOnlyPerms <- FALSE
  if(is.null(targetValue)){
    targetValue <- 0
    returnOnlyPerms <- TRUE
  }

  ts.boot1 <- boot::tsboot(tseries = stats::ts(y1), statistic = tsSame, R = Nperms, sim = sim, l = l)
  ts.boot <- ts.boot1

  if(!is.null(y2)){
    if(NROW(y2)!=NROW(y1)){stop("y1 and y2 have different lengths.")}
    ts.boot2 <- boot::tsboot(tseries = stats::ts(y2), statistic = tsSame, R = Nperms, sim = sim, l = l)
    ts.boot$t0 <- y1-y2
    ts.boot$t <- ts.boot1$t-ts.boot2$t
  } else {
    ts.boot$t0 <- y1-targetValue
    ts.boot$t  <- ts.boot1$t-targetValue
    targetValue <- 0
  }

  if(!returnOnlyPerms){

    out <- list()
    for(t in 1:NROW(y1)){
      if(ts.boot$t0[t] < targetValue){dir <- "lesser"} else {dir <- "greater"}
      if(dir%in%"lesser"){
        Ds <- sum(ts.boot$t[,t] <= ts.boot$t0[t], na.rm = TRUE) # how many <= to the observed diff <= targetValue?
      } else {
        Ds <- sum(ts.boot$t[,t] >= ts.boot$t0[t], na.rm = TRUE) # how many >= to the observed diff > targetValue?
      }
      SD <- stats::sd(c(ts.boot$t[,t],ts.boot$t0[t]), na.rm = TRUE)
      out[[t]] <- data.frame(time = t,
                             Ori  = ts.boot$t0[t],
                             targetValue = targetValue,
                             Nperms     = Nperms,
                             Ori.Diff   = dir,
                             Ori.Rank   = Ds,
                             Ori.Rank.p = Ds / (Nperms + 1),
                             alpha = alpha,
                             sig = NA,
                             Mean = mean(ts.boot$t[,t],na.rm = TRUE),
                             SD = SD,
                             SE = SD/sqrt(NROW(y1)+1)
      )
    }

    df_sig     <- plyr::ldply(out)
    df_sig$sig <- df_sig$Ori.Rank.p < alpha

    if(returnBootObject){
      return(list(df_sig = df_sig, bootOut = ts.boot))
    } else {
      return(df_sig = df_sig)
    }

  } else {
    if(!is.null(y2)){
      out <- list(y1_RND = t(ts.boot1$t),
                  y2_RND = t(ts.boot2$t))
      names(out) <- c(deparse(substitute(y1)),deparse(substitute(y2)))
      plyr::ldply(out,.id="Source")
    } else {
      out <- list(y1_RND = t(ts.boot1$t))
      names(out) <- deparse(substitute(y1))
      return(plyr::ldply(out,.id="Source"))
    }
  }

}


#' Permutation Test: Transition Matrix
#'
#' Monte Carlo resampling of a time series using a discretised version of `y`, a sequence of `bin` numbers with unique values equal to `nbins`:
#'    1. The discrete version of `y` will be used to generate a transition matrix of size `nbins X nbins`.
#'    2. This transition matrix will be used to resample values
#'
#' @inheritParams ts_permtest_block
#' @param nbins Number of bins to use (default = `ceiling(2*length(y1)^(1/3))`)
#' @param keepNA keepNA
#'
#' @return Resampled series
#' @export
#'
#' @family Time series operations
#'
#' @examples
#'
#' set.seed(4321)
#' y <- rnorm(5000)
#' ts_permtest_transmat(y)
#'
ts_permtest_transmat <- function(y1, y2 = NULL,
                                 targetValue = 0,
                                 nbins = ceiling(2*length(y1)^(1/3)),
                                 Nperms = 19,
                                 alpha = .05,
                                 keepNA = TRUE){
  if(is.null(y2)){
    y <- y1
  } else {
    y <- y2-y1
  }

  bins <- ts_discrete(y, nbins, keepNA = keepNA)

  yn   <- cbind(obin = bins, oy = y)
  tMat <- ts_transmat(yd = bins, nbins = nbins)

  y_mc <- matrix(NA,NROW(tMat),2, dimnames = list(NULL,c("rbin","ry")))
  y_mc[1,] <- yn[sample(NROW(yn),1),]

  for(t in 2:NROW(tMat)){
    y_mc[t,1] <- seq(nbins)[stats::rmultinom(n = 1, size = 1, prob = tMat[y_mc[(t-1),1],])==1]
    pset <- yn[(yn[,1] %in% y_mc[t,1]),2]
    if(length(pset)==0){
      y_mc[t,1] <- NA
    } else {
      y_mc[t,2] <- pset[sample(length(pset),1)]
    }
  }
  return(y_mc)

}


#' Find change indices
#'
#' @param y An indicator variable representing different levels of a variable or factor
#' @param returnRectdata Return a dataframe suitable for shading a `ggplot2` graph with [ggplot2::geom_rect()]
#' @param groupVar Pass a value (length 1) or variable (length of y) that can be used as a variable to join the indices by if `returnRectdata = TRUE`
#' @param labelVar If `y` is not a character vector, provide a vector of labels equal to `length(y)`
#' @param discretize If `y` is a continuous variable, setting `discretize = TRUE` will partition the values of `y` into `nbins` number of bins, each value of `y` will be replaced by its bin number.
#' @param nbins Number of bins to use to change a continuous `y` (if `discretize = TRUE`) into a variable with `nbins` levels
#'
#' @return Either a vector with the indices of change in `y`, or, a data frame with variables `xmin,xmax,ymin,ymax,label`
#'
#' @family Time series operations
#'
#' @export
#'
#' @examples
#'
#'  library(ggplot2)
#'
#'  set.seed(1234)
#'  yy     <- noise_powerlaw(standardise = TRUE, N=50, alpha = -1)
#'  tr     <- ts_levels(yy, doTreePlot = TRUE)
#'  breaks <- ts_changeindex(tr$pred$p, returnRectdata = TRUE)
#'  breaks$cols <- casnet::getColours(length(breaks$label))
#'
#'  ggplot(tr$pred) +
#'    geom_rect(data = breaks,
#'    aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, colour = label, fill = label),
#'    alpha = .3) +
#'    scale_colour_manual(values = breaks$cols) +
#'    scale_fill_manual(values = breaks$cols) +
#'    scale_x_continuous("time", expand = c(0,0)) +
#'    geom_line(aes(x=x,y=y)) +
#'    geom_step(aes(x=x,y=p), colour = "red3", size=1) +
#'    theme_bw() + theme(panel.grid.minor = element_blank())
#'
ts_changeindex <- function(y, returnRectdata=FALSE, groupVar = NULL, labelVar = NULL, discretize=FALSE, nbins = 5){

  if(is.null(names(y))){
    y <- as.numeric_discrete(y)
  }

  if(!is.null(names(y))&is.null(labelVar)){
    labelVar <- names(y)
  }

  if(length(unique(labelVar))>10){
    warning("More than 10 epochs detected!")
  }

  xmax <- c((which(diff(y)!=0))+1, NROW(y))
  names(xmax) <- labelVar[xmax]
  xmin <- c(1,xmax[1:(NROW(xmax)-1)])
  names(xmin) <- labelVar[xmin]
  group <- rep(1,length(xmin))

  if(!is.null(groupVar)){
    if(length(groupVar)==1){
      group <- groupVar
    } else {
      if(length(groupVar)==length(y)){
        group <- groupVar[xmax]
      } else {
        warning("Group values incorrect, using: groupVar = 1")
      }
    }
  }

  if(!is.null(labelVar)){
    label <- labelVar[xmax]
  } else {
    label <- y[xmax]
  }

  if(returnRectdata){
    return(data.frame(group = group, xmin = xmin, xmax = xmax, ymin = -Inf, ymax = +Inf, label = label))
  } else {
    return(unique(c(xmin,xmax)[-c(1,length(c(xmin,xmax)))]))
  }
}


#' Time series to Duration series
#'
#' @param y A time series, numeric vector, or categorical variable.
#' @param timeVec A vector, same length as `y` containing timestamps, or, sample indices.
#' @param fs Optional sampling frequency if timeVec represents sample indices. An extra column `duration.fs` will be added which represents `1/fs * duration in samples`
#' @param tolerance A number `tol` indicating a range `[y-tol,y+tol]` to consider the same value. Useful when `y` is continuous (`default = 0`)
#'
#' @return A data frame
#' @export
#'
#' @family Time series operations
#'
#' @examples
#' library(invctr)
#' # Create data with events and their timecodes
#' coder <- data.frame(beh=c("stare","stare","coffee","type","type","stare"),t=c(0,5,10,15,20,25))
#'
#' ts_duration(y = coder$beh, timeVec = coder$t)
#'
ts_duration <- function(y, timeVec = stats::time(y), fs = stats::frequency(y), tolerance = 0){


  if(plyr::is.discrete(y)){
    y.n <- as.numeric_discrete(y)
  } else {
    y.n <- as.numeric(y)
    names(y.n) <- paste(y.n)
    if(all(is.na(y.n))){
      stop("Conversion to numeric failed for all elements of y.")
    }
  }

  y <- y.n
  tID <- seq_along(y)[-1]
  y[NROW(y)+1]       <- max(y, na.rm = TRUE) + tolerance + 1
  timeVec[NROW(timeVec)+1] <- max(timeVec, na.rm = TRUE)

  same <- list()
  same[[1]] <- data.frame(y = y[1],
                          y.name  = names(y)[1],
                          t.start = timeVec[1],
                          t.end   = timeVec[1],
                          duration.time    = 0,
                          duration.samples = 1,
                          duration.fs      = fs,
                          keep = TRUE)

  for(i in tID){
    # Same as previous?
    if(y[i]%[]%c((y[i-1]-tolerance),(y[i-1]+tolerance))){
      same[[i]] <- data.frame(y = y[i],
                              y.name  = names(y[i]),
                              t.start = same[[i-1]]$t.start,
                              t.end   = timeVec[i],
                              duration.time    = 0,
                              duration.samples = (same[[i-1]]$duration.samples+1),
                              duration.fs      = fs,
                              keep = TRUE)
      same[[i-1]]$keep <- FALSE
    } else {
      same[[i-1]]$t.end <- timeVec[i]
      # Same as upcoming?
      if(y[i]%[]%c((y[i+1]-tolerance),(y[i+1]+tolerance))){
        same[[i]] <- data.frame(y = y[i],
                                y.name  = names(y[i]),
                                t.start = timeVec[i],
                                t.end   = timeVec[i+1],
                                duration.time    = 0,
                                duration.samples = 1,
                                duration.fs = fs,
                                keep = FALSE)
      } else {

        same[[i]] <- data.frame(y = y[i],
                                y.name  = names(y[i]),
                                t.start = timeVec[i],
                                t.end   = timeVec[i+1],
                                duration.time    = 0,
                                duration.samples = 1,
                                duration.fs = fs,
                                keep = TRUE)
      }
    }
  }
  same.out <- plyr::ldply(same)
  same.out <- same.out[same.out$keep,1:6]
  if(!is.null(fs)){same.out$duration.fs = (1/fs)*same.out$duration.samples}
  same.out$duration.time = (same.out$t.end - same.out$t.start)
  row.names(same.out) <- paste(seq_along(same.out$t.start))
  return(same.out)
}



#' Delay embedding of a time series
#'
#' Create a state vector based on an embedding lag and a number of embedding dimanesions.
#'
#' @param y Time series
#' @param emDim Embedding dimension
#' @param emLag Embedding lag
#' @param returnOnlyIndices Return only the index of y for each surrogate dimension, not the values (default = `FALSE`)
#' @param silent Silent-ish mode
#'
#' @note If `emLag = 0`, the assumption is the columns in `y` represent the dimensions and `y` will be returned with attributes `emLag = 0` and `emDim = NCOL(y)`. If `emLag > 0` and `NCOL(y)>1` the first column of `y` will used for embedding and a warning will be triggered.
#'
#' @return The lag embedded time series
#' @family Time series operations
#'
#' @author Fred Hasselman
#'
#' @export
#'
#' @family Time series operations
#'
ts_embed <- function (y, emDim, emLag, returnOnlyIndices = FALSE, silent = TRUE){

  id    <- ifelse(is.null(colnames(y)),ifelse(is.null(names(y)),deparse(substitute(y)),names(y)[1]),colnames(y)[1])
  y.ori <- y
  N     <- NROW(y)

  if(!is.null(dim(y))&emLag>0){
    y <- y_data <- as.numeric(y[,1])
    #if(!silent){cat("\ntaking first column...\n")}
    if(NCOL(y)>1){
      warning(cat("\nMultiple columns available, taking first column of y...\n"))
    }
  } else {
    y_data  <- y
  }

  if(any(stats::is.ts(y), zoo::is.zoo(y), xts::is.xts(y))){
    y <- stats::time(y)
    emTime <- lubridate::as_datetime(y[emLag+1])- lubridate::as_datetime(y[1])
  } else {
    y <- zoo::index(y)
    emTime <- emLag
  }


  if(emLag==0){emDim <- 1}

  #   emY   <- matrix(nrow = N, ncol = NCOL(y.ori), byrow = TRUE, dimnames = list(NULL,colnames(y.ori)))
  #   #emY   <- as.matrix(y.ori)
  #
  # } else {

    if((emDim-1) * emLag > N){
      stop(paste0("Time series length (N = ",N,") is too short to embed in ",emDim," dimensions with delay ",emLag))
    }


    if(emDim > 1){
      lag.id    <- seq(1, (emDim*emLag), emLag)
      maxN      <- (N+1) - max(lag.id)
      emY       <- matrix(nrow = maxN, ncol = emDim)

      for(tau in seq_along(lag.id)){
        emY[,tau] = y[lag.id[tau]:((N+1)-(rev(lag.id)[tau]))]
      }
      colnames(emY) <- paste0("tau.",0:(emDim-1))

    } else {
      ids <- colnames(y.ori)
      if(is.null(ids)){ids <- 0}
      emY <- as.matrix(y.ori)
      dimnames(emY) <- list(NULL,paste0("tau.",ids))
    }

    # Alternative: rollapply(y, list(-d * seq(0, k-1)), c)


  #} # if emLag = 0

  attr(emY, "embedding.dims") <- emDim
  attr(emY, "embedding.lag")  <- emLag
  attr(emY, "embedding.time") <- emTime
  attr(emY, "variable.y")     <- id

  if(returnOnlyIndices){
    attr(emY, "data.y") <- y_data
    return(emY)
  } else {
    if(emDim>1&emLag>0){
      for(c in 1:NCOL(emY)){
        emY[,c] <-  y_data[emY[,c]]
      }
    }
    #else{
    #  emY <- as.matrix(y_data)
    #}
    return(emY)
  }
}


#' Discrete representation
#'
#'  Return a discrete representation of `y` by binning the observed values and returning the transfer probabilities.
#'
#' @param y Numeric vector or time series to be discretised.
#' @param nbins Number of bins to use for calculating transfer probabilities (default = `ceiling(2*length(y)^(1/3))`)
#' @param keepNA If `TRUE`, any `NA` values will first be removed and later re-inserted into the discretised time series.
#'
#' @return A discretised version of `y`
#'
#' @family Time series operations
#'
#' @export
#'
ts_discrete <- function(y, nbins=ceiling(2*NROW(y)^(1/3)), keepNA = TRUE){

  idNA <- is.na(y)
  y <- y[!idNA]

  # Make a binvector
  binvec <- seq(min(y)-0.001, max(y),length.out=nbins+1)
  # Apply the bins
  bins <- vector("numeric",NROW(y))
  for(b in 1:nbins) {
    bins[y%(]%c(binvec[b],binvec[b+1])] <- b
  }

  if(keepNA){
    y <- rep(NA,length(idNA))
    y[!idNA] <- bins
  } else {
    y <- bins
  }

  return(y)
}

#' Symbolic representation
#'
#'  Return a discrete representation of `y` by binning the observed values and returning the transfer probabilities.
#'
#' @param y Numeric vector or time series to be discretised.
#' @param keepNA If `TRUE`, any `NA` values will first be removed and later re-inserted into the discretised time series.
#' @param usePlateaus Treat consequative "same" values after "peak" or "trough" as a "peak"/"trough".
#' @param doPlot Create a plot of the symbolized series.
#'
#' @return A discretised version of `y`
#'
#' @family Time series operations
#'
#' @export
#'
ts_symbolic <- function(y, keepNA = TRUE, usePlateaus = FALSE, doPlot = FALSE){

  cl <- class(y)
  cnames <- colnames(y)
  if(is.null(cnames)){
    cnames = paste0("Y",1:NCOL(y))
  }

  idNA <- plyr::colwise(.fun = is.na)(as.data.frame(y))
  #y <- as.data.frame(y)[stats::complete.cases(!idNA),]
  y <- as.data.frame(y)[stats::complete.cases(y),]

  # returns matrix
  sym_num <- symbolize(xy = y)

  if(is.null(dim(y))) {
    ymat <-  as.matrix(y)
  } else {
    ymat <- y
  }

  # Adjust last value
  for(i in 1:NCOL(sym_num)) {
    t = NROW(sym_num)
    # if(!is.na(ymat[t-1,i])){
    if(ymat[t-1,i]>ymat[t,i]){
      sym_num[t,i] <- 2
    }
    if(ymat[t-1,i]<ymat[t,i]){
      sym_num[t,i] <- 4
    }
    if(ymat[t-1,i]==ymat[t,i]){
      sym_num[t,i] <- 3
    }
    #}
  }

  if(keepNA&any(colSums(idNA,na.rm = TRUE)>0)){
    sym_NA <- matrix(NA,nrow=NROW(idNA),ncol = NCOL(idNA))
    for(c in 1:NCOL(y)){
      sym_NA[!idNA[,c],c] <- sym_num[,c]
    }
    sym_num <- sym_NA
    rm(sym_NA)
  }

  sym_label <- sym_num

  if(!is.null(ncol(sym_num))){
    for(c in 1:NCOL(sym_num)){
      sym_label[,c] <-  dplyr::case_when(
        sym_num[,c]==1 ~ "trough",
        sym_num[,c]==2 ~ "decrease",
        sym_num[,c]==3 ~ "same",
        sym_num[,c]==4 ~ "increase",
        sym_num[,c]==5 ~ "peak",
        is.na(sym_num[,c]) ~ NA_character_)
      if(usePlateaus){
        tmp <- symHelper(sym_label[,c],sym_num[,c])
        sym_label[,c] <- tmp[[1]]
        sym_num[,c]   <- tmp[[2]]}
    }
  } else {
    sym_label <-  dplyr::case_when(
      sym_num==1 ~ "trough",
      sym_num==2 ~ "decrease",
      sym_num==3 ~ "same",
      sym_num==4 ~ "increase",
      sym_num==5 ~ "peak",
      is.na(sym_num) ~ NA_character_)
    if(usePlateaus){
      tmp <- symHelper(sym_label,sym_num)
      sym_label <- tmp[[1]]
      sym_num   <- tmp[[2]]}
  }

  # id <- gregexpr("((increase){1}\\s(same)+\\s(trough){1})+", yc)
  #
  # yc<-paste0(y,collapse=" ")
  # substr(yc,9,9+18-1)

  if(cl%in%c("matrix","numeric")){
    out <- sym_num
    colnames(out)  <- paste0(cnames,"_sym_num")
    anames <- paste0(cnames,"_sym_label")
    for(c in 1:NCOL(out)){
      attr(out,anames[c]) <- factor(sym_label[,c],exclude=NULL)
    }
    if(cl%in%"numeric"){dim(out) <- NULL}
  } else {
    if(!is.null(ncol(y))){
      for(c in 1:NCOL(sym_label)){
        sym_label[,c]  <- factor(sym_label[,c],exclude=NULL)
      }
      colnames(sym_num) <- paste0(cnames,"_sym_num")
      colnames(sym_label) <- paste0(cnames,"_sym_label")
      out <- cbind(sym_num,sym_label) #,stringsAsFactors = FALSE)
    } else {
      out <- factor(sym_label, exclude=NULL)
      attr(out,"sym_numeric") <- sym_num
    }
  }


  if(doPlot){

    shp <- c("trough" = 25, "decrease" = 21, "same"=23,"increase" = 22, "peak" = 24,"missing" = 8)
    col <- c("trough"="darkkhaki","decrease"="orange","same"="grey60","increase"="steelblue","peak"="olivedrab","missing"="red")
    labs <- c("0" = "missing","1"="trough","2"="decrease","3"="same","4"="increase","5"="peak")

    if(!is.null(ncol(sym_num))){
      out$time <- 1:NROW(out)
      df_p1 <- out %>% dplyr::as_tibble() %>% dplyr::select(dplyr::ends_with("_sym_label"),"time") %>% tidyr::gather(key="lab_var", value="sym_label", -"time")
      #df_p1$sym_label[is.na(df_p1$sym_label)] <- "missing"
      df_p2 <- out %>% dplyr::as_tibble() %>% dplyr::select(dplyr::ends_with("_sym_num")) %>% tidyr::gather(key="num_var", value="sym_num")
      df_p2[is.na(df_p2)] <- 0
      df_plot <- cbind(df_p1,df_p2)
      df_plot$num_var <- gsub("_sym_num","",df_plot$num_var)
    } else {
      df_plot <- data.frame(time = 1:NROW(out), sym_num= attr(out,"sym_numeric"), sym_label = out)
      df_plot$sym_num[is.na(df_plot$sym_num)] <- 0
      #levels(df_plot$sym_label)[is.na(df_plot$sym_num)] <- "missing"
    }


    g <- ggplot(df_plot,aes_(x=~time,y=~sym_num)) +
      geom_line(color="grey80") +
      geom_point(aes_(shape=~sym_label, color=~sym_label, fill=~sym_label),size=2)

    if(!is.null(df_plot$num_var)){
      if(length(unique(df_plot$num_var))>1){
        g <- g + facet_grid(num_var~.)
      }
    }

    g <- g +
      scale_x_continuous("Time") +
      scale_shape_manual("Symbolic value", values = shp) +
      scale_color_manual("Symbolic value", values = col) +
      scale_fill_manual("Symbolic value", values = col) +
      scale_y_continuous("",labels = labs) +
      theme_bw() + theme(panel.grid = element_blank())

    graphics::plot(g)

    return(list(data=out,plot=invisible(g)))
  } else {
    return(out)
  }

}


#' Correct for plateaus in symbolic series
#'
#' @param sym_label Labels
#' @param sym_num Numeric labels
#'
#' @return data frame
#' @export
#'
#' @keywords internal
#'
symHelper <- function(sym_label,sym_num){
  out_num <- sym_num
  out_label <- sym_label
  i <- 1
  same <- 0
  while(i<=length(sym_label)){
    if(sym_label[i]%in%c("increase","decrease")){
      samesame <- TRUE
      r <- i+1
      same <- 0
      while(samesame){
        if(sym_label[r]%in%"same"){
          same <- same+1
          r <- r+1
        } else {
          samesame <- FALSE
        }
      }
      if(same>0){

        # sym_num[,c]==1 ~ "trough",
        # sym_num[,c]==2 ~ "decrease",
        # sym_num[,c]==3 ~ "same",
        # sym_num[,c]==4 ~ "increase",
        # sym_num[,c]==5 ~ "peak",

        if(all(!sym_label[i]%in%c("increase"),sym_label[i+same+1]%in%c("peak","increase"))){
          out_label[i:(i+same)] <- "trough"
          out_num[i:(i+same)] <- 1
        } else {
          if(all(!sym_label[i]%in%c("decrease")&sym_label[i+same+1]%in%c("trough","decrease"))){
            out_label[i:(i+same)] <- "peak"
            out_num[i:(i+same)] <- 5
          }
        }
      }
    }
    if(same>0){
      i <- (i+same)
    } else {
      i <- i+1
    }
  }
  return(list(out_label,out_num))
}


#' Convert numeric vectors to symbolic vectors.
#'
#' `symbolize` converts numeric vectors to symbolic vectors. It is a helper
#'   function for `muti`.
#'
#' @param xy An n x 2 `matrix` or `data.frame` containing the two
#'   vectors of interest.
#'
#' @return An (n-2) x 2 `matrix` of integer symbols that indicate whether
#'   the i-th value, based on the i-1 and i+1 values, is a "trough" (=1),
#'   "decrease" (=2), "same" (=3), "increase" (=4), or "peak" (=5).
#'
#' @author Mark Scheuerell (https://mdscheuerell.github.io/muti/)
#' @export
#'
#' @keywords internal
#'
symbolize <- function(xy) {
  ## of interest and converts them from numeric to symbolic.
  ## check input for errors
  if(is.null(dim(xy))) {
    xy <-  as.matrix(xy)
  }
  ## get ts length
  TT <- dim(xy)[1]
  ## init su matrix
  su <- matrix(NA, TT, NCOL(xy))
  ## convert to symbols
  ## loop over 2 vars
  for(i in 1:NCOL(xy)) {
    for(t in 2:(TT-1)) {
      ## if xy NA, also assign NA to su
      if(any(is.na(xy[(t-1):(t+1),i]))) {
        su[t,i] <- NA }
      ## else get correct symbol
      else {
        if(xy[t,i] == xy[t-1,i] | xy[t,i] == xy[t+1,i]) {
          ## same
          su[t,i] <- 3
        }
        if(xy[t,i] > xy[t-1,i]) {
          ## peak
          if(xy[t,i] > xy[t+1,i]) { su[t,i] <- 5 }
          ## increase
          else { su[t,i] <- 4 }
        }
        if(xy[t,i] < xy[t-1,i]) {
          ## trough
          if(xy[t,i] < xy[t+1,i]) { su[t,i] <- 1 }
          ## decrease
          else { su[t,i] <- 2 }
        }
      } ## end else
    } ## end t loop
  } ## end i loop
  ## return su matrix
  return(su)
}


#' Derivative of time series
#'
#' Iteratively differenced series up to `order`. The same length as the original series is recovered by calculating the mean of two vectors for each iteration: One with a duplicated first value and one with a duplicated last value.
#'
#' @param y A timeseries object or numeric vector or a matrix in which columns are variables and rows are numeric values observed over time.
#' @param order How many times should the difference iteration be applied? (default = `1`)
#' @param addColumns Should the derivative(s) be added to the input vector/matrix as columns? (default = `TRUE`)
#' @param keepDerivatives If `TRUE` and `order > 1`, all derivatives from `1:order` will be returned as a matrix )default = `FALSE`)
#' @param maskEdges Mask the values at the edges of the derivatives by any numeric type that is not `NULL` (default = `NULL`)
#' @param silent Silent-ish mode
#'
#' @return Depending on the setting of `addColumns` and the object type passed as `y`, a vector of equal length as `y` iteratively differenced by `order` times; a matrix with derivatives, or a matrix with original(s) and derivative(s).
#'
#' @note The values at the edges of the derivatives represent endpoint averages and should be excluded from any subsequent analyses. Set argument `maskEdges` to a value of your choice.
#'
#' @export
#'
#' @family Time series operations
#'
#' @examples
#'
#' # Here's an interesting numeric vector
#'y<-c(1,1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,1597,2584,4181,6765,10946,17711,28657,46368)
#'
#' # Return the first order derivative as a vector
#' ts_diff(y=y,addColumns=FALSE)
#'
#' # Return original and derivative as a matrix
#' plot(stats::ts(ts_diff(y=y, addColumns=TRUE)))
#'
#' # Works on multivariate data objects with mixed variable types
#' df <- data.frame(x=letters, y=1:26, z=sin(1:26))
#'
#' # Returns only derivatives of the numeric colunmns
#' ts_diff(y=df,addColumns=FALSE)
#'
#' # Returns original data with derivatives of the numeric columns
#' ts_diff(y=df, order=4, addColumns=TRUE)
#'
#' # Plot logistic S-curve and derivatives 1 to 3
#' S <- stats::plogis(seq(-5,5,.1))
#' plot(stats::ts(ts_diff(S, order=3, keepDerivatives = TRUE)))
#' abline(v=which(seq(-5,5,.1)==0), col= "red3", lwd=3)
#'
#' # Plot again, but with masked edges
#' (maskEdge <- ts_diff(S, order=3, keepDerivatives = TRUE, maskEdges = NA))
#' plot(stats::ts(maskEdge))
#' abline(v=which(seq(-5,5,.1)==0), col= "red3", lwd=3)
#'
ts_diff <- function(y, order= 1, addColumns = TRUE, keepDerivatives = FALSE, maskEdges = NULL, silent = TRUE){

  N <- NCOL(y)
  if(NROW(y)<N){
    warning(paste0("Assuming variables in columns and observations in rows!\nDo you really have ",N," columns with time series of length ",NROW(y),"?"))
  }

  if(!is.null(maskEdges)){
    if(is.numeric(maskEdges)|is.na(maskEdges)){
      if(length(maskEdges)==1){
        maskEdges <- as.numeric(maskEdges)
      } else {
        maskEdges = NULL
      }
    } else {
      maskEdges = NULL
    }
  }

  Nnum <- plyr::laply(y,is.numeric)
  if(length(Nnum)<1){stop("Need at least one numeric vector!")}

  if(N==1){
    Nnum <- TRUE
    yy <- matrix(y)
  } else {
    yy <- as.matrix(y[,(1:N)[Nnum]])
  }

  if(is.null(dimnames(yy)[[2]])){
    dimnames(yy) <- list(dimnames(yy)[[1]], gsub("[[:punct:]]","", paste0(deparse(substitute(y)),1:NCOL(yy))))
  }

  dynames <- dimnames(yy)[[2]]

  # dy <- matrix(y[,(1:N)[Nnum]], dimnames = dimnames(y))
  #
  # if(is.null(dimnames(dy)[[2]])){
  #    dynames <- paste0("y",1:NCOL(dy))
  #   } else {
  #    dynames <- paste(dimnames(dy)[[2]])
  #   }
  #
  # dimnames(dy) <- list(dimnames(dy)[[1]],paste0(dynames,"_d",order))

  if((N-NCOL(yy))>=0){
    if(!silent){cat(paste0("\nCalulating derivative for ",NCOL(yy)," time series.\n"))}
  }

  keepDiffs <- list() # matrix(NA, nrow = NROW(dy),ncol = length(1:order)*NCOL(dy))
  maskWhich <- list()

  dy <- yy
  # Repeat the difference operation
  for(o in 1:order){
    dif <- diff(dy)
    dy  <- cbind((rbind(dif[1,], dif)+rbind(dif,dif[NROW(dif),]))/2)
    if(!is.null(maskEdges)){
      maskWhich[[o]] <- c(1:o,(NROW(dy)-o+1):NROW(dy))
    }
    if(keepDerivatives){
      keepDiffs[[o]] <- dy
    }
  }

  if(keepDerivatives){
    dy <- as.data.frame(keepDiffs)
    colnames(dy) <- unlist(plyr::llply(1:order, function(n) paste0(dynames,"_d",n)))
  } else {
    dy <- as.data.frame(dy)
    colnames(dy) <- paste0(dynames,"_d",order)
  }

  if(!is.null(maskEdges)){
    for(c in seq_along(maskWhich)){
      dy[maskWhich[[c]],c] <- maskEdges
    }
  }

  if(addColumns){
    return(cbind(yy,dy))
  } else {
    if(N==1){
      return(dy[,1])
    } else {
      return(dy)
    }
  }
}



#' Adjust time series by summation order
#'
#' Many fluctuation analyses assume a time series' Hurst exponent is within the range of `0.2 - 1.2`. If this is not the case it is sensible to make adjustments to the time series, as well as the resutling Hurst exponent.
#'
#' @param y A time series of numeric vector
#' @param scaleS The scales to consider for `DFA1`
#' @param polyOrder Order of polynomial for detrending in DFA (default = `1`)
#' @param minData Minimum number of data points in a bin needed to calculate detrended fluctuation
#'
#' @return The input vector, possibly adjusted based on `H` with an attribute `"Hadj"` containing an integer by which a Hurst exponent calculated from the series should be adjusted.
#'
#' @details Following recommendations by <https://www.frontiersin.org/files/Articles/23948/fphys-03-00141-r2/image_m/fphys-03-00141-t001.jpg>{Ihlen (2012)}, a global Hurst exponent is estimated using DFA and `y` is adjusted accordingly:
#' \itemize{
#' \item{`1.2 < H < 1.8` first derivative of y, atribute `Hadj = 1`}
#' \item{`H > 1.8` second derivative of y, atribute `Hadj = 2`}
#' \item{`H < 0.2` y is centered and integrated, atribute `Hadj = -1`}
#' \item{`0.2 <= H <= 1.2 ` y is unaltered, atribute `Hadj = 0`}
#' }
#'
#' @references Ihlen, E. A. F. E. (2012). Introduction to multifractal detrended fluctuation analysis in Matlab. Frontiers in physiology, 3, 141.
#'
#' @family Time series operations
#'
#' @export
#'
ts_sumorder <- function(y, scaleS = NULL, polyOrder = 1, minData = 4){

  if(is.null(scaleS)){
    scaleS <- unique(round(2^(seq(2, floor(log2(NROW(y)/2)), by=((floor(log2(NROW(y)/2))-2)/30)))))
  }

  # Check global H
  TSm      <- as.matrix(cbind(t=1:NROW(y),y=ts_integrate(ts_center(y))))
  Hglobal  <- monoH(TSm = TSm, scaleS = scaleS, polyOrder = polyOrder, returnPLAW = TRUE, returnSegments = TRUE)
  rm(TSm)

  fitRange <- which(lapply(Hglobal$segments,NROW)>=minData)

  lmfit1        <- stats::lm(Hglobal$PLAW$bulk ~ Hglobal$PLAW$size.log2, na.action=stats::na.omit)
  H1  <- lmfit1$coefficients[2]
  lmfit2        <- stats::lm(Hglobal$PLAW$bulk.log2[fitRange] ~ Hglobal$PLAW$size.log2[fitRange], na.action=stats::na.omit)
  H2  <- lmfit2$coefficients[2]


  # Adjust TS by global H
  Hadj <- 0
  if(H2%(]%c(1.2,1.8)){
    y <- ts_diff(y,addColumns = FALSE)
    Hadj=1
  }
  if(H2>1.8){
    y <- ts_diff(ts_diff(y, addColumns = FALSE), addColumns = FALSE)
    Hadj <- 2
  }
  if(H2%[]%c(0.2,1.2)){
    Hadj <- 0
  }
  if(H2 < 0.2){
    y <- ts_integrate(ts_center(y))
    Hadj <- -1
  }

  attr(y,"Hglobal.full") <- H1
  attr(y,"Hglobal.excl") <- H2
  attr(y,"Hadj")         <- Hadj

  return(y)
}

#' Check and/or Fix a vector
#'
#' @param y A time series object or numeric vector
#' @param checkNumericVector is 1D numeric vector?
#' @param checkTimeVector has time vector?
#' @param checkWholeNumbers contains only wholenumbers?
#' @param checkPow2 length is power of 2?
#' @param checkScale checkScale
#' @param checkSummationOrder checkSummationOrder
#' @param checkNonStationarity checkNonStationarity
#' @param checkNonHomogeneity checkNonHomogeneity
#' @param fixNumericVector return a 1D numeric vector (WARNING: Data frames and Matrices with NCOL > 1 wil be converted to long form)
#' @param fixWholeNumbers fixWholeNumber
#' @param fixTimeVector fixTimeVector
#' @param fixPow2 foxPow2
#' @param fixNA fixNA
#' @param fixScale fixScale
#' @param fixSummationOrder fixSummationOrder
#' @param fixNonStationarity fixNonStationarity
#' @param fixNonHomogeneity fixNonHomogeneity
#'
#' @return A 'check' report and/or a 'fixed' vector y.
#'
#' @family Time series operations
#'
#' @export
#'
ts_checkfix <- function(y, checkNumericVector = TRUE, checkWholeNumbers = FALSE, checkTimeVector = FALSE, checkPow2 = FALSE, checkScale = FALSE, checkSummationOrder = FALSE, checkNonStationarity = FALSE, checkNonHomogeneity = FALSE, fixNumericVector = FALSE, fixWholeNumbers = FALSE, fixTimeVector = FALSE, fixPow2 = FALSE, fixNA = TRUE, fixScale = FALSE, fixSummationOrder = FALSE, fixNonStationarity = FALSE, fixNonHomogeneity = FALSE){

  outtext <- list()
  pre <- "\ny"
  post <- "\n"
  i <- 0

  whichClass <- class(y)
  yesNumeric <- FALSE
  yesWholeNumbers <- FALSE
  yesTime    <- FALSE
  yesPow2    <- FALSE
  yesScaled  <- FALSE
  yesNonStationary <- FALSE

  # Numeric
  tmptxt <- ""
  txt <- "NUMERIC"
  if(checkNumericVector){
    if(is.numeric(y)){
      txt <- c(txt,"is numeric")
      yesNumeric <- TRUE
    } else {
      tmptxt <- "is NOT numeric"
    }
  }

  if(fixNumericVector&!yesNumeric){
    if(is.null(dim(y))){
      suppressWarnings(y <- as.numeric(unlist(y)))
    } else {
      if(dim(y)[[2]]==1){suppressWarnings(y <- as.numeric(y[,1]))
      } else {
        if(dim(y)[[2]]>1){suppressWarnings(y <- y %>%  tidyr::gather(key="colname",value="value") %>%  as.numeric(.data$value))
        }
      }
    }
    tmptxt <- paste(tmptxt,"... FIXED by as.numeric(y)")
    yesNumeric <- TRUE
  }
  txt <- c(txt, tmptxt)
  outtext[[i%++%1]] <- paste(pre,txt,post)
  rm(txt,tmptxt)


  # Wholenumber
  tmptxt <- ""
  txt <- "WHOLENUMBER"
  if(all(is.wholenumber(y))){
    txt <- c(txt,"only contains whole numbers")
    yesWholeNumbers <- TRUE
  } else {
    tmptxt <- "does NOT only contain whole numbers"
  }
  if(fixWholeNumbers&!yesWholeNumbers){
    suppressWarnings(y <- round(y))
    tmptxt <- paste(tmptxt,"... FIXED by round(y)")
    yesWholeNumbers <- TRUE
  }
  txt <- c(txt, tmptxt)
  outtext[[i%++%1]] <- paste(pre,txt,post)
  rm(txt,tmptxt)

  # TimeVector
  tmptxt <- ""
  txt <- "TIMEVECTOR"
  if(checkTimeVector){
    if(grepl("(ts|mts|xts|zoo)",whichClass)){
      txt <-c(txt,"has a time vector:",stats::tsp(y))
      yesTime <- TRUE
    } else {
      tmptxt <- "does not have a time vector"
    }
  }

  if(is.list(fixTimeVector&!yesTime)){
    if(all(names(fixTimeVector)%in%names(formals(stats::ts)))){
      fixTimeVector$data <- NULL
      etxt <- paste0("stats::ts(data = y, ",paste0(names(fixTimeVector)," = ",fixTimeVector,collapse = ", "),")")
      y <- eval(parse(text = paste(etxt)))
      tmptxt <- paste(tmptxt,"... FIXED by",etxt)
      yesTime <- TRUE
      rm(etxt)
    } else {
      y <- stats::ts(y)
      tmptxt <- paste(tmptxt,"... FIXED by ts(y)")
      yesTime <- TRUE
    }
  }
  txt <- c(txt, tmptxt)
  outtext[[i%++%1]] <- paste(pre,txt,post)
  rm(txt,tmptxt)

  # Stationarity
  tmptxt <- ""
  txt<-"NONSTATIONARY"
  if(checkNonStationarity){
    tst1 <- suppressWarnings(tseries::kpss.test(y,null = "Trend"))$p.value
    tst2 <- 1-suppressWarnings(tseries::adf.test(y, alternative = "stationary"))$p.value
    tst3 <- 1-suppressWarnings(tseries::pp.test(y, alternative = "stationary"))$p.value
    result <- sum(c(tst1,tst2,tst3)<.05, na.rm = TRUE)
    if(result>1){
      txt <- c(txt,paste0(result,"is NOT stationary (",result,"/3 tests conclude y is stationary)"))
      yesNonStationary <- TRUE
    } else {
      txt <- c(txt,paste0(result,"is stationary (",result,"/3 tests conclude y is not stationary)"))
    }
  }
  if(fixNonStationarity&!yesNonStationary){
    y <- ts_detrend(y)
  }
  txt <- c(txt, tmptxt)
  outtext[[i%++%1]] <- paste(pre,txt,post)
  rm(txt,tmptxt)

  return(y)

}

#' Trim or Fill Vectors
#'
#' Trim the largest vector by cutting it, or filling it with `NA`.
#' Fill the shortest vector with padding.
#'
#' @param x A numeric vector
#' @param y A numeric vector
#' @param action Use `"fill"` to fill the shortest vector with `padding` (default); `"trim.cut"` to trim the longest vector to the length of the shortest; `"trim.NA"` to fill the longest vector with `NA`. This is a shortcut for running `action = "trim.cut"` with `padding=NA`, which can be useful if one wants to match the shortest series, but preserve the original length of largest vector.
#' @param type Should trimming or filling take place at the `"end"` (default), or `"front"` of the vector? The option `"center"` will try to distribute trimming by `NA` or filling by `padding` evenly across the front and end of the vector.
#' @param padding A value to use for padding (default = `0`)
#' @param silent Run silent-ish
#'
#' @return A list with two vectors of equal length.
#'
#' @seealso il_mi
#'
#' @family Time series operations
#'
#' @author Fred Hasselman
#'
#' @export
#'
#' @family Time series operations
#'
#'
ts_trimfill <- function(x,y,action=c("fill","trim.cut","trim.NA")[1],type=c("end","center","front")[1],padding=0,silent=TRUE){

  if(all(is.numeric(x),is.numeric(y))){

    l        <- lengths(list(x=x,y=y))
    oriMin   <- names(l)[which.min(l)]
    oriMax   <- names(l)[which.max(l)]
    ldiff    <- abs(diff(l))

    if(ldiff>0){

      if(any(action%in%c("trim.NA","fill"))){
        if(action=="fill"){
          if(!silent){message(paste("Padded shortest input",names(which.min(l)),"with", ldiff,"times",padding))}
        } else {
          padding <- NA
          if(!silent){message(paste("Trimmed longest input",names(which.max(l)),"with",ldiff,"times",padding))}
        } # if "filll
        if(type=="end"){
          out <- list(list(x,y)[[which.max(l)]], c(list(x,y)[[which.min(l)]],rep(padding,ldiff)))
        }
        if(type=="front"){
          out <- list(list(x,y)[[which.max(l)]], c(rep(padding,ldiff),list(x,y)[[which.min(l)]]))
        }
        if(type=="centered"){
          front<- floor(ldiff/2)
          end  <- ceiling(ldiff/2)
          out  <- list(list(x,y)[[which.max(l)]], c(rep(padding,front),list(x,y)[[which.min(l)]],rep(padding,end)) )
        }
      }

      if(action=="trim.cut"){
        if(!silent){message(paste("Trimmed longest input",names(which.max(l)),"by",ldiff))}
        if(type=="end"){
          out <- list(list(x,y)[[which.max(l)]][1:(NROW(list(x,y)[[which.max(l)]])-ldiff)], list(x,y)[[which.min(l)]])
        }
        if(type=="front"){
          out <- list(list(x,y)[[which.max(l)]][ldiff:NROW(list(x,y)[[which.max(l)]])], list(x,y)[[which.min(l)]])
        }
        if(type=="centered"){
          front<-floor(ldiff/2)
          end  <-ceiling(ldiff/2)
          out <- list(list(x,y)[[which.max(l)]][front:(NROW(list(x,y)[[which.max(l)]])-end)], list(x,y)[[which.min(l)]])
        }
      } # if "trim.cut"

    } else { # if ldiff > 0
      if(!silent){message("Vectors have equal length.")}
      out <- list(x=x,y=y)
    }

    names(out) <- c(oriMax,oriMin)
    return(out[order(c(oriMax,oriMin))])

  } else { # is.numeric
    stop("Please use 2 numeric vectors!")
  }
}

#' Get sliding window indices
#'
#' @param y A time series or numeric vector
#' @param win Size of the window to slide across `y`
#' @param step Size of steps between windows. Can be larger than `win`, but is ignored if `overlap` is not {NA}.
#' @param overlap A value between `[0 .. 1]`. If overlap is not `NA` (default), the value of `step` is ignored and set to `floor(overlap*win)`. This produces indices in which the size of `step` is always smaller than `win`, e.g. for fluctuation analyses that use binning procedures to represent time scales.
#' @param adjustY If not `NA`, or, `FALSE` a list object with fields that match one or more arguments of \link[casnet]{ts_trimfill} (except for `x,y`), e.g. `list(action="trim.NA",type="end",padding=NA,silent=TRUE)`. See `Return value` below for details.
#' @param alignment Whether to right (`"r"`), center (`"c"`), or left (`"l"`) align the window.
#'
#' @return If `adjustY = FALSE`, or, a list object with fields that represent arguments of \link[casnet]{ts_trimfill}, then the (adjusted) vector `y` is returned with an attribute `"windower"`. This is a list object with fields that contain the indices for each window that fits on `y`, given `win`, `step` or `overlap` and the settings of `adjustY`. If `adjustY = NA`, only the list object is returned.
#' @export
#'
#' @family Time series operations
#' @family Tools for windowed analyses
#'
ts_windower <- function(y, win=length(y), step=NA, overlap=NA, adjustY=NA, alignment = c("r","c","l")[1]){

  adjustOK <- FALSE
  if(is.list(adjustY)){
    if(any(names(adjustY)%in%names(formals(ts_trimfill)))){
      adjustY[!names(adjustY)%in%names(formals(ts_trimfill))[-c(1,2)]] <- NULL
      adjustY <- adjustY[!lengths(adjustY)==0]
      if(is.null(adjustY)){
        stop("No valid arguments of ts_trimfill passed to adjustY!")
      } else {
        adjustY<- plyr::llply(adjustY, function(e) if(is.character(e)){shQuote(e)} else {e})
        adjustOK <- TRUE
      }
    }
  }

  if(!is.na(win)&is.na(step)&is.na(overlap)){
    step <- win
    warning("step = NA! Stepsize was set to window size.")
  }

  if(!is.na(overlap)){
    if(overlap%[]%c(0,1)){
      if(!is.na(step)){
        warning("Step will be ignored, because overlap has a value.")
      }
      step <- floor(overlap*win)
    } else {
      stop("Overlap must be in [0,1]")
    }
  }

  wIndices <- list()
  wIndex <- seq(from = 1, to = (NROW(y) - win + step), by = step)

  iDiff    <- (dplyr::last(wIndex)+win-1)-NROW(y)
  if(abs(iDiff)>=(win/2)){
    if(iDiff>0){
      wIndex <- wIndex[1:(length(wIndex)-1)]
    } else {
      extra <- dplyr::last(wIndex)+1
      wIndex[length(wIndex)+1] <- extra
    }
  }
  wIndices <- plyr::llply(wIndex, function(i){i:min(i+win-1,length(y))})

  if(adjustOK){
    if(iDiff<0){
      if(adjustY$action%in%"fill"){
        wIndices[[length(wIndices)]] <- c(wIndices[[length(wIndices)]], seq(max(wIndices[[length(wIndices)]]), max(wIndices[[length(wIndices)]]) + abs(iDiff)))
      }
    }
  }

  wID   <- seq_along(wIndices)
  width <- nchar(paste(max(wID)))
  names(wIndices) <- paste0("window: ", formatC(wID, width = width, format = "d", flag = "0"), " | start: ",wIndex," | stop: ",wIndex+win-1)

  iDiff <- max(wIndices[[length(wIndices)]])-NROW(y)

  if(iDiff!=0|any(lengths(wIndices)!=win)){
    if(iDiff<0|any(lengths(wIndices)>win)){
      wIndices[[length(wIndices)]] <- c(wIndices[[length(wIndices)]],(dplyr::last(wIndices[[length(wIndices)]])+1):NROW(y))
      #warning("Last window was enlarged to include all data")
      warning(paste0("Window ",seq_along(wIndices)[lengths(wIndices)>win])," is larger than window size ",win,".\n")
    }
    if(iDiff>0|any(lengths(wIndices)<win)){
      wIndices[[length(wIndices)]] <- c(dplyr::first(wIndices[[length(wIndices)]]):NROW(y))
      warning(paste0("Window ",seq_along(wIndices)[lengths(wIndices)<win])," is smaller than window size ",win,".\n")
      #warning("Last window was shrunk to fit data")
    }
    names(wIndices)[length(wIndices)] <- paste0("window: ",max(seq_along(wIndices))," | start: ",max(wIndex)," | stop: ",max(wIndices[[length(wIndices)]]))
  }

  switch(alignment,
         l = tIndex <- plyr::laply(wIndices, min),#seq(from = 1,     to = (NROW(y) - win + step), by = step),
         c = tIndex <- plyr::laply(wIndices, mean),#seq(from = floor(win/2), to = (NROW(y) - floor(win/2) + step), by = step),
         r = tIndex <- plyr::laply(wIndices, max) #seq(from = win,   to = NROW(y), by = step)
  )

  attr(x = wIndices, which = "time") <- tIndex

  #names(wIndices) <- paste("stepsize",floor(step*win),"| window",1:length(wIndex))
  return(wIndices)
}



#' Detrend a time series
#'
#' @param y A time series ot numeric vector
#' @param polyOrder order Order of polynomial trend to remove
#'
#' @return Residuals after detrending polynomial of order `order`
#'
#' @export
#'
#' @family Time series operations
#'
#' @author Fred Hasselman
#'
#'
ts_detrend <- function(y, polyOrder=1){
  detR <- stats::lm(y~stats::poly(1:length(y), degree=polyOrder))$residuals
  return(detR)
}


#' Detect levels in time series
#'
#'  Use recursive partitioning function [rpart::rpart()] to perform a 'classification' of relatively stable levels in a timeseries.
#'
#' @param y A time series of numeric vector
#' @param minDataSplit An integer indicating how many datapoints should be in a segment before it will be analysed for presence of a level change (default = `12`)
#' @param minLevelDuration Minimum duration (number of datapoint) of a level (default = `round(minDataSplit/3)`)
#' @param changeSensitivity A number indicating a criterion of change that must occur before declaring a new level. Higher numbers indicate higher levels of change must occur before a new level is considered. For example, if `method = "anova"`, the overall `R^2` after a level is introduced must increase by the value of `changeSensitivity`, see the `cp` parameter in [rpart::rpart.control()].     (default = `0.01`)
#' @param maxLevels Maximum number of levels in one series (default = `30`)
#' @param method The partitioning method to use, see the manual pages of [rpart] for details.
#' @param minChange After the call to [rpart], adjust detected level changes to a minimum absolute change in `y`. If a level change is smaller than `minChange`, the previous level will be continued. Note that this is an iterative process starting at the beginning of the series and 'correcting' towards the end (default  = `sd(y, na.rm = TRUE)`)
#' @param doLevelPlot Should a plot with the original series and the levels be produced? (default = `FALSE`)
#' @param doTreePlot Should a plot of the decision tree be produced. This requires package [partykit](https://cran.r-project.org/web/packages/partykit/index.html) (default = `FALSE`)
#'
#' @return A list object with fields `tree` and `pred`. The latter is a data frame with columns `x` (time), `y` (the variable of interest) and `p` the predicted levels in `y`.
#'
#' @export
#'
#' @family Time series operations
#'
#' @author Fred Hasselman
#'
#' @examples
#'
#' # Levels in white noise?
#'
#' set.seed(4321)
#' y <- rnorm(100)
#' wn <- ts_levels(y)
#' plot(wn$pred$x,wn$pred$y, type = "l")
#' lines(wn$pred$p, col = "red3", lwd = 2)
#'
#' # This is due to the default changeSensitivity of 0.01
#'
#' wn2 <- ts_levels(y,changeSensitivity = .1)
#' lines(wn2$pred$p, col = "steelblue", lwd = 2)
#'
#'
ts_levels <- function(y, minDataSplit = 12, minLevelDuration=round(minDataSplit/3), changeSensitivity = 0.01, maxLevels=30, method=c("anova","poisson","class","exp")[1], minChange = sd(y, na.rm = TRUE), doLevelPlot = FALSE, doTreePlot = FALSE){

  checkPkg("rpart")

  x <- seq_along(y)
  dfs  <- data.frame(x=x, y=y)
  tree <- rpart::rpart(y ~ x,
                       method = method,
                       control = list(
                         minsplit  = minDataSplit,
                         minbucket = minLevelDuration,
                         maxdepth  = maxLevels,
                         cp = changeSensitivity),
                       data=dfs)

  dfs$p <- stats::predict(tree, data.frame(x=x))

  if(doLevelPlot){
   g <- ggplot2::ggplot(dfs) +
      geom_line(aes_(x=~x,y=~y)) +
      geom_step(aes_(x=~x,y=~p), colour = "red3", size=1) +
      theme_bw() + theme(panel.grid.minor = element_blank())
   print(g)
  }

  if(doTreePlot){
    checkPkg("partykit")
    plot(partykit::as.party(tree))
  }

  return(list(tree  = tree,
              pred  = dfs))
}


#' Find Peaks or Wells
#'
#' @param y A time series or numeric vector
#' @param window Window in which to look for peaks or wells
#' @param includeWells Find wells?
#' @param minPeakDist Minimum distance between peaks or wells
#' @param minPeakHeight Minimum height / depth for a peak / well
#'
#' @family Time series operations
#'
#' @author Fred Hasselman
#'
#' @return Index with peak or well coordinates
#' @export
#'
ts_peaks <- function (y,
                      window        = 3,
                      includeWells  = FALSE,
                      minPeakDist   = round(window/2),
                      minPeakHeight = .2*diff(range(y, na.rm = TRUE))){

  fp <- function(y,window){
    shape <- diff(sign(diff(y, na.pad = FALSE)))
    pksl <- sapply(which(shape < 0), FUN = function(i){
      z <- i - window + 1
      z <- ifelse(z > 0, z, 1)
      w <- i + window + 1
      w <- ifelse(w < length(y), w, length(y))
      if(all(y[c(z : i, (i + 2) : w)] <= y[i + 1])){
        return(i + 1)
      } else {
        return(numeric(0))
      }
    })
    return(unlist(pksl))
  }

  pks <- fp(y,window)
  if(includeWells){
    wls <- fp(-1*y,window)
    pks <- sort(c(pks,wls))
  }

  distOK    <- diff(c(NA,pks))>=minPeakDist
  distOK[1] <- TRUE

  heightOK <- rep(FALSE,length(pks))
  for(p in seq_along(pks)){
    if(abs(y[pks[p]] - mean(y[c(max(1,(pks[p]-window)):(pks[p]-1),(pks[p]+1):min((pks[p]+window),length(y)))])) >= minPeakHeight){
      heightOK[p] <- TRUE
    }
  }
  return(pks[heightOK&distOK])
}



#' Center a vector
#'
#' @param numvec A numeric vector
#' @param na.rm Set the `na.rm` field
#' @param type Center on the `"mean"` (default) or the `"median"` of the vector.
#'
#' @return A mean or median centered vector
#' @export
#'
#' @family Time series operations
#'
#' @author Fred Hasselman
#'
ts_center <- function(numvec, na.rm=TRUE, type = c("mean","median")[1]){
  if(!is.numeric(numvec)){
    stop("Vector must be numeric!")
  } else {
    switch(type,
           mean   = return(numvec -   mean(numvec, na.rm=na.rm)),
           median = return(numvec - stats::median(numvec, na.rm=na.rm))
    )
  }
}


#' Standardise a vector
#'
#'
#' @param y A time series or numeric vector
#' @param na.rm Set the `na.rm` field
#' @param keepNAvalues If `na.rm = TRUE` and `keepNAvalues = TRUE`, any `NA` values in `y` will be re-inserted after transformation.
#' @param type Center on the `"mean"` and divide by `sd` (default), or center on `"median"` and divide by `mad`
#' @param adjustN If `TRUE`, apply Bessel's correction (divide by `N-1`) or return the unadjusted `SD` (divide by `N`) (default = `TRUE`)
#'
#' @return A standardised vector
#'
#' @export
#'
#' @family Time series operations
#'
#' @author Fred Hasselman
#'
ts_standardise <- function(y, na.rm = TRUE, keepNAvalues = TRUE, type = c("mean.sd","median.mad")[1], adjustN = TRUE){
  if(!is.numeric(y)){
    stop("Vector must be numeric!")
  } else {
    N <- NROW(y)
    if(adjustN){
      SDtype <- "Bessel"
    } else {
      SDtype <- "unadjusted"
    }
    if(keepNAvalues){
      idNA <- !logical(length = N)
    } else {
      idNA <- !is.na(y)
    }
    switch(type,
           mean.sd    = return(((y - mean(y, na.rm = na.rm)) / ts_sd(y,na.rm = na.rm, type = SDtype))[idNA]),
           median.mad = return(((y - stats::median(y, na.rm=na.rm)) / stats::mad(y, na.rm = na.rm))[idNA])
    )
  }
}

#' Standard Deviation estimates
#'
#' Calculates the population estimate of the standard deviation, or the unadjusted standard deviation.
#'
#' @param y Time series or numeric vector
#' @param na.rm Remove missing values before calculations
#' @param type Apply Bessel's correction (divide by N-1) or return unadjusted SD (divide by N)
#' @param silent Silent-ish mode (default = `TRUE`)
#'
#' @return Standard deviation of `y`
#' @export
#'
#' @family Time series operations
#'
ts_sd <- function(y, na.rm = TRUE, type = c("Bessel","unadjusted")[1], silent=TRUE){

  if(!is.numeric(y)){
    stop("Vector must be numeric!")
  }

  if(na.rm){
    y<-y[stats::complete.cases(y)]
  }

  N          <- NROW(y)
  Y_mean     <- mean(y)
  Y_sd       <- stats::sd(y)
  Bessel     <- sqrt((sum((y - Y_mean)^2)/(N-1)))
  unadjusted <- sqrt((sum((y - Y_mean)^2)/N))

  if(type=="Bessel"){
    if(all.equal(Y_sd,Bessel)){
      if(!silent){
        cat(paste0("\nDifference adjusted-unadjusted SD:",Bessel-unadjusted,"\n"))
      }
    } else {
      cat(paste0("\nWARNING: Computation of SD based on N = ",N," and Mean = ",Y_mean," differs from function sd(y)!!!\n"))
      cat(paste0("Difference Bessel-sd(y):",Bessel-Y_sd,"\n"))
    }
    out <- Bessel
  }

  if(type=="unadjusted"){
    if(all.equal(unadjusted,(Y_sd/sqrt(N/(N-1))))){
      if(!silent){
        cat(paste0("\nDifference Unadjusted-Adjusted:",unadjusted-Bessel,"\n"))
      }
    } else {
      cat(paste0("\nWARNING: Computation of SD based on N = ",N," and Mean = ",Y_mean," differs from function sd(y)!!!\n"))
      cat(paste0("Difference unadjusted-(Y_sd/sqrt(N/(N-1))):",unadjusted-(Y_sd/sqrt(N/(N-1))),"\n"))
    }
    out <- unadjusted
  }

  return(out)
}

#' Slice a Matrix
#'
#' Slices rows of a matrix into a list of matrices representing epochs of length `epochSz`.
#'
#' @param y A matrix with timeseries as columns
#' @param epochSz Epoch size
#'
#' @return A list with epochs
#'
#' @export
#'
#' @family Time series operations
#'
#' @author Fred Hasselman
#'
ts_slice <- function(y , epochSz=4){
  if(!is.matrix(y)){
    yy <- as.matrix(y)
  } else {
    yy <- y
  }

  N <- dim(yy)
  wIndex <- plyr::llply(seq(1,N[1],epochSz),function(i) yy[i:min(i+epochSz-1,N[1]),1:N[2]])
  delID <- which(lengths(wIndex)!=epochSz)%00%NA
  for(del in delID){
    if(!any(is.na(delID))){wIndex <- wIndex[-del]}
  }
  names(wIndex) <- paste("epoch",1:length(wIndex),"| size",lengths(wIndex))
  return(wIndex)
}



#' Create a timeseries profile
#'
#' @param y A 1D timeseries
#'
#' @return The profile
#' @export
#'
#' @family Time series operations
#'
#' @examples
#' y <- runif(1000,-3,3)
#' plot(ts(y))
#' y_i <- ts_integrate(y)
#' plot(ts(y_i))
#'
ts_integrate <-function(y){
  #require(zoo)
  if(!all(is.numeric(y),is.null(dim(y)))){
    warning("Need a 1D numeric vector! Trying to convert...")
    y <- as.numeric(as.vector(y))
  }
  idOK <- !is.na(y)
  t    <- zoo::index(y)
  t[!idOK] <- NA
  yy  <- c(y[idOK][1],y[idOK])
  tt  <- c(t[idOK][1],t[idOK])
  yc  <- y
  yc[t[idOK]] <- cumsum(yy[-length(yy)]*diff(tt))
  # recover initial diff value
  yc <- y[idOK][1] + yc
  return(yc)
}



#' Turn a 1D time series vector into a 2D curve
#'
#' @param y A 1D time series object or numeric vector.
#' @param unitSquare Convert the series to a unit square? (default = `FALSE`)
#' @param toSparse Convert to sparse Matrix (default = `FALSE`)
#' @param resolution Factor by which dimensions will be multiplied (default = `2`)
#'
#' @return A (sparse) matrix representing the time series as a curve in 2D space
#' @export
#'
#' @family Time series operations
#'
#' @examples
#'
#' \donttest{
#' y <- rnorm(100)
#' plot(ts(y))
#'
#' y_img <- ts_rasterize(y)
#' image(y_img,col=c("white","black"))}
#'
#'
ts_rasterize <- function(y, unitSquare = FALSE, toSparse = TRUE, resolution = 2){

  if(!is.vector(y, mode = "numeric")){
    stop('Expecting a numeric vector')
  }

  N <- NROW(y)
  x <- stats::time(y)

  if(unitSquare){
    y <- elascer(as.numeric(y))
    x <- elascer(as.numeric(x))
    r <- raster::raster(ncols = N*resolution, nrows = N*resolution, xmn = 0, xmx = 1, ymn = 0, ymx = 1)
  } else {
    r  <- raster::raster(ncols = N*resolution, nrows = length(unique(y))*resolution, xmn = min(x, na.rm = TRUE), xmx = max(x, na.rm = TRUE), ymn = min(y, na.rm = TRUE), ymx = max(y, na.rm = TRUE))
  }

  xy <- cbind(xc=x,yc=y)
  lns <- raster::spLines(xy)
  rm(x,y,xy)

  spm <- t(raster::as.matrix(raster::rasterize(lns, r, background = 0)))
  rm(lns,r)

  if(toSparse){
    spm <- Matrix::as.matrix(spm, sparse = TRUE)
  }

  return(spm)

}


# Help lme4 get a better convergence
# nlopt <- function(par, fn, lower, upper, control) {
#   # Add to call: control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE
#   .nloptr <<- res <- nloptr(par, fn, lb = lower, ub = upper,
#                             opts = list(algorithm = "NLOPT_LN_BOBYQA", print_level = 1,
#                                         maxeval = 1000, xtol_abs = 1e-6, ftol_abs = 1e-6))
#   list(par = res$solution,
#        fval = res$objective,
#        conv = if (res$status > 0) 0 else res$status,
#        message = res$message
#   )
# }


#' Transition matrix
#'
#' Create a transition matrix from a discrete time series, e.g. to generate Monte Carlo simulations.
#'
#' @param yd A discrete numeric vector or time series, e.g. transformed using [ts_discrete()], or, [ts_symbolic()].
#' @param nbins The number of bins used to transform a continuous time series, or, the number of expected (given `nbins`, or, theoretically possible) values for a discrete series (default = `length(unique(yd))`)
#'
#' @return A transition probability matrix
#' @export
#'
#' @examples
#'
#' set.seed(4321)
#'
#' # Random uniform numbers
#' y  <- runif(10,0,20)
#'
#' # Discrete version
#' yd <- ts_discrete(y, nbins = 10)
#'
#' # Transition probabilities
#' ts_transmat(yd = yd, nbins = 10)
#'
#' # Note: The number of 'observed' bins differs from 'expected' bins
#' table(yd)
#'
#' # Not specifying the expected bins will generate a different matrix!
#' ts_transmat(yd = yd, nbins = length(unique(yd)))
#'
ts_transmat <- function(yd, nbins = length(unique(yd))){

  if(!all(is.wholenumber(yd))){
    stop("yd is not a vector of discrete numbers!")
  }

  transmat <- matrix(0,nbins,nbins)
  for(i in 1:nbins){
    for(j in 1:nbins){
      transmat[i,j] <- sum(yd[-NROW(yd)]==i & yd[-1]==j, na.rm = TRUE)
    }
    if(sum(transmat[i,])>0) {
      transmat[i,] <- transmat[i,]/sum(transmat[i,], na.rm = TRUE)
    } else {
      transmat[i,] <- 1/nbins
    }
  }
  return(transmat)
}


#' Change Profile
#'
#' @param y Time series
#' @param win A window in which
#' @param align Alignment of the window, see [zoo::rollapply()] (default = `"right"`)
#' @param keepNA Remove `NA` or return `y` with `NA` (default = `TRUE`)
#'
#' @return Transformed time series
#'
#' @export
#'
#' @references Hasselman, F. and Bosman, A.M.T. (2020). Studying Complex Adaptive Systems With Internal States: A Recurrence Network Approach to the Analysis of Multivariate Time-Series Data Representing Self-Reports of Human Experience. Frontiers of Applied Mathematics and Statistics, 6:9. doi: 10.3389/fams.2020.00009
#'
ts_cp <- function(y, win, align = "right", keepNA = TRUE){

 idOK     <- !is.na(y)

 x <- y[idOK]
 WINmean <- zoo::rollapply(x, width = win, mean, na.rm=TRUE, fill=NA, align=align)

 y[idOK] <- ts_integrate(x-WINmean[1:length(x)])

 return(y)
 #df_se[idP,paste(cn,"cp")%ci%df_se] <- ts_integrate(y)
}



#' Course grain a time series
#'
#' Generate a course grained version of a time series by summarising values into bins.
#'
#' @param y A numeric vector
#' @param grain The bin size in which to summarise the values (default = `2`)
#' @param summaryFunction How should the data be summarized in the bins?
#' @param retainLength Return only the bin values (`FALSE`), or retain the length of the original series? (default = `TRUE`)
#'
#' @return A coarse grained version of `y`.
#'
#' @export
#'
#' @examples
#'
#' set.seed(1234)
#' y <- rnorm(100)
#' y1 <- ts_coarsegrain(y, grain = 3)
#' y2 <- ts_coarsegrain(y, grain = 3, retainLength = TRUE)
#' y3 <- ts_coarsegrain(y, grain = 3, retainLength = TRUE, summaryFunction = "max")
#'
#' t1 <- seq(1,length(y), by = 3)
#'
#' plot(t1+1, y1, col = "red3", type = "l", ylim = c(-3,3), xlab = "time", ylab = "Y")
#' lines(y, col = "grey70")
#' lines(y2, col = "steelblue")
#' lines(y3, col = "green3")
#' legend(60, -1.3, legend=c("Original", "Mean", "Mean + Retain Length", "Max + Retain Length"),
#' lty = 1, col=c("grey70", "red3", "steelblue","green3"), cex = 0.7)
#'
#'
ts_coarsegrain <- function(y, grain = 2, summaryFunction = c("mean","median","min","max")[1], retainLength = FALSE){

  y <- as.numeric(y)

  if(grain>length(y)){
    stop("Cannot coarse grain beyond series length!")
  }

  if(is.na(grain%00%NA)|grain<1){
    grain <- 1
    u   <- length(y)
  } else {
    u <- ceiling(length(y)/grain)
  }
  yy <- matrix(0, nrow = grain, ncol = u)
  yy[1:length(y)] <- y

  if(NCOL(yy) != NCOL(y)){

    yy <-  switch(summaryFunction,
                  min  =  plyr::colwise(min, na.rm = TRUE)(data.frame(yy)),
                  max  =  plyr::colwise(max, na.rm = TRUE)(data.frame(yy)),
                  mean =  plyr::colwise(mean, na.rm = TRUE)(data.frame(yy)),
                  median = plyr::colwise(median, na.rm = TRUE)(data.frame(yy))
    )

    if(retainLength){
      yy <- rep(yy, each = grain)
      if(NCOL(yy) != NCOL(y)){
        warning("Coarse grained series is not exactly the same length as the original.")
      }
      #yy <- yy[1:length(y)]
    }
  }

  attr(yy,"grain") <- grain

  return(as.numeric(yy))
}
