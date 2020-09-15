# package casnet ----
#
# Dynamic Complexity ----
#
# dc functions


#' Dynamic Complexity
#'
#' Calculates Dynamic Complexity, a complexity index for short and coarse grained time series.
#'
#' @param df A data frame containing multivariate time series data from 1 person. Rows should indicate time, columns should indicate the time series variables. All time series in `df` should be on the same scale, an error will be thrown if the range of the time series in`df` is not `[scale_min,scale_max]`.
#' @param win Size of window in which to calculate Dynamic Complexity. If `win < NROW(df)` the window will move along the time series with a stepsize of `1` (default = `NROW(df)`)
#' @param scale_min The theoretical minimum value of the scale. Used to calculate expected values, so it is important to set this to the correct value.
#' @param scale_max The theoretical maximum value of the scale. Used to calculate expected values, so it is important to set this to rhe correct value.
#' @param doPlot If `TRUE` shows a Complexity Resonance Diagram of the Dynamic Complexity and returns an invisible [ggplot2::ggplot()] object. (default = `FALSE`)
#' @param doPlotF If `TRUE` shows a Complexity Resonance Diagram of the Fluctuation Intensity and returns an invisible [ggplot2::ggplot()] object. (default = `FALSE`)
#' @param doPlotD If `TRUE` shows a Complexity Resonance Diagram of the Distribution Uniformity and returns an invisible [ggplot2::ggplot()] object. (default = `FALSE`)
#' @param returnFandD Returns a list object containing the dynamic complexity series as well as the `F` and `D` series. (default = `FALSE`)
#' @param useVarNames Use the column names of `df` as variable names in the Complexity Resonance Diagram (default = `TRUE`)
#' @param colOrder If `TRUE`, the order of the columns in `df` determines the of variables on the y-axis. Use `FALSE` for alphabetic/numeric order. Use `NA` to sort by by mean value of Dynamic Complexity (default = `TRUE`)
#' @param useTimeVector Parameter used for plotting. A vector of length `NROW(df)`, containing date/time information (default = `NA`)
#' @param timeStamp If `useTimeVector` is not `NA`, a character string that can be passed to [lubridate::stamp()] to format the the dates/times passed in `useTimeVector` (default = `"01-01-1999"`)
#'
#' @return If `doPlot = TRUE`, a list object containing a data frame of Dynamic Complexity values and a `ggplot2` object of the dynamic complexity resonance diagram. If `doPlot = FALSE` the data frame with Dynamic Complexity series is returned.
#'
#' @export
#'
#' @author Merlijn Olthof
#' @author Fred Hasselman
#'
#' @family Dynamic Complexity functions
#'
#'
dc_win <- function(df, win=NROW(df), scale_min, scale_max, doPlot = FALSE, doPlotF = FALSE, doPlotD = FALSE, returnFandD = FALSE, useVarNames = TRUE, colOrder = TRUE, useTimeVector = NA, timeStamp = "01-01-1999"){

  if(any(stats::is.ts(df),xts::is.xts(df),zoo::is.zoo(df))){
    time_vec <- stats::time(df)
  } else {
    time_vec <- NA
  }

  if(is.null(dim(df))){
    if(is.numeric(df)){
      df <- data.frame(df)
    } else {
      df <- data.frame(as.numeric_discrete(df))
    }
  }


  if(any(df<scale_min,df>scale_max)){
    stop("Range of values in df is outside [scale_min,scale_max]!")
  }

  if(win<=0){stop("Need a window > 0")}

  tsNames <- colnames(df)

  data_f <- dc_f(df = df, win = win, scale_min = scale_min, scale_max = scale_max, doPlot = doPlotF, useVarNames = useVarNames, colOrder = colOrder, useTimeVector = useTimeVector, timeStamp = timeStamp)
  data_d <- dc_d(df = df, win = win, scale_min = scale_min, scale_max = scale_max, doPlot = doPlotD, useVarNames = useVarNames, colOrder = colOrder, useTimeVector = useTimeVector, timeStamp = timeStamp)

  if(doPlotF){
    graph_f <- data_f$data
    data_f  <- data_f$data
    #graphics::plot(graph_f)
  }

  if(doPlotD){
    graph_d <- data_d$data
    data_d  <- data_d$data
    # graphics::plot(graph_d)
  }

  df_win <- data_f*data_d

  # if(logDC){
  #   for(c in 1:NCOL(df_win)){
  #     idNA <- is.na(df_win[,c])
  #     id0  <- df_win[,c]<=0
  #     if(sum(!idNA&id0,na.rm = TRUE)>0){
  #    df_win[!idNA&id0,c] <- df_win[!idNA&id0,c] + .Machine$double.eps
  #    }
  #   }
  #   df_win <- log(df_win)
  # }


  colnames(df_win) <- tsNames
  attr(df_win, "time")      <- attr(data_f,"time")
  attr(df_win, "scale_min") <- scale_min
  attr(df_win, "scale_max") <- scale_max
  attr(df_win, "win")       <- win
  attr(df_win, "dataType")  <- "dc_win"

  g <- NULL
  if(doPlot){
    g <- plotDC_res(df_win = df_win, win = win, useVarNames = useVarNames, colOrder = colOrder, timeStamp = timeStamp, doPlot = doPlot)
  }

  if(returnFandD){
    return(list(dynamic_complexity =  df_win, F_data = data_f, D_data = data_d, plot = g))
  } else {
    return(df_win)
  }

}

#' Cumulative Complexity Peaks (CCP)
#'
#' Computes significant peaks in the dynamic complexity time series.
#'
#' @param df_win A data frame containing series of Dynamic Complexity values obtained by running function [dc_win()]
#' @param alpha_item The significance level of the one-sided Z-test used to determine which peaks are `> 0`.
#' @param alpha_time The significance level of the one-sided Z-test used to determine if the number of significant peaks (as determined  by `alpha_item`) at a specific time stamp are `> 0`.
#' @inheritParams dc_win
#'
#' @return A list with a dataframe of binary complexity peak indices and a cumulative complexity peak index, a CCP diagram.
#'
#' @export
#'
#' @author Merlijn Olthof
#' @author Fred Hasselman
#'
#'
#' @family Dynamic Complexity functions
#'
#'
dc_ccp = function(df_win , alpha_item = 0.05, alpha_time = 0.05, doPlot = FALSE, useVarNames = TRUE, colOrder = TRUE, useTimeVector = NA, timeStamp = "01-01-1999"){

  if(attr(df_win,"dataType")!="dc_win"){
    stop("Argument df_win must be generated by function dc_win()!")
  } else {
    attr(df_win, "time")      -> timeStamp
    attr(df_win, "from_scale_min") -> from_scale_min
    attr(df_win, "from_scale_max") -> from_scale_max
    attr(df_win, "to_scale_min") -> to_scale_min
    attr(df_win, "to_scale_max") -> to_scale_max
    attr(df_win, "col_first") -> col_first
    attr(df_win, "col_last")  -> col_last
    attr(df_win, "win")       -> win
  }

  Zitem <- stats::qnorm(1-alpha_item)
  Ztime <- stats::qnorm(1-alpha_time)


  df_ccp <- data.frame(matrix(data=NA, nrow=nrow(df_win), ncol=ncol(df_win)))
  colnames(df_ccp) <- colnames(df_win)

  for (c in (1:NCOL(df_win))){
    item_z <- ts_standardise(df_win[,c])
    df_ccp[,c] <- as.numeric(item_z > Zitem)
  }

  numPeaks <- rowSums(df_ccp,na.rm = TRUE)

  df_ccp$sig.peaks <- as.numeric(ts_standardise(rowSums(df_ccp,na.rm = TRUE)) > Ztime) * 10

  if(doPlot){

    g <- plotDC_ccp(df_ccp = df_ccp, win = win, useVarNames = useVarNames, colOrder = colOrder, useTimeVector = useTimeVector, timeStamp = timeStamp, doPlot = doPlot)

    return(invisible(list(df   = df_ccp,
                          plot = g)))
  } else {
    return(df_ccp)
  }
}


#' Fluctuation Intensity
#'
#' Fluctuation intensity is one of two components of which the product is the Dynamic Complexity measure.
#'
#' @inheritParams dc_win
#'
#' @return dataframe
#'
#' @export
#'
#' @seealso Use [dc_win()] to get the dynamic complexity measure.
#'
#' @family Dynamic Complexity functions
#'
dc_f <- function(df, win=NROW(df), scale_min, scale_max, doPlot = FALSE, useVarNames = TRUE, colOrder = TRUE, useTimeVector = NA, timeStamp = "01-01-1999"){

  if(any(stats::is.ts(df),xts::is.xts(df),zoo::is.zoo(df))){
    time_vec <- stats::time(df)
  } else {
    time_vec <- NA
  }

  if(is.null(dim(df))){
    if(is.numeric(df)){
      df <- data.frame(df)
    } else {
      df <- data.frame(as.numeric_discrete(df))
    }
  }

  if(any(df<scale_min,df>scale_max)){
    stop("Range of values in df is outside [scale_min,scale_max]!")
  }


  if(win<=0){stop("Need a window > 0")}

  tsNames <- colnames(df)

  ew_data_F <- matrix(NA, nrow=NROW(df), ncol=NCOL(df))
  ew_data_F <- data.frame(ew_data_F)
  data <- rbind(df, matrix(0,nrow=2,ncol=NCOL(df), dimnames = list(NULL,colnames(df))))

  s <- scale_max-scale_min
  length_ts <- nrow(data)

  for (column in 1:NCOL(data)){

    distance    <- 1
    fluctuation <- NA
    tsy         <- as.numeric(data[,column])

    # start sliding window
    for(i in (1:(nrow(data)-win-1))){

      y <- NA
      fluct <- NA
      dist_next <- 1
      k <- NA

      # yd2 <- diff(diff(tsy))
      # k <- rep(0,length(yd2))
      # k[yd2!=0] <- 1

      for (j in (0:(win-2))){

        # Start 2nd sliding window
        # k[j+1] <- eval_f(ts[i+j],ts[i+j+1],ts[i+j+2])

        if((tsy[i+j+1] >= tsy[i+j]) & (tsy[i+j+1] > tsy[i+j+2])){
          (k[j+1]<- 1)
        }  else if((tsy[i+j+1] <= tsy[i+j]) & (tsy[i+j+1]<tsy[i+j+2])){
          (k[j+1]<- 1)
        }  else if ((tsy[i+j+1]>tsy[i+j]) & (tsy[i+j+1] == tsy[i+j+2])){
          (k[j+1]<- 1)
        }  else if ((tsy[i+j+1]<tsy[i+j]) & (tsy[i+j+1] == tsy[i+j+2])){
          (k[j+1] <-1)
        }  else if ((tsy[i+j+1]==tsy[i+j]) & (tsy[i+j+1] > tsy[i+j+2])){
          (k[j+1] <-1)
        }  else if ((tsy[i+j+1]==tsy[i+j]) & (tsy[i+j+1] < tsy[i+j+2])){
          (k[j+1] <-1)
        }  else {
          (k[j+1] <- 0)}
      }

      k[win-1] <- 1
      k <- k[1:(win-1)]

      for (g in (1:length(k))){
        if(k[g]==1){
          y[g] <- abs(tsy[i+g]-tsy[i+g-dist_next])
          fluct[g] = (y[g]/((i+g)-(i+g-dist_next)))
          dist_next <- distance
        } else if(k[g]==0){
          y[g]=0
          fluct[g]=0
          dist_next <- dist_next+1}
      }
      ew_data_F[(i+win-1),(column)]<-sum(fluct/(s*(win-1)), na.rm=TRUE)
    }
  }
  ew_data_F <- ew_data_F[1:nrow(df),]
  attr(ew_data_F,"time") <- time_vec
  if(is.null(dim(ew_data_F))){
    ew_data_F <- data.frame(ew_data_F)
  }

  colnames(ew_data_F) <- tsNames
  g <- NULL
  if(doPlot){
    g <- plotDC_res(df_win = ew_data_F, win = win, useVarNames = useVarNames, colOrder = colOrder, timeStamp = timeStamp, doPlot = doPlot, title = "Fluctuation Intensity Diagram")
    return(list(data = ew_data_F, graph = g))
  } else {
    return(ew_data_F)
  }
}


#' Get fluctuations (internal)
#'
#' the fluctuations
#'
#' @param y_win a time series
#' @param s scale range
#'
#'
#' @return fluctuation intensity
#' @export
#'
#' @keywords internal
#'
get_fluct <- function(y_win,s){
  # Get k-points
  signs <- sign(diff(y_win))
  signs <- c(signs,signs[length(signs)])
  dur   <- ts_duration(signs)
  return(sum(abs(y_win[dur$t.end]-y_win[dur$t.start])/dur$duration.time, na.rm = TRUE)/(s*(win-1)))
}


#' @title Distribution Uniformity
#'
#' Distribution Uniformity is one of two components of which the product is the Dynamic Complexity measure.
#'
#' @inheritParams dc_win
#'
#' @return a dataframe
#' @export
#'
#' @seealso Use [dc_win()] to get the Dynamic Complexity measure.
#'
#' @family Dynamic Complexity functions
#'
dc_d <- function (df, win=NROW(df), scale_min, scale_max, doPlot = FALSE, useVarNames = TRUE, colOrder = TRUE, useTimeVector = NA, timeStamp = "01-01-1999"){

  if(any(stats::is.ts(df),xts::is.xts(df),zoo::is.zoo(df))){
    time_vec <- stats::time(df)
  } else {
    time_vec <- NA
  }

  if(is.null(dim(df))){
    if(is.numeric(df)){
      df <- data.frame(df)
    } else {
      df <- data.frame(as.numeric_discrete(df))
    }
  }

  if(any(df<scale_min,df>scale_max)){
    stop("Range of values in df is outside [scale_min,scale_max]!")
  }

  if(win<=0){stop("Need a window > 0")}

  tsNames <- colnames(df)

  ew_data_D <- matrix(NA, nrow=(nrow(df)+2), ncol=NCOL(df))
  ew_data_D <- data.frame(ew_data_D)
  for (column in 1:NCOL(df)){
    tsy <- as.numeric(df[,column])
    dens <- NA
    x <- NA
    y <- NA
    s <- scale_max-scale_min
    y <- pracma::linspace(scale_min, scale_max, win)
    for (i in (1:(length(tsy)-(win-1)))){
      x <- tsy[i:(i+win-1)]
      x <- sort(x)
      r=0
      g=0
      for (e in (1:(win-1))){
        for (d in ((e+1):win)){
          for (a in (e:(d-1))){
            for (b in (((a+1):d))){
              h <- Heaviside((y[b]-y[a])-(x[b]-x[a]))
              if(h==1){
                r <- r + (((y[b]-y[a]) - (x[b]-x[a])))
                g <- g + (y[b]-y[a])
              } else {
                r <- r
                g <- g + (y[b]-y[a])
              }
            }
          }
        }
      }
      ew_data_D[(i+win-1),(column)] <- 1-(r/g)
    }
  }
  ew_data_D <- ew_data_D[(1:nrow(df)),]
  attr(ew_data_D,"time") <- time_vec
  if(is.null(dim(ew_data_D))){
    ew_data_D <- data.frame(ew_data_D)
  }

  colnames(ew_data_D) <- tsNames
  g <- NULL
  if(doPlot){
    g <- plotDC_res(df_win = ew_data_D, win = win, useVarNames = useVarNames, colOrder = colOrder, timeStamp = timeStamp, doPlot = doPlot, title = "Distribution Uniformity Diagram")
    return(list(data = ew_data_D, graph = g))
  } else {
    return(ew_data_D)
  }
}

