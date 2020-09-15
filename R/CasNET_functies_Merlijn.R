#' @title Dynamic Complexity
#'
#' @description Calculate DynComp (deprecated, use [dc_win] instead)
#'
#' @param df df multivariate time series data - rows=time 1 person  cols = vars
#' @param col_first col_first
#' @param col_last col_last
#' @param col_time if col time use it otherwise assume ....
#' @param win win
#' @param fromScaleMin dat scale xmin
#' @param fromScaleMax data scale xmax
#' @param toScaleMin theoretical scale min
#' @param toScaleMax theoretical scale min
#' @param doPlot Show a plot
#'
#' @return A list object
#' @name dyn_comp-deprecated
#' @seealso \code{\link{casnet-deprecated}}
#' @keywords internal
#'
#' @author Merlijn Olthof
#'
#' @references Guenther
#'
NULL

#' @rdname casnet-deprecated
#' @section \code{dyn_comp}:
#' For \code{dyn_comp}, use \code{\link{dc_win}}.
#'
#' @export

dyn_comp = function(df, col_first, col_last, win=NROW(df),  fromScaleMin=NULL, fromScaleMax=NULL, toScaleMin=NULL, toScaleMax=NULL, rescale = TRUE, doPlot = FALSE){


if(win>0){ok=TRUE}

  if(rescale){
    df[,col_first:col_last] <- elascer(x = df[,col_first:col_last], mn = fromScaleMin, mx = fromScaleMax, lo = toScaleMin, hi = toScaleMax)

  }


  data_f <- DC_F(df=df, win=win, xmin=fromScaleMin, xmax=fromScaleMax, col_first =col_first, col_last = col_last)

  data_d <- DC_D(df=df, win=win, xmin=fromScaleMin, xmax=fromScaleMax, col_first =col_first, col_last = col_last)


  df.comp <- matrix(NA, nrow=nrow(data_f), ncol=ncol(data_f))
  for (column in (1:ncol(data_f))){
    df.comp[,column] <- data_f[,column]*data_d[,column]
  }

  mat.dc <- data.matrix(df.comp)
  colnames(mat.dc) <- c(1:ncol(mat.dc))

  plot.dc <- ggplot2::ggplot(reshape2::melt(mat.dc), aes_(x=~Var1, y=~Var2, fill=~value)) + ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low='blue', high='red', mid='yellow', midpoint=0.15, na.value='white') +
    ggplot2::theme_bw() +
    ggplot2::xlab('Days') +
    ggplot2::ylab('Items')+
    ggplot2::theme(panel.grid.minor = element_line(colour="NA", size=0.5)) +
    ggplot2::scale_y_continuous(minor_breaks=seq(0.5, ncol(mat.dc), 1), breaks=seq(0.5, ncol(mat.dc)+0.5, 1), label=c(0:(ncol(mat.dc)))) +
    ggplot2::scale_x_continuous(minor_breaks=seq(0.5,nrow(mat.dc),1), breaks=seq(0.5,nrow(mat.dc)+0.5,1), label=c(0:nrow(mat.dc))) +
    ggplot2::ggtitle('Complexity Resonance Diagram') +
    ggplot2::coord_cartesian(ylim=c(0.5,(ncol(mat.dc)+0.5)), xlim=c(win-0.5,nrow(mat.dc)+0.5), expand=F)


  return(list(df.comp = df.comp,
              plot.dc = plot.dc))
}



#' Fluctuation Intensity
#'
#' Used to calculate dynamic complexity
#'
#' @param df df
#' @param win win
#' @param xmin xmin
#' @param xmax xmax
#' @param col_first col_first
#' @param col_last col_last
#'
#' @return dataframe
#' @name DC_F-deprecated
#' @seealso \code{\link{casnet-deprecated}}
#' @keywords internal
NULL

#' @rdname casnet-deprecated
#' @section \code{DC_F}:
#' For \code{DC_F}, use \code{\link{dc_f}}.
#'
#' @export
DC_F = function(df, win, xmin, xmax, col_first, col_last){

  ew_data_F <- matrix(NA, nrow=nrow(df), ncol=(col_last-col_first+1))
  ew_data_F <- data.frame(ew_data_F)
  newrows = data.frame(subset(df[1,]))
  newrows[,] <- 0
  newrows[,1] <- NA
  data <- rbind(df, newrows)
  data <- rbind(data, newrows)

  s <- xmax-xmin
  length_ts <- nrow(data)

  for (column in (col_first:col_last)){
    distance<-1
    fluctuation <- NA
    ts <- data[,column]

    for (i in (1:(nrow(data)-win-1))){
      y <- NA
      fluct <- NA
      dist_next <- 1
      k <- NA

      for (j in (0:(win-2))){
        if((ts[i+j+1] >= ts[i+j]) & (ts[i+j+1] > ts[i+j+2])){
          (k[j+1]<- 1)
        }  else if((ts[i+j+1] <= ts[i+j]) & (ts[i+j+1]<ts[i+j+2])){
          (k[j+1]<- 1)
        }  else if ((ts[i+j+1]>ts[i+j]) & (ts[i+j+1] == ts[i+j+2])){
          (k[j+1]<- 1)
        }  else if ((ts[i+j+1]<ts[i+j]) & (ts[i+j+1] == ts[i+j+2])){
          (k[j+1] <-1)
        }  else if ((ts[i+j+1]==ts[i+j]) & (ts[i+j+1] >ts[i+j+2])){
          (k[j+1] <-1)
        }  else if ((ts[i+j+1]==ts[i+j]) & (ts[i+j+1] <ts[i+j+2])){
          (k[j+1] <-1)
        }  else {
          (k[j+1] <- 0)}
      }
      k[win-1] = 1
      k<-k[1:(win-1)]
      for (g in (1:length(k))){
        if(k[g]==1){
          y[g] <- abs(ts[i+g]-ts[i+g-dist_next])
          fluct[g] = (y[g]/((i+g)-(i+g-dist_next)))
          dist_next <- distance
        } else if(k[g]==0){
          y[g]=0
          fluct[g]=0
          dist_next <- dist_next+1}
      }
      ew_data_F[(i+win-1),(column-col_first+1)]<-sum(fluct/(s*(win-1)), na.rm=TRUE)
    }
  }
  ew_data_F <- ew_data_F[1:nrow(df),]
  return(ew_data_F)
}


#' Distribution Uniformity
#'
#' @param df data
#' @param win win
#' @param xmin xmin
#' @param xmax xmax
#' @param col_first col_first
#' @param col_last col_last
#'
#' @return a dataframe
#' @name DC_D-deprecated
#' @seealso \code{\link{casnet-deprecated}}
#' @keywords internal
NULL

#' @rdname casnet-deprecated
#' @section \code{DC_D}:
#' For \code{DC_D}, use \code{\link{dc_d}}.
#'
#' @export
DC_D = function (df, win, xmin, xmax, col_first, col_last){

  ew_data_D <- matrix(NA, nrow=(nrow(df)+2), ncol=(col_last-col_first+1))
  ew_data_D <- data.frame(ew_data_D)
  for (column in (col_first:col_last)){
    ts <- df[,column]
    dens <- NA
    x <- NA
    y <- NA
    s <- xmax-xmin
    y <- pracma::linspace(xmin, xmax, win)
    for (i in (1:(length(ts)-(win-1)))){
      x <- ts[i:(i+win-1)]
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
      ew_data_D[(i+win-1),(column-col_first+1)] <- 1-(r/g)
    }
  }
  ew_data_D <- ew_data_D[(1:nrow(df)),]
  return(ew_data_D)
}

# function to compute critical instability


#' Critical Instability
#'
#' @param df Data frame with output fro `dyn_comp()``
#' @param win Window used
#'
#' @return Critical instability
#'
#' @name crit_in-deprecated
#' @seealso \code{\link{casnet-deprecated}}
#' @keywords internal
NULL

#' @rdname casnet-deprecated
#' @section \code{crit_in}:
#' For \code{crit_in}, use \code{\link{dc_ccp}}.
#'
#' @export
crit_in = function(df, win){
  df.ci <- matrix(data=NA, nrow=nrow(df), ncol=ncol(df))
  df.ci <- data.frame(df.ci)
  colnames(df.ci) <- colnames(df)
  newcol <- 1

  for (column in (1:ncol(df))){

    item_z <- NA
    item_z <- scale(df[,column])

    ci_index <- which(item_z > 1.645)
    item_z_ci <- rep(0, length(item_z))
    item_z_ci[ci_index] <- 1
    df.ci[,newcol] <- item_z_ci
    newcol <- newcol+1
  }

  sum_ci <- rowSums(df.ci)
  sum_ci_z <- scale(sum_ci)
  sum_ci_index <- which(sum_ci_z > 1.645)
  CriticalIn <- rep(0, nrow(df))
  CriticalIn[sum_ci_index] <- 1
  df.ci$CritInSum <- CriticalIn

  mat.ci <- data.matrix(df.ci)
  colnames(mat.ci) <- c(1:ncol(mat.ci))
  mat.ci[,ncol(mat.ci)] <- mat.ci[,ncol(mat.ci)]*10 #the last column, showing Critical Instability of all items has values of 10 instead of 1

  plot.ci <- ggplot2::ggplot(reshape2::melt(mat.ci), aes_(x=~Var1, y=~Var2, fill=as.factor(~value))) + ggplot2::geom_tile() +
    ggplot2::scale_fill_manual(values = c("NA", "grey", "black")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab('Days') +
    ggplot2::ylab('Items')+
    ggplot2::theme(panel.grid.minor = element_line(colour="NA", size=0.5)) +
    ggplot2::scale_y_continuous(minor_breaks=seq(0.5, ncol(mat.ci), 1), breaks=seq(0.5, ncol(mat.ci)+0.5, 1), label=c(0:(ncol(mat.ci)-1),"Sum")) +
    ggplot2::scale_x_continuous(minor_breaks=seq(0.5,nrow(mat.ci),1), breaks=seq(0.5,nrow(mat.ci)+0.5,1), label=c(0:nrow(mat.ci))) +
    ggplot2::ggtitle('Critical Instability Plot') +
    ggplot2::coord_cartesian(ylim=c(0.5,(ncol(mat.ci)+0.5)), xlim=c(win-0.5,nrow(mat.ci)+0.5), expand=F)

  return(list(df.ci = df.ci,
              plot.ci = plot.ci))
}