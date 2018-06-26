# CasNET Functies Merlijn 03-2018
#
# Hoofd Functie 1: dyn_comp()
# input: dataframe, window size, minimum value of measurement scale, maximum value of maximum scale, first column to be evaluated, last column to be evaluated.
# output: dataset met dynamic complexity values, complexity resonance diagram
# deze functie gebruikt DC_F(), DC_D(), Heaviside()
#
# Hoofd functie 2: crit_in()
# input: dataframe, window size, minimum value of measurement scale, maximum value of maximum scale, first column to be evaluated, last column to be evaluated.
# output: dataset met critical instability values per item en de som, critical instability diagram
#


#' @title Dynamic Complexity
#'
#' @description Calculates Dynamic Complexity, a complexity index for short and coarse-grained time series (Schiepek & Strunk, 2010; Schiepek, 2003; Haken & Schiepek 2006).
#'
#' @param df A dataframe containing multivariate time series data from 1 person. Rows should indicate time,  columns should indicate variables. The multivariate time series should be on the same scale. If nescessary, rescale variables.
#' @param col_first The first column of the dataframe that should be analyzed.
#' @param col_last The last column of the dataframe that should be analyzed.
#' @param win Window size (default = \code{NROW(df)}
#' @param scale_min Theoretical minimum value of thescale.
#' @param scale_max Theoretical maximum value of the scale.
#' @param doPlot If \code{TRUE} shows a plot and returns an invisible \code{\link[ggplot2]{ggplot}} object. If \code{FALSE}, plot will be drawn and no object will be returned (default = \code{FALSE})
#'
#' @return If \code{doPlot = TRUE}, a list object containing a dataframe of dynamic complexity values and a \code{ggplot2} object of the dynamic complexity resonance diagram (e.g. Schiepek et al., 2016). If \code{doPlot = FALSE} the dataframe is returned.
#'
#' @export
#'
#' @author Merlijn Olthof
#'
#' @references Schiepek, G., & Strunk, G. (2010). The identification of critical fluctuations and phase transitions in short term and coarse-grained time series-a method for the real-time monitoring of human change processes. Biological cybernetics, 102(3), 197-207.
#' @references Schiepek, G. (2003). A Dynamic Systems Approach to Clinical Case Formulation. European Journal of Psychological Assessment, 19, 175-184.
#' @references Haken, H. & Schiepek, G. (2006, 2. Aufl. 2010). Synergetik in der Psychologie. Selbstorganisation verstehen und gestalten. G?ttingen: Hogrefe.
#' @references Schiepek, G. K., St?ger-Schmidinger, B., Aichhorn, W., Sch?ller, H., & Aas, B. (2016). Systemic case formulation, individualized process monitoring, and state dynamics in a case of dissociative identity disorder. Frontiers in psychology, 7, 1545.
#'
dyn_comp = function(df, col_first, col_last, scale_min, scale_max, win=NROW(df),  doPlot = FALSE){


if(win>0){ok=TRUE}

  #if(df>scale_max){ok=FALSE}

  data_f <- DC_F(df=df, win=win, xmin=scale_min, xmax=scale_max, col_first =col_first, col_last = col_last)

  data_d <- DC_D(df=df, win=win, xmin=scale_min, xmax=scale_max, col_first =col_first, col_last = col_last)


  df.comp <- matrix(NA, nrow=nrow(data_f), ncol=ncol(data_f))
  for (column in (1:ncol(data_f))){
    df.comp[,column] <- data_f[,column]*data_d[,column]
  }

  mat.dc <- data.matrix(df.comp)
  colnames(mat.dc) <- c(1:ncol(mat.dc))


if(doPlot){
  plot.dc <- ggplot2::ggplot(reshape2::melt(mat.dc), aes_(x=~Var1, y=~Var2, fill=~value)) + ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low='blue', high='red', mid='yellow', midpoint=(max(mat.dc, na.rm=TRUE)/2), na.value='white') +
    ggplot2::theme_bw() +
    ggplot2::xlab('Days') +
    ggplot2::ylab('Items')+
    ggplot2::theme(panel.grid.minor = element_line(colour="NA", size=0.5)) +
    ggplot2::scale_y_continuous(minor_breaks=seq(0.5, ncol(mat.dc), 1), breaks=seq(0.5, ncol(mat.dc)+0.5, 1), label=c(0:(ncol(mat.dc)))) +
    ggplot2::scale_x_continuous(minor_breaks=seq(0.5,nrow(mat.dc),1), breaks=seq(0.5,nrow(mat.dc)+0.5,1), label=c(0:nrow(mat.dc))) +
    ggplot2::ggtitle('Complexity Resonance Diagram') +
    ggplot2::coord_cartesian(ylim=c(0.5,(ncol(mat.dc)+0.5)), xlim=c(win-0.5,nrow(mat.dc)+0.5), expand=FALSE)

  graphics::plot.new()
  graphics::plot(plot.dc)
  return(list(df.comp = df.comp,
              plot.dc = plot.dc))
} else {
    return(df.comp)
  }
}



#' Fluctuation Intensity
#'
#' Used to calculate dynamic complexity
#' @param df df
#' @param win win
#' @param xmin xmin
#' @param xmax xmax
#' @param col_first col_first
#' @param col_last col_last
#'
#' @return dataframe
#' @export
#' @keywords internal
#'
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


###Heaviside step function###
#' Heaviside step function
#'
#' @param value value
#'
#' @return heaviside change
#' @export
#' @keywords internal
#'
Heaviside = function(value){

  if (value>0){
    h=1
  }
  else {h=0}
  return(h)
}



#' Distribution Uniformity
#'
#' @param data data
#' @param win win
#' @param xmin xmin
#' @param xmax xmax
#' @param col_first col_first
#' @param col_last col_last
#'
#' @return a dataframe
#' @export
#' @keywords internal
#'
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

#' @title Critical Instability
#' @description Computes significant peaks in the dynamic complexity time series. Example: Schiepek, Tominschek & Heinzel, 2014.
#'
#' @param df A dataframe of dynamic complexity values, e.g. output data from dyn_comp
#' @param win Window size over which dynamic complexity was calculated (neccesary for plot)
#' @param doPlot If \code{TRUE} shows a plot and returns a \code{\link[ggplot2]{ggplot}} object in a list. If \code{FALSE}, no plot will be drawn and no object will be returned (default = \code{FALSE})
#'
#' @return A list with a dataframe of binary critical instability indices and a summed critical intability index, a ggplot object containing a critical instability diagram.
#'
#' @export
#'
#' @author Merlijn Olthof
#'
#' @references Schiepek, G., & Strunk, G. (2010). The identification of critical fluctuations and phase transitions in short term and coarse-grained time series-a method for the real-time monitoring of human change processes. Biological cybernetics, 102(3), 197-207.
#' @references Schiepek, G. (2003). A Dynamic Systems Approach to Clinical Case Formulation. European Journal of Psychological Assessment, 19, 175-184.
#' @references Haken, H. & Schiepek, G. (2006, 2. Aufl. 2010). Synergetik in der Psychologie. Selbstorganisation verstehen und gestalten. G?ttingen: Hogrefe.
#' @references Schiepek, G. K., Tominschek, I., & Heinzel, S. (2014). Self-organization in psychotherapy: testing the synergetic model of change processes. Frontiers in psychology, 5, 1089.
#'
crit_in = function(df, win, doPlot = TRUE){

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

if(doPlot){
  plot.ci <- ggplot2::ggplot(reshape2::melt(mat.ci), aes_(x=~Var1, y=~Var2, fill=as.factor(~value))) +
    ggplot2::geom_tile() +
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

   graphics::plot.new()
   graphics::plot(plot.ci)

  return(list(df.ci = df.ci,
              plot.ci = plot.ci))
} else {
    return(df.ci)
  }
}
