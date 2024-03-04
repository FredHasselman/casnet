# package casnet ----
#
# Long Range Dependence - (Multi-) Fractal Estimators ----
#
# fd functions


FDrel <- function(g){
  d<-igraph::degree(g,mode="all")
  nbreaks <- round(length(igraph::V(g))/2)-1
  y<-graphics::hist(d,breaks=nbreaks,plot=FALSE)$density
  y<-y[y>0]
  return(FD <- -sum(y*log2(y))/-(log2(1/length(y))))
}


#' Informed Dimension estimate from Spectral Slope (aplha)
#'
#' @description Conversion formula: From periodogram based self-affinity parameter estimate (`sa`) to an informed estimate of the (fractal) dimension (FD).
#' @param sa Self-Affinity parameter estimate based on PSD slope (e.g., [fd_psd()]))
#' @param ... Other arguments
#'
#' @return An informed estimate of the Fractal Dimension, see Hasselman(2013) for details.
#' @export
#'
#' @details The spectral slope will be converted to a dimension estimate using:
#'
#' \deqn{D_{PSD}\approx\frac{3}{2}+\frac{14}{33}*\tanh\left(Slope * \ln(1+\sqrt{2})\right)}
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. https://doi.org/10.3389/fphys.2013.00075
#'
#'
sa2fd_psd <- function(sa, ...){return(round(3/2 + ((14/33)*tanh(sa*log(1+sqrt(2)))), digits = 2))}


#' Informed Dimension estimate from DFA slope (H)
#'
#' @description Conversion formula: Detrended Fluctuation Analysis (DFA) estimate of the Hurst exponent (a self-affinity parameter `sa`) to an informed estimate of the (fractal) dimension (FD).
#'
#' @param sa Self-Afinity parameter estimate based on DFA slope (e.g., [fd_sda()])).
#' @param ... Other arguments
#'
#' @return An informed estimate of the Fractal Dimension, see Hasselman(2013) for details.
#'
#' @export
#'
#' @details The DFA slope (H) will be converted to a dimension estimate using:
#'
#' \deqn{D_{DFA}\approx 2-(\tanh(\log(3)*sa)) }{D_{DFA} ~ 2-(tanh(log(3)*sa)) }
#'
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. https://doi.org/10.3389/fphys.2013.00075
#'
sa2fd_dfa <- function(sa, ...){return(round(2-(tanh(log(3)*sa)), digits = 2))}


#' Informed Dimension estimate from SDA slope.
#'
#' @description Conversion formula: Standardised Dispersion Analysis (SDA) estimate of self-affinity parameter (`SA`) to an informed estimate of the fractal dimension (FD).
#'
#' @param sa Self-afinity parameter estimate based on SDA slope (e.g., [fd_sda()])).
#' @param ... Other arguments
#'
#' @details
#'
#'
#' Note that for some signals different PSD slope values project to a single SDA slope. That is, SDA cannot distinguish dplyr::between all variaties of power-law scaling in the frequency domain.
#'
#' @return An informed estimate of the Fractal Dimension, see Hasselman(2013) for details.
#' @export
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. https://doi.org/10.3389/fphys.2013.00075
#'
#'
sa2fd_sda <- function(sa, ...){return(1-sa)}



#' Relative Roughness
#'
#' Relative Rougness is a ratio of local variance (autocovariance at lag-1) to global variance (autocovariance at lag-0) that can be used to classify different 'noises'.
#'
#'\deqn{RR = 2 * \left[1 - \frac{\gamma(y)}{Var(y)}\right]}{RR = 2 * [ 1 - autoCov-1(y) / Var(y) ]}
#'
#' @param y A numeric vector.
#'
#' @return The Relative Roughness of y, the values of local and global variance are returned as attributes
#' @export
#'
#' @family Fluctuation Analyses
#'
#' @references
#' \itemize{
#' \item{Marmelat, V., Torre, K., & Delignieres, D. (2012). Relative roughness: an index for testing the suitability of the monofractal model. *Frontiers in Physiology, 3*, 208.}}
#'
#'
fd_RR <- function(y){
  # lag.max = n gives autocovariance of lags 0 ... n,
  VAR  <- stats::acf(y, lag.max = 1, type = 'covariance', plot=FALSE)
  # RR formula
  RelR   <- 2*(1 - VAR$acf[2] / VAR$acf[1])
  # Add some attributes to the output
  attributes(RelR) <- list(localAutoCoVariance = VAR$acf[2], globalAutoCoVariance = VAR$acf[1])

  return(RelR)
}


#' @title Power Spectral Density Slope (PSD).

#' @description Estimate Alpha, Hurst Exponent and Fractal Dimension through log-log slope.
#'
#' @param y    A numeric vector or time series object.
#' @param fs Sample rate (default = `NULL`)
#' @param standardise    standardise the series (default = `TRUE`).
#' @param detrend    Subtract linear trend from the series (default = `TRUE`).
#' @param fitMethod Method to decide on a frequency range for log-log fit. Can be one of: "lowest25","Wijnants","Hurvich-Deo" (default). See details for more info.
#' @param doPlot    Return the log-log spectrum with linear fit (default = `TRUE`).
#' @param returnPlot Return ggplot2 object (default = `FALSE`)
#' @param returnPLAW Return the power law data (default = `FALSE`)
#' @param returnInfo Return all the data used in SDA (default = `FALSE`)
#' @param silent Run in silent-ish mode (default = `TRUE)`
#' @param noTitle Do not generate a title (only the subtitle)
#' @param tsName Name of y added as a subtitle to the plot
#'
#' @author Fred Hasselman
#'
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. https://doi.org/10.3389/fphys.2013.00075
#' @references Hurvich, C.M., & Deo, R.R. (1999). Plug-in Selection of the Number of Frequencies in Regression Estimates of the Memory Parameter of a Long Memory Time Series. *Journal of Time Series Analysis, 20(3)*, 331–341.
#'
#' @return A list object containing:
#' \itemize{
#' \item A data matrix `PLAW` with columns `freq.norm`, `size` and `bulk`.
#' \item Estimate of scaling exponent `alpha` based on a fit over the lowest 25\% frequencies (`low25`), or using the HD estimate `HD`.
#' \item Estimate of the the Fractal Dimension (`FD`) using conversion formula's reported in Hasselman(2013).
#' \item Information output by various functions.
#' }
#'
#' @family Fluctuation Analyses
#'
#' @export
#'
#' @details Calls function [stats::spec.pgram()] to estimate the scaling exponent of a timeseries based on the periodogram frequency spectrum. After detrending and normalizing the signal (if requested), [stats::spec.pgram()] is called using a cosine taper = 0.5.
#'
#' A line is fitted on the periodogram in log-log coordinates. The full range is fitted as well as one of three fit-ranges:
#'
#' * `lowest25` - The 25\% lowest frequencies
#' * `Wijnants` - The 50 lowest frequencies (Wijnants et al., 2012)
#' * `Hurvich-Deo` - The Hurvich-Deo estimate (Hurvich & Deo, 1999)
#'
#'
fd_psd <- function(y,
                   fs = NULL,
                   removeTrend = c("no","poly","adaptive","bridge")[2],
                   polyOrder=1,
                   standardise = c("none","mean.sd","median.mad")[2],
                   fitMethod = c("lowest25","Wijnants","Hurvich-Deo")[3],
                   doPlot = FALSE,
                   returnPlot = FALSE,
                   returnPLAW = FALSE,
                   returnInfo = FALSE,
                   silent = FALSE,
                   noTitle = FALSE,
                   tsName="y"){


  y <- fd_prepSeries(y = y,
                     fs = fs,
                     removeTrend = removeTrend,
                     polyOrder = polyOrder,
                     standardise = standardise,
                     adjustSumOrder = FALSE,
                     scaleMin = NA,
                     scaleMax = NA,
                     scaleResolution = NA,
                     dataMin = NA,
                     scaleS = NULL,
                     overlap = NA,
                     silent = silent
  )

  scaleS  <- attr(y, "scaleS")
  dataMin <- attr(y, "dataMin")

  N <- NROW(y)

  # Number of frequencies estimated cannot be set! (defaults to Nyquist)
  # Use Tukey window: cosine taper with r = 0.5

  # fast = TRUE ensures padding with zeros to optimize FFT to highly composite number.
  # However, we just pad to nextPow2, except if length already is a power of 2.
  # npad <- 1+(stats::nextn(N,factors=2)-N)/N
  npad <- stats::nextn(N, factors = 2)

  if(N==npad){npad <- 0}
  psd  <- stats::spec.pgram(y, fast = FALSE, demean=FALSE, detrend=FALSE, plot=FALSE, pad=npad, taper=0.5)

  # Tukey <- sapa::taper(type="raised cosine", flatness = 0.5, n.sample = N)
  # psd   <- sapa::SDF(y, taper. = Tukey, npad = npad)

 Nfreq <- NROW(psd$freq)
 freq.norm <- psd$freq/stats::frequency(y)
 size <- psd$freq
 bulk <- 2*psd$spec
 #plot(x=log2(psd$freq), y=log2(psd$spec*2),pch=".")

  powspec <- cbind.data.frame(freq.norm = freq.norm, size = size, bulk = as.matrix(bulk))

  # First check the global slope for anti-persistent noise (GT +0.20)
  # If so, fit the line starting from the highest frequency
  nr     <- length(powspec[,1])
  lsfit  <- stats::lm(log(powspec$bulk[1:nr]) ~ log(powspec$size[1:nr]))
  glob   <- stats::coef(lsfit)[2]

  # General guideline: fit over 25% frequencies
  # If signal is continuous (sampled) consider Wijnants et al. (2013) log-log fitting procedure
  nr <- switch(fitMethod,
               "lowest25" = which(powspec$size>=0.25)[1],
               "Hurvich-Deo" = HurvichDeo(nr = nr, spec = as.vector(powspec$bulk)),
               "Wijnants" = 50
  )

 if(nr>=length(powspec$freq.norm)){nr <- length(powspec$freq.norm)-1}

  ifelse((glob > 0.2), {
    lmfit1 <- stats::lm(log(rev(powspec$bulk)) ~ log(rev(powspec$size)))
    lmfit2 <- stats::lm(log(rev(powspec$bulk[1:nr])) ~ log(rev(powspec$size[1:nr])))
  },{
    lmfit1 <- stats::lm(log(powspec$bulk) ~ log(powspec$size))
    lmfit2 <- stats::lm(log(powspec$bulk[1:nr]) ~ log(powspec$size[1:nr]))
  })


  if(N>4*50){
    d <- 50
  } else {
    d <- round(N/4)
  }
  exp1 <- pracma::hurstexp(cumsum(y), d = d, display = FALSE)
  # #fractal::hurstSpec(y, sdf.method="direct", freq.max = powspec$freq.norm[length(powspec$freq.norm)-1], taper.=Tukey )
  exp2 <- pracma::hurstexp(cumsum(y), d = nr+1, display = FALSE)
  # # fractal::hurstSpec(y, sdf.method="direct", freq.max = powspec$freq.norm[nr], taper.=Tukey)


  outList <- list(
    PLAW  = powspec,
    fullRange = list(sap = stats::coef(lmfit1)[2],
                     H = exp1[1],
                     FD = sa2fd_psd(stats::coef(lmfit1)[2]),
                     fitlm1 = lmfit1,
                     method = paste0("All frequencies (n = ",length(powspec$freq.norm),")\nSlope = ",round(stats::coef(lmfit1)[2],2)," | FD = ",sa2fd_psd(stats::coef(lmfit1)[2]))),
    fitRange  = list(sap = stats::coef(lmfit2)[2],
                     H = exp2[1],
                     FD = sa2fd_psd(stats::coef(lmfit2)[2]),
                     fitlm2 = lmfit2,
                     method = paste0(fitMethod," (n = ",nr,")\nSlope = ",round(stats::coef(lmfit2)[2],2)," | FD = ",sa2fd_psd(stats::coef(lmfit2)[2]))),
    info = psd,
    plot = NA,
    analysis = list(
      name = "Power Spectral Density Slope",
      logBaseFit = "log",
      logBasePlot = 10)
  )

  if(doPlot|returnPlot){
    if(noTitle){
      title <- " "
    } else {
      title <- "log-log regression (PSD)"
    }
    g <- plotFD_loglog(fd.OUT = outList, title = title, subtitle = tsName, logBase = "10",xlabel = "Normalised Frequency", ylabel = "Power", doPlot = doPlot, returnPlot = returnPlot)


      if(returnPlot){
        outList$plot <- g
      }
    }


  if(returnInfo){returnPLAW<-TRUE}

  if(silent==FALSE){
    cat("\n~~~o~~o~~casnet~~o~~o~~~\n")
    cat(paste("\n",outList$analysis$name,"\n\n",outList$fullRange$method,"\n\n",outList$fitRange$method))
    cat("\n\n~~~o~~o~~casnet~~o~~o~~~\n")
  }

  return(invisible(outList[c(returnPLAW,TRUE,TRUE,returnInfo,returnPlot,TRUE)]))
}


#' Prepare time series for fluctuation analysis
#'
#' @param y y
#' @param fs fs
#' @param removeTrend rt
#' @param polyOrder po
#' @param standardise st
#' @param adjustSumOrder ao
#' @param scaleMin sm
#' @param scaleMax smm
#' @param scaleResolution sr
#' @param dataMin dm
#' @param scaleS sc
#' @param overlap ov
#' @param silent si
#'
#' @keywords internal
#'
#' @return transformed series
#'
#' @export
#'
fd_prepSeries <- function(y,
                          fs = NULL,
                          removeTrend = c("no","poly","adaptive","bridge")[2],
                          polyOrder=1,
                          standardise = c("none","mean.sd","median.mad")[2],
                          adjustSumOrder = FALSE,
                          scaleMin = NA,
                          scaleMax = NA,
                          scaleResolution = NA,
                          dataMin = NA,
                          scaleS = NA,
                          overlap = NA,
                          silent = TRUE){

  if(any(is.na(y))){
    warning("Missing values in y! Analysis may throw an error.")
  }

  y_ori <- y

  if(is.na(scaleMin)){
    scaleMin <- 16
  }

  if(is.na(scaleMax)){
    scaleMax <- floor(NROW(y)/2)
  }

  if(is.na(scaleResolution)){
    scaleResolution <- stats::nextn(log2(scaleMax)-log2(scaleMin), factors = 2)
  }

  if(any(is.na(scaleS))){
    exponents <- seq(log2(scaleMin), log2(scaleMax), length.out = scaleResolution+1)
    scaleS    <- round(2^exponents) #unique(round(seq(scaleMin,scaleMax, length.out= (scaleMax-scaleMin)/(scaleResolution-1))))
  } else {
    if(is.null(scaleS)){
      scaleS <- NA
    }
  }


  if(!stats::is.ts(y)){
    if(is.null(fs)){
      fs <- 1
      if(!silent){cat("\n\n(mf)dfa:\tSample rate was set to 1.\n\n")}
    }
    y <- stats::ts(y, frequency = fs)
  }


  if(!all(is.na(scaleS))){
    # Nyquist-ish
    if(max(scaleS)>floor(NROW(y)/2)){
      warning("The maximum bin should be smaller than floor(NROW(y)/2) data points")
    }


  if(!all(is.numeric(scaleS),length(scaleS)>0,scaleS%[]%c(2,(NROW(y)/2)))){
    stop("Something wrong with vector passed to scaleS.... \nUsing defaults: (scaleMax-scaleMin)/scaleResolution")
  }
  }

  if(!removeTrend%in%"no"){
    if(removeTrend%in%"adaptive"){
      adaptive <- TRUE
    } else {
      adaptive <- FALSE
    }
    y <- ts_detrend(y, polyOrder = polyOrder, adaptive = adaptive)
  }

  # Standardise by N
  if(any(standardise%in%c("mean.sd","median.mad"))){
    y <- ts_standardise(y, type = standardise,  adjustN = FALSE)
  }

  if(is.na(dataMin)){
    dataMin <- ceiling(NROW(y)/max(scaleS, na.rm = TRUE))+1
  }

  # if(adjustSumOrder){
  #   y       <- ts_sumorder(y_ori, scaleS = scaleS, polyOrder = polyOrder, dataMin = dataMin)
  #   Hadj    <- attr(y,"Hadj")
  #   Hglobal <- attr(y,"Hglobal.excl")
  # } else {
  #   Hadj    <- 0
  #   Hglobal <- NA
  # }

  attr(y,"Original time series") <- y_ori

  attr(y,"removeTrend") <- removeTrend
  attr(y,"polyOrder")   <- polyOrder
  attr(y,"standardise") <- standardise
  attr(y,"scaleMin")    <- scaleMin
  attr(y,"scaleMax")    <- scaleMax
  attr(y,"scaleResolution") <- scaleResolution
  attr(y,"dataMin")     <- dataMin
  attr(y,"scaleS")      <- scaleS
  attr(y,"overlap")     <- overlap

  attr(y,"adjustSumOrder") <- adjustSumOrder
  if(adjustSumOrder){
    y <- ts_sumorder(y, scaleS = scaleS, polyOrder = polyOrder)
    Hglobal.full <- attr(y,"Hglobal.full")
    Hglobal.excl <- attr(y,"Hglobal.excl")
    Hadj <- attr(y,"Hadj")
  } else {
    Hglobal.full <- attr(y,"Hglobal.full") <- NA
    Hglobal.excl <- attr(y,"Hglobal.excl") <- NA
    Hadj <- attr(y,"Hadj") <- 0
  }

  return(y)
}


#' fd_sda
#'
#' @title Standardised Dispersion Analysis (SDA).
#'
#' @inheritParams fd_dfa
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. https://doi.org/10.3389/fphys.2013.00075
#'
#' @return A list object containing:
#' \itemize{
#' \item A data matrix `PLAW` with columns `freq.norm`, `size` and `bulk`.
#' \item Estimate of scaling exponent `sap` based on a fit over the standard range (`fullRange`), or on a user defined range `fitRange`.
#' \item Estimate of the the Fractal Dimension (`FD`) using conversion formula's reported in Hasselman(2013).
#' \item Information output by various functions.
#' }
#'
#' @export
#'
#' @family Fluctuation Analyses
#'
#'
#'
fd_sda <- function(y,
                   fs = NULL,
                   removeTrend = c("no","poly","adaptive","bridge")[2],
                   polyOrder=1,
                   standardise = c("none","mean.sd","median.mad")[2],
                   adjustSumOrder = FALSE,
                   scaleMin = 4,
                   scaleMax = stats::nextn(floor(NROW(y)/2), factors = 2),
                   scaleResolution = log2(scaleMax)-log2(scaleMin),
                   dataMin = NA,
                   scaleS = NA,
                   overlap = 0,
                   doPlot = FALSE,
                   returnPlot = FALSE,
                   returnPLAW = FALSE,
                   returnInfo = FALSE,
                   silent = FALSE,
                   noTitle = FALSE,
                   tsName="y"){


  y <- fd_prepSeries(y = y,
                     fs = fs,
                     removeTrend = removeTrend,
                     polyOrder = polyOrder,
                     standardise = standardise,
                     adjustSumOrder = adjustSumOrder,
                     scaleMin = scaleMin,
                     scaleMax = scaleMax,
                     scaleResolution = scaleResolution,
                     dataMin = dataMin,
                     scaleS = scaleS,
                     overlap = overlap,
                     silent = silent
  )


  scaleS  <- attr(y, "scaleS")
  dataMin <- attr(y, "dataMin")

  Hadj    <- attr(y,"Hadj")
  Hglobal <- attr(y,"Hglobal.excl")

  # if(!stats::is.ts(y)){
  #   if(is.null(fs)){fs <- 1}
  #   y <- stats::ts(y, frequency = fs)
  #   if(!silent){cat("\n\nfd_sda:\tSample rate was set to 1.\n\n")}
  # }
  #
  #
  # N             <- length(y)
  #
  # if(is.na(scaleS)){
  #   scaleS <- unique(round(2^(seq(scaleMin, scaleMax, by=((scaleMax-scaleMin)/scaleResolution)))))
  # }

  if(max(scaleS)>(NROW(y)/2)){
    scaleS <- scaleS[scaleS<=(NROW(y)/2)]
  }

  # if(!all(is.numeric(scaleS),length(scaleS)>0,scaleS%[]%c(2,(NROW(y)/2)))){
  #   message("Something wrong with vector passed to scaleS.... \nUsing default: (scaleMax-scaleMin)/scaleResolution")
  # }
  #
  # # Standardise by N
  # if(any(standardise%in%c("mean.sd","median.mad"))){
  #   y <- ts_standardise(y, type = standardise,  adjustN = FALSE)
  # }
  #
  # if(adjustSumOrder){
  #   y       <- ts_sumorder(y, scaleS = scaleS, polyOrder = polyOrder, dataMin = dataMin)
  #   Hadj    <- attr(y,"Hadj")
  #   Hglobal <- attr(y,"Hglobal.excl")
  # } else {
  #   Hadj    <- 0
  #   Hglobal <- NA
  # }

  out <- SDA(y, front = FALSE, scaleS)

  scale <- out$scale[out$scale%[]%c(scaleMin,scaleMax)]
  sd    <- out$sd[out$scale%[]%c(scaleMin,scaleMax)]

  fitRange1 <- which(lengths(lapply(scale, function(s){ts_slice(y,s)}))<NROW(y))
  fitRange2 <- which(lengths(lapply(scale, function(s){ts_slice(y,s)}))%[)%c(dataMin,NROW(y)))

  lmfit1        <- stats::lm(log(sd[fitRange1]) ~ log(scale[fitRange1]))
  lmfit2        <- stats::lm(log(sd[fitRange2]) ~ log(scale[fitRange2]))

  outList <- list(
    PLAW  =  cbind.data.frame(freq.norm = stats::frequency(y)/scale[fitRange1],
                              size = scale[fitRange1],
                              bulk = sd[fitRange1]),
    fullRange = list(sap = stats::coef(lmfit1)[2],
                     H = 1+stats::coef(lmfit1)[2] + Hadj,
                     FD = sa2fd_sda(stats::coef(lmfit1)[2]),
                     fitlm1 = lmfit1,
                     method = paste0("Full range (n = ",length(scale[fitRange1]), ")\nSlope = ", round(stats::coef(lmfit1)[2],2)," | FD = ", round(sa2fd_sda(stats::coef(lmfit1)[2]),2))),
    fitRange  = list(sap = stats::coef(lmfit2)[2],
                     H = 1+stats::coef(lmfit2)[2] + Hadj,
                     FD = sa2fd_sda(stats::coef(lmfit2)[2]),
                     fitlm2 = lmfit2,
                     method = paste0("Fit range (n = ", length(scale[fitRange2]),")\nSlope = ",round(stats::coef(lmfit2)[2],2)," | FD = ", round(sa2fd_sda(stats::coef(lmfit2)[2]),2))),
    info = out,
    plot = NA,
    analysis = list(
      name = "Standardised Dispersion Analysis",
      logBaseFit = "log",
      logBasePlot = "e")
  )

  if(doPlot|returnPlot){
    if(noTitle){
      title <- " "
    } else {
      title <- "log-log regression (SDA)"
    }
    g <- plotFD_loglog(fd.OUT = outList, title = title, subtitle = tsName, logBase = "e", ylabel = "Standardised Dispersion", doPlot = doPlot, returnPlot = returnPlot)

      if(returnPlot){
        outList$plot <- g
      }
    }


  if(returnInfo){returnPLAW<-TRUE}

  if(!silent){
    cat("\n~~~o~~o~~casnet~~o~~o~~~\n")
    cat(paste("\n",outList$analysis$name,"\n\n",outList$fullRange$method,"\n\n",outList$fitRange$method))
    cat("\n\n~~~o~~o~~casnet~~o~~o~~~\n")
  }

  return(invisible(outList[c(returnPLAW,TRUE,TRUE,returnInfo,returnPlot, TRUE)]))

}


#' fd_dfa
#'
#' @title Detrended Fluctuation Analysis (DFA)
#'
#' @param y    A numeric vector or time series object.
#' @param fs   Sample rate
#' @param removeTrend Method to use for global detrending (default = `"poly"`)
#' @param polyOrder Order of global polynomial trend to remove if `removeTrend = "poly"`. If `removeTrend = "adaptive"` polynomials `1` to `polyOrder` will be evaluated and the best fitting curve (R squared) will be removed (default = `1`)
#' @param standardise Standardise the series using [casnet::ts_standardise()] with `adjustN = FALSE` (default = "mean.sd")
#' @param adjustSumOrder  Adjust the time series (summation or difference), based on the global scaling exponent, see e.g. <https://www.frontiersin.org/files/Articles/23948/fphys-03-00141-r2/image_m/fphys-03-00141-t001.jpg>{Ihlen (2012)} (default = `FALSE`)
#' @param removeTrendSegment Method to use for detrending in the bins (default = `"poly"`)
#' @param polyOrderSegment The DFA order, the order of polynomial trend to remove from the bin if `removeTrendSegment = "poly"`. If `removeTrendSegment = "adaptive"` polynomials `1` to `polyOrder` will be evaluated and the best fitting polynomial (R squared) will be removed (default = `1`)
#' @param scaleMin   Minimum scale (in data points) to use for log-log regression (default = `4`)
#' @param scaleMax   Maximum scale (in data points) to use for log-log regression. This value will be ignored if `dataMin` is not `NA`, in which case bins of size `< dataMin` will be removed (default = `stats::nextn(floor(NROW(y)/4), factors = 2)`)
#' @param scaleResolution  The scales at which detrended fluctuation will be evaluated are calculated as: `seq(scaleMin, scaleMax, length.out = scaleResolution)` (default =  `round(log2(scaleMax-scaleMin))`).
#' #' @param dataMin Minimum number of data points in a bin required for inclusion in calculation of the scaling relation. For example if `length(y) = 1024` and `dataMin = 4`, the maximum scale used to calculate the slope will be `1024 / 4 = 256`. This value will take precedence over the `scaleMax` (default = `NA`)
#' @param scaleS If not `NA`, it should be a numeric vector listing the scales on which to evaluate the detrended fluctuations. Arguments `scaleMax, scaleMin, scaleResolution` and `dataMin` will be ignored (default = `NA`)
#' @param overlap A number in `[0 ... 1]` representing the amount of 'bin overlap' when calculating the fluctuation. This reduces impact of arbitrary time series begin and end points. If `length(y) = 1024` and overlap is `.5`, a scale of `4` will be considered a sliding window of size `4` with step-size `floor(.5 * 4) = 2`, so for scale `128` step-size will be `64` (default = `NA`)
#' @param y    A numeric vector or time series object.
#' @param doPlot   Output the log-log scale versus fluctuation plot with linear fit by calling function `plotFD_loglog()` (default = `TRUE`)
#' @param returnPlot Return ggplot2 object (default = `FALSE`)
#' @param returnPLAW Return the power law data (default = `FALSE`)
#' @param returnInfo Return all the data used in SDA (default = `FALSE`)
#' @param silent Silent-ish mode (default = `FALSE`)
#' @param noTitle Do not generate a title (only the subtitle) (default = `FALSE`)
#' @param tsName Name of y added as a subtitle to the plot (default = `"y"`)
#'
#' @return Estimate of Hurst exponent (slope of `log(bin)` vs. `log(RMSE))` and an FD estimate based on Hasselman (2013)
#' A list object containing:
#' \itemize{
#' \item A data matrix `PLAW` with columns `freq.norm`, `size` and `bulk`.
#' \item Estimate of scaling exponent `sap` based on a fit over the standard range (`fullRange`), or on a user defined range `fitRange`.
#' \item Estimate of the the Fractal Dimension (`FD`) using conversion formula's reported in Hasselman(2013).
#' \item Information output by various functions.
#' }
#'
#' @export
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. https://doi.org/10.3389/fphys.2013.00075
#'
#' @family Fluctuation Analyses
#'
#' @examples
#'
#' set.seed(1234)
#'
#' # Brownian noise
#' fd_dfa(cumsum(rnorm(512)))
#'
#' # Brownian noise with overlapping bins
#' fd_dfa(cumsum(rnorm(512)), overlap = 0.5)
#'
#' # Brownian noise to white noise - windowed analysis
#' y <- rnorm(1024)
#' y[1:512] <- cumsum(y[1:512])
#'
#' id <- ts_windower(y, win = 256, step = 1)
#'
#' DFAseries <- plyr::ldply(id, function(w){
#' fd <- fd_dfa(y[w], silent = TRUE)
#' return(fd$fitRange$FD)
#' })
#'
#' op <- par(mfrow=c(2,1))
#' plot(ts(y))
#' plot(ts(DFAseries[,2]))
#' lines(c(0,770),c(1.5,1.5), col = "red3")
#' lines(c(0,770),c(1.1,1.1), col = "steelblue")
#' par(op)
#'
fd_dfa <- function(y,
                   fs = NULL,
                   removeTrend = c("no","poly","adaptive","bridge")[2],
                   polyOrder=1,
                   standardise = c("none","mean.sd","median.mad")[2],
                   adjustSumOrder = FALSE,
                   removeTrendSegment = c("no","poly","adaptive","bridge")[2],
                   polyOrderSegment = 1,
                   scaleMin = 4,
                   scaleMax = stats::nextn(floor(NROW(y)/2), factors = 2),
                   scaleResolution = log2(scaleMax)-log2(scaleMin),
                   dataMin = NA,
                   scaleS = NA,
                   overlap = NA,
                   doPlot = FALSE,
                   returnPlot = FALSE,
                   returnPLAW = FALSE,
                   returnInfo = FALSE,
                   silent = FALSE,
                   noTitle = FALSE,
                   tsName="y"){


   y <- fd_prepSeries(y = y,
                 fs = fs,
                 removeTrend = "no",
                 polyOrder = polyOrder,
                 standardise = standardise,
                 adjustSumOrder = adjustSumOrder,
                 scaleMin = scaleMin,
                 scaleMax = scaleMax,
                 scaleResolution = scaleResolution,
                 dataMin = dataMin,
                 scaleS = scaleS,
                 overlap = overlap,
                 silent = silent
                 )

  # Integrate the series
  if(standardise%in%"none"){
    y <- ts_integrate(ts_center(y)) # need negative values in profile
  } else {
    y <- ts_integrate(y)
  }

  scaleS  <- attr(y, "scaleS")
  dataMin <- attr(y, "dataMin")

  if(removeTrendSegment%in%"adaptive"){
    adaptive <- TRUE
  } else {
    adaptive <- FALSE
  }

  TSm    <- as.matrix(cbind(t=1:NROW(y),y=y))
  DFAout <- monoH(TSm = TSm, scaleS = scaleS, removeTrend = removeTrendSegment, polyOrder = polyOrderSegment, overlap = overlap, returnPLAW = TRUE, returnSegments = TRUE)

  fitRange <- which(lapply(DFAout$segments,NROW)>=dataMin)

  lmfit1 <- stats::lm(DFAout$PLAW$bulk.log2 ~ DFAout$PLAW$size.log2, na.action=stats::na.omit)
  H1     <- lmfit1$coefficients[2] + attr(y,"Hadj")
  lmfit2 <- stats::lm(DFAout$PLAW$bulk.log2[fitRange] ~ DFAout$PLAW$size.log2[fitRange], na.action=stats::na.omit)
  H2     <- lmfit2$coefficients[2] + attr(y,"Hadj")

  outList <- list(
    PLAW  =  DFAout$PLAW,
    fullRange = list(y = y,
                     sap = lmfit1$coefficients[2],
                     Hadj = attr(y,"Hglobal.full"),
                     H = H1,
                     FD = sa2fd_dfa(lmfit1$coefficients[2]),
                     fitlm1 = lmfit1,
                     method = paste0("Full range (n = ",NROW(DFAout$PLAW$size),")\nSlope = ", round(stats::coef(lmfit1)[2],2)," | FD = ", sa2fd_dfa(stats::coef(lmfit1)[2]))),
    fitRange  = list(y = y,
                     sap = lmfit2$coefficients[2],
                     H = H2,
                     Hadj = attr(y,"Hglobal.excl"),
                     FD = sa2fd_dfa(lmfit2$coefficients[2]),
                     fitlm2 = lmfit2,
                     method = paste0("Exclude large bin sizes (n = ",NROW(fitRange),")\nSlope = ", round(stats::coef(lmfit2)[2],2)," | FD = ", sa2fd_dfa(stats::coef(lmfit2)[2]))),
    info = list(fullRange=lmfit1,fitRange=lmfit2,segments=DFAout$segments),
    plot = NA,
    analysis = list(
      name = "Detrended Fluctuation Analysis",
      logBaseFit = "log2",
      logBasePlot = "2")
  )

  if(doPlot|returnPlot){
    if(noTitle){
      title <- " "
    } else {
      title <- "log-log regression (DFA)"
    }

    g <- plotFD_loglog(fd.OUT = outList, title = title, subtitle = tsName, logBase = "2",xlabel = "Scale", ylabel = "Detrended Fluctuation", doPlot = doPlot, returnPlot = returnPlot)

    if(returnPlot){
        outList$plot <- g
      }
    }

  if(returnInfo){returnPLAW<-TRUE}

  if(!silent){
    cat("\n~~~o~~o~~casnet~~o~~o~~~\n")
    cat(paste("\n",outList$analysis$name,"\n\n",outList$fullRange$method,"\n\n",outList$fitRange$method,"\n\n Detrending:",removeTrendSegment))
    cat("\n\n~~~o~~o~~casnet~~o~~o~~~\n")
  }

  return(invisible(outList[c(returnPLAW,TRUE,TRUE,returnInfo,returnPlot,TRUE)]))
}


# switch dim
#
# case 1        #------------------- 1D boxcount ---------------------#
#
# n(p+1) = sum(c);
# for g=(p-1):-1:0
# siz = 2^(p-g);
# siz2 = round(siz/2);
# for i=1:siz:(width-siz+1)
# c(i) = ( c(i) || c(i+siz2));
# end
# n(g+1) = sum(c(1:siz:(width-siz+1)));
# end
#
# case 3         #------------------- 3D boxcount ---------------------#
#
# n(p+1) = sum(c(:));
# for g=(p-1):-1:0
# siz = 2^(p-g);
# siz2 = round(siz/2);
# for i=1:siz:(width-siz+1),
# for j=1:siz:(width-siz+1),
# for k=1:siz:(width-siz+1),
# c(i,j,k)=( c(i,j,k) || c(i+siz2,j,k) || c(i,j+siz2,k) ...
# || c(i+siz2,j+siz2,k) || c(i,j,k+siz2) || c(i+siz2,j,k+siz2) ...
# || c(i,j+siz2,k+siz2) || c(i+siz2,j+siz2,k+siz2));
# end
# end
# end
# n(g+1) = sum(sum(sum(c(1:siz:(width-siz+1),1:siz:(width-siz+1),1:siz:(width-siz+1)))));
# end
#
# end

#' 2D Boxcount for 1D signal
#'
#' @param y A numeric vector or time series object.
#' @param unitSquare Create unit square image of `y`? This is required for estimating FD of time series (default = `TRUE`)
#' @param image2D A matrix representing a 2D image, argument `y` and `unitSquare` will be ignored (default = `NA`)
#' @param resolution The resolution used to embed the timeseries in 2D, a factor by which the dimensions the matrix will be multiplied (default = `1`)
#' @param removeTrend If `TRUE`, will call \link[casnet]{ts_detrend} on `y` (default = `FALSE`)
#' @param polyOrder Order of polynomial trend to remove if `removeTrend = `TRUE``
#' @param standardise Standardise `y` using [casnet::ts_standardise()] with `adjustN = FALSE` (default = `none`)
#' @param adjustSumOrder Adjust the order of the time series (by summation or differencing), based on the global scaling exponent, see e.g. <https://www.frontiersin.org/files/Articles/23948/fphys-03-00141-r2/image_m/fphys-03-00141-t001.jpg>{Ihlen (2012)} (default = `FALSE``)
#' @param scaleMax Maximum scale value (as `2^scale`) to use (default = `max` of `log2(nrows)` and `log2(ncols)`)
#' @param scaleMin Minimium scale value (as `2^scale`) to use (default = `0`)
#' @param scaleS If not `NA`, pass a numeric vector listing the scales (as a power of `2`) on which to evaluate the boxcount. Arguments `scaleMax`, `scaleMin`, and `scaleResolution` will be ignored (default = `NA`)
#' @param dataMin Minimum number of time/data points inside a box for it to be included in the slope estimation (default = `2^scaleMin`)
#' @param maxData Maximum number of time/data points inside a box for it to be included in the slope estimation (default = `2^scaleMax`)
#' @param doPlot Return the log-log scale versus bulk plot with linear fit (default = `TRUE`).
#' @param returnPlot Return ggplot2 object (default = `FALSE`)
#' @param returnPLAW Return the power law data (default = `FALSE`)
#' @param returnInfo Return all the data used in DFA (default = `FALSE`)
#' @param returnLocalScaling Return estimates of FD for each scale
#' @param silent Silent-ish mode (default = `TRUE`)
#' @param noTitle Do not generate a title (only the subtitle)
#' @param tsName Name of y added as a subtitle to the plot (default = `y`)
#'
#' @return The boxcount fractal dimension and the 'local' boxcount fractal dimension
#'
#' @note This function was inspired by the `Matlab` function `boxcount.m` [written by F. Moisy](http://www.fast.u-psud.fr/~moisy/ml/boxcount/html/demo.html). `Fred Hasselman` adapted the function for `R` for the purpose of the unit square boxcount analysis for 1D time series. The original Matlab toolbox has more options and contains more functions (e.g. `1D` and `3D` boxcount).
#'
#' @export
#'
#' @examples
#'
#' \donttest{fd_boxcount2D(y = rnorm(100))}
#'
#'
fd_boxcount2D <- function(y = NA,
                          unitSquare = TRUE,
                          image2D = NA,
                          resolution = 1,
                          removeTrend = FALSE,
                          polyOrder = 1,
                          standardise = c("none","mean.sd","median.mad")[1],
                          adjustSumOrder = FALSE,
                          scaleMin = 0,
                          scaleMax = floor(log2(NROW(y)*resolution)),
                          scaleS = NA,
                          dataMin = 2^(scaleMin+1),
                          maxData = 2^(scaleMax-1),
                          doPlot = FALSE,
                          returnPlot = FALSE,
                          returnPLAW = FALSE,
                          returnInfo = FALSE,
                          returnLocalScaling = FALSE,
                          silent = FALSE,
                          noTitle = FALSE,
                          tsName="y"
){

  # scaleResolution The range of scales at which the boxcount will be evaluated are calculatd as: `(scaleMax-scaleMin)/scaleResolution`.  (default = `(scaleMax-scaleMin)`)

  scaleResolution <- (scaleMax-scaleMin)

  if(is.na(image2D)){

    if(removeTrend){
      y <- ts_detrend(y, polyOrder = polyOrder)
    }

    if(standardise%in%(c("mean.sd","median.mad"))){
      y <- ts_standardise(y, type=standardise)
    }

    if(adjustSumOrder){
      y       <- ts_sumorder(y, scaleS = 2^(4:floor(log2(NROW(y)/2))), polyOrder = polyOrder, dataMin = 4)
      Hadj    <- attr(y,"Hadj")
      Hglobal <- attr(y,"Hglobal.excl")
    } else {
      Hadj    <- 0
      Hglobal <- NA
    }

    cat("\n\nRaterizing time series... ")
    image2D <- ts_rasterize(y = y, unitSquare = unitSquare, toSparse = TRUE)

  } else {

    unitSquare <- FALSE

  }

  cat("Done!\n")
  # Just to be sure we have a 2D matrix
  image2D <- methods::as(drop(image2D), 'lgCMatrix')
  N <- dim(image2D)
  if(length(N)!=2){
    stop('Matrix dimension must be 2!')
  }

  if(unitSquare){
    if(identical(N[1],N[2])){
      N <- N[1]
    } else {
      stop("Need a square matrix, because unitSquare = TRUE")
    }
  } else {
    N <- max(N)
  }

  if(!any(is.na(scaleS))&is.numeric(scaleS)){
    scaleS   <- unique(sort(scaleS))
    scaleMin <- min(scaleS, na.rm = TRUE)
    if(scaleMin < 0){
      scaleMin  <- 0
      scaleS[1] <- 0
    }
    scaleMax <- max(scaleS, na.rm = TRUE)
  } else {
    scaleS <- unique(round(2^(seq(scaleMin, scaleMax, by=((scaleMax-scaleMin)/scaleResolution)))))
  }

  Nscales  <- length(scaleS)

  if(scaleMax > stats::nextn(N,2)){warning("scaleMax should not be larger than the next power of 2 of largest dimension!")}

  # Need padding?
  p <- log(N)/log(2)
  if(p!=round(p) || any(dim(image2D)!=2^scaleMax)){
    pu     <- ceiling(p)
    width  <- 2^pu
    mz    <- pracma::zeros(width, width)
    mz[1:dim(image2D)[1], 1:dim(image2D)[2]] <- Matrix::as.matrix(image2D)
    image2D  <- methods::as(mz, 'lgCMatrix')
    scaleMax <- pu
    N        <- max(dim(image2D))
    scaleS   <- unique(sort(c(scaleS, 2^floor(p), 2^scaleMax)))
    Nscales  <- length(scaleS)
    rm(mz,p,pu)
  }

  #unique(round(2^(seqscaleMax-scaleMin)/scaleResolution rev(seq(0,p-1,1))

  # if(!all(is.numeric(scaleS),length(scaleS)>0,scaleS%[]%c(2,(NROW(y)/2)))){
  #   message("Something wrong with vector passed to scaleS.... \nUsing defaults: (scaleMax-scaleMin)/scaleResolution")
  # }

  #image2D <- matrix(as.logical(image2D),ncol=N,nrow=N)
  boxeS   <- unique(rev(signif(log2(scaleS[1:Nscales-1]),2)))
  Nscales <- length(boxeS)

  # 2D boxcount

  cat("Performing 2D boxcount...")

  # pre-allocate the number of boxes of size r
  n            <- numeric(Nscales+1)
  n[Nscales+1] <- sum(image2D, na.rm = TRUE)
  width        <- 2^scaleMax

  for(g in boxeS){
    siz  <- 2^(scaleMax-g)
    siz2 <- round(siz/2)
    for(i in seq(1,(width-siz+1),siz)){
      for(j in seq(1,(width-siz+1),siz)){
        image2D[i,j] <- ( image2D[i,j] || image2D[i+siz2,j] || image2D[i,j+siz2] || image2D[i+siz2,j+siz2])
      }
    }
    n[g+1] <- sum(sum(image2D[seq(1,(width-siz+1),siz),seq(1,(width-siz+1),siz)]))
  }

  n <- rev(n)
  r <- scaleS

  lmfit1 <- stats::lm(-log(n)~log(r)) #lmfit1 <- mean(-diff(log(n))/diff(log(r)))

  fitRange <- which(r%[]%c(dataMin,maxData))
  lmfit2   <- stats::lm(-log(n[fitRange])~log(r[fitRange])) #mean(-diff(log(n[fitRange]))/diff(log(r[fitRange])))

  localFit <- -pracma::gradient(log(n))/pracma::gradient(log(r))

  cat("Done!\n")

  outList <- list(
    PLAW  =  data.frame(size = r, bulk = n),
    fullRange = list(y = y,
                     sap = lmfit1$coefficients[2],
                     Hadj = Hadj,
                     H = NA,
                     FD = lmfit1$coefficients[2],
                     fitlm1 = lmfit1,
                     method = paste0("Full range (n = ",NROW(r),")\nFD = ",signif(lmfit1$coefficients[2],3)),
                     fitRange = r),
    fitRange  = list(y = y,
                     sap = lmfit2$coefficients[2],
                     H = NA,
                     Hadj = Hadj,
                     FD = lmfit2$coefficients[2],
                     fitlm2 = lmfit2,
                     method = paste0("Exclude scales (n = ",NROW(fitRange),")\nFD = ",signif(lmfit2$coefficients[2],3)),
                     fitRange = r[fitRange]),
    info = list(fullRange=lmfit1,fitRange=lmfit2, localFit = localFit),
    plot = NA,
    analysis = list(
      name = "2D boxcount of 1D curve",
      logBaseFit = "log",
      logBasePlot = "2")
  )


  if(doPlot|returnPlot){
    if(noTitle){
      title <- " "
    } else {
      title <- "log-log regression (2D Boxcount)"
    }
    if(doPlot){
      g <- plotFD_loglog(fd.OUT = outList, title = title, subtitle = tsName, logBase = "2",xlabel = "Box Size (log)", ylabel = "Box Count (-log)", doPlot = doPlot)

      if(returnLocalScaling){

        logFormat <- "log2"
        logBaseNum <- 2

        yAcc <- .01
        xAcc <- 1

        breaksX <- unique(outList$PLAW$size[1:NROW(outList[[2]]$fitlm1$fitted.values)])
        breaksY <- unique(localFit)

        if(length(breaksX)>10){breaksX <- breaksX[seq.int(1,length(breaksX),length.out = 10)]}
        if(length(breaksY)>10){breaksY <- breaksY[seq.int(1,length(breaksY),length.out = 10)]}

        gl <- ggplot2::ggplot(data.frame(size=r, bulk=localFit), ggplot2::aes_(x=~size,y=~bulk), na.rm=TRUE) +
          ggplot2::geom_point(colour = "red3") +
          ggplot2::geom_line(colour = "steelblue") +
          ggplot2::ggtitle(label = title, subtitle = "Local Dimension") +
          ggplot2::scale_x_continuous(name = paste0("Box Size (",logFormat,")"),
                                      breaks = breaksX,
                                      labels = scales::number_format(accuracy = xAcc),
                                      trans = scales::log_trans(base = logBaseNum)) +
          ggplot2::scale_y_continuous(name = "Local Dimension [ - d log(count) / d log(size) ]",
                                      breaks = breaksY,
                                      labels = scales::number_format(accuracy = yAcc),
                                      limits = c(1,2)) +
          ggplot2::scale_color_manual(values = c("red3","steelblue")) +
          ggplot2::theme_bw() +
          ggplot2::theme(panel.grid.minor =  ggplot2::element_blank(),
                         legend.text = ggplot2::element_text(margin = ggplot2::margin(t = 5,b = 5, unit = "pt"), vjust = .5),
                         plot.margin = ggplot2::margin(t = 5,b = 5, r = 5,l = 5, unit = "pt"))

        graphics::plot(gl)

      }
      if(returnPlot){
        if(returnLocalScaling){
          outList$plot <- list(global = g, local = gl)
        } else {
          outList$plot <- g
        }
      }
    }
  }

  if(returnInfo){returnPLAW <- TRUE}

  if(!silent){
    cat("\n~~~o~~o~~casnet~~o~~o~~~\n")
    cat(paste("\n",outList$analysis$name,"\n\n",outList$fullRange$method,"\n\n",outList$fitRange$method))
    cat("\n\n~~~o~~o~~casnet~~o~~o~~~\n")
  }

  return(invisible(outList[c(returnPLAW,TRUE,TRUE,returnInfo,returnPlot,TRUE)]))


}


#' Calculate FD using Sevcik's method
#'
#' @param y A time series or numeric vector
#' @param detrend Subtract linear trend from the series (default = `TRUE`).
#' @param adjustSumOrder Adjust the time series (summation or differencing), based on the global scaling exponent, see e.g. <https://www.frontiersin.org/files/Articles/23948/fphys-03-00141-r2/image_m/fphys-03-00141-t001.jpg>{Ihlen (2012)} (default = `TRUE`)
#' @param smallNapprox Force use of small sample approximation (default for N < 128)
#' @param doPlot   Return the log-log scale versus fluctuation plot with linear fit (default = `TRUE`).
#' @param returnPlot Return ggplot2 object (default = `FALSE`)
#' @param returnPLAW Return the power law data (default = `FALSE`)
#' @param returnInfo Return all the data used in DFA (default = `FALSE`)
#' @param silent Silent-ish mode
#' @param noTitle Do not generate a title (only the subtitle)
#' @param tsName Name of y added as a subtitle to the plot
#'
#' @author Fred Hasselman
#'
#' @return An FD estimate
#'
#' @export
#'
#' @family Fluctuation Analyses
#'
#' @references Sevcik, C. (1998). A procedure to Estimate the Fractal Dimension of Waveforms. Paper available at http://arxiv.org/pdf/1003.5266.pdf
#'
fd_sev <- function(y,
                   detrend = FALSE,
                   adjustSumOrder = FALSE,
                   smallNapprox = FALSE,
                   doPlot = FALSE,
                   returnPlot = FALSE,
                   returnPLAW = FALSE,
                   returnInfo = FALSE,
                   silent = FALSE,
                   noTitle = FALSE,
                   tsName="y"){

  Hadj<-0
  if(detrend){y <- ts_detrend(y)}
  if(adjustSumOrder){
    y <- ts_sumorder(y)
    Hadj <- attr(y,"Hadj")
  }

  N <- NROW(y)

  D <- data.frame(L = numeric(N),FD = numeric(N),var = numeric(N), n = numeric(N))

  # This is an implementation of Sevcik's program, which I  optimised for speed...
  #
  # ymin <- y[1]
  # ymax <- y[1]
  #
  # L <- 0
  # for(n in 2:N){
  #   if(y[n]>ymax){ymax <- y[n]}
  #   if(y[n]<ymin){ymin <- y[n]}
  #
  #   if(n>1){
  #     L <- 0
  #     for(i in 1:n){
  #       yy = (y[i] - ymin) / (ymax - ymin)
  #       if(i>1){
  #         L <- L + sqrt((yy -yant)^2 + (1 /(n-1))^2)
  #       }
  #       yant <- y
  #     }
  #   }
  #
  #   #L <- cumsum(sqrt(diff(elascer(y[1:n]))^2 + (1/(n-1)^2)))
  #   #l <- dplyr::last(L)
  #   if(L==0){L <- L +.Machine$double.eps}
  #   if(N<128|smallNapprox){
  #     D$FD[n] <- 1 + ((log(L) - log(2)) / log(2*(n-1)))
  #   } else {
  #     D$FD[n] <- 1 + (log(L) / log(2*(n-1)))
  #   }
  #  D$var[n] <- (stats::var(diff(y[1:n]), na.rm = TRUE) * (n-1)) / (L^2 * log(2*(n-1)^2))
  #  D$L[n] <- L
  #  D$n[n] <- n
  #  rm(L)
  # }

  # Also store each cumulative sum profile
  D.L <- list()

  for(n in 2:N){
    L <- cumsum(sqrt(diff(elascer(y[1:n]))^2 + (1/(n-1)^2)))
    l <- dplyr::last(L)
    if(l==0){l <- l +.Machine$double.eps}
    if(N<128|smallNapprox){
      D$FD[n] <- 1 + ((log(l) - log(2)) / log(2*(n-1)))
    } else {
      D$FD[n] <- 1 + (log(l) / log(2*(n-1)))
    }
    D$var[n] <- (stats::var(diff(y[1:n]), na.rm = TRUE) * (n-1)) / (l^2 * log(2*(n-1)^2))
    D$L[n] <- l
    D$n[n] <- n
    D.L[[n]] <- L
    rm(L)
  }

  D.sd <- ts_sd(D$FD, type = "unadjusted")
  D$SE <- sqrt(D$var/(D$n-1))

  #D.H <- lm(D$FD[D$n>=N*.25]~log(D$n[D$n>=N*.25]))

  fitlm1 <- fitlm2 <- list()
  fitlm1$fitted.values <- D$n
  fitlm2$fitted.values <- D$n[1:round(N/2)]

  outList <- list(PLAW = data.frame(bulk = D$FD, size = D$n),
                  fullRange  = list(sap = NA, H = NA,
                                    FD = dplyr::last(D$FD),
                                    fitlm1 = fitlm1,
                                    method = paste0("Lengths (n = ",N),")\n FD = ", round(dplyr::last(D$FD),2)),
                  fitRange   = list(sap = NA, H = NA, FD = dplyr::last(D$FD[round(D$FD/2)]),
                                    fitlm2 = fitlm2,
                                    method = paste0("Lengths (n = ",round(N/2),")\n FD = ", round(D$FD[round(N/2)],2))),
                  plot = NA,
                  info = list(D = D, L = D.L),
                  analysis = list(
                    name = "Planar Extent - Sevcik's method",
                    logBaseFit = "no",
                    logBasePlot = "no")
  )


  # if(doPlot){

  # dfU <- data.frame(time=seq(0,1,length.out = N),y=elascer(y))
  # ggU <- ggplot(dfU,aes_(x=time,y=y)) + geom_line() + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + coord_equal() + theme_bw() + theme(panel.grid = element_blank())
  #
  # dfOri <- data.frame(time=time(ts(y)),y=y)
  # ggOri <- ggplot(dfOri,aes_(x=time,y=y)) + geom_line() + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + theme_bw() + theme(panel.grid = element_blank())
  #
  # dfFD <- data.frame(time=time(ts(y)),FD=D.FD,ci_lo=D.FD-(1.96*D.SE),ci_hi=D.FD+(1.96*D.SE))
  # ggFD <- ggplot(dfFD,aes_(x=time,y=FD)) +  geom_ribbon(aes_(ymin=ci_lo, ymax=ci_hi),colour="grey70", fill = "grey70")+ geom_line() +  scale_x_continuous(expand = c(0,0)) + scale_y_continuous("Fractal Dimension",breaks = c(0.8,1,1.1,1.2,1.5,1.8),expand = c(0,0), limits = c(.8,2)) + theme_bw() + theme(panel.grid.major.x  = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())

  # # graphics::plot.new()
  # old <- ifultools::splitplot(2,1,1)
  # #graphics::plot(y,ylab = "Y", main = paste0('FD: ', round(dplyr::last(D.FD), digits=2)))
  # graphics::plot(elascer(y)~ seq(0,1,length.out = N), ylab = "", xlab="", main = paste0('FD: ', round(dplyr::last(D.FD), digits=2)),xlim=c(0,1),ylim=c(0,1),pty="s", type="l")
  # ifultools::splitplot(2,1,2)
  # graphics::plot(D.FD ~ seq(0,1,length.out = N), xlab="Normalised time", ylab = "", type="l")
  # #graphics::legend("topleft",c(paste0("Full (n = ",length(DFAout$PLAW$scaleS),")"), paste0("Range (n = ",length(DFAout$PLAW$scaleS[fitRange]),")")), lwd=c(3,3),col=c("darkred","darkblue"), cex = .8)
  # graphics::par(old)


  if(doPlot|returnPlot){
    if(noTitle){
      title <- " "
    } else {
      title <- "Sevcik's method (planar extent)"
    }
    g <- plotFD_loglog(fd.OUT = outList, title = title, subtitle = tsName, logBase = "no",xlabel = "Time series length", ylabel = "Fractal Dimension")
    if(doPlot){
      grid::grid.newpage()
      grid::grid.draw(g)
    }
    if(returnPlot){
      outList$plot <- g
    }
  }

  if(returnInfo){returnPLAW<-TRUE}

  if(!silent){
    cat("\n~~~o~~o~~casnet~~o~~o~~~\n")
    cat(paste("\n",outList$analysis$name,"\n\n",outList$fullRange$method,"\n\n",outList$fitRange$method))
    cat("\n\n~~~o~~o~~casnet~~o~~o~~~\n")
  }

  return(invisible(outList[c(returnPLAW,TRUE,TRUE,returnInfo,returnPlot,TRUE)]))
}


#' Allan Variance Analysis
#'
#' @param y A numeric vector or time series object
#' @param fs Sample frequency in Hz
#' @param useSD Use the standarddeviation instead of variance?
#' @param doPlot   Return the log-log scale versus fluctuation plot with linear fit (default = `TRUE`).
#' @param returnPlot Return ggplot2 object (default = `FALSE`)
#' @param returnPLAW Return the power law data (default = `FALSE`)
#' @param returnInfo Return all the data used in DFA (default = `FALSE`)
#' @param silent Silent-ish mode
#' @param noTitle Do not generate a title (only the subtitle)
#' @param tsName Name of y added as a subtitle to the plot
#'
#' @return A dataframe with the Allan Factor (variance), Alan standard deviation and error due to bin size
#' @export
#' @family Fluctuation Analyses
#'
fd_allan <- function(y,
                     fs = stats::tsp(stats::hasTsp(y))[3],
                     useSD=FALSE,
                     doPlot = FALSE,
                     returnPlot = FALSE,
                     returnPLAW = FALSE,
                     returnInfo = FALSE,
                     silent = FALSE,
                     noTitle = FALSE,
                     tsName="y"){

  N <- NROW(y)
  tau <- 1/fs
  n <- ceiling((N - 1)/2)
  p <- floor(log10(n)/log10(2))
  av <- rep(0, p + 1)
  time <- rep(0, p + 1)
  error <- rep(0, p + 1)
  for (i in 0:(p)) {
    omega <- rep(0, floor(N/(2^i)))
    TT <- (2^i) * tau
    l <- 1
    k <- 1
    while (k <= floor(N/(2^i))) {
      omega[k] <- sum(y[l:(l + ((2^i) - 1))])/(2^i)
      l <- l + (2^i)
      k <- k + 1
    }
    sumvalue <- 0
    for (k in 1:(length(omega) - 1)) {
      sumvalue <- sumvalue + (omega[k + 1] - omega[k])^2
    }
    av[i + 1] <- sumvalue/(2 * (length(omega) - 1))
    time[i + 1] <- TT
    error[i + 1] <- 1/sqrt(2 * ((N/(2^i)) - 1))
  }

  if(useSD){
    df <- data.frame(bulk = time, size = av, error=error)
  } else {
    df <- data.frame(bulk = time, size = sqrt(av), error=error)
  }

  lmfit1        <- stats::lm(log10(df$bulk) ~ log10(df$size), na.action=stats::na.omit)
  H1  <- lmfit1$coefficients[2]
  lmfit2        <- stats::lm(log10(df$bulk) ~ log10(df$size), na.action=stats::na.omit)
  H2  <- lmfit2$coefficients[2]



  outList <- list(
    PLAW  =  df,
    fullRange = list(y = y,
                     sap = lmfit1$coefficients[2],
                     AF = H1,
                     FD = NA,
                     fitlm1 = lmfit1,
                     method = paste0("Full range (n = ",NROW(df),")\nSlope = ",round(stats::coef(lmfit1)[2],2))),
    fitRange  = list(y = y,
                     sap = lmfit2$coefficients[2],
                     AF = H2,
                     FD = NA,
                     fitlm2 = lmfit2,
                     method = paste0("Full range (n = ",NROW(df),")\nSlope = ",round(stats::coef(lmfit2)[2],2))),
    info = df,
    plot = NA,
    analysis = list(
      name = "Allan Variance/Deviation",
      logBaseFit = "log10",
      logBasePlot = 10)
  )

  if(doPlot|returnPlot){
    if(noTitle){
      title <- " "
    } else {
      if(useSD){
        title <- "log-log regression (Allan Deviation)"
      } else {
        title <- "log-log regression (Allan Variance)"
      }
    }
    g <- plotFD_loglog(fd.OUT = outList, title = title, subtitle = tsName, logBase = "10",xlabel = "Scale", ylabel = "Allan Variance", doPlot = doPlot)
    if(doPlot){
      grid::grid.newpage()
      grid::grid.draw(g)
    }
    if(returnPlot){
      outList$plot <- g
    }
  }

  if(returnInfo){returnPLAW<-TRUE}

  if(!silent){
    cat("\n~~~o~~o~~casnet~~o~~o~~~\n")
    cat(paste("\n",outList$analysis$name,"\n\n",outList$fullRange$method))
    cat("\n\n~~~o~~o~~casnet~~o~~o~~~\n")
  }
  return(invisible(outList[c(returnPLAW,TRUE,TRUE,returnInfo,returnPlot,TRUE)]))

  #
  # if(doPlot){
  #   if(useSD){
  #
  #  g <-  ggplot(df,aes_(x=~Tcluster,y=~ATsd)) +
  #     geom_path() +
  #     geom_pointrange(aes_(ymin=~ATsd-(~ATsd*~error),ymax=~ATsd+(~ATsd*~error))) +
  #    scale_y_continuous(trans = scales::log10_trans(),
  #                       breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                       labels = scales::trans_format("log10", scales::math_format())
  #                       ) +
  #    scale_x_continuous(trans = scales::log10_trans(),
  #                       breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                       labels = scales::trans_format("log10", scales::math_format())) +
  #     xlab("Cluster Times (T)") +
  #     ylab("Allan Standard Deviation") +
  #     theme_bw()
  #
  #  } else {
  #    g <-  ggplot(df,aes_(x=~Tcluster,y=~AT)) +
  #      geom_path() +
  #      geom_pointrange(aes_(ymin=~AT-(~AT*~error),ymax=~AT+(~AT*~error))) +
  #    scale_y_continuous(trans = scales::log10_trans(),
  #                       breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                       labels = scales::trans_format("log10", scales::math_format())) +
  #      scale_x_continuous(trans = scales::log10_trans(),
  #                         breaks = scales::trans_breaks("log10", function(x) 10^x),
  #                         labels = scales::trans_format("log10", scales::math_format())) +
  #      xlab("Cluster Times (T)") +
  #      ylab("Allan Factor")  +
  #      theme_bw()
  #  }
  #   graphics::plot(g)
  # }
  # return(df)
}


#' Multi-fractal Detrended Fluctuation Analysis
#'
#' @inheritParams fd_dfa
#' @param qq A vector containing a range of values for the order of fluctuation `q` (default = `seq(-5, 5,length.out=101)`)
#'
#' @return A dataframe with values of `q`,`H(q)`, `t(q)`, `h(q)`, `D(q)``
#' @export
#'
#' @family Fluctuation Analyses
#'
#' @examples
#'
#' set.seed(33)
#'
#' # White noise
#' fd_mfdfa(rnorm(4096), doPlot = TRUE)
#'
#' # Pink noise
#' fd_mfdfa(noise_powerlaw(N=4096), doPlot = TRUE)
#'
#' # 'multi' fractal
#' N <- 2048
#' y <- rowSums(data.frame(elascer(noise_powerlaw(N=N, alpha = -2)), elascer(noise_powerlaw(N=N, alpha = -.5))*c(rep(.2,512),rep(.5,512),rep(.7,512),rep(1,512))))
#' fd_mfdfa(y=y, doPlot = TRUE)
#'
fd_mfdfa <- function(y,
                     fs = NULL,
                     removeTrend = c("no","poly","adaptive","bridge")[2],
                     polyOrder=1,
                     standardise = c("none","mean.sd","median.mad")[1],
                     adjustSumOrder = FALSE,
                     removeTrendSegment = c("no","poly","adaptive","bridge")[2],
                     polyOrderSegment=1,
                     scaleMin = 16,
                     scaleMax = stats::nextn(floor(NROW(y)/4), factors = 2),
                     scaleResolution = round(log2(scaleMax-scaleMin)),
                     dataMin = NA,
                     scaleS = NA,
                     overlap = NA,
                     qq = seq(-5, 5,length.out=101),
                     doPlot = FALSE,
                     returnPlot = FALSE,
                     returnInfo = FALSE,
                     silent = FALSE){

  yy <- fd_prepSeries(y = y,
                      fs = fs,
                      removeTrend = "no",
                      polyOrder = polyOrder,
                      standardise = standardise,
                      adjustSumOrder = adjustSumOrder,
                      scaleMin = scaleMin,
                      scaleMax = scaleMax,
                      scaleResolution = scaleResolution,
                      dataMin = dataMin,
                      scaleS = scaleS,
                      overlap = overlap,
                      silent = silent)

  # Integrate the series
  if(standardise%in%"none"){
    y <- cumsum(ts_center(yy)) # need negative values in profile
  } else {
    y <- cumsum(yy)
  }
  y <- rp_copy_attributes(source = yy, target = y)

  scaleS  <- attr(y, "scaleS")
  dataMin <- attr(y, "dataMin")

  qq <- c(qq,(qq[length(qq)]+.1))

  scale     <- scaleS #round(2^(seq(scaleMin,scaleMax,by=((scaleMax-scaleMin)/scaleResolution))))
  segv      <- numeric(length(scale))
  RMS_scale <- vector("list",length(scale))
  qRMSns    <- vector("list",length(scaleS))
  Fq        <- matrix(nrow = length(qq), ncol = length(scaleS))
  qRegLine  <- vector("list",length(qq))
  Hq        <- numeric(length(qq))

  TSm      <- as.matrix(cbind(t=1:NROW(y),y=y))

  if(removeTrendSegment%in%"adaptive"){
    adaptive <- TRUE
  } else {
    adaptive <- FALSE
  }
  for(ns in seq_along(scale)){
    RMS_scale[[ns]] <- plyr::ldply(ts_slice(y = TSm[,2], epochSz = scale[ns], overlap = overlap, removeUnequal = TRUE), function(sv){return(sqrt(mean(ts_detrend(sv,polyOrder = polyOrder, adaptive = adaptive)^2, na.rm = TRUE)))})
    qRMS      <- vector("list",length(qq))
    for(nq in seq_along(qq)){
      qRMS[[nq]] <- RMS_scale[[ns]]$V1^qq[nq]
      Fq[nq,ns] <- mean(qRMS[[nq]])^(1/qq[nq])
      #if(is.infinite(log2(Fq[nq,ns]))){Fq[nq,ns]<-NA
    }
    qRMSns[[ns]] <- qRMS[[nq]]
    rm(qRMS)
    Fq[which(qq==0),ns] <- exp(0.5*mean(log(RMS_scale[[ns]]$V1^2)))%00%NA
    #if(is.infinite(log2(Fq[[which(qq==0)]][ns]))){Fq[[which(qq==0)]][ns]<-NA}
  }

  # fmin<-1
  # fmax<-which(scale==max(scale))
  #for(nq in seq_along(qq)){Hq[nq] <- stats::lm(log2(Fq[[nq]])~log2(scale))$coefficients[2]}
  Hq <- plyr::ldply(seq_along(qq), function(nq){stats::lm(log2(Fq[nq,])~log2(scale),na.action=stats::na.omit)$coefficients[2]})

  tq <- Hq[,1]*qq-1
  hq <- diff(c(tq,tq[length(tq)]))/diff(c(qq,qq[length(qq)])) #[-c(1,length(tq))]
  Dq <- (qq*hq) - tq

  tq <- tq[1:(length(qq)-1)]
  hq <- hq[1:(length(qq)-1)]
  Dq <- Dq[1:(length(qq)-1)]
  Hq <- Hq[1:(length(qq)-1),1]
  qq <- qq[1:(length(qq)-1)]

  # Dq <- (qq[1:(length(qq)-1)]*hq) - tq[1:(length(qq)-1)]
  # Dh=(qq*hq)-tq
  #if(reload==TRUE){library(y,verbose=FALSE,quietly=TRUE)}

  # qRegLine <- as.data.frame(Reduce("cbind", qRegLine))
  # tq <- Hq * q - 1
  # hq <- diff(tq)/(q[2] - q[1])
  # Dq <- (q[1:(length(q) - 1)] * hq) - tq[1:(length(tq) - 1)]

  # tau=(q.*Hq)-1;
  #
  # hh=diff(tau)./(q(2)-q(1));
  # Dh=(q(1:(end-1)).*hh)-tau(1:(end-1));
  # h=hh-1;
  #plot(hq,Dq, type="l")

  out <- data.frame(q=qq,Hq=Hq,tq=tq,hq=hq,Dq=Dq)

  id <- order(out$hq)
  Spec_AUC <- sum(diff(out$hq[id])*zoo::rollmean(out$Dq[id],2), na.rm = TRUE)

  Spec_Width  <- max(out$hq, na.rm = TRUE) - min(out$hq, na.rm = TRUE)
  Spec_CVplus <- sd(out$Dq[which(out$q==0):which.max(out$q)], na.rm = TRUE)/mean(out$Dq[which(out$q==0):which.max(out$q)], na.rm = TRUE)
  Spec_CVmin  <- sd(out$Dq[which.min(out$q):which(out$q==0)], na.rm = TRUE)/mean(out$Dq[which.min(out$q):which(out$q==0)], na.rm = TRUE)
  Spec_CVtot  <- sd(out$Dq, na.rm = TRUE)/mean(out$Dq, na.rm = TRUE)
  Spec_CVasymm <- (Spec_CVplus-Spec_CVmin)/(Spec_CVplus+Spec_CVmin)


  outList <- list(MFDFA = out,
                  MFSpectrum = data.frame(Spec_AUC = Spec_AUC,
                                          Spec_Width = Spec_Width,
                                          Spec_CVplus = Spec_CVplus,
                                          Spec_CVmin = Spec_CVmin,
                                          Spec_CVtot = Spec_CVtot,
                                          Spec_CVasymm = Spec_CVasymm),
                  ScalingFunction = Fq
                  )

  # if(doPlot|returnPlot){
  #   if(noTitle){
  #     title <- " "
  #   } else {
  #     title <- "log-log regression (DFA)"
  #   }
  #   g <- plotFD_loglog(fd.OUT = outList, title = title, subtitle = tsName, logBase = "2",xlabel = "Scale", ylabel = "Detrended Fluctuation", doPlot = doPlot)
  #   if(doPlot){
  #     print(g)
  #   }
  #   if(returnPlot){
  #     outList$plot <- g
  #   }
  # }

  if(doPlot|returnPlot){

    qmin  <- which(out$q==min(out$q,na.rm = TRUE))
    qzero <- which(out$q==0)
    qmax  <- which(out$q==max(out$q,na.rm = TRUE))

    df_q <- data.frame(q = out$q[c(qmin, qzero,qmax)],
                       label = paste0("q=",out$q[c(qmin, qzero,qmax)]),
                       labelH = paste0("H(q=",out$q[c(qmin, qzero,qmax)],") = ",round(out$Hq[c(qmin, qzero,qmax)],2)),
                       hq = out$hq[c(qmin, qzero,qmax)],
                       Hq = out$Hq[c(qmin, qzero,qmax)],
                       tq = out$tq[c(qmin, qzero,qmax)],
                       Dq = out$Dq[c(qmin, qzero,qmax)])

    df_lines <- data.frame(Scale_label = paste(c(scale,scale,scale)),
                           Scale = c(scale,scale,scale),
                           labelH = rep(df_q$labelH,each = length(scale)),
                           Fq    = log2(c(Fq[qmin,],Fq[qzero,],Fq[qmax,])),
                           q     = rep(out$q[c(qmin, qzero,qmax)], each = length(scale)))
    df_lines$label <- factor(paste0("q=",df_lines$q))

    cols <- getColours(3)
    names(cols) <- levels(df_lines$labelH)

    g1 <-  ggplot(data = df_lines, aes(x=Scale,y= Fq, group = q)) +
      geom_smooth(aes(colour = labelH), method = "lm", se = FALSE) +
      geom_point(aes(fill = labelH), pch=21,colour="black") +
      scale_x_continuous("Scale", limits = c(min(df_lines$Scale, na.rm = TRUE),
                                             max(df_lines$Scale, na.rm = TRUE)), trans = "log2") +
      scale_y_continuous(expression(log[2](F[q])), limits = c(min(df_lines$Fq, na.rm = TRUE)-2,
                                                              max(df_lines$Fq, na.rm = TRUE)+2), expand = c(0,0)) +
      scale_color_manual("Slope", values = cols) +
      scale_fill_manual("Slope", values = cols) +
      ggtitle("Scaling function (Slope = Hurst exponent)") +
      theme_bw()  +
      theme(legend.position = c(.85,.25),
            legend.text = element_text(size=rel(.6)),
            legend.title = element_text(size=rel(.6)),
            legend.key.size = unit(.8,"lines"),
            legend.background = element_rect(colour = "black"))

    g2 <- ggplot(out,aes(x=q,y=Hq)) + geom_line() +
      geom_point(data = df_q, aes(x=q,y=Hq, fill = labelH), pch=21,colour="black", show.legend = FALSE, size = 2) +
      scale_x_continuous("q", limits = c(min(out$q, na.rm = TRUE)-.2,
                                         max(out$q, na.rm = TRUE)+.2), expand = c(0,0)) +
      scale_y_continuous(expression(H[q]), limits = c(min(out$Hq, na.rm = TRUE)-.1,
                                                      max(out$Hq, na.rm = TRUE)+.1), expand = c(0,0)) +
      #scale_color_discrete("q-order") +
      scale_fill_manual("Slope", values = cols) +
      theme_bw()  + ggtitle("q-order Hurst exponent")


    g3 <- ggplot(out,aes(x=q,y=tq)) +
      geom_line() +
      #      geom_point(data = df_q, aes(x=q,y=tq, colour = label),show.legend = FALSE) +
      geom_point(data = df_q, aes(x=q,y=tq, fill = label), pch=21,colour="black", show.legend = FALSE, size = 2) +
      scale_x_continuous("q", limits = c(min(out$q, na.rm = TRUE)-.2,
                                         max(out$q, na.rm = TRUE)+.2), expand = c(0,0)) +
      scale_y_continuous(expression(t[q]), limits = c(min(out$tq, na.rm = TRUE)-.2,
                                                      max(out$tq, na.rm = TRUE)+.2), expand = c(0,0)) +
      scale_fill_manual("Slope", values = cols) +
      theme_bw()  + ggtitle("q-order Mass exponent")


    g4 <- ggplot(out,aes(x=hq,y=Dq)) + geom_line() +
      geom_point(data = df_q, aes(x=hq,y=Dq, fill = label), pch=21, colour="black", show.legend = TRUE, size = 2) +
      scale_x_continuous(expression(h[q]), limits = c(min(out$hq, na.rm = TRUE)-.1,
                                                      max(out$hq, na.rm = TRUE)+.1), expand = c(0,0)) +
      scale_y_continuous(expression(D[q]), limits = c(min(c(out$Dq,0), na.rm = TRUE)-.05,
                                                      max(c(out$Dq,1), na.rm = TRUE)+.05), expand = c(0,0)) +
      scale_color_discrete("q-order") +
      scale_fill_manual("q-order", values = cols) +
      ggtitle("Multifractal spectrum") +
      theme_bw()  +
      theme(legend.position = c(.4,.5),
            legend.text = element_text(size=rel(.6)),
            legend.title = element_text(size=rel(.6)),
            legend.key.size = unit(.8,"lines"),
            legend.background = element_rect(colour = "black"))

    g <- cowplot::plot_grid(g1,g2,g3,g4, ncol = 2)

    if(doPlot){
      print(g)
    }

    if(returnPlot){
      outList$plot <- g
    } else {
      outList$plot <- NA
    }
  }

  outList$analysis = list(
      name = "Multifractal Detrended FLuctuation Analysis",
      logBaseFit = "log2",
      logBasePlot = "2")

  if(!silent){
    cat("\n~~~o~~o~~casnet~~o~~o~~~\n")
    cat(paste("\n",outList$analysis$name,"\n\n"))
    print(outList$MFSpectrum, digits = 3)
    cat("\n\n~~~o~~o~~casnet~~o~~o~~~\n")
  }

  return(invisible(outList[c(TRUE,TRUE,returnInfo,returnPlot,TRUE)]))

}


#' mono Hurst
#'
#' @inheritParams fd_dfa
#' @param TSm TS matrix with 2 columns `t` (1st) and `y` (second)
#'
#' @export
#' @keywords internal
#'
monoH <- function(TSm, scaleS, removeTrend = c("no","poly","adaptive","bridge")[2], polyOrder = 1, overlap = overlap, returnPLAW = FALSE, returnSegments = FALSE, removeRMSbelow = .Machine$double.eps){

  dfaRMS_scale <- vector("list",length(scaleS))

  if(max(scaleS)>NROW(TSm)/2){
    scaleS <- scaleS[scaleS<=NROW(TSm)/2]
  }

  F2 <- numeric(length(scaleS))

  for(ns in seq_along(scaleS)){
    dfaRMS <- plyr::ldply(ts_slice(y = TSm[,2], epochSz = scaleS[ns], overlap = overlap), function(sv){return(sqrt(mean(ts_detrend(sv,polyOrder=polyOrder)^2,na.rm = TRUE)))})
    dfaRMS$V1[dfaRMS$V1 <= removeRMSbelow^2] <- NA
    dfaRMS_scale[[ns]] <- dfaRMS$V1
    F2[ns] <- mean(dfaRMS_scale[[ns]]^2,na.rm = TRUE)^(1/2)
    if(is.infinite(log2(F2[ns]))){F2[ns] <- NA}
    rm(dfaRMS)
  }

  PLAW <- data.frame(size=scaleS%00%NA, bulk = F2%00%NA,
                     size.log2=log2(scaleS)%00%NA, bulk.log2 = log2(F2)%00%NA)
  H <- stats::lm(bulk.log2~size.log2, data = PLAW, na.action=stats::na.omit)$coefficients[2]
  if(any(returnPLAW,returnSegments)){
    return(list(H = H,
                PLAW = PLAW,
                segments = dfaRMS_scale)[c(TRUE,returnPLAW,returnSegments)])
  } else {
    return(H)
  }
}

#
#
# ts_sampEn <- function (y, edim = 2, r = 0.2 * sd(ts), tau = 1){
#
#   stopifnot(is.numeric(ts), is.numeric(edim))
#   if (tau > 1) {
#     s <- seq(1, length(ts), by = tau)
#     ts <- ts[s]
#   }
#   N <- length(ts)
#   correl <- numeric(2)
#   datamat <- zeros(edim + 1, N - edim)
#   for (i in 1:(edim + 1)) datamat[i, ] <- ts[i:(N - edim +
#                                                   i - 1)]
#   for (m in edim:(edim + 1)) {
#     count <- zeros(1, N - edim)
#     tempmat <- datamat[1:m, ]
#     for (i in 1:(N - m - 1)) {
#       X <- abs(tempmat[, (i + 1):(N - edim)] - repmat(tempmat[,
#                                                               i, drop = FALSE], 1, N - edim - i))
#       dst <- apply(X, 2, max)
#       d <- (dst < r)
#       count[i] <- sum(d)/(N - edim)
#     }
#     correl[m - edim + 1] <- sum(count)/(N - edim)
#   }
#   return(log(correl[1]/correl[2]))
# }


#' Sample Entropy
#'
#'
#' @inheritParams fd_dfa
#' @param m The size of the window in which tho evaluate whether a pattern repeats (default = `2`)
#' @param r A factor that will determine the threshold for similarity of values, calculated as r x D (default = `0.2`)
#' @param D Commonly the standard deviation of the time series, the similarity threshold will be calculated as r x D. Note that if the series is detrended and/or standardised and `D = NA` the standard deviation will be calculated after the transformations (default = `NA`)
#' @param relativeEntropy The relative entropy, SampEn / (-1 * log(1/length(y))) will be returned (default = `FALSE`)
#' @param transformBefore Detrend/standardise before coarse graining. If set to `FALSE`, each coarsegrained series will be detrended/standardised separately (default = `TRUE`)
#'
#'
#' @return The sample entropy (SampEn) of the time series y.
#'
#' @export
#'
#' @seealso info_MSE
#'
#' @family Information based complexity measures
#'
#' @examples
#'
#'
#' y <- rnorm(100)
#'
#' # Similarity threshold is r x D = 0.2 * sd(y)
#' inf_SampEn(y)
#'
#' # Similarity threshold is r = 0.2
#' inf_SampEn(y, D = 1)
#'
inf_SampEn <- function(y,
                      m = 2,
                      r = 0.2,
                      D = NA,
                      fs = NULL,
                      standardise = c("none","mean.sd","median.mad")[1],
                      transformBefore = TRUE,
                      removeTrend = c("no","poly","adaptive","bridge")[1],
                      polyOrder=1,
                      relativeEntropy = FALSE,
                      returnInfo = FALSE,
                      silent = FALSE){


  if(!stats::is.ts(y)){
    if(is.null(fs)){fs <- 1}
    y <- stats::ts(y, frequency = fs)
    if(!silent){cat("\n\nfd_SampEn:\tSample rate was set to 1.\n\n")}
  }

  N <- length(y)
  # Simple linear detrending.
  if(removeTrend%in%"poly"){
    y <- ts_detrend(y, polyOrder = polyOrder)
  }

  # Standardise by N
  if(any(standardise%in%c("mean.sd","median.mad"))){
    y <- ts_standardise(y, type = standardise,  adjustN = FALSE)
  }

  if(is.na(D)){
    D <- sd(y, na.rm = TRUE)
  }

  # (m+1)-element sequences
  N <- length(y)
  X <- matrix(0, nrow = N-m, ncol = m+1)
  for(i in 1:(m + 1)){
    X[,i] = y[((i-1)+1):(N-m+(i-1))]
  }

  # matching (m+1)-element sequences
  A <- sum(proxy::dist(X, method = "Chebyshev") < (r * D))

  # matching m-element sequences
  X <- X[, 1:m]
  B = sum(proxy::dist(X, method = "Chebyshev") < (r * D))

  # take log
  if(A==0 | B==0){
  e <- NA
  } else {
  e <- log(B / A)
  }

  if(!relativeEntropy){
    attr(e,"Relative Entropy") <- e/(-1 * log(1/N))
    return(e)
  } else {
    attr(e,"Sample Entropy") <- e
    return(e/(-1 * log(1/N)))
  }
}


#' Multi-Scale Entropy
#'
#' Calculate the Multi-Scale (Sample) Entropy of a time series.
#'
#' @inheritParams ts_coarsegrain
#' @inheritParams inf_SampEn
#' @inheritParams fd_dfa
#' @param transformBefore Detrend/standardise before coarse graining. If set to `FALSE`, each coarsegrained series will be detrended/standardised separately (default = `TRUE`)
#'
#' @note The transformation settings (detrending and/or normalisation), `transformBefore` will determine the value of the product r x D if `D = NA`.
#'
#' @return MSE
#'
#' @export
#'
#' @family Information based complexity measures
#'
#' @examples
#'
#' y <- rnorm(100)
#'
#' out <- inf_MSE(y, r=.5, scaleMin=1, scaleMax=10, returnInfo = TRUE, doPlot = FALSE)
#'
#'
inf_MSE <- function(y,
                    summaryFunction = c("mean","median","min","max")[1],
                    retainLength = FALSE,
                    m = 2,
                    r = 0.2,
                    D = NA,
                    fs = NULL,
                    transformBefore = TRUE,
                    removeTrend = c("no","poly","adaptive","bridge")[1],
                    polyOrder=1,
                    standardise = c("none","mean.sd","median.mad")[1],
                    adjustSumOrder = FALSE,
                    scaleMin = 1,
                    scaleMax = floor(NROW(y)/10),
                    scaleS = NA,
                    overlap = 0,
                    dataMin = 4,
                    relativeEntropy = FALSE,
                    doPlot = FALSE,
                    returnPlot = FALSE,
                    returnPLAW = FALSE,
                    returnInfo = FALSE,
                    silent = FALSE,
                    noTitle = FALSE,
                    tsName="y"
){

  if(transformBefore){

    # Simple linear detrending.
    if(removeTrend=="poly"){
      y <- ts_detrend(y, polyOrder = polyOrder)
    } else {
      if(removeTrend!="no"){
        message("Only polynomial detrending implemented, use removeTrend = 'poly' ")
      }
    }

    # Standardise by N
    if(any(standardise%in%c("mean.sd","median.mad"))){
      y <- ts_standardise(y, type = standardise,  adjustN = FALSE)
    }

    if(is.na(D)){
      D <- sd(y, na.rm = TRUE)
    }

    standardise <- "none"
    detrend <- FALSE

  }

  #N <- length(y)

  scaleFactor <- scaleMin:scaleMax
  out <- data.frame(scale = scaleFactor, SampEn = NA, SampEn_Relative = NA)

  for(g in seq_along(scaleFactor)){

    tmp <- ts_coarsegrain(y = y,
                          grain = scaleFactor[g],
                          summaryFunction = summaryFunction,
                          retainLength = retainLength)

    ent <- inf_SampEn(y = tmp,
                      m = m,
                      r = r,
                      D = D,
                      fs = fs,
                      standardise = "none",
                      removeTrend = "no",
                      polyOrder = polyOrder,
                      relativeEntropy = relativeEntropy,
                      returnInfo = returnInfo,
                      transformBefore = FALSE,
                      silent = TRUE)

    out$SampEn[g] <- ent

    if(!relativeEntropy){
      out$SampEn_Relative[g] <- attr(ent,"Relative Entropy")
    } else {
      out$SampEn_Relative[g] <- attr(ent,"Sample Entropy")
    }

    rm(tmp, ent)
  }


  fitRange <- 2:(length(scaleFactor)-2)

  lmfit1        <- stats::lm(out$SampEn ~ out$scale)
  lmfit2        <- stats::lm(out$SampEn[fitRange] ~ out$scale[fitRange])

  outList <- list(
    PLAW  =  cbind.data.frame(freq.norm = stats::frequency(y)/out$scale,
                              size = out$scale,
                              bulk = out$SampEn),
    fullRange = list(sap = stats::coef(lmfit1)[2],
                     H = NA, #1+stats::coef(lmfit1)[2] + Hadj,
                     FD = NA, #sa2fd_sda(stats::coef(lmfit1)[2]),
                     fitlm1 = lmfit1,
                     method = paste0("Full range (n = ",length(out$scale),")\nSlope = ",round(stats::coef(lmfit1)[2],2))), #," | FD = ",round(sa2fd_sda(stats::coef(lmfit1)[2]),2))),
    fitRange  = list(sap = stats::coef(lmfit2)[2],
                     H = NA, #1+stats::coef(lmfit2)[2] + Hadj,
                     FD = NA, #sa2fd_sda(stats::coef(lmfit2)[2]),
                     fitlm2 = lmfit2,
                     method = paste0("Fit range (n = ",length(out$scale[fitRange]),")\nSlope = ",round(stats::coef(lmfit2)[2],2))), #," | FD = ",round(sa2fd_sda(stats::coef(lmfit2)[2]),2))),
    info = out,
    plot = NA,
    analysis = list(
      name = "Multi Scale Entropy",
      logBaseFit = "no",
      logBasePlot = "no")
  )

  if(doPlot|returnPlot){

    xlabel <- "Scale Factor"
    ylabel <- "Sample Entropy"

    if(noTitle){
      title <- " "
    } else {
      title <- "regression (MSE)"
    }

    logSlopes <- data.frame(x = c(outList$PLAW$size[1:NROW(outList[[2]]$fitlm1$fitted.values)],
                                  outList$PLAW$size[1:NROW(outList[[3]]$fitlm2$fitted.values)]),
                            outList = c(outList$PLAW$bulk[1:NROW(outList[[2]]$fitlm1$fitted.values)],
                                        outList$PLAW$bulk[1:NROW(outList[[3]]$fitlm2$fitted.values)]),
                            Method = c(rep(outList[[2]]$method,NROW(outList[[2]]$fitlm1$fitted.values)),
                                       rep(outList[[3]]$method,NROW(outList[[3]]$fitlm2$fitted.values))))

    g <- ggplot2::ggplot(data.frame(outList$PLAW), ggplot2::aes_(x=~size,y=~bulk), na.rm=TRUE) +
      ggplot2::geom_point() +
      ggplot2::ggtitle(label = title, subtitle = "Multi Scale Entropy") +
      ggplot2::geom_smooth(data = logSlopes,  ggplot2::aes_(x=~x,y=~y, colour = ~Method, fill = ~Method), method="lm", alpha = .2) +
      ggplot2::scale_x_continuous(name = xlabel) + #, breaks = breaksX) +
      ggplot2::scale_y_continuous(name = ylabel)  #, breaks = breaksY) +
      ggplot2::scale_color_manual(values = c("red3","steelblue")) +
      ggplot2::scale_fill_manual(values = c("red3","steelblue")) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.minor =  ggplot2::element_blank(),
                     legend.text = ggplot2::element_text(margin = ggplot2::margin(t = 5,b = 5, unit = "pt"), vjust = .5),
                     plot.margin = ggplot2::margin(t = 5,b = 5, r = 5,l = 5, unit = "pt"))


    #g <- plotFD_loglog(fd.OUT = outList, title = title, subtitle = tsName, logBase = "no", ylabel = "Multi Scale Entropy", doPlot = doPlot)

    if(doPlot){
      print(g)
    }
    if(returnPlot){
      outList$plot <- g
    }
  }


  if(returnInfo){returnPLAW<-TRUE}

  if(!silent){
    cat("\n~~~o~~o~~casnet~~o~~o~~~\n")
    cat(paste("\n",outList$analysis$name,"\n\n",outList$fullRange$method,"\n\n",outList$fitRange$method))
    cat("\n\n~~~o~~o~~casnet~~o~~o~~~\n")
  }

  return(invisible(outList[c(returnPLAW,TRUE,TRUE,returnInfo,returnPlot, TRUE)]))

}




#
# set.seed(100)
# z <- dispersion(rnorm(1024))
# graphics::plot(log(z$scale),log(z$sd))
# #

# trace(ts_detrend,edit=T)
# seq(1,length(X),by=4096)
#
# z<-ts_slice(TSm,scale[1])
# z[[1]][,2]
#
# Hglobal <-
#
# segments <- plyr::laply(scale,function(s) floor(length(X)/s))
# IDv <- llply(segments,slice.index,)
# segv <- function(X,segments){
#
#
#
#   seq((((v-1)*scale[ns])+1),(v*scale[ns]),length=scale[ns])
#
# }
# for(ns in seq_along(scale)){
#   segv[ns] <- floor(length(X)/scale[ns])
#   for(v in 1:segv[ns]){
#     Index <- seq((((v-1)*scale[ns])+1),(v*scale[ns]),length=scale[ns])
#     Cslope<- polyfit(Index,X[Index],m)
#     fit   <- polyval(Cslope,Index)
#     RMS_scale[[ns]][v] <- sqrt(mean((X[Index]-fit)^2))
#     rm(Cslope, fit, Index)
#   }
#   for(nq in seq_along(qq)){
#     qRMS[[nq]][1:segv[ns]] <- RMS_scale[[ns]]^qq[nq]
#     Fq[[nq]][ns] <- mean(qRMS[[nq]][1:segv[ns]])^(1/qq[nq])
#   }
#   Fq[[which(qq==0)]][ns] <- exp(0.5*mean(log(RMS_scale[[ns]]^2)))
# }
#
# for(nq in seq_along(qq)){
#   Cslope <- polyfit(log2(scale),log2(Fq[[nq]]),1)
#   Hq[nq] <- Cslope[1]
#   qRegLine[[nq]] <- polyval(Cslope,log2(scale))
#   rm(Cslope)
# }
#
# tq <- (Hq*qq)-1
# hq <- diff(tq)/diff(qq)
# Dq <- (qq[1:(length(qq)-1)]*hq) - (tq[1:(length(qq)-1)])
#
# graphics::plot(hq,Dq,type="l")
#
#
# qq<-c(-10,-5,seq(-2,2,.1),5,10)
