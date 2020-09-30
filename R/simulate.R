

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


# Generate noise ----

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
#' @note This R code was adapted from a Matlab script called `powernoise.m` by Max Little. The script contained the following commented text:
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

  checkPkg("stabledist")

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

  checkPkg("stabledist")

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

  checkPkg("stabledist")

  if(!all(alpha%(]%c(0,2),beta%[]%c(-1,1),gamma%()%c(0,Inf),location%()%c(-Inf,Inf))){
    stop(paste("\nValid parameter values:\n alpha in (0,2]\n beta  in [-1,1]\n scale in (0,Inf)\n location in (-Inf,Inf)"))}
  if(alpha%)(%c(1,2)){warning("Not a Levy-Pareto distribution if alpha <= 1 or alpha >= 2")}

  df <- data.frame(matrix(rep(NA,N*ndims),ncol = ndims, dimnames = list(NULL,paste0("dim",1:ndims))))
  for(c in 1:NCOL(df)){df[,c] <- cumsum(stabledist::rstable(n=N, alpha = alpha, beta = beta, gamma = scale, delta = location, pm = 2))}
  return(df)
}


# Smallworld test ------

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

