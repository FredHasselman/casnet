
# Load required libraries --------------------------------

init <- function(){
  library(devtools)
  library(rio)
  library(plyr)
  library(tidyverse)
  library(nonlinearTseries)
  library(zoo)
  library(xts)
  library(Matrix)
  library(flexclust)
  #library(NetComp)
  library(igraph)
  library(lattice)
  library(latticeExtra)
  library(grid)
  library(gridExtra)
  source('~/Documents/GitHub/manylabRs/manylabRs/R/fRedsRutils.R')
}


#' Rose tinted infix
#'
#'
#' @param x If \code{x} is any of \code{Inf,-Inf,NA,NaN,NULL,length(x)==0}, it will return \code{y}; otherwise it returns \code{x}.
#' @param y The value to return in case of catastrophy \code{>00<}
#'
#' @export
#' @author Fred Hasselman
#' @description When your functions wear these rose tinted glasses, the world will appear to be a nicer, fluffier place.
#' @seealso purrr::%||%
#' @examples
#' Inf %00% NA
#'
#' numeric(0) %00% ''
#'
#' NA %00% 0
#'
#' NaN %00% NA
#'
#' NULL %00% NA
`%00%` <- function(x,y){
  l0<-isna<-isnan<-isinf<-isnll<-FALSE
  if(length(x)==0){
    l0=TRUE
  } else {
    if(all(is.na(x)))       isna =TRUE
    if(all(is.nan(x)))      isnan=TRUE
    if(all(is.infinite(x))) isinf=TRUE
    if(all(is.null(x)))     isnll=TRUE
  }
  if(any(l0,isna,isnan,isinf,isnll)){x<-y}
  return(x)
}


#' Create a timeseries profile
#'
#' @param y A 1D timeseries
#'
#' @return The profile
#' @export
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
  yc <-yc + y[1]
  return(yc)
}

#' Find Peaks or Wells
#'
#' @param x A time series
#' @param window Window in whcih to look for peaks or wells
#' @param includeWells Find wells?
#' @param minPeakDist Minimum distance between peaks or wells
#' @param minPeakHeight Minimum height / depth for a peak / well
#'
#' @return Index with peak or well coordinates
#' @export
#'
find_peaks <- function (x,
                        window        = 3,
                        includeWells  = FALSE,
                        minPeakDist   = round(window/2),
                        minPeakHeight = .2*diff(range(x, na.rm = TRUE))){

  fp <- function(x,window){
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pksl <- sapply(which(shape < 0), FUN = function(i){
      z <- i - window + 1
      z <- ifelse(z > 0, z, 1)
      w <- i + window + 1
      w <- ifelse(w < length(x), w, length(x))
      if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])){
        return(i + 1)
      } else {
        return(numeric(0))
      }
    })
    return(unlist(pksl))
  }

  pks <- fp(x,window)
  if(includeWells){
    wls <- fp(-1*x,window)
    pks <- sort(c(pks,wls))
  }

  distOK    <- diff(c(NA,pks))>=minPeakDist
  distOK[1] <- TRUE

  heightOK <- rep(FALSE,length(pks))
  for(p in seq_along(pks)){
    if(abs(x[pks[p]] - mean(x[c(max(1,(pks[p]-window)):(pks[p]-1),(pks[p]+1):min((pks[p]+window),length(x)))])) >= minPeakHeight){
      heightOK[p] <- TRUE
    }
  }
  return(pks[heightOK&distOK])
}


#' Repeat Copies of a Matrix
#'
#' @param X A matrix
#' @param m Multiply dim(X)[1] m times
#' @param n Multiply dim(X)[2] n times
#'
#' @return A repeated matrix
#' @export
#'
repmat <- function(X,m,n){
  # % REPMAT R equivalent of repmat (matlab)
  # % FORMAT
  # % DESC
  # % description not available.

  if (is.matrix(X))
  {
    mx = dim(X)[1]
    nx = dim(X)[2]
    out <- matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  }
  else if (is.vector(X)) {
    mx = 1
    nx = length(X)
    out <- matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  }
  else if (length(X) == 1)
  {
    out <- matrix(X,m, n)
  }
  return (out)
}



# spdiags <- function(B, d, m, n){
#   # % SPDIAGS description not available.
#   # % FORMAT
#   # % DESC
#   # % description not available
#
#
#   if (is.matrix(d))
#     d <- matrix(d,dim(d)[1]*dim(d)[2],1)
#   else
#     d <- arg2
#
#   p <- length(d)
#
#   A <- sparseMatrix(i = 1:arg3, j = 1:arg3, x = 0, dims= c(arg3,arg4))
#   cat("in spdiags.R ")
#   print(A)
#   cat(" arg3 ")
#   print(arg3)
#   cat(" arg4 ")
#   pring(arg4)
#
#   m <- dim(A)[1]
#   n <- dim(A)[2]
#
#   len<-matrix(0,p+1,1)
#   for (k in 1:p)
#     len[k+1] <- len[k]+length(max(1,1-d[k]):min(m,n-d[k]))
#
#   a <- matrix(0, len[p+1],3)
#   for (k in 1:p)
#   {
#     # % Append new d(k)-th diagonal to compact form
#     i <- t(max(1,1-d[k]):min(m,n-d[k]))
#     a[(len[k]+1):len[k+1],] <- c(i, i+d[k], B[(i+(m>=n)*d[k]),k])
#   }
#
#   res1 <- sparseMatrix(i = a[,1], j = a[,2], x = a[,3], dims = c(m,n))
#
#   return (res1)
# }



# Catch parameters
# params <- list(eDim  = eDim,
#                eLag  = eLag,
#                eRad  = eRad,
#                DLmin = DLmin,
#                VLmin = VLmin,
#                theiler = theiler,
#                JRP     = JRP,
#                distNorm       = c("EUCLIDEAN", "MAX", "MIN", "OP")[[1]],
#                returnMeasures = returnMeasures,
#                returnRPvector       = returnRPvector,
#                returnLineDist     = returnLineDist,
#                path_to_rp = path_to_rp,
#                saveOut = saveOut,
#                path_out = path_out,
#                file_ID  = file_ID)

# y2    = NULL
# eDim  = 1
# eLag  = 1
# eRad  = 1
# DLmin = 2
# VLmin = 2
# theiler = 0
# win     = min(length(y1),ifelse(is.null(y2),(length(y1)+1), length(y2)), na.rm = TRUE)
# step    = win
# JRP     = FALSE
# distNorm       = c("EUCLIDEAN", "MAX", "MIN", "OP")[[1]]
# returnMeasures = TRUE
# returnRPvector       = TRUE
# returnLineDist = FALSE
# path_to_rp = getOption("casnet.path_to_rp")
# saveOut    = FALSE
# path_out   = NULL
# file_ID    = NULL


# y1, y2, eDim, eLag, eRad, DLmin, VLmin, theiler, JRP, returnMeasures, returnRPvector, returnLineDist, path_to_rp, saveOut, path_out, file_ID = 0){
# -i <string>  	data filename (input)
# -j <string>  	filename of additional data for CRP/JRP (input)
# -J 	compute JRP instead of default CRP (only if 2nd data is given)
# -r <string>	filename recurrence plot (output)
# -o <string>	filename RQA measures (output)
# -p <string>	filename histogramme diagonal line lengths (output)
# -c	cummulative histogramme
# -f <string>	format for recurrence plot file (ASCII, TIF), default=ASCII
# -n <string>	distance norm (EUCLIDEAN, MAX, MIN, OP), default=EUCLIDEAN
# -m <number>	embedding dimension, default=1
# -t <number>	embedding delay, default=1
# -e <number>	threshold, default=1 (if negative, a distance plot will be created)
# -l <number>	minimal diagonal line, default=2
# -v <number>	minimal vertical line, default=2
# -w <number>	Theiler corrector, default=1
# -D <string>	delimiter in data file, default=TAB
# -d <string>	delimiter for RQA file, default=TAB
# -s	silent (no messages displayed)
# -V	print version information
# -h	print this help text



# (C)RQA ------------------------------
#' crqa_cl_main [internal]
#'
#' @param y1 y1
#' @param y2 y1
#' @param eDim y1
#' @param eLag y1
#' @param eRad y1
#' @param DLmin y1
#' @param VLmin y1
#' @param theiler y1
#' @param win y1
#' @param step y1
#' @param JRP y1
#' @param distNorm y1
#' @param returnMeasures y1
#' @param returnRPvector y1
#' @param returnLineDist y1
#' @param returnDistance y1
#' @param path_to_rp y1
#' @param saveOut y1
#' @param path_out y1
#' @param file_ID y1
#' @param silent y1
#' @param ...
#'
#' @keywords internal
#'
#' @export
crqa_cl_main <- function(y1,
                    y2    = NULL,
                    eDim  = 1,
                    eLag  = 1,
                    eRad  = 1,
                    DLmin = 2,
                    VLmin = 2,
                    theiler = 0,
                    win     = min(length(y1),ifelse(is.null(y2),(length(y1)+1), length(y2)), na.rm = TRUE),
                    step    = win,
                    JRP     = FALSE,
                    distNorm       = c("EUCLIDEAN", "MAX", "MIN", "OP")[[1]],
                    returnMeasures = TRUE,
                    returnRPvector       = TRUE,
                    returnLineDist = FALSE,
                    plot_recmat = c("noplot","recmat","distmat")[[1]],
                    path_to_rp = getOption("casnet.path_to_rp"),
                    saveOut    = FALSE,
                    path_out   = NULL,
                    file_ID    = NULL,
                    silent     = TRUE, ...){

  RQAmeasures <- list(
    RR      = 'Recurrence rate',
    DET     = 'Determinism',
    DET_RR  = 'Ratio DET/RR',
    LAM     = 'Laminarity',
    LAM_DET = 'Ratio LAM/DET',
    L_max   = 'maximal diagonal line length',
    L_mean  = 'mean diagonal line length',
    L_entr  = 'Entropy of diagonal line length distribution',
    DIV     = 'Divergence (1/L_max)',
    V_max   = 'maximal vertical line length',
    TT      = 'Trapping time',
    V_entr  = 'Entropy of vertical line length distribution',
    T1      = 'Recurrence times 1st type',
    T2      = 'Recurrence times 2nd type',
    W_max   = 'Max interval',
    W_mean  = 'Mean of interval',
    W_entr  = 'Entropy interval distribution',
    W_prob  = 'Prob.of intervals',
    F_min   = 'F min'
  )

  if(plot_recmat%in%"distmat"){-1 * eRad}

  if(is.null(file_ID)){
    fnames <- list.files(path=path_out)
    rqaID  <- grepl("(RQAplot|RQAhist|RQAmeasures)",fnames)
    if(any(rqaID)){
      file_ID <- max(as.numeric(gsub("(RQAplot|RQAhist|RQAmeasures|txt|[[:punct:]])+","",fnames[rqaID]),"[_]"), na.rm = TRUE)
    } else {
      file_ID <- 0
    }
  }

  if(any(is.infinite(file_ID),is.null(file_ID))){file_ID <-0}
  if(any(is.null(y2))){df <- df[,1]}


  tmpd  <- tempdir()
  tmpf1 <- tempfile(tmpdir = tmpd, fileext = ".dat")
  write.table(as.data.frame(y1), tmpf1, col.names = FALSE, row.names = FALSE, sep = "\t")

  # if(grepl("Error",try.CATCH(normalizePath(path_out,mustWork = TRUE))$value)){
  #
  # }

  fileSep <- ifelse(any(path_out%in%"/"),"/","\\")

  file_ID <- file_ID + 1
  plotOUT     <- normalizePath(file.path(path_out,paste0("RQAplot_",     file_ID, ".txt"),fsep = fileSep))
  measOUT     <- normalizePath(file.path(path_out,paste0("RQAmeasures_", file_ID, ".txt"),fsep = fileSep))
  histOUTdiag <- normalizePath(file.path(path_out,paste0("RQAhist_diag_",file_ID, ".txt"),fsep = fileSep))
  histOUThori <- normalizePath(file.path(path_out,paste0("RQAhist_hori_",file_ID, ".txt"),fsep = fileSep))


  if(any(is.null(y2))|any(is.na(y2%00%NA))){
    opts <- paste("-i", shQuote(tmpf1),
                  "-r", shQuote(plotOUT),
                  "-o", shQuote(measOUT),
                  "-p", shQuote(histOUTdiag),
                  "-q", shQuote(histOUThori),
                  "-m", eDim,
                  "-t", eLag,
                  "-e", eRad,
                  "-l", DLmin,
                  "-v", VLmin,
                  "-w", theiler,
                  "-n", shQuote(distNorm),
                  ifelse(silent,"-s",""))
  } else {
    tmpf2 <- tempfile(tmpdir = tmpd, fileext = ".dat")
    write.table(as.data.frame(y2), tmpf2, col.names = FALSE, row.names = FALSE, sep = "\t")
    opts <- paste("-i", shQuote(tmpf1),
                  "-j", shQuote(tmpf2),
                  "-r", shQuote(plotOUT),
                  "-o", shQuote(measOUT),
                  "-p", shQuote(histOUTdiag),
                  "-q", shQuote(histOUThori),
                  "-m", eDim,
                  "-t", eLag,
                  "-e", eRad,
                  "-l", DLmin,
                  "-v", VLmin,
                  "-w", theiler,
                  "-n", shQuote(distNorm),
                  ifelse(silent,"-s",""))
  }

  #closeAllConnections()

  # RCMD
  devtools::RCMD(shQuote(paste0(getOption("casnet.rp_prefix"),"rp")), options = opts, path = normalizePath(path.expand(path_to_rp), mustWork = FALSE), quiet = FALSE)

  measures     = try.CATCH(rio::import(normalizePath(gsub("[']+","",measOUT))))
  rpMAT        = try.CATCH(rio::import(normalizePath(gsub("[']+","",plotOUT))))
  disthistDiag = try.CATCH(rio::import(normalizePath(gsub("[']+","",histOUTdiag))))
  disthistHori = try.CATCH(rio::import(normalizePath(gsub("[']+","",histOUThori))))

  if(all(is.null(measures$warning),is.data.frame(measures$value))){
    measures <- measures$value
  } else {
    measures <- rbind.data.frame(rep(NA,length(RQAmeasures)))
  }
  colnames(measures) <-names(RQAmeasures)

  if(all(is.null(rpMAT$warning),is.data.frame(rpMAT$value))){
    rpMAT <- rpMAT$value
  } else {
    rpMAT <- data.frame(y1=NA,y2=NA,dist=NA)
  }
  colnames(rpMAT) <-c('y1','y2','dist')

  if(all(is.null(disthistDiag$warning),is.data.frame(grepl("Error",paste(disthistDiag$value))))){
    disthistDiag <- disthistDiag$value
  } else {
    disthistDiag <- data.frame(line.length=NA,freq=NA)
  }
  colnames(disthistDiag) <-c('diag.line.length','freq')

  if(all(is.null(disthistHori$warning),is.data.frame(grepl("Error",paste(disthistHori$value))))){
    disthistHori <- disthistHori$value
  } else {
    disthistHori <- data.frame(line.length=NA,freq=NA)
  }
  colnames(disthistHori) <-c('hori.line.length','freq')


  cat(paste0("[ID ",file_ID,"] Analysis completed... "))
  if(!saveOut){
    devtools::RCMD(getOption("casnet.sysdel"), options = paste("-f",measOUT),     quiet = TRUE)
    devtools::RCMD(getOption("casnet.sysdel"), options = paste("-f",plotOUT),     quiet = TRUE)
    devtools::RCMD(getOption("casnet.sysdel"), options = paste("-f",histOUTdiag), quiet = TRUE)
    devtools::RCMD(getOption("casnet.sysdel"), options = paste("-f",histOUThori), quiet = TRUE)
    cat("temporary files removed ...\n")
  } else {
    cat("files saved ...\n")
  }
  return(list(measures      = measures,
              rpMAT         = rpMAT,
              diag_disthist = disthistDiag,
              hori_disthist = disthistHori))
}


# devtools::RCMD(shQuote(paste0(getOption("casnet.rp_prefix"),"rp")), options = "-V", path = normalizePath(getOption("casnet.path_to_rp"), mustWork = FALSE), quiet = FALSE)

#' fast (C)RQA
#'
#' This function will run the \href{http://tocsy.pik-potsdam.de/commandline-rp.php}{commandline Recurrence Plots} executable provided by Norbert Marwan.
#'
#' @param y1 Time series 1
#' @param y2 Time series 2 for Cross Recurrence Analysis (default = \code{NULL})
#' @param eDim Embedding dimensions (default = \code{1})
#' @param eLag Embedding lag (default = \code{1})
#' @param eRad Radius on distance matrix (default = \code{1})
#' @param DLmin Minimum length of diagonal structure to be considered a line (default = \code{2})
#' @param VLmin Minimum length of vertical structure to be considered a line (default = \code{2})
#' @param theiler Theiler window (default = \code{0})
#' @param win Window to calculate the (C)RQA (default = minimum of length of \code{y1} or \code{y2})
#' @param step Stepsize for sliding windows (default = size of \code{win}, so no sliding window)
#' @param JRP Wether to calculate a Joint Recurrence Plot(default = \code{FALSE})
#' @param distNorm One of "EUCLIDEAN" (default), \code{"MAX", "MIN"}, or \code{"OP"} for an Order Pattern recurrence matrix
#' @param returnMeasures Return the (C)RQA measures? (default = \code{TRUE})
#' @param returnRPvector Return the recurrent points in a dataframe? (default = \code{TRUE})
#' @param returnLineDist Return the distribution of diagonal and horizontal line length distances (default = \code{FALSE})
#' @param plot_recmat Produce a plot of the recurrence matrix by calling \code{\link{recmat_plot}}, values can "recmat" (the thresholded recurrence matrix),"distmat" (the unthresholded recurrence matrix) or "noplot" (default = \code{"noplot"})
#' @param path_to_rp Path to the command line executable (default = path set during installation, use \code{getOption("casnet.path_to_rp")} to see)
#' @param saveOut Save the output to files? If \code{TRUE} and \code{pat_out = NA}, the current working directory will be used (default = \code{FALSE})
#' @param path_out Path to save output if \code{saveOut = TRUE} (default = \code{NULL})
#' @param file_ID A file ID which will be a prefix to to the filename if \code{saveOut = TRUE} (default = \code{NULL}, an integer will be added tot the file name to ensure unique files)
#' @param silent  Do not display any messages (default = \code{TRUE})
#' @param ... Additional parameters (currently not used)
#'
#'
#' @details The \code{rp} executable is installed when the function is called for the first time and is renamed to \code{rp}, from a platform specific filename located in the directory: \code{...\\casnet\\commandline_rp\\}.
#' The file is copied to the directory: \code{...\\casnet\\exec\\}
#' The latter location is stored as an option and can be read by calling \code{getOption("casnet.path_to_rp")}.
#'
#'
#' @section Troubleshooting:
#' Some notes on resolving errors with \code{rp}.
#'
#' \itemize{
#'  \item \emph{Copy failed} - Every time the function \code{crqa_cl()} is called it will check whether a log file \code{rp_instal_log.txt} is present in the \code{...\\casnet\\exec\\} directory. If you delete the log file, and call the function, another copy of the executable will be attempted.
#' \item \emph{Copy still fails and/or no permission to copy} - You can copy the approrpiate executable to any directory you have access to, be sure to rename it to \code{rp} (\code{rp.exe} on Windows OS). Then, either pass the path to \code{rp} as the argument \code{path_to_rp} in the \code{fast_crqa} function call, or, as a more permanent solution, set the \code{path_to_rp} option by calling \code{options(casnet.path_to_rp="YOUR_PATH")}. If you cannot acces the directory \code{...\\casnet\\commandline_rp\\}, download the appropriate executable from the \href{http://tocsy.pik-potsdam.de/commandline-rp.php}{commandline Recurrence Plots} page and copy to a directory you have access to. Then follow the instruction to set \code{path_to_rp}.
#' \item \emph{Error in execution of \code{rp}} - This can have a variety of causes, the \code{rp} executable is called using \code{\link[devtools]{RCMD}} and makes use of the \code{\link{normalizePath}} function with argument \code{mustWork = FALSE}. Problems caused by specific OS, machine, or, locale problems (e.g. the \code{winslash} can be reported as an \href{https://github.com/FredHasselman/casnet/issues}{issue on Github}). One execution error that occurs when the OS is not recognised properly can be resolved by chekcing \code{getOption("casnet.rp_prefix")}. On Windows OS this should return an empty character vector, on Linux or macOS it should return \code{"./"}. You can manually set the correct prefix by calling \code{options(casnet.rp_prefix="CORRECT OS PREFIX")} and fill in the prefix that is correct for your OS
#' }
#'
#'
#'
#' @return A list object containing 1-3 elements, depending on arguments requesting output.
#'
#' \enumerate{
#' \item \code{rqa_measures} - A list of the (C)RQA measures returned if \code{returnMeasures = TRUE}:
#' \itemize{
#'  \item RR = 'Recurrence rate'
#'  \item DET = 'Determinism'
#'  \item DET_RR = 'Ratio DET/RR'
#'  \item LAM = 'Laminarity'
#'  \item LAM_DET = 'Ratio LAM/DET'
#'  \item L_max = 'maximal diagonal line length'
#'  \item L_mean = 'mean diagonal line length'
#'  \item L_entr = 'Entropy of diagonal line length distribution'
#'  \item DIV =  'Divergence (1/L_max)'
#'  \item V_max = 'maximal vertical line length'
#'  \item TT = 'Trapping time'
#'  \item V_entr = 'Entropy of vertical line length distribution'
#'  \item T1 = 'Recurrence times 1st type'
#'  \item T2 = 'Recurrence times 2nd type'
#'  \item W_max = 'Max interval length'
#'  \item W_mean = 'Mean of interval lengths'
#'  \item W_entr = 'Entropy of interval length distribution'
#'  \item W_prob = 'Probability of interval'
#'  \item F_min = 'F min'
#' }
#' \item \code{rqa_rpvector} - The radius thresholded distance matrix (recurrence matrix), which can be visualised as a recurrence plot by calling \code{\link{recmat_plot}}. If a sliding window analysis is conducted this will be a list of matrices and could potentially grow too large to handle. It is recommended you save the output to disk by setting \code{saveOut = TRUE}.
#' \item \code{rqa_diagdist} - The distribution of diagonal line lengths
#' }
#' @note The platform specific \code{rp} command line executables were created by Norbert Marwan and obtained under a Creative Commons License from the website of the Potsdam Institute for Climate Impact Research at \url{http://tocsy.pik-potsdam.de/}.
#'
#' The full copyright statement on the website is as follows:
#'
#' © 2004-2017 SOME RIGHTS RESERVED
#'
#' University of Potsdam, Interdisciplinary Center for Dynamics of Complex Systems, Germany
#'
#' Potsdam Institute for Climate Impact Research, Transdisciplinary Concepts and Methods, Germany
#'
#' This work is licensed under a \href{https://creativecommons.org/licenses/by-nc-nd/2.0/de/}{Creative Commons Attribution-NonCommercial-NoDerivs 2.0 Germany License}.
#'
#' More information about recurrence analysis can be found on the \href{http://www.recurrence-plot.tk}{Recurrence Plot} website.
#'
#' @export
#'
crqa_cl <- function(y1,
                      y2    = NULL,
                      eDim  = 1,
                      eLag  = 1,
                      eRad  = 1,
                      DLmin = 2,
                      VLmin = 2,
                      theiler = 0,
                      win     = min(length(y1),ifelse(is.null(y2),(length(y1)+1), length(y2)), na.rm = TRUE),
                      step    = win,
                      JRP     = FALSE,
                      distNorm       = c("EUCLIDEAN", "MAX", "MIN", "OP")[[1]],
                      returnMeasures = TRUE,
                      returnRPvector       = TRUE,
                      returnLineDist = FALSE,
                      plot_recmat = c("noplot","recmat","distmat")[[1]],
                      path_to_rp = getOption("casnet.path_to_rp"),
                      saveOut    = FALSE,
                      path_out   = NULL,
                      file_ID    = NULL,
                      silent     = TRUE, ...){
require(plyr)
require(dplyr)

  if(!file.exists(normalizePath(paste0(getOption("casnet.path_to_rp"),"\rp_instal_log.txt"), mustWork = FALSE))){
    set_command_line_rp()
  }

  if(any(is.na(y1))){
    y1 <- y1[!is.na(y1)]
    warning("Removed NAs from timeseries y1 before (C)RQA")
  }

  if(any(is.na(y2%00%0))){
    y2 <- y2[!is.na(y2)]
    warning("Removed NAs from timeseries y2 before (C)RQA")
  }

  # Begin input checks

  if(is.null(path_out)){
    if(saveOut){
      path_out <- getwd()
    } else {
      path_out <- tempdir()
    }
  } else {
    path_out <- normalizePath(path_out, mustWork = TRUE)
  }

  if(!is.null(y2)){
    y2 <- as.zoo(y2)
    N1 <- NROW(y1)
    N2 <- NROW(y2)
    if(N1>N2){
      params$y2[N2:(N1+(N1-N2))] <- 0
    }
    if(N1<N2){
      y1[N1:(N2+(N2-N1))] <- 0
    }
    df    <- zoo::zoo(cbind(y1=y1,y2=y2))
  } else {
    df    <- zoo::zoo(cbind(y1=y1,y2=NA))
  }

  #if((win==min(length(y1),length(y2), na.rm = TRUE))&(step== 1)){
  if(is.null(y2)){
    cat("\nPerforming auto-RQA\n")
    if(theiler<1){theiler=1}
  } else {cat("\nPerforming cross-RQA\n")}
  #}
  if((win<min(length(y1),length(y2), na.rm = TRUE))){
    cat(paste("Calculating recurrence measures in a window of size",win,"taking steps of",step,"data points.\n\n"))
  }

  wlist <- zoo::rollapply(df,
                     width = win,
                     FUN = function(df){crqa_cl_main(y1             = df[,1],
                                                y2             = df[,2],
                                                eDim           = eDim,
                                                eLag           = eLag,
                                                eRad           = eRad,
                                                DLmin          = DLmin,
                                                VLmin          = VLmin,
                                                theiler        = theiler,
                                                JRP            = JRP,
                                                distNorm       = distNorm,
                                                returnMeasures = returnMeasures,
                                                returnRPvector = returnRPvector,
                                                returnLineDist = returnLineDist,
                                                plot_recmat    = plot_recmat,
                                                path_to_rp     = path_to_rp,
                                                saveOut        = saveOut,
                                                path_out       = path_out,
                                                file_ID        = file_ID,
                                                silent         = silent)},
                     by = step,
                     by.column = FALSE,
                     align = "left",
                     coredata = TRUE)


  wlist        <- unclass(wlist)
  rqa_measures <- ldply(wlist[,1]) %>% mutate(win = win, step = step, index = attr(wlist, "index"))
  rqa_rpvector <- ldply(wlist[,2]) %>% mutate(win = win, step = step, index = attr(wlist, "index"))
  rqa_diagdist <- ldply(wlist[,3]) %>% mutate(win = win, step = step, index = attr(wlist, "index"))
  rqa_horidist <- ldply(wlist[,4]) %>% mutate(win = win, step = step, index = attr(wlist, "index"))


  out <- list(rqa_measures = rqa_measures,
              rqa_rpvector = rqa_rpvector,
              rqa_diagdist = rqa_diagdist,
              rqa_horidist = rqa_horidist)
  if(saveOut){saveRDS(out,paste0(path_out,"CRQA_out",file_ID,".rds"))}
  return(out[c(returnMeasures,returnRPvector,returnLineDist)])
}

# emDims A range of embedding dimensions (default= 1 - maxDims)
# nSize The neighbourhood size used to estimate the number of nearest neighbours of a coordinate (default: 10\% of all points)

#' Find optimal (C)RQA parameters
#'
#' A wrapper for various algorithms used to find optimal embedding delay, number of embedding dimensions and radius.
#'
#' @param y Timeseries
#' @param maxDim Maximum number of embedding dimensions
#' @param maxLag Maximum embedding lag
#' @param emLag Optimal embedding lag (delay), e.g., provided by optimising algorithm.
#' @param emLags A range of embedding lags to consider (default)
#' @param ami.method Method
#' @param nnSizes  Neighbourhood size
#' @param nnThres  Threshold
#' @param diagPlot Plot the results
#'
#' @return A list object containing the optimal values and iteration history.
#'
#' @export
#'
crqa_parameters <- function(y,
                            maxDim   = 5,
                            maxLag   = round(length(y)/(maxDim+1)),
                            emLag    = tau(y,
                                           selection.methods = c("first.e.decay",
                                                                 "first.zero",
                                                                 "first.minimum"),
                                           maxLag =round(length(y)/(maxDim+1))
                            ),
                            emLags     = tau(y, maxLag = round(length(y)/(maxDim+1))),
                            ami.method = c("first.minimum"),
                            nnSizes  = c(.1,.3,.6,.9),
                            nnThres  = .1,
                            diagPlot = TRUE){

  if(!is.null(dim(y))){stop("y must be a 1D numeric vector!")}

  emDims  <-  1:maxDim

  y <- y[!is.na(y)]

  lagList <- list()
  cnt = 0
  for(N in seq_along(nnSizes)){
    for(L in seq_along(emLags$selection.method)){
      Nn.max <- NA
      for(D in seq_along(emDims)){
        cnt = cnt+1
        surrDims <- nonlinearTseries::buildTakens(as.numeric(y),
                                                  emDims[[D]],
                                                  emLags$opLag[L])
        rmat <- recmat(surrDims,surrDims)
        allN <- nonlinearTseries::findAllNeighbours(surrDims, radius = nnSizes[N]*max(rmat, na.rm = TRUE))
        Nn <- sum(laply(allN, length), na.rm = TRUE)
        if(D==1){Nn.max <- Nn}
        lagList[[cnt]] <- data.frame(Nsize        = nnSizes[N],
                                     emLag.method = emLags$selection.method[[L]],
                                     emLag = emLags$opLag[L],
                                     emDim = emDims[D],
                                     Nn    = Nn,
                                     Nn.max = Nn.max)
      }
    }
  }

  df        <- plyr::ldply(lagList)
  df$Nn.pct <- df$Nn/df$Nn.max

  opt <- plyr::ldply(unique(df$emLag), function(n){
    id <- which((df$Nn.pct<=nnThres)&(df$emLag==n)&(!(df$emLag.method%in%"maximum.lag")))
    if(length(id)>0){
      idmin <- id[df$emDim[id]==min(df$emDim[id], na.rm = TRUE)]
      if(length(idmin)>0){
        return(df[idmin,])
      }
    } else {
      return(df[which.min(df$Nn.pct[(df$emLag==n)&(!(df$emLag.method%in%"maximum.lag"))]),])
    }
  }
  )

  opt <- opt[opt$emDim==min(opt$emDim),][1,]

  opDim <- opt$emDim
  opLag <- opt$emLag

  #opDim <- min(df$emDim[df$Nn.pct<nnThres], na.rm = TRUE)
  # opLag <- tau(y,
  #              selection.methods = ami.method,
  #              maxLag =maxLag)$opLag[1]
  #opLag <- df$emLag[which.min(min(df$emDim[df$Nn.pct<nnThres], na.rm = TRUE))]
  #opRad = NULL

  df$emLag <- factor(df$emLag)
  df$Nsize <- factor(df$Nsize, labels = paste("nn size:",nnSizes))

  if(diagPlot){

    dfs <- data.frame(start= c(.5, hist(emDims)$mids),
                      stop = c(hist(emDims)$mids,max(emDims)+.5),
                      f=factor(seq_along(c(.5, hist(emDims)$mids))%%2))

    # use: alpha()
    #myPal <- RColorBrewer::brewer.pal(length(emLag),"Dark2")
    gNdims <- ggplot2::ggplot(df, aes(y = Nn.pct, x = emDim, colour = emLag)) +
      geom_rect(aes(xmin = start, xmax = stop, fill = f),
                ymin = -Inf, ymax = Inf, data = dfs, inherit.aes = FALSE) +
      scale_fill_manual(values = alpha(c("grey", "white"),.2), guide=FALSE) +
      geom_hline(yintercept = nnThres, linetype = 2, colour = "grey60") +
      geom_hline(yintercept = c(0,1),   colour = "grey60") +
      geom_hline(yintercept = 0.5, colour = "grey90") +
      geom_line(position  = position_dodge(.4)) +
      geom_point(position = position_dodge(.4)) +
      annotate("text",x=maxDim/3,y=nnThres, label="threshold", size = .8) +
      xlab("Embedding Dimension") + ylab("Nearest neigbours (% of max.)") +
      facet_wrap(~Nsize, ncol=2) +
      scale_x_continuous(breaks=emDims) +
      scale_y_continuous(breaks = c(nnThres,.5,1)) +
      scale_color_brewer("Lag",palette = "Set2") +
      theme_minimal() +
      theme(strip.background = element_rect(colour = "grey90", fill = "grey90"),
            strip.text.x = element_text(colour = "black", face = "bold"),
            panel.spacing = unit(1, "lines"),
            legend.background = element_rect(colour = "grey90",fill = "grey90"),
            legend.title = element_text(face = "bold"),
            legend.key = element_rect(colour = "grey90", fill = "grey90"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank()
      )

    tmi <-  nonlinearTseries::mutualInformation(y,
                                                lag.max = maxLag,
                                                do.plot = TRUE)

    dfMI <- data.frame(emDelay = tmi$time.lag[-1],
                       ami     = tmi$mutual.information[-1])

    gDelay <- ggplot2::ggplot(dfMI, aes(y = ami, x = emDelay)) +
      geom_line() +
      geom_vline(data = emLags, aes(colour=factor(selection.method),
                                    xintercept = opLag),
                 alpha = .3) +
      geom_point(data = emLags, aes(x = opLag,
                                    y=ami,
                                    colour=factor(selection.method)),
                 size=2) +
      xlab("Embedding Lag") +
      ylab("Average Mututal Information") +
      scale_color_brewer("Method",palette = "Set2") +
      theme_bw() +
      theme(legend.position = c(.95, .95),
            legend.justification = c("right", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6),
            legend.background = element_rect(colour = "grey90",fill = "grey90"),
            legend.title = element_text(face = "bold"),
            legend.key = element_rect(colour = "grey90", fill = "grey90"),
            #panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank()
      )

    gridExtra::grid.arrange(gDelay, gNdims, ncol=1, nrow=2)

  }

  return(list(optimLag  = opLag,
              optimDim  = opDim,
              #optimRad  = opRad,
              optimRow  = opt,
              optimData = df)
  )
}


#' Find a radius (fixed RR)
#'
#' @param RM Unthresholded Recurrence Matrix
#' @param startRadius Starting value for the radius (default = 10\% of the maxumum distance)
#' @param targetRR Target recurrence rate (default = \code{.05})
#' @param tol Tolerance for achieving \code{targetRR} (default = \code{0.1})
#' @param maxIter Maximum number of iterations to reach target (default = \code{100})
#' @param theiler Size of theiler window (default \code{0})
#' @param histIter Return iteration history? (defsault = \code{FALSE})
#'
#' @return A dataframe listing settings ussed to search for the radius, the radius found given the settings and the recurrence rate produced by the radius (either 1 row or the entire iteration history)
#' @export
#'
crqa_radius <- function(RM,
                        startRadius = .10*max(RM, na.rm = TRUE),
                        targetRR    = .05,
                        tol         = 0.1,
                        maxIter     = 100,
                        theiler     = 0,
                        histIter    = FALSE){

  if(!dplyr::between(tol,0,1)){stop("Argument tol must be dplyr::between 0 and 1.")}

  AUTO        <- ifelse(identical(as.vector(Matrix::tril(rmat,-1)),as.vector(Matrix::tril(t(rmat),-1))),TRUE,FALSE)
  recmat_size <- recmat_size(RM,AUTO,theiler)
  if(AUTO){
    cat(paste0("\nAuto-recurrence: Setting diagonal to 0 (matrix size: ",recmat_size,")\n"))
  } else {
    theiler = 0
  }

  RM <- bandReplace(RM,-theiler,theiler,0)

  tryRadius <- startRadius
  RR        <- 0
  iter      <- 0
  Converged <- FALSE
  seqIter <- 1:maxIter

  iterList <- data.frame(iter        = seqIter,
                         RR          = RR,
                         Radius      = tryRadius,
                         targetRR    = targetRR,
                         tolRRlo     = targetRR*(1-tol),
                         tolRRhi     = targetRR*(tol+1),
                         startRadius = startRadius,
                         recmat_size      = recmat_size,
                         AUTO        = AUTO,
                         Converged   = Converged, check.names = FALSE)

  exitIter <- FALSE

  while(!exitIter){

    iter <- iter+1

    RR = Matrix::nnzero(di2bi(RM,tryRadius),na.counted = FALSE)/recmat_size

    iterList[iter,] <-    cbind.data.frame(iter,
                                           RR,
                                           tryRadius,
                                           targetRR,
                                           targetRR*(1-tol),
                                           targetRR*(tol+1),
                                           startRadius,
                                           recmat_size,
                                           AUTO,
                                           Converged)

    if(any((dplyr::between(RR,(targetRR*(1-tol)),(targetRR*(tol+1)))),(iter>=maxIter))){
      exitIter <- TRUE
    }

    if(round(RR,digits = 2)>round(targetRR,digits = 2)){
      tryRadius <- tryRadius*(min(0.8,tol*2))
    } else {
      tryRadius <- tryRadius*(min(1.8,1+(tol*2)))
    }
  }

  iterList$Converged[iter] <- TRUE
  if(iter>=maxIter){
    warning("Max. iterations reached.")
    iterList$Converged[iter] <- FALSE
  }
  if(!dplyr::between(RR,(targetRR*(1-tol)),(targetRR*(tol+1)))){
    warning("Target RR not found, try increasing tolerance, or change starting value.")
    iterList$Converged[iter] <- FALSE
  }

  ifelse(histIter,id<-c(1:iter),id<-iter)
  return(iterList[id,])
}


#' Get (C)RQA measures from Recurrence Matrix
#'
#' Use `crqa_mat`
#'
#' @param RM A binary recurrence matrix
#' @param DLmin Minimal diagonal line length
#' @param VLmin Minimal vertical line length
#' @param HLmin Minimal horizontal line length
#' @param DLmax Maximal diagonal line length
#' @param VLmax Maximal vertical line length
#' @param HLmax Maximal horizontal line length
#' @param AUTO Is this an AUTO RQA?
#' @param chromatic Chromatic RQA?
#' @param matrices Return Matrices?
#' @param doHalf Analyse half of the matrix?
#' @param Nboot How many bootstraps?
#' @param CL Confidence Limit for bootstrap results
#'
#'
crqa_mat_measures <- function(RM,
                              DLmin = 2,
                              VLmin = 2,
                              HLmin = 2,
                              DLmax = length(diag(rmat))-1,
                              VLmax = length(diag(rmat))-1,
                              HLmax = length(diag(rmat))-1,
                              AUTO      = NULL,
                              chromatic = FALSE,
                              matrices  = FALSE,
                              doHalf    = FALSE,
                              Nboot     = NULL,
                              CL        = .95){

  # DLmin = 2
  # VLmin = 2
  # HLmin = 2
  # DLmax = length(diag(rmat))-1
  # VLmax = length(diag(rmat))-1
  # HLmax = length(diag(rmat))-1
  # chromatic = FALSE
  # matrices  = FALSE
  # CL        = .95

  require(dplyr)

  if(is.null(Nboot)){Nboot = 1}

  NRows <- NROW(RM)
  NCols <- NCOL(RM)
  mc.cores <- detectCores()
  if(Nboot<mc.cores) mc.cores <- Nboot

  tstart <- Sys.time()
  out    <- mclapply(1, function(i){
    crp_prep(matrix(RM[ceiling(NCols*NRows*runif(NCols*NRows))], ncol=NCols, nrow = NRows),
             radius= radius,
             DLmin = DLmin,
             VLmin = VLmin,
             HLmin = HLmin,
             DLmax = DLmax,
             VLmax = VLmax,
             HLmax = HLmax,
             AUTO  = AUTO,
             chromatic = chromatic,
             matrices  = matrices,
             doHalf    = doHalf)
  },
  mc.cores = mc.cores
  )
  dfori <- gather(as.data.frame(out), key = measure, value = value)
  tend  <- Sys.time()

  if(Nboot>1){

    cat(paste0("Bootstrapping Recurrence Matrix... ",Nboot," iterations...\n"))
    cat(paste0("Estimated duration: ", round((difftime(tend,tstart, unit = "mins")*Nboot)/max((round(mc.cores/2)-1),1), digits=1)," min.\n"))

    tstart <-Sys.time()
    bootOut <-  mclapply(1:Nboot, function(i){
      replicate <- as.data.frame(crp_prep(matrix(RM[ceiling(NCols*NRows*runif(NCols*NRows))],
                                                 ncol=NCols, nrow = NRows),
                                          radius= radius,
                                          DLmin = DLmin,
                                          VLmin = VLmin,
                                          HLmin = HLmin,
                                          DLmax = DLmax,
                                          VLmax = VLmax,
                                          HLmax = HLmax,
                                          AUTO  = AUTO,
                                          chromatic = chromatic,
                                          matrices  = matrices,
                                          doHalf    = doHalf))
      replicate$replicate = i
      return(replicate)
    },
    mc.cores = mc.cores
    )
    tend <- Sys.time()
    cat(paste0("Actual duration: ", round(difftime(tend,tstart, unit = "mins"), digits=1)," min.\n"))

    dfrepl <- ldply(bootOut)
    dfrepl <- gather(dfrepl, key = measure, value = value, -replicate)

    if(length(CL)==1){
      ci.lo <- (1-CL)/2
      ci.hi <- CL + ci.lo
    } else {
      ci.lo <- CL[1]
      ci.hi <- CL[2]
    }

    rqout <-  dfrepl %>% group_by(measure) %>%
      summarize(
        val      = NA,
        ci.lower = quantile(value, ci.lo, na.rm = TRUE),
        ci.upper = quantile(value, ci.hi, na.rm = TRUE),
        mean     = mean(value, na.rm = TRUE),
        sd       = sd(value, na.rm = TRUE),
        var      = var(value, na.rm = TRUE),
        N        = n(),
        se       = sd(value, na.rm = TRUE)/sqrt(n()),
        median   = mean(value, na.rm = TRUE),
        mad      = mad(value, na.rm = TRUE)
      )

    for(m in 1:nrow(rqout)){
      rqout$val[rqout$measure%in%dfori$measure[m]] <- dfori$value[m]
    }
  } else {
    rqout <- dfori
  }
  return(rqout)
}


# sizeIn <- dim(RM)
# mc.cores <- detectCores()
# if(Nboot<mc.cores) mc.cores <- Nboot
#
# rp.rg   <- function(data,mle){return(matrix(data[ceiling(mle$NCols*mle$NRows*runif(mle$NCols*mle$NRows))], ncol= mle$NCols, nrow = mle$NRows))}
#
# rp.mle <- list(NCols=sizeIn[1], NRows=sizeIn[2])
#
# rp.fun  <- function(data, ...){crp_prep(RP    = data,
#                                         radius= radius,
#                                         DLmin = DLmin,
#                                         VLmin = VLmin,
#                                         HLmin = HLmin,
#                                         DLmax = DLmax,
#                                         VLmax = VLmax,
#                                         HLmax = HLmax,
#                                         AUTO  = AUTO,
#                                         chromatic = chromatic,
#                                         matrices  = matrices,
#                                         doHalf    = doHalf)}
# rp.boot <- function(...){
#   boot(data      = RM,
#        statistic = rp.fun,
#        R         = Nboot,
#        ran.gen   = rp.rg,
#        mle       = rp.mle,
#        parallel = "multicore",
#        ncpus = detectCores()
#        )
#   #boot(RM, rp.fun, R = Nboot, sim = "parametric", ran.gen = rp.rg, mle = rp.mle)
# }

# replicates <-  boot(data      = RM,
#                     statistic = rp.fun,
#                     R         = Nboot,
#                     sim = "parametric",
#                     ran.gen   = rp.rg,
#                     mle       = rp.mle,
#                     parallel  = "multicore",
#                     ncpus     = detectCores()
#                     )
#
# cis <-  ldply(1:102, function(c){
#   ci <- boot.ci(replicates, type = c("norm"), conf = 0.95, t0 = c, t=1:9)
#   return(cbind.data.frame(measure = dimnames(replicates$t0)[[2]][c],
#                           value = replicates$t0[c],
#                           ci.lo = ci$normal[,2],
#                           ci.up = ci$normal[,3])
#          )})
# dfrepl <- llply(1:Nboot, function(i) rp.boot(i))

# bootout <- list() %>%
#   bootstrap(Nboot) %>%
#    do(as.data.frame( )
#       )


#' Get bootsrapped (C)RQA measures based on a recurrence matrix
#'
#' Measures based on horizontal line structures will be caluclated.
#'
#' @param rmat A distance matrix, or a matrix of zeroes and ones (you must set \code{radius = NULL})
#' @param radius Threshold for distance value that counts as a recurrence
#' @param DLmin Minimal diagonal line length (default = \code{2})
#' @param VLmin Minimal vertical line length (default = \code{2})
#' @param HLmin Minimal horizontal line length (default = \code{2})
#' @param DLmax Maximal diagonal line length (default = length of diagonal -1)
#' @param VLmax Maximal vertical line length (default = length of diagonal -1)
#' @param HLmax Maximal horizontal line length (default = length of diagonal -1)
#' @param AUTO Auto-recurrence? (default = \code{FALSE})
#' @param doHalf Analyse half of the matrix? (default = \code{FALSE})
#' @param chromatic Force chromatic RQA? (default = \code{FALSE})
#' @param matrices Return matrices? (default = \code{FALSE})
#' @param Nboot How many bootstrap replications? (default = \code{NULL})
#' @param CL Confidence limit for bootstrap results (default = \code{.95})

#' @return A list object containing (C)RQA measures (and matrices if requested)
#'
#' @export
#'
crqa_mat <- function(rmat,
                     radius= NULL,
                     DLmin = 2,
                     VLmin = 2,
                     HLmin = 2,
                     DLmax = length(diag(rmat))-1,
                     VLmax = length(diag(rmat))-1,
                     HLmax = length(diag(rmat))-1,
                     AUTO      = FALSE,
                     doHalf    = FALSE,
                     chromatic = FALSE,
                     matrices  = FALSE,
                     Nboot     = NULL,
                     CL        = .95){


  # Input should be a distance matrix, or a matrix of zeroes and ones with radius = NULL, output is a list
  # Fred Hasselman (me@fredhasselman.com) - August 2013

  if(is.null(AUTO)){
    AUTO <- ifelse(identical(as.vector(Matrix::tril(rmat,-1)),as.vector(Matrix::tril(t(rmat),-1))),TRUE,FALSE)
  }

  #uval <- unique(as.vector(rmat))
  if(all(as.vector(rmat)==0|as.vector(rmat)==1)){
    RM <- rmat
  } else {
    if(!is.null(radius)){
      RM <- di2bi(rmat,radius)
    } else{
      if(!chromatic){
        stop("Expecting a binary (0,1) matrix.\nUse 'crqa_radius()', or set 'chromatic = TRUE'")
      } else {
        stop("Chromatic RQA not implemented yet.")
      }
    }
  }
  rm(rmat)

  out <- crqa_mat_measures(RM,
                           DLmin = DLmin,
                           VLmin = VLmin,
                           HLmin = HLmin,
                           DLmax = DLmax,
                           VLmax = VLmax,
                           HLmax = HLmax,
                           AUTO  = AUTO,
                           chromatic = chromatic,
                           matrices  = matrices,
                           doHalf = doHalf,
                           Nboot  = Nboot,
                           CL     = CL)

  #
  #   if(is.null(Nboot)){Nboot = 1}
  #
  #   out <- crqa_mat_measures(RM,
  #                            radius= radius,
  #                            DLmin = DLmin,
  #                            VLmin = VLmin,
  #                            HLmin = HLmin,
  #                            DLmax = DLmax,
  #                            VLmax = VLmax,
  #                            HLmax = HLmax,
  #                            AUTO  = AUTO,
  #                            chromatic = chromatic,
  #                            matrices  = matrices,
  #                            doHalf = doHalf,
  #                            Nboot  = Nboot,
  #                            CL     = CL)
  #
  #   dfori <- gather(out, key = measure, value = value)
  #
  #   col.ind <- tbl_df(index(RM))
  #   row.ind <- tbl_df(sample(index(RM),size=nrow(RM)))
  #
  #   if(Nboot>1){cat(paste0("Bootstrapping Recurrence Matrix... ",Nboot," iterations.\n"))
  #     bootout <- col.ind  %>%
  #       bootstrap(Nboot) %>%
  #       do(crqa_mat_measures(RM[row.ind,unlist(.)],
  #                            radius= radius,
  #                            DLmin = DLmin,
  #                            VLmin = VLmin,
  #                            HLmin = HLmin,
  #                            DLmax = DLmax,
  #                            VLmax = VLmax,
  #                            HLmax = HLmax,
  #                            AUTO  = AUTO,
  #                            chromatic = chromatic,
  #                            matrices  = matrices,
  #                            doHalf = doHalf))
  #
  #     dfrepl <- gather(bootout, key = measure, value = value, -replicate)
  #
  #     if(length(CL)==1){
  #       ci.lo <- (1-CL)/2
  #       ci.hi <- CL + ci.lo
  #     } else {
  #       ci.lo <- CL[1]
  #       ci.hi <- CL[2]
  #     }
  #
  #     rqout <-  dfrepl %>% group_by(measure) %>%
  #       summarize(
  #         val     = NA,
  #         ci.low  = quantile(value, ci.lo, na.rm = TRUE),
  #         ci.high = quantile(value, ci.hi, na.rm = TRUE),
  #         mean    = mean(value, na.rm = TRUE),
  #         sd      = sd(value, na.rm = TRUE),
  #         var     = var(value, na.rm = TRUE),
  #         N       = n(),
  #         se      = sd(value, na.rm = TRUE)/sqrt(n()),
  #         median  = mean(value, na.rm = TRUE),
  #         mad     = mad(value, na.rm = TRUE)
  #       )
  #
  #     for(m in 1:nrow(rqout)){
  #       rqout$val[rqout$measure%in%dfori$measure[m]] <- dfori$value[m]
  #     }
  #   } else {
  #     rqout <- dfori
  #   }
  return(out)
}



#' Distance 2 binary
#'
#' Distance matrix to binary matrix based on threshold value
#'
#' @param distmat Distance matrix
#' @param radius The radius or threshold value
#' @param convMat Should the matrix be converted to an object of type Matrix
#'
#' @return A matrix with only 0s and 1s
#'
#' @export
di2bi <- function(distmat, radius, convMat = FALSE){

  if(grepl("Matrix",class(distmat))){
    distmat   <- as.matrix(distmat)
    convMat <- TRUE
  }

  # RP <- NetComp::matrix_threshold(distmat,threshold = radius, minval = 1, maxval = 0)

  RP <- matrix(0,dim(distmat)[1],dim(distmat)[2])
  RP[distmat <= radius] <- 1
  mostattributes(RP) <- attributes(distmat)

  if(!all(as.vector(RP)==0|as.vector(RP)==1)){warning("Matrix did not convert to a binary (0,1) matrix!!")}


  if(convMat){RP <- Matrix::Matrix(RP)}

  return(RP)
}

#' tau
#'
#' A wrapper for nonlinearTseries::timeLag
#'
#' @param y Time series
#' @param selection.methods Selecting an optimal embedding lag (default: Return "first.e.decay", "first.zero", "first.minimum", "first.value", where value is 1/exp(1))
#' @param value Threshold value for selection method first.value (default: 1/exp(1))
#' @param maxLag Maximal lag to consider (default: 1/5 of timeseries length)
#' @param nbins Number of bins for average mutual information function (default: 1/3 of number of unique values in timeseries)
#' @param ... Additional parameters
#'
#' @return The ami function with minima
#'
#' @export
#'
tau <- function(y,
                selection.methods = c("first.minimum"),
                maxLag = length(y)/4,
                ...){

  y   <- y[!is.na(y)]
  tmi <- nonlinearTseries::mutualInformation(y,
                                             lag.max = maxLag,
                                             do.plot = FALSE)

  lags <- numeric(length=length(selection.methods))
  cnt <- 0
  for(s.m in selection.methods){
    cnt <- cnt + 1
    lag <-try.CATCH(nonlinearTseries::timeLag(y,
                                              technique = "ami",
                                              selection.method = s.m,
                                              lag.max = maxLag,
                                              do.plot = FALSE))
    if(any(grepl("Error",lag$value))){lags[cnt] <- 1} else {lags[cnt] <- lag$value}
  }

  #id.peaks <- find_peaks(tmi$mutual.information, m = 3, wells = TRUE)
  id.min   <- tmi$time.lag[which.min(tmi$mutual.information)]
  #as.numeric(tmi$mutual.information[id.peaks[id.min]])

  out <-   cbind.data.frame(selection.method = c(selection.methods,
                                                 "global.minimum",
                                                 "maximum.lag"),
                            opLag = c(lags, id.min, maxLag)
  )

  tmi$mutual.information[tmi$time.lag%in%out$opLag]
  out$ami <- sapply(out$opLag, function(r) tmi$mutual.information[r])

  #tout <- summarise(group_by(out, opLag, ami), selection.method = paste0(unique(selection.method), collapse="|")

  return(out)
}


#' eDim
#'
#' @param y Time series
#' @param delay Embeding lag
#' @param maxDim Maximum number of embedding dimensions
#' @param ... Other arguments (not in use)
#'
#' @description A wrapper for nonlinearTseries::estimateEmbeddingDim
#'
#' @return Embedding dimensions
#' @export
#'
est_eDim <- function(y, delay = tau(y), maxDim = 20, threshold = .95, max.relative.change = .1, ...){
  cbind.data.frame(EmbeddingLag   = delay,
                   EmbeddingDim   = estimateEmbeddingDim(y,
                                                         time.lag  = delay,
                                                         threshold = threshold,
                                                         max.relative.change = max.relative.change,
                                                         max.embedding.dim = maxDim)
  )
}

#' nzdiags
#'
#' @description Get all nonzero diagonals of a binary matrix, or, diagonals specified as a vector by argument \code{d}. Output mimics Matlab's \code{[B,d] = spdiags(A)} function.
#'
#' @param A A binary (0,1) matrix.
#' @param d An optional vecotor of diagonals to extract.
#'
#' @author Fred Hasselman
#'
#' @return A list object with nonzero diagonals
#'
#' @export
#'
nzdiags <- function(A=NULL, d=NULL){
  # Loosely based on MATLAB function spdiags() by Rob Schreiber - Copyright 1984-2014 The MathWorks, Inc.
  #require(Matrix)

  if(grepl("matrix",class(A),ignore.case = TRUE)){

    if(all(A>0)){stop("All matrix elements are nonzero.")}

    # create an indicator for all diagonals in the matrix
    ind   <- col(A)-row(A)

    # Split the matrix!
    spd <- split(A, ind)

    if(is.null(d)){

      # Get diagonals which have nonzero elements
      keepID <- ldply(spd, function(di) any(di>0))
      nzdiag <- spd[keepID$V1]
      # Indices of nonzero diagonals
      dvec      <- as.numeric(keepID$.id[keepID$V1])

    } else {

      # Diagonals are specified
      d <- sort(as.vector(d))
      keepID <-  names(spd)%in%d
      nzdiag <- spd[keepID]
    }

    # Deal with uneven rows and cols
    m  <- NROW(A)
    n  <- NCOL(A)
    p  <- length(d)
    if(is.logical(A)){
      B <- matrix(FALSE,nrow = min(c(m,n), na.rm = TRUE), ncol = p)
    } else {
      B <- matrix(0, nrow = min(c(m,n), na.rm = TRUE), ncol = p)
    }

    for(i in seq_along(d)){
      B[zoo::index(nzdiag[[i]]),i] <- nzdiag[[i]]
    }
    colnames(B) <- paste(d)

    return(list(B = B,
                d = dvec))
  }
}


nzdiags.boot <- function(RP,d=NULL){
  # Loosely based on MATLAB function spdiags() by Rob Schreiber - Copyright 1984-2014 The MathWorks, Inc.
  #require(Matrix)

  if(grepl("matrix",class(RP),ignore.case = TRUE)){

    A <- RP

    if(all(A>0)){stop("All matrix elements are nonzero.")}

    # create an indicator for all diagonals in the matrix
    ind   <- col(A)-row(A)

    # Split the matrix!
    spd <- split(A, ind)

    if(is.null(d)){

      # Get diagonals which have nonzero elements
      keepID <- ldply(spd, function(di) any(di>0))
      nzdiag <- spd[keepID$V1]
      # Indices of nonzero diagonals
      d      <- as.numeric(keepID$.id[keepID$V1])

    } else {

      # Diagonals are specified
      d <- sort(as.vector(d))
      keepID <-  names(spd)%in%d
      nzdiag <- spd[keepID]
    }

    # Deal with uneven rows and cols
    m  <- NROW(A)
    n  <- NCOL(A)
    p  <- length(d)
    if(is.logical(A)){
      B <- matrix(FALSE,nrow = min(c(m,n), na.rm = TRUE), ncol = p)
    } else {
      B <- matrix(0, nrow = min(c(m,n), na.rm = TRUE), ncol = p)
    }

    for(i in seq_along(d)){
      B[zoo::index(nzdiag[[i]]),i] <- nzdiag[[i]]
    }
    colnames(B) <- paste(d)

    return(B)
  }
}


nzdiags.chroma <- function(RP, d=NULL){
  # Loosely based on MATLAB function spdiags() by Rob Schreiber - Copyright 1984-2014 The MathWorks, Inc.
  #require(Matrix)

  if(grepl("matrix",class(RP),ignore.case = TRUE)){

    A <- RP

    if(all(A>0)){stop("All matrix elements are nonzero.")}

    # create an indicator for all diagonals in the matrix
    ind   <- col(A)-row(A)

    # Split the matrix!
    spd <- split(A, ind)

    if(is.null(d)){

      # Get diagonals which have nonzero elements
      keepID <- ldply(spd, function(di) any(di>0))
      nzdiag <- spd[keepID$V1]
      # Indices of nonzero diagonals
      d      <- as.numeric(keepID$.id[keepID$V1])

    } else {

      # Diagonals are specified
      d <- sort(as.vector(d))
      keepID <-  names(spd)%in%d
      nzdiag <- spd[keepID]
    }

    # Deal with uneven rows and cols
    m  <- NROW(A)
    n  <- NCOL(A)
    p  <- length(d)
    if(is.logical(A)){
      B <- matrix(FALSE,nrow = min(c(m,n), na.rm = TRUE), ncol = p)
    } else {
      B <- matrix(0, nrow = min(c(m,n), na.rm = TRUE), ncol = p)
    }

    for(i in seq_along(d)){
      B[zoo::index(nzdiag[[i]]),i] <- nzdiag[[i]]
    }
    colnames(B) <- paste(d)

    return(B)
  }
}


#' lineDists
#'
#' @param RP A thresholded recurrence matrix (binary: 0 - 1)
#' @param d Vector of diagonals to be extracted from matrix \code{RP} before line length distributions are calculated. A one element vector will be interpreted as a windowsize, e.g., \code{d = 50} will extract the diagonal band \code{-50:50}. A two element vector will be interpreted as a band, e.g. \code{d = c(-50,100)} will extract diagonals \code{-50:100}. If \code{length(d) > 2}, the numbers will be interpreted to refer to individual diagonals, \code{d = c(-50,50,100)} will extract diagonals \code{-50,50,100}.
#' @param theiler Size of the theiler window, e.g. \code{theiler = 1} removes diagonal bands -1,0,1 from the matrix. If \code{length(d)} is \code{NULL}, 1 or 2, the theiler window is applied before diagonals are extracted. The theiler window is ignored if \code{length(d)>2}, or if it is larger than the matrix or band indicated by parameter \code{d}.
#' @param invert Relevant for Recurrence Time analysis: Return the distribution of 0 valued segments in nonzero diagonals/verticals/horizontals. This indicates the time between subsequent line structures.
#'
#' @description Extract lengths of diagonal and vertical line segments from a recurrence matrix.
#' @details Based on the Matlab function \code{lineDists} by Stefan Schinkel, Copyright (C) 2009 Stefan Schinkel, University of Potsdam, http://www.agnld.uni-potsdam.de
#'
#' References:
#' S. Schinkel, N. Marwan, O. Dimigen & J. Kurths (2009):
#' "Confidence Bounds of recurrence-based complexity measures
#' Physics Letters A,  373(26), pp. 2245-2250
#'
#' Copyright (C) 2009 Stefan Schinkel, University of Potsdam
#' \url{http://www.agnld.uni-potsdam.de}
#'
#' @author Fred Hasselman
#' @return A list object with distribution of line lengths.
#' @export
#'
lineDists <- function(RP,
                      d         = NULL,
                      theiler   = NULL,
                      Chromatic = FALSE,
                      invert    = FALSE,
                      AUTO      = NULL,
                      chromatic = FALSE,
                      matrices  = FALSE,
                      doHalf    = FALSE){

  # For boot()
  # RP <- RP[indices,]

  if(!all(as.vector(RP)==0|as.vector(RP)==1)){stop("Matrix should be a binary (0,1) matrix!!")}

  if(!is.null(d)){
    if(length(d)==1){d <- -d:d}
    if(length(d)==2){d <-  d:d}
  }
  if(!is.null(theiler)){
    if(length(d)<length(-theiler:theiler)){warning("Ignoring theiler window...")}
    RP <- double(band(x=RP,k1=-theiler,k2=theiler))
  }

  B <- nzdiags.boot(RP)

  # Get diagonal lines & pad with zeros
  diagonals   <- rbind.data.frame(rep(0,dim(B)[2]),
                                  B,
                                  rep(0,dim(B)[2])
  )
  #colnames(diagonals) <- paste(spdiagonals$d)

  # get vertical Lines & pad with zeros
  verticals <- rbind.data.frame(rep(0,dim(as.matrix(RP))[2]),
                                as.matrix(RP),
                                rep(0,dim(as.matrix(RP))[2])
  )
  colnames(verticals) <- paste(1:ncol(verticals))

  # get horizontal Lines & pad with zeros
  horizontals <- rbind.data.frame(rep(0,dim(t(as.matrix(RP)))[2]),
                                  t(as.matrix(RP)),
                                  rep(0,dim(t(as.matrix(RP)))[2])
  )
  colnames(horizontals) <- paste(1:ncol(horizontals))

  # Get indices of line lengths
  diagonals.ind   <- tidyr::gather(diagonals,   key = diagonal,  value = segment)
  verticals.ind   <- tidyr::gather(verticals,   key = vertical,  value = segment)
  horizontals.ind <- tidyr::gather(horizontals, key = horizontal,value = segment)

  # Get consecutive nonzero segments from indices, their difference is the segment lentgh
  diagonals.dist   <- which(diff(diagonals.ind$segment)==-1)-which(diff(diagonals.ind$segment)==1)
  verticals.dist   <- which(diff(verticals.ind$segment)==-1)-which(diff(verticals.ind$segment)==1)
  horizontals.dist <- which(diff(horizontals.ind$segment)==-1)-which(diff(horizontals.ind$segment)==1)

  return(list(diagonals.dist   = diagonals.dist,
              verticals.dist   = verticals.dist,
              horizontals.dist = horizontals.dist)
  )
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
dist.hamming <- function(X, Y=NULL, embedded=TRUE) {
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
  nrow(X) - H
}


#' bandReplace
#'
#' Sets a band of matrix diagonals to any given value
#'
#' @param mat A Matrix
#' @param lower Lower diagonal to be included in the band (should be $\\leq 0$)
#' @param upper Upper diagonal to be included in the band (should be $\\geq 0$)
#' @param value A single value to replace all values in the selected band
#'
#' @return A matrix in which the values in the selected diagonals have been replaced
#'
#' @export
#'
#' @examples
#' # Create a 10 by 10 matrix
#' library(Matrix)
#' m <- Matrix(rnorm(10),10,10)
#'
#' bandReplace(m,-1,1,0)   # Replace diagonal and adjacent bands with 0 (Theiler window of 1)
bandReplace <- function(mat, lower, upper, value = NA){
  if(lower>0){lower=-1*lower
  warning("lower > 0 ...\n using: -1*lower")
  }
  if(upper<0){upper=abs(upper)
  warning("upper > 0 ...\n using: abs(upper)")
  }
  if(all(lower==0,upper==0)){
    diag(mat) <- value
    message(paste0("lower and upper are both 0...\n using: diag(mat) <- ",value))
  }

  delta <- col(mat)-row(mat)
  mat[delta >= lower & delta <= upper] <- value

  return(mat)
}



#' Recurrence Matrix
#'
#' @param y1 y1
#' @param y2 y2
#' @param emDim emDim
#' @param emLag emLag
#' @param to.ts to.ts
#' @param to.sparse to.sparse
#' @param order.by order.by
#' @param method method
#' @param ... other
#'
#' @return RM
#' @export
#'
recmat <- function(y1, y2=NULL,
                   emDim = 1,
                   emLag = 1,
                   to.ts = NULL,
                   to.sparse = FALSE,
                   order.by = NULL,
                   method = "Euclidean", ...){

  if(is.null(y2)){y2 <- y1}

  if(!all(is.data.frame(y1),is.data.frame(y2))){
    y1 <- as.data.frame(y1)
    y2 <- as.data.frame(y2)
  }

  et1 <- lagEmbed(y1, emDim, emLag)
  et2 <- lagEmbed(y2, emDim, emLag)

  dmat <- proxy::dist(x = et1,
                      y = et2,
                      method = method,
                      diag = ifelse(identical(et1,et2),FALSE,TRUE))

  #dmat <- as.data.frame(dmat)
  dmat <- unclass(dmat)
  #rmat <- scales::rescale(rmat)
  if(all(!is.null(to.ts),!is.null(order.by))){

    dmat <-  switch(to.ts,
                    "xts" =  xts::xts(dmat, order.by = as_datetime(order.by)),
                    "zoo" =  zoo::zoo(dmat, order.by = as_datetime(order.by))
    )
  }
  if(!is.null(order.by)){
    colnames(dmat) <- paste(order.by)
    rownames(dmat) <- paste(order.by)
  }
  # if(is.null(to.ts)){
  #   dmat <- Matrix(dmat)
  # }

  if(to.sparse){
    dmat <- Matrix::Matrix(dmat)
  }

  attr(dmat,"eDims1") <- et1
  attr(dmat,"eDims2") <- et2
  attr(dmat,"AUTO")   <- ifelse(identical(et1,et2),TRUE,FALSE)

  # RM <- structure(.Data  = dmat,
  #                 emDims1.name = colnames(y1),
  #                 emDims2.name = colnames(y2),
  #                 eDims1 = et1,
  #                 eDims2 = et2,
  #                 AUTO   = ifelse(identical(et1,et2),TRUE,FALSE))
  # #class(RM) <- "recmat"

  return(dmat)
}


#' Plot (thresholded) distance matrix
#'
#' @param rmat A distance matrix or recurrence matrix
#' @param PhaseSpaceScale For display purposes: Rescale the data? (default = \code{TRUE})
#' @param doPlot Draw the plot? (default = \code{TRUE})
#' @param plotSurrogate Should a 2-panel comparison plot based on surrogate time series be added? If \code{rmat} has attributes \code{y1} and \code{y2} containing the time series data (i.e. it was created by a call to \code{\link{recmat}}), the following options are available: "RS" (random shuffle), "RP" (randomised phases), "AAFT" (amplitude adjusted fourier transform). If no timeseries data is included, the columns will be shuffled.  NOTE: This is not a surrogate test, just 1 surrogate is created from \code{y1}.
#' @param plotMeasures Print common (C)RQA measures in the plot
#'
#' @return A nice plot of the recurrence matrix.
#' @export
#'
recmat_plot <- function(rmat, PhaseSpaceScale= TRUE, title = "", doPlot = TRUE, plotSurrogate = NA, plotMeasures = FALSE){
  # require(grid)
  # require(ggplot2)
  # require(gtable)
  # require(cowplot)

  AUTO <-ifelse(identical(as.vector(Matrix::tril(rmat,-1)),as.vector(Matrix::tril(t(rmat),-1))),TRUE,FALSE)

  meltRP <- reshape2::melt(rmat)

  if(!all(as.vector(meltRP$value[!is.na(meltRP$value)])==0|as.vector(meltRP$value[!is.na(meltRP$value)])==1)){
    unthresholded = TRUE

    # For presentation purpose only
    if(PhaseSpaceScale){
      meltRP$value <- scales::rescale(meltRP$value)
    }
  } else {
    unthresholded = FALSE
    meltRP$value <- factor(meltRP$value)
  }


  gRP <-  ggplot2::ggplot(aes(x=Var1, y=Var2, fill = value), data= meltRP) +
    geom_raster()

  if(unthresholded){
    gRP <- gRP + scale_fill_gradient2(low      = "red3",
                                      high     = "steelblue",
                                      mid      = "white",
                                      na.value = scales::muted("slategray4"),
                                      midpoint = .5,
                                      limit    = c(min(meltRP$value),max(meltRP$value)),
                                      space    = "Lab",
                                      name     = "")
    rptheme <-     theme(panel.background = element_rect("grey99"),
                         panel.grid.major  = element_blank(),
                         panel.grid.minor  = element_blank(),
                         legend.background = element_blank(),
                         legend.key = element_blank(),
                         panel.border = element_blank(),
                         axis.ticks = element_blank(),
                         axis.text = element_blank(),
                         axis.title.x =element_blank(),
                         axis.title.y =element_blank(),
                         plot.margin = margin(0,0,0,0))

  } else {
    gRP <- gRP + scale_fill_manual(values = c("white","black"),
                                   na.value = scales::muted("slategray4"),
                                   name  = "", guide = "none")

    rptheme <-  theme(
      panel.background = element_rect("grey50"),
      panel.border = element_rect("grey50",fill=NA),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title.x =element_blank(),
      axis.title.y =element_blank(),
      plot.margin = margin(3,3,3,3))
  }

  gRP <- gRP +
    geom_abline(slope = 1,colour = "grey50", size = 1) +
    ggtitle(label=title, subtitle = ifelse(AUTO,"Auto-recurrence plot","Cross-recurrence plot")) +
    rptheme +
    coord_fixed(expand = FALSE)

  if(!is.null(attr(rmat,"eDims1"))){

    y1 <- data.frame(t1=attr(rmat,"eDims1"))
    y2 <- data.frame(t2=attr(rmat,"eDims2"))

    y1[,1] <- scales::rescale(y1[,1])
    gy1 <- ggplot2::ggplot(y1, aes(y=y1[,1], x= zoo::index(y1))) +
      geom_line() +  xlab(colnames(y1)) + ylab("") +
      geom_vline(xintercept = zoo::index(y1)[is.na(y1[,1])],
                 colour = scales::muted("slategray4"),alpha=.1, size=.5) +
      theme(panel.background = element_rect("grey99"),
            panel.grid.major  = element_blank(),
            panel.grid.minor  = element_blank(),
            legend.background = element_blank(),
            legend.key = element_blank(),
            panel.border = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x =element_text(colour = "black",angle = 0, vjust = +3),
            axis.title.y =element_blank(),
            plot.margin = margin(0,0,0,0)
      ) +
      coord_cartesian(expand = FALSE)  # +  coord_fixed(1/10)

    y2[,1] <- scales::rescale(y2[,1])
    gy2 <- ggplot2::ggplot(y2, aes(y=y2[,1], x= zoo::index(y2))) +
      geom_line() + xlab(colnames(y2)) + ylab("") +
      geom_vline(xintercept = zoo::index(y2)[is.na(y2[,1])],
                 colour = scales::muted("slategray4"),alpha=.1, size=.5) +
      theme(panel.background = element_rect("grey99"),
            panel.grid.major  = element_blank(),
            panel.grid.minor  = element_blank(),
            legend.background = element_blank(),
            legend.key = element_blank(),
            panel.border = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x =element_blank(),
            axis.title.y =element_text(colour = "black",angle = 90, vjust = -2),
            plot.margin = margin(0,0,0,0)) +
      coord_flip(expand = FALSE) +
      scale_y_reverse() #coord_fixed(1/2) +
  } else {
    gy1 <- gg.plotHolder()
    gy2 <- gg.plotHolder()
  }

  gr <- ggplot2::ggplotGrob(gRP)

  gindex <- subset(gr$layout, name == "panel")
  g <- gtable::gtable_add_cols(gr, grid::unit(2, "cm"),0)
  g <- gtable::gtable_add_rows(g, grid::unit(2,"cm"))
  g <- gtable::gtable_add_grob(g, ggplot2::ggplotGrob(gy2), t=gindex$t, l=1, b=gindex$b, r=gindex$l)
  g <- gtable::gtable_add_grob(g, ggplot2::ggplotGrob(gy1), t=9, l=2, b=11, r=6)

  if(doPlot){
  grid::grid.newpage()
  grid::grid.draw(g)
  }

  return(g)
}


#' recmat_size
#'
#' @param mat A Matrix object
#' @param AUTO Is the Matrix an Auto Recurrence Matrix? If so, the length of the diagonal will be subtracted from the matrix size, pass \code{FALSE} to prevent this behaviour. If \code{NULL} (default) \code{AUTO} will take on the value of \code{isSymmetric(mat)}.
#' @param theiler Should a Theiler window be applied?
#'
#' @return Matrix size for computation of recurrence measures.
#' @export
#'
#' @examples
#' # Create a 10 by 10 matrix
#' library(Matrix)
#' m <- Matrix(rnorm(10),10,10)
#'
#' recmat_size(m,TRUE,0)   # Subtract diagonal
#' recmat_size(m,FALSE,0)  # Do not subtract diagonal
#' recmat_size(m,NULL,0)   # Matrix is symmetrical, AUTO is set to TRUE
#' recmat_size(m,NULL,1)   # Subtract a Theiler window of 1 around and including the diagonal
recmat_size <- function(mat, AUTO=NULL, theiler = 0){
  if(is.null(AUTO)){
    AUTO <- Matrix::isSymmetric(unname(mat))
  }
  return(cumprod(dim(mat))[2] - ifelse((AUTO&theiler==0),length(Matrix::diag(mat)),
                                       ifelse(theiler>0,Matrix::nnzero(Matrix::band(mat,-theiler,theiler)),0)))
}

crp_empty <- function(){
  data.frame(
    Radius   = NA,
    RT       = NA,
    RR       = NA,
    DET      = NA,
    MEAN_dl  = NA,
    MAX_dl   = NA,
    ENT_dl   = NA,
    ENTrel_dl= NA,
    REP_tot  = NA,
    CoV_dl   = NA,
    DIV_dl   = NA,
    SING_dl  = NA,
    N_dl     = NA,
    ANI      = NA,
    LAM_vl   = NA,
    TT_vl    = NA,
    MAX_vl   = NA,
    ENT_vl   = NA,
    ENTrel_vl= NA,
    CoV_vl   = NA,
    REP_vl   = NA,
    DIV_vl   = NA,
    SING_vl  = NA,
    N_vl     = NA,
    LAM_hl   = NA,
    TT_hl    = NA,
    MAX_hl   = NA,
    ENT_hl   = NA,
    ENTrel_hl= NA,
    CoV_hl   = NA,
    REP_hl   = NA,
    DIV_hl   = NA,
    SING_hl  = NA,
    N_hl     = NA)
}

crp_calc <- function(RM,
                     radius= NULL,
                     DLmin = 2,
                     VLmin = 2,
                     HLmin = 2,
                     DLmax = length(diag(rmat))-1,
                     VLmax = length(diag(rmat))-1,
                     HLmax = length(diag(rmat))-1,
                     AUTO      = FALSE,
                     chromatic = FALSE,
                     matrices  = FALSE){

  recmat_size <- Msize(RM, auto=AUTO)

  #Total nr. recurrent points
  RT <- Matrix::nnzero(RM, na.counted = FALSE)
  #Proportion recurrence / Recurrence Rate
  RR <- RT/recmat_size
  #Get line segments
  lineSegments <- lineDists(RM)

  dlines <- lineSegments$diagonals.dist
  vlines <- lineSegments$verticals.dist
  hlines <- lineSegments$horizontals.dist

  #Frequency tables of line lengths
  freq_dl <- table(dlines)
  freq_vl <- table(vlines)
  freq_hl <- table(hlines)

  freqvec_dl <- as.numeric(names(freq_dl))
  freqvec_vl <- as.numeric(names(freq_vl))
  freqvec_hl <- as.numeric(names(freq_hl))

  #Number of recurrent points on diagonal, vertical and horizontal lines
  N_dl <- sum(freq_dl[dplyr::between(freqvec_dl,DLmin, DLmax)], na.rm = TRUE)
  N_vl <- sum(freq_vl[dplyr::between(freqvec_vl,VLmin, VLmax)], na.rm = TRUE)
  N_hl <- sum(freq_hl[dplyr::between(freqvec_hl,HLmin, HLmax)], na.rm = TRUE)

  #Determinism / Horizontal and Vertical Laminarity
  DET    <- N_dl/RT
  LAM_vl <- N_vl/RT
  LAM_hl <- N_hl/RT

  #anisotropy ratio
  ANI    <- (N_vl-N_hl)/N_dl

  # Singularities
  SING_dl <- sum(freq_dl[dplyr::between(as.numeric(names(freq_dl)),1,(DLmin-1))],na.rm = TRUE)/RT
  SING_vl <- sum(freq_vl[dplyr::between(as.numeric(names(freq_vl)),1,(VLmin-1))],na.rm = TRUE)/RT
  SING_hl <- sum(freq_hl[dplyr::between(as.numeric(names(freq_hl)),1,(HLmin-1))],na.rm = TRUE)/RT

  #Array of probabilities that a certain line length will occur (all >1)
  P_dl <- freq_dl[dplyr::between(as.numeric(names(freq_dl)), DLmin, DLmax)]/N_dl
  P_vl <- freq_vl[dplyr::between(as.numeric(names(freq_vl)), VLmin, VLmax)]/N_vl
  P_hl <- freq_hl[dplyr::between(as.numeric(names(freq_hl)), HLmin, HLmax)]/N_hl

  #Entropy of line length distributions
  ENT_dl <- -1 * sum(P_dl * log(P_dl))
  ENT_vl <- -1 * sum(P_vl * log(P_vl))
  ENT_hl <- -1 * sum(P_hl * log(P_hl))

  #Relative Entropy (Entropy / Max entropy)
  ENTrel_dl = ENT_dl/(-1 * log(1/DLmax))
  ENTrel_vl = ENT_vl/(-1 * log(1/VLmax))
  ENTrel_hl = ENT_hl/(-1 * log(1/HLmax))

  #Meanline
  MEAN_dl = mean(freq_dl)
  MEAN_vl = mean(freq_vl)
  MEAN_hl = mean(freq_hl)

  #Maxline
  MAX_dl = max(freq_dl)
  MAX_vl = max(freq_vl)
  MAX_hl = max(freq_hl)

  # REPetetiveness
  REP_tot <-(N_hl+N_vl-N_dl)/N_dl
  REP_hl  <- N_hl/N_dl
  REP_vl  <- N_vl/N_dl

  #Coefficient of determination
  CoV_dl = sd(freq_dl)/mean(freq_dl)
  CoV_vl = sd(freq_vl)/mean(freq_vl)
  CoV_hl = sd(freq_hl)/mean(freq_hl)

  #Divergence
  DIV_dl = 1/MAX_dl
  DIV_vl = 1/MAX_vl
  DIV_hl = 1/MAX_hl

  crqaMatrices  <- list()
  crqaMatrices <- list()

  #Output
  out <- data.frame(
    Radius   = radius,
    RT       = RT,
    RR       = RR,
    DET      = DET,
    MEAN_dl  = MEAN_dl,
    MAX_dl   = MAX_dl,
    ENT_dl   = ENT_dl,
    ENTrel_dl= ENTrel_dl,
    REP_tot  = REP_tot,
    CoV_dl   = CoV_dl,
    DIV_dl   = DIV_dl,
    SING_dl  = SING_dl,
    N_dl     = N_dl,
    ANI      = ANI,
    LAM_vl   = LAM_vl,
    TT_vl    = MEAN_vl,
    MAX_vl   = MAX_vl,
    ENT_vl   = ENT_vl,
    ENTrel_vl= ENTrel_vl,
    CoV_vl   = CoV_vl,
    REP_vl   = REP_vl,
    DIV_vl   = DIV_vl,
    SING_vl  = SING_vl,
    N_vl     = N_vl,
    LAM_hl   = LAM_hl,
    TT_hl    = MEAN_hl,
    MAX_hl   = MAX_hl,
    ENT_hl   = ENT_hl,
    ENTrel_hl= ENTrel_hl,
    CoV_hl   = CoV_hl,
    REP_hl   = REP_hl,
    DIV_hl   = DIV_hl,
    SING_hl  = SING_hl,
    N_hl     = N_hl)

  if(matrices){
    return(list(
      crqaMeasures = out,
      crqaMatrices = list(RM     = RM,
                          dlines = dlines,
                          vlines = vlines,
                          hlines = hlines,
                          freq_dl = freq_dl,
                          freq_vl = freq_vl,
                          freq_hl = freq_hl)
    )
    )
  } else {
    return(out)
  }
}


#' Prepare matrix
#'
#' @param RP Recurrence plot
#' @param radius Radiuc
#' @param DLmin Minimal diagonal line length
#' @param VLmin Minimal vertical line length
#' @param HLmin Minimal horizontal line length
#' @param DLmax Maximal diagonal line length
#' @param VLmax Maximal vertical line length
#' @param HLmax Maximal horizontal line length
#' @param AUTO Is this an AUTO RQA?
#' @param chromatic Chromatic RQA?
#' @param matrices Return Matrices?
#' @param doHalf Analyse half of the matrix?
#'
#' @return A prepped matrix
#' @export
#'
crp_prep <- function(RP,
                     radius= NULL,
                     DLmin = 2,
                     VLmin = 2,
                     HLmin = 2,
                     DLmax = length(diag(rmat))-1,
                     VLmax = length(diag(rmat))-1,
                     HLmax = length(diag(rmat))-1,
                     AUTO      = FALSE,
                     chromatic = FALSE,
                     matrices  = FALSE,
                     doHalf    = FALSE){

  out<-crp_calc(RP,
                radius= radius,
                DLmin = DLmin,
                VLmin = VLmin,
                HLmin = HLmin,
                DLmax = DLmax,
                VLmax = VLmax,
                HLmax = HLmax,
                AUTO  = AUTO,
                chromatic = chromatic,
                matrices  = matrices)

  if(doHalf){
    if(!AUTO){
      outLo <- crp_calc(Matrix::tril(RP,-1),
                        radius= radius,
                        DLmin = DLmin,
                        VLmin = VLmin,
                        HLmin = HLmin,
                        DLmax = DLmax,
                        VLmax = VLmax,
                        HLmax = HLmax,
                        AUTO  = AUTO,
                        chromatic = chromatic,
                        matrices  = matrices)

      outUp <- crp_calc(Matrix::triu(RP,-1),
                        radius= radius,
                        DLmin = DLmin,
                        VLmin = VLmin,
                        HLmin = HLmin,
                        DLmax = DLmax,
                        VLmax = VLmax,
                        HLmax = HLmax,
                        AUTO  = AUTO,
                        chromatic = chromatic,
                        matrices  = matrices)
      out <- cbind.data.frame(full  = out,
                              lower = outLo,
                              upper = outUp)
    } else {
      out<- cbind.data.frame(full  = out,
                             lower = crp_empty(),
                             upper = crp_empty())
    }
  }
  return(out)
}




# plotRM.crqa <- function(RM,
#                         useLattice = FALSE,
#                         plotVars   = NULL,
#                         plotTime   = NULL,
#                         plotMatrix = TRUE){
#   require(scales)
#
#   # RP dimensions
#   nr <- dim(RM)[1]
#   nc <- dim(RM)[2]
#
#   # Do some checks
#   plotVarsOK <- FALSE
#   plotTimeOK <- FALSE
#
#   if(any(class(RM)=="xts"|class(RM)=="zoo")){
#     if(!is.null(attr(RM,"eDims1"))){
#       plotVarsOK <- TRUE
#       plotTimeOK <- TRUE
#       plotVars   <- list(attributes(RM)$eDims1,attributes(RM)$eDims2)
#       matNAmes   <- c(attr(RM,"emDims1.name"),attr(RM,"emDims1.name"))
#       AUTO       <- attr(RM,"AUTO")
#       RM         <- unclass(RM)
#     }
#   } else {
#
#    AUTO <-  ifelse(isSymmetric(RM),TRUE,FALSE)
#
#     if(is.list(plotVars)){
#       if(all(lengths(plotVars)==c(nr,nc))){
#         yr <- plotVars[[1]]
#         yc <- plotVars[[2]]
#         plotVarsOK <- TRUE
#         if(any(ncol(yr)>10,ncol(yc)>10)){
#           yr <- yr[,1:10]
#           yc <- yc[,1:10]
#           warning("Detected more than 10 embedding dims, will plot 1-10.")}
#       } else {warning("Expecting: length(plotVars[[1]]) == dims(RM)[1] & length(plotVars[[2]]) == dims(RM)[2]\nFound: nrows = ",paste(lengths(plotVars),collapse ="; cols = "))}
#     } else {warning("plotVars not a list.\n Not plotting dimensions.")}
#
#     # Plot Time?
#     if(!is.null(plotTime)){
#       if(is.POSIXct(as_datetime(dimnames(RM)[[2]]))){
#         plotTimeOK <- TRUE
#       }
#     }
#   }
#
#   #   if(!is.list(plotTime)){
#   #     if(nc!=nr){
#   #       warning("Unequal rows and columns, expecting a list with 2 time variables.")
#   #     } else {
#   #       if(NROW(plotTime==nr)){
#   #         plotTimeRow <- plotTimeCol <- plotTime
#   #         plotTimeOK <- TRUE
#   #       } else {
#   #         warning("Time variable length should equal nrows(RM) & ncols(RM).")
#   #       } # if rows are not equal to length
#   #     } #if row and columns equal
#   #   } else {
#   #     if(all(lengths(plotTime)==c(nr,nc))){
#   #       plotTimeOK <- TRUE
#   #       plotTimeRow <- plotTime[[1]]
#   #       plotTimeCol <- plotTime[[2]]
#   #     } else {
#   #       warning("Expecting: length(plotTime[[1]]) == dims(RM)[1] & length(plotTime[[2]]) == dims(RM)[2]\nFound: nrows = ",paste(lengths(plotTime),collapse ="; cols = "))
#   #     }
#   #   }
#   #   if(plotTimeOK){
#   #     if(is.POSIXt(plotTimeRow)&is.POSIXt(plotTimeCol)){
#   #       plyr::amv_dimnames(RM) <-list(plotTimeRow, plotTimeCol)
#   #     } else {
#   #       warning("plotTime is not a POSIXct or POSIXlt object.")
#   #     }
#   #   }
#   # }
#
#   # AUTO or CROSS?
#   if(AUTO){
#     RM   <- Matrix::t(RM)
#   }
#
#   if(nc*nr>10^6){
#     message(paste("\nLarge RM with",nc*nr,"elements... \n >> This could take some time to plot!\n"))
#   }
#
#   if(useLattice){
#     require(lattice)
#     require(latticeExtra)
#
#     ifelse(AUTO,
#            Titles <- list(main="Auto Recurrence Plot", X = expression(Y[t]), Y = expression(Y[t])),
#            Titles <- list(main="Cross Recurrence Plot", X = expression(X[t]), Y = expression(Y[t]))
#     )
#
#     #binned <- ftable(as.vector(RM))
#     dimnames(RM) <- list(NULL,NULL)
#
#     #Thresholded or Distance matrix?
#     ifelse((all(as.vector(RM)==0|as.vector(RM)==1)),
#            distPallette <- colorRampPalette(c("white", "black"))(2),
#            distPallette <- colorRampPalette(colors = c("red3", "snow", "steelblue"),
#                                             space = "Lab",
#                                             interpolate = "spline")
#     )
#
#     lp <- levelplot(as.matrix(RM),
#                     colorkey = !AUTO,
#                     region = TRUE,
#                     col.regions = distPallette,
#                     useRaster = TRUE,
#                     aspect = nc/nr,
#                     main = Titles$main,
#                     xlab = Titles$X,
#                     ylab = Titles$Y)
#
#
#     if(plotVarsOK){
#       lp
#     }
#
#     if(plotMatrix){
#       lp
#       #grid.arrange(lp, txt, ncol=2, nrow=1, widths=c(4, 1), heights=c(4))
#     }
#
#     return(lp)
#
#   } else {
#
#     require(gridExtra)
#     require(ggplot2)
#     require(reshape2)
#
#     meltRP  <- reshape2::melt((RM))
#
#     if(!all(dplyr::between(meltRP$value[!is.na(meltRP$value)],0,1))){
#       meltRP$value <- scales::rescale(meltRP$value)
#     }
#
#   # #  if(plotTimeOK){
#   #     #if( (year(last(paste0(meltRP$Var1)))- year(meltRP$Var1[1]))>=1){
#   #       meltRP$Var1 <- parse_date_time(as.character(meltRP$Var1), c("ymd", "ymd HM", "ymd HMS"))
#   #       meltRP$Var2 <- parse_date_time(as.character(meltRP$Var2), c("ymd", "ymd HM", "ymd HMS"))
#   #     #}
#   #  # } else {
#   #     meltRP$Var1 <- index(meltRP$value)
#   #     meltRP$Var2 <- index(meltRP$value)
#   #   #}
#
#   ##  colnames(meltRP)[1:2] <- matNAmes
#
#     grp <- ggplot(data = meltRP, aes(x = meltRP[,2],
#                                      y = meltRP[,1],
#                                      fill = value)) + geom_raster()
#
#     if(all(as.vector(RM)==0|as.vector(RM)==1)){
#       grp <- grp + scale_fill_grey(name="Manhattan")
#     } else {
#       grp <- grp +  scale_fill_gradient2(low      = "red3",
#                                          high     = "steelblue",
#                                          mid      = "white",
#                                          na.value = scales::muted("slategray4"),
#                                          midpoint = .5,
#                                          limit    = c(0,1),
#                                          space    = "Lab",
#                                          name     = "Manhattan")
#     }
#
#     if(AUTO){
#       grp <- grp + labs(title = "Auto-Recurrence Plot\n")
#     } else {
#       grp <- grp + labs("Cross-Recurrence Plot\n")
#     }
#
#     # if(plotTime){
#     grp <- grp +
#       scale_x_discrete(limits = c(first(meltRP$Var1),last(meltRP$Var2))) + #breaks = meltRP$Var2[seq(1,nrow(meltRP),round(nrow(meltRP)/5))]) +
#       scale_y_discrete(limits = c(first(meltRP$Var1),last(meltRP$Var2))) #breaks = meltRP$Var1[seq(1,nrow(meltRP),round(nrow(meltRP)/5))])
#
#     grp <- grp + theme_bw() +  coord_fixed()
#
#       # theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 8, hjust = 1),
#       #       axis.text.y = element_text(size = 8)) +
#       #
#
#     if(plotVarsOK){
#       grp
#     }
#
#
#     if(plotMatrix){
#       grp
#       #grid.arrange(lp, txt, ncol=2, nrow=1, widths=c(4, 1), heights=c(4))
#     }
#
#     return(grp)
#   }
#
#   # } else {
#   #   warning("\nInput is not list output from function crqa().\n")
#   # }
#
# }


#' plotRP
#'
#' @param crqaOutput    List output from \code{\link[crqa]{crqa}}
#'
#' @return A list
#' @export
#' @author Fred Hasselman
#' @description Creates a recurrence plot from the sparse matrix output generated by \code{\link[crqa]{crqa}}.
plotRP.crqa <- function(crqaOutput){

  AUTO <- FALSE

  # Is is crqa output?
  if(all(names(crqaOutput)%in%c("RR","DET","NRLINE","maxL","L","ENTR","rENTR","LAM","TT","RP"))){

    # RP dimensions
    nr <- dim(crqaOutput$RP)[1]
    nc <- dim(crqaOutput$RP)[2]

    RP <- crqaOutput$RP

    # AUTO or CROSS?
    if(class(RP)=="dtCMatrix"){
      message("\nRP in crqa output is a Triangular Sparse Matrix, this implies auto-recurrence...\n")
      AUTO <- TRUE
    }
    if(isSymmetric(RP)){
      message("\nRP in crqa output is a Symmetric Sparse Matrix, this implies auto-recurrence: \n >> RP Measures will include the Line of Identity, unless:\n >> crqa() with side = 'upper' or 'lower' was used.\n >> crqa() with Theiler window (tw) of 1 was used.")
      AUTO <- TRUE}

    if(nc*nr>10^6){
      message(paste("\nLarge RP with",nc*nr,"elements... \n >> This could take some time to plot!\n"))
    }

    ifelse(AUTO,
           Titles <- list(main="Auto Recurrence Plot", X = expression(Y[t]), Y = expression(Y[t])),
           Titles <- list(main="Cross Recurrence Plot", X = expression(X[t]), Y = expression(Y[t]))
    )


    # Thresholded or Distance matrix?
    ifelse(is.logical(as.array(RP)),
           distPallette <- colorRampPalette(c("white", "black"))(2),
           distPallette <- colorRampPalette(c("red", "white", "blue")) #( length(unique(as.vector(RP))))
    )

    #distPallette <- colorRampPalette(c("white", "black"))(2)

    lp <- levelplot(Matrix::as.matrix(RP),
                    #colorkey = TRUE,
                    region = TRUE, col.regions = distPallette, useRaster = TRUE, aspect = nc/nr,
                    main = Titles$main,
                    xlab = Titles$X,
                    ylab = Titles$Y)

    txt <- grid.text(paste0("\n   RR: ", round(crqaOutput$RR, digits = 1),
                            "\n  DET: ", round(crqaOutput$DET, digits = 1),
                            "\nmeanL: ", round(crqaOutput$maxL, digits = 1),
                            "\n maxL: ", round(crqaOutput$NRLINE, digits = 1),
                            "\n  LAM: ", round(crqaOutput$LAM, digits = 1),
                            "\n   TT: ", round(crqaOutput$TT, digits = 1),
                            "\n  ENT: ", round(crqaOutput$ENT, digits = 1),
                            "\n rENT: ", round(crqaOutput$rENTR, digits = 1)),
                     .02,.7, just = "left", draw = FALSE, gp = gpar(fontfamily = "mono", cex = .7))
    # ,)

    grid.arrange(lp, txt, ncol=2, nrow=1, widths=c(4, 1), heights=c(4))

  } else {
    warning("\nInput is not list output from function crqa().\n")
  }

}

plotRP.fnn <- function(FNNoutput){
  plot(FNNoutput["combined",],type="b",pch=16, cex=2, col="grey80", ylim=c(0,100), xaxt="n",
       xlab = "Embedding Dimension", ylab = "False Nearest Neighbours")
  lines(FNNoutput["atol",],type="b",pch="a",col="grey30", lty=2)
  lines(FNNoutput["rtol",],type="b",pch="r", col="grey30",lty=2)
  Axis(side=1,at=seq_along(FNNoutput[1,]),labels = dimnames(FNNoutput)[[2]])
  legend("topright",c("Combined","atol","rtol"), pch = c(16,97,114), lty = c(1,2,2), col = c("grey80","grey30","grey30"), pt.cex=c(2,1,1))
}

# Toy models ------------------------------------------------------------------------------------------------------

#' Autocatlytic Growth: Iterating differential equations (maps)
#'
#' @title Autocatlytic Growth: Iterating differential equations (maps)
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
  return(ts(Y))
}

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
#' # The function such that it can take a set of conditional rules and apply them sequentially during the iterations.
#' # The conditional rules are passed as a `data.frame`
#' (cond <- cbind.data.frame(Y = c(0.2, 0.6), par = c("r", "r"), val = c(0.5, 0.1)))
#' xyplot(growth_ac_cond(cond=cond))
#'
#' # Combine a change of `r` and a change of `k`
#' (cond <- cbind.data.frame(Y = c(0.2, 1.99), par = c("r", "k"), val = c(0.5, 3)))
#' xyplot(growth_ac_cond(cond=cond))
#'
#' # A fantasy growth process
#' (cond <- cbind.data.frame(Y = c(0.1, 1.99, 1.999, 2.5, 2.9), par = c("r", "k", "r", "r","k"), val = c(0.3, 3, 0.9, 0.1, 1.3)))
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
  return(ts(Y))
}

# Complex Networks---------------------------------------------------------------
graph2svg <- function(TDM,pname){

  # Create weighted Term-Term matrix
  tTM <- as.matrix(TDM)
  TTM <- tTM %*% t(tTM)
  TTM <- log1p(TTM)

  g <- graph.adjacency(TTM,weighted=T,mode="undirected",diag=F)
  g <- simplify(g)

  # Remove vertices that were used in the search query
  Vrem <- which(V(g)$name %in% c("~dev~","~dys~","~sld~","development","children","dyslexia"))
  g <- (g - V(g)$name[Vrem])

  # Set colors and sizes for vertices
  V(g)$degree <- degree(g)
  rev         <- scaleRange(log1p(V(g)$degree))
  rev[rev<=0.3]<-0.3

  V(g)$color       <- rgb(scaleRange(V(g)$degree), 1-scaleRange(V(g)$degree),  0, rev)
  V(g)$size        <- 10*scaleRange(V(g)$degree)
  V(g)$frame.color <- NA

  # set vertex labels and their colors and sizes
  V(g)$label       <- V(g)$name
  V(g)$label.color <- rgb(0, 0, 0, rev)
  V(g)$label.cex   <- scaleRange(V(g)$degree)+.1

  # set edge width and color
  rew <- scaleRange(E(g)$weight)
  rew[rew<=0.3]<-0.3

  E(g)$width <- 2*scaleRange(E(g)$weight)
  E(g)$color <- rgb(.5, .5, 0, rew)
  set.seed(958)

  svg(paste(pname,sep=""),width=8,height=8)
  plot(g, layout=layout.fruchterman.reingold(g))
  dev.off()

  return(g)
}

# Plot vertex neighbourhood
hoodGraph2svg <- function(TDM,Vname,pname){

  # Create weighted Term-Term matrix
  tTM <- as.matrix(TDM)
  TTM <- tTM %*% t(tTM)
  TTM <- log1p(TTM)

  ig <- graph.adjacency(TTM,weighted=T,mode="undirected",diag=F)
  ig <- simplify(ig)

  # Remove vertices that were used in the search query
  Vrem <- which(V(ig)$name %in% c("~dev~","~dys~","~sld~","development","children","dyslexia"))
  ig <- (ig - V(ig)$name[Vrem])

  # This is a deletion specific for the Neighbourhood graphs
  Vrem <- which(V(ig)$name %in% c("~rdsp~","~imp~","~som~","~bod~","~mlt~"))
  ig   <- ig - V(ig)$name[Vrem]

  idx <- which(V(ig)$name==Vname)
  sg  <- graph.neighborhood(ig, order = 1, nodes=V(ig)[idx], mode = 'all')[[1]]

  # set colors and sizes for vertices
  V(sg)$degree <- degree(sg)

  rev<-scaleRange(log1p(V(sg)$degree))
  rev[rev<=0.3]<-0.3

  V(sg)$color <- rgb(scaleRange(V(sg)$degree), 1-scaleRange(log1p(V(sg)$degree*V(sg)$degree)),  0, rev)

  V(sg)$size        <- 35*scaleRange(V(sg)$degree)
  V(sg)$frame.color <- NA

  # set vertex labels and their colors and sizes
  V(sg)$label       <- V(sg)$name
  V(sg)$label.color <- rgb(0, 0, 0, rev)
  V(sg)$label.cex   <- scaleRange(V(sg)$degree)

  # set edge width and color
  rew<-scaleRange(E(sg)$weight)
  rew[rew<=0.3]<-0.3

  E(sg)$width <- 6*scaleRange(E(sg)$weight)
  E(sg)$color <- rgb(.5, .5, 0, rew)

  idV <- which(V(sg)$name==Vname)
  idE <- incident(sg,V(sg)[[idV]])
  E(sg)$color[idE] <- rgb(0, 0, 1 ,0.8)

  set.seed(958)

  idx <- which(V(sg)$name==Vname)
  svg(paste(pname,sep=""),width=8,height=8)
  plot(sg,layout=layout.star(sg,center=V(sg)[idx]))
  dev.off()

  return(sg)
}

#' lag Embed a time series
#'
#' @param y Time series
#' @param emDim Embedding dimension
#' @param emLag Embedding lag
#'
#' @return The lag embedded time series
#' @export
#'
lagEmbed <- function (y, emDim, emLag){

  if(any(is.ts(y), zoo::is.zoo(y), xts::is.xts(y))){
    timeVec <- time(y)
    emTime <- lubridate::as_datetime(timeVec[emLag+1])- lubridate::as_datetime(timeVec[1])
  } else {
    timeVec <- zoo::index(y)
    emTime <- emLag
  }

  y <- as.numeric(unlist(y))
  N  <- NROW(y)
  if(emDim > 1){
    lag.id    <- seq(1, (emDim*emLag), emLag)
    maxN      <- N - max(lag.id)
    emY       <- matrix(nrow = maxN, ncol = emDim)
    for(tau in seq_along(lag.id)){
      emY[,tau] = y[lag.id[tau]:(N-(rev(lag.id)[tau]))]
    }
    colnames(emY) <- paste0("tau.",0:(emDim-1))
  } else {
    emY <-  matrix(nrow = N, ncol = 1, byrow = TRUE, dimnames = list(NULL,"tau.0"))
    emY <- y
  }

  id <- deparse(substitute(y))
  attr(emY, "embedding.dims") = emDim
  attr(emY, "embedding.lag")  = emLag
  attr(emY, "embedding.time") = emTime
  attr(emY, "id") = id
  return(as.matrix(emY))
}


sliceTS<-function(TSmat,epochSz=1) {
  # Slice columns of TSmat in epochs of size = epochSz
  require(plyr)

  N<-dim(TSmat)
  return(llply(seq(1,N[1],epochSz),function(i) TSmat[i:min(i+epochSz-1,N[1]),1:N[2]]))
}

fltrIT <- function(TS,f){
  # Apply filtfilt to TS using f (filter settings)
  require("signal")

  return(filtfilt(f=f,x=TS))

}

SWtest0 <- function(g){
  Nreps <- 10;
  histr  <- vector("integer",Nreps)
  target<- round(mean(degree(g)))
  now   <- target/2
  for(i in 1:Nreps){
    gt      <- watts.strogatz.game(dim=1, size=length(degree(g)), nei=now, 0)
    histr[i] <- round(mean(degree(gt)))
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
#  return(list(cp=transitivity(g,type="global"),cpR=transitivity(rewire(g,mode=c("simple"),niter=N),type="global"),lp=average.path.length(g), lpR=average.path.length(rewire(g,mode=c("simple"),niter=N))))
# }

SWtestE <- function(g,p=1,N=20){
  values <- matrix(nrow=N,ncol=6,dimnames=list(c(1:N),c("cp","cpR","cp0","lp","lpR","lp0")))

  for(n in 1:N) {
    gt<-SWtest0(g)
    values[n,] <- c(transitivity(g,type="localaverage"),transitivity(rewire(g,each_edge(p=p)),type="localaverage"),transitivity(gt,type="localaverage"),average.path.length(g),average.path.length(rewire(g,each_edge(p=p))),average.path.length(gt))}
  values[n,values[n,]==0] <- NA #values[n,values[n,]==0]+1e-8}

  values   <- cbind(values,(values[,1]/values[,2])/(values[,4]/values[,5]),(values[,1]/values[,3]),(values[,4]/values[,6]),((values[,1]/values[,3])/values[,2])/((values[,4]/values[,6])/values[,5]))
  valuesSD <- data.frame(matrix(apply(values[,1:10],2,sd,na.rm=T),nrow=1,ncol=10,dimnames=list(c(1),c("cp","cpR","cp0","lp","lpR","lp0","SWI","cp:cp0","lp:lp0","SWIn"))))
  valuesAV <- data.frame(matrix(colMeans(values[,1:10],na.rm=T),nrow=1,ncol=10,dimnames=list(c(1),c("cp","cpR","cp0","lp","lpR","lp0","SWI","cp:cp0","lp:lp0","SWIn"))))
  return(list(valuesAV=valuesAV,valuesSD=valuesSD,valuesSE=valuesSD/sqrt(N)))
}

PLFsmall <- function(g){

  if(length(V(g))>100){stop("Vertices > 100, no need to use PLFsmall, use a binning procedure")}

  d <- degree(g)

  y <- hist(d,breaks=0.5:(max(d)+0.5),plot=FALSE)$counts
  if(length(y)<2){
    warning("Less than 2 points in Log-Log regression... alpha=0")
    alpha <- 0
  } else {
    if(length(y)==2){
      warning("Caution... Log-Log slope is a bridge (2 points)")
      chop <- 0
    } else {
      chop <- 1
    }
    alpha <- coef(lm(rev(log1p(y)[1:(length(y)-chop)]) ~ log1p(1:(length(y)-chop))))[2]
  }

  return(alpha)
}

plotSW <- function(n,k,p){

  g <- watts.strogatz.game(1, n, k, p)

  V(g)$degree <- degree(g)

  # set colors and sizes for vertices
  rev<-elascer(log1p(V(g)$degree))
  rev[rev<=0.2]<-0.2
  rev[rev>=0.9]<-0.9
  V(g)$rev <- rev$x

  V(g)$color       <- rgb(V(g)$rev, 1-V(g)$rev,  0, 1)
  V(g)$size        <- 25*V(g)$rev

  # set vertex labels and their colors and sizes
  V(g)$label       <- ""

  E(g)$width <- 1
  E(g)$color <- rgb(0.5, 0.5, 0.5, 1)

  return(g)
}

plotBA <- function(n,pwr,out.dist){
  #require("Cairo")

  g <- barabasi.game(n,pwr,out.dist=out.dist,directed=FALSE)
  V(g)$degree <- degree(g)

  # set colors and sizes for vertices
  rev<-elascer(log1p(V(g)$degree))
  rev[rev<=0.2] <- 0.2
  rev[rev>=0.9] <- 0.9
  V(g)$rev <- rev$x

  V(g)$color    <- rgb(V(g)$rev, 1-V(g)$rev,  0, 1)
  V(g)$size     <- 25*V(g)$rev
  # V(g)$frame.color <- rgb(.5, .5,  0, .4)

  # set vertex labels and their colors and sizes
  V(g)$label <- ""

  E(g)$width <- 1
  E(g)$color <- rgb(0.5, 0.5, 0.5, 1)

  return(g)
}


FDrel <- function(g){
  d<-degree(g,mode="all")
  nbreaks <- round(length(V(g))/2)-1
  y<-hist(d,breaks=nbreaks,plot=F)$density
  y<-y[y>0]
  return(FD <- -sum(y*log2(y))/-(log2(1/length(y))))
}

sa2fd <- function(sa, ...) UseMethod("sa2fd")

#' sa2df.default
#'
#' @param sa Self-affinity parameter
#' @param ... Other argumentd
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @family SA to FD converters
#' @keywords internal
#' @export
#'
sa2fd.default <- function(sa, ...){
  cat("No type specified.")
}


#' Informed Dimension estimate from Spectral Slope (aplha)
#'
#' @description Conversion formula: From periodogram based self-affinity parameter estimate (\code{sa}) to an informed estimate of the (fractal) dimension (FD).
#' @param sa Self-Affinity parameter estimate based on PSD slope (e.g., \code{\link{fd.psd}})).
#'
#' @return An informed estimate of the Fractal Dimension, see Hasselman(2013) for details.
#' @export
#'
#' @details The spectral slope will be converted to a dimension estimate using:
#'
#' \deqn{D_{PSD}\approx\frac{3}{2}+\frac{14}{33}*\tanh\left(Slope * \ln(1+\sqrt{2})\right)}
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @family SA to FD converters
#'
#' @examples
#' # Informed FD of white noise
#' sa2fd.psd(0)
#'
#' # Informed FD of Brownian noise
#' sa2fd.psd(-2)
#'
#' # Informed FD of blue noise
#' sa2fd.psd(2)
sa2fd.psd <- function(sa){return(round(3/2 + ((14/33)*tanh(sa*log(1+sqrt(2)))), digits = 2))}


#' Informed Dimension estimate from DFA slope (H)
#'
#' @description Conversion formula: Detrended Fluctuation Analysis (DFA) estimate of the Hurst exponent (a self-affinity parameter \code{sa}) to an informed estimate of the (fractal) dimension (FD).
#'
#' @param sa Self-Afinity parameter estimate based on DFA slope (e.g., \code{\link{fd.sda}})).
#'
#' @return An informed estimate of the Fractal Dimension, see Hasselman(2013) for details.
#'
#' @export
#'
#' @details The DFA slope (H) will be converted to a dimension estimate using:
#'
#' \deqn{D_{DFA}\approx 2-(\tanh(\log(3)*sa)) }{D_{DFA} ≈ 2-(tanh(log(3)*sa)) }
#'
#' @family SA to FD converters
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @examples
#' # Informed FD of white noise
#' sa2fd.dfa(0.5)
#'
#' # Informed FD of Pink noise
#' sa2fd.dfa(1)
#'
#' # Informed FD of blue noise
#' sa2fd.dfa(0.1)
sa2fd.dfa <- function(sa){return(round(2-(tanh(log(3)*sa)), digits = 2))}


#' Informed Dimension estimate from SDA slope.
#'
#' @description Conversion formula: Standardised Dispersion Analysis (SDA) estimate of self-affinity parameter (\code{SA}) to an informed estimate of the fractal dimension (FD).
#'
#' @param sa Self-afinity parameter estimate based on SDA slope (e.g., \code{\link{fd.sda}})).
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
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @family SA to FD converters
#'
#' @examples
#' # Informed FD of white noise
#' sa2fd.sda(-0.5)
#'
#' # Informed FD of Brownian noise
#' sa2fd.sda(-1)
#'
#' # Informed FD of blue noise
#' sa2fd.sda(-0.9)
sa2fd.sda <- function(sa){return(1-sa)}



# fd estimators ----------------------------------------------

fd <- function(y, ...) UseMethod("fd")

#' fd.default
#'
#' @param y timeseries
#' @param ... Other arguments
#'
#' @keywords internal
#' @family FD estimators
#' @export
#'
fd.default <- function(y, ...){

  cat("No type specified.\nReturning exponential growth power law.")

  r = 1.01
  y <- growth_ac(Y0=0.001, r=r, N=2048, type = "driving")
  tsp(y) <-c(1/500,2048/500,500)
  bulk <- log1p(hist(y,plot = F, breaks = seq(0,max(y),length.out = 129))$counts)
  size <- log1p(seq(0,2047,length.out = 128))
  id<-bulk==0

  lmfit <- lm(bulk[!id] ~ size[!id])

  old <- ifultools::splitplot(2,1,1)
  plot(y, ylab = "Y", main = paste0('Exponential growth  sap: ', round(coef(lmfit)[2],digits=2), ' | r:', r))
  ifultools::splitplot(2,1,2)
  plot(size[!id],bulk[!id], xlab="Size = log(bin(Time))", ylab = "Bulk = logbin(Y)", pch=21, bg="grey60", pty="s")
  lines(size[!id], predict(lmfit),lwd=4,col="darkred")
  #legend("bottomleft",c(paste0("Range (n = ",sum(powspec$size<=0.25),")"), paste0("Hurvic-Deo estimate (n = ",nr,")")), lwd=c(3,3),col=c("darkred","darkblue"), cex = .8)
  par(old)
}


# PSD -------------------------------------------------------------------------------------------------------------

#' @title Power Spectral Density Slope (PSD).

#' @description Estimate Alpha, Hurst Exponent and Fractal Dimension through log-log slope.
#'
#' @param y    A numeric vector or time series object.
#' @param normalize    Normalize the series (default).
#' @param detrend    Subtract linear trend from the series (default).
#' @param plot    Return the log-log spectrum with linear fit (default).
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @return A list object containing:
#' \itemize{
#' \item A data matrix \code{PLAW} with columns \code{freq.norm}, \code{size} and \code{bulk}.
#' \item Estimate of scaling exponent \code{alpha} based on a fit over the lowest 25\% frequencies (\code{low25}), or using the HD estimate \code{HD}.
#' \item Estimate of the the Fractal Dimension (\code{FD}) using conversion formula's reported in Hasselman(2013).
#' \item Information output by various functions.
#' }
#'
#' @family FD estimators
#'
#' @export
#'
#' @details Calls function \code{\link[sapa]{SDF}} to estimate the scaling exponent of a timeseries based on the periodogram frequency spectrum. After detrending and normalizing the signal (if requested), \code{SDF} is called using a Tukey window (\code{raised cosine \link[sapa]{taper}}).
#'
#' A line is fitted on the periodogram in log-log coordinates. Two fit-ranges are used: The 25\% lowest frequencies and the Hurvich-Deo estimate (\code{\link[fractal]{HDEst}}).
#'
fd.psd <- function(y, fs = NULL, normalize = TRUE, dtrend = TRUE, plot = FALSE){
  # require(pracma)
  # require(fractal)
  # require(sapa)
  # require(ifultools)

  if(!is.ts(y)){
    if(is.null(fs)){fs <- 1}
    y <- ts(y, frequency = fs)
    cat("\n\nfd.psd:\tSample rate was set to 1.\n\n")
  }

  N             <- length(y)
  # Simple linear detrending.
  if(dtrend)    y <- ts(pracma::detrend(as.vector(y), tt = 'linear'), frequency = fs)
  # Normalize using N instead of N-1.
  if(normalize) y <- (y - mean(y, na.rm = TRUE)) / (sd(y, na.rm = TRUE)*sqrt((N-1)/N))

  # Number of frequencies estimated cannot be set! (defaults to Nyquist)
  # Use Tukey window: cosine taper with r = 0.5

  # fast = TRUE ensures padding with zeros to optimize FFT to highly composite number.
  # However, we just pad to nextPow2, except if length already is a power of 2.
  npad <- 1+(stats::nextn(N,factors=2)-N)/N
  npad <- stats::nextn(N)

  # if(N==npad) npad = 0
  # psd  <- stats::spec.pgram(y, fast = FALSE, demean=FALSE, detrend=FALSE, plot=FALSE, pad=npad, taper=0.5)

  Tukey <- sapa::taper(type="raised cosine", flatness = 0.5, n.sample = npad)
  psd   <- sapa::SDF(y, taper. = Tukey, npad = npad)

  powspec <- cbind.data.frame(freq.norm = attr(psd, "frequency")[-1], size = attr(psd, "frequency")[-1]*frequency(y), bulk = as.matrix(psd)[-1])

  # First check the global slope for anti-persistent noise (GT +0.20)
  # If so, fit the line starting from the highest frequency
  nr     <- length(powspec[,1])
  lsfit  <- lm(log(powspec$bulk[1:nr]) ~ log(powspec$size[1:nr]))
  glob   <- coef(lsfit)[2]

  # General guideline: fit over 25% frequencies
  # If signal is continuous (sampled) consider Wijnants et al. (2013) log-log fitting procedure
  nr <- fractal::HDEst(NFT = length(powspec[,1]), sdf = psd)

  exp1 <- fractal::hurstSpec(y, sdf.method = "direct", freq.max = 0.25, taper. = Tukey )
  exp2 <- fractal::hurstSpec(y, sdf.method = "direct", freq.max = powspec$freq.norm[nr], taper. = Tukey)

  ifelse((glob > 0.2), {
    lmfit1 <- stats::lm(log(rev(powspec$bulk[powspec$size<=0.25])) ~ log(rev(powspec$size[powspec$size<=0.25])))
    lmfit2 <- stats::lm(log(rev(powspec$bulk[1:nr])) ~ log(rev(powspec$size[1:nr])))
  },{
    lmfit1 <- stats::lm(log(powspec$bulk[powspec$size<=0.25]) ~ log(powspec$size[powspec$size<=0.25]))
    lmfit2 <- stats::lm(log(powspec$bulk[1:nr]) ~ log(powspec$size[1:nr]))
  })

  if(plot){
    old<- ifultools::splitplot(2,1,1)
    plot(y,ylab = "Y", main = paste0('Lowest 25%    sap: ', round(coef(lmfit1)[2],digits=2), ' | H:', round(exp1,digits=2), ' | FD:',round(sa2fd.psd(coef(lmfit1)[2]),digits=2),'\nHurvic-Deo    sap: ', round(coef(lmfit2)[2],digits=2), ' | H:', round(exp2,digits=2), ' | FD:',round(sa2fd.psd(coef(lmfit2)[2]),digits=2)))
    ifultools::splitplot(2,1,2)
    plot(log(powspec$bulk) ~ log(powspec$size), xlab="log(Frequency)", ylab = "log(Power)")
    lines(log(powspec$size[powspec$size<=0.25]), predict(lmfit1),lwd=3,col="darkred")
    lines(log(powspec$size[1:nr]), predict(lmfit2),lwd=3,col="darkblue")
    legend("bottomleft",c(paste0("lowest 25% (n = ",sum(powspec$size<=0.25),")"), paste0("Hurvic-Deo estimate (n = ",nr,")")), lwd=c(3,3),col=c("darkred","darkblue"), cex = .8)
    par(old)
  }

  return(list(
    PLAW  = powspec,
    low25 = list(sap = coef(lmfit1)[2], H = exp1, FD = sa2fd.psd(coef(lmfit1)[2]), fitlm1 = lmfit1),
    HD    = list(sap = coef(lmfit2)[2], H = exp2, FD = sa2fd.psd(coef(lmfit2)[2]), fitlm2 = lmfit2),
    info  = psd)
  )
}


# SDA -------------------------------------------------

#' fd.sda
#'
#' @title Standardised Dispersion Analysis (SDA).
#'
#' @param y    A numeric vector or time series object.
#' @param normalize    Normalize the series (default).
#' @param plot    Return the log-log spectrum with linear fit (default).
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @return A list object containing:
#' \itemize{
#' \item A data matrix \code{PLAW} with columns \code{freq.norm}, \code{size} and \code{bulk}.
#' \item Estimate of scaling exponent \code{sap} based on a fit over the standard range (\code{fullRange}), or on a user defined range \code{fitRange}.
#' \item Estimate of the the Fractal Dimension (\code{FD}) using conversion formula's reported in Hasselman(2013).
#' \item Information output by various functions.
#' }
#'
#' @export
#'
#' @family FD estimators
#'
fd.sda <- function(y, fs = NULL, normalize = TRUE, dtrend = FALSE, scales = dispersion(y)$scale, fitRange = c(scales[1], scales[length(scales)-2]), plot = FALSE){
  require(pracma)
  require(fractal)

  if(!is.ts(y)){
    if(is.null(fs)){fs <- 1}
    y <- ts(y, frequency = fs)
    cat("\n\nfd.sda:\tSample rate was set to 1.\n\n")
  }

  N             <- length(y)
  # Simple linear detrending.
  if(dtrend)    y <- ts(pracma::detrend(as.vector(y), tt = 'linear'), frequency = fs)
  # Normalize using N instead of N-1.
  if(normalize) y <- (y - mean(y, na.rm = TRUE)) / (sd(y, na.rm = TRUE)*sqrt((N-1)/N))

  bins          <- which(fitRange[1]==scales):which(fitRange[2]==scales)
  out           <- dispersion(y, front = FALSE)
  lmfit1        <- lm(log(out$sd) ~ log(out$scale))
  lmfit2        <- lm(log(out$sd[bins]) ~ log(out$scale[bins]))

  if(plot){
    old<- ifultools::splitplot(2,1,1)
    plot(y,ylab = "Y", main = paste0('Full    sap: ', round(coef(lmfit1)[2],digits=2), ' | H:', round(1+coef(lmfit1)[2],digits=2), ' | FD:',round(sa2fd.sda(coef(lmfit1)[2]),digits=2),'\nRange    sap: ', round(coef(lmfit2)[2],digits=2), ' | H:', round(1+coef(lmfit1)[2],digits=2), ' | FD:',round(sa2fd.sda(coef(lmfit2)[2]),digits=2)))
    ifultools::splitplot(2,1,2)
    plot(log(out$sd) ~ log(out$scale), xlab="log(Bin Size)", ylab = "log(SD)")
    lines(log(out$scale), predict(lmfit1),lwd=3,col="darkred")
    lines(log(out$scale[bins]), predict(lmfit2),lwd=3,col="darkblue")
    legend("bottomleft",c(paste0("Full (n = ",length(out$scale),")"), paste0("Range (n = ",length(bins),")")), lwd=c(3,3),col=c("darkred","darkblue"), cex = .8)
    par(old)
  }

  return(list(
    PLAW  =  cbind.data.frame(freq.norm = frequency(y)/scales, size = out$scale, bulk = out$sd),
    fullRange = list(sap = coef(lmfit1)[2], H = 1+coef(lmfit1)[2], FD = sa2fd.sda(coef(lmfit1)[2]), fitlm1 = lmfit1),
    fitRange  = list(sap = coef(lmfit2)[2], H = 1+coef(lmfit2)[2], FD = sa2fd.sda(coef(lmfit2)[2]), fitlm2 = lmfit2),
    info = out)
  )
}


# DFA ---------------------------------------------

#' fd.dfa
#'
#' @title Detrended Fluctuation Analysis (DFA)
#'
#' @param y    A numeric vector or time series object.
#' @param normalize    Normalize the series (default).
#' @param detrend    Subtract linear trend from the series (default).
#' @param dmethod     Method to use for detrending, see \code{\link[fractal]{DFA}}.
#' @param plot    Return the log-log spectrum with linear fit (default).
#'
#'
#' @return Estimate of Hurst exponent (slope of \code{log(bin)} vs. \code{log(RMSE))} and an FD estimate based on Hasselman(2013)
#' A list object containing:
#' \itemize{
#' \item A data matrix \code{PLAW} with columns \code{freq.norm}, \code{size} and \code{bulk}.
#' \item Estimate of scaling exponent \code{sap} based on a fit over the standard range (\code{fullRange}), or on a user defined range \code{fitRange}.
#' \item Estimate of the the Fractal Dimension (\code{FD}) using conversion formula's reported in Hasselman(2013).
#' \item Information output by various functions.
#' }
#'
#' @export
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. \url{http://doi.org/10.3389/fphys.2013.00075}
#'
#' @family FD estimators
#'
fd.dfa <- function(y, fs = NULL, dtrend = "poly1", normalize = FALSE, sum.order = 1, scale.max=trunc(length(y)/4), scale.min=4, elasceratio=2^(1/4), overlap = 0, plot = FALSE){
  require(pracma)
  require(fractal)

  reload <- FALSE
  if("signal" %in% .packages()){
    warning("signal:poly is loaded and stats:poly is needed... will unload package:signal, compute slope, and reload...")
    reload <- TRUE
    detach("package:signal", unload=TRUE)
  }



  if(!is.ts(y)){
    if(is.null(fs)){fs <- 1}
    y <- ts(y, frequency = fs)
    cat("\n\nfd.dfa:\tSample rate was set to 1.\n\n")
  }

  N             <- length(y)
  # Normalize using N instead of N-1.
  if(normalize) y <- (y - mean(y, na.rm = TRUE)) / (sd(y, na.rm = TRUE)*sqrt((N-1)/N))

  out1 <- fractal::DFA(y, detrend=dtrend, sum.order=sum.order, scale.max=trunc(length(y)/2), scale.min=2, elasceratio=2, overlap = 0, verbose=FALSE)
  out2 <- fractal::DFA(y, detrend=dtrend, sum.order=sum.order, scale.max=scale.max, scale.min=scale.min, elasceratio=elasceratio, overlap = overlap, verbose=FALSE)

  lmfit1        <- lm(log(attributes(out1)$stat) ~ log(attributes(out1)$scale))
  lmfit2        <- lm(log(attributes(out2)$stat) ~ log(attributes(out2)$scale))

  if(plot){
    plot.new()
    old <- ifultools::splitplot(2,1,1)
    plot(y,ylab = "Y", main = paste0('Full    sap: ', round(coef(lmfit1)[2],digits=2), ' | H:',
                                     round(attributes(out1)$logfit[]$coefficients['x'] ,digits=2), ' | FD:',
                                     round(sa2fd.dfa(coef(lmfit1)[2]),digits=2),'\nRange    sap: ',
                                     round(coef(lmfit2)[2],digits=2), ' | H:',
                                     round(attributes(out2)$logfit[]$coefficients['x'] ,digits=2), ' | FD:',
                                     round(sa2fd.dfa(coef(lmfit2)[2]),digits=2)
    )
    )
    ifultools::splitplot(2,1,2)
    plot(log(attributes(out1)$stat) ~ log(attributes(out1)$scale), xlab="log(Bin Size)", ylab = "log(RMSE)")
    lines(log(attributes(out1)$scale), predict(lmfit1),lwd=3,col="darkred")
    lines(log(attributes(out2)$scale), predict(lmfit2),lwd=3,col="darkblue")
    legend("topleft",c(paste0("Full (n = ",length(attributes(out1)$scale),")"), paste0("Range (n = ",length(attributes(out2)$scale),")")), lwd=c(3,3),col=c("darkred","darkblue"), cex = .8)
    par(old)
  }

  if(reload==TRUE){library(signal,verbose=FALSE,quietly=TRUE)}

  return(list(
    PLAW  =  cbind.data.frame(freq.norm = elascer(attributes(out1)$scale*frequency(y)), size = attributes(out1)$scale, bulk = attributes(out1)$stat),
    fullRange = list(sap = coef(lmfit1)[2], H = attributes(out1)$logfit[]$coefficients['x'] , FD = sa2fd.dfa(coef(lmfit1)[2]), fitlm1 = lmfit1),
    fitRange  = list(sap = coef(lmfit2)[2], H = coef(lmfit2)[2], FD = sa2fd.dfa(coef(lmfit2)[2]), fitlm2 = lmfit2),
    info = list(out1,out2))
  )
}



#' Detrended Fluctuation Analysis
#'
#' @param signal    An input signal.
#' @param mins    Minimum scale to consider.
#' @param maxs    Maximum scale to consider.
#' @param ressc ressc
#' @param m m
#'
#' @return  output
#' @export
#'
DFA1 <- function(signal,mins=4,maxs=12,ressc=30,m=1){

  #   reload <- FALSE
  #   if("signal" %in% .packages()){
  #     warning("signal:poly is loaded and stats:poly is needed... will unload package:signal, compute slope, and reload...")
  #     reload <- TRUE
  #     detach("package:signal", unload=TRUE)
  #   }
  scale     <- round(2^(seq(mins,maxs,by=((maxs-mins)/ressc))))
  segv      <- numeric(length(scale))
  RMS_scale <- vector("list",length(scale))
  # qRMS      <- vector("list",length(qq))
  # Fq        <- vector("list",length(qq))
  # qRegLine  <- vector("list",length(qq))
  # Hq        <- numeric(length(qq))

  Y        <- cumsum(signal-mean(signal))
  TSm      <- as.matrix(cbind(t=1:length(Y),y=Y))
  Hglobal  <- monoH(TSm,scale)
}


# Multi-Fractal DFA -----------------------------------------------------------------------------------------------

#' Multi-fractal Detrended Fluctuation Analysis
#'
#' @param signal    An input signal.
#' @param qq    A vector containing a range of values for the order of fluctuation \code{q}.
#' @param mins    Minimum scale to consider.
#' @param maxs    Maximum scale to consider.
#' @param ressc rescc
#' @param m m
#'
#' @return output
#' @export
#'
MFDFA <- function(signal,qq=c(-10,-5:5,10),mins=6,maxs=12,ressc=30,m=1){

  #   reload <- FALSE
  #   if("signal" %in% .packages()){
  #     warning("signal:poly is loaded and stats:poly is needed... will unload package:signal, compute slope, and reload...")
  #     reload <- TRUE
  #     detach("package:signal", unload=TRUE)
  #   }
  scale     <- round(2^(seq(mins,maxs,by=((maxs-mins)/ressc))))
  segv      <- numeric(length(scale))
  RMS_scale <- vector("list",length(scale))
  qRMS      <- vector("list",length(qq))
  Fq        <- vector("list",length(qq))
  qRegLine  <- vector("list",length(qq))
  Hq        <- numeric(length(qq))

  Y        <- cumsum(signal-mean(signal))
  TSm      <- as.matrix(cbind(t=1:length(Y),y=Y))
  Hglobal  <- monoH(TSm,scale)

  Hadj <- 0
  if((Hglobal>1.2)&(Hglobal<1.8)){
    Y <- diff(signal)
    Hadj=1}
  if(Hglobal>1.8){
    Y <- diff(diff(signal))
    Hadj <- 2}
  if(Hglobal<0.2){
    Y <- cumsum(signal-mean(signal))
    Hadj <- -1}
  if(Hadj!=0){TSm  <- as.matrix(cbind(t=1:length(Y),y=cumsum(Y-mean(Y))))}

  for(ns in seq_along(scale)){
    RMS_scale[[ns]] <- ldply(sliceTS(TSm,scale[ns]),function(sv){return(sqrt(mean(detRend(sv[,2]))^2))})
    for(nq in seq_along(qq)){
      qRMS[[nq]][1:length(RMS_scale[[ns]]$V1)] <- RMS_scale[[ns]]$V1^qq[nq]
      Fq[[nq]][ns] <- mean(qRMS[[nq]][1:length(RMS_scale[[ns]]$V1)])^(1/qq[nq])
      if(is.inf(log2(Fq[[nq]][ns]))){Fq[[nq]][ns]<-NA}
    }
    Fq[[which(qq==0)]][ns] <- exp(0.5*mean(log(RMS_scale[[ns]]^2)))
    if(is.inf(log2(Fq[[which(qq==0)]][ns]))){Fq[[which(qq==0)]][ns]<-NA}
  }

  fmin<-1
  fmax<-which(scale==max(scale))
  #for(nq in seq_along(qq)){Hq[nq] <- lm(log2(Fq[[nq]])~log2(scale))$coefficients[2]}
  Hq <- ldply(Fq,function(Fqs){lm(log2(Fqs[fmin:fmax])~log2(scale[fmin:fmax]),na.action=na.omit)$coefficients[2]})

  tq <- (Hq[,1]*qq)-1
  hq <- diff(tq)/diff(qq)
  Dq <- (qq[1:(length(qq)-1)]*hq) - (tq[1:(length(qq)-1)])

  if(reload==TRUE){library(signal,verbose=FALSE,quietly=TRUE)}

  return(list(q=qq,Hq=Hq,tq=tq,hq=hq,Dq=Dq,Hglobal=Hglobal,Hadj=Hadj))
}


monoH <- function(TSm,scale){
  dfaRMS_scale <- vector("list",length(scale))
  F2 <- numeric(length(scale))
  for(ns in seq_along(scale)){
    dfaRMS_scale[[ns]] <- ldply(sliceTS(TSm,scale[ns]),function(sv){return(sqrt(mean(detRend(sv[,2]))^2))})
    F2[ns] <- mean(dfaRMS_scale[[ns]]$V1^2)^(1/2)
    if(is.inf(log2(F2[ns]))){F2[ns] <- NA}
  }
  return(lm(log2(F2)~log2(scale),na.action=na.omit)$coefficients[2])
}

detRend <- function(TS, Order=1){
  detR <- lm(TS~stats::poly(1:length(TS), degree=Order))$residuals
  return(detR)
}


#
# set.seed(100)
# z <- dispersion(rnorm(1024))
# plot(log(z$scale),log(z$sd))
# #

# trace(detRend,edit=T)
# seq(1,length(X),by=4096)
#
# z<-sliceTS(TSm,scale[1])
# z[[1]][,2]
#
# Hglobal <-
#
# segments <- laply(scale,function(s) floor(length(X)/s))
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
# plot(hq,Dq,type="l")
#
#
# qq<-c(-10,-5,seq(-2,2,.1),5,10)


# PLOTS -------------------------------------------------------------------
#
#
#' gg.theme
#'
#' @param type      One of \code{"clean"}, or \code{"noax"}
#' @param useArial    Use the Arial font (requires \code{.afm} font files in the \code{afmPath})
#' @param afmPATH    Path to Arial \code{.afm} font files.
#'
#' @details Will generate a \code{"clean"} ggplot theme, or a theme without any axes (\code{"noax"}).
#'
#' Some scientific journals explicitly request the Arial font should be used in figures. This can be achieved by using \code{.afm} font format (see, e.g. http://www.pure-mac.com/font.html).
#'
#' @return A theme for \code{ggplot2}.
#' @export
#'
#' @examples
#' library(ggplot2)
#' g <- ggplot(data.frame(x = rnorm(n = 100), y = rnorm(n = 100)), aes(x = x, y = y)) + geom_point()
#' g + gg.theme()
#' g + gg.theme("noax")
gg.theme <- function(type=c("clean","noax"),useArial = F, afmPATH="~/Dropbox"){

  if(length(type)>1){type <- type[1]}

  if(useArial){
    set.Arial(afmPATH)
    bf_font="Arial"
  } else {bf_font="Helvetica"}

  switch(type,
         clean = theme_bw(base_size = 16, base_family=bf_font) +
           theme(axis.text.x     = element_text(size = 14),
                 axis.title.y    = element_text(vjust = +1.5),
                 panel.grid.major  = element_blank(),
                 panel.grid.minor  = element_blank(),
                 legend.background = element_blank(),
                 legend.key = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 axis.line  = element_line(colour = "black")),
         noax = theme(line = element_blank(),
                      text  = element_blank(),
                      title = element_blank(),
                      plot.background = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank())
  )
}

#' gg.plotHolder
#'
#' @param useArial    Use the Arial font (requires \code{.afm} font files in the \code{afmPath})
#' @param afmPATH    Path to Arial \code{.afm} font files.
#'
#' @return A blank \code{ggplot2} object that can be used in concordance with \code{grid.arrange}.
#' @export
#'
#' @examples
#' # Create a plot with marginal distributions.
#' library(ggplot2)
#' library(scales)
#'
#' df <- data.frame(x = rnorm(n = 100), y = rnorm(n = 100), group = factor(sample(x=c(0,1), size = 100, replace = TRUE)))
#'
#' scatterP <- ggplot(df, aes(x = x, y =y, colour = group)) + geom_point() + gg.theme()
#' xDense <- ggplot(df, aes(x = x, fill = group)) + geom_density(aes(y= ..count..),trim=FALSE, alpha=.5) + gg.theme("noax") + theme(legend.position = "none")
#' yDense <- ggplot(df, aes(x = y, fill = group)) + geom_density(aes(y= ..count..),trim=FALSE, alpha=.5) + coord_flip() + gg.theme("noax") + theme(legend.position = "none")
#'
#' library(gridExtra)
#' grid.arrange(xDense, gg.plotHolder(), scatterP, yDense, ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
gg.plotHolder <- function(useArial = F,afmPATH="~/Dropbox"){
  ggplot() +
    geom_blank(aes(1,1)) +
    theme(line = element_blank(),
          text  = element_blank(),
          title = element_blank(),
          plot.background = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()
    )
}

set.Arial <- function(afmPATH="~/Dropbox"){
  # Set up PDF device on MAC OSX to use Arial as a font in Graphs
  if(nchar(afmPATH>0)){
    if(file.exists(paste0(afmPATH,"/Arial.afm"))){
      Arial <- Type1Font("Arial",
                         c(paste(afmPATH,"/Arial.afm",sep=""),
                           paste(afmPATH,"/Arial Bold.afm",sep=""),
                           paste(afmPATH,"/Arial Italic.afm",sep=""),
                           paste(afmPATH,"/Arial Bold Italic.afm",sep="")))
      if(!"Arial" %in% names(pdfFonts())){pdfFonts(Arial=Arial)}
      if(!"Arial" %in% names(postscriptFonts())){postscriptFonts(Arial=Arial)}
      return()
    } else {disp(header='useArial=TRUE',message='The directory did not contain the *.afm version of the Arial font family')}
  } else {disp(header='useArial=TRUE',message='Please provide the path to the *.afm version of the Arial font family')}
}


netGroupCol <- function(g,grp){
  if(length(grp)<=12){
    groupColours <-  brewer_pal(palette="Set3")(length(grp))
  } else {
    groupColours <- gradient_n_pal(brewer_pal(palette="Set3")(12))(seq(0, 1, length.out = length(grp)))
  }

  E(g)$alpha <- scales::rescale(E(g)$weight)
  E(g)$color <- "#D9D9D9"
  for(c in seq_along(grp)){
    if(length(grp[[c]])>0){
      V(g)[grp[[c]]]$color   <-  groupColours[c]

      if(length(E(g)[from(V(g)[grp[[c]]])])>0){
        id<-E(g)[from(V(g)[grp[[c]]])]$color%in%"#D9D9D9"
        if(any(id)){
          E(g)[from(V(g)[grp[[c]]])[id]]$color <- add_alpha(groupColours[c],alpha = E(g)[from(V(g)[grp[[c]]])[id]]$alpha)
        }
      }
    }
  }
  return(g)
}


plot.loglog <- function(fd.OUT){
  require(ggplot2)
  require(scales)
  g <- ggplot2::ggplot(fd.OUT$PLAW, aes(x=size,y=bulk), na.rm=T) +
    scale_x_log10(breaks = log_breaks(n=abs(diff(range(round(log10(fd.OUT$PLAW$size)))+c(-1,1))),base=10),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = range(round(log10(fd.OUT$PLAW$size)))+c(-1,1)) +
    scale_y_log10(breaks = log_breaks(n=abs(diff(range(round(log10(fd.OUT$PLAW$bulk)))+c(-1,1))),base=10),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = range(round(log10(fd.OUT$PLAW$bulk)))+c(-1,1)) +
    geom_point() +
    geom_abline(intercept = fd.OUT[[2]]$fitlm1$coefficients[[1]], slope = fd.OUT[[2]]$fitlm1$coefficients[[2]], colour = "red", size = 2) +
    ggtitle(paste("Regression over ",length(fd.OUT[[2]]$fitlm1$fitted.values)," frequencies/bins",sep=""))+
    xlab("Frequency (log10)")+ylab("Power (log10)") +
    annotation_logticks() +
    annotate("text",x=10^-2,y=10^5,label=paste("Slope = ",round(fd.OUT[[2]]$alpha,digits=2),sep="")) +
    gg.theme("clean")
  return(g)
}

## Add an alpha value to a colour
add_alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}

# Variability Analyses --------------------------------------------------------------------------------------------------------------------------

#
# #' PSDslope
# #'
# #' @param y    A time series object, or a vector that can be converted to a time series object.
# #' @param fs    Sample frequency (defults to 1).
# #' @param nfft    Number of frequencies to estimate (defaults to next power of 2)
# #' @param fitRange    Vector of length 2 with range of frequencies to perform log-log fit.
# #' @param plot    Plot the log-log spectrum and slope.
# #'
# #' @return
# #' @export
# #'
# #' @examples
# #'
# PSDslope <- function(y  = ts(rnorm(n = 1024), frequency = 1),
#                      fs = frequency(y),
#                      nfft = 2^(nextpow2(length(y)/2)),
#                      fitRange = c(1,round(.1*nfft)),
#                      plot = FALSE){
#   require(oce)
#   require(signal)
#   if(!is.ts(y)){ts(y, frequency = fs)}
#
#   win <- signal::hamming(n=nfft)
#
#   perioGram <- oce::pwelch(x = y, window = win, fs = frequency(y), nfft = nfft, plot = FALSE)
#   spec <- data.frame(Frequency = perioGram$freq, Power = perioGram$spec)
#   spec[1,1:2] <- NA
#   fit <- lm(log10(spec$Power[fitRange[1]:fitRange[2]])~log10(spec$Power[fitRange[1]:fitRange[2]]))
#   return(list(spec = spec,
#               slope = fit)
#   )
# }

#' elascer
#'
#' @description The 'elastic scaler'will rescale numeric vectors (1D, or columns in a matrix or data.frame) to a user defined minimum and maximum, either based on the extrema in the data, or, a minimum and maximum defined by the user.
#'
#' @param x     Input vector or data frame.
#' @param mn     Minimum value of original, defaults to \code{min(x, na.rm = TRUE)}.
#' @param mx     Maximum value of original, defaults to \code{max(x, na.rm = TRUE)}.
#' @param hi     Minimum value to rescale to, defaults to \code{0}.
#' @param lo     Maximum value to rescale to, defaults to \code{1}.
#'
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
#' elascer(somenumbers)
#' elascer(somenumbers,mn=-100)
#
#' # Values < mn will return < lo (default=0)
#' # Values > mx will return > hi (default=1)
#' elascer(somenumbers,mn=-1,mx=99)
#'
#' elascer(somenumbers,lo=-1,hi=1)
#' elascer(somenumbers,mn=-10,mx=101,lo=-1,hi=4)
elascer <- function(x,mn=min(x,na.rm=T),mx=max(x,na.rm=T),lo=0,hi=1){
  x <- as.data.frame(x)
  u <- x
  for(i in 1:dim(x)[2]){
    mn=min(x[,i],na.rm=T)
    mx=max(x[,i],na.rm=T)
    if(mn>mx){warning("Minimum (mn) >= maximum (mx).")}
    if(lo>hi){warning("Lowest scale value (lo) >= highest scale value (hi).")}
    ifelse(mn==mx,{u[,i]<-rep(mx,length(x[,i]))},{
      u[,i]<-(((x[i]-mn)*(hi-lo))/(mx-mn))+lo
      id<-complete.cases(u[,i])
      u[!id,i]<-0
    })
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

# MULTIPLOT FUNCTION ------------------------------------------------------------------------------------------------------------------
#
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
multi.PLOT <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
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
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
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
try.CATCH <- function(expr){
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}
