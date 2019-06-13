# (C)RQA ------------------------

#' crqa_cl_main
#'
#' @inheritParams crqa_cl
#'
#' @keywords internal
#'
#' @export
#'
crqa_cl_main <- function(data,
                         emDim  = 1,
                         emLag  = 1,
                         emRad  = NA,
                         DLmin = 2,
                         VLmin = 2,
                         theiler = 0,
                         win     = min(length(y1),ifelse(is.null(y2),(length(y1)+1), length(y2)), na.rm = TRUE),
                         step    = step,
                         JRP     = FALSE,
                         distNorm       = c("EUCLIDEAN", "MAX", "MIN", "OP")[[1]],
                         returnMeasures = TRUE,
                         returnRPvector = FALSE,
                         returnLineDist = FALSE,
                         doPlot = c("noplot","rp","distmat")[[1]],
                         path_to_rp = getOption("casnet.path_to_rp"),
                         saveOut    = FALSE,
                         path_out   = NULL,
                         file_ID    = NULL,
                         silent     = TRUE,
                         targetValue  = .05,
                         useParallel = FALSE,
                         ...){


  if(useParallel){
    #data <- unnest(data)
    returnRPvector = FALSE
    returnLineDist = FALSE
  }
  y1     <- data[,1]
  y2     <- data[,2]

  fixedRR <- FALSE
  if(is.na(emRad)){fixedRR=TRUE}

  file_ID <- file_ID%00%0

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

  if(doPlot%in%"distmat"){-1 * emRad}

  tmpd  <- tempdir()
  tmpf1 <- tempfile(tmpdir = tmpd, fileext = ".dat")
  utils::write.table(as.data.frame(y1), tmpf1, col.names = FALSE, row.names = FALSE, sep = "\t")

#  fileSep <- ifelse(any(path_out%in%"/"),"/","\\")
file_ID=1
  plotOUT     <- file.path(normalizePath(path_out,mustWork = TRUE),paste0("RQAplot_",     file_ID, ".txt"))
  measOUT     <- normalizePath(file.path(path_out,paste0("RQAmeasures_", file_ID, ".txt")), mustWork = FALSE)
  histOUTdiag <- normalizePath(file.path(path_out,paste0("RQAhist_diag_",file_ID, ".txt")), mustWork = FALSE)
  histOUThori <- normalizePath(file.path(path_out,paste0("RQAhist_hori_",file_ID, ".txt")), mustWork = FALSE)

  if(any(is.null(y2))|any(is.na(y2%00%NA))){

   if(is.na(emRad)){
      if(!is.na(targetValue)){
         emRad <- crqa_radius(y1 = y1, emDim = emDim, emLag = emLag, targetValue = targetValue, radiusOnFail = "percentile", tol = .2, silent = silent)
         if(emRad$Converged){
           emRad <- emRad$Radius
         } else {
             emRad <-  emRad <- stats::sd(c(y1,y2),na.rm = T)
             }
       }
    }
    opts <- paste("-i", shQuote(tmpf1),
                  "-r", shQuote(plotOUT),
                  "-o", shQuote(measOUT),
                  "-p", shQuote(histOUTdiag),
                  "-q", shQuote(histOUThori),
                  "-m", emDim,
                  "-t", emLag,
                  "-e", emRad,
                  "-l", DLmin,
                  "-v", VLmin,
                  "-w", theiler,
                  "-n", shQuote(distNorm),
                  ifelse(silent,"-s",""))
  } else {
    if(is.na(emRad)){
      if(!is.na(targetValue)){
        emRad <- crqa_radius(y1 = y1, y2 = y2, emDim = emDim, emLag = emLag, targetValue = targetValue, tol = .2, radiusOnFail = "percentile", silent = silent)
        if(emRad$Converged){
          emRad <- emRad$Radius
        } else {
            emRad <- stats::sd(c(y1,y2),na.rm = T)
            }
      }
    }
    tmpf2 <- tempfile(tmpdir = tmpd, fileext = ".dat")
    utils::write.table(as.data.frame(y2), tmpf2, col.names = FALSE, row.names = FALSE, sep = "\t")
    opts <- paste("-i", shQuote(tmpf1),
                  "-j", shQuote(tmpf2),
                  "-r", shQuote(plotOUT),
                  "-o", shQuote(measOUT),
                  "-p", shQuote(histOUTdiag),
                  "-q", shQuote(histOUThori),
                  "-m", emDim,
                  "-t", emLag,
                  "-e", emRad,
                  "-l", DLmin,
                  "-v", VLmin,
                  "-w", theiler,
                  "-n", shQuote(distNorm),
                  ifelse(silent,"-s",""))
  }

  #closeAllConnections()

  ## RCMD
   #callr::rcmd_copycat(cmd = paste0(getOption("casnet.rp_command")), cmdargs = opts, wd = file.path(normalizePath(path_to_rp, mustWork = FALSE)), show = silent, echo = TRUE)

  system2(command = file.path(path_to_rp,getOption("casnet.rp_command")), args = opts)

  measures     <- try_CATCH(utils::read.delim(normalizePath(gsub("[']+","",measOUT)),header=TRUE))
  rpMAT        <- try_CATCH(utils::read.delim(normalizePath(gsub("[']+","",plotOUT)),header=TRUE))
  disthistDiag <- try_CATCH(utils::read.delim(normalizePath(gsub("[']+","",histOUTdiag)), header=FALSE, sep = " "))
  disthistHori <- try_CATCH(utils::read.delim(normalizePath(gsub("[']+","",histOUThori)), header=FALSE, sep = " "))

  if(all(is.null(measures$warning),is.data.frame(measures$value))){
    measures <- measures$value
  } else {
    measures <- rbind.data.frame(rep(NA,length(RQAmeasures)))
  }
  colnames(measures) <- gsub("(#RR)|(X.RR)","RR",colnames(measures))
  measures <- cbind.data.frame(measures, emDim=emDim, emLag=emLag, emRad=emRad, DLmin=DLmin, VLmin=VLmin, distNorm =distNorm)

  if(all(is.null(rpMAT$warning),is.data.frame(rpMAT$value))){
    rpMAT <- rpMAT$value
  } else {
    rpMAT <- data.frame(y1=NA,y2=NA,dist=NA)
  }

  if(any(is.null(disthistDiag$warning),is.data.frame(grepl("Error",paste(disthistDiag$value))))){
    disthistDiag <- disthistDiag$value
  } else {
    disthistDiag <- data.frame(line.length=NA,freq=NA)
  }

  if(any(is.null(disthistHori$warning),is.data.frame(grepl("Error",paste(disthistHori$value))))){
    disthistHori <- disthistHori$value
  } else {
    disthistHori <- data.frame(line.length=NA,freq=NA)
  }

  if(!silent){cat(paste0("[ID ",file_ID,"] Analysis completed... "))}
  if(!saveOut){
    file.remove(measOUT)
    file.remove(plotOUT)
    file.remove(histOUTdiag)
    file.remove(histOUThori)
   if(!silent){cat("temporary files removed ...\n")}
  } else {
    if(!silent){
      cat("files saved ...\n")
      cat(measOUT,"\n",plotOUT,"\n",histOUTdiag,"\n",histOUThori,"\n")
      }
  }

  if(useParallel){
    return(tibble::as_tibble(measures))
  } else{
  return(list(measures = measures,
                rpMAT    = rpMAT,
                diag_disthist = disthistDiag,
                hori_disthist = disthistHori))
  }
}


#' Fast (C)RQA (command line crp)
#'
#' This function will run the \href{http://tocsy.pik-potsdam.de/commandline-rp.php}{commandline Recurrence Plots} executable provided by Norbert Marwan.
#'
#' @param y1 Time series 1
#' @param y2 Time series 2 for Cross Recurrence Analysis (default = `NULL`)
#' @param emDim Embedding dimensions (default = `1`)
#' @param emLag Embedding lag (default = `1`)
#' @param emRad Radius on distance matrix (default = `1`)
#' @param DLmin Minimum length of diagonal structure to be considered a line (default = `2`)
#' @param VLmin Minimum length of vertical structure to be considered a line (default = `2`)
#' @param theiler Theiler window (default = `0`)
#' @param win Window to calculate the (C)RQA (default = minimum of length of `y1` or `y2`)
#' @param step Stepsize for sliding windows (default = size of `win`, so no sliding window)
#' @param JRP Wether to calculate a Joint Recurrence Plot (default = `FALSE`)
#' @param distNorm One of "EUCLIDEAN" (default), `"MAX", "MIN"`, or `"OP"` for an Order Pattern recurrence matrix
#' @param standardise Standardise data: `"none"` (default), `"mean.sd"`, or `"median.mad"`
#' @param returnMeasures Return the (C)RQA measures? (default = `TRUE`)
#' @param returnRPvector Return the recurrent points in a dataframe? (default = `FALSE`)
#' @param returnLineDist Return the distribution of diagonal and horizontal line length distances (default = `FALSE`)
#' @param doPlot Produce a plot of the recurrence matrix by calling [rp_plot()], values can be `"rp"` (the thresholded recurrence matrix),`"distmat"` (the unthresholded recurrence matrix) or `"noplot"` (default = `"noplot"`)
#' @param path_to_rp Path to the command line executable (default = path set during installation, use `getOption("casnet.path_to_rp")` to see)
#' @param saveOut Save the output to files? If `TRUE` and `path_out = NA`, the current working directory will be used (default = `FALSE`)
#' @param path_out Path to save output if `saveOut = TRUE` (default = `NULL`)
#' @param file_ID A file ID which will be a prefix to to the filename if `saveOut = TRUE` (default = `NULL`, an integer will be added tot the file name to ensure unique files)
#' @param silent Do not display any messages (default = `TRUE`)
#' @param surrogateTest Perform surrogate tests. If `TRUE`, will run surrogate tests using default settings for a two-sided test of \eqn{H_0: The data generating process is a rescaled linear Gaussian process} at \eqn{\alpha = .05} (arguments `ns = 39, fft = TRUE, amplitude = TRUE`)
#' @param targetValue A value passed to `crqa_radius(...,type="fixed", targetMeasure="RR")` if `is.na(emRad)==TRUE`. This is useful for windowed analysis, it will estimate a new radius for each window.
#' @param useParallel Speed up calculations by using the parallel processing options provided by `parallel` to assign a seperate process/core for each window in windowed (C)RQA analysis using [purrr::map2()] to assign data and [parallel::detectCores()] with  `logical = TRUE` to decide on the available cores (default = `FALSE`)
#' @param ... Additional parameters (currently not used)
#'
#' @details The `rp` executable is installed when the function is called for the first time and is renamed to `rp`, from a platform specific filename downloaded from <http://tocsy.pik-potsdam.de/commandline-rp.php> or extracted from an archive located in the directory: `...\\casnet\\commandline_rp\\`.
#' The file is copied to the directory: `...\\casnet\\exec\\`
#' The latter location is stored as an option and can be read by calling `getOption("casnet.path_to_rp")`.
#'
#' @section Troubleshooting:
#' Some notes on resolving errors with `rp`.The script will first try to download the correct executable, if that fails it will try to extract the file from a .zip archive in `...\\casnet\\commandline_rp\\crp_cl.zip`. If that fails, the copy will have failed. It should be relatively easy to get `crqa_cl()` working using custom settings:
#'
#' \itemize{
#' \item *Copy failed* - Every time the function `crqa_cl()` is called it will check whether a log file `rp_instal_log.txt` is present in the `...\\casnet\\exec\\` directory. If you delete the log file, and call the function, another copy of the executable will be attempted.
#' \item *Copy still fails and/or no permission to copy* - You can copy the approrpiate executable to any directory you have access to, be sure to rename it to `rp` (`rp.exe` on Windows OS). Then, either pass the path to `rp` as the argument `path_to_rp` in the `crqa_cl` function call, or, as a more permanent solution, set the `path_to_rp` option by calling `options(casnet.path_to_rp="YOUR_PATH")`. If you cannot acces the directory `...\\casnet\\commandline_rp\\`, download the appropriate executable from the \href{http://tocsy.pik-potsdam.de/commandline-rp.php}{commandline Recurrence Plots} page and copy to a directory you have access to. Then follow the instruction to set `path_to_rp`.
#' \item *Error in execution of `rp`* - This can have a variety of causes, the `rp` executable is called using [callr::rcmd()] and makes use of the [normalizePath()] function with argument `mustWork = FALSE`. Problems caused by specific OS, machine, or, locale problems (e.g. the `winslash` can be reported as an \href{https://github.com/FredHasselman/casnet/issues}{issue on Github}). One execution error that occurs when the OS is not recognised properly can be resolved by chekcing `getOption("casnet.rp_prefix")`. On Windows OS this should return an empty character vector, on Linux or macOS it should return `"./"`. You can manually set the correct prefix by calling `options(casnet.rp_prefix="CORRECT OS PREFIX")` and fill in the prefix that is correct for your OS
#' }
#'
#' @return A list object containing 1-3 elements, depending on arguments requesting output.
#'
#' \enumerate{
#' \item `rqa_measures` - A list of the (C)RQA measures returned if `returnMeasures = TRUE`:
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
#' \item `rqa_rpvector` - The radius thresholded distance matrix (recurrence matrix), which can be visualised as a recurrence plot by calling [rp_plot()]. If a sliding window analysis is conducted this will be a list of matrices and could potentially grow too large to handle. It is recommended you save the output to disk by setting `saveOut = TRUE`.
#' \item `rqa_diagdist` - The distribution of diagonal line lengths
#' }
#'
#' @note The platform specific `rp` command line executables were created by Norbert Marwan and obtained under a Creative Commons License from the website of the Potsdam Institute for Climate Impact Research at <http://tocsy.pik-potsdam.de/>.
#'
#' The full copyright statement on the website is as follows:
#'
#' (C) 2004-2017 SOME RIGHTS RESERVED
#'
#' University of Potsdam, Interdisciplinary Center for Dynamics of Complex Systems, Germany
#'
#' Potsdam Institute for Climate Impact Research, Transdisciplinary Concepts and Methods, Germany
#'
#' This work is licensed under a \href{https://creativecommons.org/licenses/by-nc-nd/2.0/de/}{Creative Commons Attribution-NonCommercial-NoDerivs 2.0 Germany License}.
#'
#' More information about recurrence analysis can be found on the \href{http://www.recurrence-plot.tk}{Recurrence Plot} website.
#'
#' @family Recurrence Quantification Analysis

#' @export
#'
crqa_cl <- function(y1,
                    y2      = NULL,
                    emDim   = 1,
                    emLag   = 1,
                    emRad   = NA,
                    DLmin   = 2,
                    VLmin   = 2,
                    theiler = 0,
                    win            = min(length(y1),ifelse(is.null(y2),(length(y1)+1), length(y2)), na.rm = TRUE),
                    step           = win,
                    JRP            = FALSE,
                    distNorm       = c("EUCLIDEAN", "MAX", "MIN", "OP")[[1]],
                    standardise    = c("none","mean.sd","median.mad")[1],
                    returnMeasures = TRUE,
                    returnRPvector = FALSE,
                    returnLineDist = FALSE,
                    doPlot         = c("noplot","rp","distmat")[[1]],
                    path_to_rp     = getOption("casnet.path_to_rp"),
                    saveOut        = FALSE,
                    path_out       = NULL,
                    file_ID        = NULL,
                    silent         = TRUE,
                    surrogateTest  = FALSE,
                    targetValue    = .05,
                    useParallel    = FALSE,
                    ...){

  if(!file.exists(normalizePath(file.path(getOption("casnet.path_to_rp"),"/rp_install_log.txt"), mustWork = FALSE))){
    set_command_line_rp()
  }

  cat("\n~~~o~~o~~casnet~~o~~o~~~\n")

  dotArgs <- list(...)
  if(is.na(dotArgs%00%NA)){dotArgs<-list(nothing=NA)}

  if(!is.null(y2)){
    y2 <- zoo::as.zoo(y2)
    N1 <- NROW(y1)
    N2 <- NROW(y2)
    if(N1>N2){
      y2[N2:(N1+(N1-N2))] <- 0
    }
    if(N1<N2){
      y1[N1:(N2+(N2-N1))] <- 0
    }
    df <- cbind.data.frame(y1=unclass(y1),y2=unclass(y2))
    df <- df[stats::complete.cases(df),]
  } else {
    if(any(is.na(y1))){
      y1 <- y1[!is.na(y1)]
      if(!silent){message("Removed NAs from timeseries y1 before (C)RQA")}
    }
    df <- cbind.data.frame(y1=unclass(y1),y2=NA)
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

  if(is.null(y2)){
    cat("\nPerforming auto-RQA\n")
    if(theiler<1){theiler=1}
  } else {
    cat("\nPerforming cross-RQA\n")
  }

  if(surrogateTest){
    surrogateSeries <- tseries::surrogate(y1,ns=39, fft = TRUE, amplitude = TRUE)
  }

  windowedAnalysis <- FALSE
  if(win==NROW(df)|step==(NROW(df))){
    win <- step <- NROW(df)
    wIndex <- seq(1,NROW(df))
    if(is.na(emRad)){cat(paste("Radius will be calulated to yield targetValue =",targetValue,"RR...\n"))}
  } else {
    windowedAnalysis = TRUE
    wIndex <- seq(1,NROW(df)-win, by = step)
    cat(paste("\nCalculating recurrence measures in a window of size",win,"taking steps of",step,"data points. \n"))
    if(is.na(emRad)){cat(paste("Radius will be re-calculated to yield targetValue =",targetValue,"RR for each window...\n\n"))}
  }


  if(windowedAnalysis){

    # Check window length vs. embedding
    if(win < (emLag*emDim)){

      stop(paste0("The size of win = ", win, " must be larger than the product of emLag = ",emLag, " and emDim = ", emDim))

    } else {

  # Adjust time series lengths
  if((dplyr::last(wIndex)+win-NROW(df))>0){

    yy1 <- ts_trimfill(x=seq(1,dplyr::last(wIndex)+win),y=df[,1])
    yy2 <- rep(NA,length(yy1))
    if(!all(is.na(df[,2]))){
    yy2 <- ts_trimfill(x=seq(1,dplyr::last(wIndex)+win),y=df[,2])
    }
    df <- cbind.data.frame(y1=yy1,y2=yy2)
    if(!silent){message("\nAdded 0s to vector to make window and stepsize combination fit to time series length\n\n ")}
  } else {
    if(!silent){message("\nWindow and stepsize combination do not fit to the full length of the time series!\n\n ")}
  }

  if(silent){cat("\nsst! [it's a silent-ish mode]\n\n")}

  wIndices <- plyr::llply(wIndex, function(w){seq(w,w+(win))})
  names(wIndices) <- paste0("window: ",seq_along(wIndices)," | start: ",wIndex," | stop: ",wIndex+win)

  }
  } else {

    wIndices <- list(wIndex)
    names(wIndices) <- paste0("window: 1 | start: ",wIndex[1]," | stop: ",wIndex[NROW(df)])

  } # If windowed

  #cl <- parallel::makeCluster(mc.cores)
  #parallel::clusterExport(cl = cl, c("crqa_cl_main","wIndices","df","y1","y2","emDim","emLag","emRad","DLmin","VLmin","theiler", "win","step","JRP","distNorm","returnMeasures","returnRPvector","returnLineDist","doPlot","path_to_rp","saveOut","path_out","file_ID","silent","..."))

  # Create the data, for all windows,need this for parallel, but also speeds up processing in general.
  dfList <- plyr::llply(wIndices, function(ind){cbind(y1 = df[ind,1],y2 = df[ind,2])})

  if(useParallel){

    if(names(dotArgs)%in%"logical"){logical <- logical} else {logical <- TRUE}

   cat("\n...using parallel processing...\n")

   # Only return measures
   returnRPvector <- FALSE
   returnLineDist <- FALSE

   # How many cores is optimal?
   cores_available <- parallel::detectCores(logical = logical)-1
   if(cores_available>1){
    if(length(wIndices)%%2==0){odd_windows <- FALSE} else {odd_windows <- TRUE}
    if(cores_available%%2==0){odd_cores <- FALSE} else {odd_cores <- TRUE}
    if(odd_windows&!odd_cores){cores_available<-cores_available-1}
    if(!odd_windows&odd_cores){cores_available<-cores_available-1}
   } else {
     cores_available <- 1
   }

    cl       <- parallel::makeCluster(cores_available)

    parallel::clusterEvalQ(cl, library(callr))
    parallel::clusterEvalQ(cl,library(utils))
    parallel::clusterEvalQ(cl,library(plyr))
    parallel::clusterEvalQ(cl,library(tidyverse))
    parallel::clusterEvalQ(cl,library(pROC))
    parallel::clusterEvalQ(cl,library(Matrix))
    parallel::clusterEvalQ(cl,library(casnet))

   # parallel::clusterExport(cl, varlist = c("data","emDim","emLag","emRad","DLmin","VLmin","theiler","win","step","JRP","distNorm","returnMeasures","returnRPvector","returnLineDist","doPlot","path_to_rp", "saveOut","path_out","file_ID","silent","targetValue", "useParallel"))

      # cluster_library(c("devtools","utils","plyr","dplyr","tidyr","Matrix","pROC")) %>%
      # cluster_assign_value("crqa_cl_main", crqa_cl_main) %>%
      # cluster_assign_value("crqa_radius", crqa_radius) %>%
      # cluster_assign_value("emDim", emDim)  %>%
      # cluster_assign_value("emLag", emLag)  %>%
      # cluster_assign_value("emRad", emRad)  %>%
      # cluster_assign_value("DLmin", DLmin)  %>%
      # cluster_assign_value("VLmin", VLmin)  %>%
      # cluster_assign_value("theiler", theiler)  %>%
      # cluster_assign_value("JRP", JRP)  %>%
      # cluster_assign_value("distNorm", distNorm)  %>%
      # cluster_assign_value("returnMeasures", returnMeasures)  %>%
      # cluster_assign_value("returnRPvector", returnRPvector)  %>%
      # cluster_assign_value("returnLineDist", returnLineDist)  %>%
      # cluster_assign_value("doPlot", doPlot)  %>%
      # cluster_assign_value("path_to_rp", path_to_rp)  %>%
      # cluster_assign_value("saveOut", saveOut)  %>%
      # cluster_assign_value("path_out", path_out)  %>%
      # cluster_assign_value("file_ID", file_ID)  %>%
      # cluster_assign_value("silent", silent)  %>%
      # cluster_assign_value("targetValue", targetValue) %>%
      # cluster_assign_value("useParallel", useParallel)

 #  parallel::clusterExport(cl = cl, c("crqa_cl_main","wIndices","df","y1","y2","emDim","emLag","emRad","DLmin","VLmin","theiler", "win","step","JRP","distNorm","returnMeasures","returnRPvector","returnLineDist","doPlot","path_to_rp","saveOut","path_out","file_ID","silent","..."))


  start <- proc.time()
  wList <- parallel::parLapply(cl,dfList,function(df){crqa_cl_main(data = df,
                                                                       emDim          = emDim,
                                                                       emLag          = emLag,
                                                                       emRad          = emRad,
                                                                       DLmin          = DLmin,
                                                                       VLmin          = VLmin,
                                                                       theiler        = theiler,
                                                                       JRP            = JRP,
                                                                       distNorm       = distNorm,
                                                                       returnMeasures = returnMeasures,
                                                                       returnRPvector = returnRPvector,
                                                                       returnLineDist = returnLineDist,
                                                                       doPlot    = doPlot,
                                                                       path_to_rp     = path_to_rp,
                                                                       saveOut        = saveOut,
                                                                       path_out       = path_out,
                                                                       file_ID        = file_ID,
                                                                       silent         = silent,
                                                                       targetValue    = targetValue,
                                                                       useParallel    = useParallel)})

  parallel::stopCluster(cl)
  time_elapsed_parallel <- proc.time() - start # End clock

  cat("\nCompleted in:\n")
  print(time_elapsed_parallel)

  } else {

    cat("\n...using sequential processing...\n")

    start <- proc.time()

    wList <- plyr::llply(dfList, function(data){
      crqa_cl_main(
        data           = data,
        emDim          = emDim,
        emLag          = emLag,
        emRad          = emRad,
        DLmin          = DLmin,
        VLmin          = VLmin,
        theiler        = theiler,
        JRP            = JRP,
        distNorm       = distNorm,
        returnMeasures = returnMeasures,
        returnRPvector = returnRPvector,
        returnLineDist = returnLineDist,
        doPlot    = doPlot,
        path_to_rp     = path_to_rp,
        saveOut        = saveOut,
        path_out       = path_out,
        file_ID        = file_ID,
        silent         = silent,
        targetValue    = targetValue)},
      .progress = plyr::progress_text(char = "o~o"))

    time_elapsed_parallel <- proc.time() - start # End clock

    cat(paste("\nCompleted in:\n"))
    print(time_elapsed_parallel)


    } # useParallel
  rm(dfList)

     rqa_measures <- rqa_rpvector <- rqa_diagdist <- rqa_horidist <- NA
  if(useParallel){
    if(is.list(wList)){
      rqa_measures <-  plyr::ldply(wList)
    } else {
      rqa_measures <-  tidyr::unnest(wList)
    }
  } else {
    rqa_measures <-plyr::ldply(wList, function(l) l$measures) # %>% dplyr::mutate(win = win, step = step, index = attr(wlist, "index"))
    rqa_rpvector <-plyr::ldply(wList, function(l) l$rpMAT) # %>% dplyr::mutate(win = win, step = step, index = attr(wlist, "index"))
    rqa_diagdist <-plyr::ldply(wList, function(l) l$diag_disthist) # %>% dplyr::mutate(win = win, step = step, index = attr(wlist, "index"))
    rqa_horidist <-plyr::ldply(wList, function(l) l$hori_disthist) # %>% dplyr::mutate(win = win, step = step, index = attr(wlist, "index"))
  }

  wPlot <- which(doPlot%in%c("noplot","rp","distmat"))

  if(wPlot>1){
    if(wPlot==2){
      if(windowedAnalysis){
      plotList <- plyr::llply(wIndices, function(ind) rp(y1 = df[ind,1], y2 = df[ind,2], emDim = emDim, emLag = emLag, emRad = emRad, doPlot = TRUE))
      multi_PLOT(plotList)
      } else {
        rp(y1 = df[,1], y2 = df[,2], emDim = emDim, emLag = emLag, emRad = emRad, doPlot = TRUE)
      }
    }
    if(wPlot==3){
      if(windowedAnalysis){
       plotList <- plyr::llply(wIndices, function(ind) rp(y1 = df[ind,1],y2 = df[ind,2], emDim = emDim, emLag = emLag, doPlot = TRUE))
       multi_PLOT(plotList)
      } else {
        rp(y1 = df[,1],y2 = df[,2], emDim = emDim, emLag = emLag, doPlot = TRUE)
      }
    }

  }

  justData <- FALSE
  if(win!=NROW(df)){
    if(!silent){message("Windowed (c)rqa: Returning only (c)rqa measures, use saveOut = TRUE, to get the distribution files")}
    justData <- TRUE
  }

  if(returnMeasures&!returnRPvector&!returnLineDist){
    justData<- TRUE
  }

  out <- list(rqa_measures = rqa_measures,
              rqa_rpvector = rqa_rpvector,
              rqa_diagdist = rqa_diagdist,
              rqa_horidist = rqa_horidist)
  if(saveOut){saveRDS(out,paste0(path_out,"CRQA_out",file_ID,".rds"))}

  cat("\n~~~o~~o~~casnet~~o~~o~~~\n")

  if(justData){
    return(rqa_measures)
  } else {
    return(out[c(returnMeasures,returnRPvector,returnLineDist,returnLineDist)])
  }
}


#' Find optimal (C)RQA parameters
#'
#' A wrapper for various algorithms used to find optimal value pair for the embedding delay and the number of embedding dimensions
#'
#' @param y A numeric vector or time series
#' @param emLag Optimal embedding lag (delay), e.g., provided by an optimising algorithm. If `NULL` the lags based on the mutual information in `lagMethods` will be reported. If a numeric value representing a valid lag is passed, this value will be used to estimate the number of dimensions (default = `NULL`)
#' @param lagMethods A character vector with one or more of the following strings: `"first.minimum","global.minimum","max.lag"`. If `emLag` represents a valid lag this value will  be reported as "user.lag" (default = `c("first.minimum","global.minimum","max.lag")`)
#' @param maxLag Maximum embedding lag to consider. Default value is: `floor(length(y)/(maxDim+1))`
#' @param estimateDimensions Decide on an optimal embedding dimension relative to the values in `maxDim` and `lagMethods`, according to a number of preferences passed as a character vector. The order in which the preferences appear in the vector affects the selection procedure, with index `1` being most important preference. The following options are available: \itemize{
#' \item{`preferNone` - No optimal number will be picked all other preferences will be ignored}
#' \item{`preferSmallestDim` - Pick smallest number of dimensions associated with a percentage NN below `nnThres`}
#' \item{`preferSmallestNN` - Pick the number of dimensions that is associated with the smallest percentage NN below `nnThres`}
#' \item{`preferSmallestLag` - If the value of `nnThres` does not lead to a unique preference for a pair of dimension and lag values, use the pair with the smallest lag}
#' \item{`preferSmallestInLargestHood` - The default option: If no unique pair can be found, prefer pairs with smallest values for lag, dimensions, percentage NN for the largest NN size}
#' }
#' @param maxDim Maximum number of embedding dimensions to consider (default = `10`)
#' @param nnSizes Points whose distance is `nnSize` times further apart than the estimated size of the attractor will be declared false neighbours. See the argument `atol` in [fractal::FNN()] (default = `c(2,5,10,15)`)
#' @param nnRadius If the ratio of the distance between two points in successive dimensions is larger than `nnRadius`, the points are declared false neighbours. See the argument `rtol` in [fractal::FNN()] (default = `5`)
#' @param nnThres Threshold value representing the percentage of Nearest Neighbours that would be acceptable when using N surrogate dimensions. The smallest number of surrogate dimensions that yield a value below the threshold will be considered optimal (default = `10`)
#' @param theiler Theiler window on distance matrix (default = `0`)
#' @param diagPlot Produce a diagnostic plot the results (default = `TRUE`)
#' @param silent Silent-ish mode
#' @param ... Other parameters passed to [nonlinearTseries::timeLag()]
#'
#' @return A list object containing the optimal values (as indicated by the user) and iteration history.
#'
#' @details A number of functions are called to determie optimal parameters for delay embedding a time series:
#'
#' \itemize{
#' \item{Embedding lag (\eqn{\tau}, `emLag`): The default is to call [casnet::est_emLag()], which is a wrapper around [nonlinearTseries::timeLag()] with `technique="ami"` to get lags based on the mutual information function.}
#' \item{Embedding dimension (`m`, `emDim`): The default is to call [casnet::est_emDim()], which is a wrapper around [fractal::FNN()]}
#' }
#'
#' @family Recurrence Quantification Analysis
#'
#' @export
#'
crqa_parameters <- function(y,
                            emLag     = NULL,
                            maxLag   = floor(length(y)/(maxDim+1)),
                            lagMethods = c("first.minimum","global.minimum","max.lag"),
                            estimateDimensions = "preferSmallestInLargestHood",
                            maxDim   = 10,
                            nnSizes  = c(2,5,10,15),
                            nnRadius = 5,
                            nnThres  = 10,
                            theiler  = 0,
                            diagPlot = TRUE,
                            silent   = TRUE,
                            ...){

  if(!is.null(dim(y))){stop("y must be a 1D numeric vector!")}

  if(length(nnRadius)!=1){stop("nnRadius must have 1 numeric value")}
  if(length(nnSizes)!=4){stop("nnSizes must have 4 numeric values")}
  y <- y[!is.na(y)]
  #y <- ts_standardise(y, adjustN = FALSE)

  emDims  <-  1:maxDim

  doLags <- c(1:maxLag)
  if(!is.null(emLag)){
    if(NROW(emLag==1)){
      lagMethods <- c(lagMethods, "user.lag") #  lag = emLag, ami = 0)
      doLags <- unique(sort(c(1:maxLag,emLag)))
    } else {
      stop("emLag must have 1 numeric value")
    }
  }

    if(nchar(estimateDimensions)>1){
      mi <- mif(data.frame(y),lags = doLags)
      #est_emLag(y,selection.methods = lagMethods, maxLag = maxLag)
      emLags <- cbind.data.frame(selection.methods = lagMethods, lag = NA)
      for(m in seq_along(emLags$selection.methods)){
        if(emLags$selection.methods[m]=="first.minimum"){
          emLags$lag[m] <- which(ts_symbolic(data.frame(mi))[,2]%in%"trough")[1]%00%NA
        }
        if(emLags$selection.methods[m]=="global.minimum"){
          emLags$lag[m] <- as.numeric(which.min(mi))
        }
        if(emLags$selection.methods[m]=="max.lag"){
          emLags$lag[m] <- maxLag
        }
        if(emLags$selection.methods[m]=="user.lag"){
          emLags$lag[m] <- emLag
        }
        emLags$ami[m] <- mi[emLags$lag[m]]
      }
    } else {
      emLags <- cbind.data.frame(selection.methods = "Not estimated", lag = NA, ami = NA)
    }


  # (fn.out <- tseriesChaos::false.nearest(lx, m=10, d=17, t=0, eps=sd(lx)/10, rt=20))
  # plot(fn.out[1,],type="b")

if(any(estimateDimensions%in%c("preferNone","preferSmallestDim", "preferSmallestNN", "preferSmallestLag", "preferSmallestInLargestHood"))){

  lagList <- list()
  cnt = 0

  for(N in seq_along(nnSizes)){
    for(R in seq_along(nnRadius)){
      for(L in seq_along(emLags$selection.method)){

     if(!is.na(emLags$lag[L])){
        cnt = cnt+1

        Nn.max <- Nn.mean <- Nn.sd <- Nn.min <- numeric(maxDim)
      for(D in seq_along(emDims)){
        RM <- rp(y,y, emDim = emDims[[D]], emLag = emLags$lag[L])
        RM <- bandReplace(RM,-theiler,theiler,0,silent = silent)
        Nn.min[D] <- min(RM, na.rm = TRUE)
        Nn.max[D] <- max(RM, na.rm = TRUE)
        Nn.sd[D]   <- stats::sd(RM, na.rm = TRUE)
        Nn.mean[D] <- mean(RM, na.rm = TRUE)
        rm(RM)
      }

        #surrDims <- nonlinearTseries::buildTakens(time.series =  as.numeric(y), embedding.dim =  emDims[[D]], time.lag = emLags$optimal.lag[L])

        # (fn.out <- false.nearest(rnorm(1024), m=6, d=1, t=1, rt=3))
        # plot(fn.out)

        fnnSeries <- fractal::FNN(x = as.numeric(y),
                                  dimension = maxDim,
                                  tlag = emLags$lag[L],
                                  rtol = nnRadius[R],
                                  atol = nnSizes[N],
                                  olag = 1
        )
        }

        #allN <- nonlinearTseries::findAllNeighbours(surrDims, radius = nnSizes[N]*sd(y))
        #Nn <- sum(plyr::laply(allN, length), na.rm = TRUE)
      #   if(D==1){Nn.max <- Nn}
        lagList[[cnt]] <- data.frame(Nn.pct = as.numeric(fnnSeries[1,]),
                                     Nsize = nnSizes[N],
                                     Nradius = nnRadius[R],
                                     emLag.method = emLags$selection.method[[L]],
                                     emLag = emLags$lag[L],
                                     emDim = emDims,
                                     Nn.mean = Nn.mean,
                                     Nn.sd  = Nn.sd,
                                     Nn.min = Nn.min,
                                     Nn.max = Nn.max)
       } # L
      } # R
    } # N


  df        <- plyr::ldply(lagList)
 # df$Nn.pct <- df$Nn/df$Nn.max

  opt <-plyr::ldply(unique(df$emLag), function(n){
    id <- which((df$Nn.pct<=nnThres)&(df$emLag==n)&(!(df$emLag.method%in%"maximum.lag")))
    if(length(id)>0){
      #idmin <- id[df$emDim[id]==min(df$emDim[id], na.rm = TRUE)]
      # if(length(idmin)>0){
      #   return(df[idmin,])
      #}
      df <- df[id,]
      return(df[!duplicated(df),])
    } else {
      return(df[which.min(df$Nn.pct[(df$emLag==n)&(!(df$emLag.method%in%"maximum.lag"))]),])
    }
  }
  )

 opt <- switch(estimateDimensions,
         preferNone = opt[!duplicated(opt),],
         preferSmallestDim = opt[min(opt$emDim, na.rm=TRUE),],
         preferSmallestNN = opt[min(opt$NN.pct, na.rm=TRUE),],
         preferSmallestLag = opt[min(opt$emLag, na.rm=TRUE),],
         preferSmallestInLargestHood = opt[(min(opt$emLag, na.rm=TRUE)&min(opt$emDim, na.rm=TRUE)&max(opt$Nsize, na.rm = TRUE)),]
         )

  #opt <- opt[opt$emDim==min(opt$emDim),][1,]
  opDim <- min(unique(opt$emDim), na.rm = TRUE)
  #opt <- opt[!duplicated(opt),]

  } else { # if estimateDim
    opDim <- NA
  }


  #opDim <- min(df$emDim[df$Nn.pct<nnThres], na.rm = TRUE)
  # opLag <- tau(y,
  #              selection.methods = ami.method,
  #              maxLag =maxLag)$opLag[1]
  opLag <- min(unique(opt$emLag), na.rm = TRUE)
  #opRad = NULL

  opt <- opt[all(opt$emDim==opDim,opt$emLag==opLag),]

  df$emLag <- factor(df$emLag)
  df$Nns   <- interaction(df$Nsize,df$Nradius)
  df$Nns.f  <-  factor(df$Nns, levels=levels(df$Nns), labels = paste0("size=",nnSizes," | radius=",nnRadius))

  if(diagPlot){

    dfs <- data.frame(startAt= c(.5, graphics::hist(emDims,plot=FALSE)$mids),
                      stopAt = c(graphics::hist(emDims,plot=FALSE)$mids,max(emDims)+.5),
                      f=factor(seq_along(c(.5, graphics::hist(emDims,plot=FALSE)$mids))%%2))

    #tmi <-  nonlinearTseries::mutualInformation(y, lag.max = maxLag, n.partitions = , do.plot = FALSE)

    tmi <- mif(data.frame(y),lags = 1:maxLag)

    dfMI <- data.frame(emDelay = as.numeric(names(tmi)),
                       ami     = as.numeric(tmi))


    #  RColorBrewer::brewer.pal(n = length(unique(df$emLag) name = "Set2")

    # use: alpha()

   Ncol <- length(emLags$selection.method[!is.na(emLags$lag)])
   myPal <- RColorBrewer::brewer.pal(Ncol,"Set2")
  myPalLag <- myPal
  names(myPalLag) <- emLags$selection.method[!is.na(emLags$lag)]
  myPalNn <- myPal
  names(myPalNn) <- emLags$lag[!is.na(emLags$lag)]

    gNdims <- ggplot2::ggplot(df, ggplot2::aes_(y = ~Nn.pct, x = ~emDim, colour = ~emLag)) +
      ggplot2::geom_rect(ggplot2::aes_(xmin = ~startAt, xmax = ~stopAt, fill = ~f), ymin = -Inf, ymax = Inf, data = dfs, inherit.aes = FALSE) +
      ggplot2::scale_fill_manual(values = scales::alpha(c("grey", "white"),.2), guide=FALSE) +
      ggplot2::geom_hline(yintercept = nnThres, linetype = 2, colour = "grey60") +
      ggplot2::geom_hline(yintercept = c(0,100),   colour = "grey60") +
      ggplot2::geom_hline(yintercept = 50, colour = "grey90") +
      ggplot2::geom_line(position  = ggplot2::position_dodge(.4)) +
      ggplot2::geom_point(position = ggplot2::position_dodge(.4)) +
      ggplot2::annotate("text",x=maxDim/3,y=nnThres, label="threshold", size = .8) +
      ggplot2::xlab("Embedding Dimension") +
      ggplot2::ylab("Nearest neigbours (% of max.)") +
      ggplot2::facet_wrap(~Nns.f, ncol=2) +
      ggplot2::scale_x_continuous(breaks=emDims) +
      ggplot2::scale_y_continuous(breaks = c(nnThres,50,100)) +
      ggplot2::scale_color_manual("Lag", values = myPalNn ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(strip.background = element_rect(colour = "grey90", fill = "grey90"),
            strip.text.x = element_text(colour = "black", face = "bold"),
            panel.spacing = ggplot2::unit(1, "lines"),
            legend.background = element_rect(colour = "grey90",fill = "grey90"),
            legend.title = element_text(face = "bold"),
            legend.key = element_rect(colour = "grey90", fill = "grey90"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank()
      )

    gDelay <- ggplot2::ggplot(dfMI, ggplot2::aes_(y = ~ami, x = ~emDelay)) +
      ggplot2::geom_line() +
      ggplot2::geom_vline(data = emLags,  ggplot2::aes_(colour=factor(emLags$selection.method),
                                    xintercept = ~lag), alpha = .3) +
      ggplot2::geom_point(data = emLags,  ggplot2::aes_(x = ~lag, y = ~ami, colour = factor(emLags$selection.method)), size = 2) +
      ggplot2::xlab("Embedding Lag") +
      ggplot2::ylab("Average Mututal Information") +
      ggplot2::scale_color_manual("Lag", values = myPalLag) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = c(.95, .95),
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

    g <- gridExtra::grid.arrange(gDelay, gNdims, ncol=1, nrow=2)
    grid::grid.newpage()
    grid::grid.draw(g)

  } else {
    g <- NA
  }

  return(list(optimLag  = opLag,
              optimDim  = opDim,
              optimRow  = opt,
              optimData = df,
              diagPlot = invisible(g)))
}


#' Find fixed or optimal radius
#'
#' @param RM Unthresholded Recurrence Matrix
#' @param y1  A numeric vector or time series
#' @param y2  A numeric vector or time series
#' @param emLag Delay to use for embedding
#' @param emDim Number of embedding dimensions
#' @param type Either `"fixed"` (default) or `"optimal"`, `"fixed"` will search for a radius that is close to the value for the `targetMeasure` in `targetValue`, `"optimal"` will optimise the radius for the `targetMeasure`, `targetValue` is ignored.
#' @param startRadius If `type = "fixed"` this is the starting value for the radius (default = percentile of unique distances in RM given by `targetValue`). If `type = "optimal"` this will be a range of radius values (in normalised SD units) that will be considered (default = `seq(0,2,by=.01)`)
#' @param eachRadius If `type = "optimal"` this is the number of signal and noise series that will be generated for each level in `startRadius` (default = `1`)
#' @param targetMeasure If `type = "optimal"`, it must be a character vector indicating which recurrence measure to optimise the radius for, options are "RR" (default), "DET", "LAM", "T1", and "all". The option `targetMeasure = "all"` will report all the optimal values obtained from one realisation of `startRadius * eachRadius` signal and noise series.
#' @param targetValue When argument `type` is set to "fixed", the value represents the target value for the measure in `targetMeasure` (default = `RR = .05`).
#' @param tol Tolerance for achieving `targetValue` for `targetMeasure` (default = `0.1`)
#' @param maxIter If `type = "fixed"`: Maximum number of iterations to reach targetValue.
#' @param theiler Size of theiler window (default `0`)
#' @param histIter Return iteration history? (default = `FALSE`)
#' @param noiseLevel Noise level to construct the `signal + noiseLevel *` \eqn{N(\mu=0,\sigma=1)} (default = `0.75`)
#' @param noiseType Type
#' @param plotROC Generates an ROC plot if `type = "optimal"`
#' @param standardise Standardise
#' @param radiusOnFail Radius to return when search fails `"tiny" = 0 + ,Machine.double.eps`, this will likely cause a matrix full of zeros. `"huge" = 1 + max. distance in RM`, which will give a matrix full of ones, `"percentile" = quantile(RM, prob = targetValue) of distances greater than 0`.
#' @param silent Silent-ish
#'
#' @family Recurrence Quantification Analysis
#'
#' @return A dataframe listing settings ussed to search for the radius, the radius found given the settings and the recurrence rate produced by the radius (either 1 row or the entire iteration history)
#' @export
#'
crqa_radius <- function(RM = NULL,
                        y1 = NULL,
                        y2 = NULL,
                        emLag = 1,
                        emDim = 1,
                        type           = c("fixed","optimal")[1],
                        startRadius    = NULL,
                        eachRadius     = 1,
                        targetMeasure  = c("RR","DET","LAM","T1","all")[1],
                        targetValue    = 0.05,
                        tol            = 0.1,
                        maxIter        = 100,
                        theiler        = -1,
                        histIter       = FALSE,
                        noiseLevel     = 0.75,
                        noiseType      = c("normal","uniform")[1],
                        plotROC        = FALSE,
                        standardise  = c("mean.sd","median.mad","none")[3],
                        radiusOnFail   = c("tiny","huge","percentile")[1],
                        silent         = FALSE){

  optimOK <- FALSE
  if(is.null(RM)&!is.null(y1)){
    optimOK <- TRUE
    RM <- rp(y1=y1, y2=y2,emDim=emDim, emLag=emLag)
  }

  # check auto-recurrence
  RM   <- rp_checkfix(RM, checkAUTO = TRUE)
  AUTO <- attr(RM,"AUTO")

  if(AUTO){
    if(!silent){cat(paste0("\nAuto-recurrence: Setting diagonal to (1 + max. distance) for analyses\n"))}
    if(theiler < 0) theiler <- 0
  }

  if(is.null(startRadius)){
    if(type=="fixed"){
      startRadius <- as.numeric(stats::quantile(unique(as.vector(RM[lower.tri(RM)])),probs = ifelse(targetValue>=1,.05,targetValue)))
    } else {
      startRadius <- seq(0,1.5,by=0.001)
    }
  }

  recmatsize <- rp_size(RM,AUTO,theiler)

  if(theiler>=0){
    RM <- bandReplace(RM,-theiler,theiler,1+max(RM),silent = silent)
  }

  if(type%in%"fixed"){

    if(tol%][%c(0,1)){stop("Argument tol must be between 0 and 1.")}

    tryRadius <- startRadius
    Measure   <- 0
    iter      <- 0
    Converged <- FALSE
    seqIter   <- 1:maxIter

    iterList <- data.frame(iter        = seqIter,
                           Measure     = Measure,
                           Radius      = tryRadius,
                           targetValue = targetValue,
                           tollo       = targetValue*(1-tol),
                           tolhi       = targetValue*(tol+1),
                           startRadius = startRadius,
                           rp.size = recmatsize,
                           AUTO        = AUTO,
                           Converged   = Converged, check.names = FALSE)

    exitIter <- FALSE
    if(!silent){cat(paste("\nSearching for a radius that will yield",targetValue, "for", targetMeasure,"\n"))}

   # p <- dplyr::progress_estimated(maxIter)
    while(!exitIter){


      iter <- iter+1
      #p$tick()$print()

      RMs <- di2bi(RM,emRad = tryRadius, convMat = TRUE)
      RT  <- Matrix::nnzero(RMs)
      rp.size <- length(RMs)
      Measure <- RT/rp.size

      # crpOut <- crqa_rp(RM = RMs, emRad = tryRadius, AUTO=AUTO)
      #Measure  <-  crpOut[[targetMeasure]]
      #crpOut <- data.frame(RR = RR, RT = RT, size = length(RMs))
      #Measure <- RR

      iterList[iter,] <-    cbind.data.frame(iter        = iter,
                                             Measure     = Measure,
                                             Radius      = tryRadius,
                                             targetValue = targetValue,
                                             tollo       = targetValue*(1-tol),
                                             tolhi       = targetValue*(tol+1),
                                             startRadius = startRadius,
                                             rp.size = recmatsize,
                                             AUTO        = AUTO,
                                             Converged   = Converged)

      if(any(Measure%[]%c(targetValue*(1-tol),targetValue*(tol+1)),(iter>=maxIter))){
        if(Measure%[]%c(targetValue*(1-tol),targetValue*(tol+1))){
          Converged <- TRUE
          if(!silent){message("\nConverged! Found an appropriate radius...")}
        }
        iterList$Converged[iter] <- Converged
        exitIter <- TRUE
      }

      if(round(Measure,digits = 2)>round(targetValue,digits = 2)){
        tryRadius <- tryRadius*(min(0.8,tol*2))
      } else {
        tryRadius <- tryRadius*(min(1.8,1+(tol*2)))
      }
    } # While ....

    if(iter>=maxIter){
      warning("Max. iterations reached!")
    }
    if(Measure %][% c(targetValue*(1-tol),targetValue*(tol+1))){
      iterList$Radius[iter] <- dplyr::case_when(
        radiusOnFail%in%"tiny" ~ 0 + .Machine$double.eps,
        radiusOnFail%in%"huge" ~ 1 + max(RM),
        radiusOnFail%in%"percentile" ~ as.numeric(stats::quantile(unique(as.vector(Matrix::tril(RM,-1))),probs = ifelse(targetValue>=1,.05,targetValue)))
      )
      warning(paste0("\nTarget not found, try increasing tolerance, max. iterations, or, change value of startRadius.\nreturning radius: ",iterList$Radius[iter]))
    }

    ifelse(histIter,id<-c(1:iter),id<-iter)
    return(iterList[id,])

  } # if "fixed"

  if(type%in%"optimal"){

    if(optimOK){

      if(!silent){cat(paste0("\nNormalisation set to: ",standardise,"!!\n"))}

      startRadius <- rep(startRadius, each = eachRadius)

      dfREC  <-  plyr::ldply(startRadius, function(r){roc_noise(y = y1,
                                                               emDim = emDim,
                                                               emLag = emLag,
                                                               emRad = r,
                                                               noiseLevel = noiseLevel,
                                                               standardise = standardise,
                                                               noiseType = noiseType)},
                             .progress = plyr::progress_text(char = "o~o"))


      dfREC      <-  dplyr::arrange(dfREC,dfREC$response,dfREC$emRad)
      caseID    <- dfREC$response%in%"signal+noise"
      controlID <- dfREC$response%in%"noise"

       rocREC <- list(RR  = pROC::roc(cases=dfREC$RR[caseID], controls=dfREC$RR[controlID]),
                      DET = pROC::roc(cases=dfREC$DET[caseID], controls=dfREC$DET[controlID]),
                      LAM = pROC::roc(cases=dfREC$LAM[caseID], controls=dfREC$LAM[controlID]),
                      T1  = pROC::roc(cases=dfREC$T1[caseID], controls=dfREC$T1[controlID]))


       optimal.radius <- list(RR  = dfREC$emRad[(dfREC$RR>=min(pROC::coords(rocREC$RR,
                                                                                "b",best.method="c",
                                                                                ret="t")))][1],
                              DET = dfREC$emRad[(dfREC$DET  >=min(pROC::coords(rocREC$DET,
                                                                                "b",best.method="c",
                                                                                ret="t")))][1],
                              LAM = dfREC$emRad[(dfREC$LAM  >=min(pROC::coords(rocREC$LAM,
                                                                                "b",best.method="c",
                                                                                ret="t")))][1],
                              T1  = dfREC$emRad[(dfREC$T1   >=min(pROC::coords(rocREC$T1,
                                                                                "b",best.method="c",
                                                                                ret="t")))][1]
                              )

       optimal.value <- list(RR  = crqa_cl(y1, emRad = optimal.radius$RR, emDim=emDim, emLag=emLag)[["RR"]],
                             DET = crqa_cl(y1, emRad = optimal.radius$DET, emDim=emDim, emLag=emLag)[["DET"]],
                             LAM = crqa_cl(y1, emRad = optimal.radius$LAM, emDim=emDim, emLag=emLag)[["LAM"]],
                             T1  = crqa_cl(y1, emRad = optimal.radius$T1, emDim=emDim, emLag=emLag)[["T1"]]
                             )

    if(targetMeasure=="all"){

       if(plotROC){

         op <- graphics::par(pty="s", mfrow=c(2,2))

         pROC::plot.roc(rocREC$RR,
                        print.thres="best",
                        print.thres.best.method="closest.topleft",
                        print.thres.cex=.6,
                        print.auc=TRUE,
                        print.auc.x=.9,
                        print.auc.y=.1,
                        auc.polygon=TRUE,
                        main = paste("RR =",round(optimal.value$RR,3),"- Radius =",optimal.radius$RR))

         pROC::plot.roc(rocREC$DET,
                        print.thres="best",
                        print.thres.best.method="closest.topleft",
                        print.thres.cex=.6,
                        print.auc=TRUE,
                        print.auc.x=.9,
                        print.auc.y=.1,
                        auc.polygon=TRUE,
                        main = paste("DET =",round(optimal.value$DET,3),"- Radius =",optimal.radius$DET))

         pROC::plot.roc(rocREC$LAM,
                        print.thres="best",
                        print.thres.best.method="closest.topleft",
                        print.thres.cex=.6,
                        print.auc=TRUE,
                        print.auc.x=.9,
                        print.auc.y=.1,
                        auc.polygon=TRUE,
                        main = paste("LAM =",round(optimal.value$LAM,3),"- Radius =",optimal.radius$LAM))

         pROC::plot.roc(rocREC$T1,
                        print.thres="best",
                        print.thres.best.method="closest.topleft",
                        print.thres.cex=.6,
                        print.auc=TRUE,
                        print.auc.x=.9,
                        print.auc.y=.1,
                        auc.polygon=TRUE,
                        main = paste("T1 =",round(optimal.value$T1,3),"- Radius =",optimal.radius$T1))

         #    plot(precision ~ recall, t(pROC::coords(rocREC, "all", ret = c("recall", "precision"))), type="l")
         #  graphics::text(x=-5,y=0,labels=paste("Optimal radius for", targetMeasure,"=",optimal.value,"is:",optimal.radius))
         graphics::par(op)
       }

       return(data.frame(measure = c("RR","DET","LAM","T1"),
                         optimal = plyr::ldply(optimal.value)[-1],
                         radius  = plyr::ldply(optimal.radius)[-1]
                         )
              )

     } else {

     #  rocREC <- dplyr::case_when(
     #    targetMeasure == "RR"  ~ list(rocREC=pROC::roc(controls=dfREC$RR[dfREC$response%in%"signal+noise"],
     #                                                   cases=dfREC$RR[dfREC$response=="noise"])),
     #    targetMeasure == "DET" ~ list(rocREC=pROC::roc(controls=dfREC$DET[dfREC$response%in%"signal+noise"],
     #                                                   cases=dfREC$DET[dfREC$response=="noise"])),
     #    targetMeasure == "LAM" ~ list(rocREC=pROC::roc(controls=dfREC$LAM[dfREC$response%in%"signal+noise"],
     #                                                   cases=dfREC$LAM[dfREC$response=="noise"])),
     #    targetMeasure == "T1"  ~ list(rocREC=pROC::roc(controls=dfREC$T1[dfREC$response%in%"signal+noise"],
     #                                                   cases=dfREC$T1[dfREC$response=="noise"]))
     #  )


      optimal.radius <- optimal.radius[[targetMeasure]]
      optimal.value <- optimal.value[[targetMeasure]]
      rocREC        <- rocREC[[targetMeasure]]

      if(plotROC){
       op<-graphics::par(pty="s")
        pROC::plot.roc(rocREC,
                       print.thres="best",
                       print.thres.best.method="closest.topleft",
                       print.thres.cex=.6,
                       print.auc=TRUE,
                       print.auc.x=.9,
                       print.auc.y=.1,
                       auc.polygon=TRUE,
                       main = paste(targetMeasure,"=",round(optimal.value,3),"- Radius =",optimal.radius))

    #    plot(precision ~ recall, t(pROC::coords(rocREC, "all", ret = c("recall", "precision"))), type="l")
      #  graphics::text(x=-5,y=0,labels=paste("Optimal radius for", targetMeasure,"=",optimal.value,"is:",optimal.radius))
        graphics::par(op)
      }

     }

      return(cbind.data.frame(measure=targetMeasure,
                        optimal.value=optimal.value[[1]]%00%NA,
                        optimal.radius=optimal.radius%00%NA))

    } else {

      stop("Need time series vector(s) to perform optimal Radius search.")

    }
  } # if "optimal"
}


#' roc_noise
#'
#' @param y y
#' @param emRad radius
#' @param emDim embedding Dims
#' @param emLag embedding Lag
#' @param noiseLevel noise Level
#' @param standardise Standardise y? Choose from "mean.sd","median.mad","none".
#' @param noiseType Use a Normal distribution of uniform distribution for noiselevels
#'
#' @return data frame for ROC
#' @export
#'
#' @keywords internal
roc_noise <- function(y, emRad, emDim=1, emLag=1, noiseLevel=.75, standardise = c("mean.sd","median.mad","none")[3], noiseType = c("normal","uniform")[1]){
  y <- dplyr::case_when(
    standardise == "mean.sd"   ~ ts_standardise(y, type="mean.sd"),
    standardise == "median.sd" ~ ts_standardise(y, type="median.mad"),
    standardise == "none"      ~ y
  )

  yn <- dplyr::case_when(
    noiseType == "normal"  ~ stats::rnorm(NROW(y), mean=round(mean(y, na.rm = TRUE),3), sd=stats::sd(y,na.rm = TRUE)),
    noiseType == "uniform" ~ sign(stats::rnorm(1))*stats::runif(NROW(y), min=floor(min(y, na.rm = TRUE)), max = ceiling(max(y,na.rm = TRUE)))
    )

  noise_out   <- crqa_cl(yn, emRad = emRad, emDim=emDim, emLag=emLag)
  measure_out <- crqa_cl((y  + noiseLevel * yn), emRad = emRad, emDim=emDim, emLag=emLag)

  return(cbind.data.frame(radius   = emRad,
                          response = c("signal+noise","noise"),
                          rbind(measure_out,noise_out)))
}



#' Get (C)RQA measures from a Recurrence Matrix
#'
#' Use `crqa_rp`
#'
#' @inheritParams crqa_rp
#'
#' @family Recurrence Quantification Analysis
#'
#' @export
#'
#' @keywords internal
#'
crqa_rp_measures <- function(RM,
                              emRad = NULL,
                              DLmin = 2,
                              VLmin = 2,
                              HLmin = 2,
                              DLmax = length(Matrix::diag(RM)),
                              VLmax = length(Matrix::diag(RM)),
                              HLmax = length(Matrix::diag(RM)),
                              AUTO      = NULL,
                              theiler = NULL,
                              chromatic = FALSE,
                              matrices  = FALSE,
                              doHalf    = FALSE,
                              Nboot     = NULL,
                              CL        = .95,
                              doParallel = FALSE){


  if(is.null(Nboot)){Nboot = 1}

  NRows <- NROW(RM)
  NCols <- NCOL(RM)
  mc.cores <- parallel::detectCores()-1
  if(Nboot<mc.cores) mc.cores <- Nboot


  if(doParallel){
  tstart <- Sys.time()
  cl <- parallel::makeCluster(mc.cores)
  out    <- parallel::mclapply(1, function(i){
    #crqa_rp_prep(matrix(RM[ceiling(NCols*NRows*stats::runif(NCols*NRows))], ncol=NCols, nrow = NRows),
    crqa_rp_prep(RP = RM,
             emRad= emRad,
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
  parallel::stopCluster(cl)
  tend  <- Sys.time()
  } else {
    out <- crqa_rp_prep(RP = RM,
               emRad= emRad,
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
  }


  dfori <- tidyr::gather(as.data.frame(out), key = "measure", value = "value")

  if(Nboot>1){


    cat(paste0("Bootstrapping Recurrence Matrix... ",Nboot," iterations...\n"))
    cat(paste0("Estimated duration: ", round((difftime(tend,tstart, units = "mins")*Nboot)/max((round(mc.cores/2)-1),1), digits=1)," min.\n"))

    tstart <-Sys.time()
    cl <- parallel::makeCluster(mc.cores)
    bootOut <-  parallel::mclapply(1:Nboot, function(i){
      replicate <- as.data.frame(crqa_rp_prep(matrix(RM[ceiling(NCols*NRows*stats::runif(NCols*NRows))],
                                                 ncol=NCols, nrow = NRows),
                                          emRad= emRad,
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
    parallel::stopCluster(cl)
    tend <- Sys.time()
    cat(paste0("Actual duration: ", round(difftime(tend,tstart, units = "mins"), digits=1)," min.\n"))

    dfrepl <- plyr::ldply(bootOut)
    dfrepl <- tidyr::gather(dfrepl, key = "measure", value = "value", -replicate)

    if(length(CL)==1){
      ci.lo <- (1-CL)/2
      ci.hi <- CL + ci.lo
    } else {
      ci.lo <- CL[1]
      ci.hi <- CL[2]
    }

    rqout <-  dfrepl %>%
      dplyr::group_by(dfrepl$measure) %>%
      dplyr::summarise(
        val          = NA,
        ci.lower     = stats::quantile(dfrepl$value%00%NA, ci.lo, na.rm = TRUE),
        ci.upper     = stats::quantile(dfrepl$value%00%NA, ci.hi, na.rm = TRUE),
        BOOTmean     = mean(dfrepl$value%00%NA, na.rm = TRUE),
        BOOTsd       = stats::sd(dfrepl$value%00%NA, na.rm = TRUE),
        BOOTse       = stats::sd(dfrepl$value%00%NA, na.rm = TRUE)/sqrt(Nboot),
        BOOTvar      = stats::var(dfrepl$value%00%NA, na.rm = TRUE),
        BOOTmedian   = stats::median(dfrepl$value%00%NA, na.rm = TRUE),
        BOOTmad      = stats::mad(dfrepl$value%00%NA, na.rm = TRUE),
        BOOTn        = Nboot
      )

    for(m in 1:nrow(rqout)){
      rqout$val[rqout$measure%in%dfori$measure[m]] <- dfori$value[m]

    }
  } else {
    rqout <- dfori %>% tidyr::spread(dfori$measure,dfori$value)
  }
  return(rqout)
}


#' Get bootsrapped (C)RQA measures based on a recurrence matrix
#'
#' A zoo of measures based on singular recurrent points, diagonal, vertical and horizontal line structures will be caluclated.
#'
#' @param RM A distance matrix, or a matrix of zeroes and ones (you must set `emRad = NA`)
#' @param emRad Threshold for distance value that counts as a recurrence
#' @param DLmin Minimal diagonal line length (default = `2`)
#' @param VLmin Minimal vertical line length (default = `2`)
#' @param HLmin Minimal horizontal line length (default = `2`)
#' @param DLmax Maximal diagonal line length (default = length of diagonal -1)
#' @param VLmax Maximal vertical line length (default = length of diagonal -1)
#' @param HLmax Maximal horizontal line length (default = length of diagonal -1)
#' @param AUTO Auto-recurrence? (default = `FALSE`)
#' @param theiler = Use a theiler window around the line of identity / synchronisation to remove high auto-correlation at short time-lags (default = `0`)
#' @param chromatic Force chromatic RQA? (default = `FALSE`)
#' @param matrices Return matrices? (default = `FALSE`)
#' @param doHalf Analyse half of the matrix? (default = `FALSE`)
#' @param Nboot How many bootstrap replications? (default = `NULL`)
#' @param CL Confidence limit for bootstrap results (default = `.95`)
#' @param targetValue A value passed to `crqa_radius(...,type="fixed", targetMeasure="RR", tol = .2)` if `is.na(emRad)==TRUE`, it will estimate a radius (default = `.05`).
#' @param doParallel Speed up calculations by using the parallel processing options provided by `parallel` to assign a seperate process/core for each window in windowed (C)RQA analysis using [purrr::map2()] to assign data and [parallel::detectCores()] with  `logical = TRUE` to decide on the available cores (default = `FALSE`)
#' @param silent Do not display any messages (default = `TRUE`)

#' @return A list object containing (C)RQA measures (and matrices if requested)
#'
#' @export
#'
#' @family Recurrence Quantification Analysis
#'
#'
crqa_rp <- function(RM,
                    emRad = NA,
                    DLmin = 2,
                    VLmin = 2,
                    HLmin = 2,
                    DLmax = length(Matrix::diag(RM))-1,
                    VLmax = length(Matrix::diag(RM))-1,
                    HLmax = length(Matrix::diag(RM))-1,
                    AUTO      = NULL,
                    theiler   = NULL,
                    chromatic = FALSE,
                    matrices  = FALSE,
                    doHalf    = FALSE,
                    Nboot     = NULL,
                    CL        = .95,
                    targetValue = .05,
                    doParallel = FALSE,
                    silent     = TRUE){


  # Input should be a distance matrix, or a matrix of zeroes and ones with emRad = NULL, output is a list
  # Fred Hasselman - August 2013

  #require(parallel)

  # check auto-recurrence and make sure Matrix has sparse triplet representation
  RM <- rp_checkfix(RM, checkAUTO = TRUE, fixAUTO = TRUE, checkTSPARSE = TRUE, fixTSPARSE = TRUE)

  if(is.null(AUTO)){
    AUTO <- attr(RM,"AUTO")
  }

  if(is.na(emRad)){
    if(!is.null(attributes(RM)$emRad)){
      emRad <- attributes(RM)$emRad
    } else {
      # Check for attributes
      if(is.na(targetValue)){
        targetValue <- .05
      }
      if(!is.null(attributes(RM)$emDim)){
        emDim <- attributes(RM)$emDim
      } else {
        emdim <- 1
      }
      if(!is.null(attributes(RM)$emLag)){
        emLag <- attributes(RM)$emLag
      } else {
        emLag <- 1
      }

      emRad <- crqa_radius(RM = RM,
                           emDim = emDim, emLag = emLag, targetValue = targetValue, tol = .2, radiusOnFail = "percentile", silent = silent)
      if(emRad$Converged){
        emRad <- emRad$Radius
      } else {
        emRad <- stats::sd(RM,na.rm = TRUE)
      }
    }
  }

  #uval <- unique(as.vector(RM))
  if(!all(as.vector(RM)==0|as.vector(RM)==1)){
    if(!is.null(emRad)){
      RM <- di2bi(RM,emRad)
    } else{
      if(!chromatic){
        stop("Expecting a binary (0,1) matrix.\nUse 'crqa_radius()', or set 'chromatic = TRUE'")
      } else {
        stop("Chromatic RQA not implemented yet.")
      }
    }
  # } else {
  #
   }
  #rm(RM)

  out <- crqa_rp_calc(RM,
                      emRad = emRad,
                      DLmin = DLmin,
                      VLmin = VLmin,
                      HLmin = HLmin,
                      DLmax = DLmax,
                      VLmax = VLmax,
                      HLmax = HLmax,
                      theiler = theiler,
                      AUTO  = AUTO,
                      chromatic = chromatic,
                      matrices  = matrices)

  #
  #   if(is.null(Nboot)){Nboot = 1}
  #
  #   out <- crqa_rp_measures(RM,
  #                            emRad= emRad,
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
  #   dfori <- dplyr::gather(out, key = measure, value = value)
  #
  #   col.ind <- dplyr::tbl_df(index(RM))
  #   row.ind <- dplyr::tbl_df(sample(index(RM),size=nrow(RM)))
  #
  #   if(Nboot>1){cat(paste0("Bootstrapping Recurrence Matrix... ",Nboot," iterations.\n"))
  #     bootout <- col.ind  %>%
  #       bootstrap(Nboot) %>%
  #       do(crqa_rp_measures(RM[row.ind,unlist(.)],
  #                            emRad= emRad,
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
  #     rqout <-  dfrepl %>% dplyr::group_by(measure) %>%
  #       summarize(
  #         val     = NA,
  #         ci.low  = stats::quantile(value, ci.lo, na.rm = TRUE),
  #         ci.high = stats::quantile(value, ci.hi, na.rm = TRUE),
  #         mean    = mean(value, na.rm = TRUE),
  #         sd      = stats::sd(value, na.rm = TRUE),
  #         var     = stats::var(value, na.rm = TRUE),
  #         N       = n(),
  #         se      = stats::sd(value, na.rm = TRUE)/sqrt(n()),
  #         median  = mean(value, na.rm = TRUE),
  #         mad     = stats::mad(value, na.rm = TRUE)
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



#' Diagonal Recurrence Profile
#'
#' @param RM A binary recurrence matrix
#' @param diagWin Window around the line of synchrony
#' @param xname Label for x-axis
#' @param yname Label for y-axis
#' @param DLmin Minimal diagonal line length (default = `2`)
#' @param VLmin Minimal vertical line length (default = `2`)
#' @param HLmin Minimal horizontal line length (default = `2`)
#' @param DLmax Maximal diagonal line length (default = length of diagonal -1)
#' @param VLmax Maximal vertical line length (default = length of diagonal -1)
#' @param HLmax Maximal horizontal line length (default = length of diagonal -1)
#' @param doShuffle Should a shuffled baseline be calculated (default = `FALSE`)
#' @param y1 The original `y1` time series
#' @param y2 The original `y2` time series
#' @param Nshuffle How many shuffled versions to make up the baseline? The default is `19`, which is the minimum for a one-sided surrogate test.
#' @param AUTO Auto-recurrence? (default = `FALSE`)
#' @param chromatic Force chromatic RQA? (default = `FALSE`)
#' @param matrices Return matrices? (default = `FALSE`)
#' @param doPlot Plot
#'
#' @return A plot and/or the data for the plot
#'
#' @export
#'
crqa_diagProfile <- function(RM,
                             diagWin = NULL,
                             xname = "X-axis",
                             yname = "Y-axis",
                             DLmin = 2,
                             VLmin = 2,
                             HLmin = 2,
                             DLmax = length(Matrix::diag(RM))-1,
                             VLmax = length(Matrix::diag(RM))-1,
                             HLmax = length(Matrix::diag(RM))-1,
                             doShuffle = FALSE,
                             y1        = NA,
                             y2        = NA,
                             Nshuffle  = 19,
                             AUTO      = NULL,
                             chromatic = FALSE,
                             matrices  = FALSE,
                             doPlot    = TRUE){


  if(!all(as.vector(RM)==0|as.vector(RM)==1)){
    stop("Expecting a binary (0,1) matrix.")
  }

  # check auto-recurrence and make sure Matrix has sparse triplet representation
  RM <- rp_checkfix(RM, checkAUTO = TRUE, fixAUTO = TRUE, checkTSPARSE = TRUE, fixTSPARSE = TRUE)

  if(is.null(AUTO)){
    AUTO <- attr(RM,"AUTO")
  }

  if(is.numeric(diagWin)){
    if(diagWin>0){
      diagWin <- seq(-diagWin,diagWin)
    } else {
      diagWin <- seq(-1*NCOL(RM),NCOL(RM))
    }
  } else {
    stop("Diagonal window must be 1 positive integer!!")
  }

  if(doShuffle){
    if(any(is.na(y1),is.na(y2))){
      stop("Need series y1 and y2 in order to do the shuffle!!")
    } else {
      cat("Calculating diagonal recurrence profiles... \n")
      TSrnd <- tseries::surrogate(x=y2, ns=Nshuffle, fft=FALSE, amplitude = FALSE)
      if(Nshuffle==1){
        TSrnd <- matrix(TSrnd,ncol=1)
      }
    }
    emDim <- attr(RM,"emDim")
    emLag <- attr(RM,"emLag")
    emRad <- attr(RM,"emRad")
  } else {
    Nshuffle <- 0
  }

  out <- vector(mode = "list", length = Nshuffle+1)

  if(doShuffle){
    names(out) <- c("obs",paste0("Shuffled", 1:Nshuffle))
  } else {
    names(out) <- "obs"
  }

  for(r in seq_along(out)){

    if(r==1){
      RMd <- RM
      rm(RM)
    } else {
      RMd <- rp(y1 = y1, y2 = TSrnd[,(r-1)],emDim = emDim, emLag = emLag, emRad = emRad, to.sparse = TRUE)
    }
    #rp_lineDist(RM,d = diagWin, matrices = TRUE)
    B <- rp_nzdiags(RMd, removeNZ = FALSE)
    rm(RMd)

    diagID <- 1:NCOL(B)
    names(diagID) <- colnames(B)
    if(length(diagWin)<NCOL(B)){
      cID <- which(colnames(B)%in%diagWin)
      B <- B[,cID]
      diagID <- seq_along(cID)
      names(diagID) <- colnames(B)
    }

    #winRR <- sum(B==1, na.rm = TRUE)
    df <- plyr::ldply(diagID, function(i){
      data.frame(index = i, RR = sum(B[,i]==1, na.rm = TRUE)/ (NROW(B)-abs(as.numeric(colnames(B)[i]))))
    }, .id = "Diagonal")

    df$group  <- 1
    df$labels <- paste(df$Diagonal)
   # df$labels[df$Diagonal==0] <- ifelse(AUTO,"LOI","LOS")
    df$labels <- factor(df$labels,levels = df$labels,ordered = TRUE)

    out[[r]] <- df
    rm(df,B,cID,diagID)

    cat(paste("\nProfile"),r)
  }

  dy      <- plyr::ldply(out)
  if(doShuffle){
    df_shuf <- dplyr::filter(dy, .data$.id!="obs")
  } else {
    df_shuf <- dy
  }

  dy_m <- df_shuf %>%
    dplyr::group_by(.data$Diagonal,.data$labels) %>%
    dplyr::summarise(meanRRrnd = mean(.data$RR), sdRRrnd = stats::sd(.data$RR))

    if(Nshuffle==1){
    dy_m$sdRRrnd <- dy_m$meanRRrnd
  }

  dy_m$ciHI <- dy_m$meanRRrnd + 1.96*(dy_m$sdRRrnd/sqrt(Nshuffle))
  dy_m$ciLO <- dy_m$meanRRrnd - 1.96*(dy_m$sdRRrnd/sqrt(Nshuffle))

  obs <- dy$RR[dy$.id=="obs"]
  if(length(obs)!=(2*length(y1))){
    obs <-  ts_trimfill(obs,c(y1,y2),action = 'trim')$y
  }
  dy_m$y_obs <- obs

  df <- tidyr::gather(dy_m, key = "variable", value = "RR", -c(.data$Diagonal,.data$sdRRrnd, .data$labels, .data$ciLO, .data$ciHI))
  df$Diagonal <- as.numeric(df$Diagonal)

  if(doPlot){

    Diags <- as.numeric(levels(df$labels))
    Diags[is.na(Diags)] <- 0
    if(length(diagWin)>21){
      ext <- max(min(abs(Diags),na.rm = TRUE),abs(max(Diags,na.rm = TRUE)))
      breaks <- which(Diags%in%(c(seq(-ext,-1,length.out = 10),0,seq(1,ext,length.out = 10))))
      labels <- sort(unique(c(Diags[breaks],0)))
      breaks <- sort(unique(c(breaks,stats::median(breaks))))
    } else {
      breaks <- seq_along(labels)
      labels <- sort(unique(c(Diags[breaks],0)))

    }

    x1<-(which.min(as.numeric(paste(df$Diagonal))))+(length(diagWin)*.1)
    x2<-(which.max(as.numeric(paste(df$Diagonal))))-(length(diagWin)*.1)
    yL<-max(as.numeric(paste(df$RR)),na.rm = TRUE)+0.1
    col <- c("ciHI" = "grey70", "ciLO" = "grey70", "meanRRrnd" = "grey40","y_obs" = "black")
    siz <- c("ciHI" = .5, "ciLO" = .5, "meanRRrnd" = .5,"y_obs" = 1)

    g <- ggplot2::ggplot(df, ggplot2::aes_(x=~Diagonal)) +
      ggplot2::geom_ribbon(ggplot2::aes_(ymin=~ciLO, ymax=~ciHI), alpha=0.3) +
      ggplot2::geom_line(ggplot2::aes_(y=~RR, colour = ~variable), size = .5) +
      #ggplot2::geom_line(ggplot2::aes_(y=~y_obs), colour = "black", size = 1) +
      ggplot2::geom_vline(xintercept = which(df$labels%in%c("LOS","LOI")), size=1, colour = "grey50") +
      ggplot2::scale_y_continuous("Recurrence Rate",limits = c(0,yL)) +
      ggplot2::scale_x_discrete("Diagonals in recurrence Matrix", breaks = breaks, labels = labels) +
      ggplot2::geom_label(x=x1,y=yL,label=paste0("Recurrences due to\n ",xname),hjust="left", inherit.aes = FALSE) +
      ggplot2::geom_label(x=x2,y=yL,label=paste0("Recurrences due to\n ",yname),hjust="right", inherit.aes = FALSE) +
       ggplot2::scale_colour_manual(values = col) +
      # ggplot2::scale_size_manual(values = siz) +
      ggplot2::theme_bw()

    print(g)

    if(doShuffle){
      df <- tidyr::spread(df,key = .data$variable, value = .data$RR)
    }

    return(invisible(list(plot = g, data = df)))

  } else {

    if(doShuffle){df <- tidyr::spread(df,key = .data$variable, value = .data$RR)}
    return(df)

  }
}




#' Estimate embedding lag (tau)
#'
#' A wrapper for nonlinearTseries::timemLag
#'
#' @param y Time series or numeric vector
#' @param selection.methods Selecting an optimal embedding lag (default: Return "first.e.decay", "first.zero", "first.minimum", "first.value", where value is 1/exp(1))
#' @param maxLag Maximal lag to consider (default: 1/4 of timeseries length)
#' @param ... Additional parameters
#'
#' @return The ami function with requested minima
#'
#' @export
#'
#' @family Estimate CRQA parameters
#'
est_emLag <- function(y,
                selection.methods = "first.minimum",
                maxLag = length(y)/4,
                ...){

  y   <- y[!is.na(y)]
   tmi <- nonlinearTseries::mutualInformation(y,
                                              lag.max = maxLag,
                                              #n.partitions = floor((max(y)-min(y))/(2*stats::IQR(y)*length(y)^(-1/3))),
                                              do.plot = FALSE
                                             )

  lags <- numeric(length=length(selection.methods))
  cnt <- 0
  for(sm in selection.methods){
    cnt <- cnt + 1
    lag <-try_CATCH(nonlinearTseries::timeLag(y,
                                              technique = "ami",
                                              selection.method = sm,
                                              lag.max = maxLag,
                                              do.plot = FALSE,
                                              #n.partitions = floor((max(y)-min(y))/(2*stats::IQR(y)*length(y)^(-1/3)))
                                              ))
    if(any(grepl("Error",lag$value))){lags[cnt] <- NA} else {lags[cnt] <- lag$value}
  }

  #id.peaks <- find_peaks(tmi$mutual.information, m = 3, wells = TRUE)
  id.min   <- tmi$time.lag[which.min(tmi$mutual.information)]
  #as.numeric(tmi$mutual.information[id.peaks[id.min]])

  out <-   cbind.data.frame(selection.method = c(selection.methods,
                                                 "global.minimum",
                                                 "maximum.lag"),
                            optimal.lag = c(lags, id.min, maxLag)
  )

  #tmi$mutual.information[tmi$time.lag%in%out$optimal.lag]
  out$ami <- sapply(out$optimal.lag, function(r){
    if(is.na(r)){
      NA
      } else {
        tmi$mutual.information[r]
        }
    })

  #tout <- summarise(dplyr::group_by(out, opLag, ami), selection.method = paste0(unique(selection.method), collapse="|")

  return(out)
}


#' Estimate number of embedding dimensions
#'
#' @param y Time series or numeric vector
#' @param delay Embedding lag
#' @param maxDim Maximum number of embedding dimensions
#' @param threshold See [nonlinearTseries::estimateEmbeddingDim()]
#' @param max.relative.change See [nonlinearTseries::estimateEmbeddingDim()]
#' @param doPlot Plot
#' @param ... Other arguments (not in use)
#'
#' @description A wrapper for nonlinearTseries::estimateEmbeddingDim
#'
#' @return Embedding dimensions
#' @export
#'
est_emDim <- function(y, delay = est_emLag(y), maxDim = 15, threshold = .95, max.relative.change = .1, doPlot = FALSE, ...){
  cbind.data.frame(EmbeddingLag   = delay,
                   EmbeddingDim   = nonlinearTseries::estimateEmbeddingDim(y,
                                                         time.lag  = delay,
                                                         threshold = threshold,
                                                         max.relative.change = max.relative.change,
                                                         max.embedding.dim = maxDim,
                                                         do.plot = doPlot)
  )
}

#' rp_nzdiags
#'
#' @description Get all nonzero diagonals of a binary matrix, or, diagonals specified as a vector by argument `d`.
#'
#' @param RM A binary (0,1) matrix.
#' @param d An optional vector of diagonals to extract.
#' @param returnVectorList Return list
#' @param returnNZtriplets Return a dataframe with coordinates of only nonzero elements in diagonals (default = `FALSE`)
#' @param removeNZ Remove nonzero diagonals if `TRUE`. If `FALSE` returns the full diagonals matrix. Use e.g. to plot diagonal recurrence profiles (default = `TRUE`)
#' @param silent Silent-ish mode
#'
#' @author Fred Hasselman
#'
#' @return A matrix object with nonzero diagonals as columns and/or a dataframe with coordinates of nonzero diagonal elements
#'
#' @export
#'
#' @family Distance matrix operations (recurrence plot)
#'
rp_nzdiags <- function(RM=NULL, d=NULL, returnVectorList=TRUE, returnNZtriplets=FALSE, removeNZ = TRUE,silent = TRUE){
  # Loosely based on MATLAB function spdiags() by Rob Schreiber - Copyright 1984-2014 The MathWorks, Inc.

  if(grepl("matrix",class(RM),ignore.case = TRUE)){

    if(all(RM>0)){warning("All matrix elements are nonzero.")}

    s  <- Sys.time()

    nzdiagsM <- methods::as(RM, "dgTMatrix")
    nzdiags  <- data.frame(row   = nzdiagsM@i,
                           col   = nzdiagsM@j,
                           value = nzdiagsM@x,
                           ndiag = (nzdiagsM@j)-(nzdiagsM@i))
    nzdiags <- dplyr::arrange(nzdiags,nzdiags$ndiag)

    if(removeNZ){
      if(!is.null(d)){
        nd <- unique(nzdiags$ndiag)
        # Get diagonals which have nonzero elements
        d <- nd[nd%in%sort(as.vector(d))]
      } else {
        d <- unique(nzdiags$ndiag)
      }
    } else {
      d <-  c(-1*rev(1:(NCOL(nzdiagsM)-1)),0,1:(NCOL(nzdiagsM)-1))
    }

    indMat <- col(RM)-row(RM)
    indMatV <- as.vector(indMat)

    #cat("\nStart extraction of nonzero diagonals\n")
    m  <- NROW(RM)
    n  <- NCOL(RM)
    p  <- length(d)
    if(is.logical(RM)){
      B <- matrix(FALSE,nrow = min(c(m,n), na.rm = TRUE), ncol = p)
    } else {
      B <- matrix(0, nrow = min(c(m,n), na.rm = TRUE), ncol = p)
    }
    colnames(B) <- paste(d)

    for(i in seq_along(d)){
      B[(nzdiags$row[nzdiags$ndiag==d[i]])+1,i] <- nzdiags$value[nzdiags$ndiag==d[i]]
    }

    zID <- which(Matrix::colSums(Matrix::as.matrix(B))==0)
    if(length(zID)>0){
      if(removeNZ){
        B  <- B[,-zID]
        nzdiags <- nzdiags[nzdiags$ndiag%in%zID,]
      }
    }
    #dl <- plyr::llply(seq_along(toList),function(dd){RM[indMat==toList[[dd]]]}, .progress = plyr::progress_text(char = "o~o"))
    e  <- Sys.time()

    if(!silent){cat(paste0("\nTime: ", round(difftime(e,s, units = "mins"), digits=1)," min.\n"))}

    if(returnNZtriplets){
      return(list(B = B, nzdiags = nzdiags))
    } else {
      return(B)
    }
  } else {
    stop("Input to rp_nzdiags is not a matrix-like object.")
  }
}



rp_nzdiags_matlab <- function(RP,d=NULL){
  # Loosely based on MATLAB function spdiags() by Rob Schreiber - Copyright 1984-2014 The MathWorks, Inc.
  #require(Matrix)

  if(grepl("matrix",class(RP),ignore.case = TRUE)){

    A <- RP

    if(all(A>0)){warning("All matrix elements are nonzero.")}

    # create an indicator for all diagonals in the matrix
    ind   <- col(A)-row(A)

    # Split the matrix!
    spd <- split(A, ind)

    if(is.null(d)){

      # Get diagonals which have nonzero elements
      keepID <-plyr::ldply(spd, function(di) any(di>0))
      nzdiag <- spd[keepID$V1]
      # Indices of nonzero diagonals
      d      <- as.numeric(keepID$.id[keepID$V1])

    } else {

      # Diagonals are specified
      d <- sort(as.vector(d))
      keepID <- names(spd)%in%d
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


rp_nzdiags_chroma <- function(RP, d=NULL){
  # Loosely based on MATLAB function spdiags() by Rob Schreiber - Copyright 1984-2014 The MathWorks, Inc.
  #require(Matrix)

  if(grepl("matrix",class(RP),ignore.case = TRUE)){

    A <- RP

    if(all(A>0)){warning("All matrix elements are nonzero.")}

    # create an indicator for all diagonals in the matrix
    ind   <- col(A)-row(A)

    # Split the matrix!
    spd <- split(A, ind)

    if(is.null(d)){

      # Get diagonals which have nonzero elements
      keepID <-plyr::ldply(spd, function(di) any(di>0))
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


#' Line length distributions
#'
#' @param RM A thresholded recurrence matrix (binary: 0 - 1)
#' @param DLmin Minimal diagonal line length (default = `2`)
#' @param VLmin Minimal vertical line length (default = `2`)
#' @param HLmin Minimal horizontal line length (default = `2`)
#' @param DLmax Maximal diagonal line length (default = length of diagonal -1)
#' @param VLmax Maximal vertical line length (default = length of diagonal -1)
#' @param HLmax Maximal horizontal line length (default = length of diagonal -1)
#' @param d Vector of diagonals to be extracted from matrix `RP` before line length distributions are calculated. A one element vector will be interpreted as a windowsize, e.g., `d = 50` will extract the diagonal band `-50:50`. A two element vector will be interpreted as a band, e.g. `d = c(-50,100)` will extract diagonals `-50:100`. If `length(d) > 2`, the numbers will be interpreted to refer to individual diagonals, `d = c(-50,50,100)` will extract diagonals `-50,50,100`.
#' @param theiler Size of the theiler window, e.g. `theiler = 1` removes diagonal bands -1,0,1 from the matrix. If `length(d)` is `NULL`, 1 or 2, the theiler window is applied before diagonals are extracted. The theiler window is ignored if `length(d)>2`, or if it is larger than the matrix or band indicated by parameter `d`.
#' @param invert Relevant for Recurrence Time analysis: Return the distribution of 0 valued segments in nonzero diagonals/verticals/horizontals. This indicates the time between subsequent line structures.
#' @param AUTO Is this an AUTO RQA?
#' @param chromatic Chromatic RQA?
#' @param matrices Return the matrices ?
#'
#'
#' @description Extract lengths of diagonal, vertical and horizontal line segments from a recurrence matrix.
#'
#' @details Based on the Matlab function `linedists` by Stefan Schinkel, Copyright (C) 2009 Stefan Schinkel, University of Potsdam, http://www.agnld.uni-potsdam.de
#'
#' References:
#' S. Schinkel, N. Marwan, O. Dimigen & J. Kurths (2009):
#' "Confidence Bounds of recurrence-based complexity measures
#' Physics Letters A,  373(26), pp. 2245-2250
#'
#' Copyright (C) 2009 Stefan Schinkel, University of Potsdam
#' <http://www.agnld.uni-potsdam.de>
#'
#' @author Fred Hasselman
#' @return A list object with distributions of line lengths. If `matrices = TRUE` datafr are returned whose columns represent the nonzero diagonals, verticals, or, horizontals.
#'
#' @export
#'
#' @family Distance matrix operations (recurrence plot)
#'
rp_lineDist <- function(RM,
                        DLmin = 2,
                        VLmin = 2,
                        HLmin = 2,
                        DLmax = length(Matrix::diag(RM))-1,
                        VLmax = length(Matrix::diag(RM))-1,
                        HLmax = length(Matrix::diag(RM))-1,
                        d         = NULL,
                        theiler   = NULL,
                        invert    = FALSE,
                        AUTO      = NULL,
                        chromatic = FALSE,
                        matrices  = FALSE){

  # For boot()
  # RP <- RP[indices,]

  if(!all(as.vector(RM)==0|as.vector(RM)==1)){stop("Matrix should be a binary (0,1) matrix!!")}

  if(invert){
    RM <- Matrix::Matrix(1-RM,sparse = TRUE)
    if(Matrix::isSymmetric(Matrix::unname(RM))){
      RM <- bandReplace(RM,0,0,1)
    }
  }


  if(!is.null(d)){
    if(length(d)==1){d <- -d:d}
    if(length(d)==2){d <-  d[1]:d[2]}
  }
  if(!is.null(theiler)){
    if(length(d)<length(-theiler:theiler)){warning("Ignoring theiler window...")}
    RM <- bandReplace(RM,-theiler,theiler,0)
  }

  if(Matrix::isSymmetric(Matrix::unname(RM))){
    if(all(Matrix::diag(RM)==1)){
      RP <- bandReplace(RM,0,0,0)
    }
  } else {
    RP <- RM
  }

  B <- rp_nzdiags(RP)
  V <- Matrix::as.matrix(RP)[,colSums(Matrix::as.matrix(RP))>0]

  if(Matrix::isSymmetric(Matrix::unname(RM))){
    if(all(Matrix::diag(RM)==1)){
      RP <- bandReplace(Matrix::t(RM),0,0,0)
    }
  } else {
    RP <- Matrix::t(RM)
  }

  H <- Matrix::as.matrix(RP)[,colSums(Matrix::as.matrix(RP))>0]
  rm(RP)

  # Get diagonal lines & pad with zeros
  diagonals   <- rbind.data.frame(rep(0,dim(B)[2]),
                                  B,
                                  rep(0,dim(B)[2])
  )

  # get nonzero vertical Lines & pad with zeros
  verticals <- rbind.data.frame(rep(0,dim(V)[2]),
                                V,
                                rep(0,dim(V)[2])
  )
  colnames(verticals) <- paste(1:ncol(verticals))

  # get nonzero horizontal Lines & pad with zeros
  horizontals <- rbind.data.frame(rep(0,dim(H)[2]),
                                  H,
                                  rep(0,dim(H)[2])
  )
  colnames(horizontals) <- paste(1:ncol(horizontals))

  # Get indices of line lengths
  diagonals.ind   <- tidyr::gather(diagonals,   key = "diagonal",   value = "segment")
  verticals.ind   <- tidyr::gather(verticals,   key = "vertical",   value = "segment")
  horizontals.ind <- tidyr::gather(horizontals, key = "horizontal", value = "segment")

  D <- diagonals.ind$segment
  names(D) <- paste0(diagonals.ind$diagonal,ifelse(invert,"DT","D"))
  V <- verticals.ind$segment
  names(V) <- paste0(verticals.ind$vertical,ifelse(invert,"VT","V"))
  H <- horizontals.ind$segment
  names(H) <- paste0(horizontals.ind$horizontal,ifelse(invert,"HT","H"))

  # Get consecutive nonzero segments from indices, their difference is the segment length
  # We added a row of 0s so we'll get sequences of -1, 1, -1
  diagonals.dist   <- sort(which(diff(D)==-1)-which(diff(D)==1))
  verticals.dist   <- sort(which(diff(V)==-1)-which(diff(V)==1))
  horizontals.dist <- sort(which(diff(H)==-1)-which(diff(H)==1))

  diagonals.dist   <- diagonals.dist[diagonals.dist%[]%c(DLmin,DLmax)]
  verticals.dist   <- verticals.dist[verticals.dist%[]%c(VLmin,VLmax)]
  horizontals.dist <- horizontals.dist[horizontals.dist%[]%c(HLmin,HLmax)]

  if(length(diagonals.dist)==0){diagonals.dist <- NA}
  if(length(verticals.dist)==0){verticals.dist <- NA}
  if(length(horizontals.dist)==0){horizontals.dist <- NA}

  return(list(diagonals.dist   = diagonals.dist,
              verticals.dist   = verticals.dist,
              horizontals.dist = horizontals.dist,
              diagonals.mat = diagonals[-c(1,NROW(diagonals)),],
              verticals.mat = verticals[-c(1,NROW(verticals)),],
              horizontals.mat = horizontals[-c(1,NROW(horizontals)),])[c(TRUE,TRUE,TRUE,matrices,matrices,matrices)]
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


#' Create a Distance Matrix
#'
#' @param y1 A numeric vector or time series
#' @param y2 A numeric vector or time series for cross recurrence
#' @param emDim The embedding dimensions
#' @param emLag The embedding lag
#' @param emRad The threshold (emRad) to apply to the distance matrix to create a binary or weighted matrix. If `NULL`, an unthresholded matrix will be created (default = `NULL`)
#' @param to.ts Should `y1` and `y2` be converted to time series objects?
#' @param order.by If `to.ts = TRUE`, pass a vector of the same length as `y1` and `y2`. It will be used as the time index, if `NA` the vector indices will be used to represent time.
#' @param to.sparse Should sparse matrices be used?
#' @param weighted If `FALSE` a binary matrix will be returned. If `TRUE` every value larger than `emRad` will be `0`, but values smaller than `emRad` will be retained (default = `FALSE`)
#' @param method Distance measure to use. Any option that is valid for argument `method` of [proxy::dist()]. Type `proxy::pr_DB$get_entries()` to se a list of all the options. Common methods are: "Euclidean", "Manhattan", "Minkowski", "Chebysev" (or the same but shorter: "L2","L1","Lp" and "max" distance) (default = `"Euclidean"`)
#' @param targetValue A value passed to `crqa_radius(...,type="fixed", targetMeasure="RR")` if `is.na(emRad)==TRUE`.
#' @param doPlot Plot the matrix by calling [rp_plot()] with defult settings
#' @param silent Silent-ish mode
#' @param ... Any paramters to pass to [rp_plot()] if `doPlot = TRUE`
#'
#' @return A (Coss-) Recurrence matrix with attributes:
#' \enumerate{
#' \item `emdims1` and `emdims2` - A matrix of surrogate dimensions
#' \item `emdims1.name` and `emdims2.name` - Names of surrogate dimensions
#' \item `method` and `call` - The distance `method` used by [proxy::dist()]
#' \item `weighetd` - Whether a weighted matrix is returned
#' \item `emDim`, `emLag` and `emRad` - The embedding parameters
#' \item `AUTO` - Whether the matrix represents AUTO recurrence
#'  }
#'
#'
#' @export
#'
#' @family Distance matrix operations (recurrence plot)
#'
#'
rp <- function(y1, y2 = NULL,
                   emDim = 1,
                   emLag = 1,
                   emRad = NULL,
                   to.ts = NULL,
                   order.by = NULL,
                   to.sparse = FALSE,
                   weighted = FALSE,
                   method = "Euclidean",
                   targetValue  = .05,
                   doPlot = FALSE,
                   silent = TRUE,
                   ...){

  if(is.null(y2)){
    y2 <- y1
    attributes(y2) <- attributes(y1)
    }

  if(!is.data.frame(y1)){
    id <- deparse(substitute(y1))
    y1 <- as.data.frame(y1)
    if(length(id)==NCOL(y1)){
     colnames(y1) <- id
    }
  }

  if(!is.data.frame(y2)){
    id <- deparse(substitute(y2))
    y2 <- as.data.frame(y2)
    if(length(id)==NCOL(y2)){
    colnames(y2) <- id
    }
  }

  et1 <- ts_embed(y1, emDim, emLag, silent = silent)
  et2 <- ts_embed(y2, emDim, emLag, silent = silent)

  dist_method <- try_CATCH(proxy::pr_DB$get_entry(method))
  if("error"%in%class(dist_method$value)){
    stop("Unknown distance metric!\nUse proxy::pr_DB$get_entries() to see a list of valid options.")
  } else {
  dmat <- proxy::dist(x = et1,
                      y = et2,
                      method = method,
                      diag = ifelse(identical(et1,et2),FALSE,TRUE))
  }

  # Remove proxy class
  dmat <- unclass(dmat)

  # Check a time index was requested
  if(!is.null(to.ts)){
    if(is.null(order.by)){
      order.by <- lubridate::as_datetime(1:NROW(dmat), origin = lubridate::ymd_hms(Sys.time()))
    }
    dmat <-  switch(to.ts,
                    "xts" =  xts::xts(dmat, order.by = lubridate::as_datetime(order.by)),
                    "zoo" =  zoo::zoo(dmat, order.by = lubridate::as_datetime(order.by))
                    )
  }
  # Add names if ordered by time vector
  if(!is.null(order.by)){
    colnames(dmat) <- paste(order.by)
    rownames(dmat) <- paste(order.by)
  }

    if(!is.null(emRad)){
      if(is.na(emRad)){
        emRad <- crqa_radius(RM = dmat, emDim = emDim, emLag = emLag, targetValue = targetValue)$Radius
      }
      if(weighted){
        dmat <- di2we(dmat, emRad = emRad, convMat = to.sparse)
      } else {
        dmat <- di2bi(dmat, emRad = emRad, convMat = to.sparse)
      }
    }

  if(to.sparse){
    attributes(dmat)$emDims1  <- et1
    attributes(dmat)$emDims2  <- et2
    attributes(dmat)$emDims1.name <- colnames(y1)
    attributes(dmat)$emDims2.name <- colnames(y2)
    attributes(dmat)$weighted <- weighted
    attributes(dmat)$emLag <- emLag
    attributes(dmat)$emDim <- emDim
    attributes(dmat)$emRad <- emRad%00%NA
  } else {
    attr(dmat,"emDims1") <- et1
    attr(dmat,"emDims2") <- et2
    attr(dmat,"emDims1.name") <- colnames(y1)
    attr(dmat,"emDims2.name") <- colnames(y2)
    attr(dmat,"weighted") <- weighted
    attr(dmat,"emLag") <- emLag
    attr(dmat,"emDim") <- emDim
    attr(dmat,"emRad") <- emRad%00%NA
  }

  dmat <- rp_checkfix(dmat, checkAUTO = TRUE, fixAUTO = TRUE)

  if(doPlot){
    dotArgs  <- formals(rp_plot)
    nameOk   <- rep(TRUE,length(dotArgs))
    if(...length()>0){
      dotArgs <- list(...)
      nameOK  <- names(dotArgs)%in%methods::formalArgs(rp_plot)
      # Plot with defaults
      if(!all(nameOK)){
        dotArgs    <- formals(rp_plot)
        nameOk <- rep(TRUE,length(dotArgs))
      }
    }
    dotArgs$RM <- dmat
    do.call(rp_plot, dotArgs[nameOk])
  }

  return(dmat)
}


#' Copy Matrix Attributes
#'
#' Simple attribute copy used in `casnet` to convert between `matrix` and `Matrix` classes and back.
#'
#' @param source Source matrix
#' @param target Target matrix
#' @param source_remove Remove these attribute fields from the source before copying.
#'
#' @return The target matrix with attributes copied deom the source matrix.
#' @export
#'
rp_copy_attributes <- function(source, target, source_remove = c("names", "row.names", "class","dim", "dimnames","x")){
  attrList <- attributes(source)
  attrList[source_remove] <- NULL
  attributes(target) <- c(attributes(target), attrList)
  return(target)
}

#' Check a Recurrence Matrix
#'
#' @param RM RM
#' @param checkS4 checkS4
#' @param checkAUTO checkAUTO
#' @param checkSPARSE checkSPARSE
#' @param fixS4 fixS4
#' @param fixAUTO fixAUTO
#' @param fixSPARSE fixSPARSE
#'
#' @return A checked and/or fixed matrix
#' @export
#' @keywords internal
#'
rp_checkfix <- function(RM, checkS4 = TRUE, checkAUTO = TRUE, checkSPARSE = FALSE, checkTSPARSE = FALSE, fixS4 = FALSE, fixAUTO = TRUE, fixSPARSE = FALSE, fixTSPARSE = FALSE){

  dummy <- Matrix::Matrix(matrix(c(0,0)))
  dummy <- rp_copy_attributes(source = RM, target =dummy)

  # Always check S4
  if(!checkS4){checkS4 <- TRUE}

  yesS4      <- FALSE
  yesAUTO    <- FALSE
  yesSPARSE  <- FALSE
  yesTSPARSE <- FALSE

  if(checkS4){
    yesS4 <- isS4(RM)
  }
  if(!yesS4&fixS4){
    RM <- Matrix::Matrix(RM)
  }

  # check auto-recurrence
  if(checkAUTO){
    if(yesS4){
      yesAUTO <- Matrix::isSymmetric(RM)
    } else {
      yesAUTO <- identical(as.vector(RM[lower.tri(RM)]),as.vector(t(RM)[lower.tri(RM)]))
    }
  }
  if(fixAUTO){
    attributes(RM)$AUTO    <- yesAUTO
    attributes(dummy)$AUTO <- yesAUTO
  }

  if(checkTSPARSE){
    if(class(RM)%in%names(methods::getClass("TsparseMatrix")@subclasses)){
      yesTSPARSE <- TRUE
    }
  }
  if(fixTSPARSE&!yesTSPARSE){
      if(!yesS4){
        RM <- Matrix::Matrix(RM,sparse = TRUE)
      }
    Mtype <- gsub("CMatrix","TMatrix",class(RM))
    eval(parse(text=paste0("RM <- as(RM,'",Mtype,"')")))
  }


 RM <- rp_copy_attributes(source = dummy, target = RM, source_remove = c("Dimnames", "i", "class","Dim", "p","x","factors"))

 return(RM)
}


#' Plot (thresholded) distance matrix as a recurrence plot
#'
#' @param RM A distance matrix or recurrence matrix
#' @param plotDimensions Should the state vectors be plotted if they are available as attributes of RM (default = `TRUE`)
#' @param plotMeasures Print common (C)RQA measures in the plot if the matrix is binary
#' @param plotRadiusRRbar The `Radius-RR-bar` is a colour-bar guide plotted with an unthresholded distance matrix indicating a number of `RR` values one would get if a certain distance threshold were chosen (`default = TRUE`)
#' @param drawGrid Draw a grid on the recurrence plot (`default = FALSE`)
#' @param markEpochsLOI Pass a factor whose levels indicate different epochs or phases in the time series and use the line of identity to represent the levels by different colours (`default = NULL`)
#' @param Chromatic If `TRUE` and there are more than two discrete values in `RM`, give recurrent points a distinct colour. If `RM` was returned by `crqa_rp(..., chromatic = TRUE)`, the recurrence plot will colour-code recurrent points according to the category values in `attributes(RM)$chromaticRP` (`default = FALSE`)
#' @param radiusValue If `plotMeasures = TRUE` and RM is an unthresholded matrix, this value will be used to calculate recurrence measures. If `plotMeasures = TRUE` and RM is already a binary recurence matrix, pass the radius that was used as a threshold to create the matrix for display purposes. If `plotMeasures = TRUE` and `radiusValue = NA`, function `crqa_radius()` will be called with default settings (find a radius that yields .05 recurrence rate). If `plotMeasures = FALSE` this setting will be ignored.
#' @param title A title for the plot
#' @param xlab An x-axis label
#' @param ylab An y-axis label
#' @param plotSurrogate Should a 2-panel comparison plot based on surrogate time series be added? If `RM` has attributes `y1` and `y2` containing the time series data (i.e. it was created by a call to [rp()]), the following options are available: "RS" (random shuffle), "RP" (randomised phases), "AAFT" (amplitude adjusted fourier transform). If no timeseries data is included, the columns will be shuffled.  NOTE: This is not a surrogate test, just 1 surrogate is created from `y1`.
#' @param returnOnlyObject Return the ggplot object only, do not draw the plot (default = `TRUE`)
#'
#' @return A nice plot of the recurrence matrix.
#' @export
#'
#' @family Distance matrix operations (recurrence plot)
#'
rp_plot <- function(RM,
                    plotDimensions = FALSE,
                    plotMeasures   = FALSE,
                    plotRadiusRRbar = TRUE,
                    drawGrid = FALSE,
                    markEpochsLOI = NULL,
                    Chromatic = NULL,
                    radiusValue = NA,
                    title = "", xlab = "", ylab="",
                    plotSurrogate = NA,
                    returnOnlyObject = FALSE){

  useGtable <- TRUE

  # # check patchwork
  # if(!length(find.package("patchwork",quiet = TRUE))>0){
  #   warning("Package patchwork is not installed...\n1. Install Xcode from App Store (MacOS) or rwintools.exe from CRAN (Windows) \n2. Install patchwork: devtools::install_github('thomasp85/patchwork')\n3. Install casnet: devtools::install_github('FredHasselman/casnet')\n....Using gtable instead, with limited options\n")
  #   useGtable=TRUE
  # }


  colvec <- c("#FFFFFF","#000000")
  names(colvec) <- c("0","1")

  # check auto-recurrence and make sure Matrix has sparse triplet representation
  RM   <- rp_checkfix(RM, checkAUTO = TRUE)
  AUTO <- attr(RM,"AUTO")

  # prepare data
  if(attr(RM,"package")%00%""%in%"Matrix"){
    RM     <- rp_checkfix(RM, checkTSPARSE = TRUE, fixTSPARSE = TRUE)
    meltRP <- data.frame(Var1 = (RM@i+1), Var2 = (RM@j+1), value = as.numeric(RM@x))
  } else {
    meltRP <- reshape2::melt(as.matrix(RM))
  }

  hasNA <- FALSE
  if(any(is.na(meltRP$value))){
    hasNA <- TRUE
    meltRP$value[is.na(meltRP$value)] <- max(meltRP$value, na.rm = TRUE) + 1
  }

  # check unthresholded
  showL <- FALSE
  if(!all(as.vector(meltRP$value[!is.na(meltRP$value)])==0|as.vector(meltRP$value[!is.na(meltRP$value)])==1)){
    unthresholded <- TRUE

    if(!is.null(markEpochsLOI)){
      warning("Can't show epochs on an unthresholded Recurrence Plot!")
    }

  } else {

    unthresholded <- FALSE

    if(!is.null(attr(RM,"emRad"))){
      radiusValue <- attr(RM,"emRad")
    }

    # Check epochs
    if(!is.null(markEpochsLOI)){
      if(is.factor(markEpochsLOI)&length(markEpochsLOI)==max(c(NROW(RM),NCOL(RM)))){

        start <- max(meltRP$value, na.rm = TRUE) + 1

        cpal <- paletteer::paletteer_d(package = "rcartocolor",palette = "Safe", n = nlevels(markEpochsLOI))
        #cpal <- paletteer::paletteer_d(package = "ggthemes",palette = "tableau_colorblind10",n = nlevels(markEpochsLOI),direction = 1)

        if(hasNA){
        colvec <- c("#FFFFFF","#000000","#FF0000", cpal)
        names(colvec) <- c("0","1","NA",levels(markEpochsLOI))

        } else {
          colvec <- c("#FFFFFF","#000000", cpal) #viridis::viridis_pal()(nlevels(markEpochsLOI)))
          names(colvec) <- c("0","1",levels(markEpochsLOI))
        }

        N <- max(c(NROW(RM),NCOL(RM)))
        for(i in 1:N){
          j <- i
          meltRP$value[meltRP$Var1==i&meltRP$Var2==j] <- start + as.numeric_character(markEpochsLOI)[i]
        }

        meltRP$value <- factor(meltRP$value, levels = sort(unique(meltRP$value)), labels = names(colvec))
        showL <- TRUE

      } else {
        warning("Variable passed to 'markEpochsLOI' is not a factor or doesn't have correct length.")
      }
    } else {
      colvec <- c("#FFFFFF","#000000")
      names(colvec) <- c("0","1")
      meltRP$value <- factor(meltRP$value, levels = c(0,1), labels = c("0","1"))
    }
  }

  # Get CRQA measures
  if(plotMeasures){
    if(is.na(radiusValue)){
      if(!is.null(attributes(RM)$emRad)){
        radiusValue <- attr(RM,"emRad")
      } else {
        radiusValue <- crqa_radius(RM,silent = TRUE)$Radius
      }
    }
    if(unthresholded){
      rpOUT   <- crqa_rp(RM, emRad = radiusValue, AUTO = AUTO)
    } else {
      rpOUT   <- crqa_rp(RM, AUTO = AUTO)
    }
  }

  # main plot ----
  gRP <-  ggplot2::ggplot(ggplot2::aes_(x=~Var1, y=~Var2, fill = ~value), data= meltRP) +
    ggplot2::geom_raster(hjust = 0, vjust=0, show.legend = showL) +
    ggplot2::geom_abline(slope = 1,colour = "grey50", size = 1)

  if(unthresholded){
    gRP <- gRP + ggplot2::scale_fill_gradient2(low = "red3",
                                      high     = "steelblue",
                                      mid      = "white",
                                      na.value = scales::muted("slategray4"),
                                      midpoint = mean(meltRP$value, na.rm = TRUE),
                                      limit    = c(min(meltRP$value, na.rm = TRUE),max(meltRP$value, na.rm = TRUE)),
                                      space    = "Lab",
                                      name     = "")

    if(plotRadiusRRbar){
      # Create a custom legend ---
      distrange  <- round(seq(0,max(RM,na.rm = TRUE),length.out=7),2)
      resol      <- sort(unique(round(as.vector(RM),2)))
      if(length(resol)<7){
        resol <- distrange
      }
      if(length(resol)>100){
        resol <- round(seq(0,max(RM,na.rm = TRUE),length.out=100),2)
      }
      resol <- resol %>% tibble::as.tibble() %>% dplyr::mutate(y= seq(exp(0),exp(1),length.out=NROW(resol)), x=0.5)
      #resol <- resol[-1,]

      distrange <- plyr::ldply(c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5), function(t){
        suppressWarnings(crqa_radius(RM,targetValue = t,silent = TRUE, maxIter = 100, radiusOnFail = "percentile"))
      })
      #ldply(distrange[2:6],function(d) cbind(epsilon=d,RR=crqa_rp(RM = RM, emRad = d)$RR))

      RecScale <- data.frame(RR=distrange$Measure,epsilon=distrange$Radius)
      RecScale <- RecScale %>%
        dplyr::add_row(epsilon=mean(c(0,distrange$Radius[1])),RR=mean(c(0,distrange$Measure[1])),.before = 1) %>%
        dplyr::add_row(epsilon=max(RM),RR=1)

      resol$y <- elascer(x = resol$y,lo = min(log(RecScale$RR),na.rm = TRUE), hi = max(log(RecScale$RR),na.rm = TRUE))
      #resol$value <- log(resol$value)
      resol <- resol[-1,]

      gDist <-  ggplot2::ggplot(resol,ggplot2::aes_(x=~x,y=~y,fill=~value)) +
        ggplot2::geom_tile(show.legend = FALSE) +
        ggplot2::scale_y_continuous(name = "Recurrence Rate", breaks = log(RecScale$RR), labels = paste(round(RecScale$RR,3)), sec.axis = dup_axis(name=expression(paste("recurrence hreshold",~ epsilon)), labels = paste(round(RecScale$epsilon,2)))) +
        ggplot2::scale_fill_gradient2(low      = "red3",
                             high     = "steelblue",
                             mid      = "white",
                             na.value = scales::muted("slategray4"),
                             midpoint = mean(resol$value, na.rm = TRUE),
                             #limit    = c(min(meltRP$value, na.rm = TRUE),max(meltRP$value, na.rm = TRUE)),
                             space    = "Lab",
                             name     = "") +
        ggplot2::coord_equal(1, expand = FALSE) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.background = element_blank(),
              panel.grid.major  = element_blank(),
              panel.grid.minor  = element_blank(),
              legend.background = element_blank(),
              legend.key = element_blank(),
              panel.border = element_blank(),
              axis.text.y  =  element_text(size = 8),
              # axis.text.y.left  =  element_text(size = 8),
              # axis.text.y.right =  element_text(size = 8),
              axis.text.x  = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_text(hjust = 0, size = 10),
              # axis.title.y.left = element_text(hjust = 0, size = 10),
              # axis.title.y.right = element_text(hjust = 0, size = 10),
              plot.margin = margin(0,5,0,5, unit = "pt"))
    }

  } else { # unthresholded
  }

  rptheme <- ggplot2::theme_bw() + ggplot2::theme(panel.background = element_blank(),
                                panel.grid.minor  = element_blank(),
                                panel.border = element_rect("grey50",fill=NA),
                                legend.key = element_rect(colour = "grey90"),
                                axis.ticks = element_blank(),
                                axis.text = element_blank(),
                                plot.margin = margin(0,0,0,0))

  if(drawGrid){
    rptheme <- rptheme +  theme(panel.grid.major  = element_line("grey70",size = .1),
                                panel.ontop = TRUE)
  } else {
    rptheme <- rptheme +  theme(panel.grid.major  = element_blank(),
                                panel.ontop = FALSE)
  }

  if(showL){
    rptheme <- rptheme +  theme(legend.position = c(1.1,1),
                                legend.direction = "vertical",
                                legend.background = element_rect(fill="grey90"),
                                legend.justification = "top")
  }

  if(plotDimensions){
    rptheme <- rptheme + theme(axis.title.y =element_blank(),
                               axis.title.x =element_blank())
  }

  # Expand main plot ----
  if(!is.null(markEpochsLOI)){
    if(!unthresholded){
        gRP <- gRP +  ggplot2::scale_fill_manual(name  = "Key:",
                                    values = colvec,
                                    na.translate = TRUE ,
                                    na.value = scales::muted("slategray4"),
                                    guide = "legend",
                                    limits = levels(meltRP$value))
        }
          # theme(panel.ontop = TRUE,
          #       legend.position = "top",
          #       legend.background = element_rect(colour =  "grey50"))
    #geom_line(data = Edata, aes_(x=~x,y=~y,colour = ~value), size = 2, show.legend = TRUE) +

  } else {
    if(!unthresholded){
      gRP <- gRP +  ggplot2::scale_fill_manual(name  = "", breaks = c(0,1),
                                    values = colvec,
                                    na.translate = TRUE ,
                                    na.value = scales::muted("slategray4"),
                                    guide = "none")
      }
    }

  # if(!is.null(markEpochsGrid)){
  #   if(is.list(markEpochsGrid)&length(markEpochsGrid)==2){
  #
  #     markEpochsGrid[[1]] <- factor(markEpochsGrid[[1]])
  #     markEpochsGrid[[2]] <- factor(markEpochsGrid[[2]])
  #
  #     cpalV <- paletteer::paletteer_d(package = "rcartocolor",palette = "Safe", n = nlevels(markEpochsGrid[[1]]))
  #     names(cpalV) <- levels(markEpochsGrid[[1]])
  #
  #     dataV <- data.frame(t= as.numeric_character(markEpochsGrid[[1]]), xb=markEpochsGrid[[1]],col=NA)
  #     for(c in unique(dataV$t)){
  #       dataV$col[dataV$t%in%c] <- cpalV[c]
  #     }
  #
  #     cpalH <- paletteer::paletteer_d(package = "rcartocolor",palette = "Vivid", n = nlevels(markEpochsGrid[[2]]))
  #     names(cpalH) <- levels(markEpochsGrid[[2]])
  #
  #     dataH <- data.frame(t= as.numeric_character(markEpochsGrid[[2]]), yb=markEpochsGrid[[2]],col=NA)
  #     for(c in unique(dataV$t)){
  #       dataH$col[dataH$t%in%c] <- cpalH[c]
  #     }
#
#       gRP <- gRP +
#         #geom_vline(data = dataV, aes_(xintercept = diff(c((max(~t, na.rm = TRUE)+1),~xb)!=0)), colour = dataV$col) +
#         #geom_hline(data = dataH, aes_(yintercept = diff(c((max(~t, na.rm = TRUE)+1),~yb)!=0)), colour = dataH$col) +
#         ggplot2::geom_vline(data = dataV, aes(xintercept = diff(c((max(t)+1),t))!=0), colour = dataV$col) +
#         ggplot2::geom_hline(data = dataH, aes(yintercept = diff(c((max(t)+1),t))!=0), colour = dataH$col) +
#         ggplot2::theme(panel.ontop = TRUE)
#
#     } else {
#       warning("Variable passed to 'markEpochsGrid' is not a list, and/or is not of length 2.")
#     }
#   }

  if(plyr::is.discrete(meltRP$Var1)){
    gRP <- gRP + ggplot2::scale_x_discrete(breaks=meltRP$Var1,expand = c(0,0))
  } else {
    gRP <- gRP + ggplot2::scale_x_continuous(breaks=meltRP$Var1,expand = c(0,0))
  }
  if(plyr::is.discrete(meltRP$Var2)){
    gRP <- gRP + ggplot2::scale_y_discrete(breaks=meltRP$Var2,expand = c(0,0))
  } else {
    gRP <- gRP + ggplot2::scale_y_continuous(breaks=meltRP$Var2,expand = c(0,0))
  }
  gRP <- gRP + rptheme + ggplot2::coord_equal(dim(RM)[1]/dim(RM)[2]) #coord_fixed(expand = FALSE)

  gy1 <- gg_plotHolder()
  gy2 <- gg_plotHolder()

  xdims <- ""
  ydims <- ""

  if(!is.null(attr(RM,"emDims1.name"))){
    xdims <- ifelse(nchar(xlab)>0,xlab,attr(RM,"emDims1.name"))
    if(AUTO){
      ydims <- xdims
    } else {
      ydims <- ifelse(nchar(ylab)>0,ylab,attr(RM,"emDims2.name"))
    }
  }

  if(nchar(xlab)>0){
    xdims <- xlab
  }
  if(nchar(ylab)>0){
    ydims <- ylab
  }


  if(plotDimensions){

    if(!is.null(attr(RM,"emDims1"))){

      gRP <- gRP + ylab("") + xlab("")

      y1 <- data.frame(t1=attr(RM,"emDims1"))
      y2 <- data.frame(t2=attr(RM,"emDims2"))

      # Y1

      colnames(y1) <- paste0("X",1:NCOL(y1))
      y1$tm  <- 1:NROW(y1)
      y1$tmna <- 0
      y1$tmna[is.na(y1[,1])] <- y1$tm[is.na(y1[,1])]
      y1 <- tidyr::gather(y1,key="Dimension",value = "Value", -c("tm","tmna"))
      y1$Value <-  elascer(y1$Value)

      # Y2
      colnames(y2) <- paste0("Y",1:NCOL(y2))
      y2$tm <- 1:NROW(y2)
      y2$tmna <- 0
      y2$tmna[is.na(y2[,1])] <- y2$tm[is.na(y2[,1])]
      y2 <- tidyr::gather(y2,key="Dimension",value = "Value", -c("tm","tmna"))
      y2$Value <-  elascer(y2$Value)


    } else {
      if(unthresholded){
          y1 <- data.frame(Value=diag(RM),x=1:length(diag(RM)), Dimension = rep("LOS",length(diag(RM))))
          y2 <- data.frame(Value=diag(RM),x=1:length(diag(RM)), Dimension = rep("LOS",length(diag(RM))))
          xdims <- paste("LOS",xdims)
          ydims <- paste("LOS",ydims)
      } else {
        plotDimensions <- FALSE
      }
      }
  } else {
    gRP <- gRP + ylab(ydims) + xlab(xdims)
  }

  if(plotDimensions){

      gy1 <- ggplot2::ggplot(y1, ggplot2::aes_(y=~Value, x= ~tm,  group= ~Dimension)) +
        ggplot2::geom_line(aes_(colour=~Dimension), show.legend = FALSE) +
        ggplot2::xlab(xdims) + ggplot2::ylab(" ") +
        ggplot2::geom_vline(ggplot2::aes_(xintercept = ~tmna), colour = scales::muted("slategray4"),alpha=.1, size=.5) +
        ggplot2::scale_color_grey() +
        ggplot2::scale_x_continuous(expand = c(0,0)) +
        ggplot2::scale_y_continuous(expand = c(0,0)) +
        ggplot2::theme(panel.background = element_blank(),
              panel.grid.major  = element_blank(),
              panel.grid.minor  = element_blank(),
              legend.background = element_blank(),
              legend.key = element_blank(),
              panel.border = element_blank(),
              axis.text = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.title.x =element_text(colour = "black",angle = 0, vjust = +3),
              axis.title.y =element_blank(),
              plot.margin = margin(0,0,0,0, unit = "pt")) +
        ggplot2::coord_cartesian(expand = FALSE)  # +  coord_fixed(1/10)

      gy2 <- ggplot2::ggplot(y2, ggplot2::aes_(y=~Value, x=~tm, group=~Dimension)) +
        ggplot2::geom_line(ggplot2::aes_(colour=~Dimension), show.legend = FALSE) +
        ggplot2::ylab(" ") + ggplot2::xlab(ydims) +
        ggplot2::geom_vline(ggplot2::aes_(xintercept = ~tmna), colour = scales::muted("slategray4"),alpha=.1, size=.5) +
        ggplot2::scale_color_grey() +
        ggplot2::scale_x_continuous(expand = c(0,0)) +
        ggplot2::theme(panel.background = element_blank(),
              panel.grid.major  = element_blank(),
              panel.grid.minor  = element_blank(),
              legend.background = element_blank(),
              legend.key = element_blank(),
              panel.border = element_blank(),
              axis.text = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.title.x =element_blank(),
              axis.title.y =element_text(colour = "black",angle = 90, vjust = -2),
              plot.margin = margin(0,0,0,0, unit = "pt")) +
        ggplot2::coord_flip(expand = FALSE) +
        ggplot2::scale_y_reverse(expand = c(0,0))

  } # plotdimensions


  if(plotMeasures){

    rpOUT    <- round(rpOUT,3)
    if(is.na(rpOUT$emRad)){
      rpOUT$Radius <- round(radiusValue,3)
    } else {
      rpOUT$Radius <- rpOUT$emRad
      }


    rpOUTdat <- rpOUT %>%
      dplyr::select(dplyr::one_of(c("Radius","RP_N","RR","DET","MEAN_dl","ENT_dl","LAM_vl","TT_vl","ENT_vl"))) %>%
      tidyr::gather(key="measure",value="value") %>%
      dplyr::mutate(x=rep(0,9),y=9:1)

    #rpOUTdat <- cbind(rpOUTdat,rpOUTdat)
    rpOUTdat$label <-  paste0(rpOUTdat$measure,":\n",rpOUTdat$value)

    gA <-ggplot2::ggplot(rpOUTdat,ggplot2::aes_(x=~x,y=~y)) +
      ggplot2::geom_text(ggplot2::aes_(label=~label), family="mono", hjust="left", vjust="center", size=3, parse = FALSE) +
      #scale_x_continuous(limits = c(0,.3)) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = margin(0,5,0,5, unit = "pt"))

    #geom="text", label = paste("Radius:",rpOUT$Radius,"\nRec points:",rpOUT$RT,"\nRR",rpOUT$RR,"\nDET:",rpOUT$DET,"\nMEAN_dl:",rpOUT$MEAN_dl,"\nENT_dl:",rpOUT$ENT_dl,"\nLAM_vl:",rpOUT$LAM_vl, "\nTT_vl:",rpOUT$TT_vl,"\nENTR_vl:",rpOUT$ENT_vl)) + theme_minimal() + theme(text = element_text(family = "mono"))
    # ,"\nLAM_hl:",rpOUT$LAM_vl, "| TT_hl:",rpOUT$TT_vl,"| ENTR_hl:",rpOUT$ENT_hl))
  }

  if(useGtable){

    #gRP <- gRP #+ theme(panel.background = element_rect(colour="white"))

    g <- ggplot2::ggplotGrob(gRP)

    if(plotDimensions){
      gry2<-ggplot2::ggplotGrob(gy2)
      gry1<-ggplot2::ggplotGrob(gy1)
    }

    if(unthresholded&plotRadiusRRbar){
      grDist <- ggplot2::ggplotGrob(gDist)
    }

    if(plotMeasures){
      grA <- ggplot2::ggplotGrob(gA)
    }


    if(plotDimensions&!plotMeasures&unthresholded&plotRadiusRRbar){
      mat <- matrix(list(gry2, grid::nullGrob(),g, gry1, grDist, grid::nullGrob()),nrow = 2)
      gt  <- gtable::gtable_matrix("di_rp_dim", mat, widths = unit(c(.25, 1,.5), "null"), heights =  unit(c(1,.25), "null"),respect = TRUE)
    }

    if(plotDimensions&!plotMeasures&unthresholded&!plotRadiusRRbar){
      mat <- matrix(list(gry2, grid::nullGrob(),g, gry1),nrow = 2)
      gt  <- gtable::gtable_matrix("di_rp_dim", mat, widths = unit(c(.25, 1), "null"), heights =  unit(c(1,.25), "null"),respect = TRUE)
    }

    if(plotDimensions&!plotMeasures&!unthresholded){
      mat <- matrix(list(gry2, grid::nullGrob(),g, gry1),nrow = 2)
      gt  <- gtable::gtable_matrix("bi_rp_dim", mat, widths = unit(c(.25, 1), "null"), heights =  unit(c(1, .25), "null"),respect = TRUE)
    }

    if(plotDimensions&plotMeasures&unthresholded&plotRadiusRRbar){
      mat <- matrix(list(grA, grid::nullGrob(), gry2, grid::nullGrob(),g, gry1, grDist, grid::nullGrob()),nrow = 2)
      gt  <- gtable::gtable_matrix("di_rp_dim_meas", mat, widths = unit(c(.35,.25, 1,.5), "null"), heights =  unit(c(1,.25), "null"),respect = TRUE)
    }

    if(plotDimensions&plotMeasures&unthresholded&!plotRadiusRRbar){
      mat <- matrix(list(grA, grid::nullGrob(), gry2, grid::nullGrob(),g, gry1),nrow = 2)
      gt  <- gtable::gtable_matrix("di_rp_dim_meas", mat, widths = unit(c(.35,.25, 1), "null"), heights =  unit(c(1,.25), "null"),respect = TRUE)
    }

    if(plotDimensions&plotMeasures&!unthresholded){
      mat <- matrix(list(grA, grid::nullGrob(), gry2, grid::nullGrob(),g, gry1),nrow = 2)
      gt  <- gtable::gtable_matrix("bi_rp_dim_meas", mat, widths = unit(c(.35,.25, 1), "null"), heights =  unit(c(1,.25), "null"),respect = TRUE)
    }

    if(!plotDimensions&plotMeasures&unthresholded&plotRadiusRRbar){
      mat <- matrix(list(grA, g, grDist),nrow = 1)
      gt  <- gtable::gtable_matrix("di_rp_meas", mat, widths = unit(c(.35, 1,.5), "null"), heights =  unit(c(1), "null"),respect = TRUE)
    }

    if(!plotDimensions&plotMeasures&unthresholded&!plotRadiusRRbar){
      mat <- matrix(list(grA, g),nrow = 1)
      gt  <- gtable::gtable_matrix("di_rp_meas", mat, widths = unit(c(.35, 1), "null"), heights =  unit(c(1), "null"),respect = TRUE)
    }

    if(!plotDimensions&plotMeasures&!unthresholded){
      mat <- matrix(list(grA, g),nrow = 1)
      gt  <- gtable::gtable_matrix("bi_rp_meas", mat, widths = unit(c(.35, 1), "null"), heights =  unit(c(1), "null"),respect = TRUE)
    }

    if(!plotDimensions&!plotMeasures&unthresholded&plotRadiusRRbar){
      mat <- matrix(list(g, grDist),nrow = 1)
      gt  <- gtable::gtable_matrix("di_rp", mat, widths = unit(c(1,.5), "null"), heights =  unit(c(1), "null"),respect = TRUE)
    }

    if(!plotDimensions&!plotMeasures&unthresholded&!plotRadiusRRbar){
      mat <- matrix(list(g),nrow = 1)
      gt  <- gtable::gtable_matrix("di_rp", mat, widths = unit(c(1), "null"), heights =  unit(c(1), "null"),respect = TRUE)
    }

    if(!plotDimensions&!plotMeasures&!unthresholded){
      mat <- matrix(list(g),nrow = 1)
      gt  <- gtable::gtable_matrix("bi_rp", mat, widths = unit(c(1), "null"), heights =  unit(c(1), "null"),respect = TRUE)
    }

    # gindex <- subset(g$layout, name == "layout")
    # g <- gtable::gtable_add_cols(g, grid::unit(.5, "grobwidth", data = g),0)
    # g <- gtable::gtable_add_grob(g, ggplot2::ggplotGrob(gA), t=gindex$t, l=1, b=gindex$b, r=gindex$l)
    #rm(gindex)
    #}

    if(nchar(title)>0){
      grT <- ggplot2::ggplot(data.frame(x=1,y=1)) +
        ggplot2::geom_text(ggplot2::aes_(x=~x,y=~y), label=title) +
        theme_void() +
        theme(plot.margin = margin(0,0,0,0, unit = "pt"))
      gt  <- gtable::gtable_add_rows(x = gt, heights =  unit(c(.1),"null"), pos=0)
      l <- sum(c(plotDimensions,plotMeasures))+1
      gt  <- gtable::gtable_add_grob(x = gt, grobs = ggplot2::ggplotGrob(grT), name = "Title", t=1, l=l)
    }

  #  gt <- gtable::gtable_add_col_space(gt,)

    g <- gtable::gtable_add_padding(gt, unit(5, "pt"))

  }


  # else {
  #
  #   if(plotDimensions){
  #
  #     if(unthresholded){
  #       g <- (gy2 + gRP + gDist + gg_plotHolder() + gy1 + gg_plotHolder() +
  #               patchwork::plot_layout(nrow = 2, ncol = 3, widths = c(1,10,1), heights = c(10,1)) + patchwork::plot_annotation(title = title, caption = ifelse(AUTO,"Auto-recurrence plot","Cross-recurrence plot")))
  #
  #     } else {
  #
  #       if(plotMeasures){
  #         g <- (gy2 + gRP + gA + gg_plotHolder() + gy1 + gg_plotHolder() +
  #                 patchwork::plot_layout(nrow = 2, ncol = 3, widths = c(1,9,2), heights = c(10,1)) + patchwork::plot_annotation(title = title, caption = ifelse(AUTO,"Auto-recurrence plot","Cross-recurrence plot")))
  #       } else {
  #         g <- (gy2 + gRP + gg_plotHolder() + gg_plotHolder() + gy1 + gg_plotHolder() +
  #                 patchwork::plot_layout(nrow = 2, ncol = 3, widths = c(1,9,2), heights = c(10,1)) + patchwork::plot_annotation(title = title, caption = ifelse(AUTO,"Auto-recurrence plot","Cross-recurrence plot")))
  #       }
  #     }
  #
  #   } else {
  #     g <- gRP
  #   }
  # } # use gtable

  if(!returnOnlyObject){
    if(useGtable){
      grid::grid.newpage()
      grid::grid.draw(g)
    } else {
      # graphics::plot.new()
      graphics::plot(g)
    }
  }
  return(invisible(g))
}


#' rp_size
#'
#' @param mat A Matrix object
#' @param AUTO Is the Matrix an Auto Recurrence Matrix? If so, the length of the diagonal will be subtracted from the matrix size, pass `FALSE` to prevent this behaviour. If `NULL` (default) `AUTO` will take on the value of `isSymmetric(mat)`.
#' @param theiler Should a Theiler window be applied?
#'
#' @return Matrix size for computation of recurrence measures.
#' @export
#'
#' @family Distance matrix operations (recurrence plot)
#'
#' @examples
#' # Create a 10 by 10 matrix
#' library(Matrix)
#' m <- Matrix(rnorm(10),10,10)
#'
#' rp_size(m,TRUE,0)   # Subtract diagonal
#' rp_size(m,FALSE,0)  # Do not subtract diagonal
#' rp_size(m,NULL,0)   # Matrix is symmetrical, AUTO is set to TRUE
#' rp_size(m,NULL,1)   # Subtract a Theiler window of 1 around and including the diagonal
rp_size <- function(mat, AUTO=NULL, theiler = NULL){
  if(is.null(AUTO)){
    AUTO <- Matrix::isSymmetric(Matrix::unname(mat))
  }
  if(is.null(theiler)){
    if(!is.null(attributes(mat)$theiler)){
      theiler <- attr(mat,"theiler")
    } else {
      theiler <- 0
    }
  }
  return(cumprod(dim(mat))[2] - ifelse((AUTO&theiler==0),length(Matrix::diag(mat)),
                                       ifelse(theiler>0,Matrix::nnzero(Matrix::band(mat,-theiler,theiler)),0)))
}


#' Empty results vector
#'
#' @return an empty crqa_rp
#' @keywords internal
#' @export
#'
crqa_rp_empty <- function(){
  data.frame(
    emRad   = NA,
    RP_N     = NA,
    RR       = NA,
    DET      = NA,
    MEAN_dl  = NA,
    MAX_dl   = NA,
    ENT_dl   = NA,
    ENTrel_dl= NA,
    REP_av  = NA,
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

#' crqa_rp_calc
#'
#' @param RM RM
#' @param emRad r
#' @param DLmin d
#' @param VLmin v
#' @param HLmin h
#' @param DLmax dd
#' @param VLmax vv
#' @param HLmax hh
#' @param theiler t
#' @param AUTO a
#' @param chromatic c
#' @param matrices m
#'
#' @return crqa measures
#' @export
#' @keywords internal
#'
crqa_rp_calc <- function(RM,
                     emRad = NULL,
                     DLmin = 2,
                     VLmin = 2,
                     HLmin = 2,
                     DLmax = length(Matrix::diag(RM))-1,
                     VLmax = length(Matrix::diag(RM))-1,
                     HLmax = length(Matrix::diag(RM))-1,
                     theiler = 0,
                     AUTO      = NULL,
                     chromatic = FALSE,
                     matrices  = FALSE){

  recmatsize <- rp_size(RM, AUTO=AUTO)

    if(!is.null(attributes(RM)$emRad)){
      emRad <- attributes(RM)$emRad
    } else {
    emRad <- NA
  }


  if(!all(as.vector(RM)==0|as.vector(RM)==1)){
    if(!chromatic){
          stop("Need a thresholded recurrence matrix")
    } else {
      stop("Chromatic not yet implemented")
        }
    }

  #Total nr. recurrent points
  RP_N <- Matrix::nnzero(RM, na.counted = FALSE)

  minDiag <- 0
  if(Matrix::isSymmetric(Matrix::unname(RM))){
    if(all(Matrix::diag(RM)==1)){
      minDiag <- length(Matrix::diag(RM))
      }
  }

  RP_N <- RP_N-minDiag

  #Proportion recurrence / Recurrence Rate
  RR <- RP_N/recmatsize

  if(length(RR)==0){RR<-0}

  if(RR==1){
    warning("Everything is recurring!\nReturning empty vector")
    return(crqa_rp_empty())
  }


  #Get line segments
  lineSegments <- rp_lineDist(RM,
                              DLmin = DLmin, DLmax = DLmax,
                              VLmin = VLmin, VLmax = VLmax,
                              HLmin = HLmin, HLmax = HLmax,
                              theiler = theiler, AUTO = AUTO)

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

  # Number of lines
  N_dl <- sum(freq_dl, na.rm = TRUE)
  N_vl <- sum(freq_vl, na.rm = TRUE)
  N_hl <- sum(freq_hl, na.rm = TRUE)

  #Number of recurrent points on diagonal, vertical and horizontal lines
  N_dlp <- sum(freqvec_dl*freq_dl, na.rm = TRUE)
  N_vlp <- sum(freqvec_vl*freq_vl, na.rm = TRUE)
  N_hlp <- sum(freqvec_hl*freq_hl, na.rm = TRUE)

  #Determinism / Horizontal and Vertical Laminarity
  DET    <- N_dlp/RP_N
  LAM_vl <- N_vlp/RP_N
  LAM_hl <- N_hlp/RP_N

  #anisotropy ratio
  #ANI    <- ((N_vlp-N_hlp)/N_dlp)%00%NA
  ANI    <- (N_vlp/N_hlp)%00%NA

  # Singularities
  SING <- rp_lineDist(RM,
                      DLmin = 1, DLmax = DLmax,
                      VLmin = 1, VLmax = VLmax,
                      HLmin = 1, HLmax = HLmax,
                      theiler = theiler, AUTO = AUTO)

  # table(SING_N$verticals.dist)
  # table(SING_N$horizontals.dist)
  SING_N  <-  table(SING$diagonals.dist)[1]
  SING_rate <- SING_N / RP_N

  #Array of probabilities that a certain line length will occur (all >1)
  P_dl <- freq_dl/N_dl
  P_vl <- freq_vl/N_vl
  P_hl <- freq_hl/N_hl

  #Entropy of line length distributions
  ENT_dl <- -1 * sum(P_dl * log(P_dl))
  ENT_vl <- -1 * sum(P_vl * log(P_vl))
  ENT_hl <- -1 * sum(P_hl * log(P_hl))

  #Relative Entropy (Entropy / Max entropy)
  ENTrel_dl = ENT_dl/(-1 * log(1/DLmax))
  ENTrel_vl = ENT_vl/(-1 * log(1/VLmax))
  ENTrel_hl = ENT_hl/(-1 * log(1/HLmax))

  #Meanline
  MEAN_dl = mean(dlines, na.rm = TRUE)
  MEAN_vl = mean(vlines, na.rm = TRUE)
  MEAN_hl = mean(hlines, na.rm = TRUE)

  #Maxline
  MAX_dl = max(freqvec_dl, na.rm = TRUE)
  MAX_vl = max(freqvec_vl, na.rm = TRUE)
  MAX_hl = max(freqvec_hl, na.rm = TRUE)

  # REPetetiveness
  REP_av <- ((N_hlp/N_dlp) + (N_vlp/N_dlp))/2
  REP_hl  <-  N_hlp/N_dlp
  REP_vl  <-  N_vlp/N_dlp

  #Coefficient of determination
  CoV_dl = stats::sd(dlines)/mean(dlines)
  CoV_vl = stats::sd(vlines)/mean(vlines)
  CoV_hl = stats::sd(hlines)/mean(hlines)

  #Divergence
  DIV_dl = 1/MAX_dl
  DIV_vl = 1/MAX_vl
  DIV_hl = 1/MAX_hl

  #Output
  out <- data.frame(
    emRad     = emRad,
    RP_N      = RP_N,
    RR        = RR,
    SING_N    = SING_N,
    SING_rate = SING_rate,
    DIV_dl    = DIV_dl,
    REP_av    = REP_av,
    ANI       = ANI,
    N_dl      = N_dl,
    N_dlp     = N_dlp,
    DET       = DET,
    MEAN_dl   = MEAN_dl,
    MAX_dl    = MAX_dl,
    ENT_dl    = ENT_dl,
    ENTrel_dl = ENTrel_dl,
    CoV_dl    = CoV_dl,
    N_vl      = N_vl,
    N_vlp     = N_vlp,
    LAM_vl    = LAM_vl,
    TT_vl     = MEAN_vl,
    MAX_vl    = MAX_vl,
    ENT_vl    = ENT_vl,
    ENTrel_vl = ENTrel_vl,
    CoV_vl    = CoV_vl,
    REP_vl    = REP_vl,
    #    DIV_vl   = DIV_vl,
    N_hlp     = N_hlp,
    N_hl      = N_hl,
    LAM_hl    = LAM_hl,
    TT_hl     = MEAN_hl,
    MAX_hl    = MAX_hl,
    ENT_hl    = ENT_hl,
    ENTrel_hl = ENTrel_hl,
    CoV_hl    = CoV_hl,
    REP_hl    = REP_hl)
  #    DIV_hl   = DIV_hl)

  if(matrices){
    return(list(
      crqaMeasures = out,
      crqaMatrices = list(RM     = RM,
                          dlines = dlines,
                          vlines = vlines,
                          hlines = hlines,
                          freq_dl = freq_dl,
                          freq_vl = freq_vl,
                          freq_hl = freq_hl))
      )
  } else {
    return(out)
  }
}


#' Prepare matrix
#'
#' @param RP Recurrence plot
#' @param emRad Radiuc
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
#' @keywords internal
#'
#' @export
#'
crqa_rp_prep <- function(RP,
                     emRad = NULL,
                     DLmin = 2,
                     VLmin = 2,
                     HLmin = 2,
                     DLmax = length(Matrix::diag(RP)),
                     VLmax = length(Matrix::diag(RP)),
                     HLmax = length(Matrix::diag(RP)),
                     AUTO      = FALSE,
                     chromatic = FALSE,
                     matrices  = FALSE,
                     doHalf    = FALSE){

  out<-crqa_rp_calc(RP,
                emRad = emRad,
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
      outLo <- crqa_rp_calc(Matrix::tril(RP,-1),
                        emRad = emRad,
                        DLmin = DLmin,
                        VLmin = VLmin,
                        HLmin = HLmin,
                        DLmax = DLmax,
                        VLmax = VLmax,
                        HLmax = HLmax,
                        AUTO  = AUTO,
                        chromatic = chromatic,
                        matrices  = matrices)

      outUp <- crqa_rp_calc(Matrix::triu(RP,-1),
                        emRad= emRad,
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
      warning("doHalf = TRUE and AUTO = TRUE. Results would be the same for lower and upper triangle!")
      out<- cbind.data.frame(full  = out,
                             lower = crqa_rp_empty(),
                             upper = crqa_rp_empty())
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
#     if(!is.null(attr(RM,"emDims1"))){
#       plotVarsOK <- TRUE
#       plotTimeOK <- TRUE
#       plotVars   <- list(attributes(RM)$emDims1,attributes(RM)$emDims2)
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
#            distPallette <- grDevices::colorRampPalette(c("white", "black"))(2),
#            distPallette <- grDevices::colorRampPalette(colors = c("red3", "snow", "steelblue"),
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
#     grp <- ggplot(data = meltRP, aes_(x = meltRP[,2],
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


# Networks ----

#' Create a Recurrence Network Matrix
#'
#' This function serves as a wrapper for function `rp()`, it will add some attributes to the matrix related to network representation. These attributes will be used to decide which network type to generate (e.g. undirected, directed, weighted, etc.)
#'
#' @inheritParams rp
#' @param directed Should the matrix be considered to represent a directed network? (default = `FALSE`)
#' @param weighted Should the matrix be considered to represent a weighted network? (default = `FALSE`)
#' @param includeDiagonal Should the diagonal of the matrix be included when creating the network (default = `FALSE`)
#' @param returnGraph Return an [igraph::igraph()] object (default = `TRUE`)
#' @param ... Any paramters to pass to [rn_plot()] if `doPlot = TRUE`
#'
#' @return A (Coss-) Recurrence matrix that can be interpreted as an adjacency (or incidence) matrix.
#'
#' @export
#'
#' @family Distance matrix operations (recurrence network)
#'
#'
rn <- function(y1, y2 = NULL,
               emDim = 1,
               emLag = 1,
               emRad = NULL,
               directed = FALSE,
               weighted = FALSE,
               includeDiagonal = FALSE,
               to.ts = NULL,
               order.by = NULL,
               to.sparse = FALSE,
               method = "Euclidean",
               targetValue = .05,
               returnGraph = TRUE,
               doPlot = FALSE,
               silent = TRUE,
               ...){

  dmat <- rp(y1 = y1, y2 = y2,
           emDim = emDim, emLag = emLag, emRad = emRad,
           to.ts = to.ts, order.by = order.by, to.sparse = to.sparse,
           weighted = weighted,targetValue = targetValue,
           method = method, doPlot = FALSE, silent = silent)

  if(to.sparse){
    attributes(dmat)$directed <- directed
    attributes(dmat)$includeDiagonal <- includeDiagonal
  } else {
    attr(dmat,"directed") <- directed
    attr(dmat,"includeDiagonal") <- includeDiagonal
  }


  if(doPlot){
    dotArgs <- list(...)

    if(is.null(dotArgs)){
      dotArgs<- formals(rn_plot)
      }

    nameOk  <- names(dotArgs)%in%methods::formalArgs(rn_plot)
    # Plot with defaults
    if(!all(nameOk)){
      dotArgs <- formals(rn_plot)
      nameOk  <- rep(TRUE,length(dotArgs))
    }

    dotArgs$RN <- dmat

    do.call(rn_plot, dotArgs[nameOk])
  }

  return(dmat)
}


#' Plot (thresholded) distance matrix as a network
#'
#' @param RN A distance matrix or recurrence matrix
#' @inheritParams rp_plot
#'
#' @return A nice plot of the recurrence network
#' @export
#'
#' @family Distance matrix operations (recurrence network)
#'
rn_plot <- function(RN,
                    plotDimensions = FALSE,
                    plotMeasures   = FALSE,
                    drawGrid = FALSE,
                    markEpochsLOI = NULL,
                    Chromatic = NULL,
                    radiusValue = NA,
                    title = "", xlab = "", ylab="",
                    plotSurrogate = NA,
                    returnOnlyObject = FALSE){

  rp_plot(RM = RN,
          plotDimensions  = plotDimensions,
          plotMeasures    = plotMeasures,
          plotRadiusRRbar = FALSE,
          drawGrid = drawGrid,
          markEpochsLOI = markEpochsLOI,
          Chromatic = Chromatic,
          radiusValue = radiusValue,
          title = title, xlab = xlab, ylab = ylab,
          plotSurrogate = plotSurrogate,
          returnOnlyObject = returnOnlyObject)

}


#' Recurrence Time Spectrum
#'
#' Get the recurrence time distribution from a recurrence network.
#'
#' @param RN A thresholded recurrence matrix generated by function `rn()`
#' @param fitRange If `NULL` the entire range will be used for log-log slope. If a 2-element vector of integers, this will represent the range of recurrence times to use for fitting the log=log slope (e.g. `c(1,50)` would fit the first 50 recurrence times).
#' @param fs Sample rate (default = `1`)
#' @param doPlot Should a plot of the recurrence time spectrum be produced?
#' @param returnPlot Return ggplot2 object (default = `FALSE`)
#' @param returnPLAW Return the power law data (default = `FALSE`)
#' @param returnInfo Return all the data used in SDA (default = `FALSE`)
#' @param silent Silent-ish mode
#' @param noTitle Do not generate a title (only the subtitle)
#' @param tsName Name of y added as a subtitle to the plot
#'
#' @return A vector of frequencies of recurrence times and a plot (if requested)
#'
#' @export
#'
#' @family Distance matrix operations (recurrence network)
#'
rn_recSpec <- function(RN,
                       fitRange = NULL,
                       fs = 1,
                       doPlot = TRUE,
                       returnPlot = FALSE,
                       returnPLAW = FALSE,
                       returnInfo = FALSE,
                       silent = TRUE,
                       noTitle = FALSE,
                       tsName="y"){

  if(is.null(attributes(RN)$weighted)){
    stop("Wrong RN format: Create a thresholded recurrence matrix using function rn()")
  }

  diagonal <- attributes(RN)$includeDiagonal
  weighted <- NULL
  if(attributes(RN)$weighted){weighted <- TRUE}

  g1 <- igraph::graph_from_adjacency_matrix(RN, mode="upper", diag = diagonal, weighted = weighted)

  edgeFrame <- igraph::as_data_frame(g1,"edges")
  edgeFrame$rectime <- edgeFrame$to-edgeFrame$from

  #d <- density(edgeFrame$rectime, kernel = "cosine")
  #plot(d)
  #d <- table(edgeFrame$rectime)

  RT <- edgeFrame$rectime
  f1 = seq(0,1, by = 1/(1000*fs))
  P1 = graphics::hist(x = 1/(RT+1), breaks = f1)

  # f = seq(1,2*max(RT))
  # P = graphics::hist(RT, f)
  # f2 = fs/(P$mids+.5)


  yl <- c(smooth(P1$counts)+1)
  xl <- c(fs*P1$mids)
  #lreg <- stats::lm(log2(yl[yl>0] ~ xl[yl>0])
  #pred <- predict(object = lreg)

  # d <- graphics::hist(edgeFrame$rectime,breaks = (0:NROW(RN))+0.5, plot = FALSE)$count
  # names(d) <- paste(1:NROW(RN))
  # d <- d[d>0]

  ddata <- data.frame(size = xl, bulk = as.numeric(yl))
  ddata <- dplyr::arrange(ddata,size)
  rownames(ddata) <- paste(1:NROW(ddata))

  if(is.null(fitRange)){
    nr <- 1:round(NROW(ddata)/4)
  } else {
    if(length(fitRange)==2&is.numeric(fitRange)){
      if(fitRange[1]<0){fitRange[1] <- 1}
      if(fitRange[2]>NROW(ddata)){fitRange[2] <- NROW(ddata)}
      nr <- fitRange[1]:fitRange[2]
    } else {
      stop("Wrong fitRange format: Either NULL or a 2-integer vector.")
    }
  }

  nr <- which(nr%in%rownames(ddata))

  lmfit1  <- stats::lm(log(ddata$bulk) ~ log(ddata$size))
  alpha1 <- stats::coef(lmfit1)[2]

  lmfit2  <- stats::lm(log(ddata$bulk[nr]) ~ log(ddata$size[nr]))
  alpha2 <- stats::coef(lmfit2)[2]

  outList <- list(
    PLAW  =  ddata,
    fullRange = list(sap = alpha1,
                     # H = 1+stats::coef(lmfit1)[2] + Hadj,
                     # FD = sa2fd_sda(stats::coef(lmfit1)[2]),
                     fitlm1 = lmfit1,
                     method = paste0("Full range (n = ",length(ddata$size),")\nSlope = ",round(stats::coef(lmfit1)[2],2))), #" | FD = ",round(sa2fd_sda(stats::coef(lmfit1)[2]),2))),
    fitRange  = list(sap = stats::coef(lmfit2)[2],
                     # H = 1+stats::coef(lmfit2)[2] + Hadj,
                     # FD = sa2fd_sda(stats::coef(lmfit2)[2]),
                     fitlm2 = lmfit2,
                     method = paste0("Fit range (n = ",length(ddata$size[nr]),")\nSlope = ",round(stats::coef(lmfit2)[2],2))), # | FD = ",round(sa2fd_sda(stats::coef(lmfit2)[2]),2))),
    info = list(edgelist=edgeFrame,fitdata=ddata),
    plot = NA,
    analysis = list(
    name = "Recurrence Time Spectrum",
    logBaseFit = "log2",
    logBasePlot = 2))

  if(doPlot|returnPlot){
    if(noTitle){
      title <- ""
    } else {
      title <- "log-log regression (Recurrence Time Spectrum)"
    }
    g <- plotFD_loglog(fd.OUT = outList, title = title, subtitle = tsName, logBase = "2", xlabel = "Frequency", ylabel = "Power")
    if(doPlot){
      grid::grid.newpage()
      grid::grid.draw(g)
    }
    if(returnPlot){
      outList$plot <- invisible(g)
    }
  }

  if(returnInfo){returnPLAW<-TRUE}

  return(outList[c(returnPLAW,TRUE,TRUE,returnInfo,returnPlot)])


#if(doPlot){
#
#   breaks_x <- scales::log_breaks(n=abs(diff(range(round(log2(ddata$x)))+c(-1,1))),base=2)(ddata$x)
#   labels_x <- eval(quote(scales::trans_format("log2", scales::math_format(2^.x,format=scales::number_format(accuracy = .1)))(breaks_x)))
#
#   breaks_y <- scales::log_breaks(n=abs(diff(range(round(log2(ddata$y)))+c(-1,1))),base=2)(ddata$y)
#   labels_y <- eval(quote(scales::trans_format("log2", scales::math_format(2^.x,format=scales::number_format(accuracy = .1)))(breaks_y)))
#
#   g <- ggplot2::ggplot(ddata, ggplot2::aes_(x=~x_l,y=~y_l)) +
#     scale_x_continuous(breaks = log2(breaks_x),
#                        labels = labels_x,  expand = c(0,0)) + # limits = lims_x) +
#     scale_y_continuous(breaks = log2(breaks_y),
#                        labels = labels_y, expand = c(0,0)) + # limits = lims_y) +
#     geom_point() +
#     geom_line(ggplot2::aes_(x=~x_p,y=~y_p), colour = "red3") +
#     ylab("Recurrence Times (log2)")+xlab("Frequency of Occurrence (log2)") +
#     #geom_label(aes(x=ddata$x_l[25],y=ddata$y_l[25],label=round(alpha,2))) +
#     annotation_logticks() +
#     theme_bw()

#   grid::grid.newpage()
#   grid::grid.draw(g)
# }
#
# if(returnOnlyObject){
#   return(invisible(g))
# }
#
#   if(returnPowerLaw){
#     return(ddata)
#   }

}


#' Recurrence Time Scaleogram
#'
#' Display a recurrence network in a space representing Time (x-axis) x Scale (y-axis). The scale axis will be determined by the latency between the occurence of a value in the (embedded) time series vector and its recurrences in the future (i.e. only the upper triangle of the recurrence matrix will be displayed, excluding the diagonal).
#'
#' @param RN A thresholded recurrence matrix generated by function `rn()`
#' @param returnOnlyObject Return the `ggplot` / `ggraph` object only, do not draw the plot (default = `FALSE`)
#'
#' @return A `ggraph` graph object
#' @export
#'
#' @family Distance matrix operations (recurrence network)
#'
rn_scaleoGram <- function(RN, returnOnlyObject = FALSE){

  if(is.null(attributes(RN)$emRad)){
    stop("Wrong RN format: Create a thresholded recurrence matrix using function rn()")
  }

  g1 <- igraph::graph_from_adjacency_matrix(RN, mode="upper", diag = FALSE, weighted = NULL)

  edgeFrame <- igraph::as_data_frame(g1,"edges")
  edgeFrame$rectime <- edgeFrame$to-edgeFrame$from
  edgeFrame$rectime_bin <- cut(x = edgeFrame$rectime, ordered_result = TRUE, breaks = round(seq.int(1, max(edgeFrame$rectime), length.out = 11)), dig.lab = 3, include.lowest = TRUE, labels = FALSE)

  #edgeFrame$rectime_bin <- as.numeric_character(edgeFrame$rectime_bin)
  # table(edgeFrame$rectime_bin)
  # range(edgeFrame$rectime[edgeFrame$rectime_bin==1])
  #
  # edgeFrame$vcol <- "#000000"
  # edgeFrame <- arrange(edgeFrame,edgeFrame$from,edgeFrame$rectime)
  #edgeFrame$scale.vertex <- interaction(edgeFrame$rectime,edgeFrame$from)

  timepoints <- data.frame(name = as.numeric(igraph::V(g1)), degree =  igraph::degree(g1))
  vertexFrame <- dplyr::left_join(x = expand.grid(name=timepoints$name, rt_scale = c(0,sort(unique(edgeFrame$rectime)))), y = timepoints, by = "name")

  #Change Names
  vertexFrame$name <- paste0(vertexFrame$name,".",vertexFrame$rt_scale)
  edgeFrame$from <- paste0(edgeFrame$from,".0")
  edgeFrame$to   <- paste0(edgeFrame$to,".",edgeFrame$rectime)

  colvec <- paletteer::paletteer_c(package = "scico",palette = "vik",n = length(unique(edgeFrame$rectime)))
  names(colvec) <- sort(unique(edgeFrame$rectime))
  edgeFrame$colour <- NA
  for(c in names(colvec)){
    edgeFrame$colour[edgeFrame$rectime%in%c] <- colvec[names(colvec)%in%c]
  }

  g2 <- igraph::graph_from_data_frame(edgeFrame, directed = FALSE, vertices = vertexFrame)
  igraph::E(g2)$color <- edgeFrame$colour

  scaleogram <- ggraph::create_layout(g2, layout = "linear")
  scaleogram$x <- as.numeric(plyr::laply(as.character(scaleogram$name),function(s) strsplit(s,"[.]")[[1]][1]))
  scaleogram$y <- as.numeric(plyr::laply(as.character(scaleogram$name),function(s) strsplit(s,"[.]")[[1]][2]))
  #scaleogram$rt_scale <- factor(scaleogram$y)
  #attributes(scaleogram)

  y1 <- data.frame(t1=attr(RN,"emDims1"))

  # Y1
  colnames(y1) <- paste0("X",1:NCOL(y1))
  y1$tm  <- 1:NROW(y1)
  y1$tmna <- 0
  y1$tmna[is.na(y1[,1])] <- y1$tm[is.na(y1[,1])]
  y1 <- tidyr::gather(y1,key="Dimension",value = "Value", -c("tm","tmna"))
  y1$Value <-  elascer(y1$Value,lo = -round(max(igraph::E(g2)$rectime)/8),hi = 0)

  g <- ggraph::ggraph(scaleogram) +
    #geom_node_point(colour = "grey60", size=.1, alpha = .3) +
    ggraph::geom_edge_diagonal(aes_(edge_colour = ~rectime, group = ~rectime), edge_alpha=.1,lineend = "round", linejoin = "round") +
    ggplot2::geom_line(data = y1, ggplot2::aes_(y=~Value, x= ~tm, colour = ~Dimension, group = ~Dimension), show.legend = FALSE, size = .1) +
    ggplot2::scale_color_grey() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
    ggraph::scale_edge_color_gradient2("Recurrence Time", low = paste(colvec[1]), mid = paste(colvec[round(length(colvec)/2)]), high = paste(colvec[length(colvec)]), midpoint = round(length(colvec)/2)) +
    ggplot2::scale_x_continuous("Time", expand = c(0,0),
                       breaks = round(seq.int(1,max(igraph::V(g1)),length.out = 10)),
                       limits = c(0,max(as.numeric(igraph::V(g1))))) +
    ggplot2::scale_y_continuous("Recurrence Time", expand = c(0,0),
                       breaks = c(-round(max(igraph::E(g2)$rectime)/8),sort(unique(igraph::E(g2)$rectime))[round(seq.int(1,length(unique(igraph::E(g2)$rectime)), length.out = 10))]),
                       labels = c("", paste(sort(unique(igraph::E(g2)$rectime))[round(seq.int(1,length(unique(igraph::E(g2)$rectime)), length.out = 10))])),
                       limits = c(-round(max(igraph::E(g2)$rectime)/8),max(igraph::E(g2)$rectime))) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank())

  if(!returnOnlyObject){
      #grid::grid.newpage()
      grid::grid.draw(g)
  }
  return(invisible(g))
}


#' Mutual Information Function
#'
#'  Calculate the lagged mutual information fucntion within (auto-mif) or between (cross-mif) time series, or, conditional on another time series (conditional-cross-mif). Alternatively, calculate the total information of a multivariate dataset for different lags.
#'
#' @param y A `Nx1` matrix for auto-mif, a `Nx2` matrix or data frame for cross-mif, a `Nx3` matrix or data frame for mif between col 1 and 2 conditional on col 3; or a `NxM` matrix or data frame for the multi-information function. Mutual information for each lag will be calculated using functions in package [infotheo::infotheo()] for `lags` lagged versions of the time series.
#' @param lags The lags to evaluate mutual information.
#' @param nbins The number of bins passed to [infotheo::discretize()] if y is a matrix or [casnet::ts_discrete()]
#' @param doPlot Produce a plot of the symbolic time series by calling [plotRED_mif()] (default = `FALSE`)
#' @param surTest If `TRUE`, a surrogate will be conducted using simple surrogates. The surrogates will be created from the transition probabilities of the discretised time series, i.e. the probability of observing bin `j` when the current value is in bin `j`. The number of surrogates needed will be computed based on the value of the `alpha` parameter, conceived as a one-sided test: `mi > 0`.
#' @param alpha The alpha level for the surrogate test (default = `0.05`)
#'
#' @return The auto- or cross-mi function
#' @export
#'
#' @family Redundancy measures (mutual information)
#'
#' @examples
#'
#' # Lags to evaluate mututal information
#' lags <- -10:30
#'
#' # Auto-mutual information
#' y1 <- sin(seq(0, 100, by = 1/8)*pi)
#'
#' (mif(data.frame(y1),lags = lags))
#'
#' # Cross-mututal information, y2 is a lagged version y1
#' y2 <- y1[10:801]
#'
#' y <- data.frame(ts_trimfill(y1, y2, action = "trim.cut"))
#' (mif(y,lags = lags))
#'
#' # Conditional mutual information, add some noise to y2 and add it as a 3rd column
#' y$s <- y2+rnorm(NROW(y2))
#' (mif(y,lags = lags))
#'
#' # Multi-information, the information of the entire multivariate series at each lag
#' y$y3 <- cumsum(rnorm(NROW(y)))
#' (mif(y,lags = lags))
#'
#'
mif <- function(y, lags=-10:10, nbins = ceiling(2*NROW(y)^(1/3)), doPlot = FALSE, surTest = FALSE, alpha = 0.05){

  if(is.null(dim(y))){
    y <- as.matrix(y,ncol=1)
  }

  cnt <- 0
  N <- NROW(y)
  mif_out <- numeric(length=length(lags))
  for(i in lags) {
    cnt <- cnt + 1
    # shift y2 to left for neg lags
    if(i < 0) {

      ID2 <- -i:N
      ID1 <- 1:(N-abs(i)+1)
    }
    if(i == 0) {
      ID2 <- c(1:N)
      ID1 <- c(1:N)
    }
    # shift y2 to right for pos lags
    if(i > 0) {
      ID2 <- 1:(N-i+1)
      ID1 <- i:N
    }
    mif_out[cnt] <- mi_mat(y, ID1, ID2, discreteBins = nbins)
  }

  names(mif_out) <- paste(lags)
  if(NCOL(y)==1){miType <- "I(X;X)"}
  if(NCOL(y)==2){miType <- "I(X;Y)"}
  if(NCOL(y)==3){miType <- "I(X;Y|Z)"}
  if(NCOL(y)> 3){miType <- "I(X;Y;Z;...;N)"}
  attr(mif_out,"miType") <- miType

  if(doPlot){
    plotRED_mif(mif.OUT = mif_out, lags = lags, nbins = nbins)
  }

  return(mif_out)
}


#' Mutual Information variations
#'
#' @param y A matrix with time series in columns
#' @param ID1 ids
#' @param ID2 ids
#' @param discreteBins Number of bins to use when discretizing the time series
#'
#' @return mi in nats
#' @export
#'
#' @family Redundancy measures (mutual information)
#'
mi_mat <- function(y, ID1, ID2, discreteBins = ceiling(2*NROW(ID1)^(1/3))){
  Nc <- NCOL(y)
  if(!is.null(dim(y))){
    if(Nc == 1){out <- infotheo::mutinformation(X = infotheo::discretize(y[ID1,1], nbins = discreteBins),
                                                Y = infotheo::discretize(y[ID2,1], nbins = discreteBins))}
    if(Nc == 2){out <- infotheo::mutinformation(X = infotheo::discretize(y[ID1,1], nbins = discreteBins), #length(seq(-.5,(NROW(ID1)-.5)))),
                                                Y = infotheo::discretize(y[ID2,2], nbins = discreteBins))} #length(seq(-.5,(NROW(ID2)-.5)))))}
    if(Nc == 3){out <- infotheo::condinformation(X = infotheo::discretize(y[ID1,1], nbins = discreteBins),
                                                 Y = infotheo::discretize(y[ID2,2], nbins = discreteBins),
                                                 S = infotheo::discretize(y[ID1,3], nbins = discreteBins))}
    if(Nc >  3){out <- infotheo::multiinformation(X = infotheo::discretize(y[ID1,], nbins = discreteBins))}
  return(out)
  } else {
    warning("Input should be a matrix or data frame.")
  }
}


mi_ts <- function(y1,y2=NULL, nbins=NA){

  if(is.null(y2)){y2 <- y1}

    # y1 <-  ts_checkfix(y1,fixNumericVector = TRUE)
    # y2 <-  ts_checkfix(y2,fixNumericVector = TRUE)

    equal <- data.frame(ts_trimfill(y1,y2,action = "trim.cut"))
    idNA1 <- is.na(equal[,1])
    idNA2 <- is.na(equal[,2])

    equal <- equal[!apply(apply(equal,1,is.na),2,any),]

    if(is.na(nbins)){ nbins <- ceiling(2*NROW(equal)^(1/3))}

    equal[,1]  <- ts_discrete(equal[,1], nbins = nbins)
    equal[,2]  <- ts_discrete(equal[,2], nbins = nbins)

    ## init entropies
    H_s <- H_u <- H_su <- 0

     ## get ts length
     TT <- NROW(equal)
     bb <- max(equal)
     ## get ts length

     for(i in 1:bb) {
       ## calculate entropy for 1st series
       p_s <- sum(equal[,1] == i)/TT
       if(p_s != 0) { H_s <- H_s - p_s*log(p_s, 2) }
       ## calculate entropy for 2nd series
       p_u <- sum(equal[,2] == i)/TT
       if(p_u != 0) { H_u <- H_u - p_u*log(p_u, 2) }
       for(j in 1:bb) {
         ## calculate joint entropy
         p_su <- sum(equal[,1]==i & equal[,2]==j)/TT
         if(p_su != 0) { H_su <- H_su - p_su*log(p_su, 2) }
       } ## end j loop
     } ## end i loop
     ## calc mutual info
    return(MI <- H_s + H_u - H_su)
   #  if(!normal) { return(MI) } else { return(MI/sqrt(H_s*H_u)) }
}


#' Inter-layer mutual information
#'
#' @param g0 An igraph object representing a layer in a multiplex graph
#' @param g1 An igraph object representing a layer in a multiplex graph
#' @param probTable Option to return the table with marginal and joint degree distribution probabilities (default = `TRUE`)
#'
#' @return The inter-layer mutual information between `g1` and `g2`. If `probTable=TRUE`, a list object with two fields, the inter-layer mutual information and the table with marginal and joint degree distributions
#' @export
#'
#' @family Redundancy measures (mutual information)
#' @family Multiplex Recurrence Networks
#'
mif_interlayer <- function(g0,g1, probTable=FALSE){
  d0    <- igraph::degree_distribution(g0)
  d1    <- igraph::degree_distribution(g1)

  equal <- data.frame(ts_trimfill(d0,d1))

  p01 <- graphics::hist(x = igraph::degree(g0),breaks = seq(-.5,(NROW(equal)-.5)),plot=FALSE)$counts
  p10 <- graphics::hist(x = igraph::degree(g1),breaks = seq(-.5,(NROW(equal)-.5)),plot=FALSE)$counts

  equal$joint <-  (p01+p10) / sum(p01,p10,na.rm = TRUE)
  equal$degree <- seq_along(equal[,1])-1

  # imiAB <- sum(equal$joint * log(equal$joint/(equal[,1]*equal[,2]+.Machine$double.eps))%00%0, na.rm = TRUE)
  imiAB <- infotheo::mutinformation(p01,p10)
  if(probTable){
    attributes(imiAB) <- list(miType = "inter-layer mutual information", probTable = equal)
  } else {
    attributes(imiAB) <- list(miType = "inter-layer mutual information")
  }
  return(imiAB)
}


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
  if(grepl("Matrix",class(distmat))){
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
#' @return A matrix with 0s and leaves the values < threshold distance value
#'
#' @export
#'
#' @family Distance matrix operations (recurrence plot)
#' @family Distance matrix operations (recurrence network)
#'
di2we <- function(distmat, emRad, convMat = FALSE){

  matPack <- FALSE
  if(grepl("Matrix",class(distmat))){
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



#' Layout a graph on a spiral
#'
#' @param g An igraph object. If (`rev = FALSE`) the vertex with the lowest index will be placed in the centre of the spiral, the highest index will be most outer vertex,
#' @param type Spiral type, one of `"Archimedean"`,`"Bernoulli"`,`"Fermat"`, or, `"Euler"` (default = `"Archimedean"`)
#' @param arcs The number of arcs (half circles/ovals) that make up the spiral (default = `10`)
#' @param a Parameter controlling the distance between spiral arms, however, the effect will vary for different spiral types (default = `0.5`)
#' @param b Parameter controlling where the spiral originates. A value of 1 will generally place the origin in the center. The default value places the origin right of the center, however, the effect will vary for different spiral types (default = `0.1`)
#' @param rev If `TRUE` the vertex with the highest index will be placed in the centre of the spiral (default = `FALSE`)
#'
#' @return An igraph layout
#'
#' @export
#'
#' @examples
#'
#' g  <- sample_gnp(100, 1/100)
#'
#' # Equiangular spiral: Any line from the origin cuts at the same angle.
#' plot(g, layout = layout_as_spiral(g, type = "Bernoulli", arcs = 5))
#'
#' # The arms of Fermat's spiral diverge quadratically.
#' plot(g, layout = layout_as_spiral(g, type = "Fermat", arcs = 5))
#'
#' # Equidistance of intersection points along a line through the origin.
#' plot(g, layout = layout_as_spiral(g, type = "Archimedean", arcs = 5))
#'
layout_as_spiral <- function(g,
                             type = c("Archimedean","Bernoulli","Fermat","Euler"),
                             arcs = 6,
                             a = 1,
                             b = NULL,
                             rev= FALSE){

  N <- igraph::vcount(g)
  if(length(unique(V(g)$size))>1){
    res <- round(max(V(g)$size)*N)
  } else {
    res <- N*2
  }

  theta <- seq(0,arcs*pi,length.out = res)

  if(type=="Archimedean"){
    if(is.null(b)){
      b <- 1
    }
    r  <- a + b*theta
    df <- matrix(cbind(r*cos(theta),r*sin(theta)),ncol = 2)
  }
  if(type=="Bernoulli"){
    if(is.null(b)){
      b <- 0.1
    }
    r  <- a + exp(b*theta)
    df <- matrix(cbind(r*cos(theta),r*sin(theta)),ncol = 2)
  }
  if(type=="Fermat"){
    if(is.null(b)){
      b <- 1
    }
    r  <- a + b*(theta*theta)
    df <- matrix(cbind(r*cos(theta),r*sin(theta)),ncol = 2)
  }
  if(type=="Euler"){
    res1  <- ceiling(res/2)
    res2  <- floor(res/2)
    theta <- seq(0,pi/2,length.out = res1)
    df1   <- data.frame(t=theta,x=NA,y=NA)
    df1$x[1] <-0
    df1$y[1] <-0

    dt <-arcs/res1

    for(i in 2:res1){
      df1$x[i] <- df1$x[i-1] + cos(df1$t[i-1]*df1$t[i-1]) * dt
      df1$y[i] <- df1$y[i-1] + sin(df1$t[i-1]*df1$t[i-1]) * dt
      df1$t[i] <- df1$t[i-1]+dt
    }
    df <- matrix(cbind(c(-rev(df1$x[1:res2]),df1$x),c(-rev(df1$y[1:res2]),df1$y)),ncol = 2)
  }

  df     <- df[seq(1, res, length.out = N),]
  df[,1] <- elascer(df[,1],lo = 0.1, hi = 0.9)
  df[,2] <- elascer(df[,2],lo = 0.1, hi = 0.9)

  if(rev){
    df[,1] <- rev(df[,1])
    df[,2] <- rev(df[,2])
  }
  return(df)
}

#' Make Spiral Graph
#'
#' Turn an [igraph] object into a spiral graph returning a [ggplot2] object.
#'
#' @note To keep the igraph object, use the layout function [layout_as_spiral(g)] when plotting the graph.
#'
#' @inheritParams layout_as_spiral
#' @param markTimeBy Include a vector that indicates time. The time will be displayed on the plot. Pass `TRUE` to generate auto labels (experimental)
#' @param curvature The `curvature` parameter for edges see [geom_curve()] (default = `-0.7`)
#' @param angle The `angle` parameter for edges see [geom_curve()] (default = `90`)
#' @param title A title for the plot
#' @param subtitle A subtitle for the plot
#' @param showEpochLegend Should a legend be shown for the epoch colours? (default = `TRUE`)
#' @param markEpochsBy A vector of length `vcount(g)` indicating epochs or groups (default = `NULL`)
#' @param epochColours A vector of length `vcount(g)` with colour codes (default = `NULL`)
#' @param epochLabel A title for the epoch legend (default = `"Epoch"`)
#' @param showSizeLegend Should a legend be shown for the size of the nodes? (default = `FALSE`)
#' @param sizeLabel Use to indicate if `V(g)$size` represents some measure, e.g. [igraph::degree()], or, [igraph::hubscore()] (default = `"Size"`)
#' @param scaleVertexSize Scale the size of the vertices by setting a range for [ggplot2::scale_size()]. This will not affect the numbers on the size legend (default = `c(1,6)`)
#' @param vertexBorderColour Draw a border around the vertices. Pass `NULL` to use the same colour as the fill colour (default = `"black"`)
#' @param scaleEdgeSize Scale the size of the edges by a constant: `E(g)$width * scaleEdgeSize` (default = `1/5`)
#' @param defaultEdgeColour Colour of edges that do not connect to the same epoch (default = `"grey70"`)
#' @param doPlot Produce a plot? (default = `TRUE`)
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#'
#' g  <- sample_gnp(200, 1/20)
#' V(g)$size <- degree(g)
#' make_spiral_graph(g, markTimeBy = TRUE, showSizeLegend = TRUE, sizeLabel = "Node degree")
#'
make_spiral_graph <- function(g,
                              type = "Archimedean",
                              arcs = 6,
                              a = .5,
                              b = .1,
                              rev= FALSE,
                              curvature = -0.6,
                              angle = 90,
                              markTimeBy = NULL,
                              alphaV = 1,
                              alphaE = .6,
                              title = "",
                              subtitle = "",
                              showEpochLegend = TRUE,
                              markEpochsBy = NULL,
                              epochColours = NULL,
                              epochLabel = "Epoch",
                              showSizeLegend = FALSE,
                              sizeLabel = "Size",
                              scaleVertexSize = c(1,6),
                              vertexBorderColour = "black",
                              scaleEdgeSize = 1/5,
                              defaultEdgeColour = "grey70",
                              doPlot = TRUE){

  g$layout <- layout_as_spiral(g, type = type, arcs = arcs, a = a, b = b, rev = rev)

  if(is.null(markTimeBy)){
    if(!is.null(markEpochsBy)){
      grIDs <- ts_changeindex(markEpochsBy)
      tbreaks <- unique(sort(c(grIDs$xmax,grIDs$xmax)))
      tlabels <- ceiling(tbreaks)
    } else {
      x <- 1:vcount(g)
      v <- seq(1,vcount(g),by=vcount(g)/arcs)
      tbreaks <- c(1,which(diff(findInterval(x, v))!=0),vcount(g))
      #tbreaks <- which(diff(c(g$layout[1,1],g$layout[,2],g$layout[1,2])>=g$layout[1,2])!=0)
      if(max(tbreaks)!=vcount(g)){
        tbreaks[which.max(tbreaks)]<-vcount(g)
        tlabels <- paste(tbreaks)
      }
      if(min(tbreaks)>1){
        tbreaks<- c(1,tbreaks)
        tlabels <- paste(tbreaks)
      }
    }
  } else {
    if(is.numeric(markTimeBy)){
      if(all(markTimeBy%in%1:vcount(g))){
        tbreaks <- unique(markTimeBy)
        if(!is.null(names(markTimeBy))){
          tlabels <- names(unique(markTimeBy))
        } else {
          tlabels <- paste(tbreaks)
        }
      }
    }
    if(markTimeBy){
      tbreaks <- which(diff(c(g$layout[1,1],g$layout[,2],g$layout[1,2])>=g$layout[1,2])!=0)
      if(max(tbreaks)>vcount(g)){tbreaks[which.max(tbreaks)]<-vcount(g)}
      if(min(tbreaks)>1){tbreaks<- c(1,tbreaks)}
      tlabels <- paste(tbreaks)
    }
  }

  if(max(tbreaks)!=vcount(g)){
    tbreaks[which.max(tbreaks)]<-vcount(g)
    tlabels <- paste(tbreaks)
    }
  if(min(tbreaks)>1){
    tbreaks<- c(1,tbreaks)
    tlabels <- paste(tbreaks)
    }

  if(is.null(markEpochsBy)){
    markEpochsBy <- character(vcount(g))
    for(i in 1:(length(tbreaks)-1)){
    markEpochsBy[tbreaks[i]:tbreaks[i+1]] <- rep(paste0(tbreaks[i],"-",tbreaks[i+1]),length(tbreaks[i]:tbreaks[i+1]))
    }
    if(!is.null(epochColours)){
      if(length(unique(markEpochsBy))>length(unique(epochColours))){
        warning("Number of unique epochs is unequal to number of unique colours!\nUsing default colour scheme.")
        epochColours <- NULL
      }
    }
  }

  g <- plotNET_groupColour(g,
                           groups = markEpochsBy,
                           colourV = TRUE,
                           colourE = TRUE,
                           groupColours = epochColours,
                           defaultEdgeColour = defaultEdgeColour,
                           doPlot = FALSE)

  size <- 1
  if(!is.null(V(g)$size)){
    size <- V(g)$size
  }
  gNodes        <- as.data.frame(g$layout)
  gNodes$ID     <- as.numeric(V(g))
  gNodes$colour <- V(g)$colour
  gNodes$labels <- factor(V(g)$groupnum, levels = unique(V(g)$groupnum), labels = unique(V(g)$group))
  gNodes$size   <- size
  gNodes$alpha  <- V(g)$alpha

  width <- 1
  if(!is.null(E(g)$width)){
    width <- E(g)$width
  }
  gEdges        <- igraph::get.data.frame(g)
  gEdges$from.x <- gNodes$V1[match(gEdges$from, gNodes$ID)]
  gEdges$from.y <- gNodes$V2[match(gEdges$from, gNodes$ID)]
  gEdges$to.x   <- gNodes$V1[match(gEdges$to, gNodes$ID)]
  gEdges$to.y   <- gNodes$V2[match(gEdges$to, gNodes$ID)]
  gEdges$width  <- width

  if(is.null(vertexBorderColour))(
    vBc <- gNodes$colour
  ) else {
    if(length(vertexBorderColour)==1|length(vertexBorderColour)==NROW(gNodes)){
    vBc <- vertexBorderColour
    } else {
      warning("Invalif value(s) for vertexBorderColour, using default.")
      vertexBorderColour <- "black"
    }
  }

  gg <- ggplot(gNodes,aes(x=V1,y=V2)) +
    geom_curve(data=gEdges, aes(x = from.x, xend = to.x, y = from.y, yend = to.y),
               curvature = curvature,
               angle = angle,
               size = gEdges$width * scaleEdgeSize,
               colour= gEdges$color,
               alpha = alphaE) +
    geom_point(aes(fill = labels, size = size), pch=21, colour = vBc, alpha = alphaV) +
    ggtitle(label = title, subtitle = subtitle) +
    scale_fill_manual(epochLabel, values = unique(gNodes$colour)) +
    scale_size(sizeLabel, range = scaleVertexSize)

  if(showEpochLegend){
    gg <- gg + guides(fill = guide_legend(title.position = "top",
                                          byrow = TRUE,
                                          override.aes = list(size=5, order = 0)))
  } else {
    gg <- gg + guides(fill = "none")
  }

  if(showSizeLegend){
    gg <- gg + guides(size = guide_legend(title.position = "top",
                                          byrow = TRUE,
                                          override.aes = list(legend.key.size = unit(1.2,"lines"), order = 1)))
  } else {
    gg <- gg + guides(size = "none")
  }

    if(!is.null(markTimeBy)){
      gg <- gg + annotate("label", x=gNodes$V1[tbreaks], y=gNodes$V2[tbreaks], label = tlabels)
    }

  gg <- gg +
    coord_fixed() +
    theme_void() +
    theme(legend.title = element_text(face="bold"),
          legend.position =  "top",
          legend.margin = margin(t = 0,r = 1,l = 1,0))

    if(doPlot){
      print(gg)
    }

  return(invisible(gg))
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


# (M)FD estimators ----------------------------------------------


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
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. <http://doi.org/10.3389/fphys.2013.00075>
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
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. <http://doi.org/10.3389/fphys.2013.00075>
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
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. <http://doi.org/10.3389/fphys.2013.00075>
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
  RelR   <- 2*(1-VAR$acf[2] / VAR$acf[1])
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
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. <http://doi.org/10.3389/fphys.2013.00075>
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
#' @details Calls function [sapa::SDF()] to estimate the scaling exponent of a timeseries based on the periodogram frequency spectrum. After detrending and normalizing the signal (if requested), `SDF` is called using a Tukey window (\code{raised cosine \link[sapa]{taper}}).
#'
#' A line is fitted on the periodogram in log-log coordinates. The full ramge is fitted as well as one of three fit-ranges:
#' \itemize{
#' \item{`lowest25` - The 25\% lowest frequencies}
#' \item{`Wijnants` - The 50 lowest frequencies (Wijnants et al., 2012)}
#' \item{`HurvichDeo` - The Hurvich-Deo estimate, see ([fractal::HDEst()])}
#' }
#'
fd_psd <- function(y,
                   fs = NULL,
                   standardise = TRUE,
                   detrend = TRUE,
                   fitMethod = c("lowest25","Wijnants","Hurvich-Deo")[3],
                   doPlot = FALSE,
                   returnPlot = FALSE,
                   returnPLAW = FALSE,
                   returnInfo = FALSE,
                   silent = FALSE,
                   noTitle = FALSE,
                   tsName="y"){

  if(!stats::is.ts(y)){
    if(is.null(fs)){fs <- 1}
    y <- stats::ts(y, frequency = fs)
   if(!silent){cat("\n\nfd_psd:\tSample rate was set to 1.\n\n")}
  }

  N <- length(y)
  # Simple linear detrending.
  if(detrend){
    y <- stats::ts(ts_detrend(as.vector(y)), frequency = fs)
    }
  # standardise using N instead of N-1.
  if(standardise){
    y <- stats::ts(ts_standardise(y,adjustN=FALSE), frequency = fs)
    }

  # Number of frequencies estimated cannot be set! (defaults to Nyquist)
  # Use Tukey window: cosine taper with r = 0.5

  # fast = TRUE ensures padding with zeros to optimize FFT to highly composite number.
  # However, we just pad to nextPow2, except if length already is a power of 2.
  npad <- 1+(stats::nextn(N,factors=2)-N)/N
  npad <- stats::nextn(N)

  # if(N==npad) npad = 0
  # psd  <- stats::spec.pgram(y, fast = FALSE, demean=FALSE, detrend=FALSE, plot=FALSE, pad=npad, taper=0.5)

  Tukey <- sapa::taper(type="raised cosine", flatness = 0.5, n.sample = N)
  psd   <- sapa::SDF(y, taper. = Tukey, npad = npad)

  powspec <- cbind.data.frame(freq.norm = attr(psd, "frequency")[-1], size = attr(psd, "frequency")[-1]*stats::frequency(y), bulk = as.matrix(psd)[-1])

  # First check the global slope for anti-persistent noise (GT +0.20)
  # If so, fit the line starting from the highest frequency
  nr     <- length(powspec[,1])
  lsfit  <- stats::lm(log(powspec$bulk[1:nr]) ~ log(powspec$size[1:nr]))
  glob   <- stats::coef(lsfit)[2]

  # General guideline: fit over 25% frequencies
  # If signal is continuous (sampled) consider Wijnants et al. (2013) log-log fitting procedure
  nr <- switch(fitMethod,
    "lowest25" = which.min(powspec$size>=0.25),
    "Hurvich-Deo" = fractal::HDEst(NFT = length(powspec$bulk), sdf = as.vector(powspec$bulk)),
    "Wijnants" = 50
  )

  if(nr>=length(powspec$freq.norm)){nr <- length(powspec$freq.norm)-1}
  exp1 <- fractal::hurstSpec(y, sdf.method="direct", freq.max = powspec$freq.norm[length(powspec$freq.norm)-1], taper.=Tukey )
  exp2 <- fractal::hurstSpec(y, sdf.method="direct", freq.max = powspec$freq.norm[nr], taper.=Tukey)

  ifelse((glob > 0.2), {
    lmfit1 <- stats::lm(log(rev(powspec$bulk)) ~ log(rev(powspec$size)))
    lmfit2 <- stats::lm(log(rev(powspec$bulk[1:nr])) ~ log(rev(powspec$size[1:nr])))
  },{
    lmfit1 <- stats::lm(log(powspec$bulk) ~ log(powspec$size))
    lmfit2 <- stats::lm(log(powspec$bulk[1:nr]) ~ log(powspec$size[1:nr]))
  })

  outList <- list(
    PLAW  = powspec,
    fullRange = list(sap = stats::coef(lmfit1)[2],
                     H = exp1,
                     FD = sa2fd_psd(stats::coef(lmfit1)[2]),
                     fitlm1 = lmfit1,
                     method = paste0("All frequencies (n = ",length(powspec$freq.norm),")\nSlope = ",round(stats::coef(lmfit1)[2],2)," | FD = ",sa2fd_psd(stats::coef(lmfit1)[2]))),
    fitRange  = list(sap = stats::coef(lmfit2)[2],
                     H = exp2,
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
        title <- ""
      } else {
        title <- "log-log regression (PSD)"
      }
    if(doPlot){

      g <- plotFD_loglog(fd.OUT = outList, title = title, subtitle = tsName, logBase = "10",xlabel = "Normalised Frequency", ylabel = "Power")

      if(returnPlot){
        outList$plot <- g
      }
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


#' fd_sda
#'
#' @title Standardised Dispersion Analysis (SDA).
#'
#' @param y    A numeric vector or time series object.
#' @param fs Sample rate (default = NULL)
#' @param standardise standardise the series (default = "mean.sd")
#' @param detrend Subtract linear trend from the series (default = FALSE)
#' @param polyOrder Order of detrending polynomial
#' @param adjustSumOrder  Adjust the time series (summation or differencing), based on the global scaling exponent, see e.g. <https://www.frontiersin.org/files/Articles/23948/fphys-03-00141-r2/image_m/fphys-03-00141-t001.jpg>{Ihlen (2012)} (default = `FALSE`)
#' @param scaleMax   Maximum scale to use
#' @param scaleMin   Minimium scale to use
#' @param scaleResolution  The scales at which the standardised fluctuations are calculated as: `(scaleMax-scaleMin)/scaleResolution`
#' @param scaleS If not `NA`, it should be a numeric vector listing the scales on which to evaluate the fluctuations. Arguments `scaleMax, scaleMin, scaleResolution` will be ignored.
#' @param overlap Turn SDA into a sliding window analysis. A number in `[0 ... 1]` representing the amount of 'bin overlap'. If `length(y) = 1024` and overlap is `.5`, a scale of `4` will be considered a sliding window of size `4` with stepsize `floor(.5 * 4) = 2` (default = `0`)
#' @param minData Minimum number of data points in a bin needed to calculate standardised dispersion
#' @param doPlot   Output the log-log scale versus fluctuation plot with linear fit by calling function `plotFD_loglog()` (default = `TRUE`)
#' @param returnPlot Return ggplot2 object (default = `FALSE`)
#' @param returnPLAW Return the power law data (default = `FALSE`)
#' @param returnInfo Return all the data used in SDA (default = `FALSE`)
#' @param silent Silent-ish mode
#' @param noTitle Do not generate a title (only the subtitle)
#' @param tsName Name of y added as a subtitle to the plot
#'
#' @author Fred Hasselman
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. <http://doi.org/10.3389/fphys.2013.00075>
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
fd_sda <- function(y,
                   fs = NULL,
                   standardise = c("mean.sd","median.mad")[1],
                   detrend = FALSE,
                   polyOrder=1,
                   adjustSumOrder = FALSE,
                   scaleMin = 2,
                   scaleMax = floor(log2(NROW(y)/2)),
                   scaleResolution = 30,
                   scaleS = NA,
                   overlap = 0,
                   minData = 4,
                   doPlot = FALSE,
                   returnPlot = FALSE,
                   returnPLAW = FALSE,
                   returnInfo = FALSE,
                   silent = FALSE,
                   noTitle = FALSE,
                   tsName="y"){


  if(!stats::is.ts(y)){
    if(is.null(fs)){fs <- 1}
    y <- stats::ts(y, frequency = fs)
    if(!silent){cat("\n\nfd_sda:\tSample rate was set to 1.\n\n")}
  }

  N             <- length(y)
  # Simple linear detrending.
  if(detrend){y <- ts_detrend(y,polyOrder = polyOrder)} # y <- stats::ts(pracma::detrend(as.vector(y), tt = 'linear'), frequency = fs)
  # standardise using N instead of N-1.
  if(is.na(scaleS)){
    scaleS <- unique(round(2^(seq(scaleMin, scaleMax, by=((scaleMax-scaleMin)/scaleResolution)))))
  }

  if(max(scaleS)>(NROW(y)/2)){
    scaleS <- scaleS[scaleS<=(NROW(y)/2)]
  }

  if(!all(is.numeric(scaleS),length(scaleS)>0,scaleS%[]%c(2,(NROW(y)/2)))){
    message("Something wrong with vector passed to scaleS.... \nUsing default: (scaleMax-scaleMin)/scaleResolution")
  }

  # Standardise by N
  if(any(standardise%in%c("mean.sd","median.mad"))){
    y <- ts_standardise(y, type = standardise,  adjustN = FALSE)
  }

  if(adjustSumOrder){
    y       <- ts_sumorder(y, scaleS = scaleS, polyOrder = polyOrder, minData = minData)
    Hadj    <- attr(y,"Hadj")
    Hglobal <- attr(y,"Hglobal.excl")
  } else {
    Hadj    <- 0
    Hglobal <- NA
  }

  out <- fractal::dispersion(y, front = FALSE)

  fitRange <- which(lengths(lapply(out$scale, function(s){ts_slice(y,s)}))>=minData)

  lmfit1        <- stats::lm(log(out$sd) ~ log(out$scale))
  lmfit2        <- stats::lm(log(out$sd[fitRange]) ~ log(out$scale[fitRange]))

  outList <- list(
    PLAW  =  cbind.data.frame(freq.norm = stats::frequency(y)/out$scale,
                              size = out$scale,
                              bulk = out$sd),
    fullRange = list(sap = stats::coef(lmfit1)[2],
                     H = 1+stats::coef(lmfit1)[2] + Hadj,
                     FD = sa2fd_sda(stats::coef(lmfit1)[2]),
                     fitlm1 = lmfit1,
                     method = paste0("Full range (n = ",length(out$scale),")\nSlope = ",round(stats::coef(lmfit1)[2],2)," | FD = ",round(sa2fd_sda(stats::coef(lmfit1)[2]),2))),
    fitRange  = list(sap = stats::coef(lmfit2)[2],
                     H = 1+stats::coef(lmfit2)[2] + Hadj,
                     FD = sa2fd_sda(stats::coef(lmfit2)[2]),
                     fitlm2 = lmfit2,
                     method = paste0("Fit range (n = ",length(out$scale[fitRange]),")\nSlope = ",round(stats::coef(lmfit2)[2],2)," | FD = ",round(sa2fd_sda(stats::coef(lmfit2)[2]),2))),
    info = out,
    plot = NA,
    analysis = list(
      name = "Standardised Dispersion Analysis",
      logBaseFit = "log",
      logBasePlot = "e")
    )

  if(doPlot|returnPlot){
    if(noTitle){
      title <- ""
    } else {
      title <- "log-log regression (SDA)"
    }

    if(doPlot){
    g <- plotFD_loglog(fd.OUT = outList, title = title, subtitle = tsName, logBase = "e", ylabel = "Standardised Dispersion")
    if(returnPlot){
      outList$plot <- g
    }
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
#' @param removeTrend Method to use for detrending, see [fractal::DFA()] (default = "poly")
#' @param polyOrder Order of polynomial trend to remove if `removeTrend = "poly"`
#' @param standardise Standardise by the series using [casnet::ts_standardise()] with `adjustN = FALSE` (default = "mean.sd")
#' @param adjustSumOrder  Adjust the time series (summation or differencing), based on the global scaling exponent, see e.g. <https://www.frontiersin.org/files/Articles/23948/fphys-03-00141-r2/image_m/fphys-03-00141-t001.jpg>{Ihlen (2012)} (default = `FALSE`)
#' @param scaleMax   Maximum scale (as a power of 2) to use
#' @param scaleMin   Minimium scale (as a power of 2) to use
#' @param scaleResolution  The scales at which detrended fluctuation will be evaluated are calculatd as: `(scaleMax-scaleMin)/scaleResolution`. The default value yields no resolution of scales: `(scaleMax-scaleMin)`. Common values
#' @param scaleS If not `NA`, it should be a numeric vector listing the scales on which to evaluate the detrended fluctuations. Arguments `scaleMax, scaleMin, scaleResolution` will be ignored.
#' @param overlap Turn DFA into a sliding window analysis. A number in `[0 ... 1]` representing the amount of 'bin overlap'. If `length(y) = 1024` and overlap is `.5`, a scale of `4` will be considered a sliding window of size `4` with stepsize `floor(.5 * 4) = 2`. The detrended fluctuation in   For scale `128` this will be  (default = `0`)
#' @param minData Minimum number of data points in a bin needed to calculate detrended fluctuation
#' @param doPlot   Return the log-log scale versus fluctuation plot with linear fit (default = `TRUE`).
#' @param returnPlot Return ggplot2 object (default = `FALSE`)
#' @param returnPLAW Return the power law data (default = `FALSE`)
#' @param returnInfo Return all the data used in DFA (default = `FALSE`)
#' @param silent Silent-ish mode
#' @param noTitle Do not generate a title (only the subtitle)
#' @param tsName Name of y added as a subtitle to the plot
#'
#'
#' @return Estimate of Hurst exponent (slope of `log(bin)` vs. `log(RMSE))` and an FD estimate based on Hasselman(2013)
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
#' @references Hasselman, F. (2013). When the blind curve is finite: dimension estimation and model inference based on empirical waveforms. Frontiers in Physiology, 4, 75. <http://doi.org/10.3389/fphys.2013.00075>
#'
#' @family Fluctuation Analyses
#'
#'
fd_dfa <- function(y,
                   fs = NULL,
                   removeTrend = c("no","poly","adaptive","bridge")[2],
                   polyOrder=1,
                   standardise = c("none","mean.sd","median.mad")[2],
                   adjustSumOrder = FALSE,
                   scaleMin = 2,
                   scaleMax = floor(log2(NROW(y)/2)),
                   scaleResolution = (scaleMax-scaleMin),
                   scaleS = NA,
                   overlap = 0,
                   minData = 4,
                   doPlot = FALSE,
                   returnPlot = FALSE,
                   returnPLAW = FALSE,
                   returnInfo = FALSE,
                   silent = FALSE,
                   noTitle = FALSE,
                   tsName="y"){

  y_ori <- y

  if(!stats::is.ts(y)){
    if(is.null(fs)){fs <- 1}
    y <- stats::ts(y, frequency = fs)
    if(!silent){cat("\n\nfd_dfa:\tSample rate was set to 1.\n\n")}
  }

  if(is.na(scaleS)){
    scaleS <- unique(round(2^(seq(scaleMin, scaleMax, by=((scaleMax-scaleMin)/scaleResolution)))))
  }

  if(max(scaleS)>(NROW(y)/2)){
    scaleS <- scaleS[scaleS<=(NROW(y)/2)]
  }

  if(!all(is.numeric(scaleS),length(scaleS)>0,scaleS%[]%c(2,(NROW(y)/2)))){
    message("Something wrong with vector passed to scaleS.... \nUsing defaults: (scaleMax-scaleMin)/scaleResolution")
  }

  # Standardise by N
  if(any(standardise%in%c("mean.sd","median.mad"))){
    y <- ts_standardise(y, type = standardise,  adjustN = FALSE)
  }


  if(adjustSumOrder){
    y       <- ts_sumorder(y_ori, scaleS = scaleS, polyOrder = polyOrder, minData = minData)
    Hadj    <- attr(y,"Hadj")
    Hglobal <- attr(y,"Hglobal.excl")
  } else {
    Hadj    <- 0
    Hglobal <- NA
  }

  # Integrate the series
  if(standardise%in%"none"){
    y <- ts_integrate(ts_center(y)) # need negative values in profile
  } else {
    y <- ts_integrate(y)
  }

  TSm    <- as.matrix(cbind(t=1:NROW(y),y=y))
  DFAout <- monoH(TSm = TSm, scaleS = scaleS, polyOrder = polyOrder, returnPLAW = TRUE, returnSegments = TRUE)

  fitRange <- which(lapply(DFAout$segments,NROW)>=minData)

  lmfit1 <- stats::lm(DFAout$PLAW$bulk.log2 ~ DFAout$PLAW$size.log2, na.action=stats::na.omit)
  H1     <- lmfit1$coefficients[2] + Hadj
  lmfit2 <- stats::lm(DFAout$PLAW$bulk.log2[fitRange] ~ DFAout$PLAW$size.log2[fitRange], na.action=stats::na.omit)
  H2     <- lmfit2$coefficients[2] + Hadj

  outList <- list(
    PLAW  =  DFAout$PLAW,
    fullRange = list(y = y,
                     sap = lmfit1$coefficients[2],
                     Hadj = Hadj,
                     H = H1,
                     FD = sa2fd_dfa(lmfit1$coefficients[2]),
                     fitlm1 = lmfit1,
                     method = paste0("Full range (n = ",NROW(DFAout$PLAW$size),")\nSlope = ",round(stats::coef(lmfit1)[2],2)," | FD = ",sa2fd_dfa(stats::coef(lmfit1)[2]))),
    fitRange  = list(y = y,
                     sap = lmfit2$coefficients[2],
                     H = H2,
                     Hadj = Hadj,
                     FD = sa2fd_dfa(lmfit2$coefficients[2]),
                     fitlm2 = lmfit2,
                     method = paste0("Exclude large bin sizes (n = ",NROW(fitRange),")\nSlope = ",round(stats::coef(lmfit2)[2],2)," | FD = ",sa2fd_dfa(stats::coef(lmfit2)[2]))),
    info = list(fullRange=lmfit1,fitRange=lmfit2,segments=DFAout$segments),
    plot = NA,
    analysis = list(
      name = "Detrended FLuctuation Analysis",
      logBaseFit = "log2",
      logBasePlot = "2")
    )


  if(doPlot|returnPlot){
    if(noTitle){
      title <- ""
    } else {
      title <- "log-log regression (DFA)"
    }
    if(doPlot){
    g <- plotFD_loglog(fd.OUT = outList, title = title, subtitle = tsName, logBase = "2",xlabel = "Scale", ylabel = "Detrended Fluctuation")
    if(returnPlot){
      outList$plot <- g
      }
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
#' @param minData Minimum number of time/data points inside a box for it to be included in the slope estimation (default = `2^scaleMin`)
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
#' @return The boxcount fractal dimension and the 'local' boscount fractal dimension
#'
#' @note This function was inspired by the `Matlab` function `boxcount.m` [written by F. Moisy](http://www.fast.u-psud.fr/~moisy/ml/boxcount/html/demo.html). `Fred Hasselman` adapted the function for `R` for the purpose of the unit square boxcount analysis for 1D time series. The original Matlab toolbox has more options and contains more functions (e.g. `1D` and `3D` boxcount).
#'
#' @export
#'
#' @examples
#'
#' fd_boxcount2D(y = rnorm(100))
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
                          minData = 2^(scaleMin+1),
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
      y       <- ts_sumorder(y, scaleS = 2^(4:floor(log2(NROW(y)/2))), polyOrder = polyOrder, minData = 4)
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

  fitRange <- which(r%[]%c(minData,maxData))
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
      title <- ""
    } else {
      title <- "log-log regression (2D Boxcount)"
    }
    if(doPlot){
      g <- plotFD_loglog(fd.OUT = outList, title = title, subtitle = tsName, logBase = "2",xlabel = "Box Size (log)", ylabel = "Box Count (-log)")

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
        title <- ""
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
        title <- ""
      } else {
        if(useSD){
          title <- "log-log regression (Allan Deviation)"
        } else {
          title <- "log-log regression (Allan Variance)"
        }
      }
      g <- plotFD_loglog(fd.OUT = outList, title = title, subtitle = tsName, logBase = "10",xlabel = "Scale", ylabel = "Allan Variance")
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
#' @param signal An input signal.
#' @param qq A vector containing a range of values for the order of fluctuation `q`.
#' @param mins Minimum scale to consider.
#' @param maxs Maximum scale to consider.
#' @param ressc rescc
#' @param m m
#'
#' @return output
#' @export
#'
#' @family Fluctuation Analyses
#'
fd_mfdfa <- function(signal,qq=c(-10,-5:5,10),mins=6,maxs=12,ressc=30,m=1){

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

  y        <- cumsum(signal-mean(signal))
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
    RMS_scale[[ns]] <-plyr::ldply(ts_slice(TSm,scale[ns]),function(sv){return(sqrt(mean(ts_detrend(sv[,2]))^2))})
    for(nq in seq_along(qq)){
      qRMS[[nq]][1:length(RMS_scale[[ns]]$V1)] <- RMS_scale[[ns]]$V1^qq[nq]
      Fq[[nq]][ns] <- mean(qRMS[[nq]][1:length(RMS_scale[[ns]]$V1)])^(1/qq[nq])
      if(is.infinite(log2(Fq[[nq]][ns]))){Fq[[nq]][ns]<-NA}
    }
    Fq[[which(qq==0)]][ns] <- exp(0.5*mean(log(RMS_scale[[ns]]^2)))
    if(is.infinite(log2(Fq[[which(qq==0)]][ns]))){Fq[[which(qq==0)]][ns]<-NA}
  }

  fmin<-1
  fmax<-which(scale==max(scale))
  #for(nq in seq_along(qq)){Hq[nq] <- stats::lm(log2(Fq[[nq]])~log2(scale))$coefficients[2]}
  Hq <-plyr::ldply(Fq,function(Fqs){stats::lm(log2(Fqs[fmin:fmax])~log2(scale[fmin:fmax]),na.action=stats::na.omit)$coefficients[2]})

  tq <- (Hq[,1]*qq)-1
  hq <- diff(tq)/diff(qq)
  Dq <- (qq[1:(length(qq)-1)]*hq) - (tq[1:(length(qq)-1)])

  #if(reload==TRUE){library(signal,verbose=FALSE,quietly=TRUE)}

  return(list(q=qq,Hq=Hq,tq=tq,hq=hq,Dq=Dq,Hglobal=Hglobal,Hadj=Hadj))
}


#' mono Hurst
#'
#' @param TSm TS matrix with 2 columns `t` (1st) and `y` (second)
#' @param scaleS scales to evaluate
#' @param polyOrder If numeric: order to use for polynomial detrendiing, if "adaptive" will use the best fitting polynomial.
#'
#' @export
#' @keywords internal
#'
monoH <- function(TSm, scaleS, polyOrder = 1, returnPLAW = FALSE, returnSegments = FALSE, removeRMSbelow = .Machine$double.eps){

  dfaRMS_scale <- vector("list",length(scaleS))

  if(max(scaleS)>NROW(TSm)/2){
    scaleS <- scaleS[scaleS<=NROW(TSm)/2]
  }

  F2 <- numeric(length(scaleS))

  for(ns in seq_along(scaleS)){
    dfaRMS <- plyr::ldply(ts_slice(TSm[,2],scaleS[ns]), function(sv){return(sqrt(mean(ts_detrend(sv,polyOrder=polyOrder)^2,na.rm = TRUE)))})
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

# Dynamic Complexity ----


#' @title Dynamic Complexity
#'
#' @description Calculates Dynamic Complexity, a complexity index for short and coarse-grained time series (Schiepek & Strunk, 2010; Schiepek, 2003; Haken & Schiepek 2006).
#'
#' @param df A dataframe containing multivariate time series data from 1 person. Rows should indicate time, columns should indicate the time series variables. All time series in `df` should be on the same scale, an error will be thrown if the range of the time series in`df` is not `[scale_min,scale_max]`.
#' @param win Size of window in which to calculate Dynamic Complexity. If `win < NROW(df)` the window will move along the time series with a stepsize of `1` (default = `NROW(df)`)
#' @param scale_min The theoretical minimum value of the scale. Used to calculate expected values, so it is important to set this to the correct value.
#' @param scale_max The theoretical maximum value of the scale. Used to calculate expected values, so it is important to set this to rhe correct value.
#' @param doPlot If `TRUE` shows a Complexity Resonance Diagram of the Dynamic Complexity and returns an invisible [ggplot2::ggplot()] object. (default = `FALSE`)
#' @param doPlotF If `TRUE` shows a Complexity Resonance Diagram of the Fluctuation Intensity and returns an invisible [ggplot2::ggplot()] object. (default = `FALSE`)
#' #' @param doPlotD If `TRUE` shows a Complexity Resonance Diagram of the Distribution Uniformity and returns an invisible [ggplot2::ggplot()] object. (default = `FALSE`)
#' @param returnFandD Returns a list object containing the dynamic complexity series as well as the `F` and `D` series. (default = `FALSE`)
#' @param useVarNames Use the column names of `df` as variable names in the Complexity Resonance Diagram (default = `TRUE`)
#' @param colOrder If `TRUE`, the order of the columns in `df` determines the of variables on the y-axis. Use `FALSE` for alphabetic/numeric order. Use `NA` to sort by by mean value of Dynamic Complexity (default = `TRUE`)
#' @param useTimeVector Parameter used for plotting. A vector of length `NROW(df)`, containing date/time information (default = `NA`)
#' @param timeStamp If `useTimeVector` is not `NA`, a character string that can be passed to [lubridate::stamp()] to format the the dates/times passed in `useTimeVector` (default = `"01-01-1999"`)
#'
#' @return If `doPlot = TRUE`, a list object containing a data frame of Dynamic Complexity values and a `ggplot2` object of the dynamic complexity resonance diagram (e.g. Schiepek et al., 2016). If `doPlot = FALSE` the data frame with Dynamic Complexity series is returned.
#'
#' @export
#'
#' @author Merlijn Olthof
#' @author Fred Hasselman
#'
#' @family Dynamic Complexity functions
#'
#' @references Schiepek, G., & Strunk, G. (2010). The identification of critical fluctuations and phase transitions in short term and coarse-grained time series-a method for the real-time monitoring of human change processes. Biological cybernetics, 102(3), 197-207.
#' @references Schiepek, G. (2003). A Dynamic Systems Approach to Clinical Case Formulation. European Journal of Psychological Assessment, 19, 175-184.
#' @references Haken, H. & Schiepek, G. (2006, 2. Aufl. 2010). Synergetik in der Psychologie. Selbstorganisation verstehen und gestalten. G?ttingen: Hogrefe.
#' @references Schiepek, G. K., St?ger-Schmidinger, B., Aichhorn, W., Sch?ller, H., & Aas, B. (2016). Systemic case formulation, individualized process monitoring, and state dynamics in a case of dissociative identity disorder. Frontiers in psychology, 7, 1545.
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

#' @title Cumulative Complexity Peaks (CCP)
#'
#' @description Computes significant peaks in the dynamic complexity time series. Example: Schiepek, Tominschek & Heinzel, 2014.
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
#' @references Schiepek, G., & Strunk, G. (2010). The identification of critical fluctuations and phase transitions in short term and coarse-grained time series-a method for the real-time monitoring of human change processes. Biological cybernetics, 102(3), 197-207.
#' @references Schiepek, G. (2003). A Dynamic Systems Approach to Clinical Case Formulation. European Journal of Psychological Assessment, 19, 175-184.
#' @references Haken, H. & Schiepek, G. (2006, 2. Aufl. 2010). Synergetik in der Psychologie. Selbstorganisation verstehen und gestalten. G?ttingen: Hogrefe.
#' @references Schiepek, G. K., Tominschek, I., & Heinzel, S. (2014). Self-organization in psychotherapy: testing the synergetic model of change processes. Frontiers in psychology, 5, 1089.
#'
dc_ccp = function(df_win, alpha_item = 0.05, alpha_time = 0.05, doPlot = FALSE, useVarNames = TRUE, colOrder = TRUE, useTimeVector = NA, timeStamp = "01-01-1999"){

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
    distance<-1
    fluctuation <- NA
    tsy <- as.numeric(data[,column])

    for (i in (1:(nrow(data)-win-1))){
      y <- NA
      fluct <- NA
      dist_next <- 1
      k <- NA

      for (j in (0:(win-2))){
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
      k[win-1] = 1
      k<-k[1:(win-1)]
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


#' Distribution Uniformity
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


# PLOTS -------------------------------------------------------------------
#
#
#' gg_theme
#'
#' @param type      One of `"clean"`, or `"noax"`
#'
#' @details Will generate a `"clean"` ggplot theme, or a theme without any axes (`"noax"`).
#'
#' Some scientific journals explicitly request the Arial font should be used in figures. This can be achieved by using `.afm` font format (see, e.g. http://www.pure-mac.com/font.html).
#'
#' @return A theme for `ggplot2`.
#' @export
#'
#' @examples
#' library(ggplot2)
#' g <- ggplot(data.frame(x = rnorm(n = 100), y = rnorm(n = 100)), aes(x = x, y = y)) + geom_point()
#' g + gg_theme()
#' g + gg_theme("noax")
gg_theme <- function(type=c("clean","noax")){

  if(length(type)>1){type <- type[1]}

  switch(type,
         clean = theme_bw(base_size = 16, base_family="sans") +
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

#' gg_plotHolder
#'
#' @return A blank `ggplot2` object that can be used in concordance with `grid.arrange`.
#' @export
#'
#' @examples
#' # Create a plot with marginal distributions.
#' library(ggplot2)
#' library(scales)
#'
#' df <- data.frame(x = rnorm(n = 100),
#'                  y = rnorm(n = 100),
#'                  group = factor(sample(x=c(0,1),
#'                  size = 100, replace = TRUE)))
#'
#' scatterP <- ggplot(df, aes(x = x, y =y, colour = group)) +
#'                    geom_point() +
#'                    gg_theme()
#'
#' xDense <- ggplot(df, aes(x = x, fill = group)) +
#'                  geom_density(aes(y= ..count..),trim=FALSE, alpha=.5) +
#'                  gg_theme("noax") +
#'                  theme(legend.position = "none")
#'
#' yDense <- ggplot(df, aes(x = y, fill = group)) +
#'                  geom_density(aes(y= ..count..),trim=FALSE, alpha=.5) +
#'                  coord_flip() +
#'                  gg_theme("noax") +
#'                  theme(legend.position = "none")
#'
#' library(gridExtra)
#' grid.arrange(xDense,
#'              gg_plotHolder(),
#'              scatterP,
#'              yDense,
#'              ncol=2, nrow=2,
#'              widths=c(4, 1.4),
#'              heights=c(1.4, 4))
gg_plotHolder <- function(){
  return(ggplot2::ggplot() +
    ggplot2::geom_blank(ggplot2::aes_(1,1)) +
      ggplot2::theme(line = element_blank(),
          text  = element_blank(),
          title = element_blank(),
          plot.background = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  )
}



#' Plot Multivariate Time Series Data
#'
#' @param df A data frame with time series in columns.
#' @param timeVec If numeric, the number of the column in `df` which contains a time=keeping variable. If `NA`, the time vector will be `1:NROW(df)` (default = `NA`)
#' @param groupVec A vector indicating the names of the time series in the columns of `df`. If `NA`, the column names of `df` will be used, excluding the `timeVec`, if present. (default = `NA`)
#' @param doPlot Plot the [ggplot] object (default = `TRUE`)
#' @param returnPlotData Return the restructered data frame used to create the plot (default = `FALSE`)
#'
#' @return A [ggplot] object.
#' @export
#'
#' @examples
#'
#' # Generate some coloured noise
#' N <- 512
#' noises <- seq(-3,3,by=.5)
#' y <- data.frame(matrix(rep(NA,length(noises)*N), ncol=length(noises)))
#'
#' for(c in seq_along(noises)){
#'  y[,c] <- noise_powerlaw(N=N, alpha = noises[c])
#'  }
#'  colnames(y) <- paste0(noises)
#'
#'  plotTS_multi(y)
#'
plotTS_multi <- function(df,timeVec = NA, groupVec = NA, doPlot = TRUE, returnPlotData = FALSE){

  if(is.na(timeVec)){
    df$time <- 1:NROW(df)
  } else {
    if(is.numeric(timeVec)&length(timeVec)==1){
      colnames(df)[timeVec] <- "time"
    } else {
      stop("Argument timeVec is not interpretable as a column number.")
    }
  }
  if(is.na(groupVec)){
    groupVec <- colnames(df)[-c("time"%ci%df)]
  } else {
    if(length(groupVec)!=(NCOL(df)-1)){
      stop("Length of groupVec does not match number of time series in df.")
    }
  }

  tmp <- tidyr::gather(df, key = timeSeries, value = y, -time, factor_key = TRUE)
  tmp$timeSeries <- ordered(tmp$timeSeries)

  yOrder <- groupVec
  names(yOrder) <- paste(groupVec)
  offsets     <- names(yOrder) %>% {setNames(0:(length(.) - 1), .)}
  tmp$offsets <- unlist(llply(seq_along(yOrder), function(n) rep(offsets[n],sum(tmp$timeSeries%in%names(offsets)[n]))))

  # Calculate and scale group densities
  pdat <- tmp %>%
    group_by(timeSeries) %>%
    mutate(y = elascer(y,lo = -.45,hi = .45)) %>%
    ungroup()
  pdat$y_offset <- pdat$y + pdat$offsets

  if(doPlot){
    g <- ggplot(pdat, aes(x=time, y = y_offset, group = timeSeries)) +
      geom_path() +
      scale_y_continuous("Spectral Slope", breaks = offsets, labels = names(offsets), expand = c(0,0)) +
      scale_x_continuous("Time",expand = c(0,0)) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8, face = "plain"),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            plot.margin = margin(1,1,1,1,"line"))
    print(g)
  }

  if(returnPlotData){
    return(list(data = df,
                plot = invisible(g)))
  } else {
    return(invisible(g))
  }

}


#' Set Edge weights by group
#'
#'  Use a layout which takes a `weights`
#'
#' @param g  An igraph object whose edges (`get.edgelist(g)`) will be re-weighted according to the `membership` argument.
#' @param groups A named numeric vector with `length(V(g))` integers representing each group, or, a named character vector describing each group. If `names(groups)==NULL` then the names of the vector will be set as `names(groups) == V(g)$name`. If `V(g)$name==NULL`, the names of the vector will be set by the Vertex index
#' @param weigth.within The weight within a group (`default = 100`)
#' @param weight.between The weight within a group (`default = 1`)
#' @param preserve.weight.within If `E(g)$weights` is not `NULL`, try to preserve edge weigths within a group
#' @param preserve.weight.between If `E(g)$weights` is not `NULL`, try to preserve edge weigths between a groups
#' @param doPlot Plot the igraph object
#' @param returnOnlyWeights Do not return the graph, just the weights. If `FALSE` this will return the graph object, otherwis it returns `E(g)$weights`
#'
#' @return A numeric vector with `length(get.edgelist(g))` edge weights that will cluster groups defined in `membership` if a layout is used that can handle edge weights as a parameter (see examples).
#'
#' @export
#'
#' @family tools for plotting networks
#'
#' @examples
#' # Make a star graph and let the odd numbers cluster together
#' library(igraph)
#' g <-make_full_graph(10, directed=FALSE)
#' E(g)$width <- 3
#' V(g)$name <- paste(1:10)
#' membership <- rep(c(1,2),5)
#' names(membership) <- V(g)$name
#' E(g)$weight <- plotNET_groupWeight(g,membership,1000,10)
#' g$layout=layout.fruchterman.reingold(g,weights=E(g)$weight)
#' plot(g)
#'
#' # Make 3 groups by changing the 'membership' vector
#' membership[3:6] <- 3
#' names(membership) <- V(g)$name
#' E(g)$weight <- plotNET_groupWeight(g,membership,1000,10)
#' g$layout=layout.fruchterman.reingold(g,weights=E(g)$weight)
#' plot(g)
#'
#' # Use plotNET_groupColour for Vertex and Edge group colours
#' g <- plotNET_groupColour(g, membership, colourE=TRUE)
#' plot(g)
#'
plotNET_groupWeight <- function(g, groups, weigth.within=100, weight.between=1, preserve.weight.within=FALSE, preserve.weight.between=FALSE, doPlot = FALSE, returnOnlyWeights = TRUE){
  edgl <- igraph::get.edgelist(g)

  for(r in seq_along(edgl[,1])){
    row <- edgl[r,]
    if(as.numeric(groups[which(names(groups)==row[1])])==as.numeric(groups[which(names(groups)==row[2])])){
      igraph::E(g)$weight[r] <- weigth.within + ifelse(preserve.weight.within, igraph::E(g)$weight[r]%00%0, 0)
    } else {
      igraph::E(g)$weight[r] <- weight.between + ifelse(preserve.weight.between, igraph::E(g)$weight[r]%00%0, 0)
    }
  }
  if(doPlot){
    # graphics::plot.new()
    graphics::plot(g)
  }
  if(returnOnlyWeights){
    return(igraph::E(g)$weight)
  } else {
    return(invisible(g))
  }
}


#' Plot Network Based on RQA
#'
#' @param g An igraph object
#' @param labels Vertex labels
#' @param nodesize Set nodesizes by `degree(g, normalised = TRUE)` (default), `hubscore(g)$vector`, or, `strength(g)`, `eccentricity(g)`, `coreness(g)`. If a numeric value is passed all vertex sizes will be set to that value.
#' @param labelsize Set labelsize: "asnodesize" sets the `cex` for the labels to coincide with nodesize (with min of .4 and max of 1.1). A single numeric value sets the `cex` of all labels to that value. A numeric vector of length two, `c(min,max)` wil scale the node sizes to `min` and `max` which
#' @param edgeweight Set size of edges to `"E(g)$weight"` by passing "weight". If a single numeric value is provided all edges will be set to that value.
#' @param removeZeroDegree Remove vertices with `degree(g) == 0` (default = `TRUE`)
#' @param removeSelfLoops Calls `simplify(g)` (default = `TRUE`)
#' @param doPlot Plot the igraph object.
#'
#' @return an igraph object
#' @export
#'
#' @family tools for plotting networks
#'
plotNET_prep <- function(g,
                         labels     = NA,
                         nodesize   = c("degree","hubscore","strength","eccentricity","coreness")[1],
                         labelsize  = "asnodesize",
                         edgeweight = "weight",
                         removeZeroDegree = TRUE,
                         removeSelfLoops  = TRUE,
                         doPlot     = TRUE){

  if(removeSelfLoops){
    g <- simplify(g)
  }
  if(removeZeroDegree){
    g <- delete.vertices(g,degree(g)==0)
  }

  rev <- NA
  if(is.character(nodesize)){
  switch(nodesize,
         # degree   = rev <- elascer(log1p(igraph::degree(g,normalized = TRUE))),
         # hubscore = rev <- elascer(log1p(igraph::hub_score(g)$vector))
         # strength = rev <- elascer(log1p(igraph::strength(g)$vector)))
         degree   = rev <- igraph::degree(g),
         hubscore = rev <- igraph::hub_score(g)$vector,
         strength = rev <- igraph::strength(g),
         eccentricity = rev <- igraph::eccentricity(g),
         coreness = rev <- igraph::coreness(g),
           )
  } else {
    rev <- rep(as.numeric(nodesize),length(igraph::V(g)))
    }

  # set colors and sizes for vertices
  #rev<-elascer(log1p(igraph::V(g)$degree))

  igraph::V(g)$size        <- rev

  rev <- rev/max(rev, na.rm = TRUE)
  rev[rev<=0.2]<-0.2
  rev[rev>=0.9]<-0.9
  igraph::V(g)$rev <- rev

  igraph::V(g)$color       <- grDevices::rgb(igraph::V(g)$rev, 1-igraph::V(g)$rev,  0, 1)

  # set vertex labels and their colors and sizes
  if(all(is.na(labels))){
  igraph::V(g)$label       <- ""
  } else {
    igraph::V(g)$label       <- labels

   if(labelsize == "asnodesize"){igraph::V(g)$label.cex <- elascer(igraph::V(g)$size,lo = .4, hi = 1.1)}
   if(is.numeric(labelsize)&length(labelsize)==1){igraph::V(g)$label.cex <- labelsize}
   if(is.numeric(labelsize)&length(labelsize)==2){igraph::V(g)$label.cex <-  elascer(igraph::V(g)$size, lo = labelsize[1], hi = labelsize[2])}

    igraph::V(g)$label.color <- "black"

    igraph::V(g)$label.family = "Helvetica"
  }

  if(igraph::ecount(g)>0){
  if(edgeweight%in%"weight"){
    if(igraph::is_weighted(g)){
      igraph::E(g)$width <- elascer(igraph::E(g)$weight,lo = .8, hi = 5)
    }
  } else {
    if(is.numeric(edgeweight)){
      igraph::E(g)$width <- as.numeric(edgeweight)
    } else {
      igraph::E(g)$width <- 1
    }
  }
  igraph::E(g)$color <- grDevices::rgb(0.5, 0.5, 0.5, 1)
  }

  if(doPlot){
    # graphics::plot.new()
    graphics::plot(g)
  }

  return(invisible(g))
}


#' Example of Strogatz-Watts small-world network
#'
#' A wrapper around [igraph::sample_smallworld()] with `dim=1`
#'
#' @param n Size of the lattice (integer)
#' @param k Neighbourhood size (integer)
#' @param p Rewiring probability (between `0` and `1`)
#' @param doPlot PLot the igraph object
#'
#' @return A Strogatz-Watts small-world igraph object
#'
#' @export
#'
#' @family tools for plotting networks
#'
#' @seealso [igraph::sample_smallworld()]
#'
plotNET_SW <- function(n=100,k=5,p=0.05, doPlot = TRUE){

  g <- igraph::sample_smallworld(1, n, k, p)
  g <- plotNET_prep(g,doPlot = FALSE)

  # igraph::V(g)$degree <- igraph::degree(g)
  #
  # # set colors and sizes for vertices
  # rev<-elascer(log1p(igraph::V(g)$degree))
  # rev[rev<=0.2]<-0.2
  # rev[rev>=0.9]<-0.9
  # igraph::V(g)$rev <- rev$x
  #
  # igraph::V(g)$color       <- grDevices::rgb(igraph::V(g)$rev, 1-igraph::V(g)$rev,  0, 1)
  # igraph::V(g)$size        <- 25*igraph::V(g)$rev
  #
  # # set vertex labels and their colors and sizes
  # igraph::V(g)$label       <- ""
  #
  # igraph::E(g)$color <- grDevices::rgb(0.5, 0.5, 0.5, 1)

  if(doPlot){
    # graphics::plot.new()
    graphics::plot(g)
  }
  return(invisible(g))
}

#' Example of Barabasi scale-free network
#'
#' A wrapper around [igraph::sample_pa()]
#'
#' @param n Number of vertices
#' @param pwr Power of preferential attachment
#' @param out.dist Degree distribution
#' @param doPlot Plot the igraph object
#'
#' @return A Barabasi scale-free igraph object
#' @export
#'
#' @family tools for plotting networks
#'
#' @seealso [igraph::sample_pa()]
#'
plotNET_BA <- function(n=100, pwr=1, out.dist=NULL, doPlot = TRUE){

  g <- igraph::sample_pa(n, power = pwr, out.dist=out.dist, directed=FALSE)
  g <- plotNET_prep(g,doPlot = FALSE)


  # igraph::V(g)$degree <- igraph::degree(g)

  # # set colors and sizes for vertices
  # rev<-elascer(log1p(igraph::V(g)$degree))
  # rev[rev<=0.2] <- 0.2
  # rev[rev>=0.9] <- 0.9
  # igraph::V(g)$rev <- rev$x
  #
  # igraph::V(g)$color    <- grDevices::rgb(igraph::V(g)$rev, 1-igraph::V(g)$rev,  0, 1)
  # igraph::V(g)$size     <- 25*igraph::V(g)$rev
  # # igraph::V(g)$frame.color <- grDevices::rgb(.5, .5,  0, .4)
  #
  # # set vertex labels and their colors and sizes
  # igraph::V(g)$label <- ""
  # igraph::E(g)$width <- 1
  # igraph::E(g)$color <- grDevices::rgb(0.5, 0.5, 0.5, 1)

  if(doPlot){
    # graphics::plot.new()
    graphics::plot(g)
  }
  return(invisible(g))
}


#' Vertex Group Colours
#'
#' Identify Vertex and/or Edge groups by colour.
#'
#' @param g An igraph object
#' @param groups A named numeric vector with `length(V(g))` integers representing each group, or, a named character vector describing each group. If `names(groups)==NULL` then the names of the vector will be set as `names(groups) == V(g)$name`. If `V(g)$name==NULL`, the names of the vector will be set by the Vertex index
#' @param colourV Colour Vertices based on `groups` (default = `TRUE`)
#' @param alphaV Set transparency for Vertices (default = `1`)
#' @param colourE Colour Edges based on `groups`. Edges connecting to vertices of the same group will be coloured as the group (default = `FALSE`)
#' @param alphaE Set transparency for Edges. A single numeric, or a vector of length `ecount(g)` (default = `0.8`)
#' @param groupColours A list of length `groups` with valid colour codes
#' @param defaultEdgeColour Default edge colour
#' @param doPlot Plot the igraph object
#'
#' @return An igraph object with vertices and/or edges coloured by groups listed in `groups`
#'
#' @export
#'
#' @family tools for plotting networks
#'
plotNET_groupColour <- function(g, groups, colourV=TRUE, alphaV=1, colourE=FALSE, alphaE=.8, groupColours=NULL, defaultEdgeColour = "grey70", doPlot = TRUE){

  unicolours <- unique(groupColours)
  unigroups  <- unique(groups)

  if(length(groups)==igraph::gorder(g)){
    if(is.null(names(groups))){
      if(is.character(groups)){
        names(groups) <- groups
      } else {
        names(groups) <- paste0(1:igraph::gorder(g))
      }
    }
  } else {
    if(length(unique(groups))>length(unique(groupColours))){
      stop("length(groups) must be equal to number of Vertices: gorder(g)")
    }
  }

  if(is.null(groupColours)){
    if(length(unigroups)<=12){
      groupColours <-  scales::brewer_pal(palette="Paired")(length(unigroups))
    } else {
      groupColours <- scales::gradient_n_pal(scales::brewer_pal(palette="Paired")(12))(seq(0, 1, length.out = length(unigroups)))
    }
    unicolours <- groupColours
  } else {
    if(length(groups)==length(groupColours)){
      if(length(unique(groupColours))<=length(unigroups)){
        unicolours <- unique(groupColours)
      } else {
        stop("Number of groups does not match number of colours!")
      }
    } else {
      if(length(unique(groups))>length(unique(groupColours))){
        stop("Length of groups vector does not match length of groupColour vector!")
      }
    }
  }

  # Add alpha to edges
  if(all(alphaE%[]%c(0,1))){
    if(length(alphaE)==1|length(alphaE)==ecount(g)){
      igraph::E(g)$alpha <- alphaE
    } else {
      stop("Length of vector alphaE is not equal to 1 or ecount(g)")
      }
  } else {
      stop("All alphaE must be in [0,1]")
    }

  if(all(alphaV%[]%c(0,1))){
    if(length(alphaV)==1|length(alphaV)==vcount(g)){
      igraph::V(g)$alpha <- alphaV
    } else {
      stop("Length of vector alphaV is not equal to 1 or vcount(g)")
    }
  } else {
    stop("All alphaV must be in [0,1]")
  }

  igraph::E(g)$color  <- defaultEdgeColour
  igraph::E(g)$colour <- igraph::E(g)$color

  for(c in seq_along(unigroups)){
    Vid <- groups==unigroups[c]
    if(sum(Vid)>0){

      igraph::V(g)[Vid]$group      <- unigroups[c]
      igraph::V(g)[Vid]$groupnum   <- c

      if(colourV){
        igraph::V(g)[Vid]$color  <- add_alpha(unicolours[c], alpha = igraph::V(g)[Vid]$alpha)
        igraph::V(g)[Vid]$colour <- igraph::V(g)[Vid]$color
      }

      # if(alphaV){
      #   igraph::V(g)[Vid]$color <- add_alpha(igraph::V(g)[Vid]$color, alpha = igraph::V(g)[Vid]$alpha)
      #   igraph::V(g)[Vid]$colour <- igraph::V(g)[Vid]$color
      # }

      # Get ids for the edges that connect this group
      Eid <- which(igraph::E(g)%in%igraph::E(g)[igraph::V(g)[Vid]%--%igraph::V(g)[Vid]])


      if(length(Eid)>0){

        igraph::E(g)[Eid]$group             <- unigroups[c]
        igraph::E(g)[Eid]$groupnum          <- c

        if(colourE){
          igraph::E(g)[Eid]$color   <- add_alpha(unicolours[c], alpha = igraph::E(g)[Eid]$alpha)
          igraph::E(g)[Eid]$colour  <- igraph::E(g)[Eid]$color
        } else {
          # Add a default colour and alphac
          igraph::E(g)[Eid]$color  <- add_alpha(defaultEdgeColour, alpha = igraph::E(g)[Eid]$alpha)
          igraph::E(g)[Eid]$colour <- igraph::E(g)[Eid]$color
        }
        # if(alphaE){
        #   igraph::E(g)[Eid]$color <- add_alpha(groupColours[c], alpha = igraph::E(g)[Eid]$alpha)
        # }

        } # edge IDs > 0
      } # group IDs > 0
  } # group loop

  if(doPlot){
    # graphics::plot.new()
    graphics::plot(g)
  }
  return(invisible(g))
}



#' Plot output from fluctuation analyses based on log-log regression
#'
#' @param fd.OUT Output from one of the `fd_` functions that use log-log regression to get scaling exponents.
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param xlabel x label
#' @param ylabel y label
#' @param logBase base of the log used
#'
#' @return A ggplot object
#'
#' @export
#'
plotFD_loglog <- function(fd.OUT, title="", subtitle="", xlabel="Bin size", ylabel="Fluctuation", logBase=NA){

  if(!all(c("PLAW","fullRange","fitRange")%in%names(fd.OUT))){
    stop("Object fd.OUT should have 3 fields: PLAW, fullRange and fitRange")
  }

  if(nchar(title)==0){
    title <- fd.OUT$analysis$name
  }
  if(is.na(logBase)){
    logBase <- fd.OUT$analysis$logBasePlot
  }

  if(logBase%in%"e"){
    logFormat <- "log"
    logBaseNum <- exp(1)
    yAcc <- .1
    xAcc <- 1
  }

  if(logBase%in%"2"){
    logFormat <- paste0("log",logBase)
    logBaseNum <- as.numeric(logBase)
    yAcc <- .1
    xAcc <- 1
  }

  if(logBase%in%"10"){
    logFormat <- paste0("log",logBase)
    logBaseNum <- as.numeric(logBase)
    yAcc <- .1
    xAcc <- .01
  }


 logLabels <- function (expr, format = force)
  {
    quoted <- substitute(expr)
    subs <- function(x) {
      do.call("substitute", list(quoted, list(.x = as.name(x))))
    }
    function(x) {
      x <- format(x)
      lapply(x, subs)
    }
  }

 if(fd.OUT$analysis$name%in%"2D boxcount of 1D curve"){

  logSlopes <- data.frame(x = c(fd.OUT$PLAW$size[1:NROW(fd.OUT[[2]]$fitlm1$fitted.values)],
                                fd.OUT$ullRange$fitRange,fd.OUT$fitRange$fitRange),
                          y = c(fd.OUT$PLAW$bulk[1:NROW(fd.OUT[[2]]$fitlm1$fitted.values)],
                                fd.OUT$PLAW$bulk[1:NROW(fd.OUT[[3]]$fitlm2$fitted.values)]),
                          Method = c(rep(fd.OUT[[2]]$method,NROW(fd.OUT[[2]]$fitlm1$fitted.values)),
                                     rep(fd.OUT[[3]]$method,NROW(fd.OUT[[3]]$fitlm2$fitted.values))))
  } else {

  logSlopes <- data.frame(x = c(fd.OUT$PLAW$size[1:NROW(fd.OUT[[2]]$fitlm1$fitted.values)],
                                fd.OUT$PLAW$size[1:NROW(fd.OUT[[3]]$fitlm2$fitted.values)]),
                          y = c(fd.OUT$PLAW$bulk[1:NROW(fd.OUT[[2]]$fitlm1$fitted.values)],
                                fd.OUT$PLAW$bulk[1:NROW(fd.OUT[[3]]$fitlm2$fitted.values)]),
                          Method = c(rep(fd.OUT[[2]]$method,NROW(fd.OUT[[2]]$fitlm1$fitted.values)),
                                     rep(fd.OUT[[3]]$method,NROW(fd.OUT[[3]]$fitlm2$fitted.values))))
  }

  g <- ggplot2::ggplot(data.frame(fd.OUT$PLAW), ggplot2::aes_(x=~size,y=~bulk), na.rm=TRUE) +
    ggplot2::geom_point() +
    ggplot2::ggtitle(label = title, subtitle = subtitle)

  # if(logBase=="e"){
  # evalT<-
  #   paste0('g <- g + ggplot2::geom_smooth(data = logSlopes,  ggplot2::aes_(x=~x,y=~y, colour = ~Method, fill = ~Method), method="lm", alpha = .2) + ggplot2::scale_x_continuous(name = "',paste0(xlabel," (",logFormat,")"),'", breaks = scales::trans_breaks(',logFormat,', function(x) ',logBaseNum,'^x), labels =  scales::trans_format(',logFormat,',scales::math_format(e^.x)), trans = scales::log_trans(base = ',logBaseNum,')) +
  #     ggplot2::scale_y_continuous(name = "',paste0(ylabel," (",logFormat,")"),'", breaks = scales::trans_breaks(',logFormat,', function(x) ',logBaseNum,'^x), labels =  scales::trans_format(',logFormat,',scales::math_format(e^.x)), trans = scales::log_trans(base = ',logBaseNum,')) + ggplot2::annotation_logticks()')
  #
  # eval(parse(text = evalT))
  # }
  #
  # if(logBase=="2"){
  #   evalT<-
  #     paste0('g <- g + ggplot2::geom_smooth(data = logSlopes,  ggplot2::aes_(x=~x,y=~y, colour = ~Method, fill = ~Method), method="lm", alpha = .2) + ggplot2::scale_x_continuous(name = "',paste0(xlabel," (",logFormat,")"),'", breaks = scales::trans_breaks(',logFormat,', function(x) ',logBaseNum,'^x), labels =  scales::trans_format(',logFormat,',scales::math_format(2^.x)), trans = scales::log_trans(base = ',logBaseNum,')) +
  #            ggplot2::scale_y_continuous(name = "',paste0(ylabel," (",logFormat,")"),'", breaks = scales::trans_breaks(',logFormat,', function(x) ',logBaseNum,'^x), labels =  scales::trans_format(',logFormat,',scales::math_format(2^.x)), trans = scales::log_trans(base = ',logBaseNum,')) + ggplot2::annotation_logticks()')
  #
  #   eval(parse(text = evalT))
  # }
  #

  breaksX <- unique(fd.OUT$PLAW$size[1:NROW(fd.OUT[[2]]$fitlm1$fitted.values)])
  breaksY <- unique(fd.OUT$PLAW$bulk[1:NROW(fd.OUT[[2]]$fitlm1$fitted.values)])

  if(length(breaksX)>10){breaksX <- breaksX[seq.int(1,length(breaksX),length.out = 10)]}
  if(length(breaksY)>10){breaksY <- breaksY[seq.int(1,length(breaksY),length.out = 10)]}

  if(logBase!="no"){
    g <- g +
      ggplot2::geom_smooth(data = logSlopes,  ggplot2::aes_(x=~x,y=~y, colour = ~Method, fill = ~Method), method="lm", alpha = .2) +
      ggplot2::scale_x_continuous(name = paste0(xlabel," (",logFormat,")"),
                                  breaks = breaksX,
                                  labels = scales::number_format(accuracy = xAcc),
                                  trans = scales::log_trans(base = logBaseNum)) +
      ggplot2::scale_y_continuous(name = paste0(ylabel," (",logFormat,")"),
                                  breaks = breaksY,
                                  labels = scales::number_format(accuracy = yAcc),
                                  trans = scales::log_trans(base = logBaseNum))
  } else {

 # if(logBase=="no"){

  g <- g +
      ggplot2::geom_vline(data = data.frame(x = c(NROW(fd.OUT[[2]]$fitlm1$fitted.values),NROW(fd.OUT[[3]]$fitlm2$fitted.values)),Method = c(fd.OUT[[2]]$method,fd.OUT[[3]]$method)),  ggplot2::aes_(xintercept=~x, colour = ~Method)) +
      ggplot2::scale_x_continuous(name = xlabel, breaks = breaksX) +
      ggplot2::scale_y_continuous(name = ylabel, breaks = breaksY)
  }

  g <- g +
    ggplot2::scale_color_manual(values = c("red3","steelblue")) +
    ggplot2::scale_fill_manual(values = c("red3","steelblue")) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor =  ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(margin = ggplot2::margin(t = 5,b = 5, unit = "pt"), vjust = .5),
                   plot.margin = ggplot2::margin(t = 5,b = 5, r = 5,l = 5, unit = "pt"))


    # graphics::plot.new()
    graphics::plot(g)

  return(invisible(g))
}


#' Surrogate Test
#'
#' @param surrogateValues Vector of measures based on surrogate time series
#' @param observedValue The measure obtained from the observed value
#' @param sides Is this a 1 or 2-sided test (default = `1`)
#' @param binWidth The size of the histogram bins. The default is to look for the max. number of digits and set the width to `1/10^(Ndigits-1)`. If integers are detectec width will be set to 1.
#' @param measureName Label for x-axis
#' @param title A title for the plot
#' @param doPlot Plot a histogram of the distribution (default = `TRUE`)
#' @param returnOnlyPvalue Do not return the graph, just the point p-value (default = `FALSE`)
#'
#' alpha Significance threshold for the test. This value is currently calculated from the data as \eqn{\frac{1}{rank}*Nsides}, setting it will not have an effect.
#'
#' @return A point p-value for the observed value, and/or a histogram of the distribution (`ggplot2` object).
#' @export
#'
plotSUR_hist <- function(surrogateValues,
                         observedValue,
                         sides = c("two.sided","greater","less")[1],
                         binWidth = NULL,
                         measureName = "",
                         title="",
                         doPlot = TRUE,
                         returnOnlyPvalue = FALSE){

  if(any(sides%in%c("two.sided","greater","less"))){
    nsides <- dplyr::case_when(sides == "two.sided" ~ 2,
                        sides %in% c("greater","less") ~ 1)
  } else {
    stop("Use one of: 'two.sided', 'greater' or 'less' for argument 'sides'")
  }

  vec <- sort(c(surrogateValues, as.vector(observedValue)))

  if(is.null(binWidth)){
    if(all(is.wholenumber(vec))){
      binWidth<-1
    } else {
      binWidth<-1/10^(max(nchar(gsub("(\\d[.])+","",signif(vec))),na.rm = TRUE))
      if(binWidth<mean(diff(vec))){
        binWidth<-round(mean(diff(vec)),3)
      }
    }
  }

  alpha <- signif(nsides/length(vec),1)

  rank_dist_bin_lab  <- ggplot2::cut_width(vec, width=binWidth, center = 0, closed = "left")
  if(length(levels(rank_dist_bin_lab))>10*length(vec)){warning("Too many levels? Change binWidth!")}
  rank_dist_bin      <- ggplot2::cut_interval(vec, n=length(levels(rank_dist_bin_lab)), labels=FALSE)
  #ggplot2::cut_interval(rank_dist_bin,n = length(unique(rank_dist_bin)),labels=F)
  low <- seq(vec[1],vec[length(vec)],length.out = length(vec))[floor(length(vec)/2)]
  hi  <- seq(vec[1],vec[length(vec)],length.out = length(vec))[ceiling(length(vec)/2)]
  rank_dist_mid     <- mean(low,hi)
  rank_dist_bin_mid <- rank_dist_bin[which(vec>=rank_dist_mid)[1]]
  rank_dist         <- rank(vec,ties.method = "max")
  rank_obs         <- unique(rank(vec,ties.method = "max")[which(vec==observedValue)])

  rank_obs_p <- NULL

  if(sides%in%"less"){
    rank_obs_p     <- (1/((length(vec)+1)-rank_obs)*nsides)%00%NA #*table(rank(rev(vec),ties.method = "max"))[rank_obs]
    #rank_obs_p     <- table(rank(rev(vec),ties.method = "max"))[rank_obs]
    #rank_obs       <- unique(rank(vec,ties.method = "min")[which(vec==observedValue)])
    dist_lines     <- c(min(rank_dist_bin),rank_dist_bin_mid)
  }

  if(sides%in%"two.sided"){
    if(observedValue<rank_dist_mid){
      rank_obs_p     <- (1/((length(vec)+1)-rank_obs)*nsides)%00%NA
    } else {
    rank_obs_p     <- ((1/rank_obs)*nsides)%00%NA
    }
    dist_lines    <- c(min(rank_dist_bin),rank_dist_bin_mid, max(rank_dist_bin))
  }

  if(sides%in%"greater"){
    rank_obs_p     <- ((1/rank_obs)*nsides)%00%NA #*table(rank(rev(vec),ties.method = "max"))[rank_obs]
    dist_lines     <- c(rank_dist_bin_mid, max(rank_dist_bin))
  }

  if(!is.null(rank_obs)){

  rank_dist_prob <- 1/rank_dist
  df.dist <- data.frame(value      = vec,
                   rank_dist       = rank_dist,
                   rank_dist_prob  = rank_dist_prob,
                   rank_dist_bin   = rank_dist_bin,
                   rank_dist_bin_label = rank_dist_bin_lab
                   )

  ornd   <- round(observedValue,3)
  prnd   <- round(rank_obs_p,3)
  pLabel <- list(bquote(plain(P)(X==.(ornd)) ==.(prnd) ~ ~ alpha == .(alpha)))
  obs    <- cbind.data.frame(obs = "obs",
                          x = unique(df.dist$rank_dist_bin[rank_obs]),
                          y = 0)

  breaks <- seq(1,max(df.dist$rank_dist_bin),by=5)
  labels <- levels(df.dist$rank_dist_bin_label)[breaks]
  # if(length(breaks)>=length(vec)/2){
  #   breaks <- breaks[round(seq(1,max(breaks),length.out=max(breaks)/2))]
  #   labels <- paste(signif(unique(vec)[breaks],3))
  #   if(!all(diff(diff(breaks))!=0)){
  #     breaks <- seq(1,max(df.dist$rank_dist_bin),by=2)
  #     if(max(df.dist$rank_dist_bin)!=max(breaks)){
  #     breaks <- c(breaks, max(breaks)+min(diff(breaks)))
  #     labels <- c(paste(signif(na.exclude(unique(vec)[breaks],3))),paste(signif(max(na.exclude(unique(vec)[breaks],3))+min(diff(na.exclude(unique(vec)[breaks]))),3)))
  #     }
  #   }
  # }


  g <- ggplot2:: ggplot(df.dist,aes_(x=~rank_dist_bin)) +
    geom_histogram(breaks= seq(1,(max(df.dist$rank_dist_bin)+1))-.5,colour="white") +
    geom_vline(xintercept = dist_lines, colour="steelblue") +
    geom_point(data=obs, aes_(x=~x, y=~y,colour=~obs), size=5) +
    scale_x_continuous(name = measureName, breaks = breaks, labels = labels) +
    scale_color_manual("Observed",values="red", breaks="obs", labels = pLabel) +
    ggtitle(label = title, subtitle = paste0(nsides,"-sided test with ",length(surrogateValues)," surrogate values. The observed value has (max) rank ",rank_obs,".")) +
    theme_bw() +
    theme(legend.position = c(1,.90),
          legend.justification = "right",
          legend.background = element_rect(colour = "black"),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle=90, vjust = 0.5))

    if(doPlot){
      grid::grid.newpage()
      grid::grid.draw(g)
    }

  attr(rank_obs_p,"sides")         <- sides
  attr(rank_obs_p,"rank_observed") <- rank_obs
  attr(rank_obs_p,"N_values")      <- length(vec)
  attr(rank_obs_p,"alpha")         <- alpha

  if(returnOnlyPvalue){
    return(rank_obs_p)
  }
  return(invisible(list(pvalue=rank_obs_p,surrogates_plot=g,plot_data=df.dist)))
  } else {
    stop("Could not determine a rank order probability. Did you provide a correct argument for sides?")
  }
}



#' Add transparency to a colour
#'
#' @param col A colour name, hexadecimal string or positive integer `i`, such that `palette()[i]`
#' @param alpha Alpha transparency value
#'
#' @return An rgb colour with transparency
#' @export
#'
add_alpha <- function(col, alpha=1){
  if(missing(col)){stop("Please provide a vector of colours.")}
  if(!is.null(dim(col))){stop("Please provide a 1D vector of hex colours")}
  if(length(unique(alpha))==1){
    cols <- apply(sapply(col, grDevices::col2rgb, USE.NAMES = FALSE)/255, 2, function(x){grDevices::rgb(x[1], x[2], x[3], alpha=unique(alpha))})
  } else {
    if(length(col)==length(alpha)){
    colsL <- list()
    uniAlpha <- unique(alpha)
    for(a in seq_along(uniAlpha)){
      IDs <- alpha%in%uniAlpha[a]
      colsL[[a]] <- apply(sapply(col[IDs], grDevices::col2rgb, USE.NAMES = FALSE)/255, 2, function(x){grDevices::rgb(x[1], x[2], x[3], alpha=uniAlpha[a])})
    }
    cols <- unlist(colsL)
    } else {
      stop("If length(alpha)>1 it must be the same as length(col)!")
    }
  }
  return(cols)
}


#' FisherZ for PACF
#'
#' @param r r
#' @param n n
#' @param lag lag
#' @param siglevel siglevel
#' @param sides sides
#' @param type type
#'
#' @return dataframe
#' @export
#'
#' @keywords internal
#'
pacf_fisherZ <-function(r, n, lag, siglevel=.05,sides=2,type=""){
  z <- atanh(r)/sqrt(1/(n-3))
  conf.low   <- tanh(atanh(r) - (stats::qnorm(1-(siglevel/sides)) * sqrt((1/(n-3)))))
  conf.high  <- tanh(atanh(r) + (stats::qnorm(1-(siglevel/sides)) * sqrt((1/(n-3)))))
  p<-2*(1-stats::pnorm(abs(z)))
  if(p<siglevel){
    sig <-"yes"
  } else {
    sig <-"no"
  }
  return(data.frame(r=r,ciL=conf.low,ciU=conf.high,fZ=z,p=p,alpha=siglevel,sig=sig,lag=lag,n=n,type=type, stringsAsFactors = FALSE))
}


#' Plot ACF and PACF
#'
#' @param y A time series or numeric vector
#' @param Lmax Maximum number of lags
#' @param alpha Significance level
#' @param doPlot Plot output
#' @param returnCorFun Return the data
#'
#' @return Either an invisible ggplot2 object r a list containing the plot and the data
#'
#' @family Plot redundancy functions
#'
#' @export
#'
plotRED_acf <- function(y, Lmax = max(round(NROW(y)/4),10),alpha=.05 ,doPlot = TRUE, returnCorFun = FALSE){

  siglevel <- alpha
  df.acf <- stats::acf(y,plot=FALSE, lag.max = Lmax)
  df.pacf <- stats::pacf(y,plot=FALSE, lag.max = Lmax)

  dfN <- c(NROW(y), plyr::laply(1:Lmax, function(l) NROW(ts_embed(y,2,l))+1))

  corfunACF  <- plyr::ldply(seq_along(df.acf$acf), function(cc){pacf_fisherZ(r=df.acf$acf[cc],n=dfN[cc],lag=df.acf$lag[cc],type="acf")})
  corfunPACF <- plyr::ldply(seq_along(df.pacf$acf), function(cc){pacf_fisherZ(r=df.pacf$acf[cc],n=dfN[cc],lag=df.pacf$lag[cc],type="pacf")})
  corfun     <- rbind(corfunACF,corfunPACF)

  groupColours <-  scales::brewer_pal(palette="RdBu")(11)
  cols <- c("yes"=groupColours[9],"no"=groupColours[3])

  g <- ggplot2::ggplot(corfun,ggplot2::aes_(x=~lag,y=~r)) +
    ggplot2::geom_hline(yintercept = 0, colour="grey",size=1) +
    ggplot2::geom_line(data = data.frame(x=c(0,corfun$lag[1]),y=c(1,corfun$r[1])),ggplot2::aes_(x=~x,y=~y),colour="grey50") +
    ggplot2::geom_point(x=0,y=1,colour=groupColours[10],fill=groupColours[9],size=2,pch=21) +
    ggplot2::geom_ribbon(aes_(ymin=~ciL,ymax=~ciU),fill="grey70",colour="grey50") +
    ggplot2::geom_path(colour="grey50") +
    ggplot2::geom_point(aes_(fill = ~sig, colour=~sig),pch=21, cex=(1 + .01*(NROW(y)/Lmax))) +
    ggplot2::facet_grid(type ~.) +
    ggplot2::scale_fill_manual(bquote(p < .(siglevel)),values = cols,
                      labels =  list("yes"= expression(rho != 0),
                                     "no" = expression(rho == 0))) +
    ggplot2::scale_colour_manual(bquote(p < .(siglevel)),values = cols,
                        labels =  list("yes"= expression(rho != 0),
                                       "no" = expression(rho == 0))) +
    ggplot2::scale_x_continuous(limits = c(0,Lmax),expand = c(0.01,0), breaks = seq(0,Lmax,by = round(Lmax/10))) +
    ggplot2::scale_y_continuous(limits = c(-1,1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank())

  if(doPlot){
    # graphics::plot.new()
    graphics::plot(g)
  }

  if(returnCorFun){
    return(list(corfun=corfun,
                plot=invisible(g)))
  } else {
    return(invisible(g))
  }
}

#' Plot various MI functions
#'
#'
#' @inheritParams mif
#' @param mif.OUT Output from function [mif()]
#' @param returnMIFun Return the data
#'
#' @return Either an invisible ggplot2 object r a list containing the plot and the data
#'
#' @family Plot redundancy functions
#'
#' @export
#'
plotRED_mif <- function(mif.OUT = NULL, lags = 0:max(round(NROW(y)/4),10), nbins = ceiling(2*NROW(y)^(1/3)), surTest = FALSE, alpha=.05 ,doPlot = TRUE, returnMIFun = TRUE){

  if(is.null(mif.OUT)){
      stop("No data.")
  } else {
    if(is.null(attr(mif.OUT,"miType"))){
      stop("Argument mif.OUT is not the output of function mif().")
    }
  }


  # mifunMIF  <- mif(y = y, lags = lags, nbins = nbins)
  # #ldply(seq_along(df.acf$acf), function(cc){pacf_fisherZ(r=df.acf$acf[cc],n=dfN[cc],lag=df.acf$lag[cc],type="acf")})
  # mifunPMIF <- mif(y = cbind(y,y[,1],y[,1]), lags = lags, nbins = nbins)
  # #ldply(seq_along(df.pacf$acf), function(cc){pacf_fisherZ(r=df.pacf$acf[cc],n=dfN[cc],lag=df.pacf$lag[cc],type="pacf")})
  # mif.OUT   <- rbind(mifunMIF,mifunPMIF)

  groupColours <-  scales::brewer_pal(palette="RdBu")(11)
  cols <- c("yes"=groupColours[9],"no"=groupColours[3])

  mifun_long <- data.frame(lag =  c(as.numeric(names(mif.OUT)),as.numeric(names(mifunPMIF))),
                           mi = c(mifunMIF,mifunPMIF),
                           type = c(rep(attributes(mifunMIF)$miType,NROW(mifunMIF)),rep(attributes(mifunPMIF)$miType,NROW(mifunPMIF))))

  g <- ggplot2::ggplot(mifun_long,ggplot2::aes_(x=~lag,y=~mi)) +
    ggplot2::geom_hline(yintercept = 0, colour="grey",size=1) +
    ggplot2::geom_line(data = data.frame(x=c(0,mifun_long$lag[1]),y=c(1,mifun_long$mi[1])),ggplot2::aes_(x=~x,y=~y),colour="grey50")

  if(length(lags)<=50){
   g <- g + ggplot2::geom_point(x=0,y=1,colour=groupColours[10],fill=groupColours[9],size=2,pch=21)
  }
    #geom_ribbon(aes_(ymin=~ciL,ymax=~ciU),fill="grey70",colour="grey50") +
   g <- g +  ggplot2::geom_path(colour="grey50") +
    ggplot2::geom_point(pch=21, cex=(1 + .01*(NROW(y)/nbins))) +
    ggplot2::facet_grid(type ~.) +
    # scale_fill_manual(bquote(p < .(siglevel)),values = cols,
    #                   labels =  list("yes"= expression(rho != 0),
    #                                  "no" = expression(rho == 0))) +
    # scale_colour_manual(bquote(p < .(siglevel)),values = cols,
    #                     labels =  list("yes"= expression(rho != 0),
    #                                    "no" = expression(rho == 0))) +
    ggplot2::scale_x_continuous(limits = c(0,length(lags)),expand = c(0.01,0), breaks = seq(0,length(lags),by = round(length(lags)/10))) +
    ggplot2::scale_y_continuous(limits = c(-1,1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank())

  if(doPlot){
    # graphics::plot.new()
    graphics::plot(g)
  }

  if(returnMIFun){
    return(list(mifun=mifun,
                plot=invisible(g)))
  } else {
    return(invisible(g))
  }
}



#' Plot Complexity Resonance Diagram
#'
#' @inheritParams dc_ccp
#' @param title A title for the plot.
#' @param subtitle A subtitle for the plot.
#' @param xlabel A label for the x-axis.
#' @param ylabel A label for the y-axis.
#'
#' @return An invisible ggplot2 object.
#' @export
#'
#' @family Dynamic Complexity functions
#'
plotDC_res <-  function(df_win, win, useVarNames = TRUE, colOrder = TRUE, useTimeVector = NA, timeStamp = "01-01-1999", doPlot = TRUE, title = 'Complexity Resonance Diagram', subtitle = "", xlabel = "Time", ylabel = ""){

  if(!useVarNames){
    df_win <- data.matrix(df_win)
    colnames(df_win) <- c(1:ncol(df_win))
  }

  if(win>=NROW(df_win)){
    stop("Only one value to display.")
  }

  if(!all(is.na(useTimeVector))){
    if(length(useTimeVector)==1){
      if(useTimeVector%in%colnames(df_win)){
        labels <- lubridate::stamp(timeStamp)(df_win[,useTimeVector%ci%df_win])
      }
    } else {
      if(length(useTimeVector)==NROW(df_win)){
        labels <- lubridate::stamp(timeStamp)(useTimeVector)
      }
    }
  } else {
    labels <- paste(1:NROW(df_win))
  }

  df_win$time <- 1:NROW(df_win)

  minorBreaks <- df_win$time[seq(win, NROW(df_win))]
  if(NROW(df_win)>50){
    labels <- labels[minorBreaks]
    #  by = round(length(breaks)/25))]
    #labels[!labels%in%labels[seq(2,length(minorBreaks), by = round(length(minorBreaks)/25))]] <- ""
    labels <- labels[c(seq(2,length(minorBreaks), by = round(length(minorBreaks)/25)))]
    majorBreaks <- minorBreaks[c(seq(2,length(minorBreaks), by = round(length(minorBreaks)/25)))]
  }

  # breaks <- seq(win, NROW(df_win))
  # if(NROW(df_win)>50){
  #   labels <- paste(breaks)
  #   labels[!labels%in%labels[seq(2,length(breaks), by = round(length(breaks)/25))]]<- ""
  #   labels[1] <- ""
  #   labels[2] <- paste(win+1)
  # }

  if(is.na(colOrder)){
    subtitle <- 'Variables ordered by mean Dynamic Complexity'
  } else {
    subtitle <- ifelse(colOrder,
                       'Variables ordered by position in data source',
                       'Variables ordered by variable name')
  }

  if(is.na(colOrder)){
    df_win <- dplyr::select(df_win,names(sort(colMeans(df_win, na.rm = TRUE))))
    colOrder <- TRUE
  }

  df_win$time <- 1:NROW(df_win)
  dfp <- tidyr::gather(df_win, key = "variable", value = "value", c(-.data$time), factor_key = colOrder)
  dfp$time <- as.numeric(dfp$time)


  # if(subtitle=="Variables ordered by ..."){
  #   subt <- paste0("Variables ordered by ",corder)
  # }
  #max(dfp$value, na.rm=TRUE)/2
  g <- ggplot2::ggplot(dfp, ggplot2::aes_(x=~time, y=~variable, fill=~value)) +
    ggplot2::geom_raster(interpolate = FALSE) +
    ggplot2::scale_fill_gradient2("Dynamic Complexity",low='steelblue', high='red3', mid='whitesmoke', midpoint=(max(dfp$value, na.rm=TRUE)/2), na.value='white') +
    ggplot2::geom_vline(xintercept=minorBreaks-.5, colour="steelblue", alpha=.9, size=.1) +
    ggplot2::geom_hline(yintercept=1:NROW(dfp)-.5, colour="steelblue", alpha=.9, size=.1) +
    ggplot2::scale_y_discrete(ylabel, expand = c(0,0)) +
    ggplot2::scale_x_continuous(xlabel,
                                breaks = majorBreaks,
                                minor_breaks = minorBreaks-.5,
                                labels = labels,
                                expand = c(0,0), limits = c((win+1)-.5, max(majorBreaks)+.5)) +
   # ggplot2::scale_x_continuous(xlabel, breaks = breaks, labels = labels, expand = c(0,0), limits = c((win+1)-.5, max(breaks)+.5)) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_bw() +
    ggplot2::theme(#axis.text.x = element_text(vjust = 0, hjust = 1, angle = 90),
                   #axis.text.y = element_text(vjust = 1),
                   legend.position = "bottom",
                   legend.key.size = unit(.05,"npc"),
                   axis.text.x = element_text(size=8, angle = 90, vjust =0.5),
                   axis.text.y = element_text(size=8),
                   axis.line.x.top = element_line(),
                   axis.ticks.x.top = element_line(),
                   panel.grid.major.x = element_blank(),
                   panel.grid.minor.x = element_line(),
                   panel.grid.minor.y = element_line(),
                   panel.grid.major.y = element_blank(),
                   plot.subtitle = element_text(size=10),
                   plot.title = element_text(face = "bold"),
                   legend.title = element_text(vjust=.75,size = 10))


  if(doPlot){
   graphics::plot(g)
  }
  return(invisible(g))
}



#' Plot Cumulative Complexity Peaks
#'
#' @param df_ccp A dataframe generated by `dc_ccp()`
#' @inheritParams plotDC_res
#'
#' @return An invisible ggplot2 object.
#' @export
#'
#' @family Dynamic Complexity functions
#'
plotDC_ccp <-  function(df_ccp, win, useVarNames = TRUE, colOrder = TRUE, useTimeVector = NA, timeStamp = "31-01-1999", doPlot = TRUE, title = 'Critical Instability Plot', subtitle = "", xlabel = "Time", ylabel = ""){

  if(!useVarNames){
    df_ccp <- data.matrix(df_ccp)
    colnames(df_ccp) <- c(1:ncol(df_ccp))
  }

  lcol  <- df_ccp$sig.peaks
  lcol[lcol==0] <- 5
  df_ccp <- df_ccp[,-c("sig.peaks"%ci%df_ccp)]
  if(is.na(colOrder)){
    df_ccp <- dplyr::select(df_ccp,names(sort(colSums(df_ccp, na.rm = TRUE))))
    colOrder <- TRUE
  }

   if(!all(is.na(useTimeVector))){
     if(length(useTimeVector)==1){
       if(useTimeVector%in%colnames(df_ccp)){
         labels <- lubridate::stamp(timeStamp)(df_ccp[,useTimeVector%ci%df_ccp])
         }
       } else {
         if(length(useTimeVector)==NROW(df_ccp)){
           labels <- lubridate::stamp(timeStamp)(useTimeVector)
         }
       }
   } else {
     labels <- paste(1:NROW(df_ccp))
   }

  df_ccp$time <- 1:NROW(df_ccp)

  minorBreaks <- df_ccp$time[seq(win, NROW(df_ccp))]
  if(NROW(df_ccp)>50){
    labels <- labels[minorBreaks]
    labels <- labels[c(seq(2,length(minorBreaks), by = round(length(minorBreaks)/25)))]
    majorBreaks <- minorBreaks[c(seq(2,length(minorBreaks), by = round(length(minorBreaks)/25)))]
  }

 if(!colOrder){
   df_ccp <- df_ccp[,sort(colnames(df_ccp))]
 }

  df_ccp$sig.peaks <- lcol
  colnames(df_ccp)[colnames(df_ccp)%in%"sig.peaks"]<-"Sig. CCP"
  dfp <- tidyr::gather(df_ccp, key = "variable", value = "value", -(time), factor_key = TRUE)
  dfp$value <- factor(dfp$value,levels = c(0,1,5,10),
                      labels = c("0","Sig. DC level","5","Sig. CCP"))


  g <- ggplot2::ggplot(dfp[stats::complete.cases(dfp),], ggplot2::aes_(x=~time, y=~variable, fill=~value)) +
    ggplot2::geom_raster(interpolate = FALSE) +
    ggplot2::geom_vline(xintercept=minorBreaks-.5, colour="grey90", alpha=1, size=.1) +
    ggplot2::geom_hline(yintercept=1:NROW(dfp)-.5, colour="grey90", alpha=1, size=.1) +
    ggplot2::scale_y_discrete(ylabel, expand = c(0,0)) +
    ggplot2::scale_x_continuous(xlabel,
                                breaks = majorBreaks,
                                minor_breaks = minorBreaks-.5,
                                labels = labels,
                                expand = c(0,0), limits = c((win+1)-.5, max(majorBreaks)+.5)) +
    ggplot2::scale_fill_manual("Critical Insability", breaks = c("Sig. DC level","Sig. CCP"), values = c("0"="white","Sig. DC level"="grey","5"="whitesmoke","Sig. CCP"="black")) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.size = unit(.05,"npc"),
      axis.text.x = element_text(size=8, angle = 90, vjust =0.5),
      axis.text.y = element_text(size=8),
      axis.line.x.top = element_line(),
      axis.ticks.x.top = element_line(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_line(),
      panel.grid.minor.y = element_line(),
      panel.grid.major.y = element_blank(),
      plot.subtitle = element_text(size=10),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(vjust=.75,size = 10))

  if(doPlot){
    # graphics::plot.new()
    graphics::plot(g)
  }

  return(invisible(g))
}



#' Plot Peaks versus Levels
#'
#' Produce a plot in which the output of `dc_win()` and `dc_ccp()` on the same multivariate timeseries data is combined with the output of `ts_level()` on a state variable of the same length as the multivariate data.
#'
#' @param df_lvl A dataframe generated by `ts_level()` of a variable that is considered a state variable.
#' @inheritParams plotDC_ccp
#' @inheritParams plotDC_res
#' @param levelName A name for the state variable.
#'
#' @return An invisible ggplot2 object.
#' @export
#'
#' @family Dynamic Complexity functions
#'
plotDC_lvl <-  function(df_win, df_ccp, df_lvl, win, useVarNames = TRUE, colOrder = TRUE, useTimeVector = NA, timeStamp = "31-01-1999", doPlot = TRUE, title = 'Peaks versus Levels Plot', subtitle = "", xlabel = "Time", ylabel = "", levelName = "State variable"){

  if(!(all(NROW(df_win)==NROW(df_ccp) & NROW(df_win)==NROW(df_lvl$pred) & NROW(df_ccp)==NROW(df_lvl$pred)))){
    stop("The time series must all have equal lengths:\n all(NROW(df_win)==NROW(df_ccp) & NROW(df_win)==NROW(df_lvl$pred) & NROW(df_ccp)==NROW(df_lvl$pred))")
  }

  if(!(NCOL(df_win)==(NCOL(df_ccp)-1))){
    stop("There must be an equal number of time series variables in 'df_win' and 'df_ccp'")
  }

  ylabel <- paste0(ylabel," (arbitrary units)")

  if(!useVarNames){
    df_ccp <- data.matrix(df_ccp)
    colnames(df_ccp) <- c(1:ncol(df_ccp))

    df_win <- data.matrix(df_win)
    colnames(df_win) <- c(1:ncol(df_win))
  }

  lcol  <- df_ccp$sig.peaks
  lcol[lcol==0] <- 5
  df_ccp <- df_ccp[,-c("sig.peaks"%ci%df_ccp)]
  if(is.na(colOrder)){
    df_ccp <- dplyr::select(df_ccp,names(sort(colSums(df_ccp, na.rm = TRUE))))
    colOrder <- TRUE
  }

  if(!all(is.na(useTimeVector))){
    if(length(useTimeVector)==1){
      if(useTimeVector%in%colnames(df_ccp)){
        labels <- lubridate::stamp(timeStamp)(df_ccp[,useTimeVector%ci%df_ccp])
      }
    } else {
      if(length(useTimeVector)==NROW(df_ccp)){
        labels <- lubridate::stamp(timeStamp)(useTimeVector)
      }
    }
  } else {
    labels <- paste(1:NROW(df_ccp))
  }

  df_ccp$time <- 1:NROW(df_ccp)

  minorBreaks <- df_ccp$time[seq(win, NROW(df_ccp))]
  if(NROW(df_ccp)>50){
    labels <- labels[minorBreaks]
    labels <- labels[c(seq(2,length(minorBreaks), by = round(length(minorBreaks)/25)))]
    majorBreaks <- minorBreaks[c(seq(2,length(minorBreaks), by = round(length(minorBreaks)/25)))]
  }

  if(!colOrder){
    df_ccp <- df_ccp[,sort(colnames(df_ccp))]
  }

  dc.mean <- data.frame(time   = df_ccp$time,
                        meanDC = plyr::laply(1:NROW(df_win), function(r){
                          suppressWarnings(
                            sum(as.numeric(df_win[r,1:NCOL(df_win)]*df_ccp[r,1:(NCOL(df_ccp)-1)]), na.rm = TRUE)%00%NA
                            )
                          })
                        )
  dc.mean$meanDC[dc.mean$meanDC==0] <- NA

  # Make sure the CCP's turn up as max value
  #dc.mean$meanDC[df_ccp$sig.peaks==10] <- max(dc.mean$meanDC, na.rm = TRUE)
  dc.mean$meanDC <- elascer(dc.mean$meanDC)
  colnames(dc.mean)[2] <- "Sum of Sig. DC peaks"

  dc.lvl     <- data.frame(lvl = rep(NA,NROW(df_lvl$pred)))
  dc.lvl$lvl <- elascer(df_lvl$pred$p, lo = 0.25, hi = 0.75)
  colnames(dc.lvl) <- levelName
  lvl <- cbind(dc.mean, dc.lvl) %>% gather(key = "variable", value = "DC", -.data$time)
  lvl$y1 <- NA
  lvl$y2 <- NA

  ccp <- data.frame(ccp = df_ccp$time[lcol==10])

  lvl <- dplyr::add_row(lvl, time = ccp$ccp, variable = "Sig. CCP", DC = NA, y1 = 1.1, y2 = -0.1)

  cols <- eval(parse(text=paste0("c('Sum of Sig. DC peaks' = 'black','",levelName,"' = 'red3', 'Sig. CCP' = 'steelblue')")))

  g <- ggplot(lvl, aes_(x = ~time, y = ~DC, colour = ~variable)) +
    #geom_path(data = lvl, aes(x= time, y = lvl), colour = "red3", size = 1)  +
    geom_path(size=1) +
    geom_vline(data = ccp, aes_(xintercept=~ccp), colour = "steelblue") +
    geom_point(aes_(y=~y1), pch = 25, fill = "steelblue", size =3, show.legend = FALSE) +
    geom_point(aes_(y=~y2), pch = 24, fill = "steelblue", size =3, show.legend = FALSE) +
    scale_y_discrete(ylabel, breaks = NULL, labels = "", expand = c(0,0)) +
    scale_x_continuous(xlabel,
                       breaks = majorBreaks,
                       minor_breaks = minorBreaks-.5,
                       labels = labels,
                       expand = c(0,0), limits = c((win+1)-.5, max(majorBreaks)+.5)) +
    scale_colour_manual("", values = cols) +
    labs(title = title, subtitle = subtitle) +
    theme_bw() + theme(
      legend.position = "bottom",
      legend.key.size = unit(.05,"npc"),
      axis.text.x = element_text(size=8, angle = 90, vjust =0.5),
      axis.line.x.top = element_line(),
      axis.ticks.x.top = element_line(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_line(),
      panel.grid.minor.y = element_line(),
      panel.grid.major.y = element_blank(),
      plot.subtitle = element_text(size=10),
      plot.title = element_text(face = "bold"),
      legend.title = element_text(vjust=.75,size = 10))


  # df_ccp$sig.peaks <- lcol
  # colnames(df_ccp)[colnames(df_ccp)%in%"sig.peaks"]<-"Sig. CCP"
  # dfp <- tidyr::gather(df_ccp, key = variable, value = value, -(time), factor_key = TRUE)
  # dfp$value <- factor(dfp$value,levels = c(0,1,5,10),
  #                     labels = c("0","Sig. DC level","5","Sig. CCP"))
  #
  #
  # g <- ggplot2::ggplot(dfp[complete.cases(dfp),], ggplot2::aes_(x=~time, y=~variable, fill=~value)) +
  #   ggplot2::geom_raster(interpolate = FALSE) +
  #   ggplot2::geom_vline(xintercept=minorBreaks-.5, colour="grey90", alpha=1, size=.1) +
  #   ggplot2::geom_hline(yintercept=1:NROW(dfp)-.5, colour="grey90", alpha=1, size=.1) +
  #   ggplot2::scale_y_discrete(ylabel, expand = c(0,0)) +
  #   ggplot2::scale_x_continuous(xlabel,
  #                               breaks = majorBreaks,
  #                               minor_breaks = minorBreaks-.5,
  #                               labels = labels,
  #                               expand = c(0,0), limits = c((win+1)-.5, max(majorBreaks)+.5)) +
  #   ggplot2::scale_fill_manual("Critical Insability", breaks = c("Sig. DC level","Sig. CCP"), values = c("0"="white","Sig. DC level"="grey","5"="whitesmoke","Sig. CCP"="black")) +


  if(doPlot){
    # graphics::plot.new()
    graphics::plot(g)
  }

  return(invisible(g))
}




# Complex Networks# graph2svg <- function(TDM,pname){
#
#   # Create weighted Term-Term matrix
#   tTM <- as.matrix(TDM)
#   TTM <- tTM %*% t(tTM)
#   TTM <- log1p(TTM)
#
#   g <- graph.adjacency(TTM,weighted=T,mode="undirected",diag=F)
#   g <- simplify(g)
#
#   # Remove vertices that were used in the search query
#   Vrem <- which(igraph::V(g)$name %in% c("~dev~","~dys~","~sld~","development","children","dyslexia"))
#   g <- (g - igraph::V(g)$name[Vrem])
#
#   # Set colors and sizes for vertices
#   igraph::V(g)$degree <- igraph::degree(g)
#   rev         <- scaleRange(log1p(igraph::V(g)$degree))
#   rev[rev<=0.3]<-0.3
#
#   igraph::V(g)$color       <- grDevices::rgb(scaleRange(igraph::V(g)$degree), 1-scaleRange(igraph::V(g)$degree),  0, rev)
#   igraph::V(g)$size        <- 10*scaleRange(igraph::V(g)$degree)
#   igraph::V(g)$frame.color <- NA
#
#   # set vertex labels and their colors and sizes
#   igraph::V(g)$label       <- igraph::V(g)$name
#   igraph::V(g)$label.color <- grDevices::rgb(0, 0, 0, rev)
#   igraph::V(g)$label.cex   <- scaleRange(igraph::V(g)$degree)+.1
#
#   # set edge width and color
#   rew <- scaleRange(igraph::E(g)$weight)
#   rew[rew<=0.3]<-0.3
#
#   igraph::E(g)$width <- 2*scaleRange(igraph::E(g)$weight)
#   igraph::E(g)$color <- grDevices::rgb(.5, .5, 0, rew)
#   set.seed(958)
#
#   svg(paste(pname,sep=""),width=8,height=8)
#   graphics::plot(g, layout=layout.fruchterman.reingold(g))
#   dev.off()
#
#   return(g)
# }
#
# # Plot vertex neighbourhood
# hoodGraph2svg <- function(TDM,Vname,pname){
#
#   # Create weighted Term-Term matrix
#   tTM <- as.matrix(TDM)
#   TTM <- tTM %*% t(tTM)
#   TTM <- log1p(TTM)
#
#   ig <- graph.adjacency(TTM,weighted=T,mode="undirected",diag=F)
#   ig <- simplify(ig)
#
#   # Remove vertices that were used in the search query
#   Vrem <- which(igraph::V(ig)$name %in% c("~dev~","~dys~","~sld~","development","children","dyslexia"))
#   ig <- (ig - igraph::V(ig)$name[Vrem])
#
#   # This is a deletion specific for the Neighbourhood graphs
#   Vrem <- which(igraph::V(ig)$name %in% c("~rdsp~","~imp~","~som~","~bod~","~mlt~"))
#   ig   <- ig - igraph::V(ig)$name[Vrem]
#
#   idx <- which(igraph::V(ig)$name==Vname)
#   sg  <- graph.neighborhood(ig, order = 1, nodes=igraph::V(ig)[idx], mode = 'all')[[1]]
#
#   # set colors and sizes for vertices
#   igraph::V(sg)$igraph::degree <- igraph::degree(sg)
#
#   rev<-scaleRange(log1p(igraph::V(sg)$igraph::degree))
#   rev[rev<=0.3]<-0.3
#
#   igraph::V(sg)$color <- grDevices::rgb(scaleRange(igraph::V(sg)$igraph::degree), 1-scaleRange(log1p(igraph::V(sg)$igraph::degree*igraph::V(sg)$igraph::degree)),  0, rev)
#
#   igraph::V(sg)$size        <- 35*scaleRange(igraph::V(sg)$igraph::degree)
#   igraph::V(sg)$frame.color <- NA
#
#   # set vertex labels and their colors and sizes
#   igraph::V(sg)$label       <- igraph::V(sg)$name
#   igraph::V(sg)$label.color <- grDevices::rgb(0, 0, 0, rev)
#   igraph::V(sg)$label.cex   <- scaleRange(igraph::V(sg)$igraph::degree)
#
#   # set edge width and color
#   rew<-scaleRange(igraph::E(sg)$weight)
#   rew[rew<=0.3]<-0.3
#
#   igraph::E(sg)$width <- 6*scaleRange(igraph::E(sg)$weight)
#   igraph::E(sg)$color <- grDevices::rgb(.5, .5, 0, rew)
#
#   idV <- which(igraph::V(sg)$name==Vname)
#   idE <- incident(sg,igraph::V(sg)[[idV]])
#   igraph::E(sg)$color[idE] <- grDevices::rgb(0, 0, 1 ,0.8)
#
#   set.seed(958)
#
#   idx <- which(igraph::V(sg)$name==Vname)
#   svg(paste(pname,sep=""),width=8,height=8)
#   graphics::plot(sg,layout=layout.star(sg,center=igraph::V(sg)[idx]))
#   dev.off()
#
#   return(sg)
# }


#
# #' PSDslope
# #'
# #' @param y    A time series object, or a vector that can be converted to a time series object.
# #' @param fs    Sample frequency (defults to 1).
# #' @param nfft    Number of frequencies to estimate (defaults to next power of 2)
# #' @param fitRange    Vector of length 2 with range of frequencies to perform log-log fit.
# #' @param doPlot    Plot the log-log spectrum and slope.
# #'
# #' @return
# #' @export
# #'
# #' @examples
# #'
# PSDslope <- function(y  = stats::ts(rnorm(n = 1024), frequency = 1),
#                      fs = stats::frequency(y),
#                      nfft = 2^(nextpow2(length(y)/2)),
#                      fitRange = c(1,round(.1*nfft)),
#                      doPlot = FALSE){
#   require(oce)
#   require(signal)
#   if(!stats::is.ts(y)){stats::ts(y, frequency = fs)}
#
#   win <- signal::hamming(n=nfft)
#
#   perioGram <- oce::pwelch(x = y, window = win, fs = stats::frequency(y), nfft = nfft, doPlot = FALSE)
#   spec <- data.frame(Frequency = perioGram$freq, Power = perioGram$spec)
#   spec[1,1:2] <- NA
#   fit <- stats::lm(log10(spec$Power[fitRange[1]:fitRange[2]])~log10(spec$Power[fitRange[1]:fitRange[2]]))
#   return(list(spec = spec,
#               slope = fit)
#   )
# }



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


# Time Series HELPERS ----

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
#' set.seed(4321)
#' y1 <- rnorm(5000)
#' y2 <- y1-(mean(y1)+rnorm(1))
#'
#' ts_permtest_block(y1 = y1, y2 = y2)
#'
ts_permtest_block <- function(y1, y2 = NULL, targetValue = 0, Nperms = 19, sim = "geom", l = 3, alpha = .05, returnBootObject = FALSE){

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
ts_changeindex <- function(y, returnRectdata=TRUE, groupVar = NULL, labelVar = NULL, discretize=FALSE, nbins = 5){


  if(is.discrete(y)){
    y <- as.numeric_discrete(y)
    #y <- ts_checkfix(y,checkNumericVector = TRUE, fixNumericVector = TRUE, checkWholeNumbers = TRUE, fixWholeNumbers = TRUE)
  }

  if(!is.null(names(y))&is.null(labelVar)){
    labelVar <- names(y)
  }

  if(length(unique(y))>10){
    warning("More than 10 epochs detected!")
  }

  xmax <- c(which(diff(y)!=0),NROW(y))
  xmin <- c(1,xmax)[1:NROW(xmax)]
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
    return(c(1,xmin))
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
                          y.name  = names(y[1]),
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

  id <- ifelse(is.null(colnames(y)),ifelse(is.null(names(y)),deparse(substitute(y)),names(y)[1]),colnames(y)[1])
  y.ori <- y
  N  <- NROW(y)

  if((emDim-1) * emLag > N){stop(paste0("Time series length (N = ",N,") is too short to embed in ",emDim," dimensions with delay ",emLag))}

  if(!is.null(dim(y))){
    y <- y_data <- as.numeric(y[,1])
    if(!silent){cat("\ntaking first column...\n")}
  } else {
    y <- y_data <- as.numeric(y)
  }

  if(any(stats::is.ts(y), zoo::is.zoo(y), xts::is.xts(y))){
    y <- stats::time(y)
    emTime <- lubridate::as_datetime(y[emLag+1])- lubridate::as_datetime(y[1])
  } else {
    y <- zoo::index(y)
    emTime <- emLag
  }



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
    emY <-  y
  }

  # Alternative: rollapply(y, list(-d * seq(0, k-1)), c)


  attr(emY, "embedding.dims") <- emDim
  attr(emY, "embedding.lag")  <- emLag
  attr(emY, "embedding.time") <- emTime
  attr(emY, "variable.y")     <- id

  if(returnOnlyIndices){
    attr(emY, "data.y") <- y_data
    return(as.matrix(emY))
  } else {
    if(emDim>1){
      for(c in 1:NCOL(emY)){
        emY[,c] <-  y_data[emY[,c]]
      }
    } else{
      emY <- as.matrix(y_data)
    }
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
#' # Get an interesting numeric vector from package DescTools
#' y <- DescTools::Fibonacci(1:26)
#'
#' # Return the first order derivative as a vector
#' ts_diff(y=y,addColumns=FALSE)
#'
#' # Return original and derivative as a matrix
#' plot(stats::ts(ts_diff(y=y, addColumns=TRUE)))
#'
#' # Works on multivariate data objects with mixed variable types
#' df <- data.frame(x=letters, y=1:26, z=sin(y))
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
      maskWhich[[o]] <- c(1:o,(NROW(dy)-o):NROW(dy))
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
#'
#' @return If `adjustY = FALSE`, or, a list object with fields that represent arguments of \link[casnet]{ts_trimfill}, then the (adjusted) vector `y` is returned with an attribute `"windower"`. This is a list object with fields that contain the indices for each window that fits on `y`, given `win`, `step` or `overlap` and the settings of `adjustY`. If `adjustY = NA`, only the list object is returned.
#' @export
#'
#' @family Time series operations
#'
ts_windower <- function(y, win=length(y), step=round(win/2), overlap=NA, adjustY=NA){

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

  if(!is.na(overlap)){step <- floor(overlap*win)}

  wIndices <- list()
  wIndex   <- seq(1,(NROW(y)-win),step)
  iDiff    <- (dplyr::last(wIndex) + win-1)-NROW(y)

  wIndices <- plyr::llply(wIndex,function(i){i:min(i+win-1,length(y))})

  if(adjustOK){
    if(iDiff<0){
      if(adjustY$action%in%"fill"){
        wIndices[[length(wIndices)]] <- c(wIndices[[length(wIndices)]],seq(max(wIndices[[length(wIndices)]]),max(wIndices[[length(wIndices)]])+abs(iDiff)))
      }
    }
  }
  names(wIndex) <- paste("stepsize",floor(step*win),"| window",1:length(wIndex))
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
#'  Use recursive partitioning function (\link[rpart]{rpart} to perform a 'classification' of relatively stable levels in a timeseries.
#'
#' @param y A time series of numeric vector
#' @param minDataSplit An integer indicating how many datapoints should be in a segment before it will be analysed for presence of a level change (default = `12`)
#' @param minLevelDuration Minimum duration (number of datapoint) of a level (default = `round(minDataSplit/3)`)
#' @param changeSensitivity A number indicating a criterion of change that must occur before declaring a new level. Higher numbers indicate higher levels of change must occur before a new level is considered. For example, if `method = "anova"`, the overall `R^2` after a level is introduced must increase by the value of `changeSensitivity`, see the `cp` parameter in [rpart::rpart.control].     (default = `0.01`)
#' @param maxLevels Maximum number of levels in one series (default = `30`)
#' @param method The partitioning method to use, see the manual pages of [rpart] for details.
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
#' wn <- ts_levels(rnorm(100))
#' plot(wn$pred$x,wn$pred$y, type = "l")
#' lines(wn$pred$p, col = "red3", lwd = 2)
#'
#' # This is due to the default changeSensitivity of 0.01
#'
#' lines(ts_levels(rnorm(100),changeSensitivity = .1)$pred$p, col = "steelblue", lwd = 2)
#'
#'
ts_levels <- function(y, minDataSplit=12, minLevelDuration=round(minDataSplit/3), changeSensitivity = 0.01, maxLevels=30, method=c("anova","poisson","class","exp")[1]){
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

  return(list(tree  = tree,
              pred  = dfs))
}


#' Find Peaks or Wells
#'
#' @param y A time series or numeric vector
#' @param window Window in whcih to look for peaks or wells
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

#' Slice columns of a matrix in epochs
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
#'
ts_slice<-function(y,epochSz=4){
  if(!is.matrix(y)){yy <- as.matrix(y)} else {yy<-y}
  N<-dim(yy)
  wIndex <- plyr::llply(seq(1,N[1],epochSz),function(i) yy[i:min(i+epochSz-1,N[1]),1:N[2]])
  delID <- which(lengths(wIndex)!=epochSz)%00%NA
  for(del in delID){
    if(!is.na(delID)){wIndex <- wIndex[-del]}
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
#' @param y A 1D time series object or numeric verctor.
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
#' y <- rnorm(100)
#' plot(ts(y))
#'
#' y_img <- ts_rasterize(y)
#' image(y_img,col=c("white","black"))
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
#' @param nbins The number of bins used to transform a continuous time series, or, the number of expected (given `nbins`, or, theoretically possible) values for a discrete series (default = `unique(yd)`)
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
#' # Not specifying the expected bins will geberate a different matrix!
#' ts_transmat(yd = yd, nbins = length(unique(yd)))
#'
ts_transmat <- function(yd, nbins = unique(yd)){

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

# HELPERS ----

#' Get some nice colours
#'
#' @param pal The colour pallette, one of `"pe"`,`"mm"`,`"le"``"an"` (default = `"pe"`)
#' @param Ncols Number of colours
#'
#' @return A list of colours
#' @export
#'
#' @examples
#'
#' getColours(Ncol=5)
#'
getColours <- function(pal = c("pe","mm","le","an")[1],Ncols = 20){
  pl <- c("#A65141","#ABB2C4","#C89284","#7896BC","#E1CAC2","#536489","#ECE3CD","#575E7A","#BEBEC4","#C4AD75","#30261E","#BC964D","#3D434D","#D29E55","#B4B6A7","#9B7342","#E8D69D","#48211A","#EBDAB4","#3D4B68","#E8DCCF","#3E688C","#9D3D2D","#507074","#9A756C","#546959","#93725F","#858563","#E1D3B1","#A6B4BB","#78828C","#AA9C6A","#91A193","#CDB16E","#1B528D","#B8A98F","#6E432D","#4F5B71","#D9B196","#20304F","#827561","#98A39E","#8B4F31","#7E7B65","#C1B6A1","#775E45","#C0B3A2","#5A524A","#BBA27F","#3A3935","#C9BDA3","#626961","#8A4F42","#8D8A77","#947B5E","#5D3C30","#AA8470","#493A2C")

  if(Ncols<=58){
    cols <- pl[1:Ncols]
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
#  Attribution-Share Alike 2.5 Generic Creative Commons license:
#  http://creativecommons.org/licenses/by-sa/2.5/
#  If you use this work, please cite:
#  Little MA et al. (2007), "Exploiting nonlinear recurrence and fractal
#  scaling properties for voice disorder detection", Biomed Eng Online, 6:23
#
#  As of 20080323 markup
#  If you use this work, consider saying hi on comp.dsp
#  Dale B. Dalrymple
#'
noise_powerlaw <- function(y = NULL, alpha=-1, N=100, standardise = FALSE, randomPower = FALSE, seed = NA){

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
  N2 = floor(N/2)-1
  f = (2:(N2+1))
  A2 = 1/(f^(alpha/2))

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
noise_fGn <- function(H=0.5, N = 100, mu = NULL, sigma = NULL){

  # Determine whether fGn or fBn should be produced.
  if(H%)[%c(0,1)){stop("H must be in (0,1] for fGn!")}
  fBn = 0

  # Calculate the fGn.
  if (H == 0.5){
    y <- stats::rnorm(N)  # If H=0.5, then fGn is equivalent to white Gaussian noise.
    } else {
    # If this function was already in memory before being called this time,
    # AND the values for N and H are the same as the last time it was
    # called, then the following (persistent) variables do not need to be
    # recalculated.  This was done to improve the speed of this function,
    # especially when many samples of a single fGn (or fBn) process are
    # needed by the calling function.
    # persistent Zmag Nfft Nlast Hlast
    #    if isempty(Zmag) | isempty(Nfft) | isempty(Nlast) |isempty(Hlast) | N ~= Nlast #| H ~= Hlast
    # The persistent variables must be (re-)calculated.
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
noise_fBm <- function(H=1.5, N = 100, mu = NULL, sigma = NULL){

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


#' Numeric factor to numeric vector
#'
#' Converts a factor with numeric levels to a numeric vector, using the values of the levels.
#'
#' @param x A factor based on numeric values.
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
as.numeric_factor <- function(x, keepNA = FALSE){
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
#' Converts a factor with numeric levels, or, a character vector with numeric values to a numeric vector using [as.numeric_factor], or, [as.numeric_character] respectively. If an unnamed numeric vector is passed, it will be returned as a named numeric vector.
#'
#' @param x A factor with levels that are numeric, or, a character vector representing numbers.
#' @param keepNA Keep NA values (`TRUE`), or remove them (default = `FALSE`)
#'
#' @return A numeric vector with factor levels / numeric character values as names.
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
#'
as.numeric_discrete <- function(x, keepNA = FALSE){

  if(plyr::is.discrete(x)){
    if(is.factor(x)){
      if(suppressWarnings(all(is.na(as.numeric(levels(x)))))){
        x <- as.character(x)
      } else {
        y <- as.numeric_factor(x, keepNA = keepNA)
      }
    }
    if(is.character(x)){
      y <- as.numeric_character(x, keepNA = keepNA)
    }
  } else {
    if(is.numeric(x)){
      if(is.null(names(x))){
        names(x) <- paste(x)
        y <- x
      } else {
        stop("Variable is not a factor, character vector, or, unnamed numeric vector.")
      }
    }
  }
  return(y)
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
#' # Values < mn will return < lo (default=0)
#' # Values > mx will return > hi (default=1)
#' elascer(somenumbers,mn=-1,mx=99)
#'
#' elascer(somenumbers,lo=-1,hi=1)
#' elascer(somenumbers,lo=-1,hi=1, groupwise = TRUE)
#'
#' elascer(somenumbers,mn=-10,mx=100,lo=-1,hi=4)
#' elascer(somenumbers,mn= NA,mx=100,lo=-1,hi=4, groupwise = TRUE)
#'
elascer <- function(x,mn=NA,mx=NA,lo=0,hi=1,groupwise = FALSE, keepNA = TRUE){

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
