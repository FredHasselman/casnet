# package casnet ----
#
# Recurrence Plots ----
#
# rp functions

#' Create a Distance Matrix
#'
#' @param y1 A numeric vector or time series
#' @param y2 A numeric vector or time series for cross recurrence
#' @param emDim The embedding dimensions
#' @param emLag The embedding lag
#' @param emRad The threshold (emRad) to apply to the distance matrix to create a binary or weighted matrix. If `NULL`, an unthresholded matrix will be created (default = `NULL`)
#' @param theiler Use a `theiler` window around the main diagonal (Line of Identity/Synchronisation) to remove auto-correlations at short time-lags:
#' * `0` will include the main diagonal in all RQA measure calculations.
#' * `1` will remove the main diagonal from all RQA measure calculations.
#' * `NA` (default), will check if the matrix is symmetrical , if so, it will remove the diagonal by setting `theiler = 0` (Line of Identity, Auto-RQA), if it is not symmetrical (Line of Synchronisation, Cross-RQA) it will set `theiler = 1`.
#' * A value greater than `1` will remove that many diagonals around and including the diagonal from all RQA measure calculations. So `theiler = 2` means exclude `2` diagonals around the main diagonal, including the main diagonal itself: `[-1,0,1]`.
#' If `theiler` is a numeric vector of `length(theiler) == 2` it is possible to exclude an asymmetrical window. The values are interpreted as end points in a sequence of diagonal ID's, e.g. `theiler = c(-1,5)` will exclude `[-1,0,1,2,3,4,5]`. If `length(theiler) > 2`, the values will be considered individual diagonal ID's, e.g. `theiler = c(-3,-1,0,2,5)`, will exclude only those specific ID's. Also see the note.
#' @param includeDiagonal Use to force inclusion of the diagonal when calculating RQA measures. The default `NA` behaves the same way as parameter `theiler`. If set to `TRUE` it will include the diagonal even in the case of Auto-RQA (default = `NA`)
#' @param to.ts Should `y1` and `y2` be converted to time series objects?
#' @param order.by If `to.ts = TRUE`, pass a vector of the same length as `y1` and `y2`. It will be used as the time index, if `NA` the vector indices will be used to represent time.
#' @param to.sparse Should sparse matrices be used?
#' @param weighted If `FALSE` a binary matrix will be returned. If `TRUE` every value larger than `emRad` will be `0`, but values smaller than `emRad` will be retained (default = `FALSE`)
#' @param weightedBy After setting values smaller than `emRad` to `0`, what should the recurrent values represent? The default is to use the state space similarity (distance/proximity) values as weights (`"si"`). Other option are `"rt"` for *recurrence time* and `"rf"` for *recurrence time frequency* (default = `"si"`)
#' @param method Distance measure to use. Any option that is valid for argument `method` of [proxy::dist()]. Type `proxy::pr_DB$get_entries()` to see a list of all the options. Common methods are: "Euclidean", "Manhattan", "Minkowski", "Chebysev" (or the same but shorter: "L2","L1","Lp" and "max" distance) (default = `"Euclidean"`)
#' @param rescaleDist Should the distance matrix be rescaled? Options are "none", "maxDist" to create a unit scale, "meanScale" to creat z-scores based on the mean distance. (default = `"none"`)
#' @param targetValue A value passed to `est_radius(...,type="fixed", targetMeasure="RR")` if `is.na(emRad)==TRUE`.
#' @param chromatic Perform a chromatic RQA. This assumes the recurring values represent the labels of an unordered categorical variable (default = `FALSE`)
#' @param returnMeasures Should the output of [rp_measures()] be returned as an attribute `"measures"` to the matrix? If `silent = FALSE` results will also be output to the console. (default = `FALSE`)
#' @param doPlot Plot the matrix by calling [rp_plot()] with default settings
#' @param doEmbed If `FALSE`, a distance matrix will be returned that is not embedded by `emDim` and `emLag` (Multidimensional RQA). If `y1` and/or `y2` are data frames, the columns will be used as the state space dimensions (default = `TRUE`)
#' @param silent Silent-ish mode
#' @param ... Any parameters to pass to [rp_plot()] if `doPlot = TRUE`
#'
#' @note The calculation of the (C)RQA measures in [casnet] can be different from other packages. For example, depending on the value of `theiler` the main diagonal can be included or excluded from the calculations, whereas some software will always include the diagonal.
#'
#' @return A (Coss-) Recurrence matrix with attributes:
#'
#' * `emdims1` and `emdims2` - A matrix of surrogate dimensions
#' * `emdims1.name` and `emdims2.name` - Names of surrogate dimensions
#' * `method` and `call` - The distance `method` used by [proxy::dist()]
#' * `weighted` - Whether a weighted matrix is returned
#' * `emDim`, `emLag` and `emRad` - The embedding parameters
#' * `AUTO` - Whether the matrix represents AUTO recurrence
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
               theiler = NA,
               includeDiagonal = NA,
               to.ts = NULL,
               order.by = NULL,
               to.sparse = TRUE,
               weighted = FALSE,
               weightedBy = "si",
               method = "Euclidean",
               rescaleDist = c("none","maxDist","meanDist")[1],
               targetValue  = .05,
               chromatic = FALSE,
               returnMeasures = FALSE,
               doPlot = FALSE,
               doEmbed = TRUE,
               silent = TRUE,
               ...){


  if(is.null(y2)){
    y2 <- y1
    attributes(y2) <- attributes(y1)
  }

  atlist <- attributes(y1)
  if(any(names(atlist) %in% c("emDims1","emDims2","emDims1.name","emDims2.name"))){
    stop("Input is already a recurrence matrix created by 'rp()'. To manually create a binary matrix from a distance matrix use 'mat_di2bi()'. To create a weighted matrix use 'mat_di2we()'")
  }


  if(chromatic){
    if(any(is.data.frame(y1),is.data.frame(y2))){
      stop("Multidimensional Chromatic RQA has not been invented (yet)!")
    }
    uniqueID <- as.numeric_discrete(x = c(y1,y2), sortUnique = TRUE, keepNA = TRUE)
    y1  <- uniqueID[1:length(y1)]
    y2  <- uniqueID[(length(y2)+1):length(uniqueID)]

    chromaDims <- list(y1 = uniqueID[1:length(y1)],
                       y2 = uniqueID[(length(y2)+1):length(uniqueID)])

    chromaNames <- unique(uniqueID)
    nList <- list()
    for(n in seq_along(chromaNames)){
      nList[[n]] <- names(uniqueID)[which(uniqueID%in%chromaNames[n])]
    }
    names(chromaNames) <- unique(unlist(nList))

  } else {
    chromaNames <- NA
    chromaDims  <- NA
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

  if(any(!doEmbed,chromatic)){
    emDim <- 1
    emLag <- 0
    if(chromatic){
      emRad <- 0
    }
  }

  et1 <- ts_embed(y1, emDim=emDim, emLag=emLag, silent = silent)
  et2 <- ts_embed(y2, emDim=emDim, emLag=emLag, silent = silent)

  dist_method <- return_error(proxy::pr_DB$get_entry(method))
  if("error"%in%class(dist_method$value)){
    stop("Unknown distance metric!\nUse proxy::pr_DB$get_entries() to see a list of valid options.")
  } else {
    dmat <- proxy::dist(x = et1,
                        y = et2,
                        method = method,
                        diag = TRUE)
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

  if(!chromatic){
    if(rescaleDist=="maxDist"){dmat <- dmat/max(dmat,na.rm = TRUE)}
    if(rescaleDist=="meanDist"){dmat <- (dmat-mean(dmat, na.rm=TRUE))/sd(dmat, na.rm = TRUE)}
  }

  if(!is.null(emRad)){
    if(is.na(emRad)){
      suppressMessages(emRad <- est_radius(RM = dmat, emDim = emDim, emLag = emLag, targetValue = targetValue, silent = TRUE))
      if(emRad$Converged){
        emRad <- emRad$Radius
      } else {
        emRad <- est_radius(RM = dmat, emDim = emDim, emLag = emLag, targetValue = targetValue, tol = 0.2, silent = TRUE, radiusOnFail = "percentile")$Radius
      }
    }
    if(chromatic){
      dmat <- mat_di2ch(dmat, y = et1, emRad = emRad, convMat = to.sparse)
    } else {
      if(weighted){
        dmat <- mat_di2we(dmat, emRad = emRad, convMat = to.sparse)
      } else {
        dmat <- mat_di2bi(dmat, emRad = emRad, convMat = to.sparse)
      }
    }
  }

  dmat <- rp_checkfix(dmat, checkAUTO = TRUE, fixAUTO = TRUE, checkSPARSE = TRUE, fixSPARSE = TRUE, checkS4 = TRUE, fixS4 = TRUE)
  suppressMessages(dmat <- setTheiler(RM = dmat, theiler = theiler, chromatic = chromatic))

  if(returnMeasures){
    rpOut <-  rp_measures(RM = dmat,
                          emRad = emRad%00%NA,
                          silent = silent,
                          theiler = theiler,
                          includeDiagonal = includeDiagonal,
                          chromatic = chromatic)
  } else {
    rpOut <- NA
  }

  if(to.sparse){
    attributes(dmat)$emDims1  <- et1
    attributes(dmat)$emDims2  <- et2
    attributes(dmat)$emDims1.name <- colnames(y1)
    attributes(dmat)$emDims2.name <- colnames(y2)
    attributes(dmat)$embedded <- doEmbed
    attributes(dmat)$emLag <- emLag
    attributes(dmat)$emDim <- emDim
    attributes(dmat)$emRad <- emRad%00%NA
    attributes(dmat)$measures <- rpOut
    attributes(dmat)$weighted <- weighted
    attributes(dmat)$weightedBy <- weightedBy
    attributes(dmat)$chromatic <- chromatic
    attributes(dmat)$chromaNames <- chromaNames
    attributes(dmat)$chromaDims <- chromaDims
  } else {
    attr(dmat,"emDims1") <- et1
    attr(dmat,"emDims2") <- et2
    attr(dmat,"emDims1.name") <- colnames(y1)
    attr(dmat,"emDims2.name") <- colnames(y2)
    attr(dmat,"weighted") <- weighted
    attr(dmat,"embedded") <- doEmbed
    attr(dmat,"emLag") <- emLag
    attr(dmat,"emDim") <- emDim
    attr(dmat,"emRad") <- emRad%00%NA
    attr(dmat,"measures") <- rpOut
    attr(dmat,"weighted") <- weighted
    attr(dmat,"weightedBy") <- weightedBy
    attr(dmat,"chromatic") <- chromatic
    attr(dmat,"chromaNames") <- chromaNames
    attr(dmat,"chromaDims") <- chromaDims
  }

  if(is.null(attributes(dmat))){stop("lost attributes")}

  if(doPlot){

    rp_plot(RM = dmat)

    # if(...length()>0){
    #   dotArgs <- list(...)
    #   nameOK  <- names(dotArgs)%in%methods::formalArgs(rp_plot)
    #   dotArgs <- dotArgs[nameOk]
    # } else {
    #   # Plot with defaults
    #   dotArgs  <- formals(rp_plot)
    #   #nameOk   <- rep(TRUE,length(dotArgs))
    # }
    #
    # if(!"RM"%in%names(dotArgs)|length(dotArgs$RM)==1){
    #   dotArgs$RM <- dmat
    # }
    # if(!is.na(emRad%00%NA)&(!"radiusValue"%in%names(dotArgs))){
    #   dotArgs$radiusValue <- emRad
    # }
    # do.call(rp_plot, as.list(dotArgs))
  }

  return(dmat)
}


#' Get (C)RQA measures based on a binary matrix
#'
#' A zoo of measures based on singular recurrent points, diagonal, vertical and horizontal line structures (anisotropic) will be caluclated.
#'
#' @aliases crqa_rp
#' @inheritParams rp
#' @param RM A distance matrix (set `emRad = NA` to estimate a radius), or a matrix of zeroes and ones
#' @param emRad Threshold for distance value that counts as a recurrence (ignored is `RM` is a binary matrix)
#' @param DLmin Minimal diagonal line length (default = `2`)
#' @param VLmin Minimal vertical line length (default = `2`)
#' @param HLmin Minimal horizontal line length (default = `2`)
#' @param DLmax Maximal diagonal line length (default = length of diagonal -1)
#' @param VLmax Maximal vertical line length (default = length of diagonal -1)
#' @param HLmax Maximal horizontal line length (default = length of diagonal -1)
#' @param AUTO Auto-recurrence? (default = `FALSE`)
#' @param includeDiagonal Should the diagonal be included in the calculation of the Recurrence Rate? If `NA` the value will be decided by the symmetry of the matrix, the diagonal will be removed for Auto RQA (`AUTO = TRUE`) but not for Cross RQA (`AUTO = FALSE`)  (default = `NA`)
#' @param chromatic Force chromatic RQA? If `NA` the value of the `RM` attribute `"chromatic"` will be used, if present (default = `NA`)
#' @param anisotropyHV Return anisotropy ratio measures based on Horizontal and Vertical lines. The ratios are calculated as `(horizontal - vertical) / (horizontal + vertical)`. So a value of 0 means no anisotropy, negative ratios indicate the measures based on vertical lines had  higher values, positive ratios indicate the measures based on horizontal lines had higher values  (default = `FALSE`)
#' @param asymmetryUL Return asymmetry ratio measures based on Upper and Lower triangles. The ratios are calculated as `(upper - lower) / (upper + lower)`. So a value of 0 means no asymmetry, negative ratios indicate the measures based on the lower triangle had the higher values, positive ratios indicate measures based on the upper triangle had higher values (default = `FALSE`)
#' @param recurrenceTimes Return measures based on 'white lines', the recurrence times (default = `FALSE`)
#' @param matrices Return matrices? (default = `FALSE`)
#' @param Nboot How many bootstrap replications? (default = `NULL`)
#' @param CL Confidence limit for bootstrap results (default = `.95`)
#' @param targetValue A value passed to `est_radius(...,type="fixed", targetMeasure="RR", tol = .2)` if `is.na(emRad)==TRUE`, it will estimate a radius (default = `.05`).
#' @param doParallel Speed up calculations by using the parallel processing options provided by `parallel` to assign a separate process/core for each window in windowed (C)RQA analysis and [parallel::detectCores()] with  `logical = TRUE` to decide on the available cores (default = `FALSE`)
#' @param silent Do not display any messages (default = `TRUE`)

#' @return A list object containing (C)RQA measures (and matrices if requested)
#'
#' @export
#'
#' @family Recurrence Quantification Analysis
#'
#'
rp_measures <- function(RM,
                        emRad = NA,
                        DLmin = 2,
                        VLmin = 2,
                        HLmin = 2,
                        DLmax = length(Matrix::diag(RM)),
                        VLmax = length(Matrix::diag(RM)),
                        HLmax = length(Matrix::diag(RM)),
                        AUTO      = NULL,
                        theiler   = NA,
                        includeDiagonal = NA,
                        chromatic = NA,
                        anisotropyHV = FALSE,
                        asymmetryUL = FALSE,
                        recurrenceTimes = FALSE,
                        matrices  = FALSE,
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

  if(is.na(chromatic)){
    if(!is.null(attr(RM,"chromatic"))){
      chromatic <- attr(RM,"chromatic")
    } else {
      chromatic <- FALSE
    }
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

      emRad <- est_radius(RM = RM, emDim = emDim, emLag = emLag, targetValue = targetValue, tol = .2, radiusOnFail = "percentile", silent = silent)
      if(emRad$Converged){
        emRad <- emRad$Radius
      } else {
        emRad <- stats::sd(RM,na.rm = TRUE)
      }
    } # not emRad attr
  }

  #uval <- unique(as.vector(RM))
  if(!all(as.vector(RM)==0|as.vector(RM)==1)){

    if(chromatic){

      # prepare data
      if(any(grepl("Matrix",class(RM)))){
        RM     <- rp_checkfix(RM, checkTSPARSE = TRUE, fixTSPARSE = TRUE)
        #meltRP <- data.frame(Var1 = (RM@i+1), Var2 = (RM@j+1), value = as.numeric(RM@x))
        meltRP <-  mat_mat2ind(Matrix::as.matrix(RM))
      } else {
        meltRP <-  mat_mat2ind(as.matrix(RM))
      }

      if(!is.null(attr(RM,"chromaNames"))){
        chromaNumbers <- attr(RM,"chromaNames")
        Cind <- sort(chromaNumbers, index.return = TRUE)
        chromaNames <- names(chromaNumbers)[Cind$ix]
        chromaNumbers <- chromaNumbers[Cind$ix]
      } else {
        chromas <- meltRP %>% dplyr::filter(.data$value!=0)
        chromaNumbers <- as.numeric_discrete(sort(unique(chromas$value)))
      }

      chromaNumbers <- stats::na.exclude(chromaNumbers)

      rm(meltRP)

      if(NROW(chromaNumbers)>=(NROW(RM)/2)){
        warning(paste("Chromatic RQA will have a large number of categories:",NROW(chromaNumbers)))
      }

      chromaList <- list()
      matrixList <- list()

      for(i in seq_along(chromaNumbers)){

        RMtmp <- RM
        RMtmp[RM!=i] <- 0
        RMtmp[RM==i] <- 1

        suppressMessages(tmpOut <- rp_calc(RMtmp,
                                           emRad = emRad,
                                           DLmin = DLmin,
                                           VLmin = VLmin,
                                           HLmin = HLmin,
                                           DLmax = DLmax,
                                           VLmax = VLmax,
                                           HLmax = HLmax,
                                           theiler = theiler,
                                           AUTO  = AUTO,
                                           includeDiagonal = includeDiagonal,
                                           chromatic = chromatic,
                                           anisotropyHV = anisotropyHV,
                                           asymmetryUL = asymmetryUL,
                                           recurrenceTimes = recurrenceTimes,
                                           matrices  = matrices)
        )

        if(matrices){
          chromaList[[i]] <- tmpOut$crqaMeasures
          matrixList[[i]] <- tmpOut$crqaMatrices
        } else {
          chromaList[[i]] <- tmpOut
        }
        rm(RMtmp, tmpOut)
      }

      names(chromaList) <- names(chromaNumbers)
      if(matrices){
        names(matrixList) <- names(chromaNumbers)
        out <- list(crqaMeasures = plyr::ldply(chromaList, .id = "chroma"),
                    crqaMatrices =  matrixList)
      } else {
        out <- plyr::ldply(chromaList, .id = "chroma")
      }

    } else {
      if(!is.null(emRad)){
        RM <- mat_di2bi(RM,emRad)
      } else {
        if(!chromatic){
          stop("Expecting a binary (0,1) matrix.\nUse 'est_radius()', or set 'chromatic = TRUE'")
        }
      }
    }
  }
  #rm(RM)

  if(!chromatic){
    suppressMessages(out <- rp_calc(RM,
                                    emRad = emRad,
                                    DLmin = DLmin,
                                    VLmin = VLmin,
                                    HLmin = HLmin,
                                    DLmax = DLmax,
                                    VLmax = VLmax,
                                    HLmax = HLmax,
                                    theiler = theiler,
                                    AUTO  = AUTO,
                                    includeDiagonal = includeDiagonal,
                                    chromatic = chromatic,
                                    anisotropyHV = anisotropyHV,
                                    asymmetryUL = asymmetryUL,
                                    recurrenceTimes = recurrenceTimes,
                                    matrices  = matrices))
    }

    if(matrices){
      tab <- out$crqaMeasures
    } else {
      tab <- out
    }

  if(chromatic){
    Nrows <- length(chromaNames)
  } else {
    Nrows <- 1
  }

    outTable <- list(`Global Measures` = data.frame(Global = "Matrix",
                                                    `Max points` = tab$RP_max,
                                                    `N points` = tab$RP_N,
                                                    RR = tab$RR,
                                                    Singular = tab$SING_N%00%NA,
                                                    Divergence = tab$DIV_dl,
                                                    Repetitiveness = tab$REP_av),
                     `Line-based Measures` = data.frame(`Lines`    = rep(c("Diagonal", "Vertical", "Horizontal", "V+H Total"), each = Nrows),
                                                        `N lines`  = c(tab$N_dl,tab$N_vl, tab$N_hl, tab$N_hv),
                                                        `N points` = c(tab$N_dlp%00%NA,tab$N_vlp%00%NA,tab$N_hlp%00%NA, tab$N_hvp%00%NA),
                                                        `Measure`  = rep(c("DET","V LAM","H LAM", "V+H LAM"), each = Nrows),
                                                        `Rate`     = c(tab$DET, tab$LAM_vl, tab$LAM_hl, tab$LAM_hv),
                                                        `Mean`     = c(tab$MEAN_dl, tab$TT_vl, tab$TT_vl, tab$TT_hv),
                                                        `Max.`     = c(tab$MAX_dl, tab$MAX_vl, tab$MAX_hl, tab$MAX_hv),
                                                        `ENT`      = c(tab$ENT_dl, tab$ENT_vl, tab$ENT_hl, tab$ENT_hv),
                                                        `ENT_rel` = c(tab$ENTrel_dl, tab$ENTrel_vl, tab$ENTrel_hl, tab$ENTrel_hv),
                                                        `CoV`   = c(tab$CoV_dl, tab$CoV_vl, tab$CoV_hl, tab$CoV_hv)))

    if(chromatic){
      rownames(outTable$`Global Measures`) <- chromaNames
      rownames(outTable$`Line-based Measures`) <- paste(1:NROW(outTable$`Line-based Measures`), rep(chromaNames, 4))
    }

    if(anisotropyHV){

      outTable$`Horizontal/Vertical line anisotropy` <-  data.frame(`Ratio` = "H/V line measures",
                                                                    `N lines`  = tab$Nlines_ani,
                                                                    `N points` = tab$N_hlp%00%NA/tab$N_vlp%00%NA,
                                                                    `Measure` = "LAM",
                                                                    `Rate`    = tab$LAM_ani,
                                                                    `Mean`    = tab$MEAN_hvl_ani,
                                                                    `Max`     = tab$MAX_hvl_ani,
                                                                    `ENT`     = tab$ENT_hvl_ani)
      if(chromatic){
      rownames(outTable$`Horizontal/Vertical line anisotropy`) <- chromaNames
      }
    }

    if(asymmetryUL){

      outTable$`Upper/Lower triangle asymmetry` <- list(
        `Global Measures` =  data.frame(`Global Ratio` = "U/L of points",
                                        `N points` = tab$ratios.ul.Npoints_ul_ani,
                                        RR = tab$ratios.ul.RR_ul_ani,
                                        Singular = tab$ratios.ul.SING_N_ul_ani,
                                        Divergence = tab$ratios.ul.DIV_ul_ani,
                                        Repetetiveness = tab$ratios.ul.REP_ul_ani),
        `Line-based Measures` = data.frame(
          `Line ratio` = rep(c("D lines", "V lines", "H lines"), each = Nrows),
          `N lines`  = c(tab$ratios.ul.NDlines_ul_ani,
                         tab$ratios.ul.NVlines_ul_ani,
                         tab$ratios.ul.NHlines_ul_ani),
          `N points` = c(tab$upper.tri.N_dlp/tab$lower.tri.N_dlp,
                                  tab$upper.tri.N_vlp/tab$lower.tri.N_vlp,
                                  tab$upper.tri.N_hlp/tab$lower.tri.N_hlp),
          `Measure` = rep(c("DET","V LAM", "H LAM"), each = Nrows),
          `Rate`    = c(tab$ratios.ul.DET_ul_ani, tab$ratios.ul.LAM_vl_ul_ani, tab$ratios.ul.LAM_hl_ul_ani),
          `Mean`    = c(tab$ratios.ul.MEAN_dl_ul_ani, tab$ratios.ul.MEAN_vl_ul_ani, tab$ratios.ul.MEAN_hl_ul_ani),
          `Max`     = c(tab$ratios.ul.MAX_dl_ul_ani, tab$ratios.ul.MAX_vl_ul_ani, tab$ratios.ul.MAX_hl_ul_ani),
          `ENT`     = c(tab$ratios.ul.ENT_dl_ul_ani,
                                   tab$ratios.ul.ENT_vl_ul_ani,
                                   tab$ratios.ul.ENT_hl_ul_ani))
      )

      if(chromatic){
        id <- which(names(outTable)%in%"Upper/Lower triangle asymmetry")
        rownames(outTable[[id]]$`Global Measures`) <- chromaNames
        rownames(outTable[[id]]$`Line-based Measures`) <- paste(1:NROW(outTable[[id]]$`Line-based Measures`), rep(chromaNames, 3))
      }
    }


    if(!silent){
      cat("\n~~~o~~o~~casnet~~o~~o~~~\n")
      if(chromatic){
        cat(paste0("Chromatic RQA with categories: ", paste(chromaNames, collapse = ", "),"\n"))
      }
      cat(" Global Measures\n")
      print(format(outTable$`Global Measures`, digits = 3))
      cat(paste("\n\n Line-based Measures\n"))
      print(format(outTable$`Line-based Measures`, digits = 3))

      if(anisotropyHV){
        id <- which(names(outTable)%in%"Horizontal/Vertical line anisotropy")
        cat(paste0("\n\n",names(outTable)[id],"\n\n"))
        print(format(outTable$`Horizontal/Vertical line anisotropy`, digits = 3))
      }

      if(asymmetryUL){
        id <- which(names(outTable)%in%"Upper/Lower triangle asymmetry")
        cat(paste0("\n\n",names(outTable)[id],"\n\n"))
        cat(" Global Measures\n")
        print(format(outTable[[id]]$`Global Measures`, digits = 3))
        cat(paste("\n\n Line-based Measures\n"))
        print(format(outTable[[id]]$`Line-based Measures`, digits = 3))
      }
      cat("\n~~~o~~o~~casnet~~o~~o~~~\n")
    }

    attr(out,"measuresTable") <- outTable

    return(invisible(out))
}

  # tb <- c("| Max rec. points (post theiler) | N rec. points | RR | Singular points | Divergence | Repetiteveness | Anisotropy |",
  #         "|--------------------------------|---------------|----|-----------------|------------|----------------|------------|".
  #         )
  #
  #cat(c("Global Recurrence Measures\n",paste0(tb,"\n")))


  #   if(is.null(Nboot)){Nboot = 1}
  #
  #   out <- rp_measures_measures(RM,
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
  #       do(rp_measures_measures(RM[row.ind,unlist(.)],
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


#' Get (C)RQA measures from a Recurrence Matrix
#'
#' Use `rp_measures`
#'
#' @inheritParams rp_measures
#'
#' @family Recurrence Quantification Analysis
#'
#' @export
#'
#' @keywords internal
#'
rp_measures_main <- function(RM,
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


  if(is.null(Nboot)){
    Nboot <- 1
  }

  if(Nboot>1|doParallel){
    checkPkg("parallel")
    mc.cores <- parallel::detectCores()-1
    if(Nboot<mc.cores) mc.cores <- Nboot
  }

  NRows <- NROW(RM)
  NCols <- NCOL(RM)

  if(doParallel){

    tstart <- Sys.time()
    cl  <- parallel::makeCluster(mc.cores)
    out <- parallel::mclapply(1, function(i){
      #rp_prep(matrix(RM[ceiling(NCols*NRows*stats::runif(NCols*NRows))], ncol=NCols, nrow = NRows),
      rp_prep(RP = RM,
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
    out <- rp_prep(RP = RM,
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
    cl      <- parallel::makeCluster(mc.cores)
    bootOut <-  parallel::mclapply(1:Nboot, function(i){
      replicate <- as.data.frame(rp_prep(matrix(RM[ceiling(NCols*NRows*stats::runif(NCols*NRows))],
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




#' Diagonal Recurrence Profile
#'
#' @aliases crqa_diagprofile
#' @inheritParams rp_measures
#' @param RM A binary recurrence matrix
#' @param diagWin Window around the line of synchrony
#' @param xname Label for x-axis
#' @param yname Label for y-axis
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
rp_diagProfile <- function(RM,
                           diagWin = NULL,
                           xname = "X-axis",
                           yname = "Y-axis",
                           theiler = 0,
                           DLmin = 2,
                           VLmin = 2,
                           HLmin = 2,
                           DLmax = length(Matrix::diag(RM)),
                           VLmax = length(Matrix::diag(RM)),
                           HLmax = length(Matrix::diag(RM)),
                           doShuffle = FALSE,
                           y1        = NA,
                           y2        = NA,
                           Nshuffle  = 19,
                           AUTO      = NULL,
                           chromatic = FALSE,
                           matrices  = FALSE,
                           doPlot    = TRUE){


  # crqa_results_xy <– crqa(ts1 = lorData$x, ts2 = lorData$y, delay = 9, embed = 4, rescale = 2, radius = 20, normalize = 2, mindiagline = 2, minvertline = 2, tw = 0, whiteline = FALSE, recpt = FALSE, side = “both”) # compute cross-recurrence plot
  #

  if(!all(as.vector(RM)==0|as.vector(RM)==1)){
    stop("Expecting a binary (0,1) matrix.")
  }

  # check auto-recurrence and make sure Matrix has sparse triplet representation
  RM <- rp_checkfix(RM, checkAUTO = TRUE, fixAUTO = TRUE, checkTSPARSE = TRUE, fixTSPARSE = TRUE)

  if(is.null(AUTO)){
    AUTO <- attr(RM,"AUTO")
  }


  RM <- setTheiler(RM = RM, theiler = theiler, silent = TRUE)

  # if(is.na(theiler%00%NA)){
  #   if(!is.null(attributes(RM)$theiler)){
  #     if(attributes(RM)$theiler>0){
  #       message(paste0("Value found in attribute 'theiler'... assuming a theiler window of size:",attributes(RM)$theiler,"was already removed."))
  #     }
  #   }
  #   theiler <- 0
  # } else {
  #   if(theiler < 0){
  #     theiler <- 0
  #   }
  # }

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
      stop("Need time series (y1 and y2) in order to do the shuffle!!")
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
    } else {
      RMd <- rp(y1 = y1, y2 = TSrnd[,(r-1)],emDim = emDim, emLag = emLag, emRad = emRad, to.sparse = TRUE)
    }
    #rp_lineDist(RM,d = diagWin, matrices = TRUE)
    B <- rp_nzdiags(RMd, removeNZ = FALSE, d = diagWin)
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
    df$labels <- factor(df$labels,levels = df$labels, ordered = TRUE)

    out[[r]] <- df
    #rm(df,B,cID,diagID)

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
    dplyr::summarise(`Mean Shuffled` = mean(.data$RR), sdRRrnd = stats::sd(.data$RR)) %>%
    dplyr::ungroup()

  if(Nshuffle==1){
    dy_m$sdRRrnd <- 0
  }

  dy_m$ciHI <- dy_m$`Mean Shuffled` + 1.96*(dy_m$sdRRrnd/sqrt(Nshuffle))
  dy_m$ciLO <- dy_m$`Mean Shuffled` - 1.96*(dy_m$sdRRrnd/sqrt(Nshuffle))

  obs <- dy[dy$.id=="obs",]
  if(doShuffle){
    if(!all(obs$Diagonal%in%dy_m$Diagonal)){
      notInd <- which(!obs$Diagonal%in%dy_m$Diagonal)
      tmp <- matrix(nrow = NROW(dy_m), ncol = NCOL(dy_m), dimnames = list(NULL,colnames(dy_m)))
      for(r in 1:NROW(dy_m)){
        if(!r%in%notInd){
          tmp[r,] <- obs[r,]
        }
      }
      obs <- tmp
      rm(tmp)
    }
  }
  dy_m$Observed <- obs$RR

  df <- tidyr::gather(dy_m, key = "variable", value = "RR", -c(.data$Diagonal,.data$sdRRrnd, .data$labels, .data$ciLO, .data$ciHI))
  df$Diagonal <- as.numeric(df$Diagonal)

  if(doPlot){

    # Diags <- as.numeric(levels(df$labels))
    # Diags[is.na(Diags)] <- 0

    # if(length(diagWin)>21){
    #   ext <- max(min(abs(Diags),na.rm = TRUE),abs(max(Diags,na.rm = TRUE)))-1
    #   breaks <- which(Diags%in%(c(seq(-ext,-1,length.out = 10),0,seq(1,ext,length.out = 10))))
    #   labels <- sort(unique(c(Diags[breaks],0)))
    #   breaks <- sort(unique(c(breaks,stats::median(breaks))))
    # } else {
    #   breaks <- seq_along(Diags)
    #   labels <- sort(unique(c(Diags[breaks],0)))
    # }


    # if(which.min(c(length(breaks),length(labels)))==1){
    #   labels <- labels[1:length(breaks)]
    # } else {
    #   breaks <- breaks[1:length(labels)]
    # }

    breaks <- seq_along(as.numeric(levels(df$labels)))
    labels <- sort(unique(c(as.numeric(levels(df$labels))[breaks],0)))

    x1 <- round((which.min(as.numeric(paste(df$Diagonal))))+(length(diagWin)*.1))
    x2 <- round((which.max(as.numeric(paste(df$Diagonal))))-(length(diagWin)*.1))
    yL <- max(as.numeric(paste(df$RR)),na.rm = TRUE)+0.1
    # col <- c("ciHI" = "grey70", "ciLO" = "grey70", "meanRRrnd" = "grey40","y_obs" = "black")
    if(Nshuffle>0){
      col <- c("Mean Shuffled" = "grey40","Observed" = "black")
      leg <- paste0("Diagonal Profile\n(N surrogates = ",Nshuffle,")")
    } else {
      col <- c("Observed" = "black")
      leg <- paste0("Diagonal Profile")
    }

    #siz <- c("ciHI" = .5, "ciLO" = .5, "meanRRrnd" = .5,"y_obs" = 1)
    #g<- ggplot(df, aes_(x = ~labels, y = ~RR, colour = ~variable))
    #   geom_line() + theme_bw()

    g <- ggplot2::ggplot(df, ggplot2::aes_(x=~Diagonal)) +
      ggplot2::geom_vline(xintercept = df$Diagonal[df$labels=="0"][1], size=1, colour = "grey80")
      if(doShuffle){
        g <- g + ggplot2::geom_ribbon(ggplot2::aes_(ymin=~ciLO, ymax=~ciHI), alpha=0.3)
      }
    g <- g + ggplot2::geom_line(ggplot2::aes_(y=~RR, colour = ~variable), size = .5) +
      #ggplot2::geom_line(ggplot2::aes_(y=~y_obs), colour = "black", size = 1) +
      ggplot2::scale_y_continuous("Recurrence Rate",limits = c(0,1)) +
      ggplot2::scale_x_continuous("Diagonals in recurrence Matrix", breaks = breaks, labels = labels) +
      ggplot2::geom_label(x=x1,y=.8,label=paste0("Recurrences due to\n ",xname),hjust="left", inherit.aes = FALSE, size = 3) +
      ggplot2::geom_label(x=x2,y=.8,label=paste0("Recurrences due to\n ",yname),hjust="right", inherit.aes = FALSE, size = 3) +
      ggplot2::scale_colour_manual(leg, values = col) +
      ggplot2::theme_bw() +
      theme(panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = rel(.6)))
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


#' Fast (C)RQA (command line crp)
#'
#' This function will run the [commandline Recurrence Plots](http://tocsy.pik-potsdam.de/commandline-rp.php) executable provided by Norbert Marwan.
#'
#' @aliases crqa_cl
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
#' @param targetValue A value passed to `est_radius(...,type="fixed", targetMeasure="RR")` if `is.na(emRad)==TRUE`. This is useful for windowed analysis, it will estimate a new radius for each window.
#' @param useParallel Speed up calculations by using the parallel processing options provided by `parallel` to assign a seperate process/core for each window in windowed (C)RQA analysis and [parallel::detectCores()] with`logical = TRUE` to decide on the available cores (default = `FALSE`)
#' @param ... Additional parameters (currently not used)
#'
#' @details The `rp` executable is installed when the function is called for the first time and is renamed to `rp`, from a platform specific filename downloaded from http://tocsy.pik-potsdam.de/commandline-rp.php or extracted from an archive located in the directory:
#' `...\\casnet\\commandline_rp\\`
#' The file is copied to the directory: `...\\casnet\\exec\\`
#' The latter location is stored as an option and can be read by calling `getOption("casnet.path_to_rp")`.
#'
#' @section Troubleshooting:
#' Some notes on resolving errors with `rp`.The script will first try to download the correct executable, if that fails, the copy will have... failed. It should be relatively easy to get `rp_cl()` working though, by using some custom settings:
#'
#'
#' + *Copy failed* - Every time the function `rp_cl()` is called it will check whether a log file `rp_instal_log.txt` is present in the `...\\casnet\\exec\\` directory. If you delete the `rp_instal_log.txt` file, and call the function, another attempt will be made to download a copy of the executable.
#'
#' + *Copy still fails and/or no permission to copy* - If you cannot acces the directory `...\\casnet\\commandline_rp\\`, download the appropriate executable from the [Commandline Recurrence Plots](http://tocsy.pik-potsdam.de/commandline-rp.php) page and copy to a directory you do have the rights to: *execute* commands, *write* and *read* files. Make sure you rename the file to `rp` (`rp.exe` on Windows OS). Then, either pass the path to `rp` as the argument `path_to_rp` in the `rp_cl(.., path_to_rp = "YOUR_PATH")` function call, or, as a more permanent solution, set the `path_to_rp` option by calling `options(casnet.path_to_rp="YOUR_PATH")`.
#'
#' + *Error in execution of `rp`* - This can have a variety of causes, the `rp` executable is called using [system2()] and makes use of the [normalizePath()] function with argument `mustWork = FALSE`. Problems caused by specific OS, machine, or, locale problems (e.g. the `winslash` can be reported as an \href{https://github.com/FredHasselman/casnet/issues}{issue on Github}). One execution error that occurs when the OS is not recognised properly can be resolved by chekcing `getOption("casnet.rp_prefix")`. On Windows OS this should return an empty character vector, on Linux or macOS it should return `"./"`. You can manually set the correct prefix by calling `options(casnet.rp_prefix="CORRECT OS PREFIX")` and fill in the prefix that is correct for your OS
#'
#'
#' @return A list object containing 1-3 elements, depending on arguments requesting output.
#'
#' * `rqa_measures` - A list of the (C)RQA measures returned if `returnMeasures = TRUE`:
#'      1. RR = Recurrence rate
#'      2. DET = Determinism
#'      3. DET_RR = Ratio DET/RR
#'      4. LAM = Laminarity
#'      5. LAM_DET = Ratio LAM/DET
#'      6. L_max = maximal diagonal line length
#'      7. L_mean = mean diagonal line length
#'      8. L_entr = Entropy of diagonal line length distribution
#'      9. DIV =  Divergence (1/L_max)
#'      10. V_max = maximal vertical line length
#'      11. TT = Trapping time
#'      12. V_entr = Entropy of vertical line length distribution
#'      13. T1 = Recurrence times 1st type
#'      14. T2 = Recurrence times 2nd type
#'      15. W_max = Max interval length
#'      16. W_mean = Mean of interval lengths
#'      17. W_entr = Entropy of interval length distribution
#'      18. W_prob = Probability of interval
#'      19. F_min = F min
#' * `rqa_rpvector` - The radius thresholded distance matrix (recurrence matrix), which can be visualised as a recurrence plot by calling [rp_plot()]. If a sliding window analysis is conducted this will be a list of matrices and could potentially grow too large to handle. It is recommended you save the output to disk by setting `saveOut = TRUE`.
#' * `rqa_diagdist` - The distribution of diagonal line lengths
#'
#'
#' @note The platform specific `rp` command line executables were created by Norbert Marwan and obtained under a Creative Commons License from the website of the Potsdam Institute for Climate Impact Research at [http://tocsy.pik-potsdam.de/](http://tocsy.pik-potsdam.de/).
#'
#' The full copyright statement on the website is as follows:
#'
#' (C) 2004-2017 SOME RIGHTS RESERVED
#'
#' University of Potsdam, Interdisciplinary Center for Dynamics of Complex Systems, Germany
#'
#' Potsdam Institute for Climate Impact Research, Transdisciplinary Concepts and Methods, Germany
#'
#' This work is licensed under a [Creative Commons Attribution-NonCommercial-NoDerivs 2.0 Germany License](https://creativecommons.org/licenses/by-nc-nd/2.0/de/).
#'
#' More information about recurrence analysis can be found on the [Recurrence Plot](http://www.recurrence-plot.tk) website.
#'
#' @family Recurrence Quantification Analysis

#' @export
#'
rp_cl            <- function(y1,
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

  if(!file.exists(normalizePath(file.path(getOption("casnet.path_to_rp"),"rp_install_log.txt"), mustWork = FALSE))){
    set_command_line_rp()
  }

  os <- get_os()
  sysinf <- Sys.info()
  if(os%in%"osx"&as.numeric(strsplit(sysinf[which(names(sysinf)%in%"release")],"[.]")[[1]][1])>=19){
    stop("As of macOS Catalina, 32-bit applications are no longer supported, please use rp_measures() to conduct Recurrence Quantification Analysis.")
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
      path_out <- tempdir(check = TRUE)
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
  #parallel::clusterExport(cl = cl, c("rp_cl_main","wIndices","df","y1","y2","emDim","emLag","emRad","DLmin","VLmin","theiler", "win","step","JRP","distNorm","returnMeasures","returnRPvector","returnLineDist","doPlot","path_to_rp","saveOut","path_out","file_ID","silent","..."))

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

    # parallel::clusterEvalQ(cl, library(callr))
    parallel::clusterEvalQ(cl,library(utils))
    parallel::clusterEvalQ(cl,library(plyr))
    parallel::clusterEvalQ(cl,library(tidyverse))
    parallel::clusterEvalQ(cl,library(pROC))
    parallel::clusterEvalQ(cl,library(Matrix))
    parallel::clusterEvalQ(cl,library(casnet))

    # parallel::clusterExport(cl, varlist = c("data","emDim","emLag","emRad","DLmin","VLmin","theiler","win","step","JRP","distNorm","returnMeasures","returnRPvector","returnLineDist","doPlot","path_to_rp", "saveOut","path_out","file_ID","silent","targetValue", "useParallel"))

    # cluster_library(c("devtools","utils","plyr","dplyr","tidyr","Matrix","pROC")) %>%
    # cluster_assign_value("rp_cl_main", rp_cl_main) %>%
    # cluster_assign_value("est_radius", est_radius) %>%
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

    #  parallel::clusterExport(cl = cl, c("rp_cl_main","wIndices","df","y1","y2","emDim","emLag","emRad","DLmin","VLmin","theiler", "win","step","JRP","distNorm","returnMeasures","returnRPvector","returnLineDist","doPlot","path_to_rp","saveOut","path_out","file_ID","silent","..."))


    start <- proc.time()
    wList <- parallel::parLapply(cl,dfList,function(df){rp_cl_main(data = df,
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
      rp_cl_main(
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
        cowplot::plot_grid(plotList)
      } else {
        rp(y1 = df[,1], y2 = df[,2], emDim = emDim, emLag = emLag, emRad = emRad, doPlot = TRUE)
      }
    }
    if(wPlot==3){
      if(windowedAnalysis){
        plotList <- plyr::llply(wIndices, function(ind) rp(y1 = df[ind,1],y2 = df[ind,2], emDim = emDim, emLag = emLag, doPlot = TRUE))
        cowplot::plot_grid(plotList)
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


#' rp_cl_main
#'
#' @inheritParams rp_cl
#'
#' @keywords internal
#'
#' @export
#'
rp_cl_main <- function(data,
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
        emRad <- est_radius(y1 = y1,
                            emDim = emDim,
                            emLag = emLag,
                            targetValue = targetValue,
                            radiusOnFail = "percentile", tol = .2, silent = silent)
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
        emRad <- est_radius(y1 = y1, y2 = y2, emDim = emDim, emLag = emLag, targetValue = targetValue, tol = .2, radiusOnFail = "percentile", silent = silent)
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

  measures     <- return_error(utils::read.delim(normalizePath(gsub("[']+","",measOUT)),header=TRUE))
  rpMAT        <- return_error(utils::read.delim(normalizePath(gsub("[']+","",plotOUT)),header=TRUE))
  disthistDiag <- return_error(utils::read.delim(normalizePath(gsub("[']+","",histOUTdiag)), header=FALSE, sep = " "))
  disthistHori <- return_error(utils::read.delim(normalizePath(gsub("[']+","",histOUThori)), header=FALSE, sep = " "))

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
rp_nzdiags <- function(RM=NULL, d=NULL, returnVectorList=TRUE, returnNZtriplets=FALSE, removeNZ = TRUE, silent = TRUE){
  # Loosely based on MATLAB function spdiags() by Rob Schreiber - Copyright 1984-2014 The MathWorks, Inc.

  # if(!is.na(win)){
  #   if(length(win)==1){
  #     win <- c(-win,win)
  #   }
  #   if(length(win)==2){
  #       win <- sort(win)
  #   }
  #   if(!length(win)%[]%c(1,2)){
  #     stop("Windowsize must be NA, or a 1, or 2 element vector!")
  #     }
  # }

  if(grepl("matrix",class(RM)[1],ignore.case = TRUE)){

    if(all(RM>0)){warning("All matrix elements are nonzero.")}

    s  <- Sys.time()

    nzdiagsM <- methods::as(RM, "dgTMatrix")
    nzdiags  <- data.frame(row   = nzdiagsM@i,
                           col   = nzdiagsM@j,
                           value = nzdiagsM@x,
                           ndiag = (nzdiagsM@j)-(nzdiagsM@i))
    nzdiags <- dplyr::arrange(nzdiags,nzdiags$ndiag)

    if(is.null(d)){
      d <-  c(-1*rev(1:(NCOL(nzdiagsM)-1)),0,1:(NCOL(nzdiagsM)-1))
    }
    rm(nzdiagsM)

    if(removeNZ){
      nd <- unique(nzdiags$ndiag)
      # Get diagonals which have nonzero elements
      d <- nd[nd%in%sort(as.vector(d))]
    }

    indMat  <- col(RM)-row(RM)
    indMatV <- as.vector(indMat)

    #cat("\nStart extraction of nonzero diagonals\n")
    m  <- NROW(RM)
    n  <- NCOL(RM)
    p  <- length(d)
    if(is.logical(RM)){
      B <- matrix(FALSE,nrow = max(c(m,n), na.rm = TRUE), ncol = p)
    } else {
      B <- matrix(0, nrow = max(c(m,n), na.rm = TRUE), ncol = p)
    }
    colnames(B) <- paste(d)

    for(i in seq_along(d)){
      B[(nzdiags$row[nzdiags$ndiag==d[i]]+1), i] <- nzdiags$value[nzdiags$ndiag==d[i]]
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
#' @inheritParams rp
#' @inheritParams rp_measures
#' @param d Vector of diagonals to be extracted from matrix `RP` before line length distributions are calculated. A one element vector will be interpreted as a windowsize, e.g., `d = 50` will extract the diagonal band `-50:50`. A two element vector will be interpreted as a band, e.g. `d = c(-50,100)` will extract diagonals `-50:100`. If `length(d) > 2`, the numbers will be interpreted to refer to individual diagonals, `d = c(-50,50,100)` will extract diagonals `-50,50,100`. If `length(d)` is `NULL`, 1 or 2, the theiler window is applied before diagonals are extracted. The theiler window is ignored if `length(d)>2`, or if it is larger than the matrix or band indicated by parameter `d`. A warning will be given is a theiler window was already applied to the matrix.
#' @param invert Relevant for Recurrence Time analysis: Return the distribution of 0 valued segments in nonzero diagonals/verticals/horizontals. This indicates the time between subsequent line structures.
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

  if(length(d)<=2){
    suppressMessages(RM <- setTheiler(RM, theiler))
  } else {
    if(!is.null(attributes(RM)$theiler)){
      message(paste0("Value found in attribute 'theiler'... assuming a theiler window was already applied to the matrix."))
    }
  }

  if(invert){
    RM <- Matrix::Matrix(1-RM,sparse = TRUE)
  }

  if(!is.null(d)){
    if(length(d)==1){d <- seq(-d,d)}
    if(length(d)==2){d <- seq(min(d),max(d))}
  }


  B <- rp_nzdiags(RM)
  V <- Matrix::as.matrix(RM)[,colSums(Matrix::as.matrix(RM))>0]

  RPt <- Matrix::t(RM)

  H <- Matrix::as.matrix(RPt)[,colSums(Matrix::as.matrix(RPt))>0]
  rm(RPt)

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


#' Copy Matrix Attributes
#'
#' Simple attribute copy used in `casnet` to convert between `matrix` and `Matrix` classes and back.
#'
#' @param source Source matrix
#' @param target Target matrix
#' @param source_remove Remove these attribute fields from the source before copying.
#'
#' @return The target matrix with attributes copied from the source matrix.
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
    RM <- Matrix::Matrix(RM, sparse = TRUE)
    yesS4 <- isS4(RM)
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
    if(class(RM)[1]%in%names(methods::getClass("TsparseMatrix")@subclasses)){
      yesTSPARSE <- TRUE
    }
  }
  if(fixTSPARSE&!yesTSPARSE){
    if(!yesS4){
      RM <- Matrix::Matrix(RM, sparse = TRUE)
    }
    Mtype <- gsub("CMatrix","TMatrix",class(RM)[1])
    eval(parse(text=paste0("RM <- as(RM,'",Mtype,"')")))
  }


  RM <- rp_copy_attributes(source = dummy, target = RM, source_remove = c("Dimnames", "i", "class","Dim", "p","x","factors"))

  return(RM)
}


#' Plot (thresholded) distance matrix as a recurrence plot
#'
#' @param RM A distance matrix or recurrence matrix, preferably generated by function [rp] or [rn].
#' @param plotDimensions Should the state vectors be plotted if they are available as attributes of RM (default = `TRUE`)
#' @param plotDimensionLegend If `plotDimensions = TRUE` plot a simple legend (default = `FALSE`)
#' @param plotMeasures Print common (C)RQA measures in the plot if the matrix is binary (default = `FALSE`)
#' @param plotRadiusRRbar The `Radius-RR-bar` is a colour-bar guide plotted with an unthresholded distance matrix indicating a number of `RR` values one would get if a certain distance threshold were chosen (default = `TRUE`)
#' @param drawGrid Draw a grid on the recurrence plot (default = `FALSE`)
#' @param drawDiagonal One usually omits the main diagonal, the Line Of Incidence (LOI) from the calculation of Auto-RQA measures., it is however common to plot it. Set to `FALSE` to omit the LOI from an Auto-Recurrence Plot (default = `TRUE`)
#' @param drawNA Draw `NA` values in the recurrence plot? If `FALSE` then `NA` values will be considered non-recurring (default = `FALSE`)
#' @param markEpochsLOI Pass a factor whose levels indicate different epochs or phases in the time series and use the line of identity to represent the levels by different colours (default = `NULL`)
#' @param radiusValue If `plotMeasures = TRUE` and RM is an unthresholded matrix, this value will be used to calculate recurrence measures. If `plotMeasures = TRUE` and RM is already a binary recurrence matrix, pass the radius that was used as a threshold to create the matrix for display purposes. If `plotMeasures = TRUE` and `radiusValue = NA`, function `est_radius()` will be called with default settings (find a radius that yields `.05` recurrence rate). If `plotMeasures = FALSE` this setting will be ignored.
#' @param courseGrain Reduce the size of the matrix before plotting (greatly improves speed). If `plotMeasures = TRUE` coursegraining takes place after calculation of the CRQA measures. See [mat_coursegrain()] for details (default = `TRUE`)
#' @param maxSize The maximum size of the matrix (cols x rows) above which the coursegraining function [mat_coursegrain] will be implemented if `courseGrain = TRUE`. This will speed up the plotting process
#' @param title A title for the plot
#' @param xlabel An x-axis label
#' @param ylabel An y-axis label
#' @param plotSurrogate Should a 2-panel comparison plot based on surrogate time series be added? If `RM` has attributes `y1` and `y2` containing the time series data (i.e. it was created by a call to [rp]), the following options are available: "RS" (random shuffle), "RP" (randomised phases), "AAFT" (amplitude adjusted fourier transform). If no timeseries data is included, the columns will be shuffled.  NOTE: This is not a surrogate test, just 1 surrogate is created from `y1`. (default = `FALSE`)
#' @param returnOnlyObject Return the ggplot object only, do not draw the plot (default = `TRUE`)
#'
#' @return A nice plot of the recurrence matrix.
#' @export
#'
#' @family Distance matrix operations (recurrence plot)
#'
rp_plot <- function(RM,
                    plotDimensions = FALSE,
                    plotDimensionLegend = FALSE,
                    plotMeasures = FALSE,
                    plotRadiusRRbar = TRUE,
                    drawGrid = FALSE,
                    drawDiagonal = TRUE,
                    drawNA = FALSE,
                    markEpochsLOI = NULL,
                    radiusValue = NA,
                    courseGrain = TRUE,
                    maxSize = 500^2,
                    title = "",
                    xlabel = "",
                    ylabel="",
                    plotSurrogate = NA,
                    returnOnlyObject = FALSE){

  useGtable <- TRUE

  if(is.null(attr(RM,"chromatic"))){
    chromatic <- FALSE
    attr(RM,"chromatic") <- chromatic
  } else {
    chromatic <- attr(RM,"chromatic")
  }

  if(!all(stats::na.exclude(as.vector(RM))%in%c(0,1))&!chromatic){
    unthresholded <- TRUE
  } else {
    unthresholded <- FALSE
  }

  #  else {
  #   # This is just to set the attribute in case
  #   if(unthresholded|plotMeasures){
  #     radiusValue <- est_radius(RM,silent = TRUE)$Radius
  #     attr(RM,"emRad") <- radiusValue
  #   }
  # }

  # Check if we need to display CRQA measures for the original matrix

  # Get Radius attribute (if present), otherwise set it (if not present)
  if(is.na(radiusValue%00%NA)){
    if(!is.na(attr(RM,"emRad")%00%NA)){
      radiusValue <- attr(RM,"emRad")
    }
  } else {
    if(is.numeric(radiusValue)&radiusValue>=0){
      if(is.na(attr(RM,"emRad")%00%NA)){
        attr(RM,"emRad") <- radiusValue
      }
    } else {
      radiusValue <- NA
    }
  } # is.na(radiusvalue)

  # check auto-recurrence and make sure Matrix has sparse triplet representation
  RM   <- rp_checkfix(RM, checkAUTO = TRUE, fixAUTO = TRUE, checkSPARSE = TRUE)
  AUTO <- attr(RM,"AUTO")


  # Get CRQA measures if requested
  if(plotMeasures){
      if(!is.na(attr(RM,"emRad")%00%NA)){
        radiusValue <- attr(RM,"emRad")
      } else {
        radiusValue <- est_radius(RM, silent = TRUE)$Radius
        attr(RM,"emRad") <- radiusValue
      }

    if(unthresholded){
      rpOUT   <- rp_measures(RM, emRad = radiusValue, AUTO = AUTO)
    } else {
      rpOUT   <- rp_measures(RM, AUTO = AUTO)
    }
  } # plotMeasures


  # Size check
  reduced <- FALSE
  if((rp_size(RM)$rp_size_total>=maxSize)&courseGrain){
    message("NOTE: To speed up the plotting process, the RP will represent a coursegrained matrix. Set argument 'courseGrain = FALSE' to see the full matrix.")
    if((rp_size(RM)$rp_size_total/2)>=maxSize){
      target_height <- round(sqrt(maxSize))
      target_width  <- round(sqrt(maxSize))
    } else {
      target_height <- NROW(RM)/2
      target_width  <- NCOL(RM)/2
    }


    RM <- mat_coursegrain(RM, target_height = target_height, target_width = target_width)
    reduced <- TRUE
  }

  colvec <- c("#FFFFFF00","#000000")
  names(colvec) <- c("0","1")


  # Make sure we can draw a diagonal if requested
  maxDist <- max(RM, na.rm = TRUE)
  if(drawDiagonal&AUTO){
    if(all(stats::na.exclude(as.vector(RM))%in%c(0,1))){
      RM <- bandReplace(RM,0,0,1)
    } else {
      if(!chromatic){
        RM <- bandReplace(RM,0,0,(maxDist+1))
      }
    }
  }


  # prepare data
  #if(attr(RM,"package")%00%""%in%"Matrix"){
  if(any(grepl("Matrix",class(RM)))){
    RM     <- rp_checkfix(RM, checkTSPARSE = TRUE, fixTSPARSE = TRUE)
    #meltRP <- data.frame(Var1 = (RM@i+1), Var2 = (RM@j+1), value = as.numeric(RM@x))
    meltRP <-  mat_mat2ind(Matrix::as.matrix(RM))
  } else {
    meltRP <-  mat_mat2ind(as.matrix(RM))
  }

  hasNA <- FALSE

  if(chromatic){
    if(is.null(attr(RM,"chromaNames"))){
      chromaNumbers <- sort(unique(meltRP$value))[sort(unique(meltRP$value))>0]
      names(chromaNumbers) <- c("No recurrence", paste("Cat.", chromaNames))
      attr(RM,"chromaNames") <- chromaNumbers
      chromaNames <- names(chromaNumbers)
    } else {
      chromaNumbers <- attr(RM,"chromaNames")
      Cind <- sort(chromaNumbers, index.return = TRUE)
      chromaNames <- names(chromaNumbers)[Cind$ix]
      chromaNumbers <- chromaNumbers[Cind$ix]
    }
    if(NA%in%chromaNames){
      hasNA <- TRUE
    }
  }


  if(drawNA){
    if(!is.null(attr(RM,"NAij"))){
      hasNA <- TRUE
      NAid <- data.frame(attr(RM,"NAij"))
      NAid$value <- NA
      colnames(NAid)[1:2] <- c("Var1","Var2")
      #meltRP$value[meltRP$Var1==NAid[,2]&meltRP$Var2==NAid[,1]] <- NA
      meltRP <- rbind(meltRP,NAid)
    } else {
      message("Cannot find any NA to draw!")
      drawNA <- FALSE
      hasNA <- FALSE
    }
  } else {
    if(any(is.na(meltRP$value))){
      hasNA <- FALSE
      meltRP$value[is.na(meltRP$value)] <- (maxDist + 1)
    }
  }


  # Unthresholded START ----
  showL <- FALSE
  if(unthresholded){

    if(!is.null(attr(RM,"weighted"))){
      if(attr(RM,"weighted")%00%FALSE){
        warning("Can't show the radius vs. RR bar on a weighted Recurrence Plot!")
        plotRadiusRRbar <- FALSE
      }
    }

    if(!is.null(markEpochsLOI)){
      warning("Can't show epochs on an unthresholded Recurrence Plot!")
    }

    # if(chromatic){
    #   plotRadiusRRbar <- FALSE
    # }

  } else { # unthresholded

    # THRESHOLDED = TRUE
    unthresholded <- FALSE
    plotRadiusRRbar <- FALSE

    if(!is.null(markEpochsLOI)){
      if(is.factor(markEpochsLOI)&length(markEpochsLOI)==max(c(NROW(RM),NCOL(RM)))){
        start <- max(meltRP$value, na.rm = TRUE) + 1
        #cpal <- paletteer::paletteer_d(package = "rcartocolor",palette = "Safe", n = nlevels(markEpochsLOI))
        cpal <- getColours(Ncols = nlevels(markEpochsLOI))
        #cpal <- paletteer::paletteer_d(package = "ggthemes",palette = "tableau_colorblind10",n = nlevels(markEpochsLOI),direction = 1)

        if(hasNA){
          colvec <- c("#FFFFFF00","#000000","#FF0000", cpal)
          names(colvec) <- c("0","1","NA",levels(markEpochsLOI))
        } else {
          colvec <- c("#FFFFFF00","#000000", cpal) #viridis::viridis_pal()(nlevels(markEpochsLOI)))
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

      if(chromatic){

        if("No recurrence"%in%chromaNames){
          chromaNames[which("No recurrence"%in%chromaNames)] <- "0"
        } else {
          if(!"0"%in%chromaNames){
          chromaNames <- c("0",chromaNames)
          }
        }

        if(!0%in%chromaNumbers){
          chromaNumbers <- c(0,chromaNumbers)
          names(chromaNumbers)[1] <- "0"
        }

        colvec <- getColours(length(chromaNames))
        colvec[which(chromaNames%in%"0")] <- "#FFFFFF00"

        if(hasNA){
          if(NA%in%chromaNames){
            colvec[which(chromaNames%in%NA)] <-"#FF0000"
            names(colvec) <- chromaNames
          } else {
            colvec <- c(colvec,"#FF0000")
            names(colvec) <- c(chromaNames,NA)
            if(!NA%in%chromaNumbers){
              chromaNumbers <- c(chromaNumbers, NA)
            }
          }
        } else {
          names(colvec) <- chromaNames
        }

        if(!drawNA){
          chromaNames   <- stats::na.exclude(chromaNames)
          chromaNumbers <- stats::na.exclude(chromaNumbers)
        }

        # CHECK VECTOR LENGTHS
        allEqual <- FALSE
        if(length(unique(meltRP$value))==length(chromaNames)){
          if(length(unique(meltRP$value))==length(chromaNumbers)){
            allEqual <- TRUE
          }
        }

        if(allEqual){
          meltRP$value <- factor(meltRP$value, levels = as.numeric(chromaNumbers), labels = chromaNames, exclude = NULL)
        } else {
          #message("Check values and chromaNames...")

          excludeVals <- unique(meltRP$value)[!unique(meltRP$value)%in%chromaNumbers]
          meltRP$value <- factor(meltRP$value,
                                 levels = c(as.numeric(chromaNumbers), excludeVals),
                                 labels = c(chromaNames, as.character(excludeVals)),
                                 exclude = c(as.character(excludeVals),""))
        }

        showL <- FALSE
        tmp <- data.frame(x = 1:length(chromaNames), y = 1:length(chromaNames), cols = chromaNames)

        grL <- ggplot(tmp, aes_(x = ~x, y = ~y)) +
          geom_point(aes_(colour = ~cols, fill = ~cols), pch=22, size =10) +
          scale_fill_manual(name  = "Phases", breaks = chromaNames,
                            values = colvec,
                            na.translate = TRUE,
                            na.value = "#FF0000") +
          scale_colour_manual(name  = "Phases", breaks = chromaNames,
                              values = colvec,
                              na.translate = TRUE,
                              na.value = "#FF0000")
        grLegend <- cowplot::get_legend(grL)
        rm(grL)

      } else {

      colvec <- c("#FFFFFF00","#000000")
      names(colvec) <- c("0","1")
      meltRP$value <- factor(meltRP$value, levels = c(0,1), labels = c("0","1"))
      }
    }
  } # unthresholded = FALSE

  ## Unthresholded END ##

  # Main plot ----
  gRP <-  ggplot2::ggplot(ggplot2::aes_(x=~Var1, y=~Var2, fill = ~value), data= meltRP) +
    ggplot2::geom_raster(hjust = 0, vjust=0, show.legend = showL)

  if(drawDiagonal){
    gRP <- gRP + ggplot2::geom_abline(slope = 1,colour = "grey30", size = 1)
  }

  ## Unthresholded START ----
  if(unthresholded){

    # #barValue <- 0.05
    # if(is.na(radiusValue)){
    #   barValue <- mean(meltRP$value, na.rm = TRUE)
    #   } else {
    #   barValue <- est_radius(RM, targetValue = radiusValue, silent = TRUE, maxIter = 100, radiusOnFail = "percentile")$Radius
    # }

    if(!is.na(radiusValue)){
      barValue <- radiusValue
    } else {
      barValue <- mean(meltRP$value, na.rm = TRUE)
    }

    if(attr(RM,"weighted")){

      gRP <- gRP + ggplot2::scale_fill_gradient(low      = "white",
                                                high     = "red3",
                                                na.value = "#FF0000",
                                                space    = "Lab",
                                                name     = "")

    } else {




      gRP <- gRP + ggplot2::scale_fill_gradient2(low      = "red3",
                                                 high     = "steelblue",
                                                 mid      = "white",
                                                 na.value = "#FF0000",
                                                 midpoint = barValue*1.1, #mean(meltRP$value, na.rm = TRUE),
                                                 limit    = c(min(meltRP$value, na.rm = TRUE),max(meltRP$value, na.rm = TRUE)),
                                                 space    = "Lab",
                                                 name     = "")

    } # Weighted


    ## RadiusRRbar START ----
    if(plotRadiusRRbar){
      # Create a custom legend
      distrange  <- round(seq(0,max(RM,na.rm = TRUE),length.out=7),2)
      resol      <- sort(unique(round(as.vector(RM),2)))
      if(length(resol)<7){
        resol <- distrange
      }
      if(length(resol)>100){
        resol <- round(seq(0,max(RM,na.rm = TRUE),length.out=100),2)
      }
      resol <- resol %>% tibble::as_tibble() %>% dplyr::mutate(y= seq(exp(0),exp(1),length.out=NROW(resol)), x=0.5)

      distrange <- plyr::ldply(c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5), function(t){
        suppressMessages(est_radius(RM,targetValue = t, silent = TRUE, maxIter = 100, radiusOnFail = "percentile"))
      })

      if(AUTO){
      maxD <- suppressMessages(max(RM[lower.tri(RM)], na.rm = TRUE))
      } else {
      maxD <- suppressMessages(max(RM, na.rm = TRUE))
      }
      RecScale <- data.frame(RR=distrange$Measure,epsilon=distrange$Radius)
      RecScale <- RecScale %>%
        dplyr::add_row(epsilon=mean(c(0,distrange$Radius[1])),RR=mean(c(0,distrange$Measure[1])),.before = 1) %>%
        dplyr::add_row(epsilon=maxD,RR=1)

      resol$y <- elascer(x = resol$y,lo = min(log(RecScale$RR),na.rm = TRUE), hi = max(log(RecScale$RR),na.rm = TRUE))
      resol <- resol[-1,]

      if(!is.na(radiusValue)){
        barValue <- round(RecScale$RR[which(round(RecScale$epsilon,4)>=radiusValue)[1]],4)
        barValue <- resol$value[which(resol$y>=log(barValue))[1]]
      } else {
        barValue <- mean(meltRP$value, na.rm = TRUE)
      }

      gDist <-  ggplot2::ggplot(resol,ggplot2::aes_(x=~x,y=~y,fill=~value)) +
        ggplot2::geom_tile(show.legend = FALSE) +
        ggplot2::scale_y_continuous(name = "Recurrence Rate", breaks = log(RecScale$RR), labels = paste(round(RecScale$RR,3)), sec.axis = dup_axis(name=expression(paste("recurrence threshold",~ epsilon)), labels = paste(round(RecScale$epsilon,2)))) +
        ggplot2::scale_fill_gradient2(low      = "red3",
                                      high     = "steelblue",
                                      mid      = "white",
                                      na.value = "#FF0000",
                                      midpoint =  barValue * 1.1,#mean(meltRP$value, na.rm = TRUE),
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

  }
  # else { # unthresholded
  #
  # warning("This is a binary matrix. Cannot plot the Radius vs. RR bar!")
  #
  # }
  ##RadiusRRbar END##


  ## Create Theme START ----
  rptheme <- ggplot2::theme_bw() + ggplot2::theme(panel.background = element_blank(),
                                                  panel.grid.minor  = element_blank(),
                                                  panel.border = element_rect("grey50",fill=NA),
                                                  legend.key = element_rect(colour = "grey90"),
                                                  axis.ticks = element_blank(),
                                                  axis.text = element_blank(),
                                                  plot.margin = margin(0,0,0,0))

  if(drawGrid){
    rptheme <- rptheme +  theme(panel.grid.minor  = element_blank(),
                                panel.grid.major  = element_line("grey80",size = .1),
                                panel.ontop = unthresholded)
  } else {
    rptheme <- rptheme +  theme(panel.grid.major  = element_blank(),
                                panel.grid.minor  = element_blank(),
                                panel.ontop = FALSE)
  } #drawGrid


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

  if(!is.null(markEpochsLOI)){
    if(!unthresholded){
      gRP <- gRP +  ggplot2::scale_fill_manual(name  = "Key:",
                                               values = colvec,
                                               na.translate = TRUE ,
                                               na.value = scales::muted("slategray4"),
                                               guide = "legend",
                                               limits = levels(meltRP$value))
    }

  } else {
    if(!unthresholded){

      if(chromatic){

      gRP <- gRP +  ggplot2::scale_fill_manual(name  = "", breaks = sort(unique(meltRP$value)),
                                               values = colvec,
                                               na.translate = TRUE ,
                                               na.value = scales::muted("slategray4"))

      } else {

      gRP <- gRP +  ggplot2::scale_fill_manual(name  = "", breaks = c(0,1),
                                               values = colvec,
                                               na.translate = TRUE ,
                                               na.value = scales::muted("slategray4"),
                                               guide = "none")
      }
    }
  } # markEpochsLOI


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

  if(!is.null(attr(RM,"emDims1.name"))|nchar(xlabel)>0){
    xdims <- ifelse(nchar(xlabel)>0, xlabel, attr(RM,"emDims1.name"))
  }
  if(!is.null(attr(RM,"emDims2.name"))|nchar(ylabel)>0){
    ydims <- ifelse(nchar(ylabel)>0, ylabel, attr(RM,"emDims2.name"))
    if(AUTO){
      ydims <- xdims
    }
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
      warning("No atrribute with dimensions!")
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


    if(plotDimensionLegend){

      # y1 <- data.frame(t1=attr(RM,"emDims1"))
      # y2 <- data.frame(t2=attr(RM,"emDims2"))

      y1l <- data.frame(Value     = rep(0,length(unique(y1$Dimension)),each=2),
                        tm        = rep(0,length(unique(y1$Dimension)),each=2),
                        Dimension =  rep(unique(y1$Dimension),each=2))

      gdl <- ggplot2::ggplot(y1l, ggplot2::aes_(y=~Value, x=~tm, colour=~Dimension)) +
        ggplot2::geom_line(show.legend = TRUE, size = 2) + ggplot2::scale_color_grey() +
        scale_x_continuous() +
        theme_void() + theme(legend.text = element_text(size=rel(2)),
                             legend.title = element_text(size=rel(2)),
                             legend.position = c(.5,.5))


    }

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

  # Build Graph using gtable

  g <- ggplot2::ggplotGrob(gRP)

  if(plotDimensions){
    gry2<-ggplot2::ggplotGrob(gy2)
    gry1<-ggplot2::ggplotGrob(gy1)
    if(plotDimensionLegend){
      grDimL <- ggplot2::ggplotGrob(gdl)
    }
  }

  if(unthresholded&plotRadiusRRbar){
    grDist <- ggplot2::ggplotGrob(gDist)
  }

  if(plotMeasures){
    grA <- ggplot2::ggplotGrob(gA)
  }

  # cat(paste("\n\nplotMeasures =", plotMeasures,"\n"))
  # cat(paste("plotDimensions =", plotDimensions,"\n"))
  # cat(paste("plotDimensionLegend =", plotDimensionLegend,"\n"))
  # cat(paste("unthresholded =", unthresholded,"\n"))
  # cat(paste("plotRadiusRRbar =", plotRadiusRRbar,"\n"))

  # Build list for the gtable
  Lmat <- list()
  if(plotMeasures){Lmat[[1]] <- grA} else {Lmat[[1]] <- NA}
  if(plotMeasures){Lmat[[2]] <- grid::nullGrob()} else {Lmat[[2]] <- NA}
  if(plotDimensions){Lmat[[3]] <- gry2} else {Lmat[[3]] <- NA}
  if(plotDimensionLegend){Lmat[[4]] <- grDimL} else {if(plotDimensions){Lmat[[4]] <- grid::nullGrob()} else {Lmat[[4]] <- NA}}
  Lmat[[5]] <- g
  if(plotDimensions){Lmat[[6]] <- gry1} else {Lmat[[6]] <- NA}
  if((unthresholded&plotRadiusRRbar)|chromatic){
    if(chromatic){
      Lmat[[7]] <- grLegend
    } else {
      Lmat[[7]] <- grDist
    }
  } else {
    Lmat[[7]] <- NA
  }
  if(plotDimensions&(unthresholded&plotRadiusRRbar|chromatic)){Lmat[[8]] <- grid::nullGrob()} else {Lmat[[8]] <- NA}

  # Remove unused fields
  Lmat[is.na(Lmat)] <- NULL

  w <- c(.35,.25, 1,.5,.35)[c(plotMeasures,(plotMeasures|plotDimensions),TRUE,plotRadiusRRbar,chromatic)]
  h <- c(1,.25)[c(TRUE,plotDimensions)]

  widths  <- unit(w,"null")  #ifelse(plotMeasures, unit(c(.35,.25, 1,.5),  unit(c(.25, 1), "null")))
  heights <- unit(h,"null") # ifelse(plotMeasures, unit(c(1), "null"), unit(c(1,.25), "null"))

  if(plotDimensions){
    mrows <- 2
  } else {
    mrows <- 1
  }

  mat <- matrix(Lmat,nrow = mrows)
  gt  <- gtable::gtable_matrix("rp", mat, widths = widths, heights =  heights, respect = TRUE)

  #  gtable::gtable_matrix("bi_rp", mat, widths = unit(c(1), "null"), heights =  unit(c(1), "null"), respect = TRUE)

  # if(plotDimensions&plotDimensionLegend&!plotMeasures&unthresholded&!plotRadiusRRbar){
  #   mat <- matrix(list(gry2, grDimL, g, gry1),nrow = 2)
  #   gt  <- gtable::gtable_matrix("di_rp_dim", mat, widths = unit(c(.25, 1), "null"), heights =  unit(c(1,.25), "null"),respect = TRUE)
  # }
  #
  # if(plotDimensions&plotDimensionLegend&!plotMeasures&unthresholded&plotRadiusRRbar){
  #   mat <- matrix(list(gry2, grDimL, g, gry1, grDist, grid::nullGrob()),nrow = 2)
  #   gt  <- gtable::gtable_matrix("di_rp_dim", mat, widths = unit(c(.25, 1,.5), "null"), heights =  unit(c(1,.25), "null"),respect = TRUE)
  # }
  #
  # if(plotDimensions&!plotDimensionLegend&!plotMeasures&unthresholded&plotRadiusRRbar){
  #   mat <- matrix(list(gry2, grid::nullGrob(),g, gry1, grDist, grid::nullGrob()),nrow = 2)
  #   gt  <- gtable::gtable_matrix("di_rp_dim", mat, widths = unit(c(.25, 1,.5), "null"), heights =  unit(c(1,.25), "null"),respect = TRUE)
  # }
  #
  # if(plotDimensions&!plotMeasures&!plotDimensionLegend&unthresholded&!plotRadiusRRbar){
  #   mat <- matrix(list(gry2, grid::nullGrob(),g, gry1),nrow = 2)
  #   gt  <- gtable::gtable_matrix("di_rp_dim", mat, widths = unit(c(.25, 1), "null"), heights =  unit(c(1,.25), "null"),respect = TRUE)
  # }
  #
  # if(plotDimensions&!plotDimensionLegend&!plotMeasures&!unthresholded){
  #   mat <- matrix(list(gry2, grid::nullGrob(),g, gry1),nrow = 2)
  #   gt  <- gtable::gtable_matrix("bi_rp_dim", mat, widths = unit(c(.25, 1), "null"), heights =  unit(c(1, .25), "null"),respect = TRUE)
  # }
  #
  # if(plotDimensions&!plotDimensionLegend&plotMeasures&unthresholded&plotRadiusRRbar){
  #   mat <- matrix(list(grA, grid::nullGrob(), gry2, grid::nullGrob(),g, gry1, grDist, grid::nullGrob()),nrow = 2)
  #   gt  <- gtable::gtable_matrix("di_rp_dim_meas", mat, widths = unit(c(.35,.25, 1,.5), "null"), heights =  unit(c(1,.25), "null"),respect = TRUE)
  # }
  #
  # if(plotDimensions&!plotDimensionLegend&plotMeasures&unthresholded&!plotRadiusRRbar){
  #   mat <- matrix(list(grA, grid::nullGrob(), gry2, grid::nullGrob(),g, gry1),nrow = 2)
  #   gt  <- gtable::gtable_matrix("di_rp_dim_meas", mat, widths = unit(c(.35,.25, 1), "null"), heights =  unit(c(1,.25), "null"),respect = TRUE)
  # }
  #
  # if(plotDimensions&!plotDimensionLegend&plotMeasures&!unthresholded){
  #   mat <- matrix(list(grA, grid::nullGrob(), gry2, grid::nullGrob(),g, gry1),nrow = 2)
  #   gt  <- gtable::gtable_matrix("bi_rp_dim_meas", mat, widths = unit(c(.35,.25, 1), "null"), heights =  unit(c(1,.25), "null"),respect = TRUE)
  # }
  #
  # if(!plotDimensions&!plotDimensionLegend&plotMeasures&unthresholded&plotRadiusRRbar){
  #   mat <- matrix(list(grA, g, grDist),nrow = 1)
  #   gt  <- gtable::gtable_matrix("di_rp_meas", mat, widths = unit(c(.35, 1,.5), "null"), heights =  unit(c(1), "null"),respect = TRUE)
  # }
  #
  # if(!plotDimensions&!plotDimensionLegend&plotMeasures&unthresholded&!plotRadiusRRbar){
  #   mat <- matrix(list(grA, g),nrow = 1)
  #   gt  <- gtable::gtable_matrix("di_rp_meas", mat, widths = unit(c(.35, 1), "null"), heights =  unit(c(1), "null"),respect = TRUE)
  # }
  #
  # if(!plotDimensions&!plotDimensionLegend&plotMeasures&!unthresholded){
  #   mat <- matrix(list(grA, g),nrow = 1)
  #   gt  <- gtable::gtable_matrix("bi_rp_meas", mat, widths = unit(c(.35, 1), "null"), heights =  unit(c(1), "null"),respect = TRUE)
  # }
  #
  # if(!plotDimensions&!plotDimensionLegend&!plotDimensionLegend&!plotMeasures&unthresholded&plotRadiusRRbar){
  #   mat <- matrix(list(g, grDist),nrow = 1)
  #   gt  <- gtable::gtable_matrix("di_rp", mat, widths = unit(c(1,.5), "null"), heights =  unit(c(1), "null"),respect = TRUE)
  # }
  #
  # if(!plotDimensions&!plotDimensionLegend&!plotMeasures&unthresholded&!plotRadiusRRbar){
  #   mat <- matrix(list(g),nrow = 1)
  #   gt  <- gtable::gtable_matrix("di_rp", mat, widths = unit(c(1), "null"), heights =  unit(c(1), "null"),respect = TRUE)
  # }
  #
  # if(!plotDimensions&!plotDimensionLegend&!plotMeasures&!unthresholded){
  #   mat <- matrix(list(g),nrow = 1)
  #   gt  <- gtable::gtable_matrix("bi_rp", mat, widths = unit(c(1), "null"), heights =  unit(c(1), "null"),respect = TRUE)
  # }

  if(reduced){
    title <- paste(title,"(coursegrained matrix)")
  }

  if(nchar(title)>0){
    grT <- ggplot2::ggplot(data.frame(x=1,y=1)) +
      ggplot2::geom_text(ggplot2::aes_(x=~x,y=~y), label=title) +
      theme_void() +
      theme(plot.margin = margin(0,0,0,0, unit = "pt"))
    gt  <- gtable::gtable_add_rows(x = gt, heights =  unit(c(.1),"null"), pos=0)
    l <- sum(c(plotDimensions,plotMeasures))+1
    gt  <- gtable::gtable_add_grob(x = gt, grobs = ggplot2::ggplotGrob(grT), name = "Title", t=1, l=l)
  }

  g <- gtable::gtable_add_padding(gt, unit(5, "pt"))


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
#' Calculate the maximum possible number of recurrent points in a recurrence matrix.
#'
#' This function can take into account the presence of a `theiler` window, that is the points in the window will be excluded from the calculation. For example, some scholars will exclude the main diagonal from the calculation of the recurrence rate.
#'
#' @param RM A Matrix object
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
rp_size <- function(RM, AUTO=NULL, theiler = NULL){

  if(is.null(AUTO)){
      AUTO <- Matrix::isSymmetric(RM)
    }

  if(is.null(theiler)){
    if(is.null(attributes(RM)$theiler)){
      if(AUTO){
        theiler <- 0
      } else {
        theiler <- 1
      }
    } else {
      theiler <- attributes(RM)$theiler
    }
  }

  R_N <- cumprod(dim(RM))[2]
  minDiag <- 0

  if(is.na(theiler)){
    if(AUTO){
      theiler <- 1
    } else {
      theiler <- 0
    }
  }

  if(length(theiler)==1){
    if(theiler==1){
      minDiag <- length(Matrix::diag(RM))
    }
    if(theiler>1){
      minDiag <- sum(rep(length(Matrix::diag(RM)),length(seq(-theiler,theiler))) - abs(seq(-theiler,theiler)))
    }
  }

  if(length(theiler)==2){
    minDiag <- sum(rep(length(Matrix::diag(RM)),length(seq(min(theiler),max(theiler)))) - abs(seq(min(theiler),max(theiler))))
  }

  if(length(theiler)>2){
    minDiag <- sum(rep(length(Matrix::diag(RM)),length(theiler)) - abs(theiler))
  }

  return(list(rp_size_total = R_N, rp_size_theiler = R_N - minDiag))
  #cumprod(dim(mat))[2] - ifelse((includeDiag&theiler==0),length(Matrix::diag(mat)),ifelse(theiler>0,Matrix::nnzero(Matrix::band(mat,-theiler,theiler)),0)))
}


#' Empty results vector
#'
#' @return an empty rp_measures
#' @keywords internal
#' @export
#'
rp_empty <- function(){
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

#' rp_calc
#'
#' @inheritParams rp_measures
#'
#' @return CRQA measures and matrices of line distributions (if requested)
#' @export
#' @keywords internal
#'
rp_calc <- function(RM,
                    emRad = NULL,
                    DLmin = 2,
                    VLmin = 2,
                    HLmin = 2,
                    DLmax = length(Matrix::diag(RM)),
                    VLmax = length(Matrix::diag(RM)),
                    HLmax = length(Matrix::diag(RM)),
                    theiler   = NA,
                    AUTO      = NULL,
                    includeDiagonal = NA,
                    chromatic = FALSE,
                    anisotropyHV = FALSE,
                    asymmetryUL = FALSE,
                    recurrenceTimes = FALSE,
                    matrices  = FALSE){



  RM <- rp_checkfix(RM, checkAUTO = TRUE, fixAUTO = TRUE)

  if(is.null(AUTO)){
    AUTO <- attributes(RM)$AUTO
  }

  if(is.na(theiler)){
    if(is.na(attributes(RM)$theiler)){
      RM <- setTheiler(RM, theiler = theiler)
    }
    theiler <- attributes(RM)$theiler
  }

  if(!is.na(includeDiagonal)){
    AUTO <- includeDiagonal
  }
  recmatsize <- rp_size(RM, AUTO = AUTO, theiler = theiler)

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

  # Total nr. recurrent points
  RP_N <- Matrix::nnzero(RM, na.counted = FALSE)

  #Proportion recurrence / Recurrence Rate
  RR <- RP_N/recmatsize$rp_size_theiler

  if(length(RR)==0){RR<-0}

  if(RR==1){
    warning("Everything is recurring!\nReturning empty vector")
    return(rp_empty())
  }

  lineMeasures_global <- rp_calc_lineMeasures(RM = RM,
                                              RP_N = RP_N,
                                              DLmin = DLmin, DLmax = DLmax,
                                              VLmin = VLmin, VLmax = VLmax,
                                              HLmin = HLmin, HLmax = HLmax,
                                              theiler = theiler,
                                              AUTO = AUTO,
                                              chromatic = chromatic,
                                              matrices = matrices)

  if(matrices){
    lineMatrices_global <- lineMeasures_global$crqaMatrices
    lineMeasures_global <- lineMeasures_global$crqaMeasures
  }

  # Singularities
  SING <- rp_lineDist(RM,
                      DLmin = 1, DLmax = DLmax,
                      VLmin = 1, VLmax = VLmax,
                      HLmin = 1, HLmax = HLmax,
                      theiler = theiler, AUTO = AUTO)

  # table(SING_N$verticals.dist)
  # table(SING_N$horizontals.dist)
  SING_N <- table(SING$diagonals.dist)[1]
  SING_rate <- SING_N / RP_N

  # H/V Anisotropy ratios
  if(anisotropyHV){
    out_hv_ani <- data.frame(
      Nlines_ani   = ((lineMeasures_global$N_hlp - lineMeasures_global$N_vlp)  / (lineMeasures_global$N_hlp  + lineMeasures_global$N_vlp))%00%NA,
      LAM_ani      = (lineMeasures_global$LAM_hl - lineMeasures_global$LAM_vl) / (lineMeasures_global$LAM_hl + lineMeasures_global$LAM_vl),
      MEAN_hvl_ani = (lineMeasures_global$TT_hl  - lineMeasures_global$TT_vl)  / (lineMeasures_global$TT_hl  + lineMeasures_global$TT_vl),
      MAX_hvl_ani  = (lineMeasures_global$MAX_hl - lineMeasures_global$MAX_vl) / (lineMeasures_global$MAX_hl + lineMeasures_global$MAX_vl),
      ENT_hvl_ani  = (lineMeasures_global$ENT_hl - lineMeasures_global$ENT_vl) / (lineMeasures_global$ENT_hl + lineMeasures_global$ENT_vl)
    )
  } else {
    out_hv_ani <- NA
  }


  # U/L Anisotropy ratios
  if(asymmetryUL){
    RMu <- RM
    RMu[lower.tri(RMu)] <- 0

    recmatsize_u <- rp_size(RMu, AUTO = AUTO, theiler = theiler)

    lineMeasures_upper <- rp_calc_lineMeasures(RM = RMu,
                                               RP_N = RP_N,
                                               DLmin = DLmin, DLmax = DLmax,
                                               VLmin = VLmin, VLmax = VLmax,
                                               HLmin = HLmin, HLmax = HLmax,
                                               theiler = theiler,
                                               AUTO = AUTO,
                                               chromatic = chromatic,
                                               matrices = FALSE)

    # Total nr. recurrent points in upper
    lineMeasures_upper$RP_N <- Matrix::nnzero(RMu, na.counted = FALSE)

    #Proportion recurrence / Recurrence Rate un upper
    lineMeasures_upper$RR <- lineMeasures_upper$RP_N / recmatsize_u$rp_size_theiler

    if(length(lineMeasures_upper$RR)==0){lineMeasures_upper$RR<-0}


    SING_u <- rp_lineDist(RMu,
                          DLmin = 1, DLmax = DLmax,
                          VLmin = 1, VLmax = VLmax,
                          HLmin = 1, HLmax = HLmax,
                          theiler = theiler, AUTO = AUTO)

    lineMeasures_upper$SING_N  <-  table(SING_u$diagonals.dist)[1]

    rm(RMu, recmatsize_u, SING_u)

    RMl <- RM
    RMl[upper.tri(RMl)] <- 0

    recmatsize_l <- rp_size(RMl, AUTO = AUTO, theiler = theiler)

    lineMeasures_lower <- rp_calc_lineMeasures(RM = RMl,
                                               RP_N = RP_N,
                                               DLmin = DLmin, DLmax = DLmax,
                                               VLmin = VLmin, VLmax = VLmax,
                                               HLmin = HLmin, HLmax = HLmax,
                                               theiler = theiler,
                                               AUTO = AUTO,
                                               chromatic = chromatic,
                                               matrices = FALSE)


    # Total nr. recurrent points in lower
    lineMeasures_lower$RP_N <- Matrix::nnzero(RMl, na.counted = FALSE)

    #Proportion recurrence / Recurrence Rate in lower
    lineMeasures_lower$RR <- lineMeasures_lower$RP_N / recmatsize_l$rp_size_theiler

    if(length(lineMeasures_lower$RR)==0){lineMeasures_lower$RR<-0}

    SING_l <- rp_lineDist(RMl,
                          DLmin = 1, DLmax = DLmax,
                          VLmin = 1, VLmax = VLmax,
                          HLmin = 1, HLmax = HLmax,
                          theiler = theiler, AUTO = AUTO)

    lineMeasures_lower$SING_N  <-  table(SING_l$diagonals.dist)[1]
    rm(RMl,SING_l)

    # RATIOs
    out_ul_ani <- data.frame(
      Npoints_ul_ani = (lineMeasures_upper$RP_N   - lineMeasures_lower$RP_N)   / (lineMeasures_upper$RP_N   + lineMeasures_lower$RP_N),
      NDlines_ul_ani = ((lineMeasures_upper$N_dlp - lineMeasures_lower$N_dlp)  / (lineMeasures_upper$N_dlp  + lineMeasures_lower$N_dlp))%00%NA,
      NHlines_ul_ani = ((lineMeasures_upper$N_hlp - lineMeasures_lower$N_hlp)  / (lineMeasures_upper$N_hlp  + lineMeasures_lower$N_hlp))%00%NA,
      NVlines_ul_ani = ((lineMeasures_upper$N_vlp - lineMeasures_lower$N_vlp)  / (lineMeasures_upper$N_vlp  + lineMeasures_lower$N_vlp))%00%NA,
      RR_ul_ani      = ((lineMeasures_upper$RR    - lineMeasures_lower$RR)     / (lineMeasures_upper$RR     + lineMeasures_lower$RR))%00%NA,
      SING_N_ul_ani  = (lineMeasures_upper$SING_N - lineMeasures_lower$SING_N) / (lineMeasures_upper$SING_N + lineMeasures_lower$SING_N),
      DIV_ul_ani     = (lineMeasures_upper$DIV_dl - lineMeasures_lower$DIV_dl) / (lineMeasures_upper$DIV_dl + lineMeasures_lower$DIV_dl),
      REP_ul_ani     = (lineMeasures_upper$REP_av - lineMeasures_lower$REP_av) / (lineMeasures_upper$REP_av + lineMeasures_lower$REP_av),
      DET_ul_ani     = (lineMeasures_upper$DET    - lineMeasures_lower$DET)    / (lineMeasures_upper$DET    + lineMeasures_lower$DET),
      LAM_hl_ul_ani  = (lineMeasures_upper$LAM_hl - lineMeasures_lower$LAM_hl) / (lineMeasures_upper$LAM_hl + lineMeasures_lower$LAM_hl),
      LAM_vl_ul_ani  = (lineMeasures_upper$LAM_vl - lineMeasures_lower$LAM_vl) / (lineMeasures_upper$LAM_vl + lineMeasures_lower$LAM_vl),
      MEAN_dl_ul_ani = (lineMeasures_upper$MEAN_dl- lineMeasures_lower$MEAN_dl)/ (lineMeasures_upper$MEAN_dl+ lineMeasures_lower$MEAN_dl),
      MEAN_hl_ul_ani = (lineMeasures_upper$TT_hl  - lineMeasures_lower$TT_hl)  / (lineMeasures_upper$TT_hl  + lineMeasures_lower$TT_hl),
      MEAN_vl_ul_ani = (lineMeasures_upper$TT_vl  - lineMeasures_lower$TT_vl)  / (lineMeasures_upper$TT_vl  + lineMeasures_lower$TT_vl),
      MAX_dl_ul_ani  = (lineMeasures_upper$MAX_dl - lineMeasures_lower$MAX_dl) / (lineMeasures_upper$MAX_dl + lineMeasures_lower$MAX_dl),
      MAX_hl_ul_ani  = (lineMeasures_upper$MAX_hl - lineMeasures_lower$MAX_hl) / (lineMeasures_upper$MAX_hl + lineMeasures_lower$MAX_hl),
      MAX_vl_ul_ani  = (lineMeasures_upper$MAX_vl - lineMeasures_lower$MAX_vl) / (lineMeasures_upper$MAX_vl + lineMeasures_lower$MAX_vl),
      ENT_dl_ul_ani  = (lineMeasures_upper$ENT_dl - lineMeasures_lower$ENT_dl) / (lineMeasures_upper$ENT_dl + lineMeasures_lower$ENT_dl),
      ENT_hl_ul_ani  = (lineMeasures_upper$ENT_hl - lineMeasures_lower$ENT_hl) / (lineMeasures_upper$ENT_hl + lineMeasures_lower$ENT_hl),
      ENT_vl_ul_ani  = (lineMeasures_upper$ENT_vl - lineMeasures_lower$ENT_vl) / (lineMeasures_upper$ENT_vl + lineMeasures_lower$ENT_vl)
    )


  }

  #Output
  out_global <- data.frame(emRad      = emRad,
                           RP_max     = recmatsize$rp_size_theiler,
                           RR         = RR,
                           SING_N     = SING_N,
                           SING_rate  = SING_rate)


  if(asymmetryUL){
    out_UL <- data.frame(upper.tri = lineMeasures_upper,
                         lower.tri = lineMeasures_lower,
                         ratios.ul = out_ul_ani)
  } else {
    out_UL <- NA
  }

  if(anisotropyHV){
    out_HL <- rbind(ratios.hl = out_hv_ani)
  } else {
    out_HL <- NA
  }

  out <- as.data.frame(list(out_global, lineMeasures_global, out_UL, out_HL)[c(TRUE, TRUE, asymmetryUL, anisotropyHV)])

  if(matrices){
    return(list(crqaMeasures = out,
                crqaMatrices = lineMatrices_global)
    )
  } else {
    return(out)
  }
}


rp_calc_lineMeasures <- function(RM,
                                 RP_N,
                                 DLmin = 2,
                                 VLmin = 2,
                                 HLmin = 2,
                                 DLmax = length(Matrix::diag(RM)),
                                 VLmax = length(Matrix::diag(RM)),
                                 HLmax = length(Matrix::diag(RM)),
                                 d         = NULL,
                                 theiler   = NULL,
                                 invert    = FALSE,
                                 AUTO      = NULL,
                                 chromatic = FALSE,
                                 matrices  = FALSE){

  #Get line segments
  # if(Matrix::nnzero(RM)>0)
  lineSegments <- rp_lineDist(RM,
                              DLmin = DLmin, DLmax = DLmax,
                              VLmin = VLmin, VLmax = VLmax,
                              HLmin = HLmin, HLmax = HLmax,
                              d = d, theiler = theiler,
                              invert = invert, AUTO = AUTO,
                              chromatic = chromatic,
                              matrices = matrices)

  dlines <- lineSegments$diagonals.dist%00%0
  vlines <- lineSegments$verticals.dist%00%0
  hlines <- lineSegments$horizontals.dist%00%0

  #Frequency tables of line lengths
  freq_dl <- table(dlines)
  freq_vl <- table(vlines)
  freq_hl <- table(hlines)
  freq_hv <- table(c(hlines,vlines))

  freqvec_dl <- as.numeric(names(freq_dl))
  freqvec_vl <- as.numeric(names(freq_vl))
  freqvec_hl <- as.numeric(names(freq_hl))
  freqvec_hv <- as.numeric(names(freq_hv))

  # Number of lines
  N_dl <- sum(freq_dl, na.rm = TRUE)%00%0
  N_vl <- sum(freq_vl, na.rm = TRUE)%00%0
  N_hl <- sum(freq_hl, na.rm = TRUE)%00%0
  N_hv <- sum(freq_hv, na.rm = TRUE)%00%0

  #Number of recurrent points on diagonal, vertical and horizontal lines
  N_dlp <- sum(freqvec_dl*freq_dl, na.rm = TRUE)
  N_vlp <- sum(freqvec_vl*freq_vl, na.rm = TRUE)
  N_hlp <- sum(freqvec_hl*freq_hl, na.rm = TRUE)
  N_hvp <- sum(freqvec_hv*freq_hv, na.rm = TRUE)

  #Determinism / Horizontal and Vertical Laminarity
  DET    <- N_dlp/RP_N
  LAM_vl <- N_vlp/RP_N
  LAM_hl <- N_hlp/RP_N
  LAM_hv <- N_hvp/(RP_N*2)

  #Array of probabilities that a certain line length will occur (all >1)
  P_dl <- freq_dl/N_dl
  P_vl <- freq_vl/N_vl
  P_hl <- freq_hl/N_hl
  P_hv <- freq_hv/N_hv

  #Entropy of line length distributions
  ENT_dl <- -1 * sum(P_dl * log(P_dl))
  ENT_vl <- -1 * sum(P_vl * log(P_vl))
  ENT_hl <- -1 * sum(P_hl * log(P_hl))
  ENT_hv <- -1 * sum(P_hv * log(P_hv))

  #Relative Entropy (Entropy / Max entropy)
  ENTrel_dl = ENT_dl/(-1 * log(1/DLmax))
  ENTrel_vl = ENT_vl/(-1 * log(1/VLmax))
  ENTrel_hl = ENT_hl/(-1 * log(1/HLmax))
  ENTrel_hv = ENT_hv/(-1 * log(1/max(c(HLmax,VLmax), na.rm = TRUE)))

  #Meanline
  MEAN_dl = mean(dlines, na.rm = TRUE)%00%0
  MEAN_vl = mean(vlines, na.rm = TRUE)%00%0
  MEAN_hl = mean(hlines, na.rm = TRUE)%00%0
  MEAN_hv = mean(c(hlines,vlines), na.rm = TRUE)%00%0

  #Maxline
  MAX_dl = max(freqvec_dl, na.rm = TRUE)%00%0
  MAX_vl = max(freqvec_vl, na.rm = TRUE)%00%0
  MAX_hl = max(freqvec_hl, na.rm = TRUE)%00%0
  MAX_hv = max(freqvec_hv, na.rm = TRUE)%00%0

  # REPetetiveness
  REP_av  <- (N_hlp+N_vlp) / N_dlp
  REP_hl  <-  N_hlp/N_dlp
  REP_vl  <-  N_vlp/N_dlp

  #Coefficient of determination
  CoV_dl = stats::sd(dlines)/mean(dlines)
  CoV_vl = stats::sd(vlines)/mean(vlines)
  CoV_hl = stats::sd(hlines)/mean(hlines)
  CoV_hv = stats::sd(c(hlines,vlines))/mean(c(hlines,vlines))

  #Divergence
  DIV_dl = 1/MAX_dl
  DIV_vl = 1/MAX_vl
  DIV_hl = 1/MAX_hl
  DIV_hv = 1/MAX_hv

  out <- data.frame(
    RP_N      = RP_N,
    DIV_dl    = DIV_dl,
    REP_av    = REP_av,
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
    N_hlp     = N_hlp,
    N_hl      = N_hl,
    LAM_hl    = LAM_hl,
    TT_hl     = MEAN_hl,
    MAX_hl    = MAX_hl,
    ENT_hl    = ENT_hl,
    ENTrel_hl = ENTrel_hl,
    CoV_hl    = CoV_hl,
    REP_hl    = REP_hl,
    N_hvp     = N_hvp,
    N_hv      = N_hv,
    LAM_hv    = LAM_hv,
    TT_hv     = MEAN_hv,
    MAX_hv    = MAX_hv,
    ENT_hv    = ENT_hv,
    ENTrel_hv = ENTrel_hv,
    CoV_hv    = CoV_hv)

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
rp_prep <- function(RP,
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

  out<-rp_calc(RP,
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
      outLo <- rp_calc(Matrix::tril(RP,-1),
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

      outUp <- rp_calc(Matrix::triu(RP,-1),
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
                             lower = rp_empty(),
                             upper = rp_empty())
    }
  }
  return(out)
}


#' False Nearest Neighbours
#'
#' Search for FNN to get an optimal Embedding Dimension using by using [nonlinearTseries::findAllNeighbours()] in a loop.
#'
#' @inheritParams est_parameters
#' @param radius Size of the neighbourhood: Every point smaller than the radius will be considered a near neighbour, see [nonlinearTseries::findAllNeighbours()] (default = `sd(y)/10`).
#' @param number.boxes Integer representing number of boxes to to speed up neighbour search, if `NULL` an optimal number will be chosen [nonlinearTseries::findAllNeighbours()] (default = `NULL`).
#'
#' @return FNN curve
#' @export
#'
fnn <- function(y, emLag = 1, maxDim = 10, radius = sd(y)/10, number.boxes = NULL){

  out <- matrix(NA, nrow=maxDim, ncol = 1)

  # for(emL in 1:length(emLag)){
  for(emD in 1: maxDim){

    yy <- nonlinearTseries::buildTakens(time.series = y, embedding.dim = emD, time.lag = emLag)
    nn <- nonlinearTseries::findAllNeighbours(as.matrix(yy[,1:emD]), radius = radius, number.boxes = number.boxes)
    out[emD,1] <- sum(lengths(nn))

  }
  # }
  return(out)
}
