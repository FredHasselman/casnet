# package casnet ----
#
# Massively Parallel Recurrence Analysis ----
#
# rqa functions


#' Massively Parallel RQA analysis
#'
#' @description Calculate (C)RQA measures without creating a recurrence matrix. Can handle very large time series and requires package [parallel] to be installed.
#'
#' @inheritParams rp
#' @inheritParams rp_measures
#'
#' @export
#'
#' @return
#'
#' @examples
#'
rqa <- function(y1, y2 = NULL,
                emDim = 1,
                emLag = 1,
                emRad = NULL,
                theiler = NA,
                includeDiagonal = NA,
                AUTO = NULL,
                DLmin = 2,
                VLmin = 2,
                HLmin = 2,
                DLmax = NROW(y1),
                VLmax = NROW(y1),
                HLmax = NROW(y1),
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


  if(is.null(AUTO)){
    stop("Auto-RQA or Cross-RQA? (Provide a value for AUTO)")
  }

  if(!rescaleDist%in%"none"){
    warning("Cannot rescale distance matrix, please rescale the time series to standardeviation (z-score) or unit scale (min/max)")
  }

  if(is.na(chromatic)){
    # if(!is.null(attr(RM,"chromatic"))){
    #   chromatic <- attr(RM,"chromatic")
    # } else {
    chromatic <- FALSE
  }
  #}

  if(chromatic){
    if(abs(diff(range(y1)))==1){
      message("There is only 1 ctategory, chromatic (C)RQA is not sensible. Setting chromatic to FALSE...")
      chromatic <- FALSE
    }
  }

  if(any(!doEmbed,chromatic)){
    emDim <- 1
    emLag <- 0
    if(chromatic){
      emRad <- 0
    }
  }

  tmp <- rqa_getSeries(y1 = y1,
                y2 = y2,
                emDim = emDim,
                emLag = emLag,
                chromatic = chromatic,
                doEmbed = doEmbed)

  et1 <- tmp$et1
  et2 <- tmp$et2

  chromaDims  <- attributes(tmp)$chromaDims
  chromaNames <- attributes(tmp)$chromaNames

  rm(tmp)

  # Estimate radius?
  if(is.na(emRad%00%NA)){
    # if(!is.null(attributes(RM)$emRad)){
    #   emRad <- attributes(RM)$emRad
    # } else {
      # Check for attributes
      if(is.na(targetValue)){
        targetValue <- .05
      }
      if(is.null(emDim)){
        emDim <- 1
      }
      if(is.null(emLag)){
        emLag <- 1
      }

      emRad <- est_radius_rqa(y1 = et1, y2 = et2, AUTO = AUTO, emDim = emDim, emLag = emLag, targetValue = targetValue, radiusOnFail = "minimum", theiler = theiler, method = method, silent = silent)

      if(emRad$Converged){
        emRad <- emRad$Radius
      } else {
        emRad <- stats::sd(RM,na.rm = TRUE)
      }
  }

  dist_method <- return_error(proxy::pr_DB$get_entry(method))
  if("error"%in%class(dist_method$value)){
    stop("Unknown distance metric!\nUse proxy::pr_DB$get_entries() to see a list of valid options.")
  } else {

    RQAout <- rqa_measures(y1 = et1, y2 = et2, emRad = emRad, DLmin = DLmin, VLmin = VLmin, HLmin = HLmin, DLmax = DLmax, VLmax = VLmax, HLmax = HLmax, AUTO = AUTO, theiler = theiler, method = method, silent = silent)

  }

  #require(parallel)
#
#   if(to.sparse){
#     attributes(dmat)$emDims1  <- et1
#     attributes(dmat)$emDims2  <- et2
#     attributes(dmat)$emDims1.name <- colnames(y1)
#     attributes(dmat)$emDims2.name <- colnames(y2)
#     attributes(dmat)$embedded <- doEmbed
#     attributes(dmat)$emLag <- emLag
#     attributes(dmat)$emDim <- emDim
#     attributes(dmat)$emRad <- emRad%00%NA
#     attributes(dmat)$measures <- rpOut
#     attributes(dmat)$weighted <- weighted
#     attributes(dmat)$weightedBy <- weightedBy
#     attributes(dmat)$chromatic <- chromatic
#     attributes(dmat)$chromaNames <- chromaNames
#     attributes(dmat)$chromaDims <- chromaDims
#   } else {
#     attr(dmat,"emDims1") <- et1
#     attr(dmat,"emDims2") <- et2
#     attr(dmat,"emDims1.name") <- colnames(y1)
#     attr(dmat,"emDims2.name") <- colnames(y2)
#     attr(dmat,"weighted") <- weighted
#     attr(dmat,"embedded") <- doEmbed
#     attr(dmat,"emLag") <- emLag
#     attr(dmat,"emDim") <- emDim
#     attr(dmat,"emRad") <- emRad%00%NA
#     attr(dmat,"measures") <- rpOut
#     attr(dmat,"weighted") <- weighted
#     attr(dmat,"weightedBy") <- weightedBy
#     attr(dmat,"chromatic") <- chromatic
#     attr(dmat,"chromaNames") <- chromaNames
#     attr(dmat,"chromaDims") <- chromaDims
#   }
}


#' Get embedded series for rp() and rqa()
#'
#' @param y1 y1
#' @param y2 y2
#' @param emDim emDim
#' @param emLag emLag
#' @param chromatic chromatic
#' @param doEmbed doEmbed
#'
#' @return embedded series
#'
#' @export
#'
#' @keywords internal
#'
rqa_getSeries <- function(y1, y2 = NULL,
                          emDim = 1,
                          emLag = 1,
                          chromatic = FALSE,
                          doEmbed = TRUE){

  y1 <- float::as.float(y1)

  if(is.null(y2)){
    y2 <- y1
    attributes(y2) <- attributes(y1)
  }

  atlist <- attributes(y1)
  if(all(any(names(atlist) %in% c("emDims1","emDims2","emDims1.name","emDims2.name")),is.matrix(y1))){
    stop("Input is a recurrence matrix created by 'rp()'. Please provide time series (numeric vectors)")
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

  out <- list(et1 = ts_embed(y1, emDim=emDim, emLag=emLag, silent = silent, doEmbed = doEmbed),
              et2 = ts_embed(y2, emDim=emDim, emLag=emLag, silent = silent, doEmbed = doEmbed))

  attr(out,"chromatic")   <- chromatic
  attr(out,"chromaNames") <- chromaNames
  attr(out,"chromaDims")  <- chromaDims

  return(out)
}


#' Fast rqa
#'
#' @param y1 y1
#' @param y2 y2
#' @param index row
#' @param theiler theiler
#' @param emRad threshold
#' @param symmetrical symmetrical
#' @param diagonal diagonals
#' @param method distance method
#'
#' @return
#'
#' @export
#'
#' @keywords internal
#'
rqa_fast <- function(y1, y2 = NA, index, theiler, emRad, symmetrical = TRUE, diagonal = FALSE, method = "Euclidean"){

  checkPkg("bit")
  checkPkg("float")

  Xlength <- NROW(y1)
  if(all(is.na(y2))){
    Ylength <- Xlength
  } else {
    Ylength <- NROW(y2)
  }

  # if(all(symmetrical,theiler==0,!diagonal)){
  #   symmetrical <- FALSE
  # }

  horizontal <- FALSE
  if(!symmetrical&!diagonal){
    horizontal <- TRUE
  }

  inds <- mat_ind(Xlength = Xlength, Ylength = Ylength, index = index, diagonal = diagonal)

  # if(all(is.na(y2))){
  #   Y <- matrix(y1[inds$c,], ncol = 2)
  # } else {
  #   Y <- matrix(y2[inds$c,], ncol = 2)
  # }
  # X <- matrix(y1[inds$r,], ncol = 2)

  if(all(is.na(y2))){
    Y <- matrix(y1[inds$c,], nrow = NROW(inds$c))
  } else {
    Y <- matrix(y2[inds$c,], nrow = NROW(inds$c))
  }
  X <- matrix(y1[inds$r,], nrow = NROW(inds$r))
  rm(inds, y1, y2)

  if(!diagonal){
    if((index+theiler)<Ylength){
      if(theiler>0){
        Y[index+(theiler-1),] <- NA
      }
    }
  }

  if(horizontal){
    D   <- float::as.float(proxy::dist(X, Y, pairwise = diagonal, diag = TRUE, method = method))
  } else {
    D   <- float::as.float(proxy::dist(Y, X, pairwise = diagonal, diag = TRUE, method = method))
  }
  rm(X,Y)

  minDist  <- min(D[D>float::as.float(0)], na.rm = TRUE)
  minDistN <- sum((D==float::as.float(minDist)), na.rm = TRUE)
  maxDist  <- max(D[D>float::as.float(0)], na.rm = TRUE)
  meanDist <- float::colMeans(D[D>float::as.float(0)], na.rm = TRUE)
  idx <- bit::as.bit(D<=float::as.float(emRad))
  rm(D)

  attr(idx,"minDist")  <- minDist
  attr(idx,"minDistN") <- minDistN
  attr(idx,"maxDist")  <- maxDist
  attr(idx,"meanDist") <- meanDist
  #return(bit::as.bit(idx))
  return(idx)
}


# # Create Y based on row
# if(!diagonals){
#
#   if(any(symmetrical,all(is.na(Y)))){
#     Yy <- X[row,]
#     Xx <- X[(row+1):NROW(X),]
#   } else {
#     Yy <- Y[row,]
#     Xx <- X
#   }
#
# } else {# diagonals
#
#   symmetrical <- FALSE
#
#   if(row == 0){
#     Xx <- X[seq(1,NROW(X)),]
#     Yy <- Y[seq(1,NROW(Y)),]
#   }
#   if(row > 0){
#     Xx <- X[seq(1,(NROW(X)-row)),]
#     Yy <- Y[seq((row+1),NROW(Y)),]
#   }
#   if(row < 0){
#     row <- abs(row)
#     Xx <- Y[seq((row+1),NROW(Y)),]
#     Yy <- X[seq(1,(NROW(X)-row)),]
#   }
# }

#rm(X,Y)

#' Fast line dist
#'
#' @param idx idx
#'
#' @return
#' @export
#'
#' @examples
rqa_lineDist <- function(idx){
  sort(which(diff(c(0,as.numeric(idx),0))==-1)-which(diff(c(0,as.numeric(idx),0))==1))
}

#' Fast RQA
#'
#' @description Calculates RQA measures without building the recurrence matrix. Requires package [parallel].
#'
#' @inheritParams rp
#' @inheritParams rp_measures
#'
#' @return
#'
#' @export
#'
#' @examples
#'
#'
#'
rqa_measures <- function(y1,
                         y2 = NA,
                         emRad = NA,
                         DLmin = 2,
                         VLmin = 2,
                         HLmin = 2,
                         DLmax = NROW(y1),
                         VLmax = NROW(y1),
                         HLmax = NROW(y1),
                         AUTO      = NULL,
                         theiler   = NA,
                         method = "Euclidean",
                         chromatic = FALSE,
                         anisotropyHV = FALSE,
                         asymmetryUL = FALSE,
                         returnUL = FALSE,
                         recurrenceTimes = FALSE,
                         distributions = FALSE,
                         silent = TRUE){

  if(is.na(theiler%00%NA)){
    if(AUTO){
      theiler <- 1
      message("Default behaviour for Auto-RQA is to exclude the line of incidence (main diagonal), set theiler = 0 to include it.")
    } else {
      theiler <- 0
    }
  }

  if(AUTO){
    y2 <- NA
  }

  lineMeasures_global <- rqa_calc(y1 = y1, y2 = y2,
                                  emRad = emRad,
                                  DLmin = DLmin, DLmax = DLmax,
                                  VLmin = VLmin, VLmax = VLmax,
                                  HLmin = HLmin, HLmax = HLmax,
                                  theiler = theiler,
                                  recurrenceTimes = recurrenceTimes,
                                  AUTO = AUTO,
                                  method = method,
                                  chromatic = chromatic,
                                  anisotropyHV = anisotropyHV,
                                  asymmetryUL = asymmetryUL,
                                  returnUL = returnUL,
                                  distributions = distributions)

  if(distributions){
    lineDistributions_global <- lineMeasures_global$crqaMatrices
    lineMeasures_global      <- lineMeasures_global$crqaMeasures
  }


  # H/V Anisotropy ratios
  if(anisotropyHV){
    if(AUTO){
      message("Matrix is symmetrical so Horizontal/Vertical anisotropic ratios will be 0.")
    }
    out_hv_ani <- data.frame(
      Nlines_ani   = ((lineMeasures_global$N_hlp - lineMeasures_global$N_vlp)  / (lineMeasures_global$N_hlp  + lineMeasures_global$N_vlp))%00%NA,
      LAM_ani  = (lineMeasures_global$LAM_hl - lineMeasures_global$LAM_vl) / (lineMeasures_global$LAM_hl + lineMeasures_global$LAM_vl),
      MEAN_ani = (lineMeasures_global$TT_hl  - lineMeasures_global$TT_vl)  / (lineMeasures_global$TT_hl  + lineMeasures_global$TT_vl),
      MAX_ani  = (lineMeasures_global$MAX_hl - lineMeasures_global$MAX_vl) / (lineMeasures_global$MAX_hl + lineMeasures_global$MAX_vl),
      ENT_ani  = (lineMeasures_global$ENT_hl - lineMeasures_global$ENT_vl) / (lineMeasures_global$ENT_hl + lineMeasures_global$ENT_vl)
    )
  } else {
    out_hv_ani <- NA
  }


  # U/L Anisotropy ratios
  if(asymmetryUL){
    if(AUTO){
      message("Matrix is symmetrical so Upper/Lower assymetry ratios will be 0.")
    }

  #   recmatsize_u <- rp_size(dims = dims, AUTO = FALSE, theiler = theiler)
  #
  #   lineMeasures_upper <- rqa_calc_lineMeasures(RM = RMu,
  #                                              RP_N = RP_N,
  #                                              DLmin = DLmin, DLmax = DLmax,
  #                                              VLmin = VLmin, VLmax = VLmax,
  #                                              HLmin = HLmin, HLmax = HLmax,
  #                                              theiler = theiler,
  #                                              recurrenceTimes = recurrenceTimes,
  #                                              AUTO = AUTO,
  #                                              chromatic = chromatic,
  #                                              matrices = FALSE)
  #
  #   # Total nr. recurrent points in upper
  #   lineMeasures_upper$RP_N <- Matrix::nnzero(RMu, na.counted = FALSE)
  #
  #   #Proportion recurrence / Recurrence Rate un upper
  #   lineMeasures_upper$RR <- lineMeasures_upper$RP_N / recmatsize_u$rp_size_theiler
  #
  #   if(length(lineMeasures_upper$RR)==0){lineMeasures_upper$RR<-0}
  #
  #
  #   SING_u <- rp_lineDist(RMu,
  #                         DLmin = 1, DLmax = DLmax,
  #                         VLmin = 1, VLmax = VLmax,
  #                         HLmin = 1, HLmax = HLmax,
  #                         theiler = theiler, AUTO = AUTO)
  #
  #   lineMeasures_upper$SING_N  <-  table(SING_u$diagonals.dist)[1]
  #
  #   rm(RMu, recmatsize_u, SING_u)
  #
  #   RMl <- RM
  #   RMl[upper.tri(RMl)] <- 0
  #
  #   recmatsize_l <- rp_size(RM = RMl, AUTO = AUTO, theiler = theiler)
  #
  #   lineMeasures_lower <- rp_calc_lineMeasures(RM = RMl,
  #                                              RP_N = RP_N,
  #                                              DLmin = DLmin, DLmax = DLmax,
  #                                              VLmin = VLmin, VLmax = VLmax,
  #                                              HLmin = HLmin, HLmax = HLmax,
  #                                              theiler = theiler,
  #                                              AUTO = AUTO,
  #                                              recurrenceTimes = recurrenceTimes,
  #                                              chromatic = chromatic,
  #                                              matrices = FALSE)
  #
  #
  #   # Total nr. recurrent points in lower
  #   lineMeasures_lower$RP_N <- Matrix::nnzero(RMl, na.counted = FALSE)
  #
  #   #Proportion recurrence / Recurrence Rate in lower
  #   lineMeasures_lower$RR <- lineMeasures_lower$RP_N / recmatsize_l$rp_size_theiler
  #
  #   if(length(lineMeasures_lower$RR)==0){lineMeasures_lower$RR<-0}
  #
  #   SING_l <- rp_lineDist(RMl,
  #                         DLmin = 1, DLmax = DLmax,
  #                         VLmin = 1, VLmax = VLmax,
  #                         HLmin = 1, HLmax = HLmax,
  #                         theiler = theiler, AUTO = AUTO)
  #
  #   lineMeasures_lower$SING_N  <-  table(SING_l$diagonals.dist)[1]
  #   rm(RMl,SING_l)
  #
  #   # RATIOs
  #   out_ul_asym <- data.frame(
  #     Npoints_asym = (lineMeasures_upper$RP_N   - lineMeasures_lower$RP_N)/ (lineMeasures_upper$RP_N  + lineMeasures_lower$RP_N),
  #     NDlines_asym = ((lineMeasures_upper$N_dl - lineMeasures_lower$N_dl) / (lineMeasures_upper$N_dl  + lineMeasures_lower$N_dl))%00%NA,
  #     NHlines_asym = ((lineMeasures_upper$N_hl - lineMeasures_lower$N_hl) / (lineMeasures_upper$N_hl  + lineMeasures_lower$N_hl))%00%NA,
  #     NVlines_asym = ((lineMeasures_upper$N_vl - lineMeasures_lower$N_vl) / (lineMeasures_upper$N_vl  + lineMeasures_lower$N_vl))%00%NA,
  #     NDpoints_asym = ((lineMeasures_upper$N_dlp - lineMeasures_lower$N_dlp)/ (lineMeasures_upper$N_dlp  + lineMeasures_lower$N_dlp))%00%NA,
  #     NHpoints_asym = ((lineMeasures_upper$N_hlp - lineMeasures_lower$N_hlp)/ (lineMeasures_upper$N_hlp  + lineMeasures_lower$N_hlp))%00%NA,
  #     NVpoints_asym = ((lineMeasures_upper$N_vlp - lineMeasures_lower$N_vlp)/ (lineMeasures_upper$N_vlp  + lineMeasures_lower$N_vlp))%00%NA,
  #     RR_asym      = ((lineMeasures_upper$RR    - lineMeasures_lower$RR)    / (lineMeasures_upper$RR     + lineMeasures_lower$RR))%00%NA,
  #     SING_N_asym  = (lineMeasures_upper$SING_N - lineMeasures_lower$SING_N)/ (lineMeasures_upper$SING_N + lineMeasures_lower$SING_N),
  #     DIV_asym     = (lineMeasures_upper$DIV_dl - lineMeasures_lower$DIV_dl)/ (lineMeasures_upper$DIV_dl + lineMeasures_lower$DIV_dl),
  #     REP_asym     = (lineMeasures_upper$REP_av - lineMeasures_lower$REP_av)/ (lineMeasures_upper$REP_av + lineMeasures_lower$REP_av),
  #     DET_asym     = (lineMeasures_upper$DET    - lineMeasures_lower$DET)   / (lineMeasures_upper$DET    + lineMeasures_lower$DET),
  #     LAM_hl_asym  = (lineMeasures_upper$LAM_hl - lineMeasures_lower$LAM_hl)/ (lineMeasures_upper$LAM_hl + lineMeasures_lower$LAM_hl),
  #     LAM_vl_asym  = (lineMeasures_upper$LAM_vl - lineMeasures_lower$LAM_vl)/ (lineMeasures_upper$LAM_vl + lineMeasures_lower$LAM_vl),
  #     MEAN_dl_asym = (lineMeasures_upper$MEAN_dl- lineMeasures_lower$MEAN_dl)/ (lineMeasures_upper$MEAN_dl+ lineMeasures_lower$MEAN_dl),
  #     MEAN_hl_asym = (lineMeasures_upper$TT_hl  - lineMeasures_lower$TT_hl) / (lineMeasures_upper$TT_hl  + lineMeasures_lower$TT_hl),
  #     MEAN_vl_asym = (lineMeasures_upper$TT_vl  - lineMeasures_lower$TT_vl) / (lineMeasures_upper$TT_vl  + lineMeasures_lower$TT_vl),
  #     MAX_dl_asym  = (lineMeasures_upper$MAX_dl - lineMeasures_lower$MAX_dl)/ (lineMeasures_upper$MAX_dl + lineMeasures_lower$MAX_dl),
  #     MAX_hl_asym  = (lineMeasures_upper$MAX_hl - lineMeasures_lower$MAX_hl)/ (lineMeasures_upper$MAX_hl + lineMeasures_lower$MAX_hl),
  #     MAX_vl_asym  = (lineMeasures_upper$MAX_vl - lineMeasures_lower$MAX_vl)/ (lineMeasures_upper$MAX_vl + lineMeasures_lower$MAX_vl),
  #     ENT_dl_asym  = (lineMeasures_upper$ENT_dl - lineMeasures_lower$ENT_dl)/ (lineMeasures_upper$ENT_dl + lineMeasures_lower$ENT_dl),
  #     ENT_hl_asym  = (lineMeasures_upper$ENT_hl - lineMeasures_lower$ENT_hl)/ (lineMeasures_upper$ENT_hl + lineMeasures_lower$ENT_hl),
  #     ENT_vl_asym  = (lineMeasures_upper$ENT_vl - lineMeasures_lower$ENT_vl)/ (lineMeasures_upper$ENT_vl + lineMeasures_lower$ENT_vl)
  #   )
  #
  #
   }

  #Output
  out_global <- data.frame(emRad      = emRad)

  if(asymmetryUL){
    if(returnUL){
      out_UL <- data.frame(upper.tri = lineMeasures_upper,
                           lower.tri = lineMeasures_lower,
                           ratios.ul = out_ul_asym)
    } else {
      out_UL <- data.frame(ratios.ul = out_ul_asym)
    }
  } else {
     out_UL <- NA
  }

  if(anisotropyHV){
    out_HL <- data.frame(ratios.hv = out_hv_ani)
  } else {
    out_HL <- NA
  }


  out <- as.data.frame(list(out_global, lineMeasures_global, out_UL, out_HL)[c(TRUE, TRUE, asymmetryUL, anisotropyHV)])


  # if(matrices){
  #   tab <- out$crqaMeasures
  # } else {
    tab <- out
  #}

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
                                                                  `N lines`  = tab$ratios.hv.Nlines_ani,
                                                                  `N points` = tab$ratios.hv.N_hlp%00%NA/tab$ratios.hv.N_vlp%00%NA,
                                                                  `Measure` = "LAM",
                                                                  `Rate`    = tab$ratios.hv.LAM_ani,
                                                                  `Mean`    = tab$ratios.hv.MEAN_ani,
                                                                  `Max`     = tab$ratios.hv.MAX_ani,
                                                                  `ENT`     = tab$ratios.hv.ENT_ani)
    if(chromatic){
      rownames(outTable$`Horizontal/Vertical line anisotropy`) <- chromaNames
    }
  }

  if(asymmetryUL){

    outTable$`Upper/Lower triangle asymmetry` <- list(
      `Global Measures` =  data.frame(`Global Ratio` = "U/L of points",
                                      `N points` = tab$ratios.ul.Npoints_asym,
                                      RR = tab$ratios.ul.RR_asym,
                                      Singular = tab$ratios.ul.SING_N_asym,
                                      Divergence = tab$ratios.ul.DIV_asym,
                                      Repetetiveness = tab$ratios.ul.REP_asym),
      `Line-based Measures` = data.frame(
        `Line ratio` = rep(c("D lines", "V lines", "H lines"), each = Nrows),
        `N lines`  = c(tab$ratios.ul.NDlines_asym,
                       tab$ratios.ul.NVlines_asym,
                       tab$ratios.ul.NHlines_asym),
        `N points` = c(tab$ratios.ul.NDpoints_asym,
                       tab$ratios.ul.NVpoints_asym,
                       tab$ratios.ul.NHpoints_asym),
        `Measure` = rep(c("DET","V LAM", "H LAM"), each = Nrows),
        `Rate`    = c(tab$ratios.ul.DET_asym, tab$ratios.ul.LAM_vl_asym, tab$ratios.ul.LAM_hl_asym),
        `Mean`    = c(tab$ratios.ul.MEAN_dl_asym, tab$ratios.ul.MEAN_vl_asym, tab$ratios.ul.MEAN_hl_asym),
        `Max`     = c(tab$ratios.ul.MAX_dl_asym, tab$ratios.ul.MAX_vl_asym, tab$ratios.ul.MAX_hl_asym),
        `ENT`     = c(tab$ratios.ul.ENT_dl_asym, tab$ratios.ul.ENT_vl_asym, tab$ratios.ul.ENT_hl_asym))
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





  return(out)
  # return(list(out = out,
  #             diagonals.dist   = diagonals.dist,
  #             verticals.dist   = verticals.dist,
  #             horizontals.dist = horizontals.dist)
  #)
}




#' rqa_calc
#'
#' @inheritParams rqa_measures
#' @param distributions Return line length distributions.
#'
#' @return CRQA measures and matrices of line distributions (if requested) based on massively parallel analysis, without building the recurrence matrix.
#' @export
#' @keywords internal
#'
rqa_calc <- function(y1,
                     y2 = NA,
                     emRad = NULL,
                     DLmin = 2,
                     VLmin = 2,
                     HLmin = 2,
                     DLmax = NROW(y1),
                     VLmax = NROW(y1),
                     HLmax = NROW(y1),
                     theiler   = NA,
                     AUTO      = NULL,
                     method = "Euclidean",
                     includeDiagonal = NA,
                     chromatic = FALSE,
                     anisotropyHV = FALSE,
                     asymmetryUL = FALSE,
                     returnUL = FALSE,
                     recurrenceTimes = FALSE,
                     distributions = FALSE){

  checkPkg("parallel")
  numCores <- parallel::detectCores()

  if(is.null(AUTO)){
    if(all(is.na(y2))){
      AUTO <- TRUE
    }
  }

  # lineSegments <- rqa_calc_lineMeasures(y1 = y1,
  #                                       y2 = y2,
  #                                       emRad = emRad,
  #                                       DLmin = DLmin, DLmax = DLmax,
  #                                       VLmin = VLmin, VLmax = VLmax,
  #                                       HLmin = HLmin, HLmax = HLmax,
  #                                       d = d, theiler = theiler,
  #                                       recurrenceTimes = recurrenceTimes, AUTO = AUTO,
  #                                       chromatic = chromatic)
  #
  # RP_N <- lineSegments$RP_N
  # diagonals.dist <- lineSegments$diagonals.dist%00%0
  # verticals.dist <- lineSegments$verticals.dist%00%0
  # horizontals.dist <- lineSegments$horizontals.dist%00%0
  # singularities.dist <- lineSegments$singularities.distD%00%0


  # Diagonals ----

  # Includes LOS (main diagonal)
  if(theiler==0){
    rows <- c(-(NROW(y1)-1):(NROW(y1)-1))
  } else {
    rows <- c(-(NROW(y1)-1):(-theiler), (theiler):(NROW(y1)-1))
  }

  if(is.na(asymmetryUL%00%NA)){
    asymmetryUL <- FALSE
  }

  addDiag <- 0
  if(!asymmetryUL){
    if(AUTO){ # If symmetric, we can do just half the matrix
      if(theiler == 0){
        rows <- rows[rows>=1]
        addDiag <- NROW(y1) # Add the diagonal as a line length when theiler = 0
      } else {
        rows <- rows[rows>=theiler]
      }
    }
  } else {
    if(theiler == 0){
      message("NOTE: 'theiler' is set to 0 >> The Upper vs. Lower triangle asymmetry ratio excludes the main diagonal.")
      theiler <- 1
    }
    if(asymmetryUL%in%"upper"){
      rows <- rows[rows<=-theiler]
    }
    if(asymmetryUL%in%"lower"){
      rows <- rows[rows>=theiler]
    }
  }

  outD <- parallel::mcmapply(FUN = rqa_fast, index = rows, MoreArgs = list(y1 = y1, y2 = y2, emRad = emRad, theiler = theiler, symmetrical = AUTO, diagonal = TRUE, method = method), mc.cores = numCores)

  #outD <- lapply(rows, function(r) rqa_fast(index = r, y1 = y1, y2 = y2, emRad = emRad, theiler = theiler, symmetrical = AUTO, diagonal = TRUE, method = method))

  diagonals.distX <- sort(unlist(parallel::mclapply(outD,rqa_lineDist)))

  singularities.dist <- diagonals.distX[diagonals.distX%[]%c(1,1)]

  diagonals.distX <- diagonals.distX[diagonals.distX%[]%c(DLmin,DLmax)]

  if(AUTO){
    RP_N <- (2*sum(sapply(outD,sum))+addDiag)
    diagonals.distY <- diagonals.distX
    singularities.dist <- c(singularities.dist,singularities.dist)
    if(theiler==0){
      diagonals.distX <- c(diagonals.distX, addDiag)
    }
  } else {
    RP_N <- sum(sapply(outD,sum))
    diagonals.distY <- NULL
  }

  diagonals.dist <- sort(c(diagonals.distX,diagonals.distY))

  rm(outD, diagonals.distX, diagonals.distY)

  # Verticals ----
  outV <- parallel::mcmapply(FUN = rqa_fast, index = 1:NROW(y1), MoreArgs = list(y1 = y1, y2 = y2, emRad = emRad, symmetrical = AUTO, diagonal = FALSE, theiler = theiler, method = method), mc.cores = numCores)

  verticals.dist <- sort(unlist(parallel::mclapply(outV,rqa_lineDist)))
  verticals.dist <- verticals.dist[verticals.dist%[]%c(VLmin,VLmax)]

  rm(outV)

  # Horizontals ----

  if(AUTO){
    horizontals.dist <- verticals.dist
  } else {

    outH <- parallel::mcmapply(FUN = rqa_fast, index = 1:NROW(y1), MoreArgs = list(y1 = y1, y2 = y2, emRad = emRad, theiler = theiler, symmetrical = FALSE, diagonal = FALSE, method = method), mc.cores = numCores)

    horizontals.dist <- sort(unlist(parallel::mclapply(outH,rqa_lineDist)))
    horizontals.dist <- horizontals.dist[horizontals.dist%[]%c(HLmin,HLmax)]
    rm(outH)
  }

  # Recurrence Rate ----
  if(all(is.na(y2))&AUTO){
    dims <- c(NROW(y1),NROW(y1))
  } else {
    dims <- c(NROW(y1),NROW(y2))
  }

  if(asymmetryUL){
   RP_max = ((dims[1]-theiler) * (dims[2]-theiler-1))/2
  } else {
   RP_max <- rp_size(dims = dims, AUTO = AUTO, theiler = theiler)$rp_size_theiler
  }

  RR <- RP_N / RP_max

  # Singularities ----
  SING_N <- sum(singularities.dist)
  SING_rate <- SING_N / RP_N

  #Frequency tables of line lengths
  freq_dl <- table(diagonals.dist)
  freq_vl <- table(verticals.dist)
  freq_hl <- table(horizontals.dist)
  freq_hv <- table(c(horizontals.dist,verticals.dist))

  freqvec_dl <- as.numeric(names(freq_dl))
  freqvec_vl <- as.numeric(names(freq_vl))
  freqvec_hl <- as.numeric(names(freq_hl))
  freqvec_hv <- as.numeric(names(freq_hv))

  # Number of lines ----
  N_dl <- sum(freq_dl, na.rm = TRUE)%00%0
  N_vl <- sum(freq_vl, na.rm = TRUE)%00%0
  N_hl <- sum(freq_hl, na.rm = TRUE)%00%0
  N_hv <- sum(freq_hv, na.rm = TRUE)%00%0

  #Number of recurrent points on diagonal, vertical and horizontal lines
  N_dlp <- sum(freqvec_dl*freq_dl, na.rm = TRUE)
  N_vlp <- sum(freqvec_vl*freq_vl, na.rm = TRUE)
  N_hlp <- sum(freqvec_hl*freq_hl, na.rm = TRUE)
  N_hvp <- sum(freqvec_hv*freq_hv, na.rm = TRUE)

  #Determinism / Horizontal and Vertical Laminarity ----
  DET    <- N_dlp/RP_N
  LAM_vl <- N_vlp/RP_N
  LAM_hl <- N_hlp/RP_N
  LAM_hv <- N_hvp/(RP_N*2)

  #Array of probabilities that a certain line length will occur (all >1)
  P_dl <- freq_dl/N_dl
  P_vl <- freq_vl/N_vl
  P_hl <- freq_hl/N_hl
  P_hv <- freq_hv/N_hv

  #Entropy of line length distributions ----
  ENT_dl <- -1 * sum(P_dl * log(P_dl))
  ENT_vl <- -1 * sum(P_vl * log(P_vl))
  ENT_hl <- -1 * sum(P_hl * log(P_hl))
  ENT_hv <- -1 * sum(P_hv * log(P_hv))

  #Relative Entropy (Entropy / Max entropy)
  ENTrel_dl = ENT_dl/(-1 * log(1/DLmax))
  ENTrel_vl = ENT_vl/(-1 * log(1/VLmax))
  ENTrel_hl = ENT_hl/(-1 * log(1/HLmax))
  ENTrel_hv = ENT_hv/(-1 * log(1/max(c(HLmax,VLmax), na.rm = TRUE)))

  #Meanline ----
  MEAN_dl = mean(diagonals.dist, na.rm = TRUE)%00%0
  MEAN_vl = mean(verticals.dist, na.rm = TRUE)%00%0
  MEAN_hl = mean(horizontals.dist, na.rm = TRUE)%00%0
  MEAN_hv = mean(c(horizontals.dist,verticals.dist), na.rm = TRUE)%00%0

  #Maxline ----
  MAX_dl = max(freqvec_dl, na.rm = TRUE)%00%0
  MAX_vl = max(freqvec_vl, na.rm = TRUE)%00%0
  MAX_hl = max(freqvec_hl, na.rm = TRUE)%00%0
  MAX_hv = max(freqvec_hv, na.rm = TRUE)%00%0

  # REPetetiveness ----
  REP_av  <- (N_hlp+N_vlp) / N_dlp
  REP_hl  <-  N_hlp/N_dlp
  REP_vl  <-  N_vlp/N_dlp

  #Coefficient of determination ----
  CoV_dl = stats::sd(diagonals.dist)/mean(diagonals.dist)
  CoV_vl = stats::sd(verticals.dist)/mean(verticals.dist)

  CoV_hl = stats::sd(horizontals.dist)/mean(horizontals.dist)
  CoV_hv = stats::sd(c(horizontals.dist,verticals.dist))/mean(c(horizontals.dist,verticals.dist))

  #Divergence ----
  DIV_dl = 1/MAX_dl
  DIV_vl = 1/MAX_vl
  DIV_hl = 1/MAX_hl
  DIV_hv = 1/MAX_hv

  out <- data.frame(
    RP_max    = RP_max,
    RP_N      = RP_N,
    RR        = RR,
    SING_N    = SING_N,
    SING_rate = SING_rate,
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

  if(distributions){
    return(list(
      crqaMeasures = out,
      crqaMatrices = list(dlines = diagonals.dist,
                          vlines = verticals.dist,
                          hlines = horizontals.dist,
                          freq_dl = freq_dl,
                          freq_vl = freq_vl,
                          freq_hl = freq_hl))
    )
  } else {
    return(out)
  }

}


#' rqa_lineMeasures
#'
#' @inheritParams rqa
#'
#' @return
#' @export
#'
#' @keywords internal
#'
rqa_calc_lineMeasures <- function(y1,
                                  y2 = NA,
                                  emRad = NA,
                                  DLmin = 2,
                                  VLmin = 2,
                                  HLmin = 2,
                                  DLmax = NROW(y1),
                                  VLmax = NROW(y1),
                                  HLmax = NROW(y1),
                                  d         = NULL,
                                  theiler   = NA,
                                  recurrenceTimes    = FALSE,
                                  AUTO      = NULL,
                                  chromatic = FALSE,
                                  matrices  = FALSE){


  # Diagonal ----

  # Includes LOS (main diagonal)
  if(theiler==0){
    rows <- c(-(NROW(y1)-1):(NROW(y1)-1))
  } else {
    rows <- c(-(NROW(y1)-1):(-theiler), (theiler):(NROW(y1)-1))
  }

  if(AUTO){
    addDiag <- 0
    if(theiler == 0){
      rows <- rows[rows>=1]
      addDiag <- NROW(y1)
    } else {
      rows <- rows[rows>=theiler]
    }
  }

  outD <- parallel::mcmapply(FUN = rqa_fast, index = rows, MoreArgs = list(y1 = y1, y2 = y2, emRad = emRad, theiler = theiler, symmetrical = AUTO, diagonal = TRUE, method = method), mc.cores = numCores)

  diagonals.distX <- sort(unlist(parallel::mclapply(outD,rqa_lineDist)))

  singular.distD <- diagonals.distX[diagonals.distX%[]%c(1,1)]

  diagonals.distX <- diagonals.distX[diagonals.distX%[]%c(DLmin,DLmax)]

  if(AUTO){
    RP_N <- (2*sum(sapply(outD,sum))+addDiag)
    diagonals.distY <- diagonals.distX
    singular.distD <- c(singular.distD,singular.distD)
    if(theiler==0){
      diagonals.distX <- c(diagonals.distX, addDiag)
    }
  } else {
    RP_N <- sum(sapply(outD,sum))
    diagonals.distY <- NULL
  }

  diagonals.dist <- sort(c(diagonals.distX,diagonals.distY))

  rm(outD, diagonals.distX, diagonals.distY)


  # Vertical ----
  outV <- parallel::mcmapply(FUN = rqa_fast, index = 1:NROW(y1), MoreArgs = list(y1 = y1, y2 = y2, emRad = emRad, symmetrical = AUTO, diagonal = FALSE, theiler = theiler, method = method), mc.cores = numCores)

  verticals.dist <- sort(unlist(parallel::mclapply(outV,rqa_lineDist)))
  verticals.dist <- verticals.dist[verticals.dist%[]%c(VLmin,VLmax)]

  rm(outV)

  # Horizontal ----

  if(AUTO){
    horizontals.dist <- verticals.dist
  } else {

    outH <- parallel::mcmapply(FUN = rqa_fast, index = 1:NROW(y1), MoreArgs = list(y1 = y1, y2 = y2, emRad = emRad, theiler = theiler, symmetrical = FALSE, diagonal = FALSE, method = method), mc.cores = numCores)

    horizontals.dist <- sort(unlist(parallel::mclapply(outH,rqa_lineDist)))
    horizontals.dist <- horizontals.dist[horizontals.dist%[]%c(HLmin,HLmax)]
    rm(outH)
  }

  return(list(RP_N = RP_N,
              diagonals.dist = diagonals.dist,
              verticals.dist = verticals.dist,
              horizontas.dist = horizontals.dits,
              singular.distD = singular.distD))

}

#' Stitch rows
#'
#' @param rowA rowA
#' @param rowB rowB
#' @param revA revA
#'
#' @return
#'
#' @export
#'
rqa_stitchRows <- function(rowA, rowB, revA = FALSE){
  if(revA){
    nA <- length(rowA)
    nB <- length(rowB)
    return(bit::as.bit(c(rev(rowA[-nA]),rowB[1:nB])))
  } else {
    return(bit::as.bit(c(rowA,rowB)))
  }
}

