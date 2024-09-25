# package casnet ----
#
# Estimate Parameters ----
#
# est functions

#' Estimate Radius.
#'
#' Find a fixed or optimal radius.
#'
#' @aliases crqa_radius
#' @inheritParams rp
#' @param RM Unthresholded Recurrence Matrix
#' @param type Either `"fixed"` (default) or `"optimal"`, `"fixed"` will search for a radius that is close to the value for the `targetMeasure` in `targetValue`, `"optimal"` will optimise the radius for the `targetMeasure`, `targetValue` is ignored.
#' @param startRadius If `type = "fixed"` this is the starting value for the radius (default = percentile of unique distances in RM given by `targetValue`). If `type = "optimal"` this will be a range of radius values (in normalised SD units) that will be considered (default = `seq(0,2,by=.01)`)
#' @param eachRadius If `type = "optimal"` this is the number of signal and noise series that will be generated for each level in `startRadius` (default = `1`)
#' @param targetMeasure If `type = "optimal"`, it must be a character vector indicating which recurrence measure to optimise the radius for, options are "RR" (default), "DET", "LAM", "T1", and "all". The option `targetMeasure = "all"` will report all the optimal values obtained from one realisation of `startRadius * eachRadius` signal and noise series.
#' @param targetValue When argument `type` is set to "fixed", the value represents the target value for the measure in `targetMeasure` (default = `RR = .05`).
#' @param tol Tolerance for achieving `targetValue` for `targetMeasure` (default = `0.01`)
#' @param maxIter If `type = "fixed"`: Maximum number of iterations to reach targetValue.
#' @param theiler Size of theiler window (default `0`)
#' @param histIter Return iteration history? (default = `FALSE`)
#' @param noiseLevel Noise level to construct the `signal + noiseLevel *` \eqn{N(\mu=0,\sigma=1)} (default = `0.75`)
#' @param noiseType Type
#' @param plotROC Generates an ROC plot if `type = "optimal"`
#' @param standardise Standardise `y` if `type == "optimal"`
#' @param radiusOnFail Radius to return when search fails `"tiny" = 0 + ,Machine.double.eps`, this will likely cause a matrix full of zeros. `"huge" = 1 + max. distance in RM`, which will give a matrix full of ones, `"percentile" = quantile(RM, prob = targetValue) of distances greater than 0`.
#' @param silent Silent-ish
#'
#' @family Estimate Recurrence Parameters
#'
#' @return A dataframe listing settings used to search for the radius, the radius found given the settings and the recurrence rate produced by the radius (either 1 row or the entire iteration history)
#' @export
#'
est_radius <- function(RM = NULL,
                       y1 = NULL,
                       y2 = NULL,
                       emLag = 1,
                       emDim = 1,
                       method = "Euclidean",
                       type           = c("fixed","optimal")[1],
                       startRadius    = NULL,
                       eachRadius     = 1,
                       targetMeasure  = c("RR","DET","LAM","T1","all")[1],
                       targetValue    = 0.05,
                       tol            = 0.01,
                       maxIter        = 100,
                       theiler        = NA,
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
    RM <- rp(y1=y1, y2=y2,emDim=emDim, emLag=emLag, method = method)
  }

  if(any(is.na(RM))){
    NAid <- which(is.na(Matrix::as.matrix(RM)))
    RM[NAid] <- max(Matrix::as.matrix(RM), na.rm = TRUE)+1
  }

  if(all(as.vector(RM)%in%c(0,1))){
    stop("Matrix RM is already a binary matrix!")
  }

  # check auto-recurrence
  RM   <- rp_checkfix(RM, checkAUTO = TRUE)
  AUTO <- attr(RM,"AUTO")

  RM <- setTheiler(RM = RM, theiler = theiler, silent = TRUE)

  if(is.null(startRadius)){
    if(type=="fixed"){
      if(AUTO){
        startRadius <- as.numeric(stats::quantile(unique(as.vector(RM[lower.tri(RM)])),probs = ifelse(targetValue>=1,.05,targetValue)))
      } else {
        startRadius <- as.numeric(stats::quantile(unique(as.vector(RM)),probs = ifelse(targetValue>=1,.05,targetValue)))
        }
    } else {
      startRadius <- seq(0,1.5,by=0.001)
    }
  }


  if(type%in%"fixed"){

    if(tol%][%c(0,1)){stop("Argument tol must be between 0 and 1.")}

    tryRadius <- startRadius
    Measure   <- 0
    iter      <- 0
    Converged <- FALSE
    minRRfound <- FALSE
    seqIter    <- 1:maxIter
    stopRadius <- NA
    RP_N  <- NA
    tollo <- targetValue-tol #(1-tol),
    tolhi <- targetValue+tol #(tol+1),
    rp.size <- rp_size(RM=RM,AUTO = AUTO,theiler = theiler)$rp_size_theiler
    minDist <- suppressMessages(min(RM[RM>0], na.rm = TRUE))
    minRR   <- (sum((as.vector(RM)>0)&(as.vector(RM)==minDist)))/rp.size



    iterList <- data.frame(iter        = seqIter,
                           Measure     = Measure,
                           Radius      = tryRadius,
                           targetValue = targetValue,
                           tollo       = tollo, #(1-tol),
                           tolhi       = tolhi, #(tol+1),
                           startRadius = startRadius,
                           stopRadius  = stopRadius,
                           rp.size     = rp.size,
                           rp.points   = RP_N,
                           AUTO        = AUTO,
                           Converged   = Converged, check.names = FALSE)

    exitIter <- FALSE
    if(!silent){cat(paste("\nSearching for a radius that will yield",targetValue,"??",tol,"for", targetMeasure,"\n"))}

    if(tryRadius<=minDist){
      warning(paste("The minimum RR possible for this matrix is", round(Measure,3),
                    "because the minimum distance is:", round(minDist,3)))
      minRRfound <- TRUE
      exitIter   <- TRUE
    }

    # if(theiler == 0){
    #   if((NROW(y1)/RP_max)>targetValue){
    #     stop(paste0("The minimum RR possible including the diagonal is: ",round(NROW(y1)/RP_max,3),", which is larger than the targetValue for RR: ",targetValue))
    #   }
    # }

    #RM <- float::as.float(Matrix::as.matrix(RM))


    # p <- dplyr::progress_estimated(maxIter)
    while(!exitIter){

      iter <- iter+1
      #p$tick()$print()

      if(!silent){cat(paste("Iteration",iter,"\n"))}

      RP_N <- sum((as.vector(RM)>0)&(as.vector(RM)<tryRadius))

      # RMs <- mat_di2bi(RM, emRad = tryRadius, convMat = TRUE)
      # #Total nr. recurrent points
      # RP_N <- Matrix::nnzero(RMs, na.counted = FALSE)
      # rm(RMs)

      #RP_N <- RP_N-minDiag
      #rp.size <- rp_size(RMs) #length(RMs)

      Measure <- RP_N/rp.size

      # crpOut <- rp_measures(RM = RMs, emRad = tryRadius, AUTO=AUTO)
      #Measure  <-  crpOut[[targetMeasure]]
      #crpOut <- data.frame(RR = RR, RT = RT, size = length(RMs))
      #Measure <- RR

      iterList[iter,] <-    cbind.data.frame(iter        = iter,
                                             Measure     = Measure,
                                             Radius      = tryRadius,
                                             targetValue = targetValue,
                                             tollo       = tollo, #(1-tol),
                                             tolhi       = tolhi, #(tol+1),
                                             startRadius = startRadius,
                                             stopRadius  = stopRadius,
                                             rp.size     = rp.size,
                                             rp.points   = RP_N,
                                             AUTO        = AUTO,
                                             Converged   = Converged)

      if(tryRadius<=minDist){
        warning(paste("The minimum RR possible for this matrix is",round(minRR,3), "because the minimum distance is:",round(minDist,3)))
        minRRfound <- TRUE
        exitIter <- TRUE
      }

      if(any(Measure%[]%c(tollo,tolhi),(iter>=maxIter))){
        if(Measure%[]%c(tollo,tolhi)){
          Converged <- TRUE
          if(!silent){
            message("\nConverged! Found an appropriate radius...")
            }
        }
        iterList$Converged[iter] <- Converged
        exitIter <- TRUE
      }

      if(round(Measure,digits = 2)>round(targetValue,digits = 2)){
        tryRadius <- tryRadius*(min(c(0.8,1-abs(round(Measure,digits = 2)-round(targetValue,digits = 2))))) # tol*2
      } else {
        tryRadius <- tryRadius*(min(c(1.8,1+abs(round(Measure,digits = 2)-round(targetValue,digits = 2))))) #1+(tol*2)
      }


    } # While ....

    if(iter>=maxIter){
      warning("Max. iterations reached!")
      iterList$stopRadius[iter] <- tryRadius
      #iterlist$Measure[iter] <- Measure
    }
    if(!minRRfound){
    if(Measure %][% c(tollo,tolhi)){
      iterList$Radius[iter] <- dplyr::case_when(
        radiusOnFail%in%"tiny" ~ 0 + .Machine$double.eps,
        radiusOnFail%in%"huge" ~ 1 + max(RM),
        radiusOnFail%in%"percentile" ~ as.numeric(stats::quantile(unique(as.vector(Matrix::tril(RM,-1))), probs = ifelse(targetValue>=1,.05,targetValue)))
      )
      warning(paste0("\nTarget not found, try increasing tolerance, max. iterations, or, change value of startRadius.\nreturning radius: ",iterList$Radius[iter]))
      iterList$stopRadius[iter] <- tryRadius
      #iterlist$Measure[iter] <- Measure
    }
    } else {
      iterList$Radius[iter] <- minDist
      }

    ifelse(histIter,id<-c(1:iter),id<-iter)
    return(iterList[id,])

  } # if "fixed"

  if(type%in%"optimal"){

    if(optimOK){

      if(!silent){cat(paste0("\nNormalisation set to: ",standardise,"!!\n"))}

      startRadius <- rep(startRadius, each = eachRadius)

      dfREC  <-  plyr::ldply(startRadius, function(r){est_parameters_roc(y = y1,
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

      optimal.value <- list(RR  = rp_cl(y1, emRad = optimal.radius$RR, emDim=emDim, emLag=emLag)[["RR"]],
                            DET = rp_cl(y1, emRad = optimal.radius$DET, emDim=emDim, emLag=emLag)[["DET"]],
                            LAM = rp_cl(y1, emRad = optimal.radius$LAM, emDim=emDim, emLag=emLag)[["LAM"]],
                            T1  = rp_cl(y1, emRad = optimal.radius$T1, emDim=emDim, emLag=emLag)[["T1"]]
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

#' Estimate Radius without building a recurrence matrix
#'
#' @description Find a fixed radius without building the recurrence matrix.
#'
#' @inheritParams rqa_par
#' @param startRadius The starting value for the radius (default = SD of time series values)
#' @param targetValue When argument `type` is set to "fixed", the value represents the target value for the measure in `targetMeasure` (default = `RR = .05`).
#' @param tol Tolerance for achieving `targetValue` for `targetMeasure` (default = `0.01`)
#' @param maxIter If `type = "fixed"`: Maximum number of iterations to reach targetValue.
#' @param theiler Size of theiler window (default `0`)
#' @param histIter Return iteration history? (default = `FALSE`)
#' @param noiseLevel Noise level to construct the `signal + noiseLevel *` \eqn{N(\mu=0,\sigma=1)} (default = `0.75`)
#' @param noiseType Type
#' @param plotROC Generates an ROC plot if `type = "optimal"`
#' @param standardise Standardise `y` if `type == "optimal"`
#' @param radiusOnFail Radius to return when search fails `"tiny" = 0 + ,Machine.double.eps`, this will likely cause a matrix full of zeros. `"huge" = 1 + max. distance`, which will give a matrix full of ones, `"minimum" = minimum distance in matrix`.
#' @param silent Silent-ish
#' @param useParallel Should evaluation run using package parallel? This is will only be beneficial if the time series contains more than 10k data points (default = `TRUE`)
#'
#' @family Estimate Recurrence Parameters
#'
#' @return A dataframe listing settings used to search for the radius, the radius found given the settings and the recurrence rate produced by the radius (either 1 row or the entire iteration history)
#' @export
#'
est_radius_rqa <- function(y1 = NULL,
                       y2 = NULL,
                       AUTO = NULL,
                       method = "Euclidean",
                       startRadius    = NULL,
                       targetValue    = 0.05,
                       tol            = 0.01,
                       maxIter        = 100,
                       theiler        = NA,
                       histIter       = FALSE,
                       standardise  = c("mean.sd","median.mad","none")[3],
                       radiusOnFail   = c("tiny","huge","percentile")[3],
                       silent         = FALSE,
                       useParallel    = TRUE,
                       doEmbed        = TRUE){

  if(useParallel){
    checkPkg("future.apply")
    future::plan("multisession")
    #numCores <- parallel::detectCores()
  }

  if(is.null(AUTO)){
    if(is.null(y2)|all(is.na(y2))){
    message("One time series, so assuming Auto-RQA...")
    AUTO <- TRUE
    } else {
      stop("Auto-RQA or Cross-RQA? (Provide a value for AUTO)")
    }
  }


  if(!doEmbed){
    emDim <- 1
    emLag <- 0
  }

  # Embed series
  tmp <- rqa_getSeries(y1 = y1,
                       y2 = y2,
                       emDim = emDim,
                       emLag = emLag,
                       chromatic = FALSE,
                       doEmbed = doEmbed)

  y1 <- tmp$et1
  y2 <- tmp$et2

  rm(tmp)

  if(is.na(theiler)){
    if(AUTO){
      theiler <- 1
    } else {
      theiler <- 0
    }
  }

  if(AUTO){
    y2 <- NA
  }

  if(is.null(attributes(y1)$embedding.lag)){
    stop("Please use `ts_embed()` to generate `y1` and/or `y2`.")
  }

  if(all(is.na(y2%00%NA))){
    dims <- c(NROW(y1), NROW(y1))
  } else {
    dims <- c(NROW(y1), NROW(y2))
  }

  RP_max <- rp_size(dims = dims, AUTO = AUTO, theiler = theiler)$rp_size_theiler

  if(theiler == 0){
    if((NROW(y1)/RP_max)>targetValue){
      stop(paste0("The minimum RR possible including the diagonal is: ",round(NROW(y1)/RP_max,3),", which is larger than the targetValue for RR: ",targetValue))
    }
  }

  dist_method <- return_error(proxy::pr_DB$get_entry(method))
  if("error"%in%class(dist_method$value)){
    stop("Unknown distance metric!\nUse proxy::pr_DB$get_entries() to see a list of valid options.")
  }
  if(method == "SBD"){
    message("Make sure all phase space dimensions are on the same scale when using Shape Based Distance.")
  }

  if(is.null(startRadius)){
    #startRadius <- mean(c(as.numeric(minDist),as.numeric(maxDist)), na.rm = TRUE)
    if(!AUTO){
      startRadius <- as.numeric(stats::quantile(unique(proxy::dist(x = y1, y = y2, pairwise = FALSE, method = method)),
                                                probs = ifelse(targetValue>=1, .05, targetValue)))
    } else {
      startRadius <- as.numeric(stats::quantile(unique(proxy::dist(x = y1, y = y1, pairwise = FALSE, method = method)),
                                                probs = ifelse(targetValue>=1, .05, targetValue)))
    }
  }

    # as.numeric(quantile(proxy::dist(y1, pairwise = FALSE), probs = targetValue))
    tryRadius <- startRadius
    Measure   <- 0
    iter      <- 0
    Converged <- FALSE
    minRRfound <- FALSE
    seqIter    <- 1:maxIter
    stopRadius <- NA
    RP_N  <- NA

    if(tol%][%c(0,1)){stop("Argument tol must be between 0 and 1.")}
    tollo <- targetValue-tol #(1-tol),
    tolhi <- targetValue+tol #(tol+1),

    iterList <- data.frame(iter        = seqIter,
                           Measure     = Measure,
                           Radius      = tryRadius,
                           targetValue = targetValue,
                           tollo       = tollo, #(1-tol),
                           tolhi       = tolhi, #(tol+1),
                           startRadius = startRadius,
                           stopRadius  = stopRadius,
                           rp.size     = RP_max,
                           rp.points   = RP_N,
                           AUTO        = AUTO,
                           Converged   = Converged, check.names = FALSE)

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

    # while ----

    exitIter <- 0

    while(!exitIter){

      if(useParallel){

       # outD <- parallel::mcmapply(FUN = rqa_fast, index = rows, MoreArgs = list(y1 = y1, y2 = y2, emRad = tryRadius, theiler = theiler, symmetrical = AUTO, diagonal = TRUE, method = method), mc.cores = numCores)
        #parallel::stopCluster(numCores)

        outD <- future.apply::future_mapply(FUN = rqa_fast, index = rows, MoreArgs = list(y1 = y1, y2 = y2, emRad = tryRadius, theiler = theiler, symmetrical = AUTO, diagonal = TRUE, method = method))

      } else {

        outD <- lapply(rows, function(i) rqa_fast(y1 = y1, y2 = y2, index = i, emRad = tryRadius, theiler = theiler, symmetrical = AUTO, diagonal = TRUE, method = method))



      }

      if(iter == 0){

        if(useParallel){

          minDist    <- float::as.float(future.apply::future_sapply(outD, function(d) as.numeric(attributes(d)$minDist)))
          minDistN   <- float::as.float(future.apply::future_sapply(outD, function(d) as.numeric(attributes(d)$minDistN)))
          minRP_dist <- float::as.float(min(minDist, na.rm = TRUE))
          minRP_N    <- sum(minDistN[minDist==minRP_dist], na.rm = TRUE)
          maxDist    <- float::as.float(future.apply::future_sapply(outD, function(d) as.numeric(attributes(d)$maxDist)))
          maxRP_dist <- float::as.float(max(maxDist, na.rm = TRUE))

        } else {

          minDist    <- float::as.float(sapply(outD, function(d) as.numeric(attributes(d)$minDist)))
          minDistN   <- float::as.float(sapply(outD, function(d) as.numeric(attributes(d)$minDistN)))
          minRP_dist <- float::as.float(min(minDist, na.rm = TRUE))
          minRP_N    <- sum(minDistN[minDist==minRP_dist], na.rm = TRUE)
          maxDist    <- float::as.float(sapply(outD, function(d) as.numeric(attributes(d)$maxDist)))
          maxRP_dist <- float::as.float(max(maxDist, na.rm = TRUE))

        }

        if(AUTO){
          minRP_N <- (minRP_N*2) + addDiag
        }

        minRR <- minRP_N/RP_max

        # if(any(is.na(RM))){
        #   NAid <- which(is.na(Matrix::as.matrix(RM)))
        #   RM[NAid] <- max(Matrix::as.matrix(RM), na.rm = TRUE)+1
        # }
        #
        # if(!rescaleDist%in%"none"){
        #   warning("Cannot rescale distance matrix, please rescale the time series to standardeviation (z-score) or unit scale (min/max)")
        # }

        if(targetValue<=minRR){
          warning(paste("The minimum RR possible for this matrix is", round(minRR,3),
                        "because the minimum distance is:", round(as.numeric(minRP_dist),3)))
          minRRfound <- TRUE
          exitIter   <- TRUE
        }

      }

      iter <- iter+1
      if(!silent){cat(paste("Iteration",iter,"\n"))}

      if(useParallel){

      if(AUTO){
        RP_N <- (2*sum( future.apply::future_sapply(outD,sum)))+addDiag
      } else {
        RP_N <- sum( future.apply::future_sapply(outD,sum))
      }

      } else {

        if(AUTO){
          RP_N <- (2*sum(sapply(outD,sum)))+addDiag
        } else {
          RP_N <- sum(sapply(outD,sum))
        }
      }


      Measure <- RP_N/RP_max

      iterList[iter,] <-    cbind.data.frame(iter        = iter,
                                             Measure     = Measure,
                                             Radius      = tryRadius,
                                             targetValue = targetValue,
                                             tollo       = tollo, #(1-tol),
                                             tolhi       = tolhi, #(tol+1),
                                             startRadius = startRadius,
                                             stopRadius  = stopRadius,
                                             rp.size     = RP_max,
                                             rp.points   = RP_N,
                                             AUTO        = AUTO,
                                             Converged   = Converged)

      if(any(Measure%[]%c(tollo,tolhi),(iter>=maxIter),(Measure<=minRR))){
        if(Measure%[]%c(tollo,tolhi)){
          Converged <- TRUE
          if(!silent){
            message("\nConverged! Found an appropriate radius...")
          }
        }
        iterList$Converged[iter] <- Converged
        exitIter <- TRUE
        if(Measure<=minRR){
          minRRfound <- TRUE
        }
      }

      if(round(Measure,digits = 2)>round(targetValue,digits = 2)){
        tryRadius <- tryRadius*(min(c(0.8,1-abs(round(Measure,digits = 2)-round(targetValue,digits = 2))))) # tol*2
      } else {
        tryRadius <- tryRadius*(min(c(1.8,1+abs(round(Measure,digits = 2)-round(targetValue,digits = 2))))) #1+(tol*2)
      }


    } # While ....

    if(useParallel){

      future::plan("sequential")

    }

    if(iter>=maxIter){
      warning("Max. iterations reached!")
      iterList$stopRadius[iter] <- tryRadius
      #iterlist$Measure[iter] <- Measure
    }

    if(any(minRRfound,!Converged)){
      if(Measure %)(% c(tollo,tolhi)){
        iterList$Radius[iter] <- dplyr::case_when(
          radiusOnFail%in%"tiny" ~ 0 + .Machine$double.eps,
          radiusOnFail%in%"huge" ~ 1 + round(as.numeric(maxRP_dist),3),
          radiusOnFail%in%"percentile" ~ round(as.numeric(startRadius),3)
        )
        #iterlist$Measure[iter] <- Measure
    } else {
      iterList$Radius[iter] <- minDist
    }
      warning(paste0("\nTarget not found, try increasing tolerance, max. iterations, or, change value of startRadius.\nreturning radius: ",iterList$Radius[iter]))
      iterList$stopRadius[iter] <- tryRadius
    }

    ifelse(histIter,id<-c(1:iter),id<-iter)
    return(iterList[id,])
}




#' Estimate RQA parameters
#'
#' Find optimal parameters for constructing a Recurrence Matrix. A wrapper for various algorithms used to find optimal values for the embedding delay and the number of embedding dimensions.
#'
#' @param y A numeric vector or time series
#' @param emLag Optimal embedding lag (delay), e.g., provided by an optimising algorithm. If `NULL` the lags based on the mutual information in `lagMethods` will be reported. If a numeric value representing a valid lag is passed, this value will be used to estimate the number of dimensions (default = `NULL`)
#' @param lagMethods A character vector with one or more of the following strings: `"first.minimum","global.minimum","max.lag"`. If `emLag` represents a valid lag this value will be reported as `"user.lag"` (default = `c("first.minimum","global.minimum","max.lag")`)
#' @param maxLag Maximum embedding lag to consider. If `NA` then the value is caclulated as `floor(length(y)/(maxDim+1))` (default = `NA`)
#' @param estimateDimensions Decide on an optimal embedding dimension relative to the values in `maxDim` and `lagMethods`, according to a number of preferences passed as a character vector. The order in which the preferences appear in the vector affects the selection procedure, with index `1` being most important preference. The following options are available:
#'
#' * `preferNone` - No optimal number will be picked all other preferences will be ignored
#' * `preferSmallestDim` - Pick smallest number of dimensions associated with a percentage NN below `nnThres`
#' * `preferSmallestNN` - Pick the number of dimensions that is associated with the smallest percentage NN below `nnThres`
#' * `preferSmallestLag` - If the value of `nnThres` does not lead to a unique preference for a pair of dimension and lag values, use the pair with the smallest lag
#' * `preferSmallestInLargestHood` - The default option: If no unique pair can be found, prefer pairs with smallest values for lag, dimensions, percentage NN for the largest NN size
#'
#' @param maxDim Maximum number of embedding dimensions to consider (default = `10`)
#' @param minVecLength The minimum length of state space vectors after delay-embedding. For short time series, this will affect the possible values of `maxDim` that can be used to evaluate the drop in nearest neighbours. In general it is not recommended to evaluate high dimensional state spaces, based on a small number of state soace coordinates, the default is an absolute minimum and possibly even lower than that. (default = `20`)
#' @param nnSize Neighbourhood diameter (integer, the `number.boxes` parameter of [tseriesChaos::false.nearest()]) used to speed up neighbour search. (default = `NA`)
#' @param nnRadius Points smaller than the radius are considered neighbours. If `NA` the value will be `sd(y)/10`  (default = `NA`)
#' @param nnThres Threshold value (in percentage 0-100) representing the percentage of Nearest Neighbours that would be acceptable when using N surrogate dimensions. The smallest number of surrogate dimensions that yield a value below the threshold will be considered optimal (default = `10`)
#' @param theiler Theiler window on distance matrix (default = `0`)
#' @param doPlot Produce a diagnostic plot the results (default = `TRUE`)
#' @param silent Silent-ish mode
#' @param ... Other parameters passed to [nonlinearTseries::timeLag()]
#'
#' @return A list object containing the optimal values (as indicated by the user) and iteration history.
#'
#' @details A number of functions are called to determine optimal parameters for delay embedding a time series:
#'
#' * Embedding lag (`emLag`): The default is to call [casnet::est_emLag()], which is a wrapper around [nonlinearTseries::timeLag()] with `technique=ami` to get lags based on the mutual information function.
#' * Embedding dimension (`m`, `emDim`): The default is to call [casnet::est_emDim()], which is a wrapper around [tseriesChaos::false.nearest()]
#'
#'
#' @family Estimate Recurrence Parameters
#'
#' @export
#'
#' @examples
#'
#' set.seed(4321)
#' est_parameters(y=rnorm(100))
#'
est_parameters <- function(y,
                           lagMethods = c("first.minimum","global.minimum","max.lag"),
                           estimateDimensions = "preferSmallestInLargestHood",
                           maxDim   = 10,
                           emLag    = NULL,
                           maxLag   = NA,
                           minVecLength = 20,
                           nnSize  = NA,
                           nnRadius = NA,
                           nnThres  = 10,
                           theiler  = 0,
                           doPlot   = TRUE,
                           silent   = TRUE,
                           ...){

  if(!is.null(dim(y))){stop("y must be a 1D numeric vector!")}

  if(length(nnRadius)!=1){stop("nnRadius must have 1 numeric value")}
  if(length(nnSize)!=1){stop("nnSize must have 1 numeric value")}

  if(any(is.na(y))){warning("Removing NA in y before estimation!")}
  y <- y[!is.na(y)]

  if(is.na(maxLag)){
    maxLag <- floor(NROW(y)/(maxDim+1))
  }
  if(is.na(nnRadius)){
    nnRadius <- sd(y, na.rm = TRUE)/10
  }

  #y <- ts_standardise(y, adjustN = FALSE)

  if(minVecLength<20){
    stop("Please collect more data!")
  }

  if(!is.na(maxDim%00%NA)&!is.na(maxLag%00%NA)){
    if((NROW(y)-(maxDim*maxLag))<minVecLength){
      maxDim <- (1:maxDim)[max(which(NROW(y)-(1:maxDim*maxLag)>=minVecLength), na.rm = TRUE)]
      message(paste("Changed value of maxDim to",maxDim,"\n"))
    }
  }

  emDims <-  1:maxDim

  doLags <- c(1:maxLag)
  if(!is.null(emLag)){
    if(NROW(emLag)==1){
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
      if(emLags$selection.methods[m]%in%"first.minimum"){
        if(length(doLags)>2){
          emLags$lag[m] <- which(attributes(ts_symbolic(data.frame(mi)))$mi_sym_label%in%"trough")[1]%00%NA
          if(is.na(emLags$lag[m])){
            emLags$lag[m] <- as.numeric(which.min(mi))
          }
        } else {
          warning("Only 2 lags to evaluate, setting 'first.minimum' to 1.")
          emLags$lag[m] <- 1
        }
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
      if(is.na(emLags$lag[m])){
        emLags$ami[m] <- NA
      } else {
        emLags$ami[m] <- mi[emLags$lag[m]]
      }
    }
  } else {
    emLags <- cbind.data.frame(selection.methods = "Not estimated", lag = NA, ami = NA)
  }

  # (fn.out <- tseriesChaos::false.nearest(lx, m=10, d=17, t=0, eps=sd(lx)/10, rt=20))
  # plot(fn.out[1,],type="b")

  if(any(estimateDimensions%in%c("preferNone","preferSmallestDim", "preferSmallestNN", "preferSmallestLag", "preferSmallestInLargestHood"))){

    lagList <- list()
    cnt = 0

    if(is.na(nnSize)){
      number.boxes <- NULL
      nnSize <- "Auto"
    } else {
      number.boxes <- nnSize
    }

    #for(N in seq_along(nnSize)){
      for(R in seq_along(nnRadius)){
        for(L in seq_along(emLags$selection.method)){
          cnt = cnt+1
          if(!is.na(emLags$lag[L])){
            #cnt = cnt+1

            Nn.max <- Nn.mean <- Nn.sd <- Nn.min <- numeric(maxDim)
            for(D in seq_along(emDims)){
              RM <- rp(y,y, emDim = emDims[D], emLag = emLags$lag[L],returnMeasures = FALSE)
              RM <- bandReplace(RM,-theiler,theiler,0,silent = silent)
              Nn.min[D]  <- min(RM, na.rm = TRUE)
              Nn.max[D]  <- max(RM, na.rm = TRUE)
              Nn.sd[D]   <- stats::sd(RM, na.rm = TRUE)
              Nn.mean[D] <- mean(Matrix::as.matrix(RM), na.rm = TRUE)
              rm(RM)
            }

            #surrDims <- nonlinearTseries::buildTakens(time.series =  as.numeric(y), embedding.dim =  emDims[[D]], time.lag = emLags$optimal.lag[L])
            #fnnSeries <- false.nearest(as.numeric(y), m=maxDim, d=emLags$lag[L], t=theiler, rt= nnRadius[R], eps = nnSize)
            # plot(fn.out)

            fnnSeries <- fnn(y = y,
                             emLag = emLags$lag[L],
                             maxDim = maxDim,
                             radius = nnRadius[R],
                             number.boxes = number.boxes)

            lagList[[cnt]] <- data.frame(Nn.pct = fnnSeries[,1]/fnnSeries[1,1]*100,
                                         Nsize = nnSize,
                                         Nradius = nnRadius[R],
                                         emLag.method = emLags$selection.method[[L]],
                                         emLag = emLags$lag[L],
                                         emDim = emDims,
                                         Nn.mean = Nn.mean,
                                         Nn.sd  = Nn.sd,
                                         Nn.min = Nn.min,
                                         Nn.max = Nn.max)



          } else {

            lagList[[cnt]] <- data.frame(Nn.pct = NA,
                                         Nsize = nnSize,
                                         Nradius = nnRadius[R],
                                         emLag.method = emLags$selection.method[[L]],
                                         emLag = emLags$lag[L],
                                         emDim = emDims,
                                         Nn.mean = NA,
                                         Nn.sd  = NA,
                                         Nn.min = NA,
                                         Nn.max = NA)

          }

          #allN <- nonlinearTseries::findAllNeighbours(surrDims, radius = nnSize*sd(y))
          #Nn <- sum(plyr::laply(allN, length), na.rm = TRUE)
          #   if(D==1){Nn.max <- Nn}


        } #R
      } #L
   # } #N

    df <- plyr::ldply(lagList)
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
                  preferSmallestDim = opt[which.min(opt$emDim),],
                  preferSmallestNN = opt[which.min(opt$NN.pct),],
                  preferSmallestLag = opt[which.min(opt$emLag),],
                  preferSmallestInLargestHood = opt[which(min(opt$emLag, na.rm=TRUE)&min(opt$emDim, na.rm=TRUE)&max(opt$Nradius, na.rm = TRUE)),]
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
  # df$Nns   <- interaction(df$Nsize,df$Nradius)
  df$Nns.f  <-  factor(paste("NN radius =",round(df$Nradius, digits = 4)))

  if(doPlot){

    dfs <- data.frame(startAt= c(.5, graphics::hist(emDims,plot=FALSE)$mids),
                      stopAt = c(graphics::hist(emDims,plot=FALSE)$mids,max(emDims)+.5),
                      f=factor(seq_along(c(.5, graphics::hist(emDims,plot=FALSE)$mids))%%2))

    #tmi <-  nonlinearTseries::mutualInformation(y, lag.max = maxLag, n.partitions = , do.plot = FALSE)

    tmi <- mif(data.frame(y),lags = 1:maxLag)

    dfMI <- data.frame(emDelay = as.numeric(names(tmi)),
                       ami     = as.numeric(tmi))

    Ncol <- length(emLags$selection.method[!is.na(emLags$lag)])
    myPal <- getColours(Ncol)
    myPalLag <- myPal
    names(myPalLag) <- emLags$selection.method[!is.na(emLags$lag)]
    myPalNn <- myPal
    names(myPalNn) <- emLags$lag[!is.na(emLags$lag)]

    gNdims <- ggplot2::ggplot(df, ggplot2::aes_(y = ~Nn.pct, x = ~emDim, colour = ~emLag)) +
      ggplot2::geom_rect(ggplot2::aes_(xmin = ~startAt, xmax = ~stopAt, fill = ~f), ymin = -Inf, ymax = Inf, data = dfs, inherit.aes = FALSE) +
      ggplot2::scale_fill_manual(values = scales::alpha(c("grey", "white"),.2), guide="none") +
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
                     legend.position = c(.9, .8),
                     legend.background = element_rect(colour = "grey10",fill = "grey90"),
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
                     legend.background = element_rect(colour = "grey10",fill = "grey90"),
                     legend.title = element_text(face = "bold"),
                     legend.key = element_rect(colour = "grey90", fill = "grey90"),
                     #panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     panel.grid.major.y = element_blank(),
                     panel.grid.minor.y = element_blank()

      )

    g <- gridExtra::grid.arrange(gDelay, gNdims, ncol=1, nrow=2)
    #grid::grid.newpage()
    grid::grid.draw(g)

  } else {
    g <- NA
  }

  return(invisible(list(optimLag  = opLag,
                        optimDim  = opDim,
                        optimRow  = opt,
                        optimData = df,
                        diagPlot = g)))
}



#' Estimate ROC radius
#'
#'  Experimental.
#'
#' @param y y
#' @param emRad radius
#' @param emDim embedding Dims
#' @param emLag embedding Lag
#' @param noiseLevel noise Level
#' @param standardise Standardise y? Choose from "mean.sd","median.mad","none".
#' @param noiseType Use a Normal distribution of uniform distribution for noiselevels
#'
#' @family Estimate Recurrence Parameters
#'
#' @return data frame for ROC
#' @export
#'
#' @keywords internal
est_parameters_roc <- function(y, emRad, emDim=1, emLag=1, noiseLevel=.75, standardise = c("mean.sd","median.mad","none")[3], noiseType = c("normal","uniform")[1]){
  y <- dplyr::case_when(
    standardise == "mean.sd"   ~ ts_standardise(y, type="mean.sd"),
    standardise == "median.sd" ~ ts_standardise(y, type="median.mad"),
    standardise == "none"      ~ y
  )

  yn <- dplyr::case_when(
    noiseType == "normal"  ~ stats::rnorm(NROW(y), mean=round(mean(y, na.rm = TRUE),3), sd=stats::sd(y,na.rm = TRUE)),
    noiseType == "uniform" ~ sign(stats::rnorm(1))*stats::runif(NROW(y), min=floor(min(y, na.rm = TRUE)), max = ceiling(max(y,na.rm = TRUE)))
  )

  noise_out   <- rp_cl(yn, emRad = emRad, emDim=emDim, emLag=emLag)
  measure_out <- rp_cl((y  + noiseLevel * yn), emRad = emRad, emDim=emDim, emLag=emLag)

  return(cbind.data.frame(radius   = emRad,
                          response = c("signal+noise","noise"),
                          rbind(measure_out,noise_out)))
}


#' Estimate embedding lag (tau)
#'
#' A wrapper for [nonlinearTseries::timeLag]
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
#' @family Estimate Recurrence Parameters
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
    lag <-return_error(nonlinearTseries::timeLag(y,
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
#' A wrapper for [nonlinearTseries::estimateEmbeddingDim]
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
#' @family Estimate Recurrence Parameters
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

#
#
# est_RRvsRad <- function(y, RM = NA, emLag = NA, emDim = NA){
#   if(is.na(RM)){
#     RM <- rp(y1 = y, emLag = emLag, emDim = emDim)
#   }
#   df_RR <- plyr::ldply(seq(0,1, by = .05), function(t){
#     suppressMessages(out <- est_radius(RM = RM, targetValue = t, silent = TRUE, maxIter = 500, radiusOnFail = NA))
#     return(data.frame(target = t, RR = out$Measure, emRad = out$Radius))
#   })
#
# }
#


#' Estimate the maximum number of Phases
#'
#' @description Parameter sweep of function [rn_phases] for argument `maxPhases`. Use to check at which value of `maxPhases` no additional phases will be detected.
#'
#' @inheritParams rn_phases
#' @param RN Recurrence matrix
#' @param range Two element vector with minimum and maximum `c(min,max)` number of phases to check.
#'
#' @return Data frame with maxPhases by detectedPhases
#' @export
#'
est_maxPhases <- function(RN,
                          range = 2:10,
                          minStatesinPhase = 1,
                          maxStatesinPhase = NROW(RN),
                          useDegree = FALSE,
                          inverseWeight = TRUE,
                          cleanUp = TRUE,
                          removeSingularities = TRUE){

  maxPhaseOut <- plyr::ldply(range[1]:max(range), function(p){
    tmp <- rn_phases(RN, minStatesinPhase = minStatesinPhase, maxStatesinPhase = maxStatesinPhase, inverseWeight = inverseWeight, cleanUp = cleanUp, removeSingularities = removeSingularities, maxPhases = p)
    data.frame(Phases = max(tmp$phase_number), maxPhases = p)
    })

 g <- ggplot(maxPhaseOut, aes(y = Phases, x = maxPhases)) +
    geom_line() +
    geom_point() +
    theme_bw()

 print(g)


  return(invisible(maxPhaseOut))
}
