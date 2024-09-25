# package casnet ----
#
# Recurrence Networks ----
#
# rn functions

#' Create a Recurrence Network Matrix
#'
#' This function serves as a wrapper for function `rp()`, it will add some attributes to the matrix related to network representation. These attributes will be used to decide which network type to generate (e.g. undirected, directed, weighted, etc.)
#'
#' @inheritParams rp
#' @param directed Should the matrix be considered to represent a directed network? (default = `FALSE`)
#' @param cumulative To make the network represent cumulative time, set `directed = TRUE` and `cumulative = TRUE`. This will set the upper triangle of the recurrence matrix to `0` and ensures that the network edges represent recurrent values that have occurred in the `past` relative to the current observed value (node). If `directed = FALSE` the argument is ignored (default = `TRUE`).
#' @param weighted Should the matrix be considered to represent a weighted network? (default = `FALSE`)
#' @param weightedBy After setting values smaller than `emRad` to `0`, what should the recurrent values represent? The default is to use the state space similarity (distance/proximity) values as weights (`"si"`). Other option are `"rt"` for *recurrence time* and `"rf"` for *recurrence time frequency*, Because vertices represent time points in \eqn{\epsilon}-thresholded recurrence networks, a difference of two vertex-indices represents duration. If an edge `e1` connects `v1` and `v10` then the *recurrence time* will be the difference of the vertex indices, `9`, and the *recurrence time frequency* will be `1/9`.
#' @param rescaleWeights If set to `TRUE` and `weighted = TRUE`, all weight values will be rescaled to `[0,1]`, where `0` means no recurrence relation and `1` the maximum weight value.
#' @param fs Sample frequency: A numeric value interpreted as the `number of observed samples per unit of time`. If the weights represent recurrence times (`"rt"`), they will be divided by the value in `fs`. If the weights represent recurrence time frequencies (`"rf"`), they will be multiplied by the value of `fs` (default = `NA`)
#' @param returnGraph Return an [igraph::igraph()] object (default = `FALSE`)
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
               theiler = 0,
               directed = FALSE,
               cumulative = TRUE,
               weighted = FALSE,
               weightedBy = c("none","si","rt","rf")[1],
               rescaleWeights = FALSE,
               fs = NA,
               to.ts = NULL,
               order.by = NULL,
               to.sparse = FALSE,
               method = c("Euclidean","max","SBD")[1],
               targetValue = .05,
               returnGraph = FALSE,
               doPlot = FALSE,
               doEmbed = TRUE,
               silent = TRUE,
               ...){


  if(weightedBy %in% c("si","rt","rf")){
    weighted <- TRUE
    if(!silent){
      message("`weighted` was set to TRUE due to the value of `weightedBy`")
    }
  }

  if(weighted){
    if(is.na(emRad%00%NA)){
      emRad <- NA
      if(!silent){
        message("`emRad` was set to NA due to the value of `weighted`")
      }
    }
    if(!weightedBy %in% c("si","rt","rf")){
      weightedBy <- "si"
      message("Set `weightedBy` to 'si' due to the value of `weighted`")
    }
  }

  dmat <- rp(y1 = y1, y2 = y2,
             emDim = emDim, emLag = emLag, emRad = emRad, theiler = theiler,
             to.ts = to.ts, order.by = order.by, to.sparse = to.sparse,
             weighted = weighted, targetValue = targetValue,
             method = method, doPlot = FALSE, doEmbed = doEmbed, silent = silent)

  if(!doEmbed){
    emRad <- NA
  }

  if(is.null(y2)){
    y2 <- y1
    attributes(y2) <- attributes(y1)
  }

  if(to.sparse){
    attributes(dmat)$directed <- directed
    attributes(dmat)$cumulative <- cumulative
    #attributes(dmat)$weighted <- weighted
  } else {
    attr(dmat,"directed") <- directed
    attr(dmat,"cumulative") <- cumulative
    #attr(dmat,"weighted") <- weighted
  }

  if(attr(dmat,"AUTO")){
    mode <- "upper"
  }

  if(directed){
    mode <- "directed"
    # if(cumulative){
    #   dmat[upper.tri(dmat)] <- 0
    # }
  }

  #
  if(!is.null(emRad)){
    if(weighted){
      if(!weightedBy%in%c("si","rt","rf")){
        stop("Invalid string in argument weightedBy!")
      }

      grW <- igraph::graph_from_adjacency_matrix(dmat, weighted = weighted, mode = mode, diag = TRUE) #includeDiagonal)
      edgeFrame             <- igraph::as_data_frame(grW,"edges")
      E(grW)$rec_distance   <- E(grW)$weight
      E(grW)$rec_time       <- abs(edgeFrame$to-edgeFrame$from)
      E(grW)$rec_timefreq   <- 1/(edgeFrame$weight+1) #.Machine$double.eps)

      if(is.numeric(fs)){
        if(weightedBy=="rf"){edgeFrame$weight <- edgeFrame$rectime * fs}
        if(weightedBy=="rt"){edgeFrame$weight <- edgeFrame$rectime / fs}
      }

      switch(weightedBy,
             si = E(grW)$weight <- E(grW)$rec_distance,
             rt = E(grW)$weight <- E(grW)$rec_time,
             rf = E(grW)$weight <- E(grW)$rec_timefreq
      )

      # for(r in 1:NROW(edgeFrame)){
      #   dmat[edgeFrame$from[r],edgeFrame$to[r]] <- 1/edgeFrame$rectime[r]
      #   if(attr(dmat,"AUTO")){
      #     dmat[edgeFrame$to[r],edgeFrame$from[r]] <- 1/edgeFrame$rectime[r]
      #     }
      #   }

      tmp <- as_adjacency_matrix(grW, type = "both", sparse = to.sparse, attr = "weight")
      #tmp <- bandReplace(tmp,0,0,0)
      dmat <- rp_copy_attributes(source = dmat, target = tmp)
      rm(tmp)

      if(rescaleWeights==TRUE){
        dmat <- dmat/max(dmat, na.rm = TRUE)
      }
    } else {
      weighted <- NULL
    }
  } #if is.null(emRad)

  if(doPlot){

    if(doPlot){
      dotArgs  <- formals(rn_plot)
      nameOk   <- rep(TRUE,length(dotArgs))
      if(...length()>0){
        dotArgs <- list(...)
        nameOK  <- names(dotArgs)%in%methods::formalArgs(rn_plot)
        # Plot with defaults
        if(!all(nameOK)){
          dotArgs    <- formals(rn_plot)
          nameOk <- rep(TRUE,length(dotArgs))
        }
      }
      dotArgs$RN <- dmat
      do.call(rn_plot, dotArgs[nameOk])
    }

    # dotArgs <- list(...)
    #
    # if(is.null(dotArgs)){
    #   dotArgs<- formals(rn_plot)
    #   }
    #
    # nameOk  <- names(dotArgs)%in%methods::formalArgs(rn_plot)
    # # Plot with defaults
    # if(!all(nameOk)){
    #   dotArgs <- formals(rn_plot)
    #   nameOk  <- rep(TRUE,length(dotArgs))
    # }
    #
    # dotArgs$RN <- dmat
    #
    # do.call(rn_plot, dotArgs[nameOk])
  }

  if(directed){
    if(cumulative){
      dmat[upper.tri(dmat)] <- 0
    }
  }

  if(returnGraph){
    g  <- igraph::graph_from_adjacency_matrix(dmat, weighted = weighted, mode = mode) #, diag = includeDiagonal)
    # V(g)$y1 <- y1[1:igraph::vcount(g)]
    # V(g)$y2 <- y2[1:igraph::vcount(g)]

    return(list(RN = dmat,
                g  = g)
    )
  } else {
    return(dmat)
  }
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
                    radiusValue = NA,
                    title = "", xlabel = "", ylabel="",
                    plotSurrogate = NA,
                    returnOnlyObject = FALSE){

  if(attributes(RN)$directed&attributes(RN)$cumulative&is.na(attributes(RN)$emRad)){
    warning("Plotting an unthresholded, directed, distance plot may not yield a sensible result!")
  }

  rp_plot(RM = RN,
          plotDimensions  = plotDimensions,
          plotMeasures    = plotMeasures,
          plotRadiusRRbar = FALSE,
          drawGrid = drawGrid,
          markEpochsLOI = markEpochsLOI,
          radiusValue = radiusValue,
          title = title, xlabel = xlabel, ylabel = ylabel,
          plotSurrogate = plotSurrogate,
          returnOnlyObject = returnOnlyObject)

}




#' Strength versus Degree scaling relation
#'
#' Calculates the Recurrence Rate versus Recurrence Time power-law
#'
#' @param g an igraph object representing a weighted Recurrence  Network
#' @param mode Evaluate the "in", "out" degree, or "all" edges (default = `"all"`)
#' @param doPlot Plot the scaling relation? (default = TRUE)
#'
#' @return A data frame with local vertex strength and vertex degree, including
#'
#' @export
#'
#' @examples
#'
#' y  <- rnorm(100)
#' RN <- rn(y, emLag=1, emDim=3, emRad=NA, weighted = TRUE, weightedBy = "rt", returnGraph = TRUE)
#' rn_strengthDist(RN$g)
#'
rn_strengthDist <- function(g, mode = c("in","out","all")[3], doPlot = TRUE){
  xD <- igraph::degree(g, mode = mode, loops = FALSE)
  yS <- igraph::strength(g, mode = mode, loops = FALSE)

  out <- data.frame(
    xDegree   = xD[xD>0],
    yStrength = yS[xD>0],
    xDegree_log10   = log10(xD[xD>0]),
    yStrength_log10 = log10(yS[xD>0])
  )

  PL <- lm(out$yStrength_log10~out$xDegree_log10)

  out$PowerLaw <- predict(PL)
  out$PowerLawExponent <-  PL$coefficients[2]

  if(doPlot){
    gPL <- ggplot(out, aes_(y = ~yStrength_log10, x = ~xDegree_log10)) +
      geom_point() +
      geom_line(aes_(x = ~xDegree_log10, y = ~PowerLaw), colour = "steelblue3", size = 1) +
      scale_x_continuous("Vertex Degree (log10)") +
      scale_y_continuous("Vertex Strength (log10)") +
      theme_bw()
    print(gPL)
  }


  return(out)
}


#' Recurrence Network Measures
#'
#' @param g An igraph object. If V(g)$name is set the labels will be returned in a column.
#' @param cumulative Only consider out-degree.
#' @param silent Siletn(ish) mode
#'
#' @return A list with data frames with common vertex, edge amnd global network measures.
#'
#' @export
#'
#'
rn_measures <- function(g, cumulative = TRUE, silent = TRUE){

  checkPkg("DirectedClustering")
  checkPkg("brainGraph")

  if(is.null(igraph::V(g)$name)){
    V(g)$name <- paste(1:igraph::vcount(g))
  }

  vertex_prop <- data.frame(vertex_name = V(g)$name,
                            LocalClustering=rep(NA,igraph::vcount(g)),
                            DegreeDensity=rep(NA,igraph::vcount(g)),
                            StrengthDensity=rep(NA,igraph::vcount(g)),
                            LocalCloseness=rep(NA,igraph::vcount(g)),
                            LocalEfficiency=rep(NA,igraph::vcount(g)),
                            LocalBetweenness=rep(NA,igraph::vcount(g))
  )


  edge_prop   <- list()

  graph_prop  <- data.frame(EdgeDensity=NA,
                            MeanStrengthDensity = NA,
                            GlobalClustering=NA,
                            NetworkTransitivity=NA,
                            AveragePathLength=NA,
                            GlobalEfficiency=NA
  )


  efType  <- "nodal"
  degType <- "total"

  if(igraph::is_weighted(g)){
    clType <- "barrat"
    vertex_prop$StrengthDensity   <- igraph::strength(g, loops = FALSE, mode = degType)/(igraph::vcount(g)-1)
  } else{
    clType <- "local"
    vertex_prop <- vertex_prop[,-("StrengthDensity"%ci%vertex_prop)]
  }

  if(igraph::is_directed(g)){
    if(cumulative){
      dmat <- igraph::as_adjacency_matrix(g, sparse = FALSE)
      dmat[upper.tri(dmat)] <- 0
      vertex_prop$LocalClustering  <- DirectedClustering::ClustBCG(dmat, "directed", isolates = "zero")$outCC
      degType <- "out"
      efType <- "local"
    } else {
      vertex_prop$LocalClustering  <- DirectedClustering::ClustBCG(igraph::as_adjacency_matrix(g, sparse = FALSE), "directed", isolates = "zero")$totalCC
    }
  } else {
    vertex_prop$LocalClustering  <- igraph::transitivity(g, type = clType, isolates = "zero")
  }

  vertex_prop$DegreeDensity    <- igraph::degree(g, loops = FALSE)/(igraph::vcount(g)-1)
  vertex_prop$LocalCloseness   <- igraph::closeness(g, normalized = TRUE, mode = degType)
  vertex_prop$LocalEfficiency  <- brainGraph::efficiency(g,type = efType, use.parallel = FALSE, A = igraph::as_adjacency_matrix(g, sparse = FALSE))
  vertex_prop$LocalBetweenness <- igraph::betweenness(g, directed = igraph::is_directed(g), normalized = TRUE)


  # igraph::transitivity(as.undirected(GsWD[[g]]), type = "barrat", isolates = "zero")

  # outListWD$locClose[[g]]  <- igraph::estimate_closeness(GsWD[[g]], mode = "out", normalized = TRUE, cutoff = 2)
  #outListWD$locEfficiency[[g]] <- brainGraph::efficiency(GsWD[[g]], type = "local", use.parallel = FALSE, A = igraph::as_adjacency_matrix(GsWD[[g]],edges = TRUE, sparse = FALSE))
  #outListWD$locBetween[[g]] <- igraph::betweenness(GsWD[[g]], directed = TRUE, normalized = TRUE)

  graph_prop$EdgeDensity      <- mean(vertex_prop$DegreeDensity, na.rm = TRUE)
  if(igraph::is_weighted(g)){
    graph_prop$MeanStrengthDensity  <- mean(vertex_prop$StrengthDensity, na.rm=TRUE)
  }
  graph_prop$GlobalClustering    <- mean(vertex_prop$LocalClustering)
  graph_prop$NetworkTransitivity <- igraph::transitivity(g, type = "global", isolates = "zero")
  graph_prop$AveragePathLength   <- igraph::mean_distance(g, directed = FALSE, unconnected = FALSE)
  #mean(igraph::closeness(g, normalized = TRUE)^-1)
  graph_prop$GlobalEfficiency    <-(mean(brainGraph::efficiency(g, type = "nodal", use.parallel = FALSE)))^-1

  if(!silent){
    cat("\n~~~o~~o~~casnet~~o~~o~~~\n")
    cat("\nGlobal Network Measures\n\n")
    print(format(graph_prop))
    cat("\n~~~o~~o~~casnet~~o~~o~~~\n")
  }

  return(list(vertex_measures = vertex_prop,
              edge_measures  = edge_prop,
              graph_measures = graph_prop)
  )
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

  #diagonal <- attributes(RN)$includeDiagonal
  weighted <- NULL
  if(attributes(RN)$weighted){weighted <- TRUE}

  g1 <- igraph::graph_from_adjacency_matrix(RN, mode="directed", diag = FALSE, weighted = weighted)

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
  ddata <- dplyr::arrange_(ddata,~size)
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

#' Extract Phases from weighted RN
#'
#'  This function will extract phases (regions of attraction in phase space) based on a weighted `RN` object created with function [rn].
#'  The assumption is that coordinates in state space that are either close in terms of distance, or are re-visited with short recurrence times,
#'  or with high-frequency, are regions of attraction for the system.
#'
#'  The method used for the identification of phases is on the properties of the `RN` object:
#'
#'  - If weighted by distance `"si"`, the inverse distance will be used, which means higher weights correspond to closer states.
#'  - If weighted by recurrence time `"rt"`, the inverse time will be used, which means higher weights correspond to faster recurrence times.
#'
#'  The procedure is as follows:
#'
#'  1. Identify the node with the highest strength
#'  2. Identify the nodes that connect to this node
#'  3. Identify the node with highest strength that does not connect to the node identified in step 1.
#'  4. Repeat until criteria set in `maxPhases`, `minStatesinPhase` and `maxStatesinPhase` are triggered.
#'  5. Clean up. If `cleanUp = TRUE` the remaining states that connect to nodes of one of the identified phases, but not to the node with highest strength identified in step 1 will be added to that phase. Otherwise these nodes will end up in category "Other".
#'
#'
#' @inheritParams make_spiral_graph
#' @inheritParams fd_dfa
#' @param RN A matrix produced by the function [rn]
#' @param maxPhases The maximum number of phases to extract. These will be the phases associated with the highest node degree or node strength. All other recurrent point will be labelled with "Other". If `NA`, the value will be set to `NROW(RN)`, this will return all the potential phases in the data irrespective of their frequency/strength of recurrence  (default = `5`)
#' @param minStatesinPhase A parameter applied after the extraction of phases (limited by `maxPhases`). If any extracted phases do not have a minimum number of `minStatesinPhase` + `1` (= the state that was selected based on node strength), the phase will be removed from the result (default = `1`)
#' @param maxStatesinPhase A parameter applied after the extraction of phases (limited by `maxPhases`). If any extracted phases exceeds a maximum number of `maxStatesinPhase` + `1` (= the state that was selected based on node strength), the phase will be removed from the result (default = `NROW(RN)`)
#' @param useDegree By default, node strength will be used to determine the phases. Set to `TRUE` to use the node degree (default = `FALSE`)
#' @param removeSingularities Will remove states that recur only once (nodes with `degree(g) == 1`) (default = `FALSE`)
#' @param inverseWeight Whether to perform the operation `1/weight` on the edge weights. The default is `TRUE`, if the matrix was weighted by a distance metric (`weightedBy = "si"`) edges with smaller distances (recurring coordinates closer to the current coordinate) have greater impact on the node strength calculation used to select the phases. If the matrix was weighted by recurrence time (`weightedBy = "rt"`) and `inverseWeight = TRUE`, recurrent points with shorter recurrence times will have greater impact on the strength calculation. If `weightedBy = "rf"`, lower frequencies will end up having more impact. (default = `TRUE`)
#' @param cleanUp Try to assign states to phases that were not assigned by the algorithm. If `FALSE`, these states will be added to category "Other" (default = `TRUE`)
#' @param returnCentroid Values can be `"no"`, `"mean.sd"`, `"median.mad"`, `"centroid"`. Any other value than `"no"` will return a data frame with the central tendency and deviation measures for each phase (default = `"no"`)
#' @param doPhaseProfilePlot Produce a profile plot of the extracted phases (default = `TRUE`)
#' @param plotCentroid Plot the centroid requested in `returnCentroid`? (default = `FALSE`)
#' @param space_dims A vector of titles to use for the dimensions. If `NULL` the values will be read from the attribute of `RN`.
#' @param colOrder Should the order of the dimensions reflect the order of the columns in the dataset? If `FALSE` the order will be based on the values of the dimensions observed in the first extracted phase (default = `FALSE`)
#' @param phaseColours Colours for the different phases in the phase plot. If `epochColours` also has a value, `phaseColours` will be used instead (default = `NULL`)
#' @param doSpiralPlot Produce a plot of the recurrence network with the nodes coloured by phases (default = `FALSE`)
#' @param doPhaseSeriesPlot Produce a time series of the phases as they occur with a marginal histogram of their frequency (default = `FALSE`)
#' @param doPhaseSpaceProjectionPlot produce a 2D `umap` projection of the phases (default = `FALSE`)
#' @param excludeOther Should the phase "Other" be excluded from plots? (default = `FALSE`)
#' @param excludeNorec Should the category "No recurrence" be excluded from plots? (default = `TRUE`)
#' @param returnGraph Returns all the graph object objects of the plots that have been produced (default = `FALSE`)
#'
#' @export
#' @return A data frame with information about the phases or a list object with data and graph objects (if requested).
#'
#' @family Distance matrix operations (recurrence network)
#'
#' @examples
#'
#' # Use the ManyAnalysts dataset to create a phase plot with default settings
#' data("manyAnalystsESM")
#' df <- manyAnalystsESM[4:10]
#' RN <- rn(y1 = df, doEmbed = FALSE, weighted = TRUE, weightedBy = "si", emRad = NA)
#'
#' # This returns 6 phases which have minimally 1 state
#' rn_phases(RN, maxPhases = 10, doPhaseProfilePlot = TRUE)
#'
#' # Use min. number of states as the extraction criterion
#' rn_phases(RN, maxPhases = NA, minStatesinPhase = 7, doPhaseProfilePlot = TRUE)
#'
rn_phases <- function(RN,
                      maxPhases = 5,
                      minStatesinPhase = 1,
                      maxStatesinPhase = NROW(RN),
                      useDegree = FALSE,
                      inverseWeight = TRUE,
                      cleanUp = TRUE,
                      returnCentroid = c("no","mean.sd","median.mad","centroid")[1],
                      removeSingularities = FALSE,
                      standardise = c("none","mean.sd","median.mad","unit")[4],
                      returnGraph = FALSE,
                      doPhaseProfilePlot = FALSE,
                      plotCentroid = FALSE,
                      space_dims = NULL,
                      colOrder = FALSE,
                      phaseColours = NULL,
                      doSpiralPlot = FALSE,
                      doPhaseSeriesPlot = FALSE,
                      doPhaseSpacePojectionPlot = FALSE,
                      showEpochLegend = TRUE,
                      epochColours = NULL,
                      epochLabel = "Phase",
                      excludeOther = FALSE,
                      excludeNorec = TRUE){

  checkPkg("igraph")
  checkPkg("invctr")

  if(!attributes(RN)$weighted%00%FALSE){
    stop("RN has to be a weighted recurrence network.")
  }

  if(!attributes(RN)$weightedBy%in%c("si","rt","rf")){
    stop("RN has to be a distance weighted recurrence network.")
  }

  # if(!attributes(RN)$cumulative%00%FALSE){
  #   stop("RN has to be a cumulative recurrence network.")
  # }

  weighted <- attributes(RN)$weighted
  directed <- attributes(RN)$directed%00%FALSE
 if(!directed){
   directed <- "undirected"
 } else {
   directed <- "directed"
 }

  # Creat graph
  gRN <- igraph::graph_from_adjacency_matrix(RN,
                                             weighted = weighted,
                                             mode = directed)

  if(inverseWeight){
    igraph::E(gRN)$weight <- (1/igraph::E(gRN)$weight)%00%0
    }

  if(is.null(igraph::V(gRN)$size)){
    igraph::V(gRN)$size <- igraph::strength(gRN)
  }

  RN_nodes <- data.frame(time = as.numeric(igraph::V(gRN)), strength = igraph::strength(gRN))
  # Separate degree 1 (Singularities)
  Singularities <- RN_nodes %>% dplyr::filter(igraph::degree(gRN)==1)

  nonRecurring <- sum(RN_nodes$strength==0, na.rm=TRUE)

  # Remove degree 0
  RN_nodes <- RN_nodes %>% dplyr::filter(strength>0)

  # Edges
  RN_edges <- igraph::as_data_frame(gRN)

  last <- FALSE
  i <- 0
  nodeID <- list()
  phases <- list()
  strengths <- list()
  igraph::V(gRN)$phase <- NA

  if(is.na(maxPhases)){
    maxPhases <- NROW(RN)
  }

  while(!last){

    #i%++%1
    i <- i+1

    if(i > 1){
      tmp_nodes   <- RN_nodes %>% dplyr::filter(!(time%in%unique(unlist(phases[1:(i-1)]))))
      if(NROW(tmp_nodes)==0){
        last <- TRUE
        break
      }
    } else {
      tmp_nodes   <- RN_nodes
    }

    nodeID[[i]]    <- tmp_nodes$time[which.max(tmp_nodes$strength)]
    strengths[[i]] <- tmp_nodes$strength[which.max(tmp_nodes$strength)]

    tmp_edges      <- RN_edges %>% dplyr::filter(.data$from==nodeID[[i]]|.data$to==nodeID[[i]])
    if(i > 1){
      tmp_edges   <- tmp_edges %>% dplyr::filter(!((.data$from%in%as.numeric(unlist(phases[1:(i-1)])))|(.data$to%in%as.numeric(unlist(phases[1:(i-1)])))))
    }

    phases[[i]] <- unique(c(tmp_edges$from,tmp_edges$to))

    igraph::V(gRN)$phase[phases[[i]]] <- i
    names(phases)[i] <- ifelse(i<10, paste0("Phase 0",i), paste("Phase",i))

    if(i > 1){
      if(any(i>=NROW(RN_nodes), i==maxPhases, NROW(tmp_edges)==0, nodeID[[i]]==nodeID[[i-1]])){
        last <- TRUE
        break
      }
    }
    rm(tmp_edges,tmp_nodes)
  }

  if(length(unique(unlist(phases)))!=length(unlist(phases))){
    warning("Detected duplicate nodes in different phases!")
  }

  # Phases <-  V(gRN)$phase
  # Phases[is.na(Phases)] <- max(Phases, na.rm = TRUE)+1

  phases    <- phases[lengths(phases)!=0]
  nodeID    <- nodeID[which(lengths(phases)!=0)]
  strengths <- strengths[which(lengths(phases)!=0)]

  out <- tidyr::unnest(dplyr::tibble(phase_name = names(phases),
                                     phase_number = as.numeric_discrete(names(phases)),
                                     phase_size = lengths(phases),
                                     maxState_time = as.numeric(nodeID),
                                     maxState_strength = as.numeric(strengths),
                                     states_time = phases,
                                     states_singularity = 0),
                       cols = 6)

  out$states_strength <- igraph::strength(gRN)[out$states_time]

  tmp     <- plyr::ldply(seq_along(phases), function(p) data.frame(phase_name = names(phases)[[p]],
                                                                   states_time = phases[[p]],
                                                                   maxState_time = nodeID[[p]],
                                                                   states_dist2maxState = RN[phases[[p]], nodeID[[p]]],
                                                                   states_singularity = 0))

  # Check nodes and phases
  if(identical(out$states_time,tmp$states_time)){
    out$states_singularity <- tmp$states_singularity
    out$states_dist2maxState <- tmp$states_dist2maxState
  } else {
    warning("Wrong order?")
  }

  if(cleanUp){

    remains <- seq(1,igraph::vcount(gRN))[!(seq(1,igraph::vcount(gRN))%in%out$states_time)]
    remains <- remains[igraph::degree(gRN)[remains]>0]

    remainsList <- list()
    i <- 0

    for(p in unique(out$phase_name)){

      tmp <- out %>% dplyr::filter(phase_name == p)
      tmp <- dplyr::arrange(tmp, dplyr::desc(tmp$states_strength))

      dist_remain <- strength_remain <- nodes_remain <- list()

      for(n in 1:NROW(tmp)){
        nodes_remain[[n]] <- remains[remains%in%unique(c(RN_edges$to[RN_edges$from==tmp$states_time[n]|RN_edges$to==tmp$states_time[n]], RN_edges$from[RN_edges$to==tmp$states_time[n]|RN_edges$from==tmp$states_time[n]]))]
        strength_remain[[n]] <- igraph::strength(gRN)[nodes_remain[[n]]]
      }
      names(nodes_remain) <- tmp$states_time # paste("Node:",tmp$states_time,"| Strength:",tmp$states_strength)

      dist_remain <- plyr::llply(seq_along(nodes_remain), function(d){
        if(length(nodes_remain[[d]])>0){
          RN[nodes_remain[[d]], as.numeric(names(nodes_remain)[d])]
        } else {
          NA
        }
      })

      i <- i+1

      remainsList[[i]] <- tidyr::unnest(dplyr::tibble(phase_name = p,
                                                      phase_number = tmp$phase_number[tmp$phase_name==p][1],
                                                      phase_size = NA,
                                                      maxState_time = as.numeric(names(nodes_remain)),
                                                      maxState_strength = tmp$states_strength,
                                                      states_time = nodes_remain,
                                                      states_strength = strength_remain,
                                                      states_dist2maxState = dist_remain,
                                                      states_singularity = 0),
                                        cols = 6:8)

      rm(tmp)
    }

    remainsOut <- plyr::ldply(remainsList)
    remainsOut <- dplyr::arrange(remainsOut, dplyr::desc(.data$maxState_strength))
    uniIND <- plyr::laply(as.numeric(remains), function(r) which(r==remainsOut$states_time)[1])
    out <- rbind(out, remainsOut[uniIND[!is.na(uniIND)],])
    out <- dplyr::arrange(out,.data$phase_name)

    rm(remainsOut, remainsList, uniIND)
  }

  if(minStatesinPhase>maxStatesinPhase){
    maxStatesinPhase <- minStatesinPhase
  }

  if(!is.null(space_dims)){
    if(length(space_dims)!=NCOL(attributes(RN)$emDims1)){
      warning("Argument space_dims is not equal to number of phase space dimensions.")
    }
    tmp <- attributes(RN)$emDims1[out$states_time, ,drop = FALSE]
    colnames(tmp) <- space_dims
    space_dims <- tmp
  } else {
  if(!is.null(attributes(RN)$emDims1)){
      space_dims <- attributes(RN)$emDims1[out$states_time, ,drop = FALSE]
      colnames(space_dims) <- paste0("dim_", plyr::laply(strsplit(colnames(space_dims),split = "[.]"), function(s) s[[2]]))
      } else {
        stop("Could not find the attribute 'emDims1' on the Recurrence Matrix.")
      }
    }
  out <- cbind(out,space_dims)
  rm(space_dims)

  # Check if we missed some recurrent states due to the stopping parameters and create an "Other" category
  ltime <- seq(1,igraph::vcount(gRN))[!seq(1,igraph::vcount(gRN))%in%sort(out$states_time)]
  #ltime <- ltime[!ltime%in%Singularities$time]
  # These time points could contain non-recurring points, i.e. distance 0 to all other points
  ltime <- ltime[plyr::laply(ltime, function(p) sum(RN[1:NCOL(RN),p]))>0]

  if(NROW(ltime)!=0){
    lost_time <- data.frame(matrix(NA,nrow=NROW(ltime),ncol=NCOL(out),dimnames = list(NULL,colnames(out))))
    lost_time$states_time <- ltime
    lost_time[,(("states_dist2maxState"%ci%lost_time)+1):NCOL(lost_time)] <- attributes(RN)$emDims1[ltime,]
    lost_time$phase_name <- "Other"
    lost_time$phase_number <- max(out$phase_number, na.rm = TRUE)+1
    lost_time$phase_size <- length(ltime)
    lost_time$states_strength <- igraph::strength(gRN)[lost_time$states_time]
    lost_time$maxState_strength <- max(lost_time$states_strength, na.rm = TRUE)
    lost_time$maxState_time <-  lost_time$states_time[which.max(lost_time$states_strength)]
    lost_time$states_singularity <- 0
    lost_time$states_dist2maxState <- NA

    out <- rbind(out,lost_time)
  }

  # Finally add the nonrecurring states
  norec <- seq(1,igraph::vcount(gRN))[!seq(1,igraph::vcount(gRN))%in%sort(out$states_time)]
  #norec <- norec[!norec%in%Singularities$time]

  if(NROW(norec)!=0){
    empty <- data.frame(matrix(NA,nrow=NROW(norec),ncol=NCOL(out),dimnames = list(NULL,colnames(out))))
    empty$states_time <- norec
    empty[,(("states_dist2maxState"%ci%empty)+1):NCOL(empty)] <- attributes(RN)$emDims1[norec,]
    empty$phase_name <- "No recurrence"
    empty$phase_number <- max(out$phase_number, na.rm = TRUE)+1
    empty$phase_size <- length(norec)
    empty$states_strength <- .Machine$double.eps

    out <- rbind(out,empty)
  }


  if(plotCentroid&(returnCentroid=="no")){
    message("Specify which centroid type to return in argument 'returnCentroid'")
    plotCentroid <- FALSE
  }

  if(returnCentroid%in%"centroid"&!(standardise%in%"mean.sd")){
    standardise <- "mean.sd"
    warning("Changed value of standardise to 'mean.sd' which is required if returnCentroid is set to 'centroid'.")
  }

  if(any(standardise%in%c("mean.sd","median.mad","unit"))){

    if(standardise%in%"unit"){

      space_dims <- signif(elascer(out[,(("states_dist2maxState"%ci%out)+1):NCOL(out)]))

    } else {

      if(returnCentroid%in%"centroid"){
        checkPkg("dtwclust")
        NAind      <- which(stats::complete.cases(out[,(("states_dist2maxState"%ci%out)+1):NCOL(out)]))
        # #tmp        <- t(out[NAind,(("states_dist2maxState"%ci%out)+1):NCOL(out)])
        # tmp_dims   <- signif(plyr::colwise(dtwclust::zscore)(data.frame(out[,(("states_dist2maxState"%ci%out)+1):NCOL(out)])))
        # space_dims <- matrix(NA, nrow = NROW(out), ncol = length((("states_dist2maxState"%ci%out)+1):NCOL(out)))
        # space_dims[NAind,] <- tmp_dims
        # colnames(space_dims) <- rownames(tmp)

        space_dims <- signif(dtwclust::zscore(data.frame(out[,(("states_dist2maxState"%ci%out)+1):NCOL(out)]), multivariate = TRUE))
        warning("Centroid calculation by dtwclust: Any NA values will be set to 0!")

      } else {
        space_dims <- signif(plyr::colwise(ts_standardise, type = standardise, adjustN = TRUE)(data.frame(out[,(("states_dist2maxState"%ci%out)+1):NCOL(out)])))
      }
    }
    out[,(("states_dist2maxState"%ci%out)+1):NCOL(out)] <- space_dims
    rm(space_dims)
  }

  out <- out %>% dplyr::group_by(.data$phase_name) %>% dplyr::mutate(phase_size = dplyr::n()) %>% dplyr::ungroup()
  out <- out %>% dplyr::arrange(.data$states_time)
  out$states_singularity[out$states_time %in% Singularities$time] <- 1

  # Check boundaries
  keep <- names(phases)[lengths(phases)%[]%c(minStatesinPhase,maxStatesinPhase)]
  out  <- out[out$phase_size%[]%c(minStatesinPhase,maxStatesinPhase),]

  out$phaseName <- paste(out$phase_name," |  N =",out$phase_size)
  #out$phaseName <- factor(out$phaseName,unique(out$phaseName), ordered = TRUE)

  #colIND <- which(arr.ind = TRUE, sort(colIND, decreasing = TRUE))


  # CENTROID ----
  if((returnCentroid)%in%c("mean.sd","median.mad","centroid")){

    if((returnCentroid)%in%c("mean.sd","median.mad")){

      outMeans <- plyr::ldply(unique(out$phase_name), function(p){
        tmp <- out %>% dplyr::ungroup() %>% dplyr::filter(.data$phase_name==p) %>% dplyr::select(dplyr::starts_with("dim_"))
        tmp$phase_name <- p
        if(returnCentroid=="mean.sd"){
          Mean <- data.frame(phase_name = p, Mean = colMeans(tmp[,-NCOL(tmp)], na.rm = TRUE))
          SD <- plyr::colwise(sd, na.rm=TRUE)(tmp[,-NCOL(tmp)])
          Mean$SD      <- as.numeric(t(SD))
        }
        if(returnCentroid=="median.mad"){
          Mean <- data.frame(phase_name = p, Median = t(plyr::colwise(stats::median, na.rm=TRUE)(tmp[,-NCOL(tmp)])))
          SD <- plyr::colwise(stats::mad, na.rm=TRUE)(tmp[,-NCOL(tmp)])
          Mean$MAD        <- as.numeric(t(SD))
        }
        Mean$Dimension  <- rownames(Mean)
        Mean$phase_size <- unique(out$phase_size[out$phase_name==p])
        return(Mean)
      })

    }

    if(returnCentroid%in%"centroid"){
      #checkPkg("dtwclust")

      NAind <- which(!stats::complete.cases(out))
      tmpPhases <- out[stats::complete.cases(out),]

      outMeansList <- list()
      pnames <- unique(tmpPhases$phase_name)
      for(p in seq_along(pnames)){
        tmp <- tmpPhases %>%
          dplyr::ungroup() %>%
          dplyr::filter(.data$phase_name==pnames[p]) %>%
          dplyr::select(dplyr::starts_with("dim_"))

        centroid <- dtwclust::shape_extraction(data.frame(tmp), znorm = FALSE)
        outMeansList[[p]] <- data.frame(phase_name = pnames[p],
                                        Mean = centroid,
                                        SD = 0,
                                        Dimension = colnames(tmp),
                                        phase_size = unique(tmpPhases$phase_size[tmpPhases$phase_name==pnames[p]]))
      }

      outMeans <- plyr::ldply(outMeansList)
    }


    if(excludeOther){
      outMeans <- outMeans %>% dplyr::filter(.data$phase_name != "Other")
    }
    if(excludeNorec){
      outMeans <- outMeans %>% dplyr::filter(.data$phase_name != "No recurrence")
    }

    outMeans$ymin <- outMeans[,c("Median","Mean")%ci%outMeans] - outMeans[,c("MAD","SD")%ci%outMeans]
    outMeans$ymax <- outMeans[,c("Median","Mean")%ci%outMeans] + outMeans[,c("MAD","SD")%ci%outMeans]
    #outMeans$states_time <- outMeans
    outMeans$Dimension <- forcats::fct_inorder(outMeans$Dimension)
    outMeans$phaseName <- paste(outMeans$phase_name," |  N =",outMeans$phase_size)
    #outMeans$phaseName <- factor(outMeans$phaseName,unique(outMeans$phaseName), ordered = TRUE)
  }

   # spiral plot ----
  if(doSpiralPlot){

      # Phases <-  V(gRN)$phase
      # Phases[is.na(Phases)] <- max(Phases, na.rm = TRUE)+1

      if(!is.na(epochColours%00%NA)){
        markEpochsBy <- out$phase_name
        epochColours <- getColours(Ncols = length(unique(markEpochsBy)))
        names(epochColours) <- markEpochsBy
      }

      gg <- casnet::make_spiral_graph(g = gRN,
                                      type = "Euler",
                                      arcs = 4,
                                      markTimeBy = TRUE,
                                      showEpochLegend = TRUE,
                                      markEpochsBy = markEpochsBy,
                                      epochColours = epochColours,
                                      epochLabel = "Phase")
      print(gg)
    }


  bCentroid <- FALSE
  if(returnCentroid!="no"){
    bCentroid <- TRUE
  }

  if(bCentroid){
    listOut <- list(phaseSeries = out, phaseCentroids = outMeans)
  } else {
    listOut <- out
  }

  attr(listOut,"rn_phases") <- list(maxPhases = maxPhases,
                                    minStatesinPhase = minStatesinPhase,
                                    maxStatesinPhase = maxStatesinPhase,
                                    cleanUp = cleanUp,
                                    returnCentroid = returnCentroid,
                                    removeSingularities = removeSingularities,
                                    standardise = standardise)

  # phaseprofile ----
  if(doPhaseProfilePlot){
    pp <- plotRN_phaseProfile(listOut,
                              plotCentroid = plotCentroid,
                              excludeOther = excludeOther,
                              excludeNorec = excludeNorec,
                              returnGraph = TRUE)
  }

  if(doPhaseSeriesPlot){
    ps <- plotRN_phaseSeries(listOut,
                             excludeOther = excludeOther,
                             excludeNorec = excludeNorec,
                             returnGraph = TRUE)
  }


  if(returnGraph){
    listOut <- list(phases = list(phaseSeries = out, phaseCentroids = outMeans)[TRUE, Centroid],
                    plots  = list(phaseProfile = pp, phaseSeries = ps, spiralPlot = gg)[c(doPhaseProfilePlot, doPhaseSeriesPlot, doSpiralPlot)])
  }

  return(invisible(listOut))
}


#' Create transition network
#'
#' @inheritParams rn_phases
#' @param phaseSequence A vector with names or numbers that represent a sequence of phases or order parameter dynamics. If a named numeric vector is passed, the names attribute will be used to represent the phases. If the variable `phase_name`, generated by [rn_phases] is used, the arguments `excludeOther` and `excludeNorec` wil be evaluated.
#' @param threshold Provide a threshold for the relative frequencies. Values below the threshold will be set to `0`. Pass `NA` to not use a threshold (default = `NA`)
#' @param doMatrixPlot default(`TRUE`)
#' @param doNetworkPlot default(`TRUE`)
#' @param returnGraph Return an [igraph::igraph()] object (default = `FALSE`)
#'
#' @return A transition matrix based on relative frequencies
#'
#' @export
#'
#' @examples
#'
#' # This will output the transition matrix, a plot of the matrix and a transition network plot.
#' y <- c("Happy", "Happy", "Sad", "Neutral", "Neutral", "Angry", "Sad", "Sad", "Neutral")
#' TM <- rn_transition(y)
#'
#' @author Fred Hasselman
#' @author Matti Heino
#'
rn_transition <- function(phaseSequence, threshold = NA, doMatrixPlot = TRUE, doNetworkPlot = TRUE, excludeOther = FALSE, excludeNorec = TRUE, returnGraph = FALSE){

  if(excludeOther){
    phaseSequence <- phaseSequence %>% tibble::as_tibble() %>% dplyr::filter(phaseSequence != "Other")
    phaseSequence <- phaseSequence$value
  }
  if(excludeNorec){
    phaseSequence <- phaseSequence %>% tibble::as_tibble() %>% dplyr::filter(phaseSequence != "No recurrence")
    phaseSequence <- phaseSequence$value
  }

  phaseSequence <- as.numeric_discrete(phaseSequence, sortUnique = TRUE)
  phase_names <- names(phaseSequence)

  df <- data.frame(from = dplyr::lag(phaseSequence,1), to =  phaseSequence,
                   from_name = dplyr::lag(phase_names,1), to_name = phase_names) %>%
    dplyr::slice(-1) %>%
    dplyr::group_by(.data$from, .data$to) %>%
    dplyr::summarise(N = dplyr::n(),
                     from_name = dplyr::first(.data$from_name),
                     to_name = dplyr::first(.data$to_name)) %>%
    dplyr::mutate(freq = (.data$N / sum(.data$N, na.rm = TRUE))) %>%
    dplyr::select(-.data$N) %>%
    dplyr::ungroup()

  if(!is.na(threshold)){
    df$freq[df$freq<=threshold] <- 0
  }

  TM <- Matrix::sparseMatrix(x = df$freq,
                             i = df$from,
                             j = df$to,
                             dims = rep(length(unique(phase_names)),2),
                             dimnames = list(sort(unique(phase_names)), sort(unique(phase_names))))
  attr(TM,"rows") <- "from"
  attr(TM,"cols") <- "to"

  mp <- NA
  if(doMatrixPlot){

    df$freq <- signif(df$freq,3)
    df$labels <- ifelse(is.na(df$freq),"0", df$freq)
    mp <-  ggplot(df, aes_(x = ~to_name,
                           y = ~from_name)) +
      geom_tile(aes_(fill = ~freq),colour = "black",size = 0.4) +
      geom_text(aes_(label = ~labels)) +
      scale_fill_gradient(low = "#ECE3CD",
                          high = "#A65141",
                          na.value = "grey",
                          guide = "none") +
      theme_bw() +
      scale_y_discrete(expand = c(0, 0)) +
      scale_x_discrete(expand = c(0, 0)) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            axis.text.y = element_text(angle = 30, hjust = 1),
            legend.position = "right",
            legend.title = element_blank(),
            panel.grid = element_blank(),
            panel.background = element_rect(fill="grey80")) +
      coord_equal() +
      labs(x = "To ...", y = "From ...")

    print(mp)
  }

  gg <- NA
  if(doNetworkPlot){

    gg <- igraph::graph_from_adjacency_matrix(TM, mode = "directed", weighted = TRUE, diag = TRUE)
    gg$layout <- igraph::layout_as_tree(gg,mode = "all", circular = TRUE)

    gg <- plotNET_prep(gg, nodeSize = "strength", rescaleSize = c(10,20), edgeColour = TRUE, doPlot = FALSE, removeSelfLoops = FALSE)

    E(gg)$curved <- .3
    E(gg)$loop.angle <- pi/12
    E(gg)$arrow.width <- 1
    E(gg)$arrow.size <- .6

    if(!is.na(threshold)){
      gg <-  delete_edges(gg, E(gg)[E(gg)$weight == 0])
    }

    plot(gg)

  }


  if(returnGraph){
    list(TransitionMatrix = TM,
         MatrixPlot = mp,
         TransitionNetwork = gg)[TRUE,doMatrixPlot,doNetworkPlot]
  } else {
    return(TM)
  }
}

