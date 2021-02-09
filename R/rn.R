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
#' @param cumulative To make the network represent cumulative time, set `directed = TRUE` and `cumulative = TRUE`. This will set the upper triangle of the recurrence matrix to `0` and ensures that the network edges represent recurrenct values that have occurred in the `past` relative to the current opbserved value. If `directed = FALSE` the argument is ignored (default = `TRUE`).
#' @param weighted Should the matrix be considered to represent a weighted network? (default = `FALSE`)
#' @param weightedBy After setting values smaller than `emRad` to `0`, what should the recurrent values represent? The default is to use the state space similarity (distance/proximity) values as weights (`"si"`). Other option are `"rt"` for *recurrence time* and `"rf"` for *recurrence time frequency*, Because vertices represent time points in \eqn{\epsilon}-thresholded recurrence networks, a difference of two vertex-indices represents duration. If an edge `e1` connects `v1` and `v10` then the *recurrence time* will be the difference of the vertex indices, `9`, and the *recurrence time frequency* will be `1/9`.
#' @param rescaleWeights If set to `TRUE` and `weighted = TRUE`, all weight values will be rescaled to `[0,1]`, where `0` means no recurrence relation and `1` the maximum weight value.
#' @param fs Sample frequency: A numeric value interpreted as the `number of observed samples per unit of time`. If the weights represent recurrence times (`"rt"`), they will be divided by the value in `fs`. If the weights represent recurrence time frequencies (`"rf"`), they will be multiplied by the value of `fs` (default = `NA`)
#' @param includeDiagonal Should the diagonal of the matrix be included when creating the network (default = `FALSE`)
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
               directed = FALSE,
               cumulative = TRUE,
               weighted = FALSE,
               weightedBy = c("si","rt","rf")[1],
               rescaleWeights = FALSE,
               fs = NA,
               includeDiagonal = FALSE,
               to.ts = NULL,
               order.by = NULL,
               to.sparse = FALSE,
               method = "Euclidean",
               targetValue = .05,
               returnGraph = FALSE,
               doPlot = FALSE,
               doEmbed = TRUE,
               silent = TRUE,
               ...){

  dmat <- rp(y1 = y1, y2 = y2,
             emDim = emDim, emLag = emLag, emRad = emRad,
             to.ts = to.ts, order.by = order.by, to.sparse = to.sparse,
             weighted = weighted,targetValue = targetValue,
             method = method, doPlot = FALSE, doEmbed = doEmbed, silent = silent)

  if(is.null(y2)){
    y2 <- y1
    attributes(y2) <- attributes(y1)
  }

  if(to.sparse){
    attributes(dmat)$directed <- directed
    attributes(dmat)$cumulative <- cumulative
    attributes(dmat)$includeDiagonal <- includeDiagonal
    #attributes(dmat)$weighted <- weighted
  } else {
    attr(dmat,"directed") <- directed
    attr(dmat,"cumulative") <- cumulative
    attr(dmat,"includeDiagonal") <- includeDiagonal
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

      grW <- igraph::graph_from_adjacency_matrix(dmat, weighted = TRUE, mode = mode, diag = TRUE) #includeDiagonal)
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
    g  <- igraph::graph_from_adjacency_matrix(dmat, weighted = weighted, mode = mode, diag = includeDiagonal)
    V(g)$y1 <- y1[1:igraph::vcount(g)]
    V(g)$y2 <- y2[1:igraph::vcount(g)]

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
                    Chromatic = NULL,
                    radiusValue = NA,
                    title = "", xlabel = "", ylabel="",
                    plotSurrogate = NA,
                    returnOnlyObject = FALSE){

  if(attributes(RN)$directed&attributes(RN)$cumulative&is.na(attributes(RN)$emRad)){
    warning("Plotting an unthresholded, directed,distance plot may not yield a sensible result!")
  }

  rp_plot(RM = RN,
          plotDimensions  = plotDimensions,
          plotMeasures    = plotMeasures,
          plotRadiusRRbar = FALSE,
          drawGrid = drawGrid,
          markEpochsLOI = markEpochsLOI,
          Chromatic = Chromatic,
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
#' @param g An igraph object.
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

  vertex_prop <- data.frame(LocalClustering=rep(NA,vcount(g)),
                            DegreeDensity=rep(NA,vcount(g)),
                            StrengthDensity=rep(NA,vcount(g)),
                            LocalCloseness=rep(NA,vcount(g)),
                            LocalEfficiency=rep(NA,vcount(g)),
                            LocalBetweenness=rep(NA,vcount(g))
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
      vertex_prop$LocalClustering  <- DirectedClustering::ClustBCG(igraph::as_adjacency_matrix(g, edges = TRUE, sparse = FALSE), "directed", isolates = "zero")$outCC
      degType <- "out"
      efType <- "local"
    } else {
      vertex_prop$LocalClustering  <- DirectedClustering::ClustBCG(igraph::as_adjacency_matrix(g, edges = TRUE, sparse = FALSE), "directed", isolates = "zero")$totalCC
    }
  }

  vertex_prop$LocalClustering  <- igraph::transitivity(g, type = clType, isolates = "zero")
  vertex_prop$DegreeDensity    <- igraph::degree(g, loops = FALSE)/(igraph::vcount(g)-1)
  vertex_prop$LocalCloseness   <- igraph::closeness(g, normalized = TRUE, mode = degType)
  vertex_prop$LocalEfficiency  <- brainGraph::efficiency(g,type = efType, use.parallel = FALSE, A = igraph::as_adjacency_matrix(g,edges = TRUE, sparse = FALSE))
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
  graph_prop$AveragePathLength   <- igraph::mean_distance(g,directed = FALSE, unconnected = FALSE)
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


