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
    gPL <- ggplot(out, aes(y = yStrength_log10, x = xDegree_log10)) +
      geom_point() +
      geom_line(aes(x = xDegree_log10, y = PowerLaw), colour = "steelblue3", size = 1) +
      scale_x_continuous("Vertex Degree (log10)") +
      scale_y_continuous("Vertex Strength (log10)") +
      theme_bw()
    print(gPL)
  }


  return(out)
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

  colvec <- heat.colors(n = length(unique(edgeFrame$rectime)), alpha = .8)

    #paletteer::paletteer_c(package = "scico",palette = "vik",n = length(unique(edgeFrame$rectime)))
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



#' Mutliplex Recurrence Network (DEPRECATED)
#'
#' @description This function will create a Multiplex Recurrence Network from a list of [igraph] objects that can be considered the layers of a network. The layers must have the same number of nodes. There are two modes of operation: *Layer similarity* (`MRNweightedBy` is set to `"InterLayerMI"` or `"EdgeOvelap"`) and *Layer importance* (`MRNweightedBy` is `"RankorderDC"`). The former generates weighted MRN based on Interlayer Mutual Information or Edge Overlap, the latter examines the relative importance of each layer by assigning a rank to each vertex (time point), based on a vertex measure passed in argument `MRNrankedBy`.
#'
#' @inheritParams make_spiral_graph
#' @inheritParams plotNET_groupColour
#' @param layers A list of igraph objects representing the layers of the multiplex network. The layer networks must all have the same number of vertices.
#' @param MRNweightedBy The measure to be used to evaluate the average structural similarities between the layers of the network. Valid options are: `"InterLayerMI"` (Mutual information based on similarity of the vertex degree across layers), `"EdgeOverlap"` (proportion of vertices sharing the same edges across layers), `"RankorderDC"` (Dynamic Complexity of the inter layer vertex rank order based on the vertex property/measure in `MRNrankedBy`). Choosing `"InterLayerMI"` or `"EdgeOverlap"` will decide which measure is displayed in the plot of the Multiplex RN, both measures will always be returned in the numerical output.
#' @param windowedWeights If a windowed analysis is conducted and the edges of the graphs in `layers` have a weight property, there are a number of different ways to handle the weights depending on the value of `windowedWeights`: `"none"`, `"local"`, and `"cumulative"`. Value `"none"` will ignore the weights, `"local"` will limit the range of edges to those edges connecting the vertices contained within the window, `"cumulative"` will consider all edges connecting to vertices in the current window, including edges from vertices connecting from previous windows. (default = `none`)
#' @param MRNrankedBy If `MRNweightedBy = "RankorderDC"`, then `MRNrankedBy` must be a valid [igraph] command that returns vertex properties, for example: `"degree"`, `"strength"`,`"hub_score"`,`"centr_degree"`, `"transitivity"`, `"betweenness"`. The appropriate measure type, e.g. for `"directed"`, or `"weighted"` graphs, will be inferred from the graph properties of the 1st graph object in the `layers` list. For best results with weighted measures, assign a value to `E(g)$weight` for each `g` in `layers`. (default = `"degree"`)
#' @param win The window size passed to [casnet::ts_windower()] in which to evaluate `"InterLayerMI"` or `"EdgeOvelap"`. If `MRNweightedBy = "RankorderDC"`, it will be the size of the right aligned window in which Dynamic Complexity will be computed using [casnet::dc_win()] (default = `NA`)
#' @param step The stepsize for the sliding window (default = `NA`)
#' @param overlap The window overlap passed to [casnet::ts_windower()] if `MRNweightedBy` is `"InterLayerMI"` or `"EdgeOvelap"`. The value of `step` will be ignored if `overlap` is not `NA`. (default = `NA`).
#' @param alignment Alignment of the window "l","c","r".
#' @param doPlot Plot the multiplex recurrence network (default = `TRUE`).
#' @param doSave Save the plots.
#' @param coords A data frame with layout coordinastes generated by calling any of the [igraph] layout functions. If `NA` a circle layout will; be generated (default = `NA`)
#' @param RNnodes Should the vertices represent the RN of the layers? This is recommended only for a small numbers of vertices. (default = `FALSE``)
#' @param vertexColour A vector of colours for the vertices. If this is a named list, names will be displayed in the legend.
#' @param showVertexLegend Show the vertex colour legend?
#' @param createAnimation If `createAnimation = TRUE` *and* `doPlot = TRUE` *and* a windowed analysis is conducted, an animation will be produced using either package `gganimate` (if `useImageMagick = FALSE`) or `animation` (if `useImageMagick = FALSE`). The main difference is that `gganimate` has nice animation transition features, but plots the MRN using [ggplot2], which does not have great options for displaying the nodes as images. With package `animation` a sequence of [igraph] plots will be converted to an animation. If `doSave = TRUE` the animation will be saved in `imageDir` as an animated gif by calling either [gganimate::anim_save()], or [animation::saveGIF()] (default = `FALSE`)
#' @param useImageMagick Should [ImageMagick](https://imagemagick.org/index.php) be used to create the animation. **NOTE:** ImageMagick has to be installed on your system, see [animation::saveGIF()] (default = `FALSE`)
#' @param loopAnimation Should the animation loop? (default = `TRUE`)`
#' @param transitionLength Length of each transition in the animation, ignored if `useImageMagick = TRUE` (default = `3`)
#' @param stateLength Value of `state_length` if `gganimate` is used, or the  `interval` in seconds for [animation::ani.pause()] (default = `1`)
#' @param gifWidth Width of the animated `gif` in pixels. The default width will be `600/150 = 4 in` or `10.16 cm` (default = `600`)
#' @param gifRes Resolution of the animated `gif` in `ppi` (default =`150`)
#' @param noParts Do not plot the individual graphs that make up the animation to the current `dev` (default = `TRUE`)
#' @param imageDir Directory to save the layer images and windowed MRN plots. If `NA`, the value returned by `getwd()` will be used, if `NULL` no windowed images will be saved (default = `NA`)
#' @param silent Silent-ish mode
#'
#' @return A matrix with edge weights between layers that represent the measure `MRNweightedBy`.
#'
#' @export
#'
#' @examples
#'
#'
rn_multiplex <- function(layers,
                         MRNweightedBy = c("InterlayerMI", "Edgeoverlap","RankDC")[1],
                         MRNrankedBy = "degree",
                         win = NA,
                         step = NA,
                         overlap = NA,
                         alignment = "l",
                         windowedWeights = c("none","local","cumulative")[1],
                         doPlot = FALSE,
                         doSave = FALSE,
                         coords = NA,
                         RNnodes = FALSE,
                         scaleVertexSize = c(.01,5),
                         vertexColour = getColours(length(layers)),
                         vertexBorderColour = "black",
                         showVertexLegend = TRUE,
                         showSizeLegend = FALSE,
                         alphaV = .7,
                         scaleEdgeSize = 1/5,
                         alphaE = .5,
                         showEdgeColourLegend = FALSE,
                         curvature = -0.6,
                         createAnimation = FALSE,
                         useImageMagick = FALSE,
                         loopAnimation = TRUE,
                         transitionLength = 3,
                         stateLength = 1,
                         gifWidth = 600,
                         gifRes = 150,
                         noParts = TRUE,
                         imageDir = NA,
                         silent = TRUE){

  cat("\n~~~o~~o~~casnet~~o~~o~~~\n")

  if(!all(plyr::laply(layers,function(g) igraph::is.igraph(g)))){
    stop("All elements of the layers list have to be an igraph object!")
  }

  if(!all(plyr::laply(layers, function(g) igraph::vcount(g))==igraph::vcount(layers[[1]]))){
    stop("In a Multiplex Recurrence Network, the layer networks must all have the same number of vertices!")
  }
  Nsize <- igraph::vcount(layers[[1]])

  if(is.null(names(layers))){
    names(layers) <- paste0("layer",1:length(layers))
  }

  doSave <- TRUE
  if(is.null(imageDir)){
    doSave <- FALSE
  }

  if(is.na(imageDir)){
    imageDir <- getwd()
  } else {
    if(!dir.exists(normalizePath(imageDir))){
      stop(paste0("The directory '",imageDir,"' does not exist."))
    }
  }

  if(is.na(win)){
    win <- Nsize
  }
  if(is.na(step)){
    step <- win
  }

  # if(win==Nsize|step==Nsize){
  #   wIndexList <- list(win)
  #   names(wIndexList) <- paste0("window: 1 | start: ",1," | stop: ",Nsize)
  # } else {
  wIndexList <- ts_windower(y = 1:Nsize, win = win, step = step, overlap = overlap, adjustY = NA, alignment = alignment)
  #}

  func <- "mi_interlayer"
  weighted <- igraph::is.weighted(layers[[1]])
  directed <- igraph::is.directed(layers[[1]])
  if(MRNweightedBy%in%"RankDC"){
    if(MRNrankedBy%in%ls(getNamespace("igraph"))){
      args <- formals(MRNrankedBy)
      addArgs <- ""
      if("directed" %in% names(args)){
        addArgs <- c(addArgs,",directed = directed")
      }
      if("weighted" %in% names(args)){
        addArgs <- c(addArgs,",weighted = weighted")
      }
      func <- paste0("igraph::",MRNrankedBy,"(layers[[l]]",addArgs,")")
    } else {
      stop("Value of MRNrankedBy is not a function of package igraph.")
    }
  }

  if(directed){
    combis <- DescTools::CombPairs(seq_along(layers),seq_along(layers))
    mode <- "directed"
  } else {
    combis <- DescTools::CombPairs(seq_along(layers))
    mode <- "upper"
  }

  if(func%in%"mi_interlayer"){

    cat(paste("\nWelcome to the multiplex... in layer similarity mode!\n\n"))

    out_mi <- out_eo <- out_means <- list()

    for(w in seq_along(wIndexList)){

      mp <- eo <- matrix(nrow=length(layers), ncol=length(layers), dimnames = list(names(layers),names(layers)))

      for(i in seq_along(combis$X1)){

        edgeFrame1 <- igraph::as_data_frame(layers[[combis$X1[i]]],"edges")
        edgeFrame1 <- edgeFrame1 %>% filter(edgeFrame1$from%[]%range(wIndexList[[w]])&edgeFrame1$to%[]%range(wIndexList[[w]]))
        ga  <- graph_from_data_frame(edgeFrame1, directed = directed)

        edgeFrame2 <- igraph::as_data_frame(layers[[combis$X2[i]]],"edges")
        edgeFrame2 <- edgeFrame2 %>% filter(edgeFrame2$from%[]%range(wIndexList[[w]])&edgeFrame2$to%[]%range(wIndexList[[w]]))
        gb  <- graph_from_data_frame(edgeFrame2, directed = directed)

        mp[combis$X1[i],combis$X2[i]] <- mi_interlayer(ga,gb)
        eo[combis$X1[i],combis$X2[i]] <- igraph::ecount(ga %s% gb) / (igraph::ecount(ga) + igraph::ecount(gb))
        #igraph::ecount(layers[[combis$X1[i]]] %s% layers[[combis$X2[i]]]) / (igraph::ecount(layers[[combis$X1[i]]]) + igraph::ecount(layers[[combis$X2[i]]]))
      }

      mi_ave   <- mean(mp[upper.tri(mp)&mp>0], na.rm = TRUE)
      mi_sd    <- stats::sd(mp[upper.tri(mp)&mp>0], na.rm = TRUE)
      eo_ave   <- mean(eo[upper.tri(eo)&eo>0], na.rm = TRUE)
      eo_sd    <- stats::sd(eo[upper.tri(eo)&eo>0], na.rm = TRUE)
      eo_sum   <- sum(plyr::laply(layers, function(g) length(igraph::E(g))),na.rm = TRUE)
      eo_all   <- eval(parse(text = paste0("igraph::ecount(intersection(",paste0("layers[[",1:NROW(layers),"]]",collapse=","),"))")))
      eo_joint <- eo_all / eo_sum

      out_mi[[w]] <- mp
      out_eo[[w]] <- eo
      out_means[[w]] <- data.frame(mi_mean = mi_ave, mi_sd = mi_sd,
                                   eo_mean = eo_ave, eo_sd = eo_sd,
                                   eo_sum  = eo_sum, eo_all = eo_all, eo_joint = eo_joint)

    } # wIndex

    names(out_mi)    <- names(wIndexList)
    names(out_eo)    <- names(wIndexList)
    names(out_means) <- names(wIndexList)

    if(doPlot){

      gList <- MRN <- list()

      if(MRNweightedBy%in%"InterLayerMI"){
        RN <- out_mi
      } else {
        RN <- out_eo
      }

      if(RNnodes){
        gSpiro  <- plyr::llply(layers, function(g) {
          gg <- make_spiral_graph(g                  = g,
                                  arcs               = 4,
                                  curvature          = curvature,
                                  a                  = .1,
                                  b                  = 1,
                                  alphaE             = alphaE,
                                  alphaV             = alphaV,
                                  scaleVertexSize    = scaleVertexSize,
                                  scaleEdgeSize      = scaleEdgeSize,
                                  showEpochLegend    = FALSE,
                                  epochColours       = getColours(Ncols = 20),
                                  showSizeLegend     = FALSE,
                                  defaultEdgeColour  = "grey80",
                                  vertexBorderColour = "gray99",
                                  edgeColourByEpoch  = TRUE,
                                  doPlot             = FALSE)
          gg <- gg + theme(
            panel.background      = element_rect(fill = "transparent"),
            plot.background       = element_rect(fill = "transparent", color = NA),
            panel.grid.major      = element_blank(),
            panel.grid.minor      = element_blank(),
            legend.background     = element_rect(fill = "transparent"),
            legend.box.background = element_rect(fill = "transparent")
          )
        })

        g_rast <- list()
        if(useImageMagick){
          bg <- "transparent"
        } else {
          bg <- "white"
        }
        for(f in seq_along(gSpiro)){
          ggplot2::ggsave(gSpiro[[f]],
                          filename = file.path(imageDir,paste0(names(gSpiro)[f],".png")),
                          device = "png", width = 100, height = 100, units = "mm", bg = bg, dpi = gifRes)
          g_rast[[f]] <- png::readPNG(file.path(imageDir,paste0(names(gSpiro)[f],".png")), native=TRUE) #magick::image_read(file.path(imageDir,paste0(names(gSpiro)[f],".png")))


        }

      } # RNnodes

      for(w in seq_along(RN)){

        mp_net <- graph_from_adjacency_matrix(adjmatrix = RN[[w]],
                                              mode      = mode,
                                              diag      = FALSE,
                                              weighted  = weighted)
        if(is.na(coord)){
          coord <- layout_in_circle(mp_net)
        } else {
          coord <- coord
        }

        if(RNnodes){
          V(mp_net)$raster  <- g_rast
          V(mp_net)$shape   <- "raster"
          V(mp_net)$size    <- 10
          V(mp_net)$size2   <- 10
        } else {
          V(mp_net)$size    <- 1
        }

        V(mp_net)$name    <- names(layers)
        E(mp_net)$width   <- elascer(E(mp_net)$weight,lo = 1,hi = 5)
        E(mp_net)$color   <- scales::gradient_n_pal(colours = c("steelblue","grey99","red3"), values = c(0,.5,1))(elascer(E(mp_net)$weight))
        V(mp_net)$label   <- V(mp_net)$name

        V(mp_net)$label.family <- "sans"
        V(mp_net)$label.font <- 2
        mp_net$layout <- coord


        if(useImageMagick){

          if(!noParts){
            plot(mp_net)
          }

          MRN[[w]] <- mp_net

        } else {

          gNodes        <- as.data.frame(mp_net$layout,what = "nodes")
          gNodes$ID     <- as.numeric(igraph::V(mp_net))
          gNodes$name   <- V(mp_net)$name
          gNodes$colour <- V(mp_net)$colour
          gNodes$alpha  <- V(mp_net)$alpha
          gNodes$size   <- V(mp_net)$size

          if(RNnodes){
            gNodes$image  <- paste0(file.path(imageDir,paste0(gNodes$name,".png")))
          }

          gEdges        <- igraph::get.data.frame(mp_net)
          if(any(is.character(gEdges$from))){
            cName <- "name"
          } else {
            cName <- "ID"
          }
          gEdges$from.x <- gNodes$V1[match(gEdges$from, gNodes[[cName]])]
          gEdges$from.y <- gNodes$V2[match(gEdges$from, gNodes[[cName]])]
          gEdges$to.x   <- gNodes$V1[match(gEdges$to, gNodes[[cName]])]
          gEdges$to.y   <- gNodes$V2[match(gEdges$to, gNodes[[cName]])]


          gEdges$width <- 1

          # Fix same coords
          # sameID <- which(gEdges$from.x==gEdges$to.x)
          # if(length(sameID)>0){
          #   gNodes$V2[gEdges$from[sameID]] <- gNodes$V2[gEdges$from[sameID]]+mean(diff(gNodes$V2))/2
          #   gEdges$to.x[sameID] <- gEdges$to.x[sameID]+mean(diff(gEdges$to.x))/2
          #   gEdges$to.y[sameID] <- gEdges$to.y[sameID]+mean(diff(gEdges$to.y))/2
          # }

          gg <- ggplot(gNodes,aes(x=V1,y=V2)) +
            geom_curve(data=gEdges, aes(x = from.x,
                                        xend = to.x,
                                        y = from.y,
                                        yend = to.y,
                                        colour = weight),
                       size = gEdges$width*scaleEdgeSize,
                       alpha = alphaE,
                       curvature = 0)

          if(RNnodes){
            gg <- gg + ggimage::geom_image(aes(image=image), size = .2)
          } else {
            gg <- gg + geom_point(aes(fill = vertexColour, size = size), pch=21, colour = vertexBorderColour, alpha = alphaV)
          }

          gg <- gg + scale_size(range = scaleVertexSize) +
            scale_colour_gradient2(low      = "steelblue",
                                   high     = "red3",
                                   mid      = "white",
                                   na.value = scales::muted("slategray4"),
                                   midpoint = median(gEdges$weight)) +
            labs(title = names(RN)[w]) +
            coord_cartesian(clip="off") +
            theme_void() + theme(plot.margin = margin(50,50,50,50, unit = "pt"),
                                 legend.margin = margin(l = 20, unit = "pt"),
                                 plot.title = element_text(margin = margin(b=20)))

          if(showVertexLegend){
            gg <- gg + guides(fill = guide_legend(title.position = "top",
                                                  byrow = TRUE,
                                                  nrow=2,
                                                  override.aes = list(size=5, order = 0)))
          } else {
            gg <- gg + guides(fill = "none")
          }

          if(showSizeLegend){
            gg <- gg + guides(size = guide_legend(title.position = "top",
                                                  byrow = TRUE,
                                                  nrow=2,
                                                  override.aes = list(legend.key.size = unit(1.2,"lines"), order = 1)))
          } else {
            gg <- gg + guides(size = "none")
          }

          if(showEdgeColourLegend){
            gg <- gg + guides(colour = guide_legend(title.position = "top",
                                                    byrow = TRUE,
                                                    nrow=2,
                                                    override.aes = list(size=5, order = 3)))
          } else {
            gg <- gg + guides(colour = "none")
          }

          if(!noParts){
            plot(gg)
          }

          MRN[[w]] <- gg
          gList[[w]] <- list(gNodes = gNodes,
                             gEdges = gEdges)

        } # useImageMagick

      } # w

      names(MRN)   <- names(wIndexList)


      if(!createAnimation){
        return(list(MRN           = MRN,
                    interlayerMI  = out_mi,
                    edgeOverlap   = out_eo,
                    meanValues    = out_means))
      } else {

        if(useImageMagick){

          animation::ani.options(interval = stateLength, imgdir = imageDir, loop = loopAnimation, nmax = length(MRN),
                                 ani.width = gifWidth, ani.res = gifRes)

          if(doSave){

            animation::saveGIF(
              for(i in seq_along(MRN)){
                plot(MRN[[i]], xlab = names(MRN)[i])
                animation::ani.pause()
              }, img.name = paste0(names(MRN)[i]), movie.name = paste0("MRN_animation_win",win,"_step",step,".gif")
            )

          } else {
            for(i in seq_along(MRN)){
              plot(MRN[[i]], xlab = names(MRN)[i])
              animation::ani.pause()
            }
          }

          return(list(MRN           = MRN,
                      interlayerMI  = out_mi,
                      edgeOverlap   = out_eo,
                      meanValues    = out_means))

        } else {

          names(gList) <- names(wIndexList)

          gNodes <- plyr::ldply(gList , function(g) g$gNodes)
          gEdges <- plyr::ldply(gList , function(g) g$gEdges)

          gg <- ggplot(gNodes,aes(x=V1,y=V2)) +
            geom_curve(data=gEdges, aes(x = from.x,
                                        xend = to.x,
                                        y = from.y,
                                        yend = to.y,
                                        size = weight,
                                        colour = weight), curvature = 0)
          if(RNnodes){
            gg <- gg + ggimage::geom_image(aes(image=image), size = .2)
          }
          gg <- gg + scale_size(range = c(.01,5)) +
            scale_colour_gradient2(low      = "steelblue",
                                   high     = "red3",
                                   mid      = "grey99",
                                   na.value = scales::muted("slategray4"),
                                   midpoint = mean(gEdges$weight)) +
            labs(caption = "{closest_state}") +
            gganimate::transition_states(factor(.id),
                                         transition_length = transitionLength,
                                         state_length = stateLength, wrap = loopAnimation) +
            gganimate::enter_fade() +
            gganimate::exit_fade() +
            coord_cartesian(clip="off") +
            theme_void() + theme(plot.margin = margin(50,50,50,50, unit = "pt"),
                                 legend.margin = margin(l = 20, unit = "pt"),
                                 plot.title = element_text(margin = margin(b=20)))

          plot(gg)

          if(doSave){
            #if(file.exists())
            gganimate::anim_save(filename =  paste0("MRN_animation_win",win,"_step",step,".gif"),
                                 animation = gg, path = file.path(imageDir))
          }

          return(list(MRN              = MRN,
                      MRNanimationData = gList,
                      MRNanimationGG   = invisible(gg),
                      interlayerMI  = out_mi,
                      edgeOverlap   = out_eo,
                      meanValues    = out_means))
        }
      }
    } else { #doPlot
      return(list(interlayerMI  = out_mi,
                  edgeOverlap   = out_eo,
                  meanValues    = out_means))
    }

  } else {    # NOT mi_interlayer

    cat(paste("\nWelcome to the multiplex... in layer importance mode!\n\n"))

    mr <- mr_rank <- matrix(ncol = length(layers),
                            nrow = Nsize,
                            dimnames = list(NULL, names(layers)))

    for(l in seq_along(layers)){
      mr[,l] <- eval(parse(text = func))
    }
    mr <- data.frame(mr)

    for(i in 1:NROW(mr_rank)){
      mr_rank[i,] <- as.numeric(rank(mr[i,]))
    }

    mr_rank           <- data.frame(mr_rank)
    colnames(mr_rank) <- paste0(colnames(mr),"_rank")

    mr_rankDC         <- dc_win(df = mr_rank,
                                win = win,
                                scale_min = min(mr_rank, na.rm = TRUE),
                                scale_max = max(mr_rank, na.rm = TRUE))
    mr_rankDC$meanDC    <- rowMeans(mr_rankDC)
    mr_rankDC$sdDC      <- as.numeric(colwise(sd, na.rm = TRUE)(data.frame(t(mr_rankDC))))
    colnames(mr_rankDC) <- c(paste0(colnames(mr),"_DC"),"mean_DC","sd_DC")

    mr$sum       <- rowSums(mr[,1:length(layers)],na.rm = TRUE)
    mr$mean      <- rowMeans(mr[,1:length(layers)], na.rm = TRUE)
    mr$sd        <- as.numeric(colwise(sd, na.rm = TRUE)(data.frame(t(mr[,1:length(layers)]))))
    colnames(mr) <- paste0(c(paste0(names(layers),"_"),"sum_","mean_","sd_"),MRNrankedBy)

    out_rank <- cbind(time = 1:NROW(mr), mr, mr_rank, mr_rankDC)
    attr(out_rank,"MRNrankedBy") <- MRNrankedBy

    if(doPlot){

      mr_rank$Time <- 1:NROW(mr_rank)
      df_ranks <- tidyr::gather(mr_rank, key = "Layer", value = "Rank", -c("Time"))
      df_ranks$Layer <- gsub("_rank","",df_ranks$Layer)

      g1 <- ggplot2::ggplot(df_ranks, aes(x = Time, y = Rank)) +
        geom_step(aes(colour = Layer)) +
        facet_grid(Layer~.) +
        scale_y_continuous(name = paste0(MRNrankedBy," (rank)"), expand = c(0,0)) +
        scale_x_continuous(name = "Time", expand = c(0,0)) +
        scale_colour_manual("Layer", values = getColours()) +
        guides(colour = "none") +
        theme_bw()

      g2 <- plotDC_res(mr_rankDC[1:length(layers)],
                       win=win,
                       subtitle = paste0("Dynamic Complexity of rank order of ",MRNrankedBy),
                       doPlot = FALSE)

      tmp <- data.frame(Time = 1:NROW(mr_rankDC))
      tmp$mean_DC  <- mr_rankDC$mean_DC
      tmp$upper_CI <- mr_rankDC$mean_DC+(mr_rankDC$sd_DC/sqrt(length(layers)))*1.96
      tmp$lower_CI <- mr_rankDC$mean_DC-(mr_rankDC$sd_DC/sqrt(length(layers)))*1.96
      df_DC <- tidyr::gather(tmp, key = "Mean", value = "DC", -c("Time"))
      df_DC$Mean <- gsub("_"," ",df_DC$Mean)
      cols <- getColours()[c(3,6,4)]
      names(cols) <- c("upper CI","mean DC","lower CI")

      g3 <- ggplot2::ggplot(df_DC, aes(x=Time,y=DC)) +
        geom_line(aes(colour = Mean)) +
        scale_y_continuous(name = paste0("Mean DC of ",MRNrankedBy," (rank)"), expand = c(0,0)) +
        scale_x_continuous(name = "Time", expand = c(0,0), limits = c(win,Nsize)) +
        scale_colour_manual(name = "", values = cols) +
        theme_bw()

      top_row <- cowplot::plot_grid(g1,g2,greedy = FALSE,nrow = 1)
      gs      <- cowplot::plot_grid(top_row,g3, ncol = 1)

      return(list(plot = invisible(gs), rankData = out_rank))

    } else {

      return(out_rank)

    } # doPlot

  } # rankDC
  cat("\n\n~~~o~~o~~casnet~~o~~o~~~\n")
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
              edge_meassures  = edge_prop,
              graph_measures = graph_prop)
  )
}


