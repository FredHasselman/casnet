# package casnet ----
#
# Plotting ----
#
# plot functions


#' Plot Multivariate Time Series Data
#'
#' @param df A data frame with time series in columns.
#' @inheritParams plotDC_res
#' @param timeVec If numeric, the number of the column in `df` which contains a time keeping variable. If `NA`, the time vector will be `1:NROW(df)` (default = `NA`)
#' @param groupVec A vector indicating the names of the time series in the columns of `df`. If `NA`, the column names of `df` will be used, excluding the `timeVec`, if present (default = `NA`)
#' @param returnPlotData Return the restructured data frame used to create the plot (default = `FALSE`)
#' @param useRibbon Neat for distributions
#' @param overlap Multiplier for scaling the series around the y-offset. Default is `offset + elascer(y, lo = -.45*overlap, hi = .45*overlap)` and if `useRibbon = TRUE` it is `offset + elascer(y, lo = 0*overlap, hi = .95*overlap)`. (default = `1`)
#'
#' @return A [ggplot] object.
#'
#' @export
#'
#' @examples
#'
#' # Use the coloured noise data set.
#' data(ColouredNoise)
#' plotTS_multi(ColouredNoise)
#'
plotTS_multi <- function(df,timeVec = NA, groupVec = NA, useVarNames = TRUE, colOrder = TRUE, doPlot = TRUE, title = "", subtitle = "", xlabel = "Time", ylabel = "", returnPlotData = FALSE, useRibbon = FALSE, overlap = 1){

  cname <- "time"
  i <-0
  while(cname%in%colnames(df)){
    i%++%1
    cname <- paste0(cname,i)
  }
  if(is.na(timeVec)){
    df[[cname]] <- 1:NROW(df)
  } else {
    if(is.numeric(timeVec)&length(timeVec)==1){
      colnames(df)[timeVec] <- cname
    } else {
      stop("Argument timeVec is not interpretable as a column number.")
    }
  }

  time <- df[[cname]]
  df <- df %>% dplyr::select(where(is.numeric))
  eval(parse(text=paste("df$",cname," <- time")))


  if(is.na(groupVec)){
    if(colOrder){
      groupVec <- colnames(df)[-c(cname%ci%df)]
    } else {
      groupVec <- sort(colnames(df)[-c(cname%ci%df)])
    }
  } else {
    if(length(groupVec)!=(NCOL(df)-1)){
      stop("Length of groupVec does not match number of time series in df.")
    }
  }

  tmp <- tidyr::gather(df, key = "timeSeries", value = "y", -(.data[[cname]]), factor_key = colOrder)
  #tmp$timeSeries <- ordered(tmp$timeSeries)

  colnames(tmp)[cname%ci%tmp] <- "time"

  yOrder <- groupVec
  names(yOrder) <- paste(groupVec)
  offsets     <- names(yOrder)
  offsets <- {stats::setNames(0:(length(offsets) - 1), offsets)}
  tmp$offsets <- unlist(plyr::llply(seq_along(yOrder), function(n) rep(offsets[n], sum(tmp$timeSeries%in%names(offsets)[n]))))

  ymin <- -.45 * overlap
  ymax <-  .45 * overlap

  if(useRibbon){
    ymin <- 0
    ymax <- .95 * overlap
  }

  xOrder <- as.numeric_discrete(time)

  tmp$time    <- as.numeric(xOrder)
  xBreaks <- hist(xOrder, plot = FALSE)$mids
  xLabels <- names(xOrder)[xBreaks]

  if(max(nchar(xLabels), na.rm = TRUE)>4){
    ang <- 45
    hj  <- 1
  } else {
    ang <- 0
    hj  <- 0.5
  }

  # Calculate and scale group densities
  pdat <- tmp %>%
    dplyr::group_by(.data$timeSeries) %>%
    dplyr::mutate(y = elascer(.data$y,lo = ymin, hi = ymax)) %>%
    dplyr::ungroup()
  pdat$y_offset <- pdat$y + pdat$offsets
  pdat$ymin <- pdat$offsets
  pdat$ycol <- factor(pdat$offsets%%2)
  pdat$timeSeries <- rev(tmp$timeSeries)



  g <- ggplot2::ggplot(pdat, ggplot2::aes_(x=~time, y = ~y_offset, group = ~timeSeries))

  if(useRibbon){
    g <- g + geom_ribbon(aes_(ymin = ~ymin, ymax = ~y_offset, fill = ~ycol, colour = ~ycol)) +
      scale_fill_manual(values = c("0"="grey30","1"="grey70"), guide = "none") +
      scale_color_manual(values = c("0"="grey90","1"="grey90"), guide = "none")
  } else {
    g <- g +  geom_path()
  }

  g <- g +
    scale_y_continuous(ylabel, breaks = offsets, labels = names(offsets), expand = c(0,0)) +
    scale_x_continuous(xlabel, breaks = xBreaks, labels = xLabels, expand = c(0,0)) +
    ggtitle(label = title, subtitle = subtitle) +
    theme_minimal()

  if(useRibbon){
    g <- g +  theme(axis.text.y = element_text(size = 8, face = "plain"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.text.x = element_text(angle = ang, vjust = 1, hjust = hj),
                    plot.margin = margin(0,1,1,1,"line"))
  } else{
    g <- g +   theme(axis.text.y = element_text(size = 8, face = "plain"),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     axis.text.x = element_text(angle = ang,  vjust = 1, hjust = hj),
                     plot.margin = margin(0,1,1,1,"line"))
  }

  if(doPlot){
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
#' @param labels Vertex labels. If `NA` is passed then first `V(g)$name` will be checked, for node labels. To create a plot without labels pass `NULL` (default = `NA`)
#' @param nodeSize Set node sizes by `degree(g, normalised = TRUE)` (default), `hubscore(g)$vector`, or, `strength(g)`, `eccentricity(g)`, `coreness(g)`. If a numeric value is passed all vertex sizes will be set to that value.
#' @param rescaleSize Use to rescale the measure indicated under `nodeSize` to `c(min,max)` for better node visibility. Pass `c(1,1)` for no rescaling (default = `c(3,10)`)
#' @param labelSize Set labelsize: `"asnodesize"` sets the `cex` for the labels to coincide with nodesize (with min of .4 and max of 1.1). A single numeric value sets the `cex` of all labels to that value. A numeric vector of length two, `c(min,max)` will scale the label sizes to `min` and `max`
#' @param nodeColour Set to `TRUE` to colour nodes using the size of the node (default = `TRUE`)
#' @param edgeWeight Set size of edges to `"E(g)$weight"` by passing "weight". If a single numeric value is provided all edges will be set to that value.
#' @param edgeColour Set to `TRUE` to colour edges using the weight of the edge (default = `FALSE`)
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
                         nodeSize   = c("degree","hubscore","strength","eccentricity","coreness")[1],
                         rescaleSize = c(3,10),
                         nodeColour = TRUE,
                         labelSize  = "asnodesize",
                         edgeWeight = "weight",
                         edgeColour = FALSE,
                         removeZeroDegree = TRUE,
                         removeSelfLoops  = TRUE,
                         doPlot     = TRUE){

  checkPkg("igraph")

  if(removeSelfLoops){
    g <- simplify(g)
  }
  if(removeZeroDegree){
    g <- delete.vertices(g,degree(g)==0)
  }

  rev <- NA
  if(is.character(nodeSize)){
    switch(nodeSize,
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
    rev <- rep(as.numeric(nodeSize),length(igraph::V(g)))
  }

  # set colors and sizes for vertices
  #rev<-elascer(log1p(igraph::V(g)$degree))

  igraph::V(g)$size        <- elascer(rev, lo = rescaleSize[1], hi = rescaleSize[2])

  rev <- rev/max(rev, na.rm = TRUE)
  rev[rev<=0.2]<-0.2
  rev[rev>=0.9]<-0.9
  igraph::V(g)$rev <- rev

  gradFunc <- getColours(Ncols = c(2,7,1), continuous = TRUE, Dcols = c(0,.5,1))
  igraph::V(g)$color <- gradFunc(elascer(igraph::V(g)$size))

  # set vertex labels and their colors and sizes
  if(all(is.na(labels))){
    if(!is.null(V(g)$name)){
      igraph::V(g)$label <- igraph::V(g)$name
    } else {
        V(g)$label <- ""
      }
    } else {
      if(is.null(labels)){
        V(g)$label <- ""
      } else {
       if(is.character(labels)&(length(labels)==igraph::vcount(g))){
         V(g)$label <- labels
       } else
         warning("Length of labels does not equal number of nodes.")
      }
    }

    if(labelSize == "asnodesize"){igraph::V(g)$label.cex <- elascer(igraph::V(g)$size,lo = .4, hi = 1.1)}
    if(is.numeric(labelSize)&length(labelSize)==1){igraph::V(g)$label.cex <- labelSize}
    if(is.numeric(labelSize)&length(labelSize)==2){igraph::V(g)$label.cex <- elascer(igraph::V(g)$size, lo = labelSize[1], hi = labelSize[2])}

    igraph::V(g)$label.color <- "black"
    igraph::V(g)$label.family = "Helvetica"

  if(igraph::ecount(g)>0){
    if(edgeWeight%in%"weight"){
      if(igraph::is_weighted(g)){
        igraph::E(g)$width <- elascer(igraph::E(g)$weight,lo = .8, hi = 5)
      }
    } else {
      if(is.numeric(edgeWeight)){
        igraph::E(g)$width <- as.numeric(edgeWeight)
      } else {
        igraph::E(g)$width <- 1
      }
    }
  }

    if(edgeColour){
      gradFunc <- getColours(Ncols = c(21,5,1), continuous = TRUE, Dcols = c(0,.5,1))
      igraph::E(g)$color <- gradFunc(elascer(igraph::E(g)$width))
    } else {
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


#' Vertex and Edge Group Colours
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
    if(length(unigroups)<=58){
      groupColours <-  getColours(Ncols = length(unigroups))
    } else {
      groupColours <- getColours(Ncols = 1:58, continuous = TRUE, Dcols = seq(0, 1, length.out = 58))(seq(0,1,length.out=length(unigroups)))
        # scales::gradient_n_pal(scales::brewer_pal(palette="Paired")(12))(seq(0, 1, length.out = length(unigroups)))
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
    if(length(alphaV)==1|length(alphaV)==igraph::vcount(g)){
      igraph::V(g)$alpha <- alphaV
    } else {
      stop("Length of vector alphaV is not equal to 1 or vcount(g)")
    }
  } else {
    stop("All alphaV must be in [0,1]")
  }

  if(colourE){
    #init
    igraph::E(g)$color  <- defaultEdgeColour
  }

  if(length(igraph::E(g)$color)>0){
    igraph::E(g)$colour <- igraph::E(g)$color
  }


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
        }
        else {
          #  # Add a default colour and alphac
          igraph::E(g)[Eid]$color  <- add_alpha(igraph::E(g)[Eid]$color, alpha = igraph::E(g)[Eid]$alpha)
          #  igraph::E(g)[Eid]$colour <- igraph::E(g)[Eid]$color
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
#' @param doPlot Display the plot (A plot object is always returned invisibly)
#'
#' @return A ggplot object
#'
#' @export
#'
plotFD_loglog <- function(fd.OUT, title="", subtitle="", xlabel="Bin size", ylabel="Fluctuation", logBase=NA, doPlot = TRUE){

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
  if(doPlot){
  graphics::plot(g)
    }

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
plotRED_mif <- function(mif.OUT = NULL,
                        lags  = attr(mif.OUT,"lags"),
                        nbins = attr(mif.OUT,"nbins"),
                        surTest = FALSE,
                        alpha=.05,
                        doPlot = TRUE,
                        returnMIFun = TRUE){

  if(is.null(mif.OUT)){
    stop("No data.")
  } else {
    if(is.null(attr(mif.OUT,"miType"))){
      stop("Argument mif.OUT is not the output of function mif().")
    }
  }


  #mifunMIF  <- mif(y = y, lags = lags, nbins = nbins)
  # #ldply(seq_along(df.acf$acf), function(cc){pacf_fisherZ(r=df.acf$acf[cc],n=dfN[cc],lag=df.acf$lag[cc],type="acf")})
  #mifunPMIF <- mif(y = cbind(y,y[,1],y[,1]), lags = lags, nbins = nbins)
  # #ldply(seq_along(df.pacf$acf), function(cc){pacf_fisherZ(r=df.pacf$acf[cc],n=dfN[cc],lag=df.pacf$lag[cc],type="pacf")})
  # mif.OUT   <- rbind(mifunMIF,mifunPMIF)

  groupColours <-  scales::brewer_pal(palette="RdBu")(11)
  cols <- c("yes"=groupColours[9],"no"=groupColours[3])

  mifun_long <- data.frame(lag =  c(as.numeric(names(mif.OUT))), #,as.numeric(names(mifunPMIF))),
                           mi = c(mif.OUT),
                           type = c(rep(attributes(mif.OUT)$miType,NROW(mif.OUT)))) #,rep(attributes(mifunPMIF)$miType,NROW(mifunPMIF))))

  g <- ggplot2::ggplot(mifun_long,ggplot2::aes_(x=~lag,y=~mi)) +
    ggplot2::geom_hline(yintercept = 0, colour="grey",size=1) +
    #  ggplot2::geom_line(data = data.frame(x=c(0,mifun_long$lag[1]),y=c(1,mifun_long$mi[1])),ggplot2::aes_(x=~x,y=~y),colour="grey50")

    # if(length(lags)<=50){
    #  g <- g + ggplot2::geom_point(x=0,y=1,colour=groupColours[10],fill=groupColours[9],size=2,pch=21)
    # }
    #geom_ribbon(aes_(ymin=~ciL,ymax=~ciU),fill="grey70",colour="grey50") +
    ggplot2::geom_path(colour="grey50") +
    ggplot2::geom_point(pch=21, cex=(1 + .01*(NROW(mifun_long$mi)/nbins))) +
    ggplot2::facet_grid(type ~.) +
    # scale_fill_manual(bquote(p < .(siglevel)),values = cols,
    #                   labels =  list("yes"= expression(rho != 0),
    #                                  "no" = expression(rho == 0))) +
    # scale_colour_manual(bquote(p < .(siglevel)),values = cols,
    #                     labels =  list("yes"= expression(rho != 0),
    #                                    "no" = expression(rho == 0))) +
    ggplot2::scale_x_continuous(limits = c(min(lags),max(lags)),expand = c(0.01,0), breaks = seq(min(lags),max(lags), by = round(length(lags)/10))) +
    ggplot2::scale_y_continuous(limits = c(min(mifun_long$mi,na.rm = TRUE),max(mifun_long$mi,na.rm = TRUE))) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank())

  if(doPlot){
    # graphics::plot.new()
    graphics::plot(g)
  }

  if(returnMIFun){
    return(list(mifun=mif.OUT,
                plot=invisible(g)))
  } else {
    return(invisible(g))
  }
}



#' Plot Complexity Resonance Diagram
#'
#' @inheritParams dc_ccp
#' @inheritParams dc_win
#' @param resVariable Variable displayed in the plot.
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
plotDC_res <-  function(df_win, win, useVarNames = TRUE, colOrder = TRUE, useTimeVector = NA, timeStamp = "01-01-1999", doPlot = TRUE, title = 'Complexity Resonance Diagram', resVariable = "Dynamic Complexity", subtitle = "", xlabel = "Time", ylabel = ""){

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
    labels <- labels[c(seq(2,length(minorBreaks), by = round(length(minorBreaks)/25)))]
    majorBreaks <- minorBreaks[c(seq(2,length(minorBreaks), by = round(length(minorBreaks)/25)))]
  } else {
        #minorBreaks <- NULL
        majorBreaks <- NULL
  }

  #
  #   if(NROW(df_win)>50){
  #     labels <- labels[minorBreaks]
  #     #  by = round(length(breaks)/25))]
  #     #labels[!labels%in%labels[seq(2,length(minorBreaks), by = round(length(minorBreaks)/25))]] <- ""
  #     labels <- labels[c(seq(2,length(minorBreaks), by = round(length(minorBreaks)/25)))]
  #     majorBreaks <- minorBreaks[c(seq(2,length(minorBreaks), by = round(length(minorBreaks)/25)))]
  #   } else {
  #     minorBreaks <- NULL
  #     majorBreaks <- NULL
  #   }

  # breaks <- seq(win, NROW(df_win))
  # if(NROW(df_win)>50){
  #   labels <- paste(breaks)
  #   labels[!labels%in%labels[seq(2,length(breaks), by = round(length(breaks)/25))]]<- ""
  #   labels[1] <- ""
  #   labels[2] <- paste(win+1)
  # }

  if(is.na(colOrder)&subtitle==""){
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
    ggplot2::scale_fill_gradient2(resVariable,low='steelblue', high='red3', mid='whitesmoke', midpoint=(max(dfp$value, na.rm=TRUE)/2), na.value='white') +

    #   if(is.null(majorBreaks)){
    #   g <- g +
    #       ggplot2::scale_y_discrete(ylabel, expand = c(0,0)) +
    #       ggplot2::scale_x_continuous(xlabel, expand = c(0,0), limits = c((win+1)-.5, NROW(dfp)+.5))
    #
    #   } else {
    # g <- g +
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
#' @inheritParams dc_win
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
  dfp <- tidyr::gather(df_ccp, key = "variable", value = "value", -c(.data$time), factor_key = TRUE)
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
#' @param df_ccp If an object generated by [dc_ccp()], the levels shown in the plot will only be displayed if there is an cumulative complexity peak at that time point (default = `NA` )
#' @inheritParams plotDC_ccp
#' @inheritParams plotDC_res
#' @inheritParams dc_win
#' @param levelName A name for the state variable.
#'
#' @return An invisible ggplot2 object.
#' @export
#'
#' @family Dynamic Complexity functions
#'
plotDC_lvl <-  function(df_win, df_ccp = NA, df_lvl, win, useVarNames = TRUE, colOrder = TRUE, useTimeVector = NA, timeStamp = "31-01-1999", doPlot = TRUE, title = 'Peaks versus Levels Plot', subtitle = "", xlabel = "Time", ylabel = "", levelName = "State variable"){


  if(!is.na(dc_ccp)){

    if(!(all(NROW(df_win)==NROW(df_ccp) & NROW(df_win)==NROW(df_lvl$pred) & NROW(df_ccp)==NROW(df_lvl$pred)))){
      stop("The time series must all have equal lengths:\n all(NROW(df_win)==NROW(df_ccp) & NROW(df_win)==NROW(df_lvl$pred) & NROW(df_ccp)==NROW(df_lvl$pred))")
    }

    if(!(NCOL(df_win)==(NCOL(df_ccp)-1))){
      stop("There must be an equal number of time series variables in 'df_win' and 'df_ccp'")
    }

  }


  ylabel <- paste0(ylabel," (arbitrary units)")

  if(!useVarNames){
    if(!is.na(dc_ccp)){
      df_ccp <- data.matrix(df_ccp)
      colnames(df_ccp) <- c(1:ncol(df_ccp))
    }

    df_win <- data.matrix(df_win)
    colnames(df_win) <- c(1:ncol(df_win))
  }

  if(!is.na(dc_ccp)){
    lcol  <- df_ccp$sig.peaks
    lcol[lcol==0] <- 5
    df_ccp <- df_ccp[,-c("sig.peaks"%ci%df_ccp)]
    if(is.na(colOrder)){
      df_ccp <- dplyr::select(df_ccp,names(sort(colSums(df_ccp, na.rm = TRUE))))
      colOrder <- TRUE
    }
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


#' Plot windowed Multiplex Recurrence Network measures
#'
#' @param df_mrn Output from function [mrn()] with arguments set for a windowed analysis.
#' @param measures Character vector indicating which measures should be plotted. Valid elements in the vector are: `"mi"` for inter-layer mutual information,`"eo"` for edge overlap, or any function name from package [igraph] that can be applied to extract vertex-based measures from a weighted network. This is the multiplex recurrence network with a number of nodes equal to the number of layers and edge weights set by argument `mrnWeights`. (default = `"mi"`)
#' @param mrnWeights Which measure should be used for the MRN edge-weights? Valid options are `"mi"` for inter-layer mutual information,`"eo"` for edge overlap. (default = `"mi"`)
#' @param showSeries Should the time series that constitute the layers be plotted below the windowed MRN measures? This will only work if the [mrn()] was called using graphs that were generated by function [rn()] with `returnGraph = TRUE`, in which case they will have a property `V(g)$y1` (and if applicable `V(g)$y2`), which will be plotted. If more than 1 measure is requested in `measures` and `showSeries = TRUE` a seperate plot for each measure will be produced. (default = `FALSE`)
#' @param doPlot Plot the igraph object.
#'
#' @return a ggplot object
#' @export
#'
#' @family Tools for plotting multiplex recurrence networks
#' @family Tools for windowed analyses
#'
plotMRN_win <- function(df_mrn,
                        measures = "mi",
                        mrnWeights = "mi",
                        showSeries = FALSE,
                        doPlot     = TRUE){


  if(is.list(df_mrn)&(length(df_mrn)==3)&("meanValues"%in%names(df_mrn))){

    df   <- plyr::ldply(df_mrn$meanValues)
    df$t <- attr(df_mrn,"time")
    df_long <- df %>% dplyr::select(dplyr::one_of(c("t","mi_mean","eo_mean"))) %>% tidyr::pivot_longer(-t,names_to = "measure",values_to = "y")



    g <-  ggplot(df_long, aes_(x= ~t, y = ~y, colour = ~measure)) +
      geom_line() +
      facet_grid(measure~., scales="free_y") +
      scale_x_continuous(limits = c(0,max(df$t)),expand = c(0.1,0.1)) +
      theme_bw()

    if(doPlot){
      # graphics::plot.new()
      graphics::plot(g)
    }

    return(invisible(g))

  }


  #
  # rev <- NA
  # if(is.character(nodeSize)){
  #   switch(nodeSize,
  #          # degree   = rev <- elascer(log1p(igraph::degree(g,normalized = TRUE))),
  #          # hubscore = rev <- elascer(log1p(igraph::hub_score(g)$vector))
  #          # strength = rev <- elascer(log1p(igraph::strength(g)$vector)))
  #          degree   = rev <- igraph::degree(g),
  #          hubscore = rev <- igraph::hub_score(g)$vector,
  #          strength = rev <- igraph::strength(g),
  #          eccentricity = rev <- igraph::eccentricity(g),
  #          coreness = rev <- igraph::coreness(g),
  #   )
  # } else {
  #   rev <- rep(as.numeric(nodeSize),length(igraph::V(g)))
  # }
  #
  # # set colors and sizes for vertices
  # #rev<-elascer(log1p(igraph::V(g)$degree))
  #
  # igraph::V(g)$size        <- rev
  #
  # rev <- rev/max(rev, na.rm = TRUE)
  # rev[rev<=0.2]<-0.2
  # rev[rev>=0.9]<-0.9
  # igraph::V(g)$rev <- rev
  #
  # igraph::V(g)$color       <- grDevices::rgb(igraph::V(g)$rev, 1-igraph::V(g)$rev,  0, 1)
  #
  # # set vertex labels and their colors and sizes
  # if(all(is.na(labels))){
  #   igraph::V(g)$label       <- ""
  # } else {
  #   igraph::V(g)$label       <- labels
  #
  #   if(labelSize == "asnodesize"){igraph::V(g)$label.cex <- elascer(igraph::V(g)$size,lo = .4, hi = 1.1)}
  #   if(is.numeric(labelSize)&length(labelSize)==1){igraph::V(g)$label.cex <- labelSize}
  #   if(is.numeric(labelSize)&length(labelSize)==2){igraph::V(g)$label.cex <-  elascer(igraph::V(g)$size, lo = labelSize[1], hi = labelSize[2])}
  #
  #   igraph::V(g)$label.color <- "black"
  #
  #   igraph::V(g)$label.family = "Helvetica"
  # }
  #
  # if(igraph::ecount(g)>0){
  #   if(edgeweight%in%"weight"){
  #     if(igraph::is_weighted(g)){
  #       igraph::E(g)$width <- elascer(igraph::E(g)$weight,lo = .8, hi = 5)
  #     }
  #   } else {
  #     if(is.numeric(edgeweight)){
  #       igraph::E(g)$width <- as.numeric(edgeweight)
  #     } else {
  #       igraph::E(g)$width <- 1
  #     }
  #   }
  #   igraph::E(g)$color <- grDevices::rgb(0.5, 0.5, 0.5, 1)
  # }

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


# Themes, layouts ----

#' Layout a graph on a spiral
#'
#' @param g An igraph object. If (`rev = FALSE`) the vertex with the lowest index will be placed in the centre of the spiral, the highest index will be most outer vertex,
#' @param type Spiral type, one of `"Archimedean"`,`"Bernoulli"`,`"Fermat"`, or, `"Euler"` (default = `"Archimedean"`)
#' @param arcs The number of arcs (half circles/ovals) that make up the spiral (default = `10`)
#' @param a Parameter controlling the distance between spiral arms, however, the effect will vary for different spiral types (default = `0.5`)
#' @param b Parameter controlling where the spiral originates. A value of 1 will generally place the origin in the center. The default `NULL` will choose a value based on the different spiral types (default = `NULL`)
#' @param rev If `TRUE` the vertex with the highest index will be placed in the centre of the spiral (default = `FALSE`)
#'
#' @return An igraph layout
#'
#' @export
#'
#' @examples
#'
#'library(igraph)
#'
#' g  <- igraph::sample_gnp(100, 1/100)
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
  if(length(which(type%in%c("Archimedean","Bernoulli","Fermat","Euler")))==0){
    stop("Type must be one of: Archimedean, Bernoulli, Fermat, Euler")
  }

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
    if(is.null(b)){
      b = .5
    }
    res1  <- ceiling(res/2)
    res2  <- floor(res/2)
    theta <- seq(0,pi/2,length.out = res1)
    df1   <- data.frame(t=theta,x=NA,y=NA)
    df1$x[1] <-0
    df1$y[1] <-0

    dt <- arcs/res1

    for(i in 2:res1){
      df1$x[i] <- df1$x[i-1] + cos(df1$t[i-1]^(b+1)/(b+1)) * dt #*df1$t[i-1]) * dt
      df1$y[i] <- df1$y[i-1] + sin(df1$t[i-1]^(b+1)/(b+1)) * dt #*df1$t[i-1]) * dt
      df1$t[i] <- df1$t[i-1]+dt
    }
    df <- matrix(cbind(c(-rev(df1$x[1:res2]),df1$x),c(-rev(df1$y[1:res2]),df1$y)),ncol = 2)
  }

  df     <- df[seq(1, res, length.out = N),]
  df[,1] <- elascer(df[,1],lo = 0.1, hi = 0.9,boundaryPrecision = NA)
  df[,2] <- elascer(df[,2],lo = 0.1, hi = 0.9,boundaryPrecision = NA)

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
#' @inheritParams layout_as_spiral
#' @inheritParams plotNET_groupColour
#' @param markTimeBy Include a vector that indicates time. The time will be displayed on the plot. Pass `TRUE` to generate auto labels (experimental)
#' @param labelSize The size of text in the annotation labels (default = `3`)
#' @param curvature The `curvature` parameter for edges see [geom_curve()] (default = `-0.7`)
#' @param angle The `angle` parameter for edges see [geom_curve()] (default = `90`)
#' @param showArrows Show arrows at the end of the edges? (default = `FALSE`)
#' @param title A title for the plot
#' @param subtitle A subtitle for the plot
#' @param showEpochLegend Should a legend be shown for the epoch colours? (default = `TRUE`)
#' @param markEpochsBy A vector of length `vcount(g)` indicating epochs or groups (default = `NULL`)
#' @param epochColours A vector of length `vcount(g)` with colour codes (default = `NULL`)
#' @param epochLabel A title for the epoch legend (default = `"Epoch"`)
#' @param showSizeLegend Should a legend be shown for the size of the nodes? (default = `FALSE`)
#' @param sizeLabel Guide label, use it to indicate if `V(g)$size` represents some measure, e.g. [igraph::degree()], or, [igraph::hub_score()], [igraph::strength()] (default = `"Size"`)
#' @param scaleVertexSize Scale the size of the vertices by setting a range for [ggplot2::scale_size()]. This will not affect the numbers on the size legend (default = `c(1,6)`)
#' @param vertexBorderColour Draw a border around the vertices. Pass `NULL` to use the same colour as the fill colour (default = `"black"`)
#' @param edgeColourLabel Use to indicate if `E(g)$color` represents color coding based on some property. (default = `"Weight"`)
#' @param scaleEdgeSize Scale the size of the edges by a constant: `E(g)$width * scaleEdgeSize` (default = `1/5`)
#' @param showEdgeColourLegend Should a legend be shown for the colour of the edges? (default = `FALSE`)
#' @param edgeColourByEpoch Should edges that connect to the same epoch be assigned the epoch colour? This will ignore edge colour info in `E(g)$color`. (default = `TRUE`)
#' @param defaultEdgeColour Colour of edges that do not connect to the same epoch (default = `"grey70"`)
#' @param doPlot Produce a plot? (default = `TRUE`)
#'
#' @note To keep the igraph object, use the layout function [layout_as_spiral()] when plotting the graph.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#'
#' library(igraph)
#'
#' g  <- igraph::sample_gnp(200, 1/20)
#' V(g)$size <- degree(g)
#' make_spiral_graph(g, markTimeBy = TRUE, showSizeLegend = TRUE, sizeLabel = "Node degree")
#'
make_spiral_graph <- function(g,
                              type = "Archimedean",
                              arcs = 6,
                              a = 1,
                              b = NULL,
                              rev= FALSE,
                              curvature = -0.6,
                              angle = 90,
                              markTimeBy = NULL,
                              labelSize = 3,
                              alphaV = 1,
                              alphaE = .6,
                              showArrows = FALSE,
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
                              edgeColourLabel = "Weight",
                              showEdgeColourLegend = FALSE,
                              edgeColourByEpoch = TRUE,
                              defaultEdgeColour = "grey70",
                              doPlot = TRUE){

  g$layout <- layout_as_spiral(g, type = type, arcs = arcs, a = a, b = b, rev = rev)

  if(is.null(markTimeBy)){
    if(!is.null(markEpochsBy)){
      grIDs <- ts_changeindex(markEpochsBy)
      tbreaks <- unique(sort(c(grIDs$xmax,grIDs$xmax)))
      tlabels <- ceiling(tbreaks)
      #markTimeBy <- TRUE
    } else {
      x <- 1:igraph::vcount(g)
      v <- seq(1,igraph::vcount(g),by=igraph::vcount(g)/arcs)
      tbreaks <- c(1,which(diff(findInterval(x, v))!=0),igraph::vcount(g))
      #tbreaks <- which(diff(c(g$layout[1,1],g$layout[,2],g$layout[1,2])>=g$layout[1,2])!=0)
      if(max(tbreaks)!=igraph::vcount(g)){
        tbreaks[which.max(tbreaks)]<-igraph::vcount(g)
        tlabels <- paste(tbreaks)
      }
      if(min(tbreaks)>1){
        tbreaks<- c(1,tbreaks)
        tlabels <- paste(tbreaks)
      }
    }
  } else {
    if(is.numeric(markTimeBy)){
      if(all(markTimeBy%in%1:igraph::vcount(g))){
        tbreaks <- unique(markTimeBy)
      }
      if(!is.null(names(markTimeBy))){
        tlabels <- unique(names(markTimeBy))
      } else {
        tlabels <- paste(tbreaks)
      }
    } else {
      markTimeBy <- TRUE
    }
  }

  if(is.logical(markTimeBy)){
    if(markTimeBy){
      if(type == "Euler"){
        tbreaks <- which(diff(as.numeric(cut(g$layout[,2], breaks = arcs,include.lowest = TRUE)))!=0)
        tbreaks[1] <- 1
        tbreaks[which.max(tbreaks)] <- igraph::vcount(g)
      } else {
        tbreaks <- which(diff(c(g$layout[1,1],g$layout[,2],g$layout[1,2])>=g$layout[1,2])!=0)
      }
      if(max(tbreaks)>igraph::vcount(g)){tbreaks[which.max(tbreaks)]<-igraph::vcount(g)}
      if(max(tbreaks)<igraph::vcount(g)){tbreaks<- c(tbreaks,igraph::vcount(g))}
      if(min(tbreaks)>1){tbreaks<- c(1,tbreaks)}
      tlabels <- paste(tbreaks)
    }
  }

  if(max(tbreaks)!=igraph::vcount(g)){
    tbreaks[which.max(tbreaks)]<-igraph::vcount(g)
    #tlabels <- paste(tbreaks)
    tlabels <- tlabels[tbreaks]
  }

  if(min(tbreaks)>1){
    tbreaks<- c(1,tbreaks)
    #tlabels <- paste(tbreaks)
    tlabels <- tlabels[tbreaks]
  }

  if(length(which(diff(tbreaks)==1))>0){
    tbreaks <- tbreaks[-(which(diff(tbreaks)==1)+1)]
    #tlabels <- paste(tbreaks)
    tlabels <- tlabels[tbreaks]
  }


  if(is.null(markEpochsBy)){
    markEpochsBy <- character(igraph::vcount(g))
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
                           colourE = edgeColourByEpoch,
                           # alphaV = alphaV,
                           # aplhaE = alphaE,
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
  gEdges$colorVar <- as.numeric_discrete(gEdges$color)

  if(is.null(vertexBorderColour))(
    vBc <- gNodes$colour
  ) else {
    if(length(vertexBorderColour)==1|length(vertexBorderColour)==NROW(gNodes)){
      vBc <- vertexBorderColour
    } else {
      warning("Invalid value(s) for vertexBorderColour, using default.")
      vertexBorderColour <- "black"
    }
  }

  # Fix same coords
  sameID <- which(gEdges$from.x==gEdges$to.x)
  if(length(sameID)>0){
    gNodes$V2[gEdges$from[sameID]] <- gNodes$V2[gEdges$from[sameID]]+mean(diff(gNodes$V2))/2
    gEdges$to.x[sameID] <- gEdges$to.x[sameID]+mean(diff(gEdges$to.x))/2
    gEdges$to.y[sameID] <- gEdges$to.y[sameID]+mean(diff(gEdges$to.y))/2
  }


  gEdges$curvature <- curvature
  if(type=="Euler"){
    if(curvature<0){
      gEdges$curvature[gEdges$from.y<.5|gEdges$from.x>.5] <- -1*curvature
    } else {
      gEdges$curvature[gEdges$from.y>.5|gEdges$from.x<.5] <- -1*curvature
    }
  }

  ar <- NULL
  if(igraph::is_directed(g)){
    if(showArrows){
      ar <- ggplot2::arrow(type = "closed", ends = "last", length = ggplot2::unit(mean(gNodes$size, na.rm = TRUE)*scaleVertexSize[2], "mm"))
    }
  }

  gg <- ggplot(gNodes,aes(x=V1,y=V2)) +
    geom_curve(data=gEdges, aes(x = from.x,
                                xend = to.x,
                                y = from.y,
                                yend = to.y,
                                colour = colorVar),
               curvature = curvature,
               arrow = ar,
               angle = angle,
               size = gEdges$width * scaleEdgeSize,
               alpha = alphaE) +
    geom_point(aes(fill = labels, size = size), pch=21, colour = vBc, alpha = alphaV) +
    ggtitle(label = title, subtitle = subtitle) +
    scale_fill_manual(epochLabel, values = unique(gNodes$colour)) +
    scale_size(sizeLabel, range = scaleVertexSize) +
    scale_color_gradientn(edgeColourLabel, colours = unique(gEdges$color))

  if(showEpochLegend){
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

  if(!is.null(markTimeBy)){
    gg <- gg + annotate("label", x=gNodes$V1[tbreaks], y=gNodes$V2[tbreaks], label = tlabels, size=labelSize)
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

#' Get some nice colours
#'
#' @param Ncols Either an integer representing the number of requested colours, or, a numeric vector of integers between 1-58 to select specific colours. Run [getColours()] without arguments to see a plot of the colours that are available.
#' @param continuous Return a discrete vector of colours, or, a function that represents a gradient between 2 or more colours? If `TRUE` then argument `Ncols` must be a numeric vector greater than length `2`. **NOTE:** The input to the gradient function must be (re-)scaled to fall between 0 and 1, e.g. using [elascer()]. (default = `FALSE`)
#' @param Dcols If `continuous` is set to `TRUE`, this should be a vector of the same length as `Ncols` representing the relative distances between the colours in the gradient function using values between 0 and 1. (default = `c(0,1)`)
#'
#' @return A list of colours
#' @export
#'
#' @examples
#'
#' # This will plot all available colours with their numbers
#' getColours()
#'
#' # Get a specific number of colours
#' getColours(5)
#'
#' # Get specific colours
#' getColours(c(4,7,1,40))
#'
#' # Make a gradient from colour number 4 to 44 via 7
#' gradFunc <- getColours(Ncols = c(4,7,44), continuous = TRUE, Dcols = c(0,.5,1))
#'
#' df <- data.frame(x=1:50, y=sort(rnorm(50)))
#' # Make sure the input is on a scale of 0-1
#' df$ycol <- elascer(df$y)
#'
#' library(ggplot2)
#' ggplot(df, aes(x=x,y=y,colour=ycol)) +
#' geom_point() +
#' scale_colour_gradientn("Gradient",colours = gradFunc(df$ycol)) +
#' theme_bw()
#'
#' # Make a gradient from colour number 4, to 9, to 7, to 36, to 44
#' gradFunc <- getColours(Ncols = c(4,9,7,36,44), continuous = TRUE, Dcols = c(0,.33,.5,.66,1))
#'
#' ggplot(df, aes(x=x,y=y,colour=ycol)) +
#' geom_point() +
#' scale_colour_gradientn("Gradient",colours = gradFunc(df$ycol)) +
#' theme_bw()
#'
getColours <- function(Ncols, continuous = FALSE, Dcols = c(0,1)){

  pl <- c("#A65141","#ABB2C4","#C89284","#7896BC","#E1CAC2","#536489","#ECE3CD","#575E7A","#BEBEC4","#C4AD75","#30261E","#BC964D","#3D434D","#D29E55","#B4B6A7","#9B7342","#E8D69D","#48211A","#EBDAB4","#3D4B68","#E8DCCF","#3E688C","#9D3D2D","#507074","#9A756C","#546959","#93725F","#858563","#E1D3B1","#A6B4BB","#78828C","#AA9C6A","#91A193","#CDB16E","#1B528D","#B8A98F","#6E432D","#4F5B71","#D9B196","#20304F","#827561","#98A39E","#8B4F31","#7E7B65","#C1B6A1","#775E45","#C0B3A2","#5A524A","#BBA27F","#3A3935","#C9BDA3","#626961","#8A4F42","#8D8A77","#947B5E","#5D3C30","#AA8470","#493A2C")

  cols <- NULL
  if(!continuous){
    if(missing(Ncols)){
      df <- expand.grid(x=1:8,y=1:8)[1:58,]
      plot(df$x,df$y,pch=15,col=getColours(Ncols=58),cex=5, axes = FALSE, ann = FALSE)
      graphics::text(df$x,df$y,paste(1:58),col="white")
    } else {
      if(all(Ncols%[]%c(1,58))){
        if(length(Ncols)==1){
          cols <- pl[1:Ncols]
        } else {
          if(length(Ncols)%[]%c(2,58)){
            cols <- pl[Ncols]
          }
        } # Ncols >= 2

      } else {
        stop("Currently the max. number of colours is 58.")
      }
      return(cols)
    } # missing(Ncols)
  } else { # continuous
    if(length(Ncols)>=2 & all(invctr::is.wholenumber(Ncols)) & all(Ncols%[]%c(1,58))){
      if(length(Dcols)==length(Ncols)&all(Dcols%[]%c(0,1))){
        cols <- scales::gradient_n_pal(colours = pl[Ncols], values = Dcols)
      } else {
        stop("If 'continuous = TRUE' then 'Dcols' must have the same length as 'Ncols' and contain values between 0 and 1.")
      }
    } else {
      stop("If 'continuous = TRUE' then 'Ncols' must be a vector of ast least two integers between 1 and 58.")
    }
  }
  return(cols)
}


# Empty plot
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

