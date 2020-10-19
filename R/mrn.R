# package casnet ----
#
# Multiplex Recurrence Networks ----
#
# mrn functions


#' Multiplex Recurrence Network
#'
#' @description This function will create a Multiplex Recurrence Network from a list of [igraph] objects that can be considered the layers of a network. The layers must have the same number of nodes. There are two modes of operation: *Layer similarity* (`MRNweightedBy` is set to `"InterLayerMI"`, `"InterLayerCor"`, or `"EdgeOvelap"`) and *Layer importance* (not implemented yet). The former generates weighted MRN based on _Interlayer Mutual Information_, _Interlayer Correlation_, or _Edge Overlap_.
#'
#' @inheritParams ts_windower
#' @inheritParams rn
#' @param layers A list of igraph objects representing the layers of the multiplex network. The layer networks must all have the same number of vertices.
#' @param MRNweightedBy The measure to be used to evaluate the average structural similarities between the layers of the network. Valid options are: `"InterLayerMI"` (Mutual information based on similarity of the vertex degree across layers), `"EdgeOverlap"` (proportion of vertices sharing the same edges across layers). Choosing `"InterLayerMI"`, `"InterlayerCor"`, or `"EdgeOverlap"` will decide which measure is displayed in the plot of the Multiplex RN, all measures will always be returned in the numerical output.
#' @param win The window size passed to [casnet::ts_windower()] in which to evaluate `"InterLayerMI"`, `"InterLayerCor"`, or `"EdgeOvelap"`. (default = `NA`).
#' @param step The stepsize for the sliding window (default = `NA`).
#' @param overlap The window overlap passed to [casnet::ts_windower()] if `MRNweightedBy` is `"InterLayerMI"` or `"EdgeOvelap"`. The value of `step` will be ignored if `overlap` is not `NA`. (default = `NA`).
#' @param doPlot Plot the multiplex recurrence network (default = `FALSE`).
#' @param silent Silent-ish mode. (default = FALSE).
#'
#' @return A list object with fields:
#' * _interlayerMI_ - One or more matrices with edge weights between layers that represent the interlayer Mutual Information.
#' * _interlayerCor_ - One or more matrices with edge weights between layers that represent the Pearson correlation between vertex degrees of layers.
#' * _edgeOverlap_ - One or more matrices with edge weights between layers that represent the overlapping edges between layers.
#' * _meanValues_ - One or more matrices that represent the means and SDs of the interlayer Mutual Information, absolute interlayer correlation and edge overlap. Ther measure `eo_joint` refers to the number of edges shared among _all_ layers of the MRN.
#'
#' @export
#'
#' @examples
#'
#' # Create some layers
#' library(igraph)
#'
#' layers <- list(g1 = igraph::sample_smallworld(1, 100, 5, 0.05),
#' g2 = igraph::sample_smallworld(1, 100, 5, 0.5),
#' g3 = igraph::sample_smallworld(1, 100, 5, 1))
#'
#' mrn(layers = layers)
#'
#'
mrn <- function(layers,
                MRNweightedBy = c("InterlayerMI", "InterlayerCor", "Edgeoverlap")[1],
                win = NA,
                step = NA,
                overlap = NA,
                alignment = "r",
                cumulative = FALSE,
                doPlot = FALSE,
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

  if(is.na(win)){
    win <- Nsize
  }
  if(is.na(step)){
    if(win==Nsize){
      step <- 1
    } else {
      step <- win
    }
  }

  wIndexList <- ts_windower(y = 1:Nsize, win = win, step = step, overlap = overlap, adjustY = NA)

  func <- "mi_interlayer"
  weighted <- igraph::is.weighted(layers[[1]])
  directed <- igraph::is.directed(layers[[1]])

  if(directed){
    #mode <- "lower"
    combis <- getPairs(seq_along(layers))
    layers <- plyr::llply(layers, function(g) igraph::as.undirected(g))
    combis <- getPairs(seq_along(layers),seq_along(layers))
    mode <- "directed"
    if(cumulative){
      mode <- "upper"
    }
  } else {
    combis <- getPairs(seq_along(layers))
    mode <- "upper"
  }

  if(func%in%"mi_interlayer"){

    cat(paste("\nWelcome to the multiplex... in layer similarity mode!\n\n"))

    out_mi <- out_eo <- out_means <- list()

    pb <- progress::progress_bar$new(total = length(wIndexList), force = TRUE)

    for(w in seq_along(wIndexList)){

      pb$tick()
      mp <- eo <- matrix(nrow=length(layers), ncol=length(layers), dimnames = list(names(layers),names(layers)))

      for(i in seq_along(combis$X1)){

        edgeFrame1 <- igraph::as_data_frame(layers[[combis$X1[i]]],"edges")
        edgeFrame1 <- edgeFrame1 %>% filter(edgeFrame1$from%[]%range(wIndexList[[w]])&edgeFrame1$to%[]%range(wIndexList[[w]]))
        ga  <- graph_from_data_frame(edgeFrame1, directed = directed)

        edgeFrame2 <- igraph::as_data_frame(layers[[combis$X2[i]]],"edges")
        edgeFrame2 <- edgeFrame2 %>% filter(edgeFrame2$from%[]%range(wIndexList[[w]])&edgeFrame2$to%[]%range(wIndexList[[w]]))
        gb  <- graph_from_data_frame(edgeFrame2, directed = directed)

        if(any(igraph::ecount(ga)==0,igraph::ecount(gb)==0)){
          mp[combis$X1[i],combis$X2[i]] <- NA
          eo[combis$X1[i],combis$X2[i]] <- NA
        } else {
        mp[combis$X1[i],combis$X2[i]] <- mi_interlayer(ga,gb)
        eo[combis$X1[i],combis$X2[i]] <- igraph::ecount(ga %s% gb) / (igraph::ecount(ga) + igraph::ecount(gb))
        }

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

  } else {
    stop("Not implemented yet!")
  }

  cat("\n~~~o~~o~~casnet~~o~~o~~~\n")

  return(list(interlayerMI  = out_mi,
              edgeOverlap   = out_eo,
              meanValues    = out_means,
              Nsize         = Nsize))
}



#' Mutliplex Recurrence Network Plot
#'
#' @description This function will plot a Multiplex Recurrence Network from a list of [igraph] objects that can be considered the layers of a network. or based on the output of function [mrn()]. The layers must have the same number of nodes.
#'
#' @inheritParams make_spiral_graph
#' @inheritParams plotNET_groupColour
#' @inheritParams mrn
#' @param MRN The output from function [mrn()]
#' @param MRNweightedBy The measure to be used to evaluate the average structural similarities between the layers of the network. Valid options are: `"InterLayerMI"` (Mutual information based on similarity of the vertex degree across layers), `"EdgeOverlap"` (proportion of vertices sharing the same edges across layers). Choosing `"InterLayerMI"` or `"EdgeOverlap"` will decide which measure is displayed in the plot of the Multiplex RN, both measures will always be returned in the numerical output.
#' @param doPlot Plot the multiplex recurrence network (default = `TRUE`).
#' @param doSave Save the plots.
#' @param coords A data frame with layout coordinastes generated by calling any of the [igraph] layout functions. If `NA` a circle layout will; be generated (default = `NA`)
#' @param RNnodes Should the vertices of the MRN represent a plot of the RN of the layers? This is recommended only for a small numbers of vertices. (default = `FALSE``)
#' @param vertexSizeBy A valid [igraph] function that calculates node based measures, or a numeric constant. (default = `"degree"`)
#' @param vertexColour A vector of colours for the vertices. If this is a named list, names will be displayed in the legend.
#' @param showVertexLegend Show the vertex colour legend?
#' @param createAnimation If `createAnimation = TRUE` *and* `doPlot = TRUE` *and* a windowed analysis is conducted, an animation will be produced using either package `gganimate` (if `useImageMagick = FALSE`) or `animation` (if `useImageMagick = TRUE`). The main difference is that `gganimate` has nice animation transition features, but plots the MRN using [ggplot2], which does not have great options for displaying the nodes as images. With package `animation` a sequence of [igraph] plots will be converted to an animation. If `doSave = TRUE` the animation will be saved in `imageDir` as an animated gif by calling either [gganimate::anim_save()], or [animation::saveGIF()] (default = `FALSE`)
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
#' #' # Create some layers
#' library(igraph)
#'
#' layers <- list(g1 = igraph::sample_smallworld(1, 100, 5, 0.05),
#' g2 = igraph::sample_smallworld(1, 100, 5, 0.5),
#' g3 = igraph::sample_smallworld(1, 100, 5, 1))
#'
#' mrn_plot(layers = layers,showEdgeColourLegend=TRUE)
#'
mrn_plot     <- function(layers = NA,
                         MRN = NA,
                         MRNweightedBy = c("InterLayerMI", "EdgeOverlap")[1],
                         win = NA,
                         step = NA,
                         overlap = NA,
                         alignment =  "r",
                         cumulative = FALSE,
                         doPlot = TRUE,
                         doSave = FALSE,
                         coords = NA,
                         RNnodes = FALSE,
                         vertexSizeBy = "degree",
                         scaleVertexSize = c(.01,5),
                         vertexColour = NA,
                         vertexBorderColour = "black",
                         showVertexLegend = TRUE,
                         showSizeLegend = FALSE,
                         alphaV = .7,
                         scaleEdgeSize = 1/5,
                         alphaE = .5,
                         showEdgeColourLegend = FALSE,
                         curvature = -0.2,
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

  if(is.na(MRN)){

    if(missing(layers)){
      stop("Need either a list of graphs (layers), or output from function mrn()!")
    }

    MRN <- mrn(layers = layers,
               MRNweightedBy = MRNweightedBy,
               win = win,
               step = step,
               overlap = overlap,
               alignment =  alignment,
               cumulative = cumulative,
               doPlot = FALSE,
               silent = silent)

  } else {
    if(all(names(MRN)%in%c("interlayerMI","edgeOverlap","meanValues","Nsize"))){
      cat("\n~~~o~~o~~casnet~~o~~o~~~\n")
      cat(paste("\nWelcome to the multiplex... in layer similarity mode!\n\n"))
    } else {
      stop("MRN is not an object output by function mrn()")
    }
  }

  Nsize <- MRN$Nsize

  if(is.na(win)){
    win <- Nsize
  }
  if(is.na(step)){
    if(win==Nsize){
      step <- 1
    } else {
      step <- win
    }
  }

  if(missing(layers)){
    vertexColour <- getColours(NCOL(MRN$interlayerMI[[1]]))
    names(vertexColour) <- colnames(MRN$interlayerMI[[1]])
  } else {
    vertexColour <- getColours(length(layers))
    names(vertexColour) <- names(layers)
  }

  wIndexList <- ts_windower(y = 1:Nsize, win = win, step = step, overlap = overlap, adjustY = NA)

  func <- "mi_interlayer"

  if(func%in%"mi_interlayer"){

    if(doPlot){

      gList <- MRNlist <- list()

      if(MRNweightedBy%in%"InterLayerMI"){
        MRNlist <- MRN$interlayerMI
      } else {
        MRNlist <-MRN$edgeOverlap
      }

      if(RNnodes){
        checkPkg("png")

        if(!missing(layers)){
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

            gg <- gg + ggimage::theme_transparent() +
              theme(
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
            if(useImageMagick){
              g_rast[[f]] <- png::readPNG(file.path(imageDir,paste0(names(gSpiro)[f],".png")), native=TRUE)
            } else {
              g_rast[[f]] <- ggimage::image_read2(file.path(imageDir,paste0(names(gSpiro)[f],".png")))
            }
          }
        } # RNnodes
      } # missing

      for(w in seq_along(MRNlist)){

        mp_net <- graph_from_adjacency_matrix(adjmatrix = MRNlist[[w]],
                                              mode      = "upper",
                                              diag      = FALSE,
                                              weighted  = TRUE)

        if(is.na(coords)){
          coord <- layout_in_circle(mp_net)
        } else {
          coord <- coords
        }

        if(RNnodes){
          V(mp_net)$raster  <- g_rast
          V(mp_net)$shape   <- "raster"
          V(mp_net)$size    <- 10
          V(mp_net)$size2   <- 10
        } else {
          if(is.numeric(vertexSizeBy)){
            V(mp_net)$size    <- vertexSizeBy
          } else {
            if(vertexSizeBy%in%lsf.str("package:igraph")){
              V(mp_net)$size  <- eval(parse(text=paste(vertexSizeBy,"(mp_net)")))
            } else {
              V(mp_net)$size    <- 1
            }
          }
        }

        V(mp_net)$name  <- names(layers)
        V(mp_net)$label <- V(mp_net)$name
        V(mp_net)$alpha <- alphaV
        V(mp_net)$color <- vertexColour
        V(mp_net)$label.family <- "sans"
        V(mp_net)$label.font <- 2
        mp_net$layout <- coord

        E(mp_net)$width   <- elascer(E(mp_net)$weight,lo = 1,hi = 5)
        E(mp_net)$color   <- getColours(Ncols = c(3,7,4), continuous = TRUE, Dcols = c(0,.5,1))(E(mp_net)$width/5)

        if(useImageMagick){

          if(!noParts){
            plot(mp_net)
          }

          MRNlist[[w]] <- mp_net

        } else {

          gNodes        <- as.data.frame(mp_net$layout,what = "vertices")
          gNodes$ID     <- as.numeric(igraph::V(mp_net))
          gNodes$name   <- V(mp_net)$name
          gNodes$color  <- V(mp_net)$color
          gNodes$alpha  <- V(mp_net)$alpha
          gNodes$size   <- V(mp_net)$size

          if(RNnodes){
            gNodes$image  <- paste0(file.path(imageDir,paste0(gNodes$name,".png")))
          }

          gEdges        <- igraph::get.data.frame(mp_net,what = "edges")
          if(any(is.character(gEdges$from))){
            cName <- "name"
          } else {
            cName <- "ID"
          }
          gEdges$from.x <- gNodes$V1[match(gEdges$from, gNodes[[cName]])]
          gEdges$from.y <- gNodes$V2[match(gEdges$from, gNodes[[cName]])]
          gEdges$to.x   <- gNodes$V1[match(gEdges$to, gNodes[[cName]])]
          gEdges$to.y   <- gNodes$V2[match(gEdges$to, gNodes[[cName]])]

          # scales::gr colour_gradient2(low = "steelblue",
          #                        high     = "red3",
          #                        mid      = "grey90",
          #                        na.value = scales::muted("slategray4"),
          #                        midpoint = median(gEdges$weight))

          edgeCols        <- gEdges$color
          gEdges$size     <- gEdges$weight * scaleEdgeSize
          gEdges$weight   <- factor(round(gEdges$weight, digits = 4))
          names(edgeCols) <- paste(gEdges$weight)



          gg <- ggplot(gNodes,aes_(x=~V1, y=~V2)) +
            geom_curve(data=gEdges, aes_(x = ~from.x,
                                        xend = ~to.x,
                                        y = ~from.y,
                                        yend = ~to.y,
                                        colour = ~weight,
                                        size = ~size),
                       alpha = alphaE,
                       curvature = curvature)

          if(RNnodes){
            gg <- gg + ggimage::geom_image(aes(image=image), size = .2)

          } else {
            gg <- gg + geom_point(aes_(size = ~size, fill = ~name), pch=21, colour = vertexBorderColour, alpha = alphaV) +
              scale_fill_manual("Layers", values = vertexColour)
          }

          gg <- gg +
            scale_size_binned(range = scaleVertexSize) +
            scale_colour_manual(paste("Weighted by",MRNweightedBy), values = edgeCols) +
            labs(title = names(MRNlist)[w]) +
            coord_cartesian(clip="off") +
            theme_void() + theme(plot.margin = margin(50,50,50,50, unit = "pt"),
                                 legend.margin = margin(l = 20, unit = "pt"),
                                 plot.title = element_text(margin = margin(b=20)))

          if(showVertexLegend){
            if(RNnodes){

              gg <- gg + geom_text(aes_(label=~name))

             } else {
               gg <- gg + guides(fill = guide_legend(title.position = "top",
                                                     byrow = TRUE,
                                                     nrow=2,
                                                     override.aes = list(size=5,order = 0)))
            }
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

          MRNlist[[w]] <- gg
          gList[[w]] <- list(gNodes = gNodes,
                             gEdges = gEdges)

        } # useImageMagick

      } # w

      if(!createAnimation){
        return(list(MRN           = MRNlist,
                    interlayerMI  = MRN$interlayerMI,
                    edgeOverlap   = MRN$edgeOverlap,
                    meanValues    = MRN$out_means))
      } else {

        if(useImageMagick){

          checkPkg("animation")

          animation::ani.options(interval = stateLength, imgdir = imageDir, loop = loopAnimation, nmax = length(MRNlist), ani.width = gifWidth, ani.res = gifRes)

          if(doSave){

            animation::saveGIF(
              for(i in seq_along(MRNlist)){
                plot(MRNlist[[i]], xlab = names(MRNlist)[i])
                animation::ani.pause()
              }, img.name = paste0(names(MRNlist)[i]), movie.name = paste0("MRN_animation_win",win,"_step",step,".gif")
            )

          } else {
            for(i in seq_along(MRNlist)){
              plot(MRNlist[[i]], xlab = names(MRNlist)[i])
              animation::ani.pause()
            }
          }

          return(list(MRN           = MRNlist,
                      interlayerMI  = MRN$interlayerMI,
                      edgeOverlap   =  MRN$edgeOverlap,
                      meanValues    =  MRN$out_means))

        } else {

          checkPkg("gganimate")

          names(gList) <- names(wIndexList)

          gNodes <- plyr::ldply(gList , function(g) g$gNodes)
          gEdges <- plyr::ldply(gList , function(g) g$gEdges)

          gg <- ggplot(gNodes, aes_(x = ~V1, y = ~V2)) +
            geom_curve(data = gEdges, aes_(x = ~from.x,
                                    xend = ~to.x,
                                    y = ~from.y,
                                    yend = ~to.y,
                                    size = ~weight,
                                    colour = ~weight), curvature = 0)
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
            gganimate::transition_states(factor(~.id),
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

          return(list(MRN              = MRNlist,
                      MRNanimationData = gList,
                      MRNanimationGG   = invisible(gg),
                      interlayerMI  =  MRN$interlayerMI,
                      edgeOverlap   =  MRN$edgeOverlap,
                      meanValues    =  MRN$out_means))
        }
      }

    } else { #doPlot
      return(MRN)
    }

  } else {    # NOT mi_interlayer

    cat(paste("\nWelcome to the multiplex... in layer importance mode!\n\n"))
    cat(paste("\n\nNOT IMPLEMENTED\n\n"))
  } # rankDC
  cat("\n\n~~~o~~o~~casnet~~o~~o~~~\n")
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

  checkPkg("infotheo")

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
  attr(mif_out,"lags")   <- lags
  attr(mif_out,"nbins")   <- nbins

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

  checkPkg("infotheo")

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
#' @note If the networks are weighted the strength distribution will be used instead of the the degree distribution.
#'
#' @family Redundancy measures (mutual information)
#' @family Multiplex Recurrence Networks
#'
mi_interlayer <- function(g0,g1, probTable=FALSE){

  checkPkg("infotheo")

  if(any(E(g0)$weight<0)){E(g0)$weight <- abs(E(g0)$weight)}
  if(any(E(g1)$weight<0)){E(g1)$weight <- abs(E(g1)$weight)}
  if(igraph::is.weighted(g0)&igraph::is.weighted(g1)){
    s0  <- strength(g0)
    s0m <- max(s0, na.rm = TRUE)
    d0  <- graphics::hist(s0,seq(0,(ceiling(s0m)+1)), include.lowest = TRUE, plot = FALSE, warn.unused = FALSE)$density
    s1  <- strength(g1)
    s1m <- max(s1, na.rm = TRUE)
    d1  <- graphics::hist(s1,seq(0,(ceiling(s1m)+1)), include.lowest = TRUE, plot = FALSE, warn.unused = FALSE)$density
    equal <- data.frame(ts_trimfill(d0,d1))
    p01 <- graphics::hist(x = igraph::strength(g0),breaks = seq(-.5,((NROW(equal))-.5)),plot=FALSE, warn.unused = FALSE)$counts
    p10 <- graphics::hist(x = igraph::strength(g1),breaks = seq(-.5,((NROW(equal))-.5)),plot=FALSE, warn.unused = FALSE)$counts
  } else {
    d0    <- igraph::degree_distribution(g0)
    d1    <- igraph::degree_distribution(g1)
    equal <- data.frame(ts_trimfill(d0,d1))
    p01 <- graphics::hist(x = igraph::degree(g0),breaks = seq(-.5,((NROW(equal))-.5)),plot=FALSE, warn.unused = FALSE)$counts
    p10 <- graphics::hist(x = igraph::degree(g1),breaks = seq(-.5,((NROW(equal))-.5)),plot=FALSE, warn.unused = FALSE)$counts
  }

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
