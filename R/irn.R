# Inter-system recurrence networks


#' Inter system recurrence networks
#'
#' @param RMxy A cross recurrence matrix of x and y
#' @param gx An igraph object of recurrence network x
#' @param gy An igraph object of recurrence network y
#' @param merge_seed Seed
#' @inheritParams make_spiral_graph
#' @param colorExy Edge colour x >> y
#' @param colorEyx Edge colour y >> x
#' @param sizeEbyCC Node size by Clustering Coefficient
#'
#' @return A plot
#'
#' @export
#'
irn_plot <- function(RMxy,
                     gx = NULL,
                     gy = NULL,
                     merge_seed = 5,
                     curvature = 0,
                     angle = 0,
                     scaleEdgeSize = 1,
                     alphaE = .1,
                     alphaV = 1,
                     colorExy = "grey30",
                     colorEyx = "grey70",
                     sizeEbyCC = FALSE,
                     doPlot = TRUE,
                     ggplotReturn = FALSE,
                     igraphReturn = FALSE){

  if(doPlot){
    if(any(is.null(c(gx,gy)))){
      stop("Need gx and gy to create inter-system plot.")
    }
  }

  if(any(is.null(c(gx$layout,gy$layout)))){
    stop("One (or both) of the graphs do not have a `$layout` attribute!")
  }

  RMyx <- RMxy
  RMxy[lower.tri(RMxy)] <- 0
  Exy <- mat_mat2ind(as.matrix(RMxy)) %>% filter(value!=0)
  Exy$toY <- Exy$Var2 + NROW(RMyx)

  RMyx[upper.tri(RMyx)] <- 0
  Eyx <- mat_mat2ind(as.matrix(RMyx)) %>% filter(value!=0)
  Eyx$fromY <- Eyx$Var2 + NROW(RMxy)
  rm(RMxy,RMyx)

  # Graphs
  graphs <- list(gx, gy)
  set.seed(merge_seed)
  lay <- igraph::merge_coords(graphs, list(gx$layout, gy$layout))
  g  <- disjoint_union(graphs)
  V(g)$labels <- NA
  g$layout <- lay
  #plot(gd, vertex.size=2,vertex.label = NA,edge.arrow.size = 0.1)

  gNodes        <- as.data.frame(g$layout)
  gNodes$ID     <- as.numeric(V(g))
  gNodes$color  <- V(g)$color
  #gNodes$labels <- factor(V(g)$groupnum, levels = unique(V(g)$groupnum), labels = unique(V(g)$group))
  gNodes$alpha  <- alphaV
  gNodes$colorVarV <- as.numeric_discrete(gNodes$color)
  gNodes$size   <- V(g)$size%00%.1

  gEdges        <- igraph::get.data.frame(g)
  gEdges$alpha  <- alphaE
  gEdges$size   <- E(g)$size%00%.1

  if(sizeEbyCC){
    sizeXYe <- CCxy
    sizeYXe <- CCyx
  } else {
    sizeXYe <- min(E(g)$size, na.rm = TRUE)
    sizeYXe <- min(E(g)$size, na.rm = TRUE)
  }

  gEdges <- gEdges %>%
    add_row(from = Exy$Var1, to = Exy$toY, weight = Exy$value, color = colorExy, alpha = alphaE/2, size = sizeXYe)
  gEdges <- gEdges %>%
    add_row(from = Eyx$fromY, to = Eyx$Var1, weight = Eyx$value, color = colorEyx, alpha = alphaE/2, size = sizeYXe)

  gEdges$colorVar <- as.numeric_discrete(gEdges$color) #gEdges$color)

  gEdges$from.x <- gNodes$V1[match(gEdges$from, gNodes$ID)]
  gEdges$from.y <- gNodes$V2[match(gEdges$from, gNodes$ID)]
  gEdges$to.x   <- gNodes$V1[match(gEdges$to,   gNodes$ID)]
  gEdges$to.y   <- gNodes$V2[match(gEdges$to,   gNodes$ID)]

  # Fix same coords
  sameID <- which(gEdges$from.x==gEdges$to.x)
  if(length(sameID)>0){
    gNodes$V2[gEdges$from[sameID]] <- gNodes$V2[gEdges$from[sameID]]+mean(diff(gNodes$V2))/2
    gEdges$to.x[sameID] <- gEdges$to.x[sameID]+mean(diff(gEdges$to.x))/2
    gEdges$to.y[sameID] <- gEdges$to.y[sameID]+mean(diff(gEdges$to.y))/2
  }


  IR <- graph_from_edgelist(as.matrix(cbind(gEdges$from,gEdges$to)))
  gNodes$size <- degree(IR)

  ggISRN <- ggplot(data = gNodes,aes(x=V1,y=V2)) +
    geom_point(aes(fill = colorVarV, alpha = alpha, size = size), pch = 21, colour = "black") +
    geom_curve(data=gEdges, aes(x      = from.x,
                                xend   = to.x,
                                y      = from.y,
                                yend   = to.y,
                                colour = colorVar,
                                alpha  = alpha), #,linewidth = size),
               curvature = curvature,
               angle = angle) +
    scale_fill_gradientn(colors = unique(gNodes$color), guide = "none") +
    scale_colour_gradientn(colors = unique(gEdges$color), guide = "none") +
    scale_alpha(guide = "none") +
    scale_size(guide = "none") +
    theme_void()

  if(doPlot){
    plot(ggISRN)
  }

  if(ggplotReturn|igraphReturn){
    return(list(ggISRN,IR)[c(ggplotReturn,igraphReturn)])
  }

}


#' Cross Triples
#'
#' @param Ax Adjacency matrix x
#' @param Ay Adjacency matrix y
#' @param Axy Adjacency matrix xy
#'
#' @return list object with triasngles and triads
#'
#' @export
#'
irn_crossTriples <- function(Ax,Ay,Axy){

  Ayx <- Matrix::t(Axy)

  # p != q
  Matrix::diag(Ax) <- 0
  Matrix::diag(Ay) <- 0

  # x -> y
  Kxy <- cross_degree(Axy)

  v_ind <- which(Kxy>0)
  vp_ind <- pq_ind <- qv_ind <- Triads <- Triangles <- vector("list",max(v_ind))

  Triangles[1:max(v_ind)] <- 0
  Triads[1:max(v_ind)] <- 0

  for(v in v_ind){
    vp_ind[[v]] <- which(Ay[v,]>0)%00%NA
    #cat(paste(v,"\n"))

    if(any(is.na(vp_ind[[v]]))){
      pq_ind[[v]] <- 0
      Triads[[v]] <- NA
    } else {
      pq_ind[[v]] <- llply(vp_ind[[v]], function(p){
        ind <- which(Ay[p,]>0)%00%0
        return(ind[ind!=p])
      })
      Triads[[v]] <- length(unlist(pq_ind[[v]]))%00%0
    }

    if(is.na(Triads[[v]])){
      Triads[[v]] <- 0
      Triangles[[v]] <- 0
    } else {
      Triangles[[v]] <- sum(unlist(llply(pq_ind[[v]], function(q){
        qtmp <- llply(q, function(qv){
          ind <- which(Ayx[qv,]>0)%00%NULL
          ind[ind==v]
        })
        return(sum(lengths(qtmp)))
      })), na.rm = TRUE)
    }
  }
  return(list(triads = unlist(Triads),
              triangles = unlist(Triangles)))
}


#' Cross Degree
#'
#' @param Axy Adjacency matrix xy
#'
#' @return list object wih local and global cross-degree
#' @export
#'
irn_crossDegree <- function(Axy){
  Kxy <- degree(graph_from_adjacency_matrix(as.matrix(Axy)), mode = "out", loops = TRUE)
}

#' Cross CLustering Coefficient
#'
#' @param Ax Adjacency matrix x
#' @param Ay Adjacency matrix y
#' @param Axy Adjacency matrix xy
#' @param local local CC
#' @param global global CC
#'
#' @return list object with local and global CC
#'
#' @export
#'
irn_crossClustering <- function(Ax, Ay, Axy, local = FALSE, global = TRUE){

  Ayx <- Matrix::t(Axy)

  # p != q
  Matrix::diag(Ax) <- 0
  Matrix::diag(Ay) <- 0

  # x -> y
  Kxy <- irn_crossDegree(Axy)
  ind_xy <- Kxy>1

  # y -> x
  Kyx <- irn_crossDegree(Ayx)
  ind_yx <- Kyx>1

  Cedge_density_XY <- sum(Axy) * (1/prod(dim(Axy)))
  Cedge_density_YX <- sum(Ayx) * (1/prod(dim(Ayx)))

  # x -> y
  Cxy_local <-  (Matrix::diag(Ayx %*% Ay %*% Axy) / (Kyx * (Kyx-1)))%00%0
  Cxy_local[!ind_xy] <- 0

  Cxy_global <- mean(Cxy_local, na.rm = TRUE)

  # y -> x
  Cyx_local <-  (Matrix::diag(Axy %*% Ax %*% Ayx) / (Kxy * (Kxy-1)))%00%0
  Cyx_local[!ind_yx] <- 0

  Cyx_global <- mean(Cyx_local, na.rm = TRUE)

  return(list(local = cbind(Cxy_local, Cyx_local),
              global = data.frame(direction = c("XY","YX"),
                                  cross_edge_density = c(Cedge_density_XY,Cedge_density_YX),
                                  cross_cc_global = c(Cxy_global, Cyx_global)))[c(local,global)]
  )
}
