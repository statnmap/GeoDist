#' Class for distref.raster
#'
#' @slot r.ref.path path of the raster
#' @slot r.ref.dim dimensions of the raster
#' @slot r.ref.N cells numerotation
#' @slot r.ref.NA.nb numerotation of NA cells
#' @slot adj.ref.graph adjacency graph
#' @slot path.w weight (distance) of paths
#'
#' @export
#'
setClass("distref.raster",
         slots = c(r.ref.path = "character", r.ref.dim = "numeric",
                   r.ref.N = "numeric", r.ref.NA.nb = "numeric",
                   adj.ref.graph = "numeric", path.w = "numeric"),
         validity = function(object) {
           if (!length(object@r.ref.dim) %in% c(2,3))
             return("r.ref.dim should have 2 or 3 values")
           return(TRUE)
         }
)


#' Class for distref.data
#'
#' @slot data.pt.N numeric.
#' @slot NA.pos numeric
#' @slot NA.dist numeric
#'
#' @export
#'
setClass("distref.data",
         slots = c(data.pt.N = "numeric",
                   NA.pos = "numeric",
                   NA.dist = "numeric",
                   data = "ANY")
)

#' Transfrom reference raster as a network of paths for distance calculation
#'
#' @param r.ref Raster (1-layer) that will be used to calculate distances. Its resolution
#' impacts precision on distances.
#' @param step.ang Step for angles classes breaks (in degrees). This means that
#' a distance error is allowed in future calculations.
#' Distance error (\%) = 100 * (1 - cos((step.ang/2)*pi/180)).
#' With step.ang = 5 degrees, error < 0.1\%. Default set to 5.
#' @param filename Path where to save the dref.r issued from reference raster.
#' This may be usefull for re-using the raster for different datasets. Default
#' to paste0(tempdir(), "/dref.r.RData")
#' @param r.ref3D Logical. Wether the reference raster should be used as a 3D
#' surface to calculate distance accounting for elevation.
#' @param n.cores number of cores to run parallel calculation. Default to all
#' cores using parallel::detectCores()
#' @param save.rdists optional. Path to file where to save the distances matrix
#' between all non-NA raster cells. Saved as csv. If not provided, distances matrix
#' will not be calculated as it may be to big to be stored in memory.
#'
#' @details
#' Not all direct possible paths are calculated. Direct paths are retained in a
#' focal window around cells. The size of the focal window is determined by the
#' step.ang value.
#' All paths retained are stored in a undirected graph (as from library igraph).
#' All paths real distances are also stored to be weights for future points to
#' points distances calculations.
#'
#' @importFrom magrittr '%>%'
#'
#' @return
#' coords.r coordinates of all cells of the raster
#' r.ref.NA.nb cell that are NA in the reference raster (obstacles)
#' cell.Relpos.ngb row and col of neighbors retained after the angle calculation,
#'  relative to the cell of origin (moving window)
#' cell.Relpos.cross row and col of cells crossed by lines between origin and
#' neighbors, relative to the cell of origin (moving window)
#' @export

distref.raster <- function(r.ref, step.ang = 5,
                           filename = paste0(tempdir(), "/dref.r.RData"),
                           r.ref3D = FALSE, n.cores,
                           save.rdists)
{

  message(paste0("step.ang was set to ", step.ang,
                 " degrees. This conducts to potential error of ",
                 round(100 * (1 - cos((step.ang/2)*pi/180)), digits = 6),
                 " % of distances calculated."
          ))

  if (raster::nlayers(r.ref) != 1) {
    r.ref <- raster::raster(r.ref, 1)
    message("Only the first layer of the raster is used as reference.")
  }

  if (!raster::fromDisk(r.ref)) {
    raster::writeRaster(r.ref, filename = paste0(tempdir(), "/r.ref.grd"),
                overwrite = TRUE)
    r.ref <- raster::raster(paste0(tempdir(), "/r.ref.grd"))
    message(paste("r.ref has been stored in temporary directory:",
                  paste0(tempdir(), "/r.ref.grd")))
  }

  if (missing(n.cores)) {
    cl <- snow::makeCluster(parallel::detectCores())
  } else {
    cl <- snow::makeCluster(n.cores)
  }

  longlat <- raster::isLonLat(r.ref)
	coords.r <- sp::coordinates(r.ref)

	# Separate NA values in the reference raster (obstacles)
	r.ref.N <- c(1:raster::ncell(r.ref))
  r.ref.NA.nb <- which(is.na(raster::values(r.ref)))
	if (length(r.ref.NA.nb) != 0) {r.ref.N <- c(1:raster::ncell(r.ref))[-r.ref.NA.nb]}

	# Determine which cells will be retained for ngb of interest ==== ------------
	# = define the focal rectangle for distance calculations
	# Find all cells that have a different angle with x-axis
	# Define a step of angle that do not change distance a lot
	# Choose the central cell as origin
	pt.ref.nb <- raster::cellFromRowCol(
	  r.ref, row = round(nrow(r.ref) / 2),
	  col = round(ncol(r.ref) / 2))
	pt.ref <- coords.r[pt.ref.nb,]
	# angle between point and ref point
	angles <- atan2(y = (coords.r[,2] - pt.ref[2]),
	                x = (coords.r[,1] - pt.ref[1]))*360/(2*pi)
	# Only cells with angle in ]-90;90] are necessary to be kept
	# As if 1 is ngb with 2, so 2 is ngb with 1
	angles[which(angles <= (-90 + 0.5*step.ang) |
	               angles > (90 + 0.5*step.ang))] <- NA
	# Choose step for angles classes breaks
	# This means that an distance error is allowed
	brk <- seq(-90 + 0.5*step.ang, 90 + 0.5*step.ang, by = step.ang)
	ang.cut <- cut(angles, breaks = brk, labels = 1:(length(brk) - 1),
	               include.lowest = TRUE, right = TRUE)
	# fisrt angle is in the same class than last angle as -180=180
	ang.cut <- as.numeric(as.character(ang.cut))
	# distances between point and ref point
	distances <- raster::values(raster::distanceFromPoints(r.ref, pt.ref))
	distances[distances == 0] <- NA

	# Find the smallest grid around the pixels including the points of interest
	# and their position in the grid
	cell.which <- c(by(cbind(1:length(distances), distances), ang.cut,
	                   function(x) x[which(x[,2] == min(x[,2], na.rm = TRUE)),1]))
	# row and column number of cells
	RowCol.cells <- raster::rowColFromCell(r.ref, cell.which)
	RowCol.center <- raster::rowColFromCell(r.ref, pt.ref.nb)

	# equivalent matrix of weights (the moving window) for focal calculations
	# Cell positions relative to the focal cell
	cell.Relpos.ngb <- t(apply(RowCol.cells, 1, function(x) x - RowCol.center))

	# Build a little raster around the reference point
	# It will be used as moving window around each cell in future calculations
	block.r <- raster::raster(
	  raster::getValuesBlock(r.ref,
	                 row = min(cell.Relpos.ngb[,1] + RowCol.center[1]),
	                 nrows = 2*max(abs(cell.Relpos.ngb[,1])) + 1,
	                 col = min(cell.Relpos.ngb[,2] + RowCol.center[2]),
	                 ncols = 2*max(abs(cell.Relpos.ngb[,2])) + 1,
	                 format = 'matrix'))
	raster::extent(block.r) <- c(0.5, ncol(block.r) + 0.5, 0.5, nrow(block.r) + 0.5)
	# center of the grid
	block.w.center <- c(max(abs(cell.Relpos.ngb[,1])) + 1,
	                    max(abs(cell.Relpos.ngb[,2])) + 1)
	block.w.center.nb <- raster::cellFromRowCol(block.r, block.w.center[1],
	                                    block.w.center[2])

	# Find cell numbers that are crossed by direct line between 2 points
	ngb.i <- raster::cellFromRowCol(block.r,
	                        rownr = cell.Relpos.ngb[,1] + block.w.center[1],
	                        colnr = cell.Relpos.ngb[,2] + block.w.center[2])
	# Create lines
	L.i <- lapply(ngb.i, function(x) {
	  sp::Lines(list(sp::Line(list(rbind(sp::coordinates(block.r)[block.w.center.nb,],
	                             sp::coordinates(block.r)[x,])))), x)})
	ABline <- sp::SpatialLines(L.i)

	# Find cells crossed by Lines
	cellLines <- raster::cellFromLine(block.r, ABline)

	# row and column number of cells except original one (block.w.center.nb)
	RowCol.cells.i <- lapply(cellLines, function(x) {
	  raster::rowColFromCell(block.r, x[-which(x == block.w.center.nb)])
	  })

	# Standardize by cell positions relative to the focal cell
	cell.Relpos.cross <- lapply(RowCol.cells.i, function(x) {
	  t(apply(x, 1, function(x) x - block.w.center))
	  })

	# Calculate adjacency matrix from unique paths ==== --------------------------
	# Determine paths that do not cross an obstacle if any
  	adj.list <- snow::parLapply(cl, r.ref.N,
                          function(i, r.ref, coords.r, r.ref.NA.nb,
                                   cell.Relpos.ngb,
                                   cell.Relpos.cross) {
                            Find.ngb.wo.obstacle(i = i, r.ref, coords.r,
                                                  r.ref.NA.nb, cell.Relpos.ngb,
                                                  cell.Relpos.cross)
                          }, r.ref = r.ref, coords.r = coords.r,
                          r.ref.NA.nb = r.ref.NA.nb,
                          cell.Relpos.ngb = cell.Relpos.ngb,
                          cell.Relpos.cross = cell.Relpos.cross)
  	snow::stopCluster(cl)
  	gc()

	adj.list.null <- unlist(lapply(adj.list, function(x) is.null(x)))
	adj.ref <- t(as.data.frame(lapply(c(1:length(adj.list))[!adj.list.null],
	                                  function(x) t(adj.list[[x]]))))

	# Remove paths on same cell
	adj.ref.unique <- adj.ref[!apply(adj.ref, 1, function(x) (x[1] == x[2])),]

	# graph links as in package igraph
	adj.ref.graph <- igraph::graph.edgelist(adj.ref.unique, directed = FALSE)

	# Add vertices to the graph
	if (igraph::vcount(adj.ref.graph) <= raster::ncell(r.ref)) {
	  adj.ref.graph <- adj.ref.graph %>%
  	  igraph::add_vertices(raster::ncell(r.ref) - igraph::vcount(adj.ref.graph))
	}

	# Distance of each path to weight future routes
	if (!r.ref3D) {
	path.w <-
	  raster::pointDistance(
	    sp::coordinates(r.ref)[adj.ref.unique[,1],],
	    sp::coordinates(r.ref)[adj.ref.unique[,2],],
	    longlat = longlat)
	} else {
	 # use gam to smooth ngb focal window with 3d data
	  # divide each path from center cell to ngb into 1/2 raster res
	  # get coordinate, predict elevation
	  # calculte path distance with 3d dist
	  # be careful with longlat : calculate (pythagore)
	  # sqrt(spDist(x,y)^2 + (z-z0)^2)
	  # entre chaque point du sub-path

    # Use simple linear interpolation from package raster on the block raster with resample, method = "bilinear"

	}

	if (!missing(save.rdists)) {
  	r.ref.dists <- igraph::distances(adj.ref.graph,
                                     v = r.ref.N,
                                     to = r.ref.N,
                                     weights = path.w,
                                     mode = "all")
    dists.tbl <- dplyr::as.tbl(as.data.frame(r.ref.dists))
    if (raster::extension(save.rdists) != ".csv") {
      save.rdists <- paste0(save.rdists, ".csv")
    }
    readr::write_excel_csv(dists.tbl, save.rdists)
    warning("save.rdists only stores distances between non-obstacles cells")
	}

	res <- list(r.ref.path = r.ref@file@name, r.ref.dim = dim(r.ref),
	  r.ref.N = r.ref.N, r.ref.NA.nb = r.ref.NA.nb,
	  adj.ref.graph = adj.ref.graph, path.w = path.w)
	class(res) <- c("distref.raster")
	return(res)
}


#' Get characteristics of data regarding the reference raster
#'
#' @param data.pt Point data set. May be a 2-col matrix of (x,y) coordinates or a
#' SpatialPointDataFrame. Should be in same projection than r.ref
#' @param r.ref Reference raster
#' @param closest if a point is over an obstacle, set to TRUE to find the closest
#' non-obstacle r.ref cell
#' @param longlat Logical. if FALSE, Euclidean distance, if TRUE Great Circle
#' distance. Defined from r.ref.
#'
#' @export

distref.data <- function(data.pt, r.ref, closest = FALSE,
                         longlat = raster::isLonLat(r.ref))
{

  data.pt.N <- raster::extract(r.ref, data.pt, cellnumbers = TRUE)[,"cells"]

  if (sum(is.na(raster::values(r.ref)[data.pt.N])) > 0) {
    data.pt.N[which(is.na(raster::values(r.ref)[data.pt.N]))] <- NA
  }
  NA.pos <- which(is.na(data.pt.N))
	# For those that are NA (on obstacles due to r.ref resolution)
	# Find closest non-NA cell
  if (closest & sum(is.na(data.pt.N)) > 0) {
    if (length(NA.pos) != 1) {
      data.ptNA <- sp::coordinates(data.pt)[NA.pos,]
    } else {
      data.ptNA <- t(sp::coordinates(data.pt)[NA.pos,])
    }
    tabDist <-
      sp::spDists(data.ptNA,
              sp::coordinates(r.ref)[!is.na(raster::values(r.ref)),],
              longlat = longlat)
    w.close <- apply(tabDist, 1, which.min)
    NA.dist <- apply(t(1:nrow(tabDist)), 2, function(x) tabDist[x, w.close[x]])
    data.pt.N[NA.pos] <- which(!is.na(raster::values(r.ref)))[w.close]

    ## IF 3D, add message that distance is in 2D
  }

  if (!closest | length(NA.pos) ==  0) {
    res <- list(data.pt.N = data.pt.N, NA.pos = NA, NA.dist = NA, data = data.pt)
  } else {
    res <- list(data.pt.N = data.pt.N, NA.pos = NA.pos, NA.dist = NA.dist,
                data = data.pt)
  }
  class(res) <- c("distref.data")
  return(res)
}


#' Find neighbours without obstacles in different directions
#'
#' @param i cell number
#' @param coords.r coordinates of the reference raster
#' @param r.ref.NA.nb Numerotation of NA ref.raster cells
#' @param cell.Relpos.ngb Relative position of neighbours in the moving window
#' @param cell.Relpos.cross Numerotation of positions of neighbours of the moving window
#' @param r.ref reference raster
#'
#' @return matrix with all neighbours of each cell that are not NA
#' @export

Find.ngb.wo.obstacle <- function(i, r.ref, coords.r,
                                 r.ref.NA.nb, cell.Relpos.ngb,
                                 cell.Relpos.cross)
{

  mat.w.center <- raster::rowColFromCell(r.ref, i)
  ngb.i <-
    raster::cellFromRowCol(r.ref,
                   rownr = cell.Relpos.ngb[, 1] + mat.w.center[1],
                   colnr = cell.Relpos.ngb[, 2] + mat.w.center[2])

  if (length(r.ref.NA.nb) == 0) {
    # If no NA, all ngb are kept
    res <- cbind(i, ngb.i[!is.na(ngb.i)])
  } else {
    # Remove ngb which are NA (= obstacle)
    ngb.i[which(ngb.i %in% r.ref.NA.nb)] <- NA

    if (sum(is.na(ngb.i)) != length(ngb.i)) {
      # Get cell numbers of ngb (remove NA ngb)
      cell.Lines.Rel <-
        lapply(c(1:length(cell.Relpos.cross))[!is.na(ngb.i)], function(x)
          raster::cellFromRowCol(
            r.ref,
            cell.Relpos.cross[[x]][,1] + mat.w.center[1],
            cell.Relpos.cross[[x]][,2] + mat.w.center[2]
          ))

      # Only keep lines that do not cross an obstacle
      keepL <-
        which(unlist(lapply(cell.Lines.Rel, function(x) {
          (sum(x %in% r.ref.NA.nb) == 0)
        })))

      res <- cbind(i, ngb.i[!is.na(ngb.i)][keepL])
    } else {
      res <- NULL
    }
  }
  return(res)
}


#' Calculate distance from one point to a set of points using adjacency matrix
#'
#' @param dref.from object of class distref.data for a single point
#' @param dref.to object of class distref.data
#' @param dref.r object of class distref.raster corresponding to r.ref
#' @param r.ref reference raster (only used for its coordinates here)
#' @param longlat Logical. if FALSE, Euclidean distance, if TRUE Great Circle
#' distance. Defined from r.ref.
#' @param small.dist if TRUE, distances <= 1 paths are directly calculated.
#' @param keep.path Logical. Whether to keep SpatialLines of paths.
#' This allows real distances when 2 points are in the same raster cell or close.
#' @return Vector of distances of length(B)
#' @export
#'
#' @examples
#' \dontrun{
#' # Need gdistance ??
#' }

Pt2Pts.wo.obstacle <- function(dref.from, dref.to, dref.r,
                               r.ref, longlat,
                               small.dist = TRUE,
                               keep.path = FALSE)
{
  if (missing(longlat)) {longlat <- raster::isLonLat(r.ref)}

  distAtoB <- rep(-1, length(dref.to$data.pt.N))
  distAtoB.w <- 1:length(distAtoB)
  linesAtoB <- list()
  A <- dref.from$data.pt.N
  B <- dref.to$data.pt.N
  A.xy <- from.xy <- sp::coordinates(dref.from$data)
  B.xy <- to.xy <- sp::coordinates(dref.to$data)
  adj.ref.graph <- dref.r$adj.ref.graph
  path.w <- dref.r$path.w

  # Calculate distances when points are in the same raster cell
  dist0 <- which(B %in% A)
  if (length(dist0) != 0) {
    if (small.dist) {
      if (length(dist0) == 1) {
        distAtoB[dist0] <- sp::spDists(x = from.xy, y = t(to.xy[dist0,]), longlat)
      } else {
        distAtoB[dist0] <- sp::spDists(x = from.xy, y = to.xy[dist0,], longlat)
      }
    } else {
      distAtoB[dist0] <- 0
    }
    tmp <- lapply(1:length(dist0), function(x)
    {
      linesAtoB[[dist0[x]]] <<- sp::Lines(list(sp::Line(list(
        rbind(from.xy, to.xy[x,])
      ))), dist0[x])
    })
    rm(tmp)
    B <- B[-dist0]
    distAtoB.w <- distAtoB.w[-dist0]
    B.xy <- to.xy[-dist0,]
  }

  if (length(B) != 0) {
    shpath <- igraph::shortest_paths(adj.ref.graph, from = A, to = B,
                                     weights = path.w, mode = "all"
    )

    all.path <- 1:length(shpath$vpath)
    path.out <- integer(0)
    l.path <- lengths(shpath$vpath)

    # If path is only one Line, calculate real distance
    if (small.dist) {
      path.out <- which(l.path == 2)
      if (length(path.out) != 0) {
        if (length(path.out) == 1) {
          distAtoB[distAtoB.w][path.out] <-
            sp::spDists(x = A.xy, y = t(B.xy[path.out,]), longlat)
        } else {
          distAtoB[distAtoB.w][path.out] <-
            sp::spDists(x = A.xy, y = B.xy[path.out,], longlat)
        }
        tmp <- lapply(1:length(path.out), function(x)
        {
          linesAtoB[[distAtoB.w[path.out[x]]]] <<- sp::Lines(
            list(sp::Line(list(rbind(A.xy, B.xy[path.out[x],])
            ))), distAtoB.w[path.out[x]])
        })
        rm(tmp)
        all.path <- all.path[-path.out]
      }
    }

    path.out.0 <- all.path[which(l.path[all.path] == 0)]
    if (length(path.out.0) != 0) {
      warning("Couldn't reach some vertices, some distances were set to NA")
      if (keep.path) {
        warning(paste0("Hence, path do not correspond to distances.",
                       "Use ID of SpatialLines to find corresponding path"))
      }
      distAtoB[distAtoB.w][path.out.0] <- NA
      path.out <- c(path.out, path.out.0)
      all.path <- all.path[-which(all.path %in% path.out.0)]
      tmp <- lapply(path.out.0, function(x) {
        linesAtoB[[distAtoB.w[x]]] <<- NULL})
      rm(tmp)
    }

    if (length(all.path) != 0) {
      allLines <- lapply(all.path, function(x) {
      # Remove first and last path r.ref point to replace by location
        pathmin <- sp::coordinates(r.ref)[utils::head(shpath$vpath[[x]][-1],-1),]
        if (is.null(nrow(pathmin))) {
          pathmin <- data.frame(x = pathmin[1], y = pathmin[2])
        }
        if (nrow(pathmin) != 0) {
        coords <- rbind(A.xy,
                        pathmin,
                        B.xy[x,]
       )
        } else {
        coords <- rbind(A.xy, B.xy[x,])
        }
        rownames(coords) <- NULL
        linesAtoB[[distAtoB.w[x]]] <<- res <- sp::Lines(list(sp::Line(list(
          coords
        ))), distAtoB.w[x])
        res
      })

      AtoB <- sp::SpatialLines(allLines)
      distAtoB.tmp <- sp::SpatialLinesLengths(AtoB, longlat = longlat)

      if (length(path.out) == 0) {
        distAtoB[distAtoB.w] <- distAtoB.tmp
      } else {
        distAtoB[distAtoB.w][-path.out] <- distAtoB.tmp
      }
    }
    if (length(path.out.0) == 0) {
      LinesAB <- sp::SpatialLinesDataFrame(
        sp::SpatialLines(linesAtoB),
        data = data.frame(ID = 1:length(distAtoB)), match.ID = FALSE)
    } else {
      LinesAB <- sp::SpatialLinesDataFrame(
        sp::SpatialLines(linesAtoB),
        data = data.frame(ID = 1:length(distAtoB)[-distAtoB.w[path.out.0]]),
        match.ID = FALSE)
    }
    raster::projection(LinesAB) <- raster::projection(r.ref)
    # dist <- SpatialLinesLengths(LinesAB)
    if (!keep.path) {LinesAB <- NA}
    res <- list(Lines = LinesAB, dist = distAtoB)
  } else {
    res <- list(Lines = NA, dist = distAtoB)
  }
  return(res)
}


#' Calculate 2D distances from points to points avoiding obstacles
#'
#'
#' @param from Point data set from which to calculate distances. May be a 2-col
#' matrix of (x,y) coordinates or a SpatialPointDataFrame. Should be in same
#' projection than r.ref. This can also be an object of class distref.data.
#' @param to Point data set from which to calculate distances. May be a 2-col
#' matrix of (x,y) coordinates or a SpatialPointDataFrame. Should be in same
#' projection than r.ref. If not set to = from. This can also be an object of
#' class distref.data. If "to" is provided, it should be the bigger dataset.
#' @param r.ref reference raster used for distance without obstacles calculation
#' @param dref.r object of class distref.raster. Reference raster already
#' processed using distref.raster (e.g. usefull if previously saved using filename.)
#' @param n.cores number of cores to run parallel calculation. Default to all
#' cores using parallel::detectCores()
#' @param closest if a point is over an obstacle, set to TRUE to find the closest
#' non-obstacle r.ref cell
#' @param step.ang Step for angles classes breaks (in degrees). This means that
#' a distance error is allowed in future calculations.
#' Distance error (\%) = 100 * (1 - cos((step.ang/2)*pi/180)).
#' With step.ang = 5 degrees, error < 0.1\%. Default set to 5.
#' @param r.ref3D Logical. Wether the reference raster should be used as a 3D
#' surface to calculate distance accounting for elevation.
#' @param filename Path where to save the dref.r issued from reference raster.
#' This may be usefull for re-using the raster for different datasets. Default
#' to paste0(tempdir(),"/dref.r.RData")
#' @param small.dist if TRUE, distances <= 1 paths are directly calculated.
#' @param keep.path Logical. Whether to keep SpatialLines of paths.
#' This allows real distances when 2 points are in the same raster cell or close.
#' @param longlat Logical. if FALSE, Euclidean distance, if TRUE Great Circle
#' distance. Defined from r.ref.
#' @param tol numeric Tolerance to define distances between points to be null.
#' @param igraph Logical. Wether to calculate all distances through igraph (TRUE)
#' or use \code{\link{Pt2Pts.wo.obstacle}} to use igraph only when necessary.
#' Using igraph maybe less precise but a little quicker. Default to FALSE.
#'
#' @importFrom magrittr '%>%'
#'
#' @return
#' If keep.path is false, returns a matrix of distances with rows = "from"
#' and cols = "to".
#' If keep.path is true, returns a list where dist.mat is the matrix of distance
#' and listLines is a list of SpatialLines. Each list is a SpatialLinesDataFrame
#' corresponding to a starting point in "from".
#' @export
#'
#' @examples
#' \dontrun{
#' # For "from" and "to" being SpatialPoints
#' allpath <- dist.obstacle(from, to, r.ref, keep.path = TRUE)
#' plot(from)
#' points(to)
#' # All path from "from[1,]" to all "to"
#' lines(allpath$allLines[[1]])
#' }
dist.obstacle <- function(from, to, r.ref, dref.r,
                          n.cores, closest = FALSE,
                          step.ang = 5, r.ref3D = FALSE,
                          filename = paste0(tempdir(), "/dref.r.RData"),
                          small.dist = TRUE, keep.path = FALSE,
                          longlat = raster::isLonLat(r.ref), tol = 1e-5,
                          igraph = FALSE) {
  if (missing(n.cores)) {
    cl <- snow::makeCluster(parallel::detectCores())
  } else {
    cl <- snow::makeCluster(n.cores)
  }

  # r.ref preparation
  if (missing(dref.r)) {
    dref.r <- distref.raster(r.ref, step.ang,
      filename, r.ref3D, n.cores)
  }
  if (!inherits(dref.r, "distref.raster")) {
    stop("dref.raster must be of class distref.raster")
  }

  # data preparation
  if (inherits(from, "distref.data")) {
    dref.from <- from
  } else {
    dref.from <- distref.data(from, r.ref, closest)
  }
  if (missing(to)) {
    dref.to <- dref.from
  } else {
    if (inherits(to, "distref.data")) {
      dref.to <- to
    } else {
      dref.to <- distref.data(to, r.ref, closest)
    }
  }
  if (sum(is.na(dref.from$data.pt.N)) != 0 |
      sum(is.na(dref.to$data.pt.N)) != 0) {
   stop(paste(
     "Some of the stations are not over the raster or are over an obstacle,",
     "distances can not be calculated using 'r.ref' as a reference:\n",
     "In 'from':", paste(which(is.na(dref.from$data.pt.N)), collapse = ","), "\n",
     "In 'to':", paste(which(is.na(dref.to$data.pt.N)), collapse = ","), ".\n",
     "Remove these stations or use parameter 'closest = TRUE'.\n",
     "You can also previously run function 'distref.data' on your data to get",
     "distances from the closest non-obstacle 'r.ref' cell"
     ))
  }

  ## All as igraph -------------------------------------------------------------
  if (igraph) {
    adj.ref.graph <- dref.r$adj.ref.graph
    path.w <- dref.r$path.w
    equal.from.to <- FALSE
    if (nrow(dref.from$data) == nrow(dref.to$data)) {
      if (sum(dref.from$data != dref.to$data) == 0) {equal.from.to <- TRUE}
    } else {
      # Verify if "from" have common coords with "to"
      from.To.to <- sp::spDists(x = as.matrix(dref.from$data),
                                y = as.matrix(dref.to$data),
                                longlat = raster::isLonLat(r.ref))
      from.to.eq <- apply(from.To.to, 1, function(x) {which(x <= tol)})
      if (!is.list(from.to.eq)) {from.to.eq <- as.list(from.to.eq)}
      rm(from.To.to)
    }
    # Remove duplicates + message

    add.v <- 0
    adj.v.new <- add.coords <- NULL
    n.edge.r <- length(igraph::V(adj.ref.graph)[[]])
    for (i in c("from", "to")) {
      if (i == "from" | !equal.from.to) {
        fromto <- get(paste0("dref.", i))
        N.to.test <- 1:length(fromto$data.pt.N)
        fromto$vertices <- fromto$data.pt.N
        if (i == "to" & sum(lengths(from.to.eq)) != 0) {
          # Do not recalculate those in common with from
          tmp <- lapply(which(lengths(from.to.eq) != 0), function(x) {
            fromto$vertices[x] <<- rep(dref.from$vertices[x], length(x))
          })
          rm(tmp)
          N.to.test <- N.to.test[-which(lengths(from.to.eq) != 0)]
        }
        if (length(N.to.test) != 0) {
          # Do not recalculate if points are the same than r.ref positions
          fromto.To.r <- apply(t(N.to.test), 2, function(x) {
            sp::spDistsN1(pts = as.matrix(fromto$data[x,]),
                        pt = sp::coordinates(r.ref)[fromto$data.pt.N[x],],
                        longlat = raster::isLonLat(r.ref))
          })
          w.diff <- N.to.test[which(fromto.To.r >= tol)]
          if (length(w.diff) != 0) {
            # Add new points in the graph if not exists
            fromto$vertices[w.diff] <- (n.edge.r + 1):(n.edge.r + length(w.diff))
            add.v <- add.v + length(w.diff)
            n.edge.r <- n.edge.r + length(w.diff)
            add.coords.tmp <- fromto$data[w.diff,]
            if (is.null(add.coords)) {
              add.coords <- add.coords.tmp
            } else {
              add.coords <- rbind(add.coords, add.coords.tmp)
            }
            rm(add.coords.tmp)
            # Define new edges based on r.ref edges
            adj.v <- igraph::adjacent_vertices(adj.ref.graph,
                                       v = fromto$data.pt.N[w.diff],
                                       mode = "all")
            adj.v.new.tmp <- as.data.frame(lapply(1:length(w.diff), function(x) {
              rbind(w.diff[x], fromto$vertices[w.diff[x]], adj.v[[x]])
            }))
            if (is.null(adj.v.new)) {
              adj.v.new <- adj.v.new.tmp
            } else {
              adj.v.new <- cbind(adj.v.new, adj.v.new.tmp)
            }
            rm(adj.v.new.tmp)
          }
        }
        assign(paste0("dref.", i), fromto)
      } else {
        dref.to$vertices <- dref.from$vertices
      }
    }
    if (add.v != 0) {
      # Add vertices to the graph
      adj.ref.graph <- adj.ref.graph %>%
        igraph::add_vertices(add.v)
      # Add edges to graph
      adj.ref.graph <- adj.ref.graph %>%
        igraph::add_edges(unlist(adj.v.new[2:3,], use.names = FALSE))
      # Add distances to dists
      if (!r.ref3D) {
        path.w <- c(path.w, raster::pointDistance(
          p1 = as.matrix(fromto$data[unlist(adj.v.new[1,]),]),
          p2 = sp::coordinates(r.ref)[unlist(adj.v.new[3,]),],
          lonlat = raster::isLonLat(r.ref)))
      }
    }

    if (keep.path) {
      shpath <- igraph::all_shortest_paths(
        adj.ref.graph,
        from = dref.from$vertices, to = dref.to$vertices,
        weights = path.w, mode = "all")
      ## Verify that all paths from all "from" to all "to"
      # Add duplicates here in all.path
      if (is.null(shpath$vpath)) {shpath$vpath <- shpath$res}
      all.path <- 1:length(shpath$vpath)

      allLines <- lapply(all.path, function(x) {
        sp::Lines(list(sp::Line(list(
          rbind(sp::coordinates(r.ref), add.coords)[shpath$vpath[[x]],]
        ))), x)
      })

      AtoB <- sp::SpatialLines(allLines)
      raster::projection(AtoB) <- raster::projection(r.ref)
      # distAtoB.tmp <- sp::SpatialLinesLengths(AtoB, longlat = longlat)
    } else {
      dist.mat <- igraph::distances(
        adj.ref.graph,
        v = dref.from$vertices, to = dref.to$vertices,
        weights = path.w, mode = "all")
    }
  ## End all as igraph ---------------------------------------------------------
  } else {
    if (length(dref.from$data.pt.N) == 1) {
      res <- Pt2Pts.wo.obstacle(
        dref.from, dref.to, dref.r, r.ref,
        longlat = raster::isLonLat(r.ref),
        small.dist = small.dist,
        keep.path = keep.path)
      if (keep.path) {allLines <- res$Lines}
      dist.mat <- res$dist
      rm(res)
    } else {
      # stop if NA in data.pt.N because no closest
      # Reduce calculations if from = to: Symetric matrix
      if (keep.path) {allLines <- list()}
	    dist.mat <- snow::parApply(cl, t(1:length(dref.from$data.pt.N)), 2, function(
	      i,
	      dref.from, dref.to, dref.r, r.ref,
	      longlat, small.dist, keep.path)
	      {
	      dref.from <- list(data.pt.N = dref.from$data.pt.N[i],
	                        data = t(dref.from$data[i,]))
	      class(dref.from) <- "distref.data"
	      res <- Pt2Pts.wo.obstacle(
	        dref.from, dref.to, dref.r, r.ref,
	        longlat = raster::isLonLat(r.ref), small.dist = small.dist,
	        keep.path = keep.path)
	      if (keep.path) {
	        allLines[[i]] <<- res$Lines
	      }
	      return(res$dist)
	    }, dref.from = dref.from, dref.to = dref.to,
	    dref.r = dref.r, r.ref = r.ref,
	    longlat = longlat, small.dist = small.dist,
	    keep.path = keep.path)
	  snow::stopCluster(cl);gc()
    }
  }
  if (keep.path) {
    return(list(dist.mat = dist.mat, listLines = allLines))
  } else {
	  return(dist.mat)
  }
}
