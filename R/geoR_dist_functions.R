#' Inverse distance calculation using custom distances
#'
#' @param data vector of values to be interpolated at data coordinates
#' @param coords coordinates of observation points. Can also be an object of class
#' SpatialPoints*. Not necessary with loc.dist.
#' @param locations coordinates of points where to provide predictions. Can also be an object of class
#' SpatialPoints*. Not necessary with loc.dist.
#' @param loc.dist Matrix of distances between data (rows) and locations (cols)
#' @param nmax the number of nearest observations that should be used for a kriging prediction or simulation, where nearest is defined in terms of the space of the spatial locations. By default, all observations are used
#' @param nmin if the number of nearest observations within distance maxdist is less than nmin, a missing value will be generated; see maxdist
#' @param longlat Logical Wether coordinates are not-projected (TRUE) or are projected (FALSE).
#' @param max.dist Maximum distance radius for neighbours search in the idw interpolation.
#' maxdist is in km when longlat is TRUE (non-projected crs), in meters otherwise.
#' If loc.dist is specified, maxdist should be in same distance unit.
#' @param idp specify the inverse distance weighting power. Default to 2.
#'
#' @importFrom methods is
#'
#' @export

idw.dist <- function(data, coords, locations, loc.dist, idp = 2, max.dist,
  nmin = 1, nmax, longlat = TRUE) {

  ## Allow for user calculated distance matrix
  if (missing(loc.dist)) {
    if (missing(locations)) {
      locations <- coords
    }
    if (sum(grepl("SpatialPoints", is(locations))) >= 1) {
      locations <- sp::coordinates(locations)
      sp::proj4string(locations) <- sp::proj4string(coords)
    }
    # v0 <- geoR::loccoords(coords = coords, locations = locations)
    v0 <- sp::spDists(x = coords, y = locations, longlat = longlat)
  } else {
    v0 <- loc.dist
  }

  new <- TRUE
  if (new) {
    loc.dist.prep <- Prep_loc_dist(loc.dist = v0,
      nmin = nmin, nmax = nmax, max.dist = max.dist)

    idw <- idw.dist.special(data = data,
                     loc.dist.order = loc.dist.prep$loc.dist.order,
                     loc.dist.nmax = loc.dist.prep$loc.dist.nmax,
                     idp = idp)
  } else
    {
    if (missing(max.dist)) {max.dist <- NULL}

    if (!is.null(max.dist)) {
      v0[which(v0 > max.dist)] <- NA
    }
    if (missing(nmax)) {nmax <- NULL}

    calc_idw <- function(x, data, idp, nmin, nmax) {
      if (!is.null(nmax)) {
        ord.x <- order(x, na.last = TRUE)
        x <- x[ord.x]
        data <- data[ord.x]
        x[-c(1:nmax)] <- NA
        data[-c(1:nmax)] <- NA
      }
      w.noNA <- which(!is.na(data) & !is.na(x))
      if (length(w.noNA) < nmin) {
        res <- NA
      } else {
        res <- sum(data[w.noNA] * 1/(x[w.noNA]^idp), na.rm = TRUE) / sum(1/x[w.noNA]^idp, na.rm = TRUE)
      }
      return(res)
    }
      idw <- apply(v0, 2, function(x, data, idp, nmin, nmax)
                    calc_idw(x, data, idp, nmin, nmax),
                   data = data, idp = idp, nmin = nmin, nmax = nmax)
  }
  # names(idw) <- "idw"
  return(idw)
}

#' Prepare loc.dist for idw special case
#'
#' @param loc.dist Matrix of distances between data (rows) and locations (cols)
#' @param nmax the number of nearest observations that should be used for a kriging prediction or simulation, where nearest is defined in terms of the space of the spatial locations. By default, all observations are used
#' @param nmin if the number of nearest observations within distance maxdist is less than nmin, a missing value will be generated; see maxdist
#' @param max.dist Maximum distance radius for neighbours search in the idw interpolation.
#'
#' @importFrom magrittr '%>%'
#' @importFrom rlang .data
#'
#' @export

Prep_loc_dist <- function(loc.dist, nmin = 1, nmax, max.dist)
{
  if (!dplyr::is.tbl(loc.dist)) {
    loc.dist.tbl <- dplyr::as.tbl(as.data.frame(loc.dist))
  } else {
    loc.dist.tbl <- loc.dist
  }
  # NA data over max.dist
  if (missing(max.dist)) {max.dist <- NULL}
  if (!is.null(max.dist)) {
    loc.dist.tbl[loc.dist.tbl >= max.dist] <- NA
  }
  # All data if no nmax
  if (missing(nmax)) {nmax <- nrow(loc.dist.tbl)}
  # Find lines data necessary according to nmax
  loc.dist.order <- loc.dist.tbl %>%
    dplyr::mutate_all(order) %>%
    dplyr::slice(1:nmax)

  # Get distances lines necessary according to nmax
  # and corrected by nmin
  loc.dist.nmax <- loc.dist.tbl %>%
    dplyr::mutate_all(.data, function(x) sort(x, na.last = TRUE)) %>%
    dplyr::slice(1:nmax) %>%
    dplyr::mutate_all(.data, function(x) {
      if (sum(!is.na(x)) < nmin) {rep(NA, length(x))}else{x}})

  return(list(loc.dist.order = loc.dist.order, loc.dist.nmax = loc.dist.nmax))
}

#' Special case idw for reduced time calculation
#'
#' @param data vector of data to interpolate on points. Refer to output of Prep_loc_dist
#' @param loc.dist.order matrix of nmax lines and nb points columns with position
#' of data to keep
#' @param loc.dist.nmax matrix of nmax lines and nb points columns with distances
#' between kept data and points position
#' @param idp specify the inverse distance weighting power. Default to 2.
#'
#' @export
idw.dist.special <- function(data, loc.dist.order, loc.dist.nmax, idp)
{
  # Couple with data (Do for each layers)
  data.mat <- matrix(data[unlist(loc.dist.order)], ncol = ncol(loc.dist.order))

  # Calculate idw: 5sec
  idw <- apply(t(1:ncol(data.mat)), 2,
               function(x) { sum(data.mat[,x] * 1/(loc.dist.nmax[,x]^idp), na.rm = TRUE) /
                   sum(1/(loc.dist.nmax[,x]^idp), na.rm = TRUE)})
  return(idw)
}

#' Variogram calculation. Modified from geoR.
#'
#' @inheritParams geoR::variog
#' @param dist.mat Square matrix of distances between points of the dataset
#'
#' @import geoR
#'
#' @export
#'
variog.dist <- function(geodata, coords = geodata$coords, data = geodata$data, uvec = "default",
  breaks = "default", trend = "cte", lambda = 1, option = c("bin", "cloud", "smooth"), estimator.type = c("classical",
    "modulus"), nugget.tolerance, max.dist, dist.mat, pairs.min = 2, bin.cloud = FALSE,
  direction = "omnidirectional", tolerance = pi/8, unit.angle = c("radians", "degrees"),
  angles = FALSE, messages, ...)
  {
  ## If dist.mat, bin.cloud <- 'diff'
  if (!missing(dist.mat)) {
    bin.cloud <- "diff"
  }

  if (missing(geodata))
    geodata <- list(coords = coords, data = data)
  call.fc <- match.call()
  if (missing(messages))
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages"))) else messages.screen <- messages
  keep <- list(...)
  if (is.null(keep$keep.NA))
    keep.NA <- FALSE else keep.NA <- keep$keep.NA
  unit.angle <- match.arg(unit.angle)
  if (mode(direction) == "numeric") {
    if (length(direction) > 1)
      stop("only one direction is allowed")
    if (length(tolerance) > 1)
      stop("only one tolerance value is allowed")
    if (unit.angle == "degrees") {
      ang.deg <- direction
      ang.rad <- (ang.deg * pi)/180
      tol.deg <- tolerance
      tol.rad <- (tol.deg * pi)/180
    } else {
      ang.rad <- direction
      ang.deg <- (ang.rad * 180)/pi
      tol.rad <- tolerance
      tol.deg <- (tol.rad * 180)/pi
    }
    if (ang.rad > pi | ang.rad < 0)
      stop("direction must be an angle in the interval [0,pi[ radians")
    if (tol.rad > pi/2 | tol.rad < 0)
      stop("tolerance must be an angle in the interval [0,pi/2] radians")
    if (tol.deg >= 90) {
      direction <- "omnidirectional"
      cat("variog: computing omnidirectional variogram\n")
    } else {
      if (messages.screen) {
        cat(paste("variog: computing variogram for direction = ", round(ang.deg, digits = 3),
          " degrees (", round(ang.rad, digits = 3), " radians)\n", sep = ""))
        cat(paste("        tolerance angle = ", round(tol.deg, digits = 3), " degrees (",
          round(tol.rad, digits = 3), " radians)\n", sep = ""))
      }
    }
  } else if (messages.screen)
    cat("variog: computing omnidirectional variogram\n")
  coords <- as.matrix(coords)
  data <- as.matrix(data)
  if (nrow(coords) != nrow(data))
    stop("coords and data have incompatible dimensions")
  data.var <- apply(data, 2, stats::var)
  n.data <- nrow(coords)
  n.datasets <- ncol(data)
  data <- drop(data)
  option <- match.arg(option)
  estimator.type <- match.arg(estimator.type)
  if (abs(lambda - 1) > 1e-04) {
    if (abs(lambda) < 1e-04)
      data <- log(data) else data <- ((data^lambda) - 1)/lambda
  }
  xmat <- unclass(geoR::trend.spatial(trend = trend, geodata = geodata))
  if (nrow(xmat) != n.data)
    stop("coords and trend have incompatible sizes")
  if (trend != "cte") {
    if (is.vector(data)) {
      temp.fit <- stats::lm(data ~ xmat + 0)
      beta.ols <- temp.fit$coeff
      data <- temp.fit$residuals
      temp.fit <- NULL
      names(data) <- NULL
    } else {
      only.res <- function(y, x) stats::lm(y ~ xmat + 0)$residuals
      data <- apply(data, 2, only.res, x = xmat)
      only.beta <- function(y, x) stats::lm(y ~ xmat + 0)$coef
      beta.ols <- apply(data, 2, only.beta, x = xmat)
    }
  } else beta.ols <- colMeans(as.matrix(data))
  ## Allow to add your own distance matrix
  if (missing(dist.mat)) {
    u <- as.vector(stats::dist(as.matrix(coords)))
  } else {
    dist.vec <- dist.mat[which(col(dist.mat) < row(dist.mat))]
    u <- dist.vec
  }
  if (missing(nugget.tolerance) || nugget.tolerance < 1e-11) {
    nugget.tolerance <- 1e-12
    nt.ind <- FALSE
  } else {
    if (mode(nugget.tolerance) != "numeric")
      stop("nugget.tolerance must be numeric")
    nt.ind <- TRUE
  }
  min.dist <- min(u)
  if (min.dist < nugget.tolerance)
    nt.ind <- TRUE
  if (direction != "omnidirectional" | angles) {
    u.ang <- .C("tgangle", as.double(as.vector(coords[, 1])), as.double(as.vector(coords[,
      2])), as.integer(dim(coords)[1]), res = as.double(rep(0, length(u))), PACKAGE = "geoR")$res
    if (any(is.na(u.ang)))
      stop("NA returned in angle calculations maybe due to co-located data")
    u.ang <- atan(u.ang)
    u.ang[u.ang < 0] <- u.ang[u.ang < 0] + pi
  }
  if (option == "bin" && bin.cloud == FALSE && direction == "omnidirectional") {
    if (missing(max.dist))
      umax <- max(u) else umax <- max(u[u < max.dist])
    dbins <- geoR:::.define.bins(max.dist = umax, uvec = uvec, breaks = breaks, nugget.tolerance = nugget.tolerance)
    uvec <- dbins$uvec
    bins.lim <- dbins$bins.lim
    nbins <- length(bins.lim) - 1
    if (missing(max.dist))
      max.dist <- max(bins.lim)
    if (bins.lim[1] < 1e-16)
      bins.lim[1] <- -1
    bin.f <- function(data) {
      cbin <- vbin <- sdbin <- rep(0, nbins)
      .C("binit", as.integer(n.data), as.double(as.vector(coords[, 1])), as.double(as.vector(coords[,
        2])), as.double(as.vector(data)), as.integer(nbins), as.double(as.vector(bins.lim)),
        as.integer(estimator.type == "modulus"), as.double(max.dist), cbin = as.integer(cbin),
        vbin = as.double(vbin), as.integer(TRUE), sdbin = as.double(sdbin), PACKAGE = "geoR")[c("vbin",
        "cbin", "sdbin")]
    }
    result <- array(unlist(lapply(as.data.frame(data), bin.f)), dim = c(nbins, 3, n.datasets))
    indp <- (result[, 2, 1] >= pairs.min)
    result[!indp, 1, ] <- NA
    if (bins.lim[1] < 0)
      bins.lim[1] <- 0
    if (!nt.ind) {
      uvec <- uvec[-1]
      indp <- indp[-1]
      bins.lim <- bins.lim[-1]
      result <- result[-1, , , drop = FALSE]
    }
    if (keep.NA)
      result <- list(u = uvec, v = result[, 1, ], n = result[, 2, 1], sd = result[, 3,
        ], bins.lim = bins.lim, ind.bin = indp) else result <- list(u = uvec[indp], v = result[indp, 1, ], n = result[indp, 2, 1],
      sd = result[indp, 3, ], bins.lim = bins.lim, ind.bin = indp)
  } else {
    data <- as.matrix(data)
    v <- matrix(0, nrow = length(u), ncol = n.datasets)
    for (i in 1:n.datasets) {
      v[, i] <- as.vector(stats::dist(data[, i]))
      if (estimator.type == "modulus")
        v[, i] <- v[, i, drop = FALSE]^(0.5) else v[, i] <- (v[, i, drop = FALSE]^2)/2
    }
    if (!missing(max.dist)) {
      v <- v[u <= max.dist, , drop = FALSE]
      if (direction != "omnidirectional")
        u.ang <- u.ang[u <= max.dist]
      u <- u[u <= max.dist]
    }
    if (direction != "omnidirectional") {
      ang.lower <- ang.rad - tol.rad
      ang.upper <- ang.rad + tol.rad
      if (ang.lower >= 0 & ang.upper < pi)
        ang.ind <- (!is.na(u.ang) & ((u.ang >= ang.lower) & (u.ang <= ang.upper)))
      if (ang.lower < 0)
        ang.ind <- (!is.na(u.ang) & ((u.ang < ang.upper) | (u.ang > (pi + ang.lower))))
      if (ang.upper >= pi)
        ang.ind <- (!is.na(u.ang) & ((u.ang > ang.lower) | (u.ang < (ang.upper - pi))))
      v <- v[ang.ind, , drop = FALSE]
      u <- u[ang.ind]
    }
    data <- drop(data)
    v <- drop(v)
    if (option == "cloud") {
      result <- list(u = u, v = v)
      if (angles)
        result$angles <- u.ang
    }
    if (option == "bin") {
      if (missing(max.dist))
        umax <- max(u) else umax <- max(u[u < max.dist])
      if (bin.cloud == "diff")
        dd <- geoR::diffpairs(coords, data)$diff else dd <- 0
      result <- geoR:::.rfm.bin(cloud = list(u = u, v = v, d = dd), estimator.type = estimator.type,
        uvec = uvec, breaks = breaks, nugget.tolerance = nugget.tolerance, bin.cloud = bin.cloud,
        max.dist = umax, keep.NA = keep.NA)
      if (keep.NA) {
        if (pairs.min > 0) {
          indp <- (result$n < pairs.min)
          if (!nt.ind) {
          for (i in 1:5) result[[i]] <- result[[i]][-1]
          indp <- indp[-1]
          }
          if (is.matrix(result$v)) {
          result$v[indp, ] <- result$sd[indp, ] <- NA
          } else {
          result$v[indp] <- result$sd[indp] <- NA
          }
        }
        result$ind.bin <- indp
      } else {
        if (pairs.min > 0) {
          if (!nt.ind) {
          for (i in 1:5) result[[i]] <- result[[i]][-1]
          }
          indp <- (result$n >= pairs.min)
          if (is.matrix(result$v)) {
          result$v <- result$v[indp, ]
          result$sd <- result$sd[indp, ]
          } else {
          result$v <- result$v[indp]
          result$sd <- result$sd[indp]
          }
          result$u <- result$u[indp]
          result$n <- result$n[indp]
        }
        result$ind.bin <- indp
      }
    }
    if (option == "smooth") {
      if (is.matrix(v))
        stop("smooth not yet available for more than one data-set")
      temp <- stats::ksmooth(u, v, ...)
      result <- list(u = temp[[1]], v = temp[[2]])
    }
    if (missing(max.dist))
      max.dist <- max(u)
  }
  if (nt.ind) {
    if (!exists(".variog4.nomessage", where = 1))
      cat("variog: co-locatted data found, adding one bin at the origin\n")
    if (all(result$u[1:2] < 1e-11))
      result$u[2] <- sum(result$bins.lim[2:3])/2
  }
  result <- c(result, list(var.mark = data.var, beta.ols = beta.ols, output.type = option,
    max.dist = max.dist, estimator.type = estimator.type, n.data = n.data, lambda = lambda,
    trend = trend, pairs.min = pairs.min))
  result$nugget.tolerance <- nugget.tolerance
  if (direction != "omnidirectional")
    result$direction <- ang.rad else result$direction <- "omnidirectional"
  if (direction != "omnidirectional")
    result$tolerance <- tol.rad else result$tolerance <- "none"
  result$uvec <- uvec
  result$call <- call.fc
  oldClass(result) <- "variogram"
  return(result)
}


#' Fit of gaussian field. Modified from function likfit in geoR.
#'
#' @inheritParams geoR::likfit
#' @param dist.mat Square matrix of distances between data points
#' @export

likfit.dist <- function(geodata, coords = geodata$coords, data = geodata$data, trend = "cte",
  ini.cov.pars, fix.nugget = FALSE, nugget = 0, fix.kappa = TRUE, kappa = 0.5, fix.lambda = TRUE,
  lambda = 1, fix.psiA = TRUE, psiA = 0, fix.psiR = TRUE, psiR = 1, cov.model, dist.mat,
  realisations, lik.method = "ML", components = TRUE, nospatial = TRUE, limits = pars.limits(),
  print.pars = FALSE, messages, ...)
  {
  name.geodata <- deparse(substitute(geodata))
  call.fc <- match.call()
  ldots <- list(...)
  temp.list <- list()
  temp.list$print.pars <- print.pars
  if (missing(messages))
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages"))) else messages.screen <- messages
  if (!missing(ini.cov.pars)) {
    if (any(class(ini.cov.pars) == "eyefit")) {
      ini.cov.pars <- ini.cov.pars[[1]]
    }
    if (any(class(ini.cov.pars) == "variomodel")) {
      cov.model <- ini.cov.pars$cov.model
      kappa <- ini.cov.pars$kappa
    }
  }
  if (missing(cov.model))
    cov.model <- "matern"
  cov.model <- match.arg(cov.model, choices = geoR:::.geoR.cov.models)
  if (cov.model == "stable")
    cov.model <- "powered.exponential"
  if (any(cov.model == c("power", "gneiting.matern", "gencauchy")))
    stop(paste("parameter estimation for", cov.model, "is not yet implemented"))
  fixed.pars <- list(cov.model = cov.model)
  if (fix.nugget)
    fixed.pars$nugget <- nugget
  if (fix.kappa)
    fixed.pars$kappa <- kappa
  if (fix.psiA)
    fixed.pars$psiA <- psiA
  if (fix.psiR)
    fixed.pars$psiR <- psiR
  geoR:::.check.geoRparameters.values(list = fixed.pars, messages = messages.screen)
  if (cov.model == "matern" & all(kappa == 0.5))
    cov.model <- "exponential"
  temp.list$cov.model <- cov.model
  if (cov.model == "powered.exponential")
    if (limits$kappa["upper"] > 2)
      limits$kappa["upper"] <- 2
  if (cov.model == "gencauchy")
    if (limits$kappa2["upper"] > 2)
      limits$kappa2["upper"] <- 2
  lik.MET <- c("ML", "ml", "RML", "REML", "rml", "reml")
  MET <- pmatch(names(ldots), "method") == 1
  if (!is.na(MET) && any(MET) && (ldots[[which(MET)]] %in% lik.MET)) {
    warning("argument \"method\" has changed and is now used as an argument to be passed to optim(). Use \"lik.method\" to define the likelihood method")
    lik.method <- lik.MET[pmatch(ldots[[which(MET)]], lik.MET)]
    ldots[which(as.logical(pmatch(names(ldots), "method", nomatch = 0)))] <- NULL
  }
  method.lik <- lik.method
  if (method.lik %in% c("REML", "reml", "rml", "RML"))
    method.lik <- "RML"
  if (method.lik %in% c("ML", "ml"))
    method.lik <- "ML"
  if (method.lik == "ML" & cov.model == "power")
    stop("\n\"power\" model can only be used with method.lik=\"RML\".\nBe sure that what you want is not \"powered.exponential\"")
  temp.list$method.lik <- method.lik
  coords <- as.matrix(coords)
  data <- as.vector(data)
  n <- length(data)
  if ((nrow(coords) != n) | (2 * n) != length(coords))
    stop("\nnumber of locations does not match with number of data")
  if (missing(geodata))
    xmat <- trend.spatial(trend = trend, geodata = list(coords = coords, data = data)) else xmat <- unclass(trend.spatial(trend = trend, geodata = geodata))
  xmat.contrasts <- attr(xmat, "contrasts")
  xmat <- unclass(xmat)
  if (nrow(xmat) != n)
    stop("trend matrix has dimension incompatible with the data")
  .solve.geoR(crossprod(xmat))
  beta.size <- temp.list$beta.size <- dim(xmat)[2]
  if (missing(realisations))
    realisations <- as.factor(rep(1, n)) else {
    if (!missing(geodata)) {
      real.name <- deparse(substitute(realisations))
      if (all(isTRUE(as.logical(real.name))))
        if (is.null(geodata$realisations))
          stop("element realisation not available in the geodata object") else realisations <- geodata$realisations else {
        if (!is.null(geodata[[real.name]]))
          realisations <- geodata[[real.name]]
      }
    }
    if (length(realisations) != n)
      stop("realisations must be a vector with the same length of the data")
    realisations <- as.factor(realisations)
  }
  temp.list$realisations <- realisations
  nrep <- temp.list$nrep <- length(levels(realisations))
  ind.rep <- split(1:n, realisations)
  vecdist <- function(x) {
    as.vector(stats::dist(x))
  }
  if (any(class(ini.cov.pars) == "eyefit")) {
    init <- nugget <- kappa <- NULL
    for (i in 1:length(ini.cov.pars)) {
      init <- drop(rbind(init, ini.cov.pars[[i]]$cov.pars))
      nugget <- c(nugget, ini.cov.pars[[i]]$nugget)
      if (cov.model == "gneiting.matern")
        kappa <- drop(rbind(kappa, ini.cov.pars[[i]]$kappa)) else kappa <- c(kappa, ini.cov.pars[[i]]$kappa)
    }
    ini.cov.pars <- init
  }
  if (any(class(ini.cov.pars) == "variomodel")) {
    nugget <- ini.cov.pars$nugget
    kappa <- ini.cov.pars$kappa
    ini.cov.pars <- ini.cov.pars$cov.pars
  }
  if (is.matrix(ini.cov.pars) | is.data.frame(ini.cov.pars)) {
    ini.cov.pars <- as.matrix(ini.cov.pars)
    if (nrow(ini.cov.pars) == 1)
      ini.cov.pars <- as.vector(ini.cov.pars) else {
      if ((cov.model != "pure.nugget") & (ncol(ini.cov.pars) != 2))
        stop("\nini.cov.pars must be a matrix or data.frame with 2 components: \ninitial values for sigmasq and phi")
    }
  }
  if (is.vector(ini.cov.pars)) {
    if ((cov.model != "pure.nugget") & (length(ini.cov.pars) != 2))
      stop("\nini.cov.pars must be a vector with 2 components: \ninitial values for sigmasq and phi")
  }
  if (is.matrix(ini.cov.pars) | (length(nugget) > 1) | (length(kappa) > 1) | (length(lambda) >
    1) | (length(psiR) > 1) | (length(psiA) > 1)) {
    if (messages.screen)
      cat("likfit: searching for best initial value ...")
    ini.temp <- matrix(ini.cov.pars, ncol = 2)
    ## HERE disable expand.grid. Inits are those proposed by the user
#     grid.ini <- as.matrix(expand.grid(sigmasq = unique(ini.temp[, 1]), phi = unique(ini.temp[,
#       2]), tausq = unique(nugget), kappa = unique(kappa), lambda = unique(lambda), psiR = unique(psiR),
#       psiA = unique(psiA)))
    grid.ini <- cbind(sigmasq = ini.temp[, 1], phi = ini.temp[,
      2], tausq = nugget, kappa = kappa, lambda = lambda, psiR = psiR,
      psiA = psiA)
    ## HERE dist.mat ##
    if (missing(dist.mat)) {
      assign(".likGRF.dists.vec", lapply(split(as.data.frame(coords), realisations),
        vecdist), pos = 1)
    } else {
      dist.vec <- dist.mat[which(col(dist.mat) < row(dist.mat))]
      assign(".likGRF.dists.vec", list(dist.vec), pos = 1)
    }

    temp.f <- function(parms, coords, data, temp.list) {
        res <- NA
        try(res <- loglik.GRF(geodata = geodata,
        coords = coords, data = as.vector(data), cov.model = temp.list$cov.model, cov.pars = parms[1:2],
        nugget = parms["tausq"], kappa = parms["kappa"], lambda = parms["lambda"], psiR = parms["psiR"],
        psiA = parms["psiA"], trend = trend, method.lik = temp.list$method.lik, compute.dists = FALSE,
        realisations = realisations), silent = TRUE)
        return(res)
      }
    grid.lik <- apply(grid.ini, 1, temp.f, coords = coords, data = data, temp.list = temp.list)
    grid.ini <- grid.ini[(grid.lik != Inf) & (grid.lik != -Inf) & !is.na(grid.lik) & !is.nan(grid.lik),
      , drop = FALSE]
    grid.lik <- grid.lik[(grid.lik != Inf) & (grid.lik != -Inf) & !is.na(grid.lik) & !is.nan(grid.lik)]
    ini.temp <- grid.ini[which(grid.lik == max(grid.lik)), , drop = FALSE]
    if (all(ini.temp[, "phi"] == 0))
      ini.temp <- ini.temp[1, , drop = FALSE]
    rownames(ini.temp) <- "initial.value"
    if (messages.screen) {
      cat(" selected values:\n")
      print(rbind(format(ini.temp, digits = 2), status = ifelse(c(FALSE, FALSE, fix.nugget,
        fix.kappa, fix.lambda, fix.psiR, fix.psiA), "fix", "est")))
      cat(paste("likelihood value:", max(grid.lik), "\n"))
    }
    dimnames(ini.temp) <- NULL
    ini.cov.pars <- ini.temp[1:2]
    nugget <- ini.temp[3]
    kappa <- ini.temp[4]
    lambda <- ini.temp[5]
    psiR <- ini.temp[6]
    psiA <- ini.temp[7]
    grid.ini <- NULL
    remove(".likGRF.dists.vec", pos = 1)
  }
  tausq <- nugget
  if (fix.lambda) {
    if (abs(lambda - 1) < 1e-04) {
      temp.list$log.jacobian <- 0
      temp.list$z <- as.vector(data)
    } else {
      if (any(data <= 0))
        stop("Transformation option not allowed when there are zeros or negative data")
      Jdata <- data^(lambda - 1)
      if (any(Jdata <= 0))
        temp.list$log.jacobian <- log(prod(Jdata)) else temp.list$log.jacobian <- sum(log(Jdata))
      Jdata <- NULL
      if (abs(lambda) < 1e-04)
        temp.list$z <- log(data) else temp.list$z <- ((data^lambda) - 1)/lambda
    }
  } else {
    temp.list$z <- as.vector(data)
    temp.list$log.jacobian <- NULL
  }
  if (fix.psiR & fix.psiA) {
    if (psiR != 1 | psiA != 0)
      coords <- coords.aniso(coords, aniso.pars = c(psiA, psiR))
    ## HERE dist.mat ##
    if (missing(dist.mat)) {
      assign(".likGRF.dists.vec", lapply(split(as.data.frame(coords), realisations),
        vecdist), pos = 1)
    } else {
      dist.vec <- dist.mat[which(col(dist.mat) < row(dist.mat))]
      assign(".likGRF.dists.vec", list(dist.vec), pos = 1)
    }

    range.dist <- range(get(".likGRF.dists.vec", pos = 1))
    max.dist <- max(range.dist)
    min.dist <- min(range.dist)
  }
  ini <- ini.cov.pars[2]
  lower.optim <- c(limits$phi["lower"])
  upper.optim <- c(limits$phi["upper"])
  fixed.values <- list()
  if (fix.nugget) {
    fixed.values$tausq <- nugget
  } else {
    ini <- c(ini, nugget/ini.cov.pars[1])
    lower.optim <- c(lower.optim, limits$tausq.rel["lower"])
    upper.optim <- c(upper.optim, limits$tausq.rel["upper"])
  }
  if (fix.kappa) {
    fixed.values$kappa <- kappa
  } else {
    ini <- c(ini, kappa)
    lower.optim <- c(lower.optim, limits$kappa["lower"])
    upper.optim <- c(upper.optim, limits$kappa["upper"])
  }
  if (fix.lambda) {
    fixed.values$lambda <- lambda
  } else {
    ini <- c(ini, lambda)
    lower.optim <- c(lower.optim, limits$lambda["lower"])
    upper.optim <- c(upper.optim, limits$lambda["upper"])
  }
  if (fix.psiR) {
    fixed.values$psiR <- psiR
  } else {
    ini <- c(ini, psiR)
    lower.optim <- c(lower.optim, limits$psiR["lower"])
    upper.optim <- c(upper.optim, limits$psiR["upper"])
  }
  if (fix.psiA) {
    fixed.values$psiA <- psiA
  } else {
    ini <- c(ini, psiA)
    lower.optim <- c(lower.optim, limits$psiA["lower"])
    upper.optim <- c(upper.optim, limits$psiA["upper"])
  }
  if (fix.nugget & nugget > 0) {
    ini <- c(ini, ini.cov.pars[1])
    lower.optim <- c(lower.optim, limits$sigmasq["lower"])
    upper.optim <- c(upper.optim, limits$sigmasq["upper"])
  }
  names(ini) <- NULL
  if (length(ini) == 1)
    justone <- TRUE else justone <- FALSE
  ip <- list(f.tausq = fix.nugget, f.kappa = fix.kappa, f.lambda = fix.lambda, f.psiR = fix.psiR,
    f.psiA = fix.psiA)
  npars <- beta.size + 2 + sum(unlist(ip) == FALSE)
  temp.list$coords <- coords
  temp.list$xmat <- split(as.data.frame(unclass(xmat)), realisations)
  temp.list$xmat <- lapply(temp.list$xmat, as.matrix)
  temp.list$n <- as.vector(unlist(lapply(temp.list$xmat, nrow)))
  temp.list$loglik.cte <- rep(0, nrep)
  for (i in 1:nrep) {
    if (method.lik == "ML") {
      if (ip$f.tausq & (tausq > 0))
        temp.list$loglik.cte[i] <- (temp.list$n[i]/2) * (-log(2 * pi)) else temp.list$loglik.cte[i] <- (temp.list$n[i]/2) * (-log(2 * pi) + log(temp.list$n[i]) -
        1)
    }
    if (method.lik == "RML") {
      xx.eigen <- eigen(crossprod(temp.list$xmat[[i]]), symmetric = TRUE, only.values = TRUE)
      if (ip$f.tausq & (tausq > 0))
        temp.list$loglik.cte[i] <- -((temp.list$n[i] - beta.size)/2) * (log(2 * pi)) +
          0.5 * sum(log(xx.eigen$values)) else temp.list$loglik.cte[i] <- -((temp.list$n[i] - beta.size)/2) * (log(2 * pi)) +
        ((temp.list$n[i] - beta.size)/2) * (log(temp.list$n[i] - beta.size)) - ((temp.list$n[i] -
        beta.size)/2) + 0.5 * sum(log(xx.eigen$values))
    }
  }
  if (messages.screen) {
    cat("---------------------------------------------------------------\n")
    cat("likfit: likelihood maximisation using the function ")
    if (is.R()) {
      if (justone)
        cat("optimize.\n") else cat("optim.\n")
    } else cat("nlminb.\n")
    cat("likfit: Use control() to pass additional\n         arguments for the maximisation function.")
    cat("\n        For further details see documentation for ")
    if (is.R()) {
      if (justone)
        cat("optimize.\n") else cat("optim.\n")
    } else cat("nlminb.\n")
    cat("likfit: It is highly advisable to run this function several\n        times with different initial values for the parameters.\n")
    cat("likfit: WARNING: This step can be time demanding!\n")
    cat("---------------------------------------------------------------\n")
  }
  if (length(ini) == 1) {
    if (upper.optim == Inf)
      upper.optim <- 50 * max.dist
    lik.minim <- do.call("optimize", c(list(geoR:::.negloglik.GRF, lower = lower.optim,
      upper = upper.optim, fp = fixed.values, ip = ip, temp.list = temp.list), ldots))
    lik.minim <- list(par = lik.minim$minimum, value = lik.minim$objective, convergence = 0,
      message = "function optimize used")
  } else {
    MET <- pmatch(names(ldots), names(formals(stats::optim)))
    if (is.na(MET) || all(names(formals(stats::optim))[MET] != "method"))
      ldots$method <- "L-BFGS-B"
    if (!is.null(names(ldots))) {
      names(ldots)[which(as.logical(pmatch(names(ldots), "method", nomatch = 0)))] <- "method"
    }
    if (!is.null(ldots$method) && ldots$method == "L-BFGS-B") {
      ldots$lower <- lower.optim
      ldots$upper <- upper.optim
    }
    lik.minim <- do.call("optim", c(list(par = ini, fn = geoR:::.negloglik.GRF, fp = fixed.values,
      ip = ip, temp.list = temp.list), ldots))
  }
  if (messages.screen)
    cat("likfit: end of numerical maximisation.\n")
  par.est <- lik.minim$par
  if (any(par.est < 0))
    par.est <- round(par.est, digits = 12)
  phi <- par.est[1]
  if (is.R())
    loglik.max <- -lik.minim$value else loglik.max <- -lik.minim$objective
  if (ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA) {
    psiA <- par.est[2]
  }
  if (ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA) {
    psiR <- par.est[2]
  }
  if (ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA) {
    psiR <- par.est[2]
    psiA <- par.est[3]
  }
  if (ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA) {
    lambda <- par.est[2]
  }
  if (ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA) {
    lambda <- par.est[2]
    psiA <- par.est[3]
  }
  if (ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA) {
    lambda <- par.est[2]
    psiR <- par.est[3]
  }
  if (ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA) {
    lambda <- par.est[2]
    psiR <- par.est[3]
    psiA <- par.est[4]
  }
  if (ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA) {
    kappa <- par.est[2]
  }
  if (ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA) {
    kappa <- par.est[2]
    psiA <- par.est[3]
  }
  if (ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA) {
    kappa <- par.est[2]
    psiR <- par.est[3]
  }
  if (ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA) {
    kappa <- par.est[2]
    psiR <- par.est[3]
    psiA <- par.est[4]
  }
  if (ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA) {
    kappa <- par.est[2]
    lambda <- par.est[3]
  }
  if (ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA) {
    kappa <- par.est[2]
    lambda <- par.est[3]
    psiA <- par.est[4]
  }
  if (ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA) {
    kappa <- par.est[2]
    lambda <- par.est[3]
    psiR <- par.est[4]
  }
  if (ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA) {
    kappa <- par.est[2]
    lambda <- par.est[3]
    psiR <- par.est[4]
    psiA <- par.est[5]
  }
  if (!ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA) {
    tausq <- par.est[2]
  }
  if (!ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA) {
    tausq <- par.est[2]
    psiA <- par.est[3]
  }
  if (!ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA) {
    tausq <- par.est[2]
    psiR <- par.est[3]
  }
  if (!ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA) {
    tausq <- par.est[2]
    psiR <- par.est[3]
    psiA <- par.est[4]
  }
  if (!ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA) {
    tausq <- par.est[2]
    lambda <- par.est[3]
  }
  if (!ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA) {
    tausq <- par.est[2]
    lambda <- par.est[3]
    psiA <- par.est[4]
  }
  if (!ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA) {
    tausq <- par.est[2]
    lambda <- par.est[3]
    psiR <- par.est[4]
  }
  if (!ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA) {
    tausq <- par.est[2]
    lambda <- par.est[3]
    psiR <- par.est[4]
    psiA <- par.est[5]
  }
  if (!ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA) {
    tausq <- par.est[2]
    kappa <- par.est[3]
  }
  if (!ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA) {
    tausq <- par.est[2]
    kappa <- par.est[3]
    psiA <- par.est[4]
  }
  if (!ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA) {
    tausq <- par.est[2]
    kappa <- par.est[3]
    psiR <- par.est[4]
  }
  if (!ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA) {
    tausq <- par.est[2]
    kappa <- par.est[3]
    psiR <- par.est[4]
    psiA <- par.est[5]
  }
  if (!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA) {
    tausq <- par.est[2]
    kappa <- par.est[3]
    lambda <- par.est[4]
  }
  if (!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA) {
    tausq <- par.est[2]
    kappa <- par.est[3]
    lambda <- par.est[4]
    psiA <- par.est[5]
  }
  if (!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA) {
    tausq <- par.est[2]
    kappa <- par.est[3]
    lambda <- par.est[4]
    psiR <- par.est[5]
  }
  if (!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA) {
    tausq <- par.est[2]
    kappa <- par.est[3]
    lambda <- par.est[4]
    psiR <- par.est[5]
    psiA <- par.est[6]
  }
  if (fix.nugget & nugget > 0) {
    sigmasq <- par.est[length(par.est)]
    if (sigmasq > 1e-12)
      tausq <- nugget/sigmasq
    check.sigmasq <- TRUE
  } else check.sigmasq <- FALSE
  if (!fix.lambda) {
    if (abs(lambda - 1) < 1e-04) {
      log.jacobian.max <- 0
    } else {
      if (any(data^(lambda - 1) <= 0))
        log.jacobian.max <- log(prod(data^(lambda - 1))) else log.jacobian.max <- sum(log(data^(lambda - 1)))
      temp.list$z <- ((data^lambda) - 1)/lambda
    }
  } else {
    log.jacobian.max <- temp.list$log.jacobian
  }
  data.rep <- split(temp.list$z, realisations)
  coords.rep <- split(as.data.frame(coords), realisations)
  coords.rep <- lapply(coords.rep, as.matrix)
  if (fix.psiR & fix.psiA)
    remove(".likGRF.dists.vec", pos = 1) else {
    if (round(psiR, digits = 6) != 1 | round(psiA, digits = 6) != 0)
      coords <- coords.aniso(coords, aniso.pars = c(psiA, psiR))
    rangevecdist <- function(x) {
      range(as.vector(stats::dist(x)))
    }
    ## HERE dist.mat ##
    if (missing(dist.mat)) {
      range.dist <- lapply(split(as.data.frame(coords), realisations), rangevecdist)
    } else {
      range.dist <- list(range(dist.mat))
    }

    range.dist <- range(as.vector(unlist(range.dist)))
    max.dist <- max(range.dist)
    min.dist <- min(range.dist)
  }
  xivx <- matrix(0, ncol = beta.size, nrow = beta.size)
  xivy <- matrix(0, ncol = 1, nrow = beta.size)
  yivy <- 0
  for (i in 1:nrep) {
    ni <- temp.list$n[i]
    if ((phi < 1e-12))
      V <- diag(x = (1 + tausq), ni) else {
      if (check.sigmasq) {
        if (sigmasq < 1e-12) {
          if (!fix.nugget)
          V <- diag(x = (1 + tausq), ni) else V <- diag(x = sqrt(tausq), ni)
        } else V <- varcov.spatial(coords = coords.rep[[i]], cov.model = cov.model, kappa = kappa,
          nugget = tausq, cov.pars = c(1, phi))$varcov
      } else V <- varcov.spatial(coords = coords.rep[[i]], cov.model = cov.model, kappa = kappa,
        nugget = tausq, cov.pars = c(1, phi))$varcov
    }
    ivyx <- solve(V, cbind(data.rep[[i]], temp.list$xmat[[i]]))
    xivx <- xivx + crossprod(ivyx[, -1], temp.list$xmat[[i]])
    xivy <- xivy + crossprod(ivyx[, -1], data.rep[[i]])
    yivy <- yivy + crossprod(data.rep[[i]], ivyx[, 1])
  }
  betahat <- .solve.geoR(xivx, xivy)
  res <- as.vector(temp.list$z - xmat %*% betahat)
  if (!fix.nugget | (nugget < 1e-12)) {
    ssres <- as.vector(yivy - 2 * crossprod(betahat, xivy) + crossprod(betahat, xivx) %*%
      betahat)
    if (method.lik == "ML")
      sigmasq <- ssres/n else sigmasq <- ssres/(n - beta.size)
  }
  if (fix.nugget) {
    if (nugget > 0)
      tausq <- nugget
  } else tausq <- tausq * sigmasq
  betahat.var <- .solve.geoR(xivx)
  if (sigmasq > 1e-12)
    betahat.var <- sigmasq * betahat.var
  if ((phi < 0.001 * min.dist)) {
    tausq <- tausq + sigmasq
    sigmasq <- 0
  }
  if ((sigmasq < 1e-12))
    phi <- 0
  n.model.pars <- beta.size + 7
  par.su <- data.frame(status = rep(-9, n.model.pars))
  ind.par.su <- c(rep(0, beta.size), ip$f.tausq, 0, 0, ip$f.kappa, ip$f.psiR, ip$f.psiA,
    ip$f.lambda)
  par.su$status <- ifelse(ind.par.su, "fixed", "estimated")
  par.su$values <- round(c(betahat, tausq, sigmasq, phi, kappa, psiR, psiA, lambda), digits = 4)
  if (beta.size == 1)
    beta.name <- "beta" else beta.name <- paste("beta", 0:(beta.size - 1), sep = "")
  row.names(par.su) <- c(beta.name, "tausq", "sigmasq", "phi", "kappa", "psiR", "psiA", "lambda")
  par.su <- par.su[c((1:(n.model.pars - 3)), n.model.pars - 1, n.model.pars - 2, n.model.pars),
    ]
  lik.results <- list(cov.model = cov.model, nugget = tausq, cov.pars = c(sigmasq, phi),
    sigmasq = sigmasq, phi = phi, kappa = kappa, beta = as.vector(betahat), beta.var = betahat.var,
    lambda = lambda, aniso.pars = c(psiA = psiA, psiR = psiR), tausq = tausq, practicalRange = practicalRange(cov.model = cov.model,
      phi = phi, kappa = kappa), method.lik = method.lik, trend = trend, loglik = loglik.max,
    npars = npars, AIC = -2 * (loglik.max - npars), BIC = -2 * (loglik.max - 0.5 * log(n) *
      npars), parameters.summary = par.su, info.minimisation.function = lik.minim, max.dist = max.dist,
    trend = trend, trend.matrix = xmat, transform.info = list(fix.lambda = fix.lambda,
      log.jacobian = log.jacobian.max))
  if (nospatial) {
    if (fix.lambda) {
      beta.ns <- .solve.geoR(crossprod(xmat), crossprod(xmat, temp.list$z))
      ss.ns <- sum((as.vector(temp.list$z - xmat %*% beta.ns))^2)
      if (method.lik == "ML") {
        nugget.ns <- ss.ns/n
        loglik.ns <- (n/2) * ((-log(2 * pi)) - log(nugget.ns) - 1) + temp.list$log.jacobian
      }
      if (method.lik == "RML") {
        nugget.ns <- ss.ns/(n - beta.size)
        loglik.ns <- ((n - beta.size)/2) * ((-log(2 * pi)) - log(nugget.ns) - 1) +
          temp.list$log.jacobian
      }
      npars.ns <- beta.size + 1 + (!fix.lambda)
      lambda.ns <- lambda
    } else {
      if (is.R())
        lik.lambda.ns <- stats::optim(par = 1, fn = geoR:::.negloglik.boxcox, method = "L-BFGS-B",
          lower = limits$lambda["lower"], upper = limits$lambda["upper"], data = data,
          xmat = xmat, lik.method = method.lik) else lik.lambda.ns <- stats::nlminb(par = 1, fn = geoR:::.negloglik.boxcox, lower = limits$lambda["lower"],
        upper = limits$lambda["upper"], data = data, xmat = xmat, lik.method = method.lik)
      lambda.ns <- lik.lambda.ns$par
      if (abs(lambda) < 1e-04)
        tdata.ns <- log(data) else tdata.ns <- ((data^lambda.ns) - 1)/lambda.ns
      beta.ns <- .solve.geoR(crossprod(xmat), crossprod(xmat, tdata.ns))
      ss.ns <- sum((as.vector(tdata.ns - xmat %*% beta.ns))^2)
      if (is.R())
        value.min.ns <- lik.lambda.ns$value else value.min.ns <- lik.lambda.ns$objective
      if (method.lik == "ML") {
        loglik.ns <- (-value.min.ns) + (n/2) * ((-log(2 * pi)) + log(n) - 1)
        nugget.ns <- ss.ns/n
      }
      if (method.lik == "RML") {
        nugget.ns <- ss.ns/(n - beta.size)
        loglik.ns <- (-value.min.ns) + ((n - beta.size)/2) * ((-log(2 * pi)) + log(n -
          beta.size) - 1)
      }
      npars.ns <- beta.size + 1 + (!fix.lambda)
    }
    lik.results$nospatial <- list(beta.ns = beta.ns, variance.ns = nugget.ns, loglik.ns = loglik.ns,
      npars.ns = npars.ns, lambda.ns = lambda.ns, AIC.ns = -2 * (loglik.ns - npars.ns),
      BIC.ns = -2 * (loglik.ns - 0.5 * log(n) * npars.ns))
  }
  if (length(lik.results$beta.var) == 1)
    lik.results$beta.var <- as.vector(lik.results$beta.var)
  if (length(lik.results$beta) > 1) {
    if (inherits(trend, "formula") || (length(class(trend)) > 0 && any(class(trend) ==
      "trend.spatial")))
      beta.names <- c("intercept", paste("covar", 1:(ncol(xmat) - 1), sep = "")) else if (trend == "1st")
      beta.names <- c("intercept", "x", "y") else if (trend == "2nd")
      beta.names <- c("intercept", "x", "y", "x2", "xy", "y2")
    names(lik.results$beta) <- beta.names
  }
  if (components) {
    if (!fix.psiR & !fix.psiA)
      if (psiR != 1 | psiA != 0)
        coords <- coords.aniso(coords, aniso.pars = c(psiA, psiR))
    trend.comp <- temp.list$z - res
    spatial.comp <- list()
    for (i in 1:nrep) {
      spatial.comp[[i]] <- as.vector(varcov.spatial(coords = coords[ind.rep[[i]], ],
        cov.model = cov.model, kappa = kappa, nugget = 0, cov.pars = c(sigmasq, phi))$varcov %*%
        varcov.spatial(coords = coords[ind.rep[[i]], ], cov.model = cov.model, kappa = kappa,
          nugget = tausq, cov.pars = c(sigmasq, phi), inv = TRUE)$inverse %*% res[ind.rep[[i]]])
    }
    spatial.comp <- as.vector(unlist(spatial.comp))[as.vector(unlist(ind.rep))]
    predict.comp <- trend.comp + spatial.comp
    residual.comp <- as.vector(temp.list$z - predict.comp)
    lik.results$model.components <- data.frame(trend = trend.comp, spatial = spatial.comp,
      residuals = residual.comp)
  }
  lik.results$contrasts <- xmat.contrasts
  lik.results$call <- call.fc
  oldClass(lik.results) <- c("likGRF", "variomodel")
  if (messages.screen) {
    if ((lik.results$cov.pars[1] < (0.01 * (lik.results$nugget + lik.results$cov.pars[1]))) &
      lik.results$cov.pars[2] > 0)
      cat("\nWARNING: estimated sill is less than 1 hundredth of the total variance. Consider re-examine the model excluding spatial dependence\n")
    if ((lik.results$cov.pars[2] > (10 * max.dist)) & lik.results$cov.pars[1] > 0)
      cat("\nWARNING: estimated range is more than 10 times bigger than the biggest distance between two points. Consider re-examine the model:\n 1) excluding spatial dependence if estimated sill is too low and/or \n 2) taking trends (covariates) into account\n")
    if (((lik.results$cov.pars[2] < (0.1 * min.dist)) & (lik.results$cov.pars[1] > 0)) &
      lik.results$cov.pars[2] > 0)
      cat("\nWARNING: estimated range is less than 1 tenth of the minimum distance between two points. Consider re-examine the model excluding spatial dependence\n")
  }
  attr(lik.results, "geodata") <- name.geodata
  return(lik.results)
}

#' Eyefit function of geoR with higher range for phi and sigma
#'
#' @inheritParams geoR::eyefit
#' @param max.phi Numeric as multiplied by max(vario$u). Default to 2.
#' @param max.sigma Numeric as multiplied by max(vario$v). Default to 2.
#'
#' @import tcltk
#'
#' @export

eyefit.large <- function(vario, max.phi = 2, max.sigma = 2, silent = FALSE)
{
    if (!requireNamespace("tcltk", quietly = TRUE))
        stop("package tcltk is required to run eyefit()")
    # library(tcltk)
    geterrmessage()
    done <- tclVar(0)
    eyefit.env <- new.env()
    assign("eyefit.tmp", list(), envir = eyefit.env)
    dmax <- max(vario$u)
    kappa1 <- tclVar("0.5")
    kappa2 <- tclVar("1.0")
    kernel <- tclVar("exponential")
    mdist <- tclVar(max(vario$u))
    nugget <- tclVar(0.1 * (max(vario$v)))
    sill <- tclVar(0.8 * (max(vario$v)))
    range <- tclVar(dmax/3)
    replot <- function(...) {
        k <- as.character(tclObj(kernel))
        kp1 <- as.numeric(tclObj(kappa1))
        kp2 <- as.numeric(tclObj(kappa2))
        maxdist <- as.numeric(tclObj(mdist))
        sigmasq <- as.numeric(tclObj(sill))
        phi <- as.numeric(tclObj(range))
        tausq <- as.numeric(tclObj(nugget))
        eval(substitute(graphics::plot(vario)))
        fit <- get("eyefit.tmp", envir = eyefit.env)
        lapply(fit, function(x) geoR::lines.variomodel(seq(0, maxdist,
            l = 100), cov.model = x$cov.model, kappa = x$kappa,
            cov.pars = x$cov.pars, nug = x$nug, max.dist = x$max.dist))
        if (k == "gneiting.matern" | k == "gencauchy") {
            geoR::lines.variomodel(x = seq(0, maxdist, l = 100), cov.model = k,
                kappa = c(kp1, kp2), cov.pars = c(sigmasq, phi),
                nug = tausq, max.dist = maxdist)
        }
        else if (k == "powered.exponential" || k == "cauchy" ||
            k == "matern") {
            geoR::lines.variomodel(x = seq(0, maxdist, l = 100), cov.model = k,
                kappa = kp1, cov.pars = c(sigmasq, phi), nug = tausq,
                max.dist = maxdist)
        }
        else geoR::lines.variomodel(x = seq(0, maxdist, l = 100), cov.model = k,
            cov.pars = c(sigmasq, phi), nug = tausq, max.dist = maxdist)
    }
    redraw <- function(...) {
        var <- as.character(tclObj(kernel))
        if (var == "gneiting.matern" | var == "gencauchy") {
          tkconfigure(entry.kappa1, state = "normal")
          tkconfigure(ts5, state = "normal")
          tkfocus(entry.kappa1)
          tkconfigure(entry.kappa2, state = "normal")
          tkconfigure(ts6, state = "normal")
        }
        else if (var == "powered.exponential" || var == "cauchy" ||
            var == "matern") {
          tkconfigure(entry.kappa1, state = "normal")
          tkconfigure(ts5, state = "normal")
          tkfocus(entry.kappa1)
          tkconfigure(entry.kappa2, state = "disabled")
          tkconfigure(ts6, state = "disabled")
        }
        else {
          tkconfigure(ts5, state = "disabled")
          tkconfigure(ts6, state = "disabled")
          tkconfigure(entry.kappa1, state = "disabled")
          tkconfigure(entry.kappa2, state = "disabled")
        }
        replot()
    }
    base <- tktoplevel()
    tkwm.title(base, "Eyefit 1.0")
    spec.frm <- tkframe(base, borderwidth = 2)
    left.frm <- tkframe(spec.frm)
    right.frm <- tkframe(spec.frm)
    frame1 <- tkframe(left.frm, relief = "groove", borderwidth = 2)
    tkpack(tklabel(frame1, text = "Parameters"), fill = "both",
        side = "top")
    entry.mdist <- tkentry(frame1, width = "4", textvariable = mdist)
    tkpack(ts1 <- tkscale(frame1, label = "Max. Distance", command = replot,
        from = 0, to = dmax, showvalue = 0, variable = mdist,
        resolution = 0.01, orient = "horiz", relief = "groove"),
        fill = "both", expand = 1, padx = 3, pady = 2, ipadx = 3,
        ipady = 2, side = "left")
    tkpack(entry.mdist, fill = "none", side = "right")
    frame3 <- tkframe(left.frm, relief = "groove", borderwidth = 2)
    tkpack(tklabel(frame3, text = "Cov. Parameters"), fill = "both",
        side = "top")
    entry.sill <- tkentry(frame3, width = "4", textvariable = sill)
    tkpack(ts2 <- tkscale(frame3, label = "Sill (sigmasq):",
        command = replot, from = 0, to = max.sigma * max(vario$v), showvalue = 0,
        variable = sill, resolution = 0.01, orient = "horiz",
        relief = "groove"), fill = "none", expand = 1, padx = 3,
        pady = 2, ipadx = 3, ipady = 2, side = "left")
    tkpack(entry.sill, side = "right")
    frame4 <- tkframe(left.frm, relief = "groove", borderwidth = 2)
    entry.range <- tkentry(frame4, width = "4", textvariable = range)
    tkpack(ts3 <- tkscale(frame4, label = "Range (phi):", command = replot,
        from = 0, to = max.phi * dmax, showvalue = 1, variable = range,
        resolution = 0.01, orient = "horiz", relief = "groove"),
        fill = "x", expand = 1, padx = 3, pady = 2, ipadx = 3,
        ipady = 2, side = "left")
    tkpack(entry.range, side = "right")
    frame5 <- tkframe(left.frm, relief = "groove", borderwidth = 2)
    tkpack(tklabel(frame5, text = "Nugget"), fill = "both", side = "top")
    entry.nugget <- tkentry(frame5, width = "4", textvariable = nugget)
    tkpack(ts4 <- tkscale(frame5, label = "Nugget (tausq):",
        command = replot, from = 0, to = 2 * max(vario$v), showvalue = 1,
        variable = nugget, resolution = 0.01, orient = "horiz",
        relief = "groove"), fill = "x", expand = 1, padx = 3,
        pady = 2, ipadx = 3, ipady = 2, side = "left")
    tkpack(entry.nugget, side = "right")
    frame6 <- tkframe(left.frm, relief = "groove", borderwidth = 2)
    tkpack(tklabel(frame6, text = "Kappa"), fill = "both", side = "top")
    entry.kappa1 <- tkentry(frame6, width = "4", textvariable = kappa1,
        state = "disabled")
    tkpack(ts5 <- tkscale(frame6, label = "Kappa 1:", command = replot,
        from = 0, to = 10, showvalue = 1, variable = kappa1,
        state = "disabled", resolution = 0.01, orient = "horiz",
        relief = "groove"), fill = "x", expand = 1, padx = 3,
        pady = 2, ipadx = 3, ipady = 2, side = "left")
    tkpack(entry.kappa1, side = "right", fill = "none")
    frame7 <- tkframe(left.frm, relief = "groove", borderwidth = 2)
    entry.kappa2 <- tkentry(frame7, width = "4", textvariable = kappa2,
        state = "disabled")
    tkpack(ts6 <- tkscale(frame7, label = "Kappa 2:", command = replot,
        from = 0, to = 10, showvalue = 1, variable = kappa2,
        state = "disabled", resolution = 0.01, orient = "horiz",
        relief = "groove"), fill = "x", expand = 1, padx = 3,
        pady = 2, ipadx = 3, ipady = 2, side = "left")
    tkpack(entry.kappa2, side = "right", fill = "none")
    frame2 <- tkframe(right.frm, relief = "groove", borderwidth = 2)
    tkpack(tklabel(frame2, text = "Function"))
    for (i in c("cauchy", "gencauchy", "circular", "cubic", "exponential",
        "gaussian", "gneiting", "gneiting.matern", "linear",
        "matern", "power", "powered.exponential", "pure.nugget",
        "spherical", "wave")) {
        tmp <- tkradiobutton(frame2, command = redraw, text = i,
            value = i, variable = kernel)
        tkpack(tmp, anchor = "w")
    }
    OnOK <- function() {
        replot()
    }
    OnQuit <- function() {
        tclvalue(done) <- 2
    }
    OnClear <- function(aux = vario) {
        assign("eyefit.tmp", list(), envir = eyefit.env)
        graphics::plot(aux)
    }
    OnSave <- function() {
        k <- as.character(tclObj(kernel))
        kp1 <- as.numeric(tclObj(kappa1))
        if (k == "gneiting.matern")
            kp2 <- as.numeric(tclObj(kappa2))
        else kp2 <- NULL
        maxdist <- as.numeric(tclObj(mdist))
        sigmasq <- as.numeric(tclObj(sill))
        phi <- as.numeric(tclObj(range))
        tausq <- as.numeric(tclObj(nugget))
        aux <- list(cov.model = k, cov.pars = c(sigmasq, phi),
            nugget = tausq, kappa = c(kp1, kp2), lambda = vario$lambda,
            trend = vario$trend, practicalRange = practicalRange(cov.model = k,
                phi = phi, kappa = kp1), max.dist = maxdist)
        oldClass(aux) <- "variomodel"
        assign("eyefit.tmp", c(get("eyefit.tmp", envir = eyefit.env),
            list(aux)), envir = eyefit.env)
        replot()
    }
    tkpack(frame1, frame3, frame4, frame5, frame6, frame7, fill = "x")
    tkpack(frame2, fill = "x")
    tkpack(left.frm, right.frm, side = "left", anchor = "n")
    c.but <- tkbutton(base, text = "Clear", command = function() {
        OnClear(vario)
    })
    q.but <- tkbutton(base, text = "Quit", command = OnQuit)
    save.but <- tkbutton(base, text = "Save", command = OnSave)
    tkpack(spec.frm)
    tkpack(q.but, side = "right")
    tkpack(c.but, side = "left")
    tkpack(save.but, side = "right")
    replot()
    tkbind(entry.kappa1, "<Return>", function() {
        replot()
    })
    tkbind(entry.kappa2, "<Return>", function() {
        replot()
    })
    tkbind(base, "<Destroy>", function() tclvalue(done) <- 2)
    tkwait.variable(done)
    tkdestroy(base)
    if (!silent) {
        fit <- get("eyefit.tmp", envir = eyefit.env)
        oldClass(fit) <- "eyefit"
        return(fit)
    }
    else return(invisible())
}


#' Automatically fitting a variogram
#'
#' Automatically fitting a variogram to the data on which it is applied. The automatic fitting is done through fit.variogram. In fit.variogram the user had to supply an initial estimate for the sill, range etc. autofitVariogram provides this estimate based on the data and then calls fit.variogram.
#' For details, see \code{\link[automap]{autofitVariogram}}
#'
#' @inheritParams automap::autofitVariogram
#' @param cutoff Parameter included in variogram (maximum distance of variogram).
#' Default to diagonal * 0.35.
#' @param dist.mat Square matrix of distances between points of the dataset ## Not implemented yet
#'
#' @export

autofitVariogram.dist <- function(formula, input_data, model = c("Sph", "Exp", "Gau",
    "Ste"), kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10), fix.values = c(NA,
    NA, NA), verbose = FALSE, GLS.model = NA, start_vals = c(NA,
    NA, NA), cutoff, miscFitOptions = list(), dist.mat = NULL, ...)
{
    if ("alpha" %in% names(list(...)))
        warning("Anisotropic variogram model fitting not supported, see the documentation of autofitVariogram for more details.")
    miscFitOptionsDefaults = list(merge.small.bins = TRUE, min.np.bin = 5)
    miscFitOptions = utils::modifyList(miscFitOptionsDefaults, miscFitOptions)
    longlat = !sp::is.projected(input_data)
    if (is.na(longlat)) { longlat = FALSE }
      diagonal = sp::spDists(t(sp::bbox(input_data)), longlat = longlat)[1, 2]
      # Add possibility to define cutoff (maximum distance of variogram)
      if (missing(cutoff)) {cutoff <- diagonal * 0.35}
      boundaries = c(2, 4, 6, 9, 12, 15, 25, 35, 50, 65, 80, 100) *
        # diagonal * 0.35/100 # Change to cutoff
        cutoff/100
    if (!is(GLS.model, "variogramModel")) {
        experimental_variogram = gstat::variogram(formula, input_data,
            cutoff = cutoff, boundaries = boundaries, ...) # Added cutoff
    }
    else {
        if (verbose)
            cat("Calculating GLS sample variogram\n")
        g = gstat::gstat(NULL, "bla", formula, input_data, model = GLS.model,
            set = list(gls = 1))
        experimental_variogram = gstat::variogram(g, cutoff = cutoff, # Added cutoff
            boundaries = boundaries,
            ...)
    }
    if (miscFitOptions[["merge.small.bins"]]) {
        if (verbose)
            cat("Checking if any bins have less than 5 points, merging bins when necessary...\n\n")
        while (TRUE) {
            if (length(experimental_variogram$np[experimental_variogram$np <
                miscFitOptions[["min.np.bin"]]]) == 0 | length(boundaries) ==
                1)
                break
            boundaries = boundaries[2:length(boundaries)]
            if (!is(GLS.model, "variogramModel")) {
                experimental_variogram = gstat::variogram(formula, input_data,
                  cutoff = cutoff, boundaries = boundaries, ...) # Added cutoff
            }
            else {
                experimental_variogram = gstat::variogram(g, cutoff = cutoff, # Added cutoff
                  boundaries = boundaries,
                  ...)
            }
        }
    }
    if (is.na(start_vals[1])) {
        initial_nugget = min(experimental_variogram$gamma)
    }
    else {
        initial_nugget = start_vals[1]
    }
    if (is.na(start_vals[2])) {
        initial_range = 0.1 * diagonal
    }
    else {
        initial_range = start_vals[2]
    }
    if (is.na(start_vals[3])) {
        initial_sill = mean(c(max(experimental_variogram$gamma),
            stats::median(experimental_variogram$gamma)))
    }
    else {
        initial_sill = start_vals[3]
    }
    if (!is.na(fix.values[1])) {
        fit_nugget = FALSE
        initial_nugget = fix.values[1]
    }
    else fit_nugget = TRUE
    if (!is.na(fix.values[2])) {
        fit_range = FALSE
        initial_range = fix.values[2]
    }
    else fit_range = TRUE
    if (!is.na(fix.values[3])) {
        fit_sill = FALSE
        initial_sill = fix.values[3]
    }
    else fit_sill = TRUE
    getModel = function(psill, model, range, kappa, nugget, fit_range,
        fit_sill, fit_nugget, verbose) {
        if (verbose)
            debug.level = 1
        else debug.level = 0
        if (model == "Pow") {
            warning("Using the power model is at your own risk, read the docs of autofitVariogram for more details.")
            if (is.na(start_vals[1]))
                nugget = 0
            if (is.na(start_vals[2]))
                range = 1
            if (is.na(start_vals[3]))
                sill = 1
        }
        obj = try(gstat::fit.variogram(experimental_variogram, model = gstat::vgm(psill = psill,
            model = model, range = range, nugget = nugget, kappa = kappa),
            fit.ranges = c(fit_range), fit.sills = c(fit_nugget,
                fit_sill), debug.level = 0), TRUE)
        if ("try-error" %in% class(obj)) {
            warning("An error has occured during variogram fitting. Used:\n",
                "\tnugget:\t", nugget, "\n\tmodel:\t", model,
                "\n\tpsill:\t", psill, "\n\trange:\t", range,
                "\n\tkappa:\t", ifelse(kappa == 0, NA, kappa),
                "\n  as initial guess. This particular variogram fit is not taken into account. \nGstat error:\n",
                obj)
            return(NULL)
        }
        else return(obj)
    }
    test_models = model
    SSerr_list = c()
    vgm_list = list()
    counter = 1
    for (m in test_models) {
        if (m != "Mat" && m != "Ste") {
            model_fit = getModel(initial_sill - initial_nugget,
                m, initial_range, kappa = 0, initial_nugget,
                fit_range, fit_sill, fit_nugget, verbose = verbose)
            if (!is.null(model_fit)) {
                vgm_list[[counter]] = model_fit
                SSerr_list = c(SSerr_list, attr(model_fit, "SSErr"))
            }
            counter = counter + 1
        }
        else {
            for (k in kappa) {
                model_fit = getModel(initial_sill - initial_nugget,
                  m, initial_range, k, initial_nugget, fit_range,
                  fit_sill, fit_nugget, verbose = verbose)
                if (!is.null(model_fit)) {
                  vgm_list[[counter]] = model_fit
                  SSerr_list = c(SSerr_list, attr(model_fit,
                    "SSErr"))
                }
                counter = counter + 1
            }
        }
    }
    strange_entries = sapply(vgm_list, function(v) any(c(v$psill,
        v$range) < 0) | is.null(v))
    if (any(strange_entries)) {
        if (verbose) {
            print(vgm_list[strange_entries])
            cat("^^^ ABOVE MODELS WERE REMOVED ^^^\n\n")
        }
        warning("Some models where removed for being either NULL or having a negative sill/range/nugget, \n\tset verbose == TRUE for more information")
        SSerr_list = SSerr_list[!strange_entries]
        vgm_list = vgm_list[!strange_entries]
    }
    if (verbose) {
        cat("Selected:\n")
        print(vgm_list[[which.min(SSerr_list)]])
        cat("\nTested models, best first:\n")
        tested = data.frame(`Tested models` = sapply(vgm_list,
            function(x) as.character(x[2, 1])), kappa = sapply(vgm_list,
            function(x) as.character(x[2, 4])), SSerror = SSerr_list)
        tested = tested[order(tested$SSerror), ]
        print(tested)
    }
    result = list(exp_var = experimental_variogram, var_model = vgm_list[[which.min(SSerr_list)]],
        sserr = min(SSerr_list))
    class(result) = c("autofitVariogram", "list")
    return(result)
}
