#' Kriging. Modified from function krige.conv in geoR library.
#'
#' @inheritParams geoR::krige.conv
#' @param dist.mat Square matrix of distances between points of the dataset
#' @param loc.dist is a n (data) x N (locations) matrix with distances between data points and prediction locations.
#' @param loc.loc.dist is a N (locations) x N (locations) matrix with distances between prediction locations.
#'
#' @import geoR
#'
#' @export

krige.conv.dist <- function(geodata, coords = geodata$coords, data = geodata$data, locations,
                            borders, krige, output, dist.mat, loc.dist, loc.loc.dist)
{
    if (missing(geodata))
      geodata <- list(coords = coords, data = data)
    if (missing(borders))
      borders <- geodata$borders
    locations <- geoR::.check.locations(locations)
    krige1d <- FALSE
    if (length(locations) > 2) {
      if (length(unique(locations[, 1])) == 1 | length(unique(locations[, 2])) == 1)
        krige1d <- TRUE
    }
    call.fc <- match.call()
    base.env <- sys.frame(sys.nframe())
    if (missing(krige))
      krige <- geoR::krige.control() else {
        if (length(class(krige)) == 0 || class(krige) != "krige.geoR") {
          if (!is.list(krige))
            stop("krige.conv: the argument krige only takes a list or an output of the function krige.control") else {
              krige.names <- c("type.krige", "trend.d", "trend.l", "obj.model", "beta", "cov.model",
                               "cov.pars", "kappa", "nugget", "micro.scale", "dist.epsilon", "lambda", "aniso.pars")
              krige.user <- krige
              krige <- list()
              if (length(krige.user) > 0) {
                for (i in 1:length(krige.user)) {
                  n.match <- match.arg(names(krige.user)[i], krige.names)
                  krige[[n.match]] <- krige.user[[i]]
                }
              }
              if (is.null(krige$type.krige))
                krige$type.krige <- "ok"
              if (is.null(krige$trend.d))
                krige$trend.d <- "cte"
              if (is.null(krige$trend.l))
                krige$trend.l <- "cte"
              if (is.null(krige$obj.model))
                krige$obj.model <- NULL
              if (is.null(krige$beta))
                krige$beta <- NULL
              if (is.null(krige$cov.model))
                krige$cov.model <- "matern"
              if (is.null(krige$cov.pars))
                stop("covariance parameters (sigmasq and phi) should be provided in cov.pars")
              if (is.null(krige$kappa))
                krige$kappa <- 0.5
              if (is.null(krige$nugget))
                krige$nugget <- 0
              if (is.null(krige$micro.scale))
                krige$micro.scale <- 0
              if (is.null(krige$dist.epsilon))
                krige$dist.epsilon <- 1e-10
              if (is.null(krige$aniso.pars))
                krige$aniso.pars <- NULL
              if (is.null(krige$lambda))
                krige$lambda <- 1
              krige <- geoR::krige.control(type.krige = krige$type.krige, trend.d = krige$trend.d,
                                           trend.l = krige$trend.l, obj.model = krige$obj.model, beta = krige$beta,
                                           cov.model = krige$cov.model, cov.pars = krige$cov.pars, kappa = krige$kappa,
                                           nugget = krige$nugget, micro.scale = krige$micro.scale, dist.epsilon = krige$dist.epsilon,
                                           aniso.pars = krige$aniso.pars, lambda = krige$lambda)
            }
        }
      }
    cov.model <- krige$cov.model
    kappa <- krige$kappa
    lambda <- krige$lambda
    beta <- krige$beta
    cov.pars <- krige$cov.pars
    nugget <- krige$nugget
    micro.scale <- krige$micro.scale
    aniso.pars <- krige$aniso.pars
    if (missing(output))
      output <- geoR::output.control() else {
        if (length(class(krige)) == 0 || class(output) != "output.geoR") {
          if (!is.list(output))
            stop("krige.conv: the argument output can take only a list or an output of the function output.control") else {
              output.names <- c("n.posterior", "n.predictive", "moments", "n.back.moments",
                                "simulations.predictive", "mean.var", "quantile", "threshold", "signal",
                                "messages.screen")
              output.user <- output
              output <- list()
              if (length(output.user) > 0) {
                for (i in 1:length(output.user)) {
                  n.match <- match.arg(names(output.user)[i], output.names)
                  output[[n.match]] <- output.user[[i]]
                }
              }
              if (is.null(output$n.posterior))
                output$n.posterior <- 1000
              if (is.null(output$n.predictive))
                output$n.predictive <- NULL
              if (is.null(output$moments))
                output$moments <- TRUE
              if (is.null(output$n.back.moments))
                output$n.back.moments <- 1000
              if (is.null(output$simulations.predictive)) {
                if (is.null(output$n.predictive))
                  output$simulations.predictive <- NULL else output$simulations.predictive <- ifelse(output$n.predictive > 0, TRUE,
                                                                                                     FALSE)
              }
              if (is.null(output$mean.var))
                output$mean.var <- NULL
              if (is.null(output$quantile))
                output$quantile <- NULL
              if (is.null(output$threshold))
                output$threshold <- NULL
              if (is.null(output$sim.means))
                output$sim.means <- NULL
              if (is.null(output$sim.vars))
                output$sim.vars <- NULL
              if (is.null(output$signal))
                output$signal <- NULL
              if (is.null(output$messages.screen))
                output$messages.screen <- TRUE
              output <- geoR::output.control(n.posterior = output$n.posterior, n.predictive = output$n.predictive,
                                             moments = output$moments, n.back.moments = output$n.back.moments, simulations.predictive = output$simulations.predictive,
                                             mean.var = output$mean.var, quantile = output$quantile, threshold = output$threshold,
                                             sim.means = output$sim.means, sim.vars = output$sim.vars, signal = output$signal,
                                             messages = output$messages.screen)
            }
        }
      }
    signal <- ifelse(is.null(output$signal), FALSE, output$signal)
    messages.screen <- output$messages.screen
    n.predictive <- output$n.predictive
    n.back.moments <- output$n.back.moments
    n.predictive <- ifelse(is.null(n.predictive), 0, n.predictive)
    simulations.predictive <- ifelse(is.null(output$simulations.predictive), FALSE, TRUE)
    mean.estimator <- output$mean.estimator
    sim.means <- output$sim.means
    if (is.null(sim.means))
      sim.means <- ifelse(simulations.predictive, TRUE, FALSE)
    sim.vars <- output$sim.vars
    if (is.null(sim.vars))
      sim.vars <- FALSE
    if (is.null(mean.estimator) & simulations.predictive)
      mean.estimator <- TRUE
    quantile.estimator <- output$quantile.estimator
    probability.estimator <- output$probability.estimator
    if (!is.null(probability.estimator)) {
      if (length(probability.estimator) > 1 & length(probability.estimator) != nrow(locations))
        stop("krige.conv: probability.estimator must either have length 1, or have length = nrow(locations)\n")
    }
    if (simulations.predictive & n.predictive == 0)
      n.predictive <- 1000
    if (krige$type.krige == "ok")
      beta.prior <- "flat"
    if (krige$type.krige == "sk")
      beta.prior <- "deg"
    if (is.vector(coords)) {
      coords <- cbind(coords, 0)
      warning("krige.conv: coordinates provided as a vector, assuming one spatial dimension")
    }
    coords <- as.matrix(coords)
    if (!is.null(borders)) {
      nloc0 <- nrow(locations)
      ind.loc0 <- geoR::.geoR_inout(locations, borders)
      if (nrow(locations) == 1)
        locations <- locations[ind.loc0, , drop = FALSE] else locations <- locations[ind.loc0, , drop = TRUE]
      if (nrow(locations) == 0)
        stop("\nkrige.conv: there are no prediction locations inside the borders")
      if (messages.screen)
        cat("krige.conv: results will be returned only for prediction locations inside the borders\n")
    }
    dimnames(coords) <- list(NULL, NULL)
    dimnames(locations) <- list(NULL, NULL)
    if (messages.screen) {
      if (mode(krige$trend.d) == "numeric")
        cat("krige.conv: model with covariates matrix provided by the user") else cat(switch(as.character(krige$trend.d)[1], cte = "krige.conv: model with constant mean",
                                                                                             `1st` = "krige.conv: model with mean given by a 1st order polynomial on the coordinates",
                                                                                             `2nd` = "krige.conv: model with mean given by a 2nd order polynomial on the coordinates",
                                                                                             "krige.conv: model with mean defined by covariates provided by the user"))
      cat("\n")
    }
    if (class(krige$trend.d) == "trend.spatial")
      trend.d <- unclass(krige$trend.d) else trend.d <- unclass(geoR::trend.spatial(trend = krige$trend.d, geodata = geodata))
    if (nrow(trend.d) != nrow(coords))
      stop("coords and trend.d have incompatible sizes")
    beta.size <- ncol(trend.d)
    if (beta.prior == "deg")
      if (beta.size != length(beta))
        stop("size of mean vector is incompatible with trend specified")
    if (class(krige$trend.l) == "trend.spatial")
      trend.l <- unclass(krige$trend.l) else trend.l <- unclass(trend.spatial(trend = krige$trend.l, geodata = list(coords = locations)))
    if (!is.null(borders))
      if (nrow(trend.l) == nloc0)
        trend.l <- trend.l[ind.loc0, , drop = FALSE]
    if (nrow(trend.l) != nrow(locations))
      stop("locations and trend.l have incompatible sizes")
    if (beta.size > 1)
      beta.names <- paste("beta", (0:(beta.size - 1)), sep = "") else beta.names <- "beta"
    if (!is.null(aniso.pars)) {
      if (abs(aniso.pars[2] - 1) > 1e-04) {
        if (messages.screen)
          cat("krige.conv: anisotropy correction performed\n")
        coords <- geoR::coords.aniso(coords = coords, aniso.pars = aniso.pars)
        locations <- geoR::coords.aniso(coords = locations, aniso.pars = aniso.pars)
      }
    }
    if (!isTRUE(all.equal(lambda, 1))) {
      if (messages.screen)
        cat("krige.conv: performing the Box-Cox data transformation\n")
      data <- geoR::BCtransform(x = data, lambda = lambda)$data
    }
    if (is.vector(cov.pars)) {
      sigmasq <- cov.pars[1]
      phi <- cov.pars[2]
      cpars <- c(1, phi)
    } else {
      stop("current version of krige.conv does not accept nested covariance models\n")
    }
    sill.partial <- sum(sigmasq)
    if (sill.partial < 1e-16) {
      tausq.rel <- 0
      tausq.rel.micro <- 0
    } else {
      tausq.rel <- nugget/sum(sigmasq)
      tausq.rel.micro <- micro.scale/sum(sigmasq)
    }
    n <- length(data)
    ni <- nrow(trend.l)
    kc <- list()
    #Vcov <- geoR::varcov.spatial(coords = coords, cov.model = cov.model, kappa = kappa, nugget = tausq.rel,
    #                             cov.pars = cpars)$varcov
    Vcov <- geoR::varcov.spatial(dists.lowertri = as.vector(as.dist(dist.mat)),
                                 cov.model = cov.model, kappa = kappa, nugget = tausq.rel,
                                 cov.pars = cpars)$varcov
    ivtt <- solve(Vcov, trend.d)
    ttivtt <- crossprod(ivtt, trend.d)
    if (beta.prior == "flat")
      beta.flat <- drop(solve(ttivtt, crossprod(ivtt, as.vector(data))))
    remove("ivtt")
    ## Allow for personal distance matrix
    if (missing(loc.dist)) {
      v0 <- geoR::loccoords(coords = coords, locations = locations)
    } else {
      v0 <- loc.dist
    }
    if (n.predictive > 0) {
      loc.coincide <- apply(v0, 2, function(x, min.dist) {
        any(x < min.dist)
      }, min.dist = krige$dist.epsilon)
      if (any(loc.coincide))
        loc.coincide <- (1:ni)[loc.coincide]
      else loc.coincide <- NULL
      if (!is.null(loc.coincide)) {
          temp.f <- function(x, data, dist.eps) {
            return(data[x < dist.eps])
          }
          data.coincide <- apply(v0[, loc.coincide, drop = FALSE], 2, temp.f, data = data,
                                 dist.eps = krige$dist.epsilon)
      } else data.coincide <- NULL
    } else remove("locations")
    nug.factor <- ifelse(signal, tausq.rel.micro, tausq.rel)
    ind.v0 <- which(v0 < krige$dist.epsilon)
    v0 <- cov.spatial(obj = v0, cov.model = cov.model, kappa = kappa, cov.pars = cpars)
    v0[ind.v0] <- 1 + nug.factor
    ivv0 <- solve(Vcov, v0)
    tv0ivv0 <- colSums(v0 * ivv0)
    remove("Vcov", "ind.v0")
    b <- crossprod(cbind(data, trend.d), ivv0)
    if (n.predictive == 0)
      remove("v0", "ivv0")
    tv0ivdata <- drop(b[1, ])
    b <- t(trend.l) - b[-1, , drop = FALSE]
    if (beta.prior == "deg") {
      kc$predict <- tv0ivdata + drop(crossprod(b, beta))
      kc$krige.var <- sill.partial * drop(1 + nug.factor - tv0ivv0)
      kc$beta.est <- "Simple kriging performed (beta provided by user)"
    }
    if (beta.prior == "flat") {
      kc$predict <- tv0ivdata + drop(crossprod(b, beta.flat))
      bitb <- colSums(b * solve(ttivtt, b))
      kc$krige.var <- sill.partial * drop(1 + nug.factor - tv0ivv0 + bitb)
      kc$beta.est <- beta.flat
      names(kc$beta.est) <- beta.names
      if (n.predictive == 0)
        remove("bitb")
    }
    if (n.predictive == 0)
      remove("b", "tv0ivv0")
    remove("tv0ivdata")
    kc$distribution <- "normal"
    kc$krige.var[kc$krige.var < 1e-08] <- 0
    if (any(kc$krige.var < 0))
      cat("krige.conv: negative kriging variance found! Investigate why this is happening.\n")
    if (n.predictive > 0) {
      if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        warning(".Random.seed not initialised. Creating it with runif(1)")
        runif(1)
      }
      seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      if (messages.screen)
        cat("krige.conv: sampling from the predictive distribution (conditional simulations)\n")
      if (length(cov.pars) > 2) {
        if (beta.prior == "flat") {
          ok.add.var <- drop(bitb)
          tv0ivv0 <- tv0ivv0 + ok.add.var
        }
        # varcov <- (varcov.spatial(coords = locations, cov.model = cov.model, cov.pars = cov.pars,
        #                          kappa = kappa, nugget = nugget)$varcov) - tv0ivv0
        varcov <- (varcov.spatial(dists.lowertri = as.vector(as.dist(dist.mat)),
                                  cov.model = cov.model, cov.pars = cov.pars,
                                  kappa = kappa, nugget = nugget)$varcov) - tv0ivv0
        remove("tv0ivv0")
        kc$simulations <- kc$predict + crossprod(chol(varcov), matrix(rnorm(ni * n.predictive),
                                                                      ncol = n.predictive))
      } else {
        coincide.cond <- ((isTRUE(all.equal(nugget, 0)) | !signal) & (!is.null(loc.coincide)))
        if (coincide.cond) {
          nloc <- ni - length(loc.coincide)
          ind.not.coincide <- -(loc.coincide)
          v0 <- v0[, ind.not.coincide, drop = FALSE]
          b <- b[, ind.not.coincide, drop = FALSE]
        } else {
          nloc <- ni
          ind.not.coincide <- TRUE
        }
        Dval <- 1 + nug.factor
        if (beta.prior == "deg")
          vbetai <- matrix(0, ncol = beta.size, nrow = beta.size)
        else vbetai <- matrix(.solve.geoR(ttivtt), ncol = beta.size, nrow = beta.size)
        df.model <- ifelse(beta.prior == "deg", n, n - beta.size)
        kc$simulations <- matrix(NA, nrow = ni, ncol = n.predictive)
        if (nloc > 0) {
          # invcov <- varcov.spatial(coords = coords, cov.model = cov.model, kappa = kappa,
          #                          nugget = tausq.rel, cov.pars = cpars, inv = TRUE, only.inv.lower.diag = TRUE)
          invcov <- varcov.spatial(dists.lowertri = as.vector(as.dist(dist.mat)),
                                   cov.model = cov.model, kappa = kappa,
                                   nugget = tausq.rel, cov.pars = cpars, inv = TRUE, only.inv.lower.diag = TRUE)
          kc$simulations[ind.not.coincide, ] <- .cond.sim(env.loc = base.env, env.iter = base.env,
                                                          loc.coincide = loc.coincide, coincide.cond = coincide.cond, tmean = kc$predict[ind.not.coincide],
                                                          Rinv = invcov, mod = list(beta.size = beta.size, nloc = nloc, Nsims = n.predictive,
                                                                                    n = n, Dval = Dval, df.model = df.model, s2 = sill.partial, cov.model.number = .cor.number(cov.model),
                                                                                    phi = phi, kappa = kappa, nugget = nugget), vbetai = vbetai, fixed.sigmasq = TRUE)
          remove("invcov")
        }
        remove("b", "locations")
        if (coincide.cond)
          kc$simulations[loc.coincide, ] <- rep(data.coincide, n.predictive)
      }
      if (!isTRUE(all.equal(lambda, 1))) {
        if (messages.screen)
          cat("krige.conv: back-transforming the simulated values\n")
        if (any(kc$simulations < -1/lambda))
          warning("Truncation in the back-transformation: there are simulated values less than (- 1/lambda) in the normal scale.")
        kc$simulations <- BCtransform(x = kc$simulations, lambda = lambda, inverse = TRUE)$data
      }
      if (!is.null(mean.estimator) | !is.null(quantile.estimator) | !is.null(probability.estimator) |
          !is.null(sim.means) | !is.null(sim.vars)) {
        kc <- c(kc, statistics.predictive(simuls = kc$simulations, mean.var = mean.estimator,
                                          quantile = quantile.estimator, threshold = probability.estimator, sim.means = sim.means,
                                          sim.vars = sim.vars))
      }
      kc$.Random.seed <- seed
    }
    if (!isTRUE(all.equal(lambda, 1))) {
      if (messages.screen) {
        cat("krige.conv: back-transforming the predicted mean and variance\n")
        if (!isTRUE(all.equal(lambda, 0)) & !isTRUE(all.equal(lambda, 0.5)))
          cat("krige.conv: back-transforming by simulating from the predictive.\n           (run the function a few times and check stability of the results.\n")
      }
      kc[c("predict", "krige.var")] <- backtransform.moments(lambda = lambda, mean = kc$predict,
                                                             variance = kc$krige.var, distribution = "normal", n.simul = n.back.moments)[c("mean",
                                                                                                                                           "variance")]
    }
    message <- "krige.conv: Kriging performed using global neighbourhood"
    if (messages.screen)
      cat(paste(message, "\n"))
    kc$message <- message
    kc$call <- call.fc
    attr(kc, "sp.dim") <- ifelse(krige1d, "1d", "2d")
    attr(kc, "prediction.locations") <- call.fc$locations
    attr(kc, "parent.env") <- parent.frame()
    if (!is.null(call.fc$coords))
      attr(kc, "data.locations") <- call.fc$coords else attr(kc, "data.locations") <- substitute(a$coords, list(a = substitute(geodata)))
    if (!is.null(call.fc$borders))
      attr(kc, "borders") <- call.fc$borders
    oldClass(kc) <- "kriging"
    return(kc)
}
