require(Rcpp)
require(parallel)
require(dirichletprocess)
require(abind)
require(pracma)
require(ggplot2)

sourceCpp('mdpolya.cpp')

#' Marginal MDP Sampler
#'
#' @param data A numeric vector of `n` observation values.
#' @param k The number of sampling iterations to record, will be truncated if a
#'  non-integer.
#' @param alpha The concentration parameter for the Dirichlet Process prior.
#' @param mu The mean parameter for the Normal-Inverse-Gamma prior.
#' @param tau The variance parameter for the Normal-Inverse-Gamma prior.
#' @param s The shape parameter for the Normal-Inverse-Gamma prior.
#' @param S The scale parameter for the Normal-Inverse-Gamma prior.
#' @param c The shape parameter for the Gamma prior on alpha.
#' @param C The scale parameter for the Gamma prior on alpha.
#' @param a The mean parameter for the Normal prior on mu.
#' @param A The variance parameter for the Normal prior on mu.
#' @param w The shape parameter for the Inverse-Gamma prior on tau.
#' @param W The scale parameter for the Inverse-Gamma prior on tau.
#' @param fix_a A logical value indicating whether or not to fix alpha at its
#'  initial value.
#' @param fix_m A logical value indicating whether or not to fix mu at its
#'  initial value.
#' @param fix_t A logical value indicating whether or not to fix tau at its
#'  initial value.
#' @param burn The number of initial sampling iterations to discard, will be
#'  truncated if a non-integer.
#' @param thin The number of sampling iterations to discard between records,
#'  will be truncated if a non-integer.
#'
#' @return
#' An `mdpolya_result` object. A list with four entries:
#'  * `theta`: An array of dimension [`k`, `n`, 3] encoding the component label,
#'   mean, and standard deviation for each data point for each iteration. This
#'   represents the samples from the Polya posterior distribution of the
#'   marginal MDP model.
#'  * `eta`: A matrix of dimension ['k', 5] encoding the hyperparameter values
#'   for each iteration.
#'  * `args`: A list of input arguments.
#'  * `phi`: A list of matrices encoding the unique values from `theta` and
#'   associated weights for each iteration.
mdp <- function(data, k, alpha = 1, mu = 21, tau = 25, s = 4, S = 2,
                c = 2, C = 4, a = 21, A = 21, w = 1, W = 100,
                fix_a = FALSE, fix_m = FALSE, fix_t = FALSE,
                burn = 1000, thin = 150) {
  out <- mdp_cpp(data, k, alpha, mu, tau, s, S, c, C, a, A, w, W,
                fix_a, fix_m, fix_t, burn, thin)
  out[[3]] <- list(data = data, k = k, alpha = alpha, mu = mu, tau = tau,
                   s = S, S = S, c = c, C = C, a = a, A = A, w = w, W = W,
                   fix_a = fix_a, fix_m = fix_m, fix_t = fix_t,
                   burn = burn, thin = thin)
  colnames(out[[1]]) <- names(data)
  dimnames(out[[1]])[[3]] <- c('z', 'mean', 'var')
  colnames(out[[2]]) <- c('alpha', 'mu', 'tau', 's', 'S')
  names(out) <- c('theta', 'eta', 'args')
  phi <- apply(out$theta, 1, function(t) {
    z <- t[, 'z']
    z_uq <- sort(unique(z))
    z_w <- table(z) / length(z)
    mat <- matrix(c(z_w, t[z_uq, 'mean'], t[z_uq, 'var']), ncol = 3)
    colnames(mat) <- c('w', 'mean', 'var')
    return(mat)
  })
  out$phi <- phi
  class(out) <- c('mdpolya_result', 'mdp')
  return(out)
}

#' Polya Completion for Marginal MDP Samples
#'
#' @param res_mdp Samples from a marginal MDP model.
#' @param epsilon The desired maximum weight associated with the final
#'  remainder component.
#' @param upsilon The portion of samples which do not meet the desired epsilon.
#' @param nthreads UNSTABLE: The number of parallel threads to launch with OpenMP,
#'  not recommended due to induced instability.
#'
#' @return If `res_mdp` was an `mdpolya_result` object, returns another
#'  `mdpolya_result` object with `phi`, `eta` and `args` entries as in [mdp()].
#'  If `res_mdp` was a `dirichletprocess` object, returns another
#'  `dirichletprocess` object with new components and altered weights.
polya <- function(res_mdp, epsilon = 0.01, upsilon = 0.01,
                  nthreads = 1) {
  UseMethod('polya')
}

#' @describeIn polya Polya extension to a `mdpolya_result` object.
polya.mdp <- function(res_mdp, epsilon = 0.01, upsilon = 0.01,
                      nthreads = 1) {
  stopifnot(any(class(res_mdp) == 'mdp') &
              epsilon > 0 & epsilon < 1 &
              upsilon > 0 & upsilon < 1)
  out <- list(phi = NULL, eta = res_mdp$eta, args = res_mdp$args);
  out$phi <- polyaurn_cpp(res_mdp[[1]], res_mdp[[2]],
                          epsilon, upsilon, nthreads)
  sapply(1:length(out$phi), function(i) {
    out$phi[[i]] <<- out$phi[[i]][is.finite(out$phi[[i]][, 2]), , drop = FALSE];
    if (!is.matrix(out$phi[[i]])) { # This may be unnecessary now due to drop
      out$phi[[i]] <<- matrix(out$phi[[i]], nrow = 1)
    }
    colnames(out$phi[[i]]) <<- c('w', 'mean', 'var')
    })
  class(out) <- c('mdpolya_result', 'polya')
  return(out)
}

#' @describeIn polya Polya extension to a `dirichletprocess` object.
polya.dirichletprocess <- function(dp, epsilon = 0.01, upsilon = 0.01,
                                   nthreads = 1) {
  phi <- polyaurngen_cpp(t(sapply(dp$labelsChain, identity)),
                         dp$alphaChain, epsilon, upsilon, nthreads)
  sapply(1:length(phi), function(i) {
    phi[[i]] <<- phi[[i]][is.finite(phi[[i]][, 2]), ];
    if (!is.matrix(phi[[i]])) {
      phi[[i]] <<- matrix(phi[[i]], nrow = 1)
    }
    colnames(phi[[i]]) <<- c('w', 'z')
    phi[[i]][, 'z'] <<- phi[[i]][, 'z'] + 1
    # phi[[i]] <<- phi[[i]][order(phi[[i]][, 'z']), ]
  })
  dp$weightsChain <- lapply(phi, function(p) p[, 'w'])
  dp$clusterParametersChain <- lapply(1:length(phi), function(i) {
    ns <- phi[[i]][, 'z'] > max(dp$labelsChain[[i]])
    if(sum(ns) != 0) {
      pars <- dirichletprocess::PriorDraw(dp$mixingDistribution, sum(ns))
    } else {
      pars <- list(array(dim = c(0, 0, 0)))
    }
    out <- lapply(1:length(dp$clusterParametersChain[[i]]), function(j) {
      if (!all(!ns)) {
        if (prod(dim(pars[[j]])) != 0) {
          abind::abind(dp$clusterParametersChain[[i]][[j]][, ,
                                                           phi[[i]][, 'z'][!ns],
                                                           drop = FALSE],
                       pars[[j]])
        } else {
          dp$clusterParametersChain[[i]][[j]][, , which(!ns),drop = FALSE]
        }
      } else {
        dp$clusterParametersChain[[i]][[j]][, , which(!ns),drop = FALSE]
      }
    })
    return(out)
  })
  return(dp)
}

#' Subset method for `mdpolya_result` objects
#'
#' @param obj An `mdpolya_result` object.
#' @param i A numeric vector of sample indices.
#'
#' @return An `mdpolya_object` with only the samples indexed by `i`.
`[.mdpolya_result` <- function(obj, i) {
  if (!is.null(obj$theta)) {
    obj$theta <- obj$theta[i, , ]
  }
  obj$eta <- obj$eta[i, ]
  obj$phi <- obj$phi[i]
  return(obj)
}

#' Print method for `mdpolya_result` objects
#'
#' @param obj An `mdpolya_result` object.
print.mdpolya_result <- function(obj) {
  obj$args[[1]] <- '...'
  cat(paste0('\nCall: polya(mdp(',
             paste0(paste(names(obj$args),
                          obj$args,
                          sep = ' = '), collapse = ', '),
             '))\n\n'))
  cat('Eta:\n')
  eta <- apply(obj$eta, 2, fivenum)
  rownames(eta) <- c('min', 'low-hinge', 'median','up-hinge', 'max')
  print(t(eta))
  cat('\nPhi:\n')
  phi <- matrix(
    fivenum(sapply(obj$phi, nrow)
    )
  )
  colnames(phi) <- 'z'
  rownames(phi) <- c('min', 'low-hinge', 'median','up-hinge', 'max')
  print(t(phi))
  cat('\n')
}

#' Grid evaluation of `mdpolya_result` values
#'
#' @param obj An `mdpolya_result` object.
#' @param grd A numeric vector of `m` grid points.
#' @param func Either 'distribution', 'density', or 'gradient'. The function of
#'  the mixture distribution to be evaluated.
#' @param nthreads The number of parallel threads to launch with OpenMP.
#'
#' @return An `mdpolya_result_gridded` object, which is a matrix with dimension
#'  [`k`, `m`] of evaluated sample functions, with the following attributes:
#'  * `func`: The evaluated function.
#'  * `grid`: The grid points on which each of the `k` rows was evaluated.
#'  * `args`: A copy of the `args` entry from `obj`.
gridify <- function(obj, grd = NULL, func = 'density',
                                   nthreads = 1) {
  if (func == 'distribution') {
    f <- 0
  } else if (func == 'density') {
    f <- 1
  } else if (func == 'gradient') {
    f <- 2
  } else {
    stop('Unrecognized `func`. Choose either \'density\' or \'distribution\'.')
  }
  if (is.null(grd)) {
    r_x <- range(obj$args$data)
    rr_x <- diff(r_x) / 10
    grd <- seq(r_x[1] - rr_x, r_x[2] + rr_x, length = 1000)
  }
  out <- evalmdpolya_cpp(obj$phi, grd, f, nthreads)
  attr(out, 'func') <- func
  attr(out, 'grid') <- grd
  attr(out, 'args') <- obj$args
  class(out) <- 'mdpolya_result_gridded'
  return(out)
}

#' Attribute access method for `mdpolya_result_gridded` objects
#'
#' @param obj An `mdpolya_result_gridded` object.
#' @param name The name of the attribute to access (i.e. `func`, `grid`, or
#'  `args`).
#'
#' @return The value of the attribute in `obj`.
`$.mdpolya_result_gridded` <- function(obj, name) {
  attr(obj, name)
}

#' Subset method for `mdpolya_result_gridded` objects
#'
#' @param obj An `mdpolya_result_gridded` object.
#' @param i A numeric vector of sample indices.
#'
#' @return An `mdpolya_object_gridded` with only the samples indexed by `i`.
`[.mdpolya_result_gridded` <- function(obj, i) {
  out <- matrix(obj, dim(obj))[i, ]
  class(out) <- class(obj)
  attr(out, 'func') <- attr(obj, 'func')
  attr(out, 'grid') <- attr(obj, 'grid')
  attr(out, 'args') <- attr(obj, 'args')
  attr(out, 'args')$k <- length(i)
  return(out)
}

#' Plotting method for `mdpolya_result` objects
#'
#' @inheritParams gridify
#' @param confint A decimal value indicating the confidence interval width (e.g.
#'  0.95 for a 95% confidence interval). Defaults to `NULL`, in which case no
#'  confidence intervals will be drawn.
#'
#' @return A `ggplot` object.
plot.mdpolya_result <- function(obj, grd = NULL, func = 'density',
                                confint = NULL, nthreads = 1) {
  plot(gridify(obj, grd, func, nthreads), confint)
}

#' Plotting method for `mdpolya_result_gridded` objects
#'
#' @param obj An `mdpolya_result` object.
#' @param confint A decimal value indicating the confidence interval width (e.g.
#'  0.95 for a 95% confidence interval). Defaults to `NULL`, in which case no
#'  confidence intervals will be drawn.
#'
#' @return A `ggplot` object.
plot.mdpolya_result_gridded <- function(obj, confint = NULL) {
  grd <- obj$grid
  df <- data.frame(Value = rep(grd, each = nrow(obj)),
                   K = rep(1:nrow(obj), length(grd)),
                   X = as.numeric(obj))
  p <- ggplot(df, aes(x = Value, y = X, group = K)) +
    ylab(paste0(toupper(substring(obj$func, 1,1)), substring(obj$func, 2))) +
    geom_line(alpha = 0.25, color = 'grey') +
    geom_line(data = data.frame(Value = grd, X = apply(obj, 2, mean)),
              aes(group = 0),
              color = 'red') +
    theme_bw()
  data <- sort(obj$args$data)
  n <- length(data)
  if(obj$func == 'density') {
    p <- p + geom_point(data = data.frame(Value = data +
                                            runif(n, -0.001, 0.001),
                                          X = runif(n, -max(df$X) / 50, 0),
                                          K = 0),
                       shape = 16, size = 0.5, alpha = 0.5)
  } else if (obj$func == 'distribution') {
    e_cdf <- ecdf(data)
    p <- p + stat_function(fun = e_cdf, aes(group = 0), geom = 'step', n = 1001)
    if (!is.null(confint)) {
      err_int <- 1 - confint
      eps_dkw <- sqrt(log(2 / err_int) / (2 * n))
      upper_dkw <- stepfun(data, c(0, pmin(e_cdf(data) + eps_dkw, 1)))
      lower_dkw <- stepfun(data, c(0, pmax(e_cdf(data) - eps_dkw, 0)))
      eps_clt <- qnorm(1 - (err_int / 2)) *
        sqrt(e_cdf(data) * (1 - e_cdf(data)) / n)
      upper_clt <- stepfun(data, c(0, e_cdf(data) + eps_clt))
      lower_clt <- stepfun(data, c(0, e_cdf(data) - eps_clt))
      p <- p +
        stat_function(fun = upper_clt, aes(group = 0), geom = 'step', n = 1001,
                      size = 0.25, alpha = 0.5) +
        stat_function(fun = lower_clt, aes(group = 0), geom = 'step', n = 1001,
                      size = 0.25, alpha = 0.5) +
        stat_function(fun = lower_dkw, aes(group = 0), geom = 'step', n = 1001,
                      size = 0.25, color = 'grey50', linetype = 'longdash',
                      alpha = 0.5) +
        stat_function(fun = upper_dkw, aes(group = 0), geom = 'step', n = 1001,
                      size = 0.25, color = 'grey50', linetype = 'longdash',
                      alpha = 0.5)
    }
  }
  return(p)
}

#' Mode-counting generic function
#'
#' @param obj An object on which to count modes.
#' @param grd A numeric vector over an interval in which to conduct the mode
#'  search.
#'
#' @return A vector of values in `grd` corresponding to a mode.
modes <- function(obj, grd = NULL) {
  UseMethod('modes')
}

#' @describeIn modes Mode-counting method for `mdpolya_result` objects.
modes.mdpolya_result <- function(obj, grd = NULL) {
  modes(gridify(obj, grd = grd, func = 'density'))
}

#' @describeIn modes Mode-counting method for `mdpolya_result_gridded` objects.
modes.mdpolya_result_gridded <- function(obj, grd = NULL) {
  return(apply(obj, 1, function(p) obj$grid[pracma::findpeaks(p)[, 2]]))
}

#' Moment calculation generic function
#'
#' @param obj The object for which a moment will be calculated.
#' @param mom A numeric scalar indicating the moment to calculate.
#' @param cntrl A logical value indicating whether the moment should be central
#'  or not. Defaults to `TRUE`.
#' @param grd A numeric vector of grid values on which the density function
#'  samples in `obj` should be calculated for trapezoidal integration.
#'
#' @return A vector of moment values for each sampled distribution in `obj`.
moments <- function(obj, mom, cntrl = TRUE, grd = NULL) {
  UseMethod('moments')
}

#' @describeIn moments Moment calculation method for `mdpolya_result` objects.
moments.mdpolya_result <- function(obj, mom, cntrl = TRUE, grd = NULL) {
  moments(gridify(obj, grd = grd), mom = mom, cntrl = cntrl)
}

#' @describeIn moments Moment calculation method for `mdpolya_result_gridded`
#'  objects.
moments.mdpolya_result_gridded <- function(obj, mom, cntrl = TRUE, grd = NULL) {
  if (obj$func != 'density') {
    stop('`obj` is not a density evaluation, regridify the `mdpolya_result`.')
  }
  x <- obj$grid
  x_cntr <- x
  if (cntrl & (mom != 1)) {
    x_cntr <- obj$grid - apply(obj, 1, function(f)
      pracma::trapz(obj$grid, obj$grid * f))
  }
  return(apply(obj, 1, function(f) pracma::trapz(x, x_cntr ^ mom * f)))
}
