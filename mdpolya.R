# TODO: Functional methods, mode finders, replace _gridded classes with one
#  mdpolya class since they act the same anyway, just add the epsilon argument
#  to polya

require(Rcpp)
require(parallel)
require(dirichletprocess)
require(abind)
require(pracma)
require(ggplot2)

sourceCpp('mdpolya.cpp')

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

polya <- function(res_mdp, epsilon = 0.01, upsilon = 0.01,
                  nthreads = 1) {
  UseMethod('polya')
}

polya.mdp <- function(res_mdp, epsilon = 0.01, upsilon = 0.01,
                      nthreads = 1) {
  stopifnot(any(class(res_mdp) == 'mdp') &
              epsilon > 0 & epsilon < 1 &
              upsilon > 0 & upsilon < 1)
  out <- list(phi = NULL, eta = res_mdp$eta, args = res_mdp$args);
  out$phi <- polyaurn_cpp(res_mdp[[1]], res_mdp[[2]], epsilon, upsilon, nthreads)
  sapply(1:length(out$phi), function(i) {
    out$phi[[i]] <<- out$phi[[i]][is.finite(out$phi[[i]][, 2]), , drop = FALSE];
    if (!is.matrix(out$phi[[i]])) { # This may be unnecessary no due to drop
      out$phi[[i]] <<- matrix(out$phi[[i]], nrow = 1)
    }
    colnames(out$phi[[i]]) <<- c('w', 'mean', 'var')
    })
  class(out) <- c('mdpolya_result', 'polya')
  return(out)
}

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

print.mdp_result <- function(obj) {
  obj$args[[1]] <- '...'
  cat(paste0('\nCall: mdp(',
             paste0(paste(names(obj$args),
                          obj$args,
                          sep = ' = '), collapse = ', '),
             ')\n\n'))
  cat('Eta:\n')
  eta <- apply(obj$eta, 2, fivenum)
  rownames(eta) <- c('min', 'low-hinge', 'median','up-hinge', 'max')
  print(t(eta))
  cat('\nTheta:\n')
  theta <- matrix(
    fivenum(apply(obj$theta[, , 1], 1,
                  function(z) length(unique(z)))
            )
    )
  colnames(theta) <- 'z'
  rownames(theta) <- c('min', 'low-hinge', 'median','up-hinge', 'max')
  print(t(theta))
  cat('\n')
}

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

gridify <- function(obj, grd = NULL, func = 'density',
                    nthreads = 1) {
  UseMethod('gridify')
}

gridify.mdpolya_result <- function(obj, grd = NULL, func = 'density',
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

`$.mdpolya_result_gridded` <- function(obj, name) {
  attr(obj, name)
}

`[.mdpolya_result_gridded` <- function(obj, i) {
  out <- matrix(obj, dim(obj))[i, ]
  class(out) <- class(obj)
  attr(out, 'func') <- attr(obj, 'func')
  attr(out, 'grid') <- attr(obj, 'grid')
  attr(out, 'args') <- attr(obj, 'args')
  attr(out, 'args')$k <- length(i)
  return(out)
}

plot.mdpolya_result <- function(obj, grd = NULL, func = 'density',
                                confint = NULL,
                                nthreads = 1) {
  plot(gridify(obj, grd, func, nthreads), confint)
}

plot.mdpolya_result_gridded <- function(obj, confint = NULL, ...) {
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

modes <- function(obj, grd = NULL) {
  UseMethod('modes')
}

modes.mdpolya_result <- function(obj, grd = NULL) {
  modes(gridify(obj, grd = grd, func = 'density'))
}

modes.mdpolya_result_gridded <- function(obj, grd = NULL) {
  return(apply(obj, 1, function(p) obj$grid[pracma::findpeaks(p)[, 2]]))
}

moments <- function(obj, mom, cntrl = TRUE, grd = NULL) {
  UseMethod('moments')
}

moments.mdpolya_result <- function(obj, mom, cntrl = TRUE, grd = NULL) {
  moments(gridify(obj, grd = grd), mom = mom, cntrl = cntrl)
}

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
