# TODO: Functional methods, mode finders, replace _gridded classes with one
#  mdpolya class since they act the same anyway, just add the epsilon argument
#  to polya

require(Rcpp)
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

polya <- function(obj, epsilon) {
  stopifnot(any(class(obj) == 'mdp'))
  out <- list(phi = NULL, eta = obj$eta, args = obj$args);
  out$phi <- polyaurn_cpp(obj, epsilon)
  sapply(1:length(out$phi), function(i)
    colnames(out$phi[[i]]) <<- c('w', 'mean', 'var')
    )
  class(out) <- c('mdpolya_result', 'polya')
  return(out)
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

gridify <- function(obj, grd = NULL, func = 'density') {
  UseMethod('gridify')
}

gridify.mdpolya_result <- function(obj, grd = NULL, func = 'density') {
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
  out <- evalmdpolya_cpp(obj$phi, grd, f)
  attr(out, 'func') <- func
  attr(out, 'grid') <- grd
  attr(out, 'args') <- obj$args
  class(out) <- 'mdpolya_result_gridded'
  return(out)
}

`$.mdpolya_result_gridded` <- function(obj, name) {
  attr(obj, name)
}

plot.mdpolya_result <- function(obj, grd = NULL, func = 'density') {
  plot(gridify(obj, grd, func))
}

plot.mdpolya_result_gridded <- function(obj, ...) {
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
  if(obj$func == 'density') {
    n <- length(obj$args$data)
    p <- p + geom_line(data = data.frame(Value = rep(obj$args$data, 2),
                                         X = rep(c(0, 1 / n), each = n),
                                         K = rep(1:n, 2)),
                       alpha = 0.5)
  } else if (obj$func == 'distribution') {
    p <- p + geom_function(fun = ecdf(obj$args$data),
                           aes(group = 0))
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

