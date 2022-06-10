require(MASS)
require(microbenchmark)
require(ggplot2)
source('mdpolya.R')
set.seed(517)

gal <- sample(MASS::galaxies / 1000)
k1 <- 100

res_mdp <- mdp(gal, k1)
plot(res_mdp, func = 'density') +
  ggtitle('MDP: Galaxies density function')
ggsave('figs/mdp_pdf.png', width = 6, height = 4)
plot(res_mdp, func = 'distribution') +
  ggtitle('MDP: Galaxies distribution function')
ggsave('figs/mdp_cdf.png', width = 6, height = 4)
plot(res_mdp, func = 'distribution', confint = 0.95) +
  ggtitle('MDP: Galaxies distribution function')
ggsave('figs/mdp_cdf95.png', width = 6, height = 4)

res_pol <- polya(res_mdp)
plot(res_pol, func = 'density') +
  ggtitle('MDP + Polya Urn: Galaxies density function')
ggsave('figs/polya_pdf.png', width = 6, height = 4)
plot(res_pol, func = 'distribution') +
  ggtitle('MDP + Polya Urn: Galaxies distribution function')
ggsave('figs/polya_cdf.png', width = 6, height = 4)
plot(res_pol, func = 'distribution', confint = 0.95) +
  ggtitle('MDP + Polya Urn: Galaxies distribution function')
ggsave('figs/polya_cdf95.png', width = 6, height = 4)


k2 <- 1000

res_mdp <- mdp(gal, k2)
res_pol <- polya(res_mdp)

mdp_mom <- data.frame(Mean = moments(res_mdp, 1),
                      Variance = moments(res_mdp, 2))
mdp_mom_ <- mdp_mom
while (nrow(mdp_mom_) >= 0.95 * nrow(mdp_mom)) {
  mdp_mom_h <- chull(mdp_mom_)
  mdp_mom_ <- mdp_mom_[-c(mdp_mom_h), ]
}
mdp_mom_h <- mdp_mom_[chull(mdp_mom_), ]
pol_mom <- data.frame(Mean = moments(res_pol, 1),
                      Variance = moments(res_pol, 2))
pol_mom_ <- pol_mom
while (nrow(pol_mom_) >= 0.95 * nrow(pol_mom)) {
  pol_mom_h <- chull(pol_mom_)
  pol_mom_ <- pol_mom_[-c(pol_mom_h), ]
}
pol_mom_h <- pol_mom_[chull(pol_mom_), ]
ggplot(mdp_mom, aes(x = Mean, y = Variance)) +
  geom_point(alpha = 0.5, pch = 16, color = 'deeppink3') +
  geom_polygon(data = mdp_mom_h, color = 'deeppink3', fill = NA) +
  geom_point(data = pol_mom, alpha = 0.5, pch = 16, color = 'darkorange3') +
  geom_polygon(data = pol_mom_h, color = 'darkorange3', fill = NA) +
  geom_point(data = data.frame(Mean = mean(gal),
                               Variance = var(gal)),
             pch = 19, size = 2) +
  theme_bw()

nm_mdp <- sapply(modes(res_mdp), length)
nm_pol <- sapply(modes(res_pol), length)
nm_min <- min(c(nm_mdp, nm_pol))
nm_max <- max(c(nm_mdp, nm_pol))
qplot(nm_mdp, xmin = nm_min, xmax = nm_max, binwidth = 1) +
  theme_bw() +
  xlab('Modes') +
  scale_x_continuous(breaks = nm_min:nm_max) +
  ggtitle('MDP: Distribution of the number of modes')
ggsave('figs/mdp_nmodes.png', width = 6, height = 3)
qplot(nm_pol, xmin = nm_min, xmax = nm_max, binwidth = 1) +
  theme_bw() +
  xlab('Modes') +
  scale_x_continuous(breaks = nm_min:nm_max) +
  ggtitle('MDP + Polya Urn: Distribution of the number of modes')
ggsave('figs/polya_nmodes.png', width = 6, height = 3)
nk_mdp <- sapply(res_mdp$phi, nrow)
nk_pol <- sapply(res_pol$phi, nrow)
qplot(nk_mdp, binwidth = 1) +
  theme_bw() +
  xlab('k') +
  scale_x_continuous(breaks = min(nk_mdp, nk_pol):max(nk_mdp, nk_pol)) +
  ggtitle('MDP: Distribution of the number unique components')
ggsave('figs/mdp_ncomps.png', width = 6, height = 3)
qplot(nk_pol, binwidth = 1) +
  theme_bw() +
  xlab('k') +
  scale_x_continuous(breaks = min(nk_mdp, nk_pol):max(nk_mdp, nk_pol)) +
  ggtitle('MDP + Polya Urn: Distribution of the number unique components')
ggsave('figs/polya_ncomps.png', width = 6, height = 3)


rDPnorm <- function(n, alpha = 1, mu = 21, tau = 25, s = 4, S = 2,
                 c = 2, C = 4, a = 21, A = 21, w = 1, W = 100,
                 fix_a = FALSE, fix_m = FALSE, fix_t = FALSE,
                 eps = 0.05) {
  if (fix_a) {
    alpha <- rep(alpha, n)
  } else {
    alpha <- rgamma(n, c / 2, scale = C / 2)
  }
  if (fix_m) {
    mu <- rep(mu, n)
  } else {
    mu <- rnorm(n, a, sqrt(A))
  }
  if (fix_t) {
    tau <- rep(tau, n)
  } else {
    tau <- rgamma(n, w / 2, scale = W / 2)
  }
  rG0 <- sapply(1:n, function(i)
    function(nn) {
      return(matrix(c(rnorm(nn, mu[i], tau[i]),
                      rgamma(nn, s / 2, scale = S / 2)), ncol = 2))
      })
  M <- sapply(alpha,
              function(alph) 1.0 + qpois(0.95, -(alph) * log(eps)))
  phi <- function(m, alph, rG) {
    v <- rbeta(m, 1, alph)
    w <- v
    if (m > 1) {
      w[2:m] <- v[2:m] * cumprod(1 - v[1:(m - 1)])
    }
    mat <- cbind(w, rG(m))
    mat <- rbind(mat, c(1 - sum(mat[, 1]), rG(1)))
    colnames(mat) <- c('w', 'mean', 'var')
    return(mat)
  }
  return(sapply(1:n, function(nn) phi(M[nn], alpha[nn], rG0[[nn]])))
}

evalDPnorm <- function(obj, grd, func = 'density',
                       nthreads = parallel::detectCores()) {
  if (func == 'distribution') {
    f <- 0
  } else if (func == 'density') {
    f <- 1
  } else if (func == 'gradient') {
    f <- 2
  } else {
    stop('Unrecognized `func`. Choose either \'density\' or \'distribution\'.')
  }
  return(evalmdpolya_cpp(obj, grd, f, nthreads))
}

rmixnorm <- function(n, w, mean, var) {
  k <- sample(1:length(w), n, replace = TRUE, prob = w)
  return(list(y = rnorm(n, mean = mean[k], sd = sqrt(var[k])), k = k))
}

m <- 10
n <- 82
k <- 100
mdps <- rDPnorm(m)
synth <- lapply(mdps, function(x) rmixnorm(n, x[, 1], x[, 2], x[, 3]))
mdps_trunc <- lapply(1:m, function(mm) {
  out <- mdps[[mm]][sort(unique(synth[[mm]]$k)), ]
  if (!is.matrix(out)) {
    out <- t(as.matrix(out))
  }
  out[, 1] <- table(synth[[mm]]$k) / n
  return(out)
})
res_mdp <- lapply(synth, function(s) {
  mdp(s$y, k)
})
res_mdp_g <- lapply(res_mdp, gridify, func = 'distribution')
sapply(1:m, function(mm) dir.create(paste0('figs/', mm)))
for (mm in 1:m) {
  obj <- res_mdp_g[[mm]]
  grd <- obj$grid
  df <- data.frame(Value = rep(grd, 2),
                   K = rep(c('True', 'Truncated'), each = length(grd)),
                   X = c(evalDPnorm(mdps[mm], grd),
                         evalDPnorm(mdps_trunc[mm], grd)))
  p <- ggplot(df, aes(x = Value, y = X, group = K)) +
    ylab('Density') +
    geom_line(data = df[df$K == 'True', ],
              color = 'deepskyblue4') +
    geom_line(data = df[df$K == 'Truncated', ],
              color = 'deepskyblue2',
              linetype = 'dashed') +
    theme_bw() +
    theme(legend.position = 'none')
  n <- length(obj$args$data)
  p <- p + geom_point(data = data.frame(Value = rep(obj$args$data, 2) +
                                          runif(n, -0.001, 0.001),
                                        X = runif(n, -max(df$X) / 50, 0),
                                        K = 0),
                      shape = 16, size = 0.5, alpha = 0.5)
  ggsave(paste0('figs/', mm, '/density.png'), width = 6, height = 4)
  df <- data.frame(Value = rep(grd, 2),
                   K = rep(c('True', 'Truncated'), each = length(grd)),
                   X = c(evalDPnorm(mdps[mm], grd, func = 'distribution'),
                         evalDPnorm(mdps_trunc[mm], grd,
                                    func = 'distribution')))
  p <- ggplot(df, aes(x = Value, y = X, group = K)) +
    ylab('Distribtion') +
    geom_line(data = df[df$K == 'True', ],
              color = 'deepskyblue4') +
    geom_line(data = df[df$K == 'Truncated', ],
              color = 'deepskyblue2',
              linetype = 'dashed') +
    theme_bw() +
    theme(legend.position = 'none')
  p <- p + stat_function(fun = ecdf(obj$args$data), aes(group = 0),
                         geom = 'step', n = 1001)
  err_int <- 1 - 0.95
  eps_dkw <- sqrt(log(2 / err_int) / (2 * n))
  upper_dkw <- approxfun(obj$grid, c(pmin(df$X[1:length(obj$grid)] + eps_dkw, 1)))
  lower_dkw <- approxfun(obj$grid, c(pmax(df$X[1:length(obj$grid)] - eps_dkw, 0)))
  eps_clt <- qnorm(1 - (err_int / 2)) *
    sqrt(df$X[1:length(obj$grid)] * (1 - df$X[1:length(obj$grid)]) / n)
  upper_clt <- approxfun(obj$grid, c(df$X[1:length(obj$grid)] + eps_clt))
  lower_clt <- approxfun(obj$grid, c(df$X[1:length(obj$grid)] - eps_clt))
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
  ggsave(paste0('figs/', mm, '/distribution.png'), width = 6, height = 4)
}
res_pol <- lapply(res_mdp, polya)
res_pol_g <- lapply(res_pol, gridify, func = 'distribution')
for (mm in 1:m) {
  p <- plot(res_mdp_g[[mm]][1:100], confint = 0.95) +
    geom_line(data = data.frame(Value = res_mdp_g[[mm]]$grid,
                                X = c(evalDPnorm(mdps[mm],
                                                 res_mdp_g[[mm]]$grid,
                                                 func = 'distribution')),
                                K = 0),
              color = 'deepskyblue4') +
    geom_line(data = data.frame(Value = res_mdp_g[[mm]]$grid,
                                X = c(evalDPnorm(mdps_trunc[mm],
                                                 res_mdp_g[[mm]]$grid,
                                                 func = 'distribution')),
                                K = 0),
              color = 'deepskyblue2',
              linetype = 'dashed')
  ggsave(paste0('figs/', mm, '/mdp.png'), width = 6, height = 4)
  p <- plot(res_pol_g[[mm]][1:100], confint = 0.95) +
    geom_line(data = data.frame(Value = res_pol_g[[mm]]$grid,
                                X = c(evalDPnorm(mdps[mm],
                                                 res_pol_g[[mm]]$grid,
                                                 func = 'distribution')),
                                K = 0),
              color = 'deepskyblue4') +
    geom_line(data = data.frame(Value = res_pol_g[[mm]]$grid,
                                X = c(evalDPnorm(mdps_trunc[mm],
                                                 res_pol_g[[mm]]$grid,
                                                 func = 'distribution')),
                                K = 0),
              color = 'deepskyblue2',
              linetype = 'dashed')
  ggsave(paste0('figs/', mm, '/pol.png'), width = 6, height = 4)
  mdp_mom <- data.frame(Mean = moments(res_mdp[[mm]], 1),
                        Variance = moments(res_mdp[[mm]], 2))
  mdp_mom_ <- mdp_mom
  while (nrow(mdp_mom_) >= 0.95 * nrow(mdp_mom)) {
    mdp_mom_h <- chull(mdp_mom_)
    mdp_mom_ <- mdp_mom_[-c(mdp_mom_h), ]
  }
  mdp_mom_h <- mdp_mom_[chull(mdp_mom_), ]
  pol_mom <- data.frame(Mean = moments(res_pol[[mm]], 1),
                        Variance = moments(res_pol[[mm]], 2))
  pol_mom_ <- pol_mom
  while (nrow(pol_mom_) >= 0.95 * nrow(pol_mom)) {
    pol_mom_h <- chull(pol_mom_)
    pol_mom_ <- pol_mom_[-c(pol_mom_h), ]
  }
  pol_mom_h <- pol_mom_[chull(pol_mom_), ]
  grd <- res_pol_g[[mm]]$grid
  full <- evalDPnorm(mdps[mm], grd, func = 'density')
  full_m <- pracma::trapz(grd, grd * full)
  full <- data.frame(Mean = full_m,
                     Variance = pracma::trapz(grd, (grd - full_m) ^ 2 * full))
  trnc <- evalDPnorm(mdps_trunc[mm], grd, func = 'density')
  trnc_m <- pracma::trapz(grd, grd * trnc)
  trnc <- data.frame(Mean = trnc_m,
                     Variance = pracma::trapz(grd, (grd - trnc_m) ^ 2 * trnc))
  p <- ggplot(mdp_mom_h, aes(x = Mean, y = Variance)) +
    # geom_point(data = pol_mom, alpha = 0.5, pch = 16, color = 'springgreen2') +
    geom_polygon(data = pol_mom_h, alpha = 0.25, size = 0.5,
                 color = 'springgreen2', fill = 'springgreen2') +
    # geom_point(alpha = 0.5, pch = 16, color = 'springgreen4') +
    geom_polygon(data = mdp_mom_h, alpha = 0.25, size = 0.5,
                 color = 'springgreen4', fill = 'springgreen4') +
    geom_point(data = full, pch = 17, size = 3, col = 'deepskyblue4') +
    geom_point(data = trnc, pch = 18, size = 3, col = 'deepskyblue2') +
    geom_point(data = data.frame(Mean = mean(synth[[mm]]$y),
                                 Variance = var(synth[[mm]]$y)),
               pch = 19, size = 2) +
    theme_bw()
  ggsave(paste0('figs/', mm, '/moms.png'), width = 6, height = 6)
}

mbm_mdp <- microbenchmark("MDP k100n10" = { mdp(sample(gal, 10, replace = TRUE),
                                            100, burn = 0, thin = 1) },
                          "MDP k100n50" = { mdp(sample(gal, 50, replace = TRUE),
                                            100, burn = 0, thin = 1) },
                          "MDP k100n100" = { mdp(sample(gal, 100,
                                                        replace = TRUE),
                                            100, burn = 0, thin = 1) })
autoplot(mbm_mdp) + theme_bw()
ggsave('figs/mbm_mdp.png', width = 12, height = 5)

res_mdp <- mdp(sample(gal, 1000, replace = TRUE), 1000, burn = 0, thin = 1)
mbm_pol <- microbenchmark("Polya e0.05u0.05" = { polya(res_mdp, 0.05, 0.05) },
                          "Polya e0.01u0.05" = { polya(res_mdp, 0.01, 0.05) },
                          "Polya e0.05u0.01" = { polya(res_mdp, 0.05, 0.01) },
                          "Polya e0.01u0.01" = { polya(res_mdp, 0.01, 0.01) })
# mbm_pol <- microbenchmark("Polya e0.05u0.05" = { polya(res_mdp, 0.05, 0.05, nthreads = 1) },
#                           "Polya6e0.05u0.05" = { polya(res_mdp, 0.05, 0.05, nthreads = 6) },
#                           "Polya e0.01u0.01" = { polya(res_mdp, 0.01, 0.01, nthreads = 1) },
#                           "Polya6e0.01u0.01" = { polya(res_mdp, 0.01, 0.01, nthreads = 6) },
#                           times = 50L)
autoplot(mbm_pol) + theme_bw()
ggsave('figs/mbm_pol.png', width = 12, height = 5)
