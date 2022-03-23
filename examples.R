require(MASS)
source('mdpolya.R')

gal <- sample(MASS::galaxies / 1000)
k <- 100
eps <- 0.05

res_mdp <- mdp(gal, k)
plot(res_mdp, func = 'density') +
  ggtitle('MDP: Galaxies density function')
ggsave('figs/mdp_pdf.png')
plot(res_mdp, func = 'distribution') +
  ggtitle('MDP: Galaxies distribution function')
ggsave('figs/mdp_cdf.png')

res_polya <- polya(res_mdp, eps)
plot(res_polya, func = 'density') +
  ggtitle('MDP + Polya Urn: Galaxies density function')
ggsave('figs/polya_pdf.png')
plot(res_polya, func = 'distribution') +
  ggtitle('MDP + Polya Urn: Galaxies distribution function')
ggsave('figs/polya_cdf.png')

nm_mdp <- sapply(modes(res_mdp), length)
nm_polya <- sapply(modes(res_polya), length)
nm_min <- min(c(nm_mdp, nm_polya))
nm_max <- max(c(nm_mdp, nm_polya))
qplot(nm_mdp, xmin = nm_min, xmax = nm_max, binwidth = 1) +
  theme_bw() +
  scale_x_continuous(breaks = nm_min:nm_max) +
  ggtitle('MDP: Distribution of the number of modes')
ggsave('figs/mdp_nmodes.png')
qplot(nm_polya, xmin = nm_min, xmax = nm_max, binwidth = 1) +
  theme_bw() +
  scale_x_continuous(breaks = nm_min:nm_max) +
  ggtitle('MDP + Polya Urn: Distribution of the number of modes')
ggsave('figs/polya_nmodes.png')
nk_mdp <- apply(res_mdp$theta, 1, function(k) length(unique(k[, 1])))
qplot(nk_mdp, binwidth = 1) +
  theme_bw() +
  scale_x_continuous(breaks = min(nk_mdp):max(nk_mdp)) +
  ggtitle('MDP: Distribution of the number unique components')
ggsave('figs/mdp_ncomps.png')
