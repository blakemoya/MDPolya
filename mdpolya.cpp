#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;


vec rG(double s, double S, double mu, double tau) {
  vec theta(2);
  theta(1) = 1.0 / R::rgamma(s / 2.0, 2.0 / S);
  theta(0) = R::rnorm(mu, sqrt(tau));
  return theta;
}

vec rG(double y, double s, double S, double mu, double tau) {
  vec theta(2);
  theta(1) = 1.0 / R::rgamma((1.0 + s) / 2.0,
        2.0 / (S + pow(y - mu, 2.0) / (1.0 + tau)));
  double x = (mu + tau * y) / (1 + tau);
  double X = tau / (1 + tau);
  theta(0) = R::rnorm(x, sqrt(X * theta(1)));
  return theta;
}

double q(double y, double alpha, double s, double S, double mu, double tau) {
  double q0 = alpha * tgamma((1.0 + s) / 2.0) / (tgamma(s / 2.0) * sqrt(s));
  q0 *= pow(1.0 + pow(y - mu, 2.0) / ((1 + tau) * S), -(1 + s) / 2.0);
  q0 /= sqrt((1 + tau) * S / s);
  return q0;
}

double q(double y, double mu, double v) {
  return exp(-pow(y - mu, 2.0) / (2 * v)) / sqrt(2 * v);
}

// [[Rcpp::export]]
Rcpp::List mdp_cpp(vec y, uword k,
                   double alpha, double mu, double tau, double s, double S,
                   double c, double C, double a, double A, double w, double W,
                   bool fix_a, bool fix_m, bool fix_t,
                   uword burn, uword thin) {
  uword n = y.n_elem;
  cube recs(k, n, 3);
  mat hprecs(k, 5);
  int rec_cnt = 0;
  uvec z(n, fill::ones);
  z = cumsum(z) - 1;
  uvec seq = z;
  mat theta(n, 2);
  for (uword nn = 0; nn < n; nn++) {
    theta.row(nn) = rG(y(nn), s, S, mu, tau).t();
  }
  for (uword iter = 0; iter < burn + thin * k; iter++) {
    uvec z_uq = unique(z);
    mat theta_uq = theta.rows(z_uq);
    uword n_uq = z_uq.n_elem;
    if (!fix_m) {
      double x = A / (A + tau * sum(1 / theta_uq.col(1)));
      double a_post = (1 - x) * a + x * sum(theta_uq.col(1))
        * sum(theta_uq.col(0) / theta_uq.col(1));
      double A_post = x * tau * sum(theta_uq.col(1));
      mu = R::rnorm(a_post, sqrt(A_post));
    }
    if (!fix_t) {
      double w_post = w + n_uq;
      double W_post = W + sum(pow(theta_uq.col(0) - mu, 2.0) / theta_uq.col(1));
      tau = 1.0 / R::rgamma(w_post / 2.0, 2.0 / W_post);
    }
    for (uword i = 0; i < n; i++) {
      double q0 = q(y(i), alpha, s, S, mu, tau);
      vec qj(n, fill::zeros);
      for (uword j = 0; j < n; j++) {
        if (j == i) {
          qj(j) = q0;
        } else {
          qj(j) = q(y(i), theta(z(j), 0), theta(z(j), 1));
        }
      }
      qj = qj / sum(qj);
      uword z_prop = Rcpp::RcppArmadillo::sample(seq, 1, true, qj)(0);
      if (z_prop == i) {
        theta.row(i) = rG(y(i), s, S, mu, tau).t();
        z(i) = z_prop;
      } else {
        z(i) = z(z_prop);
      }
    }
    if (!fix_a) {
      double eta = R::rbeta(alpha + 1, n);
      double u = R::runif(0, 1);
      if (u < (c + n_uq - 1) / (n * (C - log(eta)) + c + n_uq - 1)) {
        alpha = R::rgamma(c + n_uq, 1.0 / (C - log(eta)));
      } else {
        alpha = R::rgamma(c + n_uq - 1, 1.0 / (C - log(eta)));
      }
    }
    if ((iter >= burn) && (iter % thin == 0)) {
      recs.slice(0).row(rec_cnt) = conv_to<vec>::from(z + 1).t();
      recs.slice(1).row(rec_cnt) = theta.col(0).t();
      recs.slice(2).row(rec_cnt) = theta.col(1).t();
      hprecs.row(rec_cnt) = vec({alpha, mu, tau, s, S}).t();
      rec_cnt++;
    }
  }
  return Rcpp::List::create(recs, hprecs);
}

// [[Rcpp::export]]
Rcpp::List polyaurn_cpp(Rcpp::List mdp_res, double eps) {
  Rcpp::List out;
  cube theta = mdp_res[0];
  mat hyper = mdp_res[1];
  uword k = hyper.n_rows;
  uword n = theta.n_cols;
  for (uword kk = 0; kk < k; kk++) {
    double alpha = hyper(kk, 0);
    uword m = 1.0 + R::qpois(0.95, -(alpha + n) * std::log(eps), true, false);
    mat phi(m, 3, fill::zeros);
    vec v = Rcpp::rbeta(m, 1, alpha + n);
    phi.col(0) = v;
    v = cumprod(1 - v);
    for (uword mm = 0; mm < m; mm++) {
      if (mm > 0) {
        phi(mm, 0) *= v(mm - 1);
      }
      double u = R::runif(0, 1) * (alpha + n);
      if (u <= alpha) {
        phi(mm, 1) = R::rnorm(hyper(kk, 1), std::sqrt(hyper(kk, 2)));
        phi(mm, 2) = R::rgamma(hyper(kk, 3) / 2.0, 2.0 / hyper(kk, 4));
      } else {
        int i = std::floor(R::runif(0, 1) * n);
        int z = theta(kk, i, 0) - 1;
        phi(mm, 1) = theta(kk, z, 1);
        phi(mm, 2) = theta(kk, z, 2);
      }
    }
    rowvec phi_ex(3);
    phi_ex(0) = 1.0 - sum(phi.col(0));
    phi_ex(1) = R::rnorm(hyper(kk, 1), std::sqrt(hyper(kk, 2)));
    phi_ex(2) = R::rgamma(hyper(kk, 3) / 2.0, 2.0 / hyper(kk, 4));
    phi = join_cols(phi, phi_ex);
    out.push_back(phi);
  }
  return out;
}

double grad(double x, double mean, double var) {
  return R::dnorm(x, mean, var, false) * (mean - x) / var;
}

double hess(double x, double mean, double var) {
  return R::dnorm(x, mean, var, false) / var
  * (std::pow(mean - x, 2) / var - 1);
}

// [[Rcpp::export]]
mat evalmdpolya_cpp(Rcpp::List mdpolya, vec x, int f) {
  mat out(mdpolya.length(), x.n_elem, fill::zeros);
  for (uword kk = 0; kk < mdpolya.length(); kk++) {
    mat thetak = mdpolya[kk];
    for (uword zz = 0; zz < thetak.n_rows; zz++) {
      for (uword xx = 0; xx < x.n_elem; xx++) {
        if (f == 0) {
          out(kk, xx) += thetak(zz, 0) * R::pnorm(x(xx), thetak(zz, 1),
              std::sqrt(thetak(zz, 2)), true, false);
        } else if (f == 1) {
          out(kk, xx) += thetak(zz, 0) * R::dnorm(x(xx), thetak(zz, 1),
              std::sqrt(thetak(zz, 2)), false);
        } else if (f == 2) {
          out(kk, xx) += thetak(zz, 0) * grad(x(xx), thetak(zz, 1),
              thetak(zz, 2));
        }
      }
    }
  }
  return(out);
}
