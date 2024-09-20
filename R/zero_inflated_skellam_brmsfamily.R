# log likelihood
log_dskellamzi_loo <- function(y, mu, delta, zi, n_draws) {
  require(matrixStats)
  require(Bessel)
  va <- abs(mu)+delta
  mu1 <- (mu+va)/2
  mu2 <- (va-mu)/2
  output <- vector("double", length(y)*n_draws)
  for(i in seq_len(length(y))) {
    if(y[i] == 0) {
      output[(((i-1)*n_draws)+1):(i*n_draws)] <- 
        logSumExp(lx = c(log(zi), log(1-zi)+(-mu1-mu2)+
                           log(besselI(2*sqrt(mu1*mu2), 0))))
    } else {
      output[(((i-1)*n_draws)+1):(i*n_draws)] <- log(1-zi)+(-mu1-mu2)+
        (y[i]/2)*(log(mu1)-log(mu2))+
        besselI.nuAsym(2*sqrt(mu1*mu2), abs(y[i]), k.max = 5, log = TRUE)
    }
  }
  return(output)
}

# rng
rskellamzi <- function(n, mu, delta, zi) {
  va = abs(mu)+delta
  mu1 = (mu+va)/2
  mu2 = (va-mu)/2
  n_zero = n-sum(rbinom(n = n,size = 1, prob = (1-zi)))
  zeros <- rep(x = 0, times = n_zero)
  non_zeros <- rpois(n-n_zero, mu1)-rpois(n-n_zero, mu2)
  all <- c(zeros, non_zeros)
  return(sample(all))
}

# brms family
zero_inflated_skellam <- function(
    link = "identity", link_delta = "log", link_zi = "logit") {
  require(brms)
  family <- custom_family(
    name = "zero_inflated_skellam",
    dpars = c("mu", "delta", "zi"),
    links = c(link, link_delta, link_zi),
    lb = c(NA, 0, 0),
    ub = c(NA, NA, 1),
    type = "int",
    log_lik = function(i, prep) {
      mu <- get_dpar(prep, "mu", i = i)
      delta <- get_dpar(prep, "delta", i = i)
      zi <- get_dpar(prep, "zi", i = i)
      y <- prep$data$Y[i]
      log_dskellamzi_loo(y, mu, delta, zi, n)
    },
    posterior_predict = function(i, prep, ...) {
      mu <- get_dpar(prep, "mu", i = i)
      delta <- get_dpar(prep, "delta", i = i)
      zi <- get_dpar(prep, "zi", i = i)
      n <- prep$ndraws
      rskellamzi(n, mu, delta, zi)
    },
    posterior_epred = function(prep) {
      mu <- get_dpar(prep, "mu")
      zi <- get_dpar(prep, "zi")
      return(mu*(1-zi))
    }
  )
  family$stanvars <- stanvar(
    scode = "
  real zero_inflated_skellam_lpmf(int y, real mu, real delta, real zi) {
    real va = abs(mu)+delta;
    real mu1 = (mu+va)/2;
    real mu2 = (va-mu)/2;
    real total = (-mu1-mu2)+(log(mu1)-log(mu2))*y/2;
    real log_prob = total+log_modified_bessel_first_kind(abs(y), 2*sqrt(mu1*mu2));
    if (y == 0) {
      return log_sum_exp(log(zi), log1m(zi) + log_prob);
    } else {
      return log1m(zi) + log_prob;
    }
  }",
    block = "functions"
  )
  return(family)
}
