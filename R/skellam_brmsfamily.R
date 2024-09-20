# log likelihood
log_dskellam_loo <- function(y, mu, delta, n_draws) {
  require(Bessel)
  va <- abs(mu)+delta
  mu1 <- (mu+va)/2
  mu2 <- (va-mu)/2
  output <- vector("double", length(y)*n_draws)
  for(i in seq_len(length(y))) {
    if(y[i] == 0) {
      output[(((i-1)*n_draws)+1):(i*n_draws)] <- (-mu1-mu2)+
        log(besselI(2*sqrt(mu1*mu2), 0))
    } else {
      output[(((i-1)*n_draws)+1):(i*n_draws)] <- (-mu1-mu2)+(y[i]/2)*
        (log(mu1)-log(mu2))+besselI.nuAsym(2*sqrt(mu1*mu2), abs(y[i]), 
                                           k.max = 5, log = TRUE)
    }
  }
  return(output)
}

# rng
rskellam <- function(n, mu, delta) { 
  va <- abs(mu) + delta
  mu1 <- (mu + va)/2
  mu2 <- (va - mu)/2
  return(rpois(n, mu1) - rpois(n, mu2))
}

# brms family
skellam <- function(link = "identity", link_delta = "log") {
  family <- custom_family(
    name = "skellam",
    dpars = c("mu", "delta"),
    links = c(link, link_delta),
    lb = c(NA, 0),
    ub = c(NA, NA),
    type = "int",
    log_lik = function(i, prep, ...) {
      mu <- get_dpar(prep, "mu", i = i)
      delta <- get_dpar(prep, "delta", i = i)
      n <- prep$ndraws
      y <- prep$data$Y[i]
      log_dskellam_loo(y, mu, delta, n)
    },
    posterior_predict = function(i, prep, ...) {
      mu <- get_dpar(prep, "mu", i = i)
      delta <- get_dpar(prep, "delta", i = i)
      n <- prep$ndraws
      rskellam(n, mu, delta)
    },
    posterior_epred = function(prep, ...) {
      mu <- get_dpar(prep, "mu")
      return(mu)
    }
  )
  family$stanvars <- stanvar(
    scode = "real skellam_lpmf(int y, real mu, real delta) {
    real va = abs(mu)+delta;
    real mu1 = (mu+va)/2;
    real mu2 = (va-mu)/2;
    real total = (-mu1-mu2)+(log(mu1)-log(mu2))*y/2;
    real log_prob = total+log_modified_bessel_first_kind(abs(y), 2*sqrt(mu1*mu2));
    return log_prob;
  }
  ",
    block = "functions"
  )
  return(family)
}
