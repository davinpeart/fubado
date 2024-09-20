# skellam-normal mixture with random intercepts and slopes
fit_rm_skellam_normal <- function(data, dv, design, id) {
  require(rstan)
  sampling(
    object = stan_model(
      model_code = "
functions {
  // mean & variance parameterization of skellam lpmf
    real skellam_lpmf(int y, real mu, real sigma) {
      if(mu == 0) {
        return(- sigma + log_modified_bessel_first_kind(abs(y), sqrt(sigma^2)));
      } else {
        return(- sigma + ((log(sigma + mu) - log(sigma - mu)) * y / 2) +
          log_modified_bessel_first_kind(abs(y), sqrt(sigma^2 - mu^2)));
      }
    }

  // skellam rng
    int skellam_rng(real mu, real sigma) {
      return(poisson_rng((sigma + mu) / 2) - poisson_rng((sigma - mu) / 2));
    }
}
data {
  // observations and fixed effects
    int<lower=1> N;  // number of observations
    int Y[N];  // response variable
    int<lower=1> K;  // number of fixed effects
    matrix[N, K] X;  // fixed effect design matrix

  // random effects
    int<lower=1> N_I;  // number of subjects
    int<lower=1,upper=N_I> I[N];  // subject identifier
    int<lower=1> J;  // number of distributional parameters with random intercepts
}
parameters {
  // fixed effects
    vector[K] beta;  // non-standardized regression coefficients
    real phi;  // log positive real dispersion parameter = standard deviation of normal
    real delta;  // log positive real difference between mean and variance of skellam

  // random effects
    matrix[J, N_I] z_I;  // standardized subject intercepts and slopes
    vector<lower=0>[J] sigma_I;  // sd for subject intercepts and slopes
    cholesky_factor_corr[J] L_I;  // correlation matrix for subject intercepts and slopes
    vector[N] mu;  // mixture mean
}
transformed parameters {
  // random effects
    matrix[J, N_I] z; // non-centered subject intercepts and slopes
    z = diag_pre_multiply(sigma_I, L_I) * z_I;
}
model {
  // priors
    //  fixed effects
      beta ~ normal(0, 3);  // positive prior for generalization coefficients
      phi ~ normal(0, 5);  // weakly informative prior for ancillary dispersion
      delta ~ normal(0, 3);  // weakly informative prior for ancillary difference

    // random effects
      L_I ~ lkj_corr_cholesky(1);  // uniform lkj prior on cholesky factors
      sigma_I ~ normal(0, 2.5);  // half-normal on subject sd
      to_vector(z_I) ~ std_normal();  // standard normal on standardized effects

  // mixture
    for(n in 1:N) {  // normal hyperprior on means
      target += normal_lpdf(mu[n] | X[n, ] * beta + z[1, I[n]], exp(phi + z[2, I[n]]));
    }
    for(n in 1:N) {  // skellam likelihood on means
      target += skellam_lpmf(Y[n] | mu[n], abs(mu[n]) + exp(delta + z[3, I[n]]));
    }
}
generated quantities {
  // recover omega
    matrix[J, J] omega;
    omega = multiply_lower_tri_self_transpose(L_I);

  // store random intercepts as array of vectors
    array[J] vector[N_I] z_array;
    for(j in 1:J) {
      z_array[j] = to_vector(z[j, ]);
    }

  // replications for posterior predictive checks
    array[J] vector[N_I] z_rep;  // random intercept replications
    z_rep = multi_normal_rng(z_array, quad_form_diag(omega, sigma_I));
    array[N] real mu_rep; // mu replications
    for(n in 1:N) {
      mu_rep[n] = normal_rng(X[n, ] * beta + z_rep[1, I[n]], exp(phi + z_rep[2, I[n]]));
    }
    array[N] int y_rep; // Y replications
    for(n in 1:N) {
      y_rep[n] = skellam_rng(mu[n], abs(mu[n]) + exp(delta + z[3, I[n]]));
    }
}"),
    data = list(
      N = length(data[[dv]]),
      Y = data[[dv]],
      K = ncol(design),
      X = design,
      I = as.integer(data[[id]]),
      N_I = length(levels(data[[id]])),
      J = 3
    ),                      # named list of data
    chains = 2,             # number of Markov chains
    warmup = 4000,          # number of warm up iterations per chain
    iter = 8000,            # total number of iterations per chain
    cores = 2,              # number of cores (one per chain)
    refresh = 1,            # progress shown
    seed = 837085670,       # set seed
    control = list(         # adjust sampler
      adapt_delta = .9,
      stepsize = .7
    )
    )}

enframe_prop_integer <- function(y) {
  if(is.vector(y)) {
    x <- table(y)
    z <- as.vector(unname(x))
    return(data.frame(integer = as.integer(names(x)), freq = z, prop = z/length(y)))
  }
  if(is.array(y)) {
    I <- dim(y)[1]
    mx <- vector("numeric", I)
    mn <- vector("numeric", I)
    x <- vector("list", I)
    nm <- vector("list", I)
    for(i in 1:I) {
      x[[i]] <- table(y[i, ])/dim(y)[2]
      nm[[i]] <- as.integer(names(x[[i]]))
      mx[i] <- max(nm[[i]])
      mn[i] <- min(nm[[i]])
    }
    z <- seq.int(min(mn), max(mx))
    J <- length(z)
    w <- vector("list", J)
    names(w) <- as.character(z)
    for(j in 1:J) {
      for(i in 1:I) {
        w[[j]][i] <- x[[i]][names(w)[j]]
      }
    }
    op <- data.frame(integer = z)
    op$mean <- 0
    op$lower <- 0
    op$upper <- 0
    for(k in 1:J) {
      op[k, "mean"] <- mean(w[[as.character(op[k, "integer"])]], na.rm = T)
      op[k, "lower"] <- unname(quantile(w[[as.character(op[k, "integer"])]], .05, na.rm = T))
      op[k, "upper"] <- unname(quantile(w[[as.character(op[k, "integer"])]], .95, na.rm = T))
    }
    return(op)
  }
}

enframe_descriptives <- function(y, stat1 = mean, stat2 = var) {
  if(is.vector(y)) {
    return(data.frame(mean = stat1(y), var = stat2(y)))
  }
  if(is.matrix(y)) {
    I <- dim(y)[1]
    s1 <- vector("numeric", I)
    s2 <- vector("numeric", I)
    for(i in 1:I) {
      s1[i] <- stat1(y[i, ])
      s2[i] <- stat2(y[i, ])
    }
    op <- data.frame(s1, s2)
    colnames(op) <- c(as.character(substitute(stat1)), as.character(substitute(stat2)))
    return(op)
  }
}

set_theme <- function() {
  require(tidyverse)
  require(ggdist)
  theme_ggdist() +
    theme(text = element_text(family = "Times New Roman", colour = "gray30", size = 8),
          legend.text = element_text(family = "Times New Roman", colour = "gray30", size = 8),
          legend.title = element_text(family = "Times New Roman", colour = "gray30", size = 8),
          strip.background = element_rect(fill="white"),
          plot.title = element_text(hjust = 0.5))
}

pp_int <- function(y, yrep, ynew, ynew2) {
  require(tidyverse)
  ggplot(enframe_prop_integer(y), aes(x = integer)) +
    geom_bar(stat = "identity", mapping = aes(y = prop), fill = "#FF88A5") +
    geom_line(data = enframe_prop_integer(ynew), aes(y = prop), colour = "gray55") +
    geom_errorbar(data = enframe_prop_integer(yrep),
                  mapping = aes(ymin = lower, ymax = upper), colour = "#687EC9", width = 0) +
    geom_point(data = enframe_prop_integer(yrep), mapping = aes(y = mean), colour = "#687EC9") +
    geom_point(data = enframe_prop_integer(ynew2), mapping = aes(y = prop),
               shape = 21, fill = NA, colour = "gray55") +
    set_theme() + xlim(-10.5, 15.5)
}

pp_stat_dens <- function(yrep, y, ynu1, ynu2) {
  require(tidyverse)
  ramp_pallette <- colorRampPalette(c("white", "#687EC9"))
  ggplot(enframe_descriptives(yrep), aes(x = mean, y = var)) +
    geom_density2d_filled(contour_var = "ndensity", bins = 5,
                          linewidth = 0, alpha = .8, adjust = 1) +
    scale_fill_manual(values = c(ramp_pallette(5))) +
    geom_segment(data = enframe_descriptives(y),
                 aes(x = mean, xend = mean, y = -Inf,  yend = var), colour = "#D86A83", linetype = 2) +
    geom_segment(data = enframe_descriptives(y),
                 aes(x = -Inf, xend = mean, y = var,  yend = var), colour = "#D86A83", linetype = 2) +
    geom_point(data = enframe_descriptives(y), colour = "#D86A83", fill = "#FF88A5",
               alpha = 1, shape = 21, size = 2.5) +
    geom_segment(data = enframe_descriptives(ynu1),
                 aes(x = mean, xend = mean, y = -Inf,  yend = var), colour = "gray55", linetype = 3) +
    geom_segment(data = enframe_descriptives(ynu1),
                 aes(x = -Inf, xend = mean, y = var,  yend = var), colour = "gray55", linetype = 3) +
    geom_point(data = enframe_descriptives(ynu1), colour = "gray55", fill = NA,
               alpha = 1, shape = 21, size = 2.5) +
    geom_segment(data = enframe_descriptives(ynu2),
                 aes(x = mean, xend = mean, y = -Inf,  yend = var), colour = "gray55", linetype = 3) +
    geom_segment(data = enframe_descriptives(ynu2),
                 aes(x = -Inf, xend = mean, y = var,  yend = var), colour = "gray55", linetype = 3) +
    geom_point(data = enframe_descriptives(ynu2), colour = "gray55", fill = NA,
               alpha = 1, shape = 21, size = 2.5) +
    set_theme()
}


