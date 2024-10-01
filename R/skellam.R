skellam <- function(predict_ancillary, repeated_measures, priors) {
  # functions for printing
  add_fixef <- function(ancillaries, parname, dpar = F) {
    if(parname %in% ancillaries) {
      paste0(
        if(!dpar) {
          "X_full[n, ] * "
        } else {
          "vector[K_full] "
        }, "beta_", parname)
    } else {
      paste0(
        if(dpar) {
          "real "
        },
        parname
      )
    }
  }
  add_ranef <- function(logical, index, rep = F) {
    if(logical) {
      paste0(" + z", if(rep) {
        "_rep"
      }, "[", index, ", I[n]]")
    }
  }

  # functions block
  functions <- "
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
  "

  # data block
  data <- paste0(
    "
// observations and fixed effects
    int<lower=1> N;  // number of observations
    int Y[N];  // response variable
    int<lower=1> K;  // number of fixed effects
    matrix[N, K] X;  // fixed effect design matrix
", if(predict_ancillary) {"
    int<lower=1> K_full;  // number of effects in full design
    matrix[N, K_full] X_full;  // full design matrix"
},
if(repeated_measures) {
  "
  // random effects
    int<lower=1> N_I;  // number of subjects
    int<lower=1,upper=N_I> I[N];  // subject identifier
    int<lower=1> J;  // number of distributional parameters
"
}
  )

  # parameters block
  parameters <- paste0(
    "
// fixed effects
    vector[K] beta;  // regression coefficients
     ", add_fixef(predict_ancillary, "delta", T),
    ";  // log positive real difference between mean and variance of skellam
     ",
    if(repeated_measures) {
      "

// random effects
    matrix[J, N_I] z_I;  // standardized subject intercepts
    vector<lower=0>[J] sigma_I;  // sd for subject intercepts and slopes
    cholesky_factor_corr[J] L_I;  // correlation matrix for subject intercepts and slopes
"
    }
  )

  #print transformed parameters
  transformed_parameters <-
    if(repeated_measures) {
      "
  // random effects
    matrix[J, N_I] z; // non-centered subject intercepts and slopes
    z = diag_pre_multiply(sigma_I, L_I) * z_I;
"
    } else { "
      "}

  # print priors
  priors <- paste0("
// fixed effects
  beta ~ ", priors[["beta"]], ";
  ", if("delta" %in% predict_ancillary) {
    "beta_"
  }, "delta ~ ", priors[["delta"]], ";
  ",
  if(repeated_measures) {
    paste0(
      "

// random effects
  L_I ~ ", priors[["cor"]], ";
  sigma_I ~ ", priors[["sd"]], ";
  to_vector(z_I) ~ std_normal();  // standard normal on standardized effects
"
    )
  }
  )

  # print likelihood
  likelihood <- paste0("
  for(n in 1:N) {  // skellam likelihood
      target += skellam_lpmf(Y[n] | X[n, ] * beta", add_ranef(repeated_measures, 1),
                       "abs(X[n, ] * beta", add_ranef(repeated_measures, 1), ") + exp(",
                       add_fixef(predict_ancillary, "delta"),
                       add_ranef(repeated_measures, 2), "));
    }
  "
  )

  # print generated quantities
  generated_quantities <- paste0(
    if(repeated_measures) {"
  // recover omega
    matrix[J, J] omega;
    omega = multiply_lower_tri_self_transpose(L_I);

  // store random intercepts as array of vectors
    array[J] vector[N_I] z_array;
    for(j in 1:J) {
      z_array[j] = to_vector(z[j, ]);
    }

  // correlated replications of random intercepts
    array[J] vector[N_I] z_rep;  // random intercept replications
    z_rep = multi_normal_rng(z_array, quad_form_diag(omega, sigma_I));

      "
    }, "
  // replications for posterior predictive checks",
    paste0("
  array[N] int y_rep; // Y replications
  for(n in 1:N) {
      y_rep[n] = skellam_rng(X[n, ] * beta", add_ranef(repeated_measures, 1, T),
           "abs(X[n, ] * beta", add_ranef(repeated_measures, 1, T), ") + exp(",
           add_fixef(predict_ancillary, "delta"),
           add_ranef(repeated_measures, 2, T), "));
    }
  "
    )
  )
  return(list(functions = functions, data = data, parameters = parameters,
              transformed_parameters = transformed_parameters,
              priors = priors, likelihood = likelihood,
              generated_quantities = generated_quantities, npar = 3))
}
