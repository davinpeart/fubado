# skellam-normal mixture with random intercepts
fit_model <-
  function(data, dv, iv, id = NULL, full_iv = NULL, predict_ancillary = FALSE,
           family, priors, stanargs = list(
             chains = 4, warmup = 5000, iter = 10000, cores = 2, refresh = 1,
             seed = runif(n = 1, min = 1, max = 99999), control = list(
               adapt_delta = .99, stepsize = .5))) {
  require(rstan)
  # model matrices with treatment contrasts
  if(is.null(full_iv)) {
    full_iv <- iv
  }
  for(i in 1:length(full_iv)) {
    attr(data[[full_iv[i]]], "contrasts") <- contr.treatment(levels(data[[full_iv[i]]]))
  }
  full_needed <- !identical(iv, full_iv) & !is.null(full_iv)
  design <- model.matrix(formula(paste0(
    dv, "~", "1+", paste0(iv, collapse = "*"))), data)
  if(full_needed) {
    design_full <- model.matrix(formula(paste0(
      dv, "~", "1+", paste0(full_iv, collapse = "*"))), data)
  }

  # generate code
  family_list <- family(predict_ancillary = predict_ancillary,
                        repeated_measures = !is.null(id),
                        priors = priors)
  stan_code <- paste0(
    "functions {
    ", family_list[["functions"]], "
    }
    data {
    ", family_list[["data"]], "
    }
    parameters {
    ", family_list[["parameters"]], "
    }
    transformed parameters {
    ", family_list[["transformed_parameters"]], "
    }
    model {
    ", family_list[["priors"]], "
    ", family_list[["likelihood"]], "
    }
    generated quantities {
    ", family_list[["generated_quantities"]], "
    }"
  )

  stan_data <- list(
    N = length(data[[dv]]),
    Y = data[[dv]],
    K = ncol(design),
    X = design,
    I = as.integer(data[[id]]),
    N_I = length(levels(data[[id]])),
    J = family_list[["npar"]]
  )
  if(full_needed) {
    stan_data$K_full <- ncol(design_full)
    stan_data$X_full <- design_full
  }

  # fit stan model
  rstan::sampling(
    object = rstan::stan_model(
      model_code = stan_code),
    data = stan_data,                  # named list of data
    chains = stanargs[["chains"]],    # number of Markov chains
    warmup = stanargs[["warmup"]],     # number of warm up iterations per chain
    iter = stanargs[["iter"]],         # total number of iterations per chain
    cores = stanargs[["cores"]],       # number of cores (one per chain)
    refresh = stanargs[["refresh"]],   # progress shown
    seed = stanargs[["seed"]],         # set seed
    control = stanargs[["control"]]    # tune sampler
  )
  }
