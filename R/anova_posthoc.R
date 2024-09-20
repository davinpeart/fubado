# vectorized matching for ordered patterns of unique strings
grep_vec_unique <- Vectorize(FUN = grep, vectorize.args = "pattern", USE.NAMES = F)

# regex sufficient for matching nested character vectors
perl_regex_without <- function(vec) {
  require(stringr)
  J <- length(vec)
  out <- vector("character", J)
  for(j in 1:J) {
    if(isFALSE(grepl(x = vec[j], pattern = ":"))) {
      out[j] <- stringr::regex(paste("^(?=.*", vec[j], ")", "(?!.*:)", sep = ""))
    } else { 
      if(isFALSE(grepl(x = vec[j], pattern = stringr::regex("^[^:]*:[^:]*:[^:]*$"), perl = TRUE))) {
        out[j] <- stringr::regex(paste("^(?=.*", vec[j], ")(?!.*^[^:]*:[^:]*:[^:]*$)", sep = ""))
      } else {
        out[j] <- stringr::regex(paste("^(?=.*", vec[j], ")", sep = ""))
      }
    }
  }
  return(out)
}

# regex for matching model interaction terms appearing in different order
reg_eff_matches <- function(vec) {
  J <- length(vec)
  out <- vector("character", J)
  effstr <- vector("list", J)
  regeff <- vector("list", J)
  for(j in 1:J) {
    if(isFALSE(grepl(x = vec[j], pattern = ":"))) {
      out[j] <- stringr::regex(paste("^(?=.*", vec[j], ")", "(?!.*:)", sep = ""))
    } else { 
      if(isFALSE(grepl(x = vec[j], pattern = stringr::regex("^[^:]*:[^:]*:[^:]*$"), perl = TRUE))) {
        effstr[[j]] <- unlist(strsplit(vec[j], split = ":", fixed = TRUE))
        I <- length(effstr[[j]])
        for(i in 1:I)
          regeff[[j]][i] <- paste("^(?=.*", effstr[[j]][i], ")", sep = "")
        out[j] <- stringr::regex(paste0(paste0(
          stringr::regeff[[j]], collapse = ""), "(?!.*^[^:]*:[^:]*:[^:]*$)", 
                               collapse = ""))
      } else {
        effstr[[j]] <- unlist(strsplit(vec[j], split = ":", fixed = TRUE))
        I <- length(effstr[[j]])
        for(i in 1:I)
          regeff[[j]][i] <- paste("^(?=.*", effstr[[j]][i], ")", sep = "")
        out[j] <- stringr::regex(paste0(regeff[[j]], collapse = ""))
      }
    }
  }
  return(out)
}

# apply GG to ez output for non-spherical effects
merge_gg <- function(aov_table) {
  require(tidyverse)
  anova <- aov_table$ANOVA %>%
    select(-ges) %>%
    mutate(p_adj = "none")
  reg.in <- aov_table$`Mauchly's Test for Sphericity`[
    grep(x = aov_table$`Mauchly's Test for Sphericity`$`p<.05`, 
         pattern = "*", fixed = TRUE), "Effect"]
  if(length(reg.in) == 0) {
    return(mutate(anova, n2p = SSn / (SSn + SSd)) %>% relocate(n2p, .after = `F`))
  } else {
    gr <- grep_vec_unique(x = anova$Effect, pattern = perl_regex_without(reg.in), perl = TRUE)
    aov_filter <- anova[-gr,]
    aov_rev <- anova[gr,]
    gg_table <- aov_table$`Sphericity Corrections`[
      grep_vec_unique(x = aov_table$`Sphericity Corrections`$Effect,
                      pattern = perl_regex_without(reg.in), perl = TRUE),] %>%
      select(Effect:`p[GG]<.05`) %>%
      mutate(p_adj = "[GG]") %>%
      rename("p" = "p[GG]",
             "p<.05" = "p[GG]<.05") %>%
      mutate(DFn = use_series(aov_rev, DFn),
             DFd = use_series(aov_rev, DFd),
             SSn = use_series(aov_rev, SSn),
             SSd = use_series(aov_rev, SSd),
             `F` = use_series(aov_rev, `F`))
    full_table <- aov_filter %>%
      mutate(GGe = as.numeric(NA)) %>%
      union(gg_table) %>% 
      mutate(n2p = SSn / (SSn + SSd)) %>%
      relocate(p_adj, .before = GGe) %>%
      relocate(n2p, .after = `F`)
    output <- full_table[grep_vec_unique(x = full_table$Effect,
                                         pattern = perl_regex_without(anova$Effect),
                                         perl = TRUE),]
    return(output)
  }
}

# extract bf by model comparison ~ type III ss
tabulate_bf <- function(bf) {
  require(magrittr)
  nums <- unique(use_series(names(bf), numerator))
  den <- unique(use_series(names(bf), denominator))
  I <- length(nums)
  effnames <- vector("character", I)
  values <- vector("double", I)
  extract_effname <- function(num, den_string) {
    grep(x = unlist(strsplit(den_string, split = " + ", fixed = TRUE)),
         pattern = paste(perl_regex_without(unlist(strsplit(num, split = " + ", fixed = TRUE))),
                         collapse = "|"), 
         perl = TRUE, value = TRUE, invert = TRUE)
  }
  for(i in I:1) {
    values[I-(i-1)] <- bf[grep(x = nums, pattern = nums[i], fixed = TRUE)]
    effnames[I-(i-1)] <- extract_effname(nums[i], den)
  }
  return(data.frame(Effect = effnames, BF01 = values, BF10 = (1/values)))
}

# combine ez and bf tables
merge_ez_bf <- function(ez, bf) {
  require(tidyverse)
  ez <- filter(ez, Effect != "(Intercept)")
  rows <- grep_vec_unique(x = bf$Effect, pattern = reg_eff_matches(ez$Effect),
                          perl = TRUE)
  bf_n <- select(bf, -Effect)
  names_ez <- colnames(ez)
  names_bf <- colnames(bf_n)
  for(i in 1:ncol(bf_n))
    ez[, ncol(ez)+1] <- bf_n[rows, i]
  colnames(ez) <- c(names_ez, names_bf)
  return(ez)
}

# add pastable column to ez and bf output
add_pastable <- function(aov_table) {
  mutate(aov_table,
         paste = if_else(p_adj == "none",
                         if_else(p < .001,
                                 paste("[F", "(", DFn, ",", DFd, ")", "=", round(`F`, 3), ", ",
                                       "p<.001", ", ",
                                       "n2p=", round(n2p, 3), ", BF10=", signif(BF10, 4), "]",
                                       sep = ""),
                                 if_else(p < .05,
                                         paste("[F", "(", DFn, ",", DFd, ")", "=", round(`F`, 3), ", ",
                                               "p=", sub('.', '', as.character(round(p, 3))), ", ",
                                               "n2p=", round(n2p, 3), ", BF10=", signif(BF10, 4), "]",
                                               sep = ""),
                                         paste("[F", "(", DFn, ",", DFd, ")", "=", round(`F`, 3), ", ",
                                               "p=", sub('.', '', as.character(round(p, 3))), ", ",
                                               "n2p=", round(n2p, 3), 
                                               ", BF01=", signif(BF01, 4), "]", sep = ""))),
                         if_else(p < .001,
                                 paste("[F", "(", DFn, ",", DFd, ")", "=", round(`F`, 3), ", ",
                                       "p[GG]<.001", ", ",
                                       "n2p=", round(n2p, 3), ", BF10=", signif(BF10, 4), "]",
                                       sep = ""),
                                 if_else(p < .05,
                                         paste("[F", "(", DFn, ",", DFd, ")", "=", round(`F`, 3), ", ",
                                               "p[GG]=", sub('.', '', as.character(round(p, 3))), ", ",
                                               "n2p=", round(n2p, 3), ", BF10=", signif(BF10, 4), "]",
                                               sep = ""),
                                         paste("[F", "(", DFn, ",", DFd, ")", "=", round(`F`, 3), ", ",
                                               "p[GG]=", sub('.', '', as.character(round(p, 3))), ", ",
                                               "n2p=", round(n2p, 3), ", BF01=", signif(BF01, 4), "]",
                                               sep = "")
                                         )))
  )
}

# prepare output in one step
prep_table <- function(ez, bf) {
  add_pastable(merge_ez_bf(merge_gg(ez), tabulate_bf(bf)))
}

# HSD main effect using ez output and data
HSD_test_main <- function(data, dv, id, aov_table, effect) {
  require(tidyverse)
  lev_across <- data %>% pull({{effect}}) %>% levels(.)
  dv_chr <- deparse(substitute(dv))
  across_chr <- deparse(substitute(effect))
  comb_across <- combn(lev_across, m = 2)
  tib <- data.frame(across_1 = comb_across[1,], across_2 = comb_across[2,]) %>%
    as_tibble()
  a_1 <- as.list(tib$across_1)
  a_2 <- as.list(tib$across_2)
  inputs <- list(a_1, a_2)
  output <- pmap(.l = inputs,
                 .f = function(a_1, a_2) {
                   sub <- data %>% group_by({{id}})
                   values <- summarise(.data = sub,
                                       diff = mean({{dv}}[{{effect}} == a_1]) -
                                         mean({{dv}}[{{effect}} == a_2]))
                   return(values)
                 })
  mdiffs <- map(.x = output,
                .f = function(x) {
                  ungroup(x) %>% pull(diff) %>% mean(.)
                })
  tib$Mean_Difference <- unlist(mdiffs)
  across_1.name <- paste(across_chr, "1", sep = "_")
  across_2.name <- paste(across_chr, "2", sep = "_")
  colnames(tib) <- c(across_1.name, across_2.name, "Mean_Difference")
  n_input <- data %>% pull({{id}})
  HSDmmd <- qtukey(p = .95, nmeans = length(lev_across), 
                   df = filter(aov_table, Effect == across_chr) %>% pull(DFd)) *
    sqrt(((filter(aov_table, Effect == across_chr) %>% pull(SSd))/
            (filter(aov_table, Effect == across_chr) %>% pull(DFd)))/
           n_distinct(n_input))
  HSD_table <- tib %>% mutate(HSDmmd = HSDmmd) %>% rowwise() %>%
    mutate(Sig = if_else(abs(Mean_Difference) > HSDmmd, "*", "ns"))
  return(HSD_table)
}

# HSD interaction using ez output and data
HSD_test_int2 <- function(data, dv, by, across, id, aov_table, effect) {
  require(tidyverse)
  lev_across <- data %>% pull({{across}}) %>% levels(.)
  lev_by <- data %>% pull({{by}}) %>% levels(.)
  by_chr <- deparse(substitute(by))
  dv_chr <- deparse(substitute(dv))
  across_chr <- deparse(substitute(across))
  comb_across <- combn(lev_across, m = 2)
  tib <- data.frame(by_level = rep(lev_by, each = ncol(comb_across)),
                    across_1 = rep(comb_across[1,], times = length(lev_by)),
                    across_2 = rep(comb_across[2,], times = length(lev_by))) %>%
    as_tibble()
  b <- as.list(tib$by_level)
  a_1 <- as.list(tib$across_1)
  a_2 <- as.list(tib$across_2)
  inputs <- list(b, a_1, a_2)
  output <- pmap(.l = inputs,
                 .f = function(b, a_1, a_2) {
                   sub <- filter(data, {{by}} == b) %>%
                     group_by({{id}})
                   values <- summarise(.data = sub,
                                       diff = {{dv}}[{{across}} == a_1] -
                                         {{dv}}[{{across}} == a_2])
                   return(values)
                 })
  mdiffs <- map(.x = output,
                .f = function(x) {
                  ungroup(x) %>%
                    pull(diff) %>%
                    mean(.)
                })
  tib$Mean_Difference <- unlist(mdiffs)
  across_1.name <- paste(across_chr, "1", sep = "_")
  across_2.name <- paste(across_chr, "2", sep = "_")
  colnames(tib) <- c(by_chr, across_1.name, across_2.name, "Mean_Difference")
  n_input <- data %>% pull({{id}})
  effect_chr <- deparse(substitute(effect))
  HSDmmd <- qtukey(p = .95, nmeans = length(lev_across), 
                   df = filter(aov_table, Effect == across_chr) %>% pull(DFd),
                   nranges = n_distinct(lev_by)) *
    sqrt(((filter(aov_table, Effect == effect_chr) %>% pull(SSd))/
            (filter(aov_table, Effect == effect_chr) %>% pull(DFd)))/
           n_distinct(n_input))
  HSD_table <- tib %>% mutate(HSDmmd = HSDmmd) %>% rowwise() %>%
    mutate(Sig = if_else(abs(Mean_Difference) > HSDmmd, "*", "ns"))
  return(HSD_table)
}

# pairwise effect size + bf
pairwise_table <- function(data, dv, by, across, id, stat) {
  require(tidyverse)
  require(BayesFactor)
  lev_across <- data %>%
    pull({{across}}) %>%
    levels(.)
  lev_by <- data %>%
    pull({{by}}) %>%
    levels(.)
  by_chr <- deparse(substitute(by))
  dv_chr <- deparse(substitute(dv))
  across_chr <- deparse(substitute(across))
  comb_across <- combn(lev_across, m = 2)
  tib <- data.frame(by_level = rep(lev_by, each = ncol(comb_across)),
                    across_1 = rep(comb_across[1,], times = length(lev_by)),
                    across_2 = rep(comb_across[2,], times = length(lev_by))) %>%
    as_tibble()
  b <- as.list(tib$by_level)
  a_1 <- as.list(tib$across_1)
  a_2 <- as.list(tib$across_2)
  inputs <- list(b, a_1, a_2)
  output <- pmap(.l = inputs,
                 .f = function(b, a_1, a_2) {
                   sub <- filter(data, {{by}} == b) %>%
                     group_by({{id}})
                   values <- summarise(.data = sub,
                                       diff = {{dv}}[{{across}} == a_1] -
                                         {{dv}}[{{across}} == a_2])
                   return(values)})
  diff_stat <- map(.x = output,
                   .f = function(x) {
                     out <- vector("list")
                     k <- ungroup(x) %>%
                       pull(diff)
                     if(stat == "ttestBF") 
                       ttestBF(x = k, mu = 0) %>% extractBF %>% pull(bf)
                     else if(stat == "hedge's_g") {
                       a <- hedges_g(x = k, mu = 0, ci = .95, alternative = "two.sided")
                       return(data.frame(Hedges_g = a$Hedges_g, Lower = a$CI_low, Upper = a$CI_high))
                     }
                   })
  tib$Difference_Stat <- if(stat == "ttestBF") {
    unlist(diff_stat) 
  } else {
    reduce(.x = diff_stat, .f = rbind)
  }
  across_1.name <- paste(across_chr,
                         "1",
                         sep = "_")
  across_2.name <- paste(across_chr,
                         "2",
                         sep = "_")
  colnames(tib) <- c(by_chr, across_1.name, across_2.name, "Difference_Stat")
  if(stat == "ttestBF") return(tib %>% rename("BF10" = "Difference_Stat"))
  else return(tib %>%
                mutate(Hedges_g = Difference_Stat$Hedges_g,
                       Lower = Difference_Stat$Lower,
                       Upper = Difference_Stat$Upper) %>%
                select(-Difference_Stat))
}

extract_aov_gg <- function(x, intxn = TRUE) {
  require(tidyverse)
  require(maggritr)
  if(isTRUE(intxn)) {
    y <- use_series(x, ANOVA) %>%
      select(-ges)
    q <- filter(y, Effect == "Dose" | Effect == "Group:Dose")
    z <- use_series(x, `Sphericity Corrections`) %>%
      select(Effect:`p[GG]<.05`) %>%
      rename("p" = "p[GG]",
             "p<.05" = "p[GG]<.05") %>%
      mutate(DFn = use_series(q, DFn),
             DFd = use_series(q, DFd),
             SSn = use_series(q, SSn),
             SSd = use_series(q, SSd),
             `F` = use_series(q, `F`))
    a <- y %>%
      filter(Effect != "Dose" & Effect != "Group:Dose") %>%
      mutate(GGe = -999) %>%
      union(z) %>% 
      mutate(n2p = SSn / (SSn + SSd)) %>%
      mutate(p_adj = c("none", "none", "[GG]", "[GG]")) %>%
      relocate(p_adj, .before = GGe) %>%
      relocate(n2p, .after = `F`)
    return(a) 
  } else {
    y <- use_series(x, ANOVA) %>%
      select(-ges)
    q <- filter(y, Effect == "Dose")
    z <- use_series(x, `Sphericity Corrections`) %>%
      select(Effect:`p[GG]<.05`) %>%
      rename("p" = "p[GG]",
             "p<.05" = "p[GG]<.05") %>%
      mutate(DFn = use_series(q, DFn),
             DFd = use_series(q, DFd),
             SSn = use_series(q, SSn),
             SSd = use_series(q, SSd),
             `F` = use_series(q, `F`))
    a <- y %>%
      filter(Effect != "Dose") %>%
      mutate(GGe = -999) %>%
      union(z) %>% 
      mutate(n2p = SSn / (SSn + SSd)) %>%
      mutate(p_adj = c("none", "[GG]")) %>%
      relocate(p_adj, .before = GGe) %>%
      relocate(n2p, .after = `F`)
    return(a) 
  }
}

# Extract partial eta squared from ez table
ez_p_etasq <- function(model, effect) (
  (model$ANOVA["SSn"][model$ANOVA["Effect"] == effect])/(
    (model$ANOVA["SSn"][model$ANOVA["Effect"] == effect]) +
      (model$ANOVA["SSd"][model$ANOVA["Effect"] == effect]))
)

# Extract MSE from ez table
ez_MSE <- function(model, effect) (
  (model$ANOVA["SSd"][model$ANOVA["Effect"] == effect])/(
    (model$ANOVA["DFd"][model$ANOVA["Effect"] == effect])
  ))

# Calculate LSDmmd from ez table
# LSDmmd = tcrit*(sqrt(2*MSE)/sqrt(n))
ez_LSDmmd <- function(model, effect, data, id) {
  require(tidyverse)
  require(magrittr)
  n_input <- data %>%
    pull({{id}})
  effect_chr <- deparse(substitute(effect))
  LSDmmd <- qt(p = 0.05,
               df = (model$ANOVA["DFd"][model$ANOVA["Effect"] == effect_chr]))*(-1)*
    (sqrt(2*(
      (model$ANOVA["SSd"][model$ANOVA["Effect"] == effect_chr])/(
        (model$ANOVA["DFd"][model$ANOVA["Effect"] == effect_chr])
      )))/sqrt(n_distinct(n_input)))
  return(LSDmmd)
}

# Mean differences table
mDiff_table <- function(data, dv, by, across, id) {
  require(tidyverse)
  require(magrittr)
  lev_across <- data %>%
    pull({{across}}) %>%
    levels(.)
  lev_by <- data %>%
    pull({{by}}) %>%
    levels(.)
  by_chr <- deparse(substitute(by))
  dv_chr <- deparse(substitute(dv))
  across_chr <- deparse(substitute(across))
  comb_across <- combn(lev_across, m = 2)
  tib <- data.frame(by_level = rep(lev_by, each = ncol(comb_across)),
                    across_1 = rep(comb_across[1,], times = length(lev_by)),
                    across_2 = rep(comb_across[2,], times = length(lev_by))) %>%
    as_tibble()
  b <- as.list(tib$by_level)
  a_1 <- as.list(tib$across_1)
  a_2 <- as.list(tib$across_2)
  inputs <- list(b, a_1, a_2)
  output <- pmap(.l = inputs,
                 .f = function(b, a_1, a_2) {
                   sub <- filter(data, {{by}} == b) %>%
                     group_by({{id}})
                   values <- summarise(.data = sub,
                                       diff = {{dv}}[{{across}} == a_1] -
                                         {{dv}}[{{across}} == a_2])
                   return(values)})
  mdiffs <- map(.x = output,
                .f = function(x) {
                  ungroup(x) %>%
                    pull(diff) %>%
                    mean(.)
                })
  tib$Mean_Difference <- unlist(mdiffs)
  across_1.name <- paste(across_chr,
                         "1",
                         sep = "_")
  across_2.name <- paste(across_chr,
                         "2",
                         sep = "_")
  colnames(tib) <- c(by_chr, across_1.name, across_2.name, "Mean_Difference")
  return(tib)
}

# LSD test for main effect
LSD_test_main <- function(data, dv, id, model, effect) {
  require(tidyverse)
  require(magrittr)
  lev_across <- data %>%
    pull({{effect}}) %>%
    levels(.)
  
  dv_chr <- deparse(substitute(dv))
  
  across_chr <- deparse(substitute(effect))
  
  comb_across <- combn(lev_across, m = 2)
  
  tib <- data.frame(across_1 = comb_across[1,],
                    across_2 = comb_across[2,]) %>%
    as_tibble()
  
  a_1 <- as.list(tib$across_1)
  
  a_2 <- as.list(tib$across_2)
  
  inputs <- list(a_1, a_2)
  
  output <- pmap(.l = inputs,
                 .f = function(a_1, a_2) {
                   sub <- data %>%
                     group_by({{id}})
                   
                   values <- summarise(.data = sub,
                                       diff = mean({{dv}}[{{effect}} == a_1]) -
                                         mean({{dv}}[{{effect}} == a_2]))
                   
                   return(values)
                 })
  
  mdiffs <- map(.x = output,
                .f = function(x) {
                  ungroup(x) %>%
                    pull(diff) %>%
                    mean(.)
                })
  
  tib$Mean_Difference <- unlist(mdiffs)
  
  across_1.name <- paste(across_chr,
                         "1",
                         sep = "_")
  
  across_2.name <- paste(across_chr,
                         "2",
                         sep = "_")
  
  colnames(tib) <- c(across_1.name, across_2.name, "Mean_Difference")
  
  n_input <- data %>%
    pull({{id}})
  
  LSDmmd <- qt(p = 0.05/2,
               df = (model["DFd"][model["Effect"] == across_chr]),
               lower.tail = FALSE)*
    sqrt(2*model["SSd"][model["Effect"] == across_chr]/
           model["DFd"][model["Effect"] == across_chr])/
    sqrt(n_distinct(n_input))
  
  LSD_table <- tib %>%
    mutate(LSDmmd = LSDmmd) %>%
    rowwise() %>%
    mutate(Sig = if_else(abs(Mean_Difference) > LSDmmd,
                         "*",
                         "ns"))
  
  return(LSD_table)
}

# LSD test for 2-way interaction
LSD_test_int2 <- function(data, dv, by, across, id, model, effect) {
  require(tidyverse)
  require(magrittr)
  lev_across <- data %>%
    pull({{across}}) %>%
    levels(.)
  
  lev_by <- data %>%
    pull({{by}}) %>%
    levels(.)
  
  by_chr <- deparse(substitute(by))
  
  dv_chr <- deparse(substitute(dv))
  
  across_chr <- deparse(substitute(across))
  
  comb_across <- combn(lev_across, m = 2)
  
  tib <- data.frame(by_level = rep(lev_by, each = ncol(comb_across)),
                    across_1 = rep(comb_across[1,], times = length(lev_by)),
                    across_2 = rep(comb_across[2,], times = length(lev_by))) %>%
    as_tibble()
  
  b <- as.list(tib$by_level)
  
  a_1 <- as.list(tib$across_1)
  
  a_2 <- as.list(tib$across_2)
  
  inputs <- list(b, a_1, a_2)
  
  output <- pmap(.l = inputs,
                 .f = function(b, a_1, a_2) {
                   sub <- filter(data, {{by}} == b) %>%
                     group_by({{id}})
                   
                   values <- summarise(.data = sub,
                                       diff = {{dv}}[{{across}} == a_1] -
                                         {{dv}}[{{across}} == a_2])
                   
                   return(values)
                 })
  
  mdiffs <- map(.x = output,
                .f = function(x) {
                  ungroup(x) %>%
                    pull(diff) %>%
                    mean(.)
                })
  
  tib$Mean_Difference <- unlist(mdiffs)
  
  across_1.name <- paste(across_chr,
                         "1",
                         sep = "_")
  
  across_2.name <- paste(across_chr,
                         "2",
                         sep = "_")
  
  colnames(tib) <- c(by_chr, across_1.name, across_2.name, "Mean_Difference")
  
  n_input <- data %>%
    pull({{id}})
  
  effect_chr <- deparse(substitute(effect))
  
  LSDmmd <- qt(p = 0.05/2,
               df = model$ANOVA["DFd"][model$ANOVA["Effect"] == effect_chr],
               lower.tail = FALSE)*
    sqrt(2*model$ANOVA["SSd"][model$ANOVA["Effect"] == effect_chr]/
           model$ANOVA["DFd"][model$ANOVA["Effect"] == effect_chr])/
    sqrt(n_distinct(n_input))
  
  LSD_table <- tib %>%
    mutate(LSDmmd = LSDmmd) %>%
    rowwise() %>%
    mutate(Sig = if_else(abs(Mean_Difference) > LSDmmd,
                         "*",
                         "ns"))
  
  return(LSD_table)
}

# LSD test for 2-way between-subjects interaction
LSD_test_int2_btw <- function(data, dv, by, across, id, model, effect) {
  require(tidyverse)
  require(magrittr)
  lev_across <- data %>%
    pull({{across}}) %>%
    levels(.)
  
  lev_by <- data %>%
    pull({{by}}) %>%
    levels(.)
  
  by_chr <- deparse(substitute(by))
  
  dv_chr <- deparse(substitute(dv))
  
  across_chr <- deparse(substitute(across))
  
  comb_across <- combn(lev_across, m = 2)
  
  tib <- data.frame(by_level = rep(lev_by, each = ncol(comb_across)),
                    across_1 = rep(comb_across[1,], times = length(lev_by)),
                    across_2 = rep(comb_across[2,], times = length(lev_by))) %>%
    as_tibble()
  
  b <- as.list(tib$by_level)
  
  a_1 <- as.list(tib$across_1)
  
  a_2 <- as.list(tib$across_2)
  
  inputs <- list(b, a_1, a_2)
  
  output <- pmap(.l = inputs,
                 .f = function(b, a_1, a_2) {
                   sub1 <- filter(data, {{by}} == b &
                                    {{across}} == a_1) %>%
                     pull({{dv}})
                   
                   sub2 <- filter(data, {{by}} == b & 
                                    {{across}} == a_2) %>%
                     pull({{dv}})
                   
                   values <- mean(sub1) - mean(sub2)
                   
                   return(values)
                 })
  
  diffs <- unlist(output)
  
  tib$Mean_Difference <- diffs
  
  across_1.name <- paste(across_chr,
                         "1",
                         sep = "_")
  
  across_2.name <- paste(across_chr,
                         "2",
                         sep = "_")
  
  colnames(tib) <- c(by_chr, across_1.name, across_2.name, "Difference_Means")
  
  n_summary <- data %>%
    group_by({{across}}) %>%
    summarise(N = n_distinct({{id}}))
  
  n_input <- mean(n_summary$N)
  
  effect_chr <- deparse(substitute(effect))
  
  LSDmmd <- qt(p = 0.05/2,
               df = model$ANOVA["DFd"][model$ANOVA["Effect"] == effect_chr],
               lower.tail = FALSE)*
    sqrt(2*model$ANOVA["SSd"][model$ANOVA["Effect"] == effect_chr]/
           model$ANOVA["DFd"][model$ANOVA["Effect"] == effect_chr])/
    sqrt(n_input)
  
  LSD_table <- tib %>%
    mutate(LSD = LSDmmd) %>%
    rowwise() %>%
    mutate(Sig = if_else(abs(Difference_Means) > LSDmmd,
                         "*",
                         "ns"))
  
  return(LSD_table)
}

# LSD test for 3-way between-subjects interaction
LSD_test_int3_btw <- function(data, dv, by1, by2, across, id, model, effect) {
  require(tidyverse)
  require(magrittr)
  lev_across <- data %>%
    pull({{across}}) %>%
    levels(.)
  
  lev_by1 <- data %>%
    pull({{by1}}) %>%
    levels(.)
  
  lev_by2 <- data %>%
    pull({{by2}}) %>%
    levels(.)
  
  by1_chr <- deparse(substitute(by1))
  
  by2_chr <- deparse(substitute(by2))
  
  dv_chr <- deparse(substitute(dv))
  
  across_chr <- deparse(substitute(across))
  
  comb_across <- combn(lev_across, m = 2)
  
  lev_by_all <- c(lev_by1, lev_by2)
  
  comb_by <- combn(lev_by_all, m = 2) %>%
    as_tibble() %>% 
    select_if(~ any(.x %in% lev_by1) & any(.x %in% lev_by2)) %>%
    as.matrix()
  
  tib <- data.frame(by_1 = rep(comb_by[1,], times = ncol(comb_across)),
                    by_2 = rep(comb_by[2,], times = ncol(comb_across)),
                    across_1 = rep(comb_across[1,], times = ncol(comb_by)),
                    across_2 = rep(comb_across[2,], times = ncol(comb_by))) %>%
    as_tibble()
  
  b_1 <- as.list(tib$by_1)
  
  b_2 <- as.list(tib$by_2)
  
  a_1 <- as.list(tib$across_1)
  
  a_2 <- as.list(tib$across_2)
  
  n_summary <- data %>%
    group_by({{across}}) %>%
    summarise(N = n_distinct({{id}}))
  
  n_input <- mean(n_summary$N)
  
  inputs <- list(b_1, b_2, a_1, a_2)
  
  output <- pmap(.l = inputs,
                 .f = function(b_1, b_2, a_1, a_2) {
                   sub1 <- filter(data, {{by1}} == b_1 & 
                                    {{by2}} == b_2 &
                                    {{across}} == a_1) %>%
                     pull({{dv}})
                   
                   sub2 <- filter(data, {{by1}} == b_1 & 
                                    {{by2}} == b_2 &
                                    {{across}} == a_2) %>%
                     pull({{dv}})
                   
                   values <- mean(sub1) - mean(sub2)
                   
                   return(values)
                 })
  
  diffs <- unlist(output)
  
  tib$Mean_Difference <- diffs
  
  across_1.name <- paste(across_chr,
                         "1",
                         sep = "_")
  
  across_2.name <- paste(across_chr,
                         "2",
                         sep = "_")
  
  colnames(tib) <- c(by1_chr, by2_chr, 
                     across_1.name, across_2.name, 
                     "Difference_Means")
  
  effect_chr <- deparse(substitute(effect))
  
  LSDmmd <- qt(p = 0.05/2,
               df = model$ANOVA["DFd"][model$ANOVA["Effect"] == effect_chr],
               lower.tail = FALSE)*
    sqrt(2*model$ANOVA["SSd"][model$ANOVA["Effect"] == effect_chr]/
           model$ANOVA["DFd"][model$ANOVA["Effect"] == effect_chr])/
    sqrt(n_input)
  
  LSD_table <- tib %>%
    rowwise() %>%
    mutate(LSD = LSDmmd) %>%
    mutate(Sig = if_else(abs(Difference_Means) > LSDmmd,
                         "*",
                         "ns"))
  
  return(LSD_table)
}

conditional_LSD <- function(.data, .dvs, .by, .across, .id, .model, .effect) {
  .dv_string <- paste0(.dvs, collapse = "|")
  efflist <- .model[lapply(X = .model, 
                           FUN = function(x) {
                             x %>%
                               filter(Effect == .effect) %>% 
                               pull(`p<.05`)
                           }) == "*"]
  L <- length(efflist)
  .dv <- vector("character", L)
  out <- vector("list", L)
  b <- function(data, dv, by, across, id, model, effect) {
    lev_across <- levels(data[[across]])
    lev_by <- levels(data[[by]])
    lev_id <- levels(data[[id]])
    comb_across <- combn(lev_across, m = 2)
    tab <- data.frame(by_l = rep(lev_by, each = ncol(comb_across)),
                      ac_l1 = rep(comb_across[1,], times = length(lev_by)),
                      ac_l2 = rep(comb_across[2,], times = length(lev_by)))
    byi <- tab[["by_l"]]
    ac1 <- tab[["ac_l1"]]
    ac2 <- tab[["ac_l2"]]
    I <- length(byi)
    J <- length(lev_id)
    mdiffs <- vector("numeric", I)
    sub <- vector("list", I)
    for(i in 1:I) {
      sub[[i]] <- vector("list", 5)
      names(sub[[i]]) <- c("subj", "by", "ac1", "ac2", "diff")
      for(k in 1:5) {
        sub[[i]][[k]] <- vector("list", J)
      }
      for(j in 1:J) {
        sub[[i]][["subj"]][[j]] <- data[data[id] == lev_id[j],]
        sub[[i]][["by"]][[j]] <- sub[[i]][["subj"]][[j]][sub[[i]][["subj"]][[j]][by] == byi[i],]
        sub[[i]][["ac1"]][[j]] <- sub[[i]][["by"]][[j]][sub[[i]][["by"]][[j]][across] == ac1[i],]
        sub[[i]][["ac2"]][[j]] <- sub[[i]][["by"]][[j]][sub[[i]][["by"]][[j]][across] == ac2[i],]
        sub[[i]][["diff"]][[j]] <- sub[[i]][["ac1"]][[j]][[dv]] - sub[[i]][["ac2"]][[j]][[dv]]
      }
      mdiffs[i] <- mean(unlist(sub[[i]][["diff"]]))
    }
    tab$Mean_Difference <- mdiffs
    across_1.name <- paste(across, "1", sep = "_")
    across_2.name <- paste(across, "2", sep = "_")
    colnames(tab) <- c(by, across_1.name, across_2.name, "Mean_Difference")
    LSDmmd <- qt(p = 0.05/2, df = model["DFd"][model["Effect"] == effect],
                 lower.tail = FALSE)*
      sqrt(2*model["SSd"][model["Effect"] == effect]/
             model["DFd"][model["Effect"] == effect])/
      sqrt(length(unique(data[[id]])))
    tab$LSDmmd <- LSDmmd
    tab$sig <- "ns"
    for(i in 1:I) {
      if(tab$LSDmmd[i] < abs(tab$Mean_Difference[i]))
        tab$sig[i] <- "*"
    }
    return(tab)
  }
  for(l in 1:L) {
    .dv[l] <- str_extract(string = names(efflist)[l], pattern = .dv_string)
    out[[l]] <- b(data = .data, dv = .dv[l], by = .by, across = .across, id = .id, 
                  model = efflist[[l]], effect = .effect)
  }
  names(out) <- names(efflist)
  return(out)
}




