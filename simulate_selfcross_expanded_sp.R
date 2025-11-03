# This code builds on Sean's simulate_selftest_expanded_sp.R code, 
# To run this code, you need to first run spacetime_lag_EBS_2025-10-28_expanded_sp.R, which fits the basic model for the species we look at here

# Contains 3 key steps:

# 1. do_one_sim(): (adapted from simulate.R). Simulates data from the operating model (obj_oper), loops over all EM configurations
# (config_set), refits each EM to that simulated data, and returns a data.frame of parameter estimates + AIC for all EMs
# for that one simulated dataset.

# 2. run_species_sim(): loads the (AIC-selected) operating model (tmb_dat_oper + obj_oper) for that species, loops n_sim times,
# calling do_one_sim() for each iteration (sequentially). Aggregates results across n_sim replicates for that species.

# 3. future_map2(): loops over species_list in parallel (via furrr) and calls run_species_sim() for each species.

# FIXME: file paths suggest cutoff = 40, but it was actually 60 when fitted because the path was hardcoded

library(TMB)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggsidekick); theme_set(theme_sleek())
library(sf)
library(tictoc)
library(furrr)
library(purrr)
library(readr)
library(ggbeeswarm)
library(ggh4x)
library(stringr)

root_dir <- here::here()
version <- "spacetime_lag_2025_10_24"

setwd(root_dir)
compile(paste0(version, ".cpp"), framework = "TMBad")
dyn.load(dynlib(version))

species_list <- c("capelin", "pacific cod", "pacific halibut")

# number of simulation replicates per species. 
# 200 sims takes ~8h
n_sim <- 200

# define the model configurations
config_list <- list(
  space = 0:1,
  time = 0:1,
  diffusion = 0,
  quadratic = 1
)
config_set <- expand.grid(config_list)
config_set <- config_set[which(config_set$space == 1 | config_set$diffusion == 0), ]

# load the full parameter template in the global (use capelin but they are all the same)
base_dir <- file.path(root_dir, "2025-10-28_cutoff=40_noprof", "capelin")
obj_path_full <- file.path(base_dir, "1-110", "obj.RDS")
full_obj <- readRDS(obj_path_full)
par_template <- full_obj$env$parList() 

# load Sean's functions
reload_obj <- function(old_obj, tmb_dat) {
  tmb_dat$sim_gmrf <- 1L #< !!
  obj <- MakeADFun(
    data = tmb_dat,
    parameters = old_obj$env$parList(),
    random = c("omega_s", "epsilon_st"),
    map = list(),
    profile = NULL,
    silent = TRUE,
    DLL = version
  )
  obj$env$last.par.best <- old_obj$env$last.par.best
  # Evaluate objective once for MC(), etc. to work
  obj$fn(obj$env$last.par.best)
  obj
}

# grab fixed effects
get_fe <- function(object) {
  lp <- object$env$last.par.best
  p <- numeric(length(lp))
  fe <- object$env$lfixed()
  p[fe] <- lp[fe]
  out <- p[fe]
  names(out) <- names(lp)[fe]
  out
}

# adpted from simulate.R to simulate a new dataset from the AIC selected operating model,
# then refit all four estimation models to the simulated data, and record how close the
# re-estimated parameters are to the original ones, and AIC for all EMs
do_one_sim <- function(obj_oper, tmb_dat_oper, par_template, iter, seed, config_set, version) {
  
  dyn.load(dynlib(version))
  set.seed(seed)
  
  # simulate data from the operating model (the AIC selected model for a given species), obj_oper
  sim <- obj_oper$simulate()
  # add those simulated data (sim$b_i) to tmb_dat_oper
  tmb_dat <- tmb_dat_oper
  tmb_dat$b_i <- sim$b_i
  
  results_list <- list()
  
  # Now that we have a data list with simulated observations (from the AIC model):
  # loop over estimation model parameter configurations (turn off what should be turn off for a given configuration), 
  # and fit the data list with simulated response + config-list-updated par list).
  for (ci in seq_len(nrow(config_set))) {
    
    do_space     = config_set[ci, 'space']
    do_time      = config_set[ci, 'time']
    do_diffusion = config_set[ci, 'diffusion']
    do_quadratic = config_set[ci, 'quadratic']
    type = paste0(do_quadratic, "-", do_space, do_time, do_diffusion)
    
    #par = obj$env$parList()
    # Start each EM from the full parameter template (ensures EM has all possible parameters present initially).
    par <- par_template
    
    # deactivate pars depending on model configuration
    if (do_space != TRUE) {
      par$log_kappaS <- numeric(0)
    }
    if (do_time != TRUE) {
      par$kappaT <- numeric(0)
    }
    if (do_diffusion != TRUE) {
      par$kappaST <- numeric(0)
    }
    if (do_quadratic != TRUE) {
      par$gamma_j <- par$gamma_j[1]
    }
    
    # fit models
    obj_sim <- MakeADFun(
      data = tmb_dat,
      parameters = par,
      random = c("omega_s", "epsilon_st"),
      map = list(),
      profile = NULL,
      silent = TRUE,
      DLL = version
    )
    
    opt_sim <- tryCatch(
      nlminb(
        start = obj_sim$par,
        objective = obj_sim$fn,
        gradient = obj_sim$gr,
        control = list(trace = 0, eval.max = 1e4, iter.max = 1e4)
      ),
      error = function(e) NULL
    )
    if (is.null(opt_sim)) next
    
    opt_sim$AIC <- 2 * length(opt_sim$par) + 2 * opt_sim$obj
    par_true <- get_fe(obj_oper)
    par_hat  <- opt_sim$par
    
    # store results
    results_list[[ci]] <- data.frame(
      par_name = names(par_hat),
      par_true = par_true[names(par_hat)],
      par_hat  = par_hat,
      iter     = iter,
      seed     = seed,
      type     = type,
      AIC      = opt_sim$AIC
    )
  }
  
  bind_rows(results_list)
}

# function to set up everything needed for do_one_sim(), and then does that (sequentially) for n_sim replicates.
run_species_sim <- function(species, version, n_sim, seeds, config_set) {
  
  # need to load here again
  dyn.load(dynlib(version))
  
  # build file paths for the species data
  base_dir <- file.path(root_dir, "2025-10-28_cutoff=40_noprof", species)
  data_path <- file.path(base_dir, "data_sf.RDS")
  
  # load AIC selected model for each species do use in the simulation part of do_one_sim
  # config codes are: type <- paste0(do_quadratic, "-", do_space, do_time, do_diffusion)
  if (species == "capelin"){
    
    tmb_path_oper <- file.path(base_dir, "1-110", "tmbdata.RDS")
    obj_path_oper <- file.path(base_dir, "1-110", "obj.RDS")
    
  } else if (species == "pacific cod"){
    
    tmb_path_oper <- file.path(base_dir, "1-000", "tmbdata.RDS")
    obj_path_oper <- file.path(base_dir, "1-000", "obj.RDS")
    
  } else if (species == "pacific halibut"){
    
    tmb_path_oper <- file.path(base_dir, "1-010", "tmbdata.RDS")
    obj_path_oper <- file.path(base_dir, "1-010", "obj.RDS")
    
  }
  
  # load the operating model given the species-specific paths
  tmb_dat_oper <- readRDS(tmb_path_oper)
  old_obj_oper <- readRDS(obj_path_oper)
  
  # make sure the operating model is properly reloaded
  obj_oper <- reload_obj(old_obj_oper, tmb_dat_oper)
  
  # per-species simulation
  out <- map_dfr(
    seq_len(n_sim),
    \(i) do_one_sim(
      obj_oper,
      tmb_dat_oper,
      par_template,
      iter = i,
      seed = seeds[i],
      config_set = config_set,
      version = version
    )
  )
  
  out$species <- species
  out
}

# now we do the top-level future_map2, which does run_species_sim() for each species in parallel
future::plan(future::multisession, workers = length(species_list))
species_seeds <- map(species_list, ~ sample(seq_len(1e5L), size = n_sim))

tic()
all_results <- future_map2(
  species_list,
  species_seeds,
  \(species, seeds) run_species_sim(species, version, n_sim = n_sim, seeds = seeds,
                                    config_set = config_set)
)
toc()

names(all_results) <- species_list
res_df <- bind_rows(all_results)

# save!
write_csv(res_df, file = file.path(root_dir, paste0("2025-10-28_cutoff=40_noprof/simulation_selfcross_out_", n_sim, "_reps.csv")))

# plot
res_df <- read_csv(file.path(root_dir, paste0("2025-10-28_cutoff=40_noprof/simulation_selfcross_out_", n_sim, "_reps.csv")))

# specify true model for each species
correct_models <- tibble(
  species = c("capelin", "pacific cod", "pacific halibut"),
  correct_type = c("space+time", "base", "time")
)

res_df |> 
  distinct(AIC, type, iter, species) |> 
  mutate(type = case_when(
    type == "1-000" ~ "base",
    type == "1-100" ~ "space",
    type == "1-010" ~ "time",
    type == "1-110" ~ "space+time",
    TRUE ~ type
  )) |>
  mutate(delta_AIC = AIC - min(AIC), .by = c(iter, species)) |> 
  summarise(fraction_best = mean(delta_AIC == 0), .by = c(type, species)) |>
  left_join(correct_models, by = "species") |>
  mutate(true_model = type == correct_type) |>
  ggplot(aes(fraction_best, type, fill = true_model)) +
  facet_wrap(~str_to_sentence(species), ncol = 3) +
  theme(legend.position = "bottom") +
  geom_col(width = 0.9) +
  labs(
    y = "Model",
    x = "Fraction of iterations with best AIC"
  ) +
  scale_fill_brewer(palette = "Set2", name = "Operating model",
                    direction = -1)

ggsave(paste0(root_dir, "/figs/self_test_aic.png"), width = 20, height = 8, unit = "cm")


# Now compare how well the full and the AIC model can recover kappaS and kappaT (for halibut and capelin)
res_df |> 
  filter(par_name %in% c("log_kappaS", "kappaT")) |> 
  filter(species %in% c("capelin", "pacific halibut")) |> 
  drop_na(par_true) |> 
  mutate(type = case_when(
    type == "1-000" ~ "base",
    type == "1-100" ~ "space",
    type == "1-010" ~ "time",
    type == "1-110" ~ "space+time",
    TRUE ~ type
  )) |>
  mutate(par_label_parsed = case_when(
    par_name == "kappaT" ~ "kappa[T]",
    par_name == "log_kappaS" ~ "log(kappa[S])",
    TRUE ~ par_name
  ),
  species_label = str_to_sentence(species)) |>
  left_join(correct_models, by = "species") |>
  mutate(true_model = type == correct_type) |>
  ggplot(aes(type, par_hat, color = true_model)) +
  geom_hline(aes(yintercept = par_true), linetype = 2, color = "tomato") +
  geom_quasirandom(alpha = 0.175, size = 1) +
  geom_boxplot(fill = NA, width = 0.1, size = 0.4, alpha = 1,
               outlier.shape = NA, color = "gray30") +
  geom_text(data = . %>% distinct(species_label, par_label_parsed),
            aes(label = par_label_parsed, x = -Inf, y = Inf),
            hjust = -0.1, vjust = 1.2, size = 4, color = "gray30", 
            parse = TRUE) +
  geom_text(data = . %>% distinct(species_label, par_label_parsed),
            aes(label = species_label, x = Inf, y = Inf),
            hjust = 1.1, vjust = 1.2, size = 4, color = "gray30") +
  facet_wrap(~ species_label+par_label_parsed, 
             scales = "free", 
             ncol = 3) +
  labs(y = "Estimated value",
       x = "Estimation model") +
  scale_color_brewer(palette = "Set2", name = "AIC-selected model",
                     direction = -1) + 
  theme(legend.position = "bottom",
        strip.text.x = element_blank()) +
  guides(color = guide_legend(override.aes = list(alpha = 0.9)))

ggsave(paste0(root_dir, "/figs/kappa_recovery.png"), width = 20, height = 10, unit = "cm")