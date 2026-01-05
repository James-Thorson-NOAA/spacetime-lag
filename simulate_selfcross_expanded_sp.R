# This code builds on Sean's simulate_selftest_expanded_sp.R code, 
# To run this code, you need to first run spacetime_lag_EBS_2025-10-28_expanded_sp.R, which fits the basic model for the species we look at here

# Contains 3 key steps:

# 1. do_one_sim(): (adapted from simulate.R). Simulates data from the operating model (obj_oper), loops over all EM configurations
# (config_set), refits each EM to that simulated data, and returns a data.frame of parameter estimates + AIC for all EMs
# for that one simulated dataset.

# 2. run_species_sim(): loads the (AIC-selected) operating model (tmb_dat_oper + obj_oper) for that species, loops n_sim times,
# calling do_one_sim() for each iteration (sequentially). Aggregates results across n_sim replicates for that species.

# 3. future_map2(): loops over species_list in parallel (via furrr) and calls run_species_sim() for each species.

library(TMB)
library(tidyr)
library(dplyr)
library(sf)
library(tictoc)
library(furrr)
library(purrr)
library(readr)
library(stringr)
library(here)

root_dir <- here::here()
version <- "spacetime_lag_2025_10_24"

setwd(root_dir)
compile(paste0(version, ".cpp"), framework = "TMBad")
dyn.load(dynlib(version))

species_list <- c("capelin", "pacific cod", "pacific halibut")

# number of simulation replicates per species. 
# 200 sims takes ~8h with cutoff 60
# 100 sims takes 2.5 days with cutoff 40
n_sim <- 300

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
#base_dir <- file.path(root_dir, "2025-11-09", "capelin")
# for linux, we want to use this with run_species_sim()!
# base_dir <- file.path(root_dir, "2025-11-09", "capelin")
# obj_path_full <- file.path(base_dir, "1-110", "obj.RDS")
# full_obj <- readRDS(obj_path_full)
# par_template <- full_obj$env$parList() 


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
  
  # Load TMB library inside worker (avoiding linux crash)
  if (!version %in% names(getLoadedDLLs())) {
    dyn.load(dynlib(version))
  }
  
  obj_path_full <- file.path(root_dir, "2025-11-09", species, "1-110", "obj.RDS")
  full_obj_local <- readRDS(obj_path_full)
  par_template <- full_obj_local$env$parList()
  
  # build file paths for the species data
  base_dir <- file.path(root_dir, "2025-11-09", species)
  
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
future::plan(future::multisession, workers = 3)
# supposedly more linux friendly, but doesnt work on server
#future::plan(future::multicore, workers = 3)
species_seeds <- map(species_list, ~ sample(seq_len(1e5L), size = n_sim))

tic()
# On Linux server:
# Error: Future (<unnamed-1>) of class MultisessionFuture interrupted, while running on ‘localhost’ (pid 974389)
# all_results <- future_map2(
#   species_list,
#   species_seeds,
#   \(species, seeds) run_species_sim(species, version, n_sim = n_sim, seeds = seeds,
#                                     config_set = config_set)
# )
all_results <- map2(
  species_list,
  species_seeds,
  \(species, seeds) run_species_sim(species, version, n_sim = n_sim, 
                                    seeds = seeds, config_set = config_set)
)
toc()

names(all_results) <- species_list
res_df <- bind_rows(all_results)

# save!
write_csv(res_df, file = file.path(root_dir, paste0("2026-01-05/simulation_selfcross_out_", n_sim, "_reps.csv")))

