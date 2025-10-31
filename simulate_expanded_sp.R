# This code builds on Sean's simulate.R code, with the function `do_one_sim()`
# Here we add run_species_sim(), which runs do_one_sim() n_sim times sequentially for a given species
# The final code, `future_map2(species_list, species_seeds, ...)`, runs all species in parallel, calling run_species_sim() for each

# To run this code, you need to first run spacetime_lag_EBS_2025-10-28_expanded_sp.R, which fits the basic model for the species we look at here

library(TMB)
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

root_dir <- here::here()
version <- "spacetime_lag_2025_10_24"

setwd(root_dir)
compile(paste0(version, ".cpp"), framework = "TMBad")
dyn.load(dynlib(version))

species_list <- c("capelin", "pacific cod", "pacific halibut")

# number of simulation replicates per species. 
# 250 sims take 3.4 hours
n_sim <- 250

# load Sean's functions
reload_obj <- function(old_obj, tmb_dat) { # ML: add tmb_dat here because we don't load it upfront since we have more species now
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

# function to simulate a new dataset from the fitted TMB model, re-fit the model to that simulated data,
# record how close the re-estimated parameters are to the original ones.
do_one_sim <- function(obj, tmb_dat, iter, seed, is_parallel = FALSE) {
  
  dyn.load(dynlib(version))
  set.seed(seed)

  sim <- obj$simulate()
  
  # fit back to the simulated data:
  tmb_dat$b_i <- sim$b_i
  pl <- obj$env$parList()
  
  obj_sim <- MakeADFun(
    data = tmb_dat,
    parameters = pl,
    random = c("omega_s", "epsilon_st"),
    map = list(),
    profile = NULL,
    silent = TRUE
  )
  
  opt_sim <- tryCatch(nlminb(
    start = obj_sim$par,
    obj = obj_sim$fn,
    gr = obj_sim$gr,
    control = list(trace = 1L, eval.max = 1e4L, iter.max = 1e4L)
  ))
  
  par_true <- get_fe(obj)
  par_hat <- opt_sim$par
  
  data.frame(
    par_name = names(par_hat),
    par_true = par_true,
    par_hat = par_hat,
    iter = iter,
    seed = seed
  )
}

# by species, loads the species data and runs do_one_sim() sequentially
run_species_sim <- function(species, version, n_sim, seeds) {
  # need to load here again
  dyn.load(dynlib(version))
  
  # build file paths for the species data
  base_dir <- file.path(root_dir, "2025-10-28_cutoff=40_noprof", species)
  data_path <- file.path(base_dir, "data_sf.RDS")
  tmb_path  <- file.path(base_dir, "1-111", "tmbdata.RDS")
  obj_path  <- file.path(base_dir, "1-111", "obj.RDS")
  
  # load species objects
  dat <- readRDS(data_path)
  tmb_dat <- readRDS(tmb_path)
  old_obj <- readRDS(obj_path)
  obj <- reload_obj(old_obj, tmb_dat) # note tmb_dat now argument in reload
  
  # per-species parallel simulation
  out <- map_dfr(
    seq_len(n_sim),
    \(i) do_one_sim(obj, tmb_dat, seed = seeds[i], iter = i)
  )
  
  out$species <- species
  out
}

# now we do the top-level future_map2, which does run_species_sim() for each species in parallel
future::plan(future::multisession, workers = length(species_list))

# create a unique set of seeds per species
species_seeds <- map(species_list, ~ sample(seq_len(1e5L), size = n_sim))

tic()
all_results <- future_map2(
  species_list,
  species_seeds,
  \(species, seeds) run_species_sim(species, version, n_sim = n_sim, seeds = seeds)
)
toc()

names(all_results) <- species_list
res_df <- bind_rows(all_results)

# save!
write_csv(res_df, file = file.path(root_dir, paste0("2025-10-28_cutoff=40_noprof/simulation_out_", n_sim, "_reps.csv")))



# plot
res_df <- read_csv(file.path(root_dir, paste0("2025-10-28_cutoff=40_noprof/simulation_out_", n_sim, "_reps.csv")))

out_sub <- res_df |> 
  filter(par_name %in% c("gamma_j", "kappaT", "kappaST", "log_kappaS", "ln_tauE", "ln_tauO")) |> 
  # the cumsum is odd here, but if I use e.g. rownumber() it becomes gamma_j4 and 5, this works because it only looks at gamma_j not all pars ... 
  mutate(
    par_name = if_else(
      par_name == "gamma_j",
      paste0("gamma_j", cumsum(par_name == "gamma_j")),
      par_name
    ),
    .by = c(species, iter)
  )

out_sub |> 
  #filter(!par_name %in% c("ln_tauE", "ln_tauO")) |>
  # tidylog::filter(!(par_name == "gamma_j1" & par_hat < -4)) |>
  # tidylog::filter(!(par_name == "gamma_j2" & par_hat < -2)) |>
  # tidylog::filter(!(par_name == "gamma_j2" & par_hat > 2)) |>
  ggplot(aes(x = 1, par_hat)) +
  geom_hline(aes(yintercept = par_true), linetype = 2, color = "tomato") +
  geom_quasirandom(alpha = 0.1, color = "grey10", size = 0.5) +
  geom_boxplot(fill = NA, width = 0.1, size = 0.4, alpha = 1,
               outlier.shape = NA, color = "gray30") +
  #geom_violin(fill = NA, width = 0.8, color = alpha("black", 0.1)) +
  ggh4x::facet_grid2(par_name~species, scales = "free", independent = "y") +
  labs(y = "Estimated value") +
  guides(color = "none") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

ggsave(paste0(root_dir, "/figs/self_test.png"), width = 19, height = 23, unit = "cm")

out_sub |> 
  filter(par_name %in% c("kappaST", "log_kappaS")) |> 
  filter(!species == "pacific halibut") |> 
  tidyr::pivot_wider(values_from = c(par_hat, par_true), names_from = par_name) |> 
  #tidylog::filter(par_hat_kappaST > -100 & par_hat_kappaST < 750) |> 
  mutate(diff = par_hat_log_kappaS - par_true_log_kappaS) |>
  #tidylog::filter(diff > -20 & diff < 20) |> 
  ggplot(aes(x = 1, diff, fill = par_hat_kappaST)) +
  geom_hline(aes(yintercept = 0), linetype = 2, color = "tomato") +
  geom_quasirandom(size = 1, shape = 21, stroke = 0.1, color = "grey30") +
  geom_boxplot(fill = NA, width = 0.1, size = 0.4, alpha = 1,
               outlier.shape = NA, color = "gray30") +
  facet_wrap(~species, scales = "free") +
  scale_fill_gradient2(midpoint = 0) +
  labs(y = "diff (est - true)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

ggsave(paste0(root_dir, "/figs/self_test_kappaST.png"), width = 20, height = 9, unit = "cm")

out_sub |> 
  filter(par_name %in% c("kappaST", "log_kappaS")) |> 
  filter(!species == "pacific halibut") |> 
  tidyr::pivot_wider(values_from = c(par_hat, par_true), names_from = par_name) |> 
  #tidylog::filter(species != "capelin" | between(par_hat_kappaST, -10, 10)) |> 
  mutate(diff = par_hat_kappaST - par_true_kappaST) |>
  #tidylog::filter(diff > -20 & diff < 20) |> 
  ggplot(aes(x = 1, diff, fill = par_hat_log_kappaS)) +
  geom_hline(aes(yintercept = 0), linetype = 2, color = "tomato") +
  geom_quasirandom(size = 1, shape = 21, stroke = 0.05, color = "grey30") +
  geom_boxplot(fill = NA, width = 0.1, size = 0.4, alpha = 1,
               outlier.shape = NA, color = "gray30") +
  facet_wrap(~species, scales = "free") +
  scale_fill_gradient2(midpoint = 0) +
  labs(y = "diff (est - true)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

ggsave(paste0(root_dir, "/figs/self_test_logkappaS.png"), width = 20, height = 9, unit = "cm")

# Correlation between parameters?
out_sub |> 
  filter(par_name %in% c("kappaST", "log_kappaS")) |> 
  filter(!species == "pacific halibut") |> 
  tidyr::pivot_wider(values_from = c(par_hat, par_true), names_from = par_name) |> 
  tidylog::filter(species != "capelin" | between(par_hat_kappaST, -10, 10)) |> 
  ggplot(aes(par_hat_log_kappaS, par_hat_kappaST)) +
  geom_point() +
  facet_wrap(~species, scales = "free")

ggsave(paste0(root_dir, "/figs/self_test_logkappaS_kappaST_corr.png"), width = 20, height = 9, unit = "cm")
