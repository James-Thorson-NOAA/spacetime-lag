library(TMB)
library(dplyr)
library(ggplot2)
library(sf)

root_dir <- here::here()
version <- "spacetime_lag_2025_10_24"

setwd(root_dir)
compile(paste0(version, ".cpp"), framework = "TMBad")
dyn.load(dynlib(version))

dat <- readRDS("2025-10-28_cutoff=40_noprof/capelin/data_sf.RDS")
tmb_dat <- readRDS("2025-10-28_cutoff=40_noprof/capelin/1-111/tmbdata.RDS")
old_obj <- readRDS("2025-10-28_cutoff=40_noprof/capelin/1-111/obj.RDS")

# Recreate the TMB object with the loaded DLL to avoid crashing sometimes on MC()
pl <- old_obj$env$parList()
random <- c("omega_s", "epsilon_st")

reload_obj <- function(old_obj) {
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
  # Evaluate objective once for MC() to work
  obj$fn(obj$env$last.par.best)
  obj
}

obj <- reload_obj(old_obj)

one_sample_posterior <- function(object) {
  # take a sample from the random effects (as MVN) with fixed effects held at MLEs
  tmp <- object$env$MC(n = 1L, keep = TRUE, antithetic = FALSE)
  re_samp <- as.vector(attr(tmp, "samples"))
  lp <- object$env$last.par.best
  p <- numeric(length(lp))
  fe <- object$env$lfixed()
  re <- object$env$lrandom()
  p[re] <- re_samp
  p[fe] <- lp[fe]
  p
}

get_fe <- function(object) {
  lp <- object$env$last.par.best
  p <- numeric(length(lp))
  fe <- object$env$lfixed()
  p[fe] <- lp[fe]
  out <- p[fe]
  names(out) <- names(lp)[fe]
  out
}

set.seed(123)

do_one_sim <- function(obj, tmb_dat, iter, seed, is_parallel = FALSE) {

  # if running in parallel, will need to run this on each iteration
  if (is_parallel) {
    dyn.load(dynlib(version))
    obj <- reload_obj(obj)
  }
  # otherwise, can skip for speed

  set.seed(seed)
  p <- one_sample_posterior(obj)
  sim <- obj$simulate(par = p)
  # sim <- obj$simulate() # EB random effects

  # now fit back to the simulated data:
  tmb_dat$b_i <- sim$b_i
  pl <- obj$env$parList()
  random <- c("omega_s", "epsilon_st")
  map <- list()
  profile <- NULL

  obj_sim <- MakeADFun(
    data = tmb_dat,
    parameters = pl,
    random = random,
    map = map,
    profile = profile,
    silent = TRUE
  )

  opt_sim <- tryCatch(nlminb(
    start = obj_sim$par,
    obj = obj_sim$fn,
    gr = obj_sim$gr,
    control = list(trace = 1L, eval.max = 1e4L, iter.max = 1e4L)
  ))

  # do we want CI coverage? If so, we'll need to run the sdreport(), much slower
  # sd_sim <- sdreport(
  #   obj_sim,
  #   getReportCovariance = TRUE, # must be TRUE to get Cov for `gamma_j` if its included in profile
  #   skip.delta.method = FALSE, # must be FALSE to get as.list( opt$SD, what = "Std. Error", report = TRUE )
  #   ignore.parm.uncertainty = FALSE # must be FALSE to include fixed effect uncertainty
  # )

  par_true <- get_fe(obj)
  par_hat <- opt_sim$par

  par_df <- data.frame(
    par_name = names(par_hat),
    par_true = par_true,
    par_hat = par_hat,
    iter = iter,
    seed = seed
  )
}

# out <- do_one_sim(obj, tmb_dat, seed = 42, iter = 1)

set.seed(123)
seeds <- sample(seq_len(1e5L), size = 500L)

# out <- purrr::map_dfr(
#   seq_len(2), \(i)
#   do_one_sim(obj, tmb_dat, seed = seeds[i], iter = i)
# )

# swap out for furrr once happy?
# OK, this gets complicated; gets crashy
# You have to MakeADFun each time on the obj first
# see my notes in the function above
future::plan(future::multisession)
out <- furrr::future_map_dfr(
  seq_len(10L), \(i)
  do_one_sim(obj, tmb_dat, seed = seeds[i], iter = i, is_parallel = TRUE),
  .options = furrr::furrr_options(seed = TRUE)
)

ggplot(out, aes(par_hat, iter)) +
  geom_point() +
  geom_vline(aes(xintercept = par_true), lty = 2) +
  facet_wrap(vars(par_name), scales = "free_x")
