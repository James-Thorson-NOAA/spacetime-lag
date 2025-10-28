library(TMB)
library(dplyr)
library(ggplot2)
library(sf)

root_dir <- here::here()
version <- "spacetime_lag_2025_10_24"

setwd(root_dir)
compile(paste0(version, ".cpp"), framework = "TMBad")
dyn.load(dynlib(version))

obj <- readRDS("2025-10-28_cutoff=40_noprof/capelin/1-111/obj.RDS")
dat <- readRDS("2025-10-28_cutoff=40_noprof/capelin/data_sf.RDS")
tmb_dat <- readRDS("2025-10-28_cutoff=40_noprof/capelin/1-111/tmbdata.RDS")

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

# sdr <- sdreport(
#   obj
#   # getReportCovariance = TRUE, # must be TRUE to get Cov for `gamma_j` if its included in profile
#   # skip.delta.method = FALSE, # must be FALSE to get as.list( opt$SD, what = "Std. Error", report = TRUE )
#   # ignore.parm.uncertainty = FALSE # must be FALSE to include fixed effect uncertainty
# )

# obj$env$MC() is crashy sometimes!?
# new_par <- one_sample_posterior(obj)
# sim <- obj$simulate(par = new_par)

# mean(sim$b_i)
# mean(dat$catch_weight)
# max(sim$b_i)
# max(dat$catch_weight)

# if crashing, for now use EB random effects, i.e. obj$simulate()

do_one_sim <- function(obj, tmb_dat, iter, seed) {

  set.seed(seed)
  sim <- obj$simulate()

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


  # do we want CI coverage? If so, we'll need to run the sdreport()
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

out <- purrr::map_df(seq_len(2), \(i) do_one_sim(obj, tmb_dat, seed = seeds[i], iter = i))

ggplot(out, aes())
