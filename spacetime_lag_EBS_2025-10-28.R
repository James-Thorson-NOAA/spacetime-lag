root_dir <- here::here()

library(TMB)
library(surveyjoin)
library(sf)
library(fmesher)
library(Matrix)
library(rnaturalearth)
options("sf_max.plot" = 45)

#
# Date = Sys.Date()
Date <- "2025-10-28_cutoff=40_noprof"
date_dir <- file.path(root_dir, Date)
dir.create(date_dir, showWarnings = FALSE)

# species
# Not available: "yellowfin sole"
# Not enough data: "sablefish"
# On boundary: "eulachon"
print(get_species(), n = 100)
species_set <- c("pacific cod", "walleye pollock", "capelin", "pacific herring", "arrowtooth flounder", "pacific halibut")

# original:
# species_set <- c("pacific cod", "walleye pollock", "capelin", "pacific herring", "arrowtooth flounder", "pacific halibut")
# selected:
species_set <- "capelin"

#
# if( !("settings.RDS" %in% list.files(date_dir)) ){

# original:
# config_list = list(
#   space = 0:1,
#   time = 0:1,
#   diffusion = 0:1,
#   quadratic = 1
# )
#
# selected:
config_list <- list(
  space = 1,
  time = 1,
  diffusion = 1,
  quadratic = 1
)

config_set <- expand.grid(config_list)
config_set <- config_set[which(config_set$space == 1 | config_set$diffusion == 0), ]

#
rescale <- 1
cutoff1 <- 60
cutoff2 <- 60
profile <- NULL # c("gamma_j")  # , # "beta_k","log_kappaS","rho_t","rho_st")
covariate <- c("temp", "anomaly")[2]
version <- "spacetime_lag_2025_10_24"
hesscompress <- FALSE # TOO SLOW
# MSD = 4 / kappaS^2
log_kappaS_start <- log(sqrt(4 / (20 / rescale))) # MSD = 20;  log(sqrt(4 / 20))
run_sdreport <- TRUE
L1version <- c("innovations", "stationary")[2] # stationary:  Lag matrix has columns that sum to 0

# settings
settings <- list(
  rescale = rescale, cutoff1 = cutoff1, cutoff2 = cutoff2, profile = profile,
  covariate = covariate, log_kappaS_start = log_kappaS_start,
  config_set = config_set, L1version = L1version, version = version,
  hesscompress = hesscompress
)
saveRDS(settings, file = file.path(date_dir, "settings.RDS"))
capture.output(settings, file = file.path(date_dir, "settings.txt"))
# }else{
#   settings = readRDS( file = file.path(date_dir,"settings.RDS") )
#   attach(settings)
# }

# survey domain from VAST
survey_domain <- st_read("data-raw/EBSshelf/")
survey_domain <- st_transform(survey_domain, crs = "+proj=utm +zone=2 +datum=WGS84 +units=km") # st_crs(data_sf))

#
worldmap <- ne_countries(country = c("united states of america", "russia"), scale = 50)
worldmap <- st_transform(worldmap, crs = "+proj=utm +zone=2 +datum=WGS84 +units=km") # st_crs(data_sf))

# Create extrapolation grid
cellsize <- 25
grid <- st_make_grid(survey_domain, cellsize = cellsize)
grid <- st_intersection(grid, survey_domain)

# Add extra year at beginning
# year_set <- min(data_sf$year):max(data_sf$year)
# year_set = c( year_set[1] - 1, year_set )
year_set <- 1981:2024

# Compile
setwd(root_dir)
compile(paste0(version, ".cpp"), framework = "TMBad")
dyn.load(dynlib(version))

#
if (hesscompress == TRUE) {
  config(tmbad.sparse_hessian_compress = TRUE, DLL = version)
}

# Loop
for (si in seq_along(species_set)) {
  for (ci in seq_len(nrow(config_set))) {
    # for( ci in 1:4 ){
    # for( ci in 5:6 ){

    #
    species <- species_set[si]
    do_space <- config_set[ci, "space"]
    do_time <- config_set[ci, "time"]
    do_diffusion <- config_set[ci, "diffusion"] #
    do_quadratic <- config_set[ci, "quadratic"]
    type <- paste0(do_quadratic, "-", do_space, do_time, do_diffusion)

    # Error checks
    if ((do_space == 0) & (do_diffusion == 1)) stop("Combination does not make sense")

    #
    species_dir <- file.path(date_dir, species)
    run_dir <- file.path(species_dir, type)
    dir.create(run_dir, recursive = TRUE)

    if (!("opt.RDS" %in% list.files(run_dir))) {
      message("running ", species, " ", type)
      # Load data
      d <- get_data(species)
      d <- as.data.frame(d)
      data <- subset(d, (survey_name == "eastern Bering Sea") & !is.na(bottom_temp_c))
      # data = subset( data, year >= 2022 )

      #
      data_sf <- st_as_sf(data, coords = c("lon_start", "lat_start"), crs = "+proj=longlat +datum=WGS84")
      data_sf <- st_transform(data_sf, crs = st_crs("+proj=utm +zone=2 +datum=WGS84 +units=km"))
      saveRDS(data_sf, file = file.path(species_dir, "data_sf.RDS"))

      # Make triangulated mesh
      mesh <- fm_mesh_2d(st_coordinates(data_sf) / rescale, cutoff = cutoff1 / rescale, refine = TRUE)

      #
      mesh2 <- fm_mesh_2d(st_coordinates(data_sf) / rescale, cutoff = cutoff2 / rescale, refine = TRUE)
      saveRDS(mesh2, file = file.path(run_dir, "mesh2.RDS"))

      # Create matrices in INLA
      spde <- fm_fem(mesh, order = 2)
      spde2 <- fm_fem(mesh2, order = 2)
      # create projection matrix from vertices to samples
      A_is <- fm_evaluator(mesh, loc = st_coordinates(data_sf) / rescale)$proj$A
      A_is2 <- fm_evaluator(mesh2, loc = st_coordinates(data_sf) / rescale)$proj$A
      # create projection matrix from vertices to grid
      A_gs <- fm_evaluator(mesh, loc = st_coordinates(st_centroid(grid)) / rescale)$proj$A
      A_gs2 <- fm_evaluator(mesh2, loc = st_coordinates(st_centroid(grid)) / rescale)$proj$A

      #
      Temp_s2t <- array(0, dim = c(mesh2$n, length(year_set)), dimnames = list(NULL, year = year_set))
      for (t in seq_along(year_set)) {
        # At_is = A_is
        # At_is[ which(data_sf$year != year_set[t]), ] = 0
        # Temp_st[,t] = ( t(At_is) %*% data$bottom_temp_c )[,1] / colSums(At_is)

        #
        if (year_set[t] %in% data_sf$year) {
          data_t <- subset(data_sf, year == year_set[t])
          idx <- RANN::nn2(data = st_coordinates(data_t) / rescale, query = mesh2$loc[, 1:2], k = 1)$nn.idx
          Temp_s2t[, t] <- data_t$bottom_temp_c[idx]
        }
      }
      if (1981 %in% year_set) Temp_s2t[, "1981"] <- rowMeans(Temp_s2t[, c("1982", "1983")])
      if (2020 %in% year_set) Temp_s2t[, "2020"] <- rowMeans(Temp_s2t[, c("2019", "2021")])

      # Convert to anomaly
      if (covariate == "anomaly") {
        Temp_s2t <- sweep(Temp_s2t, MARGIN = 1, STAT = rowMeans(Temp_s2t), FUN = "-")
      }

      #
      formula <- ~1
      X_ik <- model.matrix(formula, data_sf)

      #
      n_t <- length(year_set)
      M_tt <- sparseMatrix(i = 2:n_t, j = 1:(n_t - 1), x = 1, dims = c(n_t, n_t)) # bandSparse( n = length(year_set), k = -1, repr = "T")

      # Make inputs
      data <- list(
        options_z = c(switch(L1version,
          "innovations" = 0,
          "stationary" = 1
        ), do_quadratic),
        n_t = n_t,
        a_g = as.numeric(st_area(grid)),
        b_i = data_sf$catch_weight,
        a_i = data_sf$effort / 100, # Convert hectares to km^2
        t_i = match(data_sf$year, year_set) - 1,
        Temp_s2t = Temp_s2t,
        X_ik = X_ik,
        A_is = A_is,
        A_is2 = A_is2,
        A_gs = A_gs,
        A_gs2 = A_gs2,
        M0_ss = spde$c0,
        M1_ss = spde$g1,
        M2_ss = spde$g2,
        invM0_ss = solve(spde$c0),
        M1_s2s2 = spde2$g1,
        invM0_s2s2 = solve(spde2$c0),
        M_tt = M_tt
      )
      par <- list(
        log_kappaS = log_kappaS_start,
        kappaT = 0,
        kappaST = 0, #
        gamma_j = c(0.1, -0.1),
        beta_k = rep(0, ncol(X_ik)),
        ln_tauO = log(10) - log(rescale),
        ln_tauE = log(10) - log(rescale),
        ln_kappa = log(sqrt(8) / (100 / rescale)), # Distance = sqrt(8)/kappa  AND start Distance = 100
        ln_phi = log(1),
        logit_rhoE = 0,
        finv_power = 0,
        omega_s = rep(0, mesh$n),
        epsilon_st = matrix(0, nrow = mesh$n, ncol = data$n_t)
      )

      #
      if (do_space != TRUE) {
        par$log_kappaS <- vector()
      }
      if (do_time != TRUE) {
        par$kappaT <- vector()
      }
      if (do_diffusion != TRUE) {
        par$kappaST <- vector()
      }
      if (do_quadratic != TRUE) {
        par$gamma_j <- par$gamma_j[1]
      }

      random <- c("omega_s", "epsilon_st")
      map <- list()

      # Build and run
      start_time <- Sys.time()
      obj <- MakeADFun(
        data = data,
        parameters = par,
        random = random,
        map = map,
        profile = profile,
        silent = TRUE
      )
      compile_time <- Sys.time() - start_time
      obj$fn(obj$par)
      cbind(obj$par, obj$gr(obj$par))

      # Run
      start_time <- Sys.time()
      opt <- tryCatch(nlminb(
        start = obj$par,
        obj = obj$fn,
        gr = obj$gr,
        control = list(trace = 1, eval.max = 1e4, iter.max = 1e4)
      ))
      if (class(opt) == "try-error") {
        saveRDS(opt, file = file.path(run_dir, "opt.RDS"))
      } else {
        opt$AIC <- 2 * length(opt$par) + 2 * opt$obj
        if (run_sdreport) {
          opt$SD <- sdreport(
            obj,
            getReportCovariance = TRUE, # must be TRUE to get Cov for `gamma_j` if its included in profile
            skip.delta.method = FALSE, # must be FALSE to get as.list( opt$SD, what = "Std. Error", report = TRUE )
            ignore.parm.uncertainty = FALSE # must be FALSE to include fixed effect uncertainty
          )
        }
        opt$compile_time <- compile_time
        opt$runtime <- Sys.time() - start_time
        rep <- obj$report()
        parhat <- obj$env$parList()

        #
        capture.output(opt, file = file.path(run_dir, "opt.txt"))
        saveRDS(opt, file = file.path(run_dir, "opt.RDS"))
        saveRDS(rep, file = file.path(run_dir, "rep.RDS"))
        saveRDS(obj, file = file.path(run_dir, "obj.RDS"))
        saveRDS(parhat, file = file.path(run_dir, "parhat.RDS"))
      }
    }
  }
}

if (FALSE) {
  ##########################
  # Compile results
  ##########################

  #
  # root_dir = R'(C:\Users\james\OneDrive\Desktop\Work files (backup)\Collab-2025\2025 -- spacetime distributed lag)'
  source(file.path(R'(C:\Users\james\OneDrive\Desktop\Work files (backup)\Collab-2025\2025 -- spacetime distributed lag)', "add_legend.R"))

  library(TMB)
  library(Matrix)
  library(sf)
  library(ggplot2)
  library(patchwork)

  # Date = "2025-10-21B"
  # date_dir = file.path(root_dir, Date)

  # Create extrapolation grid
  cellsize <- 10
  grid <- st_make_grid(survey_domain, cellsize = cellsize)
  grid <- st_intersection(grid, survey_domain)

  #
  # species_set = c( "pacific cod", "walleye pollock", "capelin", "pacific herring", "arrowtooth flounder", "pacific halibut" )
  species_names <- sapply(species_set, FUN = gsub, fixed = TRUE, pattern = " ", replacement = "\n")

  #
  if (covariate == "temp") {
    Temp_z <- seq(-2.1, 11.7, by = 0.1)
  } else {
    Temp_z <- seq(-6.1, 6.1, by = 0.1)
  }
  # config_list = list(
  #  space = 0:1,
  #  diffusion = 0,
  #  time = 0:1,
  #  quadratic = 1
  # )
  # config_set = expand.grid( config_list )
  # config_set = config_set[ which(config_set$space==1 | config_set$diffusion==0), ]
  model_names <- c("null", "S", "T", "ST")

  pars <- c("runtime", "pars", "obj", "AIC", "log_kappaS", "se_log_kappaS", "kappaT", "kappaST", "rhoT", "se_rhoT", "AIC_weight", "MSD", "RMSD", "se_RMSD") # "RMSD_low", "RMSD_high"
  results_scz <- array(NA,
    dim = c(length(species_set), nrow(config_set), length(pars)),
    dimnames = list(species = species_set, model = model_names, par = pars)
  )
  gamma_scjz <- array(NA,
    dim = c(length(species_set), nrow(config_set), 2, 2),
    dimnames = list(species_set, model_names, NULL, c("hat", "se"))
  )
  p_sczq <- array(NA,
    dim = c(length(species_set), nrow(config_set), length(Temp_z), 3),
    dimnames = list(species_set, model_names, NULL, c("min", "mid", "max"))
  )
  # rhot_scz = array( NA, dim=c(length(species_set),nrow(config_set),2),
  #                    dimnames = list(species_set,NULL,c("hat","se")) )

  for (si in seq_along(species_set)) {
    for (ci in seq_len(nrow(config_set))) {
      #
      species <- species_set[si]
      do_space <- config_set[ci, "space"]
      do_time <- config_set[ci, "time"]
      do_diffusion <- config_set[ci, "diffusion"] #
      do_quadratic <- config_set[ci, "quadratic"]
      type <- paste0(do_quadratic, "-", do_space, do_time, do_diffusion)

      #
      run_dir <- file.path(root_dir, Date, species, type)

      if ("opt.RDS" %in% list.files(run_dir)) {
        #
        opt <- readRDS(file.path(run_dir, "opt.RDS"))
        rep <- readRDS(file.path(run_dir, "rep.RDS"))
        obj <- readRDS(file.path(run_dir, "obj.RDS"))
        if (class(opt) != "try-error") {
          results_scz[si, ci, "runtime"] <- opt$runtime
          results_scz[si, ci, "pars"] <- length(opt$par)
          results_scz[si, ci, "obj"] <- opt$obj
          results_scz[si, ci, "AIC"] <- opt$AIC
          results_scz[si, ci, "log_kappaS"] <- opt$par["log_kappaS"]
          results_scz[si, ci, "kappaT"] <- opt$par["kappaT"]
          if (!is.na(opt$SD[1]) && !is.na(opt$par["kappaT"])) {
            results_scz[si, ci, "rhoT"] <- as.list(opt$SD, what = "Estimate", report = TRUE)$rhoT
            results_scz[si, ci, "se_rhoT"] <- as.list(opt$SD, what = "Std. Error", report = TRUE)$rhoT
          }
          results_scz[si, ci, "kappaST"] <- opt$par["kappaST"]
          if ("MSD" %in% names(rep)) {
            results_scz[si, ci, "MSD"] <- rep$MSD
          }
          # if( L1version == "innovations" ){
          #  results_scz[si,ci,'MSD'] = 4 / exp(results_scz[si,ci,"log_kappaS"])^2
          #  if( !is.na(results_scz[si,ci,'MSD']) ){
          #    results_scz[si,ci,'se_log_kappaS'] = as.list(opt$SD, what = "Std. Error")$log_kappaS
          #    samp = rnorm( n = 1e6, mean = results_scz[si,ci,"log_kappaS"], sd = results_scz[si,ci,'se_log_kappaS'] )
          #    results_scz[si,ci,'RMSD_low'] = sqrt( 4 / exp(quantile(samp, prob = 0.025))^2 )
          #    results_scz[si,ci,'RMSD_high'] = sqrt( 4 / exp(quantile(samp, prob = 0.975))^2 )
          #  }
          # }else{
          #  results_scz[si,ci,'MSD'] = 4 / exp(results_scz[si,ci,"log_kappaS"])^2 * ifelse( is.na(results_scz[si,ci,"rhoT"]), 1, 1-results_scz[si,ci,"rhoT"] )
          # }
          if (!is.na(opt$SD[1]) && ("RMSD" %in% names(rep))) {
            results_scz[si, ci, "RMSD"] <- rep$RMSD
            results_scz[si, ci, "se_RMSD"] <- as.list(opt$SD, what = "Std. Error", report = TRUE)$RMSD
          }
          if (!is.na(opt$SD[1])) {
            gamma_scjz[si, ci, , "hat"] <- as.list(opt$SD, what = "Estimate")$gamma_j
            gamma_scjz[si, ci, , "se"] <- as.list(opt$SD, what = "Std. Error")$gamma_j
          }
          # Covariate response
          if (!is.na(opt$SD[1])) {
            hat <- as.list(opt$SD, what = "Estimate")$gamma_j
            Cov <- opt$SD$cov[which(names(opt$SD$value) == "gamma_j"), which(names(opt$SD$value) == "gamma_j")]
            beta_rj <- mvtnorm::rmvnorm(n = 1000, mean = hat, sigma = Cov)
            # Temp_z = seq( min(obj$env$data$Temp_s2t), max(obj$env$data$Temp_s2t), by = 0.1 )
            Temp_zr <- apply(beta_rj, MARGIN = 1, FUN = \(v) Temp_z * v[1] + Temp_z^2 * v[2])
            p_sczq[si, ci, , ] <- t(apply(Temp_zr, MARGIN = 1, FUN = quantile, prob = c(0.05, 0.5, 0.95)))
          }
        }
        # Load for later
        mesh2 <- readRDS(file = file.path(run_dir, "mesh2.RDS"))
        mesh2$loc <- mesh2$loc * rescale
        mesh2_sf <- fm_as_sfc(mesh2)
        A_gs2 <- fm_evaluator(mesh2, loc = st_coordinates(st_centroid(grid)))$proj$A #  obj$env$data$A_gs2
      }
    }
  }

  results_scz[, , "AIC"] <- t(apply(results_scz[, , "AIC"], MARGIN = 1, FUN = \(v) v - min(v, na.rm = TRUE)))
  results_scz[, , "AIC_weight"] <- t(apply(results_scz[, , "AIC"], MARGIN = 1, FUN = \(v) exp(-0.5 * v) / sum(exp(-0.5 * v), na.rm = TRUE)))
  # results_scz[,,'RMSD_high'] = ifelse( results_scz[,,'RMSD_low']>1e10, NA, results_scz[,,'RMSD_high'] )
  # results_scz[,,'RMSD_low'] = ifelse( results_scz[,,'RMSD_low']>1e10, NA, results_scz[,,'RMSD_low'] )

  # Display
  (summary <- apply(results_scz, MARGIN = 1, FUN = \(m) cbind(config_set, m)))
  capture.output(summary, file = file.path(root_dir, Date, "summary.txt"))

  #
  Table <- do.call(rbind, summary)
  Table <- cbind(species = rep(species_set, each = nrow(config_set)), Table)
  Table <- Table[, c("species", "space", "time", "runtime", "AIC", "rhoT", "MSD", "RMSD")]
  Table[, c("runtime", "AIC", "rhoT")] <- round(Table[, c("runtime", "AIC", "rhoT")], 2)
  Table[, c("MSD", "RMSD")] <- round(Table[, c("MSD", "RMSD")], 0)
  Table[, c("space", "time")] <- ifelse(Table[, c("space", "time")] == 1, "X", "")
  write.csv(Table, row.names = FALSE, file = file.path(date_dir, "Table1.csv"))

  #######  Plot rho_t
  # DF = expand.grid( dimnames(results_scz)[1:2] )
  # DF$est = as.vector(results_scz[,,"rhoT"])
  # DF$se = as.vector(results_scz[,,"se_rhoT"])
  # DF$low = DF$est - 1.96 * DF$se
  # DF$high = DF$est + 1.96 * DF$se
  # DF$AICw = as.vector(results_scz[,,"AIC_weight"])
  #
  # rightcolor = c("grey50", "blue")[2]
  # ggplot(data=DF, aes(x=model, y=est)) +
  #  geom_point() +
  #  geom_point( aes(y = AICw), col = rightcolor, shape = 17 ) + # ""
  #  facet_wrap( ~ species, ncol = 2 ) +
  #  geom_errorbar( aes(ymax=as.numeric(low),ymin=as.numeric(high)), width=0.25, position=position_dodge(0.9) ) +
  #  scale_y_continuous( limits = c(-0.05,1), name = expression(rho[t]), position = "left",
  #                      sec.axis = sec_axis(~ . * 1, name = "AIC weights") ) +
  #  theme( axis.title.y.right = element_text(color = rightcolor), axis.text.y.right  = element_text(color = "black") )
  # ggsave( file = file.path(date_dir,"rho_t.png"), width = 2 * 2, height = 2 * 3 )

  #
  buffer <- \(..., extra = 0.05){
    out <- range(..., na.rm = TRUE)
    mid <- mean(out)
    out <- mid + (1 + extra) * (out - mid)
    return(out)
  }
  png(file = file.path(date_dir, "summary.png"), width = 1.1 * length(species_set), height = 1.5 * 3, res = 200, units = "in")
  par(mfcol = c(3, length(species_set)), mgp = c(2, 0.5, 0), tck = -0.02, mar = c(0.5, 0.5, 0.5, 0.5), oma = c(3, 5, 3, 0), yaxs = "i", xaxs = "i")
  ylim2 <- buffer(results_scz[, , "rhoT"] - 1.96 * results_scz[, , "se_rhoT"], results_scz[, , "rhoT"] + 1.96 * results_scz[, , "se_rhoT"], extra = 0.2)
  # ylim3 = buffer( results_scz[,,"RMSD"], extra = 0.2 )
  # ylim3 = buffer( results_scz[,,"RMSD_low"], results_scz[,,"RMSD_high"], extra = 0.2 )
  ylim3 <- buffer(results_scz[, , "RMSD"] - 1.96 * results_scz[, , "se_RMSD"], results_scz[, , "RMSD"] + 1.96 * results_scz[, , "se_RMSD"], extra = 0.2)
  xlim <- c(0.5, 4.5)
  for (si in seq_along(species_set)) {
    # AIC weights
    plot(
      x = seq_len(nrow(config_set)), y = results_scz[si, , "AIC_weight"], xlab = "", ylab = "", xlim = xlim,
      ylim = c(0, 1.05), xaxt = "n", pch = 21, bg = "black", yaxt = "n", cex = 2
    )
    mtext(side = 3, text = species_names[si], line = 0.5)
    if (si == 1) mtext(side = 2, text = "AIC weight", line = 2)
    if (si == 1) axis(2)
    # rho_t
    plot(
      x = seq_len(nrow(config_set)), y = results_scz[si, , "rhoT"], xlab = "", ylab = "", xlim = xlim, xaxt = "n",
      ylim = ylim2, pch = 21, bg = "black", yaxt = "n", cex = 2
    )
    arrows(
      x0 = seq_len(nrow(config_set)), x1 = seq_len(nrow(config_set)),
      y0 = results_scz[si, , "rhoT"] - 1.96 * results_scz[si, , "se_rhoT"],
      y1 = results_scz[si, , "rhoT"] + 1.96 * results_scz[si, , "se_rhoT"],
      angle = 90, code = 3, length = 0.05
    )
    if (si == 1) mtext(side = 2, text = expression(rho[t]), line = 2)
    if (si == 1) axis(2)
    abline(h = 0, lty = "dotted")
    #
    plot(
      x = seq_len(nrow(config_set)), y = results_scz[si, , "RMSD"], xlab = "", ylab = "", xlim = xlim,
      xaxt = "n", ylim = ylim3, pch = 21, bg = "black", yaxt = "n", cex = 2
    )
    arrows(
      x0 = seq_len(nrow(config_set)), x1 = seq_len(nrow(config_set)),
      # y0 = results_scz[si,,"RMSD_low"],
      # y1 = results_scz[si,,"RMSD_high"],
      y0 = results_scz[si, , "RMSD"] - 1.96 * results_scz[si, , "se_RMSD"],
      y1 = results_scz[si, , "RMSD"] + 1.96 * results_scz[si, , "se_RMSD"],
      angle = 90, code = 3, length = 0.05
    )
    axis(side = 1, at = seq_len(nrow(config_set)), labels = model_names)
    if (si == 1) mtext(side = 2, text = "Root-mean squared\ndisplacement (km)", line = 2)
    if (si == 1) axis(2)
    abline(h = 0, lty = "dotted")
  }
  mtext(side = 1, text = c("model"), outer = TRUE, line = 1.5)
  dev.off()

  # Marginal responses
  png(file = file.path(date_dir, "marginal_response.png"), width = 2 * 2, height = 2 * 3, res = 200, units = "in")
  par(mfrow = c(3, 2), mgp = c(2, 0.5, 0), tck = -0.02, mar = c(2, 2, 1, 1), yaxs = "i", xaxs = "i", oma = c(2, 2, 0, 0))
  for (si in seq_along(species_set)) {
    #
    ci <- which.min(summary[[si]][, "AIC"])
    y_czq <- p_sczq[si, , , ]
    # y_czq = sweep( y_czq, MARGIN = 1, FUN = "-", STAT = apply(y_czq[,,'mid'], MARGIN=1, FUN=max) )
    y_czq <- exp(y_czq)
    if (length(ci) > 0) {
      matplot(
        x = Temp_z, y = t(y_czq[c(1, ci), , "mid"]), lwd = 2, col = c("black", "blue"),
        lty = "solid", type = "l", ylim = range(y_czq[c(1, ci), , ]), ylab = "", xlab = "",
        main = species_set[si], log = "y"
      )
      polygon(x = c(Temp_z, rev(Temp_z)), y = c(y_czq[1, , "min"], rev(y_czq[1, , "max"])), col = rgb(0, 0, 0, 0.2), border = NA)
      polygon(x = c(Temp_z, rev(Temp_z)), y = c(y_czq[ci, , "min"], rev(y_czq[ci, , "max"])), col = rgb(0, 0, 1, 0.2), border = NA)
    } else {
      plot.new()
    }
  }
  mtext(side = 1:2, text = c(
    ifelse(covariate == "temp", "Seafloor temperature (C)", "Seafloor temperature anomaly (C)"),
    "Density relative to average"
  ), outer = TRUE)
  legend("bottomright", bty = "n", legend = c("null", "selected"), fill = c("black", "blue"), title = "model")
  dev.off()

  #
  if (covariate == "anomaly") {
    Temp_gt <- A_gs2 %*% obj$env$data$Temp_s2t
    dimnames(Temp_gt) <- list(NULL, year = year_set)
    maxval <- max(abs(Temp_gt))
    pal <- colorRampPalette(c("blue", "white", "red"))
    png(file = file.path(date_dir, "Temp_gt.png"), width = 10, height = 10, res = 200, units = "in")
    par(mfrow = c(7, 7), oma = c(3, 3, 0, 0))
    plot_grid <- st_sf(grid, as.matrix(Temp_gt))
    for (t in seq_len(ncol(Temp_gt))) {
      plot(plot_grid[, t],
        border = NA, key.pos = NULL, reset = FALSE,
        breaks = seq(-1 * maxval, maxval, length = 11), pal = pal(11)
      ) # , zlim = range(cbind(Temp_gt))
      box()
      plot(st_geometry(worldmap), add = TRUE, col = "grey", border = NA)
      plot(st_geometry(survey_domain), add = TRUE, border = "black")
      if (t == 1) {
        add_legend(legend = round(c(-maxval, 0, maxval), 1), col = pal(10), legend_x = c(0.25, 0.3), legend_y = c(0.05, 0.45))
      }
    }
    mtext(side = 1:2, text = c("Longitude", "Latitude"), outer = TRUE)
    dev.off()
  }

  # Plot covariate diffusion from point
  # si = 3
  xlim <- 650 + c(-1, 1) * 200
  ylim <- 6500 + c(-1, 1) * 200
  MSD_sz <- array(NA, dim = c(length(species_set), 4))
  pred_s <- rep(NA, length(species_set))
  png(file = file.path(date_dir, "point_diffusion.png"), width = 2 * 4, height = 2 * 6, res = 200, units = "in")
  par(mfrow = c(length(species_set), 4), mgp = c(2, 0.5, 0), tck = -0.02, yaxs = "i", xaxs = "i", oma = c(1, 2, 3, 1))
  for (si in seq_along(species_set)) {
    #
    ci <- which.min(summary[[si]][, "AIC"])
    if (length(ci) > 0) {
      species <- species_set[si]
      do_space <- config_set[ci, "space"]
      do_time <- config_set[ci, "time"]
      do_diffusion <- config_set[ci, "diffusion"] #
      do_quadratic <- config_set[ci, "quadratic"]
      type <- paste0(do_quadratic, "-", do_space, do_time, do_diffusion)
      run_dir <- file.path(root_dir, Date, species, type)
      rep <- readRDS(file.path(run_dir, "rep.RDS"))
      obj <- readRDS(file.path(run_dir, "obj.RDS"))
      opt <- readRDS(file.path(run_dir, "opt.RDS"))
      mesh2 <- readRDS(file = file.path(run_dir, "mesh2.RDS"))
      mesh2$loc <- mesh2$loc * rescale
      mesh2_sf <- fm_as_sfc(mesh2)
      A_gs2 <- fm_evaluator(mesh2, loc = st_coordinates(st_centroid(grid)))$proj$A #  obj$env$data$A_gs2
      loc <- sweep(mesh2$loc, STAT = colMeans(mesh2$loc), MARGIN = 2, FUN = "-")
      which_s2 <- which.min(rowSums(loc^2))
      x1_s2t <- x0_s2t <- array(0, dim = dim(rep$T_s2t))
      x0_s2t[which_s2, 1] <- 1 / sum(A_gs2[, which_s2])
      x1_s2t[] <- solve(rep$IminusP_k2k2, as.vector(x0_s2t))
      x0_gt <- A_gs2 %*% x0_s2t
      x1_gt <- A_gs2 %*% x1_s2t
      mat <- as.matrix(cbind(x0_gt[, 1], x1_gt[, 1:3]))
      # mat = ifelse( mat < 1e-3, NA, mat )
      plot_grid <- st_sf(grid, mat)
      for (i in 1:4) {
        # par( mar = c(2,2,2,1) )
        # plot(1, type = "n", xlim = xlim, ylim = ylim, xaxt = "n", yaxt = "n", xlab = "", ylab = "") # , breaks = c(-0.01, seq(0.1, 1, by = 0.1)) )
        # plot(plot_grid[,i], border = NA, key.pos=NULL, reset=FALSE, main = "", xlim = xlim, ylim = ylim, add = TRUE) # , breaks = c(-0.01, seq(0.1, 1, by = 0.1)) )
        plot(plot_grid[, i], border = NA, key.pos = NULL, reset = FALSE, main = "", xlim = xlim, ylim = ylim) # , breaks = c(-0.01, seq(0.1, 1, by = 0.1)) )
        if (i == 1) mtext(side = 2, text = species, line = 0.5)
        if (si == 1) mtext(side = 3, text = c("original", "diffused", "lag-1", "lag-2")[i], line = 1.5)
        plot(st_geometry(worldmap), add = TRUE, col = "grey")
        box()
        plot(mesh2_sf, add = TRUE, border = "grey50")
        # Calculate MSD
        locs <- st_coordinates(st_centroid(grid)) / rescale
        mean_z <- apply(locs, w = mat[, i], MARGIN = 2, FUN = weighted.mean)
        devs <- sweep(locs, MARGIN = 2, STAT = mean_z, FUN = "-")
        MSD_sz[si, i] <- sum(apply(devs^2, w = mat[, i], MARGIN = 2, FUN = weighted.mean))
        kappaS <- exp(opt$par["log_kappaS"])
        pred_s[si] <- 4 / kappaS^2
        #
        # title( paste0("total = ", round(sum(mat[,i]),3)) )
        # title( paste0(round(sum(mat[,i]),3), ", ", round(MSD_sz[si,i] * rescale,0) ) )
        title(paste0("total: ", round(sum(mat[, i]), 3), " | ", "MSD: ", round(MSD_sz[si, i] * rescale, 1)))
        if (si == length(species_set)) axis(1)
        if (i == 4) axis(4)
      }
    }
  }
  dev.off()

  # Compare spatial response
  # which_year = 2014
  T_scgt <- gamma_scgt <- array(NA,
    dim = c(length(species_set), nrow(config_set), length(grid), length(year_set)),
    dimnames = list(species_set, NULL, NULL, year_set)
  )
  for (si in seq_along(species_set)) {
    for (ci in seq_len(nrow(config_set))) {
      #
      species <- species_set[si]
      do_space <- config_set[ci, "space"]
      do_time <- config_set[ci, "time"]
      do_diffusion <- config_set[ci, "diffusion"] #
      do_quadratic <- config_set[ci, "quadratic"]
      type <- paste0(do_quadratic, "-", do_space, do_time, do_diffusion)
      run_dir <- file.path(root_dir, Date, species, type)

      if ("opt.RDS" %in% list.files(run_dir)) {
        rep <- readRDS(file.path(run_dir, "rep.RDS"))
        obj <- readRDS(file.path(run_dir, "obj.RDS"))
        # Linear predictor contribution
        gamma_s2t <- gamma_scjz[si, ci, 1, "hat"] * rep$T_s2t + gamma_scjz[si, ci, 2, "hat"] * (rep$T_s2t)^2
        gamma_gt <- as.matrix(A_gs2 %*% gamma_s2t)
        # Effective temperature
        T_gt <- as.matrix(A_gs2 %*% rep$T_s2t)
        #
        colnames(gamma_gt) <- colnames(gamma_gt) <- year_set
        gamma_scgt[si, ci, , ] <- gamma_gt
        T_scgt[si, ci, , ] <- T_gt

        # Not necessary to do each individual plot
        # png( file = file.path(run_dir,"gamma_gt.png"), width = 10, height = 10, res = 200, units = "in")
        #  plot_grid = st_sf( grid, as.matrix(gamma_gt) )
        #  plot(plot_grid, border = NA) # , zlim = range(cbind(gamma_gt,T_gt)) )
        # dev.off()
      }
    }
  }

  # Plot single year across models
  which_year <- "2014"
  for (zi in 1:2) {
    png(file = file.path(date_dir, paste0("compare_", which_year, "_", c("gamma", "T")[zi], ".png")), width = 2 * nrow(config_set), height = 2 * length(species_set), res = 200, units = "in")
    par(mfrow = c(length(species_set), nrow(config_set)), oma = c(0, 3, 2, 0))
    for (si in seq_along(species_set)) {
      for (ci in seq_len(nrow(config_set))) {
        if (zi == 1) {
          vec <- gamma_scgt[si, ci, , which_year]
          vec <- vec - max(vec)
        } else {
          vec <- T_scgt[si, ci, , which_year]
        }
        if (any(!is.na(vec))) {
          if (zi == 1) {
            breaks <- log(seq(0, 1, by = 0.1))
            pal <- viridisLite::viridis
          } else {
            maxval <- max(abs(T_scgt[si, , , which_year]), na.rm = TRUE)
            breaks <- seq(-maxval, maxval, length = 11)
            pal <- colorRampPalette(c("blue", "white", "red"))
          }
          plot_grid <- st_sf(grid, vec)
          plot(plot_grid, border = NA, key.pos = NULL, reset = FALSE, main = "", breaks = breaks, pal = pal) # , zlim = range(cbind(gamma_gt,T_gt)) )
          legend("topright", legend = round(results_scz[si, ci, "AIC"], 2), bty = "n")
          plot(st_geometry(worldmap), add = TRUE, col = "grey", border = NA)
          plot(st_geometry(survey_domain), add = TRUE, border = "black")
        } else {
          plot.new()
        }
        if (zi == 2 & ci == 1) {
          add_legend(legend = round(c(-maxval, 0, maxval), 1), col = pal(10), legend_x = c(0.2, 0.25), legend_y = c(0.05, 0.45))
        }
        if (ci == 1) mtext(side = 2, text = species_names[si])
        if (si == 1) mtext(side = 3, text = model_names[ci])
        box()
      }
    }
    mtext(side = 2:3, outer = TRUE, text = c("species", "model"), cex = 1.5, line = c(1, 0.5))
    dev.off()
  }

  # selected model across years
  # Plot single year across models
  which_years <- as.character(2013:2017)
  # which_years = as.character(2012:2016)
  for (zi in 1:2) {
    png(file = file.path(date_dir, paste0("compare_", min(which_years), "-", max(which_years), "_", c("gamma", "T")[zi], ".png")), width = 2 * length(which_years), height = 2 * length(species_set), res = 200, units = "in")
    par(mfrow = c(length(species_set), length(which_years)), oma = c(0, 3, 2, 0))
    for (si in seq_along(species_set)) {
      for (ti in seq_along(which_years)) {
        ci <- which.min(summary[[si]][, "AIC"])
        which_year <- which_years[ti]
        if (zi == 1) {
          vec <- gamma_scgt[si, ci, , which_year]
          vec <- vec - max(vec)
        } else {
          vec <- T_scgt[si, ci, , which_year]
        }
        if (any(!is.na(vec))) {
          if (zi == 1) {
            breaks <- log(seq(0, 1, by = 0.1))
            pal <- viridisLite::viridis
          } else {
            maxval <- max(abs(T_scgt[si, ci, , which_years]), na.rm = TRUE)
            breaks <- seq(-maxval, maxval, length = 11)
            pal <- colorRampPalette(c("blue", "white", "red"))
          }
          plot_grid <- st_sf(grid, vec)
          plot(plot_grid, border = NA, key.pos = NULL, reset = FALSE, main = "", breaks = breaks, pal = pal) # , zlim = range(cbind(gamma_gt,T_gt)) )
          legend("topright", legend = round(results_scz[si, ci, "AIC"], 2), bty = "n")
          plot(st_geometry(worldmap), add = TRUE, col = "grey", border = NA)
          plot(st_geometry(survey_domain), add = TRUE, border = "black")
        } else {
          plot.new()
        }
        if (zi == 2 & ti == 1) {
          add_legend(legend = round(c(-maxval, 0, maxval), 1), col = pal(10), legend_x = c(0.2, 0.25), legend_y = c(0.05, 0.45))
        }
        if (ti == 1) mtext(side = 2, text = species_names[si])
        if (si == 1) mtext(side = 3, text = which_year)
        box()
      }
    }
    mtext(side = 2:3, outer = TRUE, text = c("species", "year"), cex = 1.5, line = c(1, 0.5))
    dev.off()
  }

  # Plot effective temperature for each species-model
  # VERY SLOW
  if (FALSE) {
    for (si in seq_along(species_set)) {
      for (ci in seq_len(nrow(config_set))) {
        #
        species <- species_set[si]
        do_space <- config_set[ci, "space"]
        do_time <- config_set[ci, "time"]
        do_diffusion <- config_set[ci, "diffusion"] #
        do_quadratic <- config_set[ci, "quadratic"]
        type <- paste0(do_quadratic, "-", do_space, do_time, do_diffusion)
        run_dir <- file.path(root_dir, Date, species, type)

        if ("opt.RDS" %in% list.files(run_dir)) {
          rep <- readRDS(file.path(run_dir, "rep.RDS"))
          mesh2 <- readRDS(file = file.path(run_dir, "mesh2.RDS"))
          obj <- readRDS(file = file.path(run_dir, "obj.RDS"))
          mesh2$loc <- mesh2$loc * rescale
          A_gs2 <- fm_evaluator(mesh2, loc = st_coordinates(st_centroid(grid)))$proj$A #  obj$env$data$A_gs2
          T_gt <- as.matrix(A_gs2 %*% rep$T_s2t)
          Temp_gt <- as.matrix(A_gs2 %*% obj$env$data$Temp_s2t)
          colnames(Temp_gt) <- colnames(T_gt) <- year_set
          mean(apply(T_gt, MARGIN = 1, FUN = sd))
          mean(apply(Temp_gt, MARGIN = 1, FUN = sd))
          png(file = file.path(run_dir, "T_gt.png"), width = 10, height = 10, res = 200, units = "in")
          plot_grid <- st_sf(grid, T_gt)
          plot(plot_grid, border = NA) # , zlim = range(cbind(Temp_gt,T_gt)) )
          dev.off()
        }
      }
    }
  }
}
