#
root_dir <- R'(C:\Users\james\OneDrive\Desktop\Work files (backup)\Collab-2025\2025 -- spacetime distributed lag)'

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
dir.create(date_dir)

# species
# Not available: "yellowfin sole"
# Not enough data: "sablefish"
# On boundary: "eulachon"
print(get_species(), n = 100)
species_set <- c("pacific cod", "walleye pollock", "capelin", "pacific herring", "arrowtooth flounder", "pacific halibut")

#
if (!("settings.RDS" %in% list.files(date_dir))) {
  config_list <- list(
    space = 0:1,
    time = 0:1,
    diffusion = 0,
    quadratic = 1
  )
  config_set <- expand.grid(config_list)
  config_set <- config_set[which(config_set$space == 1 | config_set$diffusion == 0), ]

  #
  rescale <- 1
  cutoff1 <- 40
  cutoff2 <- 40
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
} else {
  settings <- readRDS(file = file.path(date_dir, "settings.RDS"))
  attach(settings)
}

# survey domain from VAST
survey_domain <- st_read(file.path(R'(C:\Users\james\OneDrive\Desktop\Git\VAST\inst\region_shapefiles\EBSshelf)', "EBSshelf.shp"))
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
