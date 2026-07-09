

library(fmesher)
library(Matrix)
library(pracma)
library(sf)

set.seed(123)

#
n_x = 3200   # Number of samples
n_t = 3

# Spatial correlation
kappa_S = 10

# Temporal autocorrelation
kappa_t = c(0, 0.5)[2]
rho_T = kappa_t / (1 + kappa_t)

# Spatio-temporal diffusion parameter
kappa_st = -kappa_t
# If rho_st > 0, then negative values at release point for transformed point-mass
# If rho_st < -rho_t, then negative values at boundaries for transformed point-mass

# simulate locations from Poisson or Poisson-disk
loc = poisson2disk( n=n_x, a = 1, b = 1 )

# Make mesh and FEM matrices
boundary = fm_segm(
  rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)),
  is.bnd = TRUE
)
mesh = fm_mesh_2d(
  loc,
  boundary = boundary,
  #refine = TRUE,
  loc.domain = cbind(c(0,0,1,1,0),c(0,1,1,0,0))
)
n_s = mesh$n
n_k = n_s * n_t

# Make FEM matrices
spde = fm_fem( mesh )
summary(spde$g1@x)
invM0 = Diagonal( n = n_s, x = 1/diag(spde$c0) )
I_ss = Diagonal( n = n_s )
I_tt = Diagonal( n = n_t )
I_kk = kronecker( I_tt, I_ss )

# Adjacency matrix
A_ss = as(spde$g1, "nsparseMatrix")

# Distance matrix (same pattern as A_ss)
Atriplet = mat2triplet(A_ss)
#
sf_grid = st_make_grid( st_multipoint(cbind(c(0,0,1,1,0),c(0,1,1,0,0))), n=100 )
coords_gz = st_coordinates( st_centroid(sf_grid) )
A_gs = fm_evaluator( mesh, loc=coords_gz )$proj$A

# Rescale so that sum(A_gs) = geoscale^2
A_gs = A_gs * ( 1 / nrow(A_gs) )

# inverse-diffusion
# rowSums(G_ss) = 0
G_ss = -1 * kappa_S^(-2) * invM0 %*% spde$g1
# invD_ss = I_ss - G_ss

# Joint operator
P_kk = Diagonal( n = n_k, x = 0 )
P_kk = P_kk + kronecker( I_tt, G_ss )
if( n_t > 1 ){
  L_tt = bandSparse( n = n_t, k = -1 ) - Diagonal( n = n_t )
  P_kk = P_kk + kappa_t * kronecker( L_tt, I_ss )
  P_kk = P_kk + kappa_st * kronecker( L_tt, G_ss )
}

# Check with midpoint of domain
which_mid = which.min( rowSums(scale(mesh$loc[,1:2])^2) )
d0_st = matrix(0, nrow = n_s, ncol = n_t )
d0_st[which_mid,1] = 1/sum(A_gs[,which_mid])
d0_gt = A_gs %*% d0_st
colnames(d0_gt) = paste( "orig", seq_len(n_t) )

#
d1_k = solve( I_kk - P_kk, as.vector(d0_st) )
d1_st = array(d1_k, dim=dim(d0_st) )
d1_gt = A_gs %*% d1_st
colnames(d1_gt) = paste( "proj", seq_len(n_t) )

# Plot initial
stuff = st_sf(
  sf_grid,
  as.matrix(d0_gt) / max(d0_gt),
  as.matrix(d1_gt) / max(d1_gt)
)   # , log(as.numeric(vec2b)))
plot( 
  stuff,
  breaks = seq(0,1,length=11), 
  #logz = TRUE,
  cex=2, pch=19, border=NA 
)

# Check for volume in original
colSums(d0_gt)

# Check for volume in STDL transform ... should decay at rho_T
colSums(d1_gt)

# Calculate centroid ... shouldn't change given diffusion
# And calculate diffusion ~propto~ MSD
Mean_tz = array(NA, dim = c(n_t,2), dimnames = list(time = seq_len(n_t),coord = c("mean_X","mean_Y")) )
MSD_t = rep(NA, n_t)
for(t in seq_len(n_t) ){
  Mean_tz[t,] = apply(coords_gz, MARGIN=2, FUN=weighted.mean, w=d1_gt[,t])
  MSD_t[t] = weighted.mean( rowSums(sweep(coords_gz, MARGIN = 2, FUN = "-", STATS = Mean_tz[t,])^2), w=d1_gt[,t] )
}
cbind( Mean_tz, MSD = MSD_t )

# Check for negative values ... shouldn't be any if -kappaT < kappa_ST < 0
apply( d1_gt, MARGIN = 2, FUN = min )


