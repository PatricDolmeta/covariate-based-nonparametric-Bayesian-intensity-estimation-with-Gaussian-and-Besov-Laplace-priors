# =============================================================================
# Naive preconditioned Crank–Nicolson (pCN) MCMC
# for Bayesian inference of an Inhomogeneous Poisson Process (IPP)
# intensity function, represented in a wavelet basis.
#
# The intensity is modelled as either:
#   exponential link:  rho(z) = exp( theta^T W(z) )
#   logistic link:     rho(z) = rho_star * sigmoid( theta^T W(z) )
#
# Wavelet coefficients theta are given a Besov-space prior:
#   theta_j = tau * j^{-d/2 - alpha} * T(xi_j),   xi_j ~ N(0,1) or Laplace
# where T is the identity (Gaussian prior) or the whitening map T_white
# (Laplace prior).  alpha controls smoothness, tau controls amplitude.
# =============================================================================

library(pracma)      # numerical math utilities
library(cubature)    # numerical integration
library(rmutil)      # rlaplace()
library(purrr)       # functional programming helpers
library(furrr)       # parallel purrr
library(pals)        # colour palettes (used elsewhere)
library(tictoc)      # timing
library(future)      # parallel backends
library(future.apply)
library(utils)

# -----------------------------------------------------------------------------
# T_white: whitening / tail-probability transform
#
# Maps a real value x to sign(x) * (-log P(|Z| > |x|)) where Z ~ N(0,1).
# This turns standard-normal white noise into Laplace white noise:
# if xi ~ N(0,1) then T_white(xi) ~ Laplace(0,1).
# Used when Laplace = TRUE to build a heavier-tailed Besov prior on theta.
# -----------------------------------------------------------------------------
T_white <- function(x) {
  # log P(Z > |x|) = log of the upper-tail normal probability
  tail_log <- pnorm(abs(x), lower.tail = FALSE, log.p = TRUE)
  # -log(2 * P(Z > |x|)) = -log2 - tail_log  (log of half-normal CDF inversion)
  sign(x) * (-log(2) - tail_log)
}

# Small helper: return 0 for zero-length vectors (avoids sum() returning numeric(0))
zero_if_empty <- function(x) if (length(x) == 0) 0 else x

# -----------------------------------------------------------------------------
# lik_ratio: Metropolis–Hastings acceptance probability
#
# Computes min(1, r) where r is the MH ratio for the joint update of:
#   - the wavelet coefficients theta  (via the pCN proposal on xi)
#   - optionally alpha and/or tau     (via a random-walk proposal)
#
# Log-likelihood of an IPP:
#   log L(theta) = sum_{events i} log rho(z_i) - integral rho(z) dz
#
# Arguments
#   theta_prop / theta_kl : proposed and current wavelet coefficient vectors
#   alpha_prop / alpha    : proposed and current smoothness hyperparameters
#   tau_prop   / tau      : proposed and current scale hyperparameters
#   Wz_mat                : (n_basis x n_pixels)  wavelet basis evaluated on
#                           the discretisation grid (used for the integral)
#   Wzi_mat               : (n_basis x n_events)  wavelet basis at event locations
#   Full_basis            : row indices of non-zero basis functions on the grid
#   Full_loc              : row indices of non-zero basis functions at events
#   rho_star              : upper bound on intensity (logistic link only)
#   link                  : "exponential" or "logistic"
#   area_per_pixels       : area of one pixel (scalar; used for Riemann sum)
# -----------------------------------------------------------------------------
lik_ratio <- function(theta_prop, theta_kl, alpha_prop, alpha, tau_prop, tau,
                      Wz_mat, Wzi_mat, Full_basis, Full_loc, rho_star, link, area_per_pixels){
  
  if(identical(link, "exponential")){
    # log-intensity at event locations: log rho(z_i) = theta^T W(z_i)
    log_rho_new <- theta_prop[Full_loc] %*% Wzi_mat
    log_rho_old <- theta_kl[Full_loc]  %*% Wzi_mat
    
    # Integrated intensity: integral rho dz ≈ sum_pixels exp(theta^T W(z_pixel)) * area
    I_new <- sum(exp(theta_prop[Full_basis] %*% Wz_mat) * area_per_pixels)
    I_old <- sum(exp(theta_kl[Full_basis]   %*% Wz_mat) * area_per_pixels)
    
  } else {
    # Logistic (sigmoid) link: rho(z) = rho_star * sigma(theta^T W(z))
    log_rho_new <- log(rho_star * 1/(1 + exp(-theta_prop[Full_loc] %*% Wzi_mat)))
    log_rho_old <- log(rho_star * 1/(1 + exp(-theta_kl[Full_loc]  %*% Wzi_mat)))
    
    I_new <- sum(rho_star * 1/(1 + exp(-theta_prop[Full_basis] %*% Wz_mat)) * area_per_pixels)
    I_old <- sum(rho_star * 1/(1 + exp(-theta_kl[Full_basis]  %*% Wz_mat)) * area_per_pixels)
  }
  
  # MH log-ratio:
  #   [log L(theta_prop) - log L(theta_kl)]          <- IPP likelihood ratio
  #   + [log pi(alpha_prop) - log pi(alpha)]          <- Exp(1) prior on alpha
  #   + [log pi(tau_prop)   - log pi(tau)]            <- Exp(1) prior on tau
  #
  # The pCN proposal is designed so that the Gaussian/Laplace prior on xi
  # cancels exactly, leaving only the likelihood terms for the xi update.
  # The prior terms for alpha and tau appear when those are being updated.
  return(min(1, exp(
    sum(log_rho_new - log_rho_old) - I_new + I_old +
      dexp(alpha_prop, log = TRUE) - dexp(alpha, log = TRUE) +
      dexp(tau_prop,   log = TRUE) - dexp(tau,   log = TRUE)
  )))
}

# =============================================================================
# MCMC: main sampler
#
# Parameters
#   run_name      : string prefix for output folder
#   Laplace       : TRUE  → Laplace/Besov prior on theta via T_white
#                   FALSE → Gaussian prior (T = identity)
#   niter         : number of MCMC iterations
#   loc           : spatstat ppp object of event locations
#   covs          : list of im (pixel image) objects, one per covariate
#   discretization: (J x d) matrix of covariate grid points
#   window        : spatstat owin observation window
#   beta          : pCN step size in (0,1);  smaller = higher acceptance
#   tau           : initial scale hyperparameter (NULL → data-driven default)
#   alpha         : initial smoothness hyperparameter
#   shape, rate   : Gamma hyperparameters for rho_star (logistic link)
#   link          : "exponential" or "logistic"
#   rho_star      : intensity upper bound (logistic link)
#   Wavelet_basis : named list of wavelet basis images (from the imwd pipeline)
#   L             : number of wavelet resolution levels to include
#   refl          : TRUE → symmetric boundary; FALSE → periodic boundary
#   adaptA        : TRUE → also sample alpha via random walk
#   b_a           : random-walk SD for alpha proposals
#   adaptT        : TRUE → also sample tau via random walk
#   b_t           : random-walk SD for tau proposals
#   index_dom     : (pre-computed) covariate→discretisation indices for the grid
#   index_loc     : (pre-computed) covariate→discretisation indices for events
# =============================================================================
MCMC <- function(run_name, Laplace = TRUE, niter,
                 loc, covs,
                 discretization, window,
                 beta, tau = NULL, alpha,
                 shape, rate, link, rho_star,
                 Wavelet_basis, L, refl = FALSE,
                 adaptA = FALSE, b_a = NULL,
                 adaptT = FALSE, b_t = NULL,
                 index_dom = NA, index_loc = NA){
  
  # ---- Dimensions -----------------------------------------------------------
  J <- dim(discretization)[1]   # total number of covariate grid points
  d <- dim(discretization)[2]   # covariate dimension
  
  # ---- Identify which wavelet levels to use ---------------------------------
  # Wavelet_basis names encode the resolution level as the leading integer,
  # e.g. "w0Lconstant", "w1L1", "w2L3", ...
  # We extract that integer and keep only levels strictly below L.
  name_vec <- names(Wavelet_basis)
  levels   <- as.integer(sub("^[A-Za-z](\\d+).*", "\\1", name_vec))
  coeff    <- levels[levels < L] + 1   # +1: coeff index starts at 1 for level 0
  L_tot    <- length(coeff)            # total number of basis functions retained
  
  # ---- Initialise hyperparameters -------------------------------------------
  # Default tau: scaling that balances the Besov-norm with the window area
  if (is.null(tau)) {
    tau <- area.owin(window)^(-d / ((2 - Laplace) * (2 * alpha + d)))
  }
  # If adapting alpha or tau, start from a draw from their Exp(1) prior
  if (adaptA) alpha <- rexp(1)
  if (adaptT) tau   <- rexp(1)
  
  # ---- Covariate bookkeeping ------------------------------------------------
  # Covariate values at event locations (one list element per covariate image)
  Z_i <- lapply(covs, function(cv) cv[loc])
  
  area           <- area.owin(window)
  pixels         <- length(covs[[1]]$v)   # total pixels in the observation window
  area_per_pixels <- area / pixels         # pixel area for Riemann-sum integration
  
  # For every pixel in the window: find the nearest covariate grid point
  # (nearest-neighbour lookup in the J × d discretisation)
  index_dom <- apply(
    do.call(cbind, lapply(covs, function(cv) na.omit(as.vector(cv$v)))),
    1,
    function(w) which.min(
      rowSums((discretization -
                 matrix(w, nrow = nrow(discretization), ncol = length(covs), byrow = TRUE))^2)
    )
  )
  
  # Same lookup but for the covariate values at observed event locations
  index_loc <- apply(
    do.call(cbind, Z_i),
    1,
    function(w) which.min(
      rowSums((discretization -
                 matrix(w, nrow = nrow(discretization), ncol = length(Z_i), byrow = TRUE))^2)
    )
  )
  
  # ---- Build wavelet design matrices ----------------------------------------
  # Wz_mat  : each row is one basis function evaluated at every grid pixel
  #           shape: (L_tot x n_pixels);  used for the integrated intensity
  Wz_list <- lapply(Wavelet_basis[1:L_tot], function(w) w[index_dom])
  Wz_mat  <- do.call(rbind, Wz_list)
  
  # Keep only basis functions that are non-zero somewhere on the grid
  # (many wavelets at fine scales have zero support in the covariate range)
  Full_basis <- which(apply(Wz_mat, 1, function(x) any(x != 0)))
  Wz_mat     <- Wz_mat[Full_basis, ]
  assign("Full_basis", Full_basis, envir = .GlobalEnv)  # expose for diagnostics
  
  # Wzi_mat : same basis functions evaluated at event locations
  #           shape: (L_tot x n_events);  used for the log-likelihood sum
  Wzi_list <- lapply(Wavelet_basis[1:L_tot], function(w) w[index_loc])
  Wzi_mat  <- do.call(rbind, Wzi_list)
  Full_loc <- which(apply(Wzi_mat, 1, function(x) any(x != 0)))
  Wzi_mat  <- Wzi_mat[Full_loc, ]
  
  # ---- Initialise the white-noise vector xi and coefficients theta ----------
  # Besov prior:  theta_j = tau * j^{-d/2 - alpha} * T(xi_j)
  #   j^{-d/2 - alpha}  encodes the Sobolev/Besov decay rate
  #   T(xi_j) is either xi_j (Gaussian) or T_white(xi_j) (Laplace)
  if (Laplace) {
    T_func <- T_white
    xi_kl  <- rlaplace(L_tot, m = 0, s = 1)   # Laplace white noise
  } else {
    T_func <- function(x) x                    # identity: Gaussian prior
    xi_kl  <- rnorm(L_tot, 0, 1)
  }
  
  theta_kl <- tau * coeff^(-d/2 - alpha) * T_func(xi_kl)
  
  N <- npoints(loc)   # number of observed events
  
  # ---- Output folder and streaming CSV files --------------------------------
  folder <- paste(run_name, Laplace, refl, N, L, sep = "_")
  if (!dir.exists(folder)) dir.create(folder)
  setwd(folder)
  
  # Open append-mode connections so we stream results without holding
  # the entire chain in memory (important for long runs)
  file.create(description = "Post_intens.csv")   # posterior theta draws
  file.create(description = "Post_r_star.csv")   # posterior rho_star draws
  file.create(description = "Post_alpha.csv")    # posterior alpha draws
  file.create(description = "Post_tau.csv")      # posterior tau draws
  
  intens <- file(description = "Post_intens.csv", open = "a")
  rstar  <- file(description = "Post_r_star.csv", open = "a")
  alf    <- file(description = "Post_alpha.csv",  open = "a")
  Ta     <- file(description = "Post_tau.csv",    open = "a")
  
  # ---- Acceptance counters --------------------------------------------------
  Acc <- 0   # pCN proposals accepted (theta / xi update)
  Bcc <- 0   # random-walk proposals accepted for alpha
  Ccc <- 0   # random-walk proposals accepted for tau
  times <- c()
  
  # ==========================================================================
  # Main MCMC loop
  # ==========================================================================
  for(i in 1:niter){
    start <- proc.time()[3]
    
    # ------------------------------------------------------------------
    # Step 1: pCN proposal for xi (→ theta)
    #
    # The preconditioned Crank–Nicolson proposal is:
    #   xi_prop = sqrt(1 - beta^2) * xi_kl + beta * eta,   eta ~ N(0,I)
    #
    # This is prior-invariant: if xi_kl ~ N(0,I) then xi_prop ~ N(0,I).
    # Consequently the Gaussian/Laplace prior on xi cancels in the MH ratio,
    # and only the likelihood ratio needs to be evaluated — giving better
    # mixing than a naive random walk in high dimensions.
    # ------------------------------------------------------------------
    xi_prop   <- sqrt(1 - beta^2) * xi_kl + beta * rnorm(L_tot, 0, 1)
    theta_prop <- tau * coeff^(-d/2 - alpha) * T_func(xi_prop)
    
    lr <- lik_ratio(theta_prop, theta_kl,
                    alpha, alpha, tau, tau,       # alpha and tau unchanged here
                    Wz_mat, Wzi_mat, Full_basis, Full_loc,
                    rho_star, link, area_per_pixels)
    
    if (runif(1, 0, 1) < lr) {
      xi_kl    <- xi_prop    # accept: update white-noise vector
      theta_kl <- theta_prop # and derived coefficients
      Acc      <- Acc + 1
    }
    
    # Write current theta to disk (one row per iteration)
    write.table(t(theta_kl), file = intens, sep = ',', append = TRUE,
                quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    # ------------------------------------------------------------------
    # Step 2 (optional): random-walk Metropolis for alpha
    #
    # alpha controls the Besov smoothness decay: larger alpha → smoother
    # intensity functions.  Prior: alpha ~ Exp(1).
    # Proposal: alpha_new = alpha + b_a * N(0,1)   (reflected at 0)
    # ------------------------------------------------------------------
    if (adaptA) {
      alpha_new  <- alpha + b_a * rnorm(1)
      theta_temp <- tau * coeff^(-d/2 - alpha_new) * T_func(xi_kl)
      
      if (alpha_new > 0) {   # enforce positivity of smoothness parameter
        lr <- lik_ratio(theta_temp, theta_kl,
                        alpha_new, alpha, tau, tau,
                        Wz_mat, Wzi_mat, Full_basis, Full_loc,
                        rho_star, link, area_per_pixels)
        if (runif(1, 0, 1) < lr) {
          alpha    <- alpha_new
          theta_kl <- theta_temp
          Bcc      <- Bcc + 1
        }
      }
    }
    
    write.table(alpha, file = alf, sep = ',', append = TRUE,
                quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    # ------------------------------------------------------------------
    # Step 3 (optional): random-walk Metropolis for tau
    #
    # tau is an overall amplitude/scale parameter.  Prior: tau ~ Exp(1).
    # Proposal: tau_new = tau + b_t * N(0,1)   (reflected at 0)
    # ------------------------------------------------------------------
    if (adaptT) {
      tau_new    <- tau + b_t * rnorm(1)
      theta_temp <- tau_new * coeff^(-d/2 - alpha) * T_func(xi_kl)
      
      if (tau_new > 0) {     # enforce positivity of scale parameter
        lr <- lik_ratio(theta_temp, theta_kl,
                        alpha, alpha, tau_new, tau,
                        Wz_mat, Wzi_mat, Full_basis, Full_loc,
                        rho_star, link, area_per_pixels)
        if (runif(1, 0, 1) < lr) {
          tau      <- tau_new
          theta_kl <- theta_temp
          Ccc      <- Ccc + 1
        }
      }
    }
    
    write.table(tau, file = Ta, sep = ',', append = TRUE,
                quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    # Progress line: iteration number and running acceptance rates
    cat(paste("\rIteration", i,
              "acc", round(Acc / i, 2),   # pCN acceptance rate
              round(Bcc / i, 2),   # alpha acceptance rate
              round(Ccc / i, 2)))  # tau   acceptance rate
    
    end      <- proc.time()[3]
    times[i] <- end - start
  }
  
  assign("times", times, envir = .GlobalEnv)
  print(sum(times))
  
  # Close streaming file connections
  close(rstar)
  close(intens)
  close(alf)
  close(Ta)
}