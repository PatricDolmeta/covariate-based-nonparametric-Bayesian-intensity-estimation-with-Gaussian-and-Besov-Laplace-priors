#################################### BAYESIAN ANALYSIS #########################################

setwd("~/Research/BezovLaplacePriors/CODE")
source(file = "MCMC.R")
folder = paste0("~/")
dir.create(folder)
setwd(folder)

output_folder <- "~/Research/BezovLaplacePriors/RESULTS"
dir.create(output_folder, showWarnings = FALSE)

discr     <- seq(0, 1, length.out = length)   # length.out = length (if previously defined )
N_exp     <- 1
n_iter    <- 2000
burnin    <- floor(n_iter / 2)          

L2_lap_exp     <- matrix(NA, 5, N_exp)
L2_lap_rel_exp <- matrix(NA, 5, N_exp)
AA             <- numeric(N_exp)           

win_sizes <- c(1, 2, 4, 8, 16)          

for (row_idx in seq_along(win_sizes)) {
  win_size <- win_sizes[row_idx]
  window   <- owin(xrange = c(0, win_size), yrange = c(0, win_size))
  offset   <- N_exp * floor(sqrt(win_size)) 
  
  for (exp in 1:N_exp) {
    setwd(output_folder)
    
    MCMC("skn_adap", n_iter, Laplace = FALSE,
         loc_list[[offset + exp]],
         list(covariate_list1[[offset + exp]]),
         as.matrix(discr), window,
         beta = 0.04, tau = 20, alpha = 1,
         shape = NA, rate = NA, "exponential", 250,
         basis_list, L, refl,
         adaptA = TRUE, b_a = 0.05,
         adaptT = TRUE, b_t = 0.05)
    
    post_alpha <- read.table("Post_alpha.csv",
                             sep = ",", fill = TRUE, header = FALSE)
    AA[exp] <- mean(post_alpha[burnin:n_iter, 1])
    
    post_coeff <- read.table("Post_intens.csv",
                             sep = ",", fill = TRUE, header = FALSE)
    post_coeff <- apply(post_coeff, 2, as.numeric)
    
    name_vec  <- names(basis_list)
    levels    <- as.integer(sub("^[A-Za-z](\\d+).*", "\\1", name_vec))
    n_basis   <- sum(levels < L)        
    basis_mat <- do.call(cbind, basis_list[1:n_basis])
    
    intensity <- exp(post_coeff %*% t(basis_mat))
    mean_int <- colMeans(intensity)
    upper = apply(intensity, 2, quantile, probs = 0.95, na.rm = TRUE)
    lower = apply(intensity, 2, quantile, probs = 0.05, na.rm = TRUE)
    
    plot(t, rho(t), col = "black", type = "l")
    lines(t, exp(manual_rec), col="purple", lty=4, lwd = 1)
    lines(discr, mean_int, col = "blue", lwd=1)
    polygon(c(discr, rev(discr)), c(upper, rev(lower)), col = rgb(0, 1, 0, alpha = 0.2), border = NA)
    
    L2_lap_exp[row_idx, exp]     <- sqrt(sum((intensity - rho(discr))^2) / length(discr))
    L2_lap_rel_exp[row_idx, exp] <- sqrt(sum((intensity - rho(discr))^2) / sum(rho(discr)^2))
  }
}

#################### GGPLOT 1D POSTERIOR #######################################

plot_list <-list()

which(lapply(loc_list, npoints)==101)
rug <- covariate_list1[[2]][loc_list[[2]]]
kernel <- rhohat(loc_list[[2]], covariate_list1[[2]], method="ratio", n = length)

plot_list[[1]] <- 
  ggplot(data.frame(x=discr, y=rho(discr)), aes(x = x, y = y)) +
  geom_line(linetype = "solid", color = "black", size = 1) +
  theme_minimal(base_family = "sans") +  # Clean white background
  theme(
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    text = element_text(color = "black")
  ) + labs(x = expression(z), y = expression(rho(z))) + 
  geom_line(data = data.frame(x=discr, y=kernel$rho), mapping = aes(x = x, y = y), 
            color = "firebrick", size = 1) + ylim(0,1.5 * max(rho(discr))) + 
  geom_rug(data = data.frame(x = rug, y = 0), sides = "b") +
  geom_line(data = data.frame(x=discr, y=mean_int), mapping = aes(x = x, y = y),
            color = "dodgerblue3", size = 1) +
  geom_ribbon(data = data.frame(x = discr, y2 = upper, y1 = lower, y=mean_int), mapping = aes(x = x, ymin = y1, y = mean_int, ymax = y2),
              fill = "dodgerblue3", alpha = 0.2) #+
  # geom_line(data = data.frame(x=t, y=mean_int_g), mapping = aes(x = x, y = y),
  #           color = "seagreen4", size = 1) #+
# geom_ribbon(data = data.frame(x = t, y2 = upper_g, y1 = lower_g, y=mean_int_g), mapping = aes(x = x, ymin = y1, y = intensity_g, ymax = y2),
#             fill = "seagreen4", alpha = 0.2)

combined_plot <- wrap_plots(plot_list, nrow = 1)
combined_plot


################################### BAYESIAN ANALYSIS — 2D #########################################
setwd("~/Research/BezovLaplacePriors/CODE")
source(file = "MCMC.R") 
output_folder <- "~/Research/BezovLaplacePriors/CODE"
dir.create(output_folder, showWarnings = FALSE)

K      <- 4  # define explicitly — if not inherited from global state
N_exp  <- 10
n_iter <- 1000
burnin <- floor(n_iter / 2)

discr <- expand.grid(x = seq(0, 1, length.out = 2^K),
                     y = seq(0, 1, length.out = 2^K))

L2_2dbay_exp     <- matrix(NA, 5, N_exp)
L2_2dbay_rel_exp <- matrix(NA, 5, N_exp)
AA               <- numeric(N_exp)

win_sizes <- c(1, 2, 4, 8, 16)

for (row_idx in seq_along(win_sizes)) {
  win_size <- win_sizes[row_idx]
  window   <- owin(xrange = c(0, win_size), yrange = c(0, win_size))
  offset   <- N_exp * floor(sqrt(win_size))  # matches list layout convention from 1D
  
  for (exper in 1:N_exp) {
    setwd(output_folder)
    
    MCMC("bs_adapt", n_iter, Laplace = TRUE,
         loc_list[[offset + exper]],
         list(covariate_list1[[offset + exper]],
              covariate_list2[[offset + exper]]),
         as.matrix(discr), window,
         beta = 0.02, tau = 20, alpha = 0.5,
         shape = NA, rate = NA, "exponential", 250,
         d2_basis_list, L, refl,
         adaptA = TRUE, b_a = 0.03,
         adaptT = TRUE, b_t = 0.05)
    
    # --- alpha trace ---
    post_alpha <- read.table("Post_alpha.csv",
                             sep = ",", fill = TRUE, header = FALSE)
    if (interactive()) plot(post_alpha[1:n_iter, 1], type = "l")
    AA[exper] <- mean(post_alpha[burnin:n_iter, 1])
    
    # --- tau trace ---
    tau_post <- read.table("Post_tau.csv",
                           sep = ",", fill = TRUE, header = FALSE)
    if (interactive()) plot(tau_post[, 1], type = "l")
    
    # --- basis matrix (2D) ---
    name_vec  <- names(d2_basis_list)
    levels    <- as.integer(sub("^[A-Za-z](\\d+).*", "\\1", name_vec))
    n_basis   <- sum(levels < L)
    basis_mat <- sapply(d2_basis_list[1:n_basis], function(M) as.vector(M))
    
    # --- posterior coefficients ---
    post_coeff <- read.table("Post_intens.csv",
                             sep = ",", fill = TRUE, header = FALSE)
    post_coeff <- apply(post_coeff, 2, as.numeric)
    post_coeff <- colMeans(post_coeff[burnin:n_iter, ])
    
    # --- intensity surface ---
    intensity    <- exp(post_coeff %*% t(basis_mat))
    intensity_im <- as.im(matrix(intensity, nrow = 2^K, byrow = TRUE),
                          xrange = c(min(exp_cov$x), max(exp_cov$x)),
                          yrange = c(min(exp_cov$y), max(exp_cov$y)))
    
    if (interactive()) {
      image(intensity_im, main = loc_list[[offset + exper]]$n)
      smooth_im <- Smooth(intensity_im, sigma = 0.01)
      image(smooth_im, main = "Posterior mean")
    }
    
    # --- L2 errors ---
    rho2d_vals <- rho2d(exp_cov$x, exp_cov$y)
    n_grid     <- length(discr[, 1])
    
    L2_2dbay_exp[row_idx, exper]     <- sqrt(sum((intensity - rho2d_vals)^2) / n_grid)
    L2_2dbay_rel_exp[row_idx, exper] <- sqrt(sum((intensity - rho2d_vals)^2) / sum(rho2d_vals^2))
  }
}

#################### Figure 3: functional estimates ################################ 
plot_list2 <- list()

pp <- ppp(x = discr[,1], y = discr[,2],
            marks = intensity[1,], window = window)
post_mean_m <- Smooth(pp, kernel = "gaussian", dimyx = c(200, 200), method = "pwl")
  
exp_cov <- expand.grid(x = seq(1.084202e-19, 1, length.out = 200), 
                       y = seq(1.084202e-19, 1, length.out = 200))
  
img_df <- data.frame(value1 = as.vector(post_mean_m$v))
img_df$x <- exp_cov$y
img_df$y <- exp_cov$x
  
plot_list2[[1]] <-
  ggplot() +
  geom_raster(data = img_df, aes(x = x, y = y, fill = value1)) +
  scale_fill_gradientn(colors = c("darkblue", "mediumblue", "deeppink3", "deeppink3", "darkgoldenrod1"),
                       name = expression(rho(z))) + # , limits = c(0, 55)
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() +
  labs(x = expression(z[1]), y = expression(z[2])) + theme(legend.position = "none")
  
# add ground truth 
img_df <- data.frame(value = as.vector(rho_imag$v))
img_df$x <- gridcentres(window, nx = rho_imag$dim[1], ny = rho_imag$dim[1])$y
img_df$y <- gridcentres(window, nx = rho_imag$dim[1], ny = rho_imag$dim[1])$x

plot_list2[[4]] <- ggplot(img_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +  
  scale_fill_gradientn(colors = c("darkblue","mediumblue", "deeppink3", "deeppink3", "darkgoldenrod1"),  
                       name = expression(rho(z))) + #, limits = c(0, 55)
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() + 
  labs(x = expression(z[1]), y = expression(z[2]))


combined_plot <- wrap_plots(plot_list2, nrow = 1)
combined_plot
