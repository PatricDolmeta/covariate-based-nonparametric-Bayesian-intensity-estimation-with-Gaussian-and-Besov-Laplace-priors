####################################### LOAD LIBRARIES #############################################

library(spatstat)
library(MASS)
library(ggplot2)
library(plotly)
library(zeallot)
library(fourPNO)
library(fMultivar)
library(akima)
library(sn)
library(Matrix)
library(profvis)
library(kernlab)
library(mvtnorm)
library(xtable)
library(pracma)
library(patchwork)
library(viridis)
library(scico)
library(RColorBrewer)
library(pals)
library(tictoc)
library(waveslim)
library(wavethresh)
set.seed(123)

# load wavelet basis construction 
# set K discretization level 
# set L truncation level 
##################################### REAL DATA BEI ######################3###############
data(bei)

slope <- bei.extra$grad
elev <- bei.extra$elev

plot(slope, axes = TRUE)
plot(elev)
contour(slope)

# define color palettes
my_col_map <- colorRampPalette(c("darkblue", "deeppink3", "darkgoldenrod1"))
first_palette <- colorRampPalette(coolwarm(11))
spectral_colors <- brewer.pal(11, "Spectral")
second_palette <- colorRampPalette(rev(spectral_colors))

r <- range(elev)
mid <- mean(r)

# comment out line 68 and "[ontario_dom, drop = FALSE], add = TRUE," at line 69 to remove canada's shadow
par(mfrow = c(1,2),        
    mar = c(0, 2, 0.5, 2),  # Bottom, left, top, right margins
    oma = c(0, 0, 0, 0))  
plot(elev, main = "",
     col = first_palette(256), ribbon = TRUE, ribargs = list(cex.axis = 1.5, cex.lab = 2, at = c(120,140,159)),
     axes = TRUE, ribsep = 0.05 ) #
points(bei, pch = 18, col = "black", cex = 0.3)
plot(slope, main = "", col = second_palette(256), ribargs = list(cex.axis = 1.5, cex.lab = 2), ribsep = 0.05, axes = TRUE)
points(bei, pch = 18, col = "black", cex = 0.3)


# GP-IPP
window <- bei$window
# [0,1] transform through CFD
covariate_p <- bei.extra$grad
covariate_p2 <- bei.extra$elev
# to unit interval with empirical cdf
F_z1 <- ecdf(covariate_p$v)
F_z2 <- ecdf(covariate_p2$v)
covariate_process1 <- eval.im(F_z1(covariate_p))
covariate_process2 <- eval.im(F_z2(covariate_p2))
# [0,1] scaling correcting for range 
grid_points <- expand.grid(x = covariate_process1$xcol, y = covariate_process1$yrow)
discr <- seq(0, 1, length.out = 2**K)

######################################### 1D ANALYSIS #########################################
setwd("~/Research/BezovLaplacePriors/CODE")
source(file = "MCMC.R")
folder = paste0("~/Research/BezovLaplacePriors/CODE")
dir.create(folder)
setwd(folder)
n_iter <- 200
exp = "elev_adapt_d8"
setwd("~/Research/BezovLaplacePriors/CODE")
MCMC(exp, n_iter, Laplace = FALSE,
     bei, list(covariate_process2),
     as.matrix(discr), window, 
     beta = 0.04, tau = 1, alpha = 0.05,
     shape = NA, rate = NA, "exponential", rho_star,
     basis_list, L, refl,
     adaptA = TRUE, b_a = 0.2,
     adaptT = TRUE, b_t = 0.1)

post_alpha = read.table("Post_alpha.csv", sep = ",", fill = TRUE, header=FALSE)
plot(post_alpha[,1], type = "l")

tau_post = read.table("Post_tau.csv", sep = ",", fill = TRUE, header=FALSE)
plot(tau_post[,1], type = "l")

plot(rhohat(bei, covariate_process2, method="ratio", n = 256), 
     main = "Intensity estimates")

# upload posterior draws for rho
post_coeff = read.table("Post_intens.csv", sep = ",", fill = TRUE, header=FALSE)
post_coeff <- apply(post_coeff, 2, as.numeric)

name_vec <- names(basis_list)
levels <- as.integer(sub("^[A-Za-z](\\d+).*", "\\1", name_vec))
L_tot <- sum(levels<L)
d <- 1
basis_mat <- do.call(cbind, basis_list[1:L_tot])
res_mat <- exp(post_coeff[,1:L_tot] %*% t(basis_mat))

intensity_g <- colMeans(res_mat[seq(5/10*n_iter,n_iter),], na.rm = TRUE)

lines(discr, intensity_g, col = "blue")


upper_g = apply(res_mat[seq(5/10*n_iter,n_iter),], 2, quantile, probs = 0.95, na.rm = TRUE)
lower_g = apply(res_mat[seq(5/10*n_iter,n_iter),], 2, quantile, probs = 0.05, na.rm = TRUE)
polygon(c(discr, rev(discr)), c(upper, rev(lower)), col = rgb(0, 0, 1, alpha = 0.2), border = NA)

intensity_est <- function(rho, mask, Z1){
  # browser()
  Z1 <- Z1[!is.na(Z1)]
  Z1 <- Z1[!is.na(Z1)]
  res = matrix(NA, nrow = mask$dim[1], ncol = mask$dim[2])
  entries <- approx(seq(0,1,length.out = length(rho)), rho, xout = Z1)$y
  res[mask$v] <- entries
  return(res)
}
mask <- im(!is.na(covariate_process1$v), xcol = covariate_process1$xcol, yrow = covariate_process1$yrow)
plot(im(intensity_est(rho = rhohat(bei, covariate_process1, method="ratio", n = 256)$rho, mask = mask, covariate_process1),
        xrange = covariate_process1$xrange, yrange = covariate_process1$yrange),
     main = "Posterior mean of the intensity function with observations")
plot(bei, add = TRUE, cex = 0.5, col = "green", pch = "+")

############################################## GGplot ######################################

all_pixels <- as.vector(bei.extra$elev$v)
xx_grid = quantile(all_pixels, probs = discr, type = 1, na.rm = TRUE)

rug <- bei.extra$elev[bei]

kern <- rhohat(bei,covariate_process2, n = length)

# Figure 5: change indexing of plot_1d list when changing covariate 
# plot_1d <- list()
plot_1d[[1]] <- ggplot(data.frame(x=xx_grid, y=intensity_g), aes(x = x, y = y)) + 
  geom_line(linetype = "solid", color = "dodgerblue3", size = 1) +
  theme_minimal(base_family = "sans") +  # Clean white background
  theme(
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    text = element_text(color = "black"),
    axis.text  = element_text(size = 11)
  ) + labs(x = "Elevation", y = expression(rho(z))) + 
  geom_line(data = data.frame(x=xx_grid, y=kern$rho), mapping = aes(x = x, y = y),
            color = "firebrick", size = 1, linetype = "solid") +
  geom_ribbon(data = data.frame(x = xx_grid, y2 = upper, 
                                y1 = lower, 
                                y= intensity_g), 
              mapping = aes(x = x, ymin = y1, y = y, ymax = y2), fill = "dodgerblue3", alpha = 0.2) + 
  geom_ribbon(data = data.frame(x = xx_grid, y2 = kern$hi, 
                                y1 = kern$lo, 
                                y= kern$rho), 
              mapping = aes(x = x, ymin = y1, y = y, ymax = y2), fill = "firebrick", alpha = 0.2) +
  geom_rug(data = data.frame(x = rug, y = 0), sides = "b") + 
  ylim(c(0,0.015)) # + 
  geom_line(data = data.frame(x=xx_grid, y=intensity_g), mapping = aes(x = x, y = y),
            color = "seagreen4", size = 1) +
  geom_ribbon(data = data.frame(x = xx_grid, y2 = upper_g, 
                                y1 = lower_g, 
                                y= intensity_g), 
              mapping = aes(x = x, ymin = y1, y = y, ymax = y2), fill = "seagreen4", alpha = 0.2)

# Figure 6: change indexing of plot_lambda_1d list when changing year, 

plot_lambda_1d <- list()
# intensity on the spatial domain 
intensity_est <- function(rho, mask, Z1){
  # browser()
  Z1 <- Z1[!is.na(Z1)]
  Z1 <- Z1[!is.na(Z1)]
  res = matrix(NA, nrow = mask$dim[1], ncol = mask$dim[2])
  entries <- approx(seq(0,1,length.out = length(rho)), rho, xout = Z1)$y
  res[mask$v] <- entries
  return(res)
}

img_df <- data.frame(value = as.vector(t((im(intensity_est(rho = intensity, mask = mask, covariate_process1),
                                             xrange = covariate_process1$xrange, yrange = covariate_process1$yrange)$v)[dim(mask$v)[1]:1,])))
img_df$x <- grid_points$x
img_df$y <- rev(grid_points$y)
points <- data.frame(x = bei$x, y =bei$y)

plot_lambda_1d[[2]] <- 
  ggplot(img_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +  
  scale_fill_gradientn(na.value = "white", colors = 
                         c("darkblue", "deeppink3", "darkgoldenrod1"),  
                       name = expression(lambda(x)),
                       limits = c(0,0.014)) +
  coord_fixed() + theme_minimal() + labs(x = expression(x[1]), 
                                         y = expression(x[2])) # +
geom_point(data = points, aes(x = x, y = y), shape = 19, color="turquoise1", size = 0.5, inherit.aes = FALSE) 

combined_plot <- wrap_plots(plot_1d, nrow = 1)
combined_plot

###################################### 2D-ANALYSIS ################################################

kern2_est <- rho2hat(bei, covariate_process1, covariate_process2, method = "ratio", dimyx = c(128,128), from = 0, to =1)
plot(kern2_est)

points(covariate_process1[bei], covariate_process2[bei])
x <- as.vector(covariate_process1$v)
y <- as.vector(covariate_process2$v)
points <- cbind(x,y)[complete.cases(cbind(x,y)),]
conv_hull <- chull(points)
conv_hull <- rev(c(conv_hull, conv_hull[1]))
poly_win <- owin(poly = list(x = points[,1][conv_hull], y = points[,2][conv_hull]))
lines(poly_win$bdry[[1]]$x, poly_win$bdry[[1]]$y, col = "white", cex = 2, pch = 16)

setwd("~/Research/BezovLaplacePriors/CODE")
source(file= "MCMC.R")
folder = paste0("~/Research/BezovLaplacePriors/CODE")
dir.create(folder)
setwd(folder)
n_iter <- 150
discr2d <- expand.grid(x = seq(0,1,length.out = 2**K), 
                       y = seq(0,1,length.out = 2**K))
exp = "2dbei_symm"
setwd("~/Research/BezovLaplacePriors/CODE")
MCMC(exp, n_iter, Laplace = FALSE,
     bei, list(covariate_process1, covariate_process2),
     as.matrix(discr2d), window, 
     beta = 0.04, tau = 1, alpha = 0.5,
     shape = NA, rate = NA, "exponential", rho_star,
     d2_basis_list, L, refl,
     adaptA = TRUE, b_a = 0.15,
     adaptT = TRUE, b_t = 0.25)

post_alpha = read.table("Post_alpha.csv", sep = ",", fill = TRUE, header=FALSE)
plot(post_alpha[c(1:n_iter),1], type = "l")

tau_post = read.table("Post_tau.csv", sep = ",", fill = TRUE, header=FALSE)
plot(tau_post[c(1:n_iter),1], type = "l")

name_vec <- names(d2_basis_list)
levels <- as.integer(sub("^[A-Za-z](\\d+).*", "\\1", name_vec))
L_tot <- sum(levels<L)
post_coeff = read.table("Post_intens.csv", sep = ",", fill = TRUE, header=FALSE)
post_coeff <- apply(post_coeff, 2, as.numeric)

d <- 2
basis_mat <- sapply(d2_basis_list[1:L_tot], function(M) as.vector(M))  
res_mat <- exp(post_coeff[,1:L_tot] %*% t(basis_mat))
intensity <- colMeans(res_mat[seq(5/10*n_iter,n_iter),])

intensity_im_4 <- as.im(matrix(intensity, nrow = sqrt(dim(discr2d)[1]), byrow = TRUE),
                      xrange = c(min(exp_cov$x),max(exp_cov$x)), yrange = c(min(exp_cov$y),max(exp_cov$y)))

# smooth_im <- Smooth(intensity_im,sigma = 0.01)

image(intensity_im_4, main="Posterior mean", zlim = c(0, 0.03))

points(covariate_process1[bei], covariate_process2[bei], add = TRUE, cex = 0.5, col = "white", pch = "+")

intensity2d_est <- function(rho, covariate_p, Z1, Z2){
  Z1 = pmax(pmin(Z1, max(rho$xcol)), min(rho$xcol))
  Z2 = pmax(pmin(Z2, max(rho$yrow)), min(rho$yrow))
  Z1 = Z1[!is.na(Z1)]
  Z2 = Z2[!is.na(Z2)]
  res = matrix(NA, nrow = covariate_p$dim[1], ncol = covariate_p$dim[2])
  entries <- interp2(rho$xcol, rho$yrow, rho$v, Z1, Z2, method = "nearest")
  res[!is.na(covariate_p$v)] <- entries
  return(res)
}
plot(as.im(intensity2d_est(rho = intensity_im, covariate_process1,
                           Z1 = as.vector(covariate_process1$v), Z2 = as.vector(covariate_process2$v)), 
           W = window),
     main = "Posterior mean of the intensity function with observations")
points(bei$x, bei$y, cex = 0.2, col = "white", pch = "+")

########################################### GGplot ##########################################
gg_2D_plot <- list()

all_pixelsX <- as.vector(bei.extra$grad$v)
all_pixelsY <- as.vector(bei.extra$elev$v)

img_df <- data.frame(value = (as.vector(intensity_im_4$v))) 
exp_cov <- expand.grid(x = seq(0, 1, length.out = 2**K), 
                       y = seq(0, 1, length.out = 2**K))

img_df$x <- exp_cov$y
img_df$y <- exp_cov$x

points <- data.frame(x = covariate_process1[bei],
                     y = covariate_process2[bei])


img_df$value2 <- as.vector(intensity_im$v) 
img_df$value3 <- as.vector(kern2_est$v) 

gg_2D_plot[[1]] <- ggplot() +
  geom_raster(data = img_df, aes(x = x, y = y, fill = value)) +  # or geom_tile() if you prefer
  # scale_fill_viridis(option = "magma", name = expression(lambda(x)), limits = c(0,40)) +
  # scale_fill_gradientn(colors = colors, name = expression(lambda(x)), limits = c(0, 120)) +
  scale_fill_gradientn(na.value = "white", colors = c("darkblue", "deeppink3", "darkgoldenrod1"),  
                       name = expression(rho(z)), limits = c(0, 0.037)) +
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() + 
  labs(x = "Gradient", y = "Elevation") + 
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsX, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE),2)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsY, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE), 2)) 


###################################### REAL DATA CANADA ###################################

# load data from Canada_data_Borrajo_analysis.R

plot(cov.precip)
plot(cov.temp)

# comment out line 68 and "[ontario_dom, drop = FALSE], add = TRUE," at line 69 to remove canada's shadow
par(mfrow = c(1,2),        
    mar = c(0, 1, 0, 2),  # Bottom, left, top, right margins
    oma = c(0, 0, 0, 0))  
plot(cov.temp, main = "",
     col = first_palette(256), ribbon = TRUE, ribargs = list(cex.axis = 1.5, cex.lab = 2, at = c(8,15,23)),
     axes = TRUE, ribsep = 0.05 ) #
points(Xijun, pch = 18, col = "black", cex = 0.4)
plot(cov.precip, main = "", col = brewer.blues(256), ribargs = list(cex.axis = 1.5, cex.lab = 2, at = c(0.5,3,6)),
     axes = TRUE, ribsep = 0.05 ) #
points(Xijun, pch = 18, col = "black", cex = 0.4)


# GP-IPP
window <- owin(Xijun$window$xrange, Xijun$window$yrange)
# [0,1] transform through CFD
covariate_p <- cov.precip
covariate_p2 <- cov.temp
# to unit interval with empirical cdf
F_z1 <- ecdf(covariate_p$v)
F_z2 <- ecdf(covariate_p2$v)
covariate_process1 <- eval.im(F_z1(covariate_p))
covariate_process2 <- eval.im(F_z2(covariate_p2))
plot(density(covariate_process1$v, na.rm = TRUE))
# [0,1] scaling correcting for range 
grid_points <- expand.grid(x = covariate_process1$xcol, y = covariate_process1$yrow)
discr <- seq(0, 1, length.out = length)
discr2d <- expand.grid(x = seq(0,1,length.out = 2**K), 
                       y = seq(0,1,length.out = 2**K))

setwd("~/Research/BezovLaplacePriors/CODE")
source(file = "MCMC")
folder = paste0("~/Research/BezovLaplacePriors/CODE")
dir.create(folder)
setwd(folder)
n_iter <- 100
exp = "temp_adapt_8la"
setwd("~/Research/BezovLaplacePriors/CODE/Canada")
MCMC(exp, n_iter, Laplace = FALSE,
     Xijun, list(covariate_process2),
     as.matrix(discr), window, 
     beta = 0.04, tau = 1, alpha = 0.05,
     shape = NA, rate = NA, "exponential", rho_star,
     basis_list, L, refl,
     adaptA = TRUE, b_a = 0.15,
     adaptT = TRUE, b_t = 0.25)

post_alpha = read.table("Post_alpha.csv", sep = ",", fill = TRUE, header=FALSE)
plot(post_alpha[,1], type = "l")

tau_post = read.table("Post_tau.csv", sep = ",", fill = TRUE, header=FALSE)
plot(tau_post[,1], type = "l")

plot(rhohat(Xijun, covariate_process2, method="ratio", n = 256), 
     main = "Intensity estimates", legend = none)

# upload posterior draws for rho
post_coeff = read.table("Post_intens.csv", sep = ",", fill = TRUE, header=FALSE)
post_coeff <- apply(post_coeff, 2, as.numeric)
# post_coeff <- apply(post_coeff[seq(5/10*n_iter,n_iter),], 2, function(x){dens <- density(x); dens$x[which.max(dens$y)]})
# post_coeff <- colMeans(post_coeff)
name_vec <- names(basis_list)
levels <- as.integer(sub("^[A-Za-z](\\d+).*", "\\1", name_vec))
L_tot <- sum(levels<L)
d <- 1
basis_mat <- do.call(cbind, basis_list[1:L_tot])
res_mat <- exp(post_coeff[,1:L_tot] %*% t(basis_mat))

intensity <- colMeans(res_mat[seq(5/10*n_iter,n_iter),], na.rm = TRUE)
# intensity <- res_mat

lines(discr, intensity, col = "blue")


upper = apply(res_mat[seq(3/10*n_iter,n_iter),], 2, quantile, probs = 0.95, na.rm = TRUE)
lower = apply(res_mat[seq(3/10*n_iter,n_iter),], 2, quantile, probs = 0.05, na.rm = TRUE)
polygon(c(discr, rev(discr)), c(upper, rev(lower)), col = rgb(0, 0, 1, alpha = 0.2), border = NA)

intensity_est <- function(rho, mask, Z1){
  # browser()
  Z1 <- Z1[!is.na(Z1)]
  Z1 <- Z1[!is.na(Z1)]
  res = matrix(NA, nrow = mask$dim[1], ncol = mask$dim[2])
  entries <- approx(seq(0,1,length.out = length(rho)), rho, xout = Z1)$y
  res[mask$v] <- entries
  return(res)
}
mask <- im(!is.na(covariate_process1$v), xcol = covariate_process1$xcol, yrow = covariate_process1$yrow)
plot(im(intensity_est(rho = intensity, mask = mask, covariate_process1),
        xrange = covariate_process1$xrange, yrange = covariate_process1$yrange),
     main = "Posterior mean of the intensity function with observations",
     zlim = c(0,11))
plot(Xijun, add = TRUE, cex = 0.5, col = "green", pch = "+")

############################# GGplot #####################

all_pixels <- as.vector(cov.precip$v)
xx_grid = quantile(all_pixels, probs = discr, type = 1, na.rm = TRUE)

rug <- cov.precip[Xijun]

kern <- rhohat(Xijun,covariate_process1, n = 256)

# Figure 5: change indexing of plot_1d list when changing covariate 
# plot_1d <- list()
ggplot(data.frame(x=xx_grid, y=intensity), aes(x = x, y = y)) + 
  geom_line(linetype = "solid", color = "dodgerblue3", size = 1) +
  theme_minimal(base_family = "sans") +  # Clean white background
  theme(
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    text = element_text(color = "black"),
    axis.text  = element_text(size = 11)
  ) + labs(x = "Precipitation", y = expression(rho(z))) + 
  geom_line(data = data.frame(x=xx_grid, y=kern$rho), mapping = aes(x = x, y = y),
            color = "firebrick", size = 1) +
  geom_ribbon(data = data.frame(x = xx_grid, y2 = upper, 
                                y1 = lower, 
                                y= intensity), 
              mapping = aes(x = x, ymin = y1, y = y, ymax = y2), fill = "dodgerblue3", alpha = 0.2) + 
  geom_ribbon(data = data.frame(x = xx_grid, y2 = kern$hi, 
                                y1 = kern$lo, 
                                y= kern$rho), 
              mapping = aes(x = x, ymin = y1, y = y, ymax = y2), fill = "firebrick", alpha = 0.2) +
  geom_rug(data = data.frame(x = rug, y = 0), sides = "b")
  # geom_line(data = data.frame(x=xx_grid, y=intensity_g), mapping = aes(x = x, y = y),
  #           color = "seagreen4", size = 1) +
  # geom_ribbon(data = data.frame(x = xx_grid, y2 = upper_g, 
  #                               y1 = lower_g, 
  #                               y= intensity_g), 
  #             mapping = aes(x = x, ymin = y1, y = y, ymax = y2), fill = "seagreen4", alpha = 0.2) +
 
# Figure 6: change indexing of plot_lambda_1d list when changing year, 

plot_lambda_1d <- list()
# intensity on the spatial domain 
intensity_est <- function(rho, mask, Z1){
  # browser()
  Z1 <- Z1[!is.na(Z1)]
  Z1 <- Z1[!is.na(Z1)]
  res = matrix(NA, nrow = mask$dim[1], ncol = mask$dim[2])
  entries <- approx(seq(0,1,length.out = length(rho)), rho, xout = Z1)$y
  res[mask$v] <- entries
  return(res)
}

img_df <- data.frame(value = as.vector(t((im(intensity_est(rho = intensity, mask = mask, covariate_process1),
                                             xrange = covariate_process1$xrange, yrange = covariate_process1$yrange)$v)[dim(mask$v)[1]:1,])))
img_df$x <- grid_points$x
img_df$y <- rev(grid_points$y)
points <- data.frame(x = Xijun$x, y =Xijun$y)

plot_lambda_1d[[2]] <- 
  ggplot(img_df, aes(x = x, y = y, fill = value)) +
  geom_raster() +  
  scale_fill_gradientn(na.value = "white", colors = 
                         c("darkblue", "deeppink3", "darkgoldenrod1"),  
                       name = expression(lambda(x))) +
  coord_fixed() + theme_minimal() + labs(x = "Longitude", y = "Latitude")# +
geom_point(data = points, aes(x = x, y = y), shape = 19, color="turquoise1", size = 0.5, inherit.aes = FALSE) 

combined_plot <- wrap_plots(plot_1d, nrow = 1)
combined_plot

# 2d

kern2_est <- rho2hat(Xijun, covariate_process1, covariate_process2, 
                     method = "ratio", dimyx = c(2**K,2**K), from = 0, to =1)
plot(kern2_est)

points(covariate_process1[Xijun], covariate_process2[Xijun], cex = 0.5, pch = 19, col = "turquoise1")
x <- as.vector(covariate_process1$v)
y <- as.vector(covariate_process2$v)
points <- cbind(x,y)[complete.cases(cbind(x,y)),]
conv_hull <- chull(points)
conv_hull <- rev(c(conv_hull, conv_hull[1]))
poly_win <- owin(poly = list(x = points[,1][conv_hull], y = points[,2][conv_hull]))
lines(poly_win$bdry[[1]]$x, poly_win$bdry[[1]]$y, col = "white", cex = 2, pch = 16)

setwd("~/Research/BezovLaplacePriors/CODE")
source(file = "MCMC.R")
folder = paste0("~/Research/BezovLaplacePriors/CODE/Xijun")
dir.create(folder)
setwd(folder)
n_iter <- 200
# for(exp in 1:1){
exp = "2dXijun_symm"
setwd("~/Research/BezovLaplacePriors/CODE")
MCMC(exp, n_iter, Laplace = FALSE,
     Xijun, list(covariate_process1, covariate_process2),
     as.matrix(discr2d), window, 
     beta = 0.05, tau = 1, alpha = 0.5,
     shape = NA, rate = NA, "exponential", rho_star,
     d2_basis_list, L, refl = FALSE,
     adaptA = TRUE, b_a = 0.2,
     adaptT = TRUE, b_t = 0.1)

post_coeff = read.table("Post_intens.csv", sep = ",", fill = TRUE, header=FALSE)
post_coeff <- apply(post_coeff, 2, as.numeric)
post_coeff <- apply(post_coeff[seq(5/10*n_iter,n_iter),], 2, function(x){dens <- density(x); dens$x[which.max(dens$y)]})


L = 3

name_vec <- names(d2_basis_list)
levels <- as.integer(sub("w(\\d+).*", "\\1", name_vec))
sum(levels<L)

d <- 2
basis_mat <- sapply(d2_basis_list[1:sum(levels<L)], function(M) as.vector(M))  
res_mat <- exp(post_coeff[1:sum(levels<L)] %*% t(basis_mat))
# intensity <- colMeans(res_mat[seq(5/10*n_iter,n_iter),])
intensity <- res_mat
post_mean <- im(matrix(intensity, nrow = 2^K, byrow = TRUE),
                xrange = c(min(exp_cov$x),max(exp_cov$x)), yrange = c(min(exp_cov$y),max(exp_cov$y)))


smooth_im <- Smooth(post_mean,sigma = 0.03)
image(smooth_im, main="Posterior mean")
lines(poly_win$bdry[[1]]$x, poly_win$bdry[[1]]$y, col = "white", cex = 2, pch = 16)

intensity2d_est <- function(rho, covariate_p, Z1, Z2){
  Z1 = pmax(pmin(Z1, max(rho$xcol)), min(rho$xcol))
  Z2 = pmax(pmin(Z2, max(rho$yrow)), min(rho$yrow))
  Z1 = Z1[!is.na(Z1)]
  Z2 = Z2[!is.na(Z2)]
  res = matrix(NA, nrow = covariate_p$dim[1], ncol = covariate_p$dim[2])
  entries <- interp2(rho$xcol, rho$yrow, rho$v, Z1, Z2, method = "nearest")
  res[!is.na(covariate_p$v)] <- entries
  return(res)
}
plot(as.im(intensity2d_est(rho = kern2_est, covariate_process1,
                           Z1 = as.vector(covariate_process1$v), Z2 = as.vector(covariate_process2$v)), 
           W = window), zlim = c(0,7),
     main = "Posterior mean of the intensity function with observations")
points(Xijun$x, Xijun$y, cex = 1, col = "white", pch = "+")

############################################ GGplot ######################################
all_pixelsX <- as.vector(cov.precip$v)
all_pixelsY <- as.vector(cov.temp$v)

img_df <- data.frame(value = (as.vector(post_mean$v))) 
exp_cov <- expand.grid(x = seq(0, 1, length.out = 2^K), 
                       y = seq(0, 1, length.out = 2^K))

img_df$x <- exp_cov$y
img_df$y <- exp_cov$x

points <- data.frame(x = covariate_process1[Xijun],
                     y = covariate_process2[Xijun])

img_df$value2 <- as.vector(kern2_est$v) 
img_df$value <- as.vector(smooth_im$v) 

ggplot() +
  geom_raster(data = img_df, aes(x = x, y = y, fill = value)) +  # or geom_tile() if you prefer
  # scale_fill_viridis(option = "magma", name = expression(lambda(x)), limits = c(0,40)) +
  # scale_fill_gradientn(colors = colors, name = expression(lambda(x)), limits = c(0, 120)) +
  scale_fill_gradientn(na.value = "white", colors =
                         c("darkblue", "mediumblue", "deeppink3", "deeppink3","darkgoldenrod1"),
                       name = expression(rho(z))) +
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() +
  theme(text = element_text(size = 14),            # Base font size
        axis.title = element_text(size = 18),      # Axis titles
        axis.text  = element_text(size = 14)) +
  labs(x = "Precipitation", y = "Temperature") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = round(quantile(all_pixelsX, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE),2)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = round(quantile(all_pixelsY, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE), 2))# +
  geom_path(data = data.frame(x= poly_win$bdry[[1]]$x, y = poly_win$bdry[[1]]$y), aes(x = x, y = y, fill = NULL),
            col = "moccasin", size = 1)

ggplot() +
  geom_raster(data = img_df, aes(x = x, y = y, fill = value2)) +  # or geom_tile() if you prefer
  # scale_fill_viridis(option = "magma", name = expression(lambda(x)), limits = c(0,40)) +
  # scale_fill_gradientn(colors = colors, name = expression(lambda(x)), limits = c(0, 120)) +
  scale_fill_gradientn(na.value = "white", colors = 
                         c("darkblue", "mediumblue", "deeppink3", "deeppink3","darkgoldenrod1"),  
                       name = expression(rho(z)), limits = c(0,16)) +
  coord_fixed() +  # Keeps aspect ratio correct
  theme_minimal() + 
  theme(text = element_text(size = 14),            # Base font size
        axis.title = element_text(size = 18),      # Axis titles
        axis.text  = element_text(size = 14)) +
  labs(x = "Precipitation", y = "Temperature") + 
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsX, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE),2)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = round(quantile(all_pixelsY, probs = c(c(0, 0.25, 0.5, 0.75, 1)), type = 1, na.rm = TRUE), 2))# +
  geom_path(data = data.frame(x= poly_win$bdry[[1]]$x, y = poly_win$bdry[[1]]$y), aes(x = x, y = y, fill = NULL),
            col = "moccasin", size = 1)

  
img_df <- data.frame(value1 = as.vector(t((im(intensity2d_est(rho = post_mean, covariate_process1,
                                                             Z1 = as.vector(covariate_process1$v), 
                                                             Z2 = as.vector(covariate_process2$v)))$v)[dim(mask$v)[1]:1,])))
img_df$x <- grid_points$x
img_df$y <- rev(grid_points$y)
points <- data.frame(x = Xijun$x, y =Xijun$y)

img_df$value2 <- as.vector(t((int_multi$v)[dim(mask$v)[1]:1,]))
img_df$value3 <- as.vector(t((im(intensity2d_est(rho = kern2_est, covariate_process1,
                          Z1 = as.vector(covariate_process1$v), 
                          Z2 = as.vector(covariate_process2$v)))$v)[dim(mask$v)[1]:1,]))

ggplot(img_df, aes(x = x, y = y, fill = value3)) +
  geom_raster() + 
  scale_fill_gradientn(na.value = "white", colors = 
                         c("darkblue", "mediumblue", "deeppink3", "deeppink3","darkgoldenrod1"),  
                       name = expression(lambda(x)),limits = c(0, 17)) +
  coord_fixed() + theme_minimal() + labs(x = "Longitude", y = "Latitude") #  +
  geom_point(data = points, aes(x = x, y = y), shape = 19, color="turquoise1", size = 0.2, inherit.aes = FALSE) 

img_df$mask <- as.vector(t((mask$v)[dim(mask$v)[1]:1,]))
  
ggplot() +
  geom_raster(data = img_df, aes(x = x, y = y, fill = mask))+
  geom_raster() + scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "gray80")) +
  coord_fixed() + theme_minimal() + labs(x = "Longitude", y = "Latitude")  +
  geom_point(data = points, aes(x = x, y = y), shape = 19, color="black", size = 0.3, inherit.aes = FALSE) +
  theme(legend.position = "none")
    

                                                                  