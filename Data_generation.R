###############################################################################################
###### SIMULATIONS IN Increasing domain asymptotics for covariate-based nonparametric #########
###### Bayesian intensity estimation with Gaussian and Besov-Laplace priors ###################
###############################################################################################

#################################### LOAD LIBRARIES #######################################

library(spatstat)
library(MASS)
library(ggplot2)
library(plotly)
library(zeallot)
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
library(sn)
################### DEFINITION OF OBSERVATION WINDOW AND COVARIATE PROCESS GENERATION ##############

# covariance function
exponential_cov <- function(d, length_scale = 1, sigma_2 = 1) {
  sigma_2 * exp(-d / (2 * length_scale))
}

################################## UNIVARIATE INTENSITY FUNCTIONS #############################

# Choose among next functions 
# negative exponential
{
  rho <- function(z){
    100/2.34 * exp(3*(1-z)-1)
  }
  integral(rho, 0 ,1)
}
# skew normal bump 
{
  rho <- function(z){
    100 * dsn(z, xi = 0.8, omega = 0.3, alpha = -5)
  }
  integral(rho, 0 ,1)
}
# negative/positive deviation from plateau
{ 
  smootherstep <- function(x) {
    x <- pmin(pmax(x, 0), 1)
    6*x^5 - 15*x^4 + 10*x^3
  }
  
  # periodic distance on [0,1] from center c
  periodic_dist <- function(z, c) {
    d <- abs(((z - c + 0.5) %% 1) - 0.5)
    d
  }
  
  # plateau builder: returns value in [0,1]; ~=1 near center, ~=0 outside
  plateau <- function(z, center = 0.5, width = 0.1) {
    # width = full plateau width (value ~1 for |d| <= width/2)
    d <- periodic_dist(z, center)
    t <- d / (width/2)       # normalized distance:
    1 - smootherstep(t)      # = 1 when d=0, ->0 as d >= width/2
  }
  
  # final function: baseline 5, min plateau at 1/4 => 0, max plateau at 3/4 => 15
  rho <- function(z, width = 3/8) {
    # ensure z in [0,1] (vectorized)
    z <- z %% 1
    baseline <- 100
    # bump amplitudes:
    amp_max <- baseline   # +10
    amp_min <- baseline   # 5 (we will subtract this plateau to reach 0)
    baseline + amp_max * plateau(z, center = 3/4, width = width) -
      amp_min * plateau(z, center = 1/4, width = width)
  }
  integral(rho, 0 ,1)
}
# V shaped 1.75 sobolev regular function (SAVVA 113)
{ 
  theta <- function(l){
    l^(-2.25) *sin(10*l)
  }
  orth_basis <- function(l,x){
    sqrt(2) * sin(pi * outer(l, x, `*`))
  }
  rho <- function(x){
    100/0.43 * (1 + 1.1*(theta(seq(1:200)) %*% orth_basis(seq(1:200), x))[1,])
  }
  integral(rho, 0 ,1)
  rho_star = 250
}
# tau shaped 1 sobolev regular function (SAVVA 117)
{ 
  theta <- function(l){
    l^(-3/2) *sin(l)
  }
  orth_basis <- function(l,x){
    sqrt(2) * cos(pi * outer((l-0.5), x, `*`))
  }
  rho <- function(x){
    0.01 + 100/0.672 * (theta(seq(1:300)) %*% orth_basis(seq(1:300), x))[1,]
  }
  integral(rho, 0, 1)
}
# Block function Savva 116
{
  K_j <- function(t){
    (1 + sign(t))/2
  }
  t_j = c(0.1,0.15,0.25,0.40,0.71,0.81) #c(0.1,0.13,0.15,0.23,0.25,0.40,0.44,0.65,0.76,0.78,0.81)
  h_j = c(3,-4,3.1,-2.2,3.1,-3) # c(4,-5,3,-4,5,-4.2,2.1,4.3,-3.1,2.1,-4.2)
  rho <- function(x){
    100/1.644 * (1.01 + h_j %*% K_j(mapply(function(x){x - t_j}, x)))[1,]
  }
  integral(rho, 0, 1)
}
# Spike function Savva 116
{
  K_j <- function(t){
    (1 + abs(t))^(-4)
  }
  t_j = c(0.1,0.13,0.15,0.23,0.25,0.40,0.44,0.65,0.76,0.78,0.81)
  h_j = c(4,5,3,4,5,4.2,2.1,4.3,3.1,5.1,4.2)
  w_j = c(0.5,0.5,0.6,1,1,3,1,1,0.5,0.8,0.5)/100
  rho <- function(x){
    100/0.28 * (h_j %*% K_j(mapply(function(x){(x - t_j)/w_j}, x)))[1,]
  }
  integral(rho, 0, 1)
}
# HeaviSine
{
  rho <- function(x){
    100/5.16*(6 + 4*sin(4*pi*x) - sign(x - 0.3) - sign(0.72 - x))
  }
  integral(rho, 0, 1)
}
# Doppler
{
  epsilon = 0.05
  rho <- function(x){
    100/0.548 * (0.5 + (x*(1-x))^(1/2) * sin(2 * pi * (1 + epsilon)/(x + epsilon)))
  }
  integral(rho, 0, 1)
}

# plot of true intensity 
plot(seq(0, 1, length.out = 512), 
     rho(seq(0, 1, length.out = 512)), 
     type ="l", main = "True intensity function rho", ylab = "intensity", xlab = "covariate values")


#################################### BIVARIATE INTENSITY FUNCTIONS #############################

# Choose among next 3 functions 

# isotropic double bump function
{
  rho1 <- function(z1, z2){
    100 * dmsn(cbind(z1, z2), xi=c(0.4, 0.6), Omega = diag(0.05, nrow = 2), alpha = c(3,-2), tau=1, dp=NULL, log=FALSE)
  }
  rho2 <- function(z1, z2){
    0 * dmsn(cbind(z1, z2), xi=c(0.3, 0.8), Omega = diag(0.03, nrow = 2), alpha = c(-1,-1), tau=0, dp=NULL, log=FALSE)
  }
  rho2d <- function(z1, z2){
    pmax(0, rho1(z1,z2) + rho2(z1,z2))
  }
  integrate2d(rho2d)
}
# monotonic slope function 
{
  rho1 <- function(z1, z2){
    6 * dmsn(cbind(z1, z2), xi=c(0.3, 0.3), Omega = diag(0.5, nrow = 2), alpha = c(-1, -1), tau=1, dp=NULL, log=FALSE)
  }
  rho2d <- function(z1, z2){
    100/0.465 * pmax(0, 2 - rho1(z1,z2))
  }
  integrate2d(rho2d)
}
# anisotropic minimum and maximum 
{
  rho2d <- function(x, y) {
    baseline <- 100
    
    # bump parameters (positive Gaussian)
    A_b <- baseline
    x_b <- 0.3; y_b <- 0.8
    sx_b <- 0.08; sy_b <- 0.50
    
    # hole parameters (negative Gaussian, deep enough to reach 0)
    A_h <- baseline   # ensures depth reaches zero
    x_h <- 0.8; y_h <- 0.3
    sx_h <- 0.08; sy_h <- 0.50   # parallel anisotropy (same orientation as bump)
    
    # helper: anisotropic Gaussian
    gauss2d <- function(x, y, A, x0, y0, sx, sy) {
      A * exp(-((x - x0)^2 / (2 * sx^2) + (y - y0)^2 / (2 * sy^2)))
    }
    
    bump <- gauss2d(x, y, A_b, x_b, y_b, sx_b, sy_b)
    hole <- gauss2d(x, y, A_h, x_h, y_h, sx_h, sy_h)
    
    f_raw <- baseline + bump - hole
    
    # clip at 0 (non-smooth floor)
    pmax(f_raw, 0)
  }
  integrate2d(rho2d)
}
# 2D spike 
{
  K2d_sp <- function(u, v, p=4) {
    (1 + sqrt(u^2 + v^2))^(-p)
  }
  
  # spike centers (x_j, y_j)
  x_j <- c(0.7)#c(0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.2, 0.6)
  y_j <- c(0.8)#c(0.2, 0.4, 0.8, 0.6, 0.8, 0.1, 0.85, 0.35)
  h_j <- c(15)#c(4, 5, 3, 6, 5, 8, 4, 6)        # heights
  w_j <- c(0.1)#c(0.05, 0.08, 0.1, 0.07, 0.11, 0.09, 0.2, 0.1)  # widths
  
  # define the function f(x,y)
  rho2d <- Vectorize(function(x, y) {
    100/0.38 * sum(sapply(1:length(x_j), function(j) {
      h_j[j] * K2d_sp((x - x_j[j]) / w_j[j], (y - y_j[j]) / w_j[j])
    }))
  })
  
  # Example: evaluate on a grid
  x <- seq(0,1,length=100)
  y <- seq(0,1,length=100)
  z <- outer(x,y, rho2d)
  # Plot
  image(x,y,z, col=terrain.colors(100))
  contour(x,y,z, add=TRUE)
  
  integrate2d(rho2d)
}
# 2D block 
{
  K2d <- function(t1, t2, s1, s2){
    (1 + sign(t1))/2 * (1 + sign(t2))/2 * (1 - sign(s1))/2 * (1 -sign(s2))/2
  }
  t1_j <- c(0.1)#c(0.1, 0.6)
  t2_j <- c(0.2)#c(0.2, 0.5)
  s1_j <- c(0.3)#c(0.3, 0.8)
  s2_j <- c(0.5)#c(0.5, 0.8)
  a_j  <- c(4)#c(4, -3)
  
  rho2d <- Vectorize(function(x, y) {
    100/3.06 * sum(sapply(1:length(h_j), function(j) {
      a_j[j] * (3.01 + K2d(x - t1_j[j], y - t2_j[j], x - s1_j[j], y - s2_j[j]))
    }))
  })
  integrate2d(rho2d)
}
#spike + block
{
  rho2d <- Vectorize(function(x, y) {
    100/0.42 * sum(sapply(1:length(x_j), function(j) {
      h_j[j] * K2d_sp((x - x_j[j]) / w_j[j], (y - y_j[j]) / w_j[j]) + 
        a_j[j] * (K2d(x - t1_j[j], y - t2_j[j], x - s1_j[j], y - s2_j[j]))
    }))
  })
}

# plot of true intensity
exp_cov <- expand.grid(x = seq(0,1, length.out = 128), 
                       y = seq(0,1, length.out = 128))
rho_imag <- im(matrix(rho2d(exp_cov$x,exp_cov$y),
                      sqrt(length(exp_cov$x)), sqrt(length(exp_cov$x)), byrow = TRUE),
               xrange = c(min(exp_cov$x),max(exp_cov$x)), yrange = c(min(exp_cov$y),max(exp_cov$y)))
plot(rho_imag, axes=TRUE, main = "Rho function", xlab = "First covariate", 
     ylab = "Second covariate")

####################################### DATA GENERATION ##########################################
covariate_list1 <- list()
covariate_list2 <- list()
loc_list <- list()
i = 0

# number of covariates 
D = 1

for(n in c(1,2,4,8,16)){
  window <- owin(xrange = c(0,n), yrange = c(0,n))
  grid_points <- gridcentres(window, nx = 50, ny = 50)
  distance_matrix <- pairdist(grid_points, squared = TRUE)
  # covariance matrix of the GP generating the covariate process, discretized over the grid
  Cov <- exponential_cov(distance_matrix, length_scale = 0.5, sigma_2 = 1) + 1e-10 * diag(1,length(grid_points$x))
  Chol <- Matrix::chol(Cov)
  dim = length(grid_points$x)
  Cov2 <- exponential_cov(distance_matrix, length_scale = 1.5, sigma_2 = 1) + 1e-10 * diag(1,length(grid_points$x))
  Chol2 <- Matrix::chol(Cov2)
  # mean vector 
  mu = rep(0, length(grid_points$x))
  for(exp in 1:50){
    covariate_process <- as.vector(rnorm(n = dim, mean = 0, sd = 1) %*% Chol)
    covariate_process <- pnorm((covariate_process - mean(covariate_process))/sd(covariate_process))
    covariate_process <- im(matrix(covariate_process, sqrt(dim), sqrt(dim))[sqrt(dim):1,],
                            xrange = window$xrange, yrange = window$yrange)
    
    
    covariate_process2 <- as.vector(rnorm(n = dim, mean = 0, sd = 1) %*% Chol2)
    covariate_process2 <- pnorm((covariate_process2 - mean(covariate_process2))/sd(covariate_process2))
    covariate_process2 <- im(matrix(covariate_process2, sqrt(dim), sqrt(dim))[sqrt(dim):1,],
                             xrange = window$xrange, yrange = window$yrange)
    
    covariate_list1[[50*i + exp]] <- covariate_process 
    covariate_list2[[50*i + exp]] <- covariate_process2
    
    # 
    if(D == 1){
      loc_list[[50*i + exp]] <- rpoispp(lambda = eval.im(rho(covariate_process)), win= window, nsim = 1)
    }else{
      loc_list[[50*i + exp]] <- rpoispp(lambda = im(matrix(rho2d(as.vector(covariate_process$v),
                                as.vector(covariate_process2$v)),
                                50,50),
                                xrange = window$xrange, yrange = window$yrange),
                                win= window, nsim = 1)
    }
    print(c(n, exp))
  }
  i = i+1
}

