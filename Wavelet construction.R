################################# WAVELET BASIS CONTRUCTION ##############################

# 1. Generate a non-stationary auxiliary signal
set.seed(123)
# discretization choice 
length = 2^8
t <- seq(0, 1, length.out = length)
# true intensity funtion rho has to be defined 
signal <- rho(t)

# 2. Perform DWT with Daubechies wavelet (e.g. d8, haar)
K <- log2(length)
# reflective or periodic boundaries 
refl = FALSE
wt <- dwt(signal, wf="d8", n.levels=K, boundary = ifelse(refl, "reflection", "periodic"))

# 3. Full set of Daubechies basis functions

basis_list <- list()

coarse <- paste0("s", K, sep = "")

# coarse basis functions
for(k in 1:length(wt[[coarse]])){
  tmp <- wt
  # zero all detail coefficients
  for(i in 1:K) tmp[[paste0("d",i)]] <- rep(0, length(tmp[[paste0("d",i)]]))
  # zero s6 then set one coefficient
  tmp[[coarse]] <- rep(0, length(tmp[[coarse]]))
  tmp[[coarse]][k] <- 1
  # reconstruct
  basis_list[[paste0("s",0, "_",k)]] <- idwt(tmp)
}

# detail basis functions
for(lvl in K:1){
  dname <- paste0("d", lvl)
  for(k in 1:length(wt[[dname]])){
    tmp <- wt
    # zero all coefficients
    for(i in 1:K) tmp[[paste0("d",i)]] <- rep(0, length(tmp[[paste0("d",i)]]))
    tmp[[coarse]] <- rep(0, length(tmp[[coarse]]))
    
    # set one coefficient
    tmp[[dname]][k] <- 1
    
    # reconstruct
    basis_list[[paste0("d",K-lvl,"_",k)]] <- idwt(tmp)
  }
}

name_vec <- names(basis_list)
levels <- as.integer(sub("^[A-Za-z](\\d+).*", "\\1", name_vec))

plot(t, basis_list$d1_1, type="l", main="Detail basis function level 6, index 60")

# Manual reconstruction using basis functions 
# truncation level of choice L <= K 
L = 6
manual_rec <- rep(0, length)

# Add coarse
for(k in 1:length(wt[[coarse]])){
  manual_rec <- manual_rec + wt[[coarse]][k] * basis_list[[paste0("s",0, "_",k)]]
}

# Add details
for(lvl in 0:(L-1)){
  dname <- paste0("d",lvl)
  for(k in 1:length(wt[[paste0("d",K-lvl)]])){
    manual_rec <- manual_rec + wt[[paste0("d",K-lvl)]][k] * basis_list[[paste0(dname,"_",k)]]
  }
}

#################### FREQUENTIST KERNEL ESTIMATE #####################
plot(t, rho(t), col = "black", type = "l")
# optimal reconstruction with selected truncation level 
lines(t, manual_rec, col="purple", lty=4, lwd = 1)

###########################################################################################
########################### 2D WAVELET BASIS CONSTRUCTION #################################

#discretization grid size 
K = 4
exp_cov <- expand.grid(x = seq(0,1, length.out = 2^K), 
                       y = seq(0,1, length.out = 2^K))
# rho2d must be available 
rho_imag <- im(matrix(rho2d(exp_cov$x,exp_cov$y),
                      sqrt(length(exp_cov$x)), sqrt(length(exp_cov$x)), byrow = TRUE),
               xrange = c(min(exp_cov$x),max(exp_cov$x)), yrange = c(min(exp_cov$y),max(exp_cov$y)))


library(wavethresh)
refl=FALSE
#family = "DaubExPhase" and "DaubLeAsymm".
# filter.number = 1 is Haar, higher stands for differentiable basis

wd2d_obj <- imwd(rho_imag$v, filter.number=1, family="DaubExPhase", type="wavelet",
                 bc= ifelse(refl, "symmetric", "periodic"), RetFather=FALSE, verbose=FALSE)

recon_img <- imwr.imwd(wd2d_obj)

image(rho_imag, col=gray(0:255/255), main="Original")
image(im(recon_img), col=gray(0:255/255), main="Reconstructed")

d2_basis_list <- list()

coarse <- "w0Lconstant"

for(k in 1:length(wd2d_obj[[coarse]])){
  tmp <- wd2d_obj
  # zero all detail coefficients
  for(i in 0:(K-1)){
    D <- length(tmp[[paste0("w", i, "L1")]])
    tmp[[paste0("w", i, "L1")]] <- rep(0, D)
    tmp[[paste0("w", i, "L2")]] <- rep(0, D)
    tmp[[paste0("w", i, "L3")]] <- rep(0, D)
  } 
  # zero coarse then set one coefficient
  tmp[[coarse]] <- rep(0, D)
  tmp[[coarse]][k] <- 1
  # reconstruct
  d2_basis_list[[paste(coarse,k, sep = "_")]] <-  imwr.imwd(tmp)
}

for(lvl in 0:(K-1)){
  dname1 <- paste0("w", lvl, "L1")
  dname2 <- paste0("w", lvl, "L2")
  dname3 <- paste0("w", lvl, "L3")
  for(k in 1:length(wd2d_obj[[dname1]])){
    tmp <- wd2d_obj
    # zero all coefficients
    for(i in 0:(K-1)){
      D <- length(tmp[[paste0("w", i, "L3")]])
      tmp[[paste0("w", i, "L1")]] <- rep(0, D)
      tmp[[paste0("w", i, "L2")]] <- rep(0, D)
      tmp[[paste0("w", i, "L3")]] <- rep(0, D)
    } 
    
    tmp[[coarse]] <- rep(0, D)
    # browser()
    # set one coefficient
    tmp1 <- tmp
    tmp2 <- tmp
    tmp3 <- tmp
    
    tmp1[[dname1]][k] <- 1
    tmp2[[dname2]][k] <- 1
    tmp3[[dname3]][k] <- 1
    
    # reconstruct
    d2_basis_list[[paste0(dname1,"_",k)]] <- imwr.imwd(tmp1)
    d2_basis_list[[paste0(dname2,"_",k)]] <- imwr.imwd(tmp2)
    d2_basis_list[[paste0(dname3,"_",k)]] <- imwr.imwd(tmp3)
  }
}

#8. Manual reconstruction using basis functions 
# truncation level L <= K 
L = 3
manual_rec <- matrix(0, 2^K, 2^K)

# Add coarse
for(k in 1:length(wd2d_obj[[coarse]])){
  manual_rec <- manual_rec + wd2d_obj[[coarse]][k] * d2_basis_list[[paste0(coarse, "_",k)]]
}

# Add details
for(lvl in 0:(L-1)){
  # browser()
  dname1 <- paste0("w", lvl, "L1")
  dname2 <- paste0("w", lvl, "L2")
  dname3 <- paste0("w", lvl, "L3")
  for(k in 1:length(wd2d_obj[[dname1]])){
    manual_rec <- manual_rec + wd2d_obj[[dname1]][k] * d2_basis_list[[paste0(dname1,"_",k)]]
    manual_rec <- manual_rec + wd2d_obj[[dname2]][k] * d2_basis_list[[paste0(dname2,"_",k)]]
    manual_rec <- manual_rec + wd2d_obj[[dname3]][k] * d2_basis_list[[paste0(dname3,"_",k)]]
  }
}

# plot of single basis 
# w0Lconstant is the Coarse (Scaling) Level
# Detail Levels: w{lvl}L{1,2,3}
# For each decomposition level lvl ∈ {0, 1, ..., K-1}, there are three detail subbands, corresponding to the three orientations produced by a 2D separable wavelet transform:
# L1 Horizontal details
# L2 Vertical details
# L3 Diagonal details 
# Spatial index: k
# k is the spatial index of a single coefficient within that subband. 
# Example
# w0L1_3 3rd horizontal-detail basis element at coarsest detail level

image(im(d2_basis_list$w0L1_1,
         xrange = c(min(exp_cov$x),max(exp_cov$x)), yrange = c(min(exp_cov$y),max(exp_cov$y)))
      , col=gray(0:255/255), main="Single wavelet basis")

# best possible reconstructionat given truncation level 
image(im(manual_rec,
         xrange = c(min(exp_cov$x),max(exp_cov$x)), yrange = c(min(exp_cov$y),max(exp_cov$y)))
      , col=gray(0:255/255), main="Manual Reconstructed")
