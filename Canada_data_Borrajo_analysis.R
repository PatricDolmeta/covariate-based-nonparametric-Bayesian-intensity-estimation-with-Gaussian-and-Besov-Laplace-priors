#Load required packages
library(RandomFields)
library(spatstat)
library(ks)
library(maptools)
library(purrr)
library(mapdata)
spatstat.options(npixel=128)

#Load data set
load("CanadaJune2015.RData")

# Kernel Functions
Kepa<-function(u) { return((0.75*(1-(u)^2))*(abs(u)<1))}
Kepaderiv<-function(u){return(-1.5*u)}
Kepaderivseg<-function(u){return(rep(-1.5,length(u)))}

Kquartic<-function(u) { return((15/16)*((1-u^2)^2)*(abs(u)<1))}
Kquarticderiv<-function(u){return((-15/4)*u*(1-u^2)*(abs(u)<1))}
Kquarticderivseg<-function(u){return((15/4)*(3*u^2-1))*(abs(u)<1)}

Kgauss<-function(u){return((1/sqrt(2*pi))*exp(-0.5*(u^2)))}
Kgaussderiv<-function(u){return((-u/sqrt(2*pi))*exp(-0.5*u^2))}
Kgaussderivseg<-function(u){return((exp(-0.5*u^2)/sqrt(2*pi))*(u^2-1))}

#Simpson's numerical integration
simpson2<-function(fxs,a,b){
    n<-length(fxs)
	h=(b-a)/(n-1)
int<-3*(fxs[1]+fxs[n])/8+7*(fxs[2]+fxs[n-1])/6+23*(fxs[3]+fxs[n-2])/24+sum(fxs[4:(n-3)])
return(int*h)}

#Spatial cumulative distribution function associated to the covariate
gstar<-function(Xi,covariate){
w=Xi$window
z=as.numeric(covariate[w, drop=FALSE])
ghat=density(z)
ghatfun=approxfun(ghat$x,ghat$y)
gstarfun=function(x){ghatfun(x)*area(w)}
return(gstarfun)
}

#First order derivative of the spatial cumulative distribution function
gstard<-function(Xi,covariate){
w=Xi$window
z=as.numeric(covariate[w, drop=FALSE])
h=hns(z,deriv.order=1)
zrej=seq(min(z)-3*h,max(z)+3*h,len=512)
u<-sapply(1:length(z),function(i) (zrej-z[i])/h)
if (K=='epa') {Kud<-Kepaderiv}; if (K=='quartic') {Ku<-Kquarticderiv}; if (K=='gauss') {Ku<-Kgaussderiv}
ghatd=(h^(-2))*rowMeans(Ku(u),na.rm=T)
ghatdfun=approxfun(zrej,ghatd)
gstardfun=function(x){ghatdfun(x)*area(w)}
return(gstardfun)}
  
#Second order derivative of the spatial cumulative distribution function
gstardseg<-function(Xi,covariate){
w=Xi$window
z=as.numeric(covariate[w, drop=FALSE])
h=hns(z,deriv.order=2)
zrej=seq(min(z)-3*h,max(z)+3*h,len=512)
u<-sapply(1:length(z),function(i) (zrej-z[i])/h)
if (K=='epa') {Ku<-Kepaderivseg}; if (K=='quartic') {Ku<-Kquarticderivseg}; if (K=='gauss') {Ku<-Kgaussderivseg}
ghatd2=(h^(-3))*rowMeans(Ku(u),na.rm=T)
ghatd2fun=approxfun(zrej,ghatd2)
gstard2fun=function(x){ghatd2fun(x)*area(w)}
return(gstard2fun)}

#Intensity estimator
lambda.est.data<-function(Xi,covariate,K,h){
w=Xi$window
Zi=interp.im(covariate,Xi$x,Xi$y)
poss=is.na(Zi)
Zi=na.exclude(Zi)
z=as.numeric(covariate[w]);ngrid=length(z)
n=Xi$n
if (ngrid>1) {u<-sapply(1:length(Zi),function(i) (Zi[i]-z)/h); u<-t(u)}else {u<-(Zi-z)/h}
if (K=='epa') {Ku<-Kepa}; if (K=='quartic') {Ku<-Kquartic}; if (K=='gauss') {Ku<-Kgauss}
gstar.z=gstar(Xi,covariate)(z)
gstar.Zi=gstar(Xi,covariate)(Zi)
aux=1/gstar.Zi; Lu<-aux%*%Ku(u)
if (ngrid>1)  rho.z<-h^(-1)* Lu else {rho.z<-h^(-1)* Ku(u) }
rho.z=as.numeric(rho.z)
f.z<-(gstar.z/n)*rho.z
rhof=approxfun(z,rho.z)
lambda=eval.im(rhof(covariate))
return(list('f.est'=as.vector(f.z),'rho.est'=as.vector(rho.z),'lambda.est'=lambda))
}

# Bandwidth selectors
#Rule-of-thumb
#auxiliary function computing gstar and its derivates that only depend on the covariate
aux.RT.data<-function(Xi,covariate){ 
w=Xi$window

z=as.numeric(covariate[w])

zrej=seq(min(z),max(z),len=128)
gstarz=gstar(Xi,covariate)(zrej)
gstardz=gstard(Xi,covariate)(zrej)
gstard2z=gstardseg(Xi,covariate)(zrej)

return(list('gz'=gstarz,'gdz'=gstardz,'gd2z'=gstard2z))
}
h.RT.data=function(Xi,covariate,K,gstarz,gstardz,gstard2z){
  if(K=='epa'){RK=3/5;mu2K2=1/5;Kud=Kepaderiv;Kud2=Kepaderivseg}
  if(K=='quartic'){RK=5/7; mu2K2=1/7;Kud=Kquarticderiv;Kud2=Kquarticderivseg}
  if(K=='gauss'){RK<-sqrt(2*pi);mu2K2=1;Kud=Kgaussderiv;Kud2=Kgaussderivseg}
  n=Xi$n
  Zi=interp.im(covariate, Xi$x, Xi$y)
  Zi=na.exclude(Zi)
  w=Xi$window
  z=as.numeric(covariate[w])
  zrej=seq(min(z),max(z),len=128)  
  Am=1/n
  m=n
  mu.est=mean(Zi); sigma.est=sd(Zi)
  ff=dnorm(zrej,mu.est,sigma.est)
  fnd=function(z,mu,sigma){Kgaussderiv((z-mu)/sigma)}
  fndz=fnd(zrej,mu.est,sigma.est)
  fnd2=function(z,mu,sigma){Kgaussderivseg((z-mu)/sigma)}
  fnd2z=fnd2(zrej,mu.est,sigma.est)
  aux.f=fnd2z-2*(fndz*gstardz/gstarz)-(ff*gstard2z/gstarz)+2*(ff*(gstardz^2)/(gstarz^2))
  Rseg=simpson2(aux.f^2,zrej[1],zrej[length(zrej)])
  h=(Am*RK/((1-exp(-m))^2*Rseg*mu2K2))^(1/5)
  return(h.rt=h)
}

#Bootstrap
h.boot.data<-function(Xi,covariate, L, h.rt){
  if(L=='epa'){RL=3/5;mu2L2=1/5; Lseg=Kepaderivseg}
  if(L=='quartic'){RL=5/7; mu2L2=1/7; Lseg=Kquarticderivseg}
  if(L=='gauss'){RL<-sqrt(2*pi);mu2L2=1;Lseg=Kgaussderivseg}
m=Xi$n
g=((m^(-1/7))/(m^(-1/5)))*h.rt #pilot bandwidth rescaled from the rule-of-thumb
Zi<-interp.im(covariate, Xi$x, Xi$y); Zi=na.exclude(Zi)
n=length(Zi)
w=Xi$window
z=as.numeric(covariate[w])
hatrhog=lambda.est.data(Xi,covariate,L,g)$rho.est
gstar.z=gstar(Xi,covariate)(z)
func=hatrhog*gstar.z
auxfunc=cbind(z,func)
auxfunc=auxfunc[order(auxfunc[,1]),]
hatm=sum((auxfunc[-1,1]-auxfunc[-length(auxfunc[,1]),1])*auxfunc[-length(auxfunc[,2]),2])
ahatm=1/hatm
ngrid=length(z)
if (ngrid>1) {u<-sapply(1:n,function(i) (Zi[i]-z)/g); u<-t(u)}else {u<-(Zi-z)/g}
gstar.Zi=gstar(Xi,covariate)(Zi)
aux=1/gstar.Zi; Lu<-aux%*%Lseg(u)
if (ngrid>1)  rhoseg.z<-g^(-3)* Lu else {rhoseg.z<-g^(-3)* Lseg(u) }
rhoseg.z=as.numeric(rhoseg.z)
ff=(rhoseg.z*gstar.z/hatm)^2
auxff=cbind(z,ff)
auxff=auxff[order(auxff[,1]),]
Rrhoseg=sum((auxff[-1,1]-auxff[-length(auxff[,1]),1])*auxff[-length(auxff[,2]),2])
h=((ahatm*RL)/(mu2L2*((1-exp(-hatm))^2)*Rrhoseg))^(1/5)
return(h)
}

#Non-model-based bandwidth with covariates
bw.NM<-function(Xi,covariate,K,ns){
W <- Window(Xi)
areaW <- area.owin(W)
nnd <- nndist(Xi)
srange <- c(min(nnd[nnd > 0]), diameter(W)/2)
sigma<-seq(from=srange[1],to=srange[2],length=ns)
cv <- numeric(ns)
for (i in 1:ns) {
        si <- sigma[i]
        lam <- lambda.est.data(Xi,covariate,K,h=si)
		Zi=interp.im(covariate,Xi$x,Xi$y)
		z=as.numeric(covariate[W])
		rhofun<-approxfun(z,lam$rho.est)
		lamx<-rhofun(Zi)
		lamx<-lamx[!is.na(lamx)]
		cv[i] <- (sum(1/lamx) - areaW)^2
    }
iopt<-which.min(cv)
h<-sigma[iopt]
return(h)
}



#Application to Canada data
map("worldHires","Canada", xlim=c(-141,-53), ylim=c(40,85), col="gray90", fill=TRUE, main="Wildfires June 2015")
points(Xijun,pch=3,cex=0.5)

#Plot of the covariates 
windows()
par(mfrow = c(2,1), 
    mar = c(0, 1, 0, 2),  # Bottom, left, top, right margins
    oma = c(0, 0, 0, 0))

plot(temp_can, main = "")
plot(prec_can, main = "")


plot(cov.temp,box=FALSE,main="")
points(Xijun$x, Xijun$y, 
       cex = 0.2, col = "white", pch = "+")

plot(cov.precip,box=FALSE,main="")
points(Xijun$x, Xijun$y, 
       cex = 0.2, col = "white", pch = "+")

#Computation of the different methods
Xi=Xijun
covariate=cov.temp
K="gauss"
#RT bandwidth
aux=aux.RT.data(Xi,covariate)
gstarz=aux$gz; gstardz=aux$gdz; gstard2z=aux$gd2z
h_rt=h.RT.data(Xi,covariate,'gauss',gstarz,gstardz,gstard2z)
int_RT=lambda.est.data(Xi,covariate,'gauss',h_rt)$lambda
#Bootstrap bandwidth
hboot=h.boot.data(Xi,covariate, K, h_rt)
int_Boot=lambda.est.data(Xi,covariate,'gauss',hboot)$lambda
#non-model-based
bw_NM<-bw.NM(Xijun,cov.temp,K,16)
int_NM<-lambda.est.data(Xijun,cov.temp,'gauss',bw_NM)$lambda
#Baddeley (2012)
int_badd<-predict(rhohat(Xijun,cov.temp,method="reweight"))
#Diggle (1985)
int_dig<-density(Xijun)
#Cronie and van Lieshout(2018)
b_CvL<-bw.CvL(Xijun)
int_CvL<-density(Xijun,b_CvL)

#Representation of the resulting intensities
#limits for common plots
liminf<-min(range(int_dig)[1],range(int_CvL)[1],range(int_badd)[1],range(int_RT)[1],range(int_Boot)[1],range(int_NM)[1])
limsup<-max(range(int_dig)[2],range(int_CvL)[2],range(int_badd)[2],range(int_RT)[2],range(int_Boot)[2],range(int_NM)[2])
zlim=c(liminf,limsup)
windows()
par(mfrow=c(3,2))
plot(int_dig,zlim=zlim,main="Diggle (1985)")
plot(int_CvL,zlim=zlim,main="Cronie and van Lieshout (2018)")
plot(int_badd,zlim=zlim, main="Baddeley (2012)")
plot(int_RT,zlim=zlim, main="Our rule-of-thumb")
plot(int_Boot,zlim=zlim, main="Our bootstrap bandwidth")
plot(int_NM,zlim=zlim, main="Our non-model-based bandwidth")

#Application with two covariates (temerature and precipitation)
gstar.multi<-function(Xi,covariate1,covariate2){
if(!is.ppp(Xi)){print('Xi must be a point process');stop}
if(!is.im(covariate1)){print('The covariates must be images');stop}
if(!is.null(covariate2)){if(!is.im(covariate2)){print('The covariates must be images');stop}
if((covariate1$xrange[1]!=covariate2$xrange[1])|(covariate1$xrange[2]!=covariate2$xrange[2])|(covariate1$yrange[1]!=covariate2$yrange[1])|(covariate1$yrange[2]!=covariate2$yrange[2])){print('The covariates have different ranges');stop}}

w=Xi$window

#model with two covariates
Zi1=interp.im(covariate1,Xi$x,Xi$y)
poss1=is.na(Zi1)
Zi1=na.exclude(Zi1)
Zi1=as.vector(Zi1)
Zi2=interp.im(covariate2,Xi$x,Xi$y)
poss2=is.na(Zi2)
Zi2=na.exclude(Zi2)
Zi2=as.vector(Zi2)
Zi=cbind(Zi1,Zi2)
x.coord<-as.im((raster.x(as.mask(w), drop=FALSE)),W=w)
z1=as.numeric(x.coord[w, drop=FALSE])
y.coord<-as.im((raster.y(as.mask(w), drop=FALSE)),W=w)
z2=as.numeric(y.coord[w, drop=FALSE])
z1=as.numeric(covariate1[w, drop=FALSE])
z2=as.numeric(covariate2[w, drop=FALSE])
z=cbind(z1,z2)
ghat.z=kde(z,eval.points=z)
ghat.Zi=kde(z,eval.points=Zi)
return(list('ghat.z'=ghat.z$estimate,'ghat.Zi'=ghat.Zi$estimate))
}
lambda.est.multi<-function(Xi,covariate1,covariate2){
if(!is.ppp(Xi)){print('Xi must be a point process');stop}
if(!is.im(covariate1)){print('The covariates must be images');stop}
if(!is.null(covariate2)){if(!is.im(covariate2)){print('The covariates must be images');stop}
if((covariate1$xrange[1]!=covariate2$xrange[1])|(covariate1$xrange[2]!=covariate2$xrange[2])|(covariate1$yrange[1]!=covariate2$yrange[1])|(covariate1$yrange[2]!=covariate2$yrange[2])){print('The covariates have different ranges');stop}}

w=Xi$window
#model with two covariates
		Zi1=interp.im(covariate1,Xi$x,Xi$y)
		poss1=is.na(Zi1)
		Zi1=na.exclude(Zi1)
		Zi2=interp.im(covariate2,Xi$x,Xi$y)
		poss2=is.na(Zi2)
		Zi2=na.exclude(Zi2)
		Zi=cbind(Zi1,Zi2)
		n=dim(Zi)[1]

		z1=as.numeric(covariate1[w, drop=FALSE])
		z2=as.numeric(covariate2[w, drop=FALSE])
		z=cbind(z1,z2)

		gstar=gstar.multi(Xi,covariate1,covariate2)
		ghat.z=gstar$ghat.z
		ghat.Zi=gstar$ghat.Zi

		kk=kde(Zi,eval.points=z,w=(1/ghat.Zi)/(sum(1/ghat.Zi)/n))
		rho.z=kk$estimate*(sum(1/ghat.Zi)/n)
		rho.z=as.numeric(rho.z)
		f.z<-(ghat.z/n)*rho.z

pos.notna=which(!is.na(covariate1$v))
rho.m=matrix(NA,nr=128,nc=128)
rho.m[pos.notna]=rho.z
rho.im=as.im(rho.m,xrange=covariate1$xrange,yrange=covariate1$yrange)
lambda=rho.im

return(list('f.est'=as.vector(f.z),'rho.est'=as.vector(rho.z),'lambda.est'=lambda))}

int_multi<-lambda.est.multi(Xi,cov.temp,cov.precip)$lambda
windows()
plot(int_multi,main="Intensity estimation with TWO covariates")


############################ GP analysis for Borrajo data ###############################

cov1 <- rescaling(cov.temp)

rho_hat <- (rhohat(Xijun, cov.temp))

setwd("~/Research/IPP")
source("pCN_adaptive_MCMC_GP.R")

discr <- seq(0, 1, length.out = 200)
discr_orig <- seq(min(cov.temp), max(cov.temp), length.out = 200)
grid_points <- expand.grid(x = cov1$xcol, y = cov1$yrow)

n_iter = 20000

adapt_MCMC("canada_temp_BORRAJO_2015_ind_adapt", n_iter, beta = 0.1, 
           K = length(discr), Xijun, grid_points, 
           cov.temp, discr_orig, 
           matern52_kernel, alpha1 = 1, alpha2 = 1, sigma_2 = 1, nu = 5/2, 
           exp_param = c (2,2), shape = 1, rate = 1, link = "lambda")

lambda_post = read.table("Post_l_star.csv", sep = ",")
plot(t(lambda_post)[1:n_iter], type = "l")

Ls = read.table("Post_ls.csv", sep = ",")
plot(t(Ls)[1:n_iter], type = "l")

gs = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
gs <- apply(gs, 2, as.numeric)
post_mean_dom <- colMeans(gs[seq(8/10 * dim(gs)[1],dim(gs)[1]),],na.rm = TRUE)

par(mfrow = c(2,1), 
    mar = c(0, 0, 0, 0),  # Bottom, left, top, right margins
    oma = c(0, 0, 0, 0))


plot(rho_hat,
     xlab = expression(Z(x)),
     ylab = expression(rho(z)), main = "",
     legend = FALSE)
lines(discr_orig, post_mean_dom, col = "red", lwd = 2)
upper = apply(gs[seq( 8/10 *dim(gs)[1], dim(gs)[1]),], 2, quantile, probs = 0.975, na.rm = TRUE)
lower = apply(gs[seq( 8/10 *dim(gs)[1], dim(gs)[1]),], 2, quantile, probs = 0.025, na.rm = TRUE)
polygon(c(discr_orig, rev(discr_orig)), c(upper, rev(lower)), col = rgb(1, 0, 0, alpha = 0.2), border = NA)

# intensity on the spatial domain 
intensity_est <- function(rho, mask, Z1){
  # browser()
  Z1 <- Z1[!is.na(Z1)]
  res = matrix(NA, nrow = dim(mask)[1], ncol = dim(mask)[2])
  entries <- approx(seq(min(Z1),max(Z1),length.out = length(rho)), rho, xout = Z1)$y
  res[mask] <- entries
  return(res)
}

plot(as.im(intensity_est(rho = rho_hat$rho , !is.na(cov.temp$v), cov.temp$v), W = window),
     main = "")
points(Xijun$x, Xijun$y, 
       cex = 0.2, col = "white", pch = "+")

plot(as.im(intensity_est(rho = post_mean_dom, !is.na(cov.temp$v), cov.temp$v), W = window),
     main = "")
points(Xijun$x, Xijun$y, 
       cex = 0.2, col = "white", pch = "+")

# incorporating also precipitation 

cov2 = rescaling(cov.precip)

rho_2_hat <- rho2hat(Xijun, rescaling(cov.temp), rescaling(cov.precip))
rho_2_hat <- rho2hat(Xijun, cov1, cov2)

plot(rho_2_hat)

intensity2d_est <- function(rho, mask, Z1, Z2){
  # browser()
  Z1 <- Z1[!is.na(Z1)]
  Z2 <- Z2[!is.na(Z2)]
  Z1 <- Z1[!is.na(Z1)]
  Z2 <- Z2[!is.na(Z2)]
  res = matrix(NA, nrow = dim(mask)[1], ncol = dim(mask)[2])
  entries <- interp2(rho$xcol, rho$yrow, rho$v, Z1, Z2, method = "nearest")
  entries[is.na(entries)] <- 0
  res[mask] <- entries
  return(res)
}
plot(as.im(intensity2d_est(rho = rho_2_hat, !is.na(cov.temp$v),
                           Z1 = as.vector(cov1), Z2 = as.vector(cov2)), 
           W = Xijun$window),
     main = "Posterior mean of the intensity function with observations")

# triangula mesh on unit cube
{
points <- matrix(c(0,0,0,1,1,1,1,0), ncol = 2, byrow = T)
conv_hull <- chull(points)
conv_hull <- rev(c(conv_hull, conv_hull[1]))
poly_win <- owin(poly = list(x = points[,1][conv_hull], y = points[,2][conv_hull]))

plot(poly_win, main = "Covariate domain", xlab = "Slope", ylab = "Elevetion")

# triangular tesseletion of the convex hull, independent of points 

tri_tess <- RTriangle::pslg(P = cbind(poly_win$bdry[[1]]$x, poly_win$bdry[[1]]$y))
mesh <- RTriangle::triangulate(tri_tess, a = 0.003, j = TRUE)  #Y = TRUE, D = TRUE # smaller 'a' = finer mesh
triangle_centroids <- t(apply(mesh$T, 1, function(tri_idx) {
  colMeans(mesh$P[tri_idx, ])
}))

triangle_area <- function(a, b, c) {
  0.5 * abs((b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2]))
}
triangle_area <- apply(mesh$T, 1, function(tri) {
  a <- mesh$P[tri[1], ]
  b <- mesh$P[tri[2], ]
  c <- mesh$P[tri[3], ]
  triangle_area(a, b, c)
})

plot(mesh, asp = 1, main = "Uniform Triangular Tessellation of Convex Hull")
points(triangle_centroids, col = "blue", pch = 16, cex = .3)
vertices <- mesh$P
extremes <- mesh$P[mesh$PB==1,]
colnames(extremes) <- c("x", "y")
colnames(vertices) <- c("x", "y")

inside.owin(x = extremes[,1], y = extremes[,2], w = poly_win)
}

# triangular mesh on original covariates
{
points <- matrix(c(min(cov.temp$v, na.rm = TRUE), min(cov.precip$v, na.rm = TRUE), 
                   max(cov.temp$v, na.rm = TRUE), min(cov.precip$v, na.rm = TRUE),
                   max(cov.temp$v, na.rm = TRUE), max(cov.precip$v, na.rm = TRUE),
                   min(cov.temp$v, na.rm = TRUE), max(cov.precip$v, na.rm = TRUE)), ncol = 2, byrow = T)

conv_hull <- chull(points)
conv_hull <- rev(c(conv_hull, conv_hull[1]))
# hull_points <- cbind((y-min(slope))/(max(slope)-min(slope)), (x-min(elev))/(max(elev)-min(elev)))[conv_hull,]
poly_win <- owin(poly = list(x = points[,1][conv_hull], y = points[,2][conv_hull]))

Range = 1 #diff(range(poly_win$bdry[[1]]$y)) / diff(range(poly_win$bdry[[1]]$x))
N_triag = 400
tri_tess <- RTriangle::pslg(P = cbind(poly_win$bdry[[1]]$x*Range, poly_win$bdry[[1]]$y), PB = rep(1,length(poly_win$bdry[[1]]$x)))
mesh <- RTriangle::triangulate(tri_tess, a = area.owin(poly_win)*Range/N_triag, j = TRUE, D = TRUE)#  q = 45 # smaller 'a' = finer mesh
triangle_centroids <- t(apply(mesh$T, 1, function(tri_idx) {
  colMeans(mesh$P[tri_idx, ])
}))

triangle_area <- function(a, b, c) {
  0.5 * abs((b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2]))
}
triangle_area <- apply(mesh$T, 1, function(tri) {
  a <- mesh$P[tri[1], ]
  b <- mesh$P[tri[2], ]
  c <- mesh$P[tri[3], ]
  triangle_area(a, b, c)
})

plot(mesh, main = "Uniform Triangular Tessellation of Convex Hull",asp =1)

vertices <- mesh$P
extremes <- mesh$P[mesh$PB==1,]
colnames(extremes) <- c("x", "y")
colnames(vertices) <- c("x", "y")

center <- c(mean(extremes[,1]), mean(extremes[,2]))
points(center[1], center[2], pch = 16, col = "green")
# Compute the angle of each point relative to the centroid
angles <- atan2(extremes[,2] - center[2], extremes[,1] - center[1])
# Order indices by the angle
order_indices <- order(angles)
# Optionally, close the polygon by repeating the first point at the end
order_indices <- c(order_indices)
poly_win <- owin(poly = list(x = extremes[,1][order_indices], y = extremes[,2][order_indices]))
poly_win$bdry[[1]]$x <- extremes[,1][order_indices]
poly_win$bdry[[1]]$y <- extremes[,2][order_indices]

inside.owin(x = extremes[,1], y = extremes[,2], w = poly_win)
}

setwd("~/Research/IPP/")
source(file = "pCN_adaptive_MCMC_3mesh_GP.R")
# unlink("~/Research/IPP/2dAdaptive_MCMC_post_canadian_fires_july15", recursive = TRUE)

n_iter = 30000

MCMC("_canadian_fires_borrajo_transf", n_iter, beta = 0.08, K = dim(vertices)[1], 
     Xijun, grid_points, Xijun$window,
     cov1, cov2, 
     as.data.frame(vertices), triangle_centroids,
     matern52_kernel, alpha1=1, alpha2=2, sigma_2 = 1, nu = 5/2,
     exp_param = c(2,2), shape = 1, rate = 2, "lambda", "isotropic")


lambda_post = read.table("Post_l_star.csv", sep = ",")
plot(t(lambda_post)[1:dim(lambda_post)[1]], type = "l")

Ls = read.table("Post_ls.csv", sep = ",")
plot(t(Ls)[2,1:dim(Ls)[1]], type = "l")

gs_m = read.table("Post_gp.csv", sep = ",", fill = TRUE, header=FALSE)
gs_m <- apply(gs_m, 2, as.numeric)
pp <- ppp(x = vertices[,1], y = vertices[,2],
          marks = colMeans(gs_m[seq(8/10 * (dim(gs_m)[1]), dim(gs_m)[1]),],na.rm = TRUE), window = owin(poly_win))
post_mean_m <- Smooth(pp, kernel = "gaussian", dimyx = c(dim(mesh$P)[1], dim(mesh$P)[1]), method = "pwl")
plot(post_mean_m, main = "")
plot(rho_2_hat, xaxt = "n")
axis(side = 1, at = seq(0,1,0.2),labels = round(seq(min(cov.temp$v, na.rm = TRUE), max(cov.temp$v, na.rm = TRUE), length.out = 6)))
axis(side = 2, at = seq(0,1,0.2),labels = round(seq(min(cov.precip$v, na.rm = TRUE), max(cov.precip$v, na.rm = TRUE), length.out = 6)))

pp_uncert <- ppp(x = vertices[,1], y = vertices[,2],
                 marks = apply(gs_m[seq(8/10 * (dim(gs_m)[1]), dim(gs_m)[1]),], 2, quantile, probs = 0.95) - 
                   apply(gs_m[seq(8/10 * (dim(gs_m)[1]), dim(gs_m)[1]),], 2, quantile, probs = 0.05), window = owin(poly_win))

post_uncert <- Smooth(pp_uncert, kernel = "gaussian", dimyx = c(dim(mesh$P)[1], dim(mesh$P)[1]), method = "pwl")
plot(post_uncert, 
     main = "")
points(as.vector(cov1$v), as.vector(cov2$v), cex = 0.2, col = "white", pch = "+")


plot(as.im(intensity2d_est(rho = post_mean_m, !is.na(cov.temp$v),
                           Z1 = as.vector(rescaling(cov.temp)$v), Z2 = as.vector(rescaling(cov.precip)$v)), 
           W = Xijun$window),
     main = "")

points(Xijun$x, Xijun$y, 
       cex = 0.2, col = "white", pch = "+")

plot(as.im(intensity2d_est(rho = rho_2_hat, !is.na(cov.temp$v),
                           Z1 = as.vector(rescaling(cov.temp)$v), Z2 = as.vector(rescaling(cov.precip)$v)), 
           W = Xijun$window),
     main = "")

points(Xijun$x, Xijun$y, 
       cex = 0.2, col = "white", pch = "+")
