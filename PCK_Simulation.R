#--- Load libraries
library(geoR)
library(mvtnorm)
library(sp)
library(RColorBrewer)
library(gridExtra)
library(gstat)
library(viridis)
library(proxy)
library(matrixcalc)
library(latticeExtra)
library(sf)
seed = 123


# ---------------------------------- 1 -----------------------------------------
# ----------- Generating Poisson Multivariate Spatial Data ---------------------
# ------------------------------------------------------------------------------

# ----- Bivariate Linear Model of Coregionalization

# --- Define region
size = 20*20
x = seq(-10,10,length.out = sqrt(size))
y = seq(-10,10,length.out = sqrt(size))
coordenadas = expand.grid(x = x, y = y)
dist.matrix = as.matrix(dist(coordenadas))

# --- Define covariance structures
Cexp = function(h,a,b){b * exp(-h/a)}
Csph = function(h,a,b){ifelse(h <= a, b * (1-1.5*(h/a)+0.5*(h/a)^3), 0)}

# --- Define coefficient matrices (LMC)
B1 = matrix(c(0.81, 0.36, 0.36, 0.16),nrow=2,byrow=T)
B2 = matrix(c(0.0081, 0.054, 0.054 , 0.360),nrow=2,byrow=T)
is.positive.semi.definite(B1)
is.positive.semi.definite(B2)

# --- Define processes covariance matrix
cov11 = B1[1,1] * Cexp(dist.matrix,2,1) + B2[1,1] * Csph(dist.matrix,5,1)
cov22 = B1[2,2] * Cexp(dist.matrix,2,1) + B2[2,2] * Csph(dist.matrix,5,1)
cov12 = B1[1,2] * Cexp(dist.matrix,2,1) + B2[1,2] * Csph(dist.matrix,5,1)
cov21 = B1[2,1] * Cexp(dist.matrix,2,1) + B2[2,1] * Csph(dist.matrix,5,1)
cov.total =rbind(cbind(cov11,cov12),cbind(cov21,cov22))
is.positive.definite(cov.total)

# Plot LMC
dum = seq(1, 10, length.out = 1000)
par(mfrow = c(2,2))
plot(dum, (B1[1,1] * Cexp(0,2,1) + B2[1,1] * Csph(0,5,1)) -(B1[1,1] * Cexp(dum,2,1) + B2[1,1] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ Y[alpha]/n[alpha] - Y[alpha]/n[alpha]))
plot(dum, (B1[1,2] * Cexp(0,2,1) + B2[1,2] * Csph(0,5,1)) -(B1[1,2] * Cexp(dum,2,1) + B2[1,2] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ Y[alpha]/n[alpha] - Y[beta]/n[beta]))
plot(dum, (B1[2,1] * Cexp(0,2,1) + B2[2,1] * Csph(0,5,1)) -(B1[2,1] * Cexp(dum,2,1) + B2[2,1] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ Y[beta]/n[beta] - Y[alpha]/n[alpha]))
plot(dum, (B1[2,2] * Cexp(0,2,1) + B2[2,2] * Csph(0,5,1)) -(B1[2,2] * Cexp(dum,2,1) + B2[2,2] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance',main = expression(~ Y[beta]/n[beta] - Y[beta]/n[beta]))

# --- Generate multivariate log-Gaussian field
set.seed(seed)
sim = rmvnorm(1,mean=rep(0,2*nrow(coordenadas)), sigma=cov.total)
datos = cbind(coordenadas,Ra = sim[1:size], Rb= sim[(size+1):(2*size)])
cor(datos$Ra, datos$Rb)

# --- Determine shared risk (proportion)
rho = 0.6
rab = rho*sqrt(exp(datos$Ra) * exp(datos$Rb))
datos$Rab <- log(rab)

# --- Create dataframe
R = SpatialPixelsDataFrame(points = coordenadas, data = datos[,3:5])

# --- Plot Dataframe
plot.ra = spplot(R, 'Ra' , col.regions = rev(brewer.pal(n = 11, name = "Spectral")), cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1))
plot.rb = spplot(R, 'Rb' ,  col.regions = rev(brewer.pal(n = 11, name = "Spectral"))  , cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1))
plot.rab = spplot(R, 'Rab' , col.regions = rev(brewer.pal(n = 11, name = "Spectral")) , cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1))
grid.arrange(plot.ra, plot.rb, plot.rab, ncol = 3)


# --- Generate Poisson data
# Population sizes
effort <- "poisson"
if(effort == "uniform"){
  n = floor(runif(size, 500,2000))
}else if(effort == "normal"){
  n = floor(rnorm(size, 1000,500))
  n = n - min(n) + 50
}else if(effort == "poisson"){
  n = floor(rpois(size, 2000))
}else if(effort == "beta"){
  n = floor(1000*rbeta(size,0.1,5))
  n = n + 50
}
hist(n, breaks = 15, main = '')

# Count data realisations (Poisson log-normal)
Ya = rpois(size, exp(log(n) + R$Ra))
Yb = rpois(size, exp(log(n) + R$Rb))
Yab <- rpois(size, exp(log(n) + R$Rab))

# Add data to dataframe
R$Ya = Ya
R$Yb = Yb
R$Yab = Yab
R$Na = n
R$Nb = n

#----- Plot dataframe
color.pal = rev(brewer.pal(n = 11, name = "Spectral"))
color.pal = rev(viridis(20))
plot.ya = spplot(R, 'Ya' , col.regions = color.pal, cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1), main = '')
plot.yb = spplot(R, 'Yb' , col.regions = color.pal, cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1), main = '')
plot.yab = spplot(R, 'Yab' , col.regions = color.pal, cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1), main = '')
grid.arrange(plot.ya, plot.yb, plot.yab, ncol = 3)

# --- Split data frame
set.seed(512)
sample <- sample.int(n = nrow(R), size = floor(.3 * nrow(R)), replace = F)

# --- Smoothing (isotopic) or prediction/smoothing (heterotopic)
iso <- F
if(iso){
  train <- R
  test <-R
}else{
  train <- R[sample,]
  test <- R[-sample,]
}

# --- Plot split datasets
color.pal = rev(brewer.pal(n = 11, name = "Spectral"))
plot.ra = spplot(test, 'Ra' , col.regions = color.pal, cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1), colNa = 'black', main = 'Training Ra')
plot.rb = spplot(R, 'Rb' , col.regions = color.pal, cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1), main = 'Training Rb')
grid.arrange(plot.ra, plot.rb, ncol = 2)


# ------ Define parameters
sizea = nrow(train)
sizeb = nrow(R)
Ya <- train$Ya
Yb <- R$Yb
Yab <- train$Yab
na <- train$Na
nb <- R$Nb
rab <- train$rab
coords.a <- matrix(c(train$x, train$y), ncol= 2)
coords.b <- matrix(c(R$x, R$y), ncol= 2)

# ---------------------------------- 2 -----------------------------------------
# -----------------   Estimating the direct variogram  -------------------------
# ------------------------------------------------------------------------------

# --- Direct Variogram

# Define maximum distance
maxd <- 8
# Define semivariogram bins
nbins <- 10
# Define population rate
poprate <- 1
# Calculate number of observations
No = sizea
# Define population sizes
ni <- na
# Define process
alpha <- T
if(alpha){yi <- Ya}else{yi <- Yb}
rabi <- rab

# Create indexes to use in building of empirical semivariograms
idx1 = rep(1:(No - 1), times = (No - 1):1)
idx2 = unlist(sapply(1:(No - 1), function(i) (1:No)[-(1:i)]))

# distances for unique pairs of points
d = c(dist(coords.a, diag = T))

# Calculate difference of unique pairs of points (adjusted rates - log-normal)
diff = (poprate * (log(yi[idx1] / ni[idx1]) - log(yi[idx2] / ni[idx2]))) ^ 2

# Calculate weighted term
w = (ni[idx1] * ni[idx2]) / (ni[idx1] + ni[idx2])
# bias term
bias = (sum(yi) / sum(ni)) * poprate

# Determine maximum distance of estimator
if (is.null(maxd)) {
  maxd = max(d) / 2
}
# Determine distances of tolerance regions for empirical semivariogram estimator
bindist = seq(0, maxd, len = nbins + 1)
# Determine classification cuts for unique pairs of points by distance
dcuts = cut(d, breaks = bindist)

# Calculate population-weighted empirical semivariogram
np = unlist(lapply(split(d, dcuts), length), use.names = FALSE)
middle = unlist(lapply(split(d, dcuts), mean), use.names = FALSE)
wbinds = split(w , dcuts)
rbinds = split((w * diff - bias), dcuts)
semivariance = unlist(lapply(1:nbins, function(x) {
  sum(rbinds[[x]]) / (2 * sum(wbinds[[x]]))
}), use.names = FALSE)
# Plot population-weighted semivariogram
plot(1:nbins, semivariance)

# Transform varioram to object for gstat (easy manipulation of variogram)
varPka <- data.frame(as.numeric(np) , as.numeric(middle) ,semivariance,rep(0, nbins),
                     rep(0, nbins),rep(as.factor("var1"), nbins))
names(varPka) <-  c("np", "dist", "gamma", "dir.hor", "dir.ver", "id")
varPka <- varPka[!is.na(varPka$gamma),]
rownames(varPka) <- NULL
class(varPka) <- c("gstatVariogram", "data.frame")
plot(varPka)

# Variogram estimation
varPka.mv <- fit.variogram(varPka, vgm(psill = B2[1,1], "Sph", range = 5, add.to = vgm(psill = B1[1,1], "Exp", range = 2)), fit.method = T, fit.ranges = F,warn.if.neg = FALSE)
varioa <- plot(varPka, varPka.mv, col= 'black' , pch = 20, lwd = 2, main = expression('Auto-variogram '~ Y[alpha]/n[alpha]), xlab = 'h', ylab = expression( gamma(h)))
varioa
attr(varPka.mv, "SSErr")
attr(varPka.mv, "singular")

# --- beta

# Define process
alpha <- F
if(alpha){yi <- Ya}else{yi <- Yb}
No = sizeb
ni <- nb
rabi <- rab

# Create indexes to use in building of empirical semivariograms
idx1 = rep(1:(No - 1), times = (No - 1):1)
idx2 = unlist(sapply(1:(No - 1), function(i) (1:No)[-(1:i)]))

# distances for unique pairs of points
d = c(dist(coords.b, diag = T))

# Calculate difference of unique pairs of points
diff = (poprate * (log(yi[idx1] / ni[idx1]) - log(yi[idx2] / ni[idx2]))) ^ 2

# Calculate weighted term
w = (ni[idx1] * ni[idx2]) / (ni[idx1] + ni[idx2])
# bias term
bias = (sum(yi) / sum(ni)) * poprate

# Determine maximum distance of estimator
if (is.null(maxd)) {
  maxd = max(d) / 2
}
# Determine distances of tolerance regions for empirical semivariogram estimator
bindist = seq(0, maxd, len = nbins + 1)
# Determine classification cuts for unique pairs of points by distance
dcuts = cut(d, breaks = bindist)

# Calculate population-weighted empirical semivariogram
np = unlist(lapply(split(d, dcuts), length), use.names = FALSE)
middle = unlist(lapply(split(d, dcuts), mean), use.names = FALSE)
wbinds = split(w , dcuts)
rbinds = split((w * diff - bias), dcuts)
semivariance = unlist(lapply(1:nbins, function(x) {
  sum(rbinds[[x]]) / (2 * sum(wbinds[[x]]))
}), use.names = FALSE)
# Plot population-weighted semivariogram
plot(1:nbins, semivariance)


# Transform varioram to object for gstat (easy manipulation of variogram)
varPkb <- data.frame(as.numeric(np) , as.numeric(middle) , semivariance,rep(0, nbins),
                     rep(0, nbins),rep(as.factor("var2"), nbins))
names(varPkb) <-  c("np", "dist", "gamma", "dir.hor", "dir.ver", "id")
varPkb <- varPkb[!is.na(varPkb$gamma),]
rownames(varPkb) <- NULL
class(varPkb) <- c("gstatVariogram", "data.frame")
plot(varPkb)

# Variogram estimation
varPkb.mv <- fit.variogram(varPkb, vgm(psill = B2[2,2], "Sph", range = 5, add.to = vgm(psill = B1[2,2], "Exp", range = 2)), fit.method = 6, fit.ranges = T)
variob <- plot(varPkb, varPkb.mv, col= 'black' , pch = 20, lwd = 2, 
               main = expression('Auto-variogram '~ Y[beta]/n[beta]), xlab = 'h', ylab = expression( gamma(h)))
variob
attr(varPkb.mv, "SSErr")
attr(varPkb.mv, "singular")


# ---------------------------------- 3 -----------------------------------------
#------------------   Estimating the cross-variograms  -------------------------
#-------------------------------------------------------------------------------

# --- Normal cross-variogram

# Global parameters
No = nrow(train)

# Define population sizes
na <- train$Na
nb <- train$Nb
# Define process
ya <- train$Ya
yb <- train$Yb
# Define total means
Mab <- sum(train$Yab) / sum(train$Na)

# Create indexes to use in building of empirical semivariograms
idx1 = rep(1:No, times = No:1)
idx2 = c(1,unlist(sapply(1:No, function(i) (1:No)[-(1:(i - 1))])))

# distances for unique pairs of points
dist = as.matrix(dist(coords.a))
d = c(dist[lower.tri(dist, diag = T)])

# Calculate difference of unique pairs of points (centered and standardized)
diff = (poprate * ((log(ya[idx1] / na[idx1]) - log(ya[idx2] / na[idx2]))*(log(yb[idx1] / nb[idx1]) - log(yb[idx2] / nb[idx2]))))

# Calculate weighted term
w = (na[idx1] * nb[idx2]) / (nb[idx2] + na[idx1])
# bias term
bias = log(Mab)

# Determine maximum distance of estimator
if (is.null(maxd)) {
  maxd = max(d) / 2
}
# Determine distances of tolerance regions for empirical semivariogram estimator
bindist = seq(0, maxd, len = nbins + 1)
# Determine classification cuts for unique pairs of points by distance
dcuts = cut(d, breaks = bindist)


# Calculate population-weighted empirical semivariogram
np = unlist(lapply(split(d, dcuts), length), use.names = FALSE)
middle = unlist(lapply(split(d, dcuts), mean), use.names = FALSE)
wbinds = split(w , dcuts)
rbinds = split((w * diff - bias), dcuts)
semivariance = unlist(lapply(1:nbins, function(x) {
  sum(rbinds[[x]]) / (2 * sum(wbinds[[x]]))
}), use.names = FALSE)
# Plot population-weighted semivariogram
plot(1:nbins, semivariance)

# Transform varioram to object for gstat (easy manipulation of variogram)
crossvarPk <- data.frame(as.numeric(np) , as.numeric(middle) ,semivariance,rep(0, nbins),
                         rep(0, nbins),rep(as.factor("var1.var2"), nbins))
names(crossvarPk) <-  c("np", "dist", "gamma", "dir.hor", "dir.ver", "id")
crossvarPk <- crossvarPk[!is.na(crossvarPk$gamma),]
rownames(crossvarPk) <- NULL
class(crossvarPk) <- c("gstatVariogram", "data.frame")
plot(crossvarPk)

# Variogram estimation
crossvarPk.mv <- fit.variogram(crossvarPk, vgm(psill = B2[1,2], "Sph", range = 5, add.to = vgm(psill = B1[1,2], "Exp", range = 2)), fit.method = 6, fit.ranges = F)
cross <- plot(crossvarPk, crossvarPk.mv, col= 'black' , pch = 20, lwd = 2, 
              main = expression('Normal cross-variogram '~ Y[alpha]/n[alpha] - Y[beta]/n[beta]), xlab = 'h', ylab = expression( gamma[alpha ~beta](h)))
cross
attr(crossvarPk.mv, "SSErr")
attr(varPkb.mv, "singular")

# Plot adjusted variogarms
grid.arrange(varioa, cross, cross, variob, ncol = 2)

# ---------------------------------- 4 -----------------------------------------
#--------------------------- Poisson Cokriging ---------------------------------


# ---- Auto-variogram and cross-variogram estimation

# define auto-variogram and cross-variograms estimation
variograms <- rbind(crossvarPk,varPkb, varPka)
# create gstat object to handle the variables
g <- gstat(id = "var1", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = train$x ,y = train$y, c = train$Ya, p = train$Na))
g <- gstat(g, id = "var2", formula = c/p ~ 1, loc = ~ x + y, data = data.frame(x = train$x ,y = train$y, c = train$Yb, p = train$Nb))
# Define global model for LMC
g <- gstat(g, id = c("var1","var2"), model = vgm(psill = 1, "Sph", range = 4, add.to = vgm(psill = 1, "Exp", range = 1)))
fittedLMC <- fit.lmc(variograms, g ,model = vgm(psill = 1, "Sph", range = 5, add.to = vgm(psill = 1, "Exp", range = 2)), fit.method = 6, fit.ranges = F)
fittedLMC
plot(variograms,fittedLMC)


# ---- VARIOGRAMS PARAMETERS
para.exp <- data.frame(sill = fittedLMC$model$var1$psill[1],       range = fittedLMC$model$var1$range[1],      nugget = 0)
parb.exp <- data.frame(sill = fittedLMC$model$var2$psill[1],       range = fittedLMC$model$var2$range[1],      nugget = 0)
parab.exp <- data.frame(sill = fittedLMC$model$var1.var2$psill[1], range = fittedLMC$model$var1.var2$range[1], nugget = 0)
para.sph <- data.frame(sill = fittedLMC$model$var1$psill[2],       range = fittedLMC$model$var1$range[2],      nugget = 0)
parb.sph <- data.frame(sill = fittedLMC$model$var2$psill[2],       range = fittedLMC$model$var2$range[2],      nugget = 0)
parab.sph <- data.frame(sill = fittedLMC$model$var1.var2$psill[2], range = fittedLMC$model$var1.var2$range[2], nugget = 0)

# ---- Covariance function
Cmean <- 0

# Verify adjustment
par(mfrow = c(2,2))
plot(dum, (B1[1,1] * Cexp(dum,2,1) + B2[1,1] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ Y[alpha]/n[alpha] - Y[alpha]/n[alpha]))
par(new=TRUE)
plot(dum, Cexp(dum, para.exp$range, para.exp$sill) + Csph(dum, para.sph$range, para.sph$sill), type = 'l', col = 'red', xlab = '', ylab = '')
plot(dum, (B1[1,2] * Cexp(dum,2,1) + B2[1,1] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ Y[alpha]/n[alpha] - Y[alpha]/n[alpha]))
par(new=TRUE)
plot(dum, Cexp(dum, parab.exp$range, parab.exp$sill) + Csph(dum, parab.sph$range, parab.sph$sill), type = 'l', col = 'red', xlab = '', ylab = '')
plot(dum, (B1[2,1] * Cexp(dum,2,1) + B2[1,1] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ Y[alpha]/n[alpha] - Y[alpha]/n[alpha]))
par(new=TRUE)
plot(dum, Cexp(dum, parab.exp$range, parab.exp$sill) + Csph(dum, parab.sph$range, parab.sph$sill), type = 'l', col = 'red', xlab = '', ylab = '')
plot(dum, (B1[2,2] * Cexp(dum,2,1) + B2[1,1] * Csph(dum,5,1)), type = 'l', xlab = 'h' , ylab = 'semivariance', main = expression(~ Y[alpha]/n[alpha] - Y[alpha]/n[alpha]))
par(new=TRUE)
plot(dum, Cexp(dum, parb.exp$range, parb.exp$sill) + Csph(dum, parb.sph$range, parb.sph$sill), type = 'l', col = 'red', xlab = '', ylab = '')

# distances
distMa = as.matrix(dist(coords.a))
distMb = as.matrix(dist(coords.b))
distMab = as.matrix(dist(coords.a, coords.b))
distMba = as.matrix(dist(coords.b, coords.a))

# ---- Error term W_{lk} 

Waa <- diag((sum(train$Ya)/sum(train$Na))/train$Na)
Wbb <- diag((sum(R$Yb)/sum(R$Nb))/R$Nb)

indexab <- which(distMab == 0, arr.ind = TRUE)[,1]
Wab <- matrix(0, dim(distMab)[1], dim(distMab)[2])
Wab[which(distMab == 0, arr.ind = TRUE)] <- (sum(train$Yab)/sum(train$Na))/train$Na[indexab]

indexba <- which(distMba == 0, arr.ind = TRUE)[,2]
Wba <- matrix(0, dim(distMba)[1], dim(distMba)[2])
Wba[which(distMba == 0, arr.ind = TRUE)] <- (sum(train$Yab)/sum(train$Na))/train$Na[indexba]


# ---- Covariance matrices with W_{lk}
Caa <-  (Cmean + Cexp(distMa,  para.exp$range,  para.exp$sill) + Csph(distMa, para.sph$range, para.sph$sill))   + Waa
Cab <-  (Cmean + Cexp(distMab, parab.exp$range, parab.exp$sill) + Csph(distMab, parab.sph$range, parab.sph$sill))  + Wab
Cba <-  (Cmean + Cexp(distMba, parab.exp$range, parab.exp$sill)+ Csph(distMba, parab.sph$range, parab.sph$sill) )  + Wba
Cbb <-  (Cmean + Cexp(distMb,  parb.exp$range,  parb.exp$sill)+ Csph(distMb, parb.sph$range, parb.sph$sill))   + Wbb 

# ----- Unbiasedness constrains
Ctotal <- rbind(cbind(Caa,Cab),cbind(Cba,Cbb))
Ctotal <- cbind(Ctotal, c(rep(1, (sizea)), rep(0, sizeb)), c(rep(0, sizea), rep(1, sizeb)))
Ctotal <- rbind(Ctotal, c(rep(1, sizea), rep(0, sizeb), 0, 0), c(rep(0, sizea), rep(1, sizeb), 0, 0))

# ---- Inverse matrix
Cinv <- solve(Ctotal)

# ---- Prediction at Ra0
pos <- 1
xa0 <- coords.a[pos,]

# --- Distances to Xa0
distFun <- function(xi, x0){sqrt((x0[1]- xi[1])^2 + (x0[2]- xi[2])^2)}
distXoa <- unlist(apply(coords.a,1, distFun, x0 = xa0))
distXob <- unlist(apply(coords.b,1, distFun, x0 = xa0))

#--- Covariances to prediction
covaao <-  (Cmean + Cexp(distXoa, para.exp$range, para.exp$sill)  + Csph(distXoa, para.sph$range, para.sph$sill))
covaao[pos] <- covaao[pos] + para.exp$nugget + para.sph$nugget
covbao <-  (Cmean + Cexp(distXob, parab.exp$range, parab.exp$sill) + Csph(distXob, parab.sph$range, parab.sph$sill))
covbao[pos] <- covbao[pos] + parab.exp$nugget + parab.sph$nugget
coTotal <- c(covaao, covbao, 1, 0)

# --- Calculate weights
lambda <- Cinv%*%coTotal
max(lambda)
# --- Predict one location

predictions <- matrix(nrow = 1 , ncol = 4)
colnames(predictions) <- c('x','y','pred','var')

predictions[1,1] <- c(xa0[1])
predictions[1,2] <- c(xa0[2])
predictions[1,3] <- sum(lambda[1:(length(lambda)-2)]* c(train$Ya/train$Na, R$Yb/R$Nb) * poprate)
predictions[1,4] <- (Cmean + Cexp(0, para.exp$range, para.exp$sill) + Csph(0, para.sph$range, para.sph$sill) + para.exp$nugget) - sum(lambda*coTotal)


# ----- Predict multiple locations
ya <- train$Ya
yb <- R$Yb
na <- train$Na
nb <- R$Nb

# For all locations
poissonCokrige <- function(xao, xta, xtb, Cinv){
  
  distXoa <- unlist(apply(xta,1, distFun, x0 = xao))
  distXob <- unlist(apply(xtb,1, distFun, x0 = xao))
  posa <- which(distXoa == 0)
  posb <- which(distXob == 0)
  
  #--- Covariances to prediction
  covaao <-  (Cmean + Cexp(distXoa, para.exp$range, para.exp$sill)   + Csph(distXoa, para.sph$range, para.sph$sill))
  covbao <-  (Cmean + Cexp(distXob, parab.exp$range, parab.exp$sill) + Csph(distXob, parab.sph$range, parab.sph$sill))
  if(length(pos) > 0){
    covaao[posa] <- covaao[posa] + para.exp$nugget + para.sph$nugget
    covbao[posb] <- covbao[posb] + parab.exp$nugget + parab.sph$nugget
  }
  coTotal <- c(covaao, covbao, 1, 0)
  
  # --- Calculate weights
  lambda <- Cinv%*%coTotal
  
  # predict
  predictions = c(c(xao[1]), c(xao[2]) , sum(lambda[1:(length(lambda)-2)] * c(ya/na, yb/nb) * poprate), (Cmean + Cexp(0, para.exp$range, para.exp$sill) + Csph(0, para.sph$range, para.sph$sill) + para.exp$nugget) - sum(lambda*coTotal))
  names(predictions) <- c('x','y','pred','var')
  return(predictions)
}

preds <- as.data.frame(t(apply(coords.a, 1, poissonCokrige, xta = coords.a, xtb = coords.b, Cinv = Cinv)))
min(preds$pred)
min(preds$var)

# ------- Prediction (CV) -----------

# data
data.test = test
coords.test = matrix(c(data.test$x, data.test$y), ncol= 2)

# distances
distMa = as.matrix(dist(coords.a))
distMb = as.matrix(dist(coords.b))
distMab = as.matrix(dist(coords.a, coords.b))
distMba = as.matrix(dist(coords.b, coords.a))

# ---- Error term W_{lk} 

Waa <- diag((sum(train$Ya)/sum(train$Na))/train$Na)
Wbb <- diag((sum(R$Yb)/sum(R$Nb))/R$Nb)

indexab <- which(distMab == 0, arr.ind = TRUE)[,1]
Wab <- matrix(0, dim(distMab)[1], dim(distMab)[2])
Wab[which(distMab == 0, arr.ind = TRUE)] <- (sum(train$Yab)/sum(train$Na))/train$Na[indexab]

indexba <- which(distMba == 0, arr.ind = TRUE)[,2]
Wba <- matrix(0, dim(distMba)[1], dim(distMba)[2])
Wba[which(distMba == 0, arr.ind = TRUE)] <- (sum(train$Yab)/sum(train$Na))/train$Na[indexba]


# ---- Covariance matrices with W_{lk}
Caa <-  (Cmean + Cexp(distMa,  para.exp$range,  para.exp$sill) + Csph(distMa, para.sph$range, para.sph$sill))   + Waa
Cab <-  (Cmean + Cexp(distMab, parab.exp$range, parab.exp$sill) + Csph(distMab, parab.sph$range, parab.sph$sill))  + Wab
Cba <-  (Cmean + Cexp(distMba, parab.exp$range, parab.exp$sill)+ Csph(distMba, parab.sph$range, parab.sph$sill) )  + Wba
Cbb <-  (Cmean + Cexp(distMb,  parb.exp$range,  parb.exp$sill)+ Csph(distMb, parb.sph$range, parb.sph$sill))   + Wbb 


# ----- Unbiasedness constrains
Ctotal <- rbind(cbind(Caa,Cab),cbind(Cba,Cbb))
Ctotal <- cbind(Ctotal, c(rep(1, (sizea)), rep(0, sizeb)), c(rep(0, sizea), rep(1, sizeb)))
Ctotal <- rbind(Ctotal, c(rep(1, sizea), rep(0, sizeb), 0, 0), c(rep(0, sizea), rep(1, sizeb), 0, 0))

# ---- Inverse matrix
Cinv <- solve(Ctotal)

# prediction
ya = train$Ya
na = train$Na
yb = R$Yb
nb = R$Nb
pred.target = as.data.frame(t(apply(coords.test, 1, poissonCokrige, xta = coords.a, xtb =  coords.b, Cinv = Cinv)))
min(pred.target$var)

# empty dataframe
df = data.frame(matrix(as.numeric(NA), nrow(data.test), 5))
pred.test = SpatialPointsDataFrame(coords.test, df)
names(pred.test) = c('rate.pred','rate.var', 'observed', 'residual', 'zscore')

# predictions
pred.test$rate.pred = pred.target$pred * poprate
pred.test$rate.var = pred.target$var
pred.test$observed = data.test$Ya/data.test$Na * poprate
pred.test$residual = as.numeric(pred.target$pred - pred.test$observed)
pred.test$zscore = as.numeric(pred.test$residual/sqrt(pred.target$var))
pred.test$risk.real = exp(test$Ra)

# plots
predix <- as(pred.test, 'SpatialPixelsDataFrame')
spplot(predix, 'risk.real', at = seq(min(pred.test$risk.real),max(pred.test$risk.real), length.out = 11), col.regions = rev(brewer.pal(n = 11, name = "Spectral")), cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1))
spplot(predix, 'rate.pred', at = seq(min(pred.test$risk.real),max(pred.test$risk.real), length.out = 11),col.regions = rev(brewer.pal(n = 11, name = "Spectral")), cuts = 10, par.settings = list(strip.background = list(col =  'transparent')), colorkey = list(space = "right", height = 1))
spplot(predix, c('rate.var'), col.regions = rev(magma(20)), colorkey = list(space = "right"))



