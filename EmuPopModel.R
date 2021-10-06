##########################################################################################################################################
## emu (Dromaius novaehollandiae) demographic model
## 
## Corey Bradshaw
## corey.bradshaw@flinders.edu.au
## Flinders University, September 2021
##########################################################################################################################################

## functions
# beta distribution shape parameter estimator function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

## source
source("matrixOperators.r")


##############################
## DROMAIUS (novaehollandiae) (DN)

# mass
DN.mass <- 55 # Dromaius novaehollandiae (Sales et al. 2007 Avian and Poultry Biology Reviews 18:1–20)

## large, flightless birds (medium) log10(D) = 3.65 – 0.82×log10(M in g) (https://doi.org/10.1111/ecog.04917)
DN.D.pred <- (10^(3.65 - 0.82*(log10(DN.mass*1000))))/2

## max age
## non-volant birds & mammals (Healy K et al. 2014 PRSB)
## log10ls = 0.89 + 0.13log10m (mass in grams; ls = years)
DN.age.max <- round(10^(0.89 + (0.13*log10(DN.mass*1000))), 0)
DN.age.max <- 17 # adjusted downward

## age vector
DN.age.vec <- 0:DN.age.max

## fertility
DN.eggs <- 6.7 * 3.4 / 2 # /2 for females
DN.hatch <- 0.406 * 0.419
DN.F.pred <- DN.eggs * DN.hatch

## age at primiparity
## alpha (https://dx.doi.org/10.1093%2Fgerona%2F62.2.149)
DN.alpha <- ceiling(0.214*(DN.mass*1000)^0.303) # for birds | 4-5 years for ostriches (https://doi.org/10.1007/s11250-009-9428-2)
DN.alpha <- 3 # (http://www.veterinaryworld.org/Vol.2/November/Behavior%20of%20Emu%20bird%20(Dromaius%20novaehollandiae).pdf)

## define m function with age
DN.m.vec <- c(rep(0, DN.alpha-1), rep(0.5*DN.F.pred, round(DN.alpha/2,0)), rep(DN.F.pred, (DN.age.max+1-((DN.alpha-1+round(DN.alpha/2,0))))))
DN.m.sd.vec <- 0.05*DN.m.vec
plot(DN.age.vec, DN.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")

# fit sigmoidal function
# logistic power function y = a / (1+(x/b)^c)
DN.m.dat <- data.frame(DN.age.vec, DN.m.vec)
param.init <- c(0.5, 4, -4)
DN.fit.logp <- nls(DN.m.vec ~ a / (1+(DN.age.vec/b)^c), 
                   data = DN.m.dat,
                   algorithm = "port",
                   start = c(a = param.init[1], b = param.init[2], c = param.init[3]),
                   trace = TRUE,      
                   nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
DN.fit.logp.summ <- summary(DN.fit.logp)
plot(DN.age.vec, DN.m.vec, type="b", pch=19, xlab="age (yrs)", ylab="m")
DN.age.vec.cont <- seq(0,max(DN.age.vec),1)
DN.pred.p.mm <- coef(DN.fit.logp)[1] / (1+(DN.age.vec.cont/coef(DN.fit.logp)[2])^coef(DN.fit.logp)[3])
#DN.pred.p.mm <- ifelse(DN.pred.p.m > 1, 1, DN.pred.p.m)
lines(DN.age.vec.cont, DN.pred.p.mm,lty=2,lwd=3,col="red")

## survival
## mean adult survival (McCarthy et al. 2008 Am Nat)
## ln{-ln[s(t)]} = ln(a) + bln(M) + ln (t)
ln.a.s <- -1.78; b.s <- -0.21 # for birds
DN.s.tran <- ln.a.s + b.s*log(DN.mass*1000) + log(1)
DN.s.ad.yr <- exp(-exp(DN.s.tran))

## calculate lmax from Dillingham et al. 2016 Ecol Appl 26:322-333
# lmax_DIM occurs when lmax.fun returns 0
lmax.fun <- function(lmax, alpha, s, ar.aT) {
  abs(lmax - exp(ar.aT*(alpha + s/(lmax-s) )^(-1) ))
}
# get.lmax calculates lmax_DIM for scalar alpha and scalar/vector s
get.lmax <- function(alpha, s, ar.aT) {
  # Fixed alpha, allows a single value or vector for s
  lmax <- rep(NA, length(s))
  for (i in 1:length(s)) {
    lmax[i] <- optimize(lmax.fun, c(1, 5), tol = 0.0001, alpha = alpha, s=s[i], ar.aT=ar.aT)$minimum
  }
  return(list(lmax=lmax, alpha=alpha, s=s, ar.aT=ar.aT))
}

DN.lm.pred <- (get.lmax(alpha=DN.alpha, s=DN.s.ad.yr, ar.aT=1.107))$lmax # for birds
DN.rm.pred <- log(DN.lm.pred)

# Siler hazard h(x) (Gurven et al. 2007)
a1 <- 1 - (0.55*DN.s.ad.yr) # initial infant mortality rate (also known as αt)
b1 <- 0.105 # rate of mortality decline (also known as bt)
a2 <- 1 - DN.s.ad.yr # age-independent mortality (exogenous mortality due to environment); also known as ct
a3 <- 2.1e-04 # initial adult mortality rate (also known as βt)
b3 <- 0.43 # rate of mortality increase
longev <- DN.age.max
x <- seq(0,longev,1) # age vector
h.x <- a1 * exp(-b1*x) + a2 + a3 * exp(b3 * x) # Siler's hazard model
plot(x,h.x,pch=19,type="l")
plot(x,log(h.x),pch=19,type="l")
l.x <- exp((-a1/b1) * (1 - exp(-b1*x))) * exp(-a2 * x) * exp(a3/b3 * (1 - exp(b3 * x))) # Siler's survival (proportion surviving) model
init.pop <- 10000
lx <- round(init.pop*l.x,0)
len.lx <- length(lx)
dx <- lx[1:(len.lx-1)]-lx[2:len.lx]
qx <- dx/lx[1:(length(lx)-1)]
DN.Sx <- c(0.55*DN.s.ad.yr, 1 - qx)
plot(x, DN.Sx, pch=19, type="l", xlab="age (years)", ylab="Sx")
DN.s.sd.vec <- 0.05*DN.Sx

## create matrix
DN.popmat <- matrix(data = 0, nrow=DN.age.max+1, ncol=DN.age.max+1)
diag(DN.popmat[2:(DN.age.max+1),]) <- DN.Sx[-(DN.age.max+1)]
DN.popmat[DN.age.max+1,DN.age.max+1] <- DN.Sx[DN.age.max+1]
DN.popmat[1,] <- DN.pred.p.mm
colnames(DN.popmat) <- c(0:DN.age.max)
rownames(DN.popmat) <- c(0:DN.age.max)
DN.popmat.orig <- DN.popmat ## save original matrix

## matrix properties
max.lambda(DN.popmat.orig) ## 1-yr lambda
DN.lm.pred
max.r(DN.popmat.orig) # rate of population change, 1-yr
DN.ssd <- stable.stage.dist(DN.popmat.orig) ## stable stage distribution
plot(DN.age.vec, DN.ssd, type="l", pch=19, xlab="age (yrs)", ylab="ssd")
R.val(DN.popmat.orig, DN.age.max) # reproductive value
DN.gen.l <- G.val(DN.popmat.orig, DN.age.max) # mean generation length

## initial population vector
area <- 500*500 # km × km
DN.pop.found <- round(area*DN.D.pred, 0) # founding population size (estimated density * 100 × 100 km region [10,000 km2])
DN.init.vec <- DN.ssd * DN.pop.found

#################
## project
## set time limit for projection in 1-yr increments
yr.st <- 1
#************************
yr.end <- round(40*DN.gen.l, 0) # set projection end date
#************************
t <- (yr.end - yr.st)

DN.tot.F <- sum(DN.popmat.orig[1,])
DN.popmat <- DN.popmat.orig
yr.vec <- seq(yr.st,yr.end)

## set population storage matrices
DN.n.mat <- matrix(0, nrow=DN.age.max+1,ncol=(t+1))
DN.n.mat[,1] <- DN.init.vec

## set up projection loop
for (i in 1:t) {
  DN.n.mat[,i+1] <- DN.popmat %*% DN.n.mat[,i]
}

DN.n.pred <- colSums(DN.n.mat)
yrs <- seq(yr.st, yr.end, 1)
plot(yrs, log10(DN.n.pred),type="l",lty=2,pch=19,xlab="year",ylab="log10 N")

# compensatory density feedback
DN.K.max <- 1*DN.pop.found
DN.K.vec <- c(1, DN.K.max/2, 0.75*DN.K.max, DN.K.max) 
DN.red.vec <- c(1,0.97,0.932,0.885)
plot(DN.K.vec, DN.red.vec,pch=19,type="b")
DN.Kred.dat <- data.frame(DN.K.vec, DN.red.vec)

# logistic power function a/(1+(x/b)^c)
DN.param.init <- c(1, 2*DN.K.max, 2)
DN.fit.lp <- nls(DN.red.vec ~ a/(1+(DN.K.vec/b)^c), 
                 data = DN.Kred.dat,
                 algorithm = "port",
                 start = c(a = DN.param.init[1], b = DN.param.init[2], c = DN.param.init[3]),
                 trace = TRUE,      
                 nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/1024))
DN.fit.lp.summ <- summary(DN.fit.lp)
plot(DN.K.vec, DN.red.vec, pch=19,xlab="N",ylab="reduction factor")
DN.K.vec.cont <- seq(1,2*DN.pop.found,1)
DN.pred.lp.fx <- coef(DN.fit.lp)[1]/(1+(DN.K.vec.cont/coef(DN.fit.lp)[2])^coef(DN.fit.lp)[3])
lines(DN.K.vec.cont, DN.pred.lp.fx, lty=3,lwd=3,col="red")

DN.a.lp <- coef(DN.fit.lp)[1]
DN.b.lp <- coef(DN.fit.lp)[2]
DN.c.lp <- coef(DN.fit.lp)[3]

## compensatory density-feedback deterministic model
## set population storage matrices
DN.n.mat <- matrix(0, nrow=DN.age.max+1, ncol=(t+1))
DN.n.mat[,1] <- DN.init.vec
DN.popmat <- DN.popmat.orig

## set up projection loop
for (i in 1:t) {
  DN.totN.i <- sum(DN.n.mat[,i])
  DN.pred.red <- as.numeric(DN.a.lp/(1+(DN.totN.i/DN.b.lp)^DN.c.lp))
  diag(DN.popmat[2:(DN.age.max+1),]) <- (DN.Sx[-(DN.age.max+1)])*DN.pred.red
  DN.popmat[DN.age.max+1,DN.age.max+1] <- (DN.Sx[DN.age.max+1])*DN.pred.red
  DN.popmat[1,] <- DN.pred.p.mm
  DN.n.mat[,i+1] <- DN.popmat %*% DN.n.mat[,i]
}

DN.n.pred <- colSums(DN.n.mat)
plot(yrs, DN.n.pred, type="l",lty=2,pch=19,xlab="year",ylab="N")
abline(h=DN.pop.found, lty=2, col="red", lwd=2)

## stochatic projection with density feedback
## set storage matrices & vectors
iter <- 100
itdiv <- iter/10

DN.n.sums.mat <- matrix(data=NA, nrow=iter, ncol=(t+1))
DN.s.arr <- DN.m.arr <- array(data=NA, dim=c(t+1, DN.age.max+1, iter))

for (e in 1:iter) {
  DN.popmat <- DN.popmat.orig
  
  DN.n.mat <- matrix(0, nrow=DN.age.max+1,ncol=(t+1))
  DN.n.mat[,1] <- DN.init.vec
  
  for (i in 1:t) {
    # stochastic survival values
    DN.s.alpha <- estBetaParams(DN.Sx, DN.s.sd.vec^2)$alpha
    DN.s.beta <- estBetaParams(DN.Sx, DN.s.sd.vec^2)$beta
    DN.s.stoch <- rbeta(length(DN.s.alpha), DN.s.alpha, DN.s.beta)
    
    # stochastic fertilty sampler (gaussian)
    DN.fert.stch <- rnorm(length(DN.popmat[,1]), DN.pred.p.mm, DN.m.sd.vec)
    DN.m.arr[i,,e] <- ifelse(DN.fert.stch < 0, 0, DN.fert.stch)
    
    DN.totN.i <- sum(DN.n.mat[,i], na.rm=T)
    DN.pred.red <- DN.a.lp/(1+(DN.totN.i/DN.b.lp)^DN.c.lp)
    
    diag(DN.popmat[2:(DN.age.max+1),]) <- (DN.s.stoch[-(DN.age.max+1)])*DN.pred.red
    DN.popmat[DN.age.max+1,DN.age.max+1] <- (DN.s.stoch[DN.age.max+1])*DN.pred.red
    DN.popmat[1,] <- DN.m.arr[i,,e]
    DN.n.mat[,i+1] <- DN.popmat %*% DN.n.mat[,i]

    DN.s.arr[i,,e] <- DN.s.stoch * DN.pred.red
    
  } # end i loop
  
  DN.n.sums.mat[e,] <- ((as.vector(colSums(DN.n.mat))/DN.pop.found))
  
  if (e %% itdiv==0) print(e) 
  
} # end e loop

DN.n.md <- apply(DN.n.sums.mat, MARGIN=2, median, na.rm=T) # mean over all iterations
DN.n.up <- apply(DN.n.sums.mat, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
DN.n.lo <- apply(DN.n.sums.mat, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

par(mfrow=c(1,3))
plot(yrs,DN.n.md,type="l", main = "", xlab="year", ylab="pN1", lwd=2, ylim=c(0.95*min(DN.n.lo),1.05*max(DN.n.up)))
lines(yrs,DN.n.lo,lty=2,col="red",lwd=1.5)
lines(yrs,DN.n.up,lty=2,col="red",lwd=1.5)

DN.s.add <- DN.m.add  <- rep(0, DN.age.max+1)
for (m in 1:iter) {
  DN.s.add <- rbind(DN.s.add, DN.s.arr[ceiling(DN.gen.l):(t+1),,m])
  DN.m.add <- rbind(DN.m.add, DN.m.arr[ceiling(DN.gen.l):(t+1),,m])
}
DN.s.add <- DN.s.add[-1,]
DN.m.add <- DN.m.add[-1,]

DN.s.md <- apply(DN.s.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
DN.s.up <- apply(DN.s.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
DN.s.lo <- apply(DN.s.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(DN.age.vec,DN.s.md,type="l", main = "", xlab="age", ylab="s", lwd=2, ylim=c(0.95*min(DN.s.lo),1.05*max(DN.s.up)))
lines(DN.age.vec,DN.s.lo,lty=2,col="red",lwd=1.5)
lines(DN.age.vec,DN.s.up,lty=2,col="red",lwd=1.5)

DN.m.md <- apply(DN.m.add, MARGIN=2, median, na.rm=T) # mean s over all iterations
DN.m.up <- apply(DN.m.add, MARGIN=2, quantile, probs=0.975, na.rm=T) # upper over all iterations
DN.m.lo <- apply(DN.m.add, MARGIN=2, quantile, probs=0.025, na.rm=T) # lower over all iterations

plot(DN.age.vec,DN.m.md,type="l", main = "", xlab="age", ylab="m", lwd=2, ylim=c(0.95*min(DN.m.lo),1.05*max(DN.m.up)))
lines(DN.age.vec,DN.m.lo,lty=2,col="red",lwd=1.5)
lines(DN.age.vec,DN.m.up,lty=2,col="red",lwd=1.5)
par(mfrow=c(1,1))
