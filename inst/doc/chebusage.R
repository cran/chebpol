## ----igrid---------------------------------------------------------------
library(chebpol)
par(mfrow=c(2,2))
grid <- list(x=seq(-1,1,length.out=12), y=seq(-1,1,length.out=10))
plot(y~x,data=expand.grid(grid),typ='p',col='blue',main='Uniform 2d-grid')
grid <- chebknots(c(x=12,y=10))
plot(y ~ x, data=expand.grid(grid), typ='p',col='blue',main='Chebyshev 2d-grid')
grid <- list(x=seq(-1,1,length.out=12)^3L, y=sort(c(0,runif(10,-1,1),1)))
plot(y ~ x, data=expand.grid(grid), typ='p',col='blue',main='Irregular 2d-grid')
data <- cbind(x=runif(120,-1,1),y=runif(120,-1,1))
plot(y ~ x, data=data, typ='p', col='blue', main='Scattered data')

## ----fundef, fig.dim=c(4,4)----------------------------------------------
f <- function(x) 1/mean(log1p(0.5 + sin(0.8+2.3*pi*c(0.6,1.4)*(x+0.09)^3)^2))
s <- seq(-1, 1, length.out=1000)
plot(s, sapply(s,f), typ='l')

## ----cheb, fig.dim=c(4,4)------------------------------------------------
ch30 <- ipol(f, dims=30, method='cheb')
ch60 <- ipol(f, dims=60, method='cheb')
plot(s, sapply(s,f), typ='l')
lines(s, ch30(s), col='blue', lty=2)
lines(s, ch60(s), col='red', lty=2)
legend('topleft',c('fun','30 knots','60 knots'),fill=c('black','blue','red'))

## ----uni, fig.dim=c(4,4)-------------------------------------------------
plot(s, sapply(s,f), typ='l')
uni <- ipol(f, dims=20, method='uniform')
grid <- seq(-1,1,len=20)
ml <- ipol(f, grid=grid, method='multilinear')
fh <- ipol(f, grid=grid, method='fh', k=3)
lines(s, uni(s), col='blue', lty=2)
lines(s, ml(s), col='green', lty=2)
lines(s, fh(s, threads=4), col='red', lty=2)
legend('topleft',c('fun','uni','ml','fh'),fill=c('black','blue','green','red'))

## ----echo=FALSE----------------------------------------------------------
set.seed(46)

## ----fig.align='center',fig.dim=c(7,7)-----------------------------------
f1 <- function(x) 1.5/log(5+sin(pi/2*(x[1]^2-2*x[2]^2)))
ch1 <- ipol(f1, dims=c(9,9), method='cheb')
igrid <- list(x=seq(-1,1,len=9), y=seq(-1,1,len=9))
ml1 <- ipol(f1, grid=igrid, method='multilinear')
fh1 <- ipol(f1, grid=igrid, method='fh', k=3)
y <- x <- seq(-1,1,len=200)
testset <- expand.grid(list(x=x,y=y))
data <- cbind(testset,fun= apply(testset,1,f1), cheb=ch1(t(testset)), 
              m.lin=ml1(t(testset)), F.H.=fh1(t(testset)))

lattice::levelplot(m.lin+F.H.+fun+cheb ~ x+y, data=data, cuts=10,
          col.regions=grey.colors(11,0.1,0.9), layout=c(2,2), 
          main='Level plots of function and interpolations')

## ----poly,echo=FALSE-----------------------------------------------------
set.seed(43)

## ----fig.align='center', fig.width=7, fig.height=7-----------------------
f2 <- function(x) 1.2/min(log(3.1 + sin(0.7+1.8*pi*(x+0.39*x[1]^2-x[2]^2))^2))
ph1 <- ipol(f1, knots=matrix(runif(100,-1,1),2), method='polyharmonic', k=3)
ph2 <- ipol(f2, knots=matrix(runif(4000,-1,1),2), method='polyharmonic', k=3)
data <- cbind(testset,fun1= apply(testset,1,f1), poly1=ph1(t(testset)), 
              fun2=apply(testset,1,f2), poly2=ph2(t(testset)))

lattice::levelplot(fun2+poly2+fun1+poly1 ~ x+y, data=data, cuts=10,
          col.regions=grey.colors(11,0.1,0.9), layout=c(2,2), 
          main='Level plots of functions and interpolations')
# compute L2-error (rmse)
sqrt(cubature::hcubature(function(x) as.matrix((ph2(x)-apply(x,2,f2))^2), rep(-1,2), rep(1,2),
          vectorInterface=TRUE, absError=1e-6)$integral/4)

## ----para----------------------------------------------------------------
m <- matrix(runif(2*6,-1,1),2)
print(m)
print(ph2(m, threads=3))

## ----fast----------------------------------------------------------------
f <- ipol(sin, grid=seq(0,1,length.out=1000), method='multilinear')
a <- runif(1e6)
system.time(sin(a))
system.time(f(a,threads=4))

