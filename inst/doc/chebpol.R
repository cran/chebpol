### R code from vignette source 'chebpol.Rnw'

###################################################
### code chunk number 1: runge
###################################################
library(chebpol)


###################################################
### code chunk number 2: runge
###################################################
f <- function(x) cos(3*pi*x)/(1+25*(x-0.25)^2)
ch <- ipol(f,dims=15,method='chebyshev')
s <- seq(-1,1,length.out=401)
plot(s, f(s),type='l')
lines(s, ch(s), col='red')
kn <- chebknots(15)[[1]]
points(kn,f(kn))


###################################################
### code chunk number 3: runge
###################################################
pdf('ucrunge.pdf')
plot(s,f(s),type='l')
lines(s,ch(s), col='red')
points(kn,f(kn))


###################################################
### code chunk number 4: runge
###################################################
uc <- ipol(f, dims=15, method='uniform')
lines(s,uc(s),col='blue')


###################################################
### code chunk number 5: runge
###################################################
invisible(dev.off())


###################################################
### code chunk number 6: multi
###################################################
library(chebpol)
f <- function(x) log(x[[1]])*sqrt(x[[2]])/log(sum(x))
ch <- ipol(f, dims=c(5,8), intervals=list(c(1,2), c(15,20)), method='chebyshev')
uc <- ipol(f, dims=c(5,8), intervals=list(c(1,2), c(15,20)), method='uniform')
tp <- c(runif(1,1,2), runif(1,15,20))
cat('arg:',tp,'true:', f(tp), 'ch:', ch(tp), 'uc:',uc(tp),'\n')


###################################################
### code chunk number 7: sqrt
###################################################
f <- function(x) cos(3*pi*x)/(1+25*(x-0.25)^2)
gr <- log(seq(exp(-1),exp(1),length=15))
chg <- ipol(f, grid=gr, method='general')
plot(s, f(s), col='black', type='l')
lines(s, chg(s), col='blue')
points(gr,f(gr))


###################################################
### code chunk number 8: lagrange
###################################################
f <- function(x) 1/(1+25*x^2)
# Uniform grid:
unigrid <- list(seq(-3,2,length.out=20))
uni <- ipol(f, grid=unigrid, k=2, method='fh')
ch <- ipol(f, dims=20, intervals=c(-3,2), method='chebyshev')
s <- seq(-3,2,length.out=1000)
plot(s,f(s),ylim=range(uni(s),f(s),ch(s)),typ='l')
lines(s,uni(s),col='blue')
points(unigrid[[1]],uni(unigrid[[1]]),col='blue',pch=20)
lines(s,ch(s),col='green')
legend('topleft',c('function','FH','chebyshev'),fill=c('black','blue','green'))


###################################################
### code chunk number 9: ml
###################################################
f <- function(x) sign(sum(x^3)-0.1)*
                  sqrt(abs(25*prod(x)-4))/
                  (1+25*sum(x)^2)
grid <- replicate(4,list(seq(-1,1,length=15)))
ml <- ipol(f, grid=grid, method='multilinear')
s <- seq(-1,1,length=400)
curve <- function(x) c(cos(1.2*pi*x),
                       sin(1.5*pi*x^3), 
                       x^2, -x/(1+x^2))
wf <- sapply(s,function(x) f(curve(x)))
wml <- sapply(s,function(x) ml(curve(x)))
plot(s,wf,typ='l')  # function
lines(s,wml,col='blue') # multilinear interpolation


###################################################
### code chunk number 10: poly
###################################################
set.seed(425)  # make sure we are reproducible


###################################################
### code chunk number 11: poly
###################################################
r <- runif(20)
r <- r/sum(r)
f <- function(x) 1/mean(log1p(r*x))
knots <- matrix(runif(60000), 20)
phs <- ipol(f, knots=knots, k=3, method='polyharmonic')
rr <- runif(20)
curve <- function(x) abs(cos(5*pi*rr*x))
s <- seq(0,1,length.out=1000)
plot(s,sapply(s,function(x) f(curve(x))),typ='l')
lines(s,sapply(s,function(x) phs(curve(x))),col='blue',lty=2)


