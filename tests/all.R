library(chebpol)

cat("Has FFTW:",.Call(chebpol:::C_hasfftw),'\n')
set.seed(42)

a <- array(rnorm(24),c(2,3,4))
chebcoef(a)

f <- function(x) exp(-sum(x^2))

dims <- c(8,7,6)
ch <- chebappxf(f,dims)
s <- runif(3,-1,1)
ch(s)-f(s)
iv <- list(c(1,2),c(1,4),c(1,3))
ch <- chebappxf(f,dims,iv)
s <- c(1.4,2.3,1.9)
ch(s) - f(s)

sum(evalongrid(f,dims))

chebknots(17)

# test chebappxg
## evenly spaced grid-points
su <- seq(0,1,length.out=10)
## irregularly spaced grid-points
s <- su^3
## create approximation on the irregularly spaced grid
ch <- Vectorize(chebappxg(exp(s),list(s)))
## test it:
r <- runif(1); cat('true:',exp(r),'appx:',ch(r),'\n')

#multivariate chebappxg
su <- seq(0,1,length.out=11)
grid <- list(su,su^2,su^3)
dims <- lapply(grid,length)
fv <- structure(apply(expand.grid(grid),1,f),dim=lapply(grid,length))
ch <- chebappxg(fv,grid)
s <- runif(3)
cat('true:',f(s),'appx:',ch(s),'\n')
ch <- chebappxg(fv,grid,mapdim=7)
cat('true:',f(s),'appx:',ch(s),'\n')
ch <- chebappxg(fv,grid,mapdim=lapply(grid,length))
cat('true:',f(s),'appx:',ch(s),'\n')


