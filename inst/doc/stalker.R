## ----fig.dim=c(4,4),echo=FALSE-------------------------------------------
library(chebpol)
set.seed(42)
N <- 10
pts <- as.numeric(seq(1,N))
val <- runif(N)
val[5] <- val[4]
val[6] <- val[6] - 0.3
val[7] <- val[8] + 0.00001
s <- seq(1,N,len=1000)
plot(pts,val,pch=20,xlab='x',ylab='y')
ml <- ipol(val,grid=pts,method='multi')
lines(s,ml(s))

## ----out.width='0.5\\linewidth', echo=FALSE------------------------------
plot(pts,val,pch=20,xlab='x',ylab='y')
for(i in seq_along(pts)) {
    lines(pts[i]+c(-1,1),c(val[i],val[i]))
}

## ----fig.dim=c(5,5), echo=FALSE------------------------------------------
plot(pts,val,pch=20,xlab='x',ylab='y',ylim=c(0,1.1))
for(i in seq_along(pts)[c(-1,-length(pts))]) {
    a <- val[i]
    b <- (val[i+1] - val[i-1])/2
    c <- val[i+1]-a-b
    ss <- seq(-1,1,len=20)
    lines(ss+i,a+b*ss+c*ss^2,col=c('blue','green','red')[i%%3+1])
}

## ----basis,echo=FALSE,fig.dim=c(7,5), fig.align='center', fig.cap='"Basis functions"'----
s <- seq(-1,1,len=100)
plot(range(s),c(-0.7,1.7),pch='',xlab='x',ylab='y')
v0 <- seq(-0.7,1.7,len=41)
f0spline <- function(f0) {
  b <- 0.5
  c <- 0.5*(1-2*f0)
  if(b < 2*abs(c) && abs(c) < b)
    r <- abs(b/c)
  else if (abs(c) < 2*b && abs(b) < abs(c))
    r <- abs(c/b)
  else
    r <- 2
  f0 + b*s + c*abs(s)^r
}
for(f0 in v0) {
  col <- if( (f0 < 0.75 && f0 > 0.25) || (f0 > 1.5 || f0 < -0.5) ) 
           'darkblue' 
         else if(f0 < 1 && f0 > 0)
           'blue'
         else
           'lightblue'
  lines(s,f0spline(f0),col=col)
  points(0,f0,pch=20)
}
# illustrate overshoot
clip(0,0.5,-2,2)
abline(h=-0.7,lty=2)
clip(-0.5/(2*(0.5+0.7)),0.5,-2,2)
abline(h=-0.7-0.25/(4*(0.5+0.7)), lty=2)
clip(-2,2,-2,2)
text(x=0.7,y=-0.63,'} overshoot',pos=1)

## ----failure,echo=FALSE,fig.dim=c(7,5), fig.align='center', fig.cap='"Basis functions"'----
s <- seq(-1,1,len=100)
vseq <- seq(-1,1,len=41)
rr <- m <- v <- list()
for(i in seq_along(vseq)) {
  #  double iD = 1.0/(dmin*pow(dplus,r) + pow(dmin,r)*dplus);
  #  double b = (vplus*pow(dmin,r) - vmin*pow(dplus,r))*iD;
  #  double c = (vplus*dmin + vmin*dplus)*iD;
  b <- 0.5*(1-vseq[[i]])
  c <- 0.5*(1+vseq[[i]])
  if(b < 2*abs(c) && abs(c) <= b)
    r <- abs(b/c)
  else if (abs(c) < 2*b && abs(b) < abs(c))
    r <- abs(c/b)
  else
    r <- 2

  rr[[i]] <- r
  m[[i]] <- abs(c) <= b
  v[[i]] <- b*s + c*abs(s)^r
}
plot(range(s),do.call(range,c(v,list(finite=TRUE))),pch='',xlab='x',ylab='y')
for(i in seq_along(v)) {
  cl <- if(identical(rr[[i]],2)) 'darkblue' else if(isTRUE(m[[i]])) 'blue' else 'lightblue'
  lines(s,v[[i]],col=cl)
}
points(rep(-1,length(vseq)),vseq,pch=20)

## ----blending, fig.dim=c(4,4),fig.cap='Blending functions'---------------
sigmoid <- function(t) ifelse(t<0.5, 0.5*exp(2-1/t), 1-0.5*exp(2-1/(1-t)))
cubic <- function(t) -2*t^3 + 3*t^2
linear <- function(t) t
s <- seq(0,1,length=100)
plot(s,sigmoid(s),typ='l',ylab='y')
lines(s,linear(s), col='blue')
lines(s,cubic(s), col='green')
legend('topleft',legend=c('sigmoid','linear','cubic'),fill=c('black','blue','green'))

## ----echo=FALSE,fig.dim=c(6,4)-------------------------------------------
pts <- pts-min(pts)
pts <- pts/max(pts)
plot(pts,val,pch=20,xlab='x',ylab='y',ylim=c(0,1.1))
ns <- splinefun(pts,val,method='natural')
ms <- splinefun(pts,val,method='mono')
st <- ipol(val,grid=pts,method='stalker',k=NA)
s <- seq(0,1,len=1000)
lines(s,ns(s),col='blue')
lines(s,ms(s),col='green')
lines(s, st(s), col='magenta')
lines(s, st(s,degree=0,blend='cub'),col='red')
legend('topright',legend=c('stalker','mono','natural','hyperbolic'),
       fill=c('magenta','green','blue','red'))

## ----fig.dim=c(6,4),echo=FALSE-------------------------------------------
val <- sort(val)
plot(pts,val,pch=20,xlab='x',ylab='y',ylim=c(0,1.1))
ns <- splinefun(pts,val,method='natural')
ms <- splinefun(pts,val,method='mono')
st <- ipol(val,grid=pts, method='stalker', k=NA)
s <- seq(0,1,len=1000)
lines(s,ns(s),col='blue')
lines(s,ms(s),col='green')
lines(s, st(s), col='magenta')
lines(s, st(s,deg=0,blend='cub'), col='red')
legend('topleft',legend=c('stalker','mono','natural','hyperbolic'),
       fill=c('magenta','green','blue','red'))

## ----fig.dim=c(6,4),echo=FALSE-------------------------------------------
val <- rep(runif(N/2),each=2) + c(0,1e-16)
plot(pts,val,pch=20,xlab='x',ylab='y',ylim=c(0,1.1))
ns <- splinefun(pts,val,method='natural')
ms <- splinefun(pts,val,method='mono')
st <- ipol(val,grid=pts, method='stalker',k=NA)
s <- seq(0,1,len=1000)
lines(s,ns(s),col='blue')
lines(s,ms(s),col='green')
lines(s, st(s), col='magenta')
lines(s, st(s,deg=0,blend='cub'), col='red')
legend('topleft',legend=c('stalker','mono','natural','hyperbolic'),
       fill=c('magenta','green','blue','red'))

## ----volcano, fig.dim=c(4,4), fig.align='center', fig.pos='!ht', fig.cap='Maungawhau', out.width='.37\\linewidth', fig.ncol=2, fig.subcap=c('low resolution','multilinear','stalker','thin plate spline')----
data(volcano)
volc <- volcano[seq(1,nrow(volcano),3),seq(1,ncol(volcano),3)]/10 #low res volcano
grid <- list(x=as.numeric(seq_len(nrow(volc))), y=as.numeric(seq_len(ncol(volc))))
ph <- ipol(volc, grid=grid, method='polyharmonic',k=2)
st <- ipol(volc, grid=grid, method='stalker',k=NA)
ml <- ipol(volc, grid=grid, method='multilinear')
g <- list(x=seq(1,nrow(volc), len=71), y=seq(1,ncol(volc),len=71))
par(mar=rep(0,4)); col <- 'green'
light <- list(specular=0.2,ambient=0.0,diffuse=0.6)
plot3D::persp3D(grid$x, grid$y, volc, colvar=NULL, lighting=light,
        theta=45, ltheta=0, lphi=40, col=col, axes=FALSE, bty='n',scale=FALSE)
for(f in list(ml, st, ph)) {
  plot3D::persp3D(g$x, g$y, evalongridV(f,grid=g), colvar=NULL, lighting=light,
        theta=45, ltheta=0, lphi=40, col=col, axes=FALSE, bty='n', scale=FALSE)
}

## ----random, fig.dim=c(3.5,3.5), fig.pos='!ht', fig.align='center',fig.cap='Random surface',out.width='.4\\linewidth',fig.ncol=2,fig.subcap=c('stalker', 'thin plate spline','hyperbolic stalker','Floater-Hormann')----
set.seed(42); N <- 8
grid <- list(x=seq(0,1,length=N)+c(0,rnorm(N-2,sd=0.3/N),0), 
             y=seq(0,1,length=N)+c(0,rnorm(N-2,sd=0.3/N),0))
val <- matrix(runif(N*N,0,0.3),N)
st <- ipol(val,grid=grid, method='stalker',k=NA)
ph <- ipol(val,grid=grid, method='polyharmonic', k=2)
fh <- ipol(val,grid=grid, method='fh', k=0)
sthyp <- function(x) st(x,deg=0,blend='cubic')
g <- list(x=seq(0,1, len=70), y=seq(0,1,len=70))
par(mar=rep(0,4))
for(f in list(st, ph, sthyp, fh)) {
  plot3D::persp3D(g$x, g$y, evalongridV(f,grid=g), colvar=NULL, lighting=light,
         theta=60, ltheta=30, lphi=45, col='green', axes=FALSE, bty='n', scale=FALSE,zlim=c(0,1))
  pts <- evalongridV(f,grid=grid)+0.00
  plot3D::points3D(rep(grid$x,N),rep(grid$y,each=N),pts,add=TRUE,colvar=NULL,pch=20)
}

## ----torsion, fig.dim=c(5,5), out.width='\\linewidth', fig.cap='Torsion problem in the wet paper spline', fig.align='center'----
g <- list(x=seq(grid$x[5],grid$x[8],length=100),y=seq(grid$y[5],grid$y[8],length=100))
par(mar=rep(1,4))
plot3D::persp3D(g$x, g$y, evalongridV(st,grid=g), colvar=NULL, lighting=light,
       theta=120, ltheta=100, lphi=45, col='green', axes=TRUE)
subgrid <- list(x=grid$x[5:8],y=grid$y[5:8])
zval <- evalongridV(f,grid=subgrid)+0.005
plot3D::points3D(rep(subgrid$x,4),rep(subgrid$y,each=4),zval,add=TRUE,colvar=NULL,pch=20)

## ----tordet, fig.cap='Torsion details', fig.align='center'---------------
x <- seq(grid$x[[5]],grid$x[[8]],len=100)
pts <- c(grid$y[5:8])
pt <- mean(grid$y[6:7])
col <- c('green','magenta','blue','red')
plot(x,st(rbind(x,pt)),typ='l',ylim=c(0,0.3),ylab='z')
lines(x,st(rbind(x,pt),deg=c(1.5,2)), col='grey')
points(grid$x,st(rbind(grid$x,pt)),pch=20)
for(i in seq_along(pts)) {
  v <- st(rbind(x,pts[i]))
  lines(x,v,col=col[i],lty=2)
  points(grid$x,st(rbind(grid$x,pts[i])),pch=20,col=col[i])
}
legend('topright',title='y=',
    legend=round(c(pts[1:2],pt,pts[3:4]),3),
    fill=c(col[1:2],'black',col[3:4]))

## ----constant, fig.dim=c(5,5), out.width='\\linewidth',fig.cap='Constant degree 1.5', fig.align='center'----
par(mar=rep(1,4))
plot3D::persp3D(g$x, g$y, evalongridV(st,grid=g,degree=c(0,1.5)), colvar=NULL, lighting=light,
       theta=120, ltheta=100, lphi=45, col='green', axes=TRUE)
plot3D::points3D(rep(subgrid$x,4),rep(subgrid$y,each=4),zval,add=TRUE,colvar=NULL,pch=20)

## ----eval=FALSE----------------------------------------------------------
#  st <- ipol(val,grid=grid,method='stalker',k=1.5)
#  st(x,degree=1.2,blend='linear')

