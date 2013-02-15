\name{chebappxg}
\alias{chebappxg}
\title{Approximate on a non-Chebyshev grid}
\concept{Chebyshev polynomial}
\description{
  A poor-man's approximation on non-Chebyshev grids.
  If you for some reason can't evaluate your function on a
  Chebyshev-grid, but instead have some other grid which still is a
  Cartesian product of one-dimensional grids, you may use this function
  to create an interpolation.
}

\usage{
chebappxg(val,grid=NULL,mapdim=NULL)
}

\arguments{
\item{val}{Function values on a grid}
\item{grid}{A list. Each element is a vector of grid-points for a
  dimension. These need not be Chebyshev-knots, nor evenly spaced}
\item{mapdim}{Degree of Chebyshev-polynomial mapping grid into Chebyshev
  intervals [-1,1]}
}


\value{
A \code{fun <- function(x)} defined on the hypercube, approximating the given
function.
}

\details{
  A call \code{fun <- chebappxg(val,grid)} does the following.
  A Chebyshev interpolation \code{ch} for \code{val} is created, on the
  [-1,1] hypercube.
  For each dimension a grid-map function \code{gm} is created which maps Chebyshev knots
  into the given grid-points. When \code{fun(x)} is called, the equation
  \code{gm(y)=x} is solved for \code{y}, yielding a number in [-1,1],
  which is used as argument: \code{ch(y)}.

  If \code{mapdim} is an integer, the system \code{gm(y)=x} is solved and approximated with
  a Chebyshev-polynomial of dimension \code{mapdim}, yielding an
  inverse \code{y = igm(x))}. When \code{fun(x)} is called, it
  translates to \code{ch(igm(x))}, thus avoiding inverting \code{gm} in
  each call to \code{fun}.  \code{mapdim} may also be a vector with a
  Chebyshev-length for \code{igm} for each dimension.

  Note that even if \code{fun} is an approximation to the original
  function, and differentiable, it is not as good an
  approximation as if we had a Chebyshev-grid.
}

\examples{
## evenly spaced grid-points
su <- seq(0,1,length.out=10)
## irregularly spaced grid-points
s <- su^3
## create approximation on the irregularly spaced grid
ch <- Vectorize(chebappxg(exp(s),list(s)))
## test it:
ch(su) - exp(su)
# try one with three variables
f <- function(x) exp(-sum(x^2))
grid <- list(s,su,su^2)
ch2 <- chebappxg(evalongrid(f,grid=grid),grid)
# test it at 10 random points
replicate(10,{a<-runif(3); ch2(a)-f(a)})
}