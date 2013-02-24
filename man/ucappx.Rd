\name{ucappx}
\alias{ucappx}
\title{Create interpolation from function values on a uniform grid}
\description{
  A poor-man's approximation on uniform grids.
  If you for some reason can't evaluate your function on a
  Chebyshev-grid, but instead have a uniform grid, you may use this function
  to create an interpolation.
}

\usage{
ucappx(val, intervals=NULL)
}

\arguments{
\item{val}{Function values on a grid.}
\item{intervals}{List of vectors of length two. Specifying the hypercube
extent in each dimension}
}


\value{
A \code{function(x)} defined on the hypercube, approximating the given
function.
}

\details{
  This is a wrapper to \code{\link{chebappxg}} for unform grids.
  As the example shows, this approximation is better than the Chebyshev
  approximation for some functions.
}

\examples{
# Runge function
f <- function(x) 1/(1+25*x^2)
grid <- seq(-1,1,length.out=15)
val <- f(grid)
uc <- Vectorize(ucappx(val))
# and the Chebyshev
ch <- Vectorize(chebappxf(f,15))
# test it at 10 random points
t(replicate(10,{a<-runif(1,-1,1); c(arg=a, uc=uc(a), true=f(a), cheb=ch(a))}))
}