\name{chebappxf}
\alias{chebappxf}
\title{Create a Chebyshev interpolation for a function}
\concept{Chebyshev polynomial}
\description{
  Given a function defined on a hypercube, create a Chebyshev
  approximation to it.
}
\usage{
chebappxf(fun, dims, intervals=NULL,...)
}

\arguments{
\item{fun}{The function to be approximated.}
\item{dims}{The number of Chebyshev points in each dimension.}
\item{intervals}{A list of minimum and maximum values. One for each
  dimension of the hypercube. If NULL, assume [-1,1] in each dimension.}
\item{...}{Further arguments to \code{fun}.}
}
\value{
A function defined on the hypercube. A Chebyshev
approximation to the function \code{fun}.
}

\details{
  If \code{intervals} is not provided, it is assumed that the
  domain of the function is [-1,1] x [-1,1] x ... x [-1,1].
  Where the number of grid-points are given by \code{dims}.

  If \code{intervals} is provided, it should be a \code{list} with
  elements of length 2, providing minimum and maximum for each
  dimension. Arguments to the function will be transformed from these
  intervals into [-1,1] intervals.

  This is just a simple wrapper to \code{\link{chebappx}} in case we
  have the whole function rather than only its values in the grid points.
  I.e. it's just \code{\link{chebappx}(\link{evalongrid}(fun, dims,
  intervals, ...), intervals)}.
}

\examples{
# weird function
f <- function(x) ifelse(x==0,0,sin(1/x))
if(havefftw()) {
ch <- Vectorize(chebappxf(f,10000))
a <- runif(5,0.04,0.05); ch(a)-f(a)
}
# Make a function on [0,2] x [0,1]
f <- function(y) uniroot(function(x) x-y[[1]]*cos(pi*x^2),lower=0,upper=1)$root*sum(y^2)
# approximate it
ch <- chebappxf(f,c(12,12),intervals=list(c(0,2),c(0,1)))

# test it:
a <- c(runif(1,0,2),runif(1,0,1)); ch(a); f(a)
}