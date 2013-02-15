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
\item{fun}{The function to be approximated}
\item{dims}{The number of Chebyshev points in each dimension}
\item{intervals}{A list of minimum and maximum values. One for each
  dimension of the hypercube. If NULL, assume [-1,1] in each dimension}
\item{...}{Further arguments to \code{fun}}
}
\value{
A function defined on the hypercube. A Chebyshev
approximation to the function \code{fun}.
}

\details{
  If \code{intervals} is not provided, it is assumed that the
  domain of the function is [-1,1] x [-1,1] x ... x [-1,1].
  Where the number of grid-points are given by \code{dims}

  If \code{intervals} is provided, it should be a \code{list} with
  elements of length 2, providing minimum and maximum for each
  dimension. Arguments to the function will be transformed from these
  intervals into [-1,1] intervals.
}
\seealso{This is just a simple wrapper to \code{chebappx} in case we
  have the whole function rather than only its values in the grid points.
  I.e. it's just \code{chebappx(evalongrid(fun,dims,intervals,...),intervals)}
}

\examples{
## Make a function on [0,2]
f <- function(y) uniroot(function(x) x-y*cos(pi*x^2),lower=0,upper=1)$root
## approximate it
ch <- Vectorize(chebappxf(f,8,intervals=list(c(0,2))))
## plot it:
s <- seq(0,2,length.out=100)
plot(s,ch(s),type='l')
}