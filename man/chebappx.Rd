\name{chebappx}
\alias{chebappx}
\title{Create a Chebyshev interpolation from function values on a grid}
\concept{Chebyshev polynomial}
\description{
  Given function values on a Chebyshev grid, create an approximating function defined
  in the whole hypercube.
}
\usage{
chebappx(val,intervals=NULL)
}

\arguments{
\item{val}{The function values on the grid. \code{val} should be an
  array with appropriate dimension attribute.}
\item{intervals}{A list of minimum and maximum values. One for each
  dimension of the hypercube.}
}
\value{
A function defined on the hypercube. A Chebyshev
approximation to the function values provided in \code{val}.
}

\details{
  If \code{intervals} is not provided, it is assumed that the
  domain of the function is the Cartesian product [-1,1] x [-1,1] x
  ... x [-1,1].  Where the number of grid-points are given by
  \code{dim(val)}.

  If \code{intervals} is provided, it should be a \code{list} with
  elements of length 2, providing minimum and maximum for each
  dimension. Arguments to the function will be transformed from these
  intervals into [-1,1] intervals.

  The approximation function may be evaluated outside the hypercube, but
  be aware that it may be highly erratic there, especially if of high
  degree.
}

\examples{
f <- function(x) exp(-sum(x^2))
## we want 3 dimensions, i.e. something like
## f(x,y,z) = exp(-(x^2 + y^2 + z^2))
## 8 points in each dimension
gridsize <- list(8,8,8)
# get the function values on the Chebyshev grid
values <- evalongrid(f,gridsize)
# make an approximation
ch <- chebappx(values)
## test it:
a <- runif(3,-1,1);ch(a)-f(a)

## then one with domain [0.1,0.3] x [-1,-0.5] x [0.5,2]
intervals <- list(c(0.1,0.3),c(-1,-0.5),c(0.5,2))
# evaluate on the grid
values <- evalongrid(f,gridsize,intervals)
# make an approximation
ch2 <- chebappx(values,intervals)
a <- c(0.25,-0.68,1.43); ch2(a)-f(a)
# outside of domain:
a <- runif(3) ; ch2(a); f(a)
}