\name{evalongrid}
\alias{evalongrid}
\title{Evaluate a function on a Chebyshev grid}
\concept{Chebyshev polynomial}
\description{
  Evaluate a function on a Chebyshev grid
}

\usage{
evalongrid(fun,dims,intervals=NULL,...)
}

\arguments{
\item{fun}{Multivariate function to be evaluated}
\item{dims}{The number of grid-points in each dimension}
\item{intervals}{A list. Each entry is a vector of length 2 with
  the lower and upper end of the interval in each dimension}
\item{...}{Further arguments to fun}
}

\value{
An array with the value of \code{fun} on each grid point
}

\details{
  If \code{intervals} is not provided, it is assumed that the
  domain of the function is the hypercube [-1,1] x [-1,1] x ... x [-1,1].  
  Thus, the function is evaluated on a standard Chebyshev grid.


  If \code{intervals} is provided, it should be a \code{list} with
  elements of length 2, providing minimum and maximum for each
  dimension. 

  The grid itself may be produced by \code{expand.grid(mapply(chebknots,dims,intervals))}
}

\examples{
f <- function(x) {a <- sum(x^2); ifelse(a == 0,0,exp(-1/a))}
## Standard grid
evalongrid(f,dims=c(3,5))
## Then in certain intervals
evalongrid(f,dims=c(3,5),intervals=list(c(0,1),c(2,3)))
}