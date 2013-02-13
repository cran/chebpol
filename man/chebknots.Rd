\name{chebknots}
\alias{chebknots}
\title{Create a one-dimensional Chebyshev-grid}
\concept{Chebyshev polynomial}
\description{
  Create a one-dimensional Chebyshev grid on an interval.
}

\usage{
chebknots(n,interval=NULL)
}

\arguments{
\item{n}{The number of grid-points. For a Chebyshev-polynomial of degree \code{n-1}}
\item{interval}{A vector of length 2. The lower and upper bound of the interval}
}


\value{
A vector of length \code{n}. The Chebyshev grid-points.
}

\details{
  If \code{interval} is not provided, it is assumed that the
  domain of the function is [-1,1].  Thus, standard Chebyshev
  knots are produced.

  If \code{interval} is provided, it should be a \code{list} with
  elements of length 2, providing minimum and maximum for each
  dimension. 
}

\examples{
## Standard knots for degree 3
chebknots(4)
## Knots in the interval [2,3] for degree 3
chebknots(4,interval=c(2,3))
## Multivariate knots
g <- lapply(c(3,4,3),chebknots); print(g)
## Multivariate grid
\dontrun{
expand.grid(g)
}
}