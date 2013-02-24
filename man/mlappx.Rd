\name{mlappx}
\alias{mlappx}
\title{Multilinear interpolation on a non-Chebyshev grid}
\concept{Multilinear interpolation}
\description{
  Multilinear interpolation on an arbitrary Cartesian product.
}

\usage{
  mlappx(val,grid)
}

\arguments{
\item{val}{Function values on a grid.}
\item{grid}{A list.  Each element is a vector of grid-points for a
  dimension.  These need not be Chebyshev-knots, nor evenly spaced.}
}


\value{
A \code{function(x)} defined on the hypercube, approximating the given
function.  The function yields values for arguments outside the
hypercube as well, as a linear extension.
}

\details{
  A call \code{fun <- mlappx(val,grid)} creates a multilinear
  interpolation on the grid.  The value on the grid points will be
  exact, the value between the grid points is a convex combination of
  the values in the corners of the hypercube surrounding it.
}

\examples{
## evenly spaced grid-points
su <- seq(0,1,length.out=10)
## irregularly spaced grid-points
s <- su^3
## create approximation on the irregularly spaced grid
ml1 <- Vectorize(mlappx(exp(s),list(s)))
## test it, since exp is convex, the linear approximation lies above
## the exp between the grid points
ml1(su) - exp(su)

## multi linear approx
f <- function(x) exp(sum(x^2))
grid <- list(s,su)

ml2 <- mlappx(evalongrid(f,grid=grid),grid)
a <- runif(2); ml2(a); f(a)
# we also get an approximation outside of the domain, of disputable quality
ml2(c(1,2)); f(c(1,2))

}