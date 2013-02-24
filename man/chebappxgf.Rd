\name{chebappxgf}
\alias{chebappxgf}
\title{Create a non-Chebyshev interpolation for a function}
\concept{Chebyshev polynomial}
\description{
  Given a function defined on a hypercube, create an 
  approximation for it on a non-Chebyshev grid.
}
\usage{
chebappxgf(fun, grid, ..., mapdim=NULL) 
}

\arguments{
\item{fun}{The function to be approximated.}
\item{grid}{A list. Each element is a sorted vector of grid-points for a
  dimension. These need not be Chebyshev-knots, nor evenly spaced.}
\item{...}{Further arguments to \code{fun}.}
\item{mapdim}{Degree of Chebyshev-polynomial mapping grid into Chebyshev
  intervals [-1,1].}
}
\value{
A function defined on the hypercube. An
approximation to the function \code{fun}.
}

\details{This is just a simple wrapper to \code{\link{chebappxg}} in case we
  have the whole function rather than only its values in the grid points.
  I.e. it's just   \code{\link{chebappxg}(\link{evalongrid}(fun, ...,
  grid=grid), grid, mapdim)}.
}
\seealso{This function has the same relation to \code{\link{chebappxg}},
 as \code{\link{chebappxf}} has to \code{\link{chebappx}}.
}

\examples{
# Try Runge's function on a uniformly spaced grid.
# Ordinary polynomial fitting of high degree of Runge's function on a uniform grid
# creates large oscillations near the end of the interval. Not so with chebappxgf
f <- function(x) 1/(1+25*x^2)
chg <- Vectorize(chebappxgf(f,seq(-1,1,length.out=15)))
# also compare with Chebyshev interpolation
ch <- Vectorize(chebappxf(f,15))
\dontrun{
# plot it
s <- seq(-1,1,length.out=200)
plot(s, f(s), type='l', col='black')
lines(s, chg(s), col='blue')
lines(s, ch(s), col='red')
legend('topright',
       legend=c('Runge function','chebappxg on uniform grid','Chebyshev'),
       col=c('black','blue','red'), lty=1)
}
}