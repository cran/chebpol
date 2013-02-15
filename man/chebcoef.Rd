\name{chebcoef}
\alias{chebcoef}
\title{Create the Chebyshev-coefficients given values on a Chebyshev grid}
\concept{Chebyshev polynomial}
\description{
  Create the multivariate Chebyshev-coefficients, given values on a Chebyshev grid.
}

\usage{
chebcoef(x)
}

\arguments{
  \item{x}{An array of function values on a Chebyshev grid.
    The \code{dim}-attribute must be appropriately set.
    If not set, it is assumed to be one-dimensional.
  }
}

\value{
  An array of Chebyshev-coefficients for an interpolating Chebyshev-polynomial.
}

\details{
  If \code{x} has no \code{dim}-attribute, it is assumed to be
  one-dimensional of length the length of \code{x}.
}

\examples{
## Coefficients for a 2x3x4 grid
a <- array(rnorm(24),dim=c(2,3,4))
chebcoef(a)
}