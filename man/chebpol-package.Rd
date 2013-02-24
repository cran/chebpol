\name{chebpol-package}
\alias{chebpol-package}
\alias{chebpol}
\docType{package}
\title{Methods for creating multivariate Chebyshev interpolations on hypercubes}
\concept{Multivariate Chebyshev approximation}
\concept{Discrete Cosine Transform}
\description{
  The package contains methods for creating multivariate Chebyshev
  interpolations for real-valued functions on hypercubes.  Some methods for
  non-Chebyshev grids are also provided.
}

\details{
  Given a real-valued function on a hypercube, or hyper-rectangle,
  it may be approximated by a (multivariate) Chebyshev polynomial.

  In the one-dimensional case, the Chebyshev approximation is
  constructed by evaluating the function in certain points, a
  Chebyshev grid, and fitting a polynomial to these values.
  Thus, one needs the function values on a set of
  prespecified points.
  
  The multivariate case is similar, the grid is the Cartesian product of
  one-dimensional grids.  I.e. all combinations of all grid-points.

  The Chebyshev coefficients for the interpolating polynomial in the
  one-dimensional case is a simple linear transform of the function
  values.  The Chebyshev-transform, or Discrete Cosine Transform, being
  a variant of the Fourier transform, factors over tensor products, thus
  the multivariate transform is just a tensor product of several
  one-dimensional transforms.  If \href{http://fftw.org}{\acronym{FFTW}} was
  available at compile time, \pkg{chebpol} uses it to generate
  the Chebyshev coefficients, otherwise a slower and more
  memory-demanding matrix method is used, and a warning message is issued when
  the package is attached.

  The Chebyshev-approximation is defined on the interval [-1,1], but it
  is straightforward to map any interval into [-1,1], thus making
  Chebyshev approximation on an interval of choice.  Or, a hypercube of choice.

  The primary method of the package is \code{\link{chebappx}} which takes as
  input the function values on the grid, possibly together with a
  hypercube specification in the form of a list of intervals.  It
  produces a function which interpolates on the hypercube.

  There is also a wrapper called \code{\link{chebappxf}} which may be used if
  one has the function to be approximated rather than only its values in
  the grid-points.

  There is even an interpolation for uniform grids in
  \code{\link{ucappx}}, with a wrapper in \code{\link{ucappxf}}
  with an interesting example.  And a more general for arbitrary
  Cartesian-product grids in \code{\link{chebappxg}} with a wrapper in
  \code{\link{chebappxgf}}.
  
  A multilinear interpolation is available in \code{\link{mlappx}}.
  
  There are also functions for producing Chebyshev
  grids (\code{\link{chebknots}}) as well as a support function for evaluating a
  function on a grid (\code{\link{evalongrid}}), and a function for
  finding the Chebyshev coefficients (\code{\link{chebcoef}}).
}

\examples{
## make some function values on a 3x3x4 grid
dims <- c(x=3,y=3,z=4)
value <- array(runif(36),dims)
##fit a Chebyshev approximation to it. Note that the value-array contains the
##grid-structure.  (However, we don't really approximate the runif-function :)
ch <- chebappx(value)
## To see the full grid, use the chebknots function and expand.grid
\dontrun{
cbind(expand.grid(chebknots(dims)), value=as.numeric(value),
      appx=as.numeric(evalongrid(ch,dims)))
}
## evaluate between the grid points
ch(c(-0.3,0.2,0.16))
}
\keyword{DCT}
\keyword{Discrete Cosine Transform}
\keyword{Chebyshev approximation}
