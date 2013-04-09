
\name{polyh}
\alias{polyh}
\title{Polyharmonic splines on scattered data}
\concept{Polyharmonic splines}
\description{
Polyharmonic splines on scattered data.
}

\usage{
  polyh(val, knots, k=2, ...) 
}

\arguments{
\item{val}{array or function. Function values on scattered data, or the function
  itself.}
\item{knots}{matrix.  Each column is a point in a multidimensional space.}
\item{k}{positive integer or negative numeric. The degree of the polyharmonic spline.}
\item{...}{Further arguments to the function, if \code{is.function(val)}.}
}



\value{
A \code{function(x)} defined on the multidimensional space, approximating the given
function. 
}

\details{
  \code{polyh} fits a polyharminc spline with radial basis function \code{x^k} for odd
  \code{k}, and \code{x^k log(x)} for even \code{k}. If \code{k < 0},
  the basis \code{exp(k x^2)} is used. There are more details in a vignette.

  If \code{val} is a function it will be evaluated on the knots.
}

\examples{
# a function on a 20-dimensional space
f <- function(x) 10/(10+sum(sqrt(x)))
knots <- matrix(runif(16000), 20)
phs <- polyh(f, knots, 3)
# test it in a random point
a <- runif(20)
f(a); phs(a)
}  