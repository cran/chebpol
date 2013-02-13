# Chebyshev transformation.  I.e. coefficients for given function values in the knots.

# The Chebyshev knots of order n on an interval
chebknots <- function(n,interval=NULL) {
  kn <- cos(pi*((1:n)-0.5)/n)
  if(is.null(interval)) return(kn)
  kn*diff(interval)/2 + sum(interval)/2
}

# evaluate a function on a Chebyshev grid
evalongrid <- function(fun,dims,intervals=NULL,...) {
  if(is.null(intervals)) {
    knots <- lapply(dims,chebknots)
  } else {
    if(length(dims) == 1 && is.numeric(intervals)) intervals=list(intervals)
    if(!is.list(intervals)) stop("intervals should be a list")
    knots <- mapply(chebknots,dims,intervals,SIMPLIFY=FALSE)
  }
  structure(apply(expand.grid(knots),1,fun,...), dim=dims)
}


# Chebyshev coefficients for x, which may be an array
chebcoef <- function(x) {
  .Call(C_chebcoef,x)
}

# return a function which is a Chebyshev interpolation
chebappx <- function(val,intervals=NULL) {
  if(is.null(dim(val))) {
    dim(val) <- length(val)
    # allow for vector, e.g. intervals=c(0,1), put it inside list
  }
  if(is.numeric(intervals) && length(intervals) == 2) intervals <- list(intervals)

  cf <- chebcoef(val)
  dims <- dim(val)

  if(is.null(intervals)) {
    # it's [-1,1] intervals, so drop transformation
    cfun <- cmpfun(function(x,d) cos((0:(d-1))*acos(x)))
    imap <- cmpfun(function(x) mapply(cfun,x,dims,SIMPLIFY=FALSE))
    intervals <- replicate(length(dims),c(-1,1),simplify=FALSE)
  } else {
    # it's intervals, create mapping into [-1,1]
    if(!is.list(intervals)) stop("intervals should be a list")
    ispan <- lapply(intervals,function(x) 2/diff(x))
    mid <- lapply(intervals,function(x) sum(x)/2)
    cfun <- cmpfun(function(x,d,is,m) cos((0:(d-1))*acos((x-m)*is)))
    imap <- cmpfun(function(x) mapply(cfun,x,dims,ispan,mid,SIMPLIFY=FALSE))
  }
  structure(function(x) .Call(C_evalcheb,cf,imap(x)),arity=length(dims),domain=intervals)
}

# interpolate a function
chebappxf <- function(fun,dims,intervals=NULL,...) {
  chebappx(evalongrid(fun,dims,intervals,...),intervals)
}

# interpolate on a non-Chebyshev grid. This is useful if you for some reason
# do not have the function values on a Chebyshev-grid, but on some other grid.  It comes at a cost,
# The interpolation may not be very good compared to the Chebyshev-one.

# val are the function values on a grid, an array of appropriate dimension
# in the order of expand.grid()

# if grid is unspecified, it is assumed that it is on a Chebyshev grid in [-1,1]
# If grid is specified, it is a list of vectors. The length of the list
# is the dimension of the grid. Each vector contains grid-points in increasing order.
# val is assumed to be the function values on expand.grid(grid)
# The user-grid is mapped to [-1,1] by approximating the inverse of a one-dimensional Chebyshev
# transform in each dimension. The number of knots in this inverse approximation is
# specified in mapdim.  It defaults to 10 in each dimension.
# if mapdim is NULL, the inversion will be done on each call to the interpolation function

.onAttach <- function(libname,pkgname) {
  if(!.Call(C_hasfftw)) {
    packageStartupMessage("*** ",pkgname,": FFTW not used.\n*** You should install it from http://fftw.org and recompile ",pkgname)
  }
}

chebappxg <- function(val,grid=NULL,mapdim=NULL) {
  # grid is a list of grid points. val is the values as in expand.grid(grid)
  # if grid is null it is assumed to be a chebyshev grid. The dimensions
  # must be present in val
  if(is.null(grid)) return(chebappx(val))
  if(prod(sapply(grid,length)) != length(val)) stop('grid size must match data length')
  dim(val) <- sapply(grid,length)
  # ok, grid is something like list(c(...),c(...),c(...))
  # create piecewise linear functions which maps grid-points to chebyshev grids

  intervals <- lapply(grid,function(x) c(min(x),max(x)))

  if(is.null(mapdim)) {
    gridmaps <- mapply(
                  function(gm) {
                    gm  # force promise to be evaluated
                    function(x) {
                      uniroot(function(y) gm(y)-x, lower=-1,upper=1)$root
                    }
                  },
                  mapply(chebappx,grid)
                )
  } else {
    mapdim <- as.integer(mapdim)
    if(any(mapdim < 1)) stop('mapdim should be >= 1')
    if(length(mapdim) == 1 && length(grid) != 1) {
      mapdim <- rep(mapdim,length(grid))
    } else if(length(mapdim) != length(grid)) {
      stop('length of mapdim ',length(mapdim),' must be the same as length of grid ',length(grid))
    }
    gridmaps <- mapply(
                  chebappxf,
                  mapply(function(gm) {
                    gm  # force promise to be evaluated
                    function(x) uniroot(function(y) gm(y)-x, lower=-1,upper=1)$root
                  },
                         mapply(chebappx,grid)),
                  dims=mapdim,
                  intervals=intervals)
  }
  # now, val is the values on the grid-points 
  # create an ordinary Chebyshev interpolation, but make sure to map
  # user coordinates to [-1,1]
  # Note that gridmap may happen to map outside [-1,1] for small mapdim (perhaps for large?)
  # we simply cut it.

  cutfunc <- function(x) pmin(pmax(x,-1),1)
  gridmap <- cmpfun(function(x) mapply(function(gm,x) cutfunc(gm(x)),gridmaps,x))
  ch <- chebappx(val)
  structure(function(x) ch(gridmap(x)),arity=length(grid),domain=intervals,grid=grid)
}
