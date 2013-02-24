.onAttach <- function(libname,pkgname) {
  if(!havefftw()) {
    packageStartupMessage("*** ",pkgname,": FFTW not used.\n*** You should install it from http://fftw.org\n*** or check if your OS-distribution provides it, and recompile.",pkgname)
  }
}


# Chebyshev transformation.  I.e. coefficients for given function values in the knots.

# The Chebyshev knots of order n on an interval
chebknots1 <- function(n, interval=NULL) {
  kn <- cos(pi*((1:n)-0.5)/n)
  if((n %% 2) == 1) kn[[(n+1)/2]] = 0
  if(is.null(interval)) return(kn)
  kn*diff(interval)/2 + mean(interval)
}

chebknots <- function(dims, intervals=NULL) {
  if(is.null(intervals)) {
    res <- lapply(dims,chebknots1)
  } else {
    if(length(dims) == 1 && is.numeric(intervals)) intervals=list(intervals)
    if(!is.list(intervals)) stop("intervals should be a list")
    res <- mapply(chebknots1,dims,intervals,SIMPLIFY=FALSE)
  }

  res
}

# evaluate a function on a Chebyshev grid
evalongrid <- function(fun,dims,intervals=NULL,...,grid=NULL) {
# do a call to stuff which doesn't really expand the grid
  if(is.null(grid)) grid <- chebknots(dims,intervals)
  mf <- match.fun(fun)
  .Call(C_evalongrid,function(x) mf(x,...), grid)
#  structure(apply(expand.grid(chebknots(dims,intervals)),1,fun,...), dim=dims)
}


# Chebyshev coefficients for x, which may be an array
chebcoef <- function(val, dct=FALSE) {
  structure(.Call(C_chebcoef,as.array(val),dct),dimnames=dimnames(val))
}

chebeval <- function(x,coef,intervals=NULL) {
  if(is.null(intervals)) return(.Call(C_evalcheb,coef,x))
  # map into intervals
  .Call(C_evalcheb,mapply(function(x,i) 2*(x[[1]]-mean(i))/diff(i),x,intervals))
}

# return a function which is a Chebyshev interpolation
chebappx <- function(val,intervals=NULL) {
  if(is.null(dim(val))) {
    # allow for one-dimensional
    dim(val) <- length(val)
  }
   # allow for vector, e.g. intervals=c(0,1), put it inside list
  if(is.numeric(intervals) && length(intervals) == 2) intervals <- list(intervals)

  cf <- chebcoef(val)

  if(is.null(intervals)) {
    # it's [-1,1] intervals, so drop transformation
    fun <- structure(function(x) .Call(C_evalcheb,cf,x),arity=length(dim(val)))
  } else {
    # it's intervals, create mapping into [-1,1]
    if(!is.list(intervals)) stop("intervals should be a list")
    if(any(sapply(intervals,length) != 2)) stop("interval elements should have length 2")
    if(length(intervals) != length(dim(val))) stop("values should have the same dimension as intervals",
               length(intervals),length(dim(val)))
    ispan <- sapply(intervals,function(x) 2/diff(x))
    mid <- sapply(intervals,function(x) mean(x))
    imap <- cmpfun(function(x) (x-mid)*ispan)
    fun <- structure(function(x) .Call(C_evalcheb,cf,imap(x)),arity=length(dim(val)),domain=intervals)
  }
  fun
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


chebappxg <- function(val,grid=NULL,mapdim=NULL) {
  # grid is a list of grid points. val is the values as in expand.grid(grid)
  # if grid is null it is assumed to be a chebyshev grid. The dimensions
  # must be present in val
  if(is.null(grid)) return(chebappx(val))
  if(is.null(dim(val))) dim(val) <- length(val)
  if(!is.list(grid) && length(grid) == length(val)) grid <- list(grid)
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

chebappxgf <- function(fun, grid, ..., mapdim=NULL) {

  if(!is.list(grid)) grid <- list(grid)
  chebappxg(evalongrid(fun, ..., grid=grid),grid,mapdim)
}

ucappx <- function(val, intervals=NULL) {
  if(is.null(dim(val))) dim(val) <- length(val)
  grid <- lapply(dim(val),function(n) seq(-1,1,length.out=n))
  if(is.null(intervals)) return(chebappxg(val,grid))
  chebappxg(val,mapply(function(g,i) 0.5*g*diff(i) + mean(i),grid,intervals,SIMPLIFY=FALSE))
}

ucappxf <- function(fun, dims, intervals=NULL,...) {
  if(is.numeric(dims)) dims <- list(dims)
  grid <- lapply(dims,function(n) seq(-1,1,length.out=n))
  if(is.null(intervals)) return(chebappxgf(fun,grid))
  if(is.numeric(intervals) && length(intervals) == 2) intervals <- list(intervals)
  chebappxgf(fun,mapply(function(g,i) 0.5*g*diff(i) + mean(i),grid,intervals,SIMPLIFY=FALSE))
}

mlappx <- function(val,grid) {
  gl <- prod(sapply(grid,length))
  if(length(val) != gl)
    stop("length of values ",length(val)," do not match size of grid ",gl)
  function(x) .Call(C_evalmlip,grid,as.numeric(val),as.numeric(x))
}

havefftw <- function() .Call(C_havefftw)

