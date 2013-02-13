#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#ifdef HAVE_FFTW
#include <fftw3.h>
#endif
/*
  The xvec is an array of pointers to vectors
  Vector i is of length d=dims[i] and contains
  the 1-dimensional Chebyshev vector on coordinate i of x
  I.e. the d-vector cos((0:(d-1))*acos((x[i]))) where
  x has been interval transformed into [-1,1]

  pass dimlen by address for easier call from Fortran.
  Though, the xvec may be problematic in F77. I.e. an array of vectors of different lengths.
  Perhaps a separate Fortran interface should be made?
  However, it is probably easier to write the whole evalcheb in Fortran, if one so desires.
*/

static double evalcheb(double *cf, double *xvec[], int *dims, int *dimlen) {
  double res = 0.0;
  double *xl;
  size_t siz = 1;
  int rank = *dimlen;
  int newrank = rank-1;
  int N = dims[newrank];
  for(int i = 0; i < rank; i++) siz *= dims[i];

  if(newrank == 0) {
    for(int i = 0; i < N; i++) res += cf[i]*xvec[0][i];
    return res;
  }
  // loop over last dimension, recurse
  int len = siz/N;
  xl = xvec[newrank];
  for(int i = 0,j=0; i < N; i++,j+=len) {
    res += xl[i]*evalcheb(&cf[j],xvec,dims,&newrank);
  }
  return res;
}

static void chebcoef(double *x, int *dims, int *dimlen, double *F) {
  size_t siz = 1;
  int rank = *dimlen;
  for(int i = 0; i < *dimlen; i++) siz *= dims[i];

#ifdef HAVE_FFTW
  double isiz = 1.0/siz;
  int rdims[rank];
  /* Create a plan */
  fftw_r2r_kind kind[rank];
  for(int i = 0; i < rank; i++) kind[i] = FFTW_REDFT10;  // type II DCT

  // reverse the dimensions. fftw uses row-major order. And stuff.
  for(int i = 0; i < rank; i++) rdims[i] = dims[rank-i-1];

  // Plan and execute
  fftw_plan plan = fftw_plan_r2r(rank, rdims, x, F, kind, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
  fftw_execute(plan);
  // A new plan with the same parameters is fast to create it says, so we destroy this to
  // clean up memory
  fftw_destroy_plan(plan);

  // adjust scale to fit our setup
  for(int i = 0; i < siz; i++) F[i] *= isiz;

  // We need to adjust the first element in each dimension by 0.5
  // fftw's DCT isn't exactly equal to the Chebyshev transform
  int blocklen = 1;  // length of block
  int stride = 1; // distance between blocks
  for(int i = 0; i < rank; i++) {
    stride *= dims[i];
    for(int j = 0; j < siz; j += stride) {
      for(int s = 0; s < blocklen; s++) {
	F[j+s] *= 0.5;
      }
    }
    blocklen *= dims[i];
  }

#else
  double *src, *dest, *buf;
  double **mat;
  double beta=0;

  // Some work space
  // 
  buf = (double*) R_alloc(siz,sizeof(double));
  mat = (double**) R_alloc(rank,sizeof(double*));
  // Create the needed transform matrices
  // reuse same dimension matrices
  for(int j = 0; j < rank; j++) {
    // Is it there already?
    int N = dims[j];
    double *jmat = NULL;
    for(int k = 0; k < j; k++) {
      if(dims[k] == N) {
	jmat = mat[k];
	break;
      }
    }
    if(jmat != NULL) {
      mat[j] = jmat;
      continue;
    }
    jmat = mat[j] = (double*) R_alloc(N*N,sizeof(double));
    for(int k = 0; k < N; k++) {
      double *jkvec = &jmat[k*N];
      for(int i = 0; i < N; i++) {
	jkvec[i] = cos(M_PI*k*(i+0.5)/N);
      }
    }
  }
  // We will switch src and dest between F and buf
  // We should end up with dest=F, so if 
  // rank is odd we should start with src = F, dest=buf
  if(rank & 1 == 1) {
    src = F;
    dest = buf;
  } else {
    src = buf;
    dest = F;
  }
  memcpy(dest,x,siz*sizeof(double));
  for(int i = rank-1; i >= 0; i--) {
    // transformation of dimension i, put in front
    int N = dims[i];  // length of transform
    int stride = siz/N;
    double alpha = 2.0/N;  // Constant for transform
    // swap src and dest
    double *p = src;
    src = dest;
    dest = p;

    F77_NAME(dgemm)("TRANS","TRANS",&N,&stride,&N,&alpha,mat[i],&N,src,&stride,&beta,dest,&N);
    // Fix the first element in each vector
    for(int k = 0; k < siz; k+= N) dest[k] *= 0.5;
  }
#endif
}



static SEXP R_chebcoef(SEXP x) {

  SEXP dim;
  int rank;
  int *dims;
  size_t siz;
  SEXP resvec;
  int sdims;
  dim = getAttrib(x,R_DimSymbol);
  // If no dim-attribute, assume one-dimensional
  if(isNull(dim)) {
    sdims = LENGTH(x);
    dims = &sdims;
    rank = 1;

  } else {
    rank = LENGTH(dim);
    dims = INTEGER(dim);
  }
  siz = 1;
  for(int i = 0; i < rank; i++) siz *= dims[i];

  PROTECT(resvec = allocVector(REALSXP,siz));
  chebcoef(REAL(x),dims,&rank,REAL(resvec));
  setAttrib(resvec,R_DimSymbol,dim);
  UNPROTECT(1);

  return resvec;
}

static SEXP R_evalcheb(SEXP coef, SEXP x) {
  int rank = LENGTH(x);
  int dims[rank];
  size_t siz = 1;
  double *xvec[rank];
  double *cf = REAL(coef);
  SEXP resvec;
  for(int i = 0; i < rank; i ++) {
    xvec[i] = REAL(VECTOR_ELT(x,i));
    dims[i] = LENGTH(VECTOR_ELT(x,i));
    siz *= dims[i];
  }
  if(LENGTH(coef) != siz) error("coef length(%d) must match data length(%d)",LENGTH(coef),siz);

  PROTECT(resvec = allocVector(REALSXP,1));
  REAL(resvec)[0] = evalcheb(cf,xvec,dims,&rank);
  UNPROTECT(1);
  return resvec;
}

static SEXP R_hasfftw() {
  SEXP res = allocVector(LGLSXP,1);
#ifdef HAVE_FFTW
  LOGICAL(res)[0] = TRUE;
#else
  LOGICAL(res)[0] = FALSE;
#endif
  return res;
}

R_CallMethodDef callMethods[] = {
  {"evalcheb", (DL_FUNC) &R_evalcheb, 2},
  {"chebcoef", (DL_FUNC) &R_chebcoef, 1},
  {"hasfftw", (DL_FUNC) &R_hasfftw, 0},
  {NULL, NULL, 0}
};

void R_init_chebpol(DllInfo *info) {
  /* register our routines */
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
  //  R_RegisterCCallable("chebpol", "evalcheb", (DL_FUNC) &evalcheb);

}
void R_unload_chebpol(DllInfo *info) {
#ifdef HAVE_FFTW
  // Clean up fftw
  fftw_cleanup();
#endif
}
