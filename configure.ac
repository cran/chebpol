## 
## This file is part of chebpol
##
AC_PREREQ(2.61)

## ---[ VERSION ]-------------------------------------------------------------

AC_INIT()
AC_CONFIG_AUX_DIR(src/build/autoconf)

AC_ARG_ENABLE([fftw],
       [AS_HELP_STRING([--disable-fftw],[do not use fftw])],
       [use_fftw=$enableval],[use_fftw=yes])


# some R 
if test -n "${R_HOME}"; then
  CC=`${R_HOME}/bin/R CMD config CC`
  CFLAGS=`${R_HOME}/bin/R CMD config CFLAGS`
  CPPFLAGS=`${R_HOME}/bin/R CMD config CPPFLAGS`
fi

AC_PROG_CC

## --- FFTW 3 ----------------------------------------------------------------
AS_IF([test "x$use_fftw" = "xyes"], [
  PKG_CHECK_MODULES([FFTW], fftw3,[FFTW=-DHAVE_FFTW],[echo '*** FFTW not found. Compiling with matrix method'])
])
## --- Output ----------------------------------------------------------------
AC_SUBST(FFTW)
AC_SUBST(FFTW_CFLAGS)
AC_SUBST(FFTW_LIBS)

AC_SUBST(CPPFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(LIBS)

AC_OUTPUT(src/Makevars)
