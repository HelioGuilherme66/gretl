# AC_C_OPENMP
# -----------
# Check which options need to be passed to the C compiler to support OpenMP.
# Set the OPENMP_CFLAGS variable to these options.
# The options are necessary at compile time (so the #pragmas are understood)
# and at link time (so the appropriate library is linked with).
# This macro takes care to not produce redundant options if $CC $CFLAGS already
# supports OpenMP. It also is careful to not pass options to compilers that
# misinterpret them; for example, most compilers accept "-openmp" and create
# an output file called 'penmp' rather than activating OpenMP support.
# -----------
AC_DEFUN([AC_C_OPENMP],
[
  AC_MSG_CHECKING([whether to use OpenMP])
  AC_ARG_ENABLE(openmp,
    [AS_HELP_STRING([--enable-openmp], [use OpenMP if available [default=auto]])],
    [enable_openmp=$enableval]
  )
  ac_openmp_result=no
  OPENMP_CFLAGS=
  if test "$enable_openmp" != "no" ; then
    AC_CACHE_VAL([ac_cv_prog_cc_openmp], [
      ac_cv_prog_cc_openmp=unsupported
      AC_LINK_IFELSE([AC_LANG_SOURCE([
#ifndef _OPENMP
 choke me
#endif
#include <omp.h>
int main () { return omp_get_num_threads (); }
	])], [ac_cv_prog_cc_openmp="none needed"])
      if test "$ac_cv_prog_cc_openmp" = unsupported; then
	dnl Try these flags:
	dnl   GCC >= 4.2	   -fopenmp
	dnl   SunPRO C  	   -xopenmp
	dnl   Intel C		   -openmp
	dnl   SGI C, PGI C	   -mp
	dnl   Tru64 Compaq C	   -omp
	dnl   IBM C (AIX, Linux)   -qsmp=omp
	dnl   InteloneAPI          -qopenmp
        dnl   Apple                -Xclang -fopenmp
	for brand in apple clang GCC SunPRO Intel SGI/PGI Compaq IBM InteloneAPI; do
	  case $brand in
            apple)
              ac_conditional='defined __clang__ && defined __apple_build_version__'
              ac_option='-Xclang -fopenmp' ;;
	    clang)
	      ac_conditional='defined __clang__ && !(defined __INTEL_LLVM_COMPILER)'
	      ac_option='-fopenmp=libomp' ;;
	    GCC)
	      ac_conditional='defined __GNUC__'
	      ac_option='-fopenmp' ;;
	    SunPRO)
	      ac_conditional='defined __SUNPRO_C || defined __SUNPRO_CC'
	      ac_option='-xopenmp' ;;
	    Intel)
	      ac_conditional='defined __INTEL_COMPILER'
	      ac_option='-openmp' ;;
	    SGI/PGI)
	      ac_conditional='defined __sgi || defined __PGI || defined __PGIC__'
	      ac_option='-mp' ;;
	    Compaq)
	      ac_conditional='defined __DECC || defined __DECCXX'
	      ac_option='-omp' ;;
            IBM)
	      ac_conditional='defined __xlc__ || defined __xlC__'
	      ac_option='-qsmp=omp' ;;
        InteloneAPI)
          ac_conditional='defined __INTEL_LLVM_COMPILER'
          ac_option='-qopenmp' ;;
	  esac
	  if test $brand = GCC; then
	    if test $CC = icx || test $CC = icpx; then
	      ac_openmp_result=no
	    elif test "$GCC" = yes; then
	      ac_openmp_result=yes
	    else
	      ac_openmp_result=no
	    fi
	  else
	    AC_EGREP_CPP([Brand], [
	      #if $ac_conditional
	       Brand
	      #endif
	      ], [ac_openmp_result=yes], [ac_openmp_result=no])
	  fi
	  if test $ac_openmp_result = yes; then
            ac_save_CFLAGS=$CFLAGS
            ac_save_LIBS=$LIBS
            CFLAGS="$CFLAGS $ac_option"
            LIBS="$LIBS $OMP_LIB"
            AC_LINK_IFELSE([AC_LANG_SOURCE([
#ifndef _OPENMP
 choke me
#endif
#include <omp.h>
int main () { return omp_get_num_threads (); }
	      ])], [ac_cv_prog_cc_openmp=$ac_option])
            CFLAGS=$ac_save_CFLAGS
            LIBS=$ac_save_LIBS
	    break
	  fi
	done
      fi
      ])
    dnl AC_MSG_RESULT([$ac_cv_prog_cc_openmp])
    case $ac_cv_prog_cc_openmp in
      "none needed")
        ac_openmp_result=yes
	OPENMP_CFLAGS= ;;
      unsupported)
        ac_openmp_result=no
        OPENMP_CFLAGS= ;;
      *)
	OPENMP_CFLAGS=$ac_cv_prog_cc_openmp ;;
    esac
  fi
  AC_MSG_RESULT([$ac_openmp_result])
  AC_SUBST([OPENMP_CFLAGS])
])

