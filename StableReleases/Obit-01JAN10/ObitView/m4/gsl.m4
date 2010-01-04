# autoconfig script for gsl library
AC_DEFUN([AC_PATH_GSL], [
	AC_ARG_WITH(gsl,
                    AC_HELP_STRING([--with-gsl=DIR],
                                 [search for GSL in DIR/include and DIR/lib]),
                    [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir/include; then
      GSL_CFLAGS="$GSL_CFLAGS -I$dir/include"
    fi
    if test -d $dir/lib; then
      GSL_LDFLAGS="$GSL_LDFLAGS -L$dir/lib"
    fi
  done[]])

        AC_ARG_WITH(gsl-includes,
                    AC_HELP_STRING([--with-gsl-includes=DIR],
	                           [search for GSL includes in DIR]),
	            [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir; then
      GSL_CFLAGS="$GSL_CFLAGS -I$dir"
    fi
  done[]])

ac_gsl_saved_CFLAGS="$CFLAGS"
ac_gsl_saved_LDFLAGS="$LDFLAGS"
ac_gsl_saved_LIBS="$LIBS"
CFLAGS="$CFLAGS $GSL_CFLAGS"
LDFLAGS="$LDFLAGS $GSL_LDFLAGS"
if ! test GSL_CFLAGS; then
    GSL_CFLAGS="`gsl-config --cflags`"
fi
GSL_LIBS="`gsl-config --libs`"
ac_have_gsl=no
ac_have_gslh=no
  	touch /tmp/dummy1_gsl.h
        AC_CHECK_HEADERS([/tmp/dummy1_gsl.h], [ac_have_gslh=yes], [ac_have_gslh=no],
			[#include "gsl/gsl_blas.h"])
	rm /tmp/dummy1_gsl.h
 	if test $ac_have_gslh = yes; then
	        AC_SEARCH_LIBS(cblas_isamax, [gslcblas], [ac_have_gsl=yes], [ac_have_gsl=no], [-lm])
        	AC_SEARCH_LIBS(gsl_fit_linear, [gsl], [ac_have_gsl=yes], [ac_have_gsl=no], [-lm])
	fi
# List of places to try
testdirs="$HOME/opt/gsl $OBITINSTALL/other"
for dir in $testdirs; do
	if test $ac_have_gsl = no; then
		if  test -f $dir/include/gsl/gsl_blas.h; then
			GSL_CFLAGS="-I$dir/include"
			CPPFLAGS="$ac_gsl_saved_CPPFLAGS $GSL_CFLAGS"
			GSL_LDFLAGS="-L$dir/lib"
			LDFLAGS="$ac_gsl_saved_LDFLAGS $GSL_LDFLAGS"
  			touch /tmp/dummy3_gsl.h
	        	AC_CHECK_HEADERS(/tmp/dummy3_gsl.h, [ac_have_gslh=yes], [ac_have_gslh=no],
				[#include "gsl/gsl_blas.h"])
			rm /tmp/dummy3_gsl.h
			if test $ac_have_gslh = yes; then
				# Force check
				ac_cv_search_cblas_isamax="  "
				ac_cv_search_gsl_fit_linear="  "
			        AC_SEARCH_LIBS(cblas_isamax, [gslcblas], [ac_have_gsl=yes], [ac_have_gsl=no], [-lm])
        			AC_SEARCH_LIBS(gsl_fit_linear, [gsl], [ac_have_gsl=yes], [ac_have_gsl=no], 
					[-lm -lgslcblas])
			fi
			if test $ac_have_gsl = yes ; then
				if test $ac_have_gslh = yes ; then
					break;
				fi
			fi
		fi
	fi
done[]
if test $ac_have_gslh = no; then
	AC_MSG_WARN([cannot find GSL headers])
	ac_have_gsl=no
fi
if test $ac_have_gsl = no; then
	AC_MSG_WARN([cannot find GSL library])
	GSL_CFLAGS=""
	GSL_LDFLAGS=""
	GSL_LIBS=""
fi
if test $ac_have_gsl = yes; then
	AC_DEFINE(HAVE_GSL, 1, [Define to 1 if GSL is available.])
fi
GSL_LIBS="$LIBS"
CFLAGS="$ac_gsl_saved_CFLAGS"
LDFLAGS="$ac_gsl_saved_LDFLAGS"
LIBS="$ac_gsl_saved_LIBS"
	AC_SUBST(GSL_CFLAGS)
	AC_SUBST(GSL_LDFLAGS)
	AC_SUBST(GSL_LIBS)
])
