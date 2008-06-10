# autoconfig script for fftw library
AC_DEFUN([AC_PATH_FFTW], [
	AC_ARG_WITH(fftw,
                    AC_HELP_STRING([--with-fftw=DIR],
                                 [search for FFTW in DIR/include and DIR/lib]),
                    [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir/include; then
      FFTW_CPPFLAGS="$FFTW_CPPFLAGS -I$dir/include"
    fi
    if test -d $dir/lib; then
      FFTW_LDFLAGS="$FFTW_LDFLAGS -L$dir/lib"
    fi
  done[]])

        AC_ARG_WITH(fftw-includes,
                    AC_HELP_STRING([--with-fftw-includes=DIR],
	                           [search for FFTW includes in DIR]),
	            [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir; then
      FFTW_CPPFLAGS="$FFTW_CPPFLAGS -I$dir"
    fi
  done[]])

ac_fftw_saved_CPPFLAGS="$CPPFLAGS"
ac_fftw_saved_LDFLAGS="$LDFLAGS"
ac_fftw_saved_LIBS="$LIBS"
CPPFLAGS="$CPPFLAGS $FFTW_CPPFLAGS"
LDFLAGS="$LDFLAGS $FFTW_LDFLAGS"
ac_have_fftw=no
  	touch /tmp/dummy1_fftw.h
        AC_CHECK_HEADERS([/tmp/dummy1_fftw.h], [ac_have_fftwh=yes], [ac_have_fftwh=no],
			[#include "fftw.h"
#include "rfftw.h"])
	rm /tmp/dummy1_fftw.h
	if test $ac_have_fftwh = yes; then
        	AC_CHECK_LIB([fftw], fftw, [ac_have_fftw=yes FFTW_LIBS="$LIBS -lfftw"], [ac_have_fftw=no],
                       [-lm ])
        	AC_CHECK_LIB([rfftw], rfftw, [ac_have_fftw=yes FFTW_LIBS="$FFTW_LIBS -lrfftw"], [ac_have_fftw=no],
                       [-lm -lfftw])
	fi
# List of places to try
testdirs="$HOME/opt/fftw $OBITINSTALL/other"
for dir in $testdirs; do
	if test $ac_have_fftw = no; then
		if  test -f $dir/include/fftw.h; then
			FFTW_CPPFLAGS="-I$dir/include"
			CPPFLAGS="$ac_fftw_saved_CPPFLAGS $FFTW_CPPFLAGS"
			FFTW_LDFLAGS="-L$dir/lib"
			LDFLAGS="$ac_fftw_saved_LDFLAGS $FFTW_LDFLAGS"
  			touch /tmp/dummy3_fftw.h
	        	AC_CHECK_HEADERS(/tmp/dummy3_fftw.h, [ac_have_fftwh=yes], [ac_have_fftwh=no],
				[#include "fftw.h"
#include "rfftw.h"])
			rm /tmp/dummy3_fftw.h
			if test $ac_have_fftwh = yes; then
				# Force check
				ac_cv_search_fftw="  "
				ac_cv_search_rfftw="  "
	        		AC_CHECK_LIB([fftw], fftw, [ac_have_fftw=yes FFTW_LIBS="$LIBS -lfftw"], [ac_have_fftw=no],
        	        		[-lm ])
		        	AC_CHECK_LIB([rfftw], rfftw, [ac_have_fftw=yes FFTW_LIBS="$FFTW_LIBS -lrfftw"], [ac_have_fftw=no],
        		               [-lm -lfftw])
			fi
			if test $ac_have_fftw = yes ; then
				if test $ac_have_fftwh = yes ; then
					break;
				fi
			fi
		fi
	fi
done[]
if test $ac_have_fftwh = no; then
	AC_MSG_WARN([cannot find FFTW headers])
	ac_have_fftw=no
fi
if test $ac_have_fftw = no; then
	AC_MSG_WARN([cannot find FFTW library])
	FFTW_CPPFLAGS=""
	FFTW_LDFLAGS=""
	FFTW_LIBS=""
fi
# Use FFTW3 in preference to FFTW2
if test $ac_have_fftw3 = yes; then
	AC_MSG_WARN([Using FFTW 3 library in preference to FFTW2])
	ac_have_fftw=no
	FFTW_CPPFLAGS=""
	FFTW_CFLAGS=""
	FFTW_LIBS=""
	FFTW_LDFLAGS=""
fi
if test $ac_have_fftw = yes; then
	AC_DEFINE(HAVE_FFTW, 1, [Define to 1 if FFTW 2 is available.])
fi
CPPFLAGS="$ac_fftw_saved_CPPFLAGS"
LDFLAGS="$ac_fftw_saved_LDFLAGS"
LIBS="$ac_fftw_saved_LIBS"
	AC_SUBST(FFTW_CPPFLAGS)
	AC_SUBST(FFTW_LDFLAGS)
	AC_SUBST(FFTW_LIBS)
])
