# autoconfig script for fftw 3 library  float version
AC_DEFUN([AC_PATH_FFTW3], [
	AC_ARG_WITH(fftw3,
                    AC_HELP_STRING([--with-fftw3=DIR],
                                 [search for FFTW3 in DIR/include and DIR/lib]),
                    [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir/include; then
      FFTW3_CPPFLAGS="$FFTW3_CPPFLAGS -I$dir/include"
    fi
    if test -d $dir/lib; then
      FFTW3_LDFLAGS="$FFTW3_LDFLAGS -L$dir/lib"
    fi
  done[]])

        AC_ARG_WITH(fftw3-includes,
                    AC_HELP_STRING([--with-fftw3-includes=DIR],
	                           [search for FFTW3 includes in DIR]),
	            [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir; then
      FFTW3_CPPFLAGS="$FFTW3_CPPFLAGS -I$dir"
    fi
  done[]])

ac_fftw3_saved_CPPFLAGS="$CPPFLAGS"
ac_fftw3_saved_LDFLAGS="$LDFLAGS"
ac_fftw3_saved_LIBS="$LIBS"
CPPFLAGS="$CPPFLAGS $FFTW3_CPPFLAGS"
LDFLAGS="$LDFLAGS $FFTW3_LDFLAGS"
ac_have_fftw3=no
  	touch /tmp/dummy1_fftw3.h
        AC_CHECK_HEADERS([/tmp/dummy1_fftw3.h], [ac_have_fftw3h=yes], [ac_have_fftw3h=no],
			[#include "fftw3.h"])
	rm /tmp/dummy1_fftw3.h
	if test $ac_have_fftw3h = yes; then
        	AC_CHECK_LIB([fftw3f], fftwf_execute, [ac_have_fftw3=yes FFTW3_LIBS="$LIBS -lfftw3f"], [ac_have_fftw3=no],
                       [-lm ])
	fi
# List of places to try
testdirs="$HOME/opt/fftw3 $OBITINSTALL/other"
for dir in $testdirs; do
	if test $ac_have_fftw3 = no; then
		if  test -f $dir/include/fftw3.h; then
			FFTW3_CPPFLAGS="-I$dir/include"
			CPPFLAGS="$ac_fftw3_saved_CPPFLAGS $FFTW3_CPPFLAGS"
			FFTW3_LDFLAGS="-L$dir/lib"
			LDFLAGS="$ac_fftw3_saved_LDFLAGS $FFTW3_LDFLAGS"
  			touch /tmp/dummy3_fftw3.h
	        	AC_CHECK_HEADERS(/tmp/dummy3_fftw3.h, [ac_have_fftw3h=yes], [ac_have_fftw3h=no],
				[#include "fftw3.h"])
			rm /tmp/dummy3_fftw3.h
			if test $ac_have_fftw3h = yes; then
				# Force check
				ac_cv_search_fftw3="  "
	        		AC_CHECK_LIB([fftw3f], fftwf_execute, [ac_have_fftw3=yes FFTW3_LIBS="$LIBS -lfftw3"], [ac_have_fftw3=no],
        	        		[-lm ])
			fi
			if test $ac_have_fftw3 = yes ; then
				if test $ac_have_fftw3h = yes ; then
					break;
				fi
			fi
		fi
	fi
done[]
if test $ac_have_fftw3h = no; then
	AC_MSG_WARN([cannot find FFTW3 headers])
	ac_have_fftw3=no
fi
if test $ac_have_fftw3 = no; then
	AC_MSG_WARN([cannot find FFTW3 library])
	FFTW3_CPPFLAGS=""
	FFTW3_LDFLAGS=""
	FFTW3_LIBS=""
fi
if test $ac_have_fftw3 = yes; then
	AC_DEFINE(HAVE_FFTW3, 1, [Define to 1 if FFTW3 is available.])
fi
CPPFLAGS="$ac_fftw3_saved_CPPFLAGS"
LDFLAGS="$ac_fftw3_saved_LDFLAGS"
LIBS="$ac_fftw3_saved_LIBS"
	AC_SUBST(FFTW3_CPPFLAGS)
	AC_SUBST(FFTW3_LDFLAGS)
	AC_SUBST(FFTW3_LIBS)
])
