# autoconfig script for cfitsio library
AC_DEFUN([AC_PATH_CFITSIO], [
	AC_ARG_WITH(cfitsio,
                    AC_HELP_STRING([--with-cfitsio=DIR],
                             [search for CFITSIO in DIR/include and DIR/lib]),
                    [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir/include; then
      CFITSIO_CPPFLAGS="$CFITSIO_CPPFLAGS -I$dir/include"
    fi
    if test -d $dir/lib; then
      CFITSIO_LDFLAGS="$CFITSIO_LDFLAGS -L$dir/lib"
    fi
  done[]])

        AC_ARG_WITH(cfitsio-includes,
                    AC_HELP_STRING([--with-cfitsio-includes=DIR],
	                           [search for CFITSIO includes in DIR]),
	            [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir; then
      CFITSIO_CPPFLAGS="$CFITSIO_CPPFLAGS -I$dir"
    fi
  done[]])

ac_cfitsio_saved_CPPFLAGS="$CPPFLAGS"
ac_cfitsio_saved_LDFLAGS="$LDFLAGS"
ac_cfitsio_saved_LIBS="$LIBS"
CPPFLAGS="$CPPFLAGS $CFITSIO_CPPFLAGS"
LDFLAGS="$LDFLAGS $CFITSIO_LDFLAGS"
# Check in standard/suggested place
ac_have_cfitsioh=yes
ac_have_cfitsio=yes
	touch /tmp/dummy1_fitsio.h
        AC_CHECK_HEADER([/tmp/dummy1_fitsio.h], [], [ac_have_cfitsioh=no], [#include  <fitsio.h>])
	rm /tmp/dummy1_fitsio.h
	if test $ac_have_cfitsioh = no; then
	# Try where cfitsio likes to put includes
		if  test -f /usr/include/cfitsio/fitsio.h; then
			CFITSIO_CPPFLAGS="-I/usr/include/cfitsio"
			CPPFLAGS="$ac_cfitsio_saved_CPPFLAGS $CFITSIO_CPPFLAGS"
  			touch /tmp/dummy2_fitsio.h
      			AC_CHECK_HEADERS([/tmp/dummy2_fitsio.h], [ac_have_cfitsioh=yes], [ac_have_cfitsioh=no
	                	AC_MSG_WARN([cannot find CFITSIO headers in /usr/include/cfitsio])],
				[#include  <fitsio.h>])
			rm /tmp/dummy2_fitsio.h
		fi
	fi
 	if test $ac_have_cfitsioh = yes; then
        	AC_CHECK_LIB([cfitsio], ffinit, [], [ac_have_cfitsio=no
	             AC_MSG_WARN([cannot find CFITSIO library in default])],[-lm])
	fi
# Hack to get around caching
testdirs="$HOME/opt/cfitsio $OBITINSTALL/other"
for dir in $testdirs; do
	if test $ac_have_cfitsio = no; then
		if  test -d $dir; then
			CFITSIO_CPPFLAGS="-I$dir/include"
			CPPFLAGS="$ac_cfitsio_saved_CPPFLAGS $CFITSIO_CPPFLAGS"
			CFITSIO_LDFLAGS="-L$dir/lib"
			LDFLAGS="$ac_cfitsio_saved_LDFLAGS $CFITSIO_LDFLAGS"
  			touch /tmp/dummy3_fitsio.h
	        	AC_CHECK_HEADERS(/tmp/dummy3_fitsio.h, [], [ac_have_cfitsio=no
		        	AC_MSG_WARN([cannot find CFITSIO headers in $dir])],
				[#include  "fitsio.h"])
			rm /tmp/dummy3_fitsio.h
 			if test $ac_have_cfitsioh = yes; then
				# Force check
				ac_cv_search_ffinit="  "
		        	AC_CHECK_LIB([cfitsio], ffinit, [ac_have_cfitsio=yes], [ac_have_cfitsio=no
				     AC_MSG_WARN([cannot find CFITSIO library in $dir])],[-lm])
			fi
			if test $ac_have_fftw = yes ; then
				if test $ac_have_fftwh = yes ; then
					break;
				fi
			fi
		fi
	fi
done[]
CFITSIO_LIBS="$LIBS"
if test $ac_have_cfitsioh = no; then
	ac_have_cfitsio=no
fi
if test $ac_have_cfitsio = no; then
	ac_have_cfitsio=no
	CFITSIO_CPPFLAGS=""
	CFITSIO_LDFLAGS=""
	CFITSIO_LIBS=""
fi
if test $ac_have_cfitsio = yes; then
	AC_DEFINE(HAVE_CFITSIO, 1, [Define to 1 if you have CFITSIO.])
fi
CPPFLAGS="$ac_cfitsio_saved_CPPFLAGS"
LDFLAGS="$ac_cfitsio_saved_LDFLAGS"
LIBS="$ac_cfitsio_saved_LIBS"
	AC_SUBST(CFITSIO_CPPFLAGS)
	AC_SUBST(CFITSIO_LDFLAGS)
	AC_SUBST(CFITSIO_LIBS)
])
