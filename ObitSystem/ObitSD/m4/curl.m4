# autoconfig script for curl library
AC_DEFUN([AC_PATH_CURL], [
	AC_ARG_WITH(curl,
                    AC_HELP_STRING([--with-curl=DIR],
                                 [search for CURL in DIR/include and DIR/lib]),
                    [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir/include; then
      CURL_CFLAGS="$CURL_CFLAGS -I$dir/include/curl"
    fi
    if test -d $dir/lib; then
      CURL_LDFLAGS="$CURL_LDFLAGS -L$dir/lib"
    fi
  done[]])

        AC_ARG_WITH(curl-includes,
                    AC_HELP_STRING([--with-curl-includes=DIR],
	                           [search for CURL includes in DIR]),
	            [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir; then
      CURL_CFLAGS="$CURL_CFLAGS -I$dir"
    fi
  done[]])

ac_curl_saved_CFLAGS="$CFLAGS"
ac_curl_saved_LDFLAGS="$LDFLAGS"
ac_curl_saved_LIBS="$LIBS"
if ! test CURL_CFLAGS; then
  CURL_CFLAGS="`curl-config --cflags`"
fi
CFLAGS="$CFLAGS $CURL_CFLAGS"
CURL_LIBS="`curl-config --libs`"
LDFLAGS="$LDFLAGS $CURL_LDFLAGS"
LIBS="$CURL_LIBS"
ac_have_curl=no
ac_have_curlh=no
  	touch /tmp/dummy1_curl.h
        AC_CHECK_HEADERS([/tmp/dummy1_curl.h], [ac_have_curlh=yes], [ac_have_curlh=no],
			[#include "curl.h"])
	rm /tmp/dummy1_curl.h
 	if test $ac_have_curlh = yes; then
	        AC_SEARCH_LIBS(curl_version, [curl], [ac_have_curl=yes], [ac_have_curl=no])
	fi
# List of places to try
testdirs="$HOME/opt/curl $OBITINSTALL/other"
for dir in $testdirs; do
	if test $ac_have_curl = no; then
		if  test -f $dir/include/curl/curl.h; then
			CURL_CFLAGS="-I$dir/include/curl"
			CPPFLAGS="$ac_curl_saved_CPPFLAGS $CURL_CFLAGS"
			CURL_LDFLAGS="-L$dir/lib"
			LDFLAGS="$ac_curl_saved_LDFLAGS $CURL_LDFLAGS"
  			touch /tmp/dummy3_curl.h
	        	AC_CHECK_HEADERS(/tmp/dummy3_curl.h, [ac_have_curlh=yes], [ac_have_curlh=no],
				[#include "curl.h"])
			rm /tmp/dummy3_curl.h
			if test $ac_have_curlh = yes; then
				# Force check
				ac_cv_search_curl_version=" "
			        AC_SEARCH_LIBS(curl_version, [curl], [ac_have_curl=yes], [ac_have_curl=no])
			fi
			if test $ac_have_curl = yes ; then
				if test $ac_have_curlh = yes ; then
					break;
				fi
			fi
		fi
	fi
done[]
CURL_LIBS="$LIBS"
if test $ac_have_curl = no; then
	AC_MSG_WARN([cannot find CURL library])
fi
if test $ac_have_curlh = no; then
	AC_MSG_WARN([cannot find CURL headers])
	ac_have_curl=no
	CURL_CFLAGS=""
	CURL_LDFLAGS=""
	CURL_LIBS=""
fi
if test $ac_have_curl = yes; then
	AC_DEFINE(HAVE_CURL, 1, [Define to 1 if CURL is available.])
fi
CFLAGS="$ac_curl_saved_CFLAGS"
LDFLAGS="$ac_curl_saved_LDFLAGS"
LIBS="$ac_curl_saved_LIBS"
	AC_SUBST(CURL_CFLAGS)
	AC_SUBST(CURL_LDFLAGS)
	AC_SUBST(CURL_LIBS)
])
