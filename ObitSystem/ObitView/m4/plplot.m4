# autoconfig script for plplot library
AC_DEFUN([AC_PATH_PLPLOT], [

         AC_ARG_WITH(plplot,
                     AC_HELP_STRING([--with-plplot=DIR],
                                    [search for PLPLOT in DIR]),
                     [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir; then
 	PLPLOT_CPPFLAGS="$PLPLOT_CPPFLAGS -I$dir/include/plplot"
 	PLPLOT_CFLAGS="$PLPLOT_CFLAGS -I$dir/include/plplot"
	PLPLOT_LDFLAGS="$PLPLOT_LDFLAGS -L$dir/lib"
    fi
  done[]])
ac_plplot_saved_CPPFLAGS="$CPPFLAGS"
ac_plplot_saved_CFLAGS="$CFLAGS"
ac_plplot_saved_LDFLAGS="$LDFLAGS"
ac_plplot_saved_LIBS="$LIBS"
if test "x$PLPLOT_CFLAGS" = x; then
    PLPLOT_CFLAGS="`pkg-config plplot --cflags`"
fi
if test "x$PLPLOT_CPPFLAGS" = x; then
    PLPLOT_CPPFLAGS="`pkg-config plplot-c++ --cflags`"
fi
PLPLOT_LIBS="`pkg-config plplot --libs`"
CPPFLAGS="$CPPFLAGS $PLPLOT_CPPFLAGS"
CFLAGS="$CFLAGS $PLPLOT_CFLAGS"
LDFLAGS="$LDFLAGS $PLPLOT_LDFLAGS"
LIBS="$PLPLOT_LIBS"
ac_have_plplot=no
ac_have_plploth=no
# Standard location?
  	touch /tmp/dummy1_plplot.h
	AC_CHECK_HEADER(/tmp/dummy1_plplot.h, [ac_have_plploth=yes], [ac_have_plploth=no], 
	[#include "plplot.h"])
	rm /tmp/dummy1_plplot.h
 	if test $ac_have_plploth = yes; then
		AC_CHECK_LIB(plplot, c_plinit, [ac_have_plplot=yes], [ac_have_plplot=no])
	fi
# List of places to try
testdirs="$HOME/opt/plplot $OBITINSTALL/other"
for dir in $testdirs; do
	if test $ac_have_plplot = no; then
		if  test -f $dir/include/plplot/plplot.h; then
			#echo "Trying in $dir"
			PLPLOT_CFLAGS="-I$dir/include/plplot/"
			CPPFLAGS="$ac_plplot_saved_CPPFLAGS $PLPLOT_CFLAGS"
			PLPLOT_LDFLAGS="-L$dir/lib"
			LDFLAGS="$ac_plplot_saved_LDFLAGS $PLPLOT_LDFLAGS"
  			touch /tmp/dummy3_plplot.h
			AC_CHECK_HEADER(/tmp/dummy3_plplot.h, [ac_have_plploth=yes], [ac_have_plploth=no], 
				[#include "plplot.h"])
			rm /tmp/dummy3_plplot.h
			if test $ac_have_plploth = yes; then
				# Force check
				ac_cv_search_c_plinit="  "
				AC_SEARCH_LIBS(c_plinit, plplot, [ac_have_plplot=yes], [ac_have_plplot=no],
					[])
			fi
			if test $ac_have_plplot = yes ; then
				if test $ac_have_plploth = yes ; then
					break;
				fi
			fi
		fi
	fi
done[]
PLPLOT_LIBS="-lplplot $PLPLOT_LIBS"
if test $ac_have_plploth = no; then
	AC_MSG_WARN([cannot find PLPLOT headers])
	ac_have_plplot=no
	PLPLOT_CPPFLAGS=""
	PLPLOT_CFLAGS=""
fi
if test $ac_have_plplot = no; then
	AC_MSG_WARN([cannot find PLPLOT C library])
	PLPLOT_LIBS=""
	PLPLOT_LDFLAGS=""
fi
if test $ac_have_plplot = yes; then
	 AC_DEFINE(HAVE_PLPLOT, 1, [Define to 1 if you have PLPLOT.])
fi
CPPFLAGS="$ac_plplot_saved_CPPFLAGS"
CFLAGS="$ac_plplot_saved_CFLAGS"
LDFLAGS="$ac_plplot_saved_LDFLAGS"
LIBS="$ac_plplot_saved_LIBS"
	 AC_SUBST(PLPLOT_CPPFLAGS)
	 AC_SUBST(PLPLOT_CFLAGS)
	 AC_SUBST(PLPLOT_LDFLAGS)
	 AC_SUBST(PLPLOT_LIBS)
])
