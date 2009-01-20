# autoconf script for pgplot library
AC_DEFUN([AM_PATH_PGPLOT], [
         AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])dnl
         AC_REQUIRE([AC_PATH_XTRA])dnl

         AC_ARG_WITH(pgplot,
                     AC_HELP_STRING([--with-pgplot=DIR],
                                    [search for PGPLOT in DIR]),
                     [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir; then
      PGPLOT_CPPFLAGS="$PGPLOT_CPPFLAGS -I$dir"
      PGPLOT_LDFLAGS="$PGPLOT_LDFLAGS -L$dir"
    fi
  done[]])
ac_pgplot_saved_CPPFLAGS="$CPPFLAGS"
ac_pgplot_saved_CFLAGS="$CFLAGS"
ac_pgplot_saved_LDFLAGS="$LDFLAGS"
ac_pgplot_saved_LIBS="$LIBS"
PGPLOT_CFLAGS="$PGPLOT_CFLAGS $X_CFLAGS"
PGPLOT_LIBS="$PGPLOT_LIBS $X_LIBS $X_PRE_LIBS -lX11 $X_EXTRA_LIBS"
CPPFLAGS="$CPPFLAGS $PGPLOT_CPPFLAGS"
CFLAGS="$CFLAGS $PGPLOT_CFLAGS"
LDFLAGS="$LDFLAGS $PGPLOT_LDFLAGS"
LIBS="$PGPLOT_LIBS $FLIBS"
ac_have_pgplot=no
ac_have_pgploth=no
  	touch /tmp/dummy1_pgplot.h
	AC_CHECK_HEADER(/tmp/dummy1_pgplot.h, [ac_have_pgploth=yes], [ac_have_pgploth=no], 
	[#include "cpgplot.h"])
	rm /tmp/dummy1_pgplot.h
 	if test $ac_have_pgploth = yes; then
		AC_CHECK_LIB(pgplot, main, [ac_have_pgplot=yes], [ac_have_pgplot=no])
		AC_CHECK_LIB(cpgplot, cpgenv, [ac_have_pgplot=yes], [ac_have_pgplot=no], [-lpgplot])
	fi
# List of places to try
testdirs="$HOME/opt/pgplot $OBITINSTALL/other"
for dir in $testdirs; do
	if test $ac_have_pgplot = no; then
		if  test -f $dir/include/cpgplot.h; then
			PGPLOT_CFLAGS="-I$dir/include"
			CPPFLAGS="$ac_pgplot_saved_CPPFLAGS $PGPLOT_CFLAGS"
			PGPLOT_LDFLAGS="-L$dir/lib"
			LDFLAGS="$ac_pgplot_saved_LDFLAGS $PGPLOT_LDFLAGS"
  			touch /tmp/dummy3_pgplot.h
			AC_CHECK_HEADER(/tmp/dummy3_pgplot.h, [ac_have_pgploth=yes], [ac_have_pgploth=no], 
				[#include "cpgplot.h"])
			rm /tmp/dummy3_pgplot.h
			if test $ac_have_pgploth = yes; then
				# Force check
				ac_cv_search_main="  "
				ac_cv_search_cpgenv="  "
				AC_CHECK_LIB(pgplot, main, [ac_have_pgplot=yes], [ac_have_pgplot=no])
				AC_CHECK_LIB(cpgplot, cpgenv, [ac_have_pgplot=yes], [ac_have_pgplot=no], [-lpgplot])
			fi
			if test $ac_have_pgplot = yes ; then
				if test $ac_have_pgploth = yes ; then
					break;
				fi
			fi
		fi
	fi
done[]
PGPLOT_LIBS="-lcpgplot -lpgplot $PGPLOT_LIBS"
if test $ac_have_pgploth = no; then
	AC_MSG_WARN([cannot find PGPLOT headers])
	ac_have_pgplot=no
	PGPLOT_CPPFLAGS=""
	PGPLOT_CFLAGS=""
fi
if test $ac_have_pgplot = no; then
	AC_MSG_WARN([cannot find PGPLOT C library])
	PGPLOT_LIBS=""
	PGPLOT_LDFLAGS=""
fi
# Use plplot in preference to pgplot
if test $ac_have_plplot = yes; then
	AC_MSG_WARN([Using PLPLOT C library in preference to PGPLOT])
	ac_have_pgplot=no
	PGPLOT_CPPFLAGS=""
	PGPLOT_CFLAGS=""
	PGPLOT_LIBS=""
	PGPLOT_LDFLAGS=""
fi
if test $ac_have_pgplot = yes; then
	 AC_DEFINE(HAVE_PGPLOT, 1, [Define to 1 if you have PGPLOT.])
fi
CPPFLAGS="$ac_pgplot_saved_CPPFLAGS"
CFLAGS="$ac_pgplot_saved_CFLAGS"
LDFLAGS="$ac_pgplot_saved_LDFLAGS"
LIBS="$ac_pgplot_saved_LIBS"
	 AC_SUBST(PGPLOT_CPPFLAGS)
	 AC_SUBST(PGPLOT_CFLAGS)
	 AC_SUBST(PGPLOT_LDFLAGS)
	 AC_SUBST(PGPLOT_LIBS)
])
