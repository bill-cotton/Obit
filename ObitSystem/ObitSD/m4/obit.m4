AC_DEFUN([AC_PATH_OBIT], [
# Default root of Obit directory is $OBIT
	OBIT_DIR="$OBIT"
	if test -d $OBIT_DIR/include; then
	     OBIT_CPPFLAGS="$OBIT_CPPFLAGS -I$OBIT_DIR/include $GLIB_CFLAGS"
	fi
	if test -d $OBIT_DIR/lib/; then
	     OBIT_LDFLAGS="$OBIT_LDFLAGS -L$OBIT_DIR/lib/ $GLIB_LIBS"
	fi

	AC_ARG_WITH(obit,
                    AC_HELP_STRING([--with-obit=DIR],
                             [search for OBIT in DIR/include and DIR/lib]),
                    [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir/include; then
      OBIT_CPPFLAGS="$OBIT_CPPFLAGS -I$dir/include $GLIB_CFLAGS"
      OBIT_DIR="$dir"
      OBIT="$dir"
    fi
    if test -d $dir/lib/$ARCH; then
      OBIT_LDFLAGS="$OBIT_LDFLAGS -L$dir/lib/ $GLIB_LIBS"
    fi
  done[]])

        AC_ARG_WITH(obit-includes,
                    AC_HELP_STRING([--with-obit-includes=DIR],
	                           [search for OBIT includes in DIR]),
	            [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir; then
      OBIT_CPPFLAGS="$OBIT_CPPFLAGS -I$dir"
    fi
  done[]])

ac_obit_saved_CPPFLAGS="$CPPFLAGS"
ac_obit_saved_LDFLAGS="$LDFLAGS"
ac_obit_saved_LIBS="$LIBS"
CPPFLAGS="$CPPFLAGS $OBIT_CPPFLAGS"
LDFLAGS="$LDFLAGS $OBIT_LDFLAGS"
ac_have_obit=yes
        AC_CHECK_LIB(Obit, newObit, [], [ac_have_obit=no
	             AC_MSG_ERROR([cannot find OBIT library])]
			$GLIB_LIBS -lglib-2.0)
        AC_CHECK_HEADER(Obit.h, [], [ac_have_obit=no
	                AC_MSG_ERROR([cannot find OBIT headers])])
if test $ac_have_obit = yes; then
	AC_DEFINE(HAVE_OBIT, 1, [Define to 1 if you have OBIT.])
fi
OBIT_LIBS="$LIBS"
CPPFLAGS="$ac_obit_saved_CPPFLAGS"
LDFLAGS="$ac_obit_saved_LDFLAGS"
LIBS="$ac_obit_saved_LIBS"
	AC_SUBST(OBIT_CPPFLAGS)
	AC_SUBST(OBIT_LDFLAGS)
	AC_SUBST(OBIT_LIBS)
	AC_SUBST(OBIT_DIR)
	AC_SUBST(OBIT)
])
