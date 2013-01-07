# autoconfig script for python2.5 library
# If explicitly told - look for python2.5
AC_DEFUN([AC_PATH_PYTHON2_5], [
	AC_ARG_WITH(python,
                    AC_HELP_STRING([--with-python=DIR],
                             [search for PYTHON in DIR/include and DIR/lib]),
                    [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir/include/; then
      PYTHON_CPPFLAGS="$PYTHON_CPPFLAGS -I$dir/include/"
    fi
    if test -d $dir/lib; then
      PYTHON_LDFLAGS="$PYTHON_LDFLAGS -L$dir/lib"
    fi
  done[]])

        AC_ARG_WITH(python-includes,
                    AC_HELP_STRING([--with-python-includes=DIR],
	                           [search for PYTHON includes in DIR]),
	            [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir; then
      PYTHON_CPPFLAGS="$PYTHON_CPPFLAGS -I$dir"
    fi
  done[]])

# If not given, look in standard places
# Includes
if test "x$PYTHON_CPPFLAGS" = x; then
    if test "x$PYTHON" = x; then
        AC_PATH_PROG(PYTHON, python,, `pwd`/../../bin$PATH_SEPARATOR$PATH)
        # old PYTHON=`pwd`/../../bin/python
    fi
cat <<_ACEOF >conftest.py
import distutils.sysconfig
print distutils.sysconfig.get_python_inc()
_ACEOF
    ac_python_inc=`$PYTHON conftest.py`
    # Use python installed in other to 
    PYTHON_CPPFLAGS="-I$ac_python_inc"
fi

# Python libs
if test "x$PYTHON_LD_FLAGS" = x; then
    if test "x$PYTHON" = x; then
        AC_PATH_PROG(PYTHON, python,, `pwd`/../../bin$PATH_SEPARATOR$PATH)
        # old PYTHON=`pwd`/../../bin/python
    fi
cat <<_ACEOF >conftest.py
import distutils.sysconfig
print distutils.sysconfig.get_python_lib()
_ACEOF
    ac_python_lib=`$PYTHON conftest.py`
    PYTHON_LDFLAGS="-L$ac_python_lib"
fi

# Found it?
ac_have_pythonh=yes
if test "x$PYTHON_CPPFLAGS" = x; then
  AC_MSG_WARN([cannot find PYTHON headers])
  ac_have_pythonh=no
fi
ac_have_python=yes
if test "x$PYTHON_LDFLAGS" = x; then
  AC_MSG_WARN([cannot find PYTHON library])
  ac_have_python=no
fi

ac_python_saved_CPPFLAGS="$CPPFLAGS"
ac_python_saved_LDFLAGS="$LDFLAGS"
ac_python_saved_LIBS="$LIBS"
CPPFLAGS="$CPPFLAGS $PYTHON_CPPFLAGS"
LDFLAGS="$LDFLAGS $PYTHON_LDFLAGS"
PYTHON_LIBS="$LIBS"

if test $ac_have_pythonh = no; then
	ac_have_python=no
fi
if test $ac_have_python = no; then
	ac_have_python=no
	PYTHON_CPPFLAGS=""
	PYTHON_LDFLAGS=""
	PYTHON_LIBS=""
fi
if test $ac_have_python = yes; then
	AC_DEFINE(HAVE_PYTHON, 1, [Define to 1 if you have PYTHON.])
fi
CPPFLAGS="$ac_python_saved_CPPFLAGS"
LDFLAGS="$ac_python_saved_LDFLAGS"
LIBS="$ac_python_saved_LIBS"
	AC_SUBST(PYTHON_CPPFLAGS)
	AC_SUBST(PYTHON_LDFLAGS)
	AC_SUBST(PYTHON_LIBS)
])
