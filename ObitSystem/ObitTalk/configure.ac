# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)

AC_INIT(ObitTalk, 0.1, bcotton@nrao.edu)
AM_INIT_AUTOMAKE(1.6)
AC_CONFIG_SRCDIR(python/ObitTalk.py)

## --------- ##
## Scripts.  ##
## --------- ##

AC_CONFIG_FILES(bin/Makefile)

## -------- ##
## Python.  ##
## -------- ##

AM_PATH_PYTHON(2.2)
AC_CONFIG_FILES([python/Makefile python/Proxy/Makefile
                 python/Wizardry/Makefile test/Makefile])

## ------ ##
## Obit.  ##
## ------ ##

OBIT_PYTHONPATH=
AC_ARG_WITH(obit,
            AC_HELP_STRING([--with-obit=DIR],
                           [search for Obit in DIR]),
            [if test -f "$withval/python/Obit.py"; then
    OBIT_PYTHONPATH="$withval/python"; export ac_obit_pythonpath
  fi[]])
AC_CACHE_CHECK([for Obit], [ac_cv_obit],
               [
cat <<_ACEOF >conftest.py
try:
    import OSystem, OErr
except:
    print "missing"
else:
    print "found"
_ACEOF
PYTHONPATH=$PYTHONPATH:$OBIT_PYTHONPATH; export PYTHONPATH
ac_cv_obit=`$PYTHON conftest.py 2>/dev/null`
])
if test "$ac_cv_obit" != found; then
  AC_MSG_ERROR([Obit is not available])
fi
AC_SUBST(OBIT_PYTHONPATH)

## ------------ ##
## Conclusion.  ##
## ------------ ##

AC_CONFIG_FILES([Makefile doc/Makefile])
AC_OUTPUT
