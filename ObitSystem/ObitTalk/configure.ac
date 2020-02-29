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

# need better replacement AM_PATH_PYTHON(2.7)


        if test -n "$PYTHON"; then
      # If the user set $PYTHON, use it and don't search something else.
      echo "$as_me:$LINENO: checking whether $PYTHON version >= 2.7" >&5
echo $ECHO_N "checking whether $PYTHON version >= 2.7... $ECHO_C" >&6
      prog="import sys
sys.exit(sys.version_info < (2, 7, 0))"
  if { echo "$as_me:$LINENO: $PYTHON -c "$prog"" >&5
   ($PYTHON -c "$prog") >&5 2>&5
   ac_status=$?
   echo "$as_me:$LINENO: \$? = $ac_status" >&5
   (exit $ac_status); }; then
  echo "$as_me:$LINENO: result: yes" >&5
echo "${ECHO_T}yes" >&6
else
  { { echo "$as_me:$LINENO: error: too old" >&5
echo "$as_me: error: too old" >&2;}
   { (exit 1); exit 1; }; }
fi

      am_display_PYTHON=$PYTHON
    else
      # Otherwise, try each interpreter until we find one that satisfies
      # VERSION.
      echo "$as_me:$LINENO: checking for a Python interpreter with version >= 2.7" >&5
echo $ECHO_N "checking for a Python interpreter with version >= 2.7... $ECHO_C" >&6
if test "${am_cv_pathless_PYTHON+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else

	for am_cv_pathless_PYTHON in python python2 python3 none; do
	  test "$am_cv_pathless_PYTHON" = none && break
	  prog="import sys
sys.exit(sys.version_info < (2, 7, 0))"
  if { echo "$as_me:$LINENO: $am_cv_pathless_PYTHON -c "$prog"" >&5
   ($am_cv_pathless_PYTHON -c "$prog") >&5 2>&5
   ac_status=$?
   echo "$as_me:$LINENO: \$? = $ac_status" >&5
   (exit $ac_status); }; then
  break
fi

	done
fi
echo "$as_me:$LINENO: result: $am_cv_pathless_PYTHON" >&5
echo "${ECHO_T}$am_cv_pathless_PYTHON" >&6
      # Set $PYTHON to the absolute path of $am_cv_pathless_PYTHON.
      if test "$am_cv_pathless_PYTHON" = none; then
	PYTHON=:
      else
        # Extract the first word of "$am_cv_pathless_PYTHON", so it can be a program name with args.
set dummy $am_cv_pathless_PYTHON; ac_word=$2
echo "$as_me:$LINENO: checking for $ac_word" >&5
echo $ECHO_N "checking for $ac_word... $ECHO_C" >&6
if test "${ac_cv_path_PYTHON+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  case $PYTHON in
  [\\/]* | ?:[\\/]*)
  ac_cv_path_PYTHON="$PYTHON" # Let the user override the test with a path.
  ;;
  *)
  as_save_IFS=$IFS; IFS=$PATH_SEPARATOR
for as_dir in $PATH
do
  IFS=$as_save_IFS
  test -z "$as_dir" && as_dir=.
  for ac_exec_ext in '' $ac_executable_extensions; do
  if $as_executable_p "$as_dir/$ac_word$ac_exec_ext"; then
    ac_cv_path_PYTHON="$as_dir/$ac_word$ac_exec_ext"
    echo "$as_me:$LINENO: found $as_dir/$ac_word$ac_exec_ext" >&5
    break 2
  fi
done
done

  ;;
esac
fi
PYTHON=$ac_cv_path_PYTHON

if test -n "$PYTHON"; then
  echo "$as_me:$LINENO: result: $PYTHON" >&5
echo "${ECHO_T}$PYTHON" >&6
else
  echo "$as_me:$LINENO: result: no" >&5
echo "${ECHO_T}no" >&6
fi

      fi
      am_display_PYTHON=$am_cv_pathless_PYTHON
    fi


  if test "$PYTHON" = :; then
      { { echo "$as_me:$LINENO: error: no suitable Python interpreter found" >&5
echo "$as_me: error: no suitable Python interpreter found" >&2;}
   { (exit 1); exit 1; }; }
  else


  echo "$as_me:$LINENO: checking for $am_display_PYTHON version" >&5
echo $ECHO_N "checking for $am_display_PYTHON version... $ECHO_C" >&6
if test "${am_cv_python_version+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  am_cv_python_version=`$PYTHON -c "import sys; print sys.version[:3]"`
fi
echo "$as_me:$LINENO: result: $am_cv_python_version" >&5
echo "${ECHO_T}$am_cv_python_version" >&6
  PYTHON_VERSION=$am_cv_python_version



  PYTHON_PREFIX='${prefix}'

  PYTHON_EXEC_PREFIX='${exec_prefix}'



  echo "$as_me:$LINENO: checking for $am_display_PYTHON platform" >&5
echo $ECHO_N "checking for $am_display_PYTHON platform... $ECHO_C" >&6
if test "${am_cv_python_platform+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  am_cv_python_platform=`$PYTHON -c "import sys; print sys.platform"`
fi
echo "$as_me:$LINENO: result: $am_cv_python_platform" >&5
echo "${ECHO_T}$am_cv_python_platform" >&6
  PYTHON_PLATFORM=$am_cv_python_platform




                echo "$as_me:$LINENO: checking for $am_display_PYTHON script directory" >&5
echo $ECHO_N "checking for $am_display_PYTHON script directory... $ECHO_C" >&6
if test "${am_cv_python_pythondir+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  am_cv_python_pythondir=`$PYTHON -c "from distutils import sysconfig; print sysconfig.get_python_lib(0,0,prefix='$PYTHON_PREFIX')" 2>/dev/null ||
     echo "$PYTHON_PREFIX/lib/python$PYTHON_VERSION/site-packages"`
fi
echo "$as_me:$LINENO: result: $am_cv_python_pythondir" >&5
echo "${ECHO_T}$am_cv_python_pythondir" >&6
  pythondir=$am_cv_python_pythondir



  pkgpythondir=\${pythondir}/$PACKAGE


            echo "$as_me:$LINENO: checking for $am_display_PYTHON extension module directory" >&5
echo $ECHO_N "checking for $am_display_PYTHON extension module directory... $ECHO_C" >&6
if test "${am_cv_python_pyexecdir+set}" = set; then
  echo $ECHO_N "(cached) $ECHO_C" >&6
else
  am_cv_python_pyexecdir=`$PYTHON -c "from distutils import sysconfig; print sysconfig.get_python_lib(1,0,prefix='$PYTHON_EXEC_PREFIX')" 2>/dev/null ||
     echo "${PYTHON_EXEC_PREFIX}/lib/python${PYTHON_VERSION}/site-packages"`
fi
echo "$as_me:$LINENO: result: $am_cv_python_pyexecdir" >&5
echo "${ECHO_T}$am_cv_python_pyexecdir" >&6
  pyexecdir=$am_cv_python_pyexecdir



  pkgpyexecdir=\${pyexecdir}/$PACKAGE



  fi
## end hacked version of checking python for >= 2.2

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
