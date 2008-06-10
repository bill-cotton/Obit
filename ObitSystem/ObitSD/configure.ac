# $Id$
#-----------------------------------------------------------------------
#;  Copyright (C) 2007,2008
#;  Associated Universities, Inc. Washington DC, USA.
#;
#;  This program is free software; you can redistribute it and/or
#;  modify it under the terms of the GNU General Public License as
#;  published by the Free Software Foundation; either version 2 of
#;  the License, or (at your option) any later version.
#;
#;  This program is distributed in the hope that it will be useful,
#;  but WITHOUT ANY WARRANTY; without even the implied warranty of
#;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#;  GNU General Public License for more details.
#;
#
#;  You should have received a copy of the GNU General Public
#;  License along with this program; if not, write to the Free
#;  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
#;  MA 02139, USA.
#;
#;  Correspondence concerning this software should be addressed as follows:
#;         Internet email: bcotton@nrao.edu
#;         Postal address: W. D. Cotton
#;                         National Radio Astronomy Observatory
#;                         520 Edgemont Road
#;                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

AC_INIT(ObitSD, 1.3, bcotton@nrao.edu)
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR(src/ObitOTF.c)

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB
AM_PATH_PYTHON(2.2)
AC_PATH_PROG(SWIG, swig)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T

# Compiler options
# hack for -fPIC
AC_CHECK_FPIC

# Checks for libraries.
AM_PATH_GLIB_2_0
AC_CHECK_LIB(m, pow strtol)
AC_PATH_CFITSIO
AC_PATH_FFTW3
AC_PATH_FFTW
AC_PATH_OBIT
AC_PATH_GSL
AC_PATH_XMLRPC
AC_PATH_ZLIB
AC_PATH_PLPLOT
AM_PATH_PGPLOT

# Python stuff - standard python found?
cat <<_ACEOF >conftest.py
import distutils.sysconfig
print distutils.sysconfig.get_python_inc()
_ACEOF
ac_python_inc=`$PYTHON conftest.py`
PYTHON_CPPFLAGS="-I$ac_python_inc"
AC_SUBST(PYTHON_CPPFLAGS)

# Python libs
cat <<_ACEOF >conftest.py
import distutils.sysconfig
print distutils.sysconfig.get_python_lib()
_ACEOF
ac_python_lib=`$PYTHON conftest.py`
PYTHON_LDFLAGS="-I$ac_python_inc"
AC_SUBST(PYTHON_LDFLAGS)

# Checks for library functions
AC_FUNC_FSEEKO

# Checks for system services.
AC_SYS_LARGEFILE

# Use dummy version of xmlrpc if not available
if test $ac_have_xmlrpc = no; then
  XMLRPC_CPPFLAGS="-I../dummy_xmlrpc"
  XMLRPC_CLIENT_CPPFLAGS="-I../dummy_xmlrpc"
  XMLRPC_SERVER_CPPFLAGS="-I../dummy_xmlrpc"
  XMLRPC_LDFLAGS="-L../dummy_xmlrpc"
  XMLRPC_CLIENT_LDFLAGS="-L../dummy_xmlrpc"
  XMLRPC_SERVER_LDFLAGS="-L../dummy_xmlrpc"
  XMLRPC_LIBS="../dummy_xmlrpc/libxmlrpc.a"
  XMLRPC_CLIENT_LIBS="../dummy_xmlrpc/libxmlrpc.a"
  XMLRPC_SERVER_LIBS="../dummy_xmlrpc/libxmlrpc.a"
  WWWLIB_WL_RPATH=" "
  WWW_CPPFLAGS=" "
  WWW_LIBS=" "
fi

AC_CONFIG_FILES([Makefile python/Makefile src/Makefile test/Makefile tasks/Makefile])
AC_OUTPUT
