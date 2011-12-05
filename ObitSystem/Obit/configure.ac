# $Id$
#-----------------------------------------------------------------------
#;  Copyright (C) 2007-2009
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

AC_INIT(Obit, 1.0, bcotton@nrao.edu)
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_SRCDIR(src/Obit.c)

# Checks for programs.
CFLAGS="-g -O3"
AC_PROG_CC
AC_PROG_RANLIB
#AM_PATH_PYTHON(2.5)
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
# Check for SSE
AC_CHECK_SSE

# Checks for libraries.
AM_PATH_GLIB_2_0
AC_CHECK_LIB(m, pow strtol)
AC_PATH_CFITSIO
AC_PATH_FFTW3
AC_PATH_FFTW
AC_PATH_GSL
AC_PATH_XMLRPC
AC_PATH_ZLIB
AC_PATH_PLPLOT
AM_PATH_PGPLOT
AC_PATH_PYTHON2_5
AC_PATH_WVR

# Checks for library functions
AC_FUNC_FSEEKO

# Checks for system services.
AC_SYS_LARGEFILE

if test $ac_have_cfitsio = no; then
  CFITSIO_CPPFLAGS="-I../dummy_cfitsio"
  CFITSIO_LIBS="../dummy_cfitsio/libcfitsio.a"
fi
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

# GThread stuff - horrible hack for pkg-config bug
GTHREAD_LIBS=""
TGTHREAD_LIBS=`pkg-config --libs gthread-2.0 | sed -e 's/-pthread/-lpthread/g'`
AC_CHECK_LIB([gthread-2.0], g_thread_init, [ac_gthread=yes], [ac_gthread=no
             AC_MSG_WARN([cannot find GTHREAD library])],[$TGTHREAD_LIBS])
if test $ac_gthread = yes; then
	AC_DEFINE(OBIT_THREADS_ENABLED, 1, [Define to 1 if you have gthreads.])
	GTHREAD_LIBS=$TGTHREAD_LIBS
	AC_SUBST(GTHREAD_LIBS)
fi

AC_CONFIG_FILES([dummy_cfitsio/Makefile])
AC_CONFIG_FILES([dummy_xmlrpc/Makefile])
AC_CONFIG_FILES([Makefile python/Makefile src/Makefile tasks/Makefile])
AC_OUTPUT
