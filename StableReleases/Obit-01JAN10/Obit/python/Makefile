# $Id: Makefile.in 15 2008-08-12 10:49:55Z bill.cotton $
#-----------------------------------------------------------------------
#;  Copyright (C) 2004-2008
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
#    Makefile to build library routines
#    Obit software
#
#   Expects CC, CFLAGS, to be set by either environment variable 
#   or Make argument
#     CC is the compiler to use
#     CFLAGS is compile options
#     INCDIR directories with include files
#     LIB defines libraries and search paths,
#
#------------------------------------------------------------------------

top_srcdir = ..

CC = /usr/bin/gcc
#CPPFLAGS = 
CPPFLAGS = 
CFLAGS = -g -O3 -fPIC -msse
LDFLAGS =   -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib

ALL_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include \
	 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include/ -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include/python2.5   -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include/ -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include/python2.5 \
	 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include/plplot  -DPACKAGE_NAME=\"Obit\" -DPACKAGE_TARNAME=\"obit\" -DPACKAGE_VERSION=\"1.0\" -DPACKAGE_STRING=\"Obit\ 1.0\" -DPACKAGE_BUGREPORT=\"bcotton@nrao.edu\" -DPACKAGE=\"obit\" -DVERSION=\"1.0\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_SSE=1 -DHAVE_LIBCFITSIO=1 -DHAVE_CFITSIO=1 -DHAVE__TMP_DUMMY1_FFTW3_H=1 -DHAVE_FFTW3=1 -DHAVE__TMP_DUMMY3_FFTW_H=1 -DHAVE__TMP_DUMMY1_GSL_H=1 -DHAVE_GSL=1 -DHAVE__TMP_DUMMY1_XMLRPC_H=1 -DHAVE_XMLRPC=1 -DHAVE__TMP_DUMMY1_ZLIB_H=1 -DHAVE_ZLIB=1 -DHAVE_PLPLOT=1 -DHAVE_PYTHON=1 -D_LARGEFILE_SOURCE=1 -DHAVE_FSEEKO=1 -D_FILE_OFFSET_BITS=64 -DOBIT_THREADS_ENABLED=1 
ALL_CFLAGS = $(CFLAGS) -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include     -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include/plplot \
	 
ALL_LDFLAGS = $(LDFLAGS)  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib    -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib   \
	 -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -L/export/data_1/users/bcotton/Software.dir/SVN/ObitInstall/other/lib/python2.5/site-packages 

LIBS = ../lib/libObit.a -lcfitsio    -lfftw3f -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lglib-2.0   \
	-lgsl -lgslcblas  -lplplot -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib /export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/libplplot.so /usr/lib/libfreetype.so -lz /export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/libcsirocsa.so -ldl -lm -Wl,--rpath -Wl,/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -Wl,--rpath -Wl,/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib    -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib -L/usr/lib/gcc/i386-redhat-linux/3.4.6 -L/usr/lib/gcc/i386-redhat-linux/3.4.6/../../.. -lfrtbegin -lg2c -lm -lgcc_s -lpthread -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lgthread-2.0 -lglib-2.0  

CLIENT_CPPFLAGS = $(ALL_CPPFLAGS)   -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include 
CLIENT_CFLAGS = $(ALL_CFLAGS)  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include
CLIENT_LDFLAGS = $(ALL_LDFLAGS)     -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib 
CLIENT_LIBS = $(LIBS) -lz   -lxmlrpc_client -lcurl -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib -lcurl -L/usr/kerberos/lib -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib/lib -lxmlrpc_client -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok -Wl,--rpath -Wl,/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib     

SERVER_CPPFLAGS =   -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include
SERVER_CFLAGS = 
SERVER_LDFLAGS = $(LDFLAGS)   
SERVER_LIBS =   -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -lxmlrpc_server_abyss -lxmlrpc_server -lxmlrpc_abyss  -lpthread -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok   

SWIG = /export/users/bcotton/bin/swig

# Libraries in case they've changed
MYLIBS := $(wildcard ../lib/lib*.a)

# Do everything in one big module
TARGETS := Obit.so 

all: $(TARGETS)

# Build shared library for python interface
$(TARGETS): %.so:setupdata.py $(MYLIBS)
	python makesetup.py
	python setup.py build install --install-lib=.

# Build python/Obit interface
interface: Obit_wrap.c
	echo "rebuild Obit/python interface"

# Use Swig to build python/c interface if necessary
Obit_wrap.c: *.inc 
	echo "%module Obit" > Obit.i
	cat ObitTypeMaps.swig >> Obit.i
	cat *.inc >> Obit.i
	$(SWIG) -python  $(INCDIR) Obit.i

# Write compile/link info to setupdata.py
setupdata.py: Makefile
	echo "CFLAGS='$(ALL_CFLAGS)'">setupdata.py
	echo "CPPFLAGS='$(CLIENT_CPPFLAGS) $(SERVER_CPPFLAGS)'">>setupdata.py
	echo "LDFLAGS='$(CLIENT_LDFLAGS) $(SERVER_LDFLAGS)'">>setupdata.py
	echo "LIBS='$(CLIENT_LIBS) $(SERVER_LIBS)'">>setupdata.py

clean:
	rm -f Obit.i *.o *.so *.pyc	
	rm -rf build
