# $Id: Makefile.in,v 1.10 2006/01/26 12:14:48 bcotton Exp $
#-----------------------------------------------------------------------
#;  Copyright (C) 2002-2006
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
#   Expects CC, CFLAGS, DOCS to be set by either environment variable 
#   or Make argument
#     CC is the compiler to use
#     CFLAGS is compile options
#     LIB defines libraries and search paths,
#
# Use Local if ARCH not externally defined
ARCH ?= Local
#------------------------------------------------------------------------

top_srcdir = ..

CC = /usr/bin/gcc
CPPFLAGS = 
CFLAGS = -g -O2 -fPIC
LDFLAGS =   -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib  -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib

CLIENT_CPPFLAGS =   -I/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/include  -I/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/include
CLIENT_CFLAGS =  -I/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/include
CLIENT_LDFLAGS =     -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib 
CLIENT_LIBS =  -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib  -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok  -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib  -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib -lcurl -L/usr/kerberos/lib -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib/lib -lxmlrpc_client -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok -Wl,--rpath -Wl,/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib      -lgsl -lgslcblas  -lz   \
         -lxmlrpc_client -lcurl -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib  -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok  -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib  -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib -lcurl -L/usr/kerberos/lib -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib/lib -lxmlrpc_client -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok -Wl,--rpath -Wl,/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib      

ALL_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/include  -I/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/include \
	  $(CLIENT_CPPFLAGS)  -DPACKAGE_NAME=\"Obit\" -DPACKAGE_TARNAME=\"obit\" -DPACKAGE_VERSION=\"1.0\" -DPACKAGE_STRING=\"Obit\ 1.0\" -DPACKAGE_BUGREPORT=\"bcotton@nrao.edu\" -DPACKAGE=\"obit\" -DVERSION=\"1.0\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIBCFITSIO=1 -DHAVE_CFITSIO=1 -DHAVE__TMP_DUMMY1_FFTW3_H=1 -DHAVE_FFTW3=1 -DHAVE__TMP_DUMMY3_FFTW_H=1 -DHAVE__TMP_DUMMY1_GSL_H=1 -DHAVE_GSL=1 -DHAVE__TMP_DUMMY1_XMLRPC_H=1 -DHAVE_XMLRPC=1 -DHAVE__TMP_DUMMY1_ZLIB_H=1 -DHAVE_ZLIB=1 -DHAVE_PLPLOT=1 -D_LARGEFILE_SOURCE=1 -DHAVE_FSEEKO=1 -D_FILE_OFFSET_BITS=64 
ALL_CFLAGS = $(CFLAGS) -I/export/data_2/fits/ObitTest/test/ObitInstall/other//include/glib-2.0 -I/export/data_2/fits/ObitTest/test/ObitInstall/other//lib/glib-2.0/include    -I/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/include  \
	 $(CLIENT_CFLAGS)
ALL_LDFLAGS = $(LDFLAGS)  -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib   \
	  -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib  $(CLIENT_LDFLAGS)

LIBS = ../lib/$(ARCH)/libObit.a -lcfitsio   -L/export/data_2/fits/ObitTest/test/ObitInstall/other//lib -lglib-2.0   \
	-lgsl -lgslcblas   $(CLIENT_LIBS)   -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib -L/usr/lib/gcc/i386-redhat-linux/3.4.6 -L/usr/lib/gcc/i386-redhat-linux/3.4.6/../../.. -lfrtbegin -lg2c -lm -lgcc_s

# get list of all c source files (*.c) files
AllC := $(wildcard *.c)
EXECU   := $(patsubst %.c,%, $(AllC))
OBJECTS := $(wildcard ../lib/LINUX/*.o)
TARGETS := testCleanVis testClean testFeather testUVSub

# get list of all c header files (*.h) files
HEADERS := $(wildcard ../include/*.h)

all: $(TARGETS)

# For specific executables
$(EXECU): % : %.c ../lib/$(ARCH)/libObit.a  
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $* $(LIBS)

# generic C compile/link
# $(TARGETS): $*.c $(HEADERS) ../lib/linux/libObit.a
# 	$(CC) $(CFLAGS)  %.c -o $< $(LIB)


clean:
	rm -f $(TARGETS)
	rm -f *.o

