# $Id: Makefile.in 2 2008-06-10 15:32:27Z bill.cotton $
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

CC = gcc
CPPFLAGS = 
CFLAGS = -g -O3 -fPIC -msse
LDFLAGS =  -L/export/users/bcotton/opt/plplot/lib -L/export/users/bcotton/opt/plplot/lib

CLIENT_CPPFLAGS =  -I/export/users/bcotton/opt/xmlrpc/include 
CLIENT_CFLAGS = 
CLIENT_LDFLAGS = -L/export/users/bcotton/opt/xmlrpc/lib   -L/export/users/bcotton/opt/xmlrpc/lib  
CLIENT_LIBS =  -L/export/users/bcotton/opt/xmlrpc/lib  -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok  -L/export/users/bcotton/opt/xmlrpc/lib  -L/usr/lib -lwwwxml -lxmltok -lxmlparse -lwwwzip -lwwwinit -lwwwapp -lwwwhtml -lwwwtelnet -lwwwnews -lwwwhttp -lwwwmime -lwwwgopher -lwwwftp -lwwwfile -lwwwdir -lwwwcache -lwwwstream -lwwwmux -lwwwtrans -lwwwcore -lwwwutils -lmd5 -ldl -lz -L/usr/lib -lcurl -L/usr/kerberos/lib -lssl -lcrypto -lgssapi_krb5 -lkrb5 -lcom_err -lk5crypto -lresolv -ldl -lz -lgssapi_krb5 -lkrb5 -lk5crypto -lcom_err -lresolv -L/usr/kerberos/lib -lidn -lssl -lcrypto -lssl -lcrypto -lgssapi_krb5 -lkrb5 -lcom_err -lk5crypto -lresolv -ldl -lz -lz -lxmlrpc_client -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok  -Wl,--rpath -Wl,/usr/lib      -lgsl -lgslcblas  -lz   \
         -lxmlrpc_client -lcurl -L/export/users/bcotton/opt/xmlrpc/lib  -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok  -L/export/users/bcotton/opt/xmlrpc/lib  -L/usr/lib -lwwwxml -lxmltok -lxmlparse -lwwwzip -lwwwinit -lwwwapp -lwwwhtml -lwwwtelnet -lwwwnews -lwwwhttp -lwwwmime -lwwwgopher -lwwwftp -lwwwfile -lwwwdir -lwwwcache -lwwwstream -lwwwmux -lwwwtrans -lwwwcore -lwwwutils -lmd5 -ldl -lz -L/usr/lib -lcurl -L/usr/kerberos/lib -lssl -lcrypto -lgssapi_krb5 -lkrb5 -lcom_err -lk5crypto -lresolv -ldl -lz -lgssapi_krb5 -lkrb5 -lk5crypto -lcom_err -lresolv -L/usr/kerberos/lib -lidn -lssl -lcrypto -lssl -lcrypto -lgssapi_krb5 -lkrb5 -lcom_err -lk5crypto -lresolv -ldl -lz -lz -lxmlrpc_client -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok  -Wl,--rpath -Wl,/usr/lib      

ALL_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/include -I/usr/include/cfitsio \
	-I/export/users/bcotton/opt/fftw/include  $(CLIENT_CPPFLAGS)  -DPACKAGE_NAME=\"Obit\" -DPACKAGE_TARNAME=\"obit\" -DPACKAGE_VERSION=\"1.0\" -DPACKAGE_STRING=\"Obit\ 1.0\" -DPACKAGE_BUGREPORT=\"bcotton@nrao.edu\" -DPACKAGE=\"obit\" -DVERSION=\"1.0\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_SSE=1 -DHAVE__TMP_DUMMY2_FITSIO_H=1 -DHAVE_LIBCFITSIO=1 -DHAVE_CFITSIO=1 -DHAVE__TMP_DUMMY1_FFTW3_H=1 -DHAVE__TMP_DUMMY3_FFTW_H=1 -DHAVE_FFTW=1 -DHAVE__TMP_DUMMY1_GSL_H=1 -DHAVE_GSL=1 -DHAVE__TMP_DUMMY3_XMLRPC_H=1 -DHAVE_XMLRPC=1 -DHAVE__TMP_DUMMY1_ZLIB_H=1 -DHAVE_ZLIB=1 -DHAVE_PLPLOT=1 -DHAVE_PYTHON=1 -D_LARGEFILE_SOURCE=1 -DHAVE_FSEEKO=1 -D_FILE_OFFSET_BITS=64 -DOBIT_THREADS_ENABLED=1 
ALL_CFLAGS = $(CFLAGS) -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include     \
	 $(CLIENT_CFLAGS)
ALL_LDFLAGS = $(LDFLAGS)  -L/export/users/bcotton/opt/fftw/lib  \
	   $(CLIENT_LDFLAGS)

LIBS = ../lib/$(ARCH)/libObit.a -lcfitsio   -lfftw -lrfftw -lglib-2.0   \
	-lgsl -lgslcblas   $(CLIENT_LIBS)   -L/export/users/bcotton/opt/plplot/lib -L/usr/lib/gcc/i386-redhat-linux/3.4.6 -L/usr/lib/gcc/i386-redhat-linux/3.4.6/../../.. -lfrtbegin -lg2c -lm -lgcc_s

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

