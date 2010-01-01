# $Id: Makefile.in 15 2008-08-12 10:49:55Z bill.cotton $
#-----------------------------------------------------------------------
#;  Copyright (C) 2002-2008
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
#------------------------------------------------------------------------

top_srcdir = ..

# directory for architecture dependent executable files
BINDIR = $(top_srcdir)/bin/

CPPFLAGS = 
CFLAGS = -g -O2 -fPIC -Wall -fno-strict-aliasing -DFASTOBITMEM
CFLAGS = -g  -fPIC -Wall -fno-strict-aliasing -DFASTOBITMEM
LDFLAGS =   -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib
CC = /usr/bin/gcc

# xmlrpc client stuff
CLIENT_CPPFLAGS =   -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include
CLIENT_CFLAGS =  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include
CLIENT_LDFLAGS =     -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib 
CLIENT_LIBS =   -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib -lcurl -L/usr/kerberos/lib -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib/lib -lxmlrpc_client -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok -Wl,--rpath -Wl,/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib     -lz    -lxmlrpc_client -lcurl -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib -lcurl -L/usr/kerberos/lib -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib/lib -lxmlrpc_client -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok -Wl,--rpath -Wl,/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib     
SERVER_CPPFLAGS =   -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include
SERVER_CFLAGS = 
SERVER_LDFLAGS = $(LDFLAGS)   
SERVER_LIBS =   -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -lxmlrpc_server_abyss -lxmlrpc_server -lxmlrpc_abyss  -lpthread -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok   

ALL_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include \
	 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/include -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include   -I/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/include -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include     -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include   -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include/plplot \
	$(CLIENT_CPPFLAGS) $(SERVER_CPPFLAGS) -DPACKAGE_NAME=\"ObitSD\" -DPACKAGE_TARNAME=\"obitsd\" -DPACKAGE_VERSION=\"1.3\" -DPACKAGE_STRING=\"ObitSD\ 1.3\" -DPACKAGE_BUGREPORT=\"bcotton@nrao.edu\" -DPACKAGE=\"obitsd\" -DVERSION=\"1.3\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIBCFITSIO=1 -DHAVE_CFITSIO=1 -DHAVE__TMP_DUMMY1_FFTW3_H=1 -DHAVE_FFTW3=1 -DHAVE__TMP_DUMMY3_FFTW_H=1 -DHAVE__TMP_DUMMY1_GSL_H=1 -DHAVE_GSL=1 -DHAVE__TMP_DUMMY1_XMLRPC_H=1 -DHAVE_XMLRPC=1 -DHAVE__TMP_DUMMY1_ZLIB_H=1 -DHAVE_ZLIB=1 -DHAVE_PLPLOT=1 -D_LARGEFILE_SOURCE=1 -DHAVE_FSEEKO=1 -D_FILE_OFFSET_BITS=64 -DOBIT_THREADS_ENABLED=1 -DHAVE_LIBOBIT=1 -DHAVE_OBIT=1 
ALL_CFLAGS = $(CFLAGS) -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include     $(CLIENT_CFLAGS) $(SERVER_CFLAGS) \
	  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include/plplot  
ALL_LDFLAGS = $(LDFLAGS)  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/lib/ -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lglib-2.0   -lpthread -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lgthread-2.0 -lglib-2.0    -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  \
	  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib   $(CLIENT_LDFLAGS) $(SERVER_LDFLAGS)

LIBS = ../lib/libObitSD.a -lObit  -lcfitsio    -lfftw3f \
	-L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lglib-2.0   -lgsl -lgslcblas  $(CLIENT_LIBS) $(SERVER_LIBS) -lplplot -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib /export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/libplplot.so /usr/lib/libfreetype.so -lz /export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/libcsirocsa.so -ldl -lm -Wl,--rpath -Wl,/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -Wl,--rpath -Wl,/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib  \
	  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib -L/usr/lib/gcc/i386-redhat-linux/3.4.6 -L/usr/lib/gcc/i386-redhat-linux/3.4.6/../../.. -lfrtbegin -lg2c -lm -lgcc_s -lpthread -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lgthread-2.0 -lglib-2.0  

# get list of all c source files (*.c) files
AllC    := $(wildcard *.c)
EXECU   := $(patsubst %.c,%, $(AllC))
TARGETS := $(addprefix $(BINDIR),$(EXECU))

# get list of all c header files (*.h) files
HEADERS := $(wildcard ../include/*.h)

all: $(TARGETS)

# generic C compile/link
$(TARGETS): $(BINDIR)% : %.c ../lib/libObitSD.a  
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $* $(LIBS)
	mv $* $@

# For specific executables
$(EXECU): % : %.c ../lib/libObitSD.a  
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $* $(LIBS)
	mv $* $(BINDIR)

# special case, junk - ignore
# junk:

clean:
	rm -f $(TARGETS)
	rm -f *.o

