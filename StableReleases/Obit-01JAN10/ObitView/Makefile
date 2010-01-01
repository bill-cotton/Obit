# $Id: Makefile.in 15 2008-08-12 10:49:55Z bill.cotton $
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
#
#    ObitView: Image viewer for Obit
#
#----------------------------------------------------
CPPFLAGS = 
CC = /usr/bin/gcc

# DEBUG is the debugger/optimizer flag, 
#-------------------------------------
# O is some general level of optimization
DEBUG = -O2
OSFLAGS = 
#DEBUG = -g   # debug

# RANLIB = library "randomizer", ranlib if needed, else blank
#-----------------------------------------------------------
RANLIB = ranlib

# name of archive file
ARCHIVE = obitviewBack1.0
DISTRIB = obitview1.0
# Name of this directory
DIRN = ObitView

#------------------------------------------------------------------------
TARGETS =  srcupdate libupdate ObitView  ObitMess

top_srcdir = .

CPPFLAGS = 
CFLAGS = -g -O2 -fPIC  -Wall -DFASTOBITMEM -fno-strict-aliasing
LDFLAGS =   -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib
OBIT = /export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit
CC = /usr/bin/gcc

SERVER_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include \
           -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/include -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include   -I/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/include -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include   \
          -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include -DPACKAGE_NAME=\"ObitView\" -DPACKAGE_TARNAME=\"obitview\" -DPACKAGE_VERSION=\"1.0\" -DPACKAGE_STRING=\"ObitView\ 1.0\" -DPACKAGE_BUGREPORT=\"bcotton@nrao.edu\" -DPACKAGE=\"obitview\" -DVERSION=\"1.0\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_MOTIF=1 -DHAVE_PLPLOT=1 -DHAVE_LIBCFITSIO=1 -DHAVE_CFITSIO=1 -DHAVE__TMP_DUMMY1_FFTW3_H=1 -DHAVE_FFTW3=1 -DHAVE__TMP_DUMMY3_FFTW_H=1 -DHAVE__TMP_DUMMY1_GSL_H=1 -DHAVE_GSL=1 -DHAVE__TMP_DUMMY1_XMLRPC_H=1 -DHAVE_XMLRPC=1 -DHAVE__TMP_DUMMY1_ZLIB_H=1 -DHAVE_ZLIB=1 -D_FILE_OFFSET_BITS=64 -DOBIT_THREADS_ENABLED=1 -DHAVE_LIBOBIT=1 -DHAVE_OBIT=1 
SERVER_CFLAGS = $(CFLAGS)  -I/usr/X11R6/include -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include     -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include
SERVER_LDFLAGS = $(LDFLAGS)  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/lib/ -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lglib-2.0   -lpthread -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lgthread-2.0 -lglib-2.0   -L/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/lib/ -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lglib-2.0   -lpthread -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lgthread-2.0 -lglib-2.0    -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib   -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib \
          -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib    -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib   -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib 

SERVER_LIBS = lib/libObitView.a -lXm  -lXt -lX11  -L/usr/X11R6/lib -lObit  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lglib-2.0   \
        -lcfitsio    -lfftw3f  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib -lcurl -L/usr/kerberos/lib -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib/lib -lxmlrpc_client -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok -Wl,--rpath -Wl,/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib      -lSM -lICE \
        -lgsl -lgslcblas  -lz    -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -lxmlrpc_server_abyss -lxmlrpc_server -lxmlrpc_abyss  -lpthread -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok   -lpthread -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lgthread-2.0 -lglib-2.0  

CLIENT_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include \
           -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/include -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include   -I/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/include -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include   \
          -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include -DPACKAGE_NAME=\"ObitView\" -DPACKAGE_TARNAME=\"obitview\" -DPACKAGE_VERSION=\"1.0\" -DPACKAGE_STRING=\"ObitView\ 1.0\" -DPACKAGE_BUGREPORT=\"bcotton@nrao.edu\" -DPACKAGE=\"obitview\" -DVERSION=\"1.0\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_MOTIF=1 -DHAVE_PLPLOT=1 -DHAVE_LIBCFITSIO=1 -DHAVE_CFITSIO=1 -DHAVE__TMP_DUMMY1_FFTW3_H=1 -DHAVE_FFTW3=1 -DHAVE__TMP_DUMMY3_FFTW_H=1 -DHAVE__TMP_DUMMY1_GSL_H=1 -DHAVE_GSL=1 -DHAVE__TMP_DUMMY1_XMLRPC_H=1 -DHAVE_XMLRPC=1 -DHAVE__TMP_DUMMY1_ZLIB_H=1 -DHAVE_ZLIB=1 -D_FILE_OFFSET_BITS=64 -DOBIT_THREADS_ENABLED=1 -DHAVE_LIBOBIT=1 -DHAVE_OBIT=1 
CLIENT_CFLAGS = $(CFLAGS)  -I/usr/X11R6/include -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include     -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include
CLIENT_LDFLAGS = $(LDFLAGS)  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/lib/ -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lglib-2.0   -lpthread -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lgthread-2.0 -lglib-2.0   -L/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/lib/ -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lglib-2.0   -lpthread -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lgthread-2.0 -lglib-2.0    -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib   -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib \
          -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib    -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib   -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib 

CLIENT_LIBS = lib/libObitView.a -lXm  -lXt -lX11  -L/usr/X11R6/lib -lObit   -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lglib-2.0   \
        -lcfitsio    -lfftw3f  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib -lcurl -L/usr/kerberos/lib -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib/lib -lxmlrpc_client -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok -Wl,--rpath -Wl,/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib      -lSM -lICE \
        -lgsl -lgslcblas  -lz   \
         -lxmlrpc_client -lcurl -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib -lcurl -L/usr/kerberos/lib -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib/lib -lxmlrpc_client -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok -Wl,--rpath -Wl,/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib     -lpthread -lpthread -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib -lgthread-2.0 -lglib-2.0  


all:  $(TARGETS)

# update source/object directory
srcupdate: 
	cd src; $(MAKE)

# update library directory
libupdate: 
	cd lib; $(MAKE) RANLIB="$(RANLIB)"

# Link Obitview
ObitView: src/*.c  srcupdate ObitView.c 
	$(CC) ObitView.c -o ObitView  $(SERVER_CFLAGS) $(SERVER_CPPFLAGS) \
	$(SERVER_LDFLAGS) $(SERVER_LIBS) $(CLIENT_LIBS) 

# Link ObitMess
ObitMess: src/*.c  srcupdate ObitMess.c 
	$(CC) ObitMess.c -o ObitMess  $(SERVER_CFLAGS) $(SERVER_CPPFLAGS) \
	$(SERVER_LDFLAGS) $(SERVER_LIBS) $(CLIENT_LIBS) 

# Link clientTest
clientTest:   clientTest.c 
	$(CC) clientTest.c -o clientTest  -Iinclude $(CLIENT_CFLAGS) \
	$(CLIENT_CPPFLAGS) $(CLIENT_LDFLAGS) \
	$(CLIENT_LIBS)

# Link clientFCopy
clientFCopy:   clientFCopy.c 
	$(CC) clientFCopy.c -o clientFCopy  -Iinclude $(CLIENT_CFLAGS) \
	$(CLIENT_CPPFLAGS) $(CLIENT_LDFLAGS) \
	$(CLIENT_LIBS)

# test run ObitView
testObitView: 
	ObitView aaaSomeFile.fits

# test run ObitMess
testObitMess: 
	ObitMess &
	python testObitMess.py

# Copy to where it should go
install: /export/users/bcotton/Software.dir/SVN/ObitInstall/bin
	/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/ObitView/install-sh -s ObitView /export/users/bcotton/Software.dir/SVN/ObitInstall/bin/
	/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/ObitView/install-sh -s ObitMess /export/users/bcotton/Software.dir/SVN/ObitInstall/bin/

/export/users/bcotton/Software.dir/SVN/ObitInstall/bin:
	mkdir /export/users/bcotton/Software.dir/SVN/ObitInstall/bin

# make gzipped tar backup of directory
backup:	
	cd ..;rm -f $(ARCHIVE).tgz
	cd ..;gtar czvf $(ARCHIVE).tgz \
		$(DIRN)/README $(DIRN)/LICENSE $(DIRN)/Makefile*  \
		$(DIRN)/*.c  $(DIRN)/ObitView.hlp $(DIRN)/someFile.fits* \
		$(DIRN)/config*  $(DIRN)/aclocal.m4  \
		$(DIRN)/install-sh  $(DIRN)/missing $(DIRN)/mkinstalldirs  \
		$(DIRN)/changes  $(DIRN)/Notes.text  $(DIRN)/notes.text  \
		$(DIRN)/src/*.c  $(DIRN)/src/Makefile.in   $(DIRN)/src/Makefile \
		$(DIRN)/include/*.h   \
		$(DIRN)/m4/*.m4  \
		$(DIRN)/CVS/* $(DIRN)/m4/CVS/* $(DIRN)/src/CVS/* \
		$(DIRN)/include/CVS/*

# make gzipped tar distribution
distrib:	
	cd ..;rm -f $(DISTRIB).tgz
	cd ..;gtar czvf $(DISTRIB).tgz \
		$(DIRN)/README $(DIRN)/LICENSE $(DIRN)/Makefile.in  \
		$(DIRN)/ObitView.c  $(DIRN)/ObitView.hlp  $(DIRN)/someFile.fits* \
		$(DIRN)/config*  $(DIRN)/aclocal.m4  \
		$(DIRN)/install-sh  $(DIRN)/missing $(DIRN)/mkinstalldirs  \
		$(DIRN)/changes  $(DIRN)/Notes.text  $(DIRN)/notes.text  \
		$(DIRN)/src/*.c  $(DIRN)/src/Makefile.in \
		$(DIRN)/m4/*.m4  \
		$(DIRN)/include/*.h 
clean:
	cd src;      $(MAKE) clean
	rm -f ObitView ObitMess clientTest
	rm -f *.o
