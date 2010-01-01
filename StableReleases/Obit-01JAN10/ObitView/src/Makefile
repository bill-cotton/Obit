# $Id: Makefile.in 2 2008-06-10 15:32:27Z bill.cotton $
#-----------------------------------------------------------------------
#;  Copyright (C) 2005
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
#    Makefile to compile .c files and leave object files locally
#
#------------------------------------------------------------------------

top_srcdir = ..
LIBDIR = $(top_srcdir)/lib/
CC = /usr/bin/gcc
CPPFLAGS = 
CFLAGS = -g -O2 -fPIC -Wall -DFASTOBITMEM -fno-strict-aliasing
#CFLAGS = -g -fPIC -Wall -DFASTOBITMEM -fno-strict-aliasing
DEFS = -DPACKAGE_NAME=\"ObitView\" -DPACKAGE_TARNAME=\"obitview\" -DPACKAGE_VERSION=\"1.0\" -DPACKAGE_STRING=\"ObitView\ 1.0\" -DPACKAGE_BUGREPORT=\"bcotton@nrao.edu\" -DPACKAGE=\"obitview\" -DVERSION=\"1.0\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_MOTIF=1 -DHAVE_PLPLOT=1 -DHAVE_LIBCFITSIO=1 -DHAVE_CFITSIO=1 -DHAVE__TMP_DUMMY1_FFTW3_H=1 -DHAVE_FFTW3=1 -DHAVE__TMP_DUMMY3_FFTW_H=1 -DHAVE__TMP_DUMMY1_GSL_H=1 -DHAVE_GSL=1 -DHAVE__TMP_DUMMY1_XMLRPC_H=1 -DHAVE_XMLRPC=1 -DHAVE__TMP_DUMMY1_ZLIB_H=1 -DHAVE_ZLIB=1 -D_FILE_OFFSET_BITS=64 -DOBIT_THREADS_ENABLED=1 -DHAVE_LIBOBIT=1 -DHAVE_OBIT=1 

ALL_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include  \
          -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/include -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include   -I/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/include -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include     -I/usr/X11R6/include -DPACKAGE_NAME=\"ObitView\" -DPACKAGE_TARNAME=\"obitview\" -DPACKAGE_VERSION=\"1.0\" -DPACKAGE_STRING=\"ObitView\ 1.0\" -DPACKAGE_BUGREPORT=\"bcotton@nrao.edu\" -DPACKAGE=\"obitview\" -DVERSION=\"1.0\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_MOTIF=1 -DHAVE_PLPLOT=1 -DHAVE_LIBCFITSIO=1 -DHAVE_CFITSIO=1 -DHAVE__TMP_DUMMY1_FFTW3_H=1 -DHAVE_FFTW3=1 -DHAVE__TMP_DUMMY3_FFTW_H=1 -DHAVE__TMP_DUMMY1_GSL_H=1 -DHAVE_GSL=1 -DHAVE__TMP_DUMMY1_XMLRPC_H=1 -DHAVE_XMLRPC=1 -DHAVE__TMP_DUMMY1_ZLIB_H=1 -DHAVE_ZLIB=1 -D_FILE_OFFSET_BITS=64 -DOBIT_THREADS_ENABLED=1 -DHAVE_LIBOBIT=1 -DHAVE_OBIT=1 
ALL_CFLAGS = $(CFLAGS) -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include     -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include

CC = /usr/bin/gcc

SERVER_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include \
           -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/include -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include   -I/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/include -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include   \
          -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include -DPACKAGE_NAME=\"ObitView\" -DPACKAGE_TARNAME=\"obitview\" -DPACKAGE_VERSION=\"1.0\" -DPACKAGE_STRING=\"ObitView\ 1.0\" -DPACKAGE_BUGREPORT=\"bcotton@nrao.edu\" -DPACKAGE=\"obitview\" -DVERSION=\"1.0\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_MOTIF=1 -DHAVE_PLPLOT=1 -DHAVE_LIBCFITSIO=1 -DHAVE_CFITSIO=1 -DHAVE__TMP_DUMMY1_FFTW3_H=1 -DHAVE_FFTW3=1 -DHAVE__TMP_DUMMY3_FFTW_H=1 -DHAVE__TMP_DUMMY1_GSL_H=1 -DHAVE_GSL=1 -DHAVE__TMP_DUMMY1_XMLRPC_H=1 -DHAVE_XMLRPC=1 -DHAVE__TMP_DUMMY1_ZLIB_H=1 -DHAVE_ZLIB=1 -D_FILE_OFFSET_BITS=64 -DOBIT_THREADS_ENABLED=1 -DHAVE_LIBOBIT=1 -DHAVE_OBIT=1 
SERVER_CFLAGS = $(CFLAGS)  -I/usr/X11R6/include -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include     -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include 

CLIENT_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include \
           -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/include -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include   -I/export/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/include -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include   \
          -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include -DPACKAGE_NAME=\"ObitView\" -DPACKAGE_TARNAME=\"obitview\" -DPACKAGE_VERSION=\"1.0\" -DPACKAGE_STRING=\"ObitView\ 1.0\" -DPACKAGE_BUGREPORT=\"bcotton@nrao.edu\" -DPACKAGE=\"obitview\" -DVERSION=\"1.0\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_MOTIF=1 -DHAVE_PLPLOT=1 -DHAVE_LIBCFITSIO=1 -DHAVE_CFITSIO=1 -DHAVE__TMP_DUMMY1_FFTW3_H=1 -DHAVE_FFTW3=1 -DHAVE__TMP_DUMMY3_FFTW_H=1 -DHAVE__TMP_DUMMY1_GSL_H=1 -DHAVE_GSL=1 -DHAVE__TMP_DUMMY1_XMLRPC_H=1 -DHAVE_XMLRPC=1 -DHAVE__TMP_DUMMY1_ZLIB_H=1 -DHAVE_ZLIB=1 -D_FILE_OFFSET_BITS=64 -DOBIT_THREADS_ENABLED=1 -DHAVE_LIBOBIT=1 -DHAVE_OBIT=1 
CLIENT_CFLAGS = $(CFLAGS)  -I/usr/X11R6/include -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include     -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include 


# teach how to compile
.c.o:	
	$(CC) -c $(ALL_CFLAGS) $(ALL_CPPFLAGS) $<
	mv $@.o $(LIBDIR)

# get list of all c source files (*.c) files
AllC := $(wildcard *.c)
OBJECTS := $(patsubst %.c,%.o, $(AllC))

# get list of all c header files (*.h) files
HEADERS := $(wildcard ../include/*.h)

CTARGETS := $(addprefix $(LIBDIR),$(OBJECTS))


all:  $(CTARGETS) 

# generic C compile
$(CTARGETS): $(LIBDIR)%.o: %.c $(wildcard ../include/%*.h)
	$(CC) -c $(ALL_CPPFLAGS) $(ALL_CFLAGS) $*.c
	mv ./$*.o $@

# RPC server interface
../lib/XMLRPCserver.$(LIBDIR)%.o: %.c $(wildcard ../include/%*.h)
	$(CC) -c $(SERVER_CPPFLAGS) $(SERVER_CFLAGS) XMLRPCserver.c
	mv XMLRPCserver.o $(LIBDIR)

# need server interface
../lib/requestbox.o: $(LIBDIR)%.o: %.c $(wildcard ../include/%*.h)
	$(CC) -c $(SERVER_CPPFLAGS) $(SERVER_CFLAGS) requestbox.c
	mv requestbox.o $(LIBDIR)

clean:
	rm -f $(CTARGETS)
	rm -f *.o
