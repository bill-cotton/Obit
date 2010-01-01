# $Id: Makefile.in 2 2008-06-10 15:32:27Z bill.cotton $
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
# library directory for object files
LIBDIR = $(top_srcdir)/lib/

CC = /usr/bin/gcc
CPPFLAGS =  
# Turn off memory usage checking
CFLAGS = -Wall -fno-strict-aliasing -g -O3 -fPIC -msse -DFASTOBITMEM 
CFLAGS = -Wall -fno-strict-aliasing -g -fPIC -msse -DFASTOBITMEM  # debug
LDFLAGS =   -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib  -L/export/users/bcotton/Software.dir/SVN/ObitInstall/other/lib
LIBS = 
DEFS = -DPACKAGE_NAME=\"Obit\" -DPACKAGE_TARNAME=\"obit\" -DPACKAGE_VERSION=\"1.0\" -DPACKAGE_STRING=\"Obit\ 1.0\" -DPACKAGE_BUGREPORT=\"bcotton@nrao.edu\" -DPACKAGE=\"obit\" -DVERSION=\"1.0\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_SSE=1 -DHAVE_LIBCFITSIO=1 -DHAVE_CFITSIO=1 -DHAVE__TMP_DUMMY1_FFTW3_H=1 -DHAVE_FFTW3=1 -DHAVE__TMP_DUMMY3_FFTW_H=1 -DHAVE__TMP_DUMMY1_GSL_H=1 -DHAVE_GSL=1 -DHAVE__TMP_DUMMY1_XMLRPC_H=1 -DHAVE_XMLRPC=1 -DHAVE__TMP_DUMMY1_ZLIB_H=1 -DHAVE_ZLIB=1 -DHAVE_PLPLOT=1 -DHAVE_PYTHON=1 -D_LARGEFILE_SOURCE=1 -DHAVE_FSEEKO=1 -D_FILE_OFFSET_BITS=64 -DOBIT_THREADS_ENABLED=1 

ALL_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/include $(DEFS) \
	 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include   -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include   -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include 
ALL_CFLAGS = $(CFLAGS)  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//include/glib-2.0 -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other//lib/glib-2.0/include    -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include/plplot \
	  -I/export/users/bcotton/Software.dir/SVN/ObitInstall/other/include 

# teach how to compile
.c.o:	
	$(CC) -c $(ALL_CFLAGS) $(ALL_CPPFLAGS) $<
	mv $@.o $(LIBDIR)

# get list of all c source files (*.c) files
AllC := $(wildcard *.c)
OBJECTS := $(patsubst %.c,%.o, $(AllC))

CTARGETS := $(addprefix $(LIBDIR),$(OBJECTS))


all:  $(CTARGETS) 

# generic C compile
$(CTARGETS): $(LIBDIR)%.o: %.c $(wildcard ../include/%*.h)
	$(CC) -c $(ALL_CPPFLAGS) $(ALL_CFLAGS) $*.c
	mv ./$*.o $@

clean:
	rm -f $(CTARGETS)
	rm -f *.o
