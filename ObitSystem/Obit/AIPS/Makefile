# $Id: Makefile,v 1.1.1.1 2004/07/19 16:42:44 bcotton Exp $
#-----------------------------------------------------------------------
#;  Copyright (C) 2003
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
#;  You should have received a copy of the GNU General Public
#;  License along with this program; if not, write to the Free
#;  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
#;  MA 02139, USA.
#;
#;  Correspondence concerning 4MASSAIPS should be addressed as follows:
#;         Internet email: bcotton@nrao.edu
#;         Postal address: W. D. Cotton
#;                         National Radio Astronomy Observatory
#;                         520 Edgemont Road
#;                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------
#    Makefile to build application programs in work directory
#    Obit/AIPS
#
#   Environment variable ARCH must be set prior to call
#   PP, FC, FFLAGS, LINK may be set by either environment variable 
#   or Make argument
#     PP is the preprocessor	
#     FC is the compiler to use
#     FFLAGS is compile options
#     LINK is the AIPS link procedure
PP = PP
FC = FC
# DEBUG is the debugger flag, 
#-------------------------------------
DEBUG =                 # no debug
DEBUG = DEBUG NOOPT NOPURGE   # debug

# AIPS Compiler flags
#-------------------------------------
FFLAGS = $(DEBUG)
ifeq ($(ARCH),SGI)
   FFLAGS = $(DEBUG) -Ofast -fullwarn
endif

# Obit libraries
# For Linux (lusus)
#LIB =  -L$(CFITSIO)/lib -L../lib/$(ARCH) -L/usr/lib  -L/usr/GLib/lib\
#	-L$(FFTW)/rfftw/.libs -L$(FFTW)/fftw/.libs \
#	 -lObit -lcfitsio -lm -lnsl -lfftw -lrfftw -lnsl -lglib-2.0
# For Linux (gollum)
LIB =  -L$(CFITSIO)/lib -L../lib/$(ARCH) -L/usr/lib \
	 -lObit -lcfitsio -lm -lnsl -lfftw -lrfftw -lnsl -lglib-2.0

LINK  = LINK
PP = PP
FC = FC
CC = gcc
INCDIR = -I$(CFITSIO)/include -I../include
CFLAGS = -c $(INCDIR)
#
#------------------------------------------------------------------------
# Subroutines used
AllSub = INPUT.FOR GTPARM.FOR

FTARGETS = $(patsubst %.FOR,%.o, $(AllSub))

# get list of all c source files (*.c) files
AllC := $(wildcard *.c)
AllC = ZOINTD.c
AllCObj := $(patsubst %.c,%.o, $(AllC))
# get list of all c header files (*.h) files
AllH := $(wildcard *.h)
CTARGETS := $(patsubst %.c,%, $(AllC))
CTARGETS := $(AllCObj)

# List of target programs in work directory
TARGETS = OBTST.EXE

all:  $(FTARGETS) $(CTARGETS) $(TARGETS)

one:  $(FTARGETS) $(CTARGETS) OBTST.EXE

# generic C compile
$(CTARGETS): %.o: %.c $(AllH)
	$(CC) $(CFLAGS)  $*.c


# generic Fortran compile
$(FTARGETS): %.o: %.FOR
	$(PP) $*.FOR
	$(FC) $(FFLAGS)  $*

# generic compile and link
$(TARGETS): %.EXE: %.FOR %.OPT  $(AllSub)
	$(PP) $*.FOR
	$(FC) $(FFLAGS)  $*
	$(LINK) $*.o $(FFLAGS)  $*.OPT 


test:	test.f 
	/usr/bin/g77 -c -g $(FFLAGS)  test.f
	/usr/bin/g77 -g test.o -o test $(FFLAGS)

# make gzipped tar backup of directory 
backup:	
	cd ..;rm -f 4MASSwork.tgz
	cd ..;gtar czvf 4MASSwork.tgz work/Makefile  \
		 work/*.FOR work/*.OPT work/*.HLP \
		work/changed work/*.c work/*.h work/*.def

# copy backup to gorilla_2
copy:	backup
	scp ../4MASSwork.tgz ${GORILLA_2}/4MASS/

clean:
	rm -f $(TARGETS) *.o





