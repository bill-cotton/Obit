# $Id: Makefile.in,v 1.6 2005/12/06 18:40:43 bcotton Exp $
#-----------------------------------------------------------------------
#;  Copyright (C) 2002-2005
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

CPPFLAGS = 
CFLAGS = -g -O2 -fPIC
LDFLAGS =   -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib  -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib
CC = /usr/bin/gcc

# xmlrpc client stuff
CLIENT_CPPFLAGS =   -I/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/include  -I/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/include
CLIENT_CFLAGS =  -I/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/include
CLIENT_LDFLAGS =     -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib 
CLIENT_LIBS =  -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib  -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok  -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib  -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib -lcurl -L/usr/kerberos/lib -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib/lib -lxmlrpc_client -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok -Wl,--rpath -Wl,/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib     -lz   \
         -lxmlrpc_client -lcurl -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib  -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok  -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib  -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib -lcurl -L/usr/kerberos/lib -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib/lib -lxmlrpc_client -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok -Wl,--rpath -Wl,/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib      

ALL_CPPFLAGS = $(CPPFLAGS) -I$(top_srcdir)/include  -I/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/include \
	 -I/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/ObitSystem/Obit/include -I/export/data_2/fits/ObitTest/test/ObitInstall/other//include/glib-2.0 -I/export/data_2/fits/ObitTest/test/ObitInstall/other//lib/glib-2.0/include   -I/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/ObitSystem/Obit/include -I/export/data_2/fits/ObitTest/test/ObitInstall/other//include/glib-2.0 -I/export/data_2/fits/ObitTest/test/ObitInstall/other//lib/glib-2.0/include      \
	$(CLIENT_CPPFLAGS) -DPACKAGE_NAME=\"ObitSD\" -DPACKAGE_TARNAME=\"obitsd\" -DPACKAGE_VERSION=\"1.3\" -DPACKAGE_STRING=\"ObitSD\ 1.3\" -DPACKAGE_BUGREPORT=\"bcotton@nrao.edu\" -DPACKAGE=\"obitsd\" -DVERSION=\"1.3\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIBCFITSIO=1 -DHAVE_CFITSIO=1 -DHAVE__TMP_DUMMY1_FFTW3_H=1 -DHAVE_FFTW3=1 -DHAVE__TMP_DUMMY3_FFTW_H=1 -DHAVE_LIBOBIT=1 -DHAVE_OBIT=1 -DHAVE__TMP_DUMMY1_GSL_H=1 -DHAVE_GSL=1 -DHAVE__TMP_DUMMY1_XMLRPC_H=1 -DHAVE_XMLRPC=1 -DHAVE__TMP_DUMMY1_ZLIB_H=1 -DHAVE_ZLIB=1 -DHAVE_PLPLOT=1 -D_LARGEFILE_SOURCE=1 -DHAVE_FSEEKO=1 -D_FILE_OFFSET_BITS=64 
ALL_CFLAGS = $(CFLAGS) -I/export/data_2/fits/ObitTest/test/ObitInstall/other//include/glib-2.0 -I/export/data_2/fits/ObitTest/test/ObitInstall/other//lib/glib-2.0/include     -I/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/include $(CLIENT_CFLAGS) \
	
ALL_LDFLAGS = $(LDFLAGS)  -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/ObitSystem/Obit/lib/ -L/export/data_2/fits/ObitTest/test/ObitInstall/other//lib -lglib-2.0    -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib \
	  -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib  $(CLIENT_LDFLAGS)

LIBS = ../lib/$(ARCH)/libObitSD.a -lObit  -lcfitsio   \
	-L/export/data_2/fits/ObitTest/test/ObitInstall/other//lib -lglib-2.0   -lgsl -lgslcblas   $(CLIENT_LIBS)   -L/export/data_2/fits/ObitTest/Fix64/ObitInstallWork/other/lib -L/usr/lib/gcc/i386-redhat-linux/3.4.6 -L/usr/lib/gcc/i386-redhat-linux/3.4.6/../../.. -lfrtbegin -lg2c -lm -lgcc_s 


# get list of all c source files (*.c) files
AllC := $(wildcard *.c)
OBJECTS := $(wildcard ../lib/LINUX/*.o)
TARGETS := $(patsubst %.c,%, $(AllC))
TARGETS := Simulator Analyzer ImageOTF DCR2OTF EditOTF SP2OTF \
	GainCalOTF AtmCorOTF Soln2Cal ResidCalOTF SplitOTF ConcatOTF \
	CleanOTF
TARGETS := testSelfCal Analyzer 

# get list of all c header files (*.h) files
HEADERS := $(wildcard ../include/*.h)

all: $(TARGETS)

# generic C compile/link
# $(TARGETS): $*.c $(HEADERS) ../lib/linux/libObit.a
# 	$(CC) $(CFLAGS)  %.c -o $< $(LIB)

testSelfCal: testSelfCal.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

ConcatOTF: ConcatOTF.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

SplitOTF: SplitOTF.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

OTF2Cube: OTF2Cube.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

Simulator: Simulator.c $(HEADERS) $(OBJECTS)
	$(CC) $(CFLAGS)  Simulator.c -o Simulator  $(LIB)

PASim2OTF: PASim2OTF.c $(HEADERS) $(OBJECTS)
#noplot	$(CC) $(CFLAGS)  PASim2OTF.c -o PASim2OTF $(LIB)
	$(CC) $(CFLAGS)  PASim2OTF.c -c
	$(FC) $(FFLAGS)  PASim2OTF.o -o PASim2OTF $(LIB)

Analyzer: Analyzer.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

ImageOTF: ImageOTF.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

CleanOTF: CleanOTF.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

EditOTF: EditOTF.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

GainCalOTF: GainCalOTF.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

AtmCorOTF: AtmCorOTF.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

ResidCalOTF: ResidCalOTF.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)
#plot	$(CC) $(CFLAGS)  ResidCalOTF.c -c
#plot	$(FC) $(FFLAGS)  ResidCalOTF.o -o ResidCalOTF $(LIB)

DCR2OTF: DCR2OTF.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

SP2OTF: SP2OTF.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

X2OTF: X2OTF.c $(HEADERS) $(OBJECTS)

Soln2Cal: Soln2Cal.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

FitTip: FitTip.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

AvgCal: AvgCal.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

MakeImage: MakeImage.c MakeImageData.h $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

testFilt: testFilt.c $(HEADERS) $(OBJECTS)
	$(CC) $(CFLAGS)  testFilt.c -c 
	$(FC) $(FFLAGS)  testFilt.o -o testFilt  $(LIB)

testFilt2: testFilt2.c $(HEADERS) $(OBJECTS)
	$(CC) $(CFLAGS)  testFilt2.c -c 
	$(FC) $(FFLAGS)  testFilt2.o -o testFilt2  $(LIB)

testFitCal: testFitCal.c $(HEADERS) $(OBJECTS)
	$(CC) $(CFLAGS)  testFitCal.c -c 
	$(FC) $(FFLAGS)  testFitCal.o -o testFitCal  $(LIB)

test:	 test.c $(HEADERS) $(OBJECTS)
	$(CC) $(CFLAGS)  test.c -c 
	$(FC) $(FFLAGS)  test.o -o test  ../lib/LINUX/libPennArray.a $(LIB)

testClean: testClean.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

testUVSub: testUVSub.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

test: test.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

test2: test2.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

test3: test3.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

testPre: testPre.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

AIPS2FITS: AIPS2FITS.c $(HEADERS) $(OBJECTS)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(ALL_LDFLAGS) $< -o $@ $(LIBS)

# special case, junk - ignore
# junk:

clean:
	rm -f $(TARGETS)
	rm -f *.o

