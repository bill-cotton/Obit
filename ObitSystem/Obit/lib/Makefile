# $Id: Makefile,v 1.2 2005/11/03 14:14:04 bcotton Exp $
#-----------------------------------------------------------------------
#;  Copyright (C) 2002,2003
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
#    Makefile to build library
#    Obit: Software for the recently deceased
#
#   Expects RANLIB to be set by either environment variable  or Make argument
#     RANLIB is the library "ranomizer".
#
#------------------------------------------------------------------------
# targets to build
TARGETS =  libObit.a

# list of object modules
OBJECTS := $(wildcard *.o)

all:  $(TARGETS)

#  build Obit library
libObit.a: ${OBJECTS}
	ar rv libObit.a ${OBJECTS}
	${RANLIB} libObit.a

clean:
	rm -f $(TARGETS)


