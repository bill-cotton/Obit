# $Id: Makefile,v 1.3 2005/11/04 00:08:39 bcotton Exp $
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
#    Makefile to build library
#    Obit: Merx mollis mortibus nuper
#
#   Expects RANLIB to be set by either environment variable  or Make argument
#     RANLIB is the library "ranomizer".
#
#------------------------------------------------------------------------
# targets to build
TARGETS =  libObitSD.a

# list of object modules
OBJECTS := $(wildcard *.o)

all:  $(TARGETS)

#  build Obit library
libObitSD.a: ${OBJECTS}
	ar rv libObitSD.a ${OBJECTS}
	${RANLIB} libObitSD.a

clean:
	rm -f $(TARGETS)


