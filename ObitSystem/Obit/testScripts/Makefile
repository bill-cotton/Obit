# $Id: Makefile,v 1.6 2007/07/10 18:32:20 bcotton Exp $
#-----------------------------------------------------------------------
#;  Copyright (C) 2004-2007
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
#    Makefile to run tests of Obit/python software 
#
#------------------------------------------------------------------------
TARGETS := testObit testHGeom testFeather testUVImage testUVSub testCleanVis 

all: $(TARGETS)

testObit:
	echo " "; echo " "; echo " "; echo " ********   Run test of basic routines"
	time python testObit.py

testHGeom:
	echo " "; echo " "; echo " "; echo " ********   Run test of HGeom"
	time python testHGeom.py

testFeather:
	echo " "; echo " "; echo " "; echo " ********   Run test of Feather"
	time python testFeather.py

testUVImage:
	echo " "; echo " "; echo " "; echo " ********   Run test of UVImage"
	time python testUVImage.py

testUVSub:
	echo " "; echo " "; echo " "; echo " ********   Run test of UVSub"
	time python testUVSub.py


testCleanVis:
	echo " "; echo " "; echo " "; echo " ********   Run test of CleanVis"
	time python testCleanVis.py


testSquint:
	echo " "; echo " "; echo " "; echo " ********   Run test of Squint"
	testSquint.csh


clean:
	cd ../testIt; rm -f core.* *Out.fits	 MAImage.* Imager*Scr* HGeom*Scr* Feather*Scr* UVSub*Scr*
	cd ../testIt; rm -f CleanVis*Scr* *Out.uvtab Squint*.Scr
