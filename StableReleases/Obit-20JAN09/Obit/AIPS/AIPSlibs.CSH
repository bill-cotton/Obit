#!/bin/csh
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
# Usage: . AIPSlibs.csh
#-----------------------------------------------------------------------
# AIPSlibs.csh sets environment variables to deal with the irrational
# AIPS directory structure for link libraries.  
# The variables defined here can be used to create 
# system independent .OPT files, entries are the variable names.
# The following environment variables are used:
#   LIBR     Base of AIPS architecture specific library tree
#   ARCH     AIPS architectire soecified (SOL, SUL, LINUX...)
#   FOURMASS Base of 4MASS software tree
#   OBIT     AIPS directory of Obit
#   CFITSIO  Base of cfitsio software tree
#   GLIB     directory with libglib-2.0.so
#   FFTW     directory with libfftw.a
#   RFFTW    directory with librfftw.a
#
# Variables are set depending on the architecture:
# (also $LIBRDBG versions)
# $LIBAIPSUB    $LIBR/AIPSUB/SUBLIB	
# $LIBAPLSUB	$LIBR/APLSUB/SUBLIB
# $LIBAPLNOT	$LIBR/APLNOT/SUBLIB
# $LIBAPLOOP	$LIBR/APLOOP/SUBLIB
# $LIBQOOP	$LIBR/QOOP/SUBLIB
# $LIBQSUB	$LIBR/QSUB/SUBLIB
# $LIBQNOT	$LIBR/QNOT/SUBLIB
# $LIBQPSAP	$LIBR/QPSAP/SUBLIB
# $LIBYSUB	$LIBR/YSUB/SUBLIB
# $LIBYNOT	$LIBR/YNOT/SUBLIB
# $LIBYDEV  	$LIBR/YSS/SUBLIB
# $LIBAPLARCH	The really bothersome one, generally somewhing like
#		$LIBR/APL$ARCH/SUBLIB - but not always!!!
# $LIB4MASS     4MASS specific routinesk - screwed again - cfitsio
# $LIBCFITSIO   cfitsio library
#-----------------------------------------------------------------------
#  Set predictable ones
setenv LIBAIPSUB $LIBR/AIPSUB/SUBLIB
setenv LIBAPLSUB $LIBR/APLSUB/SUBLIB
setenv LIBAPLNOT $LIBR/APLNOT/SUBLIB
setenv LIBAPLOOP $LIBR/APLOOP/SUBLIB
setenv LIBQOOP $LIBR/QOOP/SUBLIB
setenv LIBQSUB $LIBR/QSUB/SUBLIB
setenv LIBQNOT $LIBR/QNOT/SUBLIB
setenv LIBQPSAP $LIBR/QPSAP/SUBLIB
setenv LIBYSUB $LIBR/YSUB/SUBLIB
setenv LIBYNOT $LIBR/YNOT/SUBLIB
setenv LIBYDEV $LIBR/YSS/SUBLIB
setenv LIB4MASS $FOURMASS/LIB/$ARCH/libFOURMASS.a
setenv LIBOBIT  $OBIT/lib/$ARCH/libObit.a
setenv LIBCFITSIO $CFITSIO/lib/$ARCH/libcfitsio.a
setenv LIBGLIB $GLIB/glib/.libs/libglib-2.0.so
setenv LIBFFTW $LFFTW/.libs/libfftw.a
setenv LIBRFFTW $LRFFTW/.libs/librfftw.a
# debug versions of AIPS libraries
setenv LIBDBGAIPSUB $LIBRDBG/AIPSUB/SUBLIB
setenv LIBDBGAPLSUB $LIBRDBG/APLSUB/SUBLIB
setenv LIBDBGAPLNOT $LIBRDBG/APLNOT/SUBLIB
setenv LIBDBGAPLOOP $LIBRDBG/APLOOP/SUBLIB
setenv LIBDBGQOOP $LIBRDBG/QOOP/SUBLIB
setenv LIBDBGQSUB $LIBRDBG/QSUB/SUBLIB
setenv LIBDBGQNOT $LIBRDBG/QNOT/SUBLIB
setenv LIBDBGQPSAP $LIBRDBG/QPSAP/SUBLIB
setenv LIBDBGYSUB $LIBRDBG/YSUB/SUBLIB
setenv LIBDBGYNOT $LIBRDBG/YNOT/SUBLIB
setenv LIBDBGYDEV $LIBRDBG/YSS/SUBLIB
 
# Sun Solaris
if ( "$ARCH" == "SOL" )  then
  setenv LIBAPLARCH $LIBR/APLSOL/SUBLIB
  setenv LIBDBGAPLARCH $LIBRDBG/APLSOL/SUBLIB
# Sun Solaris/Ultra (the screwball one)
else if ( "$ARCH" == "SUL" )  then
  setenv LIBAPLARCH $LIBR/APLSOL/SUBLIB
  setenv LIBDBGAPLARCH $LIBRDBG/APLSOL/SUBLIB
# Linux
else if ( "$ARCH" == "LINUX" )  then
  setenv LIBAPLARCH $LIBR/APLLINUX/SUBLIB
  setenv LIBDBGAPLARCH $LIBRDBG/APLLINUX/SUBLIB
# Alpha
else if ( "$ARCH" == "ALPHA" )  then
  setenv LIBAPLARCH $LIBR/APLDEC/SUBLIB
  setenv LIBDBGAPLARCH $LIBRDBG/APLDEC/SUBLIB
# HP (defunct?)
else if ( "$ARCH" == "HP" )  then
  setenv LIBAPLARCH $LIBR/APLHP/SUBLIB
  setenv LIBDBGAPLARCH $LIBRDBG/APLHP/SUBLIB
# IBM (defunct?)
else if ( "$ARCH" == "IBM" )  then
  setenv LIBAPLARCH $LIBR/APLIBM/SUBLIB
  setenv LIBDBGAPLARCH $LIBRDBG/APLIBM/SUBLIB
# SGI 
else if ( "$ARCH" == "SGI" )  then
  setenv LIBAPLARCH $LIBR/APLSGI/SUBLIB
  setenv LIBDBGAPLARCH $LIBRDBG/APLSGI/SUBLIB
endif
