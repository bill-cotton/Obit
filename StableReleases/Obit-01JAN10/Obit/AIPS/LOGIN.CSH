#!/bin/csh
#-----------------------------------------------------------------------
#;  Copyright (C) 2003,2004
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
#;  Correspondence concerning Obit should be addressed as follows:
#;         Internet email: bcotton@nrao.edu
#;         Postal address: W. D. Cotton
#;                         National Radio Astronomy Observatory
#;                         520 Edgemont Road
#;                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------
# Usage: source LOGIN.CSH
#-----------------------------------------------------------------------
# Set environment variables for using the Obit AIPS enhancements
# The basic AIPS setup should be done first, 
# especially the variable $ARCH.
# Knows about Gollum, NRAOCV, NRAOAOC, other sites must be added separately
# or FOURMASS, CFITSIO previously defined.
#-----------------------------------------------------------------------
#
# Base address of OBIT software directories
# Set OBIT if not already defined
if (! $?OBIT) then
  if ("$SITE" == "GOLLUM") then
#        Bill Cotton home
    setenv OBIT       /home/bcotton/Software.dir/CVS/Obit
  else if ("$SITE" == "NRAOAOC") then
#        NRAO AOC
  else if ("$SITE" == "COAOARN") then
#        NRAO AOC 
  else if ("$SITE" == "NRAOCV") then
#        NRAO CV
  else if ("$SITE" == "VCOARN") then
#        NRAO CV LINUX
    setenv OBIT       /export/users/bcotton/Software.dir/Obit
  else
#   Unknown site, need prior definition
    echo "LOGIN: You must define OBIT first"
  endif
endif
#
#  variable for AIPS VERSION adverb for executables and help files
setenv OBITAIPS   $OBIT/AIPS
#
# Set GLIB if not already defined
if (! $?GLIB) then
  if ("$SITE" == "GOLLUM") then
#        Bill Cotton home
    setenv OBIT       /home/bcotton/Software.dir/CVS/Obit
    setenv GLIB       /home/bcotton/OtherSoftware.dir/glib-2.0.7/glib
    setenv LFFTW      /home/bcotton/OtherSoftware.dir/FFTW/fftw-2.1.3/fftw
    setenv LRFFTW     /home/bcotton/OtherSoftware.dir/FFTW/fftw-2.1.3/rfftw
  else if ("$SITE" == "NRAOAOC") then
#        NRAO AOC
  else if ("$SITE" == "COAOARN") then
#        NRAO AOC 
  else if ("$SITE" == "NRAOCV") then
#        NRAO CV
  else if ("$SITE" == "VCOARN") then
#        NRAO CV LINUX
    setenv OBIT       /net/lusus/export/data_1/Software.dir/CVS/Obit
    setenv GLIB       /net/lusus/export/data_1/OtherSoftware.dir/glib-2.2.0
    setenv LFFTW      /net/lusus/export/data_1/OtherSoftware.dir/fftw-2.1.3/fftw
    setenv LRFFTW     /net/lusus/export/data_1/OtherSoftware.dir/fftw-2.1.3/rfftw
  else
#   Unknown site, need prior definition
    echo "LOGIN: You must define GLIB, FFTW first"
  endif
endif
#
# Set FFTW if not already defined
if (! $?FFTW) then
  if ("$SITE" == "GOLLUM") then
#        Bill Cotton home
    setenv FFTW      /home/bcotton/OtherSoftware.dir/FFTW/fftw-2.1.3/
  else if ("$SITE" == "NRAOAOC") then
#        NRAO AOC
  else if ("$SITE" == "COAOARN") then
#        NRAO AOC 
  else if ("$SITE" == "NRAOCV") then
#        NRAO CV
  else if ("$SITE" == "VCOARN") then
#        NRAO CV LINUX
    setenv FFTW      /net/lusus/export/data_1/OtherSoftware.dir/fftw-2.1.3
  else
#   Unknown site, need prior definition
    echo "LOGIN: You must define FFTW first"
  endif
endif
#
# Where are FFTW libraries?
setenv LFFTW      $FFTW/fftw
setenv LRFFTW     $FFTW/rfftw
#
# Base address of cfitsio
# Set CFITSIO if not already defined
if (! $?CFITSIO) then
  if ("$SITE" == "GOLLUM") then
#        Bill Cotton home
    setenv CFITSIO   /home/bcotton/OtherSoftware.dir/cfitsio2430
  else if ("$SITE" == "COAOARN") then
#        NRAO AOC
    setenv CFITSIO   /home/packages/fourmass/cfitsio
  else if ("$SITE" == "VCOARN") then
#        NRAO CV
    setenv CFITSIO   /net/lusus/export/data_1/OtherSoftware.dir/cfitsio
  else if ("$SITE" == "NRAOAOC") then
#        NRAO AOC
    setenv CFITSIO   /home/packages/fourmass/cfitsio
  else if ("$SITE" == "NRAOCV") then
#        NRAO CV
    setenv CFITSIO   /net/gorilla/gorilla_2/cfitsio
  else
#   Unknown site, need prior definition
    echo "LOGIN: You must define CFITSIO first"
  endif
endif
#
#  variable for AIPS VERSION adverb for executables and help files
setenv AIPSOBIT   $OBIT
#
# Set library variables
source $FOURMASS/bin/LOGIN.CSH
source $OBIT/AIPS/AIPSlibs.CSH
#
# FITS files directory
setenv FITS4MAS   $FOURMASS/FITS
#
# Use 4MASS includes
setenv INCS "$FOURMASS/include $INC"
