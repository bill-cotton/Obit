#!/bin/sh
# Script for clean up after installing 3rd party software for Obit
# Base address
BASE=`pwd`; export BASE
rm -f -r cfitsio
rm -f -r glib-2.2.0
rm -f -r pkgconfig-0.14.0
rm -f -r fftw-2.1.3
rm -f -r fftw-3.1.2
rm -f -r gsl-1.6
rm -f -r zlib-1.2.3  
rm -f -r openmotif-2.3.0
rm -f -r Python-2.5.1
rm -f -r w3c-libwww-5.4.0
rm -f -r xmlrpc-c-1.06.18
rm -f -r curl-7.17.0
rm -f -r plplot-5.8.0
# Clean out libraries
rm -f -r bin include info lib man share
