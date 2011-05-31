# Python/Obit build utillity
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2005
#  Associated Universities, Inc. Washington DC, USA.
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public
#  License along with this program; if not, write to the Free
#  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
#  MA 02139, USA.
#
#  Correspondence concerning this software should be addressed as follows:
#         Internet email: bcotton@nrao.edu.
#         Postal address: William Cotton
#                         National Radio Astronomy Observatory
#                         520 Edgemont Road
#                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------

# Create setup.py file from setupdata.py file with
# compiler/linker instructions (see Makefile)
# to build and install interface module, e.g.:
# CFLAGS='-g -O2 -fPIC  -I/usr/include/glib-2.0 ... '
# CPPFLAGS=' -I../include -I/usr/local/include/python2.3  ...'
# LDFLAGS=' -L/home/bcotton/opt/xmlrpc/lib ... '
# LIBS='../lib/LINUX/libObit.a -lcfitsio -lm  ... '
#
# Run setup.py as:
# python setup.py build install --install-lib=.
import os

# Read details left by Makefile
import setupdata

# Variables needed
packageName    = ''
packageVer     = ''
compileArgs    = []
incDirs        = []
libDirs        = []
runtimeLibDirs = []
libs           = []

# Parse input from  setupdata
tt = setupdata.CFLAGS
# Cleanup line
tt=tt.replace("\n","_"); tt=tt.replace("\\ ","_");
t = tt.split()
for x in t:
    if x[0:2]=='-I':
        incDirs.append(x[2:])
        
tt = setupdata.CPPFLAGS
tt=tt.replace("\n","_"); tt=tt.replace("\\ ","_");
t = tt.split()
for x in t:
    if x[0:2]=='-I':
        incDirs.append(x[2:])
    elif x[0:7]=='-DHELLO':
        pass
    elif x[0:2]=='-D':
        compileArgs.append(x)
    if x[0:16]=='-DPACKAGE_NAME=\"':
        packageName = x[16:len(x)-1]
    if x[0:19]=='-DPACKAGE_VERSION=\"':
        packageVer = x[19:len(x)-1]

tt = setupdata.LIBS
tt=tt.replace("\n","_"); tt=tt.replace("\\ ","_");
t = tt.split()
# assume package library name
libs = [packageName]
for x in t:
    if x[0:2]=='-l':
        libs.append(x[2:])
    elif x[0:3]=='-Wl':   # Ignore linker options
        pass
    elif x[0:2]=='-L':
        libDirs.append(x[2:])
    else:   # better be a path
        libDirs.append(os.path.dirname(x))
        libs.append(os.path.basename(x)[3:].split('.')[0])

tt = setupdata.LDFLAGS
tt=tt.replace("\n","_"); tt=tt.replace("\\ ","_");
t = tt.split()
for x in t:
    if x[0:2]=='-L':
        libDirs.append(x[2:])
    elif x[0:11]=='-Wl,-rpath,':
         runtimeLibDirs.append( x[11:] )

# Dump it out
outfile = file("setup.py","w")
outfile.write('from distutils.core import setup, Extension'+os.linesep)
outfile.write('setup( name=\"'+packageName+'\", version=\"'+packageVer+'\",'+os.linesep)
outfile.write('       ext_modules=[Extension(\"'+packageName+'\",'+os.linesep)
outfile.write('                              [\''+packageName+'_wrap.c\'],'+os.linesep)
outfile.write('                              extra_compile_args='+str(compileArgs)+','+os.linesep)
outfile.write('                              library_dirs='+str(libDirs)+','+os.linesep)
outfile.write('                              libraries='+str(libs)+','+os.linesep)
outfile.write('                              runtime_library_dirs='+str(runtimeLibDirs)+')],'+os.linesep)
outfile.write('       include_dirs='+str(incDirs)+os.linesep)
outfile.write(')')
