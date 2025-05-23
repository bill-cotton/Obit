# $Id$
#-----------------------------------------------------------------------
#;  Copyright (C) 2002-2025
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
#
#    Obit: Merx mollis mortibus nuper
#
# All compiles are one or two level down
#CC = gcc # default
#----------------------------------------------------
top_srcdir = .

CPPFLAGS = 
CFLAGS = -g -O3 -fPIC -msse  -Wall -fno-strict-aliasing
LDFLAGS = 
OBIT = $(top_srcdir)
CC = gcc

# RANLIB = library "randomizer", ranlib if needed, else blank
RANLIB = ranlib

# LIBS = locations of libraries and libraries to link programs
LIBS =  -lpthread -lgthread-2.0 -lglib-2.0  

# DOCS = LaTeX Documentation file
DOCS = $(top_srcdir)/doc/obitdoc.tex

# names of archive files
ARCHIVE = obitBack2.643M
DISTRIB = obit1.0
# Name of this directory
DIRN = Obit

#------------------------------------------------------------------------
#TARGETS = versionupdate cfitsioupdate xmlrpcupdate srcupdate libupdate \
#	pythonupdate taskupdate
TARGETS = versionupdate srcupdate libupdate pythonupdate taskupdate

all:  $(TARGETS)

# update the Obit version (svn revision)
versionupdate:
	python share/scripts/getVersion.py ${OBIT}

# update dummy_cfitsio directory
cfitsioupdate:
	cd dummy_cfitsio; $(MAKE)

# update dummy_xmlrpc directory
xmlrpcupdate:
	cd dummy_xmlrpc; $(MAKE)

# update tables from documentation directory
tableupdate: 
	perl bin/ObitTables.pl

# update documentation directory
docupdate: 
	cd doc; $(MAKE)

# update source/object directory
srcupdate: 
	cd src; $(MAKE)

# update library directory
libupdate: 
	cd lib; $(MAKE) RANLIB="$(RANLIB)"

# update test software directory
testupdate: 
	cd test; $(MAKE)

# update task software directory
taskupdate: 
	cd tasks; $(MAKE)

# update work directory
work: 
	cd src/work; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" LIB="$(LIB)"

# update python directory
pythonupdate: 
	cd python; $(MAKE)

# update from cvs repository
cvsupdate:
	cvs -d $(WDCCVS) update

# commit changes to cvs repository
cvscommit:
	cvs -d $(WDCCVS) commit

# Install
install:
	echo "Nothing to install"

test:
	echo "Testing not yet implemented"

# Backup to USB drives
backup:
	rsync -r -u -p -t src include tasks python TDF share doc m4 Makefile *.log *.text  /mnt/RAID_WDC2/ObitBackup/
	rsync -r -u -p -t ../ObitView /mnt/RAID_WDC2/ObitBackup/ObitView/
	rsync -r -u -p -t ../ObitTalk /mnt/RAID_WDC2/ObitBackup/ObitTalk/
	rsync -r -u -p -t src include tasks python TDF share doc m4 Makefile *.log *.text  /mnt/RAID_2024/ObitBackup/
	rsync -r -u -p -t ../ObitView /mnt/RAID_2024/ObitBackup/ObitView/
	rsync -r -u -p -t ../ObitTalk /mnt/RAID_2024/ObitBackup/ObitTalk/

# make gzipped tar backup of directory
backuptar:	
	cd ..;rm -f $(ARCHIVE).tgz
	cd ..;gtar czvf $(ARCHIVE).tgz \
		$(DIRN)/README $(DIRN)/LICENSE           \
		$(DIRN)/Makefile*  $(DIRN)/src/Makefile*  \
		$(DIRN)/tasks/Makefile* $(DIRN)/lib/Makefile*   \
		$(DIRN)/include/*.h $(DIRN)/src/*.c      \
		$(DIRN)/tasks/*.c  $(DIRN)/tasks/*.doc $(DIRN)/TDF/*.TDF  \
		$(DIRN)/m4/*.m4 \
		$(DIRN)/bin/*.pl \
		$(DIRN)/doc/*.tex $(DIRN)/doc/*.sty  $(DIRN)/doc/doxygen   \
		$(DIRN)/python/*.inc $(DIRN)/python/*.swig   \
		$(DIRN)/python/*.py  $(DIRN)/python/*.c \
		$(DIRN)/python/Makefile* \
		$(DIRN)/dummy_cfitsio/*.c $(DIRN)/dummy_cfitsio/*.h \
		$(DIRN)/dummy_cfitsio/Make* $(DIRN)/dummy_cfitsio/README \
		$(DIRN)/dummy_xmlrpc/*.c $(DIRN)/dummy_xmlrpc/*.h \
		$(DIRN)/dummy_xmlrpc/Make* $(DIRN)/dummy_xmlrpc/README \
		$(DIRN)/config* $(DIRN)/install-sh \
		$(DIRN)/mkinstalldirs $(DIRN)/missing $(DIRN)/m4/* \
		#$(DIRN)/CVS/* $(DIRN)/doc/CVS/* $(DIRN)/src/CVS/* \
		#$(DIRN)/include/CVS/* $(DIRN)/lib/CVS/*

# make gzipped tar distribution
distrib:	
	cd ..;rm -f $(DISTRIB).tgz
	cd ..;gtar czvf $(DISTRIB).tgz \
		$(DIRN)/README $(DIRN)/LICENSE           \
		$(DIRN)/Makefile*  $(DIRN)/src/Makefile*  \
		$(DIRN)/lib/Makefile $(DIRN)/scripts/*   \
		$(DIRN)/include/*.h  $(DIRN)/src/*.c $ \
		$(DIRN)/bin/*.pl  $(DIRN)/data/* \
		$(DIRN)/doc/Makefile $(DIRN)/doc/*.tex   \
		$(DIRN)/doc/*.sty  $(DIRN)/doc/doxygen   \
		$(DIRN)/AIPS/Makefile $(DIRN)/AIPS/*.FOR \
		$(DIRN)/AIPS/*.HLP $(DIRN)/AIPS/*.OPT  \
		$(DIRN)/AIPS/*.CSH   \
		$(DIRN)/python/*.inc $(DIRN)/python/*.swig   \
		$(DIRN)/python/*.py  $(DIRN)/python/*.c \
		$(DIRN)/python/Makefile* \
		$(DIRN)/test/Makefile* $(DIRN)/test/*.h   \
		$(DIRN)/test/*.c $(DIRN)/test/*.doc \
		$(DIRN)/testScripts/Makefile  $(DIRN)/testScripts/*.py  \
		$(DIRN)/tasks/Makefile* $(DIRN)/tasks/*.h   \
		$(DIRN)/tasks/*.c  $(DIRN)/tasks/*.doc $(DIRN)/TDF/*.TDF \
		$(DIRN)/dummy_cfitsio/*.c $(DIRN)/dummy_cfitsio/*.h \
		$(DIRN)/dummy_cfitsio/Make* $(DIRN)/dummy_cfitsio/README \
		$(DIRN)/dummy_xmlrpc/*.c $(DIRN)/dummy_xmlrpc/*.h \
		$(DIRN)/dummy_xmlrpc/Make* $(DIRN)/dummy_xmlrpc/README \
		$(DIRN)/config* $(DIRN)/aclocal.m4 $(DIRN)/install-sh \
		$(DIRN)/mkinstalldirs $(DIRN)/missing $(DIRN)/m4/* 

# make gzipped tar distribution of only software
copy:	
	cd ..;rm -f $(DISTRIB)Src.tgz
	cd ..;gtar czvf $(DISTRIB)Src.tgz \
		$(DIRN)/README $(DIRN)/LICENSE           \
		$(DIRN)/Makefile*  $(DIRN)/src/Makefile*  \
		$(DIRN)/lib/Makefile* $(DIRN)/scripts/*   \
		$(DIRN)/include/*.h $(DIRN)/src/*.c \
		$(DIRN)/bin/*.pl  \
		$(DIRN)/doc/Makefile* $(DIRN)/doc/*.tex   \
		$(DIRN)/doc/*.sty   \
		$(DIRN)/AIPS/Makefile $(DIRN)/AIPS/*.FOR \
		$(DIRN)/AIPS/*.HLP $(DIRN)/AIPS/*.OPT $(DIRN)/AIPS/*.CSH    \
		$(DIRN)/python/*.inc $(DIRN)/python/*.swig   \
		$(DIRN)/python/*.py $(DIRN)/python/*.c  \
		$(DIRN)/python/Makefile* \
		$(DIRN)/test/Makefile $(DIRN)/test/*.h   \
		$(DIRN)/test/*.c $(DIRN)/test/*.doc      \
		$(DIRN)/testScripts/Makefile  $(DIRN)/testScripts/*.py  \
		$(DIRN)/dummy_cfitsio/*.c $(DIRN)/dummy_cfitsio/*.h \
		$(DIRN)/dummy_cfitsio/Make* $(DIRN)/dummy_cfitsio/README \
		$(DIRN)/dummy_xmlrpc/*.c $(DIRN)/dummy_xmlrpc/*.h \
		$(DIRN)/dummy_xmlrpc/Make* $(DIRN)/dummy_xmlrpc/README \
		$(DIRN)/tasks/Makefile* $(DIRN)/tasks/*.h   \
		$(DIRN)/tasks/*.c  $(DIRN)/tasks/*.doc $(DIRN)/TDF/*.TDF 

clean:
	cd src;  $(MAKE) clean
	cd lib;  $(MAKE) clean
	cd tasks; $(MAKE) clean
	cd python; $(MAKE) clean
