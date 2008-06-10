## Process this file with automake to create Makefile.in.

check_SCRIPTS = ObitTest

TESTS = mandl.py imean.py jmfit.py \
	$(top_srcdir)/python/MinimalMatch.py \
	$(top_srcdir)/python/Task.py \
	$(top_srcdir)/python/AIPSTask.py
TESTS_ENVIRONMENT = ./ObitTest

EXTRA_DIST = $(TESTS) template.py ObitTest.in

# Files that should be removed, but which Automake does not know.
CLEANFILES = $(check_SCRIPTS)

## ------------- ##
## The scripts.  ##
## ------------- ##

edit = sed \
	-e 's,@top_srcdir\@,$(top_srcdir),g' \
	-e 's,@PYTHON\@,$(PYTHON),g' \
	-e 's,@OBIT_PYTHONPATH\@,$(OBIT_PYTHONPATH),g'

## All the files below depend on Makefile so that they are rebuilt
## when prefix etc. changes.
## Use chmod -w to prevent people from editing the wrong file by accident.
ObitTest: $(srcdir)/ObitTest.in Makefile
	rm -f ObitTest ObitTest.tmp
	$(edit) $(srcdir)/ObitTest.in >ObitTest.tmp
	chmod +x ObitTest.tmp
	chmod -w ObitTest.tmp
	mv -f ObitTest.tmp ObitTest
