## Process this file with automake to create Makefile.in.

bin_SCRIPTS = ObitTalk ObitTalkServer

EXTRA_DIST = ObitTalk.in ObitTalkServer.in

# Files that should be removed, but which Automake does not know.
CLEANFILES = $(bin_SCRIPTS)


## ------------- ##
## The scripts.  ##
## ------------- ##

edit = sed \
	-e 's,@datadir\@,$(pkgdatadir),g' \
	-e 's,@PYTHON\@,$(PYTHON),g' \
	-e 's,@OBIT_PYTHONPATH\@,$(OBIT_PYTHONPATH),g'

## All the files below depend on Makefile so that they are rebuilt
## when prefix etc. changes.
## Use chmod -w to prevent people from editing the wrong file by accident.
ObitTalk: $(srcdir)/ObitTalk.in Makefile
	rm -f ObitTalk ObitTalk.tmp
	$(edit) $(srcdir)/ObitTalk.in >ObitTalk.tmp
	chmod +x ObitTalk.tmp
	chmod -w ObitTalk.tmp
	mv -f ObitTalk.tmp ObitTalk

ObitTalkServer: $(srcdir)/ObitTalkServer.in Makefile
	rm -f ObitTalkServer ObitTalkServer.tmp
	$(edit) $(srcdir)/ObitTalkServer.in >ObitTalkServer.tmp
	chmod +x ObitTalkServer.tmp
	chmod -w ObitTalkServer.tmp
	mv -f ObitTalkServer.tmp ObitTalkServer
