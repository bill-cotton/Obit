# Generic Makefile to generate postscript files from LaTeX input
#
#---------------------------------------------------

# Names of documents  .tex and .ps, .pdf files
USER := ObitTalk
SOFT := ObitTalkSoft
DOCDIR = ../../doc
TARGETS := $(USER).ps  $(USER).pdf $(SOFT).ps  $(SOFT).pdf

#
# need this for the shell command syntax
# May need to change these
BIBTEX= bibtex
LATEX = latex
DVIPS = dvips -t letter
DVIPDF = dvipdf 

#
# rules for making postscript from tex files
%.dvi: %.tex          # run several times to toc etc.
	$(LATEX) $<
	$(LATEX) $<
	$(LATEX) $<

%.ps: %.dvi
	$(DVIPS) -o $*.ps $<

%.pdf: %.dvi 
	$(DVIPDF) $*.dvi $*.pdf

# rule to make bibliography
%.bbl: %.bib
	$(LATEX)  $(NAME)   # make sure references up to date.
	$(BIBTEX) $(NAME)
	$(LATEX)  $(NAME)   # make sure references up to date.
	$(BIBTEX) $(NAME)
	$(LATEX)  $(NAME)   # make sure references up to date.
	$(BIBTEX) $(NAME)

all: $(TARGETS)

install: $(DOCDIR)
	cp *.pdf $(DOCDIR)

$(DOCDIR):
	mkdir $(DOCDIR)


# clean up derived files
clean: 
	rm -f $(TARGETS)
	rm -f *.dvi *.lo* *.to* *.aux 






