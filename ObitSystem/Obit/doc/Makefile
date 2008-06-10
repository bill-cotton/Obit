# $Id: Makefile,v 1.1.1.1 2004/07/19 16:42:44 bcotton Exp $
# Generic Makefile to generate postscript files from Tex or LaTeX input
# Uses presence of "\begin{document}" string to indicate LaTeX
# If only the base of the tex file name is given then that postscript file
# is always built, e.g. for myDoc.tex
# % make myDoc
# With no argument all out of date .tex files are converted to postscript.
#
#---------------------------------------------------
#	 Process Tex/LaTeX files
# get list of all document files (*.tex) files
AllDoc = OBITdoc.ps
#
# need this for the shell command syntax
# May need to change these
LATEX = latex
DVIPS = dvips
DOXYGEN = doxygen
#
# Use doxygen to generate software documentation.
Doxygen = Doxygen

# built postscript files from .tex
all: $(AllDoc) $(Doxygen)


OBITdoc.ps: OBITdoc.sty OBITdoc.tex
	$(LATEX) OBITdoc
	$(LATEX) OBITdoc
	$(LATEX) OBITdoc
	$(DVIPS) -o OBITdoc.ps OBITdoc

# always run as the dependencies are pretty fierce
Doxygen:
	cd doxygen; $(DOXYGEN) Doxyfile
#	cd doxygen/latex; $(MAKE)

# clean up derived files
clean: 
	rm -f *.dvi *.log *.ps *.toc *.aux


