NAME=          doc 
BIBTEXSRCS=    dejun.bib
USE_PDFLATEX= 1
PDFLATEX= `which pdflatex`
BIBTEX= `which bibtex`
PDFLATEX_FLAGS=-shell-escape
VIEWPDF= evince
#USE_DVIPDFM=   0   # not clear I want this -- goes straight to pdf
				   # other option is USE_PDFLATEX or nothing (dvips)
#OTHER_PDF =    $(wildcard figures/*.pdf figures/*/*.pdf figures/*/*/*.pdf)
#OTHER_PDF =   figures/*.eps
OTHER+=        $(OTHER_PDF) *.tex
#OTHER+=        elsarticle-num.bst 
.SUFFIXES : .eps .epsi .pdf


.eps.pdf:
	epstopdf $< 
.epsi.pdf:
	epstopdf $<

include /usr/share/latex-mk/latex.gmk

remake:
	make clean 
	make

#status:
#	svn status . ../../bib_files

#commit:
#	svn commit . ../../bib_files

#update:
#	svn up . ../../bib_files

#up: update
