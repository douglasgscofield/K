LATEX=		pdflatex
BIBTEX=		bibtex

STEM=		K_manual

DOC=		$(STEM).pdf

FILES=		$(STEM).tex
#FILES=		$(STEM).tex \
#			$(STEM).bib \
#			authordatedgs2.bst


all: $(DOC)

$(DOC): $(FILES)
	$(LATEX) $(STEM)
	$(BIBTEX) $(STEM)
	$(LATEX) $(STEM)
	$(LATEX) $(STEM)
	$(MAKE) clean

clean-all: clean
	rm -f $(STEM).dvi
	rm -f $(STEM).ps
	rm -f $(DOC)

clean: .FORCE
	rm -f $(STEM).aux
	rm -f $(STEM).bbl
	rm -f $(STEM).blg
	rm -f $(STEM).log
	rm -f $(STEM).lof
	rm -f $(STEM).lot
	rm -f $(STEM).toc
	rm -f $(STEM).out
	rm -f $(STEM).fls
	rm -f $(STEM).fdb_latexmk
	rm -f texput.log

.FORCE:
