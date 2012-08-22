K_manual.pdf: K_manual.tex K_manual.bib
	pdflatex K_manual.tex
	bibtex K_manual
	pdflatex K_manual.tex
	pdflatex K_manual.tex

clean:
	rm -f K_manual.pdf K_manual.dvi K_manual.ps K_manual.aux K_manual.bbl K_manual.blg K_manual.lof K_manual.log K_manual.lot K_manual.tcp K_manual.toc K_manual.tps
