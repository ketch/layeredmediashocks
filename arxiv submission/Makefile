PARTS = stego.tex macros.tex 

stego.pdf: stego.ps
	ps2pdf stego.ps

stego.ps: stego.dvi
	dvips -o stego.ps stego.dvi

stego.dvi: stego.bbl ${PARTS}
	latex stego

stego.bbl: stego.bib stego.aux
	bibtex stego

stego.aux: stego.tex ${PARTS} force
	latex stego

force:

clean:	
	rm -f *.ps *.dvi *.pdf *.aux *.bbl *.log *.blg
