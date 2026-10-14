# refs.bib here is a tracked copy; research/refs.bib is the source of truth and
# is re-copied whenever it is newer (do not edit paper/refs.bib by hand).
TEXBIN ?= /Library/TeX/texbin
export PATH := $(TEXBIN):$(PATH)
PDFLATEX = pdflatex -halt-on-error -interaction=nonstopmode
SECTIONS = $(wildcard sections/*.tex)

all: main.pdf

refs.bib: ../research/refs.bib
	cp $< $@

main.pdf: main.tex refs.bib $(SECTIONS)
	$(PDFLATEX) main
	bibtex main
	$(PDFLATEX) main
	$(PDFLATEX) main

clean:
	rm -f *.aux *.log *.bbl *.blg *.out *.toc *.fls *.fdb_latexmk *.synctex.gz *Notes.bib sections/*.aux main.pdf

.PHONY: all clean
