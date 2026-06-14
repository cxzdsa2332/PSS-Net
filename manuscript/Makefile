# Build the PSS-Net manuscript.
#   make        -> compile main.pdf (latexmk, handles bibtex passes)
#   make clean  -> remove aux files
#   make purge  -> also remove main.pdf

MAIN = main

.PHONY: all clean purge

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex preamble.tex sections/*.tex ../ref/references.bib
	latexmk -pdf -bibtex $(MAIN).tex

clean:
	latexmk -c $(MAIN).tex

purge:
	latexmk -C $(MAIN).tex
