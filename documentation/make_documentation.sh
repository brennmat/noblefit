#!/bin/bash

# make list of NOBLEFIT tools:
octave --silent < ./make_tools_help_list.m

# run pdfLaTeX twice to get documentation in PDF format:
pdflatex noblefit_manual.tex
pdflatex noblefit_manual.tex

# open/view PDF document:
open noblefit_manual.pdf
