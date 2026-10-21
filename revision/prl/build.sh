#!/bin/sh
# Build the PRL Letter (with End Matter), the Supplemental Material and the reply to the referees.
# supp.tex resolves references to the main text through xr-hyper, so main must be compiled first.
set -e
cd "$(dirname "$0")"
pdflatex -interaction=nonstopmode main.tex >/dev/null; bibtex main >/dev/null; pdflatex -interaction=nonstopmode main.tex >/dev/null; pdflatex -interaction=nonstopmode main.tex >/dev/null
pdflatex -interaction=nonstopmode supp.tex >/dev/null; bibtex supp >/dev/null; pdflatex -interaction=nonstopmode supp.tex >/dev/null; pdflatex -interaction=nonstopmode supp.tex >/dev/null
pdflatex -interaction=nonstopmode response.tex >/dev/null; pdflatex -interaction=nonstopmode response.tex >/dev/null
echo "built: main.pdf supp.pdf response.pdf"
