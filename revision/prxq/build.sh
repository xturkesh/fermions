#!/bin/bash
# Builds the PRX Quantum resubmission package: manuscript, marked-up diff vs the
# version submitted to the journal, and the reply/appeal letter.
set -e
cd "$(dirname "$0")"
run() { pdflatex -interaction=nonstopmode "$1.tex" > "$1.log" || { tail -30 "$1.log"; exit 1; }; }
run main; bibtex main > /dev/null; run main; bibtex main > /dev/null; run main; run main
# base: the journal version with the leftover (invisible) latexdiff macros removed
latexdiff --math-markup=whole --type=UNDERLINE ../submitted/main_submitted_clean.tex main.tex 2>/dev/null > marked_up.tex
# additions: red wavy underline; deletions: red strike-through
sed -i 's|\\providecommand{\\DIFdeltex}\[1\]{.*}|\\providecommand{\\DIFdeltex}[1]{{\\protect\\color{red}\\sout{#1}}} %DIF PREAMBLE|' marked_up.tex
sed -i 's|\\providecommand{\\DIFaddtex}\[1\]{.*}|\\providecommand{\\DIFaddtex}[1]{{\\protect\\color{red}\\uwave{#1}}} %DIF PREAMBLE|' marked_up.tex
run marked_up; bibtex marked_up > /dev/null; run marked_up; bibtex marked_up > /dev/null; run marked_up; run marked_up
run response; run response  # Fig. R1 from ../checks/small_N_check.py
grep -h "undefined" main.log marked_up.log response.log | grep -v "There were" || true
