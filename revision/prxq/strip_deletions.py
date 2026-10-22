"""Remove the deleted blocks of a latexdiff output, so that only the new text remains,
with the additions marked (red, see red_preamble.tex)."""
import re, sys
src = open(sys.argv[1]).read()
k = src.index('\\begin{document}')
head, s = src[:k], src[k:]          # only the body: the preamble defines the DIF macros
s, n = re.subn(r'\\DIFdelbegin(FL)?\b.*?\\DIFdelend(FL)?\b ?', '', s, flags=re.S)
s = re.sub(r'^%DIF < .*\n', '', s, flags=re.M)
open(sys.argv[1], 'w').write(head + s)
print(f'removed {n} deleted blocks', file=sys.stderr)
