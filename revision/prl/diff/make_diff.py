"""Diff of the PRL Letter against the PRX Quantum version.
old.tex: the PRX Q main text (Secs. I-IV), reordered in the order of the Letter
(main text | End Matter | text moved to the SM), so that moved text is not shown
as cut and re-added. new.tex: Letter main text + End Matter + the SM sections
that hold text moved out of the main text. Run from revision/: python3 prl/diff/make_diff.py"""
P=open('prxq/main.tex').read().split('\n')
L=lambda a,b: '\n'.join(P[a-1:b])
line=lambda n: P[n-1]
def split_after(n, upto):
    t=line(n); k=t.index(upto)+len(upto); return t[:k], t[k:].lstrip()
t296=line(296); q=t296.index('We emphasize that this is a')
s459_main, s459_em = split_after(459, 'finding quantitative agreement.')
s481_main, s481_em = split_after(481, 'free-fermion Hamiltonians.')
ack=next(i for i,l in enumerate(P,1) if l.startswith(r'\begin{acknowledgments}'))
outl=next(i for i,l in enumerate(P,1) if l.startswith(r'\subsection{Outlook}'))
NOTE=(r"\noindent\textit{Diff with respect to the PRX Quantum version: removed text in red, added text in blue. "
      r"Text moved to the End Matter or to the Supplemental Material is shown in its new place (after the main text), not as removed and re-added. "
      r"SM Secs.~S1--S4 coincide with the former Appendices, except for the canonical form, Eq.~(S6), now written without qubit permutations. "
      r"Eqs.~(8) and (12) differ from the previous version only in the final punctuation.}"+"\n\n")
old_main='\n\n'.join([
 r'\section{Introduction}', L(111,126),
 r'\section{Methods}', L(130,159),
 r'\subsection{Non-local magic and its fermionic counterpart}', L(163,181),
 r'\section{Results}', L(184,238),
 r'\subsection{Typical Gaussian states}', L(268,269), t296[q:], L(298,305), line(307), L(319,326), line(329), L(333,334),
 r'\subsection{Ground states of XY model}', L(340,358), L(360,361), L(384,385), line(394), line(396), L(400,403),
 r'\subsubsection{Exact results off criticality via corner transfer matrices}', L(420,430), L(445,450), L(453,457), s459_main, L(461,475), line(477),
 r'\subsection{Area law for Gaussian-restricted non-local magic}', s481_main,
 r'\subsection{Quench dynamics and quasiparticle picture}', L(537,551),
 r'\subsubsection{Stationary state via dephasing}', line(578), L(580,587),
 r'\subsubsection{Quasiparticle picture}', L(615,622), line(625), L(627,629),
 r'\section{Discussion}', r'\subsection{Summary}', L(655,658),
])
old_em='\n\n'.join([r'\section*{End Matter}',
 r'\subsection{Experimental observability}', L(256,264),
 L(270,295), t296[:q].rstrip(), L(308,318), s459_em, s481_em, L(483,533), L(590,591),
 r'\subsection{Random Gaussian circuit dynamics}', L(632,650)])
old_moved='\n\n'.join([r'\section*{Moved to the Supplemental Material}',
 L(330,332),
 r'\subsection{Simple example: entanglement without non-local magic}', L(243,252),
 line(362), line(359), L(363,382), L(386,393), line(394), line(395), line(396), line(398), L(406,418), L(431,444), line(451), line(461),
 L(552,575), L(594,614), line(579), line(623), line(626),
 r'\subsection{Outlook}', L(660,664), L(outl+1, ack-1)])
abs_i=next(i for i,l in enumerate(P,1) if l.startswith(r'\begin{abstract}'))
N=open('prl/main.tex').read()
pre=N[:N.index(r'\begin{abstract}')]
stubs=(r"\makeatletter\newlabel{app:proof}{{V\,A}{}{}{}{}}\newlabel{app:proof:dpp}{{V\,A\,1}{}{}{}{}}\newlabel{app:variational}{{V\,B}{}{}{}{}}"
       r"\newlabel{app:experimental}{{V\,C}{}{}{}{}}\newlabel{Sec:generaliz}{{V\,D}{}{}{}{}}\newlabel{eq:analytic_NL_weights}{{9}{}{}{}{}}\makeatother"+"\n")
new_abs=N[N.index(r'\begin{abstract}'):N.index(r'\end{abstract}')+len(r'\end{abstract}')]
body=N[N.index(r'\maketitle')+len(r'\maketitle'):N.index(r'\textbf{Note added.}')]
em=open('prl/endmatter.tex').read()
S=open('prl/supp.tex').read()
def sec(a,b):
    i=S.index(a); j=S.index(b,i); return S[i:j].strip()
moved='\n\n'.join([r'\section*{Moved to the Supplemental Material}',
 sec(r'The SYK$_2$ model is specified by the Hamiltonian', r'\section{Simple example'),
 sec(r'\section{Simple example: entanglement without non-local magic}', r'\bibliography{ref_current}'),
]).replace('M-eq:','eq:').replace('M-thm:','thm:').replace('M-fig:','fig:')
tail="\n\\bibliography{../ref_current}\n\\end{document}\n"
new_doc=pre+new_abs+"\n\\maketitle\n"+NOTE+body+"\n\\clearpage\n"+em+"\n\\clearpage\n"+moved+tail
old_doc=pre+stubs+L(abs_i,abs_i+2)+"\n\\maketitle\n"+NOTE+old_main+"\n\\clearpage\n"+old_em+"\n\\clearpage\n"+old_moved+tail
for name,t in (('new',new_doc),('old',old_doc)):
    for f in ('fig1_typicality','fig1_equilibrium','fig2_dynamics_4panel'):
        t=t.replace('{%s.pdf}'%f, '{../%s.pdf}'%f)
    open('prl/diff/%s.tex'%name,'w').write(t)
print('ok')

# postprocess (after: latexdiff --math-markup=whole old.tex new.tex > diff_vs_prxq.tex):
#   python3 prl/diff/make_diff.py postprocess
import sys
if len(sys.argv)>1 and sys.argv[1]=='postprocess':
    d=open('prl/diff/diff_vs_prxq.tex').read()
    fix=("\\usepackage{physics}\n% for the diff only: plain operators (the physics versions clash with ulem in removed text)\n"
         "\\renewcommand{\\cos}{\\operatorname{cos}}\\renewcommand{\\sin}{\\operatorname{sin}}\\renewcommand{\\tan}{\\operatorname{tan}}\n"
         "\\renewcommand{\\tanh}{\\operatorname{tanh}}\\renewcommand{\\exp}{\\operatorname{exp}}\\renewcommand{\\ln}{\\operatorname{ln}}\n"
         "\\renewcommand{\\log}{\\operatorname{log}}\\renewcommand{\\det}{\\operatorname{det}}")
    d=d.replace("\\usepackage{physics}", fix, 1)
    i=d.index("\\begin{document}"); d=d[:i]+stubs+d[i:]
    open('prl/diff/diff_vs_prxq.tex','w').write(d)
