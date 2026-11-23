#!/usr/bin/env python3
"""Generate statements.tex: every statement of the paper (verbatim, paper numbering) with its Lean names.

Run from paper/summary/:  python3 make_statements.py
Inputs: ../sections/*.tex, ../main.aux (citation numbers), ../../lean/README.md (Lean correspondence).
"""
import re, os

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.dirname(HERE)
SECTIONS = ['model', 'onebody', 'main', 'dimer', 'twod', 'u1', 'numerics', 'discussion']
ENVS = ['theorem', 'proposition', 'lemma', 'corollary', 'conjecture', 'hypothesis', 'definition']

# ---------------------------------------------------------------- citations (paper numbering)
aux = open(os.path.join(PAPER, 'main.aux'), encoding='utf-8').read()
bibnum = dict(re.findall(r'\\bibcite\{([^}]*)\}\{([^}]*)\}', aux))

def cite_repl(m):
    opt, keys = m.group(1), m.group(2)
    nums = ', '.join(bibnum.get(k.strip(), '?') for k in keys.split(','))
    return '[' + nums + (', ' + opt if opt else '') + ']'

def fix_cites(s):
    return re.sub(r'\\cite(?:\[([^\]]*)\])?\{([^}]*)\}', cite_repl, s)

# ---------------------------------------------------------------- section-local macros
macros, seen = [], set()
for sec in SECTIONS:
    for line in open(os.path.join(PAPER, 'sections', sec + '.tex'), encoding='utf-8'):
        m = re.match(r'\s*\\providecommand\{(\\[A-Za-z]+)\}', line)
        if m and m.group(1) not in seen:
            seen.add(m.group(1))
            macros.append(line.split('%')[0].rstrip())

# ---------------------------------------------------------------- Lean correspondence from the README
readme = open(os.path.join(PAPER, '..', 'lean', 'README.md'), encoding='utf-8').read()
lean = {}   # label -> list of (parts, names, files)
for row in readme.split('\n'):
    if not row.startswith('| ') or row.startswith('| ---') or row.startswith('| Paper'):
        continue
    cols = re.split(r'(?<!\\)\|', row)
    if len(cols) < 5:
        continue
    c1, c2, c3 = cols[1], cols[2], cols[3]
    names = [n for n in re.findall(r'`([^`]+)`', c2) if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_'.]*", n)]
    files = re.findall(r'`([^`]+\.lean)`', c3)
    for m in re.finditer(r'`((?:thm|prop|lem|cor|def|hyp|conj):[A-Za-z0-9-]+)`((?:\([a-z0-9]+\)(?:\s*[,–-]\s*\([a-z0-9]+\))*)*)', c1):
        lean.setdefault(m.group(1), []).append((m.group(2), names, files))

NOT_FORMALIZED = {
    'prop:LDT': 'literature result (Furstenberg; uniform large deviations), not formalized: '
                'it enters the formalization as the hypothesis \\texttt{LDT}.',
    'prop:Zline': 'literature result (Ge--Zhao), not formalized; for Bernoulli potentials it also follows, '
                  'given \\texttt{LDT}, from the formalized Corollary~\\ref{cor:ECinfvol}.',
    'lem:hypcheck': 'verification of the hypotheses of the cited theorems in their own formalism, not formalized.',
    'conj:logent': 'conjecture, not proved (not formalized).',
    'conj:volmix': 'conjecture, not proved (not formalized).',
}

def tt(name):
    return '\\texttt{' + name.replace('_', '\\_\\allowbreak{}').replace('.', '.\\allowbreak{}').replace("'", "{\\textquotesingle}") + '}'

def lean_block(label):
    if label in NOT_FORMALIZED:
        return '\\leanline{' + NOT_FORMALIZED[label] + '}\n'
    rows = lean.get(label)
    if not rows:
        return '\\leanline{\\emph{no entry}}\n'
    out = []
    for parts, names, files in rows:
        head = ''  # the README rows say which parts they cover; we list names and files only
        out.append(head + ', '.join(tt(n) for n in names) + ' \\,(' + ', '.join(tt(f) for f in files) + ')')
    return '\\leanline{' + ';\\\\ '.join(out) + '}\n'

# ---------------------------------------------------------------- statements
def read_optarg(s, i):
    """s[i] == '[': return (content, index after the matching ']') with brace awareness."""
    depth_b, j = 0, i + 1
    while j < len(s):
        c = s[j]
        if c == '{':
            depth_b += 1
        elif c == '}':
            depth_b -= 1
        elif c == ']' and depth_b == 0:
            return s[i + 1:j], j + 1
        j += 1
    raise ValueError('unbalanced optional argument')

out = []
missing, inner_labels = [], []
for sec in SECTIONS:
    src = open(os.path.join(PAPER, 'sections', sec + '.tex'), encoding='utf-8').read()
    src = re.sub(r'(?<!\\)%.*', '', src)          # drop comments
    msec = re.search(r'\\section\{(.*?)\}\\label\{(sec:[^}]*)\}', src)
    items = []
    for m in re.finditer(r'\\begin\{(' + '|'.join(ENVS) + r')\}', src):
        env, i = m.group(1), m.end()
        title = ''
        if src[i] == '[':
            title, i = read_optarg(src, i)
        end = src.index('\\end{' + env + '}', i)
        body = src[i:end]
        lab = re.match(r'\s*\\label\{([^}]*)\}', body)
        label = lab.group(1)
        body = body[lab.end():]
        # labels inside the statement: equations keep the paper's number, others are dropped
        def lab_repl(mm):
            inner_labels.append((label, mm.group(1)))
            return '\\tag{\\ref{' + mm.group(1) + '}}' if mm.group(1).startswith('eq:') else ''
        body = re.sub(r'\\label\{([^}]*)\}', lab_repl, body)
        items.append((m.start(), env, title, label, body))
    if not items:
        continue
    out.append('\\subsection*{\\S\\,\\ref{' + msec.group(2) + '}\\quad ' + msec.group(1) + '}\n'
               '\\addcontentsline{toc}{subsection}{\\S\\,\\ref{' + msec.group(2) + '}\\quad ' + msec.group(1) + '}\n')
    for _, env, title, label, body in items:
        kind = 'stmtD' if env == 'definition' else 'stmtP'
        out.append('\\renewcommand{\\stmtname}{\\Cref{' + label + '}}\n')
        out.append('\\begin{' + kind + '}' + ('[' + fix_cites(title) + ']' if title else '') + '\n')
        out.append(fix_cites(body).strip() + '\n')
        out.append('\\end{' + kind + '}\n')
        out.append(lean_block(label))
        if label not in lean and label not in NOT_FORMALIZED:
            missing.append(label)

open(os.path.join(HERE, 'statements.tex'), 'w', encoding='utf-8').write(
    '% Generated by make_statements.py -- do not edit by hand.\n' + ''.join(out))
open(os.path.join(HERE, 'localmacros.tex'), 'w', encoding='utf-8').write(
    '% Section-local macros of the paper, collected by make_statements.py.\n' + '\n'.join(macros) + '\n')
n = sum(1 for x in out if x.startswith('\\renewcommand{\\stmtname}'))
print('statements:', n)
print('labels without Lean entry:', missing)
print('inner labels:', [l for _, l in inner_labels if not l.startswith('eq:')])
