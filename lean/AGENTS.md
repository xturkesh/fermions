# Instructions for agents working on the Lean formalization (`/home/user/fermions/lean`)

Project: Lean 4 (v4.36.0) + Mathlib (tag v4.36.0) formalization of the paper
"Rigorous disorder-free localization in a ℤ₂ lattice gauge theory" (`../paper/`, conventions in
`../paper/BRIEF.md`; the earlier mathematical notes are in
`../notes/mbl-lgt-number-entropy/theorems/lgt_theorems.md` together with the referee errata in
`referee_errata.md`).

## Layout

`DFL/Basic.lean` holds the definitions only. Every other file holds statements whose proofs are
`sorry`. Each agent OWNS some files and must replace every `sorry` in them by a complete proof.

## Hard rules

1. **Do not change any existing definition or theorem statement** (name, hypotheses, conclusion).
   They are the interface other agents rely on.
   - You may add new definitions, lemmas and instances in the files you own, placed before their
     first use.
   - If you are convinced a statement is false or unprovable as written, do NOT change it. Stop
     working on it and explain precisely in your final report (with a counterexample if possible).
2. **Edit only the files you own.** Never touch `DFL/Basic.lean`, `DFL.lean` or `lakefile.toml`.
3. **No `sorry`, no `admit`, no `axiom`, no `native_decide`** in the final versions.
   - At the end, check with a scratch file `scratch_<yourname>.lean` (in `/home/user/fermions/lean/`,
     delete it afterwards) containing `import DFL.<File>` and `#print axioms DFL.<theorem>`. Only
     `propext`, `Classical.choice`, `Quot.sound` are allowed for theorems in your files, except that
     theorems may depend on `sorryAx` *through other agents' files* while those are unfinished.
     Report which.
4. **Compilation.** Check a file with
   `cd /home/user/fermions/lean && PATH=$HOME/.elan/bin:$PATH lake env lean DFL/<File>.lean`.
   - **Never run `lake build`**: other agents edit files concurrently, and a build would recompile
     their half-edited files.
   - Imports are taken from the prebuilt `.olean` files of the skeleton, so you can use every
     statement of the imported files as it is in the skeleton (proved or not). You cannot use
     lemmas that another agent adds later.
   - If you own two files A and B with B importing A, B can only use the *skeleton statements* of
     A, so put helper lemmas needed by B in B itself.
5. **Resources.** The machine has 4 CPUs and 15 GB RAM shared by about 7 agents. Run at most ONE
   Lean process at a time, and never run Lean in the background.
   - A typical check takes 10–60 s.
   - Keep files compiling: develop large proofs incrementally, and keep a compiling version with
     `sorry` for unfinished parts rather than leaving the file broken for long.
6. **Mathlib API (v4.36.0).** The source is in `.lake/packages/mathlib/Mathlib`; use `grep -rn`.
   - `exp` is `NormedSpace.exp` with no field argument. Matrix-exponential lemmas are in
     `Mathlib/Analysis/Normed/Algebra/MatrixExponential.lean` (`Matrix.exp_blockDiagonal`,
     `Matrix.exp_conjTranspose`, `Matrix.exp_add_of_commute`, `Matrix.exp_neg`, `Matrix.exp_conj`,
     …).
   - Analytic facts (derivatives `hasDerivAt_exp_smul_const`, norms) need a normed-algebra
     instance on matrices, activated locally, e.g.
     `attribute [local instance] Matrix.linftyOpNormedRing Matrix.linftyOpNormedAlgebra` (check the
     exact names in `Mathlib/Analysis/Matrix/Normed.lean`).
   - Hermitian spectral theory: `Matrix.IsHermitian.eigenvectorBasis`, `.eigenvalues`,
     `.mulVec_eigenvectorBasis`, `.spectral_theorem`, `eigenvectorUnitary` in
     `Mathlib/Analysis/Matrix/Spectrum.lean`.
   - `EuclideanSpace` is `WithLp 2 _`, a structure: use `⇑v`, `PiLp`/`WithLp.ofLp` API.
   - `Real.negMulLog`, `Real.binEntropy` and their lemmas are in Mathlib.
7. **Style.** Short helper lemmas, `simp`/`ring`/`linarith`/`nlinarith`/`positivity`/`fin_cases`/
   `decide`/`omega` where appropriate. Write a one-line docstring for each new top-level lemma.
8. **Final report** (your last message):
   - for every theorem in your files: proved / not proved (with the exact blocker);
   - new lemmas you added;
   - anything suspicious about the statements.
