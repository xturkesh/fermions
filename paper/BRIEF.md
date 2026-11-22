# Brief for the paper "Rigorous disorder-free localization in a Z2 lattice gauge theory"

This brief is shared by all agents writing, refereeing or assembling the paper. Read all of it before starting.

## 0. Goal and standard

- **The paper.** A math-ph style paper (rigorous; theorem–proof), readable by physicists.
- **The claim.** The first rigorous proof of disorder-free localization (DFL) in a lattice gauge theory, namely the 1D Z2 LGT with free fermionic matter of Smith–Knolle–Kovrizhin–Moessner (SKKM 2017, ref. `SKKM2017`). It comes with:
  - a precise taxonomy of entropies;
  - a proof that DFL depends on the initial state, via a dimer counterexample;
  - exact 2D reductions and a strong-coupling proposition;
  - kinematic (Gauss-law) bounds on the number entropy in U(1) quantum link models.
- **Standard for proofs.** Every statement labelled Theorem / Proposition / Lemma / Corollary must have a complete proof. Alternatively, it can be explicitly conditional on a clearly stated hypothesis, or cite a precise theorem from the literature (with theorem number) whose hypotheses you checked.
- **No overclaiming.** If something cannot be proven, call it a Conjecture or Remark and say exactly what is missing.

## 1. Files

- Paper root: `/home/user/fermions/paper/`.
  - Section files: `sections/*.tex`. Each is `\input` by `main.tex`. A section file contains no preamble and no `\begin{document}`; it starts with `\section{...}\label{sec:...}`.
  - Shared macros: `macros.tex`. Use them; do not redefine them. You may add NEW macros at the top of your section file with `\providecommand`.
  - Bibliography: `refs.bib`, used with `\cite{key}`. Use the existing keys. If you need a new reference, add a bib entry to `refs.bib`, appending at the end with a `%` comment saying how it was verified. Never invent bibliographic data: if unsure, add `note = {to be checked}`.
  - Figures: `figs/`. Numerics code: `code/`.
  - Theorem environments (defined in `main.tex`): `theorem, proposition, lemma, corollary, definition, remark, hypothesis, conjecture`, with `\begin{proof}...\end{proof}`.
  - Labels are prefixed by section: `thm:`, `prop:`, `lem:`, `eq:`, `hyp:`, `sec:`.
- Earlier working notes. They contain most of the mathematics; reuse them, but **apply the referee errata**.
  - `/home/user/fermions/notes/mbl-lgt-number-entropy/theorems/lgt_theorems.md`
  - `/home/user/fermions/notes/mbl-lgt-number-entropy/theorems/referee_errata.md`: section "LGT note". Its fixes are mandatory.
  - Earlier scripts: `/home/user/fermions/notes/mbl-lgt-number-entropy/code/`. Relevant: `lgt1d_ed.py`, `ent_small.py`, `dfl_bounds.py`, `imbalance_bound.py`, `dimer_check.py`, `lgt2d_ed.py`, `qlm_gauss.py`, `ring_sectors.py`, `rw_potential.py`.
- Full texts of key papers (pdftotext output, verified this session):
  - Location: `/tmp/claude-0/-home-user-fermions/e9e4eeb3-4dfa-5769-bb28-5161ded2b166/scratchpad/papers/`.
  - Files: `1701.04748.txt` (SKKM), `1803.06574.txt` (SKKM 2D), `2110.00097.txt` (Macera–Sodin), `2307.07626.txt` (Macera 2023, contains the Ge–Zhao statement and a short proof), `math-ph_0405024.txt` (Jitomirskaya–Schulz-Baldes–Stolz), `1706.06135.txt` (Bucaj et al.), `1803.10697.txt` (Jitomirskaya–Zhu), `1705.09143.txt`.
  - You can also download more arXiv papers with `curl -L -o X.pdf https://arxiv.org/pdf/<id>` followed by `pdftotext -layout`. Put downloads in a NEW empty directory under the scratchpad, never inside the repo.
- Tools: Python 3 with numpy, scipy and matplotlib; LaTeX (`pdflatex`, `latexmk`).

## 2. Conventions (MANDATORY, identical in all sections)

**Lattice.** Sites j = 1,…,N carry spinless fermions f_j, with n_j = f_j† f_j. Links b_j = (j, j+1), j = 1,…,N, carry spins 1/2. Link b_N is a **boundary link**: no hopping crosses it. We set σ^x_{b_0} := 1.

**Hamiltonian (SKKM Eq. (1); VERIFIED).**
```
H = −J Σ_{j=1}^{N−1} σ^z_{b_j} (f_j† f_{j+1} + h.c.) − h Σ_{j=1}^{N} σ^x_{b_{j−1}} σ^x_{b_j},   J > 0, h ∈ R∖{0}.
```
SKKM write the two-link term exactly as −h Σ_i σ^x_{i−1,i} σ^x_{i,i+1}.

**Gauss operators and charges.** G_j = σ^x_{b_{j−1}} (−1)^{n_j} σ^x_{b_j}, with eigenvalues q_j ∈ {±1}. SKKM write q_j = τ^z_j (−1)^{n_j} with τ^z_j = σ^x_{j−1,j} σ^x_{j,j+1}.

**Dual one-body operator in sector q.**
```
h_q = −J Δ_N + 2h diag(q)   on ℓ²({1,…,N}),
```
where Δ_N is the adjacency matrix of the path with Dirichlet ends.
- Disorder parameter: λ := 2h/J.
- Sector Hamiltonian: U† H U = ⊕_q ( dΓ(h_q) − h Σ_j q_j ).
- This matches SKKM Eq. (2): H_{q} = −J Σ c_i† c_j + 2h Σ_j q_j (c_j† c_j − 1/2).

**Dressed fermions.** c̃_j := f_j Π_{k≥j} σ^z_{b_k}. The Wilson line runs to the boundary link.

**Initial states.**
- Ψ0 = |↑…↑⟩_z ⊗ ψ. Then the q are exactly i.i.d. Bernoulli(1/2). This is SKKM Eq. (5); our boundary-link convention reproduces their "sum over all 2^N charge configurations".
- The variants "open" (no boundary link: parity constraint, density ≤ 2 relative to the product law) and "ring" (Wilson loop W; twist = W) are remarks only.

**Region and entropies.**
- A = sites 1..ℓ plus links b_1..b_ℓ.
- S_G := E_q S(ρ^q_A): entanglement conditional on all Gauss charges.
- S_full: matter + gauge entanglement of the actual state.
- S_mix := S(E_q ρ^q_A).
- S_N: number entropy of N_A = Σ_{x≤ℓ} n_x.

**Eigenfunction correlator.** For a finite interval Λ,
```
Q^Λ(x,y) := sup{ |f(H_Λ)_{xy}| : f Borel, |f| ≤ 1 } = Σ_α |φ_α(x) φ_α(y)|   (simple spectrum).
```
Simple spectrum holds for Jacobi matrices.

**Normalization of the one-body problem.** h_q / J = −Δ + λ diag(q). This is the 1D Anderson model with Bernoulli(1/2) potential ±λ, i.e. Macera–Sodin's setting with W = 1 and L_x ≡ 1 (up to sign and scale).

## 3. Verified literature facts (quote these, with theorem numbers)

1. **SKKM 2017 (arXiv:1701.04748v2).**
   - Eqs. (1)–(5) as above; open boundary conditions; N up to 200.
   - Initial fermion states: domain wall, density wave, ground state at h = 0.
   - Observables: imbalance, domain-wall melting, light cone.
   - **No rigorous claims.**
2. **Macera–Sodin, AHP 2022 (arXiv:2110.00097v2).**
   - Setting: quasi-1D operators (Hψ)(x) = L_x ψ(x+1) + V_x ψ(x) + L_{x−1}^T ψ(x−1).
   - Assumptions, for W = 1 and L ≡ 1:
     - (A): finite fractional moment;
     - (B): holds when 1 ∈ supp L;
     - (C): supp V irreducible (trivial for W = 1), and supp V − supp V contains a rank-one matrix, i.e. a nonzero number. **Bernoulli ±λ with λ ≠ 0 satisfies (A)–(C).**
   - Q_I^Λ(x,y) := sup{‖f(H_Λ)_{x,y}‖ : supp f ⊂ I, |f| ≤ 1}, and Q_I(x,y) := sup_{a≤x,y≤b} Q_I^{[a,b]}(x,y), the supremum over ALL finite intervals with Dirichlet restriction.
   - **Theorem 2:** for any compact I, P{ limsup_{x→±∞} (1/|x|) log Q_I(x,y) ≤ −inf_{E∈I} γ_W(E) } = 1. This is an almost-sure statement.
   - **Proof of Thm 2 (Sect. 2.5):** on the event of Prop. 2.1 (probability ≥ 1 − Ce^{−cN}, with N²/10 ≤ |y| ≤ N²/5) one has Q_I^Λ(0,y) ≤ e^{−(γ−2τ)|y|} for |y| > C_0, uniformly in Λ ⊃ {0,y}.
     - CAUTION: the resolvent expansion uses boxes [x−N, x+N] ⊂ Λ; check how boundaries are handled before claiming uniformity for pairs near ∂Λ.
   - **Prop. 2.1:** P{ max_{E∈I} diam(Res(τ,E,N) ∩ [−N², N²]) > 2N } ≤ C e^{−cN}, where x is (τ,E,N)-non-resonant if ‖G_E[H_{[x−N,x+N]}](x, x±N)‖ ≤ e^{−(γ−τ)N}.
   - **Prop. 1.2:** a large-deviation estimate uniform over E ∈ I and over Lagrangian subspaces F, i.e. boundary conditions. READ IT.
   - Eq. (23): Q^Λ_I(x,y) ≤ lim_{ε→0} (ε/2) ∫_I ‖G_E[H_Λ](x,y)‖^{1−ε} dE ≤ W.
   - **Theorem 4 + Remark 4.1 (half-line):** with Dirichlet or any Lagrangian boundary condition, pure point spectrum and exponential decay of eigenfunctions. "Combining the proof of Theorem 2 with the additional argument described below, one can also prove dynamical localisation." Gorodetski–Kleptsyn treat W = 1 with Dirichlet.
3. **Macera 2023 (arXiv:2307.07626).**
   - **Theorem 1.1 (Ge–Zhao):** for H = Δ + V on ℓ²(Z), with V_0 bounded and not concentrated at one point, γ_E := −limsup_{|y−x|→∞} (1/|y−x|) log E Q(x,y) > 0. Here Q is the infinite-volume correlator.
   - **Sect. 2:** a ten-line proof. With N = ⌊y/10⌋ and the event B = {∃E : 0, y ∈ Res(τ,E,N)}, P(B) ≤ Ce^{−cN} by Claim 1 (MS Prop. 2.1). On B^c, for each E one of 0, y is non-resonant, and one resolvent step gives |G_E[H_{[−L,L]}](0,y)| ≤ e^{−τN}(…). Hence E Q_L(0,y) ≤ Ce^{−cN} + 2e^{−τN}, **uniformly in L** (for the symmetric boxes [−L, L] containing [−N, N] and [y−N, y+N]).
   - **Theorem 1.2:** for Bernoulli disorder the annealed rate γ_E stays bounded as the disorder grows (it is NOT the Lyapunov exponent).
4. **Ge–Zhao, J. Spectr. Theory 10 (2020) 887–904:** the statement is as in item 3. Cite through Macera 2023, Thm 1.1, plus the original.
5. **Jitomirskaya–Schulz-Baldes–Stolz, CMP 233 (2003) (math-ph/0405024).**
   - Setting: random polymer models (H ψ)(n) = −t(n+1)ψ(n+1) + v(n)ψ(n) − t(n)ψ(n−1).
   - **Definition 1:** E_c is critical if the polymer transfer matrices T_±^{E_c} are elliptic (|Tr| < 2) or equal to ±1, and commute.
   - **Theorem 1:** "There exists a constant C such that for EVERY configuration ω and for q ≥ 0, M_{ω,q}(T) ≥ C T^{q−1}", where M_{ω,q}(T) = ∫_0^∞ (dt/T) e^{−t/T} ⟨0|e^{iHt}|X|^q e^{−iHt}|0⟩. This is deterministic.
   - **Theorem 4:** if |⟨e^{2iη_±}⟩| < 1, then a.s. β^±_{ω,q} ≥ 1 − 1/(2q) (overdiffusive).
   - **Random dimer example:** v̂_+ = (λ, λ), v̂_− = (−λ, −λ), t ≡ 1. The critical energies are E_c = ±λ "as long as λ < 1".
   - **Remark 5:** ⟨e^{4iη}⟩ ≠ 1 iff λ ≠ 1/√2 (needed only for Thm 2).
   - Our computation for our dimer law: at E_c = λ (with t = 1), η_+ = π, and η_− = 2θ with cos θ = −λ, so θ ∈ (π/2, π). Then ⟨e^{2iη}⟩ = ½ + ½ e^{4iθ}, which has modulus < 1 for all λ ∈ (0, 1). **So JSBS Thm 4 applies for all 0 < 2h < J.** Re-verify this computation in the dimer section.
6. **Bucaj et al. (arXiv:1706.06135):** spectral and exponential dynamical localization, sup_t |⟨δ_n, e^{−itH}δ_m⟩| ≲ e^{ε|m|} e^{−β|n−m|} almost surely (their Thm 6.4); half-line versions mentioned.
7. **SKKM 2D (arXiv:1803.06574):** H = −Σ J_jk σ^z_jk f_j† f_k − Σ h_j A_j with star operators A_j = Π_k σ^x_jk (their Eqs. (1)–(2)); plaquettes B_p.

## 4. Mandatory errata (from the referee) to apply

- **Imbalance bound.** Use (J+|h|) everywhere, never (J+h); ‖h_q‖ ≤ 2J + 2|h|. Saturation in N is "slowly converging" for λ ≤ 0.5.
- **Lieb–Robinson.** The zero-velocity LR bound is only SECTOR-AVERAGED (state-dependent). Say explicitly that no operator-norm zero-velocity LR bound holds, because the clean sector q ≡ +1 is ballistic.
- **Thermodynamic limit / boundary.** The chain has a fixed left boundary and A = [1, ℓ] touches it. Infinite-line theorems (Ge–Zhao) do not directly give finite-volume uniform bounds near a boundary. Handle this explicitly:
  - (i) via Hypothesis (EC) on finite intervals, with a careful status discussion; and/or
  - (ii) via a proof (onebody section); and/or
  - (iii) via two-sided chains with bulk observables.
- **Entropy inequality.** P5(3) is false for S_G: S_G can be < S_N (numerical example: N = 10, h = J, ℓ = 5, t = 1 gives S_G = 0.61 < S_N = 0.67). Only S_full, S_mix and S_f are ≥ S_N.
- **Number-entropy bound.** Report the theorem's actual constant. The "1.09" was a non-rigorous plug-in.
- **Dimer.** State the critical-energy condition 2h < J and verify JSBS's hypotheses (see §3.5). Avoid resonant special values if a theorem needs it.
- **2D.**
  - Torus: add the non-contractible Wilson loops.
  - Strong coupling: state the smallness condition J/g ≤ c_0 for the Schrieffer–Wolff step, or give the trivial complement.
- **U(1) QLM.** The number entropy for superpositions over charge sectors adds "at most H(Σ_A Q)". Wilson-link bound: fix the hopping-norm normalization.
- **ADHH corollaries.** Either verify the ADHH statement from the arXiv text (1509.05386) or label them as remarks ("by [ADHH] one expects …").
- **Random-walk potential.** Cite the half-line Simon–Spencer statement only as a remark.

## 5. Planned structure of the paper

1. Introduction: context, DFL, MBL debate; what is new; summary of results.
2. `sec:model`: model, exact duality (Thm 2.1), sector averaging (Cor 2.2), boundary conventions, dictionary, cut-locality, Gauss-law projection remark (single sector ⇒ ballistic).
3. `sec:onebody`: the one-body input. Bernoulli localization on Z (Ge–Zhao); finite-interval uniformity (Hypothesis (EC), its proof or status); absence of critical energies; consequences (sup_t propagator bounds, time-averaged transition probabilities).
4. `sec:main`: main theorems for the LGT.
   - (a) Dressed-fermion localization and bounded transport.
   - (b) Sector-averaged zero-velocity LR bound.
   - (c) Positive infinite-time imbalance for every h ≠ 0.
   - (d) Entanglement: area law for S_G and for S_full of definite-charge states; ln 2 duality defect; entropy taxonomy and correct inequalities.
   - (e) Number entropy bounded, with explicit constant.
   - (f) Infinite-volume versions.
5. `sec:dimer`: DFL depends on the initial state. A product gauge state ⇒ random dimer sector law ⇒ for 0 < 2h < J transport at least diffusive for every configuration (JSBS Thm 1), and a.s. overdiffusive (JSBS Thm 4). Biased product states (1-block factors) as a remark.
6. `sec:2d`: 2D mapping (Bernoulli or binary π-flux), strong-coupling cluster confinement up to t ~ h/J², open problems (Ding–Smart, L. Li, π-flux).
7. `sec:u1`: number entropy in U(1) gauge theories. Gauss-law caps for quantum links; Wilson links; consequences for "number entropy" diagnostics in LGTs. The single-link Z2 model maps a single charge to a random-walk potential (remark; localization open).
8. `sec:numerics`: numerical illustrations (figures).
9. Discussion and open problems: interacting DFL = binary-disorder MBL; S_full ~ ln t conjecture; Loschmidt / spin observables; etc.
