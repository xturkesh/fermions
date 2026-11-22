# First-hand results: Lean 4 prototype built in this session (verified)

Source: compiled in this session, not literature. Location in the repository: `lean/` (README in
`lean/README.md`), committed on branch `claude/sleepy-pasteur-x7ygma`.

## Setup facts (verified)
- Lean `v4.37.0-rc1`, Mathlib commit `a495318fed6482b3c90df686de5131b46a70fd5b` (master, 21 Nov 2026).
- In the cloud container, `release.lean-lang.org` and `cache.mathlib.org` were blocked by the proxy
  (HTTP 403); the toolchain was downloaded from GitHub releases and the needed Mathlib modules
  (~5,600 build jobs) were compiled from source in about one hour on 4 cores.
- `lake build` of the repository layout succeeds. Zero errors, zero warnings, zero `sorry`.
  `#print axioms` on all 12 main theorems: only `propext`, `Classical.choice`, `Quot.sound`.
- Size: 4 files, ~730 lines. Written and debugged in roughly one hour of work once Mathlib was built,
  which supports the effort estimates in lean_formalization_landscape.md for targets (1)–(3)
  (they were done in hours, not weeks, with AI assistance).

## Theorems proven (names as in the Lean source, namespace `MonitoredFermions`)

Setting: Slater determinant on finite site set; correlation matrix C (C_jk = <c_j† c_k>) is a
Hermitian projector. Measurement of n_i: outcome 1 with probability C_ii gives
`measOcc C i = C - (C_ii)^{-1} C P_i C + P_i` (P_i = single i i 1); outcome 0 gives
`measEmp C i = 1 - measOcc (1 - C) i`.

File `Measurement.lean`:
- `idem_update`: abstract ring identity: C, P idempotent, PCP = aP, a ≠ 0 ⇒ C - a⁻¹CPC + P idempotent.
- `measOcc_mul_self`, `measEmp_mul_self`: post-measurement C is again a projector (Gaussianity and
  purity preserved).
- `measOcc_isHermitian`, `measEmp_isHermitian`.
- `measOcc_apply_self/row/col`, `measEmp_apply_self/row/col`: measured site has definite occupation
  and is decoupled (row and column i vanish off the diagonal).
- `trace_measOcc`, `trace_measEmp`: particle number conserved.
- `born_prob_nonneg`, `born_prob_le_one`: 0 ≤ C_ii ≤ 1 for Hermitian projectors.
- `born_average`: C_ii·measOcc + (1-C_ii)·measEmp = (1-P_i)C(1-P_i) + P_i C P_i, i.e. the Born
  average of the post-measurement correlation matrix is the dephased one: averaged two-point
  functions follow the linear dephasing (Lindblad) dynamics, so a measurement-induced transition
  can only appear in nonlinear quantities (entanglement, cumulants).

File `CutRank.lean` (sites split as A ⊕ B):
- `block_defect`: C_AA - C_AA² = C_AB C_BA for idempotent C.
- `rank_block_defect`: rank(C_AA - C_AA²) = rank C_AB for Hermitian projectors.
- `rank_toBlocks₁₂_measure_le`: a projective occupation measurement (any site, either outcome)
  never increases rank C_AB.
- `rank_toBlocks₁₂_blockDiag_le`: gates acting inside A and inside B never increase rank C_AB.

File `EntropyBound.lean`:
- `four_log_two_mul_le_binEntropy`: 4 ln2 · p(1-p) ≤ binEntropy p on [0,1] (not previously in
  Mathlib; proof via concavity on [0,s*] and monotonicity of the derivative on [s*,1/2],
  s*(1-s*) = 1/(8 ln 2)).
- `four_log_two_mul_var_le_entropy`: S_A ≥ 4 ln2 · Var(N_A), stated on the spectrum ν_k ∈ [0,1].
- `entropy_le_log_two_mul_card`: S_A ≤ ln2 · #{k : ν_k ∉ {0,1}}.

File `MinCut.lean`:
- `rank_sub_mul_self_eq_card`: for Hermitian M, rank(M - M²) = #eigenvalues not in {0,1}.
- `entropy_le_log_two_mul_rank`: S_A ≤ ln2 · rank C_AB for every Slater determinant (S_A defined as
  Σ_k binEntropy of the eigenvalues of C_AA).

This formally confirms the "unverified derivation" flagged in lean_formalization_landscape.md and
the inequalities 4 ln2·Var(N_A) ≤ S_A ≤ ln2·rank C_AB derived (but not found in a verified source)
in rigorous_math_physics.md.

## Not yet formalized (next steps, ordered by difficulty)
1. ν_k ∈ [0,1] for eigenvalues of C_AA (to state S_A ≥ 4 ln2 Tr C_AA(1-C_AA) directly on matrices).
2. Two-site gate across the cut raises rank C_AB by at most 2; decoupled (just-measured) sites can be
   moved across the cut without changing the rank ⇒ spacetime min-cut bound for hybrid Gaussian
   circuits.
3. Born-rule trajectory measure via `ProbabilityTheory.Kernel.traj` + Peierls counting ⇒ uniform-in-L
   area law above an explicit projective-measurement threshold (a rigorous upper bound on any
   transition point).
4. Many-body bridge (Fock space, von Neumann entropy of reduced state equals Σ binEntropy(ν_k)).
5. Existence/absence of the 1D transition: out of reach (no mathematical proof exists to formalize).

## Caution for the report writer
lean_formalization_landscape.md mentions an unreferenced, unrefereed orphan commit claiming
θ(p_c) = 0 in all dimensions d ≥ 2 (a famous open problem for 3 ≤ d ≤ 10). It is unverified; do
not rely on it. If mentioned at all, flag it explicitly as an unverified, unreviewed claim.
