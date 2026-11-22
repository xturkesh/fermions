# Lean 4 formalization landscape for monitored free fermions (status: 22 November 2026)

> **How this was verified.** On 22 Nov 2026 I made shallow clones of the primary repositories from GitHub and grepped their source at fixed commits:
> - Mathlib `a495318f` (committed 21 Nov 2026; toolchain `v4.37.0-rc1`)
> - Physlib `56430c65` (21 Nov 2026)
> - Lean-QuantumInfo `56e83a92` (5 May 2026, the final snapshot before its merge into Physlib)
> - formal-conjectures `2207c54d` (21 Nov 2026)
> - anthropics/formal-math: `main` = `16a5ce41` (26 Oct 2026), plus the unreferenced commit `795efb86` (28 Aug 2026)
> - elan `master` (18 Nov 2026) and the Lake sources on lean4 `master` (22 Nov 2026)
> - math-inc/strongpnt (10 Sep 2025) and Hayata-Yamasaki-Group/lean-quantum (28 Oct 2026)
>
> **Caveats.**
> - "Sorry-free" means only that no `sorry` token appears in the source. I could not compile anything or run `#print axioms`, because the Lean toolchain and cache hosts are blocked in this container.
> - arXiv, leanprover-community.github.io (the Zulip archive and the docs site) and IEEE Spectrum were also unreachable. Claims that rest only on search-result snippets are marked **(snippet)**.
> - Line numbers (`#Lnn`) refer to the pinned commits. Mathlib moves fast: for example, the Hermitian spectral theory moved recently from `LinearAlgebra/Matrix` to `Analysis/Matrix`. Expect paths and names to drift.

## 1. Mathlib master (late 2026): what exists for matrix analysis, entropy functions, Clifford/exterior algebras, probability, graphs and percolation?

### Takeaway
Mathlib already covers the core finite-dimensional tools: Hermitian spectral theory, positive semidefinite (PSD) and Schur-complement criteria, Woodbury and the matrix-determinant lemma, the matrix exponential, a Hermitian continuous functional calculus (CFC), a scoped Loewner order, and operator monotonicity/concavity of `log`, `sqrt` and `x^p` (p ∈ [0,1]). It also has `negMulLog` and `binEntropy` with their derivatives and concavity, Clifford and exterior algebras with contraction and a graded splitting Cl(Q₁⊕Q₂) ≅ Cl(Q₁) ⊗̂ Cl(Q₂), and a strong probability layer: Ionescu–Tulcea `Kernel.traj`, infinite product measures, `setBernoulli`, martingale convergence, the strong law of large numbers (SLLN), the central limit theorem (CLT) and Brownian motion.

It does **not** have:
- Courant–Fischer or Cauchy interlacing (Physlib has them; see §2)
- the inequality binEntropy p ≥ 4 ln2 · p(1−p)
- a Hilbert-space Fock space or the canonical anticommutation relations (CAR)
- the pointwise Birkhoff ergodic theorem (only the maximal and mean ergodic theorems), and nothing from Kingman, Furstenberg or Oseledets
- Menger's theorem, max-flow/min-cut, self-avoiding-walk counting, percolation, Ising or Peierls arguments.

### Cited Findings

#### 1a. Finite-dimensional matrix analysis

**Spectral theory.** `Mathlib/Analysis/Matrix/Spectrum.lean` contains:
- `Matrix.IsHermitian.eigenvalues : n → ℝ` [L67](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/Matrix/Spectrum.lean#L67)
- the sorted `eigenvalues₀` with `eigenvalues₀_antitone` (L60, L63)
- `eigenvectorBasis : OrthonormalBasis n 𝕜 (EuclideanSpace 𝕜 n)` (L71) and `eigenvectorUnitary` (L89)
- `spectral_theorem` [L143](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/Matrix/Spectrum.lean#L143), `eigenvalues_eq` (L148), `charpoly_eq` (L157), `det_eq_prod_eigenvalues` (L193)
- `rank_eq_rank_diagonal` (L197) and `rank_eq_card_non_zero_eigs : A.rank = Fintype.card {i // hA.eigenvalues i ≠ 0}` [L202](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/Matrix/Spectrum.lean#L202)
- `spectrum_real_eq_range_eigenvalues` (L213), `eigenvalues_eq_zero_iff` (L219), and `trace_eq_sum_eigenvalues` [L240](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/Matrix/Spectrum.lean#L240)

**Positive (semi)definite matrices.**
- The definitions live in `Mathlib/LinearAlgebra/Matrix/PosDef.lean` (783 lines).
- Eigenvalue facts are in `Mathlib/Analysis/Matrix/PosDef.lean` — [source](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/Matrix/PosDef.lean#L33):
  - `Matrix.IsHermitian.posSemidef_iff_eigenvalues_nonneg` (L33)
  - `Matrix.PosSemidef.eigenvalues_nonneg` (L41)
  - `PosSemidef.trace_eq_zero_iff` (L54)
  - `IsHermitian.posDef_iff_eigenvalues_pos` (L73)
- Compressions preserve PSD: `Matrix.PosSemidef.submatrix` [PosDef.lean L84](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/LinearAlgebra/Matrix/PosDef.lean#L84) and `Matrix.PosDef.submatrix` (L188).

**Schur complements.**
- `Matrix.PosSemidef.fromBlocks₁₁` states `(fromBlocks A B Bᴴ D).PosSemidef ↔ (D − Bᴴ A⁻¹ B).PosSemidef` [L699](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/LinearAlgebra/Matrix/PosDef.lean#L699); `fromBlocks₂₂` is at L721.
- Determinants: `Matrix.det_fromBlocks₁₁` and `det_fromBlocks₂₂` [SchurComplement.lean L368/L382](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/LinearAlgebra/Matrix/SchurComplement.lean#L368).
- Matrix determinant lemma: `Matrix.det_one_add_replicateCol_mul_replicateRow` and `det_add_replicateCol_mul_replicateRow` [L416/L426](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/LinearAlgebra/Matrix/SchurComplement.lean#L416).

**Woodbury identity.** `Matrix.add_mul_mul_inv_eq_sub` (Woodbury) and `Matrix.add_mul_mul_inv_eq_sub'` (binomial inverse theorem) — [NonsingularInverse.lean L616/L627](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/LinearAlgebra/Matrix/NonsingularInverse.lean#L616).

**Rank.**
- `Matrix.rank_mul_le` [Rank.lean L207](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/LinearAlgebra/Matrix/Rank.lean#L207), `Matrix.rank_conjTranspose_mul_self` (L527) and `Matrix.rank_self_mul_conjTranspose` [L549](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/LinearAlgebra/Matrix/Rank.lean#L549).
- I found no theorem named `rank_add_le` (rank subadditivity) in `Mathlib/LinearAlgebra/Matrix/` (my grep).

**Projectors and matrix units.**
- `IsStarProjection` (p·p = p ∧ star p = p) with `IsStarProjection.one_sub` (L73) and `pow_eq` (L64) — [StarProjection.lean](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Algebra/Star/StarProjection.lean#L73).
- The matrix unit is `Matrix.single i j a` [Data/Matrix/Basis.lean L34](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Data/Matrix/Basis.lean#L34).

**Matrix exponential.** `Mathlib/Analysis/Normed/Algebra/MatrixExponential.lean` has `Matrix.exp_diagonal` (L82), `exp_conjTranspose` (L93), `IsHermitian.exp` (L97), `exp_neg` (L175) and `exp_conj` (L185) — [source](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/Normed/Algebra/MatrixExponential.lean#L82).

**CFC for Hermitian matrices.**
- `Matrix.IsHermitian.instContinuousFunctionalCalculus` (L98) makes the general `cfc` available.
- Also `Matrix.IsHermitian.cfc (f : ℝ → ℝ)` (L133), `cfc_eq : cfc f A = hA.cfc f` (L144) and `charpoly_cfc_eq` (L153) — [HermitianFunctionalCalculus.lean](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/Matrix/HermitianFunctionalCalculus.lean#L133).
- That file contains no trace lemma such as `trace (cfc f A) = ∑ f(λᵢ)` (grep for "trace" is empty). Physlib has one (see §2).

**Loewner order.** `Mathlib/Analysis/Matrix/Order.lean` (Monica Omar, 2025) — [source](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/Matrix/Order.lean):
- defines `Matrix.instPreOrder` / `Matrix.instPartialOrder` with `A ≤ B := (B − A).PosSemidef`
- the order is scoped: you must `open scoped MatrixOrder`
- provides `Matrix.le_iff` and `Matrix.nonneg_iff_posSemidef`
- its stated purpose is to reuse C⋆-algebra results such as `CFC.sqrt` on matrices.

**Operator monotone and concave functions.**
- `ExpLog/Order.lean` — [source](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/SpecialFunctions/ContinuousFunctionalCalculus/ExpLog/Order.lean):
  - `CFC.log_monotoneOn` (log is operator monotone) and `CFC.concaveOn_log` (operator concave)
  - listed TODO: "Show that `x => x * log x` is operator convex".
- `Rpow/Order.lean` (Frédéric Dupuis, 2025) — [source](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/SpecialFunctions/ContinuousFunctionalCalculus/Rpow/Order.lean):
  - `CFC.monotone_rpow` and `CFC.monotone_nnrpow` (p ∈ [0,1]), `CFC.monotone_sqrt`
  - `CFC.concaveOn_rpow`, `CFC.concaveOn_nnrpow`, `CFC.concaveOn_sqrt`
  - TODO: rpow is operator antitone/convex on [−1,0] and operator convex on [1,2]
  - built on `Rpow/IntegralRepresentation.lean`.

**Interlacing / min–max.**
- `grep -i "Courant\|interlac"` over `Mathlib/Analysis` and `Mathlib/LinearAlgebra` finds nothing relevant.
- Only the sorted eigenvalues of symmetric operators exist, e.g. `LinearMap.IsSymmetric.eigenvalues_antitone` [InnerProductSpace/Spectrum.lean L308](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/InnerProductSpace/Spectrum.lean#L308).

**C⋆-algebra files.** `Mathlib/Analysis/CStarAlgebra/` contains `Matrix.lean`, `Projection.lean`, `CompletelyPositiveMap.lean`, `GelfandNaimarkSegal.lean`, `Fuglede.lean` and others — [directory](https://github.com/leanprover-community/mathlib4/tree/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/CStarAlgebra) (I checked file names only).

#### 1b. Entropy-related real analysis

**`negMulLog`.** `Real.negMulLog x = -x * log x` [Log/NegMulLog.lean L163](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/SpecialFunctions/Log/NegMulLog.lean#L163), with `Real.strictConcaveOn_negMulLog` and `Real.concaveOn_negMulLog` on `Set.Ici 0` (L223, L226).

**`binEntropy`.** `Mathlib/Analysis/SpecialFunctions/BinaryEntropy.lean` is in namespace `Real` and uses natural logs — [source](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/SpecialFunctions/BinaryEntropy.lean#L449):
- `binEntropy_eq_negMulLog_add_negMulLog_one_sub` (L72)
- `binEntropy_two_inv_add` (symmetry about ½, L84)
- `binEntropy_pos` and `binEntropy_nonneg` (L87, L94)
- `binEntropy_eq_zero : binEntropy p = 0 ↔ p = 0 ∨ p = 1` (L130)
- `binEntropy_le_log_two` (L158) and `binEntropy_eq_log_two` (L163)
- `deriv_binEntropy` (L191), `hasDerivAt_binEntropy` (L258), and `deriv2_binEntropy : deriv^[2] binEntropy p = -1 / (p * (1 - p))` (L366)
- `binEntropy_strictMonoOn` on Icc 0 ½ (L428) and `binEntropy_strictAntiOn` (L433)
- `strictConcave_binEntropy : StrictConcaveOn ℝ (Icc 0 1) binEntropy` (L449; note there is no "On" in the name)
- q-ary analogues (`qaryEntropy`).

**What is not in the file.**
- No lemma of the form `4 * log 2 * p * (1 - p) ≤ binEntropy p` (grep for `4 *` returns 0 matches).
- No lemma named `continuous_binEntropy`.

**Tools available for proving the missing bound.**
- Decimal bounds on log 2: `Real.log_two_gt_d9 : 0.6931471803 < log 2` and `Real.log_two_lt_d9` — [Analysis/Complex/ExponentialBounds.lean L86/L89](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/Complex/ExponentialBounds.lean#L86).
- Convexity from second derivatives: `convexOn_of_deriv2_nonneg` and `concaveOn_of_deriv2_nonpos` — [Analysis/Convex/Deriv.lean L215/L224](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/Convex/Deriv.lean#L215).

#### 1c. Clifford and exterior algebras as a CAR algebra / Fock space

**Contraction.** `Mathlib/LinearAlgebra/CliffordAlgebra/Contraction.lean` — [source](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/LinearAlgebra/CliffordAlgebra/Contraction.lean#L81):
- `CliffordAlgebra.contractLeft : Module.Dual R M →ₗ[R] CliffordAlgebra Q →ₗ[R] CliffordAlgebra Q` (L81)
- `contractLeft_ι_mul` (L122), the graded-Leibniz/anticommutation rule
- `contractLeft_contractLeft : d⌋(d⌋x) = 0` (L185)
- `changeForm` (L238) and `CliffordAlgebra.equivExterior [Invertible (2 : R)]` (L350).

**Graded splitting of a direct sum.** `CliffordAlgebra.equivProd : CliffordAlgebra (Q₁.prod Q₂) ≃ₐ[R] (evenOdd Q₁ ᵍ⊗[R] evenOdd Q₂)`, a graded tensor product. Its TODO is to upgrade it to a graded-algebra equivalence — [CliffordAlgebra/Prod.lean](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/LinearAlgebra/CliffordAlgebra/Prod.lean#L20).

**Exterior algebra.**
- `ExteriorAlgebra.ι : M →ₗ[R] ExteriorAlgebra R M` [ExteriorAlgebra/Basic.lean L69](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/LinearAlgebra/ExteriorAlgebra/Basic.lean#L69).
- An occupation-number basis: `Module.Basis.ExteriorAlgebra : Basis (Finset I) R (ExteriorAlgebra R M)` [ExteriorAlgebra/Basis.lean L32](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/LinearAlgebra/ExteriorAlgebra/Basis.lean#L32).
- `Mathlib/LinearAlgebra/ExteriorPower/` contains `BilinForm.lean`, `Pairing.lean`, `WedgePairing.lean`, `HodgeStar.lean` and others — [directory](https://github.com/leanprover-community/mathlib4/tree/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/LinearAlgebra/ExteriorPower). I did not read their contents.

**No Fock-space or CAR material.** "Fock" and "canonical anticommutation" do not occur in the CliffordAlgebra or ExteriorAlgebra directories (grep).

#### 1d. Probability

**Ionescu–Tulcea.** `Mathlib/Probability/Kernel/IonescuTulcea/{Maps,PartialTraj,Traj}.lean` (Etienne Marion, 2024) — [Traj.lean L519](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Probability/Kernel/IonescuTulcea/Traj.lean#L519):
- `ProbabilityTheory.Kernel.traj (a : ℕ) : Kernel (Π i : Iic a, X i) (Π n, X n)`, built from a family `κ n : Kernel (Π i : Iic n, X i) (X (n+1))`
- `trajMeasure μ₀ κ`
- main statements: `map_traj_succ_self`, `eq_traj`, `traj_comp_partialTraj`, `condExp_traj`, `condDistrib_trajMeasure`.

**Building kernels.** `ProbabilityTheory.Kernel.ofFunOfCountable` [Kernel/Basic.lean L236](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Probability/Kernel/Basic.lean#L236) builds a kernel from any function on a countable domain.

**Product measures.** `Measure.infinitePiNat` (L112) and `Measure.infinitePi` (L354, arbitrary index type) — [ProductMeasure.lean](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Probability/ProductMeasure.lean#L354).

**Bernoulli subsets.**
- `ProbabilityTheory.setBernoulli : Measure (Set ι)`, the product Bernoulli measure on subsets [Distributions/SetBernoulli.lean L44](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Probability/Distributions/SetBernoulli.lean#L44).
- It is used by `Mathlib/Probability/Combinatorics/BinomialRandomGraph/Defs.lean` [source](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Probability/Combinatorics/BinomialRandomGraph/Defs.lean).

**Limit theorems.**
- SLLN: `ProbabilityTheory.strong_law_ae` [StrongLaw.lean L787](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Probability/StrongLaw.lean#L787), `strong_law_ae_real` (L598) and `strong_law_Lp` (L831).
- CLT: `tendstoInDistribution_inv_sqrt_mul_sum` [CentralLimitTheorem.lean L77](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Probability/CentralLimitTheorem.lean#L77) (authors: Thomas Zhu, Etienne Marion).
- Brownian motion: `Mathlib/Probability/BrownianMotion/{Basic,GaussianProjectiveFamily}.lean` (Marion, 2025) — [directory](https://github.com/leanprover-community/mathlib4/tree/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Probability/BrownianMotion).

**Martingales.** `Mathlib/Probability/Martingale/` contains Basic, BorelCantelli, Centering, Convergence, OptionalSampling, OptionalStopping and Upcrossing. Key theorems: `Submartingale.exists_ae_tendsto_of_bdd` (L189), `Submartingale.ae_tendsto_limitProcess` (L207) and `…_of_uniformIntegrable` (L313) — [Convergence.lean](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Probability/Martingale/Convergence.lean#L207).

**Other directories (names only).** `Process/` (Filtration, Stopping, HittingTime, Kolmogorov, …), `Kernel/` (Composition, Disintegration, CondDistrib, Invariance, Irreducible, …), `IdentDistrib.lean` and `Independence/` — [directory](https://github.com/leanprover-community/mathlib4/tree/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Probability).

**Ergodic theory.**
- Maximal ergodic theorem: `setIntegral_birkhoffSumSup_nonneg` plus a weak-type bound for `birkhoffAverageSup` [Dynamics/BirkhoffSum/Maximal.lean](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Dynamics/BirkhoffSum/Maximal.lean) (Lua Viana Reis, 2025).
- Mean ergodic theorem (von Neumann): `ContinuousLinearMap.tendsto_birkhoffAverage_orthogonalProjection` [InnerProductSpace/MeanErgodic.lean L89](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/InnerProductSpace/MeanErgodic.lean#L89).
- I found no a.e. (pointwise Birkhoff) convergence theorem in `Mathlib/Dynamics`.
- Word greps for `Kingman` and `Furstenberg` are empty. `Lyapunov` occurs only in `Analysis/MeanInequalities.lean`, which is presumably Lyapunov's interpolation inequality (I did not inspect it); there is nothing on Lyapunov exponents.

#### 1e. Combinatorics, graphs, percolation

**Absent.** Word-boundary greps over all of `Mathlib/` for Menger, Ising, Peierls, percolation, self-avoiding, max-flow/maxFlow/minCut, Kingman and Furstenberg return no hits. (A substring grep for "ising" only matched British spellings such as "Randomisation".)

**What exists.**
- `SimpleGraph/` has Walk/, Paths, Connectivity/, Hall, Tutte, Matching, LapMatrix, Regularity/, Hamiltonian and Coloring/.
- There are also `Combinatorics/Graph/` (multigraphs), `Digraph/` and `Optimization/ValuedCSP.lean` — [directory](https://github.com/leanprover-community/mathlib4/tree/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Combinatorics/SimpleGraph).

**Size.** Mathlib has 8,730 `.lean` files under `Mathlib/` at `a495318f` (`git ls-tree` count); see the [repository](https://github.com/leanprover-community/mathlib4/tree/a495318fed6482b3c90df686de5131b46a70fd5b).

### Inferences
- **Correlation-matrix (single-particle) level is ready.** The matrix-level parts of targets (1)–(3) need no new theory. They can be built from `IsStarProjection`, `Matrix.single`, `submatrix`, `PosSemidef.submatrix`, the spectral theorem, `rank_eq_card_non_zero_eigs` and `rank_self_mul_conjTranspose`. The main friction points are:
  - eigenvalues of polynomials of a Hermitian matrix (use `spectral_theorem`, `rank_eq_rank_diagonal` or `charpoly_cfc_eq`)
  - working with compressions and block decompositions
  - the missing `rank_add_le`-style lemma.
- **Fock space would have to be built.** The exterior algebra gives the algebraic structure:
  - basis indexed by `Finset` (occupation numbers)
  - creation operators as left multiplication by `ι v`
  - annihilation operators as `contractLeft`, with the CAR-type rule `contractLeft_ι_mul` and nilpotency.

  Missing: an inner product on ⋀V, adjointness between creation and annihilation, a C⋆/Hilbert structure, Jordan–Wigner, and Gaussian states. `equivProd` handles the A|B mode splitting at the level of graded algebras, but not partial traces of states. A pragmatic shortcut is to represent everything as explicit 2^L×2^L matrices on `Finset (Fin L)` with Jordan–Wigner signs.
- **Trajectory measures fit existing APIs.** `Kernel.traj` and `ofFunOfCountable` fit Born-rule trajectories with finite outcome spaces directly. `setBernoulli` and `Measure.infinitePi` fit i.i.d. measurement locations.
- **No percolation reasoning is available.** No percolation, path-counting or Peierls/contour infrastructure is in Mathlib; any such argument must be built locally or ported (see §2).

### Gaps
- Not checked:
  - whether `ExteriorPower/BilinForm.lean` or `Pairing.lean` provide a Gram-determinant pairing ⟨v₁∧…∧vₖ, w₁∧…∧wₖ⟩ = det⟨vᵢ,wⱼ⟩ (a prerequisite for a Fock inner product)
  - `Mathlib/InformationTheory/` (Kullback–Leibler (KL) divergence, etc.)
  - the Rayleigh-quotient API.
- Open Mathlib pull requests (for example a pointwise ergodic theorem or Courant–Fischer) could not be listed. GitHub organization and search API paths were blocked for this session, and the docs site and Zulip were unreachable.
- Exact names for rank invariance under unitary multiplication and for the rank of a diagonal matrix were not verified.

## 2. External Lean projects: Lean-QuantumInfo / PhysLean (now Physlib), formal_conjectures, free fermions, Lieb–Robinson bounds, Peierls arguments, percolation, random matrices

### Takeaway
- **Physlib is now the main home for physics and quantum information in Lean.** Physlib (`leanprover-community/physlib`) was formed in March 2026 by merging PhysLean (formerly HepLean) with Lean-QuantumInfo. Its `QuantumInfo/` part (≈40k lines, sorry-free by token grep) contains von Neumann entropy, strong subadditivity (SSA), subadditivity, Araki–Lieb, weak monotonicity, joint convexity of relative entropy, a `DPI.lean` file for data-processing inequalities (sandwiched-Rényi machinery; I did not check the exact DPI theorem names), and the generalized quantum Stein's lemma (GQSL). A separate Courant–Fischer/Weyl/Cauchy-interlacing file sits in PhyslibAlpha.
- **Physlib's physics side** has Wick's theorem (operator-algebraic, for QFT), a single-particle tight-binding chain, a classical Ising-on-a-graph energy file, the Born rule for effect-valued measures, and Glimm–Jaffe Euclidean field axioms.
- **Nothing about free fermions, fermionic Gaussian states, Slater determinants, Jordan–Wigner, Lieb–Robinson bounds or quantum trajectories was found anywhere.**
- **Percolation in Lean** exists only as:
  - (a) an AI-written, unrefereed and now unreferenced artifact claiming θ(p_c)=0 in all d ≥ 2, which contains a 2D Peierls dual-circuit lemma
  - (b) a hierarchical-lattice development known only from a snippet.
- **formal_conjectures** has three quantum problems and nothing on percolation, Ising models or random matrices.

### Cited Findings

#### Lean-QuantumInfo → Physlib

**History and scope (Lean-QuantumInfo README).**
- "In March 2026, it merged with PhysLean (previously known as HEPLean) and formed Physlib … under the leanprover-community umbrella … new development will take place there!" — [Lean-QuantumInfo README](https://github.com/Timeroot/Lean-QuantumInfo/blob/56e83a9288a3c616285038748e273b3c0e1a36bf/README.md)
- The GQSL milestones, from the same README:
  - October 2025: formalizing the arguments of the Hayashi–Yamasaki paper while assuming standard results
  - April 2026: completion, once those standard results were also formalized
  - the final theorem is `limit_hypotesting_eq_limit_rel_entropy`.
- Stats "as of May 5th 2026": ~2143 theorems, 423 definitions, 38,105 lines; toolchain `v4.28.0`.
- Hayata Yamasaki contributed operator-inequality code "ported to this repository from his own" — [README](https://github.com/Timeroot/Lean-QuantumInfo/blob/56e83a9288a3c616285038748e273b3c0e1a36bf/README.md).

**Final Lean-QuantumInfo snapshot (5 May 2026).**
- 69 `sorry` tokens in 20 files, for example `AxiomatizedEntropy/Defs.lean` (12), `ResourceTheory/ResourceTheory.lean` (8) and `Capacity.lean` (7).
- In `Entropy/Relative.lean`, `qRelativeEnt_joint_convexity` and `qMutualInfo_as_qRelativeEnt` were proved by `sorry` — [Relative.lean L2091–2108](https://github.com/Timeroot/Lean-QuantumInfo/blob/56e83a9288a3c616285038748e273b3c0e1a36bf/QuantumInfo/Finite/Entropy/Relative.lean#L2091).
- `Entropy/SSA.lean` (© 2026 Meiburg) had none.

**The GQSL paper (snippet).** arXiv 2510.08672 (Meiburg, Lessa, Soldati) formalizes the GQSL following Hayashi–Yamasaki (2024). At that time the library had ">1000 theorems, 250 definitions, and 15,000 lines" — [arXiv 2510.08672](https://arxiv.org/abs/2510.08672).

**Physlib `56430c65` (21 Nov 2026).**
- Lean `v4.36.0`, depending on Mathlib tag `v4.36.0` — [lakefile.toml L12–14](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/lakefile.toml#L12).
- Sizes (my counts) — [repository](https://github.com/leanprover-community/physlib/tree/56430c65b76f73c7c25fc03eed46be1d888520bb):

  | Directory | `.lean` files | Lines | `sorry` tokens |
  |---|---|---|---|
  | `Physlib/` | 696 | 168,070 | 0 |
  | `PhyslibAlpha/` | 878 | 245,502 | 0 |
  | `QuantumInfo/` | 81 | 40,357 | 6 |

- All 6 `sorry` tokens in `QuantumInfo/` are inside comments or docstrings: commented-out lines in `ForMathlib/HermitianMat/CFC.lean` L337–354, a `/- … -/` block in [ForMathlib/Misc.lean L22–36](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/QuantumInfo/ForMathlib/Misc.lean#L22), and a docstring in `ClassicalInfo/Capacity.lean` L47.

**Entropy results in Physlib.** `QuantumInfo/Entropy/{VonNeumann,Relative,DPI,SSA}.lean` contain no `sorry` tokens. In particular:
- `Sᵥₙ_strong_subadditivity` [SSA.lean L1191](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/QuantumInfo/Entropy/SSA.lean#L1191) and `Sᵥₙ_subadditivity` (L1203)
- `Sᵥₙ_triangle_subaddivity` (Araki–Lieb; note the misspelling in the name, L1232), `Sᵥₙ_weak_monotonicity` (L1243), `qConditionalEnt_strong_subadditivity` (L1257) and `qMutualInfo_strong_subadditivity` (L1266)
- `qRelativeEnt_joint_convexity`, now proved [DPI.lean L1523](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/QuantumInfo/Entropy/DPI.lean#L1523)
- the GQSL main theorem `limit_hypotesting_eq_limit_rel_entropy` [ResourceTheory/SteinsLemma.lean L2114](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/QuantumInfo/ResourceTheory/SteinsLemma.lean#L2114).

**Von Neumann entropy API** — [VonNeumann.lean](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/QuantumInfo/Entropy/VonNeumann.lean#L72):
- `Sᵥₙ (ρ : MState d) : ℝ` (L72), `Sᵥₙ_eq_trace_cfc_negMulLog` (L123), and `Sᵥₙ_of_partial_eq` (a pure bipartite state has equal marginal entropies, L316)
- in the Lean-QuantumInfo snapshot also `Sᵥₙ_nonneg`, `Sᵥₙ_le_log_d`, `Sᵥₙ_of_pure_zero`, `Sᵥₙ_continuous`, `qConditionalEnt`, `qMutualInfo`, `coherentInfo` and `qcmi` — [snapshot](https://github.com/Timeroot/Lean-QuantumInfo/blob/56e83a9288a3c616285038748e273b3c0e1a36bf/QuantumInfo/Finite/Entropy/VonNeumann.lean#L68).

**Matrix helpers in Physlib QuantumInfo.**
- `HermitianMat.trace_cfc_eq` [ForMathlib/HermitianMat/CFC.lean L1278](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/QuantumInfo/ForMathlib/HermitianMat/CFC.lean#L1278)
- operator `log_concave` [HermitianMat/LogExp.lean L621](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/QuantumInfo/ForMathlib/HermitianMat/LogExp.lean#L621)
- joint-concavity results in `ForMathlib/HayataGroup/TraceInequality/`
- `ForMathlib/Majorization.lean` (singular values, `cauchyBinet`)
- `ForMathlib/HermitianMat/Peierls.lean` contains `trace_cfc_conj_unitary`, a trace lemma. In the Lean-QuantumInfo snapshot, `peierls_inequality` (a trace inequality for convex g, not the Peierls contour argument) sits in [DPI.lean L171](https://github.com/Timeroot/Lean-QuantumInfo/blob/56e83a9288a3c616285038748e273b3c0e1a36bf/QuantumInfo/Finite/Entropy/DPI.lean#L171).
- I found no theorem whose name mentions concavity of `Sᵥₙ` (grep over `QuantumInfo`).

**Courant–Fischer, Weyl and Cauchy interlacing (outside Mathlib).** [PhyslibAlpha/Mathematics/TauCeti/Analysis/InnerProductSpace/CourantFischer.lean](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/PhyslibAlpha/Mathematics/TauCeti/Analysis/InnerProductSpace/CourantFischer.lean) (645 lines) contains:
- `LinearMap.IsSymmetric.eigenvalues_eq_iSup_iInf_re_inner` (Courant–Fischer)
- `eigenvalue_mono` (monotonicity under quadratic-form order)
- `abs_eigenvalue_sub_eigenvalue_le` and `abs_eigenvalue_sub_eigenvalue_le_norm` (Weyl)
- `eigenvalue_compression_le` and `eigenvalue_le_eigenvalue_compression` (Cauchy interlacing for compressions)
- `eigenvalue_ascending_le_ritz`.

**Wick's theorem (Physlib).**
- `wicks_theorem : 𝓣(ofFieldOpList φs) = …` [WickAlgebra/WicksTheorem.lean L79](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/Physlib/QFT/PerturbationTheory/WickAlgebra/WicksTheorem.lean#L79)
- `wicks_theorem_normal_order` (WicksTheoremNormal.lean L213) and `static_wick_theorem` (StaticWickTheorem.lean L57)
- the `Physlib/QFT/PerturbationTheory/` tree (FieldSpecification, FieldStatistics, Koszul, WickContraction/Sign, FieldOpFreeAlgebra) contains 0 `sorry` tokens — [directory](https://github.com/leanprover-community/physlib/tree/56430c65b76f73c7c25fc03eed46be1d888520bb/Physlib/QFT/PerturbationTheory).

**Condensed matter.**
- [TightBindingChain/Basic.lean](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/Physlib/CondensedMatter/TightBindingChain/Basic.lean) (633 lines): a single electron on an N-site periodic chain, with `localizedState`, `hamiltonian`, `BrillouinZone`, `energyEigenstate`, `energyEigenvalue` and `hasEigenvalue_hamiltonian_iff`.
- These files are stubs:
  - `CondensedMatter/LatticeModels/Basic.lean` (18 lines; "a home for hopping Hamiltonians, tight-binding and Hubbard-type models, and spin lattice models")
  - `CondensedMatter/ManyBody/Basic.lean` (17 lines)
  - `LatticeQFT/Basic.lean` (22 lines) — [LatticeModels](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/Physlib/CondensedMatter/LatticeModels/Basic.lean).

**Ising model.** [PhyslibAlpha/CondensedMatter/IsingModel/Basic.lean](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/PhyslibAlpha/CondensedMatter/IsingModel/Basic.lean) (353 lines) covers classical Ising energy on a graph, boundary/volume bounds, the ferromagnetic ground state, and convergence of energy per site. Its docstring does not mention a phase transition or a Peierls argument.

**Other PhyslibAlpha material (file names only).** PhyslibAlpha also has `CondensedMatter/QuasiLocal/{Gibbs, GibbsSymmetry, Interaction, LimitState, …}` and `Bloch/`, `Crystal/` — [directory](https://github.com/leanprover-community/physlib/tree/56430c65b76f73c7c25fc03eed46be1d888520bb/PhyslibAlpha/CondensedMatter).

**Measurement.**
- [PhyslibAlpha/ProbabilisticTheory/Measurement/BornRule.lean](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/PhyslibAlpha/ProbabilisticTheory/Measurement/BornRule.lean): `EffectValuedMeasure.toProbabilityMeasure` and `EffectValuedMeasure.probabilityLaw`.
- Also present (names only): `Measurement/Instrument.lean`, `JordanOrderUnit/Lueders.lean` and `QuantumInfo/Measurements/POVM.lean`.

**Euclidean field theory.** [Physlib/QFT/Euclidean/GlimmJaffe/](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/Physlib/QFT/Euclidean/GlimmJaffe/Basic.lean) gives a measure-based real scalar Euclidean field with Schwinger functions and the Glimm–Jaffe axioms (analyticity, growth, Euclidean invariance, reflection positivity, ergodicity).

**Negative keyword search.** Case-insensitive greps over `Physlib/`, `PhyslibAlpha/` and `QuantumInfo/` for Jordan–Wigner, Lieb–Robinson, "correlation matrix", "free fermion", percolation and Slater gave no relevant hits — [repository](https://github.com/leanprover-community/physlib/tree/56430c65b76f73c7c25fc03eed46be1d888520bb). The two apparent hits were false positives: "Slater" matched inside the identifier `translateRegion`, and "free fermion" matched the identifier `fermionicProjFree`.

**AI policy.** [AI-POLICY.md](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/AI-POLICY.md):
- "Physlib generally welcomes AI-assisted contributions" (1.1)
- "A human must vouch that each definition, theorem statement, and proof step means what it claims" (1.7)
- reviewer communication must be by humans (3.1)
- an `AGENTS.md` file holds agent instructions.

#### Other 2026 quantum-information libraries in Lean
- **Lean-Quantum** (arXiv 2607.05492, July 2026; Kasaura, …, Yamasaki) — [arXiv](https://arxiv.org/abs/2607.05492):
  - basis-independent framework; the data-processing inequality (DPI) for the sandwiched Rényi relative entropy, with SSA as a corollary
  - claims "the last missing component needed to complete the Lean formalization of the generalized quantum Stein's lemma" (snippet)
  - repository [Hayata-Yamasaki-Group/lean-quantum](https://github.com/Hayata-Yamasaki-Group/lean-quantum): HEAD 28 Oct 2026, 27,444 Lean lines, 34 `sorry` tokens at HEAD; its README points to a `paper2` branch for the published version.
- **Lean-QIT** (arXiv 2607.09632; QuAIR / University of Hong Kong) — [arXiv](https://arxiv.org/abs/2607.09632), [news](https://quantumzeitgeist.com/qudeleap-research-hkust-lean-library-quantum/):
  - separates operational definitions (codes, rates, capacities) from analytic characterizations
  - validated by formalizing the Holevo–Schumacher–Westmoreland (HSW) classical capacity theorem
  - a follow-up benchmark, Lean-QIT-Bench, has 40 tasks (snippet)
  - the repository [QuAIR/Lean-QIT](https://github.com/QuAIR/Lean-QIT) exists.
- **Further 2026 titles (content not verified):** [QECert (2611.05047)](https://arxiv.org/abs/2611.05047), [End-to-End Formalization of Quantum Error Correction (2605.16523)](https://arxiv.org/abs/2605.16523), [AxQM benchmark (2609.05157)](https://arxiv.org/abs/2609.05157), [Benchmarking Agents for Proving Theorems in Quantum Algorithms and Quantum Information (2607.21533)](https://arxiv.org/abs/2607.21533), [Lean-QuantumAlg (2610.30216)](https://arxiv.org/abs/2610.30216).

#### formal_conjectures (Google DeepMind)
- At `2207c54d` (21 Nov 2026) the repository has 2,154 `.lean` files. Top folders: Arxiv, Books, ErdosProblems, GreensOpenProblems, HilbertProblems, Kourovka, LittProblems, Mathoverflow, Millennium, OEIS, **OpenQuantumProblems**, OptimizationConstants, Other, Paper, Subsets, Wikipedia, WrittenOnTheWallII — [repository](https://github.com/google-deepmind/formal-conjectures/tree/2207c54d0dbe960de591d280f997b195bf3f6e28).
- `OpenQuantumProblems/` contains only problems 13 (mutually unbiased bases), 23 (SIC-POVMs) and 35 (absolutely maximally entangled states) from the IQOQI Vienna list — [folder](https://github.com/google-deepmind/formal-conjectures/tree/2207c54d0dbe960de591d280f997b195bf3f6e28/FormalConjectures/OpenQuantumProblems).
- Greps find no percolation, Ising, random-matrix, Lyapunov or self-avoiding-walk content; "random walk" appears only in Erdős problems 1164–1166 — [repository](https://github.com/google-deepmind/formal-conjectures/tree/2207c54d0dbe960de591d280f997b195bf3f6e28).
- A benchmark paper describes the repository: "Formal Conjectures: An Open and Evolving Benchmark for Verified Discovery in Mathematics" (snippet) — [arXiv 2605.13171](https://arxiv.org/abs/2605.13171).

#### Percolation, Peierls, random matrices, Lieb–Robinson

**The anthropics/formal-math percolation artifact.** At commit `795efb86` (28 Aug 2026) of `anthropics/formal-math` there is a `percolation/` Lake project claiming **θ(p_c)=0 for nearest-neighbour Bernoulli bond percolation on ℤ^d for every d ≥ 2** — [README](https://github.com/anthropics/formal-math/blob/795efb86f191735c5481675763537cfb4ff37e55/percolation/README.md).
- *Claimed route:* prove the Kozma–Nitzan Conjecture 3 via an "additive gluing" inequality, then use a re-proof of their Theorem 6.
- *Classical inputs re-proved:* Harris, van den Berg–Häggström–Kahn, the four-functions theorem, Russo–Seymour–Welsh (RSW), Kesten's p_c(ℤ²)=½, and Barsky–Grimmett–Newman.
- *Size (my count):* 251 Lean files, 97,574 lines, 3 `sorry` tokens. The README says there are two deliberate placeholders in `Challenge.lean`.
- *Model:* the percolation measure uses Mathlib's `ProbabilityTheory.setBernoulli` — [Challenge.lean](https://github.com/anthropics/formal-math/blob/795efb86f191735c5481675763537cfb4ff37e55/percolation/Challenge.lean).
- *Peierls ingredient:* a 2D Peierls dual-circuit lemma `Percolation.Literature.Contour.exists_dualCircuit` in [DualContours.lean](https://github.com/anthropics/formal-math/blob/795efb86f191735c5481675763537cfb4ff37e55/percolation/Percolation/Literature/DualContours.lean) (1,082 lines). The guide says "the library runs an explicit Peierls contour estimate with ε ≤ 2^{-32}" — [summary.tex L378](https://github.com/anthropics/formal-math/blob/795efb86f191735c5481675763537cfb4ff37e55/percolation/summary.tex#L378).
- *Provenance (README):* the Lean sources "were written by an AI system (Anthropic's Claude models) working autonomously under the direction of Justin Leder; no human wrote or edited the Lean code". The work "has not yet been refereed by human mathematicians or by anyone independent of the author" — [README](https://github.com/anthropics/formal-math/blob/795efb86f191735c5481675763537cfb4ff37e55/percolation/README.md).
- *Toolchain:* Lean `v4.32.0`.

**Status of that commit.**
- As of 22 Nov 2026 it is **not on `main`**: it is absent from `main`'s history since 1 Aug 2026, and it is not the tip of any branch, tag or pull-request ref (`git ls-remote`).
- The [`main` README](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/README.md) lists only `zeta23`, `3sum-apsp` and `eth`. I found no public explanation for the removal.

**Other percolation-adjacent work.**
- "Iterated Graph Systems (II)" (arXiv 2610.10469) reports a Lean 4/Mathlib development on Bernoulli percolation on hierarchical lattices (snippet) — [arXiv](https://arxiv.org/abs/2610.10469).
- [lean-stat-learning-theory](https://github.com/YuanheZ/lean-stat-learning-theory) (ICML 2026; repository exists) claims a sorry-free library including matrix Bernstein, Lieb's concavity inequality and Hanson–Wright (snippet).
- arXiv 2609.13331 (random permutation matrices) reports a Lean 4 formalization of its finite model with no placeholders (snippet) — [arXiv](https://arxiv.org/abs/2609.13331).
- [j49liu/homogeneous-lyapunov-counterexample-lean](https://github.com/j49liu/homogeneous-lyapunov-counterexample-lean) is about Lyapunov *functions* for ODE stability, not exponents.

**Nothing found.** A web search found no Lean formalization of Lieb–Robinson bounds, Lyapunov exponents of random-matrix products, Furstenberg's theorem, or fermionic Gaussian states/Slater determinants in any proof assistant. The closest analogue was a bosonic free-field QFT formalization following Glimm–Jaffe (arXiv 2603.15770, snippet) — [search result](https://arxiv.org/abs/2603.15770).

### Inferences
- **Entropy statements should target Physlib.** The natural dependency for von Neumann entropy statements is Physlib's `QuantumInfo` (`MState d`, `Sᵥₙ`, `traceLeft`/`traceRight`). Using it for fermions requires a Jordan–Wigner identification of Fock(A∪B) with `d_A × d_B`. This is easiest when A is a contiguous block at the start of the chain.
- **Physlib's Wick theorem is not the Wick theorem needed here.** It is an operator-algebra statement (time and normal ordering in a field-operator algebra), not "Gaussian-state expectations factorize into pairings". Its contraction and sign combinatorics might be reusable, but the bridge is non-trivial.
- **Courant–Fischer and Cauchy interlacing exist in PhyslibAlpha** (at the operator level). They are directly useful for eigenvalues of C_A, but they live in an "Alpha" tier and are not in Mathlib.
- **The percolation artifact is a reuse source, not a citable result.** Its `setBernoulli`-based bond model, Harris/FKG and contour code are a possible source to port: the repository is Apache-2.0 and the folder has its own LICENSE and NOTICE files. However, it targets a different (2D dual-circuit) Peierls argument, sits on Lean `v4.32.0`, is unreviewed and is orphaned. Its headline claim must be treated as unverified.
- **Version skew means porting costs.** Toolchains range from `v4.28` (Lean-QuantumInfo) and `v4.32` (percolation) to `v4.36` (Physlib) and `v4.37-rc1` (Mathlib master); any reuse of older code needs porting.

### Gaps
- The `sorry` status is based on token grep only; no `#print axioms` was run.
- Concavity of `Sᵥₙ` may exist under a name not matching "concave".
- The reason the percolation commit is unreferenced, and whether it was retracted, is unknown.
- The contents of Lean-QIT, QECert and the related papers were not verified, because arXiv was unreachable.
- PhysLean's earlier `informal_lemma` placeholders and Physlib's TODO list were not reviewed.

## 3. Other proof assistants: Isabelle/HOL and Coq/Rocq quantum-information formalizations and entropy inequalities

### Takeaway
The Coq/Rocq and Isabelle quantum libraries found in this session (CoqQ, QWIRE, Isabelle Marries Dirac, certified quantum computation in Isabelle) target quantum programs and circuits. I found no evidence of machine-checked quantum entropy inequalities such as SSA outside Lean. Lean's Physlib/QuantumInfo appears to be ahead on quantum Shannon theory, but this is not established.

### Cited Findings
- **CoqQ** is a Coq framework for verifying quantum programs, restricted to finite-dimensional Hilbert spaces through MathComp's vector-space library, with a deeply embedded language and labelled Dirac notation (snippet) — [arXiv 2207.11350](https://arxiv.org/abs/2207.11350).
- **QWIRE** is a Coq-embedded circuit language with density-matrix semantics and no program logic of its own, which is one of the gaps CoqQ addresses (snippet) — [Formally Verified Quantum Programming](https://rand.cs.uchicago.edu/publication/rand-2018-formally/).
- **Isabelle Marries Dirac (IMD)** is an Isabelle/HOL library built on matrix representations of quantum circuits, covering no-cloning, teleportation, Deutsch and Deutsch–Jozsa. A survey notes its Dirac-notation support is limited (snippet) — [AFP entry](https://isa-afp.org/entries/Isabelle_Marries_Dirac.html); [Certified Quantum Computation in Isabelle/HOL](https://arxiv.org/abs/2012.13925).
- A survey of the area exists: "Formal Verification of Quantum Programs: Theory, Tools and Challenges" (snippet) — [arXiv 2110.01320](https://arxiv.org/abs/2110.01320).
- A targeted search found no source describing a formalized quantum entropy inequality in CoqQ, QWIRE or IMD (search summary) — [CoqQ](https://arxiv.org/abs/2207.11350).

### Inferences
- For entropy inequalities and Gaussian-state entanglement, the other assistants offer no ready infrastructure. Lean plus Physlib is the only realistic platform among those examined.

### Gaps
- Not checked, because the search budget was exhausted:
  - Isabelle AFP classical information theory
  - Coq `infotheo` (classical Shannon theory)
  - SQIR/VOQC
  - the AFP "projective measurements" and CHSH entries
  - any HOL Light or Metamath quantum work.

## 4. AI-assisted theorem proving in 2025–2026, and documented research-level physics or probability formalized with AI help

### Takeaway
- **Large AI-produced developments are now routine under human direction.** Math Inc's Gauss produced about 25k lines (strong PNT) in three weeks; secondary reports describe 200k+ lines for sphere packing in two weeks. Harmonic's Aristotle produced Lean proofs for Erdős problems. Claude models autonomously wrote a 97.6k-line percolation library.
- **Competition benchmarks are nearly saturated.** miniF2F is at ≈99%, and Seed-Prover solved 331 of 657 PutnamBench problems.
- **Quantum information saw several AI-assisted Lean libraries in 2026.**
- **The bottleneck is statement fidelity and design, not proof search.** Every project emphasizes human checking of statements; the Gauss discussion explicitly mentions misformalized lemmas.

### Cited Findings

#### Gauss and Aristotle
- **Gauss, strong PNT.** The math-inc/strongpnt README calls it "an AI-generated Lean formalization of the strong Prime Number Theorem". "Most of the statements and proofs were produced by Gauss … with targeted human scaffolding and review of key lemmas and strategies". Scale is "≈25k lines of Lean and ≈1.1k theorems/definitions", with "time to completion: three weeks". It reuses parts of PrimeNumberTheoremAnd. My count is 28,407 Lean lines and 0 `sorry` tokens; the last commit is 10 Sep 2025 — [math-inc/strongpnt](https://github.com/math-inc/strongpnt).
- **Statement problems in the Gauss output (snippet).** A Zulip discussion notes that some statements were modified (for example, an exponent change), and that the agent sometimes skipped misformalized lemmas and sometimes used them — [Zulip archive](https://leanprover-community.github.io/archive/stream/219941-Machine-Learning-for-Theorem-Proving/topic/Math.2C.20Inc.2C.20and.20Gauss.20discussion.html).
- **Gauss, sphere packing.** A March 2026 Math Inc document is titled "Completing the formal proof of higher-dimensional sphere packing" — [gwern mirror](https://gwern.net/doc/math/2026-03-mathinc-completingtheformalproofofhigherdimensionalspherepacking.html). IEEE Spectrum reportedly says Gauss autoformalized Viazovska's 24-dimensional proof, "200,000+ lines", in two weeks (snippet; unverified) — [IEEE Spectrum](https://spectrum.ieee.org/ai-proof-verification).
- **Aristotle, Erdős #728.** In January 2026, a GPT-5.2 Pro + Aristotle workflow resolved Erdős Problem #728 with a Lean proof, reported as the first Erdős problem regarded as fully resolved autonomously by AI. The problem's original statement was vague (snippet) — [arXiv 2601.07421](https://arxiv.org/abs/2601.07421).
- **Aristotle, Sidon sets.** In a paper on multiplicative Sidon sets, Theorem 1.3 was "autonomously discovered and formally verified in Lean by Aristotle" (snippet) — [arXiv 2605.02064](https://arxiv.org/abs/2605.02064).

#### Benchmarks
All numbers are self-reported by each team and were not rerun on a common setup.

| System | miniF2F-test | PutnamBench | Source |
|---|---|---|---|
| DeepSeek-Prover-V2 (671B, Apr 2025) | 88.9% | 47 (pass@1024) | [2504.21801](https://arxiv.org/abs/2504.21801) |
| Goedel-Prover-V2-32B (Aug 2025) | 88.1% pass@32; 90.4% with self-correction | 86 (pass@184, self-correction) | [2508.03613](https://arxiv.org/abs/2508.03613) |
| Seed-Prover (ByteDance) | 99.6% | 331/657 (medium setting) | [review page](https://liner.com/review/seedprover-deep-and-broad-reasoning-for-automated-theorem-proving) |
| Kimina-Prover | 92.2% (as quoted by Seed-Prover) | — | same review page |
| HILBERT | 99.2% | — | [2026 survey, 2607.07779](https://arxiv.org/abs/2607.07779) |

#### Long-horizon research formalization
- **LeanMarathon** (2026) formalized all 7 target theorems of two papers, covering 4 Erdős problems, with no `sorry` (snippet) — [arXiv 2606.05400](https://arxiv.org/abs/2606.05400).
- **AIProver** is described as "Agentic Auto-Formalization of Mathematical Research" (title only) — [arXiv 2610.05367](https://arxiv.org/abs/2610.05367).

#### Physics, quantum and probability examples
- Lean-Quantum's title is "Toward AI-Assisted Formalization of Quantum Information"; the AI tools used were not visible in the snippet — [arXiv 2607.05492](https://arxiv.org/abs/2607.05492).
- Lean-QIT presents itself as a base "for AI-assisted formalization and proof search" (snippet) — [arXiv 2607.09632](https://arxiv.org/abs/2607.09632).
- Physlib's AI policy welcomes AI-assisted pull requests as long as a human vouches for them — [AI-POLICY.md](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/AI-POLICY.md).
- The percolation artifact (§2) is AI-written: 97,574 lines, with "the only review so far … carried out by AI systems" — [README](https://github.com/anthropics/formal-math/blob/795efb86f191735c5481675763537cfb4ff37e55/percolation/README.md).

#### anthropics/formal-math `main`
- The `main` README lists three other AI-era projects. I did not evaluate them; one of them, the claim to refute the Exponential Time Hypothesis (ETH), is extraordinary and should not be relied on:
  - `zeta23`: more than two thirds of the zeros of zeta are simple and on the critical line (Alpöge–Furman, arXiv 2608.13637)
  - `3sum-apsp`
  - `eth`, which claims a refutation of the ETH, with the paper credited to "Claude"
- Its CI checks each project with `lake build` and with the Palomar registry's `leanprover/comparator`, an independent `nanoda` kernel and a `landrun` sandbox, against a trusted `Challenge.lean` statement file — [README](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/README.md).

### Inferences
- **What AI changes for our targets.** Finite-dimensional linear algebra and one-variable real inequalities (items 1–3, plus the scalar part of item 2) are exactly the profile where 2026 tools (Aristotle-class provers, agentic Claude/GPT workflows) compress proof-writing several-fold. The remaining human cost is:
  - choosing definitions (Fock space, circuit model, trajectory measure)
  - checking that statements mean what is intended.
- **The comparator/`Challenge.lean` pattern is a good model for item 6.** It separates a small, human-audited statement file from a large, possibly AI-written proof.

### Gaps
- AlphaProof's 2025 Nature paper and current Aristotle, Seed-Prover 1.5/2 and Gauss numbers were not retrieved; the search budget was exhausted.
- The IEEE Spectrum and gwern claims are unverified.
- No time or cost data per line was found for AI-assisted physics formalization.

## 5. Practical setup: which network hosts a Lean 4 + Mathlib project needs (elan, toolchain, `lake exe cache get`), and how long a from-source build takes

### Takeaway
A fresh Lean + Mathlib setup in late 2026 needs:
- `elan.lean-lang.org` (installer)
- `github.com` plus GitHub's release-download redirect hosts (elan binary and toolchain archives; the exact redirect host is not traced)
- `release.lean-lang.org` (toolchain release index read by elan)
- `cache.mathlib.org` (Mathlib build cache)
- optionally `raw.githubusercontent.com`, `api.github.com`, `reservoir.lean-lang.org` (Reservoir-style `require`) and `lake-cache.physlib.io` (if depending on Physlib).

In this container, `elan.lean-lang.org`, `release.lean-lang.org`, `cache.mathlib.org` and `reservoir.lean-lang.org` are blocked, while the GitHub hosts work.

Building Mathlib from source without the cache took about 30–50 min on CI-class machines and about 2 h or more on personal computers (2024–2025 anecdotes). Mathlib now has 8,730 files, so expect more.

### Cited Findings

#### elan
- Installer command: `curl https://elan.lean-lang.org/elan-init.sh -sSf | sh` — [elan README L37](https://github.com/leanprover/elan/blob/master/README.md#L37).
- The installer downloads elan from `ELAN_UPDATE_ROOT="https://github.com/leanprover/elan/releases"` — [elan-init.sh L17](https://github.com/leanprover/elan/blob/master/elan-init.sh#L17).
- Self-update uses `https://github.com/leanprover/elan/releases/download`.
- For the `leanprover/lean4` origin (and its nightly variant), elan fetches the release index JSON from `https://release.lean-lang.org` (`DEFAULT_ORIGIN_JSON_URL`) and downloads the asset's `browser_download_url` — [manifestation.rs L15–16](https://github.com/leanprover/elan/blob/master/src/elan-dist/src/manifestation.rs#L16).
- For other origins it uses `https://github.com/{origin}/releases/expanded_assets/{tag}` and `https://api.github.com/` — [utils.rs L608](https://github.com/leanprover/elan/blob/master/src/elan-utils/src/utils.rs#L608).
- It reads remote `lean-toolchain` files from `raw.githubusercontent.com` — [toolchain.rs L149](https://github.com/leanprover/elan/blob/master/src/elan/toolchain.rs#L149).
- All of this is from elan `master`, commit dated 18 Nov 2026.

#### Lake and Reservoir
- Lake's default Reservoir endpoint is `RESERVOIR_API_BASE_URL = https://reservoir.lean-lang.org/api` — [Lake Config/Env.lean L169](https://github.com/leanprover/lean4/blob/master/src/lake/Lake/Config/Env.lean#L169) (lean4 `master`, 22 Nov 2026).

#### Mathlib cache
- `lake exe cache get` reads from `publicCacheEndpoint = "https://cache.mathlib.org"` — [Cache/Infra.lean L120](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Cache/Infra.lean#L120).
- The read base can be overridden with `MATHLIB_CACHE_BASE_URL`, or `MATHLIB_CACHE_GET_URL` for a single flat endpoint.
- The legacy path is retired: "MATHLIB_CACHE_DEBUG_USE_LEGACY is retired; cache reads use https://cache.mathlib.org" — [Cache/Main.lean L128](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Cache/Main.lean#L128).
- If the system curl is too old, on Linux x86_64/aarch64 the cache tool downloads a static curl from `github.com/leanprover-community/static-curl/releases` — [Cache/IO.lean L302](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Cache/IO.lean#L302).

#### Mathlib dependencies and toolchain
- Mathlib's dependencies are fetched by git, per `lake-manifest.json`: plausible, LeanSearchClient, importGraph, proofwidgets, aesop, Qq, batteries and Cli — [lake-manifest.json](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/lake-manifest.json).
- The README advises running `lake exe cache get` ("Skipping this step means the next step will be very slow") — [README L58](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/README.md#L58).
- Mathlib master's toolchain is `leanprover/lean4:v4.37.0-rc1`. Physlib pins `v4.36.0`, both for Lean and for Mathlib — [Mathlib lean-toolchain](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/lean-toolchain), [Physlib lakefile](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/lakefile.toml#L12).

#### Physlib cache
- Physlib uses Lake's built-in cache (`LAKE_CONFIG=$PWD/lake-cache.toml lake cache get`). Its public read endpoint is `https://lake-cache.physlib.io/artifacts` (and `/revisions`); uploads go to a Cloudflare R2 S3 endpoint — [lake-cache.toml](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/lake-cache.toml).

#### Reachability from this container (22 Nov 2026)
- I ran `curl` probes and `git clone`s.
- Connection failed (`000`): `https://elan.lean-lang.org/elan-init.sh`, `https://release.lean-lang.org`, `https://cache.mathlib.org` and `https://reservoir.lean-lang.org`.
- Reachable:
  - `raw.githubusercontent.com` (HTTP 200)
  - `release-assets.githubusercontent.com` and `objects.githubusercontent.com` (connections succeeded; HTTP 404 at the root)
  - `git clone` from `github.com`.
- WebFetch could not resolve `arxiv.org`, `leanprover-community.github.io` or `spectrum.ieee.org`.
- These are my own observations (Bash logs), so there is no external URL. The hosts involved come from the [elan](https://github.com/leanprover/elan/blob/master/src/elan-dist/src/manifestation.rs#L16) and [Mathlib cache](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Cache/Infra.lean#L120) sources above.

#### From-source build time (snippets of Zulip threads)
- "about 30–45 minutes" on community CI machines.
- "50 minutes in GitHub Actions" (Nov 2025).
- "about 2 hours to build mathlib on my computer" (2024); in 2023 one build had not finished after two hours.
- The community advice is to use the cache instead — [Zulip: Long compilation of mathlib](https://leanprover-community.github.io/archive/stream/113488-general/topic/Long.20compilation.20of.20mathlib.html), [Zulip: How long does it take to initially build a mathlib proj?](https://leanprover-community.github.io/archive/stream/113489-new-members/topic/How.20long.20does.20it.20take.20to.20initially.20build.20a.20mathlib.20proj.3F.html), [Zulip: Cache on GitHub Actions](https://leanprover-community.github.io/archive/stream/113488-general/topic/Cache.20on.20GitHub.20Actions.html).

### Inferences
**Minimal allowlist for this container:**
- `elan.lean-lang.org` (or use the reachable `raw.githubusercontent.com/leanprover/elan/master/elan-init.sh`)
- `release.lean-lang.org`
- `cache.mathlib.org`
- `reservoir.lean-lang.org`, only if using `require "leanprover-community" / "mathlib"`; a git `require mathlib from git "https://github.com/leanprover-community/mathlib4" @ "<rev>"` avoids it
- `lake-cache.physlib.io`, if depending on Physlib
- the GitHub hosts, which already work.

**Without the cache**, budget at least 1 h on a 16+ core machine and several hours on a laptop for Mathlib's ~8.7k files. Building Physlib (~450k further lines) on top adds substantially more. This is an extrapolation, not a measurement.

### Gaps
- I did not trace which GitHub host serves release redirects in 2026 (`release-assets.githubusercontent.com` versus `objects.githubusercontent.com`).
- I did not check whether ProofWidgets downloads prebuilt JavaScript from GitHub releases.
- No current (late-2026) measured from-source build time was found.
- The Zulip threads could not be opened directly.

## 6. Effort estimates for the candidate targets: projective measurement update (1), entanglement entropy and its variance bound (2), rank bounds (3), Peierls-type area law (4), Born-rule trajectory measure (5), formal statement of the transition (6)

### Takeaway
- **Matrix-level targets (1)–(3) are realistic short projects.** Each needs roughly 1–3 person-weeks (pw) by hand, or a few days with 2026-grade AI assistance, on current Mathlib.
- **Definitional targets (5)–(6) need about 2–5 pw**, mostly design work.
- **Target (4), a uniform-in-L area law at high measurement probability via a Peierls/path-counting union bound, is a plausible multi-month formalization.** It becomes much easier if one uses a pathwise "measurement kills the row of C_AB" lemma, derived below and still to be checked, which makes the bound independent of Born outcomes.
- **The genuinely expensive layer is the many-body bridge**, i.e. proving that the correlation-matrix rules equal Fock-space quantum mechanics and that S_A equals the von Neumann entropy of the reduced state: about 2–5 person-months.
- **Proving that a transition exists or does not exist is out of reach;** only stating it is realistic.

### Cited Findings

#### Calibration data
- Strong PNT took about 25k lines and three weeks with Gauss plus human scaffolding — [strongpnt](https://github.com/math-inc/strongpnt).
- Lean-QuantumInfo: 38k lines by May 2026, starting from a 2024 repository (per its citation entry). GQSL went from its first milestone in October 2025 to completion in April 2026 — [README](https://github.com/Timeroot/Lean-QuantumInfo/blob/56e83a9288a3c616285038748e273b3c0e1a36bf/README.md).
- The AI-written percolation library has 97.6k lines; its development time is not stated — [README](https://github.com/anthropics/formal-math/blob/795efb86f191735c5481675763537cfb4ff37e55/percolation/README.md).
- Statement-only formalization of open problems is an established deliverable (formal_conjectures, including OpenQuantumProblems) — [folder](https://github.com/google-deepmind/formal-conjectures/tree/2207c54d0dbe960de591d280f997b195bf3f6e28/FormalConjectures/OpenQuantumProblems).

#### Infrastructure available per target (verified names, §1–§2)
- **(1):**
  - `IsStarProjection` — [StarProjection.lean](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Algebra/Star/StarProjection.lean#L73)
  - `Matrix.single` — [Basis.lean](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Data/Matrix/Basis.lean#L34)
  - Woodbury and the matrix-determinant lemma — [NonsingularInverse](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/LinearAlgebra/Matrix/NonsingularInverse.lean#L616).
- **(2):**
  - `Real.binEntropy` and its derivatives and concavity — [BinaryEntropy.lean](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/SpecialFunctions/BinaryEntropy.lean#L366)
  - `convexOn_of_deriv2_nonneg` and `concaveOn_of_deriv2_nonpos` — [Convex/Deriv.lean](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/Convex/Deriv.lean#L215)
  - `log_two_gt_d9` — [ExponentialBounds](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/Complex/ExponentialBounds.lean#L86)
  - `IsHermitian.eigenvalues` and `trace_eq_sum_eigenvalues` — [Spectrum.lean](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/Matrix/Spectrum.lean#L240)
  - `PosSemidef.submatrix` — [PosDef.lean](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/LinearAlgebra/Matrix/PosDef.lean#L84)
  - `Sᵥₙ` and SSA — [Physlib SSA.lean](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/QuantumInfo/Entropy/SSA.lean#L1191)
  - interlacing — [CourantFischer.lean](https://github.com/leanprover-community/physlib/blob/56430c65b76f73c7c25fc03eed46be1d888520bb/PhyslibAlpha/Mathematics/TauCeti/Analysis/InnerProductSpace/CourantFischer.lean).
- **(3):** `rank_eq_card_non_zero_eigs`, `rank_eq_rank_diagonal`, `rank_self_mul_conjTranspose` and `rank_mul_le` — [Spectrum.lean L202](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Analysis/Matrix/Spectrum.lean#L202), [Rank.lean L549](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/LinearAlgebra/Matrix/Rank.lean#L549).
- **(4):**
  - `setBernoulli` — [SetBernoulli.lean](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Probability/Distributions/SetBernoulli.lean#L44)
  - `Measure.infinitePi` — [ProductMeasure.lean](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Probability/ProductMeasure.lean#L354)
  - `SimpleGraph.Walk` — [SimpleGraph](https://github.com/leanprover-community/mathlib4/tree/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Combinatorics/SimpleGraph)
  - the 2D contour code in the orphaned artifact — [DualContours.lean](https://github.com/anthropics/formal-math/blob/795efb86f191735c5481675763537cfb4ff37e55/percolation/Percolation/Literature/DualContours.lean).
- **(5):** `Kernel.traj`, `trajMeasure` and `Kernel.ofFunOfCountable` — [Traj.lean](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Probability/Kernel/IonescuTulcea/Traj.lean#L519), [Kernel/Basic.lean](https://github.com/leanprover-community/mathlib4/blob/a495318fed6482b3c90df686de5131b46a70fd5b/Mathlib/Probability/Kernel/Basic.lean#L236).
- **(6):** the comparator and trusted-statement pattern — [formal-math README](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/README.md).

### Inferences
*(The estimates and derivations below are my own. They need checking by the physics and Lean team. Units: pw = person-weeks of a Lean user fluent in Mathlib's linear-algebra API. "AI" = with 2026 agentic provers, with a human writing and reviewing statements.)*

**Mathematical observations that set the effort**
- **(1) Projector property.** Let C be a projector and v := C eᵢ/√Cᵢᵢ. Then ‖v‖² = eᵢ†C†Ceᵢ = Cᵢᵢ = 1, and v ∈ ran C. So P := C − vv† is a projector, and P eᵢ = Ceᵢ − v·(v†eᵢ) = Ceᵢ − Ceᵢ = 0 because v†eᵢ = √Cᵢᵢ. Hence C′ = P + eᵢeᵢ† is a sum of orthogonal projectors, so it is a projector, with rank C′ = rank C and tr C′ = tr C.
  - Both Born probabilities are real and in [0,1]: Cᵢᵢ = ‖Ceᵢ‖² and 1 − Cᵢᵢ = ‖(1−C)eᵢ‖².
  - The outcome-0 update is C″ = C + (1−C)eᵢeᵢ†(1−C)/(1−Cᵢᵢ) − eᵢeᵢ†, by the particle–hole symmetry C ↦ 1−C.
- **(2) Scalar inequality.** Let g(p) = binEntropy(p) − 4 ln2·p(1−p). Then g(0) = g(½) = g(1) = 0, g is symmetric about ½, and g″ = −1/(p(1−p)) + 8 ln 2.
  - g is concave on [0, p*] and convex on [p*, ½], where p* ≈ 0.236 solves p(1−p) = 1/(8 ln 2).
  - Convexity with g′(½) = 0 gives g ≥ 0 on [p*, ½]. Concavity with g(0) = 0 and g(p*) ≥ 0 gives g ≥ 0 on [0, p*].
  - This is provable with `deriv2_binEntropy`, `concaveOn_of_deriv2_nonpos`, `convexOn_of_deriv2_nonneg` and `log_two_gt_d9`.
  - Spot checks: at p = 0.1, 0.3251 ≥ 0.2495; at p = 0.4, 0.6730 ≥ 0.6654.
  - Summing over the eigenvalues ν of C_A gives S_A ≥ 4 ln2·Σν(1−ν) = 4 ln2·tr C_A(1−C_A).
- **(3) Rank identity.** Write C in blocks over A and B. Then C² = C gives C_A − C_A² = C_AB C_AB†. Hence rank C_AB = rank(C_A(1−C_A)) = #{ν : ν(1−ν) ≠ 0} = #{ν ∈ (0,1)}, using `rank_self_mul_conjTranspose` and `rank_eq_card_non_zero_eigs` applied to C_A(1−C_A), plus the spectral theorem to identify that matrix's eigenvalues with ν(1−ν). How gates and measurements change it:
  - local unitaries on A or on B preserve rank C_AB
  - a two-site gate on (a ∈ A, b ∈ B) changes C_AB by at most two rank-one terms, so rank changes by at most 2
  - a measurement of nᵢ with i ∈ A (either outcome) gives C′_AB = (1 − C_{A,i}eᵢ†/Cᵢᵢ)·C_AB for outcome 1 (similarly for outcome 0); it zeroes row i, so rank C_AB never increases.
- **(4) Pathwise mechanism.** By the row-killing in (3), each measured site erases its row (or column) of C_AB whatever the outcome. Combining this with S_A ≤ ln2·rank C_AB (from `binEntropy_le_log_two` and `binEntropy_eq_zero`), S_A is bounded *deterministically, given the measurement locations*, by the number of A-sites connected to the cut by unmeasured backward spacetime paths.
  - At high p, a union bound over paths of length n, roughly (c(1−p))ⁿ, gives E[S_A] ≤ K(p), uniformly in L and t. This needs only i.i.d. Bernoulli measurement locations; Born outcomes do not enter. The support-propagation lemma itself must be stated and proved carefully.

**Estimates**

| Target | What is missing | Human (pw) | AI-assisted (pw) |
|---|---|---|---|
| (1) Matrix level: Born probability Cᵢᵢ ∈ [0,1]; C′ (and C″) is a projector; rank and trace preserved | Only convenience lemmas for rank-one matrices | 0.5–1.5 | 0.2–0.5 |
| (1′) Many-body derivation: the Fock-space projective measurement on a Slater determinant yields C′ | Fock space with inner product (⋀ℂ^L or 2^L-dimensional Jordan–Wigner matrices); CAR; Slater determinants; Gram-determinant/Wick factorization of expectations; post-measurement state | 6–12 | 2–5 |
| (2) S_A := Σ binEntropy(νᵢ) and S_A ≥ 4 ln2·tr C_A(1−C_A) | The scalar inequality (absent from Mathlib); ν ∈ [0,1] via `PosSemidef.submatrix` or Rayleigh quotients; tr C_A(1−C_A) = Σν(1−ν) | 2–4 (scalar part 1–2) | 0.5–1.5 |
| (2′) S_A equals `Sᵥₙ` of the reduced Fock state; Var(N_A) = tr C_A(1−C_A) | Everything in (1′), plus the A|B tensor factorization (Jordan–Wigner), the reduced state of a Gaussian state, and second quantization Γ(U) | 10–20 | 4–8 |
| (3) #{ν ∈ (0,1)} = rank C_AB, plus the gate and measurement rank lemmas | Block-decomposition plumbing; possibly a `rank_add_le` substitute | 2–5 | 0.5–1.5 |
| (4) Uniform-in-L area law at high p via Peierls/path counting | Circuit model (brickwork plus Bernoulli measurement layers); support-propagation lemma; spacetime path enumeration and counting bound (absent from Mathlib); union bound and expectation assembly; uniformity in L and t | 8–16 | 3–6 |
| (5) Born-rule trajectory measure | Finite outcome spaces, kernels from history via `ofFunOfCountable`, `Kernel.traj`/`trajMeasure` (or `Measure.pi` for a finite horizon); junk values for probability-zero branches | 1–3 (+2–4 for Haar-random continuous gates) | 0.5–1 (+1–2) |
| (6) Formal Props such as `AreaLaw p := ∃ K, ∀ L t, avgEnt L t p ≤ K`, `TransitionExists := ∃ p₁ < p₂ in (0,1), ¬AreaLaw p₁ ∧ AreaLaw p₂` (or a sharp-threshold or log-law variant) | Decisions on observable, order of limits, mean versus typical behaviour, gate ensemble and boundary conditions; depends on (5) | 1–2 | 0.5–1 |

- **Proving the formal Props in (6).** The high-p half of `AreaLaw` is target (4). A low-p failure of the area law, and hence a "transition exists" theorem, is research-level and may be physically false for generic 1D monitored free fermions. That question belongs to the physics notes; I do not assess it here.
- **A minimum viable project.** Targets (1), (2), (3), (5) and (6) at matrix level, plus (4), come to roughly 15–30 pw by hand or 6–12 pw AI-assisted, i.e. about 3–7 months for one person. Adding the many-body bridge, (1′) and (2′), adds about 2–5 person-months.

### Gaps
- The pathwise support-propagation lemma behind (4) is my own derivation and is unverified. The exact threshold it yields (via the path-counting constant c) is not computed.
- No empirical person-week data exists for comparable physics formalizations; the estimates are extrapolated from the calibration points above.
- Whether Mathlib's exterior-power pairing files already give the Gram-determinant inner product needed for (1′)/(2′) was not checked; it could cut the cost of (1′) noticeably.
