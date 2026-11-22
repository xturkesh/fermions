# DFL — Lean 4 formalization of "Rigorous disorder-free localization in a ℤ₂ lattice gauge theory"

Lean 4 (v4.36.0) + Mathlib (v4.36.0) formalization of the results of the paper in `../paper/`.

## Building

```
cd lean
lake exe cache get     # downloads the prebuilt Mathlib (needs network access to cache.mathlib.org)
lake build
```

## What is formalized

All objects are finite-dimensional and defined from scratch:
- the 1D ℤ₂ lattice gauge theory of Smith–Knolle–Kovrizhin–Moessner, with genuine Jordan–Wigner
  fermions, link spins in the `σˣ` basis, Gauss operators and a boundary link;
- matrix exponentials (`NormedSpace.exp`) for the dynamics;
- Mathlib's spectral theorem for the one-body problem.

**External input.** There is one external input, **Hypothesis (EC)** (`DFL.EC`, `DFL/OneBody.lean`):
decay of the Bernoulli-averaged eigenfunction correlator of the one-dimensional Anderson model with
binary potential, uniformly in the volume. In the paper it comes from the localization literature
(Ge–Zhao; Macera–Sodin; Furstenberg's theorem and large deviations for random matrix products),
which is far beyond current Mathlib. Every Lean theorem that needs it takes `hEC : EC J h F` as an
explicit argument, so the Lean statements are exactly as conditional as the paper's. No `axiom` is
introduced anywhere.

| Paper | Lean | File |
|---|---|---|
| Model, Gauss law (Sec. 2) | `H`, `gauss`, `ann`, `cre`, `dGamma`, `hq`, `Hsec`, `Hdual` | `Basic.lean` |
| Theorem `thm:duality` (exact duality) | `duality`, `gauss_dual`, `dGamma_hq`, `hop_apply`, `numOp_dual` | `Duality.lean` |
| Corollary `cor:sector` (exact sector averaging) | `sector_average`, `sector_average_matter`, `gaugeInv_submatrix` | `Sector.lean` |
| Lemma `lem:productlaw` (product initial states) | `sector_sum_prod`, `sector_average_prod_fock`, `sector_average_zprod`, `expect_xprod_fock` | `ProductLaw.lean` |
| Proposition `prop:variants`(i) (open chain) | `H_boundary_link`, `linkOf_last`, `sector_average_open` | `Variants.lean` |
| Unitary dynamics | `evo_unitary`, `evo_submatrix`, `evo_blockDiagonal`, `expect_evo_sub_smul_one` | `Evolution.lean` |
| Free fermions: CAR, Heisenberg evolution, Wick | `ann_mul_cre_add`, `heisenberg_ann`, `two_point`, `density`, `number_variance` | `FreeDyn.lean` |
| Lemma `lem:consequences`, Hypothesis `hyp:EC` | `corrQ`, `norm_evo_le_corrQ`, `transProb_le_corrQ`, `EC` | `OneBody.lean` |
| Theorem `thm:transport` (bounded transport) | `lgt_density`, `lgt_transport` | `Main.lean` |
| Proposition `prop:LR` (sector-averaged zero-velocity LR) | `avg_norm_evo_le`, `avg_iSup_norm_evo_le`, `avg_norm_comm_le` | `LiebRobinson.lean` |
| Theorem `thm:imbalance` (persistent imbalance) | `stag_hq_anticomm`, `eigenvalue_mul_stagWeight`, `hq_eigenvalues_injective`, `imbOne_timeAvg`, `avg_imbInf_ge`, `avg_imbInf_uniform`, `lgtImbalance_eq`, `lgtImbalance_timeAvg`, `lgtImbalance_lower` | `Imbalance.lean`, `Main.lean` |
| Theorem `thm:numberentropy` (bounded number entropy) | `lgt_numberEntropy`, `lgtNumDist_isProb` | `Main.lean` |
| Lemma `lem:intentropy` (Massey) and entropy tools | `massey`, `shannon_le_log_card`, `pvar_mixture`, `binEntropy_le_two_sqrt` | `Entropy.lean` |
| Theorem `thm:arealaw`(a),(b) (area law, `S_G`) | `vonNeumann_slater_le`, `lgt_areaLaw` | `AreaLaw.lean` |
| Theorem `thm:arealaw`(c), Lemma `lem:cutlocal`(c) (definite charges) | `vonNeumann_rdmLGT_xprod`, `lgt_areaLaw_full` | `CutLocal.lean` |
| Proposition `prop:entropies`(d) (superposition = ensemble, all multi-time correlators) | `expect_initZ_eq_trace_ensemble`, `correlator_eq_ensemble`, `gaugeInv_heis` | `Ensemble.lean` |
| Remark `rem:randomwalk` (one charge of the standard ℤ₂ LGT in a random-walk potential) | `Hstd_dual_oneP`, `partialSum_bijective` | `RandomWalk.lean` |
| Lemma `lem:duhamel` | `colNorm_evo_sub_le`, `avg_colNorm_evo_sub_le` | `Duhamel.lean` |
| Proposition `prop:dimerlaw` (dimer sector law) | `dimerInit_dual`, `card_isDimer`, `sector_average_dimer` | `Dimer.lean` |
| Lemma `lem:dimercritical`, JSBS hypotheses | `dimerTM_plus_crit`, `dimerTM_minus_elliptic`, `dimerTM_commute_crit`, `jsbs_condition`, `dimerTM_minus_not_elliptic` | `DimerCritical.lean` |
| Proposition `prop:dimerimbalance`(i),(ii) | `dimer_stagWeight_identity`, `dimerImbalance_timeAvg` | `DimerImbalance.lean` |
| Theorem `thm:2dmap` for `σᶻ`-polarized gauge states, on **any finite graph** (2D square lattice included) | `FinGraph.sector_average_matter`, `FinGraph.gaussOp_comm_ham`, `FinGraph.dGamma_hq` | `Graph.lean` |
| Proposition `prop:gausscap` (U(1) quantum links) | `countInt_gauss`, `qlm_numberEntropy_le`, `qlm_numberEntropy_le_edge` | `QLM.lean` |

## Conventions and deviations from the paper

- **Indexing.** Sites and links are 0-indexed; the boundary link is `b_{N-1}`. The density wave
  `cdw N` occupies the even 0-indexed sites. The dimer state of `Dimer.lean` uses that density wave
  and fixes the odd links; it is a different but equivalent instance of the paper's
  Lemma `lem:productlaw`(e).
- **Area law constant.** The Lean area law has constant `2` instead of the paper's `2 ln 2`. The
  proof is elementary: `-t log t ≤ 2(√t - t)`.
- **Constant `K`.** The number-entropy and area-law constants are stated as `∃ K` (depending only on
  `J`, `h`, `F`), not with the paper's explicit expressions.

## Not formalized

- The one-body localization theorems themselves (Hypothesis (EC) is an input; see above). This
  includes the paper's own proof of (EC) on finite intervals, `thm:EC`, which builds on
  Macera–Sodin.
- Infinite-volume statements: `thm:infvol`, the JSBS transport lower bounds behind `thm:dimer`, and
  part (iii) of `prop:dimerimbalance`. All of these need operators on `ℓ²(ℤ)`.
- In two dimensions, the flux sectors other than the trivial one (`σˣ`-polarized gauge states,
  i.e. the binary random π-flux), the strong-coupling proposition `prop:cluster`, and the ring
  variant (`prop:variants`(ii)).
- The pinching inequality (`lem:pinching`) and the parts of the entropy taxonomy (`prop:entropies`)
  that need concavity of the von Neumann entropy.
- Wilson links (`prop:wilson`); conjectures.
