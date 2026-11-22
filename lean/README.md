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
  fermions, link spins in the `σˣ` basis, Gauss operators and a boundary link (also the open chain
  and the ring), and its two-dimensional version on arbitrary finite graphs;
- matrix exponentials (`NormedSpace.exp`) for the dynamics;
- Mathlib's spectral theorem for the one-body problem;
- von Neumann and Shannon entropies.

### The localization input

The one-body localization input of the paper is **Hypothesis (EC)** (`DFL.EC`, `DFL/OneBody.lean`):
decay of the Bernoulli-averaged eigenfunction correlator of the one-dimensional Anderson model with
binary potential, uniformly in the volume. The main theorems are stated with `hEC : EC J h F` as an
explicit argument.

The paper's **proof of (EC)** (Theorem `thm:EC`) is formalized as well (`Transfer.lean`,
`Resonance.lean`, `LDTEntry.lean`, `ECProof.lean`): `DFL.ec_of_ldt` derives (EC) from `DFL.LDT λ`,
`λ = -2h/J`, the **uniform large-deviation estimate** for products of i.i.d. Bernoulli transfer
matrices (vector form) together with a uniformly positive Lyapunov exponent on
`I_λ = [-2-|λ|, 2+|λ|]`. This is a published theorem (Bucaj–Damanik–Fillman–Gerbuz–VandenBoom–
Wang–Zhang, Trans. AMS 2019, Prop. 3.6 with Thms. 2.1 and 2.6; also Macera–Sodin 2022, Prop. 1.2)
whose proof (Furstenberg's theory of random matrix products) is far beyond current Mathlib; the
hypotheses of Furstenberg's theorem for our model are checked in `Furstenberg.lean`. Hence
`lgt_transport_of_ldt`, `lgt_numberEntropy_of_ldt`, `lgt_areaLaw_of_ldt` and
`lgtImbalance_lower_of_ldt` state the main results with (EC) replaced by this literature input.
No `axiom` is introduced anywhere.

| Paper | Lean | File |
|---|---|---|
| Definition `def:model` (the model), Gauss law (Sec. 2) | `H`, `gauss`, `ann`, `cre`, `dGamma`, `hq`, `Hsec`, `Hdual` | `Basic.lean` |
| Lemma `lem:gauss` (Gauss operators: self-adjoint involutions commuting with each other, `N_f` and `H`; 1D Gauss law `σˣ_{b_j} = (-1)^{N_{≤j}} ∏_{k≤j} G_k`, `∏_j G_j = (-1)^{N_f} σˣ_{b_{N-1}}`; sector projections of rank `2^N`) | `gauss_isSymm`, `gauss_mul_self`, `gauss_comm`, `gauss_comm_numF`, `gauss_comm_H`, `gauss_law_1d`, `gauss_prod`, `sxL_last_comm_H`, `gauss_mul_gaussProj`, `gaussProj_mul`, `sum_gaussProj`, `rank_gaussProj`, `gaussProj_comm_H` | `GaussLaw.lean` |
| Theorem `thm:duality` (exact duality) | `duality`, `gauss_dual`, `dGamma_hq`, `hop_apply`, `numOp_dual` | `Duality.lean` |
| Theorem `thm:duality`(c),(d),(e) (dictionary `U†f_jU = c_jτᶻ_j`, `U†σᶻ_{b_j}U = τᶻ_jτᶻ_{j+1}`, Wilson line to the boundary, `σˣ`, electric and hopping terms, parity; dressed fermions `c̃_j = Uc_jU† = f_j ∏_{k≥j} σᶻ_{b_k}`: CAR, gauge invariance, oddness, `c̃_x†c̃_y` = hopping with a Wilson line, `H` in dressed variables; `U†e^{-itH}U = ⊕_q e^{iht∑q}Γ(e^{-ith_q})` and the Heisenberg evolution of `c̃_x`) | `ann_dual`, `sz_dual`, `szTail_dual`, `sx_dual`, `elec_dual`, `hopTerm_dual`, `parity_dual`, `dressed_eq`, `dressed_anticomm`, `dressed_CAR`, `dressed_comm_gauss`, `dressed_odd`, `dressed_number`, `dressed_hop`, `H_dressed`, `evo_dual`, `heisenberg_dressed` | `DualityFull.lean` |
| Theorem `thm:duality`(b), spectrum (`char dΓ(A) = ∏_T (X - ∑_{α∈T} λ_α)` via the Slater eigenbasis; the spectrum of `H` with multiplicities is `⋃_q {∑_{α∈T} E_α(q) - h∑_j q_j}`) | `charpoly_dGamma`, `charpoly_H` | `DualitySpectrum.lean` |
| Corollary `cor:sector` (exact sector averaging) | `sector_average`, `sector_average_matter`, `gaugeInv_submatrix` | `Sector.lean` |
| Lemma `lem:sectornorm` (sector-diagonal operators: `X = ∑_q P_q X P_q`, `‖X‖ = max_q ‖P_q X P_q‖`, `\|Tr(ρX)\| ≤ ∑_q w_q ‖P_q X P_q‖`, weights `2^{-N}`, `⟨Ψ₀,\|X\|Ψ₀⟩ = 4^{-N} Tr\|X\| = E_q\|x_q\|`) | `gaugeInv_eq_sum`, `norm_gaugeInv`, `norm_trace_le`, `trace_mixed_gaussProj`, `trace_initZ_gaussProj`, `norm_trace_mixed_le`, `norm_expect_initZ_le`, `sectorDiag_abs`, `expect_sectorDiag_abs`, `trace_sectorDiag_abs` | `SectorNorm.lean` |
| Lemma `lem:productlaw` (product initial states) | `sector_sum_prod`, `sector_average_prod_fock`, `sector_average_zprod`, `expect_xprod_fock` | `ProductLaw.lean` |
| Proposition `prop:variants`(i) (open chain) | `H_boundary_link`, `linkOf_last`, `sector_average_open` | `Variants.lean` |
| Proposition `prop:variants`(ii) (ring: Wilson loop, twist `w`) | `ring_duality`, `gaussRing_dual_zero`, `gaussRing_dual`, `wilson_dual`, `Hring_comm_gaussRing`, `Hring_comm_wilson`, `gaussRing_prod`, `ringVec_orthonormal`, `card_ringLabels`, `Hring_ringVec`, `sector_average_ring` | `Ring.lean` |
| Unitary dynamics | `evo_unitary`, `evo_submatrix`, `evo_blockDiagonal`, `expect_evo_sub_smul_one` | `Evolution.lean` |
| Free fermions: CAR, Heisenberg evolution, Wick | `ann_mul_cre_add`, `heisenberg_ann`, `two_point`, `density`, `number_variance` | `FreeDyn.lean` |
| Lemma `lem:consequences`, Hypothesis `hyp:EC` | `corrQ`, `norm_evo_le_corrQ`, `transProb_le_corrQ`, `EC` | `OneBody.lean` |
| Lemma `lem:normalization` (`Γ h_q Γ⁻¹ = J(Δ + λ diag q)`; `Q` independent of the hopping sign; law of `Q` invariant under `λ ↦ -λ`) | `stag_conj_hq`, `stag_conj_hq_eq_smul`, `corrQ_hq_neg`, `hq_neg_h`, `avgQ_law_neg_h` | `Normalization.lean` |
| Theorem `thm:transport` (bounded transport; constant `2Σ₁(F)`) | `lgt_density`, `lgt_transport` | `Main.lean` |
| Theorem `thm:transport`(a)–(d) **with the constants of the paper** (dressed anticommutator `{c̃_x(t), c̃_y†} = ∑_q U_q(t)_{xy} P_q` and its localization; one-particle moments `≤ 2M_p(F)`; transport `≤ Σ₁(F)`, second moment `≤ Θ(F) + Σ'₁(F)²`; dressed correlations) and Theorem `thm:numberentropy` (`S_N(t) ≤ ½ ln(2πe(Θ(F) + Σ'₁(F)² + 1/12))`, single sector with `Θ(F)`) | `dressed_anticomm_t`, `norm_sector_anticomm`, `avg_sup_anticomm_le`, `expect_abs_anticomm_le`, `one_particle_moments`, `sector_transport_le`, `avg_sup_transport_le`, `lgt_transport_sharp`, `lgt_transport_sq`, `dressed_corr_le`, `lgt_numberEntropy_sharp`, `avg_sup_sector_numberEntropy_le` | `TransportSharp.lean` |
| Proposition `prop:LR` (sector-averaged zero-velocity LR) | `avg_norm_evo_le`, `avg_iSup_norm_evo_le`, `avg_norm_comm_le` | `LiebRobinson.lean` |
| Proposition `prop:LR`(a),(b),(c) on the gauge-theory Hilbert space (dressed monomials and Majorana expansions with `8\|A\|₁`; `B` localized in `Y` (`LocB`, contains `𝔅⁺(Y)`); sector average under (EC); expectations for density matrices with sector weights `≤ K2^{-N}`; `‖{c̃_x(t), c̃_y†}‖ = ‖[c̃_x(t), n_y]‖ = max_q \|U_q(t)_{xy}\|`) | `LR_sector_mono`, `LR_sector_maj`, `LR_avg_mono`, `LR_avg_maj`, `LR_trace_mono`, `norm_dressedAnticomm_eq`, `norm_comm_number_eq` | `LRFull.lean` |
| Proposition `prop:LR`(c), rest (`‖[n_x(t), n_y]‖ = max_q \|U\|(1-\|U\|²)^{1/2}`; clean sector: `sup_t ‖{c̃_1(t), c̃_N†}‖ ≥ (3/(2(N+1)))^{1/2}`, `sup_t ‖[n_1(t), n_N]‖ ≥ (3/(4(N+1)))^{1/2}`; sector average `≤ F(N-1)`) | `numComm`, `norm_comm_nn_eq`, `sup_norm_dressedAnticomm_ge`, `sup_norm_comm_nn_ge`, `avg_sup_anticomm_corner_le` | `LRNorm.lean` |
| Theorem `thm:imbalance` (persistent imbalance) | `stag_hq_anticomm`, `eigenvalue_mul_stagWeight`, `hq_eigenvalues_injective`, `imbOne_timeAvg`, `avg_imbInf_ge`, `avg_imbInf_uniform`, `lgtImbalance_eq`, `lgtImbalance_timeAvg`, `lgtImbalance_lower` | `Imbalance.lean`, `Main.lean` |
| Corollary `cor:imbalanceopen` (persistent imbalance on the open chain of SKKM, without boundary link) | `lgtImbalanceOpen_lower`, `lgtImbalanceOpen_timeAvg`, `lgtImbalanceOpen_eq`, `avgPar_imbInf_uniform` | `OpenImbalance.lean` |
| Theorem `thm:numberentropy` (bounded number entropy, `∃ K`; see `TransportSharp.lean` for the constant of the paper) | `lgt_numberEntropy`, `lgtNumDist_isProb` | `Main.lean` |
| Lemma `lem:intentropy` (Massey) and entropy tools | `massey`, `shannon_le_log_card`, `pvar_mixture`, `binEntropy_le_two_sqrt` | `Entropy.lean` |
| Theorem `thm:arealaw`(a),(b) (area law, `S_G`, with the constant `2` in place of `2 ln 2`; see below for the sharp form) | `vonNeumann_slater_le`, `lgt_areaLaw` | `AreaLaw.lean` |
| Lemma `lem:slater`(a),(b) and `lem:quasifree`(b) for Slater determinants `Γ(e^{-itA})\|S⟩` and left segments (`S(ρ_A) = ∑_i h(ν_i)`, `ν_i ∈ [0,1]` the eigenvalues of `C_A = P_A K P_A`; `N_A` Poisson-binomial with generating function `det(1 + (z-1) C_A)`; `∑_i √(ν_i(1-ν_i)) ≤ ∑_{x∈A, y∉A} \|K_{xy}\|`) | `corrA_isHermitian`, `eig_corrA_mem`, `vonNeumann_rdm_slater`, `numDist_slater`, `genfun_slater`, `sum_sqrt_eig_le` | `SlaterExact.lean` |
| Lemma `lem:slater`(a), `lem:quasifree`(b) for **arbitrary regions** `R` (e.g. bulk intervals): `N_R` Poisson-binomial in the eigenvalues of `C_R`, generating function `det(1 + (z-1) C_R)` | `corrR_isHermitian`, `eig_corrR_mem`, `genfun_slater_gen`, `numDist_slater_gen` | `SlaterGeneral.lean` |
| Lemma `lem:slater`(a),(b) for **arbitrary regions** with the fermionic reduced state (`rdmF`, reordering sign of `𝓕 ≅ 𝓕_R ⊗ 𝓕_B`): `S(ρ_R) = ∑_i h(ν_i)`, `S(ρ_R) ≤ 2 ln 2 ∑_{x∈R, y∉R} \|K_{xy}\|`, agreement with `rdm` for left segments; bulk area law `E_q sup_t S(ρ^q_B(t)) ≤ 4 ln 2 Θ(F)` (finite-volume form of `thm:infvol`(d)) | `rdmF`, `reorderSign`, `vonNeumann_rdmF_slater`, `vonNeumann_rdmF_le`, `vonNeumann_rdmF_segment`, `lgt_areaLaw_bulk` | `SlaterBulk.lean` |
| Lemma `lem:binent` (`h(ν) ≤ 2 ln 2 √(ν(1-ν))`, equality iff `ν ∈ {0, ½, 1}`) | `binEntropy_le_two_log_two_sqrt`, `binEntropy_eq_two_log_two_sqrt_iff` | `BinEntropy.lean` |
| Lemma `lem:annealed`(i),(ii) (first bound) and Theorem `thm:arealaw`(a),(b) **with the constants of the paper**: `S(ρ_A) ≤ 2 ln 2 ∑_{x∈A, y∉A} \|K_{xy}\|`, `E_q sup_t S(ρ^q_A(t)) ≤ 2 ln 2 Θ(F) ≤ 8 ln 2 M_2(F)` | `annealed_K`, `annealed_cut`, `vonNeumann_slater_sharp`, `lgt_areaLaw_sharp`, `thetaF_le`, `lgt_areaLaw_sharp_of_ldt` | `AreaLawSharp.lean` |
| Lemma `lem:annealed`(ii) (remaining bounds) and (iii) (bulk intervals: constants `2Θ(F)`, `2Σ₁(F)`, `4Σ'₁(F)²`) | `annealed_cut_Q`, `annealed_cut_Q_sq`, `annealed_ivl_K`, `annealed_ivl_Q`, `annealed_ivl_Q_sq` | `AnnealedInterval.lean` |
| Number entropy of bulk intervals in finite volume (the finite-volume part of the proof of `thm:infvol`(e)): `⟨(N_B - s_B)²⟩ ≤ 2Θ(F) + 4Σ'₁(F)²`, `S_N(t) ≤ ½ ln(2πe(2Θ(F) + 4Σ'₁(F)² + 1/12))` | `lgt_bulk_transport_sq`, `lgt_numberEntropy_bulk` | `NumberEntropyBulk.lean` |
| Theorem `thm:arealaw`(c), Lemma `lem:cutlocal`(c) (definite charges) | `vonNeumann_rdmLGT_xprod`, `lgt_areaLaw_full` | `CutLocal.lean` |
| Theorem `thm:arealaw`(b),(c) **with the constants of the paper** (`E_ξ S_full(t) = S_G(t)`, `sup_t S_G ≤ 8 ln 2 M_2(F)`, `E_ξ sup_t S_full ≤ 8 ln 2 M_2(F)`, Markov bound `P_ξ(sup_t S_full ≥ s) ≤ 8 ln 2 M_2(F)/s`) | `avg_full_eq_SG`, `sup_SG_le`, `avg_sup_full_le`, `prob_sup_full_ge_le` | `AreaLawFull.lean` |
| Lemma `lem:pinching`; concavity and mixing bound of the von Neumann entropy | `vonNeumann_le_pinching`, `vonNeumann_pinching_le`, `vonNeumann_concave`, `vonNeumann_mix_le`, `vonNeumann_le_diag`, `vonNeumann_le_log_card`, `vonNeumann_mul_conjTranspose_comm` | `VNEntropy.lean` |
| Lemma `lem:cutlocal`(a) (in configuration form: `U⁻¹` of a split configuration is the split of the dual configurations after the controlled flip `C_ℓ` of the right links), (b),(d) (duality defect `≤ ln 2`, sharp) | `EntropiesAux.dualEquiv_symm_splitLGT`, `cutlocal_defect`, `exProd_entropy`, `exProd_U_entropy`, `exProd_Udag_entropy` | `Entropies.lean`, `CutLocalExamples.lean` |
| Proposition `prop:entropies`(a)–(c) (entropy taxonomy) | `fullEntropy_eq_numberEntropy_add`, `bareEntropy_eq_numberEntropy_add`, `vonNeumann_block_decomp`, `lgt_duality_defect`, `mixRdm_eq_marginal`, `condEntropy_le_mix`, `mix_le_cond_add`, `mix_le_log`, `sectorLawEntropy_le`, `sectorLawEntropy_zero`, `lgtNumDist_eq_avg`, `numberEntropy_le_mix`, `numberEntropy_le_bare`, `avg_sector_shannon_le_cond`, `avg_sector_shannon_le_numberEntropy`, `numberEntropy_le_avg_add`; `shannon_le_vonNeumann_of_block`, `numDist_shannon_le_vonNeumann_slater`, `lgt_numberEntropy_le_full` | `Entropies.lean`, `Taxonomy.lean`, `EntropyBlocks.lean` |
| Proposition `prop:entropies`(d) (superposition = ensemble, all multi-time correlators) | `expect_initZ_eq_trace_ensemble`, `correlator_eq_ensemble`, `gaugeInv_heis` | `Ensemble.lean` |
| Lemma `lem:duhamel` | `colNorm_evo_sub_le`, `avg_colNorm_evo_sub_le` | `Duhamel.lean` |
| Theorem `thm:EC` (Hypothesis (EC) from uniform large deviations), with Lemmas `lem:green`, `lem:entry`, `lem:badset`, `lem:double` | `ec_of_ldt`, `corrQ_decay_of_ldt`; `poisson`, `box_right`, `box_left`, `dirL_eq_dirR`, `sl2_entry_ineq`, `eigvec_recursion`, `abs_eigenvalue_hq_le`; `decay_of_not_res`, `exists_resEnd`, `card_resEnds`; `LDT`, `entry_ldt`, `single_box`, `two_box` | `Transfer.lean`, `Resonance.lean`, `LDTEntry.lean`, `ECProof.lean` |
| Theorem `thm:EC`(ii) (uniformity over the volume: `E sup_{Λ ∋ x,y} Q^Λ(x,y) ≤ Ce^{-η\|x-y\|}`, supremum over all windows, uniformly in `N`) | `ec_sup_of_ldt`, `corrQwin`, `corrQsup` | `ECSup.lean` |
| Main results with (EC) replaced by `LDT` | `lgt_transport_of_ldt`, `lgt_numberEntropy_of_ldt`, `lgt_areaLaw_of_ldt`, `lgtImbalance_lower_of_ldt` | `ECProof.lean` |
| Lemma `lem:propagation`(b), Corollary `cor:ECinfvol` (propagator form) and Theorem `thm:infvol`(b) (infinite volume: `U^ℤ_q(t)` as the exponential series of the tridiagonal `h^ℤ_q`, i.i.d. charges on `ℤ`; `E sup_t \|U^ℤ_q(t)_{xy}\| ≤ F(\|x-y\|)`, one-particle density, moments, transported charge) | `tendsto_evoZ`, `lintegral_iSup_UZ_le`, `infvol_propagator_of_ldt`, `tendsto_density_one`, `moments_UZ`, `abs_transportM_le`, `transport_tendsto` | `InfVol.lean`, `InfVolTransport.lean` |
| Theorem `thm:infvol`(e) (bulk intervals `A`: the law of `N_A` in `Ψ^{(M)}(t)` converges to `p_t = E PB(ν(K^ℤ_q(t)\|_A))`, `S_N^{(M)}(t) → H(p_t)`, and `H(p_t) ≤ ½ ln(2πe(2Θ(F) + 4Σ'₁(F)² + 1/12))`) | `KZ`, `KZA`, `pZ`, `tendsto_numDistM`, `tendsto_numberEntropyM`, `numberEntropyZ_le` | `InfVolNumber.lean` |
| Theorem `thm:infvol`(c) (Lieb–Robinson bound for volume-independent Majorana expansions on `X ⊂ ℤ` and observables localized in `Y`, uniformly in `M`, hence for the `limsup`) | `majListM`, `infvol_LR` | `InfVolLR.lean` |
| Proposition `prop:nocritical`(i),(ii),(iv) (hypotheses of Furstenberg's theorem; no critical energies) | `transferT_mul_inv`, `transferT_commutator`, `transferT_not_commute`, `unipotent_pow_unbounded`, `no_finite_invariant_lines` | `Furstenberg.lean` |
| Remark `rem:randomwalk` (one charge of the standard ℤ₂ LGT in a random-walk potential) | `Hstd_dual_oneP`, `partialSum_bijective` | `RandomWalk.lean` |
| Proposition `prop:dimerlaw` (dimer sector law; Definition `def:dimerstate`, for the CDW dimer state) | `dimerInit_dual`, `card_isDimer`, `sector_average_dimer` | `Dimer.lean` |
| Lemma `lem:dimercritical`, JSBS hypotheses | `dimerTM_plus_crit`, `dimerTM_minus_elliptic`, `dimerTM_commute_crit`, `jsbs_condition`, `dimerTM_minus_not_elliptic` | `DimerCritical.lean` |
| Lemma `lem:dimercritical`(i),(iii),(iv) and eq. `eq:dimerangles` (commutator `[D_+, D_-] = 2λ(λ²-E²)σ_x`; rotation normal form `S⁻¹D_+(λ)S = -1`, `S⁻¹D_-(λ)S = R_{2θ}`; averaged phases; `κ = ‖S‖‖S⁻¹‖ = ((1+λ)/(1-λ))^{1/2}`; `E_c = -λ`), Lemma `lem:dimertransfer` (`‖T^z(n,m)‖ ≤ K₁e^{K₂\|z-λ\|(n-m)}` for every dimer configuration) | `dimerD_commutator`, `dimerD_commute_iff`, `dimerTheta_spec`, `dimerS_conj_plus`, `dimerS_conj_tm`, `dimerS_conj_minus`, `dimer_angles`, `dimerS_kappa`, `gamma2_tmC`, `dimerS_conj_neg`, `dimer_transfer_bound` | `DimerTransfer.lean` |
| Lemma `lem:dimerTL` (finite speed uniformly in the potential; finite-volume difference, for finite `Λ'`; ballistic bound `∑_y \|y\|^p \|⟨δ_y, e^{-isH}δ_0⟩\|² ≤ (2\|s\|)^p`) | `dimerTL_i`, `dimerTL_i'`, `dimerTL_ii`, `dimerTL_iii` | `DimerTL.lean` |
| Proposition `prop:dimerlower` on finite intervals (Plancherel identity eq. `eq:plancherel`; `M_p(T;x) ≥ c_{λ,p} T^{p-1}` for every dimer configuration) | `dimerMom`, `dimerGreen`, `dimerMom_plancherel`, `dimer_lower_bound` | `DimerLower.lean` |
| Theorem `thm:dimer`(a) for every dimer sector and every volume, and (d) for `p > 1` (moments unbounded sector by sector; (EC) fails for the dimer law) | `dimer_sector_moment_lower`, `dimer_sup_moment_unbounded`, `not_EC_dimer`, `dimerAvg` | `DimerLGT.lean` |
| Proposition `prop:dimerimbalance`(i),(ii) | `dimer_stagWeight_identity`, `dimerImbalance_timeAvg` | `DimerImbalance.lean` |
| Theorem `thm:2dmap` for `σᶻ`-polarized gauge states, on **any finite graph** | `FinGraph.sector_average_matter`, `FinGraph.gaussOp_comm_ham`, `FinGraph.dGamma_hq` | `Graph.lean` |
| Lemma `lem:2dgauge`(i)–(iii), Theorem `thm:2dmap`(a) (fluxes = gauge classes, `2^β` flux sectors, orbits of size `2^{\|V\|-1}`; Wilson loops; `∏ G_i = (-1)^{N_f}`) | `FinGraph.flux_eq_iff_gauge`, `FinGraph.gaugeAct_eq_self_iff`, `FinGraph.card_gaugeOrbit`, `FinGraph.card_fluxClasses`, `FinGraph.wilsonLoop_comm_ham`, `FinGraph.wilsonLoop_comm_gaussOp`, `FinGraph.gaussOp_prod` | `GraphFlux.lean`, `GraphWilson.lean` |
| Theorem `thm:2dmap`(b) (joint eigenbasis `Ψ_{n,q,s}` of the Gauss operators and Wilson loops, gauge covariance, completeness, sector Hamiltonians `dΓ(h_{q,s}) - h ∑ q_i`) | `FinGraph.fluxVec_gauss`, `FinGraph.fluxVec_wilson`, `FinGraph.fluxVec_occ`, `FinGraph.fluxVec_gaugeAct`, `FinGraph.fluxVec_norm`, `FinGraph.fluxVec_orthogonal`, `FinGraph.card_fluxLabels`, `FinGraph.ham_fluxVec` | `TwoDMap.lean` |
| Theorem `thm:2dmap`(c),(d) (planar boxes and tori: cycle bases from the plaquettes, resp. the plaquettes but one and the two non-contractible loops; `β = (L₁-1)(L₂-1)`, resp. `L₁L₂+1`; `∏_p W_p = 1` on the torus; all `2^β` flux labels occur) | `FinGraph.IsCycleBasis`, `FinGraph.gauge_iff_flux_basis`, `FinGraph.exists_flux_basis`, `FinGraph.card_cycleBasis`, `FinGraph.box_connected`, `FinGraph.box_beta`, `FinGraph.box_plaq_basis`, `FinGraph.box_card_fluxClasses`, `FinGraph.torus_connected`, `FinGraph.torus_beta`, `FinGraph.torus_sum_plaq`, `FinGraph.torus_prod_wilson`, `FinGraph.torus_basis`, `FinGraph.torus_card_fluxClasses` | `TwoDLattices.lean` |
| Proposition `prop:2dlaws`(a),(b) (any flux sector; `σˣ`-polarized links: random `π`-flux) | `FinGraph.sector_average_flux`, `FinGraph.sector_average_xpol`, `FinGraph.xpol_summand_gauge` | `GraphFlux.lean` |
| Proposition `prop:2dlaws`(c)(i),(ii) (checkerboard: `h² = K_s² + 4h²`, chiral gap), gauge covariance of `h_{q,s}` | `FinGraph.checkerboard_sq`, `FinGraph.checkerboard_gap`, `FinGraph.checkerboard_edge_iff`, `FinGraph.Ks_anticomm`, `FinGraph.gauge_hqs`, `FinGraph.hqs_one` | `TwoD.lean` |
| Proposition `prop:2dlaws`(c)(iii),(iv) (checkerboard propagator `e^{-ith} = cos(tR) - i h R⁻¹ sin(tR)`, `R² = K_s² + 4h²`, on eigenvectors and as power series; `(K_s²)_{ii} = J² deg(i)` on simple graphs) | `FinGraph.checkerboard_evo_eigvec`, `FinGraph.checkerboard_evo`, `FinGraph.Ks_sq_diag` | `TwoDLaws.lean` |
| Proposition `prop:nocompact` (no compactly supported eigenvectors on `ℤ²`) | `no_compact_eigenvector`, `no_compact_eigenvector_hqs` | `TwoD.lean` |
| Proposition `prop:cluster` (strong-coupling cluster confinement; (c) with the percolation probability as input) | `FinGraph.hqs0_cluster`, `FinGraph.cluster_bound`, `FinGraph.cluster_leak`, `FinGraph.cluster_avg` | `Cluster.lean` |
| Proposition `prop:gausscap` (U(1) quantum links: Gauss-law caps, superpositions of sectors, any dynamics) | `countInt_gauss`, `qlm_numberEntropy_le`, `qlm_numberEntropy_le_edge`, `qlm_numberEntropy_le_super`, `qlm_numberEntropy_le_super_edge`, `qlm_sector_preserved`, `bgEntropy_invariant` | `QLM.lean`, `QLMSuper.lean` |
| Proposition `prop:wilson` (Wilson links: cut-averaged number entropy from energy conservation) | `wField_hop`, `wField_eq_count`, `wilson_field_bound`, `wilson_numberEntropy_avg`, `wilson_numberEntropy_most`, `wilson_numberEntropy_each` | `Wilson.lean` |

## Verification

`lake build` compiles all 68 modules (about 35,200 lines) without `sorry`; there is no `axiom`,
`admit`, `native_decide` or `implemented_by` in the library. The script `CheckAxioms.lean` walks
through every theorem declared in a `DFL` module (including private helper lemmas) and collects the
axioms it depends on:

```
lake env lean CheckAxioms.lean
# checked 3828 theorems in DFL modules; with non-standard axioms: 0
```

i.e. every theorem depends only on Lean's standard axioms `propext`, `Classical.choice` and
`Quot.sound`. Individual theorems can be checked with `#print axioms`, e.g.

```
cat > Check.lean <<'X'
import DFL
#print axioms DFL.ec_of_ldt
#print axioms DFL.lgt_transport_of_ldt
#print axioms DFL.lgt_numberEntropy_of_ldt
#print axioms DFL.lgtImbalance_lower_of_ldt
#print axioms DFL.cutlocal_defect
#print axioms DFL.FinGraph.sector_average_xpol
#print axioms DFL.FinGraph.cluster_bound
X
lake env lean Check.lean
```

## Conventions and deviations from the paper

- **Indexing.** Sites and links are 0-indexed; the boundary link is `b_{N-1}`. The density wave
  `cdw N` occupies the even 0-indexed sites. The dimer state of `Dimer.lean` uses that density wave
  and fixes the odd links; it is a different but equivalent instance of the paper's
  Lemma `lem:productlaw`(e). On the ring the extra bond is `(N-1, 0)`, carried by `b_{N-1}`. In
  `Wilson.lean` Lean site `j` is paper site `j + 1`.
- **Charge labels.** Gauss charges, link signs and Wilson-loop eigenvalues are stored in `ZMod 2`
  and converted to `±1` by `sgn`.
- **Two proofs of the area law.** `AreaLaw.lean` proves the area law with the constant `2` without
  computing the spectrum of `ρ_A`: `S(ρ_A) ≤ 2 log ∑_i √λ_i` by Jensen, `∑_i √λ_i ≤ ∏_k (√ν_k +
  √(1-ν_k))` by a product decomposition of the Slater determinant in canonical orbitals, and
  `log(√ν + √(1-ν)) ≤ √(ν(1-ν))`. `SlaterExact.lean`, `BinEntropy.lean` and `AreaLawSharp.lean` follow the
  paper: the Schmidt decomposition across the cut gives `S(ρ_A) = ∑_i h(ν_i)` exactly, and with
  `lem:binent` the constant of the paper, `2 ln 2 Θ(F)`. The reduced state is identified through its
  entropy and the law of `N_A` (the quasi-free density matrix `ρ_C` of `lem:quasifree`(a) is not
  defined in Lean).
- **Constants.** `TransportSharp.lean` and `AreaLawSharp.lean` state transport, number entropy and the
  area law with the explicit constants of the paper (`Σ₁(F)`, `Θ(F) + Σ'₁(F)²`, `2 ln 2 Θ(F)`). The
  earlier versions in `Main.lean` and `AreaLaw.lean` have weaker or existential constants (transport
  `2Σ₁(F)`, number entropy and area law `∃ K`). The constants `C, η` of (EC) are existential, as in
  the paper.
- **Proof routes that differ from the paper.** In `cluster_bound` (`‖e^{-ith} - e^{-ith⁰}‖ ≤ 8J/g +
  16|t|J²/g`, the propagator argument of the paper via `e^{iuh}(1-S)e^{-iuh⁰}`) the generator `S` with
  `[h⁰, S] = h¹`, `‖S‖ ≤ 4J/g` is obtained from the anticommutator equation for `Γ h⁰ ≥ g/2` instead of
  the integral formula. (The paper's constant was `48` before this argument was found during the
  formalization; it is now `16` in both.) The
  Green-function formula of `lem:green` is replaced by the Poisson identity
  `P_B(E) φ(x) = p_a(E) φ(a-1) + p_b(E) φ(b+1)` (Wronskians), and `lem:entry` uses only the vector
  form of the large-deviation estimate (the norm of the product is bounded through two vectors).

## Not formalized

- The literature theorems used as inputs: the uniform large-deviation estimate and the positivity
  and continuity of the Lyapunov exponent (Furstenberg; `prop:LDT`, entering as the hypothesis
  `LDT`), the verification of the hypotheses of the cited theorems in their own formalism
  (`lem:hypcheck`), the Zariski density of `prop:nocritical`(iii), Ge–Zhao's annealed localization
  on `ℤ` (`prop:Zline`), the transport bounds of Jitomirskaya–Schulz-Baldes–Stolz for the random
  dimer (`thm:dimer`(b),(c), and (d) for `1/2 < p ≤ 1`), and the percolation estimate of Li–Liang
  used in `prop:cluster`(c) (entering `cluster_avg` as the probability that the cluster leaves the
  region).
- Infinite-volume statements beyond the one-body propagator, the Lieb–Robinson bound and the number entropy:
  `thm:infvol`(a),(d) (the quasi-free limit state and the area-law limits),
  `cor:ECinfvol` for the full correlator `Q^X` (Borel functional calculus on `ℓ²(X)`; the
  propagator bound is formalized), part (iii) of `prop:dimerimbalance` and the density of states of
  the dimer model (`lem:dimerdos`), which need spectral measures on `ℓ²(ℤ)`; `prop:LR`(d) (no
  operator-norm zero-velocity bound, which uses the ballistic motion of the free chain on `ℤ`).
- Remarks, conjectures and numerics.
