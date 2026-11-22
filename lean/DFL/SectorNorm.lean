import DFL.GaussLaw
import DFL.OneBody

/-!
# Sector-diagonal operators (paper: Lemma `lem:sectornorm`)

Let `X` be a gauge-invariant operator on the LGT Hilbert space and `P_q` the sector projections of
`DFL.GaussLaw`.

* (a) `X = ∑_q P_q X P_q` and `‖X‖ = max_q ‖P_q X P_q‖` (`gaugeInv_eq_sum`, `norm_gaugeInv`).
* (b) For a density matrix `ρ` with sector weights `w_q = Tr(ρ P_q)`, `|Tr(ρ X)| ≤ ∑_q w_q ‖P_q X P_q‖`
  (`norm_trace_le`); `w_q = 2^{-N}` for `ρ = 4^{-N} 1` and for `ρ = |Ψ₀⟩⟨Ψ₀|`,
  `Ψ₀ = |↑⋯↑⟩_z ⊗ ψ` (`trace_mixed_gaussProj`, `trace_initZ_gaussProj`), hence
  `|Tr(ρ X)| ≤ E_q ‖P_q X P_q‖` (`norm_trace_mixed_le`, `norm_expect_initZ_le`).
* (c) If `X = ∑_q x_q P_q`, then `|X| = ∑_q |x_q| P_q` (the positive semidefinite square root of
  `X^* X`; `sectorDiag_abs`) and `⟨Ψ₀, |X| Ψ₀⟩ = 4^{-N} Tr |X| = E_q |x_q|`
  (`expect_sectorDiag_abs`, `trace_sectorDiag_abs`).

Norms of matrices are `ℓ²` operator norms (`Matrix.Norms.L2Operator`).
-/

open Matrix Finset
open scoped Kronecker Matrix.Norms.L2Operator ComplexOrder

noncomputable section

namespace DFL

variable {N : ℕ}

namespace SectorNormAux

section Analytic

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Squared Euclidean norm of a matrix-vector product, in coordinates. -/
lemma norm_toEuclideanCLM_sq (A : Matrix ι ι ℂ) (w : EuclideanSpace ℂ ι) :
    ‖toEuclideanCLM (n := ι) (𝕜 := ℂ) A w‖ ^ 2 = ∑ a, ‖(A *ᵥ WithLp.ofLp w) a‖ ^ 2 := by
  rw [EuclideanSpace.norm_sq_eq]
  simp

/-- An operator-norm bound from a bound on the coordinates of `A v`. -/
lemma norm_le_of_sq (A : Matrix ι ι ℂ) (M : ℝ) (hM : 0 ≤ M)
    (h : ∀ v : ι → ℂ, ∑ a, ‖(A *ᵥ v) a‖ ^ 2 ≤ M ^ 2 * ∑ a, ‖v a‖ ^ 2) : ‖A‖ ≤ M := by
  rw [Matrix.cstar_norm_def]
  refine ContinuousLinearMap.opNorm_le_bound _ hM fun w => ?_
  have h1 := h (WithLp.ofLp w)
  rw [← norm_toEuclideanCLM_sq] at h1
  have h2 : ∑ a, ‖WithLp.ofLp w a‖ ^ 2 = ‖w‖ ^ 2 := by
    rw [EuclideanSpace.norm_sq_eq]
  rw [h2, ← mul_pow] at h1
  exact (sq_le_sq₀ (norm_nonneg _) (by positivity)).mp h1

/-- Coordinate form of `‖A v‖ ≤ ‖A‖ ‖v‖`. -/
lemma sum_sq_le (A : Matrix ι ι ℂ) (v : ι → ℂ) :
    ∑ a, ‖(A *ᵥ v) a‖ ^ 2 ≤ ‖A‖ ^ 2 * ∑ a, ‖v a‖ ^ 2 := by
  have h := (toEuclideanCLM (n := ι) (𝕜 := ℂ) A).le_opNorm (WithLp.toLp 2 v)
  rw [← Matrix.cstar_norm_def] at h
  have h2 := pow_le_pow_left₀ (norm_nonneg _) h 2
  rw [mul_pow, norm_toEuclideanCLM_sq, EuclideanSpace.norm_sq_eq] at h2
  simpa using h2

/-- `|⟨v, Y v⟩| ≤ ‖Y‖ ‖v‖²`. -/
lemma norm_dotProduct_mulVec_le (Y : Matrix ι ι ℂ) (v : ι → ℂ) :
    ‖star v ⬝ᵥ (Y *ᵥ v)‖ ≤ ‖Y‖ * (star v ⬝ᵥ v).re := by
  set w : EuclideanSpace ℂ ι := WithLp.toLp 2 v with hw
  set T := toEuclideanCLM (n := ι) (𝕜 := ℂ) Y with hT
  have hinner : star v ⬝ᵥ (Y *ᵥ v) = inner ℂ w (T w) := by
    rw [EuclideanSpace.inner_eq_star_dotProduct, hT, ofLp_toEuclideanCLM, dotProduct_comm]
  have hnorm : (star v ⬝ᵥ v).re = ‖w‖ ^ 2 := by
    rw [EuclideanSpace.norm_sq_eq]
    simp only [dotProduct, Complex.re_sum, hw]
    refine Finset.sum_congr rfl fun a _ => ?_
    rw [Pi.star_apply, Complex.star_def, Complex.conj_mul']
    norm_cast
  rw [hinner, hnorm]
  calc ‖inner ℂ w (T w)‖ ≤ ‖w‖ * ‖T w‖ := norm_inner_le_norm _ _
    _ ≤ ‖w‖ * (‖T‖ * ‖w‖) := mul_le_mul_of_nonneg_left (T.le_opNorm w) (norm_nonneg _)
    _ = ‖Y‖ * ‖w‖ ^ 2 := by rw [hT, ← Matrix.cstar_norm_def]; ring

omit [DecidableEq ι] in
/-- `Tr(|u⟩⟨u| Y) = ⟨u, Y u⟩`. -/
lemma trace_vecMulVec_mul (u : ι → ℂ) (Y : Matrix ι ι ℂ) :
    (vecMulVec u (star u) * Y).trace = star u ⬝ᵥ (Y *ᵥ u) := by
  rw [vecMulVec_mul, trace_vecMulVec, dotProduct_comm, ← dotProduct_mulVec]

/-- `|Tr(σ Y)| ≤ Tr(σ) ‖Y‖` for positive semidefinite `σ`. -/
lemma norm_trace_mul_le (σ Y : Matrix ι ι ℂ) (hσ : σ.PosSemidef) :
    ‖(σ * Y).trace‖ ≤ σ.trace.re * ‖Y‖ := by
  obtain ⟨m, v, rfl⟩ := Matrix.posSemidef_iff_eq_sum_vecMulVec.mp hσ
  rw [Finset.sum_mul, trace_sum, trace_sum, Complex.re_sum, Finset.sum_mul]
  refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun i _ => ?_)
  rw [trace_vecMulVec_mul, trace_vecMulVec, dotProduct_comm (v i), mul_comm]
  exact norm_dotProduct_mulVec_le Y (v i)

end Analytic

/-- The complexified sector projection is the diagonal indicator of the charge `q`. -/
lemma cplx_gaussProj (q : Cfg N) :
    cplx (gaussProj q) = diagonal (fun p => if chargeOf p.1 p.2 = q then (1 : ℂ) else 0) := by
  rw [gaussProj, sec_cplx_diagonal]
  congr 1
  funext p
  split_ifs <;> simp

/-- `P_q² = P_q`. -/
lemma cplx_gaussProj_mul_self (q : Cfg N) :
    cplx (gaussProj q) * cplx (gaussProj q) = cplx (gaussProj q) := by
  rw [cplx_gaussProj, diagonal_mul_diagonal]
  congr 1
  funext p
  split_ifs <;> simp

/-- `P_q` is self-adjoint. -/
lemma cplx_gaussProj_conjTranspose (q : Cfg N) :
    (cplx (gaussProj q))ᴴ = cplx (gaussProj q) := by
  rw [cplx_gaussProj, diagonal_conjTranspose]
  congr 1
  funext p
  simp only [Pi.star_apply]
  split_ifs <;> simp

/-- Entries of `P_q X P_q`. -/
lemma block_apply (X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (q : Cfg N)
    (a b : Cfg N × Cfg N) :
    (cplx (gaussProj q) * X * cplx (gaussProj q)) a b =
      if chargeOf a.1 a.2 = q ∧ chargeOf b.1 b.2 = q then X a b else 0 := by
  rw [cplx_gaussProj, mul_diagonal, diagonal_mul]
  by_cases ha : chargeOf a.1 a.2 = q <;> by_cases hb : chargeOf b.1 b.2 = q <;> simp [ha, hb]

/-- Gauge-invariant operators have no matrix elements between different charges. -/
lemma gaugeInv_apply_eq_zero (X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hX : GaugeInv X)
    (a b : Cfg N × Cfg N) (h : chargeOf a.1 a.2 ≠ chargeOf b.1 b.2) : X a b = 0 := by
  have e := congrFun (congrFun (gaugeInv_submatrix X hX) (dualEquiv.symm a)) (dualEquiv.symm b)
  rw [submatrix_apply, Equiv.apply_symm_apply, Equiv.apply_symm_apply, blockDiagonal_apply] at e
  rw [e, ite_eq_right]
  simpa using h

/-- For gauge-invariant `X`, `(P_q X P_q) v = P_q (X v)`. -/
lemma block_mulVec (X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hX : GaugeInv X) (q : Cfg N)
    (v : Cfg N × Cfg N → ℂ) (a : Cfg N × Cfg N) :
    ((cplx (gaussProj q) * X * cplx (gaussProj q)) *ᵥ v) a =
      if chargeOf a.1 a.2 = q then (X *ᵥ v) a else 0 := by
  simp only [mulVec, dotProduct, block_apply]
  split_ifs with ha
  · refine Finset.sum_congr rfl fun b _ => ?_
    by_cases hb : chargeOf b.1 b.2 = q
    · simp [ha, hb]
    · rw [ite_eq_right (fun h => hb h.2),
        gaugeInv_apply_eq_zero X hX a b (fun h => hb (h.symm.trans ha))]
  · exact Finset.sum_eq_zero fun b _ => by rw [ite_eq_right (fun h => ha h.1), zero_mul]

/-- `P_q X P_q` only sees the component of the vector in the sector `q`. -/
lemma block_mulVec_proj (X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (q : Cfg N)
    (v : Cfg N × Cfg N → ℂ) :
    (cplx (gaussProj q) * X * cplx (gaussProj q)) *ᵥ v =
      (cplx (gaussProj q) * X * cplx (gaussProj q)) *ᵥ
        (fun b => if chargeOf b.1 b.2 = q then v b else 0) := by
  ext a
  simp only [mulVec, dotProduct, block_apply]
  refine Finset.sum_congr rfl fun b _ => ?_
  by_cases hb : chargeOf b.1 b.2 = q <;> simp [hb]

/-- `chargeOf n x = q ↔ x = linkOf n q`. -/
lemma chargeOf_eq_iff (n x q : Cfg N) : chargeOf n x = q ↔ x = linkOf n q := by
  constructor
  · rintro rfl
    exact (linkOf_chargeOf n x).symm
  · rintro rfl
    exact chargeOf_linkOf n q

/-- For fixed occupations, `x ↦ chargeOf n x` is a bijection. -/
def chargeEquiv (n : Cfg N) : Cfg N ≃ Cfg N where
  toFun := chargeOf n
  invFun := linkOf n
  left_inv x := linkOf_chargeOf n x
  right_inv q := chargeOf_linkOf n q

/-- Summing a function of the charge over the link labels is summing over the charges. -/
lemma sum_comp_chargeOf {M : Type*} [AddCommMonoid M] (n : Cfg N) (f : Cfg N → M) :
    ∑ x, f (chargeOf n x) = ∑ q, f q :=
  Equiv.sum_comp (chargeEquiv n) f

/-- For fixed occupations, exactly one link configuration has charge `q`. -/
lemma sum_ite_chargeOf {M : Type*} [AddCommMonoid M] (n q : Cfg N) (g : Cfg N → M) :
    ∑ x, (if chargeOf n x = q then g x else 0) = g (linkOf n q) := by
  simp only [chargeOf_eq_iff]
  simp

/-- The number of configurations is `2^N`. -/
lemma card_cfg : Fintype.card (Cfg N) = 2 ^ N := by
  simp

/-- `Tr P_q = 2^N`. -/
lemma trace_cplx_gaussProj (q : Cfg N) : (cplx (gaussProj q)).trace = 2 ^ N := by
  rw [cplx_gaussProj, trace_diagonal, Fintype.sum_prod_type]
  simp only [sum_ite_chargeOf]
  simp

/-- `Re (2^{-N}) = 2^{-N}`. -/
lemma re_inv_two_pow : (((2 : ℂ) ^ N)⁻¹).re = 1 / 2 ^ N := by
  rw [show ((2 : ℂ) ^ N)⁻¹ = (((2 : ℝ) ^ N)⁻¹ : ℝ) by push_cast; rfl, Complex.ofReal_re, one_div]

/-- `4^N = 2^N 2^N`. -/
lemma four_pow : (4 : ℂ) ^ N = 2 ^ N * 2 ^ N := by
  rw [← mul_pow]
  norm_num

/-- `Tr(P_q ρ P_q) = Tr(ρ P_q)`. -/
lemma trace_proj_mul_proj (ρ : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (q : Cfg N) :
    (cplx (gaussProj q) * ρ * cplx (gaussProj q)).trace = (ρ * cplx (gaussProj q)).trace := by
  rw [mul_assoc, trace_mul_comm, mul_assoc, cplx_gaussProj_mul_self]

/-- `Tr(ρ P_q X P_q) = Tr((P_q ρ P_q)(P_q X P_q))`. -/
lemma trace_mul_block (ρ X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (q : Cfg N) :
    (ρ * (cplx (gaussProj q) * X * cplx (gaussProj q))).trace =
      (cplx (gaussProj q) * ρ * cplx (gaussProj q) *
        (cplx (gaussProj q) * X * cplx (gaussProj q))).trace := by
  have hPP := cplx_gaussProj_mul_self q
  set P := cplx (gaussProj q)
  have e1 : P * ρ * P * (P * X * P) = P * (ρ * (P * X * P)) := by
    simp only [mul_assoc]
    rw [← mul_assoc P P (X * P), hPP]
  have e2 : ρ * (P * X * P) * P = ρ * (P * X * P) := by
    simp only [mul_assoc]
    rw [hPP]
  rw [e1, trace_mul_comm P, e2]

/-- `expect X v = Tr(|v⟩⟨v| X)`. -/
lemma expect_eq_trace (X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (v : Cfg N × Cfg N → ℂ) :
    expect X v = (vecMulVec v (star v) * X).trace := by
  rw [trace_vecMulVec_mul]
  rfl

/-- A linear combination of the sector projections is diagonal. -/
lemma sum_smul_gaussProj (c : Cfg N → ℂ) :
    ∑ q : Cfg N, c q • cplx (gaussProj q) = diagonal (fun p => c (chargeOf p.1 p.2)) := by
  ext a b
  simp only [Matrix.sum_apply, Matrix.smul_apply, cplx_gaussProj, diagonal_apply, smul_eq_mul]
  split_ifs with hab
  · simp
  · simp

end SectorNormAux

open SectorNormAux

/-! ### (a) -/

/-- A gauge-invariant operator is block diagonal: `X = ∑_q P_q X P_q`. -/
theorem gaugeInv_eq_sum (X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hX : GaugeInv X) :
    X = ∑ q : Cfg N, cplx (gaussProj q) * X * cplx (gaussProj q) := by
  ext a b
  rw [Matrix.sum_apply]
  simp only [block_apply]
  by_cases hab : chargeOf a.1 a.2 = chargeOf b.1 b.2
  · rw [Finset.sum_eq_single (chargeOf a.1 a.2)]
    · simp [hab]
    · intro q _ hq
      exact ite_eq_right fun h => hq h.1.symm
    · simp
  · rw [gaugeInv_apply_eq_zero X hX a b hab]
    simp

/-- **Norm of a gauge-invariant operator** (paper: Lemma `lem:sectornorm`(a)):
`‖X‖ = max_q ‖P_q X P_q‖`. -/
theorem norm_gaugeInv (X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hX : GaugeInv X) :
    ‖X‖ = univ.sup' univ_nonempty
      (fun q : Cfg N => ‖cplx (gaussProj q) * X * cplx (gaussProj q)‖) := by
  have hle : ∀ q : Cfg N, ‖cplx (gaussProj q) * X * cplx (gaussProj q)‖ ≤
      univ.sup' univ_nonempty
        (fun q : Cfg N => ‖cplx (gaussProj q) * X * cplx (gaussProj q)‖) := fun q =>
    Finset.le_sup' (fun q : Cfg N => ‖cplx (gaussProj q) * X * cplx (gaussProj q)‖) (mem_univ q)
  have hM0 : 0 ≤ univ.sup' univ_nonempty
      (fun q : Cfg N => ‖cplx (gaussProj q) * X * cplx (gaussProj q)‖) :=
    (norm_nonneg _).trans (hle 0)
  apply le_antisymm
  · refine norm_le_of_sq X _ hM0 fun v => ?_
    have hX_sum : ∀ a, ‖(X *ᵥ v) a‖ ^ 2 =
        ∑ q : Cfg N, ‖((cplx (gaussProj q) * X * cplx (gaussProj q)) *ᵥ v) a‖ ^ 2 := by
      intro a
      simp only [block_mulVec X hX]
      simp [apply_ite]
    have hv_sum : ∀ a, ‖v a‖ ^ 2 =
        ∑ q : Cfg N, ‖(fun b => if chargeOf b.1 b.2 = q then v b else 0) a‖ ^ 2 := by
      intro a
      simp [apply_ite]
    calc ∑ a, ‖(X *ᵥ v) a‖ ^ 2
        = ∑ q : Cfg N, ∑ a, ‖((cplx (gaussProj q) * X * cplx (gaussProj q)) *ᵥ
            (fun b => if chargeOf b.1 b.2 = q then v b else 0)) a‖ ^ 2 := by
          simp_rw [hX_sum, ← block_mulVec_proj]
          exact Finset.sum_comm
      _ ≤ ∑ q : Cfg N, ‖cplx (gaussProj q) * X * cplx (gaussProj q)‖ ^ 2 *
            ∑ a, ‖(fun b => if chargeOf b.1 b.2 = q then v b else 0) a‖ ^ 2 :=
          Finset.sum_le_sum fun q _ => sum_sq_le _ _
      _ ≤ ∑ q : Cfg N, (univ.sup' univ_nonempty
            (fun q : Cfg N => ‖cplx (gaussProj q) * X * cplx (gaussProj q)‖)) ^ 2 *
            ∑ a, ‖(fun b => if chargeOf b.1 b.2 = q then v b else 0) a‖ ^ 2 := by
          refine Finset.sum_le_sum fun q _ => mul_le_mul_of_nonneg_right ?_
            (Finset.sum_nonneg fun _ _ => sq_nonneg _)
          exact pow_le_pow_left₀ (norm_nonneg _) (hle q) 2
      _ = (univ.sup' univ_nonempty
            (fun q : Cfg N => ‖cplx (gaussProj q) * X * cplx (gaussProj q)‖)) ^ 2 *
            ∑ a, ‖v a‖ ^ 2 := by
          rw [← Finset.mul_sum, Finset.sum_comm]
          simp_rw [← hv_sum]
  · refine Finset.sup'_le _ _ fun q _ => norm_le_of_sq _ _ (norm_nonneg _) fun v => ?_
    calc ∑ a, ‖((cplx (gaussProj q) * X * cplx (gaussProj q)) *ᵥ v) a‖ ^ 2
        ≤ ∑ a, ‖(X *ᵥ v) a‖ ^ 2 := by
          refine Finset.sum_le_sum fun a _ => ?_
          rw [block_mulVec X hX]
          split_ifs <;> simp
      _ ≤ ‖X‖ ^ 2 * ∑ a, ‖v a‖ ^ 2 := sum_sq_le _ _

/-! ### (b) -/

/-- **Expectations of gauge-invariant operators** (paper: Lemma `lem:sectornorm`(b)):
`|Tr(ρ X)| ≤ ∑_q Tr(ρ P_q) ‖P_q X P_q‖` for every positive semidefinite `ρ`. -/
theorem norm_trace_le (ρ X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hρ : ρ.PosSemidef)
    (hX : GaugeInv X) :
    ‖(ρ * X).trace‖ ≤
      ∑ q : Cfg N, (ρ * cplx (gaussProj q)).trace.re * ‖cplx (gaussProj q) * X * cplx (gaussProj q)‖ := by
  conv_lhs => rw [gaugeInv_eq_sum X hX]
  rw [Finset.mul_sum, trace_sum]
  simp_rw [trace_mul_block ρ X]
  refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun q _ => ?_)
  have hσ : (cplx (gaussProj q) * ρ * cplx (gaussProj q)).PosSemidef := by
    have := hρ.conjTranspose_mul_mul_same (cplx (gaussProj q))
    rwa [cplx_gaussProj_conjTranspose] at this
  refine (norm_trace_mul_le _ _ hσ).trans (le_of_eq ?_)
  rw [trace_proj_mul_proj]

/-- The sector weights of the maximally mixed state `4^{-N} 1` are `2^{-N}`. -/
theorem trace_mixed_gaussProj (q : Cfg N) :
    (((4 : ℂ) ^ N)⁻¹ • (1 : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) * cplx (gaussProj q)).trace =
      ((2 : ℂ) ^ N)⁻¹ := by
  rw [smul_mul_assoc, one_mul, trace_smul, trace_cplx_gaussProj, smul_eq_mul, four_pow]
  field_simp

/-- The sector weights of `Ψ₀ = |↑⋯↑⟩_z ⊗ ψ` (`ψ` a unit vector) are `2^{-N}`. -/
theorem trace_initZ_gaussProj (ψ : Cfg N → ℂ) (hψ : star ψ ⬝ᵥ ψ = 1) (q : Cfg N) :
    (vecMulVec (initZ ψ) (star (initZ ψ)) * cplx (gaussProj q)).trace = ((2 : ℂ) ^ N)⁻¹ := by
  rw [trace_vecMulVec_mul, cplx_gaussProj]
  simp only [dotProduct, mulVec_diagonal, Fintype.sum_prod_type, initZ, Pi.star_apply]
  have e : ∀ n : Cfg N, ∑ x : Cfg N, star ((((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ * ψ n) *
      ((if chargeOf n x = q then (1 : ℂ) else 0) * ((((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ * ψ n)) =
      star ((((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹) * (((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ *
        (star (ψ n) * ψ n) := by
    intro n
    simp only [ite_mul, one_mul, zero_mul, mul_ite, mul_zero]
    rw [sum_ite_chargeOf n q (fun _ => star ((((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ * ψ n) *
      ((((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ * ψ n))]
    rw [star_mul']
    ring
  simp only [e]
  rw [← Finset.mul_sum, sec_norm_sq_prefactor]
  have hψ' : ∑ n, star (ψ n) * ψ n = 1 := hψ
  rw [hψ', mul_one, one_div]

/-- `|Tr(4^{-N} X)| ≤ E_q ‖P_q X P_q‖` for gauge-invariant `X` (paper: Lemma `lem:sectornorm`(b)). -/
theorem norm_trace_mixed_le (X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hX : GaugeInv X) :
    ‖((4 : ℂ) ^ N)⁻¹ * X.trace‖ ≤
      avgQ (fun q : Cfg N => ‖cplx (gaussProj q) * X * cplx (gaussProj q)‖) := by
  have hρ : (((4 : ℂ) ^ N)⁻¹ • (1 : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ)).PosSemidef := by
    rw [smul_one_eq_diagonal]
    refine PosSemidef.diagonal fun _ => ?_
    rw [show ((4 : ℂ) ^ N)⁻¹ = (((4 : ℝ) ^ N)⁻¹ : ℝ) by push_cast; rfl]
    exact Complex.zero_le_real.mpr (by positivity)
  have h := norm_trace_le _ X hρ hX
  simp only [trace_mixed_gaussProj, re_inv_two_pow] at h
  rw [smul_mul_assoc, one_mul, trace_smul, smul_eq_mul] at h
  refine h.trans (le_of_eq ?_)
  rw [avgQ, Finset.mul_sum]

/-- `|⟨Ψ₀, X Ψ₀⟩| ≤ E_q ‖P_q X P_q‖` for gauge-invariant `X` and `Ψ₀ = |↑⋯↑⟩_z ⊗ ψ`
(paper: Lemma `lem:sectornorm`(b)). -/
theorem norm_expect_initZ_le (ψ : Cfg N → ℂ) (hψ : star ψ ⬝ᵥ ψ = 1)
    (X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hX : GaugeInv X) :
    ‖expect X (initZ ψ)‖ ≤
      avgQ (fun q : Cfg N => ‖cplx (gaussProj q) * X * cplx (gaussProj q)‖) := by
  rw [expect_eq_trace]
  refine (norm_trace_le _ X (posSemidef_vecMulVec_self_star _) hX).trans (le_of_eq ?_)
  simp only [trace_initZ_gaussProj ψ hψ, re_inv_two_pow]
  rw [avgQ, Finset.mul_sum]

/-! ### (c) -/

/-- **Absolute value of a sector-diagonal operator** (paper: Lemma `lem:sectornorm`(c)): for
`X = ∑_q x_q P_q`, the operator `∑_q |x_q| P_q` is positive semidefinite and squares to `X^* X`; it is
therefore `|X|`. -/
theorem sectorDiag_abs (x : Cfg N → ℂ) :
    (∑ q : Cfg N, ((‖x q‖ : ℝ) : ℂ) • cplx (gaussProj q)).PosSemidef ∧
      (∑ q : Cfg N, ((‖x q‖ : ℝ) : ℂ) • cplx (gaussProj q)) *
          (∑ q : Cfg N, ((‖x q‖ : ℝ) : ℂ) • cplx (gaussProj q)) =
        (∑ q : Cfg N, x q • cplx (gaussProj q))ᴴ * (∑ q : Cfg N, x q • cplx (gaussProj q)) := by
  rw [sum_smul_gaussProj, sum_smul_gaussProj]
  refine ⟨PosSemidef.diagonal fun p => Complex.zero_le_real.mpr (norm_nonneg _), ?_⟩
  rw [diagonal_conjTranspose, diagonal_mul_diagonal, diagonal_mul_diagonal]
  congr 1
  funext p
  simp only [Pi.star_apply, Complex.star_def, Complex.conj_mul']
  ring

/-- `⟨Ψ₀, |X| Ψ₀⟩ = E_q |x_q|` (paper: Lemma `lem:sectornorm`(c)). -/
theorem expect_sectorDiag_abs (ψ : Cfg N → ℂ) (hψ : star ψ ⬝ᵥ ψ = 1) (x : Cfg N → ℂ) :
    expect (∑ q : Cfg N, ((‖x q‖ : ℝ) : ℂ) • cplx (gaussProj q)) (initZ ψ) =
      ((avgQ (fun q : Cfg N => ‖x q‖) : ℝ) : ℂ) := by
  rw [sum_smul_gaussProj]
  simp only [expect, dotProduct, mulVec_diagonal, Fintype.sum_prod_type, initZ, Pi.star_apply]
  have e : ∀ n : Cfg N, ∑ y : Cfg N, star ((((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ * ψ n) *
      (((‖x (chargeOf n y)‖ : ℝ) : ℂ) * ((((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ * ψ n)) =
      star ((((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹) * (((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ *
        (star (ψ n) * ψ n) * ∑ q : Cfg N, ((‖x q‖ : ℝ) : ℂ) := by
    intro n
    conv_rhs => rw [Finset.mul_sum, ← sum_comp_chargeOf n]
    refine Finset.sum_congr rfl fun y _ => ?_
    rw [star_mul']
    ring
  simp only [e]
  rw [← Finset.sum_mul, ← Finset.mul_sum, sec_norm_sq_prefactor]
  have hψ' : ∑ n, star (ψ n) * ψ n = 1 := hψ
  rw [hψ', avgQ]
  push_cast
  ring

/-- `4^{-N} Tr |X| = E_q |x_q|` (paper: Lemma `lem:sectornorm`(c)). -/
theorem trace_sectorDiag_abs (x : Cfg N → ℂ) :
    ((4 : ℂ) ^ N)⁻¹ * (∑ q : Cfg N, ((‖x q‖ : ℝ) : ℂ) • cplx (gaussProj q)).trace =
      ((avgQ (fun q : Cfg N => ‖x q‖) : ℝ) : ℂ) := by
  rw [sum_smul_gaussProj, trace_diagonal, Fintype.sum_prod_type]
  have e : ∀ n : Cfg N, ∑ y : Cfg N, ((‖x (chargeOf n y)‖ : ℝ) : ℂ) =
      ∑ q : Cfg N, ((‖x q‖ : ℝ) : ℂ) := fun n =>
    sum_comp_chargeOf n (fun q => ((‖x q‖ : ℝ) : ℂ))
  simp only [e]
  rw [Finset.sum_const, card_univ, card_cfg, nsmul_eq_mul, avgQ, four_pow]
  push_cast
  field_simp

end DFL
