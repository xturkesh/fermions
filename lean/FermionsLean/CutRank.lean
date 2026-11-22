import FermionsLean.Measurement
import Mathlib.LinearAlgebra.Matrix.Rank
import Mathlib.Data.Matrix.Block

/-!
# Entanglement across a cut and the rank of the off-diagonal block

Split the sites into `A ⊕ B` and let `C` be the correlation matrix of a Slater determinant
(a Hermitian projector). Each eigenvalue `ν ∈ (0, 1)` of the block `C_AA` is an "entangled mode"
contributing `binEntropy ν ≤ log 2` to the entanglement entropy `S_A`, and the number of such
modes is `rank (C_AA - C_AA²)`. Here we prove the purely algebraic facts:

* `block_defect`: `C_AA - C_AA * C_AA = C_AB * C_BA` for every idempotent `C`;
* `rank_block_defect`: for a Hermitian projector, `rank (C_AA - C_AA²) = rank C_AB`;
* `rank_toBlocks₁₂_measOcc_inl_le` & co.: a projective occupation measurement, on either side
  of the cut and with either outcome, never increases `rank C_AB`;
* `rank_toBlocks₁₂_blockDiag_le`: (block-diagonal) gates acting inside `A` and inside `B` never
  increase `rank C_AB`.

Together with `S_A ≤ log 2 · #(entangled modes)` this is a Gaussian, linear-algebraic version of
the "min-cut" bound used for hybrid circuits: entanglement across the cut can only be created by
gates that straddle it.
-/

open Matrix

namespace MonitoredFermions

variable {A B 𝕜 : Type*} [Fintype A] [Fintype B] [DecidableEq A] [DecidableEq B] [Field 𝕜]

omit [DecidableEq A] [DecidableEq B] in
lemma toBlocks₁₁_mul (M N : Matrix (A ⊕ B) (A ⊕ B) 𝕜) :
    (M * N).toBlocks₁₁ = M.toBlocks₁₁ * N.toBlocks₁₁ + M.toBlocks₁₂ * N.toBlocks₂₁ := by
  ext a a'
  simp [toBlocks₁₁, toBlocks₁₂, toBlocks₂₁, Matrix.mul_apply, Fintype.sum_sum_type]

omit [DecidableEq A] [DecidableEq B] in
lemma toBlocks₁₂_mul (M N : Matrix (A ⊕ B) (A ⊕ B) 𝕜) :
    (M * N).toBlocks₁₂ = M.toBlocks₁₁ * N.toBlocks₁₂ + M.toBlocks₁₂ * N.toBlocks₂₂ := by
  ext a b
  simp [toBlocks₁₁, toBlocks₁₂, toBlocks₂₂, Matrix.mul_apply, Fintype.sum_sum_type]

omit [Fintype A] [Fintype B] [DecidableEq A] [DecidableEq B] in
lemma toBlocks₁₂_sub (M N : Matrix (A ⊕ B) (A ⊕ B) 𝕜) :
    (M - N).toBlocks₁₂ = M.toBlocks₁₂ - N.toBlocks₁₂ := rfl

omit [Fintype A] [Fintype B] in
lemma toBlocks₁₂_one : (1 : Matrix (A ⊕ B) (A ⊕ B) 𝕜).toBlocks₁₂ = 0 := by
  ext a b
  simp only [toBlocks₁₂, Matrix.of_apply, Matrix.one_apply, Matrix.zero_apply]
  simp

omit [DecidableEq A] [DecidableEq B] in
/-- For an idempotent `C`: `C_AA - C_AA * C_AA = C_AB * C_BA`. -/
theorem block_defect {C : Matrix (A ⊕ B) (A ⊕ B) 𝕜} (hC : C * C = C) :
    C.toBlocks₁₁ - C.toBlocks₁₁ * C.toBlocks₁₁ = C.toBlocks₁₂ * C.toBlocks₂₁ := by
  have h := congrArg toBlocks₁₁ hC
  rw [toBlocks₁₁_mul] at h
  rw [sub_eq_iff_eq_add']
  exact h.symm

section Star

variable [StarRing 𝕜]

omit [Fintype A] [Fintype B] [DecidableEq A] [DecidableEq B] in
lemma toBlocks₂₁_eq_conjTranspose {C : Matrix (A ⊕ B) (A ⊕ B) 𝕜} (hH : C.IsHermitian) :
    C.toBlocks₂₁ = (C.toBlocks₁₂)ᴴ := by
  ext b a
  simp [toBlocks₂₁, toBlocks₁₂, conjTranspose_apply, hH.apply]

omit [DecidableEq A] [DecidableEq B] in
/-- **Number of entangled modes = rank of the off-diagonal block.** For a Hermitian projector
`C`, `rank (C_AA - C_AA²) = rank C_AB`. -/
theorem rank_block_defect [PartialOrder 𝕜] [StarOrderedRing 𝕜]
    {C : Matrix (A ⊕ B) (A ⊕ B) 𝕜} (hH : C.IsHermitian) (hC : C * C = C) :
    (C.toBlocks₁₁ - C.toBlocks₁₁ * C.toBlocks₁₁).rank = C.toBlocks₁₂.rank := by
  rw [block_defect hC, toBlocks₂₁_eq_conjTranspose hH, rank_self_mul_conjTranspose]

end Star

/-! ### Measurements never increase `rank C_AB` -/

/-- Off-diagonal block after measuring a site `a₀ ∈ A` (outcome `1`). -/
theorem toBlocks₁₂_measOcc_inl (C : Matrix (A ⊕ B) (A ⊕ B) 𝕜) (a₀ : A) :
    (measOcc C (Sum.inl a₀)).toBlocks₁₂ =
      (1 - (C (Sum.inl a₀) (Sum.inl a₀))⁻¹ • (C.toBlocks₁₁ * single a₀ a₀ (1 : 𝕜))) *
        C.toBlocks₁₂ := by
  rw [Matrix.sub_mul, Matrix.one_mul, Matrix.smul_mul]
  ext a b
  simp only [Matrix.sub_apply, Matrix.smul_apply, mul_single_mul_apply, measOcc,
    Matrix.add_apply, mul_siteProj_mul_apply, siteProj_apply, toBlocks₁₂, toBlocks₁₁,
    Matrix.of_apply, smul_eq_mul]
  simp

/-- Off-diagonal block after measuring a site `b₀ ∈ B` (outcome `1`). -/
theorem toBlocks₁₂_measOcc_inr (C : Matrix (A ⊕ B) (A ⊕ B) 𝕜) (b₀ : B) :
    (measOcc C (Sum.inr b₀)).toBlocks₁₂ =
      C.toBlocks₁₂ *
        (1 - (C (Sum.inr b₀) (Sum.inr b₀))⁻¹ • (single b₀ b₀ (1 : 𝕜) * C.toBlocks₂₂)) := by
  rw [Matrix.mul_sub, Matrix.mul_one, Matrix.mul_smul, ← Matrix.mul_assoc]
  ext a b
  simp only [Matrix.sub_apply, Matrix.smul_apply, mul_single_mul_apply, measOcc,
    Matrix.add_apply, mul_siteProj_mul_apply, siteProj_apply, toBlocks₁₂, toBlocks₂₂,
    Matrix.of_apply, smul_eq_mul]
  simp

theorem rank_toBlocks₁₂_measOcc_inl_le (C : Matrix (A ⊕ B) (A ⊕ B) 𝕜) (a₀ : A) :
    (measOcc C (Sum.inl a₀)).toBlocks₁₂.rank ≤ C.toBlocks₁₂.rank := by
  rw [toBlocks₁₂_measOcc_inl]
  exact rank_mul_le_right _ _

theorem rank_toBlocks₁₂_measOcc_inr_le (C : Matrix (A ⊕ B) (A ⊕ B) 𝕜) (b₀ : B) :
    (measOcc C (Sum.inr b₀)).toBlocks₁₂.rank ≤ C.toBlocks₁₂.rank := by
  rw [toBlocks₁₂_measOcc_inr]
  exact rank_mul_le_left _ _

lemma toBlocks₁₂_measEmp (C : Matrix (A ⊕ B) (A ⊕ B) 𝕜) (s : A ⊕ B) :
    (measEmp C s).toBlocks₁₂ = -(measOcc (1 - C) s).toBlocks₁₂ := by
  rw [measEmp, toBlocks₁₂_sub, toBlocks₁₂_one, zero_sub]

omit [Fintype A] [Fintype B] in
lemma toBlocks₁₂_one_sub (C : Matrix (A ⊕ B) (A ⊕ B) 𝕜) :
    (1 - C).toBlocks₁₂ = -C.toBlocks₁₂ := by
  rw [toBlocks₁₂_sub, toBlocks₁₂_one, zero_sub]

lemma rank_neg' {m o : Type*} [Fintype o] [DecidableEq o] (M : Matrix m o 𝕜) :
    (-M).rank = M.rank := by
  have h : -M = M * (-1 : Matrix o o 𝕜) := by simp
  rw [h]
  refine rank_mul_eq_left_of_isUnit_det _ _ ?_
  rw [isUnit_iff_ne_zero, det_neg, det_one, mul_one]
  exact pow_ne_zero _ (neg_ne_zero.mpr one_ne_zero)

theorem rank_toBlocks₁₂_measEmp_le (C : Matrix (A ⊕ B) (A ⊕ B) 𝕜) (s : A ⊕ B) :
    (measEmp C s).toBlocks₁₂.rank ≤ C.toBlocks₁₂.rank := by
  rw [toBlocks₁₂_measEmp, rank_neg']
  calc (measOcc (1 - C) s).toBlocks₁₂.rank ≤ (1 - C).toBlocks₁₂.rank := by
        rcases s with a₀ | b₀
        · exact rank_toBlocks₁₂_measOcc_inl_le _ a₀
        · exact rank_toBlocks₁₂_measOcc_inr_le _ b₀
    _ = C.toBlocks₁₂.rank := by rw [toBlocks₁₂_one_sub, rank_neg']

/-- **Measurements never create entanglement across the cut**: for either outcome of the
measurement of any site `s`, `rank C_AB` does not increase. -/
theorem rank_toBlocks₁₂_measure_le (C : Matrix (A ⊕ B) (A ⊕ B) 𝕜) (s : A ⊕ B) :
    (measOcc C s).toBlocks₁₂.rank ≤ C.toBlocks₁₂.rank ∧
      (measEmp C s).toBlocks₁₂.rank ≤ C.toBlocks₁₂.rank := by
  refine ⟨?_, rank_toBlocks₁₂_measEmp_le C s⟩
  rcases s with a₀ | b₀
  · exact rank_toBlocks₁₂_measOcc_inl_le C a₀
  · exact rank_toBlocks₁₂_measOcc_inr_le C b₀

/-! ### Gates acting on one side of the cut -/

omit [DecidableEq A] [DecidableEq B] in
/-- Conjugating by block-diagonal matrices (gates acting inside `A` and inside `B`) never
increases `rank C_AB`. For unitaries `U = diag(U_A, U_B)` take `N = Uᴴ`. -/
theorem rank_toBlocks₁₂_blockDiag_le (C : Matrix (A ⊕ B) (A ⊕ B) 𝕜)
    (MA NA : Matrix A A 𝕜) (MB NB : Matrix B B 𝕜) :
    (fromBlocks MA 0 0 MB * C * fromBlocks NA 0 0 NB).toBlocks₁₂.rank ≤ C.toBlocks₁₂.rank := by
  have h : (fromBlocks MA 0 0 MB * C * fromBlocks NA 0 0 NB).toBlocks₁₂ =
      MA * C.toBlocks₁₂ * NB := by
    rw [toBlocks₁₂_mul, toBlocks₁₁_mul, toBlocks₁₂_mul]
    simp
  rw [h]
  exact (rank_mul_le_left _ _).trans (rank_mul_le_right _ _)

end MonitoredFermions
