import FermionsLean.CutRank
import FermionsLean.EntropyBound
import Mathlib.Analysis.Matrix.Spectrum

/-!
# A Gaussian "min-cut" bound: `S_A ≤ log 2 · rank C_AB`

For a Slater determinant with correlation matrix `C` on `A ⊕ B`, the entanglement entropy of `A`
is `S_A = ∑ₖ binEntropy ν_k`, where `ν_k` are the eigenvalues of the Hermitian block `C_AA`.

* `rank_sub_mul_self_eq_card`: for Hermitian `M`, `rank (M - M²)` counts the eigenvalues of `M`
  different from `0` and `1`;
* `entropy_le_log_two_mul_rank`: `S_A ≤ log 2 · rank C_AB`.

Combined with `rank_toBlocks₁₂_measure_le` and `rank_toBlocks₁₂_blockDiag_le`, this says that the
bound `log 2 · rank C_AB` can only grow through gates that straddle the cut: it never grows under
projective occupation measurements (any site, any outcome) or under gates acting on one side.
-/

open Matrix Real
open scoped ComplexOrder
open Unitary

namespace MonitoredFermions

variable {A B 𝕜 : Type*} [Fintype A] [Fintype B] [DecidableEq A] [DecidableEq B] [RCLike 𝕜]

/-- For Hermitian `M`, `rank (M - M * M)` is the number of eigenvalues of `M` not in `{0, 1}`. -/
theorem rank_sub_mul_self_eq_card {M : Matrix A A 𝕜} (hM : M.IsHermitian) :
    (M - M * M).rank =
      Fintype.card {k // hM.eigenvalues k ≠ 0 ∧ hM.eigenvalues k ≠ 1} := by
  have hD : M - M * M = conjStarAlgAut 𝕜 _ hM.eigenvectorUnitary
      (diagonal fun k => ((hM.eigenvalues k : 𝕜) - (hM.eigenvalues k : 𝕜) *
        (hM.eigenvalues k : 𝕜))) := by
    conv_lhs => rw [hM.spectral_theorem]
    rw [← map_mul, ← map_sub, diagonal_mul_diagonal, diagonal_sub]
    rfl
  rw [hD, conjStarAlgAut_apply, ← Unitary.coe_star,
    rank_mul_eq_left_of_isUnit_det _ _ (UnitaryGroup.det_isUnit _),
    rank_mul_eq_right_of_isUnit_det _ _ (UnitaryGroup.det_isUnit _), rank_diagonal]
  refine Fintype.card_congr (Equiv.subtypeEquivRight fun k => ?_)
  have e : (hM.eigenvalues k : 𝕜) - (hM.eigenvalues k : 𝕜) * (hM.eigenvalues k : 𝕜) =
      ((hM.eigenvalues k * (1 - hM.eigenvalues k) : ℝ) : 𝕜) := by
    push_cast
    ring
  rw [e, Ne, RCLike.ofReal_eq_zero, mul_eq_zero, not_or, sub_eq_zero, eq_comm (a := (1 : ℝ))]

omit [Fintype A] [Fintype B] [DecidableEq A] [DecidableEq B] in
lemma isHermitian_toBlocks₁₁ {C : Matrix (A ⊕ B) (A ⊕ B) 𝕜} (hH : C.IsHermitian) :
    C.toBlocks₁₁.IsHermitian := by
  ext a a'
  simp [toBlocks₁₁, conjTranspose_apply, hH.apply]

omit [DecidableEq B] in
/-- **Gaussian min-cut bound.** For a Slater determinant (Hermitian projector `C`), the
entanglement entropy of `A`, `S_A = ∑ₖ binEntropy ν_k` with `ν_k` the eigenvalues of `C_AA`,
is at most `log 2 · rank C_AB`. -/
theorem entropy_le_log_two_mul_rank {C : Matrix (A ⊕ B) (A ⊕ B) 𝕜} (hH : C.IsHermitian)
    (hC : C * C = C) :
    ∑ k, binEntropy ((isHermitian_toBlocks₁₁ hH).eigenvalues k) ≤
      log 2 * C.toBlocks₁₂.rank := by
  set hA := isHermitian_toBlocks₁₁ hH
  have h1 := entropy_le_log_two_mul_card Finset.univ hA.eigenvalues
  have h2 : ((Finset.univ.filter fun k => hA.eigenvalues k ≠ 0 ∧ hA.eigenvalues k ≠ 1).card :
      ℝ) = C.toBlocks₁₂.rank := by
    rw [← rank_block_defect hH hC, rank_sub_mul_self_eq_card hA, Fintype.card_subtype]
  rw [h2] at h1
  exact h1

end MonitoredFermions
