import DFL.DualityFull
import DFL.GaussLaw

/-!
# Where the Wilson line ends (paper: Remark `rem:strings`)

The dressed fermion `c̃_j = f_j ∏_{k ≥ j} σᶻ_{b_k}` (`dressed`) carries its Wilson line to the boundary link.

* A left-anchored string does not work: `c'_j = (∏_{k<j} σᶻ_{b_k}) f_j = c̃_j ∏_k σᶻ_{b_k}` (`leftDressed_eq`)
  satisfies the CAR (`leftDressed_anticomm`, `leftDressed_CAR`) but anticommutes with the first Gauss
  operator `G_1` (`leftDressed_anticomm_gauss0`; 0-indexed `gauss 0`): for `j = 1` through `(-1)^{n_1}`, for
  `j ≥ 2` through `σᶻ_{b_1}`.
* The dressed fermions anticommute with the boundary flux `σˣ_{b_N}` (`dressed_anticomm_sxL_last`), consistently
  with `∏_j G_j = (-1)^{N_f} σˣ_{b_N}` (`gauss_prod`) and `dressed_odd`.
* `{c̃_j, G_j, ∏_{k ≥ j} σᶻ_{b_k}}` is a CAR ⊗ Pauli system: the strings `T_j = ∏_{k ≥ j} σᶻ_{b_k}` (`tailString`,
  SKKM's `τ̂ˣ_j`) are commuting involutions, `G_i T_j = (-1)^{δ_{ij}} T_j G_i` (`gauss_tailString`), and both
  commute with the dressed fermions (`dressed_comm_gauss`, `dressed_comm_tailString`).
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The left-anchored string fermion `c'_j = (∏_{k<j} σᶻ_{b_k}) f_j`. -/
def leftDressed (j : Fin N) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  ann j ⊗ₖ flipSet (univ.filter (fun k => k < j))

/-- The string `T_j = ∏_{k ≥ j} σᶻ_{b_k}`. -/
def tailString (j : Fin N) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flipSet (tailOf j)

/-! ### Auxiliary lemmas -/

namespace StringsAux

/-- `{k ≥ j} ∆ univ = {k < j}`. -/
lemma tailOf_symmDiff_univ (j : Fin N) :
    symmDiff (tailOf j) univ = univ.filter (fun k => k < j) := by
  ext k
  simp [tailOf, Finset.mem_symmDiff, not_le]

/-- `c'_j† = c_j† ⊗ ∏_{k<j} σᶻ_{b_k}`. -/
lemma leftDressed_transpose (j : Fin N) :
    (leftDressed j)ᵀ = cre j ⊗ₖ flipSet (univ.filter (fun k => k < j)) := by
  rw [leftDressed, DualityFullAux.kronecker_transpose', DualityFullAux.flipSet_transpose]
  rfl

/-- If `d a = c * d b` on the support of `M`, then `diag(d) M = c • (M diag(d))`. -/
lemma diagonal_mul_eq_smul_of_support {ι : Type*} [Fintype ι] [DecidableEq ι] (d : ι → ℝ)
    (M : Matrix ι ι ℝ) (c : ℝ) (h : ∀ a b, M a b ≠ 0 → d a = c * d b) :
    diagonal d * M = c • (M * diagonal d) := by
  ext a b
  rw [diagonal_mul, Matrix.smul_apply, mul_diagonal, smul_eq_mul]
  by_cases hab : M a b = 0
  · rw [hab]; ring
  · rw [h a b hab]; ring

/-- If `d a = -d b` on the support of `M`, then `M diag(d) = -(diag(d) M)`. -/
lemma mul_diagonal_eq_neg_of_support {ι : Type*} [Fintype ι] [DecidableEq ι] (d : ι → ℝ)
    (M : Matrix ι ι ℝ) (h : ∀ a b, M a b ≠ 0 → d a = -d b) :
    M * diagonal d = -(diagonal d * M) := by
  rw [diagonal_mul_eq_smul_of_support d M (-1) (fun a b hab => by rw [h a b hab]; ring),
    neg_one_smul, neg_neg]

/-- Support of `c_j ⊗ flipSet T`: `n_j = 1`, `n' = n - e_j` and `x' = x + 1_T`. -/
lemma ann_kronecker_flipSet_ne_zero {j : Fin N} {T : Finset (Fin N)} {n' x' n x : Cfg N}
    (h : (ann j ⊗ₖ flipSet T) (n', x') (n, x) ≠ 0) :
    n j = 1 ∧ n' = Function.update n j 0 ∧ x' = x + fockCfg T := by
  rw [kronecker_apply] at h
  by_cases ha : n j = 1 ∧ n' = Function.update n j 0
  · refine ⟨ha.1, ha.2, ?_⟩
    by_contra hx
    apply h
    simp only [flipSet]
    rw [ite_eq_right hx, mul_zero]
  · exact (h (by rw [FreeDynAux.ann_apply, ite_eq_right ha, zero_mul])).elim

/-- Support of `1 ⊗ flipSet T`: `n' = n` and `x' = x + 1_T`. -/
lemma one_kronecker_flipSet_ne_zero {T : Finset (Fin N)} {n' x' n x : Cfg N}
    (h : ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flipSet T) (n', x') (n, x) ≠ 0) :
    n' = n ∧ x' = x + fockCfg T := by
  rw [kronecker_apply] at h
  by_cases hn : n' = n
  · refine ⟨hn, ?_⟩
    by_contra hx
    apply h
    simp only [flipSet]
    rw [ite_eq_right hx, mul_zero]
  · exact (h (by rw [Matrix.one_apply_ne hn, zero_mul])).elim

end StringsAux

/-- `c'_j = c̃_j ∏_k σᶻ_{b_k}`. -/
theorem leftDressed_eq (j : Fin N) :
    leftDressed j = dressed j * ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flipSet univ) := by
  rw [leftDressed, dressed, ← mul_kronecker_mul, Matrix.mul_one, DualityFullAux.flipSet_mul,
    StringsAux.tailOf_symmDiff_univ]

/-- CAR: `{c'_x, c'_y} = 0`. -/
theorem leftDressed_anticomm (x y : Fin N) :
    leftDressed x * leftDressed y + leftDressed y * leftDressed x = 0 := by
  simp only [leftDressed]
  rw [← mul_kronecker_mul, ← mul_kronecker_mul, DualityFullAux.flipSet_mul,
    DualityFullAux.flipSet_mul, symmDiff_comm (univ.filter fun k => k < y), ← add_kronecker,
    ann_mul_ann_add, zero_kronecker]

/-- CAR: `{c'_x, c'_y†} = δ_{xy}`. -/
theorem leftDressed_CAR (x y : Fin N) :
    leftDressed x * (leftDressed y)ᵀ + (leftDressed y)ᵀ * leftDressed x = if x = y then 1 else 0 := by
  rw [StringsAux.leftDressed_transpose]
  simp only [leftDressed]
  rw [← mul_kronecker_mul, ← mul_kronecker_mul, DualityFullAux.flipSet_mul,
    DualityFullAux.flipSet_mul, symmDiff_comm (univ.filter fun k => k < y), ← add_kronecker,
    ann_mul_cre_add]
  split_ifs with hxy
  · subst hxy
    rw [symmDiff_self, Finset.bot_eq_empty, DualityFullAux.flipSet_empty, one_kronecker_one]
  · rw [zero_kronecker]

/-- **The left-anchored string fermion is not gauge invariant**: it anticommutes with `G_1`. -/
theorem leftDressed_anticomm_gauss0 (hN : 0 < N) (j : Fin N) :
    leftDressed j * gauss ⟨0, hN⟩ = -(gauss ⟨0, hN⟩ * leftDressed j) := by
  unfold gauss
  refine StringsAux.mul_diagonal_eq_neg_of_support _ _ ?_
  rintro ⟨n', x'⟩ ⟨n, x⟩ hab
  obtain ⟨hn, rfl, rfl⟩ := StringsAux.ann_kronecker_flipSet_ne_zero hab
  have hp : ∀ z : Cfg N, prevLink z ⟨0, hN⟩ = 0 := fun z => by simp [prevLink]
  simp only [hp, zero_add, Pi.add_apply]
  by_cases hj : j = ⟨0, hN⟩
  · subst hj
    have h0 : fockCfg (univ.filter (fun k : Fin N => k < ⟨0, hN⟩)) ⟨0, hN⟩ = 0 := by
      simp [fockCfg]
    rw [Function.update_self, h0, add_zero, hn, FreeDynAux.sgn_zero, FreeDynAux.sgn_one]
    ring
  · have hlt : (⟨0, hN⟩ : Fin N) < j := by
      rw [Fin.lt_def]
      exact Nat.pos_of_ne_zero (fun h => hj (Fin.ext h))
    have h0 : fockCfg (univ.filter (fun k : Fin N => k < j)) ⟨0, hN⟩ = 1 := by
      simp [fockCfg, hlt]
    rw [Function.update_of_ne (Ne.symm hj), h0, GaussLawAux.sgn_add, FreeDynAux.sgn_one]
    ring

/-- The dressed fermions anticommute with the boundary flux `σˣ_{b_N}`. -/
theorem dressed_anticomm_sxL_last (hN : 0 < N) (j : Fin N) :
    dressed j * sxL ⟨N - 1, by omega⟩ = -(sxL ⟨N - 1, by omega⟩ * dressed j) := by
  rw [sxL, sx, GaussLawAux.one_kronecker_diagonal]
  refine StringsAux.mul_diagonal_eq_neg_of_support _ _ ?_
  rintro ⟨n', x'⟩ ⟨n, x⟩ hab
  obtain ⟨-, -, rfl⟩ := StringsAux.ann_kronecker_flipSet_ne_zero hab
  have hle : j ≤ ⟨N - 1, by omega⟩ := by
    rw [Fin.le_def]
    show j.val ≤ N - 1
    have := j.isLt
    omega
  have hmem : fockCfg (tailOf j) ⟨N - 1, by omega⟩ = 1 := by
    simp [fockCfg, tailOf, hle]
  simp only [Pi.add_apply, hmem, GaussLawAux.sgn_add, FreeDynAux.sgn_one]
  ring

/-- `T_j² = 1`. -/
theorem tailString_mul_self (j : Fin N) : tailString j * tailString j = 1 := by
  rw [tailString, ← mul_kronecker_mul, Matrix.mul_one, DualityFullAux.flipSet_mul_self,
    one_kronecker_one]

/-- The strings commute. -/
theorem tailString_comm (i j : Fin N) : tailString i * tailString j = tailString j * tailString i := by
  simp only [tailString]
  rw [← mul_kronecker_mul, ← mul_kronecker_mul, DualityFullAux.flipSet_mul,
    DualityFullAux.flipSet_mul, symmDiff_comm]

/-- **Pauli relations**: `G_i T_j = (-1)^{δ_{ij}} T_j G_i`. -/
theorem gauss_tailString (i j : Fin N) :
    gauss i * tailString j = (if i = j then (-1 : ℝ) else 1) • (tailString j * gauss i) := by
  rw [GaussLawAux.gauss_eq]
  refine StringsAux.diagonal_mul_eq_smul_of_support _ _ _ ?_
  rintro ⟨n', x'⟩ ⟨n, x⟩ hab
  obtain ⟨hn, hx⟩ := StringsAux.one_kronecker_flipSet_ne_zero hab
  subst hn hx
  dsimp only
  rw [DualityFullAux.chargeOf_add_right, DualityFullAux.delta_tail, Pi.add_apply,
    GaussLawAux.sgn_add]
  by_cases hij : i = j
  · rw [ite_eq_left hij, show fockCfg {j} i = 1 by simp [fockCfg, hij], FreeDynAux.sgn_one]
    ring
  · rw [ite_eq_right hij, show fockCfg {j} i = 0 by simp [fockCfg, hij], FreeDynAux.sgn_zero]
    ring

/-- The dressed fermions commute with the strings. -/
theorem dressed_comm_tailString (i j : Fin N) :
    dressed i * tailString j = tailString j * dressed i := by
  rw [dressed, tailString, ← mul_kronecker_mul, ← mul_kronecker_mul, Matrix.mul_one,
    Matrix.one_mul, DualityFullAux.flipSet_mul, DualityFullAux.flipSet_mul, symmDiff_comm]

end DFL
