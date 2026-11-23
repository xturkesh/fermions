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

/-- `c'_j = c̃_j ∏_k σᶻ_{b_k}`. -/
theorem leftDressed_eq (j : Fin N) :
    leftDressed j = dressed j * ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flipSet univ) := by
  sorry

/-- CAR: `{c'_x, c'_y} = 0`. -/
theorem leftDressed_anticomm (x y : Fin N) :
    leftDressed x * leftDressed y + leftDressed y * leftDressed x = 0 := by
  sorry

/-- CAR: `{c'_x, c'_y†} = δ_{xy}`. -/
theorem leftDressed_CAR (x y : Fin N) :
    leftDressed x * (leftDressed y)ᵀ + (leftDressed y)ᵀ * leftDressed x = if x = y then 1 else 0 := by
  sorry

/-- **The left-anchored string fermion is not gauge invariant**: it anticommutes with `G_1`. -/
theorem leftDressed_anticomm_gauss0 (hN : 0 < N) (j : Fin N) :
    leftDressed j * gauss ⟨0, hN⟩ = -(gauss ⟨0, hN⟩ * leftDressed j) := by
  sorry

/-- The dressed fermions anticommute with the boundary flux `σˣ_{b_N}`. -/
theorem dressed_anticomm_sxL_last (hN : 0 < N) (j : Fin N) :
    dressed j * sxL ⟨N - 1, by omega⟩ = -(sxL ⟨N - 1, by omega⟩ * dressed j) := by
  sorry

/-- `T_j² = 1`. -/
theorem tailString_mul_self (j : Fin N) : tailString j * tailString j = 1 := by
  sorry

/-- The strings commute. -/
theorem tailString_comm (i j : Fin N) : tailString i * tailString j = tailString j * tailString i := by
  sorry

/-- **Pauli relations**: `G_i T_j = (-1)^{δ_{ij}} T_j G_i`. -/
theorem gauss_tailString (i j : Fin N) :
    gauss i * tailString j = (if i = j then (-1 : ℝ) else 1) • (tailString j * gauss i) := by
  sorry

/-- The dressed fermions commute with the strings. -/
theorem dressed_comm_tailString (i j : Fin N) :
    dressed i * tailString j = tailString j * dressed i := by
  sorry

end DFL
