import DFL.GaussLaw
import DFL.FreeDyn

/-!
# The duality dictionary, dressed fermions and their dynamics
(paper: Theorem `thm:duality`(c),(d),(e))

`U` is the permutation `|n, q⟩ ↦ |n⟩ ⊗ |x(n, q)⟩` (`dualEquiv`); `U† X U` is
`X.submatrix dualEquiv dualEquiv`. In the dual frame `c_j = ann j ⊗ 1`, `τˣ_j = diag(sgn q_j)` and `τᶻ_j`
flips the charge `q_j` (`τᶻ_N := 1`). `flipSet T` flips the labels in `T`.

* (c) **Dictionary**: `U† f_j U = c_j τᶻ_j` (`ann_dual`), `U† σᶻ_{b_j} U = τᶻ_j τᶻ_{j+1}` (`sz_dual`),
  `U† ∏_{k ≥ j} σᶻ_{b_k} U = τᶻ_j` (`szTail_dual`), `U† σˣ_{b_j} U = ∏_{k ≤ j} τˣ_k (-1)^{n_k}`
  (`sx_dual`), `U† σˣ_{b_{j-1}} σˣ_{b_j} U = τˣ_j (-1)^{n_j}` (`elec_dual`),
  `U† σᶻ_{b_j} f_j† f_{j+1} U = c_j† c_{j+1}` (`hopTerm_dual`), `U† (-1)^{N_f} U = (-1)^{N_c}`
  (`parity_dual`).
* (d) **Dressed fermions** `c̃_j = U c_j U† = f_j ∏_{k ≥ j} σᶻ_{b_k}` (`dressed`, `dressed_eq`): they
  satisfy the CAR (`dressed_anticomm`, `dressed_CAR`), commute with every `G_k` (`dressed_comm_gauss`),
  are odd (`dressed_odd`), `c̃_j† c̃_j = n_j` (`dressed_number`), `c̃_x† c̃_y = f_x† (∏_{x ≤ k < y} σᶻ_{b_k}) f_y`
  for `x < y` (`dressed_hop`), and `H = -J ∑_j (c̃_j† c̃_{j+1} + h.c.) - h ∑_j G_j (-1)^{n_j}`
  (`H_dressed`).
* (e) **Dynamics**: `U† e^{-itH} U = ⊕_q e^{iht ∑_j q_j} Γ(e^{-ith_q})` (`evo_dual`) and the Heisenberg
  evolution of the dressed fermions `e^{itH} c̃_x e^{-itH} = ∑_q P_q ∑_y (e^{-ith_q})_{xy} c̃_y`
  (`heisenberg_dressed`).
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The flip of the labels in `T`: `|y⟩ ↦ |y + 1_T⟩`. -/
def flipSet (T : Finset (Fin N)) : Matrix (Cfg N) (Cfg N) ℝ := fun y' y =>
  if y' = y + fockCfg T then 1 else 0

/-- The fermion parity `(-1)^{N_f}` on the matter Fock space. -/
def parityF : Matrix (Cfg N) (Cfg N) ℝ := diagonal (fun n => sgn (∑ k, n k))

/-- The sites `k ≥ j`. -/
def tailOf (j : Fin N) : Finset (Fin N) := univ.filter (fun k => j ≤ k)

/-- The sites `x ≤ k < y`. -/
def between (x y : Fin N) : Finset (Fin N) := univ.filter (fun k => x ≤ k ∧ k < y)

/-! ### (c) The dictionary -/

/-- `σᶻ_{b_j}` is the flip of the link `b_j`. -/
theorem sz_eq_flipSet (j : Fin N) : sz j = flipSet {j} := by
  sorry

/-- `U† f_j U = c_j τᶻ_j`. -/
theorem ann_dual (j : Fin N) :
    (ann j ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)).submatrix (dualEquiv (N := N)) dualEquiv =
      ann j ⊗ₖ flipSet {j} := by
  sorry

/-- `U† σᶻ_{b_j} U = τᶻ_j τᶻ_{j+1}` (`τᶻ_N := 1`). -/
theorem sz_dual (j : Fin N) :
    ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ sz j).submatrix (dualEquiv (N := N)) dualEquiv =
      (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flipSet (univ.filter (fun k => k = j ∨ k.val = j.val + 1)) := by
  sorry

/-- `U† (∏_{k ≥ j} σᶻ_{b_k}) U = τᶻ_j`. -/
theorem szTail_dual (j : Fin N) :
    ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flipSet (tailOf j)).submatrix (dualEquiv (N := N)) dualEquiv =
      (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flipSet {j} := by
  sorry

/-- `U† σˣ_{b_j} U = ∏_{k ≤ j} τˣ_k (-1)^{n_k}`. -/
theorem sx_dual (j : Fin N) :
    ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ sx j).submatrix (dualEquiv (N := N)) dualEquiv =
      diagonal (fun a : Cfg N × Cfg N => ∏ k ∈ univ.filter (· ≤ j), sgn (a.2 k) * sgn (a.1 k)) := by
  sorry

/-- `U† σˣ_{b_{j-1}} σˣ_{b_j} U = τˣ_j (-1)^{n_j}`. -/
theorem elec_dual (j : Fin N) :
    ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ elec j).submatrix (dualEquiv (N := N)) dualEquiv =
      diagonal (fun a : Cfg N × Cfg N => sgn (a.2 j) * sgn (a.1 j)) := by
  sorry

/-- `U† (σᶻ_{b_j} f_j† f_{j+1}) U = c_j† c_{j+1}` (`j < N - 1`). -/
theorem hopTerm_dual (j : Fin N) (hj : j.val + 1 < N) :
    ((cre j * ann ⟨j.val + 1, hj⟩) ⊗ₖ sz j).submatrix (dualEquiv (N := N)) dualEquiv =
      (cre j * ann ⟨j.val + 1, hj⟩) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) := by
  sorry

/-- `U† (-1)^{N_f} U = (-1)^{N_c}`: `U` is even. -/
theorem parity_dual :
    (parityF ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)).submatrix (dualEquiv (N := N)) dualEquiv =
      parityF ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) := by
  sorry

/-! ### (d) Dressed fermions -/

/-- The dressed fermion `c̃_j = f_j ∏_{k ≥ j} σᶻ_{b_k}` (Wilson line to the boundary). -/
def dressed (j : Fin N) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ := ann j ⊗ₖ flipSet (tailOf j)

/-- **`c̃_j = U c_j U†`** (paper: eq. `eq:dressed`). -/
theorem dressed_eq (j : Fin N) :
    dressed j =
      (ann j ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)).submatrix (dualEquiv (N := N)).symm dualEquiv.symm := by
  sorry

/-- CAR: `{c̃_x, c̃_y} = 0`. -/
theorem dressed_anticomm (x y : Fin N) : dressed x * dressed y + dressed y * dressed x = 0 := by
  sorry

/-- CAR: `{c̃_x, c̃_y†} = δ_{xy}`. -/
theorem dressed_CAR (x y : Fin N) :
    dressed x * (dressed y)ᵀ + (dressed y)ᵀ * dressed x = if x = y then 1 else 0 := by
  sorry

/-- The dressed fermions are gauge invariant. -/
theorem dressed_comm_gauss (j k : Fin N) : dressed j * gauss k = gauss k * dressed j := by
  sorry

/-- The dressed fermions are odd. -/
theorem dressed_odd (j : Fin N) :
    dressed j * (parityF ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) =
      -((parityF ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) * dressed j) := by
  sorry

/-- `c̃_j† c̃_j = n_j`. -/
theorem dressed_number (j : Fin N) :
    (dressed j)ᵀ * dressed j = numOp {j} ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) := by
  sorry

/-- **Gauge-invariant hopping** (paper: Theorem `thm:duality`(d)): for `x < y`,
`c̃_x† c̃_y = f_x† (∏_{x ≤ k < y} σᶻ_{b_k}) f_y`. -/
theorem dressed_hop (x y : Fin N) (hxy : x < y) :
    (dressed x)ᵀ * dressed y =
      (cre x ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) *
        ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flipSet (between x y)) *
          (ann y ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) := by
  sorry

/-- **`H` in terms of dressed fermions** (paper: eq. `eq:Hdressed`):
`H = -J ∑_j (c̃_j† c̃_{j+1} + c̃_{j+1}† c̃_j) - h ∑_j G_j (-1)^{n_j}`. -/
theorem H_dressed (J h : ℝ) :
    H J h =
      (-J) • ∑ j : Fin N, (if hj : j.val + 1 < N then
          (dressed j)ᵀ * dressed ⟨j.val + 1, hj⟩ + (dressed ⟨j.val + 1, hj⟩)ᵀ * dressed j else 0) -
        h • ∑ j : Fin N, gauss j * (diagonal (fun n : Cfg N => sgn (n j)) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) := by
  sorry

/-! ### (e) Dynamics -/

/-- **The dynamics in the dual frame** (paper: Theorem `thm:duality`(e)):
`U† e^{-itH} U = ⊕_q e^{iht ∑_j q_j} Γ(e^{-ith_q})`. -/
theorem evo_dual (J h : ℝ) (t : ℝ) :
    (evo (H J h) t).submatrix (dualEquiv (N := N)) dualEquiv =
      blockDiagonal (fun q : Cfg N =>
        Complex.exp (Complex.I * (h * t * ∑ j, sgn (q j) : ℝ)) • evo (dGamma (hq J h q)) t) := by
  sorry

/-- **Heisenberg evolution of the dressed fermions** (paper: eq. `eq:heisenberg`):
`e^{itH} c̃_x e^{-itH} = ∑_q P_q ∑_y (e^{-ith_q})_{xy} c̃_y`. -/
theorem heisenberg_dressed (J h : ℝ) (x : Fin N) (t : ℝ) :
    (evo (H J h) t)ᴴ * cplx (dressed x) * evo (H J h) t =
      ∑ q : Cfg N, cplx (gaussProj q) * ∑ y : Fin N, evo (hq J h q) t x y • cplx (dressed y) := by
  sorry

end DFL
