import DFL.LRNorm

/-!
# No operator-norm zero-velocity bound (paper: Proposition `prop:LR`(d))

There is no function `G ≥ 0` with `∑_d d² G(d)² < ∞` such that `sup_t ‖[n_x(t), n_y]‖ ≤ G(|x - y|)` for all `N` and
all `x, y`; the same holds with `[n_x(t), n_y]` replaced by `{c̃_x(t), c̃_y†}` or by `[c̃_x(t), n_y]`
(`no_uniform_LR_nn`, `no_uniform_LR_anticomm`, `no_uniform_LR_cn`).

The operator norms are maxima over the sectors (`norm_comm_nn_eq`, `norm_dressedAnticomm_eq`, `norm_comm_number_eq`) and
include the clean sector, where the dressed fermion is a free particle on the chain. It moves ballistically:
`∑_r r² |U(t)_{r, c}|² ≥ (√2 |J t| - ε)²` in the middle `c` of a long chain (Heisenberg equations `Ẋ = V`, `V̇ = i[H, V]`,
with `[H, V]` supported at the two ends, which the particle reaches only with exponentially small amplitude), which is
incompatible with `|U(t)_{rc}| ≤ √2 G(|r - c|)` for large `|r - c|`.
-/

open Matrix Finset
open scoped Kronecker Matrix.Norms.L2Operator

noncomputable section

namespace DFL

/-- **No operator-norm zero-velocity bound for `[n_x(t), n_y]`** (paper: Proposition `prop:LR`(d)). -/
theorem no_uniform_LR_nn (J h : ℝ) (hJ : J ≠ 0) (G : ℕ → ℝ) (hG0 : ∀ d, 0 ≤ G d)
    (hG : Summable (fun d : ℕ => (d : ℝ) ^ 2 * G d ^ 2)) :
    ¬ ∀ (N : ℕ) (x y : Fin N) (t : ℝ), ‖numComm (N := N) J h x y t‖ ≤ G (Nat.dist x.val y.val) := by
  sorry

/-- **No operator-norm zero-velocity bound for `{c̃_x(t), c̃_y†}`** (paper: Proposition `prop:LR`(d)). -/
theorem no_uniform_LR_anticomm (J h : ℝ) (hJ : J ≠ 0) (G : ℕ → ℝ) (hG0 : ∀ d, 0 ≤ G d)
    (hG : Summable (fun d : ℕ => (d : ℝ) ^ 2 * G d ^ 2)) :
    ¬ ∀ (N : ℕ) (x y : Fin N) (t : ℝ),
      ‖dressedAnticomm (N := N) J h x y t‖ ≤ G (Nat.dist x.val y.val) := by
  sorry

/-- **No operator-norm zero-velocity bound for `[c̃_x(t), n_y]`** (paper: Proposition `prop:LR`(d)). -/
theorem no_uniform_LR_cn (J h : ℝ) (hJ : J ≠ 0) (G : ℕ → ℝ) (hG0 : ∀ d, 0 ≤ G d)
    (hG : Summable (fun d : ℕ => (d : ℝ) ^ 2 * G d ^ 2)) :
    ¬ ∀ (N : ℕ) (x y : Fin N) (t : ℝ),
      ‖heisT J h t (cplx (dressed x)) * cplx (numOp {y} ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) -
          cplx (numOp {y} ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) * heisT J h t (cplx (dressed x))‖ ≤
        G (Nat.dist x.val y.val) := by
  sorry

end DFL
