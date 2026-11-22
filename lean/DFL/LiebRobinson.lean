import DFL.OneBody
import DFL.FreeDyn

/-!
# Sector-averaged zero-velocity Lieb–Robinson bounds (paper: Proposition `prop:LR`)

In the dual frame the dressed fermions of the LGT are the free fermions `c_x` of each charge
sector, evolving by `c_x(t) = ∑_y (e^{-ith_q})_{xy} c_y` (`DFL.heisenberg_ann`). Under Hypothesis
(EC) their anticommutators and their commutators with local observables decay exponentially in the
distance *on average over the sectors*, uniformly in time. (No state-independent operator-norm
bound holds: the clean sector `q ≡ 0` is ballistic.)

Norms of matrices are `ℓ²` operator norms (`Matrix.Norms.L2Operator`).
-/

open Matrix Finset
open scoped Matrix.Norms.L2Operator

noncomputable section

namespace DFL

variable {N : ℕ}

/-- **Sector-averaged anticommutator bound**: `avg_q |{c_x(t), c_y†}_q| ≤ F(|x - y|)` for every
time, where `{c_x(t), c_y†}_q = (e^{-ith_q})_{xy}`. -/
theorem avg_norm_evo_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (x y : Fin N) (t : ℝ) :
    avgQ (fun q : Cfg N => ‖evo (hq J h q) t x y‖) ≤ F (Nat.dist x.val y.val) := by
  sorry

/-- Uniform-in-time version, with the supremum over `t` inside the average. -/
theorem avg_iSup_norm_evo_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (x y : Fin N) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ, ‖evo (hq J h q) t x y‖) ≤ F (Nat.dist x.val y.val) := by
  sorry

/-- The Heisenberg-evolved annihilation operator in the sector `q`. -/
def annT (J h : ℝ) (q : Cfg N) (x : Fin N) (t : ℝ) : Matrix (Cfg N) (Cfg N) ℂ :=
  (evo (dGamma (hq J h q)) t)ᴴ * cplx (ann x) * evo (dGamma (hq J h q)) t

/-- Annihilation operators have operator norm at most one. -/
theorem norm_cplx_ann_le (y : Fin N) : ‖cplx (ann (N := N) y)‖ ≤ 1 := by
  sorry

/-- **Sector-averaged zero-velocity Lieb–Robinson bound**: if `B` commutes with every `c_y`,
`y ∉ Y` (e.g. an even observable supported in `Y`), then for every time
`avg_q ‖[c_x(t), B]‖ ≤ 2 ‖B‖ ∑_{y ∈ Y} F(|x - y|)`. -/
theorem avg_norm_comm_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (Y : Finset (Fin N))
    (B : Matrix (Cfg N) (Cfg N) ℂ) (hB : ∀ y ∉ Y, cplx (ann y) * B = B * cplx (ann y))
    (x : Fin N) (t : ℝ) :
    avgQ (fun q : Cfg N => ‖annT J h q x t * B - B * annT J h q x t‖) ≤
      2 * ‖B‖ * ∑ y ∈ Y, F (Nat.dist x.val y.val) := by
  sorry

end DFL
