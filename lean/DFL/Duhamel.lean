import DFL.OneBody

/-!
# Duhamel sensitivity to distant charges (paper: Lemma `lem:duhamel`)

Flipping the charges on a set `D` changes the one-body propagator, acting on a site `z`, by an
amount controlled by the propagation from `z` to `D`:
`‖(U_q(s) - U_{q'}(s)) δ_z‖ ≤ 4|h| |s| ∑_{x ∈ D} Q_{q'}(x, z)`, where `U_q(s) = e^{-ish_q}`.
Under Hypothesis (EC) the Bernoulli average is at most `4|h| |s| ∑_{x ∈ D} F(|x - z|)`: charges at
distance `R` from `z` influence the dynamics near `z` only after times of order `e^{ηR}`.

Vector norms are Euclidean (`EuclideanSpace ℂ (Fin N)`).
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Flipping a configuration on a set `D`. -/
def flipOn (D : Finset (Fin N)) (q : Cfg N) : Cfg N := fun x => if x ∈ D then q x + 1 else q x

/-- The Euclidean norm of the column `U δ_z` of a matrix. -/
def colNorm (U : Matrix (Fin N) (Fin N) ℂ) (z : Fin N) : ℝ :=
  ‖(WithLp.toLp 2 (fun x => U x z) : EuclideanSpace ℂ (Fin N))‖

/-- **Duhamel bound** (deterministic): if `q'` is `q` flipped on `D`, then
`‖(U_q(s) - U_{q'}(s)) δ_z‖ ≤ 4|h| |s| ∑_{x ∈ D} Q_{q'}(x, z)`. -/
theorem colNorm_evo_sub_le (J h : ℝ) (q : Cfg N) (D : Finset (Fin N)) (z : Fin N) (s : ℝ) :
    colNorm (evo (hq J h q) s - evo (hq J h (flipOn D q)) s) z ≤
      4 * |h| * |s| *
        ∑ x ∈ D, corrQ (hq J h (flipOn D q)) (hq_isHermitian J h (flipOn D q)) x z := by
  sorry

/-- **Averaged Duhamel bound** under Hypothesis (EC), uniformly for `|s| ≤ t`. -/
theorem avg_colNorm_evo_sub_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (D : Finset (Fin N))
    (z : Fin N) (t s : ℝ) (hs : |s| ≤ t) :
    avgQ (fun q : Cfg N => colNorm (evo (hq J h q) s - evo (hq J h (flipOn D q)) s) z) ≤
      4 * |h| * t * ∑ x ∈ D, F (Nat.dist x.val z.val) := by
  sorry

end DFL
