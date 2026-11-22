import Mathlib

/-!
# The binary entropy (paper: Lemma `lem:binent`)

`h(ν) ≤ 2 ln 2 √(ν(1-ν))` on `[0, 1]`, with equality exactly at `ν ∈ {0, ½, 1}`. With the exact
formula `S(ρ_A) = ∑_i h(ν_i)` for Slater determinants and the trace-norm bound, this gives the area
law with the constant `2 ln 2` of the paper.
-/

noncomputable section

namespace DFL

/-- **Binary entropy bound** (paper: Lemma `lem:binent`). -/
theorem binEntropy_le_two_log_two_sqrt (ν : ℝ) (h0 : 0 ≤ ν) (h1 : ν ≤ 1) :
    Real.binEntropy ν ≤ 2 * Real.log 2 * √(ν * (1 - ν)) := by
  sorry

/-- **Equality cases** (paper: Lemma `lem:binent`): equality holds exactly for `ν ∈ {0, ½, 1}`. -/
theorem binEntropy_eq_two_log_two_sqrt_iff (ν : ℝ) (h0 : 0 ≤ ν) (h1 : ν ≤ 1) :
    Real.binEntropy ν = 2 * Real.log 2 * √(ν * (1 - ν)) ↔ ν = 0 ∨ ν = 1 / 2 ∨ ν = 1 := by
  sorry

end DFL
