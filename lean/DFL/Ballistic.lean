import DFL.InfVol

/-!
# The clean sector is ballistic on `ℤ` (paper: Remark `rem:gaussproj`(iii))

For the clean sector on `ℤ`, `h_ℤ = -JΔ + 2h` (`UZ J h (fun _ => 0)`, the exponential series of the
tridiagonal `hZ` of `InfVol.lean`), the mean-square displacement of a particle started at `x` is exactly
`∑_y (y - x)² |⟨δ_y, e^{-ith_ℤ} δ_x⟩|² = 2J²t²` for all `t` (`ballistic_clean`): on finitely supported
vectors `e^{ith_ℤ} X e^{-ith_ℤ} = X + tV`, with `V = i[h_ℤ, X]`, because `V` commutes with `h_ℤ`, and
`‖V δ_x‖² = 2J²`. Equivalently, in the power series of `t ↦ ⟨δ_x, e^{ith} (X - x)² e^{-ith} δ_x⟩` the
coefficients are `⟨δ_x, ad_h^n((X - x)²) δ_x⟩ (i^n / n!)`, and `ad_h³((X - x)²) = 0`.

Consequently there is no bound `sup_t |⟨δ_x, e^{-ith_ℤ} δ_y⟩| ≤ F(|x - y|)` with `∑_d d² F(d)² < ∞`
(`no_uniform_decay_clean`).
-/

open Filter Topology

noncomputable section

namespace DFL

/-- **The clean sector is ballistic** (paper: `rem:gaussproj`(iii)):
`∑_y (y - x)² |⟨δ_y, e^{-ith_ℤ} δ_x⟩|² = 2J²t²`. -/
theorem ballistic_clean (J h t : ℝ) (x : ℤ) :
    ∑' y : ℤ, ((y - x : ℤ) : ℝ) ^ 2 * ‖UZ J h (fun _ => 0) t y x‖ ^ 2 = 2 * J ^ 2 * t ^ 2 := by
  sorry

/-- The series in `ballistic_clean` converges. -/
theorem summable_ballistic_clean (J h t : ℝ) (x : ℤ) :
    Summable (fun y : ℤ => ((y - x : ℤ) : ℝ) ^ 2 * ‖UZ J h (fun _ => 0) t y x‖ ^ 2) := by
  sorry

/-- **No uniform decay in the clean sector** (paper: `rem:gaussproj`(iii)): for `J ≠ 0` there is no `F` with
`∑_d d² F(d)² < ∞` and `|⟨δ_y, e^{-ith_ℤ} δ_x⟩| ≤ F(|x - y|)` for all `t`, `x`, `y`. -/
theorem no_uniform_decay_clean (J h : ℝ) (hJ : J ≠ 0) (F : ℕ → ℝ)
    (hF : Summable (fun d : ℕ => (d : ℝ) ^ 2 * F d ^ 2)) :
    ¬ ∀ (t : ℝ) (x y : ℤ), ‖UZ J h (fun _ => 0) t y x‖ ≤ F (y - x).natAbs := by
  sorry

end DFL
