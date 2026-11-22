import DFL.DimerTransfer
import DFL.DimerTL

/-!
# Every dimer configuration spreads at least diffusively (paper: Proposition `prop:dimerlower`)

For a dimer configuration `v` (values `±λ`, `0 < λ < 1`, equal on the two sites of each dimer) on a finite
interval `I = [a, b]`, `H_v = -Δ_I + diag(v)` (`schrH`), `x ∈ I` and `T ≥ 1` with `[x - T, x + T] ⊆ I`, the
time-averaged moments
`M_p(T; x) = ∫_0^∞ (dt/T) e^{-t/T} ∑_y |y - x|^p |⟨δ_y, e^{-itH_v} δ_x⟩|²` (`dimerMom`) satisfy

* (Plancherel, eq. `eq:plancherel`) `M_p(T; x) = (2πT)⁻¹ ∫_ℝ ∑_y |y - x|^p |G_{E + i/(2T)}(y)|² dE`, with
  `G_z(y) = ⟨δ_y, (H_v - z)⁻¹ δ_x⟩` (`dimerMom_plancherel`);
* `M_p(T; x) ≥ c_{λ,p} T^{p-1}`, `c_{λ,p} = c₁² / (2^{p+3} π (p+1))`, `c₁ = 1/((3 + 2λ) K₁ e^{K₂})`, with `K₁, K₂`
  the constants of Lemma `lem:dimertransfer` (`dimer_lower_bound`).
-/

open Matrix MeasureTheory

noncomputable section

namespace DFL

/-- The time-averaged moment `M_p(T; x) = ∫_0^∞ (dt/T) e^{-t/T} ∑_y |y - x|^p |⟨δ_y, e^{-itH_v} δ_x⟩|²`. -/
def dimerMom (a b : ℤ) (v : ℤ → ℝ) (p T : ℝ) (x : Set.Icc a b) : ℝ :=
  ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) *
    ∑ y : Set.Icc a b, |((y.1 - x.1 : ℤ) : ℝ)| ^ p * ‖evo (schrH a b v) t y x‖ ^ 2

/-- The Green function `G_z(y) = ⟨δ_y, (H_v - z)⁻¹ δ_x⟩`. -/
def dimerGreen (a b : ℤ) (v : ℤ → ℝ) (z : ℂ) (x y : Set.Icc a b) : ℂ :=
  (cplx (schrH a b v) - z • (1 : Matrix (Set.Icc a b) (Set.Icc a b) ℂ))⁻¹ y x

/-- `c₁ = 1/((3 + 2λ) K₁ e^{K₂})`, `K₁ = κ (2 + 2λ)²`, `K₂ = κ (3 + 4λ)/2`. -/
def dimerC1 (lam : ℝ) : ℝ :=
  1 / ((3 + 2 * lam) * (dimerKappa lam * (2 + 2 * lam) ^ 2) *
    Real.exp (dimerKappa lam * (3 + 4 * lam) / 2))

/-- `c_{λ,p} = c₁² / (2^{p+3} π (p+1))`. -/
def dimerClp (lam p : ℝ) : ℝ := dimerC1 lam ^ 2 / (2 ^ (p + 3) * Real.pi * (p + 1))

/-- **Plancherel identity** (paper: eq. `eq:plancherel`). -/
theorem dimerMom_plancherel (a b : ℤ) (v : ℤ → ℝ) (p T : ℝ) (hp : 0 < p) (hT : 0 < T)
    (x : Set.Icc a b) :
    dimerMom a b v p T x =
      (1 / (2 * Real.pi * T)) * ∫ E : ℝ, ∑ y : Set.Icc a b, |((y.1 - x.1 : ℤ) : ℝ)| ^ p *
        ‖dimerGreen a b v (E + Complex.I / (2 * T)) x y‖ ^ 2 := by
  sorry

/-- **Every dimer configuration spreads at least diffusively** (paper: Proposition `prop:dimerlower`), on a
finite interval `[a, b]`. -/
theorem dimer_lower_bound (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (p : ℝ) (hp : 0 < p) (r : ℤ)
    (a b : ℤ) (v : ℤ → ℝ) (hv : IsDimerCfg lam r v a (b - a + 1).toNat) (x : Set.Icc a b) (T : ℝ)
    (hT : 1 ≤ T) (hxT : (a : ℝ) ≤ (x.1 : ℝ) - T ∧ (x.1 : ℝ) + T ≤ (b : ℝ)) :
    dimerClp lam p * T ^ (p - 1) ≤ dimerMom a b v p T x := by
  sorry

end DFL
