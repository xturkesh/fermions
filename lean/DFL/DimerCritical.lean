import Mathlib

/-!
# Critical energies of the random dimer sector law (paper: Section `sec:dimer`)

For the dimer sector law of `DFL.Dimer`, the one-body operator (in units of `J`) is the random
dimer model `-Δ + a·(±1 on both sites of each dimer)` with `a = 2h/J`. We verify the hypotheses of
Jitomirskaya–Schulz-Baldes–Stolz (CMP 233 (2003)), Definition 1 and Theorem 4, for `0 < a < 1`:
at the energy `E = a` the `(+a, +a)` dimer has transfer matrix `-1`, the `(-a, -a)` dimer has an
elliptic transfer matrix (`|tr| < 2`), the two commute (so `E = a` is a *critical energy*), and the
average of `e^{2iη}` over the two dimers (equal probabilities; `η₊ = π`, `cos η₋ = tr/2`) has
modulus `< 1`. The transport lower bounds of JSBS (infinite volume) are not formalized.

Convention (JSBS eq. (3) with hopping `t = 1`): the eigenvalue equation
`-ψ(n+1) + v(n) ψ(n) - ψ(n-1) = E ψ(n)` is transported by `(ψ(n), ψ(n-1)) ↦ (ψ(n+1), ψ(n))`,
i.e. by the matrix `!![v - E, -1; 1, 0]`.
-/

open Matrix

noncomputable section

namespace DFL

/-- One-step transfer matrix of `-ψ(n+1) + v ψ(n) - ψ(n-1) = E ψ(n)`. -/
def tm (v E : ℝ) : Matrix (Fin 2) (Fin 2) ℝ := !![v - E, -1; 1, 0]

/-- Transfer matrix across a dimer carrying the potential `v` on both of its sites. -/
def dimerTM (v E : ℝ) : Matrix (Fin 2) (Fin 2) ℝ := tm v E * tm v E

/-- Transfer matrices have determinant one. -/
theorem det_dimerTM (v E : ℝ) : (dimerTM v E).det = 1 := by
  have h : (tm v E).det = 1 := by
    rw [tm, Matrix.det_fin_two_of]; ring
  rw [dimerTM, Matrix.det_mul, h, one_mul]

/-- At `E = a` the `(+a, +a)` dimer is transparent up to a sign: its transfer matrix is `-1`. -/
theorem dimerTM_plus_crit (a : ℝ) : dimerTM a a = -1 := by
  rw [dimerTM, tm, sub_self, Matrix.mul_fin_two]
  ext i j
  fin_cases i <;> fin_cases j <;> simp

/-- At `E = a` the `(-a, -a)` dimer has trace `4a² - 2`. -/
theorem trace_dimerTM_minus_crit (a : ℝ) : (dimerTM (-a) a).trace = 4 * a ^ 2 - 2 := by
  rw [dimerTM, tm, Matrix.mul_fin_two, Matrix.trace_fin_two_of]
  ring

/-- For `0 < a < 1` the `(-a, -a)` dimer transfer matrix at `E = a` is elliptic. -/
theorem dimerTM_minus_elliptic {a : ℝ} (h0 : 0 < a) (h1 : a < 1) :
    |(dimerTM (-a) a).trace| < 2 := by
  rw [trace_dimerTM_minus_crit, abs_lt]
  constructor <;> nlinarith

/-- The two dimer transfer matrices commute at `E = a`; with ellipticity this makes `E = a` a
critical energy in the sense of JSBS, Definition 1. -/
theorem dimerTM_commute_crit (a : ℝ) : Commute (dimerTM a a) (dimerTM (-a) a) := by
  rw [dimerTM_plus_crit]
  exact (Commute.one_left _).neg_left

/-- The hypothesis of JSBS, Theorem 4, for the dimer law with equal probabilities:
`|½ e^{2iη₊} + ½ e^{2iη₋}| < 1` with `η₊ = π` and `η₋ = arccos(tr/2) = arccos(2a² - 1)`. -/
theorem jsbs_condition {a : ℝ} (h0 : 0 < a) (h1 : a < 1) :
    ‖(1 / 2 : ℂ) * Complex.exp (2 * (Real.pi : ℂ) * Complex.I) +
        (1 / 2 : ℂ) * Complex.exp (2 * (Real.arccos (2 * a ^ 2 - 1) : ℂ) * Complex.I)‖ < 1 := by
  set c : ℝ := 2 * a ^ 2 - 1 with hc
  have hc1 : -1 < c := by nlinarith
  have hc2 : c < 1 := by nlinarith
  set θ : ℝ := 2 * Real.arccos c with hθ
  -- `cos θ = cos (2 arccos c) = 2c² - 1`
  have hcosθ : Real.cos θ = 2 * c ^ 2 - 1 := by
    rw [hθ, Real.cos_two_mul, Real.cos_arccos hc1.le hc2.le]
  have harg : 2 * (Real.arccos c : ℂ) * Complex.I = (θ : ℂ) * Complex.I := by
    rw [hθ]; push_cast; ring
  -- the average is `(1 + e^{iθ}) / 2 = (1 + cos θ)/2 + i sin θ / 2`
  have hz : (1 / 2 : ℂ) * Complex.exp (2 * (Real.pi : ℂ) * Complex.I) +
      (1 / 2 : ℂ) * Complex.exp (2 * (Real.arccos c : ℂ) * Complex.I) =
      (((1 + Real.cos θ) / 2 : ℝ) : ℂ) + ((Real.sin θ / 2 : ℝ) : ℂ) * Complex.I := by
    rw [Complex.exp_two_pi_mul_I, harg, Complex.exp_mul_I, ← Complex.ofReal_cos,
      ← Complex.ofReal_sin]
    push_cast
    ring
  -- its squared modulus is `(1 + cos θ)/2 = c²`
  have hsq : ‖(((1 + Real.cos θ) / 2 : ℝ) : ℂ) + ((Real.sin θ / 2 : ℝ) : ℂ) * Complex.I‖ ^ 2 =
      c ^ 2 := by
    rw [Complex.sq_norm, Complex.normSq_add_mul_I]
    linear_combination (1 / 4 : ℝ) * Real.sin_sq_add_cos_sq θ + (1 / 2 : ℝ) * hcosθ
  rw [hz, ← sq_lt_one_iff₀ (norm_nonneg _), hsq]
  nlinarith

/-- For `a ≥ 1` the energy `E = a` is not critical: the `(-a,-a)` dimer is not elliptic. -/
theorem dimerTM_minus_not_elliptic {a : ℝ} (h1 : 1 ≤ a) : 2 ≤ |(dimerTM (-a) a).trace| := by
  rw [trace_dimerTM_minus_crit]
  exact le_abs.mpr (Or.inl (by nlinarith))

end DFL
