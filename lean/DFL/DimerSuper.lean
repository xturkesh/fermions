import DFL.DimerSecond
import DFL.DimerState

/-!
# Superdiffusion of the random dimer model (paper: Proposition `prop:dimersuper`, Theorem `thm:dimer`(c),(d))

For the random dimer potential `v = v^{(r)}` on `ℤ` (values `±λ`, `λ = 2|h|/J ∈ (0, 1)`):

* a deterministic bound (paper: eq. `eq:superdet`): with `z = E + i/(2T)`, `L = ⌊T⌋` and
  `X_j(E) = ‖T^z(j+1, 1)‖₂² + ‖T^z(0, -j)‖₂²`,
  `M^ℤ_{v,p}(T; 0) ≥ c₀²/(4πT) ∫_{λ-δ}^{λ+δ} ∑_{j ≤ L} j^p / X_j(E) dE` (Plancherel `dimerMom_plancherel` and the
  Green-function propagation of `DimerLower.lean` on `[-L', L']`, with `‖A⁻¹‖ ≤ ‖A‖` in `SL(2, ℂ)`, and the limit
  `L' → ∞`, `tendsto_dimerMom_Z`);
* `dimer_super_annealed` (paper: `prop:dimersuper`(a)): `E M^ℤ_{v,p}(T; 0) ≥ c T^{p - 1/2}` for `T ≥ T₀`, by
  `δ = T^{-1/2}`, `lem:dimersecond` (`integral_hs2_tmProd_le`) and Jensen's inequality `E[1/X] ≥ 1/E[X]`;
* `dimer_super_ae` (paper: `prop:dimersuper`(b)): almost surely `M^ℤ_{v,p}(T; 0) ≥ T^{p - 1/2 - α}` eventually, for
  every `α > 0` (Cauchy–Schwarz, Markov, Borel–Cantelli along `T = 2^k`);
* `dimerMInf_super`, `dimerMInf_exponent` (paper: `thm:dimer`(c)): `M^{∞,r}_p(T) ≥ c (JT)^{p - 1/2}`;
* `dimerD_unbounded_half`, `dimer_LGT_unbounded_half`, `not_EC_dimer_half` (paper: `thm:dimer`(d) for every
  `p > 1/2`).
-/

open MeasureTheory Filter Topology Matrix Finset

noncomputable section

namespace DFL

/-- **Superdiffusion on average** (paper: Proposition `prop:dimersuper`(a)):
`E M^ℤ_{v^{(r)},p}(T; 0) ≥ c T^{p - 1/2}` for `T ≥ T₀`. -/
theorem dimer_super_annealed (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (r : ℤ)
    (p : ℝ) (hp : 0 < p) :
    ∃ c T₀ : ℝ, 0 < c ∧ 0 < T₀ ∧ ∀ T : ℝ, T₀ ≤ T →
      c * T ^ (p - 1 / 2) ≤ ∫ η, dimerMomZ (dimerPotZ J h r η) p T 0 ∂muZ := by
  sorry

/-- **Superdiffusion almost surely** (paper: Proposition `prop:dimersuper`(b)): almost surely, for every
`α > 0`, `M^ℤ_{v^{(r)},p}(T; 0) ≥ T^{p - 1/2 - α}` for all large `T`. -/
theorem dimer_super_ae (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (r : ℤ) (p : ℝ)
    (hp : 0 < p) :
    ∀ᵐ η ∂muZ, ∀ α : ℝ, 0 < α → ∀ᶠ T in atTop,
      T ^ (p - 1 / 2 - α) ≤ dimerMomZ (dimerPotZ J h r η) p T 0 := by
  sorry

/-- **Theorem `thm:dimer`(c)**: `M^{∞,r}_p(T) ≥ c (JT)^{p - 1/2}` for `JT ≥ T₀`. -/
theorem dimerMInf_super (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (r : ℤ) (p : ℝ)
    (hp : 0 < p) :
    ∃ c T₀ : ℝ, 0 < c ∧ 0 < T₀ ∧ ∀ T : ℝ, T₀ ≤ J * T →
      c * (J * T) ^ (p - 1 / 2) ≤ dimerMInf J h r p T := by
  sorry

/-- **Theorem `thm:dimer`(c)**, the annealed exponent: for every `α > 0`, `M^{∞,r}_p(T) ≥ T^{p - 1/2 - α}` for all
large `T`. -/
theorem dimerMInf_exponent (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (r : ℤ) (p : ℝ)
    (hp : 0 < p) (α : ℝ) (hα : 0 < α) :
    ∀ᶠ T in atTop, T ^ (p - 1 / 2 - α) ≤ dimerMInf J h r p T := by
  sorry

/-- **Theorem `thm:dimer`(d)** for every `p > 1/2`: `sup_{N even, x, t} D^{(N)}_p(t; x) = ∞`. -/
theorem dimerD_unbounded_half (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (p : ℝ)
    (hp : 1 / 2 < p) :
    ¬ BddAbove (Set.range (fun r : (Σ M : ℕ, Fin (2 * M) × ℝ) =>
      dimerD J h (2 * r.1) p r.2.2 r.2.1)) := by
  sorry

/-- **Theorem `thm:dimer`(d)** for every `p > 1/2`, for the single-fermion dimer state. -/
theorem dimer_LGT_unbounded_half (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (p : ℝ)
    (hp : 1 / 2 < p) :
    ¬ BddAbove (Set.range (fun r : (Σ M : ℕ, Fin (2 * M) × ℝ) =>
      dimerDLGT J h (2 * r.1) p r.2.2 r.2.1)) := by
  sorry

/-- **Hypothesis (EC) fails for the dimer law** (paper: `thm:dimer`(d)), for every `p > 1/2`: there is no `F ≥ 0`
with `∑_d (d+1)^p F(d) < ∞` and `E_dim Q^{[1,N]}_q(x, y) ≤ F(|x - y|)` for all even `N`. -/
theorem not_EC_dimer_half (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (p : ℝ)
    (hp : 1 / 2 < p) (F : ℕ → ℝ) (hF0 : ∀ d, 0 ≤ F d)
    (hFs : Summable (fun d : ℕ => ((d : ℝ) + 1) ^ p * F d)) :
    ¬ ∀ (N : ℕ), Even N → ∀ x y : Fin N,
      dimerAvg (fun q : Cfg N => corrQ (hq J h q) (hq_isHermitian J h q) x y) ≤
        F (Nat.dist x.val y.val) := by
  sorry

end DFL
