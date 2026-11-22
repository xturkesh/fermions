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

/-! ### Auxiliary lemmas

*Plancherel.* `H_v` is diagonalised by Mathlib's spectral theorem, `H_v = W diag(μ) W*`. Both sides of
`eq:plancherel` then become `∑_y |y - x|^p ∑_{α,β} c_α(y) conj(c_β(y)) / (1 + iT(μ_α - μ_β))`, with
`c_α(y) = W_{yα} conj(W_{xα})`: on the time side by `∫_0^∞ T⁻¹ e^{-t/T} e^{-itd} dt = 1/(1 + iTd)`, on the
energy side by `∫_ℝ ((a - E) - iε)⁻¹ ((b - E) + iε)⁻¹ dE = 2πi/((b - a) + 2iε)` (explicit antiderivative
built from `log` and `arctan`).

*Lower bound.* Steps 2–4 of the paper: the Green function extended by zero solves the three-term recursion
away from `x`, so it is propagated by the transfer matrices `tmProd`; these have determinant one and
are bounded by `K₁ e^{K₂}` on windows of length `≤ T` (Lemma `lem:dimertransfer`).
-/

namespace DimerLowerAux

open Filter Topology

/-- `((a - E)² + ε²)⁻¹` is integrable in `E`. -/
lemma integrable_lorentz (a ε : ℝ) (hε : ε ≠ 0) :
    Integrable (fun E : ℝ => ((a - E) ^ 2 + ε ^ 2)⁻¹) := by
  have h := ((integrable_inv_one_add_sq.comp_div hε).comp_sub_right a).const_mul (ε ^ 2)⁻¹
  refine h.congr (ae_of_all _ fun E => ?_)
  have hε2 : ε ^ 2 ≠ 0 := by positivity
  have h3 : (a - E) ^ 2 + ε ^ 2 ≠ 0 := by positivity
  simp only
  field_simp
  ring

/-- A complex number with nonzero imaginary part is nonzero. -/
lemma ne_zero_of_im_ne {w : ℂ} (h : w.im ≠ 0) : w ≠ 0 := fun h' => h (by simp [h'])

/-- A complex number with nonzero real part is nonzero. -/
lemma ne_zero_of_re_ne {w : ℂ} (h : w.re ≠ 0) : w ≠ 0 := fun h' => h (by simp [h'])

/-- `‖A - ε i‖² = A² + ε²`. -/
lemma norm_sq_sub_mul_I (A ε : ℝ) : ‖(A : ℂ) - ε * Complex.I‖ ^ 2 = A ^ 2 + ε ^ 2 := by
  rw [Complex.sq_norm, Complex.normSq_apply]
  simp
  ring

/-- `(A - ε i)⁻¹ = A/(A² + ε²) + i ε/(A² + ε²)`. -/
lemma inv_sub_mul_I (A ε : ℝ) (hε : ε ≠ 0) :
    ((A : ℂ) - ε * Complex.I)⁻¹ =
      ((A / (A ^ 2 + ε ^ 2) : ℝ) : ℂ) + ((ε / (A ^ 2 + ε ^ 2) : ℝ) : ℂ) * Complex.I := by
  have hden : A ^ 2 + ε ^ 2 ≠ 0 := by positivity
  have hne : (A : ℂ) - ε * Complex.I ≠ 0 := ne_zero_of_im_ne (by simpa using hε)
  have hden' : ((A ^ 2 + ε ^ 2 : ℝ) : ℂ) ≠ 0 := by exact_mod_cast hden
  rw [show ((A / (A ^ 2 + ε ^ 2) : ℝ) : ℂ) + ((ε / (A ^ 2 + ε ^ 2) : ℝ) : ℂ) * Complex.I =
      ((A : ℂ) + ε * Complex.I) / ((A ^ 2 + ε ^ 2 : ℝ) : ℂ) by push_cast; ring]
  rw [eq_div_iff hden', inv_mul_eq_iff_eq_mul₀ hne]
  push_cast
  linear_combination (ε : ℂ) ^ 2 * Complex.I_sq

/-- The resolvent product `((a - E) - εi)⁻¹ ((b - E) + εi)⁻¹` is integrable in `E`. -/
lemma integrable_resolvent_prod (a b ε : ℝ) (hε : 0 < ε) :
    Integrable (fun E : ℝ =>
      (((a - E : ℝ) : ℂ) - ε * Complex.I)⁻¹ * (((b - E : ℝ) : ℂ) + ε * Complex.I)⁻¹) := by
  have hg := ((integrable_lorentz a ε hε.ne').add (integrable_lorentz b ε hε.ne')).div_const 2
  have hX : ∀ E : ℝ, (((a - E : ℝ) : ℂ) - ε * Complex.I) ≠ 0 := fun E =>
    ne_zero_of_im_ne (by simp; exact hε.ne')
  have hY : ∀ E : ℝ, (((b - E : ℝ) : ℂ) + ε * Complex.I) ≠ 0 := fun E =>
    ne_zero_of_im_ne (by simp; exact hε.ne')
  refine hg.mono' ?_ (ae_of_all _ fun E => ?_)
  · apply Continuous.aestronglyMeasurable
    refine Continuous.mul (Continuous.inv₀ (by fun_prop) hX) (Continuous.inv₀ (by fun_prop) hY)
  · have h1 := norm_sq_sub_mul_I (a - E) ε
    have h2 : ‖((b - E : ℝ) : ℂ) + ε * Complex.I‖ ^ 2 = (b - E) ^ 2 + ε ^ 2 := by
      have := norm_sq_sub_mul_I (b - E) (-ε)
      simpa using this
    set u := ‖((a - E : ℝ) : ℂ) - ε * Complex.I‖ with hu
    set w := ‖((b - E : ℝ) : ℂ) + ε * Complex.I‖ with hw
    have hu0 : 0 < u := norm_pos_iff.mpr (hX E)
    have hw0 : 0 < w := norm_pos_iff.mpr (hY E)
    rw [norm_mul, norm_inv, norm_inv]
    show u⁻¹ * w⁻¹ ≤ (((a - E) ^ 2 + ε ^ 2)⁻¹ + ((b - E) ^ 2 + ε ^ 2)⁻¹) / 2
    rw [← h1, ← h2, ← inv_pow, ← inv_pow]
    nlinarith [sq_nonneg (u⁻¹ - w⁻¹)]

/-- An antiderivative of `((c - E) - εi)⁻¹`. -/
def antider (c ε E : ℝ) : ℂ :=
  -((Real.log ((c - E) ^ 2 + ε ^ 2) / 2 : ℝ) : ℂ) - Complex.I * (Real.arctan ((c - E) / ε) : ℂ)

/-- The derivative of `antider`. -/
lemma hasDerivAt_antider (c ε : ℝ) (hε : ε ≠ 0) (E : ℝ) :
    HasDerivAt (antider c ε) ((((c - E : ℝ) : ℂ) - ε * Complex.I)⁻¹) E := by
  have h0 : HasDerivAt (fun E : ℝ => c - E) (-1) E := by
    simpa using (hasDerivAt_id E).const_sub c
  have h1 : HasDerivAt (fun E : ℝ => (c - E) ^ 2 + ε ^ 2) (2 * (c - E) * (-1)) E := by
    have := (h0.pow 2).add_const (ε ^ 2)
    convert this using 1
    simp
  have hpos : (c - E) ^ 2 + ε ^ 2 ≠ 0 := by positivity
  have hlog := (h1.log hpos).div_const 2
  have harc := (h0.div_const ε).arctan
  have h := (hlog.ofReal_comp.neg).sub (harc.ofReal_comp.const_mul Complex.I)
  convert h using 1
  · funext y
    rfl
  rw [inv_sub_mul_I _ _ hε]
  have e1 : 2 * (c - E) * (-1) / ((c - E) ^ 2 + ε ^ 2) / 2 = -((c - E) / ((c - E) ^ 2 + ε ^ 2)) := by
    field_simp
  have e2 : 1 / (1 + ((c - E) / ε) ^ 2) * (-1 / ε) = -(ε / ((c - E) ^ 2 + ε ^ 2)) := by
    field_simp
    ring
  rw [e1, e2]
  push_cast
  ring

/-- `log((b - E)² + c) - log((a - E)² + d) → 0` as `E → ±∞`. -/
lemma tendsto_log_sub_log (a b c d : ℝ) (hc : 0 < c) (hd : 0 < d) {l : Filter ℝ}
    (hl : Tendsto (fun E : ℝ => E⁻¹) l (𝓝 0)) (hl' : ∀ᶠ E in l, E ≠ 0) :
    Tendsto (fun E => Real.log ((b - E) ^ 2 + c) - Real.log ((a - E) ^ 2 + d)) l (𝓝 0) := by
  have hgc : ContinuousAt
      (fun u : ℝ => ((b * u - 1) ^ 2 + c * u ^ 2) / ((a * u - 1) ^ 2 + d * u ^ 2)) 0 :=
    ContinuousAt.div₀ (by fun_prop) (by fun_prop) (by norm_num)
  have hgt := hgc.tendsto
  have hg0 : ((b * 0 - 1) ^ 2 + c * 0 ^ 2) / ((a * 0 - 1) ^ 2 + d * 0 ^ 2) = (1 : ℝ) := by norm_num
  rw [hg0] at hgt
  have h1 : Tendsto (fun E : ℝ => Real.log
      (((b * E⁻¹ - 1) ^ 2 + c * E⁻¹ ^ 2) / ((a * E⁻¹ - 1) ^ 2 + d * E⁻¹ ^ 2))) l (𝓝 0) := by
    have := ((Real.continuousAt_log (by norm_num : (1:ℝ) ≠ 0)).tendsto.comp hgt).comp hl
    rw [Real.log_one] at this
    exact this
  refine h1.congr' ?_
  filter_upwards [hl'] with E hE
  have hb : (b - E) ^ 2 + c ≠ 0 := by positivity
  have ha : (a - E) ^ 2 + d ≠ 0 := by positivity
  rw [← Real.log_div hb ha]
  congr 1
  field_simp

/-- `arctan((c - E)/ε) → -π/2` as `E → +∞`. -/
lemma tendsto_arctan_atTop' (c ε : ℝ) (hε : 0 < ε) :
    Tendsto (fun E : ℝ => Real.arctan ((c - E) / ε)) atTop (𝓝 (-(Real.pi / 2))) := by
  have h : Tendsto (fun E : ℝ => (c - E) / ε) atTop atBot := by
    have := (tendsto_atBot_add_const_left atTop c tendsto_neg_atTop_atBot).atBot_div_const hε
    simpa [sub_eq_add_neg] using this
  exact (Real.tendsto_arctan_atBot.mono_right nhdsWithin_le_nhds).comp h

/-- `arctan((c - E)/ε) → π/2` as `E → -∞`. -/
lemma tendsto_arctan_atBot' (c ε : ℝ) (hε : 0 < ε) :
    Tendsto (fun E : ℝ => Real.arctan ((c - E) / ε)) atBot (𝓝 (Real.pi / 2)) := by
  have h : Tendsto (fun E : ℝ => (c - E) / ε) atBot atTop := by
    have := (tendsto_atTop_add_const_left atBot c tendsto_neg_atBot_atTop).atTop_div_const hε
    simpa [sub_eq_add_neg] using this
  exact (Real.tendsto_arctan_atTop.mono_right nhdsWithin_le_nhds).comp h

/-- **The resolvent-product integral**:
`∫ ((a - E) - εi)⁻¹ ((b - E) + εi)⁻¹ dE = 2πi/((b - a) + 2εi)`. -/
lemma integral_resolvent_prod (a b ε : ℝ) (hε : 0 < ε) :
    ∫ E : ℝ, (((a - E : ℝ) : ℂ) - ε * Complex.I)⁻¹ * (((b - E : ℝ) : ℂ) + ε * Complex.I)⁻¹ =
      2 * Real.pi * Complex.I / (((b - a : ℝ) : ℂ) + 2 * ε * Complex.I) := by
  set D : ℂ := ((b - a : ℝ) : ℂ) + 2 * ε * Complex.I with hD
  have hDne : D ≠ 0 := ne_zero_of_im_ne (by simp [hD]; exact hε.ne')
  have hX : ∀ E : ℝ, (((a - E : ℝ) : ℂ) - ε * Complex.I) ≠ 0 := fun E =>
    ne_zero_of_im_ne (by simp; exact hε.ne')
  have hY : ∀ E : ℝ, (((b - E : ℝ) : ℂ) + ε * Complex.I) ≠ 0 := fun E =>
    ne_zero_of_im_ne (by simp; exact hε.ne')
  have hderiv : ∀ E : ℝ, HasDerivAt (fun E => D⁻¹ * (antider a ε E - antider b (-ε) E))
      ((((a - E : ℝ) : ℂ) - ε * Complex.I)⁻¹ * (((b - E : ℝ) : ℂ) + ε * Complex.I)⁻¹) E := by
    intro E
    have h := ((hasDerivAt_antider a ε hε.ne' E).sub
      (hasDerivAt_antider b (-ε) (neg_ne_zero.mpr hε.ne') E)).const_mul D⁻¹
    convert h using 1
    have hXE := hX E
    have hYE := hY E
    rw [show (((b - E : ℝ) : ℂ) - ((-ε : ℝ) : ℂ) * Complex.I) =
      (((b - E : ℝ) : ℂ) + ε * Complex.I) by push_cast; ring, inv_sub_inv hXE hYE]
    rw [show (((b - E : ℝ) : ℂ) + ε * Complex.I) - (((a - E : ℝ) : ℂ) - ε * Complex.I) = D by
      rw [hD]; push_cast; ring]
    field_simp
  have hdiff : ∀ E : ℝ, antider a ε E - antider b (-ε) E =
      (((Real.log ((b - E) ^ 2 + (-ε) ^ 2) - Real.log ((a - E) ^ 2 + ε ^ 2)) / 2 : ℝ) : ℂ) -
        Complex.I * ((Real.arctan ((a - E) / ε) + Real.arctan ((b - E) / ε) : ℝ) : ℂ) := by
    intro E
    simp only [antider, div_neg, Real.arctan_neg]
    push_cast
    ring
  have hlog_top := tendsto_log_sub_log a b ((-ε) ^ 2) (ε ^ 2) (by rw [neg_sq]; positivity) (by positivity)
    tendsto_inv_atTop_zero (eventually_ne_atTop 0)
  have hlog_bot := tendsto_log_sub_log a b ((-ε) ^ 2) (ε ^ 2) (by rw [neg_sq]; positivity) (by positivity)
    tendsto_inv_atBot_zero (eventually_ne_atBot 0)
  have htop : Tendsto (fun E => D⁻¹ * (antider a ε E - antider b (-ε) E)) atTop
      (𝓝 (D⁻¹ * (((0 / 2 : ℝ) : ℂ) -
        Complex.I * ((-(Real.pi / 2) + -(Real.pi / 2) : ℝ) : ℂ)))) := by
    simp_rw [hdiff]
    exact ((((hlog_top.div_const 2).ofReal).sub (((tendsto_arctan_atTop' a ε hε).add
      (tendsto_arctan_atTop' b ε hε)).ofReal.const_mul Complex.I))).const_mul D⁻¹
  have hbot : Tendsto (fun E => D⁻¹ * (antider a ε E - antider b (-ε) E)) atBot
      (𝓝 (D⁻¹ * (((0 / 2 : ℝ) : ℂ) - Complex.I * ((Real.pi / 2 + Real.pi / 2 : ℝ) : ℂ)))) := by
    simp_rw [hdiff]
    exact ((((hlog_bot.div_const 2).ofReal).sub (((tendsto_arctan_atBot' a ε hε).add
      (tendsto_arctan_atBot' b ε hε)).ofReal.const_mul Complex.I))).const_mul D⁻¹
  rw [integral_of_hasDerivAt_of_tendsto hderiv (integrable_resolvent_prod a b ε hε) hbot htop]
  push_cast
  field_simp
  ring

/-- `∫_0^∞ (1/T) e^{-t/T} e^{-itd} dt = 1/(1 + iTd)`. -/
lemma integral_exp_time (T d : ℝ) (hT : 0 < T) :
    ∫ t in Set.Ioi (0 : ℝ), ((1 / T : ℝ) : ℂ) * Complex.exp ((-(1 / T : ℝ) - d * Complex.I) * t) =
      1 / (1 + T * d * Complex.I) := by
  have hre : (-(1 / T : ℝ) - d * Complex.I : ℂ).re < 0 := by simp; positivity
  have hne : (-(1 / T : ℝ) - d * Complex.I : ℂ) ≠ 0 := fun h => by
    rw [h, Complex.zero_re] at hre; exact lt_irrefl _ hre
  have hne2 : (1 + T * d * Complex.I : ℂ) ≠ 0 := ne_zero_of_re_ne (by simp)
  rw [integral_const_mul, integral_exp_mul_complex_Ioi hre 0]
  have hT' : (T : ℂ) ≠ 0 := by exact_mod_cast hT.ne'
  have hne3 : (-1 - (T : ℂ) * d * Complex.I) ≠ 0 := ne_zero_of_re_ne (by simp)
  push_cast
  simp only [mul_zero, Complex.exp_zero]
  field_simp
  ring


section Spectral

variable {n : Type*} [Fintype n] [DecidableEq n]

omit [Fintype n] [DecidableEq n] in
/-- The complexification of a real symmetric matrix is Hermitian. -/
lemma cplx_isHermitian (A : Matrix n n ℝ) (hA : A.IsSymm) : (cplx A).IsHermitian := by
  ext i j
  simp only [Matrix.conjTranspose_apply, cplx, Matrix.map_apply, Complex.star_def,
    Complex.conj_ofReal]
  rw [hA.apply]

/-- Entries of `W D W*`. -/
lemma conj_diagonal_apply (W : Matrix n n ℂ) (d : n → ℂ) (i j : n) :
    (W * diagonal d * star W) i j = ∑ α, (W i α * star (W j α)) * d α := by
  rw [Matrix.mul_apply]
  refine Finset.sum_congr rfl fun α _ => ?_
  rw [Matrix.mul_diagonal, Matrix.star_apply]
  ring

/-- The spectral theorem for the complexification of a real matrix. -/
lemma spectral_cplx (A : Matrix n n ℝ) (hA : (cplx A).IsHermitian) :
    cplx A = (hA.eigenvectorUnitary : Matrix n n ℂ) *
      diagonal (fun α => (hA.eigenvalues α : ℂ)) * star (hA.eigenvectorUnitary : Matrix n n ℂ) := by
  have h := hA.spectral_theorem
  rw [Unitary.conjStarAlgAut_apply] at h
  exact h

/-- The propagator in the eigenbasis. -/
lemma evo_eq_spectral (A : Matrix n n ℝ) (hA : (cplx A).IsHermitian) (t : ℝ) :
    evo A t = (hA.eigenvectorUnitary : Matrix n n ℂ) *
      diagonal (fun α => Complex.exp (-(t : ℂ) * Complex.I * (hA.eigenvalues α : ℂ))) *
      star (hA.eigenvectorUnitary : Matrix n n ℂ) := by
  set W : Matrix n n ℂ := ↑hA.eigenvectorUnitary with hW
  have hWu : IsUnit W :=
    ⟨⟨W, star W, Unitary.coe_mul_star_self _, Unitary.coe_star_mul_self _⟩, rfl⟩
  have hinv : W⁻¹ = star W := Matrix.inv_eq_right_inv (Unitary.coe_mul_star_self _)
  unfold evo
  conv_lhs => rw [spectral_cplx A hA, ← hW]
  rw [← hinv]
  rw [show ((-(t : ℂ) * Complex.I) • (W * diagonal (fun α => (hA.eigenvalues α : ℂ)) * W⁻¹)) =
      W * diagonal (fun α => (-(t : ℂ) * Complex.I) * (hA.eigenvalues α : ℂ)) * W⁻¹ by
    rw [show (fun α => (-(t : ℂ) * Complex.I) * (hA.eigenvalues α : ℂ)) =
      (-(t : ℂ) * Complex.I) • (fun α => (hA.eigenvalues α : ℂ)) from rfl, Matrix.diagonal_smul,
      Matrix.mul_smul, Matrix.smul_mul]]
  have hv : NormedSpace.exp (fun α => -(t : ℂ) * Complex.I * (hA.eigenvalues α : ℂ)) =
      fun α => Complex.exp (-(t : ℂ) * Complex.I * (hA.eigenvalues α : ℂ)) := by
    funext α
    rw [Pi.coe_exp, Complex.exp_eq_exp_ℂ]
  rw [Matrix.exp_conj _ _ hWu, Matrix.exp_diagonal, hv]

/-- The resolvent in the eigenbasis: `(A - z) W (μ - z)⁻¹ W* = 1`. -/
lemma resolvent_mul_spectral (A : Matrix n n ℝ) (hA : (cplx A).IsHermitian) (z : ℂ)
    (hz : z.im ≠ 0) :
    (cplx A - z • (1 : Matrix n n ℂ)) * ((hA.eigenvectorUnitary : Matrix n n ℂ) *
      diagonal (fun α => ((hA.eigenvalues α : ℂ) - z)⁻¹) *
        star (hA.eigenvectorUnitary : Matrix n n ℂ)) = 1 := by
  set W : Matrix n n ℂ := ↑hA.eigenvectorUnitary with hW
  have h1 : W * star W = 1 := Unitary.coe_mul_star_self _
  have h2 : star W * W = 1 := Unitary.coe_star_mul_self _
  have hz1 : z • (1 : Matrix n n ℂ) = W * diagonal (fun _ => z) * star W := by
    rw [← smul_one_eq_diagonal, Matrix.mul_smul, Matrix.mul_one, Matrix.smul_mul, h1]
  have hsub : cplx A - z • (1 : Matrix n n ℂ) =
      W * diagonal (fun α => (hA.eigenvalues α : ℂ) - z) * star W := by
    conv_lhs => rw [spectral_cplx A hA, ← hW, hz1]
    rw [← Matrix.sub_mul, ← Matrix.mul_sub, diagonal_sub]
  have hne : ∀ α, (hA.eigenvalues α : ℂ) - z ≠ 0 := fun α h => hz (by
    have := congrArg Complex.im h
    simpa using this.symm)
  rw [hsub]
  calc W * diagonal (fun α => (hA.eigenvalues α : ℂ) - z) * star W *
        (W * diagonal (fun α => ((hA.eigenvalues α : ℂ) - z)⁻¹) * star W)
      = W * (diagonal (fun α => (hA.eigenvalues α : ℂ) - z) * (star W * W) *
          diagonal (fun α => ((hA.eigenvalues α : ℂ) - z)⁻¹)) * star W := by
        simp only [Matrix.mul_assoc]
    _ = W * 1 * star W := by
        rw [h2, Matrix.mul_one, diagonal_mul_diagonal]
        congr 2
        rw [← diagonal_one]
        congr 1
        funext α
        exact mul_inv_cancel₀ (hne α)
    _ = 1 := by rw [Matrix.mul_one, h1]

/-- The resolvent in the eigenbasis. -/
lemma resolvent_eq_spectral (A : Matrix n n ℝ) (hA : (cplx A).IsHermitian) (z : ℂ)
    (hz : z.im ≠ 0) :
    (cplx A - z • (1 : Matrix n n ℂ))⁻¹ = (hA.eigenvectorUnitary : Matrix n n ℂ) *
      diagonal (fun α => ((hA.eigenvalues α : ℂ) - z)⁻¹) *
        star (hA.eigenvectorUnitary : Matrix n n ℂ) :=
  Matrix.inv_eq_right_inv (resolvent_mul_spectral A hA z hz)

/-- `‖∑ f‖² = ∑_{α,β} f_α conj(f_β)` (as complex numbers). -/
lemma ofReal_norm_sq_sum {ι : Type*} (s : Finset ι) (f : ι → ℂ) :
    ((‖∑ α ∈ s, f α‖ ^ 2 : ℝ) : ℂ) = ∑ α ∈ s, ∑ β ∈ s, f α * star (f β) := by
  rw [Complex.ofReal_pow, ← Complex.mul_conj', map_sum, Finset.sum_mul_sum]
  rfl

/-- The coefficient `W_{yα} conj(W_{xα})` of the spectral expansion of `⟨δ_y, f(A) δ_x⟩`. -/
def coef {A : Matrix n n ℂ} (hA : A.IsHermitian) (x y α : n) : ℂ :=
  (hA.eigenvectorUnitary : Matrix n n ℂ) y α * star ((hA.eigenvectorUnitary : Matrix n n ℂ) x α)

/-- Entries of the propagator in the eigenbasis. -/
lemma evo_apply_spectral (A : Matrix n n ℝ) (hA : (cplx A).IsHermitian) (t : ℝ) (x y : n) :
    evo A t y x =
      ∑ α, coef hA x y α * Complex.exp (-(t : ℂ) * Complex.I * (hA.eigenvalues α : ℂ)) := by
  rw [evo_eq_spectral A hA t, conj_diagonal_apply]
  rfl

/-- Entries of the resolvent in the eigenbasis. -/
lemma resolvent_apply_spectral (A : Matrix n n ℝ) (hA : (cplx A).IsHermitian) (z : ℂ)
    (hz : z.im ≠ 0) (x y : n) :
    (cplx A - z • (1 : Matrix n n ℂ))⁻¹ y x =
      ∑ α, coef hA x y α * ((hA.eigenvalues α : ℂ) - z)⁻¹ := by
  rw [resolvent_eq_spectral A hA z hz, conj_diagonal_apply]
  rfl

/-- `e^{-t/T} e^{-ita} conj(e^{-itb}) = e^{(-1/T - i(a - b)) t}`. -/
lemma exp_mul_star_exp (t T a b : ℝ) :
    (Real.exp (-t / T) : ℂ) * (Complex.exp (-(t : ℂ) * Complex.I * a) *
      star (Complex.exp (-(t : ℂ) * Complex.I * b))) =
      Complex.exp ((-(1 / T : ℝ) - ((a - b : ℝ) : ℂ) * Complex.I) * t) := by
  rw [Complex.star_def, ← Complex.exp_conj, Complex.ofReal_exp, ← Complex.exp_add,
    ← Complex.exp_add]
  congr 1
  simp only [map_mul, map_neg, Complex.conj_ofReal, Complex.conj_I]
  push_cast
  ring

/-- One term of the time integrand. -/
lemma term_time (T t a b wy : ℝ) (ca cb : ℂ) :
    ((1 / T : ℝ) : ℂ) * (Real.exp (-t / T) : ℂ) * ((wy : ℂ) *
      ((ca * Complex.exp (-(t : ℂ) * Complex.I * a)) *
        star (cb * Complex.exp (-(t : ℂ) * Complex.I * b)))) =
    ((wy : ℂ) * (ca * star cb)) * (((1 / T : ℝ) : ℂ) *
      Complex.exp ((-(1 / T : ℝ) - ((a - b : ℝ) : ℂ) * Complex.I) * t)) := by
  rw [← exp_mul_star_exp t T a b, star_mul']
  ring

/-- The time integrand, expanded in the eigenbasis. -/
lemma time_integrand_eq (A : Matrix n n ℝ) (hA : (cplx A).IsHermitian) (w : n → ℝ) (T t : ℝ)
    (x : n) :
    ((1 / T * Real.exp (-t / T) * ∑ y, w y * ‖evo A t y x‖ ^ 2 : ℝ) : ℂ) =
      ∑ y, ∑ α, ∑ β, ((w y : ℂ) * (coef hA x y α * star (coef hA x y β))) *
        (((1 / T : ℝ) : ℂ) * Complex.exp ((-(1 / T : ℝ) -
          ((hA.eigenvalues α - hA.eigenvalues β : ℝ) : ℂ) * Complex.I) * t)) := by
  simp only [Complex.ofReal_mul, Complex.ofReal_sum]
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun y _ => ?_
  rw [evo_apply_spectral A hA t x y, ofReal_norm_sq_sum, Finset.mul_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun α _ => ?_
  rw [Finset.mul_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun β _ => ?_
  exact term_time T t _ _ (w y) _ _

/-- The energy integrand, expanded in the eigenbasis. -/
lemma green_integrand_eq (A : Matrix n n ℝ) (hA : (cplx A).IsHermitian) (w : n → ℝ) (T : ℝ)
    (hT : 0 < T) (E : ℝ) (x : n) :
    ((∑ y, w y * ‖(cplx A - ((E : ℂ) + Complex.I / (2 * T)) • (1 : Matrix n n ℂ))⁻¹ y x‖ ^ 2 :
      ℝ) : ℂ) =
      ∑ y, ∑ α, ∑ β, ((w y : ℂ) * (coef hA x y α * star (coef hA x y β))) *
        ((((hA.eigenvalues α - E : ℝ) : ℂ) - ((1 / (2 * T) : ℝ) : ℂ) * Complex.I)⁻¹ *
          (((hA.eigenvalues β - E : ℝ) : ℂ) + ((1 / (2 * T) : ℝ) : ℂ) * Complex.I)⁻¹) := by
  have hz : ((E : ℂ) + Complex.I / (2 * T)).im ≠ 0 := by
    have : ((E : ℂ) + Complex.I / (2 * T)).im = 1 / (2 * T) := by
      rw [show (Complex.I / (2 * T) : ℂ) = ((1 / (2 * T) : ℝ) : ℂ) * Complex.I by push_cast; ring]
      simp
    rw [this]
    positivity
  simp only [Complex.ofReal_sum, Complex.ofReal_mul]
  refine Finset.sum_congr rfl fun y _ => ?_
  rw [resolvent_apply_spectral A hA _ hz x y, ofReal_norm_sq_sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun α _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun β _ => ?_
  rw [star_mul']
  have e1 : ((hA.eigenvalues α : ℂ) - ((E : ℂ) + Complex.I / (2 * T))) =
      ((hA.eigenvalues α - E : ℝ) : ℂ) - ((1 / (2 * T) : ℝ) : ℂ) * Complex.I := by
    push_cast; ring
  have e2 : star (((hA.eigenvalues β : ℂ) - ((E : ℂ) + Complex.I / (2 * T)))⁻¹) =
      (((hA.eigenvalues β - E : ℝ) : ℂ) + ((1 / (2 * T) : ℝ) : ℂ) * Complex.I)⁻¹ := by
    rw [star_inv₀]
    congr 1
    rw [Complex.star_def]
    simp only [map_sub, map_add, map_div₀, Complex.conj_ofReal, Complex.conj_I, map_mul,
      map_ofNat]
    push_cast
    ring
  rw [e1, e2]
  ring

end Spectral

/-- The integral of a finite triple sum of multiples of integrable functions. -/
lemma integral_triple_sum {X : Type*} [MeasurableSpace X] {μ : Measure X} {ι κ : Type*}
    (s : Finset ι) (s' : Finset κ) (K : ι → κ → κ → ℂ) (e : κ → κ → X → ℂ)
    (he : ∀ α β, Integrable (e α β) μ) :
    ∫ t, ∑ y ∈ s, ∑ α ∈ s', ∑ β ∈ s', K y α β * e α β t ∂μ =
      ∑ y ∈ s, ∑ α ∈ s', ∑ β ∈ s', K y α β * ∫ t, e α β t ∂μ := by
  rw [integral_finsetSum _ fun y _ => integrable_finsetSum _ fun α _ =>
    integrable_finsetSum _ fun β _ => (he α β).const_mul _]
  refine Finset.sum_congr rfl fun y _ => ?_
  rw [integral_finsetSum _ fun α _ => integrable_finsetSum _ fun β _ => (he α β).const_mul _]
  refine Finset.sum_congr rfl fun α _ => ?_
  rw [integral_finsetSum _ fun β _ => (he α β).const_mul _]
  refine Finset.sum_congr rfl fun β _ => ?_
  rw [integral_const_mul]

/-- Integrability of a finite triple sum of multiples of integrable functions. -/
lemma integrable_triple_sum {X : Type*} [MeasurableSpace X] {μ : Measure X} {ι κ : Type*}
    (s : Finset ι) (s' : Finset κ) (K : ι → κ → κ → ℂ) (e : κ → κ → X → ℂ)
    (he : ∀ α β, Integrable (e α β) μ) :
    Integrable (fun t => ∑ y ∈ s, ∑ α ∈ s', ∑ β ∈ s', K y α β * e α β t) μ :=
  integrable_finsetSum _ fun _ _ => integrable_finsetSum _ fun α _ =>
    integrable_finsetSum _ fun β _ => (he α β).const_mul _

/-- `1/(1 + iT(a - b)) = (2πT)⁻¹ · 2πi/((b - a) + 2i/(2T))`. -/
lemma kernel_eq (T a b : ℝ) (hT : 0 < T) :
    (1 : ℂ) / (1 + T * ((a - b : ℝ) : ℂ) * Complex.I) =
      ((1 / (2 * Real.pi * T) : ℝ) : ℂ) *
        (2 * Real.pi * Complex.I / (((b - a : ℝ) : ℂ) + 2 * ((1 / (2 * T) : ℝ) : ℂ) * Complex.I)) := by
  have hT' : (T : ℂ) ≠ 0 := by exact_mod_cast hT.ne'
  have hpi : (Real.pi : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
  have h1 : (1 + T * ((a - b : ℝ) : ℂ) * Complex.I) ≠ 0 := ne_zero_of_re_ne (by simp)
  have h2 : (((b - a : ℝ) : ℂ) + 2 * ((1 / (2 * T) : ℝ) : ℂ) * Complex.I) ≠ 0 :=
    ne_zero_of_im_ne (by simp; positivity)
  have h3 : ((T : ℂ) * ((b : ℂ) - a) + Complex.I) ≠ 0 := ne_zero_of_im_ne (by simp)
  rw [div_eq_iff h1]
  push_cast at h2 ⊢
  field_simp
  linear_combination (-(T : ℂ) * (a - b)) * Complex.I_sq

section Plancherel

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- **Plancherel identity** for a real symmetric matrix (general form of `dimerMom_plancherel`). -/
theorem plancherel_general (A : Matrix n n ℝ) (hA : A.IsSymm) (w : n → ℝ) (T : ℝ) (hT : 0 < T)
    (x : n) :
    ∫ t in Set.Ioi (0 : ℝ), 1 / T * Real.exp (-t / T) * ∑ y, w y * ‖evo A t y x‖ ^ 2 =
      1 / (2 * Real.pi * T) * ∫ E : ℝ, ∑ y, w y *
        ‖(cplx A - ((E : ℂ) + Complex.I / (2 * T)) • (1 : Matrix n n ℂ))⁻¹ y x‖ ^ 2 := by
  have hH := cplx_isHermitian A hA
  set μ := hH.eigenvalues with hμ
  have hε : 0 < 1 / (2 * T) := by positivity
  apply Complex.ofReal_injective
  rw [← integral_complex_ofReal, Complex.ofReal_mul, ← integral_complex_ofReal]
  simp_rw [time_integrand_eq A hH w T _ x, green_integrand_eq A hH w T hT _ x]
  rw [integral_triple_sum _ _ _ _ fun α β =>
      (integrableOn_exp_mul_complex_Ioi (by simp; positivity) 0).const_mul _,
    integral_triple_sum _ _ _ _ fun α β =>
      integrable_resolvent_prod (μ α) (μ β) (1 / (2 * T)) hε]
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun y _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun α _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun β _ => ?_
  rw [integral_exp_time T _ hT, integral_resolvent_prod _ _ _ hε, kernel_eq T (μ α) (μ β) hT]
  ring

end Plancherel

/-- `Im(E + i/(2T)) = 1/(2T)`. -/
lemma im_add_I_div (E T : ℝ) : ((E : ℂ) + Complex.I / (2 * T)).im = 1 / (2 * T) := by
  rw [show (Complex.I / (2 * T) : ℂ) = ((1 / (2 * T) : ℝ) : ℂ) * Complex.I by push_cast; ring]
  simp

/-- `|z - λ| ≤ 1/T` for `z = E + i/(2T)` with `|E - λ| ≤ 1/(2T)`. -/
lemma norm_z_sub_le (E lam T : ℝ) (hT : 0 < T) (hE : |E - lam| ≤ 1 / (2 * T)) :
    ‖((E : ℂ) + Complex.I / (2 * T)) - lam‖ ≤ 1 / T := by
  have e : ((E : ℂ) + Complex.I / (2 * T)) - lam =
      ((E - lam : ℝ) : ℂ) + ((1 / (2 * T) : ℝ) : ℂ) * Complex.I := by
    push_cast
    ring
  have h2 : 1 / (2 * T) + 1 / (2 * T) = 1 / T := by
    field_simp
    ring
  rw [e]
  calc ‖((E - lam : ℝ) : ℂ) + ((1 / (2 * T) : ℝ) : ℂ) * Complex.I‖
      ≤ ‖((E - lam : ℝ) : ℂ)‖ + ‖((1 / (2 * T) : ℝ) : ℂ) * Complex.I‖ := norm_add_le _ _
    _ = |E - lam| + 1 / (2 * T) := by
        rw [Complex.norm_real, norm_mul, Complex.norm_real, Complex.norm_I, mul_one,
          Real.norm_eq_abs, Real.norm_eq_abs, abs_of_pos (by positivity : (0 : ℝ) < 1 / (2 * T))]
    _ ≤ 1 / T := by linarith

section Resolvent

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- `(A - z) (A - z)⁻¹ = 1` for a real symmetric `A` and `Im z ≠ 0`. -/
lemma resolvent_mul_inv (A : Matrix n n ℝ) (hA : A.IsSymm) (z : ℂ) (hz : z.im ≠ 0) :
    (cplx A - z • (1 : Matrix n n ℂ)) * (cplx A - z • (1 : Matrix n n ℂ))⁻¹ = 1 := by
  have hH := cplx_isHermitian A hA
  rw [resolvent_eq_spectral A hH z hz]
  exact resolvent_mul_spectral A hH z hz

/-- The energy integrand of the Plancherel identity is integrable. -/
lemma integrable_green_sum (A : Matrix n n ℝ) (hA : A.IsSymm) (w : n → ℝ) (T : ℝ) (hT : 0 < T)
    (x : n) :
    Integrable (fun E : ℝ => ∑ y, w y *
      ‖(cplx A - ((E : ℂ) + Complex.I / (2 * T)) • (1 : Matrix n n ℂ))⁻¹ y x‖ ^ 2) := by
  have hH := cplx_isHermitian A hA
  have hε : 0 < 1 / (2 * T) := by positivity
  have hc : Integrable (fun E : ℝ => ((∑ y, w y *
      ‖(cplx A - ((E : ℂ) + Complex.I / (2 * T)) • (1 : Matrix n n ℂ))⁻¹ y x‖ ^ 2 : ℝ) : ℂ)) := by
    simp_rw [green_integrand_eq A hH w T hT _ x]
    exact integrable_triple_sum _ _ _ _ fun α β => integrable_resolvent_prod _ _ _ hε
  refine hc.re.congr (ae_of_all _ fun E => ?_)
  simp only [RCLike.re_to_complex, Complex.ofReal_re]

end Resolvent

section LowerBound

open scoped Matrix.Norms.L2Operator

/-- Extension by zero of a function on `[a, b]` to `ℤ`. -/
def ext (a b : ℤ) (f : Set.Icc a b → ℂ) (k : ℤ) : ℂ :=
  if h : a ≤ k ∧ k ≤ b then f ⟨k, h⟩ else 0

/-- `ext` agrees with `f` on `[a, b]`. -/
lemma ext_val (a b : ℤ) (f : Set.Icc a b → ℂ) (y : Set.Icc a b) : ext a b f y.1 = f y := by
  simp [ext, y.2.1, y.2.2]

/-- `ext` vanishes outside `[a, b]`. -/
lemma ext_out (a b : ℤ) (f : Set.Icc a b → ℂ) (k : ℤ) (hk : ¬ (a ≤ k ∧ k ≤ b)) :
    ext a b f k = 0 := by
  simp [ext, hk]

/-- Picking out one site of `[a, b]` in a sum. -/
lemma sum_ite_val (a b : ℤ) (f : Set.Icc a b → ℂ) (c : ℂ) (k : ℤ) :
    ∑ u : Set.Icc a b, (if u.1 = k then c * f u else 0) = c * ext a b f k := by
  by_cases h : a ≤ k ∧ k ≤ b
  · rw [Finset.sum_eq_single ⟨k, h⟩]
    · simp [ext, h]
    · intro u _ hu
      rw [ite_eq_right]
      intro huk
      exact hu (Subtype.ext huk)
    · simp
  · rw [ext_out a b f k h, mul_zero]
    refine Finset.sum_eq_zero fun u _ => ?_
    rw [ite_eq_right]
    intro huk
    exact h ⟨huk ▸ u.2.1, huk ▸ u.2.2⟩

/-- Row `y` of `H - z` applied to `f`. -/
lemma schrH_sub_apply_mul (a b : ℤ) (v : ℤ → ℝ) (z : ℂ) (f : Set.Icc a b → ℂ)
    (y u : Set.Icc a b) :
    (cplx (schrH a b v) - z • (1 : Matrix (Set.Icc a b) (Set.Icc a b) ℂ)) y u * f u =
      (if u.1 = y.1 then ((v y.1 : ℂ) - z) * f u else 0) +
        (if u.1 = y.1 + 1 then (-1) * f u else 0) +
          (if u.1 = y.1 - 1 then (-1) * f u else 0) := by
  simp only [Matrix.sub_apply, Matrix.smul_apply, cplx, Matrix.map_apply, schrH,
    Matrix.one_apply, smul_eq_mul]
  by_cases h1 : y.1 = u.1
  · have hyu : y = u := Subtype.ext h1
    subst hyu
    have h3 : ¬ y.1 = y.1 + 1 := by omega
    have h4 : ¬ y.1 = y.1 - 1 := by omega
    simp [h3, h4]
  · have hyu : y ≠ u := fun e => h1 (congrArg Subtype.val e)
    rw [ite_eq_right h1, ite_eq_right hyu, ite_eq_right (Ne.symm h1)]
    by_cases h2 : |y.1 - u.1| = 1
    · rw [ite_eq_left h2]
      rcases (DimerTLAux.abs_sub_eq_one_iff y.1 u.1).mp h2 with h | h
      · rw [ite_eq_left h, ite_eq_right (by omega)]
        push_cast
        ring
      · rw [ite_eq_right (by omega), ite_eq_left h]
        push_cast
        ring
    · rw [ite_eq_right h2]
      have h3 : ¬ u.1 = y.1 + 1 := fun h => h2 (by rw [h]; simp)
      have h4 : ¬ u.1 = y.1 - 1 := fun h => h2 (by rw [h]; simp)
      rw [ite_eq_right h3, ite_eq_right h4]
      simp

/-- **The equation `(H - z) G = δ_x`** for the Green function extended by zero. -/
lemma green_eq (a b : ℤ) (v : ℤ → ℝ) (z : ℂ)
    (hR : (cplx (schrH a b v) - z • (1 : Matrix (Set.Icc a b) (Set.Icc a b) ℂ)) *
      (cplx (schrH a b v) - z • (1 : Matrix (Set.Icc a b) (Set.Icc a b) ℂ))⁻¹ = 1)
    (x y : Set.Icc a b) :
    ((v y.1 : ℂ) - z) * ext a b (dimerGreen a b v z x) y.1 -
      ext a b (dimerGreen a b v z x) (y.1 + 1) - ext a b (dimerGreen a b v z x) (y.1 - 1) =
      if y = x then 1 else 0 := by
  have h := congrFun (congrFun hR y) x
  rw [Matrix.mul_apply, Matrix.one_apply] at h
  have h2 : ∑ u, (cplx (schrH a b v) - z • (1 : Matrix (Set.Icc a b) (Set.Icc a b) ℂ)) y u *
      dimerGreen a b v z x u = if y = x then 1 else 0 := h
  simp_rw [schrH_sub_apply_mul a b v z (dimerGreen a b v z x) y] at h2
  rw [Finset.sum_add_distrib, Finset.sum_add_distrib, sum_ite_val, sum_ite_val, sum_ite_val,
    ext_val] at h2
  rw [← h2, ext_val]
  ring

/-- **Transfer matrices**: if `g` satisfies the three-term recursion on the sites `m, …, m + L - 1`,
then `(g(m + L), g(m + L - 1)) = T^z(m + L, m) (g(m), g(m - 1))`. -/
lemma transfer (v : ℤ → ℝ) (z : ℂ) (g : ℤ → ℂ) (m : ℤ) (L : ℕ)
    (hrec : ∀ j : ℕ, j < L → g (m + j + 1) = ((v (m + j) : ℂ) - z) * g (m + j) - g (m + j - 1)) :
    ![g (m + L), g (m + L - 1)] = tmProd v z m L *ᵥ ![g m, g (m - 1)] := by
  induction L with
  | zero => simp [tmProd]
  | succ L ih =>
    rw [tmProd, ← Matrix.mulVec_mulVec, ← ih (fun j hj => hrec j (by omega))]
    have h := hrec L (by omega)
    have e1 : m + ((L + 1 : ℕ) : ℤ) = m + L + 1 := by push_cast; ring
    have e2 : m + ((L + 1 : ℕ) : ℤ) - 1 = m + L := by push_cast; ring
    rw [e2, e1, h]
    ext i
    fin_cases i
    · simp [tmC, Matrix.mulVec, dotProduct, Fin.sum_univ_two]
      ring
    · simp [tmC, Matrix.mulVec, dotProduct, Fin.sum_univ_two]

/-- The squared Euclidean norm of a vector of `ℂ²`. -/
def nsq (u : Fin 2 → ℂ) : ℝ := ‖u 0‖ ^ 2 + ‖u 1‖ ^ 2

/-- `‖u‖² ≥ 0`. -/
lemma nsq_nonneg (u : Fin 2 → ℂ) : 0 ≤ nsq u := by
  unfold nsq
  positivity

/-- `‖A u‖² ≤ ‖A‖² ‖u‖²`. -/
lemma nsq_mulVec_le (A : Matrix (Fin 2) (Fin 2) ℂ) (u : Fin 2 → ℂ) :
    nsq (A *ᵥ u) ≤ ‖A‖ ^ 2 * nsq u :=
  DimerTransferAux.sq_le_norm_sq_mul A u

/-- For `det A = 1`: `‖u‖² ≤ ‖A‖² ‖A u‖²` (i.e. `‖A⁻¹‖ ≤ ‖A‖` in `SL(2, ℂ)`). -/
lemma nsq_le_of_det (A : Matrix (Fin 2) (Fin 2) ℂ) (hA : A.det = 1) (u : Fin 2 → ℂ) :
    nsq u ≤ ‖A‖ ^ 2 * nsq (A *ᵥ u) := by
  set w := A *ᵥ u with hw
  have hu : u = A.adjugate *ᵥ w := by
    rw [hw, Matrix.mulVec_mulVec, Matrix.adjugate_mul, hA, one_smul, Matrix.one_mulVec]
  set w' : Fin 2 → ℂ := ![star (w 1), -star (w 0)] with hw'
  have e0 : (Aᴴ *ᵥ w') 0 = star ((A.adjugate *ᵥ w) 1) := by
    simp [hw', Matrix.mulVec, dotProduct, Fin.sum_univ_two, Matrix.conjTranspose_apply,
      Matrix.adjugate_fin_two]
    ring
  have e1 : (Aᴴ *ᵥ w') 1 = -star ((A.adjugate *ᵥ w) 0) := by
    simp [hw', Matrix.mulVec, dotProduct, Fin.sum_univ_two, Matrix.conjTranspose_apply,
      Matrix.adjugate_fin_two]
  have h1 : nsq (A.adjugate *ᵥ w) = nsq (Aᴴ *ᵥ w') := by
    simp only [nsq, e0, e1, norm_neg, norm_star]
    ring
  have h2 : nsq w' = nsq w := by
    simp only [nsq, hw', Matrix.cons_val_zero, Matrix.cons_val_one, norm_neg, norm_star]
    ring
  rw [hu, h1, ← h2, ← Matrix.l2_opNorm_conjTranspose A]
  exact nsq_mulVec_le _ _

/-- `det T_c = 1`. -/
lemma det_tmC (c : ℂ) : (tmC c).det = 1 := by
  simp [tmC, Matrix.det_fin_two_of]

/-- `det T^z(m + L, m) = 1`. -/
lemma det_tmProd (v : ℤ → ℝ) (z : ℂ) (m : ℤ) (L : ℕ) : (tmProd v z m L).det = 1 := by
  induction L with
  | zero => simp [tmProd]
  | succ L ih => rw [tmProd, Matrix.det_mul, det_tmC, ih, one_mul]

/-- Dimer configurations restrict to subintervals. -/
lemma isDimerCfg_mono {lam : ℝ} {r : ℤ} {v : ℤ → ℝ} {m m' : ℤ} {L L' : ℕ}
    (h : IsDimerCfg lam r v m L) (hm : m ≤ m') (hL : m' + L' ≤ m + L) :
    IsDimerCfg lam r v m' L' := by
  obtain ⟨h1, h2⟩ := h
  exact ⟨fun y hy1 hy2 => h1 y (by omega) (by omega), fun k hk1 hk2 => h2 k (by omega) (by omega)⟩

/-- The transfer matrices over windows of length `≤ T` near `E_c = λ` are bounded by `K₁ e^{K₂}`. -/
lemma tmProd_norm_le (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (r : ℤ) (v : ℤ → ℝ) (m : ℤ)
    (L : ℕ) (hv : IsDimerCfg lam r v m L) (z : ℂ) (T : ℝ) (hT : 1 ≤ T)
    (hz : ‖z - lam‖ ≤ 1 / T) (hLT : (L : ℝ) ≤ T) :
    ‖tmProd v z m L‖ ≤
      dimerKappa lam * (2 + 2 * lam) ^ 2 * Real.exp (dimerKappa lam * (3 + 4 * lam) / 2) := by
  have hT0 : 0 < T := by linarith
  have hz1 : ‖z - lam‖ ≤ 1 := hz.trans (by rw [div_le_one hT0]; exact hT)
  refine (dimer_transfer_bound lam h0 h1 r v m L hv z hz1).trans ?_
  have hκ : 0 ≤ dimerKappa lam := Real.sqrt_nonneg _
  have hK2 : 0 ≤ dimerKappa lam * (3 + 4 * lam) / 2 := by positivity
  have hzL : ‖z - lam‖ * L ≤ 1 := by
    calc ‖z - lam‖ * L ≤ (1 / T) * T := by gcongr
      _ = 1 := by field_simp
  gcongr
  calc dimerKappa lam * (3 + 4 * lam) / 2 * ‖z - ↑lam‖ * ↑L
      = dimerKappa lam * (3 + 4 * lam) / 2 * (‖z - ↑lam‖ * ↑L) := by ring
    _ ≤ dimerKappa lam * (3 + 4 * lam) / 2 * 1 := by gcongr
    _ = dimerKappa lam * (3 + 4 * lam) / 2 := mul_one _

/-- `∫_0^L s^p ds ≤ ∑_{j=1}^L j^p`, i.e. `L^{p+1}/(p+1) ≤ ∑_{i<L} (i+1)^p`. -/
lemma sum_rpow_ge (p : ℝ) (hp : 0 < p) (L : ℕ) :
    (L : ℝ) ^ (p + 1) / (p + 1) ≤ ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p := by
  have hmono : MonotoneOn (fun s : ℝ => s ^ p) (Set.Icc 0 (0 + (L : ℝ))) := by
    intro s hs t ht hst
    exact Real.rpow_le_rpow hs.1 hst hp.le
  have h := hmono.integral_le_sum
  rw [integral_rpow (Or.inl (by linarith)), zero_add,
    Real.zero_rpow (by linarith), sub_zero] at h
  refine h.trans (le_of_eq (Finset.sum_congr rfl fun i _ => ?_))
  push_cast
  ring_nf

/-- Partial sums of the moment along injective families of sites are bounded by the moment. -/
lemma sum_le_moment (a b : ℤ) (x : Set.Icc a b) (p : ℝ) (f : Set.Icc a b → ℂ) (s : Finset ℕ)
    (g : ℕ → ℤ) (hg : Set.InjOn g s) :
    ∑ j ∈ s, |((g j - x.1 : ℤ) : ℝ)| ^ p * ‖ext a b f (g j)‖ ^ 2 ≤
      ∑ y : Set.Icc a b, |((y.1 - x.1 : ℤ) : ℝ)| ^ p * ‖f y‖ ^ 2 := by
  set h : ℤ → ℝ := fun k => |((k - x.1 : ℤ) : ℝ)| ^ p * ‖ext a b f k‖ ^ 2 with hh
  have h0 : ∀ k, 0 ≤ h k := fun k => by positivity
  have hout : ∀ k, ¬ (a ≤ k ∧ k ≤ b) → h k = 0 := fun k hk => by
    simp [hh, ext_out a b f k hk]
  calc ∑ j ∈ s, |((g j - x.1 : ℤ) : ℝ)| ^ p * ‖ext a b f (g j)‖ ^ 2
      = ∑ k ∈ s.image g, h k := (Finset.sum_image (f := h) fun i hi j hj e => hg hi hj e).symm
    _ ≤ ∑ k ∈ Finset.Icc a b, h k := by
        rw [← Finset.sum_filter_add_sum_filter_not (s.image g) (fun k => a ≤ k ∧ k ≤ b),
          Finset.sum_eq_zero (s := Finset.filter (fun k => ¬ (a ≤ k ∧ k ≤ b)) _)
            (fun k hk => hout k (Finset.mem_filter.mp hk).2), add_zero]
        apply Finset.sum_le_sum_of_subset_of_nonneg
        · intro k hk
          rw [Finset.mem_filter] at hk
          exact Finset.mem_Icc.mpr hk.2
        · intro k _ _
          exact h0 k
    _ = ∑ y : Set.Icc a b, h y.1 :=
        Finset.sum_subtype _ (fun k => by simp) h
    _ = ∑ y : Set.Icc a b, |((y.1 - x.1 : ℤ) : ℝ)| ^ p * ‖f y‖ ^ 2 := by
        refine Finset.sum_congr rfl fun y _ => ?_
        simp only [hh, ext_val]

/-- `(i+1)^p c ≤ (i+1)^p A + (i+2)^p B` from `c ≤ A + B`. -/
lemma weighted_le (c A B' p : ℝ) (i : ℕ) (hp : 0 < p) (hB : 0 ≤ B') (h : c ≤ A + B') :
    ((i : ℝ) + 1) ^ p * c ≤ ((i : ℝ) + 1) ^ p * A + ((i : ℝ) + 2) ^ p * B' := by
  have h1 : ((i : ℝ) + 1) ^ p ≤ ((i : ℝ) + 2) ^ p :=
    Real.rpow_le_rpow (by positivity) (by linarith) hp.le
  have h2 : 0 ≤ ((i : ℝ) + 1) ^ p := by positivity
  nlinarith

/-- `|n|^p` for a natural number `n` written as an integer difference. -/
lemma abs_cast_rpow (k : ℤ) (q : ℝ) (hk : (k : ℝ) = q) (hq : 0 ≤ q) (p : ℝ) :
    |((k : ℤ) : ℝ)| ^ p = q ^ p := by
  rw [hk, abs_of_nonneg hq]

/-- `|-q|^p = q^p` for `q ≥ 0`, for an integer cast. -/
lemma abs_cast_rpow' (k : ℤ) (q : ℝ) (hk : (k : ℝ) = -q) (hq : 0 ≤ q) (p : ℝ) :
    |((k : ℤ) : ℝ)| ^ p = q ^ p := by
  rw [hk, abs_neg, abs_of_nonneg hq]

/-- **Steps 2–3 of the proof**: the moment of the Green function at `z` with `|z - λ| ≤ 1/T`. -/
lemma green_moment_ge (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (p : ℝ) (hp : 0 < p) (r : ℤ)
    (a b : ℤ) (v : ℤ → ℝ) (hv : IsDimerCfg lam r v a (b - a + 1).toNat) (x : Set.Icc a b)
    (T : ℝ) (hT : 1 ≤ T) (hxT : (a : ℝ) ≤ (x.1 : ℝ) - T ∧ (x.1 : ℝ) + T ≤ (b : ℝ)) (z : ℂ)
    (hR : (cplx (schrH a b v) - z • (1 : Matrix (Set.Icc a b) (Set.Icc a b) ℂ)) *
      (cplx (schrH a b v) - z • (1 : Matrix (Set.Icc a b) (Set.Icc a b) ℂ))⁻¹ = 1)
    (hz : ‖z - lam‖ ≤ 1 / T) :
    dimerC1 lam ^ 2 * T ^ (p + 1) / (2 ^ (p + 2) * (p + 1)) ≤
      ∑ y : Set.Icc a b, |((y.1 - x.1 : ℤ) : ℝ)| ^ p * ‖dimerGreen a b v z x y‖ ^ 2 := by
  set G := ext a b (dimerGreen a b v z x) with hG
  set S := ∑ y : Set.Icc a b, |((y.1 - x.1 : ℤ) : ℝ)| ^ p * ‖dimerGreen a b v z x y‖ ^ 2
    with hS
  set B := dimerKappa lam * (2 + 2 * lam) ^ 2 * Real.exp (dimerKappa lam * (3 + 4 * lam) / 2)
    with hB
  set c₀ := 1 / (3 + 2 * lam) with hc₀
  set L := ⌊T⌋₊ with hL
  have hT0 : 0 < T := by linarith
  have hκ := DimerTransferAux.kappa_ge_one lam h0 h1
  have hB0 : 0 < B := by positivity
  have hc₀0 : 0 < c₀ := by positivity
  have hc₀eq : (3 + 2 * lam) * c₀ = 1 := by rw [hc₀]; field_simp
  have hc1 : dimerC1 lam = c₀ / B := by
    simp only [dimerC1, hc₀, hB, div_div, mul_assoc]
  -- the window `[x - L, x + L]`
  have hLT : (L : ℝ) ≤ T := Nat.floor_le hT0.le
  have hL1 : 1 ≤ L := Nat.le_floor (by simpa using hT)
  have hTL : T < L + 1 := Nat.lt_floor_add_one T
  have hLT2 : T / 2 ≤ L := by
    have : (1 : ℝ) ≤ L := by exact_mod_cast hL1
    linarith
  have hxa : a ≤ x.1 - (L : ℤ) := by
    have : (a : ℝ) ≤ ((x.1 - (L : ℤ) : ℤ) : ℝ) := by push_cast; linarith [hxT.1]
    exact_mod_cast this
  have hxb : x.1 + (L : ℤ) ≤ b := by
    have : ((x.1 + (L : ℤ) : ℤ) : ℝ) ≤ (b : ℝ) := by push_cast; linarith [hxT.2]
    exact_mod_cast this
  have hxa' : a ≤ x.1 := x.2.1
  have hN : a + ((b - a + 1).toNat : ℤ) = b + 1 := by
    rw [Int.toNat_of_nonneg (by omega)]
    ring
  -- the equation `(H - z) G = δ_x`
  have hrec : ∀ k, a ≤ k → k ≤ b → k ≠ x.1 →
      G (k + 1) = ((v k : ℂ) - z) * G k - G (k - 1) := by
    intro k hka hkb hkx
    have h := green_eq a b v z hR x ⟨k, hka, hkb⟩
    rw [ite_eq_right (fun e => hkx (congrArg Subtype.val e))] at h
    linear_combination -h
  have hX : ((v x.1 : ℂ) - z) * G x.1 - G (x.1 + 1) - G (x.1 - 1) = 1 := by
    have h := green_eq a b v z hR x x
    rwa [ite_eq_left rfl] at h
  -- transfer-matrix bounds
  have hTM : ∀ (m : ℤ) (j : ℕ), a ≤ m → m + j ≤ b + 1 → (j : ℝ) ≤ T →
      ‖tmProd v z m j‖ ≤ B := fun m j hm hmj hjT =>
    tmProd_norm_le lam h0 h1 r v m j (isDimerCfg_mono hv hm (by rw [hN]; exact hmj)) z T hT hz
      hjT
  have hprop : ∀ u w : Fin 2 → ℂ, c₀ ^ 2 ≤ nsq u → nsq u ≤ B ^ 2 * nsq w →
      dimerC1 lam ^ 2 ≤ nsq w := by
    intro u w hu huw
    rw [hc1, div_pow, div_le_iff₀ (by positivity)]
    linarith
  -- Step 3: one of `|G(x+1)|, |G(x)|, |G(x-1)|` is at least `c₀`
  have hvx : ‖(v x.1 : ℂ) - z‖ ≤ 2 * lam + 1 := by
    have hw : v x.1 = lam ∨ v x.1 = -lam := hv.1 x.1 x.2.1 (by rw [hN]; linarith [x.2.2])
    have h1' := DimerTransferAux.norm_w_sub_le lam h0 hw
    have h2' : 1 / T ≤ 1 := by rw [div_le_one hT0]; exact hT
    calc ‖(v x.1 : ℂ) - z‖ = ‖((v x.1 : ℂ) - lam) - (z - lam)‖ := by ring_nf
      _ ≤ ‖(v x.1 : ℂ) - lam‖ + ‖z - lam‖ := norm_sub_le _ _
      _ ≤ 2 * lam + 1 := by linarith
  have hthree : 1 ≤ (2 * lam + 1) * ‖G x.1‖ + ‖G (x.1 + 1)‖ + ‖G (x.1 - 1)‖ := by
    have h := congrArg norm hX
    rw [norm_one] at h
    calc (1 : ℝ) = ‖((v x.1 : ℂ) - z) * G x.1 - G (x.1 + 1) - G (x.1 - 1)‖ := h.symm
      _ ≤ ‖((v x.1 : ℂ) - z) * G x.1 - G (x.1 + 1)‖ + ‖G (x.1 - 1)‖ := norm_sub_le _ _
      _ ≤ ‖((v x.1 : ℂ) - z) * G x.1‖ + ‖G (x.1 + 1)‖ + ‖G (x.1 - 1)‖ := by
          gcongr
          exact norm_sub_le _ _
      _ = ‖(v x.1 : ℂ) - z‖ * ‖G x.1‖ + ‖G (x.1 + 1)‖ + ‖G (x.1 - 1)‖ := by rw [norm_mul]
      _ ≤ (2 * lam + 1) * ‖G x.1‖ + ‖G (x.1 + 1)‖ + ‖G (x.1 - 1)‖ := by gcongr
  -- the summed bound `c₁² ∑_{i<L} (i+1)^p ≤ 2 S`
  have hmain : dimerC1 lam ^ 2 * ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p ≤ 2 * S := by
    rw [Finset.mul_sum]
    by_cases hcase : c₀ ≤ ‖G (x.1 + 1)‖ ∨ c₀ ≤ ‖G x.1‖
    · -- `‖Φ(x+1)‖ ≥ c₀`: propagate to the right
      have hΦ : c₀ ^ 2 ≤ nsq ![G (x.1 + 1), G (x.1 + 1 - 1)] := by
        simp only [nsq, Matrix.cons_val_zero, Matrix.cons_val_one, add_sub_cancel_right]
        rcases hcase with h | h
        · nlinarith [sq_nonneg ‖G x.1‖]
        · nlinarith [sq_nonneg ‖G (x.1 + 1)‖]
      have hj : ∀ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p * dimerC1 lam ^ 2 ≤
          |((x.1 + 1 + ((i + 1 : ℕ) : ℤ) - 1 - x.1 : ℤ) : ℝ)| ^ p *
              ‖G (x.1 + 1 + ((i + 1 : ℕ) : ℤ) - 1)‖ ^ 2 +
            |((x.1 + 1 + ((i + 1 : ℕ) : ℤ) - x.1 : ℤ) : ℝ)| ^ p *
              ‖G (x.1 + 1 + ((i + 1 : ℕ) : ℤ))‖ ^ 2 := by
        intro i hi
        rw [Finset.mem_range] at hi
        have htr := transfer v z G (x.1 + 1) (i + 1) (fun i' hi' =>
          hrec (x.1 + 1 + i') (by omega) (by omega) (by omega))
        have hM := hTM (x.1 + 1) (i + 1) (by omega) (by push_cast; omega)
          (by push_cast; linarith [(by exact_mod_cast Nat.succ_le_of_lt hi : (i : ℝ) + 1 ≤ L)])
        have h1' := nsq_le_of_det _ (det_tmProd v z (x.1 + 1) (i + 1))
          ![G (x.1 + 1), G (x.1 + 1 - 1)]
        rw [← htr] at h1'
        have h2' : nsq ![G (x.1 + 1), G (x.1 + 1 - 1)] ≤
            B ^ 2 * nsq ![G (x.1 + 1 + ((i + 1 : ℕ) : ℤ)), G (x.1 + 1 + ((i + 1 : ℕ) : ℤ) - 1)] :=
          h1'.trans (mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (norm_nonneg _) hM 2)
            (nsq_nonneg _))
        have h3' := hprop _ _ hΦ h2'
        simp only [nsq, Matrix.cons_val_zero, Matrix.cons_val_one] at h3'
        rw [abs_cast_rpow _ ((i : ℝ) + 1) (by push_cast; ring) (by positivity),
          abs_cast_rpow _ ((i : ℝ) + 2) (by push_cast; ring) (by positivity)]
        exact weighted_le _ _ _ p i hp (by positivity) (by linarith)
      calc ∑ i ∈ Finset.range L, dimerC1 lam ^ 2 * ((i : ℝ) + 1) ^ p
          = ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p * dimerC1 lam ^ 2 :=
            Finset.sum_congr rfl fun i _ => mul_comm _ _
        _ ≤ ∑ i ∈ Finset.range L,
            (|((x.1 + 1 + ((i + 1 : ℕ) : ℤ) - 1 - x.1 : ℤ) : ℝ)| ^ p *
                ‖G (x.1 + 1 + ((i + 1 : ℕ) : ℤ) - 1)‖ ^ 2 +
              |((x.1 + 1 + ((i + 1 : ℕ) : ℤ) - x.1 : ℤ) : ℝ)| ^ p *
                ‖G (x.1 + 1 + ((i + 1 : ℕ) : ℤ))‖ ^ 2) := Finset.sum_le_sum hj
        _ = ∑ i ∈ Finset.range L,
              |(((fun i : ℕ => x.1 + 1 + ((i + 1 : ℕ) : ℤ) - 1) i - x.1 : ℤ) : ℝ)| ^ p *
                ‖G ((fun i : ℕ => x.1 + 1 + ((i + 1 : ℕ) : ℤ) - 1) i)‖ ^ 2 +
            ∑ i ∈ Finset.range L,
              |(((fun i : ℕ => x.1 + 1 + ((i + 1 : ℕ) : ℤ)) i - x.1 : ℤ) : ℝ)| ^ p *
                ‖G ((fun i : ℕ => x.1 + 1 + ((i + 1 : ℕ) : ℤ)) i)‖ ^ 2 :=
            Finset.sum_add_distrib
        _ ≤ S + S := add_le_add
            (sum_le_moment a b x p _ _ _ fun i _ j _ e => by simp at e; omega)
            (sum_le_moment a b x p _ _ _ fun i _ j _ e => by simp at e; omega)
        _ = 2 * S := by ring
    · -- otherwise `|G(x-1)| ≥ c₀`, so `‖Φ(x)‖ ≥ c₀`: propagate to the left
      rw [not_or, not_le, not_le] at hcase
      have hxm : c₀ ≤ ‖G (x.1 - 1)‖ := by nlinarith [norm_nonneg (G x.1)]
      have hΦ : c₀ ^ 2 ≤ nsq ![G x.1, G (x.1 - 1)] := by
        simp only [nsq, Matrix.cons_val_zero, Matrix.cons_val_one]
        nlinarith [sq_nonneg ‖G x.1‖]
      have hj : ∀ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p * dimerC1 lam ^ 2 ≤
          |((x.1 - ((i + 1 : ℕ) : ℤ) - x.1 : ℤ) : ℝ)| ^ p *
              ‖G (x.1 - ((i + 1 : ℕ) : ℤ))‖ ^ 2 +
            |((x.1 - ((i + 1 : ℕ) : ℤ) - 1 - x.1 : ℤ) : ℝ)| ^ p *
              ‖G (x.1 - ((i + 1 : ℕ) : ℤ) - 1)‖ ^ 2 := by
        intro i hi
        rw [Finset.mem_range] at hi
        have htr := transfer v z G (x.1 - ((i + 1 : ℕ) : ℤ)) (i + 1) (fun i' hi' =>
          hrec (x.1 - ((i + 1 : ℕ) : ℤ) + i') (by omega) (by omega) (by omega))
        have hM := hTM (x.1 - ((i + 1 : ℕ) : ℤ)) (i + 1) (by omega) (by omega)
          (by push_cast; linarith [(by exact_mod_cast Nat.succ_le_of_lt hi : (i : ℝ) + 1 ≤ L)])
        have e1 : x.1 - ((i + 1 : ℕ) : ℤ) + ((i + 1 : ℕ) : ℤ) = x.1 := by ring
        rw [e1] at htr
        have h1' := nsq_mulVec_le (tmProd v z (x.1 - ((i + 1 : ℕ) : ℤ)) (i + 1))
          ![G (x.1 - ((i + 1 : ℕ) : ℤ)), G (x.1 - ((i + 1 : ℕ) : ℤ) - 1)]
        rw [← htr] at h1'
        have h2' : nsq ![G x.1, G (x.1 - 1)] ≤
            B ^ 2 * nsq ![G (x.1 - ((i + 1 : ℕ) : ℤ)), G (x.1 - ((i + 1 : ℕ) : ℤ) - 1)] :=
          h1'.trans (mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (norm_nonneg _) hM 2)
            (nsq_nonneg _))
        have h3' := hprop _ _ hΦ h2'
        simp only [nsq, Matrix.cons_val_zero, Matrix.cons_val_one] at h3'
        rw [abs_cast_rpow' _ ((i : ℝ) + 1) (by push_cast; ring) (by positivity),
          abs_cast_rpow' _ ((i : ℝ) + 2) (by push_cast; ring) (by positivity)]
        exact weighted_le _ _ _ p i hp (by positivity) h3'
      calc ∑ i ∈ Finset.range L, dimerC1 lam ^ 2 * ((i : ℝ) + 1) ^ p
          = ∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p * dimerC1 lam ^ 2 :=
            Finset.sum_congr rfl fun i _ => mul_comm _ _
        _ ≤ ∑ i ∈ Finset.range L,
            (|((x.1 - ((i + 1 : ℕ) : ℤ) - x.1 : ℤ) : ℝ)| ^ p *
                ‖G (x.1 - ((i + 1 : ℕ) : ℤ))‖ ^ 2 +
              |((x.1 - ((i + 1 : ℕ) : ℤ) - 1 - x.1 : ℤ) : ℝ)| ^ p *
                ‖G (x.1 - ((i + 1 : ℕ) : ℤ) - 1)‖ ^ 2) := Finset.sum_le_sum hj
        _ = ∑ i ∈ Finset.range L,
              |(((fun i : ℕ => x.1 - ((i + 1 : ℕ) : ℤ)) i - x.1 : ℤ) : ℝ)| ^ p *
                ‖G ((fun i : ℕ => x.1 - ((i + 1 : ℕ) : ℤ)) i)‖ ^ 2 +
            ∑ i ∈ Finset.range L,
              |(((fun i : ℕ => x.1 - ((i + 1 : ℕ) : ℤ) - 1) i - x.1 : ℤ) : ℝ)| ^ p *
                ‖G ((fun i : ℕ => x.1 - ((i + 1 : ℕ) : ℤ) - 1) i)‖ ^ 2 :=
            Finset.sum_add_distrib
        _ ≤ S + S := add_le_add
            (sum_le_moment a b x p _ _ _ fun i _ j _ e => by simp at e; omega)
            (sum_le_moment a b x p _ _ _ fun i _ j _ e => by simp at e; omega)
        _ = 2 * S := by ring
  -- conclusion
  have hsum := sum_rpow_ge p hp L
  have hTL' : (T / 2) ^ (p + 1) ≤ (L : ℝ) ^ (p + 1) :=
    Real.rpow_le_rpow (by positivity) hLT2 (by linarith)
  have e2 : (2 : ℝ) ^ (p + 2) = 2 ^ (p + 1) * 2 := by
    rw [show p + 2 = (p + 1) + 1 by ring, Real.rpow_add two_pos, Real.rpow_one]
  have e3 : T ^ (p + 1) / (2 ^ (p + 2) * (p + 1)) = (T / 2) ^ (p + 1) / (p + 1) / 2 := by
    rw [Real.div_rpow hT0.le (by norm_num), e2]
    field_simp
  have hc1sq : 0 ≤ dimerC1 lam ^ 2 := sq_nonneg _
  calc dimerC1 lam ^ 2 * T ^ (p + 1) / (2 ^ (p + 2) * (p + 1))
      = dimerC1 lam ^ 2 * ((T / 2) ^ (p + 1) / (p + 1)) / 2 := by
        rw [mul_div_assoc, e3]
        ring
    _ ≤ dimerC1 lam ^ 2 * ((L : ℝ) ^ (p + 1) / (p + 1)) / 2 := by
        gcongr
    _ ≤ dimerC1 lam ^ 2 * (∑ i ∈ Finset.range L, ((i : ℝ) + 1) ^ p) / 2 := by
        gcongr
    _ ≤ S := by linarith

end LowerBound

end DimerLowerAux

/-- **Plancherel identity** (paper: eq. `eq:plancherel`). -/
theorem dimerMom_plancherel (a b : ℤ) (v : ℤ → ℝ) (p T : ℝ) (hp : 0 < p) (hT : 0 < T)
    (x : Set.Icc a b) :
    dimerMom a b v p T x =
      (1 / (2 * Real.pi * T)) * ∫ E : ℝ, ∑ y : Set.Icc a b, |((y.1 - x.1 : ℤ) : ℝ)| ^ p *
        ‖dimerGreen a b v (E + Complex.I / (2 * T)) x y‖ ^ 2 := by
  exact DimerLowerAux.plancherel_general (schrH a b v) (DimerTLAux.schrH_isSymm a b v)
    (fun y => |((y.1 - x.1 : ℤ) : ℝ)| ^ p) T hT x

/-- **Every dimer configuration spreads at least diffusively** (paper: Proposition `prop:dimerlower`), on a
finite interval `[a, b]`. -/
theorem dimer_lower_bound (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (p : ℝ) (hp : 0 < p) (r : ℤ)
    (a b : ℤ) (v : ℤ → ℝ) (hv : IsDimerCfg lam r v a (b - a + 1).toNat) (x : Set.Icc a b) (T : ℝ)
    (hT : 1 ≤ T) (hxT : (a : ℝ) ≤ (x.1 : ℝ) - T ∧ (x.1 : ℝ) + T ≤ (b : ℝ)) :
    dimerClp lam p * T ^ (p - 1) ≤ dimerMom a b v p T x := by
  have hT0 : 0 < T := by linarith
  rw [dimerMom_plancherel a b v p T hp hT0 x]
  set F : ℝ → ℝ := fun E => ∑ y : Set.Icc a b, |((y.1 - x.1 : ℤ) : ℝ)| ^ p *
    ‖dimerGreen a b v (E + Complex.I / (2 * T)) x y‖ ^ 2 with hF
  have hF0 : ∀ E, 0 ≤ F E := fun E => Finset.sum_nonneg fun y _ => by positivity
  have hFint : Integrable F :=
    DimerLowerAux.integrable_green_sum (schrH a b v) (DimerTLAux.schrH_isSymm a b v)
      (fun y => |((y.1 - x.1 : ℤ) : ℝ)| ^ p) T hT0 x
  set Bnd := dimerC1 lam ^ 2 * T ^ (p + 1) / (2 ^ (p + 2) * (p + 1)) with hBnd
  -- Steps 2–3: the integrand is large near `E_c = λ`
  have hwin : ∀ E ∈ Set.Icc (lam - 1 / (2 * T)) (lam + 1 / (2 * T)), Bnd ≤ F E := by
    intro E hE
    have hzim : ((E : ℂ) + Complex.I / (2 * T)).im ≠ 0 := by
      rw [DimerLowerAux.im_add_I_div]
      positivity
    have hR := DimerLowerAux.resolvent_mul_inv (schrH a b v) (DimerTLAux.schrH_isSymm a b v) _
      hzim
    have hz := DimerLowerAux.norm_z_sub_le E lam T hT0
      (abs_sub_le_iff.mpr ⟨by linarith [hE.2], by linarith [hE.1]⟩)
    exact DimerLowerAux.green_moment_ge lam h0 h1 p hp r a b v hv x T hT hxT _ hR hz
  -- Step 4: restrict the energy integral to the window `[λ - 1/(2T), λ + 1/(2T)]`
  have hint : Bnd * (1 / T) ≤ ∫ E in Set.Icc (lam - 1 / (2 * T)) (lam + 1 / (2 * T)), F E := by
    have h := setIntegral_ge_of_const_le_real measurableSet_Icc
      (by simp [Real.volume_Icc]) hwin hFint.integrableOn
    have hvol : volume.real (Set.Icc (lam - 1 / (2 * T)) (lam + 1 / (2 * T))) = 1 / T := by
      rw [Measure.real, Real.volume_Icc,
        ENNReal.toReal_ofReal (by linarith [show (0 : ℝ) < 1 / (2 * T) by positivity])]
      field_simp
      ring
    rwa [hvol] at h
  have hle : ∫ E in Set.Icc (lam - 1 / (2 * T)) (lam + 1 / (2 * T)), F E ≤ ∫ E, F E :=
    setIntegral_le_integral hFint (ae_of_all _ hF0)
  have e1 : T ^ (p + 1) = T ^ (p - 1) * T ^ 2 := by
    rw [show p + 1 = (p - 1) + 2 by ring, Real.rpow_add hT0, Real.rpow_two]
  have e2 : (2 : ℝ) ^ (p + 3) = 2 ^ (p + 2) * 2 := by
    rw [show p + 3 = (p + 2) + 1 by ring, Real.rpow_add two_pos, Real.rpow_one]
  have hpi : 0 < Real.pi := Real.pi_pos
  calc dimerClp lam p * T ^ (p - 1) = 1 / (2 * Real.pi * T) * (Bnd * (1 / T)) := by
        rw [dimerClp, hBnd, e1, e2]
        field_simp
    _ ≤ 1 / (2 * Real.pi * T) * ∫ E, F E := by
        gcongr
        exact hint.trans hle

end DFL
