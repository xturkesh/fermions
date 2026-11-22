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


variable {n : Type*} [Fintype n] [DecidableEq n]

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

end DimerLowerAux

/-- **Plancherel identity** (paper: eq. `eq:plancherel`). -/
theorem dimerMom_plancherel (a b : ℤ) (v : ℤ → ℝ) (p T : ℝ) (hp : 0 < p) (hT : 0 < T)
    (x : Set.Icc a b) :
    dimerMom a b v p T x =
      (1 / (2 * Real.pi * T)) * ∫ E : ℝ, ∑ y : Set.Icc a b, |((y.1 - x.1 : ℤ) : ℝ)| ^ p *
        ‖dimerGreen a b v (E + Complex.I / (2 * T)) x y‖ ^ 2 :=
  DimerLowerAux.plancherel_general (schrH a b v) (DimerTLAux.schrH_isSymm a b v)
    (fun y => |((y.1 - x.1 : ℤ) : ℝ)| ^ p) T hT x

/-- **Every dimer configuration spreads at least diffusively** (paper: Proposition `prop:dimerlower`), on a
finite interval `[a, b]`. -/
theorem dimer_lower_bound (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (p : ℝ) (hp : 0 < p) (r : ℤ)
    (a b : ℤ) (v : ℤ → ℝ) (hv : IsDimerCfg lam r v a (b - a + 1).toNat) (x : Set.Icc a b) (T : ℝ)
    (hT : 1 ≤ T) (hxT : (a : ℝ) ≤ (x.1 : ℝ) - T ∧ (x.1 : ℝ) + T ≤ (b : ℝ)) :
    dimerClp lam p * T ^ (p - 1) ≤ dimerMom a b v p T x := by
  sorry

end DFL
