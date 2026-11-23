import DFL.ECInfVolQ
import DFL.DimerImbalanceLim
import DFL.DimerLGT

/-!
# The density of states of the dimer model (paper: Remark `rem:dimerdos`)

On `ℓ²(ℤ)` let `H_η = Δ + λ diag(q(η))` with the dimer charges `q_{2k-1} = q_{2k} = η_k` (`dimerQZ`,
`η` i.i.d. uniform under `muZ`); conjugation by `Γ = diag((-1)^x)` turns it into `-Δ + λ diag(q)` without
changing the spectral measures at the vectors `δ_x`. The density-of-states functional is
`ν(f) = ½ E(⟨δ_1, f(H) δ_1⟩ + ⟨δ_2, f(H) δ_2⟩)` (`dosFun`, with the continuous functional calculus `cfc`
of `z ↦ f(Re z)`; the spectrum of `H` is real).

* `tendsto_dosFin` (`rem:dimerdos`(ii)): for continuous `f`, `N⁻¹ E_dim Tr f(h_q/J) → ν(f)` along even `N`
  (finite propagation for polynomials, stationarity of the dimer law under the shift by one dimer, and
  Weierstrass approximation on the spectral interval `[-2-|λ|, 2+|λ|]`).
* `dosFun_nonneg`, `dosFun_one`, `dosFun_eq_zero_of_support`, `dosFun_le_of_near`: `ν` is positive and
  normalized, it lives on `[-2-|λ|, 2+|λ|]`, and it gives mass at most `2/L` to every interval
  `[E - ε, E + ε]` with `ε L (|λ| + |E| + 1)^L < 1` (from `card_eig_near_le`, the deterministic finite-volume
  count of `lem:fewapprox`).
* `exists_dosMeasure`, `dosMeasure`, `isProbabilityMeasure_dosMeasure`, `dosMeasure_compl_Icc`,
  `dosMeasure_singleton` (`rem:dimerdos`(i)): by the Riesz–Markov–Kakutani theorem `ν` is a probability measure
  on `[-2-|λ|, 2+|λ|]` without atoms.
-/

open MeasureTheory Filter Topology Matrix
open scoped InnerProductSpace

noncomputable section

namespace DFL

/-- The dimer charges on `ℤ` with the dimers `{2k-1, 2k}`: `q_y = η_{⌊(y+1)/2⌋}`. -/
def dimerQZ (η : CfgZ) : CfgZ := fun y => η ((y + 1) / 2)

/-- The dimer Hamiltonian `H_η = Δ + λ diag(q(η))` on `ℓ²(ℤ)`. -/
def dimerHZ (lam : ℝ) (η : CfgZ) : L2Z →L[ℂ] L2Z := ECInfVolQAux.AZop lam (dimerQZ η)

/-- The density-of-states functional `ν(f) = ½ E(⟨δ_1, f(H) δ_1⟩ + ⟨δ_2, f(H) δ_2⟩)`, with `f(H)` the
continuous functional calculus of the self-adjoint `H` applied to `z ↦ f(Re z)`. -/
def dosFun (lam : ℝ) (f : ℝ → ℝ) : ℝ :=
  (1 / 2) * ∫ η, (⟪deltaZ 1, cfc (fun z : ℂ => ((f z.re : ℝ) : ℂ)) (dimerHZ lam η) (deltaZ 1)⟫_ℂ +
    ⟪deltaZ 2, cfc (fun z : ℂ => ((f z.re : ℝ) : ℂ)) (dimerHZ lam η) (deltaZ 2)⟫_ℂ).re ∂muZ

/-- The finite-volume eigenvalue average `N⁻¹ E_dim Tr f(h_q/J)`, `h_q/J = -Δ_N + λ diag(q)`, `N = 2M`. -/
def dosFin (lam : ℝ) (M : ℕ) (f : ℝ → ℝ) : ℝ :=
  (1 / (2 * M : ℝ)) * dimerAvg (fun q : Cfg (2 * M) =>
    ∑ α, f ((hq_isHermitian 1 (lam / 2) q).eigenvalues α))

namespace DimerDOSAux

/-! ### The infinite-volume operator -/

/-- A weighted relabelling has operator norm at most one. -/
lemma norm_wOp_le (e : ℤ ≃ ℤ) (c : ℤ → ℂ) (hc : ∀ z, ‖c z‖ ≤ 1) :
    ‖ECInfVolQAux.wOp e c hc‖ ≤ 1 :=
  LinearMap.mkContinuous_norm_le _ zero_le_one _

/-- `‖Δ + λ diag(sgn q)‖ ≤ 2 + |λ|`. -/
lemma norm_AZop_le (lam : ℝ) (q : CfgZ) : ‖ECInfVolQAux.AZop lam q‖ ≤ 2 + |lam| := by
  have h1 := norm_wOp_le (Equiv.addRight 1) (fun _ => (1 : ℂ)) (fun _ => by simp)
  have h2 := norm_wOp_le (Equiv.subRight 1) (fun _ => (1 : ℂ)) (fun _ => by simp)
  have h3 := norm_wOp_le (Equiv.refl ℤ) (fun z => ((sgn (q z) : ℝ) : ℂ))
    (fun z => ECInfVolQAux.norm_sgn_le (q z))
  have h4 : ‖(lam : ℂ) • ECInfVolQAux.wOp (Equiv.refl ℤ) (fun z => ((sgn (q z) : ℝ) : ℂ))
      (fun z => ECInfVolQAux.norm_sgn_le (q z))‖ ≤ |lam| := by
    rw [norm_smul, Complex.norm_real, Real.norm_eq_abs]
    exact mul_le_of_le_one_right (abs_nonneg _) h3
  unfold ECInfVolQAux.AZop
  refine le_trans (norm_add_le _ _) ?_
  refine le_trans (add_le_add (norm_add_le _ _) h4) ?_
  linarith

/-- `H_η` has the matrix elements of `Δ + λ diag(sgn q(η))`. -/
lemma isAZ_dimerHZ (lam : ℝ) (η : CfgZ) : IsAZ lam (dimerQZ η) (dimerHZ lam η) :=
  ECInfVolQAux.isAZ_AZop lam _

/-- `H_η` is self-adjoint. -/
lemma isSelfAdjoint_dimerHZ (lam : ℝ) (η : CfgZ) : IsSelfAdjoint (dimerHZ lam η) :=
  ECInfVolQAux.isSelfAdjoint_of_isAZ (isAZ_dimerHZ lam η)

/-- The spectrum of `H_η` lies in `[-2-|λ|, 2+|λ|] ⊆ ℝ`. -/
lemma mem_spectrum_dimerHZ (lam : ℝ) (η : CfgZ) {z : ℂ}
    (hz : z ∈ spectrum ℂ (dimerHZ lam η)) : z = (z.re : ℂ) ∧ |z.re| ≤ 2 + |lam| := by
  obtain ⟨h1, h2⟩ := ECInfVolQAux.mem_spectrum_re (isSelfAdjoint_dimerHZ lam η) hz
  refine ⟨h1, h2.trans ?_⟩
  have h3 : ‖dimerHZ lam η‖ ≤ 2 + |lam| := norm_AZop_le lam _
  have h4 : ‖(1 : L2Z →L[ℂ] L2Z)‖ ≤ 1 := ContinuousLinearMap.norm_id_le
  calc ‖dimerHZ lam η‖ * ‖(1 : L2Z →L[ℂ] L2Z)‖ ≤ (2 + |lam|) * 1 :=
        mul_le_mul h3 h4 (norm_nonneg _) (by positivity)
    _ = 2 + |lam| := mul_one _

/-- The diagonal matrix element `p(H_η)_{xx} = ∑_k p_k (H_η^k)_{xx}` of a real polynomial. -/
def momZ (lam : ℝ) (p : Polynomial ℝ) (η : CfgZ) (x : ℤ) : ℝ :=
  ∑ k ∈ Finset.range (p.natDegree + 1), p.coeff k * hZpow (-1) (lam / 2) (dimerQZ η) k x x

/-- `⟨δ_x, p(H_η) δ_x⟩ = p(H_η)_{xx}`. -/
lemma inner_aeval_deltaZ (lam : ℝ) (p : Polynomial ℝ) (η : CfgZ) (x : ℤ) :
    ⟪deltaZ x, Polynomial.aeval (dimerHZ lam η) (p.map Complex.ofRealHom) (deltaZ x)⟫_ℂ =
      (momZ lam p η x : ℂ) := by
  rw [Polynomial.aeval_eq_sum_range, Polynomial.natDegree_map, _root_.sum_apply,
    inner_sum, momZ, Complex.ofReal_sum]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [_root_.smul_apply, inner_smul_right, ECInfVolQAux.inner_pow_deltaZ (isAZ_dimerHZ lam η),
    Polynomial.coeff_map, Complex.ofReal_mul]
  rfl

/-- `|Re⟨δ_x, f(H_η) δ_x⟩ - p(H_η)_{xx}| ≤ δ` if `|p - f| ≤ δ` on `[-2-|λ|, 2+|λ|]`. -/
lemma abs_re_inner_cfc_sub_le (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f) (p : Polynomial ℝ)
    (δ : ℝ) (hδ : 0 ≤ δ) (hp : ∀ t ∈ Set.Icc (-(2 + |lam|)) (2 + |lam|), |p.eval t - f t| ≤ δ)
    (η : CfgZ) (x : ℤ) :
    |(⟪deltaZ x, cfc (fun z : ℂ => ((f z.re : ℝ) : ℂ)) (dimerHZ lam η) (deltaZ x)⟫_ℂ).re -
      momZ lam p η x| ≤ δ := by
  have hnormal : IsStarNormal (dimerHZ lam η) := (isSelfAdjoint_dimerHZ lam η).isStarNormal
  have hPev : ∀ t : ℝ, (p.map Complex.ofRealHom).eval (t : ℂ) = ((p.eval t : ℝ) : ℂ) :=
    fun t => by
      rw [Polynomial.eval_map]
      exact Polynomial.eval₂_hom Complex.ofRealHom t
  have hcont : Continuous fun z : ℂ => ((f z.re : ℝ) : ℂ) :=
    Complex.continuous_ofReal.comp (hf.comp Complex.continuous_re)
  have h1 : ‖cfc (fun z : ℂ => ((f z.re : ℝ) : ℂ)) (dimerHZ lam η) -
      Polynomial.aeval (dimerHZ lam η) (p.map Complex.ofRealHom)‖ ≤ δ := by
    rw [← cfc_polynomial (p.map Complex.ofRealHom) (dimerHZ lam η),
      ← cfc_sub (fun z : ℂ => ((f z.re : ℝ) : ℂ)) (fun z => (p.map Complex.ofRealHom).eval z)
        (dimerHZ lam η)]
    refine norm_cfc_le hδ fun z hz => ?_
    obtain ⟨hz1, hz2⟩ := mem_spectrum_dimerHZ lam η hz
    have e : (p.map Complex.ofRealHom).eval z = ((p.eval z.re : ℝ) : ℂ) := by
      rw [← hPev z.re, ← hz1]
    rw [e, ← Complex.ofReal_sub, Complex.norm_real, Real.norm_eq_abs, abs_sub_comm]
    exact hp z.re ⟨by linarith [(abs_le.mp hz2).1], (abs_le.mp hz2).2⟩
  have h2 := (ECInfVolQAux.norm_inner_deltaZ_le _ x x).trans h1
  rw [_root_.sub_apply, inner_sub_right, inner_aeval_deltaZ] at h2
  calc |(⟪deltaZ x, cfc (fun z : ℂ => ((f z.re : ℝ) : ℂ)) (dimerHZ lam η) (deltaZ x)⟫_ℂ).re -
        momZ lam p η x|
      = |(⟪deltaZ x, cfc (fun z : ℂ => ((f z.re : ℝ) : ℂ)) (dimerHZ lam η) (deltaZ x)⟫_ℂ -
          (momZ lam p η x : ℂ)).re| := by
        rw [Complex.sub_re, Complex.ofReal_re]
    _ ≤ ‖⟪deltaZ x, cfc (fun z : ℂ => ((f z.re : ℝ) : ℂ)) (dimerHZ lam η) (deltaZ x)⟫_ℂ -
          (momZ lam p η x : ℂ)‖ := Complex.abs_re_le_norm _
    _ ≤ δ := h2

/-! ### Measurability and integrability -/

/-- `η ↦ (H_η^k)_{xy}` is measurable. -/
lemma measurable_hZpow (J h : ℝ) (k : ℕ) :
    ∀ x y : ℤ, Measurable fun η : CfgZ => hZpow J h (dimerQZ η) k x y := by
  induction k with
  | zero => intro x y; exact measurable_const
  | succ k ih =>
    intro x y
    simp only [hZpow]
    refine Finset.measurable_sum _ fun z _ => Measurable.fun_mul ?_ (ih z y)
    have e : (fun η : CfgZ => hZ J h (dimerQZ η) x z) = fun η =>
        (if x = z + 1 ∨ z = x + 1 then -J else 0) +
          (if x = z then 2 * h * sgn (η ((x + 1) / 2)) else 0) := rfl
    rw [e]
    refine Measurable.fun_add measurable_const ?_
    by_cases hxz : x = z
    · subst hxz
      simp only [ite_true]
      exact Measurable.const_mul ((measurable_of_countable sgn).comp (measurable_pi_apply _)) _
    · simp only [hxz, ite_false]
      exact measurable_const

/-- `η ↦ p(H_η)_{xx}` is measurable. -/
lemma measurable_momZ (lam : ℝ) (p : Polynomial ℝ) (x : ℤ) :
    Measurable fun η => momZ lam p η x :=
  Finset.measurable_sum _ fun k _ => measurable_const.mul (measurable_hZpow _ _ k x x)

/-- The bound `∑_k |p_k| (2 + |λ|)^k` on the matrix elements of `p(H)`. -/
def polyBound (lam : ℝ) (p : Polynomial ℝ) : ℝ :=
  ∑ k ∈ Finset.range (p.natDegree + 1), |p.coeff k| * (2 + |lam|) ^ k

/-- `(2|-1| + 2|λ/2|)^k = (2 + |λ|)^k`. -/
lemma two_abs_add (lam : ℝ) : 2 * |(-1 : ℝ)| + 2 * |lam / 2| = 2 + |lam| := by
  rw [abs_neg, abs_one, abs_div, abs_two]
  ring

/-- `|p(H_η)_{xx}| ≤ ∑_k |p_k| (2 + |λ|)^k`. -/
lemma abs_momZ_le (lam : ℝ) (p : Polynomial ℝ) (η : CfgZ) (x : ℤ) :
    |momZ lam p η x| ≤ polyBound lam p := by
  unfold momZ polyBound
  refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun k _ => ?_)
  rw [abs_mul]
  refine mul_le_mul_of_nonneg_left ?_ (abs_nonneg _)
  have := InfVolAux.abs_hZpow_le (-1) (lam / 2) (dimerQZ η) k x x
  rwa [two_abs_add] at this

/-- `η ↦ p(H_η)_{xx}` is integrable. -/
lemma integrable_momZ (lam : ℝ) (p : Polynomial ℝ) (x : ℤ) :
    Integrable (fun η => momZ lam p η x) muZ :=
  Integrable.of_bound (measurable_momZ lam p x).aestronglyMeasurable (polyBound lam p)
    (ae_of_all _ fun η => by rw [Real.norm_eq_abs]; exact abs_momZ_le lam p η x)

/-- The integrand of `dosFun`. -/
def dosInt (lam : ℝ) (f : ℝ → ℝ) (η : CfgZ) : ℝ :=
  (⟪deltaZ 1, cfc (fun z : ℂ => ((f z.re : ℝ) : ℂ)) (dimerHZ lam η) (deltaZ 1)⟫_ℂ +
    ⟪deltaZ 2, cfc (fun z : ℂ => ((f z.re : ℝ) : ℂ)) (dimerHZ lam η) (deltaZ 2)⟫_ℂ).re

/-- `ν(f) = ½ ∫ dosInt`. -/
lemma dosFun_eq (lam : ℝ) (f : ℝ → ℝ) : dosFun lam f = (1 / 2) * ∫ η, dosInt lam f η ∂muZ :=
  rfl

/-- The integrand is approximated by the polynomial matrix elements. -/
lemma abs_dosInt_sub_le (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f) (p : Polynomial ℝ)
    (δ : ℝ) (hδ : 0 ≤ δ) (hp : ∀ t ∈ Set.Icc (-(2 + |lam|)) (2 + |lam|), |p.eval t - f t| ≤ δ)
    (η : CfgZ) : |dosInt lam f η - (momZ lam p η 1 + momZ lam p η 2)| ≤ 2 * δ := by
  have h1 := abs_re_inner_cfc_sub_le lam f hf p δ hδ hp η 1
  have h2 := abs_re_inner_cfc_sub_le lam f hf p δ hδ hp η 2
  rw [dosInt, Complex.add_re]
  rw [abs_le] at h1 h2 ⊢
  constructor <;> linarith [h1.1, h1.2, h2.1, h2.2]

/-- The integrand of `dosFun` is measurable (a pointwise limit of polynomial matrix elements). -/
lemma measurable_dosInt (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f) :
    Measurable (dosInt lam f) := by
  have hex : ∀ n : ℕ, ∃ p : Polynomial ℝ, ∀ t ∈ Set.Icc (-(2 + |lam|)) (2 + |lam|),
      |p.eval t - f t| < 1 / ((n : ℝ) + 1) := fun n =>
    exists_polynomial_near_of_continuousOn _ _ f hf.continuousOn _ (by positivity)
  choose p hp using hex
  refine measurable_of_tendsto_metrizable
    (f := fun n η => momZ lam (p n) η 1 + momZ lam (p n) η 2)
    (fun n => (measurable_momZ lam (p n) 1).add (measurable_momZ lam (p n) 2)) ?_
  rw [tendsto_pi_nhds]
  intro η
  rw [tendsto_iff_norm_sub_tendsto_zero]
  refine squeeze_zero (g := fun n : ℕ => 2 * (1 / ((n : ℝ) + 1))) (fun n => norm_nonneg _)
    (fun n => ?_) ?_
  · rw [Real.norm_eq_abs, abs_sub_comm]
    exact abs_dosInt_sub_le lam f hf (p n) _ (by positivity) (fun t ht => (hp n t ht).le) η
  · simpa using (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)).const_mul 2

/-- The integrand of `dosFun` is integrable. -/
lemma integrable_dosInt (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f) :
    Integrable (dosInt lam f) muZ := by
  obtain ⟨p, hp⟩ := exists_polynomial_near_of_continuousOn (-(2 + |lam|)) (2 + |lam|) f
    hf.continuousOn 1 one_pos
  refine Integrable.of_bound (measurable_dosInt lam f hf).aestronglyMeasurable
    (2 + 2 * polyBound lam p) (ae_of_all _ fun η => ?_)
  have h1 := abs_dosInt_sub_le lam f hf p 1 zero_le_one (fun t ht => (hp t ht).le) η
  have h2 := abs_momZ_le lam p η 1
  have h3 := abs_momZ_le lam p η 2
  rw [Real.norm_eq_abs, abs_le]
  rw [abs_le] at h1 h2 h3
  constructor <;> linarith [h1.1, h1.2, h2.1, h2.2, h3.1, h3.2]

/-- `|ν(f) - ½ E(p(H)_{11} + p(H)_{22})| ≤ δ` if `|p - f| ≤ δ` on the spectral interval. -/
lemma abs_dosFun_sub_le (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f) (p : Polynomial ℝ)
    (δ : ℝ) (hδ : 0 ≤ δ) (hp : ∀ t ∈ Set.Icc (-(2 + |lam|)) (2 + |lam|), |p.eval t - f t| ≤ δ) :
    |dosFun lam f - (1 / 2) * (∫ η, momZ lam p η 1 ∂muZ + ∫ η, momZ lam p η 2 ∂muZ)| ≤ δ := by
  have hi1 := integrable_momZ lam p 1
  have hi2 := integrable_momZ lam p 2
  have hI := integrable_dosInt lam f hf
  have hi12 : Integrable (fun η => momZ lam p η 1 + momZ lam p η 2) muZ := hi1.add hi2
  have e1 : ∫ η, (dosInt lam f η - (momZ lam p η 1 + momZ lam p η 2)) ∂muZ =
      ∫ η, dosInt lam f η ∂muZ - ∫ η, (momZ lam p η 1 + momZ lam p η 2) ∂muZ :=
    integral_sub hI hi12
  have e2 : ∫ η, (momZ lam p η 1 + momZ lam p η 2) ∂muZ =
      ∫ η, momZ lam p η 1 ∂muZ + ∫ η, momZ lam p η 2 ∂muZ := integral_add hi1 hi2
  have e : dosFun lam f - (1 / 2) * (∫ η, momZ lam p η 1 ∂muZ + ∫ η, momZ lam p η 2 ∂muZ) =
      (1 / 2) * ∫ η, (dosInt lam f η - (momZ lam p η 1 + momZ lam p η 2)) ∂muZ := by
    rw [e1, e2, dosFun_eq]
    ring
  have hb : ‖∫ η, (dosInt lam f η - (momZ lam p η 1 + momZ lam p η 2)) ∂muZ‖ ≤
      2 * δ * muZ.real Set.univ :=
    norm_integral_le_of_norm_le_const (ae_of_all _ fun η => by
      rw [Real.norm_eq_abs]; exact abs_dosInt_sub_le lam f hf p δ hδ hp η)
  rw [probReal_univ, mul_one, Real.norm_eq_abs] at hb
  rw [e, abs_mul, abs_of_pos (by norm_num : (0 : ℝ) < 1 / 2)]
  linarith

end DimerDOSAux

/-- **Convergence of the finite-volume density of states** (paper: `rem:dimerdos`(ii)). -/
theorem tendsto_dosFin (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f) :
    Tendsto (fun M : ℕ => dosFin lam M f) atTop (𝓝 (dosFun lam f)) := by
  sorry

/-- `ν` is positive. -/
theorem dosFun_nonneg (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f) (h0 : ∀ x, 0 ≤ f x) :
    0 ≤ dosFun lam f := by
  sorry

/-- `ν` is normalized. -/
theorem dosFun_one (lam : ℝ) : dosFun lam (fun _ => 1) = 1 := by
  sorry

/-- `ν` lives on the spectral interval `[-2-|λ|, 2+|λ|]`. -/
theorem dosFun_eq_zero_of_support (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f)
    (hsupp : ∀ x ∈ Set.Icc (-2 - |lam|) (2 + |lam|), f x = 0) : dosFun lam f = 0 := by
  sorry

/-- **No concentration near any energy**: `ν(f) ≤ 2/L` for `0 ≤ f ≤ 1` supported in `[E - ε, E + ε]` with
`ε L (|λ| + |E| + 1)^L < 1`. -/
theorem dosFun_le_of_near (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f)
    (h01 : ∀ x, 0 ≤ f x ∧ f x ≤ 1) (E ε : ℝ) (L : ℕ) (hL : 0 < L)
    (hε : ε * L * (|lam| + |E| + 1) ^ L < 1) (hsupp : ∀ x, f x ≠ 0 → |x - E| ≤ ε) :
    dosFun lam f ≤ 2 / L := by
  sorry

/-- **The density of states is a probability measure** (Riesz–Markov–Kakutani). -/
theorem exists_dosMeasure (lam : ℝ) :
    ∃ ν : Measure ℝ, IsProbabilityMeasure ν ∧ ∀ f : CompactlySupportedContinuousMap ℝ ℝ, ∫ x, f x ∂ν = dosFun lam f := by
  sorry

/-- The density of states `ν` of the dimer model. -/
def dosMeasure (lam : ℝ) : Measure ℝ := Classical.choose (exists_dosMeasure lam)

/-- `ν` is a probability measure (paper: `rem:dimerdos`(i)). -/
theorem isProbabilityMeasure_dosMeasure (lam : ℝ) : IsProbabilityMeasure (dosMeasure lam) :=
  (Classical.choose_spec (exists_dosMeasure lam)).1

/-- `∫ f dν = ν(f)` for compactly supported continuous `f`. -/
theorem integral_dosMeasure (lam : ℝ) (f : CompactlySupportedContinuousMap ℝ ℝ) :
    ∫ x, f x ∂(dosMeasure lam) = dosFun lam f :=
  (Classical.choose_spec (exists_dosMeasure lam)).2 f

/-- `ν` is supported in `[-2-|λ|, 2+|λ|]` (paper: `rem:dimerdos`(i)). -/
theorem dosMeasure_compl_Icc (lam : ℝ) :
    dosMeasure lam (Set.Icc (-2 - |lam|) (2 + |lam|))ᶜ = 0 := by
  sorry

/-- **`ν` has no atoms** (paper: `rem:dimerdos`(i)). -/
theorem dosMeasure_singleton (lam E : ℝ) : dosMeasure lam {E} = 0 := by
  sorry

end DFL
