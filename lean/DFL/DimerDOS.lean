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
      have hm : Measurable fun η : CfgZ => sgn (η ((x + 1) / 2)) :=
        (measurable_of_countable sgn).comp (measurable_pi_apply _)
      exact hm.const_mul _
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

/-! ### The dimer law as a marginal of `muZ` -/

/-- The dimer charges `q_i = q(η)_{i+1}` of the chain `Fin (2M)` (Lean site `i` = paper site
`i + 1`). -/
def qOf (M : ℕ) (η : CfgZ) : Cfg (2 * M) := fun i => dimerQZ η ((i.val : ℤ) + 1)

/-- `η ↦ q(η)|_{[1, 2M]}` is measurable. -/
lemma measurable_qOf (M : ℕ) : Measurable (qOf M) :=
  Measurable.of_eval fun _ => measurable_pi_apply _

/-- `q(η)|_{[1, 2M]}` is a dimer configuration. -/
lemma isDimer_qOf (M : ℕ) (η : CfgZ) : IsDimer (qOf M η) := by
  intro k hk hk1
  simp only [qOf, dimerQZ]
  congr 1
  push_cast
  omega

/-- The law of `q(η)|_{[1, 2M]}` under `muZ` is uniform on the dimer configurations. -/
lemma muZ_qOf_preimage (M : ℕ) (q : Cfg (2 * M)) :
    muZ (qOf M ⁻¹' {q}) = if IsDimer q then ((2 : ENNReal) ^ M)⁻¹ else 0 := by
  split_ifs with hq
  · have hset : qOf M ⁻¹' {q} = Set.pi (↑(Finset.Icc (1 : ℤ) M)) (fun j =>
        if h : 2 * (j.toNat - 1) < 2 * M then {q ⟨2 * (j.toNat - 1), h⟩} else Set.univ) := by
      ext η
      simp only [Set.mem_preimage, Set.mem_singleton_iff, Set.mem_pi, Finset.mem_coe,
        Finset.mem_Icc]
      constructor
      · rintro rfl j ⟨hj1, hj2⟩
        rw [dite_eq_left (by omega)]
        simp only [Set.mem_singleton_iff, qOf, dimerQZ]
        congr 1
        omega
      · intro H
        funext i
        have hi := i.isLt
        have h1 := H ((i.val / 2 : ℕ) + 1) ⟨by omega, by omega⟩
        rw [dite_eq_left (by omega)] at h1
        simp only [Set.mem_singleton_iff] at h1
        have e1 : ((i.val : ℤ) + 1 + 1) / 2 = ((i.val / 2 : ℕ) : ℤ) + 1 := by omega
        simp only [qOf, dimerQZ]
        rw [e1, h1]
        by_cases hpar : i.val % 2 = 0
        · congr 1
          ext
          simp only
          omega
        · have := hq ⟨2 * (i.val / 2), by omega⟩ (by simp only; omega) (by simp only; omega)
          have e2 : (⟨2 * ((((i.val / 2 : ℕ) : ℤ) + 1).toNat - 1), by omega⟩ : Fin (2 * M)) =
              ⟨2 * (i.val / 2), by omega⟩ := by
            ext
            simp only
            omega
          rw [e2, this]
          congr 1
          ext
          simp only
          omega
    rw [hset, muZ, Measure.infinitePi_pi _ (fun _ _ => MeasurableSet.of_discrete)]
    rw [Finset.prod_congr rfl (g := fun _ => (2 : ENNReal)⁻¹) (fun j hj => by
      rw [Finset.mem_Icc] at hj
      rw [dite_eq_left (by omega), InfVolAux.unifZ2_singleton])]
    rw [Finset.prod_const, Int.card_Icc, ENNReal.inv_pow]
    congr 2
    omega
  · have : qOf M ⁻¹' {q} = ∅ := by
      ext η
      simp only [Set.mem_preimage, Set.mem_singleton_iff, Set.mem_empty_iff_false, iff_false]
      rintro rfl
      exact hq (isDimer_qOf M η)
    rw [this, measure_empty]

/-- Functions of `q(η)|_{[1, 2M]}` are integrable. -/
lemma integrable_comp_qOf (M : ℕ) (G : Cfg (2 * M) → ℝ) :
    Integrable (fun η => G (qOf M η)) muZ :=
  (integrable_map_measure (measurable_of_countable G).aestronglyMeasurable
    (measurable_qOf M).aemeasurable).mp Integrable.of_finite

/-- **The dimer average is an expectation over `η`**: `E_dim G = E G(q(η)|_{[1, 2M]})`. -/
lemma dimerAvg_eq_integral (M : ℕ) (G : Cfg (2 * M) → ℝ) :
    dimerAvg G = ∫ η, G (qOf M η) ∂muZ := by
  rw [← integral_map (measurable_qOf M).aemeasurable
    (measurable_of_countable G).aestronglyMeasurable, integral_fintype Integrable.of_finite]
  unfold dimerAvg
  rw [card_isDimer ⟨M, by ring⟩, Finset.sum_filter, Finset.sum_div]
  refine Finset.sum_congr rfl fun q _ => ?_
  rw [Measure.real, Measure.map_apply (measurable_qOf M) (measurableSet_singleton q),
    muZ_qOf_preimage, smul_eq_mul]
  have e : 2 * M / 2 = M := by omega
  split_ifs
  · rw [e, ENNReal.toReal_inv, ENNReal.toReal_pow, ENNReal.toReal_ofNat]
    push_cast
    ring
  · simp

/-! ### Finite volume: eigenvalues and traces -/

/-- The finite-volume eigenvalues lie in `[-2-|λ|, 2+|λ|]`. -/
lemma abs_eig_le (lam : ℝ) {n : ℕ} (q : Cfg n) (α : Fin n) :
    |(hq_isHermitian 1 (lam / 2) q).eigenvalues α| ≤ 2 + |lam| := by
  set hB := hq_isHermitian 1 (lam / 2) q with hBdef
  refine ECInfVolQAux.abs_eig_le (hq 1 (lam / 2) q) _ (fun i => ?_)
    (fun j => hB.eigenvectorBasis α j) _ ?_ (NormalizationAux.eig_eq _ hB α)
  · have := InfVolAux.sum_abs_hq_le 1 (lam / 2) q i
    rw [abs_one, abs_div, abs_two] at this
    linarith
  · by_contra hcon
    have h0 : ∀ i, hB.eigenvectorBasis α i = 0 := fun i => by
      by_contra hi
      exact hcon ⟨i, hi⟩
    have := NormalizationAux.eig_orth _ hB α α
    simp [h0] at this

/-- `(Γ B Γ)^k = Γ B^k Γ`. -/
lemma pow_stag_conj {n : ℕ} (B : Matrix (Fin n) (Fin n) ℝ) (k : ℕ) :
    (stag n * B * stag n) ^ k = stag n * B ^ k * stag n := by
  induction k with
  | zero => rw [pow_zero, pow_zero, Matrix.mul_one, stag_mul_stag]
  | succ k ih =>
    rw [pow_succ, ih, pow_succ]
    simp only [Matrix.mul_assoc]
    rw [← Matrix.mul_assoc (stag n) (stag n), stag_mul_stag, Matrix.one_mul]

/-- The diagonal of the powers of `h_q` does not depend on the sign of the hopping. -/
lemma pow_hq_neg_apply_self (J h : ℝ) {n : ℕ} (q : Cfg n) (k : ℕ) (i : Fin n) :
    ((hq (-J) h q) ^ k) i i = ((hq J h q) ^ k) i i := by
  rw [← stag_conj_hq J h q, pow_stag_conj]
  simp only [stag, Matrix.diagonal_mul, Matrix.mul_diagonal]
  linear_combination (((hq J h q) ^ k) i i) * NormalizationAux.neg_one_pow_mul_self i.val

/-- The finite-volume diagonal matrix element `p(A_q)_{ii}`, `A_q = Δ + λ diag(q) = hq (-1) (λ/2) q`. -/
def momFin (lam : ℝ) (p : Polynomial ℝ) {n : ℕ} (q : Cfg n) (i : Fin n) : ℝ :=
  ∑ k ∈ Finset.range (p.natDegree + 1), p.coeff k * ((hq (-1) (lam / 2) q) ^ k) i i

/-- **Trace formula**: `∑_α p(μ_α) = ∑_i p(A_q)_{ii}`, i.e. `Tr p(h_q) = Tr p(Γ h_q Γ)`. -/
lemma sum_eval_eig (lam : ℝ) {n : ℕ} (q : Cfg n) (p : Polynomial ℝ) :
    ∑ α, p.eval ((hq_isHermitian 1 (lam / 2) q).eigenvalues α) = ∑ i, momFin lam p q i := by
  unfold momFin
  set hB := hq_isHermitian 1 (lam / 2) q with hBdef
  have hdiag : ∀ k i, ((hq (-1) (lam / 2) q) ^ k) i i =
      ∑ α, hB.eigenvalues α ^ k * (hB.eigenvectorBasis α i * hB.eigenvectorBasis α i) := by
    intro k i
    rw [pow_hq_neg_apply_self]
    exact ECInfVolQAux.pow_apply_eq_sum _ (fun α => hB.eigenvectorBasis α) _
      (NormalizationAux.eig_compl _ hB) (NormalizationAux.eig_eq _ hB) k i i
  simp_rw [hdiag, Polynomial.eval_eq_sum_range, Finset.mul_sum]
  symm
  calc ∑ i, ∑ k ∈ Finset.range (p.natDegree + 1), ∑ α,
        p.coeff k * (hB.eigenvalues α ^ k * (hB.eigenvectorBasis α i * hB.eigenvectorBasis α i))
      = ∑ k ∈ Finset.range (p.natDegree + 1), ∑ α, ∑ i,
        p.coeff k * (hB.eigenvalues α ^ k *
          (hB.eigenvectorBasis α i * hB.eigenvectorBasis α i)) := by
        rw [Finset.sum_comm]
        exact Finset.sum_congr rfl fun k _ => Finset.sum_comm
    _ = ∑ k ∈ Finset.range (p.natDegree + 1), ∑ α, p.coeff k * hB.eigenvalues α ^ k := by
        refine Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun α _ => ?_
        rw [← Finset.mul_sum, ← Finset.mul_sum, NormalizationAux.eig_orth _ hB α α]
        simp
    _ = ∑ α, ∑ k ∈ Finset.range (p.natDegree + 1), p.coeff k * hB.eigenvalues α ^ k :=
        Finset.sum_comm

/-! ### Finite propagation -/

/-- The entries of `h_q` in a shifted box `q_i = Q_{i + a}` are those of `h^ℤ_Q`. -/
lemma hq_apply_shift (J h : ℝ) {n : ℕ} (q : Cfg n) (Q : CfgZ) (a : ℤ)
    (hqQ : ∀ i : Fin n, q i = Q (i.val + a)) (i l : Fin n) :
    hq J h q i l = hZ J h Q (i.val + a) (l.val + a) := by
  have e1 : (i.val + 1 = l.val ∨ l.val + 1 = i.val) ↔
      ((i.val : ℤ) + a = l.val + a + 1 ∨ (l.val : ℤ) + a = i.val + a + 1) := by omega
  have e2 : i = l ↔ (i.val : ℤ) + a = l.val + a := by rw [Fin.ext_iff]; omega
  simp only [hq, hZ, pathAdj, Matrix.add_apply, Matrix.smul_apply, Matrix.diagonal_apply,
    smul_eq_mul, hqQ]
  by_cases c1 : i.val + 1 = l.val ∨ l.val + 1 = i.val
  · rw [ite_eq_left c1, ite_eq_left (e1.mp c1)]
    by_cases c2 : i = l
    · rw [ite_eq_left c2, ite_eq_left (e2.mp c2)]; ring
    · rw [ite_eq_right c2, ite_eq_right (fun h => c2 (e2.mpr h))]; ring
  · rw [ite_eq_right c1, ite_eq_right (fun h => c1 (e1.mpr h))]
    by_cases c2 : i = l
    · rw [ite_eq_left c2, ite_eq_left (e2.mp c2)]; ring
    · rw [ite_eq_right c2, ite_eq_right (fun h => c2 (e2.mpr h))]; ring

/-- One row of `h_q` in a shifted box, away from the boundary. -/
lemma sum_hq_mul_shift (J h : ℝ) {n : ℕ} (q : Cfg n) (Q : CfgZ) (a : ℤ)
    (hqQ : ∀ i : Fin n, q i = Q (i.val + a)) (i : Fin n) (hi1 : 1 ≤ i.val)
    (hi2 : i.val + 1 < n) (v : Fin n → ℝ) (w : ℤ → ℝ)
    (hvw : ∀ l : Fin n, i.val ≤ l.val + 1 → l.val ≤ i.val + 1 → v l = w (l.val + a)) :
    ∑ l, hq J h q i l * v l =
      ∑ z ∈ Finset.Icc ((i.val : ℤ) + a - 1) (i.val + a + 1), hZ J h Q (i.val + a) z * w z := by
  set S := Finset.univ.filter (fun l : Fin n => i.val ≤ l.val + 1 ∧ l.val ≤ i.val + 1)
    with hS
  rw [← Finset.sum_subset (Finset.subset_univ S) ?_]
  · refine Finset.sum_bij (fun l _ => (l.val : ℤ) + a) ?_ ?_ ?_ ?_
    · intro l hl
      simp only [hS, Finset.mem_filter, Finset.mem_univ, true_and] at hl
      simp only [Finset.mem_Icc]
      omega
    · intro l1 _ l2 _ h12
      exact Fin.ext (by omega)
    · intro z hz
      simp only [Finset.mem_Icc] at hz
      refine ⟨⟨(z - a).toNat, by omega⟩, ?_, ?_⟩
      · simp only [hS, Finset.mem_filter, Finset.mem_univ, true_and]
        omega
      · simp only
        omega
    · intro l hl
      simp only [hS, Finset.mem_filter, Finset.mem_univ, true_and] at hl
      rw [hq_apply_shift J h q Q a hqQ, hvw l hl.1 hl.2]
  · intro l _ hl
    simp only [hS, Finset.mem_filter, Finset.mem_univ, true_and, not_and_or, not_le] at hl
    rw [hq_apply_shift J h q Q a hqQ, InfVolAux.hZ_eq_zero J h Q _ _ ?_, zero_mul]
    simp only [Finset.mem_Icc]
    omega

/-- **Finite propagation in a shifted box**: if `q_i = Q_{i + a}`, the powers of `h_q` agree with
those of `h^ℤ_Q` as long as the paths stay in the box. -/
lemma pow_hq_apply_shift (J h : ℝ) {n : ℕ} (q : Cfg n) (Q : CfgZ) (a : ℤ)
    (hqQ : ∀ i : Fin n, q i = Q (i.val + a)) (j : Fin n) (k : ℕ) :
    ∀ i : Fin n, k ≤ i.val → i.val + k < n →
      ((hq J h q) ^ k) i j = hZpow J h Q k (i.val + a) (j.val + a) := by
  induction k with
  | zero =>
    intro i _ _
    simp only [pow_zero, Matrix.one_apply, hZpow]
    by_cases hij : i = j
    · rw [ite_eq_left hij, ite_eq_left (by rw [hij])]
    · rw [ite_eq_right hij, ite_eq_right (fun h => hij (Fin.ext (by omega)))]
  | succ k ih =>
    intro i hi1 hi2
    rw [pow_succ', Matrix.mul_apply]
    simp only [hZpow]
    exact sum_hq_mul_shift J h q Q a hqQ i (by omega) (by omega) _ _
      fun l hl1 hl2 => ih l (by omega) (by omega)

/-! ### Stationarity of the dimer law -/

/-- The shift `η ↦ η(· + m)` preserves `muZ`. -/
lemma map_shift_muZ (m : ℤ) : muZ.map (fun η : CfgZ => fun j => η (j + m)) = muZ := by
  have h := Measure.infinitePi_map_piCongrLeft (X := fun _ : ℤ => ZMod 2)
    (fun _ : ℤ => unifZ2) (Equiv.subRight m)
  have e : ⇑(MeasurableEquiv.piCongrLeft (fun _ : ℤ => ZMod 2) (Equiv.subRight m)) =
      fun η : CfgZ => fun j => η (j + m) := by
    funext η j
    simp [MeasurableEquiv.coe_piCongrLeft, Equiv.piCongrLeft_apply_eq_cast]
  rw [e] at h
  exact h

/-- Integrals are invariant under the shift. -/
lemma integral_comp_shift (m : ℤ) (G : CfgZ → ℝ) (hG : Measurable G) :
    ∫ η, G (fun j => η (j + m)) ∂muZ = ∫ η, G η ∂muZ := by
  have hmeas : Measurable fun η : CfgZ => fun j => η (j + m) :=
    Measurable.of_eval fun j => measurable_pi_apply _
  rw [← integral_map hmeas.aemeasurable hG.aestronglyMeasurable, map_shift_muZ]

/-- `q(η(· + m))_y = q(η)_{y + 2m}`. -/
lemma dimerQZ_shift (η : CfgZ) (m y : ℤ) :
    dimerQZ (fun j => η (j + m)) y = dimerQZ η (y + 2 * m) := by
  simp only [dimerQZ]
  congr 1
  omega

/-- Translation covariance of `h^ℤ`. -/
lemma hZ_shift (J h : ℝ) (Q : CfgZ) (c x z : ℤ) :
    hZ J h (fun w => Q (w + c)) x z = hZ J h Q (x + c) (z + c) := by
  simp only [hZ]
  by_cases c1 : x = z + 1 ∨ z = x + 1
  · rw [ite_eq_left c1, ite_eq_left (by omega : x + c = z + c + 1 ∨ z + c = x + c + 1)]
    by_cases c2 : x = z
    · rw [ite_eq_left c2, ite_eq_left (by omega : x + c = z + c)]
    · rw [ite_eq_right c2, ite_eq_right (by omega : ¬ (x + c = z + c))]
  · rw [ite_eq_right c1, ite_eq_right (by omega : ¬ (x + c = z + c + 1 ∨ z + c = x + c + 1))]
    by_cases c2 : x = z
    · rw [ite_eq_left c2, ite_eq_left (by omega : x + c = z + c)]
    · rw [ite_eq_right c2, ite_eq_right (by omega : ¬ (x + c = z + c))]

/-- Translation covariance of the powers of `h^ℤ`. -/
lemma hZpow_shift (J h : ℝ) (Q : CfgZ) (c : ℤ) (k : ℕ) :
    ∀ x y : ℤ, hZpow J h (fun w => Q (w + c)) k x y = hZpow J h Q k (x + c) (y + c) := by
  induction k with
  | zero =>
    intro x y
    simp only [hZpow, add_left_inj]
  | succ k ih =>
    intro x y
    simp only [hZpow]
    refine Finset.sum_nbij' (· + c) (· - c) ?_ ?_ ?_ ?_ ?_
    · intro z hz
      simp only [Finset.mem_Icc] at hz ⊢
      omega
    · intro z hz
      simp only [Finset.mem_Icc] at hz ⊢
      omega
    · intro z _
      simp
    · intro z _
      simp
    · intro z _
      rw [hZ_shift, ih]

/-- `p(H_{η(· + m)})_{xx} = p(H_η)_{x + 2m, x + 2m}`. -/
lemma momZ_shift (lam : ℝ) (p : Polynomial ℝ) (η : CfgZ) (m x : ℤ) :
    momZ lam p (fun j => η (j + m)) x = momZ lam p η (x + 2 * m) := by
  have e : dimerQZ (fun j => η (j + m)) = fun w => dimerQZ η (w + 2 * m) :=
    funext fun w => dimerQZ_shift η m w
  unfold momZ
  rw [e]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [hZpow_shift]

/-- **Stationarity**: `E p(H)_{x + 2m, x + 2m} = E p(H)_{xx}`. -/
lemma integral_momZ_shift (lam : ℝ) (p : Polynomial ℝ) (m x : ℤ) :
    ∫ η, momZ lam p η (x + 2 * m) ∂muZ = ∫ η, momZ lam p η x ∂muZ := by
  simp_rw [← momZ_shift]
  exact integral_comp_shift m (fun η => momZ lam p η x) (measurable_momZ lam p x)

/-! ### Convergence for polynomials -/

/-- Expectations of bounded functions are bounded. -/
lemma abs_integral_le {G : CfgZ → ℝ} {C : ℝ} (h : ∀ η, |G η| ≤ C) :
    |∫ η, G η ∂muZ| ≤ C := by
  have := norm_integral_le_of_norm_le_const (μ := muZ) (f := G) (C := C)
    (ae_of_all _ fun η => by rw [Real.norm_eq_abs]; exact h η)
  rwa [probReal_univ, mul_one, Real.norm_eq_abs] at this

/-- `0 ≤ ∑_k |p_k| (2 + |λ|)^k`. -/
lemma polyBound_nonneg (lam : ℝ) (p : Polynomial ℝ) : 0 ≤ polyBound lam p :=
  Finset.sum_nonneg fun k _ => mul_nonneg (abs_nonneg _) (by positivity)

/-- `|p(A_q)_{ii}| ≤ ∑_k |p_k| (2 + |λ|)^k`. -/
lemma abs_momFin_le (lam : ℝ) (p : Polynomial ℝ) {n : ℕ} (q : Cfg n) (i : Fin n) :
    |momFin lam p q i| ≤ polyBound lam p := by
  unfold momFin polyBound
  refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun k _ => ?_)
  rw [abs_mul]
  refine mul_le_mul_of_nonneg_left ?_ (abs_nonneg _)
  have := InfVolAux.abs_pow_hq_le (-1) (lam / 2) q k i i
  rwa [two_abs_add] at this

/-- Away from the boundary, the finite-volume matrix elements are the infinite-volume ones. -/
lemma momFin_qOf_eq (lam : ℝ) (p : Polynomial ℝ) (M : ℕ) (η : CfgZ) (i : Fin (2 * M))
    (h1 : p.natDegree ≤ i.val) (h2 : i.val + p.natDegree < 2 * M) :
    momFin lam p (qOf M η) i = momZ lam p η ((i.val : ℤ) + 1) := by
  unfold momFin momZ
  refine Finset.sum_congr rfl fun k hk => ?_
  rw [Finset.mem_range] at hk
  rw [pow_hq_apply_shift (-1) (lam / 2) (qOf M η) (dimerQZ η) 1 (fun _ => rfl) i k i
    (by omega) (by omega)]

/-- `N⁻¹ E_dim Tr f(h_q)` as an expectation over `η`. -/
lemma dosFin_eq_integral (lam : ℝ) (M : ℕ) (f : ℝ → ℝ) :
    dosFin lam M f = (1 / (2 * M : ℝ)) *
      ∫ η, ∑ α, f ((hq_isHermitian 1 (lam / 2) (qOf M η)).eigenvalues α) ∂muZ := by
  unfold dosFin
  rw [dimerAvg_eq_integral]

/-- The polynomial case: `N⁻¹ E_dim Tr p(h_q) = N⁻¹ ∑_i E p(A_q)_{ii}`. -/
lemma dosFin_poly_eq (lam : ℝ) (p : Polynomial ℝ) (M : ℕ) :
    dosFin lam M (fun t => p.eval t) =
      (1 / (2 * M : ℝ)) * ∑ i : Fin (2 * M), ∫ η, momFin lam p (qOf M η) i ∂muZ := by
  have hi : ∀ i ∈ (Finset.univ : Finset (Fin (2 * M))),
      Integrable (fun η => momFin lam p (qOf M η) i) muZ :=
    fun i _ => integrable_comp_qOf M (fun q => momFin lam p q i)
  rw [dosFin_eq_integral, ← integral_finsetSum _ hi]
  congr 1
  exact integral_congr_ae (ae_of_all _ fun η => sum_eval_eig lam (qOf M η) p)

/-- `∑_{i < 2M} g_i = ∑_{m < M} (g_{2m} + g_{2m+1})`. -/
lemma sum_range_two_mul (g : ℕ → ℝ) (M : ℕ) :
    ∑ i ∈ Finset.range (2 * M), g i = ∑ m ∈ Finset.range M, (g (2 * m) + g (2 * m + 1)) := by
  induction M with
  | zero => simp
  | succ M ih =>
    rw [show 2 * (M + 1) = 2 * M + 1 + 1 by ring, Finset.sum_range_succ, Finset.sum_range_succ,
      ih, Finset.sum_range_succ]
    ring

/-- At most `2d` sites of `Fin n` are within distance `d` of the boundary. -/
lemma sum_boundary_le (n d : ℕ) :
    ∑ i : Fin n, (if d ≤ i.val ∧ i.val + d < n then (0 : ℝ) else 1) ≤ 2 * d := by
  rw [Fin.sum_univ_eq_sum_range (fun i => if d ≤ i ∧ i + d < n then (0 : ℝ) else 1) n]
  have h1 : ((Finset.range n).filter (fun i => i < d)).card ≤ d :=
    (Finset.card_le_card fun i hi => by
      simp only [Finset.mem_filter, Finset.mem_range] at hi
      exact Finset.mem_range.mpr hi.2).trans (Finset.card_range d).le
  have h2 : ((Finset.range n).filter (fun i => n ≤ i + d)).card ≤ d := by
    refine (Finset.card_le_card (t := Finset.Ico (n - d) n) fun i hi => ?_).trans ?_
    · simp only [Finset.mem_filter, Finset.mem_range] at hi
      exact Finset.mem_Ico.mpr ⟨by omega, hi.1⟩
    · rw [Nat.card_Ico]
      omega
  calc ∑ i ∈ Finset.range n, (if d ≤ i ∧ i + d < n then (0 : ℝ) else 1)
      ≤ ∑ i ∈ Finset.range n, ((if i < d then (1 : ℝ) else 0) + (if n ≤ i + d then 1 else 0)) := by
        refine Finset.sum_le_sum fun i _ => ?_
        split_ifs <;> norm_num
        omega
    _ = (((Finset.range n).filter (fun i => i < d)).card : ℝ) +
          (((Finset.range n).filter (fun i => n ≤ i + d)).card : ℝ) := by
        rw [Finset.sum_add_distrib, Finset.sum_boole, Finset.sum_boole]
    _ ≤ d + d := by
        gcongr
    _ = 2 * d := by ring

/-- **Quantitative convergence for polynomials**: `|N⁻¹ E_dim Tr p(h_q) - ν(p)| ≤ 2 d B / M`. -/
lemma abs_dosFin_poly_sub_le (lam : ℝ) (p : Polynomial ℝ) (M : ℕ) (hM : 0 < M) :
    |dosFin lam M (fun t => p.eval t) -
        (1 / 2) * (∫ η, momZ lam p η 1 ∂muZ + ∫ η, momZ lam p η 2 ∂muZ)| ≤
      2 * p.natDegree * polyBound lam p / M := by
  set c₁ := ∫ η, momZ lam p η 1 ∂muZ with hc₁
  set c₂ := ∫ η, momZ lam p η 2 ∂muZ with hc₂
  set a : Fin (2 * M) → ℝ := fun i => ∫ η, momFin lam p (qOf M η) i ∂muZ with ha
  set b : Fin (2 * M) → ℝ := fun i => ∫ η, momZ lam p η ((i.val : ℤ) + 1) ∂muZ with hb
  set B := polyBound lam p with hB
  have hB0 : 0 ≤ B := polyBound_nonneg lam p
  have hsumb : ∑ i, b i = M * (c₁ + c₂) := by
    rw [hb, Fin.sum_univ_eq_sum_range (fun i : ℕ => ∫ η, momZ lam p η ((i : ℤ) + 1) ∂muZ) (2 * M),
      sum_range_two_mul]
    have e1 : ∀ m : ℕ, ∫ η, momZ lam p η (((2 * m : ℕ) : ℤ) + 1) ∂muZ = c₁ := fun m => by
      have := integral_momZ_shift lam p m 1
      rw [show (1 : ℤ) + 2 * (m : ℤ) = ((2 * m : ℕ) : ℤ) + 1 by push_cast; ring] at this
      exact this
    have e2 : ∀ m : ℕ, ∫ η, momZ lam p η (((2 * m + 1 : ℕ) : ℤ) + 1) ∂muZ = c₂ := fun m => by
      have := integral_momZ_shift lam p m 2
      rw [show (2 : ℤ) + 2 * (m : ℤ) = ((2 * m + 1 : ℕ) : ℤ) + 1 by push_cast; ring] at this
      exact this
    simp only [e1, e2, Finset.sum_const, Finset.card_range, nsmul_eq_mul]
  have hab : ∀ i : Fin (2 * M), |a i - b i| ≤
      2 * B * (if p.natDegree ≤ i.val ∧ i.val + p.natDegree < 2 * M then 0 else 1) := by
    intro i
    split_ifs with hi
    · have : a i = b i := integral_congr_ae (ae_of_all _ fun η =>
        momFin_qOf_eq lam p M η i hi.1 hi.2)
      rw [this, sub_self, abs_zero, mul_zero]
    · have h1 : |a i| ≤ B := abs_integral_le fun η => abs_momFin_le lam p _ i
      have h2 : |b i| ≤ B := abs_integral_le fun η => abs_momZ_le lam p η _
      rw [mul_one]
      calc |a i - b i| ≤ |a i| + |b i| := abs_sub _ _
        _ ≤ 2 * B := by linarith
  have hM' : (0 : ℝ) < M := by exact_mod_cast hM
  have e : dosFin lam M (fun t => p.eval t) - (1 / 2) * (c₁ + c₂) =
      (1 / (2 * M : ℝ)) * ∑ i, (a i - b i) := by
    rw [dosFin_poly_eq, Finset.sum_sub_distrib, hsumb]
    field_simp
    ring
  rw [e, abs_mul, abs_of_pos (by positivity)]
  have h3 : |∑ i, (a i - b i)| ≤ 2 * B * (2 * p.natDegree) := by
    calc |∑ i, (a i - b i)| ≤ ∑ i, |a i - b i| := Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ i : Fin (2 * M), 2 * B *
          (if p.natDegree ≤ i.val ∧ i.val + p.natDegree < 2 * M then 0 else 1) :=
          Finset.sum_le_sum fun i _ => hab i
      _ = 2 * B * ∑ i : Fin (2 * M),
          (if p.natDegree ≤ i.val ∧ i.val + p.natDegree < 2 * M then (0 : ℝ) else 1) := by
          rw [Finset.mul_sum]
      _ ≤ 2 * B * (2 * p.natDegree) := by
          gcongr
          exact sum_boundary_le _ _
  calc 1 / (2 * M : ℝ) * |∑ i, (a i - b i)| ≤ 1 / (2 * M : ℝ) * (2 * B * (2 * p.natDegree)) := by
        gcongr
    _ = 2 * p.natDegree * B / M := by
        field_simp

/-- **Convergence for polynomials** (finite propagation and stationarity). -/
lemma tendsto_dosFin_poly (lam : ℝ) (p : Polynomial ℝ) :
    Tendsto (fun M : ℕ => dosFin lam M (fun t => p.eval t)) atTop
      (𝓝 ((1 / 2) * (∫ η, momZ lam p η 1 ∂muZ + ∫ η, momZ lam p η 2 ∂muZ))) := by
  rw [tendsto_iff_norm_sub_tendsto_zero]
  refine squeeze_zero' (Eventually.of_forall fun M => norm_nonneg _)
    (eventually_atTop.mpr ⟨1, fun M hM => ?_⟩)
    (tendsto_const_div_atTop_nhds_zero_nat (2 * p.natDegree * polyBound lam p))
  rw [Real.norm_eq_abs]
  exact abs_dosFin_poly_sub_le lam p M hM

/-- `|N⁻¹ E_dim Tr f(h_q) - N⁻¹ E_dim Tr p(h_q)| ≤ δ` if `|p - f| ≤ δ` on the spectral interval. -/
lemma abs_dosFin_sub_le (lam : ℝ) (M : ℕ) (hM : 0 < M) (f : ℝ → ℝ) (p : Polynomial ℝ) (δ : ℝ)
    (hp : ∀ t ∈ Set.Icc (-(2 + |lam|)) (2 + |lam|), |p.eval t - f t| ≤ δ) :
    |dosFin lam M f - dosFin lam M (fun t => p.eval t)| ≤ δ := by
  have hi1 : Integrable (fun η =>
      ∑ α, f ((hq_isHermitian 1 (lam / 2) (qOf M η)).eigenvalues α)) muZ :=
    integrable_comp_qOf M (fun q => ∑ α, f ((hq_isHermitian 1 (lam / 2) q).eigenvalues α))
  have hi2 : Integrable (fun η =>
      ∑ α, (fun t => p.eval t) ((hq_isHermitian 1 (lam / 2) (qOf M η)).eigenvalues α)) muZ :=
    integrable_comp_qOf M
      (fun q => ∑ α, (fun t => p.eval t) ((hq_isHermitian 1 (lam / 2) q).eigenvalues α))
  rw [dosFin_eq_integral, dosFin_eq_integral, ← mul_sub, ← integral_sub hi1 hi2]
  have hM' : (0 : ℝ) < M := by exact_mod_cast hM
  have hb : ∀ η, |∑ α, f ((hq_isHermitian 1 (lam / 2) (qOf M η)).eigenvalues α) -
      ∑ α, (fun t => p.eval t) ((hq_isHermitian 1 (lam / 2) (qOf M η)).eigenvalues α)| ≤
        2 * M * δ := by
    intro η
    rw [← Finset.sum_sub_distrib]
    refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
    calc ∑ α, |f ((hq_isHermitian 1 (lam / 2) (qOf M η)).eigenvalues α) -
          (fun t => p.eval t) ((hq_isHermitian 1 (lam / 2) (qOf M η)).eigenvalues α)|
        ≤ ∑ α : Fin (2 * M), δ := Finset.sum_le_sum fun α _ => by
          have hα := abs_le.mp (abs_eig_le lam (qOf M η) α)
          rw [abs_sub_comm]
          exact hp _ ⟨by linarith [hα.1], hα.2⟩
      _ = 2 * M * δ := by simp
  rw [abs_mul, abs_of_pos (by positivity)]
  calc 1 / (2 * M : ℝ) * |∫ η, (∑ α, f ((hq_isHermitian 1 (lam / 2) (qOf M η)).eigenvalues α) -
        ∑ α, (fun t => p.eval t) ((hq_isHermitian 1 (lam / 2) (qOf M η)).eigenvalues α)) ∂muZ|
      ≤ 1 / (2 * M : ℝ) * (2 * M * δ) := by
        gcongr
        exact abs_integral_le hb
    _ = δ := by
        field_simp

end DimerDOSAux

/-- **Convergence of the finite-volume density of states** (paper: `rem:dimerdos`(ii)). -/
theorem tendsto_dosFin (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f) :
    Tendsto (fun M : ℕ => dosFin lam M f) atTop (𝓝 (dosFun lam f)) := by
  rw [Metric.tendsto_atTop]
  intro ε hε
  obtain ⟨p, hp⟩ := exists_polynomial_near_of_continuousOn (-(2 + |lam|)) (2 + |lam|) f
    hf.continuousOn (ε / 4) (by positivity)
  have hp' : ∀ t ∈ Set.Icc (-(2 + |lam|)) (2 + |lam|), |p.eval t - f t| ≤ ε / 4 :=
    fun t ht => (hp t ht).le
  have h1 := DimerDOSAux.abs_dosFun_sub_le lam f hf p (ε / 4) (by positivity) hp'
  obtain ⟨N, hN⟩ := Metric.tendsto_atTop.mp (DimerDOSAux.tendsto_dosFin_poly lam p) (ε / 4)
    (by positivity)
  refine ⟨max N 1, fun M hM => ?_⟩
  have h3 := DimerDOSAux.abs_dosFin_sub_le lam M
    (Nat.lt_of_lt_of_le Nat.one_pos (le_of_max_le_right hM)) f p (ε / 4) hp'
  have h4 := hN M (le_of_max_le_left hM)
  rw [Real.dist_eq] at h4 ⊢
  rw [abs_le] at h1 h3
  rw [abs_lt] at h4 ⊢
  constructor <;> linarith [h1.1, h1.2, h3.1, h3.2, h4.1, h4.2]

namespace DimerDOSAux

/-! ### Consequences of the convergence -/

/-- `N⁻¹ E_dim Tr c = c`. -/
lemma dosFin_const (lam : ℝ) (M : ℕ) (hM : 0 < M) (c : ℝ) : dosFin lam M (fun _ => c) = c := by
  have hM' : (0 : ℝ) < M := by exact_mod_cast hM
  rw [dosFin_eq_integral]
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, integral_const,
    probReal_univ, one_smul]
  push_cast
  field_simp

/-- `N⁻¹ E_dim Tr f(h_q)` is additive in `f`. -/
lemma dosFin_add (lam : ℝ) (M : ℕ) (f g : ℝ → ℝ) :
    dosFin lam M (fun x => f x + g x) = dosFin lam M f + dosFin lam M g := by
  unfold dosFin dimerAvg
  simp only [Finset.sum_add_distrib, add_div, mul_add]

/-- `N⁻¹ E_dim Tr f(h_q)` is homogeneous in `f`. -/
lemma dosFin_smul (lam : ℝ) (M : ℕ) (c : ℝ) (f : ℝ → ℝ) :
    dosFin lam M (fun x => c * f x) = c * dosFin lam M f := by
  unfold dosFin dimerAvg
  simp only [← Finset.mul_sum, mul_div_assoc]
  ring

/-- `ν` is additive. -/
lemma dosFun_add (lam : ℝ) (f g : ℝ → ℝ) (hf : Continuous f) (hg : Continuous g) :
    dosFun lam (fun x => f x + g x) = dosFun lam f + dosFun lam g := by
  refine tendsto_nhds_unique (tendsto_dosFin lam _ (hf.add hg)) ?_
  simp_rw [dosFin_add]
  exact (tendsto_dosFin lam f hf).add (tendsto_dosFin lam g hg)

/-- `ν` is homogeneous. -/
lemma dosFun_smul (lam : ℝ) (c : ℝ) (f : ℝ → ℝ) (hf : Continuous f) :
    dosFun lam (fun x => c * f x) = c * dosFun lam f := by
  refine tendsto_nhds_unique (tendsto_dosFin lam _ (continuous_const.mul hf)) ?_
  simp_rw [dosFin_smul]
  exact (tendsto_dosFin lam f hf).const_mul c

/-- `ν` is monotone. -/
lemma dosFun_mono (lam : ℝ) (f g : ℝ → ℝ) (hf : Continuous f) (hg : Continuous g)
    (hfg : ∀ x, f x ≤ g x) : dosFun lam f ≤ dosFun lam g := by
  refine le_of_tendsto_of_tendsto' (tendsto_dosFin lam f hf) (tendsto_dosFin lam g hg)
    fun M => ?_
  unfold dosFin
  exact mul_le_mul_of_nonneg_left (DimerLGTAux.dimerAvg_mono fun q _ =>
    Finset.sum_le_sum fun α _ => hfg _) (by positivity)

/-- `ν(f)` only depends on `f` on the spectral interval. -/
lemma dosFun_congr (lam : ℝ) (f g : ℝ → ℝ) (hf : Continuous f) (hg : Continuous g)
    (hfg : ∀ x ∈ Set.Icc (-(2 + |lam|)) (2 + |lam|), f x = g x) :
    dosFun lam f = dosFun lam g := by
  refine tendsto_nhds_unique (tendsto_dosFin lam f hf) ?_
  have e : ∀ M, dosFin lam M f = dosFin lam M g := fun M => by
    unfold dosFin
    congr 2
    funext q
    exact Finset.sum_congr rfl fun α _ => hfg _ (abs_le.mp (abs_eig_le lam q α))
  simp_rw [e]
  exact tendsto_dosFin lam g hg

/-- `ν(K) ≤ ∫ f dν` for a finite measure and `f ≥ 0` with `f = 1` on `K`. -/
lemma measureReal_le_integral (ν : Measure ℝ) [IsFiniteMeasure ν] {K : Set ℝ}
    (hK : MeasurableSet K) (f : CompactlySupportedContinuousMap ℝ ℝ) (h0 : ∀ x, 0 ≤ f x)
    (h1 : ∀ x ∈ K, f x = 1) : ν.real K ≤ ∫ x, f x ∂ν := by
  rw [← integral_indicator_one hK]
  refine integral_mono ((integrable_const (1 : ℝ)).indicator hK) f.integrable fun x => ?_
  by_cases hx : x ∈ K
  · simp [Set.indicator_of_mem hx, h1 x hx]
  · rw [Set.indicator_of_notMem hx]
    exact h0 x

end DimerDOSAux

/-- `ν` is positive. -/
theorem dosFun_nonneg (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f) (h0 : ∀ x, 0 ≤ f x) :
    0 ≤ dosFun lam f := by
  refine ge_of_tendsto' (tendsto_dosFin lam f hf) fun M => ?_
  unfold dosFin
  exact mul_nonneg (by positivity)
    (DimerLGTAux.dimerAvg_nonneg fun q => Finset.sum_nonneg fun α _ => h0 _)

/-- `ν` is normalized. -/
theorem dosFun_one (lam : ℝ) : dosFun lam (fun _ => 1) = 1 := by
  refine tendsto_nhds_unique (tendsto_dosFin lam _ continuous_const) ?_
  exact tendsto_const_nhds.congr'
    (eventually_atTop.mpr ⟨1, fun M hM => (DimerDOSAux.dosFin_const lam M hM 1).symm⟩)

/-- `ν` lives on the spectral interval `[-2-|λ|, 2+|λ|]`. -/
theorem dosFun_eq_zero_of_support (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f)
    (hsupp : ∀ x ∈ Set.Icc (-2 - |lam|) (2 + |lam|), f x = 0) : dosFun lam f = 0 := by
  have e := DimerDOSAux.dosFun_congr lam f (fun _ => 0) hf continuous_const fun x hx =>
    hsupp x ⟨by linarith [hx.1], hx.2⟩
  rw [e]
  refine tendsto_nhds_unique (tendsto_dosFin lam _ continuous_const) ?_
  exact tendsto_const_nhds.congr'
    (eventually_atTop.mpr ⟨1, fun M hM => (DimerDOSAux.dosFin_const lam M hM 0).symm⟩)

/-- **No concentration near any energy**: `ν(f) ≤ 2/L` for `0 ≤ f ≤ 1` supported in `[E - ε, E + ε]` with
`ε L (|λ| + |E| + 1)^L < 1`. -/
theorem dosFun_le_of_near (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f)
    (h01 : ∀ x, 0 ≤ f x ∧ f x ≤ 1) (E ε : ℝ) (L : ℕ) (hL : 0 < L)
    (hε : ε * L * (|lam| + |E| + 1) ^ L < 1) (hsupp : ∀ x, f x ≠ 0 → |x - E| ≤ ε) :
    dosFun lam f ≤ 2 / L := by
  have hε' : ε * L * ((2 * |lam / 2| + |E|) / |(1 : ℝ)| + 1) ^ L < |(1 : ℝ)| := by
    rw [abs_one, div_one, abs_div, abs_two, mul_div_cancel₀ _ two_ne_zero]
    exact hε
  have hbound : ∀ M : ℕ, 0 < M → dosFin lam M f ≤ 2 / L + 1 / M := by
    intro M hM
    have hM' : (0 : ℝ) < M := by exact_mod_cast hM
    have hL' : (0 : ℝ) < L := by exact_mod_cast hL
    have hq : ∀ q : Cfg (2 * M), ∑ α, f ((hq_isHermitian 1 (lam / 2) q).eigenvalues α) ≤
        2 * ((2 * M : ℝ) / L + 1) := by
      intro q
      have hc := card_eig_near_le 1 (lam / 2) one_ne_zero q E ε L hL hε'
      push_cast at hc
      calc ∑ α, f ((hq_isHermitian 1 (lam / 2) q).eigenvalues α)
          ≤ ∑ α, (if |(hq_isHermitian 1 (lam / 2) q).eigenvalues α - E| ≤ ε then (1 : ℝ)
              else 0) := Finset.sum_le_sum fun α _ => by
            split_ifs with h
            · exact (h01 _).2
            · have : f ((hq_isHermitian 1 (lam / 2) q).eigenvalues α) = 0 := by
                by_contra hne
                exact h (hsupp _ hne)
              rw [this]
        _ = ((Finset.univ.filter (fun α : Fin (2 * M) =>
              |(hq_isHermitian 1 (lam / 2) q).eigenvalues α - E| ≤ ε)).card : ℝ) := by
            rw [Finset.sum_boole]
        _ ≤ 2 * ((2 * M : ℝ) / L + 1) := hc
    have hi : Integrable (fun η =>
        ∑ α, f ((hq_isHermitian 1 (lam / 2) (DimerDOSAux.qOf M η)).eigenvalues α)) muZ :=
      DimerDOSAux.integrable_comp_qOf M
        (fun q => ∑ α, f ((hq_isHermitian 1 (lam / 2) q).eigenvalues α))
    rw [DimerDOSAux.dosFin_eq_integral]
    calc 1 / (2 * M : ℝ) *
          ∫ η, ∑ α, f ((hq_isHermitian 1 (lam / 2) (DimerDOSAux.qOf M η)).eigenvalues α) ∂muZ
        ≤ 1 / (2 * M : ℝ) * (2 * ((2 * M : ℝ) / L + 1)) := by
          gcongr
          calc ∫ η, ∑ α, f ((hq_isHermitian 1 (lam / 2) (DimerDOSAux.qOf M η)).eigenvalues α) ∂muZ
              ≤ ∫ η, 2 * ((2 * M : ℝ) / L + 1) ∂muZ :=
                integral_mono hi (integrable_const _) fun η => hq _
            _ = 2 * ((2 * M : ℝ) / L + 1) := by simp
      _ = 2 / L + 1 / M := by
          field_simp
  refine le_of_tendsto_of_tendsto (tendsto_dosFin lam f hf) ?_
    (eventually_atTop.mpr ⟨1, fun M hM => hbound M hM⟩)
  simpa using (tendsto_const_nhds (x := (2 / L : ℝ))).add
    (tendsto_one_div_atTop_nhds_zero_nat (𝕜 := ℝ))

namespace DimerDOSAux

/-- The positive linear functional `f ↦ ν(f)` on `C_c(ℝ, ℝ)`. -/
def dosPLM (lam : ℝ) : CompactlySupportedContinuousMap ℝ ℝ →ₚ[ℝ] ℝ where
  toFun f := dosFun lam f
  map_add' f g := by
    rw [CompactlySupportedContinuousMap.coe_add]
    exact dosFun_add lam f g f.continuous g.continuous
  map_smul' c f := by
    rw [CompactlySupportedContinuousMap.coe_smul, RingHom.id_apply, smul_eq_mul]
    exact dosFun_smul lam c f f.continuous
  monotone' f g hfg := dosFun_mono lam f g f.continuous g.continuous
    (CompactlySupportedContinuousMap.le_def.mp hfg)

end DimerDOSAux

/-- **The density of states is a probability measure** (Riesz–Markov–Kakutani). -/
theorem exists_dosMeasure (lam : ℝ) :
    ∃ ν : Measure ℝ, IsProbabilityMeasure ν ∧ ∀ f : CompactlySupportedContinuousMap ℝ ℝ, ∫ x, f x ∂ν = dosFun lam f := by
  set Λ := DimerDOSAux.dosPLM lam with hΛ
  have hΛf : ∀ f, Λ f = dosFun lam f := fun f => rfl
  refine ⟨RealRMK.rieszMeasure Λ, ⟨le_antisymm ?_ ?_⟩,
    fun f => (RealRMK.integral_rieszMeasure Λ f).trans (hΛf f)⟩
  · rw [isOpen_univ.measure_eq_iSup_isCompact]
    refine iSup_le fun K => iSup_le fun _ => iSup_le fun hK => ?_
    obtain ⟨g, hg1, hg2, -, hg4⟩ :=
      exists_continuousMap_one_of_isCompact_subset_isOpen hK isOpen_univ (Set.subset_univ K)
    refine (RealRMK.rieszMeasure_le_of_eq_one Λ (f := ⟨g, hg2⟩) (fun x => (hg4 x).1) hK
      fun x hx => hg1 hx).trans ?_
    rw [hΛf, ENNReal.ofReal_le_one]
    calc dosFun lam g ≤ dosFun lam (fun _ => 1) :=
          DimerDOSAux.dosFun_mono lam g _ g.continuous continuous_const fun x => (hg4 x).2
      _ = 1 := dosFun_one lam
  · obtain ⟨g, hg1, hg2, -, hg4⟩ := exists_continuousMap_one_of_isCompact_subset_isOpen
      (isCompact_Icc (a := -(2 + |lam|)) (b := 2 + |lam|)) isOpen_univ (Set.subset_univ _)
    have h1 : dosFun lam g = 1 := (DimerDOSAux.dosFun_congr lam g (fun _ => 1) g.continuous
      continuous_const fun x hx => hg1 hx).trans (dosFun_one lam)
    calc (1 : ENNReal) = ENNReal.ofReal (Λ ⟨g, hg2⟩) := by
          show (1 : ENNReal) = ENNReal.ofReal (dosFun lam g)
          rw [h1, ENNReal.ofReal_one]
      _ ≤ RealRMK.rieszMeasure Λ Set.univ :=
          RealRMK.le_rieszMeasure_tsupport_subset Λ (f := ⟨g, hg2⟩)
            (fun x => ⟨(hg4 x).1, (hg4 x).2⟩) (Set.subset_univ _)

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
  have := isProbabilityMeasure_dosMeasure lam
  have hU : IsOpen (Set.Icc (-2 - |lam|) (2 + |lam|))ᶜ := isOpen_compl_iff.mpr isClosed_Icc
  rw [hU.measure_eq_iSup_isCompact]
  refine le_antisymm (iSup_le fun K => iSup_le fun hKU => iSup_le fun hK => ?_) zero_le
  obtain ⟨g, hg1, hg2, hg3, hg4⟩ := exists_continuousMap_one_of_isCompact_subset_isOpen hK hU hKU
  have hg0 : dosFun lam g = 0 := dosFun_eq_zero_of_support lam g g.continuous fun x hx =>
    image_eq_zero_of_notMem_tsupport fun hx' => hg3 hx' hx
  have h := DimerDOSAux.measureReal_le_integral (dosMeasure lam) hK.measurableSet ⟨g, hg2⟩
    (fun x => (hg4 x).1) (fun x hx => hg1 hx)
  rw [integral_dosMeasure] at h
  have h' : (dosMeasure lam).real K = 0 :=
    le_antisymm (h.trans (le_of_eq hg0)) measureReal_nonneg
  rw [measureReal_eq_zero_iff (measure_ne_top _ _)] at h'
  rw [h']

/-- **`ν` has no atoms** (paper: `rem:dimerdos`(i)). -/
theorem dosMeasure_singleton (lam E : ℝ) : dosMeasure lam {E} = 0 := by
  have := isProbabilityMeasure_dosMeasure lam
  have hbound : ∀ L : ℕ, 0 < L → (dosMeasure lam).real {E} ≤ 2 / L := by
    intro L hL
    have hL' : (0 : ℝ) < L := by exact_mod_cast hL
    have hX : 0 < (|lam| + |E| + 1) ^ L := by positivity
    set ε : ℝ := 1 / (2 * L * (|lam| + |E| + 1) ^ L) with hεdef
    have hε0 : 0 < ε := by positivity
    have hε : ε * L * (|lam| + |E| + 1) ^ L < 1 := by
      have e : ε * L * (|lam| + |E| + 1) ^ L = 1 / 2 := by
        rw [hεdef]
        field_simp
      rw [e]
      norm_num
    have hsub : ({E} : Set ℝ) ⊆ Set.Ioo (E - ε) (E + ε) := by
      intro x hx
      rw [Set.mem_singleton_iff] at hx
      subst hx
      exact ⟨by linarith, by linarith⟩
    obtain ⟨g, hg1, hg2, hg3, hg4⟩ :=
      exists_continuousMap_one_of_isCompact_subset_isOpen isCompact_singleton isOpen_Ioo hsub
    have hsupp : ∀ x, g x ≠ 0 → |x - E| ≤ ε := fun x hx => by
      have := hg3 (subset_tsupport _ hx)
      rw [abs_le]
      constructor <;> linarith [this.1, this.2]
    have h1 := dosFun_le_of_near lam g g.continuous hg4 E ε L hL hε hsupp
    have h2 := DimerDOSAux.measureReal_le_integral (dosMeasure lam) (measurableSet_singleton E)
      ⟨g, hg2⟩ (fun x => (hg4 x).1) (fun x hx => hg1 hx)
    rw [integral_dosMeasure] at h2
    exact h2.trans h1
  have h0 : (dosMeasure lam).real {E} ≤ 0 :=
    le_of_tendsto_of_tendsto tendsto_const_nhds (tendsto_const_div_atTop_nhds_zero_nat 2)
      (eventually_atTop.mpr ⟨1, fun L hL => hbound L hL⟩)
  have h' : (dosMeasure lam).real {E} = 0 := le_antisymm h0 measureReal_nonneg
  rwa [measureReal_eq_zero_iff (measure_ne_top _ _)] at h'

end DFL
