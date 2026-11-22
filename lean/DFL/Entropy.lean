import Mathlib

/-!
# Entropy inequalities

Shannon entropy of distributions on `{0, …, M}`, the support bound, Massey's maximum-entropy bound
for integer-valued random variables (paper: proof of Theorem `thm:numberentropy`), the variance of a
mixture, and the binary-entropy bound `h(ν) ≤ 2 √(ν(1-ν))` used for the area law.
-/

open Finset Real

noncomputable section

namespace DFL

/-- Shannon entropy (in nats) of a distribution `p` on `{0, …, M}`. -/
def shannon (M : ℕ) (p : ℕ → ℝ) : ℝ := ∑ k ∈ range (M + 1), negMulLog (p k)

/-- Mean of a distribution on `{0, …, M}`. -/
def pmean (M : ℕ) (p : ℕ → ℝ) : ℝ := ∑ k ∈ range (M + 1), p k * k

/-- Variance of a distribution on `{0, …, M}`. -/
def pvar (M : ℕ) (p : ℕ → ℝ) : ℝ := ∑ k ∈ range (M + 1), p k * ((k : ℝ) - pmean M p) ^ 2

/-- `p` is a probability distribution on `{0, …, M}`. -/
structure IsProb (M : ℕ) (p : ℕ → ℝ) : Prop where
  nonneg : ∀ k, 0 ≤ p k
  sum_one : ∑ k ∈ range (M + 1), p k = 1

variable {M : ℕ} {p : ℕ → ℝ}

/-- Each probability of a distribution on `{0, …, M}` is at most `1`. -/
private lemma IsProb.le_one (hp : IsProb M p) {k : ℕ} (hk : k ∈ range (M + 1)) : p k ≤ 1 := by
  rw [← hp.sum_one]
  exact Finset.single_le_sum (fun i _ => hp.nonneg i) hk

theorem shannon_nonneg (hp : IsProb M p) : 0 ≤ shannon M p := by
  exact Finset.sum_nonneg fun k hk => negMulLog_nonneg (hp.nonneg k) (hp.le_one hk)

/-- Termwise Gibbs inequality: `-x log x ≤ w - x - x log w` for `x ≥ 0`, `w > 0`. -/
private lemma negMulLog_le_gibbs {x w : ℝ} (hx : 0 ≤ x) (hw : 0 < w) :
    negMulLog x ≤ w - x - x * Real.log w := by
  rcases hx.eq_or_lt with rfl | hx
  · simp only [negMulLog_zero, sub_zero, zero_mul]
    exact hw.le
  · have h := Real.log_le_sub_one_of_pos (div_pos hw hx)
    rw [Real.log_div hw.ne' hx.ne'] at h
    have h2 : x * (Real.log w - Real.log x) ≤ x * (w / x - 1) :=
      mul_le_mul_of_nonneg_left h hx.le
    have h3 : x * (w / x - 1) = w - x := by field_simp
    rw [negMulLog]
    nlinarith

/-- The entropy is at most the logarithm of the size of the support. -/
theorem shannon_le_log_card (hp : IsProb M p) (s : Finset ℕ) (hs : s.Nonempty)
    (hsupp : ∀ k ∈ range (M + 1), k ∉ s → p k = 0) :
    shannon M p ≤ Real.log s.card := by
  have hcard : (0 : ℝ) < s.card := by exact_mod_cast hs.card_pos
  have key : ∀ k ∈ range (M + 1), negMulLog (p k) ≤
      (if k ∈ s then (s.card : ℝ)⁻¹ else 0) - p k + p k * Real.log s.card := by
    intro k hk
    by_cases hks : k ∈ s
    · simp only [hks, ↓reduceIte]
      have := negMulLog_le_gibbs (hp.nonneg k) (inv_pos.mpr hcard)
      rw [Real.log_inv] at this
      linarith
    · simp only [hks, ↓reduceIte, hsupp k hk hks]
      simp
  have hind : ∑ k ∈ range (M + 1), (if k ∈ s then (s.card : ℝ)⁻¹ else 0) ≤ 1 := by
    rw [Finset.sum_ite_mem, Finset.sum_const, nsmul_eq_mul]
    calc ((range (M + 1) ∩ s).card : ℝ) * (s.card : ℝ)⁻¹
        ≤ s.card * (s.card : ℝ)⁻¹ := by
          gcongr
          exact Finset.inter_subset_right
      _ = 1 := mul_inv_cancel₀ hcard.ne'
  calc shannon M p ≤ ∑ k ∈ range (M + 1),
        ((if k ∈ s then (s.card : ℝ)⁻¹ else 0) - p k + p k * Real.log s.card) :=
        Finset.sum_le_sum key
    _ = ∑ k ∈ range (M + 1), (if k ∈ s then (s.card : ℝ)⁻¹ else 0) - 1 +
          Real.log s.card := by
        rw [Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.sum_mul, hp.sum_one]
        ring
    _ ≤ Real.log s.card := by linarith

/-- The variance is nonnegative. -/
private lemma pvar_nonneg (hp : IsProb M p) : 0 ≤ pvar M p :=
  Finset.sum_nonneg fun k _ => mul_nonneg (hp.nonneg k) (sq_nonneg _)

/-- `∫_k^{k+1} (A - (y - μ)² / (2v)) dy = A - ((k + ½ - μ)² + 1/12) / (2v)`. -/
private lemma integral_tangent_quad (k μ v A : ℝ) (hv : v ≠ 0) :
    ∫ y in k..(k + 1), (A - (y - μ) ^ 2 / (2 * v)) =
      A - ((k + 1 / 2 - μ) ^ 2 + 1 / 12) / (2 * v) := by
  have hderiv : ∀ y ∈ Set.uIcc k (k + 1),
      HasDerivAt (fun y => A * y - (y - μ) ^ 3 / (6 * v)) (A - (y - μ) ^ 2 / (2 * v)) y := by
    intro y _
    have h1 : HasDerivAt (fun y => A * y) A y := by
      simpa using (hasDerivAt_id y).const_mul A
    have h2 : HasDerivAt (fun y => (y - μ) ^ 3 / (6 * v)) ((y - μ) ^ 2 / (2 * v)) y := by
      have := (((hasDerivAt_id' y).sub_const μ).fun_pow 3).div_const (6 * v)
      convert this using 1
      norm_num
      field_simp
      ring
    exact h1.sub h2
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hderiv
    (by apply Continuous.intervalIntegrable; fun_prop)]
  field_simp
  ring

/-- **Massey's bound**: the Shannon entropy of an integer-valued random variable with variance `σ²`
is at most `½ ln(2πe(σ² + 1/12))`. -/
theorem massey (hp : IsProb M p) :
    shannon M p ≤ (1 / 2) * Real.log (2 * π * exp 1 * (pvar M p + 1 / 12)) := by
  -- Proof: compare `p` with the weights `w k = ∫_k^{k+1} g` of the Gaussian density `g` of mean
  -- `m + 1/2` and variance `v = σ² + 1/12`: `∑ w k ≤ ∫ g = 1`, discrete Gibbs gives
  -- `H(p) ≤ -∑ p k log (w k)`, and Jensen (tangent line of `exp`) gives
  -- `log (w k) ≥ ∫_k^{k+1} log g = log c - ((k - m)² + 1/12) / (2v)`, `c = (2πv)^{-1/2}`.
  -- notation: mean `m`, smoothed variance `v = σ² + 1/12`, Gaussian centre `μ = m + 1/2`
  obtain ⟨m, hm⟩ : ∃ m, m = pmean M p := ⟨_, rfl⟩
  obtain ⟨v, hv⟩ : ∃ v, v = pvar M p + 1 / 12 := ⟨_, rfl⟩
  have hv0 : 0 < v := by have := pvar_nonneg hp; rw [hv]; linarith
  rw [← hv]
  obtain ⟨c, hc⟩ : ∃ c, c = (Real.sqrt (2 * π * v))⁻¹ := ⟨_, rfl⟩
  have hc0 : 0 < c := by rw [hc]; positivity
  obtain ⟨μ, hμ⟩ : ∃ μ, μ = m + 1 / 2 := ⟨_, rfl⟩
  let vv : NNReal := ⟨v, hv0.le⟩
  have hvv : vv ≠ 0 := by
    intro h
    have : (vv : ℝ) = 0 := by rw [h]; rfl
    exact hv0.ne' this
  -- the Gaussian density `g = exp φ` with `φ(y) = log c - (y - μ)² / (2v)`
  obtain ⟨g, hgdef⟩ : ∃ g : ℝ → ℝ,
      g = fun y => Real.exp (Real.log c - (y - μ) ^ 2 / (2 * v)) :=
    ⟨_, rfl⟩
  have hgy : ∀ y, g y = Real.exp (Real.log c - (y - μ) ^ 2 / (2 * v)) := fun y => by rw [hgdef]
  have hgG : g = ProbabilityTheory.gaussianPDFReal μ vv := by
    rw [hgdef]
    funext y
    rw [ProbabilityTheory.gaussianPDFReal, Real.exp_sub, Real.exp_log hc0, neg_div, Real.exp_neg,
      div_eq_mul_inv, hc]
    rfl
  have hgc : Continuous g := by rw [hgdef]; fun_prop
  have hgi : MeasureTheory.Integrable g := by
    rw [hgG]; exact ProbabilityTheory.integrable_gaussianPDFReal μ vv
  -- weights `w k = ∫_k^{k+1} g` and their lower bounds `exp (B k)`
  obtain ⟨w, hw⟩ : ∃ w : ℕ → ℝ,
      w = fun k : ℕ => ∫ y in (k : ℝ)..((k : ℝ) + 1), g y :=
    ⟨_, rfl⟩
  obtain ⟨B, hB⟩ : ∃ B : ℕ → ℝ,
      B = fun k : ℕ => Real.log c - (((k : ℝ) - m) ^ 2 + 1 / 12) / (2 * v) :=
    ⟨_, rfl⟩
  have hwB : ∀ k, Real.exp (B k) ≤ w k := by
    intro k
    -- tangent line of `exp` at `B k`: `exp x ≥ exp (B k) * (1 + x - B k)`
    have hlow : ∀ y,
        Real.exp (B k) * ((1 - B k + Real.log c) - (y - μ) ^ 2 / (2 * v)) ≤ g y := by
      intro y
      have ht := Real.add_one_le_exp (Real.log c - (y - μ) ^ 2 / (2 * v) - B k)
      have he : g y = Real.exp (B k) * Real.exp (Real.log c - (y - μ) ^ 2 / (2 * v) - B k) := by
        rw [hgy, ← Real.exp_add]
        congr 1
        ring
      rw [he]
      calc Real.exp (B k) * ((1 - B k + Real.log c) - (y - μ) ^ 2 / (2 * v))
          = Real.exp (B k) * ((Real.log c - (y - μ) ^ 2 / (2 * v) - B k) + 1) := by ring
        _ ≤ _ := mul_le_mul_of_nonneg_left ht (Real.exp_pos _).le
    have hint : ∫ y in (k : ℝ)..((k : ℝ) + 1),
        Real.exp (B k) * ((1 - B k + Real.log c) - (y - μ) ^ 2 / (2 * v)) ≤
          ∫ y in (k : ℝ)..((k : ℝ) + 1), g y :=
      intervalIntegral.integral_mono_on (by linarith)
        (by apply Continuous.intervalIntegrable; fun_prop) (hgc.intervalIntegrable _ _)
        (fun y _ => hlow y)
    rw [intervalIntegral.integral_const_mul, integral_tangent_quad _ _ _ _ hv0.ne'] at hint
    have h1 : (1 - B k + Real.log c) - (((k : ℝ) + 1 / 2 - μ) ^ 2 + 1 / 12) / (2 * v) = 1 := by
      simp only [hB, hμ]
      ring
    rw [h1, mul_one] at hint
    simpa only [hw] using hint
  have hw0 : ∀ k, 0 < w k := fun k => (Real.exp_pos _).trans_le (hwB k)
  have hwsum : ∑ k ∈ range (M + 1), w k ≤ 1 := by
    have h := intervalIntegral.sum_integral_adjacent_intervals (μ := MeasureTheory.volume)
      (f := g) (a := fun k : ℕ => (k : ℝ)) (n := M + 1) (fun k _ => hgc.intervalIntegrable _ _)
    simp only [Nat.cast_add, Nat.cast_one, Nat.cast_zero] at h
    calc ∑ k ∈ range (M + 1), w k = ∫ x in (0 : ℝ)..((M : ℝ) + 1), g x := by
          simp only [hw]; exact h
      _ = ∫ x in Set.Ioc (0 : ℝ) ((M : ℝ) + 1), g x :=
          intervalIntegral.integral_of_le (by positivity)
      _ ≤ ∫ x, g x := MeasureTheory.setIntegral_le_integral hgi
          (Filter.Eventually.of_forall fun y => by rw [Pi.zero_apply, hgy]; positivity)
      _ = 1 := by rw [hgG]; exact ProbabilityTheory.integral_gaussianPDFReal_eq_one μ hvv
  -- discrete Gibbs inequality against the weights `w`
  have h1 : shannon M p ≤ ∑ k ∈ range (M + 1), (w k - p k - p k * Real.log (w k)) :=
    Finset.sum_le_sum fun k _ => negMulLog_le_gibbs (hp.nonneg k) (hw0 k)
  have h1' : ∑ k ∈ range (M + 1), (w k - p k - p k * Real.log (w k)) =
      ∑ k ∈ range (M + 1), w k - ∑ k ∈ range (M + 1), p k -
        ∑ k ∈ range (M + 1), p k * Real.log (w k) := by
    rw [Finset.sum_sub_distrib, Finset.sum_sub_distrib]
  have h2 : ∑ k ∈ range (M + 1), p k * B k ≤ ∑ k ∈ range (M + 1), p k * Real.log (w k) :=
    Finset.sum_le_sum fun k _ => mul_le_mul_of_nonneg_left
      ((Real.le_log_iff_exp_le (hw0 k)).mpr (hwB k)) (hp.nonneg k)
  have h3 : ∑ k ∈ range (M + 1), p k * B k = Real.log c - (pvar M p + 1 / 12) / (2 * v) := by
    have : ∑ k ∈ range (M + 1), p k * B k = Real.log c * ∑ k ∈ range (M + 1), p k -
        (∑ k ∈ range (M + 1), p k * ((k : ℝ) - m) ^ 2 +
          1 / 12 * ∑ k ∈ range (M + 1), p k) / (2 * v) := by
      rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib, Finset.sum_div,
        ← Finset.sum_sub_distrib]
      apply Finset.sum_congr rfl
      intro k _
      simp only [hB]
      ring
    rw [this, hp.sum_one, hm, pvar]
    ring
  have h4 : Real.log c = -(1 / 2) * Real.log (2 * π * v) := by
    rw [hc, Real.log_inv, Real.log_sqrt (by positivity)]
    ring
  have h5 : Real.log (2 * π * exp 1 * v) = Real.log (2 * π * v) + 1 := by
    rw [show 2 * π * exp 1 * v = (2 * π * v) * exp 1 by ring,
      Real.log_mul (by positivity) (by positivity), Real.log_exp]
  have h6 : (pvar M p + 1 / 12) / (2 * v) = 1 / 2 := by
    rw [← hv]
    field_simp
  rw [h5]
  rw [h6] at h3
  have hsum := hp.sum_one
  linarith

/-- The Massey bound is monotone in the variance. -/
theorem shannon_le_of_pvar_le (hp : IsProb M p) {V : ℝ} (hV : pvar M p ≤ V) :
    shannon M p ≤ (1 / 2) * Real.log (2 * π * exp 1 * (V + 1 / 12)) := by
  refine (massey hp).trans ?_
  have h0 := pvar_nonneg hp
  gcongr

/-- Shifted second moment: `∑ p (k - c)² = pvar + (mean - c)²` for a normalized `p`. -/
private lemma sum_sq_eq_pvar_add (hp1 : ∑ k ∈ range (M + 1), p k = 1) (c : ℝ) :
    ∑ k ∈ range (M + 1), p k * ((k : ℝ) - c) ^ 2 = pvar M p + (pmean M p - c) ^ 2 := by
  have expand : ∀ d : ℝ, ∑ k ∈ range (M + 1), p k * ((k : ℝ) - d) ^ 2 =
      ∑ k ∈ range (M + 1), p k * (k : ℝ) ^ 2 - 2 * d * pmean M p + d ^ 2 := by
    intro d
    have : ∑ k ∈ range (M + 1), p k * ((k : ℝ) - d) ^ 2 =
        ∑ k ∈ range (M + 1), (p k * (k : ℝ) ^ 2 - 2 * d * (p k * k) + d ^ 2 * p k) := by
      apply Finset.sum_congr rfl
      intro k _
      ring
    rw [this, Finset.sum_add_distrib, Finset.sum_sub_distrib, ← Finset.mul_sum,
      ← Finset.mul_sum, hp1, pmean]
    ring
  rw [pvar, expand c, expand (pmean M p)]
  ring

set_option linter.unusedVariables false in
/-- Variance of a mixture: `Var(∑ w_i p_i) = ∑ w_i Var(p_i) + ∑ w_i (mean_i - mean)²`. -/
theorem pvar_mixture {ι : Type*} [Fintype ι] (w : ι → ℝ) (hw0 : ∀ i, 0 ≤ w i)
    (hw1 : ∑ i, w i = 1) (P : ι → ℕ → ℝ) (hP : ∀ i, IsProb M (P i)) :
    pvar M (fun k => ∑ i, w i * P i k) =
      ∑ i, w i * pvar M (P i) +
        ∑ i, w i * (pmean M (P i) - pmean M (fun k => ∑ i, w i * P i k)) ^ 2 := by
  -- the identity holds for arbitrary weights: `hw0` and `hw1` are not needed
  set m := pmean M (fun k => ∑ i, w i * P i k)
  calc pvar M (fun k => ∑ i, w i * P i k)
      = ∑ k ∈ range (M + 1), (∑ i, w i * P i k) * ((k : ℝ) - m) ^ 2 := rfl
    _ = ∑ k ∈ range (M + 1), ∑ i, w i * (P i k * ((k : ℝ) - m) ^ 2) := by
        apply Finset.sum_congr rfl
        intro k _
        rw [Finset.sum_mul]
        apply Finset.sum_congr rfl
        intro i _
        ring
    _ = ∑ i, w i * ∑ k ∈ range (M + 1), P i k * ((k : ℝ) - m) ^ 2 := by
        rw [Finset.sum_comm]
        apply Finset.sum_congr rfl
        intro i _
        rw [Finset.mul_sum]
    _ = ∑ i, w i * (pvar M (P i) + (pmean M (P i) - m) ^ 2) := by
        apply Finset.sum_congr rfl
        intro i _
        rw [sum_sq_eq_pvar_add (hP i).sum_one]
    _ = _ := by
        rw [← Finset.sum_add_distrib]
        apply Finset.sum_congr rfl
        intro i _
        ring

/-- The variance is at most the mean square deviation from any point. -/
theorem pvar_le_sum_sq (hp : IsProb M p) (c : ℝ) :
    pvar M p ≤ ∑ k ∈ range (M + 1), p k * ((k : ℝ) - c) ^ 2 := by
  rw [sum_sq_eq_pvar_add hp.sum_one c]
  exact le_add_of_nonneg_right (sq_nonneg _)

/-- `-t log t ≤ 2 (√t - t)` for `t ≥ 0`. -/
private lemma negMulLog_le_two_sqrt_sub {t : ℝ} (ht : 0 ≤ t) :
    negMulLog t ≤ 2 * (Real.sqrt t - t) := by
  -- `-t log t = -2t log √t ≤ 2t (1/√t - 1)` by `log x ≥ 1 - 1/x`
  rcases ht.eq_or_lt with rfl | ht
  · simp
  · have hs : 0 < Real.sqrt t := Real.sqrt_pos.mpr ht
    have h1 := Real.one_sub_inv_le_log_of_pos hs
    rw [Real.log_sqrt ht.le] at h1
    have h2 : t * (Real.sqrt t)⁻¹ = Real.sqrt t := by
      rw [← div_eq_mul_inv, Real.div_sqrt]
    have h3 : t * (1 - (Real.sqrt t)⁻¹) ≤ t * (Real.log t / 2) :=
      mul_le_mul_of_nonneg_left h1 ht.le
    rw [negMulLog]
    nlinarith

/-- Binary entropy bound `h(ν) ≤ 2 √(ν(1-ν))` (used for the area law). Proof idea: Rényi
monotonicity `h(ν) ≤ 2 ln(√ν + √(1-ν)) = ln(1 + 2√(ν(1-ν)))` (Jensen for `log`), then
`ln(1+x) ≤ x`. (The paper uses the sharper constant `2 ln 2`.) -/
theorem binEntropy_le_two_sqrt {ν : ℝ} (h0 : 0 ≤ ν) (h1 : ν ≤ 1) :
    binEntropy ν ≤ 2 * Real.sqrt (ν * (1 - ν)) := by
  -- termwise `-t log t ≤ 2(√t - t)` gives `h(ν) ≤ 2(√ν + √(1-ν) - 1)`, and
  -- `√ν + √(1-ν) - 1 ≤ √ν √(1-ν)` is `(1 - √ν)(1 - √(1-ν)) ≥ 0`
  rw [binEntropy_eq_negMulLog_add_negMulLog_one_sub, Real.sqrt_mul h0]
  have ha := negMulLog_le_two_sqrt_sub h0
  have hb := negMulLog_le_two_sqrt_sub (sub_nonneg.mpr h1)
  have hx1 : Real.sqrt ν ≤ 1 := Real.sqrt_le_one.mpr h1
  have hy1 : Real.sqrt (1 - ν) ≤ 1 := Real.sqrt_le_one.mpr (by linarith)
  nlinarith [mul_nonneg (sub_nonneg.mpr hx1) (sub_nonneg.mpr hy1)]

end DFL
