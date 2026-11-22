import Mathlib

/-!
# The binary entropy (paper: Lemma `lem:binent`)

`h(ν) ≤ 2 ln 2 √(ν(1-ν))` on `[0, 1]`, with equality exactly at `ν ∈ {0, ½, 1}`. With the exact
formula `S(ρ_A) = ∑_i h(ν_i)` for Slater determinants and the trace-norm bound, this gives the area
law with the constant `2 ln 2` of the paper.

## Proof strategy

Instead of the power-series argument of the paper we use calculus: the auxiliary function
`R(ν) = h(ν)² / (ν(1-ν))` has derivative `h(ν) ((1-ν) ln(1-ν) - ν ln ν) / (ν(1-ν))²`, which is
positive on `(0, ½)` by the elementary inequality `ν ln ν < (1-ν) ln(1-ν)` there (a consequence of
`ln x ≤ x - 1`). Hence `R` is strictly increasing on `(0, ½]`, so `R(ν) < R(½) = (2 ln 2)²` for
`ν ∈ (0, ½)`; the case `ν ∈ (½, 1)` follows by the symmetry `ν ↦ 1 - ν`.
-/

noncomputable section

namespace DFL

namespace BinEntropyAux

open Real

/-- Elementary inequality: `ν ln ν < (1 - ν) ln (1 - ν)` for `0 < ν < 1/2`. -/
lemma mul_log_lt_one_sub_mul_log {ν : ℝ} (h0 : 0 < ν) (h1 : ν < 1 / 2) :
    ν * log ν < (1 - ν) * log (1 - ν) := by
  have hb : 0 < 1 - ν := by linarith
  have hu : 0 < ν / (1 - ν) := div_pos h0 hb
  have hu1 : ν / (1 - ν) ≠ 1 := by
    rw [Ne, div_eq_one_iff_eq hb.ne']
    intro h
    linarith
  have e1 := log_lt_sub_one_of_pos hu hu1
  have e2 := log_le_sub_one_of_pos (inv_pos.mpr hb)
  rw [log_div h0.ne' hb.ne'] at e1
  rw [log_inv] at e2
  have k1 : ν * (log ν - log (1 - ν)) < ν * (ν / (1 - ν) - 1) := mul_lt_mul_of_pos_left e1 h0
  have k2 : (1 - 2 * ν) * (-log (1 - ν)) ≤ (1 - 2 * ν) * ((1 - ν)⁻¹ - 1) :=
    mul_le_mul_of_nonneg_left e2 (by linarith)
  have k3 : ν * (ν / (1 - ν) - 1) + (1 - 2 * ν) * ((1 - ν)⁻¹ - 1) = 0 := by
    field_simp
    ring
  nlinarith

/-- The auxiliary function `h(ν)² / (ν (1 - ν))`. -/
def ratio (ν : ℝ) : ℝ := binEntropy ν ^ 2 / (ν * (1 - ν))

/-- Derivative of `ratio` on `(0, 1)`. -/
lemma hasDerivAt_ratio {ν : ℝ} (h0 : 0 < ν) (h1 : ν < 1) :
    HasDerivAt ratio
      (binEntropy ν * ((1 - ν) * log (1 - ν) - ν * log ν) / (ν * (1 - ν)) ^ 2) ν := by
  have hb := hasDerivAt_binEntropy h0.ne' h1.ne
  have hd : HasDerivAt (fun x : ℝ => x * (1 - x)) (1 * (1 - ν) + ν * (0 - 1)) ν :=
    (hasDerivAt_id ν).mul ((hasDerivAt_const ν (1 : ℝ)).sub (hasDerivAt_id ν))
  have hne : ν * (1 - ν) ≠ 0 := by
    have : 0 < 1 - ν := by linarith
    positivity
  refine ((hb.fun_pow 2).fun_div hd hne).congr_deriv ?_
  simp only [binEntropy, log_inv]
  norm_num
  ring

/-- `ratio` is strictly increasing on `(0, 1/2]`. -/
lemma ratio_strictMonoOn : StrictMonoOn ratio (Set.Ioc 0 (1 / 2)) := by
  apply strictMonoOn_of_deriv_pos (convex_Ioc 0 (1 / 2))
  · intro x hx
    have h0 : 0 < x := hx.1
    have h1 : 0 < 1 - x := by linarith [hx.2]
    apply ContinuousAt.continuousWithinAt
    exact ((binEntropy_continuous.continuousAt).pow 2).div (by fun_prop) (by positivity)
  · intro x hx
    rw [interior_Ioc] at hx
    obtain ⟨h0, h1⟩ := hx
    rw [(hasDerivAt_ratio h0 (by linarith)).deriv]
    have hb : 0 < 1 - x := by linarith
    have hpos := binEntropy_pos h0 (by linarith : x < 1)
    have hkey := mul_log_lt_one_sub_mul_log h0 h1
    apply div_pos
    · exact mul_pos hpos (by linarith)
    · positivity

/-- Value of `ratio` at `1/2`. -/
lemma ratio_half : ratio (1 / 2) = (2 * log 2) ^ 2 := by
  rw [ratio, one_div, binEntropy_two_inv]
  norm_num
  ring

/-- Strict inequality on `(0, 1/2)`. -/
lemma lt_of_pos_of_lt_half {ν : ℝ} (h0 : 0 < ν) (h1 : ν < 1 / 2) :
    binEntropy ν < 2 * log 2 * √(ν * (1 - ν)) := by
  have hm := ratio_strictMonoOn ⟨h0, h1.le⟩ ⟨by norm_num, le_rfl⟩ h1
  have hb : 0 < 1 - ν := by linarith
  have hd : 0 < ν * (1 - ν) := by positivity
  rw [ratio_half, ratio, div_lt_iff₀ hd] at hm
  have hs : √(ν * (1 - ν)) ^ 2 = ν * (1 - ν) := sq_sqrt hd.le
  have hl : 0 < log 2 := log_pos (by norm_num)
  have hR : 0 ≤ 2 * log 2 * √(ν * (1 - ν)) := by positivity
  have hsq : binEntropy ν ^ 2 < (2 * log 2 * √(ν * (1 - ν))) ^ 2 := by
    rw [mul_pow, hs]
    exact hm
  exact lt_of_pow_lt_pow_left₀ 2 hR hsq

/-- Strict inequality on `(0, 1) \ {1/2}`. -/
lemma lt_of_ne_half {ν : ℝ} (h0 : 0 < ν) (h1 : ν < 1) (hne : ν ≠ 1 / 2) :
    binEntropy ν < 2 * log 2 * √(ν * (1 - ν)) := by
  rcases lt_or_gt_of_ne hne with h | h
  · exact lt_of_pos_of_lt_half h0 h
  · have := lt_of_pos_of_lt_half (ν := 1 - ν) (by linarith) (by linarith)
    rwa [binEntropy_one_sub, show (1 - ν) * (1 - (1 - ν)) = ν * (1 - ν) by ring] at this

/-- Equality at `1/2`. -/
lemma eq_half : binEntropy (1 / 2) = 2 * log 2 * √((1 / 2 : ℝ) * (1 - 1 / 2)) := by
  rw [show ((1 : ℝ) / 2) * (1 - 1 / 2) = (1 / 2) ^ 2 by norm_num, sqrt_sq (by norm_num), one_div,
    binEntropy_two_inv]
  ring

end BinEntropyAux

/-- **Binary entropy bound** (paper: Lemma `lem:binent`). -/
theorem binEntropy_le_two_log_two_sqrt (ν : ℝ) (h0 : 0 ≤ ν) (h1 : ν ≤ 1) :
    Real.binEntropy ν ≤ 2 * Real.log 2 * √(ν * (1 - ν)) := by
  rcases h0.eq_or_lt with rfl | h0
  · simp
  rcases h1.eq_or_lt with rfl | h1
  · simp
  by_cases hne : ν = 1 / 2
  · subst hne
    exact BinEntropyAux.eq_half.le
  · exact (BinEntropyAux.lt_of_ne_half h0 h1 hne).le

/-- **Equality cases** (paper: Lemma `lem:binent`): equality holds exactly for `ν ∈ {0, ½, 1}`. -/
theorem binEntropy_eq_two_log_two_sqrt_iff (ν : ℝ) (h0 : 0 ≤ ν) (h1 : ν ≤ 1) :
    Real.binEntropy ν = 2 * Real.log 2 * √(ν * (1 - ν)) ↔ ν = 0 ∨ ν = 1 / 2 ∨ ν = 1 := by
  constructor
  · intro heq
    by_contra hcon
    simp only [not_or] at hcon
    obtain ⟨ha, hb, hc⟩ := hcon
    have := BinEntropyAux.lt_of_ne_half (lt_of_le_of_ne h0 (Ne.symm ha)) (lt_of_le_of_ne h1 hc) hb
    linarith
  · rintro (rfl | rfl | rfl)
    · simp
    · exact BinEntropyAux.eq_half
    · simp

end DFL
