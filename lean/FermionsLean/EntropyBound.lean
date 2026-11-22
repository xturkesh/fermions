import Mathlib.Analysis.SpecialFunctions.BinaryEntropy
import Mathlib.Analysis.Convex.Deriv
import Mathlib.Analysis.Convex.Jensen

/-!
# Entanglement entropy versus particle-number fluctuations for Gaussian states

For a Slater determinant, the entanglement entropy of a region `A` and the variance of the
particle number in `A` are both spectral functions of the reduced correlation matrix `C_A`,
whose eigenvalues `ν_k` lie in `[0, 1]`:

  `S_A = ∑ₖ binEntropy ν_k`,   `Var N_A = ∑ₖ ν_k (1 - ν_k) = Tr C_A (1 - C_A)`.

* `four_log_two_mul_le_binEntropy`: `4 log 2 · p (1 - p) ≤ binEntropy p` on `[0, 1]`
  (sharp at `p ∈ {0, 1/2, 1}`), hence
* `four_log_two_mul_var_le_entropy`: `4 log 2 · Var N_A ≤ S_A`. A log-law (or any growth) of
  the particle-number variance forces at least the same growth of the entanglement entropy.
* `entropy_le_log_two_mul_card`: `S_A ≤ log 2 · #{k : ν_k ∉ {0, 1}}`, the number of entangled
  modes (equal to `rank C_AB`, see `FermionsLean.CutRank`).
-/

open Real Set

namespace MonitoredFermions

lemma half_lt_log_two : (1 : ℝ) / 2 < log 2 := by
  have h := Real.log_lt_sub_one_of_pos (by norm_num : (0 : ℝ) < 1 / 2) (by norm_num)
  have h2 : log (1 / 2 : ℝ) = -log 2 := by rw [one_div, Real.log_inv]
  linarith

/-- `g p = binEntropy p - 4 log 2 · p (1 - p)`. -/
noncomputable def gFun (x : ℝ) : ℝ := binEntropy x - 4 * log 2 * (x * (1 - x))

/-- The derivative of `gFun` on `(0, 1)`. -/
noncomputable def gDeriv (x : ℝ) : ℝ := log (1 - x) - log x - 4 * log 2 * (1 - 2 * x)

lemma gFun_continuous : Continuous gFun := by
  unfold gFun
  fun_prop

lemma hasDerivAt_gFun {x : ℝ} (h0 : x ≠ 0) (h1 : x ≠ 1) : HasDerivAt gFun (gDeriv x) x := by
  have hb := hasDerivAt_binEntropy h0 h1
  have hq : HasDerivAt (fun y : ℝ => y * (1 - y)) (1 - 2 * x) x := by
    have := (hasDerivAt_id' x).mul ((hasDerivAt_id' x).const_sub 1)
    convert this using 1
    ring
  exact hb.sub (hq.const_mul (4 * log 2))

lemma hasDerivAt_gDeriv {x : ℝ} (h0 : x ≠ 0) (h1 : x ≠ 1) :
    HasDerivAt gDeriv (8 * log 2 - 1 / (x * (1 - x))) x := by
  have h1' : 1 - x ≠ 0 := sub_ne_zero.mpr (Ne.symm h1)
  have ha : HasDerivAt (fun y : ℝ => log (1 - y)) (-1 / (1 - x)) x :=
    ((hasDerivAt_id' x).const_sub 1).log h1'
  have hb : HasDerivAt (fun y : ℝ => log y) x⁻¹ x := Real.hasDerivAt_log h0
  have hc : HasDerivAt (fun y : ℝ => 4 * log 2 * (1 - 2 * y)) (4 * log 2 * (-(2 * 1))) x :=
    (((hasDerivAt_id' x).const_mul 2).const_sub 1).const_mul (4 * log 2)
  refine ((ha.sub hb).sub hc).congr_deriv ?_
  field_simp
  ring

/-- The point `s ∈ (0, 1/2)` where `s (1 - s) = 1 / (8 log 2)`: `g` is concave on `[0, s]`
and `g'` is increasing on `[s, 1/2]`. -/
noncomputable def sStar : ℝ := (1 - √(1 - 1 / (2 * log 2))) / 2

lemma sStar_facts : 0 < sStar ∧ sStar < 2⁻¹ ∧ sStar * (1 - sStar) = 1 / (8 * log 2) := by
  have hL : (1 : ℝ) / 2 < log 2 := half_lt_log_two
  have hLpos : 0 < log 2 := by linarith
  have hpos : 0 < 1 / (2 * log 2) := by positivity
  have hin : 0 < 1 - 1 / (2 * log 2) := by
    have : 1 / (2 * log 2) < 1 := by
      rw [div_lt_one (by linarith)]
      linarith
    linarith
  obtain ⟨r, hr⟩ : ∃ r, r = √(1 - 1 / (2 * log 2)) := ⟨_, rfl⟩
  have hr0 : 0 < r := hr ▸ Real.sqrt_pos.mpr hin
  have hr2 : r ^ 2 = 1 - 1 / (2 * log 2) := hr ▸ Real.sq_sqrt hin.le
  have hr1 : r < 1 := by nlinarith
  have hs : sStar = (1 - r) / 2 := by rw [sStar, hr]
  rw [hs]
  refine ⟨by linarith, by linarith, ?_⟩
  have e : (1 - r) / 2 * (1 - (1 - r) / 2) = (1 - r ^ 2) / 4 := by ring
  rw [e, hr2]
  field_simp
  ring

lemma gDeriv_half : gDeriv 2⁻¹ = 0 := by
  norm_num [gDeriv]

lemma gFun_half : gFun 2⁻¹ = 0 := by
  rw [gFun, binEntropy_two_inv]
  ring

lemma gFun_zero : gFun 0 = 0 := by
  simp [gFun]

/-- On `[s, 1/2]`, `g' ≤ 0`, because `g'` is increasing there and vanishes at `1/2`. -/
lemma gDeriv_nonpos_of_mem {x : ℝ} (hx : x ∈ Icc sStar 2⁻¹) : gDeriv x ≤ 0 := by
  obtain ⟨hs0, hs1, hss⟩ := sStar_facts
  have hL : (1 : ℝ) / 2 < log 2 := half_lt_log_two
  have hmono : MonotoneOn gDeriv (Icc sStar 2⁻¹) := by
    apply monotoneOn_of_deriv_nonneg (convex_Icc _ _)
    · intro y hy
      have hy0 : y ≠ 0 := by linarith [hy.1]
      have hy1 : y ≠ 1 := by linarith [hy.2]
      exact (hasDerivAt_gDeriv hy0 hy1).continuousAt.continuousWithinAt
    · intro y hy
      rw [interior_Icc] at hy
      have hy0 : y ≠ 0 := by linarith [hy.1]
      have hy1 : y ≠ 1 := by linarith [hy.2]
      exact (hasDerivAt_gDeriv hy0 hy1).differentiableAt.differentiableWithinAt
    · intro y hy
      rw [interior_Icc] at hy
      have hy0 : y ≠ 0 := by linarith [hy.1]
      have hy1 : y ≠ 1 := by linarith [hy.2]
      rw [(hasDerivAt_gDeriv hy0 hy1).deriv]
      have hyy : sStar * (1 - sStar) ≤ y * (1 - y) := by nlinarith [hy.1, hy.2]
      have hypos : 0 < y * (1 - y) := by nlinarith [hy.1, hy.2]
      rw [hss] at hyy
      have h8 : 1 ≤ y * (1 - y) * (8 * log 2) :=
        (div_le_iff₀ (by positivity : (0 : ℝ) < 8 * log 2)).mp hyy
      have : 1 / (y * (1 - y)) ≤ 8 * log 2 := by
        rw [div_le_iff₀ hypos]
        linarith
      linarith
  calc gDeriv x ≤ gDeriv 2⁻¹ := hmono hx ⟨hs1.le, le_refl _⟩ hx.2
    _ = 0 := gDeriv_half

/-- On `[s, 1/2]`, `g ≥ 0` (it decreases to `g (1/2) = 0`). -/
lemma gFun_nonneg_of_mem_right {x : ℝ} (hx : x ∈ Icc sStar 2⁻¹) : 0 ≤ gFun x := by
  obtain ⟨hs0, hs1, _⟩ := sStar_facts
  have hanti : AntitoneOn gFun (Icc sStar 2⁻¹) := by
    apply antitoneOn_of_deriv_nonpos (convex_Icc _ _) gFun_continuous.continuousOn
    · intro y hy
      rw [interior_Icc] at hy
      have hy0 : y ≠ 0 := by linarith [hy.1]
      have hy1 : y ≠ 1 := by linarith [hy.2]
      exact (hasDerivAt_gFun hy0 hy1).differentiableAt.differentiableWithinAt
    · intro y hy
      rw [interior_Icc] at hy
      have hy0 : y ≠ 0 := by linarith [hy.1]
      have hy1 : y ≠ 1 := by linarith [hy.2]
      rw [(hasDerivAt_gFun hy0 hy1).deriv]
      exact gDeriv_nonpos_of_mem ⟨hy.1.le, hy.2.le⟩
  calc 0 = gFun 2⁻¹ := gFun_half.symm
    _ ≤ gFun x := hanti hx ⟨hs1.le, le_refl _⟩ hx.2

/-- On `[0, s]`, `g` is concave, so `g ≥ min (g 0) (g s) ≥ 0`. -/
lemma gFun_nonneg_of_mem_left {x : ℝ} (hx : x ∈ Icc 0 sStar) : 0 ≤ gFun x := by
  obtain ⟨hs0, hs1, hss⟩ := sStar_facts
  have hL : (1 : ℝ) / 2 < log 2 := half_lt_log_two
  have hantiD : AntitoneOn gDeriv (Ioo 0 sStar) := by
    apply antitoneOn_of_deriv_nonpos (convex_Ioo _ _)
    · intro y hy
      have hy0 : y ≠ 0 := by linarith [hy.1]
      have hy1 : y ≠ 1 := by linarith [hy.2]
      exact (hasDerivAt_gDeriv hy0 hy1).continuousAt.continuousWithinAt
    · intro y hy
      rw [interior_Ioo] at hy
      have hy0 : y ≠ 0 := by linarith [hy.1]
      have hy1 : y ≠ 1 := by linarith [hy.2]
      exact (hasDerivAt_gDeriv hy0 hy1).differentiableAt.differentiableWithinAt
    · intro y hy
      rw [interior_Ioo] at hy
      have hy0 : y ≠ 0 := by linarith [hy.1]
      have hy1 : y ≠ 1 := by linarith [hy.2]
      rw [(hasDerivAt_gDeriv hy0 hy1).deriv]
      have hyy : y * (1 - y) ≤ sStar * (1 - sStar) := by nlinarith [hy.1, hy.2]
      have hypos : 0 < y * (1 - y) := by nlinarith [hy.1, hy.2]
      rw [hss] at hyy
      have h8 : y * (1 - y) * (8 * log 2) ≤ 1 :=
        (le_div_iff₀ (by positivity : (0 : ℝ) < 8 * log 2)).mp hyy
      have : 8 * log 2 ≤ 1 / (y * (1 - y)) := by
        rw [le_div_iff₀ hypos]
        linarith
      linarith
  have hconc : ConcaveOn ℝ (Icc 0 sStar) gFun := by
    apply AntitoneOn.concaveOn_of_deriv (convex_Icc _ _) gFun_continuous.continuousOn
    · intro y hy
      rw [interior_Icc] at hy
      have hy0 : y ≠ 0 := by linarith [hy.1]
      have hy1 : y ≠ 1 := by linarith [hy.2]
      exact (hasDerivAt_gFun hy0 hy1).differentiableAt.differentiableWithinAt
    · rw [interior_Icc]
      intro a ha b hb hab
      have ha0 : a ≠ 0 := by linarith [ha.1]
      have ha1 : a ≠ 1 := by linarith [ha.2]
      have hb0 : b ≠ 0 := by linarith [hb.1]
      have hb1 : b ≠ 1 := by linarith [hb.2]
      rw [(hasDerivAt_gFun ha0 ha1).deriv, (hasDerivAt_gFun hb0 hb1).deriv]
      exact hantiD ha hb hab
  have hs : 0 ≤ gFun sStar := gFun_nonneg_of_mem_right ⟨le_refl _, hs1.le⟩
  have hmin := hconc.min_le_of_mem_Icc (left_mem_Icc.mpr hs0.le) (right_mem_Icc.mpr hs0.le) hx
  rw [gFun_zero] at hmin
  exact le_trans (le_min (le_refl 0) hs) hmin

/-- **`4 log 2 · p (1 - p) ≤ binEntropy p`** for `p ∈ [0, 1]`. -/
theorem four_log_two_mul_le_binEntropy {p : ℝ} (hp0 : 0 ≤ p) (hp1 : p ≤ 1) :
    4 * log 2 * (p * (1 - p)) ≤ binEntropy p := by
  have key : ∀ q ∈ Icc (0 : ℝ) 2⁻¹, 4 * log 2 * (q * (1 - q)) ≤ binEntropy q := by
    intro q hq
    have h : 0 ≤ gFun q := by
      rcases le_total q sStar with h | h
      · exact gFun_nonneg_of_mem_left ⟨hq.1, h⟩
      · exact gFun_nonneg_of_mem_right ⟨h, hq.2⟩
    simp only [gFun] at h
    linarith
  rcases le_total p 2⁻¹ with h | h
  · exact key p ⟨hp0, h⟩
  · have h' := key (1 - p) ⟨by linarith, by linarith⟩
    rw [binEntropy_one_sub] at h'
    have e : (1 - p) * (1 - (1 - p)) = p * (1 - p) := by ring
    rwa [e] at h'

/-- **Entanglement ≥ 4 log 2 × particle-number variance** for Gaussian states, phrased on the
spectrum `ν` of the reduced correlation matrix. -/
theorem four_log_two_mul_var_le_entropy {ι : Type*} (s : Finset ι) (ν : ι → ℝ)
    (hν : ∀ k ∈ s, ν k ∈ Icc (0 : ℝ) 1) :
    4 * log 2 * ∑ k ∈ s, ν k * (1 - ν k) ≤ ∑ k ∈ s, binEntropy (ν k) := by
  rw [Finset.mul_sum]
  exact Finset.sum_le_sum fun k hk => four_log_two_mul_le_binEntropy (hν k hk).1 (hν k hk).2

/-- **`S_A ≤ log 2 · (number of entangled modes)`**: eigenvalues `0` and `1` do not contribute,
every other eigenvalue contributes at most `log 2`. -/
theorem entropy_le_log_two_mul_card {ι : Type*} (s : Finset ι) (ν : ι → ℝ) :
    ∑ k ∈ s, binEntropy (ν k) ≤ log 2 * (s.filter fun k => ν k ≠ 0 ∧ ν k ≠ 1).card := by
  classical
  rw [← Finset.sum_filter_add_sum_filter_not s (fun k => ν k ≠ 0 ∧ ν k ≠ 1)]
  have h2 : ∑ k ∈ s.filter (fun k => ¬(ν k ≠ 0 ∧ ν k ≠ 1)), binEntropy (ν k) = 0 := by
    apply Finset.sum_eq_zero
    intro k hk
    simp only [Finset.mem_filter, not_and_or, not_not] at hk
    rcases hk.2 with h | h <;> simp [h]
  rw [h2, add_zero]
  calc ∑ k ∈ s.filter (fun k => ν k ≠ 0 ∧ ν k ≠ 1), binEntropy (ν k)
      ≤ ∑ _k ∈ s.filter (fun k => ν k ≠ 0 ∧ ν k ≠ 1), log 2 :=
        Finset.sum_le_sum fun _ _ => binEntropy_le_log_two
    _ = log 2 * (s.filter fun k => ν k ≠ 0 ∧ ν k ≠ 1).card := by
        rw [Finset.sum_const, nsmul_eq_mul, mul_comm]

end MonitoredFermions
