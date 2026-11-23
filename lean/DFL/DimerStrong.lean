import DFL.DimerInfVol
import DFL.DimerState

/-!
# Strong coupling: the dressed fermion stays tight (paper: Remark `rem:dimerstrong`)

For `λ = 2|h|/J > 1`, `λ ≠ √2`, De Bièvre and Germinet prove dynamical localization for the random dimer
model on `ℤ`: almost surely, for every `p > 0`, `sup_t ∑_y |y|^p |⟨δ_y, e^{-itH}δ_0⟩|² < ∞`. This is a
literature input; here it enters as the hypothesis `hDL` for one `p > 0` (for the operator
`-Δ + v^{(r)}` of `DimerInfVol.lean`, which is the operator of De Bièvre–Germinet up to the diagonal gauge
transformation `diag((-1)^x)` and a translation).

* `dimerTailLGT`: the weight `⟨Ψ^x_N(t), ∑_{|y-x|>R} n_y Ψ^x_N(t)⟩` of the single-fermion dimer state outside
  the window of radius `R` (normalized as `dimerDLGT`); `dimerTailLGT_eq`: it is the dimer-sector
  average of `∑_{|y-x|>R} |⟨δ_y, e^{-ith_q}δ_x⟩|²` (eq. `eq:dimerdensity`).
* `dimer_TL_tail`: its thermodynamic limit is `dimerTailZ R t = E ∑_{|y|>R} |⟨δ_y, e^{-iJtH}δ_0⟩|²`
  (Lemma `lem:dimerTL`, as in the proof of Theorem `thm:dimer`(b)).
* `dimerTailZ_tendsto_zero`: given dynamical localization, `sup_t dimerTailZ R t → 0` as `R → ∞`
  (dominated convergence; the quantities are at most `1`).
* `dimer_tight_of_dynLoc`: the statement of the remark,
  `lim_{R→∞} sup_t lim_{N→∞} ⟨Ψ^{x_N}_N(t), ∑_{|y-x_N|>R} n_y Ψ^{x_N}_N(t)⟩ = 0`.
-/

open MeasureTheory Filter Topology Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

/-- The weight of the single-fermion dimer state `Ψ^x_N(t)` outside the window of radius `R` around
`x`, `⟨Ψ^x_N(t), ∑_{|y-x|>R} n_y Ψ^x_N(t)⟩` (normalized as `dimerDLGT`). -/
def dimerTailLGT (J h : ℝ) (N R : ℕ) (t : ℝ) (x : Fin N) : ℝ :=
  (1 / 2 ^ (N / 2) : ℝ) *
    (expect (cplx (diagonal (fun n : Cfg N => ∑ y : Fin N,
        (if (R : ℤ) < |(y.val : ℤ) - x.val| then (1 : ℝ) else 0) * occ (n y)) ⊗ₖ
      (1 : Matrix (Cfg N) (Cfg N) ℝ))) (evo (H J h) t *ᵥ dimerState (fockCfg {x}))).re

/-- The infinite-volume counterpart `E ∑_{|y|>R} |⟨δ_y, e^{-iJtH}δ_0⟩|²` for the random dimer potential
`v^{(r)}`. -/
def dimerTailZ (J h : ℝ) (r : ℤ) (R : ℕ) (t : ℝ) : ℝ :=
  ∫ η, ∑' y : ℤ, (if (R : ℤ) < |y| then ‖evoZw (dimerPotZ J h r η) (J * t) y 0‖ ^ 2 else 0) ∂muZ

/-- The tail weight of the dimer state is the dimer-sector average of the one-particle tail weight
(eq. `eq:dimerdensity`). -/
theorem dimerTailLGT_eq (J h : ℝ) {N : ℕ} (hN : Even N) (R : ℕ) (t : ℝ) (x : Fin N) :
    dimerTailLGT J h N R t x = dimerAvg (fun q : Cfg N =>
      ∑ y : Fin N, (if (R : ℤ) < |(y.val : ℤ) - x.val| then ‖evo (hq J h q) t y x‖ ^ 2 else 0)) := by
  unfold dimerTailLGT
  have e0 : (fun q : Cfg N => ∑ y : Fin N,
      (if (R : ℤ) < |(y.val : ℤ) - x.val| then ‖evo (hq J h q) t y x‖ ^ 2 else 0)) =
      fun q => ∑ y : Fin N, (if (R : ℤ) < |(y.val : ℤ) - x.val| then (1 : ℝ) else 0) *
        ‖evo (hq J h q) t y x‖ ^ 2 := by
    funext q
    refine Finset.sum_congr rfl fun y _ => ?_
    split_ifs <;> simp
  rw [e0, DimerLGTAux.dimerAvg_sum_mul (N := N)
    (fun y : Fin N => if (R : ℤ) < |(y.val : ℤ) - x.val| then (1 : ℝ) else 0)
    (fun y q => ‖evo (hq J h q) t y x‖ ^ 2)]
  simp only [← dimer_density_single J h hN x, DimerStateAux.expect_matter_re]
  simp only [Finset.sum_mul, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun y _ => Finset.sum_congr rfl fun r _ => ?_
  ring

namespace DimerStrongAux

/-- The tail weight of `h_{Φ(η)}` at `x` is that of `-Δ + v^{(r)}` on `[-x, N-1-x]` at `0`, time `Jt`. -/
lemma tail_dimerPhi (J h : ℝ) (hJ : J ≠ 0) (N : ℕ) (x : Fin N) (r : ℤ) (η : CfgZ) (R : ℕ)
    (t : ℝ) :
    ∑ y : Fin N, (if (R : ℤ) < |(y.val : ℤ) - x.val| then
        ‖evo (hq J h (DimerInfVolAux.dimerPhi N x.val r η)) t y x‖ ^ 2 else 0) =
      ∑ y : Set.Icc (-(x.val : ℤ)) ((N : ℤ) - 1 - x.val), (if (R : ℤ) < |y.1| then
        ‖evo (schrH (-(x.val : ℤ)) ((N : ℤ) - 1 - x.val) (dimerPotZ J h r η)) (J * t) y
          (origin _ _ (by omega) (by have := x.isLt; omega))‖ ^ 2 else 0) := by
  rw [DimerInfVolAux.hq_dimerPhi J h hJ N x.val r η, ← evo_submatrix, DimerLGTAux.evo_smul,
    ← (DimerInfVolAux.shiftIcc N x.val).sum_comp]
  refine Finset.sum_congr rfl fun y _ => ?_
  have ho : DimerInfVolAux.shiftIcc N x.val x =
      origin _ _ (by omega) (by have := x.isLt; omega) := by
    ext
    simp [origin]
  rw [Matrix.submatrix_apply, ho, DimerInfVolAux.shiftIcc_apply]

/-- A partial column weight of a propagator is at most `1`. -/
lemma sum_ite_norm_sq_le_one {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℝ)
    (hA : A.IsSymm) (s : ℝ) (z : ι) (P : ι → Prop) [DecidablePred P] :
    ∑ y, (if P y then ‖evo A s y z‖ ^ 2 else 0) ≤ 1 := by
  refine (Finset.sum_le_sum fun y _ => ?_).trans (DimerTLAux.sum_norm_sq_evo A hA s z).le
  split_ifs
  · exact le_rfl
  · positivity

/-- `y ↦ 4^{-|y|}` is summable on `ℤ`. -/
lemma summable_geom_Z : Summable (fun y : ℤ => (1 / 4 : ℝ) ^ y.natAbs) := by
  have hN : Summable (fun n : ℕ => (1 / 4 : ℝ) ^ n) :=
    summable_geometric_of_lt_one (by norm_num) (by norm_num)
  refine Summable.of_nat_of_neg ?_ ?_
  · refine hN.congr fun n => ?_
    simp
  · refine hN.congr fun n => ?_
    simp

/-- The partial column weights `∑_{y ∈ P} |⟨δ_y, e^{-isH_Λ} δ_0⟩|²` converge as `Λ ↑ ℤ`. -/
lemma tendsto_ite_Z (w : ℤ → ℝ) (W : ℝ) (hw : ∀ x, |w x| ≤ W) (as bs : ℕ → ℤ)
    (ha : ∀ n, as n ≤ 0) (hb : ∀ n, 0 ≤ bs n) (has : Tendsto as atTop atBot)
    (hbs : Tendsto bs atTop atTop) (P : ℤ → Prop) [DecidablePred P] (s : ℝ) :
    Tendsto (fun n => ∑ y : Set.Icc (as n) (bs n), (if P y.1 then
        ‖evo (schrH (as n) (bs n) w) s y (origin (as n) (bs n) (ha n) (hb n))‖ ^ 2 else 0))
      atTop (𝓝 (∑' y : ℤ, if P y then ‖evoZw w s y 0‖ ^ 2 else 0)) := by
  have hW : 0 ≤ W := (abs_nonneg _).trans (hw 0)
  set M := ⌈4 * Real.exp 2 * |s|⌉₊ with hM
  set G : ℕ → ℤ → ℂ := fun n y => if h : as n ≤ y ∧ y ≤ bs n then
      evo (schrH (as n) (bs n) w) s ⟨y, h⟩ (origin (as n) (bs n) (ha n) (hb n)) else 0 with hG
  have hF : ∀ n, ∑ y : Set.Icc (as n) (bs n), (if P y.1 then
      ‖evo (schrH (as n) (bs n) w) s y (origin (as n) (bs n) (ha n) (hb n))‖ ^ 2 else 0) =
      ∑' y : ℤ, (if P y then ‖G n y‖ ^ 2 else 0) := by
    intro n
    rw [DimerInfVolAux.sum_Icc_eq_tsum]
    congr 1
    funext y
    by_cases h : as n ≤ y ∧ y ≤ bs n
    · simp only [hG, dite_eq_left h]
    · simp only [hG, dite_eq_right h, norm_zero]
      simp
  simp_rw [hF]
  have hm : Tendsto (fun n => (min (bs n) (-as n)).toNat) atTop atTop := by
    rw [tendsto_atTop]
    intro K
    filter_upwards [hbs.eventually_ge_atTop (K : ℤ), has.eventually_le_atBot (-(K : ℤ))]
      with n h1 h2
    omega
  have hτ := ((DimerInfVolAux.tendsto_tail ((2 + W) * |s|) (by positivity)).comp hm).const_mul 2
  rw [mul_zero] at hτ
  refine tendsto_tsum_of_dominated_convergence
    (bound := fun y : ℤ => 4 ^ M * (1 / 4 : ℝ) ^ y.natAbs)
    (summable_geom_Z.mul_left _) (fun y => ?_) (Eventually.of_forall fun n y => ?_)
  · have hGt : Tendsto (fun n => G n y) atTop (𝓝 (evoZw w s y 0)) := by
      rw [tendsto_iff_norm_sub_tendsto_zero]
      refine squeeze_zero' (Eventually.of_forall fun n => norm_nonneg _) ?_ hτ
      filter_upwards [has.eventually_le_atBot y, hbs.eventually_ge_atTop y] with n h1 h2
      have h : as n ≤ y ∧ y ≤ bs n := ⟨h1, h2⟩
      simp only [hG, dite_eq_left h]
      exact dimerTL_ii_Z (as n) (bs n) (ha n) (hb n) w W hw s ⟨y, h⟩
    by_cases hP : P y
    · simp only [hP, ite_true]
      exact (tendsto_norm.comp hGt).pow 2
    · simp only [hP, ite_false]
      exact tendsto_const_nhds
  · by_cases h : as n ≤ y ∧ y ≤ bs n
    · have h1 := DimerInfVolAux.norm_sq_evo_le (as n) (bs n) (ha n) (hb n) w s ⟨y, h⟩
      simp only [hG, dite_eq_left h]
      split_ifs
      · rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
        exact h1
      · rw [norm_zero]
        positivity
    · simp only [hG, dite_eq_right h, norm_zero]
      split_ifs
      · simp only [ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, norm_zero]
        positivity
      · rw [norm_zero]
        positivity

end DimerStrongAux

/-- **Thermodynamic limit of the tail weight** (paper: Remark `rem:dimerstrong`, by Lemma
`lem:dimerTL` as in the proof of Theorem `thm:dimer`(b)). -/
theorem dimer_TL_tail (J h : ℝ) (hJ : 0 < J) (r : ℤ) (Ns : ℕ → ℕ) (xs : ∀ n, Fin (Ns n))
    (hN : ∀ n, Even (Ns n)) (hpar : ∀ n, ((xs n).val : ℤ) % 2 = r)
    (hx1 : Tendsto (fun n => (xs n).val) atTop atTop)
    (hx2 : Tendsto (fun n => Ns n - (xs n).val) atTop atTop) (R : ℕ) (t : ℝ) :
    Tendsto (fun n => dimerTailLGT J h (Ns n) R t (xs n)) atTop (𝓝 (dimerTailZ J h r R t)) := by
  obtain ⟨has, hbs⟩ := DimerInfVolAux.tendsto_bounds Ns xs hx1 hx2
  have ha : ∀ n, -((xs n).val : ℤ) ≤ 0 := fun n => by omega
  have hb : ∀ n, 0 ≤ (Ns n : ℤ) - 1 - (xs n).val := fun n => by have := (xs n).isLt; omega
  have hD : ∀ n, dimerTailLGT J h (Ns n) R t (xs n) = ∫ η, (fun q : Cfg (Ns n) =>
      ∑ y : Fin (Ns n), (if (R : ℤ) < |(y.val : ℤ) - (xs n).val| then
        ‖evo (hq J h q) t y (xs n)‖ ^ 2 else 0))
        (DimerInfVolAux.dimerPhi (Ns n) (xs n).val r η) ∂muZ := fun n => by
    rw [dimerTailLGT_eq J h (hN n) R t (xs n)]
    exact (DimerInfVolAux.integral_comp_dimerPhi (Ns n) (hN n) (xs n).val r (hpar n) _).symm
  simp_rw [hD]
  unfold dimerTailZ
  refine tendsto_integral_of_dominated_convergence (fun _ => (1 : ℝ))
    (fun n => ((measurable_of_countable (fun q : Cfg (Ns n) => ∑ y : Fin (Ns n),
      (if (R : ℤ) < |(y.val : ℤ) - (xs n).val| then ‖evo (hq J h q) t y (xs n)‖ ^ 2
        else 0))).comp
      (DimerInfVolAux.measurable_dimerPhi (Ns n) (xs n).val r)).aestronglyMeasurable)
    (integrable_const _) (fun n => ae_of_all _ fun η => ?_) (ae_of_all _ fun η => ?_)
  · rw [DimerStrongAux.tail_dimerPhi J h hJ.ne' (Ns n) (xs n) r η R t, Real.norm_eq_abs,
      abs_of_nonneg (Finset.sum_nonneg fun y _ => by split_ifs <;> positivity)]
    exact DimerStrongAux.sum_ite_norm_sq_le_one _ (DimerTLAux.schrH_isSymm _ _ _) _ _ _
  · simp_rw [DimerStrongAux.tail_dimerPhi J h hJ.ne' (Ns _) (xs _) r η R t]
    exact DimerStrongAux.tendsto_ite_Z (dimerPotZ J h r η) |2 * h / J|
      (DimerInfVolAux.abs_dimerPotZ_le J h r η) _ _ ha hb has hbs _ (J * t)

namespace DimerStrongAux

/-- `η ↦ v^{(r)}(y)` is measurable. -/
lemma measurable_dimerPotZ_apply (J h : ℝ) (r y : ℤ) :
    Measurable (fun η : CfgZ => dimerPotZ J h r η y) :=
  measurable_const.mul ((measurable_of_countable sgn).comp (measurable_pi_apply _))

/-- The entries of the powers of `H_{v^{(r)}}` depend measurably on the signs. -/
lemma measurable_schrZpow (J h : ℝ) (r : ℤ) (k : ℕ) (x y : ℤ) :
    Measurable (fun η : CfgZ => schrZpow (dimerPotZ J h r η) k x y) := by
  induction k generalizing x with
  | zero => exact measurable_const
  | succ k ih =>
    simp only [schrZpow]
    refine Finset.measurable_sum _ fun z _ => ?_
    have hz : Measurable (fun η : CfgZ => schrZ (dimerPotZ J h r η) x z) := by
      by_cases hxz : x = z
      · subst hxz
        simp only [DimerInfVolAux.schrZ_diag]
        exact measurable_dimerPotZ_apply J h r x
      · simp only [schrZ, hxz, ite_false]
        exact measurable_const
    exact hz.mul (ih z)

/-- The entries of `e^{-isH_{v^{(r)}}}` on `ℤ` depend measurably on the signs. -/
lemma measurable_evoZw (J h : ℝ) (r : ℤ) (s : ℝ) (x y : ℤ) :
    Measurable (fun η : CfgZ => evoZw (dimerPotZ J h r η) s x y) := by
  refine measurable_of_tendsto_metrizable (f := fun (n : ℕ) (η : CfgZ) => ∑ k ∈ Finset.range n,
    ((-(s : ℂ) * Complex.I) ^ k / (k.factorial : ℂ)) *
      ((schrZpow (dimerPotZ J h r η) k x y : ℝ) : ℂ)) (fun n => ?_) ?_
  · exact Finset.measurable_sum _ fun k _ => measurable_const.mul
      (Complex.measurable_ofReal.comp (measurable_schrZpow J h r k x y))
  · rw [tendsto_pi_nhds]
    intro η
    exact (summable_evoZw (dimerPotZ J h r η) |2 * h / J|
      (DimerInfVolAux.abs_dimerPotZ_le J h r η) s x y).hasSum.tendsto_sum_nat

/-- `s ↦ ⟨δ_x, e^{-isH_w} δ_y⟩` is continuous on `ℤ` (locally uniform convergence of the series). -/
lemma continuous_evoZw (w : ℤ → ℝ) (W : ℝ) (hw : ∀ x, |w x| ≤ W) (x y : ℤ) :
    Continuous (fun s : ℝ => evoZw w s x y) := by
  have hW : 0 ≤ W := (abs_nonneg _).trans (hw 0)
  rw [continuous_iff_continuousAt]
  intro s₀
  have hcont : ContinuousOn (fun s : ℝ => evoZw w s x y)
      (Set.Ioo (-(|s₀| + 1)) (|s₀| + 1)) := by
    unfold evoZw
    refine continuousOn_tsum (u := fun k => ((2 + W) * (|s₀| + 1)) ^ k / (k.factorial : ℝ))
      (fun k => Continuous.continuousOn (by fun_prop)) (Real.summable_pow_div_factorial _)
      (fun k s hs => ?_)
    refine (DimerTLAux.norm_term_le s k _ (2 + W)
      (DimerInfVolAux.abs_schrZpow_le w W hw k x y)).trans ?_
    have hs' : |s| ≤ |s₀| + 1 := abs_le.mpr ⟨hs.1.le, hs.2.le⟩
    have h1 : (2 + W) * |s| ≤ (2 + W) * (|s₀| + 1) :=
      mul_le_mul_of_nonneg_left hs' (by linarith)
    exact div_le_div_of_nonneg_right
      (pow_le_pow_left₀ (mul_nonneg (by linarith) (abs_nonneg s)) h1 k) (by positivity)
  exact hcont.continuousAt (Ioo_mem_nhds (by linarith [neg_abs_le s₀])
    (by linarith [le_abs_self s₀]))

/-- Partial column weights of `e^{-isH_w}` on `ℤ` are at most `1`. -/
lemma tsum_ite_le_one (w : ℤ → ℝ) (W : ℝ) (hw : ∀ x, |w x| ≤ W) (P : ℤ → Prop)
    [DecidablePred P] (s : ℝ) : ∑' y : ℤ, (if P y then ‖evoZw w s y 0‖ ^ 2 else 0) ≤ 1 := by
  have hlim := tendsto_ite_Z w W hw (fun n : ℕ => -(n : ℤ)) (fun n => (n : ℤ))
    (fun n => by simp) (fun n => by positivity)
    (tendsto_neg_atTop_atBot.comp tendsto_natCast_atTop_atTop) tendsto_natCast_atTop_atTop P s
  exact le_of_tendsto' hlim fun n =>
    sum_ite_norm_sq_le_one _ (DimerTLAux.schrH_isSymm _ _ _) _ _ _

/-- Chebyshev: the weight beyond `R ≥ 1` is at most the `p`-th moment divided by `R^p`. -/
lemma tsum_tail_le (u : ℤ → ℝ) (hu : ∀ y, 0 ≤ u y) (p : ℝ) (hp : 0 < p) (R : ℕ) (hR : 1 ≤ R)
    (B : ℝ) (hB : ∀ S : Finset ℤ, ∑ y ∈ S, |(y : ℝ)| ^ p * u y ≤ B) :
    ∑' y : ℤ, (if (R : ℤ) < |y| then u y else 0) ≤ B / (R : ℝ) ^ p := by
  have hR0 : (0 : ℝ) < R := by exact_mod_cast hR
  have hRp : 0 < (R : ℝ) ^ p := Real.rpow_pos_of_pos hR0 p
  refine Real.tsum_le_of_sum_le (fun y => ?_) fun S => ?_
  · show (0 : ℝ) ≤ if (R : ℤ) < |y| then u y else 0
    split_ifs
    · exact hu y
    · exact le_rfl
  rw [le_div_iff₀ hRp, Finset.sum_mul]
  refine le_trans (Finset.sum_le_sum fun y _ => ?_) (hB S)
  split_ifs with hy
  · have h1 : (R : ℝ) ≤ |(y : ℝ)| := by exact_mod_cast hy.le
    have h2 : (R : ℝ) ^ p ≤ |(y : ℝ)| ^ p := Real.rpow_le_rpow hR0.le h1 hp.le
    rw [mul_comm]
    exact mul_le_mul_of_nonneg_right h2 (hu y)
  · rw [zero_mul]
    exact mul_nonneg (by positivity) (hu y)

end DimerStrongAux

/-- **Tightness from dynamical localization** (paper: Remark `rem:dimerstrong`): if almost surely
`sup_t ∑_y |y|^p |⟨δ_y, e^{-itH}δ_0⟩|² < ∞` for some `p > 0`, then `sup_t dimerTailZ R t → 0`. -/
theorem dimerTailZ_tendsto_zero (J h : ℝ) (hJ : 0 < J) (r : ℤ) (p : ℝ) (hp : 0 < p)
    (hDL : ∀ᵐ η ∂muZ, ∃ C : ℝ, ∀ (t : ℝ) (S : Finset ℤ),
      ∑ y ∈ S, |(y : ℝ)| ^ p * ‖evoZw (dimerPotZ J h r η) t y 0‖ ^ 2 ≤ C) :
    Tendsto (fun R : ℕ => ⨆ t : ℝ, dimerTailZ J h r R t) atTop (𝓝 0) := by
  have _ := hJ  -- the sign of `J` is not needed here (`hDL` holds for all times)
  -- the moment supremum over rational times and finite windows
  set A : CfgZ → ENNReal := fun η => ⨆ (q : ℚ) (S : Finset ℤ),
    ENNReal.ofReal (∑ y ∈ S, |(y : ℝ)| ^ p * ‖evoZw (dimerPotZ J h r η) (q : ℝ) y 0‖ ^ 2)
    with hA
  have hAm : Measurable A := by
    refine Measurable.iSup fun q => Measurable.iSup fun S => ?_
    exact ENNReal.measurable_ofReal.comp (Finset.measurable_sum _ fun y _ =>
      measurable_const.mul (((DimerStrongAux.measurable_evoZw J h r _ y 0).norm).pow_const 2))
  -- the dominating functions `G_R = min(1, A / R^p)`
  set G : ℕ → CfgZ → ℝ := fun R η => min 1 ((A η).toReal / (R : ℝ) ^ p) with hG
  have hG0 : ∀ R η, 0 ≤ G R η := fun R η =>
    le_min zero_le_one (div_nonneg ENNReal.toReal_nonneg (Real.rpow_nonneg (Nat.cast_nonneg R) p))
  have hG1 : ∀ R η, G R η ≤ 1 := fun R η => min_le_left _ _
  have hGn : ∀ R η, ‖G R η‖ ≤ 1 := fun R η => by
    rw [Real.norm_eq_abs, abs_of_nonneg (hG0 R η)]
    exact hG1 R η
  have hGm : ∀ R, Measurable (G R) := fun R =>
    measurable_const.min ((ENNReal.measurable_toReal.comp hAm).div_const _)
  have hGi : ∀ R, Integrable (G R) muZ := fun R =>
    Integrable.mono' (integrable_const 1) (hGm R).aestronglyMeasurable (ae_of_all _ (hGn R))
  -- `∫ G_R → 0` by dominated convergence
  have hlimG : Tendsto (fun R : ℕ => ∫ η, G R η ∂muZ) atTop (𝓝 0) := by
    have hpt : ∀ η, Tendsto (fun R : ℕ => G R η) atTop (𝓝 0) := fun η => by
      have h1 : Tendsto (fun R : ℕ => (A η).toReal / (R : ℝ) ^ p) atTop (𝓝 0) :=
        tendsto_const_nhds.div_atTop ((tendsto_rpow_atTop hp).comp tendsto_natCast_atTop_atTop)
      exact squeeze_zero (fun R => hG0 R η) (fun R => min_le_right _ _) h1
    have h := tendsto_integral_of_dominated_convergence (μ := muZ) (F := G)
      (f := fun _ => (0 : ℝ)) (fun _ => (1 : ℝ)) (fun R => (hGm R).aestronglyMeasurable) (integrable_const 1)
      (fun R => ae_of_all _ (hGn R)) (ae_of_all _ hpt)
    simpa using h
  -- almost surely, the tail beyond `R ≥ 1` is at most `G_R`, for all times
  have hbound : ∀ R : ℕ, 1 ≤ R → ∀ t : ℝ, ∀ᵐ η ∂muZ, ∑' y : ℤ,
      (if (R : ℤ) < |y| then ‖evoZw (dimerPotZ J h r η) (J * t) y 0‖ ^ 2 else 0) ≤ G R η := by
    intro R hR t
    filter_upwards [hDL] with η hη
    obtain ⟨C, hC⟩ := hη
    have hW := DimerInfVolAux.abs_dimerPotZ_le J h r η
    have hAC : A η ≤ ENNReal.ofReal C :=
      iSup₂_le fun q S => ENNReal.ofReal_le_ofReal (hC q S)
    have hAt : A η ≠ ⊤ := ne_top_of_le_ne_top ENNReal.ofReal_ne_top hAC
    have hmom : ∀ (s : ℝ) (S : Finset ℤ),
        ∑ y ∈ S, |(y : ℝ)| ^ p * ‖evoZw (dimerPotZ J h r η) s y 0‖ ^ 2 ≤ (A η).toReal := by
      intro s S
      have hcl : IsClosed {s : ℝ | ∑ y ∈ S, |(y : ℝ)| ^ p *
          ‖evoZw (dimerPotZ J h r η) s y 0‖ ^ 2 ≤ (A η).toReal} :=
        isClosed_le (continuous_finsetSum _ fun y _ => continuous_const.mul
          (((DimerStrongAux.continuous_evoZw _ _ hW y 0).norm).pow 2)) continuous_const
      refine isClosed_property Rat.denseRange_cast hcl (fun q => ?_) s
      show ∑ y ∈ S, |(y : ℝ)| ^ p * ‖evoZw (dimerPotZ J h r η) (q : ℝ) y 0‖ ^ 2 ≤ (A η).toReal
      rw [← ENNReal.ofReal_le_iff_le_toReal hAt]
      exact le_iSup₂ (f := fun (q : ℚ) (S : Finset ℤ) => ENNReal.ofReal
        (∑ y ∈ S, |(y : ℝ)| ^ p * ‖evoZw (dimerPotZ J h r η) (q : ℝ) y 0‖ ^ 2)) q S
    refine le_min (DimerStrongAux.tsum_ite_le_one _ _ hW _ _) ?_
    exact DimerStrongAux.tsum_tail_le _ (fun y => by positivity) p hp R hR _ (hmom (J * t))
  -- conclusion
  have hnn : ∀ (R : ℕ) (t : ℝ) (η : CfgZ), 0 ≤ ∑' y : ℤ,
      (if (R : ℤ) < |y| then ‖evoZw (dimerPotZ J h r η) (J * t) y 0‖ ^ 2 else 0) :=
    fun R t η => tsum_nonneg fun y => by split_ifs <;> positivity
  refine squeeze_zero' (Eventually.of_forall fun R => Real.iSup_nonneg fun t =>
    integral_nonneg fun η => hnn R t η) ?_ hlimG
  filter_upwards [eventually_ge_atTop 1] with R hR
  refine ciSup_le fun t => ?_
  unfold dimerTailZ
  exact integral_mono_of_nonneg (ae_of_all _ fun η => hnn R t η) (hGi R) (hbound R hR t)

/-- **The dressed fermion stays tight** (paper: Remark `rem:dimerstrong`): given dynamical localization
of the random dimer model, `lim_{R→∞} sup_t lim_{N→∞} ⟨Ψ^{x_N}_N(t), ∑_{|y-x_N|>R} n_y Ψ^{x_N}_N(t)⟩ = 0`. -/
theorem dimer_tight_of_dynLoc (J h : ℝ) (hJ : 0 < J) (r : ℤ) (Ns : ℕ → ℕ) (xs : ∀ n, Fin (Ns n))
    (hN : ∀ n, Even (Ns n)) (hpar : ∀ n, ((xs n).val : ℤ) % 2 = r)
    (hx1 : Tendsto (fun n => (xs n).val) atTop atTop)
    (hx2 : Tendsto (fun n => Ns n - (xs n).val) atTop atTop) (p : ℝ) (hp : 0 < p)
    (hDL : ∀ᵐ η ∂muZ, ∃ C : ℝ, ∀ (t : ℝ) (S : Finset ℤ),
      ∑ y ∈ S, |(y : ℝ)| ^ p * ‖evoZw (dimerPotZ J h r η) t y 0‖ ^ 2 ≤ C) :
    ∃ L : ℕ → ℝ → ℝ, (∀ R t, Tendsto (fun n => dimerTailLGT J h (Ns n) R t (xs n)) atTop
        (𝓝 (L R t))) ∧ Tendsto (fun R : ℕ => ⨆ t : ℝ, L R t) atTop (𝓝 0) :=
  ⟨dimerTailZ J h r, fun R t => dimer_TL_tail J h hJ r Ns xs hN hpar hx1 hx2 R t,
    dimerTailZ_tendsto_zero J h hJ r p hp hDL⟩

end DFL
