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

/-- **Tightness from dynamical localization** (paper: Remark `rem:dimerstrong`): if almost surely
`sup_t ∑_y |y|^p |⟨δ_y, e^{-itH}δ_0⟩|² < ∞` for some `p > 0`, then `sup_t dimerTailZ R t → 0`. -/
theorem dimerTailZ_tendsto_zero (J h : ℝ) (hJ : 0 < J) (r : ℤ) (p : ℝ) (hp : 0 < p)
    (hDL : ∀ᵐ η ∂muZ, ∃ C : ℝ, ∀ (t : ℝ) (S : Finset ℤ),
      ∑ y ∈ S, |(y : ℝ)| ^ p * ‖evoZw (dimerPotZ J h r η) t y 0‖ ^ 2 ≤ C) :
    Tendsto (fun R : ℕ => ⨆ t : ℝ, dimerTailZ J h r R t) atTop (𝓝 0) := by
  sorry

/-- **The dressed fermion stays tight** (paper: Remark `rem:dimerstrong`): given dynamical localization
of the random dimer model, `lim_{R→∞} sup_t lim_{N→∞} ⟨Ψ^{x_N}_N(t), ∑_{|y-x_N|>R} n_y Ψ^{x_N}_N(t)⟩ = 0`. -/
theorem dimer_tight_of_dynLoc (J h : ℝ) (hJ : 0 < J) (r : ℤ) (Ns : ℕ → ℕ) (xs : ∀ n, Fin (Ns n))
    (hN : ∀ n, Even (Ns n)) (hpar : ∀ n, ((xs n).val : ℤ) % 2 = r)
    (hx1 : Tendsto (fun n => (xs n).val) atTop atTop)
    (hx2 : Tendsto (fun n => Ns n - (xs n).val) atTop atTop) (p : ℝ) (hp : 0 < p)
    (hDL : ∀ᵐ η ∂muZ, ∃ C : ℝ, ∀ (t : ℝ) (S : Finset ℤ),
      ∑ y ∈ S, |(y : ℝ)| ^ p * ‖evoZw (dimerPotZ J h r η) t y 0‖ ^ 2 ≤ C) :
    ∃ L : ℕ → ℝ → ℝ, (∀ R t, Tendsto (fun n => dimerTailLGT J h (Ns n) R t (xs n)) atTop
        (𝓝 (L R t))) ∧ Tendsto (fun R : ℕ => ⨆ t : ℝ, L R t) atTop (𝓝 0) := by
  sorry

end DFL
