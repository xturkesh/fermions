import DFL.InfVol

/-!
# Infinite volume: the transported charge (paper: Theorem `thm:infvol`(b))

For an occupation configuration `n ∈ {0,1}^ℤ` (occupied set `S`, possibly infinite) and a cut at
`m ∈ ℤ`, the transported charge of the gauge theory on `Λ_M = [-M, M]`,
`⟨Ψ^{(M)}(t), N_{[-M, m]} Ψ^{(M)}(t)⟩ - |S ∩ [-M, m]|`, converges as `M → ∞`, and under Hypothesis (EC)
the limit is bounded uniformly in time: in finite volume it is the difference of the charge flowing
across the cut in the two directions, `∑_{z ∈ S, z > m} ∑_{x ≤ m} E|U_{xz}|² - ∑_{z ∈ S, z ≤ m} ∑_{x > m} E|U_{xz}|²`,
each term converges by `tendsto_evoZ`, and the double sums are dominated by `∑_d d F(d) < ∞`.
-/

open MeasureTheory Filter Topology Matrix Finset

noncomputable section

namespace DFL

/-- The occupied sites of `n` inside `Λ_M`, relabelled in `Fin (2M+1)`. -/
def occM (M : ℕ) (n : CfgZ) : Finset (Fin (2 * M + 1)) :=
  univ.filter (fun i => n ((i.val : ℤ) - M) = 1)

/-- The sites `[-M, m]` of `Λ_M`, relabelled in `Fin (2M+1)`. -/
def leftM (M : ℕ) (m : ℤ) : Finset (Fin (2 * M + 1)) :=
  univ.filter (fun i => (i.val : ℤ) - M ≤ m)

/-- The transported charge on `Λ_M`: `⟨Ψ^{(M)}(t), N_{[-M,m]} Ψ^{(M)}(t)⟩ - |S ∩ [-M, m]|`, for the gauge
theory with Fock matter `n|_{Λ_M}` and `z`-polarized links. -/
def transportM (J h : ℝ) (n : CfgZ) (m : ℤ) (t : ℝ) (M : ℕ) : ℝ :=
  (expect (matter (numOp (leftM M m))) (lgtState J h (fockVec (occM M n)) t)).re -
    ((occM M n ∩ leftM M m).card : ℝ)

namespace InfVolTransportAux

/-! ### The cut relative to `Λ_M` -/

/-- If the cut `m` is not to the left of `Λ_M`, then `leftM M m` is the left segment
`leftOf (m + M)` of `Fin (2M+1)`. -/
lemma leftM_eq_leftOf (M : ℕ) (m : ℤ) (hm : 0 ≤ m + M) :
    leftM M m = leftOf (N := 2 * M + 1) (m + M).toNat := by
  ext i
  simp only [leftM, leftOf, mem_filter, mem_univ, true_and]
  omega

/-- If the cut `m` is to the left of `Λ_M`, then `leftM M m` is empty. -/
lemma leftM_eq_empty (M : ℕ) (m : ℤ) (hm : m + M < 0) : leftM M m = ∅ := by
  ext i
  simp only [leftM, mem_filter, mem_univ, true_and]
  simp only [Finset.notMem_empty, iff_false, not_le]
  omega

/-- The number operator of the empty region vanishes. -/
lemma numOp_empty (N : ℕ) : numOp (∅ : Finset (Fin N)) = 0 := by
  have hc : count (∅ : Finset (Fin N)) = fun _ => 0 := by
    funext n; simp [count]
  rw [numOp, hc, diagonal_zero]

end InfVolTransportAux

/-- The finite-volume transported charge is bounded uniformly in `M` and `t` under (EC). -/
theorem abs_transportM_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (n : CfgZ) (m : ℤ) (t : ℝ)
    (M : ℕ) : |transportM J h n m t M| ≤ 2 * ∑' d : ℕ, (d : ℝ) * F d := by
  have hC : 0 ≤ 2 * ∑' d : ℕ, (d : ℝ) * F d :=
    mul_nonneg (by norm_num) (tsum_nonneg fun d => mul_nonneg (Nat.cast_nonneg d) (hEC.nonneg d))
  by_cases hm : 0 ≤ m + M
  · unfold transportM
    rw [InfVolTransportAux.leftM_eq_leftOf M m hm]
    exact lgt_transport J h hEC _ _ t
  · have h0 : transportM J h n m t M = 0 := by
      unfold transportM
      rw [InfVolTransportAux.leftM_eq_empty M m (by omega), InfVolTransportAux.numOp_empty]
      simp [expect, cplx]
    rw [h0, abs_zero]
    exact hC

namespace InfVolTransportAux

variable {N : ℕ}

/-! ### Bernoulli averages (as in `DFL.Main`, where these lemmas are private) -/

/-- There are `2^N` charge configurations. -/
lemma card_cfg : Fintype.card (Cfg N) = 2 ^ N := by
  simp [Cfg, ZMod.card]

/-- Average of a constant. -/
lemma avgQ_const (c : ℝ) : avgQ (fun _ : Cfg N => c) = c := by
  simp only [avgQ, Finset.sum_const, Finset.card_univ, card_cfg, nsmul_eq_mul]
  push_cast
  field_simp

/-- The average commutes with finite sums. -/
lemma avgQ_sum {κ : Type*} (s : Finset κ) (f : κ → Cfg N → ℝ) :
    avgQ (fun q => ∑ i ∈ s, f i q) = ∑ i ∈ s, avgQ (f i) := by
  simp only [avgQ]
  rw [Finset.sum_comm, Finset.mul_sum]

/-- Monotonicity of the average. -/
lemma avgQ_mono {f g : Cfg N → ℝ} (h : ∀ q, f q ≤ g q) : avgQ f ≤ avgQ g := by
  simp only [avgQ]
  exact mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun q _ => h q) (by positivity)

/-- Real part of a complex Bernoulli average of real numbers. -/
lemma re_avg (f : Cfg N → ℝ) :
    ((1 / 2 ^ N : ℂ) * ∑ q : Cfg N, ((f q : ℝ) : ℂ)).re = avgQ f := by
  simp only [avgQ]
  rw [← Complex.ofReal_sum]
  have : (1 / 2 ^ N : ℂ) = ((1 / 2 ^ N : ℝ) : ℂ) := by push_cast; ring
  rw [this, ← Complex.ofReal_mul, Complex.ofReal_re]

/-! ### Free fermions in a fixed sector (as in `DFL.Main`) -/

section Aux

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Expectation of a complexified real diagonal matrix: `⟨v, diag(d) v⟩ = ∑ d_i |v_i|²`. -/
lemma expect_cplx_diagonal (d : ι → ℝ) (v : ι → ℂ) :
    expect (cplx (diagonal d)) v = ((∑ i, d i * Complex.normSq (v i) : ℝ) : ℂ) := by
  have hc : cplx (diagonal d) = diagonal (fun i => (d i : ℂ)) := by
    unfold cplx; exact Matrix.diagonal_map (by simp)
  rw [expect, hc]
  simp only [dotProduct, mulVec_diagonal, Pi.star_apply, Complex.ofReal_sum, Complex.ofReal_mul]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [Complex.normSq_eq_conj_mul_self]
  simp only [Complex.star_def]
  ring

/-- Columns of a unitary matrix are normalized. -/
lemma sum_normSq_col (U : Matrix ι ι ℂ) (hU : Uᴴ * U = 1) (z : ι) :
    ∑ x, Complex.normSq (U x z) = 1 := by
  have h := congrFun (congrFun hU z) z
  simp only [one_apply_eq] at h
  have : ((∑ x, Complex.normSq (U x z) : ℝ) : ℂ) = 1 := by
    rw [← h, Complex.ofReal_sum]
    refine Finset.sum_congr rfl (fun x _ => ?_)
    rw [Complex.normSq_eq_conj_mul_self]
    rfl
  exact_mod_cast this

end Aux

/-- `⟨n_x⟩` of an evolved Fock state, as a real sum over configurations. -/
lemma sum_occ_normSq (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (x : Fin N) (t : ℝ) :
    ∑ n, occ (n x) * Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) =
      ∑ z ∈ S, Complex.normSq (evo A t x z) := by
  have := density A hA S x t
  rw [expect_cplx_diagonal] at this
  exact_mod_cast this

/-- `⟨N_R⟩` of an evolved Fock state: `∑_{x ∈ R} ∑_{z ∈ S} |U_{xz}|²`. -/
lemma sum_count_normSq (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S R : Finset (Fin N))
    (t : ℝ) :
    ∑ n, count R n * Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) =
      ∑ x ∈ R, ∑ z ∈ S, Complex.normSq (evo A t x z) := by
  simp only [count, Finset.sum_mul]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl (fun x _ => sum_occ_normSq A hA S x t)

/-- LGT expectation of a diagonal matter observable, as a Bernoulli average of sector sums. -/
lemma lgt_expect_diag_re (J h : ℝ) (ψ : Cfg N → ℂ) (d : Cfg N → ℝ) (t : ℝ) :
    (expect (matter (diagonal d)) (lgtState J h ψ t)).re =
      avgQ (fun q => ∑ n, d n * Complex.normSq ((evo (dGamma (hq J h q)) t *ᵥ ψ) n)) := by
  rw [matter, sector_average_matter]
  simp_rw [expect_cplx_diagonal]
  exact re_avg _

/-- `⟨Ψ(t), N_L Ψ(t)⟩ = ∑_{x ∈ L} ∑_{z ∈ S} E|U_{xz}|²` for the LGT started from a Fock state. -/
lemma lgt_numOp_re (J h : ℝ) (L S : Finset (Fin N)) (t : ℝ) :
    (expect (matter (numOp L)) (lgtState J h (fockVec S) t)).re =
      ∑ x ∈ L, ∑ z ∈ S, avgQ (fun q => Complex.normSq (evo (hq J h q) t x z)) := by
  rw [numOp, lgt_expect_diag_re]
  have e1 : (fun q => ∑ n, count L n *
      Complex.normSq ((evo (dGamma (hq J h q)) t *ᵥ fockVec S) n)) =
      fun q => ∑ x ∈ L, ∑ z ∈ S, Complex.normSq (evo (hq J h q) t x z) := by
    funext q; exact sum_count_normSq _ (hq_isSymm J h q) S _ t
  rw [e1, avgQ_sum]
  exact Finset.sum_congr rfl fun x _ => avgQ_sum _ _

/-- The averaged transition probabilities `E|U_{xz}|²` have unit column sums. -/
lemma sum_avgQ_normSq_col (J h t : ℝ) (z : Fin N) :
    ∑ x, avgQ (fun q : Cfg N => Complex.normSq (evo (hq J h q) t x z)) = 1 := by
  rw [← avgQ_sum]
  have : (fun q : Cfg N => ∑ x ∈ univ, Complex.normSq (evo (hq J h q) t x z)) = fun _ => 1 := by
    funext q; exact sum_normSq_col _ (evo_unitary _ (hq_isSymm J h q) t) z
  rw [this, avgQ_const]

/-- With unit column sums, `∑_{i ∈ L} ∑_{k ∈ S} a_{ik} - |S ∩ L|` is the weight flowing into `L`
from `S \ L` minus the weight flowing out of `L` from `S ∩ L`. -/
lemma col_identity {ι : Type*} [Fintype ι] [DecidableEq ι] (a : ι → ι → ℝ)
    (hcol : ∀ k, ∑ i, a i k = 1) (L S : Finset ι) :
    ∑ i ∈ L, ∑ k ∈ S, a i k - ((S ∩ L).card : ℝ) =
      ∑ i, ∑ k, (if i ∈ L ∧ k ∉ L ∧ k ∈ S then a i k else 0) -
        ∑ i, ∑ k, (if i ∉ L ∧ k ∈ L ∧ k ∈ S then a i k else 0) := by
  have hsub : ∀ (T : Finset ι) (f : ι → ℝ), ∑ i ∈ T, f i = ∑ i, if i ∈ T then f i else 0 := by
    intro T f; rw [Finset.sum_ite_mem, univ_inter]
  have h1 : ∑ i ∈ L, ∑ k ∈ S, a i k = ∑ i, ∑ k, (if i ∈ L ∧ k ∈ S then a i k else 0) := by
    rw [hsub L]
    refine Finset.sum_congr rfl fun i _ => ?_
    by_cases hi : i ∈ L
    · simp only [hi, true_and, ite_true]
      exact hsub S _
    · simp [hi]
  have h2 : ((S ∩ L).card : ℝ) = ∑ i, ∑ k, (if k ∈ S ∧ k ∈ L then a i k else 0) := by
    rw [Finset.sum_comm]
    have e : ((S ∩ L).card : ℝ) = ∑ k ∈ S ∩ L, ∑ i, a i k := by
      rw [Finset.sum_congr rfl fun k _ => hcol k]
      simp
    rw [e, hsub (S ∩ L)]
    refine Finset.sum_congr rfl fun k _ => ?_
    by_cases hk : k ∈ S ∧ k ∈ L
    · simp only [Finset.mem_inter, hk, and_self, ite_true]
    · rw [ite_eq_right (by rwa [Finset.mem_inter]), Finset.sum_eq_zero fun i _ => ite_eq_right hk]
  rw [h1, h2, ← Finset.sum_sub_distrib, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [← Finset.sum_sub_distrib, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl fun k _ => ?_
  by_cases hi : i ∈ L <;> by_cases hk : k ∈ L <;> by_cases hs : k ∈ S <;> simp [hi, hk, hs]

/-! ### The law of `q|_{Λ_M}` under `muZ` is uniform -/

/-- The restriction to `Λ_M` is measurable. -/
lemma measurable_restrictZ (M : ℕ) : Measurable (restrictZ M) :=
  Measurable.of_eval fun _ => measurable_pi_apply _

/-- The cylinder `{q : q|_{Λ_M} = σ}` has probability `2^{-(2M+1)}`. -/
lemma muZ_restrictZ_preimage (M : ℕ) (σ : Cfg (2 * M + 1)) :
    muZ (restrictZ M ⁻¹' {σ}) = (2⁻¹ : ENNReal) ^ (2 * M + 1) := by
  classical
  let tt : ℤ → Set (ZMod 2) := fun j => if hj : |j| ≤ M then {σ (idxZ M j hj)} else Set.univ
  have hset : restrictZ M ⁻¹' {σ} = Set.pi (Finset.Icc (-(M : ℤ)) M : Set ℤ) tt := by
    ext q
    simp only [Set.mem_preimage, Set.mem_singleton_iff, Set.mem_pi, Finset.coe_Icc, Set.mem_Icc]
    constructor
    · rintro rfl j hj
      have hj' : |j| ≤ M := abs_le.mpr hj
      simp only [tt, dite_eq_left hj', Set.mem_singleton_iff, restrictZ, idxZ]
      congr 1
      omega
    · intro H
      funext i
      have hi : |((i.val : ℤ) - M)| ≤ M := by
        have := i.isLt
        rw [abs_le]; omega
      have := H ((i.val : ℤ) - M) (abs_le.mp hi)
      simp only [tt, dite_eq_left hi, Set.mem_singleton_iff] at this
      show q ((i.val : ℤ) - M) = σ i
      rw [this]
      congr 1
      ext
      simp [idxZ]
  rw [hset, muZ, Measure.infinitePi_pi _ (fun _ _ => MeasurableSet.of_discrete)]
  have h2 : ∀ j ∈ Finset.Icc (-(M : ℤ)) M, unifZ2 (tt j) = 2⁻¹ := by
    intro j hj
    have hj' : |j| ≤ M := abs_le.mpr (Finset.mem_Icc.mp hj)
    simp only [tt, dite_eq_left hj', unifZ2, smul_apply, Measure.count_singleton, smul_eq_mul,
      mul_one]
  rw [Finset.prod_congr rfl h2, Finset.prod_const, Int.card_Icc]
  congr 1
  omega

/-- The law of `q|_{Λ_M}` under `muZ` is uniform: `E f(q|_{Λ_M}) = avgQ f`. -/
lemma integral_restrictZ (M : ℕ) (f : Cfg (2 * M + 1) → ℝ) :
    ∫ q, f (restrictZ M q) ∂muZ = avgQ f := by
  have hmeas : ∀ σ : Cfg (2 * M + 1), MeasurableSet (restrictZ M ⁻¹' {σ}) :=
    fun σ => measurable_restrictZ M (measurableSet_singleton σ)
  have hdecomp : (fun q => f (restrictZ M q)) =
      fun q => ∑ σ, (restrictZ M ⁻¹' {σ}).indicator (fun _ => f σ) q := by
    funext q
    rw [Finset.sum_eq_single (restrictZ M q)]
    · simp
    · intro b _ hb
      rw [Set.indicator_of_notMem]
      simpa [eq_comm] using hb
    · simp
  rw [hdecomp, integral_finsetSum]
  · simp_rw [integral_indicator_const _ (hmeas _), measureReal_def, muZ_restrictZ_preimage]
    rw [avgQ, Finset.mul_sum]
    refine Finset.sum_congr rfl fun σ _ => ?_
    simp [ENNReal.toReal_pow, ENNReal.toReal_inv]
  · intro σ _
    exact (integrable_const _).indicator (hmeas σ)

/-! ### Averaged transition probabilities on `Λ_M` and on `ℤ` -/

/-- `E |U^{(M)}_q(t)_{xz}|²` (`0` if `x` or `z` is outside `Λ_M`). -/
def aM (J h t : ℝ) (M : ℕ) (x z : ℤ) : ℝ := ∫ q, Complex.normSq (evoZ J h q t M x z) ∂muZ

/-- `E |U^ℤ_q(t)_{xz}|²`. -/
def aZ (J h t : ℝ) (x z : ℤ) : ℝ := ∫ q, Complex.normSq (UZ J h q t x z) ∂muZ

/-- The finite-volume entries are bounded by one. -/
lemma norm_evoZ_le_one (J h : ℝ) (q : CfgZ) (t : ℝ) (M : ℕ) (x z : ℤ) :
    ‖evoZ J h q t M x z‖ ≤ 1 := by
  unfold evoZ
  split_ifs
  · exact norm_evo_le_one _ (hq_isSymm J h _) t _ _
  · simp

/-- The finite-volume entries depend measurably on the charges (they factor through the
restriction to `Λ_M`, a map into a finite discrete space). -/
lemma measurable_evoZ' (J h t : ℝ) (M : ℕ) (x z : ℤ) :
    Measurable (fun q : CfgZ => evoZ J h q t M x z) := by
  unfold evoZ
  split_ifs with hxz
  · exact (Measurable.of_discrete (f := fun σ : Cfg (2 * M + 1) =>
      evo (hq J h σ) t (idxZ M x hxz.1) (idxZ M z hxz.2))).comp (measurable_restrictZ M)
  · exact measurable_const

/-- `E |U^{(M)}_q(t)_{xz}|² → E |U^ℤ_q(t)_{xz}|²` (dominated convergence). -/
lemma tendsto_aM (J h t : ℝ) (x z : ℤ) :
    Tendsto (fun M => aM J h t M x z) atTop (𝓝 (aZ J h t x z)) := by
  refine tendsto_integral_of_dominated_convergence (fun _ => 1) (fun M => ?_) (integrable_const 1)
    (fun M => Filter.Eventually.of_forall fun q => ?_) (Filter.Eventually.of_forall fun q => ?_)
  · exact (Complex.continuous_normSq.measurable.comp
      (measurable_evoZ' J h t M x z)).aestronglyMeasurable
  · rw [Real.norm_of_nonneg (Complex.normSq_nonneg _), Complex.normSq_eq_norm_sq]
    have h1 := norm_evoZ_le_one J h q t M x z
    nlinarith [norm_nonneg (evoZ J h q t M x z)]
  · exact (Complex.continuous_normSq.tendsto _).comp (tendsto_evoZ J h q t x z)

/-- The sites of `Λ_M`, as integers, lie in `[-M, M]`. -/
lemma abs_val_sub_le (M : ℕ) (i : Fin (2 * M + 1)) : |((i.val : ℤ) - M)| ≤ M := by
  have := i.isLt
  rw [abs_le]; omega

/-- `idxZ` inverts the relabelling `i ↦ i - M`. -/
lemma idxZ_val_sub (M : ℕ) (i : Fin (2 * M + 1)) (hi : |((i.val : ℤ) - M)| ≤ M) :
    idxZ M ((i.val : ℤ) - M) hi = i := by
  ext; simp [idxZ]

/-- On `Λ_M`, `aM` is the Bernoulli average of the finite-volume transition probability. -/
lemma aM_eq_avgQ (J h t : ℝ) (M : ℕ) (x z : ℤ) (hx : |x| ≤ M) (hz : |z| ≤ M) :
    aM J h t M x z = avgQ (fun σ : Cfg (2 * M + 1) =>
      Complex.normSq (evo (hq J h σ) t (idxZ M x hx) (idxZ M z hz))) := by
  rw [← integral_restrictZ]
  unfold aM evoZ
  simp only [dite_eq_left (And.intro hx hz)]

/-- Outside `Λ_M`, `aM` vanishes. -/
lemma aM_eq_zero (J h t : ℝ) (M : ℕ) (x z : ℤ) (hxz : ¬ (|x| ≤ M ∧ |z| ≤ M)) :
    aM J h t M x z = 0 := by
  unfold aM evoZ
  simp [dite_eq_right hxz]

/-- `aM` is non-negative. -/
lemma aM_nonneg (J h t : ℝ) (M : ℕ) (x z : ℤ) : 0 ≤ aM J h t M x z :=
  integral_nonneg fun _ => Complex.normSq_nonneg _

/-- Under (EC), `aM J h t M x z ≤ F(|x - z|)`, uniformly in `M`. -/
lemma aM_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (t : ℝ) (M : ℕ) (x z : ℤ) :
    aM J h t M x z ≤ F (x - z).natAbs := by
  by_cases hxz : |x| ≤ M ∧ |z| ≤ M
  · rw [aM_eq_avgQ J h t M x z hxz.1 hxz.2]
    set i := idxZ M x hxz.1
    set k := idxZ M z hxz.2
    have hik : Nat.dist i.val k.val = (x - z).natAbs := by
      have := abs_le.mp hxz.1
      have := abs_le.mp hxz.2
      simp only [i, k, idxZ, Nat.dist]
      omega
    calc avgQ (fun σ : Cfg (2 * M + 1) => Complex.normSq (evo (hq J h σ) t i k))
        ≤ avgQ (fun σ => corrQ (hq J h σ) (hq_isHermitian J h σ) i k) := by
          refine avgQ_mono fun σ => ?_
          rw [Complex.normSq_eq_norm_sq]
          have h1 := norm_evo_le_one _ (hq_isSymm J h σ) t i k
          have h2 := norm_evo_le_corrQ _ (hq_isHermitian J h σ) t i k
          nlinarith [norm_nonneg (evo (hq J h σ) t i k)]
      _ ≤ F (Nat.dist i.val k.val) := hEC.bound _ i k
      _ = F (x - z).natAbs := by rw [hik]
  · rw [aM_eq_zero J h t M x z hxz]; exact hEC.nonneg _

/-! ### The two crossing sums -/

/-- Pairs `(x, z)` with `x ≤ m < z` and `z` occupied: weight carried into `(-∞, m]`. -/
def G1 (n : CfgZ) (m : ℤ) (a : ℤ → ℤ → ℝ) (p : ℤ × ℤ) : ℝ :=
  if p.1 ≤ m ∧ m < p.2 ∧ n p.2 = 1 then a p.1 p.2 else 0

/-- Pairs `(x, z)` with `z ≤ m < x` and `z` occupied: weight carried out of `(-∞, m]`. -/
def G2 (n : CfgZ) (m : ℤ) (a : ℤ → ℤ → ℝ) (p : ℤ × ℤ) : ℝ :=
  if m < p.1 ∧ p.2 ≤ m ∧ n p.2 = 1 then a p.1 p.2 else 0

/-- A family on `ℤ × ℤ` supported in `Λ_M × Λ_M` sums to the finite sum over `Fin (2M+1)²`. -/
lemma tsum_eq_sum_fin (M : ℕ) (g : ℤ × ℤ → ℝ)
    (hg : ∀ p : ℤ × ℤ, ¬ (|p.1| ≤ M ∧ |p.2| ≤ M) → g p = 0) :
    ∑' p, g p = ∑ i : Fin (2 * M + 1), ∑ k : Fin (2 * M + 1),
      g ((i.val : ℤ) - M, (k.val : ℤ) - M) := by
  have hinj : Function.Injective (fun c : Fin (2 * M + 1) × Fin (2 * M + 1) =>
      (((c.1.val : ℤ) - M), ((c.2.val : ℤ) - M))) := by
    intro c d hcd
    simp only [Prod.mk.injEq] at hcd
    ext <;> omega
  rw [← hinj.tsum_eq, tsum_fintype, Fintype.sum_prod_type]
  intro p hp
  rw [Function.mem_support] at hp
  have hp' : |p.1| ≤ M ∧ |p.2| ≤ M := by
    by_contra hcon; exact hp (hg p hcon)
  refine ⟨(idxZ M p.1 hp'.1, idxZ M p.2 hp'.2), ?_⟩
  have h1 := abs_le.mp hp'.1
  have h2 := abs_le.mp hp'.2
  ext <;> simp [idxZ] <;> omega

/-- The finite-volume transported charge is the difference of the two crossing sums. -/
lemma transportM_eq (J h : ℝ) (n : CfgZ) (m : ℤ) (t : ℝ) (M : ℕ) :
    transportM J h n m t M =
      ∑' p, G1 n m (aM J h t M) p - ∑' p, G2 n m (aM J h t M) p := by
  have hsupp : ∀ p : ℤ × ℤ, ¬ (|p.1| ≤ M ∧ |p.2| ≤ M) → aM J h t M p.1 p.2 = 0 :=
    fun p hp => aM_eq_zero J h t M p.1 p.2 hp
  rw [tsum_eq_sum_fin M _ (fun p hp => by simp [G1, hsupp p hp]),
    tsum_eq_sum_fin M _ (fun p hp => by simp [G2, hsupp p hp])]
  unfold transportM
  rw [lgt_numOp_re, col_identity _ (sum_avgQ_normSq_col J h t)]
  congr 1
  · refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun k _ => ?_
    simp only [G1]
    rw [aM_eq_avgQ J h t M _ _ (abs_val_sub_le M i) (abs_val_sub_le M k), idxZ_val_sub,
      idxZ_val_sub]
    refine if_congr ?_ rfl rfl
    simp only [leftM, occM, mem_filter, mem_univ, true_and, not_le]
  · refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun k _ => ?_
    simp only [G2]
    rw [aM_eq_avgQ J h t M _ _ (abs_val_sub_le M i) (abs_val_sub_le M k), idxZ_val_sub,
      idxZ_val_sub]
    refine if_congr ?_ rfl rfl
    simp only [leftM, occM, mem_filter, mem_univ, true_and, not_le]

/-! ### Domination -/

/-- For a non-increasing non-negative profile, `F(d) ≤ √F(a) √F(b)` when `a, b ≤ d`. -/
lemma F_le_sqrt_mul {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) {a b d : ℕ} (ha : a ≤ d)
    (hb : b ≤ d) : F d ≤ √(F a) * √(F b) := by
  calc F d = √(F d) * √(F d) := (Real.mul_self_sqrt (hEC.nonneg d)).symm
    _ ≤ √(F a) * √(F b) := mul_le_mul (Real.sqrt_le_sqrt (hEC.antitone ha))
        (Real.sqrt_le_sqrt (hEC.antitone hb)) (Real.sqrt_nonneg _) (Real.sqrt_nonneg _)

/-- `x ↦ √F(|x - m|)` is summable on `ℤ` under (EC). -/
lemma summable_sqrt_F {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (m : ℤ) :
    Summable (fun x : ℤ => √(F (x - m).natAbs)) := by
  have hs : Summable (fun d : ℕ => √(F d)) :=
    hEC.summable_sqrt.of_nonneg_of_le (fun d => Real.sqrt_nonneg _)
      (fun d => le_mul_of_one_le_left (Real.sqrt_nonneg _)
        (one_le_pow₀ (by linarith [(Nat.cast_nonneg d : (0 : ℝ) ≤ d)])))
  have h0 : Summable (fun x : ℤ => √(F x.natAbs)) :=
    Summable.of_nat_of_neg (by simpa using hs) (by simpa using hs)
  refine (Equiv.addRight m).summable_iff.mp ?_
  simpa [Function.comp_def] using h0

/-- Domination of the first crossing sum, uniformly in `M`. -/
lemma norm_G1_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (n : CfgZ) (m : ℤ) (t : ℝ) (M : ℕ)
    (p : ℤ × ℤ) :
    ‖G1 n m (aM J h t M) p‖ ≤ √(F (p.1 - m).natAbs) * √(F (p.2 - m).natAbs) := by
  unfold G1
  split_ifs with hc
  · obtain ⟨h1, h2, -⟩ := hc
    rw [Real.norm_of_nonneg (aM_nonneg J h t M _ _)]
    exact (aM_le hEC t M _ _).trans (F_le_sqrt_mul hEC (by omega) (by omega))
  · rw [norm_zero]; positivity

/-- Domination of the second crossing sum, uniformly in `M`. -/
lemma norm_G2_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (n : CfgZ) (m : ℤ) (t : ℝ) (M : ℕ)
    (p : ℤ × ℤ) :
    ‖G2 n m (aM J h t M) p‖ ≤ √(F (p.1 - m).natAbs) * √(F (p.2 - m).natAbs) := by
  unfold G2
  split_ifs with hc
  · obtain ⟨h1, h2, -⟩ := hc
    rw [Real.norm_of_nonneg (aM_nonneg J h t M _ _)]
    exact (aM_le hEC t M _ _).trans (F_le_sqrt_mul hEC (by omega) (by omega))
  · rw [norm_zero]; positivity

end InfVolTransportAux

/-- **Transported charge in the thermodynamic limit** (paper: Theorem `thm:infvol`(b)): the limit
`T_m(t) = lim_{M → ∞} transportM` exists, and under (EC) `|T_m(t)| ≤ 2 ∑_d d F(d)` for every `t`. -/
theorem transport_tendsto (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (n : CfgZ) (m : ℤ) (t : ℝ) :
    ∃ T : ℝ, Tendsto (fun M : ℕ => transportM J h n m t M) atTop (𝓝 T) ∧
      |T| ≤ 2 * ∑' d : ℕ, (d : ℝ) * F d := by
  open InfVolTransportAux in
  have hbnd : Summable (fun p : ℤ × ℤ => √(F (p.1 - m).natAbs) * √(F (p.2 - m).natAbs)) :=
    (summable_sqrt_F hEC m).mul_of_nonneg (summable_sqrt_F hEC m)
      (fun _ => Real.sqrt_nonneg _) (fun _ => Real.sqrt_nonneg _)
  open InfVolTransportAux in
  have hT1 : Tendsto (fun M => ∑' p, G1 n m (aM J h t M) p) atTop
      (𝓝 (∑' p, G1 n m (aZ J h t) p)) := by
    refine tendsto_tsum_of_dominated_convergence hbnd (fun p => ?_)
      (Filter.Eventually.of_forall fun M p => norm_G1_le hEC n m t M p)
    by_cases hc : p.1 ≤ m ∧ m < p.2 ∧ n p.2 = 1
    · simp only [G1, ite_eq_left hc]
      exact tendsto_aM J h t _ _
    · simp only [G1, ite_eq_right hc]
      exact tendsto_const_nhds
  open InfVolTransportAux in
  have hT2 : Tendsto (fun M => ∑' p, G2 n m (aM J h t M) p) atTop
      (𝓝 (∑' p, G2 n m (aZ J h t) p)) := by
    refine tendsto_tsum_of_dominated_convergence hbnd (fun p => ?_)
      (Filter.Eventually.of_forall fun M p => norm_G2_le hEC n m t M p)
    by_cases hc : m < p.1 ∧ p.2 ≤ m ∧ n p.2 = 1
    · simp only [G2, ite_eq_left hc]
      exact tendsto_aM J h t _ _
    · simp only [G2, ite_eq_right hc]
      exact tendsto_const_nhds
  have hlim := (hT1.sub hT2).congr fun M => (InfVolTransportAux.transportM_eq J h n m t M).symm
  exact ⟨_, hlim, le_of_tendsto' ((continuous_abs.tendsto _).comp hlim)
    fun M => abs_transportM_le J h hEC n m t M⟩

end DFL
