import DFL.AreaLawSharp

/-!
# Annealed one-body bounds for cuts and intervals (paper: Lemma `lem:annealed`(ii),(iii))

Under (EC), with `Σ₁(F) = ∑_{d ≥ 1} d F(d)`, `Σ'₁(F) = ∑_{d ≥ 1} d F(d)^{1/2}` and `Θ(F)` (`thetaF`):

* (ii) for the cut after the site `m`: `E ∑_{x ≤ m < z} Q_q(x, z) ≤ Σ₁(F)` (`annealed_cut_Q`) and
  `E (∑_{x ≤ m < z} Q_q(x, z))² ≤ Σ'₁(F)²` (`annealed_cut_Q_sq`); the first bound of (ii) is
  `annealed_cut`;
* (iii) for an interval `B = [a, b]`, the same bounds for the sums over `x ∈ B`, `y ∉ B`, with the
  constants `2Θ(F)`, `2Σ₁(F)` and `4Σ'₁(F)²` (`annealed_ivl_K`, `annealed_ivl_Q`, `annealed_ivl_Q_sq`).
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The interval `[a, b]` of sites. -/
def ivl (a b : ℕ) : Finset (Fin N) := univ.filter (fun x => a ≤ x.val ∧ x.val ≤ b)

namespace AnnealedIntervalAux

open AreaLawAux AreaLawSharpAux

/-! ### Elementary facts -/

/-- `leftOf m = segment (m + 1)`. -/
lemma leftOf_eq_segment_succ (m : ℕ) : leftOf (N := N) m = segment (m + 1) := by
  ext x
  simp [leftOf, segment]

/-- Additivity of the Bernoulli average. -/
lemma avgQ_add_eq (f g : Cfg N → ℝ) : avgQ (fun q => f q + g q) = avgQ f + avgQ g := by
  unfold avgQ; rw [Finset.sum_add_distrib]; ring

/-- `Q(x, z)² ≤ Q(x, z)`, since `0 ≤ Q ≤ 1`. -/
lemma corrQ_sq_le (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x z : Fin N) :
    corrQ A hA x z ^ 2 ≤ corrQ A hA x z := by
  have h0 := corrQ_nonneg A hA x z
  have h1 := ECProofAux.corrQ_le_one A hA x z
  nlinarith

/-- `y² ≤ 2u² + 2v²` whenever `0 ≤ y ≤ u + v`. -/
lemma sq_le_two_sq_add_two_sq {y u v : ℝ} (h0 : 0 ≤ y) (h : y ≤ u + v) :
    y ^ 2 ≤ 2 * u ^ 2 + 2 * v ^ 2 := by
  nlinarith [mul_nonneg h0 (sub_nonneg.2 h), sq_nonneg (u - v), sq_nonneg (u + v - y)]

/-! ### Summability under (EC) -/

/-- Under (EC), `d ↦ d F(d)` is summable. -/
lemma summable_mul_F {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) :
    Summable (fun d : ℕ => (d : ℝ) * F d) :=
  hEC.summable.of_nonneg_of_le (fun d => mul_nonneg (Nat.cast_nonneg d) (hEC.nonneg d))
    (fun d => mul_le_mul_of_nonneg_right
      (by nlinarith [(Nat.cast_nonneg d : (0 : ℝ) ≤ d)]) (hEC.nonneg d))

/-- Under (EC), `d ↦ d √F(d)` is summable. -/
lemma summable_mul_sqrt_F {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) :
    Summable (fun d : ℕ => (d : ℝ) * √(F d)) :=
  hEC.summable_sqrt.of_nonneg_of_le
    (fun d => mul_nonneg (Nat.cast_nonneg d) (Real.sqrt_nonneg _))
    (fun d => mul_le_mul_of_nonneg_right
      (by nlinarith [(Nat.cast_nonneg d : (0 : ℝ) ≤ d)]) (Real.sqrt_nonneg _))

/-! ### Counting the pairs across a cut -/

/-- Bounding `∑_{a ∈ s} g(φ a)` by the cardinalities of the fibres of `φ`. -/
lemma sum_le_of_card_fiber {α : Type*} (s : Finset α) (φ : α → ℕ) (M : ℕ)
    (hφ : ∀ a ∈ s, φ a < M) (c : ℕ → ℝ)
    (hc : ∀ d, ((s.filter (fun a => φ a = d)).card : ℝ) ≤ c d)
    (g : ℕ → ℝ) (hg : ∀ d, 0 ≤ g d) :
    ∑ a ∈ s, g (φ a) ≤ ∑ d ∈ range M, c d * g d := by
  rw [← Finset.sum_fiberwise_of_maps_to' (t := range M) (fun a ha => mem_range.2 (hφ a ha))]
  refine Finset.sum_le_sum (fun d _ => ?_)
  rw [Finset.sum_const, nsmul_eq_mul]
  exact mul_le_mul_of_nonneg_right (hc d) (hg d)

/-- At most `d` pairs `x < ℓ ≤ y` are at distance `d`. -/
lemma card_seg_fiber_le (ℓ d : ℕ) :
    (((segment (N := N) ℓ ×ˢ (univ \ segment ℓ)).filter
      (fun p : Fin N × Fin N => Nat.dist p.1.val p.2.val = d)).card : ℝ) ≤ d := by
  have h : ((segment (N := N) ℓ ×ˢ (univ \ segment ℓ)).filter
      (fun p : Fin N × Fin N => Nat.dist p.1.val p.2.val = d)).card ≤ (range d).card := by
    apply Finset.card_le_card_of_injOn (fun p : Fin N × Fin N => ℓ - 1 - p.1.val)
    · intro p hp
      simp only [segment, coe_filter, mem_product, mem_filter, mem_univ, true_and,
        Finset.mem_sdiff, Set.mem_ofPred_eq] at hp
      simp only [coe_range, Set.mem_Iio]
      unfold Nat.dist at hp
      omega
    · intro p hp p' hp' hpp'
      simp only [segment, coe_filter, mem_product, mem_filter, mem_univ, true_and,
        Finset.mem_sdiff, Set.mem_ofPred_eq] at hp hp'
      simp only at hpp'
      unfold Nat.dist at hp hp'
      ext <;> omega
  rw [card_range] at h
  exact_mod_cast h

/-- **Crossing pairs of a cut**: `∑_{x < ℓ ≤ y} g(|x - y|) ≤ ∑_d d g(d)`. -/
lemma sum_seg_cross_le (g : ℕ → ℝ) (hg : ∀ d, 0 ≤ g d)
    (hs : Summable (fun d : ℕ => (d : ℝ) * g d)) (ℓ : ℕ) :
    ∑ x ∈ segment (N := N) ℓ, ∑ y ∈ univ \ segment (N := N) ℓ, g (Nat.dist x.val y.val) ≤
      ∑' d : ℕ, (d : ℝ) * g d := by
  rw [← Finset.sum_product' (f := fun x y : Fin N => g (Nat.dist x.val y.val))]
  calc _ ≤ ∑ d ∈ range N, (d : ℝ) * g d :=
        sum_le_of_card_fiber _ (fun p : Fin N × Fin N => Nat.dist p.1.val p.2.val) N
          (fun p _ => by have := p.1.isLt; have := p.2.isLt; unfold Nat.dist; omega)
          (fun d => d) (card_seg_fiber_le ℓ) g hg
    _ ≤ _ := hs.sum_le_tsum _ (fun d _ => mul_nonneg (Nat.cast_nonneg d) (hg d))

/-! ### Annealed bounds for a cut -/

/-- `E ∑_{x < ℓ ≤ z} Q(x, z) ≤ Σ₁(F)`. -/
lemma avgQ_seg_corrQ_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (ℓ : ℕ) :
    avgQ (fun q : Cfg N => ∑ x ∈ segment ℓ, ∑ z ∈ univ \ segment ℓ,
      corrQ (hq J h q) (hq_isHermitian J h q) x z) ≤ ∑' d : ℕ, (d : ℝ) * F d := by
  calc _ = ∑ x ∈ segment (N := N) ℓ, ∑ z ∈ univ \ segment (N := N) ℓ,
        avgQ (fun q : Cfg N => corrQ (hq J h q) (hq_isHermitian J h q) x z) := by
        rw [avgQ_sum']
        exact Finset.sum_congr rfl fun x _ => avgQ_sum' _ _
    _ ≤ ∑ x ∈ segment (N := N) ℓ, ∑ z ∈ univ \ segment (N := N) ℓ, F (Nat.dist x.val z.val) :=
        Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun z _ => hEC.bound N x z
    _ ≤ _ := sum_seg_cross_le F hEC.nonneg (summable_mul_F hEC) ℓ

/-- **Minkowski-type bound** for the Bernoulli average:
`avg((∑ᵢ aᵢ)²) ≤ (∑ᵢ bᵢ)²` whenever `aᵢ ≥ 0` and `√(avg aᵢ²) ≤ bᵢ`. -/
lemma avgQ_sq_sum_le' {κ : Type*} (s : Finset κ) (a : κ → Cfg N → ℝ) (b : κ → ℝ)
    (ha : ∀ i ∈ s, ∀ q, 0 ≤ a i q) (hb : ∀ i ∈ s, √(avgQ (fun q => a i q ^ 2)) ≤ b i) :
    avgQ (fun q => (∑ i ∈ s, a i q) ^ 2) ≤ (∑ i ∈ s, b i) ^ 2 := by
  have h1 : (fun q => (∑ i ∈ s, a i q) ^ 2) = fun q => ∑ i ∈ s, ∑ j ∈ s, a i q * a j q := by
    funext q; rw [sq, Finset.sum_mul_sum]
  rw [h1, avgQ_sum', sq, Finset.sum_mul_sum]
  refine Finset.sum_le_sum (fun i hi => ?_)
  rw [avgQ_sum']
  refine Finset.sum_le_sum (fun j hj => ?_)
  have e : (fun q => a i q * a j q) = fun q => √(a i q ^ 2) * √(a j q ^ 2) := by
    funext q
    rw [Real.sqrt_sq (ha i hi q), Real.sqrt_sq (ha j hj q)]
  rw [e]
  refine (avgQ_sqrt_mul_sqrt_le _ _ (fun q => sq_nonneg _) (fun q => sq_nonneg _)).trans ?_
  exact mul_le_mul (hb i hi) (hb j hj) (Real.sqrt_nonneg _)
    ((Real.sqrt_nonneg _).trans (hb i hi))

/-- `E (∑_{x < ℓ ≤ z} Q(x, z))² ≤ Σ'₁(F)²`. -/
lemma avgQ_seg_corrQ_sq_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (ℓ : ℕ) :
    avgQ (fun q : Cfg N => (∑ x ∈ segment ℓ, ∑ z ∈ univ \ segment ℓ,
      corrQ (hq J h q) (hq_isHermitian J h q) x z) ^ 2) ≤
      (∑' d : ℕ, (d : ℝ) * √(F d)) ^ 2 := by
  have e : ∀ q : Cfg N, ∑ x ∈ segment (N := N) ℓ, ∑ z ∈ univ \ segment (N := N) ℓ,
      corrQ (hq J h q) (hq_isHermitian J h q) x z =
      ∑ p ∈ segment (N := N) ℓ ×ˢ (univ \ segment (N := N) ℓ),
        corrQ (hq J h q) (hq_isHermitian J h q) p.1 p.2 :=
    fun q => (Finset.sum_product' _ _ _).symm
  simp_rw [e]
  have hpair : ∀ x z : Fin N,
      √(avgQ (fun q : Cfg N => corrQ (hq J h q) (hq_isHermitian J h q) x z ^ 2)) ≤
        √(F (Nat.dist x.val z.val)) := fun x z =>
    Real.sqrt_le_sqrt ((avgQ_mono' fun q => corrQ_sq_le _ _ x z).trans (hEC.bound N x z))
  refine (avgQ_sq_sum_le' (segment (N := N) ℓ ×ˢ (univ \ segment (N := N) ℓ))
    (fun p q => corrQ (hq J h q) (hq_isHermitian J h q) p.1 p.2)
    (fun p => √(F (Nat.dist p.1.val p.2.val)))
    (fun p _ q => corrQ_nonneg _ _ _ _) (fun p _ => hpair p.1 p.2)).trans ?_
  refine pow_le_pow_left₀ (Finset.sum_nonneg fun p _ => Real.sqrt_nonneg _) ?_ 2
  rw [Finset.sum_product]
  exact sum_seg_cross_le (fun d => √(F d)) (fun d => Real.sqrt_nonneg _)
    (summable_mul_sqrt_F hEC) ℓ

/-! ### Splitting an interval into two cuts -/

/-- Sites outside `[a, b]` lie to the left of `a` or to the right of `b`. -/
lemma compl_ivl_subset (a b : ℕ) :
    univ \ ivl (N := N) a b ⊆ segment a ∪ (univ \ segment (b + 1)) := by
  intro y hy
  simp only [ivl, segment, Finset.mem_sdiff, Finset.mem_univ, Finset.mem_filter, true_and,
    Finset.mem_union, not_and, not_le, not_lt] at hy ⊢
  by_cases hya : y.val < a
  · exact Or.inl hya
  · exact Or.inr (by have := hy (by omega); omega)

/-- The interval `[a, b]` lies to the right of the cut at `a`. -/
lemma ivl_subset_compl_segment (a b : ℕ) : ivl (N := N) a b ⊆ univ \ segment a := by
  intro x hx
  simp only [ivl, segment, Finset.mem_sdiff, Finset.mem_univ, Finset.mem_filter, true_and,
    not_lt] at hx ⊢
  exact hx.1

/-- The interval `[a, b]` lies to the left of the cut at `b + 1`. -/
lemma ivl_subset_segment (a b : ℕ) : ivl (N := N) a b ⊆ segment (b + 1) := by
  intro x hx
  simp only [ivl, segment, Finset.mem_univ, Finset.mem_filter, true_and] at hx ⊢
  omega

/-- For a nonnegative symmetric kernel `G`, the sum over `x ∈ [a, b]`, `y ∉ [a, b]` is at most the
sum of the cut sums at `a` (pairs `y < a ≤ x`) and at `b + 1` (pairs `x ≤ b < y`). -/
lemma sum_ivl_le_two_cuts (G : Fin N → Fin N → ℝ) (hG0 : ∀ x y, 0 ≤ G x y)
    (hGs : ∀ x y, G x y = G y x) (a b : ℕ) :
    ∑ x ∈ ivl a b, ∑ y ∈ univ \ ivl a b, G x y ≤
      ∑ x ∈ segment (N := N) a, ∑ y ∈ univ \ segment a, G x y +
        ∑ x ∈ segment (N := N) (b + 1), ∑ y ∈ univ \ segment (b + 1), G x y := by
  have h1 : ∀ x, ∑ y ∈ univ \ ivl (N := N) a b, G x y ≤
      ∑ y ∈ segment (N := N) a, G x y + ∑ y ∈ univ \ segment (N := N) (b + 1), G x y := by
    intro x
    calc ∑ y ∈ univ \ ivl (N := N) a b, G x y
        ≤ ∑ y ∈ segment (N := N) a ∪ (univ \ segment (b + 1)), G x y :=
          Finset.sum_le_sum_of_subset_of_nonneg (compl_ivl_subset a b) (fun y _ _ => hG0 x y)
      _ ≤ ∑ y ∈ segment (N := N) a, G x y + ∑ y ∈ univ \ segment (N := N) (b + 1), G x y := by
          rw [← Finset.sum_union_inter]
          exact le_add_of_nonneg_right (Finset.sum_nonneg fun y _ => hG0 x y)
  calc ∑ x ∈ ivl a b, ∑ y ∈ univ \ ivl a b, G x y
      ≤ ∑ x ∈ ivl (N := N) a b, (∑ y ∈ segment (N := N) a, G x y +
          ∑ y ∈ univ \ segment (N := N) (b + 1), G x y) := Finset.sum_le_sum fun x _ => h1 x
    _ = ∑ x ∈ ivl (N := N) a b, ∑ y ∈ segment (N := N) a, G x y +
          ∑ x ∈ ivl (N := N) a b, ∑ y ∈ univ \ segment (N := N) (b + 1), G x y :=
        Finset.sum_add_distrib
    _ ≤ ∑ x ∈ segment (N := N) a, ∑ y ∈ univ \ segment a, G x y +
        ∑ x ∈ segment (N := N) (b + 1), ∑ y ∈ univ \ segment (b + 1), G x y := by
        refine add_le_add ?_ ?_
        · calc ∑ x ∈ ivl (N := N) a b, ∑ y ∈ segment (N := N) a, G x y
              ≤ ∑ x ∈ univ \ segment (N := N) a, ∑ y ∈ segment (N := N) a, G x y :=
                Finset.sum_le_sum_of_subset_of_nonneg (ivl_subset_compl_segment a b)
                  (fun x _ _ => Finset.sum_nonneg fun y _ => hG0 x y)
            _ = ∑ y ∈ segment (N := N) a, ∑ x ∈ univ \ segment (N := N) a, G x y :=
                Finset.sum_comm
            _ = ∑ x ∈ segment (N := N) a, ∑ y ∈ univ \ segment a, G x y :=
                Finset.sum_congr rfl fun y _ => Finset.sum_congr rfl fun x _ => hGs x y
        · exact Finset.sum_le_sum_of_subset_of_nonneg (ivl_subset_segment a b)
            (fun x _ _ => Finset.sum_nonneg fun y _ => hG0 x y)

end AnnealedIntervalAux

open AreaLawAux AreaLawSharpAux AnnealedIntervalAux

/-- `E ∑_{x ≤ m < z} Q_q(x, z) ≤ Σ₁(F)` (paper: Lemma `lem:annealed`(ii)). -/
theorem annealed_cut_Q (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (m : ℕ) :
    avgQ (fun q : Cfg N => ∑ x ∈ leftOf m, ∑ z ∈ univ \ leftOf m,
      corrQ (hq J h q) (hq_isHermitian J h q) x z) ≤ ∑' d : ℕ, (d : ℝ) * F d := by
  rw [leftOf_eq_segment_succ]
  exact avgQ_seg_corrQ_le hEC (m + 1)

/-- `E (∑_{x ≤ m < z} Q_q(x, z))² ≤ Σ'₁(F)²` (paper: Lemma `lem:annealed`(ii)). -/
theorem annealed_cut_Q_sq (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (m : ℕ) :
    avgQ (fun q : Cfg N => (∑ x ∈ leftOf m, ∑ z ∈ univ \ leftOf m,
      corrQ (hq J h q) (hq_isHermitian J h q) x z) ^ 2) ≤ (∑' d : ℕ, (d : ℝ) * √(F d)) ^ 2 := by
  rw [leftOf_eq_segment_succ]
  exact avgQ_seg_corrQ_sq_le hEC (m + 1)

/-- `E sup_t ∑_{x ∈ B, y ∉ B} |K_q(t)_{xy}| ≤ 2Θ(F)` for an interval `B` (paper: Lemma `lem:annealed`(iii)). -/
theorem annealed_ivl_K (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (a b : ℕ) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ,
      ∑ x ∈ ivl a b, ∑ y ∈ univ \ ivl a b, ‖corrMat (hq J h q) S t x y‖) ≤ 2 * thetaF F := by
  have hpt : ∀ q : Cfg N, (⨆ t : ℝ,
      ∑ x ∈ ivl a b, ∑ y ∈ univ \ ivl a b, ‖corrMat (hq J h q) S t x y‖) ≤
      (∑ x ∈ segment (N := N) a, ∑ y ∈ univ \ segment a,
        ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) x z *
          corrQ (hq J h q) (hq_isHermitian J h q) y z) +
      ∑ x ∈ segment (N := N) (b + 1), ∑ y ∈ univ \ segment (b + 1),
        ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) x z *
          corrQ (hq J h q) (hq_isHermitian J h q) y z := by
    intro q
    refine ciSup_le fun t => ?_
    refine (Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ =>
      norm_corrMat_le_corrQ J h q S t x y).trans ?_
    exact sum_ivl_le_two_cuts
      (fun x y => ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) x z *
        corrQ (hq J h q) (hq_isHermitian J h q) y z)
      (fun x y => Finset.sum_nonneg fun z _ =>
        mul_nonneg (corrQ_nonneg _ _ x z) (corrQ_nonneg _ _ y z))
      (fun x y => Finset.sum_congr rfl fun z _ => mul_comm _ _) a b
  refine (avgQ_mono' hpt).trans ?_
  rw [avgQ_add_eq, two_mul]
  exact add_le_add (avgQ_cut_corrQ_le hEC a) (avgQ_cut_corrQ_le hEC (b + 1))

/-- `E ∑_{x ∈ B, z ∉ B} Q_q(x, z) ≤ 2Σ₁(F)` (paper: Lemma `lem:annealed`(iii)). -/
theorem annealed_ivl_Q (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (a b : ℕ) :
    avgQ (fun q : Cfg N => ∑ x ∈ ivl a b, ∑ z ∈ univ \ ivl a b,
      corrQ (hq J h q) (hq_isHermitian J h q) x z) ≤ 2 * ∑' d : ℕ, (d : ℝ) * F d := by
  have hpt : ∀ q : Cfg N, ∑ x ∈ ivl a b, ∑ z ∈ univ \ ivl a b,
      corrQ (hq J h q) (hq_isHermitian J h q) x z ≤
      (∑ x ∈ segment (N := N) a, ∑ z ∈ univ \ segment a,
        corrQ (hq J h q) (hq_isHermitian J h q) x z) +
      ∑ x ∈ segment (N := N) (b + 1), ∑ z ∈ univ \ segment (b + 1),
        corrQ (hq J h q) (hq_isHermitian J h q) x z := fun q =>
    sum_ivl_le_two_cuts _ (fun x z => corrQ_nonneg _ _ x z)
      (fun x z => ECProofAux.corrQ_symm _ _ x z) a b
  refine (avgQ_mono' hpt).trans ?_
  rw [avgQ_add_eq, two_mul]
  exact add_le_add (avgQ_seg_corrQ_le hEC a) (avgQ_seg_corrQ_le hEC (b + 1))

/-- `E (∑_{x ∈ B, z ∉ B} Q_q(x, z))² ≤ 4Σ'₁(F)²` (paper: Lemma `lem:annealed`(iii)). -/
theorem annealed_ivl_Q_sq (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (a b : ℕ) :
    avgQ (fun q : Cfg N => (∑ x ∈ ivl a b, ∑ z ∈ univ \ ivl a b,
      corrQ (hq J h q) (hq_isHermitian J h q) x z) ^ 2) ≤
      4 * (∑' d : ℕ, (d : ℝ) * √(F d)) ^ 2 := by
  have hpt : ∀ q : Cfg N, (∑ x ∈ ivl a b, ∑ z ∈ univ \ ivl a b,
      corrQ (hq J h q) (hq_isHermitian J h q) x z) ^ 2 ≤
      2 * (∑ x ∈ segment (N := N) a, ∑ z ∈ univ \ segment a,
        corrQ (hq J h q) (hq_isHermitian J h q) x z) ^ 2 +
      2 * (∑ x ∈ segment (N := N) (b + 1), ∑ z ∈ univ \ segment (b + 1),
        corrQ (hq J h q) (hq_isHermitian J h q) x z) ^ 2 := fun q =>
    sq_le_two_sq_add_two_sq
      (Finset.sum_nonneg fun x _ => Finset.sum_nonneg fun z _ => corrQ_nonneg _ _ x z)
      (sum_ivl_le_two_cuts _ (fun x z => corrQ_nonneg _ _ x z)
        (fun x z => ECProofAux.corrQ_symm _ _ x z) a b)
  refine (avgQ_mono' hpt).trans ?_
  rw [avgQ_add_eq, avgQ_const_mul', avgQ_const_mul']
  have e1 := avgQ_seg_corrQ_sq_le (N := N) hEC a
  have e2 := avgQ_seg_corrQ_sq_le (N := N) hEC (b + 1)
  linarith

end DFL
