import DFL.Sector
import DFL.FreeDyn
import DFL.Entropy
import DFL.Imbalance

/-!
# Main results: disorder-free localization in the 1D ℤ₂ LGT (paper: Section `sec:main`)

All statements are about the actual lattice gauge theory `H J h` with initial state
`|↑…↑⟩_z ⊗ |S⟩` (all links polarized along `z`, matter in a Fock state). The duality and the sector
averaging are proven in `DFL.Duality` and `DFL.Sector`; the only external input is
Hypothesis (EC) (`DFL.EC`), an explicit argument of the theorems that need it.

Main results:
* `DFL.lgt_density`: exact sector-averaged density formula;
* `DFL.lgt_transport`: bounded particle transport across any cut, uniformly in time and volume;
* `DFL.lgt_numberEntropy`: the number entropy of a half chain is bounded uniformly in time and
  volume;
* `DFL.lgtImbalance_eq`, `DFL.lgtImbalance_timeAvg`, `DFL.lgtImbalance_lower`: the
  density-wave imbalance has a strictly positive infinite-time value for every `h ≠ 0`, uniformly in
  the volume.
-/

open Matrix Finset Filter Topology
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Sites `x ≤ m` (the region to the left of the cut after site `m`). -/
def leftOf (m : ℕ) : Finset (Fin N) := univ.filter (fun x => x.val ≤ m)

/-- Density-wave (CDW) occupied set: even sites. -/
def cdw (N : ℕ) : Finset (Fin N) := univ.filter (fun x => Even x.val)

/-- Lift of a matter observable to the LGT (identity on the links). -/
abbrev matter (D : Matrix (Cfg N) (Cfg N) ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  cplx (D ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))

/-! ### Auxiliary lemmas: expectations, unitarity, averages -/

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

omit [DecidableEq ι] in
/-- `∑ |v_i|² = star v ⬝ᵥ v`. -/
lemma sum_normSq_eq_dotProduct (w : ι → ℂ) :
    ((∑ i, Complex.normSq (w i) : ℝ) : ℂ) = star w ⬝ᵥ w := by
  simp only [Complex.ofReal_sum, dotProduct, Pi.star_apply]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [Complex.normSq_eq_conj_mul_self]
  rfl

/-- An isometry preserves `∑ |v_i|²`. -/
lemma sum_normSq_mulVec (U : Matrix ι ι ℂ) (hU : Uᴴ * U = 1) (v : ι → ℂ) :
    ∑ i, Complex.normSq ((U *ᵥ v) i) = ∑ i, Complex.normSq (v i) := by
  have h1 : star (U *ᵥ v) ⬝ᵥ (U *ᵥ v) = star v ⬝ᵥ v := by
    rw [star_mulVec, dotProduct_mulVec, vecMul_vecMul, hU, vecMul_one]
  exact_mod_cast (sum_normSq_eq_dotProduct (U *ᵥ v)).trans
    (h1.trans (sum_normSq_eq_dotProduct v).symm)

/-- Rows of a unitary matrix are normalized. -/
lemma sum_normSq_row (U : Matrix ι ι ℂ) (hU : U * Uᴴ = 1) (x : ι) :
    ∑ z, Complex.normSq (U x z) = 1 := by
  have h := congrFun (congrFun hU x) x
  simp only [one_apply_eq] at h
  have : ((∑ z, Complex.normSq (U x z) : ℝ) : ℂ) = 1 := by
    rw [← h, Complex.ofReal_sum]
    refine Finset.sum_congr rfl (fun z _ => ?_)
    rw [← Complex.mul_conj]
    rfl
  exact_mod_cast this

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

/-! #### The Bernoulli average `avgQ` -/

/-- There are `2^N` charge configurations. -/
lemma card_cfg : Fintype.card (Cfg N) = 2 ^ N := by
  simp [Cfg, ZMod.card]

/-- Average of a constant. -/
lemma avgQ_const (c : ℝ) : avgQ (fun _ : Cfg N => c) = c := by
  simp only [avgQ, Finset.sum_const, Finset.card_univ, card_cfg, nsmul_eq_mul]
  push_cast
  field_simp

/-- Additivity of the average. -/
lemma avgQ_add (f g : Cfg N → ℝ) : avgQ (fun q => f q + g q) = avgQ f + avgQ g := by
  simp only [avgQ, Finset.sum_add_distrib]
  ring

/-- The average commutes with finite sums. -/
lemma avgQ_sum {κ : Type*} (s : Finset κ) (f : κ → Cfg N → ℝ) :
    avgQ (fun q => ∑ i ∈ s, f i q) = ∑ i ∈ s, avgQ (f i) := by
  simp only [avgQ]
  rw [Finset.sum_comm, Finset.mul_sum]

/-- Homogeneity of the average. -/
lemma avgQ_const_mul (c : ℝ) (f : Cfg N → ℝ) : avgQ (fun q => c * f q) = c * avgQ f := by
  simp only [avgQ, ← Finset.mul_sum]
  ring

/-- Average of a function minus a constant. -/
lemma avgQ_sub_const (f : Cfg N → ℝ) (c : ℝ) : avgQ (fun q => f q - c) = avgQ f - c := by
  have : avgQ (fun q => f q - c) = avgQ (fun q => f q + (-c)) := by
    simp only [sub_eq_add_neg]
  rw [this, avgQ_add, avgQ_const]
  ring

/-- Monotonicity of the average. -/
lemma avgQ_mono {f g : Cfg N → ℝ} (h : ∀ q, f q ≤ g q) : avgQ f ≤ avgQ g := by
  simp only [avgQ]
  exact mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun q _ => h q) (by positivity)

/-- Positivity of the average. -/
lemma avgQ_nonneg {f : Cfg N → ℝ} (h : ∀ q, 0 ≤ f q) : 0 ≤ avgQ f := by
  simp only [avgQ]
  exact mul_nonneg (by positivity) (Finset.sum_nonneg fun q _ => h q)

/-- `|avg f| ≤ avg |f|`. -/
lemma abs_avgQ_le (f : Cfg N → ℝ) : |avgQ f| ≤ avgQ (fun q => |f q|) := by
  simp only [avgQ]
  rw [abs_mul, abs_of_pos (by positivity)]
  exact mul_le_mul_of_nonneg_left (Finset.abs_sum_le_sum_abs _ _) (by positivity)

/-- Real part of a complex Bernoulli average of real numbers. -/
lemma re_avg (f : Cfg N → ℝ) :
    ((1 / 2 ^ N : ℂ) * ∑ q : Cfg N, ((f q : ℝ) : ℂ)).re = avgQ f := by
  simp only [avgQ]
  rw [← Complex.ofReal_sum]
  have : (1 / 2 ^ N : ℂ) = ((1 / 2 ^ N : ℝ) : ℂ) := by push_cast; ring
  rw [this, ← Complex.ofReal_mul, Complex.ofReal_re]

/-- **Cauchy–Schwarz** for the Bernoulli average: `avg(√u √v) ≤ √(avg u) √(avg v)`. -/
lemma avgQ_sqrt_mul_sqrt_le {u v : Cfg N → ℝ} (hu : ∀ q, 0 ≤ u q) (hv : ∀ q, 0 ≤ v q) :
    avgQ (fun q => √(u q) * √(v q)) ≤ √(avgQ u) * √(avgQ v) := by
  simp only [avgQ]
  have hc : (0 : ℝ) ≤ 1 / 2 ^ N := by positivity
  have h1 : (1 / 2 ^ N : ℝ) * ∑ q, √(u q) * √(v q) =
      ∑ q, √((1 / 2 ^ N : ℝ) * u q) * √((1 / 2 ^ N : ℝ) * v q) := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl (fun q _ => ?_)
    rw [Real.sqrt_mul hc, Real.sqrt_mul hc,
      show √(1 / 2 ^ N : ℝ) * √(u q) * (√(1 / 2 ^ N : ℝ) * √(v q)) =
        (√(1 / 2 ^ N : ℝ) * √(1 / 2 ^ N : ℝ)) * (√(u q) * √(v q)) by ring,
      Real.mul_self_sqrt hc]
  rw [h1, Finset.mul_sum, Finset.mul_sum]
  exact Real.sum_sqrt_mul_sqrt_le _ (fun q => mul_nonneg hc (hu q)) (fun q => mul_nonneg hc (hv q))

/-- **Minkowski-type bound** for the Bernoulli average:
`avg((∑ᵢ aᵢ)²) ≤ (∑ᵢ bᵢ)²` whenever `aᵢ ≥ 0` and `√(avg aᵢ²) ≤ bᵢ`. -/
lemma avgQ_sq_sum_le {κ : Type*} (s : Finset κ) (a : κ → Cfg N → ℝ) (b : κ → ℝ)
    (ha : ∀ i ∈ s, ∀ q, 0 ≤ a i q) (hb : ∀ i ∈ s, √(avgQ (fun q => a i q ^ 2)) ≤ b i) :
    avgQ (fun q => (∑ i ∈ s, a i q) ^ 2) ≤ (∑ i ∈ s, b i) ^ 2 := by
  have h1 : (fun q => (∑ i ∈ s, a i q) ^ 2) = fun q => ∑ i ∈ s, ∑ j ∈ s, a i q * a j q := by
    funext q; rw [sq, Finset.sum_mul_sum]
  rw [h1, avgQ_sum, sq, Finset.sum_mul_sum]
  refine Finset.sum_le_sum (fun i hi => ?_)
  rw [avgQ_sum]
  refine Finset.sum_le_sum (fun j hj => ?_)
  have e : (fun q => a i q * a j q) = fun q => √(a i q ^ 2) * √(a j q ^ 2) := by
    funext q
    rw [Real.sqrt_sq (ha i hi q), Real.sqrt_sq (ha j hj q)]
  rw [e]
  refine (avgQ_sqrt_mul_sqrt_le (fun q => sq_nonneg _) (fun q => sq_nonneg _)).trans ?_
  exact mul_le_mul (hb i hi) (hb j hj) (Real.sqrt_nonneg _)
    ((Real.sqrt_nonneg _).trans (hb i hi))

/-! #### Counting lemmas for sums of a distance profile -/

/-- Bounding `∑_{a ∈ s} g(φ a)` by the cardinalities of the fibres of `φ`. -/
lemma sum_comp_le_of_card_fiber {α : Type*} (s : Finset α) (φ : α → ℕ) (M : ℕ)
    (hφ : ∀ a ∈ s, φ a < M) (c : ℕ → ℝ)
    (hc : ∀ d, ((s.filter (fun a => φ a = d)).card : ℝ) ≤ c d)
    (g : ℕ → ℝ) (hg : ∀ d, 0 ≤ g d) :
    ∑ a ∈ s, g (φ a) ≤ ∑ d ∈ range M, c d * g d := by
  rw [← Finset.sum_fiberwise_of_maps_to' (t := range M) (fun a ha => mem_range.2 (hφ a ha))]
  refine Finset.sum_le_sum (fun d _ => ?_)
  rw [Finset.sum_const, nsmul_eq_mul]
  exact mul_le_mul_of_nonneg_right (hc d) (hg d)

/-- At most two sites are at a given distance from a site. -/
lemma card_dist_fiber_le (z : Fin N) (d : ℕ) :
    ((univ.filter (fun x : Fin N => Nat.dist x.val z.val = d)).card : ℝ) ≤ 2 := by
  have h : (univ.filter (fun x : Fin N => Nat.dist x.val z.val = d)).card ≤
      ({z.val - d, z.val + d} : Finset ℕ).card := by
    apply Finset.card_le_card_of_injOn Fin.val
    · intro x hx
      simp only [coe_filter, mem_univ, true_and, Set.mem_ofPred_eq] at hx
      simp only [coe_insert, coe_singleton, Set.mem_insert_iff, Set.mem_singleton_iff]
      unfold Nat.dist at hx
      omega
    · intro x _ y _ hxy
      exact Fin.ext hxy
  have h2 : ({z.val - d, z.val + d} : Finset ℕ).card ≤ 2 := Finset.card_le_two
  exact_mod_cast h.trans h2

/-- At most `d` pairs `(x, y)` with `x ≤ m < y` are at distance `d`. -/
lemma card_cross_fiber_le (m d : ℕ) :
    (((leftOf (N := N) m ×ˢ (univ \ leftOf m)).filter
      (fun p : Fin N × Fin N => Nat.dist p.1.val p.2.val = d)).card : ℝ) ≤ d := by
  have h : ((leftOf (N := N) m ×ˢ (univ \ leftOf m)).filter
      (fun p : Fin N × Fin N => Nat.dist p.1.val p.2.val = d)).card ≤ (range d).card := by
    apply Finset.card_le_card_of_injOn (fun p : Fin N × Fin N => m - p.1.val)
    · intro p hp
      simp only [leftOf, coe_filter, mem_product, mem_filter, mem_univ, true_and, Finset.mem_sdiff,
        Set.mem_ofPred_eq] at hp
      simp only [coe_range, Set.mem_Iio]
      unfold Nat.dist at hp
      omega
    · intro p hp p' hp' hpp'
      simp only [leftOf, coe_filter, mem_product, mem_filter, mem_univ, true_and, Finset.mem_sdiff,
        Set.mem_ofPred_eq] at hp hp'
      simp only at hpp'
      unfold Nat.dist at hp hp'
      ext <;> omega
  rw [card_range] at h
  exact_mod_cast h

/-- `∑_x g(|x - z|) ≤ 2 ∑_d g(d)`. -/
lemma sum_dist_le (g : ℕ → ℝ) (hg : ∀ d, 0 ≤ g d) (hs : Summable g) (z : Fin N) :
    ∑ x : Fin N, g (Nat.dist x.val z.val) ≤ 2 * ∑' d, g d := by
  calc ∑ x : Fin N, g (Nat.dist x.val z.val) ≤ ∑ d ∈ range N, 2 * g d :=
        sum_comp_le_of_card_fiber univ (fun x : Fin N => Nat.dist x.val z.val) N
          (fun x _ => by have := x.isLt; have := z.isLt; unfold Nat.dist; omega) (fun _ => 2)
          (card_dist_fiber_le z) g hg
    _ = 2 * ∑ d ∈ range N, g d := by rw [Finset.mul_sum]
    _ ≤ 2 * ∑' d, g d := by
        gcongr
        exact hs.sum_le_tsum _ (fun d _ => hg d)

/-- **Crossing pairs**: `∑_{x ≤ m < y} g(|x - y|) ≤ ∑_d d g(d)`. -/
lemma sum_cross_le (g : ℕ → ℝ) (hg : ∀ d, 0 ≤ g d) (hs : Summable (fun d : ℕ => (d : ℝ) * g d))
    (m : ℕ) :
    ∑ x ∈ leftOf (N := N) m, ∑ y ∈ univ \ leftOf (N := N) m, g (Nat.dist x.val y.val) ≤
      ∑' d : ℕ, (d : ℝ) * g d := by
  rw [← Finset.sum_product' (f := fun x y : Fin N => g (Nat.dist x.val y.val))]
  calc _ ≤ ∑ d ∈ range N, (d : ℝ) * g d :=
        sum_comp_le_of_card_fiber _ (fun p : Fin N × Fin N => Nat.dist p.1.val p.2.val) N
          (fun p _ => by have := p.1.isLt; have := p.2.isLt; unfold Nat.dist; omega)
          (fun d => d) (card_cross_fiber_le m) g hg
    _ ≤ _ := hs.sum_le_tsum _ (fun d _ => mul_nonneg (Nat.cast_nonneg d) (hg d))

/-- **Exact density formula**: `⟨n_x(t)⟩ = 2^{-N} ∑_q ∑_{z ∈ S} |(e^{-ith_q})_{xz}|²`. -/
theorem lgt_density (J h : ℝ) (S : Finset (Fin N)) (x : Fin N) (t : ℝ) :
    (expect (matter (diagonal fun n => occ (n x))) (lgtState J h (fockVec S) t)).re =
      avgQ (fun q => ∑ z ∈ S, Complex.normSq (evo (hq J h q) t x z)) := by
  sorry

/-- **Bounded transport** (paper: Theorem `thm:transport`): under Hypothesis (EC), the number of
particles transferred across any cut is bounded uniformly in time and volume. -/
theorem lgt_transport (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ)
    (t : ℝ) :
    |(expect (matter (numOp (leftOf m))) (lgtState J h (fockVec S) t)).re -
        ((S ∩ leftOf m).card : ℝ)| ≤
      2 * ∑' d : ℕ, (d : ℝ) * F d := by
  sorry

/-- Distribution of the particle number `N_R` in the LGT state at time `t`. -/
def lgtNumDist (J h : ℝ) (R : Finset (Fin N)) (ψ : Cfg N → ℂ) (t : ℝ) (k : ℕ) : ℝ :=
  (expect (matter (numProj R k)) (lgtState J h ψ t)).re

/-- The number distribution is a probability distribution on `{0, …, N}`. -/
theorem lgtNumDist_isProb (J h : ℝ) (R S : Finset (Fin N)) (t : ℝ) :
    IsProb N (lgtNumDist J h R (fockVec S) t) := by
  sorry

/-- **Bounded number entropy** (paper: Theorem `thm:numberentropy`): under Hypothesis (EC), the
Shannon entropy of the half-chain particle-number distribution is bounded uniformly in time,
volume, position of the cut and occupied set. -/
theorem lgt_numberEntropy (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ K : ℝ, ∀ (N : ℕ) (m : ℕ) (S : Finset (Fin N)) (t : ℝ),
      shannon N (lgtNumDist J h (leftOf m) (fockVec S) t) ≤
        (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (K + 1 / 12)) := by
  sorry

/-- The density-wave imbalance `I(t) = (2/N) ∑_x (-1)^x ⟨n_x(t)⟩` of the LGT. -/
def lgtImbalance (J h : ℝ) (N : ℕ) (t : ℝ) : ℝ :=
  (2 / N : ℝ) * ∑ x : Fin N, (-1 : ℝ) ^ x.val *
    (expect (matter (diagonal fun n => occ (n x))) (lgtState J h (fockVec (cdw N)) t)).re

/-- The LGT imbalance is the Bernoulli average of the one-body imbalances (`N` even). -/
theorem lgtImbalance_eq (J h : ℝ) (hN : Even N) (t : ℝ) :
    lgtImbalance J h N t = avgQ (fun q : Cfg N => imbOne J h q t) := by
  sorry

/-- The time-averaged LGT imbalance converges to `Ī = 2^{-N} ∑_q N⁻¹ ∑_α a_α(q)²`. -/
theorem lgtImbalance_timeAvg (J h : ℝ) (hJ : J ≠ 0) (hN : Even N) :
    Tendsto (fun T : ℝ => (1 / T) * ∫ t in (0 : ℝ)..T, lgtImbalance J h N t) atTop
      (𝓝 (avgQ (fun q : Cfg N => imbInf J h q))) := by
  sorry

/-- **Persistent imbalance** (paper: Theorem `thm:imbalance`): for every `h ≠ 0`, under
Hypothesis (EC), the infinite-time imbalance of the LGT is bounded below by a positive constant
independent of the (even) volume. -/
theorem lgtImbalance_lower (J h : ℝ) (hJ : J ≠ 0) (hh : h ≠ 0) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ c : ℝ, 0 < c ∧ ∀ (N : ℕ), Even N → 2 ≤ N →
      ∀ ε > 0, ∀ᶠ T in atTop,
        c - ε ≤ (1 / T) * ∫ t in (0 : ℝ)..T, lgtImbalance J h N t := by
  sorry

end DFL
