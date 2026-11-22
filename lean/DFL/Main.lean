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

namespace MainAux

/-! ### Auxiliary lemmas: expectations, unitarity, averages -/

section Aux

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Expectation of a complexified real diagonal matrix: `⟨v, diag(d) v⟩ = ∑ d_i |v_i|²`. -/
private lemma expect_cplx_diagonal (d : ι → ℝ) (v : ι → ℂ) :
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
private lemma sum_normSq_eq_dotProduct (w : ι → ℂ) :
    ((∑ i, Complex.normSq (w i) : ℝ) : ℂ) = star w ⬝ᵥ w := by
  simp only [Complex.ofReal_sum, dotProduct, Pi.star_apply]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [Complex.normSq_eq_conj_mul_self]
  rfl

/-- An isometry preserves `∑ |v_i|²`. -/
private lemma sum_normSq_mulVec (U : Matrix ι ι ℂ) (hU : Uᴴ * U = 1) (v : ι → ℂ) :
    ∑ i, Complex.normSq ((U *ᵥ v) i) = ∑ i, Complex.normSq (v i) := by
  have h1 : star (U *ᵥ v) ⬝ᵥ (U *ᵥ v) = star v ⬝ᵥ v := by
    rw [star_mulVec, dotProduct_mulVec, vecMul_vecMul, hU, vecMul_one]
  exact_mod_cast (sum_normSq_eq_dotProduct (U *ᵥ v)).trans
    (h1.trans (sum_normSq_eq_dotProduct v).symm)

/-- Rows of a unitary matrix are normalized. -/
private lemma sum_normSq_row (U : Matrix ι ι ℂ) (hU : U * Uᴴ = 1) (x : ι) :
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
private lemma sum_normSq_col (U : Matrix ι ι ℂ) (hU : Uᴴ * U = 1) (z : ι) :
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
private lemma card_cfg : Fintype.card (Cfg N) = 2 ^ N := by
  simp [Cfg, ZMod.card]

/-- Average of a constant. -/
private lemma avgQ_const (c : ℝ) : avgQ (fun _ : Cfg N => c) = c := by
  simp only [avgQ, Finset.sum_const, Finset.card_univ, card_cfg, nsmul_eq_mul]
  push_cast
  field_simp

/-- Additivity of the average. -/
private lemma avgQ_add (f g : Cfg N → ℝ) : avgQ (fun q => f q + g q) = avgQ f + avgQ g := by
  simp only [avgQ, Finset.sum_add_distrib]
  ring

/-- The average commutes with finite sums. -/
private lemma avgQ_sum {κ : Type*} (s : Finset κ) (f : κ → Cfg N → ℝ) :
    avgQ (fun q => ∑ i ∈ s, f i q) = ∑ i ∈ s, avgQ (f i) := by
  simp only [avgQ]
  rw [Finset.sum_comm, Finset.mul_sum]

/-- Homogeneity of the average. -/
private lemma avgQ_const_mul (c : ℝ) (f : Cfg N → ℝ) : avgQ (fun q => c * f q) = c * avgQ f := by
  simp only [avgQ, ← Finset.mul_sum]
  ring

/-- Average of a function minus a constant. -/
private lemma avgQ_sub_const (f : Cfg N → ℝ) (c : ℝ) : avgQ (fun q => f q - c) = avgQ f - c := by
  have : avgQ (fun q => f q - c) = avgQ (fun q => f q + (-c)) := by
    simp only [sub_eq_add_neg]
  rw [this, avgQ_add, avgQ_const]
  ring

/-- Monotonicity of the average. -/
private lemma avgQ_mono {f g : Cfg N → ℝ} (h : ∀ q, f q ≤ g q) : avgQ f ≤ avgQ g := by
  simp only [avgQ]
  exact mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun q _ => h q) (by positivity)

/-- Positivity of the average. -/
private lemma avgQ_nonneg {f : Cfg N → ℝ} (h : ∀ q, 0 ≤ f q) : 0 ≤ avgQ f := by
  simp only [avgQ]
  exact mul_nonneg (by positivity) (Finset.sum_nonneg fun q _ => h q)

/-- `|avg f| ≤ avg |f|`. -/
private lemma abs_avgQ_le (f : Cfg N → ℝ) : |avgQ f| ≤ avgQ (fun q => |f q|) := by
  simp only [avgQ]
  rw [abs_mul, abs_of_pos (by positivity)]
  exact mul_le_mul_of_nonneg_left (Finset.abs_sum_le_sum_abs _ _) (by positivity)

/-- Real part of a complex Bernoulli average of real numbers. -/
private lemma re_avg (f : Cfg N → ℝ) :
    ((1 / 2 ^ N : ℂ) * ∑ q : Cfg N, ((f q : ℝ) : ℂ)).re = avgQ f := by
  simp only [avgQ]
  rw [← Complex.ofReal_sum]
  have : (1 / 2 ^ N : ℂ) = ((1 / 2 ^ N : ℝ) : ℂ) := by push_cast; ring
  rw [this, ← Complex.ofReal_mul, Complex.ofReal_re]

/-- **Cauchy–Schwarz** for the Bernoulli average: `avg(√u √v) ≤ √(avg u) √(avg v)`. -/
private lemma avgQ_sqrt_mul_sqrt_le {u v : Cfg N → ℝ} (hu : ∀ q, 0 ≤ u q) (hv : ∀ q, 0 ≤ v q) :
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
private lemma avgQ_sq_sum_le {κ : Type*} (s : Finset κ) (a : κ → Cfg N → ℝ) (b : κ → ℝ)
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
private lemma sum_comp_le_of_card_fiber {α : Type*} (s : Finset α) (φ : α → ℕ) (M : ℕ)
    (hφ : ∀ a ∈ s, φ a < M) (c : ℕ → ℝ)
    (hc : ∀ d, ((s.filter (fun a => φ a = d)).card : ℝ) ≤ c d)
    (g : ℕ → ℝ) (hg : ∀ d, 0 ≤ g d) :
    ∑ a ∈ s, g (φ a) ≤ ∑ d ∈ range M, c d * g d := by
  rw [← Finset.sum_fiberwise_of_maps_to' (t := range M) (fun a ha => mem_range.2 (hφ a ha))]
  refine Finset.sum_le_sum (fun d _ => ?_)
  rw [Finset.sum_const, nsmul_eq_mul]
  exact mul_le_mul_of_nonneg_right (hc d) (hg d)

/-- At most two sites are at a given distance from a site. -/
private lemma card_dist_fiber_le (z : Fin N) (d : ℕ) :
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
private lemma card_cross_fiber_le (m d : ℕ) :
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
private lemma sum_dist_le (g : ℕ → ℝ) (hg : ∀ d, 0 ≤ g d) (hs : Summable g) (z : Fin N) :
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
private lemma sum_cross_le (g : ℕ → ℝ) (hg : ∀ d, 0 ≤ g d)
    (hs : Summable (fun d : ℕ => (d : ℝ) * g d)) (m : ℕ) :
    ∑ x ∈ leftOf (N := N) m, ∑ y ∈ univ \ leftOf (N := N) m, g (Nat.dist x.val y.val) ≤
      ∑' d : ℕ, (d : ℝ) * g d := by
  rw [← Finset.sum_product' (f := fun x y : Fin N => g (Nat.dist x.val y.val))]
  calc _ ≤ ∑ d ∈ range N, (d : ℝ) * g d :=
        sum_comp_le_of_card_fiber _ (fun p : Fin N × Fin N => Nat.dist p.1.val p.2.val) N
          (fun p _ => by have := p.1.isLt; have := p.2.isLt; unfold Nat.dist; omega)
          (fun d => d) (card_cross_fiber_le m) g hg
    _ ≤ _ := hs.sum_le_tsum _ (fun d _ => mul_nonneg (Nat.cast_nonneg d) (hg d))

/-! #### Free-fermion computations in a fixed sector -/

/-- `⟨n_x⟩` of an evolved Fock state, as a real sum over configurations. -/
private lemma sum_occ_normSq (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (x : Fin N) (t : ℝ) :
    ∑ n, occ (n x) * Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) =
      ∑ z ∈ S, Complex.normSq (evo A t x z) := by
  have := density A hA S x t
  rw [expect_cplx_diagonal] at this
  exact_mod_cast this

/-- `⟨N_R⟩` of an evolved Fock state: `∑_{x ∈ R} ∑_{z ∈ S} |U_{xz}|²`. -/
private lemma sum_count_normSq (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S R : Finset (Fin N))
    (t : ℝ) :
    ∑ n, count R n * Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) =
      ∑ x ∈ R, ∑ z ∈ S, Complex.normSq (evo A t x z) := by
  simp only [count, Finset.sum_mul]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl (fun x _ => sum_occ_normSq A hA S x t)

/-- The evolved Fock state is normalized. -/
private lemma sum_normSq_evo_fock (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm)
    (S : Finset (Fin N)) (t : ℝ) : ∑ n, Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) = 1 := by
  rw [sum_normSq_mulVec _ (evo_unitary _ (dGamma_isSymm A hA) t), fockVec,
    Finset.sum_eq_single (fockCfg S) (fun n _ hn => by rw [Pi.single_eq_of_ne hn, map_zero])
      (fun h => absurd (mem_univ _) h), Pi.single_eq_same, map_one]

/-- LGT expectation of a diagonal matter observable, as a Bernoulli average of sector sums. -/
private lemma lgt_expect_diag_re (J h : ℝ) (ψ : Cfg N → ℂ) (d : Cfg N → ℝ) (t : ℝ) :
    (expect (matter (diagonal d)) (lgtState J h ψ t)).re =
      avgQ (fun q => ∑ n, d n * Complex.normSq ((evo (dGamma (hq J h q)) t *ᵥ ψ) n)) := by
  rw [matter, sector_average_matter]
  simp_rw [expect_cplx_diagonal]
  exact re_avg _

/-- The particle number of `R` in the configuration `n`, as a natural number. -/
private def countN (R : Finset (Fin N)) (n : Cfg N) : ℕ := (R.filter (fun x => n x ≠ 0)).card

/-- `count` is the natural number `countN`. -/
private lemma count_eq_countN (R : Finset (Fin N)) (n : Cfg N) : count R n = countN R n := by
  unfold count countN occ
  rw [Finset.card_filter]
  push_cast
  refine Finset.sum_congr rfl (fun x _ => ?_)
  split_ifs <;> simp_all

/-- The particle number is at most `N`. -/
private lemma countN_le (R : Finset (Fin N)) (n : Cfg N) : countN R n ≤ N :=
  (Finset.card_filter_le _ _).trans ((Finset.card_le_univ R).trans (by simp))

/-- Summing over the level sets `{count R = k}`, `k ∈ {0, …, N}`, is summing over everything. -/
private lemma sum_fiber_count (R : Finset (Fin N)) (φ : Cfg N → ℝ) :
    ∑ k ∈ range (N + 1), ∑ n ∈ univ.filter (fun n => count R n = k), φ n = ∑ n, φ n := by
  have h : ∀ k : ℕ, univ.filter (fun n : Cfg N => count R n = k) =
      univ.filter (fun n => countN R n = k) := by
    intro k; ext n; simp [count_eq_countN]
  simp_rw [h]
  exact Finset.sum_fiberwise_of_maps_to
    (fun n _ => mem_range.2 (Nat.lt_succ_of_le (countN_le R n))) φ

/-- The number distribution of a normalized vector is a probability distribution. -/
private lemma numDist_isProb (R : Finset (Fin N)) (v : Cfg N → ℂ)
    (hv : ∑ n, Complex.normSq (v n) = 1) :
    IsProb N (numDist R v) where
  nonneg k := Finset.sum_nonneg (fun n _ => Complex.normSq_nonneg _)
  sum_one := by unfold numDist; rw [sum_fiber_count]; exact hv

/-- Second moment about `c` of the number distribution. -/
private lemma sum_numDist_mul_sq (R : Finset (Fin N)) (v : Cfg N → ℂ) (c : ℝ) :
    ∑ k ∈ range (N + 1), numDist R v k * ((k : ℝ) - c) ^ 2 =
      ∑ n, Complex.normSq (v n) * (count R n - c) ^ 2 := by
  unfold numDist
  simp_rw [Finset.sum_mul]
  rw [← sum_fiber_count R]
  refine Finset.sum_congr rfl (fun k _ => Finset.sum_congr rfl (fun n hn => ?_))
  rw [(mem_filter.1 hn).2]

/-- Deviation of the mean particle number of `R` from its initial value `|S ∩ R|`, bounded by the
transition probabilities across the boundary of `R` (uses only column normalization). -/
private lemma abs_mean_sub_le (U : Matrix (Fin N) (Fin N) ℂ) (hU : Uᴴ * U = 1)
    (S R : Finset (Fin N)) :
    |∑ x ∈ R, ∑ z ∈ S, Complex.normSq (U x z) - ((S ∩ R).card : ℝ)| ≤
      ∑ z ∈ R, ∑ x ∈ univ \ R, Complex.normSq (U x z) +
        ∑ z ∈ univ \ R, ∑ x ∈ R, Complex.normSq (U x z) := by
  have ha : ∀ x z, 0 ≤ Complex.normSq (U x z) := fun x z => Complex.normSq_nonneg _
  have hcol : ∀ z,
      ∑ x ∈ R, Complex.normSq (U x z) + ∑ x ∈ univ \ R, Complex.normSq (U x z) = 1 := by
    intro z
    rw [add_comm, Finset.sum_sdiff (subset_univ R)]
    exact sum_normSq_col U hU z
  have hcard : ((S ∩ R).card : ℝ) = ∑ z ∈ S, if z ∈ R then (1 : ℝ) else 0 := by
    rw [Finset.sum_boole, Finset.filter_mem_eq_inter]
  have h1 : ∑ x ∈ R, ∑ z ∈ S, Complex.normSq (U x z) - ((S ∩ R).card : ℝ) =
      ∑ z ∈ S, (∑ x ∈ R, Complex.normSq (U x z) - if z ∈ R then 1 else 0) := by
    rw [hcard, Finset.sum_comm, ← Finset.sum_sub_distrib]
  have h2 : ∀ z, |∑ x ∈ R, Complex.normSq (U x z) - if z ∈ R then 1 else 0| ≤
      if z ∈ R then ∑ x ∈ univ \ R, Complex.normSq (U x z)
      else ∑ x ∈ R, Complex.normSq (U x z) := by
    intro z
    split_ifs with hz
    · have := hcol z
      rw [show ∑ x ∈ R, Complex.normSq (U x z) - 1 = -(∑ x ∈ univ \ R, Complex.normSq (U x z)) by
        linarith, abs_neg]
      exact le_of_eq (abs_of_nonneg (Finset.sum_nonneg fun x _ => ha x z))
    · rw [sub_zero]
      exact le_of_eq (abs_of_nonneg (Finset.sum_nonneg fun x _ => ha x z))
  have hd0 : ∀ z, 0 ≤ (if z ∈ R then ∑ x ∈ univ \ R, Complex.normSq (U x z)
      else ∑ x ∈ R, Complex.normSq (U x z)) := by
    intro z; split_ifs <;> exact Finset.sum_nonneg fun x _ => ha x z
  rw [h1]
  calc |∑ z ∈ S, (∑ x ∈ R, Complex.normSq (U x z) - if z ∈ R then 1 else 0)|
      ≤ ∑ z ∈ S, |∑ x ∈ R, Complex.normSq (U x z) - if z ∈ R then 1 else 0| :=
        Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ z ∈ S, (if z ∈ R then ∑ x ∈ univ \ R, Complex.normSq (U x z)
          else ∑ x ∈ R, Complex.normSq (U x z)) := Finset.sum_le_sum fun z _ => h2 z
    _ ≤ ∑ z, (if z ∈ R then ∑ x ∈ univ \ R, Complex.normSq (U x z)
          else ∑ x ∈ R, Complex.normSq (U x z)) :=
        Finset.sum_le_sum_of_subset_of_nonneg (subset_univ S) (fun z _ _ => hd0 z)
    _ = ∑ z ∈ R, ∑ x ∈ univ \ R, Complex.normSq (U x z) +
          ∑ z ∈ univ \ R, ∑ x ∈ R, Complex.normSq (U x z) := by
        rw [← Finset.sum_sdiff (subset_univ R), add_comm]
        congr 1
        · exact Finset.sum_congr rfl (fun z hz => ite_eq_left hz)
        · exact Finset.sum_congr rfl (fun z hz => ite_eq_right (Finset.mem_sdiff.1 hz).2)

/-! #### Consequences of Hypothesis (EC) -/

/-- Summability of `d ↦ d F(d)`. -/
private lemma ec_summable_mul {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) :
    Summable (fun d : ℕ => (d : ℝ) * F d) :=
  hEC.summable.of_nonneg_of_le (fun d => mul_nonneg (Nat.cast_nonneg d) (hEC.nonneg d))
    (fun d => mul_le_mul_of_nonneg_right
      (by nlinarith [(Nat.cast_nonneg d : (0 : ℝ) ≤ d)]) (hEC.nonneg d))

/-- Summability of `d ↦ √F(d)`. -/
private lemma ec_summable_sqrt {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) :
    Summable (fun d : ℕ => √(F d)) :=
  hEC.summable_sqrt.of_nonneg_of_le (fun d => Real.sqrt_nonneg _)
    (fun d => le_mul_of_one_le_left (Real.sqrt_nonneg _)
      (one_le_pow₀ (by linarith [(Nat.cast_nonneg d : (0 : ℝ) ≤ d)])))

/-- Summability of `d ↦ d √F(d)`. -/
private lemma ec_summable_mul_sqrt {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) :
    Summable (fun d : ℕ => (d : ℝ) * √(F d)) :=
  hEC.summable_sqrt.of_nonneg_of_le (fun d => mul_nonneg (Nat.cast_nonneg d) (Real.sqrt_nonneg _))
    (fun d => mul_le_mul_of_nonneg_right
      (by nlinarith [(Nat.cast_nonneg d : (0 : ℝ) ≤ d)]) (Real.sqrt_nonneg _))

/-- Averaged propagator bound from Hypothesis (EC). -/
private lemma avg_norm_evo_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (t : ℝ) (x z : Fin N) :
    avgQ (fun q => ‖evo (hq J h q) t x z‖) ≤ F (Nat.dist x.val z.val) :=
  (avgQ_mono (fun q => norm_evo_le_corrQ _ (hq_isHermitian J h q) t x z)).trans (hEC.bound N x z)

/-- `|U_{xz}|² ≤ |U_{xz}|` for the unitary propagator. -/
private lemma normSq_evo_le_norm (J h : ℝ) (q : Cfg N) (t : ℝ) (x z : Fin N) :
    Complex.normSq (evo (hq J h q) t x z) ≤ ‖evo (hq J h q) t x z‖ := by
  rw [Complex.normSq_eq_norm_sq]
  have h1 := norm_evo_le_one _ (hq_isSymm J h q) t x z
  have h0 := norm_nonneg (evo (hq J h q) t x z)
  nlinarith

/-- `|U_{xz}|⁴ ≤ |U_{xz}|` for the unitary propagator. -/
private lemma normSq_sq_evo_le_norm (J h : ℝ) (q : Cfg N) (t : ℝ) (x z : Fin N) :
    Complex.normSq (evo (hq J h q) t x z) ^ 2 ≤ ‖evo (hq J h q) t x z‖ := by
  rw [Complex.normSq_eq_norm_sq]
  have h1 := norm_evo_le_one _ (hq_isSymm J h q) t x z
  have h0 := norm_nonneg (evo (hq J h q) t x z)
  have h3 : ‖evo (hq J h q) t x z‖ ^ 3 ≤ 1 := pow_le_one₀ h0 h1
  nlinarith [mul_le_mul_of_nonneg_left h3 h0]

end MainAux

/-- **Exact density formula**: `⟨n_x(t)⟩ = 2^{-N} ∑_q ∑_{z ∈ S} |(e^{-ith_q})_{xz}|²`. -/
theorem lgt_density (J h : ℝ) (S : Finset (Fin N)) (x : Fin N) (t : ℝ) :
    (expect (matter (diagonal fun n => occ (n x))) (lgtState J h (fockVec S) t)).re =
      avgQ (fun q => ∑ z ∈ S, Complex.normSq (evo (hq J h q) t x z)) := by
  rw [MainAux.lgt_expect_diag_re]
  congr 1
  funext q
  exact MainAux.sum_occ_normSq _ (hq_isSymm J h q) S x t

/-- **Bounded transport** (paper: Theorem `thm:transport`): under Hypothesis (EC), the number of
particles transferred across any cut is bounded uniformly in time and volume. -/
theorem lgt_transport (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ)
    (t : ℝ) :
    |(expect (matter (numOp (leftOf m))) (lgtState J h (fockVec S) t)).re -
        ((S ∩ leftOf m).card : ℝ)| ≤
      2 * ∑' d : ℕ, (d : ℝ) * F d := by
  have hF0 := hEC.nonneg
  have hpair : ∀ x z : Fin N, avgQ (fun q => Complex.normSq (evo (hq J h q) t x z)) ≤
      F (Nat.dist x.val z.val) := fun x z =>
    (MainAux.avgQ_mono (fun q => MainAux.normSq_evo_le_norm J h q t x z)).trans
      (MainAux.avg_norm_evo_le hEC t x z)
  have hB := MainAux.sum_cross_le (N := N) F hF0 (MainAux.ec_summable_mul hEC) m
  rw [numOp, MainAux.lgt_expect_diag_re]
  have e1 : (fun q => ∑ n, count (leftOf m) n *
      Complex.normSq ((evo (dGamma (hq J h q)) t *ᵥ fockVec S) n)) =
      fun q => ∑ x ∈ leftOf m, ∑ z ∈ S, Complex.normSq (evo (hq J h q) t x z) := by
    funext q; exact MainAux.sum_count_normSq _ (hq_isSymm J h q) S _ t
  rw [e1, ← MainAux.avgQ_sub_const]
  calc |avgQ (fun q => ∑ x ∈ leftOf m, ∑ z ∈ S, Complex.normSq (evo (hq J h q) t x z) -
        ((S ∩ leftOf m).card : ℝ))|
      ≤ avgQ (fun q => |∑ x ∈ leftOf m, ∑ z ∈ S, Complex.normSq (evo (hq J h q) t x z) -
        ((S ∩ leftOf m).card : ℝ)|) := MainAux.abs_avgQ_le _
    _ ≤ avgQ (fun q => ∑ z ∈ leftOf m, ∑ x ∈ univ \ leftOf m,
          Complex.normSq (evo (hq J h q) t x z) +
          ∑ z ∈ univ \ leftOf m, ∑ x ∈ leftOf m, Complex.normSq (evo (hq J h q) t x z)) :=
        MainAux.avgQ_mono (fun q =>
          MainAux.abs_mean_sub_le _ (evo_unitary _ (hq_isSymm J h q) t) S _)
    _ = ∑ z ∈ leftOf m, ∑ x ∈ univ \ leftOf m,
          avgQ (fun q => Complex.normSq (evo (hq J h q) t x z)) +
          ∑ z ∈ univ \ leftOf m, ∑ x ∈ leftOf m,
            avgQ (fun q => Complex.normSq (evo (hq J h q) t x z)) := by
        simp only [MainAux.avgQ_add, MainAux.avgQ_sum]
    _ ≤ ∑ z ∈ leftOf m, ∑ x ∈ univ \ leftOf m, F (Nat.dist x.val z.val) +
          ∑ z ∈ univ \ leftOf m, ∑ x ∈ leftOf m, F (Nat.dist x.val z.val) :=
        add_le_add (Finset.sum_le_sum fun z _ => Finset.sum_le_sum fun x _ => hpair x z)
          (Finset.sum_le_sum fun z _ => Finset.sum_le_sum fun x _ => hpair x z)
    _ ≤ (∑' d : ℕ, (d : ℝ) * F d) + ∑' d : ℕ, (d : ℝ) * F d := by
        refine add_le_add ?_ ?_
        · calc _ = ∑ z ∈ leftOf m, ∑ x ∈ univ \ leftOf m, F (Nat.dist z.val x.val) :=
                Finset.sum_congr rfl fun z _ => Finset.sum_congr rfl fun x _ => by
                  rw [Nat.dist_comm]
            _ ≤ _ := hB
        · rw [Finset.sum_comm]
          exact hB
    _ = 2 * ∑' d : ℕ, (d : ℝ) * F d := by ring

/-- Distribution of the particle number `N_R` in the LGT state at time `t`. -/
def lgtNumDist (J h : ℝ) (R : Finset (Fin N)) (ψ : Cfg N → ℂ) (t : ℝ) (k : ℕ) : ℝ :=
  (expect (matter (numProj R k)) (lgtState J h ψ t)).re

namespace MainAux

/-- The LGT number distribution is the Bernoulli average of the sector number distributions. -/
private lemma lgtNumDist_eq_avgQ (J h : ℝ) (R : Finset (Fin N)) (ψ : Cfg N → ℂ) (t : ℝ) (k : ℕ) :
    lgtNumDist J h R ψ t k = avgQ (fun q => numDist R (evo (dGamma (hq J h q)) t *ᵥ ψ) k) := by
  rw [lgtNumDist, numProj, lgt_expect_diag_re]
  congr 1
  funext q
  rw [numDist, Finset.sum_filter]
  refine Finset.sum_congr rfl (fun n _ => ?_)
  split_ifs <;> ring

end MainAux

/-- The number distribution is a probability distribution on `{0, …, N}`. -/
theorem lgtNumDist_isProb (J h : ℝ) (R S : Finset (Fin N)) (t : ℝ) :
    IsProb N (lgtNumDist J h R (fockVec S) t) := by
  refine ⟨fun k => ?_, ?_⟩
  · rw [MainAux.lgtNumDist_eq_avgQ]
    exact MainAux.avgQ_nonneg (fun q =>
      (MainAux.numDist_isProb R _ (MainAux.sum_normSq_evo_fock _ (hq_isSymm J h q) S t)).nonneg k)
  · simp_rw [MainAux.lgtNumDist_eq_avgQ]
    rw [← MainAux.avgQ_sum]
    have : ∀ q : Cfg N, ∑ k ∈ range (N + 1),
        numDist R (evo (dGamma (hq J h q)) t *ᵥ fockVec S) k = 1 := fun q =>
      (MainAux.numDist_isProb R _ (MainAux.sum_normSq_evo_fock _ (hq_isSymm J h q) S t)).sum_one
    simp_rw [this]
    exact MainAux.avgQ_const 1

namespace MainAux

/-! #### The number variance -/

/-- `a ≤ √a` on `[0, 1]`. -/
private lemma le_sqrt_self_of_le_one {a : ℝ} (h0 : 0 ≤ a) (h1 : a ≤ 1) : a ≤ √a :=
  Real.le_sqrt_of_sq_le (by nlinarith)

/-- Entry bound for the correlation matrix `C = U P_S U†`:
`|C_{xy}|² ≤ ∑_z √|U_{xz}| √|U_{yz}|`. -/
private lemma normSq_corrMat_le (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (x y : Fin N) :
    Complex.normSq (corrMat A S t x y) ≤ ∑ z, √‖evo A t x z‖ * √‖evo A t y z‖ := by
  have hrow : ∀ x, ∑ z, ‖evo A t x z‖ ^ 2 = 1 := fun x => by
    simpa [Complex.normSq_eq_norm_sq] using sum_normSq_row (evo A t) (evo_unitary' A hA t) x
  have hle1 : ∀ x z, ‖evo A t x z‖ ≤ 1 := norm_evo_le_one A hA t
  have hC : corrMat A S t x y =
      ∑ z, evo A t x z * (if z ∈ S then 1 else 0) * star (evo A t y z) := by
    have e : corrMat A S t x y = ∑ z, (evo A t * diagonal (fun z => if z ∈ S then (1 : ℂ) else 0))
        x z * (evo A t)ᴴ z y := Matrix.mul_apply
    rw [e]
    simp only [mul_diagonal, conjTranspose_apply]
  have h1 : ‖corrMat A S t x y‖ ≤ ∑ z, ‖evo A t x z‖ * ‖evo A t y z‖ := by
    rw [hC]
    refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun z _ => ?_)
    rw [norm_mul, norm_mul, norm_star]
    have : ‖(if z ∈ S then (1 : ℂ) else 0)‖ ≤ 1 := by split_ifs <;> simp
    calc ‖evo A t x z‖ * ‖(if z ∈ S then (1 : ℂ) else 0)‖ * ‖evo A t y z‖
        ≤ ‖evo A t x z‖ * 1 * ‖evo A t y z‖ := by gcongr
      _ = ‖evo A t x z‖ * ‖evo A t y z‖ := by ring
  have h2 : ∑ z, ‖evo A t x z‖ * ‖evo A t y z‖ ≤ 1 := by
    have : ∀ z, ‖evo A t x z‖ * ‖evo A t y z‖ ≤ (‖evo A t x z‖ ^ 2 + ‖evo A t y z‖ ^ 2) / 2 :=
      fun z => by nlinarith [sq_nonneg (‖evo A t x z‖ - ‖evo A t y z‖)]
    calc ∑ z, ‖evo A t x z‖ * ‖evo A t y z‖
        ≤ ∑ z, (‖evo A t x z‖ ^ 2 + ‖evo A t y z‖ ^ 2) / 2 := Finset.sum_le_sum fun z _ => this z
      _ = (∑ z, ‖evo A t x z‖ ^ 2 + ∑ z, ‖evo A t y z‖ ^ 2) / 2 := by
        rw [← Finset.sum_add_distrib, Finset.sum_div]
      _ = 1 := by rw [hrow, hrow]; norm_num
  have h3 : ∀ z, ‖evo A t x z‖ * ‖evo A t y z‖ ≤ √‖evo A t x z‖ * √‖evo A t y z‖ := fun z =>
    mul_le_mul (le_sqrt_self_of_le_one (norm_nonneg _) (hle1 x z))
      (le_sqrt_self_of_le_one (norm_nonneg _) (hle1 y z)) (norm_nonneg _) (Real.sqrt_nonneg _)
  rw [Complex.normSq_eq_norm_sq]
  have hC0 := norm_nonneg (corrMat A S t x y)
  have hC1 := h1.trans h2
  calc ‖corrMat A S t x y‖ ^ 2 ≤ ‖corrMat A S t x y‖ := by nlinarith
    _ ≤ ∑ z, ‖evo A t x z‖ * ‖evo A t y z‖ := h1
    _ ≤ ∑ z, √‖evo A t x z‖ * √‖evo A t y z‖ := Finset.sum_le_sum fun z _ => h3 z

/-- The combinatorial bound behind the variance estimate:
`∑_{x ≤ m < y} ∑_z f(|x-z|) f(|y-z|) ≤ 4 (∑_d f(d)) (∑_d d f(d))`. -/
private lemma triple_sum_le (f : ℕ → ℝ) (hf : ∀ d, 0 ≤ f d) (hs : Summable f)
    (hs' : Summable (fun d : ℕ => (d : ℝ) * f d)) (m : ℕ) :
    ∑ x ∈ leftOf (N := N) m, ∑ y ∈ univ \ leftOf (N := N) m, ∑ z : Fin N,
        f (Nat.dist x.val z.val) * f (Nat.dist y.val z.val) ≤
      4 * (∑' d, f d) * ∑' d : ℕ, (d : ℝ) * f d := by
  have hA0 : 0 ≤ ∑' d, f d := tsum_nonneg hf
  have hB := sum_cross_le (N := N) f hf hs' m
  have hP : ∀ (z : Fin N) (s : Finset (Fin N)), ∑ x ∈ s, f (Nat.dist x.val z.val) ≤
      2 * ∑' d, f d := fun z s =>
    (Finset.sum_le_sum_of_subset_of_nonneg (subset_univ s) (fun x _ _ => hf _)).trans
      (sum_dist_le f hf hs z)
  have hnn : ∀ (z : Fin N) (s : Finset (Fin N)), 0 ≤ ∑ x ∈ s, f (Nat.dist x.val z.val) :=
    fun z s => Finset.sum_nonneg fun x _ => hf _
  have e1 : ∑ x ∈ leftOf (N := N) m, ∑ y ∈ univ \ leftOf (N := N) m, ∑ z : Fin N,
        f (Nat.dist x.val z.val) * f (Nat.dist y.val z.val) =
      ∑ z : Fin N, (∑ x ∈ leftOf (N := N) m, f (Nat.dist x.val z.val)) *
        ∑ y ∈ univ \ leftOf (N := N) m, f (Nat.dist y.val z.val) := by
    calc _ = ∑ x ∈ leftOf (N := N) m, ∑ z : Fin N, ∑ y ∈ univ \ leftOf (N := N) m,
          f (Nat.dist x.val z.val) * f (Nat.dist y.val z.val) :=
          Finset.sum_congr rfl fun x _ => Finset.sum_comm
      _ = ∑ z : Fin N, ∑ x ∈ leftOf (N := N) m, ∑ y ∈ univ \ leftOf (N := N) m,
          f (Nat.dist x.val z.val) * f (Nat.dist y.val z.val) := Finset.sum_comm
      _ = _ := Finset.sum_congr rfl fun z _ => (Finset.sum_mul_sum _ _ _ _).symm
  rw [e1, ← Finset.sum_sdiff (subset_univ (leftOf (N := N) m))]
  have h1 : ∑ z ∈ univ \ leftOf (N := N) m, (∑ x ∈ leftOf (N := N) m, f (Nat.dist x.val z.val)) *
      ∑ y ∈ univ \ leftOf (N := N) m, f (Nat.dist y.val z.val) ≤
      2 * (∑' d, f d) * ∑' d : ℕ, (d : ℝ) * f d := by
    calc _ ≤ ∑ z ∈ univ \ leftOf (N := N) m,
          (∑ x ∈ leftOf (N := N) m, f (Nat.dist x.val z.val)) * (2 * ∑' d, f d) :=
          Finset.sum_le_sum fun z _ => mul_le_mul_of_nonneg_left (hP z _) (hnn z _)
      _ = 2 * (∑' d, f d) * ∑ z ∈ univ \ leftOf (N := N) m,
          ∑ x ∈ leftOf (N := N) m, f (Nat.dist x.val z.val) := by
          rw [Finset.mul_sum]
          exact Finset.sum_congr rfl fun z _ => by ring
      _ = 2 * (∑' d, f d) * ∑ x ∈ leftOf (N := N) m,
          ∑ z ∈ univ \ leftOf (N := N) m, f (Nat.dist x.val z.val) := by
          rw [Finset.sum_comm]
      _ ≤ _ := mul_le_mul_of_nonneg_left hB (by positivity)
  have h2 : ∑ z ∈ leftOf (N := N) m, (∑ x ∈ leftOf (N := N) m, f (Nat.dist x.val z.val)) *
      ∑ y ∈ univ \ leftOf (N := N) m, f (Nat.dist y.val z.val) ≤
      2 * (∑' d, f d) * ∑' d : ℕ, (d : ℝ) * f d := by
    calc _ ≤ ∑ z ∈ leftOf (N := N) m,
          (2 * ∑' d, f d) * ∑ y ∈ univ \ leftOf (N := N) m, f (Nat.dist y.val z.val) :=
          Finset.sum_le_sum fun z _ => mul_le_mul_of_nonneg_right (hP z _) (hnn z _)
      _ = 2 * (∑' d, f d) * ∑ z ∈ leftOf (N := N) m,
          ∑ y ∈ univ \ leftOf (N := N) m, f (Nat.dist z.val y.val) := by
          rw [Finset.mul_sum]
          refine Finset.sum_congr rfl fun z _ => ?_
          congr 1
          exact Finset.sum_congr rfl fun y _ => by rw [Nat.dist_comm]
      _ ≤ _ := mul_le_mul_of_nonneg_left hB (by positivity)
  linarith

/-- Averaged half-chain number variance under (EC). -/
private lemma avg_var_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ)
    (t : ℝ) :
    avgQ (fun q => ∑ x ∈ leftOf (N := N) m, ∑ y ∈ univ \ leftOf (N := N) m,
        Complex.normSq (corrMat (hq J h q) S t x y)) ≤
      4 * (∑' d, √(F d)) * ∑' d : ℕ, (d : ℝ) * √(F d) := by
  have hpair : ∀ x y z : Fin N,
      avgQ (fun q => √‖evo (hq J h q) t x z‖ * √‖evo (hq J h q) t y z‖) ≤
        √(F (Nat.dist x.val z.val)) * √(F (Nat.dist y.val z.val)) := fun x y z =>
    (avgQ_sqrt_mul_sqrt_le (fun q => norm_nonneg _) (fun q => norm_nonneg _)).trans
      (mul_le_mul (Real.sqrt_le_sqrt (avg_norm_evo_le hEC t x z))
        (Real.sqrt_le_sqrt (avg_norm_evo_le hEC t y z)) (Real.sqrt_nonneg _) (Real.sqrt_nonneg _))
  calc _ ≤ avgQ (fun q => ∑ x ∈ leftOf (N := N) m, ∑ y ∈ univ \ leftOf (N := N) m, ∑ z,
        √‖evo (hq J h q) t x z‖ * √‖evo (hq J h q) t y z‖) :=
        avgQ_mono fun q => Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ =>
          normSq_corrMat_le _ (hq_isSymm J h q) S t x y
    _ = ∑ x ∈ leftOf (N := N) m, ∑ y ∈ univ \ leftOf (N := N) m, ∑ z,
        avgQ (fun q => √‖evo (hq J h q) t x z‖ * √‖evo (hq J h q) t y z‖) := by
        simp only [avgQ_sum]
    _ ≤ ∑ x ∈ leftOf (N := N) m, ∑ y ∈ univ \ leftOf (N := N) m, ∑ z,
        √(F (Nat.dist x.val z.val)) * √(F (Nat.dist y.val z.val)) :=
        Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ => Finset.sum_le_sum fun z _ =>
          hpair x y z
    _ ≤ _ := triple_sum_le (fun d => √(F d)) (fun d => Real.sqrt_nonneg _) (ec_summable_sqrt hEC)
        (ec_summable_mul_sqrt hEC) m

/-- Averaged squared deviation of the mean half-chain particle number under (EC). -/
private lemma avg_mean_sq_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ)
    (t : ℝ) :
    avgQ (fun q => (∑ x ∈ leftOf (N := N) m, ∑ z ∈ S, Complex.normSq (evo (hq J h q) t x z) -
        ((S ∩ leftOf (N := N) m).card : ℝ)) ^ 2) ≤
      4 * (∑' d : ℕ, (d : ℝ) * √(F d)) ^ 2 := by
  have hf0 : ∀ d, 0 ≤ √(F d) := fun d => Real.sqrt_nonneg _
  have hB := sum_cross_le (N := N) (fun d => √(F d)) hf0 (ec_summable_mul_sqrt hEC) m
  have hpair : ∀ x z : Fin N,
      √(avgQ (fun q => Complex.normSq (evo (hq J h q) t x z) ^ 2)) ≤
        √(F (Nat.dist x.val z.val)) := fun x z =>
    Real.sqrt_le_sqrt ((avgQ_mono fun q => normSq_sq_evo_le_norm J h q t x z).trans
      (avg_norm_evo_le hEC t x z))
  have hD1 : avgQ (fun q => (∑ p ∈ leftOf (N := N) m ×ˢ (univ \ leftOf (N := N) m),
      Complex.normSq (evo (hq J h q) t p.2 p.1)) ^ 2) ≤
      (∑' d : ℕ, (d : ℝ) * √(F d)) ^ 2 := by
    refine (avgQ_sq_sum_le (leftOf (N := N) m ×ˢ (univ \ leftOf (N := N) m))
      (fun p q => Complex.normSq (evo (hq J h q) t p.2 p.1))
      (fun p => √(F (Nat.dist p.2.val p.1.val)))
      (fun p _ q => Complex.normSq_nonneg _) (fun p _ => hpair p.2 p.1)).trans ?_
    refine pow_le_pow_left₀ (Finset.sum_nonneg fun p _ => hf0 _) ?_ 2
    calc ∑ p ∈ leftOf (N := N) m ×ˢ (univ \ leftOf (N := N) m), √(F (Nat.dist p.2.val p.1.val))
        = ∑ z ∈ leftOf (N := N) m, ∑ x ∈ univ \ leftOf (N := N) m,
            √(F (Nat.dist z.val x.val)) := by
          rw [Finset.sum_product]
          exact Finset.sum_congr rfl fun z _ => Finset.sum_congr rfl fun x _ => by
            simp only [Nat.dist_comm x.val z.val]
      _ ≤ _ := hB
  have hD2 : avgQ (fun q => (∑ p ∈ leftOf (N := N) m ×ˢ (univ \ leftOf (N := N) m),
      Complex.normSq (evo (hq J h q) t p.1 p.2)) ^ 2) ≤
      (∑' d : ℕ, (d : ℝ) * √(F d)) ^ 2 := by
    refine (avgQ_sq_sum_le (leftOf (N := N) m ×ˢ (univ \ leftOf (N := N) m))
      (fun p q => Complex.normSq (evo (hq J h q) t p.1 p.2))
      (fun p => √(F (Nat.dist p.1.val p.2.val)))
      (fun p _ q => Complex.normSq_nonneg _) (fun p _ => hpair p.1 p.2)).trans ?_
    refine pow_le_pow_left₀ (Finset.sum_nonneg fun p _ => hf0 _) ?_ 2
    rw [Finset.sum_product]
    exact hB
  have hpt : ∀ q : Cfg N, (∑ x ∈ leftOf (N := N) m, ∑ z ∈ S,
      Complex.normSq (evo (hq J h q) t x z) - ((S ∩ leftOf (N := N) m).card : ℝ)) ^ 2 ≤
      2 * (∑ p ∈ leftOf (N := N) m ×ˢ (univ \ leftOf (N := N) m),
        Complex.normSq (evo (hq J h q) t p.2 p.1)) ^ 2 +
      2 * (∑ p ∈ leftOf (N := N) m ×ˢ (univ \ leftOf (N := N) m),
        Complex.normSq (evo (hq J h q) t p.1 p.2)) ^ 2 := by
    intro q
    have hab := abs_mean_sub_le (evo (hq J h q) t) (evo_unitary _ (hq_isSymm J h q) t) S
      (leftOf (N := N) m)
    have e1 : ∑ p ∈ leftOf (N := N) m ×ˢ (univ \ leftOf (N := N) m),
        Complex.normSq (evo (hq J h q) t p.2 p.1) =
        ∑ z ∈ leftOf (N := N) m, ∑ x ∈ univ \ leftOf (N := N) m,
          Complex.normSq (evo (hq J h q) t x z) := Finset.sum_product _ _ _
    have e2 : ∑ p ∈ leftOf (N := N) m ×ˢ (univ \ leftOf (N := N) m),
        Complex.normSq (evo (hq J h q) t p.1 p.2) =
        ∑ z ∈ univ \ leftOf (N := N) m, ∑ x ∈ leftOf (N := N) m,
          Complex.normSq (evo (hq J h q) t x z) := by
      rw [Finset.sum_product, Finset.sum_comm]
    rw [e1, e2]
    have h1 : 0 ≤ ∑ z ∈ leftOf (N := N) m, ∑ x ∈ univ \ leftOf (N := N) m,
        Complex.normSq (evo (hq J h q) t x z) :=
      Finset.sum_nonneg fun _ _ => Finset.sum_nonneg fun _ _ => Complex.normSq_nonneg _
    have h2 : 0 ≤ ∑ z ∈ univ \ leftOf (N := N) m, ∑ x ∈ leftOf (N := N) m,
        Complex.normSq (evo (hq J h q) t x z) :=
      Finset.sum_nonneg fun _ _ => Finset.sum_nonneg fun _ _ => Complex.normSq_nonneg _
    have h3 := pow_le_pow_left₀ (abs_nonneg _) hab 2
    rw [sq_abs] at h3
    nlinarith [sq_nonneg (∑ z ∈ leftOf (N := N) m, ∑ x ∈ univ \ leftOf (N := N) m,
        Complex.normSq (evo (hq J h q) t x z) - ∑ z ∈ univ \ leftOf (N := N) m,
        ∑ x ∈ leftOf (N := N) m, Complex.normSq (evo (hq J h q) t x z))]
  calc _ ≤ _ := avgQ_mono hpt
    _ = 2 * avgQ (fun q => (∑ p ∈ leftOf (N := N) m ×ˢ (univ \ leftOf (N := N) m),
          Complex.normSq (evo (hq J h q) t p.2 p.1)) ^ 2) +
        2 * avgQ (fun q => (∑ p ∈ leftOf (N := N) m ×ˢ (univ \ leftOf (N := N) m),
          Complex.normSq (evo (hq J h q) t p.1 p.2)) ^ 2) := by
        rw [avgQ_add, avgQ_const_mul, avgQ_const_mul]
    _ ≤ 4 * (∑' d : ℕ, (d : ℝ) * √(F d)) ^ 2 := by linarith

/-- Per-sector second moment of the number distribution about `c`:
`∑_k p_q(k) (k - c)² = Var_q(N_R) + (⟨N_R⟩_q - c)²`. -/
private lemma sector_second_moment (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm)
    (S R : Finset (Fin N)) (t c : ℝ) :
    ∑ k ∈ range (N + 1), numDist R (evo (dGamma A) t *ᵥ fockVec S) k * ((k : ℝ) - c) ^ 2 =
      (∑ x ∈ R, ∑ y ∈ univ \ R, Complex.normSq (corrMat A S t x y)) +
        (∑ x ∈ R, ∑ z ∈ S, Complex.normSq (evo A t x z) - c) ^ 2 := by
  rw [sum_numDist_mul_sq]
  have hnorm := sum_normSq_evo_fock A hA S t
  have hM1 := sum_count_normSq A hA S R t
  have hvar := number_variance A hA S R t
  simp only [numOp, diagonal_mul_diagonal] at hvar
  rw [expect_cplx_diagonal, expect_cplx_diagonal] at hvar
  have hvar' : ∑ n, count R n * count R n *
        Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) -
      (∑ n, count R n * Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n)) ^ 2 =
      ∑ x ∈ R, ∑ y ∈ univ \ R, Complex.normSq (corrMat A S t x y) := by
    exact_mod_cast hvar
  have hexp : ∑ n, Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) * (count R n - c) ^ 2 =
      ∑ n, count R n * count R n * Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) -
        2 * c * ∑ n, count R n * Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) +
        c ^ 2 * ∑ n, Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) := by
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl (fun n _ => by ring)
  rw [hexp, hnorm, hM1]
  rw [hM1] at hvar'
  linear_combination hvar'

end MainAux

/-- **Bounded number entropy** (paper: Theorem `thm:numberentropy`): under Hypothesis (EC), the
Shannon entropy of the half-chain particle-number distribution is bounded uniformly in time,
volume, position of the cut and occupied set. -/
theorem lgt_numberEntropy (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ K : ℝ, ∀ (N : ℕ) (m : ℕ) (S : Finset (Fin N)) (t : ℝ),
      shannon N (lgtNumDist J h (leftOf m) (fockVec S) t) ≤
        (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (K + 1 / 12)) := by
  refine ⟨4 * (∑' d, √(F d)) * (∑' d : ℕ, (d : ℝ) * √(F d)) +
    4 * (∑' d : ℕ, (d : ℝ) * √(F d)) ^ 2, fun N m S t => ?_⟩
  have hp := lgtNumDist_isProb J h (leftOf m) S t
  refine shannon_le_of_pvar_le hp ?_
  calc pvar N (lgtNumDist J h (leftOf m) (fockVec S) t)
      ≤ ∑ k ∈ range (N + 1), lgtNumDist J h (leftOf m) (fockVec S) t k *
          ((k : ℝ) - ((S ∩ leftOf m).card : ℝ)) ^ 2 := pvar_le_sum_sq hp _
    _ = avgQ (fun q => ∑ k ∈ range (N + 1),
          numDist (leftOf m) (evo (dGamma (hq J h q)) t *ᵥ fockVec S) k *
            ((k : ℝ) - ((S ∩ leftOf m).card : ℝ)) ^ 2) := by
        rw [MainAux.avgQ_sum]
        refine Finset.sum_congr rfl (fun k _ => ?_)
        rw [MainAux.lgtNumDist_eq_avgQ, mul_comm, ← MainAux.avgQ_const_mul]
        congr 1
        funext q
        ring
    _ = avgQ (fun q => (∑ x ∈ leftOf m, ∑ y ∈ univ \ leftOf m,
            Complex.normSq (corrMat (hq J h q) S t x y)) +
          (∑ x ∈ leftOf m, ∑ z ∈ S, Complex.normSq (evo (hq J h q) t x z) -
            ((S ∩ leftOf m).card : ℝ)) ^ 2) := by
        congr 1
        funext q
        exact MainAux.sector_second_moment _ (hq_isSymm J h q) S (leftOf m) t _
    _ = _ := MainAux.avgQ_add _ _
    _ ≤ _ := add_le_add (MainAux.avg_var_le hEC S m t) (MainAux.avg_mean_sq_le hEC S m t)

/-- The density-wave imbalance `I(t) = (2/N) ∑_x (-1)^x ⟨n_x(t)⟩` of the LGT. -/
def lgtImbalance (J h : ℝ) (N : ℕ) (t : ℝ) : ℝ :=
  (2 / N : ℝ) * ∑ x : Fin N, (-1 : ℝ) ^ x.val *
    (expect (matter (diagonal fun n => occ (n x))) (lgtState J h (fockVec (cdw N)) t)).re

namespace MainAux

/-! #### The imbalance -/

/-- `∑_z (-1)^z a_z = 2 ∑_{z even} a_z - ∑_z a_z`. -/
private lemma sum_neg_one_pow_mul (a : Fin N → ℝ) :
    ∑ z : Fin N, (-1 : ℝ) ^ z.val * a z = 2 * ∑ z ∈ cdw N, a z - ∑ z, a z := by
  rw [cdw, Finset.sum_filter, Finset.mul_sum, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl (fun z _ => ?_)
  rcases Nat.even_or_odd z.val with he | ho
  · rw [he.neg_one_pow, ite_eq_left he]; ring
  · rw [ho.neg_one_pow, ite_eq_right (Nat.not_even_iff_odd.2 ho)]; ring

/-- `∑_{x < N} (-1)^x = 0` for even `N`. -/
private lemma sum_neg_one_pow_fin (hN : Even N) : ∑ x : Fin N, (-1 : ℝ) ^ x.val = 0 := by
  rw [Fin.sum_univ_eq_sum_range (fun i => (-1 : ℝ) ^ i) N]
  obtain ⟨k, rfl⟩ := hN
  induction k with
  | zero => simp
  | succ k ih =>
    have h2 : k + 1 + (k + 1) = (k + k) + 1 + 1 := by ring
    rw [h2, Finset.sum_range_succ, Finset.sum_range_succ, ih, pow_succ,
      Even.neg_one_pow ⟨k, rfl⟩]
    ring

/-- `Re Tr(Γ U Γ U†) = ∑_{x,z} (-1)^x (-1)^z |U_{xz}|²`. -/
private lemma trace_stag_re (U : Matrix (Fin N) (Fin N) ℂ) :
    (trace (cplx (stag N) * U * cplx (stag N) * Uᴴ)).re =
      ∑ x : Fin N, ∑ z : Fin N, (-1 : ℝ) ^ x.val * ((-1 : ℝ) ^ z.val * Complex.normSq (U x z)) := by
  have hs : cplx (stag N) = diagonal (fun x : Fin N => (((-1 : ℝ) ^ x.val : ℝ) : ℂ)) := by
    unfold cplx stag; exact diagonal_map (by simp)
  rw [hs, trace, Complex.re_sum]
  refine Finset.sum_congr rfl (fun x _ => ?_)
  rw [Matrix.diag_apply, Matrix.mul_apply, Complex.re_sum]
  refine Finset.sum_congr rfl (fun z _ => ?_)
  rw [mul_diagonal, diagonal_mul, conjTranspose_apply]
  have : (((-1 : ℝ) ^ x.val : ℝ) : ℂ) * U x z * (((-1 : ℝ) ^ z.val : ℝ) : ℂ) * star (U x z) =
      (((-1 : ℝ) ^ x.val * ((-1 : ℝ) ^ z.val * Complex.normSq (U x z)) : ℝ) : ℂ) := by
    rw [Complex.ofReal_mul, Complex.ofReal_mul, ← Complex.mul_conj]
    simp only [Complex.star_def]
    ring
  rw [this, Complex.ofReal_re]

/-- The one-body imbalance through the propagator (`N` even):
`imbOne = (2/N) ∑_x (-1)^x ∑_{z even} |U_{xz}|²`. -/
private lemma imbOne_eq (J h : ℝ) (hN : Even N) (q : Cfg N) (t : ℝ) :
    imbOne J h q t = (2 / N : ℝ) * ∑ x : Fin N, (-1 : ℝ) ^ x.val *
      ∑ z ∈ cdw N, Complex.normSq (evo (hq J h q) t x z) := by
  have hrow := sum_normSq_row (evo (hq J h q) t) (evo_unitary' _ (hq_isSymm J h q) t)
  have h1 : ∀ x : Fin N, ∑ z : Fin N, (-1 : ℝ) ^ x.val *
      ((-1 : ℝ) ^ z.val * Complex.normSq (evo (hq J h q) t x z)) =
      2 * ((-1 : ℝ) ^ x.val * ∑ z ∈ cdw N, Complex.normSq (evo (hq J h q) t x z)) -
        (-1 : ℝ) ^ x.val := by
    intro x
    rw [← Finset.mul_sum, sum_neg_one_pow_mul, hrow x]
    ring
  rw [imbOne, trace_stag_re, Finset.sum_congr rfl (fun x _ => h1 x), Finset.sum_sub_distrib,
    ← Finset.mul_sum, sum_neg_one_pow_fin hN]
  ring

end MainAux

/-- The LGT imbalance is the Bernoulli average of the one-body imbalances (`N` even). -/
theorem lgtImbalance_eq (J h : ℝ) (hN : Even N) (t : ℝ) :
    lgtImbalance J h N t = avgQ (fun q : Cfg N => imbOne J h q t) := by
  unfold lgtImbalance
  simp_rw [lgt_density, MainAux.imbOne_eq J h hN]
  rw [MainAux.avgQ_const_mul, MainAux.avgQ_sum]
  congr 1
  refine Finset.sum_congr rfl (fun x _ => ?_)
  rw [MainAux.avgQ_const_mul]

namespace MainAux

set_option backward.isDefEq.respectTransparency false in
open scoped Matrix.Norms.Operator in
/-- The propagator `t ↦ e^{-itA}` is continuous. -/
private lemma continuous_evo {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℝ) :
    Continuous (fun t : ℝ => evo A t) := by
  unfold evo
  exact NormedSpace.exp_continuous.comp (by fun_prop)

/-- The one-body imbalance is continuous in time. -/
private lemma continuous_imbOne (J h : ℝ) (q : Cfg N) : Continuous (fun t => imbOne J h q t) := by
  have hU := continuous_evo (hq J h q)
  unfold imbOne
  exact continuous_const.mul (Complex.continuous_re.comp
    ((((continuous_const.matrix_mul hU).matrix_mul continuous_const).matrix_mul
      hU.matrix_conjTranspose).matrix_trace))

end MainAux

/-- The time-averaged LGT imbalance converges to `Ī = 2^{-N} ∑_q N⁻¹ ∑_α a_α(q)²`. -/
theorem lgtImbalance_timeAvg (J h : ℝ) (hJ : J ≠ 0) (hN : Even N) :
    Tendsto (fun T : ℝ => (1 / T) * ∫ t in (0 : ℝ)..T, lgtImbalance J h N t) atTop
      (𝓝 (avgQ (fun q : Cfg N => imbInf J h q))) := by
  have heq : (fun T : ℝ => (1 / T) * ∫ t in (0 : ℝ)..T, lgtImbalance J h N t) =
      fun T => (1 / 2 ^ N : ℝ) * ∑ q : Cfg N, ((1 / T) * ∫ t in (0 : ℝ)..T, imbOne J h q t) := by
    funext T
    simp_rw [lgtImbalance_eq J h hN, avgQ]
    rw [intervalIntegral.integral_const_mul, intervalIntegral.integral_finsetSum
      (fun q _ => (MainAux.continuous_imbOne J h q).intervalIntegrable 0 T), Finset.mul_sum,
      Finset.mul_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl (fun q _ => ?_)
    ring
  rw [heq]
  exact (tendsto_finsetSum _ (fun q _ => imbOne_timeAvg J h hJ q)).const_mul _

/-- **Persistent imbalance** (paper: Theorem `thm:imbalance`): for every `h ≠ 0`, under
Hypothesis (EC), the infinite-time imbalance of the LGT is bounded below by a positive constant
independent of the (even) volume. -/
theorem lgtImbalance_lower (J h : ℝ) (hJ : J ≠ 0) (hh : h ≠ 0) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ c : ℝ, 0 < c ∧ ∀ (N : ℕ), Even N → 2 ≤ N →
      ∀ ε > 0, ∀ᶠ T in atTop,
        c - ε ≤ (1 / T) * ∫ t in (0 : ℝ)..T, lgtImbalance J h N t := by
  obtain ⟨c, hc, hcN⟩ := avg_imbInf_uniform J h hJ hh hEC
  refine ⟨c, hc, fun N hN hN2 ε hε => ?_⟩
  have hge : c ≤ avgQ (fun q : Cfg N => imbInf J h q) := hcN N (by omega)
  filter_upwards [(lgtImbalance_timeAvg J h hJ hN).eventually
    (lt_mem_nhds (show avgQ (fun q : Cfg N => imbInf J h q) - ε <
      avgQ (fun q : Cfg N => imbInf J h q) by linarith))] with T hT
  linarith

end DFL
