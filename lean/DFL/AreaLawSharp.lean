import DFL.SlaterExact
import DFL.BinEntropy
import DFL.ECProof

/-!
# The area law with the constants of the paper (paper: Lemma `lem:annealed`, Theorem `thm:arealaw`(a),(b))

* `vonNeumann_slater_sharp` (paper: `thm:arealaw`(a)): `S(ρ_A) ≤ 2 ln 2 ∑_{x ∈ A, y ∉ A} |K_{xy}|`
  for every Slater determinant, from `S(ρ_A) = ∑ h(ν_i)`, `h(ν) ≤ 2 ln 2 √(ν(1-ν))` and the trace-norm
  bound.
* `annealed_K`, `annealed_cut` (paper: `lem:annealed`(i),(ii)): under (EC),
  `E sup_t |K_q(t)_{xy}| ≤ 2 ∑_{2k ≥ |x-y|} F(k)` and `E sup_t ∑_{x ∈ A, y ∉ A} |K_q(t)_{xy}| ≤ Θ(F)`,
  `Θ(F) = 2 ∑_{k ≥ 1} k (2k+1) F(k)`.
* `lgt_areaLaw_sharp` (paper: `thm:arealaw`(b)): `E_q sup_t S(ρ^q_A(t)) ≤ 2 ln 2 Θ(F) ≤ 8 ln 2 M_2(F)`;
  `lgt_areaLaw_sharp_of_ldt`: the same with (EC) replaced by the large-deviation input (`ec_of_ldt`).
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- `Θ(F) = 2 ∑_{k ≥ 1} k (2k+1) F(k)`. -/
def thetaF (F : ℕ → ℝ) : ℝ := 2 * ∑' k : ℕ, (k : ℝ) * (2 * k + 1) * F k

namespace AreaLawSharpAux

open AreaLawAux

/-! ### Summability -/

/-- Under (EC) the profile `F` is summable. -/
lemma summable_F {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) : Summable F :=
  Summable.of_nonneg_of_le hEC.nonneg
    (fun d => le_mul_of_one_le_left (hEC.nonneg d)
      (by nlinarith [(Nat.cast_nonneg d : (0 : ℝ) ≤ d)])) hEC.summable

/-- A truncation `k ↦ [P k] F(k)` of a nonnegative summable profile is summable. -/
lemma summable_ite {F : ℕ → ℝ} (hF0 : ∀ d, 0 ≤ F d) (hF : Summable F) (P : ℕ → Prop)
    [DecidablePred P] : Summable (fun k => if P k then F k else 0) := by
  refine Summable.of_nonneg_of_le (fun k => ?_) (fun k => ?_) hF
  · split_ifs
    · exact hF0 k
    · exact le_rfl
  · split_ifs
    · exact le_rfl
    · exact hF0 k

/-- Under (EC), `k ↦ k (2k+1) F(k)` is summable. -/
lemma summable_theta {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) :
    Summable (fun k : ℕ => (k : ℝ) * (2 * k + 1) * F k) := by
  refine Summable.of_nonneg_of_le (fun k => ?_) (fun k => ?_) (hEC.summable.mul_left 2)
  · exact mul_nonneg (mul_nonneg (Nat.cast_nonneg k) (by positivity)) (hEC.nonneg k)
  · have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
    have hF := hEC.nonneg k
    nlinarith [mul_nonneg hk hF]

/-! ### Counting -/

/-- Each value of `z ↦ max(|x-z|, |y-z|)` is attained at most twice, and only at values `k` with
`|x-y| ≤ 2k`: `∑_z f(max(|x-z|, |y-z|)) ≤ 2 ∑_{k : |x-y| ≤ 2k} f(k)` (here `x ≤ y`). -/
lemma sum_max_dist_le (f : ℕ → ℝ) (hf0 : ∀ n, 0 ≤ f n) (hf : Summable f) (x y : ℕ)
    (hxy : x ≤ y) :
    ∑ z : Fin N, f (max (Nat.dist x z.val) (Nat.dist y z.val)) ≤
      2 * ∑' k, (if Nat.dist x y ≤ 2 * k then f k else 0) := by
  set g : ℕ → ℝ := fun k => if Nat.dist x y ≤ 2 * k then f k else 0 with hg
  have hg0 : ∀ n, 0 ≤ g n := fun n => by
    simp only [hg]
    split_ifs
    · exact hf0 n
    · exact le_rfl
  have hgs : Summable g := summable_ite hf0 hf _
  rw [← Finset.sum_filter_add_sum_filter_not univ (fun z : Fin N => 2 * z.val ≤ x + y)]
  have h1 : ∑ z ∈ univ.filter (fun z : Fin N => 2 * z.val ≤ x + y),
      f (max (Nat.dist x z.val) (Nat.dist y z.val)) ≤ ∑' k, g k := by
    have e : ∀ z ∈ univ.filter (fun z : Fin N => 2 * z.val ≤ x + y),
        f (max (Nat.dist x z.val) (Nat.dist y z.val)) = g (y - z.val) := by
      intro z hz
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hz
      have hm : max (Nat.dist x z.val) (Nat.dist y z.val) = y - z.val := by
        rw [max_eq_right (by unfold Nat.dist; omega)]
        unfold Nat.dist; omega
      have hd : Nat.dist x y ≤ 2 * (y - z.val) := by unfold Nat.dist; omega
      simp only [hg, hm, ite_eq_left hd]
    rw [Finset.sum_congr rfl e]
    refine sum_le_tsum_injOn _ (fun z : Fin N => y - z.val) ?_ g hg0 hgs
    intro a ha b hb hab
    simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq] at ha hb
    exact Fin.ext (by simp only at hab; omega)
  have h2 : ∑ z ∈ univ.filter (fun z : Fin N => ¬ 2 * z.val ≤ x + y),
      f (max (Nat.dist x z.val) (Nat.dist y z.val)) ≤ ∑' k, g k := by
    have e : ∀ z ∈ univ.filter (fun z : Fin N => ¬ 2 * z.val ≤ x + y),
        f (max (Nat.dist x z.val) (Nat.dist y z.val)) = g (z.val - x) := by
      intro z hz
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hz
      have hm : max (Nat.dist x z.val) (Nat.dist y z.val) = z.val - x := by
        rw [max_eq_left (by unfold Nat.dist; omega)]
        unfold Nat.dist; omega
      have hd : Nat.dist x y ≤ 2 * (z.val - x) := by unfold Nat.dist; omega
      simp only [hg, hm, ite_eq_left hd]
    rw [Finset.sum_congr rfl e]
    refine sum_le_tsum_injOn _ (fun z : Fin N => z.val - x) ?_ g hg0 hgs
    intro a ha b hb hab
    simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq] at ha hb
    exact Fin.ext (by simp only at hab; omega)
  linarith

/-- `∑_{s < n} (s + 1) = n (n+1) / 2`, in the form `∑_{s < 2k} (s+1) = k (2k+1)`. -/
lemma sum_range_succ_two_mul (k : ℕ) : ∑ s ∈ range (2 * k), (s + 1) = k * (2 * k + 1) := by
  have h1 := Finset.sum_range_succ' (fun i => i) (2 * k)
  have h2 := Finset.sum_range_id_mul_two (2 * k + 1)
  rw [h1, add_zero, Nat.add_sub_cancel] at h2
  have h3 : (∑ s ∈ range (2 * k), (s + 1)) * 2 = (k * (2 * k + 1)) * 2 := by rw [h2]; ring
  exact Nat.eq_of_mul_eq_mul_right (by norm_num) h3

/-- At most `k (2k+1)` pairs `x < ℓ ≤ y` have `y - x ≤ 2k`. -/
lemma card_cut_pairs_le (ℓ k : ℕ) :
    ((segment (N := N) ℓ ×ˢ (univ \ segment (N := N) ℓ)).filter
      (fun p : Fin N × Fin N => Nat.dist p.1.val p.2.val ≤ 2 * k)).card ≤ k * (2 * k + 1) := by
  have hT : ((range (2 * k)).sigma (fun s => range (s + 1))).card = k * (2 * k + 1) := by
    rw [Finset.card_sigma]
    simp only [Finset.card_range]
    exact sum_range_succ_two_mul k
  rw [← hT]
  refine Finset.card_le_card_of_injOn
    (fun p : Fin N × Fin N => (⟨p.2.val - p.1.val - 1, ℓ - 1 - p.1.val⟩ : Σ _ : ℕ, ℕ)) ?_ ?_
  · intro p hp
    simp only [Finset.coe_filter, Finset.mem_product, segment, Finset.mem_sdiff,
      Finset.mem_filter, Finset.mem_univ, true_and, not_lt, Set.mem_ofPred_eq] at hp
    obtain ⟨⟨h1, h2⟩, h3⟩ := hp
    unfold Nat.dist at h3
    simp only [Finset.coe_sigma, Set.mem_sigma_iff, Finset.coe_range, Set.mem_Iio]
    constructor <;> omega
  · intro p hp p' hp' hpp'
    simp only [Finset.coe_filter, Finset.mem_product, segment, Finset.mem_sdiff,
      Finset.mem_filter, Finset.mem_univ, true_and, not_lt, Set.mem_ofPred_eq] at hp hp'
    obtain ⟨⟨h1, h2⟩, h3⟩ := hp
    obtain ⟨⟨h1', h2'⟩, h3'⟩ := hp'
    unfold Nat.dist at h3 h3'
    simp only [Sigma.mk.inj_iff, heq_eq_eq] at hpp'
    obtain ⟨e1, e2⟩ := hpp'
    exact Prod.ext (Fin.ext (by omega)) (Fin.ext (by omega))

/-- `∑_{x<ℓ≤y} [|x-y| ≤ 2k] F(k) ≤ k (2k+1) F(k)`. -/
lemma sum_cut_ite_le (F : ℕ → ℝ) (hF0 : ∀ d, 0 ≤ F d) (ℓ k : ℕ) :
    ∑ x ∈ segment (N := N) ℓ, ∑ y ∈ univ \ segment (N := N) ℓ,
        (if Nat.dist x.val y.val ≤ 2 * k then F k else 0) ≤ (k : ℝ) * (2 * k + 1) * F k := by
  have e : ∑ x ∈ segment (N := N) ℓ, ∑ y ∈ univ \ segment (N := N) ℓ,
        (if Nat.dist x.val y.val ≤ 2 * k then F k else 0) =
      ∑ p ∈ segment (N := N) ℓ ×ˢ (univ \ segment (N := N) ℓ),
        (if Nat.dist p.1.val p.2.val ≤ 2 * k then F k else 0) :=
    (Finset.sum_product' _ _
      (fun x y : Fin N => if Nat.dist x.val y.val ≤ 2 * k then F k else 0)).symm
  rw [e, ← Finset.sum_filter, Finset.sum_const, nsmul_eq_mul]
  have hc : ((((segment (N := N) ℓ ×ˢ (univ \ segment (N := N) ℓ)).filter
      (fun p : Fin N × Fin N => Nat.dist p.1.val p.2.val ≤ 2 * k)).card : ℕ) : ℝ) ≤
      (k : ℝ) * (2 * k + 1) := by exact_mod_cast card_cut_pairs_le ℓ k
  exact mul_le_mul_of_nonneg_right hc (hF0 k)

/-- `∑_{x<ℓ≤y} 2 ∑_{k : |x-y| ≤ 2k} F(k) ≤ Θ(F)`. -/
lemma cut_count_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (ℓ : ℕ) :
    ∑ x ∈ segment (N := N) ℓ, ∑ y ∈ univ \ segment (N := N) ℓ,
        2 * ∑' k : ℕ, (if Nat.dist x.val y.val ≤ 2 * k then F k else 0) ≤ thetaF F := by
  have hs : ∀ x y : Fin N,
      Summable (fun k : ℕ => if Nat.dist x.val y.val ≤ 2 * k then F k else 0) :=
    fun x y => summable_ite hEC.nonneg (summable_F hEC) _
  have hsum : ∑ x ∈ segment (N := N) ℓ, ∑ y ∈ univ \ segment (N := N) ℓ,
        2 * ∑' k : ℕ, (if Nat.dist x.val y.val ≤ 2 * k then F k else 0) =
      2 * ∑' k : ℕ, ∑ x ∈ segment (N := N) ℓ, ∑ y ∈ univ \ segment (N := N) ℓ,
        (if Nat.dist x.val y.val ≤ 2 * k then F k else 0) := by
    rw [Summable.tsum_finsetSum (fun x _ => summable_sum fun y _ => hs x y), Finset.mul_sum]
    refine Finset.sum_congr rfl fun x _ => ?_
    rw [Summable.tsum_finsetSum (fun y _ => hs x y), Finset.mul_sum]
  rw [hsum, thetaF]
  exact mul_le_mul_of_nonneg_left (Summable.tsum_le_tsum (fun k => sum_cut_ite_le F hEC.nonneg ℓ k)
    (summable_sum fun x _ => summable_sum fun y _ => hs x y) (summable_theta hEC)) (by norm_num)

/-! ### Annealed bounds on the correlator -/

/-- `|K_{xy}(t)| ≤ ∑_z Q(x,z) Q(y,z)`, uniformly in time. -/
lemma norm_corrMat_le_corrQ (J h : ℝ) (q : Cfg N) (S : Finset (Fin N)) (t : ℝ) (x y : Fin N) :
    ‖corrMat (hq J h q) S t x y‖ ≤
      ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) x z *
        corrQ (hq J h q) (hq_isHermitian J h q) y z := by
  refine (norm_corrMat_le _ S t x y).trans (Finset.sum_le_sum fun z _ => ?_)
  exact mul_le_mul (norm_evo_le_corrQ _ _ t x z) (norm_evo_le_corrQ _ _ t y z) (norm_nonneg _)
    (corrQ_nonneg _ _ x z)

/-- `E[Q(x,z) Q(y,z)] ≤ F(max(|x-z|, |y-z|))`, from `0 ≤ Q ≤ 1` and (EC). -/
lemma avgQ_corrQ_mul_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (x y z : Fin N) :
    avgQ (fun q => corrQ (hq J h q) (hq_isHermitian J h q) x z *
        corrQ (hq J h q) (hq_isHermitian J h q) y z) ≤
      F (max (Nat.dist x.val z.val) (Nat.dist y.val z.val)) := by
  rcases le_total (Nat.dist x.val z.val) (Nat.dist y.val z.val) with hle | hle
  · rw [max_eq_right hle]
    refine (avgQ_mono' fun q => ?_).trans (hEC.bound N y z)
    exact mul_le_of_le_one_left (corrQ_nonneg _ _ y z) (ECProofAux.corrQ_le_one _ _ x z)
  · rw [max_eq_left hle]
    refine (avgQ_mono' fun q => ?_).trans (hEC.bound N x z)
    exact mul_le_of_le_one_right (corrQ_nonneg _ _ x z) (ECProofAux.corrQ_le_one _ _ y z)

/-- `E ∑_z Q(x,z) Q(y,z) ≤ 2 ∑_{k : |x-y| ≤ 2k} F(k)`. -/
lemma avgQ_sum_corrQ_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (x y : Fin N) :
    avgQ (fun q => ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) x z *
        corrQ (hq J h q) (hq_isHermitian J h q) y z) ≤
      2 * ∑' k : ℕ, (if Nat.dist x.val y.val ≤ 2 * k then F k else 0) := by
  rw [avgQ_sum']
  refine (Finset.sum_le_sum fun z _ => avgQ_corrQ_mul_le hEC x y z).trans ?_
  rcases le_total x.val y.val with hxy | hxy
  · exact sum_max_dist_le F hEC.nonneg (summable_F hEC) x.val y.val hxy
  · calc ∑ z : Fin N, F (max (Nat.dist x.val z.val) (Nat.dist y.val z.val))
        = ∑ z : Fin N, F (max (Nat.dist y.val z.val) (Nat.dist x.val z.val)) :=
          Finset.sum_congr rfl fun z _ => by rw [max_comm]
      _ ≤ 2 * ∑' k : ℕ, (if Nat.dist y.val x.val ≤ 2 * k then F k else 0) :=
          sum_max_dist_le F hEC.nonneg (summable_F hEC) y.val x.val hxy
      _ = 2 * ∑' k : ℕ, (if Nat.dist x.val y.val ≤ 2 * k then F k else 0) := by
          rw [Nat.dist_comm]

/-- `E ∑_{x<ℓ≤y} ∑_z Q(x,z) Q(y,z) ≤ Θ(F)`. -/
lemma avgQ_cut_corrQ_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (ℓ : ℕ) :
    avgQ (fun q : Cfg N => ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ,
        ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) x z *
          corrQ (hq J h q) (hq_isHermitian J h q) y z) ≤ thetaF F := by
  calc _ = ∑ x ∈ segment (N := N) ℓ, ∑ y ∈ univ \ segment (N := N) ℓ,
        avgQ (fun q : Cfg N => ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) x z *
          corrQ (hq J h q) (hq_isHermitian J h q) y z) := by
        rw [avgQ_sum']
        exact Finset.sum_congr rfl fun x _ => avgQ_sum' _ _
    _ ≤ ∑ x ∈ segment (N := N) ℓ, ∑ y ∈ univ \ segment (N := N) ℓ,
        2 * ∑' k : ℕ, (if Nat.dist x.val y.val ≤ 2 * k then F k else 0) :=
        Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ => avgQ_sum_corrQ_le hEC x y
    _ ≤ thetaF F := cut_count_le hEC ℓ

/-- `0 ≤ 2 ln 2`. -/
lemma two_log_two_nonneg : (0 : ℝ) ≤ 2 * Real.log 2 :=
  mul_nonneg (by norm_num) (Real.log_nonneg (by norm_num))

end AreaLawSharpAux

open AreaLawAux AreaLawSharpAux

/-- **Area law for Slater determinants with the constant `2 ln 2`** (paper: Theorem
`thm:arealaw`(a)). -/
theorem vonNeumann_slater_sharp (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (rdm ℓ hℓ (evo (dGamma A) t *ᵥ fockVec S)) (rdm_isHermitian ℓ hℓ _) ≤
      2 * Real.log 2 * ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ‖corrMat A S t x y‖ := by
  rw [vonNeumann_rdm_slater A hA S t ℓ hℓ]
  calc ∑ i, Real.binEntropy ((corrA_isHermitian A hA S t ℓ hℓ).eigenvalues i)
      ≤ ∑ i, 2 * Real.log 2 * √((corrA_isHermitian A hA S t ℓ hℓ).eigenvalues i *
          (1 - (corrA_isHermitian A hA S t ℓ hℓ).eigenvalues i)) :=
        Finset.sum_le_sum fun i _ => binEntropy_le_two_log_two_sqrt _
          (eig_corrA_mem A hA S t ℓ hℓ i).1 (eig_corrA_mem A hA S t ℓ hℓ i).2
    _ = 2 * Real.log 2 * ∑ i, √((corrA_isHermitian A hA S t ℓ hℓ).eigenvalues i *
          (1 - (corrA_isHermitian A hA S t ℓ hℓ).eigenvalues i)) := by
        rw [Finset.mul_sum]
    _ ≤ 2 * Real.log 2 * ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ‖corrMat A S t x y‖ :=
        mul_le_mul_of_nonneg_left (sum_sqrt_eig_le A hA S t ℓ hℓ) two_log_two_nonneg

/-- **Annealed correlation bound** (paper: Lemma `lem:annealed`(i)):
`E sup_t |K_q(t)_{xy}| ≤ 2 ∑_{k : 2k ≥ |x-y|} F(k)`. -/
theorem annealed_K (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (x y : Fin N) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ, ‖corrMat (hq J h q) S t x y‖) ≤
      2 * ∑' k : ℕ, (if Nat.dist x.val y.val ≤ 2 * k then F k else 0) := by
  refine (avgQ_mono' fun q => ?_).trans (avgQ_sum_corrQ_le hEC x y)
  exact ciSup_le fun t => norm_corrMat_le_corrQ J h q S t x y

/-- **Annealed cut bound** (paper: Lemma `lem:annealed`(ii)):
`E sup_t ∑_{x ∈ A, y ∉ A} |K_q(t)_{xy}| ≤ Θ(F)` for a left segment `A`. -/
theorem annealed_cut (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (ℓ : ℕ) :
    avgQ (fun q : Cfg N =>
      ⨆ t : ℝ, ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ‖corrMat (hq J h q) S t x y‖) ≤
      thetaF F := by
  refine (avgQ_mono' fun q => ?_).trans (avgQ_cut_corrQ_le hEC ℓ)
  exact ciSup_le fun t => Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ =>
    norm_corrMat_le_corrQ J h q S t x y

/-- **Sector-averaged area law with the constant of the paper** (paper: Theorem `thm:arealaw`(b)):
`E_q sup_t S(ρ^q_A(t)) ≤ 2 ln 2 Θ(F)`, uniformly in `N`, `ℓ` and the Fock state. -/
theorem lgt_areaLaw_sharp (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (ℓ : ℕ) (hℓ : ℓ ≤ N)
    (S : Finset (Fin N)) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ,
      vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ fockVec S)) (rdm_isHermitian ℓ hℓ _)) ≤
      2 * Real.log 2 * thetaF F := by
  have hpt : ∀ q : Cfg N, (⨆ t : ℝ,
      vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ fockVec S)) (rdm_isHermitian ℓ hℓ _)) ≤
      2 * Real.log 2 * ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ,
        ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) x z *
          corrQ (hq J h q) (hq_isHermitian J h q) y z := by
    intro q
    refine ciSup_le fun t =>
      (vonNeumann_slater_sharp (hq J h q) (hq_isSymm J h q) S t ℓ hℓ).trans ?_
    exact mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ =>
      norm_corrMat_le_corrQ J h q S t x y) two_log_two_nonneg
  refine (avgQ_mono' hpt).trans ?_
  rw [avgQ_const_mul']
  exact mul_le_mul_of_nonneg_left (avgQ_cut_corrQ_le hEC ℓ) two_log_two_nonneg

/-- `Θ(F) ≤ 4 M_2(F)`, `M_2(F) = ∑_{d ≥ 0} (d+1)² F(d)`. -/
theorem thetaF_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) :
    thetaF F ≤ 4 * ∑' d : ℕ, ((d : ℝ) + 1) ^ 2 * F d := by
  rw [thetaF, ← tsum_mul_left, ← tsum_mul_left]
  refine Summable.tsum_le_tsum (fun k => ?_) ((summable_theta hEC).mul_left 2)
    (hEC.summable.mul_left 4)
  have hk : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  have hF := hEC.nonneg k
  nlinarith [mul_nonneg hk hF]

/-- **Sharp area law from uniform large deviations** (paper: Theorems `thm:EC`, `thm:arealaw`(b)):
`E_q sup_t S(ρ^q_A(t)) ≤ 2 ln 2 Θ(F)` with `F(d) = min(1, C e^{-ηd})`, uniformly in `N`, `ℓ` and the
Fock state. -/
theorem lgt_areaLaw_sharp_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ ∀ (N : ℕ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (S : Finset (Fin N)),
      avgQ (fun q : Cfg N => ⨆ t : ℝ,
        vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ fockVec S)) (rdm_isHermitian ℓ hℓ _)) ≤
        2 * Real.log 2 * thetaF (fun d => min 1 (C * Real.exp (-η * d))) := by
  obtain ⟨C, η, hη, hEC⟩ := ec_of_ldt J h hJ hL
  exact ⟨C, η, hη, fun N ℓ hℓ S => lgt_areaLaw_sharp J h hEC ℓ hℓ S⟩

end DFL
