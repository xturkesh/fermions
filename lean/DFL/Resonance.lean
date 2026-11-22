import DFL.Transfer

/-!
# Resonant energies of a box (paper: Lemma `lem:badset`, Eq. `eq:Bad`)

For a box `B = [a, a+L-1]`, a site `x = a + i ∈ B`, a choice of ends (`endL`: the left end `a`,
`endR`: the right end `a+L-1`) and a threshold `θ = e^{-κn}`, the resonant set is
`{E ∈ I_λ : |p_e(E)| ≥ θ |P_B(E)| for some selected end e}`, with `P_B = dirL V E a (L+1)`,
`p_b = dirL V E a (i+1)` (right end) and `p_a = dirR V E (a+L-1) (L-i)` (left end).

* `decay_of_not_res`: outside the resonant set every solution decays across the box:
  `|φ(x)| ≤ θ (|φ(a-1)| + |φ(a+L)|)` (paper: Lemma `lem:badset`(a) with Eq. `eq:poisson`).
* `exists_resEnd`, `card_resEnds`: every resonant energy `E` lies in an interval `[z, E]` of resonant
  energies whose left end `z` belongs to an explicit set of at most `4L + 2` energies, determined by
  the potential on the box (paper: Lemma `lem:badset`(b)).
-/

open Matrix Finset Polynomial

noncomputable section

namespace DFL

/-- `I_λ = [-2-|λ|, 2+|λ|]`, which contains the spectrum of `Δ + λ diag(q)`. -/
def Ilam (lam : ℝ) : Set ℝ := Set.Icc (-2 - |lam|) (2 + |lam|)

/-- The resonant set of `(x = a + i, B = [a, a+L-1], ends)` at threshold `θ` (paper: Eq. `eq:Bad`). -/
def resSet (V : ℤ → ℝ) (lam : ℝ) (a : ℤ) (L i : ℕ) (endL endR : Bool) (θ : ℝ) : Set ℝ :=
  {E | E ∈ Ilam lam ∧
    ((endR = true ∧ θ * |dirL V E a (L + 1)| ≤ |dirL V E a (i + 1)|) ∨
      (endL = true ∧ θ * |dirL V E a (L + 1)| ≤ |dirR V E (a + L - 1) (L - i)|))}

/-- The candidate left endpoints of the components of the resonant set: the two ends of `I_λ` and the
real zeros of `p_b² - θ² P_B²` and of `p_a² - θ² P_B²`. -/
def resEnds (V : ℤ → ℝ) (lam : ℝ) (a : ℤ) (L i : ℕ) (θ : ℝ) : Finset ℝ :=
  {-2 - |lam|, 2 + |lam|} ∪
    ((dirLpoly V a (i + 1)) ^ 2 - C (θ ^ 2) * (dirLpoly V a (L + 1)) ^ 2).roots.toFinset ∪
    ((dirRpoly V (a + L - 1) (L - i)) ^ 2 - C (θ ^ 2) * (dirLpoly V a (L + 1)) ^ 2).roots.toFinset

namespace ResonanceAux

/-- `p² - c P²` has at most `2L` distinct roots when `deg p ≤ L` and `deg P ≤ L`. -/
theorem card_roots_toFinset_le (p P : ℝ[X]) (c : ℝ) (L : ℕ) (hp : p.natDegree ≤ L)
    (hP : P.natDegree ≤ L) : (p ^ 2 - C c * P ^ 2).roots.toFinset.card ≤ 2 * L := by
  calc (p ^ 2 - C c * P ^ 2).roots.toFinset.card
      ≤ Multiset.card (p ^ 2 - C c * P ^ 2).roots := Multiset.toFinset_card_le _
    _ ≤ (p ^ 2 - C c * P ^ 2).natDegree := Polynomial.card_roots' _
    _ ≤ max (p ^ 2).natDegree (C c * P ^ 2).natDegree := natDegree_sub_le _ _
    _ ≤ 2 * L := by
      apply max_le
      · calc (p ^ 2).natDegree ≤ 2 * p.natDegree := natDegree_pow_le
          _ ≤ 2 * L := by omega
      · calc (C c * P ^ 2).natDegree ≤ (P ^ 2).natDegree := natDegree_C_mul_le _ _
          _ ≤ 2 * P.natDegree := natDegree_pow_le
          _ ≤ 2 * L := by omega

/-- `p² - c P² ≠ 0` when `P` is monic of degree `L`, `deg p < L` and `c ≠ 0` (its coefficient of
`X^{2L}` is `-c`). -/
theorem sq_sub_C_mul_sq_ne_zero (p P : ℝ[X]) (c : ℝ) (L : ℕ) (hP : P.Monic)
    (hPdeg : P.natDegree = L) (hp : p.natDegree < L) (hc : c ≠ 0) :
    p ^ 2 - C c * P ^ 2 ≠ 0 := by
  intro h
  have h1 : (p ^ 2 - C c * P ^ 2).coeff (2 * L) = 0 := by rw [h, coeff_zero]
  rw [coeff_sub, coeff_C_mul] at h1
  have h2 : (p ^ 2).coeff (2 * L) = 0 := by
    apply coeff_eq_zero_of_natDegree_lt
    calc (p ^ 2).natDegree ≤ 2 * p.natDegree := natDegree_pow_le
      _ < 2 * L := by omega
  have h3 : (P ^ 2).coeff (2 * L) = 1 := by
    have hd : (P ^ 2).natDegree = 2 * L := by rw [hP.natDegree_pow, hPdeg]
    rw [← hd]
    exact (hP.pow 2).coeff_natDegree
  rw [h2, h3, mul_one, zero_sub, neg_eq_zero] at h1
  exact hc h1

/-- `dirL V E a k` is continuous in the energy `E`. -/
theorem continuous_dirL (V : ℤ → ℝ) (a : ℤ) (k : ℕ) : Continuous fun E : ℝ => dirL V E a k := by
  have : (fun E : ℝ => dirL V E a k) = fun E => (dirLpoly V a k).eval E :=
    funext fun E => (eval_dirLpoly V E a k).symm
  rw [this]
  exact (dirLpoly V a k).continuous

/-- `dirR V E b k` is continuous in the energy `E`. -/
theorem continuous_dirR (V : ℤ → ℝ) (b : ℤ) (k : ℕ) : Continuous fun E : ℝ => dirR V E b k := by
  have : (fun E : ℝ => dirR V E b k) = fun E => (dirRpoly V b k).eval E :=
    funext fun E => (eval_dirRpoly V E b k).symm
  rw [this]
  exact (dirRpoly V b k).continuous

/-- The degree of `p_a = dirRpoly V b (L - i)` is `L - i - 1`. -/
theorem natDegree_dirRpoly_sub (V : ℤ → ℝ) (b : ℤ) (L i : ℕ) (hi : i < L) :
    (dirRpoly V b (L - i)).natDegree = L - i - 1 := by
  obtain ⟨k, hk⟩ : ∃ k, L - i = k + 1 := ⟨L - i - 1, by omega⟩
  rw [hk, natDegree_dirRpoly]
  omega

/-- Rewriting the sites `a + L - 1 - j` of the box as `a + (L - 1 - j)`. -/
theorem box_site (a : ℤ) (L j : ℕ) (hj : j + 1 ≤ L) :
    a + (L : ℤ) - 1 - (j : ℤ) = a + ((L - 1 - j : ℕ) : ℤ) := by
  omega

end ResonanceAux

/-- **Decay outside the resonant set**: if at least one end is selected, `E ∈ I_λ` is not resonant,
`φ` solves the equation on the box, and `φ` vanishes beyond every unselected end, then
`|φ(x)| ≤ θ (|φ(a-1)| + |φ(a+L)|)`. -/
theorem decay_of_not_res (V : ℤ → ℝ) (lam : ℝ) (a : ℤ) (L i : ℕ) (endL endR : Bool) (θ : ℝ)
    (hi : i < L) (hend : endL = true ∨ endR = true) (E : ℝ) (hE : E ∈ Ilam lam)
    (hres : E ∉ resSet V lam a L i endL endR θ) (φ : ℤ → ℝ)
    (hφ : ∀ z, a ≤ z → z ≤ a + L - 1 → φ (z + 1) + φ (z - 1) + V z * φ z = E * φ z)
    (hL : endL = false → φ (a - 1) = 0) (hR : endR = false → φ (a + L) = 0) :
    |φ (a + i)| ≤ θ * (|φ (a - 1)| + |φ (a + L)|) := by
  have hpois := poisson V E φ a L i hi hφ
  set P := dirL V E a (L + 1)
  set pb := dirL V E a (i + 1)
  set pa := dirR V E (a + L - 1) (L - i)
  have hnot : ¬((endR = true ∧ θ * |P| ≤ |pb|) ∨ (endL = true ∧ θ * |P| ≤ |pa|)) :=
    fun h => hres ⟨hE, h⟩
  push Not at hnot
  obtain ⟨hR', hL'⟩ := hnot
  have hpos : 0 < θ * |P| := by
    rcases hend with h | h
    · exact lt_of_le_of_lt (abs_nonneg _) (hL' h)
    · exact lt_of_le_of_lt (abs_nonneg _) (hR' h)
  have hPpos : 0 < |P| := by
    rcases (abs_nonneg P).lt_or_eq with h | h
    · exact h
    · rw [← h, mul_zero] at hpos
      exact absurd hpos (lt_irrefl 0)
  have ha : |pa| * |φ (a - 1)| ≤ θ * |P| * |φ (a - 1)| := by
    cases hl : endL
    · rw [hL hl, abs_zero, mul_zero, mul_zero]
    · exact mul_le_mul_of_nonneg_right (hL' hl).le (abs_nonneg _)
  have hb : |pb| * |φ (a + L)| ≤ θ * |P| * |φ (a + L)| := by
    cases hr : endR
    · rw [hR hr, abs_zero, mul_zero, mul_zero]
    · exact mul_le_mul_of_nonneg_right (hR' hr).le (abs_nonneg _)
  have key : |P| * |φ (a + i)| ≤ |P| * (θ * (|φ (a - 1)| + |φ (a + L)|)) := by
    calc |P| * |φ (a + i)| = |P * φ (a + i)| := (abs_mul _ _).symm
      _ = |pa * φ (a - 1) + pb * φ (a + L)| := by rw [hpois]
      _ ≤ |pa * φ (a - 1)| + |pb * φ (a + L)| := abs_add_le _ _
      _ = |pa| * |φ (a - 1)| + |pb| * |φ (a + L)| := by rw [abs_mul, abs_mul]
      _ ≤ θ * |P| * |φ (a - 1)| + θ * |P| * |φ (a + L)| := add_le_add ha hb
      _ = |P| * (θ * (|φ (a - 1)| + |φ (a + L)|)) := by ring
  exact le_of_mul_le_mul_left key hPpos

/-- The resonant set depends only on the potential on the box. -/
theorem resSet_congr (V V' : ℤ → ℝ) (lam : ℝ) (a : ℤ) (L i : ℕ) (endL endR : Bool) (θ : ℝ)
    (hi : i < L) (h : ∀ j : ℕ, j < L → V (a + j) = V' (a + j)) :
    resSet V lam a L i endL endR θ = resSet V' lam a L i endL endR θ := by
  have hP : ∀ E, dirL V E a (L + 1) = dirL V' E a (L + 1) := fun E =>
    dirL_congr V V' E a (L + 1) fun j hj => h j (by omega)
  have hpb : ∀ E, dirL V E a (i + 1) = dirL V' E a (i + 1) := fun E =>
    dirL_congr V V' E a (i + 1) fun j hj => h j (by omega)
  have hpa : ∀ E, dirR V E (a + L - 1) (L - i) = dirR V' E (a + L - 1) (L - i) := fun E =>
    dirR_congr V V' E (a + L - 1) (L - i) fun j hj => by
      rw [ResonanceAux.box_site a L j (by omega)]
      exact h (L - 1 - j) (by omega)
  ext E
  simp only [resSet, Set.mem_ofPred_eq, hP, hpb, hpa]

/-- The candidate endpoints depend only on the potential on the box. -/
theorem resEnds_congr (V V' : ℤ → ℝ) (lam : ℝ) (a : ℤ) (L i : ℕ) (θ : ℝ) (hi : i < L)
    (h : ∀ j : ℕ, j < L → V (a + j) = V' (a + j)) :
    resEnds V lam a L i θ = resEnds V' lam a L i θ := by
  have hP : dirLpoly V a (L + 1) = dirLpoly V' a (L + 1) :=
    dirLpoly_congr V V' a (L + 1) fun j hj => h j (by omega)
  have hpb : dirLpoly V a (i + 1) = dirLpoly V' a (i + 1) :=
    dirLpoly_congr V V' a (i + 1) fun j hj => h j (by omega)
  have hpa : dirRpoly V (a + L - 1) (L - i) = dirRpoly V' (a + L - 1) (L - i) :=
    dirRpoly_congr V V' (a + L - 1) (L - i) fun j hj => by
      rw [ResonanceAux.box_site a L j (by omega)]
      exact h (L - 1 - j) (by omega)
  unfold resEnds
  rw [hP, hpb, hpa]

/-- There are at most `4L + 2` candidate endpoints. -/
theorem card_resEnds (V : ℤ → ℝ) (lam : ℝ) (a : ℤ) (L i : ℕ) (θ : ℝ) (hi : i < L) :
    (resEnds V lam a L i θ).card ≤ 4 * L + 2 := by
  have hP : (dirLpoly V a (L + 1)).natDegree ≤ L := (natDegree_dirLpoly V a L).le
  have hpb : (dirLpoly V a (i + 1)).natDegree ≤ L := by
    rw [natDegree_dirLpoly]
    omega
  have hpa : (dirRpoly V (a + L - 1) (L - i)).natDegree ≤ L := by
    rw [ResonanceAux.natDegree_dirRpoly_sub V _ L i hi]
    omega
  have h1 := ResonanceAux.card_roots_toFinset_le _ _ (θ ^ 2) L hpb hP
  have h2 := ResonanceAux.card_roots_toFinset_le _ _ (θ ^ 2) L hpa hP
  have h0 : ({-2 - |lam|, 2 + |lam|} : Finset ℝ).card ≤ 2 := Finset.card_le_two
  have h3 := card_union_le ({-2 - |lam|, 2 + |lam|} : Finset ℝ)
    ((dirLpoly V a (i + 1)) ^ 2 - C (θ ^ 2) * (dirLpoly V a (L + 1)) ^ 2).roots.toFinset
  have h4 := card_union_le (({-2 - |lam|, 2 + |lam|} : Finset ℝ) ∪
    ((dirLpoly V a (i + 1)) ^ 2 - C (θ ^ 2) * (dirLpoly V a (L + 1)) ^ 2).roots.toFinset)
    ((dirRpoly V (a + L - 1) (L - i)) ^ 2 - C (θ ^ 2) * (dirLpoly V a (L + 1)) ^ 2).roots.toFinset
  unfold resEnds
  omega

/-- **Components of the resonant set** (paper: Lemma `lem:badset`(b)): every resonant energy `E`
lies in an interval `[z, E]` of resonant energies with `z` a candidate endpoint. -/
theorem exists_resEnd (V : ℤ → ℝ) (lam : ℝ) (a : ℤ) (L i : ℕ) (endL endR : Bool) (θ : ℝ)
    (hθ : 0 < θ) (hi : i < L) (E : ℝ) (hE : E ∈ resSet V lam a L i endL endR θ) :
    ∃ z ∈ resEnds V lam a L i θ, z ≤ E ∧ Set.Icc z E ⊆ resSet V lam a L i endL endR θ := by
  classical
  set R := resSet V lam a L i endL endR θ
  have cP := ResonanceAux.continuous_dirL V a (L + 1)
  have cpb := ResonanceAux.continuous_dirL V a (i + 1)
  have cpa := ResonanceAux.continuous_dirR V (a + L - 1) (L - i)
  -- the resonant set is closed
  have hRclosed : IsClosed R := by
    have h1 : IsClosed ({t : ℝ | endR = true} ∩
        {t : ℝ | θ * |dirL V t a (L + 1)| ≤ |dirL V t a (i + 1)|}) :=
      isClosed_const.inter (isClosed_le (continuous_const.mul cP.abs) cpb.abs)
    have h2 : IsClosed ({t : ℝ | endL = true} ∩
        {t : ℝ | θ * |dirL V t a (L + 1)| ≤ |dirR V t (a + L - 1) (L - i)|}) :=
      isClosed_const.inter (isClosed_le (continuous_const.mul cP.abs) cpa.abs)
    exact (isClosed_Icc (a := -2 - |lam|) (b := 2 + |lam|)).inter (h1.union h2)
  -- the admissible left endpoints `s` (with `[s, E] ⊆ R`) and their infimum `z`
  set S : Set ℝ := {s | -2 - |lam| ≤ s ∧ s ≤ E ∧ Set.Icc s E ⊆ R}
  have hES : E ∈ S :=
    ⟨hE.1.1, le_rfl, by rw [Set.Icc_self]; exact Set.singleton_subset_iff.mpr hE⟩
  have hSne : S.Nonempty := ⟨E, hES⟩
  have hSbdd : BddBelow S := ⟨-2 - |lam|, fun s hs => hs.1⟩
  set z := sInf S
  have hmz : -2 - |lam| ≤ z := le_csInf hSne fun s hs => hs.1
  have hzE : z ≤ E := csInf_le hSbdd hES
  have hIoc : Set.Ioc z E ⊆ R := by
    rintro t ⟨ht1, ht2⟩
    obtain ⟨s, hs, hst⟩ := exists_lt_of_csInf_lt hSne ht1
    exact hs.2.2 ⟨hst.le, ht2⟩
  have hIcc : Set.Icc z E ⊆ R := by
    rcases eq_or_lt_of_le hzE with h | h
    · rw [h, Set.Icc_self]
      exact Set.singleton_subset_iff.mpr hE
    · rw [← closure_Ioc h.ne]
      exact hRclosed.closure_subset_iff.mpr hIoc
  refine ⟨z, ?_, hzE, hIcc⟩
  rcases eq_or_lt_of_le hmz with hm | hm
  · rw [← hm]
    exact mem_union_left _ (mem_union_left _ (mem_insert_self _ _))
  -- `m < z`: arbitrarily close to the left of `z` there are non-resonant energies of `I_λ`
  have key : ∀ δ > 0, ¬ ∀ t, z - δ < t → t < z → -2 - |lam| ≤ t → t ∈ R := by
    intro δ hδ hall
    have hs : max (-2 - |lam|) (z - δ / 2) ∈ S := by
      refine ⟨le_max_left _ _, (max_lt hm (by linarith)).le.trans hzE, ?_⟩
      rintro t ⟨ht1, ht2⟩
      by_cases htz : t < z
      · exact hall t (by linarith [le_max_right (-2 - |lam|) (z - δ / 2)]) htz
          ((le_max_left _ _).trans ht1)
      · exact hIcc ⟨not_lt.mp htz, ht2⟩
    have h1 : z ≤ max (-2 - |lam|) (z - δ / 2) := csInf_le hSbdd hs
    have h2 : max (-2 - |lam|) (z - δ / 2) < z := max_lt hm (by linarith)
    linarith
  have hzR : z ∈ R := hIcc ⟨le_rfl, hzE⟩
  have hz2 : z ≤ 2 + |lam| := hzE.trans hE.1.2
  rcases hzR.2 with ⟨hR, hfR⟩ | ⟨hL, hfL⟩
  · -- the right end is selected and resonant at `z`: equality holds at `z`
    have heq : |dirL V z a (i + 1)| = θ * |dirL V z a (L + 1)| := by
      by_contra hne
      have hlt : θ * |dirL V z a (L + 1)| < |dirL V z a (i + 1)| :=
        lt_of_le_of_ne hfR (Ne.symm hne)
      have hopen : IsOpen {t : ℝ | θ * |dirL V t a (L + 1)| < |dirL V t a (i + 1)|} :=
        isOpen_lt (continuous_const.mul cP.abs) cpb.abs
      obtain ⟨δ, hδ, hball⟩ := Metric.isOpen_iff.mp hopen z hlt
      refine key δ hδ fun t ht1 ht2 ht3 => ?_
      have htb : t ∈ Metric.ball z δ := by
        rw [Metric.mem_ball, Real.dist_eq, abs_lt]
        constructor <;> linarith
      have ht := hball htb
      simp only [Set.mem_ofPred_eq] at ht
      exact ⟨⟨ht3, by linarith⟩, Or.inl ⟨hR, ht.le⟩⟩
    have hne0 := ResonanceAux.sq_sub_C_mul_sq_ne_zero (dirLpoly V a (i + 1))
      (dirLpoly V a (L + 1)) (θ ^ 2) L (dirLpoly_monic V a L) (natDegree_dirLpoly V a L)
      (by rw [natDegree_dirLpoly]; exact hi) (pow_ne_zero 2 hθ.ne')
    have hroot :
        z ∈ ((dirLpoly V a (i + 1)) ^ 2 - C (θ ^ 2) * (dirLpoly V a (L + 1)) ^ 2).roots := by
      rw [mem_roots hne0, IsRoot.def]
      simp only [eval_sub, eval_pow, eval_mul, eval_C, eval_dirLpoly]
      have h2 := congrArg (fun x : ℝ => x ^ 2) heq
      simp only [mul_pow, sq_abs] at h2
      linarith
    exact mem_union_left _ (mem_union_right _ (Multiset.mem_toFinset.mpr hroot))
  · -- the left end is selected and resonant at `z`: equality holds at `z`
    have heq : |dirR V z (a + L - 1) (L - i)| = θ * |dirL V z a (L + 1)| := by
      by_contra hne
      have hlt : θ * |dirL V z a (L + 1)| < |dirR V z (a + L - 1) (L - i)| :=
        lt_of_le_of_ne hfL (Ne.symm hne)
      have hopen :
          IsOpen {t : ℝ | θ * |dirL V t a (L + 1)| < |dirR V t (a + L - 1) (L - i)|} :=
        isOpen_lt (continuous_const.mul cP.abs) cpa.abs
      obtain ⟨δ, hδ, hball⟩ := Metric.isOpen_iff.mp hopen z hlt
      refine key δ hδ fun t ht1 ht2 ht3 => ?_
      have htb : t ∈ Metric.ball z δ := by
        rw [Metric.mem_ball, Real.dist_eq, abs_lt]
        constructor <;> linarith
      have ht := hball htb
      simp only [Set.mem_ofPred_eq] at ht
      exact ⟨⟨ht3, by linarith⟩, Or.inr ⟨hL, ht.le⟩⟩
    have hne0 := ResonanceAux.sq_sub_C_mul_sq_ne_zero (dirRpoly V (a + L - 1) (L - i))
      (dirLpoly V a (L + 1)) (θ ^ 2) L (dirLpoly_monic V a L) (natDegree_dirLpoly V a L)
      (by rw [ResonanceAux.natDegree_dirRpoly_sub V _ L i hi]; omega) (pow_ne_zero 2 hθ.ne')
    have hroot : z ∈ ((dirRpoly V (a + L - 1) (L - i)) ^ 2 -
        C (θ ^ 2) * (dirLpoly V a (L + 1)) ^ 2).roots := by
      rw [mem_roots hne0, IsRoot.def]
      simp only [eval_sub, eval_pow, eval_mul, eval_C, eval_dirLpoly, eval_dirRpoly]
      have h2 := congrArg (fun x : ℝ => x ^ 2) heq
      simp only [mul_pow, sq_abs] at h2
      linarith
    exact mem_union_right _ (Multiset.mem_toFinset.mpr hroot)

end DFL
