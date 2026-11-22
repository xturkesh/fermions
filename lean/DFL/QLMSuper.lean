import DFL.QLM

/-!
# Quantum links: superpositions of background-charge sectors and arbitrary dynamics
(paper: Proposition `prop:gausscap`(c),(d))

The summed Gauss law `N_A = E_b - E_{a-1} + o_A - Q_A` holds configuration by configuration, with
`Q_A = ∑_{j ∈ A} G_j` the total background charge of `A = [a, b]`. Hence, for every normalized state
(not necessarily in one sector), `S_N ≤ ln(4S+1) + H(Q_A)` in the bulk and `S_N ≤ ln(2S+1) + H(Q_A)` at
the fixed-flux edge. Every unitary commuting with the Gauss operators preserves each sector and the
law of `Q_A`, so these bounds, and the bound of `qlm_numberEntropy_le` for states in one sector, hold
along any gauge-invariant dynamics.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The Gauss operator `G_j = E_j - E_{j-1} - ρ_j` evaluated on a configuration. -/
def gaussVal (S : ℕ) (E0 : ℤ) (c : Cfg N × (Fin N → Fin (2 * S + 1))) (j : Fin N) : ℤ :=
  efield S (c.2 j) - prevField S E0 c.2 j - stagCharge c.1 j

/-- `G_j` as a diagonal matrix. -/
def gaussQLM (S : ℕ) (E0 : ℤ) (j : Fin N) :
    Matrix (Cfg N × (Fin N → Fin (2 * S + 1))) (Cfg N × (Fin N → Fin (2 * S + 1))) ℂ :=
  diagonal (fun c => ((gaussVal S E0 c j : ℤ) : ℂ))

/-- The total background charge `Q_A = ∑_{j ∈ [a, b]} G_j` of a configuration. -/
def bgCharge (S : ℕ) (E0 : ℤ) (a b : ℕ) (c : Cfg N × (Fin N → Fin (2 * S + 1))) : ℤ :=
  ∑ j ∈ univ.filter (fun j : Fin N => a ≤ j.val ∧ j.val ≤ b), gaussVal S E0 c j

/-- The Shannon entropy `H(Q_A)` of the law of `Q_A` in the state `v`. -/
def bgEntropy (S : ℕ) (E0 : ℤ) (a b : ℕ) (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ) : ℝ :=
  ∑ y ∈ univ.image (bgCharge (N := N) S E0 a b),
    Real.negMulLog (∑ c ∈ univ.filter (fun c => bgCharge S E0 a b c = y), Complex.normSq (v c))

/-! ### Auxiliary lemmas -/

namespace QLMSuperAux

/-- Subadditivity of `η = negMulLog` on `[0, ∞)`. -/
lemma negMulLog_add_le {x y : ℝ} (hx : 0 ≤ x) (hy : 0 ≤ y) :
    Real.negMulLog (x + y) ≤ Real.negMulLog x + Real.negMulLog y := by
  have h1 : x * Real.log x ≤ x * Real.log (x + y) := by
    rcases hx.eq_or_lt with rfl | hx'
    · simp
    · exact mul_le_mul_of_nonneg_left (Real.log_le_log hx' (by linarith)) hx
  have h2 : y * Real.log y ≤ y * Real.log (x + y) := by
    rcases hy.eq_or_lt with rfl | hy'
    · simp
    · exact mul_le_mul_of_nonneg_left (Real.log_le_log hy' (by linarith)) hy
  simp only [Real.negMulLog]
  have h3 : -(x + y) * Real.log (x + y) = -(x * Real.log (x + y)) - y * Real.log (x + y) := by
    ring
  rw [h3]
  linarith

/-- Subadditivity of `η` over a finite family of nonnegative numbers (grouping). -/
lemma negMulLog_sum_le {ι : Type*} (s : Finset ι) (p : ι → ℝ) (hp : ∀ i ∈ s, 0 ≤ p i) :
    Real.negMulLog (∑ i ∈ s, p i) ≤ ∑ i ∈ s, Real.negMulLog (p i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | insert j s hj ih =>
    rw [Finset.sum_insert hj, Finset.sum_insert hj]
    have hpj : 0 ≤ p j := hp j (Finset.mem_insert_self j s)
    have hps : ∀ i ∈ s, 0 ≤ p i := fun i hi => hp i (Finset.mem_insert_of_mem hi)
    calc Real.negMulLog (p j + ∑ i ∈ s, p i)
        ≤ Real.negMulLog (p j) + Real.negMulLog (∑ i ∈ s, p i) :=
          negMulLog_add_le hpj (Finset.sum_nonneg hps)
      _ ≤ Real.negMulLog (p j) + ∑ i ∈ s, Real.negMulLog (p i) := by linarith [ih hps]

/-- Termwise Gibbs inequality: `-x log x ≤ w - x - x log w` for `x ≥ 0`, `w > 0`. -/
lemma negMulLog_le_gibbs {x w : ℝ} (hx : 0 ≤ x) (hw : 0 < w) :
    Real.negMulLog x ≤ w - x - x * Real.log w := by
  rcases hx.eq_or_lt with rfl | hx
  · simp only [Real.negMulLog_zero, sub_zero, zero_mul]
    exact hw.le
  · have h := Real.log_le_sub_one_of_pos (div_pos hw hx)
    rw [Real.log_div hw.ne' hx.ne'] at h
    have h2 : x * (Real.log w - Real.log x) ≤ x * (w / x - 1) :=
      mul_le_mul_of_nonneg_left h hx.le
    have h3 : x * (w / x - 1) = w - x := by field_simp
    rw [Real.negMulLog]
    nlinarith

/-- Conditional entropy bound: `∑_{i ∈ s} η(p_i) ≤ η(∑ p_i) + (∑ p_i) log L` for nonnegative `p`
and `|s| ≤ L`. -/
lemma sum_negMulLog_le {ι : Type*} (s : Finset ι) (p : ι → ℝ) (hp : ∀ i ∈ s, 0 ≤ p i)
    {L : ℝ} (hL : (s.card : ℝ) ≤ L) :
    ∑ i ∈ s, Real.negMulLog (p i) ≤
      Real.negMulLog (∑ i ∈ s, p i) + (∑ i ∈ s, p i) * Real.log L := by
  have hq0 : 0 ≤ ∑ i ∈ s, p i := Finset.sum_nonneg hp
  rcases hq0.eq_or_lt with hq0 | hqpos
  · -- all terms vanish
    have hz : ∀ i ∈ s, p i = 0 := (Finset.sum_eq_zero_iff_of_nonneg hp).1 hq0.symm
    rw [← hq0, Real.negMulLog_zero, zero_mul, add_zero]
    exact le_of_eq (Finset.sum_eq_zero (fun i hi => by rw [hz i hi, Real.negMulLog_zero]))
  · have hs : s.Nonempty := by
      by_contra hne
      rw [Finset.not_nonempty_iff_eq_empty] at hne
      rw [hne, Finset.sum_empty] at hqpos
      exact lt_irrefl 0 hqpos
    have hn : (0 : ℝ) < s.card := by exact_mod_cast hs.card_pos
    have hwpos : 0 < (∑ i ∈ s, p i) / s.card := div_pos hqpos hn
    have hsum : ∑ i ∈ s, Real.negMulLog (p i) ≤
        ∑ i ∈ s, ((∑ i ∈ s, p i) / s.card - p i - p i * Real.log ((∑ i ∈ s, p i) / s.card)) :=
      Finset.sum_le_sum (fun i hi => negMulLog_le_gibbs (hp i hi) hwpos)
    rw [Finset.sum_sub_distrib, Finset.sum_sub_distrib, Finset.sum_const, nsmul_eq_mul,
      ← Finset.sum_mul] at hsum
    have hcw : (s.card : ℝ) * ((∑ i ∈ s, p i) / s.card) = ∑ i ∈ s, p i := by
      field_simp
    have hlogw : Real.log ((∑ i ∈ s, p i) / s.card) = Real.log (∑ i ∈ s, p i) - Real.log s.card :=
      Real.log_div hqpos.ne' hn.ne'
    have hlogL : Real.log s.card ≤ Real.log L := Real.log_le_log hn hL
    rw [hcw, hlogw] at hsum
    have hneg : Real.negMulLog (∑ i ∈ s, p i) = -(∑ i ∈ s, p i) * Real.log (∑ i ∈ s, p i) := rfl
    rw [hneg]
    nlinarith [mul_le_mul_of_nonneg_left hlogL hq0]

/-- Variant of `sum_negMulLog_le` for `p` vanishing on `s` outside a set `t` with `|t| ≤ L`. -/
lemma sum_negMulLog_le_of_supp {ι : Type*} [DecidableEq ι] (s t : Finset ι) (p : ι → ℝ)
    (hp : ∀ i ∈ s, 0 ≤ p i) (hsupp : ∀ i ∈ s, i ∉ t → p i = 0) {L : ℝ}
    (hL : (t.card : ℝ) ≤ L) :
    ∑ i ∈ s, Real.negMulLog (p i) ≤
      Real.negMulLog (∑ i ∈ s, p i) + (∑ i ∈ s, p i) * Real.log L := by
  have h1 : ∑ i ∈ s, Real.negMulLog (p i) = ∑ i ∈ s ∩ t, Real.negMulLog (p i) := by
    symm
    apply Finset.sum_subset Finset.inter_subset_left
    intro i hi hni
    have hit : i ∉ t := fun ht => hni (Finset.mem_inter.2 ⟨hi, ht⟩)
    rw [hsupp i hi hit, Real.negMulLog_zero]
  have h2 : ∑ i ∈ s, p i = ∑ i ∈ s ∩ t, p i := by
    symm
    apply Finset.sum_subset Finset.inter_subset_left
    intro i hi hni
    exact hsupp i hi (fun ht => hni (Finset.mem_inter.2 ⟨hi, ht⟩))
  rw [h1, h2]
  apply sum_negMulLog_le
  · intro i hi
    exact hp i (Finset.mem_inter.1 hi).1
  · calc ((s ∩ t).card : ℝ) ≤ t.card := by
          exact_mod_cast Finset.card_le_card Finset.inter_subset_right
      _ ≤ L := hL

/-- The number `o_A` of odd sites of `[a, b]`. -/
def oddCount (N : ℕ) (a b : ℕ) : ℤ :=
  ∑ j ∈ univ.filter (fun j : Fin N => a ≤ j.val ∧ j.val ≤ b), (if j.val % 2 = 1 then (1 : ℤ) else 0)

/-- **Summed Gauss law, configuration by configuration**: `N_A = E_b - E_{a-1} + o_A - Q_A`. -/
lemma countInt_eq (S : ℕ) (E0 : ℤ) (c : Cfg N × (Fin N → Fin (2 * S + 1))) (a b : Fin N)
    (hab : a ≤ b) :
    countInt a.val b.val c.1 = efield S (c.2 b) - prevField S E0 c.2 a + oddCount N a.val b.val -
      bgCharge S E0 a.val b.val c := by
  have hG : GaussLaw S E0 (fun j => gaussVal S E0 c j) c.1 c.2 := by
    intro j
    unfold gaussVal
    ring
  rw [countInt_gauss S E0 _ c.1 c.2 hG a b hab, Finset.sum_sub_distrib]
  unfold oddCount bgCharge
  ring

/-- If, given `Q_A = y`, the particle number of `[a, b]` lies in a window `[L y, L y + W]`, then
`S_N ≤ ln(W + 1) + H(Q_A)`. -/
lemma shannon_le_log_add_bgEntropy (S : ℕ) (E0 : ℤ) (a b : ℕ)
    (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ) (hnorm : ∑ c, Complex.normSq (v c) = 1)
    (L : ℤ → ℤ) (W : ℕ)
    (hwin : ∀ c : Cfg N × (Fin N → Fin (2 * S + 1)), L (bgCharge S E0 a b c) ≤ countInt a b c.1 ∧
      countInt a b c.1 ≤ L (bgCharge S E0 a b c) + W) :
    shannon N (qlmNumDist S a b v) ≤ Real.log ((W : ℝ) + 1) + bgEntropy S E0 a b v := by
  classical
  -- joint masses of `(N_A, Q_A)` and the law of `Q_A`
  let m : ℕ → ℤ → ℝ := fun k y =>
    ∑ c ∈ univ.filter (fun c : Cfg N × (Fin N → Fin (2 * S + 1)) =>
      countInt a b c.1 = k ∧ bgCharge S E0 a b c = y), Complex.normSq (v c)
  let q : ℤ → ℝ := fun y =>
    ∑ c ∈ univ.filter (fun c : Cfg N × (Fin N → Fin (2 * S + 1)) =>
      bgCharge S E0 a b c = y), Complex.normSq (v c)
  let Y : Finset ℤ := univ.image (bgCharge (N := N) S E0 a b)
  have hm0 : ∀ k y, 0 ≤ m k y := fun k y =>
    Finset.sum_nonneg (fun c _ => Complex.normSq_nonneg _)
  -- the marginal of `N_A`
  have hmarg1 : ∀ k, qlmNumDist S a b v k = ∑ y ∈ Y, m k y := by
    intro k
    unfold qlmNumDist
    rw [← Finset.sum_fiberwise_of_maps_to (g := bgCharge S E0 a b) (t := Y)
      (fun c _ => Finset.mem_image_of_mem _ (Finset.mem_univ c))]
    refine Finset.sum_congr rfl (fun y _ => ?_)
    rw [Finset.filter_filter]
  -- the marginal of `Q_A`
  have hmarg2 : ∀ y, ∑ k ∈ range (N + 1), m k y = q y := by
    intro y
    have h := Finset.sum_fiberwise_of_maps_to
      (s := univ.filter (fun c : Cfg N × (Fin N → Fin (2 * S + 1)) => bgCharge S E0 a b c = y))
      (g := fun c : Cfg N × (Fin N → Fin (2 * S + 1)) => (countInt a b c.1).toNat)
      (t := range (N + 1))
      (fun c _ => by
        have := countInt_le a b c.1
        rw [Finset.mem_range]
        omega) (fun c => Complex.normSq (v c))
    refine Eq.trans ?_ h
    refine Finset.sum_congr rfl (fun k _ => ?_)
    rw [Finset.filter_filter]
    apply Finset.sum_congr _ (fun _ _ => rfl)
    ext c
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    have := countInt_nonneg a b c.1
    constructor
    · rintro ⟨h1, h2⟩
      exact ⟨h2, by omega⟩
    · rintro ⟨h1, h2⟩
      exact ⟨by omega, h1⟩
  -- given `Q_A = y`, `N_A` lies in a window of `W + 1` values
  have hwin' : ∀ y, ∀ k ∈ range (N + 1),
      k ∉ (Finset.Icc (L y) (L y + W)).image Int.toNat → m k y = 0 := by
    intro y k _ hk
    apply Finset.sum_eq_zero
    intro c hc
    exfalso
    rw [Finset.mem_filter] at hc
    obtain ⟨-, hck, hcy⟩ := hc
    obtain ⟨h1, h2⟩ := hwin c
    rw [hcy, hck] at h1 h2
    exact hk (Finset.mem_image.2 ⟨(k : ℤ), Finset.mem_Icc.2 ⟨h1, h2⟩, by simp⟩)
  have hcard : ∀ y, (((Finset.Icc (L y) (L y + W)).image Int.toNat).card : ℝ) ≤ (W : ℝ) + 1 := by
    intro y
    have h : ((Finset.Icc (L y) (L y + W)).image Int.toNat).card ≤ W + 1 := by
      calc ((Finset.Icc (L y) (L y + W)).image Int.toNat).card
          ≤ (Finset.Icc (L y) (L y + W)).card := Finset.card_image_le
        _ = W + 1 := by rw [Int.card_Icc]; omega
    exact_mod_cast h
  have hq1 : ∑ y ∈ Y, q y = 1 := by
    rw [← hnorm]
    exact Finset.sum_fiberwise_of_maps_to
      (fun c _ => Finset.mem_image_of_mem _ (Finset.mem_univ c)) _
  calc shannon N (qlmNumDist S a b v)
      = ∑ k ∈ range (N + 1), Real.negMulLog (∑ y ∈ Y, m k y) := by
        unfold shannon
        exact Finset.sum_congr rfl (fun k _ => by rw [hmarg1])
    _ ≤ ∑ k ∈ range (N + 1), ∑ y ∈ Y, Real.negMulLog (m k y) :=
        Finset.sum_le_sum (fun k _ => negMulLog_sum_le _ _ (fun y _ => hm0 k y))
    _ = ∑ y ∈ Y, ∑ k ∈ range (N + 1), Real.negMulLog (m k y) := Finset.sum_comm
    _ ≤ ∑ y ∈ Y, (Real.negMulLog (q y) + q y * Real.log ((W : ℝ) + 1)) := by
        apply Finset.sum_le_sum
        intro y _
        have h := sum_negMulLog_le_of_supp (range (N + 1))
          ((Finset.Icc (L y) (L y + W)).image Int.toNat) (fun k => m k y)
          (fun k _ => hm0 k y) (hwin' y) (hcard y)
        rw [hmarg2] at h
        exact h
    _ = bgEntropy S E0 a b v + (∑ y ∈ Y, q y) * Real.log ((W : ℝ) + 1) := by
        rw [Finset.sum_add_distrib, ← Finset.sum_mul]
        rfl
    _ = Real.log ((W : ℝ) + 1) + bgEntropy S E0 a b v := by
        rw [hq1, one_mul, add_comm]

/-- A matrix commuting with `G_j` has vanishing entries between configurations with different
values of `G_j`. -/
lemma gaussVal_eq_of_comm (S : ℕ) (E0 : ℤ) (j : Fin N)
    (U : Matrix (Cfg N × (Fin N → Fin (2 * S + 1))) (Cfg N × (Fin N → Fin (2 * S + 1))) ℂ)
    (hU : U * gaussQLM S E0 j = gaussQLM S E0 j * U) (c' c : Cfg N × (Fin N → Fin (2 * S + 1)))
    (h : U c' c ≠ 0) : gaussVal S E0 c' j = gaussVal S E0 c j := by
  have h1 := congrFun (congrFun hU c') c
  rw [gaussQLM, mul_diagonal, diagonal_mul] at h1
  have h2 : U c' c * (((gaussVal S E0 c j : ℤ) : ℂ) - ((gaussVal S E0 c' j : ℤ) : ℂ)) = 0 := by
    rw [mul_sub, h1]
    ring
  rcases mul_eq_zero.1 h2 with h3 | h3
  · exact absurd h3 h
  · exact_mod_cast (sub_eq_zero.1 h3).symm

/-- A matrix commuting with all `G_j` has vanishing entries between configurations with different
background charges `Q_A`. -/
lemma bgCharge_eq_of_comm (S : ℕ) (E0 : ℤ) (a b : ℕ)
    (U : Matrix (Cfg N × (Fin N → Fin (2 * S + 1))) (Cfg N × (Fin N → Fin (2 * S + 1))) ℂ)
    (hU : ∀ j, U * gaussQLM S E0 j = gaussQLM S E0 j * U)
    (c' c : Cfg N × (Fin N → Fin (2 * S + 1))) (h : U c' c ≠ 0) :
    bgCharge S E0 a b c' = bgCharge S E0 a b c := by
  unfold bgCharge
  exact Finset.sum_congr rfl (fun j _ => gaussVal_eq_of_comm S E0 j U (hU j) c' c h)

/-- An isometry preserves `∑ |w|²`. -/
lemma sum_normSq_mulVec {ι : Type*} [Fintype ι] [DecidableEq ι] (U : Matrix ι ι ℂ)
    (hUu : Uᴴ * U = 1) (w : ι → ℂ) :
    ∑ i, Complex.normSq ((U *ᵥ w) i) = ∑ i, Complex.normSq (w i) := by
  have h : star (U *ᵥ w) ⬝ᵥ (U *ᵥ w) = star w ⬝ᵥ w := by
    rw [Matrix.star_mulVec, ← Matrix.dotProduct_mulVec, Matrix.mulVec_mulVec, hUu,
      Matrix.one_mulVec]
  have h2 : ∀ u : ι → ℂ, star u ⬝ᵥ u = ((∑ i, Complex.normSq (u i) : ℝ) : ℂ) := by
    intro u
    simp only [dotProduct, Pi.star_apply, Complex.ofReal_sum]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [Complex.normSq_eq_conj_mul_self]
    rfl
  rw [h2, h2] at h
  exact_mod_cast h

end QLMSuperAux

/-- `GaussLaw` is the statement `G_j = Q_j` for all `j`. -/
theorem gaussLaw_iff (S : ℕ) (E0 : ℤ) (Q : Fin N → ℤ) (c : Cfg N × (Fin N → Fin (2 * S + 1))) :
    GaussLaw S E0 Q c.1 c.2 ↔ ∀ j, gaussVal S E0 c j = Q j := by
  unfold GaussLaw gaussVal
  constructor
  · intro h j
    have := h j
    linarith
  · intro h j
    have := h j
    linarith

/-- **Superpositions of sectors, bulk** (paper: Proposition `prop:gausscap`(c)):
`S_N ≤ ln(4S+1) + H(Q_A)` for every normalized state and every interval `[a, b]` with `a ≥ 1`. -/
theorem qlm_numberEntropy_le_super (S : ℕ) (E0 : ℤ) (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ)
    (hnorm : ∑ c, Complex.normSq (v c) = 1) (a b : Fin N) (hab : a ≤ b) (ha : a.val ≠ 0) :
    shannon N (qlmNumDist S a.val b.val v) ≤
      Real.log (4 * S + 1) + bgEntropy S E0 a.val b.val v := by
  have h := QLMSuperAux.shannon_le_log_add_bgEntropy S E0 a.val b.val v hnorm
    (fun y => QLMSuperAux.oddCount N a.val b.val - 2 * S - y) (4 * S) (fun c => by
      have hc := QLMSuperAux.countInt_eq S E0 c a b hab
      have hp := prevField_bounds S E0 c.2 a ha
      have he := efield_bounds S (c.2 b)
      rw [hc]
      push_cast
      constructor <;> linarith [hp.1, hp.2, he.1, he.2])
  push_cast at h
  exact h

/-- **Superpositions of sectors, fixed-flux edge** (paper: Proposition `prop:gausscap`(c)):
`S_N ≤ ln(2S+1) + H(Q_A)` for every normalized state and every interval `[0, b]`. -/
theorem qlm_numberEntropy_le_super_edge (S : ℕ) (E0 : ℤ)
    (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ) (hnorm : ∑ c, Complex.normSq (v c) = 1)
    (b : Fin N) :
    shannon N (qlmNumDist S 0 b.val v) ≤ Real.log (2 * S + 1) + bgEntropy S E0 0 b.val v := by
  have hN : 0 < N := Nat.lt_of_le_of_lt (Nat.zero_le _) b.isLt
  let a : Fin N := ⟨0, hN⟩
  have hab : a ≤ b := Nat.zero_le _
  have h := QLMSuperAux.shannon_le_log_add_bgEntropy S E0 a.val b.val v hnorm
    (fun y => QLMSuperAux.oddCount N a.val b.val - E0 - S - y) (2 * S) (fun c => by
      have hc := QLMSuperAux.countInt_eq S E0 c a b hab
      have hp : prevField S E0 c.2 a = E0 := by
        unfold prevField
        rw [dite_eq_left rfl]
      have he := efield_bounds S (c.2 b)
      rw [hc, hp]
      push_cast
      constructor <;> linarith [he.1, he.2])
  push_cast at h
  exact h

/-- **Any gauge-invariant dynamics preserves the sectors** (paper: Proposition `prop:gausscap`(d)). -/
theorem qlm_sector_preserved (S : ℕ) (E0 : ℤ) (Q : Fin N → ℤ)
    (U : Matrix (Cfg N × (Fin N → Fin (2 * S + 1))) (Cfg N × (Fin N → Fin (2 * S + 1))) ℂ)
    (hU : ∀ j, U * gaussQLM S E0 j = gaussQLM S E0 j * U)
    (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ) (hsupp : ∀ c, v c ≠ 0 → GaussLaw S E0 Q c.1 c.2) :
    ∀ c, (U *ᵥ v) c ≠ 0 → GaussLaw S E0 Q c.1 c.2 := by
  intro c hc
  rw [gaussLaw_iff]
  intro j
  obtain ⟨c', -, hc'⟩ := Finset.exists_ne_zero_of_sum_ne_zero hc
  have hU' : U c c' ≠ 0 := left_ne_zero_of_mul hc'
  have hv : v c' ≠ 0 := right_ne_zero_of_mul hc'
  rw [QLMSuperAux.gaussVal_eq_of_comm S E0 j U (hU j) c c' hU']
  exact (gaussLaw_iff S E0 Q c').1 (hsupp c' hv) j

/-- **Any gauge-invariant dynamics preserves the law of `Q_A`** (paper: Proposition
`prop:gausscap`(d)). -/
theorem bgEntropy_invariant (S : ℕ) (E0 : ℤ) (a b : ℕ)
    (U : Matrix (Cfg N × (Fin N → Fin (2 * S + 1))) (Cfg N × (Fin N → Fin (2 * S + 1))) ℂ)
    (hU : ∀ j, U * gaussQLM S E0 j = gaussQLM S E0 j * U) (hUu : Uᴴ * U = 1)
    (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ) :
    bgEntropy S E0 a b (U *ᵥ v) = bgEntropy S E0 a b v := by
  unfold bgEntropy
  refine Finset.sum_congr rfl (fun y _ => ?_)
  congr 1
  -- the projection onto `{Q_A = y}` commutes with `U`
  let P : (Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ) → (Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ) :=
    fun u c => if bgCharge S E0 a b c = y then u c else 0
  have hPU : P (U *ᵥ v) = U *ᵥ P v := by
    funext c'
    simp only [P, mulVec, dotProduct]
    split_ifs with hc'
    · refine Finset.sum_congr rfl (fun c _ => ?_)
      by_cases hc : bgCharge S E0 a b c = y
      · rw [ite_eq_left hc]
      · rw [ite_eq_right hc]
        have hU0 : U c' c = 0 := by
          by_contra hne
          exact hc ((QLMSuperAux.bgCharge_eq_of_comm S E0 a b U hU c' c hne).symm.trans hc')
        rw [hU0, zero_mul, zero_mul]
    · symm
      apply Finset.sum_eq_zero
      intro c _
      by_cases hc : bgCharge S E0 a b c = y
      · rw [ite_eq_left hc]
        have hU0 : U c' c = 0 := by
          by_contra hne
          exact hc' ((QLMSuperAux.bgCharge_eq_of_comm S E0 a b U hU c' c hne).trans hc)
        rw [hU0, zero_mul]
      · rw [ite_eq_right hc, mul_zero]
  have hsumP : ∀ u : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ,
      ∑ c ∈ univ.filter (fun c => bgCharge S E0 a b c = y), Complex.normSq (u c) =
        ∑ c, Complex.normSq (P u c) := by
    intro u
    rw [Finset.sum_filter]
    refine Finset.sum_congr rfl (fun c _ => ?_)
    simp only [P]
    split_ifs <;> simp
  rw [hsumP, hsumP, hPU, QLMSuperAux.sum_normSq_mulVec U hUu]

end DFL
