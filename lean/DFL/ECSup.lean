import DFL.ECProof

/-!
# Uniformity of (EC) over the volume (paper: Theorem `thm:EC`(ii))

For a window `Λ = [a, a + L)` of the chain `[0, N)` let `Q^Λ(x, y)` be the eigenfunction correlator of the
Bernoulli chain on `Λ` (charges `q|_Λ`), `0` unless `x, y ∈ Λ` (`corrQwin`). Under the large-deviation input:
`E sup_{Λ ∋ x, y} Q^Λ(x, y) ≤ C e^{-η|x-y|}`, with the supremum over all windows of `[0, N)`, uniformly in `N`
(`ec_sup_of_ldt`). Since the charges are i.i.d., this is the statement of the paper for all finite intervals
`Λ ⊂ ℤ` containing `x` and `y`; with windows `[0, L)` it gives `E sup_L Q^{[0,L)}(x, y) ≤ C e^{-η|x-y|}` for one
sequence of charges.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The charges of `q` on the window `[a, a + L)`, relabelled as `Fin L`. -/
def winCfg (q : Cfg N) (a L : ℕ) (h : a + L ≤ N) : Cfg L := fun i => q ⟨a + i.val, by omega⟩

/-- `Q^Λ(x, y)` for the window `Λ = [a, a + L)` (`0` unless `x, y ∈ Λ`). -/
def corrQwin (J h : ℝ) (q : Cfg N) (a L : ℕ) (hL : a + L ≤ N) (x y : Fin N) : ℝ :=
  if hxy : a ≤ x.val ∧ x.val < a + L ∧ a ≤ y.val ∧ y.val < a + L then
    corrQ (hq J h (winCfg q a L hL)) (hq_isHermitian J h _) ⟨x.val - a, by omega⟩ ⟨y.val - a, by omega⟩
  else 0

/-- `sup_{Λ ∋ x, y} Q^Λ(x, y)`, the supremum over all windows `Λ = [a, a + L)` of `[0, N)`. -/
def corrQsup (J h : ℝ) (q : Cfg N) (x y : Fin N) : ℝ :=
  (univ : Finset (Fin (N + 1) × Fin (N + 1))).sup' univ_nonempty
    (fun p => if hp : p.1.val + p.2.val ≤ N then corrQwin J h q p.1.val p.2.val hp x y else 0)

namespace ECSupAux

/-! ### Translation of boxes -/

/-- `dirL` of a translated potential. -/
lemma dirL_shift (V : ℤ → ℝ) (E : ℝ) (a s : ℤ) (k : ℕ) :
    dirL (fun z => V (z + s)) E a k = dirL V E (a + s) k := by
  induction k using Nat.twoStepInduction with
  | zero => rfl
  | one => rfl
  | more n ih0 ih1 =>
    simp only [dirL]
    rw [ih0, ih1, show a + (n : ℤ) + s = a + s + n by ring]

/-- `dirR` of a translated potential. -/
lemma dirR_shift (V : ℤ → ℝ) (E : ℝ) (b s : ℤ) (k : ℕ) :
    dirR (fun z => V (z + s)) E b k = dirR V E (b + s) k := by
  induction k using Nat.twoStepInduction with
  | zero => rfl
  | one => rfl
  | more n ih0 ih1 =>
    simp only [dirR]
    rw [ih0, ih1, show b - (n : ℤ) + s = b + s - n by ring]

/-- The resonant set of a translated potential is the translated resonant set. -/
lemma resSet_shift (V : ℤ → ℝ) (lam : ℝ) (a s : ℤ) (L i : ℕ) (endL endR : Bool) (θ : ℝ) :
    resSet (fun z => V (z + s)) lam a L i endL endR θ = resSet V lam (a + s) L i endL endR θ := by
  ext E
  simp only [resSet, Set.mem_ofPred_eq, dirL_shift, dirR_shift]
  rw [show a + (L : ℤ) - 1 + s = a + s + L - 1 by ring]

/-! ### Windows -/

/-- The potential of a window is the translated potential of the chain, inside the window. -/
lemma potZ_winCfg (lam : ℝ) (q : Cfg N) (a L : ℕ) (hL : a + L ≤ N) (z : ℤ) (hz0 : 0 ≤ z)
    (hzL : z < L) : potZ lam (winCfg q a L hL) z = potZ lam q (z + a) := by
  unfold potZ extZ
  rw [dite_eq_left ⟨hz0, hzL⟩, dite_eq_left ⟨by omega, by omega⟩]
  refine congrArg (fun i => lam * sgn (q i)) (Fin.ext ?_)
  show a + z.toNat = (z + a).toNat
  omega

/-- The resonant set of a box inside a window, in the coordinates of the chain. -/
lemma resSet_winCfg (lam : ℝ) (q : Cfg N) (a L : ℕ) (hL : a + L ≤ N) (b : ℤ) (L1 i : ℕ)
    (endL endR : Bool) (θ : ℝ) (hb : 0 ≤ b) (hbL : b + L1 ≤ L) (hi : i < L1) :
    resSet (potZ lam (winCfg q a L hL)) lam b L1 i endL endR θ =
      resSet (potZ lam q) lam (b + a) L1 i endL endR θ := by
  rw [resSet_congr (potZ lam (winCfg q a L hL)) (fun z => potZ lam q (z + a)) lam b L1 i endL endR
    θ hi (fun j hj => potZ_winCfg lam q a L hL _ (by omega) (by omega))]
  exact resSet_shift _ _ _ _ _ _ _ _ _

/-- `Q^Λ(x, y) ≤ 1`. -/
lemma corrQwin_le_one (J h : ℝ) (q : Cfg N) (a L : ℕ) (hL : a + L ≤ N) (x y : Fin N) :
    corrQwin J h q a L hL x y ≤ 1 := by
  unfold corrQwin
  split_ifs
  · exact ECProofAux.corrQ_le_one _ _ _ _
  · exact zero_le_one

/-- `Q^Λ` is symmetric. -/
lemma corrQwin_symm (J h : ℝ) (q : Cfg N) (a L : ℕ) (hL : a + L ≤ N) (x y : Fin N) :
    corrQwin J h q a L hL x y = corrQwin J h q a L hL y x := by
  unfold corrQwin
  by_cases hxy : a ≤ x.val ∧ x.val < a + L ∧ a ≤ y.val ∧ y.val < a + L
  · have hyx : a ≤ y.val ∧ y.val < a + L ∧ a ≤ x.val ∧ x.val < a + L :=
      ⟨hxy.2.2.1, hxy.2.2.2, hxy.1, hxy.2.1⟩
    rw [dite_eq_left hxy, dite_eq_left hyx]
    exact ECProofAux.corrQ_symm _ _ _ _
  · have hyx : ¬ (a ≤ y.val ∧ y.val < a + L ∧ a ≤ x.val ∧ x.val < a + L) :=
      fun h' => hxy ⟨h'.2.2.1, h'.2.2.2, h'.1, h'.2.1⟩
    rw [dite_eq_right hxy, dite_eq_right hyx]

/-- A bound on all windows bounds the supremum. -/
lemma corrQsup_le (J h : ℝ) (q : Cfg N) (x y : Fin N) (B : ℝ) (hB : 0 ≤ B)
    (hw : ∀ (a L : ℕ) (hL : a + L ≤ N), corrQwin J h q a L hL x y ≤ B) :
    corrQsup J h q x y ≤ B := by
  unfold corrQsup
  refine Finset.sup'_le _ _ fun p _ => ?_
  split_ifs with hp
  · exact hw _ _ hp
  · exact hB

/-- `sup_Λ Q^Λ(x, y) ≤ 1`. -/
lemma corrQsup_le_one (J h : ℝ) (q : Cfg N) (x y : Fin N) : corrQsup J h q x y ≤ 1 :=
  corrQsup_le J h q x y 1 zero_le_one fun a L hL => corrQwin_le_one J h q a L hL x y

/-- `sup_Λ Q^Λ` is symmetric. -/
lemma corrQsup_symm (J h : ℝ) (q : Cfg N) (x y : Fin N) :
    corrQsup J h q x y = corrQsup J h q y x := by
  unfold corrQsup
  refine Finset.sup'_congr _ rfl fun p _ => ?_
  split_ifs with hp
  · exact corrQwin_symm J h q _ _ hp x y
  · rfl

/-! ### The bad events -/

/-- The bad event of the boxes around `x` and `y` determined by the left end `ac` and the right end
`bc` (exclusive) of a window, clamped to `[x-n-1, x]` and `[y+1, y+n+2]`: the two boxes share a resonant
energy. -/
def badE (lam θ : ℝ) (x y n ac bc : ℕ) (q : Cfg N) : Prop :=
  (resSet (potZ lam q) lam ((max (x - n) ac : ℕ) : ℤ) (x + n + 1 - max (x - n) ac)
      (x - max (x - n) ac) (decide (ac + n + 1 ≤ x)) true θ ∩
    resSet (potZ lam q) lam ((y - n : ℕ) : ℤ) (n + 1 + min n (bc - 1 - y)) n true
      (decide (y + n + 2 ≤ bc)) θ).Nonempty

open Classical in
/-- The union of the bad events over the at most `(n + 2)²` clamped ends of the windows. -/
def badSet (lam θ : ℝ) (x y n : ℕ) : Finset (Cfg N) :=
  ((Finset.Icc (x - n - 1) x) ×ˢ ((Finset.Icc (y + 1) (y + n + 2)).filter (· ≤ N))).biUnion
    (fun p => univ.filter (fun q => badE lam θ x y n p.1 p.2 q))

/-- **Deterministic step for a window** (paper: Eq. `eq:Qgood` for `Λ`): if the window contains `x`
and `y`, and the boxes determined by its clamped ends share no resonant energy, then
`Q^Λ(x, y) ≤ 4θ`. -/
lemma corrQwin_le_of_good (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) (x y : Fin N) (n : ℕ) (θ : ℝ)
    (hθ : 0 ≤ θ) (h2n : x.val + 2 * n + 1 ≤ y.val) (a L : ℕ) (hL : a + L ≤ N)
    (hin : a ≤ x.val ∧ x.val < a + L ∧ a ≤ y.val ∧ y.val < a + L)
    (hgood : ¬ badE (-(2 * h / J)) θ x.val y.val n (max a (x.val - n - 1))
      (min (a + L) (y.val + n + 2)) q) :
    corrQwin J h q a L hL x y ≤ 4 * θ := by
  unfold corrQwin
  rw [dite_eq_left hin]
  obtain ⟨hax, hxL, hay, hyL⟩ := hin
  set ac := max a (x.val - n - 1) with hac
  set bc := min (a + L) (y.val + n + 2) with hbc
  set s1 := max (x.val - n) ac with hs1
  have hs1a : a ≤ s1 := by omega
  have hs1x : s1 ≤ x.val := by omega
  have hs1n : x.val - n ≤ s1 := by omega
  have hbcy : y.val + 1 ≤ bc := by omega
  have hbcL : bc ≤ a + L := by omega
  have hw1 := resSet_winCfg (-(2 * h / J)) q a L hL ((s1 : ℤ) - a) (x.val + n + 1 - s1)
    (x.val - s1) (decide (ac + n + 1 ≤ x.val)) true θ (by omega) (by omega) (by omega)
  have hw2 := resSet_winCfg (-(2 * h / J)) q a L hL (((y.val - n : ℕ) : ℤ) - a)
    (n + 1 + min n (bc - 1 - y.val)) n true (decide (y.val + n + 2 ≤ bc)) θ (by omega) (by omega)
    (by omega)
  rw [sub_add_cancel] at hw1 hw2
  refine ECProofAux.corrQ_le_of_good J h hJ (winCfg q a L hL) ⟨x.val - a, by omega⟩
    ⟨y.val - a, by omega⟩ θ hθ ((s1 : ℤ) - a) (x.val + n + 1 - s1) (x.val - s1)
    (decide (ac + n + 1 ≤ x.val)) true (((y.val - n : ℕ) : ℤ) - a) (n + 1 + min n (bc - 1 - y.val)) n
    true (decide (y.val + n + 2 ≤ bc)) (by omega) (by omega) (by omega) (by simp only; omega)
    (Or.inr rfl) (fun hd => by simp only [decide_eq_false_iff_not] at hd; omega)
    (fun hd => by simp at hd) (by omega) (by omega) (by omega) (by simp only; omega)
    (Or.inl rfl) (fun hd => by simp at hd)
    (fun hd => by simp only [decide_eq_false_iff_not] at hd; omega) ?_
  rw [hw1, hw2]
  exact hgood

open Classical in
/-- **Probabilistic step, uniformly over windows** (paper: proof of Theorem `thm:EC`(ii)): the union
of the bad events has probability at most `(n+2)² · 4(8n+6) C₁ e^{-c₁ n}`. -/
lemma card_badSet_le (lam κ C1 c1 : ℝ) (m0 : ℕ) (hEB : EntryBound lam κ C1 c1 m0) (hC1 : 0 ≤ C1)
    (x y n : ℕ) (h2n : x + 2 * n + 1 ≤ y) (hn : m0 ≤ n) :
    ((badSet (N := N) lam (Real.exp (-(κ * n))) x y n).card : ℝ) ≤
      ((n : ℝ) + 2) ^ 2 * (4 * (8 * n + 6) * C1 * Real.exp (-c1 * n)) * 2 ^ N := by
  unfold badSet
  set I := (Finset.Icc (x - n - 1) x) ×ˢ ((Finset.Icc (y + 1) (y + n + 2)).filter (· ≤ N)) with hI_def
  have hI : (I.card : ℝ) ≤ ((n : ℝ) + 2) ^ 2 := by
    have h1 : (Finset.Icc (x - n - 1) x).card ≤ n + 2 := by rw [Nat.card_Icc]; omega
    have h2 : ((Finset.Icc (y + 1) (y + n + 2)).filter (· ≤ N)).card ≤ n + 2 :=
      (Finset.card_filter_le _ _).trans (by rw [Nat.card_Icc]; omega)
    have h3 : I.card ≤ (n + 2) ^ 2 := by
      rw [hI_def, Finset.card_product, sq]
      exact Nat.mul_le_mul h1 h2
    exact_mod_cast h3
  have hterm : ∀ p ∈ I, ((univ.filter (fun q : Cfg N =>
      badE lam (Real.exp (-(κ * n))) x y n p.1 p.2 q)).card : ℝ) ≤
      4 * (8 * n + 6) * C1 * Real.exp (-c1 * n) * 2 ^ N := by
    intro p hp
    simp only [hI_def, Finset.mem_product, Finset.mem_Icc, Finset.mem_filter] at hp
    obtain ⟨⟨hp1, hp2⟩, ⟨hp3, hp4⟩, hp5⟩ := hp
    exact two_box lam κ C1 c1 m0 hEB hC1 ((max (x - n) p.1 : ℕ) : ℤ) (x + n + 1 - max (x - n) p.1)
      (x - max (x - n) p.1) (decide (p.1 + n + 1 ≤ x)) true ((y - n : ℕ) : ℤ)
      (n + 1 + min n (p.2 - 1 - y)) n true (decide (y + n + 2 ≤ p.2)) n (by positivity)
      (by omega) (by omega) (by omega) (by omega) (by omega) (by omega) hn (fun _ => by omega)
      (fun hd => by simp only [decide_eq_true_eq] at hd; omega)
      (fun hd => by simp only [decide_eq_true_eq] at hd; omega) (fun _ => rfl)
  have hB0 : 0 ≤ 4 * (8 * (n : ℝ) + 6) * C1 * Real.exp (-c1 * n) * 2 ^ N := by positivity
  calc ((I.biUnion (fun p => univ.filter (fun q : Cfg N =>
        badE lam (Real.exp (-(κ * n))) x y n p.1 p.2 q))).card : ℝ)
      ≤ ∑ p ∈ I, ((univ.filter (fun q : Cfg N =>
        badE lam (Real.exp (-(κ * n))) x y n p.1 p.2 q)).card : ℝ) := by
        exact_mod_cast Finset.card_biUnion_le
    _ ≤ ∑ _p ∈ I, 4 * (8 * (n : ℝ) + 6) * C1 * Real.exp (-c1 * n) * 2 ^ N :=
        Finset.sum_le_sum hterm
    _ = I.card * (4 * (8 * (n : ℝ) + 6) * C1 * Real.exp (-c1 * n) * 2 ^ N) := by
        rw [Finset.sum_const, nsmul_eq_mul]
    _ ≤ ((n : ℝ) + 2) ^ 2 * (4 * (8 * (n : ℝ) + 6) * C1 * Real.exp (-c1 * n) * 2 ^ N) :=
        mul_le_mul_of_nonneg_right hI hB0
    _ = _ := by ring

/-- **(EC) at distance `n`, uniformly over the windows** (paper: Eq. `eq:ECsupn`): for `x < y` and
`n = ⌊(d-1)/2⌋ ≥ m₀`, `E sup_Λ Q^Λ(x, y) ≤ 4 e^{-κn} + (n+2)² 4(8n+6) C₁ e^{-c₁n}`. -/
lemma avgQ_sup_le_of_entry (J h : ℝ) (hJ : J ≠ 0) (κ C1 c1 : ℝ) (m0 : ℕ)
    (hEB : EntryBound (-(2 * h / J)) κ C1 c1 m0) (hC1 : 0 ≤ C1) (x y : Fin N)
    (hxy : x.val < y.val) (n : ℕ) (hn_def : n = (y.val - x.val - 1) / 2) (hn : m0 ≤ n) :
    avgQ (fun q => corrQsup J h q x y) ≤
      4 * Real.exp (-(κ * n)) + ((n : ℝ) + 2) ^ 2 * (4 * (8 * n + 6) * C1 * Real.exp (-c1 * n)) := by
  have h2n : x.val + 2 * n + 1 ≤ y.val := by omega
  refine ECProofAux.avgQ_le_add _ (badSet (-(2 * h / J)) (Real.exp (-(κ * n))) x.val y.val n) _ _
    (fun q => corrQsup_le_one J h q x y) (fun q hq => ?_) (by positivity)
    (card_badSet_le _ κ C1 c1 m0 hEB hC1 x.val y.val n h2n hn)
  refine corrQsup_le J h q x y _ (by positivity) fun a L hL => ?_
  by_cases hin : a ≤ x.val ∧ x.val < a + L ∧ a ≤ y.val ∧ y.val < a + L
  · refine corrQwin_le_of_good J h hJ q x y n _ (Real.exp_pos _).le h2n a L hL hin fun hbad => hq ?_
    unfold badSet
    rw [Finset.mem_biUnion]
    refine ⟨(max a (x.val - n - 1), min (a + L) (y.val + n + 2)), ?_, ?_⟩
    · simp only [Finset.mem_product, Finset.mem_Icc, Finset.mem_filter]
      omega
    · simp only [Finset.mem_filter, Finset.mem_univ, true_and]
      exact hbad
  · unfold corrQwin
    rw [dite_eq_right hin]
    positivity

/-- `(n+2)² (8n+6) e^{-c₁ n} ≤ 384 e^{c₁} c₁⁻³ e^{-cn}` for `0 < c ≤ c₁/2`. -/
lemma poly3_exp_le (c1 c : ℝ) (hc1 : 0 < c1) (hc : c ≤ c1 / 2) (n : ℕ) :
    ((n : ℝ) + 2) ^ 2 * (8 * n + 6) * Real.exp (-c1 * n) ≤
      384 * Real.exp c1 / c1 ^ 3 * Real.exp (-(c * n)) := by
  have hn0 : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have hx : 0 ≤ c1 * ((n : ℝ) + 2) / 2 := by positivity
  have hf := Real.pow_div_factorial_le_exp _ hx 3
  have h6 : ((Nat.factorial 3 : ℕ) : ℝ) = 6 := by norm_num [Nat.factorial]
  rw [h6] at hf
  have hexp : Real.exp (c1 * ((n : ℝ) + 2) / 2) = Real.exp c1 * Real.exp (c1 * n / 2) := by
    rw [← Real.exp_add]
    congr 1
    ring
  rw [hexp] at hf
  have hc13 : 0 < c1 ^ 3 := by positivity
  have hcube : ((n : ℝ) + 2) ^ 3 ≤ 48 / c1 ^ 3 * (Real.exp c1 * Real.exp (c1 * n / 2)) := by
    rw [div_mul_eq_mul_div, le_div_iff₀ hc13]
    have e : (c1 * ((n : ℝ) + 2) / 2) ^ 3 / 6 = c1 ^ 3 * ((n : ℝ) + 2) ^ 3 / 48 := by ring
    rw [e] at hf
    linarith
  have hpoly : ((n : ℝ) + 2) ^ 2 * (8 * n + 6) ≤ 8 * ((n : ℝ) + 2) ^ 3 := by
    have : 0 ≤ ((n : ℝ) + 2) ^ 2 := by positivity
    nlinarith
  have hmono : Real.exp (-(c1 / 2 * n)) ≤ Real.exp (-(c * n)) := by
    apply Real.exp_le_exp.mpr
    nlinarith
  have hK0 : 0 ≤ 384 * Real.exp c1 / c1 ^ 3 := by positivity
  calc ((n : ℝ) + 2) ^ 2 * (8 * n + 6) * Real.exp (-c1 * n)
      ≤ 8 * (48 / c1 ^ 3 * (Real.exp c1 * Real.exp (c1 * n / 2))) * Real.exp (-c1 * n) :=
        mul_le_mul_of_nonneg_right (by linarith) (Real.exp_pos _).le
    _ = 384 * Real.exp c1 / c1 ^ 3 * (Real.exp (c1 * n / 2) * Real.exp (-c1 * n)) := by ring
    _ = 384 * Real.exp c1 / c1 ^ 3 * Real.exp (-(c1 / 2 * n)) := by
        rw [← Real.exp_add]
        congr 2
        ring
    _ ≤ 384 * Real.exp c1 / c1 ^ 3 * Real.exp (-(c * n)) := mul_le_mul_of_nonneg_left hmono hK0

end ECSupAux

/-- **(EC) uniformly over the volume** (paper: Theorem `thm:EC`(ii), eq. `eq:ECsup`):
`E sup_{Λ ∋ x, y} Q^Λ(x, y) ≤ C e^{-η|x-y|}`, uniformly in `N`. -/
theorem ec_sup_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ ∀ (N : ℕ) (x y : Fin N),
      avgQ (fun q : Cfg N => corrQsup J h q x y) ≤ C * Real.exp (-η * Nat.dist x.val y.val) := by
  obtain ⟨κ, C1, c1, m0, hκ, hc1, hC1, hEB⟩ := entry_ldt _ hL
  set c := min κ (c1 / 2) with hc_def
  have hc : 0 < c := lt_min hκ (by linarith)
  have hcκ : c ≤ κ := min_le_left _ _
  have hcc1 : c ≤ c1 / 2 := min_le_right _ _
  set K := 384 * Real.exp c1 / c1 ^ 3 with hK_def
  have hK : 0 ≤ K := by positivity
  set C' := 4 + 4 * C1 * K with hC'_def
  have hC' : 0 ≤ C' := by positivity
  set η := c / 2 with hη_def
  have hη : 0 < η := by positivity
  refine ⟨max (C' * Real.exp c) (Real.exp (η * (2 * m0 + 1))), η, hη, ?_⟩
  -- the case `x < y` with `n ≥ m₀`
  have key : ∀ (N : ℕ) (x y : Fin N), x.val < y.val → 2 * m0 < y.val - x.val →
      avgQ (fun q : Cfg N => corrQsup J h q x y) ≤
        C' * Real.exp c * Real.exp (-η * ((y.val - x.val : ℕ) : ℝ)) := by
    intro N x y hxy hd
    obtain ⟨n, hn_def⟩ : ∃ n : ℕ, n = (y.val - x.val - 1) / 2 := ⟨_, rfl⟩
    have hn : m0 ≤ n := by omega
    refine (ECSupAux.avgQ_sup_le_of_entry J h hJ κ C1 c1 m0 hEB hC1 x y hxy n hn_def hn).trans ?_
    have hn0 : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    have hnR : ((y.val - x.val : ℕ) : ℝ) ≤ 2 * (n : ℝ) + 2 := by
      have : y.val - x.val ≤ 2 * n + 2 := by omega
      exact_mod_cast this
    have e1 : Real.exp (-(κ * n)) ≤ Real.exp (-(c * n)) := by
      apply Real.exp_le_exp.mpr
      nlinarith
    have e2 := ECSupAux.poly3_exp_le c1 c hc1 hcc1 n
    have e3 : Real.exp (-(c * n)) ≤ Real.exp c * Real.exp (-η * ((y.val - x.val : ℕ) : ℝ)) := by
      rw [← Real.exp_add]
      apply Real.exp_le_exp.mpr
      rw [hη_def]
      nlinarith
    calc 4 * Real.exp (-(κ * n)) + ((n : ℝ) + 2) ^ 2 * (4 * (8 * n + 6) * C1 * Real.exp (-c1 * n))
        = 4 * Real.exp (-(κ * n)) +
            4 * C1 * (((n : ℝ) + 2) ^ 2 * (8 * n + 6) * Real.exp (-c1 * n)) := by ring
      _ ≤ 4 * Real.exp (-(c * n)) + 4 * C1 * (K * Real.exp (-(c * n))) := by
          gcongr
      _ = C' * Real.exp (-(c * n)) := by rw [hC'_def]; ring
      _ ≤ C' * (Real.exp c * Real.exp (-η * ((y.val - x.val : ℕ) : ℝ))) := by gcongr
      _ = C' * Real.exp c * Real.exp (-η * ((y.val - x.val : ℕ) : ℝ)) := by ring
  -- short distances
  have small : ∀ (N : ℕ) (x y : Fin N), Nat.dist x.val y.val ≤ 2 * m0 →
      avgQ (fun q : Cfg N => corrQsup J h q x y) ≤
        Real.exp (η * (2 * m0 + 1)) * Real.exp (-η * Nat.dist x.val y.val) := by
    intro N x y hd
    refine (ECProofAux.avgQ_le_of_le _ 1 fun q => ECSupAux.corrQsup_le_one J h q x y).trans ?_
    rw [← Real.exp_add]
    apply Real.one_le_exp
    have : (Nat.dist x.val y.val : ℝ) ≤ 2 * m0 := by exact_mod_cast hd
    nlinarith
  intro N x y
  rcases le_or_gt (Nat.dist x.val y.val) (2 * m0) with hd | hd
  · exact (small N x y hd).trans
      (mul_le_mul_of_nonneg_right (le_max_right _ _) (Real.exp_pos _).le)
  · rcases lt_trichotomy x.val y.val with hxy | hxy | hxy
    · rw [Nat.dist_eq_sub_of_le hxy.le] at hd ⊢
      exact (key N x y hxy hd).trans
        (mul_le_mul_of_nonneg_right (le_max_left _ _) (Real.exp_pos _).le)
    · rw [hxy, Nat.dist_self] at hd
      omega
    · rw [Nat.dist_eq_sub_of_le_right hxy.le] at hd ⊢
      have hsym : avgQ (fun q : Cfg N => corrQsup J h q x y) =
          avgQ (fun q : Cfg N => corrQsup J h q y x) := by
        unfold avgQ
        congr 1
        exact Finset.sum_congr rfl fun q _ => ECSupAux.corrQsup_symm J h q x y
      rw [hsym]
      exact (key N y x hxy hd).trans
        (mul_le_mul_of_nonneg_right (le_max_left _ _) (Real.exp_pos _).le)

end DFL
