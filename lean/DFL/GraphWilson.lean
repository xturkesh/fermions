import DFL.GraphFlux

/-!
# Two dimensions: Wilson loops and the counting of flux sectors
(paper: Theorem `thm:2dmap`(a), Lemma `lem:2dgauge`(i)–(iii))

* `wilsonLoop z = ∏_{e ∈ z} σᶻ_e` for a cycle `z ∈ ker ∂` (in the `σˣ` basis: it flips the labels of the
  edges of `z`). Wilson loops are involutions, commute with each other, with every Gauss operator and
  with the Hamiltonian; `∏_i G_i = (-1)^{N_f}` (paper: Theorem `thm:2dmap`(a), Eq. `eq:2dparity`).
* On a connected graph, a gauge transformation fixes a sign configuration iff it is constant
  (uniqueness up to a global sign, Lemma `lem:2dgauge`(i)); every gauge orbit has `2^{|V|-1}` elements
  (Lemma `lem:2dgauge`(iii)), and there are `2^β = 2^{|E| - |V| + 1}` orbits, i.e. flux sectors
  (Lemma `lem:2dgauge`(ii) with `flux_eq_iff_gauge`).
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

namespace FinGraph

variable {N M : ℕ} (G : FinGraph N M)

/-- The Wilson loop of `z ∈ Cfg M`: flips the labels of the edges in `z`. -/
def wilsonLoop (_G : FinGraph N M) (z : Cfg M) : Matrix (Cfg N × Cfg M) (Cfg N × Cfg M) ℝ :=
  (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ (fun x' x => if x' = x + z then 1 else 0)

/-! ### Auxiliary lemmas: `ZMod 2` and signs -/

/-- Every element of `ZMod 2` is `0` or `1`. -/
private lemma z2_cases (a : ZMod 2) : a = 0 ∨ a = 1 := by
  revert a; decide

/-- `a + a = 0` in `ZMod 2`. -/
private lemma z2_add_self (a : ZMod 2) : a + a = 0 := by
  revert a; decide

/-- `a + b = 0` forces `a = b` in `ZMod 2`. -/
private lemma z2_eq_of_add_eq_zero {a b : ZMod 2} (h : a + b = 0) : a = b := by
  revert a b; decide

/-- `a ≠ a + 1` in `ZMod 2`. -/
private lemma z2_ne_add_one (a : ZMod 2) : a ≠ a + 1 := by
  revert a; decide

/-- `z + z = 0` for configurations. -/
private lemma cfg_add_self {K : ℕ} (z : Cfg K) : z + z = 0 := by
  funext e
  exact z2_add_self (z e)

/-- `sgn` turns addition in `ZMod 2` into multiplication. -/
private lemma sgn_add (a b : ZMod 2) : sgn (a + b) = sgn a * sgn b := by
  rcases z2_cases a with rfl | rfl <;> rcases z2_cases b with rfl | rfl <;>
    simp [sgn, show (1 : ZMod 2) + 1 = 0 from by decide]

/-- A product of signs is the sign of the sum. -/
private lemma prod_sgn {ι : Type*} (s : Finset ι) (a : ι → ZMod 2) :
    ∏ i ∈ s, sgn (a i) = sgn (∑ i ∈ s, a i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [sgn]
  | insert j s hj ih => rw [Finset.prod_insert hj, Finset.sum_insert hj, ih, sgn_add]

/-! ### Auxiliary lemmas: the boundary map -/

/-- `∂` is additive. -/
private lemma bdry_add (x y : Cfg M) : G.bdry (x + y) = G.bdry x + G.bdry y := by
  funext i
  simp only [bdry, Pi.add_apply, Finset.sum_add_distrib]

/-- A sum over the two endpoints of an edge. -/
private lemma sum_ite_endpoint {R : Type*} [AddCommMonoid R] (e : Fin M) (f : Fin N → R) :
    ∑ i, (if G.src e = i ∨ G.tgt e = i then f i else 0) = f (G.src e) + f (G.tgt e) := by
  have h : ∀ i, (if G.src e = i ∨ G.tgt e = i then f i else 0) =
      (if G.src e = i then f i else 0) + (if G.tgt e = i then f i else 0) := by
    intro i
    by_cases h1 : G.src e = i
    · have h2 : ¬ G.tgt e = i := fun h2 => G.src_ne_tgt e (h1.trans h2.symm)
      rw [ite_eq_left (Or.inl h1), ite_eq_left h1, ite_eq_right h2, add_zero]
    · by_cases h2 : G.tgt e = i
      · rw [ite_eq_left (Or.inr h2), ite_eq_right h1, ite_eq_left h2, zero_add]
      · rw [ite_eq_right (not_or.2 ⟨h1, h2⟩), ite_eq_right h1, ite_eq_right h2, add_zero]
  rw [Finset.sum_congr rfl (fun i _ => h i), Finset.sum_add_distrib, Finset.sum_ite_eq,
    Finset.sum_ite_eq, ite_eq_left (Finset.mem_univ _), ite_eq_left (Finset.mem_univ _)]

/-- `∑_i (∂x)_i = 0`: every edge has exactly two endpoints. -/
private lemma sum_bdry (x : Cfg M) : ∑ i, G.bdry x i = 0 := by
  simp only [bdry, Finset.sum_filter]
  rw [Finset.sum_comm]
  apply Finset.sum_eq_zero
  intro e _
  rw [G.sum_ite_endpoint e (fun _ => x e), z2_add_self]

/-! ### Auxiliary lemmas: translations of edge labels -/

/-- The translation `x ↦ x + z` of edge labels, as a matrix. -/
private def transl (z : Cfg M) : Matrix (Cfg M) (Cfg M) ℝ :=
  fun x' x => if x' = x + z then 1 else 0

/-- Translations compose additively. -/
private lemma transl_mul (z z' : Cfg M) : transl z * transl z' = transl (z + z') := by
  ext x' x
  simp only [transl, Matrix.mul_apply, mul_ite, mul_one, mul_zero]
  rw [Finset.sum_eq_single (x + z')]
  · simp only [ite_true, add_assoc, add_comm z' z]
  · intro y _ hy
    simp [hy]
  · simp

/-- The translation by `0` is the identity. -/
private lemma transl_zero : transl (0 : Cfg M) = 1 := by
  ext x' x
  simp [transl, Matrix.one_apply]

/-- The Wilson loop is `1 ⊗ (translation by z)`. -/
private lemma wilsonLoop_eq (z : Cfg M) :
    G.wilsonLoop z = (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ transl z := rfl

/-- Flipping the label of `e` is adding `δ_e`. -/
private lemma update_add_one (x : Cfg M) (e : Fin M) :
    Function.update x e (x e + 1) = x + Pi.single e 1 := by
  funext k
  by_cases hk : k = e
  · subst hk; simp
  · simp [hk]

/-- `σᶻ_e` is the translation by `δ_e`. -/
private lemma szE_eq_transl (e : Fin M) : G.szE e = transl (Pi.single e 1) := by
  ext x' x
  simp only [szE, transl, update_add_one]

/-- The star terms commute with translations by cycles. -/
private lemma star_mul_transl (i : Fin N) (z : Cfg M) (hz : G.bdry z = 0) :
    G.star i * transl z = transl z * G.star i := by
  ext x' x
  rw [FinGraph.star, diagonal_mul, mul_diagonal]
  simp only [transl]
  by_cases hx : x' = x + z
  · subst hx
    rw [G.bdry_add, hz, add_zero]
    ring
  · simp [hx]

/-! ### Wilson loops -/

/-- Wilson loops are involutions. -/
theorem wilsonLoop_mul_self (z : Cfg M) : G.wilsonLoop z * G.wilsonLoop z = 1 := by
  rw [wilsonLoop_eq, ← mul_kronecker_mul, Matrix.one_mul, transl_mul, cfg_add_self, transl_zero,
    one_kronecker_one]

/-- Wilson loops commute with each other. -/
theorem wilsonLoop_comm (z z' : Cfg M) :
    G.wilsonLoop z * G.wilsonLoop z' = G.wilsonLoop z' * G.wilsonLoop z := by
  rw [wilsonLoop_eq, wilsonLoop_eq, ← mul_kronecker_mul, ← mul_kronecker_mul, transl_mul,
    transl_mul, add_comm]

/-- Wilson loops of cycles commute with the Gauss operators. -/
theorem wilsonLoop_comm_gaussOp (z : Cfg M) (hz : G.bdry z = 0) (i : Fin N) :
    G.gaussOp i * G.wilsonLoop z = G.wilsonLoop z * G.gaussOp i := by
  ext ⟨n', x'⟩ ⟨n, x⟩
  rw [gaussOp, diagonal_mul, mul_diagonal, wilsonLoop_eq, kronecker_apply, Matrix.one_apply]
  simp only [transl]
  by_cases hn : n' = n
  · by_cases hx : x' = x + z
    · subst hn hx
      simp only [G.bdry_add, hz, add_zero]
      ring
    · simp [hx]
  · simp [hn]

/-- Wilson loops of cycles commute with the Hamiltonian (paper: Theorem `thm:2dmap`(a)). -/
theorem wilsonLoop_comm_ham (J h : ℝ) (z : Cfg M) (hz : G.bdry z = 0) :
    G.ham J h * G.wilsonLoop z = G.wilsonLoop z * G.ham J h := by
  have h1 : ∀ e, (G.hopE e ⊗ₖ G.szE e) * G.wilsonLoop z =
      G.wilsonLoop z * (G.hopE e ⊗ₖ G.szE e) := by
    intro e
    rw [wilsonLoop_eq, ← mul_kronecker_mul, ← mul_kronecker_mul, Matrix.mul_one, Matrix.one_mul,
      szE_eq_transl, transl_mul, transl_mul, add_comm]
  have h2 : ∀ i, ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ G.star i) * G.wilsonLoop z =
      G.wilsonLoop z * ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ G.star i) := by
    intro i
    rw [wilsonLoop_eq, ← mul_kronecker_mul, ← mul_kronecker_mul, G.star_mul_transl i z hz]
  simp only [ham, Matrix.sub_mul, Matrix.mul_sub, Matrix.smul_mul, Matrix.mul_smul,
    Matrix.sum_mul, Matrix.mul_sum, h1, h2]

/-- `∏_i G_i = (-1)^{N_f}` (paper: Eq. `eq:2dparity`); the Gauss operators are diagonal and this is the
product of their entries. -/
theorem gaussOp_prod (p : Cfg N × Cfg M) : ∏ i, G.gaussOp i p p = sgn (∑ i, p.1 i) := by
  simp only [gaussOp, diagonal_apply_eq]
  rw [Finset.prod_mul_distrib, prod_sgn, prod_sgn, G.sum_bdry, mul_right_eq_self₀]
  left
  simp [sgn]

/-- Adjacency of vertices. -/
def Adj (i j : Fin N) : Prop := ∃ e, (G.src e = i ∧ G.tgt e = j) ∨ (G.tgt e = i ∧ G.src e = j)

/-- The graph is connected. -/
def Connected : Prop := ∀ i j : Fin N, Relation.ReflTransGen G.Adj i j

/-! ### Gauge orbits -/

/-- Composition of gauge transformations. -/
private lemma gaugeAct_gaugeAct (θ θ' : Cfg N) (s : Cfg M) :
    G.gaugeAct θ (G.gaugeAct θ' s) = G.gaugeAct (θ + θ') s := by
  funext e
  simp only [gaugeAct, Pi.add_apply]
  ring

/-- The trivial gauge transformation. -/
private lemma gaugeAct_zero (s : Cfg M) : G.gaugeAct 0 s = s := by
  funext e
  simp [gaugeAct]

/-- **Uniqueness of the gauge transformation up to a global sign** (paper: Lemma `lem:2dgauge`(i)): on
a connected graph, `θ` fixes `s` iff `θ` is constant. -/
theorem gaugeAct_eq_self_iff (hG : G.Connected) (θ : Cfg N) (s : Cfg M) :
    G.gaugeAct θ s = s ↔ ∀ i j, θ i = θ j := by
  constructor
  · intro h
    have hedge : ∀ e, θ (G.src e) = θ (G.tgt e) := by
      intro e
      have h1 := congrFun h e
      simp only [gaugeAct] at h1
      apply z2_eq_of_add_eq_zero
      have h2 : s e + (θ (G.src e) + θ (G.tgt e)) = s e + 0 := by
        rw [← add_assoc, h1, add_zero]
      exact add_left_cancel h2
    have hadj : ∀ i j, G.Adj i j → θ i = θ j := by
      rintro i j ⟨e, ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩⟩
      · exact hedge e
      · exact (hedge e).symm
    intro i j
    induction hG i j with
    | refl => rfl
    | tail _ hbc ih => exact ih.trans (hadj _ _ hbc)
  · intro h
    funext e
    simp only [gaugeAct]
    rw [h (G.src e) (G.tgt e), add_assoc, z2_add_self, add_zero]

/-- Two gauge transformations act in the same way iff they differ by a constant. -/
private lemma gaugeAct_eq_gaugeAct_iff (hG : G.Connected) (θ θ' : Cfg N) (s : Cfg M) :
    G.gaugeAct θ s = G.gaugeAct θ' s ↔ ∀ i j, (θ + θ') i = (θ + θ') j := by
  rw [← G.gaugeAct_eq_self_iff hG (θ + θ') (G.gaugeAct θ' s), gaugeAct_gaugeAct, add_assoc,
    cfg_add_self, add_zero]

/-- **Gauge orbits** (paper: Lemma `lem:2dgauge`(iii)): on a connected graph with `N ≥ 1` vertices,
every gauge orbit has `2^{N-1}` elements. -/
theorem card_gaugeOrbit (hG : G.Connected) (hN : 0 < N) (s : Cfg M) :
    (univ.image (fun θ : Cfg N => G.gaugeAct θ s)).card = 2 ^ (N - 1) := by
  set v0 : Fin N := ⟨0, hN⟩
  have hfib : ∀ y ∈ univ.image (fun θ : Cfg N => G.gaugeAct θ s),
      (univ.filter (fun θ : Cfg N => G.gaugeAct θ s = y)).card = 2 := by
    intro y hy
    obtain ⟨θ0, -, rfl⟩ := Finset.mem_image.1 hy
    have hset : univ.filter (fun θ : Cfg N => G.gaugeAct θ s = G.gaugeAct θ0 s) =
        {θ0, θ0 + 1} := by
      ext θ
      simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_insert,
        Finset.mem_singleton]
      rw [G.gaugeAct_eq_gaugeAct_iff hG]
      constructor
      · intro h
        rcases z2_cases ((θ + θ0) v0) with h0 | h0
        · left
          funext i
          apply z2_eq_of_add_eq_zero
          have := h i v0
          rw [h0] at this
          exact this
        · right
          funext i
          have := h i v0
          rw [h0] at this
          simp only [Pi.add_apply, Pi.one_apply] at this ⊢
          rw [← this, add_comm (θ i), ← add_assoc, z2_add_self, zero_add]
      · have h1 : ∀ a : ZMod 2, a + 1 + a = 1 := by decide
        rintro (rfl | rfl) i j
        · simp [z2_add_self]
        · simp only [Pi.add_apply, Pi.one_apply]
          rw [h1, h1]
    rw [hset, Finset.card_pair]
    intro h
    exact z2_ne_add_one (θ0 v0) (congrFun h v0)
  have h := Finset.card_eq_sum_card_image (fun θ : Cfg N => G.gaugeAct θ s) univ
  rw [Finset.sum_congr rfl hfib, Finset.sum_const, smul_eq_mul, Finset.card_univ,
    Fintype.card_fun, ZMod.card, Fintype.card_fin] at h
  have h2 : 2 ^ N = 2 ^ (N - 1) * 2 := by
    rw [← pow_succ]
    congr 1
    omega
  rw [h2] at h
  exact (Nat.eq_of_mul_eq_mul_right (by norm_num) h).symm

/-- **Number of flux sectors** (paper: Lemma `lem:2dgauge`(ii)): on a connected graph with `N ≥ 1`
vertices there are `2^β`, `β = M - N + 1`, gauge orbits of sign configurations. -/
theorem card_fluxClasses (hG : G.Connected) (hN : 0 < N) :
    (univ.image (fun s : Cfg M => univ.image (fun θ : Cfg N => G.gaugeAct θ s))).card =
      2 ^ (M + 1 - N) := by
  -- an element of an orbit has the same orbit
  have hmem : ∀ s s0 : Cfg M, s ∈ univ.image (fun θ : Cfg N => G.gaugeAct θ s0) ↔
      univ.image (fun θ : Cfg N => G.gaugeAct θ s) =
        univ.image (fun θ : Cfg N => G.gaugeAct θ s0) := by
    intro s s0
    constructor
    · intro hs
      obtain ⟨θ, -, rfl⟩ := Finset.mem_image.1 hs
      ext t
      simp only [Finset.mem_image, Finset.mem_univ, true_and]
      constructor
      · rintro ⟨θ', rfl⟩
        exact ⟨θ' + θ, (G.gaugeAct_gaugeAct θ' θ s0).symm⟩
      · rintro ⟨θ', rfl⟩
        refine ⟨θ' + θ, ?_⟩
        rw [gaugeAct_gaugeAct, add_assoc, cfg_add_self, add_zero]
    · intro h
      rw [← h]
      exact Finset.mem_image.2 ⟨0, Finset.mem_univ _, G.gaugeAct_zero s⟩
  have hfib : ∀ o ∈ univ.image (fun s : Cfg M => univ.image (fun θ : Cfg N => G.gaugeAct θ s)),
      (univ.filter (fun s : Cfg M =>
        univ.image (fun θ : Cfg N => G.gaugeAct θ s) = o)).card = 2 ^ (N - 1) := by
    intro o ho
    obtain ⟨s0, -, rfl⟩ := Finset.mem_image.1 ho
    have hset : univ.filter (fun s : Cfg M => univ.image (fun θ : Cfg N => G.gaugeAct θ s) =
        univ.image (fun θ : Cfg N => G.gaugeAct θ s0)) =
        univ.image (fun θ : Cfg N => G.gaugeAct θ s0) := by
      ext s
      simp only [Finset.mem_filter, Finset.mem_univ, true_and]
      exact (hmem s s0).symm
    rw [hset]
    exact G.card_gaugeOrbit hG hN s0
  have h := Finset.card_eq_sum_card_image
    (fun s : Cfg M => univ.image (fun θ : Cfg N => G.gaugeAct θ s)) univ
  rw [Finset.sum_congr rfl hfib, Finset.sum_const, smul_eq_mul, Finset.card_univ,
    Fintype.card_fun, ZMod.card, Fintype.card_fin] at h
  set k := (univ.image (fun s : Cfg M => univ.image (fun θ : Cfg N => G.gaugeAct θ s))).card
  have hk : 0 < k := Finset.card_pos.2 (Finset.image_nonempty.2 Finset.univ_nonempty)
  have hle : N - 1 ≤ M := by
    have h3 : 2 ^ (N - 1) ≤ 2 ^ M := by
      rw [h]
      exact Nat.le_mul_of_pos_left _ hk
    exact (Nat.pow_le_pow_iff_right (by norm_num)).1 h3
  have h2 : 2 ^ M = 2 ^ (M + 1 - N) * 2 ^ (N - 1) := by
    rw [← pow_add]
    congr 1
    omega
  rw [h2] at h
  exact (Nat.eq_of_mul_eq_mul_right (by positivity) h).symm

end FinGraph

end DFL
