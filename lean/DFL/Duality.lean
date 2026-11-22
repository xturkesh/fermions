import DFL.Basic

/-!
# The exact free-fermion duality (paper: Theorem `thm:duality`, Corollary `cor:sector`)

We prove that the 1D ℤ₂ LGT Hamiltonian `H J h` is unitarily equivalent — by a permutation of
the configuration basis `(n, x) ↔ (n, q)` — to the block-diagonal free-fermion Hamiltonian
`Hdual J h = ⊕_q (dΓ(h_q) - h ∑_j q_j)`, and that the Gauss operators become the diagonal
charge operators `τˣ_j = diag(sgn q_j)`.

Main results:
* `DFL.hop_apply`: on an open chain, nearest-neighbour fermionic hopping has no Jordan–Wigner
  sign;
* `DFL.dGamma_hq`: explicit form of `dΓ(h_q)`;
* `DFL.dualEquiv`: the bijection `(n, q) ↦ (n, linkOf n q)`;
* `DFL.duality`: `(H J h).submatrix dualEquiv dualEquiv = Hdual J h`;
* `DFL.gauss_dual`: `(gauss j).submatrix dualEquiv dualEquiv = diagonal (sgn ∘ q_j)`;
* `DFL.numOp_dual`: the matter number operators are unchanged by the duality.
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-! ### Elementary facts about `ZMod 2`, `sgn` and `occ` (private helpers) -/

/-- Every element of `ZMod 2` is `0` or `1`. -/
private lemma zmod2_eq_zero_or_one (a : ZMod 2) : a = 0 ∨ a = 1 := by
  revert a; decide

/-- `1 + 1 = 0` in `ZMod 2`. -/
private lemma zmod2_one_add_one : (1 : ZMod 2) + 1 = 0 := by decide

/-- `a + a = 0` in `ZMod 2`. -/
private lemma zmod2_add_self (a : ZMod 2) : a + a = 0 := by
  revert a; decide

/-- In `ZMod 2`, `a ≠ b` iff `b = a + 1`. -/
private lemma zmod2_ne_iff (a b : ZMod 2) : a ≠ b ↔ b = a + 1 := by
  revert a b; decide

/-- `sgn` turns addition in `ZMod 2` into multiplication. -/
private lemma sgn_add (a b : ZMod 2) : sgn (a + b) = sgn a * sgn b := by
  rcases zmod2_eq_zero_or_one a with rfl | rfl <;>
  rcases zmod2_eq_zero_or_one b with rfl | rfl <;> simp [sgn, zmod2_one_add_one]

/-- `sgn a ^ 2 = 1`. -/
private lemma sgn_mul_self (a : ZMod 2) : sgn a * sgn a = 1 := by
  rcases zmod2_eq_zero_or_one a with rfl | rfl <;> simp [sgn]

/-- `(-1)^a = 1 - 2a` for `a ∈ {0, 1}`. -/
private lemma sgn_eq_one_sub_two_occ (a : ZMod 2) : sgn a = 1 - 2 * occ a := by
  rcases zmod2_eq_zero_or_one a with rfl | rfl <;> norm_num [sgn, occ]

/-! ### Fermionic matrix elements -/

/-- Matrix elements of `c_a† c_b`. -/
private lemma cre_mul_ann_apply (a b : Fin N) (n' n : Cfg N) :
    (cre a * ann b) n' n =
      if n b = 1 ∧ n' a = 1 ∧ Function.update n b 0 = Function.update n' a 0 then
        sgn (jwPar n' a) * sgn (jwPar n b) else 0 := by
  simp only [cre, Matrix.mul_apply, Matrix.transpose_apply, ann]
  rw [Finset.sum_eq_single (Function.update n b 0)]
  · by_cases h1 : n b = 1 <;> by_cases h2 : n' a = 1 <;>
      by_cases h3 : Function.update n b 0 = Function.update n' a 0 <;> simp [h1, h2, h3]
  · intro m _ hm
    simp [hm]
  · simp

/-- Two configurations with the same value at `x` that agree after emptying `x` are equal. -/
private lemma eq_of_update_eq {x : Fin N} {n n' : Cfg N} (hx : n' x = n x)
    (h : Function.update n x 0 = Function.update n' x 0) : n' = n := by
  funext k
  by_cases hk : k = x
  · subst hk; exact hx
  · have := congrFun h k
    simpa [Function.update_of_ne hk] using this.symm

/-- `c_x† c_x` is the diagonal occupation operator. -/
theorem cre_mul_ann_self (x : Fin N) : cre x * ann x = diagonal (fun n => occ (n x)) := by
  ext n' n
  rw [cre_mul_ann_apply, diagonal_apply]
  by_cases hn : n x = 1
  · by_cases h' : n' x = 1 ∧ Function.update n x 0 = Function.update n' x 0
    · have : n' = n := eq_of_update_eq (h'.1.trans hn.symm) h'.2
      subst this
      simp [hn, sgn_mul_self, occ]
    · have : n' ≠ n := by rintro rfl; exact h' ⟨hn, rfl⟩
      simp only [hn, true_and, h', this, ite_false]
  · have h0 : n x = 0 := (zmod2_eq_zero_or_one _).resolve_right hn
    by_cases h' : n' = n <;> simp [h', h0, occ]

/-- The support condition of `c_a† c_b` for `a ≠ b`: `n_a = 0`, `n_b = 1` and
`n' = n ∘ swap a b`. -/
private lemma cre_ann_cond {a b : Fin N} (hab : a ≠ b) (n' n : Cfg N) :
    (n b = 1 ∧ n' a = 1 ∧ Function.update n b 0 = Function.update n' a 0) ↔
      (n a = 0 ∧ n b = 1 ∧ n' = n ∘ Equiv.swap a b) := by
  constructor
  · rintro ⟨hb, ha, h⟩
    have hna : n a = 0 := by
      have := congrFun h a
      simpa [Function.update_of_ne hab] using this
    refine ⟨hna, hb, ?_⟩
    funext k
    by_cases hka : k = a
    · subst hka; simp [ha, hb]
    · by_cases hkb : k = b
      · subst hkb
        have := congrFun h k
        simp [Function.update_of_ne hka] at this
        simp [← this, hna]
      · have := congrFun h k
        simp [Function.update_of_ne hka, Function.update_of_ne hkb] at this
        simp [Equiv.swap_apply_of_ne_of_ne hka hkb, this]
  · rintro ⟨ha, hb, rfl⟩
    refine ⟨hb, by simp [hb], ?_⟩
    funext k
    by_cases hka : k = a
    · subst hka; simp [Function.update_of_ne hab, ha]
    · by_cases hkb : k = b
      · subst hkb; simp [Function.update_of_ne hka, ha]
      · simp [Function.update_of_ne hka, Function.update_of_ne hkb,
          Equiv.swap_apply_of_ne_of_ne hka hkb]

/-- The Jordan–Wigner parity only depends on the modes to the left of `x`. -/
private lemma jwPar_congr {n m : Cfg N} {x : Fin N} (h : ∀ k, k < x → n k = m k) :
    jwPar n x = jwPar m x :=
  Finset.sum_congr rfl (fun k hk => h k (Finset.mem_filter.1 hk).2)

/-- Recursion for the Jordan–Wigner parity: `jwPar n (j+1) = jwPar n j + n j`. -/
private lemma jwPar_succ (n : Cfg N) {j j' : Fin N} (hj' : j'.val = j.val + 1) :
    jwPar n j' = jwPar n j + n j := by
  unfold jwPar
  have : univ.filter (fun k : Fin N => k < j') = insert j (univ.filter (fun k => k < j)) := by
    ext k; simp only [mem_filter, mem_univ, true_and, mem_insert, Fin.lt_def, Fin.ext_iff]; omega
  rw [this, Finset.sum_insert (by simp), add_comm]

/-- Swapping the modes `j, j+1` does not change the parity to the left of `j`. -/
private lemma jwPar_comp_swap (n : Cfg N) {j j' : Fin N} (hj' : j'.val = j.val + 1) :
    jwPar (n ∘ Equiv.swap j j') j = jwPar n j := by
  apply jwPar_congr
  intro k hk
  have h1 : k ≠ j := ne_of_lt hk
  have h2 : k ≠ j' := by intro h; rw [Fin.lt_def, h] at hk; omega
  simp [Equiv.swap_apply_of_ne_of_ne h1 h2]

/-- Matrix elements of the hopping term: on an open chain nearest-neighbour hopping is sign-free,
`⟨n'| c_j† c_{j+1} + c_{j+1}† c_j |n⟩ = 1` iff `n_j ≠ n_{j+1}` and `n'` is `n` with the two
occupations exchanged. -/
theorem hop_apply (j : Fin N) (n' n : Cfg N) :
    hop j n' n =
      if hj : j.val + 1 < N then
        (if n j ≠ n ⟨j.val + 1, hj⟩ ∧ n' = n ∘ Equiv.swap j ⟨j.val + 1, hj⟩ then 1 else 0)
      else 0 := by
  by_cases hj : j.val + 1 < N
  · rw [hop, dite_eq_left hj, dite_eq_left hj, Matrix.add_apply, cre_mul_ann_apply,
      cre_mul_ann_apply]
    generalize hJ : (⟨j.val + 1, hj⟩ : Fin N) = j'
    have hj' : j'.val = j.val + 1 := by rw [← hJ]
    have hjj' : j ≠ j' := by intro h; rw [h] at hj'; omega
    simp only [cre_ann_cond hjj', cre_ann_cond hjj'.symm, Equiv.swap_comm j' j]
    by_cases hs : n' = n ∘ Equiv.swap j j'
    · subst hs
      have h1 := jwPar_comp_swap n hj'
      have h2 : jwPar (n ∘ Equiv.swap j j') j' = jwPar n j + n j' := by
        rw [jwPar_succ _ hj', h1]; simp
      have h3 : jwPar n j' = jwPar n j + n j := jwPar_succ n hj'
      rcases zmod2_eq_zero_or_one (n j) with h0 | h0 <;>
        rcases zmod2_eq_zero_or_one (n j') with h1' | h1' <;>
        simp [h0, h1', h1, h2, h3, sgn_mul_self]
    · simp [hs]
  · simp [hop, hj]

/-- `dΓ` is additive. -/
private lemma dGamma_add (A B : Matrix (Fin N) (Fin N) ℝ) :
    dGamma (A + B) = dGamma A + dGamma B := by
  simp only [dGamma, Matrix.add_apply, add_smul, Finset.sum_add_distrib]

/-- `dΓ` is homogeneous. -/
private lemma dGamma_smul (c : ℝ) (A : Matrix (Fin N) (Fin N) ℝ) :
    dGamma (c • A) = c • dGamma A := by
  simp only [dGamma, Matrix.smul_apply, smul_eq_mul, mul_smul, Finset.smul_sum]

/-- Second quantization of a diagonal one-body matrix. -/
private lemma dGamma_diagonal (d : Fin N → ℝ) :
    dGamma (diagonal d) = diagonal (fun n => ∑ x, d x * occ (n x)) := by
  have e : ∀ x, ∑ y, diagonal d x y • (cre x * ann y) = d x • diagonal (fun n => occ (n x)) := by
    intro x
    rw [Finset.sum_eq_single x]
    · rw [diagonal_apply_eq, cre_mul_ann_self]
    · intro y _ hy; rw [diagonal_apply_ne _ (Ne.symm hy), zero_smul]
    · simp
  simp only [dGamma, e]
  ext n' n
  simp only [Matrix.sum_apply, Matrix.smul_apply, diagonal_apply, smul_eq_mul]
  split_ifs <;> simp

/-- Summing over the right neighbour `y = x + 1` of `x`. -/
private lemma sum_ite_succ_eq {M : Type*} [AddCommMonoid M] (x : Fin N) (f : Fin N → M) :
    ∑ y : Fin N, (if x.val + 1 = y.val then f y else 0) =
      if h : x.val + 1 < N then f ⟨x.val + 1, h⟩ else 0 := by
  by_cases h : x.val + 1 < N
  · rw [dite_eq_left h, Finset.sum_eq_single ⟨x.val + 1, h⟩]
    · simp
    · intro y _ hy
      rw [ite_eq_right]
      intro h'
      exact absurd (Fin.ext h'.symm) hy
    · simp
  · rw [dite_eq_right h]
    apply Finset.sum_eq_zero
    intro y _
    rw [ite_eq_right]
    intro h'
    have := y.isLt
    omega

/-- Second quantization of the path adjacency matrix is the sum of the hopping terms. -/
private lemma dGamma_pathAdj : dGamma (pathAdj N) = ∑ j, hop j := by
  have e1 : ∀ x y : Fin N, pathAdj N x y • (cre x * ann y) =
      (if x.val + 1 = y.val then cre x * ann y else 0) +
        (if y.val + 1 = x.val then cre x * ann y else 0) := by
    intro x y
    by_cases h1 : x.val + 1 = y.val <;> by_cases h2 : y.val + 1 = x.val <;> simp [pathAdj, h1, h2]
    omega
  have e2 : ∑ x : Fin N, ∑ y : Fin N, (if y.val + 1 = x.val then cre x * ann y else 0) =
      ∑ y : Fin N, if h : y.val + 1 < N then cre ⟨y.val + 1, h⟩ * ann y else 0 := by
    rw [Finset.sum_comm]
    simp only [sum_ite_succ_eq]
  simp only [dGamma, e1, Finset.sum_add_distrib, e2]
  simp only [sum_ite_succ_eq]
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  unfold hop
  split_ifs <;> simp

/-- Explicit second quantization of `h_q`. -/
theorem dGamma_hq (J h : ℝ) (q : Cfg N) :
    dGamma (hq J h q) =
      (-J) • ∑ j, hop j + (2 * h) • diagonal (fun n => ∑ x, sgn (q x) * occ (n x)) := by
  rw [hq, dGamma_add, dGamma_smul, dGamma_smul, dGamma_pathAdj, dGamma_diagonal]

/-! ### The duality bijection -/

/-- `{k ≤ j} = insert j {k < j}`. -/
private lemma filter_le_eq_insert (j : Fin N) :
    univ.filter (fun k => k ≤ j) = insert j (univ.filter (fun k => k < j)) := by
  ext k; simp [le_iff_lt_or_eq, or_comm]

/-- `linkOf n q j` is the strict partial sum plus the `j`-th term. -/
private lemma linkOf_eq_lt_add (n q : Cfg N) (j : Fin N) :
    linkOf n q j = ∑ k ∈ univ.filter (fun k => k < j), (q k + n k) + (q j + n j) := by
  unfold linkOf
  rw [filter_le_eq_insert, Finset.sum_insert (by simp), add_comm]

/-- The label of the link to the left of `j` is the strict partial sum. -/
private lemma prevLink_linkOf (n q : Cfg N) (j : Fin N) :
    prevLink (linkOf n q) j = ∑ k ∈ univ.filter (fun k => k < j), (q k + n k) := by
  unfold prevLink
  split_ifs with hj
  · symm
    apply Finset.sum_eq_zero
    intro k hk
    simp only [mem_filter, mem_univ, true_and, Fin.lt_def] at hk
    omega
  · unfold linkOf
    congr 1
    ext k
    simp only [mem_filter, mem_univ, true_and, Fin.lt_def, Fin.le_def]
    omega

/-- The electric term in dual variables (used before `prevLink_linkOf_add` is available). -/
private lemma prevLink_add_linkOf (n q : Cfg N) (j : Fin N) :
    prevLink (linkOf n q) j + linkOf n q j = q j + n j := by
  rw [prevLink_linkOf, linkOf_eq_lt_add, ← add_assoc, zmod2_add_self, zero_add]

theorem chargeOf_linkOf (n q : Cfg N) : chargeOf n (linkOf n q) = q := by
  funext j
  simp only [chargeOf]
  rw [prevLink_add_linkOf, add_assoc, zmod2_add_self, add_zero]

/-- For fixed occupations, `q ↦ linkOf n q` is injective. -/
private lemma linkOf_injective (n : Cfg N) : Function.Injective (linkOf n) := by
  intro q q' h
  rw [← chargeOf_linkOf n q, h, chargeOf_linkOf]

theorem linkOf_chargeOf (n x : Cfg N) : linkOf n (chargeOf n x) = x := by
  obtain ⟨q, rfl⟩ := Finite.surjective_of_injective (linkOf_injective n) x
  rw [chargeOf_linkOf]

/-- The duality bijection `(n, q) ↦ (n, x(n, q))` between (occupations × charges) and
(occupations × link labels). -/
def dualEquiv : Cfg N × Cfg N ≃ Cfg N × Cfg N where
  toFun a := (a.1, linkOf a.1 a.2)
  invFun b := (b.1, chargeOf b.1 b.2)
  left_inv a := by simp [chargeOf_linkOf]
  right_inv b := by simp [linkOf_chargeOf]

@[simp] theorem dualEquiv_apply (a : Cfg N × Cfg N) : dualEquiv a = (a.1, linkOf a.1 a.2) := rfl

@[simp] theorem dualEquiv_symm_apply (b : Cfg N × Cfg N) :
    dualEquiv.symm b = (b.1, chargeOf b.1 b.2) := rfl

/-- The electric term in dual variables: `σˣ_{b_{j-1}} σˣ_{b_j} = q_j (-1)^{n_j}`. -/
theorem prevLink_linkOf_add (n q : Cfg N) (j : Fin N) :
    prevLink (linkOf n q) j + linkOf n q j = q j + n j := by
  exact prevLink_add_linkOf n q j

/-- `linkOf n` is injective (iff form). -/
private lemma linkOf_inj (n q' q : Cfg N) : linkOf n q' = linkOf n q ↔ q' = q :=
  (linkOf_injective n).eq_iff

/-- Swapping two different occupations adds `δ_j + δ_{j'}` (in `ZMod 2`). -/
private lemma comp_swap_apply (n : Cfg N) {j j' : Fin N} (hjj' : j ≠ j') (hn : n j ≠ n j')
    (m : Fin N) :
    (n ∘ Equiv.swap j j') m = n m + (if m = j then 1 else 0) + (if m = j' then 1 else 0) := by
  rw [zmod2_ne_iff] at hn
  simp only [Function.comp_apply]
  by_cases h1 : m = j
  · subst h1; simp [hjj', hn]
  · by_cases h2 : m = j'
    · subst h2; simp [h1, hn, add_assoc, zmod2_one_add_one]
    · simp [h1, h2, Equiv.swap_apply_of_ne_of_ne h1 h2]

/-- Effect of swapping two different occupations on the dual link labels. -/
private lemma linkOf_comp_swap_apply (n q : Cfg N) {j j' : Fin N} (hjj' : j ≠ j')
    (hn : n j ≠ n j') (k : Fin N) :
    linkOf (n ∘ Equiv.swap j j') q k =
      linkOf n q k + (if j ≤ k then 1 else 0) + (if j' ≤ k then 1 else 0) := by
  simp only [linkOf]
  have e : ∀ m, q m + (n ∘ Equiv.swap j j') m =
      (q m + n m) + (if m = j then 1 else 0) + (if m = j' then 1 else 0) := by
    intro m; rw [comp_swap_apply n hjj' hn m]; ring
  rw [Finset.sum_congr rfl (fun m _ => e m), Finset.sum_add_distrib, Finset.sum_add_distrib,
    Finset.sum_ite_eq', Finset.sum_ite_eq']
  simp

/-- Moving a particle across the bond `(j, j+1)` flips exactly the label of link `b_j`. -/
theorem linkOf_swap (n q : Cfg N) (j : Fin N) (hj : j.val + 1 < N)
    (hn : n j ≠ n ⟨j.val + 1, hj⟩) :
    linkOf (n ∘ Equiv.swap j ⟨j.val + 1, hj⟩) q =
      Function.update (linkOf n q) j (linkOf n q j + 1) := by
  generalize hJ : (⟨j.val + 1, hj⟩ : Fin N) = j' at hn ⊢
  have hj' : j'.val = j.val + 1 := by rw [← hJ]
  have hjj' : j ≠ j' := by intro h; rw [h] at hj'; omega
  funext k
  rw [linkOf_comp_swap_apply n q hjj' hn k]
  by_cases hkj : k = j
  · subst hkj
    have : ¬ j' ≤ k := by rw [Fin.le_def]; omega
    simp [this]
  · rw [Function.update_of_ne hkj]
    have hkj' : k.val ≠ j.val := fun h => hkj (Fin.ext h)
    by_cases hlt : k.val < j.val
    · have h1 : ¬ j ≤ k := by rw [Fin.le_def]; omega
      have h2 : ¬ j' ≤ k := by rw [Fin.le_def]; omega
      simp [h1, h2]
    · have h1 : j ≤ k := by rw [Fin.le_def]; omega
      have h2 : j' ≤ k := by rw [Fin.le_def]; omega
      simp [h1, h2, add_assoc, zmod2_one_add_one]

/-! ### Main statements -/

/-- The hopping term in dual variables: it conserves the charges. -/
private lemma hop_mul_sz (j : Fin N) (n' n q' q : Cfg N) :
    hop j n' n * sz j (linkOf n' q') (linkOf n q) = if q' = q then hop j n' n else 0 := by
  rw [hop_apply]
  by_cases hj : j.val + 1 < N
  · rw [dite_eq_left hj]
    by_cases hc : n j ≠ n ⟨j.val + 1, hj⟩ ∧ n' = n ∘ Equiv.swap j ⟨j.val + 1, hj⟩
    · rw [ite_eq_left hc, one_mul]
      obtain ⟨hne, rfl⟩ := hc
      simp only [sz]
      rw [← linkOf_swap n q j hj hne]
      simp only [linkOf_inj]
    · simp [hc]
  · simp [hj]

/-- The electric term in dual variables: it is diagonal, with entries `sgn (q j + n j)`. -/
private lemma one_mul_elec (j : Fin N) (n' n q' q : Cfg N) :
    (1 : Matrix (Cfg N) (Cfg N) ℝ) n' n * elec j (linkOf n' q') (linkOf n q) =
      if q' = q then (if n' = n then sgn (q j + n j) else 0) else 0 := by
  rw [Matrix.one_apply, elec, diagonal_apply]
  by_cases hn : n' = n
  · subst hn
    simp only [linkOf_inj, ite_true, one_mul]
    split_ifs with hqq
    · subst hqq; rw [prevLink_linkOf_add]
    · rfl
  · simp [hn]

/-- **Exact duality** (paper: Theorem `thm:duality`). In the basis `(n, q)` the LGT Hamiltonian is
the block-diagonal free-fermion Hamiltonian `⊕_q (dΓ(h_q) - h ∑_j q_j)`. -/
theorem duality (J h : ℝ) :
    (H J h).submatrix (dualEquiv (N := N)) dualEquiv = Hdual J h := by
  ext ⟨n', q'⟩ ⟨n, q⟩
  rw [submatrix_apply, dualEquiv_apply, dualEquiv_apply, Hdual, blockDiagonal_apply]
  simp only [H, Matrix.sub_apply, Matrix.smul_apply, Matrix.sum_apply, kronecker_apply,
    smul_eq_mul, hop_mul_sz, one_mul_elec]
  by_cases hqq : q' = q
  · subst hqq
    simp only [ite_true]
    rw [Hsec, dGamma_hq]
    simp only [Matrix.sub_apply, Matrix.add_apply, Matrix.smul_apply, Matrix.sum_apply,
      diagonal_apply, Matrix.one_apply, smul_eq_mul]
    by_cases hn : n' = n
    · subst hn
      have e : ∀ j, sgn (q' j + n' j) = sgn (q' j) - 2 * (sgn (q' j) * occ (n' j)) := by
        intro j; rw [sgn_add, sgn_eq_one_sub_two_occ (n' j)]; ring
      simp only [ite_true, e, Finset.sum_sub_distrib, ← Finset.mul_sum]
      ring
    · simp [hn]
  · simp [hqq]

/-- The Gauss operators become the charge operators `τˣ_j = diag(sgn q_j)`. -/
theorem gauss_dual (j : Fin N) :
    (gauss j).submatrix (dualEquiv (N := N)) dualEquiv = diagonal (fun a => sgn (a.2 j)) := by
  rw [gauss, submatrix_diagonal_equiv]
  congr 1
  funext a
  simp only [Function.comp_apply, dualEquiv_apply]
  rw [prevLink_linkOf_add, sgn_add, mul_left_comm, sgn_mul_self, mul_one]

/-- Operators of the form `diagonal d ⊗ₖ 1` (diagonal in the matter, trivial on the gauge
field) are unchanged by the duality. -/
private lemma diagonal_kron_one_dual (d : Cfg N → ℝ) :
    (diagonal d ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)).submatrix (dualEquiv (N := N)) dualEquiv =
      diagonal d ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) := by
  rw [← diagonal_one, diagonal_kronecker_diagonal, submatrix_diagonal_equiv]
  rfl

/-- The matter number operators are unchanged by the duality. -/
theorem numOp_dual (A : Finset (Fin N)) :
    (numOp A ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)).submatrix (dualEquiv (N := N)) dualEquiv =
      numOp A ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) := by
  exact diagonal_kron_one_dual _

/-- The spectral projections of the matter number operators are unchanged by the duality. -/
theorem numProj_dual (A : Finset (Fin N)) (k : ℕ) :
    (numProj A k ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)).submatrix (dualEquiv (N := N)) dualEquiv =
      numProj A k ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) := by
  exact diagonal_kron_one_dual _

end DFL
