import DFL.Variants
import DFL.RandomWalk
import DFL.CutLocal

/-!
# The ring: Wilson loop and twisted boundary condition (paper: Proposition `prop:variants`(ii))

On the ring (0-indexed sites `0, …, N-1`, links `b_j = (j, j+1 mod N)`), the hopping sum includes
the bond `(N-1, 0)` carried by the link `b_{N-1}`, and the Gauss law at site `0` involves
`σˣ_{b_{N-1}}`. In the dual variables `(n, q)` of the open chain (`dualEquiv`), with the same map,
* `G_0^ring ↦ (-1)^{N_f} ∏_{k ≠ 0} τˣ_k`, `G_j^ring ↦ τˣ_j` (`j ≠ 0`), `W = ∏_j σᶻ_{b_j} ↦ τᶻ_0`;
* the ring bond becomes `-J τᶻ_0 (c†_{N-1} c_0 + h.c.)`.
Hence on the joint eigenspace `{G_j^ring = q_j ∀ j, W = w}` the Hamiltonian acts as
`dΓ(h_{q,w}) - h ∑_j q_j` on the matter states of parity `∑_j q_j`, where `h_{q,w}` is `h_q` with the
bond `(N-1, 0)` added with amplitude `-J w`: the twist is exactly the Wilson loop eigenvalue.

Conventions: `ω : ZMod 2` labels `w = sgn ω`; `q : Cfg N` labels `q_j = sgn (q j)`.
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ} [NeZero N]

/-- The last site `N - 1`; the ring bond `(N - 1, 0)` is carried by the link `b_{N-1}`. -/
def lastSite : Fin N := ⟨N - 1, Nat.sub_lt (NeZero.pos N) Nat.one_pos⟩

/-- The link to the left of site `j` on the ring: `b_{j-1}`, with `b_{-1} := b_{N-1}`. -/
def prevLinkRing (x : Cfg N) (j : Fin N) : ZMod 2 :=
  if j.val = 0 then x lastSite else prevLink x j

/-- Ring electric term `σˣ_{b_{j-1}} σˣ_{b_j}`. -/
def elecRing (j : Fin N) : Matrix (Cfg N) (Cfg N) ℝ :=
  diagonal (fun x => sgn (prevLinkRing x j + x j))

/-- Ring Gauss operator `G_j^ring = (-1)^{n_j} σˣ_{b_{j-1}} σˣ_{b_j}`. -/
def gaussRing (j : Fin N) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  diagonal (fun p => sgn (p.1 j) * sgn (prevLinkRing p.2 j + p.2 j))

/-- Hopping across the ring bond `(N - 1, 0)`. -/
def hopWrap : Matrix (Cfg N) (Cfg N) ℝ :=
  cre lastSite * ann 0 + cre 0 * ann lastSite

/-- The ring Hamiltonian
`H^ring = -J ∑_{j} σᶻ_{b_j} (c†_j c_{j+1} + h.c.) - h ∑_j σˣ_{b_{j-1}} σˣ_{b_j}` (indices mod `N`). -/
def Hring (J h : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  (-J) • (∑ j, hop j ⊗ₖ sz j + hopWrap ⊗ₖ sz lastSite) -
    h • ∑ j, (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ elecRing j

/-- `∏_j σᶻ_{b_j}` in the `σˣ` basis: flips every link label. -/
def flipAll : Matrix (Cfg N) (Cfg N) ℝ := fun x' x => if x' = x + 1 then 1 else 0

/-- The Wilson loop `W = ∏_j σᶻ_{b_j}`. -/
def wilson : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flipAll

/-- `τᶻ_0` in the dual variables: flips the charge label `q_0`. -/
def flip0 : Matrix (Cfg N) (Cfg N) ℝ :=
  fun q' q => if q' = Function.update q 0 (q 0 + 1) then 1 else 0

/-- The ring Hamiltonian in the dual variables `(n, q)`. -/
def HringDual (J h : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  (-J) • (∑ j, hop j ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) + hopWrap ⊗ₖ flip0) -
    h • (diagonal (fun p : Cfg N × Cfg N => sgn (p.1 0 + ∑ k, p.1 k + ∑ k ∈ univ.erase 0, p.2 k)) +
      ∑ j ∈ univ.erase (0 : Fin N), diagonal (fun p : Cfg N × Cfg N => sgn (p.1 j + p.2 j)))

/-- The adjacency matrix of the ring bond `(N - 1, 0)`. -/
def wrapAdj : Matrix (Fin N) (Fin N) ℝ :=
  fun x y => if (x = lastSite ∧ y = 0) ∨ (x = 0 ∧ y = lastSite) then 1 else 0

/-- The twisted one-body Hamiltonian `h_{q,w} = h_q - J w (|N-1⟩⟨0| + |0⟩⟨N-1|)`, `w = sgn ω`. -/
def hqw (J h : ℝ) (q : Cfg N) (ω : ZMod 2) : Matrix (Fin N) (Fin N) ℝ :=
  hq J h q + (-J * sgn ω) • wrapAdj


/-! ### Auxiliary lemmas (namespace `DFL.RingAux`) -/

namespace RingAux

/-- Every element of `ZMod 2` is `0` or `1`. -/
lemma zmod2_cases (a : ZMod 2) : a = 0 ∨ a = 1 := by
  revert a; decide

/-- `a + a = 0` in `ZMod 2`. -/
lemma zmod2_add_self (a : ZMod 2) : a + a = 0 := by
  revert a; decide

/-- `sgn` turns addition in `ZMod 2` into multiplication. -/
lemma sgn_add (a b : ZMod 2) : sgn (a + b) = sgn a * sgn b := RandomWalkAux.sgn_add a b

/-- `sgn a ^ 2 = 1`. -/
lemma sgn_mul_self (a : ZMod 2) : sgn a * sgn a = 1 := by
  rcases zmod2_cases a with rfl | rfl <;> simp [sgn]

/-- `(-1)^a = 1 - 2a` for `a ∈ {0, 1}`. -/
lemma sgn_eq_one_sub_two_occ (a : ZMod 2) : sgn a = 1 - 2 * occ a := by
  rcases zmod2_cases a with rfl | rfl <;> norm_num [sgn, occ]

/-- A product of signs is the sign of the sum. -/
lemma prod_sgn {ι : Type*} (s : Finset ι) (f : ι → ZMod 2) :
    ∏ i ∈ s, sgn (f i) = sgn (∑ i ∈ s, f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [sgn]
  | insert a s ha ih => rw [Finset.prod_insert ha, Finset.sum_insert ha, ih, sgn_add]

/-- Splitting off site `0` from a sum. -/
lemma sum_eq_zero_add (q : Cfg N) : ∑ k, q k = q 0 + ∑ k ∈ univ.erase 0, q k :=
  (Finset.add_sum_erase _ _ (mem_univ 0)).symm

/-- Two configurations that agree off site `0` and have the same total are equal. -/
lemma eq_of_eqOff_of_sum {q' q : Cfg N} (h : ∀ k, k ≠ 0 → q' k = q k)
    (hs : ∑ k, q' k = ∑ k, q k) : q' = q := by
  funext k
  by_cases hk : k = 0
  · subst hk
    rw [sum_eq_zero_add q', sum_eq_zero_add q,
      Finset.sum_congr rfl (fun k hk => h k (Finset.ne_of_mem_erase hk))] at hs
    exact add_right_cancel hs
  · exact h k hk

/-- The value of `lastSite`. -/
lemma lastSite_val : (lastSite : Fin N).val = N - 1 := rfl

/-- Every site is `≤ lastSite`. -/
lemma le_lastSite (k : Fin N) : k ≤ lastSite := by
  rw [Fin.le_def, lastSite_val]; have := k.isLt; omega

/-- `lastSite ≤ k` iff `k = lastSite`. -/
lemma lastSite_le_iff (k : Fin N) : lastSite ≤ k ↔ k = lastSite := by
  rw [Fin.le_def, Fin.ext_iff, lastSite_val]; have := k.isLt; omega

/-- For `N ≥ 2`, the last site is not the site `0`. -/
lemma lastSite_ne_zero (hN : 2 ≤ N) : (lastSite : Fin N) ≠ 0 := by
  intro h
  have := congrArg Fin.val h
  rw [lastSite_val, Fin.val_zero] at this
  omega

/-- Flipping the charge label of site `0`. -/
def flipQ (q : Cfg N) : Cfg N := Function.update q 0 (q 0 + 1)

/-- Entries of `flipQ q`. -/
lemma flipQ_apply (q : Cfg N) (k : Fin N) : flipQ q k = q k + if k = 0 then 1 else 0 := by
  by_cases hk : k = 0
  · subst hk; simp [flipQ]
  · simp [flipQ, hk]

/-- `flipQ` is an involution. -/
lemma flipQ_flipQ (q : Cfg N) : flipQ (flipQ q) = q := by
  funext k
  rw [flipQ_apply, flipQ_apply, add_assoc]
  split_ifs <;> simp [zmod2_add_self]

/-- `q' = flipQ q ↔ q = flipQ q'`. -/
lemma eq_flipQ_comm (q' q : Cfg N) : q' = flipQ q ↔ q = flipQ q' := by
  constructor
  · rintro rfl; rw [flipQ_flipQ]
  · rintro rfl; rw [flipQ_flipQ]

/-- `flipQ` does not change the labels off site `0`. -/
lemma flipQ_ne (q : Cfg N) {k : Fin N} (hk : k ≠ 0) : flipQ q k = q k := by
  rw [flipQ_apply, ite_eq_right hk, add_zero]

/-- `flip0` in terms of `flipQ`. -/
lemma flip0_apply (q' q : Cfg N) : flip0 q' q = if q' = flipQ q then 1 else 0 := rfl

/-- Flipping `q_0` flips every dual link label. -/
lemma linkOf_flipQ (n q : Cfg N) : linkOf n (flipQ q) = linkOf n q + 1 := by
  funext k
  simp only [linkOf, Pi.add_apply, Pi.one_apply, flipQ_apply]
  have e : ∀ m : Fin N, q m + (if m = 0 then 1 else 0) + n m =
      (q m + n m) + (if m = 0 then 1 else 0) := fun m => by ring
  rw [Finset.sum_congr rfl (fun m _ => e m), Finset.sum_add_distrib, Finset.sum_ite_eq']
  simp

/-- The last dual link label is the total parity. -/
lemma linkOf_lastSite (n q : Cfg N) : linkOf n q lastSite = ∑ k, q k + ∑ k, n k := by
  simp only [linkOf]
  rw [Finset.filter_true_of_mem (fun k _ => le_lastSite k), Finset.sum_add_distrib]

/-- The first dual link label. -/
lemma linkOf_zero (n q : Cfg N) : linkOf n q 0 = q 0 + n 0 := by
  simp only [linkOf]
  have : univ.filter (fun k : Fin N => k ≤ 0) = {0} := by
    ext k; simp [nonpos_iff_eq_zero]
  rw [this, Finset.sum_singleton]

/-- The ring electric term in dual variables. -/
lemma prevLinkRing_linkOf_add (n q : Cfg N) (j : Fin N) :
    prevLinkRing (linkOf n q) j + linkOf n q j =
      if j = 0 then n 0 + ∑ k, n k + ∑ k ∈ univ.erase 0, q k else q j + n j := by
  by_cases hj : j = 0
  · subst hj
    rw [ite_eq_left rfl]
    simp only [prevLinkRing, Fin.val_zero, ite_true]
    rw [linkOf_lastSite, linkOf_zero, sum_eq_zero_add q]
    linear_combination zmod2_add_self (q 0)
  · rw [ite_eq_right hj]
    have hj' : j.val ≠ 0 := fun h => hj (Fin.ext (by rw [h, Fin.val_zero]))
    simp only [prevLinkRing, hj', ite_false]
    exact prevLink_linkOf_add n q j

omit [NeZero N] in
/-- Support of `c_a† c_b`: the occupations change exactly at `a` and `b` (in `ZMod 2`). -/
lemma cre_mul_ann_support {a b : Fin N} {n' n : Cfg N} (h : (cre a * ann b) n' n ≠ 0)
    (m : Fin N) : n' m = n m + (if m = a then 1 else 0) + (if m = b then 1 else 0) := by
  simp only [cre, Matrix.mul_apply, Matrix.transpose_apply, ann] at h
  obtain ⟨p, -, hp⟩ := Finset.exists_ne_zero_of_sum_ne_zero h
  by_cases h1 : n' a = 1 ∧ p = Function.update n' a 0
  · by_cases h2 : n b = 1 ∧ p = Function.update n b 0
    · have hu : Function.update n' a 0 = Function.update n b 0 := h1.2.symm.trans h2.2
      have hm := congrFun hu m
      by_cases hma : m = a <;> by_cases hmb : m = b
      · subst hma; subst hmb
        rw [h1.1, h2.1, ite_eq_left rfl]; decide
      · subst hma
        rw [Function.update_self, Function.update_of_ne hmb] at hm
        rw [h1.1, ← hm, ite_eq_left rfl, ite_eq_right hmb]; decide
      · subst hmb
        rw [Function.update_self, Function.update_of_ne hma] at hm
        rw [hm, h2.1, ite_eq_right hma, ite_eq_left rfl]; decide
      · rw [Function.update_of_ne hma, Function.update_of_ne hmb] at hm
        rw [hm, ite_eq_right hma, ite_eq_right hmb, add_zero, add_zero]
    · simp [h2] at hp
  · simp [h1] at hp

/-- Support of the ring bond hopping. -/
lemma hopWrap_support {n' n : Cfg N} (h : hopWrap n' n ≠ 0) (m : Fin N) :
    n' m = n m + (if m = 0 then 1 else 0) + (if m = lastSite then 1 else 0) := by
  have h' : (cre lastSite * ann 0 : Matrix (Cfg N) (Cfg N) ℝ) n' n ≠ 0 ∨
      (cre 0 * ann lastSite : Matrix (Cfg N) (Cfg N) ℝ) n' n ≠ 0 := by
    by_contra hc
    rw [not_or, not_not, not_not] at hc
    exact h (by simp [hopWrap, Matrix.add_apply, hc.1, hc.2])
  rcases h' with h' | h'
  · rw [cre_mul_ann_support h' m]; ring
  · rw [cre_mul_ann_support h' m]

/-- The ring bond hopping conserves the particle number parity. -/
lemma sum_of_hopWrap {n' n : Cfg N} (h : hopWrap n' n ≠ 0) : ∑ k, n' k = ∑ k, n k := by
  rw [Finset.sum_congr rfl (fun m _ => hopWrap_support h m), Finset.sum_add_distrib,
    Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.sum_ite_eq', ite_eq_left (mem_univ _),
    ite_eq_left (mem_univ _), add_assoc, zmod2_add_self, add_zero]

omit [NeZero N] in
/-- The open-chain hopping conserves the particle number parity. -/
lemma sum_of_hop {j : Fin N} {n' n : Cfg N} (h : hop j n' n ≠ 0) : ∑ k, n' k = ∑ k, n k := by
  rw [hop_apply] at h
  split_ifs at h with hj hc
  · rw [hc.2]
    exact Equiv.sum_comp (Equiv.swap j ⟨j.val + 1, hj⟩) n
  · exact absurd rfl h
  · exact absurd rfl h

/-- Moving a particle across the ring bond flips every dual link label except the last one. -/
lemma linkOf_of_hopWrap {n' n : Cfg N} (h : hopWrap n' n ≠ 0) (q : Cfg N) (k : Fin N) :
    linkOf n' q k = linkOf n q k + 1 + (if k = lastSite then 1 else 0) := by
  simp only [linkOf]
  have e : ∀ m : Fin N, q m + n' m =
      (q m + n m) + (if m = 0 then 1 else 0) + (if m = lastSite then 1 else 0) := by
    intro m; rw [hopWrap_support h m]; ring
  rw [Finset.sum_congr rfl (fun m _ => e m), Finset.sum_add_distrib, Finset.sum_add_distrib,
    Finset.sum_ite_eq', Finset.sum_ite_eq']
  simp only [Finset.mem_filter, Finset.mem_univ, true_and, Fin.zero_le, ite_true,
    lastSite_le_iff]

/-- The ring bond in dual variables: `σᶻ_{b_{N-1}}` becomes `τᶻ_0`. -/
lemma sz_lastSite_linkOf {n' n : Cfg N} (h : hopWrap n' n ≠ 0) (q' q : Cfg N) :
    sz lastSite (linkOf n' q') (linkOf n q) = flip0 q' q := by
  have key : Function.update (linkOf n q) lastSite (linkOf n q lastSite + 1) =
      linkOf n' (flipQ q) := by
    funext k
    rw [linkOf_flipQ, Pi.add_apply, Pi.one_apply, linkOf_of_hopWrap h q k]
    by_cases hk : k = lastSite
    · subst hk
      rw [Function.update_self, ite_eq_left rfl, add_assoc, add_assoc, zmod2_add_self, add_zero]
    · rw [Function.update_of_ne hk, ite_eq_right hk, add_zero, add_assoc, zmod2_add_self, add_zero]
  simp only [sz, key, RandomWalkAux.linkOf_inj, flip0_apply]

/-- The ring electric term in dual variables: diagonal. -/
lemma one_mul_elecRing (j : Fin N) (n' n q' q : Cfg N) :
    (1 : Matrix (Cfg N) (Cfg N) ℝ) n' n * elecRing j (linkOf n' q') (linkOf n q) =
      if n' = n ∧ q' = q then sgn (prevLinkRing (linkOf n q) j + linkOf n q j) else 0 := by
  rw [Matrix.one_apply, elecRing, diagonal_apply]
  by_cases hn : n' = n
  · subst hn
    by_cases hq : q' = q
    · subst hq; simp
    · simp [hq, RandomWalkAux.linkOf_inj]
  · simp [hn]

/-- The total ring electric energy in dual variables. -/
lemma sum_sgn_elecRing_dual (n q : Cfg N) :
    ∑ j, sgn (prevLinkRing (linkOf n q) j + linkOf n q j) =
      sgn (n 0 + ∑ k, n k + ∑ k ∈ univ.erase 0, q k) + ∑ j ∈ univ.erase 0, sgn (n j + q j) := by
  rw [← Finset.add_sum_erase _ _ (mem_univ 0), prevLinkRing_linkOf_add, ite_eq_left rfl]
  congr 1
  refine Finset.sum_congr rfl (fun j hj => ?_)
  rw [prevLinkRing_linkOf_add, ite_eq_right (Finset.ne_of_mem_erase hj), add_comm]

end RingAux

/-! ### The dual frame -/

/-- **Ring duality** (paper: Proposition `prop:variants`(ii)): the same map `U` that dualizes the
open chain turns the ring Hamiltonian into `HringDual`. -/
theorem ring_duality (J h : ℝ) :
    (Hring J h).submatrix (dualEquiv (N := N)) dualEquiv = HringDual J h := by
  ext ⟨n', q'⟩ ⟨n, q⟩
  rw [submatrix_apply, dualEquiv_apply, dualEquiv_apply]
  have hw : hopWrap n' n * sz lastSite (linkOf n' q') (linkOf n q) = hopWrap n' n * flip0 q' q := by
    by_cases hz : hopWrap n' n = 0
    · simp [hz]
    · rw [RingAux.sz_lastSite_linkOf hz]
  simp only [Hring, HringDual, Matrix.sub_apply, Matrix.smul_apply, Matrix.add_apply,
    Matrix.sum_apply, kronecker_apply, smul_eq_mul, RandomWalkAux.hop_mul_sz,
    RingAux.one_mul_elecRing, diagonal_apply, hw]
  simp only [Matrix.one_apply, Prod.mk.injEq]
  by_cases hq : q' = q
  · subst hq
    by_cases hn : n' = n
    · subst hn
      simp only [and_self, ite_true, mul_one, RingAux.sum_sgn_elecRing_dual]
    · simp [hn]
  · simp [hq]

/-- `U† G_0^ring U = (-1)^{N_f} ∏_{k ≠ 0} τˣ_k`. -/
theorem gaussRing_dual_zero :
    (gaussRing (0 : Fin N)).submatrix (dualEquiv (N := N)) dualEquiv =
      diagonal (fun p : Cfg N × Cfg N => sgn (∑ k, p.1 k + ∑ k ∈ univ.erase 0, p.2 k)) := by
  rw [gaussRing, submatrix_diagonal_equiv]
  congr 1
  funext a
  simp only [Function.comp_apply, dualEquiv_apply]
  rw [RingAux.prevLinkRing_linkOf_add, ite_eq_left rfl, ← RingAux.sgn_add]
  congr 1
  linear_combination RingAux.zmod2_add_self (a.1 0)

/-- `U† G_j^ring U = τˣ_j` for `j ≠ 0`. -/
theorem gaussRing_dual (j : Fin N) (hj : j ≠ 0) :
    (gaussRing j).submatrix (dualEquiv (N := N)) dualEquiv =
      diagonal (fun p : Cfg N × Cfg N => sgn (p.2 j)) := by
  rw [gaussRing, submatrix_diagonal_equiv]
  congr 1
  funext a
  simp only [Function.comp_apply, dualEquiv_apply]
  rw [RingAux.prevLinkRing_linkOf_add, ite_eq_right hj, RingAux.sgn_add, mul_left_comm,
    RingAux.sgn_mul_self, mul_one]

/-- `U† W U = τᶻ_0`. -/
theorem wilson_dual :
    (wilson (N := N)).submatrix (dualEquiv (N := N)) dualEquiv =
      (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flip0 := by
  ext ⟨n', q'⟩ ⟨n, q⟩
  rw [submatrix_apply, dualEquiv_apply, dualEquiv_apply]
  simp only [wilson, kronecker_apply, Matrix.one_apply, flipAll, RingAux.flip0_apply]
  by_cases hn : n' = n
  · subst hn
    simp only [ite_true, one_mul]
    rw [← RingAux.linkOf_flipQ]
    simp only [RandomWalkAux.linkOf_inj]
  · simp [hn]

namespace RingAux

/-- A matrix commutes with a diagonal matrix whose entries are constant along its support. -/
lemma commute_diagonal_of_support {ι : Type*} [Fintype ι] [DecidableEq ι] (M : Matrix ι ι ℝ)
    (g : ι → ℝ) (h : ∀ a b, M a b ≠ 0 → g a = g b) : Commute M (diagonal g) := by
  show M * diagonal g = diagonal g * M
  ext a b
  rw [mul_diagonal, diagonal_mul]
  by_cases hab : M a b = 0
  · simp [hab]
  · rw [h a b hab, mul_comm]

/-- Diagonal matrices commute. -/
lemma commute_diagonal_diagonal {ι : Type*} [Fintype ι] [DecidableEq ι] (f g : ι → ℝ) :
    Commute (diagonal f) (diagonal g) := by
  show diagonal f * diagonal g = diagonal g * diagonal f
  rw [diagonal_mul_diagonal, diagonal_mul_diagonal]
  congr 1
  funext i
  ring

omit [NeZero N] in
/-- Commutation can be checked after relabelling the basis by `dualEquiv`. -/
lemma comm_of_dual (A B : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ)
    (h : Commute (A.submatrix dualEquiv dualEquiv) (B.submatrix dualEquiv dualEquiv)) :
    A * B = B * A := by
  have h' := h.eq
  rw [submatrix_mul_equiv, submatrix_mul_equiv] at h'
  have h'' := congrArg (fun M => M.submatrix dualEquiv.symm dualEquiv.symm) h'
  simpa [submatrix_submatrix] using h''

/-- `HringDual` commutes with every diagonal matrix depending only on the particle-number parity
and on the charges off site `0`. -/
lemma HringDual_comm_diagonal (J h : ℝ) (g : Cfg N × Cfg N → ℝ)
    (h1 : ∀ m' m r, ∑ k, m' k = ∑ k, m k → g (m', r) = g (m, r))
    (h2 : ∀ m' m r, ∑ k, m' k = ∑ k, m k → g (m', flipQ r) = g (m, r)) :
    Commute (HringDual J h) (diagonal g) := by
  unfold HringDual
  refine Commute.sub_left (Commute.smul_left (Commute.add_left
    (Commute.sum_left _ _ _ fun j _ => ?_) ?_) _) (Commute.smul_left (Commute.add_left
    (commute_diagonal_diagonal _ _) (Commute.sum_left _ _ _ fun j _ =>
      commute_diagonal_diagonal _ _)) _)
  · apply commute_diagonal_of_support
    rintro ⟨m', r'⟩ ⟨m, r⟩ hne
    rw [kronecker_apply, Matrix.one_apply] at hne
    by_cases hr : r' = r
    · subst hr
      have hh : hop j m' m ≠ 0 := by intro h0; apply hne; simp [h0]
      exact h1 m' m r' (sum_of_hop hh)
    · exact absurd (by simp [hr]) hne
  · apply commute_diagonal_of_support
    rintro ⟨m', r'⟩ ⟨m, r⟩ hne
    rw [kronecker_apply, flip0_apply] at hne
    by_cases hr : r' = flipQ r
    · subst hr
      have hh : hopWrap m' m ≠ 0 := by intro h0; apply hne; simp [h0]
      exact h2 m' m r (sum_of_hopWrap hh)
    · exact absurd (by simp [hr]) hne

/-- A diagonal matrix not depending on `q_0` commutes with `1 ⊗ τᶻ_0`. -/
lemma commute_diagonal_flip0 (g : Cfg N × Cfg N → ℝ) (hg : ∀ m r, g (m, flipQ r) = g (m, r)) :
    Commute (diagonal g) ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ flip0) := by
  show diagonal g * _ = _ * diagonal g
  ext ⟨m', r'⟩ ⟨m, r⟩
  rw [diagonal_mul, mul_diagonal, kronecker_apply, Matrix.one_apply, flip0_apply]
  by_cases hm : m' = m
  · subst hm
    by_cases hr : r' = flipQ r
    · subst hr; simp [hg]
    · simp [hr]
  · simp [hm]

/-- `∑_{k ≠ 0} q_k` does not depend on `q_0`. -/
lemma sum_erase_flipQ (r : Cfg N) :
    ∑ k ∈ univ.erase 0, flipQ r k = ∑ k ∈ univ.erase 0, r k :=
  Finset.sum_congr rfl (fun _ hk => flipQ_ne r (Finset.ne_of_mem_erase hk))

/-- The diagonal entries of the ring Gauss operators in the dual frame. -/
lemma gaussRing_dualEquiv_apply (j : Fin N) (a : Cfg N × Cfg N) :
    gaussRing j (dualEquiv a) (dualEquiv a) =
      if j = 0 then sgn (∑ k, a.1 k + ∑ k ∈ univ.erase 0, a.2 k) else sgn (a.2 j) := by
  by_cases hj : j = 0
  · subst hj
    have := congrFun (congrFun (gaussRing_dual_zero (N := N)) a) a
    rw [submatrix_apply, diagonal_apply_eq] at this
    rw [this, ite_eq_left rfl]
  · have := congrFun (congrFun (gaussRing_dual j hj) a) a
    rw [submatrix_apply, diagonal_apply_eq] at this
    rw [this, ite_eq_right hj]

end RingAux

/-- The Gauss operators commute with the ring Hamiltonian. -/
theorem Hring_comm_gaussRing (J h : ℝ) (j : Fin N) :
    Hring J h * gaussRing j = gaussRing j * Hring J h := by
  apply RingAux.comm_of_dual
  rw [ring_duality]
  by_cases hj : j = 0
  · subst hj
    rw [gaussRing_dual_zero]
    apply RingAux.HringDual_comm_diagonal
    · intro m' m r hs; simp only [hs]
    · intro m' m r hs; simp only [hs, RingAux.sum_erase_flipQ]
  · rw [gaussRing_dual j hj]
    apply RingAux.HringDual_comm_diagonal
    · intro m' m r _; rfl
    · intro m' m r _; simp only [RingAux.flipQ_ne r hj]

/-- The Wilson loop commutes with the ring Hamiltonian. -/
theorem Hring_comm_wilson (J h : ℝ) :
    Hring J h * wilson (N := N) = wilson (N := N) * Hring J h := by
  apply RingAux.comm_of_dual
  rw [ring_duality, wilson_dual]
  unfold HringDual
  refine Commute.sub_left (Commute.smul_left (Commute.add_left
    (Commute.sum_left _ _ _ fun j _ => ?_) ?_) _) (Commute.smul_left (Commute.add_left
    (RingAux.commute_diagonal_flip0 _ ?_) (Commute.sum_left _ _ _ fun j hj =>
      RingAux.commute_diagonal_flip0 _ ?_)) _)
  · show _ * _ = _ * _
    rw [← mul_kronecker_mul, ← mul_kronecker_mul, Matrix.mul_one, Matrix.one_mul, Matrix.one_mul,
      Matrix.mul_one]
  · show _ * _ = _ * _
    rw [← mul_kronecker_mul, ← mul_kronecker_mul, Matrix.mul_one, Matrix.one_mul]
  · intro m r; simp only [RingAux.sum_erase_flipQ]
  · intro m r; simp only [RingAux.flipQ_ne r (Finset.ne_of_mem_erase hj)]

/-- `∏_j G_j^ring = (-1)^{N_f}` (the Gauss operators are diagonal; this is the product of their
entries). -/
theorem gaussRing_prod (p : Cfg N × Cfg N) :
    ∏ j, gaussRing j p p = sgn (∑ k, p.1 k) := by
  obtain ⟨a, rfl⟩ := dualEquiv.surjective p
  have hj : ∀ j ∈ univ.erase (0 : Fin N), gaussRing j (dualEquiv a) (dualEquiv a) = sgn (a.2 j) :=
    fun j hj => by
      rw [RingAux.gaussRing_dualEquiv_apply, ite_eq_right (Finset.ne_of_mem_erase hj)]
  rw [← Finset.mul_prod_erase _ _ (mem_univ 0), RingAux.gaussRing_dualEquiv_apply, ite_eq_left rfl,
    Finset.prod_congr rfl hj, RingAux.prod_sgn, ← RingAux.sgn_add, dualEquiv_apply]
  congr 1
  linear_combination RingAux.zmod2_add_self (∑ k ∈ univ.erase 0, a.2 k)

/-! ### Joint eigenvectors and the twisted sector Hamiltonians -/

/-- The joint eigenvector of all `G_j^ring` (eigenvalues `sgn (q j)`) and of `W` (eigenvalue
`sgn ω`) with matter `|n⟩`, for `∑ n = ∑ q`. In the dual frame it is
`|n⟩ ⊗ (|q⁰⟩ + sgn ω |q¹⟩)/√2`, where `q⁰, q¹` agree with `q` off site `0` and have `q_0 = 0, 1`. -/
def ringVec (q : Cfg N) (ω : ZMod 2) (n : Cfg N) : Cfg N × Cfg N → ℂ :=
  fun p => if p.1 = n ∧ ∀ k, k ≠ 0 → chargeOf p.1 p.2 k = q k then
    ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ * ((sgn (ω * chargeOf p.1 p.2 0) : ℝ) : ℂ) else 0


namespace RingAux

/-- The charge part of `ringVec` in the dual frame: `(|q⁰⟩ + sgn ω |q¹⟩)/√2` (real-valued). -/
def chVec (q : Cfg N) (ω : ZMod 2) (r : Cfg N) : ℝ :=
  if ∀ k, k ≠ 0 → r k = q k then (Real.sqrt 2)⁻¹ * sgn (ω * r 0) else 0

/-- `ringVec` in the dual frame is `|n⟩ ⊗ chVec`. -/
lemma ringVec_dualEquiv (q : Cfg N) (ω : ZMod 2) (n : Cfg N) (a : Cfg N × Cfg N) :
    ringVec q ω n (dualEquiv a) = ((if a.1 = n then chVec q ω a.2 else 0 : ℝ) : ℂ) := by
  simp only [ringVec, dualEquiv_apply, chargeOf_linkOf, chVec]
  by_cases h1 : a.1 = n
  · by_cases h2 : ∀ k, k ≠ 0 → a.2 k = q k
    · rw [ite_eq_left ⟨h1, h2⟩, ite_eq_left h1, ite_eq_left h2, Complex.ofReal_mul,
        Complex.ofReal_inv]
    · rw [ite_eq_right (fun h => h2 h.2), ite_eq_left h1, ite_eq_right h2, Complex.ofReal_zero]
  · rw [ite_eq_right (fun h => h1 h.1), ite_eq_right h1, Complex.ofReal_zero]

/-- Summing over the charge configurations that agree with `q` off site `0`. -/
lemma sum_eqOff {M : Type*} [AddCommMonoid M] (q : Cfg N) (F : Cfg N → M) :
    ∑ r : Cfg N, (if ∀ k, k ≠ 0 → r k = q k then F r else 0) =
      F (Function.update q 0 0) + F (Function.update q 0 1) := by
  rw [← Finset.sum_filter]
  have hS : univ.filter (fun r : Cfg N => ∀ k, k ≠ 0 → r k = q k) =
      {Function.update q 0 0, Function.update q 0 1} := by
    ext r
    simp only [mem_filter, mem_univ, true_and, mem_insert, mem_singleton]
    constructor
    · intro hr
      rcases zmod2_cases (r 0) with h0 | h0
      · left
        funext k
        by_cases hk : k = 0
        · subst hk; simp [h0]
        · simp [hk, hr k hk]
      · right
        funext k
        by_cases hk : k = 0
        · subst hk; simp [h0]
        · simp [hk, hr k hk]
    · rintro (rfl | rfl) <;> intro k hk <;> simp [hk]
  rw [hS, Finset.sum_pair]
  intro heq
  have := congrFun heq 0
  simp at this

/-- `(√2)⁻¹ (√2)⁻¹ = 1/2`. -/
lemma inv_sqrt_two_mul_self : (Real.sqrt 2)⁻¹ * (Real.sqrt 2)⁻¹ = (1 / 2 : ℝ) := by
  rw [← mul_inv, Real.mul_self_sqrt (by norm_num : (0 : ℝ) ≤ 2)]
  norm_num

/-- The charge parts are orthonormal up to the label `q_0`. -/
lemma sum_chVec_mul (q' q : Cfg N) (ω' ω : ZMod 2) :
    ∑ r, chVec q' ω' r * chVec q ω r =
      if (∀ k, k ≠ 0 → q' k = q k) ∧ ω' = ω then 1 else 0 := by
  by_cases hq : ∀ k, k ≠ 0 → q' k = q k
  · have e : ∀ r, chVec q' ω' r * chVec q ω r =
        if ∀ k, k ≠ 0 → r k = q k then (1 / 2) * (sgn (ω' * r 0) * sgn (ω * r 0)) else 0 := by
      intro r
      simp only [chVec]
      by_cases hr : ∀ k, k ≠ 0 → r k = q k
      · have hr' : ∀ k, k ≠ 0 → r k = q' k := fun k hk => (hr k hk).trans (hq k hk).symm
        rw [ite_eq_left hr', ite_eq_left hr, ite_eq_left hr]
        linear_combination (sgn (ω' * r 0) * sgn (ω * r 0)) * inv_sqrt_two_mul_self
      · rw [ite_eq_right hr, ite_eq_right hr, mul_zero]
    rw [Finset.sum_congr rfl (fun r _ => e r), sum_eqOff]
    simp only [Function.update_self, mul_zero, mul_one]
    by_cases hω : ω' = ω
    · rw [ite_eq_left ⟨hq, hω⟩]
      subst hω
      rcases zmod2_cases ω' with rfl | rfl <;> norm_num [sgn]
    · rw [ite_eq_right (fun h => hω h.2)]
      rcases zmod2_cases ω with rfl | rfl <;> rcases zmod2_cases ω' with rfl | rfl <;>
        simp_all [sgn]
  · rw [ite_eq_right (fun h => hq h.1)]
    apply Finset.sum_eq_zero
    intro r _
    simp only [chVec]
    by_cases h1 : ∀ k, k ≠ 0 → r k = q' k
    · by_cases h2 : ∀ k, k ≠ 0 → r k = q k
      · exact absurd (fun k hk => (h1 k hk).symm.trans (h2 k hk)) hq
      · rw [ite_eq_right h2, mul_zero]
    · rw [ite_eq_right h1, zero_mul]

/-- The overlaps of the dual-frame vectors `|n⟩ ⊗ chVec`. -/
lemma sum_ringVecR (q' q : Cfg N) (ω' ω : ZMod 2) (n' n : Cfg N) :
    ∑ a : Cfg N × Cfg N, (if a.1 = n' then chVec q' ω' a.2 else 0) *
        (if a.1 = n then chVec q ω a.2 else 0) =
      if n' = n ∧ (∀ k, k ≠ 0 → q' k = q k) ∧ ω' = ω then 1 else 0 := by
  simp only [Fintype.sum_prod_type]
  by_cases hn : n' = n
  · subst hn
    rw [Finset.sum_eq_single n']
    · simp only [ite_true, true_and]
      exact sum_chVec_mul q' q ω' ω
    · intro m _ hm
      simp [hm]
    · simp
  · rw [ite_eq_right (fun h => hn h.1)]
    apply Finset.sum_eq_zero
    intro m _
    apply Finset.sum_eq_zero
    intro r _
    by_cases h1 : m = n'
    · subst h1; simp [hn]
    · simp [h1]

/-- Given the total charge, `∑ n + ∑_{k ≠ 0} r_k = q_0` on the support of `ringVec`. -/
lemma gauss0_sign {q n r : Cfg N} (hn : ∑ k, n k = ∑ k, q k) (hr : ∀ k, k ≠ 0 → r k = q k) :
    ∑ k, n k + ∑ k ∈ univ.erase 0, r k = q 0 := by
  rw [hn, sum_eq_zero_add q,
    Finset.sum_congr rfl (fun k hk => hr k (Finset.ne_of_mem_erase hk)), add_assoc,
    zmod2_add_self, add_zero]

omit [NeZero N] in
/-- Matrix-vector products in the dual frame. -/
lemma cplx_mulVec_dualEquiv (M : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ)
    (v : Cfg N × Cfg N → ℂ) (a : Cfg N × Cfg N) :
    (cplx M *ᵥ v) (dualEquiv a) =
      ∑ a', ((M.submatrix dualEquiv dualEquiv a a' : ℝ) : ℂ) * v (dualEquiv a') := by
  simp only [mulVec, dotProduct, cplx, Matrix.map_apply, submatrix_apply]
  exact (Equiv.sum_comp dualEquiv (fun p => ((M (dualEquiv a) p : ℝ) : ℂ) * v p)).symm

/-- In the dual frame the Wilson loop flips `q_0`. -/
lemma wilson_mulVec_dualEquiv (v : Cfg N × Cfg N → ℂ) (m r : Cfg N) :
    (cplx wilson *ᵥ v) (dualEquiv (m, r)) = v (dualEquiv (m, flipQ r)) := by
  rw [cplx_mulVec_dualEquiv, wilson_dual, Finset.sum_eq_single (m, flipQ r)]
  · simp [flip0_apply, flipQ_flipQ]
  · rintro ⟨m', r'⟩ _ hne
    rw [kronecker_apply, Matrix.one_apply, flip0_apply]
    by_cases h1 : m = m'
    · by_cases h2 : r = flipQ r'
      · exfalso
        apply hne
        rw [← h1, (eq_flipQ_comm _ _).1 h2]
      · simp [h2]
    · simp [h1]
  · simp

/-- The charge part picks up `sgn ω` under the flip of `q_0`. -/
lemma chVec_flipQ (q : Cfg N) (ω : ZMod 2) (r : Cfg N) :
    chVec q ω (flipQ r) = sgn ω * chVec q ω r := by
  simp only [chVec]
  have h1 : (∀ k, k ≠ 0 → flipQ r k = q k) ↔ (∀ k, k ≠ 0 → r k = q k) := by
    constructor <;> intro h k hk
    · rw [← flipQ_ne r hk]; exact h k hk
    · rw [flipQ_ne r hk]; exact h k hk
  by_cases hr : ∀ k, k ≠ 0 → r k = q k
  · rw [ite_eq_left (h1.2 hr), ite_eq_left hr, flipQ_apply, ite_eq_left rfl, mul_add, mul_one,
      sgn_add]
    ring
  · rw [ite_eq_right (fun h => hr (h1.1 h)), ite_eq_right hr, mul_zero]

omit [NeZero N] in
/-- `dΓ` is additive. -/
lemma dGamma_add (A B : Matrix (Fin N) (Fin N) ℝ) : dGamma (A + B) = dGamma A + dGamma B := by
  simp only [dGamma, Matrix.add_apply, add_smul, Finset.sum_add_distrib]

omit [NeZero N] in
/-- `dΓ` is homogeneous. -/
lemma dGamma_smul (c : ℝ) (A : Matrix (Fin N) (Fin N) ℝ) : dGamma (c • A) = c • dGamma A := by
  simp only [dGamma, Matrix.smul_apply, smul_eq_mul, mul_smul, Finset.smul_sum]

/-- For `N ≥ 2`, second quantization of the ring bond adjacency is the ring bond hopping. -/
lemma dGamma_wrapAdj (hN : 2 ≤ N) : dGamma (wrapAdj (N := N)) = hopWrap := by
  have hne : (lastSite : Fin N) ≠ 0 := lastSite_ne_zero hN
  have e : ∀ x y : Fin N, wrapAdj x y • (cre x * ann y) =
      (if x = lastSite then (if y = 0 then cre x * ann y else 0) else 0) +
        (if x = 0 then (if y = lastSite then cre x * ann y else 0) else 0) := by
    intro x y
    simp only [wrapAdj]
    by_cases h1 : x = lastSite <;> by_cases h2 : y = 0 <;> by_cases h3 : x = 0 <;>
      by_cases h4 : y = lastSite <;> simp_all
  simp only [dGamma, e, Finset.sum_add_distrib, Finset.sum_ite_irrel, Finset.sum_const_zero,
    Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  rfl

/-- Second quantization of the twisted one-body Hamiltonian. -/
lemma dGamma_hqw (hN : 2 ≤ N) (J h : ℝ) (q : Cfg N) (ω : ZMod 2) :
    dGamma (hqw J h q ω) = (-J) • ∑ j, hop j +
      (2 * h) • diagonal (fun n => ∑ x, sgn (q x) * occ (n x)) + (-J * sgn ω) • hopWrap := by
  rw [hqw, dGamma_add, dGamma_smul, dGamma_hq, dGamma_wrapAdj hN]

/-- Entries of `HringDual`. -/
lemma HringDual_apply (J h : ℝ) (m r n r' : Cfg N) :
    HringDual J h (m, r) (n, r') =
      -J * ((∑ j, hop j m n) * (if r = r' then 1 else 0) + hopWrap m n * flip0 r r') -
        h * (if m = n ∧ r = r' then (sgn (m 0 + ∑ k, m k + ∑ k ∈ univ.erase 0, r k) +
          ∑ j ∈ univ.erase 0, sgn (m j + r j)) else 0) := by
  simp only [HringDual, Matrix.sub_apply, Matrix.smul_apply, Matrix.add_apply, Matrix.sum_apply,
    kronecker_apply, Matrix.one_apply, diagonal_apply, smul_eq_mul, Prod.mk.injEq, Finset.sum_mul]
  by_cases hmn : m = n
  · subst hmn
    by_cases hrr : r = r'
    · subst hrr; simp
    · simp [hrr]
  · simp [hmn]

/-- `∑_{r'} ⟨r|τᶻ_0|r'⟩ f(r') = f(flipQ r)`. -/
lemma sum_flip0_mul (r : Cfg N) (f : Cfg N → ℝ) : ∑ r', flip0 r r' * f r' = f (flipQ r) := by
  simp only [flip0_apply, ite_mul, one_mul, zero_mul]
  rw [Finset.sum_eq_single (flipQ r)]
  · rw [ite_eq_left (flipQ_flipQ r).symm]
  · intro r' _ hr'
    rw [ite_eq_right]
    intro h
    apply hr'
    rw [h, flipQ_flipQ]
  · simp

/-- The electric energy on the support of `ringVec`. -/
lemma elec_sum_dual {n q r : Cfg N} (hn : ∑ k, n k = ∑ k, q k)
    (hr : ∀ k, k ≠ 0 → r k = q k) :
    sgn (n 0 + ∑ k, n k + ∑ k ∈ univ.erase 0, r k) + ∑ j ∈ univ.erase 0, sgn (n j + r j) =
      ∑ j, sgn (q j) - 2 * ∑ x, sgn (q x) * occ (n x) := by
  have h0 : n 0 + ∑ k, n k + ∑ k ∈ univ.erase 0, r k = n 0 + q 0 := by
    rw [add_assoc, gauss0_sign hn hr]
  have h1 : ∑ j ∈ univ.erase 0, sgn (n j + r j) = ∑ j ∈ univ.erase 0, sgn (n j + q j) :=
    Finset.sum_congr rfl (fun j hj => by rw [hr j (Finset.ne_of_mem_erase hj)])
  have hsplit : ∑ j, sgn (n j + q j) = sgn (n 0 + q 0) + ∑ j ∈ univ.erase 0, sgn (n j + q j) :=
    (Finset.add_sum_erase _ _ (mem_univ 0)).symm
  rw [h0, h1, ← hsplit, Finset.mul_sum, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  rw [sgn_add, sgn_eq_one_sub_two_occ (n j)]
  ring

/-- The key computation of `Hring_ringVec` in the dual frame. -/
lemma HringDual_chVec (hN : 2 ≤ N) (J h : ℝ) (q : Cfg N) (ω : ZMod 2) (n : Cfg N)
    (hn : ∑ k, n k = ∑ k, q k) (m r : Cfg N) :
    ∑ r', HringDual J h (m, r) (n, r') * chVec q ω r' =
      (if ∑ k, m k = ∑ k, q k then
        (dGamma (hqw J h q ω) - (h * ∑ j, sgn (q j)) • 1) m n else 0) * chVec q ω r := by
  set G := sgn (m 0 + ∑ k, m k + ∑ k ∈ univ.erase 0, r k) + ∑ j ∈ univ.erase 0, sgn (m j + r j)
    with hG
  have e : ∀ r', HringDual J h (m, r) (n, r') * chVec q ω r' =
      (-J * ∑ j, hop j m n) * ((if r = r' then 1 else 0) * chVec q ω r') +
      (-J * hopWrap m n) * (flip0 r r' * chVec q ω r') -
      h * ((if m = n ∧ r = r' then G else 0) * chVec q ω r') := by
    intro r'; rw [HringDual_apply]; ring
  have s1 : ∑ r', (if r = r' then (1 : ℝ) else 0) * chVec q ω r' = chVec q ω r := by
    simp [ite_mul]
  have s3 : ∑ r', (if m = n ∧ r = r' then G else 0) * chVec q ω r' =
      if m = n then G * chVec q ω r else 0 := by
    by_cases hmn : m = n <;> simp [hmn, ite_mul]
  rw [Finset.sum_congr rfl (fun r' _ => e r'), Finset.sum_sub_distrib, Finset.sum_add_distrib,
    ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum, s1, s3, sum_flip0_mul, chVec_flipQ,
    dGamma_hqw hN]
  simp only [Matrix.sub_apply, Matrix.add_apply, Matrix.smul_apply, Matrix.sum_apply,
    diagonal_apply, Matrix.one_apply, smul_eq_mul]
  by_cases hξ : chVec q ω r = 0
  · simp [hξ]
  · have hr : ∀ k, k ≠ 0 → r k = q k := by
      by_contra hc
      apply hξ
      simp only [chVec, ite_eq_right hc]
    by_cases hmn : m = n
    · subst hmn
      rw [ite_eq_left hn]
      simp only [ite_true]
      rw [hG, elec_sum_dual hn hr]
      ring
    · by_cases hs : ∑ k, m k = ∑ k, q k
      · rw [ite_eq_left hs]
        simp only [hmn, ite_false]
        ring
      · rw [ite_eq_right hs]
        have hS : ∑ j, hop j m n = 0 := Finset.sum_eq_zero (fun j _ => by
          by_contra hc
          exact hs ((sum_of_hop hc).trans hn))
        have hW : hopWrap m n = 0 := by
          by_contra hc
          exact hs ((sum_of_hopWrap hc).trans hn)
        simp [hS, hW, hmn]

end RingAux

/-- The vectors `ringVec q ω n`, `∑ n = ∑ q`, are orthonormal. -/
theorem ringVec_orthonormal (q q' : Cfg N) (ω ω' : ZMod 2) (n n' : Cfg N)
    (hn : ∑ k, n k = ∑ k, q k) (hn' : ∑ k, n' k = ∑ k, q' k) :
    star (ringVec q' ω' n') ⬝ᵥ ringVec q ω n = if q' = q ∧ ω' = ω ∧ n' = n then 1 else 0 := by
  have hre : star (ringVec q' ω' n') ⬝ᵥ ringVec q ω n =
      ∑ a : Cfg N × Cfg N, star (ringVec q' ω' n' (dualEquiv a)) * ringVec q ω n (dualEquiv a) := by
    simp only [dotProduct, Pi.star_apply]
    exact (Equiv.sum_comp dualEquiv (fun p => star (ringVec q' ω' n' p) * ringVec q ω n p)).symm
  have e : ∀ a : Cfg N × Cfg N, star (ringVec q' ω' n' (dualEquiv a)) * ringVec q ω n (dualEquiv a) =
      (((if a.1 = n' then RingAux.chVec q' ω' a.2 else 0) *
        (if a.1 = n then RingAux.chVec q ω a.2 else 0) : ℝ) : ℂ) := by
    intro a
    rw [RingAux.ringVec_dualEquiv, RingAux.ringVec_dualEquiv, Complex.star_def,
      Complex.conj_ofReal, Complex.ofReal_mul]
  rw [hre, Finset.sum_congr rfl (fun a _ => e a), ← Complex.ofReal_sum, RingAux.sum_ringVecR]
  have hiff : (n' = n ∧ (∀ k, k ≠ 0 → q' k = q k) ∧ ω' = ω) ↔ (q' = q ∧ ω' = ω ∧ n' = n) := by
    constructor
    · rintro ⟨rfl, hq, hω⟩
      exact ⟨RingAux.eq_of_eqOff_of_sum hq (hn'.symm.trans hn), hω, rfl⟩
    · rintro ⟨rfl, hω, rfl⟩
      exact ⟨rfl, fun k _ => rfl, hω⟩
  by_cases hc : q' = q ∧ ω' = ω ∧ n' = n
  · rw [ite_eq_left (hiff.2 hc), ite_eq_left hc, Complex.ofReal_one]
  · rw [ite_eq_right (fun h => hc (hiff.1 h)), ite_eq_right hc, Complex.ofReal_zero]

/-- There are `4^N = dim` labels `(q, ω, n)` with `∑ n = ∑ q`: with orthonormality, the vectors
`ringVec` form an orthonormal basis, so every joint eigenspace `{G^ring = q, W = w}` is non-zero and
spanned by the `ringVec q ω n`. -/
theorem card_ringLabels :
    (univ.filter (fun a : Cfg N × ZMod 2 × Cfg N => ∑ k, a.2.2 k = ∑ k, a.1 k)).card =
      Fintype.card (Cfg N × Cfg N) := by
  have hin : ∀ q : Cfg N, ∑ n : Cfg N, (if ∑ k, n k = ∑ k, q k then 1 else 0) = 2 ^ (N - 1) := by
    intro q
    rw [← Finset.card_filter]
    exact card_parity_sector (NeZero.pos N) _
  rw [Finset.card_filter]
  simp only [Fintype.sum_prod_type, hin, Finset.sum_const, Finset.card_univ, smul_eq_mul,
    ZMod.card, Fintype.card_prod, Fintype.card_fun, Fintype.card_fin]
  have h2 : 2 ^ N = 2 * 2 ^ (N - 1) := by
    rw [← pow_succ']
    congr 1
    have := NeZero.pos N
    omega
  rw [h2]

/-- `G_j^ring ringVec = sgn (q j) ringVec`. -/
theorem gaussRing_ringVec (q : Cfg N) (ω : ZMod 2) (n : Cfg N) (hn : ∑ k, n k = ∑ k, q k)
    (j : Fin N) :
    cplx (gaussRing j) *ᵥ ringVec q ω n = ((sgn (q j) : ℝ) : ℂ) • ringVec q ω n := by
  funext p
  obtain ⟨⟨m, r⟩, rfl⟩ := dualEquiv.surjective p
  have hG : cplx (gaussRing j) = diagonal (fun p => ((gaussRing j p p : ℝ) : ℂ)) := by
    unfold gaussRing
    rw [sec_cplx_diagonal]
    simp [diagonal_apply_eq]
  rw [hG, mulVec_diagonal, Pi.smul_apply, smul_eq_mul, RingAux.ringVec_dualEquiv,
    RingAux.gaussRing_dualEquiv_apply]
  dsimp only
  by_cases h1 : m = n
  · subst h1
    simp only [ite_true, RingAux.chVec]
    by_cases h2 : ∀ k, k ≠ 0 → r k = q k
    · rw [ite_eq_left h2]
      have hs : (if j = 0 then sgn (∑ k, m k + ∑ k ∈ univ.erase 0, r k) else sgn (r j)) =
          sgn (q j) := by
        by_cases hj : j = 0
        · subst hj
          rw [ite_eq_left rfl, RingAux.gauss0_sign hn h2]
        · rw [ite_eq_right hj, h2 j hj]
      rw [hs]
    · simp [h2]
  · simp [h1]

/-- `W ringVec = sgn ω ringVec`. -/
theorem wilson_ringVec (q : Cfg N) (ω : ZMod 2) (n : Cfg N) :
    cplx wilson *ᵥ ringVec q ω n = ((sgn ω : ℝ) : ℂ) • ringVec q ω n := by
  funext p
  obtain ⟨⟨m, r⟩, rfl⟩ := dualEquiv.surjective p
  rw [RingAux.wilson_mulVec_dualEquiv, Pi.smul_apply, smul_eq_mul, RingAux.ringVec_dualEquiv,
    RingAux.ringVec_dualEquiv, ← Complex.ofReal_mul]
  congr 1
  dsimp only
  rw [RingAux.chVec_flipQ]
  split_ifs <;> ring

/-- `ringVec q ω n` is a joint eigenvector of the occupations, with eigenvalues `n`. -/
theorem occ_ringVec (q : Cfg N) (ω : ZMod 2) (n : Cfg N) (i : Fin N) :
    cplx (diagonal (fun p : Cfg N × Cfg N => occ (p.1 i))) *ᵥ ringVec q ω n =
      ((occ (n i) : ℝ) : ℂ) • ringVec q ω n := by
  funext p
  rw [sec_cplx_diagonal, mulVec_diagonal, Pi.smul_apply, smul_eq_mul]
  simp only [ringVec]
  split_ifs with hp
  · rw [hp.1]
  · simp

/-- **Twisted sector Hamiltonians** (paper: Proposition `prop:variants`(ii), Eq. `eq:ringtwist`): the
span of the `ringVec q ω n`, `∑ n = ∑ q`, is invariant under `H^ring`, which acts there as
`dΓ(h_{q,w}) - h ∑_j q_j` in the basis `|n⟩`. -/
theorem Hring_ringVec (hN : 2 ≤ N) (J h : ℝ) (q : Cfg N) (ω : ZMod 2) (n : Cfg N)
    (hn : ∑ k, n k = ∑ k, q k) :
    cplx (Hring J h) *ᵥ ringVec q ω n =
      ∑ n' ∈ univ.filter (fun n' : Cfg N => ∑ k, n' k = ∑ k, q k),
        (((dGamma (hqw J h q ω) - (h * ∑ j, sgn (q j)) • 1) n' n : ℝ) : ℂ) • ringVec q ω n' := by
  funext p
  obtain ⟨⟨m, r⟩, rfl⟩ := dualEquiv.surjective p
  rw [RingAux.cplx_mulVec_dualEquiv, ring_duality, Finset.sum_apply]
  have hL : ∑ a', ((HringDual J h (m, r) a' : ℝ) : ℂ) * ringVec q ω n (dualEquiv a') =
      ((∑ r', HringDual J h (m, r) (n, r') * RingAux.chVec q ω r' : ℝ) : ℂ) := by
    rw [Fintype.sum_prod_type, Finset.sum_eq_single n]
    · rw [Complex.ofReal_sum]
      refine Finset.sum_congr rfl (fun r' _ => ?_)
      rw [RingAux.ringVec_dualEquiv]
      simp
    · intro m' _ hm'
      apply Finset.sum_eq_zero
      intro r' _
      rw [RingAux.ringVec_dualEquiv]
      simp [hm']
    · simp
  have hR : ∑ n' ∈ univ.filter (fun n' : Cfg N => ∑ k, n' k = ∑ k, q k),
      ((((dGamma (hqw J h q ω) - (h * ∑ j, sgn (q j)) • 1) n' n : ℝ) : ℂ) • ringVec q ω n')
        (dualEquiv (m, r)) =
      (((if ∑ k, m k = ∑ k, q k then (dGamma (hqw J h q ω) - (h * ∑ j, sgn (q j)) • 1) m n
        else 0) * RingAux.chVec q ω r : ℝ) : ℂ) := by
    simp only [Pi.smul_apply, smul_eq_mul, RingAux.ringVec_dualEquiv]
    rw [Finset.sum_filter, Finset.sum_eq_single m]
    · by_cases hm : ∑ k, m k = ∑ k, q k <;> simp [hm]
    · intro n' _ hn'
      simp [Ne.symm hn']
    · simp
  rw [hL, hR, RingAux.HringDual_chVec hN J h q ω n hn m r]

/-! ### The `z`-polarized initial state -/

/-- The `z`-polarized gauge state has Wilson loop `W = +1`. -/
theorem wilson_initZ (ψ : Cfg N → ℂ) : cplx wilson *ᵥ initZ ψ = initZ ψ := by
  funext p
  obtain ⟨⟨m, r⟩, rfl⟩ := dualEquiv.surjective p
  rw [RingAux.wilson_mulVec_dualEquiv]
  simp [initZ]

/-- **Ring, sector law**: for a matter state of definite parity `π` and `z`-polarized links, the
dynamics of every matter observable diagonal in the occupations is the uniform average over the
`2^{N-1}` charge configurations with `∑_j q_j = π` of free-fermion dynamics with the periodic
one-body Hamiltonians `h_{q,+1}`. -/
theorem sector_average_ring (hN : 2 ≤ N) (J h : ℝ) (ψ : Cfg N → ℂ) (π : ZMod 2)
    (hψ : ∀ n, ψ n ≠ 0 → ∑ j, n j = π) (d : Cfg N → ℝ) (t : ℝ) :
    expect (cplx (diagonal d ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))) (evo (Hring J h) t *ᵥ initZ ψ) =
      (1 / 2 ^ (N - 1) : ℂ) * ∑ q ∈ univ.filter (fun q : Cfg N => ∑ j, q j = π),
        expect (cplx (diagonal d)) (evo (dGamma (hqw J h q 0)) t *ᵥ ψ) := by
  sorry

end DFL
