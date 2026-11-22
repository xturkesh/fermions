import DFL.Duality

/-!
# The standard ℤ₂ gauge theory: one charge in a random-walk potential (paper: Remark
`rem:randomwalk`)

With the single-link electric term `-h ∑_j σˣ_{b_j}` (instead of the two-link term of
Smith–Knolle–Kovrizhin–Moessner) the dual model is no longer free. In the one-particle sector,
however, it is an exact single-particle problem: in the charge sector `q`, a fermion at site `y`
feels the potential
`V_q(y) = 2h ∑_{j ≥ y} Q_j - h ∑_j Q_j`, with `Q_j = ∏_{k ≤ j} q_k`,
a random walk in `y` with i.i.d. ±1 increments when the charges are i.i.d. (the map `q ↦ Q` is a
bijection).
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The standard ℤ₂ LGT Hamiltonian with single-link electric term:
`H = -J ∑_j σᶻ_{b_j} (c_j† c_{j+1} + h.c.) - h ∑_j σˣ_{b_j}`. -/
def Hstd (J h : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  (-J) • ∑ j, hop j ⊗ₖ sz j - h • ∑ j, (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ sx j

/-- The configuration with a single fermion at site `y`. -/
def oneP (y : Fin N) : Cfg N := Pi.single y 1

/-- `Q_j = ∏_{k ≤ j} q_k`, as a sign. -/
def Qsgn (q : Cfg N) (j : Fin N) : ℝ := sgn (∑ k ∈ univ.filter (fun k => k ≤ j), q k)

/-- The random-walk potential of the charge sector `q`. -/
def rwPot (h : ℝ) (q : Cfg N) (y : Fin N) : ℝ :=
  2 * h * ∑ j ∈ univ.filter (fun j => y ≤ j), Qsgn q j - h * ∑ j, Qsgn q j

namespace RandomWalkAux

/-- Every element of `ZMod 2` is `0` or `1`. -/
lemma zmod2_cases (a : ZMod 2) : a = 0 ∨ a = 1 := by
  revert a; decide

/-- `sgn` turns addition in `ZMod 2` into multiplication. -/
lemma sgn_add (a b : ZMod 2) : sgn (a + b) = sgn a * sgn b := by
  have h11 : (1 : ZMod 2) + 1 = 0 := by decide
  rcases zmod2_cases a with rfl | rfl <;> rcases zmod2_cases b with rfl | rfl <;> simp [sgn, h11]

/-- For fixed occupations, `q ↦ linkOf n q` is injective (iff form). -/
lemma linkOf_inj (n q' q : Cfg N) : linkOf n q' = linkOf n q ↔ q' = q := by
  constructor
  · intro hq
    rw [← chargeOf_linkOf n q', hq, chargeOf_linkOf]
  · rintro rfl
    rfl

/-- The hopping term in dual variables conserves the charges. -/
lemma hop_mul_sz (j : Fin N) (n' n q' q : Cfg N) :
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

/-- The electric term `1 ⊗ σˣ_{b_j}` in dual variables: diagonal, with entries `sgn (x_{b_j})`. -/
lemma one_mul_sx (j : Fin N) (n' n q' q : Cfg N) :
    (1 : Matrix (Cfg N) (Cfg N) ℝ) n' n * sx j (linkOf n' q') (linkOf n q) =
      if q' = q then (if n' = n then sgn (linkOf n q j) else 0) else 0 := by
  rw [Matrix.one_apply, sx, diagonal_apply]
  by_cases hn : n' = n
  · subst hn
    by_cases hq : q' = q
    · subst hq
      simp
    · simp [linkOf_inj, hq]
  · simp [hn]

/-- Occupations of the one-particle configuration. -/
lemma oneP_apply (y k : Fin N) : oneP y k = if k = y then 1 else 0 := by
  simp [oneP, Pi.single_apply]

/-- `y ↦ oneP y` is injective (iff form). -/
lemma oneP_inj (y' y : Fin N) : oneP (N := N) y' = oneP y ↔ y' = y :=
  Pi.single_left_inj one_ne_zero

/-- Moving the particle with a transposition of sites. -/
lemma oneP_comp_swap (y a b : Fin N) : oneP y ∘ Equiv.swap a b = oneP (Equiv.swap a b y) := by
  rw [oneP, Pi.single_comp_equiv, Equiv.symm_swap]
  rfl

/-- One-particle matrix elements of the hopping term across the bond `(j, j+1)`. -/
lemma hop_oneP (j y' y : Fin N) :
    hop j (oneP y') (oneP y) =
      if (y.val = j.val ∧ y'.val = j.val + 1) ∨ (y.val = j.val + 1 ∧ y'.val = j.val) then 1
      else 0 := by
  rw [hop_apply]
  by_cases hj : j.val + 1 < N
  · rw [dite_eq_left hj]
    apply if_congr _ rfl rfl
    rw [oneP_comp_swap, oneP_inj, oneP_apply, oneP_apply]
    have hne : (⟨j.val + 1, hj⟩ : Fin N) ≠ j := fun h => by simp [Fin.ext_iff] at h
    by_cases h1 : j = y
    · subst h1
      rw [Equiv.swap_apply_left, ite_eq_left rfl, ite_eq_right hne]
      simp [Fin.ext_iff]
    · by_cases h2 : (⟨j.val + 1, hj⟩ : Fin N) = y
      · subst h2
        rw [Equiv.swap_apply_right, ite_eq_right h1, ite_eq_left rfl]
        simp [Fin.ext_iff]
      · rw [ite_eq_right h1, ite_eq_right h2]
        simp only [ne_eq, not_true_eq_false, false_and, false_iff]
        have h1' : j.val ≠ y.val := fun h => h1 (Fin.ext h)
        have h2' : j.val + 1 ≠ y.val := fun h => h2 (Fin.ext h)
        omega
  · rw [dite_eq_right hj, ite_eq_right]
    have := y.isLt
    have := y'.isLt
    omega

/-- Summing the hopping terms in the one-particle sector gives the path adjacency matrix. -/
lemma sum_hop_oneP (y' y : Fin N) : ∑ j, hop j (oneP y') (oneP y) = pathAdj N y' y := by
  simp only [hop_oneP, pathAdj]
  by_cases h1 : y.val + 1 = y'.val
  · rw [ite_eq_left (Or.inr h1), Finset.sum_eq_single y]
    · rw [ite_eq_left (Or.inl ⟨rfl, h1.symm⟩)]
    · intro j _ hj
      have : j.val ≠ y.val := fun h => hj (Fin.ext h)
      rw [ite_eq_right]
      omega
    · simp
  · by_cases h2 : y'.val + 1 = y.val
    · rw [ite_eq_left (Or.inl h2), Finset.sum_eq_single y']
      · rw [ite_eq_left (Or.inr ⟨h2.symm, rfl⟩)]
      · intro j _ hj
        have : j.val ≠ y'.val := fun h => hj (Fin.ext h)
        rw [ite_eq_right]
        omega
      · simp
    · rw [ite_eq_right (by omega)]
      apply Finset.sum_eq_zero
      intro j _
      rw [ite_eq_right]
      omega

/-- The link labels of the one-particle configuration: `sgn x_{b_j} = Q_j (-1)^{[y ≤ j]}`. -/
lemma sgn_linkOf_oneP (q : Cfg N) (y j : Fin N) :
    sgn (linkOf (oneP y) q j) = Qsgn q j * (if y ≤ j then -1 else 1) := by
  have e : linkOf (oneP y) q j =
      (∑ k ∈ univ.filter (fun k => k ≤ j), q k) + (if y ≤ j then 1 else 0) := by
    simp only [linkOf, Finset.sum_add_distrib, oneP, Finset.sum_pi_single', Finset.mem_filter,
      Finset.mem_univ, true_and]
  rw [e, sgn_add, Qsgn]
  split_ifs <;> simp [sgn]

/-- The electric energy of the one-particle configuration is the random-walk potential. -/
lemma elec_oneP (h : ℝ) (q : Cfg N) (y : Fin N) :
    -(h * ∑ j, sgn (linkOf (oneP y) q j)) = rwPot h q y := by
  have e : ∀ j, sgn (linkOf (oneP y) q j) = Qsgn q j - 2 * (if y ≤ j then Qsgn q j else 0) := by
    intro j
    rw [sgn_linkOf_oneP]
    split_ifs <;> ring
  rw [rwPot, Finset.sum_filter, Finset.sum_congr rfl (fun j _ => e j), Finset.sum_sub_distrib,
    ← Finset.mul_sum]
  ring

end RandomWalkAux

open RandomWalkAux

/-- **One charge in a random-walk potential**: in the dual frame, the one-particle block of the
standard ℤ₂ LGT in the sector `q` is `-J Δ + diag(V_q)`. -/
theorem Hstd_dual_oneP (J h : ℝ) (y y' : Fin N) (q q' : Cfg N) :
    Hstd J h (dualEquiv (oneP y', q')) (dualEquiv (oneP y, q)) =
      if q' = q then (-J) * pathAdj N y' y + (if y' = y then rwPot h q y else 0) else 0 := by
  rw [dualEquiv_apply, dualEquiv_apply]
  simp only [Hstd, Matrix.sub_apply, Matrix.smul_apply, Matrix.sum_apply, kronecker_apply,
    smul_eq_mul, hop_mul_sz, one_mul_sx]
  by_cases hq : q' = q
  · simp only [hq, ite_true]
    rw [sum_hop_oneP]
    by_cases hy : y' = y
    · simp only [hy, ite_true]
      rw [← elec_oneP]
      ring
    · simp only [oneP_inj, hy, ite_false, Finset.sum_const_zero, mul_zero, sub_zero, add_zero]
  · simp only [hq, ite_false, Finset.sum_const_zero, mul_zero, sub_zero]

/-- The partial-sum map `q ↦ (∑_{k ≤ j} q_k)_j` is a bijection, so `Q` is i.i.d. uniform when `q` is. -/
theorem partialSum_bijective :
    Function.Bijective (fun q : Cfg N => fun j : Fin N => ∑ k ∈ univ.filter (fun k => k ≤ j), q k) := by
  have e : (fun q : Cfg N => fun j : Fin N => ∑ k ∈ univ.filter (fun k => k ≤ j), q k) =
      linkOf 0 := by
    funext q j
    simp [linkOf]
  rw [e]
  exact Function.bijective_iff_has_inverse.mpr ⟨chargeOf 0, chargeOf_linkOf 0, linkOf_chargeOf 0⟩

end DFL
