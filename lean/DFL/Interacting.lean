import DFL.Ensemble
import DFL.Entropies
import DFL.GaussLaw

/-!
# Interacting matter (paper: Remark `rem:interacting`)

Following Smith, Knolle, Moessner and Kovrizhin, add the density interaction
`Δ ∑_{j=1}^{N-1} (2n_j - 1)(2n_{j+1} - 1)` to the Hamiltonian: `H_Δ = H + Δ ∑_j (2n_j - 1)(2n_{j+1} - 1)`
(`Hint`). The interaction is a function of the occupations, so it commutes with every Gauss operator
(`Hint_comm_gauss`) and with the duality, and:

* `duality_int`: `U† H_Δ U = ⊕_q (dΓ(h_q) - h ∑_j q_j + Δ ∑_j (2n_j - 1)(2n_{j+1} - 1))`;
* `HsecInt_eq_xxzChain`: with `S^z_j = n_j - ½` and `S^±_j` the spin raising and lowering operators in the
  occupation basis (`c_j† c_{j+1} = S⁺_j S⁻_{j+1}`: no Jordan–Wigner sign for neighbours), the sector
  Hamiltonian is the XXZ chain in the binary random field `2h q_j`,
  `-J ∑ (S⁺_j S⁻_{j+1} + S⁻_j S⁺_{j+1}) + 4Δ ∑ S^z_j S^z_{j+1} + 2h ∑ q_j S^z_j`; the constants cancel exactly;
* `sector_average_int`, `sector_average_matter_int`: exact sector averaging for `|↑⋯↑⟩_z ⊗ ψ`, every
  gauge-invariant expectation value is the average over i.i.d. charges of the XXZ chains;
  `correlator_sector_average_int`: the same for products `O_1(t_1)⋯O_m(t_m)` at any number of times
  (with `Δ = 0`, `correlator_sector_average` is the multi-time statement of Corollary `cor:sector`(b));
  `correlator_eq_ensemble_int`: the superposition equals the ensemble of definite-charge states;
  `lgt_duality_defect_int`: the duality changes the entanglement of a left segment by at most `ln 2`;
* `nnInt_one_particle`: with a single particle the interaction is the constant `N - 5` plus the boundary
  potential `2(n_1 + n_N)`, so one particle still sees a binary Anderson chain; interactions matter from
  two particles on.

Localization for `Δ ≠ 0` would be many-body localization of the XXZ chain with binary disorder, an open
problem; it is not claimed here.
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The density interaction `∑_{j=1}^{N-1} (2n_j - 1)(2n_{j+1} - 1)`, diagonal in the occupations
(Lean sites `j, j+1` with `j + 1 < N`). -/
def nnInt : Matrix (Cfg N) (Cfg N) ℝ :=
  diagonal (fun n => ∑ j : Fin N, if hj : j.val + 1 < N then
    (2 * occ (n j) - 1) * (2 * occ (n ⟨j.val + 1, hj⟩) - 1) else 0)

/-- **The interacting gauge theory** `H_Δ = H + Δ ∑_j (2n_j - 1)(2n_{j+1} - 1)` (paper: Remark
`rem:interacting`). -/
def Hint (J h Δ : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  H J h + (Δ • nnInt) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)

/-- The interacting sector Hamiltonian `dΓ(h_q) - h ∑_j q_j + Δ ∑_j (2n_j - 1)(2n_{j+1} - 1)`. -/
def HsecInt (J h Δ : ℝ) (q : Cfg N) : Matrix (Cfg N) (Cfg N) ℝ :=
  Hsec J h q + Δ • nnInt

namespace InteractingAux

/-- The interaction term on the LGT Hilbert space is diagonal in `(n, x)`. -/
lemma nnInt_kron_eq (Δ : ℝ) :
    (Δ • nnInt (N := N)) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) =
      diagonal (fun p : Cfg N × Cfg N => Δ * nnInt p.1 p.1) := by
  rw [nnInt, ← diagonal_smul, sec_diagonal_kronecker_one]
  congr 1
  funext p
  simp [diagonal_apply_eq]

/-- In the dual frame the interaction term is block diagonal, with the same block in every
sector. -/
lemma nnInt_kron_dual (Δ : ℝ) :
    ((Δ • nnInt (N := N)) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)).submatrix dualEquiv dualEquiv =
      blockDiagonal (fun _ : Cfg N => Δ • nnInt (N := N)) := by
  rw [nnInt, ← diagonal_smul, sec_diagonal_kronecker_one, submatrix_diagonal_equiv,
    blockDiagonal_diagonal]
  rfl

end InteractingAux

open InteractingAux

/-- Without interaction, `H_0 = H`. -/
theorem Hint_zero (J h : ℝ) : Hint (N := N) J h 0 = H J h := by
  rw [Hint, zero_smul, zero_kronecker, add_zero]

/-- `H_Δ` is gauge invariant. -/
theorem Hint_comm_gauss (J h Δ : ℝ) (j : Fin N) : gauss j * Hint J h Δ = Hint J h Δ * gauss j := by
  have h2 : gauss j * (Δ • nnInt) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) =
      (Δ • nnInt) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) * gauss j := by
    rw [nnInt_kron_eq, gauss, diagonal_mul_diagonal, diagonal_mul_diagonal]
    congr 1
    funext p
    ring
  rw [Hint, mul_add, add_mul, gauss_comm_H, h2]

/-- **Exact duality for interacting matter** (paper: Remark `rem:interacting`):
`U† H_Δ U = ⊕_q (dΓ(h_q) - h ∑_j q_j + Δ ∑_j (2n_j - 1)(2n_{j+1} - 1))`. -/
theorem duality_int (J h Δ : ℝ) :
    (Hint J h Δ).submatrix (dualEquiv (N := N)) dualEquiv =
      blockDiagonal (fun q => HsecInt J h Δ q) := by
  have h3 : (Hint J h Δ).submatrix (dualEquiv (N := N)) dualEquiv =
      (H J h).submatrix dualEquiv dualEquiv +
        ((Δ • nnInt (N := N)) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)).submatrix dualEquiv dualEquiv :=
    rfl
  rw [h3, duality, nnInt_kron_dual, Hdual, ← blockDiagonal_add]
  rfl

/-! ### The XXZ chain -/

/-- `S⁺_j` in the occupation basis: raises `n_j` from `0` to `1` (no Jordan–Wigner string). -/
def spinPlusAt (j : Fin N) : Matrix (Cfg N) (Cfg N) ℝ :=
  fun n' n => if n j = 0 ∧ n' = Function.update n j 1 then 1 else 0

/-- `S⁻_j = (S⁺_j)ᵀ`. -/
def spinMinusAt (j : Fin N) : Matrix (Cfg N) (Cfg N) ℝ := (spinPlusAt j)ᵀ

/-- `S^z_j = n_j - ½`. -/
def spinZAt (j : Fin N) : Matrix (Cfg N) (Cfg N) ℝ := diagonal (fun n => occ (n j) - 1 / 2)

/-- The XXZ chain in the binary field `2h q_j` (SKKM 2017b, Eq. (5)):
`-J ∑_j (S⁺_j S⁻_{j+1} + S⁻_j S⁺_{j+1}) + 4Δ ∑_j S^z_j S^z_{j+1} + 2h ∑_j q_j S^z_j`. -/
def xxzChain (J h Δ : ℝ) (q : Cfg N) : Matrix (Cfg N) (Cfg N) ℝ :=
  (-J) • ∑ j : Fin N, (if hj : j.val + 1 < N then
      spinPlusAt j * spinMinusAt ⟨j.val + 1, hj⟩ + spinMinusAt j * spinPlusAt ⟨j.val + 1, hj⟩
    else 0) +
  (4 * Δ) • ∑ j : Fin N, (if hj : j.val + 1 < N then spinZAt j * spinZAt ⟨j.val + 1, hj⟩ else 0) +
  (2 * h) • ∑ j : Fin N, sgn (q j) • spinZAt j

namespace InteractingAux

/-- Every element of `ZMod 2` is `0` or `1`. -/
lemma zmod2_cases (a : ZMod 2) : a = 0 ∨ a = 1 := by
  revert a; decide

/-- Emptying the occupied site `b` and filling the empty site `a` exchanges `a` and `b`. -/
lemma update_update_eq_comp_swap {a b : Fin N} (hab : a ≠ b) (n : Cfg N) (ha : n a = 0)
    (hb : n b = 1) : Function.update (Function.update n b 0) a 1 = n ∘ Equiv.swap a b := by
  funext k
  by_cases hka : k = a
  · subst hka; simp [hb]
  · by_cases hkb : k = b
    · subst hkb; simp [Function.update_of_ne hka, ha]
    · simp [Function.update_of_ne hka, Function.update_of_ne hkb,
        Equiv.swap_apply_of_ne_of_ne hka hkb]

/-- Support condition of `S⁻_a S⁺_b` (`a ≠ b`, `n_b = 0`): the particle at `a` moves to `b`. -/
lemma update_eq_update_iff {a b : Fin N} (hab : a ≠ b) (n n' : Cfg N) (hb : n b = 0) :
    (n' a = 0 ∧ Function.update n b 1 = Function.update n' a 1) ↔
      (n a = 1 ∧ n' = n ∘ Equiv.swap a b) := by
  constructor
  · rintro ⟨ha', h⟩
    have hna : n a = 1 := by
      have := congrFun h a
      rwa [Function.update_of_ne hab, Function.update_self] at this
    refine ⟨hna, ?_⟩
    funext k
    by_cases hka : k = a
    · subst hka; simp [ha', hb]
    · by_cases hkb : k = b
      · subst hkb
        have := congrFun h k
        rw [Function.update_self, Function.update_of_ne hka] at this
        simp [← this, hna]
      · have := congrFun h k
        rw [Function.update_of_ne hkb, Function.update_of_ne hka] at this
        simp [Equiv.swap_apply_of_ne_of_ne hka hkb, this]
  · rintro ⟨hna, rfl⟩
    refine ⟨by simp [hb], ?_⟩
    funext k
    by_cases hka : k = a
    · subst hka; simp [Function.update_of_ne hab, hna]
    · by_cases hkb : k = b
      · subst hkb; simp [Function.update_of_ne hka, hna]
      · simp [Function.update_of_ne hka, Function.update_of_ne hkb,
          Equiv.swap_apply_of_ne_of_ne hka hkb]

/-- Entries of `S⁺_a S⁻_b` (`a ≠ b`): a particle hops from `b` to `a`. -/
lemma spinPlus_mul_spinMinus_apply {a b : Fin N} (hab : a ≠ b) (n' n : Cfg N) :
    (spinPlusAt a * spinMinusAt b) n' n =
      if n a = 0 ∧ n b = 1 ∧ n' = n ∘ Equiv.swap a b then 1 else 0 := by
  simp only [Matrix.mul_apply, spinMinusAt, Matrix.transpose_apply, spinPlusAt]
  rw [Finset.sum_eq_single (Function.update n b 0)]
  · rw [Function.update_of_ne hab, Function.update_self, Function.update_idem]
    by_cases hb : n b = 1
    · have hn : Function.update n b 1 = n := Function.update_eq_self_iff.2 hb.symm
      by_cases ha : n a = 0
      · rw [update_update_eq_comp_swap hab n ha hb, hn]
        simp [ha, hb]
      · simp [ha]
    · have hn : ¬ n = Function.update n b 1 := fun h => hb (by rw [h, Function.update_self])
      simp [hb, hn]
  · intro m _ hm
    rw [mul_eq_zero]
    right
    rw [ite_eq_right_iff]
    rintro ⟨hm0, rfl⟩
    exact absurd (by rw [Function.update_idem, Function.update_eq_self_iff.2 hm0.symm]) hm
  · simp

/-- Entries of `S⁻_a S⁺_b` (`a ≠ b`): a particle hops from `a` to `b`. -/
lemma spinMinus_mul_spinPlus_apply {a b : Fin N} (hab : a ≠ b) (n' n : Cfg N) :
    (spinMinusAt a * spinPlusAt b) n' n =
      if n a = 1 ∧ n b = 0 ∧ n' = n ∘ Equiv.swap a b then 1 else 0 := by
  simp only [Matrix.mul_apply, spinMinusAt, Matrix.transpose_apply, spinPlusAt]
  rw [Finset.sum_eq_single (Function.update n b 1)]
  · by_cases hb : n b = 0
    · simp only [update_eq_update_iff hab n n' hb]
      simp [hb]
    · simp [hb]
  · intro m _ hm
    rw [mul_eq_zero]
    right
    rw [ite_eq_right_iff]
    rintro ⟨_, rfl⟩
    exact absurd rfl hm
  · simp

/-- **Nearest-neighbour hopping in spin language**: `c_j† c_{j+1} + c_{j+1}† c_j =
S⁺_j S⁻_{j+1} + S⁻_j S⁺_{j+1}` (no Jordan–Wigner sign between neighbours). -/
lemma hop_eq_spin (j : Fin N) :
    hop j = if hj : j.val + 1 < N then
      spinPlusAt j * spinMinusAt ⟨j.val + 1, hj⟩ + spinMinusAt j * spinPlusAt ⟨j.val + 1, hj⟩
    else 0 := by
  ext n' n
  rw [hop_apply, Matrix.dite_apply]
  by_cases hj : j.val + 1 < N
  · rw [dite_eq_left hj, dite_eq_left hj]
    have hjj' : j ≠ ⟨j.val + 1, hj⟩ := fun h => by
      have := congrArg Fin.val h
      simp at this
    rw [Matrix.add_apply, spinPlus_mul_spinMinus_apply hjj', spinMinus_mul_spinPlus_apply hjj']
    rcases zmod2_cases (n j) with h0 | h0 <;>
      rcases zmod2_cases (n ⟨j.val + 1, hj⟩) with h1 | h1 <;>
      by_cases hs : n' = n ∘ Equiv.swap j ⟨j.val + 1, hj⟩ <;> simp [h0, h1, hs]
  · rw [dite_eq_right hj, dite_eq_right hj, Matrix.zero_apply]

end InteractingAux

/-- **The sector Hamiltonians are XXZ chains in a binary random field** (paper: Remark
`rem:interacting`); the constants cancel exactly. -/
theorem HsecInt_eq_xxzChain (J h Δ : ℝ) (q : Cfg N) : HsecInt J h Δ q = xxzChain J h Δ q := by
  have e1 : ∑ j, hop j = ∑ j : Fin N, (if hj : j.val + 1 < N then
      spinPlusAt j * spinMinusAt ⟨j.val + 1, hj⟩ + spinMinusAt j * spinPlusAt ⟨j.val + 1, hj⟩
      else 0) :=
    Finset.sum_congr rfl fun j _ => hop_eq_spin j
  have e2 : ∑ j : Fin N, (if hj : j.val + 1 < N then spinZAt j * spinZAt ⟨j.val + 1, hj⟩ else 0) =
      diagonal (fun n => ∑ j : Fin N, if hj : j.val + 1 < N then
        (occ (n j) - 1 / 2) * (occ (n ⟨j.val + 1, hj⟩) - 1 / 2) else 0) := by
    ext n' n
    simp only [Matrix.sum_apply, Matrix.dite_apply, spinZAt, diagonal_mul_diagonal, diagonal_apply,
      Matrix.zero_apply]
    by_cases hn : n' = n
    · subst hn; simp
    · simp [hn]
  have e3 : ∑ j, sgn (q j) • spinZAt j =
      diagonal (fun n => ∑ j, sgn (q j) * (occ (n j) - 1 / 2)) := by
    ext n' n
    simp only [Matrix.sum_apply, Matrix.smul_apply, spinZAt, diagonal_apply, smul_eq_mul]
    by_cases hn : n' = n
    · subst hn; simp
    · simp [hn]
  unfold HsecInt Hsec xxzChain
  rw [dGamma_hq, e1, e2, e3]
  ext n' n
  simp only [Matrix.add_apply, Matrix.sub_apply, Matrix.smul_apply, smul_eq_mul, diagonal_apply,
    Matrix.one_apply, nnInt]
  by_cases hn : n' = n
  · subst hn
    simp only [ite_true, mul_one]
    have hA : ∑ j, sgn (q j) * (occ (n' j) - 1 / 2) =
        ∑ x, sgn (q x) * occ (n' x) - 1 / 2 * ∑ j, sgn (q j) := by
      rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
      exact Finset.sum_congr rfl fun j _ => by ring
    have hB : ∑ j : Fin N, (if hj : j.val + 1 < N then
        (occ (n' j) - 1 / 2) * (occ (n' ⟨j.val + 1, hj⟩) - 1 / 2) else 0) =
        1 / 4 * ∑ j : Fin N, (if hj : j.val + 1 < N then
          (2 * occ (n' j) - 1) * (2 * occ (n' ⟨j.val + 1, hj⟩) - 1) else 0) := by
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl fun j _ => ?_
      split_ifs <;> ring
    rw [hA, hB]
    ring
  · simp [hn]

/-- **One particle**: on configurations with a single particle the interaction is the constant `N - 5`
plus the boundary potential `2(n_1 + n_N)` (Lean sites `0` and `N - 1`). -/
theorem nnInt_one_particle (hN : 2 ≤ N) (n : Cfg N) (h1 : count univ n = 1) :
    nnInt n n = ((N : ℝ) - 5) + 2 * (occ (n ⟨0, by omega⟩) + occ (n ⟨N - 1, by omega⟩)) := by
  obtain ⟨M, rfl⟩ : ∃ M, N = M + 1 := ⟨N - 1, by omega⟩
  have hocc : ∀ x, occ (n x) = 0 ∨ occ (n x) = 1 := fun x => by
    unfold occ; split_ifs <;> simp
  have hnn : ∀ x, 0 ≤ occ (n x) := fun x => by
    rcases hocc x with h | h <;> simp [h]
  have hT : ∑ x, occ (n x) = 1 := h1
  have hpair : ∀ x y : Fin (M + 1), x ≠ y → occ (n x) * occ (n y) = 0 := by
    intro x y hxy
    have hle : occ (n x) + occ (n y) ≤ 1 := by
      rw [← Finset.sum_pair (f := fun z => occ (n z)) hxy, ← hT]
      exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _)
        (fun i _ _ => hnn i)
    rcases hocc x with h | h
    · rw [h, zero_mul]
    · rcases hocc y with h' | h'
      · rw [h', mul_zero]
      · rw [h, h'] at hle
        norm_num at hle
  have hS1 : ∑ i : Fin M, occ (n (Fin.castSucc i)) = 1 - occ (n (Fin.last M)) := by
    rw [← hT, Fin.sum_univ_castSucc]; ring
  have hS2 : ∑ i : Fin M, occ (n i.succ) = 1 - occ (n 0) := by
    rw [← hT, Fin.sum_univ_succ]; ring
  have hP : ∀ i : Fin M, occ (n (Fin.castSucc i)) * occ (n i.succ) = 0 := fun i =>
    hpair _ _ Fin.castSucc_lt_succ.ne
  have e0 : (⟨0, by omega⟩ : Fin (M + 1)) = 0 := rfl
  have eL : (⟨M + 1 - 1, by omega⟩ : Fin (M + 1)) = Fin.last M := Fin.ext (by simp)
  rw [e0, eL, nnInt, diagonal_apply_eq, Fin.sum_univ_castSucc]
  have hlast : ¬ ((Fin.last M).val + 1 < M + 1) := by simp
  rw [dite_eq_right hlast, add_zero]
  have hterm : ∀ i : Fin M, (if hj : (Fin.castSucc i).val + 1 < M + 1 then
      (2 * occ (n (Fin.castSucc i)) - 1) * (2 * occ (n ⟨(Fin.castSucc i).val + 1, hj⟩) - 1)
      else 0) =
      4 * (occ (n (Fin.castSucc i)) * occ (n i.succ)) - 2 * occ (n (Fin.castSucc i)) -
        2 * occ (n i.succ) + 1 := by
    intro i
    have hi : (Fin.castSucc i).val + 1 < M + 1 := by simp
    rw [dite_eq_left hi]
    have : (⟨(Fin.castSucc i).val + 1, hi⟩ : Fin (M + 1)) = i.succ := Fin.ext (by simp)
    rw [this]
    ring
  rw [Finset.sum_congr rfl (fun i _ => hterm i), Finset.sum_add_distrib, Finset.sum_sub_distrib,
    Finset.sum_sub_distrib, ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum, hS1, hS2,
    Finset.sum_congr rfl (fun i _ => hP i)]
  simp only [Finset.sum_const_zero, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
    nsmul_eq_mul, mul_one]
  push_cast
  ring

/-! ### Exact sector averaging -/

/-- The interacting LGT state at time `t` evolved from `|↑…↑⟩_z ⊗ ψ`. -/
def lgtStateInt (J h Δ : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) : Cfg N × Cfg N → ℂ :=
  evo (Hint J h Δ) t *ᵥ initZ ψ

/-- The Heisenberg picture `O(t) = e^{itH_Δ} O e^{-itH_Δ}`. -/
def heisInt (J h Δ : ℝ) (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (t : ℝ) :
    Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  (evo (Hint J h Δ) t)ᴴ * O * evo (Hint J h Δ) t

namespace InteractingAux

/-- The interacting propagator is block diagonal in the dual frame, with the sector propagators as
blocks. -/
lemma evo_Hint_dual (J h Δ : ℝ) (t : ℝ) :
    (evo (Hint (N := N) J h Δ) t).submatrix dualEquiv dualEquiv =
      blockDiagonal (fun q => evo (HsecInt J h Δ q) t) := by
  rw [evo_submatrix, duality_int, evo_blockDiagonal]

/-- The initial state `|↑⋯↑⟩_z ⊗ ψ` in the dual frame is `2^{-N/2} ψ` in every sector. -/
lemma initZ_dual (ψ : Cfg N → ℂ) (q : Cfg N) :
    (fun n => initZ ψ (dualEquiv (n, q))) = (((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ • ψ := by
  ext n
  simp [initZ]

/-- **Dual-frame expectation formula** for the interacting dynamics: for every initial state `Ψ₀`
and every gauge-invariant `O`, the expectation is a sum over the charge sectors. -/
lemma expect_evo_dual_int (J h Δ : ℝ) (Ψ₀ : Cfg N × Cfg N → ℂ)
    (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (Hint J h Δ) t *ᵥ Ψ₀) =
      ∑ q : Cfg N, expect (sectorBlock O q)
        (evo (HsecInt J h Δ q) t *ᵥ fun n => Ψ₀ (dualEquiv (n, q))) := by
  have hstate : (evo (Hint J h Δ) t *ᵥ Ψ₀) ∘ dualEquiv =
      blockDiagonal (fun q : Cfg N => evo (HsecInt J h Δ q) t) *ᵥ (Ψ₀ ∘ dualEquiv) := by
    rw [← evo_Hint_dual, submatrix_mulVec_equiv, Function.comp_assoc, Equiv.self_comp_symm,
      Function.comp_id]
  rw [← sec_expect_submatrix O dualEquiv, gaugeInv_submatrix O hO, hstate,
    sec_expect_blockDiagonal]
  refine Finset.sum_congr rfl fun q _ => ?_
  congr 1
  ext n
  rw [sec_blockDiagonal_mulVec]
  rfl

/-- A Heisenberg-picture gauge-invariant observable is block diagonal in the dual frame, with the
sector Heisenberg-picture blocks. -/
lemma heisInt_dual (J h Δ : ℝ) (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O)
    (t : ℝ) :
    (heisInt J h Δ O t).submatrix dualEquiv dualEquiv =
      blockDiagonal (fun q => (evo (HsecInt J h Δ q) t)ᴴ * sectorBlock O q *
        evo (HsecInt J h Δ q) t) := by
  rw [heisInt, ← submatrix_mul_equiv (e₂ := dualEquiv), ← submatrix_mul_equiv (e₂ := dualEquiv),
    ← conjTranspose_submatrix, evo_Hint_dual, gaugeInv_submatrix O hO, blockDiagonal_conjTranspose,
    ← blockDiagonal_mul, ← blockDiagonal_mul]

/-- Products of Heisenberg-picture gauge-invariant observables are block diagonal in the dual
frame. -/
lemma list_prod_heisInt_dual (J h Δ : ℝ)
    (Os : List (Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ × ℝ)) (hOs : ∀ p ∈ Os, GaugeInv p.1) :
    ((Os.map (fun p => heisInt J h Δ p.1 p.2)).prod).submatrix dualEquiv dualEquiv =
      blockDiagonal (fun q => (Os.map (fun p =>
        (evo (HsecInt J h Δ q) p.2)ᴴ * sectorBlock p.1 q * evo (HsecInt J h Δ q) p.2)).prod) := by
  induction Os with
  | nil =>
    simp only [List.map_nil, List.prod_nil, submatrix_one_equiv]
    exact blockDiagonal_one.symm
  | cons p Os ih =>
    rw [List.map_cons, List.prod_cons, ← submatrix_mul_equiv (e₂ := dualEquiv),
      ih fun p' hp' => hOs p' (List.mem_cons_of_mem _ hp'),
      heisInt_dual J h Δ p.1 (hOs p List.mem_cons_self) p.2, ← blockDiagonal_mul]
    simp only [List.map_cons, List.prod_cons]

end InteractingAux

/-- **Exact sector averaging for interacting matter** (paper: Remark `rem:interacting`). -/
theorem sector_average_int (J h Δ : ℝ) (ψ : Cfg N → ℂ)
    (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) (t : ℝ) :
    expect O (lgtStateInt J h Δ ψ t) =
      (1 / 2 ^ N : ℂ) * ∑ q : Cfg N, expect (sectorBlock O q) (evo (HsecInt J h Δ q) t *ᵥ ψ) := by
  unfold lgtStateInt
  rw [expect_evo_dual_int J h Δ _ O hO t, Finset.mul_sum]
  refine Finset.sum_congr rfl fun q _ => ?_
  rw [initZ_dual, mulVec_smul, sec_expect_smul, sec_norm_sq_prefactor]

/-- For functions of the occupations: the average over i.i.d. charges of the XXZ chains in the binary
field `2h q_j`. -/
theorem sector_average_matter_int (J h Δ : ℝ) (ψ : Cfg N → ℂ) (d : Cfg N → ℝ) (t : ℝ) :
    expect (cplx (diagonal d ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))) (lgtStateInt J h Δ ψ t) =
      (1 / 2 ^ N : ℂ) * ∑ q : Cfg N,
        expect (cplx (diagonal d)) (evo (xxzChain J h Δ q) t *ᵥ ψ) := by
  rw [sector_average_int J h Δ ψ _ (gaugeInv_matter d) t]
  simp only [sectorBlock_matter, HsecInt_eq_xxzChain]

/-- **Sector averaging at any number of times** (paper: Remark `rem:interacting` and Corollary
`cor:sector`(b)): for gauge-invariant `O_1, …, O_m` and times `t_1, …, t_m`,
`⟨Ψ₀, O_1(t_1)⋯O_m(t_m) Ψ₀⟩ = 2^{-N} ∑_q ⟨ψ, O_{1,q}(t_1)⋯O_{m,q}(t_m) ψ⟩`. -/
theorem correlator_sector_average_int (J h Δ : ℝ) (ψ : Cfg N → ℂ)
    (Os : List (Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ × ℝ)) (hOs : ∀ p ∈ Os, GaugeInv p.1) :
    expect ((Os.map (fun p => heisInt J h Δ p.1 p.2)).prod) (initZ ψ) =
      (1 / 2 ^ N : ℂ) * ∑ q : Cfg N, expect ((Os.map (fun p =>
        (evo (HsecInt J h Δ q) p.2)ᴴ * sectorBlock p.1 q * evo (HsecInt J h Δ q) p.2)).prod) ψ := by
  rw [← sec_expect_submatrix _ dualEquiv, list_prod_heisInt_dual J h Δ Os hOs,
    sec_expect_blockDiagonal, Finset.mul_sum]
  refine Finset.sum_congr rfl fun q _ => ?_
  have hvec : (fun n => (initZ ψ ∘ dualEquiv) (n, q)) = (((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ • ψ :=
    initZ_dual ψ q
  rw [hvec, sec_expect_smul, sec_norm_sq_prefactor]

/-- The multi-time statement of Corollary `cor:sector`(b) for the model without interaction. -/
theorem correlator_sector_average (J h : ℝ) (ψ : Cfg N → ℂ)
    (Os : List (Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ × ℝ)) (hOs : ∀ p ∈ Os, GaugeInv p.1) :
    expect ((Os.map (fun p => heis J h p.1 p.2)).prod) (initZ ψ) =
      (1 / 2 ^ N : ℂ) * ∑ q : Cfg N, expect ((Os.map (fun p =>
        (evo (Hsec J h q) p.2)ᴴ * sectorBlock p.1 q * evo (Hsec J h q) p.2)).prod) ψ := by
  have h0 := correlator_sector_average_int J h 0 ψ Os hOs
  have e1 : ∀ (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (t : ℝ),
      heisInt J h 0 O t = heis J h O t := fun O t => by
    rw [heisInt, heis, Hint_zero]
  have e2 : ∀ q : Cfg N, HsecInt J h 0 q = Hsec J h q := fun q => by
    rw [HsecInt, zero_smul, add_zero]
  simp only [e1, e2] at h0
  exact h0

namespace InteractingAux

/-- The interacting propagator commutes with every Gauss operator. -/
lemma evo_Hint_comm_gauss (J h Δ : ℝ) (j : Fin N) (t : ℝ) :
    Commute (evo (Hint J h Δ) t) (cplx (gauss j)) := by
  have hc : Commute (cplx (Hint J h Δ)) (cplx (gauss j)) := by
    show cplx (Hint J h Δ) * cplx (gauss j) = cplx (gauss j) * cplx (Hint J h Δ)
    rw [← EnsembleAux.cplx_mul, ← EnsembleAux.cplx_mul, Hint_comm_gauss]
  exact (hc.smul_left _).exp_left

/-- The adjoint interacting propagator commutes with every Gauss operator. -/
lemma evo_Hint_conjTranspose_comm_gauss (J h Δ : ℝ) (j : Fin N) (t : ℝ) :
    Commute (evo (Hint J h Δ) t)ᴴ (cplx (gauss j)) := by
  have h2 := congrArg conjTranspose (evo_Hint_comm_gauss J h Δ j t).eq
  rw [conjTranspose_mul, conjTranspose_mul, EnsembleAux.cplx_gauss_conjTranspose] at h2
  exact h2.symm

/-- The interacting Heisenberg evolution preserves gauge invariance. -/
lemma gaugeInv_heisInt (J h Δ : ℝ) {O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ}
    (hO : GaugeInv O) (t : ℝ) : GaugeInv (heisInt J h Δ O t) := by
  intro j
  exact ((evo_Hint_conjTranspose_comm_gauss J h Δ j t).mul_left (hO j)).mul_left
    (evo_Hint_comm_gauss J h Δ j t)

end InteractingAux

/-- **Superposition equals ensemble** for interacting matter (as Proposition `prop:entropies`(d)). -/
theorem correlator_eq_ensemble_int (J h Δ : ℝ) (n : Cfg N)
    (Os : List (Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ × ℝ)) (hOs : ∀ p ∈ Os, GaugeInv p.1) :
    expect ((Os.map (fun p => heisInt J h Δ p.1 p.2)).prod) (initZ (Pi.single n 1)) =
      trace (ensemble n * (Os.map (fun p => heisInt J h Δ p.1 p.2)).prod) := by
  apply expect_initZ_eq_trace_ensemble
  induction Os with
  | nil =>
    intro j
    simp
  | cons p Os ih =>
    rw [List.map_cons, List.prod_cons]
    exact gaugeInv_mul (gaugeInv_heisInt J h Δ (hOs p List.mem_cons_self) p.2)
      (ih fun p' hp' => hOs p' (List.mem_cons_of_mem _ hp'))

namespace InteractingAux

/-- The interacting sector Hamiltonians are symmetric. -/
lemma HsecInt_isSymm (J h Δ : ℝ) (q : Cfg N) : (HsecInt J h Δ q).IsSymm := by
  unfold HsecInt Hsec nnInt
  exact ((dGamma_isSymm _ (hq_isSymm J h q)).sub (isSymm_one.smul _)).add
    ((isSymm_diagonal _).smul _)

/-- `H_Δ` is symmetric (checked in the dual frame). -/
lemma Hint_isSymm (J h Δ : ℝ) : (Hint (N := N) J h Δ).IsSymm := by
  have hB : (blockDiagonal (fun q => HsecInt (N := N) J h Δ q)).IsSymm := by
    rw [Matrix.IsSymm, blockDiagonal_transpose]
    congr 1
    funext q
    exact (HsecInt_isSymm J h Δ q).eq
  have h2 : Hint (N := N) J h Δ =
      (blockDiagonal (fun q => HsecInt (N := N) J h Δ q)).submatrix dualEquiv.symm
        dualEquiv.symm := by
    rw [← duality_int, submatrix_submatrix, Equiv.self_comp_symm, submatrix_id_id]
  rw [h2]
  exact hB.submatrix _

/-- The initial state `|↑⋯↑⟩_z ⊗ ψ` is normalized when `ψ` is. -/
lemma sum_normSq_initZ (ψ : Cfg N → ℂ) (hψ : ∑ n, Complex.normSq (ψ n) = 1) :
    ∑ p, Complex.normSq (initZ ψ p) = 1 := by
  rw [Fintype.sum_prod_type]
  simp only [initZ, Complex.normSq_mul, EntropiesAux.normSq_prefactor, Finset.sum_const,
    Finset.card_univ, EntropiesAux.card_cfg, nsmul_eq_mul]
  rw [← Finset.mul_sum, ← Finset.mul_sum, hψ]
  push_cast
  field_simp

/-- The interacting quench state is normalized in the dual frame. -/
lemma sum_normSq_lgtStateInt_dual (J h Δ : ℝ) (ψ : Cfg N → ℂ)
    (hψ : ∑ n, Complex.normSq (ψ n) = 1) (t : ℝ) :
    ∑ p, Complex.normSq (lgtStateInt J h Δ ψ t (dualEquiv p)) = 1 := by
  rw [Equiv.sum_comp dualEquiv (fun p => Complex.normSq (lgtStateInt J h Δ ψ t p)), lgtStateInt,
    EntropiesAux.sum_normSq_mulVec _ (evo_unitary _ (Hint_isSymm J h Δ) t), sum_normSq_initZ ψ hψ]

end InteractingAux

/-- **Cut locality for interacting matter** (as Proposition `prop:entropies`(a)): the entanglement of a
left segment in `Ψ(t)` and in `U†Ψ(t)` differ by at most `ln 2`. -/
theorem lgt_duality_defect_int (J h Δ : ℝ) (ψ : Cfg N → ℂ)
    (hψ : ∑ n, Complex.normSq (ψ n) = 1) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    |vonNeumann (rdmLGT ℓ hℓ (lgtStateInt J h Δ ψ t)) (rdmLGT_isHermitian ℓ hℓ _) -
      vonNeumann (rdmLGT ℓ hℓ (fun p => lgtStateInt J h Δ ψ t (dualEquiv p)))
        (rdmLGT_isHermitian ℓ hℓ _)| ≤ Real.log 2 := by
  have h := cutlocal_defect ℓ hℓ (fun p => lgtStateInt J h Δ ψ t (dualEquiv p))
    (sum_normSq_lgtStateInt_dual J h Δ ψ hψ t)
  simpa only [Equiv.apply_symm_apply] using h

end DFL
