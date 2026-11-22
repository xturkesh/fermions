import DFL.Sector

/-!
# The open-chain variant (paper: Proposition `prop:variants`(i))

The chain without boundary link is our chain with the boundary link `b_{N-1}` frozen in the state
`|+⟩` (label `0`): the Hamiltonian never flips the boundary link, so the subspace
`{x_{b_{N-1}} = 0}` is invariant. In the dual frame `x_{b_{N-1}} = ∑_j q_j + ∑_j n_j`, so this
subspace is the parity sector `∏_j q_j = (-1)^{N_f}`. Consequently, for a matter state of definite
parity `π` and links `|↑…↑⟩_z` (the boundary link in `|+⟩`), the sector law is uniform on
`{q : ∑_j q_j = π}` — density `2` with respect to the i.i.d. Bernoulli(1/2) law.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-! ### Elementary facts about `ZMod 2` (private helpers) -/

/-- `a + a = 0` in `ZMod 2`. -/
private lemma var_add_self (a : ZMod 2) : a + a = 0 := by
  revert a; decide

/-- In `ZMod 2`, `a + b = 0 ↔ a = b`. -/
private lemma var_add_eq_zero_iff (a b : ZMod 2) : a + b = 0 ↔ a = b := by
  revert a b; decide

/-- In `ZMod 2`, `a ≠ b ↔ a = b + 1`. -/
private lemma var_ne_iff (a b : ZMod 2) : a ≠ b ↔ a = b + 1 := by
  revert a b; decide

/-- The Hamiltonian never flips the boundary link `b_{N-1}`. -/
theorem H_boundary_link (J h : ℝ) (hN : 0 < N) (p p' : Cfg N × Cfg N)
    (hp : H J h p' p ≠ 0) : p'.2 ⟨N - 1, by omega⟩ = p.2 ⟨N - 1, by omega⟩ := by
  by_contra hne
  apply hp
  have hx : p'.2 ≠ p.2 := fun h' => hne (by rw [h'])
  have h1 : ∀ j : Fin N, hop j p'.1 p.1 * sz j p'.2 p.2 = 0 := by
    intro j
    by_cases hj : j.val + 1 < N
    · have hsz : sz j p'.2 p.2 = 0 := by
        simp only [sz]
        rw [ite_eq_right]
        intro h'
        apply hne
        rw [h', Function.update_of_ne]
        intro h''
        rw [Fin.ext_iff] at h''
        have : (⟨N - 1, by omega⟩ : Fin N).val = N - 1 := rfl
        omega
      rw [hsz, mul_zero]
    · have hhop : hop j = 0 := by
        unfold hop
        rw [dite_eq_right hj]
      rw [hhop, Matrix.zero_apply, zero_mul]
  have h2 : ∀ j : Fin N, (1 : Matrix (Cfg N) (Cfg N) ℝ) p'.1 p.1 * elec j p'.2 p.2 = 0 := by
    intro j
    rw [elec, diagonal_apply_ne _ hx, mul_zero]
  simp only [H, Matrix.sub_apply, Matrix.smul_apply, Matrix.sum_apply, Matrix.kroneckerMap_apply,
    h1, h2, Finset.sum_const_zero, smul_zero, sub_zero]

/-- In the dual frame the boundary-link label is the total parity of charges plus particles. -/
theorem linkOf_last (hN : 0 < N) (n q : Cfg N) :
    linkOf n q ⟨N - 1, by omega⟩ = (∑ j, q j) + ∑ j, n j := by
  simp only [linkOf]
  rw [Finset.filter_true_of_mem, Finset.sum_add_distrib]
  intro k _
  rw [Fin.le_def]
  show k.val ≤ N - 1
  have := k.isLt
  omega

/-- The open-chain initial state: links `0, …, N-2` in `|↑⟩_z`, boundary link in `|+⟩`, matter in
`ψ`: amplitude `2^{-(N-1)/2} ψ(n)` on link configurations with `x_{b_{N-1}} = 0`. -/
def openInit (hN : 0 < N) (ψ : Cfg N → ℂ) : Cfg N × Cfg N → ℂ :=
  fun p => if p.2 ⟨N - 1, by omega⟩ = 0 then (((Real.sqrt 2) ^ (N - 1) : ℝ) : ℂ)⁻¹ * ψ p.1
    else 0

/-- **Open chain, sector law** (paper: Proposition `prop:variants`(i)): for a matter state of
definite parity `π`, every gauge-invariant expectation is the average over the `2^{N-1}` charge
configurations with `∑_j q_j = π` of free-fermion expectations. -/
theorem sector_average_open (J h : ℝ) (hN : 0 < N) (ψ : Cfg N → ℂ) (π : ZMod 2)
    (hψ : ∀ n, ψ n ≠ 0 → ∑ j, n j = π) (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ)
    (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (H J h) t *ᵥ openInit hN ψ) =
      (1 / 2 ^ (N - 1) : ℂ) * ∑ q ∈ univ.filter (fun q : Cfg N => ∑ j, q j = π),
        expect (sectorBlock O q) (evo (Hsec J h q) t *ᵥ ψ) := by
  rw [sec_expect_evo_dual J h _ O hO t, Finset.sum_filter, Finset.mul_sum]
  refine Finset.sum_congr rfl fun q _ => ?_
  have hvec : (fun n => openInit hN ψ (dualEquiv (n, q))) =
      (if ∑ j, q j = π then (((Real.sqrt 2) ^ (N - 1) : ℝ) : ℂ)⁻¹ else 0) • ψ := by
    ext n
    simp only [openInit, dualEquiv_apply, linkOf_last hN, Pi.smul_apply, smul_eq_mul]
    by_cases hn : ψ n = 0
    · simp [hn]
    · rw [hψ n hn]
      by_cases hq : ∑ j, q j = π
      · rw [ite_eq_left ((var_add_eq_zero_iff _ _).2 hq), ite_eq_left hq]
      · rw [ite_eq_right (fun h' => hq ((var_add_eq_zero_iff _ _).1 h')), ite_eq_right hq,
          zero_mul]
  rw [hvec, mulVec_smul, sec_expect_smul]
  split_ifs with hq
  · rw [sec_norm_sq_prefactor]
  · simp

/-- The open-chain sector law has `2^{N-1}` elements (`N ≥ 1`). -/
theorem card_parity_sector (hN : 0 < N) (π : ZMod 2) :
    (univ.filter (fun q : Cfg N => ∑ j, q j = π)).card = 2 ^ (N - 1) := by
  set e : Cfg N := Pi.single (⟨0, hN⟩ : Fin N) 1 with he
  have hsum : ∀ q : Cfg N, ∑ j, (q + e) j = ∑ j, q j + 1 := by
    intro q
    simp only [Pi.add_apply, Finset.sum_add_distrib, he, Fintype.sum_pi_single']
  have hinv : ∀ q : Cfg N, q + e + e = q := by
    intro q
    funext j
    simp only [Pi.add_apply]
    rw [add_assoc, var_add_self, add_zero]
  have hswap : ∀ a : ZMod 2, (univ.filter (fun q : Cfg N => ∑ j, q j = a)).card =
      (univ.filter (fun q : Cfg N => ∑ j, q j = a + 1)).card := by
    intro a
    apply Finset.card_nbij' (· + e) (· + e)
    · intro q hq
      simp only [Finset.mem_coe, Finset.mem_filter, Finset.mem_univ, true_and] at hq ⊢
      rw [hsum, hq]
    · intro q hq
      simp only [Finset.mem_coe, Finset.mem_filter, Finset.mem_univ, true_and] at hq ⊢
      rw [hsum, hq, add_assoc, var_add_self, add_zero]
    · intro q _
      exact hinv q
    · intro q _
      exact hinv q
  have hfil : univ.filter (fun q : Cfg N => ∑ j, q j = π + 1) =
      univ.filter (fun q : Cfg N => ¬ ∑ j, q j = π) :=
    Finset.filter_congr (fun q _ => (var_ne_iff _ _).symm)
  have htot : (univ.filter (fun q : Cfg N => ∑ j, q j = π)).card +
      (univ.filter (fun q : Cfg N => ∑ j, q j = π + 1)).card = 2 ^ N := by
    rw [hfil, Finset.card_filter_add_card_filter_not, Finset.card_univ, Fintype.card_fun,
      ZMod.card, Fintype.card_fin]
  have h2 : 2 ^ N = 2 * 2 ^ (N - 1) := by
    rw [← pow_succ']
    congr 1
    omega
  rw [← hswap π] at htot
  omega

end DFL
