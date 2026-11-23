import DFL.QLM
import DFL.Wilson

/-!
# Remarks on the U(1) models (paper: Remarks `rem:adhh`(a), `rem:wilson`(c), `rem:gausscap`(b))

* `energy_conservation_bound` (`rem:adhh`(a)): if `H = ν 𝒩 + H'` with `ν > 0` and `|⟨H'⟩_v| ≤ C ‖v‖²`, then
  `|⟨𝒩⟩_t - ⟨𝒩⟩_0| ≤ 2C/ν` for every normalized initial state and all `t`, since `ν⟨𝒩⟩_t + ⟨H'⟩_t` is
  conserved.
* `spinPlus`, `spinPlus_opNorm`, `max_spinPlus_weight_even`, `max_spinPlus_weight_odd`, `cS_bounds`
  (`rem:wilson`(c)): for spin `S = n/2` the raising operator `S⁺` has `‖S⁺‖ = max_m ((S-m)(S+m+1))^{1/2}`; the
  maximum is `S(S+1)` for integer `S` and `(S+1/2)²` for half-integer `S`, so
  `1 ≤ c_S = ‖S⁺‖/√(S(S+1)) ≤ 2/√3`.
* `shannon_uniform`, `qlm_numberEntropy_attained` (`rem:gausscap`(b)): the equal-weight superposition of states
  of one sector with `r` distinct values of `N_A` has number entropy `ln r`, so the caps of `prop:gausscap` are
  attained whenever a sector contains such states.
-/

open Matrix Finset Real

noncomputable section

namespace DFL

/-! ### Energy conservation (paper: `rem:adhh`(a)) -/

namespace U1RemarksAux

/-- The real part of an expectation is linear in the real generator: `ν • A + B`. -/
theorem re_expect_cplx_smul_add {ι : Type*} [Fintype ι] (ν : ℝ) (A B : Matrix ι ι ℝ)
    (v : ι → ℂ) :
    (expect (cplx (ν • A + B)) v).re = ν * (expect (cplx A) v).re + (expect (cplx B) v).re := by
  rw [WilsonAux.re_expect_cplx, WilsonAux.re_expect_cplx, WilsonAux.re_expect_cplx]
  simp only [Matrix.add_apply, Matrix.smul_apply, smul_eq_mul, add_mul, Finset.sum_add_distrib,
    Finset.mul_sum, mul_assoc]

end U1RemarksAux

/-- **An elementary bound from energy conservation** (paper: `rem:adhh`(a)). -/
theorem energy_conservation_bound {ι : Type*} [Fintype ι] [DecidableEq ι]
    (Nop Hp : Matrix ι ι ℝ) (hNop : Nop.IsSymm) (hHp : Hp.IsSymm) (ν C : ℝ) (hν : 0 < ν)
    (hC : ∀ v : ι → ℂ, |(expect (cplx Hp) v).re| ≤ C * ∑ i, Complex.normSq (v i))
    (ψ : ι → ℂ) (hψ : ∑ i, Complex.normSq (ψ i) = 1) (t : ℝ) :
    |(expect (cplx Nop) (evo (ν • Nop + Hp) t *ᵥ ψ)).re - (expect (cplx Nop) ψ).re| ≤ 2 * C / ν := by
  have hA : (ν • Nop + Hp).IsSymm := (hNop.smul ν).add hHp
  -- energy conservation: `ν ⟨𝒩⟩ + ⟨H'⟩` is constant along the evolution
  have hre := congrArg Complex.re (WilsonAux.expect_evo (ν • Nop + Hp) hA t ψ)
  rw [U1RemarksAux.re_expect_cplx_smul_add, U1RemarksAux.re_expect_cplx_smul_add] at hre
  -- the evolved state stays normalized
  have hφn : ∑ i, Complex.normSq ((evo (ν • Nop + Hp) t *ᵥ ψ) i) = 1 := by
    rw [WilsonAux.sum_normSq_mulVec _ (evo_unitary _ hA t), hψ]
  have h1 := hC (evo (ν • Nop + Hp) t *ᵥ ψ)
  have h2 := hC ψ
  rw [hφn, mul_one] at h1
  rw [hψ, mul_one] at h2
  have hb : |ν * ((expect (cplx Nop) (evo (ν • Nop + Hp) t *ᵥ ψ)).re -
      (expect (cplx Nop) ψ).re)| ≤ 2 * C := by
    rw [abs_le] at h1 h2 ⊢
    constructor <;> nlinarith
  rw [abs_mul, abs_of_pos hν] at hb
  rw [le_div_iff₀ hν]
  linarith

/-! ### The norm of the spin raising operator (paper: `rem:wilson`(c)) -/

/-- The spin raising operator for spin `S = n/2` on the `n + 1` states `m = -S + k`, `k = 0, …, n`:
`S⁺ |k⟩ = ((n - k)(k + 1))^{1/2} |k + 1⟩`. -/
def spinPlus (n : ℕ) : Matrix (Fin (n + 1)) (Fin (n + 1)) ℝ := fun i j =>
  if i.val = j.val + 1 then Real.sqrt (((n - j.val : ℕ) : ℝ) * ((j.val + 1 : ℕ) : ℝ)) else 0

namespace U1RemarksAux

/-- `S⁺ᴴ S⁺` is diagonal with entries `(n - j)(j + 1)`. -/
theorem spinPlus_conjTranspose_mul_self (n : ℕ) :
    (spinPlus n)ᴴ * spinPlus n =
      diagonal (fun j : Fin (n + 1) => ((n - j.val : ℕ) : ℝ) * ((j.val + 1 : ℕ) : ℝ)) := by
  ext i j
  rw [Matrix.mul_apply, diagonal_apply]
  simp only [conjTranspose_apply, star_trivial, spinPlus]
  by_cases hij : i = j
  · subst hij
    rw [ite_eq_left rfl]
    by_cases hi : i.val < n
    · -- the only nonzero term is `k = i + 1`
      rw [Finset.sum_eq_single_of_mem (⟨i.val + 1, by omega⟩ : Fin (n + 1)) (Finset.mem_univ _)]
      · simp only [ite_true]
        rw [Real.mul_self_sqrt (by positivity)]
      · intro k _ hk
        have : k.val ≠ i.val + 1 := fun h => hk (Fin.ext h)
        simp [this]
    · -- top state `i = n`: the column of `S⁺` vanishes, and so does the weight
      have hi' : i.val = n := by omega
      rw [Finset.sum_eq_zero]
      · simp [hi']
      · intro k _
        have : k.val ≠ i.val + 1 := by omega
        simp [this]
  · rw [ite_eq_right hij]
    apply Finset.sum_eq_zero
    intro k _
    by_cases h1 : k.val = i.val + 1
    · have h2 : k.val ≠ j.val + 1 := fun h => hij (Fin.ext (by omega))
      simp [h2]
    · simp [h1]

/-- The sup norm of the weights `(n - j)(j + 1)`, `j ≤ n`, is their maximum over `j < n`. -/
theorem norm_spinWeights (n : ℕ) (hn : 0 < n) :
    ‖(fun j : Fin (n + 1) => ((n - j.val : ℕ) : ℝ) * ((j.val + 1 : ℕ) : ℝ))‖ =
      (Finset.range n).sup' (Finset.nonempty_range_iff.mpr (by omega))
        (fun k => (((n - k : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ))) := by
  have hM0 : 0 ≤ (Finset.range n).sup' (Finset.nonempty_range_iff.mpr (by omega))
        (fun k => (((n - k : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ))) :=
    le_trans (by positivity) (Finset.le_sup' (fun k => (((n - k : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ)))
      (Finset.mem_range.mpr hn))
  apply le_antisymm
  · rw [pi_norm_le_iff_of_nonneg hM0]
    intro j
    rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
    by_cases hj : j.val < n
    · exact Finset.le_sup' (fun k => (((n - k : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ)))
        (Finset.mem_range.mpr hj)
    · have : n - j.val = 0 := by omega
      rw [this]
      simpa using hM0
  · apply Finset.sup'_le
    intro k hk
    rw [Finset.mem_range] at hk
    have h := norm_le_pi_norm
      (fun j : Fin (n + 1) => ((n - j.val : ℕ) : ℝ) * ((j.val + 1 : ℕ) : ℝ)) ⟨k, by omega⟩
    rw [Real.norm_eq_abs] at h
    exact le_trans (le_abs_self _) h

end U1RemarksAux

open scoped Matrix.Norms.L2Operator in
/-- `‖S⁺‖ = max_k ((n - k)(k + 1))^{1/2}` (operator norm). -/
theorem spinPlus_opNorm (n : ℕ) (hn : 0 < n) :
    ‖spinPlus n‖ = Real.sqrt ((Finset.range n).sup' (Finset.nonempty_range_iff.mpr (by omega))
      (fun k => (((n - k : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ)))) := by
  -- C*-identity `‖S⁺‖² = ‖S⁺ᴴ S⁺‖`, and `S⁺ᴴ S⁺` is diagonal
  have h1 : ‖spinPlus n‖ * ‖spinPlus n‖ = (Finset.range n).sup'
      (Finset.nonempty_range_iff.mpr (by omega))
      (fun k => (((n - k : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ))) := by
    rw [← l2_opNorm_conjTranspose_mul_self, U1RemarksAux.spinPlus_conjTranspose_mul_self,
      l2_opNorm_diagonal, U1RemarksAux.norm_spinWeights n hn]
  rw [← h1, Real.sqrt_mul_self (norm_nonneg _)]

/-- Integer spin `S = m`: `max_k (2m - k)(k + 1) = m (m + 1) = S(S+1)`. -/
theorem max_spinPlus_weight_even (m : ℕ) (hm : 0 < m) :
    (Finset.range (2 * m)).sup' (Finset.nonempty_range_iff.mpr (by omega))
      (fun k => (((2 * m - k : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ))) = (m : ℝ) * (m + 1) := by
  apply le_antisymm
  · -- `m (m + 1) - (2m - k)(k + 1) = (m - k)(m - k - 1) ≥ 0` for integers
    apply Finset.sup'_le
    intro k hk
    rw [Finset.mem_range] at hk
    rw [Nat.cast_sub hk.le]
    push_cast
    have key : (0 : ℤ) ≤ ((m : ℤ) - k) * ((m : ℤ) - k - 1) := by
      rcases le_or_gt (k : ℤ) ((m : ℤ) - 1) with h | h
      · exact mul_nonneg (by omega) (by omega)
      · exact mul_nonneg_of_nonpos_of_nonpos (by omega) (by omega)
    have key' : (0 : ℝ) ≤ ((m : ℝ) - k) * ((m : ℝ) - k - 1) := by exact_mod_cast key
    nlinarith
  · -- attained at `k = m`
    refine le_trans ?_ (Finset.le_sup' _ (Finset.mem_range.mpr (by omega : m < 2 * m)))
    show (m : ℝ) * (m + 1) ≤ ((2 * m - m : ℕ) : ℝ) * ((m + 1 : ℕ) : ℝ)
    rw [show 2 * m - m = m by omega]
    push_cast
    exact le_rfl

/-- Half-integer spin `S = m + 1/2`: `max_k (2m + 1 - k)(k + 1) = (m + 1)² = (S + 1/2)²`. -/
theorem max_spinPlus_weight_odd (m : ℕ) :
    (Finset.range (2 * m + 1)).sup' (Finset.nonempty_range_iff.mpr (by omega))
      (fun k => (((2 * m + 1 - k : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ))) = ((m : ℝ) + 1) ^ 2 := by
  apply le_antisymm
  · -- `(m + 1)² - (2m + 1 - k)(k + 1) = (m - k)² ≥ 0`
    apply Finset.sup'_le
    intro k hk
    rw [Finset.mem_range] at hk
    rw [Nat.cast_sub hk.le]
    push_cast
    nlinarith [sq_nonneg ((m : ℝ) - k)]
  · -- attained at `k = m`
    refine le_trans ?_ (Finset.le_sup' _ (Finset.mem_range.mpr (by omega : m < 2 * m + 1)))
    show ((m : ℝ) + 1) ^ 2 ≤ ((2 * m + 1 - m : ℕ) : ℝ) * ((m + 1 : ℕ) : ℝ)
    rw [show 2 * m + 1 - m = m + 1 by omega]
    push_cast
    nlinarith

/-- `1 ≤ c_S ≤ 2/√3` for every spin `S = n/2 > 0`, where `c_S² = max_k (n - k)(k + 1) / (S(S+1))`. -/
theorem cS_bounds (n : ℕ) (hn : 0 < n) :
    1 ≤ (Finset.range n).sup' (Finset.nonempty_range_iff.mpr (by omega))
        (fun k => (((n - k : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ))) / ((n : ℝ) / 2 * ((n : ℝ) / 2 + 1)) ∧
      (Finset.range n).sup' (Finset.nonempty_range_iff.mpr (by omega))
        (fun k => (((n - k : ℕ) : ℝ) * ((k + 1 : ℕ) : ℝ))) / ((n : ℝ) / 2 * ((n : ℝ) / 2 + 1)) ≤
        4 / 3 := by
  obtain ⟨m, rfl | rfl⟩ := Nat.even_or_odd' n
  · -- integer spin: the ratio is exactly `1`
    have hm : 0 < m := by omega
    rw [max_spinPlus_weight_even m hm]
    have hden : ((2 * m : ℕ) : ℝ) / 2 * (((2 * m : ℕ) : ℝ) / 2 + 1) = (m : ℝ) * (m + 1) := by
      push_cast
      ring
    have hpos : (0 : ℝ) < (m : ℝ) * (m + 1) := by positivity
    rw [hden, div_self hpos.ne']
    norm_num
  · -- half-integer spin: the ratio is `4(m+1)² / ((2m+1)(2m+3)) ∈ [1, 4/3]`
    rw [max_spinPlus_weight_odd m]
    have hden : ((2 * m + 1 : ℕ) : ℝ) / 2 * (((2 * m + 1 : ℕ) : ℝ) / 2 + 1) =
        (2 * (m : ℝ) + 1) * (2 * m + 3) / 4 := by
      push_cast
      ring
    have hpos : (0 : ℝ) < (2 * (m : ℝ) + 1) * (2 * m + 3) / 4 := by positivity
    rw [hden]
    have hm0 : (0 : ℝ) ≤ m := Nat.cast_nonneg m
    constructor
    · rw [le_div_iff₀ hpos]
      nlinarith
    · rw [div_le_iff₀ hpos]
      nlinarith

/-! ### The caps of `prop:gausscap` are attained (paper: `rem:gausscap`(b)) -/

/-- The uniform distribution on `r` values has Shannon entropy `ln r`. -/
theorem shannon_uniform (M : ℕ) (K : Finset ℕ) (hK : ∀ k ∈ K, k ≤ M) (hne : K.Nonempty) :
    shannon M (fun k => if k ∈ K then 1 / (K.card : ℝ) else 0) = Real.log K.card := by
  unfold shannon
  have h1 : ∀ k ∈ range (M + 1), negMulLog (if k ∈ K then 1 / (K.card : ℝ) else 0) =
      if k ∈ K then negMulLog (1 / (K.card : ℝ)) else 0 := by
    intro k _
    split_ifs <;> simp
  rw [Finset.sum_congr rfl h1, ← Finset.sum_filter, Finset.sum_const]
  have hfilt : (range (M + 1)).filter (fun k => k ∈ K) = K := by
    ext k
    simp only [Finset.mem_filter, Finset.mem_range, and_iff_right_iff_imp]
    intro hk
    have := hK k hk
    omega
  rw [hfilt, nsmul_eq_mul, Real.negMulLog, one_div, Real.log_inv]
  have hr : (0 : ℝ) < K.card := by exact_mod_cast hne.card_pos
  field_simp

namespace U1RemarksAux

/-- If at most one term of a finite sum is nonzero, `normSq` of the sum is the sum of `normSq`s. -/
theorem normSq_sum_of_unique {α : Type*} (K : Finset α) (g : α → ℂ)
    (hg : ∀ k ∈ K, ∀ k' ∈ K, g k ≠ 0 → g k' ≠ 0 → k = k') :
    Complex.normSq (∑ k ∈ K, g k) = ∑ k ∈ K, Complex.normSq (g k) := by
  by_cases h : ∃ k ∈ K, g k ≠ 0
  · obtain ⟨k₀, hk₀, hne⟩ := h
    have hz : ∀ k ∈ K, k ≠ k₀ → g k = 0 := by
      intro k hk hkk
      by_contra hgk
      exact hkk (hg k hk k₀ hk₀ hgk hne)
    rw [Finset.sum_eq_single_of_mem k₀ hk₀ hz,
      Finset.sum_eq_single_of_mem k₀ hk₀
        (fun k hk hkk => by rw [hz k hk hkk, Complex.normSq_zero])]
  · push Not at h
    rw [Finset.sum_eq_zero h,
      Finset.sum_eq_zero (fun k hk => by rw [h k hk, Complex.normSq_zero]), Complex.normSq_zero]

end U1RemarksAux

/-- **The caps are attained** (paper: `rem:gausscap`(b)): if `φ_k`, `k ∈ K`, are normalized states with
`N_{[a,b]} φ_k = k φ_k`, their equal-weight superposition is normalized and has number entropy
`ln |K|`. -/
theorem qlm_numberEntropy_attained {N : ℕ} (S : ℕ) (a b : ℕ) (K : Finset ℕ)
    (hK : ∀ k ∈ K, k ≤ N) (hne : K.Nonempty)
    (φ : ℕ → Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ)
    (hφn : ∀ k ∈ K, ∑ c, Complex.normSq (φ k c) = 1)
    (hφs : ∀ k ∈ K, ∀ c, φ k c ≠ 0 → countInt a b c.1 = k) :
    (∑ c, Complex.normSq ((Real.sqrt (K.card : ℝ))⁻¹ * ∑ k ∈ K, φ k c) = 1) ∧
      shannon N (qlmNumDist S a b (fun c => (Real.sqrt (K.card : ℝ))⁻¹ * ∑ k ∈ K, φ k c)) =
        Real.log K.card := by
  have hr : (0 : ℝ) < K.card := by exact_mod_cast hne.card_pos
  have hs2 : Complex.normSq (((Real.sqrt (K.card : ℝ))⁻¹ : ℝ) : ℂ) = 1 / (K.card : ℝ) := by
    rw [Complex.normSq_ofReal, ← mul_inv, Real.mul_self_sqrt hr.le, one_div]
  -- the `φ_k` have disjoint supports (different particle numbers)
  have huniq : ∀ c, ∀ k ∈ K, ∀ k' ∈ K, φ k c ≠ 0 → φ k' c ≠ 0 → k = k' := by
    intro c k hk k' hk' h1 h2
    have e1 := hφs k hk c h1
    have e2 := hφs k' hk' c h2
    exact_mod_cast e1.symm.trans e2
  have hnormSq : ∀ c, Complex.normSq ((Real.sqrt (K.card : ℝ))⁻¹ * ∑ k ∈ K, φ k c) =
      1 / (K.card : ℝ) * ∑ k ∈ K, Complex.normSq (φ k c) := by
    intro c
    rw [Complex.normSq_mul, hs2, U1RemarksAux.normSq_sum_of_unique K _ (huniq c)]
  constructor
  · simp_rw [hnormSq]
    rw [← Finset.mul_sum, Finset.sum_comm, Finset.sum_congr rfl hφn, Finset.sum_const,
      nsmul_eq_mul, mul_one]
    field_simp
  · -- the particle-number distribution is uniform on `K`
    have hdist : qlmNumDist S a b (fun c => (Real.sqrt (K.card : ℝ))⁻¹ * ∑ k ∈ K, φ k c) =
        fun j => if j ∈ K then 1 / (K.card : ℝ) else 0 := by
      funext j
      unfold qlmNumDist
      simp_rw [hnormSq]
      rw [← Finset.mul_sum, Finset.sum_comm]
      have hk : ∀ k ∈ K, ∑ c ∈ univ.filter
          (fun c : Cfg N × (Fin N → Fin (2 * S + 1)) => countInt a b c.1 = j),
          Complex.normSq (φ k c) = if k = j then 1 else 0 := by
        intro k hk
        split_ifs with hkj
        · subst hkj
          rw [Finset.sum_filter_of_ne]
          · exact hφn k hk
          · intro c _ hc
            exact hφs k hk c (fun h => hc (by rw [h, Complex.normSq_zero]))
        · apply Finset.sum_eq_zero
          intro c hc
          rw [Finset.mem_filter] at hc
          by_contra hne'
          have hφ : φ k c ≠ 0 := fun h => hne' (by rw [h, Complex.normSq_zero])
          have e := hφs k hk c hφ
          rw [hc.2] at e
          exact hkj (by exact_mod_cast e.symm)
      rw [Finset.sum_congr rfl hk, Finset.sum_ite_eq']
      split_ifs <;> simp
    rw [hdist]
    exact shannon_uniform N K hK hne

end DFL
