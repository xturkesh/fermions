import DFL.LRNorm

/-!
# The constants of (EC) degenerate at weak disorder (paper: Remark `rem:rate`(iii))

* `corrQ_free_mid`, `tendsto_corrQ_free_mid`: for the free chain on `N = 2m + 1` sites (0-indexed), the
  eigenfunction correlator between the first site and the middle site is
  `Q_{Δ_N}(1, (N+1)/2) = 2 / ((N+1) sin(π/(N+1))) → 2/π`: with the sine eigenvectors
  `φ_k(x) = (2/(N+1))^{1/2} sin(πkx/(N+1))` only odd `k` contribute, and
  `∑_{k odd} sin(πk/(N+1)) = 1/sin(π/(N+1))`.
* `corrQ_end_lower`: **Chebyshev lower bound**, for every charge configuration `q` on `m + 1` sites,
  `Q_q(1, m + 1) ≥ ½ (|J|/(|J| + |h|))^m = ½ (1 + |λ|/2)^{-m}`: the Chebyshev polynomial `T_m(E/(2|J| + 2|h|))` is
  bounded by one on the spectrum, and its `(1, m+1)` matrix element is `2^{m-1}(-J)^m/(2|J| + 2|h|)^m`, because
  `h_q^k` vanishes beyond the `k`-th superdiagonal and equals `(-J)^k` on it.
* `ec_rate_le`: hence every rate `η` admissible in (EC) satisfies `η ≤ ln(1 + |λ|/2) ≤ |λ|/2`, and
  `no_uniform_EC_weak`: no constants `C, η > 0` work uniformly for `0 < λ ≤ λ₁`.
-/

open Filter Topology Real

noncomputable section

namespace DFL

namespace WeakDisorderAux

open Finset

variable {N : ℕ}

/-! ### The free chain: the end-to-middle correlator -/

/-- Telescoping: `2 sin θ ∑_{k<2j+1} sin((k+1)θ) |sin((k+1)π/2)| = 1 - cos((2j+2)θ)`. -/
lemma two_sin_mul_sum (θ : ℝ) (j : ℕ) :
    2 * Real.sin θ * ∑ k ∈ range (2 * j + 1),
      Real.sin (((k : ℝ) + 1) * θ) * |Real.sin (((k : ℝ) + 1) * (π / 2))| =
      1 - Real.cos ((2 * (j : ℝ) + 2) * θ) := by
  induction j with
  | zero =>
    simp only [mul_zero, zero_add, range_one, sum_singleton, Nat.cast_zero, one_mul,
      Real.sin_pi_div_two, abs_one, mul_one]
    have h2 : Real.cos (2 * θ) = 1 - 2 * Real.sin θ ^ 2 := by
      rw [Real.cos_two_mul, Real.cos_sq']; ring
    rw [h2]; ring
  | succ j ih =>
    have hs : 2 * (j + 1) + 1 = (2 * j + 1) + 1 + 1 := by ring
    have t1 : |Real.sin ((((2 * j + 1 : ℕ) : ℝ) + 1) * (π / 2))| = 0 := by
      rw [show (((2 * j + 1 : ℕ) : ℝ) + 1) * (π / 2) = ((j + 1 : ℕ) : ℝ) * π by push_cast; ring,
        Real.sin_nat_mul_pi, abs_zero]
    have t2 : |Real.sin ((((2 * j + 1 + 1 : ℕ) : ℝ) + 1) * (π / 2))| = 1 := by
      rw [show (((2 * j + 1 + 1 : ℕ) : ℝ) + 1) * (π / 2) = ((j + 1 : ℕ) : ℝ) * π + π / 2 by
          push_cast; ring,
        Real.sin_add_pi_div_two, Real.cos_nat_mul_pi, abs_pow, abs_neg, abs_one, one_pow]
    have t3 : Real.sin ((((2 * j + 1 + 1 : ℕ) : ℝ) + 1) * θ) =
        Real.sin ((2 * (j : ℝ) + 3) * θ) := by
      congr 1; push_cast; ring
    rw [hs, sum_range_succ, sum_range_succ, t1, t2, t3, mul_zero, add_zero, mul_one, mul_add, ih]
    have e4 : Real.cos ((2 * (j : ℝ) + 2) * θ) = Real.cos ((2 * (j : ℝ) + 3) * θ - θ) := by
      congr 1; ring
    have e5 : Real.cos ((2 * ((j + 1 : ℕ) : ℝ) + 2) * θ) =
        Real.cos ((2 * (j : ℝ) + 3) * θ + θ) := by
      congr 1; push_cast; ring
    rw [e4, e5, Real.cos_sub, Real.cos_add]
    ring

/-- `∑_{k<2m+1} sin((k+1)θ) |sin((k+1)π/2)| = 1 / sin θ` for `θ = π/(2m+2)`. -/
lemma sum_free_mid (m : ℕ) :
    ∑ k ∈ range (2 * m + 1), Real.sin (((k : ℝ) + 1) * (π / (2 * (m : ℝ) + 2))) *
      |Real.sin (((k : ℝ) + 1) * (π / 2))| = 1 / Real.sin (π / (2 * (m : ℝ) + 2)) := by
  have hm2 : (0 : ℝ) < 2 * (m : ℝ) + 2 := by positivity
  have hpos : 0 < Real.sin (π / (2 * (m : ℝ) + 2)) := by
    apply Real.sin_pos_of_pos_of_lt_pi (by positivity)
    rw [div_lt_iff₀ hm2]
    nlinarith [Real.pi_pos, (Nat.cast_nonneg m : (0 : ℝ) ≤ m)]
  have h := two_sin_mul_sum (π / (2 * (m : ℝ) + 2)) m
  have hc : (2 * (m : ℝ) + 2) * (π / (2 * (m : ℝ) + 2)) = π := by field_simp
  rw [hc, Real.cos_pi] at h
  rw [eq_div_iff hpos.ne']
  linear_combination h / 2

/-- The sine mode `φ_k` at the first site of the chain with `2m+1` sites. -/
lemma sineMode_first (m : ℕ) (κ : Fin (2 * m + 1)) :
    LRNormAux.sineMode κ ⟨0, by omega⟩ = √(2 / (2 * (m : ℝ) + 2)) *
      Real.sin ((((κ : ℕ) : ℝ) + 1) * (π / (2 * (m : ℝ) + 2))) := by
  have hN : ((2 * m + 1 : ℕ) : ℝ) + 1 = 2 * (m : ℝ) + 2 := by push_cast; ring
  simp only [LRNormAux.sineMode, LRNormAux.sineAngle, hN]
  congr 2
  simp only [Nat.cast_zero, zero_add, mul_one]
  ring

/-- The sine mode `φ_k` at the middle site of the chain with `2m+1` sites. -/
lemma sineMode_mid (m : ℕ) (κ : Fin (2 * m + 1)) :
    LRNormAux.sineMode κ ⟨m, by omega⟩ = √(2 / (2 * (m : ℝ) + 2)) *
      Real.sin ((((κ : ℕ) : ℝ) + 1) * (π / 2)) := by
  have hN : ((2 * m + 1 : ℕ) : ℝ) + 1 = 2 * (m : ℝ) + 2 := by push_cast; ring
  have hm2 : (2 * (m : ℝ) + 2) ≠ 0 := by positivity
  simp only [LRNormAux.sineMode, LRNormAux.sineAngle, hN]
  congr 2
  field_simp

/-! ### Polynomials of a symmetric matrix and the correlator -/

/-- Spectral form of matrix powers: `(A^k)_{xy} = ∑_α E_α^k φ_α(x) φ_α(y)`. -/
lemma pow_apply_eq_sum (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (k : ℕ)
    (x y : Fin N) :
    (A ^ k) x y = ∑ α, hA.eigenvalues α ^ k *
      (hA.eigenvectorBasis α x * hA.eigenvectorBasis α y) := by
  induction k generalizing x with
  | zero =>
    simp only [pow_zero, one_mul]
    rw [NormalizationAux.eig_compl A hA x y, Matrix.one_apply]
  | succ k ih =>
    rw [pow_succ', Matrix.mul_apply]
    simp_rw [ih, Finset.mul_sum]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun α _ => ?_
    have e := NormalizationAux.eig_eq A hA α x
    calc ∑ z, A x z * (hA.eigenvalues α ^ k *
          (hA.eigenvectorBasis α z * hA.eigenvectorBasis α y))
        = (∑ z, A x z * hA.eigenvectorBasis α z) *
            (hA.eigenvalues α ^ k * hA.eigenvectorBasis α y) := by
          rw [Finset.sum_mul]
          exact Finset.sum_congr rfl fun z _ => by ring
      _ = hA.eigenvalues α ^ (k + 1) * (hA.eigenvectorBasis α x * hA.eigenvectorBasis α y) := by
          rw [e]; ring

/-- A polynomial bounded by one on the spectrum is dominated by the correlator:
`|∑_k c_k (A^k)_{xy}| ≤ Q(x, y)` if `|∑_k c_k E_α^k| ≤ 1` for every eigenvalue `E_α`. -/
lemma abs_sum_pow_le_corrQ (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (c : ℕ → ℝ)
    (n : ℕ) (hc : ∀ α, |∑ k ∈ range n, c k * hA.eigenvalues α ^ k| ≤ 1) (x y : Fin N) :
    |∑ k ∈ range n, c k * (A ^ k) x y| ≤ corrQ A hA x y := by
  have e : ∑ k ∈ range n, c k * (A ^ k) x y = ∑ α, (∑ k ∈ range n, c k * hA.eigenvalues α ^ k) *
      (hA.eigenvectorBasis α x * hA.eigenvectorBasis α y) := by
    simp_rw [pow_apply_eq_sum A hA, Finset.mul_sum, Finset.sum_mul]
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun α _ => Finset.sum_congr rfl fun k _ => by ring
  rw [e, corrQ]
  refine (Finset.abs_sum_le_sum_abs _ _).trans (Finset.sum_le_sum fun α _ => ?_)
  rw [abs_mul, abs_mul]
  have h0 : 0 ≤ |hA.eigenvectorBasis α x| * |hA.eigenvectorBasis α y| := by positivity
  calc |∑ k ∈ range n, c k * hA.eigenvalues α ^ k| *
        (|hA.eigenvectorBasis α x| * |hA.eigenvectorBasis α y|)
      ≤ 1 * (|hA.eigenvectorBasis α x| * |hA.eigenvectorBasis α y|) :=
        mul_le_mul_of_nonneg_right (hc α) h0
    _ = |hA.eigenvectorBasis α x| * |hA.eigenvectorBasis α y| := one_mul _

/-- Powers of `h_q` vanish beyond the `k`-th superdiagonal, where they equal `(-J)^k`. -/
lemma hq_pow_apply (J h : ℝ) (q : Cfg N) (x : Fin N) (k : ℕ) (y : Fin N) :
    (x.val + k < y.val → (hq J h q ^ k) x y = 0) ∧
      (x.val + k = y.val → (hq J h q ^ k) x y = (-J) ^ k) := by
  induction k generalizing y with
  | zero =>
    refine ⟨fun hxy => ?_, fun hxy => ?_⟩
    · rw [pow_zero, Matrix.one_apply_ne]
      intro e
      rw [e] at hxy
      omega
    · have e : x = y := Fin.ext (by omega)
      rw [e, pow_zero, pow_zero, Matrix.one_apply_eq]
  | succ k ih =>
    rw [pow_succ (hq J h q) k, Matrix.mul_apply]
    refine ⟨fun hxy => Finset.sum_eq_zero fun z _ => ?_, fun hxy => ?_⟩
    · by_cases hz : z.val + 1 < y.val
      · rw [ImbalanceAux.hq_apply_far J h q z y hz, mul_zero]
      · rw [(ih z).1 (by omega), zero_mul]
    · have hz0 : x.val + k < N := by have := y.isLt; omega
      rw [Finset.sum_eq_single ⟨x.val + k, hz0⟩]
      · rw [(ih _).2 rfl, ImbalanceAux.hq_apply_succ J h q _ y (by show x.val + k + 1 = y.val; omega)]
        ring
      · intro z _ hz
        by_cases hz' : z.val + 1 < y.val
        · rw [ImbalanceAux.hq_apply_far J h q z y hz', mul_zero]
        · have hne : z.val ≠ x.val + k := fun e => hz (Fin.ext e)
          rw [(ih z).1 (by omega), zero_mul]
      · intro hmem
        exact absurd (Finset.mem_univ _) hmem

/-- **Chebyshev lower bound**: on `m + 1` sites, `Q_q(0, m) ≥ (1/2) (|J|/(|J|+|h|))^m` for every
charge configuration `q` (test function `T_m(E/(2|J| + 2|h|))`, bounded by one on the spectrum). -/
lemma corrQ_end_ge (J h : ℝ) (hJ : J ≠ 0) (m : ℕ) (hm : 1 ≤ m) (q : Cfg (m + 1)) :
    1 / 2 * (|J| / (|J| + |h|)) ^ m ≤
      corrQ (hq J h q) (hq_isHermitian J h q) ⟨0, by omega⟩ ⟨m, by omega⟩ := by
  have hJp : 0 < |J| := abs_pos.mpr hJ
  have hS : 0 < |J| + |h| := by have := abs_nonneg h; linarith
  have hRp : 0 < 2 * |J| + 2 * |h| := by linarith
  have hdeg : (Polynomial.Chebyshev.T ℝ (m : ℤ)).natDegree = m := by
    rw [Polynomial.Chebyshev.natDegree_T]; simp
  have hlead : (Polynomial.Chebyshev.T ℝ (m : ℤ)).coeff m = 2 ^ (m - 1) := by
    have h1 := Polynomial.coeff_natDegree (p := Polynomial.Chebyshev.T ℝ (m : ℤ))
    rw [hdeg, Polynomial.Chebyshev.leadingCoeff_T] at h1
    simpa using h1
  set c : ℕ → ℝ := fun k => (Polynomial.Chebyshev.T ℝ (m : ℤ)).coeff k /
    (2 * |J| + 2 * |h|) ^ k with hc
  have hspec : ∀ α, |∑ k ∈ range (m + 1), c k * (hq_isHermitian J h q).eigenvalues α ^ k| ≤ 1 := by
    intro α
    have e : ∑ k ∈ range (m + 1), c k * (hq_isHermitian J h q).eigenvalues α ^ k =
        (Polynomial.Chebyshev.T ℝ (m : ℤ)).eval
          ((hq_isHermitian J h q).eigenvalues α / (2 * |J| + 2 * |h|)) := by
      rw [Polynomial.eval_eq_sum_range' (n := m + 1) (by rw [hdeg]; omega)]
      refine Finset.sum_congr rfl fun k _ => ?_
      simp only [hc, div_pow]
      ring
    rw [e]
    apply Polynomial.Chebyshev.abs_eval_T_real_le_one
    rw [abs_div, abs_of_pos hRp, div_le_one hRp]
    exact abs_eigenvalue_hq_le J h q α
  have hsum : ∑ k ∈ range (m + 1), c k * (hq J h q ^ k) ⟨0, by omega⟩ ⟨m, by omega⟩ =
      c m * (-J) ^ m := by
    rw [Finset.sum_range_succ, (hq_pow_apply J h q ⟨0, by omega⟩ m ⟨m, by omega⟩).2 (by simp),
      Finset.sum_eq_zero, zero_add]
    intro k hk
    rw [(hq_pow_apply J h q ⟨0, by omega⟩ k ⟨m, by omega⟩).1
      (by simp only [Finset.mem_range] at hk; simp only [zero_add]; omega), mul_zero]
  have key := abs_sum_pow_le_corrQ (hq J h q) (hq_isHermitian J h q) c (m + 1) hspec
    ⟨0, by omega⟩ ⟨m, by omega⟩
  rw [hsum] at key
  refine le_trans (le_of_eq ?_) key
  simp only [hc, hlead, abs_mul, abs_div, abs_pow, abs_neg, abs_two, abs_of_pos hRp]
  obtain ⟨n, rfl⟩ : ∃ n, m = n + 1 := ⟨m - 1, by omega⟩
  rw [Nat.add_sub_cancel, show 2 * |J| + 2 * |h| = 2 * (|J| + |h|) by ring, mul_pow, div_pow,
    pow_succ (2 : ℝ) n]
  field_simp

end WeakDisorderAux

open WeakDisorderAux in
/-- **The free-chain correlator between the end and the middle** (paper: `rem:rate`(iii)): for `N = 2m+1`,
`Q_{Δ_N}(1, (N+1)/2) = 2 / ((N+1) sin(π/(N+1)))`. -/
theorem corrQ_free_mid (J : ℝ) (hJ : J ≠ 0) (m : ℕ) :
    corrQ (hq J 0 (0 : Cfg (2 * m + 1))) (hq_isHermitian J 0 0) ⟨0, by omega⟩ ⟨m, by omega⟩ =
      2 / ((2 * (m : ℝ) + 2) * Real.sin (π / (2 * (m : ℝ) + 2))) := by
  set hA := hq_isHermitian J 0 (0 : Cfg (2 * m + 1)) with hAdef
  have hinj := hq_eigenvalues_injective J 0 hJ (0 : Cfg (2 * m + 1))
  have key : ∀ κ : Fin (2 * m + 1), ∃ β, hA.eigenvalues β = LRNormAux.sineEig J 0 κ ∧
      ∀ x, |LRNormAux.sineMode κ x| = |hA.eigenvectorBasis β x| := fun κ =>
    NormalizationAux.exists_abs_eq_eigvec _ hA hinj _ _ (LRNormAux.sineMode_eig J 0 κ)
      (LRNormAux.sineMode_norm κ)
  choose σ hσ1 hσ2 using key
  have hσinj : Function.Injective σ := by
    intro κ κ' e
    exact LRNormAux.sineEig_injective J 0 hJ (by rw [← hσ1 κ, ← hσ1 κ', e])
  have hm2 : (0 : ℝ) < 2 * (m : ℝ) + 2 := by positivity
  have hterm : ∀ κ : Fin (2 * m + 1),
      |hA.eigenvectorBasis (σ κ) ⟨0, by omega⟩| * |hA.eigenvectorBasis (σ κ) ⟨m, by omega⟩| =
      2 / (2 * (m : ℝ) + 2) * (Real.sin ((((κ : ℕ) : ℝ) + 1) * (π / (2 * (m : ℝ) + 2))) *
        |Real.sin ((((κ : ℕ) : ℝ) + 1) * (π / 2))|) := by
    intro κ
    rw [← hσ2 κ ⟨0, by omega⟩, ← hσ2 κ ⟨m, by omega⟩, sineMode_first, sineMode_mid]
    have hsin : 0 ≤ Real.sin ((((κ : ℕ) : ℝ) + 1) * (π / (2 * (m : ℝ) + 2))) := by
      apply Real.sin_nonneg_of_nonneg_of_le_pi (by positivity)
      have hk : ((κ : ℕ) : ℝ) + 1 ≤ 2 * (m : ℝ) + 2 := by
        have h1 : (κ : ℕ) ≤ 2 * m := by have := κ.isLt; omega
        have h2 : ((κ : ℕ) : ℝ) ≤ 2 * (m : ℝ) := by exact_mod_cast h1
        linarith
      calc (((κ : ℕ) : ℝ) + 1) * (π / (2 * (m : ℝ) + 2))
          ≤ (2 * (m : ℝ) + 2) * (π / (2 * (m : ℝ) + 2)) := by gcongr
        _ = π := by field_simp
    have hc : √(2 / (2 * (m : ℝ) + 2)) * √(2 / (2 * (m : ℝ) + 2)) = 2 / (2 * (m : ℝ) + 2) :=
      Real.mul_self_sqrt (by positivity)
    rw [abs_mul, abs_mul, abs_of_nonneg (Real.sqrt_nonneg _), abs_of_nonneg hsin,
      mul_mul_mul_comm, hc]
  rw [corrQ, ← hσinj.bijective_of_finite.sum_comp
    (fun β => |hA.eigenvectorBasis β ⟨0, by omega⟩| * |hA.eigenvectorBasis β ⟨m, by omega⟩|)]
  simp only [hterm]
  rw [← Finset.mul_sum, Fin.sum_univ_eq_sum_range (fun k : ℕ =>
      Real.sin (((k : ℝ) + 1) * (π / (2 * (m : ℝ) + 2))) * |Real.sin (((k : ℝ) + 1) * (π / 2))|)
      (2 * m + 1), sum_free_mid m, mul_one_div, div_div]

/-- `Q_{Δ_N}(1, (N+1)/2) → 2/π` along odd `N`. -/
theorem tendsto_corrQ_free_mid (J : ℝ) (hJ : J ≠ 0) :
    Tendsto (fun m : ℕ => corrQ (hq J 0 (0 : Cfg (2 * m + 1))) (hq_isHermitian J 0 0) ⟨0, by omega⟩
      ⟨m, by omega⟩) atTop (𝓝 (2 / π)) := by
  have e : ∀ m : ℕ, corrQ (hq J 0 (0 : Cfg (2 * m + 1))) (hq_isHermitian J 0 0) ⟨0, by omega⟩
      ⟨m, by omega⟩ = 2 / (π * Real.sinc (π / (2 * (m : ℝ) + 2))) := by
    intro m
    have hm2 : (2 * (m : ℝ) + 2) ≠ 0 := by positivity
    have hπ : π ≠ 0 := Real.pi_ne_zero
    rw [corrQ_free_mid J hJ m, Real.sinc_of_ne_zero (by positivity)]
    congr 1
    field_simp
  refine Tendsto.congr (fun m => (e m).symm) ?_
  have hu : Tendsto (fun m : ℕ => π / (2 * (m : ℝ) + 2)) atTop (𝓝 0) := by
    refine tendsto_const_nhds.div_atTop ?_
    refine tendsto_atTop_mono (fun m => ?_) tendsto_natCast_atTop_atTop
    have := (Nat.cast_nonneg m : (0 : ℝ) ≤ m)
    linarith
  have hs : Tendsto (fun m : ℕ => Real.sinc (π / (2 * (m : ℝ) + 2))) atTop (𝓝 1) := by
    have := (Real.continuous_sinc.tendsto 0).comp hu
    rwa [Real.sinc_zero] at this
  have hπs := hs.const_mul π
  rw [mul_one] at hπs
  exact tendsto_const_nhds.div hπs Real.pi_ne_zero

/-- **The constants of (EC) must degenerate at weak disorder** (paper: `rem:rate`(iii)): for `J ≠ 0` there are no
`C, η > 0` with `E_q Q^{[1,N]}_q(x, y) ≤ C e^{-η|x-y|}` for all `0 < λ = 2h/J ≤ λ₁`, all `N` and all `x, y`. -/
theorem no_uniform_EC_weak (J : ℝ) (hJ : J ≠ 0) (lam1 : ℝ) (hlam1 : 0 < lam1) (C η : ℝ)
    (hη : 0 < η) :
    ¬ ∀ h : ℝ, 0 < 2 * h / J → 2 * h / J ≤ lam1 → ∀ (N : ℕ) (x y : Fin N),
      avgQ (fun q : Cfg N => corrQ (hq J h q) (hq_isHermitian J h q) x y) ≤
        C * Real.exp (-η * (Nat.dist x.val y.val : ℝ)) := by
  intro H
  have hJp : 0 < |J| := abs_pos.mpr hJ
  -- a distance `m ≥ 1` at which `C e^{-ηm} < 1/4`
  have hlim : Tendsto (fun m : ℕ => C * Real.exp (-η * (m : ℝ))) atTop (𝓝 0) := by
    have h1 : Tendsto (fun m : ℕ => -η * (m : ℝ)) atTop atBot := by
      have := tendsto_neg_atTop_atBot.comp (tendsto_natCast_atTop_atTop.const_mul_atTop hη)
      refine this.congr fun m => ?_
      simp [Function.comp]
    have h2 := (Real.tendsto_exp_atBot.comp h1).const_mul C
    simpa using h2
  obtain ⟨m, hmC, hm1⟩ := ((hlim.eventually (gt_mem_nhds (by norm_num : (0 : ℝ) < 1 / 4))).and
    (eventually_ge_atTop 1)).exists
  have hmC' : C * Real.exp (-η * (m : ℝ)) < 1 / 4 := hmC
  have hm0 : (1 : ℝ) ≤ m := by exact_mod_cast hm1
  -- the disorder strength `λ = min(λ₁, 1/m)`, `h = Jλ/2`
  set lam : ℝ := min lam1 (1 / (m : ℝ)) with hlam
  have hlam0 : 0 < lam := lt_min hlam1 (by positivity)
  have hlamle : lam ≤ lam1 := min_le_left _ _
  have hlamm : (m : ℝ) * lam ≤ 1 := by
    have : lam ≤ 1 / m := min_le_right _ _
    calc (m : ℝ) * lam ≤ m * (1 / m) := by gcongr
      _ = 1 := by field_simp
  have hlam1' : lam ≤ 1 := by nlinarith
  set h : ℝ := J * lam / 2 with hh
  have h2h : 2 * h / J = lam := by rw [hh]; field_simp
  have habsh : |h| = |J| * lam / 2 := by
    rw [hh, abs_div, abs_mul, abs_of_pos hlam0, abs_two]
  -- `(|J|/(|J|+|h|))^m ≥ 1/2`
  have hratio : 1 / 2 ≤ (|J| / (|J| + |h|)) ^ m := by
    have h1 : 1 - lam / 2 ≤ |J| / (|J| + |h|) := by
      rw [habsh, le_div_iff₀ (by nlinarith)]
      nlinarith [mul_nonneg hJp.le (sq_nonneg lam)]
    have h2 : 1 + (m : ℝ) * (-(lam / 2)) ≤ (1 + -(lam / 2)) ^ m :=
      one_add_mul_le_pow (by linarith) m
    have h3 : (1 - lam / 2) ^ m ≤ (|J| / (|J| + |h|)) ^ m :=
      pow_le_pow_left₀ (by linarith) h1 m
    rw [← sub_eq_add_neg] at h2
    nlinarith
  have hlow : ∀ q : Cfg (m + 1), 1 / 4 ≤
      corrQ (hq J h q) (hq_isHermitian J h q) ⟨0, by omega⟩ ⟨m, by omega⟩ := by
    intro q
    have := WeakDisorderAux.corrQ_end_ge J h hJ m hm1 q
    linarith
  have havg : 1 / 4 ≤ avgQ (fun q : Cfg (m + 1) =>
      corrQ (hq J h q) (hq_isHermitian J h q) ⟨0, by omega⟩ ⟨m, by omega⟩) := by
    calc (1 / 4 : ℝ) = avgQ (fun _ : Cfg (m + 1) => (1 / 4 : ℝ)) :=
          (ImbalanceAux.avgQ_const _).symm
      _ ≤ _ := ImbalanceAux.avgQ_mono hlow
  have hH := H h (by rw [h2h]; exact hlam0) (by rw [h2h]; exact hlamle) (m + 1)
    ⟨0, by omega⟩ ⟨m, by omega⟩
  have hd : ((Nat.dist (⟨0, by omega⟩ : Fin (m + 1)).val (⟨m, by omega⟩ : Fin (m + 1)).val : ℕ) : ℝ)
      = m := by
    simp [Nat.dist_zero_left]
  rw [hd] at hH
  linarith

/-- **Chebyshev lower bound on the correlator** (paper: `rem:rate`(iii)): on `m + 1` sites and for every charge
configuration, `Q_q(1, m + 1) ≥ ½ (|J|/(|J| + |h|))^m = ½ (1 + |λ|/2)^{-m}`, with the test function
`T_m(E/(2|J| + 2|h|))` (Chebyshev polynomial, bounded by one on the spectrum). -/
theorem corrQ_end_lower (J h : ℝ) (hJ : J ≠ 0) (m : ℕ) (hm : 1 ≤ m) (q : Cfg (m + 1)) :
    1 / 2 * (|J| / (|J| + |h|)) ^ m ≤
      corrQ (hq J h q) (hq_isHermitian J h q) ⟨0, by omega⟩ ⟨m, by omega⟩ :=
  WeakDisorderAux.corrQ_end_ge J h hJ m hm q

/-- **An upper bound on the rate of (EC)** (paper: `rem:rate`(iii)): if `E_q Q^{[1,N]}_q(x, y) ≤ C e^{-η|x-y|}` for
all `N`, `x`, `y`, then `η ≤ ln(1 + |h|/|J|) = ln(1 + |λ|/2)`; in particular `η → 0` as `λ → 0`. -/
theorem ec_rate_le (J h : ℝ) (hJ : J ≠ 0) (C η : ℝ)
    (H : ∀ (N : ℕ) (x y : Fin N),
      avgQ (fun q : Cfg N => corrQ (hq J h q) (hq_isHermitian J h q) x y) ≤
        C * Real.exp (-η * (Nat.dist x.val y.val : ℝ))) :
    η ≤ Real.log (1 + |h| / |J|) := by
  have hJp : 0 < |J| := abs_pos.mpr hJ
  set L := Real.log (1 + |h| / |J|) with hL
  have ha : 0 < 1 + |h| / |J| := by positivity
  have hr : |J| / (|J| + |h|) = Real.exp (-L) := by
    rw [hL, Real.exp_neg, Real.exp_log ha]
    field_simp
  by_contra hη
  push Not at hη
  have key : ∀ m : ℕ, 1 ≤ m → Real.exp ((η - L) * m) ≤ 2 * C := by
    intro m hm
    have h1 : 1 / 2 * (|J| / (|J| + |h|)) ^ m ≤ avgQ (fun q : Cfg (m + 1) =>
        corrQ (hq J h q) (hq_isHermitian J h q) ⟨0, by omega⟩ ⟨m, by omega⟩) := by
      calc 1 / 2 * (|J| / (|J| + |h|)) ^ m =
          avgQ (fun _ : Cfg (m + 1) => 1 / 2 * (|J| / (|J| + |h|)) ^ m) :=
            (ImbalanceAux.avgQ_const _).symm
        _ ≤ _ := ImbalanceAux.avgQ_mono (fun q => corrQ_end_lower J h hJ m hm q)
    have h2 := H (m + 1) ⟨0, by omega⟩ ⟨m, by omega⟩
    have hd : ((Nat.dist (⟨0, by omega⟩ : Fin (m + 1)).val (⟨m, by omega⟩ : Fin (m + 1)).val : ℕ) :
        ℝ) = m := by
      simp [Nat.dist_zero_left]
    rw [hd] at h2
    rw [hr, ← Real.exp_nat_mul] at h1
    have h3 : 1 / 2 * Real.exp (↑m * -L) ≤ C * Real.exp (-η * ↑m) := h1.trans h2
    have e : Real.exp ((η - L) * m) = Real.exp (↑m * -L) * Real.exp (η * m) := by
      rw [← Real.exp_add]; congr 1; ring
    have e2 : Real.exp (-η * ↑m) * Real.exp (η * m) = 1 := by
      rw [← Real.exp_add]; simp
    have hpos : 0 < Real.exp (η * m) := Real.exp_pos _
    calc Real.exp ((η - L) * m) = Real.exp (↑m * -L) * Real.exp (η * m) := e
      _ = 2 * (1 / 2 * Real.exp (↑m * -L)) * Real.exp (η * m) := by ring
      _ ≤ 2 * (C * Real.exp (-η * ↑m)) * Real.exp (η * m) := by gcongr
      _ = 2 * C * (Real.exp (-η * ↑m) * Real.exp (η * m)) := by ring
      _ = 2 * C := by rw [e2, mul_one]
  have hδ : 0 < η - L := by linarith
  obtain ⟨m, hm⟩ := exists_nat_gt (2 * C / (η - L) + 1)
  have hm1 : 1 ≤ m := by
    have : (0 : ℝ) < m := by
      have hC : 0 < 2 * C := by
        have := key 1 le_rfl
        exact lt_of_lt_of_le (Real.exp_pos _) this
      have : 0 ≤ 2 * C / (η - L) := div_nonneg hC.le hδ.le
      linarith
    exact_mod_cast this
  have hk := key m hm1
  have hexp := Real.add_one_le_exp ((η - L) * m)
  have : (η - L) * m > 2 * C := by
    have h1 : 2 * C / (η - L) < m := by linarith
    rw [div_lt_iff₀ hδ] at h1
    linarith
  linarith

end DFL
