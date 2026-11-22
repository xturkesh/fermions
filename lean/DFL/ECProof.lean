import DFL.LDTEntry
import DFL.Main
import DFL.AreaLaw

/-!
# Hypothesis (EC) from uniform large deviations (paper: Theorem `thm:EC`(i))

The paper's proof of Hypothesis (EC), uniformly in the volume and up to the boundary: for `x < y`,
`n = ⌊(d-1)/2⌋`, boxes of radius `n` around `x` and `y` truncated by `[0, N-1]`, and the event that no
energy is resonant for both boxes (probability `≥ 1 - 4(8n+6) C₁ e^{-c₁ n}`, `two_box`), every
eigenvector satisfies the Poisson bound at `x` or at `y` (`decay_of_not_res`), which gives
`Q(x, y) ≤ 4 e^{-κ n}` on that event.

* `corrQ_decay_of_ldt`: `E_q Q_q(x, y) ≤ C e^{-η |x-y|}` for `h_q = -J Δ + 2h diag(q)`, uniformly in
  `N`, given `LDT (-(2h/J))`.
* `ec_of_ldt`: Hypothesis (EC) holds with `F(d) = min(1, C e^{-η d})`.
* `lgt_transport_of_ldt`, `lgt_numberEntropy_of_ldt`, `lgt_areaLaw_of_ldt`: the main results of the
  paper with Hypothesis (EC) replaced by the literature input `LDT`.
-/

open Matrix Finset

noncomputable section

namespace DFL

namespace ECProofAux

variable {N : ℕ}

/-! ### Elementary facts on the correlator and on averages -/

/-- `|a| |b| ≤ (a² + b²) / 2`. -/
lemma abs_mul_abs_le (a b : ℝ) : |a| * |b| ≤ (a ^ 2 + b ^ 2) / 2 := by
  nlinarith [sq_nonneg (|a| - |b|), sq_abs a, sq_abs b]

/-- Cauchy–Schwarz for two rows of the eigenvector matrix: `∑_α |φ_α(u)| |φ_α(v)| ≤ 1`. -/
lemma sum_abs_eigvec_mul_le_one (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian)
    (u v : Fin N) :
    ∑ α, |hA.eigenvectorBasis α u| * |hA.eigenvectorBasis α v| ≤ 1 := by
  calc ∑ α, |hA.eigenvectorBasis α u| * |hA.eigenvectorBasis α v|
      ≤ ∑ α, ((hA.eigenvectorBasis α u) ^ 2 + (hA.eigenvectorBasis α v) ^ 2) / 2 :=
        Finset.sum_le_sum fun α _ => abs_mul_abs_le _ _
    _ = ((∑ α, (hA.eigenvectorBasis α u) ^ 2) + ∑ α, (hA.eigenvectorBasis α v) ^ 2) / 2 := by
        rw [← Finset.sum_div, Finset.sum_add_distrib]
    _ = 1 := by rw [OneBodyAux.sum_eigvec_sq, OneBodyAux.sum_eigvec_sq]; norm_num

/-- The eigenfunction correlator is at most one. -/
lemma corrQ_le_one (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x y : Fin N) :
    corrQ A hA x y ≤ 1 :=
  sum_abs_eigvec_mul_le_one A hA x y

/-- The eigenfunction correlator is symmetric. -/
lemma corrQ_symm (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x y : Fin N) :
    corrQ A hA x y = corrQ A hA y x := by
  unfold corrQ
  exact Finset.sum_congr rfl fun α _ => mul_comm _ _

/-- `extZ` at a site of the chain. -/
lemma extZ_coe (φ : Fin N → ℝ) (i : Fin N) : extZ φ (i.val : ℤ) = φ i := by
  have h : (0 : ℤ) ≤ (i.val : ℤ) ∧ (i.val : ℤ) < N :=
    ⟨Int.natCast_nonneg _, by exact_mod_cast i.isLt⟩
  rw [extZ, dite_eq_left h]
  congr 1

/-- `extZ` vanishes outside the chain. -/
lemma extZ_of_not (φ : Fin N → ℝ) (z : ℤ) (hz : ¬ (0 ≤ z ∧ z < N)) : extZ φ z = 0 := by
  rw [extZ, dite_eq_right hz]

/-- Cauchy–Schwarz with one row evaluated at an arbitrary integer (zero outside the chain). -/
lemma sum_abs_extZ_mul_le_one (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (u : ℤ)
    (v : Fin N) :
    ∑ α, |extZ (fun j => hA.eigenvectorBasis α j) u| * |hA.eigenvectorBasis α v| ≤ 1 := by
  by_cases hu : 0 ≤ u ∧ u < N
  · have hk : u.toNat < N := by omega
    have : ∀ α, extZ (fun j => hA.eigenvectorBasis α j) u =
        hA.eigenvectorBasis α ⟨u.toNat, hk⟩ := by
      intro α
      rw [extZ, dite_eq_left hu]
    simp only [this]
    exact sum_abs_eigvec_mul_le_one A hA _ v
  · simp [extZ_of_not _ _ hu]

/-- Cauchy–Schwarz with one row evaluated at an arbitrary integer, other order. -/
lemma sum_abs_mul_extZ_le_one (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (u : ℤ)
    (v : Fin N) :
    ∑ α, |hA.eigenvectorBasis α v| * |extZ (fun j => hA.eigenvectorBasis α j) u| ≤ 1 := by
  simp only [mul_comm |hA.eigenvectorBasis _ v|]
  exact sum_abs_extZ_mul_le_one A hA u v

/-- There are `2^N` charge configurations. -/
lemma card_cfg : (Finset.univ : Finset (Cfg N)).card = 2 ^ N := by
  simp [Finset.card_univ, ZMod.card]

/-- An average of a function bounded by `a` is bounded by `a`. -/
lemma avgQ_le_of_le (f : Cfg N → ℝ) (a : ℝ) (hf : ∀ q, f q ≤ a) : avgQ f ≤ a := by
  unfold avgQ
  have h2 : (0 : ℝ) < 2 ^ N := by positivity
  have hs : ∑ q, f q ≤ 2 ^ N * a := by
    calc ∑ q, f q ≤ ∑ _q : Cfg N, a := Finset.sum_le_sum fun q _ => hf q
      _ = 2 ^ N * a := by rw [Finset.sum_const, card_cfg, nsmul_eq_mul]; push_cast; ring
  rw [div_mul_eq_mul_div, one_mul, div_le_iff₀ h2]
  linarith

/-- Average of a function bounded by `1`, and by `a` outside a set of probability `≤ B`. -/
lemma avgQ_le_add (f : Cfg N → ℝ) (Bad : Finset (Cfg N)) (a B : ℝ) (hf1 : ∀ q, f q ≤ 1)
    (hfa : ∀ q, q ∉ Bad → f q ≤ a) (ha : 0 ≤ a) (hB : (Bad.card : ℝ) ≤ B * 2 ^ N) :
    avgQ f ≤ a + B := by
  have hpt : ∀ q, f q ≤ a + (if q ∈ Bad then 1 else 0) := by
    intro q
    by_cases hq : q ∈ Bad
    · rw [ite_eq_left hq]
      linarith [hf1 q]
    · rw [ite_eq_right hq, add_zero]
      exact hfa q hq
  have hs : ∑ q, f q ≤ 2 ^ N * a + Bad.card := by
    calc ∑ q, f q ≤ ∑ q, (a + (if q ∈ Bad then 1 else 0)) := Finset.sum_le_sum fun q _ => hpt q
      _ = 2 ^ N * a + Bad.card := by
        rw [Finset.sum_add_distrib, Finset.sum_const, card_cfg, nsmul_eq_mul, Finset.sum_boole,
          Finset.filter_mem_eq_inter, Finset.univ_inter]
        push_cast
        ring
  unfold avgQ
  have h2 : (0 : ℝ) < 2 ^ N := by positivity
  rw [div_mul_eq_mul_div, one_mul, div_le_iff₀ h2]
  linarith

/-! ### The deterministic step -/

/-- The rescaled eigenvalue `E = -E'/J` of `h_q` lies in `I_λ`, `λ = -(2h/J)`. -/
lemma eigval_mem_Ilam (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) (α : Fin N) :
    -(hq_isHermitian J h q).eigenvalues α / J ∈ Ilam (-(2 * h / J)) := by
  have hb := abs_eigenvalue_hq_le J h q α
  have hJp : 0 < |J| := abs_pos.mpr hJ
  have hlam : |-(2 * h / J)| * |J| = 2 * |h| := by
    rw [abs_neg, abs_div, abs_mul, div_mul_cancel₀ _ hJp.ne']
    norm_num
  have habs : |-(hq_isHermitian J h q).eigenvalues α / J| ≤ 2 + |-(2 * h / J)| := by
    rw [abs_div, abs_neg, div_le_iff₀ hJp, add_mul, hlam]
    linarith
  exact ⟨by linarith [(abs_le.mp habs).1], (abs_le.mp habs).2⟩

/-- **Decay at one point** (paper: deterministic step): if the rescaled eigenvalue is not resonant
for a box `[a, a+L-1] ⊆ [0, N-1]` around `x = a + i`, and the unselected ends of the box are ends
of the chain, then `|φ_α(x)| ≤ θ (|φ_α(a-1)| + |φ_α(a+L)|)`. -/
lemma eigvec_decay (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) (α x : Fin N) (θ : ℝ) (a : ℤ) (L i : ℕ)
    (endL endR : Bool) (ha : 0 ≤ a) (haL : a + L ≤ N) (hi : i < L) (hx : a + i = x.val)
    (hend : endL = true ∨ endR = true) (hL : endL = false → a - 1 < 0)
    (hR : endR = false → (N : ℤ) ≤ a + L)
    (hres : -(hq_isHermitian J h q).eigenvalues α / J ∉
      resSet (potZ (-(2 * h / J)) q) (-(2 * h / J)) a L i endL endR θ) :
    |(hq_isHermitian J h q).eigenvectorBasis α x| ≤
      θ * (|extZ (fun j => (hq_isHermitian J h q).eigenvectorBasis α j) (a - 1)| +
        |extZ (fun j => (hq_isHermitian J h q).eigenvectorBasis α j) (a + L)|) := by
  have key := decay_of_not_res (potZ (-(2 * h / J)) q) (-(2 * h / J)) a L i endL endR θ hi hend
    _ (eigval_mem_Ilam J h hJ q α) hres
    (extZ (fun j => (hq_isHermitian J h q).eigenvectorBasis α j))
    (fun z hz1 hz2 => eigvec_recursion J h hJ q α z (by omega) (by omega))
    (fun hf => extZ_of_not _ _ (by have := hL hf; omega))
    (fun hf => extZ_of_not _ _ (by have := hR hf; omega))
  rw [hx, extZ_coe] at key
  exact key

/-- **Deterministic step** (paper: Eq. `eq:Qgood`): if no energy is resonant for both boxes, then
`Q(x, y) ≤ 4θ`. -/
lemma corrQ_le_of_good (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) (x y : Fin N) (θ : ℝ) (hθ : 0 ≤ θ)
    (a1 : ℤ) (L1 i1 : ℕ) (endL1 endR1 : Bool) (a2 : ℤ) (L2 i2 : ℕ) (endL2 endR2 : Bool)
    (ha1 : 0 ≤ a1) (ha1L : a1 + L1 ≤ N) (hi1 : i1 < L1) (hx : a1 + i1 = x.val)
    (hend1 : endL1 = true ∨ endR1 = true) (hL1 : endL1 = false → a1 - 1 < 0)
    (hR1 : endR1 = false → (N : ℤ) ≤ a1 + L1)
    (ha2 : 0 ≤ a2) (ha2L : a2 + L2 ≤ N) (hi2 : i2 < L2) (hy : a2 + i2 = y.val)
    (hend2 : endL2 = true ∨ endR2 = true) (hL2 : endL2 = false → a2 - 1 < 0)
    (hR2 : endR2 = false → (N : ℤ) ≤ a2 + L2)
    (hgood : ¬ (resSet (potZ (-(2 * h / J)) q) (-(2 * h / J)) a1 L1 i1 endL1 endR1 θ ∩
        resSet (potZ (-(2 * h / J)) q) (-(2 * h / J)) a2 L2 i2 endL2 endR2 θ).Nonempty) :
    corrQ (hq J h q) (hq_isHermitian J h q) x y ≤ 4 * θ := by
  set hA := hq_isHermitian J h q with hA_def
  set φ : Fin N → ℤ → ℝ := fun α => extZ (fun j => hA.eigenvectorBasis α j) with hφ_def
  have hpt : ∀ α, |hA.eigenvectorBasis α x| * |hA.eigenvectorBasis α y| ≤
      θ * (|φ α (a1 - 1)| * |hA.eigenvectorBasis α y| +
          |φ α (a1 + L1)| * |hA.eigenvectorBasis α y|) +
        θ * (|hA.eigenvectorBasis α x| * |φ α (a2 - 1)| +
          |hA.eigenvectorBasis α x| * |φ α (a2 + L2)|) := by
    intro α
    have hx0 := abs_nonneg (hA.eigenvectorBasis α x)
    have hy0 := abs_nonneg (hA.eigenvectorBasis α y)
    have h1n := abs_nonneg (φ α (a1 - 1))
    have h2n := abs_nonneg (φ α (a1 + L1))
    have h3n := abs_nonneg (φ α (a2 - 1))
    have h4n := abs_nonneg (φ α (a2 + L2))
    by_cases h1 : -hA.eigenvalues α / J ∈
        resSet (potZ (-(2 * h / J)) q) (-(2 * h / J)) a1 L1 i1 endL1 endR1 θ
    · have h2 : -hA.eigenvalues α / J ∉
          resSet (potZ (-(2 * h / J)) q) (-(2 * h / J)) a2 L2 i2 endL2 endR2 θ :=
        fun h2 => hgood ⟨_, h1, h2⟩
      have hd := eigvec_decay J h hJ q α y θ a2 L2 i2 endL2 endR2 ha2 ha2L hi2 hy hend2 hL2 hR2 h2
      have hd' := mul_le_mul_of_nonneg_left hd hx0
      have : 0 ≤ θ * (|φ α (a1 - 1)| * |hA.eigenvectorBasis α y| +
          |φ α (a1 + L1)| * |hA.eigenvectorBasis α y|) := by positivity
      have e : |hA.eigenvectorBasis α x| * (θ * (|φ α (a2 - 1)| + |φ α (a2 + L2)|)) =
          θ * (|hA.eigenvectorBasis α x| * |φ α (a2 - 1)| +
            |hA.eigenvectorBasis α x| * |φ α (a2 + L2)|) := by ring
      linarith
    · have hd := eigvec_decay J h hJ q α x θ a1 L1 i1 endL1 endR1 ha1 ha1L hi1 hx hend1 hL1 hR1 h1
      have hd' := mul_le_mul_of_nonneg_right hd hy0
      have : 0 ≤ θ * (|hA.eigenvectorBasis α x| * |φ α (a2 - 1)| +
          |hA.eigenvectorBasis α x| * |φ α (a2 + L2)|) := by positivity
      have e : θ * (|φ α (a1 - 1)| + |φ α (a1 + L1)|) * |hA.eigenvectorBasis α y| =
          θ * (|φ α (a1 - 1)| * |hA.eigenvectorBasis α y| +
            |φ α (a1 + L1)| * |hA.eigenvectorBasis α y|) := by ring
      linarith
  have s1 := sum_abs_extZ_mul_le_one (hq J h q) hA (a1 - 1) y
  have s2 := sum_abs_extZ_mul_le_one (hq J h q) hA (a1 + L1) y
  have s3 := sum_abs_mul_extZ_le_one (hq J h q) hA (a2 - 1) x
  have s4 := sum_abs_mul_extZ_le_one (hq J h q) hA (a2 + L2) x
  calc corrQ (hq J h q) hA x y
      = ∑ α, |hA.eigenvectorBasis α x| * |hA.eigenvectorBasis α y| := rfl
    _ ≤ ∑ α, (θ * (|φ α (a1 - 1)| * |hA.eigenvectorBasis α y| +
          |φ α (a1 + L1)| * |hA.eigenvectorBasis α y|) +
        θ * (|hA.eigenvectorBasis α x| * |φ α (a2 - 1)| +
          |hA.eigenvectorBasis α x| * |φ α (a2 + L2)|)) := Finset.sum_le_sum fun α _ => hpt α
    _ = θ * ((∑ α, |φ α (a1 - 1)| * |hA.eigenvectorBasis α y|) +
          ∑ α, |φ α (a1 + L1)| * |hA.eigenvectorBasis α y|) +
        θ * ((∑ α, |hA.eigenvectorBasis α x| * |φ α (a2 - 1)|) +
          ∑ α, |hA.eigenvectorBasis α x| * |φ α (a2 + L2)|) := by
        rw [Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum, Finset.sum_add_distrib,
          Finset.sum_add_distrib]
    _ ≤ θ * (1 + 1) + θ * (1 + 1) := by
        gcongr
    _ = 4 * θ := by ring

/-! ### The probabilistic step -/

/-- **Probabilistic step** (paper: Eq. `eq:ECn`): for `x < y`, `n = ⌊(d-1)/2⌋ ≥ m₀`,
`E_q Q_q(x, y) ≤ 4 e^{-κn} + 4(8n+6) C₁ e^{-c₁n}`. -/
lemma avgQ_le_of_entry (J h : ℝ) (hJ : J ≠ 0) (κ C1 c1 : ℝ) (m0 : ℕ)
    (hEB : EntryBound (-(2 * h / J)) κ C1 c1 m0) (hC1 : 0 ≤ C1) (x y : Fin N)
    (hxy : x.val < y.val) (n : ℕ) (hn_def : n = (y.val - x.val - 1) / 2) (hn : m0 ≤ n) :
    avgQ (fun q => corrQ (hq J h q) (hq_isHermitian J h q) x y) ≤
      4 * Real.exp (-(κ * n)) + 4 * (8 * n + 6) * C1 * Real.exp (-c1 * n) := by
  have h2n : x.val + 2 * n + 1 ≤ y.val := by omega
  have hyN := y.isLt
  obtain ⟨a1, ha1⟩ : ∃ a1 : ℕ, a1 = x.val - n := ⟨_, rfl⟩
  obtain ⟨L1, hL1⟩ : ∃ L1 : ℕ, L1 = x.val + n + 1 - a1 := ⟨_, rfl⟩
  obtain ⟨i1, hi1⟩ : ∃ i1 : ℕ, i1 = x.val - a1 := ⟨_, rfl⟩
  obtain ⟨a2, ha2⟩ : ∃ a2 : ℕ, a2 = y.val - n := ⟨_, rfl⟩
  obtain ⟨L2, hL2⟩ : ∃ L2 : ℕ, L2 = n + 1 + min n (N - 1 - y.val) := ⟨_, rfl⟩
  have hB := two_box (N := N) (-(2 * h / J)) κ C1 c1 m0 hEB hC1 (a1 : ℤ) L1 i1
    (decide (n + 1 ≤ x.val)) true (a2 : ℤ) L2 n true (decide (y.val + n + 2 ≤ N)) n
    (by positivity) (by omega) (by omega)
    (by omega) (by omega) (by omega) (by omega) hn (fun _ => by omega)
    (fun hd => by simp only [decide_eq_true_eq] at hd; omega)
    (fun hd => by simp only [decide_eq_true_eq] at hd; omega) (fun _ => rfl)
  refine avgQ_le_add _ _ _ _ (fun q => corrQ_le_one _ _ x y) (fun q hq => ?_)
    (by positivity) hB
  have hgood : ¬ (resSet (potZ (-(2 * h / J)) q) (-(2 * h / J)) a1 L1 i1
        (decide (n + 1 ≤ x.val)) true (Real.exp (-(κ * n))) ∩
      resSet (potZ (-(2 * h / J)) q) (-(2 * h / J)) a2 L2 n true
        (decide (y.val + n + 2 ≤ N)) (Real.exp (-(κ * n)))).Nonempty := by
    simpa using hq
  exact corrQ_le_of_good J h hJ q x y _ (Real.exp_pos _).le (a1 : ℤ) L1 i1 _ true (a2 : ℤ) L2 n
    true _ (by positivity) (by omega) (by omega) (by omega) (Or.inr rfl)
    (fun hd => by simp only [decide_eq_false_iff_not] at hd; omega) (fun hd => by simp at hd)
    (by positivity) (by omega) (by omega) (by omega) (Or.inl rfl) (fun hd => by simp at hd)
    (fun hd => by simp only [decide_eq_false_iff_not] at hd; omega) hgood

/-! ### Constants -/

/-- `(8n+6) e^{-c₁ n} ≤ (6 + 16/c₁) e^{-cn}` for `0 < c ≤ c₁/2`. -/
lemma poly_exp_le (c1 c : ℝ) (hc1 : 0 < c1) (hc : c ≤ c1 / 2) (n : ℕ) :
    (8 * (n : ℝ) + 6) * Real.exp (-c1 * n) ≤ (6 + 16 / c1) * Real.exp (-(c * n)) := by
  have hn0 : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  have hK0 : 0 ≤ 6 + 16 / c1 := by positivity
  have hK : 8 * (n : ℝ) + 6 ≤ (6 + 16 / c1) * (1 + c1 * n / 2) := by
    have e : (6 + 16 / c1) * (1 + c1 * n / 2) = 6 + 3 * c1 * n + 16 / c1 + 8 * n := by
      field_simp
      ring
    rw [e]
    have : 0 ≤ 16 / c1 := by positivity
    nlinarith
  have hexp : 1 + c1 * n / 2 ≤ Real.exp ((c1 - c) * n) := by
    have h1 := Real.add_one_le_exp (c1 * n / 2)
    have h2 : Real.exp (c1 * n / 2) ≤ Real.exp ((c1 - c) * n) := by
      apply Real.exp_le_exp.mpr
      nlinarith
    linarith
  calc (8 * (n : ℝ) + 6) * Real.exp (-c1 * n)
      ≤ (6 + 16 / c1) * Real.exp ((c1 - c) * n) * Real.exp (-c1 * n) := by
        gcongr
        exact hK.trans (mul_le_mul_of_nonneg_left hexp hK0)
    _ = (6 + 16 / c1) * Real.exp (-(c * n)) := by
        rw [mul_assoc, ← Real.exp_add]
        congr 2
        ring

/-- `∑_d (d+1)² e^{-rd}` converges for `r > 0`. -/
lemma summable_sq_mul_exp (r : ℝ) (hr : 0 < r) :
    Summable (fun d : ℕ => ((d : ℝ) + 1) ^ 2 * Real.exp (-r * d)) := by
  have h2 := Real.summable_pow_mul_exp_neg_nat_mul 2 hr
  have h1 := Real.summable_pow_mul_exp_neg_nat_mul 1 hr
  have h0 := Real.summable_pow_mul_exp_neg_nat_mul 0 hr
  refine ((h2.add (h1.mul_left 2)).add h0).congr fun d => ?_
  ring

end ECProofAux

/-- **Uniform decay of the averaged eigenfunction correlator** (paper: Theorem `thm:EC`(i)). -/
theorem corrQ_decay_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ ∀ (N : ℕ) (x y : Fin N),
      avgQ (fun q => corrQ (hq J h q) (hq_isHermitian J h q) x y) ≤
        C * Real.exp (-η * Nat.dist x.val y.val) := by
  obtain ⟨κ, C1, c1, m0, hκ, hc1, hC1, hEB⟩ := entry_ldt _ hL
  set c := min κ (c1 / 2) with hc_def
  have hc : 0 < c := lt_min hκ (by linarith)
  have hcκ : c ≤ κ := min_le_left _ _
  have hcc1 : c ≤ c1 / 2 := min_le_right _ _
  set C' := 4 + 4 * C1 * (6 + 16 / c1) with hC'_def
  have hC' : 0 ≤ C' := by positivity
  set η := c / 2 with hη_def
  have hη : 0 < η := by positivity
  refine ⟨max (C' * Real.exp c) (Real.exp (η * (2 * m0 + 1))), η, hη, ?_⟩
  -- the case `x < y` with `n ≥ m₀`
  have key : ∀ (N : ℕ) (x y : Fin N), x.val < y.val → 2 * m0 < y.val - x.val →
      avgQ (fun q => corrQ (hq J h q) (hq_isHermitian J h q) x y) ≤
        C' * Real.exp c * Real.exp (-η * ((y.val - x.val : ℕ) : ℝ)) := by
    intro N x y hxy hd
    obtain ⟨n, hn_def⟩ : ∃ n : ℕ, n = (y.val - x.val - 1) / 2 := ⟨_, rfl⟩
    have hn : m0 ≤ n := by omega
    refine (ECProofAux.avgQ_le_of_entry J h hJ κ C1 c1 m0 hEB hC1 x y hxy n hn_def hn).trans ?_
    have hn0 : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    have hnR : ((y.val - x.val : ℕ) : ℝ) ≤ 2 * (n : ℝ) + 2 := by
      have : y.val - x.val ≤ 2 * n + 2 := by omega
      exact_mod_cast this
    have e1 : Real.exp (-(κ * n)) ≤ Real.exp (-(c * n)) := by
      apply Real.exp_le_exp.mpr
      nlinarith
    have e2 := ECProofAux.poly_exp_le c1 c hc1 hcc1 n
    have e3 : Real.exp (-(c * n)) ≤ Real.exp c * Real.exp (-η * ((y.val - x.val : ℕ) : ℝ)) := by
      rw [← Real.exp_add]
      apply Real.exp_le_exp.mpr
      rw [hη_def]
      nlinarith
    calc 4 * Real.exp (-(κ * n)) + 4 * (8 * n + 6) * C1 * Real.exp (-c1 * n)
        ≤ 4 * Real.exp (-(c * n)) + 4 * C1 * ((6 + 16 / c1) * Real.exp (-(c * n))) := by
          have := mul_le_mul_of_nonneg_left e2 (by positivity : (0 : ℝ) ≤ 4 * C1)
          nlinarith
      _ = C' * Real.exp (-(c * n)) := by rw [hC'_def]; ring
      _ ≤ C' * (Real.exp c * Real.exp (-η * ((y.val - x.val : ℕ) : ℝ))) := by gcongr
      _ = C' * Real.exp c * Real.exp (-η * ((y.val - x.val : ℕ) : ℝ)) := by ring
  -- short distances
  have small : ∀ (N : ℕ) (x y : Fin N), Nat.dist x.val y.val ≤ 2 * m0 →
      avgQ (fun q => corrQ (hq J h q) (hq_isHermitian J h q) x y) ≤
        Real.exp (η * (2 * m0 + 1)) * Real.exp (-η * Nat.dist x.val y.val) := by
    intro N x y hd
    refine (ECProofAux.avgQ_le_of_le _ 1 fun q => ECProofAux.corrQ_le_one _ _ x y).trans ?_
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
      have hsym : avgQ (fun q => corrQ (hq J h q) (hq_isHermitian J h q) x y) =
          avgQ (fun q => corrQ (hq J h q) (hq_isHermitian J h q) y x) := by
        unfold avgQ
        congr 1
        exact Finset.sum_congr rfl fun q _ => ECProofAux.corrQ_symm _ _ _ _
      rw [hsym]
      exact (key N y x hxy hd).trans
        (mul_le_mul_of_nonneg_right (le_max_left _ _) (Real.exp_pos _).le)

/-- **Hypothesis (EC) holds** given the uniform large deviations (paper: Theorem `thm:EC`). -/
theorem ec_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ EC J h (fun d => min 1 (C * Real.exp (-η * d))) := by
  obtain ⟨C, η, hη, hb⟩ := corrQ_decay_of_ldt J h hJ hL
  set C₀ := max C 1 with hC₀_def
  have hC₀ : 0 < C₀ := lt_of_lt_of_le one_pos (le_max_right _ _)
  refine ⟨C₀, η, hη, ?_⟩
  have hmaj : ∀ d : ℕ, min 1 (C₀ * Real.exp (-η * d)) ≤ C₀ * Real.exp (-η * d) :=
    fun d => min_le_right _ _
  have hnn : ∀ d : ℕ, 0 ≤ min 1 (C₀ * Real.exp (-η * d)) :=
    fun d => le_min zero_le_one (by positivity)
  refine ⟨hnn, ?_, ?_, ?_, ?_⟩
  · intro a b hab
    show min 1 (C₀ * Real.exp (-η * b)) ≤ min 1 (C₀ * Real.exp (-η * a))
    apply min_le_min le_rfl
    apply mul_le_mul_of_nonneg_left _ hC₀.le
    apply Real.exp_le_exp.mpr
    have : (a : ℝ) ≤ b := by exact_mod_cast hab
    nlinarith
  · intro N x y
    show _ ≤ min 1 (C₀ * Real.exp (-η * Nat.dist x.val y.val))
    apply le_min
    · exact ECProofAux.avgQ_le_of_le _ 1 fun q => ECProofAux.corrQ_le_one _ _ x y
    · exact (hb N x y).trans
        (mul_le_mul_of_nonneg_right (le_max_left _ _) (Real.exp_pos _).le)
  · refine Summable.of_nonneg_of_le (fun d => mul_nonneg (by positivity) (hnn d)) (fun d => ?_)
      ((ECProofAux.summable_sq_mul_exp η hη).mul_left C₀)
    calc ((d : ℝ) + 1) ^ 2 * min 1 (C₀ * Real.exp (-η * d))
        ≤ ((d : ℝ) + 1) ^ 2 * (C₀ * Real.exp (-η * d)) :=
          mul_le_mul_of_nonneg_left (hmaj d) (by positivity)
      _ = C₀ * (((d : ℝ) + 1) ^ 2 * Real.exp (-η * d)) := by ring
  · refine Summable.of_nonneg_of_le
      (fun d => mul_nonneg (by positivity) (Real.sqrt_nonneg _)) (fun d => ?_)
      ((ECProofAux.summable_sq_mul_exp (η / 2) (by positivity)).mul_left (Real.sqrt C₀))
    have hs : Real.sqrt (min 1 (C₀ * Real.exp (-η * d))) ≤
        Real.sqrt C₀ * Real.exp (-(η / 2) * d) := by
      calc Real.sqrt (min 1 (C₀ * Real.exp (-η * d)))
          ≤ Real.sqrt (C₀ * Real.exp (-η * d)) := Real.sqrt_le_sqrt (hmaj d)
        _ = Real.sqrt C₀ * Real.sqrt (Real.exp (-η * d)) := Real.sqrt_mul hC₀.le _
        _ = Real.sqrt C₀ * Real.exp (-(η / 2) * d) := by
            rw [← Real.exp_half]
            congr 2
            ring
    calc ((d : ℝ) + 1) ^ 2 * Real.sqrt (min 1 (C₀ * Real.exp (-η * d)))
        ≤ ((d : ℝ) + 1) ^ 2 * (Real.sqrt C₀ * Real.exp (-(η / 2) * d)) :=
          mul_le_mul_of_nonneg_left hs (by positivity)
      _ = Real.sqrt C₀ * (((d : ℝ) + 1) ^ 2 * Real.exp (-(η / 2) * d)) := by ring

/-! ### The main results with (EC) replaced by the literature input -/

variable {N : ℕ}

/-- **Bounded transport** (paper: Theorem `thm:transport`) from uniform large deviations. -/
theorem lgt_transport_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ K : ℝ, ∀ (N : ℕ) (S : Finset (Fin N)) (m : ℕ) (t : ℝ),
      |(expect (matter (numOp (leftOf m))) (lgtState J h (fockVec S) t)).re -
        ((S ∩ leftOf m).card : ℝ)| ≤ K := by
  obtain ⟨C, η, -, hEC⟩ := ec_of_ldt J h hJ hL
  exact ⟨_, fun N S m t => lgt_transport J h hEC S m t⟩

/-- **Bounded number entropy** (paper: Theorem `thm:numberentropy`) from uniform large deviations. -/
theorem lgt_numberEntropy_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ K : ℝ, ∀ (N : ℕ) (m : ℕ) (S : Finset (Fin N)) (t : ℝ),
      shannon N (lgtNumDist J h (leftOf m) (fockVec S) t) ≤
        (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (K + 1 / 12)) := by
  obtain ⟨C, η, -, hEC⟩ := ec_of_ldt J h hJ hL
  exact lgt_numberEntropy J h hEC

/-- **Area law** (paper: Theorem `thm:arealaw`(a),(b)) from uniform large deviations. -/
theorem lgt_areaLaw_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ K : ℝ, ∀ (N : ℕ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (S : Finset (Fin N)) (t : ℝ),
      avgQ (fun q : Cfg N =>
        vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ fockVec S))
          (rdm_isHermitian ℓ hℓ _)) ≤ K := by
  obtain ⟨C, η, -, hEC⟩ := ec_of_ldt J h hJ hL
  exact lgt_areaLaw J h hEC

/-- **Persistent imbalance** (paper: Theorem `thm:imbalance`) from uniform large deviations. -/
theorem lgtImbalance_lower_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hh : h ≠ 0)
    (hL : LDT (-(2 * h / J))) :
    ∃ c : ℝ, 0 < c ∧ ∀ (N : ℕ), Even N → 2 ≤ N →
      ∀ ε > 0, ∀ᶠ T in Filter.atTop,
        c - ε ≤ (1 / T) * ∫ t in (0 : ℝ)..T, lgtImbalance J h N t := by
  obtain ⟨C, η, -, hEC⟩ := ec_of_ldt J h hJ hL
  exact lgtImbalance_lower J h hJ hh hEC

end DFL
