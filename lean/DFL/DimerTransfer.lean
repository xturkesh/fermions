import DFL.DimerCritical

/-!
# Dimer transfer matrices: rotation normal form and the uniform bound
(paper: Lemma `lem:dimercritical`(i),(iii),(iv), eq. `eq:dimerangles`, Lemma `lem:dimertransfer`)

For `a ∈ ℂ`, `T_a = !![a, -1; 1, 0]`; the polymer transfer matrices of the random dimer model at energy `z`
are `D_±(z) = T_{±λ - z}²`.

* (i) `[D_+(E), D_-(E)] = 2λ(λ² - E²) σ_x`, so `D_±(E)` commute iff `E = ±λ` (`dimerD_commutator`,
  `dimerD_commute_iff`).
* (iii) For `0 < λ < 1`, `θ = arccos(-λ) ∈ (π/2, π)` and `S = !![sin θ, cos θ; 0, 1]`:
  `S⁻¹ D_+(λ) S = -1 = R_π`, `S⁻¹ T_{-2λ} S = R_θ`, `S⁻¹ D_-(λ) S = R_{2θ}`; the averaged phases
  `⟨e^{2iη±}⟩ = ½(1 + e^{4iθ})` have modulus `|1 - 2λ²| < 1`, and `⟨e^{4iη±}⟩ = 1 ⟺ λ = 1/√2`;
  `κ = ‖S‖ ‖S⁻¹‖ = ((1+λ)/(1-λ))^{1/2}`.
* (iv) `Γ₂ T_a Γ₂ = -T_{-a}` with `Γ₂ = diag(1, -1)`, and at `E_c = -λ` the roles of `D_±` are exchanged.
* `lem:dimertransfer`: for every dimer configuration `v` on `[m, n-1]` (of either alignment) and every
  `z ∈ ℂ` with `|z - λ| ≤ 1`, `‖T^z(n, m)‖ ≤ K₁ e^{K₂ |z - λ| (n - m)}`, `K₁ = κ (2 + 2λ)²`,
  `K₂ = κ (3 + 4λ)/2`.

Norms of matrices are `ℓ²` operator norms (`Matrix.Norms.L2Operator`).
-/

open Matrix
open scoped Matrix.Norms.L2Operator

noncomputable section

namespace DFL

/-- `T_a = !![a, -1; 1, 0]`. -/
def tmC (a : ℂ) : Matrix (Fin 2) (Fin 2) ℂ := !![a, -1; 1, 0]

/-- `D_+(z) = T_{λ - z}²`. -/
def dimerDp (lam : ℝ) (z : ℂ) : Matrix (Fin 2) (Fin 2) ℂ := tmC (lam - z) * tmC (lam - z)

/-- `D_-(z) = T_{-λ - z}²`. -/
def dimerDm (lam : ℝ) (z : ℂ) : Matrix (Fin 2) (Fin 2) ℂ := tmC (-lam - z) * tmC (-lam - z)

/-- The rotation `R_φ`. -/
def rotC (φ : ℝ) : Matrix (Fin 2) (Fin 2) ℂ :=
  !![(Real.cos φ : ℂ), -(Real.sin φ : ℂ); (Real.sin φ : ℂ), (Real.cos φ : ℂ)]

/-- `σ_x`. -/
def sigmaXC : Matrix (Fin 2) (Fin 2) ℂ := !![0, 1; 1, 0]

/-- `Γ₂ = diag(1, -1)`. -/
def gamma2C : Matrix (Fin 2) (Fin 2) ℂ := !![1, 0; 0, -1]

/-- The angle `θ ∈ (π/2, π)` with `cos θ = -λ`. -/
def dimerTheta (lam : ℝ) : ℝ := Real.arccos (-lam)

/-- `S = !![sin θ, cos θ; 0, 1]`. -/
def dimerS (lam : ℝ) : Matrix (Fin 2) (Fin 2) ℂ :=
  !![(Real.sin (dimerTheta lam) : ℂ), (Real.cos (dimerTheta lam) : ℂ); 0, 1]

/-- `κ = ((1+λ)/(1-λ))^{1/2}`. -/
def dimerKappa (lam : ℝ) : ℝ := √((1 + lam) / (1 - lam))

/-! ### (i) -/

/-- **Commutator of the polymer transfer matrices** (paper: Lemma `lem:dimercritical`(i)). -/
theorem dimerD_commutator (lam : ℝ) (z : ℂ) :
    dimerDp lam z * dimerDm lam z - dimerDm lam z * dimerDp lam z =
      (2 * lam * (lam ^ 2 - z ^ 2)) • sigmaXC := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [dimerDp, dimerDm, tmC, sigmaXC] <;> ring1

/-- `D_+(E)` and `D_-(E)` commute iff `E = ±λ` (paper: Lemma `lem:dimercritical`(i)). -/
theorem dimerD_commute_iff (lam : ℝ) (hlam : 0 < lam) (E : ℝ) :
    Commute (dimerDp lam E) (dimerDm lam E) ↔ E = lam ∨ E = -lam := by
  have key := dimerD_commutator lam E
  constructor
  · intro h
    have h2 : (2 * (lam : ℂ) * ((lam : ℂ) ^ 2 - (E : ℂ) ^ 2)) • sigmaXC = 0 := by
      rw [← key, h.eq, sub_self]
    have h3 : (2 * (lam : ℂ) * ((lam : ℂ) ^ 2 - (E : ℂ) ^ 2)) = 0 := by
      have := congrFun (congrFun h2 0) 1
      simpa [sigmaXC] using this
    have h4 : 2 * lam * (lam ^ 2 - E ^ 2) = 0 := by exact_mod_cast h3
    have h5 : (E - lam) * (E + lam) = 0 := by
      rcases mul_eq_zero.mp h4 with h | h
      · exfalso; linarith
      · linear_combination -h
    rcases mul_eq_zero.mp h5 with h | h
    · left; linarith
    · right; linarith
  · intro h
    show dimerDp lam E * dimerDm lam E = dimerDm lam E * dimerDp lam E
    rw [← sub_eq_zero, key]
    rcases h with h | h <;> subst h <;> simp

/-! ### (iii) -/

/-- `θ ∈ (π/2, π)`, `cos θ = -λ`, `sin θ = √(1 - λ²)`. -/
theorem dimerTheta_spec (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    Real.pi / 2 < dimerTheta lam ∧ dimerTheta lam < Real.pi ∧
      Real.cos (dimerTheta lam) = -lam ∧ Real.sin (dimerTheta lam) = √(1 - lam ^ 2) := by
  unfold dimerTheta
  refine ⟨Real.pi_div_two_lt_arccos.mpr (by linarith), Real.arccos_lt_pi.mpr (by linarith),
    Real.cos_arccos (by linarith) (by linarith), ?_⟩
  rw [Real.sin_arccos, neg_sq]

namespace DimerTransferAux

/-- `sin θ > 0`. -/
lemma sin_pos (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) : 0 < Real.sin (dimerTheta lam) := by
  rw [(dimerTheta_spec lam h0 h1).2.2.2]
  exact Real.sqrt_pos.mpr (by nlinarith)

/-- `sin² θ = 1 - λ²`. -/
lemma sin_sq (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    Real.sin (dimerTheta lam) ^ 2 = 1 - lam ^ 2 := by
  rw [(dimerTheta_spec lam h0 h1).2.2.2, Real.sq_sqrt (by nlinarith)]

/-- `sin² θ = 1 - λ²`, as complex numbers. -/
lemma sin_sq_C (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    ((Real.sin (dimerTheta lam) : ℝ) : ℂ) ^ 2 = 1 - (lam : ℂ) ^ 2 := by
  exact_mod_cast sin_sq lam h0 h1

/-- If `A P = P B` with `P` invertible, then `P⁻¹ A P = B`. -/
lemma conj_of_mul_eq {A B P : Matrix (Fin 2) (Fin 2) ℂ} (hP : IsUnit P.det)
    (h : A * P = P * B) : P⁻¹ * A * P = B := by
  rw [Matrix.mul_assoc, h, ← Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hP, Matrix.one_mul]

/-- If `A P = P B` with `P` invertible, then `P⁻¹ A² P = B²`. -/
lemma conj_sq_of_mul_eq {A B P : Matrix (Fin 2) (Fin 2) ℂ} (hP : IsUnit P.det)
    (h : A * P = P * B) : P⁻¹ * (A * A) * P = B * B := by
  apply conj_of_mul_eq hP
  rw [Matrix.mul_assoc, h, ← Matrix.mul_assoc, h, Matrix.mul_assoc]

/-- `R_a R_b = R_{a+b}`. -/
lemma rotC_mul (a b : ℝ) : rotC a * rotC b = rotC (a + b) := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [rotC, Real.cos_add, Real.sin_add] <;> ring1

/-- `T_{-2λ} S = S R_θ`. -/
lemma tm_mul_S (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    tmC (-2 * lam) * dimerS lam = dimerS lam * rotC (dimerTheta lam) := by
  have hc := (dimerTheta_spec lam h0 h1).2.2.1
  have hs := sin_sq_C lam h0 h1
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [tmC, dimerS, rotC, hc, -Complex.ofReal_sin, -Complex.ofReal_cos] <;>
    first | ring1 | linear_combination hs

/-- `T_{2λ} (Γ₂ S) = (Γ₂ S) (-R_θ)`. -/
lemma tm_mul_gS (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    tmC (2 * lam) * (gamma2C * dimerS lam) =
      (gamma2C * dimerS lam) * (-rotC (dimerTheta lam)) := by
  have hc := (dimerTheta_spec lam h0 h1).2.2.1
  have hs := sin_sq_C lam h0 h1
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [tmC, dimerS, rotC, gamma2C, hc, -Complex.ofReal_sin, -Complex.ofReal_cos] <;>
    first | ring1 | linear_combination -hs

end DimerTransferAux

open DimerTransferAux in
/-- `S` is invertible. -/
theorem dimerS_det (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) : (dimerS lam).det ≠ 0 := by
  rw [dimerS, Matrix.det_fin_two_of]
  have := sin_pos lam h0 h1
  simp only [mul_one, mul_zero, sub_zero]
  exact_mod_cast this.ne'

namespace DimerTransferAux

/-- `det S` is a unit. -/
lemma S_isUnit (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) : IsUnit (dimerS lam).det :=
  isUnit_iff_ne_zero.mpr (dimerS_det lam h0 h1)

/-- `det (Γ₂ S)` is a unit. -/
lemma gS_isUnit (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    IsUnit (gamma2C * dimerS lam).det := by
  rw [Matrix.det_mul]
  refine (isUnit_iff_ne_zero.mpr ?_).mul (S_isUnit lam h0 h1)
  simp [gamma2C, Matrix.det_fin_two_of]

/-- `D_+(λ) = -1`. -/
lemma dimerDp_self (lam : ℝ) : dimerDp lam lam = -1 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [dimerDp, tmC]

/-- `D_-(-λ) = -1`. -/
lemma dimerDm_neg_self (lam : ℝ) : dimerDm lam (-lam) = -1 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [dimerDm, tmC]

/-- `R_π = -1`. -/
lemma rotC_pi : rotC Real.pi = -1 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [rotC]

end DimerTransferAux

/-- **Rotation normal form at `E_c = λ`** (paper: Lemma `lem:dimercritical`(iii), eq. `eq:dimerrot`):
`S⁻¹ D_+(λ) S = -1 = R_π`. -/
theorem dimerS_conj_plus (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    (dimerS lam)⁻¹ * dimerDp lam lam * dimerS lam = -1 ∧ rotC Real.pi = -1 := by
  refine ⟨?_, DimerTransferAux.rotC_pi⟩
  rw [DimerTransferAux.dimerDp_self, mul_neg, mul_one, neg_mul,
    Matrix.nonsing_inv_mul _ (DimerTransferAux.S_isUnit lam h0 h1)]

/-- `S⁻¹ T_{-2λ} S = R_θ` (paper: eq. `eq:dimerrot`). -/
theorem dimerS_conj_tm (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    (dimerS lam)⁻¹ * tmC (-2 * lam) * dimerS lam = rotC (dimerTheta lam) :=
  DimerTransferAux.conj_of_mul_eq (DimerTransferAux.S_isUnit lam h0 h1)
    (DimerTransferAux.tm_mul_S lam h0 h1)

/-- `S⁻¹ D_-(λ) S = R_{2θ}` (paper: eq. `eq:dimerrot`). -/
theorem dimerS_conj_minus (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    (dimerS lam)⁻¹ * dimerDm lam lam * dimerS lam = rotC (2 * dimerTheta lam) := by
  have h := DimerTransferAux.conj_sq_of_mul_eq (DimerTransferAux.S_isUnit lam h0 h1)
    (DimerTransferAux.tm_mul_S lam h0 h1)
  rw [DimerTransferAux.rotC_mul, ← two_mul] at h
  rw [← h, dimerDm]
  congr 4 <;> ring1

/-- **The averaged rotation phases** (paper: eq. `eq:dimerangles`), with `η₊ = π`, `η₋ = 2θ` and equal
probabilities. -/
theorem dimer_angles (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    (1 / 2 : ℂ) * (Complex.exp (2 * Real.pi * Complex.I) +
        Complex.exp (2 * (2 * dimerTheta lam) * Complex.I)) =
      (1 + Complex.exp (4 * dimerTheta lam * Complex.I)) / 2 ∧
    ‖(1 + Complex.exp (4 * dimerTheta lam * Complex.I)) / 2‖ = |1 - 2 * lam ^ 2| ∧
    |1 - 2 * lam ^ 2| < 1 ∧
    ((1 / 2 : ℂ) * (Complex.exp (4 * Real.pi * Complex.I) +
        Complex.exp (4 * (2 * dimerTheta lam) * Complex.I)) = 1 ↔ lam = 1 / √2) := by
  obtain ⟨hθ1, hθ2, hc, -⟩ := dimerTheta_spec lam h0 h1
  set θ := dimerTheta lam with hθ
  refine ⟨?_, ?_, ?_, ?_⟩
  · rw [Complex.exp_two_pi_mul_I]
    ring_nf
  · have h4 : (4 * (θ : ℂ) * Complex.I) = ((4 * θ : ℝ) : ℂ) * Complex.I := by push_cast; ring
    rw [h4, Complex.exp_mul_I, ← Complex.ofReal_cos, ← Complex.ofReal_sin]
    have hcos4 : Real.cos (4 * θ) = 2 * (2 * lam ^ 2 - 1) ^ 2 - 1 := by
      rw [show 4 * θ = 2 * (2 * θ) by ring, Real.cos_two_mul, Real.cos_two_mul, hc]
      ring
    have hw : (1 + ((Real.cos (4 * θ) : ℂ) + (Real.sin (4 * θ) : ℂ) * Complex.I)) / 2 =
        (((1 + Real.cos (4 * θ)) / 2 : ℝ) : ℂ) + ((Real.sin (4 * θ) / 2 : ℝ) : ℂ) * Complex.I := by
      push_cast; ring
    rw [hw]
    refine (sq_eq_sq₀ (norm_nonneg _) (abs_nonneg _)).mp ?_
    rw [sq_abs, Complex.sq_norm, Complex.normSq_add_mul_I]
    linear_combination (1 / 4 : ℝ) * Real.sin_sq_add_cos_sq (4 * θ) + (1 / 2 : ℝ) * hcos4
  · exact abs_lt.mpr ⟨by nlinarith, by nlinarith⟩
  · have h4pi : Complex.exp (4 * Real.pi * Complex.I) = 1 :=
      Complex.exp_eq_one_iff.mpr ⟨2, by push_cast; ring⟩
    rw [h4pi]
    have hiff : (1 / 2 : ℂ) * (1 + Complex.exp (4 * (2 * (θ : ℂ)) * Complex.I)) = 1 ↔
        Complex.exp (4 * (2 * (θ : ℂ)) * Complex.I) = 1 := by
      constructor
      · intro h; linear_combination 2 * h
      · intro h; rw [h]; norm_num
    rw [hiff, Complex.exp_eq_one_iff]
    have hsq2 : √2 * √2 = 2 := Real.mul_self_sqrt (by norm_num)
    have hsq2pos : 0 < √2 := by positivity
    have hinv : 1 / √2 = √2 / 2 := by
      field_simp
      linarith [hsq2]
    have hc34 : Real.cos (3 * Real.pi / 4) = -(√2 / 2) := by
      rw [show 3 * Real.pi / 4 = Real.pi - Real.pi / 4 by ring, Real.cos_pi_sub,
        Real.cos_pi_div_four]
    constructor
    · rintro ⟨n, hn⟩
      have hn' : (8 * θ : ℂ) = ((2 * Real.pi * n : ℝ) : ℂ) := by
        have := mul_right_cancel₀ Complex.I_ne_zero
          (show (8 * (θ : ℂ)) * Complex.I = (2 * Real.pi * n : ℂ) * Complex.I by
            rw [show (8 * (θ : ℂ)) * Complex.I = 4 * (2 * (θ : ℂ)) * Complex.I by ring, hn]; ring)
        rw [this]; push_cast; ring
      have hn'' : 8 * θ = 2 * Real.pi * n := by exact_mod_cast hn'
      have hpi := Real.pi_pos
      have hn1 : (2 : ℝ) < n := by nlinarith
      have hn2 : (n : ℝ) < 4 := by nlinarith
      have hn3 : n = 3 := by
        have a : (2 : ℤ) < n := by exact_mod_cast hn1
        have b : n < (4 : ℤ) := by exact_mod_cast hn2
        omega
      subst hn3
      have hθ34 : θ = 3 * Real.pi / 4 := by push_cast at hn''; linarith
      rw [hθ34, hc34] at hc
      rw [hinv]; linarith
    · intro hl
      have hθ34 : θ = 3 * Real.pi / 4 := by
        rw [hθ, dimerTheta]
        apply Real.arccos_eq_of_eq_cos (by positivity) (by nlinarith [Real.pi_pos])
        rw [hc34, hl, hinv]
      refine ⟨3, ?_⟩
      rw [hθ34]; push_cast; ring

namespace DimerTransferAux

/-- The squared Euclidean norm of a vector in `ℂ²`. -/
lemma euclid_norm_sq (x : EuclideanSpace ℂ (Fin 2)) : ‖x‖ ^ 2 = ‖x 0‖ ^ 2 + ‖x 1‖ ^ 2 := by
  rw [EuclideanSpace.norm_sq_eq, Fin.sum_univ_two]

/-- Upper bound for the squared `ℓ²` operator norm of a `2 × 2` matrix. -/
lemma norm_sq_le_of_forall {A : Matrix (Fin 2) (Fin 2) ℂ} {c : ℝ} (hc : 0 ≤ c)
    (h : ∀ v : Fin 2 → ℂ, ‖(A *ᵥ v) 0‖ ^ 2 + ‖(A *ᵥ v) 1‖ ^ 2 ≤ c * (‖v 0‖ ^ 2 + ‖v 1‖ ^ 2)) :
    ‖A‖ ^ 2 ≤ c := by
  have hA : ‖A‖ ≤ √c := by
    rw [Matrix.cstar_norm_def]
    refine ContinuousLinearMap.opNorm_le_bound _ (Real.sqrt_nonneg c) fun x => ?_
    rw [← sq_le_sq₀ (norm_nonneg _) (by positivity), mul_pow, Real.sq_sqrt hc, euclid_norm_sq,
      euclid_norm_sq]
    simpa using h (WithLp.ofLp x)
  calc ‖A‖ ^ 2 ≤ (√c) ^ 2 := by gcongr
    _ = c := Real.sq_sqrt hc

/-- Lower bound for the squared `ℓ²` operator norm via a test vector. -/
lemma sq_le_norm_sq_mul (A : Matrix (Fin 2) (Fin 2) ℂ) (v : Fin 2 → ℂ) :
    ‖(A *ᵥ v) 0‖ ^ 2 + ‖(A *ᵥ v) 1‖ ^ 2 ≤ ‖A‖ ^ 2 * (‖v 0‖ ^ 2 + ‖v 1‖ ^ 2) := by
  have h := (Matrix.toEuclideanCLM (n := Fin 2) (𝕜 := ℂ) A).le_opNorm (WithLp.toLp 2 v)
  rw [← Matrix.cstar_norm_def] at h
  have h2 := pow_le_pow_left₀ (norm_nonneg _) h 2
  rw [mul_pow, euclid_norm_sq, euclid_norm_sq] at h2
  simpa using h2

/-- `‖x‖² = (re x)² + (im x)²`. -/
lemma norm_sq_re_im (x : ℂ) : ‖x‖ ^ 2 = x.re ^ 2 + x.im ^ 2 := by
  rw [Complex.sq_norm, Complex.normSq_apply]; ring

/-- A matrix with real entries acting on `ℂ²`. -/
lemma mulVec_real (a b c d : ℝ) (v : Fin 2 → ℂ) :
    (!![(a : ℂ), (b : ℂ); (c : ℂ), (d : ℂ)] *ᵥ v) 0 = a * v 0 + b * v 1 ∧
    (!![(a : ℂ), (b : ℂ); (c : ℂ), (d : ℂ)] *ᵥ v) 1 = c * v 0 + d * v 1 := by
  simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two]

/-- `‖a x + b y‖²` for real `a, b`, in terms of real and imaginary parts. -/
lemma norm_sq_real_comb (a b : ℝ) (x y : ℂ) :
    ‖(a : ℂ) * x + (b : ℂ) * y‖ ^ 2 = (a * x.re + b * y.re) ^ 2 + (a * x.im + b * y.im) ^ 2 := by
  rw [norm_sq_re_im]; simp

/-- `S⁻¹ = !![1/s, λ/s; 0, 1]`. -/
lemma S_inv_eq (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    (dimerS lam)⁻¹ = !![((1 / Real.sin (dimerTheta lam) : ℝ) : ℂ),
      ((lam / Real.sin (dimerTheta lam) : ℝ) : ℂ); ((0 : ℝ) : ℂ), ((1 : ℝ) : ℂ)] := by
  apply Matrix.inv_eq_left_inv
  have hs := (sin_pos lam h0 h1).ne'
  have hc := (dimerTheta_spec lam h0 h1).2.2.1
  have hs' : ((Real.sin (dimerTheta lam) : ℝ) : ℂ) ≠ 0 := by exact_mod_cast hs
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [dimerS, hc, hs', -Complex.ofReal_sin, -Complex.ofReal_cos]
  ring1

/-- `‖S‖² = 1 + λ`. -/
lemma norm_S_sq (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) : ‖dimerS lam‖ ^ 2 = 1 + lam := by
  have hc := (dimerTheta_spec lam h0 h1).2.2.1
  have hs := sin_sq lam h0 h1
  have hsq := (dimerTheta_spec lam h0 h1).2.2.2
  set s := Real.sin (dimerTheta lam) with hsdef
  have hS : dimerS lam = !![(s : ℂ), ((-lam : ℝ) : ℂ); ((0 : ℝ) : ℂ), ((1 : ℝ) : ℂ)] := by
    rw [dimerS, hc, Complex.ofReal_zero, Complex.ofReal_one]
  set p := √(1 + lam) with hpdef
  set q := √(1 - lam) with hqdef
  have hp : p ^ 2 = 1 + lam := Real.sq_sqrt (by linarith)
  have hq : q ^ 2 = 1 - lam := Real.sq_sqrt (by linarith)
  have hpq : p * q = s := by
    rw [hpdef, hqdef, ← Real.sqrt_mul (by linarith), hsq]
    congr 1; ring
  apply le_antisymm
  · apply norm_sq_le_of_forall (by linarith)
    intro v
    rw [hS, (mulVec_real _ _ _ _ v).1, (mulVec_real _ _ _ _ v).2, norm_sq_real_comb,
      norm_sq_real_comb, norm_sq_re_im (v 0), norm_sq_re_im (v 1)]
    set x1 := (v 0).re
    set x2 := (v 0).im
    set y1 := (v 1).re
    set y2 := (v 1).im
    have key := mul_nonneg h0.le (add_nonneg (sq_nonneg (p * x1 + q * y1))
      (sq_nonneg (p * x2 + q * y2)))
    have e : (1 + lam) * (x1 ^ 2 + x2 ^ 2 + (y1 ^ 2 + y2 ^ 2)) -
        ((s * x1 + -lam * y1) ^ 2 + (s * x2 + -lam * y2) ^ 2 +
          ((0 * x1 + 1 * y1) ^ 2 + (0 * x2 + 1 * y2) ^ 2)) =
        lam * ((p * x1 + q * y1) ^ 2 + (p * x2 + q * y2) ^ 2) := by
      linear_combination (-(x1 ^ 2 + x2 ^ 2)) * hs - lam * (x1 ^ 2 + x2 ^ 2) * hp -
        2 * lam * (x1 * y1 + x2 * y2) * hpq - lam * (y1 ^ 2 + y2 ^ 2) * hq
    linarith
  · have h := sq_le_norm_sq_mul (dimerS lam) ![(q : ℂ), -(p : ℂ)]
    have hreal : s * q + lam * p = p := by rw [← hpq]; linear_combination p * hq
    have e0 : (dimerS lam *ᵥ ![(q : ℂ), -(p : ℂ)]) 0 = (p : ℂ) := by
      rw [hS, (mulVec_real _ _ _ _ _).1]
      have := congrArg (fun x : ℝ => (x : ℂ)) hreal
      push_cast at this
      simp only [Matrix.cons_val_zero, Matrix.cons_val_one]
      push_cast
      linear_combination this
    have e1 : (dimerS lam *ᵥ ![(q : ℂ), -(p : ℂ)]) 1 = -(p : ℂ) := by
      rw [hS, (mulVec_real _ _ _ _ _).2]
      simp
    rw [e0, e1] at h
    simp only [norm_neg, Complex.norm_real, Real.norm_eq_abs, sq_abs, Matrix.cons_val_zero,
      Matrix.cons_val_one] at h
    rw [hp, hq] at h
    nlinarith [h]

/-- `‖S⁻¹‖² = 1/(1 - λ)`. -/
lemma norm_S_inv_sq (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    ‖(dimerS lam)⁻¹‖ ^ 2 = 1 / (1 - lam) := by
  have hs := sin_sq lam h0 h1
  have hspos := sin_pos lam h0 h1
  have hSi := S_inv_eq lam h0 h1
  set s := Real.sin (dimerTheta lam) with hsdef
  have hsne : s ≠ 0 := hspos.ne'
  have h1l : 1 - lam ≠ 0 := by linarith
  have h1l' : 0 < 1 - lam := by linarith
  apply le_antisymm
  · apply norm_sq_le_of_forall (by positivity)
    intro v
    rw [hSi, (mulVec_real _ _ _ _ v).1, (mulVec_real _ _ _ _ v).2, norm_sq_real_comb,
      norm_sq_real_comb, norm_sq_re_im (v 0), norm_sq_re_im (v 1)]
    set x1 := (v 0).re
    set x2 := (v 0).im
    set y1 := (v 1).re
    set y2 := (v 1).im
    have hA : (x1 + lam * y1) ^ 2 + (x2 + lam * y2) ^ 2 + (1 - lam ^ 2) * (y1 ^ 2 + y2 ^ 2) ≤
        (1 + lam) * (x1 ^ 2 + x2 ^ 2 + (y1 ^ 2 + y2 ^ 2)) := by
      nlinarith [mul_nonneg h0.le (sq_nonneg (x1 - y1)), mul_nonneg h0.le (sq_nonneg (x2 - y2))]
    calc (1 / s * x1 + lam / s * y1) ^ 2 + (1 / s * x2 + lam / s * y2) ^ 2 +
          ((0 * x1 + 1 * y1) ^ 2 + (0 * x2 + 1 * y2) ^ 2)
        = (1 / s) ^ 2 * ((x1 + lam * y1) ^ 2 + (x2 + lam * y2) ^ 2 +
            s ^ 2 * (y1 ^ 2 + y2 ^ 2)) := by
          field_simp
          ring
      _ ≤ (1 / s) ^ 2 * ((1 + lam) * (x1 ^ 2 + x2 ^ 2 + (y1 ^ 2 + y2 ^ 2))) := by
          rw [hs]; gcongr
      _ = 1 / (1 - lam) * (x1 ^ 2 + x2 ^ 2 + (y1 ^ 2 + y2 ^ 2)) := by
          rw [div_pow, one_pow, hs]
          have : 1 - lam ^ 2 = (1 - lam) * (1 + lam) := by ring
          rw [this]
          field_simp
  · have h := sq_le_norm_sq_mul (dimerS lam)⁻¹ ![1, 1]
    have e0 : ((dimerS lam)⁻¹ *ᵥ ![(1 : ℂ), 1]) 0 = (((1 + lam) / s : ℝ) : ℂ) := by
      rw [hSi, (mulVec_real _ _ _ _ _).1]
      simp only [Matrix.cons_val_zero, Matrix.cons_val_one, mul_one]
      push_cast
      ring
    have e1 : ((dimerS lam)⁻¹ *ᵥ ![(1 : ℂ), 1]) 1 = 1 := by
      rw [hSi, (mulVec_real _ _ _ _ _).2]
      simp
    rw [e0, e1] at h
    simp only [Complex.norm_real, Real.norm_eq_abs, sq_abs, norm_one, Matrix.cons_val_zero,
      Matrix.cons_val_one, one_pow] at h
    have e2 : ((1 + lam) / s) ^ 2 = (1 + lam) / (1 - lam) := by
      rw [div_pow, hs]
      have : 1 - lam ^ 2 = (1 - lam) * (1 + lam) := by ring
      rw [this]
      field_simp
    rw [e2] at h
    have e3 : 1 / (1 - lam) = ((1 + lam) / (1 - lam) + 1) / 2 := by
      field_simp
      ring
    rw [e3]
    linarith

end DimerTransferAux

/-- **The condition number of `S`** (paper: Lemma `lem:dimercritical`(iii)):
`κ = ‖S‖ ‖S⁻¹‖ = ((1+λ)/(1-λ))^{1/2}`. -/
theorem dimerS_kappa (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    ‖dimerS lam‖ * ‖(dimerS lam)⁻¹‖ = dimerKappa lam := by
  have hS := DimerTransferAux.norm_S_sq lam h0 h1
  have hSi := DimerTransferAux.norm_S_inv_sq lam h0 h1
  rw [dimerKappa, ← Real.sqrt_sq (mul_nonneg (norm_nonneg (dimerS lam))
    (norm_nonneg (dimerS lam)⁻¹)), mul_pow, hS, hSi, mul_one_div]

/-! ### (iv) -/

/-- `Γ₂ T_a Γ₂ = -T_{-a}` (paper: proof of Lemma `lem:dimercritical`(iv)). -/
theorem gamma2_tmC (a : ℂ) : gamma2C * tmC a * gamma2C = -tmC (-a) := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [gamma2C, tmC]

/-- **Rotation normal form at `E_c = -λ`** (paper: Lemma `lem:dimercritical`(iv)):
`(Γ₂S)⁻¹ D_-(-λ) Γ₂S = R_π` and `(Γ₂S)⁻¹ D_+(-λ) Γ₂S = R_{2θ}`. -/
theorem dimerS_conj_neg (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    (gamma2C * dimerS lam)⁻¹ * dimerDm lam (-lam) * (gamma2C * dimerS lam) = rotC Real.pi ∧
      (gamma2C * dimerS lam)⁻¹ * dimerDp lam (-lam) * (gamma2C * dimerS lam) =
        rotC (2 * dimerTheta lam) := by
  have hU := DimerTransferAux.gS_isUnit lam h0 h1
  refine ⟨?_, ?_⟩
  · rw [DimerTransferAux.rotC_pi, DimerTransferAux.dimerDm_neg_self, mul_neg, mul_one, neg_mul,
      Matrix.nonsing_inv_mul _ hU]
  · have h := DimerTransferAux.conj_sq_of_mul_eq hU (DimerTransferAux.tm_mul_gS lam h0 h1)
    rw [neg_mul_neg, DimerTransferAux.rotC_mul, ← two_mul] at h
    rw [← h, dimerDp]
    congr 4 <;> ring1

/-! ### Lemma `lem:dimertransfer` -/

/-- The transfer matrix `T^z(m + L, m) = T_{v(m+L-1) - z} ⋯ T_{v(m) - z}` across the sites
`m, …, m + L - 1`. -/
def tmProd (v : ℤ → ℝ) (z : ℂ) (m : ℤ) : ℕ → Matrix (Fin 2) (Fin 2) ℂ
  | 0 => 1
  | L + 1 => tmC (v (m + L) - z) * tmProd v z m L

/-- `v` is a dimer configuration with alignment `r` on the sites `m, …, m + L - 1`: `v = ±λ`, and `v`
takes equal values on the two sites of every dimer `{2k - r, 2k - r + 1}` contained in these sites. -/
def IsDimerCfg (lam : ℝ) (r : ℤ) (v : ℤ → ℝ) (m : ℤ) (L : ℕ) : Prop :=
  (∀ y : ℤ, m ≤ y → y < m + L → v y = lam ∨ v y = -lam) ∧
    ∀ k : ℤ, m ≤ 2 * k - r → 2 * k - r + 1 < m + L → v (2 * k - r) = v (2 * k - r + 1)

namespace DimerTransferAux

/-- `E₁₁ = diag(1, 0)`. -/
def e11 : Matrix (Fin 2) (Fin 2) ℂ := !![1, 0; 0, 0]

/-- `‖E₁₁‖ ≤ 1`. -/
lemma norm_e11_le : ‖e11‖ ≤ 1 := by
  have h : ‖e11‖ ^ 2 ≤ 1 := norm_sq_le_of_forall zero_le_one fun v => by
    simp [e11, Matrix.mulVec, dotProduct, Fin.sum_univ_two]
  nlinarith [norm_nonneg e11]

/-- `T_a - T_b = (a - b) E₁₁`. -/
lemma tmC_sub (a b : ℂ) : tmC a - tmC b = (a - b) • e11 := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [tmC, e11]

/-- `‖T_a - T_b‖ ≤ |a - b|`. -/
lemma norm_tmC_sub_le (a b : ℂ) : ‖tmC a - tmC b‖ ≤ ‖a - b‖ := by
  rw [tmC_sub, norm_smul]
  calc ‖a - b‖ * ‖e11‖ ≤ ‖a - b‖ * 1 := by gcongr; exact norm_e11_le
    _ = ‖a - b‖ := mul_one _

/-- `‖T_0‖ ≤ 1`. -/
lemma norm_tmC_zero_le : ‖tmC 0‖ ≤ 1 := by
  have h : ‖tmC 0‖ ^ 2 ≤ 1 := norm_sq_le_of_forall zero_le_one fun v => by
    simp [tmC, Matrix.mulVec, dotProduct, Fin.sum_univ_two]
    linarith
  nlinarith [norm_nonneg (tmC 0)]

/-- `‖T_a‖ ≤ |a| + 1`. -/
lemma norm_tmC_le (a : ℂ) : ‖tmC a‖ ≤ ‖a‖ + 1 := by
  calc ‖tmC a‖ = ‖(tmC a - tmC 0) + tmC 0‖ := by rw [sub_add_cancel]
    _ ≤ ‖tmC a - tmC 0‖ + ‖tmC 0‖ := norm_add_le _ _
    _ ≤ ‖a - 0‖ + 1 := add_le_add (norm_tmC_sub_le a 0) norm_tmC_zero_le
    _ = ‖a‖ + 1 := by rw [sub_zero]

/-- `‖R_φ‖ ≤ 1`. -/
lemma norm_rotC_le (φ : ℝ) : ‖rotC φ‖ ≤ 1 := by
  have h : ‖rotC φ‖ ^ 2 ≤ 1 := norm_sq_le_of_forall zero_le_one fun v => by
    have hm := mulVec_real (Real.cos φ) (-Real.sin φ) (Real.sin φ) (Real.cos φ) v
    have hR : rotC φ = !![((Real.cos φ : ℝ) : ℂ), ((-Real.sin φ : ℝ) : ℂ);
        ((Real.sin φ : ℝ) : ℂ), ((Real.cos φ : ℝ) : ℂ)] := by
      rw [rotC, Complex.ofReal_neg]
    rw [hR, hm.1, hm.2, norm_sq_real_comb, norm_sq_real_comb, norm_sq_re_im (v 0),
      norm_sq_re_im (v 1)]
    apply le_of_eq
    linear_combination ((v 0).re ^ 2 + (v 0).im ^ 2 + ((v 1).re ^ 2 + (v 1).im ^ 2)) *
      Real.cos_sq_add_sin_sq φ
  nlinarith [norm_nonneg (rotC φ)]

/-- `‖1‖ ≤ 1`. -/
lemma norm_one_le : ‖(1 : Matrix (Fin 2) (Fin 2) ℂ)‖ ≤ 1 := by
  have h : ‖(1 : Matrix (Fin 2) (Fin 2) ℂ)‖ ^ 2 ≤ 1 :=
    norm_sq_le_of_forall zero_le_one fun v => by simp
  nlinarith [norm_nonneg (1 : Matrix (Fin 2) (Fin 2) ℂ)]

/-- `κ ≥ 1`. -/
lemma kappa_ge_one (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) : 1 ≤ dimerKappa lam := by
  rw [dimerKappa, Real.one_le_sqrt, le_div_iff₀ (by linarith)]
  linarith

/-- `‖S⁻¹ X S‖ ≤ κ ‖X‖`. -/
lemma norm_conj_le (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (X : Matrix (Fin 2) (Fin 2) ℂ) :
    ‖(dimerS lam)⁻¹ * X * dimerS lam‖ ≤ dimerKappa lam * ‖X‖ := by
  rw [← dimerS_kappa lam h0 h1]
  calc ‖(dimerS lam)⁻¹ * X * dimerS lam‖ ≤ ‖(dimerS lam)⁻¹ * X‖ * ‖dimerS lam‖ :=
        norm_mul_le _ _
    _ ≤ ‖(dimerS lam)⁻¹‖ * ‖X‖ * ‖dimerS lam‖ := by gcongr; exact norm_mul_le _ _
    _ = ‖dimerS lam‖ * ‖(dimerS lam)⁻¹‖ * ‖X‖ := by ring

/-- `‖X‖ ≤ κ ‖S⁻¹ X S‖`. -/
lemma norm_le_conj (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (X : Matrix (Fin 2) (Fin 2) ℂ) :
    ‖X‖ ≤ dimerKappa lam * ‖(dimerS lam)⁻¹ * X * dimerS lam‖ := by
  have hU := S_isUnit lam h0 h1
  have hX : dimerS lam * ((dimerS lam)⁻¹ * X * dimerS lam) * (dimerS lam)⁻¹ = X := by
    rw [Matrix.mul_assoc (dimerS lam)⁻¹, Matrix.mul_nonsing_inv_cancel_left _ _ hU,
      Matrix.mul_nonsing_inv_cancel_right _ _ hU]
  rw [← dimerS_kappa lam h0 h1]
  calc ‖X‖ = ‖dimerS lam * ((dimerS lam)⁻¹ * X * dimerS lam) * (dimerS lam)⁻¹‖ := by rw [hX]
    _ ≤ ‖dimerS lam * ((dimerS lam)⁻¹ * X * dimerS lam)‖ * ‖(dimerS lam)⁻¹‖ := norm_mul_le _ _
    _ ≤ ‖dimerS lam‖ * ‖(dimerS lam)⁻¹ * X * dimerS lam‖ * ‖(dimerS lam)⁻¹‖ := by
        gcongr; exact norm_mul_le _ _
    _ = ‖dimerS lam‖ * ‖(dimerS lam)⁻¹‖ * ‖(dimerS lam)⁻¹ * X * dimerS lam‖ := by ring

/-- `|w - λ| ≤ 2λ` for `w = ±λ`. -/
lemma norm_w_sub_le (lam : ℝ) (h0 : 0 < lam) {w : ℝ} (hw : w = lam ∨ w = -lam) :
    ‖(w : ℂ) - lam‖ ≤ 2 * lam := by
  rw [← Complex.ofReal_sub, Complex.norm_real, Real.norm_eq_abs, abs_le]
  rcases hw with hw | hw <;> subst hw <;> constructor <;> linarith

/-- A single site: `‖T_{w - z}‖ ≤ 2 + 2λ` for `w = ±λ` and `|z - λ| ≤ 1`. -/
lemma norm_tm_single_le (lam : ℝ) (h0 : 0 < lam) {w : ℝ} (hw : w = lam ∨ w = -lam) (z : ℂ)
    (hz : ‖z - lam‖ ≤ 1) : ‖tmC (w - z)‖ ≤ 2 + 2 * lam := by
  have hwl := norm_w_sub_le lam h0 hw
  have htri : ‖(w : ℂ) - z‖ ≤ ‖(w : ℂ) - lam‖ + ‖z - lam‖ := by
    calc ‖(w : ℂ) - z‖ ≤ ‖(w : ℂ) - lam‖ + ‖(lam : ℂ) - z‖ :=
          norm_sub_le_norm_sub_add_norm_sub _ _ _
      _ = ‖(w : ℂ) - lam‖ + ‖z - lam‖ := by rw [norm_sub_rev (lam : ℂ) z]
  calc ‖tmC (w - z)‖ ≤ ‖(w : ℂ) - z‖ + 1 := norm_tmC_le _
    _ ≤ 2 + 2 * lam := by linarith

/-- Abstract perturbation step: `‖S⁻¹ A² S‖ ≤ 1 + κ (a + b) d` if `‖S⁻¹ B² S‖ ≤ 1`, `‖A‖ ≤ a`,
`‖B‖ ≤ b`, `‖A - B‖ ≤ d`. -/
lemma conj_sq_perturb {S A B : Matrix (Fin 2) (Fin 2) ℂ} {κ a b d : ℝ}
    (hκ : ∀ X : Matrix (Fin 2) (Fin 2) ℂ, ‖S⁻¹ * X * S‖ ≤ κ * ‖X‖) (hκ0 : 0 ≤ κ)
    (hB : ‖S⁻¹ * (B * B) * S‖ ≤ 1) (ha : ‖A‖ ≤ a) (hb : ‖B‖ ≤ b) (hd : ‖A - B‖ ≤ d) :
    ‖S⁻¹ * (A * A) * S‖ ≤ 1 + κ * ((a + b) * d) := by
  have e1 : A * A - B * B = A * (A - B) + (A - B) * B := by noncomm_ring
  have hAB : ‖A * A - B * B‖ ≤ (a + b) * d := by
    rw [e1]
    calc ‖A * (A - B) + (A - B) * B‖ ≤ ‖A‖ * ‖A - B‖ + ‖A - B‖ * ‖B‖ :=
          (norm_add_le _ _).trans (add_le_add (norm_mul_le _ _) (norm_mul_le _ _))
      _ ≤ a * d + d * b := by gcongr <;> linarith [norm_nonneg A, norm_nonneg (A - B)]
      _ = (a + b) * d := by ring
  have e2 : S⁻¹ * (A * A) * S = S⁻¹ * (B * B) * S + S⁻¹ * (A * A - B * B) * S := by
    noncomm_ring
  rw [e2]
  calc ‖S⁻¹ * (B * B) * S + S⁻¹ * (A * A - B * B) * S‖
      ≤ ‖S⁻¹ * (B * B) * S‖ + ‖S⁻¹ * (A * A - B * B) * S‖ := norm_add_le _ _
    _ ≤ 1 + κ * ((a + b) * d) := by
        gcongr
        exact (hκ _).trans (by gcongr)

/-- One dimer: `‖S⁻¹ T_{w-z}² S‖ ≤ exp(κ(3+4λ)|z-λ|)` for `w = ±λ` and `|z - λ| ≤ 1`. -/
lemma conj_dimer_le (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) {w : ℝ} (hw : w = lam ∨ w = -lam)
    (z : ℂ) (hz : ‖z - lam‖ ≤ 1) :
    ‖(dimerS lam)⁻¹ * (tmC (w - z) * tmC (w - z)) * dimerS lam‖ ≤
      Real.exp (dimerKappa lam * (3 + 4 * lam) * ‖z - lam‖) := by
  have hB : ‖(dimerS lam)⁻¹ * (tmC (w - lam) * tmC (w - lam)) * dimerS lam‖ ≤ 1 := by
    rcases hw with hw | hw <;> rw [hw]
    · rw [show tmC ((lam : ℂ) - lam) * tmC ((lam : ℂ) - lam) = dimerDp lam lam from rfl,
        (dimerS_conj_plus lam h0 h1).1, norm_neg]
      exact norm_one_le
    · rw [show tmC (((-lam : ℝ) : ℂ) - lam) * tmC (((-lam : ℝ) : ℂ) - lam) = dimerDm lam lam by
          rw [dimerDm, Complex.ofReal_neg], dimerS_conj_minus lam h0 h1]
      exact norm_rotC_le _
  have hwl := norm_w_sub_le lam h0 hw
  have ha := norm_tm_single_le lam h0 hw z hz
  have hb : ‖tmC (w - lam)‖ ≤ 2 * lam + 1 :=
    (norm_tmC_le _).trans (by linarith)
  have hd : ‖tmC (w - z) - tmC (w - lam)‖ ≤ ‖z - lam‖ := by
    refine (norm_tmC_sub_le _ _).trans (le_of_eq ?_)
    rw [show ((w : ℂ) - z) - (w - lam) = -(z - lam) by ring, norm_neg]
  have hκ0 : 0 ≤ dimerKappa lam := Real.sqrt_nonneg _
  have h := conj_sq_perturb (norm_conj_le lam h0 h1) hκ0 hB ha hb hd
  calc _ ≤ 1 + dimerKappa lam * ((2 + 2 * lam + (2 * lam + 1)) * ‖z - lam‖) := h
    _ = dimerKappa lam * (3 + 4 * lam) * ‖z - lam‖ + 1 := by ring
    _ ≤ _ := Real.add_one_le_exp _

/-- `T^z(m + L₁ + L₂, m) = T^z(m + L₁ + L₂, m + L₁) T^z(m + L₁, m)`. -/
lemma tmProd_add (v : ℤ → ℝ) (z : ℂ) (m : ℤ) (L₁ L₂ : ℕ) :
    tmProd v z m (L₁ + L₂) = tmProd v z (m + L₁) L₂ * tmProd v z m L₁ := by
  induction L₂ with
  | zero => simp [tmProd]
  | succ n ih =>
    rw [← add_assoc, tmProd, ih, tmProd, Matrix.mul_assoc]
    have e : m + ((L₁ + n : ℕ) : ℤ) = m + (L₁ : ℤ) + (n : ℤ) := by push_cast; ring
    rw [e]

/-- Products over complete dimers starting at the first site of a dimer:
`‖S⁻¹ T^z(m₀ + 2j, m₀) S‖ ≤ exp(κ(3+4λ)|z-λ| j)`. -/
lemma aligned_bound (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (v : ℤ → ℝ) (z : ℂ)
    (hz : ‖z - lam‖ ≤ 1) (m₀ : ℤ) :
    ∀ j : ℕ, (∀ i : ℕ, i < j → (v (m₀ + 2 * i) = lam ∨ v (m₀ + 2 * i) = -lam) ∧
        v (m₀ + 2 * i + 1) = v (m₀ + 2 * i)) →
      ‖(dimerS lam)⁻¹ * tmProd v z m₀ (2 * j) * dimerS lam‖ ≤
        Real.exp (dimerKappa lam * (3 + 4 * lam) * ‖z - lam‖ * j)
  | 0, _ => by
    simp only [tmProd, Matrix.mul_one,
      Matrix.nonsing_inv_mul _ (S_isUnit lam h0 h1), CharP.cast_eq_zero, mul_zero, Real.exp_zero]
    exact norm_one_le
  | j + 1, hv => by
    have ih := aligned_bound lam h0 h1 v z hz m₀ j (fun i hi => hv i (by omega))
    obtain ⟨hw, heq⟩ := hv j (by omega)
    have hU := S_isUnit lam h0 h1
    have hprod : tmProd v z m₀ (2 * (j + 1)) =
        tmC (v (m₀ + 2 * j) - z) * tmC (v (m₀ + 2 * j) - z) * tmProd v z m₀ (2 * j) := by
      rw [show 2 * (j + 1) = 2 * j + 1 + 1 by ring, tmProd, tmProd, ← Matrix.mul_assoc]
      have e1 : m₀ + ((2 * j + 1 : ℕ) : ℤ) = m₀ + 2 * j + 1 := by push_cast; ring
      have e2 : m₀ + ((2 * j : ℕ) : ℤ) = m₀ + 2 * j := by push_cast; ring
      rw [e1, e2, heq]
    have hsplit : (dimerS lam)⁻¹ * tmProd v z m₀ (2 * (j + 1)) * dimerS lam =
        ((dimerS lam)⁻¹ * (tmC (v (m₀ + 2 * j) - z) * tmC (v (m₀ + 2 * j) - z)) * dimerS lam) *
          ((dimerS lam)⁻¹ * tmProd v z m₀ (2 * j) * dimerS lam) := by
      rw [hprod]
      simp only [Matrix.mul_assoc, Matrix.mul_nonsing_inv_cancel_left _ _ hU]
    rw [hsplit]
    calc _ ≤ ‖(dimerS lam)⁻¹ * (tmC (v (m₀ + 2 * j) - z) * tmC (v (m₀ + 2 * j) - z)) *
            dimerS lam‖ * ‖(dimerS lam)⁻¹ * tmProd v z m₀ (2 * j) * dimerS lam‖ :=
          norm_mul_le _ _
      _ ≤ Real.exp (dimerKappa lam * (3 + 4 * lam) * ‖z - lam‖) *
            Real.exp (dimerKappa lam * (3 + 4 * lam) * ‖z - lam‖ * j) := by
          gcongr
          exact conj_dimer_le lam h0 h1 hw z hz
      _ = Real.exp (dimerKappa lam * (3 + 4 * lam) * ‖z - lam‖ * ↑(j + 1)) := by
          rw [← Real.exp_add]
          push_cast
          ring_nf

/-- `e^{K d j} = e^{(K/2) d (2j)}`. -/
lemma exp_even (K d : ℝ) (j : ℕ) :
    Real.exp (K * d * j) = Real.exp (K / 2 * d * ((2 * j : ℕ) : ℝ)) := by
  congr 1; push_cast; ring

/-- `e^{K d j} ≤ e^{(K/2) d (2j+1)}` for `K, d ≥ 0`. -/
lemma exp_odd_le {K d : ℝ} (hK : 0 ≤ K) (hd : 0 ≤ d) (j : ℕ) :
    Real.exp (K * d * j) ≤ Real.exp (K / 2 * d * ((2 * j + 1 : ℕ) : ℝ)) := by
  rw [Real.exp_le_exp]
  have : 0 ≤ K / 2 * d := by positivity
  push_cast
  nlinarith

/-- `e^{c L} ≤ e^{c (L+1)}` for `c ≥ 0`. -/
lemma exp_succ_le {c : ℝ} (hc : 0 ≤ c) (L : ℕ) :
    Real.exp (c * L) ≤ Real.exp (c * ((L + 1 : ℕ) : ℝ)) := by
  rw [Real.exp_le_exp]
  push_cast
  nlinarith

end DimerTransferAux

/-- **Uniform transfer-matrix bound near `E_c = λ`** (paper: Lemma `lem:dimertransfer`):
`‖T^z(n, m)‖ ≤ K₁ e^{K₂ |z - λ| (n - m)}` with `K₁ = κ (2 + 2λ)²`, `K₂ = κ (3 + 4λ)/2`, for every dimer
configuration and every `z` with `|z - λ| ≤ 1`. -/
theorem dimer_transfer_bound (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (r : ℤ) (v : ℤ → ℝ) (m : ℤ)
    (L : ℕ) (hv : IsDimerCfg lam r v m L) (z : ℂ) (hz : ‖z - lam‖ ≤ 1) :
    ‖tmProd v z m L‖ ≤
      dimerKappa lam * (2 + 2 * lam) ^ 2 *
        Real.exp (dimerKappa lam * (3 + 4 * lam) / 2 * ‖z - lam‖ * L) := by
  obtain ⟨hv1, hv2⟩ := hv
  have hκ1 : 1 ≤ dimerKappa lam := DimerTransferAux.kappa_ge_one lam h0 h1
  have hκ0 : 0 ≤ dimerKappa lam := by linarith
  have hd0 : 0 ≤ ‖z - (lam : ℂ)‖ := norm_nonneg _
  have hK0 : 0 ≤ dimerKappa lam * (3 + 4 * lam) := by positivity
  have hc1 : 1 ≤ 2 + 2 * lam := by linarith
  have hc2 : 2 + 2 * lam ≤ (2 + 2 * lam) ^ 2 := by nlinarith
  -- the bound for a product starting at the first site `2k - r` of a dimer
  have key : ∀ (k : ℤ) (L' : ℕ), m ≤ 2 * k - r → 2 * k - r + L' ≤ m + L →
      ‖tmProd v z (2 * k - r) L'‖ ≤
        dimerKappa lam * (2 + 2 * lam) *
          Real.exp (dimerKappa lam * (3 + 4 * lam) / 2 * ‖z - lam‖ * L') := by
    intro k L' hm hL
    obtain ⟨j, hj⟩ := Nat.even_or_odd' L'
    have hal : ‖(dimerS lam)⁻¹ * tmProd v z (2 * k - r) (2 * j) * dimerS lam‖ ≤
        Real.exp (dimerKappa lam * (3 + 4 * lam) * ‖z - lam‖ * j) := by
      apply DimerTransferAux.aligned_bound lam h0 h1 v z hz
      intro i hi
      have hi' : 2 * (i : ℤ) + 1 < L' := by omega
      refine ⟨hv1 _ (by omega) (by omega), ?_⟩
      have e : 2 * (k + (i : ℤ)) - r = 2 * k - r + 2 * i := by ring
      have h := hv2 (k + i) (by rw [e]; omega) (by rw [e]; omega)
      rw [e] at h
      exact h.symm
    have hP : ‖tmProd v z (2 * k - r) (2 * j)‖ ≤
        dimerKappa lam * Real.exp (dimerKappa lam * (3 + 4 * lam) * ‖z - lam‖ * j) :=
      (DimerTransferAux.norm_le_conj lam h0 h1 _).trans (by gcongr)
    rcases hj with rfl | rfl
    · calc _ ≤ dimerKappa lam * Real.exp (dimerKappa lam * (3 + 4 * lam) * ‖z - lam‖ * j) := hP
        _ = dimerKappa lam * 1 *
              Real.exp (dimerKappa lam * (3 + 4 * lam) / 2 * ‖z - lam‖ * ((2 * j : ℕ) : ℝ)) := by
            rw [DimerTransferAux.exp_even, mul_one]
        _ ≤ _ := by gcongr
    · have hy : v (2 * k - r + ((2 * j : ℕ) : ℤ)) = lam ∨
          v (2 * k - r + ((2 * j : ℕ) : ℤ)) = -lam := hv1 _ (by omega) (by omega)
      rw [tmProd]
      calc _ ≤ ‖tmC (v (2 * k - r + ((2 * j : ℕ) : ℤ)) - z)‖ *
              ‖tmProd v z (2 * k - r) (2 * j)‖ := norm_mul_le _ _
        _ ≤ (2 + 2 * lam) *
              (dimerKappa lam * Real.exp (dimerKappa lam * (3 + 4 * lam) * ‖z - lam‖ * j)) := by
            gcongr
            exact DimerTransferAux.norm_tm_single_le lam h0 hy z hz
        _ = dimerKappa lam * (2 + 2 * lam) *
              Real.exp (dimerKappa lam * (3 + 4 * lam) * ‖z - lam‖ * j) := by ring
        _ ≤ _ := mul_le_mul_of_nonneg_left (DimerTransferAux.exp_odd_le hK0 hd0 j)
            (by positivity)
  -- assemble, according to the parity of `m + r`
  obtain ⟨k, hk⟩ := Int.even_or_odd' (m + r)
  rcases hk with hk | hk
  · -- `m = 2k - r` is the first site of a dimer
    have hm : m = 2 * k - r := by omega
    have h := key k L (by omega) (by omega)
    rw [← hm] at h
    calc _ ≤ _ := h
      _ ≤ _ := by gcongr
  · -- `m + 1 = 2(k+1) - r` is the first site of a dimer
    cases L with
    | zero =>
      simp only [tmProd, CharP.cast_eq_zero, mul_zero, Real.exp_zero, mul_one]
      calc ‖(1 : Matrix (Fin 2) (Fin 2) ℂ)‖ ≤ 1 := DimerTransferAux.norm_one_le
        _ ≤ dimerKappa lam * (2 + 2 * lam) ^ 2 := by nlinarith
    | succ L'' =>
      have hm : m + 1 = 2 * (k + 1) - r := by omega
      have hsplit : tmProd v z m (L'' + 1) = tmProd v z (m + 1) L'' * tmProd v z m 1 := by
        rw [add_comm L'' 1, DimerTransferAux.tmProd_add]
        push_cast
        rfl
      have hP := key (k + 1) L'' (by omega) (by push_cast; omega)
      rw [← hm] at hP
      have hy : v m = lam ∨ v m = -lam := hv1 m le_rfl (by push_cast; omega)
      have h1f : ‖tmProd v z m 1‖ ≤ 2 + 2 * lam := by
        simp only [tmProd, Matrix.mul_one, CharP.cast_eq_zero, add_zero]
        exact DimerTransferAux.norm_tm_single_le lam h0 hy z hz
      rw [hsplit]
      calc _ ≤ ‖tmProd v z (m + 1) L''‖ * ‖tmProd v z m 1‖ := norm_mul_le _ _
        _ ≤ (dimerKappa lam * (2 + 2 * lam) *
              Real.exp (dimerKappa lam * (3 + 4 * lam) / 2 * ‖z - lam‖ * L'')) * (2 + 2 * lam) := by
            gcongr
        _ = dimerKappa lam * (2 + 2 * lam) ^ 2 *
              Real.exp (dimerKappa lam * (3 + 4 * lam) / 2 * ‖z - lam‖ * L'') := by ring
        _ ≤ _ := mul_le_mul_of_nonneg_left
            (DimerTransferAux.exp_succ_le (c := dimerKappa lam * (3 + 4 * lam) / 2 * ‖z - lam‖)
              (by positivity) L'') (by positivity)

end DFL
