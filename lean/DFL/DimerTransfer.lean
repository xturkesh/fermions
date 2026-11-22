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

/-- **Uniform transfer-matrix bound near `E_c = λ`** (paper: Lemma `lem:dimertransfer`):
`‖T^z(n, m)‖ ≤ K₁ e^{K₂ |z - λ| (n - m)}` with `K₁ = κ (2 + 2λ)²`, `K₂ = κ (3 + 4λ)/2`, for every dimer
configuration and every `z` with `|z - λ| ≤ 1`. -/
theorem dimer_transfer_bound (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (r : ℤ) (v : ℤ → ℝ) (m : ℤ)
    (L : ℕ) (hv : IsDimerCfg lam r v m L) (z : ℂ) (hz : ‖z - lam‖ ≤ 1) :
    ‖tmProd v z m L‖ ≤
      dimerKappa lam * (2 + 2 * lam) ^ 2 *
        Real.exp (dimerKappa lam * (3 + 4 * lam) / 2 * ‖z - lam‖ * L) := by
  sorry

end DFL
