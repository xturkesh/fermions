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
  congr 3 <;> ring1

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
  sorry

/-- **The condition number of `S`** (paper: Lemma `lem:dimercritical`(iii)):
`κ = ‖S‖ ‖S⁻¹‖ = ((1+λ)/(1-λ))^{1/2}`. -/
theorem dimerS_kappa (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    ‖dimerS lam‖ * ‖(dimerS lam)⁻¹‖ = dimerKappa lam := by
  sorry

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
    congr 3 <;> ring1

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
