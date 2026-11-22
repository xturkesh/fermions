import DFL.ECSup
import DFL.InfVolNumber
import DFL.Normalization

/-!
# Annealed localization of the full correlator on `ℤ` (paper: Corollary `cor:ECinfvol`, `X = ℤ`)

On `ℓ²(ℤ)` let `A` be a bounded operator with the matrix elements of `A_ℤ = Δ_ℤ + λ diag(q)`, `λ = -2h/J` (the normalization
of `DFL.LDTEntry`: `⟨δ_x, A δ_y⟩ = 1` for `|x - y| = 1`, `λ sgn(q_x)` for `x = y`, `0` otherwise; `IsAZ`; such an operator
exists, `exists_isAZ`). The eigenfunction correlator is
`Q^ℤ(x, y) = sup { |⟨δ_x, f(A) δ_y⟩| : |f| ≤ 1 }` (`corrQZ`), with `f(A)` the continuous functional calculus of the
self-adjoint `A` (the supremum over Borel `f` coincides with the supremum over continuous `f`, paper: proof of
Corollary `cor:ECinfvol`). Under the large-deviation input, `E Q^ℤ(x, y) ≤ C e^{-η |x - y|}` (`lintegral_corrQZ_le`).
-/

open MeasureTheory Filter Topology
open scoped InnerProductSpace

noncomputable section

namespace DFL

/-- `ℓ²(ℤ)`. -/
abbrev L2Z : Type := lp (fun _ : ℤ => ℂ) 2

/-- The basis vector `δ_x` of `ℓ²(ℤ)`. -/
def deltaZ (x : ℤ) : L2Z := lp.single 2 x 1

/-- `A` has the matrix elements of `A_ℤ = Δ_ℤ + λ diag(q)` (`q_x` read as `sgn(q_x) ∈ {±1}`). -/
def IsAZ (lam : ℝ) (q : CfgZ) (A : L2Z →L[ℂ] L2Z) : Prop :=
  ∀ x y : ℤ, ⟪deltaZ x, A (deltaZ y)⟫_ℂ =
    ((if x = y + 1 ∨ y = x + 1 then 1 else 0) + (if x = y then lam * sgn (q x) else 0) : ℝ)

namespace ECInfVolQAux

/-- The coordinates of a vector of `ℓ²(ℤ)` are its inner products with the `δ_x`. -/
lemma inner_deltaZ_left (x : ℤ) (v : L2Z) : ⟪deltaZ x, v⟫_ℂ = v x := by
  rw [deltaZ, lp.inner_single_left]
  simp

/-- The coordinates of `δ_y`. -/
lemma deltaZ_apply (y z : ℤ) : deltaZ y z = if z = y then 1 else 0 := by
  simp [deltaZ, Pi.single_apply]

/-- A weighted relabelling `z ↦ c z * ψ (e z)` (`|c| ≤ 1`) of a vector of `ℓ²(ℤ)` is in `ℓ²(ℤ)`. -/
lemma memℓp_wt (e : ℤ ≃ ℤ) (c : ℤ → ℂ) (hc : ∀ z, ‖c z‖ ≤ 1) (ψ : L2Z) :
    Memℓp (fun z => c z * ψ (e z)) 2 := by
  have hp : 0 < (2 : ENNReal).toReal := by norm_num
  have hs : Summable (fun z => ‖ψ (e z)‖ ^ (2 : ENNReal).toReal) :=
    (e.summable_iff (f := fun z => ‖ψ z‖ ^ (2 : ENNReal).toReal)).mpr ((lp.memℓp ψ).summable hp)
  refine memℓp_gen (hs.of_nonneg_of_le (fun z => by positivity) fun z => ?_)
  refine Real.rpow_le_rpow (norm_nonneg _) ?_ hp.le
  rw [norm_mul]
  exact mul_le_of_le_one_left (norm_nonneg _) (hc z)

/-- The weighted relabelling `(Tψ)(z) = c z ψ(e z)` as a bounded operator on `ℓ²(ℤ)`. -/
def wOp (e : ℤ ≃ ℤ) (c : ℤ → ℂ) (hc : ∀ z, ‖c z‖ ≤ 1) : L2Z →L[ℂ] L2Z :=
  LinearMap.mkContinuous
    { toFun := fun ψ => ⟨fun z => c z * ψ (e z), memℓp_wt e c hc ψ⟩
      map_add' := fun ψ φ => by
        ext z
        show c z * (ψ + φ) (e z) = c z * ψ (e z) + c z * φ (e z)
        rw [lp.coeFn_add, Pi.add_apply, mul_add]
      map_smul' := fun a ψ => by
        ext z
        simp only [lp.coeFn_smul, Pi.smul_apply, smul_eq_mul, RingHom.id_apply]
        ring }
    1 (fun ψ => by
      have hp : 0 < (2 : ENNReal).toReal := by norm_num
      rw [one_mul]
      refine lp.norm_le_of_tsum_le hp (norm_nonneg _) ?_
      have hs : Summable (fun z => ‖ψ (e z)‖ ^ (2 : ENNReal).toReal) :=
        (e.summable_iff (f := fun z => ‖ψ z‖ ^ (2 : ENNReal).toReal)).mpr
          ((lp.memℓp ψ).summable hp)
      calc ∑' z, ‖c z * ψ (e z)‖ ^ (2 : ENNReal).toReal
          ≤ ∑' z, ‖ψ (e z)‖ ^ (2 : ENNReal).toReal := by
            refine Summable.tsum_le_tsum (fun z => ?_) ((memℓp_wt e c hc ψ).summable hp) hs
            refine Real.rpow_le_rpow (norm_nonneg _) ?_ hp.le
            rw [norm_mul]
            exact mul_le_of_le_one_left (norm_nonneg _) (hc z)
        _ = ∑' z, ‖ψ z‖ ^ (2 : ENNReal).toReal := e.tsum_eq (fun z => ‖ψ z‖ ^ (2 : ENNReal).toReal)
        _ = ‖ψ‖ ^ (2 : ENNReal).toReal := (lp.norm_rpow_eq_tsum hp ψ).symm)

/-- `(Tψ)(z) = c z ψ(e z)`. -/
lemma wOp_apply (e : ℤ ≃ ℤ) (c : ℤ → ℂ) (hc : ∀ z, ‖c z‖ ≤ 1) (ψ : L2Z) (z : ℤ) :
    wOp e c hc ψ z = c z * ψ (e z) := rfl

/-- `|sgn a| ≤ 1`, as a complex number. -/
lemma norm_sgn_le (a : ZMod 2) : ‖((sgn a : ℝ) : ℂ)‖ ≤ 1 := by
  rw [Complex.norm_real, Real.norm_eq_abs, InfVolAux.abs_sgn]

/-- The operator `Δ_ℤ + λ diag(sgn q)` on `ℓ²(ℤ)`. -/
def AZop (lam : ℝ) (q : CfgZ) : L2Z →L[ℂ] L2Z :=
  wOp (Equiv.addRight 1) (fun _ => 1) (fun _ => by simp) +
    wOp (Equiv.subRight 1) (fun _ => 1) (fun _ => by simp) +
    (lam : ℂ) • wOp (Equiv.refl ℤ) (fun z => ((sgn (q z) : ℝ) : ℂ)) (fun z => norm_sgn_le (q z))

/-- `(Aψ)(z) = ψ(z+1) + ψ(z-1) + λ sgn(q_z) ψ(z)`. -/
lemma AZop_apply (lam : ℝ) (q : CfgZ) (ψ : L2Z) (z : ℤ) :
    AZop lam q ψ z = ψ (z + 1) + ψ (z - 1) + (lam : ℂ) * ((sgn (q z) : ℝ) : ℂ) * ψ z := by
  simp only [AZop, add_apply, smul_apply, lp.coeFn_add, lp.coeFn_smul, Pi.add_apply,
    Pi.smul_apply, wOp_apply, Equiv.coe_addRight, Equiv.subRight_apply, Equiv.refl_apply,
    smul_eq_mul, one_mul]
  ring

/-- `AZop` has the matrix elements of `Δ_ℤ + λ diag(sgn q)`. -/
lemma isAZ_AZop (lam : ℝ) (q : CfgZ) : IsAZ lam q (AZop lam q) := by
  intro x y
  rw [inner_deltaZ_left, AZop_apply, deltaZ_apply, deltaZ_apply, deltaZ_apply]
  by_cases h1 : x + 1 = y
  · have h2 : ¬ x - 1 = y := by omega
    have h3 : ¬ x = y := by omega
    simp [h1, h2, h3]
  · by_cases h2 : x - 1 = y
    · have h3 : ¬ x = y := by omega
      have h4 : x = y + 1 ∨ y = x + 1 := Or.inl (by omega)
      simp [h1, h2, h3, h4]
    · have h4 : ¬ (x = y + 1 ∨ y = x + 1) := by omega
      by_cases h3 : x = y
      · simp [h3]
      · simp [h1, h2, h3, h4]

end ECInfVolQAux

/-- **The operator exists** (a bounded tridiagonal operator, `‖A‖ ≤ 2 + |λ|`). -/
theorem exists_isAZ (lam : ℝ) (q : CfgZ) : ∃ A : L2Z →L[ℂ] L2Z, IsAZ lam q A :=
  ⟨ECInfVolQAux.AZop lam q, ECInfVolQAux.isAZ_AZop lam q⟩

/-- The eigenfunction correlator `Q^ℤ(x, y) = sup_{|f| ≤ 1} |⟨δ_x, f(A) δ_y⟩|` (`f(A) = cfc f A`). -/
def corrQZ (A : L2Z →L[ℂ] L2Z) (x y : ℤ) : ℝ :=
  ⨆ (f : ℂ → ℂ) (_ : ∀ z, ‖f z‖ ≤ 1), ‖⟪deltaZ x, cfc f A (deltaZ y)⟫_ℂ‖

namespace ECInfVolQAux

/-- A bounded operator on `ℓ²(ℤ)` is determined by its matrix elements. -/
lemma ext_deltaZ {T S : L2Z →L[ℂ] L2Z}
    (h : ∀ x y, ⟪deltaZ x, T (deltaZ y)⟫_ℂ = ⟪deltaZ x, S (deltaZ y)⟫_ℂ) : T = S := by
  have hv : ∀ y, T (deltaZ y) = S (deltaZ y) := fun y => by
    ext x
    rw [← inner_deltaZ_left, ← inner_deltaZ_left]
    exact h x y
  refine lp.ext_continuousLinearMap (by norm_num : (2 : ENNReal) ≠ ⊤) fun i => ?_
  refine ContinuousLinearMap.ext fun a => ?_
  simp only [ContinuousLinearMap.comp_apply, lp.singleContinuousLinearMap_apply]
  have e : lp.single 2 i a = a • deltaZ i := by
    rw [deltaZ, ← lp.single_smul, smul_eq_mul, mul_one]
  rw [e, map_smul, map_smul, hv i]

/-- The matrix elements of `Δ_ℤ + λ diag(sgn q)`, as those of `h^ℤ` with `J = -1`, `2h = λ`. -/
lemma isAZ_iff_hZ (lam : ℝ) (q : CfgZ) (x y : ℤ) :
    ((if x = y + 1 ∨ y = x + 1 then 1 else 0) + (if x = y then lam * sgn (q x) else 0) : ℝ) =
      hZ (-1) (lam / 2) q x y := by
  simp only [hZ, neg_neg]
  congr 1
  split_ifs <;> ring

/-- `h^ℤ` is symmetric. -/
lemma hZ_comm (J h : ℝ) (q : CfgZ) (x y : ℤ) : hZ J h q x y = hZ J h q y x := by
  unfold hZ
  by_cases hxy : x = y
  · subst hxy; rfl
  · have : ¬ y = x := fun e => hxy e.symm
    rw [ite_eq_right hxy, ite_eq_right this]
    congr 1
    simp only [or_comm]

/-- An operator with the matrix elements of `A_ℤ` is self-adjoint. -/
lemma isSelfAdjoint_of_isAZ {lam : ℝ} {q : CfgZ} {A : L2Z →L[ℂ] L2Z} (hA : IsAZ lam q A) :
    IsSelfAdjoint A := by
  rw [ContinuousLinearMap.isSelfAdjoint_iff']
  refine ext_deltaZ fun x y => ?_
  rw [ContinuousLinearMap.adjoint_inner_right, ← inner_conj_symm, hA y x, hA x y,
    Complex.conj_ofReal, isAZ_iff_hZ, isAZ_iff_hZ, hZ_comm]

/-- `A δ_x = ∑_{|z - x| ≤ 1} A_{zx} δ_z`. -/
lemma apply_deltaZ {lam : ℝ} {q : CfgZ} {A : L2Z →L[ℂ] L2Z} (hA : IsAZ lam q A) (x : ℤ) :
    A (deltaZ x) =
      ∑ z ∈ Finset.Icc (x - 1) (x + 1), ((hZ (-1) (lam / 2) q z x : ℝ) : ℂ) • deltaZ z := by
  ext w
  rw [← inner_deltaZ_left, hA w x, isAZ_iff_hZ, lp.coeFn_sum, Finset.sum_apply]
  simp only [lp.coeFn_smul, Pi.smul_apply, deltaZ_apply, smul_eq_mul, mul_ite, mul_one, mul_zero]
  rw [Finset.sum_ite_eq]
  split_ifs with hw
  · rfl
  · rw [InfVolAux.hZ_eq_zero (-1) (lam / 2) q w x ?_, Complex.ofReal_zero]
    simp only [Finset.mem_Icc] at hw ⊢
    omega

/-- **Matrix elements of the powers**: `⟨δ_x, A^k δ_y⟩ = (A_ℤ^k)_{xy}` (a finite sum over paths). -/
lemma inner_pow_deltaZ {lam : ℝ} {q : CfgZ} {A : L2Z →L[ℂ] L2Z} (hA : IsAZ lam q A) (k : ℕ) :
    ∀ x y : ℤ, ⟪deltaZ x, (A ^ k) (deltaZ y)⟫_ℂ = ((hZpow (-1) (lam / 2) q k x y : ℝ) : ℂ) := by
  have hsa := ContinuousLinearMap.isSelfAdjoint_iff'.mp (isSelfAdjoint_of_isAZ hA)
  induction k with
  | zero =>
    intro x y
    rw [pow_zero, one_apply_eq_self, inner_deltaZ_left, deltaZ_apply]
    simp only [hZpow]
    split_ifs <;> simp
  | succ k ih =>
    intro x y
    rw [pow_succ', mul_apply_eq_comp, ← ContinuousLinearMap.adjoint_inner_left, hsa,
      apply_deltaZ hA x, sum_inner]
    simp only [inner_smul_left, ih, Complex.conj_ofReal, hZpow]
    push_cast
    refine Finset.sum_congr rfl fun z _ => ?_
    rw [hZ_comm]

/-! ### Finite volume -/

/-- Spectral form of the powers of a real matrix with a complete orthonormal family of eigenvectors:
`(B^k)_{ij} = ∑_α μ_α^k φ_α(i) φ_α(j)`. -/
lemma pow_apply_eq_sum {n : ℕ} (B : Matrix (Fin n) (Fin n) ℝ) (φ : Fin n → Fin n → ℝ)
    (μ : Fin n → ℝ) (hc : ∀ i j, ∑ α, φ α i * φ α j = if i = j then 1 else 0)
    (he : ∀ α i, ∑ j, B i j * φ α j = μ α * φ α i) (k : ℕ) (i j : Fin n) :
    (B ^ k) i j = ∑ α, μ α ^ k * (φ α i * φ α j) := by
  induction k generalizing i with
  | zero =>
    rw [pow_zero, Matrix.one_apply, ← hc i j]
    simp
  | succ k ih =>
    rw [pow_succ', Matrix.mul_apply]
    simp_rw [ih, Finset.mul_sum]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun α _ => ?_
    calc ∑ l, B i l * (μ α ^ k * (φ α l * φ α j))
        = μ α ^ k * φ α j * ∑ l, B i l * φ α l := by
          rw [Finset.mul_sum]
          exact Finset.sum_congr rfl fun l _ => by ring
      _ = μ α ^ (k + 1) * (φ α i * φ α j) := by rw [he α i]; ring

/-- An eigenvalue of a matrix whose rows have `ℓ¹`-norm at most `c` has modulus at most `c`. -/
lemma abs_eig_le {n : ℕ} (B : Matrix (Fin n) (Fin n) ℝ) (c : ℝ) (hB : ∀ i, ∑ l, |B i l| ≤ c)
    (v : Fin n → ℝ) (μ : ℝ) (hv : ∃ i, v i ≠ 0) (he : ∀ i, ∑ j, B i j * v j = μ * v i) :
    |μ| ≤ c := by
  obtain ⟨i1, hi1⟩ := hv
  have : Nonempty (Fin n) := ⟨i1⟩
  obtain ⟨i0, -, hmax⟩ := Finset.exists_max_image Finset.univ (fun i => |v i|) Finset.univ_nonempty
  have hpos : 0 < |v i0| := lt_of_lt_of_le (abs_pos.mpr hi1) (hmax i1 (Finset.mem_univ _))
  have h1 : |μ| * |v i0| ≤ c * |v i0| := by
    calc |μ| * |v i0| = |∑ j, B i0 j * v j| := by rw [he i0, abs_mul]
      _ ≤ ∑ j, |B i0 j * v j| := Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ j, |B i0 j| * |v i0| := Finset.sum_le_sum fun j _ => by
          rw [abs_mul]
          exact mul_le_mul_of_nonneg_left (hmax j (Finset.mem_univ _)) (abs_nonneg _)
      _ = (∑ j, |B i0 j|) * |v i0| := by rw [Finset.sum_mul]
      _ ≤ c * |v i0| := mul_le_mul_of_nonneg_right (hB i0) hpos.le
  exact le_of_mul_le_mul_right h1 hpos

/-- `h_q = -J (Δ + λ diag(q))` with `λ = -2h/J`: `h_q(J, h) = (-J) • h_q(-1, λ/2)`. -/
lemma hq_eq_smul {n : ℕ} (J h : ℝ) (hJ : J ≠ 0) (q : Cfg n) :
    hq J h q = (-J) • hq (-1) (-(2 * h / J) / 2) q := by
  unfold hq
  rw [smul_add, smul_smul, smul_smul]
  congr 2
  · ring
  · field_simp

/-- **Finite volume**: `|p(A_Λ)_{ij}| ≤ c Q^Λ(i, j)` if `|p| ≤ c` on `[-2-|λ|, 2+|λ|]`, where
`A_Λ = Δ_Λ + λ diag(q)` and `Q^Λ` is the correlator of `h_q = -J A_Λ`. -/
lemma norm_poly_apply_le (J h : ℝ) (hJ : J ≠ 0) {n : ℕ} (q' : Cfg n) (P : Polynomial ℂ) (c : ℝ)
    (hP : ∀ t : ℝ, |t| ≤ 2 + |(-(2 * h / J))| → ‖P.eval (t : ℂ)‖ ≤ c) (i j : Fin n) :
    ‖∑ k ∈ Finset.range (P.natDegree + 1), P.coeff k *
      (((hq (-1) (-(2 * h / J) / 2) q' ^ k) i j : ℝ) : ℂ)‖ ≤
      c * corrQ (hq J h q') (hq_isHermitian J h q') i j := by
  set B := hq (-1) (-(2 * h / J) / 2) q' with hB
  set hA := hq_isHermitian J h q' with hA_def
  set μ : Fin n → ℝ := fun α => -hA.eigenvalues α / J with hμ
  have hent : ∀ i j, hq J h q' i j = -J * B i j := fun i j => by
    rw [hq_eq_smul J h hJ q']
    rfl
  have heig : ∀ α i, ∑ j, B i j * hA.eigenvectorBasis α j = μ α * hA.eigenvectorBasis α i := by
    intro α i
    calc ∑ j, B i j * hA.eigenvectorBasis α j
        = (-1 / J) * ∑ j, hq J h q' i j * hA.eigenvectorBasis α j := by
          rw [Finset.mul_sum]
          refine Finset.sum_congr rfl fun j _ => ?_
          rw [hent]
          field_simp
      _ = μ α * hA.eigenvectorBasis α i := by
          rw [NormalizationAux.eig_eq _ hA α i, hμ]
          ring
  have hpow := pow_apply_eq_sum B (fun α => hA.eigenvectorBasis α) μ
    (NormalizationAux.eig_compl _ hA) heig
  have hμle : ∀ α, |μ α| ≤ 2 + |(-(2 * h / J))| := by
    intro α
    refine abs_eig_le B _ (fun i => ?_) _ _ ?_ (heig α)
    · have := InfVolAux.sum_abs_hq_le (-1) (-(2 * h / J) / 2) q' i
      rw [abs_div, abs_neg, abs_one, abs_two] at this
      linarith
    · by_contra hcon
      have h0 : ∀ i, hA.eigenvectorBasis α i = 0 := fun i => by
        by_contra hi
        exact hcon ⟨i, hi⟩
      have := NormalizationAux.eig_orth _ hA α α
      simp [h0] at this
  have hsum : ∑ k ∈ Finset.range (P.natDegree + 1), P.coeff k * (((B ^ k) i j : ℝ) : ℂ) =
      ∑ α, P.eval (μ α : ℂ) * ((hA.eigenvectorBasis α i * hA.eigenvectorBasis α j : ℝ) : ℂ) := by
    simp_rw [hpow]
    push_cast
    simp_rw [Polynomial.eval_eq_sum_range, Finset.sum_mul, Finset.mul_sum]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun α _ => Finset.sum_congr rfl fun k _ => ?_
    ring
  rw [hsum, corrQ, Finset.mul_sum]
  refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun α _ => ?_)
  rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_mul]
  exact mul_le_mul_of_nonneg_right (hP _ (hμle α)) (by positivity)

/-! ### Approximation -/

/-- **Weierstrass** for complex-valued continuous functions on a real interval. -/
lemma exists_poly_near (F : ℂ → ℂ) (hF : Continuous F) (R ε : ℝ) (hε : 0 < ε) :
    ∃ P : Polynomial ℂ, ∀ t ∈ Set.Icc (-R) R, ‖F t - P.eval (t : ℂ)‖ < ε := by
  have hre : Continuous fun t : ℝ => (F t).re :=
    Complex.continuous_re.comp (hF.comp Complex.continuous_ofReal)
  have him : Continuous fun t : ℝ => (F t).im :=
    Complex.continuous_im.comp (hF.comp Complex.continuous_ofReal)
  obtain ⟨p1, hp1⟩ := exists_polynomial_near_of_continuousOn (-R) R _ hre.continuousOn (ε / 2)
    (by positivity)
  obtain ⟨p2, hp2⟩ := exists_polynomial_near_of_continuousOn (-R) R _ him.continuousOn (ε / 2)
    (by positivity)
  refine ⟨p1.map Complex.ofRealHom + Polynomial.C Complex.I * p2.map Complex.ofRealHom,
    fun t ht => ?_⟩
  have e : ∀ p : Polynomial ℝ, (p.map Complex.ofRealHom).eval (t : ℂ) = ((p.eval t : ℝ) : ℂ) :=
    fun p => by
      rw [Polynomial.eval_map]
      exact Polynomial.eval₂_hom Complex.ofRealHom t
  rw [Polynomial.eval_add, Polynomial.eval_mul, Polynomial.eval_C, e, e]
  have h1 := hp1 t ht
  have h2 := hp2 t ht
  have hsplit : F t - (((p1.eval t : ℝ) : ℂ) + Complex.I * ((p2.eval t : ℝ) : ℂ)) =
      (((F t).re - p1.eval t : ℝ) : ℂ) + (((F t).im - p2.eval t : ℝ) : ℂ) * Complex.I := by
    apply Complex.ext <;> simp
  rw [hsplit]
  calc ‖(((F t).re - p1.eval t : ℝ) : ℂ) + (((F t).im - p2.eval t : ℝ) : ℂ) * Complex.I‖
      ≤ ‖(((F t).re - p1.eval t : ℝ) : ℂ)‖ + ‖(((F t).im - p2.eval t : ℝ) : ℂ) * Complex.I‖ :=
        norm_add_le _ _
    _ = |(F t).re - p1.eval t| + |(F t).im - p2.eval t| := by
        rw [norm_mul, Complex.norm_I, mul_one, Complex.norm_real, Complex.norm_real,
          Real.norm_eq_abs, Real.norm_eq_abs]
    _ < ε / 2 + ε / 2 :=
        add_lt_add (by rw [abs_sub_comm]; exact h1) (by rw [abs_sub_comm]; exact h2)
    _ = ε := by ring

/-- The spectrum of a self-adjoint operator is real and bounded by `‖A‖ ‖1‖`. -/
lemma mem_spectrum_re {A : L2Z →L[ℂ] L2Z} (hsa : IsSelfAdjoint A) {z : ℂ}
    (hz : z ∈ spectrum ℂ A) :
    z = (z.re : ℂ) ∧ |z.re| ≤ ‖A‖ * ‖(1 : L2Z →L[ℂ] L2Z)‖ :=
  ⟨hsa.mem_spectrum_eq_re hz,
    (Complex.abs_re_le_norm z).trans (spectrum.norm_le_norm_mul_of_mem hz)⟩

/-- `‖δ_x‖ = 1`. -/
lemma norm_deltaZ (x : ℤ) : ‖deltaZ x‖ = 1 := by
  rw [deltaZ, lp.norm_single (by norm_num)]
  simp

/-- Matrix elements are bounded by the operator norm. -/
lemma norm_inner_deltaZ_le (T : L2Z →L[ℂ] L2Z) (x y : ℤ) : ‖⟪deltaZ x, T (deltaZ y)⟫_ℂ‖ ≤ ‖T‖ := by
  calc ‖⟪deltaZ x, T (deltaZ y)⟫_ℂ‖ ≤ ‖deltaZ x‖ * ‖T (deltaZ y)‖ := norm_inner_le_norm _ _
    _ ≤ ‖deltaZ x‖ * (‖T‖ * ‖deltaZ y‖) :=
        mul_le_mul_of_nonneg_left (T.le_opNorm _) (norm_nonneg _)
    _ = ‖T‖ := by rw [norm_deltaZ, norm_deltaZ]; ring

/-! ### Comparison with finite volume -/

/-- `Q^{Λ_M}(x, y)` as a function of the charges in `Λ_M = [-M, M]` (`0` unless `x, y ∈ Λ_M`). -/
def corrQR (J h : ℝ) (M : ℕ) (x y : ℤ) (q' : Cfg (2 * M + 1)) : ℝ :=
  if hxy : |x| ≤ M ∧ |y| ≤ M then
    corrQ (hq J h q') (hq_isHermitian J h q') (idxZ M x hxy.1) (idxZ M y hxy.2) else 0

/-- `0 ≤ Q^{Λ_M}`. -/
lemma corrQR_nonneg (J h : ℝ) (M : ℕ) (x y : ℤ) (q' : Cfg (2 * M + 1)) :
    0 ≤ corrQR J h M x y q' := by
  unfold corrQR
  split_ifs
  · exact corrQ_nonneg _ _ _ _
  · exact le_refl 0

/-- `Q^{Λ_M} ≤ 1`. -/
lemma corrQR_le_one (J h : ℝ) (M : ℕ) (x y : ℤ) (q' : Cfg (2 * M + 1)) :
    corrQR J h M x y q' ≤ 1 := by
  unfold corrQR
  split_ifs
  · exact ECProofAux.corrQ_le_one _ _ _ _
  · exact zero_le_one

/-- **Comparison with finite volume** (paper: proof of Corollary `cor:ECinfvol`): for `|f| ≤ 1`
and `ε > 0`, `|⟨δ_x, f(A) δ_y⟩| ≤ 2ε + Q^{Λ_M}(x, y)` for all large `M`. -/
lemma norm_inner_cfc_le (J h : ℝ) (hJ : J ≠ 0) (q : CfgZ) (A : L2Z →L[ℂ] L2Z)
    (hA : IsAZ (-(2 * h / J)) q A) (x y : ℤ) (f : ℂ → ℂ) (hf : ∀ z, ‖f z‖ ≤ 1) (ε : ℝ)
    (hε : 0 < ε) : ∃ M0 : ℕ, ∀ M : ℕ, M0 ≤ M →
      ‖⟪deltaZ x, cfc f A (deltaZ y)⟫_ℂ‖ ≤ 2 * ε + corrQR J h M x y (restrictZ M q) := by
  have hsa := isSelfAdjoint_of_isAZ hA
  have hnormal : IsStarNormal A := hsa.isStarNormal
  by_cases hcont : ContinuousOn f (spectrum ℂ A)
  swap
  · refine ⟨0, fun M _ => ?_⟩
    rw [cfc_apply_of_not_continuousOn A hcont, zero_apply, inner_zero_right, norm_zero]
    have := corrQR_nonneg J h M x y (restrictZ M q)
    linarith
  -- extend `f|σ(A)` to a continuous function on `ℂ` with values in the closed unit disc
  obtain ⟨F, hF1, hF2⟩ := ContinuousMap.exists_forall_mem_restrict_eq (spectrum.isClosed A)
    (⟨fun z => f z, hcont.domRestrict⟩ : C(spectrum ℂ A, ℂ))
    (t := Metric.closedBall (0 : ℂ) 1) (fun z => by
      rw [Metric.mem_closedBall, dist_zero_right]
      exact hf z)
  have hFn : ∀ z, ‖F z‖ ≤ 1 := fun z => by simpa using hF1 z
  have hfF : cfc f A = cfc (⇑F) A := cfc_congr fun z hz => by
    have := congrArg (fun g : C(spectrum ℂ A, ℂ) => g ⟨z, hz⟩) hF2
    simp only [ContinuousMap.restrict_apply] at this
    exact this.symm
  set R := max (‖A‖ * ‖(1 : L2Z →L[ℂ] L2Z)‖) (2 + |(-(2 * h / J))|) with hR
  obtain ⟨P, hP⟩ := exists_poly_near F F.continuous R ε hε
  refine ⟨x.natAbs + P.natDegree + y.natAbs, fun M hM => ?_⟩
  have hx : |x| ≤ M := by rw [Int.abs_eq_natAbs]; omega
  have hy : |y| ≤ M := by rw [Int.abs_eq_natAbs]; omega
  -- the operator part: `‖f(A) - p(A)‖ ≤ ε`
  have h1 : ‖cfc (⇑F) A - Polynomial.aeval A P‖ ≤ ε := by
    rw [← cfc_polynomial P A, ← cfc_sub (⇑F) (fun z => P.eval z) A]
    refine norm_cfc_le hε.le fun z hz => ?_
    obtain ⟨hz1, hz2⟩ := mem_spectrum_re hsa hz
    have hz3 : z.re ∈ Set.Icc (-R) R := by
      have := le_max_left (‖A‖ * ‖(1 : L2Z →L[ℂ] L2Z)‖) (2 + |(-(2 * h / J))|)
      constructor <;> linarith [abs_le.mp hz2]
    have := hP z.re hz3
    rw [← hz1] at this
    exact this.le
  -- the matrix part: `|p(A)_{xy}| ≤ (1 + ε) Q^{Λ_M}(x, y)`
  have h2 : ‖⟪deltaZ x, Polynomial.aeval A P (deltaZ y)⟫_ℂ‖ ≤
      (1 + ε) * corrQR J h M x y (restrictZ M q) := by
    have e1 : ⟪deltaZ x, Polynomial.aeval A P (deltaZ y)⟫_ℂ =
        ∑ k ∈ Finset.range (P.natDegree + 1), P.coeff k *
          (((hq (-1) (-(2 * h / J) / 2) (restrictZ M q) ^ k) (idxZ M x hx) (idxZ M y hy) : ℝ) :
            ℂ) := by
      rw [Polynomial.aeval_eq_sum_range, sum_apply, inner_sum]
      refine Finset.sum_congr rfl fun k hk => ?_
      have hk' : |x| + k ≤ M := by
        rw [Finset.mem_range] at hk
        rw [Int.abs_eq_natAbs]
        omega
      rw [smul_apply, inner_smul_right, inner_pow_deltaZ hA, InfVolAux.idxZ_eq_toFin,
        InfVolAux.idxZ_eq_toFin, InfVolAux.pow_hq_apply (-1) (-(2 * h / J) / 2) q M y hy k x hk']
    rw [e1]
    unfold corrQR
    rw [dite_eq_left ⟨hx, hy⟩]
    refine norm_poly_apply_le J h hJ _ P (1 + ε) (fun t ht => ?_) _ _
    have ht' : t ∈ Set.Icc (-R) R := by
      have := le_max_right (‖A‖ * ‖(1 : L2Z →L[ℂ] L2Z)‖) (2 + |(-(2 * h / J))|)
      constructor <;> linarith [abs_le.mp ht]
    have := hP t ht'
    calc ‖P.eval (t : ℂ)‖ = ‖F t - (F t - P.eval (t : ℂ))‖ := by rw [sub_sub_cancel]
      _ ≤ ‖F t‖ + ‖F t - P.eval (t : ℂ)‖ := norm_sub_le _ _
      _ ≤ 1 + ε := add_le_add (hFn _) this.le
  have hQ1 := corrQR_le_one J h M x y (restrictZ M q)
  have hQ0 := corrQR_nonneg J h M x y (restrictZ M q)
  calc ‖⟪deltaZ x, cfc f A (deltaZ y)⟫_ℂ‖
      = ‖⟪deltaZ x, (cfc (⇑F) A - Polynomial.aeval A P) (deltaZ y)⟫_ℂ +
          ⟪deltaZ x, Polynomial.aeval A P (deltaZ y)⟫_ℂ‖ := by
        rw [hfF, sub_apply, inner_sub_right, sub_add_cancel]
    _ ≤ ‖⟪deltaZ x, (cfc (⇑F) A - Polynomial.aeval A P) (deltaZ y)⟫_ℂ‖ +
          ‖⟪deltaZ x, Polynomial.aeval A P (deltaZ y)⟫_ℂ‖ := norm_add_le _ _
    _ ≤ ε + (1 + ε) * corrQR J h M x y (restrictZ M q) :=
        add_le_add ((norm_inner_deltaZ_le _ x y).trans h1) h2
    _ ≤ 2 * ε + corrQR J h M x y (restrictZ M q) := by nlinarith

/-- **Pointwise bound** (paper: Eq. `eq:ECinf`): `Q^ℤ(x, y) ≤ liminf_M Q^{Λ_M}(x, y)`. -/
lemma ofReal_corrQZ_le_liminf (J h : ℝ) (hJ : J ≠ 0) (q : CfgZ) (A : L2Z →L[ℂ] L2Z)
    (hA : IsAZ (-(2 * h / J)) q A) (x y : ℤ) :
    ENNReal.ofReal (corrQZ A x y) ≤
      liminf (fun M : ℕ => ENNReal.ofReal (corrQR J h M x y (restrictZ M q))) atTop := by
  set L := liminf (fun M : ℕ => ENNReal.ofReal (corrQR J h M x y (restrictZ M q))) atTop with hL
  have hL1 : L ≤ 1 := liminf_le_of_frequently_le' (Frequently.of_forall fun M =>
    ENNReal.ofReal_le_one.mpr (corrQR_le_one J h M x y (restrictZ M q)))
  have hLtop : L ≠ ⊤ := ne_top_of_le_ne_top ENNReal.one_ne_top hL1
  have hbound : corrQZ A x y ≤ L.toReal := by
    refine Real.iSup_le (fun f => Real.iSup_le (fun hf => ?_) ENNReal.toReal_nonneg)
      ENNReal.toReal_nonneg
    refine le_of_forall_pos_le_add fun ε hε => ?_
    obtain ⟨M0, hM0⟩ := norm_inner_cfc_le J h hJ q A hA x y f hf (ε / 2) (by positivity)
    have hev : ∀ᶠ M in atTop, ENNReal.ofReal (‖⟪deltaZ x, cfc f A (deltaZ y)⟫_ℂ‖ - ε) ≤
        ENNReal.ofReal (corrQR J h M x y (restrictZ M q)) :=
      eventually_atTop.mpr ⟨M0, fun M hM => ENNReal.ofReal_le_ofReal (by linarith [hM0 M hM])⟩
    have h1 := (ENNReal.ofReal_le_iff_le_toReal hLtop).mp (le_liminf_of_le (h := hev))
    linarith
  calc ENNReal.ofReal (corrQZ A x y) ≤ ENNReal.ofReal L.toReal := ENNReal.ofReal_le_ofReal hbound
    _ = L := ENNReal.ofReal_toReal hLtop

/-- Finite volume: `E Q^{Λ_M}(x, y) ≤ F(|x - y|)` under (EC). -/
lemma lintegral_corrQR_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (M : ℕ) (x y : ℤ)
    (hx : |x| ≤ M) (hy : |y| ≤ M) :
    ∫⁻ q, ENNReal.ofReal (corrQR J h M x y (restrictZ M q)) ∂muZ ≤
      ENNReal.ofReal (F (x - y).natAbs) := by
  have hxy : |x| ≤ M ∧ |y| ≤ M := ⟨hx, hy⟩
  set i := idxZ M x hx
  set j := idxZ M y hy
  have hcorr : ∀ q' : Cfg (2 * M + 1), 0 ≤ corrQ (hq J h q') (hq_isHermitian J h q') i j :=
    fun q' => corrQ_nonneg _ _ _ _
  have he : ∀ q' : Cfg (2 * M + 1),
      corrQR J h M x y q' = corrQ (hq J h q') (hq_isHermitian J h q') i j := fun q' => by
    unfold corrQR
    rw [dite_eq_left hxy]
  rw [InfVolAux.lintegral_comp_restrictZ M (fun q' => ENNReal.ofReal (corrQR J h M x y q'))]
  simp_rw [he]
  calc ∑ q' : Cfg (2 * M + 1), ENNReal.ofReal (corrQ (hq J h q') (hq_isHermitian J h q') i j) *
          ((2 : ENNReal) ^ (2 * M + 1))⁻¹
      = ENNReal.ofReal (avgQ fun q' => corrQ (hq J h q') (hq_isHermitian J h q') i j) := by
        unfold avgQ
        rw [ENNReal.ofReal_mul (by positivity), ENNReal.ofReal_sum_of_nonneg fun q' _ => hcorr q',
          Finset.mul_sum]
        refine Finset.sum_congr rfl fun q' _ => ?_
        rw [mul_comm]
        congr 1
        rw [one_div, ENNReal.ofReal_inv_of_pos (by positivity), ENNReal.ofReal_pow (by norm_num),
          ENNReal.ofReal_ofNat]
    _ ≤ ENNReal.ofReal (F (x - y).natAbs) := by
        refine ENNReal.ofReal_le_ofReal ?_
        have := hEC.bound (2 * M + 1) i j
        rwa [InfVolAux.dist_idxZ] at this

end ECInfVolQAux

/-- **Annealed localization on `ℤ`** (paper: Corollary `cor:ECinfvol`): `E Q^ℤ(x, y) ≤ C e^{-η|x-y|}`, for the operators
`A_q` (determined by their matrix elements; no measurability is assumed, the bound is for the lower integral). -/
theorem lintegral_corrQZ_le (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ ∀ (A : CfgZ → L2Z →L[ℂ] L2Z), (∀ q, IsAZ (-(2 * h / J)) q (A q)) →
      ∀ x y : ℤ, ∫⁻ q, ENNReal.ofReal (corrQZ (A q) x y) ∂muZ ≤
        ENNReal.ofReal (C * Real.exp (-η * (x - y).natAbs)) := by
  obtain ⟨C, η, hη, hEC⟩ := ec_of_ldt J h hJ hL
  refine ⟨C, η, hη, fun A hA x y => ?_⟩
  have hmeas : ∀ M : ℕ, Measurable fun q : CfgZ =>
      ENNReal.ofReal (ECInfVolQAux.corrQR J h M x y (restrictZ M q)) := fun M =>
    InfVolAux.measurable_comp_restrictZ M
      (fun q' => ENNReal.ofReal (ECInfVolQAux.corrQR J h M x y q'))
  calc ∫⁻ q, ENNReal.ofReal (corrQZ (A q) x y) ∂muZ
      ≤ ∫⁻ q, liminf (fun M : ℕ =>
          ENNReal.ofReal (ECInfVolQAux.corrQR J h M x y (restrictZ M q))) atTop ∂muZ :=
        lintegral_mono fun q => ECInfVolQAux.ofReal_corrQZ_le_liminf J h hJ q (A q) (hA q) x y
    _ ≤ liminf (fun M : ℕ =>
          ∫⁻ q, ENNReal.ofReal (ECInfVolQAux.corrQR J h M x y (restrictZ M q)) ∂muZ) atTop :=
        lintegral_liminf_le hmeas
    _ ≤ ENNReal.ofReal (min 1 (C * Real.exp (-η * (x - y).natAbs))) := by
        refine liminf_le_of_frequently_le' (Eventually.frequently ?_)
        filter_upwards [eventually_ge_atTop (x.natAbs + y.natAbs)] with M hM
        have hx : |x| ≤ M := by rw [Int.abs_eq_natAbs]; omega
        have hy : |y| ≤ M := by rw [Int.abs_eq_natAbs]; omega
        exact ECInfVolQAux.lintegral_corrQR_le J h hEC M x y hx hy
    _ ≤ ENNReal.ofReal (C * Real.exp (-η * (x - y).natAbs)) :=
        ENNReal.ofReal_le_ofReal (min_le_right _ _)

end DFL
