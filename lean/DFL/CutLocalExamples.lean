import DFL.Entropies

/-!
# Sharpness of the duality defect (paper: Lemma `lem:cutlocal`(d))

Two sites, cut after the first site (`N = 2`, `ℓ = 1`). The same function
`exProd (n, y) = 2^{-1/2} [n = 0] [y₁ = 0]`, read either in the dual variables `(n, q)` or in the
original variables `(n, x)`, is a product state across the cut. Its image under `U`
(`(UΦ)(n, x) = Φ(n, chargeOf n x)`) and under `U†` (`(U†Ψ)(n, q) = Ψ(n, linkOf n q)`) is a Bell state of
the two links, resp. of the two charges, with entanglement entropy `ln 2`. So the bound `ln 2` of
`DFL.cutlocal_defect` is attained in both directions.
-/

open Matrix Finset

noncomputable section

namespace DFL

/-- `2^{-1/2} |0, 0⟩ ⊗ (|y₀ = 0⟩ + |y₀ = 1⟩) ⊗ |y₁ = 0⟩`. -/
def exProd : Cfg 2 × Cfg 2 → ℂ :=
  fun p => if p.1 = 0 ∧ p.2 1 = 0 then ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ else 0

namespace CutLocalExamplesAux

/-- The reduced density matrix of a state supported on the graph `b = g a` of a map from the
configurations of `A` to those of the complement. -/
lemma rdmLGT_graph {N ℓ : ℕ} (hℓ : ℓ ≤ N) (Ψ : Cfg N × Cfg N → ℂ)
    (g : Cfg ℓ × Cfg ℓ → Cfg (N - ℓ) × Cfg (N - ℓ)) (f : Cfg ℓ × Cfg ℓ → ℂ)
    (hΨ : ∀ a b, Ψ ((splitLGT ℓ hℓ).symm (a, b)) = if b = g a then f a else 0)
    (a a' : Cfg ℓ × Cfg ℓ) :
    rdmLGT ℓ hℓ Ψ a a' = if g a = g a' then f a * star (f a') else 0 := by
  simp only [rdmLGT, hΨ, ite_mul, zero_mul]
  rw [Finset.sum_ite_eq' Finset.univ (g a)]
  simp only [Finset.mem_univ, ite_true]
  split_ifs <;> simp

/-- `|2^{-1/2}|² = 1/2`. -/
lemma c_mul_star :
    ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ * star (((Real.sqrt 2 : ℝ) : ℂ)⁻¹) = ((1 / 2 : ℝ) : ℂ) := by
  rw [← Complex.ofReal_inv, Complex.star_def, Complex.conj_ofReal, ← Complex.ofReal_mul]
  congr 1
  rw [← mul_inv, Real.mul_self_sqrt (by norm_num : (0 : ℝ) ≤ 2)]
  norm_num

/-- Rewriting the support condition of an indicator. -/
lemma ite_eq_of_iff {P Q R : Prop} [Decidable P] [Decidable Q] [Decidable R] (h : P ↔ Q ∧ R)
    (x : ℂ) : (if P then x else 0) = if Q then (if R then x else 0) else 0 := by
  by_cases hQ : Q <;> by_cases hR : R <;> simp [h, hQ, hR]

/-- The amplitude `2^{-1/2} [α = 0]` on the region `A`. -/
def amp (a : Cfg 1 × Cfg 1) : ℂ := if a.1 = 0 then ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ else 0

/-- Product of two amplitudes. -/
lemma amp_mul_star (a a' : Cfg 1 × Cfg 1) :
    amp a * star (amp a') = if a.1 = 0 ∧ a'.1 = 0 then ((1 / 2 : ℝ) : ℂ) else 0 := by
  unfold amp
  by_cases h : a.1 = 0
  · by_cases h' : a'.1 = 0
    · rw [ite_eq_left h, ite_eq_left h', ite_eq_left ⟨h, h'⟩, c_mul_star]
    · simp [h']
  · simp [h]

/-- The complement configuration `(0, β₀)` carried by the Bell states. -/
def gBell (a : Cfg 1 × Cfg 1) : Cfg (2 - 1) × Cfg (2 - 1) := (0, fun _ => a.2 0)

/-- Support of `exProd` in split coordinates. -/
lemma exProd_support : ∀ (α β : Cfg 1) (γ δ : Cfg (2 - 1)),
    ((splitCfg 1 (by norm_num)).symm (α, γ) = 0 ∧ (splitCfg 1 (by norm_num)).symm (β, δ) 1 = 0) ↔
      ((γ, δ) = ((0 : Cfg (2 - 1)), (0 : Cfg (2 - 1))) ∧ α = 0) := by
  decide

/-- Support of `U exProd` in split coordinates. -/
lemma exProdU_support : ∀ (α β : Cfg 1) (γ δ : Cfg (2 - 1)),
    ((splitCfg 1 (by norm_num)).symm (α, γ) = 0 ∧
      chargeOf ((splitCfg 1 (by norm_num)).symm (α, γ)) ((splitCfg 1 (by norm_num)).symm (β, δ)) 1
        = 0) ↔ ((γ, δ) = gBell (α, β) ∧ α = 0) := by
  decide

/-- Support of `U† exProd` in split coordinates. -/
lemma exProdUdag_support : ∀ (α β : Cfg 1) (γ δ : Cfg (2 - 1)),
    ((splitCfg 1 (by norm_num)).symm (α, γ) = 0 ∧
      linkOf ((splitCfg 1 (by norm_num)).symm (α, γ)) ((splitCfg 1 (by norm_num)).symm (β, δ)) 1
        = 0) ↔ ((γ, δ) = gBell (α, β) ∧ α = 0) := by
  decide

/-- On the support of the Bell states, `gBell` is injective. -/
lemma gBell_iff : ∀ a a' : Cfg 1 × Cfg 1,
    (gBell a = gBell a' ∧ (a.1 = 0 ∧ a'.1 = 0)) ↔ (a = a' ∧ a.1 = 0) := by
  decide

/-- The number of configurations of `A` with empty site. -/
lemma card_support : (univ.filter (fun a : Cfg 1 × Cfg 1 => a.1 = 0)).card = 2 := by
  decide

/-- The column vector of the amplitudes `amp`. -/
def ampCol : Matrix (Cfg 1 × Cfg 1) Unit ℂ := Matrix.of fun a _ => amp a

/-- `ampCol ampColᴴ` is the rank-one matrix of the amplitudes. -/
lemma ampCol_mul_conjTranspose_apply (a a' : Cfg 1 × Cfg 1) :
    (ampCol * ampColᴴ) a a' = amp a * star (amp a') := by
  simp [ampCol, Matrix.mul_apply]

/-- The amplitudes are normalized: `ampColᴴ ampCol = 1`. -/
lemma conjTranspose_ampCol_mul :
    ampColᴴ * ampCol = diagonal (fun _ : Unit => ((1 : ℝ) : ℂ)) := by
  have h1 : ∀ x : Cfg 1 × Cfg 1,
      star (amp x) * amp x = if x.1 = 0 then ((1 / 2 : ℝ) : ℂ) else 0 := by
    intro x
    rw [mul_comm, amp_mul_star]
    simp only [and_self]
  ext u u'
  rw [Subsingleton.elim u' u, diagonal_apply_eq, Matrix.mul_apply]
  simp only [conjTranspose_apply, ampCol, Matrix.of_apply, h1]
  rw [Finset.sum_ite, Finset.sum_const_zero, add_zero, Finset.sum_const, card_support]
  push_cast
  norm_num

/-- The entropy of the diagonal matrix `diag(½ [α = 0])` is `ln 2`. -/
lemma sum_negMulLog_half :
    ∑ a : Cfg 1 × Cfg 1, Real.negMulLog (if a.1 = 0 then 1 / 2 else 0) = Real.log 2 := by
  have h : ∀ a : Cfg 1 × Cfg 1, Real.negMulLog (if a.1 = 0 then 1 / 2 else 0) =
      if a.1 = 0 then Real.negMulLog (1 / 2) else 0 := by
    intro a
    split_ifs <;> simp
  simp only [h]
  rw [Finset.sum_ite, Finset.sum_const_zero, add_zero, Finset.sum_const, card_support,
    nsmul_eq_mul]
  simp only [Real.negMulLog, one_div, Real.log_inv]
  ring

/-- The reduced matrix of a Bell state supported on `b = gBell a` is `diag(½ [α = 0])`. -/
lemma rdmLGT_bell (Ψ : Cfg 2 × Cfg 2 → ℂ)
    (hΨ : ∀ a b, Ψ ((splitLGT 1 (by norm_num)).symm (a, b)) = if b = gBell a then amp a else 0) :
    rdmLGT 1 (by norm_num) Ψ =
      diagonal (fun a => (((if a.1 = 0 then 1 / 2 else 0 : ℝ)) : ℂ)) := by
  ext a a'
  rw [rdmLGT_graph _ Ψ gBell amp hΨ, diagonal_apply]
  have key := gBell_iff a a'
  by_cases ha : a = a'
  · subst ha
    rw [ite_eq_left rfl, ite_eq_left rfl, amp_mul_star]
    by_cases h0 : a.1 = 0 <;> simp [h0]
  · rw [ite_eq_right ha]
    split_ifs with hg
    · rw [amp_mul_star, ite_eq_right]
      intro h
      exact ha (key.1 ⟨hg, h⟩).1
    · rfl

/-- The entropy of the reduced matrix of a Bell state supported on `b = gBell a` is `ln 2`. -/
lemma vonNeumann_bell (Ψ : Cfg 2 × Cfg 2 → ℂ)
    (hΨ : ∀ a b, Ψ ((splitLGT 1 (by norm_num)).symm (a, b)) = if b = gBell a then amp a else 0) :
    vonNeumann (rdmLGT 1 (by norm_num) Ψ) (rdmLGT_isHermitian 1 (by norm_num) _) = Real.log 2 := by
  have hD : (diagonal (fun a : Cfg 1 × Cfg 1 =>
      (((if a.1 = 0 then 1 / 2 else 0 : ℝ)) : ℂ))).IsHermitian :=
    rdmLGT_bell Ψ hΨ ▸ rdmLGT_isHermitian 1 (by norm_num) Ψ
  rw [CutLocalAux.vonNeumann_congr _ hD (rdmLGT_bell Ψ hΨ), vonNeumann_diagonal]
  exact sum_negMulLog_half

end CutLocalExamplesAux

/-- `exProd` is a product state across the cut: zero entanglement entropy. -/
theorem exProd_entropy :
    vonNeumann (rdmLGT 1 (by norm_num) exProd) (rdmLGT_isHermitian 1 (by norm_num) _) = 0 := by
  have hΨ : ∀ a b, exProd ((splitLGT 1 (by norm_num)).symm (a, b)) =
      if b = ((0 : Cfg (2 - 1)), (0 : Cfg (2 - 1))) then CutLocalExamplesAux.amp a else 0 := by
    rintro ⟨α, β⟩ ⟨γ, δ⟩
    rw [CutLocalAux.splitLGT_symm_apply]
    exact CutLocalExamplesAux.ite_eq_of_iff (CutLocalExamplesAux.exProd_support α β γ δ) _
  have hEq : rdmLGT 1 (by norm_num) exProd =
      CutLocalExamplesAux.ampCol * CutLocalExamplesAux.ampColᴴ := by
    ext a a'
    rw [CutLocalExamplesAux.rdmLGT_graph _ exProd _ CutLocalExamplesAux.amp hΨ, ite_eq_left rfl,
      CutLocalExamplesAux.ampCol_mul_conjTranspose_apply]
  have hD : (diagonal (fun _ : Unit => ((1 : ℝ) : ℂ))).IsHermitian :=
    CutLocalExamplesAux.conjTranspose_ampCol_mul ▸
      isHermitian_conjTranspose_mul_self CutLocalExamplesAux.ampCol
  rw [CutLocalAux.vonNeumann_congr _
      (isHermitian_mul_conjTranspose_self CutLocalExamplesAux.ampCol) hEq,
    vonNeumann_mul_conjTranspose_comm,
    CutLocalAux.vonNeumann_congr _ hD CutLocalExamplesAux.conjTranspose_ampCol_mul,
    vonNeumann_diagonal]
  simp

/-- **`U` creates entanglement `ln 2`** (paper: Lemma `lem:cutlocal`(d), first example). -/
theorem exProd_U_entropy :
    vonNeumann (rdmLGT 1 (by norm_num) (fun p => exProd (dualEquiv.symm p)))
      (rdmLGT_isHermitian 1 (by norm_num) _) = Real.log 2 := by
  refine CutLocalExamplesAux.vonNeumann_bell _ ?_
  rintro ⟨α, β⟩ ⟨γ, δ⟩
  rw [CutLocalAux.splitLGT_symm_apply, dualEquiv_symm_apply]
  exact CutLocalExamplesAux.ite_eq_of_iff (CutLocalExamplesAux.exProdU_support α β γ δ) _

/-- **`U†` creates entanglement `ln 2`** (paper: Lemma `lem:cutlocal`(d), second example). -/
theorem exProd_Udag_entropy :
    vonNeumann (rdmLGT 1 (by norm_num) (fun p => exProd (dualEquiv p)))
      (rdmLGT_isHermitian 1 (by norm_num) _) = Real.log 2 := by
  refine CutLocalExamplesAux.vonNeumann_bell _ ?_
  rintro ⟨α, β⟩ ⟨γ, δ⟩
  rw [CutLocalAux.splitLGT_symm_apply, dualEquiv_apply]
  exact CutLocalExamplesAux.ite_eq_of_iff (CutLocalExamplesAux.exProdUdag_support α β γ δ) _

end DFL
