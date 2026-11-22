import DFL.Imbalance

/-!
# Normalization of the one-body operators (paper: Lemma `lem:normalization`)

With `Γ = diag((-1)^x)` (`stag`):

* `Γ h_q Γ⁻¹ = J (Δ + λ diag(q))`, `λ = 2h/J` (`stag_conj_hq`, `stag_conj_hq_eq_smul`): the sign of
  the hopping is irrelevant;
* consequently the eigenfunction correlator is the same for `J` and `-J` (`corrQ_hq_neg`);
* the joint law of `(Q(x, y))_{x,y}` is the same for `λ` and `-λ` (`avgQ_law_neg_h`), because
  `h_q` at `-h` is `h_{-q}` at `h` (`hq_neg_h`) and `-q` has the same law as `q`.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

namespace NormalizationAux

/-- `(-1)^n (-1)^n = 1`. -/
lemma neg_one_pow_mul_self (n : ℕ) : (-1 : ℝ) ^ n * (-1) ^ n = 1 := by
  rw [← mul_pow, neg_one_mul, neg_neg, one_pow]

/-- `sgn (a + 1) = - sgn a` in `ZMod 2`. -/
lemma sgn_add_one (a : ZMod 2) : sgn (a + 1) = - sgn a := by
  have h01 : (1 : ZMod 2) ≠ 0 := by decide
  have h11 : (1 + 1 : ZMod 2) = 0 := by decide
  rcases (by decide : ∀ b : ZMod 2, b = 0 ∨ b = 1) a with rfl | rfl
  · simp [sgn, h01]
  · simp [sgn, h01, h11]

/-- Orthonormality of the eigenvectors: `∑_x φ_α(x) φ_β(x) = δ_{αβ}`. -/
lemma eig_orth (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (α β : Fin N) :
    ∑ x, hA.eigenvectorBasis α x * hA.eigenvectorBasis β x = if α = β then 1 else 0 := by
  have h := congrFun (congrFun (Unitary.coe_star_mul_self hA.eigenvectorUnitary) α) β
  simp only [Matrix.mul_apply, Matrix.one_apply, Matrix.star_apply,
    Matrix.IsHermitian.eigenvectorUnitary_apply, star_trivial] at h
  exact h

/-- Completeness of the eigenvectors: `∑_α φ_α(x) φ_α(y) = δ_{xy}`. -/
lemma eig_compl (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x y : Fin N) :
    ∑ α, hA.eigenvectorBasis α x * hA.eigenvectorBasis α y = if x = y then 1 else 0 := by
  have h := congrFun (congrFun (Unitary.coe_mul_star_self hA.eigenvectorUnitary) x) y
  simp only [Matrix.mul_apply, Matrix.one_apply, Unitary.coe_star, Matrix.star_apply,
    Matrix.IsHermitian.eigenvectorUnitary_apply, star_trivial] at h
  exact h

/-- The eigenvalue equation, entrywise: `∑_y A_{xy} φ_α(y) = E_α φ_α(x)`. -/
lemma eig_eq (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (α x : Fin N) :
    ∑ y, A x y * hA.eigenvectorBasis α y = hA.eigenvalues α * hA.eigenvectorBasis α x := by
  have h := congrFun (hA.mulVec_eigenvectorBasis α) x
  simpa [Matrix.mulVec, dotProduct] using h

/-- For a real symmetric matrix with simple spectrum, a unit eigenvector agrees in absolute value,
entry by entry, with one of the vectors of Mathlib's eigenbasis, whose eigenvalue is the same. -/
lemma exists_abs_eq_eigvec (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian)
    (hinj : Function.Injective hA.eigenvalues) (μ : ℝ) (v : Fin N → ℝ)
    (hv : ∀ x, ∑ y, A x y * v y = μ * v x) (hnorm : ∑ x, v x ^ 2 = 1) :
    ∃ β, hA.eigenvalues β = μ ∧ ∀ x, |v x| = |hA.eigenvectorBasis β x| := by
  have hsymm : ∀ x y, A x y = A y x := fun x y => by
    have := hA.apply y x
    rwa [star_trivial] at this
  obtain ⟨c, hc⟩ : ∃ c : Fin N → ℝ, ∀ β, c β = ∑ x, hA.eigenvectorBasis β x * v x :=
    ⟨_, fun _ => rfl⟩
  -- the coefficient of `v` along `φ_β` vanishes unless `E_β = μ`
  have h1 : ∀ β, c β * (hA.eigenvalues β - μ) = 0 := by
    intro β
    have e1 : ∑ x, hA.eigenvectorBasis β x * ∑ y, A x y * v y = μ * c β := by
      rw [hc β, Finset.mul_sum]
      refine Finset.sum_congr rfl fun x _ => ?_
      rw [hv x]; ring
    have e2 : ∑ x, hA.eigenvectorBasis β x * ∑ y, A x y * v y = hA.eigenvalues β * c β := by
      have : ∀ y, ∑ x, hA.eigenvectorBasis β x * (A x y * v y) =
          hA.eigenvalues β * (hA.eigenvectorBasis β y * v y) := by
        intro y
        calc ∑ x, hA.eigenvectorBasis β x * (A x y * v y)
            = (∑ x, A y x * hA.eigenvectorBasis β x) * v y := by
              rw [Finset.sum_mul]
              refine Finset.sum_congr rfl fun x _ => ?_
              rw [hsymm x y]; ring
          _ = hA.eigenvalues β * (hA.eigenvectorBasis β y * v y) := by
              rw [eig_eq A hA β y]; ring
      simp_rw [Finset.mul_sum]
      rw [Finset.sum_comm, Finset.sum_congr rfl fun y _ => this y, ← Finset.mul_sum, hc β]
    rw [mul_sub, mul_comm, ← e2, mul_comm (c β) μ, ← e1, sub_self]
  -- expansion of `v` in the eigenbasis
  have h2 : ∀ x, v x = ∑ β, c β * hA.eigenvectorBasis β x := by
    intro x
    simp only [hc, Finset.sum_mul]
    rw [Finset.sum_comm]
    have : ∀ y, ∑ β, hA.eigenvectorBasis β y * v y * hA.eigenvectorBasis β x =
        v y * (if y = x then 1 else 0) := by
      intro y
      rw [← eig_compl A hA y x, Finset.mul_sum]
      exact Finset.sum_congr rfl fun β _ => by ring
    rw [Finset.sum_congr rfl fun y _ => this y]
    simp
  -- Parseval
  have h3 : ∑ β, c β ^ 2 = 1 := by
    rw [← hnorm]
    have : ∀ x, v x ^ 2 = ∑ β, c β * (hA.eigenvectorBasis β x * v x) := by
      intro x
      have e : ∑ β, c β * (hA.eigenvectorBasis β x * v x) =
          (∑ β, c β * hA.eigenvectorBasis β x) * v x := by
        rw [Finset.sum_mul]; exact Finset.sum_congr rfl fun β _ => by ring
      rw [e, ← h2 x, sq]
    rw [Finset.sum_congr rfl fun x _ => this x, Finset.sum_comm]
    refine Finset.sum_congr rfl fun β _ => ?_
    rw [← Finset.mul_sum, ← hc β, sq]
  obtain ⟨β, hβ⟩ : ∃ β, c β ≠ 0 := by
    by_contra hcon
    have h0 : ∀ β, c β = 0 := fun β => by
      by_contra hβ; exact hcon ⟨β, hβ⟩
    simp [h0] at h3
  have hE : hA.eigenvalues β = μ := by
    rcases mul_eq_zero.mp (h1 β) with h | h
    · exact absurd h hβ
    · linarith
  have hzero : ∀ γ, γ ≠ β → c γ = 0 := by
    intro γ hγ
    rcases mul_eq_zero.mp (h1 γ) with h | h
    · exact h
    · exact absurd (hinj (by linarith : hA.eigenvalues γ = hA.eigenvalues β)) hγ
  have hv' : ∀ x, v x = c β * hA.eigenvectorBasis β x := by
    intro x
    rw [h2 x, Finset.sum_eq_single β]
    · intro γ _ hγ; rw [hzero γ hγ, zero_mul]
    · intro h; exact absurd (Finset.mem_univ β) h
  have hcβ : c β ^ 2 = 1 := by
    rw [← h3, Finset.sum_eq_single β]
    · intro γ _ hγ; rw [hzero γ hγ]; ring
    · intro h; exact absurd (Finset.mem_univ β) h
  have habs : |c β| = 1 := by
    have h' : |c β| ^ 2 = 1 := by rw [sq_abs]; exact hcβ
    nlinarith [abs_nonneg (c β)]
  refine ⟨β, hE, fun x => ?_⟩
  rw [hv' x, abs_mul, habs, one_mul]

/-- `corrQ` depends only on the matrix, not on the proof that it is Hermitian. -/
lemma corrQ_congr {A B : Matrix (Fin N) (Fin N) ℝ} (hAB : A = B) (hA : A.IsHermitian)
    (hB : B.IsHermitian) : corrQ A hA = corrQ B hB := by
  subst hAB; rfl

end NormalizationAux

/-- `Γ² = 1`. -/
theorem stag_mul_stag : stag N * stag N = 1 := by
  rw [stag, diagonal_mul_diagonal, ← diagonal_one]
  congr 1
  funext x
  exact NormalizationAux.neg_one_pow_mul_self x.val

/-- `Γ h_q Γ = h_q` with the opposite hopping sign. -/
theorem stag_conj_hq (J h : ℝ) (q : Cfg N) : stag N * hq J h q * stag N = hq (-J) h q := by
  ext x y
  simp only [stag, hq, Matrix.diagonal_mul, Matrix.mul_diagonal, Matrix.add_apply,
    Matrix.smul_apply, smul_eq_mul, Matrix.diagonal_apply, pathAdj]
  by_cases hxy : x = y
  · subst hxy
    have : ¬ (x.val + 1 = x.val ∨ x.val + 1 = x.val) := by omega
    simp only [this, ite_false, ite_true]
    linear_combination (2 * h * sgn (q x)) * NormalizationAux.neg_one_pow_mul_self x.val
  · simp only [hxy, ite_false, mul_zero, add_zero]
    split_ifs with hadj
    · rcases hadj with e | e
      · rw [← e, pow_succ]
        linear_combination J * NormalizationAux.neg_one_pow_mul_self x.val
      · rw [← e, pow_succ]
        linear_combination J * NormalizationAux.neg_one_pow_mul_self y.val
    · ring

/-- **Normalization** (paper: Lemma `lem:normalization`): `Γ h_q Γ⁻¹ = J (Δ + λ diag(q))`,
`λ = 2h/J`. -/
theorem stag_conj_hq_eq_smul (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) :
    stag N * hq J h q * stag N =
      J • (pathAdj N + (2 * h / J) • diagonal (fun x => sgn (q x))) := by
  rw [stag_conj_hq, hq, neg_neg, smul_add, smul_smul, mul_div_cancel₀ _ hJ]

/-- The eigenfunction correlator does not depend on the sign of the hopping. -/
theorem corrQ_hq_neg (J h : ℝ) (q : Cfg N) (x y : Fin N) :
    corrQ (hq (-J) h q) (hq_isHermitian (-J) h q) x y = corrQ (hq J h q) (hq_isHermitian J h q) x y := by
  by_cases hJ : J = 0
  · subst hJ
    rw [neg_zero]
  · set hA := hq_isHermitian J h q with hAdef
    set hB := hq_isHermitian (-J) h q with hBdef
    -- entrywise form of `h_q(J) = Γ h_q(-J) Γ`
    have hconj : ∀ x y, hq J h q x y = (-1 : ℝ) ^ x.val * hq (-J) h q x y * (-1) ^ y.val := by
      intro x y
      have := congrFun (congrFun (stag_conj_hq (-J) h q) x) y
      rw [neg_neg] at this
      rw [← this]
      simp only [stag, Matrix.diagonal_mul, Matrix.mul_diagonal]
    -- `Γ ψ_α` is, up to sign, an eigenvector `φ_{σ α}` of `h_q(J)`
    have key : ∀ α, ∃ β, hA.eigenvalues β = hB.eigenvalues α ∧
        ∀ x, |hB.eigenvectorBasis α x| = |hA.eigenvectorBasis β x| := by
      intro α
      have hv : ∀ x, ∑ y, hq J h q x y * ((-1 : ℝ) ^ y.val * hB.eigenvectorBasis α y) =
          hB.eigenvalues α * ((-1 : ℝ) ^ x.val * hB.eigenvectorBasis α x) := by
        intro x
        simp only [hconj]
        rw [show hB.eigenvalues α * ((-1 : ℝ) ^ x.val * hB.eigenvectorBasis α x) =
            (-1 : ℝ) ^ x.val * ∑ y, hq (-J) h q x y * hB.eigenvectorBasis α y by
          rw [NormalizationAux.eig_eq _ hB α x]; ring]
        rw [Finset.mul_sum]
        refine Finset.sum_congr rfl fun y _ => ?_
        linear_combination ((-1 : ℝ) ^ x.val * hq (-J) h q x y * hB.eigenvectorBasis α y) *
          NormalizationAux.neg_one_pow_mul_self y.val
      have hn : ∑ x, ((-1 : ℝ) ^ x.val * hB.eigenvectorBasis α x) ^ 2 = 1 := by
        have e := NormalizationAux.eig_orth _ hB α α
        rw [ite_eq_left rfl] at e
        refine (Finset.sum_congr rfl fun x _ => ?_).trans e
        linear_combination (hB.eigenvectorBasis α x) ^ 2 *
          NormalizationAux.neg_one_pow_mul_self x.val
      obtain ⟨β, hβ1, hβ2⟩ := NormalizationAux.exists_abs_eq_eigvec (hq J h q) hA
        (hq_eigenvalues_injective J h hJ q) (hB.eigenvalues α)
        (fun x => (-1 : ℝ) ^ x.val * hB.eigenvectorBasis α x) hv hn
      refine ⟨β, hβ1, fun x => ?_⟩
      have h' := hβ2 x
      simp only [abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul] at h'
      exact h'
    choose σ hσ1 hσ2 using key
    have hσinj : Function.Injective σ := by
      intro α α' e
      apply hq_eigenvalues_injective (-J) h (neg_ne_zero.mpr hJ) q
      rw [← hσ1 α, ← hσ1 α', e]
    have hσbij := hσinj.bijective_of_finite
    unfold corrQ
    calc ∑ α, |hB.eigenvectorBasis α x| * |hB.eigenvectorBasis α y|
        = ∑ α, |hA.eigenvectorBasis (σ α) x| * |hA.eigenvectorBasis (σ α) y| :=
          Finset.sum_congr rfl fun α _ => by rw [hσ2 α x, hσ2 α y]
      _ = ∑ β, |hA.eigenvectorBasis β x| * |hA.eigenvectorBasis β y| :=
          hσbij.sum_comp (fun β => |hA.eigenvectorBasis β x| * |hA.eigenvectorBasis β y|)

/-- Flipping the sign of `h` is the same as flipping all charges. -/
theorem hq_neg_h (J h : ℝ) (q : Cfg N) : hq J (-h) q = hq J h (fun x => q x + 1) := by
  ext x y
  simp only [hq, Matrix.add_apply, Matrix.smul_apply, smul_eq_mul, Matrix.diagonal_apply,
    NormalizationAux.sgn_add_one]
  split_ifs <;> ring

/-- **Symmetry `λ ↦ -λ`** (paper: Lemma `lem:normalization`): the joint law of the correlators
`(Q(x, y))_{x,y}` under i.i.d. uniform charges is the same for `h` and `-h`. -/
theorem avgQ_law_neg_h (J h : ℝ) (f : (Fin N → Fin N → ℝ) → ℝ) :
    avgQ (fun q : Cfg N => f (fun x y => corrQ (hq J (-h) q) (hq_isHermitian J (-h) q) x y)) =
      avgQ (fun q : Cfg N => f (fun x y => corrQ (hq J h q) (hq_isHermitian J h q) x y)) := by
  unfold avgQ
  congr 1
  have e : ∀ q : Cfg N, f (fun x y => corrQ (hq J (-h) q) (hq_isHermitian J (-h) q) x y) =
      f (fun x y => corrQ (hq J h (q + 1)) (hq_isHermitian J h (q + 1)) x y) := by
    intro q
    have e2 : hq J (-h) q = hq J h (q + 1) := hq_neg_h J h q
    rw [NormalizationAux.corrQ_congr e2 (hq_isHermitian J (-h) q) (hq_isHermitian J h (q + 1))]
  rw [Finset.sum_congr rfl fun q _ => e q]
  exact Fintype.sum_equiv (Equiv.addRight (1 : Cfg N)) _ _ (fun q => rfl)

end DFL
