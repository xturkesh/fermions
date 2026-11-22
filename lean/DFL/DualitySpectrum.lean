import DFL.SlaterExact
import DFL.Imbalance

/-!
# The spectrum of the gauge theory (paper: Theorem `thm:duality`(b))

The second quantization `dΓ(A)` of a real symmetric `A` with eigenvalues `λ_1, …, λ_N` has the
eigenvalues `∑_{α ∈ T} λ_α`, `T ⊆ {1, …, N}`, with multiplicities (`charpoly_dGamma`): the Slater
determinants of the eigenvectors of `A` form an orthonormal eigenbasis. Hence, by the duality, the spectrum
of `H` (with multiplicities) is the union over `q` of the multisets `{∑_α m_α E_α(q) - h ∑_j q_j}`,
`m ∈ {0,1}^N`, where `E_α(q)` are the eigenvalues of `h_q` (`charpoly_H`); for `J ≠ 0` they are simple
(`hq_eigenvalues_injective`).
-/

open Matrix Finset Polynomial

noncomputable section

namespace DFL

variable {N : ℕ}

namespace DualitySpectrumAux

open AreaLawAux SlaterExactAux

/-! ### Linear algebra -/

/-- A matrix with an orthonormal eigenbasis `v_k` (eigenvalues `μ_k`) has characteristic
polynomial `∏_k (X - μ_k)`. -/
lemma charpoly_eq_prod_of_orthonormal {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ]
    [DecidableEq κ] (hc : Fintype.card κ = Fintype.card ι) (M : Matrix ι ι ℂ) (v : κ → ι → ℂ)
    (μ : κ → ℂ) (horth : ∀ k l, star (v k) ⬝ᵥ v l = if k = l then 1 else 0)
    (heig : ∀ k, M *ᵥ v k = μ k • v k) :
    M.charpoly = ∏ k, (X - C (μ k)) := by
  let P : Matrix ι κ ℂ := Matrix.of fun i k => v k i
  have hPP : Pᴴ * P = 1 := by
    ext k l
    rw [Matrix.mul_apply, Matrix.one_apply, ← horth k l]
    simp [P, dotProduct, Matrix.conjTranspose_apply]
  have hPP' : P * Pᴴ = 1 := (Matrix.mul_eq_one_comm_of_equiv (Fintype.equivOfCardEq hc)).mp hPP
  have hMP : M * P = P * diagonal μ := by
    ext i k
    rw [Matrix.mul_diagonal]
    have h := congrFun (heig k) i
    simp only [Matrix.mulVec, dotProduct, Pi.smul_apply, smul_eq_mul] at h
    simp only [Matrix.mul_apply, P, Matrix.of_apply]
    rw [h, mul_comm]
  have key := Matrix.charpoly_mul_comm' (M * P) Pᴴ
  rw [Matrix.mul_assoc, hPP', Matrix.mul_one, hMP, ← Matrix.mul_assoc, hPP, Matrix.one_mul,
    Matrix.charpoly_diagonal, hc] at key
  exact (isRegular_X_pow _).left key

/-- The characteristic polynomial of a block-diagonal matrix is the product of those of the
blocks. -/
lemma charpoly_blockDiagonal {m o R : Type*} [Fintype m] [DecidableEq m] [Fintype o]
    [DecidableEq o] [CommRing R] (M : o → Matrix m m R) :
    (blockDiagonal M).charpoly = ∏ k, (M k).charpoly := by
  have h : charmatrix (blockDiagonal M) = blockDiagonal (fun k => charmatrix (M k)) := by
    ext ⟨i, k⟩ ⟨j, k'⟩
    simp only [charmatrix_apply, blockDiagonal_apply, diagonal_apply, Prod.mk.injEq]
    by_cases hk : k = k' <;> by_cases hij : i = j <;> simp [hk, hij]
  rw [charpoly, h, det_blockDiagonal]
  rfl

/-! ### Complexification -/

/-- Complexification commutes with finite sums. -/
lemma cplx_sum' {ι κ : Type*} [Fintype ι] [DecidableEq ι] (s : Finset κ)
    (f : κ → Matrix ι ι ℝ) : cplx (∑ k ∈ s, f k) = ∑ k ∈ s, cplx (f k) := by
  ext i j; simp [cplx, Matrix.sum_apply]

/-- Complexification of a real multiple. -/
lemma cplx_smul' {ι : Type*} [Fintype ι] [DecidableEq ι] (r : ℝ) (M : Matrix ι ι ℝ) :
    cplx (r • M) = (r : ℂ) • cplx M := by
  ext i j; simp [cplx]

/-- Complexification of a difference. -/
lemma cplx_sub' {ι : Type*} [Fintype ι] [DecidableEq ι] (M M' : Matrix ι ι ℝ) :
    cplx (M - M') = cplx M - cplx M' := by
  ext i j; simp [cplx]

/-! ### Commutators with `dΓ(A)` -/

/-- `[c_z† c_w, c_x†] = δ_{wx} c_z†`. -/
lemma comm_cre_mul_ann_cre (z w x : Fin N) :
    cre z * ann w * cre x - cre x * (cre z * ann w) = if w = x then cre z else 0 := by
  have h1 : ann w * cre x = (if w = x then 1 else 0) - cre x * ann w :=
    eq_sub_of_add_eq (ann_mul_cre_add w x)
  have h2 : cre x * cre z = -(cre z * cre x) := eq_neg_of_add_eq_zero_left (cre_mul_cre_add x z)
  calc cre z * ann w * cre x - cre x * (cre z * ann w)
      = cre z * (ann w * cre x) - (cre x * cre z) * ann w := by simp only [Matrix.mul_assoc]
    _ = if w = x then cre z else 0 := by
        rw [h1, h2]
        split_ifs <;> noncomm_ring

/-- `[dΓ(A), c_x†] = ∑_z A_{zx} c_z†`. -/
lemma dGamma_comm_cre (A : Matrix (Fin N) (Fin N) ℝ) (x : Fin N) :
    dGamma A * cre x - cre x * dGamma A = ∑ z, A z x • cre z := by
  have key : ∀ z w, (A z w • (cre z * ann w)) * cre x - cre x * (A z w • (cre z * ann w)) =
      if w = x then A z w • cre z else 0 := by
    intro z w
    rw [Matrix.smul_mul, Matrix.mul_smul, ← smul_sub, comm_cre_mul_ann_cre]
    split_ifs <;> simp
  simp only [dGamma, Finset.sum_mul, Finset.mul_sum, ← Finset.sum_sub_distrib, key]
  simp [Finset.sum_ite_eq']

/-- `c†(c v) = c c†(v)`. -/
lemma crOp_smul (c : ℂ) (v : Fin N → ℂ) : crOp (c • v) = c • crOp v := by
  simp [crOp, Finset.smul_sum, smul_smul]

/-- `[dΓ(A), c†(v)] = c†(A v)` (complexified). -/
lemma dGamma_comm_crOp (A : Matrix (Fin N) (Fin N) ℝ) (v : Fin N → ℂ) :
    cplx (dGamma A) * crOp v - crOp v * cplx (dGamma A) = crOp (cplx A *ᵥ v) := by
  have h1 : ∀ x, cplx (dGamma A) * cplx (cre x) - cplx (cre x) * cplx (dGamma A) =
      ∑ z, (A z x : ℂ) • cplx (cre z) := by
    intro x
    rw [← cplx_mul', ← cplx_mul', ← cplx_sub', dGamma_comm_cre, cplx_sum']
    simp only [cplx_smul']
  simp only [crOp, Finset.mul_sum, Finset.sum_mul, Matrix.mul_smul, Matrix.smul_mul,
    ← Finset.sum_sub_distrib, ← smul_sub, h1, Finset.smul_sum, smul_smul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun z _ => ?_
  rw [← Finset.sum_smul]
  congr 1
  simp [Matrix.mulVec, dotProduct, cplx, mul_comm]

/-- `dΓ(A) |0⟩ = 0`. -/
lemma dGamma_vac (A : Matrix (Fin N) (Fin N) ℝ) :
    cplx (dGamma A) *ᵥ fockVec (∅ : Finset (Fin N)) = 0 := by
  have h : cplx (dGamma A) = ∑ x, ∑ y, (A x y : ℂ) • (cplx (cre x) * cplx (ann y)) := by
    simp only [dGamma, cplx_sum', cplx_smul', cplx_mul']
  rw [h]
  simp [Matrix.sum_mulVec, Matrix.smul_mulVec, ← Matrix.mulVec_mulVec, ann_mulVec_vac]

/-- Slater determinants of eigenvectors are eigenvectors of `dΓ(A)`:
`dΓ(A) ∏_{σ k} c†(u_k) |0⟩ = (∑_{σ k} λ_k) ∏_{σ k} c†(u_k) |0⟩`. -/
lemma dGamma_prodSel (A : Matrix (Fin N) (Fin N) ℝ) {M : ℕ} (u : Fin M → Fin N → ℂ)
    (lam : Fin M → ℂ) (hu : ∀ k, cplx A *ᵥ u k = lam k • u k) (σ : Fin M → Bool) :
    cplx (dGamma A) *ᵥ prodSel u σ = (∑ k, if σ k then lam k else 0) • prodSel u σ := by
  induction M with
  | zero => simp [prodSel_zero, dGamma_vac]
  | succ M ih =>
    have ih' := ih (fun k => u k.succ) (fun k => lam k.succ) (fun k => hu k.succ)
      (fun k => σ k.succ)
    rw [prodSel_succ, Fin.sum_univ_succ]
    cases h0 : σ 0
    · simp only [Bool.false_eq_true, ite_false, Matrix.one_mulVec, zero_add]
      exact ih'
    · simp only [ite_true]
      have hc := dGamma_comm_crOp A (u 0)
      rw [hu 0, sub_eq_iff_eq_add] at hc
      rw [Matrix.mulVec_mulVec, hc, Matrix.add_mulVec, ← Matrix.mulVec_mulVec, ih', crOp_smul,
        Matrix.smul_mulVec, Matrix.mulVec_smul, add_smul, add_comm]

/-! ### The Slater eigenbasis -/

/-- The complexified eigenvectors of a real symmetric matrix. -/
def eigC (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (α : Fin N) : Fin N → ℂ :=
  fun x => (hA.eigenvectorBasis α x : ℂ)

/-- Orthonormality of the complexified eigenvectors. -/
lemma eigC_orth (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (α β : Fin N) :
    star (eigC A hA α) ⬝ᵥ eigC A hA β = if α = β then 1 else 0 := by
  have h := ImbalanceAux.eig_orth A hA α β
  have e : star (eigC A hA α) ⬝ᵥ eigC A hA β =
      ((∑ x, hA.eigenvectorBasis α x * hA.eigenvectorBasis β x : ℝ) : ℂ) := by
    simp [eigC, dotProduct]
  rw [e, h]
  split_ifs <;> simp

/-- Eigenvalue equation of the complexified eigenvectors. -/
lemma eigC_eig (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (α : Fin N) :
    cplx A *ᵥ eigC A hA α = (hA.eigenvalues α : ℂ) • eigC A hA α := by
  funext x
  have h := ImbalanceAux.eig_eq A hA α x
  simp only [Matrix.mulVec, dotProduct, cplx, Matrix.map_apply, eigC, Pi.smul_apply, smul_eq_mul]
  exact_mod_cast h

/-- The Slater determinant `∏_{α ∈ T} c†(φ_α) |0⟩` of the eigenvectors. -/
def slaterT (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (T : Finset (Fin N)) :
    Cfg N → ℂ :=
  prodSel (eigC A hA) (fun α => decide (α ∈ T))

/-- The Slater determinants of the eigenvectors are orthonormal. -/
lemma slaterT_orth (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (T T' : Finset (Fin N)) :
    star (slaterT A hA T) ⬝ᵥ slaterT A hA T' = if T = T' then 1 else 0 := by
  have h := gram_prodSel (eigC A hA) (fun _ => 1) (fun k l => by rw [eigC_orth])
    (fun α => decide (α ∈ T)) (fun α => decide (α ∈ T'))
  rw [slaterT, slaterT, h]
  have hiff : ((fun α => decide (α ∈ T)) = fun α => decide (α ∈ T')) ↔ T = T' := by
    constructor
    · intro hT
      ext α
      have := congrFun hT α
      simpa using this
    · rintro rfl; rfl
  simp only [hiff, ite_self, Finset.prod_const_one]

/-- `dΓ(A) - c` acts diagonally on the Slater determinants of the eigenvectors. -/
lemma slaterT_eig (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (c : ℝ)
    (T : Finset (Fin N)) :
    cplx (dGamma A - c • 1) *ᵥ slaterT A hA T =
      ((∑ α ∈ T, hA.eigenvalues α - c : ℝ) : ℂ) • slaterT A hA T := by
  rw [cplx_sub', cplx_smul', cplx_one', Matrix.sub_mulVec, Matrix.smul_mulVec, Matrix.one_mulVec,
    slaterT, dGamma_prodSel A (eigC A hA) (fun α => (hA.eigenvalues α : ℂ)) (eigC_eig A hA),
    ← sub_smul]
  congr 1
  simp [Finset.sum_ite_mem]

/-- **Spectrum of `dΓ(A) - c`**: `char (dΓ(A) - c) = ∏_{T} (X - (∑_{α ∈ T} λ_α - c))`. -/
lemma charpoly_dGamma_sub (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (c : ℝ) :
    (dGamma A - c • 1).charpoly =
      ∏ T : Finset (Fin N), (X - C (∑ α ∈ T, hA.eigenvalues α - c)) := by
  apply Polynomial.map_injective Complex.ofRealHom Complex.ofReal_injective
  rw [← Matrix.charpoly_map, Polynomial.map_prod]
  simp only [Polynomial.map_sub, Polynomial.map_X, Polynomial.map_C, Complex.ofRealHom_eq_coe]
  have hc : Fintype.card (Finset (Fin N)) = Fintype.card (Cfg N) := by
    simp [Fintype.card_finset, ZMod.card]
  exact charpoly_eq_prod_of_orthonormal hc _ (slaterT A hA) _ (slaterT_orth A hA)
    (slaterT_eig A hA c)

end DualitySpectrumAux

open DualitySpectrumAux

/-- **Spectrum of a quadratic Hamiltonian**: `char dΓ(A) = ∏_{T} (X - ∑_{α ∈ T} λ_α)`. -/
theorem charpoly_dGamma (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) :
    (dGamma A).charpoly = ∏ T : Finset (Fin N), (X - C (∑ α ∈ T, hA.eigenvalues α)) := by
  have h := charpoly_dGamma_sub A hA 0
  simpa using h

/-- **The spectrum of `H`** (paper: Theorem `thm:duality`(b)): with multiplicities, the union over the
charge sectors `q` of `{∑_{α ∈ T} E_α(q) - h ∑_j q_j : T ⊆ {1, …, N}}`. -/
theorem charpoly_H (J h : ℝ) :
    (H (N := N) J h).charpoly =
      ∏ q : Cfg N, ∏ T : Finset (Fin N),
        (X - C (∑ α ∈ T, (hq_isHermitian J h q).eigenvalues α - h * ∑ j, sgn (q j))) := by
  rw [← Matrix.charpoly_reindex (dualEquiv (N := N)).symm (H J h)]
  have e : Matrix.reindex dualEquiv.symm dualEquiv.symm (H J h) = Hdual (N := N) J h := by
    rw [Matrix.reindex_apply, Equiv.symm_symm]
    exact duality J h
  rw [e, Hdual, charpoly_blockDiagonal]
  refine Finset.prod_congr rfl fun q _ => ?_
  rw [Hsec]
  exact charpoly_dGamma_sub (hq J h q) (hq_isHermitian J h q) _

end DFL
