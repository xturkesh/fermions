import DFL.FreeDyn
import DFL.Entropy
import DFL.OneBody

/-!
# Entanglement area law (paper: Theorem `thm:arealaw`, part on `S_G`)

For a Fock state evolved by a quadratic Hamiltonian (a Slater determinant), the von Neumann
entropy of a left segment `A = {0, …, ℓ-1}` satisfies
`S(ρ_A) = ∑_k h(ν_k) ≤ 2 ∑_k √(ν_k (1 - ν_k)) ≤ 2 ∑_{x ∈ A, y ∉ A} |C_{xy}|`,
where `ν_k` are the eigenvalues of the restricted correlation matrix `C_A`. Combined with the
duality, the sector averaging and Hypothesis (EC), this gives a uniform-in-time area law for the
sector-averaged (charge-conditioned) entanglement `S_G` of the lattice gauge theory.

For a left segment the Jordan–Wigner occupation basis factorizes canonically,
`|n⟩ = |n_A⟩ ⊗ |n_B⟩`, so the fermionic reduced density matrix is the ordinary partial trace in this
basis.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The left segment `{0, …, ℓ-1}`. -/
def segment (ℓ : ℕ) : Finset (Fin N) := univ.filter (fun x => x.val < ℓ)

/-- Splitting a configuration into its left segment `{0, …, ℓ-1}` and the rest. -/
def splitCfg (ℓ : ℕ) (hℓ : ℓ ≤ N) : Cfg N ≃ Cfg ℓ × Cfg (N - ℓ) :=
  (Equiv.arrowCongr (finCongr (by omega : N = ℓ + (N - ℓ))) (Equiv.refl _)).trans
    ((Equiv.arrowCongr finSumFinEquiv.symm (Equiv.refl _)).trans
      (Equiv.sumArrowEquivProdArrow _ _ _))

/-- Reduced density matrix of a pure state `v` on the left segment `{0, …, ℓ-1}`. -/
def rdm (ℓ : ℕ) (hℓ : ℓ ≤ N) (v : Cfg N → ℂ) : Matrix (Cfg ℓ) (Cfg ℓ) ℂ :=
  fun a a' => ∑ b : Cfg (N - ℓ),
    v ((splitCfg ℓ hℓ).symm (a, b)) * star (v ((splitCfg ℓ hℓ).symm (a', b)))

theorem rdm_isHermitian (ℓ : ℕ) (hℓ : ℓ ≤ N) (v : Cfg N → ℂ) : (rdm ℓ hℓ v).IsHermitian := by
  ext a a'
  simp only [Matrix.conjTranspose_apply, rdm, star_sum, star_mul', star_star]
  exact Finset.sum_congr rfl fun b _ => mul_comm _ _

/-- von Neumann entropy (in nats) of a Hermitian matrix: `∑_i -λ_i ln λ_i`. -/
def vonNeumann {ι : Type*} [Fintype ι] [DecidableEq ι] (ρ : Matrix ι ι ℂ) (hρ : ρ.IsHermitian) :
    ℝ :=
  ∑ i, Real.negMulLog (hρ.eigenvalues i)

namespace AreaLawAux

/-- Complexified matrices multiply like the real ones. -/
lemma cplx_mul' {ι : Type*} [Fintype ι] [DecidableEq ι] (A B : Matrix ι ι ℝ) :
    cplx (A * B) = cplx A * cplx B := by
  ext i j; simp [cplx, Matrix.mul_apply]

lemma cplx_conjTranspose' {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℝ) :
    (cplx A)ᴴ = cplx Aᵀ := by
  ext i j; simp [cplx, Matrix.conjTranspose_apply]

lemma cplx_ann_conjTranspose (x : Fin N) : (cplx (ann x))ᴴ = cplx (cre x) := by
  rw [cplx_conjTranspose']; rfl

lemma cplx_cre_conjTranspose (x : Fin N) : (cplx (cre x))ᴴ = cplx (ann x) := by
  rw [cplx_conjTranspose']; rfl

/-- CAR for creation operators. -/
lemma cre_mul_cre_add (x y : Fin N) : cre x * cre y + cre y * cre x = 0 := by
  have h := congrArg Matrix.transpose (ann_mul_ann_add y x)
  simpa [Matrix.transpose_add, Matrix.transpose_mul, cre, add_comm] using h

/-- The complex creation operator `c†(w) = ∑_x w_x c_x†`. -/
def crOp (w : Fin N → ℂ) : Matrix (Cfg N) (Cfg N) ℂ := ∑ x, w x • cplx (cre x)

/-- `c†(w)` as a linear map. -/
def crLin : (Fin N → ℂ) →ₗ[ℂ] Matrix (Cfg N) (Cfg N) ℂ where
  toFun := crOp
  map_add' v w := by simp [crOp, add_smul, Finset.sum_add_distrib]
  map_smul' c w := by simp [crOp, Finset.smul_sum, smul_smul]

lemma crLin_apply (w : Fin N → ℂ) : crLin w = crOp w := rfl

/-- `c†(w)² = 0`. -/
lemma crOp_mul_self (w : Fin N → ℂ) : crOp w * crOp w = 0 := by
  have key : ∀ x y : Fin N, cplx (cre x) * cplx (cre y) = -(cplx (cre y) * cplx (cre x)) := by
    intro x y
    rw [← cplx_mul', ← cplx_mul', eq_neg_iff_add_eq_zero]
    have := cre_mul_cre_add x y
    have h2 : cplx (cre x * cre y + cre y * cre x) = 0 := by rw [this]; simp [cplx]
    simpa [cplx, Matrix.map_add] using h2
  have hS : crOp w * crOp w = ∑ x, ∑ y, (w x * w y) • (cplx (cre x) * cplx (cre y)) := by
    simp only [crOp, Finset.sum_mul, Finset.mul_sum, smul_mul_smul_comm]
    exact Finset.sum_comm
  have hS' : crOp w * crOp w = -(crOp w * crOp w) := by
    calc crOp w * crOp w = ∑ x, ∑ y, (w x * w y) • (cplx (cre x) * cplx (cre y)) := hS
    _ = ∑ y, ∑ x, (w x * w y) • (cplx (cre x) * cplx (cre y)) := Finset.sum_comm
    _ = ∑ y, ∑ x, -((w y * w x) • (cplx (cre y) * cplx (cre x))) := by
        refine Finset.sum_congr rfl fun y _ => Finset.sum_congr rfl fun x _ => ?_
        rw [key x y, smul_neg, mul_comm]
    _ = -(crOp w * crOp w) := by rw [hS]; simp [Finset.sum_neg_distrib]
  have h2 : (2 : ℂ) • (crOp w * crOp w) = 0 := by
    rw [two_smul]; nth_rewrite 2 [hS']; simp
  exact (smul_eq_zero.mp h2).resolve_left two_ne_zero

/-- The algebra map from the exterior algebra induced by `c†`. -/
def crAlg : ExteriorAlgebra ℂ (Fin N → ℂ) →ₐ[ℂ] Matrix (Cfg N) (Cfg N) ℂ :=
  ExteriorAlgebra.lift ℂ ⟨crLin, fun w => crOp_mul_self w⟩

/-- Applying a matrix to the vacuum, as a linear map. -/
def vacLin : Matrix (Cfg N) (Cfg N) ℂ →ₗ[ℂ] (Cfg N → ℂ) where
  toFun X := X *ᵥ fockVec ∅
  map_add' X Y := Matrix.add_mulVec X Y _
  map_smul' c X := Matrix.smul_mulVec c X _

/-- The Slater determinant `c†(w_0) ⋯ c†(w_{M-1}) |0⟩` as an alternating map. -/
def slaterAlt (N M : ℕ) : (Fin N → ℂ) [⋀^Fin M]→ₗ[ℂ] (Cfg N → ℂ) :=
  ((vacLin (N := N)).comp (crAlg (N := N)).toLinearMap).compAlternatingMap
    (ExteriorAlgebra.ιMulti ℂ M)

lemma slaterAlt_apply {M : ℕ} (w : Fin M → Fin N → ℂ) :
    slaterAlt N M w = (List.ofFn fun k => crOp (w k)).prod *ᵥ fockVec ∅ := by
  simp only [slaterAlt, LinearMap.compAlternatingMap_apply, LinearMap.comp_apply,
    AlgHom.toLinearMap_apply, ExteriorAlgebra.ιMulti_apply, map_list_prod, List.map_ofFn]
  rfl


/-- Change of orbital basis multiplies a Slater determinant by the determinant. -/
lemma slaterAlt_mul {M : ℕ} (Φ : Matrix (Fin N) (Fin M) ℂ) (V : Matrix (Fin M) (Fin M) ℂ) :
    slaterAlt N M (fun k x => (Φ * V) x k) = V.det • slaterAlt N M (fun k x => Φ x k) := by
  funext n
  let g : (Fin M → ℂ) [⋀^Fin M]→ₗ[ℂ] ℂ :=
    ((LinearMap.proj n : (Cfg N → ℂ) →ₗ[ℂ] ℂ).compAlternatingMap (slaterAlt N M)).compLinearMap
      Φ.mulVecLin
  have hg := AlternatingMap.eq_smul_basis_det (Pi.basisFun ℂ (Fin M)) g
  have h1 := congrArg (fun f => f (fun k m => V m k)) hg
  simp only [smul_apply, Pi.basisFun_det_apply, smul_eq_mul] at h1
  have e1 : (fun k x => (Φ * V) x k) = fun k => Φ.mulVecLin (fun m => V m k) := by
    funext k x; simp [Matrix.mulVec, Matrix.mul_apply, dotProduct]
  have e2 : (fun k x => Φ x k) = fun k => Φ.mulVecLin (Pi.basisFun ℂ (Fin M) k) := by
    funext k x; simp [Matrix.mulVec, dotProduct, Pi.single_apply]
  have e3 : (Matrix.of fun k m => V m k).det = V.det := by
    rw [← Matrix.det_transpose]; rfl
  rw [Pi.smul_apply, smul_eq_mul, e1, e2, ← e3, mul_comm]
  exact h1


/-! ### Action of creation and annihilation operators on vectors -/

lemma cre_apply' (x : Fin N) (n m : Cfg N) :
    cre x n m = if n x = 1 ∧ m = Function.update n x 0 then sgn (jwPar n x) else 0 := rfl

/-- `(c_x† ψ)(n) = [n_x = 1] (-1)^{jw} ψ(n - e_x)`. -/
lemma cre_mulVec_apply (x : Fin N) (ψ : Cfg N → ℂ) (n : Cfg N) :
    (cplx (cre x) *ᵥ ψ) n =
      if n x = 1 then (sgn (jwPar n x) : ℂ) * ψ (Function.update n x 0) else 0 := by
  simp only [Matrix.mulVec, dotProduct, cplx, Matrix.map_apply, cre_apply']
  split_ifs with h
  · rw [Finset.sum_eq_single (Function.update n x 0)]
    · simp [h]
    · intro m _ hm; simp [hm]
    · simp
  · simp [h]

/-- Annihilating a particle at an occupied site. -/
lemma ann_mulVec_update (x : Fin N) (ψ : Cfg N → ℂ) (n : Cfg N) (hn : n x = 1) :
    (cplx (ann x) *ᵥ ψ) (Function.update n x 0) = (sgn (jwPar n x) : ℂ) * ψ n := by
  simp only [Matrix.mulVec, dotProduct, cplx, Matrix.map_apply, ann]
  rw [Finset.sum_eq_single n]
  · simp [hn]
  · intro m _ hm
    have hc : ¬ (m x = 1 ∧ Function.update n x 0 = Function.update m x 0) := by
      rintro ⟨h1, h2⟩
      apply hm
      funext k
      by_cases hk : k = x
      · subst hk; rw [h1, hn]
      · have := congrFun h2 k
        simpa [Function.update_apply, hk] using this.symm
    simp [hc]
  · simp

/-- A vector annihilated by all `c_x` is a multiple of the vacuum. -/
lemma eq_smul_vac_of_ann (ψ : Cfg N → ℂ) (h : ∀ x, cplx (ann x) *ᵥ ψ = 0) :
    ψ = ψ (fockCfg ∅) • fockVec ∅ := by
  funext n
  by_cases hn : n = fockCfg ∅
  · subst hn; simp [fockVec]
  · have : ∃ x, n x = 1 := by
      by_contra hcon
      push Not at hcon
      apply hn
      funext x
      have := hcon x
      simp only [fockCfg, Finset.notMem_empty, ite_false]
      revert this; generalize n x = a; revert a; decide
    obtain ⟨x, hx⟩ := this
    have h1 := congrFun (h x) (Function.update n x 0)
    rw [ann_mulVec_update x ψ n hx] at h1
    have hs : (sgn (jwPar n x) : ℂ) ≠ 0 := by
      unfold sgn; split_ifs <;> norm_num
    have : ψ n = 0 := by simpa [hs] using h1
    simp [this, fockVec, hn]

/-- `c_y |0⟩ = 0`. -/
lemma ann_mulVec_vac (y : Fin N) : cplx (ann y) *ᵥ fockVec ∅ = 0 := by
  funext n
  simp [fockVec, Matrix.mulVec_single, cplx, ann, fockCfg]


/-! ### Norms -/

/-- Squared `ℓ²` norm of a vector. -/
def nsq {α : Type*} [Fintype α] (v : α → ℂ) : ℝ := ∑ a, ‖v a‖ ^ 2

/-- `ℓ²` norm of a vector. -/
def l2 {α : Type*} [Fintype α] (v : α → ℂ) : ℝ := √(nsq v)

lemma nsq_nonneg {α : Type*} [Fintype α] (v : α → ℂ) : 0 ≤ nsq v :=
  Finset.sum_nonneg fun _ _ => by positivity

lemma l2_nonneg {α : Type*} [Fintype α] (v : α → ℂ) : 0 ≤ l2 v := Real.sqrt_nonneg _

lemma star_dotProduct_self {α : Type*} [Fintype α] (v : α → ℂ) :
    star v ⬝ᵥ v = ((nsq v : ℝ) : ℂ) := by
  simp only [dotProduct, Pi.star_apply, nsq, Complex.ofReal_sum, Complex.ofReal_pow]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [Complex.star_def, Complex.conj_mul']

lemma star_mulVec_dotProduct {α : Type*} [Fintype α] (M : Matrix α α ℂ) (v : α → ℂ) :
    star (M *ᵥ v) ⬝ᵥ (M *ᵥ v) = star v ⬝ᵥ ((Mᴴ * M) *ᵥ v) := by
  rw [Matrix.star_mulVec, Matrix.dotProduct_mulVec (star v ᵥ* Mᴴ) M v, Matrix.vecMul_vecMul,
    Matrix.dotProduct_mulVec (star v) (Mᴴ * M) v]

lemma nsq_mulVec_unitary {α : Type*} [Fintype α] [DecidableEq α] (W : Matrix α α ℂ) (hW : Wᴴ * W = 1)
    (v : α → ℂ) : nsq (W *ᵥ v) = nsq v := by
  have h := star_mulVec_dotProduct W v
  rw [hW, Matrix.one_mulVec, star_dotProduct_self, star_dotProduct_self] at h
  exact_mod_cast h

lemma nsq_smul {α : Type*} [Fintype α] (c : ℂ) (v : α → ℂ) : nsq (c • v) = ‖c‖ ^ 2 * nsq v := by
  simp [nsq, Finset.mul_sum, mul_pow]

lemma l2_smul {α : Type*} [Fintype α] (c : ℂ) (v : α → ℂ) : l2 (c • v) = ‖c‖ * l2 v := by
  rw [l2, nsq_smul, Real.sqrt_mul (by positivity), Real.sqrt_sq (norm_nonneg _), l2]

/-! ### The CAR bound `‖c†(w) ψ‖ ≤ ‖w‖ ‖ψ‖` -/

lemma crOp_conjTranspose (w : Fin N → ℂ) : (crOp w)ᴴ = ∑ x, star (w x) • cplx (ann x) := by
  simp [crOp, Matrix.conjTranspose_sum, Matrix.conjTranspose_smul, cplx_cre_conjTranspose]

lemma cplx_one' {ι : Type*} [Fintype ι] [DecidableEq ι] : cplx (1 : Matrix ι ι ℝ) = 1 := by
  simp [cplx]

lemma crOp_CAR (w : Fin N → ℂ) :
    (crOp w)ᴴ * crOp w + crOp w * (crOp w)ᴴ =
      ((nsq w : ℝ) : ℂ) • (1 : Matrix (Cfg N) (Cfg N) ℂ) := by
  have hcar : ∀ x y : Fin N, cplx (ann x) * cplx (cre y) + cplx (cre y) * cplx (ann x) =
      if x = y then 1 else 0 := by
    intro x y
    rw [← cplx_mul', ← cplx_mul']
    have h := congrArg cplx (ann_mul_cre_add x y)
    have hadd : cplx (ann x * cre y + cre y * ann x) =
        cplx (ann x * cre y) + cplx (cre y * ann x) := by simp [cplx, Matrix.map_add]
    rw [← hadd, h]
    split_ifs
    · exact cplx_one'
    · simp [cplx]
  rw [crOp_conjTranspose]
  calc (∑ x, star (w x) • cplx (ann x)) * crOp w + crOp w * (∑ x, star (w x) • cplx (ann x))
      = ∑ x, ∑ y, ((star (w x) * w y) • (cplx (ann x) * cplx (cre y)) +
          (w y * star (w x)) • (cplx (cre y) * cplx (ann x))) := by
        simp only [crOp, Finset.sum_mul_sum, smul_mul_smul_comm, Finset.sum_add_distrib]
        rw [Finset.sum_comm (f := fun x y => (w x * star (w y)) • (cplx (cre x) * cplx (ann y)))]
    _ = ∑ x, ∑ y, (star (w x) * w y) • (if x = y then (1 : Matrix (Cfg N) (Cfg N) ℂ) else 0) := by
        refine Finset.sum_congr rfl fun x _ => Finset.sum_congr rfl fun y _ => ?_
        rw [mul_comm (w y), ← smul_add, hcar]
    _ = ((nsq w : ℝ) : ℂ) • 1 := by
        simp_rw [smul_ite, smul_zero, Finset.sum_ite_eq, Finset.mem_univ, ite_true,
          ← Finset.sum_smul]
        simp only [nsq, Complex.ofReal_sum, Complex.ofReal_pow]
        refine congrArg (· • (1 : Matrix (Cfg N) (Cfg N) ℂ)) (Finset.sum_congr rfl fun a _ => ?_)
        rw [Complex.star_def, Complex.conj_mul']

lemma nsq_crOp_mulVec_le (w : Fin N → ℂ) (ψ : Cfg N → ℂ) :
    nsq (crOp w *ᵥ ψ) ≤ nsq w * nsq ψ := by
  have h1 := star_mulVec_dotProduct (crOp w) ψ
  have h2 := star_mulVec_dotProduct (crOp w)ᴴ ψ
  rw [Matrix.conjTranspose_conjTranspose] at h2
  have h3 : star (crOp w *ᵥ ψ) ⬝ᵥ (crOp w *ᵥ ψ) + star ((crOp w)ᴴ *ᵥ ψ) ⬝ᵥ ((crOp w)ᴴ *ᵥ ψ) =
      ((nsq w : ℝ) : ℂ) * (star ψ ⬝ᵥ ψ) := by
    rw [h1, h2, ← dotProduct_add, ← Matrix.add_mulVec, crOp_CAR, Matrix.smul_mulVec,
      Matrix.one_mulVec, dotProduct_smul, smul_eq_mul]
  rw [star_dotProduct_self, star_dotProduct_self, star_dotProduct_self] at h3
  have h4 : nsq (crOp w *ᵥ ψ) + nsq ((crOp w)ᴴ *ᵥ ψ) = nsq w * nsq ψ := by exact_mod_cast h3
  linarith [nsq_nonneg ((crOp w)ᴴ *ᵥ ψ)]

lemma l2_crOp_mulVec_le (w : Fin N → ℂ) (ψ : Cfg N → ℂ) :
    l2 (crOp w *ᵥ ψ) ≤ l2 w * l2 ψ := by
  rw [l2, l2, l2, ← Real.sqrt_mul (nsq_nonneg w)]
  exact Real.sqrt_le_sqrt (nsq_crOp_mulVec_le w ψ)

/-! ### Heisenberg picture: the evolved Fock state is a Slater determinant -/

section Heis

variable (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (t : ℝ)
include hA

lemma ann_mul_evo (x : Fin N) :
    cplx (ann x) * evo (dGamma A) t = evo (dGamma A) t * ∑ y, evo A t x y • cplx (ann y) := by
  rw [← heisenberg_ann A hA x t]
  simp only [← mul_assoc, evo_unitary' _ (dGamma_isSymm A hA) t, one_mul]

lemma evo_mul_cre (z : Fin N) :
    evo (dGamma A) t * cplx (cre z) = crOp (fun x => evo A t x z) * evo (dGamma A) t := by
  have hW1 := evo_unitary _ (dGamma_isSymm A hA) t
  have hW2 := evo_unitary' _ (dGamma_isSymm A hA) t
  set W := evo (dGamma A) t with hW
  set U := evo A t with hUdef
  have hstar : ∀ x, Wᴴ * cplx (cre x) = (∑ y, star (U x y) • cplx (cre y)) * Wᴴ := by
    intro x
    have h := congrArg Matrix.conjTranspose (ann_mul_evo A hA t x)
    simpa [Matrix.conjTranspose_mul, Matrix.conjTranspose_sum, Matrix.conjTranspose_smul,
      cplx_ann_conjTranspose] using h
  have hU : ∀ y, ∑ x, U x z * star (U x y) = if y = z then 1 else 0 := by
    intro y
    have h := congrFun (congrFun (evo_unitary A hA t) y) z
    simp only [Matrix.mul_apply, Matrix.conjTranspose_apply, Matrix.one_apply] at h
    rw [← h]
    exact Finset.sum_congr rfl fun x _ => mul_comm _ _
  have inner : ∑ x, U x z • ∑ y, star (U x y) • cplx (cre y) = cplx (cre z) := by
    simp_rw [Finset.smul_sum, smul_smul]
    rw [Finset.sum_comm]
    simp_rw [← Finset.sum_smul, hU, ite_smul, one_smul, zero_smul]
    simp
  have key : Wᴴ * crOp (fun x => U x z) = cplx (cre z) * Wᴴ := by
    simp only [crOp, Finset.mul_sum, mul_smul_comm, hstar, ← smul_mul_assoc, ← Finset.sum_mul,
      inner]
  calc W * cplx (cre z) = W * cplx (cre z) * (Wᴴ * W) := by rw [hW1, mul_one]
  _ = W * (cplx (cre z) * Wᴴ) * W := by simp only [mul_assoc]
  _ = W * (Wᴴ * crOp (fun x => U x z)) * W := by rw [key]
  _ = (W * Wᴴ) * crOp (fun x => U x z) * W := by simp only [mul_assoc]
  _ = crOp (fun x => U x z) * W := by rw [hW2, one_mul]

lemma evo_mul_prod {M : ℕ} (z : Fin M → Fin N) :
    evo (dGamma A) t * (List.ofFn fun k => cplx (cre (z k))).prod =
      (List.ofFn fun k => crOp (fun x => evo A t x (z k))).prod * evo (dGamma A) t := by
  induction M with
  | zero => simp
  | succ M ih =>
    rw [List.ofFn_succ, List.ofFn_succ, List.prod_cons, List.prod_cons, ← mul_assoc,
      evo_mul_cre A hA t, mul_assoc, ih (fun k => z k.succ), mul_assoc]

lemma evo_vac : ∃ c : ℂ, ‖c‖ = 1 ∧ evo (dGamma A) t *ᵥ fockVec ∅ = c • fockVec ∅ := by
  set ψ := evo (dGamma A) t *ᵥ fockVec (∅ : Finset (Fin N)) with hψ
  have hann : ∀ x, cplx (ann x) *ᵥ ψ = 0 := by
    intro x
    rw [hψ, Matrix.mulVec_mulVec, ann_mul_evo A hA t x, ← Matrix.mulVec_mulVec,
      Matrix.sum_mulVec]
    simp [Matrix.smul_mulVec, ann_mulVec_vac]
  refine ⟨ψ (fockCfg ∅), ?_, eq_smul_vac_of_ann ψ hann⟩
  have h1 : nsq ψ = nsq (fockVec (∅ : Finset (Fin N))) :=
    nsq_mulVec_unitary _ (evo_unitary _ (dGamma_isSymm A hA) t) _
  have hv : nsq (fockVec (∅ : Finset (Fin N))) = 1 := by
    have : ∀ x : Cfg N, ‖(fockVec (∅ : Finset (Fin N))) x‖ ^ 2 =
        if x = fockCfg ∅ then 1 else 0 := by
      intro x; simp only [fockVec, Pi.single_apply]; split_ifs <;> simp
    simp [nsq, this]
  have h2 : nsq ψ = ‖ψ (fockCfg ∅)‖ ^ 2 * nsq (fockVec (∅ : Finset (Fin N))) := by
    conv_lhs => rw [eq_smul_vac_of_ann ψ hann]
    rw [nsq_smul]
  rw [h1, hv, mul_one] at h2
  exact (pow_eq_one_iff_of_nonneg (norm_nonneg _) two_ne_zero).mp h2.symm

end Heis

lemma fockVec_insert (a : Fin N) (s : Finset (Fin N)) (h : ∀ x ∈ s, a < x) :
    fockVec (insert a s) = cplx (cre a) *ᵥ fockVec s := by
  have ha : a ∉ s := fun h' => lt_irrefl a (h a h')
  funext n
  rw [cre_mulVec_apply]
  simp only [fockVec, Pi.single_apply]
  by_cases hn : n = fockCfg (insert a s)
  · subst hn
    have h1 : fockCfg (insert a s) a = 1 := by simp [fockCfg]
    have h2 : Function.update (fockCfg (insert a s)) a 0 = fockCfg s := by
      funext k
      by_cases hk : k = a
      · subst hk; simp [fockCfg, ha]
      · simp [hk, fockCfg]
    have h3 : jwPar (fockCfg (insert a s)) a = 0 := by
      unfold jwPar
      apply Finset.sum_eq_zero
      intro k hk
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hk
      have hka : k ≠ a := ne_of_lt hk
      have hks : k ∉ s := fun hks => lt_asymm hk (h k hks)
      simp [fockCfg, hka, hks]
    simp [h1, h2, h3, sgn]
  · rw [ite_eq_right hn]
    split_ifs with h1 h2
    · exfalso; apply hn
      funext k
      by_cases hk : k = a
      · subst hk; simp [fockCfg, h1]
      · have := congrFun h2 k
        simp only [Function.update_apply, hk, ite_false] at this
        simp [fockCfg, hk, this]
    · simp
    · rfl

lemma fockVec_eq_prod (S : Finset (Fin N)) :
    ∃ (M : ℕ) (z : Fin M → Fin N), Function.Injective z ∧ (∀ x, x ∈ S ↔ ∃ k, z k = x) ∧
      fockVec S = (List.ofFn fun k => cplx (cre (z k))).prod *ᵥ fockVec ∅ := by
  induction S using Finset.induction_on_min with
  | empty => exact ⟨0, Fin.elim0, fun i => i.elim0, by simp, by simp⟩
  | insert a s has ih =>
    obtain ⟨M, z, hz, hzS, hprod⟩ := ih
    have ha : a ∉ s := fun h' => lt_irrefl a (has a h')
    refine ⟨M + 1, Fin.cons a z, ?_, ?_, ?_⟩
    · rw [Fin.cons_injective_iff]
      refine ⟨?_, hz⟩
      rintro ⟨k, hk⟩
      exact ha ((hzS a).mpr ⟨k, hk⟩)
    · intro x
      simp only [Finset.mem_insert, Fin.exists_fin_succ, Fin.cons_zero, Fin.cons_succ, hzS]
      constructor
      · rintro (h | h)
        · exact Or.inl h.symm
        · exact Or.inr h
      · rintro (h | h)
        · exact Or.inl h.symm
        · exact Or.inr h
    · rw [List.ofFn_succ, List.prod_cons, ← Matrix.mulVec_mulVec]
      simp only [Fin.cons_zero, Fin.cons_succ]
      rw [← hprod, fockVec_insert a s has]

/-- **The evolved Fock state is a Slater determinant** (up to a phase). -/
lemma evo_fockVec_eq (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (t : ℝ) (S : Finset (Fin N)) :
    ∃ (M : ℕ) (z : Fin M → Fin N) (c : ℂ), Function.Injective z ∧
      (∀ x, x ∈ S ↔ ∃ k, z k = x) ∧ ‖c‖ = 1 ∧
      evo (dGamma A) t *ᵥ fockVec S = c • slaterAlt N M (fun k x => evo A t x (z k)) := by
  obtain ⟨M, z, hz, hzS, hS⟩ := fockVec_eq_prod S
  obtain ⟨c, hc, hvac⟩ := evo_vac A hA t
  refine ⟨M, z, c, hz, hzS, hc, ?_⟩
  rw [hS, Matrix.mulVec_mulVec, evo_mul_prod A hA t, ← Matrix.mulVec_mulVec, hvac,
    Matrix.mulVec_smul, slaterAlt_apply]

end AreaLawAux

open AreaLawAux

/-- **Entanglement of Slater determinants**: the entanglement entropy of a left segment is
bounded by `2 ∑_{x ∈ A, y ∉ A} |C_{xy}|` (the paper has the sharper constant `2 ln 2`). -/
theorem vonNeumann_slater_le (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (rdm ℓ hℓ (evo (dGamma A) t *ᵥ fockVec S))
        (rdm_isHermitian ℓ hℓ _) ≤
      2 * ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ‖corrMat A S t x y‖ := by
  sorry

/-- **Area law for the charge-conditioned entanglement** `S_G` (paper: Theorem `thm:arealaw`):
under Hypothesis (EC), the Bernoulli average over charge sectors of the entanglement entropy of a
left segment is bounded uniformly in time, volume, segment length and occupied set. -/
theorem lgt_areaLaw (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ K : ℝ, ∀ (N : ℕ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (S : Finset (Fin N)) (t : ℝ),
      avgQ (fun q : Cfg N =>
        vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ fockVec S))
          (rdm_isHermitian ℓ hℓ _)) ≤ K := by
  sorry

end DFL
