import DFL.FreeDyn
import DFL.Entropy
import DFL.OneBody

/-!
# Entanglement area law (paper: Theorem `thm:arealaw`, part on `S_G`)

For a Fock state evolved by a quadratic Hamiltonian (a Slater determinant), the von Neumann
entropy of a left segment `A = {0, …, ℓ-1}` satisfies
`S(ρ_A) ≤ 2 ∑_k √(ν_k (1 - ν_k)) ≤ 2 ∑_{x ∈ A, y ∉ A} |C_{xy}|`,
where `ν_k` are the eigenvalues of the restricted correlation matrix `C_A`. (Classically
`S(ρ_A) = ∑_k h(ν_k)`; the formal proof does not compute the spectrum of `ρ_A` but bounds
`S(ρ_A) ≤ 2 log ∑_i √λ_i` by Jensen and `∑_i √λ_i ≤ ∏_k (√ν_k + √(1 - ν_k))` through a
product decomposition of the Slater determinant in canonical orbitals.) Combined with the duality,
the sector averaging and Hypothesis (EC), this gives a uniform-in-time area law for the
sector-averaged (charge-conditioned) entanglement `S_G` of the lattice gauge theory.

For a left segment the Jordan–Wigner occupation basis factorizes canonically,
`|n⟩ = |n_A⟩ ⊗ |n_B⟩`, so the fermionic reduced density matrix is the ordinary partial trace in this
basis.
-/

open Matrix Finset
open scoped ComplexOrder

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

/-!
### Proof strategy for `vonNeumann_slater_le`

1. (`evo_fockVec_eq`) By the Heisenberg relation `W c_z† = c†(U e_z) W` and `W |0⟩ = phase · |0⟩`,
   the evolved Fock state is (up to a phase) the Slater determinant
   `c†(u_1) ⋯ c†(u_M) |0⟩` of the orbitals `u_k = U e_{z_k}`.
2. (`slaterAlt`, `slaterAlt_mul`) Via the exterior algebra (`c†(w)² = 0`), Slater determinants
   form an alternating map, so a unitary change of orbitals only changes a phase.
3. (`canonical_orbitals`) Diagonalizing the overlap matrix `⟨u_k, P_A u_l⟩` gives orbitals with
   `⟨P_A u_k, P_A u_l⟩ = ν_k δ_kl` and `⟨P_B u_k, P_B u_l⟩ = (1 - ν_k) δ_kl`.
4. (`slater_decomp`) Writing `c†(u) = c†(P_A u) + c†(P_B u)` and using the factorization of the
   Jordan–Wigner basis for a left segment, the determinant is a sum of product vectors of total
   cross norm `≤ ∏_k (‖P_A u_k‖ + ‖P_B u_k‖) = ∏_k (√ν_k + √(1-ν_k))`.
5. (`vonNeumann_rdm_le`) For a unit vector `v = ∑_j X_j ⊗ Y_j`, the von Neumann entropy is at most
   the Rényi-`1/2` entropy `2 log ∑_i √λ_i`, and `∑_i √λ_i ≤ ∑_j ‖X_j‖ ‖Y_j‖` (nuclear norm bound).
6. `log(√ν + √(1-ν)) ≤ √(ν(1-ν))` and (`sum_sqrt_nu_le`) `∑_k √(ν_k(1-ν_k)) ≤ ∑_{x∈A,y∉A} |C_xy|`
   (Cauchy–Schwarz and Bessel for the normalized restricted canonical orbitals).
-/

/-- Complexified matrices multiply like the real ones. -/
lemma cplx_mul' {ι : Type*} [Fintype ι] [DecidableEq ι] (A B : Matrix ι ι ℝ) :
    cplx (A * B) = cplx A * cplx B := by
  ext i j; simp [cplx, Matrix.mul_apply]

/-- The conjugate transpose of a complexified real matrix is its transpose. -/
lemma cplx_conjTranspose' {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℝ) :
    (cplx A)ᴴ = cplx Aᵀ := by
  ext i j; simp [cplx, Matrix.conjTranspose_apply]

/-- `(c_x)† = c_x†` for the complexified operators. -/
lemma cplx_ann_conjTranspose (x : Fin N) : (cplx (ann x))ᴴ = cplx (cre x) := by
  rw [cplx_conjTranspose']; rfl

/-- `(c_x†)† = c_x` for the complexified operators. -/
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

/-- `slaterAlt` is the ordered product of creation operators applied to the vacuum. -/
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

/-- The entries of `c_x†`. -/
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

/-- The squared norm is nonnegative. -/
lemma nsq_nonneg {α : Type*} [Fintype α] (v : α → ℂ) : 0 ≤ nsq v :=
  Finset.sum_nonneg fun _ _ => by positivity

/-- The norm is nonnegative. -/
lemma l2_nonneg {α : Type*} [Fintype α] (v : α → ℂ) : 0 ≤ l2 v := Real.sqrt_nonneg _

/-- `⟨v, v⟩ = ‖v‖²`. -/
lemma star_dotProduct_self {α : Type*} [Fintype α] (v : α → ℂ) :
    star v ⬝ᵥ v = ((nsq v : ℝ) : ℂ) := by
  simp only [dotProduct, Pi.star_apply, nsq, Complex.ofReal_sum, Complex.ofReal_pow]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [Complex.star_def, Complex.conj_mul']

/-- `⟨Mv, Mv⟩ = ⟨v, M†M v⟩`. -/
lemma star_mulVec_dotProduct {α : Type*} [Fintype α] (M : Matrix α α ℂ) (v : α → ℂ) :
    star (M *ᵥ v) ⬝ᵥ (M *ᵥ v) = star v ⬝ᵥ ((Mᴴ * M) *ᵥ v) := by
  rw [Matrix.star_mulVec, Matrix.dotProduct_mulVec (star v ᵥ* Mᴴ) M v, Matrix.vecMul_vecMul,
    Matrix.dotProduct_mulVec (star v) (Mᴴ * M) v]

/-- Unitary matrices preserve the norm. -/
lemma nsq_mulVec_unitary {α : Type*} [Fintype α] [DecidableEq α] (W : Matrix α α ℂ) (hW : Wᴴ * W = 1)
    (v : α → ℂ) : nsq (W *ᵥ v) = nsq v := by
  have h := star_mulVec_dotProduct W v
  rw [hW, Matrix.one_mulVec, star_dotProduct_self, star_dotProduct_self] at h
  exact_mod_cast h

/-- `‖c v‖² = |c|² ‖v‖²`. -/
lemma nsq_smul {α : Type*} [Fintype α] (c : ℂ) (v : α → ℂ) : nsq (c • v) = ‖c‖ ^ 2 * nsq v := by
  simp [nsq, Finset.mul_sum, mul_pow]

/-- `‖c v‖ = |c| ‖v‖`. -/
lemma l2_smul {α : Type*} [Fintype α] (c : ℂ) (v : α → ℂ) : l2 (c • v) = ‖c‖ * l2 v := by
  rw [l2, nsq_smul, Real.sqrt_mul (by positivity), Real.sqrt_sq (norm_nonneg _), l2]

/-- A basis vector has norm one. -/
lemma nsq_pi_single {α : Type*} [Fintype α] [DecidableEq α] (c : α) :
    nsq (Pi.single c (1 : ℂ) : α → ℂ) = 1 := by
  have : ∀ y : α, ‖(Pi.single c (1 : ℂ) : α → ℂ) y‖ ^ 2 = if y = c then 1 else 0 := by
    intro y; simp only [Pi.single_apply]; split_ifs <;> simp
  simp [nsq, this]

/-! ### The CAR bound `‖c†(w) ψ‖ ≤ ‖w‖ ‖ψ‖` -/

/-- `c†(w)† = ∑_x w̄_x c_x`. -/
lemma crOp_conjTranspose (w : Fin N → ℂ) : (crOp w)ᴴ = ∑ x, star (w x) • cplx (ann x) := by
  simp [crOp, Matrix.conjTranspose_sum, Matrix.conjTranspose_smul, cplx_cre_conjTranspose]

/-- Complexification maps `1` to `1`. -/
lemma cplx_one' {ι : Type*} [Fintype ι] [DecidableEq ι] : cplx (1 : Matrix ι ι ℝ) = 1 := by
  simp [cplx]

/-- CAR for `c†(w)`: `c(w̄) c†(w) + c†(w) c(w̄) = ‖w‖²`. -/
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

/-- `‖c†(w) ψ‖² ≤ ‖w‖² ‖ψ‖²`. -/
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

/-- `‖c†(w) ψ‖ ≤ ‖w‖ ‖ψ‖`. -/
lemma l2_crOp_mulVec_le (w : Fin N → ℂ) (ψ : Cfg N → ℂ) :
    l2 (crOp w *ᵥ ψ) ≤ l2 w * l2 ψ := by
  rw [l2, l2, l2, ← Real.sqrt_mul (nsq_nonneg w)]
  exact Real.sqrt_le_sqrt (nsq_crOp_mulVec_le w ψ)

/-! ### Heisenberg picture: the evolved Fock state is a Slater determinant -/

section Heis

variable (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (t : ℝ)
include hA

/-- Heisenberg relation in the form `c_x W = W ∑_y U_{xy} c_y`. -/
lemma ann_mul_evo (x : Fin N) :
    cplx (ann x) * evo (dGamma A) t = evo (dGamma A) t * ∑ y, evo A t x y • cplx (ann y) := by
  rw [← heisenberg_ann A hA x t]
  simp only [← mul_assoc, evo_unitary' _ (dGamma_isSymm A hA) t, one_mul]

/-- `W c_z† = c†(U e_z) W`: the propagator maps creation operators to rotated ones. -/
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

/-- `W` intertwines products of creation operators with products of rotated ones. -/
lemma evo_mul_prod {M : ℕ} (z : Fin M → Fin N) :
    evo (dGamma A) t * (List.ofFn fun k => cplx (cre (z k))).prod =
      (List.ofFn fun k => crOp (fun x => evo A t x (z k))).prod * evo (dGamma A) t := by
  induction M with
  | zero => simp
  | succ M ih =>
    rw [List.ofFn_succ, List.ofFn_succ, List.prod_cons, List.prod_cons, ← mul_assoc,
      evo_mul_cre A hA t, mul_assoc, ih (fun k => z k.succ), mul_assoc]

/-- The propagator maps the vacuum to a phase times the vacuum. -/
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
  have hv : nsq (fockVec (∅ : Finset (Fin N))) = 1 := nsq_pi_single _
  have h2 : nsq ψ = ‖ψ (fockCfg ∅)‖ ^ 2 * nsq (fockVec (∅ : Finset (Fin N))) := by
    conv_lhs => rw [eq_smul_vac_of_ann ψ hann]
    rw [nsq_smul]
  rw [h1, hv, mul_one] at h2
  exact (pow_eq_one_iff_of_nonneg (norm_nonneg _) two_ne_zero).mp h2.symm

end Heis

/-- `|{a} ∪ s⟩ = c_a† |s⟩` when `a` lies to the left of every site of `s`. -/
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

/-- Every Fock vector is an ordered product of creation operators applied to the vacuum. -/
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

/-! ### Factorization of the Jordan–Wigner basis for a left segment -/

section Split

variable (ℓ : ℕ) (hℓ : ℓ ≤ N)

/-- The site `a` of the left segment, as a site of the chain. -/
def inL (a : Fin ℓ) : Fin N := ⟨a.val, lt_of_lt_of_le a.isLt hℓ⟩

/-- The site `b` of the right part, as a site of the chain. -/
def inR (b : Fin (N - ℓ)) : Fin N := ⟨ℓ + b.val, by omega⟩

/-- The left part of the split configuration. -/
lemma splitCfg_fst (n : Cfg N) (a : Fin ℓ) : (splitCfg ℓ hℓ n).1 a = n (inL ℓ hℓ a) := rfl

/-- The right part of the split configuration. -/
lemma splitCfg_snd (n : Cfg N) (b : Fin (N - ℓ)) : (splitCfg ℓ hℓ n).2 b = n (inR ℓ hℓ b) :=
  rfl

/-- Splitting a sum over the sites of the chain. -/
lemma sum_fin_split {M : Type*} [AddCommMonoid M] (f : Fin N → M) :
    ∑ x, f x = ∑ a : Fin ℓ, f (inL ℓ hℓ a) + ∑ b : Fin (N - ℓ), f (inR ℓ hℓ b) := by
  have h : N = ℓ + (N - ℓ) := by omega
  rw [← (finCongr h).symm.sum_comp f, Fin.sum_univ_add]
  rfl

/-- A sum over the left segment, reindexed by `Fin ℓ`. -/
lemma sum_segment {M : Type*} [AddCommMonoid M] (f : Fin N → M) :
    ∑ x ∈ segment ℓ, f x = ∑ a : Fin ℓ, f (inL ℓ hℓ a) := by
  rw [segment, Finset.sum_filter, sum_fin_split ℓ hℓ]
  have h1 : ∀ a : Fin ℓ, (inL ℓ hℓ a).val < ℓ := fun a => a.isLt
  have h2 : ∀ b : Fin (N - ℓ), ¬ (inR ℓ hℓ b).val < ℓ := fun b => by simp [inR]
  simp only [h1, h2, ite_true, ite_false, Finset.sum_const_zero, add_zero]

/-- A sum over the complement of the left segment, reindexed by `Fin (N - ℓ)`. -/
lemma sum_compl_segment {M : Type*} [AddCommMonoid M] (f : Fin N → M) :
    ∑ x ∈ univ \ segment ℓ, f x = ∑ b : Fin (N - ℓ), f (inR ℓ hℓ b) := by
  have hs : (univ \ segment ℓ : Finset (Fin N)) = univ.filter (fun x => ¬ x.val < ℓ) := by
    ext x; simp [segment]
  rw [hs, Finset.sum_filter, sum_fin_split ℓ hℓ]
  have h1 : ∀ a : Fin ℓ, (inL ℓ hℓ a).val < ℓ := fun a => a.isLt
  have h2 : ∀ b : Fin (N - ℓ), ¬ (inR ℓ hℓ b).val < ℓ := fun b => by simp [inR]
  simp only [h1, h2, not_true_eq_false, not_false_eq_true, ite_true, ite_false,
    Finset.sum_const_zero, zero_add]

/-- `inL` is injective. -/
lemma inL_injective : Function.Injective (inL ℓ hℓ) := by
  intro a a' h; simpa [inL, Fin.ext_iff] using h

/-- `inR` is injective. -/
lemma inR_injective : Function.Injective (inR ℓ hℓ) := by
  intro b b' h; simpa [inR, Fin.ext_iff] using h

/-- Left and right sites are distinct. -/
lemma inL_ne_inR (a : Fin ℓ) (b : Fin (N - ℓ)) : inL ℓ hℓ a ≠ inR ℓ hℓ b := by
  intro h; have := congrArg Fin.val h; simp [inL, inR] at this; omega

/-- Updating a left site only updates the left part of the split configuration. -/
lemma splitCfg_update_inL (n : Cfg N) (a : Fin ℓ) (v : ZMod 2) :
    splitCfg ℓ hℓ (Function.update n (inL ℓ hℓ a) v) =
      (Function.update (splitCfg ℓ hℓ n).1 a v, (splitCfg ℓ hℓ n).2) := by
  ext i
  · simp only [splitCfg_fst, Function.update_apply, (inL_injective ℓ hℓ).eq_iff]
  · simp only [splitCfg_snd, Function.update_apply, (inL_ne_inR ℓ hℓ a i).symm, ite_false]

/-- Updating a right site only updates the right part of the split configuration. -/
lemma splitCfg_update_inR (n : Cfg N) (b : Fin (N - ℓ)) (v : ZMod 2) :
    splitCfg ℓ hℓ (Function.update n (inR ℓ hℓ b) v) =
      ((splitCfg ℓ hℓ n).1, Function.update (splitCfg ℓ hℓ n).2 b v) := by
  ext i
  · simp only [splitCfg_fst, Function.update_apply, inL_ne_inR ℓ hℓ i b, ite_false]
  · simp only [splitCfg_snd, Function.update_apply, (inR_injective ℓ hℓ).eq_iff]

/-- The Jordan–Wigner string of a left site only sees the left part. -/
lemma jwPar_inL (n : Cfg N) (a : Fin ℓ) :
    jwPar n (inL ℓ hℓ a) = jwPar (splitCfg ℓ hℓ n).1 a := by
  unfold jwPar
  rw [Finset.sum_filter, Finset.sum_filter, sum_fin_split ℓ hℓ]
  have h2 : ∀ b, ¬ (inR ℓ hℓ b < inL ℓ hℓ a) := by
    intro b; simp only [Fin.lt_def, inR, inL]; omega
  have h1 : ∀ a', (inL ℓ hℓ a' < inL ℓ hℓ a) ↔ a' < a := by
    intro a'; simp only [Fin.lt_def, inL]
  simp only [h2, ite_false, Finset.sum_const_zero, add_zero, h1, splitCfg_fst]

/-- The Jordan–Wigner string of a right site is the left parity plus the right string. -/
lemma jwPar_inR (n : Cfg N) (b : Fin (N - ℓ)) :
    jwPar n (inR ℓ hℓ b) = (∑ i, (splitCfg ℓ hℓ n).1 i) + jwPar (splitCfg ℓ hℓ n).2 b := by
  unfold jwPar
  rw [Finset.sum_filter, Finset.sum_filter, sum_fin_split ℓ hℓ]
  have h1 : ∀ a, inL ℓ hℓ a < inR ℓ hℓ b := by
    intro a; simp only [Fin.lt_def, inR, inL]; omega
  have h2 : ∀ b', (inR ℓ hℓ b' < inR ℓ hℓ b) ↔ b' < b := by
    intro b'; simp only [Fin.lt_def, inR]; omega
  simp only [h1, ite_true, h2, splitCfg_fst, splitCfg_snd]

/-- `(-1)^{a+b} = (-1)^a (-1)^b`. -/
lemma sgn_add (a b : ZMod 2) : sgn (a + b) = sgn a * sgn b := by
  have h01 : ∀ c : ZMod 2, c = 0 ∨ c = 1 := by decide
  have h11 : (1 + 1 : ZMod 2) = 0 := by decide
  rcases h01 a with rfl | rfl <;> rcases h01 b with rfl | rfl <;> simp [sgn, h11]

/-- Every site is either a left site or a right site. -/
lemma fin_cases_split (x : Fin N) : (∃ a, inL ℓ hℓ a = x) ∨ (∃ b, inR ℓ hℓ b = x) := by
  by_cases h : x.val < ℓ
  · exact Or.inl ⟨⟨x.val, h⟩, Fin.ext rfl⟩
  · exact Or.inr ⟨⟨x.val - ℓ, by omega⟩, Fin.ext (by simp only [inR]; omega)⟩

/-- The product vector `X ⊗ Y` in the factorized occupation basis. -/
def tens (X : Cfg ℓ → ℂ) (Y : Cfg (N - ℓ) → ℂ) : Cfg N → ℂ :=
  fun n => X (splitCfg ℓ hℓ n).1 * Y (splitCfg ℓ hℓ n).2

/-- Parity sign `(-1)^{N}` of a configuration. -/
def parC {K : ℕ} (a : Cfg K) : ℂ := ((sgn (∑ i, a i) : ℝ) : ℂ)

/-- A left creation operator acts on the left factor of a product vector. -/
lemma cre_tens_inL (a : Fin ℓ) (X : Cfg ℓ → ℂ) (Y : Cfg (N - ℓ) → ℂ) :
    cplx (cre (inL ℓ hℓ a)) *ᵥ tens ℓ hℓ X Y = tens ℓ hℓ (cplx (cre a) *ᵥ X) Y := by
  funext n
  rw [cre_mulVec_apply]
  simp only [tens, cre_mulVec_apply, splitCfg_update_inL, jwPar_inL, splitCfg_fst]
  by_cases h : n (inL ℓ hℓ a) = 1
  · simp only [h, ite_true]; ring
  · simp only [h, ite_false, zero_mul]

/-- A right creation operator acts on the right factor, with the parity of the left factor. -/
lemma cre_tens_inR (b : Fin (N - ℓ)) (X : Cfg ℓ → ℂ) (Y : Cfg (N - ℓ) → ℂ) :
    cplx (cre (inR ℓ hℓ b)) *ᵥ tens ℓ hℓ X Y =
      tens ℓ hℓ (fun a => parC a * X a) (cplx (cre b) *ᵥ Y) := by
  funext n
  rw [cre_mulVec_apply]
  simp only [tens, cre_mulVec_apply, splitCfg_update_inR, jwPar_inR, splitCfg_snd, parC,
    sgn_add, Complex.ofReal_mul]
  by_cases h : n (inR ℓ hℓ b) = 1
  · simp only [h, ite_true]; ring
  · simp only [h, ite_false, mul_zero]

/-- `⊗` is linear in the left factor. -/
lemma tens_sum_smul_left {ι : Type*} (s : Finset ι) (c : ι → ℂ) (X : ι → Cfg ℓ → ℂ)
    (Y : Cfg (N - ℓ) → ℂ) :
    tens ℓ hℓ (∑ i ∈ s, c i • X i) Y = ∑ i ∈ s, c i • tens ℓ hℓ (X i) Y := by
  funext n; simp [tens, Finset.sum_apply, Finset.sum_mul, mul_assoc]

/-- `⊗` is linear in the right factor. -/
lemma tens_sum_smul_right {ι : Type*} (s : Finset ι) (c : ι → ℂ) (X : Cfg ℓ → ℂ)
    (Y : ι → Cfg (N - ℓ) → ℂ) :
    tens ℓ hℓ X (∑ i ∈ s, c i • Y i) = ∑ i ∈ s, c i • tens ℓ hℓ X (Y i) := by
  funext n; simp only [tens, Finset.sum_apply, Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
  exact Finset.sum_congr rfl fun i _ => by ring

/-- `c†(v) Z = ∑_x v_x c_x† Z`. -/
lemma crOp_mulVec_eq {K : ℕ} (v : Fin K → ℂ) (Z : Cfg K → ℂ) :
    crOp v *ᵥ Z = ∑ x, v x • (cplx (cre x) *ᵥ Z) := by
  simp [crOp, Matrix.sum_mulVec, Matrix.smul_mulVec]

/-- `c†(w)` acting on a product vector. -/
lemma crOp_tens (w : Fin N → ℂ) (X : Cfg ℓ → ℂ) (Y : Cfg (N - ℓ) → ℂ) :
    crOp w *ᵥ tens ℓ hℓ X Y =
      tens ℓ hℓ (crOp (fun a => w (inL ℓ hℓ a)) *ᵥ X) Y +
        tens ℓ hℓ (fun a => parC a * X a) (crOp (fun b => w (inR ℓ hℓ b)) *ᵥ Y) := by
  rw [crOp_mulVec_eq, sum_fin_split ℓ hℓ, crOp_mulVec_eq, crOp_mulVec_eq, tens_sum_smul_left,
    tens_sum_smul_right]
  simp only [cre_tens_inL, cre_tens_inR]

/-- The vacuum is the product of the vacua of the two parts. -/
lemma tens_vac : tens ℓ hℓ (fockVec ∅) (fockVec ∅) = fockVec (∅ : Finset (Fin N)) := by
  funext n
  have h0 : ∀ {K : ℕ}, fockCfg (∅ : Finset (Fin K)) = 0 := by
    intro K; funext x; simp [fockCfg]
  simp only [tens, fockVec, Pi.single_apply, h0]
  have : n = 0 ↔ (splitCfg ℓ hℓ n).1 = 0 ∧ (splitCfg ℓ hℓ n).2 = 0 := by
    constructor
    · rintro rfl; exact ⟨rfl, rfl⟩
    · rintro ⟨h1, h2⟩
      funext x
      rcases (fin_cases_split ℓ hℓ x) with ⟨a, rfl⟩ | ⟨b, rfl⟩
      · exact congrFun h1 a
      · exact congrFun h2 b
  by_cases hn : n = 0
  · subst hn
    have h00 : splitCfg ℓ hℓ (0 : Cfg N) = (0, 0) := by ext <;> rfl
    simp [h00]
  · rw [ite_eq_right hn]
    rw [this] at hn
    by_cases h1 : (splitCfg ℓ hℓ n).1 = 0
    · simp [h1, show ¬ (splitCfg ℓ hℓ n).2 = 0 from fun h2 => hn ⟨h1, h2⟩]
    · simp [h1]

/-- Multiplying by the parity sign does not change the norm. -/
lemma l2_parC_mul {K : ℕ} (X : Cfg K → ℂ) : l2 (fun a => parC a * X a) = l2 X := by
  have h : ∀ a : Cfg K, ‖parC a‖ = 1 := by
    intro a; simp only [parC, sgn]; split_ifs <;> simp
  simp [l2, nsq, h]

/-- The vacuum has norm one. -/
lemma l2_fockVec_empty {K : ℕ} : l2 (fockVec (∅ : Finset (Fin K))) = 1 := by
  rw [l2, fockVec, nsq_pi_single, Real.sqrt_one]

/-- **Product decomposition of a Slater determinant** across the cut: the determinant is a sum
of product vectors whose total cross norm is at most `∏_k (‖P_A w_k‖ + ‖P_B w_k‖)`. -/
lemma slater_decomp {M : ℕ} (ws : Fin M → Fin N → ℂ) :
    ∃ L : List ((Cfg ℓ → ℂ) × (Cfg (N - ℓ) → ℂ)),
      slaterAlt N M ws = (L.map fun p => tens ℓ hℓ p.1 p.2).sum ∧
      (L.map fun p => l2 p.1 * l2 p.2).sum ≤
        ∏ k, (l2 (fun a => ws k (inL ℓ hℓ a)) + l2 (fun b => ws k (inR ℓ hℓ b))) := by
  induction M with
  | zero =>
    refine ⟨[(fockVec ∅, fockVec ∅)], ?_, ?_⟩
    · rw [slaterAlt_apply]
      simp [tens_vac]
    · simp [l2_fockVec_empty]
  | succ M ih =>
    obtain ⟨L, hL, hn⟩ := ih (fun k => ws k.succ)
    set wA : Fin ℓ → ℂ := fun a => ws 0 (inL ℓ hℓ a) with hwA
    set wB : Fin (N - ℓ) → ℂ := fun b => ws 0 (inR ℓ hℓ b) with hwB
    refine ⟨L.map (fun p => (crOp wA *ᵥ p.1, p.2)) ++
      L.map (fun p => ((fun a => parC a * p.1 a), crOp wB *ᵥ p.2)), ?_, ?_⟩
    · have h1 : slaterAlt N (M + 1) ws = crOp (ws 0) *ᵥ slaterAlt N M (fun k => ws k.succ) := by
        rw [slaterAlt_apply, slaterAlt_apply, List.ofFn_succ, List.prod_cons,
          ← Matrix.mulVec_mulVec]
      rw [h1, hL]
      have h2 : crOp (ws 0) *ᵥ (L.map fun p => tens ℓ hℓ p.1 p.2).sum =
          (L.map fun p => crOp (ws 0) *ᵥ tens ℓ hℓ p.1 p.2).sum := by
        rw [← Matrix.mulVecLin_apply, map_list_sum, List.map_map]
        rfl
      rw [h2]
      simp only [crOp_tens, List.sum_map_add, List.map_append, List.sum_append, List.map_map]
      rfl
    · simp only [List.map_append, List.sum_append, List.map_map, Function.comp_def]
      have hA : (L.map fun p => l2 (crOp wA *ᵥ p.1) * l2 p.2).sum ≤
          (L.map fun p => l2 wA * (l2 p.1 * l2 p.2)).sum := by
        apply List.sum_le_sum
        intro p _
        rw [← mul_assoc]
        exact mul_le_mul_of_nonneg_right (l2_crOp_mulVec_le wA p.1) (l2_nonneg _)
      have hB : (L.map fun p => l2 (fun a => parC a * p.1 a) * l2 (crOp wB *ᵥ p.2)).sum ≤
          (L.map fun p => l2 wB * (l2 p.1 * l2 p.2)).sum := by
        apply List.sum_le_sum
        intro p _
        rw [l2_parC_mul, mul_left_comm]
        exact mul_le_mul_of_nonneg_left (l2_crOp_mulVec_le wB p.2) (l2_nonneg _)
      rw [List.sum_map_mul_left] at hA hB
      rw [Fin.prod_univ_succ]
      have hpos : 0 ≤ l2 wA + l2 wB := add_nonneg (l2_nonneg _) (l2_nonneg _)
      calc _ ≤ l2 wA * (L.map fun p => l2 p.1 * l2 p.2).sum +
            l2 wB * (L.map fun p => l2 p.1 * l2 p.2).sum := add_le_add hA hB
        _ = (l2 wA + l2 wB) * (L.map fun p => l2 p.1 * l2 p.2).sum := by ring
        _ ≤ _ := mul_le_mul_of_nonneg_left hn hpos

end Split

/-! ### Entropy bounds -/

/-- Rényi-`1/2` bound: `∑ -p log p ≤ 2 log ∑ √p` for a probability vector `p`. -/
lemma sum_negMulLog_le_two_log_sum_sqrt {ι : Type*} [Fintype ι] (p : ι → ℝ)
    (h0 : ∀ i, 0 ≤ p i) (h1 : ∑ i, p i = 1) :
    ∑ i, Real.negMulLog (p i) ≤ 2 * Real.log (∑ i, √(p i)) := by
  set q : ι → ℝ := fun i => if 0 < p i then (√(p i))⁻¹ else 1 with hq
  have hmem : ∀ i ∈ (univ : Finset ι), q i ∈ Set.Ioi (0 : ℝ) := by
    intro i _
    simp only [hq, Set.mem_Ioi]
    split_ifs with h
    · exact inv_pos.mpr (Real.sqrt_pos.mpr h)
    · exact one_pos
  have hJ := strictConcaveOn_log_Ioi.concaveOn.le_map_sum (t := univ) (w := p) (p := q)
    (fun i _ => h0 i) h1 hmem
  have e1 : ∑ i, p i • q i = ∑ i, √(p i) := by
    refine Finset.sum_congr rfl fun i _ => ?_
    simp only [hq, smul_eq_mul]
    split_ifs with h
    · have hs : 0 < √(p i) := Real.sqrt_pos.mpr h
      rw [mul_inv_eq_iff_eq_mul₀ hs.ne', Real.mul_self_sqrt (h0 i)]
    · have : p i = 0 := le_antisymm (not_lt.mp h) (h0 i)
      simp [this]
  have e2 : ∑ i, p i • Real.log (q i) = (1 / 2) * ∑ i, Real.negMulLog (p i) := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    simp only [hq, smul_eq_mul, Real.negMulLog]
    split_ifs with h
    · rw [Real.log_inv, Real.log_sqrt (h0 i)]; ring
    · have : p i = 0 := le_antisymm (not_lt.mp h) (h0 i)
      simp [this]
  rw [e1, e2] at hJ
  linarith

/-- Bessel's inequality for a family of vectors that is orthonormal on `K`. -/
lemma bessel_on {ι κ : Type*} [Fintype ι] [Fintype κ] [DecidableEq κ] (u : κ → ι → ℂ)
    (K : Finset κ) (horth : ∀ k ∈ K, ∀ l ∈ K, star (u k) ⬝ᵥ u l = if k = l then 1 else 0)
    (w : ι → ℂ) : ∑ k ∈ K, ‖star (u k) ⬝ᵥ w‖ ^ 2 ≤ nsq w := by
  let v : K → EuclideanSpace ℂ ι := fun k => WithLp.toLp 2 (u k.1)
  have hv : Orthonormal ℂ v := by
    rw [orthonormal_iff_ite]
    intro k l
    simp only [v, EuclideanSpace.inner_toLp_toLp]
    rw [dotProduct_comm, horth k.1 k.2 l.1 l.2]
    simp only [Subtype.ext_iff]
  have hb := hv.sum_inner_products_le (WithLp.toLp 2 w) (s := Finset.univ)
  rw [EuclideanSpace.norm_sq_eq] at hb
  simp only [v, EuclideanSpace.inner_toLp_toLp] at hb
  rw [← Finset.sum_coe_sort K, nsq]
  refine le_of_eq_of_le (Finset.sum_congr rfl fun k _ => ?_) hb
  rw [dotProduct_comm]

/-- `⟨Mv, Mw⟩ = ⟨v, M†M w⟩` for rectangular `M`. -/
lemma star_mulVec_dotProduct_mulVec {α β : Type*} [Fintype α] [Fintype β] (M : Matrix β α ℂ)
    (v w : α → ℂ) : star (M *ᵥ v) ⬝ᵥ (M *ᵥ w) = star v ⬝ᵥ ((Mᴴ * M) *ᵥ w) := by
  rw [Matrix.star_mulVec, Matrix.dotProduct_mulVec (star v ᵥ* Mᴴ) M w, Matrix.vecMul_vecMul,
    Matrix.dotProduct_mulVec (star v) (Mᴴ * M) w]

/-- `|⟨x, y⟩| = |⟨y, x⟩|`. -/
lemma norm_star_dotProduct_comm {α : Type*} [Fintype α] (x y : α → ℂ) :
    ‖star x ⬝ᵥ y‖ = ‖star y ⬝ᵥ x‖ := by
  rw [← norm_star (star x ⬝ᵥ y)]
  congr 1
  simp only [dotProduct, star_sum, star_mul', Pi.star_apply, star_star]
  exact Finset.sum_congr rfl fun a _ => mul_comm _ _

/-- Complex conjugation preserves the norm. -/
lemma nsq_star {α : Type*} [Fintype α] (x : α → ℂ) : nsq (star x) = nsq x := by
  simp [nsq]

/-- **Nuclear-norm bound**: if `Ψ = ∑_j X_j Y_jᵀ` then `∑_i √λ_i(ΨΨ†) ≤ ∑_j ‖X_j‖ ‖Y_j‖`. -/
lemma sum_sqrt_eigenvalues_le {α β J : Type*} [Fintype α] [DecidableEq α] [Fintype β]
    [Fintype J] (ρ : Matrix α α ℂ) (hρ : ρ.IsHermitian) (Ψ : Matrix α β ℂ) (hρΨ : ρ = Ψ * Ψᴴ)
    (X : J → α → ℂ) (Y : J → β → ℂ) (hΨ : ∀ a b, Ψ a b = ∑ j, X j a * Y j b) :
    ∑ i, √(hρ.eigenvalues i) ≤ ∑ j, l2 (X j) * l2 (Y j) := by
  classical
  set lam := hρ.eigenvalues with hlam
  have hpsd : ρ.PosSemidef := hρΨ ▸ Matrix.posSemidef_self_mul_conjTranspose Ψ
  have hlam0 : ∀ i, 0 ≤ lam i := fun i => hpsd.eigenvalues_nonneg i
  set p : α → α → ℂ := fun i => (hρ.eigenvectorBasis i).ofLp with hp
  have hp_orth : ∀ i k, star (p i) ⬝ᵥ p k = if i = k then 1 else 0 := by
    intro i k
    have h := orthonormal_iff_ite.mp hρ.eigenvectorBasis.orthonormal i k
    rw [EuclideanSpace.inner_eq_star_dotProduct] at h
    rw [dotProduct_comm]; exact h
  have hp_eig : ∀ i, ρ *ᵥ p i = (lam i : ℂ) • p i := by
    intro i
    rw [hp, hρ.mulVec_eigenvectorBasis i]
    ext a; simp only [Pi.smul_apply, Complex.real_smul, smul_eq_mul, hlam]
  set g : α → β → ℂ := fun i => Ψᴴ *ᵥ p i with hg
  have hgg : ∀ i k, star (g i) ⬝ᵥ g k = if i = k then (lam i : ℂ) else 0 := by
    intro i k
    rw [hg, star_mulVec_dotProduct_mulVec, Matrix.conjTranspose_conjTranspose, ← hρΨ, hp_eig,
      dotProduct_smul, hp_orth]
    split_ifs with h
    · subst h; simp
    · simp
  set K := univ.filter (fun i => 0 < lam i) with hK
  set q : α → β → ℂ := fun i => if 0 < lam i then ((√(lam i) : ℝ) : ℂ)⁻¹ • g i else 0
    with hq
  have hq_orth : ∀ i ∈ K, ∀ k ∈ K, star (q i) ⬝ᵥ q k = if i = k then 1 else 0 := by
    intro i hi k hk
    simp only [hK, Finset.mem_filter, Finset.mem_univ, true_and] at hi hk
    simp only [hq, hi, hk, ite_true, star_smul, smul_dotProduct, dotProduct_smul, hgg]
    split_ifs with h
    · subst h
      have hs : (0:ℝ) < √(lam i) := Real.sqrt_pos.mpr hi
      have hs' : ((√(lam i) : ℝ) : ℂ) ≠ 0 := by exact_mod_cast hs.ne'
      have hl : ((lam i : ℝ) : ℂ) = ((√(lam i) : ℝ) : ℂ) * ((√(lam i) : ℝ) : ℂ) := by
        rw [← Complex.ofReal_mul, Real.mul_self_sqrt (hlam0 i)]
      rw [hl, smul_eq_mul, smul_eq_mul, star_inv₀, Complex.star_def, Complex.conj_ofReal]
      field_simp
    · simp
  have hsq : ∀ i, (((√(lam i)) : ℝ) : ℂ) = star (q i) ⬝ᵥ g i := by
    intro i
    by_cases hi : 0 < lam i
    · simp only [hq, hi, ite_true, star_smul, smul_dotProduct, hgg]
      have hs : (0:ℝ) < √(lam i) := Real.sqrt_pos.mpr hi
      have hs' : ((√(lam i) : ℝ) : ℂ) ≠ 0 := by exact_mod_cast hs.ne'
      have hl : ((lam i : ℝ) : ℂ) = ((√(lam i) : ℝ) : ℂ) * ((√(lam i) : ℝ) : ℂ) := by
        rw [← Complex.ofReal_mul, Real.mul_self_sqrt (hlam0 i)]
      rw [hl, smul_eq_mul, star_inv₀, Complex.star_def, Complex.conj_ofReal]
      field_simp
    · have h0 : lam i = 0 := le_antisymm (not_lt.mp hi) (hlam0 i)
      simp [hq, h0]
  have hΨH : ∀ v : α → ℂ, Ψᴴ *ᵥ v = ∑ j, (star (X j) ⬝ᵥ v) • star (Y j) := by
    intro v; funext b
    simp only [Matrix.mulVec, dotProduct, Matrix.conjTranspose_apply, hΨ, Finset.sum_apply,
      Pi.smul_apply, smul_eq_mul, Pi.star_apply, star_sum, star_mul', Finset.sum_mul]
    rw [Finset.sum_comm]
    exact Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun a _ => by ring
  have hexp : ∀ i, star (q i) ⬝ᵥ g i =
      ∑ j, (star (X j) ⬝ᵥ p i) * (star (q i) ⬝ᵥ star (Y j)) := by
    intro i
    rw [hg]
    simp only [hΨH, dotProduct_sum, dotProduct_smul, smul_eq_mul]
  have hzero : ∀ j, ∀ i ∉ K, star (q i) ⬝ᵥ star (Y j) = 0 := by
    intro j i hi
    simp only [hK, Finset.mem_filter, Finset.mem_univ, true_and] at hi
    simp [hq, hi]
  calc ∑ i, √(lam i) = ∑ i, ‖star (q i) ⬝ᵥ g i‖ := by
        refine Finset.sum_congr rfl fun i _ => ?_
        rw [← hsq i, Complex.norm_real, Real.norm_of_nonneg (Real.sqrt_nonneg _)]
    _ ≤ ∑ i, ∑ j, ‖star (X j) ⬝ᵥ p i‖ * ‖star (q i) ⬝ᵥ star (Y j)‖ := by
        refine Finset.sum_le_sum fun i _ => ?_
        rw [hexp i]
        refine (norm_sum_le _ _).trans (le_of_eq ?_)
        simp only [norm_mul]
    _ = ∑ j, ∑ i, ‖star (X j) ⬝ᵥ p i‖ * ‖star (q i) ⬝ᵥ star (Y j)‖ := Finset.sum_comm
    _ ≤ ∑ j, √(∑ i, ‖star (X j) ⬝ᵥ p i‖ ^ 2) * √(∑ i, ‖star (q i) ⬝ᵥ star (Y j)‖ ^ 2) :=
        Finset.sum_le_sum fun j _ => Real.sum_mul_le_sqrt_mul_sqrt _ _ _
    _ ≤ ∑ j, l2 (X j) * l2 (Y j) := by
        refine Finset.sum_le_sum fun j _ => ?_
        apply mul_le_mul (Real.sqrt_le_sqrt ?_) (Real.sqrt_le_sqrt ?_) (Real.sqrt_nonneg _)
          (l2_nonneg _)
        · have hb := bessel_on p univ (fun i _ k _ => hp_orth i k) (X j)
          refine le_of_eq_of_le (Finset.sum_congr rfl fun i _ => ?_) hb
          rw [norm_star_dotProduct_comm]
        · have hb := bessel_on q K hq_orth (star (Y j))
          rw [nsq_star] at hb
          refine le_of_eq_of_le ?_ hb
          symm
          apply Finset.sum_subset (Finset.subset_univ K)
          intro i _ hi
          rw [hzero j i hi, norm_zero]
          norm_num

/-- **Entropy of a sum of product vectors**: for a unit vector `v = ∑_j X_j ⊗ Y_j`,
`S(ρ_A) ≤ 2 log B` whenever `∑_j ‖X_j‖ ‖Y_j‖ ≤ B`. -/
lemma vonNeumann_rdm_le (ℓ : ℕ) (hℓ : ℓ ≤ N) (v : Cfg N → ℂ) (hv : nsq v = 1) {J : Type*}
    [Fintype J] (X : J → Cfg ℓ → ℂ) (Y : J → Cfg (N - ℓ) → ℂ)
    (hdec : v = ∑ j, tens ℓ hℓ (X j) (Y j)) (B : ℝ) (hB : ∑ j, l2 (X j) * l2 (Y j) ≤ B) :
    vonNeumann (rdm ℓ hℓ v) (rdm_isHermitian ℓ hℓ v) ≤ 2 * Real.log B := by
  set Ψ : Matrix (Cfg ℓ) (Cfg (N - ℓ)) ℂ :=
    Matrix.of fun a b => v ((splitCfg ℓ hℓ).symm (a, b)) with hΨdef
  have hρΨ : rdm ℓ hℓ v = Ψ * Ψᴴ := by
    ext a a'; simp [rdm, Ψ, Matrix.mul_apply, Matrix.conjTranspose_apply]
  have hΨ : ∀ a b, Ψ a b = ∑ j, X j a * Y j b := by
    intro a b
    simp only [Ψ, Matrix.of_apply, hdec, Finset.sum_apply, tens, Equiv.apply_symm_apply]
  set hρ := rdm_isHermitian ℓ hℓ v
  have hpsd : (rdm ℓ hℓ v).PosSemidef := hρΨ ▸ Matrix.posSemidef_self_mul_conjTranspose Ψ
  have h0 : ∀ i, 0 ≤ hρ.eigenvalues i := fun i => hpsd.eigenvalues_nonneg i
  have htr : ∑ i, hρ.eigenvalues i = 1 := by
    have h1 := hρ.trace_eq_sum_eigenvalues
    have h2 : (rdm ℓ hℓ v).trace = ((nsq v : ℝ) : ℂ) := by
      rw [Matrix.trace, nsq, Complex.ofReal_sum, ← (splitCfg ℓ hℓ).symm.sum_comp,
        Fintype.sum_prod_type]
      refine Finset.sum_congr rfl fun a _ => ?_
      simp only [Matrix.diag_apply, rdm]
      refine Finset.sum_congr rfl fun b _ => ?_
      rw [Complex.star_def, Complex.mul_conj', Complex.ofReal_pow]
    rw [h2, hv] at h1
    apply Complex.ofReal_injective
    rw [Complex.ofReal_sum]
    exact h1.symm
  have hR := sum_negMulLog_le_two_log_sum_sqrt _ h0 htr
  have hN := sum_sqrt_eigenvalues_le (rdm ℓ hℓ v) hρ Ψ hρΨ X Y hΨ
  have hpos : 0 < ∑ i, √(hρ.eigenvalues i) := by
    by_contra hcon
    have hle : ∑ i, √(hρ.eigenvalues i) = 0 :=
      le_antisymm (not_lt.mp hcon) (Finset.sum_nonneg fun i _ => Real.sqrt_nonneg _)
    have hz := (Finset.sum_eq_zero_iff_of_nonneg (fun i _ => Real.sqrt_nonneg _)).mp hle
    have : ∑ i, hρ.eigenvalues i = 0 := Finset.sum_eq_zero fun i hi => by
      have := hz i hi; rwa [Real.sqrt_eq_zero (h0 i)] at this
    linarith
  calc vonNeumann (rdm ℓ hℓ v) hρ = ∑ i, Real.negMulLog (hρ.eigenvalues i) := rfl
    _ ≤ 2 * Real.log (∑ i, √(hρ.eigenvalues i)) := hR
    _ ≤ 2 * Real.log B :=
        mul_le_mul_of_nonneg_left (Real.log_le_log hpos (hN.trans hB)) (by norm_num)

/-! ### Canonical orbitals -/

/-- **Canonical orbitals**: a unitary change of basis of the occupied orbitals that makes the
overlap matrices of their restrictions to `R` and to its complement diagonal. -/
lemma canonical_orbitals {M : ℕ} (Φ : Matrix (Fin N) (Fin M) ℂ) (hΦ : Φᴴ * Φ = 1)
    (R : Finset (Fin N)) :
    ∃ V : Matrix (Fin M) (Fin M) ℂ, V ∈ Matrix.unitaryGroup (Fin M) ℂ ∧ ∃ ν : Fin M → ℝ,
      (∀ k l, ∑ x ∈ R, star ((Φ * V) x k) * (Φ * V) x l = if k = l then (ν k : ℂ) else 0) ∧
      (∀ k l, ∑ x ∈ univ \ R, star ((Φ * V) x k) * (Φ * V) x l =
        if k = l then ((1 - ν k : ℝ) : ℂ) else 0) := by
  classical
  set P : Matrix (Fin N) (Fin N) ℂ := diagonal (fun x => if x ∈ R then 1 else 0) with hPdef
  have hP : Pᴴ = P := by
    rw [hPdef, Matrix.diagonal_conjTranspose]
    congr 1; funext x; by_cases hx : x ∈ R <;> simp [hx]
  set G := Φᴴ * P * Φ with hGdef
  have hG : G.IsHermitian := by
    unfold Matrix.IsHermitian
    rw [hGdef, Matrix.conjTranspose_mul, Matrix.conjTranspose_mul, hP,
      Matrix.conjTranspose_conjTranspose, Matrix.mul_assoc]
  set U := hG.eigenvectorUnitary with hU
  set V : Matrix (Fin M) (Fin M) ℂ := (U : Matrix (Fin M) (Fin M) ℂ) with hV
  have hV1 : star V * V = 1 := Unitary.coe_star_mul_self U
  have hV2 : V * star V = 1 := Unitary.coe_mul_star_self U
  have hspec : G = V * diagonal (RCLike.ofReal ∘ hG.eigenvalues) * star V := by
    have := hG.spectral_theorem
    rwa [Unitary.conjStarAlgAut_apply] at this
  have hD : star V * G * V = diagonal (RCLike.ofReal ∘ hG.eigenvalues) := by
    conv_lhs => rw [hspec]
    simp only [← Matrix.mul_assoc, hV1, Matrix.one_mul]
    rw [Matrix.mul_assoc, hV1, Matrix.mul_one]
  have hstar : (Φ * V)ᴴ = star V * Φᴴ := by
    rw [Matrix.conjTranspose_mul, Matrix.star_eq_conjTranspose]
  have hA : (Φ * V)ᴴ * P * (Φ * V) = diagonal (RCLike.ofReal ∘ hG.eigenvalues) := by
    rw [hstar, ← hD, hGdef]
    simp only [Matrix.mul_assoc]
  have hB : (Φ * V)ᴴ * (Φ * V) = 1 := by
    rw [hstar, Matrix.mul_assoc, ← Matrix.mul_assoc Φᴴ, hΦ, Matrix.one_mul, hV1]
  have eA : ∀ k l, ((Φ * V)ᴴ * P * (Φ * V)) k l =
      ∑ x ∈ R, star ((Φ * V) x k) * (Φ * V) x l := by
    intro k l
    rw [Matrix.mul_apply]
    simp only [hPdef, Matrix.mul_diagonal, Matrix.conjTranspose_apply, mul_ite, mul_one,
      mul_zero, ite_mul, zero_mul]
    rw [Finset.sum_ite_mem, Finset.univ_inter]
  have eB : ∀ k l, ((Φ * V)ᴴ * (Φ * V)) k l = ∑ x, star ((Φ * V) x k) * (Φ * V) x l := by
    intro k l
    rw [Matrix.mul_apply]
    simp only [Matrix.conjTranspose_apply]
  refine ⟨V, U.2, hG.eigenvalues, fun k l => ?_, fun k l => ?_⟩
  · rw [← eA, hA, Matrix.diagonal_apply]
    rfl
  · have h1 := eB k l
    rw [hB, Matrix.one_apply] at h1
    have h2 : ∑ x ∈ R, star ((Φ * V) x k) * (Φ * V) x l = if k = l then (hG.eigenvalues k : ℂ) else 0 := by
      rw [← eA, hA, Matrix.diagonal_apply]
      rfl
    have h3 := Finset.sum_add_sum_compl R (fun x => star ((Φ * V) x k) * (Φ * V) x l)
    rw [Finset.compl_eq_univ_sdiff] at h3
    have h4 : ∑ x ∈ univ \ R, star ((Φ * V) x k) * (Φ * V) x l =
        (if k = l then 1 else 0) - (if k = l then (hG.eigenvalues k : ℂ) else 0) := by
      rw [h1, ← h2, ← h3]; ring
    rw [h4]
    split_ifs <;> push_cast <;> ring

/-- The vector `x ↦ [x ∈ T] c u(x)`. -/
def rvec (T : Finset (Fin N)) (c : ℂ) (u : Fin N → ℂ) : Fin N → ℂ :=
  fun x => if x ∈ T then c * u x else 0

/-- Inner product of two vectors restricted to `T`. -/
lemma rvec_dot (T : Finset (Fin N)) (c c' : ℂ) (u u' : Fin N → ℂ) :
    star (rvec T c u) ⬝ᵥ rvec T c' u' = star c * c' * ∑ x ∈ T, star (u x) * u' x := by
  calc star (rvec T c u) ⬝ᵥ rvec T c' u'
      = ∑ x, if x ∈ T then star c * c' * (star (u x) * u' x) else 0 := by
        refine Finset.sum_congr rfl fun x _ => ?_
        simp only [rvec, Pi.star_apply]
        split_ifs
        · rw [star_mul']; ring
        · simp
    _ = star c * c' * ∑ x ∈ T, star (u x) * u' x := by
        rw [Finset.sum_ite_mem, Finset.univ_inter, Finset.mul_sum]

/-- Inner product of a restricted vector with a basis vector. -/
lemma rvec_dot_single (T : Finset (Fin N)) (c : ℂ) (u : Fin N → ℂ) (x : Fin N) :
    star (rvec T c u) ⬝ᵥ Pi.single x 1 = star (rvec T c u x) := by
  simp [dotProduct, Pi.single_apply]

/-- `∑_k √(ν_k(1-ν_k)) ≤ ∑_{x ∈ R, y ∉ R} |C_xy|` for canonical orbitals (`C = Φ'Φ'†`). -/
lemma sum_sqrt_nu_le {M : ℕ} (Φ' : Matrix (Fin N) (Fin M) ℂ) (R : Finset (Fin N))
    (ν : Fin M → ℝ)
    (h1 : ∀ k l, ∑ x ∈ R, star (Φ' x k) * Φ' x l = if k = l then (ν k : ℂ) else 0)
    (h2 : ∀ k l, ∑ x ∈ univ \ R, star (Φ' x k) * Φ' x l =
      if k = l then ((1 - ν k : ℝ) : ℂ) else 0) :
    ∑ k, √(ν k * (1 - ν k)) ≤ ∑ x ∈ R, ∑ y ∈ univ \ R, ‖(Φ' * Φ'ᴴ) x y‖ := by
  classical
  have hdiag : ∀ (T : Finset (Fin N)) (k : Fin M),
      ∑ x ∈ T, star (Φ' x k) * Φ' x k = ((∑ x ∈ T, ‖Φ' x k‖ ^ 2 : ℝ) : ℂ) := by
    intro T k; push_cast
    refine Finset.sum_congr rfl fun x _ => ?_
    rw [Complex.star_def, Complex.conj_mul']
  have hν : ∀ k, ν k = ∑ x ∈ R, ‖Φ' x k‖ ^ 2 := by
    intro k; have := h1 k k; rw [ite_eq_left rfl, hdiag] at this; exact_mod_cast this.symm
  have hν' : ∀ k, 1 - ν k = ∑ x ∈ univ \ R, ‖Φ' x k‖ ^ 2 := by
    intro k; have := h2 k k; rw [ite_eq_left rfl, hdiag] at this; exact_mod_cast this.symm
  have hν0 : ∀ k, 0 ≤ ν k := fun k => (hν k).symm ▸ Finset.sum_nonneg fun _ _ => by positivity
  have hν1 : ∀ k, 0 ≤ 1 - ν k := fun k =>
    (hν' k).symm ▸ Finset.sum_nonneg fun _ _ => by positivity
  set K := univ.filter (fun k => 0 < ν k ∧ 0 < 1 - ν k) with hK
  set c : Fin M → ℂ := fun k => ((√(ν k) : ℝ) : ℂ)⁻¹ with hc
  set d : Fin M → ℂ := fun k => ((√(1 - ν k) : ℝ) : ℂ)⁻¹ with hd
  set a : Fin M → Fin N → ℂ := fun k => rvec R (c k) (fun x => Φ' x k) with ha
  set b : Fin M → Fin N → ℂ := fun k => rvec (univ \ R) (d k) (fun x => Φ' x k) with hb
  have hcs : ∀ k, star (c k) = c k := fun k => by
    simp only [hc, star_inv₀, Complex.star_def, Complex.conj_ofReal]
  have hds : ∀ k, star (d k) = d k := fun k => by
    simp only [hd, star_inv₀, Complex.star_def, Complex.conj_ofReal]
  have hcc : ∀ k ∈ K, c k * c k * (ν k : ℂ) = 1 := by
    intro k hk
    simp only [hK, Finset.mem_filter, Finset.mem_univ, true_and] at hk
    have hs : ((√(ν k) : ℝ) : ℂ) ≠ 0 := by exact_mod_cast (Real.sqrt_pos.mpr hk.1).ne'
    have hl : ((ν k : ℝ) : ℂ) = ((√(ν k) : ℝ) : ℂ) * ((√(ν k) : ℝ) : ℂ) := by
      rw [← Complex.ofReal_mul, Real.mul_self_sqrt (hν0 k)]
    rw [hl, hc]; field_simp
  have hdd : ∀ k ∈ K, d k * d k * ((1 - ν k : ℝ) : ℂ) = 1 := by
    intro k hk
    simp only [hK, Finset.mem_filter, Finset.mem_univ, true_and] at hk
    have hs : ((√(1 - ν k) : ℝ) : ℂ) ≠ 0 := by exact_mod_cast (Real.sqrt_pos.mpr hk.2).ne'
    have hl : ((1 - ν k : ℝ) : ℂ) = ((√(1 - ν k) : ℝ) : ℂ) * ((√(1 - ν k) : ℝ) : ℂ) := by
      rw [← Complex.ofReal_mul, Real.mul_self_sqrt (hν1 k)]
    rw [hl, hd]; field_simp
  have ha_orth : ∀ k ∈ K, ∀ l ∈ K, star (a k) ⬝ᵥ a l = if k = l then 1 else 0 := by
    intro k hk l _
    simp only [ha, rvec_dot, h1, hcs]
    split_ifs with h
    · subst h; rw [hcc k hk]
    · simp
  have hb_orth : ∀ k ∈ K, ∀ l ∈ K, star (b k) ⬝ᵥ b l = if k = l then 1 else 0 := by
    intro k hk l _
    simp only [hb, rvec_dot, h2, hds]
    split_ifs with h
    · subst h; rw [hdd k hk]
    · simp
  have ha_b : ∀ x, ∑ k ∈ K, ‖a k x‖ ^ 2 ≤ 1 := by
    intro x
    have := bessel_on a K ha_orth (Pi.single x 1)
    rw [nsq_pi_single] at this
    refine le_of_eq_of_le (Finset.sum_congr rfl fun k _ => ?_) this
    rw [ha, rvec_dot_single, norm_star]
  have hb_b : ∀ y, ∑ k ∈ K, ‖b k y‖ ^ 2 ≤ 1 := by
    intro y
    have := bessel_on b K hb_orth (Pi.single y 1)
    rw [nsq_pi_single] at this
    refine le_of_eq_of_le (Finset.sum_congr rfl fun k _ => ?_) this
    rw [hb, rvec_dot_single, norm_star]
  have hC : ∀ x y, (Φ' * Φ'ᴴ) x y = ∑ m, Φ' x m * star (Φ' y m) := by
    intro x y; simp [Matrix.mul_apply, Matrix.conjTranspose_apply]
  have hsum : ∀ k, ∑ x ∈ R, ∑ y ∈ univ \ R, star (Φ' x k) * (Φ' * Φ'ᴴ) x y * Φ' y k =
      (ν k : ℂ) * ((1 - ν k : ℝ) : ℂ) := by
    intro k
    calc ∑ x ∈ R, ∑ y ∈ univ \ R, star (Φ' x k) * (Φ' * Φ'ᴴ) x y * Φ' y k
        = ∑ x ∈ R, ∑ y ∈ univ \ R, ∑ m, (star (Φ' x k) * Φ' x m) * (star (Φ' y m) * Φ' y k) := by
          refine Finset.sum_congr rfl fun x _ => Finset.sum_congr rfl fun y _ => ?_
          rw [hC, Finset.mul_sum, Finset.sum_mul]
          exact Finset.sum_congr rfl fun m _ => by ring
      _ = ∑ x ∈ R, ∑ m, ∑ y ∈ univ \ R, (star (Φ' x k) * Φ' x m) * (star (Φ' y m) * Φ' y k) :=
          Finset.sum_congr rfl fun x _ => Finset.sum_comm
      _ = ∑ m, ∑ x ∈ R, ∑ y ∈ univ \ R, (star (Φ' x k) * Φ' x m) * (star (Φ' y m) * Φ' y k) :=
          Finset.sum_comm
      _ = ∑ m, (∑ x ∈ R, star (Φ' x k) * Φ' x m) * (∑ y ∈ univ \ R, star (Φ' y m) * Φ' y k) :=
          Finset.sum_congr rfl fun m _ => (Finset.sum_mul_sum _ _ _ _).symm
      _ = ∑ m, (if k = m then (ν k : ℂ) else 0) * (if m = k then ((1 - ν m : ℝ) : ℂ) else 0) := by
          simp only [h1, h2]
      _ = (ν k : ℂ) * ((1 - ν k : ℝ) : ℂ) := by
          rw [Finset.sum_eq_single k]
          · simp
          · intro m _ hm; simp [Ne.symm hm, hm]
          · simp
  have hkey : ∀ k ∈ K, ∑ x ∈ R, ∑ y ∈ univ \ R, star (a k x) * (Φ' * Φ'ᴴ) x y * b k y =
      ((√(ν k * (1 - ν k)) : ℝ) : ℂ) := by
    intro k hk
    have e : ∀ x ∈ R, ∀ y ∈ univ \ R, star (a k x) * (Φ' * Φ'ᴴ) x y * b k y =
        (c k * d k) * (star (Φ' x k) * (Φ' * Φ'ᴴ) x y * Φ' y k) := by
      intro x hx y hy
      simp only [ha, hb, rvec, hx, hy, ite_true, star_mul', hcs]
      ring
    rw [Finset.sum_congr rfl fun x hx => Finset.sum_congr rfl fun y hy => e x hx y hy]
    simp_rw [← Finset.mul_sum]
    rw [hsum k]
    simp only [hK, Finset.mem_filter, Finset.mem_univ, true_and] at hk
    have hs1 : ((√(ν k) : ℝ) : ℂ) ≠ 0 := by exact_mod_cast (Real.sqrt_pos.mpr hk.1).ne'
    have hs2 : ((√(1 - ν k) : ℝ) : ℂ) ≠ 0 := by exact_mod_cast (Real.sqrt_pos.mpr hk.2).ne'
    have hl1 : ((ν k : ℝ) : ℂ) = ((√(ν k) : ℝ) : ℂ) * ((√(ν k) : ℝ) : ℂ) := by
      rw [← Complex.ofReal_mul, Real.mul_self_sqrt (hν0 k)]
    have hl2 : ((1 - ν k : ℝ) : ℂ) = ((√(1 - ν k) : ℝ) : ℂ) * ((√(1 - ν k) : ℝ) : ℂ) := by
      rw [← Complex.ofReal_mul, Real.mul_self_sqrt (hν1 k)]
    rw [Real.sqrt_mul (hν0 k), Complex.ofReal_mul, hl1, hl2, hc, hd]
    field_simp
  have hzero : ∀ k ∉ K, √(ν k * (1 - ν k)) = 0 := by
    intro k hk
    simp only [hK, Finset.mem_filter, Finset.mem_univ, true_and, not_and_or, not_lt] at hk
    rcases hk with h | h
    · rw [le_antisymm h (hν0 k), zero_mul, Real.sqrt_zero]
    · rw [le_antisymm h (hν1 k), mul_zero, Real.sqrt_zero]
  calc ∑ k, √(ν k * (1 - ν k)) = ∑ k ∈ K, √(ν k * (1 - ν k)) := by
        symm; apply Finset.sum_subset (Finset.subset_univ K)
        intro k _ hk; exact hzero k hk
    _ = ∑ k ∈ K, ‖∑ x ∈ R, ∑ y ∈ univ \ R, star (a k x) * (Φ' * Φ'ᴴ) x y * b k y‖ := by
        refine Finset.sum_congr rfl fun k hk => ?_
        rw [hkey k hk, Complex.norm_real, Real.norm_of_nonneg (Real.sqrt_nonneg _)]
    _ ≤ ∑ k ∈ K, ∑ x ∈ R, ∑ y ∈ univ \ R, ‖a k x‖ * ‖(Φ' * Φ'ᴴ) x y‖ * ‖b k y‖ := by
        refine Finset.sum_le_sum fun k _ => (norm_sum_le _ _).trans
          (Finset.sum_le_sum fun x _ => (norm_sum_le _ _).trans (le_of_eq ?_))
        refine Finset.sum_congr rfl fun y _ => ?_
        rw [norm_mul, norm_mul, norm_star]
    _ = ∑ x ∈ R, ∑ y ∈ univ \ R, ‖(Φ' * Φ'ᴴ) x y‖ * ∑ k ∈ K, ‖a k x‖ * ‖b k y‖ := by
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl fun x _ => ?_
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl fun y _ => ?_
        rw [Finset.mul_sum]
        exact Finset.sum_congr rfl fun k _ => by ring
    _ ≤ ∑ x ∈ R, ∑ y ∈ univ \ R, ‖(Φ' * Φ'ᴴ) x y‖ * 1 := by
        refine Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ =>
          mul_le_mul_of_nonneg_left ?_ (norm_nonneg _)
        have h3 : ∑ k ∈ K, ‖a k x‖ * ‖b k y‖ ≤
            ∑ k ∈ K, (‖a k x‖ ^ 2 + ‖b k y‖ ^ 2) / 2 :=
          Finset.sum_le_sum fun k _ => by nlinarith [sq_nonneg (‖a k x‖ - ‖b k y‖)]
        rw [← Finset.sum_div, Finset.sum_add_distrib] at h3
        linarith [ha_b x, hb_b y]
    _ = ∑ x ∈ R, ∑ y ∈ univ \ R, ‖(Φ' * Φ'ᴴ) x y‖ := by simp

/-- `log(√ν + √(1-ν)) ≤ √(ν(1-ν))`. -/
lemma log_sqrt_add_sqrt_le (ν : ℝ) (h0 : 0 ≤ ν) (h1 : 0 ≤ 1 - ν) :
    Real.log (√ν + √(1 - ν)) ≤ √(ν * (1 - ν)) := by
  have hsq : (√ν + √(1 - ν)) ^ 2 = 1 + 2 * (√ν * √(1 - ν)) := by
    nlinarith [Real.sq_sqrt h0, Real.sq_sqrt h1]
  have hs : 0 < √ν + √(1 - ν) := by
    by_contra hneg
    have h : √ν + √(1 - ν) = 0 :=
      le_antisymm (not_lt.mp hneg) (add_nonneg (Real.sqrt_nonneg _) (Real.sqrt_nonneg _))
    rw [h] at hsq
    nlinarith [mul_nonneg (Real.sqrt_nonneg ν) (Real.sqrt_nonneg (1 - ν))]
  have hlog := Real.log_le_sub_one_of_pos (pow_pos hs 2)
  rw [Real.log_pow, hsq] at hlog
  rw [Real.sqrt_mul h0]
  push_cast at hlog
  linarith

/-- `√ν + √(1-ν) > 0` for `ν ≥ 0`. -/
lemma sqrt_add_sqrt_pos (ν : ℝ) (h0 : 0 ≤ ν) : 0 < √ν + √(1 - ν) := by
  rcases eq_or_lt_of_le h0 with h | h
  · rw [← h]; simp
  · have := Real.sqrt_pos.mpr h
    exact add_pos_of_pos_of_nonneg this (Real.sqrt_nonneg _)

/-- `(c X) ⊗ Y = c (X ⊗ Y)`. -/
lemma tens_smul_left (ℓ : ℕ) (hℓ : ℓ ≤ N) (c : ℂ) (X : Cfg ℓ → ℂ) (Y : Cfg (N - ℓ) → ℂ) :
    tens ℓ hℓ (c • X) Y = c • tens ℓ hℓ X Y := by
  funext n; simp [tens, mul_assoc]

/-- A list sum as a sum over the positions of the list. -/
lemma list_sum_eq_sum_fin {α β : Type*} [AddCommMonoid β] (L : List α) (f : α → β) :
    (L.map f).sum = ∑ j : Fin L.length, f (L.get j) := by
  conv_lhs => rw [← List.ofFn_get L]
  rw [List.map_ofFn, List.sum_ofFn]
  rfl

end AreaLawAux

open AreaLawAux

/-- **Entanglement of Slater determinants**: the entanglement entropy of a left segment is
bounded by `2 ∑_{x ∈ A, y ∉ A} |C_{xy}|` (the paper has the sharper constant `2 ln 2`). -/
theorem vonNeumann_slater_le (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (rdm ℓ hℓ (evo (dGamma A) t *ᵥ fockVec S))
        (rdm_isHermitian ℓ hℓ _) ≤
      2 * ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ‖corrMat A S t x y‖ := by
  classical
  -- the evolved Fock state is a Slater determinant of the orbitals `U e_z`, `z ∈ S`
  obtain ⟨M, z, c, hz, hzS, hc, hv⟩ := evo_fockVec_eq A hA t S
  set Φ : Matrix (Fin N) (Fin M) ℂ := Matrix.of fun x k => evo A t x (z k) with hΦdef
  have hΦ : Φᴴ * Φ = 1 := by
    ext k l
    have h := congrFun (congrFun (evo_unitary A hA t) (z k)) (z l)
    simp only [Matrix.mul_apply, Matrix.conjTranspose_apply, Matrix.one_apply] at h ⊢
    simp only [hΦdef, Matrix.of_apply]
    rw [h]
    exact if_congr hz.eq_iff rfl rfl
  have hC : ∀ x y, corrMat A S t x y = (Φ * Φᴴ) x y := by
    intro x y
    rw [corrMat, Matrix.mul_apply, Matrix.mul_apply]
    simp only [Matrix.mul_diagonal, Matrix.conjTranspose_apply, mul_ite, mul_one, mul_zero,
      ite_mul, zero_mul, hΦdef, Matrix.of_apply]
    rw [Finset.sum_ite_mem, Finset.univ_inter]
    have hS : S = Finset.univ.image z := by
      ext x; simp [hzS x]
    rw [hS, Finset.sum_image (fun k _ l _ h => hz h)]
  -- canonical orbitals
  obtain ⟨V, hVu, ν, h1, h2⟩ := canonical_orbitals Φ hΦ (segment ℓ)
  have hV1 : star V * V = 1 := Matrix.mem_unitaryGroup_iff'.mp hVu
  have hV2 : V * star V = 1 := Matrix.mem_unitaryGroup_iff.mp hVu
  have hCC : (Φ * V) * (Φ * V)ᴴ = Φ * Φᴴ := by
    rw [Matrix.conjTranspose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc V,
      ← Matrix.star_eq_conjTranspose, hV2, Matrix.one_mul]
  have hdet : ‖V.det‖ = 1 := by
    have h := congrArg Matrix.det hV1
    rw [Matrix.det_mul, Matrix.star_eq_conjTranspose, Matrix.det_conjTranspose, Matrix.det_one,
      Complex.star_def, Complex.conj_mul'] at h
    have h' : ‖V.det‖ ^ 2 = 1 := by exact_mod_cast h
    exact (pow_eq_one_iff_of_nonneg (norm_nonneg _) two_ne_zero).mp h'
  have hdet0 : V.det ≠ 0 := by
    intro h0; rw [h0, norm_zero] at hdet; exact zero_ne_one hdet
  set c' : ℂ := c * (V.det)⁻¹ with hc'def
  have hc' : ‖c'‖ = 1 := by rw [hc'def, norm_mul, norm_inv, hc, hdet]; norm_num
  have hv' : evo (dGamma A) t *ᵥ fockVec S = c' • slaterAlt N M (fun k x => (Φ * V) x k) := by
    rw [hv, slaterAlt_mul Φ V, smul_smul, hc'def, mul_assoc, inv_mul_cancel₀ hdet0, mul_one]
    rfl
  -- product decomposition across the cut
  obtain ⟨L, hL, hLn⟩ := slater_decomp ℓ hℓ (fun k x => (Φ * V) x k)
  have hvn : nsq (evo (dGamma A) t *ᵥ fockVec S) = 1 := by
    rw [nsq_mulVec_unitary _ (evo_unitary _ (dGamma_isSymm A hA) t), fockVec, nsq_pi_single]
  have hdec : evo (dGamma A) t *ᵥ fockVec S =
      ∑ j : Fin L.length, tens ℓ hℓ (c' • (L.get j).1) (L.get j).2 := by
    rw [hv', hL, list_sum_eq_sum_fin, Finset.smul_sum]
    exact Finset.sum_congr rfl fun j _ => (tens_smul_left ℓ hℓ c' _ _).symm
  have hsq : ∀ (T : Finset (Fin N)) (k : Fin M), ∑ x ∈ T, star ((Φ * V) x k) * (Φ * V) x k =
      ((∑ x ∈ T, ‖(Φ * V) x k‖ ^ 2 : ℝ) : ℂ) := by
    intro T k
    push_cast
    exact Finset.sum_congr rfl fun x _ => by rw [Complex.star_def, Complex.conj_mul']
  have hνA : ∀ k, ν k = ∑ x ∈ segment ℓ, ‖(Φ * V) x k‖ ^ 2 := by
    intro k
    have h := h1 k k
    rw [ite_eq_left rfl, hsq] at h
    exact_mod_cast h.symm
  have hνB : ∀ k, 1 - ν k = ∑ x ∈ univ \ segment ℓ, ‖(Φ * V) x k‖ ^ 2 := by
    intro k
    have h := h2 k k
    rw [ite_eq_left rfl, hsq] at h
    exact_mod_cast h.symm
  have hlA : ∀ k, l2 (fun a => (Φ * V) (inL ℓ hℓ a) k) = √(ν k) := by
    intro k
    rw [l2, nsq, hνA k, sum_segment ℓ hℓ (fun x => ‖(Φ * V) x k‖ ^ 2)]
  have hlB : ∀ k, l2 (fun b => (Φ * V) (inR ℓ hℓ b) k) = √(1 - ν k) := by
    intro k
    rw [l2, nsq, hνB k, sum_compl_segment ℓ hℓ (fun x => ‖(Φ * V) x k‖ ^ 2)]
  have hν0 : ∀ k, 0 ≤ ν k := fun k => by
    rw [hνA k]; exact Finset.sum_nonneg fun _ _ => by positivity
  have hν1 : ∀ k, 0 ≤ 1 - ν k := fun k => by
    rw [hνB k]; exact Finset.sum_nonneg fun _ _ => by positivity
  have hB : ∑ j : Fin L.length, l2 (c' • (L.get j).1) * l2 (L.get j).2 ≤
      ∏ k, (√(ν k) + √(1 - ν k)) := by
    simp only [l2_smul, hc', one_mul]
    rw [← list_sum_eq_sum_fin L (fun p => l2 p.1 * l2 p.2)]
    refine hLn.trans (le_of_eq ?_)
    exact Finset.prod_congr rfl fun k _ => by rw [hlA k, hlB k]
  have hent := vonNeumann_rdm_le ℓ hℓ _ hvn (fun j => c' • (L.get j).1) (fun j => (L.get j).2)
    hdec _ hB
  have hlog : Real.log (∏ k, (√(ν k) + √(1 - ν k))) ≤ ∑ k, √(ν k * (1 - ν k)) := by
    rw [Real.log_prod (fun k _ => (sqrt_add_sqrt_pos (ν k) (hν0 k)).ne')]
    exact Finset.sum_le_sum fun k _ => log_sqrt_add_sqrt_le (ν k) (hν0 k) (hν1 k)
  have hnu := sum_sqrt_nu_le (Φ * V) (segment ℓ) ν h1 h2
  rw [hCC] at hnu
  calc _ ≤ 2 * Real.log (∏ k, (√(ν k) + √(1 - ν k))) := hent
    _ ≤ 2 * ∑ k, √(ν k * (1 - ν k)) := mul_le_mul_of_nonneg_left hlog (by norm_num)
    _ ≤ 2 * ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ‖(Φ * Φᴴ) x y‖ :=
        mul_le_mul_of_nonneg_left hnu (by norm_num)
    _ = 2 * ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ‖corrMat A S t x y‖ := by
        simp only [hC]

namespace AreaLawAux

/-! ### Averages over charge sectors -/

/-- Averaging is monotone. -/
lemma avgQ_mono' {f g : Cfg N → ℝ} (h : ∀ q, f q ≤ g q) : avgQ f ≤ avgQ g := by
  unfold avgQ
  exact mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun q _ => h q) (by positivity)

/-- Averaging commutes with finite sums. -/
lemma avgQ_sum' {ι : Type*} (s : Finset ι) (f : ι → Cfg N → ℝ) :
    avgQ (fun q => ∑ i ∈ s, f i q) = ∑ i ∈ s, avgQ (f i) := by
  unfold avgQ; rw [Finset.sum_comm, Finset.mul_sum]

/-- Averaging commutes with multiplication by constants. -/
lemma avgQ_const_mul' (c : ℝ) (f : Cfg N → ℝ) : avgQ (fun q => c * f q) = c * avgQ f := by
  unfold avgQ; rw [← Finset.mul_sum]; ring

/-- Cauchy–Schwarz for averages: `E[√a √b] ≤ √(E a) √(E b)`. -/
lemma avgQ_sqrt_mul_sqrt_le (a b : Cfg N → ℝ) (ha : ∀ q, 0 ≤ a q) (hb : ∀ q, 0 ≤ b q) :
    avgQ (fun q => √(a q) * √(b q)) ≤ √(avgQ a) * √(avgQ b) := by
  unfold avgQ
  have h := Real.sum_sqrt_mul_sqrt_le univ ha hb
  have hc : (0 : ℝ) ≤ 1 / 2 ^ N := by positivity
  calc (1 / 2 ^ N : ℝ) * ∑ q, √(a q) * √(b q)
      ≤ (1 / 2 ^ N : ℝ) * (√(∑ q, a q) * √(∑ q, b q)) := mul_le_mul_of_nonneg_left h hc
    _ = √((1 / 2 ^ N : ℝ) * ∑ q, a q) * √((1 / 2 ^ N : ℝ) * ∑ q, b q) := by
        rw [Real.sqrt_mul hc, Real.sqrt_mul hc, mul_mul_mul_comm, Real.mul_self_sqrt hc]

/-! ### Counting pairs across the cut -/

/-- A finite sum of values of `f` at injectively chosen points is at most `∑' f`. -/
lemma sum_le_tsum_injOn {ι : Type*} (s : Finset ι) (k : ι → ℕ) (hk : Set.InjOn k s)
    (f : ℕ → ℝ) (hf0 : ∀ n, 0 ≤ f n) (hf : Summable f) : ∑ i ∈ s, f (k i) ≤ ∑' n, f n := by
  classical
  rw [← Finset.sum_image hk]
  exact hf.sum_le_tsum _ (fun n _ => hf0 n)

/-- `∑_x f(|x - z|) ≤ 2 ∑_d f(d)`. -/
lemma sum_dist_le (f : ℕ → ℝ) (hf0 : ∀ n, 0 ≤ f n) (hf : Summable f) (z : Fin N) :
    ∑ x : Fin N, f (Nat.dist x.val z.val) ≤ 2 * ∑' n, f n := by
  rw [← Finset.sum_filter_add_sum_filter_not univ (fun x : Fin N => x.val ≤ z.val)]
  have h1 : ∑ x ∈ univ.filter (fun x : Fin N => x.val ≤ z.val), f (Nat.dist x.val z.val) ≤
      ∑' n, f n := by
    have e : ∀ x ∈ univ.filter (fun x : Fin N => x.val ≤ z.val),
        f (Nat.dist x.val z.val) = f (z.val - x.val) := by
      intro x hx
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hx
      rw [Nat.dist_eq_sub_of_le hx]
    rw [Finset.sum_congr rfl e]
    refine sum_le_tsum_injOn _ (fun x : Fin N => z.val - x.val) ?_ f hf0 hf
    intro x hx y hy hxy
    simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq] at hx hy
    exact Fin.ext (by simp only at hxy; omega)
  have h2 : ∑ x ∈ univ.filter (fun x : Fin N => ¬ x.val ≤ z.val), f (Nat.dist x.val z.val) ≤
      ∑' n, f n := by
    have e : ∀ x ∈ univ.filter (fun x : Fin N => ¬ x.val ≤ z.val),
        f (Nat.dist x.val z.val) = f (x.val - z.val) := by
      intro x hx
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hx
      rw [Nat.dist_eq_sub_of_le_right (by omega)]
    rw [Finset.sum_congr rfl e]
    refine sum_le_tsum_injOn _ (fun x : Fin N => x.val - z.val) ?_ f hf0 hf
    intro x hx y hy hxy
    simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq] at hx hy
    exact Fin.ext (by simp only at hxy; omega)
  linarith

/-- Summability of `1/(d+1)²`. -/
lemma summable_inv_succ_sq : Summable (fun d : ℕ => 1 / ((d : ℝ) + 1) ^ 2) := by
  have h := (summable_nat_add_iff 1).mpr (Real.summable_one_div_nat_pow.mpr one_lt_two)
  refine h.congr fun d => ?_
  push_cast; ring

/-- Uniform bound on the cut sum `∑_{x<ℓ≤y} ∑_z g(|x-z|) g(|y-z|)`. -/
lemma cut_sum_le (g : ℕ → ℝ) (hg0 : ∀ d, 0 ≤ g d)
    (hφ : Summable (fun d : ℕ => ((d : ℝ) + 1) ^ 2 * g d)) (ℓ : ℕ) :
    ∑ x ∈ segment (N := N) ℓ, ∑ y ∈ univ \ segment (N := N) ℓ, ∑ z : Fin N,
        g (Nat.dist x.val z.val) * g (Nat.dist y.val z.val) ≤
      8 * (∑' d, g d) * (∑' d : ℕ, ((d : ℝ) + 1) ^ 2 * g d) *
        (∑' d : ℕ, 1 / ((d : ℝ) + 1) ^ 2) := by
  set φ : ℕ → ℝ := fun d => ((d : ℝ) + 1) ^ 2 * g d with hφdef
  have hφ0 : ∀ d, 0 ≤ φ d := fun d => mul_nonneg (by positivity) (hg0 d)
  have hgφ : ∀ d, g d ≤ φ d := fun d => by
    have : (1 : ℝ) ≤ ((d : ℝ) + 1) ^ 2 := by nlinarith [(Nat.cast_nonneg d : (0 : ℝ) ≤ d)]
    simpa [hφdef] using mul_le_mul_of_nonneg_right this (hg0 d)
  have hgs : Summable g := Summable.of_nonneg_of_le hg0 hgφ hφ
  set Sg := ∑' d, g d with hSg
  set Sφ := ∑' d, φ d with hSφ
  set Z := ∑' d : ℕ, 1 / ((d : ℝ) + 1) ^ 2 with hZ
  have hSg0 : 0 ≤ Sg := tsum_nonneg hg0
  have hSφ0 : 0 ≤ Sφ := tsum_nonneg hφ0
  -- weighted decay: `g d ≤ φ d / (m+1)^2` for `m ≤ d`
  have hdecay : ∀ m d : ℕ, m ≤ d → g d ≤ φ d / ((m : ℝ) + 1) ^ 2 := by
    intro m d hmd
    rw [le_div_iff₀ (by positivity)]
    have : ((m : ℝ) + 1) ^ 2 ≤ ((d : ℝ) + 1) ^ 2 := by
      have : (m : ℝ) ≤ d := by exact_mod_cast hmd
      nlinarith
    simpa [hφdef, mul_comm] using mul_le_mul_of_nonneg_left this (hg0 d)
  -- rewrite as a sum over `z`
  have e : ∑ x ∈ segment (N := N) ℓ, ∑ y ∈ univ \ segment (N := N) ℓ, ∑ z : Fin N,
        g (Nat.dist x.val z.val) * g (Nat.dist y.val z.val) =
      ∑ z : Fin N, (∑ x ∈ segment (N := N) ℓ, g (Nat.dist x.val z.val)) *
        (∑ y ∈ univ \ segment (N := N) ℓ, g (Nat.dist y.val z.val)) := by
    calc _ = ∑ x ∈ segment (N := N) ℓ, ∑ z : Fin N, ∑ y ∈ univ \ segment (N := N) ℓ,
          g (Nat.dist x.val z.val) * g (Nat.dist y.val z.val) :=
          Finset.sum_congr rfl fun x _ => Finset.sum_comm
      _ = ∑ z : Fin N, ∑ x ∈ segment (N := N) ℓ, ∑ y ∈ univ \ segment (N := N) ℓ,
          g (Nat.dist x.val z.val) * g (Nat.dist y.val z.val) := Finset.sum_comm
      _ = _ := Finset.sum_congr rfl fun z _ => (Finset.sum_mul_sum _ _ _ _).symm
  rw [e, ← Finset.sum_filter_add_sum_filter_not univ (fun z : Fin N => z.val < ℓ)]
  -- the full sums are bounded by `2 Sg` resp. `2 Sφ`
  have hA : ∀ z : Fin N, ∑ x ∈ segment (N := N) ℓ, g (Nat.dist x.val z.val) ≤ 2 * Sg := fun z =>
    (Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _) (fun _ _ _ => hg0 _)).trans
      (sum_dist_le g hg0 hgs z)
  have hB : ∀ z : Fin N, ∑ y ∈ univ \ segment (N := N) ℓ, g (Nat.dist y.val z.val) ≤ 2 * Sg := fun z =>
    (Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _) (fun _ _ _ => hg0 _)).trans
      (sum_dist_le g hg0 hgs z)
  -- left `z`: the `y`-sum decays
  have hBl : ∀ z : Fin N, z.val < ℓ → ∑ y ∈ univ \ segment (N := N) ℓ, g (Nat.dist y.val z.val) ≤
      2 * Sφ * (1 / (((ℓ - z.val : ℕ) : ℝ) + 1) ^ 2) := by
    intro z hz
    calc ∑ y ∈ univ \ segment (N := N) ℓ, g (Nat.dist y.val z.val)
        ≤ ∑ y ∈ univ \ segment (N := N) ℓ, φ (Nat.dist y.val z.val) / (((ℓ - z.val : ℕ) : ℝ) + 1) ^ 2 := by
          refine Finset.sum_le_sum fun y hy => hdecay _ _ ?_
          simp only [segment, Finset.mem_sdiff, Finset.mem_univ, Finset.mem_filter, true_and,
            not_lt] at hy
          rw [Nat.dist_eq_sub_of_le_right (by omega)]
          omega
      _ = (∑ y ∈ univ \ segment (N := N) ℓ, φ (Nat.dist y.val z.val)) *
            (1 / (((ℓ - z.val : ℕ) : ℝ) + 1) ^ 2) := by
          rw [Finset.sum_mul]; exact Finset.sum_congr rfl fun y _ => by ring
      _ ≤ 2 * Sφ * (1 / (((ℓ - z.val : ℕ) : ℝ) + 1) ^ 2) := by
          refine mul_le_mul_of_nonneg_right ?_ (by positivity)
          exact (Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _)
            (fun _ _ _ => hφ0 _)).trans (sum_dist_le φ hφ0 hφ z)
  -- right `z`: the `x`-sum decays
  have hAr : ∀ z : Fin N, ¬ z.val < ℓ → ∑ x ∈ segment (N := N) ℓ, g (Nat.dist x.val z.val) ≤
      2 * Sφ * (1 / (((z.val - ℓ + 1 : ℕ) : ℝ) + 1) ^ 2) := by
    intro z hz
    calc ∑ x ∈ segment (N := N) ℓ, g (Nat.dist x.val z.val)
        ≤ ∑ x ∈ segment (N := N) ℓ, φ (Nat.dist x.val z.val) / (((z.val - ℓ + 1 : ℕ) : ℝ) + 1) ^ 2 := by
          refine Finset.sum_le_sum fun x hx => hdecay _ _ ?_
          simp only [segment, Finset.mem_filter, Finset.mem_univ, true_and] at hx
          rw [Nat.dist_eq_sub_of_le (by omega)]
          omega
      _ = (∑ x ∈ segment (N := N) ℓ, φ (Nat.dist x.val z.val)) *
            (1 / (((z.val - ℓ + 1 : ℕ) : ℝ) + 1) ^ 2) := by
          rw [Finset.sum_mul]; exact Finset.sum_congr rfl fun y _ => by ring
      _ ≤ 2 * Sφ * (1 / (((z.val - ℓ + 1 : ℕ) : ℝ) + 1) ^ 2) := by
          refine mul_le_mul_of_nonneg_right ?_ (by positivity)
          exact (Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _)
            (fun _ _ _ => hφ0 _)).trans (sum_dist_le φ hφ0 hφ z)
  have hZl : ∑ z ∈ univ.filter (fun z : Fin N => z.val < ℓ),
      (1 / (((ℓ - z.val : ℕ) : ℝ) + 1) ^ 2) ≤ Z := by
    refine sum_le_tsum_injOn _ (fun z : Fin N => ℓ - z.val) ?_ (fun d => 1 / ((d : ℝ) + 1) ^ 2)
      (fun d => by positivity) summable_inv_succ_sq
    intro x hx y hy hxy
    simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq] at hx hy
    exact Fin.ext (by simp only at hxy; omega)
  have hZr : ∑ z ∈ univ.filter (fun z : Fin N => ¬ z.val < ℓ),
      (1 / (((z.val - ℓ + 1 : ℕ) : ℝ) + 1) ^ 2) ≤ Z := by
    refine sum_le_tsum_injOn _ (fun z : Fin N => z.val - ℓ + 1) ?_
      (fun d => 1 / ((d : ℝ) + 1) ^ 2) (fun d => by positivity) summable_inv_succ_sq
    intro x hx y hy hxy
    simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq] at hx hy
    exact Fin.ext (by simp only at hxy; omega)
  have hl : ∑ z ∈ univ.filter (fun z : Fin N => z.val < ℓ),
      (∑ x ∈ segment (N := N) ℓ, g (Nat.dist x.val z.val)) *
        (∑ y ∈ univ \ segment (N := N) ℓ, g (Nat.dist y.val z.val)) ≤ 4 * Sg * Sφ * Z := by
    calc _ ≤ ∑ z ∈ univ.filter (fun z : Fin N => z.val < ℓ),
          (2 * Sg) * (2 * Sφ * (1 / (((ℓ - z.val : ℕ) : ℝ) + 1) ^ 2)) := by
          refine Finset.sum_le_sum fun z hz => ?_
          simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hz
          exact mul_le_mul (hA z) (hBl z hz) (Finset.sum_nonneg fun _ _ => hg0 _) (by positivity)
      _ = 4 * Sg * Sφ * ∑ z ∈ univ.filter (fun z : Fin N => z.val < ℓ),
          (1 / (((ℓ - z.val : ℕ) : ℝ) + 1) ^ 2) := by
          rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun z _ => by ring
      _ ≤ 4 * Sg * Sφ * Z := mul_le_mul_of_nonneg_left hZl (by positivity)
  have hr : ∑ z ∈ univ.filter (fun z : Fin N => ¬ z.val < ℓ),
      (∑ x ∈ segment (N := N) ℓ, g (Nat.dist x.val z.val)) *
        (∑ y ∈ univ \ segment (N := N) ℓ, g (Nat.dist y.val z.val)) ≤ 4 * Sg * Sφ * Z := by
    calc _ ≤ ∑ z ∈ univ.filter (fun z : Fin N => ¬ z.val < ℓ),
          (2 * Sφ * (1 / (((z.val - ℓ + 1 : ℕ) : ℝ) + 1) ^ 2)) * (2 * Sg) := by
          refine Finset.sum_le_sum fun z hz => ?_
          simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hz
          exact mul_le_mul (hAr z hz) (hB z) (Finset.sum_nonneg fun _ _ => hg0 _) (by positivity)
      _ = 4 * Sg * Sφ * ∑ z ∈ univ.filter (fun z : Fin N => ¬ z.val < ℓ),
          (1 / (((z.val - ℓ + 1 : ℕ) : ℝ) + 1) ^ 2) := by
          rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun z _ => by ring
      _ ≤ 4 * Sg * Sφ * Z := mul_le_mul_of_nonneg_left hZr (by positivity)
  linarith

/-- Entries of the correlation matrix: `|C_xy| ≤ ∑_z |U_xz| |U_yz|`. -/
lemma norm_corrMat_le (A : Matrix (Fin N) (Fin N) ℝ) (S : Finset (Fin N)) (t : ℝ) (x y : Fin N) :
    ‖corrMat A S t x y‖ ≤ ∑ z, ‖evo A t x z‖ * ‖evo A t y z‖ := by
  rw [corrMat, Matrix.mul_apply]
  simp only [Matrix.mul_diagonal, Matrix.conjTranspose_apply]
  refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun z _ => ?_)
  rw [norm_mul, norm_mul, norm_star]
  have : ‖(if z ∈ S then (1 : ℂ) else 0)‖ ≤ 1 := by split_ifs <;> simp
  rw [mul_right_comm]
  exact mul_le_of_le_one_right (mul_nonneg (norm_nonneg _) (norm_nonneg _)) this

/-- `|U_xz| ≤ √Q(x,z)` (entries are bounded by `1` and by the eigenfunction correlator). -/
lemma norm_evo_le_sqrt_corrQ (J h : ℝ) (q : Cfg N) (t : ℝ) (x z : Fin N) :
    ‖evo (hq J h q) t x z‖ ≤ √(corrQ (hq J h q) (hq_isHermitian J h q) x z) := by
  have h1 := norm_evo_le_corrQ (hq J h q) (hq_isHermitian J h q) t x z
  have h2 := norm_evo_le_one (hq J h q) (hq_isSymm J h q) t x z
  rw [← Real.sqrt_sq (norm_nonneg (evo (hq J h q) t x z))]
  refine Real.sqrt_le_sqrt ?_
  nlinarith [norm_nonneg (evo (hq J h q) t x z)]

end AreaLawAux

/-- **Area law for the charge-conditioned entanglement** `S_G` (paper: Theorem `thm:arealaw`):
under Hypothesis (EC), the Bernoulli average over charge sectors of the entanglement entropy of a
left segment is bounded uniformly in time, volume, segment length and occupied set. -/
theorem lgt_areaLaw (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ K : ℝ, ∀ (N : ℕ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (S : Finset (Fin N)) (t : ℝ),
      avgQ (fun q : Cfg N =>
        vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ fockVec S))
          (rdm_isHermitian ℓ hℓ _)) ≤ K := by
  set g : ℕ → ℝ := fun d => √(F d) with hg
  have hg0 : ∀ d, 0 ≤ g d := fun d => Real.sqrt_nonneg _
  refine ⟨2 * (8 * (∑' d, g d) * (∑' d : ℕ, ((d : ℝ) + 1) ^ 2 * g d) *
    (∑' d : ℕ, 1 / ((d : ℝ) + 1) ^ 2)), fun N ℓ hℓ S t => ?_⟩
  set Q : Cfg N → Fin N → Fin N → ℝ := fun q x z => corrQ (hq J h q) (hq_isHermitian J h q) x z
    with hQ
  have hQ0 : ∀ q x z, 0 ≤ Q q x z := fun q x z => corrQ_nonneg _ _ x z
  have hpt : ∀ q : Cfg N, vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ fockVec S))
      (rdm_isHermitian ℓ hℓ _) ≤
      2 * ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ∑ z, √(Q q x z) * √(Q q y z) := by
    intro q
    refine (vonNeumann_slater_le (hq J h q) (hq_isSymm J h q) S t ℓ hℓ).trans ?_
    refine mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ => ?_)
      (by norm_num)
    refine (norm_corrMat_le _ S t x y).trans (Finset.sum_le_sum fun z _ => ?_)
    exact mul_le_mul (norm_evo_le_sqrt_corrQ J h q t x z) (norm_evo_le_sqrt_corrQ J h q t y z)
      (norm_nonneg _) (Real.sqrt_nonneg _)
  calc _ ≤ avgQ (fun q : Cfg N =>
        2 * ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ∑ z, √(Q q x z) * √(Q q y z)) :=
        avgQ_mono' hpt
    _ = 2 * ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ∑ z,
          avgQ (fun q : Cfg N => √(Q q x z) * √(Q q y z)) := by
        rw [avgQ_const_mul', avgQ_sum']
        simp_rw [avgQ_sum']
    _ ≤ 2 * ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ∑ z : Fin N,
          g (Nat.dist x.val z.val) * g (Nat.dist y.val z.val) := by
        refine mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun x _ => Finset.sum_le_sum
          fun y _ => Finset.sum_le_sum fun z _ => ?_) (by norm_num)
        refine (avgQ_sqrt_mul_sqrt_le _ _ (fun q => hQ0 q x z) (fun q => hQ0 q y z)).trans ?_
        exact mul_le_mul (Real.sqrt_le_sqrt (hEC.bound N x z))
          (Real.sqrt_le_sqrt (hEC.bound N y z)) (Real.sqrt_nonneg _) (Real.sqrt_nonneg _)
    _ ≤ _ := mul_le_mul_of_nonneg_left (cut_sum_le g hg0 hEC.summable_sqrt ℓ) (by norm_num)


end DFL
