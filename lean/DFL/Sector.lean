import DFL.Duality
import DFL.Evolution

/-!
# Exact sector averaging (paper: Corollary `cor:sector`)

For the initial state `|↑…↑⟩_z ⊗ ψ` (all links polarized along `z`), the dual state is the uniform
superposition `2^{-N/2} ∑_q ψ ⊗ |q⟩` over all `2^N` charge sectors. Hence the expectation of every
gauge-invariant observable at every time is an exact average over i.i.d. Bernoulli(1/2) charges of
free-fermion expectations.
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The LGT initial state `|↑…↑⟩_z ⊗ ψ` written in the `σˣ` basis of the links:
`Ψ₀(n, x) = 2^{-N/2} ψ(n)` for every link configuration `x`. -/
def initZ (ψ : Cfg N → ℂ) : Cfg N × Cfg N → ℂ :=
  fun p => (((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ * ψ p.1

/-- The LGT state at time `t` evolved from `|↑…↑⟩_z ⊗ ψ`. -/
def lgtState (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) : Cfg N × Cfg N → ℂ :=
  evo (H J h) t *ᵥ initZ ψ

/-- Gauge invariance: `O` commutes with every Gauss operator. -/
def GaugeInv (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) : Prop :=
  ∀ j : Fin N, O * cplx (gauss j) = cplx (gauss j) * O

/-- The block of `O` in the charge sector `q`, in the dual frame. -/
def sectorBlock (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (q : Cfg N) :
    Matrix (Cfg N) (Cfg N) ℂ :=
  fun n' n => O (dualEquiv (n', q)) (dualEquiv (n, q))

/-! ### Auxiliary lemmas -/

/-- `ZMod 2` has exactly the two elements `0` and `1`. -/
theorem sec_zmod2_cases (a : ZMod 2) : a = 0 ∨ a = 1 := by
  revert a; decide

/-- Distinct `ZMod 2` labels have distinct signs. -/
theorem sec_sgn_ne {a b : ZMod 2} (hab : a ≠ b) : sgn a ≠ sgn b := by
  have h10 : (1 : ZMod 2) ≠ 0 := by decide
  rcases sec_zmod2_cases a with rfl | rfl <;> rcases sec_zmod2_cases b with rfl | rfl <;>
    simp_all [sgn] <;> norm_num

/-- Complexification of a diagonal matrix. -/
theorem sec_cplx_diagonal {κ : Type*} [DecidableEq κ] (d : κ → ℝ) :
    cplx (diagonal d) = diagonal (fun k => (d k : ℂ)) := by
  unfold cplx
  exact diagonal_map (by simp)

/-- Entrywise form of commuting with a diagonal matrix. -/
theorem sec_comm_diagonal_apply {κ : Type*} [Fintype κ] [DecidableEq κ] (O : Matrix κ κ ℂ)
    (g : κ → ℂ) (hc : O * diagonal g = diagonal g * O) (x y : κ) :
    O x y * g y = g x * O x y := by
  have := congrFun (congrFun hc x) y
  rwa [mul_diagonal, diagonal_mul] at this

/-- The (complexified) Gauss operator is diagonal. -/
theorem sec_cplx_gauss (j : Fin N) :
    cplx (gauss j) = diagonal (fun p => ((gauss j p p : ℝ) : ℂ)) := by
  unfold gauss
  rw [sec_cplx_diagonal]
  simp [diagonal_apply_eq]

/-- The diagonal entries of a Gauss operator in the dual frame are the charge signs. -/
theorem sec_gauss_dual_apply (j : Fin N) (a : Cfg N × Cfg N) :
    gauss j (dualEquiv a) (dualEquiv a) = sgn (a.2 j) := by
  have h := congrFun (congrFun (gauss_dual j) a) a
  rwa [submatrix_apply, diagonal_apply_eq] at h

/-- A diagonal matter observable `D ⊗ 1` is diagonal in `(n, x)`. -/
theorem sec_diagonal_kronecker_one (d : Cfg N → ℝ) :
    diagonal d ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) = diagonal (fun p => d p.1) := by
  rw [← diagonal_one, diagonal_kronecker_diagonal]
  simp

/-- Gauge-invariant observables are block diagonal in the charges (dual frame). -/
theorem gaugeInv_submatrix (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) :
    O.submatrix dualEquiv dualEquiv = blockDiagonal (sectorBlock O) := by
  ext ⟨n', q'⟩ ⟨n, q⟩
  rw [submatrix_apply, blockDiagonal_apply']
  split_ifs with hq
  · subst hq; rfl
  · obtain ⟨j, hj⟩ := Function.ne_iff.mp hq
    have key := sec_comm_diagonal_apply O _ (by rw [← sec_cplx_gauss]; exact hO j)
      (dualEquiv (n', q')) (dualEquiv (n, q))
    simp only [sec_gauss_dual_apply] at key
    have hne : ((sgn (q j) : ℝ) : ℂ) ≠ ((sgn (q' j) : ℝ) : ℂ) := by
      exact_mod_cast (sec_sgn_ne hj).symm
    have hzero : O (dualEquiv (n', q')) (dualEquiv (n, q)) *
        (((sgn (q j) : ℝ) : ℂ) - ((sgn (q' j) : ℝ) : ℂ)) = 0 := by
      rw [mul_sub, key]; ring
    rcases mul_eq_zero.mp hzero with h0 | h0
    · exact h0
    · exact absurd (sub_eq_zero.mp h0) hne

/-- Diagonal matter observables `D ⊗ 1` (densities, number operators and their spectral
projections) are gauge invariant. -/
theorem gaugeInv_matter (d : Cfg N → ℝ) :
    GaugeInv (cplx (diagonal d ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))) := by
  intro j
  rw [sec_diagonal_kronecker_one, sec_cplx_gauss, sec_cplx_diagonal, diagonal_mul_diagonal,
    diagonal_mul_diagonal]
  congr 1
  ext p
  ring

/-- The sector blocks of a diagonal matter observable `D ⊗ 1` are all equal to `D`. -/
theorem sectorBlock_matter (d : Cfg N → ℝ) (q : Cfg N) :
    sectorBlock (cplx (diagonal d ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))) q = cplx (diagonal d) := by
  ext n' n
  rw [sec_diagonal_kronecker_one, sec_cplx_diagonal, sec_cplx_diagonal]
  unfold sectorBlock
  by_cases hn : n' = n
  · subst hn; simp only [diagonal_apply_eq, dualEquiv_apply]
  · have hne : dualEquiv (n', q) ≠ dualEquiv (n, q) := fun h' =>
      hn (congrArg Prod.fst (dualEquiv.injective h'))
    rw [diagonal_apply_ne _ hne, diagonal_apply_ne _ hn]

/-- Expectations are invariant under relabelling the basis by a bijection. -/
theorem sec_expect_submatrix {κ ι : Type*} [Fintype κ] [Fintype ι] (M : Matrix ι ι ℂ)
    (e : κ ≃ ι) (v : ι → ℂ) : expect (M.submatrix e e) (v ∘ e) = expect M v := by
  unfold expect
  rw [submatrix_mulVec_equiv, Function.comp_assoc, Equiv.self_comp_symm, Function.comp_id]
  exact comp_equiv_dotProduct_comp_equiv (star v) (M *ᵥ v) e

/-- Matrix–vector product with a block-diagonal matrix, block by block. -/
theorem sec_blockDiagonal_mulVec {m o : Type*} [Fintype m] [Fintype o] [DecidableEq o]
    (M : o → Matrix m m ℂ) (v : m × o → ℂ) (i : m) (k : o) :
    (blockDiagonal M *ᵥ v) (i, k) = (M k *ᵥ fun j => v (j, k)) i := by
  simp only [mulVec, dotProduct, blockDiagonal_apply', Fintype.sum_prod_type]
  refine Finset.sum_congr rfl fun j _ => ?_
  simp [ite_mul]

/-- Expectation of a block-diagonal observable is the sum of the blockwise expectations. -/
theorem sec_expect_blockDiagonal {m o : Type*} [Fintype m] [Fintype o] [DecidableEq o]
    (M : o → Matrix m m ℂ) (v : m × o → ℂ) :
    expect (blockDiagonal M) v = ∑ k, expect (M k) (fun j => v (j, k)) := by
  unfold expect
  simp only [dotProduct, Fintype.sum_prod_type, sec_blockDiagonal_mulVec, Pi.star_apply]
  exact Finset.sum_comm

/-- Multiplying the state by a scalar `z` multiplies expectations by `|z|²`. -/
theorem sec_expect_smul {ι : Type*} [Fintype ι] (M : Matrix ι ι ℂ) (z : ℂ) (v : ι → ℂ) :
    expect M (z • v) = (star z * z) * expect M v := by
  unfold expect
  rw [mulVec_smul, star_smul, smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul]
  ring

/-- **Dual-frame expectation formula**: for every initial state `Ψ₀` and every gauge-invariant
`O`, the expectation in the evolved LGT state is the sum over charge sectors `q` of free-fermion
expectations, with initial sector states `n ↦ Ψ₀ (dualEquiv (n, q))`. -/
theorem sec_expect_evo_dual (J h : ℝ) (Ψ₀ : Cfg N × Cfg N → ℂ)
    (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (H J h) t *ᵥ Ψ₀) =
      ∑ q : Cfg N, expect (sectorBlock O q)
        (evo (Hsec J h q) t *ᵥ fun n => Ψ₀ (dualEquiv (n, q))) := by
  have hU : blockDiagonal (fun q : Cfg N => evo (Hsec J h q) t) =
      (evo (H (N := N) J h) t).submatrix dualEquiv dualEquiv := by
    rw [evo_submatrix, duality, ← evo_blockDiagonal]
    rfl
  have hstate : (evo (H J h) t *ᵥ Ψ₀) ∘ dualEquiv =
      blockDiagonal (fun q : Cfg N => evo (Hsec J h q) t) *ᵥ (Ψ₀ ∘ dualEquiv) := by
    rw [hU, submatrix_mulVec_equiv, Function.comp_assoc, Equiv.self_comp_symm, Function.comp_id]
  rw [← sec_expect_submatrix O dualEquiv, gaugeInv_submatrix O hO, hstate,
    sec_expect_blockDiagonal]
  refine Finset.sum_congr rfl fun q _ => ?_
  congr 1
  ext n
  rw [sec_blockDiagonal_mulVec]
  rfl

/-- `|2^{-N/2}|² = 2^{-N}`. -/
theorem sec_norm_sq_prefactor :
    star ((((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹) * (((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ = (1 / 2 ^ N : ℂ) := by
  rw [star_inv₀, Complex.star_def, Complex.conj_ofReal, ← mul_inv, ← Complex.ofReal_mul,
    ← mul_pow, Real.mul_self_sqrt (by norm_num : (0 : ℝ) ≤ 2)]
  push_cast
  ring

/-- **Exact sector averaging** (paper: Corollary `cor:sector`): for every gauge-invariant `O` and
every time `t`, `⟨Ψ(t)|O|Ψ(t)⟩ = 2^{-N} ∑_q ⟨ψ_q(t)|O_q|ψ_q(t)⟩` with
`ψ_q(t) = e^{-it (dΓ(h_q) - h∑q)} ψ`. -/
theorem sector_average (J h : ℝ) (ψ : Cfg N → ℂ) (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ)
    (hO : GaugeInv O) (t : ℝ) :
    expect O (lgtState J h ψ t) =
      (1 / 2 ^ N : ℂ) * ∑ q : Cfg N, expect (sectorBlock O q) (evo (Hsec J h q) t *ᵥ ψ) := by
  unfold lgtState
  rw [sec_expect_evo_dual J h _ O hO t, Finset.mul_sum]
  refine Finset.sum_congr rfl fun q _ => ?_
  have hvec : (fun n => initZ ψ (dualEquiv (n, q))) =
      (((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ • ψ := by
    ext n
    simp [initZ]
  rw [hvec, mulVec_smul, sec_expect_smul, sec_norm_sq_prefactor]

/-- Sector averaging for diagonal matter observables, with the constant phase removed:
`⟨Ψ(t)| D ⊗ 1 |Ψ(t)⟩ = 2^{-N} ∑_q ⟨e^{-itdΓ(h_q)}ψ| D |e^{-itdΓ(h_q)}ψ⟩`. -/
theorem sector_average_matter (J h : ℝ) (ψ : Cfg N → ℂ) (d : Cfg N → ℝ) (t : ℝ) :
    expect (cplx (diagonal d ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))) (lgtState J h ψ t) =
      (1 / 2 ^ N : ℂ) * ∑ q : Cfg N,
        expect (cplx (diagonal d)) (evo (dGamma (hq J h q)) t *ᵥ ψ) := by
  rw [sector_average J h ψ _ (gaugeInv_matter d) t]
  refine congrArg (fun s => (1 / 2 ^ N : ℂ) * s) (Finset.sum_congr rfl fun q _ => ?_)
  rw [sectorBlock_matter]
  unfold Hsec
  exact expect_evo_sub_smul_one _ _ _ _ _

end DFL
