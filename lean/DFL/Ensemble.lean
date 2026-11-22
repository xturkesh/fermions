import DFL.Sector

/-!
# Superposition versus ensemble (paper: Proposition `prop:entropies`(d))

The initial state `|↑…↑⟩_z ⊗ |n⟩` is a coherent superposition of all `2^N` charge sectors, while
the ensemble `ρ^ens = 2^{-N} 1_links ⊗ |n⟩⟨n|` is the incoherent mixture of the definite-charge
states `|ξ⟩_σ ⊗ |n⟩`. Every gauge-invariant observable, and every product of gauge-invariant
observables at arbitrary times (multi-time correlators), has the same expectation in the two. Only
entanglement-type quantities, which are not expectations of gauge-invariant observables, can
distinguish them (the "superposition vs sampling" point of the Google experiment).
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Heisenberg-picture observable `O(t) = e^{itH} O e^{-itH}`. -/
def heis (J h : ℝ) (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (t : ℝ) :
    Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  (evo (H J h) t)ᴴ * O * evo (H J h) t

/-- The ensemble `ρ^ens = 2^{-N} 1_links ⊗ |n⟩⟨n|` (the dephased initial state). -/
def ensemble (n : Cfg N) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  (1 / 2 ^ N : ℂ) • diagonal (fun p => if p.1 = n then 1 else 0)

namespace EnsembleAux

/-- Complexification is multiplicative. -/
lemma cplx_mul {ι : Type*} [Fintype ι] (A B : Matrix ι ι ℝ) : cplx (A * B) = cplx A * cplx B := by
  ext i k
  simp [cplx, Matrix.mul_apply]

/-- Two matrices commute if their relabellings by a bijection commute. -/
lemma comm_of_submatrix {ι κ : Type*} [Fintype ι] [Fintype κ] (A B : Matrix ι ι ℝ) (e : κ ≃ ι)
    (h : A.submatrix e e * B.submatrix e e = B.submatrix e e * A.submatrix e e) :
    A * B = B * A := by
  rw [submatrix_mul_equiv, submatrix_mul_equiv] at h
  have := congrArg (fun M => M.submatrix e.symm e.symm) h
  simpa using this

/-- The LGT Hamiltonian commutes with every Gauss operator (checked in the dual frame). -/
lemma H_mul_gauss (J h : ℝ) (j : Fin N) : H J h * gauss j = gauss j * H J h := by
  apply comm_of_submatrix _ _ dualEquiv
  rw [duality, gauss_dual, Hdual]
  ext ⟨n', q'⟩ ⟨n, q⟩
  rw [mul_diagonal, diagonal_mul, blockDiagonal_apply']
  split_ifs with hq
  · subst hq; ring
  · simp

/-- The complexified Gauss operators are self-adjoint. -/
lemma cplx_gauss_conjTranspose (j : Fin N) : (cplx (gauss j))ᴴ = cplx (gauss j) := by
  rw [sec_cplx_gauss, diagonal_conjTranspose]
  refine congrArg diagonal (funext fun p => ?_)
  simp

/-- The propagator commutes with every Gauss operator. -/
lemma evo_comm_gauss (J h : ℝ) (j : Fin N) (t : ℝ) :
    Commute (evo (H J h) t) (cplx (gauss j)) := by
  have hc : Commute (cplx (H J h)) (cplx (gauss j)) := by
    show cplx (H J h) * cplx (gauss j) = cplx (gauss j) * cplx (H J h)
    rw [← cplx_mul, ← cplx_mul, H_mul_gauss]
  exact (hc.smul_left _).exp_left

/-- The adjoint propagator commutes with every Gauss operator. -/
lemma evo_conjTranspose_comm_gauss (J h : ℝ) (j : Fin N) (t : ℝ) :
    Commute (evo (H J h) t)ᴴ (cplx (gauss j)) := by
  have h2 := congrArg conjTranspose (evo_comm_gauss J h j t).eq
  rw [conjTranspose_mul, conjTranspose_mul, cplx_gauss_conjTranspose] at h2
  exact h2.symm

/-- The expectation in a basis vector is a diagonal entry. -/
lemma expect_single {ι : Type*} [Fintype ι] [DecidableEq ι] (M : Matrix ι ι ℂ) (i : ι) :
    expect M (Pi.single i 1) = M i i := by
  simp [expect, dotProduct, mulVec, Pi.single_apply]

end EnsembleAux

open EnsembleAux

/-- Products of gauge-invariant observables are gauge invariant. -/
theorem gaugeInv_mul {O P : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ} (hO : GaugeInv O)
    (hP : GaugeInv P) : GaugeInv (O * P) := by
  intro j
  exact Commute.mul_left (hO j) (hP j)

/-- Heisenberg evolution preserves gauge invariance. -/
theorem gaugeInv_heis (J h : ℝ) {O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ}
    (hO : GaugeInv O) (t : ℝ) : GaugeInv (heis J h O t) := by
  intro j
  exact ((evo_conjTranspose_comm_gauss J h j t).mul_left (hO j)).mul_left
    (evo_comm_gauss J h j t)

/-- **Superposition = ensemble** for one gauge-invariant observable. -/
theorem expect_initZ_eq_trace_ensemble (n : Cfg N)
    (X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hX : GaugeInv X) :
    expect X (initZ (Pi.single n 1)) = trace (ensemble n * X) := by
  rw [← sec_expect_submatrix X dualEquiv, gaugeInv_submatrix X hX, sec_expect_blockDiagonal]
  have hvec : ∀ q : Cfg N, (fun n' => (initZ (Pi.single n 1) ∘ dualEquiv) (n', q)) =
      (((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ • Pi.single n 1 := by
    intro q; ext n'; simp [initZ]
  simp only [hvec, sec_expect_smul, sec_norm_sq_prefactor, expect_single]
  rw [← Finset.mul_sum, ensemble, Matrix.smul_mul, trace_smul, smul_eq_mul]
  congr 1
  have hbij : Function.Bijective (linkOf n) :=
    Function.bijective_iff_has_inverse.mpr ⟨chargeOf n, chargeOf_linkOf n, linkOf_chargeOf n⟩
  have hL : ∑ q, sectorBlock X q n n = ∑ x, X (n, x) (n, x) :=
    hbij.sum_comp (fun x => X (n, x) (n, x))
  have hR : ((diagonal fun p : Cfg N × Cfg N => if p.1 = n then (1 : ℂ) else 0) * X).trace =
      ∑ x, X (n, x) (n, x) := by
    rw [trace, Fintype.sum_prod_type, Fintype.sum_eq_single n]
    · simp [diagonal_mul]
    · intro n' hn'
      simp [diagonal_mul, hn']
  rw [hL, hR]

/-- **Superposition = ensemble for all multi-time correlators of gauge-invariant observables**
(paper: Proposition `prop:entropies`(d), eq. `eq:main-ens`). -/
theorem correlator_eq_ensemble (J h : ℝ) (n : Cfg N)
    (Os : List (Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ × ℝ)) (hOs : ∀ p ∈ Os, GaugeInv p.1) :
    expect ((Os.map (fun p => heis J h p.1 p.2)).prod) (initZ (Pi.single n 1)) =
      trace (ensemble n * (Os.map (fun p => heis J h p.1 p.2)).prod) := by
  apply expect_initZ_eq_trace_ensemble
  induction Os with
  | nil =>
    intro j
    simp
  | cons p Os ih =>
    rw [List.map_cons, List.prod_cons]
    exact gaugeInv_mul (gaugeInv_heis J h (hOs p List.mem_cons_self) p.2)
      (ih fun p' hp' => hOs p' (List.mem_cons_of_mem _ hp'))

end DFL
