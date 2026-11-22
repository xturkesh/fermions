import DFL.Basic
import DFL.Evolution

/-!
# One-body input: the eigenfunction correlator and Hypothesis (EC)

The localization input of the paper (Section `sec:onebody`) is a statement about the
one-dimensional Anderson model with Bernoulli potential `h_q = -J Δ_N + 2h diag(q)`: the averaged
eigenfunction correlator decays, uniformly in the volume (Hypothesis `hyp:EC`). In the paper its
status is discussed in detail (Ge–Zhao on `ℤ`; Macera–Sodin and Macera for finite intervals). Here
it is an explicit hypothesis: every result that uses it takes `EC J h F` as an argument, so the
Lean statements are exactly as conditional as the paper's.

Main results:
* `DFL.hq_isHermitian`;
* `DFL.corrQ`: the eigenfunction correlator `Q(x,y) = ∑_α |φ_α(x)| |φ_α(y)|`;
* `DFL.norm_evo_le_corrQ`: `|(e^{-itA})_{xy}| ≤ Q(x,y)` for every `t`;
* `DFL.transProb_le_corrQ`: the time-averaged transition probability is `≤ Q(x,y)`;
* `DFL.EC`: Hypothesis (EC).
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

theorem pathAdj_isSymm (N : ℕ) : (pathAdj N).IsSymm := by
  refine Matrix.IsSymm.ext fun x y => ?_
  simp only [pathAdj, or_comm]

theorem hq_isSymm (J h : ℝ) (q : Cfg N) : (hq J h q).IsSymm := by
  exact ((pathAdj_isSymm N).smul _).add ((Matrix.isSymm_diagonal _).smul _)

theorem hq_isHermitian (J h : ℝ) (q : Cfg N) : (hq J h q).IsHermitian := by
  unfold Matrix.IsHermitian
  rw [Matrix.conjTranspose_eq_transpose_of_trivial]
  exact hq_isSymm J h q

/-- The eigenfunction correlator `Q(x, y) = ∑_α |φ_α(x)| |φ_α(y)|` of a real symmetric matrix, for
the orthonormal eigenbasis provided by Mathlib. For a Jacobi matrix (simple spectrum) this equals
`sup {|f(A)_{xy}| : |f| ≤ 1}`, the eigenfunction correlator of the literature. -/
def corrQ (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x y : Fin N) : ℝ :=
  ∑ α, |hA.eigenvectorBasis α x| * |hA.eigenvectorBasis α y|

/-- The time-averaged transition probability `R(x, y) = ∑_α φ_α(x)² φ_α(y)²`. -/
def transProb (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x y : Fin N) : ℝ :=
  ∑ α, (hA.eigenvectorBasis α x) ^ 2 * (hA.eigenvectorBasis α y) ^ 2

theorem corrQ_nonneg (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x y : Fin N) :
    0 ≤ corrQ A hA x y := by
  exact Finset.sum_nonneg fun _ _ => mul_nonneg (abs_nonneg _) (abs_nonneg _)

namespace OneBodyAux

/-- The rows of the eigenvector matrix are unit vectors: `∑_α φ_α(x)² = 1`. -/
lemma sum_eigvec_sq (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x : Fin N) :
    ∑ α, (hA.eigenvectorBasis α x) ^ 2 = 1 := by
  have h := congrFun (congrFun (Unitary.coe_mul_star_self hA.eigenvectorUnitary) x) x
  simp only [Matrix.mul_apply, Matrix.one_apply_eq, Unitary.coe_star, Matrix.star_apply,
    Matrix.IsHermitian.eigenvectorUnitary_apply, star_trivial] at h
  simpa [sq] using h

/-- Each eigenvector entry is bounded by one in absolute value. -/
lemma abs_eigvec_le_one (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (α x : Fin N) :
    |hA.eigenvectorBasis α x| ≤ 1 := by
  have h1 : (hA.eigenvectorBasis α x) ^ 2 ≤ 1 := by
    rw [← sum_eigvec_sq A hA x]
    exact Finset.single_le_sum (f := fun β => (hA.eigenvectorBasis β x) ^ 2)
      (fun _ _ => sq_nonneg _) (Finset.mem_univ α)
  rw [abs_le_one_iff_mul_self_le_one]
  nlinarith

/-- Spectral form of the propagator: `(e^{-itA})_{xy} = ∑_α e^{-itE_α} φ_α(x) φ_α(y)`. -/
lemma evo_apply_eq_sum (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (t : ℝ)
    (x y : Fin N) :
    evo A t x y = ∑ α, Complex.exp (-(t : ℂ) * Complex.I * (hA.eigenvalues α : ℂ)) *
      ((hA.eigenvectorBasis α x : ℂ) * (hA.eigenvectorBasis α y : ℂ)) := by
  set U : Matrix (Fin N) (Fin N) ℝ := (hA.eigenvectorUnitary : Matrix (Fin N) (Fin N) ℝ)
    with hU
  have hspec : A = U * diagonal hA.eigenvalues * star U := by
    have h3 := hA.spectral_theorem
    rw [Unitary.conjStarAlgAut_apply] at h3
    simpa using h3
  have hUU : U * star U = 1 := Unitary.coe_mul_star_self hA.eigenvectorUnitary
  have hcplx : ∀ M : Matrix (Fin N) (Fin N) ℝ, cplx M = (Complex.ofRealHom).mapMatrix M :=
    fun M => rfl
  have hinv : (cplx U)⁻¹ = cplx (star U) := by
    apply Matrix.inv_eq_right_inv
    rw [hcplx, hcplx, ← map_mul, hUU, map_one]
  have hunit : IsUnit (cplx U) := by
    refine (Matrix.isUnit_iff_isUnit_det _).mpr ?_
    have : cplx U * cplx (star U) = 1 := by rw [hcplx, hcplx, ← map_mul, hUU, map_one]
    exact Matrix.isUnit_det_of_right_inverse this
  set c : ℂ := -(t : ℂ) * Complex.I with hc
  have hdiag : c • cplx (diagonal hA.eigenvalues) =
      diagonal (fun α => c * (hA.eigenvalues α : ℂ)) := by
    ext i j
    by_cases hij : i = j
    · subst hij; simp [cplx]
    · simp [cplx, diagonal_apply_ne _ hij]
  have hconj : c • cplx A = cplx U * diagonal (fun α => c * (hA.eigenvalues α : ℂ)) *
      (cplx U)⁻¹ := by
    have hA' : cplx A = cplx U * cplx (diagonal hA.eigenvalues) * cplx (star U) := by
      simp only [hcplx, ← map_mul]
      rw [← hspec]
    rw [hinv, ← hdiag, hA', ← smul_mul_assoc, ← mul_smul_comm]
  have hexp : evo A t = cplx U * diagonal (fun α => Complex.exp (c * (hA.eigenvalues α : ℂ))) *
      cplx (star U) := by
    unfold evo
    rw [← hc, hconj, Matrix.exp_conj _ _ hunit, Matrix.exp_diagonal, hinv, Pi.exp_def]
    simp [Complex.exp_eq_exp_ℂ]
  rw [hexp, Matrix.mul_apply]
  refine Finset.sum_congr rfl fun α _ => ?_
  simp [cplx, Matrix.mul_diagonal, hU, Matrix.star_apply]
  ring

end OneBodyAux

/-- The propagator is dominated by the eigenfunction correlator, uniformly in time. -/
theorem norm_evo_le_corrQ (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (t : ℝ)
    (x y : Fin N) : ‖evo A t x y‖ ≤ corrQ A hA x y := by
  rw [OneBodyAux.evo_apply_eq_sum A hA t x y, corrQ]
  refine (norm_sum_le _ _).trans (le_of_eq (Finset.sum_congr rfl fun α _ => ?_))
  have hE : ‖Complex.exp (-(t : ℂ) * Complex.I * (hA.eigenvalues α : ℂ))‖ = 1 := by
    have : -(t : ℂ) * Complex.I * (hA.eigenvalues α : ℂ) =
        ((-(t * hA.eigenvalues α) : ℝ) : ℂ) * Complex.I := by push_cast; ring
    rw [this, Complex.norm_exp_ofReal_mul_I]
  rw [norm_mul, norm_mul, hE, one_mul, Complex.norm_real, Complex.norm_real, Real.norm_eq_abs,
    Real.norm_eq_abs]

/-- Matrix entries of a unitary propagator have modulus at most one. -/
theorem norm_evo_le_one (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (t : ℝ) (x y : Fin N) :
    ‖evo A t x y‖ ≤ 1 := by
  have h := congrFun (congrFun (evo_unitary' A hA t) x) x
  simp only [Matrix.mul_apply, Matrix.one_apply_eq, Matrix.conjTranspose_apply] at h
  have h2 : ∑ j, ‖evo A t x j‖ ^ 2 = 1 := by
    have : ∀ j, evo A t x j * star (evo A t x j) = ((‖evo A t x j‖ ^ 2 : ℝ) : ℂ) := by
      intro j
      rw [Complex.star_def, Complex.mul_conj, Complex.normSq_eq_norm_sq]
    simp only [this] at h
    exact_mod_cast h
  have h3 : ‖evo A t x y‖ ^ 2 ≤ 1 := h2 ▸ Finset.single_le_sum
    (f := fun j => ‖evo A t x j‖ ^ 2) (fun _ _ => sq_nonneg _) (Finset.mem_univ y)
  nlinarith [norm_nonneg (evo A t x y)]

/-- The time-averaged transition probability is dominated by the eigenfunction correlator. -/
theorem transProb_le_corrQ (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x y : Fin N) :
    transProb A hA x y ≤ corrQ A hA x y := by
  refine Finset.sum_le_sum fun α _ => ?_
  have hx := OneBodyAux.abs_eigvec_le_one A hA α x
  have hy := OneBodyAux.abs_eigvec_le_one A hA α y
  rw [← sq_abs (hA.eigenvectorBasis α x), ← sq_abs (hA.eigenvectorBasis α y)]
  have ha := abs_nonneg (hA.eigenvectorBasis α x)
  have hb := abs_nonneg (hA.eigenvectorBasis α y)
  have h1 : |hA.eigenvectorBasis α x| ^ 2 ≤ |hA.eigenvectorBasis α x| := by nlinarith
  have h2 : |hA.eigenvectorBasis α y| ^ 2 ≤ |hA.eigenvectorBasis α y| := by nlinarith
  exact mul_le_mul h1 h2 (sq_nonneg _) ha

/-- The uniform (Bernoulli(1/2)) average over charge configurations. -/
def avgQ (f : Cfg N → ℝ) : ℝ := (1 / 2 ^ N : ℝ) * ∑ q : Cfg N, f q

/-- **Hypothesis (EC)** (paper: `hyp:EC`): the averaged eigenfunction correlator of the Bernoulli
chain `h_q` on `{0, …, N-1}` (Dirichlet ends) is bounded by a non-increasing summable profile `F`,
uniformly in the volume `N`. Expected (and proven on `ℤ` by Ge–Zhao): `F(d) = C e^{-ηd}`. -/
structure EC (J h : ℝ) (F : ℕ → ℝ) : Prop where
  nonneg : ∀ d, 0 ≤ F d
  antitone : Antitone F
  bound : ∀ (N : ℕ) (x y : Fin N),
    avgQ (fun q => corrQ (hq J h q) (hq_isHermitian J h q) x y) ≤ F (Nat.dist x.val y.val)
  summable : Summable (fun d : ℕ => ((d : ℝ) + 1) ^ 2 * F d)
  summable_sqrt : Summable (fun d : ℕ => ((d : ℝ) + 1) ^ 2 * Real.sqrt (F d))

end DFL
