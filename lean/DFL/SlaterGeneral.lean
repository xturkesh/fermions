import DFL.SlaterExact

/-!
# The law of `N_R` for Slater determinants and arbitrary regions (paper: Lemmas `lem:slater`(a), `lem:quasifree`(b))

For the evolved Fock state `ψ = Γ(e^{-itA}) |S⟩` (a Slater determinant with one-body projection
`K = corrMat A S t`) and an **arbitrary** set of sites `R` (e.g. a bulk interval), let `C_R = P_R K P_R`
(`corrR`), with eigenvalues `ν_i ∈ [0, 1]`. Then the law of `N_R` is the Poisson-binomial law of the `ν_i`
(`numDist_slater_gen`), with generating function `det(1 + (z - 1) C_R)` (`genfun_slater_gen`). For left
segments this is `numDist_slater`, `genfun_slater` of `DFL.SlaterExact`.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- `C_R = P_R K P_R`: the correlation matrix restricted to the sites of `R`. -/
def corrR (A : Matrix (Fin N) (Fin N) ℝ) (S : Finset (Fin N)) (t : ℝ) (R : Finset (Fin N)) :
    Matrix R R ℂ :=
  (corrMat A S t).submatrix (fun x : R => (x : Fin N)) (fun x : R => (x : Fin N))

/-- `C_R` is Hermitian. -/
theorem corrR_isHermitian (A : Matrix (Fin N) (Fin N) ℝ) (S : Finset (Fin N)) (t : ℝ)
    (R : Finset (Fin N)) : (corrR A S t R).IsHermitian := by
  sorry

/-- The eigenvalues of `C_R` lie in `[0, 1]`. -/
theorem eig_corrR_mem (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N)) (t : ℝ)
    (R : Finset (Fin N)) (i : R) :
    0 ≤ (corrR_isHermitian A S t R).eigenvalues i ∧ (corrR_isHermitian A S t R).eigenvalues i ≤ 1 := by
  sorry

/-- **Generating function of `N_R`** (paper: Lemma `lem:quasifree`(b)), any region `R`:
`∑_k P(N_R = k) z^k = det(1 + (z - 1) C_R)`. -/
theorem genfun_slater_gen (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N)) (t : ℝ)
    (R : Finset (Fin N)) (z : ℂ) :
    ∑ k ∈ range (N + 1), (numDist R (evo (dGamma A) t *ᵥ fockVec S) k : ℂ) * z ^ k =
      (1 + (z - 1) • corrR A S t R).det := by
  sorry

/-- **The law of `N_R` is Poisson-binomial** (paper: Lemmas `lem:slater`(a), `lem:quasifree`(b)), any
region `R`. -/
theorem numDist_slater_gen (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N)) (t : ℝ)
    (R : Finset (Fin N)) (k : ℕ) :
    numDist R (evo (dGamma A) t *ᵥ fockVec S) k =
      poissonBinomial (corrR_isHermitian A S t R).eigenvalues k := by
  sorry

end DFL
