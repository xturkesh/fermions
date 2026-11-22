import DFL.AreaLaw
import DFL.VNEntropy

/-!
# Reduced states of Slater determinants (paper: Lemmas `lem:quasifree`, `lem:slater`)

For the evolved Fock state `ψ = Γ(e^{-itA}) |S⟩` (a Slater determinant with one-body projection
`K = corrMat A S t`) and the left segment `{0, …, ℓ-1}`, let `C_A = P_A K P_A` (`corrA`), with
eigenvalues `ν_1, …, ν_ℓ ∈ [0, 1]`. The canonical orbitals give a Schmidt decomposition of `ψ` across
the cut, with Schmidt coefficients `∏_{j ∈ T} ν_j ∏_{j ∉ T} (1 - ν_j)`, whence

* `vonNeumann_rdm_slater`: `S(ρ_A) = ∑_i h(ν_i)`, `h` the binary entropy (paper: `lem:slater`(a));
* `numDist_slater`: the law of `N_A` is the Poisson-binomial law of the `ν_i`;
* `genfun_slater`: its generating function is `det(1 + (z-1) C_A)` (paper: `lem:quasifree`(b));
* `sum_sqrt_eig_le`: `∑_i √(ν_i (1 - ν_i)) ≤ ∑_{x ∈ A, y ∉ A} |K_{xy}|` (paper: `lem:slater`(b)).
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- `C_A = P_A K P_A`: the correlation matrix restricted to the left segment `{0, …, ℓ-1}`. -/
def corrA (A : Matrix (Fin N) (Fin N) ℝ) (S : Finset (Fin N)) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    Matrix (Fin ℓ) (Fin ℓ) ℂ :=
  (corrMat A S t).submatrix (AreaLawAux.inL ℓ hℓ) (AreaLawAux.inL ℓ hℓ)

/-- The Poisson-binomial law of the parameters `ν`. -/
def poissonBinomial {ι : Type*} [Fintype ι] [DecidableEq ι] (ν : ι → ℝ) (k : ℕ) : ℝ :=
  ∑ T ∈ (univ : Finset (Finset ι)).filter (fun T => T.card = k),
    (∏ i ∈ T, ν i) * ∏ i ∈ univ \ T, (1 - ν i)

/-- `C_A` is Hermitian. -/
theorem corrA_isHermitian (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) : (corrA A S t ℓ hℓ).IsHermitian := by
  sorry

/-- The eigenvalues of `C_A` lie in `[0, 1]`. -/
theorem eig_corrA_mem (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (i : Fin ℓ) :
    0 ≤ (corrA_isHermitian A hA S t ℓ hℓ).eigenvalues i ∧
      (corrA_isHermitian A hA S t ℓ hℓ).eigenvalues i ≤ 1 := by
  sorry

/-- **Entanglement entropy of a Slater determinant** (paper: Lemma `lem:slater`(a)):
`S(ρ_A) = ∑_i h(ν_i)`. -/
theorem vonNeumann_rdm_slater (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (rdm ℓ hℓ (evo (dGamma A) t *ᵥ fockVec S)) (rdm_isHermitian ℓ hℓ _) =
      ∑ i, Real.binEntropy ((corrA_isHermitian A hA S t ℓ hℓ).eigenvalues i) := by
  sorry

/-- **The law of `N_A` is Poisson-binomial** (paper: Lemmas `lem:slater`(a), `lem:quasifree`(b)). -/
theorem numDist_slater (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (k : ℕ) :
    numDist (segment ℓ) (evo (dGamma A) t *ᵥ fockVec S) k =
      poissonBinomial (corrA_isHermitian A hA S t ℓ hℓ).eigenvalues k := by
  sorry

/-- **Generating function** (paper: Lemma `lem:quasifree`(b)):
`∑_k P(N_A = k) z^k = det(1 + (z - 1) C_A)`. -/
theorem genfun_slater (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (z : ℂ) :
    ∑ k ∈ range (N + 1), (numDist (segment ℓ) (evo (dGamma A) t *ᵥ fockVec S) k : ℂ) * z ^ k =
      (1 + (z - 1) • corrA A S t ℓ hℓ).det := by
  sorry

/-- **Trace-norm bound** (paper: Lemma `lem:slater`(b)):
`∑_i √(ν_i (1 - ν_i)) = ‖P_B K P_A‖₁ ≤ ∑_{x ∈ A, y ∉ A} |K_{xy}|`. -/
theorem sum_sqrt_eig_le (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    ∑ i, √((corrA_isHermitian A hA S t ℓ hℓ).eigenvalues i *
        (1 - (corrA_isHermitian A hA S t ℓ hℓ).eigenvalues i)) ≤
      ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ‖corrMat A S t x y‖ := by
  sorry

end DFL
