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

/-- **Spectrum of a quadratic Hamiltonian**: `char dΓ(A) = ∏_{T} (X - ∑_{α ∈ T} λ_α)`. -/
theorem charpoly_dGamma (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) :
    (dGamma A).charpoly = ∏ T : Finset (Fin N), (X - C (∑ α ∈ T, hA.eigenvalues α)) := by
  sorry

/-- **The spectrum of `H`** (paper: Theorem `thm:duality`(b)): with multiplicities, the union over the
charge sectors `q` of `{∑_{α ∈ T} E_α(q) - h ∑_j q_j : T ⊆ {1, …, N}}`. -/
theorem charpoly_H (J h : ℝ) :
    (H (N := N) J h).charpoly =
      ∏ q : Cfg N, ∏ T : Finset (Fin N),
        (X - C (∑ α ∈ T, (hq_isHermitian J h q).eigenvalues α - h * ∑ j, sgn (q j))) := by
  sorry

end DFL
