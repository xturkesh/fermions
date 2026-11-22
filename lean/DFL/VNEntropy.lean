import DFL.CutLocal

/-!
# Von Neumann entropy: Schur bound, concavity, mixing bound, two-block pinching
(paper: Lemma `lem:pinching`, and the entropy inequalities used in Proposition `prop:entropies`)

`vonNeumann ρ hρ = ∑_i η(λ_i)`, with `η(x) = -x log x` (`Real.negMulLog`) and `λ_i` the eigenvalues
of the Hermitian matrix `ρ` (defined in `DFL.AreaLaw`).

* `vonNeumann_le_diag`: `S(ρ) ≤ ∑_a η(ρ_aa)` for `ρ ≥ 0`. The diagonal of `ρ` is a doubly
  stochastic average of its eigenvalues, `ρ_aa = ∑_i |U_ai|² λ_i`, and `η` is concave.
* `vonNeumann_diagonal`: the entropy of a real diagonal matrix is the Shannon entropy of its diagonal.
* `vonNeumann_le_log_card`: `S(ρ) ≤ log dim` for density matrices.
* `vonNeumann_mul_conjTranspose_comm`: `S(A Aᴴ) = S(Aᴴ A)` (equal nonzero spectra).
* `vonNeumann_concave`: concavity, `∑ w_i S(σ_i) ≤ S(∑ w_i σ_i)`.
* `vonNeumann_mix_le`: the mixing bound `S(∑ w_i σ_i) ≤ ∑ w_i S(σ_i) + H(w)` for density matrices.
* `vonNeumann_le_pinching`, `vonNeumann_pinching_le`: the two-block pinching lemma
  `S(σ) ≤ S(PσP + P'σP') ≤ S(σ) + ln 2`, `P' = 1 - P`.
-/

open Matrix Finset
open scoped ComplexOrder

noncomputable section

namespace DFL

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- **Schur bound**: the von Neumann entropy of a positive semidefinite matrix is at most the
Shannon entropy of its diagonal. -/
theorem vonNeumann_le_diag (ρ : Matrix ι ι ℂ) (hρ : ρ.PosSemidef) :
    vonNeumann ρ hρ.isHermitian ≤ ∑ a, Real.negMulLog (ρ a a).re := by
  sorry

/-- The entropy of a real diagonal matrix is the Shannon entropy of its diagonal. -/
theorem vonNeumann_diagonal (d : ι → ℝ) (hd : (diagonal (fun a => (d a : ℂ))).IsHermitian) :
    vonNeumann (diagonal (fun a => (d a : ℂ))) hd = ∑ a, Real.negMulLog (d a) := by
  sorry

/-- The entropy of a density matrix is at most `log` of the dimension. -/
theorem vonNeumann_le_log_card (ρ : Matrix ι ι ℂ) (hρ : ρ.PosSemidef) (htr : ρ.trace = 1) :
    vonNeumann ρ hρ.isHermitian ≤ Real.log (Fintype.card ι) := by
  sorry

/-- `A Aᴴ` and `Aᴴ A` have the same nonzero eigenvalues, hence the same entropy. -/
theorem vonNeumann_mul_conjTranspose_comm {κ : Type*} [Fintype κ] [DecidableEq κ]
    (A : Matrix ι κ ℂ) :
    vonNeumann (A * Aᴴ) (isHermitian_mul_conjTranspose_self A) =
      vonNeumann (Aᴴ * A) (isHermitian_conjTranspose_mul_self A) := by
  sorry

/-- **Concavity** of the von Neumann entropy (finite convex combinations of positive semidefinite
matrices). -/
theorem vonNeumann_concave {κ : Type*} [Fintype κ] (w : κ → ℝ) (hw0 : ∀ i, 0 ≤ w i)
    (hw1 : ∑ i, w i = 1) (σ : κ → Matrix ι ι ℂ) (hσ : ∀ i, (σ i).PosSemidef)
    (hmix : (∑ i, (w i : ℂ) • σ i).IsHermitian) :
    ∑ i, w i * vonNeumann (σ i) (hσ i).isHermitian ≤ vonNeumann (∑ i, (w i : ℂ) • σ i) hmix := by
  sorry

/-- **Mixing bound**: `S(∑ w_i σ_i) ≤ ∑ w_i S(σ_i) + H(w)` for density matrices `σ_i` and weights
`w_i ≥ 0`. -/
theorem vonNeumann_mix_le {κ : Type*} [Fintype κ] (w : κ → ℝ) (hw0 : ∀ i, 0 ≤ w i)
    (σ : κ → Matrix ι ι ℂ) (hσ : ∀ i, (σ i).PosSemidef) (htr : ∀ i, (σ i).trace = 1)
    (hmix : (∑ i, (w i : ℂ) • σ i).IsHermitian) :
    vonNeumann (∑ i, (w i : ℂ) • σ i) hmix ≤
      ∑ i, w i * vonNeumann (σ i) (hσ i).isHermitian + ∑ i, Real.negMulLog (w i) := by
  sorry

/-- The two-block pinching `P σ P + (1 - P) σ (1 - P)`. -/
def pinching (P σ : Matrix ι ι ℂ) : Matrix ι ι ℂ := P * σ * P + (1 - P) * σ * (1 - P)

/-- The pinching of a positive semidefinite matrix by a Hermitian `P` is positive semidefinite. -/
theorem pinching_posSemidef (P σ : Matrix ι ι ℂ) (hP : P.IsHermitian) (hσ : σ.PosSemidef) :
    (pinching P σ).PosSemidef := by
  sorry

/-- **Two-block pinching, lower bound** (paper: Lemma `lem:pinching`): for an orthogonal projection
`P`, pinching does not decrease the entropy. -/
theorem vonNeumann_le_pinching (P σ : Matrix ι ι ℂ) (hP : P.IsHermitian) (hPP : P * P = P)
    (hσ : σ.PosSemidef) (h' : (pinching P σ).IsHermitian) :
    vonNeumann σ hσ.isHermitian ≤ vonNeumann (pinching P σ) h' := by
  sorry

/-- **Two-block pinching, upper bound** (paper: Lemma `lem:pinching`): for a density matrix `σ` and
an orthogonal projection `P`, pinching increases the entropy by at most `ln 2`. -/
theorem vonNeumann_pinching_le (P σ : Matrix ι ι ℂ) (hP : P.IsHermitian) (hPP : P * P = P)
    (hσ : σ.PosSemidef) (htr : σ.trace = 1) (h' : (pinching P σ).IsHermitian) :
    vonNeumann (pinching P σ) h' ≤ vonNeumann σ hσ.isHermitian + Real.log 2 := by
  sorry

end DFL
