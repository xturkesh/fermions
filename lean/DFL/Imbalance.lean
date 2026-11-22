import DFL.OneBody

/-!
# The density-wave imbalance (paper: Theorem `thm:imbalance`)

One-body part of the imbalance theorem: the chiral anticommutation identity, simplicity of the
spectrum of the Jacobi matrix `h_q`, the infinite-time (time-averaged) imbalance
`Ī_q = N⁻¹ ∑_α a_α²` with `a_α = ⟨φ_α, Γ φ_α⟩`, and lower bounds on its Bernoulli average:
a finite-`N` bound for every `h ≠ 0`, and a bound uniform in `N` under Hypothesis (EC).
-/

open Matrix Finset Filter Topology

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Staggered weight `a_α = ⟨φ_α, Γ φ_α⟩ = ∑_x (-1)^x φ_α(x)²` of the eigenvector `φ_α` of `h_q`. -/
def stagWeight (J h : ℝ) (q : Cfg N) (α : Fin N) : ℝ :=
  ∑ x, (-1 : ℝ) ^ x.val * ((hq_isHermitian J h q).eigenvectorBasis α x) ^ 2

/-- The one-body imbalance in sector `q`: `N⁻¹ Re Tr(Γ U(t) Γ U(t)†)`, `U(t) = e^{-ith_q}`. -/
def imbOne (J h : ℝ) (q : Cfg N) (t : ℝ) : ℝ :=
  (1 / N : ℝ) *
    (trace (cplx (stag N) * evo (hq J h q) t * cplx (stag N) * (evo (hq J h q) t)ᴴ)).re

/-- The infinite-time imbalance in sector `q`: `N⁻¹ ∑_α a_α²`. -/
def imbInf (J h : ℝ) (q : Cfg N) : ℝ := (1 / N : ℝ) * ∑ α, stagWeight J h q α ^ 2

/-- `Γ` anticommutes with the path adjacency matrix (the path is bipartite). -/
theorem stag_pathAdj_anticomm (N : ℕ) : stag N * pathAdj N + pathAdj N * stag N = 0 := by
  sorry

/-- Chiral identity `Γ h_q + h_q Γ = 4h V_q Γ`, with `V_q = diag(q)`. -/
theorem stag_hq_anticomm (J h : ℝ) (q : Cfg N) :
    stag N * hq J h q + hq J h q * stag N =
      (4 * h) • (diagonal (fun x => sgn (q x)) * stag N) := by
  sorry

/-- The eigenvalue identity `E_α a_α = 2h ⟨φ_α, V_q Γ φ_α⟩`. -/
theorem eigenvalue_mul_stagWeight (J h : ℝ) (q : Cfg N) (α : Fin N) :
    (hq_isHermitian J h q).eigenvalues α * stagWeight J h q α =
      2 * h * ∑ x, sgn (q x) * (-1 : ℝ) ^ x.val *
        ((hq_isHermitian J h q).eigenvectorBasis α x) ^ 2 := by
  sorry

/-- Jacobi matrices with non-zero off-diagonal entries have simple spectrum. -/
theorem hq_eigenvalues_injective (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) :
    Function.Injective (hq_isHermitian J h q).eigenvalues := by
  sorry

/-- The time average of the one-body imbalance converges to `Ī_q = N⁻¹ ∑_α a_α²`. -/
theorem imbOne_timeAvg (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) :
    Tendsto (fun T : ℝ => (1 / T) * ∫ t in (0 : ℝ)..T, imbOne J h q t) atTop
      (𝓝 (imbInf J h q)) := by
  sorry

/-- **Finite-volume lower bound** for every `h ≠ 0`:
`2^{-N} ∑_q Ī_q ≥ h² / (N (|J| + |h|)²)`. -/
theorem avg_imbInf_ge (J h : ℝ) (hJ : J ≠ 0) (hh : h ≠ 0) (hN : 1 ≤ N) :
    h ^ 2 / (N * (|J| + |h|) ^ 2) ≤ avgQ (fun q : Cfg N => imbInf J h q) := by
  sorry

/-- **Uniform lower bound** under Hypothesis (EC): the averaged infinite-time imbalance is bounded
below by a positive constant independent of the volume. -/
theorem avg_imbInf_uniform (J h : ℝ) (hJ : J ≠ 0) (hh : h ≠ 0) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ c : ℝ, 0 < c ∧ ∀ (N : ℕ), 1 ≤ N → c ≤ avgQ (fun q : Cfg N => imbInf J h q) := by
  sorry

end DFL
