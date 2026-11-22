import DFL.Imbalance
import DFL.Dimer

/-!
# The imbalance of the dimer state (paper: Proposition `prop:dimerimbalance`(i),(ii))

For dimer charge configurations (`q_{2k} = q_{2k+1}`, `N` even), every eigenvector `φ` of `h_q`
with eigenvalue `E` satisfies the identity
`((2h)² - E²) ∑_x (-1)^x φ(x)² = J² (φ(0)² - φ(N-1)²)`
(0-indexed sites; the paper writes it with 1-indexed sites and `λ = 2h/J`). Hence the staggered
weight of eigenvectors away from the critical energies `±2h` is a boundary term, which is the
mechanism behind the vanishing of the dimer imbalance. We also identify the infinite-time
imbalance of the dimer LGT state with the dimer average of `N⁻¹ ∑_α a_α²`. (Part (iii) of the
paper's proposition, `Ī_N → 0`, uses the infinite-volume density of states of the dimer model and is
not formalized.)
-/

open Matrix Finset Filter Topology
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- **Dimer identity** (paper: Proposition `prop:dimerimbalance`(i), 0-indexed). -/
theorem dimer_stagWeight_identity (J h : ℝ) (hN : Even N) (hN0 : 0 < N) (q : Cfg N)
    (hq : IsDimer q) (α : Fin N) :
    ((2 * h) ^ 2 - (hq_isHermitian J h q).eigenvalues α ^ 2) * stagWeight J h q α =
      J ^ 2 * (((hq_isHermitian J h q).eigenvectorBasis α ⟨0, hN0⟩) ^ 2 -
        ((hq_isHermitian J h q).eigenvectorBasis α ⟨N - 1, by omega⟩) ^ 2) := by
  sorry

/-- The density-wave imbalance of the LGT started from the dimer state (unnormalized
`dimerInit`, whose squared norm is `2^{N/2}`). -/
def dimerImbalance (J h : ℝ) (N : ℕ) (t : ℝ) : ℝ :=
  (2 / N : ℝ) * (1 / 2 ^ (N / 2) : ℝ) * ∑ x : Fin N, (-1 : ℝ) ^ x.val *
    (expect (cplx (diagonal (fun n : Cfg N => occ (n x)) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)))
      (evo (H J h) t *ᵥ dimerInit N)).re

/-- **Infinite-time imbalance of the dimer state** (paper: Proposition `prop:dimerimbalance`(ii)):
the time average converges to the dimer average of `N⁻¹ ∑_α a_α²`. -/
theorem dimerImbalance_timeAvg (J h : ℝ) (hJ : J ≠ 0) (hN : Even N) :
    Tendsto (fun T : ℝ => (1 / T) * ∫ t in (0 : ℝ)..T, dimerImbalance J h N t) atTop
      (𝓝 ((1 / 2 ^ (N / 2) : ℝ) *
        ∑ q ∈ univ.filter (fun q : Cfg N => IsDimer q), imbInf J h q)) := by
  sorry

end DFL
