import DFL.AnnealedInterval
import DFL.TransportSharp

/-!
# Number entropy of bulk intervals in finite volume (paper: proof of Theorem `thm:infvol`(e))

For an interval `B = [a, b]` of sites and a Fock state with occupied set `S`, `s_B = |S ∩ B|`:

* `⟨Ψ(t), (N_B - s_B)² Ψ(t)⟩ ≤ 2Θ(F) + 4Σ'₁(F)²` (`lgt_bulk_transport_sq`): the law of total variance,
  `Var_q(N_B) = ∑_{x ∈ B, y ∉ B} |K_{xy}|² ≤ ∑_{x ∈ B, y ∉ B} |K_{xy}|` and `|⟨N_B⟩_q - s_B| ≤ ∑_{x ∈ B, z ∉ B} Q(x, z)`,
  with Lemma `lem:annealed`(iii);
* `S_N(t) ≤ ½ ln(2πe(2Θ(F) + 4Σ'₁(F)² + 1/12))` (`lgt_numberEntropy_bulk`), uniformly in `N`, `t`, `B`, `S`.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- **Second moment of the charge in a bulk interval**: `⟨Ψ(t), (N_B - s_B)² Ψ(t)⟩ ≤ 2Θ(F) + 4Σ'₁(F)²`. -/
theorem lgt_bulk_transport_sq (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (a b : ℕ)
    (t : ℝ) :
    (expect (matter ((numOp (ivl a b) - ((S ∩ ivl a b).card : ℝ) • 1) *
        (numOp (ivl a b) - ((S ∩ ivl a b).card : ℝ) • 1))) (lgtState J h (fockVec S) t)).re ≤
      2 * thetaF F + 4 * sigma1' F ^ 2 := by
  sorry

/-- **Bounded number entropy of a bulk interval**:
`S_N(t) ≤ ½ ln(2πe(2Θ(F) + 4Σ'₁(F)² + 1/12))`. -/
theorem lgt_numberEntropy_bulk (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (a b : ℕ)
    (t : ℝ) :
    shannon N (lgtNumDist J h (ivl a b) (fockVec S) t) ≤
      (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (2 * thetaF F + 4 * sigma1' F ^ 2 + 1 / 12)) := by
  sorry

end DFL
