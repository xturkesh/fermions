import DFL.AreaLawSharp

/-!
# Annealed one-body bounds for cuts and intervals (paper: Lemma `lem:annealed`(ii),(iii))

Under (EC), with `Σ₁(F) = ∑_{d ≥ 1} d F(d)`, `Σ'₁(F) = ∑_{d ≥ 1} d F(d)^{1/2}` and `Θ(F)` (`thetaF`):

* (ii) for the cut after the site `m`: `E ∑_{x ≤ m < z} Q_q(x, z) ≤ Σ₁(F)` (`annealed_cut_Q`) and
  `E (∑_{x ≤ m < z} Q_q(x, z))² ≤ Σ'₁(F)²` (`annealed_cut_Q_sq`); the first bound of (ii) is
  `annealed_cut`;
* (iii) for an interval `B = [a, b]`, the same bounds for the sums over `x ∈ B`, `y ∉ B`, with the
  constants `2Θ(F)`, `2Σ₁(F)` and `4Σ'₁(F)²` (`annealed_ivl_K`, `annealed_ivl_Q`, `annealed_ivl_Q_sq`).
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The interval `[a, b]` of sites. -/
def ivl (a b : ℕ) : Finset (Fin N) := univ.filter (fun x => a ≤ x.val ∧ x.val ≤ b)

/-- `E ∑_{x ≤ m < z} Q_q(x, z) ≤ Σ₁(F)` (paper: Lemma `lem:annealed`(ii)). -/
theorem annealed_cut_Q (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (m : ℕ) :
    avgQ (fun q : Cfg N => ∑ x ∈ leftOf m, ∑ z ∈ univ \ leftOf m,
      corrQ (hq J h q) (hq_isHermitian J h q) x z) ≤ ∑' d : ℕ, (d : ℝ) * F d := by
  sorry

/-- `E (∑_{x ≤ m < z} Q_q(x, z))² ≤ Σ'₁(F)²` (paper: Lemma `lem:annealed`(ii)). -/
theorem annealed_cut_Q_sq (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (m : ℕ) :
    avgQ (fun q : Cfg N => (∑ x ∈ leftOf m, ∑ z ∈ univ \ leftOf m,
      corrQ (hq J h q) (hq_isHermitian J h q) x z) ^ 2) ≤ (∑' d : ℕ, (d : ℝ) * √(F d)) ^ 2 := by
  sorry

/-- `E sup_t ∑_{x ∈ B, y ∉ B} |K_q(t)_{xy}| ≤ 2Θ(F)` for an interval `B` (paper: Lemma `lem:annealed`(iii)). -/
theorem annealed_ivl_K (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (a b : ℕ) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ,
      ∑ x ∈ ivl a b, ∑ y ∈ univ \ ivl a b, ‖corrMat (hq J h q) S t x y‖) ≤ 2 * thetaF F := by
  sorry

/-- `E ∑_{x ∈ B, z ∉ B} Q_q(x, z) ≤ 2Σ₁(F)` (paper: Lemma `lem:annealed`(iii)). -/
theorem annealed_ivl_Q (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (a b : ℕ) :
    avgQ (fun q : Cfg N => ∑ x ∈ ivl a b, ∑ z ∈ univ \ ivl a b,
      corrQ (hq J h q) (hq_isHermitian J h q) x z) ≤ 2 * ∑' d : ℕ, (d : ℝ) * F d := by
  sorry

/-- `E (∑_{x ∈ B, z ∉ B} Q_q(x, z))² ≤ 4Σ'₁(F)²` (paper: Lemma `lem:annealed`(iii)). -/
theorem annealed_ivl_Q_sq (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (a b : ℕ) :
    avgQ (fun q : Cfg N => (∑ x ∈ ivl a b, ∑ z ∈ univ \ ivl a b,
      corrQ (hq J h q) (hq_isHermitian J h q) x z) ^ 2) ≤
      4 * (∑' d : ℕ, (d : ℝ) * √(F d)) ^ 2 := by
  sorry

end DFL
