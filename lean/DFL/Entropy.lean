import Mathlib

/-!
# Entropy inequalities

Shannon entropy of distributions on `{0, …, M}`, the support bound, Massey's maximum-entropy bound
for integer-valued random variables (paper: proof of Theorem `thm:numberentropy`), the variance of a
mixture, and the binary-entropy bound `h(ν) ≤ 2 √(ν(1-ν))` used for the area law.
-/

open Finset Real

noncomputable section

namespace DFL

/-- Shannon entropy (in nats) of a distribution `p` on `{0, …, M}`. -/
def shannon (M : ℕ) (p : ℕ → ℝ) : ℝ := ∑ k ∈ range (M + 1), negMulLog (p k)

/-- Mean of a distribution on `{0, …, M}`. -/
def pmean (M : ℕ) (p : ℕ → ℝ) : ℝ := ∑ k ∈ range (M + 1), p k * k

/-- Variance of a distribution on `{0, …, M}`. -/
def pvar (M : ℕ) (p : ℕ → ℝ) : ℝ := ∑ k ∈ range (M + 1), p k * ((k : ℝ) - pmean M p) ^ 2

/-- `p` is a probability distribution on `{0, …, M}`. -/
structure IsProb (M : ℕ) (p : ℕ → ℝ) : Prop where
  nonneg : ∀ k, 0 ≤ p k
  sum_one : ∑ k ∈ range (M + 1), p k = 1

variable {M : ℕ} {p : ℕ → ℝ}

theorem shannon_nonneg (hp : IsProb M p) : 0 ≤ shannon M p := by
  sorry

/-- The entropy is at most the logarithm of the size of the support. -/
theorem shannon_le_log_card (hp : IsProb M p) (s : Finset ℕ) (hs : s.Nonempty)
    (hsupp : ∀ k ∈ range (M + 1), k ∉ s → p k = 0) :
    shannon M p ≤ Real.log s.card := by
  sorry

/-- **Massey's bound**: the Shannon entropy of an integer-valued random variable with variance `σ²`
is at most `½ ln(2πe(σ² + 1/12))`. -/
theorem massey (hp : IsProb M p) :
    shannon M p ≤ (1 / 2) * Real.log (2 * π * exp 1 * (pvar M p + 1 / 12)) := by
  sorry

/-- The Massey bound is monotone in the variance. -/
theorem shannon_le_of_pvar_le (hp : IsProb M p) {V : ℝ} (hV : pvar M p ≤ V) :
    shannon M p ≤ (1 / 2) * Real.log (2 * π * exp 1 * (V + 1 / 12)) := by
  sorry

/-- Variance of a mixture: `Var(∑ w_i p_i) = ∑ w_i Var(p_i) + ∑ w_i (mean_i - mean)²`. -/
theorem pvar_mixture {ι : Type*} [Fintype ι] (w : ι → ℝ) (hw0 : ∀ i, 0 ≤ w i)
    (hw1 : ∑ i, w i = 1) (P : ι → ℕ → ℝ) (hP : ∀ i, IsProb M (P i)) :
    pvar M (fun k => ∑ i, w i * P i k) =
      ∑ i, w i * pvar M (P i) +
        ∑ i, w i * (pmean M (P i) - pmean M (fun k => ∑ i, w i * P i k)) ^ 2 := by
  sorry

/-- The variance is at most the mean square deviation from any point. -/
theorem pvar_le_sum_sq (hp : IsProb M p) (c : ℝ) :
    pvar M p ≤ ∑ k ∈ range (M + 1), p k * ((k : ℝ) - c) ^ 2 := by
  sorry

/-- Binary entropy bound `h(ν) ≤ 2 √(ν(1-ν))` (used for the area law). Proof idea: Rényi
monotonicity `h(ν) ≤ 2 ln(√ν + √(1-ν)) = ln(1 + 2√(ν(1-ν)))` (Jensen for `log`), then
`ln(1+x) ≤ x`. (The paper uses the sharper constant `2 ln 2`.) -/
theorem binEntropy_le_two_sqrt {ν : ℝ} (h0 : 0 ≤ ν) (h1 : ν ≤ 1) :
    binEntropy ν ≤ 2 * Real.sqrt (ν * (1 - ν)) := by
  sorry

end DFL
