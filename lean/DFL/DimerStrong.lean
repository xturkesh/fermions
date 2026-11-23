import DFL.DimerInfVol
import DFL.DimerState

/-!
# Strong coupling: the dressed fermion stays tight (paper: Remark `rem:dimerstrong`)

For `λ = 2|h|/J > 1`, `λ ≠ √2`, De Bièvre and Germinet prove dynamical localization for the random dimer
model on `ℤ`: almost surely, for every `p > 0`, `sup_t ∑_y |y|^p |⟨δ_y, e^{-itH}δ_0⟩|² < ∞`. This is a
literature input; here it enters as the hypothesis `hDL` for one `p > 0` (for the operator
`-Δ + v^{(r)}` of `DimerInfVol.lean`, which is the operator of De Bièvre–Germinet up to the diagonal gauge
transformation `diag((-1)^x)` and a translation).

* `dimerTailLGT`: the weight `⟨Ψ^x_N(t), ∑_{|y-x|>R} n_y Ψ^x_N(t)⟩` of the single-fermion dimer state outside
  the window of radius `R` (normalized as `dimerDLGT`); `dimerTailLGT_eq`: it is the dimer-sector
  average of `∑_{|y-x|>R} |⟨δ_y, e^{-ith_q}δ_x⟩|²` (eq. `eq:dimerdensity`).
* `dimer_TL_tail`: its thermodynamic limit is `dimerTailZ R t = E ∑_{|y|>R} |⟨δ_y, e^{-iJtH}δ_0⟩|²`
  (Lemma `lem:dimerTL`, as in the proof of Theorem `thm:dimer`(b)).
* `dimerTailZ_tendsto_zero`: given dynamical localization, `sup_t dimerTailZ R t → 0` as `R → ∞`
  (dominated convergence; the quantities are at most `1`).
* `dimer_tight_of_dynLoc`: the statement of the remark,
  `lim_{R→∞} sup_t lim_{N→∞} ⟨Ψ^{x_N}_N(t), ∑_{|y-x_N|>R} n_y Ψ^{x_N}_N(t)⟩ = 0`.
-/

open MeasureTheory Filter Topology Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

/-- The weight of the single-fermion dimer state `Ψ^x_N(t)` outside the window of radius `R` around
`x`, `⟨Ψ^x_N(t), ∑_{|y-x|>R} n_y Ψ^x_N(t)⟩` (normalized as `dimerDLGT`). -/
def dimerTailLGT (J h : ℝ) (N R : ℕ) (t : ℝ) (x : Fin N) : ℝ :=
  (1 / 2 ^ (N / 2) : ℝ) *
    (expect (cplx (diagonal (fun n : Cfg N => ∑ y : Fin N,
        (if (R : ℤ) < |(y.val : ℤ) - x.val| then (1 : ℝ) else 0) * occ (n y)) ⊗ₖ
      (1 : Matrix (Cfg N) (Cfg N) ℝ))) (evo (H J h) t *ᵥ dimerState (fockCfg {x}))).re

/-- The infinite-volume counterpart `E ∑_{|y|>R} |⟨δ_y, e^{-iJtH}δ_0⟩|²` for the random dimer potential
`v^{(r)}`. -/
def dimerTailZ (J h : ℝ) (r : ℤ) (R : ℕ) (t : ℝ) : ℝ :=
  ∫ η, ∑' y : ℤ, (if (R : ℤ) < |y| then ‖evoZw (dimerPotZ J h r η) (J * t) y 0‖ ^ 2 else 0) ∂muZ

/-- The tail weight of the dimer state is the dimer-sector average of the one-particle tail weight
(eq. `eq:dimerdensity`). -/
theorem dimerTailLGT_eq (J h : ℝ) {N : ℕ} (hN : Even N) (R : ℕ) (t : ℝ) (x : Fin N) :
    dimerTailLGT J h N R t x = dimerAvg (fun q : Cfg N =>
      ∑ y : Fin N, (if (R : ℤ) < |(y.val : ℤ) - x.val| then ‖evo (hq J h q) t y x‖ ^ 2 else 0)) := by
  sorry

/-- **Thermodynamic limit of the tail weight** (paper: Remark `rem:dimerstrong`, by Lemma
`lem:dimerTL` as in the proof of Theorem `thm:dimer`(b)). -/
theorem dimer_TL_tail (J h : ℝ) (hJ : 0 < J) (r : ℤ) (Ns : ℕ → ℕ) (xs : ∀ n, Fin (Ns n))
    (hN : ∀ n, Even (Ns n)) (hpar : ∀ n, ((xs n).val : ℤ) % 2 = r)
    (hx1 : Tendsto (fun n => (xs n).val) atTop atTop)
    (hx2 : Tendsto (fun n => Ns n - (xs n).val) atTop atTop) (R : ℕ) (t : ℝ) :
    Tendsto (fun n => dimerTailLGT J h (Ns n) R t (xs n)) atTop (𝓝 (dimerTailZ J h r R t)) := by
  sorry

/-- **Tightness from dynamical localization** (paper: Remark `rem:dimerstrong`): if almost surely
`sup_t ∑_y |y|^p |⟨δ_y, e^{-itH}δ_0⟩|² < ∞` for some `p > 0`, then `sup_t dimerTailZ R t → 0`. -/
theorem dimerTailZ_tendsto_zero (J h : ℝ) (hJ : 0 < J) (r : ℤ) (p : ℝ) (hp : 0 < p)
    (hDL : ∀ᵐ η ∂muZ, ∃ C : ℝ, ∀ (t : ℝ) (S : Finset ℤ),
      ∑ y ∈ S, |(y : ℝ)| ^ p * ‖evoZw (dimerPotZ J h r η) t y 0‖ ^ 2 ≤ C) :
    Tendsto (fun R : ℕ => ⨆ t : ℝ, dimerTailZ J h r R t) atTop (𝓝 0) := by
  sorry

/-- **The dressed fermion stays tight** (paper: Remark `rem:dimerstrong`): given dynamical localization
of the random dimer model, `lim_{R→∞} sup_t lim_{N→∞} ⟨Ψ^{x_N}_N(t), ∑_{|y-x_N|>R} n_y Ψ^{x_N}_N(t)⟩ = 0`. -/
theorem dimer_tight_of_dynLoc (J h : ℝ) (hJ : 0 < J) (r : ℤ) (Ns : ℕ → ℕ) (xs : ∀ n, Fin (Ns n))
    (hN : ∀ n, Even (Ns n)) (hpar : ∀ n, ((xs n).val : ℤ) % 2 = r)
    (hx1 : Tendsto (fun n => (xs n).val) atTop atTop)
    (hx2 : Tendsto (fun n => Ns n - (xs n).val) atTop atTop) (p : ℝ) (hp : 0 < p)
    (hDL : ∀ᵐ η ∂muZ, ∃ C : ℝ, ∀ (t : ℝ) (S : Finset ℤ),
      ∑ y ∈ S, |(y : ℝ)| ^ p * ‖evoZw (dimerPotZ J h r η) t y 0‖ ^ 2 ≤ C) :
    ∃ L : ℕ → ℝ → ℝ, (∀ R t, Tendsto (fun n => dimerTailLGT J h (Ns n) R t (xs n)) atTop
        (𝓝 (L R t))) ∧ Tendsto (fun R : ℕ => ⨆ t : ℝ, L R t) atTop (𝓝 0) := by
  sorry

end DFL
