import DFL.LDTEntry
import DFL.Main
import DFL.AreaLaw

/-!
# Hypothesis (EC) from uniform large deviations (paper: Theorem `thm:EC`(i))

The paper's proof of Hypothesis (EC), uniformly in the volume and up to the boundary: for `x < y`,
`n = ⌊(d-1)/2⌋`, boxes of radius `n` around `x` and `y` truncated by `[0, N-1]`, and the event that no
energy is resonant for both boxes (probability `≥ 1 - 4(8n+6) C₁ e^{-c₁ n}`, `two_box`), every
eigenvector satisfies the Poisson bound at `x` or at `y` (`decay_of_not_res`), which gives
`Q(x, y) ≤ 4 e^{-κ n}` on that event.

* `corrQ_decay_of_ldt`: `E_q Q_q(x, y) ≤ C e^{-η |x-y|}` for `h_q = -J Δ + 2h diag(q)`, uniformly in
  `N`, given `LDT (-(2h/J))`.
* `ec_of_ldt`: Hypothesis (EC) holds with `F(d) = min(1, C e^{-η d})`.
* `lgt_transport_of_ldt`, `lgt_numberEntropy_of_ldt`, `lgt_areaLaw_of_ldt`: the main results of the
  paper with Hypothesis (EC) replaced by the literature input `LDT`.
-/

open Matrix Finset

noncomputable section

namespace DFL

/-- **Uniform decay of the averaged eigenfunction correlator** (paper: Theorem `thm:EC`(i)). -/
theorem corrQ_decay_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ ∀ (N : ℕ) (x y : Fin N),
      avgQ (fun q => corrQ (hq J h q) (hq_isHermitian J h q) x y) ≤
        C * Real.exp (-η * Nat.dist x.val y.val) := by
  sorry

/-- **Hypothesis (EC) holds** given the uniform large deviations (paper: Theorem `thm:EC`). -/
theorem ec_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ EC J h (fun d => min 1 (C * Real.exp (-η * d))) := by
  sorry

/-! ### The main results with (EC) replaced by the literature input -/

variable {N : ℕ}

/-- **Bounded transport** (paper: Theorem `thm:transport`) from uniform large deviations. -/
theorem lgt_transport_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ K : ℝ, ∀ (N : ℕ) (S : Finset (Fin N)) (m : ℕ) (t : ℝ),
      |(expect (matter (numOp (leftOf m))) (lgtState J h (fockVec S) t)).re -
        ((S ∩ leftOf m).card : ℝ)| ≤ K := by
  obtain ⟨C, η, -, hEC⟩ := ec_of_ldt J h hJ hL
  exact ⟨_, fun N S m t => lgt_transport J h hEC S m t⟩

/-- **Bounded number entropy** (paper: Theorem `thm:numberentropy`) from uniform large deviations. -/
theorem lgt_numberEntropy_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ K : ℝ, ∀ (N : ℕ) (m : ℕ) (S : Finset (Fin N)) (t : ℝ),
      shannon N (lgtNumDist J h (leftOf m) (fockVec S) t) ≤
        (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (K + 1 / 12)) := by
  obtain ⟨C, η, -, hEC⟩ := ec_of_ldt J h hJ hL
  exact lgt_numberEntropy J h hEC

/-- **Area law** (paper: Theorem `thm:arealaw`(a),(b)) from uniform large deviations. -/
theorem lgt_areaLaw_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ K : ℝ, ∀ (N : ℕ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (S : Finset (Fin N)) (t : ℝ),
      avgQ (fun q : Cfg N =>
        vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ fockVec S))
          (rdm_isHermitian ℓ hℓ _)) ≤ K := by
  obtain ⟨C, η, -, hEC⟩ := ec_of_ldt J h hJ hL
  exact lgt_areaLaw J h hEC

/-- **Persistent imbalance** (paper: Theorem `thm:imbalance`) from uniform large deviations. -/
theorem lgtImbalance_lower_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hh : h ≠ 0)
    (hL : LDT (-(2 * h / J))) :
    ∃ c : ℝ, 0 < c ∧ ∀ (N : ℕ), Even N → 2 ≤ N →
      ∀ ε > 0, ∀ᶠ T in Filter.atTop,
        c - ε ≤ (1 / T) * ∫ t in (0 : ℝ)..T, lgtImbalance J h N t := by
  obtain ⟨C, η, -, hEC⟩ := ec_of_ldt J h hJ hL
  exact lgtImbalance_lower J h hJ hh hEC

end DFL
