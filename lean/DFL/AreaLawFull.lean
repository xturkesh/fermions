import DFL.AreaLawSharp
import DFL.CutLocal

/-!
# Area law for definite-charge initial states with the constants of the paper
(paper: Theorem `thm:arealaw`(b),(c))

For `Ψ₀ = |ξ⟩_σ ⊗ |n⟩` (links in `σˣ` eigenstates) the full matter-plus-gauge entanglement is the free-fermion
entanglement in the sector `q = q(n, ξ)` (`vonNeumann_rdmLGT_xprod`). For uniformly distributed `ξ`:

* `E_ξ S_full(t) = S_G(t) = E_q S(ρ^q_A(t))` for every `t` (`avg_full_eq_SG`);
* `sup_t S_G(t) ≤ 8 ln 2 M_2(F)` (`sup_SG_le`);
* `E_ξ sup_t S_full(t) ≤ 8 ln 2 M_2(F)` (`avg_sup_full_le`);
* `P_ξ(sup_t S_full(t) ≥ s) ≤ 8 ln 2 M_2(F)/s` for `s > 0` (`prob_sup_full_ge_le`).

Here `M_2(F) = ∑_{d ≥ 0} (d+1)² F(d)`.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- **`E_ξ S_full(t) = S_G(t)`** (paper: Theorem `thm:arealaw`(c)). -/
theorem avg_full_eq_SG (J h : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (n : Cfg N) (t : ℝ) :
    (1 / 2 ^ N : ℝ) * ∑ ξ : Cfg N,
        vonNeumann (rdmLGT ℓ hℓ (evo (H J h) t *ᵥ xprodInit ξ n)) (rdmLGT_isHermitian ℓ hℓ _) =
      avgQ (fun q : Cfg N => vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1))
        (rdm_isHermitian ℓ hℓ _)) := by
  sorry

/-- **`sup_t S_G(t) ≤ 8 ln 2 M_2(F)`** (paper: Theorem `thm:arealaw`(b)). -/
theorem sup_SG_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (ℓ : ℕ) (hℓ : ℓ ≤ N) (n : Cfg N) (t : ℝ) :
    avgQ (fun q : Cfg N => vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1))
      (rdm_isHermitian ℓ hℓ _)) ≤ 8 * Real.log 2 * ∑' d : ℕ, ((d : ℝ) + 1) ^ 2 * F d := by
  sorry

/-- **`E_ξ sup_t S_full(t) ≤ 8 ln 2 M_2(F)`** (paper: Theorem `thm:arealaw`(c)). -/
theorem avg_sup_full_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (ℓ : ℕ) (hℓ : ℓ ≤ N) (n : Cfg N) :
    (1 / 2 ^ N : ℝ) * ∑ ξ : Cfg N, ⨆ t : ℝ,
        vonNeumann (rdmLGT ℓ hℓ (evo (H J h) t *ᵥ xprodInit ξ n)) (rdmLGT_isHermitian ℓ hℓ _) ≤
      8 * Real.log 2 * ∑' d : ℕ, ((d : ℝ) + 1) ^ 2 * F d := by
  sorry

/-- **Markov bound** (paper: Theorem `thm:arealaw`(c)):
`P_ξ(sup_t S_full(t) ≥ s) ≤ 8 ln 2 M_2(F)/s`. -/
theorem prob_sup_full_ge_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (ℓ : ℕ) (hℓ : ℓ ≤ N) (n : Cfg N)
    (s : ℝ) (hs : 0 < s) :
    ((univ.filter (fun ξ : Cfg N => s ≤ ⨆ t : ℝ,
        vonNeumann (rdmLGT ℓ hℓ (evo (H J h) t *ᵥ xprodInit ξ n)) (rdmLGT_isHermitian ℓ hℓ _))).card : ℝ) /
        2 ^ N ≤
      8 * Real.log 2 * (∑' d : ℕ, ((d : ℝ) + 1) ^ 2 * F d) / s := by
  sorry

end DFL
