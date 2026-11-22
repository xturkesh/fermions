import DFL.SlaterExact
import DFL.BinEntropy
import DFL.ECProof

/-!
# The area law with the constants of the paper (paper: Lemma `lem:annealed`, Theorem `thm:arealaw`(a),(b))

* `vonNeumann_slater_sharp` (paper: `thm:arealaw`(a)): `S(ρ_A) ≤ 2 ln 2 ∑_{x ∈ A, y ∉ A} |K_{xy}|`
  for every Slater determinant, from `S(ρ_A) = ∑ h(ν_i)`, `h(ν) ≤ 2 ln 2 √(ν(1-ν))` and the trace-norm
  bound.
* `annealed_K`, `annealed_cut` (paper: `lem:annealed`(i),(ii)): under (EC),
  `E sup_t |K_q(t)_{xy}| ≤ 2 ∑_{2k ≥ |x-y|} F(k)` and `E sup_t ∑_{x ∈ A, y ∉ A} |K_q(t)_{xy}| ≤ Θ(F)`,
  `Θ(F) = 2 ∑_{k ≥ 1} k (2k+1) F(k)`.
* `lgt_areaLaw_sharp` (paper: `thm:arealaw`(b)): `E_q sup_t S(ρ^q_A(t)) ≤ 2 ln 2 Θ(F) ≤ 8 ln 2 M_2(F)`.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- `Θ(F) = 2 ∑_{k ≥ 1} k (2k+1) F(k)`. -/
def thetaF (F : ℕ → ℝ) : ℝ := 2 * ∑' k : ℕ, (k : ℝ) * (2 * k + 1) * F k

/-- **Area law for Slater determinants with the constant `2 ln 2`** (paper: Theorem
`thm:arealaw`(a)). -/
theorem vonNeumann_slater_sharp (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (rdm ℓ hℓ (evo (dGamma A) t *ᵥ fockVec S)) (rdm_isHermitian ℓ hℓ _) ≤
      2 * Real.log 2 * ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ‖corrMat A S t x y‖ := by
  sorry

/-- **Annealed correlation bound** (paper: Lemma `lem:annealed`(i)):
`E sup_t |K_q(t)_{xy}| ≤ 2 ∑_{k : 2k ≥ |x-y|} F(k)`. -/
theorem annealed_K (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (x y : Fin N) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ, ‖corrMat (hq J h q) S t x y‖) ≤
      2 * ∑' k : ℕ, (if Nat.dist x.val y.val ≤ 2 * k then F k else 0) := by
  sorry

/-- **Annealed cut bound** (paper: Lemma `lem:annealed`(ii)):
`E sup_t ∑_{x ∈ A, y ∉ A} |K_q(t)_{xy}| ≤ Θ(F)` for a left segment `A`. -/
theorem annealed_cut (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (ℓ : ℕ) :
    avgQ (fun q : Cfg N =>
      ⨆ t : ℝ, ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ‖corrMat (hq J h q) S t x y‖) ≤
      thetaF F := by
  sorry

/-- **Sector-averaged area law with the constant of the paper** (paper: Theorem `thm:arealaw`(b)):
`E_q sup_t S(ρ^q_A(t)) ≤ 2 ln 2 Θ(F)`, uniformly in `N`, `ℓ` and the Fock state. -/
theorem lgt_areaLaw_sharp (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (ℓ : ℕ) (hℓ : ℓ ≤ N)
    (S : Finset (Fin N)) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ,
      vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ fockVec S)) (rdm_isHermitian ℓ hℓ _)) ≤
      2 * Real.log 2 * thetaF F := by
  sorry

/-- `Θ(F) ≤ 4 M_2(F)`, `M_2(F) = ∑_{d ≥ 0} (d+1)² F(d)`. -/
theorem thetaF_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) :
    thetaF F ≤ 4 * ∑' d : ℕ, ((d : ℝ) + 1) ^ 2 * F d := by
  sorry

end DFL
