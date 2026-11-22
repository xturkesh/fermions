import DFL.Transfer

/-!
# Resonant energies of a box (paper: Lemma `lem:badset`, Eq. `eq:Bad`)

For a box `B = [a, a+L-1]`, a site `x = a + i ∈ B`, a choice of ends (`endL`: the left end `a`,
`endR`: the right end `a+L-1`) and a threshold `θ = e^{-κn}`, the resonant set is
`{E ∈ I_λ : |p_e(E)| ≥ θ |P_B(E)| for some selected end e}`, with `P_B = dirL V E a (L+1)`,
`p_b = dirL V E a (i+1)` (right end) and `p_a = dirR V E (a+L-1) (L-i)` (left end).

* `decay_of_not_res`: outside the resonant set every solution decays across the box:
  `|φ(x)| ≤ θ (|φ(a-1)| + |φ(a+L)|)` (paper: Lemma `lem:badset`(a) with Eq. `eq:poisson`).
* `exists_resEnd`, `card_resEnds`: every resonant energy `E` lies in an interval `[z, E]` of resonant
  energies whose left end `z` belongs to an explicit set of at most `4L + 2` energies, determined by
  the potential on the box (paper: Lemma `lem:badset`(b)).
-/

open Matrix Finset Polynomial

noncomputable section

namespace DFL

/-- `I_λ = [-2-|λ|, 2+|λ|]`, which contains the spectrum of `Δ + λ diag(q)`. -/
def Ilam (lam : ℝ) : Set ℝ := Set.Icc (-2 - |lam|) (2 + |lam|)

/-- The resonant set of `(x = a + i, B = [a, a+L-1], ends)` at threshold `θ` (paper: Eq. `eq:Bad`). -/
def resSet (V : ℤ → ℝ) (lam : ℝ) (a : ℤ) (L i : ℕ) (endL endR : Bool) (θ : ℝ) : Set ℝ :=
  {E | E ∈ Ilam lam ∧
    ((endR = true ∧ θ * |dirL V E a (L + 1)| ≤ |dirL V E a (i + 1)|) ∨
      (endL = true ∧ θ * |dirL V E a (L + 1)| ≤ |dirR V E (a + L - 1) (L - i)|))}

/-- The candidate left endpoints of the components of the resonant set: the two ends of `I_λ` and the
real zeros of `p_b² - θ² P_B²` and of `p_a² - θ² P_B²`. -/
def resEnds (V : ℤ → ℝ) (lam : ℝ) (a : ℤ) (L i : ℕ) (θ : ℝ) : Finset ℝ :=
  {-2 - |lam|, 2 + |lam|} ∪
    ((dirLpoly V a (i + 1)) ^ 2 - C (θ ^ 2) * (dirLpoly V a (L + 1)) ^ 2).roots.toFinset ∪
    ((dirRpoly V (a + L - 1) (L - i)) ^ 2 - C (θ ^ 2) * (dirLpoly V a (L + 1)) ^ 2).roots.toFinset

/-- **Decay outside the resonant set**: if at least one end is selected, `E ∈ I_λ` is not resonant,
`φ` solves the equation on the box, and `φ` vanishes beyond every unselected end, then
`|φ(x)| ≤ θ (|φ(a-1)| + |φ(a+L)|)`. -/
theorem decay_of_not_res (V : ℤ → ℝ) (lam : ℝ) (a : ℤ) (L i : ℕ) (endL endR : Bool) (θ : ℝ)
    (hi : i < L) (hend : endL = true ∨ endR = true) (E : ℝ) (hE : E ∈ Ilam lam)
    (hres : E ∉ resSet V lam a L i endL endR θ) (φ : ℤ → ℝ)
    (hφ : ∀ z, a ≤ z → z ≤ a + L - 1 → φ (z + 1) + φ (z - 1) + V z * φ z = E * φ z)
    (hL : endL = false → φ (a - 1) = 0) (hR : endR = false → φ (a + L) = 0) :
    |φ (a + i)| ≤ θ * (|φ (a - 1)| + |φ (a + L)|) := by
  sorry

/-- The resonant set depends only on the potential on the box. -/
theorem resSet_congr (V V' : ℤ → ℝ) (lam : ℝ) (a : ℤ) (L i : ℕ) (endL endR : Bool) (θ : ℝ)
    (hi : i < L) (h : ∀ j : ℕ, j < L → V (a + j) = V' (a + j)) :
    resSet V lam a L i endL endR θ = resSet V' lam a L i endL endR θ := by
  sorry

/-- The candidate endpoints depend only on the potential on the box. -/
theorem resEnds_congr (V V' : ℤ → ℝ) (lam : ℝ) (a : ℤ) (L i : ℕ) (θ : ℝ) (hi : i < L)
    (h : ∀ j : ℕ, j < L → V (a + j) = V' (a + j)) :
    resEnds V lam a L i θ = resEnds V' lam a L i θ := by
  sorry

/-- There are at most `4L + 2` candidate endpoints. -/
theorem card_resEnds (V : ℤ → ℝ) (lam : ℝ) (a : ℤ) (L i : ℕ) (θ : ℝ) (hi : i < L) :
    (resEnds V lam a L i θ).card ≤ 4 * L + 2 := by
  sorry

/-- **Components of the resonant set** (paper: Lemma `lem:badset`(b)): every resonant energy `E`
lies in an interval `[z, E]` of resonant energies with `z` a candidate endpoint. -/
theorem exists_resEnd (V : ℤ → ℝ) (lam : ℝ) (a : ℤ) (L i : ℕ) (endL endR : Bool) (θ : ℝ)
    (hθ : 0 < θ) (hi : i < L) (E : ℝ) (hE : E ∈ resSet V lam a L i endL endR θ) :
    ∃ z ∈ resEnds V lam a L i θ, z ≤ E ∧ Set.Icc z E ⊆ resSet V lam a L i endL endR θ := by
  sorry

end DFL
