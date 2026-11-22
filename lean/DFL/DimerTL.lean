import DFL.Evolution

/-!
# Finite propagation speed for one-dimensional Schrödinger operators (paper: Lemma `lem:dimerTL`)

For an interval `Λ = [a, b] ∋ 0` and a real potential `w`, let `H_Λ = -Δ_Λ + diag(w)` (Dirichlet ends).

* (i) `|⟨δ_y, e^{-isH_Λ} δ_0⟩| ≤ ∑_{k ≥ |y|} (2|s|)^k / k! ≤ (2|s|)^{|y|}/|y|! · e^{2|s|}`, uniformly in the
  potential (`dimerTL_i`, `dimerTL_i'`);
* (ii) for `Λ ⊆ Λ' = [a', b']` and `|w| ≤ W` on `Λ'`,
  `|⟨δ_y, (e^{-isH_Λ} - e^{-isH_{Λ'}}) δ_0⟩| ≤ 2 ∑_{k ≥ R} ((2 + W)|s|)^k / k!` with
  `R = dist(0, ℤ \ Λ)` (`dimerTL_ii`; here `Λ'` is finite);
* (iii) for `0 < p ≤ 2`, `∑_{y ∈ Λ} |y|^p |⟨δ_y, e^{-isH_Λ} δ_0⟩|² ≤ (2|s|)^p` (`dimerTL_iii`).
-/

open Matrix

noncomputable section

namespace DFL

/-- `H_Λ = -Δ_Λ + diag(w)` on `ℓ²([a, b])`. -/
def schrH (a b : ℤ) (w : ℤ → ℝ) : Matrix (Set.Icc a b) (Set.Icc a b) ℝ := fun x y =>
  if x.1 = y.1 then w x.1 else if |x.1 - y.1| = 1 then -1 else 0

/-- The origin as a site of `[a, b]`. -/
def origin (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) : Set.Icc a b := ⟨0, ha, hb⟩

/-- **Finite speed, uniformly in the potential** (paper: Lemma `lem:dimerTL`(i)). -/
theorem dimerTL_i (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (w : ℤ → ℝ) (s : ℝ) (y : Set.Icc a b) :
    ‖evo (schrH a b w) s y (origin a b ha hb)‖ ≤
      ∑' k : ℕ, if (y.1).natAbs ≤ k then (2 * |s|) ^ k / (k.factorial : ℝ) else 0 := by
  sorry

/-- The tail bound in (i): `∑_{k ≥ m} c^k/k! ≤ (c^m / m!) e^c`. -/
theorem dimerTL_i' (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (w : ℤ → ℝ) (s : ℝ) (y : Set.Icc a b) :
    ‖evo (schrH a b w) s y (origin a b ha hb)‖ ≤
      (2 * |s|) ^ (y.1).natAbs / ((y.1).natAbs.factorial : ℝ) * Real.exp (2 * |s|) := by
  sorry

/-- **Thermodynamic limit** (paper: Lemma `lem:dimerTL`(ii), finite `Λ'`). -/
theorem dimerTL_ii (a b a' b' : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (ha' : a' ≤ a) (hb' : b ≤ b')
    (w : ℤ → ℝ) (W : ℝ) (hw : ∀ x, a' ≤ x → x ≤ b' → |w x| ≤ W) (s : ℝ) (y : Set.Icc a b) :
    ‖evo (schrH a b w) s y (origin a b ha hb) -
        evo (schrH a' b' w) s ⟨y.1, ha'.trans y.2.1, y.2.2.trans hb'⟩
          (origin a' b' (ha'.trans ha) (hb.trans hb'))‖ ≤
      2 * ∑' k : ℕ, if (min b (-a)).toNat + 1 ≤ k then ((2 + W) * |s|) ^ k / (k.factorial : ℝ) else 0 := by
  sorry

/-- **Ballistic bound** (paper: Lemma `lem:dimerTL`(iii)). -/
theorem dimerTL_iii (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (w : ℤ → ℝ) (s : ℝ) (p : ℝ) (hp0 : 0 < p)
    (hp2 : p ≤ 2) :
    ∑ y : Set.Icc a b, |((y.1 : ℤ) : ℝ)| ^ p * ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2 ≤
      (2 * |s|) ^ p := by
  sorry

end DFL
