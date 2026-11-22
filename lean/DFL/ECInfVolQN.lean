import DFL.ECInfVolQ

/-!
# Annealed localization of the full correlator on the half-line (paper: Corollary `cor:ECinfvol`, `X = ℕ`)

On `ℓ²(ℕ)` (sites `0, 1, 2, …`, Dirichlet condition at `-1`; the paper's `ℕ = {1, 2, …}` shifted by one) let `A` be the
bounded operator with the matrix elements of `A_ℕ = Δ_ℕ + λ diag(q)` (`IsAN`; it exists, `exists_isAN`), with i.i.d.
uniform charges `q : ℕ → ZMod 2` (`muN`). The correlator `Q^ℕ(x, y) = sup_{|f| ≤ 1} |⟨δ_x, f(A) δ_y⟩|` (`corrQN`, continuous
functional calculus) satisfies `E Q^ℕ(x, y) ≤ C e^{-η|x - y|}` under the large-deviation input (`lintegral_corrQN_le`). The
finite windows `[0, L)` are exactly the chains of `DFL.Basic` with the restricted charges.
-/

open MeasureTheory Filter Topology
open scoped InnerProductSpace

noncomputable section

namespace DFL

/-- Charge configurations on the half-line. -/
abbrev CfgN := ℕ → ZMod 2

/-- I.i.d. uniform charges on the half-line. -/
def muN : Measure CfgN := Measure.infinitePi (fun _ : ℕ => unifZ2)

/-- `ℓ²(ℕ)`. -/
abbrev L2N : Type := lp (fun _ : ℕ => ℂ) 2

/-- The basis vector `δ_x` of `ℓ²(ℕ)`. -/
def deltaN (x : ℕ) : L2N := lp.single 2 x 1

/-- `A` has the matrix elements of `A_ℕ = Δ_ℕ + λ diag(q)` (Dirichlet at `-1`). -/
def IsAN (lam : ℝ) (q : CfgN) (A : L2N →L[ℂ] L2N) : Prop :=
  ∀ x y : ℕ, ⟪deltaN x, A (deltaN y)⟫_ℂ =
    ((if x = y + 1 ∨ y = x + 1 then 1 else 0) + (if x = y then lam * sgn (q x) else 0) : ℝ)

/-- **The operator exists**. -/
theorem exists_isAN (lam : ℝ) (q : CfgN) : ∃ A : L2N →L[ℂ] L2N, IsAN lam q A := by
  sorry

/-- The eigenfunction correlator `Q^ℕ(x, y) = sup_{|f| ≤ 1} |⟨δ_x, f(A) δ_y⟩|`. -/
def corrQN (A : L2N →L[ℂ] L2N) (x y : ℕ) : ℝ :=
  ⨆ (f : ℂ → ℂ) (_ : ∀ z, ‖f z‖ ≤ 1), ‖⟪deltaN x, cfc f A (deltaN y)⟫_ℂ‖

/-- **Annealed localization on the half-line** (paper: Corollary `cor:ECinfvol`, `X = ℕ`):
`E Q^ℕ(x, y) ≤ C e^{-η|x-y|}`. -/
theorem lintegral_corrQN_le (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ ∀ (A : CfgN → L2N →L[ℂ] L2N), (∀ q, IsAN (-(2 * h / J)) q (A q)) →
      ∀ x y : ℕ, ∫⁻ q, ENNReal.ofReal (corrQN (A q) x y) ∂muN ≤
        ENNReal.ofReal (C * Real.exp (-η * Nat.dist x y)) := by
  sorry

end DFL
