import DFL.ECSup
import DFL.InfVolNumber
import DFL.Normalization

/-!
# Annealed localization of the full correlator on `ℤ` (paper: Corollary `cor:ECinfvol`, `X = ℤ`)

On `ℓ²(ℤ)` let `A` be a bounded operator with the matrix elements of `A_ℤ = Δ_ℤ + λ diag(q)`, `λ = -2h/J` (the normalization
of `DFL.LDTEntry`: `⟨δ_x, A δ_y⟩ = 1` for `|x - y| = 1`, `λ sgn(q_x)` for `x = y`, `0` otherwise; `IsAZ`; such an operator
exists, `exists_isAZ`). The eigenfunction correlator is
`Q^ℤ(x, y) = sup { |⟨δ_x, f(A) δ_y⟩| : |f| ≤ 1 }` (`corrQZ`), with `f(A)` the continuous functional calculus of the
self-adjoint `A` (the supremum over Borel `f` coincides with the supremum over continuous `f`, paper: proof of
Corollary `cor:ECinfvol`). Under the large-deviation input, `E Q^ℤ(x, y) ≤ C e^{-η |x - y|}` (`lintegral_corrQZ_le`).
-/

open MeasureTheory Filter Topology
open scoped InnerProductSpace

noncomputable section

namespace DFL

/-- `ℓ²(ℤ)`. -/
abbrev L2Z : Type := lp (fun _ : ℤ => ℂ) 2

/-- The basis vector `δ_x` of `ℓ²(ℤ)`. -/
def deltaZ (x : ℤ) : L2Z := lp.single 2 x 1

/-- `A` has the matrix elements of `A_ℤ = Δ_ℤ + λ diag(q)` (`q_x` read as `sgn(q_x) ∈ {±1}`). -/
def IsAZ (lam : ℝ) (q : CfgZ) (A : L2Z →L[ℂ] L2Z) : Prop :=
  ∀ x y : ℤ, ⟪deltaZ x, A (deltaZ y)⟫_ℂ =
    ((if x = y + 1 ∨ y = x + 1 then 1 else 0) + (if x = y then lam * sgn (q x) else 0) : ℝ)

/-- **The operator exists** (a bounded tridiagonal operator, `‖A‖ ≤ 2 + |λ|`). -/
theorem exists_isAZ (lam : ℝ) (q : CfgZ) : ∃ A : L2Z →L[ℂ] L2Z, IsAZ lam q A := by
  sorry

/-- The eigenfunction correlator `Q^ℤ(x, y) = sup_{|f| ≤ 1} |⟨δ_x, f(A) δ_y⟩|` (`f(A) = cfc f A`). -/
def corrQZ (A : L2Z →L[ℂ] L2Z) (x y : ℤ) : ℝ :=
  ⨆ (f : ℂ → ℂ) (_ : ∀ z, ‖f z‖ ≤ 1), ‖⟪deltaZ x, cfc f A (deltaZ y)⟫_ℂ‖

/-- **Annealed localization on `ℤ`** (paper: Corollary `cor:ECinfvol`): `E Q^ℤ(x, y) ≤ C e^{-η|x-y|}`, for the operators
`A_q` (determined by their matrix elements; no measurability is assumed, the bound is for the lower integral). -/
theorem lintegral_corrQZ_le (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ ∀ (A : CfgZ → L2Z →L[ℂ] L2Z), (∀ q, IsAZ (-(2 * h / J)) q (A q)) →
      ∀ x y : ℤ, ∫⁻ q, ENNReal.ofReal (corrQZ (A q) x y) ∂muZ ≤
        ENNReal.ofReal (C * Real.exp (-η * (x - y).natAbs)) := by
  sorry

end DFL
