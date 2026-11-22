import DFL.Entropies

/-!
# Sharpness of the duality defect (paper: Lemma `lem:cutlocal`(d))

Two sites, cut after the first site (`N = 2`, `ℓ = 1`). The same function
`exProd (n, y) = 2^{-1/2} [n = 0] [y₁ = 0]`, read either in the dual variables `(n, q)` or in the
original variables `(n, x)`, is a product state across the cut. Its image under `U`
(`(UΦ)(n, x) = Φ(n, chargeOf n x)`) and under `U†` (`(U†Ψ)(n, q) = Ψ(n, linkOf n q)`) is a Bell state of
the two links, resp. of the two charges, with entanglement entropy `ln 2`. So the bound `ln 2` of
`DFL.cutlocal_defect` is attained in both directions.
-/

open Matrix Finset

noncomputable section

namespace DFL

/-- `2^{-1/2} |0, 0⟩ ⊗ (|y₀ = 0⟩ + |y₀ = 1⟩) ⊗ |y₁ = 0⟩`. -/
def exProd : Cfg 2 × Cfg 2 → ℂ :=
  fun p => if p.1 = 0 ∧ p.2 1 = 0 then ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ else 0

/-- `exProd` is a product state across the cut: zero entanglement entropy. -/
theorem exProd_entropy :
    vonNeumann (rdmLGT 1 (by norm_num) exProd) (rdmLGT_isHermitian 1 (by norm_num) _) = 0 := by
  sorry

/-- **`U` creates entanglement `ln 2`** (paper: Lemma `lem:cutlocal`(d), first example). -/
theorem exProd_U_entropy :
    vonNeumann (rdmLGT 1 (by norm_num) (fun p => exProd (dualEquiv.symm p)))
      (rdmLGT_isHermitian 1 (by norm_num) _) = Real.log 2 := by
  sorry

/-- **`U†` creates entanglement `ln 2`** (paper: Lemma `lem:cutlocal`(d), second example). -/
theorem exProd_Udag_entropy :
    vonNeumann (rdmLGT 1 (by norm_num) (fun p => exProd (dualEquiv p)))
      (rdmLGT_isHermitian 1 (by norm_num) _) = Real.log 2 := by
  sorry

end DFL
