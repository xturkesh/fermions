import DFL.GaussLaw
import DFL.OneBody

/-!
# Sector-diagonal operators (paper: Lemma `lem:sectornorm`)

Let `X` be a gauge-invariant operator on the LGT Hilbert space and `P_q` the sector projections of
`DFL.GaussLaw`.

* (a) `X = ∑_q P_q X P_q` and `‖X‖ = max_q ‖P_q X P_q‖` (`gaugeInv_eq_sum`, `norm_gaugeInv`).
* (b) For a density matrix `ρ` with sector weights `w_q = Tr(ρ P_q)`, `|Tr(ρ X)| ≤ ∑_q w_q ‖P_q X P_q‖`
  (`norm_trace_le`); `w_q = 2^{-N}` for `ρ = 4^{-N} 1` and for `ρ = |Ψ₀⟩⟨Ψ₀|`,
  `Ψ₀ = |↑⋯↑⟩_z ⊗ ψ` (`trace_mixed_gaussProj`, `trace_initZ_gaussProj`), hence
  `|Tr(ρ X)| ≤ E_q ‖P_q X P_q‖` (`norm_trace_mixed_le`, `norm_expect_initZ_le`).
* (c) If `X = ∑_q x_q P_q`, then `|X| = ∑_q |x_q| P_q` (the positive semidefinite square root of
  `X^* X`; `sectorDiag_abs`) and `⟨Ψ₀, |X| Ψ₀⟩ = 4^{-N} Tr |X| = E_q |x_q|`
  (`expect_sectorDiag_abs`, `trace_sectorDiag_abs`).

Norms of matrices are `ℓ²` operator norms (`Matrix.Norms.L2Operator`).
-/

open Matrix Finset
open scoped Kronecker Matrix.Norms.L2Operator ComplexOrder

noncomputable section

namespace DFL

variable {N : ℕ}

/-! ### (a) -/

/-- A gauge-invariant operator is block diagonal: `X = ∑_q P_q X P_q`. -/
theorem gaugeInv_eq_sum (X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hX : GaugeInv X) :
    X = ∑ q : Cfg N, cplx (gaussProj q) * X * cplx (gaussProj q) := by
  sorry

/-- **Norm of a gauge-invariant operator** (paper: Lemma `lem:sectornorm`(a)):
`‖X‖ = max_q ‖P_q X P_q‖`. -/
theorem norm_gaugeInv (X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hX : GaugeInv X) :
    ‖X‖ = univ.sup' univ_nonempty
      (fun q : Cfg N => ‖cplx (gaussProj q) * X * cplx (gaussProj q)‖) := by
  sorry

/-! ### (b) -/

/-- **Expectations of gauge-invariant operators** (paper: Lemma `lem:sectornorm`(b)):
`|Tr(ρ X)| ≤ ∑_q Tr(ρ P_q) ‖P_q X P_q‖` for every positive semidefinite `ρ`. -/
theorem norm_trace_le (ρ X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hρ : ρ.PosSemidef)
    (hX : GaugeInv X) :
    ‖(ρ * X).trace‖ ≤
      ∑ q : Cfg N, (ρ * cplx (gaussProj q)).trace.re * ‖cplx (gaussProj q) * X * cplx (gaussProj q)‖ := by
  sorry

/-- The sector weights of the maximally mixed state `4^{-N} 1` are `2^{-N}`. -/
theorem trace_mixed_gaussProj (q : Cfg N) :
    (((4 : ℂ) ^ N)⁻¹ • (1 : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) * cplx (gaussProj q)).trace =
      ((2 : ℂ) ^ N)⁻¹ := by
  sorry

/-- The sector weights of `Ψ₀ = |↑⋯↑⟩_z ⊗ ψ` (`ψ` a unit vector) are `2^{-N}`. -/
theorem trace_initZ_gaussProj (ψ : Cfg N → ℂ) (hψ : star ψ ⬝ᵥ ψ = 1) (q : Cfg N) :
    (vecMulVec (initZ ψ) (star (initZ ψ)) * cplx (gaussProj q)).trace = ((2 : ℂ) ^ N)⁻¹ := by
  sorry

/-- `|Tr(4^{-N} X)| ≤ E_q ‖P_q X P_q‖` for gauge-invariant `X` (paper: Lemma `lem:sectornorm`(b)). -/
theorem norm_trace_mixed_le (X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hX : GaugeInv X) :
    ‖((4 : ℂ) ^ N)⁻¹ * X.trace‖ ≤
      avgQ (fun q : Cfg N => ‖cplx (gaussProj q) * X * cplx (gaussProj q)‖) := by
  sorry

/-- `|⟨Ψ₀, X Ψ₀⟩| ≤ E_q ‖P_q X P_q‖` for gauge-invariant `X` and `Ψ₀ = |↑⋯↑⟩_z ⊗ ψ`
(paper: Lemma `lem:sectornorm`(b)). -/
theorem norm_expect_initZ_le (ψ : Cfg N → ℂ) (hψ : star ψ ⬝ᵥ ψ = 1)
    (X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hX : GaugeInv X) :
    ‖expect X (initZ ψ)‖ ≤
      avgQ (fun q : Cfg N => ‖cplx (gaussProj q) * X * cplx (gaussProj q)‖) := by
  sorry

/-! ### (c) -/

/-- **Absolute value of a sector-diagonal operator** (paper: Lemma `lem:sectornorm`(c)): for
`X = ∑_q x_q P_q`, the operator `∑_q |x_q| P_q` is positive semidefinite and squares to `X^* X`; it is
therefore `|X|`. -/
theorem sectorDiag_abs (x : Cfg N → ℂ) :
    (∑ q : Cfg N, ((‖x q‖ : ℝ) : ℂ) • cplx (gaussProj q)).PosSemidef ∧
      (∑ q : Cfg N, ((‖x q‖ : ℝ) : ℂ) • cplx (gaussProj q)) *
          (∑ q : Cfg N, ((‖x q‖ : ℝ) : ℂ) • cplx (gaussProj q)) =
        (∑ q : Cfg N, x q • cplx (gaussProj q))ᴴ * (∑ q : Cfg N, x q • cplx (gaussProj q)) := by
  sorry

/-- `⟨Ψ₀, |X| Ψ₀⟩ = E_q |x_q|` (paper: Lemma `lem:sectornorm`(c)). -/
theorem expect_sectorDiag_abs (ψ : Cfg N → ℂ) (hψ : star ψ ⬝ᵥ ψ = 1) (x : Cfg N → ℂ) :
    expect (∑ q : Cfg N, ((‖x q‖ : ℝ) : ℂ) • cplx (gaussProj q)) (initZ ψ) =
      ((avgQ (fun q : Cfg N => ‖x q‖) : ℝ) : ℂ) := by
  sorry

/-- `4^{-N} Tr |X| = E_q |x_q|` (paper: Lemma `lem:sectornorm`(c)). -/
theorem trace_sectorDiag_abs (x : Cfg N → ℂ) :
    ((4 : ℂ) ^ N)⁻¹ * (∑ q : Cfg N, ((‖x q‖ : ℝ) : ℂ) • cplx (gaussProj q)).trace =
      ((avgQ (fun q : Cfg N => ‖x q‖) : ℝ) : ℂ) := by
  sorry

end DFL
