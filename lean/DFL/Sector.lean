import DFL.Duality
import DFL.Evolution

/-!
# Exact sector averaging (paper: Corollary `cor:sector`)

For the initial state `|↑…↑⟩_z ⊗ ψ` (all links polarized along `z`), the dual state is the uniform
superposition `2^{-N/2} ∑_q ψ ⊗ |q⟩` over all `2^N` charge sectors. Hence the expectation of every
gauge-invariant observable at every time is an exact average over i.i.d. Bernoulli(1/2) charges of
free-fermion expectations.
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The LGT initial state `|↑…↑⟩_z ⊗ ψ` written in the `σˣ` basis of the links:
`Ψ₀(n, x) = 2^{-N/2} ψ(n)` for every link configuration `x`. -/
def initZ (ψ : Cfg N → ℂ) : Cfg N × Cfg N → ℂ :=
  fun p => (((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ * ψ p.1

/-- The LGT state at time `t` evolved from `|↑…↑⟩_z ⊗ ψ`. -/
def lgtState (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) : Cfg N × Cfg N → ℂ :=
  evo (H J h) t *ᵥ initZ ψ

/-- Gauge invariance: `O` commutes with every Gauss operator. -/
def GaugeInv (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) : Prop :=
  ∀ j : Fin N, O * cplx (gauss j) = cplx (gauss j) * O

/-- The block of `O` in the charge sector `q`, in the dual frame. -/
def sectorBlock (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (q : Cfg N) :
    Matrix (Cfg N) (Cfg N) ℂ :=
  fun n' n => O (dualEquiv (n', q)) (dualEquiv (n, q))

/-- Gauge-invariant observables are block diagonal in the charges (dual frame). -/
theorem gaugeInv_submatrix (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) :
    O.submatrix dualEquiv dualEquiv = blockDiagonal (sectorBlock O) := by
  sorry

/-- Diagonal matter observables `D ⊗ 1` (densities, number operators and their spectral
projections) are gauge invariant. -/
theorem gaugeInv_matter (d : Cfg N → ℝ) :
    GaugeInv (cplx (diagonal d ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))) := by
  sorry

/-- The sector blocks of a diagonal matter observable `D ⊗ 1` are all equal to `D`. -/
theorem sectorBlock_matter (d : Cfg N → ℝ) (q : Cfg N) :
    sectorBlock (cplx (diagonal d ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))) q = cplx (diagonal d) := by
  sorry

/-- **Exact sector averaging** (paper: Corollary `cor:sector`): for every gauge-invariant `O` and
every time `t`, `⟨Ψ(t)|O|Ψ(t)⟩ = 2^{-N} ∑_q ⟨ψ_q(t)|O_q|ψ_q(t)⟩` with
`ψ_q(t) = e^{-it (dΓ(h_q) - h∑q)} ψ`. -/
theorem sector_average (J h : ℝ) (ψ : Cfg N → ℂ) (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ)
    (hO : GaugeInv O) (t : ℝ) :
    expect O (lgtState J h ψ t) =
      (1 / 2 ^ N : ℂ) * ∑ q : Cfg N, expect (sectorBlock O q) (evo (Hsec J h q) t *ᵥ ψ) := by
  sorry

/-- Sector averaging for diagonal matter observables, with the constant phase removed:
`⟨Ψ(t)| D ⊗ 1 |Ψ(t)⟩ = 2^{-N} ∑_q ⟨e^{-itdΓ(h_q)}ψ| D |e^{-itdΓ(h_q)}ψ⟩`. -/
theorem sector_average_matter (J h : ℝ) (ψ : Cfg N → ℂ) (d : Cfg N → ℝ) (t : ℝ) :
    expect (cplx (diagonal d ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))) (lgtState J h ψ t) =
      (1 / 2 ^ N : ℂ) * ∑ q : Cfg N,
        expect (cplx (diagonal d)) (evo (dGamma (hq J h q)) t *ᵥ ψ) := by
  sorry

end DFL
