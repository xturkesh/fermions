import DFL.Sector

/-!
# Superposition versus ensemble (paper: Proposition `prop:entropies`(d))

The initial state `|↑…↑⟩_z ⊗ |n⟩` is a coherent superposition of all `2^N` charge sectors, while
the ensemble `ρ^ens = 2^{-N} 1_links ⊗ |n⟩⟨n|` is the incoherent mixture of the definite-charge
states `|ξ⟩_σ ⊗ |n⟩`. Every gauge-invariant observable, and every product of gauge-invariant
observables at arbitrary times (multi-time correlators), has the same expectation in the two. Only
entanglement-type quantities, which are not expectations of gauge-invariant observables, can
distinguish them (the "superposition vs sampling" point of the Google experiment).
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Heisenberg-picture observable `O(t) = e^{itH} O e^{-itH}`. -/
def heis (J h : ℝ) (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (t : ℝ) :
    Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  (evo (H J h) t)ᴴ * O * evo (H J h) t

/-- The ensemble `ρ^ens = 2^{-N} 1_links ⊗ |n⟩⟨n|` (the dephased initial state). -/
def ensemble (n : Cfg N) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  (1 / 2 ^ N : ℂ) • diagonal (fun p => if p.1 = n then 1 else 0)

/-- Products of gauge-invariant observables are gauge invariant. -/
theorem gaugeInv_mul {O P : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ} (hO : GaugeInv O)
    (hP : GaugeInv P) : GaugeInv (O * P) := by
  sorry

/-- Heisenberg evolution preserves gauge invariance. -/
theorem gaugeInv_heis (J h : ℝ) {O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ}
    (hO : GaugeInv O) (t : ℝ) : GaugeInv (heis J h O t) := by
  sorry

/-- **Superposition = ensemble** for one gauge-invariant observable. -/
theorem expect_initZ_eq_trace_ensemble (n : Cfg N)
    (X : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hX : GaugeInv X) :
    expect X (initZ (Pi.single n 1)) = trace (ensemble n * X) := by
  sorry

/-- **Superposition = ensemble for all multi-time correlators of gauge-invariant observables**
(paper: Proposition `prop:entropies`(d), eq. `eq:main-ens`). -/
theorem correlator_eq_ensemble (J h : ℝ) (n : Cfg N)
    (Os : List (Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ × ℝ)) (hOs : ∀ p ∈ Os, GaugeInv p.1) :
    expect ((Os.map (fun p => heis J h p.1 p.2)).prod) (initZ (Pi.single n 1)) =
      trace (ensemble n * (Os.map (fun p => heis J h p.1 p.2)).prod) := by
  sorry

end DFL
