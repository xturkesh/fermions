import DFL.Sector

/-!
# The open-chain variant (paper: Proposition `prop:variants`(i))

The chain without boundary link is our chain with the boundary link `b_{N-1}` frozen in the state
`|+⟩` (label `0`): the Hamiltonian never flips the boundary link, so the subspace
`{x_{b_{N-1}} = 0}` is invariant. In the dual frame `x_{b_{N-1}} = ∑_j q_j + ∑_j n_j`, so this
subspace is the parity sector `∏_j q_j = (-1)^{N_f}`. Consequently, for a matter state of definite
parity `π` and links `|↑…↑⟩_z` (the boundary link in `|+⟩`), the sector law is uniform on
`{q : ∑_j q_j = π}` — density `2` with respect to the i.i.d. Bernoulli(1/2) law.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The Hamiltonian never flips the boundary link `b_{N-1}`. -/
theorem H_boundary_link (J h : ℝ) (hN : 0 < N) (p p' : Cfg N × Cfg N)
    (hp : H J h p' p ≠ 0) : p'.2 ⟨N - 1, by omega⟩ = p.2 ⟨N - 1, by omega⟩ := by
  sorry

/-- In the dual frame the boundary-link label is the total parity of charges plus particles. -/
theorem linkOf_last (hN : 0 < N) (n q : Cfg N) :
    linkOf n q ⟨N - 1, by omega⟩ = (∑ j, q j) + ∑ j, n j := by
  sorry

/-- The open-chain initial state: links `0, …, N-2` in `|↑⟩_z`, boundary link in `|+⟩`, matter in
`ψ`: amplitude `2^{-(N-1)/2} ψ(n)` on link configurations with `x_{b_{N-1}} = 0`. -/
def openInit (hN : 0 < N) (ψ : Cfg N → ℂ) : Cfg N × Cfg N → ℂ :=
  fun p => if p.2 ⟨N - 1, by omega⟩ = 0 then (((Real.sqrt 2) ^ (N - 1) : ℝ) : ℂ)⁻¹ * ψ p.1
    else 0

/-- **Open chain, sector law** (paper: Proposition `prop:variants`(i)): for a matter state of
definite parity `π`, every gauge-invariant expectation is the average over the `2^{N-1}` charge
configurations with `∑_j q_j = π` of free-fermion expectations. -/
theorem sector_average_open (J h : ℝ) (hN : 0 < N) (ψ : Cfg N → ℂ) (π : ZMod 2)
    (hψ : ∀ n, ψ n ≠ 0 → ∑ j, n j = π) (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ)
    (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (H J h) t *ᵥ openInit hN ψ) =
      (1 / 2 ^ (N - 1) : ℂ) * ∑ q ∈ univ.filter (fun q : Cfg N => ∑ j, q j = π),
        expect (sectorBlock O q) (evo (Hsec J h q) t *ᵥ ψ) := by
  sorry

/-- The open-chain sector law has `2^{N-1}` elements (`N ≥ 1`). -/
theorem card_parity_sector (hN : 0 < N) (π : ZMod 2) :
    (univ.filter (fun q : Cfg N => ∑ j, q j = π)).card = 2 ^ (N - 1) := by
  sorry

end DFL
