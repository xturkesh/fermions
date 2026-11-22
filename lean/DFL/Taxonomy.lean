import DFL.CutLocal
import DFL.Main

/-!
# The number entropy is a lower bound to the entanglement entropy
(paper: Proposition `prop:entropies`(c))

If a density matrix commutes with a "charge" (here the particle number of the region), its von
Neumann entropy is at least the Shannon entropy of the charge distribution:
`S(ρ) = S_N + ∑_k p_k S(ρ^{(k)}) ≥ S_N`. We prove the inequality directly (each eigenvalue with
weight in the charge block `k` is at most the block weight `p_k`) and apply it to
* the entanglement of an evolved Fock state (Slater determinant) in every charge sector:
  `S(ρ^q_A(t)) ≥ H(p^q_t)`;
* the full matter + gauge entanglement of the LGT state started from `|↑…↑⟩_z ⊗ |n⟩`:
  `S_full(t) ≥ S_N(t)`.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- **Block entropy bound**: if the density matrix `ρ` is block diagonal with respect to a charge
`c : ι → ℕ` (with values `≤ M`), its von Neumann entropy is at least the Shannon entropy of the
distribution of the charge, `p_k = ∑_{a : c a = k} ρ_{aa}`. -/
theorem shannon_le_vonNeumann_of_block {ι : Type*} [Fintype ι] [DecidableEq ι]
    (ρ : Matrix ι ι ℂ) (hρ : ρ.PosSemidef) (htr : ρ.trace = 1) (c : ι → ℕ) (M : ℕ)
    (hcM : ∀ a, c a ≤ M) (hblock : ∀ a b, c a ≠ c b → ρ a b = 0) :
    shannon M (fun k => (∑ a ∈ univ.filter (fun a => c a = k), ρ a a).re) ≤
      vonNeumann ρ hρ.isHermitian := by
  sorry

/-- **Slater determinants** (every charge sector): the entanglement entropy of the left segment
`{0, …, ℓ-1}` is at least the number entropy of the segment. -/
theorem numDist_shannon_le_vonNeumann_slater (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm)
    (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    shannon N (numDist (segment ℓ) (evo (dGamma A) t *ᵥ Pi.single n 1)) ≤
      vonNeumann (rdm ℓ hℓ (evo (dGamma A) t *ᵥ Pi.single n 1)) (rdm_isHermitian ℓ hℓ _) := by
  sorry

/-- **The LGT** (paper: Proposition `prop:entropies`(c) for `S_full`): for the state evolved from
`|↑…↑⟩_z ⊗ |n⟩`, the full matter + gauge entanglement entropy of the region
`A` = (sites `< ℓ`, links `< ℓ`) is at least the number entropy of the matter in `A`. -/
theorem lgt_numberEntropy_le_full (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    shannon N (lgtNumDist J h (segment ℓ) (Pi.single n 1) t) ≤
      vonNeumann (rdmLGT ℓ hℓ (lgtState J h (Pi.single n 1) t)) (rdmLGT_isHermitian ℓ hℓ _) := by
  sorry

end DFL
