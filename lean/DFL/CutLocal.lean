import DFL.AreaLaw
import DFL.ProductLaw

/-!
# Full entanglement of definite-charge states (paper: Lemma `lem:cutlocal`(c),
Theorem `thm:arealaw`(c))

Let the links start in a `σˣ`-product state `|ξ⟩_σ` and the matter in a Fock state `|n⟩`; then the
LGT state stays in the single charge sector `q = chargeOf n ξ` (`DFL.expect_xprod_fock`). For the
region `A` = (sites `< ℓ`, links `< ℓ`), the full matter + gauge entanglement entropy of the LGT
state equals the entanglement entropy of the free-fermion state `e^{-it dΓ(h_q)} |n⟩` for the modes
`< ℓ`. Reason: in the LGT state the link labels in `A` are a function of the occupations in `A`
(Gauss law), and the link labels outside `A` only depend on the particle-number parity of `A`, which
is the same for all configurations in the support (particle-number conservation). Combined with
the area law of `DFL.AreaLaw`, the average over uniformly random `ξ` of the full entanglement is
bounded uniformly in time and volume.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Splitting an LGT configuration (occupations, link labels) into the region
`A` = (sites `< ℓ`, links `< ℓ`) and its complement. -/
def splitLGT (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    Cfg N × Cfg N ≃ (Cfg ℓ × Cfg ℓ) × (Cfg (N - ℓ) × Cfg (N - ℓ)) :=
  (Equiv.prodCongr (splitCfg ℓ hℓ) (splitCfg ℓ hℓ)).trans (Equiv.prodProdProdComm _ _ _ _)

/-- Reduced density matrix of an LGT pure state `Ψ` on the region `A`. -/
def rdmLGT (ℓ : ℕ) (hℓ : ℓ ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    Matrix (Cfg ℓ × Cfg ℓ) (Cfg ℓ × Cfg ℓ) ℂ :=
  fun a a' => ∑ b : Cfg (N - ℓ) × Cfg (N - ℓ),
    Ψ ((splitLGT ℓ hℓ).symm (a, b)) * star (Ψ ((splitLGT ℓ hℓ).symm (a', b)))

theorem rdmLGT_isHermitian (ℓ : ℕ) (hℓ : ℓ ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    (rdmLGT ℓ hℓ Ψ).IsHermitian := by
  sorry

/-- The `σˣ`-product gauge state `|ξ⟩_σ` with Fock matter `|n⟩`. -/
def xprodInit (ξ n : Cfg N) : Cfg N × Cfg N → ℂ :=
  prodInit (fun j x => if x = ξ j then 1 else 0) (Pi.single n 1)

/-- **Definite-charge states** (paper: Lemma `lem:cutlocal`(c), Theorem `thm:arealaw`(c)): the full
matter + gauge entanglement of the evolved LGT state equals the entanglement of the free-fermion
state of the sector `q = chargeOf n ξ`. -/
theorem vonNeumann_rdmLGT_xprod (J h : ℝ) (ξ n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (rdmLGT ℓ hℓ (evo (H J h) t *ᵥ xprodInit ξ n)) (rdmLGT_isHermitian ℓ hℓ _) =
      vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h (chargeOf n ξ))) t *ᵥ Pi.single n 1))
        (rdm_isHermitian ℓ hℓ _) := by
  sorry

/-- **Area law for the full entanglement of definite-charge states** (paper: Theorem
`thm:arealaw`(c)): under Hypothesis (EC), the average over uniformly random link labels `ξ` of the
full entanglement entropy is bounded uniformly in time, volume, `ℓ` and the Fock state. -/
theorem lgt_areaLaw_full (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ K : ℝ, ∀ (N : ℕ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (n : Cfg N) (t : ℝ),
      (1 / 2 ^ N : ℝ) * ∑ ξ : Cfg N,
        vonNeumann (rdmLGT ℓ hℓ (evo (H J h) t *ᵥ xprodInit ξ n)) (rdmLGT_isHermitian ℓ hℓ _) ≤
          K := by
  sorry

end DFL
