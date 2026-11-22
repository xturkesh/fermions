import DFL.VNEntropy
import DFL.Taxonomy

/-!
# Entropy taxonomy (paper: Lemma `lem:cutlocal`(b), Proposition `prop:entropies`(a)–(c))

For the initial state `|↑⋯↑⟩_z ⊗ ψ` and the region `A` = (sites `< ℓ`, links `< ℓ`):

* `cutlocal_defect`: the duality `U` is local across the cut up to one controlled gate, so for every
  normalized state `Φ` of the dual variables, `|S_A(UΦ) - S_{A'}(Φ)| ≤ ln 2` (pinching argument);
  `lgt_duality_defect` is the case `Φ = Φ(t) = U†Ψ(t)` (paper: `prop:entropies`(a)).
* `ρ^q_A(t)` (`secRdm`) is the reduced state of the sector-`q` free-fermion state
  `ψ_q(t) = Γ(e^{-ith_q}) ψ`, `ρ_mix(t) = E_q ρ^q_A(t)` (`mixRdm`), `S_G(t) = E_q S(ρ^q_A(t))`
  (`condEntropy`), and `H_t` (`sectorLawEntropy`) is the Shannon entropy of the law of the random
  matrix `q ↦ ρ^q_A(t)`. Then `S_G ≤ S_mix ≤ min(S_G + H_t, ℓ ln 2)`, `H_t ≤ N ln 2`, `H_0 = 0`, and
  `ρ_mix` is the marginal of `ρ_{A'}(Φ(t))` on the modes (paper: `prop:entropies`(b)).
* For Fock `ψ`: `p_t = E_q p^q_t`; the number entropy `S_N` is a lower bound for `S_mix` and for the
  entropy `S_f` of the bare fermion modes (and for `S_full`, `DFL.lgt_numberEntropy_le_full`);
  `S_G ≥ E_q H(p^q_t)`, and `E_q H(p^q_t) ≤ S_N ≤ E_q H(p^q_t) + H_t` (paper: `prop:entropies`(c)).
-/

open Matrix Finset
open scoped ComplexOrder

noncomputable section

namespace DFL

variable {N : ℕ}

/-! ### Cut locality of the duality -/

/-- **Duality defect** (paper: Lemma `lem:cutlocal`(b)): for every normalized state `Φ` of the dual
variables `(n, q)`, the entanglement entropy of `A` in `UΦ`, `(UΦ)(n, x) = Φ(n, chargeOf n x)`, and
that of `A'` = (modes `< ℓ`, charges `< ℓ`) in `Φ` differ by at most `ln 2`. -/
theorem cutlocal_defect (ℓ : ℕ) (hℓ : ℓ ≤ N) (Φ : Cfg N × Cfg N → ℂ)
    (hΦ : ∑ p, Complex.normSq (Φ p) = 1) :
    |vonNeumann (rdmLGT ℓ hℓ (fun p => Φ (dualEquiv.symm p))) (rdmLGT_isHermitian ℓ hℓ _) -
      vonNeumann (rdmLGT ℓ hℓ Φ) (rdmLGT_isHermitian ℓ hℓ _)| ≤ Real.log 2 := by
  sorry

/-- **Duality defect for the quench** (paper: Proposition `prop:entropies`(a)):
`|S_full(t) - S_{A'}(Φ(t))| ≤ ln 2`, with `Φ(t) = U†Ψ(t)`. -/
theorem lgt_duality_defect (J h : ℝ) (ψ : Cfg N → ℂ) (hψ : ∑ n, Complex.normSq (ψ n) = 1)
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    |vonNeumann (rdmLGT ℓ hℓ (lgtState J h ψ t)) (rdmLGT_isHermitian ℓ hℓ _) -
      vonNeumann (rdmLGT ℓ hℓ (fun p => lgtState J h ψ t (dualEquiv p)))
        (rdmLGT_isHermitian ℓ hℓ _)| ≤ Real.log 2 := by
  sorry

/-! ### Conditional versus mixed entropy -/

/-- The sector-`q` reduced state `ρ^q_A(t)` of `Γ(e^{-ith_q}) ψ` on the modes `< ℓ`. -/
def secRdm (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (q : Cfg N) :
    Matrix (Cfg ℓ) (Cfg ℓ) ℂ :=
  rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ ψ)

/-- `ρ_mix(t) = E_q ρ^q_A(t)`. -/
def mixRdm (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) : Matrix (Cfg ℓ) (Cfg ℓ) ℂ :=
  ((1 / 2 ^ N : ℝ) : ℂ) • ∑ q, secRdm J h ψ t ℓ hℓ q

/-- `ρ_mix(t)` is Hermitian. -/
theorem mixRdm_isHermitian (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    (mixRdm J h ψ t ℓ hℓ).IsHermitian := by
  sorry

/-- `S_G(t) = E_q S(ρ^q_A(t))`, the entanglement conditional on all Gauss charges. -/
def condEntropy (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) : ℝ :=
  avgQ (fun q => vonNeumann (secRdm J h ψ t ℓ hℓ q) (rdm_isHermitian ℓ hℓ _))

open Classical in
/-- `H_t`: the Shannon entropy of the law of the random density matrix `q ↦ ρ^q_A(t)` under the
uniform law of `q`. -/
def sectorLawEntropy (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) : ℝ :=
  ∑ σ ∈ univ.image (secRdm J h ψ t ℓ hℓ),
    Real.negMulLog (((univ.filter (fun q => secRdm J h ψ t ℓ hℓ q = σ)).card : ℝ) / 2 ^ N)

/-- `ρ_mix(t)` is the marginal on the modes `< ℓ` of the reduced state `ρ_{A'}(Φ(t))` of the dual
state `Φ(t) = U†Ψ(t)` (paper: Proposition `prop:entropies`(b)). -/
theorem mixRdm_eq_marginal (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N)
    (a a' : Cfg ℓ) :
    mixRdm J h ψ t ℓ hℓ a a' =
      ∑ c : Cfg ℓ, rdmLGT ℓ hℓ (fun p => lgtState J h ψ t (dualEquiv p)) (a, c) (a', c) := by
  sorry

/-- `S_G ≤ S_mix` (paper: Proposition `prop:entropies`(b)). -/
theorem condEntropy_le_mix (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    condEntropy J h ψ t ℓ hℓ ≤ vonNeumann (mixRdm J h ψ t ℓ hℓ) (mixRdm_isHermitian J h ψ t ℓ hℓ) := by
  sorry

/-- `S_mix ≤ S_G + H_t` (paper: Proposition `prop:entropies`(b)). -/
theorem mix_le_cond_add (J h : ℝ) (ψ : Cfg N → ℂ) (hψ : ∑ n, Complex.normSq (ψ n) = 1) (t : ℝ)
    (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (mixRdm J h ψ t ℓ hℓ) (mixRdm_isHermitian J h ψ t ℓ hℓ) ≤
      condEntropy J h ψ t ℓ hℓ + sectorLawEntropy J h ψ t ℓ hℓ := by
  sorry

/-- `S_mix ≤ ℓ ln 2` (paper: Proposition `prop:entropies`(b)). -/
theorem mix_le_log (J h : ℝ) (ψ : Cfg N → ℂ) (hψ : ∑ n, Complex.normSq (ψ n) = 1) (t : ℝ)
    (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (mixRdm J h ψ t ℓ hℓ) (mixRdm_isHermitian J h ψ t ℓ hℓ) ≤ ℓ * Real.log 2 := by
  sorry

/-- `H_t ≤ N ln 2` (paper: Proposition `prop:entropies`(b)). -/
theorem sectorLawEntropy_le (J h : ℝ) (ψ : Cfg N → ℂ) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    sectorLawEntropy J h ψ t ℓ hℓ ≤ N * Real.log 2 := by
  sorry

/-- `H_0 = 0` (paper: Proposition `prop:entropies`(b)). -/
theorem sectorLawEntropy_zero (J h : ℝ) (ψ : Cfg N → ℂ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    sectorLawEntropy J h ψ 0 ℓ hℓ = 0 := by
  sorry

/-! ### Number entropy as a lower bound -/

/-- `p_t = E_q p^q_t`: the LGT number distribution is the average of the sector distributions
(paper: Corollary `cor:sector`(c), used in Proposition `prop:entropies`(c)). -/
theorem lgtNumDist_eq_avg (J h : ℝ) (R : Finset (Fin N)) (ψ : Cfg N → ℂ) (t : ℝ) (k : ℕ) :
    lgtNumDist J h R ψ t k = avgQ (fun q => numDist R (evo (dGamma (hq J h q)) t *ᵥ ψ) k) := by
  sorry

/-- `S_N ≤ S_mix` for Fock initial states (paper: Proposition `prop:entropies`(c)). -/
theorem numberEntropy_le_mix (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    shannon N (lgtNumDist J h (segment ℓ) (Pi.single n 1) t) ≤
      vonNeumann (mixRdm J h (Pi.single n 1) t ℓ hℓ)
        (mixRdm_isHermitian J h (Pi.single n 1) t ℓ hℓ) := by
  sorry

/-- The reduced state `ρ^f_A` of the bare fermion modes `< ℓ`, with all links traced out. -/
def rdmBare (ℓ : ℕ) (hℓ : ℓ ≤ N) (Ψ : Cfg N × Cfg N → ℂ) : Matrix (Cfg ℓ) (Cfg ℓ) ℂ :=
  fun a a' => ∑ c : Cfg ℓ, rdmLGT ℓ hℓ Ψ (a, c) (a', c)

/-- `ρ^f_A` is Hermitian. -/
theorem rdmBare_isHermitian (ℓ : ℕ) (hℓ : ℓ ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    (rdmBare ℓ hℓ Ψ).IsHermitian := by
  sorry

/-- `S_N ≤ S_f` for Fock initial states (paper: Proposition `prop:entropies`(c)). -/
theorem numberEntropy_le_bare (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    shannon N (lgtNumDist J h (segment ℓ) (Pi.single n 1) t) ≤
      vonNeumann (rdmBare ℓ hℓ (lgtState J h (Pi.single n 1) t))
        (rdmBare_isHermitian ℓ hℓ _) := by
  sorry

/-- `S_G ≥ E_q H(p^q_t)` for Fock initial states (paper: Proposition `prop:entropies`(c)). -/
theorem avg_sector_shannon_le_cond (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    avgQ (fun q => shannon N (numDist (segment ℓ) (evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1))) ≤
      condEntropy J h (Pi.single n 1) t ℓ hℓ := by
  sorry

/-- `E_q H(p^q_t) ≤ S_N(t)` (concavity; paper: Proposition `prop:entropies`(c)). -/
theorem avg_sector_shannon_le_numberEntropy (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    avgQ (fun q => shannon N (numDist (segment ℓ) (evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1))) ≤
      shannon N (lgtNumDist J h (segment ℓ) (Pi.single n 1) t) := by
  sorry

/-- `S_N(t) ≤ E_q H(p^q_t) + H_t` (classical mixing bound, sectors grouped by the value of
`ρ^q_A(t)`, which determines `p^q_t`; paper: Proposition `prop:entropies`(c)). -/
theorem numberEntropy_le_avg_add (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    shannon N (lgtNumDist J h (segment ℓ) (Pi.single n 1) t) ≤
      avgQ (fun q => shannon N (numDist (segment ℓ) (evo (dGamma (hq J h q)) t *ᵥ Pi.single n 1))) +
        sectorLawEntropy J h (Pi.single n 1) t ℓ hℓ := by
  sorry

end DFL
