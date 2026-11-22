import DFL.TwoD

/-!
# Two dimensions: flux sectors (paper: Lemma `lem:2dgauge`(i), Proposition `prop:2dlaws`(a),(b))

The ℤ₂ gauge theory `FinGraph.ham` on an arbitrary finite graph, written in the `σˣ` basis of the
edges (`DFL/Graph.lean`). Link signs `s ∈ {±1}^E` are encoded as `s : Cfg M` (`sgn (s e)`); the
corresponding one-body operator is `h_{q,s} = K_s + 2h diag(q)` (`FinGraph.hqs`, `DFL/TwoD.lean`).

* `flux_eq_iff_gauge` (paper: `lem:2dgauge`(i)): two sign configurations have the same fluxes
  (the same `±1` product around every cycle, i.e. the same values `∑_e s_e z_e` on every
  `z ∈ ker ∂`) iff they differ by a gauge transformation `s ↦ s + δθ`, `(δθ)_e = θ_{s(e)} + θ_{t(e)}`.
* `sector_average_flux` (paper: `prop:2dlaws`(a)): for the `σᶻ`-product gauge state `|s⟩_z` and Fock
  matter `|n⟩`, the flux is that of `s`, and the dynamics of every matter observable diagonal in the
  occupations is the uniform average over the reachable charge sectors `q ∈ n + im ∂` of free-fermion
  dynamics with `h_{q,s}` (for `s = +` this is `FinGraph.sector_average_matter`).
* `sector_average_xpol` (paper: `prop:2dlaws`(b)): for the `σˣ`-polarized gauge state `|+⋯+⟩` the
  charges are fixed (`q⁰ = n`, i.e. `G_i = (-1)^{n_i}`) and the dynamics is the uniform average over
  all sign configurations `s` of free-fermion dynamics with `h_{q⁰,s}`; by `xpol_summand_gauge` the
  summand depends on `s` only through its gauge class, so this is the uniform average over the flux
  configurations: random `π`-flux.
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

namespace FinGraph

variable {N M : ℕ} (G : FinGraph N M)

/-- The `σᶻ`-product gauge state `|s⟩_z` (edge `e` in the `σᶻ` eigenstate `sgn (s e)`) times the
matter state `ψ`, in the `σˣ` basis: amplitude `2^{-M/2} (-1)^{∑_e s_e x_e} ψ(n)`. -/
def initS (_G : FinGraph N M) (s : Cfg M) (ψ : Cfg N → ℂ) : Cfg N × Cfg M → ℂ :=
  fun p => (((Real.sqrt 2) ^ M : ℝ) : ℂ)⁻¹ * ((sgn (∑ e, s e * p.2 e) : ℝ) : ℂ) * ψ p.1

/-- The gauge transformation `θ` acting on link labels: `(s + δθ)_e = s_e + θ_{s(e)} + θ_{t(e)}`. -/
def gaugeAct (θ : Cfg N) (s : Cfg M) : Cfg M := fun e => s e + θ (G.src e) + θ (G.tgt e)

/-- **Fluxes are gauge classes** (paper: Lemma `lem:2dgauge`(i)): `s` and `s'` take the same values on
every cycle `z ∈ ker ∂` iff `s' = s + δθ` for some gauge transformation `θ`. -/
theorem flux_eq_iff_gauge (s s' : Cfg M) :
    (∀ z : Cfg M, G.bdry z = 0 → ∑ e, s e * z e = ∑ e, s' e * z e) ↔
      ∃ θ : Cfg N, s' = G.gaugeAct θ s := by
  sorry

/-- The `σᶻ`-product state with `s = 0` (all edges `|↑⟩_z`) is `initZ`. -/
theorem initS_zero (ψ : Cfg N → ℂ) : G.initS 0 ψ = G.initZ ψ := by
  sorry

/-- **`σᶻ`-product links: Anderson–Bernoulli models in a fixed flux sector** (paper: Proposition
`prop:2dlaws`(a)). -/
theorem sector_average_flux (J h : ℝ) (s : Cfg M) (n : Cfg N) (d : Cfg N → ℝ) (t : ℝ) :
    expect (cplx (diagonal d ⊗ₖ (1 : Matrix (Cfg M) (Cfg M) ℝ)))
        (evo (G.ham J h) t *ᵥ G.initS s (Pi.single n 1)) =
      (1 / ((G.sectors n).card : ℂ)) * ∑ q ∈ G.sectors n,
        expect (cplx (diagonal d))
          (evo (dGamma (G.hqs J h (fun i => sgn (q i)) (fun e => sgn (s e)))) t *ᵥ
            Pi.single n 1) := by
  sorry

/-- **`σˣ`-polarized links: random `π`-flux** (paper: Proposition `prop:2dlaws`(b)). The state
`|n⟩ ⊗ |+⋯+⟩` is the basis vector `(n, 0)` of the `σˣ` basis. -/
theorem sector_average_xpol (J h : ℝ) (n : Cfg N) (d : Cfg N → ℝ) (t : ℝ) :
    expect (cplx (diagonal d ⊗ₖ (1 : Matrix (Cfg M) (Cfg M) ℝ)))
        (evo (G.ham J h) t *ᵥ Pi.single (n, 0) 1) =
      (1 / 2 ^ M : ℂ) * ∑ s : Cfg M,
        expect (cplx (diagonal d))
          (evo (dGamma (G.hqs J h (fun i => sgn (n i)) (fun e => sgn (s e)))) t *ᵥ
            Pi.single n 1) := by
  sorry

/-- The free-fermion expectations depend on the signs only through their gauge class (hence only
through the fluxes, by `flux_eq_iff_gauge`). -/
theorem xpol_summand_gauge (J h : ℝ) (q n : Cfg N) (d : Cfg N → ℝ) (t : ℝ) (s : Cfg M)
    (θ : Cfg N) :
    expect (cplx (diagonal d))
        (evo (dGamma (G.hqs J h (fun i => sgn (q i)) (fun e => sgn (G.gaugeAct θ s e)))) t *ᵥ
          Pi.single n 1) =
      expect (cplx (diagonal d))
        (evo (dGamma (G.hqs J h (fun i => sgn (q i)) (fun e => sgn (s e)))) t *ᵥ
          Pi.single n 1) := by
  sorry

end FinGraph

end DFL
