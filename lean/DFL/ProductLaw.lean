import DFL.Sector

/-!
# Sector laws of product initial states (paper: Lemma `lem:productlaw`)

For a product state `χ = ⊗_j χ_j` of the links (given by its `σˣ` amplitudes `χ_j(x_j)`) and any
matter state `ψ`, the dual state is `∑_q ψ^{(q)} ⊗ |q⟩` with `ψ^{(q)}(n) = χ(x(n,q)) ψ(n)`.
Consequences:
* (b) Fock matter `|n⟩`: the sector law is `μ(q) = |χ(x(n,q))|²`, the image of the product law of
  the link labels under the bijection `x ↦ q`;
* (c) `σᶻ`-product gauge states `|s⟩_z`: the charges are i.i.d. Bernoulli(1/2), with matter state
  `W_D ψ` (a sign `(-1)^{N_{≤j}}` for every link `j` in the state `|↓⟩_z`);
* (d) `σˣ`-product gauge states with Fock matter: a single, definite-charge sector.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- A product state of the links, given by the `σˣ` amplitudes `χ j (x_j)` of each link. -/
def prodLink (χ : Fin N → ZMod 2 → ℂ) (x : Cfg N) : ℂ := ∏ j, χ j (x j)

/-- The product initial state `χ ⊗ ψ` (links in the product state `χ`, matter in `ψ`). -/
def prodInit (χ : Fin N → ZMod 2 → ℂ) (ψ : Cfg N → ℂ) : Cfg N × Cfg N → ℂ :=
  fun p => ψ p.1 * prodLink χ p.2

/-- (a) The dual form of a product initial state. -/
theorem prodInit_dual (χ : Fin N → ZMod 2 → ℂ) (ψ : Cfg N → ℂ) (n q : Cfg N) :
    prodInit χ ψ (dualEquiv (n, q)) = prodLink χ (linkOf n q) * ψ n := by
  sorry

/-- (a) Sector decomposition of the dynamics for a product initial state: every gauge-invariant
expectation is a sum over sectors of free-fermion expectations with initial matter states
`ψ^{(q)} = χ(x(·, q)) ψ`. -/
theorem sector_sum_prod (J h : ℝ) (χ : Fin N → ZMod 2 → ℂ) (ψ : Cfg N → ℂ)
    (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (H J h) t *ᵥ prodInit χ ψ) =
      ∑ q : Cfg N, expect (sectorBlock O q)
        (evo (Hsec J h q) t *ᵥ fun n => prodLink χ (linkOf n q) * ψ n) := by
  sorry

/-- (b) Fock matter: the sector law is `μ(q) = |χ(x(n,q))|²`. -/
theorem sector_average_prod_fock (J h : ℝ) (χ : Fin N → ZMod 2 → ℂ) (n : Cfg N)
    (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (H J h) t *ᵥ prodInit χ (Pi.single n 1)) =
      ∑ q : Cfg N, (Complex.normSq (prodLink χ (linkOf n q)) : ℂ) *
        expect (sectorBlock O q) (evo (Hsec J h q) t *ᵥ Pi.single n 1) := by
  sorry

/-- (b) The sector law is the image of the product law of the link labels under `x ↦ q`. -/
theorem sum_sectorLaw_prod (χ : Fin N → ZMod 2 → ℂ) (n : Cfg N) (f : Cfg N → ℂ) :
    ∑ q : Cfg N, (Complex.normSq (prodLink χ (linkOf n q)) : ℂ) * f q =
      ∑ x : Cfg N, (Complex.normSq (prodLink χ x) : ℂ) * f (chargeOf n x) := by
  sorry

/-- `σˣ` amplitudes of the `σᶻ` eigenstates: `⟨x|↑⟩_z = 2^{-1/2}`, `⟨x|↓⟩_z = 2^{-1/2} (-1)^x`
(label `s = 0` for `|↑⟩_z`, `s = 1` for `|↓⟩_z`). -/
def zAmp (s x : ZMod 2) : ℂ := ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ * (if s = 1 ∧ x = 1 then -1 else 1)

/-- The matter sign `W_D = ∏_{j : s_j = 1} (-1)^{N_{≤ j}}`. -/
def signWD (s : Cfg N) (n : Cfg N) : ℂ :=
  ∏ j, if s j = 1 then ((sgn (∑ k ∈ univ.filter (fun k => k ≤ j), n k) : ℝ) : ℂ) else 1

/-- (c) `σᶻ`-product gauge states: the charges are i.i.d. Bernoulli(1/2) and every sector carries
the matter state `W_D ψ`. -/
theorem sector_average_zprod (J h : ℝ) (s : Cfg N) (ψ : Cfg N → ℂ)
    (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (H J h) t *ᵥ prodInit (fun j => zAmp (s j)) ψ) =
      (1 / 2 ^ N : ℂ) * ∑ q : Cfg N,
        expect (sectorBlock O q) (evo (Hsec J h q) t *ᵥ fun n => signWD s n * ψ n) := by
  sorry

/-- (d) `σˣ`-product gauge states `|ξ⟩_σ` with Fock matter `|n⟩` are definite-charge states:
the dynamics stays in the single sector `q = chargeOf n ξ`. -/
theorem expect_xprod_fock (J h : ℝ) (ξ n : Cfg N)
    (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (H J h) t *ᵥ prodInit (fun j x => if x = ξ j then 1 else 0) (Pi.single n 1)) =
      expect (sectorBlock O (chargeOf n ξ))
        (evo (Hsec J h (chargeOf n ξ)) t *ᵥ Pi.single n 1) := by
  sorry

end DFL
