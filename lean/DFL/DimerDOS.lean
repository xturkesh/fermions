import DFL.ECInfVolQ
import DFL.DimerImbalanceLim
import DFL.DimerLGT

/-!
# The density of states of the dimer model (paper: Remark `rem:dimerdos`)

On `ℓ²(ℤ)` let `H_η = Δ + λ diag(q(η))` with the dimer charges `q_{2k-1} = q_{2k} = η_k` (`dimerQZ`,
`η` i.i.d. uniform under `muZ`); conjugation by `Γ = diag((-1)^x)` turns it into `-Δ + λ diag(q)` without
changing the spectral measures at the vectors `δ_x`. The density-of-states functional is
`ν(f) = ½ E(⟨δ_1, f(H) δ_1⟩ + ⟨δ_2, f(H) δ_2⟩)` (`dosFun`, with the continuous functional calculus `cfc`
of `z ↦ f(Re z)`; the spectrum of `H` is real).

* `tendsto_dosFin` (`rem:dimerdos`(ii)): for continuous `f`, `N⁻¹ E_dim Tr f(h_q/J) → ν(f)` along even `N`
  (finite propagation for polynomials, stationarity of the dimer law under the shift by one dimer, and
  Weierstrass approximation on the spectral interval `[-2-|λ|, 2+|λ|]`).
* `dosFun_nonneg`, `dosFun_one`, `dosFun_eq_zero_of_support`, `dosFun_le_of_near`: `ν` is positive and
  normalized, it lives on `[-2-|λ|, 2+|λ|]`, and it gives mass at most `2/L` to every interval
  `[E - ε, E + ε]` with `ε L (|λ| + |E| + 1)^L < 1` (from `card_eig_near_le`, the deterministic finite-volume
  count of `lem:fewapprox`).
* `exists_dosMeasure`, `dosMeasure`, `isProbabilityMeasure_dosMeasure`, `dosMeasure_compl_Icc`,
  `dosMeasure_singleton` (`rem:dimerdos`(i)): by the Riesz–Markov–Kakutani theorem `ν` is a probability measure
  on `[-2-|λ|, 2+|λ|]` without atoms.
-/

open MeasureTheory Filter Topology Matrix
open scoped InnerProductSpace

noncomputable section

namespace DFL

/-- The dimer charges on `ℤ` with the dimers `{2k-1, 2k}`: `q_y = η_{⌊(y+1)/2⌋}`. -/
def dimerQZ (η : CfgZ) : CfgZ := fun y => η ((y + 1) / 2)

/-- The dimer Hamiltonian `H_η = Δ + λ diag(q(η))` on `ℓ²(ℤ)`. -/
def dimerHZ (lam : ℝ) (η : CfgZ) : L2Z →L[ℂ] L2Z := ECInfVolQAux.AZop lam (dimerQZ η)

/-- The density-of-states functional `ν(f) = ½ E(⟨δ_1, f(H) δ_1⟩ + ⟨δ_2, f(H) δ_2⟩)`, with `f(H)` the
continuous functional calculus of the self-adjoint `H` applied to `z ↦ f(Re z)`. -/
def dosFun (lam : ℝ) (f : ℝ → ℝ) : ℝ :=
  (1 / 2) * ∫ η, (⟪deltaZ 1, cfc (fun z : ℂ => ((f z.re : ℝ) : ℂ)) (dimerHZ lam η) (deltaZ 1)⟫_ℂ +
    ⟪deltaZ 2, cfc (fun z : ℂ => ((f z.re : ℝ) : ℂ)) (dimerHZ lam η) (deltaZ 2)⟫_ℂ).re ∂muZ

/-- The finite-volume eigenvalue average `N⁻¹ E_dim Tr f(h_q/J)`, `h_q/J = -Δ_N + λ diag(q)`, `N = 2M`. -/
def dosFin (lam : ℝ) (M : ℕ) (f : ℝ → ℝ) : ℝ :=
  (1 / (2 * M : ℝ)) * dimerAvg (fun q : Cfg (2 * M) =>
    ∑ α, f ((hq_isHermitian 1 (lam / 2) q).eigenvalues α))

/-- **Convergence of the finite-volume density of states** (paper: `rem:dimerdos`(ii)). -/
theorem tendsto_dosFin (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f) :
    Tendsto (fun M : ℕ => dosFin lam M f) atTop (𝓝 (dosFun lam f)) := by
  sorry

/-- `ν` is positive. -/
theorem dosFun_nonneg (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f) (h0 : ∀ x, 0 ≤ f x) :
    0 ≤ dosFun lam f := by
  sorry

/-- `ν` is normalized. -/
theorem dosFun_one (lam : ℝ) : dosFun lam (fun _ => 1) = 1 := by
  sorry

/-- `ν` lives on the spectral interval `[-2-|λ|, 2+|λ|]`. -/
theorem dosFun_eq_zero_of_support (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f)
    (hsupp : ∀ x ∈ Set.Icc (-2 - |lam|) (2 + |lam|), f x = 0) : dosFun lam f = 0 := by
  sorry

/-- **No concentration near any energy**: `ν(f) ≤ 2/L` for `0 ≤ f ≤ 1` supported in `[E - ε, E + ε]` with
`ε L (|λ| + |E| + 1)^L < 1`. -/
theorem dosFun_le_of_near (lam : ℝ) (f : ℝ → ℝ) (hf : Continuous f)
    (h01 : ∀ x, 0 ≤ f x ∧ f x ≤ 1) (E ε : ℝ) (L : ℕ) (hL : 0 < L)
    (hε : ε * L * (|lam| + |E| + 1) ^ L < 1) (hsupp : ∀ x, f x ≠ 0 → |x - E| ≤ ε) :
    dosFun lam f ≤ 2 / L := by
  sorry

/-- **The density of states is a probability measure** (Riesz–Markov–Kakutani). -/
theorem exists_dosMeasure (lam : ℝ) :
    ∃ ν : Measure ℝ, IsProbabilityMeasure ν ∧ ∀ f : CompactlySupportedContinuousMap ℝ ℝ, ∫ x, f x ∂ν = dosFun lam f := by
  sorry

/-- The density of states `ν` of the dimer model. -/
def dosMeasure (lam : ℝ) : Measure ℝ := Classical.choose (exists_dosMeasure lam)

/-- `ν` is a probability measure (paper: `rem:dimerdos`(i)). -/
theorem isProbabilityMeasure_dosMeasure (lam : ℝ) : IsProbabilityMeasure (dosMeasure lam) :=
  (Classical.choose_spec (exists_dosMeasure lam)).1

/-- `∫ f dν = ν(f)` for compactly supported continuous `f`. -/
theorem integral_dosMeasure (lam : ℝ) (f : CompactlySupportedContinuousMap ℝ ℝ) :
    ∫ x, f x ∂(dosMeasure lam) = dosFun lam f :=
  (Classical.choose_spec (exists_dosMeasure lam)).2 f

/-- `ν` is supported in `[-2-|λ|, 2+|λ|]` (paper: `rem:dimerdos`(i)). -/
theorem dosMeasure_compl_Icc (lam : ℝ) :
    dosMeasure lam (Set.Icc (-2 - |lam|) (2 + |lam|))ᶜ = 0 := by
  sorry

/-- **`ν` has no atoms** (paper: `rem:dimerdos`(i)). -/
theorem dosMeasure_singleton (lam E : ℝ) : dosMeasure lam {E} = 0 := by
  sorry

end DFL
