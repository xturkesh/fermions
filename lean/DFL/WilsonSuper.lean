import DFL.Wilson
import DFL.QLMSuper

/-!
# Wilson links: superpositions of background-charge sectors (paper: Remark `rem:wilson`(a))

A finite superposition of sectors is described by a finite family of background charges
`Q : ι → (Fin L → ℤ)` (distinct sectors: `Q` injective; the results hold for any family) and the
components `Ψ₀ a : Cfg L → ℂ` of the state in the sectors `Q a` (a Fock vector `|n⟩` of sector `a`
stands for `|n⟩ ⊗ |E(n)⟩`). The Hamiltonian is block diagonal, so the component in sector `a`
evolves with the sector Hamiltonian `schwinger … (Q a)` (`wStateS`), the energy is the sum of the
sector energies (`energyDensityS`), and every quantity diagonal in the occupations and fields is the
sum over the sectors.

* `wilson_field_bound_super`: Proposition `prop:wilson`(b) holds for superpositions;
* `wilson_numberEntropy_le_super`: `S_N(Ψ(t); [1,ℓ]) ≤ H_{Ψ(t)}(E_ℓ) + H_{Ψ₀}(UQ_{[1,ℓ]})`, by the summed
  Gauss law and subadditivity; the second term is independent of time (`bgEntropyS_wStateS`);
* `lawEntropy_le_massey`: Lemma `lem:intentropy` for integer-valued laws on a finite space;
* `wilson_numberEntropy_avg_super`: the time-dependent part `S_N - H_{Ψ₀}(UQ)` is bounded on average
  over the cuts, uniformly in `t` and `L`;
* `bgEntropyS_bare_vacuum`, `bgEntropyS_le_log_three`: for bare-vacuum initial states
  (`ρ_j Ψ₀ = 0`) the second term is the entropy of the initial field on the cut link, at most `ln 3`
  if the initial fields take values in `{-1, 0, 1}`.
-/

open Matrix Finset Real

noncomputable section

namespace DFL

variable {L : ℕ} {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The Shannon entropy of the law of an integer-valued function `f` on a finite configuration
space in the state `v` (weights `|v c|²`). -/
def lawEntropy {κ : Type*} [Fintype κ] (f : κ → ℤ) (v : κ → ℂ) : ℝ :=
  ∑ y ∈ univ.image f, Real.negMulLog (∑ c ∈ univ.filter (fun c => f c = y), Complex.normSq (v c))

/-- The state at time `t` of a superposition of the sectors `Q a`: each component evolves with its
sector Hamiltonian. -/
def wStateS (w m g θ : ℝ) (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ₀ : ι → Cfg L → ℂ) (t : ℝ) :
    ι → Cfg L → ℂ :=
  fun a => wState w m g θ e₀ (Q a) (Ψ₀ a) t

/-- The energy density `⟨Ψ₀, HΨ₀⟩ / L = ∑_a ⟨Ψ₀^a, H_{Q_a} Ψ₀^a⟩ / L` of a superposition. -/
def energyDensityS (w m g θ : ℝ) (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ₀ : ι → Cfg L → ℂ) : ℝ :=
  (∑ a, (expect (cplx (schwinger w m g θ e₀ (Q a))) (Ψ₀ a)).re) / L

/-- `⟨Ψ, (E_j + ϑ)² Ψ⟩` for a superposition of sectors. -/
def fieldSqS (θ : ℝ) (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ : ι → Cfg L → ℂ) (j : Fin L) : ℝ :=
  ∑ a, fieldSq θ e₀ (Q a) (Ψ a) j

/-- The law of `N_R` in a superposition of sectors. -/
def numDistS (R : Finset (Fin L)) (Ψ : ι → Cfg L → ℂ) (k : ℕ) : ℝ :=
  ∑ a, numDist R (Ψ a) k

/-- The background charge `UQ_{[1,ℓ]} = ∑_{k ≤ ℓ} Q_k` (paper indexing; Lean sites `< ℓ`). -/
def bgSum (Q : Fin L → ℤ) (ℓ : ℕ) : ℤ :=
  ∑ k ∈ univ.filter (fun k : Fin L => k.val < ℓ), Q k

/-- `H_Ψ(E_j)`: the Shannon entropy of the law of the field on the link `j` in a superposition. -/
def fieldEntropyS (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ : ι → Cfg L → ℂ) (j : Fin L) : ℝ :=
  lawEntropy (fun c : ι × Cfg L => wField e₀ (Q c.1) c.2 j) (fun c => Ψ c.1 c.2)

/-- `H_Ψ(UQ_{[1,ℓ]})`: the Shannon entropy of the law of the background charge of `[1, ℓ]`. -/
def bgEntropyS (Q : ι → Fin L → ℤ) (Ψ : ι → Cfg L → ℂ) (ℓ : ℕ) : ℝ :=
  lawEntropy (fun c : ι × Cfg L => bgSum (Q c.1) ℓ) (fun c => Ψ c.1 c.2)

/-! ### Integer-valued laws -/

/-- **Lemma `lem:intentropy` for integer-valued laws**: `H(X) ≤ ½ ln(2πe(E(X + ϑ)² + 1/12))`. -/
theorem lawEntropy_le_massey {κ : Type*} [Fintype κ] (f : κ → ℤ) (v : κ → ℂ)
    (hv : ∑ c, Complex.normSq (v c) = 1) (θ : ℝ) :
    lawEntropy f v ≤ 1 / 2 * Real.log (2 * π * exp 1 *
      (∑ c, Complex.normSq (v c) * ((f c : ℝ) + θ) ^ 2 + 1 / 12)) := by
  sorry

/-- A law supported on a set `T` of values has entropy at most `ln |T|`. -/
theorem lawEntropy_le_log_card {κ : Type*} [Fintype κ] (f : κ → ℤ) (v : κ → ℂ)
    (hv : ∑ c, Complex.normSq (v c) = 1) (T : Finset ℤ) (hT : ∀ c, v c ≠ 0 → f c ∈ T) :
    lawEntropy f v ≤ Real.log T.card := by
  sorry

/-! ### Energy conservation for superpositions -/

/-- **Proposition `prop:wilson`(b) for superpositions of sectors** (paper: Remark `rem:wilson`(a)):
`ε₀ + w + |m| ≥ 0` and `(1/L) ∑_j ⟨(E_j + ϑ)²⟩_t ≤ (2/g²)(ε₀ + w + |m|)`. -/
theorem wilson_field_bound_super (w m g θ : ℝ) (hw : 0 ≤ w) (hg : 0 < g) (e₀ : ℤ)
    (Q : ι → Fin L → ℤ) (hL : 0 < L) (Ψ₀ : ι → Cfg L → ℂ)
    (hΨ₀ : ∑ a, ∑ n, Complex.normSq (Ψ₀ a n) = 1) (t : ℝ) :
    0 ≤ energyDensityS w m g θ e₀ Q Ψ₀ + w + |m| ∧
      (1 / (L : ℝ)) * ∑ j, fieldSqS θ e₀ Q (wStateS w m g θ e₀ Q Ψ₀ t) j ≤
        2 / g ^ 2 * (energyDensityS w m g θ e₀ Q Ψ₀ + w + |m|) := by
  sorry

/-! ### Number entropy of superpositions -/

/-- The law of the background charge does not change in time. -/
theorem bgEntropyS_wStateS (w m g θ : ℝ) (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ₀ : ι → Cfg L → ℂ)
    (t : ℝ) (ℓ : ℕ) : bgEntropyS Q (wStateS w m g θ e₀ Q Ψ₀ t) ℓ = bgEntropyS Q Ψ₀ ℓ := by
  sorry

/-- **Number entropy of superpositions** (paper: Remark `rem:wilson`(a)):
`S_N(Ψ(t); [1,ℓ]) ≤ H_{Ψ(t)}(E_ℓ) + H_{Ψ₀}(UQ_{[1,ℓ]})`. -/
theorem wilson_numberEntropy_le_super (w m g θ : ℝ) (e₀ : ℤ) (Q : ι → Fin L → ℤ)
    (Ψ₀ : ι → Cfg L → ℂ) (hΨ₀ : ∑ a, ∑ n, Complex.normSq (Ψ₀ a n) = 1) (t : ℝ) (ℓ : ℕ)
    (hℓ : 0 < ℓ) (hℓL : ℓ ≤ L) :
    shannon L (numDistS (segment ℓ) (wStateS w m g θ e₀ Q Ψ₀ t)) ≤
      fieldEntropyS e₀ Q (wStateS w m g θ e₀ Q Ψ₀ t) ⟨ℓ - 1, by omega⟩ + bgEntropyS Q Ψ₀ ℓ := by
  sorry

/-- **The time-dependent part of `S_N` is bounded on average over the cuts** (paper: Remark
`rem:wilson`(a)): `(1/L) ∑_ℓ (S_N(Ψ(t); [1,ℓ]) - H_{Ψ₀}(UQ_{[1,ℓ]})) ≤ ½ ln(2πe((2/g²)(ε₀ + w + |m|) + 1/12))`,
uniformly in `t` and `L`. -/
theorem wilson_numberEntropy_avg_super (w m g θ : ℝ) (hw : 0 ≤ w) (hg : 0 < g) (e₀ : ℤ)
    (Q : ι → Fin L → ℤ) (hL : 0 < L) (Ψ₀ : ι → Cfg L → ℂ)
    (hΨ₀ : ∑ a, ∑ n, Complex.normSq (Ψ₀ a n) = 1) (t : ℝ) :
    (1 / (L : ℝ)) * ∑ ℓ ∈ Icc 1 L,
        (shannon L (numDistS (segment ℓ) (wStateS w m g θ e₀ Q Ψ₀ t)) - bgEntropyS Q Ψ₀ ℓ) ≤
      1 / 2 * Real.log (2 * π * exp 1 *
        (2 / g ^ 2 * (energyDensityS w m g θ e₀ Q Ψ₀ + w + |m|) + 1 / 12)) := by
  sorry

/-! ### Bare-vacuum initial states -/

/-- For bare-vacuum initial states (`ρ_j Ψ₀ = 0`), `UQ_{[1,ℓ]} = E_ℓ - e₀` on `Ψ₀`, so the background
term is the entropy of the initial field on the cut link. -/
theorem bgEntropyS_bare_vacuum (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ₀ : ι → Cfg L → ℂ)
    (hvac : ∀ a n, Ψ₀ a n ≠ 0 → ∀ j, stagRho n j = 0) (ℓ : ℕ) (hℓ : 0 < ℓ) (hℓL : ℓ ≤ L) :
    bgEntropyS Q Ψ₀ ℓ = fieldEntropyS e₀ Q Ψ₀ ⟨ℓ - 1, by omega⟩ := by
  sorry

/-- For bare-vacuum initial states whose field on the cut link takes values in `{-1, 0, 1}`, the
background term is at most `ln 3`. -/
theorem bgEntropyS_le_log_three (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ₀ : ι → Cfg L → ℂ)
    (hΨ₀ : ∑ a, ∑ n, Complex.normSq (Ψ₀ a n) = 1)
    (hvac : ∀ a n, Ψ₀ a n ≠ 0 → ∀ j, stagRho n j = 0) (ℓ : ℕ) (hℓ : 0 < ℓ) (hℓL : ℓ ≤ L)
    (h3 : ∀ a n, Ψ₀ a n ≠ 0 → wField e₀ (Q a) n ⟨ℓ - 1, by omega⟩ ∈ ({-1, 0, 1} : Finset ℤ)) :
    bgEntropyS Q Ψ₀ ℓ ≤ Real.log 3 := by
  sorry

end DFL
