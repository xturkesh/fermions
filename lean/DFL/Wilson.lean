import DFL.FreeDyn
import DFL.Entropy
import DFL.AreaLaw

/-!
# Wilson links: cut-averaged number entropy from energy conservation
(paper: Proposition `prop:wilson`)

The lattice Schwinger model with Wilson links (paper: Eq. `eq:HU1`),
`H = -w ∑_j (f†_j U_j f_{j+1} + h.c.) + m ∑_j (-1)^j n_j + (g²/2) ∑_j (E_j + ϑ)²`,
in a fixed background-charge sector `Q`. The Gauss law `G_j = E_j - E_{j-1} - ρ_j = Q_j` (with
`E_0 = e₀` fixed) determines the electric fields from the occupations,
`E_j = e₀ + ∑_{k ≤ j} (ρ_k + Q_k)`, and the hopping `f†_j U_j f_{j+1}` raises exactly the field it
crosses (`wField_hop`). The sector is therefore the space of functions of the occupations, and `H`
acts on it as the `2^L × 2^L` real symmetric matrix `schwinger` (a Fock basis vector `|n⟩` stands for
`|n⟩ ⊗ |E(n)⟩`).

Indexing: Lean sites are `0, …, L-1`; Lean site `j` is paper site `j + 1`. So the paper's
`(-1)^j` is `(-1)^{j+1}` here, the staggered charge is `ρ_j = n_j - [j even]`, the paper's
interval `[1, ℓ]` is `segment ℓ`, and the paper's field `E_ℓ` (link `b_ℓ`) is `wField … ⟨ℓ - 1, _⟩`.
-/

open Matrix Finset Real

noncomputable section

namespace DFL

variable {L : ℕ}

/-- Staggered charge `ρ_j = n_j - [j even]` (0-indexed; paper: `n_j - ½(1 - (-1)^j)`, 1-indexed). -/
def stagRho (n : Cfg L) (j : Fin L) : ℤ :=
  (if n j = 1 then 1 else 0) - (if j.val % 2 = 0 then 1 else 0)

/-- The electric field on the link `b_j` in the sector `Q`: the Gauss law integrated from the left
edge, `E_j = e₀ + ∑_{k ≤ j} (ρ_k + Q_k)`. -/
def wField (e₀ : ℤ) (Q : Fin L → ℤ) (n : Cfg L) (j : Fin L) : ℤ :=
  e₀ + ∑ k ∈ univ.filter (fun k => k ≤ j), (stagRho n k + Q k)

/-- The lattice Schwinger Hamiltonian in the background-charge sector `Q`. -/
def schwinger (w m g θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) : Matrix (Cfg L) (Cfg L) ℝ :=
  (-w) • ∑ j, hop j +
    diagonal (fun n => m * ∑ j : Fin L, (-1 : ℝ) ^ (j.val + 1) * occ (n j) +
      g ^ 2 / 2 * ∑ j, ((wField e₀ Q n j : ℝ) + θ) ^ 2)

/-- The state at time `t`. -/
def wState (w m g θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) (ψ₀ : Cfg L → ℂ) (t : ℝ) : Cfg L → ℂ :=
  evo (schwinger w m g θ e₀ Q) t *ᵥ ψ₀

/-- The energy density `ε₀ = ⟨ψ₀, H ψ₀⟩ / L`. -/
def energyDensity (w m g θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) (ψ₀ : Cfg L → ℂ) : ℝ :=
  (expect (cplx (schwinger w m g θ e₀ Q)) ψ₀).re / L

/-- `⟨ψ, (E_j + ϑ)² ψ⟩`. -/
def fieldSq (θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) (ψ : Cfg L → ℂ) (j : Fin L) : ℝ :=
  ∑ n, Complex.normSq (ψ n) * ((wField e₀ Q n j : ℝ) + θ) ^ 2

/-- The constant in the summed Gauss law: `E_ℓ - N_{[1,ℓ]} = e₀ - #{odd paper sites in [1,ℓ]} + ∑_{k ≤ ℓ} Q_k`. -/
def wConst (e₀ : ℤ) (Q : Fin L → ℤ) (ℓ : ℕ) : ℝ :=
  (e₀ : ℝ) - ((univ.filter (fun k : Fin L => k.val < ℓ ∧ k.val % 2 = 0)).card : ℝ) +
    ∑ k ∈ univ.filter (fun k : Fin L => k.val < ℓ), (Q k : ℝ)

/-- The hopping across the bond `(j, j+1)` changes no field other than `E_j` (so `schwinger` is the
restriction of `H` to the sector). -/
theorem wField_hop (e₀ : ℤ) (Q : Fin L → ℤ) (j k : Fin L) (n' n : Cfg L) (hjk : k ≠ j)
    (hhop : hop j n' n ≠ 0) : wField e₀ Q n' k = wField e₀ Q n k := by
  sorry

/-- **Summed Gauss law** (paper: Proposition `prop:wilson`(a), via Eq. `eq:summedgauss`): in the
sector, `E_ℓ` is an injective affine function of `N_{[1,ℓ]}`, so their laws have the same Shannon
entropy. -/
theorem wField_eq_count (e₀ : ℤ) (Q : Fin L → ℤ) (n : Cfg L) (ℓ : ℕ) (hℓ : 0 < ℓ) (hℓL : ℓ ≤ L) :
    (wField e₀ Q n ⟨ℓ - 1, by omega⟩ : ℝ) = count (segment ℓ) n + wConst e₀ Q ℓ := by
  sorry

/-- **Energy bound on the fields** (paper: Proposition `prop:wilson`(b)): `ε₀ + w + |m| ≥ 0` and
`(1/L) ∑_ℓ ⟨(E_ℓ + ϑ)²⟩_t ≤ (2/g²)(ε₀ + w + |m|)` for all `t`. -/
theorem wilson_field_bound (w m g θ : ℝ) (hw : 0 ≤ w) (hg : 0 < g) (e₀ : ℤ) (Q : Fin L → ℤ)
    (hL : 0 < L) (ψ₀ : Cfg L → ℂ) (hψ₀ : ∑ n, Complex.normSq (ψ₀ n) = 1) (t : ℝ) :
    0 ≤ energyDensity w m g θ e₀ Q ψ₀ + w + |m| ∧
      (1 / (L : ℝ)) * ∑ j, fieldSq θ e₀ Q (wState w m g θ e₀ Q ψ₀ t) j ≤
        2 / g ^ 2 * (energyDensity w m g θ e₀ Q ψ₀ + w + |m|) := by
  sorry

/-- **Cut-averaged number entropy** (paper: Proposition `prop:wilson`(c)). -/
theorem wilson_numberEntropy_avg (w m g θ : ℝ) (hw : 0 ≤ w) (hg : 0 < g) (e₀ : ℤ)
    (Q : Fin L → ℤ) (hL : 0 < L) (ψ₀ : Cfg L → ℂ) (hψ₀ : ∑ n, Complex.normSq (ψ₀ n) = 1)
    (t : ℝ) :
    (1 / (L : ℝ)) * ∑ ℓ ∈ Icc 1 L, shannon L (numDist (segment ℓ) (wState w m g θ e₀ Q ψ₀ t)) ≤
      1 / 2 * Real.log (2 * π * exp 1 *
        (2 / g ^ 2 * (energyDensity w m g θ e₀ Q ψ₀ + w + |m|) + 1 / 12)) := by
  sorry

/-- **Most cuts** (paper: Proposition `prop:wilson`(d)): for `δ ∈ (0, 1]`, the bound with `δ g²` in
place of `g²` fails for at most `δ L` of the cuts. -/
theorem wilson_numberEntropy_most (w m g θ : ℝ) (hw : 0 ≤ w) (hg : 0 < g) (e₀ : ℤ)
    (Q : Fin L → ℤ) (hL : 0 < L) (ψ₀ : Cfg L → ℂ) (hψ₀ : ∑ n, Complex.normSq (ψ₀ n) = 1)
    (t δ : ℝ) (hδ0 : 0 < δ) (hδ1 : δ ≤ 1) :
    (((Icc 1 L).filter (fun ℓ =>
        1 / 2 * Real.log (2 * π * exp 1 *
          (2 * (energyDensity w m g θ e₀ Q ψ₀ + w + |m|) / (δ * g ^ 2) + 1 / 12)) <
        shannon L (numDist (segment ℓ) (wState w m g θ e₀ Q ψ₀ t)))).card : ℝ) ≤ δ * L := by
  sorry

/-- **Every cut** (paper: Proposition `prop:wilson`(e)). -/
theorem wilson_numberEntropy_each (w m g θ : ℝ) (hw : 0 ≤ w) (hg : 0 < g) (e₀ : ℤ)
    (Q : Fin L → ℤ) (hL : 0 < L) (ψ₀ : Cfg L → ℂ) (hψ₀ : ∑ n, Complex.normSq (ψ₀ n) = 1)
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ∈ Icc 1 L) :
    shannon L (numDist (segment ℓ) (wState w m g θ e₀ Q ψ₀ t)) ≤
      1 / 2 * Real.log (2 * π * exp 1 *
        (2 * L * (energyDensity w m g θ e₀ Q ψ₀ + w + |m|) / g ^ 2 + 1 / 12)) := by
  sorry

end DFL
