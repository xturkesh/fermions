import DFL.SlaterGeneral
import DFL.AnnealedInterval

/-!
# Entanglement of Slater determinants for arbitrary regions (paper: Lemma `lem:slater`(a),(b))

For a region `R` (e.g. a bulk interval) the paper identifies `𝓕_Λ ≅ 𝓕_R ⊗ 𝓕_B`, `B = Λ \ R`, by
`c(a) ↦ c(a) ⊗ 1`, `c(b) ↦ (-1)^{N_R} ⊗ c(b)`. In the Jordan–Wigner occupation basis this maps `|n⟩` to
`ε(α, β) |α⟩ ⊗ |β⟩`, where `α = n|_R`, `β = n|_B` and `ε(α, β) = (-1)^{#{(a, b) : a ∈ R, b ∉ R, b < a, α_a = β_b = 1}}`
(`reorderSign`). The reduced density matrix of the modes in `R` (`rdmF`) is therefore
`ρ_R(α, α') = ∑_β ε(α, β) ε(α', β) ψ(α ⊔ β) conj ψ(α' ⊔ β)`.

* `vonNeumann_rdmF_slater`: for `ψ = Γ(e^{-itA})|S⟩`, `S(ρ_R) = ∑_i h(ν_i)`, `ν_i` the eigenvalues of `C_R`;
* `vonNeumann_rdmF_le`: `S(ρ_R) ≤ 2 ln 2 ∑_{x ∈ R, y ∉ R} |K_{xy}|`;
* `vonNeumann_rdmF_segment`: for a left segment this is the reduced state `rdm` of `DFL.AreaLaw`;
* `lgt_areaLaw_bulk`: for an interval `B = [a, b]`, `E_q sup_t S(ρ^q_B(t)) ≤ 4 ln 2 Θ(F)` under (EC)
  (finite-volume form of Theorem `thm:infvol`(d)).
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Splitting a configuration into its restrictions to `R` and to the complement of `R`. -/
def splitR (R : Finset (Fin N)) : Cfg N ≃ (R → ZMod 2) × ({x : Fin N // x ∉ R} → ZMod 2) :=
  Equiv.piEquivPiSubtypeProd (fun x : Fin N => x ∈ R) (fun _ => ZMod 2)

/-- The reordering sign `ε(α, β) = (-1)^{#{(a, b) : a ∈ R, b ∉ R, b < a, α_a = β_b = 1}}`. -/
def reorderSign (R : Finset (Fin N)) (α : R → ZMod 2) (β : {x : Fin N // x ∉ R} → ZMod 2) : ℝ :=
  (-1 : ℝ) ^ ((univ : Finset (R × {x : Fin N // x ∉ R})).filter
    (fun p => (p.2 : Fin N) < (p.1 : Fin N) ∧ α p.1 = 1 ∧ β p.2 = 1)).card

/-- The reduced density matrix of the modes in `R`:
`ρ_R(α, α') = ∑_β ε(α, β) ε(α', β) ψ(α ⊔ β) conj ψ(α' ⊔ β)`. -/
def rdmF (R : Finset (Fin N)) (v : Cfg N → ℂ) : Matrix (R → ZMod 2) (R → ZMod 2) ℂ :=
  fun α α' => ∑ β : {x : Fin N // x ∉ R} → ZMod 2,
    ((reorderSign R α β * reorderSign R α' β : ℝ) : ℂ) *
      v ((splitR R).symm (α, β)) * star (v ((splitR R).symm (α', β)))

/-- `ρ_R` is Hermitian. -/
theorem rdmF_isHermitian (R : Finset (Fin N)) (v : Cfg N → ℂ) : (rdmF R v).IsHermitian := by
  ext α α'
  simp only [Matrix.conjTranspose_apply, rdmF, star_sum]
  refine Finset.sum_congr rfl fun β _ => ?_
  rw [star_mul', star_mul', star_star, Complex.star_def, Complex.conj_ofReal]
  push_cast
  ring

/-- **Entanglement entropy of a Slater determinant for any region** (paper: Lemma `lem:slater`(a)):
`S(ρ_R) = ∑_i h(ν_i)`. -/
theorem vonNeumann_rdmF_slater (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (R : Finset (Fin N)) :
    vonNeumann (rdmF R (evo (dGamma A) t *ᵥ fockVec S)) (rdmF_isHermitian R _) =
      ∑ i, Real.binEntropy ((corrR_isHermitian A S t R).eigenvalues i) := by
  sorry

/-- **Area law for Slater determinants and any region** (paper: Lemma `lem:slater`(b), Lemma `lem:binent`):
`S(ρ_R) ≤ 2 ln 2 ∑_{x ∈ R, y ∉ R} |K_{xy}|`. -/
theorem vonNeumann_rdmF_le (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (R : Finset (Fin N)) :
    vonNeumann (rdmF R (evo (dGamma A) t *ᵥ fockVec S)) (rdmF_isHermitian R _) ≤
      2 * Real.log 2 * ∑ x ∈ R, ∑ y ∈ univ \ R, ‖corrMat A S t x y‖ := by
  sorry

/-- For a left segment, `ρ_R` is the reduced state `rdm` of `DFL.AreaLaw` (no reordering signs). -/
theorem vonNeumann_rdmF_segment (ℓ : ℕ) (hℓ : ℓ ≤ N) (v : Cfg N → ℂ) :
    vonNeumann (rdmF (segment ℓ) v) (rdmF_isHermitian _ v) =
      vonNeumann (rdm ℓ hℓ v) (rdm_isHermitian ℓ hℓ v) := by
  sorry

/-- **Area law for bulk intervals** (paper: finite-volume form of Theorem `thm:infvol`(d)):
`E_q sup_t S(ρ^q_B(t)) ≤ 2 ln 2 · 2Θ(F) = 4 ln 2 Θ(F)` for every interval `B = [a, b]`. -/
theorem lgt_areaLaw_bulk (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (a b : ℕ) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ,
      vonNeumann (rdmF (ivl a b) (evo (dGamma (hq J h q)) t *ᵥ fockVec S)) (rdmF_isHermitian _ _)) ≤
      4 * Real.log 2 * thetaF F := by
  sorry

end DFL
