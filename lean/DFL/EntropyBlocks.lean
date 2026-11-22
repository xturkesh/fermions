import DFL.Entropies

/-!
# Block decomposition of the entropy (paper: Proposition `prop:entropies`(c), the identity)

For a Hermitian matrix that is block diagonal with respect to a labelling `lab : ι → κ`, the von
Neumann entropy is the sum of the entropies of the diagonal blocks (`vonNeumann_block`), and for a
positive semidefinite block `B` of trace `p`, `S(B) = p S(B/p) + η(p)` (`vonNeumann_smul_trace`).
Hence `S(ρ) = H(p) + ∑_k p_k S(ρ^{(k)})`, `ρ^{(k)} = B_k / p_k` (`vonNeumann_block_decomp`). Applied to the
reduced states of the quench, block diagonal in the particle number `N_A` of `A` with weights the
number distribution `p_t`, this is the identity `S_X(t) = S_N(t) + ∑_k p_t(k) S(ρ_X^{(k)}(t))` for
`X ∈ {full, f}` (paper: Proposition `prop:entropies`(c)).
-/

open Matrix Finset
open scoped ComplexOrder

noncomputable section

namespace DFL

variable {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ]

/-- The diagonal block of `ρ` on the labels `{a | lab a = k}`. -/
def blockOf (ρ : Matrix ι ι ℂ) (lab : ι → κ) (k : κ) :
    Matrix {a // lab a = k} {a // lab a = k} ℂ :=
  ρ.submatrix (fun a => a.1) (fun a => a.1)

/-- The normalized block `B_k / Tr B_k` (the zero matrix if `Tr B_k = 0`). -/
def normBlock (ρ : Matrix ι ι ℂ) (lab : ι → κ) (k : κ) :
    Matrix {a // lab a = k} {a // lab a = k} ℂ :=
  ((blockOf ρ lab k).trace)⁻¹ • blockOf ρ lab k

/-- Blocks of a Hermitian matrix are Hermitian. -/
theorem blockOf_isHermitian (ρ : Matrix ι ι ℂ) (hρ : ρ.IsHermitian) (lab : ι → κ) (k : κ) :
    (blockOf ρ lab k).IsHermitian := by
  sorry

/-- Normalized blocks of a Hermitian matrix are Hermitian. -/
theorem normBlock_isHermitian (ρ : Matrix ι ι ℂ) (hρ : ρ.IsHermitian) (lab : ι → κ) (k : κ) :
    (normBlock ρ lab k).IsHermitian := by
  sorry

/-- **Additivity over blocks**: the entropy of a block-diagonal Hermitian matrix is the sum of the
entropies of its blocks. -/
theorem vonNeumann_block (ρ : Matrix ι ι ℂ) (hρ : ρ.IsHermitian) (lab : ι → κ)
    (hblock : ∀ a b, lab a ≠ lab b → ρ a b = 0) :
    vonNeumann ρ hρ = ∑ k, vonNeumann (blockOf ρ lab k) (blockOf_isHermitian ρ hρ lab k) := by
  sorry

/-- **Scaling**: for a positive semidefinite `B` of trace `p`, `S(B) = p S(B/p) + η(p)`. -/
theorem vonNeumann_smul_trace {μ : Type*} [Fintype μ] [DecidableEq μ] (B : Matrix μ μ ℂ)
    (hB : B.PosSemidef) (hB' : ((B.trace)⁻¹ • B).IsHermitian) :
    vonNeumann B hB.isHermitian =
      (B.trace).re * vonNeumann ((B.trace)⁻¹ • B) hB' + Real.negMulLog (B.trace).re := by
  sorry

/-- **Block decomposition**: for a positive semidefinite block-diagonal `ρ`,
`S(ρ) = ∑_k η(p_k) + ∑_k p_k S(ρ^{(k)})` with `p_k = Tr B_k`. -/
theorem vonNeumann_block_decomp (ρ : Matrix ι ι ℂ) (hρ : ρ.PosSemidef) (lab : ι → κ)
    (hblock : ∀ a b, lab a ≠ lab b → ρ a b = 0) :
    vonNeumann ρ hρ.isHermitian =
      ∑ k, Real.negMulLog ((blockOf ρ lab k).trace).re +
        ∑ k, ((blockOf ρ lab k).trace).re *
          vonNeumann (normBlock ρ lab k) (normBlock_isHermitian ρ hρ.isHermitian lab k) := by
  sorry

variable {N : ℕ}

/-- The matter particle number of `A` as a label in `{0, …, N}` of the configurations of `A`. -/
def cntA (ℓ : ℕ) (hℓ : ℓ ≤ N) (a : Cfg ℓ × Cfg ℓ) : Fin (N + 1) :=
  ⟨(univ.filter (fun i => a.1 i = 1)).card, by
    have := Finset.card_le_univ (univ.filter (fun i => a.1 i = 1))
    simp only [Fintype.card_fin] at this
    omega⟩

/-- The particle number of the bare modes of `A` as a label in `{0, …, N}`. -/
def cntBare (ℓ : ℕ) (hℓ : ℓ ≤ N) (a : Cfg ℓ) : Fin (N + 1) :=
  ⟨(univ.filter (fun i => a i = 1)).card, by
    have := Finset.card_le_univ (univ.filter (fun i => a i = 1))
    simp only [Fintype.card_fin] at this
    omega⟩

/-- **Entropy taxonomy, the identity for `S_full`** (paper: Proposition `prop:entropies`(c)):
`S_full(t) = S_N(t) + ∑_k p_t(k) S(ρ_full^{(k)}(t))` for Fock initial states. -/
theorem fullEntropy_eq_numberEntropy_add (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (rdmLGT ℓ hℓ (lgtState J h (Pi.single n 1) t)) (rdmLGT_isHermitian ℓ hℓ _) =
      shannon N (lgtNumDist J h (segment ℓ) (Pi.single n 1) t) +
        ∑ k : Fin (N + 1), lgtNumDist J h (segment ℓ) (Pi.single n 1) t k *
          vonNeumann (normBlock (rdmLGT ℓ hℓ (lgtState J h (Pi.single n 1) t)) (cntA ℓ hℓ) k)
            (normBlock_isHermitian _ (rdmLGT_isHermitian ℓ hℓ _) _ k) := by
  sorry

/-- **Entropy taxonomy, the identity for `S_f`** (paper: Proposition `prop:entropies`(c)):
`S_f(t) = S_N(t) + ∑_k p_t(k) S(ρ_f^{(k)}(t))` for Fock initial states. -/
theorem bareEntropy_eq_numberEntropy_add (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (rdmBare ℓ hℓ (lgtState J h (Pi.single n 1) t)) (rdmBare_isHermitian ℓ hℓ _) =
      shannon N (lgtNumDist J h (segment ℓ) (Pi.single n 1) t) +
        ∑ k : Fin (N + 1), lgtNumDist J h (segment ℓ) (Pi.single n 1) t k *
          vonNeumann (normBlock (rdmBare ℓ hℓ (lgtState J h (Pi.single n 1) t)) (cntBare ℓ hℓ) k)
            (normBlock_isHermitian _ (rdmBare_isHermitian ℓ hℓ _) _ k) := by
  sorry

end DFL
