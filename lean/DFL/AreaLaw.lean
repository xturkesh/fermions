import DFL.FreeDyn
import DFL.Entropy
import DFL.OneBody

/-!
# Entanglement area law (paper: Theorem `thm:arealaw`, part on `S_G`)

For a Fock state evolved by a quadratic Hamiltonian (a Slater determinant), the von Neumann
entropy of a left segment `A = {0, …, ℓ-1}` satisfies
`S(ρ_A) = ∑_k h(ν_k) ≤ 2 ∑_k √(ν_k (1 - ν_k)) ≤ 2 ∑_{x ∈ A, y ∉ A} |C_{xy}|`,
where `ν_k` are the eigenvalues of the restricted correlation matrix `C_A`. Combined with the
duality, the sector averaging and Hypothesis (EC), this gives a uniform-in-time area law for the
sector-averaged (charge-conditioned) entanglement `S_G` of the lattice gauge theory.

For a left segment the Jordan–Wigner occupation basis factorizes canonically,
`|n⟩ = |n_A⟩ ⊗ |n_B⟩`, so the fermionic reduced density matrix is the ordinary partial trace in this
basis.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The left segment `{0, …, ℓ-1}`. -/
def segment (ℓ : ℕ) : Finset (Fin N) := univ.filter (fun x => x.val < ℓ)

/-- Splitting a configuration into its left segment `{0, …, ℓ-1}` and the rest. -/
def splitCfg (ℓ : ℕ) (hℓ : ℓ ≤ N) : Cfg N ≃ Cfg ℓ × Cfg (N - ℓ) :=
  (Equiv.arrowCongr (finCongr (by omega : N = ℓ + (N - ℓ))) (Equiv.refl _)).trans
    ((Equiv.arrowCongr finSumFinEquiv.symm (Equiv.refl _)).trans
      (Equiv.sumArrowEquivProdArrow _ _ _))

/-- Reduced density matrix of a pure state `v` on the left segment `{0, …, ℓ-1}`. -/
def rdm (ℓ : ℕ) (hℓ : ℓ ≤ N) (v : Cfg N → ℂ) : Matrix (Cfg ℓ) (Cfg ℓ) ℂ :=
  fun a a' => ∑ b : Cfg (N - ℓ),
    v ((splitCfg ℓ hℓ).symm (a, b)) * star (v ((splitCfg ℓ hℓ).symm (a', b)))

theorem rdm_isHermitian (ℓ : ℕ) (hℓ : ℓ ≤ N) (v : Cfg N → ℂ) : (rdm ℓ hℓ v).IsHermitian := by
  sorry

/-- von Neumann entropy (in nats) of a Hermitian matrix: `∑_i -λ_i ln λ_i`. -/
def vonNeumann {ι : Type*} [Fintype ι] [DecidableEq ι] (ρ : Matrix ι ι ℂ) (hρ : ρ.IsHermitian) :
    ℝ :=
  ∑ i, Real.negMulLog (hρ.eigenvalues i)

/-- **Entanglement of Slater determinants**: the entanglement entropy of a left segment is
bounded by `2 ∑_{x ∈ A, y ∉ A} |C_{xy}|` (the paper has the sharper constant `2 ln 2`). -/
theorem vonNeumann_slater_le (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    vonNeumann (rdm ℓ hℓ (evo (dGamma A) t *ᵥ fockVec S))
        (rdm_isHermitian ℓ hℓ _) ≤
      2 * ∑ x ∈ segment ℓ, ∑ y ∈ univ \ segment ℓ, ‖corrMat A S t x y‖ := by
  sorry

/-- **Area law for the charge-conditioned entanglement** `S_G` (paper: Theorem `thm:arealaw`):
under Hypothesis (EC), the Bernoulli average over charge sectors of the entanglement entropy of a
left segment is bounded uniformly in time, volume, segment length and occupied set. -/
theorem lgt_areaLaw (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ K : ℝ, ∀ (N : ℕ) (ℓ : ℕ) (hℓ : ℓ ≤ N) (S : Finset (Fin N)) (t : ℝ),
      avgQ (fun q : Cfg N =>
        vonNeumann (rdm ℓ hℓ (evo (dGamma (hq J h q)) t *ᵥ fockVec S))
          (rdm_isHermitian ℓ hℓ _)) ≤ K := by
  sorry

end DFL
