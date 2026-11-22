import DFL.SlaterBulk

/-!
# Quasi-free density matrices (paper: Lemma `lem:quasifree`, and `lem:slater`(a): `ρ_A = ρ[C_A]`)

For a matrix `C` on `ℓ²(Fin L)` let `Γ(C)` be its second quantization on the Fock space `𝓕_L`,
`Γ(C) c†_{x_1} ⋯ c†_{x_k} |0⟩ = c†(C δ_{x_1}) ⋯ c†(C δ_{x_k}) |0⟩` (`secondQ`), and (eq. `eq:main-rhoformula`)
`ρ[C] = ∑_{j=0}^{L} ((-1)^j / j!) ∑_{x_1, …, x_j} c_{x_1} ⋯ c_{x_j} Γ(C) c†_{x_j} ⋯ c†_{x_1}` (`qfState`).

* (a) For Hermitian `C` with orthonormal eigenvectors `v_i` and eigenvalues `ν_i`,
  `ρ[C] = ∏_i (ν_i c†(v_i)c(v_i) + (1 - ν_i) c(v_i)c†(v_i))` (`qfState_eq_prod`): `ρ[C]` does not depend on the
  eigenbasis; its matrix elements are polynomials in the entries of `C` (`continuous_qfState`), and `C ↦ ρ[C]` is
  uniformly continuous on `{0 ≤ C ≤ 1}` (`uniformContinuousOn_qfState`).
* (b) For `0 ≤ C ≤ 1`, `ρ[C]` is a density matrix (`qfState_posSemidef`, `trace_qfState`) with
  `S(ρ[C]) = ∑_i h(ν_i)` (`vonNeumann_qfState`); its generating function is
  `Tr(ρ[C] z^{N}) = det(1 + (z - 1) C)` (`trace_qfState_zN`).
* `lem:slater`(a): the reduced state of a Slater determinant on any region `R` is the quasi-free state with symbol
  `C_R` (`rdmF_eq_qfState`), after relabelling the modes of `R` increasingly by `Fin |R|`.
-/

open Matrix Finset
open scoped ComplexOrder

noncomputable section

namespace DFL

open AreaLawAux

variable {L : ℕ}

/-- The occupied sites of `n`. -/
def occSet (n : Cfg L) : Finset (Fin L) := univ.filter (fun x => n x = 1)

/-- **Second quantization** `Γ(C)`: `Γ(C)|n⟩ = c†(C δ_{x_1}) ⋯ c†(C δ_{x_k}) |0⟩` for the occupied sites
`x_1 < ⋯ < x_k` of `n`. -/
def secondQ (C : Matrix (Fin L) (Fin L) ℂ) : Matrix (Cfg L) (Cfg L) ℂ := fun n' n =>
  slaterAlt L (occSet n).card
    (fun j : Fin (occSet n).card => fun x => C x ((occSet n).orderEmbOfFin rfl j)) n'

/-- **The quasi-free density matrix** `ρ[C]` (paper: eq. `eq:main-rhoformula`). -/
def qfState (C : Matrix (Fin L) (Fin L) ℂ) : Matrix (Cfg L) (Cfg L) ℂ :=
  ∑ j ∈ range (L + 1), ((-1 : ℂ) ^ j / (j.factorial : ℂ)) •
    ∑ x : Fin j → Fin L,
      (List.ofFn fun i => cplx (ann (x i))).prod * secondQ C *
        (List.ofFn fun i => cplx (cre (x (Fin.rev i)))).prod

/-- The operator `z^{N}` on `𝓕_L`. -/
def zPowN (z : ℂ) : Matrix (Cfg L) (Cfg L) ℂ := diagonal (fun n => z ^ (occSet n).card)

/-- **Eigenmode form** (paper: Lemma `lem:quasifree`(a)): `ρ[C] = ∏_i (ν_i c†(v_i)c(v_i) + (1 - ν_i) c(v_i)c†(v_i))`
for Hermitian `C` (the factors commute); in particular `ρ[C]` does not depend on the choice of eigenbasis. -/
theorem qfState_eq_prod (C : Matrix (Fin L) (Fin L) ℂ) (hC : C.IsHermitian) :
    qfState C = (List.ofFn fun i : Fin L =>
      ((hC.eigenvalues i : ℝ) : ℂ) • (crOp (hC.eigenvectorBasis i) * (crOp (hC.eigenvectorBasis i))ᴴ) +
        ((1 - hC.eigenvalues i : ℝ) : ℂ) •
          ((crOp (hC.eigenvectorBasis i))ᴴ * crOp (hC.eigenvectorBasis i))).prod := by
  sorry

/-- The matrix elements of `ρ[C]` are polynomials in the entries of `C`: `C ↦ ρ[C]` is continuous. -/
theorem continuous_qfState : Continuous (qfState : Matrix (Fin L) (Fin L) ℂ → Matrix (Cfg L) (Cfg L) ℂ) := by
  sorry

/-- **Uniform continuity** (paper: Lemma `lem:quasifree`(a)) on `{0 ≤ C ≤ 1}`. -/
theorem uniformContinuousOn_qfState :
    UniformContinuousOn (qfState : Matrix (Fin L) (Fin L) ℂ → Matrix (Cfg L) (Cfg L) ℂ)
      {C | C.PosSemidef ∧ (1 - C).PosSemidef} := by
  sorry

/-- `ρ[C] ≥ 0` for `0 ≤ C ≤ 1` (paper: Lemma `lem:quasifree`(b)). -/
theorem qfState_posSemidef (C : Matrix (Fin L) (Fin L) ℂ) (hC0 : C.PosSemidef) (hC1 : (1 - C).PosSemidef) :
    (qfState C).PosSemidef := by
  sorry

/-- `Tr ρ[C] = 1` for Hermitian `C` (paper: Lemma `lem:quasifree`(b)). -/
theorem trace_qfState (C : Matrix (Fin L) (Fin L) ℂ) (hC : C.IsHermitian) : (qfState C).trace = 1 := by
  sorry

/-- **Entropy of a quasi-free state** (paper: Lemma `lem:quasifree`(b)): `S(ρ[C]) = ∑_i h(ν_i)`. -/
theorem vonNeumann_qfState (C : Matrix (Fin L) (Fin L) ℂ) (hC0 : C.PosSemidef) (hC1 : (1 - C).PosSemidef) :
    vonNeumann (qfState C) (qfState_posSemidef C hC0 hC1).isHermitian =
      ∑ i, Real.binEntropy (hC0.isHermitian.eigenvalues i) := by
  sorry

/-- **Generating function** (paper: Lemma `lem:quasifree`(b)): `Tr(ρ[C] z^{N}) = det(1 + (z - 1) C)`. -/
theorem trace_qfState_zN (C : Matrix (Fin L) (Fin L) ℂ) (hC : C.IsHermitian) (z : ℂ) :
    (qfState C * zPowN z).trace = (1 + (z - 1) • C).det := by
  sorry

/-- **Reduced states of Slater determinants are quasi-free** (paper: Lemma `lem:slater`(a), `ρ_A = ρ[C_A]`): for any
region `R`, relabelled increasingly by `Fin |R|`. -/
theorem rdmF_eq_qfState {N : ℕ} (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N)) (t : ℝ)
    (R : Finset (Fin N)) :
    (rdmF R (evo (dGamma A) t *ᵥ fockVec S)).submatrix
        (fun n : Cfg R.card => fun a : R => n ((R.orderIsoOfFin rfl).symm a))
        (fun n : Cfg R.card => fun a : R => n ((R.orderIsoOfFin rfl).symm a)) =
      qfState ((corrR A S t R).submatrix (fun i => (R.orderIsoOfFin rfl) i) (fun i => (R.orderIsoOfFin rfl) i)) := by
  sorry

end DFL
