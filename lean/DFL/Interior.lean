import DFL.Entropies

/-!
# Interior intervals (paper: Remark `rem:interior`)

For an interval `I = [ℓ₁, ℓ₁ + m)` (sites and links; Lean indexing, `m ≥ 1`) of the open chain, the
duality `U`, `(UΦ)(n, x) = Φ(n, chargeOf n x)`, changes the entanglement entropy of `I` by at most
`2 ln 2`: `|S_I(UΦ) - S_{I'}(Φ)| ≤ 2 ln 2` (`interior_defect`), with the partial trace in the
occupation basis (the Jordan–Wigner partial trace) on both sides. In configuration language:
* the right gate (control: the last link of `I`, target: the links right of `I`) is removed by a
  substitution in the sum over the right part, up to the pinching by the label of the last link of
  `I`; this costs at most `ln 2` (as in `cutlocal_defect`);
* the left gate (control: the last link left of `I`, target: all links right of it) splits the
  reduced state as `ω₀ + F ω₁ F†`, where `ω_c` is the part with control label `c` and `F` flips all
  link labels of `I`; concavity and the mixing bound give a change of at most `ln 2`.
The value `2 ln 2` is attained for `N = 3`, `I` = site `1` and link `1` (paper: site 2 and link `b_2`):
`exInt` is a product state across `I'` whose image under `U` has maximally mixed reduced state on
`I` (`exInt_entropy`, `exInt_U_entropy`). For states of definite total fermion parity the reduced
state of `I` is block diagonal with respect to the parity of `I` (`rdmMid_parity_block`).
-/

open Matrix Finset
open scoped ComplexOrder

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Split an LGT configuration `(n, x)` into the interval `I = [ℓ₁, ℓ₁ + m)` (sites and links), the
left part `[0, ℓ₁)` and the right part `[ℓ₁ + m, N)`. -/
def splitMid (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) :
    Cfg N × Cfg N ≃
      (Cfg m × Cfg m) × ((Cfg ℓ₁ × Cfg ℓ₁) × (Cfg (N - ℓ₁ - m) × Cfg (N - ℓ₁ - m))) :=
  (splitLGT ℓ₁ (by omega)).trans <|
    (Equiv.prodCongr (Equiv.refl _) (splitLGT (N := N - ℓ₁) m (by omega))).trans
      { toFun := fun p => (p.2.1, (p.1, p.2.2))
        invFun := fun p => (p.2.1, (p.1, p.2.2))
        left_inv := fun _ => rfl
        right_inv := fun _ => rfl }

/-- The reduced density matrix of the interval `[ℓ₁, ℓ₁ + m)`: partial trace over the left and the
right part in the occupation basis (the Jordan–Wigner partial trace). -/
def rdmMid (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    Matrix (Cfg m × Cfg m) (Cfg m × Cfg m) ℂ :=
  fun a a' => ∑ b, Ψ ((splitMid ℓ₁ m h).symm (a, b)) * star (Ψ ((splitMid ℓ₁ m h).symm (a', b)))

/-- The reduced density matrix of an interval is Hermitian. -/
theorem rdmMid_isHermitian (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    (rdmMid ℓ₁ m h Ψ).IsHermitian := by
  sorry

/-- The reduced density matrix of an interval is positive semidefinite. -/
theorem rdmMid_posSemidef (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    (rdmMid ℓ₁ m h Ψ).PosSemidef := by
  sorry

/-- The trace of the reduced density matrix is the squared norm of the state. -/
theorem trace_rdmMid (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Ψ : Cfg N × Cfg N → ℂ) :
    (rdmMid ℓ₁ m h Ψ).trace = ((∑ p, Complex.normSq (Ψ p) : ℝ) : ℂ) := by
  sorry

/-- **Interior intervals** (paper: Remark `rem:interior`): for every normalized state `Φ` of the dual
variables `(n, q)`, the entanglement entropy of the interval `[ℓ₁, ℓ₁ + m)` (sites and links) in
`UΦ`, `(UΦ)(n, x) = Φ(n, chargeOf n x)`, and that of the interval (modes and charges) in `Φ` differ
by at most `2 ln 2`. -/
theorem interior_defect (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (hm : 0 < m) (Φ : Cfg N × Cfg N → ℂ)
    (hΦ : ∑ p, Complex.normSq (Φ p) = 1) :
    |vonNeumann (rdmMid ℓ₁ m h (fun p => Φ (dualEquiv.symm p))) (rdmMid_isHermitian ℓ₁ m h _) -
      vonNeumann (rdmMid ℓ₁ m h Φ) (rdmMid_isHermitian ℓ₁ m h _)| ≤ 2 * Real.log 2 := by
  sorry

/-- **Parity blocks** (paper: Remark `rem:interior`): if `Ψ` has definite total fermion parity, its
reduced state on an interval vanishes between configurations of the interval of different fermion
parity. -/
theorem rdmMid_parity_block (ℓ₁ m : ℕ) (h : ℓ₁ + m ≤ N) (Ψ : Cfg N × Cfg N → ℂ) (π : ZMod 2)
    (hΨ : ∀ p, Ψ p ≠ 0 → ∑ k, p.1 k = π) (a a' : Cfg m × Cfg m)
    (ha : ∑ k, a.1 k ≠ ∑ k, a'.1 k) : rdmMid ℓ₁ m h Ψ a a' = 0 := by
  sorry

/-! ### Sharpness -/

/-- `Φ = ½ |n₀ = 0, n₂ = 0⟩ ⊗ (|n₁ = 0⟩ + |n₁ = 1⟩) ⊗ (|q₀ = 0⟩ + |q₀ = 1⟩) ⊗ |q₁ = 0, q₂ = 0⟩`
(`N = 3`), a product state across the interval `I'` = (mode `1`, charge `1`). -/
def exInt : Cfg 3 × Cfg 3 → ℂ :=
  fun p => if p.1 0 = 0 ∧ p.1 2 = 0 ∧ p.2 1 = 0 ∧ p.2 2 = 0 then 1 / 2 else 0

/-- `exInt` is normalized. -/
theorem exInt_norm : ∑ p, Complex.normSq (exInt p) = 1 := by
  sorry

/-- `S_{I'}(Φ) = 0` for `Φ = exInt`, `I'` = (mode `1`, charge `1`). -/
theorem exInt_entropy :
    vonNeumann (rdmMid 1 1 (by norm_num) exInt) (rdmMid_isHermitian 1 1 (by norm_num) _) = 0 := by
  sorry

/-- `S_I(UΦ) = 2 ln 2` for `Φ = exInt`, `I` = (site `1`, link `1`): the reduced state is maximally
mixed. -/
theorem exInt_U_entropy :
    vonNeumann (rdmMid 1 1 (by norm_num) (fun p => exInt (dualEquiv.symm p)))
      (rdmMid_isHermitian 1 1 (by norm_num) _) = 2 * Real.log 2 := by
  sorry

/-- **The bound `2 ln 2` of `interior_defect` is attained.** -/
theorem interior_defect_attained :
    |vonNeumann (rdmMid 1 1 (by norm_num) (fun p => exInt (dualEquiv.symm p)))
        (rdmMid_isHermitian 1 1 (by norm_num) _) -
      vonNeumann (rdmMid 1 1 (by norm_num) exInt) (rdmMid_isHermitian 1 1 (by norm_num) _)| =
      2 * Real.log 2 := by
  sorry

end DFL
