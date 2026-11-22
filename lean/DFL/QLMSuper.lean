import DFL.QLM

/-!
# Quantum links: superpositions of background-charge sectors and arbitrary dynamics
(paper: Proposition `prop:gausscap`(c),(d))

The summed Gauss law `N_A = E_b - E_{a-1} + o_A - Q_A` holds configuration by configuration, with
`Q_A = ∑_{j ∈ A} G_j` the total background charge of `A = [a, b]`. Hence, for every normalized state
(not necessarily in one sector), `S_N ≤ ln(4S+1) + H(Q_A)` in the bulk and `S_N ≤ ln(2S+1) + H(Q_A)` at
the fixed-flux edge. Every unitary commuting with the Gauss operators preserves each sector and the
law of `Q_A`, so these bounds, and the bound of `qlm_numberEntropy_le` for states in one sector, hold
along any gauge-invariant dynamics.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The Gauss operator `G_j = E_j - E_{j-1} - ρ_j` evaluated on a configuration. -/
def gaussVal (S : ℕ) (E0 : ℤ) (c : Cfg N × (Fin N → Fin (2 * S + 1))) (j : Fin N) : ℤ :=
  efield S (c.2 j) - prevField S E0 c.2 j - stagCharge c.1 j

/-- `G_j` as a diagonal matrix. -/
def gaussQLM (S : ℕ) (E0 : ℤ) (j : Fin N) :
    Matrix (Cfg N × (Fin N → Fin (2 * S + 1))) (Cfg N × (Fin N → Fin (2 * S + 1))) ℂ :=
  diagonal (fun c => ((gaussVal S E0 c j : ℤ) : ℂ))

/-- The total background charge `Q_A = ∑_{j ∈ [a, b]} G_j` of a configuration. -/
def bgCharge (S : ℕ) (E0 : ℤ) (a b : ℕ) (c : Cfg N × (Fin N → Fin (2 * S + 1))) : ℤ :=
  ∑ j ∈ univ.filter (fun j : Fin N => a ≤ j.val ∧ j.val ≤ b), gaussVal S E0 c j

/-- The Shannon entropy `H(Q_A)` of the law of `Q_A` in the state `v`. -/
def bgEntropy (S : ℕ) (E0 : ℤ) (a b : ℕ) (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ) : ℝ :=
  ∑ y ∈ univ.image (bgCharge (N := N) S E0 a b),
    Real.negMulLog (∑ c ∈ univ.filter (fun c => bgCharge S E0 a b c = y), Complex.normSq (v c))

/-- `GaussLaw` is the statement `G_j = Q_j` for all `j`. -/
theorem gaussLaw_iff (S : ℕ) (E0 : ℤ) (Q : Fin N → ℤ) (c : Cfg N × (Fin N → Fin (2 * S + 1))) :
    GaussLaw S E0 Q c.1 c.2 ↔ ∀ j, gaussVal S E0 c j = Q j := by
  sorry

/-- **Superpositions of sectors, bulk** (paper: Proposition `prop:gausscap`(c)):
`S_N ≤ ln(4S+1) + H(Q_A)` for every normalized state and every interval `[a, b]` with `a ≥ 1`. -/
theorem qlm_numberEntropy_le_super (S : ℕ) (E0 : ℤ) (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ)
    (hnorm : ∑ c, Complex.normSq (v c) = 1) (a b : Fin N) (hab : a ≤ b) (ha : a.val ≠ 0) :
    shannon N (qlmNumDist S a.val b.val v) ≤
      Real.log (4 * S + 1) + bgEntropy S E0 a.val b.val v := by
  sorry

/-- **Superpositions of sectors, fixed-flux edge** (paper: Proposition `prop:gausscap`(c)):
`S_N ≤ ln(2S+1) + H(Q_A)` for every normalized state and every interval `[0, b]`. -/
theorem qlm_numberEntropy_le_super_edge (S : ℕ) (E0 : ℤ)
    (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ) (hnorm : ∑ c, Complex.normSq (v c) = 1)
    (b : Fin N) :
    shannon N (qlmNumDist S 0 b.val v) ≤ Real.log (2 * S + 1) + bgEntropy S E0 0 b.val v := by
  sorry

/-- **Any gauge-invariant dynamics preserves the sectors** (paper: Proposition `prop:gausscap`(d)). -/
theorem qlm_sector_preserved (S : ℕ) (E0 : ℤ) (Q : Fin N → ℤ)
    (U : Matrix (Cfg N × (Fin N → Fin (2 * S + 1))) (Cfg N × (Fin N → Fin (2 * S + 1))) ℂ)
    (hU : ∀ j, U * gaussQLM S E0 j = gaussQLM S E0 j * U)
    (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ) (hsupp : ∀ c, v c ≠ 0 → GaussLaw S E0 Q c.1 c.2) :
    ∀ c, (U *ᵥ v) c ≠ 0 → GaussLaw S E0 Q c.1 c.2 := by
  sorry

/-- **Any gauge-invariant dynamics preserves the law of `Q_A`** (paper: Proposition
`prop:gausscap`(d)). -/
theorem bgEntropy_invariant (S : ℕ) (E0 : ℤ) (a b : ℕ)
    (U : Matrix (Cfg N × (Fin N → Fin (2 * S + 1))) (Cfg N × (Fin N → Fin (2 * S + 1))) ℂ)
    (hU : ∀ j, U * gaussQLM S E0 j = gaussQLM S E0 j * U) (hUu : Uᴴ * U = 1)
    (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ) :
    bgEntropy S E0 a b (U *ᵥ v) = bgEntropy S E0 a b v := by
  sorry

end DFL
