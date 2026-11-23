import DFL.QLM

/-!
# Number entropy at the free end with fixed fermion number (paper: Remark `rem:gausscap`(c))

If the total fermion number is fixed, then the electric field at the right end is fixed in each Gauss-law sector, so
an interval `[a, N-1]` ending at the right end obeys the edge bound `ln(2S + 1)` of `qlm_numberEntropy_le_edge`: the
number of particles in `[a, N-1]` is the fixed total minus the number in `[0, a-1]`, whose distribution has entropy at
most `ln(2S + 1)`.
-/

open Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- **The edge bound at the right end** (paper: `rem:gausscap`(c)): for states of one background-charge sector with
a fixed total number of fermions, every interval `[a, N-1]` has number entropy at most `ln(2S + 1)`. -/
theorem qlm_numberEntropy_le_right_edge (S : ℕ) (E0 : ℤ) (Q : Fin N → ℤ)
    (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ)
    (hnorm : ∑ c, Complex.normSq (v c) = 1)
    (hsupp : ∀ c, v c ≠ 0 → GaussLaw S E0 Q c.1 c.2) (ntot : ℤ)
    (hfix : ∀ c, v c ≠ 0 → countInt 0 (N - 1) c.1 = ntot) (a : Fin N) :
    shannon N (qlmNumDist S a.val (N - 1) v) ≤ Real.log (2 * S + 1) := by
  sorry

end DFL
