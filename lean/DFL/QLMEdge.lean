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
  have hN : 0 < N := Nat.lt_of_le_of_lt (Nat.zero_le _) a.isLt
  obtain ⟨b, hb⟩ : ∃ b : Fin N, b.val = N - 1 := ⟨⟨N - 1, by omega⟩, rfl⟩
  obtain ⟨z, hz⟩ : ∃ z : Fin N, z.val = 0 := ⟨⟨0, hN⟩, rfl⟩
  have hab : a ≤ b := by
    rw [Fin.le_iff_val_le_val, hb]
    have := a.isLt
    omega
  have hzb : z ≤ b := by
    rw [Fin.le_iff_val_le_val, hz]
    exact Nat.zero_le _
  rw [← hb]
  by_cases ha : a.val = 0
  · -- `a = 0`: the particle number of `[0, N-1]` is the fixed total `ntot`
    have h := shannon_qlmNumDist_le_of_window S a.val b.val v hnorm (ntot - S) (2 * S)
      (fun c hc => by
        rw [ha, hb, hfix c hc]
        constructor <;> push_cast <;> omega)
    push_cast at h
    exact h
  · -- `a > 0`: the field at the right end is fixed by `ntot`, the one left of `a` is in `[-S, S]`
    set Ca : ℤ := ∑ j ∈ univ.filter (fun j : Fin N => a.val ≤ j.val ∧ j.val ≤ b.val),
      ((if j.val % 2 = 1 then 1 else 0) - Q j) with hCa
    set C0 : ℤ := ∑ j ∈ univ.filter (fun j : Fin N => z.val ≤ j.val ∧ j.val ≤ b.val),
      ((if j.val % 2 = 1 then 1 else 0) - Q j) with hC0
    have h := shannon_qlmNumDist_le_of_window S a.val b.val v hnorm
      (ntot + E0 - C0 + Ca - S) (2 * S)
      (fun c hc => by
        have hGa := countInt_gauss S E0 Q c.1 c.2 (hsupp c hc) a b hab
        have hG0 := countInt_gauss S E0 Q c.1 c.2 (hsupp c hc) z b hzb
        have hp0 : prevField S E0 c.2 z = E0 := by
          unfold prevField
          rw [dite_eq_left hz]
        have hpa := prevField_bounds S E0 c.2 a ha
        have hf : countInt z.val b.val c.1 = ntot := by
          rw [hz, hb]
          exact hfix c hc
        rw [← hCa] at hGa
        rw [← hC0, hp0, hf] at hG0
        rw [hGa]
        push_cast
        constructor <;> linarith [hpa.1, hpa.2])
    push_cast at h
    exact h

end DFL
