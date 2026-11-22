import DFL.Entropy
import DFL.Basic

/-!
# Number entropy in U(1) quantum link models: kinematic bound (paper: Proposition `prop:gausscap`)

In a U(1) quantum link model with spin-`S` links, in a fixed sector of background charges, the
summed Gauss law expresses the matter particle number of an interval through the electric fields
on its two boundary links. Hence, for *any* state of the sector (and therefore for any dynamics),
the particle-number distribution of an interval is supported on at most `4S + 1` values, and its
Shannon entropy is at most `ln(4S + 1)`.

Conventions: staggered fermions with charge `ρ_j = n_j - [j odd]`; electric fields
`E_b ∈ {-S, …, S}` on the links `b_j = (j, j+1)`, stored as `E_b + S ∈ Fin (2S+1)`; a fixed
boundary field `E_{-1}`; integer background charges `Q_j`; Gauss law
`E_{b_j} - E_{b_{j-1}} = ρ_j + Q_j`.
-/

open Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Electric field `E ∈ {-S, …, S}` encoded by `e ∈ Fin (2S+1)`. -/
def efield (S : ℕ) (e : Fin (2 * S + 1)) : ℤ := (e.val : ℤ) - S

/-- Electric field on the link to the left of site `j` (the boundary field `E0` for `j = 0`). -/
def prevField (S : ℕ) (E0 : ℤ) (e : Fin N → Fin (2 * S + 1)) (j : Fin N) : ℤ :=
  if hj : j.val = 0 then E0 else efield S (e ⟨j.val - 1, by have := j.isLt; omega⟩)

/-- Staggered matter charge `ρ_j = n_j - [j odd]`. -/
def stagCharge (n : Cfg N) (j : Fin N) : ℤ :=
  (if n j = 0 then 0 else 1) - (if j.val % 2 = 1 then 1 else 0)

/-- The Gauss law `E_{b_j} - E_{b_{j-1}} = ρ_j + Q_j` for all sites. -/
def GaussLaw (S : ℕ) (E0 : ℤ) (Q : Fin N → ℤ) (n : Cfg N) (e : Fin N → Fin (2 * S + 1)) :
    Prop :=
  ∀ j : Fin N, efield S (e j) - prevField S E0 e j = stagCharge n j + Q j

/-- Number of particles of `n` in the interval `[a, b]`, as an integer. -/
def countInt (a b : ℕ) (n : Cfg N) : ℤ :=
  ∑ j ∈ univ.filter (fun j : Fin N => a ≤ j.val ∧ j.val ≤ b), (if n j = 0 then 0 else 1)

/-- **Summed Gauss law**: the particle number of `[a, b]` is fixed by the boundary fields. -/
theorem countInt_gauss (S : ℕ) (E0 : ℤ) (Q : Fin N → ℤ) (n : Cfg N)
    (e : Fin N → Fin (2 * S + 1)) (hG : GaussLaw S E0 Q n e) (a b : Fin N) (hab : a ≤ b) :
    countInt a.val b.val n =
      efield S (e b) - prevField S E0 e a +
        ∑ j ∈ univ.filter (fun j : Fin N => a.val ≤ j.val ∧ j.val ≤ b.val),
          ((if j.val % 2 = 1 then 1 else 0) - Q j) := by
  sorry

/-- Distribution of the particle number of `[a, b]` in a state `v` of the quantum link model. -/
def qlmNumDist (S : ℕ) (a b : ℕ) (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ) (k : ℕ) : ℝ :=
  ∑ c ∈ univ.filter (fun c : Cfg N × (Fin N → Fin (2 * S + 1)) => countInt a b c.1 = k),
    Complex.normSq (v c)

/-- **Kinematic cap on the number entropy** (paper: Proposition `prop:gausscap`): for every
normalized state supported on a fixed background-charge sector, the number entropy of an interval
is at most `ln(4S + 1)`. -/
theorem qlm_numberEntropy_le (S : ℕ) (E0 : ℤ) (Q : Fin N → ℤ)
    (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ)
    (hnorm : ∑ c, Complex.normSq (v c) = 1)
    (hsupp : ∀ c, v c ≠ 0 → GaussLaw S E0 Q c.1 c.2) (a b : Fin N) (hab : a ≤ b) :
    shannon N (qlmNumDist S a.val b.val v) ≤ Real.log (4 * S + 1) := by
  sorry

/-- At a fixed-flux edge (`a = 0`) the cap improves to `ln(2S + 1)`. -/
theorem qlm_numberEntropy_le_edge (S : ℕ) (E0 : ℤ) (Q : Fin N → ℤ)
    (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ)
    (hnorm : ∑ c, Complex.normSq (v c) = 1)
    (hsupp : ∀ c, v c ≠ 0 → GaussLaw S E0 Q c.1 c.2) (b : Fin N) :
    shannon N (qlmNumDist S 0 b.val v) ≤ Real.log (2 * S + 1) := by
  sorry

end DFL
