import DFL.Duality

/-!
# The standard ℤ₂ gauge theory: one charge in a random-walk potential (paper: Remark
`rem:randomwalk`)

With the single-link electric term `-h ∑_j σˣ_{b_j}` (instead of the two-link term of
Smith–Knolle–Kovrizhin–Moessner) the dual model is no longer free. In the one-particle sector,
however, it is an exact single-particle problem: in the charge sector `q`, a fermion at site `y`
feels the potential
`V_q(y) = 2h ∑_{j ≥ y} Q_j - h ∑_j Q_j`, with `Q_j = ∏_{k ≤ j} q_k`,
a random walk in `y` with i.i.d. ±1 increments when the charges are i.i.d. (the map `q ↦ Q` is a
bijection).
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The standard ℤ₂ LGT Hamiltonian with single-link electric term:
`H = -J ∑_j σᶻ_{b_j} (c_j† c_{j+1} + h.c.) - h ∑_j σˣ_{b_j}`. -/
def Hstd (J h : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  (-J) • ∑ j, hop j ⊗ₖ sz j - h • ∑ j, (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ sx j

/-- The configuration with a single fermion at site `y`. -/
def oneP (y : Fin N) : Cfg N := Pi.single y 1

/-- `Q_j = ∏_{k ≤ j} q_k`, as a sign. -/
def Qsgn (q : Cfg N) (j : Fin N) : ℝ := sgn (∑ k ∈ univ.filter (fun k => k ≤ j), q k)

/-- The random-walk potential of the charge sector `q`. -/
def rwPot (h : ℝ) (q : Cfg N) (y : Fin N) : ℝ :=
  2 * h * ∑ j ∈ univ.filter (fun j => y ≤ j), Qsgn q j - h * ∑ j, Qsgn q j

/-- **One charge in a random-walk potential**: in the dual frame, the one-particle block of the
standard ℤ₂ LGT in the sector `q` is `-J Δ + diag(V_q)`. -/
theorem Hstd_dual_oneP (J h : ℝ) (y y' : Fin N) (q q' : Cfg N) :
    Hstd J h (dualEquiv (oneP y', q')) (dualEquiv (oneP y, q)) =
      if q' = q then (-J) * pathAdj N y' y + (if y' = y then rwPot h q y else 0) else 0 := by
  sorry

/-- The partial-sum map `q ↦ (∑_{k ≤ j} q_k)_j` is a bijection, so `Q` is i.i.d. uniform when `q` is. -/
theorem partialSum_bijective :
    Function.Bijective (fun q : Cfg N => fun j : Fin N => ∑ k ∈ univ.filter (fun k => k ≤ j), q k) := by
  sorry

end DFL
