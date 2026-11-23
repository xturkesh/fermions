import DFL.ECProof
import DFL.Ring
import DFL.OpenImbalance

/-!
# Hypothesis (EC) on the open chain and on the ring (paper: Remark `rem:ringEC`)

* `ec_open_chain` (paper: Remark `rem:ringEC`(i)): the sector law of the open chain of Proposition
  `prop:variants`(i) is uniform on a parity hyperplane, so it has density at most `2` with respect to
  the i.i.d. law, and (EC) holds for it with the extra factor `2`.
* `ringEC_of_ldt` (paper: Remark `rem:ringEC`(ii)): the proof of Theorem `thm:EC` applies to the twisted
  ring operators `h_{q,w}` for pairs at ring distance `d < N/2`: the two boxes are disjoint arcs, both
  ends of each box are active, and an eigenvector of `h_{q,w}`, unwrapped to a function on `ℤ`
  (multiplied by `w` at every turn), solves the untwisted eigenvalue equation with the periodic
  potential, so the Poisson bound of `decay_of_not_res` applies; after a rotation of the ring the two
  arcs lie in `[0, N-1]` and the probabilistic estimate `two_box` applies.
  `ringEC_par_of_ldt`: the same for the ring sector law (uniform on a parity hyperplane, factor `2`).
-/

open Matrix Finset

noncomputable section

namespace DFL

/-- The ring distance `min(|x - y|, N - |x - y|)` on `ℤ/N`. -/
def ringDist {N : ℕ} (x y : Fin N) : ℕ := min (Nat.dist x.val y.val) (N - Nat.dist x.val y.val)

/-- The twisted ring operator `h_{q,w}` is symmetric. -/
theorem hqw_isHermitian {N : ℕ} [NeZero N] (J h : ℝ) (q : Cfg N) (ω : ZMod 2) :
    (hqw J h q ω).IsHermitian := by
  sorry

/-- **The open chain** (paper: Remark `rem:ringEC`(i)): under (EC), the parity-hyperplane average of
the eigenfunction correlator is at most `2 F(|x - y|)`. -/
theorem ec_open_chain (J h : ℝ) (F : ℕ → ℝ) (hEC : EC J h F) {N : ℕ} (hN : 0 < N) (π : ZMod 2)
    (x y : Fin N) :
    avgPar π (fun q : Cfg N => corrQ (hq J h q) (hq_isHermitian J h q) x y) ≤
      2 * F (Nat.dist x.val y.val) := by
  sorry

/-- **(EC) on the ring** (paper: Remark `rem:ringEC`(ii)): given the large-deviation estimate, the
sector-averaged eigenfunction correlator of the twisted ring operators `h_{q,w}` decays exponentially
in the ring distance, for pairs at ring distance `d < N/2`, uniformly in `N` and in the twist. -/
theorem ringEC_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ ∀ (N : ℕ) [NeZero N] (ω : ZMod 2) (x y : Fin N), 2 * ringDist x y < N →
      avgQ (fun q : Cfg N => corrQ (hqw J h q ω) (hqw_isHermitian J h q ω) x y) ≤
        C * Real.exp (-η * ringDist x y) := by
  sorry

/-- **(EC) for the ring sector law** (paper: Remark `rem:ringEC`(ii) and Proposition
`prop:variants`(ii)): the uniform law on a parity hyperplane costs a factor `2`. -/
theorem ringEC_par_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ ∀ (N : ℕ) [NeZero N] (ω π : ZMod 2) (x y : Fin N), 2 * ringDist x y < N →
      avgPar π (fun q : Cfg N => corrQ (hqw J h q ω) (hqw_isHermitian J h q ω) x y) ≤
        C * Real.exp (-η * ringDist x y) := by
  sorry

end DFL
