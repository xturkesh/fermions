import DFL.Main
import DFL.Variants

/-!
# Persistent imbalance on the open chain (paper: Corollary `cor:imbalanceopen`)

The model of Smith–Knolle–Kovrizhin–Moessner without the boundary link is the open chain of
`DFL/Variants.lean` (`openInit`: links `0, …, N-2` in `|↑⟩_z`, boundary link in `|+⟩`). Its sector
law is uniform on the parity hyperplane `{q : ∑_j q_j = π}`, `π` the parity of the matter state, so
the imbalance is the parity-hyperplane average `E'_q` of the sector imbalances.

Three properties of `E'` replace the independence of the charges in `avg_imbInf_uniform`: the
charges of any `N - 1` sites are i.i.d. uniform under `E'`; `E' f ≤ 2 E f` for `f ≥ 0`; and the
charges are pairwise independent (`N ≥ 3`). Hence, under (EC), the infinite-time imbalance of the
open chain is bounded below uniformly in the volume.
-/

open Matrix Finset Filter Topology

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The average over the parity hyperplane, `E'_q f = 2^{-(N-1)} ∑_{∑ q = π} f(q)`. -/
def avgPar (π : ZMod 2) (f : Cfg N → ℝ) : ℝ :=
  (1 / 2 ^ (N - 1) : ℝ) * ∑ q ∈ univ.filter (fun q : Cfg N => ∑ j, q j = π), f q

/-- The density-wave imbalance of the open chain. -/
def lgtImbalanceOpen (J h : ℝ) (N : ℕ) (hN : 0 < N) (t : ℝ) : ℝ :=
  (2 / N : ℝ) * ∑ x : Fin N, (-1 : ℝ) ^ x.val *
    (expect (matter (diagonal fun n => occ (n x)))
      (evo (H J h) t *ᵥ openInit hN (fockVec (cdw N)))).re

/-- The imbalance of the open chain is the parity-hyperplane average of the sector imbalances. -/
theorem lgtImbalanceOpen_eq (J h : ℝ) (hN0 : 0 < N) (hN : Even N) (t : ℝ) :
    lgtImbalanceOpen J h N hN0 t =
      avgPar (∑ j, fockCfg (cdw N) j) (fun q : Cfg N => imbOne J h q t) := by
  sorry

/-- The time average of the imbalance of the open chain converges to the parity-hyperplane average of
the infinite-time sector imbalances. -/
theorem lgtImbalanceOpen_timeAvg (J h : ℝ) (hJ : J ≠ 0) (hN0 : 0 < N) (hN : Even N) :
    Tendsto (fun T : ℝ => (1 / T) * ∫ t in (0 : ℝ)..T, lgtImbalanceOpen J h N hN0 t) atTop
      (𝓝 (avgPar (∑ j, fockCfg (cdw N) j) (fun q : Cfg N => imbInf J h q))) := by
  sorry

/-- **Uniform lower bound for the parity-hyperplane average** under Hypothesis (EC). -/
theorem avgPar_imbInf_uniform (J h : ℝ) (hJ : J ≠ 0) (hh : h ≠ 0) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ c : ℝ, 0 < c ∧ ∀ (N : ℕ), 4 ≤ N → ∀ π : ZMod 2,
      c ≤ avgPar π (fun q : Cfg N => imbInf J h q) := by
  sorry

/-- **Persistent imbalance on the open chain** (paper: Corollary `cor:imbalanceopen`): under
Hypothesis (EC), the time-averaged density-wave imbalance of the model of Smith–Knolle–Kovrizhin–
Moessner without boundary link is eventually above a positive constant independent of the volume. -/
theorem lgtImbalanceOpen_lower (J h : ℝ) (hJ : J ≠ 0) (hh : h ≠ 0) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ c : ℝ, 0 < c ∧ ∀ (N : ℕ) (hN0 : 0 < N), Even N → 4 ≤ N →
      ∀ ε > 0, ∀ᶠ T in atTop,
        c - ε ≤ (1 / T) * ∫ t in (0 : ℝ)..T, lgtImbalanceOpen J h N hN0 t := by
  sorry

end DFL
