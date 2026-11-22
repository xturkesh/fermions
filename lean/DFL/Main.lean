import DFL.Sector
import DFL.FreeDyn
import DFL.Entropy
import DFL.Imbalance

/-!
# Main results: disorder-free localization in the 1D ℤ₂ LGT (paper: Section `sec:main`)

All statements are about the actual lattice gauge theory `H J h` with initial state
`|↑…↑⟩_z ⊗ |S⟩` (all links polarized along `z`, matter in a Fock state). The duality and the sector
averaging are proven in `DFL.Duality` and `DFL.Sector`; the only external input is
Hypothesis (EC) (`DFL.EC`), an explicit argument of the theorems that need it.

Main results:
* `DFL.lgt_density`: exact sector-averaged density formula;
* `DFL.lgt_transport`: bounded particle transport across any cut, uniformly in time and volume;
* `DFL.lgt_numberEntropy`: the number entropy of a half chain is bounded uniformly in time and
  volume;
* `DFL.lgtImbalance_eq`, `DFL.lgtImbalance_timeAvg`, `DFL.lgtImbalance_lower`: the
  density-wave imbalance has a strictly positive infinite-time value for every `h ≠ 0`, uniformly in
  the volume.
-/

open Matrix Finset Filter Topology
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Sites `x ≤ m` (the region to the left of the cut after site `m`). -/
def leftOf (m : ℕ) : Finset (Fin N) := univ.filter (fun x => x.val ≤ m)

/-- Density-wave (CDW) occupied set: even sites. -/
def cdw (N : ℕ) : Finset (Fin N) := univ.filter (fun x => Even x.val)

/-- Lift of a matter observable to the LGT (identity on the links). -/
abbrev matter (D : Matrix (Cfg N) (Cfg N) ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  cplx (D ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))

/-- **Exact density formula**: `⟨n_x(t)⟩ = 2^{-N} ∑_q ∑_{z ∈ S} |(e^{-ith_q})_{xz}|²`. -/
theorem lgt_density (J h : ℝ) (S : Finset (Fin N)) (x : Fin N) (t : ℝ) :
    (expect (matter (diagonal fun n => occ (n x))) (lgtState J h (fockVec S) t)).re =
      avgQ (fun q => ∑ z ∈ S, Complex.normSq (evo (hq J h q) t x z)) := by
  sorry

/-- **Bounded transport** (paper: Theorem `thm:transport`): under Hypothesis (EC), the number of
particles transferred across any cut is bounded uniformly in time and volume. -/
theorem lgt_transport (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ)
    (t : ℝ) :
    |(expect (matter (numOp (leftOf m))) (lgtState J h (fockVec S) t)).re -
        ((S ∩ leftOf m).card : ℝ)| ≤
      2 * ∑' d : ℕ, (d : ℝ) * F d := by
  sorry

/-- Distribution of the particle number `N_R` in the LGT state at time `t`. -/
def lgtNumDist (J h : ℝ) (R : Finset (Fin N)) (ψ : Cfg N → ℂ) (t : ℝ) (k : ℕ) : ℝ :=
  (expect (matter (numProj R k)) (lgtState J h ψ t)).re

/-- The number distribution is a probability distribution on `{0, …, N}`. -/
theorem lgtNumDist_isProb (J h : ℝ) (R S : Finset (Fin N)) (t : ℝ) :
    IsProb N (lgtNumDist J h R (fockVec S) t) := by
  sorry

/-- **Bounded number entropy** (paper: Theorem `thm:numberentropy`): under Hypothesis (EC), the
Shannon entropy of the half-chain particle-number distribution is bounded uniformly in time,
volume, position of the cut and occupied set. -/
theorem lgt_numberEntropy (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ K : ℝ, ∀ (N : ℕ) (m : ℕ) (S : Finset (Fin N)) (t : ℝ),
      shannon N (lgtNumDist J h (leftOf m) (fockVec S) t) ≤
        (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (K + 1 / 12)) := by
  sorry

/-- The density-wave imbalance `I(t) = (2/N) ∑_x (-1)^x ⟨n_x(t)⟩` of the LGT. -/
def lgtImbalance (J h : ℝ) (N : ℕ) (t : ℝ) : ℝ :=
  (2 / N : ℝ) * ∑ x : Fin N, (-1 : ℝ) ^ x.val *
    (expect (matter (diagonal fun n => occ (n x))) (lgtState J h (fockVec (cdw N)) t)).re

/-- The LGT imbalance is the Bernoulli average of the one-body imbalances (`N` even). -/
theorem lgtImbalance_eq (J h : ℝ) (hN : Even N) (t : ℝ) :
    lgtImbalance J h N t = avgQ (fun q : Cfg N => imbOne J h q t) := by
  sorry

/-- The time-averaged LGT imbalance converges to `Ī = 2^{-N} ∑_q N⁻¹ ∑_α a_α(q)²`. -/
theorem lgtImbalance_timeAvg (J h : ℝ) (hJ : J ≠ 0) (hN : Even N) :
    Tendsto (fun T : ℝ => (1 / T) * ∫ t in (0 : ℝ)..T, lgtImbalance J h N t) atTop
      (𝓝 (avgQ (fun q : Cfg N => imbInf J h q))) := by
  sorry

/-- **Persistent imbalance** (paper: Theorem `thm:imbalance`): for every `h ≠ 0`, under
Hypothesis (EC), the infinite-time imbalance of the LGT is bounded below by a positive constant
independent of the (even) volume. -/
theorem lgtImbalance_lower (J h : ℝ) (hJ : J ≠ 0) (hh : h ≠ 0) {F : ℕ → ℝ} (hEC : EC J h F) :
    ∃ c : ℝ, 0 < c ∧ ∀ (N : ℕ), Even N → 2 ≤ N →
      ∀ ε > 0, ∀ᶠ T in atTop,
        c - ε ≤ (1 / T) * ∫ t in (0 : ℝ)..T, lgtImbalance J h N t := by
  sorry

end DFL
