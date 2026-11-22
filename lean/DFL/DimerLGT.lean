import DFL.DimerLower
import DFL.Dimer

/-!
# Delocalization for the dimer sector law (paper: Theorem `thm:dimer`(a), and (d) for `p > 1`)

Let `0 < 2|h| < J` and `λ = 2|h|/J ∈ (0, 1)`. In every dimer sector `q` (`IsDimer q`: `q_{2k} = q_{2k+1}`) the
one-body operator is `h_q = J(-Δ + (2h/J) diag(q))`, a random-dimer operator with values `±λ`.

* (a) **Every sector, every volume**: for `p > 0` and `1 ≤ JT ≤ min(x, N-1-x)` (0-indexed sites),
  `∫_0^∞ (dt/T) e^{-t/T} ∑_y |y - x|^p |⟨δ_y, e^{-ith_q} δ_x⟩|² ≥ c_{λ,p} (JT)^{p-1}`
  (`dimer_sector_moment_lower`), by Proposition `prop:dimerlower` after the time change `s = Jt`.
* (d) for `p > 1`: **failure of uniform localization**, sector by sector: for every choice of dimer sectors
  `q^{(N)}`, `sup_{N, x, t} ∑_y |y - x|^p |⟨δ_y, e^{-ith_{q^{(N)}}} δ_x⟩|² = ∞` (`dimer_sup_moment_unbounded`);
  consequently there is no summable profile `F` with `∑_d (d+1)^p F(d) < ∞` such that
  `E_dim Q^{[1,N]}_q(x, y) ≤ F(|x - y|)` for all even `N` (`not_EC_dimer`), i.e. Hypothesis (EC) fails for the
  dimer law. (For `1/2 < p ≤ 1` the paper uses the transport bounds of Jitomirskaya–Schulz-Baldes–Stolz,
  which are not formalized.)
-/

open Matrix Finset MeasureTheory

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The uniform average over the dimer charge configurations (the dimer sector law `μ_dim`). -/
def dimerAvg (f : Cfg N → ℝ) : ℝ :=
  (∑ q ∈ univ.filter (fun q : Cfg N => IsDimer q), f q) /
    ((univ.filter (fun q : Cfg N => IsDimer q)).card : ℝ)

/-- **Every dimer sector spreads at least diffusively, in every volume** (paper: Theorem `thm:dimer`(a)). -/
theorem dimer_sector_moment_lower (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J)
    (q : Cfg N) (hq : IsDimer q) (x : Fin N) (p T : ℝ) (hp : 0 < p) (hT1 : 1 ≤ J * T)
    (hT2 : J * T ≤ x.val) (hT3 : J * T ≤ (N : ℝ) - 1 - x.val) :
    dimerClp (2 * |h| / J) p * (J * T) ^ (p - 1) ≤
      ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) *
        ∑ y : Fin N, |((y.val : ℝ) - x.val)| ^ p * ‖evo (hq J h q) t y x‖ ^ 2 := by
  sorry

/-- **Failure of uniform localization, sector by sector** (paper: Theorem `thm:dimer`(d), `p > 1`). -/
theorem dimer_sup_moment_unbounded (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J)
    (p : ℝ) (hp : 1 < p) (qs : ∀ N : ℕ, Cfg N) (hqs : ∀ N, IsDimer (qs N)) :
    ¬ BddAbove (Set.range (fun r : (Σ N : ℕ, Fin N × ℝ) =>
      ∑ y : Fin r.1, |((y.val : ℝ) - r.2.1.val)| ^ p * ‖evo (hq J h (qs r.1)) r.2.2 y r.2.1‖ ^ 2)) := by
  sorry

/-- **Hypothesis (EC) fails for the dimer law** (paper: Theorem `thm:dimer`(d), `p > 1`): there is no profile
`F ≥ 0` with `∑_d (d+1)^p F(d) < ∞` such that `E_dim Q^{[1,N]}_q(x, y) ≤ F(|x - y|)` for all even `N`. -/
theorem not_EC_dimer (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (p : ℝ) (hp : 1 < p)
    (F : ℕ → ℝ) (hF0 : ∀ d, 0 ≤ F d) (hFs : Summable (fun d : ℕ => ((d : ℝ) + 1) ^ p * F d)) :
    ¬ ∀ (N : ℕ), Even N → ∀ x y : Fin N,
      dimerAvg (fun q : Cfg N => corrQ (hq J h q) (hq_isHermitian J h q) x y) ≤ F (Nat.dist x.val y.val) := by
  sorry

end DFL
