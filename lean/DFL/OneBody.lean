import DFL.Basic
import DFL.Evolution

/-!
# One-body input: the eigenfunction correlator and Hypothesis (EC)

The localization input of the paper (Section `sec:onebody`) is a statement about the
one-dimensional Anderson model with Bernoulli potential `h_q = -J Δ_N + 2h diag(q)`: the averaged
eigenfunction correlator decays, uniformly in the volume (Hypothesis `hyp:EC`). In the paper its
status is discussed in detail (Ge–Zhao on `ℤ`; Macera–Sodin and Macera for finite intervals). Here
it is an explicit hypothesis: every result that uses it takes `EC J h F` as an argument, so the
Lean statements are exactly as conditional as the paper's.

Main results:
* `DFL.hq_isHermitian`;
* `DFL.corrQ`: the eigenfunction correlator `Q(x,y) = ∑_α |φ_α(x)| |φ_α(y)|`;
* `DFL.norm_evo_le_corrQ`: `|(e^{-itA})_{xy}| ≤ Q(x,y)` for every `t`;
* `DFL.transProb_le_corrQ`: the time-averaged transition probability is `≤ Q(x,y)`;
* `DFL.EC`: Hypothesis (EC).
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

theorem pathAdj_isSymm (N : ℕ) : (pathAdj N).IsSymm := by
  sorry

theorem hq_isSymm (J h : ℝ) (q : Cfg N) : (hq J h q).IsSymm := by
  sorry

theorem hq_isHermitian (J h : ℝ) (q : Cfg N) : (hq J h q).IsHermitian := by
  sorry

/-- The eigenfunction correlator `Q(x, y) = ∑_α |φ_α(x)| |φ_α(y)|` of a real symmetric matrix, for
the orthonormal eigenbasis provided by Mathlib. For a Jacobi matrix (simple spectrum) this equals
`sup {|f(A)_{xy}| : |f| ≤ 1}`, the eigenfunction correlator of the literature. -/
def corrQ (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x y : Fin N) : ℝ :=
  ∑ α, |hA.eigenvectorBasis α x| * |hA.eigenvectorBasis α y|

/-- The time-averaged transition probability `R(x, y) = ∑_α φ_α(x)² φ_α(y)²`. -/
def transProb (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x y : Fin N) : ℝ :=
  ∑ α, (hA.eigenvectorBasis α x) ^ 2 * (hA.eigenvectorBasis α y) ^ 2

theorem corrQ_nonneg (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x y : Fin N) :
    0 ≤ corrQ A hA x y := by
  sorry

/-- The propagator is dominated by the eigenfunction correlator, uniformly in time. -/
theorem norm_evo_le_corrQ (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (t : ℝ)
    (x y : Fin N) : ‖evo A t x y‖ ≤ corrQ A hA x y := by
  sorry

/-- Matrix entries of a unitary propagator have modulus at most one. -/
theorem norm_evo_le_one (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (t : ℝ) (x y : Fin N) :
    ‖evo A t x y‖ ≤ 1 := by
  sorry

/-- The time-averaged transition probability is dominated by the eigenfunction correlator. -/
theorem transProb_le_corrQ (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (x y : Fin N) :
    transProb A hA x y ≤ corrQ A hA x y := by
  sorry

/-- The uniform (Bernoulli(1/2)) average over charge configurations. -/
def avgQ (f : Cfg N → ℝ) : ℝ := (1 / 2 ^ N : ℝ) * ∑ q : Cfg N, f q

/-- **Hypothesis (EC)** (paper: `hyp:EC`): the averaged eigenfunction correlator of the Bernoulli
chain `h_q` on `{0, …, N-1}` (Dirichlet ends) is bounded by a non-increasing summable profile `F`,
uniformly in the volume `N`. Expected (and proven on `ℤ` by Ge–Zhao): `F(d) = C e^{-ηd}`. -/
structure EC (J h : ℝ) (F : ℕ → ℝ) : Prop where
  nonneg : ∀ d, 0 ≤ F d
  antitone : Antitone F
  bound : ∀ (N : ℕ) (x y : Fin N),
    avgQ (fun q => corrQ (hq J h q) (hq_isHermitian J h q) x y) ≤ F (Nat.dist x.val y.val)
  summable : Summable (fun d : ℕ => ((d : ℝ) + 1) ^ 2 * F d)
  summable_sqrt : Summable (fun d : ℕ => ((d : ℝ) + 1) ^ 2 * Real.sqrt (F d))

end DFL
