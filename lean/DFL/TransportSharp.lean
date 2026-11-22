import DFL.DualityFull
import DFL.SectorNorm
import DFL.AreaLawSharp

/-!
# Bounded transport and number entropy with the constants of the paper
(paper: Theorems `thm:transport`(a)–(d), `thm:numberentropy`)

Constants: `Σ₁(F) = ∑_{d ≥ 1} d F(d)` (`sigma1`), `Σ'₁(F) = ∑_{d ≥ 1} d F(d)^{1/2}` (`sigma1'`),
`Θ(F) = 2 ∑_{k ≥ 1} k (2k+1) F(k)` (`thetaF`), `M_p(F) = ∑_{d ≥ 0} (d+1)^p F(d)` (`momentF`).

* (a) **Propagator**: `{c̃_x(t), c̃_y†} = ∑_q U_q(t)_{xy} P_q` (`dressed_anticomm_t`), so
  `‖P_q {c̃_x(t), c̃_y†} P_q‖ = |U_q(t)_{xy}|` (`norm_sector_anticomm`),
  `E_q sup_t ‖P_q {c̃_x(t), c̃_y†} P_q‖ ≤ F(|x-y|)` (`avg_sup_anticomm_le`) and
  `⟨Ψ₀, |{c̃_x(t), c̃_y†}| Ψ₀⟩ ≤ F(|x-y|)` (`expect_abs_anticomm_le`).
* (b) **One particle**: `∑_x (|x-y|+1)^p ⟨n_x⟩_t ≤ E_q sup_t ∑_x (|x-y|+1)^p |U_q(t)_{xy}|² ≤ 2 M_p(F)`
  for `0 ≤ p ≤ 2` (`one_particle_moments`).
* (c) **Transport across a cut**: `|⟨N_{≤m}⟩_{q,t} - s_m| ≤ ∑_{x ≤ m < z} |U_q(t)_{xz}|²`
  (`sector_transport_le`), `E_q sup_t |⟨N_{≤m}⟩_{q,t} - s_m| ≤ Σ₁(F)` (`avg_sup_transport_le`),
  `|⟨Ψ(t), N_{≤m} Ψ(t)⟩ - s_m| ≤ Σ₁(F)` (`lgt_transport_sharp`) and
  `⟨Ψ(t), (N_{≤m} - s_m)² Ψ(t)⟩ ≤ Θ(F) + Σ'₁(F)²` (`lgt_transport_sq`).
* (d) **Dressed correlations**: `|⟨Ψ(t), c̃_x† c̃_y Ψ(t)⟩| ≤ 2 ∑_{k ≥ ⌈|x-y|/2⌉} F(k)` (`dressed_corr_le`).
* **Number entropy**: `S_N(t) ≤ ½ ln(2πe(Θ(F) + Σ'₁(F)² + 1/12))` (`lgt_numberEntropy_sharp`) and, in a
  single sector, `E_q sup_t H(p^q_t) ≤ ½ ln(2πe(Θ(F) + 1/12))` (`avg_sup_sector_numberEntropy_le`).

Sites are 0-indexed: `leftOf m = {0, …, m}` is the paper's `[1, m+1]`.
-/

open Matrix Finset
open scoped Kronecker Matrix.Norms.L2Operator ComplexOrder

noncomputable section

namespace DFL

variable {N : ℕ}

/-- `Σ₁(F) = ∑_{d ≥ 1} d F(d)`. -/
def sigma1 (F : ℕ → ℝ) : ℝ := ∑' d : ℕ, (d : ℝ) * F d

/-- `Σ'₁(F) = ∑_{d ≥ 1} d F(d)^{1/2}`. -/
def sigma1' (F : ℕ → ℝ) : ℝ := ∑' d : ℕ, (d : ℝ) * √(F d)

/-- `M_p(F) = ∑_{d ≥ 0} (d+1)^p F(d)`. -/
def momentF (p : ℝ) (F : ℕ → ℝ) : ℝ := ∑' d : ℕ, ((d : ℝ) + 1) ^ p * F d

/-- The dressed fermion at time `t`, `c̃_x(t) = e^{itH} c̃_x e^{-itH}`. -/
def dressedT (J h : ℝ) (x : Fin N) (t : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  (evo (H J h) t)ᴴ * cplx (dressed x) * evo (H J h) t

/-- The anticommutator `{c̃_x(t), c̃_y†}`. -/
def dressedAnticomm (J h : ℝ) (x y : Fin N) (t : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  dressedT J h x t * cplx (dressed y)ᵀ + cplx (dressed y)ᵀ * dressedT J h x t

/-! ### (a) The dressed propagator -/

/-- **The dressed anticommutator is the one-body propagator in each sector** (paper: eq. `eq:main-anticomm`). -/
theorem dressed_anticomm_t (J h : ℝ) (x y : Fin N) (t : ℝ) :
    dressedAnticomm J h x y t = ∑ q : Cfg N, evo (hq J h q) t x y • cplx (gaussProj q) := by
  sorry

/-- `‖P_q {c̃_x(t), c̃_y†} P_q‖ = |U_q(t)_{xy}|`. -/
theorem norm_sector_anticomm (J h : ℝ) (x y : Fin N) (t : ℝ) (q : Cfg N) :
    ‖cplx (gaussProj q) * dressedAnticomm J h x y t * cplx (gaussProj q)‖ = ‖evo (hq J h q) t x y‖ := by
  sorry

/-- **Localization of the dressed propagator** (paper: eq. `eq:main-prop`). -/
theorem avg_sup_anticomm_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (x y : Fin N) :
    avgQ (fun q : Cfg N =>
      ⨆ t : ℝ, ‖cplx (gaussProj q) * dressedAnticomm J h x y t * cplx (gaussProj q)‖) ≤
      F (Nat.dist x.val y.val) := by
  sorry

/-- **`⟨Ψ₀, |{c̃_x(t), c̃_y†}| Ψ₀⟩ ≤ F(|x-y|)`** (paper: Theorem `thm:transport`(a)): the operator
`∑_q |U_q(t)_{xy}| P_q` is positive semidefinite with square `{c̃_x(t), c̃_y†}^* {c̃_x(t), c̃_y†}` (so it is the
absolute value), and its expectation in `Ψ₀ = |↑⋯↑⟩_z ⊗ ψ` is at most `F(|x-y|)`. -/
theorem expect_abs_anticomm_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (ψ : Cfg N → ℂ)
    (hψ : star ψ ⬝ᵥ ψ = 1) (x y : Fin N) (t : ℝ) :
    (∑ q : Cfg N, ((‖evo (hq J h q) t x y‖ : ℝ) : ℂ) • cplx (gaussProj q)).PosSemidef ∧
      (∑ q : Cfg N, ((‖evo (hq J h q) t x y‖ : ℝ) : ℂ) • cplx (gaussProj q)) *
          (∑ q : Cfg N, ((‖evo (hq J h q) t x y‖ : ℝ) : ℂ) • cplx (gaussProj q)) =
        (dressedAnticomm J h x y t)ᴴ * dressedAnticomm J h x y t ∧
      (expect (∑ q : Cfg N, ((‖evo (hq J h q) t x y‖ : ℝ) : ℂ) • cplx (gaussProj q)) (initZ ψ)).re ≤
        F (Nat.dist x.val y.val) := by
  sorry

/-! ### (b) One particle -/

/-- **Moments of a single particle** (paper: Theorem `thm:transport`(b)), `0 ≤ p ≤ 2`. -/
theorem one_particle_moments (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (y : Fin N) (p : ℝ)
    (hp0 : 0 ≤ p) (hp2 : p ≤ 2) (t : ℝ) :
    ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p *
        (expect (matter (diagonal fun n => occ (n x))) (lgtState J h (fockVec {y}) t)).re ≤
      avgQ (fun q : Cfg N => ⨆ s : ℝ,
        ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p * ‖evo (hq J h q) s x y‖ ^ 2) ∧
    avgQ (fun q : Cfg N => ⨆ s : ℝ,
        ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p * ‖evo (hq J h q) s x y‖ ^ 2) ≤
      2 * momentF p F := by
  sorry

/-! ### (c) Transport across a cut -/

/-- **Transport in one sector** (paper: Theorem `thm:transport`(c)). -/
theorem sector_transport_le (J h : ℝ) (q : Cfg N) (S : Finset (Fin N)) (m : ℕ) (t : ℝ) :
    |(expect (cplx (numOp (leftOf m))) (evo (dGamma (hq J h q)) t *ᵥ fockVec S)).re -
        ((S ∩ leftOf m).card : ℝ)| ≤
      ∑ x ∈ leftOf m, ∑ z ∈ univ \ leftOf m, ‖evo (hq J h q) t x z‖ ^ 2 := by
  sorry

/-- `E_q sup_t |⟨N_{≤m}⟩_{q,t} - s_m| ≤ Σ₁(F)` (paper: eq. `eq:main-transport`). -/
theorem avg_sup_transport_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ,
      |(expect (cplx (numOp (leftOf m))) (evo (dGamma (hq J h q)) t *ᵥ fockVec S)).re -
        ((S ∩ leftOf m).card : ℝ)|) ≤ sigma1 F := by
  sorry

/-- **Bounded transport with the constant of the paper** (paper: eq. `eq:main-transport`):
`|⟨Ψ(t), N_{≤m} Ψ(t)⟩ - s_m| ≤ Σ₁(F)`. -/
theorem lgt_transport_sharp (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ)
    (t : ℝ) :
    |(expect (matter (numOp (leftOf m))) (lgtState J h (fockVec S) t)).re -
        ((S ∩ leftOf m).card : ℝ)| ≤ sigma1 F := by
  sorry

/-- **Second moment of the transported charge** (paper: eq. `eq:main-transport2`):
`⟨Ψ(t), (N_{≤m} - s_m)² Ψ(t)⟩ ≤ Θ(F) + Σ'₁(F)²`. -/
theorem lgt_transport_sq (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ)
    (t : ℝ) :
    (expect (matter ((numOp (leftOf m) - ((S ∩ leftOf m).card : ℝ) • 1) *
        (numOp (leftOf m) - ((S ∩ leftOf m).card : ℝ) • 1))) (lgtState J h (fockVec S) t)).re ≤
      thetaF F + sigma1' F ^ 2 := by
  sorry

/-! ### (d) Dressed correlations -/

/-- **Dressed correlations** (paper: Theorem `thm:transport`(d)): for a Fock state and `x ≠ y`,
`|⟨Ψ(t), c̃_x† c̃_y Ψ(t)⟩| ≤ 2 ∑_{k ≥ ⌈|x-y|/2⌉} F(k)`. -/
theorem dressed_corr_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (x y : Fin N)
    (hxy : x ≠ y) (t : ℝ) :
    ‖expect (cplx ((dressed x)ᵀ * dressed y)) (lgtState J h (fockVec S) t)‖ ≤
      2 * ∑' k : ℕ, (if Nat.dist x.val y.val ≤ 2 * k then F k else 0) := by
  sorry

/-! ### Number entropy -/

/-- **Bounded number entropy with the constants of the paper** (paper: Theorem `thm:numberentropy`):
`S_N(t) ≤ ½ ln(2πe(K₁ + K₂ + 1/12))`, `K₁ = Θ(F)`, `K₂ = Σ'₁(F)²`. -/
theorem lgt_numberEntropy_sharp (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ)
    (t : ℝ) :
    shannon N (lgtNumDist J h (leftOf m) (fockVec S) t) ≤
      (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (thetaF F + sigma1' F ^ 2 + 1 / 12)) := by
  sorry

/-- **Number entropy in a single sector** (paper: Theorem `thm:numberentropy`):
`E_q sup_t H(p^q_t) ≤ ½ ln(2πe(Θ(F) + 1/12))`. -/
theorem avg_sup_sector_numberEntropy_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N))
    (m : ℕ) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ,
      shannon N (numDist (leftOf m) (evo (dGamma (hq J h q)) t *ᵥ fockVec S))) ≤
      (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (thetaF F + 1 / 12)) := by
  sorry

end DFL
