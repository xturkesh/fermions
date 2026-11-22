import DFL.InfVolTransport
import DFL.AreaLawSharp

/-!
# Number entropy of bulk intervals in the thermodynamic limit (paper: Theorem `thm:infvol`(e))

For a finite interval `A = [a, b] ⊂ ℤ`, a Fock configuration `n` of the matter on `ℤ` and the gauge theory
on `Λ_M = [-M, M]` with Fock matter `n|_{Λ_M}` and `z`-polarized links, let `p^{(M)}_t` be the law of `N_A`
in `Ψ^{(M)}(t)` (`numDistM`). With `K^ℤ_q(t) = U^ℤ_q(t) P_S U^ℤ_q(t)^*` (`KZ`, `S = {s : n_s = 1}`) and
`p_t := E PB(ν(K^ℤ_q(t)|_A))` (`pZ`), the Poisson-binomial law of the eigenvalues of `K^ℤ_q(t)|_A` averaged
over i.i.d. charges on `ℤ`:

* `p^{(M)}_t(k) → p_t(k)` for every `k` (`tendsto_numDistM`);
* `S_N^{(M)}(t) → S_N(t) := H(p_t)` (`tendsto_numberEntropyM`);
* `S_N(t) ≤ ½ ln(2πe(2Θ(F) + 4Σ'₁(F)² + 1/12))` for every `t` under (EC) (`numberEntropyZ_le`).
-/

open Matrix Finset MeasureTheory Filter Topology

noncomputable section

namespace DFL

/-- The bulk interval `[a, b]` inside `Λ_M`, relabelled in `Fin (2M+1)`. -/
def bulkM (M : ℕ) (a b : ℤ) : Finset (Fin (2 * M + 1)) :=
  univ.filter (fun i => a ≤ (i.val : ℤ) - M ∧ (i.val : ℤ) - M ≤ b)

/-- The law of `N_A`, `A = [a, b]`, in the gauge theory on `Λ_M` with Fock matter `n|_{Λ_M}` and
`z`-polarized links. -/
def numDistM (J h : ℝ) (n : CfgZ) (a b : ℤ) (t : ℝ) (M : ℕ) (k : ℕ) : ℝ :=
  lgtNumDist J h (bulkM M a b) (fockVec (occM M n)) t k

/-- `K^ℤ_q(t)_{xy} = ∑_{s : n_s = 1} U^ℤ_q(t)_{xs} conj(U^ℤ_q(t)_{ys})`. -/
def KZ (J h : ℝ) (q n : CfgZ) (t : ℝ) (x y : ℤ) : ℂ :=
  ∑' s : ℤ, if n s = 1 then UZ J h q t x s * star (UZ J h q t y s) else 0

/-- `K^ℤ_q(t)|_A`, `A = [a, b]`. -/
def KZA (J h : ℝ) (q n : CfgZ) (t : ℝ) (a b : ℤ) : Matrix (Set.Icc a b) (Set.Icc a b) ℂ :=
  fun x y => KZ J h q n t x y

/-- `K^ℤ_q(t)|_A` is Hermitian. -/
theorem KZA_isHermitian (J h : ℝ) (q n : CfgZ) (t : ℝ) (a b : ℤ) : (KZA J h q n t a b).IsHermitian := by
  ext x y
  simp only [Matrix.conjTranspose_apply, KZA, KZ]
  rw [tsum_star]
  refine tsum_congr fun s => ?_
  split_ifs
  · rw [star_mul', star_star, mul_comm]
  · exact star_zero _

/-- `p_t(k) = E PB(ν(K^ℤ_q(t)|_A))(k)`. -/
def pZ (J h : ℝ) (n : CfgZ) (a b : ℤ) (t : ℝ) (k : ℕ) : ℝ :=
  ∫ q, poissonBinomial (KZA_isHermitian J h q n t a b).eigenvalues k ∂muZ

/-- **Convergence of the law of `N_A`** (paper: Theorem `thm:infvol`(e)). -/
theorem tendsto_numDistM (J h : ℝ) (n : CfgZ) (a b : ℤ) (t : ℝ) (k : ℕ) :
    Tendsto (fun M : ℕ => numDistM J h n a b t M k) atTop (𝓝 (pZ J h n a b t k)) := by
  sorry

/-- **Convergence of the number entropy** (paper: Theorem `thm:infvol`(e)):
`S_N^{(M)}(t) → H(p_t)`. -/
theorem tendsto_numberEntropyM (J h : ℝ) (n : CfgZ) (a b : ℤ) (hab : a ≤ b) (t : ℝ) :
    Tendsto (fun M : ℕ => shannon (2 * M + 1) (numDistM J h n a b t M)) atTop
      (𝓝 (shannon (b - a + 1).toNat (pZ J h n a b t))) := by
  sorry

/-- **Bounded number entropy in the thermodynamic limit** (paper: Theorem `thm:infvol`(e)):
`S_N(t) ≤ ½ ln(2πe(2Θ(F) + 4Σ'₁(F)² + 1/12))` for every `t`. -/
theorem numberEntropyZ_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (n : CfgZ) (a b : ℤ) (hab : a ≤ b)
    (t : ℝ) :
    shannon (b - a + 1).toNat (pZ J h n a b t) ≤
      (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 *
        (2 * thetaF F + 4 * (∑' d : ℕ, (d : ℝ) * √(F d)) ^ 2 + 1 / 12)) := by
  sorry

end DFL
