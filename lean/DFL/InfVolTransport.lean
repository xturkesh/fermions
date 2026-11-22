import DFL.InfVol

/-!
# Infinite volume: the transported charge (paper: Theorem `thm:infvol`(b))

For an occupation configuration `n ∈ {0,1}^ℤ` (occupied set `S`, possibly infinite) and a cut at
`m ∈ ℤ`, the transported charge of the gauge theory on `Λ_M = [-M, M]`,
`⟨Ψ^{(M)}(t), N_{[-M, m]} Ψ^{(M)}(t)⟩ - |S ∩ [-M, m]|`, converges as `M → ∞`, and under Hypothesis (EC)
the limit is bounded uniformly in time: in finite volume it is the difference of the charge flowing
across the cut in the two directions, `∑_{z ∈ S, z > m} ∑_{x ≤ m} E|U_{xz}|² - ∑_{z ∈ S, z ≤ m} ∑_{x > m} E|U_{xz}|²`,
each term converges by `tendsto_evoZ`, and the double sums are dominated by `∑_d d F(d) < ∞`.
-/

open MeasureTheory Filter Topology Matrix Finset

noncomputable section

namespace DFL

/-- The occupied sites of `n` inside `Λ_M`, relabelled in `Fin (2M+1)`. -/
def occM (M : ℕ) (n : CfgZ) : Finset (Fin (2 * M + 1)) :=
  univ.filter (fun i => n ((i.val : ℤ) - M) = 1)

/-- The sites `[-M, m]` of `Λ_M`, relabelled in `Fin (2M+1)`. -/
def leftM (M : ℕ) (m : ℤ) : Finset (Fin (2 * M + 1)) :=
  univ.filter (fun i => (i.val : ℤ) - M ≤ m)

/-- The transported charge on `Λ_M`: `⟨Ψ^{(M)}(t), N_{[-M,m]} Ψ^{(M)}(t)⟩ - |S ∩ [-M, m]|`, for the gauge
theory with Fock matter `n|_{Λ_M}` and `z`-polarized links. -/
def transportM (J h : ℝ) (n : CfgZ) (m : ℤ) (t : ℝ) (M : ℕ) : ℝ :=
  (expect (matter (numOp (leftM M m))) (lgtState J h (fockVec (occM M n)) t)).re -
    ((occM M n ∩ leftM M m).card : ℝ)

/-- The finite-volume transported charge is bounded uniformly in `M` and `t` under (EC). -/
theorem abs_transportM_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (n : CfgZ) (m : ℤ) (t : ℝ)
    (M : ℕ) : |transportM J h n m t M| ≤ 2 * ∑' d : ℕ, (d : ℝ) * F d := by
  sorry

/-- **Transported charge in the thermodynamic limit** (paper: Theorem `thm:infvol`(b)): the limit
`T_m(t) = lim_{M → ∞} transportM` exists, and under (EC) `|T_m(t)| ≤ 2 ∑_d d F(d)` for every `t`. -/
theorem transport_tendsto (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (n : CfgZ) (m : ℤ) (t : ℝ) :
    ∃ T : ℝ, Tendsto (fun M : ℕ => transportM J h n m t M) atTop (𝓝 T) ∧
      |T| ≤ 2 * ∑' d : ℕ, (d : ℝ) * F d := by
  sorry

end DFL
