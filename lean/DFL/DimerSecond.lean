import DFL.DimerInfVol

/-!
# Second moment of the dimer transfer matrices near the critical energy (paper: Lemma `lem:dimersecond`)

Let `0 < λ < 1` and `D_±(z) = T_{±λ - z}²` (`dimerDp`, `dimerDm`). In the basis `S` of `lem:dimercritical`(iii)
(`dimerS`), `M_±(z) = S⁻¹ D_±(z) S` are, at `z = λ`, the rotations `-1` and `R_{2θ}`. The averaged map
`𝓛_z(X) = ½ ∑_± M_±(z)* X M_±(z)` on Hermitian `2 × 2` matrices, in the coordinates
`X = x 1 + y σ_y + w₁ σ_z + w₂ σ_x`, is the identity on `(x, y)` and the contraction `½(1 + R_{4θ})` (norm
`|cos 2θ| = |2λ² - 1| < 1`) on `w` at `z = λ`; for real `z` the `y`-coordinate decouples (the `M_±` are real
symplectic), and the first-order change of the `(x, x)`-entry vanishes (`det M_± ≡ 1`). In the weighted norm
`|x| + |y| + c|w|`, `c ≍ |z - λ|`, this gives `‖𝓛_z‖ ≤ 1 + K((Re z - λ)² + |Im z|)`, hence
`E ‖D_{ε_n}(z) ⋯ D_{ε_1}(z)‖₂² ≤ C e^{K n ((Re z - λ)² + |Im z|)}` (`dimer_second_moment`).

`integral_hs2_tmProd_le` is the consequence for the transfer matrices `T^z(m + L, m)` (`tmProd`) of the random dimer
potential `v^{(r)}` over any window: complete dimers are i.i.d., and at most two further sites contribute a bounded
factor.
-/

open Matrix Finset MeasureTheory

noncomputable section

namespace DFL

/-- The squared Hilbert–Schmidt norm of a `2 × 2` complex matrix. -/
def hs2 (A : Matrix (Fin 2) (Fin 2) ℂ) : ℝ := ∑ i, ∑ j, ‖A i j‖ ^ 2

/-- **Second moment near the critical energy** (paper: Lemma `lem:dimersecond`): the average over the `2^n` dimer
words of `‖D_{ε_1}(z) ⋯ D_{ε_n}(z)‖₂²` is at most `C e^{K n ((Re z - λ)² + |Im z|)}` for `|z - λ| ≤ δ₀`. -/
theorem dimer_second_moment (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) :
    ∃ C K δ₀ : ℝ, 0 < C ∧ 0 ≤ K ∧ 0 < δ₀ ∧ ∀ z : ℂ, ‖z - (lam : ℂ)‖ ≤ δ₀ → ∀ n : ℕ,
      (∑ s : Fin n → Bool,
          hs2 (List.ofFn fun i => if s i then dimerDp lam z else dimerDm lam z).prod) / 2 ^ n ≤
        C * Real.exp (K * n * ((z.re - lam) ^ 2 + |z.im|)) := by
  sorry

/-- **Second moment over windows of the random dimer potential**: for `λ = 2|h|/J ∈ (0, 1)`, every alignment `r`,
every window `[m, m + L - 1]` and `z` near `λ`, `E ‖T^z(m + L, m)‖₂² ≤ C e^{K L ((Re z - λ)² + |Im z|)}`. -/
theorem integral_hs2_tmProd_le (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (r : ℤ) :
    ∃ C K δ₀ : ℝ, 0 < C ∧ 0 ≤ K ∧ 0 < δ₀ ∧ ∀ z : ℂ, ‖z - ((2 * |h| / J : ℝ) : ℂ)‖ ≤ δ₀ →
      ∀ (m : ℤ) (L : ℕ),
        Integrable (fun η => hs2 (tmProd (dimerPotZ J h r η) z m L)) muZ ∧
        ∫ η, hs2 (tmProd (dimerPotZ J h r η) z m L) ∂muZ ≤
          C * Real.exp (K * L * ((z.re - 2 * |h| / J) ^ 2 + |z.im|)) := by
  sorry

end DFL
