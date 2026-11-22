import DFL.TwoD
import DFL.Evolution

/-!
# Checkerboard propagator and the diagonal of `K_s²` (paper: Proposition `prop:2dlaws`(c)(iii),(iv))

For a bipartite graph with sublattice signs `ε` and the checkerboard charges `q⁰ = -ε`, the sector
operator is `h = h_{q⁰,s} = K_s - 2h Γ_ε` with `h² = K_s² + 4h²` (`FinGraph.checkerboard_sq`). With
`R = (K_s² + 4h²)^{1/2}`:

* (iii) `e^{-ith} = cos(tR) - i h R⁻¹ sin(tR)`: on every eigenvector `v` of `K_s²` with eigenvalue `μ`,
  `e^{-ith} v = cos(tω) v - i ω⁻¹ sin(tω) h v`, `ω = √(μ + 4h²)` (`checkerboard_evo_eigvec`); as an
  operator identity, with `cos(tR) = ∑_k (-1)^k t^{2k} (R²)^k / (2k)!` and
  `R⁻¹ sin(tR) = ∑_k (-1)^k t^{2k+1} (R²)^k / (2k+1)!` (`checkerboard_evo`);
* (iv) on a simple graph and for signs `s ∈ {±1}^E`, `(K_s²)_{ii} = J² deg(i)` (`Ks_sq_diag`).
-/

open Matrix Finset

noncomputable section

namespace DFL

namespace FinGraph

variable {N M : ℕ} (G : FinGraph N M)

/-- **Checkerboard propagator on eigenvectors** (paper: Proposition `prop:2dlaws`(c)(iii)): if
`K_s² v = μ v`, then `e^{-ith} v = cos(tω) v - i ω⁻¹ sin(tω) h v` with `ω = √(μ + 4h²)`. -/
theorem checkerboard_evo_eigvec (J h : ℝ) (hh : h ≠ 0) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε1 : ∀ i, ε i = 1 ∨ ε i = -1) (hε : ∀ e, ε (G.src e) = -ε (G.tgt e))
    (t μ : ℝ) (v : Fin N → ℂ) (hv : cplx (G.Ks J s * G.Ks J s) *ᵥ v = (μ : ℂ) • v) :
    evo (G.hqs J h (fun i => -ε i) s) t *ᵥ v =
      ((Real.cos (t * √(μ + 4 * h ^ 2)) : ℝ) : ℂ) • v -
        (Complex.I * ((Real.sin (t * √(μ + 4 * h ^ 2)) / √(μ + 4 * h ^ 2) : ℝ) : ℂ)) •
          (cplx (G.hqs J h (fun i => -ε i) s) *ᵥ v) := by
  sorry

/-- **Checkerboard propagator** (paper: Proposition `prop:2dlaws`(c)(iii)):
`e^{-ith} = cos(tR) - i h R⁻¹ sin(tR)`, `R² = K_s² + 4h²`, with `cos(tR)` and `R⁻¹ sin(tR)` given by
their (even) power series in `R²`. -/
theorem checkerboard_evo (J h : ℝ) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε1 : ∀ i, ε i = 1 ∨ ε i = -1) (hε : ∀ e, ε (G.src e) = -ε (G.tgt e)) (t : ℝ) :
    evo (G.hqs J h (fun i => -ε i) s) t =
      (∑' k : ℕ, (((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k) / ((2 * k).factorial : ℂ))) •
          cplx (G.Ks J s * G.Ks J s + (4 * h ^ 2) • (1 : Matrix (Fin N) (Fin N) ℝ)) ^ k) -
        Complex.I • (cplx (G.hqs J h (fun i => -ε i) s) *
          ∑' k : ℕ, (((-1 : ℂ) ^ k * (t : ℂ) ^ (2 * k + 1) / ((2 * k + 1).factorial : ℂ))) •
            cplx (G.Ks J s * G.Ks J s + (4 * h ^ 2) • (1 : Matrix (Fin N) (Fin N) ℝ)) ^ k) := by
  sorry

/-- **Diagonal of `K_s²`** (paper: Proposition `prop:2dlaws`(c)(iv)): on a simple graph (no parallel
edges) and for signs `s_e ∈ {±1}`, `(K_s²)_{ii} = J² deg(i)`. -/
theorem Ks_sq_diag (J : ℝ) (s : Fin M → ℝ) (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hsimple : ∀ e e', ((G.src e = G.src e' ∧ G.tgt e = G.tgt e') ∨
      (G.src e = G.tgt e' ∧ G.tgt e = G.src e')) → e = e')
    (i : Fin N) :
    (G.Ks J s * G.Ks J s) i i =
      J ^ 2 * ((univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card : ℝ) := by
  sorry

end FinGraph

end DFL
