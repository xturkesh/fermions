import DFL.LRFull
import DFL.Imbalance

/-!
# The operator norm sees the clean sector (paper: Proposition `prop:LR`(c))

* `‖[n_x(t), n_y]‖ = max_q |U_q(t)_{xy}| (1 - |U_q(t)_{xy}|²)^{1/2}` (`norm_comm_nn_eq`): in the sector `q` the
  commutator is `dΓ(Z)` with `Z = [|w⟩⟨w|, |δ_y⟩⟨δ_y|]`, `w = U_q(t)^* δ_x`, whose eigenvalues are `±iμ, 0`,
  `μ = |U_{xy}| (1 - |U_{xy}|²)^{1/2}`;
* in the clean sector `q = (1, …, 1)` the eigenvectors of `h_q = -JΔ + 2h` are the sine modes, and
  `sup_t |U(t)_{1N}|² ≥ ∑_k φ_k(1)⁴ = 3/(2(N+1))` for `N ≥ 2`, whence
  `sup_t ‖{c̃_1(t), c̃_N†}‖ ≥ (3/(2(N+1)))^{1/2}` and `sup_t ‖[n_1(t), n_N]‖ ≥ (3/(4(N+1)))^{1/2}`
  (`sup_norm_dressedAnticomm_ge`, `sup_norm_comm_nn_ge`), while the sector averages are at most `F(N-1)`.

Sites are 0-indexed: the paper's sites `1` and `N` are `0` and `N - 1`.
-/

open Matrix Finset
open scoped Kronecker Matrix.Norms.L2Operator

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The commutator `[n_x(t), n_y]`. -/
def numComm (J h : ℝ) (x y : Fin N) (t : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  heisT J h t (cplx (numOp {x} ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))) *
      cplx (numOp {y} ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) -
    cplx (numOp {y} ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) *
      heisT J h t (cplx (numOp {x} ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)))

/-- **`‖[n_x(t), n_y]‖ = max_q |U_q(t)_{xy}| (1 - |U_q(t)_{xy}|²)^{1/2}`** (paper: eq. `eq:main-LRnorm`). -/
theorem norm_comm_nn_eq (J h : ℝ) (x y : Fin N) (t : ℝ) :
    ‖numComm J h x y t‖ =
      univ.sup' univ_nonempty (fun q : Cfg N =>
        ‖evo (hq J h q) t x y‖ * √(1 - ‖evo (hq J h q) t x y‖ ^ 2)) := by
  sorry

/-- **The clean sector forbids a uniform operator-norm bound** (paper: Proposition `prop:LR`(c)):
`sup_t ‖{c̃_1(t), c̃_N†}‖ ≥ (3/(2(N+1)))^{1/2}` for `N ≥ 2`. -/
theorem sup_norm_dressedAnticomm_ge (J h : ℝ) (hJ : J ≠ 0) (hN : 2 ≤ N) :
    √(3 / (2 * ((N : ℝ) + 1))) ≤
      ⨆ t : ℝ, ‖dressedAnticomm (N := N) J h ⟨0, by omega⟩ ⟨N - 1, by omega⟩ t‖ := by
  sorry

/-- `sup_t ‖[n_1(t), n_N]‖ ≥ (3/(4(N+1)))^{1/2}` for `N ≥ 2` (paper: Proposition `prop:LR`(c)). -/
theorem sup_norm_comm_nn_ge (J h : ℝ) (hJ : J ≠ 0) (hN : 2 ≤ N) :
    √(3 / (4 * ((N : ℝ) + 1))) ≤ ⨆ t : ℝ, ‖numComm (N := N) J h ⟨0, by omega⟩ ⟨N - 1, by omega⟩ t‖ := by
  sorry

/-- The sector average is small: `E_q sup_t ‖P_q {c̃_1(t), c̃_N†} P_q‖ ≤ F(N-1)` (paper: Proposition
`prop:LR`(c)). -/
theorem avg_sup_anticomm_corner_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (hN : 1 ≤ N) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ, ‖cplx (gaussProj q) *
      dressedAnticomm J h ⟨0, by omega⟩ ⟨N - 1, by omega⟩ t * cplx (gaussProj q)‖) ≤ F (N - 1) := by
  sorry

end DFL
