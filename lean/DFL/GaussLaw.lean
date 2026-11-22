import DFL.Sector

/-!
# The Gauss law (paper: Lemma `lem:gauss`)

* (i) The Gauss operators `G_j` are self-adjoint involutions; they commute with each other, with the
  fermion number `N_f` and with `H` (`gauss_isSymm`, `gauss_mul_self`, `gauss_comm`,
  `gauss_comm_numF`, `gauss_comm_H`).
* (ii) The one-dimensional `ℤ₂` Gauss law `σˣ_{b_j} = (-1)^{N_{≤ j}} ∏_{k ≤ j} G_k`
  (`gauss_law_1d`); in particular `∏_j G_j = (-1)^{N_f} σˣ_{b_N}` (`gauss_prod`), and the boundary flux
  `σˣ_{b_N}` commutes with `H` (`sxL_last_comm_H`).
* (iii) The joint spectral projections `P_q` of the Gauss operators (`gaussProj`) have rank `2^N`,
  are mutually orthogonal, sum to `1` and commute with `H`: `2^N` non-trivial superselection sectors,
  each invariant under `H`.

The link of site `j` to its right is `b_j` (0-indexed sites and links); `b_{N-1}` is the boundary link.
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The total fermion number `N_f = ∑_x n_x` on the LGT Hilbert space. -/
def numF : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  numOp univ ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)

/-- `σˣ` on the link `b_j`, as an operator on the LGT Hilbert space. -/
def sxL (j : Fin N) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ sx j

/-- The parity `(-1)^{N_T}` of the fermions in `T`, on the LGT Hilbert space. -/
def parL (T : Finset (Fin N)) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  diagonal (fun p => sgn (∑ k ∈ T, p.1 k))

/-- The joint spectral projection `P_q` of the Gauss operators onto `{G_j = sgn (q_j) ∀ j}`. -/
def gaussProj (q : Cfg N) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  diagonal (fun p => if chargeOf p.1 p.2 = q then 1 else 0)

/-! ### (i) -/

/-- `G_j` is self-adjoint (a real symmetric matrix). -/
theorem gauss_isSymm (j : Fin N) : (gauss j).IsSymm := by
  sorry

/-- `G_j² = 1`. -/
theorem gauss_mul_self (j : Fin N) : gauss j * gauss j = 1 := by
  sorry

/-- The Gauss operators commute with each other. -/
theorem gauss_comm (i j : Fin N) : gauss i * gauss j = gauss j * gauss i := by
  sorry

/-- The Gauss operators commute with the fermion number. -/
theorem gauss_comm_numF (j : Fin N) : gauss j * numF = numF * gauss j := by
  sorry

/-- **Gauge invariance of `H`** (paper: Lemma `lem:gauss`(i)). -/
theorem gauss_comm_H (J h : ℝ) (j : Fin N) : gauss j * H J h = H J h * gauss j := by
  sorry

/-! ### (ii) -/

/-- **The one-dimensional `ℤ₂` Gauss law** (paper: Lemma `lem:gauss`(ii), eq. `eq:gausslaw1d`):
`σˣ_{b_j} = (-1)^{N_{≤ j}} ∏_{k ≤ j} G_k`. -/
theorem gauss_law_1d (j : Fin N) :
    sxL j = parL (univ.filter (· ≤ j)) *
      (((List.finRange N).filter (fun k => decide (k ≤ j))).map gauss).prod := by
  sorry

/-- `∏_j G_j = (-1)^{N_f} σˣ_{b_{N-1}}` (paper: Lemma `lem:gauss`(ii)). -/
theorem gauss_prod (hN : 0 < N) :
    ((List.finRange N).map gauss).prod = parL univ * sxL ⟨N - 1, by omega⟩ := by
  sorry

/-- The boundary flux `σˣ_{b_{N-1}}` commutes with `H` (paper: Lemma `lem:gauss`(ii)). -/
theorem sxL_last_comm_H (J h : ℝ) (hN : 0 < N) :
    sxL (N := N) ⟨N - 1, by omega⟩ * H J h = H J h * sxL (N := N) ⟨N - 1, by omega⟩ := by
  sorry

/-! ### (iii) -/

/-- `P_q` is the joint eigenprojection: `G_j P_q = sgn (q_j) P_q`. -/
theorem gauss_mul_gaussProj (j : Fin N) (q : Cfg N) :
    gauss j * gaussProj q = sgn (q j) • gaussProj q := by
  sorry

/-- The `P_q` are mutually orthogonal projections. -/
theorem gaussProj_mul (q q' : Cfg N) :
    gaussProj q * gaussProj q' = if q = q' then gaussProj q else 0 := by
  sorry

/-- `∑_q P_q = 1`. -/
theorem sum_gaussProj : ∑ q : Cfg N, gaussProj q = 1 := by
  sorry

/-- **Every sector is non-trivial** (paper: Lemma `lem:gauss`(iii)): `rank P_q = 2^N`. -/
theorem rank_gaussProj (q : Cfg N) : (gaussProj q).rank = 2 ^ N := by
  sorry

/-- Each sector is invariant under `H`. -/
theorem gaussProj_comm_H (J h : ℝ) (q : Cfg N) : gaussProj q * H J h = H J h * gaussProj q := by
  sorry

end DFL
