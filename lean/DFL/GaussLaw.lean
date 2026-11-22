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

namespace GaussLawAux

/-- Every element of `ZMod 2` is `0` or `1`. -/
lemma zmod2_cases (a : ZMod 2) : a = 0 ∨ a = 1 := by
  revert a; decide

/-- `1 + 1 = 0` in `ZMod 2`. -/
lemma zmod2_one_add_one : (1 : ZMod 2) + 1 = 0 := by decide

/-- `sgn` turns addition in `ZMod 2` into multiplication. -/
lemma sgn_add (a b : ZMod 2) : sgn (a + b) = sgn a * sgn b := by
  rcases zmod2_cases a with rfl | rfl <;>
  rcases zmod2_cases b with rfl | rfl <;> simp [sgn, zmod2_one_add_one]

/-- `sgn a ^ 2 = 1`. -/
lemma sgn_mul_self (a : ZMod 2) : sgn a * sgn a = 1 := by
  rcases zmod2_cases a with rfl | rfl <;> simp [sgn]

/-- `sgn` of a sum is the product of the signs. -/
lemma sgn_sum {ι : Type*} (s : Finset ι) (f : ι → ZMod 2) :
    sgn (∑ i ∈ s, f i) = ∏ i ∈ s, sgn (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [sgn]
  | insert a s ha ih => rw [Finset.sum_insert ha, Finset.prod_insert ha, sgn_add, ih]

/-- The Gauss operator is the diagonal matrix of the signs of the charges. -/
lemma gauss_eq (j : Fin N) : gauss j = diagonal (fun p => sgn (chargeOf p.1 p.2 j)) := by
  unfold gauss chargeOf
  congr 1
  funext p
  rw [sgn_add (prevLink p.2 j + p.2 j) (p.1 j), mul_comm]

/-- A diagonal matrix commutes with every matrix whose support only connects indices on which
the diagonal is constant. -/
lemma diagonal_comm_of_support {ι : Type*} [Fintype ι] [DecidableEq ι] (d : ι → ℝ)
    (M : Matrix ι ι ℝ) (h : ∀ a b, M a b ≠ 0 → d a = d b) :
    diagonal d * M = M * diagonal d := by
  ext a b
  rw [diagonal_mul, mul_diagonal]
  by_cases hab : M a b = 0
  · rw [hab, mul_zero, zero_mul]
  · rw [h a b hab, mul_comm]

/-- `H` conserves the Gauss charges: its matrix elements between different charges vanish. -/
lemma H_apply_ne_zero (J h : ℝ) (p p' : Cfg N × Cfg N) (hne : H J h p p' ≠ 0) :
    chargeOf p.1 p.2 = chargeOf p'.1 p'.2 := by
  have e := congrFun (congrFun (duality (N := N) J h) (dualEquiv.symm p)) (dualEquiv.symm p')
  rw [submatrix_apply, Equiv.apply_symm_apply, Equiv.apply_symm_apply, Hdual,
    blockDiagonal_apply] at e
  rw [e] at hne
  by_contra hc
  apply hne
  rw [ite_eq_right]
  simpa using hc

/-- `1 ⊗ diag(d)` is diagonal in `(n, x)`. -/
lemma one_kronecker_diagonal (d : Cfg N → ℝ) :
    (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ diagonal d = diagonal (fun p => d p.2) := by
  rw [← diagonal_one, diagonal_kronecker_diagonal]
  simp

/-- The list product of diagonal matrices is the diagonal matrix of the products. -/
lemma list_prod_diagonal {ι α : Type*} [Fintype ι] [DecidableEq ι] (l : List α)
    (f : α → ι → ℝ) :
    (l.map (fun k => diagonal (f k))).prod = diagonal (fun i => (l.map (fun k => f k i)).prod) := by
  induction l with
  | nil => simp
  | cons a l ih => rw [List.map_cons, List.prod_cons, ih, diagonal_mul_diagonal]; simp

/-- A product of Gauss operators over a list without repetitions. -/
lemma list_prod_gauss (l : List (Fin N)) (hl : l.Nodup) :
    (l.map gauss).prod = diagonal (fun p => sgn (∑ k ∈ l.toFinset, chargeOf p.1 p.2 k)) := by
  have e : l.map gauss = l.map (fun k => diagonal (fun p => sgn (chargeOf p.1 p.2 k))) := by
    congr 1
    funext k
    exact gauss_eq k
  rw [e, list_prod_diagonal]
  congr 1
  funext p
  rw [sgn_sum, List.prod_toFinset _ hl]

/-- `H` does not change the label of the boundary link `b_{N-1}`. -/
lemma H_apply_last (J h : ℝ) (hN : 0 < N) (p p' : Cfg N × Cfg N) (hne : H J h p p' ≠ 0) :
    p.2 ⟨N - 1, by omega⟩ = p'.2 ⟨N - 1, by omega⟩ := by
  obtain ⟨n', x'⟩ := p
  obtain ⟨n, x⟩ := p'
  by_contra hc
  apply hne
  have h1 : ∀ j, (hop j ⊗ₖ sz j) (n', x') (n, x) = 0 := by
    intro j
    rw [kronecker_apply]
    by_cases hj : j = ⟨N - 1, by omega⟩
    · subst hj
      have : hop (⟨N - 1, by omega⟩ : Fin N) = 0 := by
        unfold hop
        rw [dite_eq_right]
        simp only [not_lt]
        omega
      rw [this, Matrix.zero_apply, zero_mul]
    · have : sz j x' x = 0 := by
        unfold sz
        rw [ite_eq_right]
        intro heq
        apply hc
        simp only
        rw [heq, Function.update_of_ne (Ne.symm hj)]
      rw [this, mul_zero]
  have h2 : ∀ j, ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ elec j) (n', x') (n, x) = 0 := by
    intro j
    rw [kronecker_apply, elec, diagonal_apply_ne]
    · rw [mul_zero]
    · intro heq
      exact hc (by rw [heq])
  simp only [H, Matrix.sub_apply, Matrix.smul_apply, Matrix.sum_apply, h1, h2,
    Finset.sum_const_zero, smul_zero, sub_zero]

end GaussLawAux

open GaussLawAux

/-! ### (i) -/

/-- `G_j` is self-adjoint (a real symmetric matrix). -/
theorem gauss_isSymm (j : Fin N) : (gauss j).IsSymm := by
  unfold gauss
  exact isSymm_diagonal _

/-- `G_j² = 1`. -/
theorem gauss_mul_self (j : Fin N) : gauss j * gauss j = 1 := by
  rw [gauss_eq, diagonal_mul_diagonal, ← diagonal_one]
  congr 1
  funext p
  exact sgn_mul_self _

/-- The Gauss operators commute with each other. -/
theorem gauss_comm (i j : Fin N) : gauss i * gauss j = gauss j * gauss i := by
  rw [gauss_eq, gauss_eq, diagonal_mul_diagonal, diagonal_mul_diagonal]
  congr 1
  funext p
  ring

/-- The Gauss operators commute with the fermion number. -/
theorem gauss_comm_numF (j : Fin N) : gauss j * numF = numF * gauss j := by
  rw [numF, numOp, sec_diagonal_kronecker_one, gauss_eq, diagonal_mul_diagonal,
    diagonal_mul_diagonal]
  congr 1
  funext p
  ring

/-- **Gauge invariance of `H`** (paper: Lemma `lem:gauss`(i)). -/
theorem gauss_comm_H (J h : ℝ) (j : Fin N) : gauss j * H J h = H J h * gauss j := by
  rw [gauss_eq]
  exact diagonal_comm_of_support _ _ fun a b hab => by
    simp only [H_apply_ne_zero J h a b hab]

/-! ### (ii) -/

/-- **The one-dimensional `ℤ₂` Gauss law** (paper: Lemma `lem:gauss`(ii), eq. `eq:gausslaw1d`):
`σˣ_{b_j} = (-1)^{N_{≤ j}} ∏_{k ≤ j} G_k`. -/
theorem gauss_law_1d (j : Fin N) :
    sxL j = parL (univ.filter (· ≤ j)) *
      (((List.finRange N).filter (fun k => decide (k ≤ j))).map gauss).prod := by
  rw [list_prod_gauss _ ((List.nodup_finRange N).filter _), sxL, sx, one_kronecker_diagonal, parL,
    diagonal_mul_diagonal]
  congr 1
  funext p
  obtain ⟨n, x⟩ := p
  have hfs : ((List.finRange N).filter (fun k => decide (k ≤ j))).toFinset =
      univ.filter (· ≤ j) := by
    ext k
    simp
  have hx : linkOf n (chargeOf n x) j = x j := congrFun (linkOf_chargeOf n x) j
  rw [hfs]
  simp only
  rw [← hx, linkOf, Finset.sum_add_distrib, sgn_add, mul_comm]

/-- `∏_j G_j = (-1)^{N_f} σˣ_{b_{N-1}}` (paper: Lemma `lem:gauss`(ii)). -/
theorem gauss_prod (hN : 0 < N) :
    ((List.finRange N).map gauss).prod = parL univ * sxL ⟨N - 1, by omega⟩ := by
  rw [list_prod_gauss _ (List.nodup_finRange N), sxL, sx, one_kronecker_diagonal, parL,
    diagonal_mul_diagonal, List.toFinset_finRange]
  congr 1
  funext p
  obtain ⟨n, x⟩ := p
  have hx : linkOf n (chargeOf n x) ⟨N - 1, by omega⟩ = x ⟨N - 1, by omega⟩ :=
    congrFun (linkOf_chargeOf n x) _
  have hfilter : univ.filter (fun k : Fin N => k ≤ ⟨N - 1, by omega⟩) = univ := by
    ext k
    simp only [mem_filter, mem_univ, true_and, iff_true, Fin.le_def]
    omega
  simp only
  rw [← hx, linkOf, hfilter, Finset.sum_add_distrib, sgn_add, mul_left_comm, sgn_mul_self,
    mul_one]

/-- The boundary flux `σˣ_{b_{N-1}}` commutes with `H` (paper: Lemma `lem:gauss`(ii)). -/
theorem sxL_last_comm_H (J h : ℝ) (hN : 0 < N) :
    sxL (N := N) ⟨N - 1, by omega⟩ * H J h = H J h * sxL (N := N) ⟨N - 1, by omega⟩ := by
  rw [sxL, sx, one_kronecker_diagonal]
  exact diagonal_comm_of_support _ _ fun a b hab => by
    simp only [H_apply_last J h hN a b hab]

/-! ### (iii) -/

/-- `P_q` is the joint eigenprojection: `G_j P_q = sgn (q_j) P_q`. -/
theorem gauss_mul_gaussProj (j : Fin N) (q : Cfg N) :
    gauss j * gaussProj q = sgn (q j) • gaussProj q := by
  rw [gauss_eq, gaussProj, diagonal_mul_diagonal, ← diagonal_smul]
  congr 1
  funext p
  simp only [Pi.smul_apply, smul_eq_mul]
  split_ifs with hp
  · rw [hp]
  · simp

/-- The `P_q` are mutually orthogonal projections. -/
theorem gaussProj_mul (q q' : Cfg N) :
    gaussProj q * gaussProj q' = if q = q' then gaussProj q else 0 := by
  unfold gaussProj
  rw [diagonal_mul_diagonal]
  split_ifs with hqq
  · subst hqq
    congr 1
    funext p
    split_ifs <;> simp
  · rw [← diagonal_zero]
    congr 1
    funext p
    split_ifs with h1 h2
    · exact absurd (h1.symm.trans h2) hqq
    · simp
    · simp
    · simp

/-- `∑_q P_q = 1`. -/
theorem sum_gaussProj : ∑ q : Cfg N, gaussProj q = 1 := by
  ext a b
  simp only [gaussProj, Matrix.sum_apply, diagonal_apply, Matrix.one_apply]
  split_ifs with hab
  · simp [Finset.sum_ite_eq]
  · simp

/-- **Every sector is non-trivial** (paper: Lemma `lem:gauss`(iii)): `rank P_q = 2^N`. -/
theorem rank_gaussProj (q : Cfg N) : (gaussProj q).rank = 2 ^ N := by
  rw [gaussProj, rank_diagonal]
  have e : {p : Cfg N × Cfg N // (if chargeOf p.1 p.2 = q then (1 : ℝ) else 0) ≠ 0} ≃ Cfg N :=
    { toFun := fun p => p.1.1
      invFun := fun n => ⟨(n, linkOf n q), by simp [chargeOf_linkOf]⟩
      left_inv := by
        rintro ⟨⟨n, x⟩, hp⟩
        have hq : chargeOf n x = q := by
          by_contra hc
          exact hp (by simp [hc])
        subst hq
        simp [linkOf_chargeOf]
      right_inv := fun n => rfl }
  rw [Fintype.card_congr e]
  simp

/-- Each sector is invariant under `H`. -/
theorem gaussProj_comm_H (J h : ℝ) (q : Cfg N) : gaussProj q * H J h = H J h * gaussProj q := by
  unfold gaussProj
  exact diagonal_comm_of_support _ _ fun a b hab => by
    simp only [H_apply_ne_zero J h a b hab]

end DFL
