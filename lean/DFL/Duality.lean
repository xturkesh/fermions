import DFL.Basic

/-!
# The exact free-fermion duality (paper: Theorem `thm:duality`, Corollary `cor:sector`)

We prove that the 1D ℤ₂ LGT Hamiltonian `H J h` is unitarily equivalent — by a permutation of
the configuration basis `(n, x) ↔ (n, q)` — to the block-diagonal free-fermion Hamiltonian
`Hdual J h = ⊕_q (dΓ(h_q) - h ∑_j q_j)`, and that the Gauss operators become the diagonal
charge operators `τˣ_j = diag(sgn q_j)`.

Main results:
* `DFL.hop_apply`: on an open chain, nearest-neighbour fermionic hopping has no Jordan–Wigner
  sign;
* `DFL.dGamma_hq`: explicit form of `dΓ(h_q)`;
* `DFL.dualEquiv`: the bijection `(n, q) ↦ (n, linkOf n q)`;
* `DFL.duality`: `(H J h).submatrix dualEquiv dualEquiv = Hdual J h`;
* `DFL.gauss_dual`: `(gauss j).submatrix dualEquiv dualEquiv = diagonal (sgn ∘ q_j)`;
* `DFL.numOp_dual`: the matter number operators are unchanged by the duality.
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-! ### Fermionic matrix elements -/

/-- `c_x† c_x` is the diagonal occupation operator. -/
theorem cre_mul_ann_self (x : Fin N) : cre x * ann x = diagonal (fun n => occ (n x)) := by
  sorry

/-- Matrix elements of the hopping term: on an open chain nearest-neighbour hopping is sign-free,
`⟨n'| c_j† c_{j+1} + c_{j+1}† c_j |n⟩ = 1` iff `n_j ≠ n_{j+1}` and `n'` is `n` with the two
occupations exchanged. -/
theorem hop_apply (j : Fin N) (n' n : Cfg N) :
    hop j n' n =
      if hj : j.val + 1 < N then
        (if n j ≠ n ⟨j.val + 1, hj⟩ ∧ n' = n ∘ Equiv.swap j ⟨j.val + 1, hj⟩ then 1 else 0)
      else 0 := by
  sorry

/-- Explicit second quantization of `h_q`. -/
theorem dGamma_hq (J h : ℝ) (q : Cfg N) :
    dGamma (hq J h q) =
      (-J) • ∑ j, hop j + (2 * h) • diagonal (fun n => ∑ x, sgn (q x) * occ (n x)) := by
  sorry

/-! ### The duality bijection -/

theorem chargeOf_linkOf (n q : Cfg N) : chargeOf n (linkOf n q) = q := by
  sorry

theorem linkOf_chargeOf (n x : Cfg N) : linkOf n (chargeOf n x) = x := by
  sorry

/-- The duality bijection `(n, q) ↦ (n, x(n, q))` between (occupations × charges) and
(occupations × link labels). -/
def dualEquiv : Cfg N × Cfg N ≃ Cfg N × Cfg N where
  toFun a := (a.1, linkOf a.1 a.2)
  invFun b := (b.1, chargeOf b.1 b.2)
  left_inv a := by simp [chargeOf_linkOf]
  right_inv b := by simp [linkOf_chargeOf]

@[simp] theorem dualEquiv_apply (a : Cfg N × Cfg N) : dualEquiv a = (a.1, linkOf a.1 a.2) := rfl

@[simp] theorem dualEquiv_symm_apply (b : Cfg N × Cfg N) :
    dualEquiv.symm b = (b.1, chargeOf b.1 b.2) := rfl

/-- The electric term in dual variables: `σˣ_{b_{j-1}} σˣ_{b_j} = q_j (-1)^{n_j}`. -/
theorem prevLink_linkOf_add (n q : Cfg N) (j : Fin N) :
    prevLink (linkOf n q) j + linkOf n q j = q j + n j := by
  sorry

/-- Moving a particle across the bond `(j, j+1)` flips exactly the label of link `b_j`. -/
theorem linkOf_swap (n q : Cfg N) (j : Fin N) (hj : j.val + 1 < N)
    (hn : n j ≠ n ⟨j.val + 1, hj⟩) :
    linkOf (n ∘ Equiv.swap j ⟨j.val + 1, hj⟩) q =
      Function.update (linkOf n q) j (linkOf n q j + 1) := by
  sorry

/-! ### Main statements -/

/-- **Exact duality** (paper: Theorem `thm:duality`). In the basis `(n, q)` the LGT Hamiltonian is
the block-diagonal free-fermion Hamiltonian `⊕_q (dΓ(h_q) - h ∑_j q_j)`. -/
theorem duality (J h : ℝ) :
    (H J h).submatrix (dualEquiv (N := N)) dualEquiv = Hdual J h := by
  sorry

/-- The Gauss operators become the charge operators `τˣ_j = diag(sgn q_j)`. -/
theorem gauss_dual (j : Fin N) :
    (gauss j).submatrix (dualEquiv (N := N)) dualEquiv = diagonal (fun a => sgn (a.2 j)) := by
  sorry

/-- The matter number operators are unchanged by the duality. -/
theorem numOp_dual (A : Finset (Fin N)) :
    (numOp A ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)).submatrix (dualEquiv (N := N)) dualEquiv =
      numOp A ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) := by
  sorry

/-- The spectral projections of the matter number operators are unchanged by the duality. -/
theorem numProj_dual (A : Finset (Fin N)) (k : ℕ) :
    (numProj A k ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)).submatrix (dualEquiv (N := N)) dualEquiv =
      numProj A k ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) := by
  sorry

end DFL
