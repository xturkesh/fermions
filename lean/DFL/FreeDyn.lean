import DFL.Basic
import DFL.Evolution

/-!
# Free-fermion dynamics

Canonical anticommutation relations, the Heisenberg evolution of the annihilation operators under
a quadratic Hamiltonian `dΓ(A)`, and the density and number-variance formulas for Fock initial
states (Wick's theorem for Fock states).

Main results:
* `DFL.ann_mul_ann_add`, `DFL.ann_mul_cre_add`: CAR;
* `DFL.dGamma_comm_ann`: `[dΓ(A), c_x] = -∑_y A_{xy} c_y`;
* `DFL.heisenberg_ann`: `W(t)† c_x W(t) = ∑_y (e^{-itA})_{xy} c_y`, with `W(t) = e^{-itdΓ(A)}`;
* `DFL.two_point`: `⟨c_x† c_y⟩_t = C(t)_{yx}` with `C(t) = U(t) P_S U(t)†`;
* `DFL.density`: `⟨n_x⟩_t = ∑_{z ∈ S} |U(t)_{xz}|²`;
* `DFL.number_variance`: `Var N_R = ∑_{x ∈ R, y ∉ R} |C(t)_{xy}|²`.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-! ### Canonical anticommutation relations -/

theorem ann_mul_ann_add (x y : Fin N) : ann x * ann y + ann y * ann x = 0 := by
  sorry

theorem ann_mul_cre_add (x y : Fin N) :
    ann x * cre y + cre y * ann x = if x = y then 1 else 0 := by
  sorry

/-- `[dΓ(A), c_x] = -∑_y A_{xy} c_y`. -/
theorem dGamma_comm_ann (A : Matrix (Fin N) (Fin N) ℝ) (x : Fin N) :
    dGamma A * ann x - ann x * dGamma A = -∑ y, A x y • ann y := by
  sorry

/-- `dΓ` of a symmetric matrix is symmetric. -/
theorem dGamma_isSymm (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) : (dGamma A).IsSymm := by
  sorry

/-! ### Heisenberg evolution -/

/-- **Heisenberg evolution of the annihilation operators** under a quadratic Hamiltonian:
`e^{itdΓ(A)} c_x e^{-itdΓ(A)} = ∑_y (e^{-itA})_{xy} c_y`. -/
theorem heisenberg_ann (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (x : Fin N) (t : ℝ) :
    (evo (dGamma A) t)ᴴ * cplx (ann x) * evo (dGamma A) t = ∑ y, evo A t x y • cplx (ann y) := by
  sorry

/-! ### Fock initial states -/

/-- The Fock basis vector with occupied set `S`. -/
def fockVec (S : Finset (Fin N)) : Cfg N → ℂ := Pi.single (fockCfg S) 1

/-- The one-body correlation matrix `C(t) = U(t) P_S U(t)†` of the evolved Fock state, where
`U(t) = e^{-itA}` and `P_S` is the projection onto `S`. -/
def corrMat (A : Matrix (Fin N) (Fin N) ℝ) (S : Finset (Fin N)) (t : ℝ) :
    Matrix (Fin N) (Fin N) ℂ :=
  evo A t * diagonal (fun z => if z ∈ S then (1 : ℂ) else 0) * (evo A t)ᴴ

/-- The two-point function of an evolved Fock state: `⟨c_x† c_y⟩_t = C(t)_{yx}`. -/
theorem two_point (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (x y : Fin N) (t : ℝ) :
    expect (cplx (cre x * ann y)) (evo (dGamma A) t *ᵥ fockVec S) = corrMat A S t y x := by
  sorry

/-- The density of an evolved Fock state: `⟨n_x⟩_t = ∑_{z ∈ S} |U(t)_{xz}|²`. -/
theorem density (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N)) (x : Fin N)
    (t : ℝ) :
    expect (cplx (diagonal fun n => occ (n x))) (evo (dGamma A) t *ᵥ fockVec S) =
      ((∑ z ∈ S, Complex.normSq (evo A t x z) : ℝ) : ℂ) := by
  sorry

/-- The variance of the particle number of a region `R` in an evolved Fock state:
`⟨N_R²⟩ - ⟨N_R⟩² = ∑_{x ∈ R} ∑_{y ∉ R} |C(t)_{xy}|²`. -/
theorem number_variance (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S R : Finset (Fin N))
    (t : ℝ) :
    expect (cplx (numOp R * numOp R)) (evo (dGamma A) t *ᵥ fockVec S) -
        expect (cplx (numOp R)) (evo (dGamma A) t *ᵥ fockVec S) ^ 2 =
      ((∑ x ∈ R, ∑ y ∈ univ \ R, Complex.normSq (corrMat A S t x y) : ℝ) : ℂ) := by
  sorry

/-! ### The distribution of the particle number -/

/-- Distribution of the particle number `N_R` in the state `v`: `p(k) = ∑_{n : N_R(n) = k} |v(n)|²`. -/
def numDist (R : Finset (Fin N)) (v : Cfg N → ℂ) (k : ℕ) : ℝ :=
  ∑ n ∈ univ.filter (fun n => count R n = k), Complex.normSq (v n)

/-- The distribution `numDist` is the spectral measure of `N_R`: its moments are expectations. -/
theorem numDist_eq_expect (R : Finset (Fin N)) (v : Cfg N → ℂ) (k : ℕ) :
    ((numDist R v k : ℝ) : ℂ) = expect (cplx (numProj R k)) v := by
  sorry

end DFL
