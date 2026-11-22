import DFL.OneBody

/-!
# Transfer matrices and Dirichlet solutions of the one-dimensional Schrödinger equation
(paper: Section `sec:onebody:EC`, Lemma `lem:green` and the deterministic step of the proof of
Theorem `thm:EC`)

Deterministic ingredients of the proof of Hypothesis (EC), for an arbitrary real potential
`V : ℤ → ℝ` and the equation `u(z+1) + u(z-1) + V z u(z) = E u(z)` (the operator `Δ + diag V`).

* `dirL V E a k = u₋(a - 1 + k)` is the solution with `u₋(a-1) = 0`, `u₋(a) = 1`; for the box
  `B = [a, a+L-1]`, `dirL V E a (L+1) = det(E - A_B)` (`P_B(E)`), and `dirL V E a (i+1)` is the
  determinant of `[a, x-1]`, `x = a + i` (`p_b(E)`).
* `dirR V E b k = u₊(b + 1 - k)` is the solution with `u₊(b+1) = 0`, `u₊(b) = 1`;
  `dirR V E b (L-i)` is the determinant of `[x+1, b]` (`p_a(E)`).
* `poisson`: for every solution `φ` on `B`, `P_B(E) φ(x) = p_a(E) φ(a-1) + p_b(E) φ(b+1)` (this
  replaces the Green-function formula of Lemma `lem:green`(b) and Eq. `eq:poisson`).
* `box_right`, `box_left`: `P_B(E)` through one product of transfer matrices applied to a vector
  determined by the potential on the other side (Lemma `lem:green`(c),(d)).
* `sl2_entry_ineq`: the deterministic inequality behind Lemma `lem:entry`.
* `eigvec_recursion`, `abs_eigenvalue_hq_le`: eigenvectors of `h_q` solve the equation with
  `V = -(2h/J) sgn q`, `E = -E'/J`, and `|E'| ≤ 2|J| + 2|h|`.
-/

open Matrix Finset Polynomial

noncomputable section

namespace DFL

/-- One-step transfer matrix of `u(z+1) + u(z-1) + v u(z) = E u(z)`:
`(u(z+1), u(z)) = T_E(v) (u(z), u(z-1))`. -/
def transferT (E v : ℝ) : Matrix (Fin 2) (Fin 2) ℝ := !![E - v, -1; 1, 0]

/-- The Euclidean norm on `ℝ²`. -/
def enorm2 (u : Fin 2 → ℝ) : ℝ := Real.sqrt (u 0 ^ 2 + u 1 ^ 2)

/-- Left Dirichlet solution: `dirL V E a k = u(a - 1 + k)`, where `u(a-1) = 0`, `u(a) = 1` and
`u(z+1) = (E - V z) u(z) - u(z-1)`. -/
def dirL (V : ℤ → ℝ) (E : ℝ) (a : ℤ) : ℕ → ℝ
  | 0 => 0
  | 1 => 1
  | (k + 2) => (E - V (a + k)) * dirL V E a (k + 1) - dirL V E a k

/-- Right Dirichlet solution: `dirR V E b k = u(b + 1 - k)`, where `u(b+1) = 0`, `u(b) = 1` and
`u(z-1) = (E - V z) u(z) - u(z+1)`. -/
def dirR (V : ℤ → ℝ) (E : ℝ) (b : ℤ) : ℕ → ℝ
  | 0 => 0
  | 1 => 1
  | (k + 2) => (E - V (b - k)) * dirR V E b (k + 1) - dirR V E b k

/-- `dirL` as a polynomial in `E`. -/
def dirLpoly (V : ℤ → ℝ) (a : ℤ) : ℕ → ℝ[X]
  | 0 => 0
  | 1 => 1
  | (k + 2) => (X - C (V (a + k))) * dirLpoly V a (k + 1) - dirLpoly V a k

/-- `dirR` as a polynomial in `E`. -/
def dirRpoly (V : ℤ → ℝ) (b : ℤ) : ℕ → ℝ[X]
  | 0 => 0
  | 1 => 1
  | (k + 2) => (X - C (V (b - k))) * dirRpoly V b (k + 1) - dirRpoly V b k

/-- `T_E(V(s+m-1)) ⋯ T_E(V(s))`: the product over the sites `s, …, s+m-1`, later sites to the left. -/
def tprodL (V : ℤ → ℝ) (E : ℝ) (s : ℤ) : ℕ → Matrix (Fin 2) (Fin 2) ℝ
  | 0 => 1
  | (m + 1) => transferT E (V (s + m)) * tprodL V E s m

/-- `T_E(V(s)) ⋯ T_E(V(s+m-1))`: the product over the sites `s, …, s+m-1`, earlier sites to the left. -/
def tprodR (V : ℤ → ℝ) (E : ℝ) (s : ℤ) : ℕ → Matrix (Fin 2) (Fin 2) ℝ
  | 0 => 1
  | (m + 1) => tprodR V E s m * transferT E (V (s + m))

/-! ### Polynomials -/

theorem eval_dirLpoly (V : ℤ → ℝ) (E : ℝ) (a : ℤ) (k : ℕ) :
    (dirLpoly V a k).eval E = dirL V E a k := by
  sorry

theorem eval_dirRpoly (V : ℤ → ℝ) (E : ℝ) (b : ℤ) (k : ℕ) :
    (dirRpoly V b k).eval E = dirR V E b k := by
  sorry

theorem dirLpoly_monic (V : ℤ → ℝ) (a : ℤ) (k : ℕ) : (dirLpoly V a (k + 1)).Monic := by
  sorry

theorem natDegree_dirLpoly (V : ℤ → ℝ) (a : ℤ) (k : ℕ) :
    (dirLpoly V a (k + 1)).natDegree = k := by
  sorry

theorem dirRpoly_monic (V : ℤ → ℝ) (b : ℤ) (k : ℕ) : (dirRpoly V b (k + 1)).Monic := by
  sorry

theorem natDegree_dirRpoly (V : ℤ → ℝ) (b : ℤ) (k : ℕ) :
    (dirRpoly V b (k + 1)).natDegree = k := by
  sorry

/-! ### Locality -/

/-- `dirL V E a k` depends only on `V (a), …, V (a + k - 2)`. -/
theorem dirL_congr (V V' : ℤ → ℝ) (E : ℝ) (a : ℤ) (k : ℕ)
    (h : ∀ j : ℕ, j + 1 < k → V (a + j) = V' (a + j)) : dirL V E a k = dirL V' E a k := by
  sorry

/-- `dirR V E b k` depends only on `V (b), …, V (b - k + 2)`. -/
theorem dirR_congr (V V' : ℤ → ℝ) (E : ℝ) (b : ℤ) (k : ℕ)
    (h : ∀ j : ℕ, j + 1 < k → V (b - j) = V' (b - j)) : dirR V E b k = dirR V' E b k := by
  sorry

theorem dirLpoly_congr (V V' : ℤ → ℝ) (a : ℤ) (k : ℕ)
    (h : ∀ j : ℕ, j + 1 < k → V (a + j) = V' (a + j)) : dirLpoly V a k = dirLpoly V' a k := by
  sorry

theorem dirRpoly_congr (V V' : ℤ → ℝ) (b : ℤ) (k : ℕ)
    (h : ∀ j : ℕ, j + 1 < k → V (b - j) = V' (b - j)) : dirRpoly V b k = dirRpoly V' b k := by
  sorry

/-! ### Identities -/

/-- Consecutive values of a Dirichlet solution do not both vanish. -/
theorem dirL_pair_ne_zero (V : ℤ → ℝ) (E : ℝ) (a : ℤ) (k : ℕ) :
    dirL V E a (k + 1) ≠ 0 ∨ dirL V E a k ≠ 0 := by
  sorry

/-- Consecutive values of a Dirichlet solution do not both vanish. -/
theorem dirR_pair_ne_zero (V : ℤ → ℝ) (E : ℝ) (b : ℤ) (k : ℕ) :
    dirR V E b (k + 1) ≠ 0 ∨ dirR V E b k ≠ 0 := by
  sorry

/-- Both Dirichlet solutions give the determinant of the box `[a, a+L-1]` (constancy of the
Wronskian; paper: Lemma `lem:green`(a)). -/
theorem dirL_eq_dirR (V : ℤ → ℝ) (E : ℝ) (a : ℤ) (L : ℕ) :
    dirL V E a (L + 1) = dirR V E (a + L - 1) (L + 1) := by
  sorry

/-- **Poisson identity**: if `φ` solves `φ(z+1) + φ(z-1) + V z φ(z) = E φ(z)` on the box
`[a, a+L-1]`, then at `x = a + i`, `P_B(E) φ(x) = p_a(E) φ(a-1) + p_b(E) φ(a+L)`. -/
theorem poisson (V : ℤ → ℝ) (E : ℝ) (φ : ℤ → ℝ) (a : ℤ) (L i : ℕ) (hi : i < L)
    (hφ : ∀ z, a ≤ z → z ≤ a + L - 1 → φ (z + 1) + φ (z - 1) + V z * φ z = E * φ z) :
    dirL V E a (L + 1) * φ (a + i) =
      dirR V E (a + L - 1) (L - i) * φ (a - 1) + dirL V E a (i + 1) * φ (a + L) := by
  sorry

/-- **Right end** (paper: Lemma `lem:green`(c)): `P_B(E) = (T_b ⋯ T_{x+1} Y)₀` with
`Y = (u₋(x+1), u₋(x))`, which depends only on `V` on `[a, x]`; note `p_b(E) = Y₁`. -/
theorem box_right (V : ℤ → ℝ) (E : ℝ) (a : ℤ) (L i : ℕ) (hi : i < L) :
    dirL V E a (L + 1) =
      (tprodL V E (a + i + 1) (L - 1 - i) *ᵥ ![dirL V E a (i + 2), dirL V E a (i + 1)]) 0 := by
  sorry

/-- **Left end** (paper: Lemma `lem:green`(d)): `P_B(E) = (T_a ⋯ T_{x-1} Y')₀` with
`Y' = (u₊(x-1), u₊(x))`, which depends only on `V` on `[x, b]`; note `p_a(E) = Y'₁`. -/
theorem box_left (V : ℤ → ℝ) (E : ℝ) (a : ℤ) (L i : ℕ) (hi : i < L) :
    dirL V E a (L + 1) =
      (tprodR V E a i *ᵥ
        ![dirR V E (a + L - 1) (L - i + 1), dirR V E (a + L - 1) (L - i)]) 0 := by
  sorry

/-- Transfer matrices have determinant one. -/
theorem det_transferT (E v : ℝ) : (transferT E v).det = 1 := by
  sorry

/-- `T_E(v)ᵀ = D T_E(v) D`, `D = diag(1, -1)`. -/
theorem transferT_transpose (E v : ℝ) :
    (transferT E v)ᵀ = !![1, 0; 0, -1] * transferT E v * !![1, 0; 0, -1] := by
  sorry

/-- **The `SL(2)` entry inequality** (behind Lemma `lem:entry`): for `det M = 1`, unit vectors `e, w`
and `a = ⟨w, M e⟩`, `(‖Mᵀ w‖² - a²)(‖M e‖² - a²) ≤ (|a| ‖M‖_F + 1)²`. -/
theorem sl2_entry_ineq (M : Matrix (Fin 2) (Fin 2) ℝ) (hM : M.det = 1) (e w : Fin 2 → ℝ)
    (he : e 0 ^ 2 + e 1 ^ 2 = 1) (hw : w 0 ^ 2 + w 1 ^ 2 = 1) :
    ((Mᵀ *ᵥ w) ⬝ᵥ (Mᵀ *ᵥ w) - (w ⬝ᵥ (M *ᵥ e)) ^ 2) * ((M *ᵥ e) ⬝ᵥ (M *ᵥ e) - (w ⬝ᵥ (M *ᵥ e)) ^ 2) ≤
      (|w ⬝ᵥ (M *ᵥ e)| * Real.sqrt (∑ i, ∑ j, M i j ^ 2) + 1) ^ 2 := by
  sorry

/-! ### The finite chain -/

variable {N : ℕ}

/-- Extension by zero of a function on `{0, …, N-1}` to `ℤ`. -/
def extZ (φ : Fin N → ℝ) (z : ℤ) : ℝ := if h : 0 ≤ z ∧ z < N then φ ⟨z.toNat, by omega⟩ else 0

/-- The potential `λ sgn(q)` on `{0, …, N-1}`, extended by `0`. -/
def potZ (lam : ℝ) (q : Cfg N) : ℤ → ℝ := extZ (fun i => lam * sgn (q i))

/-- `potZ` at a site of the chain. -/
theorem potZ_apply (lam : ℝ) (q : Cfg N) (i : Fin N) :
    potZ lam q (i : ℤ) = lam * sgn (q i) := by
  sorry

/-- Eigenvectors of `h_q = -J Δ + 2h diag(q)`, extended by zero, solve
`φ(z+1) + φ(z-1) + V z φ(z) = E φ(z)` on `{0, …, N-1}` with `V = -(2h/J) sgn q` and
`E = -E'/J` (`E'` the eigenvalue). -/
theorem eigvec_recursion (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) (α : Fin N) (z : ℤ) (hz0 : 0 ≤ z)
    (hzN : z < N) :
    extZ (fun i => (hq_isHermitian J h q).eigenvectorBasis α i) (z + 1) +
        extZ (fun i => (hq_isHermitian J h q).eigenvectorBasis α i) (z - 1) +
        potZ (-(2 * h / J)) q z * extZ (fun i => (hq_isHermitian J h q).eigenvectorBasis α i) z =
      (-(hq_isHermitian J h q).eigenvalues α / J) *
        extZ (fun i => (hq_isHermitian J h q).eigenvectorBasis α i) z := by
  sorry

/-- Every eigenvalue of `h_q` has modulus at most `2|J| + 2|h|`. -/
theorem abs_eigenvalue_hq_le (J h : ℝ) (q : Cfg N) (α : Fin N) :
    |(hq_isHermitian J h q).eigenvalues α| ≤ 2 * |J| + 2 * |h| := by
  sorry

end DFL
