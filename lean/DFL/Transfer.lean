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

namespace TransferAux

/-- A monic polynomial step: `(X - c) p - r` is monic of degree `n + 2` when `p` is monic of degree
`n + 1` and `r` has degree `n`. -/
lemma monic_step (c : ℝ) (p r : ℝ[X]) (n : ℕ) (hp : p.Monic) (hpd : p.natDegree = n + 1)
    (hrd : r.natDegree = n) :
    ((X - C c) * p - r).Monic ∧ ((X - C c) * p - r).natDegree = n + 2 := by
  have hm : ((X - C c) * p).Monic := (monic_X_sub_C c).mul hp
  have hd : ((X - C c) * p).natDegree = n + 2 := by
    rw [(monic_X_sub_C c).natDegree_mul hp, natDegree_X_sub_C, hpd]; ring
  have hlt : r.natDegree < ((X - C c) * p).natDegree := by rw [hd, hrd]; omega
  exact ⟨hm.sub_of_left (degree_lt_degree hlt),
    by rw [natDegree_sub_eq_left_of_natDegree_lt hlt, hd]⟩

/-- `dirLpoly V a (k+1)` is monic of degree `k`. -/
lemma dirLpoly_monic_natDegree (V : ℤ → ℝ) (a : ℤ) (k : ℕ) :
    (dirLpoly V a (k + 1)).Monic ∧ (dirLpoly V a (k + 1)).natDegree = k := by
  induction k using Nat.twoStepInduction with
  | zero => simp [dirLpoly]
  | one =>
    have h : dirLpoly V a (1 + 1) = X - C (V a) := by simp [dirLpoly]
    rw [h]
    exact ⟨monic_X_sub_C _, natDegree_X_sub_C _⟩
  | more n ih0 ih1 =>
    have h : dirLpoly V a (n + 2 + 1) =
        (X - C (V (a + ((n + 1 : ℕ) : ℤ)))) * dirLpoly V a (n + 1 + 1) - dirLpoly V a (n + 1) :=
      rfl
    rw [h]
    exact monic_step _ _ _ n ih1.1 ih1.2 ih0.2

/-- `dirRpoly V b (k+1)` is monic of degree `k`. -/
lemma dirRpoly_monic_natDegree (V : ℤ → ℝ) (b : ℤ) (k : ℕ) :
    (dirRpoly V b (k + 1)).Monic ∧ (dirRpoly V b (k + 1)).natDegree = k := by
  induction k using Nat.twoStepInduction with
  | zero => simp [dirRpoly]
  | one =>
    have h : dirRpoly V b (1 + 1) = X - C (V b) := by simp [dirRpoly]
    rw [h]
    exact ⟨monic_X_sub_C _, natDegree_X_sub_C _⟩
  | more n ih0 ih1 =>
    have h : dirRpoly V b (n + 2 + 1) =
        (X - C (V (b - ((n + 1 : ℕ) : ℤ)))) * dirRpoly V b (n + 1 + 1) - dirRpoly V b (n + 1) :=
      rfl
    rw [h]
    exact monic_step _ _ _ n ih1.1 ih1.2 ih0.2

end TransferAux

/-! ### Polynomials -/

theorem eval_dirLpoly (V : ℤ → ℝ) (E : ℝ) (a : ℤ) (k : ℕ) :
    (dirLpoly V a k).eval E = dirL V E a k := by
  induction k using Nat.twoStepInduction with
  | zero => simp [dirLpoly, dirL]
  | one => simp [dirLpoly, dirL]
  | more n ih0 ih1 => simp [dirLpoly, dirL, ih0, ih1]

theorem eval_dirRpoly (V : ℤ → ℝ) (E : ℝ) (b : ℤ) (k : ℕ) :
    (dirRpoly V b k).eval E = dirR V E b k := by
  induction k using Nat.twoStepInduction with
  | zero => simp [dirRpoly, dirR]
  | one => simp [dirRpoly, dirR]
  | more n ih0 ih1 => simp [dirRpoly, dirR, ih0, ih1]

theorem dirLpoly_monic (V : ℤ → ℝ) (a : ℤ) (k : ℕ) : (dirLpoly V a (k + 1)).Monic := by
  exact (TransferAux.dirLpoly_monic_natDegree V a k).1

theorem natDegree_dirLpoly (V : ℤ → ℝ) (a : ℤ) (k : ℕ) :
    (dirLpoly V a (k + 1)).natDegree = k := by
  exact (TransferAux.dirLpoly_monic_natDegree V a k).2

theorem dirRpoly_monic (V : ℤ → ℝ) (b : ℤ) (k : ℕ) : (dirRpoly V b (k + 1)).Monic := by
  exact (TransferAux.dirRpoly_monic_natDegree V b k).1

theorem natDegree_dirRpoly (V : ℤ → ℝ) (b : ℤ) (k : ℕ) :
    (dirRpoly V b (k + 1)).natDegree = k := by
  exact (TransferAux.dirRpoly_monic_natDegree V b k).2

/-! ### Locality -/

/-- `dirL V E a k` depends only on `V (a), …, V (a + k - 2)`. -/
theorem dirL_congr (V V' : ℤ → ℝ) (E : ℝ) (a : ℤ) (k : ℕ)
    (h : ∀ j : ℕ, j + 1 < k → V (a + j) = V' (a + j)) : dirL V E a k = dirL V' E a k := by
  revert h
  induction k using Nat.twoStepInduction with
  | zero => intro _; rfl
  | one => intro _; rfl
  | more n ih0 ih1 =>
    intro h
    simp only [dirL]
    rw [ih0 (fun j hj => h j (by omega)), ih1 (fun j hj => h j (by omega)), h n (by omega)]

/-- `dirR V E b k` depends only on `V (b), …, V (b - k + 2)`. -/
theorem dirR_congr (V V' : ℤ → ℝ) (E : ℝ) (b : ℤ) (k : ℕ)
    (h : ∀ j : ℕ, j + 1 < k → V (b - j) = V' (b - j)) : dirR V E b k = dirR V' E b k := by
  revert h
  induction k using Nat.twoStepInduction with
  | zero => intro _; rfl
  | one => intro _; rfl
  | more n ih0 ih1 =>
    intro h
    simp only [dirR]
    rw [ih0 (fun j hj => h j (by omega)), ih1 (fun j hj => h j (by omega)), h n (by omega)]

theorem dirLpoly_congr (V V' : ℤ → ℝ) (a : ℤ) (k : ℕ)
    (h : ∀ j : ℕ, j + 1 < k → V (a + j) = V' (a + j)) : dirLpoly V a k = dirLpoly V' a k := by
  revert h
  induction k using Nat.twoStepInduction with
  | zero => intro _; rfl
  | one => intro _; rfl
  | more n ih0 ih1 =>
    intro h
    simp only [dirLpoly]
    rw [ih0 (fun j hj => h j (by omega)), ih1 (fun j hj => h j (by omega)), h n (by omega)]

theorem dirRpoly_congr (V V' : ℤ → ℝ) (b : ℤ) (k : ℕ)
    (h : ∀ j : ℕ, j + 1 < k → V (b - j) = V' (b - j)) : dirRpoly V b k = dirRpoly V' b k := by
  revert h
  induction k using Nat.twoStepInduction with
  | zero => intro _; rfl
  | one => intro _; rfl
  | more n ih0 ih1 =>
    intro h
    simp only [dirRpoly]
    rw [ih0 (fun j hj => h j (by omega)), ih1 (fun j hj => h j (by omega)), h n (by omega)]

/-! ### Identities -/

namespace TransferAux

/-- One step of the constancy of the Wronskian of `dirL` and `dirR`. -/
lemma wronsk_step (V : ℤ → ℝ) (E : ℝ) (a b : ℤ) (i j : ℕ)
    (hij : a + (i : ℤ) = b - (j : ℤ)) :
    dirL V E a (i + 2) * dirR V E b (j + 1) - dirL V E a (i + 1) * dirR V E b j =
      dirL V E a (i + 1) * dirR V E b (j + 2) - dirL V E a i * dirR V E b (j + 1) := by
  simp only [dirL, dirR]
  rw [hij]
  ring

/-- The Wronskian of `dirL` and `dirR` is constant along the box `[a, b]`, `b = a + L - 1`. -/
lemma wronsk_const (V : ℤ → ℝ) (E : ℝ) (a b : ℤ) (L : ℕ) (hab : b = a + L - 1) :
    ∀ i j : ℕ, i + j = L →
      dirL V E a (i + 1) * dirR V E b (j + 1) - dirL V E a i * dirR V E b j =
        dirR V E b (L + 1) := by
  intro i
  induction i with
  | zero =>
    intro j hj
    obtain rfl : j = L := by omega
    simp [dirL]
  | succ n ih =>
    intro j hj
    rw [wronsk_step V E a b n j (by omega)]
    exact ih (j + 1) (by omega)

/-- Wronskian of the left Dirichlet solution with a solution `φ` of the equation on the box. -/
lemma wronskL_phi (V : ℤ → ℝ) (E : ℝ) (φ : ℤ → ℝ) (a : ℤ) (L : ℕ)
    (hφ : ∀ z, a ≤ z → z ≤ a + L - 1 → φ (z + 1) + φ (z - 1) + V z * φ z = E * φ z) :
    ∀ i : ℕ, i ≤ L →
      dirL V E a (i + 1) * φ (a + i - 1) - dirL V E a i * φ (a + i) = φ (a - 1) := by
  intro i
  induction i with
  | zero => intro _; simp [dirL]
  | succ n ih =>
    intro hn
    have h1 := ih (by omega)
    have h2 := hφ (a + n) (by omega) (by omega)
    rw [show a + ((n + 1 : ℕ) : ℤ) - 1 = a + n by push_cast; ring,
      show a + ((n + 1 : ℕ) : ℤ) = a + n + 1 by push_cast; ring]
    simp only [dirL]
    linear_combination h1 - dirL V E a (n + 1) * h2

/-- Wronskian of the right Dirichlet solution with a solution `φ` of the equation on the box. -/
lemma wronskR_phi (V : ℤ → ℝ) (E : ℝ) (φ : ℤ → ℝ) (a b : ℤ) (L : ℕ)
    (hab : b = a + L - 1)
    (hφ : ∀ z, a ≤ z → z ≤ a + L - 1 → φ (z + 1) + φ (z - 1) + V z * φ z = E * φ z) :
    ∀ j : ℕ, j ≤ L →
      dirR V E b j * φ (b - j) - dirR V E b (j + 1) * φ (b - j + 1) = -φ (b + 1) := by
  intro j
  induction j with
  | zero => intro _; simp [dirR]
  | succ n ih =>
    intro hn
    have h1 := ih (by omega)
    have h2 := hφ (b - n) (by omega) (by omega)
    rw [show b - ((n + 1 : ℕ) : ℤ) + 1 = b - n by push_cast; ring,
      show b - ((n + 1 : ℕ) : ℤ) = b - n - 1 by push_cast; ring]
    simp only [dirR]
    linear_combination h1 + dirR V E b (n + 1) * h2

/-- Action of a transfer matrix on a vector. -/
lemma transferT_mulVec (E v x y : ℝ) : transferT E v *ᵥ ![x, y] = ![(E - v) * x - y, x] := by
  ext k
  fin_cases k <;> simp [transferT, Matrix.mulVec, dotProduct, Fin.sum_univ_two]
  ring

/-- Transfer matrices propagate the left Dirichlet solution to the right. -/
lemma tprodL_dirL (V : ℤ → ℝ) (E : ℝ) (a : ℤ) (n m : ℕ) :
    tprodL V E (a + n) m *ᵥ ![dirL V E a (n + 1), dirL V E a n] =
      ![dirL V E a (n + m + 1), dirL V E a (n + m)] := by
  induction m with
  | zero => simp [tprodL]
  | succ m ih =>
    rw [tprodL, ← Matrix.mulVec_mulVec, ih, transferT_mulVec]
    have h1 : dirL V E a (n + (m + 1) + 1) =
        (E - V (a + n + m)) * dirL V E a (n + m + 1) - dirL V E a (n + m) := by
      rw [show n + (m + 1) + 1 = (n + m) + 2 by omega]
      simp only [dirL]
      rw [show a + ((n + m : ℕ) : ℤ) = a + n + m by push_cast; ring]
    rw [h1, show n + (m + 1) = n + m + 1 by omega]

/-- Transfer matrices propagate the right Dirichlet solution to the left. -/
lemma tprodR_dirR (V : ℤ → ℝ) (E : ℝ) (b s : ℤ) (m : ℕ) :
    ∀ n : ℕ, s + m + n = b + 1 →
      tprodR V E s m *ᵥ ![dirR V E b (n + 1), dirR V E b n] =
        ![dirR V E b (n + m + 1), dirR V E b (n + m)] := by
  induction m with
  | zero => intro n _; simp [tprodR]
  | succ m ih =>
    intro n hn
    push_cast at hn
    rw [tprodR, ← Matrix.mulVec_mulVec, transferT_mulVec]
    have h1 : (E - V (s + m)) * dirR V E b (n + 1) - dirR V E b n = dirR V E b (n + 2) := by
      simp only [dirR]
      rw [show b - (n : ℤ) = s + m by omega]
    rw [h1, ih (n + 1) (by push_cast; omega), show n + 1 + m = n + (m + 1) by omega]

end TransferAux

/-- Consecutive values of a Dirichlet solution do not both vanish. -/
theorem dirL_pair_ne_zero (V : ℤ → ℝ) (E : ℝ) (a : ℤ) (k : ℕ) :
    dirL V E a (k + 1) ≠ 0 ∨ dirL V E a k ≠ 0 := by
  induction k with
  | zero => left; simp [dirL]
  | succ n ih =>
    by_cases h1 : dirL V E a (n + 1) = 0
    · left
      have h0 : dirL V E a n ≠ 0 := by
        rcases ih with h | h
        · exact absurd h1 h
        · exact h
      simp only [dirL, h1, mul_zero, zero_sub, neg_ne_zero]
      exact h0
    · right; exact h1

/-- Consecutive values of a Dirichlet solution do not both vanish. -/
theorem dirR_pair_ne_zero (V : ℤ → ℝ) (E : ℝ) (b : ℤ) (k : ℕ) :
    dirR V E b (k + 1) ≠ 0 ∨ dirR V E b k ≠ 0 := by
  induction k with
  | zero => left; simp [dirR]
  | succ n ih =>
    by_cases h1 : dirR V E b (n + 1) = 0
    · left
      have h0 : dirR V E b n ≠ 0 := by
        rcases ih with h | h
        · exact absurd h1 h
        · exact h
      simp only [dirR, h1, mul_zero, zero_sub, neg_ne_zero]
      exact h0
    · right; exact h1

/-- Both Dirichlet solutions give the determinant of the box `[a, a+L-1]` (constancy of the
Wronskian; paper: Lemma `lem:green`(a)). -/
theorem dirL_eq_dirR (V : ℤ → ℝ) (E : ℝ) (a : ℤ) (L : ℕ) :
    dirL V E a (L + 1) = dirR V E (a + L - 1) (L + 1) := by
  have h := TransferAux.wronsk_const V E a (a + L - 1) L rfl L 0 (by simp)
  simpa [dirR] using h

/-- **Poisson identity**: if `φ` solves `φ(z+1) + φ(z-1) + V z φ(z) = E φ(z)` on the box
`[a, a+L-1]`, then at `x = a + i`, `P_B(E) φ(x) = p_a(E) φ(a-1) + p_b(E) φ(a+L)`. -/
theorem poisson (V : ℤ → ℝ) (E : ℝ) (φ : ℤ → ℝ) (a : ℤ) (L i : ℕ) (hi : i < L)
    (hφ : ∀ z, a ≤ z → z ≤ a + L - 1 → φ (z + 1) + φ (z - 1) + V z * φ z = E * φ z) :
    dirL V E a (L + 1) * φ (a + i) =
      dirR V E (a + L - 1) (L - i) * φ (a - 1) + dirL V E a (i + 1) * φ (a + L) := by
  obtain ⟨j, rfl⟩ : ∃ j, L = i + 1 + j := ⟨L - i - 1, by omega⟩
  have h1 := TransferAux.wronskL_phi V E φ a (i + 1 + j) hφ (i + 1) (by omega)
  have h2 := TransferAux.wronskR_phi V E φ a (a + ((i + 1 + j : ℕ) : ℤ) - 1) (i + 1 + j) rfl
    hφ j (by omega)
  have h3 := TransferAux.wronsk_const V E a (a + ((i + 1 + j : ℕ) : ℤ) - 1) (i + 1 + j) rfl
    (i + 1) j rfl
  rw [← dirL_eq_dirR] at h3
  rw [show i + 1 + j - i = j + 1 by omega]
  rw [show a + ((i + 1 : ℕ) : ℤ) - 1 = a + i by push_cast; ring,
    show a + ((i + 1 : ℕ) : ℤ) = a + i + 1 by push_cast; ring] at h1
  rw [show a + ((i + 1 + j : ℕ) : ℤ) - 1 - j + 1 = a + i + 1 by push_cast; ring,
    show a + ((i + 1 + j : ℕ) : ℤ) - 1 - j = a + i by push_cast; ring,
    show a + ((i + 1 + j : ℕ) : ℤ) - 1 + 1 = a + ((i + 1 + j : ℕ) : ℤ) by ring] at h2
  linear_combination dirR V E (a + ((i + 1 + j : ℕ) : ℤ) - 1) (j + 1) * h1 -
    dirL V E a (i + 1) * h2 - φ (a + i) * h3

/-- **Right end** (paper: Lemma `lem:green`(c)): `P_B(E) = (T_b ⋯ T_{x+1} Y)₀` with
`Y = (u₋(x+1), u₋(x))`, which depends only on `V` on `[a, x]`; note `p_b(E) = Y₁`. -/
theorem box_right (V : ℤ → ℝ) (E : ℝ) (a : ℤ) (L i : ℕ) (hi : i < L) :
    dirL V E a (L + 1) =
      (tprodL V E (a + i + 1) (L - 1 - i) *ᵥ ![dirL V E a (i + 2), dirL V E a (i + 1)]) 0 := by
  have h := TransferAux.tprodL_dirL V E a (i + 1) (L - 1 - i)
  rw [show a + ((i + 1 : ℕ) : ℤ) = a + i + 1 by push_cast; ring] at h
  rw [show i + 2 = i + 1 + 1 from rfl, h]
  simp only [Matrix.cons_val_zero]
  congr 1
  omega

/-- **Left end** (paper: Lemma `lem:green`(d)): `P_B(E) = (T_a ⋯ T_{x-1} Y')₀` with
`Y' = (u₊(x-1), u₊(x))`, which depends only on `V` on `[x, b]`; note `p_a(E) = Y'₁`. -/
theorem box_left (V : ℤ → ℝ) (E : ℝ) (a : ℤ) (L i : ℕ) (hi : i < L) :
    dirL V E a (L + 1) =
      (tprodR V E a i *ᵥ
        ![dirR V E (a + L - 1) (L - i + 1), dirR V E (a + L - 1) (L - i)]) 0 := by
  rw [TransferAux.tprodR_dirR V E (a + L - 1) a i (L - i) (by omega)]
  simp only [Matrix.cons_val_zero]
  rw [show L - i + i + 1 = L + 1 by omega]
  exact dirL_eq_dirR V E a L

/-- Transfer matrices have determinant one. -/
theorem det_transferT (E v : ℝ) : (transferT E v).det = 1 := by
  simp [transferT, Matrix.det_fin_two_of]

/-- `T_E(v)ᵀ = D T_E(v) D`, `D = diag(1, -1)`. -/
theorem transferT_transpose (E v : ℝ) :
    (transferT E v)ᵀ = !![1, 0; 0, -1] * transferT E v * !![1, 0; 0, -1] := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [transferT, Matrix.mul_apply, Fin.sum_univ_two]

namespace TransferAux

/-- The `SL(2)` entry inequality in coordinates. -/
lemma sl2_coord (p q r s e0 e1 w0 w1 : ℝ) (hM : p * s - q * r = 1)
    (he : e0 ^ 2 + e1 ^ 2 = 1) (hw : w0 ^ 2 + w1 ^ 2 = 1) :
    (((p * w0 + r * w1) ^ 2 + (q * w0 + s * w1) ^ 2) -
        (w0 * (p * e0 + q * e1) + w1 * (r * e0 + s * e1)) ^ 2) *
      (((p * e0 + q * e1) ^ 2 + (r * e0 + s * e1) ^ 2) -
        (w0 * (p * e0 + q * e1) + w1 * (r * e0 + s * e1)) ^ 2) ≤
      (|w0 * (p * e0 + q * e1) + w1 * (r * e0 + s * e1)| *
        Real.sqrt (p ^ 2 + q ^ 2 + r ^ 2 + s ^ 2) + 1) ^ 2 := by
  set A := w0 * (p * e0 + q * e1) + w1 * (r * e0 + s * e1) with hA
  set B := w0 * (-p * e1 + q * e0) + w1 * (-r * e1 + s * e0) with hB
  set C := -w1 * (p * e0 + q * e1) + w0 * (r * e0 + s * e1) with hC
  set D := -w1 * (-p * e1 + q * e0) + w0 * (-r * e1 + s * e0) with hD
  set F := p ^ 2 + q ^ 2 + r ^ 2 + s ^ 2 with hF
  have h1 : (p * w0 + r * w1) ^ 2 + (q * w0 + s * w1) ^ 2 - A ^ 2 = B ^ 2 := by
    linear_combination (-((p * w0 + r * w1) ^ 2 + (q * w0 + s * w1) ^ 2)) * he
  have h2 : (p * e0 + q * e1) ^ 2 + (r * e0 + s * e1) ^ 2 - A ^ 2 = C ^ 2 := by
    linear_combination (-((p * e0 + q * e1) ^ 2 + (r * e0 + s * e1) ^ 2)) * hw
  have h3 : A * D - B * C = 1 := by
    linear_combination ((w0 ^ 2 + w1 ^ 2) * (e0 ^ 2 + e1 ^ 2)) * hM +
      (w0 ^ 2 + w1 ^ 2) * he + hw
  have h4 : A ^ 2 + B ^ 2 + C ^ 2 + D ^ 2 = F := by
    linear_combination (F * (e0 ^ 2 + e1 ^ 2)) * hw + F * he
  have hDF : D ^ 2 ≤ F := by nlinarith [sq_nonneg A, sq_nonneg B, sq_nonneg C]
  have hD' : |D| ≤ Real.sqrt F := Real.abs_le_sqrt hDF
  have key : |A * D - 1| ≤ |A| * Real.sqrt F + 1 := by
    calc |A * D - 1| ≤ |A * D| + |(1 : ℝ)| := abs_sub _ _
      _ = |A| * |D| + 1 := by rw [abs_mul, abs_one]
      _ ≤ |A| * Real.sqrt F + 1 := by gcongr
  rw [h1, h2]
  calc B ^ 2 * C ^ 2 = (A * D - 1) ^ 2 := by
        have : B * C = A * D - 1 := by linarith
        rw [← mul_pow, this]
    _ = |A * D - 1| ^ 2 := (sq_abs _).symm
    _ ≤ (|A| * Real.sqrt F + 1) ^ 2 := pow_le_pow_left₀ (abs_nonneg _) key 2

end TransferAux

/-- **The `SL(2)` entry inequality** (behind Lemma `lem:entry`): for `det M = 1`, unit vectors `e, w`
and `a = ⟨w, M e⟩`, `(‖Mᵀ w‖² - a²)(‖M e‖² - a²) ≤ (|a| ‖M‖_F + 1)²`. -/
theorem sl2_entry_ineq (M : Matrix (Fin 2) (Fin 2) ℝ) (hM : M.det = 1) (e w : Fin 2 → ℝ)
    (he : e 0 ^ 2 + e 1 ^ 2 = 1) (hw : w 0 ^ 2 + w 1 ^ 2 = 1) :
    ((Mᵀ *ᵥ w) ⬝ᵥ (Mᵀ *ᵥ w) - (w ⬝ᵥ (M *ᵥ e)) ^ 2) * ((M *ᵥ e) ⬝ᵥ (M *ᵥ e) - (w ⬝ᵥ (M *ᵥ e)) ^ 2) ≤
      (|w ⬝ᵥ (M *ᵥ e)| * Real.sqrt (∑ i, ∑ j, M i j ^ 2) + 1) ^ 2 := by
  rw [Matrix.det_fin_two] at hM
  have hX : (Mᵀ *ᵥ w) ⬝ᵥ (Mᵀ *ᵥ w) =
      (M 0 0 * w 0 + M 1 0 * w 1) ^ 2 + (M 0 1 * w 0 + M 1 1 * w 1) ^ 2 := by
    simp [dotProduct, Matrix.mulVec, Fin.sum_univ_two]; ring
  have hY : (M *ᵥ e) ⬝ᵥ (M *ᵥ e) =
      (M 0 0 * e 0 + M 0 1 * e 1) ^ 2 + (M 1 0 * e 0 + M 1 1 * e 1) ^ 2 := by
    simp [dotProduct, Matrix.mulVec, Fin.sum_univ_two]; ring
  have hA : w ⬝ᵥ (M *ᵥ e) =
      w 0 * (M 0 0 * e 0 + M 0 1 * e 1) + w 1 * (M 1 0 * e 0 + M 1 1 * e 1) := by
    simp [dotProduct, Matrix.mulVec, Fin.sum_univ_two]
  have hF : ∑ i, ∑ j, M i j ^ 2 = M 0 0 ^ 2 + M 0 1 ^ 2 + M 1 0 ^ 2 + M 1 1 ^ 2 := by
    simp [Fin.sum_univ_two]; ring
  rw [hX, hY, hA, hF]
  exact TransferAux.sl2_coord _ _ _ _ _ _ _ _ hM he hw

/-! ### The finite chain -/

variable {N : ℕ}

/-- Extension by zero of a function on `{0, …, N-1}` to `ℤ`. -/
def extZ (φ : Fin N → ℝ) (z : ℤ) : ℝ := if h : 0 ≤ z ∧ z < N then φ ⟨z.toNat, by omega⟩ else 0

/-- The potential `λ sgn(q)` on `{0, …, N-1}`, extended by `0`. -/
def potZ (lam : ℝ) (q : Cfg N) : ℤ → ℝ := extZ (fun i => lam * sgn (q i))

namespace TransferAux

/-- `extZ` at a site of the chain. -/
lemma extZ_natCast (φ : Fin N → ℝ) (i : Fin N) : extZ φ (i : ℤ) = φ i := by
  have h : (0 : ℤ) ≤ ((i : ℕ) : ℤ) ∧ ((i : ℕ) : ℤ) < N :=
    ⟨by positivity, by exact_mod_cast i.isLt⟩
  rw [extZ, dite_eq_left h]
  congr 1

/-- `extZ φ t` as a sum over the sites of the chain. -/
lemma sum_ite_extZ (φ : Fin N → ℝ) (t : ℤ) :
    ∑ y : Fin N, (if ((y : ℕ) : ℤ) = t then φ y else 0) = extZ φ t := by
  by_cases ht : 0 ≤ t ∧ t < N
  · rw [Finset.sum_eq_single ⟨t.toNat, by omega⟩]
    · rw [ite_eq_left (by simp; omega), extZ, dite_eq_left ht]
    · intro y _ hy
      rw [ite_eq_right]
      intro h
      apply hy
      ext
      simp
      omega
    · intro h
      exact absurd (Finset.mem_univ _) h
  · rw [extZ, dite_eq_right ht]
    apply Finset.sum_eq_zero
    intro y _
    rw [ite_eq_right]
    intro h
    apply ht
    have := y.isLt
    omega

/-- The path adjacency matrix acting on a vector, via the extension by zero. -/
lemma pathAdj_mulVec_apply (φ : Fin N → ℝ) (x : Fin N) :
    (pathAdj N *ᵥ φ) x = extZ φ ((x : ℤ) + 1) + extZ φ ((x : ℤ) - 1) := by
  have hsplit : ∀ y : Fin N, pathAdj N x y * φ y =
      (if ((y : ℕ) : ℤ) = (x : ℤ) + 1 then φ y else 0) +
        (if ((y : ℕ) : ℤ) = (x : ℤ) - 1 then φ y else 0) := by
    intro y
    simp only [pathAdj]
    by_cases h1 : x.val + 1 = y.val
    · rw [ite_eq_left (Or.inl h1), ite_eq_left (by omega), ite_eq_right (by omega)]
      ring
    · by_cases h2 : y.val + 1 = x.val
      · rw [ite_eq_left (Or.inr h2), ite_eq_right (by omega), ite_eq_left (by omega)]
        ring
      · rw [ite_eq_right (by tauto), ite_eq_right (by omega), ite_eq_right (by omega)]
        ring
  simp only [Matrix.mulVec, dotProduct, hsplit, Finset.sum_add_distrib, sum_ite_extZ]

/-- The matrix `h_q` acting on a vector, via the extension by zero. -/
lemma hq_mulVec_apply (J h : ℝ) (q : Cfg N) (φ : Fin N → ℝ) (x : Fin N) :
    (hq J h q *ᵥ φ) x =
      -J * (extZ φ ((x : ℤ) + 1) + extZ φ ((x : ℤ) - 1)) + 2 * h * sgn (q x) * φ x := by
  rw [hq, Matrix.add_mulVec, Matrix.smul_mulVec, Matrix.smul_mulVec]
  simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, Matrix.mulVec_diagonal,
    pathAdj_mulVec_apply]
  ring

/-- A bound on `extZ`. -/
lemma abs_extZ_le (φ : Fin N → ℝ) (M : ℝ) (hM : 0 ≤ M) (h : ∀ y, |φ y| ≤ M)
    (t : ℤ) :
    |extZ φ t| ≤ M := by
  unfold extZ
  split_ifs
  · exact h _
  · simpa using hM

/-- `|sgn a| = 1`. -/
lemma abs_sgn (a : ZMod 2) : |sgn a| = 1 := by
  unfold sgn
  split_ifs <;> simp

end TransferAux

/-- `potZ` at a site of the chain. -/
theorem potZ_apply (lam : ℝ) (q : Cfg N) (i : Fin N) :
    potZ lam q (i : ℤ) = lam * sgn (q i) := by
  exact TransferAux.extZ_natCast _ i

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
  obtain ⟨x, rfl⟩ : ∃ x : Fin N, ((x : ℕ) : ℤ) = z :=
    ⟨⟨z.toNat, by omega⟩, by simp; omega⟩
  have hev := congrFun ((hq_isHermitian J h q).mulVec_eigenvectorBasis α) x
  rw [TransferAux.hq_mulVec_apply] at hev
  simp only [Pi.smul_apply, smul_eq_mul] at hev
  rw [potZ_apply, TransferAux.extZ_natCast]
  have hJJ : J * J⁻¹ = 1 := mul_inv_cancel₀ hJ
  linear_combination (-J⁻¹) * hev -
    (extZ (fun i => (hq_isHermitian J h q).eigenvectorBasis α i) ((x : ℤ) + 1) +
      extZ (fun i => (hq_isHermitian J h q).eigenvectorBasis α i) ((x : ℤ) - 1)) * hJJ

/-- Every eigenvalue of `h_q` has modulus at most `2|J| + 2|h|`. -/
theorem abs_eigenvalue_hq_le (J h : ℝ) (q : Cfg N) (α : Fin N) :
    |(hq_isHermitian J h q).eigenvalues α| ≤ 2 * |J| + 2 * |h| := by
  have : Nonempty (Fin N) := ⟨α⟩
  obtain ⟨x, -, hxmax⟩ := Finset.exists_max_image Finset.univ
    (fun y => |(hq_isHermitian J h q).eigenvectorBasis α y|) Finset.univ_nonempty
  have hnorm : ‖(hq_isHermitian J h q).eigenvectorBasis α‖ = 1 :=
    (hq_isHermitian J h q).eigenvectorBasis.orthonormal.1 α
  have hpos : 0 < |(hq_isHermitian J h q).eigenvectorBasis α x| := by
    refine lt_of_le_of_ne (abs_nonneg _) (fun hcon => ?_)
    have h0 : ∀ y, (hq_isHermitian J h q).eigenvectorBasis α y = 0 := fun y =>
      abs_nonpos_iff.mp ((hxmax y (Finset.mem_univ y)).trans hcon.symm.le)
    rw [EuclideanSpace.norm_eq] at hnorm
    simp [h0] at hnorm
  have hev := congrFun ((hq_isHermitian J h q).mulVec_eigenvectorBasis α) x
  rw [TransferAux.hq_mulVec_apply] at hev
  simp only [Pi.smul_apply, smul_eq_mul] at hev
  set φ : Fin N → ℝ := fun i => (hq_isHermitian J h q).eigenvectorBasis α i with hφ
  have hb : ∀ t : ℤ, |extZ φ t| ≤ |φ x| := fun t =>
    TransferAux.abs_extZ_le φ _ (abs_nonneg _) (fun y => hxmax y (Finset.mem_univ y)) t
  have key : |(hq_isHermitian J h q).eigenvalues α| * |φ x| ≤ (2 * |J| + 2 * |h|) * |φ x| := by
    rw [← abs_mul, ← hev]
    calc |-J * (extZ φ ((x : ℤ) + 1) + extZ φ ((x : ℤ) - 1)) + 2 * h * sgn (q x) * φ x|
        ≤ |-J * (extZ φ ((x : ℤ) + 1) + extZ φ ((x : ℤ) - 1))| +
            |2 * h * sgn (q x) * φ x| := abs_add_le _ _
      _ = |J| * |extZ φ ((x : ℤ) + 1) + extZ φ ((x : ℤ) - 1)| + 2 * |h| * |φ x| := by
          rw [abs_mul, abs_mul, abs_mul, abs_mul, abs_neg, TransferAux.abs_sgn, abs_two, mul_one]
      _ ≤ |J| * (|φ x| + |φ x|) + 2 * |h| * |φ x| := by
          gcongr
          exact (abs_add_le _ _).trans (add_le_add (hb _) (hb _))
      _ = (2 * |J| + 2 * |h|) * |φ x| := by ring
  exact le_of_mul_le_mul_right key hpos

end DFL
