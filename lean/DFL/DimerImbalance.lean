import DFL.Imbalance
import DFL.Dimer
import DFL.FreeDyn

/-!
# The imbalance of the dimer state (paper: Proposition `prop:dimerimbalance`(i),(ii))

For dimer charge configurations (`q_{2k} = q_{2k+1}`, `N` even), every eigenvector `φ` of `h_q`
with eigenvalue `E` satisfies the identity
`((2h)² - E²) ∑_x (-1)^x φ(x)² = J² (φ(0)² - φ(N-1)²)`
(0-indexed sites; the paper writes it with 1-indexed sites and `λ = 2h/J`). Hence the staggered
weight of eigenvectors away from the critical energies `±2h` is a boundary term, which is the
mechanism behind the vanishing of the dimer imbalance. We also identify the infinite-time
imbalance of the dimer LGT state with the dimer average of `N⁻¹ ∑_α a_α²`. Part (iii) of the
paper's proposition, `Ī_N → 0` with an explicit rate, is proved in `DimerImbalanceLim.lean`.
-/

open Matrix Finset Filter Topology
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-! ### The dimer identity for eigenvectors of `h_q` -/

namespace DimerImbalanceAux

/-- Extension by zero of a vector on `Fin N` to a sequence indexed by `ℕ`. -/
def ext0 (φ : Fin N → ℝ) (i : ℕ) : ℝ := if h : i < N then φ ⟨i, h⟩ else 0

/-- The left neighbour `f (i - 1)` of a sequence, with value `0` at `i = 0`. -/
def prev0 (f : ℕ → ℝ) : ℕ → ℝ
  | 0 => 0
  | i + 1 => f i

/-- The extension by zero agrees with the vector on the chain. -/
lemma ext0_val (φ : Fin N → ℝ) (x : Fin N) : ext0 φ x.val = φ x := by
  simp [ext0, x.isLt]

/-- The extension by zero vanishes beyond the chain. -/
lemma ext0_of_le (φ : Fin N → ℝ) {i : ℕ} (hi : N ≤ i) : ext0 φ i = 0 := by
  simp [ext0, not_lt.mpr hi]

/-- The path adjacency matrix acts as `(Δφ)(x) = φ(x-1) + φ(x+1)` (Dirichlet ends). -/
lemma pathAdj_mulVec_apply (φ : Fin N → ℝ) (x : Fin N) :
    (pathAdj N *ᵥ φ) x = prev0 (ext0 φ) x.val + ext0 φ (x.val + 1) := by
  have h1 : (pathAdj N *ᵥ φ) x = ∑ j ∈ range N,
      ((if x.val + 1 = j then ext0 φ j else 0) + (if j + 1 = x.val then ext0 φ j else 0)) := by
    rw [← Fin.sum_univ_eq_sum_range (fun j => (if x.val + 1 = j then ext0 φ j else 0) +
      (if j + 1 = x.val then ext0 φ j else 0)) N]
    simp only [mulVec, dotProduct, pathAdj]
    refine sum_congr rfl fun y _ => ?_
    rw [ext0_val]
    by_cases h1 : x.val + 1 = y.val
    · have h2 : ¬ y.val + 1 = x.val := by omega
      simp [h1, h2]
    · by_cases h2 : y.val + 1 = x.val <;> simp [h1, h2]
  rw [h1, sum_add_distrib, sum_ite_eq, add_comm]
  congr 1
  · obtain ⟨xv, hx⟩ := x
    cases xv with
    | zero => simp [prev0]
    | succ k =>
      simp only [prev0, Nat.add_right_cancel_iff]
      rw [sum_ite_eq']
      simp only [mem_range, ite_eq_left_iff]
      omega
  · split_ifs with h2
    · rfl
    · rw [ext0_of_le]
      simp only [mem_range, not_lt] at h2
      exact h2

/-- `∑_{i<2m} g i = ∑_{k<m} (g (2k) + g (2k+1))`. -/
lemma sum_range_two_mul (g : ℕ → ℝ) (m : ℕ) :
    ∑ i ∈ range (2 * m), g i = ∑ k ∈ range m, (g (2 * k) + g (2 * k + 1)) := by
  induction m with
  | zero => simp
  | succ m ih =>
    rw [show 2 * (m + 1) = 2 * m + 1 + 1 by ring, sum_range_succ, sum_range_succ, ih,
      sum_range_succ]
    ring

/-- Bipartiteness of the path: `∑_x (-1)^x φ(x) (φ(x-1) + φ(x+1)) = 0`. -/
lemma sum_stag_cross (f : ℕ → ℝ) (M : ℕ) (hf : f (M + 1) = 0) :
    ∑ i ∈ range (M + 1), (-1 : ℝ) ^ i * f i * (prev0 f i + f (i + 1)) = 0 := by
  have h1 : ∑ i ∈ range (M + 1), (-1 : ℝ) ^ i * f i * prev0 f i =
      -∑ i ∈ range M, (-1 : ℝ) ^ i * f i * f (i + 1) := by
    rw [sum_range_succ', ← sum_neg_distrib]
    simp only [prev0, mul_zero, add_zero]
    refine sum_congr rfl fun i _ => ?_
    ring
  have h2 : ∑ i ∈ range (M + 1), (-1 : ℝ) ^ i * f i * f (i + 1) =
      ∑ i ∈ range M, (-1 : ℝ) ^ i * f i * f (i + 1) := by
    rw [sum_range_succ, hf, mul_zero, add_zero]
  simp only [mul_add, sum_add_distrib, h1, h2]
  ring

/-- `∑_{x<M+1} (-1)^x φ(x-1)² = -∑_{x<M} (-1)^x φ(x)²`. -/
lemma sum_stag_prev_sq (f : ℕ → ℝ) (M : ℕ) :
    ∑ i ∈ range (M + 1), (-1 : ℝ) ^ i * prev0 f i ^ 2 =
      -∑ i ∈ range M, (-1 : ℝ) ^ i * f i ^ 2 := by
  rw [sum_range_succ', ← sum_neg_distrib]
  simp only [prev0, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, mul_zero, add_zero]
  refine sum_congr rfl fun i _ => ?_
  ring

/-- `∑_{x<M+1} (-1)^x φ(x+1)² = φ(0)² - ∑_{x<M+1} (-1)^x φ(x)²` when `φ(M+1) = 0`. -/
lemma sum_stag_next_sq (f : ℕ → ℝ) (M : ℕ) (hf : f (M + 1) = 0) :
    ∑ i ∈ range (M + 1), (-1 : ℝ) ^ i * f (i + 1) ^ 2 =
      f 0 ^ 2 - ∑ i ∈ range (M + 1), (-1 : ℝ) ^ i * f i ^ 2 := by
  have h1 := sum_range_succ (fun i => (-1 : ℝ) ^ i * f i ^ 2) (M + 1)
  have h2 := sum_range_succ' (fun i => (-1 : ℝ) ^ i * f i ^ 2) (M + 1)
  have h3 : ∑ i ∈ range (M + 1), (-1 : ℝ) ^ (i + 1) * f (i + 1) ^ 2 =
      -∑ i ∈ range (M + 1), (-1 : ℝ) ^ i * f (i + 1) ^ 2 := by
    rw [← sum_neg_distrib]
    refine sum_congr rfl fun i _ => ?_
    ring
  simp only [hf, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, mul_zero,
    add_zero] at h1
  simp only [pow_zero, one_mul] at h2
  rw [h3] at h2
  linarith

/-- The dimer identity for sequences: if `f` vanishes at `2m` and solves the eigenvalue equation
`J (f(i-1) + f(i+1)) = (2h v_i - E) f(i)` on `[0, 2m)` with `v_i = ±1` constant on dimers, then
`((2h)² - E²) ∑_{i<2m} (-1)^i f(i)² = J² (f(0)² - f(2m-1)²)`. -/
lemma dimer_identity_seq (f v : ℕ → ℝ) (J h E : ℝ) (m : ℕ)
    (hf : f (2 * m) = 0)
    (heig : ∀ i < 2 * m, J * (prev0 f i + f (i + 1)) = (2 * h * v i - E) * f i)
    (hv : ∀ i < 2 * m, v i ^ 2 = 1)
    (hd : ∀ k < m, v (2 * k + 1) = v (2 * k)) :
    ((2 * h) ^ 2 - E ^ 2) * ∑ i ∈ range (2 * m), (-1 : ℝ) ^ i * f i ^ 2 =
      J ^ 2 * (f 0 ^ 2 - f (2 * m - 1) ^ 2) := by
  rcases Nat.eq_zero_or_pos m with rfl | hm
  · simp
  obtain ⟨M, hM⟩ : ∃ M, 2 * m = M + 1 := ⟨2 * m - 1, by omega⟩
  have hMo : Odd M := ⟨m - 1, by omega⟩
  have hf' : f (M + 1) = 0 := hM ▸ hf
  -- bipartiteness
  have hP : ∑ i ∈ range (2 * m), (-1 : ℝ) ^ i * f i * (prev0 f i + f (i + 1)) = 0 := by
    rw [hM]; exact sum_stag_cross f M hf'
  -- the staggered norm of `Δφ`
  have hW : ∑ i ∈ range (2 * m), (-1 : ℝ) ^ i * (prev0 f i + f (i + 1)) ^ 2 =
      -2 * ∑ i ∈ range (2 * m), (-1 : ℝ) ^ i * f i ^ 2 + f 0 ^ 2 - f (2 * m - 1) ^ 2 +
        2 * ∑ i ∈ range (2 * m), (-1 : ℝ) ^ i * prev0 f i * f (i + 1) := by
    have e1 : ∀ i, (-1 : ℝ) ^ i * (prev0 f i + f (i + 1)) ^ 2 =
        (-1 : ℝ) ^ i * prev0 f i ^ 2 + (-1 : ℝ) ^ i * f (i + 1) ^ 2 +
          2 * ((-1 : ℝ) ^ i * prev0 f i * f (i + 1)) := fun i => by ring
    have hlast : ∑ i ∈ range (M + 1), (-1 : ℝ) ^ i * f i ^ 2 =
        ∑ i ∈ range M, (-1 : ℝ) ^ i * f i ^ 2 - f M ^ 2 := by
      rw [sum_range_succ, hMo.neg_one_pow]
      ring
    simp only [e1, sum_add_distrib, ← mul_sum]
    rw [hM, sum_stag_prev_sq, sum_stag_next_sq f M hf', hlast, Nat.add_sub_cancel]
    ring
  -- the dimer pairing
  have hC : J * ∑ i ∈ range (2 * m), (-1 : ℝ) ^ i * prev0 f i * f (i + 1) =
      J * ∑ i ∈ range (2 * m), (-1 : ℝ) ^ i * f i ^ 2 := by
    rw [mul_sum, mul_sum, sum_range_two_mul, sum_range_two_mul]
    refine sum_congr rfl fun k hk => ?_
    rw [mem_range] at hk
    have e1 := heig (2 * k) (by omega)
    have e2 := heig (2 * k + 1) (by omega)
    rw [hd k hk] at e2
    have hp : prev0 f (2 * k + 1) = f (2 * k) := rfl
    rw [hp] at e2 ⊢
    have p1 : (-1 : ℝ) ^ (2 * k) = 1 := by rw [pow_mul]; norm_num
    have p2 : (-1 : ℝ) ^ (2 * k + 1) = -1 := by rw [pow_succ, pow_mul]; norm_num
    rw [p1, p2]
    linear_combination f (2 * k + 1) * e1 - f (2 * k) * e2
  -- the eigenvalue equation, summed
  have h1 : ((2 * h) ^ 2 - E ^ 2) * ∑ i ∈ range (2 * m), (-1 : ℝ) ^ i * f i ^ 2 =
      J ^ 2 * ∑ i ∈ range (2 * m), (-1 : ℝ) ^ i * (prev0 f i + f (i + 1)) ^ 2 +
        2 * E * J * ∑ i ∈ range (2 * m), (-1 : ℝ) ^ i * f i * (prev0 f i + f (i + 1)) := by
    rw [mul_sum, mul_sum, mul_sum, ← sum_add_distrib]
    refine sum_congr rfl fun i hi => ?_
    have e := heig i (mem_range.mp hi)
    have ev := hv i (mem_range.mp hi)
    linear_combination (-(-1 : ℝ) ^ i * (J * (prev0 f i + f (i + 1)) +
      (2 * h * v i - E) * f i + 2 * E * f i)) * e + (-(-1 : ℝ) ^ i * 4 * h ^ 2 * f i ^ 2) * ev
  linear_combination h1 + 2 * J * hC + 2 * E * J * hP + J ^ 2 * hW

/-- The sign `sgn a` squares to one. -/
lemma sgn_sq (a : ZMod 2) : sgn a ^ 2 = 1 := by
  unfold sgn; split_ifs <;> norm_num

/-- The eigenvalue equation of `h_q` in coordinates: `J (φ(x-1) + φ(x+1)) = (2h q_x - E) φ(x)`. -/
lemma hq_eig_apply (J h : ℝ) (q : Cfg N) (φ : Fin N → ℝ) (E : ℝ) (hφ : hq J h q *ᵥ φ = E • φ)
    (x : Fin N) :
    J * (prev0 (ext0 φ) x.val + ext0 φ (x.val + 1)) = (2 * h * sgn (q x) - E) * ext0 φ x.val := by
  have hE := congrFun hφ x
  rw [DFL.hq, add_mulVec, smul_mulVec, smul_mulVec, Pi.add_apply, Pi.smul_apply, Pi.smul_apply,
    mulVec_diagonal, pathAdj_mulVec_apply, Pi.smul_apply, smul_eq_mul, smul_eq_mul,
    smul_eq_mul] at hE
  rw [ext0_val]
  linear_combination -hE

end DimerImbalanceAux

/-- **Dimer identity** (paper: Proposition `prop:dimerimbalance`(i), 0-indexed). -/
theorem dimer_stagWeight_identity (J h : ℝ) (hN : Even N) (hN0 : 0 < N) (q : Cfg N)
    (hq : IsDimer q) (α : Fin N) :
    ((2 * h) ^ 2 - (hq_isHermitian J h q).eigenvalues α ^ 2) * stagWeight J h q α =
      J ^ 2 * (((hq_isHermitian J h q).eigenvectorBasis α ⟨0, hN0⟩) ^ 2 -
        ((hq_isHermitian J h q).eigenvectorBasis α ⟨N - 1, by omega⟩) ^ 2) := by
  obtain ⟨m, hm⟩ := hN.two_dvd
  set φ : Fin N → ℝ := ⇑((hq_isHermitian J h q).eigenvectorBasis α) with hφ
  have key := DimerImbalanceAux.dimer_identity_seq (DimerImbalanceAux.ext0 φ)
    (DimerImbalanceAux.ext0 (fun x => sgn (q x))) J h ((hq_isHermitian J h q).eigenvalues α) m
    (DimerImbalanceAux.ext0_of_le _ (by omega))
    (fun i hi => by
      have h1 := DimerImbalanceAux.hq_eig_apply J h q φ _
        ((hq_isHermitian J h q).mulVec_eigenvectorBasis α) ⟨i, by omega⟩
      have h2 : DimerImbalanceAux.ext0 (fun x => sgn (q x)) i = sgn (q ⟨i, by omega⟩) :=
        DimerImbalanceAux.ext0_val (fun x => sgn (q x)) ⟨i, by omega⟩
      rw [h2]
      exact h1)
    (fun i hi => by
      have h2 : DimerImbalanceAux.ext0 (fun x => sgn (q x)) i = sgn (q ⟨i, by omega⟩) :=
        DimerImbalanceAux.ext0_val (fun x => sgn (q x)) ⟨i, by omega⟩
      rw [h2]
      exact DimerImbalanceAux.sgn_sq _)
    (fun k hk => by
      have h1 : DimerImbalanceAux.ext0 (fun x => sgn (q x)) (2 * k + 1) =
          sgn (q ⟨2 * k + 1, by omega⟩) :=
        DimerImbalanceAux.ext0_val (fun x => sgn (q x)) ⟨2 * k + 1, by omega⟩
      have h2 : DimerImbalanceAux.ext0 (fun x => sgn (q x)) (2 * k) =
          sgn (q ⟨2 * k, by omega⟩) :=
        DimerImbalanceAux.ext0_val (fun x => sgn (q x)) ⟨2 * k, by omega⟩
      rw [h1, h2, hq ⟨2 * k, by omega⟩ (by simp) (by simp only; omega)])
  rw [← hm] at key
  have e0 : DimerImbalanceAux.ext0 φ 0 = φ ⟨0, hN0⟩ := DimerImbalanceAux.ext0_val φ ⟨0, hN0⟩
  have e1 : DimerImbalanceAux.ext0 φ (N - 1) = φ ⟨N - 1, by omega⟩ :=
    DimerImbalanceAux.ext0_val φ ⟨N - 1, by omega⟩
  have es : ∑ i ∈ range N, (-1 : ℝ) ^ i * DimerImbalanceAux.ext0 φ i ^ 2 =
      stagWeight J h q α := by
    rw [stagWeight, ← Fin.sum_univ_eq_sum_range]
    simp only [DimerImbalanceAux.ext0_val, hφ]
  rw [e0, e1, es] at key
  exact key

/-- The density-wave imbalance of the LGT started from the dimer state (unnormalized
`dimerInit`, whose squared norm is `2^{N/2}`). -/
def dimerImbalance (J h : ℝ) (N : ℕ) (t : ℝ) : ℝ :=
  (2 / N : ℝ) * (1 / 2 ^ (N / 2) : ℝ) * ∑ x : Fin N, (-1 : ℝ) ^ x.val *
    (expect (cplx (diagonal (fun n : Cfg N => occ (n x)) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)))
      (evo (H J h) t *ᵥ dimerInit N)).re

/-! ### The time-averaged imbalance of the dimer state -/

namespace DimerImbalanceAux

/-- Rows of a unitary matrix are normalized: `∑_z |U_{xz}|² = 1`. -/
lemma sum_normSq_row {ι : Type*} [Fintype ι] [DecidableEq ι] (U : Matrix ι ι ℂ)
    (hU : U * Uᴴ = 1) (x : ι) : ∑ z, Complex.normSq (U x z) = 1 := by
  have h := congrFun (congrFun hU x) x
  simp only [one_apply_eq] at h
  have : ((∑ z, Complex.normSq (U x z) : ℝ) : ℂ) = 1 := by
    rw [← h, Complex.ofReal_sum, Matrix.mul_apply]
    refine Finset.sum_congr rfl (fun z _ => ?_)
    rw [conjTranspose_apply, ← Complex.mul_conj]
    rfl
  exact_mod_cast this

/-- `∑_z (-1)^z a_z = 2 ∑_{z even} a_z - ∑_z a_z`. -/
lemma sum_neg_one_pow_mul (a : Fin N → ℝ) :
    ∑ z : Fin N, (-1 : ℝ) ^ z.val * a z =
      2 * ∑ z ∈ univ.filter (fun x : Fin N => Even x.val), a z - ∑ z, a z := by
  rw [Finset.sum_filter, Finset.mul_sum, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl (fun z _ => ?_)
  rcases Nat.even_or_odd z.val with he | ho
  · rw [he.neg_one_pow, ite_eq_left he]; ring
  · rw [ho.neg_one_pow, ite_eq_right (Nat.not_even_iff_odd.2 ho)]; ring

/-- `∑_{x < N} (-1)^x = 0` for even `N`. -/
lemma sum_neg_one_pow_fin (hN : Even N) : ∑ x : Fin N, (-1 : ℝ) ^ x.val = 0 := by
  rw [Fin.sum_univ_eq_sum_range (fun i => (-1 : ℝ) ^ i) N]
  obtain ⟨k, rfl⟩ := hN
  induction k with
  | zero => simp
  | succ k ih =>
    rw [show k + 1 + (k + 1) = (k + k) + 1 + 1 by ring, Finset.sum_range_succ,
      Finset.sum_range_succ, ih, pow_succ, Even.neg_one_pow ⟨k, rfl⟩]
    ring

/-- `Re Tr(Γ U Γ U†) = ∑_{x,z} (-1)^x (-1)^z |U_{xz}|²`. -/
lemma trace_stag_re (U : Matrix (Fin N) (Fin N) ℂ) :
    (trace (cplx (stag N) * U * cplx (stag N) * Uᴴ)).re =
      ∑ x : Fin N, ∑ z : Fin N,
        (-1 : ℝ) ^ x.val * ((-1 : ℝ) ^ z.val * Complex.normSq (U x z)) := by
  have hs : cplx (stag N) = diagonal (fun x : Fin N => (((-1 : ℝ) ^ x.val : ℝ) : ℂ)) := by
    unfold cplx stag; exact diagonal_map (by simp)
  rw [hs, trace, Complex.re_sum]
  refine Finset.sum_congr rfl (fun x _ => ?_)
  rw [Matrix.diag_apply, Matrix.mul_apply, Complex.re_sum]
  refine Finset.sum_congr rfl (fun z _ => ?_)
  rw [mul_diagonal, diagonal_mul, conjTranspose_apply]
  have : (((-1 : ℝ) ^ x.val : ℝ) : ℂ) * U x z * (((-1 : ℝ) ^ z.val : ℝ) : ℂ) * star (U x z) =
      (((-1 : ℝ) ^ x.val * ((-1 : ℝ) ^ z.val * Complex.normSq (U x z)) : ℝ) : ℂ) := by
    rw [Complex.ofReal_mul, Complex.ofReal_mul, ← Complex.mul_conj]
    simp only [Complex.star_def]
    ring
  rw [this, Complex.ofReal_re]

/-- The one-body imbalance through the propagator (`N` even):
`imbOne = (2/N) ∑_x (-1)^x ∑_{z even} |U_{xz}|²`. -/
lemma imbOne_eq (J h : ℝ) (hN : Even N) (q : Cfg N) (t : ℝ) :
    imbOne J h q t = (2 / N : ℝ) * ∑ x : Fin N, (-1 : ℝ) ^ x.val *
      ∑ z ∈ univ.filter (fun x : Fin N => Even x.val),
        Complex.normSq (evo (hq J h q) t x z) := by
  have hrow := sum_normSq_row (evo (hq J h q) t) (evo_unitary' _ (hq_isSymm J h q) t)
  have h1 : ∀ x : Fin N, ∑ z : Fin N, (-1 : ℝ) ^ x.val *
      ((-1 : ℝ) ^ z.val * Complex.normSq (evo (hq J h q) t x z)) =
      2 * ((-1 : ℝ) ^ x.val * ∑ z ∈ univ.filter (fun x : Fin N => Even x.val),
        Complex.normSq (evo (hq J h q) t x z)) - (-1 : ℝ) ^ x.val := by
    intro x
    rw [← Finset.mul_sum, sum_neg_one_pow_mul, hrow x]
    ring
  rw [imbOne, trace_stag_re, Finset.sum_congr rfl (fun x _ => h1 x), Finset.sum_sub_distrib,
    ← Finset.mul_sum, sum_neg_one_pow_fin hN]
  ring

/-- The density of the dimer LGT state at site `x`: the dimer sum of the densities of the evolved
density-wave Fock state, `∑_{q dimer} ∑_{z even} |(e^{-ith_q})_{xz}|²`. -/
lemma dimer_density_re (J h : ℝ) (hN : Even N) (x : Fin N) (t : ℝ) :
    (expect (cplx (diagonal (fun n : Cfg N => occ (n x)) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)))
      (evo (H J h) t *ᵥ dimerInit N)).re =
      ∑ q ∈ univ.filter (fun q : Cfg N => IsDimer q),
        ∑ z ∈ univ.filter (fun x : Fin N => Even x.val),
          Complex.normSq (evo (hq J h q) t x z) := by
  rw [sector_average_dimer hN J h _ (gaugeInv_matter _) t, Complex.re_sum]
  refine sum_congr rfl fun q _ => ?_
  rw [sectorBlock_matter, Hsec, expect_evo_sub_smul_one]
  have hv : (Pi.single (fockCfg (univ.filter (fun x : Fin N => Even x.val))) 1 : Cfg N → ℂ) =
      fockVec (univ.filter (fun x : Fin N => Even x.val)) := rfl
  rw [hv, density _ (hq_isSymm J h q), Complex.ofReal_re]

/-- The dimer imbalance is the dimer average of the one-body imbalances (`N` even). -/
lemma dimerImbalance_eq (J h : ℝ) (hN : Even N) (t : ℝ) :
    dimerImbalance J h N t = (1 / 2 ^ (N / 2) : ℝ) *
      ∑ q ∈ univ.filter (fun q : Cfg N => IsDimer q), imbOne J h q t := by
  simp only [dimerImbalance, dimer_density_re J h hN, imbOne_eq J h hN, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine sum_congr rfl fun q _ => sum_congr rfl fun x _ => sum_congr rfl fun z _ => ?_
  ring

set_option backward.isDefEq.respectTransparency false in
open scoped Matrix.Norms.Operator in
/-- The propagator `t ↦ e^{-itA}` is continuous. -/
lemma continuous_evo {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℝ) :
    Continuous (fun t : ℝ => evo A t) := by
  unfold evo
  exact NormedSpace.exp_continuous.comp (by fun_prop)

/-- The one-body imbalance is continuous in time. -/
lemma continuous_imbOne (J h : ℝ) (q : Cfg N) : Continuous (fun t => imbOne J h q t) := by
  have hU := continuous_evo (hq J h q)
  unfold imbOne
  exact continuous_const.mul (Complex.continuous_re.comp
    ((((continuous_const.matrix_mul hU).matrix_mul continuous_const).matrix_mul
      hU.matrix_conjTranspose).matrix_trace))

end DimerImbalanceAux

/-- **Infinite-time imbalance of the dimer state** (paper: Proposition `prop:dimerimbalance`(ii)):
the time average converges to the dimer average of `N⁻¹ ∑_α a_α²`. -/
theorem dimerImbalance_timeAvg (J h : ℝ) (hJ : J ≠ 0) (hN : Even N) :
    Tendsto (fun T : ℝ => (1 / T) * ∫ t in (0 : ℝ)..T, dimerImbalance J h N t) atTop
      (𝓝 ((1 / 2 ^ (N / 2) : ℝ) *
        ∑ q ∈ univ.filter (fun q : Cfg N => IsDimer q), imbInf J h q)) := by
  have heq : (fun T : ℝ => (1 / T) * ∫ t in (0 : ℝ)..T, dimerImbalance J h N t) =
      fun T => (1 / 2 ^ (N / 2) : ℝ) * ∑ q ∈ univ.filter (fun q : Cfg N => IsDimer q),
        ((1 / T) * ∫ t in (0 : ℝ)..T, imbOne J h q t) := by
    funext T
    simp_rw [DimerImbalanceAux.dimerImbalance_eq J h hN]
    rw [intervalIntegral.integral_const_mul, intervalIntegral.integral_finsetSum
      (fun q _ => (DimerImbalanceAux.continuous_imbOne J h q).intervalIntegrable 0 T),
      Finset.mul_sum, Finset.mul_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl (fun q _ => ?_)
    ring
  rw [heq]
  exact (tendsto_finsetSum _ (fun q _ => imbOne_timeAvg J h hJ q)).const_mul _

end DFL
