import DFL.Imbalance
import DFL.Normalization
import DFL.DimerTransfer

/-!
# Single sectors (paper: Remarks `rem:gaussproj`(i),(ii), `rem:main-imbalance`(ii), `rem:dimersigns`)

The clean sector is `q ≡ 1` (all Gauss charges `+1`; the label `0 : ZMod 2`), with
`h_{1_N} = -J Δ_N + 2h`.

* `corrQ_clean_ends` (`rem:gaussproj`(i)): no decay of the eigenfunction correlator in the clean sector,
  `Q(1, N) = ∑_α |φ_α(1) φ_α(N)| = 1` for every `N`, the largest possible value. (The eigenvectors of the
  reflection-symmetric `h_{1_N}` are even or odd, since the spectrum is simple.)
* `stagWeight_clean`, `imbInf_clean`, `imbOne_clean_timeAvg` (`rem:gaussproj`(ii)): for even `N`,
  `Γ h + h Γ = 4h Γ` and `2h` is not an eigenvalue, so every eigenvector has `a_α = ⟨φ_α, Γ φ_α⟩ = 0` and the
  infinite-time imbalance of the clean sector vanishes.
* `stagWeight_h_zero` (`rem:main-imbalance`(ii)): for `h = 0` and even `N`, the chiral symmetry
  `Γ h_q Γ = -h_q` forces `a_α = 0` in every sector.
* `hq_flip_all`, `norm_evo_flip_all` (`rem:dimersigns`(i)): the global reversal `h_{-q} = -Γ h_q Γ` is harmless,
  `|e^{-ith_{-q}}_{xy}| = |e^{-ith_q}_{xy}|`; `trace_hq_dimer_example`, `not_similar_dimer_example`: reversing the
  charges of one dimer out of two is not a unitary equivalence (`Tr h_{(1,1,1,1)} = 8h ≠ 0 = Tr h_{(1,1,-1,-1)}`).
* `antidimer_commutator`, `antidimer_comm_ne` (`rem:dimersigns`(ii)): the anti-dimer transfer matrices
  `T_a T_b`, `T_b T_a` (`a = -λ - E`, `b = λ - E`) satisfy
  `[T_a T_b, T_b T_a] = [[b² - a², (a-b)ab], [-(a-b)ab, a² - b²]]` and never commute for `λ ≠ 0`.
-/

open Matrix Finset Filter Topology

noncomputable section

namespace DFL

variable {N : ℕ}

namespace CleanSectorAux

/-! ### The kernel of the path adjacency matrix -/

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

/-- The odd entries of a sequence annihilated by the Dirichlet path adjacency on `[0, M)` vanish. -/
lemma ker_odd (f : ℕ → ℝ) (M : ℕ) (hf : ∀ i, M ≤ i → f i = 0)
    (heq : ∀ i < M, prev0 f i + f (i + 1) = 0) (k : ℕ) : f (2 * k + 1) = 0 := by
  induction k with
  | zero =>
    by_cases hM : 0 < M
    · have e := heq 0 hM
      have hp : prev0 f 0 = 0 := rfl
      rw [hp, zero_add] at e
      exact e
    · exact hf 1 (by omega)
  | succ k ih =>
    by_cases hM : 2 * k + 2 < M
    · have e := heq (2 * k + 2) hM
      have hp : prev0 f (2 * k + 2) = f (2 * k + 1) := rfl
      rw [hp, ih, zero_add] at e
      rw [show 2 * (k + 1) + 1 = 2 * k + 2 + 1 by ring]
      exact e
    · exact hf _ (by omega)

/-- On a chain of even length `2m`, the even entries of a kernel sequence vanish. -/
lemma ker_even (f : ℕ → ℝ) (m : ℕ) (hf : ∀ i, 2 * m ≤ i → f i = 0)
    (heq : ∀ i < 2 * m, prev0 f i + f (i + 1) = 0) (k : ℕ) :
    ∀ j, j + k + 1 = m → f (2 * j) = 0 := by
  induction k with
  | zero =>
    intro j hj
    have e := heq (2 * j + 1) (by omega)
    have hp : prev0 f (2 * j + 1) = f (2 * j) := rfl
    rw [hp, hf (2 * j + 1 + 1) (by omega), add_zero] at e
    exact e
  | succ k ih =>
    intro j hj
    have e := heq (2 * j + 1) (by omega)
    have hp : prev0 f (2 * j + 1) = f (2 * j) := rfl
    have h2 : f (2 * j + 1 + 1) = 0 := by
      rw [show 2 * j + 1 + 1 = 2 * (j + 1) by ring]
      exact ih (j + 1) (by omega)
    rw [hp, h2, add_zero] at e
    exact e

/-- For even `N`, the path adjacency matrix has trivial kernel. -/
lemma pathAdj_ker_even (hN : Even N) (φ : Fin N → ℝ) (hφ : pathAdj N *ᵥ φ = 0) :
    φ = 0 := by
  obtain ⟨m, hm⟩ := hN
  have hf : ∀ i, 2 * m ≤ i → ext0 φ i = 0 := fun i hi => ext0_of_le φ (by omega)
  have heq : ∀ i < 2 * m, prev0 (ext0 φ) i + ext0 φ (i + 1) = 0 := by
    intro i hi
    have e := congrFun hφ ⟨i, by omega⟩
    rw [pathAdj_mulVec_apply] at e
    exact e
  funext x
  rw [← ext0_val φ x]
  obtain ⟨j, hj | hj⟩ := Nat.even_or_odd' x.val
  · rw [hj]
    by_cases hjm : j < m
    · exact ker_even _ m hf heq (m - 1 - j) j (by omega)
    · exact hf _ (by omega)
  · rw [hj]
    exact ker_odd _ (2 * m) hf heq j

/-- For even `N` and `J ≠ 0`, an eigenvector of `h_q` with `2h q_x = E` for all `x` vanishes. -/
lemma eigvec_zero_of_const (J h : ℝ) (hJ : J ≠ 0) (hN : Even N) (q : Cfg N) (E : ℝ)
    (φ : Fin N → ℝ) (hφ : hq J h q *ᵥ φ = E • φ) (hc : ∀ x, 2 * h * sgn (q x) = E) :
    φ = 0 := by
  apply pathAdj_ker_even hN
  funext x
  have e := congrFun hφ x
  simp only [hq, Matrix.add_mulVec, Matrix.smul_mulVec, Pi.add_apply, Pi.smul_apply, smul_eq_mul,
    Matrix.mulVec_diagonal] at e
  have e2 : -J * (pathAdj N *ᵥ φ) x = 0 := by linear_combination e - (φ x) * hc x
  rcases mul_eq_zero.mp e2 with h' | h'
  · exact absurd (neg_eq_zero.mp h') hJ
  · simpa using h'

/-! ### Reflection symmetry of the clean sector -/

/-- The clean Hamiltonian `h_{1_N}` is invariant under the reflection `x ↦ N-1-x`. -/
lemma hq_zero_rev (J h : ℝ) (x y : Fin N) :
    hq J h (0 : Cfg N) (Fin.rev x) (Fin.rev y) = hq J h (0 : Cfg N) x y := by
  have hx := x.isLt
  have hy := y.isLt
  have e : ((Fin.rev x).val + 1 = (Fin.rev y).val ∨ (Fin.rev y).val + 1 = (Fin.rev x).val) ↔
      (x.val + 1 = y.val ∨ y.val + 1 = x.val) := by
    simp only [Fin.val_rev]; omega
  simp only [hq, Matrix.add_apply, Matrix.smul_apply, pathAdj, Matrix.diagonal_apply, Fin.rev_inj,
    Pi.zero_apply, e]

/-- In the clean sector (`J ≠ 0`), the eigenvectors satisfy `|φ_α(N-1-x)| = |φ_α(x)|`. -/
lemma abs_eigvec_rev (J h : ℝ) (hJ : J ≠ 0) (α x : Fin N) :
    |(hq_isHermitian J h (0 : Cfg N)).eigenvectorBasis α (Fin.rev x)| =
      |(hq_isHermitian J h (0 : Cfg N)).eigenvectorBasis α x| := by
  set hA := hq_isHermitian J h (0 : Cfg N) with hAdef
  have hv : ∀ x, ∑ y, hq J h 0 x y * hA.eigenvectorBasis α (Fin.rev y) =
      hA.eigenvalues α * hA.eigenvectorBasis α (Fin.rev x) := by
    intro x
    rw [← NormalizationAux.eig_eq _ hA α (Fin.rev x)]
    refine Fintype.sum_equiv Fin.revPerm _ _ fun y => ?_
    rw [Fin.revPerm_apply, hq_zero_rev]
  have hn : ∑ x, (hA.eigenvectorBasis α (Fin.rev x)) ^ 2 = 1 := by
    rw [← ImbalanceAux.sum_eig_sq' _ hA α]
    exact Fintype.sum_equiv Fin.revPerm _ _ fun y => by rw [Fin.revPerm_apply]
  obtain ⟨β, hβ1, hβ2⟩ := NormalizationAux.exists_abs_eq_eigvec _ hA
    (hq_eigenvalues_injective J h hJ 0) _ _ hv hn
  have hβα : β = α := hq_eigenvalues_injective J h hJ 0 hβ1
  subst hβα
  exact hβ2 x

/-! ### Conjugation of the propagator by `Γ` -/

/-- For a real symmetric `A`, `(e^{itA})_{xy}` is the complex conjugate of `(e^{-itA})_{xy}`. -/
lemma evo_neg_apply (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (t : ℝ) (x y : Fin N) :
    evo A (-t) x y = starRingEnd ℂ (evo A t x y) := by
  rw [OneBodyAux.evo_apply_eq_sum A hA (-t), OneBodyAux.evo_apply_eq_sum A hA t, map_sum]
  refine Finset.sum_congr rfl fun α _ => ?_
  rw [map_mul, map_mul, Complex.conj_ofReal, Complex.conj_ofReal, ← Complex.exp_conj]
  congr 2
  simp only [map_mul, map_neg, Complex.conj_ofReal, Complex.conj_I]
  push_cast
  ring

/-- `cplx Γ` is the diagonal matrix of the signs `(-1)^x`. -/
lemma cplx_stag : cplx (stag N) = diagonal (fun x => (((-1 : ℝ) ^ x.val : ℝ) : ℂ)) := by
  unfold cplx stag
  rw [Matrix.diagonal_map (by simp)]

/-- `cplx Γ * cplx Γ = 1`. -/
lemma cplx_stag_mul_self : cplx (stag N) * cplx (stag N) = 1 := by
  have hc : ∀ M : Matrix (Fin N) (Fin N) ℝ, cplx M = (Complex.ofRealHom).mapMatrix M :=
    fun M => rfl
  rw [hc, ← map_mul, stag_mul_stag, map_one]

/-- `e^{-it(-ΓAΓ)} = Γ e^{itA} Γ`. -/
lemma evo_neg_conj_stag (A : Matrix (Fin N) (Fin N) ℝ) (t : ℝ) :
    evo (-(stag N * A * stag N)) t = cplx (stag N) * evo A (-t) * cplx (stag N) := by
  have hc : ∀ M : Matrix (Fin N) (Fin N) ℝ, cplx M = (Complex.ofRealHom).mapMatrix M :=
    fun M => rfl
  have hinv : (cplx (stag N))⁻¹ = cplx (stag N) := Matrix.inv_eq_right_inv cplx_stag_mul_self
  have hunit : IsUnit (cplx (stag N)) :=
    (Matrix.isUnit_iff_isUnit_det _).mpr (Matrix.isUnit_det_of_right_inverse cplx_stag_mul_self)
  have hneg : cplx (-(stag N * A * stag N)) = -(cplx (stag N) * cplx A * cplx (stag N)) := by
    simp only [hc, map_neg, map_mul]
  have hgen : (-(t : ℂ) * Complex.I) • cplx (-(stag N * A * stag N)) =
      cplx (stag N) * ((-((-t : ℝ) : ℂ) * Complex.I) • cplx A) * (cplx (stag N))⁻¹ := by
    rw [hinv, hneg, mul_smul_comm, smul_mul_assoc, smul_neg, ← neg_smul]
    congr 1
    push_cast
    ring
  unfold evo
  rw [hgen, Matrix.exp_conj _ _ hunit, hinv]

end CleanSectorAux

/-! ### The clean sector (paper: `rem:gaussproj`) -/

/-- **No decay of the eigenfunction correlator in the clean sector** (paper: `rem:gaussproj`(i)):
`Q_{1_N}(1, N) = 1`. -/
theorem corrQ_clean_ends (J h : ℝ) (hJ : J ≠ 0) (hN : 0 < N) :
    corrQ (hq J h (0 : Cfg N)) (hq_isHermitian J h 0) ⟨0, hN⟩ ⟨N - 1, by omega⟩ = 1 := by
  have hlast : (⟨N - 1, by omega⟩ : Fin N) = Fin.rev ⟨0, hN⟩ := Fin.ext (by simp)
  rw [corrQ, hlast]
  simp_rw [CleanSectorAux.abs_eigvec_rev J h hJ]
  rw [← OneBodyAux.sum_eigvec_sq _ (hq_isHermitian J h 0) ⟨0, hN⟩]
  refine Finset.sum_congr rfl fun α _ => ?_
  rw [abs_mul_abs_self, sq]

/-- **No memory of a density wave in the clean sector** (paper: `rem:gaussproj`(ii)): for even `N`, every
eigenvector of `h_{1_N}` has staggered weight `a_α = 0`. -/
theorem stagWeight_clean (J h : ℝ) (hJ : J ≠ 0) (hN : Even N) (α : Fin N) :
    stagWeight J h (0 : Cfg N) α = 0 := by
  have key := eigenvalue_mul_stagWeight J h (0 : Cfg N) α
  have hs : ∑ x, sgn ((0 : Cfg N) x) * (-1 : ℝ) ^ x.val *
      ((hq_isHermitian J h (0 : Cfg N)).eigenvectorBasis α x) ^ 2 =
        stagWeight J h (0 : Cfg N) α := by
    simp [stagWeight, sgn]
  rw [hs] at key
  by_contra hne
  have hE : (hq_isHermitian J h (0 : Cfg N)).eigenvalues α = 2 * h := mul_right_cancel₀ hne key
  have h0 := CleanSectorAux.eigvec_zero_of_const J h hJ hN 0 _ _
    ((hq_isHermitian J h (0 : Cfg N)).mulVec_eigenvectorBasis α) (fun x => by simp [sgn, hE])
  have h0' : ∀ x, (hq_isHermitian J h (0 : Cfg N)).eigenvectorBasis α x = 0 :=
    fun x => congrFun h0 x
  have hn := ImbalanceAux.sum_eig_sq' _ (hq_isHermitian J h (0 : Cfg N)) α
  simp [h0'] at hn

/-- The infinite-time imbalance of the clean sector is `N⁻¹ ∑_α a_α² = 0` (paper: `rem:gaussproj`(ii)). -/
theorem imbInf_clean (J h : ℝ) (hJ : J ≠ 0) (hN : Even N) : imbInf J h (0 : Cfg N) = 0 := by
  simp [imbInf, stagWeight_clean J h hJ hN]

/-- The time-averaged imbalance of the clean sector tends to `0` (paper: `rem:gaussproj`(ii)). -/
theorem imbOne_clean_timeAvg (J h : ℝ) (hJ : J ≠ 0) (hN : Even N) :
    Tendsto (fun T : ℝ => (1 / T) * ∫ t in (0 : ℝ)..T, imbOne J h (0 : Cfg N) t) atTop (𝓝 0) := by
  have := imbOne_timeAvg J h hJ (0 : Cfg N)
  rwa [imbInf_clean J h hJ hN] at this

/-- **Chiral symmetry at `h = 0`** (paper: `rem:main-imbalance`(ii)): for even `N`, `a_α = 0` in every
sector. -/
theorem stagWeight_h_zero (J : ℝ) (hJ : J ≠ 0) (hN : Even N) (q : Cfg N) (α : Fin N) :
    stagWeight J 0 q α = 0 := by
  have key := eigenvalue_mul_stagWeight J 0 q α
  rw [mul_zero, zero_mul] at key
  by_contra hne
  have hE : (hq_isHermitian J 0 q).eigenvalues α = 0 := by
    rcases mul_eq_zero.mp key with h' | h'
    · exact h'
    · exact absurd h' hne
  have h0 := CleanSectorAux.eigvec_zero_of_const J 0 hJ hN q _ _
    ((hq_isHermitian J 0 q).mulVec_eigenvectorBasis α) (fun x => by rw [hE]; ring)
  have h0' : ∀ x, (hq_isHermitian J 0 q).eigenvectorBasis α x = 0 := fun x => congrFun h0 x
  have hn := ImbalanceAux.sum_eig_sq' _ (hq_isHermitian J 0 q) α
  simp [h0'] at hn

/-! ### Deterministic signs (paper: `rem:dimersigns`) -/

/-- The global reversal of the charges: `h_{-q} = -Γ h_q Γ` (paper: `rem:dimersigns`(i)). -/
theorem hq_flip_all (J h : ℝ) (q : Cfg N) :
    hq J h (fun x => q x + 1) = -(stag N * hq J h q * stag N) := by
  rw [stag_conj_hq]
  ext x y
  simp only [hq, Matrix.add_apply, Matrix.smul_apply, Matrix.neg_apply, smul_eq_mul,
    Matrix.diagonal_apply, NormalizationAux.sgn_add_one]
  split_ifs <;> ring

/-- `|e^{-ith_{-q}}_{xy}| = |e^{-ith_q}_{xy}|` (paper: `rem:dimersigns`(i)). -/
theorem norm_evo_flip_all (J h : ℝ) (q : Cfg N) (t : ℝ) (x y : Fin N) :
    ‖evo (hq J h (fun x => q x + 1)) t x y‖ = ‖evo (hq J h q) t x y‖ := by
  rw [hq_flip_all J h q, CleanSectorAux.evo_neg_conj_stag, CleanSectorAux.cplx_stag,
    Matrix.mul_diagonal, Matrix.diagonal_mul,
    CleanSectorAux.evo_neg_apply _ (hq_isHermitian J h q), norm_mul, norm_mul, Complex.norm_conj]
  simp

/-- `Tr h_{(1,1,1,1)} = 8h` and `Tr h_{(1,1,-1,-1)} = 0` (paper: `rem:dimersigns`(i)). -/
theorem trace_hq_dimer_example (J h : ℝ) :
    trace (hq J h (0 : Cfg 4)) = 8 * h ∧ trace (hq J h (![0, 0, 1, 1] : Cfg 4)) = 0 := by
  have h1 : (1 : ZMod 2) ≠ 0 := by decide
  constructor
  · simp [hq, Matrix.trace, sgn, pathAdj]
    ring
  · simp [hq, Matrix.trace, Fin.sum_univ_four, sgn, pathAdj, h1]

/-- Reversing the charges of one dimer out of two is not a similarity transformation for `h ≠ 0`
(paper: `rem:dimersigns`(i)). -/
theorem not_similar_dimer_example (J h : ℝ) (hh : h ≠ 0) (U : Matrix (Fin 4) (Fin 4) ℝ)
    (hU : IsUnit U) : U * hq J h (0 : Cfg 4) * U⁻¹ ≠ hq J h (![0, 0, 1, 1] : Cfg 4) := by
  intro heq
  have ht := congrArg Matrix.trace heq
  rw [Matrix.trace_mul_cycle, Matrix.nonsing_inv_mul U ((Matrix.isUnit_iff_isUnit_det U).mp hU),
    Matrix.one_mul, (trace_hq_dimer_example J h).1, (trace_hq_dimer_example J h).2] at ht
  exact hh (by linarith)

/-- The commutator of the anti-dimer transfer matrices (paper: `rem:dimersigns`(ii)). -/
theorem antidimer_commutator (a b : ℂ) :
    tmC a * tmC b * (tmC b * tmC a) - tmC b * tmC a * (tmC a * tmC b) =
      !![b ^ 2 - a ^ 2, (a - b) * a * b; -((a - b) * a * b), a ^ 2 - b ^ 2] := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [tmC] <;> ring

/-- **The anti-dimer law has no critical energy** (paper: `rem:dimersigns`(ii)): for `λ ≠ 0` the two
anti-dimer transfer matrices never commute. -/
theorem antidimer_comm_ne (lam : ℝ) (hl : lam ≠ 0) (E : ℂ) :
    tmC (-lam - E) * tmC (lam - E) * (tmC (lam - E) * tmC (-lam - E)) ≠
      tmC (lam - E) * tmC (-lam - E) * (tmC (-lam - E) * tmC (lam - E)) := by
  intro heq
  have hc := antidimer_commutator (-lam - E) (lam - E)
  rw [heq, sub_self] at hc
  have e00 := congrFun (congrFun hc 0) 0
  have e01 := congrFun (congrFun hc 0) 1
  simp only [Matrix.zero_apply, Matrix.of_apply, Matrix.cons_val', Matrix.cons_val_zero,
    Matrix.cons_val_one, Matrix.empty_val', Matrix.cons_val_fin_one] at e00 e01
  have hl' : (lam : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hl
  have hE : E = 0 := by
    have hmul : (lam : ℂ) * E = 0 := by linear_combination (1 / 4 : ℂ) * e00
    rcases mul_eq_zero.mp hmul with h' | h'
    · exact absurd h' hl'
    · exact h'
  subst hE
  have h3 : (lam : ℂ) ^ 3 = 0 := by linear_combination (-1 / 2 : ℂ) * e01
  exact hl' (pow_eq_zero_iff (by norm_num) |>.mp h3)

end DFL
