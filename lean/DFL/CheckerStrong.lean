import DFL.Cluster
import DFL.TwoDLaws
import DFL.TwoDLattices

/-!
# The checkerboard state at strong coupling (paper: Remark `rem:checkerstrong`)

For the checkerboard charges `q⁰ = -ε` of a bipartite graph (sublattice signs `ε`, every edge joins
opposite signs):

* every equal-charge cluster is a single site (`checker_inCluster_iff`), `h⁰_{q⁰,s} = 2h diag(q⁰)`
  (`hqs0_checker`), and Proposition `prop:cluster`(b) gives the deterministic return amplitude
  `|⟨δ_x, e^{-ith_{q⁰,s}} δ_x⟩| ≥ 1 - 8J/g - 16|t|J²/g` in every flux sector (`checker_return`);
* `R_s = (K_s² + 4h²)^{1/2}` (`checkerR`, continuous functional calculus of the Hermitian matrix
  `K_s`) satisfies `R_s² = K_s² + 4h²` (`checkerR_sq`) and, by the scalar inequality
  `0 ≤ 1 + y/2 - √(1+y) ≤ y²/8` (`one_add_half_sub_sqrt_bounds`) and `‖K_s‖ ≤ 4J`
  (`norm_cplx_Ks_le`),
  `‖R_s - 2|h| - K_s²/(4|h|)‖ ≤ ‖K_s‖⁴/(64|h|³) ≤ 4J⁴/|h|³` (`norm_checkerR_sub_le`,
  `norm_checkerR_sub_le_J`); hence `‖e^{-itR_s} - e^{-it(2|h| + K_s²/(4|h|))}‖ ≤ |t| ‖K_s‖⁴/(64|h|³)`
  (`norm_exp_checkerR_sub_le`, `norm_exp_checkerR_sub_le_J`);
* `K_s²` connects only sites of the same sublattice (`Ks_sq_sublattice`); between the two ends of a
  diagonal of a plaquette it equals `J² s s (1 + W_p)` (`Ks_sq_plaquette`), between the ends of a
  straight path of length two `J² s s` (`Ks_sq_line`), and it vanishes between distinct sites without
  a common neighbour (`Ks_sq_far`). On the planar box (`box_simple`, `box_deg_le_four`) this gives
  the effective sublattice hopping of the remark: `±2J²` (resp. `0`) on the plaquette diagonals if
  `W_p = +1` (resp. `-1`), and `±J²` on the straight next-nearest-neighbour bonds
  (`box_Ks_sq_diag`, `box_Ks_sq_antidiag`, `box_Ks_sq_horiz`, `box_Ks_sq_vert`).
-/

open Matrix Finset NormedSpace
open scoped Matrix.Norms.L2Operator

noncomputable section

namespace DFL

namespace CheckerStrongAux

variable {N M : ℕ} (G : FinGraph N M)

/-! ### Auxiliary lemmas: graphs and `K_s` -/

/-- Entries of `h_{q,s}`. -/
lemma hqs_apply (J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) (i j : Fin N) :
    G.hqs J h q s i j = G.Ks J s i j + 2 * h * (if i = j then q i else 0) := by
  simp only [FinGraph.hqs, Matrix.add_apply, Matrix.smul_apply, smul_eq_mul,
    Matrix.diagonal_apply]

/-- `K_s` vanishes between vertices that are not joined by an edge. -/
lemma Ks_eq_zero_of_not_adj (J : ℝ) (s : Fin M → ℝ) (i j : Fin N)
    (hij : ¬ ∃ e, (G.src e = i ∧ G.tgt e = j) ∨ (G.tgt e = i ∧ G.src e = j)) :
    G.Ks J s i j = 0 := by
  show -J * ∑ e ∈ univ.filter (fun e => (G.src e = i ∧ G.tgt e = j) ∨
    (G.tgt e = i ∧ G.src e = j)), s e = 0
  rw [Finset.filter_false_of_mem, Finset.sum_empty, mul_zero]
  intro e _ he
  exact hij ⟨e, he⟩

/-- On a bipartite graph, the two ends of an edge have opposite sublattice signs. -/
lemma eps_of_adj (ε : Fin N → ℝ) (hε : ∀ e, ε (G.src e) = -ε (G.tgt e)) {i j : Fin N}
    (hij : ∃ e, (G.src e = i ∧ G.tgt e = j) ∨ (G.tgt e = i ∧ G.src e = j)) : ε i = -ε j := by
  obtain ⟨e, ⟨h1, h2⟩ | ⟨h1, h2⟩⟩ := hij
  · rw [← h1, ← h2]
    exact hε e
  · have := hε e
    rw [h1, h2] at this
    linarith

/-- On a bipartite graph, no edge joins two vertices of the same sublattice. -/
lemma not_adj_of_eps_eq (ε : Fin N → ℝ) (hε1 : ∀ i, ε i = 1 ∨ ε i = -1)
    (hε : ∀ e, ε (G.src e) = -ε (G.tgt e)) {i j : Fin N} (hij : ε i = ε j) :
    ¬ ∃ e, (G.src e = i ∧ G.tgt e = j) ∨ (G.tgt e = i ∧ G.src e = j) := by
  intro h
  have h1 := eps_of_adj G ε hε h
  rcases hε1 i with h2 | h2 <;> linarith

/-- On a simple graph, the only edge joining the two ends of an edge `e₀` is `e₀`. -/
lemma filter_join_eq_singleton
    (hsimple : ∀ e e', ((G.src e = G.src e' ∧ G.tgt e = G.tgt e') ∨
      (G.src e = G.tgt e' ∧ G.tgt e = G.src e')) → e = e')
    {i j : Fin N} {e₀ : Fin M}
    (h₀ : (G.src e₀ = i ∧ G.tgt e₀ = j) ∨ (G.tgt e₀ = i ∧ G.src e₀ = j)) :
    univ.filter (fun e => (G.src e = i ∧ G.tgt e = j) ∨ (G.tgt e = i ∧ G.src e = j)) = {e₀} := by
  ext e
  simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_singleton]
  constructor
  · intro he
    apply hsimple
    rcases he with ⟨h1, h2⟩ | ⟨h1, h2⟩ <;> rcases h₀ with ⟨h3, h4⟩ | ⟨h3, h4⟩
    · exact Or.inl ⟨h1.trans h3.symm, h2.trans h4.symm⟩
    · exact Or.inr ⟨h1.trans h3.symm, h2.trans h4.symm⟩
    · exact Or.inr ⟨h2.trans h4.symm, h1.trans h3.symm⟩
    · exact Or.inl ⟨h2.trans h4.symm, h1.trans h3.symm⟩
  · rintro rfl
    exact h₀

/-- On a simple graph, `(K_s)_{ij} = -J s_{e₀}` if the edge `e₀` joins `i` and `j`. -/
lemma Ks_of_edge (J : ℝ) (s : Fin M → ℝ)
    (hsimple : ∀ e e', ((G.src e = G.src e' ∧ G.tgt e = G.tgt e') ∨
      (G.src e = G.tgt e' ∧ G.tgt e = G.src e')) → e = e')
    {i j : Fin N} {e₀ : Fin M}
    (h₀ : (G.src e₀ = i ∧ G.tgt e₀ = j) ∨ (G.tgt e₀ = i ∧ G.src e₀ = j)) :
    G.Ks J s i j = -J * s e₀ := by
  show -J * ∑ e ∈ univ.filter (fun e => (G.src e = i ∧ G.tgt e = j) ∨
    (G.tgt e = i ∧ G.src e = j)), s e = -J * s e₀
  rw [filter_join_eq_singleton G hsimple h₀, Finset.sum_singleton]

/-- A term `(K_s)_{ij} (K_s)_{jk}` vanishes unless `j` is a common neighbour of `i` and `k`. -/
lemma Ks_mul_Ks_eq_zero (J : ℝ) (s : Fin M → ℝ) (i j k : Fin N)
    (h : ¬ ((∃ e, (G.src e = i ∧ G.tgt e = j) ∨ (G.tgt e = i ∧ G.src e = j)) ∧
      (∃ e, (G.src e = j ∧ G.tgt e = k) ∨ (G.tgt e = j ∧ G.src e = k)))) :
    G.Ks J s i j * G.Ks J s j k = 0 := by
  by_cases h1 : ∃ e, (G.src e = i ∧ G.tgt e = j) ∨ (G.tgt e = i ∧ G.src e = j)
  · rw [Ks_eq_zero_of_not_adj G J s j k (fun h2 => h ⟨h1, h2⟩), mul_zero]
  · rw [Ks_eq_zero_of_not_adj G J s i j h1, zero_mul]

/-! ### Auxiliary lemmas: scalar inequalities -/

/-- `0 ≤ 1 + y/2 - √(1+y) ≤ y²/8` for `y ≥ 0`. -/
lemma sqrt_bounds (y : ℝ) (hy : 0 ≤ y) :
    0 ≤ 1 + y / 2 - √(1 + y) ∧ 1 + y / 2 - √(1 + y) ≤ y ^ 2 / 8 := by
  have hs : √(1 + y) ^ 2 = 1 + y := Real.sq_sqrt (by linarith)
  have hs0 : 0 ≤ √(1 + y) := Real.sqrt_nonneg _
  have hs1 : 1 ≤ √(1 + y) := by nlinarith
  constructor
  · nlinarith [sq_nonneg (√(1 + y) - 1)]
  · have hy' : y = √(1 + y) ^ 2 - 1 := by linarith
    have key : 0 ≤ (√(1 + y) - 1) ^ 3 * (√(1 + y) + 3) :=
      mul_nonneg (pow_nonneg (by linarith) 3) (by linarith)
    have e : y ^ 2 = (√(1 + y) ^ 2 - 1) ^ 2 := by rw [← hy']
    nlinarith [key, e]

/-! ### Auxiliary lemmas: matrices -/

/-- An entry of a matrix is bounded by its operator norm. -/
lemma norm_entry_le {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℂ) (i j : ι) :
    ‖A i j‖ ≤ ‖A‖ := by
  have h1 := ClusterAux.col_sum_le A j
  have h2 : ‖A i j‖ ^ 2 ≤ ∑ y, ‖A y j‖ ^ 2 :=
    Finset.single_le_sum (f := fun y => ‖A y j‖ ^ 2) (fun y _ => by positivity)
      (Finset.mem_univ i)
  exact (sq_le_sq₀ (norm_nonneg _) (norm_nonneg _)).mp (h2.trans h1)

/-- The propagator of a real diagonal generator has unimodular diagonal entries. -/
lemma norm_evo_diagonal {ι : Type*} [Fintype ι] [DecidableEq ι] (d : ι → ℝ) (t : ℝ) (x : ι) :
    ‖evo (diagonal d) t x x‖ = 1 := by
  have hc : cplx (diagonal d) = diagonal (fun i => (d i : ℂ)) := by
    ext i j
    by_cases hij : i = j <;> simp [cplx, hij]
  unfold evo
  rw [hc, ← diagonal_smul, Matrix.exp_diagonal, diagonal_apply_eq, Pi.coe_exp, Pi.smul_apply,
    smul_eq_mul, ← Complex.exp_eq_exp_ℂ, Complex.norm_exp]
  simp

/-- Complexification of a real multiple of a matrix. -/
lemma cplx_smul {ι : Type*} (c : ℝ) (A : Matrix ι ι ℝ) : cplx (c • A) = (c : ℂ) • cplx A := by
  ext i j
  simp [cplx]

/-- `algebraMap` of a real scalar into the complex matrices. -/
lemma algebraMap_eq {ι : Type*} [Fintype ι] [DecidableEq ι] (c : ℝ) :
    algebraMap ℝ (Matrix ι ι ℂ) c = (c : ℂ) • (1 : Matrix ι ι ℂ) := by
  rw [Algebra.algebraMap_eq_smul_one, Complex.coe_smul]

/-- The complexified `K_s` is Hermitian. -/
lemma cplx_Ks_isHermitian (J : ℝ) (s : Fin M → ℝ) : (cplx (G.Ks J s)).IsHermitian :=
  TwoDAux.conjTranspose_map_ofReal _ (TwoDLawsAux.Ks_transpose G J s)

/-- The real spectrum of a self-adjoint matrix lies in `[-‖A‖, ‖A‖]`. -/
lemma abs_le_norm_of_mem_spectrum {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℂ)
    (hA : IsSelfAdjoint A) {x : ℝ} (hx : x ∈ spectrum ℝ A) : |x| ≤ ‖A‖ := by
  have := norm_apply_le_norm_cfc (fun x : ℝ => x) A hx (by fun_prop) hA
  rwa [cfc_id' ℝ A hA, Real.norm_eq_abs] at this

/-- `|√(x² + 4h²) - 2|h| - x²/(4|h|)| ≤ k⁴/(64|h|³)` for `|x| ≤ k` (from `sqrt_bounds` with
`y = x²/(4h²)`). -/
lemma sqrt_sub_le (x k h : ℝ) (hh : h ≠ 0) (hx : |x| ≤ k) :
    ‖√(x ^ 2 + 4 * h ^ 2) - (2 * |h| + 1 / (4 * |h|) * x ^ 2)‖ ≤ k ^ 4 / (64 * |h| ^ 3) := by
  have ha0 : 0 < |h| := abs_pos.mpr hh
  have hh2 : h ^ 2 = |h| ^ 2 := (sq_abs h).symm
  obtain ⟨b1, b2⟩ := sqrt_bounds (x ^ 2 / (4 * |h| ^ 2)) (by positivity)
  have hsq : √(x ^ 2 + 4 * h ^ 2) = 2 * |h| * √(1 + x ^ 2 / (4 * |h| ^ 2)) := by
    rw [show x ^ 2 + 4 * h ^ 2 = (2 * |h|) ^ 2 * (1 + x ^ 2 / (4 * |h| ^ 2)) by
      rw [hh2]; field_simp; ring]
    rw [Real.sqrt_mul (by positivity), Real.sqrt_sq (by positivity)]
  have e1 : 2 * |h| + 1 / (4 * |h|) * x ^ 2 = 2 * |h| * (1 + x ^ 2 / (4 * |h| ^ 2) / 2) := by
    field_simp
  rw [hsq, e1, Real.norm_eq_abs, abs_sub_comm, ← mul_sub, abs_mul,
    abs_of_pos (by positivity : 0 < 2 * |h|), abs_of_nonneg b1]
  have hx4 : x ^ 4 ≤ k ^ 4 := by
    have := pow_le_pow_left₀ (abs_nonneg x) hx 4
    rwa [pow_abs, abs_of_nonneg (by positivity : (0 : ℝ) ≤ x ^ 4)] at this
  calc 2 * |h| * (1 + x ^ 2 / (4 * |h| ^ 2) / 2 - √(1 + x ^ 2 / (4 * |h| ^ 2)))
      ≤ 2 * |h| * ((x ^ 2 / (4 * |h| ^ 2)) ^ 2 / 8) := by gcongr
    _ = x ^ 4 / (64 * |h| ^ 3) := by field_simp; ring
    _ ≤ k ^ 4 / (64 * |h| ^ 3) := by gcongr

/-- Spectral form `U diag(f(λ)) U*` of the continuous functional calculus of a Hermitian
matrix. -/
lemma cfc_eq_conj {ι : Type*} [Fintype ι] [DecidableEq ι] {A : Matrix ι ι ℂ}
    (hA : A.IsHermitian) (f : ℝ → ℝ) :
    cfc f A = (hA.eigenvectorUnitary : Matrix ι ι ℂ) *
      diagonal (fun i => ((f (hA.eigenvalues i) : ℝ) : ℂ)) *
        star (hA.eigenvectorUnitary : Matrix ι ι ℂ) := by
  rw [hA.cfc_eq, Matrix.IsHermitian.cfc, Unitary.conjStarAlgAut_apply]
  rfl

/-- The exponential commutes with unitary conjugation. -/
lemma exp_unitary_conj {ι : Type*} [Fintype ι] [DecidableEq ι] (U : unitary (Matrix ι ι ℂ))
    (D : Matrix ι ι ℂ) (c : ℂ) :
    exp (c • ((U : Matrix ι ι ℂ) * D * star (U : Matrix ι ι ℂ))) =
      (U : Matrix ι ι ℂ) * exp (c • D) * star (U : Matrix ι ι ℂ) := by
  have hinv : (U : Matrix ι ι ℂ)⁻¹ = star (U : Matrix ι ι ℂ) :=
    Matrix.inv_eq_left_inv (Unitary.coe_star_mul_self U)
  rw [← hinv, ← Matrix.exp_conj _ _ Unitary.isUnit_coe, Matrix.mul_smul, Matrix.smul_mul]

/-- The exponential of a multiple of a diagonal matrix. -/
lemma exp_smul_diagonal {ι : Type*} [Fintype ι] [DecidableEq ι] (c : ℂ) (v : ι → ℂ) :
    exp (c • diagonal v) = diagonal (fun i => exp (c * v i)) := by
  rw [← diagonal_smul, Matrix.exp_diagonal]
  congr 1
  funext i
  rw [Pi.coe_exp]
  rfl

/-- Unitary conjugation preserves the operator norm. -/
lemma norm_unitary_conj {ι : Type*} [Fintype ι] [DecidableEq ι] (U : unitary (Matrix ι ι ℂ))
    (X : Matrix ι ι ℂ) : ‖(U : Matrix ι ι ℂ) * X * star (U : Matrix ι ι ℂ)‖ = ‖X‖ := by
  rw [← Unitary.coe_star, CStarRing.norm_mul_coe_unitary, CStarRing.norm_coe_unitary_mul]

/-- `|e^{-itr} - e^{-itb}| ≤ |t| |r - b|` for real `t`, `r`, `b`. -/
lemma norm_exp_sub_exp_le (t r b : ℝ) :
    ‖exp ((-(t : ℂ) * Complex.I) * (r : ℂ)) - exp ((-(t : ℂ) * Complex.I) * (b : ℂ))‖ ≤
      |t| * |r - b| := by
  rw [← Complex.exp_eq_exp_ℂ]
  have e : Complex.exp ((-(t : ℂ) * Complex.I) * (r : ℂ)) -
      Complex.exp ((-(t : ℂ) * Complex.I) * (b : ℂ)) =
      Complex.exp ((-(t : ℂ) * Complex.I) * (b : ℂ)) *
        (Complex.exp (Complex.I * ((-(t * (r - b)) : ℝ) : ℂ)) - 1) := by
    rw [mul_sub, mul_one, ← Complex.exp_add]
    congr 2
    push_cast
    ring
  rw [e, norm_mul, Complex.norm_exp]
  have hre : ((-(t : ℂ) * Complex.I) * (b : ℂ)).re = 0 := by simp
  rw [hre, Real.exp_zero, one_mul]
  refine Real.norm_exp_I_mul_ofReal_sub_one_le.trans ?_
  rw [Real.norm_eq_abs, abs_neg, abs_mul]

end CheckerStrongAux

namespace FinGraph

variable {N M : ℕ} (G : FinGraph N M)

/-! ### Clusters and the return amplitude -/

/-- For the checkerboard charges `q⁰ = -ε` every equal-charge cluster is a single site. -/
theorem checker_inCluster_iff (ε : Fin N → ℝ) (hε1 : ∀ i, ε i = 1 ∨ ε i = -1)
    (hε : ∀ e, ε (G.src e) = -ε (G.tgt e)) (x y : Fin N) :
    G.InCluster (fun i => -ε i) x y ↔ x = y := by
  have hnone : ∀ b, ¬ G.eqChargeAdj (fun i => -ε i) x b := by
    rintro b ⟨hadj, hq⟩
    exact CheckerStrongAux.not_adj_of_eps_eq G ε hε1 hε (neg_inj.mp hq) hadj
  unfold InCluster
  rw [Relation.reflTransGen_iff_eq hnone]
  exact eq_comm

/-- `h⁰_{q⁰,s} = 2h diag(q⁰)` for the checkerboard charges. -/
theorem hqs0_checker (J h : ℝ) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε1 : ∀ i, ε i = 1 ∨ ε i = -1) (hε : ∀ e, ε (G.src e) = -ε (G.tgt e)) :
    G.hqs0 J h (fun i => -ε i) s = (2 * h) • diagonal (fun i => -ε i) := by
  ext i j
  simp only [hqs0, Matrix.smul_apply, smul_eq_mul, Matrix.diagonal_apply]
  by_cases hij : i = j
  · subst hij
    simp [CheckerStrongAux.hqs_apply, ClusterAux.Ks_diag]
  · rw [ite_eq_right hij, mul_zero]
    split_ifs with hq
    · rw [CheckerStrongAux.hqs_apply, ite_eq_right hij, mul_zero, add_zero]
      exact CheckerStrongAux.Ks_eq_zero_of_not_adj G J s i j
        (CheckerStrongAux.not_adj_of_eps_eq G ε hε1 hε (neg_inj.mp hq))
    · rfl

/-- **Return amplitude of the checkerboard state** (paper: Remark `rem:checkerstrong`): in every
flux sector, `|⟨δ_x, e^{-ith_{q⁰,s}} δ_x⟩| ≥ 1 - 8J/g - 16|t|J²/g`, `g = 4|h| - 8J`. -/
theorem checker_return (J h : ℝ) (hJ : 0 ≤ J) (hh : 2 * J < |h|) (s : Fin M → ℝ)
    (hs : ∀ e, s e = 1 ∨ s e = -1) (ε : Fin N → ℝ) (hε1 : ∀ i, ε i = 1 ∨ ε i = -1)
    (hε : ∀ e, ε (G.src e) = -ε (G.tgt e))
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4) (t : ℝ)
    (x : Fin N) :
    1 - 8 * J / (4 * |h| - 8 * J) - 16 * |t| * J ^ 2 / (4 * |h| - 8 * J) ≤
      ‖evo (G.hqs J h (fun i => -ε i) s) t x x‖ := by
  have hq : ∀ i, -ε i = 1 ∨ -ε i = -1 := fun i => by
    rcases hε1 i with h1 | h1
    · right; linarith
    · left; linarith
  have hb := G.cluster_bound J h hJ hh (fun i => -ε i) hq s hs hdeg t
  rw [hqs0_checker G J h s ε hε1 hε, ← diagonal_smul] at hb
  have h0 := CheckerStrongAux.norm_evo_diagonal ((2 * h) • fun i => -ε i) t x
  have h1 := CheckerStrongAux.norm_entry_le (evo (G.hqs J h (fun i => -ε i) s) t -
    evo (diagonal ((2 * h) • fun i => -ε i)) t) x x
  rw [Matrix.sub_apply] at h1
  have h2 := norm_sub_norm_le (evo (diagonal ((2 * h) • fun i => -ε i)) t x x)
    (evo (G.hqs J h (fun i => -ε i) s) t x x)
  rw [norm_sub_rev] at h2
  linarith

/-! ### The slow generator -/

/-- The scalar inequality `0 ≤ 1 + y/2 - √(1+y) ≤ y²/8` for `y ≥ 0`. -/
theorem one_add_half_sub_sqrt_bounds (y : ℝ) (hy : 0 ≤ y) :
    0 ≤ 1 + y / 2 - √(1 + y) ∧ 1 + y / 2 - √(1 + y) ≤ y ^ 2 / 8 :=
  CheckerStrongAux.sqrt_bounds y hy

/-- `‖K_s‖ ≤ 4J` on graphs of degree at most four (Schur test). -/
theorem norm_cplx_Ks_le (J : ℝ) (hJ : 0 ≤ J) (s : Fin M → ℝ) (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4) :
    ‖cplx (G.Ks J s)‖ ≤ 4 * J :=
  ClusterAux.norm_cplx_le G J hJ s hs hdeg (G.Ks J s) (fun _ _ => le_rfl)

/-- `R_s = (K_s² + 4h²)^{1/2}`, by the continuous functional calculus of the Hermitian matrix
`K_s` (complexified). -/
def checkerR (J h : ℝ) (s : Fin M → ℝ) : Matrix (Fin N) (Fin N) ℂ :=
  cfc (fun x : ℝ => √(x ^ 2 + 4 * h ^ 2)) (cplx (G.Ks J s))

/-- The effective slow generator `2|h| + K_s²/(4|h|)`. -/
def checkerEff (J h : ℝ) (s : Fin M → ℝ) : Matrix (Fin N) (Fin N) ℝ :=
  (2 * |h|) • (1 : Matrix (Fin N) (Fin N) ℝ) + (1 / (4 * |h|)) • (G.Ks J s * G.Ks J s)

namespace CheckerStrongAux

/-- The effective generator `2|h| + K_s²/(4|h|)` is the function `x ↦ 2|h| + x²/(4|h|)` of
`K_s`. -/
lemma cplx_checkerEff (J h : ℝ) (s : Fin M → ℝ) :
    cplx (G.checkerEff J h s) =
      cfc (fun x : ℝ => 2 * |h| + 1 / (4 * |h|) * x ^ 2) (cplx (G.Ks J s)) := by
  have hA := (CheckerStrongAux.cplx_Ks_isHermitian G J s).isSelfAdjoint
  rw [cfc_const_add _ _ _ (by fun_prop) hA, cfc_const_mul _ _ _ (by fun_prop), cfc_pow_id _ 2 hA,
    checkerEff, TwoDLawsAux.cplx_add, TwoDLawsAux.cplx_smul_one, CheckerStrongAux.cplx_smul,
    TwoDLawsAux.cplx_mul, CheckerStrongAux.algebraMap_eq, sq]
  simp only [Complex.coe_smul]

end CheckerStrongAux

/-- `R_s` is Hermitian. -/
theorem checkerR_isHermitian (J h : ℝ) (s : Fin M → ℝ) : (G.checkerR J h s).IsHermitian :=
  (cfc_predicate _ _ : IsSelfAdjoint (G.checkerR J h s)).isHermitian

/-- `R_s² = K_s² + 4h²`. -/
theorem checkerR_sq (J h : ℝ) (s : Fin M → ℝ) :
    G.checkerR J h s * G.checkerR J h s =
      cplx (G.Ks J s * G.Ks J s + (4 * h ^ 2) • (1 : Matrix (Fin N) (Fin N) ℝ)) := by
  have hA := (CheckerStrongAux.cplx_Ks_isHermitian G J s).isSelfAdjoint
  unfold checkerR
  rw [← cfc_mul _ _ _ (by fun_prop) (by fun_prop)]
  rw [cfc_congr (g := fun x : ℝ => x ^ 2 + 4 * h ^ 2) (fun x _ => by
    simp only
    exact Real.mul_self_sqrt (by positivity))]
  rw [cfc_add_const _ _ _ (by fun_prop) hA, cfc_pow_id _ 2 hA, TwoDLawsAux.cplx_add,
    TwoDLawsAux.cplx_mul, TwoDLawsAux.cplx_smul_one, CheckerStrongAux.algebraMap_eq, sq]

/-- **Slow generator** (paper: Remark `rem:checkerstrong`):
`‖R_s - 2|h| - K_s²/(4|h|)‖ ≤ ‖K_s‖⁴/(64|h|³)`. -/
theorem norm_checkerR_sub_le (J h : ℝ) (hh : h ≠ 0) (s : Fin M → ℝ) :
    ‖G.checkerR J h s - cplx (G.checkerEff J h s)‖ ≤ ‖cplx (G.Ks J s)‖ ^ 4 / (64 * |h| ^ 3) := by
  have hA := (CheckerStrongAux.cplx_Ks_isHermitian G J s).isSelfAdjoint
  rw [CheckerStrongAux.cplx_checkerEff, checkerR, ← cfc_sub _ _ _ (by fun_prop) (by fun_prop)]
  refine norm_cfc_le (by positivity) fun x hx => ?_
  exact CheckerStrongAux.sqrt_sub_le x _ h hh
    (CheckerStrongAux.abs_le_norm_of_mem_spectrum _ hA hx)

/-- `‖R_s - 2|h| - K_s²/(4|h|)‖ ≤ 4J⁴/|h|³` on graphs of degree at most four. -/
theorem norm_checkerR_sub_le_J (J h : ℝ) (hJ : 0 ≤ J) (hh : h ≠ 0) (s : Fin M → ℝ)
    (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4) :
    ‖G.checkerR J h s - cplx (G.checkerEff J h s)‖ ≤ 4 * J ^ 4 / |h| ^ 3 := by
  have hK := G.norm_cplx_Ks_le J hJ s hs hdeg
  have ha : 0 < |h| := abs_pos.mpr hh
  calc ‖G.checkerR J h s - cplx (G.checkerEff J h s)‖
      ≤ ‖cplx (G.Ks J s)‖ ^ 4 / (64 * |h| ^ 3) := G.norm_checkerR_sub_le J h hh s
    _ ≤ (4 * J) ^ 4 / (64 * |h| ^ 3) := by gcongr
    _ = 4 * J ^ 4 / |h| ^ 3 := by field_simp; ring

/-- **The slow dynamics** (paper: Remark `rem:checkerstrong`): up to `|t| ‖K_s‖⁴/(64|h|³)`, the
dynamics generated by `R_s` is that generated by `2|h| + K_s²/(4|h|)`. -/
theorem norm_exp_checkerR_sub_le (J h : ℝ) (hh : h ≠ 0) (s : Fin M → ℝ) (t : ℝ) :
    ‖exp ((-(t : ℂ) * Complex.I) • G.checkerR J h s) - evo (G.checkerEff J h s) t‖ ≤
      |t| * (‖cplx (G.Ks J s)‖ ^ 4 / (64 * |h| ^ 3)) := by
  have hA := CheckerStrongAux.cplx_Ks_isHermitian G J s
  -- both generators are diagonal in the eigenbasis of `K_s`
  have hR := CheckerStrongAux.cfc_eq_conj hA (fun x : ℝ => √(x ^ 2 + 4 * h ^ 2))
  have hE := CheckerStrongAux.cfc_eq_conj hA (fun x : ℝ => 2 * |h| + 1 / (4 * |h|) * x ^ 2)
  rw [← CheckerStrongAux.cplx_checkerEff] at hE
  rw [← checkerR] at hR
  unfold evo
  rw [hR, hE, CheckerStrongAux.exp_unitary_conj, CheckerStrongAux.exp_unitary_conj,
    CheckerStrongAux.exp_smul_diagonal, CheckerStrongAux.exp_smul_diagonal, ← Matrix.sub_mul,
    ← Matrix.mul_sub, CheckerStrongAux.norm_unitary_conj, diagonal_sub, Matrix.l2_opNorm_diagonal]
  refine (pi_norm_le_iff_of_nonneg (by positivity)).mpr fun i => ?_
  refine (CheckerStrongAux.norm_exp_sub_exp_le t _ _).trans ?_
  gcongr
  rw [← Real.norm_eq_abs]
  exact CheckerStrongAux.sqrt_sub_le _ _ h hh (CheckerStrongAux.abs_le_norm_of_mem_spectrum _
    hA.isSelfAdjoint (hA.eigenvalues_mem_spectrum_real i))

/-- The slow dynamics up to `4|t|J⁴/|h|³` on graphs of degree at most four. -/
theorem norm_exp_checkerR_sub_le_J (J h : ℝ) (hJ : 0 ≤ J) (hh : h ≠ 0) (s : Fin M → ℝ)
    (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4) (t : ℝ) :
    ‖exp ((-(t : ℂ) * Complex.I) • G.checkerR J h s) - evo (G.checkerEff J h s) t‖ ≤
      |t| * (4 * J ^ 4 / |h| ^ 3) := by
  have hK := G.norm_cplx_Ks_le J hJ s hs hdeg
  have ha : 0 < |h| := abs_pos.mpr hh
  calc ‖exp ((-(t : ℂ) * Complex.I) • G.checkerR J h s) - evo (G.checkerEff J h s) t‖
      ≤ |t| * (‖cplx (G.Ks J s)‖ ^ 4 / (64 * |h| ^ 3)) := G.norm_exp_checkerR_sub_le J h hh s t
    _ ≤ |t| * ((4 * J) ^ 4 / (64 * |h| ^ 3)) := by gcongr
    _ = |t| * (4 * J ^ 4 / |h| ^ 3) := by field_simp; ring

/-! ### Structure of `K_s²` -/

/-- `K_s²` connects only sites of the same sublattice. -/
theorem Ks_sq_sublattice (J : ℝ) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε : ∀ e, ε (G.src e) = -ε (G.tgt e)) (i k : Fin N) (hik : ε i ≠ ε k) :
    (G.Ks J s * G.Ks J s) i k = 0 := by
  rw [Matrix.mul_apply]
  refine Finset.sum_eq_zero fun j _ => CheckerStrongAux.Ks_mul_Ks_eq_zero G J s i j k ?_
  rintro ⟨h1, h2⟩
  have e1 := CheckerStrongAux.eps_of_adj G ε hε h1
  have e2 := CheckerStrongAux.eps_of_adj G ε hε h2
  exact hik (by linarith)

/-- `K_s²` vanishes between sites without a common neighbour. -/
theorem Ks_sq_far (J : ℝ) (s : Fin M → ℝ) (i k : Fin N)
    (hfar : ∀ j, ¬ ((∃ e, (G.src e = i ∧ G.tgt e = j) ∨ (G.tgt e = i ∧ G.src e = j)) ∧
      (∃ e, (G.src e = j ∧ G.tgt e = k) ∨ (G.tgt e = j ∧ G.src e = k)))) :
    (G.Ks J s * G.Ks J s) i k = 0 := by
  rw [Matrix.mul_apply]
  exact Finset.sum_eq_zero fun j _ => CheckerStrongAux.Ks_mul_Ks_eq_zero G J s i j k (hfar j)

/-- **Plaquette diagonal**: on a simple graph, if `j₁ ≠ j₂` are the only common neighbours of `i` and
`k`, joined by the edges `e₁ = {i, j₁}`, `e₂ = {j₁, k}`, `e₃ = {k, j₂}`, `e₄ = {j₂, i}`, then
`(K_s²)_{ik} = J² s_{e₁} s_{e₂} (1 + W)` with the plaquette flux `W = s_{e₁} s_{e₂} s_{e₃} s_{e₄}`. -/
theorem Ks_sq_plaquette (J : ℝ) (s : Fin M → ℝ) (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hsimple : ∀ e e', ((G.src e = G.src e' ∧ G.tgt e = G.tgt e') ∨
      (G.src e = G.tgt e' ∧ G.tgt e = G.src e')) → e = e')
    (i k j₁ j₂ : Fin N) (hj : j₁ ≠ j₂) (e₁ e₂ e₃ e₄ : Fin M)
    (h₁ : (G.src e₁ = i ∧ G.tgt e₁ = j₁) ∨ (G.tgt e₁ = i ∧ G.src e₁ = j₁))
    (h₂ : (G.src e₂ = j₁ ∧ G.tgt e₂ = k) ∨ (G.tgt e₂ = j₁ ∧ G.src e₂ = k))
    (h₃ : (G.src e₃ = k ∧ G.tgt e₃ = j₂) ∨ (G.tgt e₃ = k ∧ G.src e₃ = j₂))
    (h₄ : (G.src e₄ = j₂ ∧ G.tgt e₄ = i) ∨ (G.tgt e₄ = j₂ ∧ G.src e₄ = i))
    (hcommon : ∀ j, (∃ e, (G.src e = i ∧ G.tgt e = j) ∨ (G.tgt e = i ∧ G.src e = j)) →
      (∃ e, (G.src e = j ∧ G.tgt e = k) ∨ (G.tgt e = j ∧ G.src e = k)) → j = j₁ ∨ j = j₂) :
    (G.Ks J s * G.Ks J s) i k = J ^ 2 * s e₁ * s e₂ * (1 + s e₁ * s e₂ * s e₃ * s e₄) := by
  have hK1 : G.Ks J s i j₁ = -J * s e₁ := CheckerStrongAux.Ks_of_edge G J s hsimple h₁
  have hK2 : G.Ks J s j₁ k = -J * s e₂ := CheckerStrongAux.Ks_of_edge G J s hsimple h₂
  have hK3 : G.Ks J s j₂ k = -J * s e₃ := by
    rw [ClusterAux.Ks_symm G J s j₂ k]
    exact CheckerStrongAux.Ks_of_edge G J s hsimple h₃
  have hK4 : G.Ks J s i j₂ = -J * s e₄ := by
    rw [ClusterAux.Ks_symm G J s i j₂]
    exact CheckerStrongAux.Ks_of_edge G J s hsimple h₄
  rw [Matrix.mul_apply, Fintype.sum_eq_add j₁ j₂ hj]
  · rw [hK1, hK2, hK3, hK4]
    rcases hs e₁ with h1 | h1 <;> rcases hs e₂ with h2 | h2 <;> rw [h1, h2] <;> ring
  · rintro j ⟨hj1, hj2⟩
    refine CheckerStrongAux.Ks_mul_Ks_eq_zero G J s i j k ?_
    rintro ⟨ha, hb⟩
    rcases hcommon j ha hb with h | h
    · exact hj1 h
    · exact hj2 h

/-- **Straight path of length two**: on a simple graph, if `j` is the only common neighbour of `i`
and `k`, joined by the edges `e₁ = {i, j}` and `e₂ = {j, k}`, then `(K_s²)_{ik} = J² s_{e₁} s_{e₂}`. -/
theorem Ks_sq_line (J : ℝ) (s : Fin M → ℝ)
    (hsimple : ∀ e e', ((G.src e = G.src e' ∧ G.tgt e = G.tgt e') ∨
      (G.src e = G.tgt e' ∧ G.tgt e = G.src e')) → e = e')
    (i k j : Fin N) (e₁ e₂ : Fin M)
    (h₁ : (G.src e₁ = i ∧ G.tgt e₁ = j) ∨ (G.tgt e₁ = i ∧ G.src e₁ = j))
    (h₂ : (G.src e₂ = j ∧ G.tgt e₂ = k) ∨ (G.tgt e₂ = j ∧ G.src e₂ = k))
    (hcommon : ∀ j', (∃ e, (G.src e = i ∧ G.tgt e = j') ∨ (G.tgt e = i ∧ G.src e = j')) →
      (∃ e, (G.src e = j' ∧ G.tgt e = k) ∨ (G.tgt e = j' ∧ G.src e = k)) → j' = j) :
    (G.Ks J s * G.Ks J s) i k = J ^ 2 * s e₁ * s e₂ := by
  rw [Matrix.mul_apply, Fintype.sum_eq_single j]
  · rw [CheckerStrongAux.Ks_of_edge G J s hsimple h₁, CheckerStrongAux.Ks_of_edge G J s hsimple h₂]
    ring
  · intro j' hj'
    refine CheckerStrongAux.Ks_mul_Ks_eq_zero G J s i j' k ?_
    rintro ⟨ha, hb⟩
    exact hj' (hcommon j' ha hb)

/-! ### The planar box -/

open TwoDLatticesAux

/-- The corner `(x + a, y + b)` of the box plaquette `p = (x, y)` (`a, b ∈ {0, 1}`). -/
def plaqCorner (L₁ L₂ : ℕ) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) (a b : Fin 2) : Fin (L₁ * L₂) :=
  boxVert ⟨p.1.val + a.val, by have := p.1.isLt; have := a.isLt; omega⟩
    ⟨p.2.val + b.val, by have := p.2.isLt; have := b.isLt; omega⟩

/-- The plaquette flux `W_p = ∏_{e ∈ ∂p} s_e` of a box plaquette. -/
def boxFlux (L₁ L₂ : ℕ) (s : Fin (boxM L₁ L₂) → ℝ) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) : ℝ :=
  s (boxBot L₁ L₂ p) * s (boxRight L₁ L₂ p) * s (boxTop L₁ L₂ p) * s (boxLeft L₁ L₂ p)

namespace CheckerStrongAux

/-- First coordinate of a box vertex. -/
def bx {L₁ L₂ : ℕ} (v : Fin (L₁ * L₂)) : ℕ := (finProdFinEquiv.symm v).1.val

/-- Second coordinate of a box vertex. -/
def bY {L₁ L₂ : ℕ} (v : Fin (L₁ * L₂)) : ℕ := (finProdFinEquiv.symm v).2.val

/-- First coordinate of `(x, y)`. -/
@[simp] lemma bx_boxVert {L₁ L₂ : ℕ} (x : Fin L₁) (y : Fin L₂) : bx (boxVert x y) = x.val := by
  simp [bx, boxVert]

/-- Second coordinate of `(x, y)`. -/
@[simp] lemma bY_boxVert {L₁ L₂ : ℕ} (x : Fin L₁) (y : Fin L₂) : bY (boxVert x y) = y.val := by
  simp [bY, boxVert]

/-- A box vertex is determined by its coordinates. -/
lemma vert_ext {L₁ L₂ : ℕ} {v w : Fin (L₁ * L₂)} (h1 : bx v = bx w) (h2 : bY v = bY w) :
    v = w :=
  finProdFinEquiv.symm.injective (Prod.ext (Fin.ext h1) (Fin.ext h2))

/-- First coordinate of a plaquette corner. -/
@[simp] lemma bx_plaqCorner (L₁ L₂ : ℕ) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) (a b : Fin 2) :
    bx (plaqCorner L₁ L₂ p a b) = p.1.val + a.val := by
  simp [plaqCorner]

/-- Second coordinate of a plaquette corner. -/
@[simp] lemma bY_plaqCorner (L₁ L₂ : ℕ) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) (a b : Fin 2) :
    bY (plaqCorner L₁ L₂ p a b) = p.2.val + b.val := by
  simp [plaqCorner]

/-- Every box edge is a unit step to the right or upwards. -/
lemma box_edge_coords (L₁ L₂ : ℕ) (e : Fin (boxM L₁ L₂)) :
    (bx ((box L₁ L₂).tgt e) = bx ((box L₁ L₂).src e) + 1 ∧
      bY ((box L₁ L₂).tgt e) = bY ((box L₁ L₂).src e)) ∨
    (bx ((box L₁ L₂).tgt e) = bx ((box L₁ L₂).src e) ∧
      bY ((box L₁ L₂).tgt e) = bY ((box L₁ L₂).src e) + 1) := by
  obtain ⟨a, rfl⟩ := (boxEdge L₁ L₂).surjective e
  rw [box_src, box_tgt]
  rcases a with ⟨x, y⟩ | ⟨x, y⟩
  · left
    simp [boxSrc, boxTgt]
  · right
    simp [boxSrc, boxTgt]

/-- A box edge is determined by its endpoints. -/
lemma box_edge_ext (L₁ L₂ : ℕ) {e e' : Fin (boxM L₁ L₂)}
    (h1 : (box L₁ L₂).src e = (box L₁ L₂).src e')
    (h2 : (box L₁ L₂).tgt e = (box L₁ L₂).tgt e') : e = e' := by
  obtain ⟨a, rfl⟩ := (boxEdge L₁ L₂).surjective e
  obtain ⟨b, rfl⟩ := (boxEdge L₁ L₂).surjective e'
  rw [box_src, box_src] at h1
  rw [box_tgt, box_tgt] at h2
  have c1 := congrArg bx h1
  have c2 := congrArg bY h1
  have c3 := congrArg bx h2
  have c4 := congrArg bY h2
  congr 1
  rcases a with ⟨x, y⟩ | ⟨x, y⟩ <;> rcases b with ⟨x', y'⟩ | ⟨x', y'⟩ <;>
    simp only [boxSrc, boxTgt, bx_boxVert, bY_boxVert] at c1 c2 c3 c4
  · simp only [Sum.inl.injEq, Prod.mk.injEq, Fin.ext_iff]
    omega
  · exfalso
    omega
  · exfalso
    omega
  · simp only [Sum.inr.injEq, Prod.mk.injEq, Fin.ext_iff]
    omega

/-- Adjacent box vertices differ by a unit step. -/
lemma box_adj_coords (L₁ L₂ : ℕ) {i j : Fin (L₁ * L₂)}
    (h : ∃ e, ((box L₁ L₂).src e = i ∧ (box L₁ L₂).tgt e = j) ∨
      ((box L₁ L₂).tgt e = i ∧ (box L₁ L₂).src e = j)) :
    (bx j = bx i + 1 ∧ bY j = bY i) ∨ (bx i = bx j + 1 ∧ bY j = bY i) ∨
      (bx j = bx i ∧ bY j = bY i + 1) ∨ (bx j = bx i ∧ bY i = bY j + 1) := by
  obtain ⟨e, ⟨h1, h2⟩ | ⟨h1, h2⟩⟩ := h
  · have := box_edge_coords L₁ L₂ e
    rw [h1, h2] at this
    omega
  · have := box_edge_coords L₁ L₂ e
    rw [h1, h2] at this
    omega

end CheckerStrongAux

/-- The box is a simple graph. -/
theorem box_simple (L₁ L₂ : ℕ) (e e' : Fin (boxM L₁ L₂))
    (h : ((box L₁ L₂).src e = (box L₁ L₂).src e' ∧ (box L₁ L₂).tgt e = (box L₁ L₂).tgt e') ∨
      ((box L₁ L₂).src e = (box L₁ L₂).tgt e' ∧ (box L₁ L₂).tgt e = (box L₁ L₂).src e')) :
    e = e' := by
  rcases h with ⟨h1, h2⟩ | ⟨h1, h2⟩
  · exact CheckerStrongAux.box_edge_ext L₁ L₂ h1 h2
  · exfalso
    have k1 := CheckerStrongAux.box_edge_coords L₁ L₂ e
    have k2 := CheckerStrongAux.box_edge_coords L₁ L₂ e'
    rw [h1, h2] at k1
    omega

/-- Every vertex of the box has degree at most four. -/
theorem box_deg_le_four (L₁ L₂ : ℕ) (i : Fin (L₁ * L₂)) :
    (univ.filter (fun e => (box L₁ L₂).src e = i ∨ (box L₁ L₂).tgt e = i)).card ≤ 4 := by
  -- an edge at `i` is determined by its direction and by whether it starts at `i`
  let f : Fin (boxM L₁ L₂) → Bool × Bool := fun e =>
    (decide (CheckerStrongAux.bx ((box L₁ L₂).tgt e) =
      CheckerStrongAux.bx ((box L₁ L₂).src e) + 1), decide ((box L₁ L₂).src e = i))
  have hinj : Set.InjOn f (univ.filter (fun e => (box L₁ L₂).src e = i ∨
      (box L₁ L₂).tgt e = i) : Finset (Fin (boxM L₁ L₂))) := by
    intro e he e' he' hf
    simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq] at he he'
    simp only [f, Prod.mk.injEq, decide_eq_decide] at hf
    obtain ⟨hd, hs⟩ := hf
    have k1 := CheckerStrongAux.box_edge_coords L₁ L₂ e
    have k2 := CheckerStrongAux.box_edge_coords L₁ L₂ e'
    by_cases hsrc : (box L₁ L₂).src e = i
    · have hsrc' : (box L₁ L₂).src e' = i := hs.mp hsrc
      have hx : CheckerStrongAux.bx ((box L₁ L₂).src e) =
          CheckerStrongAux.bx ((box L₁ L₂).src e') := by rw [hsrc, hsrc']
      have hy : CheckerStrongAux.bY ((box L₁ L₂).src e) =
          CheckerStrongAux.bY ((box L₁ L₂).src e') := by rw [hsrc, hsrc']
      refine CheckerStrongAux.box_edge_ext L₁ L₂ (hsrc.trans hsrc'.symm)
        (CheckerStrongAux.vert_ext ?_ ?_) <;> omega
    · have htgt : (box L₁ L₂).tgt e = i := he.resolve_left hsrc
      have htgt' : (box L₁ L₂).tgt e' = i := he'.resolve_left (fun h => hsrc (hs.mpr h))
      have hx : CheckerStrongAux.bx ((box L₁ L₂).tgt e) =
          CheckerStrongAux.bx ((box L₁ L₂).tgt e') := by rw [htgt, htgt']
      have hy : CheckerStrongAux.bY ((box L₁ L₂).tgt e) =
          CheckerStrongAux.bY ((box L₁ L₂).tgt e') := by rw [htgt, htgt']
      refine CheckerStrongAux.box_edge_ext L₁ L₂ (CheckerStrongAux.vert_ext ?_ ?_)
        (htgt.trans htgt'.symm) <;> omega
  calc (univ.filter (fun e => (box L₁ L₂).src e = i ∨ (box L₁ L₂).tgt e = i)).card
      ≤ (univ : Finset (Bool × Bool)).card :=
        Finset.card_le_card_of_injOn f (fun _ _ => Finset.mem_coe.mpr (Finset.mem_univ _)) hinj
    _ = 4 := rfl

/-- **Plaquette diagonals of the box** (paper: Remark `rem:checkerstrong`):
`(K_s²)_{(x,y),(x+1,y+1)} = J² s_bot s_right (1 + W_p)`. -/
theorem box_Ks_sq_diag (L₁ L₂ : ℕ) (J : ℝ) (s : Fin (boxM L₁ L₂) → ℝ)
    (hs : ∀ e, s e = 1 ∨ s e = -1) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) :
    ((box L₁ L₂).Ks J s * (box L₁ L₂).Ks J s) (plaqCorner L₁ L₂ p 0 0) (plaqCorner L₁ L₂ p 1 1) =
      J ^ 2 * s (boxBot L₁ L₂ p) * s (boxRight L₁ L₂ p) * (1 + boxFlux L₁ L₂ s p) := by
  rw [boxFlux]
  refine (box L₁ L₂).Ks_sq_plaquette J s hs (box_simple L₁ L₂) (plaqCorner L₁ L₂ p 0 0)
    (plaqCorner L₁ L₂ p 1 1) (plaqCorner L₁ L₂ p 1 0) (plaqCorner L₁ L₂ p 0 1) ?_
    (boxBot L₁ L₂ p) (boxRight L₁ L₂ p) (boxTop L₁ L₂ p) (boxLeft L₁ L₂ p)
    (Or.inl ⟨?_, ?_⟩) (Or.inl ⟨?_, ?_⟩) (Or.inr ⟨?_, ?_⟩) (Or.inr ⟨?_, ?_⟩) ?_
  · intro h
    have := congrArg CheckerStrongAux.bx h
    simp at this
  all_goals try
    exact CheckerStrongAux.vert_ext
      (by simp [boxBot, boxTop, boxLeft, boxRight, box_src, box_tgt, boxSrc, boxTgt])
      (by simp [boxBot, boxTop, boxLeft, boxRight, box_src, box_tgt, boxSrc, boxTgt])
  -- the common neighbours of `(x, y)` and `(x + 1, y + 1)`
  intro j hj1 hj2
  have a1 := CheckerStrongAux.box_adj_coords L₁ L₂ hj1
  have a2 := CheckerStrongAux.box_adj_coords L₁ L₂ hj2
  simp only [CheckerStrongAux.bx_plaqCorner, CheckerStrongAux.bY_plaqCorner, Fin.isValue,
    Fin.val_zero, Fin.val_one, add_zero] at a1 a2
  rcases (by omega : (CheckerStrongAux.bx j = p.1.val + 1 ∧ CheckerStrongAux.bY j = p.2.val) ∨
      (CheckerStrongAux.bx j = p.1.val ∧ CheckerStrongAux.bY j = p.2.val + 1)) with
    ⟨h1, h2⟩ | ⟨h1, h2⟩
  · left
    exact CheckerStrongAux.vert_ext (by simp [h1]) (by simp [h2])
  · right
    exact CheckerStrongAux.vert_ext (by simp [h1]) (by simp [h2])

/-- The other diagonal: `(K_s²)_{(x+1,y),(x,y+1)} = J² s_bot s_left (1 + W_p)`. -/
theorem box_Ks_sq_antidiag (L₁ L₂ : ℕ) (J : ℝ) (s : Fin (boxM L₁ L₂) → ℝ)
    (hs : ∀ e, s e = 1 ∨ s e = -1) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) :
    ((box L₁ L₂).Ks J s * (box L₁ L₂).Ks J s) (plaqCorner L₁ L₂ p 1 0) (plaqCorner L₁ L₂ p 0 1) =
      J ^ 2 * s (boxBot L₁ L₂ p) * s (boxLeft L₁ L₂ p) * (1 + boxFlux L₁ L₂ s p) := by
  have e : boxFlux L₁ L₂ s p =
      s (boxBot L₁ L₂ p) * s (boxLeft L₁ L₂ p) * s (boxTop L₁ L₂ p) * s (boxRight L₁ L₂ p) := by
    rw [boxFlux]
    ring
  rw [e]
  refine (box L₁ L₂).Ks_sq_plaquette J s hs (box_simple L₁ L₂) (plaqCorner L₁ L₂ p 1 0)
    (plaqCorner L₁ L₂ p 0 1) (plaqCorner L₁ L₂ p 0 0) (plaqCorner L₁ L₂ p 1 1) ?_
    (boxBot L₁ L₂ p) (boxLeft L₁ L₂ p) (boxTop L₁ L₂ p) (boxRight L₁ L₂ p)
    (Or.inr ⟨?_, ?_⟩) (Or.inl ⟨?_, ?_⟩) (Or.inl ⟨?_, ?_⟩) (Or.inr ⟨?_, ?_⟩) ?_
  · intro h
    have := congrArg CheckerStrongAux.bx h
    simp at this
  all_goals try
    exact CheckerStrongAux.vert_ext
      (by simp [boxBot, boxTop, boxLeft, boxRight, box_src, box_tgt, boxSrc, boxTgt])
      (by simp [boxBot, boxTop, boxLeft, boxRight, box_src, box_tgt, boxSrc, boxTgt])
  -- the common neighbours of `(x + 1, y)` and `(x, y + 1)`
  intro j hj1 hj2
  have a1 := CheckerStrongAux.box_adj_coords L₁ L₂ hj1
  have a2 := CheckerStrongAux.box_adj_coords L₁ L₂ hj2
  simp only [CheckerStrongAux.bx_plaqCorner, CheckerStrongAux.bY_plaqCorner, Fin.isValue,
    Fin.val_zero, Fin.val_one, add_zero] at a1 a2
  rcases (by omega : (CheckerStrongAux.bx j = p.1.val ∧ CheckerStrongAux.bY j = p.2.val) ∨
      (CheckerStrongAux.bx j = p.1.val + 1 ∧ CheckerStrongAux.bY j = p.2.val + 1)) with
    ⟨h1, h2⟩ | ⟨h1, h2⟩
  · left
    exact CheckerStrongAux.vert_ext (by simp [h1]) (by simp [h2])
  · right
    exact CheckerStrongAux.vert_ext (by simp [h1]) (by simp [h2])

/-- **Straight next-nearest-neighbour bonds of the box**, horizontal:
`(K_s²)_{(x,y),(x+2,y)} = J² s_{(x,y)—(x+1,y)} s_{(x+1,y)—(x+2,y)}`. -/
theorem box_Ks_sq_horiz (L₁ L₂ : ℕ) (J : ℝ) (s : Fin (boxM L₁ L₂) → ℝ)
    (x : Fin (L₁ - 2)) (y : Fin L₂) :
    ((box L₁ L₂).Ks J s * (box L₁ L₂).Ks J s)
        (boxVert ⟨x.val, by have := x.isLt; omega⟩ y)
        (boxVert ⟨x.val + 2, by have := x.isLt; omega⟩ y) =
      J ^ 2 * s (boxEdge L₁ L₂ (Sum.inl (⟨x.val, by have := x.isLt; omega⟩, y))) *
        s (boxEdge L₁ L₂ (Sum.inl (⟨x.val + 1, by have := x.isLt; omega⟩, y))) := by
  refine (box L₁ L₂).Ks_sq_line J s (box_simple L₁ L₂) _ _
    (boxVert ⟨x.val + 1, by have := x.isLt; omega⟩ y) _ _ (Or.inl ⟨?_, ?_⟩) (Or.inl ⟨?_, ?_⟩) ?_
  all_goals try
    exact CheckerStrongAux.vert_ext (by simp [box_src, box_tgt, boxSrc, boxTgt])
      (by simp [box_src, box_tgt, boxSrc, boxTgt])
  -- the only common neighbour of `(x, y)` and `(x + 2, y)` is `(x + 1, y)`
  intro j hj1 hj2
  have a1 := CheckerStrongAux.box_adj_coords L₁ L₂ hj1
  have a2 := CheckerStrongAux.box_adj_coords L₁ L₂ hj2
  simp only [CheckerStrongAux.bx_boxVert, CheckerStrongAux.bY_boxVert] at a1 a2
  exact CheckerStrongAux.vert_ext (by simp only [CheckerStrongAux.bx_boxVert]; omega)
    (by simp only [CheckerStrongAux.bY_boxVert]; omega)

/-- Vertical: `(K_s²)_{(x,y),(x,y+2)} = J² s_{(x,y)—(x,y+1)} s_{(x,y+1)—(x,y+2)}`. -/
theorem box_Ks_sq_vert (L₁ L₂ : ℕ) (J : ℝ) (s : Fin (boxM L₁ L₂) → ℝ)
    (x : Fin L₁) (y : Fin (L₂ - 2)) :
    ((box L₁ L₂).Ks J s * (box L₁ L₂).Ks J s)
        (boxVert x ⟨y.val, by have := y.isLt; omega⟩)
        (boxVert x ⟨y.val + 2, by have := y.isLt; omega⟩) =
      J ^ 2 * s (boxEdge L₁ L₂ (Sum.inr (x, ⟨y.val, by have := y.isLt; omega⟩))) *
        s (boxEdge L₁ L₂ (Sum.inr (x, ⟨y.val + 1, by have := y.isLt; omega⟩))) := by
  refine (box L₁ L₂).Ks_sq_line J s (box_simple L₁ L₂) _ _
    (boxVert x ⟨y.val + 1, by have := y.isLt; omega⟩) _ _ (Or.inl ⟨?_, ?_⟩) (Or.inl ⟨?_, ?_⟩) ?_
  all_goals try
    exact CheckerStrongAux.vert_ext (by simp [box_src, box_tgt, boxSrc, boxTgt])
      (by simp [box_src, box_tgt, boxSrc, boxTgt])
  -- the only common neighbour of `(x, y)` and `(x, y + 2)` is `(x, y + 1)`
  intro j hj1 hj2
  have a1 := CheckerStrongAux.box_adj_coords L₁ L₂ hj1
  have a2 := CheckerStrongAux.box_adj_coords L₁ L₂ hj2
  simp only [CheckerStrongAux.bx_boxVert, CheckerStrongAux.bY_boxVert] at a1 a2
  exact CheckerStrongAux.vert_ext (by simp only [CheckerStrongAux.bx_boxVert]; omega)
    (by simp only [CheckerStrongAux.bY_boxVert]; omega)

end FinGraph

end DFL
