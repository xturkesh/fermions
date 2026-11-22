import DFL.TwoD

/-!
# Strong coupling: cluster confinement up to times of order `|h|/J²`
(paper: Proposition `prop:cluster`)

On a finite graph with at most four edges at every vertex (e.g. a box in `ℤ²`), with link signs
`s = ±1` and charges `q = ±1`, let `h = h_{q,s} = K_s + 2h diag(q)` and let `h⁰` be `h` with the
hopping across all bonds joining opposite charges removed (`hqs0`, the pinching of `h` by the
projections `P_± = 1_{q = ±1}`). The equal-charge clusters are the connected components of the
subgraphs induced on `{q = +1}` and `{q = -1}`.

* `hqs0_cluster` (paper: (a)): `h⁰` leaves `ℓ²(C)` invariant for every equal-charge cluster `C`.
* `cluster_bound` (paper: (b)): for `|h| > 2J` and `g = 4|h| - 8J`,
  `‖e^{-ith} - e^{-ith⁰}‖ ≤ 8J/g + 48|t|J²/g` (operator norm), uniformly in the graph, `q` and `s`.
* `cluster_leak` (paper: (b)): the weight of `e^{-ith} δ_x` outside the cluster of `x` is at most
  `(8J/g + 48|t|J²/g)²`.
* `cluster_avg` (paper: (c), up to the percolation input): for any law of `q` with density at most
  `K` with respect to the uniform law, the averaged weight outside a region `B ∋ x` is at most
  `(8J/g + 48|t|J²/g)² + K · P_unif(C_q(x) ⊄ B)`. The exponential bound on `P_unif(C_q(x) ⊄ B_R(x))`
  for i.i.d. charges on `ℤ²` is the percolation estimate of Li–Liang (paper: proof of (c)); it is
  not formalized here.

## Proof of (b)

We follow the paper with two simplifications, both finite-dimensional and with better constants.
* With `Γ = ε diag(q)` (`ε = sgn h`), `A = Γ h⁰ = 2|h| + B` with `‖B‖ ≤ 4J` (Schur test), and
  `V = h - h⁰` anticommutes with `Γ`. The map `S ↦ AS + SA` is bounded below by
  `2(2|h| - 4J) = g`, hence invertible; its preimage `S` of `ΓV` anticommutes with `Γ` and solves
  `[h⁰, S] = V` with `‖S‖ ≤ ‖ΓV‖/g ≤ 4J/g` (`ClusterAux.exists_sylvester`).
* Instead of the rotation `e^S`, the mean value inequality for `Z(u) = e^{iuh} (1 - S) e^{-iuh⁰}`,
  whose derivative is `-i e^{iuh} V S e^{-iuh⁰}`, gives
  `‖e^{-ith} - e^{-ith⁰}‖ ≤ 2‖S‖ + |t| ‖V‖ ‖S‖ ≤ 8J/g + 16|t|J²/g` (`ClusterAux.norm_evo_sub_le`).
-/

open Matrix Finset
open scoped Matrix.Norms.L2Operator

noncomputable section

namespace DFL

namespace ClusterAux

open NormedSpace

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-! ### Operator-norm estimates -/

/-- **Schur test**: if `‖M i j‖ ≤ a i j` and all row and column sums of `a` are at most `R`, then
`‖M‖ ≤ R` (`ℓ²` operator norm). -/
lemma norm_le_of_entry_le (M : Matrix ι ι ℂ) (a : ι → ι → ℝ) (R : ℝ) (hR : 0 ≤ R)
    (ha : ∀ i j, ‖M i j‖ ≤ a i j) (hrow : ∀ i, ∑ j, a i j ≤ R)
    (hcol : ∀ j, ∑ i, a i j ≤ R) : ‖M‖ ≤ R := by
  rw [Matrix.cstar_norm_def]
  refine ContinuousLinearMap.opNorm_le_bound _ hR fun x => ?_
  have ha0 : ∀ i j, 0 ≤ a i j := fun i j => (norm_nonneg _).trans (ha i j)
  have key : ‖Matrix.toEuclideanCLM (n := ι) (𝕜 := ℂ) M x‖ ^ 2 ≤ (R * ‖x‖) ^ 2 := by
    rw [EuclideanSpace.norm_sq_eq, mul_pow, EuclideanSpace.norm_sq_eq]
    simp only [Matrix.ofLp_toEuclideanCLM, Matrix.mulVec, dotProduct]
    calc ∑ i, ‖∑ j, M i j * x j‖ ^ 2 ≤ ∑ i, (∑ j, a i j * ‖x j‖) ^ 2 := by
          gcongr with i
          refine (norm_sum_le _ _).trans (Finset.sum_le_sum fun j _ => ?_)
          rw [norm_mul]
          exact mul_le_mul_of_nonneg_right (ha i j) (norm_nonneg _)
      _ ≤ ∑ i, (∑ j, a i j) * (∑ j, a i j * ‖x j‖ ^ 2) := by
          gcongr with i
          exact Finset.sum_sq_le_sum_mul_sum_of_sq_le_mul _ (fun j _ => ha0 i j)
            (fun j _ => by have := ha0 i j; positivity) (fun j _ => by ring_nf; rfl)
      _ ≤ ∑ i, R * (∑ j, a i j * ‖x j‖ ^ 2) := by
          gcongr with i
          · exact Finset.sum_nonneg fun j _ => by have := ha0 i j; positivity
          · exact hrow i
      _ = R * ∑ j, (∑ i, a i j) * ‖x j‖ ^ 2 := by
          rw [← Finset.mul_sum, Finset.sum_comm]
          simp only [Finset.sum_mul]
      _ ≤ R * ∑ j, R * ‖x j‖ ^ 2 := by
          gcongr with j
          exact hcol j
      _ = R ^ 2 * ∑ j, ‖x j‖ ^ 2 := by rw [← Finset.mul_sum]; ring
  exact (sq_le_sq₀ (norm_nonneg _) (by positivity)).mp key

/-- The real-linear map `M ↦ M v`. -/
def mulVecLin (v : EuclideanSpace ℂ ι) : Matrix ι ι ℂ →ₗ[ℝ] EuclideanSpace ℂ ι where
  toFun M := WithLp.toLp 2 (M *ᵥ WithLp.ofLp v)
  map_add' M M' := by
    rw [Matrix.add_mulVec, WithLp.toLp_add]
  map_smul' c M := by
    rw [RingHom.id_apply, ← WithLp.toLp_smul, Matrix.smul_mulVec]

/-- The continuous real-linear map `M ↦ M v`. -/
def mulVecL (v : EuclideanSpace ℂ ι) : Matrix ι ι ℂ →L[ℝ] EuclideanSpace ℂ ι :=
  LinearMap.toContinuousLinearMap (mulVecLin v)

omit [DecidableEq ι] in
/-- Evaluation of `mulVecL`. -/
lemma mulVecL_apply (v : EuclideanSpace ℂ ι) (M : Matrix ι ι ℂ) :
    mulVecL v M = WithLp.toLp 2 (M *ᵥ WithLp.ofLp v) := rfl

/-- `‖M v‖ ≤ ‖M‖ ‖v‖`. -/
lemma norm_mulVecL_le (v : EuclideanSpace ℂ ι) (M : Matrix ι ι ℂ) :
    ‖mulVecL v M‖ ≤ ‖M‖ * ‖v‖ := by
  rw [Matrix.cstar_norm_def]
  exact (Matrix.toEuclideanCLM (n := ι) (𝕜 := ℂ) M).le_opNorm v

/-- An operator-norm bound from bounds on vectors. -/
lemma norm_le_of_mulVecL (M : Matrix ι ι ℂ) (C : ℝ) (hC : 0 ≤ C)
    (h : ∀ v, ‖mulVecL v M‖ ≤ C * ‖v‖) : ‖M‖ ≤ C := by
  rw [Matrix.cstar_norm_def]
  exact ContinuousLinearMap.opNorm_le_bound _ hC h

/-- The squared moduli of a column add up to at most `‖M‖²`. -/
lemma col_sum_le (M : Matrix ι ι ℂ) (x : ι) : ∑ y, ‖M y x‖ ^ 2 ≤ ‖M‖ ^ 2 := by
  have h1 := norm_mulVecL_le (EuclideanSpace.single x (1 : ℂ)) M
  rw [PiLp.norm_single, norm_one, mul_one] at h1
  have h2 : ‖mulVecL (EuclideanSpace.single x (1 : ℂ)) M‖ ^ 2 = ∑ y, ‖M y x‖ ^ 2 := by
    rw [EuclideanSpace.norm_sq_eq]
    refine Finset.sum_congr rfl fun y _ => ?_
    simp [mulVecL_apply]
  rw [← h2]
  exact pow_le_pow_left₀ (norm_nonneg _) h1 2

/-- The identity matrix has norm at most one. -/
lemma norm_one_le : ‖(1 : Matrix ι ι ℂ)‖ ≤ 1 := by
  rw [← Matrix.diagonal_one, Matrix.l2_opNorm_diagonal]
  exact (pi_norm_le_iff_of_nonneg zero_le_one).mpr fun i => by simp

/-- A unitary matrix has norm at most one. -/
lemma norm_le_one_of_unitary (U : Matrix ι ι ℂ) (hU : Uᴴ * U = 1) : ‖U‖ ≤ 1 := by
  have h1 : ‖Uᴴ * U‖ = ‖U‖ * ‖U‖ := Matrix.l2_opNorm_conjTranspose_mul_self U
  rw [hU] at h1
  nlinarith [norm_nonneg U, norm_one_le (ι := ι)]

/-- `‖E X F‖ ≤ ‖X‖` for unitary `E`, `F`. -/
lemma norm_unitary_mul_mul_le (E X F : Matrix ι ι ℂ) (hE : Eᴴ * E = 1) (hF : Fᴴ * F = 1) :
    ‖E * X * F‖ ≤ ‖X‖ := by
  calc ‖E * X * F‖ ≤ ‖E‖ * ‖X‖ * ‖F‖ := by
        refine (norm_mul_le _ _).trans ?_
        exact mul_le_mul_of_nonneg_right (norm_mul_le _ _) (norm_nonneg _)
    _ ≤ 1 * ‖X‖ * 1 := by
        gcongr
        · exact norm_le_one_of_unitary E hE
        · exact norm_le_one_of_unitary F hF
    _ = ‖X‖ := by ring

/-! ### Complexification -/

omit [Fintype ι] [DecidableEq ι] in
/-- Complexification is additive. -/
lemma cplx_sub (A B : Matrix ι ι ℝ) : cplx (A - B) = cplx A - cplx B := by
  ext i j
  simp [cplx]

omit [Fintype ι] [DecidableEq ι] in
/-- Complexification commutes with negation. -/
lemma cplx_neg (A : Matrix ι ι ℝ) : cplx (-A) = -cplx A := by
  ext i j
  simp [cplx]

omit [DecidableEq ι] in
/-- Complexification is multiplicative. -/
lemma cplx_mul (A B : Matrix ι ι ℝ) : cplx (A * B) = cplx A * cplx B :=
  Matrix.map_mul (f := Complex.ofRealHom)

omit [Fintype ι] in
/-- Complexification of a multiple of the identity. -/
lemma cplx_smul_one (c : ℝ) : cplx (c • (1 : Matrix ι ι ℝ)) = (c : ℂ) • (1 : Matrix ι ι ℂ) := by
  ext i j
  by_cases hij : i = j <;> simp [cplx, hij]

omit [Fintype ι] in
/-- Complexification of the identity. -/
lemma cplx_one : cplx (1 : Matrix ι ι ℝ) = 1 := by
  simpa using cplx_smul_one (ι := ι) 1

/-! ### Propagators -/

/-- The generator of `evo` as a real multiple of a fixed matrix. -/
lemma evo_eq_exp_smul (A : Matrix ι ι ℝ) (t : ℝ) :
    evo A t = exp (t • ((-Complex.I) • cplx A)) := by
  unfold evo
  congr 1
  ext i j
  simp only [Matrix.smul_apply, smul_eq_mul, Complex.real_smul]
  ring

/-- `evo A 0 = 1`. -/
lemma evo_zero (A : Matrix ι ι ℝ) : evo A 0 = 1 := by
  simp [evo]

/-- `e^{-itA} e^{itA} = 1`. -/
lemma evo_mul_evo_neg (A : Matrix ι ι ℝ) (t : ℝ) : evo A t * evo A (-t) = 1 := by
  unfold evo
  rw [← Matrix.exp_add_of_commute]
  · rw [← add_smul]
    simp
  · exact ((Commute.refl (cplx A)).smul_left _).smul_right _

section Deriv

open scoped Matrix.Norms.Operator

set_option backward.isDefEq.respectTransparency false

/-- Derivative of the propagator (generator on the right). -/
lemma hasDerivAt_evo (A : Matrix ι ι ℝ) (t : ℝ) :
    HasDerivAt (fun u : ℝ => evo A u) (evo A t * ((-Complex.I) • cplx A)) t := by
  simp only [evo_eq_exp_smul]
  exact hasDerivAt_exp_smul_const _ t

/-- Derivative of the propagator (generator on the left). -/
lemma hasDerivAt_evo' (A : Matrix ι ι ℝ) (t : ℝ) :
    HasDerivAt (fun u : ℝ => evo A u) (((-Complex.I) • cplx A) * evo A t) t := by
  simp only [evo_eq_exp_smul]
  exact hasDerivAt_exp_smul_const' _ t

/-- Derivative of `u ↦ e^{iuA} R e^{-iuA₀}`. -/
lemma hasDerivAt_conj (A A0 : Matrix ι ι ℝ) (R : Matrix ι ι ℂ) (u : ℝ) :
    HasDerivAt (fun u : ℝ => evo A (-u) * R * evo A0 u)
      (evo A (-u) * (Complex.I • (cplx A * R - R * cplx A0)) * evo A0 u) u := by
  have h1 := (hasDerivAt_evo A (-u)).scomp u (hasDerivAt_neg u)
  have h2 := hasDerivAt_evo' A0 u
  have h3 := (h1.mul_const R).mul h2
  convert h3 using 1
  · rfl
  · rfl
  · simp only [Function.comp_apply, neg_smul, one_smul, smul_neg, Matrix.mul_neg, Matrix.neg_mul,
      Matrix.mul_smul, Matrix.smul_mul, smul_sub, Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_assoc,
      neg_neg]
    abel

/-- Derivative of `u ↦ e^{iuA} R e^{-iuA₀} v`. -/
lemma hasDerivAt_conj_vec (A A0 : Matrix ι ι ℝ) (R : Matrix ι ι ℂ) (v : EuclideanSpace ℂ ι)
    (u : ℝ) :
    HasDerivAt (fun u : ℝ => mulVecL v (evo A (-u) * R * evo A0 u))
      (mulVecL v (evo A (-u) * (Complex.I • (cplx A * R - R * cplx A0)) * evo A0 u)) u :=
  (mulVecL v).hasFDerivAt.comp_hasDerivAt u (hasDerivAt_conj A A0 R u)

end Deriv

/-- **Propagator comparison**: if `[A₀, S] = A - A₀`, then
`‖e^{-itA} - e^{-itA₀}‖ ≤ 2‖S‖ + |t| ‖A - A₀‖ ‖S‖` (mean value inequality for
`Z(u) = e^{iuA} (1 - S) e^{-iuA₀}`, whose derivative is `-i e^{iuA} (A - A₀) S e^{-iuA₀}`). -/
lemma norm_evo_sub_le (A A0 : Matrix ι ι ℝ) (hA : A.IsSymm) (hA0 : A0.IsSymm)
    (S : Matrix ι ι ℂ) (hS : cplx A0 * S - S * cplx A0 = cplx A - cplx A0) (t : ℝ) :
    ‖evo A t - evo A0 t‖ ≤ 2 * ‖S‖ + |t| * (‖cplx A - cplx A0‖ * ‖S‖) := by
  set V := cplx A - cplx A0 with hV
  have hC : 0 ≤ ‖V‖ * ‖S‖ := by positivity
  -- the generator of `u ↦ e^{iuA} (1 - S) e^{-iuA₀}`
  have hgen : cplx A * (1 - S) - (1 - S) * cplx A0 = -(V * S) := by
    have h1 : cplx A * S - S * cplx A0 = V + V * S := by
      calc cplx A * S - S * cplx A0 = (cplx A0 * S - S * cplx A0) + (cplx A - cplx A0) * S := by
            rw [Matrix.sub_mul]; abel
        _ = V + V * S := by rw [hS, ← hV]
    calc cplx A * (1 - S) - (1 - S) * cplx A0
        = (cplx A - cplx A0) - (cplx A * S - S * cplx A0) := by
          simp only [Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_one, Matrix.one_mul]
          abel
      _ = -(V * S) := by rw [h1, ← hV]; abel
  -- mean value inequality for `Z(u) = e^{iuA} (1 - S) e^{-iuA₀}`
  have hZ : ‖evo A (-t) * (1 - S) * evo A0 t - (1 - S)‖ ≤ ‖V‖ * ‖S‖ * |t| := by
    refine norm_le_of_mulVecL _ _ (by positivity) fun v => ?_
    have key := Convex.norm_image_sub_le_of_norm_hasDerivWithin_le
      (f := fun u : ℝ => mulVecL v (evo A (-u) * (1 - S) * evo A0 u)) (s := Set.univ)
      (x := 0) (y := t) (C := ‖V‖ * ‖S‖ * ‖v‖)
      (fun u _ => (hasDerivAt_conj_vec A A0 (1 - S) v u).hasDerivWithinAt)
      (fun u _ => by
        refine (norm_mulVecL_le v _).trans (mul_le_mul_of_nonneg_right ?_ (norm_nonneg v))
        rw [hgen]
        refine (norm_unitary_mul_mul_le _ _ _ (evo_unitary A hA _) (evo_unitary A0 hA0 _)).trans ?_
        rw [smul_neg, norm_neg, norm_smul, Complex.norm_I, one_mul]
        exact norm_mul_le _ _)
      convex_univ trivial trivial
    simp only [neg_zero, evo_zero, Matrix.one_mul, Matrix.mul_one, sub_zero] at key
    rw [map_sub]
    calc ‖mulVecL v (evo A (-t) * (1 - S) * evo A0 t) - mulVecL v (1 - S)‖
        ≤ ‖V‖ * ‖S‖ * ‖v‖ * ‖t‖ := key
      _ = ‖V‖ * ‖S‖ * |t| * ‖v‖ := by rw [Real.norm_eq_abs]; ring
  -- assemble
  have hW : evo A t - evo A0 t = -(evo A t * ((evo A (-t) * (1 - S) * evo A0 t - (1 - S)) - S +
      evo A (-t) * S * evo A0 t)) := by
    have h1 := evo_mul_evo_neg A t
    simp only [Matrix.mul_sub, Matrix.sub_mul, Matrix.mul_add, Matrix.mul_one, Matrix.one_mul,
      ← Matrix.mul_assoc, h1]
    abel
  rw [hW, norm_neg]
  calc ‖evo A t * ((evo A (-t) * (1 - S) * evo A0 t - (1 - S)) - S + evo A (-t) * S * evo A0 t)‖
      ≤ ‖evo A t‖ * ‖(evo A (-t) * (1 - S) * evo A0 t - (1 - S)) - S +
          evo A (-t) * S * evo A0 t‖ := norm_mul_le _ _
    _ ≤ 1 * ‖(evo A (-t) * (1 - S) * evo A0 t - (1 - S)) - S + evo A (-t) * S * evo A0 t‖ :=
        mul_le_mul_of_nonneg_right (norm_le_one_of_unitary _ (evo_unitary A hA t)) (norm_nonneg _)
    _ ≤ ‖V‖ * ‖S‖ * |t| + ‖S‖ + ‖S‖ := by
        rw [one_mul]
        refine (norm_add_le _ _).trans (add_le_add ((norm_sub_le _ _).trans
          (add_le_add hZ le_rfl)) ?_)
        exact norm_unitary_mul_mul_le _ _ _ (evo_unitary A hA _) (evo_unitary A0 hA0 _)
    _ = 2 * ‖S‖ + |t| * (‖V‖ * ‖S‖) := by ring

/-- **Sylvester equation** (with a reflection `Γ`): if `Γ² = 1`, `Γ` commutes with `H₀` and
anticommutes with `V`, and `‖ΓH₀ - c‖ ≤ b < c`, then some `S` solves `[H₀, S] = V` with
`2(c - b)‖S‖ ≤ ‖ΓV‖`. -/
lemma exists_sylvester (Γ H0 V : Matrix ι ι ℂ) (c b : ℝ)
    (hB : ‖Γ * H0 - (c : ℂ) • (1 : Matrix ι ι ℂ)‖ ≤ b) (hbc : b < c)
    (hΓ : Γ * Γ = 1) (hΓH : Γ * H0 = H0 * Γ) (hΓV : Γ * V = -(V * Γ)) :
    ∃ S : Matrix ι ι ℂ, H0 * S - S * H0 = V ∧ 2 * (c - b) * ‖S‖ ≤ ‖Γ * V‖ := by
  obtain ⟨A, hA⟩ : ∃ A, A = Γ * H0 := ⟨_, rfl⟩
  obtain ⟨B, hBdef⟩ : ∃ B, B = A - (c : ℂ) • (1 : Matrix ι ι ℂ) := ⟨_, rfl⟩
  rw [← hA, ← hBdef] at hB
  have hb0 : 0 ≤ b := (norm_nonneg _).trans hB
  have hc0 : 0 < c := lt_of_le_of_lt hb0 hbc
  -- lower bound for `S ↦ AS + SA`
  have hlow : ∀ S : Matrix ι ι ℂ, 2 * (c - b) * ‖S‖ ≤ ‖A * S + S * A‖ := by
    intro S
    have e : A * S + S * A = ((2 * c : ℝ) : ℂ) • S + (B * S + S * B) := by
      have hAB : A = B + (c : ℂ) • (1 : Matrix ι ι ℂ) := by rw [hBdef]; abel
      rw [hAB]
      simp only [Matrix.add_mul, Matrix.mul_add, Matrix.smul_mul, Matrix.mul_smul,
        Matrix.one_mul, Matrix.mul_one]
      rw [show ((2 * c : ℝ) : ℂ) = (c : ℂ) + (c : ℂ) by push_cast; ring, add_smul]
      abel
    have h1 : ‖B * S + S * B‖ ≤ 2 * b * ‖S‖ := by
      calc ‖B * S + S * B‖ ≤ ‖B * S‖ + ‖S * B‖ := norm_add_le _ _
        _ ≤ ‖B‖ * ‖S‖ + ‖S‖ * ‖B‖ := add_le_add (norm_mul_le _ _) (norm_mul_le _ _)
        _ ≤ b * ‖S‖ + ‖S‖ * b := by gcongr
        _ = 2 * b * ‖S‖ := by ring
    have h2 : ‖((2 * c : ℝ) : ℂ) • S‖ = 2 * c * ‖S‖ := by
      rw [norm_smul, Complex.norm_real, Real.norm_eq_abs, abs_of_pos (by positivity)]
    have h3 : ‖((2 * c : ℝ) : ℂ) • S‖ ≤ ‖((2 * c : ℝ) : ℂ) • S + (B * S + S * B)‖ +
        ‖B * S + S * B‖ := by
      simpa using norm_sub_le (((2 * c : ℝ) : ℂ) • S + (B * S + S * B)) (B * S + S * B)
    rw [e]
    linarith
  -- the linear map `S ↦ AS + SA` is injective, hence surjective
  let L : Matrix ι ι ℂ →ₗ[ℂ] Matrix ι ι ℂ := LinearMap.mulLeft ℂ A + LinearMap.mulRight ℂ A
  have hL : ∀ S, L S = A * S + S * A := fun S => rfl
  have hinj : Function.Injective L := by
    rw [← LinearMap.ker_eq_bot, LinearMap.ker_eq_bot']
    intro S hS
    have h := hlow S
    rw [← hL, hS, norm_zero] at h
    have h' : ‖S‖ ≤ 0 := by nlinarith [norm_nonneg S]
    exact norm_le_zero_iff.mp h'
  obtain ⟨S, hS⟩ := (LinearMap.injective_iff_surjective.mp hinj) (Γ * V)
  rw [hL] at hS
  have hAΓ : A * Γ = H0 := by
    rw [hA, Matrix.mul_assoc, ← hΓH, ← Matrix.mul_assoc, hΓ, Matrix.one_mul]
  have hΓA : Γ * A = H0 := by
    rw [hA, ← Matrix.mul_assoc, hΓ, Matrix.one_mul]
  -- `S` anticommutes with `Γ`, by uniqueness
  have hanti : -(Γ * S * Γ) = S := by
    apply hinj
    rw [hL, hL, hS]
    have e1 : A * (Γ * S * Γ) + Γ * S * Γ * A = Γ * (A * S + S * A) * Γ := by
      calc A * (Γ * S * Γ) + Γ * S * Γ * A = (A * Γ) * S * Γ + Γ * S * (Γ * A) := by
            simp only [Matrix.mul_assoc]
        _ = (Γ * A) * S * Γ + Γ * S * (A * Γ) := by rw [hAΓ, hΓA]
        _ = Γ * (A * S + S * A) * Γ := by
            simp only [Matrix.mul_add, Matrix.add_mul, Matrix.mul_assoc]
    calc A * -(Γ * S * Γ) + -(Γ * S * Γ) * A = -(A * (Γ * S * Γ) + Γ * S * Γ * A) := by
          rw [Matrix.mul_neg, Matrix.neg_mul, neg_add]
      _ = -(Γ * (Γ * V) * Γ) := by rw [e1, hS]
      _ = -(V * Γ) := by rw [← Matrix.mul_assoc Γ Γ V, hΓ, Matrix.one_mul]
      _ = Γ * V := hΓV.symm
  have hSΓ : S * Γ = -(Γ * S) := by
    calc S * Γ = -(Γ * S * Γ) * Γ := by rw [hanti]
      _ = -(Γ * S * (Γ * Γ)) := by rw [Matrix.neg_mul, Matrix.mul_assoc (Γ * S)]
      _ = -(Γ * S) := by rw [hΓ, Matrix.mul_one]
  refine ⟨S, ?_, ?_⟩
  · calc H0 * S - S * H0 = (Γ * A) * S - S * (Γ * A) := by rw [hΓA]
      _ = Γ * (A * S) - (S * Γ) * A := by simp only [Matrix.mul_assoc]
      _ = Γ * (A * S) + Γ * (S * A) := by
          rw [hSΓ, Matrix.neg_mul, sub_neg_eq_add, Matrix.mul_assoc]
      _ = Γ * (A * S + S * A) := by rw [Matrix.mul_add]
      _ = Γ * (Γ * V) := by rw [hS]
      _ = V := by rw [← Matrix.mul_assoc, hΓ, Matrix.one_mul]
  · calc 2 * (c - b) * ‖S‖ ≤ ‖A * S + S * A‖ := hlow S
      _ = ‖Γ * V‖ := by rw [hS]

/-! ### The hopping operator `K_s` -/

variable {N M : ℕ} (G : FinGraph N M)

/-- The edges joining `i` to `j` are the edges at `i` whose other endpoint is `j`. -/
lemma filter_join_eq (i j : Fin N) :
    (univ.filter (fun e => (G.src e = i ∧ G.tgt e = j) ∨ (G.tgt e = i ∧ G.src e = j))) =
      (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).filter
        (fun e => (if G.src e = i then G.tgt e else G.src e) = j) := by
  ext e
  simp only [mem_filter, mem_univ, true_and]
  have hne := G.src_ne_tgt e
  by_cases h1 : G.src e = i
  · simp only [h1, true_or, ite_true, true_and]
    constructor
    · rintro (h | ⟨h2, _⟩)
      · exact h
      · exact absurd (h1.trans h2.symm) hne
    · intro h; exact Or.inl h
  · simp only [h1, false_and, false_or, ite_false]

/-- Row sums of `|K_s|` are at most `4J` (at most four edges at every vertex, `|s_e| = 1`). -/
lemma Ks_row_sum (J : ℝ) (hJ : 0 ≤ J) (s : Fin M → ℝ) (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4) (i : Fin N) :
    ∑ j, |G.Ks J s i j| ≤ 4 * J := by
  have h1 : ∀ j, |G.Ks J s i j| ≤ J * ((univ.filter (fun e => (G.src e = i ∧ G.tgt e = j) ∨
      (G.tgt e = i ∧ G.src e = j))).card : ℝ) := by
    intro j
    simp only [FinGraph.Ks, abs_mul, abs_neg, abs_of_nonneg hJ]
    refine mul_le_mul_of_nonneg_left ?_ hJ
    refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
    rw [Finset.card_eq_sum_ones, Nat.cast_sum]
    refine Finset.sum_le_sum fun e _ => ?_
    rcases hs e with h | h <;> simp [h]
  have h2 : ∑ j, ((univ.filter (fun e => (G.src e = i ∧ G.tgt e = j) ∨
      (G.tgt e = i ∧ G.src e = j))).card : ℝ) =
      ((univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card : ℝ) := by
    rw [Finset.card_eq_sum_card_fiberwise (f := fun e => if G.src e = i then G.tgt e else G.src e)
      (t := univ) (fun _ _ => Finset.mem_coe.mpr (mem_univ _))]
    push_cast
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [filter_join_eq]
  calc ∑ j, |G.Ks J s i j| ≤ ∑ j, J * ((univ.filter (fun e => (G.src e = i ∧ G.tgt e = j) ∨
      (G.tgt e = i ∧ G.src e = j))).card : ℝ) := Finset.sum_le_sum fun j _ => h1 j
    _ = J * ((univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card : ℝ) := by
      rw [← Finset.mul_sum, h2]
    _ ≤ J * 4 := by
      refine mul_le_mul_of_nonneg_left ?_ hJ
      exact_mod_cast hdeg i
    _ = 4 * J := by ring

/-- `K_s` is symmetric. -/
lemma Ks_symm (J : ℝ) (s : Fin M → ℝ) (i j : Fin N) : G.Ks J s i j = G.Ks J s j i := by
  simp only [FinGraph.Ks]
  congr 2
  ext e
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]
  tauto

/-- `K_s` vanishes on the diagonal. -/
lemma Ks_diag (J : ℝ) (s : Fin M → ℝ) (i : Fin N) : G.Ks J s i i = 0 := by
  simp only [FinGraph.Ks]
  rw [Finset.filter_false_of_mem, Finset.sum_empty, mul_zero]
  intro e _
  have := G.src_ne_tgt e
  rintro (⟨h1, h2⟩ | ⟨h1, h2⟩)
  · exact this (h1.trans h2.symm)
  · exact this (h2.trans h1.symm)

/-- A real matrix dominated entrywise by `|K_s|` has complexified norm at most `4J`. -/
lemma norm_cplx_le (J : ℝ) (hJ : 0 ≤ J) (s : Fin M → ℝ) (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4)
    (A : Matrix (Fin N) (Fin N) ℝ) (hA : ∀ i j, |A i j| ≤ |G.Ks J s i j|) :
    ‖cplx A‖ ≤ 4 * J := by
  refine norm_le_of_entry_le (cplx A) (fun i j => |G.Ks J s i j|) (4 * J) (by linarith)
    (fun i j => ?_) (fun i => Ks_row_sum G J hJ s hs hdeg i) (fun j => ?_)
  · simp only [cplx, Matrix.map_apply, Complex.norm_real, Real.norm_eq_abs]
    exact hA i j
  · simp_rw [Ks_symm G J s _ j]
    exact Ks_row_sum G J hJ s hs hdeg j

end ClusterAux

namespace FinGraph

variable {N M : ℕ} (G : FinGraph N M)

/-- `h⁰_{q,s}`: the hopping across bonds joining opposite charges removed. -/
def hqs0 (J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) : Matrix (Fin N) (Fin N) ℝ :=
  fun i j => if q i = q j then G.hqs J h q s i j else 0

/-- Adjacent vertices of equal charge. -/
def eqChargeAdj (q : Fin N → ℝ) (i j : Fin N) : Prop :=
  (∃ e, (G.src e = i ∧ G.tgt e = j) ∨ (G.tgt e = i ∧ G.src e = j)) ∧ q i = q j

/-- `y` lies in the equal-charge cluster `C_q(x)` of `x`. -/
def InCluster (q : Fin N → ℝ) (x y : Fin N) : Prop :=
  Relation.ReflTransGen (G.eqChargeAdj q) x y

/-- Entries of `h_{q,s}`. -/
private lemma hqs_apply' (J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) (i j : Fin N) :
    G.hqs J h q s i j = G.Ks J s i j + 2 * h * (if i = j then q i else 0) := by
  simp only [hqs, Matrix.add_apply, Matrix.smul_apply, smul_eq_mul, Matrix.diagonal_apply]

/-- `h_{q,s}` is symmetric. -/
private lemma hqs_isSymm' (J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) :
    (G.hqs J h q s).IsSymm := by
  refine Matrix.IsSymm.ext fun i j => ?_
  rw [hqs_apply', hqs_apply', ClusterAux.Ks_symm G J s j i]
  by_cases hij : i = j
  · subst hij; rfl
  · rw [ite_eq_right hij, ite_eq_right (Ne.symm hij)]

/-- `h⁰_{q,s}` is symmetric (entrywise). -/
private lemma hqs0_symm (J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) (i j : Fin N) :
    G.hqs0 J h q s j i = G.hqs0 J h q s i j := by
  simp only [hqs0]
  by_cases hq : q i = q j
  · rw [ite_eq_left hq.symm, ite_eq_left hq, (hqs_isSymm' G J h q s).apply i j]
  · rw [ite_eq_right (Ne.symm hq), ite_eq_right hq]

/-- `h⁰_{q,s}` is symmetric. -/
private lemma hqs0_isSymm' (J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) :
    (G.hqs0 J h q s).IsSymm :=
  Matrix.IsSymm.ext fun i j => hqs0_symm G J h q s i j

/-- **Clusters are invariant** (paper: Proposition `prop:cluster`(a)): if `z ∈ C_q(x)` and
`h⁰_{yz} ≠ 0`, then `y ∈ C_q(x)`. -/
theorem hqs0_cluster (J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) (x y z : Fin N)
    (hz : G.InCluster q x z) (hyz : G.hqs0 J h q s y z ≠ 0) : G.InCluster q x y := by
  by_cases hyz' : y = z
  · rw [hyz']; exact hz
  · have hq : q y = q z := by
      by_contra hne
      exact hyz (by simp [hqs0, hne])
    have hK : G.hqs J h q s y z ≠ 0 := by
      intro h0; exact hyz (by simp [hqs0, h0])
    have hne : (univ.filter (fun e => (G.src e = y ∧ G.tgt e = z) ∨
        (G.tgt e = y ∧ G.src e = z))).Nonempty := by
      by_contra hem
      rw [Finset.not_nonempty_iff_eq_empty] at hem
      apply hK
      rw [hqs_apply', ite_eq_right hyz', mul_zero, add_zero, Ks, hem, Finset.sum_empty, mul_zero]
    obtain ⟨e, he⟩ := hne
    rw [Finset.mem_filter] at he
    refine Relation.ReflTransGen.tail hz ⟨⟨e, ?_⟩, hq.symm⟩
    rcases he.2 with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · exact Or.inr ⟨h2, h1⟩
    · exact Or.inl ⟨h2, h1⟩

/-- The reflection `Γ = ε diag(q)` squares to one. -/
private lemma refl_sq (ε : ℝ) (hε : ε = 1 ∨ ε = -1) (q : Fin N → ℝ)
    (hq : ∀ i, q i = 1 ∨ q i = -1) :
    diagonal (fun i => ε * q i) * diagonal (fun i => ε * q i) =
      (1 : Matrix (Fin N) (Fin N) ℝ) := by
  rw [diagonal_mul_diagonal, ← diagonal_one]
  congr 1
  funext i
  rcases hε with h1 | h1 <;> rcases hq i with h2 | h2 <;> rw [h1, h2] <;> norm_num

/-- `Γ` commutes with `h⁰`. -/
private lemma refl_comm (ε J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) :
    diagonal (fun i => ε * q i) * G.hqs0 J h q s =
      G.hqs0 J h q s * diagonal (fun i => ε * q i) := by
  ext i j
  rw [diagonal_mul, mul_diagonal]
  simp only [hqs0]
  split_ifs with hij
  · rw [hij]; ring
  · ring

/-- `Γ` anticommutes with `h¹ = h - h⁰`. -/
private lemma refl_anti (ε J h : ℝ) (q : Fin N → ℝ) (hq : ∀ i, q i = 1 ∨ q i = -1)
    (s : Fin M → ℝ) :
    diagonal (fun i => ε * q i) * (G.hqs J h q s - G.hqs0 J h q s) =
      -((G.hqs J h q s - G.hqs0 J h q s) * diagonal (fun i => ε * q i)) := by
  ext i j
  rw [Matrix.neg_apply, diagonal_mul, mul_diagonal]
  simp only [Matrix.sub_apply, hqs0]
  split_ifs with hij
  · ring
  · have hqj : q j = -q i := by
      rcases hq i with h1 | h1 <;> rcases hq j with h2 | h2
      · exact absurd (h1.trans h2.symm) hij
      · rw [h1, h2]
      · rw [h1, h2]; norm_num
      · exact absurd (h1.trans h2.symm) hij
    rw [hqj]; ring

/-- `|ε q_i| = 1`. -/
private lemma abs_refl (ε : ℝ) (hε : ε = 1 ∨ ε = -1) (q : Fin N → ℝ)
    (hq : ∀ i, q i = 1 ∨ q i = -1) (i : Fin N) : |ε * q i| = 1 := by
  rcases hε with h1 | h1 <;> rcases hq i with h2 | h2 <;> rw [h1, h2] <;> norm_num

/-- The entries of `Γh⁰ - 2|h|` are dominated by those of `K_s`. -/
private lemma entry_refl_hqs0 (ε J h : ℝ) (hε : ε = 1 ∨ ε = -1) (hεh : ε * h = |h|)
    (q : Fin N → ℝ) (hq : ∀ i, q i = 1 ∨ q i = -1) (s : Fin M → ℝ) (i j : Fin N) :
    |(diagonal (fun i => ε * q i) * G.hqs0 J h q s -
      (2 * |h|) • (1 : Matrix (Fin N) (Fin N) ℝ)) i j| ≤ |G.Ks J s i j| := by
  rw [Matrix.sub_apply, diagonal_mul, Matrix.smul_apply, Matrix.one_apply, smul_eq_mul]
  simp only [hqs0, hqs_apply']
  by_cases hij : i = j
  · subst hij
    simp only [ite_true, ClusterAux.Ks_diag, zero_add, mul_one]
    have h0 : ε * q i * (2 * h * q i) - 2 * |h| = 0 := by
      rcases hq i with h1 | h1 <;> rw [h1] <;> linear_combination (2 : ℝ) * hεh
    simp [h0]
  · simp only [ite_eq_right hij, mul_zero, add_zero, sub_zero]
    split_ifs
    · rw [abs_mul, abs_refl ε hε q hq i, one_mul]
    · simp

/-- The entries of `h¹ = h - h⁰` are dominated by those of `K_s`. -/
private lemma entry_V (J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) (i j : Fin N) :
    |(G.hqs J h q s - G.hqs0 J h q s) i j| ≤ |G.Ks J s i j| := by
  rw [Matrix.sub_apply]
  simp only [hqs0]
  split_ifs with hq'
  · simp
  · have hij : i ≠ j := by rintro rfl; exact hq' rfl
    rw [hqs_apply', ite_eq_right hij, mul_zero, add_zero, sub_zero]

/-- The entries of `Γh¹` are dominated by those of `K_s`. -/
private lemma entry_refl_V (ε J h : ℝ) (hε : ε = 1 ∨ ε = -1) (q : Fin N → ℝ)
    (hq : ∀ i, q i = 1 ∨ q i = -1) (s : Fin M → ℝ) (i j : Fin N) :
    |(diagonal (fun i => ε * q i) * (G.hqs J h q s - G.hqs0 J h q s)) i j| ≤ |G.Ks J s i j| := by
  rw [diagonal_mul, abs_mul, abs_refl ε hε q hq i, one_mul]
  exact entry_V G J h q s i j

/-- **Cluster confinement, operator-norm bound** (paper: Proposition `prop:cluster`(b)). -/
theorem cluster_bound (J h : ℝ) (hJ : 0 ≤ J) (hh : 2 * J < |h|) (q : Fin N → ℝ)
    (hq : ∀ i, q i = 1 ∨ q i = -1) (s : Fin M → ℝ) (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4) (t : ℝ) :
    ‖evo (G.hqs J h q s) t - evo (G.hqs0 J h q s) t‖ ≤
      8 * J / (4 * |h| - 8 * J) + 48 * |t| * J ^ 2 / (4 * |h| - 8 * J) := by
  have hg : 0 < 4 * |h| - 8 * J := by linarith
  obtain ⟨ε, hε, hεh⟩ : ∃ ε : ℝ, (ε = 1 ∨ ε = -1) ∧ ε * h = |h| := by
    rcases le_or_gt 0 h with h0 | h0
    · exact ⟨1, Or.inl rfl, by rw [abs_of_nonneg h0, one_mul]⟩
    · exact ⟨-1, Or.inr rfl, by rw [abs_of_neg h0]; ring⟩
  -- the three norm bounds (Schur test)
  have N1 : ‖cplx (diagonal (fun i => ε * q i)) * cplx (G.hqs0 J h q s) -
      ((2 * |h| : ℝ) : ℂ) • (1 : Matrix (Fin N) (Fin N) ℂ)‖ ≤ 4 * J := by
    rw [← ClusterAux.cplx_mul, ← ClusterAux.cplx_smul_one, ← ClusterAux.cplx_sub]
    exact ClusterAux.norm_cplx_le G J hJ s hs hdeg _ (entry_refl_hqs0 G ε J h hε hεh q hq s)
  have N2 : ‖cplx (G.hqs J h q s) - cplx (G.hqs0 J h q s)‖ ≤ 4 * J := by
    rw [← ClusterAux.cplx_sub]
    exact ClusterAux.norm_cplx_le G J hJ s hs hdeg _ (entry_V G J h q s)
  have N3 : ‖cplx (diagonal (fun i => ε * q i)) *
      (cplx (G.hqs J h q s) - cplx (G.hqs0 J h q s))‖ ≤ 4 * J := by
    rw [← ClusterAux.cplx_sub, ← ClusterAux.cplx_mul]
    exact ClusterAux.norm_cplx_le G J hJ s hs hdeg _ (entry_refl_V G ε J h hε q hq s)
  -- the Sylvester equation `[h⁰, S] = h¹`
  obtain ⟨S, hS, hSn⟩ := ClusterAux.exists_sylvester (cplx (diagonal (fun i => ε * q i)))
    (cplx (G.hqs0 J h q s)) (cplx (G.hqs J h q s) - cplx (G.hqs0 J h q s)) (2 * |h|) (4 * J)
    N1 (by linarith)
    (by rw [← ClusterAux.cplx_mul, refl_sq ε hε q hq, ClusterAux.cplx_one])
    (by rw [← ClusterAux.cplx_mul, ← ClusterAux.cplx_mul, refl_comm G ε J h q s])
    (by rw [← ClusterAux.cplx_sub, ← ClusterAux.cplx_mul, ← ClusterAux.cplx_mul,
      refl_anti G ε J h q hq s, ClusterAux.cplx_neg])
  -- propagators
  have hP := ClusterAux.norm_evo_sub_le _ _ (hqs_isSymm' G J h q s) (hqs0_isSymm' G J h q s)
    S hS t
  have hS0 := norm_nonneg S
  have hSle : ‖S‖ ≤ 4 * J / (4 * |h| - 8 * J) := by
    rw [le_div_iff₀ hg]
    nlinarith [hSn.trans N3]
  calc ‖evo (G.hqs J h q s) t - evo (G.hqs0 J h q s) t‖
      ≤ 2 * ‖S‖ + |t| * (‖cplx (G.hqs J h q s) - cplx (G.hqs0 J h q s)‖ * ‖S‖) := hP
    _ ≤ 2 * (4 * J / (4 * |h| - 8 * J)) + |t| * (4 * J * (4 * J / (4 * |h| - 8 * J))) := by
        gcongr
    _ = 8 * J / (4 * |h| - 8 * J) + 16 * |t| * J ^ 2 / (4 * |h| - 8 * J) := by ring
    _ ≤ 8 * J / (4 * |h| - 8 * J) + 48 * |t| * J ^ 2 / (4 * |h| - 8 * J) := by
        gcongr
        norm_num

/-- By (a), `e^{-ith⁰} δ_x` vanishes outside the cluster of `x`. -/
private lemma evo_hqs0_eq_zero (J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) (t : ℝ) (x y : Fin N)
    (hy : ¬ G.InCluster q x y) : evo (G.hqs0 J h q s) t y x = 0 := by
  classical
  let P : Matrix (Fin N) (Fin N) ℂ := diagonal (fun i => if G.InCluster q x i then 1 else 0)
  have hcomm : Commute P (cplx (G.hqs0 J h q s)) := by
    show P * cplx (G.hqs0 J h q s) = cplx (G.hqs0 J h q s) * P
    ext i j
    rw [diagonal_mul, mul_diagonal]
    by_cases h0 : G.hqs0 J h q s i j = 0
    · simp [cplx, h0]
    · have hiff : G.InCluster q x i ↔ G.InCluster q x j := by
        constructor
        · intro hi
          exact G.hqs0_cluster J h q s x j i hi (by rwa [hqs0_symm])
        · intro hj
          exact G.hqs0_cluster J h q s x i j hj h0
      by_cases hi : G.InCluster q x i
      · simp [hi, hiff.mp hi]
      · simp [hi, mt hiff.mpr hi]
  have hcomm' : Commute P (evo (G.hqs0 J h q s) t) := by
    unfold evo
    exact (hcomm.smul_right _).exp_right
  have key := congrFun (congrFun hcomm'.eq y) x
  rw [diagonal_mul, mul_diagonal] at key
  have hx : G.InCluster q x x := Relation.ReflTransGen.refl
  simpa [hy, hx] using key.symm

open Classical in
/-- **Leakage out of the cluster** (paper: Proposition `prop:cluster`(b)): the weight of
`e^{-ith} δ_x` outside `C_q(x)` is at most `(8J/g + 48|t|J²/g)²`. -/
theorem cluster_leak (J h : ℝ) (hJ : 0 ≤ J) (hh : 2 * J < |h|) (q : Fin N → ℝ)
    (hq : ∀ i, q i = 1 ∨ q i = -1) (s : Fin M → ℝ) (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4) (t : ℝ)
    (x : Fin N) :
    ∑ y ∈ univ.filter (fun y => ¬ G.InCluster q x y), ‖evo (G.hqs J h q s) t y x‖ ^ 2 ≤
      (8 * J / (4 * |h| - 8 * J) + 48 * |t| * J ^ 2 / (4 * |h| - 8 * J)) ^ 2 := by
  have hb := G.cluster_bound J h hJ hh q hq s hs hdeg t
  calc ∑ y ∈ univ.filter (fun y => ¬ G.InCluster q x y), ‖evo (G.hqs J h q s) t y x‖ ^ 2
      = ∑ y ∈ univ.filter (fun y => ¬ G.InCluster q x y),
          ‖(evo (G.hqs J h q s) t - evo (G.hqs0 J h q s) t) y x‖ ^ 2 := by
        refine Finset.sum_congr rfl fun y hy => ?_
        rw [Matrix.sub_apply, G.evo_hqs0_eq_zero J h q s t x y (Finset.mem_filter.mp hy).2,
          sub_zero]
    _ ≤ ∑ y, ‖(evo (G.hqs J h q s) t - evo (G.hqs0 J h q s) t) y x‖ ^ 2 :=
        Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _)
          (fun _ _ _ => by positivity)
    _ ≤ ‖evo (G.hqs J h q s) t - evo (G.hqs0 J h q s) t‖ ^ 2 := ClusterAux.col_sum_le _ _
    _ ≤ _ := pow_le_pow_left₀ (norm_nonneg _) hb 2

open Classical in
/-- **Averaged confinement** (paper: Proposition `prop:cluster`(c), with the percolation estimate
left as the input `P_unif(C_q(x) ⊄ B)`): charges `q ∈ {±1}^V` encoded as `sgn ∘ q`, `q : Cfg N`,
distributed with weights `w ≤ K · 2^{-N}`. -/
theorem cluster_avg (J h : ℝ) (hJ : 0 ≤ J) (hh : 2 * J < |h|) (s : Fin M → ℝ)
    (hs : ∀ e, s e = 1 ∨ s e = -1)
    (hdeg : ∀ i, (univ.filter (fun e => G.src e = i ∨ G.tgt e = i)).card ≤ 4) (t : ℝ)
    (x : Fin N) (B : Finset (Fin N)) (w : Cfg N → ℝ) (K : ℝ) (hw0 : ∀ q, 0 ≤ w q)
    (hw1 : ∑ q, w q = 1) (hwK : ∀ q, w q ≤ K / 2 ^ N) :
    ∑ q, w q * ∑ y ∈ univ.filter (fun y => y ∉ B),
        ‖evo (G.hqs J h (fun i => sgn (q i)) s) t y x‖ ^ 2 ≤
      (8 * J / (4 * |h| - 8 * J) + 48 * |t| * J ^ 2 / (4 * |h| - 8 * J)) ^ 2 +
        K * ((univ.filter (fun q : Cfg N =>
          ∃ y, y ∉ B ∧ G.InCluster (fun i => sgn (q i)) x y)).card : ℝ) / 2 ^ N := by
  set E := (8 * J / (4 * |h| - 8 * J) + 48 * |t| * J ^ 2 / (4 * |h| - 8 * J)) ^ 2 with hE
  have hsgn : ∀ (q : Cfg N) i, sgn (q i) = 1 ∨ sgn (q i) = -1 := by
    intro q i
    unfold sgn
    split_ifs <;> simp
  have hpt : ∀ q : Cfg N, w q * ∑ y ∈ univ.filter (fun y => y ∉ B),
      ‖evo (G.hqs J h (fun i => sgn (q i)) s) t y x‖ ^ 2 ≤
      w q * E + (if (∃ y, y ∉ B ∧ G.InCluster (fun i => sgn (q i)) x y) then K / 2 ^ N
        else 0) := by
    intro q
    have hcol : ∑ y ∈ univ.filter (fun y => y ∉ B),
        ‖evo (G.hqs J h (fun i => sgn (q i)) s) t y x‖ ^ 2 ≤ 1 := by
      have hU := ClusterAux.norm_le_one_of_unitary _
        (evo_unitary _ (hqs_isSymm' G J h (fun i => sgn (q i)) s) t)
      refine (Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _)
        (fun _ _ _ => by positivity)).trans ((ClusterAux.col_sum_le _ x).trans ?_)
      exact pow_le_one₀ (norm_nonneg _) hU
    split_ifs with hbad
    · have h1 : w q * ∑ y ∈ univ.filter (fun y => y ∉ B),
          ‖evo (G.hqs J h (fun i => sgn (q i)) s) t y x‖ ^ 2 ≤ w q :=
        mul_le_of_le_one_right (hw0 q) hcol
      have h2 : 0 ≤ w q * E := mul_nonneg (hw0 q) (sq_nonneg _)
      linarith [hwK q]
    · have hsub : univ.filter (fun y => y ∉ B) ⊆
          univ.filter (fun y => ¬ G.InCluster (fun i => sgn (q i)) x y) := by
        intro y hy
        rw [Finset.mem_filter] at hy ⊢
        exact ⟨mem_univ _, fun hc => hbad ⟨y, hy.2, hc⟩⟩
      have h1 := G.cluster_leak J h hJ hh (fun i => sgn (q i)) (hsgn q) s hs hdeg t x
      rw [add_zero]
      exact mul_le_mul_of_nonneg_left ((Finset.sum_le_sum_of_subset_of_nonneg hsub
        (fun _ _ _ => by positivity)).trans h1) (hw0 q)
  calc ∑ q, w q * ∑ y ∈ univ.filter (fun y => y ∉ B),
        ‖evo (G.hqs J h (fun i => sgn (q i)) s) t y x‖ ^ 2
      ≤ ∑ q : Cfg N, (w q * E + (if (∃ y, y ∉ B ∧ G.InCluster (fun i => sgn (q i)) x y) then
          K / 2 ^ N else 0)) := Finset.sum_le_sum fun q _ => hpt q
    _ = E * ∑ q, w q + ((univ.filter (fun q : Cfg N =>
          ∃ y, y ∉ B ∧ G.InCluster (fun i => sgn (q i)) x y)).card : ℝ) * (K / 2 ^ N) := by
        rw [Finset.sum_add_distrib, ← Finset.sum_mul, Finset.sum_ite, Finset.sum_const_zero,
          add_zero, Finset.sum_const, nsmul_eq_mul]
        ring
    _ = E + K * ((univ.filter (fun q : Cfg N =>
          ∃ y, y ∉ B ∧ G.InCluster (fun i => sgn (q i)) x y)).card : ℝ) / 2 ^ N := by
        rw [hw1]
        ring

end FinGraph

end DFL
