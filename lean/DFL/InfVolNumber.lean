import DFL.InfVolTransport
import DFL.AreaLawSharp
import DFL.SlaterGeneral
import DFL.NumberEntropyBulk

/-!
# Number entropy of bulk intervals in the thermodynamic limit (paper: Theorem `thm:infvol`(e))

For a finite interval `A = [a, b] ⊂ ℤ`, a Fock configuration `n` of the matter on `ℤ` and the gauge theory
on `Λ_M = [-M, M]` with Fock matter `n|_{Λ_M}` and `z`-polarized links, let `p^{(M)}_t` be the law of `N_A`
in `Ψ^{(M)}(t)` (`numDistM`). With `K^ℤ_q(t) = U^ℤ_q(t) P_S U^ℤ_q(t)^*` (`KZ`, `S = {s : n_s = 1}`) and
`p_t := E PB(ν(K^ℤ_q(t)|_A))` (`pZ`), the Poisson-binomial law of the eigenvalues of `K^ℤ_q(t)|_A` averaged
over i.i.d. charges on `ℤ`:

* `p^{(M)}_t(k) → p_t(k)` for every `k` (`tendsto_numDistM`);
* `S_N^{(M)}(t) → S_N(t) := H(p_t)` (`tendsto_numberEntropyM`);
* `S_N(t) ≤ ½ ln(2πe(2Θ(F) + 4Σ'₁(F)² + 1/12))` for every `t` under (EC) (`numberEntropyZ_le`).

## Proof outline (auxiliary lemmas in `DFL.InfVolNumberAux`)

1. The Poisson-binomial law of the eigenvalues of a Hermitian `K` is read off the generating polynomial
   `det(1 + (X - 1) K)` (`pb_eq_pbC`, via the spectral theorem); its coefficients are polynomials in the
   entries of `K` (`continuous_pbC`, via the universal matrix over `MvPolynomial`), and it is invariant
   under relabelling the rows and columns (`pbC_submatrix`).
2. Finite propagation: the entries of `h_q^k` vanish beyond distance `k`, so
   `|U^{(M)}_q(t)_{xy}| ≤ c^{|x-y|}/|x-y|! e^c`, `c = |t|(2|J| + 2|h|)`, uniformly in `M` and `q`
   (`norm_evoZ_le_wt`).
3. Hence `K^{(M)}_q(t)_{xy} → K^ℤ_q(t)_{xy}` by dominated convergence for series (`tendsto_KM`).
4. For `[a, b] ⊆ Λ_M`, `p^{(M)}_t(k) = E PB(ν(K^{(M)}_q(t)|_A))(k)` (`numDistM_eq_integral`,
   `corrR_submatrix_eqvM`), and dominated convergence for the integral gives the limit.
-/

open Matrix Finset MeasureTheory Filter Topology

noncomputable section

namespace DFL

/-- The bulk interval `[a, b]` inside `Λ_M`, relabelled in `Fin (2M+1)`. -/
def bulkM (M : ℕ) (a b : ℤ) : Finset (Fin (2 * M + 1)) :=
  univ.filter (fun i => a ≤ (i.val : ℤ) - M ∧ (i.val : ℤ) - M ≤ b)

/-- The law of `N_A`, `A = [a, b]`, in the gauge theory on `Λ_M` with Fock matter `n|_{Λ_M}` and
`z`-polarized links. -/
def numDistM (J h : ℝ) (n : CfgZ) (a b : ℤ) (t : ℝ) (M : ℕ) (k : ℕ) : ℝ :=
  lgtNumDist J h (bulkM M a b) (fockVec (occM M n)) t k

/-- `K^ℤ_q(t)_{xy} = ∑_{s : n_s = 1} U^ℤ_q(t)_{xs} conj(U^ℤ_q(t)_{ys})`. -/
def KZ (J h : ℝ) (q n : CfgZ) (t : ℝ) (x y : ℤ) : ℂ :=
  ∑' s : ℤ, if n s = 1 then UZ J h q t x s * star (UZ J h q t y s) else 0

/-- `K^ℤ_q(t)|_A`, `A = [a, b]`. -/
def KZA (J h : ℝ) (q n : CfgZ) (t : ℝ) (a b : ℤ) : Matrix (Set.Icc a b) (Set.Icc a b) ℂ :=
  fun x y => KZ J h q n t x y

/-- `K^ℤ_q(t)|_A` is Hermitian. -/
theorem KZA_isHermitian (J h : ℝ) (q n : CfgZ) (t : ℝ) (a b : ℤ) : (KZA J h q n t a b).IsHermitian := by
  ext x y
  simp only [Matrix.conjTranspose_apply, KZA, KZ]
  rw [tsum_star]
  refine tsum_congr fun s => ?_
  split_ifs
  · rw [star_mul', star_star, mul_comm]
  · exact star_zero _

/-- `p_t(k) = E PB(ν(K^ℤ_q(t)|_A))(k)`. -/
def pZ (J h : ℝ) (n : CfgZ) (a b : ℤ) (t : ℝ) (k : ℕ) : ℝ :=
  ∫ q, poissonBinomial (KZA_isHermitian J h q n t a b).eigenvalues k ∂muZ

namespace InfVolNumberAux

/-! ### The Poisson-binomial law as a continuous function of the matrix -/

section PB

variable {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ] [DecidableEq κ]

/-- The generating polynomial `det(1 + (X - 1) K)` of a square matrix `K`. -/
def gpoly (K : Matrix ι ι ℂ) : Polynomial ℂ :=
  (1 + (Polynomial.X - 1 : Polynomial ℂ) • K.map Polynomial.C).det

/-- The (real part of the) `k`-th coefficient of `det(1 + (X - 1) K)`. -/
def pbC (K : Matrix ι ι ℂ) (k : ℕ) : ℝ := ((gpoly K).coeff k).re

/-- `det(1 + c P D Q) = det(1 + c D)` if `P Q = 1`. -/
lemma det_one_add_smul_conj {R : Type*} [CommRing R] (P Q D : Matrix ι ι R) (hPQ : P * Q = 1)
    (c : R) : (1 + c • (P * D * Q)).det = (1 + c • D).det := by
  have e : 1 + c • (P * D * Q) = P * (1 + c • D) * Q := by
    rw [Matrix.mul_add, Matrix.add_mul, Matrix.mul_one, hPQ, Matrix.mul_smul, Matrix.smul_mul]
  rw [e, Matrix.det_mul, Matrix.det_mul, mul_right_comm, ← Matrix.det_mul, hPQ, Matrix.det_one,
    one_mul]

/-- For a Hermitian `K` with eigenvalues `ν`, `det(1 + (X - 1) K) = ∏_i (ν_i X + 1 - ν_i)`. -/
lemma gpoly_eq (K : Matrix ι ι ℂ) (hK : K.IsHermitian) :
    gpoly K = (SlaterExactAux.pbPoly hK.eigenvalues).map Complex.ofRealHom := by
  set U : Matrix ι ι ℂ := (hK.eigenvectorUnitary : Matrix ι ι ℂ) with hUdef
  have hU : U * star U = 1 := Unitary.mul_star_self_of_mem hK.eigenvectorUnitary.2
  have hK' : K = U * diagonal (fun i => (hK.eigenvalues i : ℂ)) * star U :=
    hK.spectral_theorem.trans rfl
  have hPQ : U.map Polynomial.C * (star U).map Polynomial.C = (1 : Matrix ι ι (Polynomial ℂ)) := by
    rw [← Matrix.map_mul, hU, Matrix.map_one _ (map_zero _) (map_one _)]
  have hmap : K.map Polynomial.C = U.map Polynomial.C *
      (diagonal (fun i => (hK.eigenvalues i : ℂ))).map Polynomial.C *
        (star U).map Polynomial.C := by
    conv_lhs => rw [hK']
    rw [Matrix.map_mul, Matrix.map_mul]
  have hD : (1 : Matrix ι ι (Polynomial ℂ)) + (Polynomial.X - 1 : Polynomial ℂ) •
      (diagonal (fun i => (hK.eigenvalues i : ℂ))).map Polynomial.C =
      diagonal (fun i => 1 + (Polynomial.X - 1) * Polynomial.C (hK.eigenvalues i : ℂ)) := by
    rw [Matrix.diagonal_map (map_zero _)]
    ext i j
    by_cases hij : i = j
    · subst hij; simp
    · simp [hij]
  rw [gpoly, hmap, det_one_add_smul_conj _ _ _ hPQ, hD, Matrix.det_diagonal, SlaterExactAux.pbPoly,
    Polynomial.map_prod]
  refine Finset.prod_congr rfl fun i _ => ?_
  simp only [Polynomial.map_add, Polynomial.map_mul, Polynomial.map_sub, Polynomial.map_one,
    Polynomial.map_C, Polynomial.map_X, Complex.ofRealHom_eq_coe, map_sub, map_one]
  ring

/-- For a Hermitian `K`, the Poisson-binomial law of its eigenvalues is `pbC K`. -/
lemma pb_eq_pbC (K : Matrix ι ι ℂ) (hK : K.IsHermitian) (k : ℕ) :
    poissonBinomial hK.eigenvalues k = pbC K k := by
  rw [pbC, gpoly_eq K hK, Polynomial.coeff_map, SlaterExactAux.coeff_pbPoly]
  simp

/-- `det(1 + (X - 1) K)` is invariant under relabelling rows and columns. -/
lemma gpoly_submatrix (K : Matrix ι ι ℂ) (e : κ ≃ ι) : gpoly (K.submatrix e e) = gpoly K := by
  unfold gpoly
  rw [← Matrix.det_submatrix_equiv_self e
    (1 + (Polynomial.X - 1 : Polynomial ℂ) • K.map Polynomial.C)]
  congr 1
  ext i j
  simp [Matrix.one_apply]

/-- `pbC` is invariant under relabelling rows and columns. -/
lemma pbC_submatrix (K : Matrix ι ι ℂ) (e : κ ≃ ι) (k : ℕ) :
    pbC (K.submatrix e e) k = pbC K k := by
  rw [pbC, pbC, gpoly_submatrix]

/-- The coefficients of `det(1 + (X - 1) K)` are polynomials in the entries of `K`, hence continuous. -/
lemma continuous_coeff_gpoly (k : ℕ) : Continuous (fun K : Matrix ι ι ℂ => (gpoly K).coeff k) := by
  set Xu : Matrix ι ι (MvPolynomial (ι × ι) ℂ) := Matrix.of fun i j => MvPolynomial.X (i, j)
    with hXu
  set Du : Polynomial (MvPolynomial (ι × ι) ℂ) :=
    (1 + (Polynomial.X - 1 : Polynomial (MvPolynomial (ι × ι) ℂ)) • Xu.map Polynomial.C).det
    with hDu
  have key : ∀ K : Matrix ι ι ℂ,
      (gpoly K).coeff k = MvPolynomial.eval (fun p : ι × ι => K p.1 p.2) (Du.coeff k) := by
    intro K
    have h1 : Polynomial.mapRingHom (MvPolynomial.eval (fun p : ι × ι => K p.1 p.2)) Du =
        gpoly K := by
      rw [hDu, RingHom.map_det]
      unfold gpoly
      congr 1
      ext i j
      simp [hXu, Matrix.one_apply, apply_ite]
    rw [← h1, Polynomial.coe_mapRingHom, Polynomial.coeff_map]
  rw [show (fun K : Matrix ι ι ℂ => (gpoly K).coeff k) =
      fun K => MvPolynomial.eval (fun p : ι × ι => K p.1 p.2) (Du.coeff k) from funext key]
  exact (MvPolynomial.continuous_eval _).comp
    (continuous_pi fun p => continuous_id.matrix_elem p.1 p.2)

/-- `K ↦ pbC K k` is continuous. -/
lemma continuous_pbC (k : ℕ) : Continuous (fun K : Matrix ι ι ℂ => pbC K k) :=
  Complex.continuous_re.comp (continuous_coeff_gpoly k)

/-- Poisson-binomial probabilities are non-negative. -/
lemma pb_nonneg (ν : ι → ℝ) (hν : ∀ i, 0 ≤ ν i ∧ ν i ≤ 1) (k : ℕ) : 0 ≤ poissonBinomial ν k :=
  Finset.sum_nonneg fun _ _ => mul_nonneg (Finset.prod_nonneg fun i _ => (hν i).1)
    (Finset.prod_nonneg fun i _ => sub_nonneg.mpr (hν i).2)

/-- Poisson-binomial probabilities are at most one. -/
lemma pb_le_one (ν : ι → ℝ) (hν : ∀ i, 0 ≤ ν i ∧ ν i ≤ 1) (k : ℕ) : poissonBinomial ν k ≤ 1 := by
  have h1 : ∏ i, (ν i + (1 - ν i)) = 1 := by simp
  rw [Finset.prod_add, Finset.powerset_univ] at h1
  rw [← h1, poissonBinomial]
  exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _) fun T _ _ =>
    mul_nonneg (Finset.prod_nonneg fun i _ => (hν i).1)
      (Finset.prod_nonneg fun i _ => sub_nonneg.mpr (hν i).2)

/-- The Poisson-binomial law of `m` parameters vanishes beyond `m`. -/
lemma pb_eq_zero (ν : ι → ℝ) {k : ℕ} (hk : Fintype.card ι < k) : poissonBinomial ν k = 0 := by
  rw [poissonBinomial]
  refine Finset.sum_eq_zero fun T hT => ?_
  have h1 := (Finset.mem_filter.mp hT).2
  have h2 := T.card_le_univ
  omega

end PB

/-! ### Finite propagation in finite volume -/

/-- Nonzero entries of `h_q` join sites at distance at most one. -/
lemma hq_apply_eq_zero {N : ℕ} (J h : ℝ) (q : Cfg N) (i l : Fin N)
    (hil : 1 < Nat.dist i.val l.val) : hq J h q i l = 0 := by
  have h1 : ¬ (i.val + 1 = l.val ∨ l.val + 1 = i.val) := by
    unfold Nat.dist at hil; omega
  have h2 : i ≠ l := by
    rintro rfl; simp [Nat.dist_self] at hil
  simp [hq, pathAdj, h1, h2]

/-- **Finite propagation**: the entries of `h_q^k` vanish beyond distance `k`. -/
lemma pow_hq_apply_eq_zero {N : ℕ} (J h : ℝ) (q : Cfg N) (k : ℕ) :
    ∀ i j : Fin N, k < Nat.dist i.val j.val → ((hq J h q) ^ k) i j = 0 := by
  induction k with
  | zero =>
    intro i j hij
    rw [pow_zero, Matrix.one_apply_ne]
    rintro rfl
    simp [Nat.dist_self] at hij
  | succ k ih =>
    intro i j hij
    rw [pow_succ', Matrix.mul_apply]
    refine Finset.sum_eq_zero fun l _ => ?_
    by_cases hil : 1 < Nat.dist i.val l.val
    · rw [hq_apply_eq_zero J h q i l hil, zero_mul]
    · rw [ih l j (by unfold Nat.dist at hij hil ⊢; omega), mul_zero]

/-- The majorant series of the exponential tail is summable. -/
lemma summable_tail (x : ℝ) (hx : 0 ≤ x) (m : ℕ) :
    Summable (fun k : ℕ => if m ≤ k then x ^ k / (k.factorial : ℝ) else 0) := by
  refine Summable.of_nonneg_of_le (fun k => ?_) (fun k => ?_) (Real.summable_pow_div_factorial x)
  · split_ifs <;> positivity
  · split_ifs
    · exact le_rfl
    · positivity

/-- Tail of the exponential series: `∑_{k ≥ m} x^k/k! ≤ (x^m/m!) e^x`. -/
lemma tail_exp_le (x : ℝ) (hx : 0 ≤ x) (m : ℕ) :
    ∑' k : ℕ, (if m ≤ k then x ^ k / (k.factorial : ℝ) else 0) ≤
      x ^ m / (m.factorial : ℝ) * Real.exp x := by
  set f : ℕ → ℝ := fun k => if m ≤ k then x ^ k / (k.factorial : ℝ) else 0 with hf
  have hfs : Summable f := summable_tail x hx m
  rw [← hfs.sum_add_tsum_nat_add m]
  have h0 : ∑ i ∈ Finset.range m, f i = 0 := Finset.sum_eq_zero fun i hi => by
    rw [Finset.mem_range] at hi
    simp only [hf]
    rw [ite_eq_right (by omega)]
  rw [h0, zero_add]
  have hexp : Real.exp x = ∑' j : ℕ, x ^ j / (j.factorial : ℝ) := by
    rw [Real.exp_eq_exp_ℝ, NormedSpace.exp_eq_tsum_div]
  rw [hexp, ← tsum_mul_left]
  refine ((summable_nat_add_iff m).mpr hfs).tsum_le_tsum (fun j => ?_)
    ((Real.summable_pow_div_factorial x).mul_left _)
  simp only [hf]
  rw [ite_eq_left (by omega)]
  have h1 : ((m.factorial : ℝ) * j.factorial) ≤ ((j + m).factorial : ℝ) := by
    have := Nat.factorial_mul_factorial_dvd_factorial_add m j
    rw [add_comm m j] at this
    exact_mod_cast Nat.le_of_dvd (Nat.factorial_pos _) this
  rw [pow_add, div_mul_div_comm, mul_comm (x ^ j) (x ^ m)]
  exact div_le_div_of_nonneg_left (by positivity) (by positivity) h1

/-- The propagation weight `w_c(d) = c^d/d! e^c`. -/
def wt (c : ℝ) (d : ℕ) : ℝ := c ^ d / (d.factorial : ℝ) * Real.exp c

/-- The propagation weight is non-negative. -/
lemma wt_nonneg {c : ℝ} (hc : 0 ≤ c) (d : ℕ) : 0 ≤ wt c d := by
  unfold wt; positivity

/-- **Propagation bound** in finite volume: `|(e^{-ith_q})_{ij}| ≤ w_c(|i - j|)` with
`c = |t|(2|J| + 2|h|)`, uniformly in the volume and in the charges. -/
lemma norm_evo_hq_le {N : ℕ} (J h t : ℝ) (q : Cfg N) (i j : Fin N) :
    ‖evo (hq J h q) t i j‖ ≤ wt (|t| * (2 * |J| + 2 * |h|)) (Nat.dist i.val j.val) := by
  have hc : 0 ≤ |t| * (2 * |J| + 2 * |h|) := by positivity
  refine le_trans ?_ (tail_exp_le _ hc (Nat.dist i.val j.val))
  refine (InfVolAux.hasSum_evo_apply (hq J h q) t i j).norm_le_of_bounded
    (summable_tail _ hc _).hasSum fun k => ?_
  by_cases hk : Nat.dist i.val j.val ≤ k
  · rw [ite_eq_left hk]
    exact InfVolAux.norm_term_le t k _ _ (InfVolAux.abs_pow_hq_le J h q k i j)
  · rw [ite_eq_right hk, pow_hq_apply_eq_zero J h q k i j (by omega)]
    simp

/-- **Propagation bound** for the finite-volume entries embedded in `ℤ`, uniformly in `M`. -/
lemma norm_evoZ_le_wt (J h : ℝ) (q : CfgZ) (t : ℝ) (M : ℕ) (x s : ℤ) :
    ‖evoZ J h q t M x s‖ ≤ wt (|t| * (2 * |J| + 2 * |h|)) (x - s).natAbs := by
  unfold evoZ
  split_ifs with hxs
  · rw [← InfVolAux.dist_idxZ M x s hxs.1 hxs.2]
    exact norm_evo_hq_le J h t _ _ _
  · rw [norm_zero]
    exact wt_nonneg (by positivity) _

/-- `s ↦ w_c(|x - s|)` is summable over `ℤ`. -/
lemma summable_wt_int (c : ℝ) (x : ℤ) : Summable (fun s : ℤ => wt c (x - s).natAbs) := by
  have hs : Summable (fun d : ℕ => wt c d) := (Real.summable_pow_div_factorial c).mul_right _
  have h0 : Summable (fun s : ℤ => wt c s.natAbs) :=
    Summable.of_nat_of_neg (by simpa using hs) (by simpa using hs)
  have e : (fun s : ℤ => wt c (x - s).natAbs) = (fun s : ℤ => wt c s.natAbs) ∘ (Equiv.subLeft x) := by
    funext s
    simp [Equiv.subLeft_apply]
  rw [e]
  exact (Equiv.subLeft x).summable_iff.mpr h0

/-! ### Convergence of the correlation matrices -/

/-- `K^{(M)}_q(t)_{xy}`: the correlation matrix of the gauge theory on `Λ_M`, embedded in `ℤ`. -/
def KM (J h : ℝ) (q n : CfgZ) (t : ℝ) (M : ℕ) (x y : ℤ) : ℂ :=
  ∑' s : ℤ, if n s = 1 then evoZ J h q t M x s * star (evoZ J h q t M y s) else 0

/-- `K^{(M)}_q(t)|_A` on the sites of `A = [a, b]`. -/
def KMA (J h : ℝ) (q n : CfgZ) (t : ℝ) (M : ℕ) (a b : ℤ) : Matrix (Set.Icc a b) (Set.Icc a b) ℂ :=
  fun x y => KM J h q n t M x y

/-- **Convergence of the correlation matrices** (paper: Step 1 of the proof of `thm:infvol`):
`K^{(M)}_q(t)_{xy} → K^ℤ_q(t)_{xy}`. -/
lemma tendsto_KM (J h : ℝ) (q n : CfgZ) (t : ℝ) (x y : ℤ) :
    Tendsto (fun M => KM J h q n t M x y) atTop (𝓝 (KZ J h q n t x y)) := by
  have hc : 0 ≤ |t| * (2 * |J| + 2 * |h|) := by positivity
  unfold KM KZ
  refine tendsto_tsum_of_dominated_convergence
    (bound := fun s => wt (|t| * (2 * |J| + 2 * |h|)) (x - s).natAbs) (summable_wt_int _ x)
    (fun s => ?_) (Eventually.of_forall fun M s => ?_)
  · by_cases hs : n s = 1
    · simp only [hs, ite_true]
      exact (tendsto_evoZ J h q t x s).mul (tendsto_evoZ J h q t y s).star
    · simp only [hs, ite_false]
      exact tendsto_const_nhds
  · split_ifs
    · rw [norm_mul, norm_star]
      calc ‖evoZ J h q t M x s‖ * ‖evoZ J h q t M y s‖
          ≤ wt (|t| * (2 * |J| + 2 * |h|)) (x - s).natAbs * 1 :=
            mul_le_mul (norm_evoZ_le_wt J h q t M x s) (InfVolAux.norm_evoZ_le_one J h q t M y s)
              (norm_nonneg _) (wt_nonneg hc _)
        _ = _ := mul_one _
    · rw [norm_zero]
      exact wt_nonneg hc _

/-- Convergence of `K^{(M)}_q(t)|_A` to `K^ℤ_q(t)|_A`. -/
lemma tendsto_KMA (J h : ℝ) (q n : CfgZ) (t : ℝ) (a b : ℤ) :
    Tendsto (fun M => KMA J h q n t M a b) atTop (𝓝 (KZA J h q n t a b)) :=
  tendsto_pi_nhds.2 fun x => tendsto_pi_nhds.2 fun y => tendsto_KM J h q n t x y

/-! ### The finite-volume law of `N_A` -/

/-- `[a, b] ⊆ Λ_M`. -/
lemma abs_le_of_mem {a b : ℤ} {M : ℕ} (hM : a.natAbs ≤ M ∧ b.natAbs ≤ M) {x : ℤ}
    (hx : x ∈ Set.Icc a b) : |x| ≤ M := by
  obtain ⟨h1, h2⟩ := hx
  obtain ⟨h3, h4⟩ := hM
  rw [abs_le]
  constructor <;> omega

/-- For `[a, b] ⊆ Λ_M`, the relabelling `x ↦ x + M` identifies `[a, b]` with `bulkM M a b`. -/
def eqvM (M : ℕ) (a b : ℤ) (hM : a.natAbs ≤ M ∧ b.natAbs ≤ M) : Set.Icc a b ≃ ↥(bulkM M a b) where
  toFun x := ⟨idxZ M x.1 (abs_le_of_mem hM x.2), by
    obtain ⟨h1, h2⟩ := x.2
    obtain ⟨h3, h4⟩ := hM
    simp only [bulkM, Finset.mem_filter, Finset.mem_univ, true_and, idxZ]
    omega⟩
  invFun i := ⟨(i.1.val : ℤ) - M, (Finset.mem_filter.mp i.2).2⟩
  left_inv x := by
    obtain ⟨h1, h2⟩ := x.2
    obtain ⟨h3, h4⟩ := hM
    apply Subtype.ext
    simp only [idxZ]
    omega
  right_inv i := by
    apply Subtype.ext
    apply Fin.ext
    simp only [idxZ]
    omega

/-- The finite-volume entries at the sites `i - M` of `Λ_M`. -/
lemma evoZ_val_sub (J h : ℝ) (q : CfgZ) (t : ℝ) (M : ℕ) (x : ℤ) (hx : |x| ≤ M)
    (i : Fin (2 * M + 1)) :
    evoZ J h q t M x ((i.val : ℤ) - M) = evo (hq J h (restrictZ M q)) t (idxZ M x hx) i := by
  have hi := InfVolTransportAux.abs_val_sub_le M i
  unfold evoZ
  rw [dite_eq_left ⟨hx, hi⟩, InfVolTransportAux.idxZ_val_sub]

/-- A series over `ℤ` supported in `Λ_M` is the finite sum over `Fin (2M+1)`. -/
lemma tsum_eq_sum_fin (M : ℕ) (g : ℤ → ℂ) (hg : ∀ s : ℤ, ¬ |s| ≤ M → g s = 0) :
    ∑' s, g s = ∑ i : Fin (2 * M + 1), g ((i.val : ℤ) - M) := by
  have hinj : Function.Injective (fun i : Fin (2 * M + 1) => ((i.val : ℤ) - M)) := by
    intro i j hij
    simp only at hij
    ext
    omega
  rw [← hinj.tsum_eq, tsum_fintype]
  intro s hs
  rw [Function.mem_support] at hs
  have hs' : |s| ≤ M := by
    by_contra hc
    exact hs (hg s hc)
  refine ⟨idxZ M s hs', ?_⟩
  have := abs_le.mp hs'
  simp only [idxZ]
  omega

/-- For `[a, b] ⊆ Λ_M`, `C_R` of the bulk interval, relabelled by `[a, b]`, is `K^{(M)}_q(t)|_A`. -/
lemma corrR_submatrix_eqvM (J h : ℝ) (q n : CfgZ) (t : ℝ) (M : ℕ) (a b : ℤ)
    (hM : a.natAbs ≤ M ∧ b.natAbs ≤ M) :
    (corrR (hq J h (restrictZ M q)) (occM M n) t (bulkM M a b)).submatrix (eqvM M a b hM)
      (eqvM M a b hM) = KMA J h q n t M a b := by
  ext x y
  have hx := abs_le_of_mem hM x.2
  have hy := abs_le_of_mem hM y.2
  simp only [Matrix.submatrix_apply, corrR, KMA, KM]
  rw [FreeDynAux.corrMat_apply, tsum_eq_sum_fin M _ ?_, occM, Finset.sum_filter]
  · refine Finset.sum_congr rfl fun i _ => ?_
    rw [evoZ_val_sub J h q t M x hx i, evoZ_val_sub J h q t M y hy i]
    rfl
  · intro s hs
    simp [evoZ, hs]

/-- The finite-volume law of `N_A` as an integral over the charges on `ℤ`. -/
lemma numDistM_eq_integral (J h : ℝ) (n : CfgZ) (a b : ℤ) (t : ℝ) (M k : ℕ) :
    numDistM J h n a b t M k = ∫ q, poissonBinomial
      (corrR_isHermitian (hq J h (restrictZ M q)) (occM M n) t (bulkM M a b)).eigenvalues k ∂muZ := by
  calc numDistM J h n a b t M k
      = avgQ (fun q' : Cfg (2 * M + 1) => poissonBinomial
          (corrR_isHermitian (hq J h q') (occM M n) t (bulkM M a b)).eigenvalues k) := by
        unfold numDistM
        rw [TransportSharpAux.lgtNumDist_eq_avgQ]
        congr 1
        funext q'
        exact numDist_slater_gen _ (hq_isSymm J h q') _ t _ k
    _ = _ := (InfVolAux.integral_comp_restrictZ M _).symm

end InfVolNumberAux

open InfVolNumberAux in
/-- **Convergence of the law of `N_A`** (paper: Theorem `thm:infvol`(e)). -/
theorem tendsto_numDistM (J h : ℝ) (n : CfgZ) (a b : ℤ) (t : ℝ) (k : ℕ) :
    Tendsto (fun M : ℕ => numDistM J h n a b t M k) atTop (𝓝 (pZ J h n a b t k)) := by
  refine Tendsto.congr (fun M => (numDistM_eq_integral J h n a b t M k).symm) ?_
  unfold pZ
  refine tendsto_integral_of_dominated_convergence (fun _ => (1 : ℝ)) (fun M => ?_)
    (integrable_const _) (fun M => ae_of_all _ fun q => ?_) (ae_of_all _ fun q => ?_)
  · exact (InfVolAux.measurable_comp_restrictZ M (fun q' => poissonBinomial
      (corrR_isHermitian (hq J h q') (occM M n) t (bulkM M a b)).eigenvalues k)).aestronglyMeasurable
  · have hν := eig_corrR_mem _ (hq_isSymm J h (restrictZ M q)) (occM M n) t (bulkM M a b)
    rw [Real.norm_eq_abs, abs_le]
    exact ⟨by linarith [pb_nonneg _ hν k], pb_le_one _ hν k⟩
  · rw [pb_eq_pbC _ (KZA_isHermitian J h q n t a b)]
    refine (((continuous_pbC (ι := Set.Icc a b) k).tendsto (KZA J h q n t a b)).comp
      (tendsto_KMA J h q n t a b)).congr' ?_
    filter_upwards [eventually_ge_atTop (max a.natAbs b.natAbs)] with M hM
    have hM' : a.natAbs ≤ M ∧ b.natAbs ≤ M := ⟨le_of_max_le_left hM, le_of_max_le_right hM⟩
    simp only [Function.comp_apply]
    rw [pb_eq_pbC _ (corrR_isHermitian _ _ _ _),
      ← pbC_submatrix (corrR (hq J h (restrictZ M q)) (occM M n) t (bulkM M a b)) (eqvM M a b hM'),
      corrR_submatrix_eqvM]

open InfVolNumberAux in
/-- **Convergence of the number entropy** (paper: Theorem `thm:infvol`(e)):
`S_N^{(M)}(t) → H(p_t)`. -/
theorem tendsto_numberEntropyM (J h : ℝ) (n : CfgZ) (a b : ℤ) (hab : a ≤ b) (t : ℝ) :
    Tendsto (fun M : ℕ => shannon (2 * M + 1) (numDistM J h n a b t M)) atTop
      (𝓝 (shannon (b - a + 1).toNat (pZ J h n a b t))) := by
  -- the hypothesis `a ≤ b` is not needed: for `a > b` both laws are `δ_0`
  have _ := hab
  set L := (b - a + 1).toNat with hL
  have hcard : Fintype.card (Set.Icc a b) = L := by
    rw [Fintype.card_Icc, Int.card_Icc, hL]
    congr 1
    ring
  have hnum : ∀ M : ℕ, max a.natAbs b.natAbs ≤ M → ∀ k, L < k →
      numDistM J h n a b t M k = 0 := by
    intro M hM k hk
    have hM' : a.natAbs ≤ M ∧ b.natAbs ≤ M := ⟨le_of_max_le_left hM, le_of_max_le_right hM⟩
    have hc : Fintype.card ↥(bulkM M a b) = L := by
      rw [← Fintype.card_congr (eqvM M a b hM'), hcard]
    rw [numDistM_eq_integral]
    have h0 : ∀ q : CfgZ, poissonBinomial (corrR_isHermitian (hq J h (restrictZ M q)) (occM M n) t
        (bulkM M a b)).eigenvalues k = 0 :=
      fun q => pb_eq_zero _ (by rw [hc]; exact hk)
    simp only [h0, integral_zero]
  have hlim : Tendsto (fun M : ℕ => ∑ k ∈ range (L + 1),
      Real.negMulLog (numDistM J h n a b t M k)) atTop (𝓝 (shannon L (pZ J h n a b t))) :=
    tendsto_finsetSum _ fun k _ =>
      (Real.continuous_negMulLog.tendsto _).comp (tendsto_numDistM J h n a b t k)
  refine hlim.congr' ?_
  filter_upwards [eventually_ge_atTop (max (max a.natAbs b.natAbs) L)] with M hM
  have hM1 : max a.natAbs b.natAbs ≤ M := le_of_max_le_left hM
  have hM2 : L ≤ M := le_of_max_le_right hM
  unfold shannon
  refine Finset.sum_subset (Finset.range_subset_range.mpr (by omega)) fun k hk1 hk2 => ?_
  rw [hnum M hM1 k (by simp only [Finset.mem_range] at hk2; omega), Real.negMulLog_zero]

/-- **Bounded number entropy in the thermodynamic limit** (paper: Theorem `thm:infvol`(e)):
`S_N(t) ≤ ½ ln(2πe(2Θ(F) + 4Σ'₁(F)² + 1/12))` for every `t`. -/
theorem numberEntropyZ_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (n : CfgZ) (a b : ℤ) (hab : a ≤ b)
    (t : ℝ) :
    shannon (b - a + 1).toNat (pZ J h n a b t) ≤
      (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 *
        (2 * thetaF F + 4 * (∑' d : ℕ, (d : ℝ) * √(F d)) ^ 2 + 1 / 12)) := by
  refine le_of_tendsto (tendsto_numberEntropyM J h n a b hab t) ?_
  filter_upwards [eventually_ge_atTop (max a.natAbs b.natAbs)] with M hM
  have hM' : a.natAbs ≤ M ∧ b.natAbs ≤ M := ⟨le_of_max_le_left hM, le_of_max_le_right hM⟩
  have hb : bulkM M a b = ivl (a + M).toNat (b + M).toNat := by
    obtain ⟨h3, h4⟩ := hM'
    ext i
    simp only [bulkM, ivl, Finset.mem_filter, Finset.mem_univ, true_and]
    omega
  have := lgt_numberEntropy_bulk J h hEC (occM M n) (a + M).toNat (b + M).toNat t
  rw [← hb] at this
  exact this

end DFL
