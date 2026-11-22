import DFL.ECProof

/-!
# Infinite volume: the one-body propagator on `ℤ` (paper: Lemma `lem:propagation`, Corollary
`cor:ECinfvol`, Theorem `thm:infvol`(b))

Two-sided chains `Λ_M = [-M, M] ≅ Fin (2M+1)` (relabelling `x ↦ x + M`) carry the finite-volume models
of the rest of the library. The charges `q ∈ {±1}^ℤ` are i.i.d. uniform (`muZ`, the infinite product of
the uniform law on `ZMod 2`), and the sector of `Λ_M` is the restriction `q|_{Λ_M}`.

The infinite-volume one-body operator `h^ℤ_q = -J Δ_ℤ + 2h diag(q)` is a tridiagonal infinite matrix;
its powers are finite sums (`hZpow`), and the matrix entries of `e^{-ith^ℤ_q}` are given by the
absolutely convergent exponential series (`UZ`). This avoids operator theory on `ℓ²(ℤ)`.

* `tendsto_evoZ` (paper: `lem:propagation`(b)): the finite-volume entries converge to those of `UZ`.
* `lintegral_iSup_UZ_le` (paper: `cor:ECinfvol`, `thm:infvol`(b)): under (EC),
  `E sup_t |U^ℤ_q(t)_{xy}| ≤ F(|x-y|)`; `infvol_propagator_of_ldt`: the same with an exponential
  profile, from the uniform large deviations.
* `tendsto_density_one`: the one-particle density converges to `E |U^ℤ_q(t)_{xy}|²`.
* `moments_UZ`: `∑_x (|x-y|+1)^p E|U^ℤ_q(t)_{xy}|² ≤ 2 M_p(F)` for `0 ≤ p ≤ 2`.
-/

open MeasureTheory Filter Topology Matrix Finset

noncomputable section

namespace DFL

/-- Charge configurations on `ℤ`. -/
abbrev CfgZ := ℤ → ZMod 2

instance : MeasurableSpace (ZMod 2) := ⊤

instance : DiscreteMeasurableSpace (ZMod 2) := ⟨fun _ => trivial⟩

/-- The uniform law on `ZMod 2`. -/
def unifZ2 : Measure (ZMod 2) := (2 : ENNReal)⁻¹ • Measure.count

instance : IsProbabilityMeasure unifZ2 := by
  constructor
  rw [unifZ2, smul_apply, Measure.count_univ, smul_eq_mul]
  have h2 : ((ENat.card (ZMod 2) : ℕ∞) : ENNReal) = 2 := by
    rw [ENat.card_eq_coe_fintype_card, ZMod.card]; rfl
  rw [h2]
  exact ENNReal.inv_mul_cancel (by norm_num) (by norm_num)

/-- i.i.d. uniform charges on `ℤ`. -/
def muZ : Measure CfgZ := Measure.infinitePi (fun _ : ℤ => unifZ2)

/-- The restriction of `q` to `Λ_M = [-M, M]`, relabelled as `Fin (2M+1)`. -/
def restrictZ (M : ℕ) (q : CfgZ) : Cfg (2 * M + 1) := fun i => q ((i.val : ℤ) - M)

/-- The site `x ∈ Λ_M` as an element of `Fin (2M+1)`. -/
def idxZ (M : ℕ) (x : ℤ) (hx : |x| ≤ M) : Fin (2 * M + 1) :=
  ⟨(x + M).toNat, by
    have h1 := (abs_le.mp hx).1
    have h2 := (abs_le.mp hx).2
    omega⟩

/-- `h^ℤ_q = -J Δ_ℤ + 2h diag(q)` as an infinite tridiagonal matrix. -/
def hZ (J h : ℝ) (q : CfgZ) (x y : ℤ) : ℝ :=
  (if x = y + 1 ∨ y = x + 1 then -J else 0) + (if x = y then 2 * h * sgn (q x) else 0)

/-- The powers of `h^ℤ_q` (finite sums: the matrix is tridiagonal). -/
def hZpow (J h : ℝ) (q : CfgZ) : ℕ → ℤ → ℤ → ℝ
  | 0 => fun x y => if x = y then 1 else 0
  | (k + 1) => fun x y => ∑ z ∈ Finset.Icc (x - 1) (x + 1), hZ J h q x z * hZpow J h q k z y

/-- The matrix entries of `U^ℤ_q(t) = e^{-ith^ℤ_q}`: `∑_k (-it)^k/k! (h^ℤ_q)^k_{xy}`. -/
def UZ (J h : ℝ) (q : CfgZ) (t : ℝ) (x y : ℤ) : ℂ :=
  ∑' k : ℕ, ((-(t : ℂ) * Complex.I) ^ k / (k.factorial : ℂ)) * ((hZpow J h q k x y : ℝ) : ℂ)

/-- The finite-volume entries `(e^{-ith^{(M)}_q})_{xy}` on `Λ_M` (`0` if `x` or `y` is outside). -/
def evoZ (J h : ℝ) (q : CfgZ) (t : ℝ) (M : ℕ) (x y : ℤ) : ℂ :=
  if hxy : |x| ≤ M ∧ |y| ≤ M then
    evo (hq J h (restrictZ M q)) t (idxZ M x hxy.1) (idxZ M y hxy.2) else 0

/-- `muZ` is a probability measure. -/
instance : IsProbabilityMeasure muZ := by
  unfold muZ; infer_instance

namespace InfVolAux

/-- The site `z ∈ ℤ` as an element of `Fin (2M+1)` (clamped to the box). -/
def toFin (M : ℕ) (z : ℤ) : Fin (2 * M + 1) := ⟨min (z + M).toNat (2 * M), by omega⟩

/-- `toFin` is the relabelling `z ↦ z + M` inside the box. -/
lemma toFin_val (M : ℕ) (z : ℤ) (hz : |z| ≤ M) : ((toFin M z).val : ℤ) = z + M := by
  have := abs_le.mp hz
  simp only [toFin]
  omega

/-- `idxZ` agrees with `toFin`. -/
lemma idxZ_eq_toFin (M : ℕ) (z : ℤ) (hz : |z| ≤ M) : idxZ M z hz = toFin M z := by
  have := abs_le.mp hz
  ext
  simp only [idxZ, toFin]
  omega

/-- `toFin` inverts the relabelling `i ↦ i - M`. -/
lemma toFin_sub (M : ℕ) (i : Fin (2 * M + 1)) : toFin M ((i.val : ℤ) - M) = i := by
  ext
  have := i.isLt
  simp only [toFin]
  omega

/-- `|sgn a| = 1`. -/
lemma abs_sgn (a : ZMod 2) : |sgn a| = 1 := by
  unfold sgn; split_ifs <;> simp

/-- At most one index of `Fin N` has a given value. -/
lemma sum_ite_val_le_one {N : ℕ} (m : ℕ) :
    ∑ l : Fin N, (if l.val = m then (1 : ℝ) else 0) ≤ 1 := by
  rw [Finset.sum_boole]
  exact_mod_cast Finset.card_le_one.mpr fun a ha b hb => by
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at ha hb
    exact Fin.ext (ha.trans hb.symm)

/-- The rows of `h_q` have `ℓ¹`-norm at most `2|J| + 2|h|`. -/
lemma sum_abs_hq_le {N : ℕ} (J h : ℝ) (q : Cfg N) (i : Fin N) :
    ∑ l, |hq J h q i l| ≤ 2 * |J| + 2 * |h| := by
  have h1 : ∀ l, |hq J h q i l| ≤ |J| * ((if l.val = i.val + 1 then (1 : ℝ) else 0) +
      (if l.val + 1 = i.val then (1 : ℝ) else 0)) + 2 * |h| * (if i = l then (1 : ℝ) else 0) := by
    intro l
    simp only [hq, pathAdj, Matrix.add_apply, Matrix.smul_apply, Matrix.diagonal_apply,
      smul_eq_mul]
    refine (abs_add_le _ _).trans (add_le_add ?_ ?_)
    · rw [abs_mul, abs_neg]
      refine mul_le_mul_of_nonneg_left ?_ (abs_nonneg _)
      split_ifs <;> simp_all
    · rw [abs_mul, abs_mul, abs_two]
      refine mul_le_mul_of_nonneg_left ?_ (by positivity)
      split_ifs <;> simp [abs_sgn]
  calc ∑ l, |hq J h q i l|
      ≤ ∑ l, (|J| * ((if l.val = i.val + 1 then (1 : ℝ) else 0) +
          (if l.val + 1 = i.val then (1 : ℝ) else 0)) + 2 * |h| * (if i = l then (1 : ℝ) else 0)) :=
        Finset.sum_le_sum fun l _ => h1 l
    _ = |J| * (∑ l : Fin N, (if l.val = i.val + 1 then (1 : ℝ) else 0) +
          ∑ l : Fin N, (if l.val + 1 = i.val then (1 : ℝ) else 0)) +
          2 * |h| * ∑ l : Fin N, (if i = l then (1 : ℝ) else 0) := by
        rw [Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum, Finset.sum_add_distrib]
    _ ≤ |J| * (1 + 1) + 2 * |h| * 1 := by
        have e1 := sum_ite_val_le_one (N := N) (i.val + 1)
        have e2 : ∑ l : Fin N, (if l.val + 1 = i.val then (1 : ℝ) else 0) ≤ 1 := by
          rcases Nat.eq_zero_or_pos i.val with h0 | h0
          · refine le_of_eq_of_le (Finset.sum_eq_zero fun l _ => ?_) zero_le_one
            rw [ite_eq_right (by omega)]
          · have := sum_ite_val_le_one (N := N) (i.val - 1)
            refine le_of_eq_of_le (Finset.sum_congr rfl fun l _ => ?_) this
            congr 1
            apply propext
            omega
        have e3 : ∑ l : Fin N, (if i = l then (1 : ℝ) else 0) = 1 := by
          rw [Finset.sum_ite_eq]; simp
        rw [e3]
        gcongr
    _ = 2 * |J| + 2 * |h| := by ring

/-- Entries of powers of a matrix whose rows have `ℓ¹`-norm at most `c`. -/
lemma abs_pow_apply_le {N : ℕ} (A : Matrix (Fin N) (Fin N) ℝ) (c : ℝ) (hc : 0 ≤ c)
    (hA : ∀ i, ∑ l, |A i l| ≤ c) (k : ℕ) (i j : Fin N) : |(A ^ k) i j| ≤ c ^ k := by
  induction k generalizing i with
  | zero =>
    rw [pow_zero, pow_zero, Matrix.one_apply]
    split_ifs <;> simp
  | succ k ih =>
    rw [pow_succ', Matrix.mul_apply, pow_succ']
    calc |∑ l, A i l * (A ^ k) l j| ≤ ∑ l, |A i l * (A ^ k) l j| :=
          Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ l, |A i l| * c ^ k := by
          refine Finset.sum_le_sum fun l _ => ?_
          rw [abs_mul]
          exact mul_le_mul_of_nonneg_left (ih l) (abs_nonneg _)
      _ = (∑ l, |A i l|) * c ^ k := by rw [Finset.sum_mul]
      _ ≤ c * c ^ k := mul_le_mul_of_nonneg_right (hA i) (by positivity)

/-- Reindexing a sum over `Fin (2M+1)` as a sum over `[-M, M]`. -/
lemma sum_fin_eq_sum_Icc {α : Type*} [AddCommMonoid α] (M : ℕ) (F : Fin (2 * M + 1) → α) :
    ∑ l, F l = ∑ z ∈ Finset.Icc (-(M : ℤ)) M, F (toFin M z) := by
  refine Finset.sum_nbij' (fun l => (l.val : ℤ) - M) (toFin M) ?_ ?_ ?_ ?_ ?_
  · intro l _
    have := l.isLt
    simp only [Finset.mem_Icc]
    omega
  · intro z _
    exact Finset.mem_coe.mpr (Finset.mem_univ _)
  · intro l _
    exact toFin_sub M l
  · intro z hz
    simp only [Finset.mem_Icc] at hz
    have := toFin_val M z (abs_le.mpr hz)
    omega
  · intro l _
    rw [toFin_sub]

/-- The entries of `h_{q|Λ_M}` are those of `h^ℤ_q` inside the box. -/
lemma hq_restrictZ_apply (J h : ℝ) (q : CfgZ) (M : ℕ) (x z : ℤ) (hx : |x| ≤ M) (hz : |z| ≤ M) :
    hq J h (restrictZ M q) (toFin M x) (toFin M z) = hZ J h q x z := by
  have hx' := toFin_val M x hx
  have hz' := toFin_val M z hz
  have e1 : ((toFin M x).val + 1 = (toFin M z).val ∨ (toFin M z).val + 1 = (toFin M x).val) ↔
      (x = z + 1 ∨ z = x + 1) := by omega
  have e2 : toFin M x = toFin M z ↔ x = z := by
    rw [Fin.ext_iff]; omega
  have e3 : ((toFin M x).val : ℤ) - M = x := by omega
  simp only [hq, hZ, pathAdj, Matrix.add_apply, Matrix.smul_apply, Matrix.diagonal_apply,
    smul_eq_mul, restrictZ, e3]
  by_cases c1 : x = z + 1 ∨ z = x + 1
  · rw [ite_eq_left (e1.mpr c1), ite_eq_left c1]
    by_cases c2 : x = z
    · rw [ite_eq_left (e2.mpr c2), ite_eq_left c2]; ring
    · rw [ite_eq_right (fun h => c2 (e2.mp h)), ite_eq_right c2]; ring
  · rw [ite_eq_right (fun h => c1 (e1.mp h)), ite_eq_right c1]
    by_cases c2 : x = z
    · rw [ite_eq_left (e2.mpr c2), ite_eq_left c2]; ring
    · rw [ite_eq_right (fun h => c2 (e2.mp h)), ite_eq_right c2]; ring

/-- `h^ℤ_q` is tridiagonal. -/
lemma hZ_eq_zero (J h : ℝ) (q : CfgZ) (x z : ℤ) (hz : z ∉ Finset.Icc (x - 1) (x + 1)) :
    hZ J h q x z = 0 := by
  simp only [Finset.mem_Icc, not_and_or, not_le] at hz
  simp only [hZ]
  rw [ite_eq_right (by omega), ite_eq_right (by omega)]
  ring

/-- One row of `h_{q|Λ_M}` acting on a vector, away from the boundary of the box. -/
lemma sum_hq_mul (J h : ℝ) (q : CfgZ) (M : ℕ) (x : ℤ) (hx : |x| + 1 ≤ M)
    (v : Fin (2 * M + 1) → ℝ) :
    ∑ l, hq J h (restrictZ M q) (toFin M x) l * v l =
      ∑ z ∈ Finset.Icc (x - 1) (x + 1), hZ J h q x z * v (toFin M z) := by
  have hx1 := abs_le.mp (show |x| ≤ (M : ℤ) by linarith [abs_nonneg x])
  have hx2 : -(M : ℤ) + 1 ≤ x ∧ x + 1 ≤ M := by
    rw [Int.abs_eq_natAbs] at hx; omega
  rw [sum_fin_eq_sum_Icc]
  have hsub : Finset.Icc (x - 1) (x + 1) ⊆ Finset.Icc (-(M : ℤ)) M := by
    intro z hz
    simp only [Finset.mem_Icc] at hz ⊢
    omega
  rw [← Finset.sum_subset hsub]
  · refine Finset.sum_congr rfl fun z hz => ?_
    have hz' : |z| ≤ M := abs_le.mpr (by simp only [Finset.mem_Icc] at hz; omega)
    rw [hq_restrictZ_apply J h q M x z (abs_le.mpr hx1) hz']
  · intro z hz hz2
    have hz' : |z| ≤ M := abs_le.mpr (by simp only [Finset.mem_Icc] at hz; omega)
    rw [hq_restrictZ_apply J h q M x z (abs_le.mpr hx1) hz', hZ_eq_zero J h q x z hz2, zero_mul]

/-- **Finite propagation**: the powers of `h_{q|Λ_M}` agree with those of `h^ℤ_q` as long as the
paths stay in the box. -/
lemma pow_hq_apply (J h : ℝ) (q : CfgZ) (M : ℕ) (y : ℤ) (hy : |y| ≤ M) (k : ℕ) :
    ∀ x : ℤ, |x| + k ≤ M →
      ((hq J h (restrictZ M q)) ^ k) (toFin M x) (toFin M y) = hZpow J h q k x y := by
  induction k with
  | zero =>
    intro x hx
    have hx' : |x| ≤ M := by simpa using hx
    have e : toFin M x = toFin M y ↔ x = y := by
      rw [Fin.ext_iff]
      have := toFin_val M x hx'
      have := toFin_val M y hy
      omega
    simp only [pow_zero, Matrix.one_apply, hZpow]
    by_cases hxy : x = y
    · rw [ite_eq_left (e.mpr hxy), ite_eq_left hxy]
    · rw [ite_eq_right (fun h => hxy (e.mp h)), ite_eq_right hxy]
  | succ k ih =>
    intro x hx
    rw [pow_succ', Matrix.mul_apply, sum_hq_mul J h q M x (by push_cast at hx; linarith)]
    simp only [hZpow]
    refine Finset.sum_congr rfl fun z hz => ?_
    rw [ih z ?_]
    simp only [Finset.mem_Icc] at hz
    push_cast at hx
    rw [Int.abs_eq_natAbs] at hx ⊢
    omega

open scoped Matrix.Norms.Operator in
/-- Entrywise exponential series of the propagator. -/
lemma hasSum_evo_apply {N : ℕ} (A : Matrix (Fin N) (Fin N) ℝ) (t : ℝ) (i j : Fin N) :
    HasSum (fun k : ℕ => ((-(t : ℂ) * Complex.I) ^ k / (k.factorial : ℂ)) *
      (((A ^ k) i j : ℝ) : ℂ)) (evo A t i j) := by
  have h := NormedSpace.exp_series_hasSum_exp' (𝕂 := ℂ) ((-(t : ℂ) * Complex.I) • cplx A)
  have h2 := Pi.hasSum.mp (Pi.hasSum.mp h i) j
  convert h2 using 1
  swap
  · rfl
  funext k
  have hc : (cplx A) ^ k = cplx (A ^ k) := by
    unfold cplx
    exact (Matrix.map_pow A Complex.ofRealHom k).symm
  rw [smul_pow, hc]
  simp only [Matrix.smul_apply, cplx, Matrix.map_apply, smul_eq_mul]
  ring

/-- The coefficients of the exponential series: `‖(-it)^k/k!‖ = |t|^k/k!`. -/
lemma norm_coef (t : ℝ) (k : ℕ) :
    ‖(-(t : ℂ) * Complex.I) ^ k / (k.factorial : ℂ)‖ = |t| ^ k / k.factorial := by
  rw [norm_div, norm_pow, norm_mul, norm_neg, Complex.norm_real, Complex.norm_I, mul_one,
    Complex.norm_natCast, Real.norm_eq_abs]

/-- A term of the exponential series with an entry bounded by `c^k`. -/
lemma norm_term_le (t : ℝ) (k : ℕ) (a c : ℝ) (ha : |a| ≤ c ^ k) :
    ‖((-(t : ℂ) * Complex.I) ^ k / (k.factorial : ℂ)) * ((a : ℝ) : ℂ)‖ ≤
      (|t| * c) ^ k / k.factorial := by
  rw [norm_mul, norm_coef, Complex.norm_real, Real.norm_eq_abs, mul_pow]
  calc |t| ^ k / k.factorial * |a| ≤ |t| ^ k / k.factorial * c ^ k :=
        mul_le_mul_of_nonneg_left ha (by positivity)
    _ = |t| ^ k * c ^ k / k.factorial := by ring

/-- Entries of the powers of `h_q`, uniformly in the volume. -/
lemma abs_pow_hq_le {N : ℕ} (J h : ℝ) (q : Cfg N) (k : ℕ) (i j : Fin N) :
    |((hq J h q) ^ k) i j| ≤ (2 * |J| + 2 * |h|) ^ k :=
  abs_pow_apply_le _ _ (by positivity) (sum_abs_hq_le J h q) k i j

/-- Entries of the powers of `h^ℤ_q`. -/
lemma abs_hZpow_le (J h : ℝ) (q : CfgZ) (k : ℕ) (x y : ℤ) :
    |hZpow J h q k x y| ≤ (2 * |J| + 2 * |h|) ^ k := by
  have hy : |y| ≤ ((x.natAbs + k + y.natAbs : ℕ) : ℤ) := by
    rw [Int.abs_eq_natAbs]; omega
  have hx : |x| + k ≤ ((x.natAbs + k + y.natAbs : ℕ) : ℤ) := by
    rw [Int.abs_eq_natAbs]; omega
  rw [← pow_hq_apply J h q _ y hy k x hx]
  exact abs_pow_hq_le J h _ k _ _

/-- The propagator entries have modulus at most one in every finite volume. -/
lemma norm_evoZ_le_one (J h : ℝ) (q : CfgZ) (t : ℝ) (M : ℕ) (x y : ℤ) :
    ‖evoZ J h q t M x y‖ ≤ 1 := by
  unfold evoZ
  split_ifs
  · exact norm_evo_le_one _ (hq_isSymm J h _) t _ _
  · rw [norm_zero]; exact zero_le_one

/-- The finite-volume entry as a function of the charges in the box. -/
def evoR (J h t : ℝ) (M : ℕ) (x y : ℤ) (q' : Cfg (2 * M + 1)) : ℂ :=
  if hxy : |x| ≤ M ∧ |y| ≤ M then evo (hq J h q') t (idxZ M x hxy.1) (idxZ M y hxy.2) else 0

/-- The restriction map is measurable. -/
lemma measurable_restrictZ (M : ℕ) : Measurable (restrictZ M) :=
  Measurable.of_eval fun _ => measurable_pi_apply _

/-- Every function of the restricted charges is measurable. -/
lemma measurable_comp_restrictZ {β : Type*} [MeasurableSpace β] (M : ℕ)
    (g : Cfg (2 * M + 1) → β) : Measurable fun q => g (restrictZ M q) :=
  (measurable_of_countable g).comp (measurable_restrictZ M)

/-- `unifZ2` gives mass `1/2` to each point. -/
lemma unifZ2_singleton (a : ZMod 2) : unifZ2 {a} = 2⁻¹ := by
  rw [unifZ2, smul_apply, Measure.count_singleton, smul_eq_mul, mul_one]

/-- The law of `q|_{Λ_M}` under `muZ` is uniform. -/
lemma muZ_restrictZ_preimage (M : ℕ) (q' : Cfg (2 * M + 1)) :
    muZ (restrictZ M ⁻¹' {q'}) = ((2 : ENNReal) ^ (2 * M + 1))⁻¹ := by
  have hset : restrictZ M ⁻¹' {q'} =
      Set.pi (↑(Finset.Icc (-(M : ℤ)) M)) (fun z => {q' (toFin M z)}) := by
    ext q
    simp only [Set.mem_preimage, Set.mem_singleton_iff, Set.mem_pi, Finset.mem_coe,
      Finset.mem_Icc]
    constructor
    · rintro rfl z hz
      have hz' : |z| ≤ M := abs_le.mpr hz
      simp only [restrictZ, toFin_val M z hz']
      ring_nf
    · intro H
      funext i
      have := i.isLt
      have h1 := H ((i.val : ℤ) - M) ⟨by omega, by omega⟩
      rw [toFin_sub] at h1
      exact h1
  rw [hset, muZ, Measure.infinitePi_pi _ (fun _ _ => MeasurableSet.of_discrete)]
  simp only [unifZ2_singleton]
  rw [Finset.prod_const, Int.card_Icc, ENNReal.inv_pow]
  congr 1
  omega

/-- Lower integrals of functions of `q|_{Λ_M}`. -/
lemma lintegral_comp_restrictZ (M : ℕ) (g : Cfg (2 * M + 1) → ENNReal) :
    ∫⁻ q, g (restrictZ M q) ∂muZ = ∑ q', g q' * ((2 : ENNReal) ^ (2 * M + 1))⁻¹ := by
  rw [← lintegral_map (measurable_of_countable g) (measurable_restrictZ M), lintegral_fintype]
  refine Finset.sum_congr rfl fun q' _ => ?_
  rw [Measure.map_apply (measurable_restrictZ M) (measurableSet_singleton q'),
    muZ_restrictZ_preimage]

/-- Integrals of functions of `q|_{Λ_M}` are the uniform averages. -/
lemma integral_comp_restrictZ (M : ℕ) (g : Cfg (2 * M + 1) → ℝ) :
    ∫ q, g (restrictZ M q) ∂muZ = avgQ g := by
  rw [← integral_map (measurable_restrictZ M).aemeasurable
    (measurable_of_countable g).aestronglyMeasurable, integral_fintype Integrable.of_finite]
  unfold avgQ
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun q' _ => ?_
  rw [Measure.real, Measure.map_apply (measurable_restrictZ M) (measurableSet_singleton q'),
    muZ_restrictZ_preimage, smul_eq_mul, ENNReal.toReal_inv, ENNReal.toReal_pow]
  norm_num

/-- The relabelling `Λ_M ≅ Fin (2M+1)` preserves distances. -/
lemma dist_idxZ (M : ℕ) (x y : ℤ) (hx : |x| ≤ M) (hy : |y| ≤ M) :
    Nat.dist (idxZ M x hx).val (idxZ M y hy).val = (x - y).natAbs := by
  have := abs_le.mp hx
  have := abs_le.mp hy
  simp only [idxZ, Nat.dist]
  omega

/-- `|z| ≤ ofReal c` in `ENNReal` from `‖z‖ ≤ c`. -/
lemma coe_nnnorm_le_ofReal {z : ℂ} {c : ℝ} (h : ‖z‖ ≤ c) : (‖z‖₊ : ENNReal) ≤ ENNReal.ofReal c := by
  rw [← ENNReal.ofReal_coe_nnreal, coe_nnNorm]
  exact ENNReal.ofReal_le_ofReal h

/-- `ofReal a ≤ |z|` in `ENNReal` from `a ≤ ‖z‖`. -/
lemma ofReal_le_coe_nnnorm {z : ℂ} {a : ℝ} (h : a ≤ ‖z‖) : ENNReal.ofReal a ≤ (‖z‖₊ : ENNReal) := by
  rw [← ENNReal.ofReal_coe_nnreal, coe_nnNorm]
  exact ENNReal.ofReal_le_ofReal h

/-- Finite-volume bound: `E sup_t |(e^{-ith_{q|Λ_M}})_{xy}| ≤ F(|x-y|)` under (EC). -/
lemma lintegral_iSup_evoZ_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (M : ℕ) (x y : ℤ)
    (hx : |x| ≤ M) (hy : |y| ≤ M) :
    ∫⁻ q, ⨆ t : ℝ, (‖evoZ J h q t M x y‖₊ : ENNReal) ∂muZ ≤ ENNReal.ofReal (F (x - y).natAbs) := by
  have hxy : |x| ≤ M ∧ |y| ≤ M := ⟨hx, hy⟩
  set i := idxZ M x hx
  set j := idxZ M y hy
  have hcorr : ∀ q' : Cfg (2 * M + 1), 0 ≤ corrQ (hq J h q') (hq_isHermitian J h q') i j :=
    fun q' => corrQ_nonneg _ _ _ _
  calc ∫⁻ q, ⨆ t : ℝ, (‖evoZ J h q t M x y‖₊ : ENNReal) ∂muZ
      = ∫⁻ q, (fun q' => ⨆ t : ℝ, (‖evoR J h t M x y q'‖₊ : ENNReal)) (restrictZ M q) ∂muZ := rfl
    _ = ∑ q', (⨆ t : ℝ, (‖evoR J h t M x y q'‖₊ : ENNReal)) * ((2 : ENNReal) ^ (2 * M + 1))⁻¹ :=
        lintegral_comp_restrictZ M (fun q' => ⨆ t : ℝ, (‖evoR J h t M x y q'‖₊ : ENNReal))
    _ ≤ ∑ q' : Cfg (2 * M + 1), ENNReal.ofReal (corrQ (hq J h q') (hq_isHermitian J h q') i j) *
          ((2 : ENNReal) ^ (2 * M + 1))⁻¹ := by
        refine Finset.sum_le_sum fun q' _ => mul_le_mul_left ?_ _
        refine iSup_le fun t => coe_nnnorm_le_ofReal ?_
        simp only [evoR, dite_eq_left hxy]
        exact norm_evo_le_corrQ _ _ t _ _
    _ = ENNReal.ofReal (avgQ fun q' => corrQ (hq J h q') (hq_isHermitian J h q') i j) := by
        unfold avgQ
        rw [ENNReal.ofReal_mul (by positivity), ENNReal.ofReal_sum_of_nonneg fun q' _ => hcorr q',
          Finset.mul_sum]
        refine Finset.sum_congr rfl fun q' _ => ?_
        rw [mul_comm]
        congr 1
        rw [one_div, ENNReal.ofReal_inv_of_pos (by positivity), ENNReal.ofReal_pow (by norm_num),
          ENNReal.ofReal_ofNat]
    _ ≤ ENNReal.ofReal (F (x - y).natAbs) := by
        refine ENNReal.ofReal_le_ofReal ?_
        have := hEC.bound (2 * M + 1) i j
        rwa [dist_idxZ] at this

end InfVolAux

open InfVolAux in
/-- The exponential series defining `UZ` converges absolutely. -/
theorem summable_UZ (J h : ℝ) (q : CfgZ) (t : ℝ) (x y : ℤ) :
    Summable (fun k : ℕ =>
      ((-(t : ℂ) * Complex.I) ^ k / (k.factorial : ℂ)) * ((hZpow J h q k x y : ℝ) : ℂ)) := by
  exact Summable.of_norm_bounded (Real.summable_pow_div_factorial (|t| * (2 * |J| + 2 * |h|)))
    fun k => norm_term_le t k _ _ (abs_hZpow_le J h q k x y)

open InfVolAux in
/-- **Finite propagation speed** (paper: Lemma `lem:propagation`(b)): the finite-volume propagator
entries converge to the entries of `U^ℤ_q(t)`. -/
theorem tendsto_evoZ (J h : ℝ) (q : CfgZ) (t : ℝ) (x y : ℤ) :
    Tendsto (fun M : ℕ => evoZ J h q t M x y) atTop (𝓝 (UZ J h q t x y)) := by
  let f : ℕ → ℕ → ℂ := fun M k => if |x| ≤ M ∧ |y| ≤ M then
      ((-(t : ℂ) * Complex.I) ^ k / (k.factorial : ℂ)) *
        ((((hq J h (restrictZ M q)) ^ k) (toFin M x) (toFin M y) : ℝ) : ℂ) else 0
  have hf : (fun M : ℕ => evoZ J h q t M x y) = fun M => ∑' k, f M k := by
    funext M
    by_cases hxy : |x| ≤ M ∧ |y| ≤ M
    · simp only [evoZ, dite_eq_left hxy, f, ite_eq_left hxy, idxZ_eq_toFin]
      exact (hasSum_evo_apply _ t _ _).tsum_eq.symm
    · simp only [evoZ, dite_eq_right hxy, f, ite_eq_right hxy, tsum_zero]
  rw [hf]
  refine tendsto_tsum_of_dominated_convergence
    (Real.summable_pow_div_factorial (|t| * (2 * |J| + 2 * |h|))) (fun k => ?_)
    (Eventually.of_forall fun M k => ?_)
  · apply tendsto_const_nhds.congr'
    filter_upwards [eventually_ge_atTop (x.natAbs + k + y.natAbs)] with M hM
    have hy : |y| ≤ M := by rw [Int.abs_eq_natAbs]; omega
    have hx : |x| + k ≤ M := by rw [Int.abs_eq_natAbs]; omega
    have hxy : |x| ≤ M ∧ |y| ≤ M := ⟨by linarith [(Nat.cast_nonneg k : (0 : ℤ) ≤ k)], hy⟩
    simp only [f, ite_eq_left hxy, pow_hq_apply J h q M y hy k x hx]
  · simp only [f]
    split_ifs
    · exact norm_term_le t k _ _ (abs_pow_hq_le J h _ k _ _)
    · rw [norm_zero]; positivity

/-- `|U^ℤ_q(t)_{xy}| ≤ 1`. -/
theorem norm_UZ_le_one (J h : ℝ) (q : CfgZ) (t : ℝ) (x y : ℤ) : ‖UZ J h q t x y‖ ≤ 1 := by
  exact le_of_tendsto' (tendsto_norm.comp (tendsto_evoZ J h q t x y)) fun M =>
    InfVolAux.norm_evoZ_le_one J h q t M x y

/-- The finite-volume entries depend measurably on the charges. -/
theorem measurable_evoZ (J h : ℝ) (t : ℝ) (M : ℕ) (x y : ℤ) :
    Measurable (fun q : CfgZ => evoZ J h q t M x y) := by
  exact InfVolAux.measurable_comp_restrictZ M (InfVolAux.evoR J h t M x y)

/-- The infinite-volume entries depend measurably on the charges. -/
theorem measurable_UZ (J h : ℝ) (t : ℝ) (x y : ℤ) :
    Measurable (fun q : CfgZ => UZ J h q t x y) := by
  exact measurable_of_tendsto_metrizable (fun M => measurable_evoZ J h t M x y)
    (tendsto_pi_nhds.mpr fun q => tendsto_evoZ J h q t x y)

/-- **Annealed dynamical localization on `ℤ`** (paper: Corollary `cor:ECinfvol`, Theorem
`thm:infvol`(b)): under Hypothesis (EC), `E sup_t |U^ℤ_q(t)_{xy}| ≤ F(|x-y|)`. -/
theorem lintegral_iSup_UZ_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (x y : ℤ) :
    ∫⁻ q, ⨆ t : ℝ, (‖UZ J h q t x y‖₊ : ENNReal) ∂muZ ≤ ENNReal.ofReal (F (x - y).natAbs) := by
  set G : ℕ → CfgZ → ENNReal := fun M q => ⨆ t : ℝ, (‖evoZ J h q t M x y‖₊ : ENNReal)
  have hpt : ∀ q, ⨆ t : ℝ, (‖UZ J h q t x y‖₊ : ENNReal) ≤ liminf (fun M => G M q) atTop := by
    intro q
    refine iSup_le fun t => ?_
    have hlim : Tendsto (fun M => (‖evoZ J h q t M x y‖₊ : ENNReal)) atTop
        (𝓝 (‖UZ J h q t x y‖₊ : ENNReal)) :=
      ((ENNReal.continuous_coe.comp continuous_nnNorm).tendsto _).comp (tendsto_evoZ J h q t x y)
    rw [← hlim.liminf_eq]
    exact liminf_le_liminf (Eventually.of_forall fun M =>
      le_iSup (fun t => (‖evoZ J h q t M x y‖₊ : ENNReal)) t)
  have hmeas : ∀ M, Measurable (G M) := fun M =>
    InfVolAux.measurable_comp_restrictZ M
      (fun q' => ⨆ t : ℝ, (‖InfVolAux.evoR J h t M x y q'‖₊ : ENNReal))
  calc ∫⁻ q, ⨆ t : ℝ, (‖UZ J h q t x y‖₊ : ENNReal) ∂muZ
      ≤ ∫⁻ q, liminf (fun M => G M q) atTop ∂muZ := lintegral_mono hpt
    _ ≤ liminf (fun M => ∫⁻ q, G M q ∂muZ) atTop := lintegral_liminf_le hmeas
    _ ≤ ENNReal.ofReal (F (x - y).natAbs) := by
        refine liminf_le_of_frequently_le' (Eventually.frequently ?_)
        filter_upwards [eventually_ge_atTop (x.natAbs + y.natAbs)] with M hM
        have hx : |x| ≤ M := by rw [Int.abs_eq_natAbs]; omega
        have hy : |y| ≤ M := by rw [Int.abs_eq_natAbs]; omega
        exact InfVolAux.lintegral_iSup_evoZ_le J h hEC M x y hx hy

/-- **Annealed dynamical localization on `ℤ` from uniform large deviations**: there are `C` and
`η > 0` with `E sup_t |U^ℤ_q(t)_{xy}| ≤ C e^{-η |x-y|}`. -/
theorem infvol_propagator_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ ∀ x y : ℤ,
      ∫⁻ q, ⨆ t : ℝ, (‖UZ J h q t x y‖₊ : ENNReal) ∂muZ ≤
        ENNReal.ofReal (C * Real.exp (-η * (x - y).natAbs)) := by
  obtain ⟨C, η, hη, hEC⟩ := ec_of_ldt J h hJ hL
  exact ⟨C, η, hη, fun x y => (lintegral_iSup_UZ_le J h hEC x y).trans
    (ENNReal.ofReal_le_ofReal (min_le_right _ _))⟩

/-- The one-particle density `⟨Ψ^{(M)}(t), n_x Ψ^{(M)}(t)⟩` of the gauge theory on `Λ_M` with one
fermion at `y` (`0` if `x` or `y` is outside `Λ_M`). -/
def densOne (J h : ℝ) (t : ℝ) (M : ℕ) (x y : ℤ) : ℝ :=
  if hxy : |x| ≤ M ∧ |y| ≤ M then
    (expect (matter (diagonal fun n => occ (n (idxZ M x hxy.1))))
      (lgtState J h (fockVec {idxZ M y hxy.2}) t)).re else 0

/-- **One-particle density in the thermodynamic limit** (paper: Theorem `thm:infvol`(b)):
`⟨Ψ^{(M)}(t), n_x Ψ^{(M)}(t)⟩ → E |U^ℤ_q(t)_{xy}|²`. -/
theorem tendsto_density_one (J h : ℝ) (t : ℝ) (x y : ℤ) :
    Tendsto (fun M : ℕ => densOne J h t M x y) atTop
      (𝓝 (∫ q, Complex.normSq (UZ J h q t x y) ∂muZ)) := by
  have hlim : Tendsto (fun M : ℕ => ∫ q, Complex.normSq (evoZ J h q t M x y) ∂muZ) atTop
      (𝓝 (∫ q, Complex.normSq (UZ J h q t x y) ∂muZ)) := by
    refine tendsto_integral_of_dominated_convergence (fun _ => (1 : ℝ)) (fun M => ?_)
      (integrable_const _) (fun M => ae_of_all _ fun q => ?_) (ae_of_all _ fun q => ?_)
    · exact (Complex.continuous_normSq.measurable.comp
        (measurable_evoZ J h t M x y)).aestronglyMeasurable
    · rw [Real.norm_eq_abs, abs_of_nonneg (Complex.normSq_nonneg _), Complex.normSq_eq_norm_sq]
      have := InfVolAux.norm_evoZ_le_one J h q t M x y
      nlinarith [norm_nonneg (evoZ J h q t M x y)]
    · exact (Complex.continuous_normSq.tendsto _).comp (tendsto_evoZ J h q t x y)
  refine hlim.congr' ?_
  filter_upwards [eventually_ge_atTop (x.natAbs + y.natAbs)] with M hM
  have hx : |x| ≤ M := by rw [Int.abs_eq_natAbs]; omega
  have hy : |y| ≤ M := by rw [Int.abs_eq_natAbs]; omega
  have hxy : |x| ≤ M ∧ |y| ≤ M := ⟨hx, hy⟩
  rw [densOne, dite_eq_left hxy, lgt_density, ← InfVolAux.integral_comp_restrictZ]
  refine integral_congr_ae (ae_of_all _ fun q => ?_)
  simp only [evoZ, dite_eq_left hxy, Finset.sum_singleton]

/-- **Moments** (paper: Theorem `thm:infvol`(b)): under (EC), for `0 ≤ p ≤ 2`,
`∑_x (|x-y|+1)^p E |U^ℤ_q(t)_{xy}|² ≤ 2 M_p(F)`. -/
theorem moments_UZ (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (t : ℝ) (y : ℤ) (p : ℝ) (hp0 : 0 ≤ p)
    (hp2 : p ≤ 2) :
    Summable (fun x : ℤ => (((x - y).natAbs : ℝ) + 1) ^ p * ∫ q, Complex.normSq (UZ J h q t x y) ∂muZ) ∧
      ∑' x : ℤ, (((x - y).natAbs : ℝ) + 1) ^ p * ∫ q, Complex.normSq (UZ J h q t x y) ∂muZ ≤
        2 * ∑' d : ℕ, ((d : ℝ) + 1) ^ p * F d := by
  set φ : ℕ → ℝ := fun d => ((d : ℝ) + 1) ^ p * F d
  have hφnn : ∀ d, 0 ≤ φ d := fun d =>
    mul_nonneg (Real.rpow_nonneg (by positivity) _) (hEC.nonneg d)
  have hφsum : Summable φ := by
    refine Summable.of_nonneg_of_le hφnn (fun d => ?_) hEC.summable
    refine mul_le_mul_of_nonneg_right ?_ (hEC.nonneg d)
    have h1 : (1 : ℝ) ≤ (d : ℝ) + 1 := by linarith [Nat.cast_nonneg (α := ℝ) d]
    calc ((d : ℝ) + 1) ^ p ≤ ((d : ℝ) + 1) ^ (2 : ℝ) := Real.rpow_le_rpow_of_exponent_le h1 hp2
      _ = ((d : ℝ) + 1) ^ 2 := by norm_cast
  -- the integrated densities
  have hu : ∀ x : ℤ, 0 ≤ ∫ q, Complex.normSq (UZ J h q t x y) ∂muZ ∧
      ∫ q, Complex.normSq (UZ J h q t x y) ∂muZ ≤ F (x - y).natAbs := by
    intro x
    refine ⟨integral_nonneg fun q => Complex.normSq_nonneg _, ?_⟩
    rw [integral_eq_lintegral_of_nonneg_ae (ae_of_all _ fun q => Complex.normSq_nonneg _)
      ((Complex.continuous_normSq.measurable.comp (measurable_UZ J h t x y)).aestronglyMeasurable)]
    refine ENNReal.toReal_le_of_le_ofReal (hEC.nonneg _)
      ((lintegral_mono fun q => ?_).trans (lintegral_iSup_UZ_le J h hEC x y))
    refine le_trans ?_ (le_iSup (fun t' => (‖UZ J h q t' x y‖₊ : ENNReal)) t)
    refine InfVolAux.ofReal_le_coe_nnnorm ?_
    rw [Complex.normSq_eq_norm_sq]
    have := norm_UZ_le_one J h q t x y
    nlinarith [norm_nonneg (UZ J h q t x y)]
  -- the sum over `ℤ` of the profile
  have hv : Summable (fun x : ℤ => φ (x - y).natAbs) ∧
      ∑' x : ℤ, φ (x - y).natAbs ≤ 2 * ∑' d, φ d := by
    have hs1 : HasSum (fun n : ℕ => (fun z : ℤ => φ z.natAbs) (n : ℤ)) (∑' d, φ d) := by
      simpa only [Int.natAbs_natCast] using hφsum.hasSum
    have hφsum1 : Summable (fun n : ℕ => φ (n + 1)) := (summable_nat_add_iff 1).mpr hφsum
    have hs2 : HasSum (fun n : ℕ => (fun z : ℤ => φ z.natAbs) (-((n : ℤ) + 1)))
        (∑' n : ℕ, φ (n + 1)) := by
      have e : ∀ n : ℕ, (-((n : ℤ) + 1)).natAbs = n + 1 := fun n => by omega
      simpa only [e] using hφsum1.hasSum
    have hs := HasSum.of_nat_of_neg_add_one (f := fun z : ℤ => φ z.natAbs) hs1 hs2
    have hcomp : (fun x : ℤ => φ (x - y).natAbs) =
        (fun z : ℤ => φ z.natAbs) ∘ (Equiv.subRight y) := by
      funext x
      simp only [Function.comp_apply, Equiv.subRight_apply]
    have hsy : HasSum (fun x : ℤ => φ (x - y).natAbs) (∑' d, φ d + ∑' n : ℕ, φ (n + 1)) := by
      rw [hcomp]
      exact (Equiv.subRight y).hasSum_iff.mpr hs
    refine ⟨hsy.summable, ?_⟩
    rw [hsy.tsum_eq, hφsum.tsum_eq_zero_add]
    linarith [hφnn 0]
  have hle : ∀ x : ℤ, (((x - y).natAbs : ℝ) + 1) ^ p * ∫ q, Complex.normSq (UZ J h q t x y) ∂muZ ≤
      φ (x - y).natAbs := fun x =>
    mul_le_mul_of_nonneg_left (hu x).2 (Real.rpow_nonneg (by positivity) _)
  have hnn : ∀ x : ℤ, 0 ≤ (((x - y).natAbs : ℝ) + 1) ^ p *
      ∫ q, Complex.normSq (UZ J h q t x y) ∂muZ := fun x =>
    mul_nonneg (Real.rpow_nonneg (by positivity) _) (hu x).1
  have hs : Summable (fun x : ℤ => (((x - y).natAbs : ℝ) + 1) ^ p *
      ∫ q, Complex.normSq (UZ J h q t x y) ∂muZ) := Summable.of_nonneg_of_le hnn hle hv.1
  exact ⟨hs, (Summable.tsum_le_tsum hle hs hv.1).trans hv.2⟩

end DFL
