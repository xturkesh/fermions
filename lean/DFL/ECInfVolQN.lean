import DFL.ECInfVolQ

/-!
# Annealed localization of the full correlator on the half-line (paper: Corollary `cor:ECinfvol`, `X = ℕ`)

On `ℓ²(ℕ)` (sites `0, 1, 2, …`, Dirichlet condition at `-1`; the paper's `ℕ = {1, 2, …}` shifted by one) let `A` be the
bounded operator with the matrix elements of `A_ℕ = Δ_ℕ + λ diag(q)` (`IsAN`; it exists, `exists_isAN`), with i.i.d.
uniform charges `q : ℕ → ZMod 2` (`muN`). The correlator `Q^ℕ(x, y) = sup_{|f| ≤ 1} |⟨δ_x, f(A) δ_y⟩|` (`corrQN`, continuous
functional calculus) satisfies `E Q^ℕ(x, y) ≤ C e^{-η|x - y|}` under the large-deviation input (`lintegral_corrQN_le`). The
finite windows `[0, L)` are exactly the chains of `DFL.Basic` with the restricted charges.
-/

open MeasureTheory Filter Topology
open scoped InnerProductSpace

noncomputable section

namespace DFL

/-- Charge configurations on the half-line. -/
abbrev CfgN := ℕ → ZMod 2

/-- I.i.d. uniform charges on the half-line. -/
def muN : Measure CfgN := Measure.infinitePi (fun _ : ℕ => unifZ2)

/-- `ℓ²(ℕ)`. -/
abbrev L2N : Type := lp (fun _ : ℕ => ℂ) 2

/-- The basis vector `δ_x` of `ℓ²(ℕ)`. -/
def deltaN (x : ℕ) : L2N := lp.single 2 x 1

/-- `A` has the matrix elements of `A_ℕ = Δ_ℕ + λ diag(q)` (Dirichlet at `-1`). -/
def IsAN (lam : ℝ) (q : CfgN) (A : L2N →L[ℂ] L2N) : Prop :=
  ∀ x y : ℕ, ⟪deltaN x, A (deltaN y)⟫_ℂ =
    ((if x = y + 1 ∨ y = x + 1 then 1 else 0) + (if x = y then lam * sgn (q x) else 0) : ℝ)

namespace ECInfVolQNAux

/-- The coordinates of a vector of `ℓ²(ℕ)` are its inner products with the `δ_x`. -/
lemma inner_deltaN_left (x : ℕ) (v : L2N) : ⟪deltaN x, v⟫_ℂ = v x := by
  rw [deltaN, lp.inner_single_left]
  simp

/-- The coordinates of `δ_y`. -/
lemma deltaN_apply (y z : ℕ) : deltaN y z = if z = y then 1 else 0 := by
  simp [deltaN, Pi.single_apply]

/-- `∑_z ‖c z ψ(e z)‖² ≤ ∑_w ‖ψ w‖²` for an injective `e` and `|c| ≤ 1`. -/
lemma summable_wtN {e : ℕ → ℕ} (he : Function.Injective e) (c : ℕ → ℂ) (hc : ∀ z, ‖c z‖ ≤ 1)
    (ψ : L2N) :
    Summable (fun z => ‖c z * ψ (e z)‖ ^ (2 : ENNReal).toReal) ∧
      ∑' z, ‖c z * ψ (e z)‖ ^ (2 : ENNReal).toReal ≤ ∑' w, ‖ψ w‖ ^ (2 : ENNReal).toReal := by
  have hp : 0 < (2 : ENNReal).toReal := by norm_num
  have hs0 := (lp.memℓp ψ).summable hp
  have hle : ∀ z, ‖c z * ψ (e z)‖ ^ (2 : ENNReal).toReal ≤ ‖ψ (e z)‖ ^ (2 : ENNReal).toReal :=
    fun z => by
      refine Real.rpow_le_rpow (norm_nonneg _) ?_ hp.le
      rw [norm_mul]
      exact mul_le_of_le_one_left (norm_nonneg _) (hc z)
  have hs : Summable (fun z => ‖c z * ψ (e z)‖ ^ (2 : ENNReal).toReal) :=
    (hs0.comp_injective he).of_nonneg_of_le (fun z => by positivity) hle
  exact ⟨hs, hs.tsum_le_tsum_of_inj e he (fun _ _ => by positivity) hle hs0⟩

/-- An injective weighted relabelling `z ↦ c z * ψ (e z)` (`|c| ≤ 1`) of a vector of `ℓ²(ℕ)` is in
`ℓ²(ℕ)`. -/
lemma memℓp_wtN {e : ℕ → ℕ} (he : Function.Injective e) (c : ℕ → ℂ) (hc : ∀ z, ‖c z‖ ≤ 1)
    (ψ : L2N) : Memℓp (fun z => c z * ψ (e z)) 2 :=
  memℓp_gen (summable_wtN he c hc ψ).1

/-- The weighted relabelling `(Tψ)(z) = c z ψ(e z)` (`e` injective) as a bounded operator on
`ℓ²(ℕ)`. -/
def wOpN {e : ℕ → ℕ} (he : Function.Injective e) (c : ℕ → ℂ) (hc : ∀ z, ‖c z‖ ≤ 1) :
    L2N →L[ℂ] L2N :=
  LinearMap.mkContinuous
    { toFun := fun ψ => ⟨fun z => c z * ψ (e z), memℓp_wtN he c hc ψ⟩
      map_add' := fun ψ φ => by
        ext z
        show c z * (ψ + φ) (e z) = c z * ψ (e z) + c z * φ (e z)
        rw [lp.coeFn_add, Pi.add_apply, mul_add]
      map_smul' := fun a ψ => by
        ext z
        simp only [lp.coeFn_smul, Pi.smul_apply, smul_eq_mul, RingHom.id_apply]
        ring }
    1 (fun ψ => by
      have hp : 0 < (2 : ENNReal).toReal := by norm_num
      rw [one_mul]
      refine lp.norm_le_of_tsum_le hp (norm_nonneg _) ?_
      exact (summable_wtN he c hc ψ).2.trans_eq (lp.norm_rpow_eq_tsum hp ψ).symm)

/-- `(Tψ)(z) = c z ψ(e z)`. -/
lemma wOpN_apply {e : ℕ → ℕ} (he : Function.Injective e) (c : ℕ → ℂ) (hc : ∀ z, ‖c z‖ ≤ 1)
    (ψ : L2N) (z : ℕ) : wOpN he c hc ψ z = c z * ψ (e z) := rfl

/-- The left shift `(Sψ)(z) = ψ(z + 1)` on `ℓ²(ℕ)`. -/
def shiftN : L2N →L[ℂ] L2N := wOpN (add_left_injective 1) (fun _ => 1) (fun _ => by simp)

/-- The multiplication operator `diag(sgn q)` on `ℓ²(ℕ)`. -/
def diagN (q : CfgN) : L2N →L[ℂ] L2N :=
  wOpN Function.injective_id (fun z => ((sgn (q z) : ℝ) : ℂ))
    (fun z => ECInfVolQAux.norm_sgn_le (q z))

/-- The operator `Δ_ℕ + λ diag(sgn q) = S + S* + λ diag(sgn q)` on `ℓ²(ℕ)`. -/
def ANop (lam : ℝ) (q : CfgN) : L2N →L[ℂ] L2N :=
  shiftN + ContinuousLinearMap.adjoint shiftN + (lam : ℂ) • diagN q

/-- `⟨δ_x, S δ_y⟩ = [x + 1 = y]`. -/
lemma inner_shiftN (x y : ℕ) :
    ⟪deltaN x, shiftN (deltaN y)⟫_ℂ = if x + 1 = y then 1 else 0 := by
  rw [inner_deltaN_left, shiftN, wOpN_apply, deltaN_apply, one_mul]

/-- `ANop` has the matrix elements of `Δ_ℕ + λ diag(sgn q)`. -/
lemma isAN_ANop (lam : ℝ) (q : CfgN) : IsAN lam q (ANop lam q) := by
  intro x y
  rw [ANop, add_apply, add_apply, smul_apply, inner_add_right, inner_add_right, inner_smul_right,
    ContinuousLinearMap.adjoint_inner_right, ← inner_conj_symm (shiftN (deltaN x)) (deltaN y),
    inner_shiftN, inner_shiftN,
    inner_deltaN_left, diagN, wOpN_apply, deltaN_apply, id]
  by_cases h1 : x + 1 = y
  · have h2 : ¬ y + 1 = x := by omega
    have h3 : ¬ x = y := by omega
    simp [h1, h2, h3]
  · by_cases h2 : y + 1 = x
    · have h3 : ¬ x = y := by omega
      simp [h1, h2, h3]
    · have h4 : ¬ (x = y + 1 ∨ y = x + 1) := by omega
      by_cases h3 : x = y
      · subst h3
        simp
      · simp [h1, h2, h3, h4]

end ECInfVolQNAux

/-- **The operator exists**. -/
theorem exists_isAN (lam : ℝ) (q : CfgN) : ∃ A : L2N →L[ℂ] L2N, IsAN lam q A :=
  ⟨ECInfVolQNAux.ANop lam q, ECInfVolQNAux.isAN_ANop lam q⟩

/-- The eigenfunction correlator `Q^ℕ(x, y) = sup_{|f| ≤ 1} |⟨δ_x, f(A) δ_y⟩|`. -/
def corrQN (A : L2N →L[ℂ] L2N) (x y : ℕ) : ℝ :=
  ⨆ (f : ℂ → ℂ) (_ : ∀ z, ‖f z‖ ≤ 1), ‖⟪deltaN x, cfc f A (deltaN y)⟫_ℂ‖

namespace ECInfVolQNAux

/-! ### The operator `A` -/

/-- The matrix entries `(A_ℕ)_{xy}` of `Δ_ℕ + λ diag(sgn q)`. -/
def entN (lam : ℝ) (q : CfgN) (x y : ℕ) : ℝ :=
  (if x = y + 1 ∨ y = x + 1 then 1 else 0) + (if x = y then lam * sgn (q x) else 0)

/-- `IsAN` in terms of `entN`. -/
lemma inner_eq_entN {lam : ℝ} {q : CfgN} {A : L2N →L[ℂ] L2N} (hA : IsAN lam q A) (x y : ℕ) :
    ⟪deltaN x, A (deltaN y)⟫_ℂ = ((entN lam q x y : ℝ) : ℂ) := hA x y

/-- `A_ℕ` is symmetric. -/
lemma entN_comm (lam : ℝ) (q : CfgN) (x y : ℕ) : entN lam q x y = entN lam q y x := by
  unfold entN
  by_cases hxy : x = y
  · subst hxy; rfl
  · rw [ite_eq_right hxy, ite_eq_right (Ne.symm hxy)]
    congr 1
    simp only [or_comm]

/-- `A_ℕ` is tridiagonal. -/
lemma entN_eq_zero (lam : ℝ) (q : CfgN) (x y : ℕ) (hxy : x + 2 ≤ y) : entN lam q x y = 0 := by
  unfold entN
  rw [ite_eq_right (by omega), ite_eq_right (by omega)]
  ring

/-- A bounded operator on `ℓ²(ℕ)` is determined by its matrix elements. -/
lemma ext_deltaN {T S : L2N →L[ℂ] L2N}
    (h : ∀ x y, ⟪deltaN x, T (deltaN y)⟫_ℂ = ⟪deltaN x, S (deltaN y)⟫_ℂ) : T = S := by
  have hv : ∀ y, T (deltaN y) = S (deltaN y) := fun y => by
    ext x
    rw [← inner_deltaN_left, ← inner_deltaN_left]
    exact h x y
  refine lp.ext_continuousLinearMap (by norm_num : (2 : ENNReal) ≠ ⊤) fun i => ?_
  refine ContinuousLinearMap.ext fun a => ?_
  simp only [ContinuousLinearMap.comp_apply, lp.singleContinuousLinearMap_apply]
  have e : lp.single 2 i a = a • deltaN i := by
    rw [deltaN, ← lp.single_smul, smul_eq_mul, mul_one]
  rw [e, map_smul, map_smul, hv i]

/-- An operator with the matrix elements of `A_ℕ` is self-adjoint. -/
lemma isSelfAdjoint_of_isAN {lam : ℝ} {q : CfgN} {A : L2N →L[ℂ] L2N} (hA : IsAN lam q A) :
    IsSelfAdjoint A := by
  rw [ContinuousLinearMap.isSelfAdjoint_iff']
  refine ext_deltaN fun x y => ?_
  rw [ContinuousLinearMap.adjoint_inner_right, ← inner_conj_symm, inner_eq_entN hA y x,
    inner_eq_entN hA x y, Complex.conj_ofReal, entN_comm]

/-- `A δ_x = ∑_{z ≤ x + 1} A_{zx} δ_z`. -/
lemma apply_deltaN {lam : ℝ} {q : CfgN} {A : L2N →L[ℂ] L2N} (hA : IsAN lam q A) (x : ℕ) :
    A (deltaN x) = ∑ z ∈ Finset.range (x + 2), ((entN lam q z x : ℝ) : ℂ) • deltaN z := by
  ext w
  rw [← inner_deltaN_left, inner_eq_entN hA w x, lp.coeFn_sum, Finset.sum_apply]
  simp only [lp.coeFn_smul, Pi.smul_apply, deltaN_apply, smul_eq_mul, mul_ite, mul_one, mul_zero]
  rw [Finset.sum_ite_eq]
  split_ifs with hw
  · rfl
  · rw [entN_comm, entN_eq_zero lam q x w (by simp only [Finset.mem_range] at hw; omega),
      Complex.ofReal_zero]

/-! ### Finite windows `[0, L)` -/

/-- The restriction of `q` to the window `[0, L)`. -/
def restrictN (L : ℕ) (q : CfgN) : Cfg L := fun i => q i.val

/-- The window matrix `h_{q|[0,L)}(-1, λ/2) = Δ_L + λ diag(q)` has the entries of `A_ℕ`. -/
lemma hq_restrictN_apply (lam : ℝ) (q : CfgN) (L i j : ℕ) (hi : i < L) (hj : j < L) :
    hq (-1) (lam / 2) (restrictN L q) ⟨i, hi⟩ ⟨j, hj⟩ = entN lam q i j := by
  simp only [hq, pathAdj, entN, restrictN, Matrix.add_apply, Matrix.smul_apply,
    Matrix.diagonal_apply, smul_eq_mul, Fin.mk.injEq]
  split_ifs <;> first | omega | ring

/-- **Finite propagation** (from the right end of the window only): `⟨δ_x, A^k δ_y⟩ = (A_L^k)_{xy}`
for the window `[0, L)` with `x + k < L` and `y < L`. -/
lemma inner_pow_deltaN {lam : ℝ} {q : CfgN} {A : L2N →L[ℂ] L2N} (hA : IsAN lam q A) (L y : ℕ)
    (hy : y < L) (k : ℕ) : ∀ (x : ℕ) (hx : x < L), x + k < L →
      ⟪deltaN x, (A ^ k) (deltaN y)⟫_ℂ =
        (((hq (-1) (lam / 2) (restrictN L q) ^ k) ⟨x, hx⟩ ⟨y, hy⟩ : ℝ) : ℂ) := by
  have hsa := ContinuousLinearMap.isSelfAdjoint_iff'.mp (isSelfAdjoint_of_isAN hA)
  induction k with
  | zero =>
    intro x hx _
    rw [pow_zero, one_apply_eq_self, inner_deltaN_left, deltaN_apply, pow_zero, Matrix.one_apply]
    simp only [Fin.mk.injEq]
    split_ifs <;> simp
  | succ k ih =>
    intro x hx hxk
    rw [pow_succ', mul_apply_eq_comp, ← ContinuousLinearMap.adjoint_inner_left, hsa,
      apply_deltaN hA x, sum_inner]
    simp only [inner_smul_left, Complex.conj_ofReal]
    rw [pow_succ', Matrix.mul_apply]
    set g : ℕ → ℝ := fun z => if hz : z < L then
      hq (-1) (lam / 2) (restrictN L q) ⟨x, hx⟩ ⟨z, hz⟩ *
        (hq (-1) (lam / 2) (restrictN L q) ^ k) ⟨z, hz⟩ ⟨y, hy⟩ else 0 with hg
    have e2 : ∑ l : Fin L, hq (-1) (lam / 2) (restrictN L q) ⟨x, hx⟩ l *
        (hq (-1) (lam / 2) (restrictN L q) ^ k) l ⟨y, hy⟩ = ∑ z ∈ Finset.range L, g z := by
      rw [← Fin.sum_univ_eq_sum_range]
      refine Finset.sum_congr rfl fun l _ => ?_
      simp only [hg, dite_eq_left l.isLt, Fin.eta]
    have e3 : ∑ z ∈ Finset.range L, g z = ∑ z ∈ Finset.range (x + 2), g z := by
      refine (Finset.sum_subset (Finset.range_subset_range.mpr (by omega)) fun z hz hz2 => ?_).symm
      simp only [Finset.mem_range, not_lt] at hz hz2
      simp only [hg, dite_eq_left hz, hq_restrictN_apply, entN_eq_zero lam q x z hz2, zero_mul]
    rw [e2, e3, Complex.ofReal_sum]
    refine Finset.sum_congr rfl fun z hz => ?_
    simp only [Finset.mem_range] at hz
    have hz' : z < L := by omega
    rw [ih z hz' (by omega)]
    simp only [hg, dite_eq_left hz', hq_restrictN_apply, entN_comm lam q x z]
    push_cast
    ring

/-- `Q^{[0,L)}(x, y)` as a function of the charges in the window (`0` unless `x, y < L`). -/
def corrQRN (J h : ℝ) (L x y : ℕ) (q' : Cfg L) : ℝ :=
  if hxy : x < L ∧ y < L then
    corrQ (hq J h q') (hq_isHermitian J h q') ⟨x, hxy.1⟩ ⟨y, hxy.2⟩ else 0

/-- `0 ≤ Q^{[0,L)}`. -/
lemma corrQRN_nonneg (J h : ℝ) (L x y : ℕ) (q' : Cfg L) : 0 ≤ corrQRN J h L x y q' := by
  unfold corrQRN
  split_ifs
  · exact corrQ_nonneg _ _ _ _
  · exact le_refl 0

/-- `Q^{[0,L)} ≤ 1`. -/
lemma corrQRN_le_one (J h : ℝ) (L x y : ℕ) (q' : Cfg L) : corrQRN J h L x y q' ≤ 1 := by
  unfold corrQRN
  split_ifs
  · exact ECProofAux.corrQ_le_one _ _ _ _
  · exact zero_le_one

/-! ### Comparison with finite volume -/

/-- The spectrum of a self-adjoint operator is real and bounded by `‖A‖ ‖1‖`. -/
lemma mem_spectrum_reN {A : L2N →L[ℂ] L2N} (hsa : IsSelfAdjoint A) {z : ℂ}
    (hz : z ∈ spectrum ℂ A) :
    z = (z.re : ℂ) ∧ |z.re| ≤ ‖A‖ * ‖(1 : L2N →L[ℂ] L2N)‖ :=
  ⟨hsa.mem_spectrum_eq_re hz,
    (Complex.abs_re_le_norm z).trans (spectrum.norm_le_norm_mul_of_mem hz)⟩

/-- `‖δ_x‖ = 1`. -/
lemma norm_deltaN (x : ℕ) : ‖deltaN x‖ = 1 := by
  rw [deltaN, lp.norm_single (by norm_num)]
  simp

/-- Matrix elements are bounded by the operator norm. -/
lemma norm_inner_deltaN_le (T : L2N →L[ℂ] L2N) (x y : ℕ) :
    ‖⟪deltaN x, T (deltaN y)⟫_ℂ‖ ≤ ‖T‖ := by
  calc ‖⟪deltaN x, T (deltaN y)⟫_ℂ‖ ≤ ‖deltaN x‖ * ‖T (deltaN y)‖ := norm_inner_le_norm _ _
    _ ≤ ‖deltaN x‖ * (‖T‖ * ‖deltaN y‖) :=
        mul_le_mul_of_nonneg_left (T.le_opNorm _) (norm_nonneg _)
    _ = ‖T‖ := by rw [norm_deltaN, norm_deltaN]; ring

/-- **Comparison with finite volume** (paper: proof of Corollary `cor:ECinfvol`): for `|f| ≤ 1`
and `ε > 0`, `|⟨δ_x, f(A) δ_y⟩| ≤ 2ε + Q^{[0,L)}(x, y)` for all large `L`. -/
lemma norm_inner_cfc_le (J h : ℝ) (hJ : J ≠ 0) (q : CfgN) (A : L2N →L[ℂ] L2N)
    (hA : IsAN (-(2 * h / J)) q A) (x y : ℕ) (f : ℂ → ℂ) (hf : ∀ z, ‖f z‖ ≤ 1) (ε : ℝ)
    (hε : 0 < ε) : ∃ L0 : ℕ, ∀ L : ℕ, L0 ≤ L →
      ‖⟪deltaN x, cfc f A (deltaN y)⟫_ℂ‖ ≤ 2 * ε + corrQRN J h L x y (restrictN L q) := by
  have hsa := isSelfAdjoint_of_isAN hA
  have hnormal : IsStarNormal A := hsa.isStarNormal
  by_cases hcont : ContinuousOn f (spectrum ℂ A)
  swap
  · refine ⟨0, fun L _ => ?_⟩
    rw [cfc_apply_of_not_continuousOn A hcont, zero_apply, inner_zero_right, norm_zero]
    have := corrQRN_nonneg J h L x y (restrictN L q)
    linarith
  -- extend `f|σ(A)` to a continuous function on `ℂ` with values in the closed unit disc
  obtain ⟨F, hF1, hF2⟩ := ContinuousMap.exists_forall_mem_restrict_eq (spectrum.isClosed A)
    (⟨fun z => f z, hcont.domRestrict⟩ : C(spectrum ℂ A, ℂ))
    (t := Metric.closedBall (0 : ℂ) 1) (fun z => by
      rw [Metric.mem_closedBall, dist_zero_right]
      exact hf z)
  have hFn : ∀ z, ‖F z‖ ≤ 1 := fun z => by simpa using hF1 z
  have hfF : cfc f A = cfc (⇑F) A := cfc_congr fun z hz => by
    have := congrArg (fun g : C(spectrum ℂ A, ℂ) => g ⟨z, hz⟩) hF2
    simp only [ContinuousMap.restrict_apply] at this
    exact this.symm
  set R := max (‖A‖ * ‖(1 : L2N →L[ℂ] L2N)‖) (2 + |(-(2 * h / J))|) with hR
  obtain ⟨P, hP⟩ := ECInfVolQAux.exists_poly_near F F.continuous R ε hε
  refine ⟨x + P.natDegree + y + 1, fun L hL => ?_⟩
  have hx : x < L := by omega
  have hy : y < L := by omega
  -- the operator part: `‖f(A) - p(A)‖ ≤ ε`
  have h1 : ‖cfc (⇑F) A - Polynomial.aeval A P‖ ≤ ε := by
    rw [← cfc_polynomial P A, ← cfc_sub (⇑F) (fun z => P.eval z) A]
    refine norm_cfc_le hε.le fun z hz => ?_
    obtain ⟨hz1, hz2⟩ := mem_spectrum_reN hsa hz
    have hz3 : z.re ∈ Set.Icc (-R) R := by
      have := le_max_left (‖A‖ * ‖(1 : L2N →L[ℂ] L2N)‖) (2 + |(-(2 * h / J))|)
      constructor <;> linarith [abs_le.mp hz2]
    have := hP z.re hz3
    rw [← hz1] at this
    exact this.le
  -- the matrix part: `|p(A)_{xy}| ≤ (1 + ε) Q^{[0,L)}(x, y)`
  have h2 : ‖⟪deltaN x, Polynomial.aeval A P (deltaN y)⟫_ℂ‖ ≤
      (1 + ε) * corrQRN J h L x y (restrictN L q) := by
    have e1 : ⟪deltaN x, Polynomial.aeval A P (deltaN y)⟫_ℂ =
        ∑ k ∈ Finset.range (P.natDegree + 1), P.coeff k *
          (((hq (-1) (-(2 * h / J) / 2) (restrictN L q) ^ k) ⟨x, hx⟩ ⟨y, hy⟩ : ℝ) : ℂ) := by
      rw [Polynomial.aeval_eq_sum_range, sum_apply, inner_sum]
      refine Finset.sum_congr rfl fun k hk => ?_
      have hk' : x + k < L := by
        rw [Finset.mem_range] at hk
        omega
      rw [smul_apply, inner_smul_right, inner_pow_deltaN hA L y hy k x hx hk']
    rw [e1]
    unfold corrQRN
    rw [dite_eq_left ⟨hx, hy⟩]
    refine ECInfVolQAux.norm_poly_apply_le J h hJ _ P (1 + ε) (fun t ht => ?_) _ _
    have ht' : t ∈ Set.Icc (-R) R := by
      have := le_max_right (‖A‖ * ‖(1 : L2N →L[ℂ] L2N)‖) (2 + |(-(2 * h / J))|)
      constructor <;> linarith [abs_le.mp ht]
    have := hP t ht'
    calc ‖P.eval (t : ℂ)‖ = ‖F t - (F t - P.eval (t : ℂ))‖ := by rw [sub_sub_cancel]
      _ ≤ ‖F t‖ + ‖F t - P.eval (t : ℂ)‖ := norm_sub_le _ _
      _ ≤ 1 + ε := add_le_add (hFn _) this.le
  have hQ1 := corrQRN_le_one J h L x y (restrictN L q)
  have hQ0 := corrQRN_nonneg J h L x y (restrictN L q)
  calc ‖⟪deltaN x, cfc f A (deltaN y)⟫_ℂ‖
      = ‖⟪deltaN x, (cfc (⇑F) A - Polynomial.aeval A P) (deltaN y)⟫_ℂ +
          ⟪deltaN x, Polynomial.aeval A P (deltaN y)⟫_ℂ‖ := by
        rw [hfF, sub_apply, inner_sub_right, sub_add_cancel]
    _ ≤ ‖⟪deltaN x, (cfc (⇑F) A - Polynomial.aeval A P) (deltaN y)⟫_ℂ‖ +
          ‖⟪deltaN x, Polynomial.aeval A P (deltaN y)⟫_ℂ‖ := norm_add_le _ _
    _ ≤ ε + (1 + ε) * corrQRN J h L x y (restrictN L q) :=
        add_le_add ((norm_inner_deltaN_le _ x y).trans h1) h2
    _ ≤ 2 * ε + corrQRN J h L x y (restrictN L q) := by nlinarith

/-- **Pointwise bound** (paper: Eq. `eq:ECinf`): `Q^ℕ(x, y) ≤ liminf_L Q^{[0,L)}(x, y)`. -/
lemma ofReal_corrQN_le_liminf (J h : ℝ) (hJ : J ≠ 0) (q : CfgN) (A : L2N →L[ℂ] L2N)
    (hA : IsAN (-(2 * h / J)) q A) (x y : ℕ) :
    ENNReal.ofReal (corrQN A x y) ≤
      liminf (fun L : ℕ => ENNReal.ofReal (corrQRN J h L x y (restrictN L q))) atTop := by
  set M := liminf (fun L : ℕ => ENNReal.ofReal (corrQRN J h L x y (restrictN L q))) atTop
    with hM
  have hM1 : M ≤ 1 := liminf_le_of_frequently_le' (Frequently.of_forall fun L =>
    ENNReal.ofReal_le_one.mpr (corrQRN_le_one J h L x y (restrictN L q)))
  have hMtop : M ≠ ⊤ := ne_top_of_le_ne_top ENNReal.one_ne_top hM1
  have hbound : corrQN A x y ≤ M.toReal := by
    refine Real.iSup_le (fun f => Real.iSup_le (fun hf => ?_) ENNReal.toReal_nonneg)
      ENNReal.toReal_nonneg
    refine le_of_forall_pos_le_add fun ε hε => ?_
    obtain ⟨L0, hL0⟩ := norm_inner_cfc_le J h hJ q A hA x y f hf (ε / 2) (by positivity)
    have hev : ∀ᶠ L in atTop, ENNReal.ofReal (‖⟪deltaN x, cfc f A (deltaN y)⟫_ℂ‖ - ε) ≤
        ENNReal.ofReal (corrQRN J h L x y (restrictN L q)) :=
      eventually_atTop.mpr ⟨L0, fun L hL => ENNReal.ofReal_le_ofReal (by linarith [hL0 L hL])⟩
    have h1 := (ENNReal.ofReal_le_iff_le_toReal hMtop).mp (le_liminf_of_le (h := hev))
    linarith
  calc ENNReal.ofReal (corrQN A x y) ≤ ENNReal.ofReal M.toReal := ENNReal.ofReal_le_ofReal hbound
    _ = M := ENNReal.ofReal_toReal hMtop

/-! ### The marginal of `muN` on a window -/

/-- The restriction map is measurable. -/
lemma measurable_restrictN (L : ℕ) : Measurable (restrictN L) :=
  Measurable.of_eval fun _ => measurable_pi_apply _

/-- Every function of the restricted charges is measurable. -/
lemma measurable_comp_restrictN {β : Type*} [MeasurableSpace β] (L : ℕ)
    (g : Cfg L → β) : Measurable fun q => g (restrictN L q) :=
  (measurable_of_countable g).comp (measurable_restrictN L)

/-- The law of `q|_{[0,L)}` under `muN` is uniform. -/
lemma muN_restrictN_preimage (L : ℕ) (q' : Cfg L) :
    muN (restrictN L ⁻¹' {q'}) = ((2 : ENNReal) ^ L)⁻¹ := by
  have hset : restrictN L ⁻¹' {q'} =
      Set.pi (↑(Finset.range L)) (fun z => {if hz : z < L then q' ⟨z, hz⟩ else 0}) := by
    ext q
    simp only [Set.mem_preimage, Set.mem_singleton_iff, Set.mem_pi, Finset.mem_coe,
      Finset.mem_range]
    constructor
    · rintro rfl z hz
      simp only [restrictN, dite_eq_left hz]
    · intro H
      funext i
      have h1 := H i.val i.isLt
      rw [dite_eq_left i.isLt] at h1
      exact h1
  rw [hset, muN, Measure.infinitePi_pi _ (fun _ _ => MeasurableSet.of_discrete)]
  simp only [InfVolAux.unifZ2_singleton]
  rw [Finset.prod_const, Finset.card_range, ENNReal.inv_pow]

/-- Lower integrals of functions of `q|_{[0,L)}`. -/
lemma lintegral_comp_restrictN (L : ℕ) (g : Cfg L → ENNReal) :
    ∫⁻ q, g (restrictN L q) ∂muN = ∑ q', g q' * ((2 : ENNReal) ^ L)⁻¹ := by
  rw [← lintegral_map (measurable_of_countable g) (measurable_restrictN L), lintegral_fintype]
  refine Finset.sum_congr rfl fun q' _ => ?_
  rw [Measure.map_apply (measurable_restrictN L) (measurableSet_singleton q'),
    muN_restrictN_preimage]

/-- Finite volume: `E Q^{[0,L)}(x, y) ≤ F(|x - y|)` under (EC). -/
lemma lintegral_corrQRN_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (L x y : ℕ)
    (hx : x < L) (hy : y < L) :
    ∫⁻ q, ENNReal.ofReal (corrQRN J h L x y (restrictN L q)) ∂muN ≤
      ENNReal.ofReal (F (Nat.dist x y)) := by
  have hxy : x < L ∧ y < L := ⟨hx, hy⟩
  set i : Fin L := ⟨x, hx⟩
  set j : Fin L := ⟨y, hy⟩
  have hcorr : ∀ q' : Cfg L, 0 ≤ corrQ (hq J h q') (hq_isHermitian J h q') i j :=
    fun q' => corrQ_nonneg _ _ _ _
  have he : ∀ q' : Cfg L,
      corrQRN J h L x y q' = corrQ (hq J h q') (hq_isHermitian J h q') i j := fun q' => by
    unfold corrQRN
    rw [dite_eq_left hxy]
  rw [lintegral_comp_restrictN L (fun q' => ENNReal.ofReal (corrQRN J h L x y q'))]
  simp_rw [he]
  calc ∑ q' : Cfg L, ENNReal.ofReal (corrQ (hq J h q') (hq_isHermitian J h q') i j) *
          ((2 : ENNReal) ^ L)⁻¹
      = ENNReal.ofReal (avgQ fun q' => corrQ (hq J h q') (hq_isHermitian J h q') i j) := by
        unfold avgQ
        rw [ENNReal.ofReal_mul (by positivity), ENNReal.ofReal_sum_of_nonneg fun q' _ => hcorr q',
          Finset.mul_sum]
        refine Finset.sum_congr rfl fun q' _ => ?_
        rw [mul_comm]
        congr 1
        rw [one_div, ENNReal.ofReal_inv_of_pos (by positivity), ENNReal.ofReal_pow (by norm_num),
          ENNReal.ofReal_ofNat]
    _ ≤ ENNReal.ofReal (F (Nat.dist x y)) := ENNReal.ofReal_le_ofReal (hEC.bound L i j)

end ECInfVolQNAux

/-- **Annealed localization on the half-line** (paper: Corollary `cor:ECinfvol`, `X = ℕ`):
`E Q^ℕ(x, y) ≤ C e^{-η|x-y|}`. -/
theorem lintegral_corrQN_le (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ ∀ (A : CfgN → L2N →L[ℂ] L2N), (∀ q, IsAN (-(2 * h / J)) q (A q)) →
      ∀ x y : ℕ, ∫⁻ q, ENNReal.ofReal (corrQN (A q) x y) ∂muN ≤
        ENNReal.ofReal (C * Real.exp (-η * Nat.dist x y)) := by
  obtain ⟨C, η, hη, hEC⟩ := ec_of_ldt J h hJ hL
  refine ⟨C, η, hη, fun A hA x y => ?_⟩
  have hmeas : ∀ L : ℕ, Measurable fun q : CfgN =>
      ENNReal.ofReal (ECInfVolQNAux.corrQRN J h L x y (ECInfVolQNAux.restrictN L q)) := fun L =>
    ECInfVolQNAux.measurable_comp_restrictN L
      (fun q' => ENNReal.ofReal (ECInfVolQNAux.corrQRN J h L x y q'))
  calc ∫⁻ q, ENNReal.ofReal (corrQN (A q) x y) ∂muN
      ≤ ∫⁻ q, liminf (fun L : ℕ =>
          ENNReal.ofReal (ECInfVolQNAux.corrQRN J h L x y (ECInfVolQNAux.restrictN L q)))
            atTop ∂muN :=
        lintegral_mono fun q => ECInfVolQNAux.ofReal_corrQN_le_liminf J h hJ q (A q) (hA q) x y
    _ ≤ liminf (fun L : ℕ =>
          ∫⁻ q, ENNReal.ofReal (ECInfVolQNAux.corrQRN J h L x y (ECInfVolQNAux.restrictN L q))
            ∂muN) atTop :=
        lintegral_liminf_le hmeas
    _ ≤ ENNReal.ofReal (min 1 (C * Real.exp (-η * Nat.dist x y))) := by
        refine liminf_le_of_frequently_le' (Eventually.frequently ?_)
        filter_upwards [eventually_ge_atTop (x + y + 1)] with L hLxy
        exact ECInfVolQNAux.lintegral_corrQRN_le J h hEC L x y (by omega) (by omega)
    _ ≤ ENNReal.ofReal (C * Real.exp (-η * Nat.dist x y)) :=
        ENNReal.ofReal_le_ofReal (min_le_right _ _)

end DFL
