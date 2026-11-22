import DFL.DualityFull
import DFL.SectorNorm
import DFL.AreaLawSharp

/-!
# Bounded transport and number entropy with the constants of the paper
(paper: Theorems `thm:transport`(a)–(d), `thm:numberentropy`)

Constants: `Σ₁(F) = ∑_{d ≥ 1} d F(d)` (`sigma1`), `Σ'₁(F) = ∑_{d ≥ 1} d F(d)^{1/2}` (`sigma1'`),
`Θ(F) = 2 ∑_{k ≥ 1} k (2k+1) F(k)` (`thetaF`), `M_p(F) = ∑_{d ≥ 0} (d+1)^p F(d)` (`momentF`).

* (a) **Propagator**: `{c̃_x(t), c̃_y†} = ∑_q U_q(t)_{xy} P_q` (`dressed_anticomm_t`), so
  `‖P_q {c̃_x(t), c̃_y†} P_q‖ = |U_q(t)_{xy}|` (`norm_sector_anticomm`),
  `E_q sup_t ‖P_q {c̃_x(t), c̃_y†} P_q‖ ≤ F(|x-y|)` (`avg_sup_anticomm_le`) and
  `⟨Ψ₀, |{c̃_x(t), c̃_y†}| Ψ₀⟩ ≤ F(|x-y|)` (`expect_abs_anticomm_le`).
* (b) **One particle**: `∑_x (|x-y|+1)^p ⟨n_x⟩_t ≤ E_q sup_t ∑_x (|x-y|+1)^p |U_q(t)_{xy}|² ≤ 2 M_p(F)`
  for `0 ≤ p ≤ 2` (`one_particle_moments`).
* (c) **Transport across a cut**: `|⟨N_{≤m}⟩_{q,t} - s_m| ≤ ∑_{x ≤ m < z} |U_q(t)_{xz}|²`
  (`sector_transport_le`), `E_q sup_t |⟨N_{≤m}⟩_{q,t} - s_m| ≤ Σ₁(F)` (`avg_sup_transport_le`),
  `|⟨Ψ(t), N_{≤m} Ψ(t)⟩ - s_m| ≤ Σ₁(F)` (`lgt_transport_sharp`) and
  `⟨Ψ(t), (N_{≤m} - s_m)² Ψ(t)⟩ ≤ Θ(F) + Σ'₁(F)²` (`lgt_transport_sq`).
* (d) **Dressed correlations**: `|⟨Ψ(t), c̃_x† c̃_y Ψ(t)⟩| ≤ 2 ∑_{k ≥ ⌈|x-y|/2⌉} F(k)` (`dressed_corr_le`).
* **Number entropy**: `S_N(t) ≤ ½ ln(2πe(Θ(F) + Σ'₁(F)² + 1/12))` (`lgt_numberEntropy_sharp`) and, in a
  single sector, `E_q sup_t H(p^q_t) ≤ ½ ln(2πe(Θ(F) + 1/12))` (`avg_sup_sector_numberEntropy_le`).

Sites are 0-indexed: `leftOf m = {0, …, m}` is the paper's `[1, m+1]`.
-/

open Matrix Finset
open scoped Kronecker Matrix.Norms.L2Operator ComplexOrder

noncomputable section

namespace DFL

variable {N : ℕ}

/-- `Σ₁(F) = ∑_{d ≥ 1} d F(d)`. -/
def sigma1 (F : ℕ → ℝ) : ℝ := ∑' d : ℕ, (d : ℝ) * F d

/-- `Σ'₁(F) = ∑_{d ≥ 1} d F(d)^{1/2}`. -/
def sigma1' (F : ℕ → ℝ) : ℝ := ∑' d : ℕ, (d : ℝ) * √(F d)

/-- `M_p(F) = ∑_{d ≥ 0} (d+1)^p F(d)`. -/
def momentF (p : ℝ) (F : ℕ → ℝ) : ℝ := ∑' d : ℕ, ((d : ℝ) + 1) ^ p * F d

/-- The dressed fermion at time `t`, `c̃_x(t) = e^{itH} c̃_x e^{-itH}`. -/
def dressedT (J h : ℝ) (x : Fin N) (t : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  (evo (H J h) t)ᴴ * cplx (dressed x) * evo (H J h) t

/-- The anticommutator `{c̃_x(t), c̃_y†}`. -/
def dressedAnticomm (J h : ℝ) (x y : Fin N) (t : ℝ) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ :=
  dressedT J h x t * cplx (dressed y)ᵀ + cplx (dressed y)ᵀ * dressedT J h x t

namespace TransportSharpAux

/-! ### Auxiliary lemmas -/

section Cplx

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

omit [DecidableEq ι] in
/-- `cplx` is multiplicative. -/
lemma cplx_mul' (A B : Matrix ι ι ℝ) : cplx (A * B) = cplx A * cplx B := by
  ext i j; simp [cplx, Matrix.mul_apply]

omit [Fintype ι] [DecidableEq ι] in
/-- `cplx` is additive. -/
lemma cplx_add' (A B : Matrix ι ι ℝ) : cplx (A + B) = cplx A + cplx B := by
  ext i j; simp [cplx]

omit [Fintype ι] in
/-- `cplx 1 = 1`. -/
lemma cplx_one' : cplx (1 : Matrix ι ι ℝ) = 1 := by
  ext i j; by_cases h : i = j <;> simp [cplx, Matrix.one_apply, h]

omit [Fintype ι] [DecidableEq ι] in
/-- `cplx 0 = 0`. -/
lemma cplx_zero' : cplx (0 : Matrix ι ι ℝ) = 0 := by
  ext i j; simp [cplx]

/-- `e^{-itA}` is a symmetric matrix when `A` is symmetric. -/
lemma evo_transpose (A : Matrix ι ι ℝ) (hA : A.IsSymm) (t : ℝ) : (evo A t)ᵀ = evo A t := by
  unfold evo
  rw [← Matrix.exp_transpose, transpose_smul]
  congr 2
  unfold cplx
  rw [← Matrix.transpose_map, hA.eq]

end Cplx

/-- A real matrix commuting with every Gauss operator only connects configurations with the same
charges. -/
lemma charge_eq_of_comm_gauss (M : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ)
    (hM : ∀ k, M * gauss k = gauss k * M) (a b : Cfg N × Cfg N) (hab : M a b ≠ 0) :
    chargeOf a.1 a.2 = chargeOf b.1 b.2 := by
  funext k
  have h := congrFun (congrFun (hM k) a) b
  rw [GaussLawAux.gauss_eq, mul_diagonal, diagonal_mul] at h
  by_contra hne
  have hs := sec_sgn_ne hne
  apply hab
  have h0 : M a b * (sgn (chargeOf b.1 b.2 k) - sgn (chargeOf a.1 a.2 k)) = 0 := by
    rw [mul_sub, h]; ring
  rcases mul_eq_zero.mp h0 with h1 | h1
  · exact h1
  · exact absurd (sub_eq_zero.mp h1).symm hs

/-- A real matrix commuting with every Gauss operator commutes with every sector projection. -/
lemma gaussProj_comm_of_comm_gauss (M : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ)
    (hM : ∀ k, M * gauss k = gauss k * M) (q : Cfg N) :
    gaussProj q * M = M * gaussProj q := by
  unfold gaussProj
  exact GaussLawAux.diagonal_comm_of_support _ _ fun a b hab => by
    rw [charge_eq_of_comm_gauss M hM a b hab]

/-- `c̃_j†` commutes with every Gauss operator. -/
lemma dressed_transpose_comm_gauss (j k : Fin N) :
    (dressed j)ᵀ * gauss k = gauss k * (dressed j)ᵀ := by
  have h := congrArg Matrix.transpose (dressed_comm_gauss j k)
  rw [transpose_mul, transpose_mul, (gauss_isSymm k).eq] at h
  exact h.symm

/-- `P_q` commutes with `c̃_y†`. -/
lemma cplx_gaussProj_comm_dressed_transpose (q : Cfg N) (y : Fin N) :
    cplx (gaussProj q) * cplx (dressed y)ᵀ = cplx (dressed y)ᵀ * cplx (gaussProj q) := by
  rw [← cplx_mul', ← cplx_mul',
    gaussProj_comm_of_comm_gauss _ (dressed_transpose_comm_gauss y) q]

/-- The CAR `{c̃_z, c̃_y†} = δ_{zy}`, complexified. -/
lemma cplx_dressed_CAR (z y : Fin N) :
    cplx (dressed z) * cplx (dressed y)ᵀ + cplx (dressed y)ᵀ * cplx (dressed z) =
      if z = y then 1 else 0 := by
  rw [← cplx_mul', ← cplx_mul', ← cplx_add', dressed_CAR]
  split_ifs
  · exact cplx_one'
  · exact cplx_zero'

/-- The Bernoulli average is monotone. -/
lemma avgQ_mono {f g : Cfg N → ℝ} (h : ∀ q, f q ≤ g q) : avgQ f ≤ avgQ g := by
  unfold avgQ
  exact mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun q _ => h q) (by positivity)

/-- Averaged propagator bound from Hypothesis (EC). -/
lemma avgQ_norm_evo_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (x y : Fin N) (t : ℝ) :
    avgQ (fun q : Cfg N => ‖evo (hq J h q) t x y‖) ≤ F (Nat.dist x.val y.val) :=
  (avgQ_mono fun q => norm_evo_le_corrQ (hq J h q) (hq_isHermitian J h q) t x y).trans
    (hEC.bound N x y)

/-- Uniform-in-time averaged propagator bound from Hypothesis (EC). -/
lemma avgQ_iSup_norm_evo_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (x y : Fin N) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ, ‖evo (hq J h q) t x y‖) ≤ F (Nat.dist x.val y.val) :=
  (avgQ_mono fun q => ciSup_le fun t =>
    norm_evo_le_corrQ (hq J h q) (hq_isHermitian J h q) t x y).trans (hEC.bound N x y)

/-- There are `2^N` charge configurations. -/
lemma card_cfg : Fintype.card (Cfg N) = 2 ^ N := by
  simp [Cfg, ZMod.card]

/-- Average of a constant. -/
lemma avgQ_const (c : ℝ) : avgQ (fun _ : Cfg N => c) = c := by
  simp only [avgQ, Finset.sum_const, Finset.card_univ, card_cfg, nsmul_eq_mul]
  push_cast
  field_simp

/-- Additivity of the average. -/
lemma avgQ_add (f g : Cfg N → ℝ) : avgQ (fun q => f q + g q) = avgQ f + avgQ g := by
  simp only [avgQ, Finset.sum_add_distrib]
  ring

/-- The average commutes with finite sums. -/
lemma avgQ_sum {κ : Type*} (s : Finset κ) (f : κ → Cfg N → ℝ) :
    avgQ (fun q => ∑ i ∈ s, f i q) = ∑ i ∈ s, avgQ (f i) := by
  simp only [avgQ]
  rw [Finset.sum_comm, Finset.mul_sum]

/-- Homogeneity of the average. -/
lemma avgQ_const_mul (c : ℝ) (f : Cfg N → ℝ) : avgQ (fun q => c * f q) = c * avgQ f := by
  simp only [avgQ, ← Finset.mul_sum]
  ring

/-- Average of a function minus a constant. -/
lemma avgQ_sub_const (f : Cfg N → ℝ) (c : ℝ) : avgQ (fun q => f q - c) = avgQ f - c := by
  have : avgQ (fun q => f q - c) = avgQ (fun q => f q + (-c)) := by
    simp only [sub_eq_add_neg]
  rw [this, avgQ_add, avgQ_const]
  ring

/-- `|avg f| ≤ avg |f|`. -/
lemma abs_avgQ_le (f : Cfg N → ℝ) : |avgQ f| ≤ avgQ (fun q => |f q|) := by
  simp only [avgQ]
  rw [abs_mul, abs_of_pos (by positivity)]
  exact mul_le_mul_of_nonneg_left (Finset.abs_sum_le_sum_abs _ _) (by positivity)

/-- Real part of a complex Bernoulli average of real numbers. -/
lemma re_avg (f : Cfg N → ℝ) :
    ((1 / 2 ^ N : ℂ) * ∑ q : Cfg N, ((f q : ℝ) : ℂ)).re = avgQ f := by
  simp only [avgQ]
  rw [← Complex.ofReal_sum]
  have : (1 / 2 ^ N : ℂ) = ((1 / 2 ^ N : ℝ) : ℂ) := by push_cast; ring
  rw [this, ← Complex.ofReal_mul, Complex.ofReal_re]

/-- `‖U_{xz}‖² ≤ ‖U_{xz}‖ ≤ Q(x, z)` for the propagator of `h_q`. -/
lemma norm_evo_sq_le_corrQ (J h : ℝ) (q : Cfg N) (t : ℝ) (x z : Fin N) :
    ‖evo (hq J h q) t x z‖ ^ 2 ≤ corrQ (hq J h q) (hq_isHermitian J h q) x z := by
  have h1 := norm_evo_le_one _ (hq_isSymm J h q) t x z
  have h0 := norm_nonneg (evo (hq J h q) t x z)
  have h2 := norm_evo_le_corrQ (hq J h q) (hq_isHermitian J h q) t x z
  nlinarith

section Unitary

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Expectation of a complexified real diagonal matrix: `⟨v, diag(d) v⟩ = ∑ d_i |v_i|²`. -/
lemma expect_cplx_diagonal (d : ι → ℝ) (v : ι → ℂ) :
    expect (cplx (diagonal d)) v = ((∑ i, d i * Complex.normSq (v i) : ℝ) : ℂ) := by
  have hc : cplx (diagonal d) = diagonal (fun i => (d i : ℂ)) := by
    unfold cplx; exact Matrix.diagonal_map (by simp)
  rw [expect, hc]
  simp only [dotProduct, mulVec_diagonal, Pi.star_apply, Complex.ofReal_sum, Complex.ofReal_mul]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [Complex.normSq_eq_conj_mul_self]
  simp only [Complex.star_def]
  ring

omit [DecidableEq ι] in
/-- `∑ |v_i|² = star v ⬝ᵥ v`. -/
lemma sum_normSq_eq_dotProduct (w : ι → ℂ) :
    ((∑ i, Complex.normSq (w i) : ℝ) : ℂ) = star w ⬝ᵥ w := by
  simp only [Complex.ofReal_sum, dotProduct, Pi.star_apply]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [Complex.normSq_eq_conj_mul_self]
  rfl

/-- An isometry preserves `∑ |v_i|²`. -/
lemma sum_normSq_mulVec (U : Matrix ι ι ℂ) (hU : Uᴴ * U = 1) (v : ι → ℂ) :
    ∑ i, Complex.normSq ((U *ᵥ v) i) = ∑ i, Complex.normSq (v i) := by
  have h1 : star (U *ᵥ v) ⬝ᵥ (U *ᵥ v) = star v ⬝ᵥ v := by
    rw [star_mulVec, dotProduct_mulVec, vecMul_vecMul, hU, vecMul_one]
  exact_mod_cast (sum_normSq_eq_dotProduct (U *ᵥ v)).trans
    (h1.trans (sum_normSq_eq_dotProduct v).symm)

/-- Rows of a unitary matrix are normalized. -/
lemma sum_normSq_row (U : Matrix ι ι ℂ) (hU : U * Uᴴ = 1) (x : ι) :
    ∑ z, Complex.normSq (U x z) = 1 := by
  have h := congrFun (congrFun hU x) x
  simp only [one_apply_eq] at h
  have : ((∑ z, Complex.normSq (U x z) : ℝ) : ℂ) = 1 := by
    rw [← h, Complex.ofReal_sum]
    refine Finset.sum_congr rfl (fun z _ => ?_)
    rw [← Complex.mul_conj]
    rfl
  exact_mod_cast this

/-- Columns of a unitary matrix are normalized. -/
lemma sum_normSq_col (U : Matrix ι ι ℂ) (hU : Uᴴ * U = 1) (z : ι) :
    ∑ x, Complex.normSq (U x z) = 1 := by
  have h := congrFun (congrFun hU z) z
  simp only [one_apply_eq] at h
  have : ((∑ x, Complex.normSq (U x z) : ℝ) : ℂ) = 1 := by
    rw [← h, Complex.ofReal_sum]
    refine Finset.sum_congr rfl (fun x _ => ?_)
    rw [Complex.normSq_eq_conj_mul_self]
    rfl
  exact_mod_cast this

end Unitary

/-- `⟨n_x⟩` of an evolved Fock state, as a real sum over configurations. -/
lemma sum_occ_normSq (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (x : Fin N) (t : ℝ) :
    ∑ n, occ (n x) * Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) =
      ∑ z ∈ S, Complex.normSq (evo A t x z) := by
  have := density A hA S x t
  rw [expect_cplx_diagonal] at this
  exact_mod_cast this

/-- `⟨N_R⟩` of an evolved Fock state: `∑_{x ∈ R} ∑_{z ∈ S} |U_{xz}|²`. -/
lemma sum_count_normSq (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S R : Finset (Fin N))
    (t : ℝ) :
    ∑ n, count R n * Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) =
      ∑ x ∈ R, ∑ z ∈ S, Complex.normSq (evo A t x z) := by
  simp only [count, Finset.sum_mul]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl (fun x _ => sum_occ_normSq A hA S x t)

/-- The evolved Fock state is normalized. -/
lemma sum_normSq_evo_fock (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm)
    (S : Finset (Fin N)) (t : ℝ) :
    ∑ n, Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) = 1 := by
  rw [sum_normSq_mulVec _ (evo_unitary _ (dGamma_isSymm A hA) t), fockVec,
    Finset.sum_eq_single (fockCfg S) (fun n _ hn => by rw [Pi.single_eq_of_ne hn, map_zero])
      (fun h => absurd (mem_univ _) h), Pi.single_eq_same, map_one]

/-- LGT expectation of a diagonal matter observable, as a Bernoulli average of sector sums. -/
lemma lgt_expect_diag_re (J h : ℝ) (ψ : Cfg N → ℂ) (d : Cfg N → ℝ) (t : ℝ) :
    (expect (matter (diagonal d)) (lgtState J h ψ t)).re =
      avgQ (fun q => ∑ n, d n * Complex.normSq ((evo (dGamma (hq J h q)) t *ᵥ ψ) n)) := by
  rw [matter, sector_average_matter]
  simp_rw [expect_cplx_diagonal]
  exact re_avg _

/-- Sector expectation of `N_R` in an evolved Fock state. -/
lemma expect_numOp_re (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S R : Finset (Fin N))
    (t : ℝ) :
    (expect (cplx (numOp R)) (evo (dGamma A) t *ᵥ fockVec S)).re =
      ∑ x ∈ R, ∑ z ∈ S, Complex.normSq (evo A t x z) := by
  rw [numOp, expect_cplx_diagonal, Complex.ofReal_re, sum_count_normSq A hA S R t]

/-- **Deviation of the mean particle number** for a symmetric unitary `U`:
`|∑_{x ∈ R, z ∈ S} |U_{xz}|² - |S ∩ R|| ≤ ∑_{x ∈ R, z ∉ R} |U_{xz}|²`. -/
lemma abs_mean_sub_le_symm (U : Matrix (Fin N) (Fin N) ℂ) (hU : Uᴴ * U = 1) (hT : Uᵀ = U)
    (S R : Finset (Fin N)) :
    |∑ x ∈ R, ∑ z ∈ S, Complex.normSq (U x z) - ((S ∩ R).card : ℝ)| ≤
      ∑ x ∈ R, ∑ z ∈ univ \ R, Complex.normSq (U x z) := by
  have ha : ∀ x z, 0 ≤ Complex.normSq (U x z) := fun x z => Complex.normSq_nonneg _
  have hsym : ∀ x z, U x z = U z x := fun x z => by
    have := congrFun (congrFun hT z) x
    rwa [transpose_apply] at this
  have hcol : ∀ z,
      ∑ x ∈ R, Complex.normSq (U x z) + ∑ x ∈ univ \ R, Complex.normSq (U x z) = 1 := by
    intro z
    rw [add_comm, Finset.sum_sdiff (subset_univ R)]
    exact sum_normSq_col U hU z
  set D := ∑ x ∈ R, ∑ z ∈ univ \ R, Complex.normSq (U x z)
  set a := ∑ z ∈ S \ R, ∑ x ∈ R, Complex.normSq (U x z)
  set b := ∑ z ∈ S ∩ R, ∑ x ∈ univ \ R, Complex.normSq (U x z) with hb_def
  have hab : ∑ x ∈ R, ∑ z ∈ S, Complex.normSq (U x z) - ((S ∩ R).card : ℝ) = a - b := by
    rw [Finset.sum_comm, ← Finset.sum_inter_add_sum_sdiff S R]
    have h1 : ∑ z ∈ S ∩ R, ∑ x ∈ R, Complex.normSq (U x z) = ((S ∩ R).card : ℝ) - b := by
      rw [show ((S ∩ R).card : ℝ) = ∑ z ∈ S ∩ R, (1 : ℝ) by simp, hb_def,
        ← Finset.sum_sub_distrib]
      exact Finset.sum_congr rfl fun z _ => by linarith [hcol z]
    rw [h1]
    ring
  have ha0 : 0 ≤ a := Finset.sum_nonneg fun z _ => Finset.sum_nonneg fun x _ => ha x z
  have hb0 : 0 ≤ b := Finset.sum_nonneg fun z _ => Finset.sum_nonneg fun x _ => ha x z
  have haD : a ≤ D := by
    calc a ≤ ∑ z ∈ univ \ R, ∑ x ∈ R, Complex.normSq (U x z) :=
          Finset.sum_le_sum_of_subset_of_nonneg (Finset.sdiff_subset_sdiff (subset_univ S) le_rfl)
            (fun z _ _ => Finset.sum_nonneg fun x _ => ha x z)
      _ = D := Finset.sum_comm
  have hbD : b ≤ D := by
    calc b ≤ ∑ z ∈ R, ∑ x ∈ univ \ R, Complex.normSq (U x z) :=
          Finset.sum_le_sum_of_subset_of_nonneg Finset.inter_subset_right
            (fun z _ _ => Finset.sum_nonneg fun x _ => ha x z)
      _ = D := Finset.sum_congr rfl fun z _ => Finset.sum_congr rfl fun x _ => by rw [hsym]
  rw [hab, abs_sub_le_iff]
  constructor <;> linarith

/-! #### Counting lemmas -/

/-- Bounding `∑_{a ∈ s} g(φ a)` by the cardinalities of the fibres of `φ`. -/
lemma sum_comp_le_of_card_fiber {α : Type*} (s : Finset α) (φ : α → ℕ) (M : ℕ)
    (hφ : ∀ a ∈ s, φ a < M) (c : ℕ → ℝ)
    (hc : ∀ d, ((s.filter (fun a => φ a = d)).card : ℝ) ≤ c d)
    (g : ℕ → ℝ) (hg : ∀ d, 0 ≤ g d) :
    ∑ a ∈ s, g (φ a) ≤ ∑ d ∈ range M, c d * g d := by
  rw [← Finset.sum_fiberwise_of_maps_to' (t := range M) (fun a ha => mem_range.2 (hφ a ha))]
  refine Finset.sum_le_sum (fun d _ => ?_)
  rw [Finset.sum_const, nsmul_eq_mul]
  exact mul_le_mul_of_nonneg_right (hc d) (hg d)

/-- At most `d` pairs `(x, y)` with `x ≤ m < y` are at distance `d`. -/
lemma card_cross_fiber_le (m d : ℕ) :
    (((leftOf (N := N) m ×ˢ (univ \ leftOf m)).filter
      (fun p : Fin N × Fin N => Nat.dist p.1.val p.2.val = d)).card : ℝ) ≤ d := by
  have h : ((leftOf (N := N) m ×ˢ (univ \ leftOf m)).filter
      (fun p : Fin N × Fin N => Nat.dist p.1.val p.2.val = d)).card ≤ (range d).card := by
    apply Finset.card_le_card_of_injOn (fun p : Fin N × Fin N => m - p.1.val)
    · intro p hp
      simp only [leftOf, coe_filter, mem_product, mem_filter, mem_univ, true_and, Finset.mem_sdiff,
        Set.mem_ofPred_eq] at hp
      simp only [coe_range, Set.mem_Iio]
      unfold Nat.dist at hp
      omega
    · intro p hp p' hp' hpp'
      simp only [leftOf, coe_filter, mem_product, mem_filter, mem_univ, true_and, Finset.mem_sdiff,
        Set.mem_ofPred_eq] at hp hp'
      simp only at hpp'
      unfold Nat.dist at hp hp'
      ext <;> omega
  rw [card_range] at h
  exact_mod_cast h

/-- **Crossing pairs**: `∑_{x ≤ m < y} g(|x - y|) ≤ ∑_d d g(d)`. -/
lemma sum_cross_le (g : ℕ → ℝ) (hg : ∀ d, 0 ≤ g d)
    (hs : Summable (fun d : ℕ => (d : ℝ) * g d)) (m : ℕ) :
    ∑ x ∈ leftOf (N := N) m, ∑ y ∈ univ \ leftOf (N := N) m, g (Nat.dist x.val y.val) ≤
      ∑' d : ℕ, (d : ℝ) * g d := by
  rw [← Finset.sum_product' (f := fun x y : Fin N => g (Nat.dist x.val y.val))]
  calc _ ≤ ∑ d ∈ range N, (d : ℝ) * g d :=
        sum_comp_le_of_card_fiber _ (fun p : Fin N × Fin N => Nat.dist p.1.val p.2.val) N
          (fun p _ => by have := p.1.isLt; have := p.2.isLt; unfold Nat.dist; omega)
          (fun d => d) (card_cross_fiber_le m) g hg
    _ ≤ _ := hs.sum_le_tsum _ (fun d _ => mul_nonneg (Nat.cast_nonneg d) (hg d))

/-! #### Consequences of Hypothesis (EC) -/

/-- Summability of `d ↦ d F(d)`. -/
lemma ec_summable_mul {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) :
    Summable (fun d : ℕ => (d : ℝ) * F d) :=
  hEC.summable.of_nonneg_of_le (fun d => mul_nonneg (Nat.cast_nonneg d) (hEC.nonneg d))
    (fun d => mul_le_mul_of_nonneg_right
      (by nlinarith [(Nat.cast_nonneg d : (0 : ℝ) ≤ d)]) (hEC.nonneg d))

/-- Summability of `d ↦ d √F(d)`. -/
lemma ec_summable_mul_sqrt {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) :
    Summable (fun d : ℕ => (d : ℝ) * √(F d)) :=
  hEC.summable_sqrt.of_nonneg_of_le (fun d => mul_nonneg (Nat.cast_nonneg d) (Real.sqrt_nonneg _))
    (fun d => mul_le_mul_of_nonneg_right
      (by nlinarith [(Nat.cast_nonneg d : (0 : ℝ) ≤ d)]) (Real.sqrt_nonneg _))

/-- `E ∑_{x ≤ m < z} Q(x, z) ≤ Σ₁(F)`. -/
lemma avgQ_cross_corrQ_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (m : ℕ) :
    avgQ (fun q : Cfg N => ∑ x ∈ leftOf m, ∑ z ∈ univ \ leftOf m,
      corrQ (hq J h q) (hq_isHermitian J h q) x z) ≤ sigma1 F := by
  calc _ = ∑ x ∈ leftOf (N := N) m, ∑ z ∈ univ \ leftOf (N := N) m,
        avgQ (fun q : Cfg N => corrQ (hq J h q) (hq_isHermitian J h q) x z) := by
        rw [avgQ_sum]
        exact Finset.sum_congr rfl fun x _ => avgQ_sum _ _
    _ ≤ ∑ x ∈ leftOf (N := N) m, ∑ z ∈ univ \ leftOf (N := N) m, F (Nat.dist x.val z.val) :=
        Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun z _ => hEC.bound N x z
    _ ≤ sigma1 F := sum_cross_le F hEC.nonneg (ec_summable_mul hEC) m

/-! #### The second moment of the transported charge -/

/-- `leftOf m = segment (m + 1)`. -/
lemma leftOf_eq_segment (m : ℕ) : leftOf (N := N) m = segment (m + 1) := by
  ext x
  simp [leftOf, segment]

/-- **Cauchy–Schwarz** for the Bernoulli average: `avg(√u √v) ≤ √(avg u) √(avg v)`. -/
lemma avgQ_sqrt_mul_sqrt_le {u v : Cfg N → ℝ} (hu : ∀ q, 0 ≤ u q) (hv : ∀ q, 0 ≤ v q) :
    avgQ (fun q => √(u q) * √(v q)) ≤ √(avgQ u) * √(avgQ v) := by
  simp only [avgQ]
  have hc : (0 : ℝ) ≤ 1 / 2 ^ N := by positivity
  have h1 : (1 / 2 ^ N : ℝ) * ∑ q, √(u q) * √(v q) =
      ∑ q, √((1 / 2 ^ N : ℝ) * u q) * √((1 / 2 ^ N : ℝ) * v q) := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl (fun q _ => ?_)
    rw [Real.sqrt_mul hc, Real.sqrt_mul hc,
      show √(1 / 2 ^ N : ℝ) * √(u q) * (√(1 / 2 ^ N : ℝ) * √(v q)) =
        (√(1 / 2 ^ N : ℝ) * √(1 / 2 ^ N : ℝ)) * (√(u q) * √(v q)) by ring,
      Real.mul_self_sqrt hc]
  rw [h1, Finset.mul_sum, Finset.mul_sum]
  exact Real.sum_sqrt_mul_sqrt_le _ (fun q => mul_nonneg hc (hu q)) (fun q => mul_nonneg hc (hv q))

/-- **Minkowski-type bound** for the Bernoulli average:
`avg((∑ᵢ aᵢ)²) ≤ (∑ᵢ bᵢ)²` whenever `aᵢ ≥ 0` and `√(avg aᵢ²) ≤ bᵢ`. -/
lemma avgQ_sq_sum_le {κ : Type*} (s : Finset κ) (a : κ → Cfg N → ℝ) (b : κ → ℝ)
    (ha : ∀ i ∈ s, ∀ q, 0 ≤ a i q) (hb : ∀ i ∈ s, √(avgQ (fun q => a i q ^ 2)) ≤ b i) :
    avgQ (fun q => (∑ i ∈ s, a i q) ^ 2) ≤ (∑ i ∈ s, b i) ^ 2 := by
  have h1 : (fun q => (∑ i ∈ s, a i q) ^ 2) = fun q => ∑ i ∈ s, ∑ j ∈ s, a i q * a j q := by
    funext q; rw [sq, Finset.sum_mul_sum]
  rw [h1, avgQ_sum, sq, Finset.sum_mul_sum]
  refine Finset.sum_le_sum (fun i hi => ?_)
  rw [avgQ_sum]
  refine Finset.sum_le_sum (fun j hj => ?_)
  have e : (fun q => a i q * a j q) = fun q => √(a i q ^ 2) * √(a j q ^ 2) := by
    funext q
    rw [Real.sqrt_sq (ha i hi q), Real.sqrt_sq (ha j hj q)]
  rw [e]
  refine (avgQ_sqrt_mul_sqrt_le (fun q => sq_nonneg _) (fun q => sq_nonneg _)).trans ?_
  exact mul_le_mul (hb i hi) (hb j hj) (Real.sqrt_nonneg _)
    ((Real.sqrt_nonneg _).trans (hb i hi))

/-- `|U_{xz}|⁴ ≤ |U_{xz}|` for the unitary propagator. -/
lemma normSq_sq_evo_le_norm (J h : ℝ) (q : Cfg N) (t : ℝ) (x z : Fin N) :
    Complex.normSq (evo (hq J h q) t x z) ^ 2 ≤ ‖evo (hq J h q) t x z‖ := by
  rw [Complex.normSq_eq_norm_sq]
  have h1 := norm_evo_le_one _ (hq_isSymm J h q) t x z
  have h0 := norm_nonneg (evo (hq J h q) t x z)
  have h3 : ‖evo (hq J h q) t x z‖ ^ 3 ≤ 1 := pow_le_one₀ h0 h1
  nlinarith [mul_le_mul_of_nonneg_left h3 h0]

/-- `|K_{xy}| ≤ 1` for the correlation matrix `K = U P_S U†`. -/
lemma norm_corrMat_le_one (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N))
    (t : ℝ) (x y : Fin N) : ‖corrMat A S t x y‖ ≤ 1 := by
  have hrow : ∀ x, ∑ z, ‖evo A t x z‖ ^ 2 = 1 := fun x => by
    simpa [Complex.normSq_eq_norm_sq] using sum_normSq_row (evo A t) (evo_unitary' A hA t) x
  refine (AreaLawAux.norm_corrMat_le A S t x y).trans ?_
  have hz : ∀ z, ‖evo A t x z‖ * ‖evo A t y z‖ ≤ (‖evo A t x z‖ ^ 2 + ‖evo A t y z‖ ^ 2) / 2 :=
    fun z => by nlinarith [sq_nonneg (‖evo A t x z‖ - ‖evo A t y z‖)]
  calc ∑ z, ‖evo A t x z‖ * ‖evo A t y z‖
      ≤ ∑ z, (‖evo A t x z‖ ^ 2 + ‖evo A t y z‖ ^ 2) / 2 := Finset.sum_le_sum fun z _ => hz z
    _ = (∑ z, ‖evo A t x z‖ ^ 2 + ∑ z, ‖evo A t y z‖ ^ 2) / 2 := by
      rw [← Finset.sum_add_distrib, Finset.sum_div]
    _ = 1 := by rw [hrow, hrow]; norm_num

/-- `|K_{xy}|² ≤ ∑_z Q(x,z) Q(y,z)`, uniformly in time. -/
lemma normSq_corrMat_le_corrQ (J h : ℝ) (q : Cfg N) (S : Finset (Fin N)) (t : ℝ) (x y : Fin N) :
    Complex.normSq (corrMat (hq J h q) S t x y) ≤
      ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) x z *
        corrQ (hq J h q) (hq_isHermitian J h q) y z := by
  rw [Complex.normSq_eq_norm_sq]
  have h1 := norm_corrMat_le_one _ (hq_isSymm J h q) S t x y
  have h0 := norm_nonneg (corrMat (hq J h q) S t x y)
  have h2 := AreaLawSharpAux.norm_corrMat_le_corrQ J h q S t x y
  nlinarith

/-- Per-sector second moment of `N_R` about `c`:
`∑_n (N_R(n) - c)² |v_n|² = Var_q(N_R) + (⟨N_R⟩_q - c)²`. -/
lemma sector_second_moment (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm)
    (S R : Finset (Fin N)) (t c : ℝ) :
    ∑ n, (count R n - c) ^ 2 * Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) =
      (∑ x ∈ R, ∑ y ∈ univ \ R, Complex.normSq (corrMat A S t x y)) +
        (∑ x ∈ R, ∑ z ∈ S, Complex.normSq (evo A t x z) - c) ^ 2 := by
  have hnorm := sum_normSq_evo_fock A hA S t
  have hM1 := sum_count_normSq A hA S R t
  have hvar := number_variance A hA S R t
  simp only [numOp, diagonal_mul_diagonal] at hvar
  rw [expect_cplx_diagonal, expect_cplx_diagonal] at hvar
  have hvar' : ∑ n, count R n * count R n *
        Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) -
      (∑ n, count R n * Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n)) ^ 2 =
      ∑ x ∈ R, ∑ y ∈ univ \ R, Complex.normSq (corrMat A S t x y) := by
    exact_mod_cast hvar
  have hexp : ∑ n, (count R n - c) ^ 2 * Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) =
      ∑ n, count R n * count R n * Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) -
        2 * c * ∑ n, count R n * Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) +
        c ^ 2 * ∑ n, Complex.normSq ((evo (dGamma A) t *ᵥ fockVec S) n) := by
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
    exact Finset.sum_congr rfl (fun n _ => by ring)
  rw [hexp, hnorm, hM1]
  rw [hM1] at hvar'
  linear_combination hvar'

/-- The LGT second moment of `N_R` about `c` is the Bernoulli average of the sector second
moments. -/
lemma lgt_second_moment_eq (J h : ℝ) (S R : Finset (Fin N)) (c t : ℝ) :
    (expect (matter ((numOp R - c • (1 : Matrix (Cfg N) (Cfg N) ℝ)) *
        (numOp R - c • (1 : Matrix (Cfg N) (Cfg N) ℝ)))) (lgtState J h (fockVec S) t)).re =
      avgQ (fun q : Cfg N => ∑ n, (count R n - c) ^ 2 *
        Complex.normSq ((evo (dGamma (hq J h q)) t *ᵥ fockVec S) n)) := by
  have e : (numOp R - c • (1 : Matrix (Cfg N) (Cfg N) ℝ)) *
      (numOp R - c • (1 : Matrix (Cfg N) (Cfg N) ℝ)) =
      diagonal (fun n => (count R n - c) ^ 2) := by
    rw [numOp, smul_one_eq_diagonal, diagonal_sub, diagonal_mul_diagonal]
    congr 1
    funext n
    ring
  rw [e, lgt_expect_diag_re]

/-- **Averaged second moment of the transported charge**:
`E_q ∑_n (N_{≤m}(n) - s_m)² |ψ_q(t)(n)|² ≤ Θ(F) + Σ'₁(F)²`. -/
lemma avgQ_second_moment_le {J h : ℝ} {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N))
    (m : ℕ) (t : ℝ) :
    avgQ (fun q : Cfg N => ∑ n, (count (leftOf m) n - ((S ∩ leftOf m).card : ℝ)) ^ 2 *
      Complex.normSq ((evo (dGamma (hq J h q)) t *ᵥ fockVec S) n)) ≤
      thetaF F + sigma1' F ^ 2 := by
  have e1 : (fun q : Cfg N => ∑ n, (count (leftOf m) n - ((S ∩ leftOf m).card : ℝ)) ^ 2 *
      Complex.normSq ((evo (dGamma (hq J h q)) t *ᵥ fockVec S) n)) =
      fun q : Cfg N => (∑ x ∈ leftOf (N := N) m, ∑ y ∈ univ \ leftOf (N := N) m,
        Complex.normSq (corrMat (hq J h q) S t x y)) +
        (∑ x ∈ leftOf (N := N) m, ∑ z ∈ S, Complex.normSq (evo (hq J h q) t x z) -
          ((S ∩ leftOf m).card : ℝ)) ^ 2 := by
    funext q
    exact sector_second_moment _ (hq_isSymm J h q) S _ t _
  rw [e1, avgQ_add]
  refine add_le_add ?_ ?_
  · refine (avgQ_mono fun q => Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ =>
      normSq_corrMat_le_corrQ J h q S t x y).trans ?_
    rw [leftOf_eq_segment]
    exact AreaLawSharpAux.avgQ_cut_corrQ_le hEC (m + 1)
  · have hpt : ∀ q : Cfg N, (∑ x ∈ leftOf (N := N) m, ∑ z ∈ S,
        Complex.normSq (evo (hq J h q) t x z) - ((S ∩ leftOf m).card : ℝ)) ^ 2 ≤
        (∑ p ∈ leftOf (N := N) m ×ˢ (univ \ leftOf (N := N) m),
          Complex.normSq (evo (hq J h q) t p.1 p.2)) ^ 2 := by
      intro q
      have hab := abs_mean_sub_le_symm (evo (hq J h q) t) (evo_unitary _ (hq_isSymm J h q) t)
        (evo_transpose _ (hq_isSymm J h q) t) S (leftOf m)
      have e : ∑ x ∈ leftOf (N := N) m, ∑ z ∈ univ \ leftOf (N := N) m,
          Complex.normSq (evo (hq J h q) t x z) =
          ∑ p ∈ leftOf (N := N) m ×ˢ (univ \ leftOf (N := N) m),
            Complex.normSq (evo (hq J h q) t p.1 p.2) := (Finset.sum_product' _ _ _).symm
      rw [e] at hab
      rw [← sq_abs]
      exact pow_le_pow_left₀ (abs_nonneg _) hab 2
    have hpair : ∀ x z : Fin N,
        √(avgQ (fun q : Cfg N => Complex.normSq (evo (hq J h q) t x z) ^ 2)) ≤
          √(F (Nat.dist x.val z.val)) := fun x z =>
      Real.sqrt_le_sqrt ((avgQ_mono fun q => normSq_sq_evo_le_norm J h q t x z).trans
        (avgQ_norm_evo_le J h hEC x z t))
    refine (avgQ_mono hpt).trans ?_
    refine (avgQ_sq_sum_le (leftOf (N := N) m ×ˢ (univ \ leftOf (N := N) m))
      (fun p q => Complex.normSq (evo (hq J h q) t p.1 p.2))
      (fun p => √(F (Nat.dist p.1.val p.2.val)))
      (fun p _ q => Complex.normSq_nonneg _) (fun p _ => hpair p.1 p.2)).trans ?_
    refine pow_le_pow_left₀ (Finset.sum_nonneg fun p _ => Real.sqrt_nonneg _) ?_ 2
    rw [Finset.sum_product]
    exact sum_cross_le (fun d => √(F d)) (fun d => Real.sqrt_nonneg _)
      (ec_summable_mul_sqrt hEC) m

/-! #### Number distributions -/

/-- The particle number of `R` in the configuration `n`, as a natural number. -/
def countNat (R : Finset (Fin N)) (n : Cfg N) : ℕ := (R.filter (fun x => n x ≠ 0)).card

/-- `count` is the natural number `countNat`. -/
lemma count_eq_countNat (R : Finset (Fin N)) (n : Cfg N) : count R n = countNat R n := by
  unfold count countNat occ
  rw [Finset.card_filter]
  push_cast
  refine Finset.sum_congr rfl (fun x _ => ?_)
  split_ifs <;> simp_all

/-- The particle number is at most `N`. -/
lemma countNat_le (R : Finset (Fin N)) (n : Cfg N) : countNat R n ≤ N :=
  (Finset.card_filter_le _ _).trans ((Finset.card_le_univ R).trans (by simp))

/-- Summing over the level sets `{count R = k}`, `k ∈ {0, …, N}`, is summing over everything. -/
lemma sum_fiber_count (R : Finset (Fin N)) (φ : Cfg N → ℝ) :
    ∑ k ∈ range (N + 1), ∑ n ∈ univ.filter (fun n => count R n = k), φ n = ∑ n, φ n := by
  have h : ∀ k : ℕ, univ.filter (fun n : Cfg N => count R n = k) =
      univ.filter (fun n => countNat R n = k) := by
    intro k; ext n; simp [count_eq_countNat]
  simp_rw [h]
  exact Finset.sum_fiberwise_of_maps_to
    (fun n _ => mem_range.2 (Nat.lt_succ_of_le (countNat_le R n))) φ

/-- The number distribution of a normalized vector is a probability distribution. -/
lemma numDist_isProb (R : Finset (Fin N)) (v : Cfg N → ℂ)
    (hv : ∑ n, Complex.normSq (v n) = 1) :
    IsProb N (numDist R v) where
  nonneg k := Finset.sum_nonneg (fun n _ => Complex.normSq_nonneg _)
  sum_one := by unfold numDist; rw [sum_fiber_count]; exact hv

/-- Second moment about `c` of the number distribution. -/
lemma sum_numDist_mul_sq (R : Finset (Fin N)) (v : Cfg N → ℂ) (c : ℝ) :
    ∑ k ∈ range (N + 1), numDist R v k * ((k : ℝ) - c) ^ 2 =
      ∑ n, (count R n - c) ^ 2 * Complex.normSq (v n) := by
  unfold numDist
  simp_rw [Finset.sum_mul]
  rw [← sum_fiber_count R]
  refine Finset.sum_congr rfl (fun k _ => Finset.sum_congr rfl (fun n hn => ?_))
  rw [(mem_filter.1 hn).2]
  ring

/-- The LGT number distribution is the Bernoulli average of the sector number distributions. -/
lemma lgtNumDist_eq_avgQ (J h : ℝ) (R : Finset (Fin N)) (ψ : Cfg N → ℂ) (t : ℝ) (k : ℕ) :
    lgtNumDist J h R ψ t k = avgQ (fun q => numDist R (evo (dGamma (hq J h q)) t *ᵥ ψ) k) := by
  rw [lgtNumDist, numProj, lgt_expect_diag_re]
  congr 1
  funext q
  rw [numDist, Finset.sum_filter]
  refine Finset.sum_congr rfl (fun n _ => ?_)
  split_ifs <;> ring

/-- **Jensen's inequality** for the logarithm and the Bernoulli average. -/
lemma avgQ_log_le {f : Cfg N → ℝ} (hf : ∀ q, 0 < f q) :
    avgQ (fun q => Real.log (f q)) ≤ Real.log (avgQ f) := by
  have hw : ∑ _q : Cfg N, (1 / 2 ^ N : ℝ) = 1 := by
    rw [Finset.sum_const, Finset.card_univ, card_cfg, nsmul_eq_mul]
    push_cast
    field_simp
  have h := strictConcaveOn_log_Ioi.concaveOn.le_map_sum (t := (univ : Finset (Cfg N)))
    (w := fun _ => (1 / 2 ^ N : ℝ)) (p := f) (fun _ _ => by positivity) hw
    (fun q _ => Set.mem_Ioi.mpr (hf q))
  simp only [smul_eq_mul] at h
  simpa only [avgQ, Finset.mul_sum] using h

/-! #### Dressed hopping -/

/-- `c̃_x† c̃_y` is gauge invariant. -/
lemma gaugeInv_dressed_hop (x y : Fin N) : GaugeInv (cplx ((dressed x)ᵀ * dressed y)) := by
  intro j
  rw [← cplx_mul', ← cplx_mul']
  congr 1
  rw [Matrix.mul_assoc, dressed_comm_gauss, ← Matrix.mul_assoc, dressed_transpose_comm_gauss,
    Matrix.mul_assoc]

/-- In the dual frame `U† c̃_x† c̃_y U = c_x† c_y ⊗ 1`, so every sector block of `c̃_x† c̃_y` is the
free hopping `c_x† c_y`. -/
lemma sectorBlock_dressed_hop (x y : Fin N) (q : Cfg N) :
    sectorBlock (cplx ((dressed x)ᵀ * dressed y)) q = cplx (cre x * ann y) := by
  have hk : (ann x ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))ᵀ * (ann y ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) =
      (cre x * ann y) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ) := by
    rw [← kroneckerMap_transpose, transpose_one, ← mul_kronecker_mul, Matrix.one_mul]
    rfl
  have hd : (dressed x)ᵀ * dressed y =
      ((cre x * ann y) ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)).submatrix (dualEquiv (N := N)).symm
        dualEquiv.symm := by
    rw [dressed_eq x, dressed_eq y, transpose_submatrix, submatrix_mul_equiv, hk]
  ext n' n
  simp only [sectorBlock, hd, cplx, Matrix.map_apply, submatrix_apply, Equiv.symm_apply_apply,
    kronecker_apply, Matrix.one_apply_eq, mul_one]

end TransportSharpAux

open TransportSharpAux

/-! ### (a) The dressed propagator -/

/-- **The dressed anticommutator is the one-body propagator in each sector** (paper: eq. `eq:main-anticomm`). -/
theorem dressed_anticomm_t (J h : ℝ) (x y : Fin N) (t : ℝ) :
    dressedAnticomm J h x y t = ∑ q : Cfg N, evo (hq J h q) t x y • cplx (gaussProj q) := by
  unfold dressedAnticomm dressedT
  rw [heisenberg_dressed]
  calc (∑ q : Cfg N, cplx (gaussProj q) * ∑ z : Fin N, evo (hq J h q) t x z • cplx (dressed z)) *
        cplx (dressed y)ᵀ + cplx (dressed y)ᵀ *
        (∑ q : Cfg N, cplx (gaussProj q) * ∑ z : Fin N, evo (hq J h q) t x z • cplx (dressed z))
      = ∑ q : Cfg N, cplx (gaussProj q) * ∑ z : Fin N, evo (hq J h q) t x z •
          (cplx (dressed z) * cplx (dressed y)ᵀ + cplx (dressed y)ᵀ * cplx (dressed z)) := by
        rw [Finset.sum_mul, Finset.mul_sum, ← Finset.sum_add_distrib]
        refine Finset.sum_congr rfl fun q _ => ?_
        rw [Matrix.mul_assoc, ← Matrix.mul_assoc (cplx (dressed y)ᵀ) (cplx (gaussProj q)),
          ← cplx_gaussProj_comm_dressed_transpose q y, Matrix.mul_assoc, ← Matrix.mul_add]
        congr 1
        rw [Finset.sum_mul, Finset.mul_sum, ← Finset.sum_add_distrib]
        refine Finset.sum_congr rfl fun z _ => ?_
        rw [Matrix.smul_mul, Matrix.mul_smul, ← smul_add]
    _ = ∑ q : Cfg N, evo (hq J h q) t x y • cplx (gaussProj q) := by
        refine Finset.sum_congr rfl fun q _ => ?_
        simp only [cplx_dressed_CAR, smul_ite, smul_zero, Finset.sum_ite_eq', Finset.mem_univ,
          ite_true]
        rw [Matrix.mul_smul, Matrix.mul_one]

/-- `‖P_q {c̃_x(t), c̃_y†} P_q‖ = |U_q(t)_{xy}|`. -/
theorem norm_sector_anticomm (J h : ℝ) (x y : Fin N) (t : ℝ) (q : Cfg N) :
    ‖cplx (gaussProj q) * dressedAnticomm J h x y t * cplx (gaussProj q)‖ = ‖evo (hq J h q) t x y‖ := by
  have hdiag : cplx (gaussProj q) * dressedAnticomm J h x y t * cplx (gaussProj q) =
      diagonal (fun p : Cfg N × Cfg N =>
        if chargeOf p.1 p.2 = q then evo (hq J h q) t x y else 0) := by
    rw [dressed_anticomm_t, SectorNormAux.sum_smul_gaussProj, SectorNormAux.cplx_gaussProj,
      diagonal_mul_diagonal, diagonal_mul_diagonal]
    congr 1
    funext p
    split_ifs with hp
    · rw [hp]; ring
    · ring
  rw [hdiag, l2_opNorm_diagonal]
  apply le_antisymm
  · refine (pi_norm_le_iff_of_nonneg (norm_nonneg _)).mpr fun p => ?_
    split_ifs
    · exact le_rfl
    · simp
  · set v : Cfg N × Cfg N → ℂ := fun p =>
      if chargeOf p.1 p.2 = q then evo (hq J h q) t x y else 0 with hv
    have h1 : v ((0 : Cfg N), linkOf 0 q) = evo (hq J h q) t x y := by
      simp only [hv]
      exact ite_eq_left (chargeOf_linkOf 0 q)
    calc ‖evo (hq J h q) t x y‖ = ‖v ((0 : Cfg N), linkOf 0 q)‖ := by rw [h1]
      _ ≤ ‖v‖ := norm_le_pi_norm v _

/-- **Localization of the dressed propagator** (paper: eq. `eq:main-prop`). -/
theorem avg_sup_anticomm_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (x y : Fin N) :
    avgQ (fun q : Cfg N =>
      ⨆ t : ℝ, ‖cplx (gaussProj q) * dressedAnticomm J h x y t * cplx (gaussProj q)‖) ≤
      F (Nat.dist x.val y.val) := by
  simp_rw [norm_sector_anticomm]
  exact avgQ_iSup_norm_evo_le J h hEC x y

/-- **`⟨Ψ₀, |{c̃_x(t), c̃_y†}| Ψ₀⟩ ≤ F(|x-y|)`** (paper: Theorem `thm:transport`(a)): the operator
`∑_q |U_q(t)_{xy}| P_q` is positive semidefinite with square `{c̃_x(t), c̃_y†}^* {c̃_x(t), c̃_y†}` (so it is the
absolute value), and its expectation in `Ψ₀ = |↑⋯↑⟩_z ⊗ ψ` is at most `F(|x-y|)`. -/
theorem expect_abs_anticomm_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (ψ : Cfg N → ℂ)
    (hψ : star ψ ⬝ᵥ ψ = 1) (x y : Fin N) (t : ℝ) :
    (∑ q : Cfg N, ((‖evo (hq J h q) t x y‖ : ℝ) : ℂ) • cplx (gaussProj q)).PosSemidef ∧
      (∑ q : Cfg N, ((‖evo (hq J h q) t x y‖ : ℝ) : ℂ) • cplx (gaussProj q)) *
          (∑ q : Cfg N, ((‖evo (hq J h q) t x y‖ : ℝ) : ℂ) • cplx (gaussProj q)) =
        (dressedAnticomm J h x y t)ᴴ * dressedAnticomm J h x y t ∧
      (expect (∑ q : Cfg N, ((‖evo (hq J h q) t x y‖ : ℝ) : ℂ) • cplx (gaussProj q)) (initZ ψ)).re ≤
        F (Nat.dist x.val y.val) := by
  obtain ⟨h1, h2⟩ := sectorDiag_abs (fun q : Cfg N => evo (hq J h q) t x y)
  refine ⟨h1, ?_, ?_⟩
  · rw [h2, dressed_anticomm_t]
  · rw [expect_sectorDiag_abs ψ hψ (fun q : Cfg N => evo (hq J h q) t x y), Complex.ofReal_re]
    exact avgQ_norm_evo_le J h hEC x y t

/-! ### (b) One particle -/

/-- **Moments of a single particle** (paper: Theorem `thm:transport`(b)), `0 ≤ p ≤ 2`. -/
theorem one_particle_moments (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (y : Fin N) (p : ℝ)
    (hp0 : 0 ≤ p) (hp2 : p ≤ 2) (t : ℝ) :
    ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p *
        (expect (matter (diagonal fun n => occ (n x))) (lgtState J h (fockVec {y}) t)).re ≤
      avgQ (fun q : Cfg N => ⨆ s : ℝ,
        ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p * ‖evo (hq J h q) s x y‖ ^ 2) ∧
    avgQ (fun q : Cfg N => ⨆ s : ℝ,
        ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p * ‖evo (hq J h q) s x y‖ ^ 2) ≤
      2 * momentF p F := by
  -- the lower bound `0 ≤ p` is not needed: `(d + 1)^p ≥ 0` for every real `p`
  have _ := hp0
  have hw0 : ∀ x : Fin N, 0 ≤ ((Nat.dist x.val y.val : ℝ) + 1) ^ p := fun x =>
    Real.rpow_nonneg (by positivity) p
  constructor
  · have e : ∀ x : Fin N,
        (expect (matter (diagonal fun n => occ (n x))) (lgtState J h (fockVec {y}) t)).re =
          avgQ (fun q : Cfg N => ‖evo (hq J h q) t x y‖ ^ 2) := by
      intro x
      rw [lgt_density]
      simp only [Finset.sum_singleton, Complex.normSq_eq_norm_sq]
    simp_rw [e, ← avgQ_const_mul]
    rw [← avgQ_sum]
    refine avgQ_mono fun q => ?_
    refine le_ciSup (f := fun s : ℝ => ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p *
      ‖evo (hq J h q) s x y‖ ^ 2) ⟨∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p, ?_⟩ t
    rintro _ ⟨s, rfl⟩
    refine Finset.sum_le_sum fun x _ => mul_le_of_le_one_right (hw0 x) ?_
    exact pow_le_one₀ (norm_nonneg _) (norm_evo_le_one _ (hq_isSymm J h q) s x y)
  · have hpt : ∀ q : Cfg N, (⨆ s : ℝ, ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p *
        ‖evo (hq J h q) s x y‖ ^ 2) ≤ ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p *
        corrQ (hq J h q) (hq_isHermitian J h q) x y := fun q =>
      ciSup_le fun s => Finset.sum_le_sum fun x _ =>
        mul_le_mul_of_nonneg_left (norm_evo_sq_le_corrQ J h q s x y) (hw0 x)
    have hg0 : ∀ d : ℕ, 0 ≤ ((d : ℝ) + 1) ^ p * F d := fun d =>
      mul_nonneg (Real.rpow_nonneg (by positivity) p) (hEC.nonneg d)
    have hgs : Summable (fun d : ℕ => ((d : ℝ) + 1) ^ p * F d) := by
      refine hEC.summable.of_nonneg_of_le hg0 fun d => ?_
      refine mul_le_mul_of_nonneg_right ?_ (hEC.nonneg d)
      have h1 : (1 : ℝ) ≤ (d : ℝ) + 1 := by linarith [(Nat.cast_nonneg d : (0 : ℝ) ≤ d)]
      calc ((d : ℝ) + 1) ^ p ≤ ((d : ℝ) + 1) ^ (2 : ℝ) := Real.rpow_le_rpow_of_exponent_le h1 hp2
        _ = ((d : ℝ) + 1) ^ 2 := Real.rpow_two _
    calc avgQ (fun q : Cfg N => ⨆ s : ℝ, ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p *
          ‖evo (hq J h q) s x y‖ ^ 2)
        ≤ avgQ (fun q : Cfg N => ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p *
          corrQ (hq J h q) (hq_isHermitian J h q) x y) := avgQ_mono hpt
      _ = ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p *
          avgQ (fun q : Cfg N => corrQ (hq J h q) (hq_isHermitian J h q) x y) := by
          rw [avgQ_sum]
          exact Finset.sum_congr rfl fun x _ => avgQ_const_mul _ _
      _ ≤ ∑ x : Fin N, ((Nat.dist x.val y.val : ℝ) + 1) ^ p * F (Nat.dist x.val y.val) :=
          Finset.sum_le_sum fun x _ => mul_le_mul_of_nonneg_left (hEC.bound N x y) (hw0 x)
      _ ≤ 2 * momentF p F :=
          AreaLawAux.sum_dist_le (fun d : ℕ => ((d : ℝ) + 1) ^ p * F d) hg0 hgs y

/-! ### (c) Transport across a cut -/

/-- **Transport in one sector** (paper: Theorem `thm:transport`(c)). -/
theorem sector_transport_le (J h : ℝ) (q : Cfg N) (S : Finset (Fin N)) (m : ℕ) (t : ℝ) :
    |(expect (cplx (numOp (leftOf m))) (evo (dGamma (hq J h q)) t *ᵥ fockVec S)).re -
        ((S ∩ leftOf m).card : ℝ)| ≤
      ∑ x ∈ leftOf m, ∑ z ∈ univ \ leftOf m, ‖evo (hq J h q) t x z‖ ^ 2 := by
  rw [expect_numOp_re _ (hq_isSymm J h q)]
  refine (abs_mean_sub_le_symm _ (evo_unitary _ (hq_isSymm J h q) t)
    (evo_transpose _ (hq_isSymm J h q) t) S _).trans (le_of_eq ?_)
  simp only [Complex.normSq_eq_norm_sq]

/-- `E_q sup_t |⟨N_{≤m}⟩_{q,t} - s_m| ≤ Σ₁(F)` (paper: eq. `eq:main-transport`). -/
theorem avg_sup_transport_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ,
      |(expect (cplx (numOp (leftOf m))) (evo (dGamma (hq J h q)) t *ᵥ fockVec S)).re -
        ((S ∩ leftOf m).card : ℝ)|) ≤ sigma1 F := by
  refine (avgQ_mono fun q => ciSup_le fun t => ?_).trans (avgQ_cross_corrQ_le hEC m)
  exact (sector_transport_le J h q S m t).trans (Finset.sum_le_sum fun x _ =>
    Finset.sum_le_sum fun z _ => norm_evo_sq_le_corrQ J h q t x z)

/-- **Bounded transport with the constant of the paper** (paper: eq. `eq:main-transport`):
`|⟨Ψ(t), N_{≤m} Ψ(t)⟩ - s_m| ≤ Σ₁(F)`. -/
theorem lgt_transport_sharp (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ)
    (t : ℝ) :
    |(expect (matter (numOp (leftOf m))) (lgtState J h (fockVec S) t)).re -
        ((S ∩ leftOf m).card : ℝ)| ≤ sigma1 F := by
  rw [numOp, lgt_expect_diag_re, ← avgQ_sub_const]
  refine (abs_avgQ_le _).trans ((avgQ_mono fun q => ?_).trans (avgQ_cross_corrQ_le hEC m))
  have h1 := sector_transport_le J h q S m t
  rw [numOp, expect_cplx_diagonal, Complex.ofReal_re] at h1
  exact h1.trans (Finset.sum_le_sum fun x _ =>
    Finset.sum_le_sum fun z _ => norm_evo_sq_le_corrQ J h q t x z)

/-- **Second moment of the transported charge** (paper: eq. `eq:main-transport2`):
`⟨Ψ(t), (N_{≤m} - s_m)² Ψ(t)⟩ ≤ Θ(F) + Σ'₁(F)²`. -/
theorem lgt_transport_sq (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ)
    (t : ℝ) :
    (expect (matter ((numOp (leftOf m) - ((S ∩ leftOf m).card : ℝ) • 1) *
        (numOp (leftOf m) - ((S ∩ leftOf m).card : ℝ) • 1))) (lgtState J h (fockVec S) t)).re ≤
      thetaF F + sigma1' F ^ 2 := by
  rw [lgt_second_moment_eq]
  exact avgQ_second_moment_le hEC S m t

/-! ### (d) Dressed correlations -/

/-- **Dressed correlations** (paper: Theorem `thm:transport`(d)): for a Fock state and `x ≠ y`,
`|⟨Ψ(t), c̃_x† c̃_y Ψ(t)⟩| ≤ 2 ∑_{k ≥ ⌈|x-y|/2⌉} F(k)`. -/
theorem dressed_corr_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (x y : Fin N)
    (hxy : x ≠ y) (t : ℝ) :
    ‖expect (cplx ((dressed x)ᵀ * dressed y)) (lgtState J h (fockVec S) t)‖ ≤
      2 * ∑' k : ℕ, (if Nat.dist x.val y.val ≤ 2 * k then F k else 0) := by
  -- the bound also holds for `x = y`
  have _ := hxy
  rw [sector_average J h (fockVec S) _ (gaugeInv_dressed_hop x y) t]
  have e : ∀ q : Cfg N, expect (sectorBlock (cplx ((dressed x)ᵀ * dressed y)) q)
      (evo (Hsec J h q) t *ᵥ fockVec S) = corrMat (hq J h q) S t y x := by
    intro q
    rw [sectorBlock_dressed_hop, Hsec, expect_evo_sub_smul_one, two_point _ (hq_isSymm J h q)]
  simp_rw [e]
  calc ‖(1 / 2 ^ N : ℂ) * ∑ q : Cfg N, corrMat (hq J h q) S t y x‖
      ≤ avgQ (fun q : Cfg N => ‖corrMat (hq J h q) S t y x‖) := by
        rw [norm_mul, avgQ]
        have h2 : ‖(1 / 2 ^ N : ℂ)‖ = 1 / 2 ^ N := by
          rw [norm_div, norm_one, norm_pow, Complex.norm_ofNat]
        rw [h2]
        exact mul_le_mul_of_nonneg_left (norm_sum_le _ _) (by positivity)
    _ ≤ avgQ (fun q : Cfg N => ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) y z *
          corrQ (hq J h q) (hq_isHermitian J h q) x z) :=
        avgQ_mono fun q => AreaLawSharpAux.norm_corrMat_le_corrQ J h q S t y x
    _ ≤ 2 * ∑' k : ℕ, (if Nat.dist y.val x.val ≤ 2 * k then F k else 0) :=
        AreaLawSharpAux.avgQ_sum_corrQ_le hEC y x
    _ = 2 * ∑' k : ℕ, (if Nat.dist x.val y.val ≤ 2 * k then F k else 0) := by
        rw [Nat.dist_comm]

/-! ### Number entropy -/

/-- **Bounded number entropy with the constants of the paper** (paper: Theorem `thm:numberentropy`):
`S_N(t) ≤ ½ ln(2πe(K₁ + K₂ + 1/12))`, `K₁ = Θ(F)`, `K₂ = Σ'₁(F)²`. -/
theorem lgt_numberEntropy_sharp (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N)) (m : ℕ)
    (t : ℝ) :
    shannon N (lgtNumDist J h (leftOf m) (fockVec S) t) ≤
      (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (thetaF F + sigma1' F ^ 2 + 1 / 12)) := by
  have hp := lgtNumDist_isProb J h (leftOf m) S t
  refine shannon_le_of_pvar_le hp ?_
  calc pvar N (lgtNumDist J h (leftOf m) (fockVec S) t)
      ≤ ∑ k ∈ range (N + 1), lgtNumDist J h (leftOf m) (fockVec S) t k *
          ((k : ℝ) - ((S ∩ leftOf m).card : ℝ)) ^ 2 := pvar_le_sum_sq hp _
    _ = avgQ (fun q => ∑ k ∈ range (N + 1),
          numDist (leftOf m) (evo (dGamma (hq J h q)) t *ᵥ fockVec S) k *
            ((k : ℝ) - ((S ∩ leftOf m).card : ℝ)) ^ 2) := by
        rw [avgQ_sum]
        refine Finset.sum_congr rfl (fun k _ => ?_)
        rw [lgtNumDist_eq_avgQ, mul_comm, ← avgQ_const_mul]
        congr 1
        funext q
        ring
    _ = avgQ (fun q : Cfg N => ∑ n, (count (leftOf m) n - ((S ∩ leftOf m).card : ℝ)) ^ 2 *
          Complex.normSq ((evo (dGamma (hq J h q)) t *ᵥ fockVec S) n)) := by
        congr 1
        funext q
        exact sum_numDist_mul_sq _ _ _
    _ ≤ thetaF F + sigma1' F ^ 2 := avgQ_second_moment_le hEC S m t

/-- **Number entropy in a single sector** (paper: Theorem `thm:numberentropy`):
`E_q sup_t H(p^q_t) ≤ ½ ln(2πe(Θ(F) + 1/12))`. -/
theorem avg_sup_sector_numberEntropy_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (S : Finset (Fin N))
    (m : ℕ) :
    avgQ (fun q : Cfg N => ⨆ t : ℝ,
      shannon N (numDist (leftOf m) (evo (dGamma (hq J h q)) t *ᵥ fockVec S))) ≤
      (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (thetaF F + 1 / 12)) := by
  set V : Cfg N → ℝ := fun q => ∑ x ∈ leftOf (N := N) m, ∑ y ∈ univ \ leftOf (N := N) m,
    ∑ z, corrQ (hq J h q) (hq_isHermitian J h q) x z *
      corrQ (hq J h q) (hq_isHermitian J h q) y z with hV
  have hV0 : ∀ q, 0 ≤ V q := fun q =>
    Finset.sum_nonneg fun x _ => Finset.sum_nonneg fun y _ => Finset.sum_nonneg fun z _ =>
      mul_nonneg (corrQ_nonneg _ _ x z) (corrQ_nonneg _ _ y z)
  have hVavg : avgQ V ≤ thetaF F := by
    rw [hV, leftOf_eq_segment]
    exact AreaLawSharpAux.avgQ_cut_corrQ_le hEC (m + 1)
  have hpt : ∀ q : Cfg N, (⨆ t : ℝ,
      shannon N (numDist (leftOf m) (evo (dGamma (hq J h q)) t *ᵥ fockVec S))) ≤
      (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (V q + 1 / 12)) := by
    intro q
    refine ciSup_le fun t => ?_
    have hp := numDist_isProb (leftOf m) _ (sum_normSq_evo_fock _ (hq_isSymm J h q) S t)
    refine shannon_le_of_pvar_le hp ?_
    set μ : ℝ := ∑ x ∈ leftOf (N := N) m, ∑ z ∈ S, Complex.normSq (evo (hq J h q) t x z)
    calc pvar N (numDist (leftOf m) (evo (dGamma (hq J h q)) t *ᵥ fockVec S))
        ≤ ∑ k ∈ range (N + 1), numDist (leftOf m) (evo (dGamma (hq J h q)) t *ᵥ fockVec S) k *
            ((k : ℝ) - μ) ^ 2 := pvar_le_sum_sq hp μ
      _ = (∑ x ∈ leftOf (N := N) m, ∑ y ∈ univ \ leftOf (N := N) m,
            Complex.normSq (corrMat (hq J h q) S t x y)) + (μ - μ) ^ 2 := by
          rw [sum_numDist_mul_sq]
          exact sector_second_moment _ (hq_isSymm J h q) S _ t μ
      _ ≤ V q := by
          rw [sub_self, sq, mul_zero, add_zero]
          exact Finset.sum_le_sum fun x _ => Finset.sum_le_sum fun y _ =>
            normSq_corrMat_le_corrQ J h q S t x y
  have hpos : ∀ q, 0 < 2 * Real.pi * Real.exp 1 * (V q + 1 / 12) := fun q => by
    have := hV0 q
    positivity
  have hlin : avgQ (fun q => 2 * Real.pi * Real.exp 1 * (V q + 1 / 12)) =
      2 * Real.pi * Real.exp 1 * (avgQ V + 1 / 12) := by
    rw [avgQ_const_mul, avgQ_add, avgQ_const]
  have hA0 : 0 ≤ avgQ V := by
    have := avgQ_mono (f := fun _ => (0 : ℝ)) (g := V) hV0
    rwa [avgQ_const] at this
  calc avgQ (fun q : Cfg N => ⨆ t : ℝ,
        shannon N (numDist (leftOf m) (evo (dGamma (hq J h q)) t *ᵥ fockVec S)))
      ≤ avgQ (fun q => (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (V q + 1 / 12))) :=
        avgQ_mono hpt
    _ = (1 / 2) * avgQ (fun q => Real.log (2 * Real.pi * Real.exp 1 * (V q + 1 / 12))) :=
        avgQ_const_mul _ _
    _ ≤ (1 / 2) * Real.log (avgQ (fun q => 2 * Real.pi * Real.exp 1 * (V q + 1 / 12))) :=
        mul_le_mul_of_nonneg_left (avgQ_log_le hpos) (by norm_num)
    _ = (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (avgQ V + 1 / 12)) := by rw [hlin]
    _ ≤ (1 / 2) * Real.log (2 * Real.pi * Real.exp 1 * (thetaF F + 1 / 12)) := by
        refine mul_le_mul_of_nonneg_left (Real.log_le_log (by positivity) ?_) (by norm_num)
        have : 0 < 2 * Real.pi * Real.exp 1 := by positivity
        nlinarith

end DFL
