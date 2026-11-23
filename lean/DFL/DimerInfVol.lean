import DFL.DimerLGT
import DFL.InfVol

/-!
# The thermodynamic limit of the dimer dynamics (paper: Lemma `lem:dimerTL`(ii) for `Λ' = ℤ`,
Proposition `prop:dimerlower` for `I = ℤ`, Theorem `thm:dimer`(b) and the lower bound in (c))

On `ℤ`, `H_w = -Δ + diag(w)` is an infinite tridiagonal matrix. Its powers are finite sums (`schrZpow`), and
the matrix entries of `e^{-isH_w}` are given by the absolutely convergent exponential series (`evoZw`), as for
`UZ` in `InfVol.lean`. No operator theory on `ℓ²(ℤ)` is needed.

* `dimerTL_ii_Z` (paper: `lem:dimerTL`(ii) with `Λ' = ℤ`):
  `|⟨δ_y, (e^{-isH_Λ} - e^{-isH_ℤ}) δ_0⟩| ≤ 2 ∑_{k ≥ R} ((2+W)|s|)^k / k!`.
* `moment_poly_bound`: `∑_y |y|^p |⟨δ_y, e^{-isH_Λ} δ_0⟩|² ≤ (4e²|s|)^p + ∑_{y ∈ ℤ} |y|^p 4^{-|y|}` for every
  `p > 0`, uniformly in the interval `Λ` and the potential.
* `tendsto_moment_Z`, `tendsto_dimerMom_Z`: the moments and the time-averaged moments converge as `Λ ↑ ℤ`.
* `dimer_lower_bound_Z` (paper: `prop:dimerlower` for `I = ℤ`): `M^ℤ_{v,p}(T; x) ≥ c_{λ,p} T^{p-1}` for every
  dimer configuration `v` on `ℤ`, every `x ∈ ℤ`, `p > 0` and `T ≥ 1`.
* `dimerMN_lower` (paper: `thm:dimer`(a), averaged form): `M^{(N)}_p(T; x) ≥ c_{λ,p} (JT)^{p-1}`.
* `dimer_TL_D`, `dimer_TL_M` (paper: `thm:dimer`(b)): the sector-averaged moments `D^{(N)}_p(t; x_N)` and their
  time averages `M^{(N)}_p(T; x_N)` converge to those of the random dimer model `v^{(r)}` on `ℤ`.
* `dimerMInf_lower` (paper: `thm:dimer`(c), first statement): `M^{∞,r}_p(T) ≥ c_{λ,p} (JT)^{p-1}` for `JT ≥ 1`.
* `dimerD_unbounded` (paper: `thm:dimer`(d), eq. `eq:dimerLGTd`, for `p > 1`).

The dimer law: for even `N` and a site `x ≡ r (mod 2)` (0-indexed; the paper's `x_N ≡ r + 1` with 1-indexed
sites), translating the chain `[0, N-1]` by `-x` maps the dimers `{2k, 2k+1}` onto the blocks
`{2k' - r, 2k' - r + 1}` of `v^{(r)}`, and the uniform law on the dimer sectors becomes the law of
`v^{(r)}|_{[-x, N-1-x]}` under `muZ`.
-/

open MeasureTheory Filter Topology Matrix Finset

noncomputable section

namespace DFL

/-- `H_w = -Δ_ℤ + diag(w)` as an infinite tridiagonal matrix (same convention as `schrH`). -/
def schrZ (w : ℤ → ℝ) (x y : ℤ) : ℝ := if x = y then w x else if |x - y| = 1 then -1 else 0

/-- The powers of `H_w` (finite sums: the matrix is tridiagonal). -/
def schrZpow (w : ℤ → ℝ) : ℕ → ℤ → ℤ → ℝ
  | 0 => fun x y => if x = y then 1 else 0
  | (k + 1) => fun x y => ∑ z ∈ Finset.Icc (x - 1) (x + 1), schrZ w x z * schrZpow w k z y

/-- The matrix entries of `e^{-isH_w}` on `ℓ²(ℤ)`: `∑_k (-is)^k/k! (H_w^k)_{xy}`. -/
def evoZw (w : ℤ → ℝ) (s : ℝ) (x y : ℤ) : ℂ :=
  ∑' k : ℕ, ((-(s : ℂ) * Complex.I) ^ k / (k.factorial : ℂ)) * ((schrZpow w k x y : ℝ) : ℂ)

/-- The time-averaged moment on `ℤ` (paper: eq. `eq:Mom` with `I = ℤ`):
`M^ℤ_{w,p}(T; x) = ∫_0^∞ (dt/T) e^{-t/T} ∑_y |y - x|^p |⟨δ_y, e^{-itH_w} δ_x⟩|²`. -/
def dimerMomZ (w : ℤ → ℝ) (p T : ℝ) (x : ℤ) : ℝ :=
  ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) *
    ∑' y : ℤ, |((y - x : ℤ) : ℝ)| ^ p * ‖evoZw w t y x‖ ^ 2

/-- Dimer configurations on `ℤ` with values `±λ` and alignment `r`: `v(2k - r) = v(2k - r + 1)`. -/
def IsDimerCfgZ (lam : ℝ) (r : ℤ) (v : ℤ → ℝ) : Prop :=
  (∀ y, v y = lam ∨ v y = -lam) ∧ ∀ k : ℤ, v (2 * k - r) = v (2 * k - r + 1)

/-- The random dimer potential `v^{(r)}` on `ℤ` built from the signs `η` (i.i.d. uniform under `muZ`):
`v(y) = (2h/J) sgn η_k` for `y ∈ {2k - r, 2k - r + 1}`. -/
def dimerPotZ (J h : ℝ) (r : ℤ) (η : CfgZ) (y : ℤ) : ℝ := 2 * h / J * sgn (η ((y + r) / 2))

/-- The sector-averaged moment `D^{(N)}_p(t; x) = E_dim ∑_y |y - x|^p |⟨δ_y, e^{-ith_q} δ_x⟩|²` (by eq.
`eq:dimerdensity` this is the moment of the single-fermion dimer state; see `DimerState.lean`). -/
def dimerD (J h : ℝ) (N : ℕ) (p t : ℝ) (x : Fin N) : ℝ :=
  dimerAvg (fun q : Cfg N =>
    ∑ y : Fin N, |((y.val : ℝ) - x.val)| ^ p * ‖evo (hq J h q) t y x‖ ^ 2)

/-- `M^{(N)}_p(T; x) = ∫_0^∞ (dt/T) e^{-t/T} D^{(N)}_p(t; x)` (paper: eq. `eq:LGTmoments`). -/
def dimerMN (J h : ℝ) (N : ℕ) (p T : ℝ) (x : Fin N) : ℝ :=
  ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) * dimerD J h N p t x

/-- `M^{∞,r}_p(T) = E M^ℤ_{v^{(r)},p}(JT; 0)`. -/
def dimerMInf (J h : ℝ) (r : ℤ) (p T : ℝ) : ℝ :=
  ∫ η, dimerMomZ (dimerPotZ J h r η) p (J * T) 0 ∂muZ

/-! ### The infinite-volume propagator -/

namespace DimerInfVolAux

/-- The three neighbours `x - 1, x, x + 1`. -/
lemma Icc_three (x : ℤ) : Finset.Icc (x - 1) (x + 1) = {x - 1, x, x + 1} := by
  ext z
  simp only [Finset.mem_Icc, Finset.mem_insert, Finset.mem_singleton]
  omega

/-- A sum over the three neighbours. -/
lemma sum_Icc_three {M : Type*} [AddCommMonoid M] (x : ℤ) (f : ℤ → M) :
    ∑ z ∈ Finset.Icc (x - 1) (x + 1), f z = f (x - 1) + f x + f (x + 1) := by
  rw [Icc_three, Finset.sum_insert (by simp only [Finset.mem_insert, Finset.mem_singleton]; omega),
    Finset.sum_insert (by simp only [Finset.mem_singleton]; omega), Finset.sum_singleton,
    add_assoc]

/-- `H_w(x, x - 1) = -1`. -/
lemma schrZ_left (w : ℤ → ℝ) (x : ℤ) : schrZ w x (x - 1) = -1 := by
  unfold schrZ
  rw [ite_eq_right (by omega), ite_eq_left (by rw [sub_sub_cancel, abs_one])]

/-- `H_w(x, x) = w(x)`. -/
lemma schrZ_diag (w : ℤ → ℝ) (x : ℤ) : schrZ w x x = w x := by
  unfold schrZ
  rw [ite_eq_left rfl]

/-- `H_w(x, x + 1) = -1`. -/
lemma schrZ_right (w : ℤ → ℝ) (x : ℤ) : schrZ w x (x + 1) = -1 := by
  unfold schrZ
  rw [ite_eq_right (by omega), ite_eq_left (by rw [sub_add_cancel_left, abs_neg, abs_one])]

/-- `H_w` is tridiagonal. -/
lemma schrZ_eq_zero (w : ℤ → ℝ) (x z : ℤ) (hz : z ∉ Finset.Icc (x - 1) (x + 1)) :
    schrZ w x z = 0 := by
  simp only [Finset.mem_Icc, not_and_or, not_le] at hz
  unfold schrZ
  rw [ite_eq_right (by omega), ite_eq_right (by rw [abs_eq zero_le_one]; omega)]

/-- The recursion for the powers of `H_w`. -/
lemma schrZpow_succ (w : ℤ → ℝ) (k : ℕ) (x y : ℤ) :
    schrZpow w (k + 1) x y =
      -schrZpow w k (x - 1) y + w x * schrZpow w k x y - schrZpow w k (x + 1) y := by
  simp only [schrZpow]
  rw [sum_Icc_three, schrZ_left, schrZ_diag, schrZ_right]
  ring

/-- Entries of the powers of `H_w`: `|H_w^k(x, y)| ≤ (2 + W)^k`. -/
lemma abs_schrZpow_le (w : ℤ → ℝ) (W : ℝ) (hw : ∀ x, |w x| ≤ W) (k : ℕ) (x y : ℤ) :
    |schrZpow w k x y| ≤ (2 + W) ^ k := by
  have hW : 0 ≤ W := (abs_nonneg _).trans (hw 0)
  induction k generalizing x with
  | zero =>
    simp only [schrZpow, pow_zero]
    split_ifs <;> simp
  | succ k ih =>
    rw [schrZpow_succ, pow_succ]
    have h1 := abs_le.mp (ih (x - 1))
    have h3 := abs_le.mp (ih (x + 1))
    have hp : |w x * schrZpow w k x y| ≤ W * (2 + W) ^ k := by
      rw [abs_mul]
      exact mul_le_mul (hw x) (ih x) (abs_nonneg _) hW
    have h2 := abs_le.mp hp
    rw [abs_le]
    constructor <;> nlinarith

/-- Finite propagation on `ℤ`: `H_w^k(z, x) = 0` for `|z - x| > k`. -/
lemma schrZpow_eq_zero (w : ℤ → ℝ) (k : ℕ) (z x : ℤ) (h : (k : ℤ) < |z - x|) :
    schrZpow w k z x = 0 := by
  induction k generalizing z with
  | zero =>
    simp only [schrZpow]
    rw [ite_eq_right]
    intro e
    rw [e, sub_self, abs_zero] at h
    simp at h
  | succ k ih =>
    rw [lt_abs] at h
    push_cast at h
    rw [schrZpow_succ, ih (z - 1) (by rw [lt_abs]; omega), ih z (by rw [lt_abs]; omega),
      ih (z + 1) (by rw [lt_abs]; omega)]
    ring

/-- Translation covariance of the powers of `H_w`. -/
lemma schrZpow_shift (w : ℤ → ℝ) (x : ℤ) (k : ℕ) (y z : ℤ) :
    schrZpow w k (y + x) (z + x) = schrZpow (fun u => w (u + x)) k y z := by
  induction k generalizing y with
  | zero =>
    simp only [schrZpow, add_left_inj]
  | succ k ih =>
    rw [schrZpow_succ, schrZpow_succ, show y + x - 1 = (y - 1) + x by ring,
      show y + x + 1 = (y + 1) + x by ring, ih, ih, ih]

/-- Paths of length `k ≤ min(b, -a)` ending at `0` stay in `[a, b]`: the column `0` of `H_Λ^k` is
that of `H_ℤ^k`. -/
lemma pow_schrH_origin (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (w : ℤ → ℝ) (k : ℕ)
    (hk : (k : ℤ) ≤ min b (-a)) (y : Set.Icc a b) :
    ((schrH a b w) ^ k) y (origin a b ha hb) = schrZpow w k y.1 0 := by
  induction k generalizing y with
  | zero =>
    simp only [pow_zero, Matrix.one_apply, schrZpow]
    by_cases h : y.1 = 0
    · rw [ite_eq_left (Subtype.ext h), ite_eq_left h]
    · rw [ite_eq_right (fun e => h (congrArg Subtype.val e)), ite_eq_right h]
  | succ k ih =>
    have hk' : (k : ℤ) ≤ min b (-a) := by push_cast at hk; linarith
    have hkb : (k : ℤ) + 1 ≤ b := by push_cast at hk; exact hk.trans (min_le_left _ _)
    have hka : (k : ℤ) + 1 ≤ -a := by push_cast at hk; exact hk.trans (min_le_right _ _)
    rw [pow_succ', Matrix.mul_apply]
    simp_rw [ih hk']
    have e1 : ∑ u : Set.Icc a b, schrH a b w y u * schrZpow w k u.1 0 =
        ∑ z ∈ Finset.Icc a b, schrZ w y.1 z * schrZpow w k z 0 :=
      (Finset.sum_subtype (Finset.Icc a b) (p := fun z => z ∈ Set.Icc a b) (fun z => by simp)
        (fun z => schrZ w y.1 z * schrZpow w k z 0)).symm
    rw [e1]
    have hy1 := y.2.1
    have hy2 := y.2.2
    have s1 : ∀ z ∈ Finset.Icc (a - 1) (b + 1), z ∉ Finset.Icc a b →
        schrZ w y.1 z * schrZpow w k z 0 = 0 := by
      intro z hz hz'
      simp only [Finset.mem_Icc] at hz hz'
      rw [schrZpow_eq_zero w k z 0 (by rw [sub_zero, lt_abs]; omega), mul_zero]
    have s2 : ∀ z ∈ Finset.Icc (a - 1) (b + 1), z ∉ Finset.Icc (y.1 - 1) (y.1 + 1) →
        schrZ w y.1 z * schrZpow w k z 0 = 0 := by
      intro z _ hz
      rw [schrZ_eq_zero w y.1 z hz, zero_mul]
    rw [Finset.sum_subset (Finset.Icc_subset_Icc (show a - 1 ≤ a by omega)
      (show b ≤ b + 1 by omega)) s1, ← Finset.sum_subset (Finset.Icc_subset_Icc
      (show a - 1 ≤ y.1 - 1 by omega) (show y.1 + 1 ≤ b + 1 by omega)) s2]
    rfl

/-- `|⟨δ_y, e^{-isH_Λ} δ_x⟩| ≤ 1`. -/
lemma norm_evo_schrH_le_one (a b : ℤ) (w : ℤ → ℝ) (s : ℝ) (y x : Set.Icc a b) :
    ‖evo (schrH a b w) s y x‖ ≤ 1 := by
  have h1 := DimerTLAux.sum_norm_sq_evo _ (DimerTLAux.schrH_isSymm a b w) s x
  have h2 : ‖evo (schrH a b w) s y x‖ ^ 2 ≤ 1 :=
    h1 ▸ Finset.single_le_sum (f := fun u => ‖evo (schrH a b w) s u x‖ ^ 2)
      (fun _ _ => by positivity) (Finset.mem_univ y)
  nlinarith [norm_nonneg (evo (schrH a b w) s y x)]

/-- The tails of the exponential series tend to zero. -/
lemma tendsto_tail (c : ℝ) (hc : 0 ≤ c) :
    Tendsto (fun m : ℕ => ∑' k : ℕ, if m + 1 ≤ k then c ^ k / (k.factorial : ℝ) else 0) atTop
      (𝓝 0) := by
  have h1 := (FloorSemiring.tendsto_pow_div_factorial_atTop c).comp (tendsto_add_atTop_nat 1)
  have h2 := h1.mul_const (Real.exp c)
  rw [zero_mul] at h2
  refine squeeze_zero (fun m => ?_) (fun m => ?_) h2
  · exact tsum_nonneg fun k => by split_ifs <;> positivity
  · exact DimerTLAux.tail_exp_le c hc (m + 1)

end DimerInfVolAux

/-- The exponential series defining `evoZw` converges absolutely for a bounded potential. -/
theorem summable_evoZw (w : ℤ → ℝ) (W : ℝ) (hw : ∀ x, |w x| ≤ W) (s : ℝ) (x y : ℤ) :
    Summable (fun k : ℕ =>
      ((-(s : ℂ) * Complex.I) ^ k / (k.factorial : ℂ)) * ((schrZpow w k x y : ℝ) : ℂ)) := by
  exact Summable.of_norm_bounded (Real.summable_pow_div_factorial ((2 + W) * |s|))
    fun k => DimerTLAux.norm_term_le s k _ (2 + W) (DimerInfVolAux.abs_schrZpow_le w W hw k x y)

/-- Translation covariance: `⟨δ_{y+x}, e^{-isH_w} δ_x⟩ = ⟨δ_y, e^{-isH_{w(·+x)}} δ_0⟩`. -/
theorem evoZw_shift (w : ℤ → ℝ) (s : ℝ) (x y : ℤ) :
    evoZw w s (y + x) x = evoZw (fun z => w (z + x)) s y 0 := by
  unfold evoZw
  congr 1
  funext k
  have h := DimerInfVolAux.schrZpow_shift w x k y 0
  rw [zero_add] at h
  rw [h]

/-- **Thermodynamic limit** (paper: Lemma `lem:dimerTL`(ii) with `Λ' = ℤ`). -/
theorem dimerTL_ii_Z (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (w : ℤ → ℝ) (W : ℝ)
    (hw : ∀ x, |w x| ≤ W) (s : ℝ) (y : Set.Icc a b) :
    ‖evo (schrH a b w) s y (origin a b ha hb) - evoZw w s y.1 0‖ ≤
      2 * ∑' k : ℕ, if (min b (-a)).toNat + 1 ≤ k then ((2 + W) * |s|) ^ k / (k.factorial : ℝ)
        else 0 := by
  set m := (min b (-a)).toNat with hm
  have hW : 0 ≤ W := (abs_nonneg _).trans (hw 0)
  have h3 := (DimerTLAux.hasSum_evo_apply (schrH a b w) s y (origin a b ha hb)).sub
    (summable_evoZw w W hw s y.1 0).hasSum
  set f : ℕ → ℝ := fun k => if m + 1 ≤ k then ((2 + W) * |s|) ^ k / (k.factorial : ℝ) else 0
    with hf
  have hfs : Summable f := DimerTLAux.summable_tail _ (by positivity) _
  have hrow : ∀ u : Set.Icc a b, ∑ l, |schrH a b w u l| ≤ 2 + W := fun u =>
    (DimerTLAux.sum_abs_schrH_le a b w u).trans (by linarith [hw u.1])
  refine (h3.norm_le_of_bounded (hfs.mul_left 2).hasSum fun k => ?_).trans
    (le_of_eq tsum_mul_left)
  by_cases hk : m + 1 ≤ k
  · have e1 := DimerTLAux.norm_term_le s k _ (2 + W)
      (DimerTLAux.abs_pow_apply_le _ _ (by linarith) hrow k y (origin a b ha hb))
    have e2 := DimerTLAux.norm_term_le s k _ (2 + W)
      (DimerInfVolAux.abs_schrZpow_le w W hw k y.1 0)
    simp only [hf]
    rw [ite_eq_left hk]
    refine (norm_sub_le _ _).trans ?_
    linarith
  · have hk' : (k : ℤ) ≤ min b (-a) := by
      have h0 : (0 : ℤ) ≤ min b (-a) := le_min hb (by linarith)
      have : (k : ℤ) ≤ m := by exact_mod_cast (by omega : k ≤ m)
      rwa [hm, Int.toNat_of_nonneg h0] at this
    rw [DimerInfVolAux.pow_schrH_origin a b ha hb w k hk' y, sub_self, norm_zero]
    positivity

/-- `|⟨δ_y, e^{-isH_ℤ} δ_0⟩| ≤ 1`. -/
theorem norm_evoZw_le_one (w : ℤ → ℝ) (W : ℝ) (hw : ∀ x, |w x| ≤ W) (s : ℝ) (y : ℤ) :
    ‖evoZw w s y 0‖ ≤ 1 := by
  have hW : 0 ≤ W := (abs_nonneg _).trans (hw 0)
  have hτ := DimerInfVolAux.tendsto_tail ((2 + W) * |s|) (by positivity)
  have hlim := (hτ.const_mul 2).const_add 1
  rw [mul_zero, add_zero] at hlim
  refine ge_of_tendsto hlim ?_
  filter_upwards [eventually_ge_atTop y.natAbs] with n hn
  have hy : -(n : ℤ) ≤ y ∧ y ≤ n := by omega
  have h1 := dimerTL_ii_Z (-(n : ℤ)) n (by omega) (by omega) w W hw s ⟨y, hy⟩
  have h2 := DimerInfVolAux.norm_evo_schrH_le_one (-(n : ℤ)) n w s ⟨y, hy⟩
    (origin (-(n : ℤ)) n (by omega) (by omega))
  have hmn : (min (n : ℤ) (-(-(n : ℤ)))).toNat = n := by simp
  rw [hmn] at h1
  exact (norm_le_norm_add_norm_sub _ _).trans (add_le_add h2 h1)

/-! ### Moments -/

/-- **Ballistic bound for all moments**, uniformly in the interval and the potential:
`∑_y |y|^p |⟨δ_y, e^{-isH_Λ} δ_0⟩|² ≤ (4e²|s|)^p + ∑_{y ∈ ℤ} |y|^p 4^{-|y|}`. -/
theorem moment_poly_bound (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (w : ℤ → ℝ) (s p : ℝ) (hp : 0 < p) :
    ∑ y : Set.Icc a b, |((y.1 : ℤ) : ℝ)| ^ p * ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2 ≤
      (4 * Real.exp 2 * |s|) ^ p + ∑' y : ℤ, |(y : ℝ)| ^ p * (1 / 4 : ℝ) ^ y.natAbs := by
  sorry

/-- The moments converge in the thermodynamic limit `[a_n, b_n] ↑ ℤ`. -/
theorem tendsto_moment_Z (w : ℤ → ℝ) (W : ℝ) (hw : ∀ x, |w x| ≤ W) (as bs : ℕ → ℤ)
    (ha : ∀ n, as n ≤ 0) (hb : ∀ n, 0 ≤ bs n) (has : Tendsto as atTop atBot)
    (hbs : Tendsto bs atTop atTop) (p s : ℝ) (hp : 0 < p) :
    Tendsto (fun n => ∑ y : Set.Icc (as n) (bs n), |((y.1 : ℤ) : ℝ)| ^ p *
        ‖evo (schrH (as n) (bs n) w) s y (origin (as n) (bs n) (ha n) (hb n))‖ ^ 2) atTop
      (𝓝 (∑' y : ℤ, |(y : ℝ)| ^ p * ‖evoZw w s y 0‖ ^ 2)) := by
  sorry

/-- The time-averaged moments converge in the thermodynamic limit `[a_n, b_n] ↑ ℤ`. -/
theorem tendsto_dimerMom_Z (w : ℤ → ℝ) (W : ℝ) (hw : ∀ x, |w x| ≤ W) (as bs : ℕ → ℤ)
    (ha : ∀ n, as n ≤ 0) (hb : ∀ n, 0 ≤ bs n) (has : Tendsto as atTop atBot)
    (hbs : Tendsto bs atTop atTop) (p T : ℝ) (hp : 0 < p) (hT : 0 < T) :
    Tendsto (fun n => dimerMom (as n) (bs n) w p T (origin (as n) (bs n) (ha n) (hb n))) atTop
      (𝓝 (dimerMomZ w p T 0)) := by
  sorry

/-- **Every dimer configuration on `ℤ` spreads at least diffusively** (paper: Proposition `prop:dimerlower`
for `I = ℤ`). -/
theorem dimer_lower_bound_Z (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (p : ℝ) (hp : 0 < p) (r : ℤ)
    (v : ℤ → ℝ) (hv : IsDimerCfgZ lam r v) (x : ℤ) (T : ℝ) (hT : 1 ≤ T) :
    dimerClp lam p * T ^ (p - 1) ≤ dimerMomZ v p T x := by
  sorry

/-! ### The gauge theory -/

/-- `v^{(r)}` is a dimer configuration with values `±2|h|/J` and alignment `r`. -/
theorem isDimerCfgZ_dimerPotZ (J h : ℝ) (hJ : 0 < J) (r : ℤ) (η : CfgZ) :
    IsDimerCfgZ (2 * |h| / J) r (dimerPotZ J h r η) := by
  sorry

/-- **Every dimer sector, averaged** (paper: Theorem `thm:dimer`(a), second statement, in sector form):
`M^{(N)}_p(T; x) ≥ c_{λ,p} (JT)^{p-1}` for `1 ≤ JT ≤ min(x, N-1-x)` (0-indexed sites). -/
theorem dimerMN_lower (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) {N : ℕ}
    (x : Fin N) (p T : ℝ) (hp : 0 < p) (hT1 : 1 ≤ J * T) (hT2 : J * T ≤ x.val)
    (hT3 : J * T ≤ (N : ℝ) - 1 - x.val) :
    dimerClp (2 * |h| / J) p * (J * T) ^ (p - 1) ≤ dimerMN J h N p T x := by
  sorry

/-- **Thermodynamic limit of the dimer moments** (paper: Theorem `thm:dimer`(b), first limit). -/
theorem dimer_TL_D (J h : ℝ) (hJ : 0 < J) (r : ℤ) (Ns : ℕ → ℕ) (xs : ∀ n, Fin (Ns n))
    (hN : ∀ n, Even (Ns n)) (hpar : ∀ n, ((xs n).val : ℤ) % 2 = r)
    (hx1 : Tendsto (fun n => (xs n).val) atTop atTop)
    (hx2 : Tendsto (fun n => Ns n - (xs n).val) atTop atTop) (p t : ℝ) (hp0 : 0 < p)
    (hp2 : p ≤ 2) :
    Tendsto (fun n => dimerD J h (Ns n) p t (xs n)) atTop
      (𝓝 (∫ η, ∑' y : ℤ, |(y : ℝ)| ^ p * ‖evoZw (dimerPotZ J h r η) (J * t) y 0‖ ^ 2 ∂muZ)) := by
  sorry

/-- **Thermodynamic limit of the time-averaged dimer moments** (paper: Theorem `thm:dimer`(b), second
limit). -/
theorem dimer_TL_M (J h : ℝ) (hJ : 0 < J) (r : ℤ) (Ns : ℕ → ℕ) (xs : ∀ n, Fin (Ns n))
    (hN : ∀ n, Even (Ns n)) (hpar : ∀ n, ((xs n).val : ℤ) % 2 = r)
    (hx1 : Tendsto (fun n => (xs n).val) atTop atTop)
    (hx2 : Tendsto (fun n => Ns n - (xs n).val) atTop atTop) (p T : ℝ) (hp0 : 0 < p)
    (hp2 : p ≤ 2) (hT : 0 < T) :
    Tendsto (fun n => dimerMN J h (Ns n) p T (xs n)) atTop (𝓝 (dimerMInf J h r p T)) := by
  sorry

/-- **Superdiffusion, lower bound** (paper: Theorem `thm:dimer`(c), first statement):
`M^{∞,r}_p(T) ≥ c_{λ,p} (JT)^{p-1}` for `JT ≥ 1`. -/
theorem dimerMInf_lower (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (r : ℤ)
    (hr : r = 0 ∨ r = 1) (p T : ℝ) (hp0 : 0 < p) (hp2 : p ≤ 2) (hT : 1 ≤ J * T) :
    dimerClp (2 * |h| / J) p * (J * T) ^ (p - 1) ≤ dimerMInf J h r p T := by
  sorry

/-- **Failure of uniform localization for the dimer law** (paper: Theorem `thm:dimer`(d),
eq. `eq:dimerLGTd`, for `p > 1`): `sup_{N even, x, t} D^{(N)}_p(t; x) = ∞`. -/
theorem dimerD_unbounded (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (p : ℝ)
    (hp : 1 < p) :
    ¬ BddAbove (Set.range (fun r : (Σ M : ℕ, Fin (2 * M) × ℝ) =>
      dimerD J h (2 * r.1) p r.2.2 r.2.1)) := by
  sorry

end DFL
