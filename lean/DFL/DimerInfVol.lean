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

namespace DimerInfVolAux

/-- `2^n (2ce)^n ≤ n!` for `n ≥ 4e²c` (from `n^n/n! ≤ e^n`). -/
lemma pow_le_factorial (c : ℝ) (hc : 0 ≤ c) (n : ℕ) (hn : 4 * Real.exp 2 * c ≤ n) :
    2 ^ n * (2 * c * Real.exp 1) ^ n ≤ (n.factorial : ℝ) := by
  have he : 0 < Real.exp 1 := Real.exp_pos 1
  have hf : (0 : ℝ) < n.factorial := by positivity
  have h1 : (n : ℝ) ^ n / n.factorial ≤ Real.exp n :=
    Real.pow_div_factorial_le_exp _ (Nat.cast_nonneg n) n
  have h2 : ((n : ℝ) / Real.exp 1) ^ n ≤ n.factorial := by
    rw [div_pow, Real.exp_one_pow, div_le_iff₀ (Real.exp_pos _)]
    rw [div_le_iff₀ hf] at h1
    linarith [mul_comm (Real.exp n) (n.factorial : ℝ)]
  have h3 : 2 * (2 * c * Real.exp 1) ≤ (n : ℝ) / Real.exp 1 := by
    rw [le_div_iff₀ he]
    have : Real.exp 2 = Real.exp 1 * Real.exp 1 := by rw [← Real.exp_add]; norm_num
    nlinarith
  rw [← mul_pow]
  exact (pow_le_pow_left₀ (by positivity) h3 n).trans h2

/-- `(2c)^n e^{2c}/n! ≤ 2^{-n}` for `n ≥ 4e²c`. -/
lemma decay_le (c : ℝ) (hc : 0 ≤ c) (n : ℕ) (hn : 4 * Real.exp 2 * c ≤ n) :
    (2 * c) ^ n / (n.factorial : ℝ) * Real.exp (2 * c) ≤ (1 / 2) ^ n := by
  have h1 := pow_le_factorial c hc n hn
  have hf : (0 : ℝ) < n.factorial := by positivity
  have hc2 : 2 * c ≤ n := by
    have : 1 ≤ Real.exp 2 := Real.one_le_exp (by norm_num)
    nlinarith
  have h2 : Real.exp (2 * c) ≤ Real.exp 1 ^ n := by
    rw [Real.exp_one_pow]
    exact Real.exp_le_exp.mpr hc2
  have h3 : (1 / 2 : ℝ) ^ n * 2 ^ n = 1 := by rw [← mul_pow]; norm_num
  have h4 : (2 * c * Real.exp 1) ^ n ≤ (1 / 2) ^ n * n.factorial := by
    calc (2 * c * Real.exp 1) ^ n = (1 / 2) ^ n * (2 ^ n * (2 * c * Real.exp 1) ^ n) := by
          rw [← mul_assoc, h3, one_mul]
      _ ≤ (1 / 2) ^ n * n.factorial := by gcongr
  calc (2 * c) ^ n / (n.factorial : ℝ) * Real.exp (2 * c)
      ≤ (2 * c) ^ n / (n.factorial : ℝ) * Real.exp 1 ^ n := by gcongr
    _ = (2 * c * Real.exp 1) ^ n / n.factorial := by ring
    _ ≤ (1 / 2) ^ n * n.factorial / n.factorial := div_le_div_of_nonneg_right h4 hf.le
    _ = (1 / 2) ^ n := by field_simp

/-- Far from the origin, `|⟨δ_y, e^{-isH_Λ} δ_0⟩|² ≤ 4^{-|y|}`. -/
lemma norm_sq_evo_le_far (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (w : ℤ → ℝ) (s : ℝ)
    (y : Set.Icc a b) (hy : 4 * Real.exp 2 * |s| ≤ (y.1.natAbs : ℝ)) :
    ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2 ≤ (1 / 4 : ℝ) ^ y.1.natAbs := by
  have h3 := (dimerTL_i' a b ha hb w s y).trans (decay_le |s| (abs_nonneg s) _ hy)
  calc ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2 ≤ ((1 / 2) ^ y.1.natAbs) ^ 2 :=
        pow_le_pow_left₀ (norm_nonneg _) h3 2
    _ = (1 / 4 : ℝ) ^ y.1.natAbs := by rw [← pow_mul, mul_comm, pow_mul]; norm_num

/-- `|⟨δ_y, e^{-isH_Λ} δ_0⟩|² ≤ 4^M 4^{-|y|}` with `M = ⌈4e²|s|⌉`, uniformly in `Λ` and `w`. -/
lemma norm_sq_evo_le (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (w : ℤ → ℝ) (s : ℝ)
    (y : Set.Icc a b) :
    ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2 ≤
      4 ^ ⌈4 * Real.exp 2 * |s|⌉₊ * (1 / 4 : ℝ) ^ y.1.natAbs := by
  by_cases hy : 4 * Real.exp 2 * |s| ≤ (y.1.natAbs : ℝ)
  · exact (norm_sq_evo_le_far a b ha hb w s y hy).trans
      (le_mul_of_one_le_left (by positivity) (one_le_pow₀ (by norm_num)))
  · push Not at hy
    have hn : y.1.natAbs ≤ ⌈4 * Real.exp 2 * |s|⌉₊ := by
      have := Nat.le_ceil (4 * Real.exp 2 * |s|)
      exact_mod_cast hy.le.trans this
    have h1 := norm_evo_schrH_le_one a b w s y (origin a b ha hb)
    calc ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2 ≤ 1 := by
          nlinarith [norm_nonneg (evo (schrH a b w) s y (origin a b ha hb))]
      _ = 4 ^ ⌈4 * Real.exp 2 * |s|⌉₊ * (1 / 4 : ℝ) ^ ⌈4 * Real.exp 2 * |s|⌉₊ := by
          rw [← mul_pow]; norm_num
      _ ≤ 4 ^ ⌈4 * Real.exp 2 * |s|⌉₊ * (1 / 4 : ℝ) ^ y.1.natAbs :=
          mul_le_mul_of_nonneg_left (pow_le_pow_of_le_one (by norm_num) (by norm_num) hn)
            (by positivity)

/-- `y ↦ |y|^p 4^{-|y|}` is summable on `ℤ`. -/
lemma summable_weight (p : ℝ) (hp : 0 < p) :
    Summable (fun y : ℤ => |(y : ℝ)| ^ p * (1 / 4 : ℝ) ^ y.natAbs) := by
  have hN : Summable (fun n : ℕ => (n : ℝ) ^ p * (1 / 4 : ℝ) ^ n) := by
    refine Summable.of_nonneg_of_le (fun n => by positivity) (fun n => ?_)
      (summable_pow_mul_geometric_of_norm_lt_one ⌈p⌉₊ (r := (1 / 4 : ℝ)) (by norm_num))
    refine mul_le_mul_of_nonneg_right ?_ (by positivity)
    rcases Nat.eq_zero_or_pos n with rfl | hn
    · rw [Nat.cast_zero, Real.zero_rpow hp.ne']
      positivity
    · rw [← Real.rpow_natCast]
      exact Real.rpow_le_rpow_of_exponent_le (by exact_mod_cast hn) (Nat.le_ceil p)
  refine Summable.of_nat_of_neg ?_ ?_
  · refine hN.congr fun n => ?_
    simp
  · refine hN.congr fun n => ?_
    simp

/-- A sum over `[a, b]` as a `tsum` over `ℤ` of the extension by zero. -/
lemma sum_Icc_eq_tsum (a b : ℤ) (f : Set.Icc a b → ℝ) :
    ∑ y, f y = ∑' y : ℤ, (if h : a ≤ y ∧ y ≤ b then f ⟨y, h⟩ else 0) := by
  rw [tsum_eq_sum (s := Finset.Icc a b) (fun y hy => by
      simp only [Finset.mem_Icc] at hy
      simp [hy]),
    Finset.sum_subtype (Finset.Icc a b) (p := fun z => z ∈ Set.Icc a b) (fun z => by simp)]
  refine Finset.sum_congr rfl fun y _ => ?_
  simp [y.2.1, y.2.2]

end DimerInfVolAux

/-- **Ballistic bound for all moments**, uniformly in the interval and the potential:
`∑_y |y|^p |⟨δ_y, e^{-isH_Λ} δ_0⟩|² ≤ (4e²|s|)^p + ∑_{y ∈ ℤ} |y|^p 4^{-|y|}`. -/
theorem moment_poly_bound (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (w : ℤ → ℝ) (s p : ℝ) (hp : 0 < p) :
    ∑ y : Set.Icc a b, |((y.1 : ℤ) : ℝ)| ^ p * ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2 ≤
      (4 * Real.exp 2 * |s|) ^ p + ∑' y : ℤ, |(y : ℝ)| ^ p * (1 / 4 : ℝ) ^ y.natAbs := by
  set R := 4 * Real.exp 2 * |s| with hR
  set g : ℤ → ℝ := fun y => |(y : ℝ)| ^ p * (1 / 4 : ℝ) ^ y.natAbs with hg
  have hg0 : ∀ y, 0 ≤ g y := fun y => by positivity
  have hgs : Summable g := DimerInfVolAux.summable_weight p hp
  have hpt : ∀ y : Set.Icc a b,
      |((y.1 : ℤ) : ℝ)| ^ p * ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2 ≤
        R ^ p * ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2 + g y.1 := by
    intro y
    have habs : |((y.1 : ℤ) : ℝ)| = (y.1.natAbs : ℝ) := by rw [Nat.cast_natAbs, Int.cast_abs]
    by_cases hy : R ≤ (y.1.natAbs : ℝ)
    · have h1 := DimerInfVolAux.norm_sq_evo_le_far a b ha hb w s y hy
      have h2 : |((y.1 : ℤ) : ℝ)| ^ p * ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2 ≤ g y.1 :=
        mul_le_mul_of_nonneg_left h1 (by positivity)
      have h3 : 0 ≤ R ^ p * ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2 := by positivity
      linarith
    · push Not at hy
      have h1 : |((y.1 : ℤ) : ℝ)| ^ p ≤ R ^ p := by
        rw [habs]
        exact Real.rpow_le_rpow (by positivity) hy.le hp.le
      have h2 := mul_le_mul_of_nonneg_right h1
        (sq_nonneg ‖evo (schrH a b w) s y (origin a b ha hb)‖)
      linarith [hg0 y.1]
  have hsub : ∑ y : Set.Icc a b, g y.1 ≤ ∑' y, g y := by
    rw [← Finset.sum_subtype (Finset.Icc a b) (p := fun z => z ∈ Set.Icc a b) (fun z => by simp) g]
    exact hgs.sum_le_tsum _ (fun y _ => hg0 y)
  calc ∑ y : Set.Icc a b, |((y.1 : ℤ) : ℝ)| ^ p * ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2
      ≤ ∑ y : Set.Icc a b, (R ^ p * ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2 + g y.1) :=
        Finset.sum_le_sum fun y _ => hpt y
    _ = R ^ p * ∑ y : Set.Icc a b, ‖evo (schrH a b w) s y (origin a b ha hb)‖ ^ 2 +
          ∑ y : Set.Icc a b, g y.1 := by
        rw [Finset.sum_add_distrib, Finset.mul_sum]
    _ = R ^ p + ∑ y : Set.Icc a b, g y.1 := by
        rw [DimerTLAux.sum_norm_sq_evo _ (DimerTLAux.schrH_isSymm a b w) s _, mul_one]
    _ ≤ R ^ p + ∑' y, g y := by linarith [hsub]

/-- The moments converge in the thermodynamic limit `[a_n, b_n] ↑ ℤ`. -/
theorem tendsto_moment_Z (w : ℤ → ℝ) (W : ℝ) (hw : ∀ x, |w x| ≤ W) (as bs : ℕ → ℤ)
    (ha : ∀ n, as n ≤ 0) (hb : ∀ n, 0 ≤ bs n) (has : Tendsto as atTop atBot)
    (hbs : Tendsto bs atTop atTop) (p s : ℝ) (hp : 0 < p) :
    Tendsto (fun n => ∑ y : Set.Icc (as n) (bs n), |((y.1 : ℤ) : ℝ)| ^ p *
        ‖evo (schrH (as n) (bs n) w) s y (origin (as n) (bs n) (ha n) (hb n))‖ ^ 2) atTop
      (𝓝 (∑' y : ℤ, |(y : ℝ)| ^ p * ‖evoZw w s y 0‖ ^ 2)) := by
  have hW : 0 ≤ W := (abs_nonneg _).trans (hw 0)
  set M := ⌈4 * Real.exp 2 * |s|⌉₊ with hM
  set G : ℕ → ℤ → ℂ := fun n y => if h : as n ≤ y ∧ y ≤ bs n then
      evo (schrH (as n) (bs n) w) s ⟨y, h⟩ (origin (as n) (bs n) (ha n) (hb n)) else 0 with hG
  have hF : ∀ n, ∑ y : Set.Icc (as n) (bs n), |((y.1 : ℤ) : ℝ)| ^ p *
      ‖evo (schrH (as n) (bs n) w) s y (origin (as n) (bs n) (ha n) (hb n))‖ ^ 2 =
      ∑' y : ℤ, |(y : ℝ)| ^ p * ‖G n y‖ ^ 2 := by
    intro n
    rw [DimerInfVolAux.sum_Icc_eq_tsum]
    congr 1
    funext y
    by_cases h : as n ≤ y ∧ y ≤ bs n
    · simp only [hG, dite_eq_left h]
    · simp only [hG, dite_eq_right h, norm_zero]
      simp
  simp_rw [hF]
  have hm : Tendsto (fun n => (min (bs n) (-as n)).toNat) atTop atTop := by
    rw [tendsto_atTop]
    intro K
    filter_upwards [hbs.eventually_ge_atTop (K : ℤ), has.eventually_le_atBot (-(K : ℤ))]
      with n h1 h2
    omega
  have hτ := ((DimerInfVolAux.tendsto_tail ((2 + W) * |s|) (by positivity)).comp hm).const_mul 2
  rw [mul_zero] at hτ
  refine tendsto_tsum_of_dominated_convergence
    (bound := fun y : ℤ => 4 ^ M * (|(y : ℝ)| ^ p * (1 / 4 : ℝ) ^ y.natAbs))
    ((DimerInfVolAux.summable_weight p hp).mul_left _) (fun y => ?_)
    (Eventually.of_forall fun n y => ?_)
  · have hGt : Tendsto (fun n => G n y) atTop (𝓝 (evoZw w s y 0)) := by
      rw [tendsto_iff_norm_sub_tendsto_zero]
      refine squeeze_zero' (Eventually.of_forall fun n => norm_nonneg _) ?_ hτ
      filter_upwards [has.eventually_le_atBot y, hbs.eventually_ge_atTop y] with n h1 h2
      have h : as n ≤ y ∧ y ≤ bs n := ⟨h1, h2⟩
      simp only [hG, dite_eq_left h]
      exact dimerTL_ii_Z (as n) (bs n) (ha n) (hb n) w W hw s ⟨y, h⟩
    exact ((tendsto_norm.comp hGt).pow 2).const_mul _
  · have h0 : 0 ≤ |(y : ℝ)| ^ p * ‖G n y‖ ^ 2 := by positivity
    rw [Real.norm_eq_abs, abs_of_nonneg h0]
    by_cases h : as n ≤ y ∧ y ≤ bs n
    · have h1 := DimerInfVolAux.norm_sq_evo_le (as n) (bs n) (ha n) (hb n) w s ⟨y, h⟩
      simp only [hG, dite_eq_left h]
      calc |(y : ℝ)| ^ p * ‖evo (schrH (as n) (bs n) w) s ⟨y, h⟩
            (origin (as n) (bs n) (ha n) (hb n))‖ ^ 2
          ≤ |(y : ℝ)| ^ p * (4 ^ M * (1 / 4 : ℝ) ^ y.natAbs) :=
            mul_le_mul_of_nonneg_left h1 (by positivity)
        _ = 4 ^ M * (|(y : ℝ)| ^ p * (1 / 4 : ℝ) ^ y.natAbs) := by ring
    · simp only [hG, dite_eq_right h, norm_zero]
      simp only [ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, mul_zero]
      positivity

namespace DimerInfVolAux

set_option backward.isDefEq.respectTransparency false in
open scoped Matrix.Norms.Operator in
/-- The propagator `t ↦ e^{-itA}` is continuous. -/
lemma continuous_evo {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℝ) :
    Continuous (fun t : ℝ => evo A t) := by
  unfold evo
  exact NormedSpace.exp_continuous.comp (by fun_prop)

/-- The time integrand of a moment is continuous. -/
lemma continuous_integrand {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℝ)
    (c : ι → ℝ) (x : ι) (T : ℝ) :
    Continuous (fun t : ℝ => (1 / T) * Real.exp (-t / T) * ∑ y, c y * ‖evo A t y x‖ ^ 2) := by
  have hU := continuous_evo A
  refine (continuous_const.mul (Real.continuous_exp.comp (continuous_neg.div_const T))).mul ?_
  exact continuous_finsetSum _ fun y _ => continuous_const.mul ((hU.matrix_elem y x).norm.pow 2)

/-- The dominating function `(1/T) e^{-t/T} ((4e²|t|)^p + ∑_y |y|^p 4^{-|y|})`. -/
def momBound (p T t : ℝ) : ℝ :=
  (1 / T) * Real.exp (-t / T) *
    ((4 * Real.exp 2 * |t|) ^ p + ∑' y : ℤ, |(y : ℝ)| ^ p * (1 / 4 : ℝ) ^ y.natAbs)

/-- The dominating function is integrable on `(0, ∞)`. -/
lemma integrableOn_momBound (p T : ℝ) (hp : 0 < p) (hT : 0 < T) :
    IntegrableOn (momBound p T) (Set.Ioi 0) := by
  set S := ∑' y : ℤ, |(y : ℝ)| ^ p * (1 / 4 : ℝ) ^ y.natAbs with hS
  have h1 : IntegrableOn (fun t : ℝ => t ^ p * Real.exp (-(1 / T) * t ^ (1 : ℝ))) (Set.Ioi 0) :=
    integrableOn_rpow_mul_exp_neg_mul_rpow (by linarith) one_pos (by positivity)
  have h2 := (h1.const_mul ((1 / T) * (4 * Real.exp 2) ^ p)).add
    ((DimerLGTAux.integrableOn_weight T hT).mul_const S)
  refine IntegrableOn.congr_fun (f := _) h2 (fun t ht => ?_) measurableSet_Ioi
  have ht' : 0 < t := ht
  simp only [momBound, Pi.add_apply, Real.rpow_one]
  rw [abs_of_pos ht', Real.mul_rpow (by positivity) ht'.le,
    show -(1 / T) * t = -t / T by ring]
  ring

/-- The integrand of `dimerMom` at the origin is dominated by `momBound`. -/
lemma integrand_le (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (w : ℤ → ℝ) (p T : ℝ) (hp : 0 < p)
    (hT : 0 < T) (t : ℝ) :
    ‖(1 / T) * Real.exp (-t / T) * ∑ y : Set.Icc a b, |((y.1 - (origin a b ha hb).1 : ℤ) : ℝ)| ^ p *
        ‖evo (schrH a b w) t y (origin a b ha hb)‖ ^ 2‖ ≤ momBound p T t := by
  have e : ∀ y : Set.Icc a b, y.1 - (origin a b ha hb).1 = y.1 := fun y => sub_zero _
  simp only [e]
  have h0 : 0 ≤ ∑ y : Set.Icc a b, |((y.1 : ℤ) : ℝ)| ^ p *
      ‖evo (schrH a b w) t y (origin a b ha hb)‖ ^ 2 :=
    Finset.sum_nonneg fun y _ => by positivity
  rw [Real.norm_eq_abs, abs_of_nonneg (by positivity), momBound]
  exact mul_le_mul_of_nonneg_left (moment_poly_bound a b ha hb w t p hp) (by positivity)

/-- `M_p(T; 0) ≤ ∫ momBound`, uniformly in the interval and the potential. -/
lemma norm_dimerMom_le (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (w : ℤ → ℝ) (p T : ℝ) (hp : 0 < p)
    (hT : 0 < T) :
    ‖dimerMom a b w p T (origin a b ha hb)‖ ≤ ∫ t in Set.Ioi (0 : ℝ), momBound p T t :=
  norm_integral_le_of_norm_le (integrableOn_momBound p T hp hT)
    (ae_of_all _ fun t => integrand_le a b ha hb w p T hp hT t)

end DimerInfVolAux

/-- The time-averaged moments converge in the thermodynamic limit `[a_n, b_n] ↑ ℤ`. -/
theorem tendsto_dimerMom_Z (w : ℤ → ℝ) (W : ℝ) (hw : ∀ x, |w x| ≤ W) (as bs : ℕ → ℤ)
    (ha : ∀ n, as n ≤ 0) (hb : ∀ n, 0 ≤ bs n) (has : Tendsto as atTop atBot)
    (hbs : Tendsto bs atTop atTop) (p T : ℝ) (hp : 0 < p) (hT : 0 < T) :
    Tendsto (fun n => dimerMom (as n) (bs n) w p T (origin (as n) (bs n) (ha n) (hb n))) atTop
      (𝓝 (dimerMomZ w p T 0)) := by
  unfold dimerMom dimerMomZ
  refine tendsto_integral_of_dominated_convergence (DimerInfVolAux.momBound p T)
    (fun n => (DimerInfVolAux.continuous_integrand _ _ _ T).aestronglyMeasurable)
    (DimerInfVolAux.integrableOn_momBound p T hp hT)
    (fun n => ae_of_all _ fun t =>
      DimerInfVolAux.integrand_le (as n) (bs n) (ha n) (hb n) w p T hp hT t)
    (ae_of_all _ fun t => ?_)
  have e : ∀ n, ∀ y : Set.Icc (as n) (bs n),
      y.1 - (origin (as n) (bs n) (ha n) (hb n)).1 = y.1 := fun n y => sub_zero _
  simp only [e, sub_zero]
  exact (tendsto_moment_Z w W hw as bs ha hb has hbs p t hp).const_mul _

namespace DimerInfVolAux

/-- Shifting a dimer configuration shifts its alignment. -/
lemma isDimerCfgZ_shift (lam : ℝ) (r : ℤ) (v : ℤ → ℝ) (hv : IsDimerCfgZ lam r v) (x : ℤ) :
    IsDimerCfgZ lam (r + x) (fun z => v (z + x)) := by
  refine ⟨fun y => hv.1 (y + x), fun k => ?_⟩
  have h := hv.2 k
  simp only
  rw [show 2 * k - (r + x) + x = 2 * k - r by ring,
    show 2 * k - (r + x) + 1 + x = 2 * k - r + 1 by ring]
  exact h

/-- Translation covariance of the time-averaged moments on `ℤ`. -/
lemma dimerMomZ_shift (v : ℤ → ℝ) (p T : ℝ) (x : ℤ) :
    dimerMomZ v p T x = dimerMomZ (fun z => v (z + x)) p T 0 := by
  unfold dimerMomZ
  refine setIntegral_congr_fun measurableSet_Ioi fun t _ => ?_
  rw [← (Equiv.addRight x).tsum_eq]
  simp only [Equiv.coe_addRight, add_sub_cancel_right, sub_zero, evoZw_shift]

/-- `prop:dimerlower` on `ℤ`, at the origin. -/
lemma lower_bound_origin (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (p : ℝ) (hp : 0 < p) (r : ℤ)
    (v : ℤ → ℝ) (hv : IsDimerCfgZ lam r v) (T : ℝ) (hT : 1 ≤ T) :
    dimerClp lam p * T ^ (p - 1) ≤ dimerMomZ v p T 0 := by
  have hW : ∀ z, |v z| ≤ lam := fun z => by
    rcases hv.1 z with h | h <;> rw [h] <;> simp [abs_of_pos h0]
  have hlim := tendsto_dimerMom_Z v lam hW (fun n : ℕ => -(n : ℤ)) (fun n => (n : ℤ))
    (fun n => by simp) (fun n => by positivity)
    (tendsto_neg_atTop_atBot.comp tendsto_natCast_atTop_atTop) tendsto_natCast_atTop_atTop p T hp
    (by linarith)
  refine ge_of_tendsto hlim ?_
  filter_upwards [eventually_ge_atTop ⌈T⌉₊] with n hn
  have hTn : T ≤ n := (Nat.le_ceil T).trans (by exact_mod_cast hn)
  exact dimer_lower_bound lam h0 h1 p hp r (-(n : ℤ)) n v
    ⟨fun y _ _ => hv.1 y, fun k _ _ => hv.2 k⟩ _ T hT
    ⟨by simp only [origin, Int.cast_neg, Int.cast_natCast, Int.cast_zero]; linarith,
      by simp only [origin, Int.cast_natCast, Int.cast_zero]; linarith⟩

end DimerInfVolAux

/-- **Every dimer configuration on `ℤ` spreads at least diffusively** (paper: Proposition `prop:dimerlower`
for `I = ℤ`). -/
theorem dimer_lower_bound_Z (lam : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (p : ℝ) (hp : 0 < p) (r : ℤ)
    (v : ℤ → ℝ) (hv : IsDimerCfgZ lam r v) (x : ℤ) (T : ℝ) (hT : 1 ≤ T) :
    dimerClp lam p * T ^ (p - 1) ≤ dimerMomZ v p T x := by
  rw [DimerInfVolAux.dimerMomZ_shift v p T x]
  exact DimerInfVolAux.lower_bound_origin lam h0 h1 p hp (r + x) _
    (DimerInfVolAux.isDimerCfgZ_shift lam r v hv x) T hT

/-! ### The gauge theory -/

/-- `v^{(r)}` is a dimer configuration with values `±2|h|/J` and alignment `r`. -/
theorem isDimerCfgZ_dimerPotZ (J h : ℝ) (hJ : 0 < J) (r : ℤ) (η : CfgZ) :
    IsDimerCfgZ (2 * |h| / J) r (dimerPotZ J h r η) := by
  refine ⟨fun y => ?_, fun k => ?_⟩
  · simp only [dimerPotZ, sgn]
    rcases abs_choice h with h' | h'
    · rw [h']
      split_ifs
      · left; ring
      · right; ring
    · rw [h']
      split_ifs
      · right; ring
      · left; ring
  · simp only [dimerPotZ]
    rw [show (2 * k - r + r) / 2 = k by omega, show (2 * k - r + 1 + r) / 2 = k by omega]

namespace DimerInfVolAux

/-- The moment of `h_q` is integrable against the time weight. -/
lemma integrable_moment_hq (J h : ℝ) {N : ℕ} (x : Fin N) (p T : ℝ) (hT : 0 < T) (q : Cfg N) :
    Integrable (fun t => (1 / T) * Real.exp (-t / T) *
      ∑ y : Fin N, |((y.val : ℝ) - x.val)| ^ p * ‖evo (hq J h q) t y x‖ ^ 2)
      (volume.restrict (Set.Ioi 0)) := by
  set C := ∑ y : Fin N, |((y.val : ℝ) - x.val)| ^ p with hC
  refine Integrable.mono' ((DimerLGTAux.integrableOn_weight T hT).mul_const C)
    (continuous_integrand _ _ _ T).aestronglyMeasurable (ae_of_all _ fun t => ?_)
  have h0 : 0 ≤ ∑ y : Fin N, |((y.val : ℝ) - x.val)| ^ p * ‖evo (hq J h q) t y x‖ ^ 2 :=
    Finset.sum_nonneg fun y _ => by positivity
  rw [Real.norm_eq_abs, abs_of_nonneg (by positivity)]
  refine mul_le_mul_of_nonneg_left (Finset.sum_le_sum fun y _ => ?_) (by positivity)
  have h1 := norm_evo_le_one (hq J h q) (hq_isSymm J h q) t y x
  have h2 : ‖evo (hq J h q) t y x‖ ^ 2 ≤ 1 := by
    nlinarith [norm_nonneg (evo (hq J h q) t y x)]
  calc |((y.val : ℝ) - x.val)| ^ p * ‖evo (hq J h q) t y x‖ ^ 2
      ≤ |((y.val : ℝ) - x.val)| ^ p * 1 := mul_le_mul_of_nonneg_left h2 (by positivity)
    _ = |((y.val : ℝ) - x.val)| ^ p := mul_one _

/-- The time average commutes with the dimer average. -/
lemma dimerMN_eq (J h : ℝ) (N : ℕ) (x : Fin N) (p T : ℝ) (hT : 0 < T) :
    dimerMN J h N p T x = dimerAvg (fun q => ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) *
        ∑ y : Fin N, |((y.val : ℝ) - x.val)| ^ p * ‖evo (hq J h q) t y x‖ ^ 2) := by
  unfold dimerMN dimerD dimerAvg
  rw [← integral_finsetSum _ (fun q _ => integrable_moment_hq J h x p T hT q), ← integral_div]
  congr 1
  funext t
  rw [mul_div_assoc', Finset.mul_sum]

end DimerInfVolAux

/-- **Every dimer sector, averaged** (paper: Theorem `thm:dimer`(a), second statement, in sector form):
`M^{(N)}_p(T; x) ≥ c_{λ,p} (JT)^{p-1}` for `1 ≤ JT ≤ min(x, N-1-x)` (0-indexed sites). -/
theorem dimerMN_lower (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) {N : ℕ}
    (x : Fin N) (p T : ℝ) (hp : 0 < p) (hT1 : 1 ≤ J * T) (hT2 : J * T ≤ x.val)
    (hT3 : J * T ≤ (N : ℝ) - 1 - x.val) :
    dimerClp (2 * |h| / J) p * (J * T) ^ (p - 1) ≤ dimerMN J h N p T x := by
  have hT : 0 < T := by nlinarith
  rw [DimerInfVolAux.dimerMN_eq J h N x p T hT]
  have hcard : 0 < ((univ.filter (fun q : Cfg N => IsDimer q)).card : ℝ) := by
    have h0 : (0 : Cfg N) ∈ univ.filter (fun q : Cfg N => IsDimer q) :=
      Finset.mem_filter.mpr ⟨Finset.mem_univ _, fun k _ hk => rfl⟩
    exact_mod_cast Finset.card_pos.mpr ⟨0, h0⟩
  unfold dimerAvg
  rw [le_div_iff₀ hcard]
  calc dimerClp (2 * |h| / J) p * (J * T) ^ (p - 1) *
        ((univ.filter (fun q : Cfg N => IsDimer q)).card : ℝ)
      = ∑ q ∈ univ.filter (fun q : Cfg N => IsDimer q),
          dimerClp (2 * |h| / J) p * (J * T) ^ (p - 1) := by
        rw [Finset.sum_const, nsmul_eq_mul, mul_comm]
    _ ≤ _ := Finset.sum_le_sum fun q hq =>
        dimer_sector_moment_lower J h hJ hh0 hh1 q (Finset.mem_filter.mp hq).2 x p T hp hT1 hT2
          hT3

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
  rintro ⟨C, hC⟩
  rw [mem_upperBounds] at hC
  have hl0 : 0 < 2 * |h| / J := by positivity
  have hl1 : 2 * |h| / J < 1 := by rw [div_lt_one hJ]; exact hh1
  have hc0 := DimerLGTAux.dimerClp_pos (2 * |h| / J) p hl0 hl1 (by linarith)
  obtain ⟨s, hs1, hsC⟩ := DimerLGTAux.exists_large p hp _ hc0 C
  have hsM : s ≤ ⌈s⌉₊ := Nat.le_ceil s
  have hJT : J * (s / J) = s := by field_simp
  set x : Fin (2 * (⌈s⌉₊ + 1)) := ⟨⌈s⌉₊ + 1, by omega⟩ with hx
  have hlow := dimerMN_lower J h hJ hh0 hh1 x p (s / J) (by linarith) (by rw [hJT]; exact hs1)
    (by rw [hJT]; show s ≤ ((⌈s⌉₊ + 1 : ℕ) : ℝ); push_cast; linarith)
    (by rw [hJT]; show s ≤ ((2 * (⌈s⌉₊ + 1) : ℕ) : ℝ) - 1 - ((⌈s⌉₊ + 1 : ℕ) : ℝ); push_cast
        linarith)
  have hup : dimerMN J h (2 * (⌈s⌉₊ + 1)) p (s / J) x ≤ C :=
    DimerLGTAux.integral_weight_le (s / J) C (by positivity)
      (fun t => dimerD J h (2 * (⌈s⌉₊ + 1)) p t x)
      (fun t => DimerLGTAux.dimerAvg_nonneg fun q => Finset.sum_nonneg fun y _ => by positivity)
      (fun t => hC _ ⟨⟨⌈s⌉₊ + 1, x, t⟩, rfl⟩)
  rw [hJT] at hlow
  linarith

end DFL
