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

/-- The exponential series defining `evoZw` converges absolutely for a bounded potential. -/
theorem summable_evoZw (w : ℤ → ℝ) (W : ℝ) (hw : ∀ x, |w x| ≤ W) (s : ℝ) (x y : ℤ) :
    Summable (fun k : ℕ =>
      ((-(s : ℂ) * Complex.I) ^ k / (k.factorial : ℂ)) * ((schrZpow w k x y : ℝ) : ℂ)) := by
  sorry

/-- Translation covariance: `⟨δ_{y+x}, e^{-isH_w} δ_x⟩ = ⟨δ_y, e^{-isH_{w(·+x)}} δ_0⟩`. -/
theorem evoZw_shift (w : ℤ → ℝ) (s : ℝ) (x y : ℤ) :
    evoZw w s (y + x) x = evoZw (fun z => w (z + x)) s y 0 := by
  sorry

/-- **Thermodynamic limit** (paper: Lemma `lem:dimerTL`(ii) with `Λ' = ℤ`). -/
theorem dimerTL_ii_Z (a b : ℤ) (ha : a ≤ 0) (hb : 0 ≤ b) (w : ℤ → ℝ) (W : ℝ)
    (hw : ∀ x, |w x| ≤ W) (s : ℝ) (y : Set.Icc a b) :
    ‖evo (schrH a b w) s y (origin a b ha hb) - evoZw w s y.1 0‖ ≤
      2 * ∑' k : ℕ, if (min b (-a)).toNat + 1 ≤ k then ((2 + W) * |s|) ^ k / (k.factorial : ℝ)
        else 0 := by
  sorry

/-- `|⟨δ_y, e^{-isH_ℤ} δ_0⟩| ≤ 1`. -/
theorem norm_evoZw_le_one (w : ℤ → ℝ) (W : ℝ) (hw : ∀ x, |w x| ≤ W) (s : ℝ) (y : ℤ) :
    ‖evoZw w s y 0‖ ≤ 1 := by
  sorry

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
