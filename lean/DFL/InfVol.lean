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

/-- The exponential series defining `UZ` converges absolutely. -/
theorem summable_UZ (J h : ℝ) (q : CfgZ) (t : ℝ) (x y : ℤ) :
    Summable (fun k : ℕ =>
      ((-(t : ℂ) * Complex.I) ^ k / (k.factorial : ℂ)) * ((hZpow J h q k x y : ℝ) : ℂ)) := by
  sorry

/-- **Finite propagation speed** (paper: Lemma `lem:propagation`(b)): the finite-volume propagator
entries converge to the entries of `U^ℤ_q(t)`. -/
theorem tendsto_evoZ (J h : ℝ) (q : CfgZ) (t : ℝ) (x y : ℤ) :
    Tendsto (fun M : ℕ => evoZ J h q t M x y) atTop (𝓝 (UZ J h q t x y)) := by
  sorry

/-- `|U^ℤ_q(t)_{xy}| ≤ 1`. -/
theorem norm_UZ_le_one (J h : ℝ) (q : CfgZ) (t : ℝ) (x y : ℤ) : ‖UZ J h q t x y‖ ≤ 1 := by
  sorry

/-- The finite-volume entries depend measurably on the charges. -/
theorem measurable_evoZ (J h : ℝ) (t : ℝ) (M : ℕ) (x y : ℤ) :
    Measurable (fun q : CfgZ => evoZ J h q t M x y) := by
  sorry

/-- The infinite-volume entries depend measurably on the charges. -/
theorem measurable_UZ (J h : ℝ) (t : ℝ) (x y : ℤ) :
    Measurable (fun q : CfgZ => UZ J h q t x y) := by
  sorry

/-- **Annealed dynamical localization on `ℤ`** (paper: Corollary `cor:ECinfvol`, Theorem
`thm:infvol`(b)): under Hypothesis (EC), `E sup_t |U^ℤ_q(t)_{xy}| ≤ F(|x-y|)`. -/
theorem lintegral_iSup_UZ_le (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (x y : ℤ) :
    ∫⁻ q, ⨆ t : ℝ, (‖UZ J h q t x y‖₊ : ENNReal) ∂muZ ≤ ENNReal.ofReal (F (x - y).natAbs) := by
  sorry

/-- **Annealed dynamical localization on `ℤ` from uniform large deviations**: there are `C` and
`η > 0` with `E sup_t |U^ℤ_q(t)_{xy}| ≤ C e^{-η |x-y|}`. -/
theorem infvol_propagator_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ ∀ x y : ℤ,
      ∫⁻ q, ⨆ t : ℝ, (‖UZ J h q t x y‖₊ : ENNReal) ∂muZ ≤
        ENNReal.ofReal (C * Real.exp (-η * (x - y).natAbs)) := by
  sorry

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
  sorry

/-- **Moments** (paper: Theorem `thm:infvol`(b)): under (EC), for `0 ≤ p ≤ 2`,
`∑_x (|x-y|+1)^p E |U^ℤ_q(t)_{xy}|² ≤ 2 M_p(F)`. -/
theorem moments_UZ (J h : ℝ) {F : ℕ → ℝ} (hEC : EC J h F) (t : ℝ) (y : ℤ) (p : ℝ) (hp0 : 0 ≤ p)
    (hp2 : p ≤ 2) :
    Summable (fun x : ℤ => (((x - y).natAbs : ℝ) + 1) ^ p * ∫ q, Complex.normSq (UZ J h q t x y) ∂muZ) ∧
      ∑' x : ℤ, (((x - y).natAbs : ℝ) + 1) ^ p * ∫ q, Complex.normSq (UZ J h q t x y) ∂muZ ≤
        2 * ∑' d : ℕ, ((d : ℝ) + 1) ^ p * F d := by
  sorry

end DFL
