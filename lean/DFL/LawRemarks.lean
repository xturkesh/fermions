import DFL.ProductLaw
import DFL.InfVol

/-!
# Sector laws of biased product states; transfer of upper bounds; the dimer law on `ℤ`
(paper: Remarks `rem:productlaws`(b),(c) and `rem:biased`)

Biased product states `χ_j = cos θ |+⟩ + sin θ |-⟩` on all links with Fock matter `|n⟩` have the sector law
`μ_θ(q) = |χ(x(n, q))|²` (`lem:productlaw`(b), `sector_average_prod_fock`): the `σˣ` labels `x_j` of the links
are i.i.d. with `P(x_j = 0) = ϱ = cos²θ`, and `q_j = x_{j-1} + x_j + n_j` (0-indexed, `x_{-1} = 0`). With
`m = 2ϱ - 1 = cos 2θ` and the signs `s_j = (-1)^{n_j} q_j`:

* `biasedE_first`, `biasedE_bulk`, `biasedE_pair01`, `biased_cov01`: `E s_0 = m`, `E s_j = m²` (`j ≥ 1`),
  `E s_0 s_1 = m`, so `Cov(s_0, s_1) = m(1 - m²)`: the charges are independent only if
  `cos²θ ∈ {0, 1/2, 1}` (`rem:productlaws`(b));
* `biased_cov_bulk`: `Cov(q_j, q_{j+1}) = ±m²(1 - m²)` for `j ≥ 1` (`rem:biased`);
* `biased_indep_far`: `q_j` and `q_k` are independent for `|j - k| ≥ 2` (a block factor of an i.i.d.
  sequence);
* `biasedLaw_le`, `exists_biasedLaw_eq`, `biasedLaw_pos`, `biasedLaw_half`, `biasedLaw_pure`: the density with
  respect to the uniform law is at most `(2 max(ϱ, 1 - ϱ))^N`, with equality somewhere; full support for
  `0 < ϱ < 1`; the uniform law for `ϱ = 1/2`; a point mass for `ϱ = 1`;
* `upper_bound_transfer` (`rem:productlaws`(c)): upper bounds on averages of non-negative sector functionals
  transfer to every law with density at most `K`;
* `muZ_dimerEvent`, `muDimZ_dimerEvent_compl`, `muZ_mutuallySingular_muDimZ` (`rem:productlaws`(c)): on `ℤ`
  the i.i.d. law and the dimer law are mutually singular.
-/

open MeasureTheory Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- `σˣ` amplitudes of the biased link state `cos θ |+⟩ + sin θ |-⟩` (label `0` for `|+⟩`). -/
def biasAmp (θ : ℝ) (x : ZMod 2) : ℂ := if x = 0 then (Real.cos θ : ℂ) else (Real.sin θ : ℂ)

/-- The sector law of the biased product state with Fock matter `n`: `μ_θ(q) = |χ(x(n, q))|²`. -/
def biasedLaw (θ : ℝ) (n q : Cfg N) : ℝ :=
  Complex.normSq (prodLink (fun _ => biasAmp θ) (linkOf n q))

/-- The expectation under the biased sector law. -/
def biasedE (θ : ℝ) (n : Cfg N) (f : Cfg N → ℝ) : ℝ := ∑ q, biasedLaw θ n q * f q

/-- The dimer law on `ℤ` (alignment `0`): the image of the i.i.d. law under `η ↦ (y ↦ η ⌊y/2⌋)`. -/
def muDimZ : Measure CfgZ := muZ.map (fun η : CfgZ => fun y : ℤ => η (y / 2))

/-! ### The biased law -/

/-- `μ_θ(q) = ∏_j (cos²θ or sin²θ)` according to the link label `x_j(n, q)`. -/
theorem biasedLaw_eq (θ : ℝ) (n q : Cfg N) :
    biasedLaw θ n q =
      ∏ j, (if linkOf n q j = 0 then Real.cos θ ^ 2 else Real.sin θ ^ 2) := by
  sorry

/-- `μ_θ` is a probability law. -/
theorem sum_biasedLaw (θ : ℝ) (n : Cfg N) : ∑ q, biasedLaw θ n q = 1 := by
  sorry

/-- `E s_0 = m = cos 2θ`. -/
theorem biasedE_first (θ : ℝ) (n : Cfg N) (hN : 0 < N) :
    biasedE θ n (fun q => sgn (n ⟨0, hN⟩) * sgn (q ⟨0, hN⟩)) = Real.cos (2 * θ) := by
  sorry

/-- `E s_j = m² = cos² 2θ` for `j ≥ 1`. -/
theorem biasedE_bulk (θ : ℝ) (n : Cfg N) (j : Fin N) (hj : 0 < j.val) :
    biasedE θ n (fun q => sgn (n j) * sgn (q j)) = Real.cos (2 * θ) ^ 2 := by
  sorry

/-- `E s_0 s_1 = m`. -/
theorem biasedE_pair01 (θ : ℝ) (n : Cfg N) (hN : 1 < N) :
    biasedE θ n (fun q => sgn (n ⟨0, by omega⟩) * sgn (q ⟨0, by omega⟩) *
      (sgn (n ⟨1, hN⟩) * sgn (q ⟨1, hN⟩))) = Real.cos (2 * θ) := by
  sorry

/-- `Cov(s_0, s_1) = m (1 - m²)`, nonzero unless `cos²θ ∈ {0, 1/2, 1}` (paper: `rem:productlaws`(b)). -/
theorem biased_cov01 (θ : ℝ) (n : Cfg N) (hN : 1 < N) :
    biasedE θ n (fun q => sgn (n ⟨0, by omega⟩) * sgn (q ⟨0, by omega⟩) *
        (sgn (n ⟨1, hN⟩) * sgn (q ⟨1, hN⟩))) -
      biasedE θ n (fun q => sgn (n ⟨0, by omega⟩) * sgn (q ⟨0, by omega⟩)) *
        biasedE θ n (fun q => sgn (n ⟨1, hN⟩) * sgn (q ⟨1, hN⟩)) =
      Real.cos (2 * θ) * (1 - Real.cos (2 * θ) ^ 2) := by
  sorry

/-- `Cov(q_j, q_{j+1}) = (-1)^{n_j + n_{j+1}} m² (1 - m²)` for `j ≥ 1` (paper: `rem:biased`). -/
theorem biased_cov_bulk (θ : ℝ) (n : Cfg N) (j : ℕ) (hj : 1 ≤ j) (hjN : j + 1 < N) :
    biasedE θ n (fun q => sgn (q ⟨j, by omega⟩) * sgn (q ⟨j + 1, hjN⟩)) -
      biasedE θ n (fun q => sgn (q ⟨j, by omega⟩)) * biasedE θ n (fun q => sgn (q ⟨j + 1, hjN⟩)) =
      sgn (n ⟨j, by omega⟩) * sgn (n ⟨j + 1, hjN⟩) *
        (Real.cos (2 * θ) ^ 2 * (1 - Real.cos (2 * θ) ^ 2)) := by
  sorry

/-- `q_j` and `q_k` are independent for `k ≥ j + 2` (paper: `rem:productlaws`(b)). -/
theorem biased_indep_far (θ : ℝ) (n : Cfg N) (j k : Fin N) (hjk : j.val + 2 ≤ k.val)
    (f g : ZMod 2 → ℝ) :
    biasedE θ n (fun q => f (q j) * g (q k)) =
      biasedE θ n (fun q => f (q j)) * biasedE θ n (fun q => g (q k)) := by
  sorry

/-- The density with respect to the uniform law is at most `(2 max(ϱ, 1 - ϱ))^N` (paper: `rem:biased`). -/
theorem biasedLaw_le (θ : ℝ) (n q : Cfg N) :
    biasedLaw θ n q ≤ max (Real.cos θ ^ 2) (Real.sin θ ^ 2) ^ N := by
  sorry

/-- The maximal density is attained (paper: `rem:biased`). -/
theorem exists_biasedLaw_eq (θ : ℝ) (n : Cfg N) :
    ∃ q, biasedLaw θ n q = max (Real.cos θ ^ 2) (Real.sin θ ^ 2) ^ N := by
  sorry

/-- Full support for `0 < ϱ < 1` (paper: `rem:biased`). -/
theorem biasedLaw_pos (θ : ℝ) (h1 : Real.cos θ ≠ 0) (h2 : Real.sin θ ≠ 0) (n q : Cfg N) :
    0 < biasedLaw θ n q := by
  sorry

/-- For `ϱ = 1/2` the law is i.i.d. uniform (paper: `rem:productlaws`(b)). -/
theorem biasedLaw_half (θ : ℝ) (hθ : Real.cos θ ^ 2 = 1 / 2) (n q : Cfg N) :
    biasedLaw θ n q = (1 / 2) ^ N := by
  sorry

/-- For `ϱ = 1` the law is a point mass: a single sector (paper: `rem:productlaws`(b)). -/
theorem biasedLaw_pure (θ : ℝ) (hθ : Real.sin θ = 0) (n q : Cfg N) :
    biasedLaw θ n q = if q = chargeOf n 0 then 1 else 0 := by
  sorry

/-! ### Transfer of upper bounds (paper: `rem:productlaws`(c)) -/

/-- An upper bound on a uniform average of a non-negative sector functional transfers to every law with
density at most `K` with respect to the uniform law. -/
theorem upper_bound_transfer (μ : Cfg N → ℝ) (K : ℝ) (hμ0 : ∀ q, 0 ≤ μ q)
    (hμK : ∀ q, μ q ≤ K * (1 / 2) ^ N) (F : Cfg N → ℝ) (hF : ∀ q, 0 ≤ F q) :
    ∑ q, μ q * F q ≤ K * avgQ F := by
  sorry

/-! ### The dimer law on `ℤ` (paper: `rem:productlaws`(c)) -/

/-- `muDimZ` is a probability measure. -/
theorem isProbabilityMeasure_muDimZ : IsProbabilityMeasure muDimZ := by
  sorry

/-- Under the i.i.d. law, almost surely not all pairs `(2k, 2k+1)` carry equal charges. -/
theorem muZ_dimerEvent : muZ {q : CfgZ | ∀ k : ℤ, q (2 * k) = q (2 * k + 1)} = 0 := by
  sorry

/-- Under the dimer law, almost surely all pairs `(2k, 2k+1)` carry equal charges. -/
theorem muDimZ_dimerEvent_compl :
    muDimZ {q : CfgZ | ∀ k : ℤ, q (2 * k) = q (2 * k + 1)}ᶜ = 0 := by
  sorry

/-- **In infinite volume the i.i.d. law and the dimer law are mutually singular** (paper:
`rem:productlaws`(c)). -/
theorem muZ_mutuallySingular_muDimZ : muZ ⟂ₘ muDimZ := by
  sorry

end DFL
