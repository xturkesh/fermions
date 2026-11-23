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

namespace LawRemarksAux

/-- The weight `|⟨a|χ_j⟩|² ∈ {cos²θ, sin²θ}` of a single link label `a`. -/
def linkWt (θ : ℝ) (a : ZMod 2) : ℝ := if a = 0 then Real.cos θ ^ 2 else Real.sin θ ^ 2

/-- `|biasAmp θ a|² = linkWt θ a`. -/
lemma normSq_biasAmp (θ : ℝ) (a : ZMod 2) : Complex.normSq (biasAmp θ a) = linkWt θ a := by
  unfold biasAmp linkWt
  split_ifs <;> rw [Complex.normSq_ofReal, sq]

end LawRemarksAux

/-- `μ_θ(q) = ∏_j (cos²θ or sin²θ)` according to the link label `x_j(n, q)`. -/
theorem biasedLaw_eq (θ : ℝ) (n q : Cfg N) :
    biasedLaw θ n q =
      ∏ j, (if linkOf n q j = 0 then Real.cos θ ^ 2 else Real.sin θ ^ 2) := by
  unfold biasedLaw prodLink
  rw [map_prod]
  exact Finset.prod_congr rfl fun j _ => LawRemarksAux.normSq_biasAmp θ _

namespace LawRemarksAux

/-- A sum over `ZMod 2`. -/
lemma sum_zmod2 (F : ZMod 2 → ℝ) : ∑ a : ZMod 2, F a = F 0 + F 1 := by
  rw [show (Finset.univ : Finset (ZMod 2)) = {0, 1} by decide, Finset.sum_pair (by decide)]

/-- The single-link weights sum to one. -/
lemma sum_linkWt (θ : ℝ) : ∑ a : ZMod 2, linkWt θ a = 1 := by
  rw [sum_zmod2]
  simp only [linkWt, ite_eq_right (show (1 : ZMod 2) ≠ 0 by decide), ite_true]
  exact Real.cos_sq_add_sin_sq θ

/-- The mean sign of a link label is `cos 2θ`. -/
lemma sum_linkWt_sgn (θ : ℝ) : ∑ a : ZMod 2, linkWt θ a * sgn a = Real.cos (2 * θ) := by
  rw [sum_zmod2, Real.cos_two_mul']
  simp only [linkWt, ite_eq_right (show (1 : ZMod 2) ≠ 0 by decide), ite_true,
    FreeDynAux.sgn_zero, FreeDynAux.sgn_one]
  ring

/-- Expectations under `μ_θ` are expectations under the i.i.d. law of the link labels. -/
lemma biasedE_eq_sum (θ : ℝ) (n : Cfg N) (f : Cfg N → ℝ) :
    biasedE θ n f = ∑ x : Cfg N, (∏ j, linkWt θ (x j)) * f (chargeOf n x) := by
  unfold biasedE
  let e : Cfg N ≃ Cfg N :=
    ⟨linkOf n, chargeOf n, fun q => chargeOf_linkOf n q, fun x => linkOf_chargeOf n x⟩
  refine Fintype.sum_equiv e _ _ fun q => ?_
  simp only [e, Equiv.coe_fn_mk, chargeOf_linkOf, biasedLaw_eq]
  rfl

/-- The product weight factorizes against products of one-link functions. -/
lemma sum_linkWt_mul_prod (θ : ℝ) (g : Fin N → ZMod 2 → ℝ) :
    ∑ x : Cfg N, (∏ j, linkWt θ (x j)) * ∏ j, g j (x j) = ∏ j, ∑ a, linkWt θ a * g j a := by
  rw [Fintype.prod_sum]
  exact Finset.sum_congr rfl fun x _ => Finset.prod_mul_distrib.symm

/-- Sign moments of the link labels: `E ∏_{i ∈ S} (-1)^{x_i} = cos(2θ)^{|S|}`. -/
lemma sum_linkWt_mul_prod_sgn (θ : ℝ) (S : Finset (Fin N)) :
    ∑ x : Cfg N, (∏ j, linkWt θ (x j)) * ∏ i ∈ S, sgn (x i) = Real.cos (2 * θ) ^ S.card := by
  have h := sum_linkWt_mul_prod θ (fun i a => if i ∈ S then sgn a else 1)
  simp only [Finset.prod_ite_mem, Finset.univ_inter] at h
  rw [h]
  have h2 : ∀ i, ∑ a, linkWt θ a * (if i ∈ S then sgn a else 1) =
      if i ∈ S then Real.cos (2 * θ) else 1 := by
    intro i
    split_ifs
    · exact sum_linkWt_sgn θ
    · simp only [mul_one]; exact sum_linkWt θ
  simp only [h2, Finset.prod_ite_mem, Finset.univ_inter, Finset.prod_const]

/-- Expectation of a sign monomial in the link labels. -/
lemma biasedE_of_monomial (θ : ℝ) (n : Cfg N) (f : Cfg N → ℝ) (c : ℝ) (S : Finset (Fin N))
    (hf : ∀ x, f (chargeOf n x) = c * ∏ i ∈ S, sgn (x i)) :
    biasedE θ n f = c * Real.cos (2 * θ) ^ S.card := by
  rw [biasedE_eq_sum, ← sum_linkWt_mul_prod_sgn θ S, Finset.mul_sum]
  exact Finset.sum_congr rfl fun x _ => by rw [hf]; ring

/-- The previous link of a site `j` with left neighbour `i` is the link `i`. -/
lemma prevLink_of_succ (x : Cfg N) (i j : Fin N) (h : i.val + 1 = j.val) :
    prevLink x j = x i := by
  unfold prevLink
  rw [dite_eq_right (by omega)]
  exact congrArg x (Fin.ext (by show j.val - 1 = i.val; omega))

/-- `(-1)^{q_0} = (-1)^{n_0} (-1)^{x_0}` for `q = chargeOf n x`. -/
lemma sgn_chargeOf_zero (n x : Cfg N) (j : Fin N) (h : j.val = 0) :
    sgn (chargeOf n x j) = sgn (n j) * sgn (x j) := by
  have hp : prevLink x j = 0 := by unfold prevLink; rw [dite_eq_left h]
  rw [show chargeOf n x j = prevLink x j + x j + n j from rfl, hp, zero_add, FreeDynAux.sgn_add,
    mul_comm]

/-- `(-1)^{q_j} = (-1)^{n_j} (-1)^{x_{j-1}} (-1)^{x_j}` for `q = chargeOf n x` and `j ≥ 1`. -/
lemma sgn_chargeOf_succ (n x : Cfg N) (i j : Fin N) (h : i.val + 1 = j.val) :
    sgn (chargeOf n x j) = sgn (n j) * (sgn (x i) * sgn (x j)) := by
  rw [show chargeOf n x j = prevLink x j + x j + n j from rfl, prevLink_of_succ x i j h,
    FreeDynAux.sgn_add, FreeDynAux.sgn_add]
  ring

end LawRemarksAux

/-- `μ_θ` is a probability law. -/
theorem sum_biasedLaw (θ : ℝ) (n : Cfg N) : ∑ q, biasedLaw θ n q = 1 := by
  have h := LawRemarksAux.biasedE_of_monomial θ n (fun _ => 1) 1 ∅ (fun x => by simp)
  simpa [biasedE] using h

/-- `E s_0 = m = cos 2θ`. -/
theorem biasedE_first (θ : ℝ) (n : Cfg N) (hN : 0 < N) :
    biasedE θ n (fun q => sgn (n ⟨0, hN⟩) * sgn (q ⟨0, hN⟩)) = Real.cos (2 * θ) := by
  rw [LawRemarksAux.biasedE_of_monomial θ n _ 1 {⟨0, hN⟩}, Finset.card_singleton, pow_one,
    one_mul]
  intro x
  rw [LawRemarksAux.sgn_chargeOf_zero n x ⟨0, hN⟩ rfl, Finset.prod_singleton]
  linear_combination sgn (x ⟨0, hN⟩) * FreeDynAux.sgn_mul_self (n ⟨0, hN⟩)

/-- `E s_j = m² = cos² 2θ` for `j ≥ 1`. -/
theorem biasedE_bulk (θ : ℝ) (n : Cfg N) (j : Fin N) (hj : 0 < j.val) :
    biasedE θ n (fun q => sgn (n j) * sgn (q j)) = Real.cos (2 * θ) ^ 2 := by
  obtain ⟨i, hi⟩ : ∃ i : Fin N, i.val + 1 = j.val :=
    ⟨⟨j.val - 1, by have := j.isLt; omega⟩, by show j.val - 1 + 1 = j.val; omega⟩
  have hij : i ≠ j := fun e => by rw [e] at hi; omega
  rw [LawRemarksAux.biasedE_of_monomial θ n _ 1 {i, j}, Finset.card_pair hij, one_mul]
  intro x
  rw [LawRemarksAux.sgn_chargeOf_succ n x i j hi, Finset.prod_pair hij]
  linear_combination sgn (x i) * sgn (x j) * FreeDynAux.sgn_mul_self (n j)

/-- `E s_0 s_1 = m`. -/
theorem biasedE_pair01 (θ : ℝ) (n : Cfg N) (hN : 1 < N) :
    biasedE θ n (fun q => sgn (n ⟨0, by omega⟩) * sgn (q ⟨0, by omega⟩) *
      (sgn (n ⟨1, hN⟩) * sgn (q ⟨1, hN⟩))) = Real.cos (2 * θ) := by
  rw [LawRemarksAux.biasedE_of_monomial θ n _ 1 {⟨1, hN⟩}, Finset.card_singleton, pow_one,
    one_mul]
  intro x
  rw [LawRemarksAux.sgn_chargeOf_zero n x ⟨0, by omega⟩ rfl,
    LawRemarksAux.sgn_chargeOf_succ n x ⟨0, by omega⟩ ⟨1, hN⟩ rfl, Finset.prod_singleton]
  have h0 := FreeDynAux.sgn_mul_self (n ⟨0, by omega⟩)
  have h1 := FreeDynAux.sgn_mul_self (n ⟨1, hN⟩)
  have h2 := FreeDynAux.sgn_mul_self (x ⟨0, by omega⟩)
  linear_combination (sgn (n ⟨1, hN⟩) * sgn (n ⟨1, hN⟩) * sgn (x ⟨0, by omega⟩) *
      sgn (x ⟨0, by omega⟩) * sgn (x ⟨1, hN⟩)) * h0 +
    (sgn (x ⟨0, by omega⟩) * sgn (x ⟨0, by omega⟩) * sgn (x ⟨1, hN⟩)) * h1 +
    sgn (x ⟨1, hN⟩) * h2

/-- `Cov(s_0, s_1) = m (1 - m²)`, nonzero unless `cos²θ ∈ {0, 1/2, 1}` (paper: `rem:productlaws`(b)). -/
theorem biased_cov01 (θ : ℝ) (n : Cfg N) (hN : 1 < N) :
    biasedE θ n (fun q => sgn (n ⟨0, by omega⟩) * sgn (q ⟨0, by omega⟩) *
        (sgn (n ⟨1, hN⟩) * sgn (q ⟨1, hN⟩))) -
      biasedE θ n (fun q => sgn (n ⟨0, by omega⟩) * sgn (q ⟨0, by omega⟩)) *
        biasedE θ n (fun q => sgn (n ⟨1, hN⟩) * sgn (q ⟨1, hN⟩)) =
      Real.cos (2 * θ) * (1 - Real.cos (2 * θ) ^ 2) := by
  rw [biasedE_pair01 θ n hN, biasedE_first θ n (by omega), biasedE_bulk θ n ⟨1, hN⟩ Nat.one_pos]
  ring

namespace LawRemarksAux

/-- `E (-1)^{q_j} = (-1)^{n_j} cos² 2θ` for a site `j` with left neighbour `i`. -/
lemma biasedE_sgn_succ (θ : ℝ) (n : Cfg N) (i j : Fin N) (h : i.val + 1 = j.val) :
    biasedE θ n (fun q => sgn (q j)) = sgn (n j) * Real.cos (2 * θ) ^ 2 := by
  have hij : i ≠ j := fun e => by rw [e] at h; omega
  rw [biasedE_of_monomial θ n _ (sgn (n j)) {i, j}, Finset.card_pair hij]
  intro x
  rw [sgn_chargeOf_succ n x i j h, Finset.prod_pair hij]

/-- `E (-1)^{q_j + q_k} = (-1)^{n_j + n_k} cos² 2θ` for consecutive sites `i, j, k`. -/
lemma biasedE_sgn_pair_succ (θ : ℝ) (n : Cfg N) (i j k : Fin N) (hij : i.val + 1 = j.val)
    (hjk : j.val + 1 = k.val) :
    biasedE θ n (fun q => sgn (q j) * sgn (q k)) =
      sgn (n j) * sgn (n k) * Real.cos (2 * θ) ^ 2 := by
  have hik : i ≠ k := fun e => by rw [e] at hij; omega
  rw [biasedE_of_monomial θ n _ (sgn (n j) * sgn (n k)) {i, k}, Finset.card_pair hik]
  intro x
  rw [sgn_chargeOf_succ n x i j hij, sgn_chargeOf_succ n x j k hjk, Finset.prod_pair hik]
  linear_combination sgn (n j) * sgn (n k) * sgn (x i) * sgn (x k) * FreeDynAux.sgn_mul_self (x j)

end LawRemarksAux

/-- `Cov(q_j, q_{j+1}) = (-1)^{n_j + n_{j+1}} m² (1 - m²)` for `j ≥ 1` (paper: `rem:biased`). -/
theorem biased_cov_bulk (θ : ℝ) (n : Cfg N) (j : ℕ) (hj : 1 ≤ j) (hjN : j + 1 < N) :
    biasedE θ n (fun q => sgn (q ⟨j, by omega⟩) * sgn (q ⟨j + 1, hjN⟩)) -
      biasedE θ n (fun q => sgn (q ⟨j, by omega⟩)) * biasedE θ n (fun q => sgn (q ⟨j + 1, hjN⟩)) =
      sgn (n ⟨j, by omega⟩) * sgn (n ⟨j + 1, hjN⟩) *
        (Real.cos (2 * θ) ^ 2 * (1 - Real.cos (2 * θ) ^ 2)) := by
  rw [LawRemarksAux.biasedE_sgn_pair_succ θ n ⟨j - 1, by omega⟩ ⟨j, by omega⟩ ⟨j + 1, hjN⟩
      (by show j - 1 + 1 = j; omega) rfl,
    LawRemarksAux.biasedE_sgn_succ θ n ⟨j - 1, by omega⟩ ⟨j, by omega⟩
      (by show j - 1 + 1 = j; omega),
    LawRemarksAux.biasedE_sgn_succ θ n ⟨j, by omega⟩ ⟨j + 1, hjN⟩ rfl]
  ring

namespace LawRemarksAux

/-- The links `b_{j-1}, b_j` adjacent to the site `j` (only `b_0` for `j = 0`). -/
def nbrLinks (j : Fin N) : Finset (Fin N) := Finset.univ.filter (fun i => i.val + 1 = j.val ∨ i = j)

/-- Membership in `nbrLinks`. -/
lemma mem_nbrLinks {i j : Fin N} : i ∈ nbrLinks j ↔ i.val + 1 = j.val ∨ i = j := by
  simp [nbrLinks]

/-- `(-1)^{q_j} = (-1)^{n_j} ∏_{i ∈ nbrLinks j} (-1)^{x_i}` for `q = chargeOf n x`. -/
lemma sgn_chargeOf (n x : Cfg N) (j : Fin N) :
    sgn (chargeOf n x j) = sgn (n j) * ∏ i ∈ nbrLinks j, sgn (x i) := by
  by_cases hj : j.val = 0
  · have hS : nbrLinks j = {j} := by
      ext i
      simp only [mem_nbrLinks, Finset.mem_singleton, Fin.ext_iff]
      omega
    rw [hS, Finset.prod_singleton, sgn_chargeOf_zero n x j hj]
  · obtain ⟨i, hi⟩ : ∃ i : Fin N, i.val + 1 = j.val :=
      ⟨⟨j.val - 1, by have := j.isLt; omega⟩, by show j.val - 1 + 1 = j.val; omega⟩
    have hij : i ≠ j := fun e => by rw [e] at hi; omega
    have hS : nbrLinks j = {i, j} := by
      ext k
      simp only [mem_nbrLinks, Finset.mem_insert, Finset.mem_singleton, Fin.ext_iff]
      omega
    rw [hS, Finset.prod_pair hij, sgn_chargeOf_succ n x i j hi]

/-- Sites at distance at least two have disjoint sets of adjacent links. -/
lemma disjoint_nbrLinks (j k : Fin N) (h : j.val + 2 ≤ k.val) :
    Disjoint (nbrLinks j) (nbrLinks k) := by
  rw [Finset.disjoint_left]
  intro i hj hk
  rw [mem_nbrLinks, Fin.ext_iff] at hj hk
  omega

/-- Every function on `ZMod 2` is affine in the sign. -/
lemma affine_zmod2 (f : ZMod 2 → ℝ) (a : ZMod 2) :
    f a = (f 0 + f 1) / 2 + (f 0 - f 1) / 2 * sgn a := by
  rcases FreeDynAux.zmod2_cases a with rfl | rfl
  · rw [FreeDynAux.sgn_zero]; ring
  · rw [FreeDynAux.sgn_one]; ring

/-- Expectation of `f(q_j)` for `f` affine in the sign. -/
lemma biasedE_affine (θ : ℝ) (n : Cfg N) (j : Fin N) (f : ZMod 2 → ℝ) (α β : ℝ)
    (hf : ∀ a, f a = α + β * sgn a) :
    biasedE θ n (fun q => f (q j)) =
      α + β * (sgn (n j) * Real.cos (2 * θ) ^ (nbrLinks j).card) := by
  have h0 := sum_linkWt_mul_prod_sgn θ (∅ : Finset (Fin N))
  have h1 := sum_linkWt_mul_prod_sgn θ (nbrLinks j)
  simp only [Finset.prod_empty, mul_one, Finset.card_empty, pow_zero] at h0
  rw [biasedE_eq_sum]
  calc ∑ x : Cfg N, (∏ i, linkWt θ (x i)) * f (chargeOf n x j)
      = ∑ x : Cfg N, (α * (∏ i, linkWt θ (x i)) +
          β * sgn (n j) * ((∏ i, linkWt θ (x i)) * ∏ i ∈ nbrLinks j, sgn (x i))) :=
        Finset.sum_congr rfl fun x _ => by rw [hf, sgn_chargeOf]; ring
    _ = α * ∑ x : Cfg N, (∏ i, linkWt θ (x i)) +
          β * sgn (n j) * ∑ x : Cfg N, ((∏ i, linkWt θ (x i)) * ∏ i ∈ nbrLinks j, sgn (x i)) := by
        rw [Finset.sum_add_distrib, Finset.mul_sum, Finset.mul_sum]
    _ = _ := by rw [h0, h1]; ring

/-- Expectation of `f(q_j) g(q_k)` for affine `f, g` and disjoint sets of adjacent links. -/
lemma biasedE_affine₂ (θ : ℝ) (n : Cfg N) (j k : Fin N) (hd : Disjoint (nbrLinks j) (nbrLinks k))
    (f g : ZMod 2 → ℝ) (α β γ δ : ℝ) (hf : ∀ a, f a = α + β * sgn a)
    (hg : ∀ a, g a = γ + δ * sgn a) :
    biasedE θ n (fun q => f (q j) * g (q k)) =
      (α + β * (sgn (n j) * Real.cos (2 * θ) ^ (nbrLinks j).card)) *
        (γ + δ * (sgn (n k) * Real.cos (2 * θ) ^ (nbrLinks k).card)) := by
  have h0 := sum_linkWt_mul_prod_sgn θ (∅ : Finset (Fin N))
  have hj := sum_linkWt_mul_prod_sgn θ (nbrLinks j)
  have hk := sum_linkWt_mul_prod_sgn θ (nbrLinks k)
  have hjk := sum_linkWt_mul_prod_sgn θ (nbrLinks j ∪ nbrLinks k)
  simp only [Finset.prod_empty, mul_one, Finset.card_empty, pow_zero] at h0
  rw [Finset.card_union_of_disjoint hd, pow_add] at hjk
  rw [biasedE_eq_sum]
  calc ∑ x : Cfg N, (∏ i, linkWt θ (x i)) * (f (chargeOf n x j) * g (chargeOf n x k))
      = ∑ x : Cfg N, (α * γ * (∏ i, linkWt θ (x i)) +
          α * δ * sgn (n k) * ((∏ i, linkWt θ (x i)) * ∏ i ∈ nbrLinks k, sgn (x i)) +
          β * γ * sgn (n j) * ((∏ i, linkWt θ (x i)) * ∏ i ∈ nbrLinks j, sgn (x i)) +
          β * δ * (sgn (n j) * sgn (n k)) *
            ((∏ i, linkWt θ (x i)) * ∏ i ∈ nbrLinks j ∪ nbrLinks k, sgn (x i))) :=
        Finset.sum_congr rfl fun x _ => by
          rw [hf, hg, sgn_chargeOf, sgn_chargeOf, Finset.prod_union hd]; ring
    _ = α * γ * ∑ x : Cfg N, (∏ i, linkWt θ (x i)) +
          α * δ * sgn (n k) *
            ∑ x : Cfg N, ((∏ i, linkWt θ (x i)) * ∏ i ∈ nbrLinks k, sgn (x i)) +
          β * γ * sgn (n j) *
            ∑ x : Cfg N, ((∏ i, linkWt θ (x i)) * ∏ i ∈ nbrLinks j, sgn (x i)) +
          β * δ * (sgn (n j) * sgn (n k)) *
            ∑ x : Cfg N, ((∏ i, linkWt θ (x i)) * ∏ i ∈ nbrLinks j ∪ nbrLinks k, sgn (x i)) := by
        simp only [Finset.sum_add_distrib, Finset.mul_sum]
    _ = _ := by rw [h0, hj, hk, hjk]; ring

end LawRemarksAux

/-- `q_j` and `q_k` are independent for `k ≥ j + 2` (paper: `rem:productlaws`(b)). -/
theorem biased_indep_far (θ : ℝ) (n : Cfg N) (j k : Fin N) (hjk : j.val + 2 ≤ k.val)
    (f g : ZMod 2 → ℝ) :
    biasedE θ n (fun q => f (q j) * g (q k)) =
      biasedE θ n (fun q => f (q j)) * biasedE θ n (fun q => g (q k)) := by
  rw [LawRemarksAux.biasedE_affine₂ θ n j k (LawRemarksAux.disjoint_nbrLinks j k hjk) f g _ _ _ _
      (LawRemarksAux.affine_zmod2 f) (LawRemarksAux.affine_zmod2 g),
    LawRemarksAux.biasedE_affine θ n j f _ _ (LawRemarksAux.affine_zmod2 f),
    LawRemarksAux.biasedE_affine θ n k g _ _ (LawRemarksAux.affine_zmod2 g)]

/-- The density with respect to the uniform law is at most `(2 max(ϱ, 1 - ϱ))^N` (paper: `rem:biased`). -/
theorem biasedLaw_le (θ : ℝ) (n q : Cfg N) :
    biasedLaw θ n q ≤ max (Real.cos θ ^ 2) (Real.sin θ ^ 2) ^ N := by
  rw [biasedLaw_eq]
  calc ∏ j, (if linkOf n q j = 0 then Real.cos θ ^ 2 else Real.sin θ ^ 2)
      ≤ ∏ _j : Fin N, max (Real.cos θ ^ 2) (Real.sin θ ^ 2) :=
        Finset.prod_le_prod₀ (fun j _ => by split_ifs <;> positivity)
          (fun j _ => by split_ifs; exacts [le_max_left _ _, le_max_right _ _])
    _ = _ := by rw [Finset.prod_const, Finset.card_univ, Fintype.card_fin]

/-- The maximal density is attained (paper: `rem:biased`). -/
theorem exists_biasedLaw_eq (θ : ℝ) (n : Cfg N) :
    ∃ q, biasedLaw θ n q = max (Real.cos θ ^ 2) (Real.sin θ ^ 2) ^ N := by
  by_cases h : Real.sin θ ^ 2 ≤ Real.cos θ ^ 2
  · refine ⟨chargeOf n 0, ?_⟩
    rw [biasedLaw_eq, linkOf_chargeOf, max_eq_left h]
    simp
  · refine ⟨chargeOf n (fun _ => 1), ?_⟩
    rw [biasedLaw_eq, linkOf_chargeOf, max_eq_right (le_of_not_ge h)]
    simp

/-- Full support for `0 < ϱ < 1` (paper: `rem:biased`). -/
theorem biasedLaw_pos (θ : ℝ) (h1 : Real.cos θ ≠ 0) (h2 : Real.sin θ ≠ 0) (n q : Cfg N) :
    0 < biasedLaw θ n q := by
  rw [biasedLaw_eq]
  exact Finset.prod_pos fun j _ => by split_ifs <;> positivity

/-- For `ϱ = 1/2` the law is i.i.d. uniform (paper: `rem:productlaws`(b)). -/
theorem biasedLaw_half (θ : ℝ) (hθ : Real.cos θ ^ 2 = 1 / 2) (n q : Cfg N) :
    biasedLaw θ n q = (1 / 2) ^ N := by
  have hs : Real.sin θ ^ 2 = 1 / 2 := by have := Real.sin_sq_add_cos_sq θ; linarith
  rw [biasedLaw_eq]
  simp only [hθ, hs, ite_self, Finset.prod_const, Finset.card_univ, Fintype.card_fin]

/-- For `ϱ = 1` the law is a point mass: a single sector (paper: `rem:productlaws`(b)). -/
theorem biasedLaw_pure (θ : ℝ) (hθ : Real.sin θ = 0) (n q : Cfg N) :
    biasedLaw θ n q = if q = chargeOf n 0 then 1 else 0 := by
  have hc : Real.cos θ ^ 2 = 1 := by
    have := Real.sin_sq_add_cos_sq θ; rw [hθ] at this; linarith
  have hs : Real.sin θ ^ 2 = 0 := by rw [hθ]; ring
  have key : (∀ j, linkOf n q j = 0) ↔ q = chargeOf n 0 := by
    constructor
    · intro h
      rw [← chargeOf_linkOf n q]
      congr 1
      funext j
      exact h j
    · rintro rfl j
      rw [linkOf_chargeOf]
      rfl
  rw [biasedLaw_eq]
  simp only [hc, hs]
  by_cases hq : q = chargeOf n 0
  · rw [ite_eq_left hq]
    exact Finset.prod_eq_one fun j _ => ite_eq_left (key.mpr hq j)
  · rw [ite_eq_right hq]
    obtain ⟨j, hj⟩ : ∃ j, linkOf n q j ≠ 0 := by
      by_contra hcon
      exact hq (key.mp fun j => not_not.mp fun hj => hcon ⟨j, hj⟩)
    exact Finset.prod_eq_zero (Finset.mem_univ j) (ite_eq_right hj)

/-! ### Transfer of upper bounds (paper: `rem:productlaws`(c)) -/

/-- An upper bound on a uniform average of a non-negative sector functional transfers to every law with
density at most `K` with respect to the uniform law. -/
theorem upper_bound_transfer (μ : Cfg N → ℝ) (K : ℝ) (hμ0 : ∀ q, 0 ≤ μ q)
    (hμK : ∀ q, μ q ≤ K * (1 / 2) ^ N) (F : Cfg N → ℝ) (hF : ∀ q, 0 ≤ F q) :
    ∑ q, μ q * F q ≤ K * avgQ F := by
  unfold avgQ
  calc ∑ q, μ q * F q ≤ ∑ q, K * (1 / 2) ^ N * F q :=
        Finset.sum_le_sum fun q _ => mul_le_mul_of_nonneg_right (hμK q) (hF q)
    _ = K * (1 / 2 ^ N * ∑ q, F q) := by
        rw [Finset.mul_sum, Finset.mul_sum]
        exact Finset.sum_congr rfl fun q _ => by rw [one_div_pow]; ring

/-! ### The dimer law on `ℤ` (paper: `rem:productlaws`(c)) -/

/-- `muDimZ` is a probability measure. -/
theorem isProbabilityMeasure_muDimZ : IsProbabilityMeasure muDimZ := by
  unfold muDimZ
  infer_instance

namespace LawRemarksAux

/-- The map `η ↦ (y ↦ η ⌊y/2⌋)` defining the dimer law is measurable. -/
lemma measurable_dimerMap : Measurable (fun η : CfgZ => fun y : ℤ => η (y / 2)) :=
  Measurable.of_eval fun y => measurable_pi_apply (y / 2)

/-- The dimer event `{q | ∀ k, q (2k) = q (2k+1)}` is measurable. -/
lemma measurableSet_dimerEvent :
    MeasurableSet {q : CfgZ | ∀ k : ℤ, q (2 * k) = q (2 * k + 1)} := by
  have : {q : CfgZ | ∀ k : ℤ, q (2 * k) = q (2 * k + 1)} =
      ⋂ k : ℤ, {q : CfgZ | q (2 * k) = q (2 * k + 1)} := by
    ext q; simp
  rw [this]
  exact MeasurableSet.iInter fun k =>
    measurableSet_eq_fun (measurable_pi_apply _) (measurable_pi_apply _)

/-- The index map `(k, b) ↦ 2k + b` of the first `K` pairs of sites. -/
def pairIdx (K : ℕ) (p : Fin K × Fin 2) : ℤ := 2 * (p.1.val : ℤ) + (p.2.val : ℤ)

/-- `pairIdx` is injective. -/
lemma pairIdx_injective (K : ℕ) : Function.Injective (pairIdx K) := by
  rintro ⟨⟨k, hk⟩, ⟨b, hb⟩⟩ ⟨⟨k', hk'⟩, ⟨b', hb'⟩⟩ h
  simp only [pairIdx] at h
  have h1 : k = k' := by omega
  have h2 : b = b' := by omega
  subst h1 h2
  rfl

/-- The diagonal `{v | v 0 = v 1}` has probability `1/2` under the uniform law on `(ZMod 2)²`. -/
lemma infinitePi_diag :
    Measure.infinitePi (fun _ : Fin 2 => unifZ2) {v : Fin 2 → ZMod 2 | v 0 = v 1} = 2⁻¹ := by
  have hset : {v : Fin 2 → ZMod 2 | v 0 = v 1} = {fun _ => 0} ∪ {fun _ => 1} := by
    ext v
    simp only [Set.mem_ofPred_eq, Set.mem_union, Set.mem_singleton_iff, funext_iff,
      Fin.forall_fin_two]
    rcases FreeDynAux.zmod2_cases (v 0) with h0 | h0 <;>
      rcases FreeDynAux.zmod2_cases (v 1) with h1 | h1 <;> simp [h0, h1]
  have hdisj : Disjoint ({fun _ => 0} : Set (Fin 2 → ZMod 2)) {fun _ => 1} := by
    rw [Set.disjoint_singleton]
    intro h
    exact absurd (congrFun h 0) (by decide)
  rw [hset, measure_union hdisj (measurableSet_singleton _),
    Measure.infinitePi_singleton_of_fintype, Measure.infinitePi_singleton_of_fintype]
  simp only [InfVolAux.unifZ2_singleton, Finset.prod_const, Finset.card_univ, Fintype.card_fin]
  rw [pow_two, ← add_mul, ENNReal.inv_two_add_inv_two, one_mul]

/-- Under the i.i.d. law the first `K` pairs `(2k, 2k+1)` all carry equal charges with
probability `2^{-K}`. -/
lemma muZ_pairs (K : ℕ) :
    muZ {q : CfgZ | ∀ k : Fin K, q (2 * (k.val : ℤ)) = q (2 * (k.val : ℤ) + 1)} = 2⁻¹ ^ K := by
  let G : CfgZ → (Fin K → Fin 2 → ZMod 2) := fun q k b => q (pairIdx K (k, b))
  have hG : Measurable G :=
    Measurable.of_eval fun k => Measurable.of_eval fun b => measurable_pi_apply _
  have hR : Measurable (fun (ω : CfgZ) (p : Fin K × Fin 2) => ω (pairIdx K p)) :=
    Measurable.of_eval fun p => measurable_pi_apply _
  have hmap : muZ.map G =
      Measure.infinitePi (fun _ : Fin K => Measure.infinitePi (fun _ : Fin 2 => unifZ2)) := by
    calc muZ.map G = (muZ.map (fun (ω : CfgZ) (p : Fin K × Fin 2) => ω (pairIdx K p))).map
          (MeasurableEquiv.curry (Fin K) (Fin 2) (ZMod 2)) := by
          rw [Measure.map_map (MeasurableEquiv.measurable _) hR]; rfl
      _ = (Measure.infinitePi (fun _ : Fin K × Fin 2 => unifZ2)).map
          (MeasurableEquiv.curry (Fin K) (Fin 2) (ZMod 2)) := by
          rw [muZ, Measure.map_infinitePi_infinitePi_of_inj (pairIdx_injective K)]
      _ = _ := Measure.infinitePi_map_curry (fun (_ : Fin K) (_ : Fin 2) => unifZ2)
  have hset : {q : CfgZ | ∀ k : Fin K, q (2 * (k.val : ℤ)) = q (2 * (k.val : ℤ) + 1)} =
      G ⁻¹' Set.pi Set.univ (fun _ => {v : Fin 2 → ZMod 2 | v 0 = v 1}) := by
    ext q
    simp [G, pairIdx]
  rw [hset, ← Measure.map_apply hG MeasurableSet.of_discrete, hmap, Measure.infinitePi_eq_pi,
    Measure.pi_pi]
  simp only [infinitePi_diag, Finset.prod_const, Finset.card_univ, Fintype.card_fin]

end LawRemarksAux

/-- Under the i.i.d. law, almost surely not all pairs `(2k, 2k+1)` carry equal charges. -/
theorem muZ_dimerEvent : muZ {q : CfgZ | ∀ k : ℤ, q (2 * k) = q (2 * k + 1)} = 0 := by
  have hle : ∀ K : ℕ, muZ {q : CfgZ | ∀ k : ℤ, q (2 * k) = q (2 * k + 1)} ≤ 2⁻¹ ^ K :=
    fun K => by
      rw [← LawRemarksAux.muZ_pairs K]
      exact measure_mono fun q hq k => hq _
  exact le_antisymm
    (ge_of_tendsto' (ENNReal.tendsto_pow_atTop_nhds_zero_of_lt_one (by norm_num)) hle)
    zero_le

/-- Under the dimer law, almost surely all pairs `(2k, 2k+1)` carry equal charges. -/
theorem muDimZ_dimerEvent_compl :
    muDimZ {q : CfgZ | ∀ k : ℤ, q (2 * k) = q (2 * k + 1)}ᶜ = 0 := by
  unfold muDimZ
  rw [Measure.map_apply LawRemarksAux.measurable_dimerMap
    LawRemarksAux.measurableSet_dimerEvent.compl]
  have hempty : (fun η : CfgZ => fun y : ℤ => η (y / 2)) ⁻¹'
      {q : CfgZ | ∀ k : ℤ, q (2 * k) = q (2 * k + 1)}ᶜ = ∅ := by
    ext η
    simp only [Set.mem_preimage, Set.mem_compl_iff, Set.mem_ofPred_eq, Set.mem_empty_iff_false,
      iff_false, not_not]
    intro k
    rw [show 2 * k / 2 = k by omega, show (2 * k + 1) / 2 = k by omega]
  rw [hempty, measure_empty]

/-- **In infinite volume the i.i.d. law and the dimer law are mutually singular** (paper:
`rem:productlaws`(c)). -/
theorem muZ_mutuallySingular_muDimZ : muZ ⟂ₘ muDimZ :=
  ⟨_, LawRemarksAux.measurableSet_dimerEvent, muZ_dimerEvent, muDimZ_dimerEvent_compl⟩

end DFL
#print axioms DFL.biasedLaw_eq
#print axioms DFL.sum_biasedLaw
#print axioms DFL.biasedE_first
#print axioms DFL.biasedE_bulk
#print axioms DFL.biasedE_pair01
#print axioms DFL.biased_cov01
#print axioms DFL.biased_cov_bulk
#print axioms DFL.biased_indep_far
#print axioms DFL.biasedLaw_le
#print axioms DFL.exists_biasedLaw_eq
#print axioms DFL.biasedLaw_pos
#print axioms DFL.biasedLaw_half
#print axioms DFL.biasedLaw_pure
#print axioms DFL.upper_bound_transfer
#print axioms DFL.isProbabilityMeasure_muDimZ
#print axioms DFL.muZ_dimerEvent
#print axioms DFL.muDimZ_dimerEvent_compl
#print axioms DFL.muZ_mutuallySingular_muDimZ
