import DFL.Sector

/-!
# Sector laws of product initial states (paper: Lemma `lem:productlaw`)

For a product state `χ = ⊗_j χ_j` of the links (given by its `σˣ` amplitudes `χ_j(x_j)`) and any
matter state `ψ`, the dual state is `∑_q ψ^{(q)} ⊗ |q⟩` with `ψ^{(q)}(n) = χ(x(n,q)) ψ(n)`.
Consequences:
* (b) Fock matter `|n⟩`: the sector law is `μ(q) = |χ(x(n,q))|²`, the image of the product law of
  the link labels under the bijection `x ↦ q`;
* (c) `σᶻ`-product gauge states `|s⟩_z`: the charges are i.i.d. Bernoulli(1/2), with matter state
  `W_D ψ` (a sign `(-1)^{N_{≤j}}` for every link `j` in the state `|↓⟩_z`);
* (d) `σˣ`-product gauge states with Fock matter: a single, definite-charge sector.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- A product state of the links, given by the `σˣ` amplitudes `χ j (x_j)` of each link. -/
def prodLink (χ : Fin N → ZMod 2 → ℂ) (x : Cfg N) : ℂ := ∏ j, χ j (x j)

/-- The product initial state `χ ⊗ ψ` (links in the product state `χ`, matter in `ψ`). -/
def prodInit (χ : Fin N → ZMod 2 → ℂ) (ψ : Cfg N → ℂ) : Cfg N × Cfg N → ℂ :=
  fun p => ψ p.1 * prodLink χ p.2

/-- (a) The dual form of a product initial state. -/
theorem prodInit_dual (χ : Fin N → ZMod 2 → ℂ) (ψ : Cfg N → ℂ) (n q : Cfg N) :
    prodInit χ ψ (dualEquiv (n, q)) = prodLink χ (linkOf n q) * ψ n := by
  simp only [prodInit, dualEquiv_apply]
  exact mul_comm _ _

/-- (a) Sector decomposition of the dynamics for a product initial state: every gauge-invariant
expectation is a sum over sectors of free-fermion expectations with initial matter states
`ψ^{(q)} = χ(x(·, q)) ψ`. -/
theorem sector_sum_prod (J h : ℝ) (χ : Fin N → ZMod 2 → ℂ) (ψ : Cfg N → ℂ)
    (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (H J h) t *ᵥ prodInit χ ψ) =
      ∑ q : Cfg N, expect (sectorBlock O q)
        (evo (Hsec J h q) t *ᵥ fun n => prodLink χ (linkOf n q) * ψ n) := by
  rw [sec_expect_evo_dual J h _ O hO t]
  simp only [prodInit_dual]

/-- (b) Fock matter: the sector law is `μ(q) = |χ(x(n,q))|²`. -/
theorem sector_average_prod_fock (J h : ℝ) (χ : Fin N → ZMod 2 → ℂ) (n : Cfg N)
    (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (H J h) t *ᵥ prodInit χ (Pi.single n 1)) =
      ∑ q : Cfg N, (Complex.normSq (prodLink χ (linkOf n q)) : ℂ) *
        expect (sectorBlock O q) (evo (Hsec J h q) t *ᵥ Pi.single n 1) := by
  rw [sector_sum_prod J h χ _ O hO t]
  refine Finset.sum_congr rfl fun q _ => ?_
  have hvec : (fun n' => prodLink χ (linkOf n' q) * (Pi.single n 1 : Cfg N → ℂ) n') =
      prodLink χ (linkOf n q) • (Pi.single n 1 : Cfg N → ℂ) := by
    ext n'
    by_cases hn : n' = n
    · subst hn; simp
    · simp [hn]
  rw [hvec, mulVec_smul, sec_expect_smul, Complex.normSq_eq_conj_mul_self]
  rfl

/-- (b) The sector law is the image of the product law of the link labels under `x ↦ q`. -/
theorem sum_sectorLaw_prod (χ : Fin N → ZMod 2 → ℂ) (n : Cfg N) (f : Cfg N → ℂ) :
    ∑ q : Cfg N, (Complex.normSq (prodLink χ (linkOf n q)) : ℂ) * f q =
      ∑ x : Cfg N, (Complex.normSq (prodLink χ x) : ℂ) * f (chargeOf n x) := by
  let e : Cfg N ≃ Cfg N :=
    ⟨linkOf n, chargeOf n, fun q => chargeOf_linkOf n q, fun x => linkOf_chargeOf n x⟩
  exact Fintype.sum_equiv e _ _ (fun q => by simp only [e, Equiv.coe_fn_mk, chargeOf_linkOf])

/-- `σˣ` amplitudes of the `σᶻ` eigenstates: `⟨x|↑⟩_z = 2^{-1/2}`, `⟨x|↓⟩_z = 2^{-1/2} (-1)^x`
(label `s = 0` for `|↑⟩_z`, `s = 1` for `|↓⟩_z`). -/
def zAmp (s x : ZMod 2) : ℂ := ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ * (if s = 1 ∧ x = 1 then -1 else 1)

/-- The matter sign `W_D = ∏_{j : s_j = 1} (-1)^{N_{≤ j}}`. -/
def signWD (s : Cfg N) (n : Cfg N) : ℂ :=
  ∏ j, if s j = 1 then ((sgn (∑ k ∈ univ.filter (fun k => k ≤ j), n k) : ℝ) : ℂ) else 1

/-! ### Auxiliary lemmas for `σᶻ`-product gauge states -/

/-- `sgn` turns addition in `ZMod 2` into multiplication. -/
private lemma pl_sgn_add (a b : ZMod 2) : sgn (a + b) = sgn a * sgn b := by
  rcases sec_zmod2_cases a with rfl | rfl <;> rcases sec_zmod2_cases b with rfl | rfl <;>
    simp [sgn, show (1 : ZMod 2) + 1 = 0 from by decide, show (1 : ZMod 2) ≠ 0 from by decide]

/-- `sgn a ^ 2 = 1`. -/
private lemma pl_sgn_mul_self (a : ZMod 2) : sgn a * sgn a = 1 := by
  rcases sec_zmod2_cases a with rfl | rfl <;>
    simp [sgn, show (1 : ZMod 2) ≠ 0 from by decide]

/-- The amplitudes `zAmp` in terms of `sgn`. -/
private lemma zAmp_eq_sgn (s x : ZMod 2) :
    zAmp s x = ((Real.sqrt 2 : ℝ) : ℂ)⁻¹ * (if s = 1 then ((sgn x : ℝ) : ℂ) else 1) := by
  unfold zAmp
  rcases sec_zmod2_cases s with rfl | rfl <;> rcases sec_zmod2_cases x with rfl | rfl <;>
    simp [sgn, show (1 : ZMod 2) ≠ 0 from by decide, show (0 : ZMod 2) ≠ 1 from by decide]

/-- `signWD` squares to one. -/
private lemma signWD_mul_self (s n : Cfg N) : signWD s n * signWD s n = 1 := by
  unfold signWD
  rw [← Finset.prod_mul_distrib]
  refine Finset.prod_eq_one fun j _ => ?_
  split_ifs
  · rw [← Complex.ofReal_mul, pl_sgn_mul_self, Complex.ofReal_one]
  · exact one_mul 1

/-- `signWD` is real. -/
private lemma star_signWD (s n : Cfg N) : star (signWD s n) = signWD s n := by
  unfold signWD
  rw [star_prod]
  refine Finset.prod_congr rfl fun j _ => ?_
  split_ifs
  · exact Complex.conj_ofReal _
  · exact star_one ℂ

/-- In the dual frame, the `σᶻ`-product amplitude factorizes into the constant `2^{-N/2}`, a charge
sign `ε(q) = signWD s q` and the matter sign `W_D = signWD s n`. -/
private lemma prodLink_zAmp_linkOf (s n q : Cfg N) :
    prodLink (fun j => zAmp (s j)) (linkOf n q) =
      (((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ * (signWD s q * signWD s n) := by
  unfold prodLink signWD
  simp only [zAmp_eq_sgn]
  rw [Finset.prod_mul_distrib, Finset.prod_const, Finset.card_univ, Fintype.card_fin,
    ← Finset.prod_mul_distrib]
  congr 1
  · push_cast
    rw [inv_pow]
  · refine Finset.prod_congr rfl fun j _ => ?_
    split_ifs
    · simp only [linkOf, Finset.sum_add_distrib, pl_sgn_add, Complex.ofReal_mul]
    · exact (one_mul 1).symm

/-- (c) `σᶻ`-product gauge states: the charges are i.i.d. Bernoulli(1/2) and every sector carries
the matter state `W_D ψ`. -/
theorem sector_average_zprod (J h : ℝ) (s : Cfg N) (ψ : Cfg N → ℂ)
    (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (H J h) t *ᵥ prodInit (fun j => zAmp (s j)) ψ) =
      (1 / 2 ^ N : ℂ) * ∑ q : Cfg N,
        expect (sectorBlock O q) (evo (Hsec J h q) t *ᵥ fun n => signWD s n * ψ n) := by
  rw [sector_sum_prod J h _ ψ O hO t, Finset.mul_sum]
  refine Finset.sum_congr rfl fun q _ => ?_
  have hvec : (fun n => prodLink (fun j => zAmp (s j)) (linkOf n q) * ψ n) =
      ((((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ * signWD s q) • (fun n => signWD s n * ψ n) := by
    ext n
    rw [prodLink_zAmp_linkOf, Pi.smul_apply, smul_eq_mul]
    ring
  have hz : star ((((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ * signWD s q) *
      ((((Real.sqrt 2) ^ N : ℝ) : ℂ)⁻¹ * signWD s q) = (1 / 2 ^ N : ℂ) := by
    rw [star_mul', star_signWD, mul_mul_mul_comm, sec_norm_sq_prefactor, signWD_mul_self,
      mul_one]
  rw [hvec, mulVec_smul, sec_expect_smul, hz]

/-- The product of the `σˣ`-eigenstate amplitudes is the indicator of `x = ξ`. -/
private lemma prodLink_indicator (ξ x : Cfg N) :
    prodLink (fun j x => if x = ξ j then (1 : ℂ) else 0) x = if x = ξ then 1 else 0 := by
  unfold prodLink
  rw [Fintype.prod_boole]
  by_cases hx : x = ξ
  · subst hx; simp
  · rw [ite_eq_right hx, ite_eq_right]
    exact fun h => hx (funext h)

/-- (d) `σˣ`-product gauge states `|ξ⟩_σ` with Fock matter `|n⟩` are definite-charge states:
the dynamics stays in the single sector `q = chargeOf n ξ`. -/
theorem expect_xprod_fock (J h : ℝ) (ξ n : Cfg N)
    (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (H J h) t *ᵥ prodInit (fun j x => if x = ξ j then 1 else 0) (Pi.single n 1)) =
      expect (sectorBlock O (chargeOf n ξ))
        (evo (Hsec J h (chargeOf n ξ)) t *ᵥ Pi.single n 1) := by
  rw [sector_sum_prod J h _ _ O hO t, Finset.sum_eq_single (chargeOf n ξ)]
  · have hvec : (fun n' => prodLink (fun j x => if x = ξ j then (1 : ℂ) else 0)
        (linkOf n' (chargeOf n ξ)) * (Pi.single n 1 : Cfg N → ℂ) n') =
        (Pi.single n 1 : Cfg N → ℂ) := by
      ext n'
      rw [prodLink_indicator]
      by_cases hn : n' = n
      · subst hn; simp [linkOf_chargeOf]
      · simp [hn]
    rw [hvec]
  · intro q _ hq
    have hvec : (fun n' => prodLink (fun j x => if x = ξ j then (1 : ℂ) else 0) (linkOf n' q) *
        (Pi.single n 1 : Cfg N → ℂ) n') = 0 := by
      ext n'
      rw [prodLink_indicator]
      by_cases hn : n' = n
      · subst hn
        have hne : linkOf n' q ≠ ξ := fun h' => hq (by rw [← h', chargeOf_linkOf])
        simp [hne]
      · simp [hn]
    rw [hvec, mulVec_zero]
    simp [expect]
  · intro h'
    exact absurd (Finset.mem_univ _) h'

end DFL
