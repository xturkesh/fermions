import DFL.TwoD

/-!
# Two dimensions: flux sectors (paper: Lemma `lem:2dgauge`(i), Proposition `prop:2dlaws`(a),(b))

The ℤ₂ gauge theory `FinGraph.ham` on an arbitrary finite graph, written in the `σˣ` basis of the
edges (`DFL/Graph.lean`). Link signs `s ∈ {±1}^E` are encoded as `s : Cfg M` (`sgn (s e)`); the
corresponding one-body operator is `h_{q,s} = K_s + 2h diag(q)` (`FinGraph.hqs`, `DFL/TwoD.lean`).

* `flux_eq_iff_gauge` (paper: `lem:2dgauge`(i)): two sign configurations have the same fluxes
  (the same `±1` product around every cycle, i.e. the same values `∑_e s_e z_e` on every
  `z ∈ ker ∂`) iff they differ by a gauge transformation `s ↦ s + δθ`, `(δθ)_e = θ_{s(e)} + θ_{t(e)}`.
* `sector_average_flux` (paper: `prop:2dlaws`(a)): for the `σᶻ`-product gauge state `|s⟩_z` and Fock
  matter `|n⟩`, the flux is that of `s`, and the dynamics of every matter observable diagonal in the
  occupations is the uniform average over the reachable charge sectors `q ∈ n + im ∂` of free-fermion
  dynamics with `h_{q,s}` (for `s = +` this is `FinGraph.sector_average_matter`).
* `sector_average_xpol` (paper: `prop:2dlaws`(b)): for the `σˣ`-polarized gauge state `|+⋯+⟩` the
  charges are fixed (`q⁰ = n`, i.e. `G_i = (-1)^{n_i}`) and the dynamics is the uniform average over
  all sign configurations `s` of free-fermion dynamics with `h_{q⁰,s}`; by `xpol_summand_gauge` the
  summand depends on `s` only through its gauge class, so this is the uniform average over the flux
  configurations: random `π`-flux.
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

namespace FinGraph

variable {N M : ℕ} (G : FinGraph N M)

/-- The `σᶻ`-product gauge state `|s⟩_z` (edge `e` in the `σᶻ` eigenstate `sgn (s e)`) times the
matter state `ψ`, in the `σˣ` basis: amplitude `2^{-M/2} (-1)^{∑_e s_e x_e} ψ(n)`. -/
def initS (_G : FinGraph N M) (s : Cfg M) (ψ : Cfg N → ℂ) : Cfg N × Cfg M → ℂ :=
  fun p => (((Real.sqrt 2) ^ M : ℝ) : ℂ)⁻¹ * ((sgn (∑ e, s e * p.2 e) : ℝ) : ℂ) * ψ p.1

/-- The gauge transformation `θ` acting on link labels: `(s + δθ)_e = s_e + θ_{s(e)} + θ_{t(e)}`. -/
def gaugeAct (θ : Cfg N) (s : Cfg M) : Cfg M := fun e => s e + θ (G.src e) + θ (G.tgt e)

/-! ### Auxiliary lemmas: `ZMod 2` and signs -/

/-- Every element of `ZMod 2` is `0` or `1`. -/
private lemma z2_cases (a : ZMod 2) : a = 0 ∨ a = 1 := by
  revert a; decide

/-- `1 + 1 = 0` in `ZMod 2`. -/
private lemma z2_one_add_one : (1 : ZMod 2) + 1 = 0 := by decide

/-- `a + a = 0` in `ZMod 2`. -/
private lemma z2_add_self (a : ZMod 2) : a + a = 0 := by
  revert a; decide

/-- `sgn` turns addition in `ZMod 2` into multiplication. -/
private lemma sgn_add (a b : ZMod 2) : sgn (a + b) = sgn a * sgn b := by
  rcases z2_cases a with rfl | rfl <;> rcases z2_cases b with rfl | rfl <;>
    simp [sgn, z2_one_add_one]

/-- `sgn a ^ 2 = 1`. -/
private lemma sgn_mul_self (a : ZMod 2) : sgn a * sgn a = 1 := by
  rcases z2_cases a with rfl | rfl <;> norm_num [sgn]

/-- `sgn a = ±1`. -/
private lemma sgn_cases (a : ZMod 2) : sgn a = 1 ∨ sgn a = -1 := by
  rcases z2_cases a with rfl | rfl <;> simp [sgn]

/-- `sgn 0 = 1`. -/
private lemma sgn_zero : sgn 0 = 1 := by simp [sgn]

/-- `(-1)^a = 1 - 2a` for `a ∈ {0, 1}`. -/
private lemma sgn_eq_one_sub_two_occ (a : ZMod 2) : sgn a = 1 - 2 * occ a := by
  rcases z2_cases a with rfl | rfl <;> norm_num [sgn, occ]

/-! ### Auxiliary lemmas: the boundary map -/

/-- `∂` is additive. -/
private lemma bdry_add (x y : Cfg M) : G.bdry (x + y) = G.bdry x + G.bdry y := by
  funext i
  simp only [bdry, Pi.add_apply, Finset.sum_add_distrib]

/-- `∂ 0 = 0`. -/
private lemma bdry_zero : G.bdry 0 = 0 := by
  funext i
  simp [bdry]

/-- The boundary map as an additive homomorphism. -/
private def bdryHom : Cfg M →+ Cfg N where
  toFun := G.bdry
  map_zero' := bdry_zero G
  map_add' := bdry_add G

/-- `∂` commutes with subtraction. -/
private lemma bdry_sub (x y : Cfg M) : G.bdry (x - y) = G.bdry x - G.bdry y :=
  map_sub (bdryHom G) x y

/-- `∂` of a single edge. -/
private lemma bdry_single (e : Fin M) (i : Fin N) :
    G.bdry (Pi.single e 1) i = if G.src e = i ∨ G.tgt e = i then 1 else 0 := by
  unfold bdry
  rw [Finset.sum_pi_single']
  simp

/-- A sum over the two endpoints of an edge. -/
private lemma sum_ite_endpoint {R : Type*} [AddCommMonoid R] (e : Fin M) (f : Fin N → R) :
    ∑ i, (if G.src e = i ∨ G.tgt e = i then f i else 0) = f (G.src e) + f (G.tgt e) := by
  have h : ∀ i, (if G.src e = i ∨ G.tgt e = i then f i else 0) =
      (if G.src e = i then f i else 0) + (if G.tgt e = i then f i else 0) := by
    intro i
    by_cases h1 : G.src e = i
    · have h2 : ¬ G.tgt e = i := fun h2 => G.src_ne_tgt e (h1.trans h2.symm)
      rw [ite_eq_left (Or.inl h1), ite_eq_left h1, ite_eq_right h2, add_zero]
    · by_cases h2 : G.tgt e = i
      · rw [ite_eq_left (Or.inr h2), ite_eq_right h1, ite_eq_left h2, zero_add]
      · rw [ite_eq_right (not_or.2 ⟨h1, h2⟩), ite_eq_right h1, ite_eq_right h2, add_zero]
  rw [Finset.sum_congr rfl (fun i _ => h i), Finset.sum_add_distrib, Finset.sum_ite_eq,
    Finset.sum_ite_eq, ite_eq_left (Finset.mem_univ _), ite_eq_left (Finset.mem_univ _)]

/-- `∂` of a single edge is the sum of the two endpoints. -/
private lemma bdry_single_eq (e : Fin M) :
    G.bdry (Pi.single e 1) = Pi.single (G.src e) 1 + Pi.single (G.tgt e) 1 := by
  funext i
  rw [bdry_single, Pi.add_apply, Pi.single_apply, Pi.single_apply]
  by_cases h1 : G.src e = i
  · have h2 : ¬ G.tgt e = i := fun h2 => G.src_ne_tgt e (h1.trans h2.symm)
    rw [ite_eq_left (Or.inl h1), ite_eq_left h1.symm, ite_eq_right (Ne.symm h2), add_zero]
  · by_cases h2 : G.tgt e = i
    · rw [ite_eq_left (Or.inr h2), ite_eq_right (Ne.symm h1), ite_eq_left h2.symm, zero_add]
    · rw [ite_eq_right (not_or.2 ⟨h1, h2⟩), ite_eq_right (Ne.symm h1), ite_eq_right (Ne.symm h2), add_zero]

/-- Duality of `δ` and `∂`: `⟨δθ, z⟩ = ⟨θ, ∂z⟩`. -/
private lemma sum_coboundary_mul (θ : Cfg N) (z : Cfg M) :
    ∑ e, (θ (G.src e) + θ (G.tgt e)) * z e = ∑ i, θ i * G.bdry z i := by
  simp only [bdry, Finset.sum_filter, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun e _ => ?_)
  have h := G.sum_ite_endpoint e (fun i => θ i * z e)
  simp only [mul_ite, mul_zero]
  rw [h, add_mul]

/-- The boundary map as a `ZMod 2`-linear map. -/
private def bdryLin : Cfg M →ₗ[ZMod 2] Cfg N where
  toFun := G.bdry
  map_add' := G.bdry_add
  map_smul' := by
    intro c x
    funext i
    simp [bdry, Finset.mul_sum]

/-- **Fluxes are gauge classes** (paper: Lemma `lem:2dgauge`(i)): `s` and `s'` take the same values on
every cycle `z ∈ ker ∂` iff `s' = s + δθ` for some gauge transformation `θ`. -/
theorem flux_eq_iff_gauge (s s' : Cfg M) :
    (∀ z : Cfg M, G.bdry z = 0 → ∑ e, s e * z e = ∑ e, s' e * z e) ↔
      ∃ θ : Cfg N, s' = G.gaugeAct θ s := by
  constructor
  · intro hz
    -- the functional `z ↦ ⟨s' - s, z⟩` annihilates `ker ∂`, hence factors through `∂`
    let φ : Module.Dual (ZMod 2) (Cfg M) :=
      { toFun := fun z => ∑ e, (s' e - s e) * z e
        map_add' := by
          intro x y
          simp only [Pi.add_apply, mul_add, Finset.sum_add_distrib]
        map_smul' := by
          intro c x
          simp only [Pi.smul_apply, smul_eq_mul, RingHom.id_apply, Finset.mul_sum]
          refine Finset.sum_congr rfl (fun e _ => ?_)
          ring }
    have hφ : φ ∈ (LinearMap.ker G.bdryLin).dualAnnihilator := by
      rw [Submodule.mem_dualAnnihilator]
      intro z hz0
      have h0 : G.bdry z = 0 := hz0
      have := hz z h0
      change ∑ e, (s' e - s e) * z e = 0
      simp only [sub_mul, Finset.sum_sub_distrib, this, sub_self]
    rw [← LinearMap.range_dualMap_eq_dualAnnihilator_ker] at hφ
    obtain ⟨ψ, hψ⟩ := hφ
    refine ⟨fun i => ψ (Pi.single i 1), ?_⟩
    funext e
    have h1 : φ (Pi.single e 1) = s' e - s e := by
      change ∑ k, (s' k - s k) * (Pi.single e 1 : Cfg M) k = s' e - s e
      simp [Pi.single_apply]
    have h2 : φ (Pi.single e 1) = ψ (Pi.single (G.src e) 1) + ψ (Pi.single (G.tgt e) 1) := by
      rw [← hψ, LinearMap.dualMap_apply]
      change ψ (G.bdry (Pi.single e 1)) = _
      rw [bdry_single_eq, map_add]
    rw [h1] at h2
    simp only [gaugeAct]
    rw [add_assoc, ← h2]
    ring
  · rintro ⟨θ, rfl⟩ z hz
    have h := G.sum_coboundary_mul θ z
    rw [hz] at h
    simp only [Pi.zero_apply, mul_zero, Finset.sum_const_zero] at h
    simp only [gaugeAct, add_assoc, add_mul, Finset.sum_add_distrib] at h ⊢
    rw [h, add_zero]

/-- The `σᶻ`-product state with `s = 0` (all edges `|↑⟩_z`) is `initZ`. -/
theorem initS_zero (ψ : Cfg N → ℂ) : G.initS 0 ψ = G.initZ ψ := by
  funext p
  simp [initS, initZ, sgn]

/-! ### Auxiliary lemmas: fermionic matrix elements -/

/-- Matrix elements of `c_a† c_b`. -/
private lemma cre_mul_ann_apply (a b : Fin N) (n' n : Cfg N) :
    (cre a * ann b) n' n =
      if n b = 1 ∧ n' a = 1 ∧ Function.update n b 0 = Function.update n' a 0 then
        sgn (jwPar n' a) * sgn (jwPar n b) else 0 := by
  simp only [cre, Matrix.mul_apply, Matrix.transpose_apply, ann]
  rw [Finset.sum_eq_single (Function.update n b 0)]
  · by_cases h1 : n b = 1 <;> by_cases h2 : n' a = 1 <;>
      by_cases h3 : Function.update n b 0 = Function.update n' a 0 <;> simp [h1, h2, h3]
  · intro m _ hm
    simp [hm]
  · simp

/-- Two configurations with the same value at `x` that agree after emptying `x` are equal. -/
private lemma eq_of_update_eq {x : Fin N} {n n' : Cfg N} (hx : n' x = n x)
    (h : Function.update n x 0 = Function.update n' x 0) : n' = n := by
  funext k
  by_cases hk : k = x
  · subst hk; exact hx
  · have := congrFun h k
    simpa [Function.update_of_ne hk] using this.symm

/-- `c_x† c_x` is the diagonal occupation operator. -/
private lemma cre_mul_ann_self (x : Fin N) : cre x * ann x = diagonal (fun n => occ (n x)) := by
  ext n' n
  rw [cre_mul_ann_apply, diagonal_apply]
  by_cases hn : n x = 1
  · by_cases h' : n' x = 1 ∧ Function.update n x 0 = Function.update n' x 0
    · have : n' = n := eq_of_update_eq (h'.1.trans hn.symm) h'.2
      subst this
      have hs : sgn (jwPar n' x) * sgn (jwPar n' x) = 1 := by
        unfold sgn; split_ifs <;> norm_num
      simp [hn, hs, occ]
    · have : n' ≠ n := by rintro rfl; exact h' ⟨hn, rfl⟩
      simp only [hn, true_and, h', this, ite_false]
  · have h0 : n x = 0 := (z2_cases _).resolve_right hn
    by_cases h' : n' = n <;> simp [h', h0, occ]

/-- Support of `c_a† c_b` for `a ≠ b`: a nonzero entry `(n', n)` has `n' = n + δ_a + δ_b`. -/
private lemma cre_mul_ann_support {a b : Fin N} (hab : a ≠ b) {n' n : Cfg N}
    (h : (cre a * ann b) n' n ≠ 0) (i : Fin N) :
    n' i = n i + if a = i ∨ b = i then 1 else 0 := by
  rw [cre_mul_ann_apply] at h
  split_ifs at h with hc
  · obtain ⟨hb, ha, hu⟩ := hc
    have hua := congrFun hu a
    have hub := congrFun hu b
    rw [Function.update_of_ne hab, Function.update_self] at hua
    rw [Function.update_self, Function.update_of_ne hab.symm] at hub
    by_cases hia : a = i
    · subst hia
      simp [hua, ha]
    · by_cases hib : b = i
      · subst hib
        simp [← hub, hb, z2_one_add_one]
      · have hui := congrFun hu i
        rw [Function.update_of_ne (Ne.symm hib), Function.update_of_ne (Ne.symm hia)] at hui
        simp [hia, hib, hui]
  · exact absurd rfl h

/-- Support of the edge hopping: a nonzero entry `(n', n)` has `n' = n + ∂ δ_e`. -/
private lemma hopE_support (e : Fin M) {n' n : Cfg N} (h : G.hopE e n' n ≠ 0) :
    n' = n + G.bdry (Pi.single e 1) := by
  funext i
  rw [Pi.add_apply, bdry_single]
  rw [hopE, Matrix.add_apply] at h
  by_cases h1 : (cre (G.src e) * ann (G.tgt e)) n' n = 0
  · rw [h1, zero_add] at h
    rw [cre_mul_ann_support (G.src_ne_tgt e).symm h i]
    exact congrArg _ (if_congr or_comm rfl rfl)
  · exact cre_mul_ann_support (G.src_ne_tgt e) h1 i

/-! ### Auxiliary lemmas: the edge terms -/

/-- Flipping the label of `e` is adding `δ_e`. -/
private lemma update_add_one (x : Cfg M) (e : Fin M) :
    Function.update x e (x e + 1) = x + Pi.single e 1 := by
  funext k
  by_cases hk : k = e
  · subst hk; simp
  · simp [hk]

/-- Flipping twice is the identity. -/
private lemma add_single_add_single (x : Cfg M) (e : Fin M) :
    x + Pi.single e 1 + Pi.single e 1 = x := by
  rw [add_assoc, ← Pi.single_add, z2_one_add_one, Pi.single_zero, add_zero]

/-- Matrix elements of `σᶻ_e`. -/
private lemma szE_apply (e : Fin M) (x' x : Cfg M) :
    G.szE e x' x = if x = x' + Pi.single e 1 then 1 else 0 := by
  unfold szE
  rw [update_add_one]
  refine if_congr ⟨?_, ?_⟩ rfl rfl
  · rintro rfl; rw [add_single_add_single]
  · rintro rfl; rw [add_single_add_single]

/-- The star terms `1 ⊗ ∏ σˣ` are diagonal. -/
private lemma one_kron_star (i : Fin N) :
    (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ G.star i =
      diagonal (fun p : Cfg N × Cfg M => sgn (G.bdry p.2 i)) := by
  rw [FinGraph.star, ← diagonal_one, diagonal_kronecker_diagonal]
  simp

/-- `⟨s, x + δ_e⟩ = ⟨s, x⟩ + s_e`. -/
private lemma pair_add_single (s x : Cfg M) (e : Fin M) :
    ∑ k, s k * (x + Pi.single e 1 : Cfg M) k = ∑ k, s k * x k + s e := by
  simp [mul_add, Finset.sum_add_distrib, Pi.single_apply]

/-! ### Auxiliary lemmas: second quantization -/

/-- `dΓ` is additive. -/
private lemma dGamma_add (A B : Matrix (Fin N) (Fin N) ℝ) :
    dGamma (A + B) = dGamma A + dGamma B := by
  simp only [dGamma, Matrix.add_apply, add_smul, Finset.sum_add_distrib]

/-- `dΓ` is homogeneous. -/
private lemma dGamma_smul (c : ℝ) (A : Matrix (Fin N) (Fin N) ℝ) :
    dGamma (c • A) = c • dGamma A := by
  simp only [dGamma, Matrix.smul_apply, smul_eq_mul, mul_smul, Finset.smul_sum]

/-- Second quantization of a diagonal one-body matrix. -/
private lemma dGamma_diagonal (d : Fin N → ℝ) :
    dGamma (diagonal d) = diagonal (fun n => ∑ x, d x * occ (n x)) := by
  have e : ∀ x, ∑ y, diagonal d x y • (cre x * ann y) = d x • diagonal (fun n => occ (n x)) := by
    intro x
    rw [Finset.sum_eq_single x]
    · rw [diagonal_apply_eq, cre_mul_ann_self]
    · intro y _ hy; rw [diagonal_apply_ne _ (Ne.symm hy), zero_smul]
    · simp
  simp only [dGamma, e]
  ext n' n
  simp only [Matrix.sum_apply, Matrix.smul_apply, diagonal_apply, smul_eq_mul]
  split_ifs <;> simp

/-- The entries of `K_s` as a sum over all edges. -/
private lemma Ks_entry_eq_sum (J : ℝ) (s : Fin M → ℝ) (i j : Fin N) :
    G.Ks J s i j = -J * ∑ e, ((if G.src e = i ∧ G.tgt e = j then s e else 0) +
      (if G.tgt e = i ∧ G.src e = j then s e else 0)) := by
  simp only [Ks, Finset.sum_filter]
  congr 1
  refine Finset.sum_congr rfl (fun e _ => ?_)
  by_cases h1 : G.src e = i ∧ G.tgt e = j
  · have h2 : ¬ (G.tgt e = i ∧ G.src e = j) := by
      rintro ⟨h3, -⟩
      exact G.src_ne_tgt e (h1.1.trans h3.symm)
    rw [ite_eq_left (Or.inl h1), ite_eq_left h1, ite_eq_right h2, add_zero]
  · by_cases h2 : G.tgt e = i ∧ G.src e = j
    · rw [ite_eq_left (Or.inr h2), ite_eq_right h1, ite_eq_left h2, zero_add]
    · rw [ite_eq_right (not_or.2 ⟨h1, h2⟩), ite_eq_right h1, ite_eq_right h2, add_zero]

/-- Second quantization of the signed hopping: `dΓ(K_s) = -J ∑_e s_e (c†_{s(e)} c_{t(e)} + h.c.)`. -/
private lemma dGamma_Ks_eq (J : ℝ) (s : Fin M → ℝ) :
    dGamma (G.Ks J s) = (-J) • ∑ e, s e • G.hopE e := by
  have e1 : ∀ i j, G.Ks J s i j • (cre i * ann j) =
      (-J) • ∑ e, ((if G.src e = i ∧ G.tgt e = j then s e • (cre i * ann j) else 0) +
        (if G.tgt e = i ∧ G.src e = j then s e • (cre i * ann j) else 0)) := by
    intro i j
    rw [Ks_entry_eq_sum, mul_smul, Finset.sum_smul]
    congr 1
    refine Finset.sum_congr rfl (fun e _ => ?_)
    rw [add_smul]
    split_ifs <;> simp
  have e2 : ∀ e, s e • G.hopE e = ∑ i, ∑ j,
      ((if G.src e = i ∧ G.tgt e = j then s e • (cre i * ann j) else 0) +
        (if G.tgt e = i ∧ G.src e = j then s e • (cre i * ann j) else 0)) := by
    intro e
    simp only [Finset.sum_add_distrib, ite_and, Finset.sum_ite_irrel, Finset.sum_const_zero,
      Finset.sum_ite_eq, Finset.mem_univ, ite_true, hopE, smul_add]
  simp only [dGamma, e1, e2]
  simp only [← Finset.smul_sum]
  congr 1
  symm
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [Finset.sum_comm]

/-- `dΓ(h_{q,s})` in terms of the edge hoppings. -/
private lemma dGamma_hqs (J h : ℝ) (q : Cfg N) (s : Fin M → ℝ) :
    dGamma (G.hqs J h (fun i => sgn (q i)) s) =
      (-J) • ∑ e, s e • G.hopE e + (2 * h) • diagonal (fun n => ∑ i, sgn (q i) * occ (n i)) := by
  rw [FinGraph.hqs, dGamma_add, dGamma_smul, dGamma_Ks_eq, dGamma_diagonal]

/-! ### The signed intertwiner `V_s` -/

/-- The (real, rectangular) signed intertwiner `V_s |n, q⟩ = ∑_{x : n + ∂x = q} (-1)^{⟨s,x⟩} |n, x⟩`:
the entry `((n, x), b)` is `(-1)^{⟨s,x⟩}` iff `b = (n, n + ∂x)`. -/
private def VSmat (s : Cfg M) : Matrix (Cfg N × Cfg M) (Cfg N × Cfg N) ℝ :=
  fun a b => if b = (a.1, a.1 + G.bdry a.2) then sgn (∑ e, s e * a.2 e) else 0

/-- Right multiplication by `V_s`. -/
private lemma mul_VSmat_apply (s : Cfg M) (A : Matrix (Cfg N × Cfg M) (Cfg N × Cfg M) ℝ)
    (a : Cfg N × Cfg M) (n0 q : Cfg N) :
    (A * G.VSmat s) a (n0, q) =
      ∑ x1, if q = n0 + G.bdry x1 then A a (n0, x1) * sgn (∑ e, s e * x1 e) else 0 := by
  rw [Matrix.mul_apply, Fintype.sum_prod_type]
  simp only [VSmat, Prod.mk.injEq, mul_ite, mul_zero]
  rw [Finset.sum_eq_single n0]
  · simp
  · intro n1 _ hn1
    apply Finset.sum_eq_zero
    intro x1 _
    rw [ite_eq_right]
    rintro ⟨h1, -⟩
    exact hn1 h1.symm
  · simp

/-- Left multiplication by `V_s`. -/
private lemma VSmat_mul_apply (s : Cfg M) (B : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ)
    (n : Cfg N) (x : Cfg M) (b : Cfg N × Cfg N) :
    (G.VSmat s * B) (n, x) b = sgn (∑ e, s e * x e) * B (n, n + G.bdry x) b := by
  rw [Matrix.mul_apply]
  simp only [VSmat, ite_mul, zero_mul]
  rw [Finset.sum_ite_eq']
  simp

/-- The hopping terms intertwine with `V_s`, picking up the link sign `(-1)^{s_e}`. -/
private lemma hop_VSmat (s : Cfg M) (e : Fin M) :
    (G.hopE e ⊗ₖ G.szE e) * G.VSmat s =
      G.VSmat s * (sgn (s e) • (G.hopE e ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ))) := by
  ext ⟨n, x⟩ ⟨n0, q⟩
  rw [mul_VSmat_apply, VSmat_mul_apply]
  rw [Finset.sum_eq_single (x + Pi.single e 1)]
  · simp only [kronecker_apply, szE_apply, Matrix.smul_apply, Matrix.one_apply, smul_eq_mul,
      ite_true, mul_one]
    by_cases hn : G.hopE e n n0 = 0
    · simp [hn]
    · have hs := G.hopE_support e hn
      subst hs
      rw [pair_add_single, sgn_add, bdry_add]
      have e1 : n0 + (G.bdry x + G.bdry (Pi.single e 1)) =
          n0 + G.bdry (Pi.single e 1) + G.bdry x := by abel
      rw [e1]
      by_cases hq : q = n0 + G.bdry (Pi.single e 1) + G.bdry x
      · rw [ite_eq_left hq, ite_eq_left hq.symm]
        ring
      · rw [ite_eq_right hq, ite_eq_right (Ne.symm hq)]
        ring
  · intro x1 _ hx1
    simp [szE_apply, hx1]
  · simp

/-- The star terms intertwine with `V_s`: on the fibre `n + ∂x = q`,
`∏_{e ∋ i} σˣ_e = sgn(q_i + n_i)`. -/
private lemma star_VSmat (s : Cfg M) (i : Fin N) :
    ((1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ G.star i) * G.VSmat s =
      G.VSmat s * diagonal (fun b : Cfg N × Cfg N => sgn (b.2 i + b.1 i)) := by
  ext ⟨n, x⟩ ⟨n0, q⟩
  rw [mul_VSmat_apply, VSmat_mul_apply, one_kron_star]
  simp only [diagonal_apply, Prod.mk.injEq]
  rw [Finset.sum_eq_single x]
  · by_cases hn : n = n0
    · subst n0
      by_cases hq : q = n + G.bdry x
      · subst q
        simp only [and_self, ite_true, Pi.add_apply]
        rw [add_comm (n i), add_assoc, z2_add_self, add_zero]
        ring
      · simp [hq, Ne.symm hq]
    · simp [hn]
  · intro x1 _ hx1
    simp [Ne.symm hx1]
  · simp

/-- The dual Hamiltonian with link signs `s`, on (occupations × charges). -/
private def dualHamS (J h : ℝ) (s : Cfg M) : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  (-J) • ∑ e, sgn (s e) • (G.hopE e ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) -
    h • ∑ i, diagonal (fun b : Cfg N × Cfg N => sgn (b.2 i + b.1 i))

/-- **Intertwining**: `H V_s = V_s H_dual,s`. -/
private lemma ham_mul_VSmat (J h : ℝ) (s : Cfg M) :
    G.ham J h * G.VSmat s = G.VSmat s * G.dualHamS J h s := by
  simp only [ham, dualHamS, Matrix.sub_mul, Matrix.mul_sub, Matrix.smul_mul, Matrix.mul_smul,
    Matrix.sum_mul, Matrix.mul_sum, hop_VSmat, star_VSmat]

/-- The signed dual Hamiltonian is block diagonal, with blocks `dΓ(h_{q,s}) - h ∑ q_i`. -/
private lemma dualHamS_eq (J h : ℝ) (s : Cfg M) :
    G.dualHamS J h s = blockDiagonal (fun q =>
      dGamma (G.hqs J h (fun i => sgn (q i)) (fun e => sgn (s e))) -
        (h * ∑ i, sgn (q i)) • (1 : Matrix (Cfg N) (Cfg N) ℝ)) := by
  ext ⟨n, q⟩ ⟨n0, q0⟩
  rw [blockDiagonal_apply, dGamma_hqs]
  simp only [dualHamS, Matrix.sub_apply, Matrix.add_apply, Matrix.smul_apply, Matrix.sum_apply,
    kronecker_apply, diagonal_apply, Matrix.one_apply, Prod.mk.injEq, smul_eq_mul]
  by_cases hq : q = q0
  · subst q0
    by_cases hn : n = n0
    · subst n0
      simp only [and_self, ite_true, mul_one]
      have e : ∀ i, sgn (q i + n i) = sgn (q i) - 2 * (sgn (q i) * occ (n i)) := by
        intro i; rw [sgn_add, sgn_eq_one_sub_two_occ (n i)]; ring
      simp only [e, Finset.sum_sub_distrib, ← Finset.mul_sum]
      ring
    · simp [hn]
  · simp [hq]

/-! ### Admissible pairs -/

/-- Indicator of the admissible pairs `(n, q)`, i.e. `q ∈ n + im ∂`. -/
private def Pmat : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℝ :=
  diagonal (fun b => if b.2 ∈ G.sectors b.1 then 1 else 0)

/-- Changing the occupations by a boundary does not change the reachable sectors. -/
private lemma sectors_add_bdry (n : Cfg N) (y : Cfg M) :
    G.sectors (n + G.bdry y) = G.sectors n := by
  ext q
  simp only [sectors, Finset.mem_image, Finset.mem_univ, true_and]
  constructor
  · rintro ⟨x, rfl⟩
    exact ⟨y + x, by rw [bdry_add, add_assoc]⟩
  · rintro ⟨x, rfl⟩
    exact ⟨x - y, by rw [bdry_sub]; abel⟩

/-- The hopping terms preserve the admissible pairs. -/
private lemma hop_comm_Pmat (e : Fin M) :
    (G.hopE e ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) * G.Pmat =
      G.Pmat * (G.hopE e ⊗ₖ (1 : Matrix (Cfg N) (Cfg N) ℝ)) := by
  ext ⟨n, q⟩ ⟨n0, q0⟩
  rw [Pmat, mul_diagonal, diagonal_mul, kronecker_apply, Matrix.one_apply]
  by_cases hn : G.hopE e n n0 = 0
  · simp [hn]
  · have hs := G.hopE_support e hn
    subst hs
    by_cases hq : q = q0
    · subst q0
      simp only [sectors_add_bdry]
      ring
    · simp [hq]

/-- Diagonal matrices commute with `P`. -/
private lemma diagonal_comm_Pmat (f : Cfg N × Cfg N → ℝ) :
    diagonal f * G.Pmat = G.Pmat * diagonal f := by
  rw [Pmat, diagonal_mul_diagonal, diagonal_mul_diagonal]
  congr 1
  funext b
  ring

/-- The signed dual Hamiltonian preserves the admissible pairs. -/
private lemma dualHamS_comm_Pmat (J h : ℝ) (s : Cfg M) :
    G.dualHamS J h s * G.Pmat = G.Pmat * G.dualHamS J h s := by
  simp only [dualHamS, Matrix.sub_mul, Matrix.mul_sub, Matrix.smul_mul, Matrix.mul_smul,
    Matrix.sum_mul, Matrix.mul_sum, hop_comm_Pmat, diagonal_comm_Pmat]

/-! ### Exponentials of intertwined matrices -/

/-- The exponential series of a complex matrix. -/
private lemma hasSum_exp_matrix {ι : Type*} [Fintype ι] [DecidableEq ι] (A : Matrix ι ι ℂ) :
    HasSum (fun k : ℕ => ((k.factorial : ℂ)⁻¹) • A ^ k) (NormedSpace.exp A) := by
  open scoped Matrix.Norms.Operator in exact NormedSpace.exp_series_hasSum_exp' A

/-- If `A W = W B` then `exp A W = W exp B` (for a rectangular `W`). -/
private lemma exp_mul_eq_mul_exp {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ]
    [DecidableEq κ] {A : Matrix ι ι ℂ} {B : Matrix κ κ ℂ} {W : Matrix ι κ ℂ}
    (h : A * W = W * B) : NormedSpace.exp A * W = W * NormedSpace.exp B := by
  have hpow : ∀ k : ℕ, A ^ k * W = W * B ^ k := by
    intro k
    induction k with
    | zero => simp
    | succ k ih =>
      rw [pow_succ, pow_succ, Matrix.mul_assoc, h, ← Matrix.mul_assoc, ih, Matrix.mul_assoc]
  have h1 := (hasSum_exp_matrix A).map (Matrix.addMonoidHomMulRight W)
    (continuous_id.matrix_mul continuous_const)
  have h2 := (hasSum_exp_matrix B).map (Matrix.addMonoidHomMulLeft W)
    (continuous_const.matrix_mul continuous_id)
  have e : (⇑(Matrix.addMonoidHomMulRight W) ∘ fun k : ℕ => ((k.factorial : ℂ)⁻¹) • A ^ k) =
      (⇑(Matrix.addMonoidHomMulLeft W) ∘ fun k : ℕ => ((k.factorial : ℂ)⁻¹) • B ^ k) := by
    funext k
    simp only [Function.comp_apply, Matrix.addMonoidHomMulRight_apply,
      Matrix.addMonoidHomMulLeft_apply, Matrix.smul_mul, Matrix.mul_smul, hpow]
  rw [e] at h1
  exact h1.unique h2

/-- If `A W = W B` (real matrices) then `e^{-itA} W = W e^{-itB}`. -/
private lemma evo_mul_eq_mul_evo {ι κ : Type*} [Fintype ι] [DecidableEq ι] [Fintype κ]
    [DecidableEq κ] {A : Matrix ι ι ℝ} {B : Matrix κ κ ℝ} {W : Matrix ι κ ℝ}
    (h : A * W = W * B) (t : ℝ) :
    evo A t * W.map (fun r => (r : ℂ)) = W.map (fun r => (r : ℂ)) * evo B t := by
  unfold evo
  apply exp_mul_eq_mul_exp
  rw [Matrix.smul_mul, Matrix.mul_smul]
  congr 1
  have := congrArg (fun X => X.map Complex.ofRealHom) h
  simp only [Matrix.map_mul] at this
  exact this

/-! ### States and expectations -/

/-- Applying the complexified `V_s`. -/
private lemma VSc_mulVec_apply (s : Cfg M) (φ : Cfg N × Cfg N → ℂ) (n : Cfg N) (x : Cfg M) :
    ((G.VSmat s).map (fun r => (r : ℂ)) *ᵥ φ) (n, x) =
      (sgn (∑ e, s e * x e) : ℂ) * φ (n, n + G.bdry x) := by
  simp only [mulVec, dotProduct, Matrix.map_apply, VSmat]
  rw [Finset.sum_eq_single (n, n + G.bdry x)]
  · simp
  · intro b _ hb
    simp [hb]
  · simp

/-- Applying the complexified `P`. -/
private lemma Pc_mulVec_apply (φ : Cfg N × Cfg N → ℂ) (b : Cfg N × Cfg N) :
    (G.Pmat.map (fun r => (r : ℂ)) *ᵥ φ) b = if b.2 ∈ G.sectors b.1 then φ b else 0 := by
  rw [Pmat, diagonal_map (by simp), mulVec_diagonal]
  split_ifs <;> simp

/-- Matrix–vector product with a block-diagonal matrix, block by block. -/
private lemma blockDiagonal_mulVec_apply {m o : Type*} [Fintype m] [Fintype o] [DecidableEq o]
    (Mb : o → Matrix m m ℂ) (v : m × o → ℂ) (i : m) (k : o) :
    (blockDiagonal Mb *ᵥ v) (i, k) = (Mb k *ᵥ fun j => v (j, k)) i := by
  simp only [mulVec, dotProduct, blockDiagonal_apply', Fintype.sum_prod_type]
  refine Finset.sum_congr rfl fun j _ => ?_
  simp [ite_mul]

/-- Multiplying the state by `z` multiplies expectations by `|z|²`. -/
private lemma expect_smul {ι : Type*} [Fintype ι] (Mx : Matrix ι ι ℂ) (z : ℂ) (v : ι → ℂ) :
    expect Mx (z • v) = (Star.star z * z) * expect Mx v := by
  unfold expect
  rw [mulVec_smul, star_smul, smul_dotProduct, dotProduct_smul, smul_eq_mul, smul_eq_mul]
  ring

/-- Expectation of a (complexified) real diagonal observable. -/
private lemma expect_cplx_diagonal {ι : Type*} [Fintype ι] [DecidableEq ι] (f : ι → ℝ)
    (v : ι → ℂ) : expect (cplx (diagonal f)) v = ∑ a, Star.star (v a) * ((f a : ℂ) * v a) := by
  unfold expect cplx
  rw [diagonal_map (by simp)]
  simp only [dotProduct, Pi.star_apply, mulVec_diagonal]

/-- `|2^{-M/2}|² = 2^{-M}`. -/
private lemma norm_sq_prefactor :
    Star.star ((((Real.sqrt 2) ^ M : ℝ) : ℂ)⁻¹) * (((Real.sqrt 2) ^ M : ℝ) : ℂ)⁻¹ =
      (1 / 2 ^ M : ℂ) := by
  rw [star_inv₀, Complex.star_def, Complex.conj_ofReal, ← mul_inv, ← Complex.ofReal_mul,
    ← mul_pow, Real.mul_self_sqrt (by norm_num : (0 : ℝ) ≤ 2)]
  push_cast
  ring

/-- A real sign does not change `|z|²`-type expressions: `(±1) z̄ d (±1) z = z̄ d z`. -/
private lemma star_sgn_mul (a : ZMod 2) (z w : ℂ) :
    Star.star ((sgn a : ℂ) * z) * (w * ((sgn a : ℂ) * z)) = Star.star z * (w * z) := by
  have hs : Star.star (sgn a : ℂ) = (sgn a : ℂ) := Complex.conj_ofReal _
  rw [star_mul', hs]
  have h1 : (sgn a : ℂ) * (sgn a : ℂ) = 1 := by
    rw [← Complex.ofReal_mul, sgn_mul_self, Complex.ofReal_one]
  linear_combination (Star.star z * (w * z)) * h1

/-- The number `|ker ∂|` of edge configurations with zero boundary. -/
private def kerCard : ℕ := (univ.filter (fun x : Cfg M => G.bdry x = 0)).card

/-- All nonempty fibres of `x ↦ n + ∂x` have `|ker ∂|` elements. -/
private lemma fibre_card {n q : Cfg N} (hq : q ∈ G.sectors n) :
    (univ.filter (fun x : Cfg M => n + G.bdry x = q)).card = G.kerCard := by
  obtain ⟨x0, -, rfl⟩ := Finset.mem_image.1 hq
  have h1 : univ.filter (fun x : Cfg M => n + G.bdry x = n + G.bdry x0) =
      univ.filter (fun x => G.bdryHom x = G.bdryHom x0) := by
    apply Finset.filter_congr
    intro x _
    simp [bdryHom]
  have h2 : univ.filter (fun x : Cfg M => G.bdry x = 0) =
      univ.filter (fun x => G.bdryHom x = G.bdryHom 0) := by
    apply Finset.filter_congr
    intro x _
    simp [bdryHom, bdry_zero]
  rw [kerCard, h1, h2]
  exact AddMonoidHom.card_fiber_eq_of_mem_range _ ⟨x0, rfl⟩ ⟨0, rfl⟩

/-- `|n + im ∂| · |ker ∂| = 2^M`: the fibres partition the edge configurations. -/
private lemma card_sectors_mul_kerCard (n : Cfg N) :
    (G.sectors n).card * G.kerCard = 2 ^ M := by
  have h := Finset.card_eq_sum_card_image (fun x : Cfg M => n + G.bdry x) univ
  rw [Finset.card_univ, Fintype.card_fun, ZMod.card, Fintype.card_fin,
    show univ.image (fun x : Cfg M => n + G.bdry x) = G.sectors n from rfl] at h
  rw [h, Finset.sum_congr rfl (fun q hq => G.fibre_card hq), Finset.sum_const, smul_eq_mul]

/-- Expectation of `D ⊗ 1` in the state `V_s φ`, for `φ` supported on the admissible pairs. -/
private lemma expect_VSc (s : Cfg M) (d : Cfg N → ℝ) (φ : Cfg N × Cfg N → ℂ)
    (hsupp : ∀ n' q, q ∉ G.sectors n' → φ (n', q) = 0) :
    expect (cplx (diagonal d ⊗ₖ (1 : Matrix (Cfg M) (Cfg M) ℝ)))
        ((G.VSmat s).map (fun r => (r : ℂ)) *ᵥ φ) =
      (G.kerCard : ℂ) * ∑ q, ∑ n', Star.star (φ (n', q)) * ((d n' : ℂ) * φ (n', q)) := by
  have hD : diagonal d ⊗ₖ (1 : Matrix (Cfg M) (Cfg M) ℝ) =
      diagonal (fun p : Cfg N × Cfg M => d p.1) := by
    rw [← diagonal_one, diagonal_kronecker_diagonal]
    simp
  rw [hD, expect_cplx_diagonal, Fintype.sum_prod_type]
  conv_rhs => rw [Finset.sum_comm, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun n' _ => ?_)
  simp only [VSc_mulVec_apply, star_sgn_mul]
  have hcomp : ∑ x : Cfg M,
        Star.star (φ (n', n' + G.bdry x)) * ((d n' : ℂ) * φ (n', n' + G.bdry x)) =
      ∑ q ∈ G.sectors n', (univ.filter (fun x : Cfg M => n' + G.bdry x = q)).card •
        (Star.star (φ (n', q)) * ((d n' : ℂ) * φ (n', q))) :=
    Finset.sum_comp (fun q => Star.star (φ (n', q)) * ((d n' : ℂ) * φ (n', q)))
      (fun x => n' + G.bdry x)
  rw [hcomp, Finset.mul_sum, ← Finset.sum_subset (Finset.subset_univ (G.sectors n'))]
  · refine Finset.sum_congr rfl (fun q hq => ?_)
    rw [G.fibre_card hq, nsmul_eq_mul]
  · intro q _ hq
    rw [hsupp n' q hq]
    simp

/-- **`σᶻ`-product links: Anderson–Bernoulli models in a fixed flux sector** (paper: Proposition
`prop:2dlaws`(a)). -/
theorem sector_average_flux (J h : ℝ) (s : Cfg M) (n : Cfg N) (d : Cfg N → ℝ) (t : ℝ) :
    expect (cplx (diagonal d ⊗ₖ (1 : Matrix (Cfg M) (Cfg M) ℝ)))
        (evo (G.ham J h) t *ᵥ G.initS s (Pi.single n 1)) =
      (1 / ((G.sectors n).card : ℂ)) * ∑ q ∈ G.sectors n,
        expect (cplx (diagonal d))
          (evo (dGamma (G.hqs J h (fun i => sgn (q i)) (fun e => sgn (s e)))) t *ᵥ
            Pi.single n 1) := by
  -- the dual initial state `w = ∑_{q ∈ n + im ∂} |n, q⟩` and its evolution `φ`
  set w : Cfg N × Cfg N → ℂ :=
    fun b => if b.2 ∈ G.sectors n then (Pi.single n 1 : Cfg N → ℂ) b.1 else 0 with hw
  set φ : Cfg N × Cfg N → ℂ := evo (G.dualHamS J h s) t *ᵥ w with hφ
  -- the LGT initial state is `2^{-M/2} V_s w`
  have hinit : G.initS s (Pi.single n 1) =
      (((Real.sqrt 2) ^ M : ℝ) : ℂ)⁻¹ • ((G.VSmat s).map (fun r => (r : ℂ)) *ᵥ w) := by
    funext ⟨n', x⟩
    rw [Pi.smul_apply, VSc_mulVec_apply, smul_eq_mul]
    simp only [initS, hw]
    by_cases hn : n' = n
    · subst n'
      have hmem : n + G.bdry x ∈ G.sectors n := Finset.mem_image_of_mem _ (Finset.mem_univ x)
      rw [ite_eq_left hmem]
      ring
    · simp [hn]
  -- the LGT state at time `t` is `2^{-M/2} V_s φ`
  have hstate : evo (G.ham J h) t *ᵥ G.initS s (Pi.single n 1) =
      (((Real.sqrt 2) ^ M : ℝ) : ℂ)⁻¹ • ((G.VSmat s).map (fun r => (r : ℂ)) *ᵥ φ) := by
    rw [hinit, mulVec_smul, mulVec_mulVec, evo_mul_eq_mul_evo (G.ham_mul_VSmat J h s) t,
      ← mulVec_mulVec]
  -- `φ` is supported on the admissible pairs
  have hsupp : ∀ n' q, q ∉ G.sectors n' → φ (n', q) = 0 := by
    have hPw : G.Pmat.map (fun r => (r : ℂ)) *ᵥ w = w := by
      funext ⟨n', q⟩
      rw [Pc_mulVec_apply]
      simp only [hw]
      by_cases hn : n' = n
      · subst n'
        split_ifs <;> rfl
      · simp [hn]
    have hPφ : G.Pmat.map (fun r => (r : ℂ)) *ᵥ φ = φ := by
      rw [hφ, mulVec_mulVec, ← evo_mul_eq_mul_evo (G.dualHamS_comm_Pmat J h s) t,
        ← mulVec_mulVec, hPw]
    intro n' q hq
    rw [← hPφ, Pc_mulVec_apply, ite_eq_right hq]
  -- `φ` in terms of the sector evolutions
  have hφq : ∀ n' q, φ (n', q) = if q ∈ G.sectors n then
      (evo (dGamma (G.hqs J h (fun i => sgn (q i)) (fun e => sgn (s e))) -
        (h * ∑ i, sgn (q i)) • (1 : Matrix (Cfg N) (Cfg N) ℝ)) t *ᵥ Pi.single n 1) n'
      else 0 := by
    intro n' q
    rw [hφ, dualHamS_eq, evo_blockDiagonal]
    simp only [blockDiagonal_mulVec_apply]
    split_ifs with hq
    · have h1 : (fun j => w (j, q)) = Pi.single n 1 := by
        funext j
        simp [hw, hq]
      rw [h1]
    · have h0 : (fun j => w (j, q)) = 0 := by
        funext j
        simp [hw, hq]
      rw [h0, mulVec_zero]
      rfl
  have hsum : ∑ q, ∑ n', Star.star (φ (n', q)) * ((d n' : ℂ) * φ (n', q)) =
      ∑ q ∈ G.sectors n, expect (cplx (diagonal d))
        (evo (dGamma (G.hqs J h (fun i => sgn (q i)) (fun e => sgn (s e)))) t *ᵥ
          Pi.single n 1) := by
    rw [← Finset.sum_subset (Finset.subset_univ (G.sectors n))]
    · refine Finset.sum_congr rfl (fun q hq => ?_)
      rw [← expect_evo_sub_smul_one
        (dGamma (G.hqs J h (fun i => sgn (q i)) (fun e => sgn (s e)))) (h * ∑ i, sgn (q i)),
        expect_cplx_diagonal]
      refine Finset.sum_congr rfl (fun n' _ => ?_)
      rw [hφq, ite_eq_left hq]
    · intro q _ hq
      apply Finset.sum_eq_zero
      intro n' _
      rw [hφq, ite_eq_right hq]
      simp
  rw [hstate, expect_smul, norm_sq_prefactor, G.expect_VSc s d φ hsupp, hsum, ← mul_assoc]
  congr 1
  -- the prefactor: `2^{-M} |ker ∂| = 1 / |n + im ∂|`
  have hcard : ((G.sectors n).card : ℂ) * (G.kerCard : ℂ) = 2 ^ M := by
    exact_mod_cast G.card_sectors_mul_kerCard n
  have h2 : ((G.sectors n).card : ℂ) * (G.kerCard : ℂ) ≠ 0 := by
    rw [hcard]
    exact pow_ne_zero _ two_ne_zero
  have hK : (G.kerCard : ℂ) ≠ 0 := right_ne_zero_of_mul h2
  rw [← hcard, one_div, mul_inv, mul_assoc, inv_mul_cancel₀ hK, mul_one, one_div]

/-- **`σˣ`-polarized links: random `π`-flux** (paper: Proposition `prop:2dlaws`(b)). The state
`|n⟩ ⊗ |+⋯+⟩` is the basis vector `(n, 0)` of the `σˣ` basis. -/
theorem sector_average_xpol (J h : ℝ) (n : Cfg N) (d : Cfg N → ℝ) (t : ℝ) :
    expect (cplx (diagonal d ⊗ₖ (1 : Matrix (Cfg M) (Cfg M) ℝ)))
        (evo (G.ham J h) t *ᵥ Pi.single (n, 0) 1) =
      (1 / 2 ^ M : ℂ) * ∑ s : Cfg M,
        expect (cplx (diagonal d))
          (evo (dGamma (G.hqs J h (fun i => sgn (n i)) (fun e => sgn (s e)))) t *ᵥ
            Pi.single n 1) := by
  sorry

/-! ### Gauge covariance of the free-fermion dynamics -/

/-- The second quantization `Γ(D_θ) = diag((-1)^{∑_i θ_i m_i})` of `D_θ = diag((-1)^{θ_i})`. -/
private def gaugeSign (θ : Cfg N) : Matrix (Cfg N) (Cfg N) ℝ :=
  diagonal (fun m => sgn (∑ i, θ i * m i))

/-- `Γ(D_θ)² = 1`. -/
private lemma gaugeSign_mul_self (θ : Cfg N) : gaugeSign θ * gaugeSign θ = 1 := by
  rw [gaugeSign, diagonal_mul_diagonal, ← diagonal_one]
  congr 1
  funext m
  exact sgn_mul_self _

/-- Removing the particle at `y` changes the gauge sign by `(-1)^{θ_y}`. -/
private lemma sgn_pair_update (θ n : Cfg N) (y : Fin N) (hy : n y = 1) :
    sgn (∑ i, θ i * Function.update n y 0 i) * sgn (∑ i, θ i * n i) = sgn (θ y) := by
  have hn : Function.update n y 0 + Pi.single y 1 = n := by
    funext k
    by_cases hk : k = y
    · subst hk; simp [hy]
    · simp [hk]
  have h2 : ∑ i, θ i * n i = ∑ i, θ i * Function.update n y 0 i + θ y := by
    have := pair_add_single θ (Function.update n y 0) y
    rw [hn] at this
    exact this
  rw [h2, sgn_add, ← mul_assoc, sgn_mul_self, one_mul]

/-- `Γ(D_θ) c_y Γ(D_θ) = (-1)^{θ_y} c_y`. -/
private lemma gaugeSign_ann (θ : Cfg N) (y : Fin N) :
    gaugeSign θ * ann y * gaugeSign θ = sgn (θ y) • ann y := by
  ext n' n
  rw [gaugeSign, mul_diagonal, diagonal_mul, Matrix.smul_apply, smul_eq_mul]
  simp only [ann]
  split_ifs with hc
  · obtain ⟨hy, rfl⟩ := hc
    have := sgn_pair_update θ n y hy
    linear_combination (sgn (jwPar n y)) * this
  · simp

/-- `Γ(D_θ) c_y† Γ(D_θ) = (-1)^{θ_y} c_y†`. -/
private lemma gaugeSign_cre (θ : Cfg N) (y : Fin N) :
    gaugeSign θ * cre y * gaugeSign θ = sgn (θ y) • cre y := by
  have h := congrArg Matrix.transpose (gaugeSign_ann θ y)
  rw [Matrix.transpose_mul, Matrix.transpose_mul, Matrix.transpose_smul] at h
  have hT : (gaugeSign θ)ᵀ = gaugeSign θ := diagonal_transpose _
  rw [hT, ← Matrix.mul_assoc] at h
  exact h

/-- `Γ(D_θ) dΓ(A) Γ(D_θ) = dΓ(D_θ A D_θ)`. -/
private lemma gaugeSign_dGamma (θ : Cfg N) (A : Matrix (Fin N) (Fin N) ℝ) :
    gaugeSign θ * dGamma A * gaugeSign θ =
      dGamma (diagonal (fun i => sgn (θ i)) * A * diagonal (fun i => sgn (θ i))) := by
  have key : ∀ x y, gaugeSign θ * (cre x * ann y) * gaugeSign θ =
      (sgn (θ x) * sgn (θ y)) • (cre x * ann y) := by
    intro x y
    have e1 : gaugeSign θ * (cre x * ann y) * gaugeSign θ =
        (gaugeSign θ * cre x * gaugeSign θ) * (gaugeSign θ * ann y * gaugeSign θ) := by
      have e2 : (gaugeSign θ * cre x * gaugeSign θ) * (gaugeSign θ * ann y * gaugeSign θ) =
          gaugeSign θ * cre x * (gaugeSign θ * gaugeSign θ) * ann y * gaugeSign θ := by
        simp only [Matrix.mul_assoc]
      rw [e2, gaugeSign_mul_self, Matrix.mul_one]
      simp only [Matrix.mul_assoc]
    rw [e1, gaugeSign_cre, gaugeSign_ann, Matrix.smul_mul, Matrix.mul_smul, smul_smul]
  simp only [dGamma, Matrix.mul_sum, Matrix.sum_mul, Matrix.mul_smul, Matrix.smul_mul, key,
    smul_smul, Matrix.mul_diagonal, Matrix.diagonal_mul]
  refine Finset.sum_congr rfl (fun x _ => Finset.sum_congr rfl (fun y _ => ?_))
  congr 1
  ring

/-- **Gauge covariance** (as in `gauge_hqs`): `h_{q,θ·s} = D_θ h_{q,s} D_θ`. -/
private lemma hqs_gaugeAct (J h : ℝ) (q : Fin N → ℝ) (s : Cfg M) (θ : Cfg N) :
    G.hqs J h q (fun e => sgn (G.gaugeAct θ s e)) =
      diagonal (fun i => sgn (θ i)) * G.hqs J h q (fun e => sgn (s e)) *
        diagonal (fun i => sgn (θ i)) := by
  ext i j
  simp only [hqs, Matrix.mul_diagonal, Matrix.diagonal_mul, Matrix.add_apply, Matrix.smul_apply,
    smul_eq_mul, Ks]
  have hK : ∑ e ∈ univ.filter (fun e => (G.src e = i ∧ G.tgt e = j) ∨
      (G.tgt e = i ∧ G.src e = j)), sgn (G.gaugeAct θ s e) =
      sgn (θ i) * sgn (θ j) * ∑ e ∈ univ.filter (fun e => (G.src e = i ∧ G.tgt e = j) ∨
      (G.tgt e = i ∧ G.src e = j)), sgn (s e) := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl (fun e he => ?_)
    rcases (Finset.mem_filter.mp he).2 with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · simp only [gaugeAct, sgn_add, h1, h2]; ring
    · simp only [gaugeAct, sgn_add, h1, h2]; ring
  have hD : sgn (θ i) * (2 * h * diagonal q i j) * sgn (θ j) = 2 * h * diagonal q i j := by
    by_cases hij : i = j
    · subst hij
      have := sgn_mul_self (θ i)
      linear_combination (2 * h * diagonal q i i) * this
    · rw [diagonal_apply_ne _ hij]; ring
  rw [hK]
  linear_combination (-1 : ℝ) * hD

/-- Expectations of diagonal observables are invariant under the diagonal sign `Γ(D_θ)`. -/
private lemma expect_diag_gaugeSign (θ : Cfg N) (d : Cfg N → ℝ) (u : Cfg N → ℂ) :
    expect (cplx (diagonal d)) ((gaugeSign θ).map (fun r => (r : ℂ)) *ᵥ u) =
      expect (cplx (diagonal d)) u := by
  rw [expect_cplx_diagonal, expect_cplx_diagonal, gaugeSign, diagonal_map (by simp)]
  refine Finset.sum_congr rfl (fun m _ => ?_)
  rw [mulVec_diagonal]
  exact star_sgn_mul _ _ _

/-- The free-fermion expectations depend on the signs only through their gauge class (hence only
through the fluxes, by `flux_eq_iff_gauge`). -/
theorem xpol_summand_gauge (J h : ℝ) (q n : Cfg N) (d : Cfg N → ℝ) (t : ℝ) (s : Cfg M)
    (θ : Cfg N) :
    expect (cplx (diagonal d))
        (evo (dGamma (G.hqs J h (fun i => sgn (q i)) (fun e => sgn (G.gaugeAct θ s e)))) t *ᵥ
          Pi.single n 1) =
      expect (cplx (diagonal d))
        (evo (dGamma (G.hqs J h (fun i => sgn (q i)) (fun e => sgn (s e)))) t *ᵥ
          Pi.single n 1) := by
  have hA : dGamma (G.hqs J h (fun i => sgn (q i)) (fun e => sgn (G.gaugeAct θ s e))) =
      gaugeSign θ * dGamma (G.hqs J h (fun i => sgn (q i)) (fun e => sgn (s e))) *
        gaugeSign θ := by
    rw [hqs_gaugeAct, gaugeSign_dGamma]
  have hAW : gaugeSign θ * dGamma (G.hqs J h (fun i => sgn (q i)) (fun e => sgn (s e))) *
      gaugeSign θ * gaugeSign θ =
      gaugeSign θ * dGamma (G.hqs J h (fun i => sgn (q i)) (fun e => sgn (s e))) := by
    rw [Matrix.mul_assoc _ (gaugeSign θ) (gaugeSign θ), gaugeSign_mul_self, Matrix.mul_one]
  have hevo := evo_mul_eq_mul_evo hAW t
  -- `|n⟩ = ±Γ(D_θ)|n⟩`
  have hsingle : (Pi.single n 1 : Cfg N → ℂ) =
      (sgn (∑ i, θ i * n i) : ℂ) • ((gaugeSign θ).map (fun r => (r : ℂ)) *ᵥ Pi.single n 1) := by
    funext m
    rw [Pi.smul_apply, gaugeSign, diagonal_map (by simp), mulVec_diagonal, smul_eq_mul]
    by_cases hm : m = n
    · subst hm
      simp only [Pi.single_eq_same, mul_one]
      rw [← Complex.ofReal_mul, sgn_mul_self, Complex.ofReal_one]
    · simp [hm]
  rw [hA]
  conv_lhs => rw [hsingle]
  rw [mulVec_smul, mulVec_mulVec, hevo, ← mulVec_mulVec, expect_smul, expect_diag_gaugeSign]
  have h1 : Star.star (sgn (∑ i, θ i * n i) : ℂ) * (sgn (∑ i, θ i * n i) : ℂ) = 1 := by
    have hs : Star.star (sgn (∑ i, θ i * n i) : ℂ) = (sgn (∑ i, θ i * n i) : ℂ) :=
      Complex.conj_ofReal _
    rw [hs, ← Complex.ofReal_mul, sgn_mul_self, Complex.ofReal_one]
  rw [h1, one_mul]

end FinGraph

end DFL
