import DFL.GraphWilson

/-!
# The exact mapping on graphs: joint eigenvectors of the Gauss operators and Wilson loops
(paper: Theorem `thm:2dmap`(b))

For a connected finite graph, charges `q ∈ Cfg N`, link signs `s ∈ Cfg M` and occupations `n` with
`∑ n = ∑ q`, the vector
`Ψ_{n,q,s} = |ker ∂|^{-1/2} ∑_{x : n + ∂x = q} (-1)^{⟨s,x⟩} |n, x⟩`   (`fluxVec q s n`)
is a normalized joint eigenvector of the occupations, of every Gauss operator (eigenvalue `q_i`) and
of every Wilson loop `W_z`, `z ∈ ker ∂` (eigenvalue `(-1)^{⟨s,z⟩}`, i.e. the flux of `s`). Vectors with
different labels `(q, n)` or different fluxes are orthogonal, gauge-equivalent signs give the same
vector up to the sign `q(θ)(-1)^{⟨θ,n⟩}` (paper: `V_{q,θ·s} = q(θ) V_{q,s} Γ(D_θ)`), and there are
`2^N · 2^β · 2^{N-1} = dim` labels, so one representative `s` per flux class gives an orthonormal
basis. In this basis `H` acts as `dΓ(h_{q,s}) - h ∑ q_i` (paper: Eqs. `eq:2ddict`, `eq:2dsector`).
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

namespace FinGraph

variable {N M : ℕ} (G : FinGraph N M)

/-- The number of cycles, `|ker ∂|`. -/
def cycleCard : ℕ := (univ.filter (fun z : Cfg M => G.bdry z = 0)).card

/-- The flux classes (gauge orbits of sign configurations). -/
def fluxClasses : Finset (Finset (Cfg M)) :=
  univ.image (fun s : Cfg M => univ.image (fun θ : Cfg N => G.gaugeAct θ s))

/-- The joint eigenvector `Ψ_{n,q,s}`. -/
def fluxVec (q : Cfg N) (s : Cfg M) (n : Cfg N) : Cfg N × Cfg M → ℂ :=
  fun p => if p.1 = n ∧ n + G.bdry p.2 = q then
    ((Real.sqrt G.cycleCard : ℝ) : ℂ)⁻¹ * ((sgn (∑ e, s e * p.2 e) : ℝ) : ℂ) else 0

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

/-- `(-1)^a = 1 - 2a` for `a ∈ {0, 1}`. -/
private lemma sgn_eq_one_sub_two_occ (a : ZMod 2) : sgn a = 1 - 2 * occ a := by
  rcases z2_cases a with rfl | rfl <;> norm_num [sgn, occ]

/-- Different elements of `ZMod 2` have opposite signs. -/
private lemma sgn_mul_sgn_of_ne {a b : ZMod 2} (h : a ≠ b) : sgn a * sgn b = -1 := by
  rcases z2_cases a with rfl | rfl <;> rcases z2_cases b with rfl | rfl
  · exact absurd rfl h
  · simp [sgn]
  · simp [sgn]
  · exact absurd rfl h

/-- `z + z = 0` for configurations. -/
private lemma cfg_add_self {K : ℕ} (z : Cfg K) : z + z = 0 := by
  funext e
  exact z2_add_self (z e)

/-- `x + c + c = x` for configurations. -/
private lemma cfg_add_add_self {K : ℕ} (x c : Cfg K) : x + c + c = x := by
  rw [add_assoc, cfg_add_self, add_zero]

/-- `a + b + (b + c) = a + c` for configurations. -/
private lemma cfg_add_add_cancel {K : ℕ} (a b c : Cfg K) : a + b + (b + c) = a + c := by
  rw [add_assoc, ← add_assoc b b c, cfg_add_self, zero_add]

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

/-- The boundary map as a `ZMod 2`-linear map. -/
private def bdryLin : Cfg M →ₗ[ZMod 2] Cfg N where
  toFun := G.bdry
  map_add' := G.bdry_add
  map_smul' := by
    intro c x
    funext i
    simp [bdry, Finset.mul_sum]

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
    · rw [ite_eq_right (not_or.2 ⟨h1, h2⟩), ite_eq_right (Ne.symm h1),
        ite_eq_right (Ne.symm h2), add_zero]

/-- Duality of `δ` and `∂`: `⟨δθ, z⟩ = ⟨θ, ∂z⟩`. -/
private lemma sum_coboundary_mul (θ : Cfg N) (z : Cfg M) :
    ∑ e, (θ (G.src e) + θ (G.tgt e)) * z e = ∑ i, θ i * G.bdry z i := by
  simp only [bdry, Finset.sum_filter, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun e _ => ?_)
  have h := G.sum_ite_endpoint e (fun i => θ i * z e)
  simp only [mul_ite, mul_zero]
  rw [h, add_mul]

/-- `⟨s, x + y⟩ = ⟨s, x⟩ + ⟨s, y⟩`. -/
private lemma pair_add (s x y : Cfg M) :
    ∑ e, s e * (x + y) e = ∑ e, s e * x e + ∑ e, s e * y e := by
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl (fun e _ => ?_)
  rw [Pi.add_apply, mul_add]

/-! ### Auxiliary lemmas: the image of `∂` on a connected graph, and the fibres -/

/-- Along a path from `i` to `j` there is an edge set with boundary `δ_i + δ_j`. -/
private lemma exists_bdry_pair {i j : Fin N} (h : Relation.ReflTransGen G.Adj i j) :
    ∃ x : Cfg M, G.bdry x = Pi.single i 1 + Pi.single j 1 := by
  induction h with
  | refl => exact ⟨0, by rw [bdry_zero, cfg_add_self]⟩
  | tail _ hbc ih =>
    obtain ⟨x, hx⟩ := ih
    obtain ⟨e, he⟩ := hbc
    refine ⟨x + Pi.single e 1, ?_⟩
    rw [bdry_add, hx, bdry_single_eq]
    rcases he with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · rw [h1, h2]
      exact cfg_add_add_cancel _ _ _
    · rw [add_comm (Pi.single (G.src e) 1), h1, h2]
      exact cfg_add_add_cancel _ _ _

/-- On a connected graph every `y` with `∑ y = 0` is a boundary. -/
private lemma exists_bdry_eq (hG : G.Connected) (y : Cfg N) (hy : ∑ i, y i = 0) :
    ∃ x : Cfg M, G.bdry x = y := by
  rcases Nat.eq_zero_or_pos N with hN | hN
  · subst hN
    exact ⟨0, Subsingleton.elim _ _⟩
  · choose xp hxp using fun i => G.exists_bdry_pair (hG i ⟨0, hN⟩)
    refine ⟨∑ i, y i • xp i, ?_⟩
    change G.bdryLin (∑ i, y i • xp i) = y
    rw [map_sum]
    simp only [map_smul]
    have h1 : ∀ i, G.bdryLin (xp i) = Pi.single i 1 + Pi.single ⟨0, hN⟩ 1 := hxp
    simp only [h1, smul_add, Finset.sum_add_distrib, ← Finset.sum_smul, hy, zero_smul, add_zero]
    funext j
    simp [Finset.sum_apply, Pi.single_apply]

/-- All nonempty fibres of `x ↦ n + ∂x` have `|ker ∂|` elements. -/
private lemma fibre_card {n q : Cfg N} {x0 : Cfg M} (hx0 : n + G.bdry x0 = q) :
    (univ.filter (fun x : Cfg M => n + G.bdry x = q)).card = G.cycleCard := by
  subst hx0
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
  rw [cycleCard, h1, h2]
  exact AddMonoidHom.card_fiber_eq_of_mem_range _ ⟨x0, rfl⟩ ⟨0, rfl⟩

/-- `|ker ∂| > 0`. -/
private lemma cycleCard_pos : 0 < G.cycleCard :=
  Finset.card_pos.2 ⟨0, Finset.mem_filter.2 ⟨Finset.mem_univ _, G.bdry_zero⟩⟩

/-! ### Auxiliary lemmas: translations by cycles -/

/-- Translating the edge labels by a cycle `z` multiplies `Ψ_{n,q,s}` by `(-1)^{⟨s,z⟩}`. -/
private lemma fluxVec_add_cycle (q : Cfg N) (s : Cfg M) (n : Cfg N) {z : Cfg M}
    (hz : G.bdry z = 0) (m : Cfg N) (x : Cfg M) :
    G.fluxVec q s n (m, x + z) = ((sgn (∑ e, s e * z e) : ℝ) : ℂ) * G.fluxVec q s n (m, x) := by
  simp only [fluxVec, G.bdry_add, hz, add_zero]
  split_ifs
  · rw [pair_add, sgn_add]
    push_cast
    ring
  · simp

/-- The translation `x ↦ x + z` of edge labels, as a matrix. -/
private def translM (z : Cfg M) : Matrix (Cfg M) (Cfg M) ℝ :=
  fun x' x => if x' = x + z then 1 else 0

/-- The Wilson loop is `1 ⊗ (translation by z)`. -/
private lemma wilsonLoop_eq_transl (z : Cfg M) :
    G.wilsonLoop z = (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ translM z := rfl

/-- `(W_z v)(m, x) = v(m, x + z)`. -/
private lemma wilsonLoop_mulVec_apply (z : Cfg M) (v : Cfg N × Cfg M → ℂ) (m : Cfg N)
    (x : Cfg M) : (cplx (G.wilsonLoop z) *ᵥ v) (m, x) = v (m, x + z) := by
  simp only [mulVec, dotProduct]
  rw [Fintype.sum_prod_type]
  simp only [cplx, Matrix.map_apply, wilsonLoop_eq_transl, kronecker_apply, translM,
    Matrix.one_apply]
  rw [Finset.sum_eq_single m]
  · rw [Finset.sum_eq_single (x + z)]
    · simp [cfg_add_add_self]
    · intro x' _ hx'
      have hne : ¬ x = x' + z := by
        rintro rfl
        exact hx' (cfg_add_add_self x' z).symm
      simp [hne]
    · simp
  · intro m' _ hm'
    simp [Ne.symm hm']
  · simp

/-- `Ψ_{n,q,s}` is an eigenvector of the occupations with eigenvalues `n`. -/
theorem fluxVec_occ (q : Cfg N) (s : Cfg M) (n : Cfg N) (i : Fin N) :
    cplx (diagonal (fun p : Cfg N × Cfg M => occ (p.1 i))) *ᵥ G.fluxVec q s n =
      ((occ (n i) : ℝ) : ℂ) • G.fluxVec q s n := by
  funext ⟨m, x⟩
  rw [cplx, diagonal_map (by simp), mulVec_diagonal, Pi.smul_apply, smul_eq_mul]
  simp only [fluxVec]
  split_ifs with hc
  · rw [hc.1]
  · simp

/-- `G_i Ψ_{n,q,s} = q_i Ψ_{n,q,s}`. -/
theorem fluxVec_gauss (q : Cfg N) (s : Cfg M) (n : Cfg N) (i : Fin N) :
    cplx (G.gaussOp i) *ᵥ G.fluxVec q s n = ((sgn (q i) : ℝ) : ℂ) • G.fluxVec q s n := by
  funext ⟨m, x⟩
  rw [gaussOp, cplx, diagonal_map (by simp), mulVec_diagonal, Pi.smul_apply, smul_eq_mul]
  simp only [fluxVec]
  split_ifs with hc
  · rw [← sgn_add, ← Pi.add_apply, hc.1, hc.2]
  · simp

/-- `W_z Ψ_{n,q,s} = (-1)^{⟨s,z⟩} Ψ_{n,q,s}` for every cycle `z`. -/
theorem fluxVec_wilson (q : Cfg N) (s : Cfg M) (n : Cfg N) (z : Cfg M) (hz : G.bdry z = 0) :
    cplx (G.wilsonLoop z) *ᵥ G.fluxVec q s n =
      ((sgn (∑ e, s e * z e) : ℝ) : ℂ) • G.fluxVec q s n := by
  funext ⟨m, x⟩
  rw [wilsonLoop_mulVec_apply, G.fluxVec_add_cycle q s n hz, Pi.smul_apply, smul_eq_mul]

/-- **Gauge covariance** (paper: `V_{q,θ·s} = q(θ) V_{q,s} Γ(D_θ)`). -/
theorem fluxVec_gaugeAct (q : Cfg N) (s : Cfg M) (n : Cfg N) (θ : Cfg N) :
    G.fluxVec q (G.gaugeAct θ s) n = ((sgn (∑ i, θ i * (q i + n i)) : ℝ) : ℂ) • G.fluxVec q s n := by
  funext ⟨m, x⟩
  simp only [fluxVec, Pi.smul_apply, smul_eq_mul]
  split_ifs with hc
  · obtain ⟨-, h2⟩ := hc
    have hb : G.bdry x = q + n := by
      rw [← h2, add_comm n, add_assoc, cfg_add_self, add_zero]
    have h3 : ∑ e, G.gaugeAct θ s e * x e = ∑ e, s e * x e + ∑ i, θ i * (q i + n i) := by
      have h4 : ∑ i, θ i * (q i + n i) = ∑ i, θ i * G.bdry x i := by
        rw [hb]
        rfl
      rw [h4, ← G.sum_coboundary_mul, ← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl (fun e _ => ?_)
      simp only [gaugeAct]
      ring
    rw [h3, sgn_add]
    push_cast
    ring
  · simp

/-- `Ψ_{n,q,s}` is normalized when `∑ n = ∑ q` (connected graph). -/
theorem fluxVec_norm (hG : G.Connected) (q : Cfg N) (s : Cfg M) (n : Cfg N)
    (hn : ∑ i, n i = ∑ i, q i) :
    Star.star (G.fluxVec q s n) ⬝ᵥ G.fluxVec q s n = 1 := by
  -- the fibre `{x | n + ∂x = q}` is nonempty
  obtain ⟨x0, hx0⟩ := G.exists_bdry_eq hG (q - n) (by
    simp only [Pi.sub_apply, Finset.sum_sub_distrib, hn, sub_self])
  have hx0' : n + G.bdry x0 = q := by
    rw [hx0]
    abel
  have hfib := G.fibre_card hx0'
  have hK : (0 : ℝ) < G.cycleCard := by exact_mod_cast G.cycleCard_pos
  -- each term of the sum
  have hterm : ∀ p : Cfg N × Cfg M, Star.star (G.fluxVec q s n p) * G.fluxVec q s n p =
      if p.1 = n ∧ n + G.bdry p.2 = q then ((G.cycleCard : ℂ))⁻¹ else 0 := by
    intro p
    simp only [fluxVec]
    split_ifs
    · simp only [star_mul', star_inv₀, Complex.star_def, Complex.conj_ofReal]
      have hσ : ((sgn (∑ e, s e * p.2 e) : ℝ) : ℂ) * ((sgn (∑ e, s e * p.2 e) : ℝ) : ℂ) = 1 := by
        rw [← Complex.ofReal_mul, sgn_mul_self, Complex.ofReal_one]
      have hr : ((Real.sqrt G.cycleCard : ℝ) : ℂ) * ((Real.sqrt G.cycleCard : ℝ) : ℂ) =
          (G.cycleCard : ℂ) := by
        rw [← Complex.ofReal_mul, Real.mul_self_sqrt hK.le, Complex.ofReal_natCast]
      calc ((Real.sqrt G.cycleCard : ℝ) : ℂ)⁻¹ * ((sgn (∑ e, s e * p.2 e) : ℝ) : ℂ) *
            (((Real.sqrt G.cycleCard : ℝ) : ℂ)⁻¹ * ((sgn (∑ e, s e * p.2 e) : ℝ) : ℂ)) =
            (((Real.sqrt G.cycleCard : ℝ) : ℂ) * ((Real.sqrt G.cycleCard : ℝ) : ℂ))⁻¹ *
              (((sgn (∑ e, s e * p.2 e) : ℝ) : ℂ) * ((sgn (∑ e, s e * p.2 e) : ℝ) : ℂ)) := by
            ring
        _ = ((G.cycleCard : ℂ))⁻¹ := by rw [hr, hσ, mul_one]
    · simp
  rw [show Star.star (G.fluxVec q s n) ⬝ᵥ G.fluxVec q s n =
    ∑ p, Star.star (G.fluxVec q s n p) * G.fluxVec q s n p from rfl]
  rw [Fintype.sum_prod_type, Finset.sum_eq_single n]
  · simp only [hterm, true_and]
    rw [← Finset.sum_filter, Finset.sum_const, hfib, nsmul_eq_mul]
    exact mul_inv_cancel₀ (by exact_mod_cast hK.ne')
  · intro m _ hm
    apply Finset.sum_eq_zero
    intro x _
    rw [hterm, ite_eq_right (fun h => hm h.1)]
  · simp

/-- Vectors with different charges, occupations or fluxes are orthogonal. -/
theorem fluxVec_orthogonal (q q' : Cfg N) (s s' : Cfg M) (n n' : Cfg N)
    (h : q' ≠ q ∨ n' ≠ n ∨ ¬ ∃ θ : Cfg N, s' = G.gaugeAct θ s) :
    Star.star (G.fluxVec q' s' n') ⬝ᵥ G.fluxVec q s n = 0 := by
  rw [show Star.star (G.fluxVec q' s' n') ⬝ᵥ G.fluxVec q s n =
    ∑ p, Star.star (G.fluxVec q' s' n' p) * G.fluxVec q s n p from rfl]
  rcases h with hq | hn | hθ
  · -- disjoint supports
    apply Finset.sum_eq_zero
    rintro ⟨m, x⟩ -
    simp only [fluxVec]
    split_ifs with h1 h2
    · exfalso
      apply hq
      rw [← h1.2, ← h2.2, ← h1.1, ← h2.1]
    all_goals simp
  · -- disjoint supports
    apply Finset.sum_eq_zero
    rintro ⟨m, x⟩ -
    simp only [fluxVec]
    split_ifs with h1 h2
    · exact absurd (h1.1.symm.trans h2.1) hn
    all_goals simp
  · -- a cycle `z` on which the fluxes differ
    obtain ⟨z, hz, hne⟩ : ∃ z : Cfg M, G.bdry z = 0 ∧ ∑ e, s e * z e ≠ ∑ e, s' e * z e := by
      by_contra hcon
      apply hθ
      rw [← G.flux_eq_iff_gauge]
      intro z hz
      by_contra hne
      exact hcon ⟨z, hz, hne⟩
    have hσ : ((sgn (∑ e, s e * z e) : ℝ) : ℂ) * ((sgn (∑ e, s' e * z e) : ℝ) : ℂ) = -1 := by
      rw [← Complex.ofReal_mul, sgn_mul_sgn_of_ne hne]
      push_cast
      ring
    -- the involution `x ↦ x + z` reverses the sign of each term
    have key : ∀ p : Cfg N × Cfg M,
        Star.star (G.fluxVec q' s' n' (p + (0, z))) * G.fluxVec q s n (p + (0, z)) =
          -(Star.star (G.fluxVec q' s' n' p) * G.fluxVec q s n p) := by
      rintro ⟨m, x⟩
      rw [Prod.mk_add_mk, add_zero, G.fluxVec_add_cycle q' s' n' hz,
        G.fluxVec_add_cycle q s n hz, star_mul', Complex.star_def, Complex.conj_ofReal]
      linear_combination ((starRingEnd ℂ) (G.fluxVec q' s' n' (m, x)) *
        G.fluxVec q s n (m, x)) * hσ
    have hsum := Equiv.sum_comp (Equiv.addRight ((0 : Cfg N), z))
      (fun p => Star.star (G.fluxVec q' s' n' p) * G.fluxVec q s n p)
    simp only [Equiv.coe_addRight, key, Finset.sum_neg_distrib] at hsum
    linear_combination -hsum / 2

/-- **Completeness**: `2^N` charge configurations, `2^β` flux classes and `2^{N-1}` occupations of the
right parity give `2^{N+M}` labels, the dimension of the Hilbert space (connected graph, `N ≥ 1`). -/
theorem card_fluxLabels (hG : G.Connected) (hN : 0 < N) :
    2 ^ N * G.fluxClasses.card * 2 ^ (N - 1) = Fintype.card (Cfg N × Cfg M) := by
  have h1 : G.fluxClasses.card = 2 ^ (M + 1 - N) := G.card_fluxClasses hG hN
  have h2 : Fintype.card (Cfg N × Cfg M) = 2 ^ N * 2 ^ M := by
    simp only [Fintype.card_prod, Fintype.card_fun, ZMod.card, Fintype.card_fin]
  have hle : N - 1 ≤ M := by
    have h3 : 2 ^ (N - 1) ≤ 2 ^ M := by
      rw [← G.card_gaugeOrbit hG hN 0]
      calc _ ≤ Fintype.card (Cfg M) := Finset.card_le_univ _
        _ = 2 ^ M := by simp only [Fintype.card_fun, ZMod.card, Fintype.card_fin]
    exact (Nat.pow_le_pow_iff_right (by norm_num)).1 h3
  rw [h1, h2, mul_assoc, ← pow_add, show M + 1 - N + (N - 1) = M by omega]

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

/-- `Ψ_{n,q,s}` is the normalized image of the basis vector `|n, q⟩` under `V_s`. -/
private lemma fluxVec_eq_VSc (q : Cfg N) (s : Cfg M) (n : Cfg N) :
    G.fluxVec q s n = ((Real.sqrt G.cycleCard : ℝ) : ℂ)⁻¹ •
      ((G.VSmat s).map (fun r => (r : ℂ)) *ᵥ Pi.single (n, q) 1) := by
  funext ⟨m, x⟩
  rw [Pi.smul_apply, VSc_mulVec_apply, smul_eq_mul]
  simp only [fluxVec, Pi.single_apply, Prod.mk.injEq]
  by_cases hm : m = n
  · subst hm
    by_cases hq : m + G.bdry x = q
    · simp [hq]
    · simp [hq]
  · simp [hm]

/-- Applying a matrix to a basis vector picks a column. -/
private lemma mulVec_single_one_apply {ι κ : Type*} [Fintype κ] [DecidableEq κ]
    (A : Matrix ι κ ℂ) (j : κ) (i : ι) : (A *ᵥ Pi.single j 1) i = A i j := by
  simp [mulVec, dotProduct, Pi.single_apply]

/-- **Sector Hamiltonians** (paper: Eqs. `eq:2ddict`, `eq:2dsector`): `H Ψ_{n,q,s}` is the combination
of the `Ψ_{n',q,s}` with the matrix elements of `dΓ(h_{q,s}) - h ∑_i q_i`. -/
theorem ham_fluxVec (J h : ℝ) (q : Cfg N) (s : Cfg M) (n : Cfg N) :
    cplx (G.ham J h) *ᵥ G.fluxVec q s n =
      ∑ n', (((dGamma (G.hqs J h (fun i => sgn (q i)) (fun e => sgn (s e))) -
        (h * ∑ i, sgn (q i)) • 1) n' n : ℝ) : ℂ) • G.fluxVec q s n' := by
  -- `H V_s = V_s H_dual,s`, complexified
  have hint : cplx (G.ham J h) * (G.VSmat s).map (fun r => (r : ℂ)) =
      (G.VSmat s).map (fun r => (r : ℂ)) * cplx (G.dualHamS J h s) := by
    have := congrArg (fun X => X.map Complex.ofRealHom) (G.ham_mul_VSmat J h s)
    simp only [Matrix.map_mul] at this
    exact this
  simp only [fluxVec_eq_VSc]
  rw [mulVec_smul, mulVec_mulVec, hint, ← mulVec_mulVec]
  funext ⟨m, x⟩
  simp only [Pi.smul_apply, Finset.sum_apply, VSc_mulVec_apply, smul_eq_mul]
  simp only [mulVec_single_one_apply, cplx, Matrix.map_apply, dualHamS_eq, blockDiagonal_apply,
    Pi.single_apply, Prod.mk.injEq]
  rw [Finset.sum_eq_single m]
  · by_cases hq : m + G.bdry x = q
    · subst hq
      simp only [and_self, ite_true, mul_one]
      ring
    · simp [hq]
  · intro n' _ hn'
    simp [Ne.symm hn']
  · simp

end FinGraph

end DFL
