import DFL.TwoDMap

/-!
# Planar boxes and tori (paper: Theorem `thm:2dmap`(c),(d))

On a connected graph the flux sectors are the gauge classes of the link signs (`FinGraph.flux_eq_iff_gauge`,
`FinGraph.card_fluxClasses`); a basis of the cycle space `ker ∂ ⊆ (ℤ₂)^E` therefore gives a complete set of
independent flux labels (`IsCycleBasis`, `gauge_iff_flux_basis`, `exists_flux_basis`, `card_cycleBasis`).

* (c) **Planar box** `{0,…,L₁-1} × {0,…,L₂-1}` (`box`): `β = (L₁-1)(L₂-1)` (`box_beta`), the plaquette
  boundaries form a basis of the cycle space (`box_plaq_basis`); hence the flux sectors are labelled by
  the independent plaquette fluxes `W_p ∈ ℤ₂`, and all `2^β` labels occur (`box_card_fluxClasses`).
* (d) **Torus** `(ℤ/L₁) × (ℤ/L₂)`, `L₁, L₂ ≥ 3` (`torus`): `β = L₁L₂ + 1` (`torus_beta`),
  `∑_p ∂p = 0`, i.e. `∏_p W_p = 1` (`torus_sum_plaq`, `torus_prod_wilson`), and for every plaquette `p₀`
  the boundaries `∂p`, `p ≠ p₀`, together with the two non-contractible loops `γ₁` (the horizontal edges of
  the row `y = 0`) and `γ₂` (the vertical edges of the column `x = 0`) form a basis of the cycle space
  (`torus_basis`); all `2^β` labels occur (`torus_card_fluxClasses`).

The sector Hamiltonians in each labelled sector are given by `FinGraph.ham_fluxVec` (Theorem `thm:2dmap`(b)).
-/

open Matrix Finset

noncomputable section

namespace DFL

namespace FinGraph

variable {N M : ℕ}

/-! ### Cycle bases -/

/-- A basis (over `ℤ₂`) of the cycle space `ker ∂`. -/
def IsCycleBasis (G : FinGraph N M) {ι : Type*} (c : ι → Cfg M) : Prop :=
  (∀ i, G.bdry (c i) = 0) ∧ LinearIndependent (ZMod 2) c ∧
    ∀ z : Cfg M, G.bdry z = 0 → z ∈ Submodule.span (ZMod 2) (Set.range c)

/-- The flux of the sign configuration `s` through `z`, `W_z(s) = ∑_{e ∈ z} s_e` in `ℤ₂`
(`(-1)^{W_z(s)} = ∏_{e ∈ z} s_e` in multiplicative notation). -/
def fluxThrough (s z : Cfg M) : ZMod 2 := ∑ e, s e * z e

namespace TwoDLatticesAux

/-- `a + a = 0` in `ZMod 2`. -/
lemma z2_add_self (a : ZMod 2) : a + a = 0 := by
  revert a; decide

/-- `a + b = 0` forces `a = b` in `ZMod 2`. -/
lemma z2_eq_of_add_eq_zero {a b : ZMod 2} (h : a + b = 0) : a = b := by
  revert a b; decide

/-- `z + z = 0` for configurations. -/
lemma cfg_add_self {K : ℕ} (z : Cfg K) : z + z = 0 := by
  funext e
  exact z2_add_self (z e)

/-- `∂` is additive. -/
lemma bdry_add (G : FinGraph N M) (x y : Cfg M) : G.bdry (x + y) = G.bdry x + G.bdry y := by
  funext i
  simp only [bdry, Pi.add_apply, Finset.sum_add_distrib]

/-- The boundary map as a `ZMod 2`-linear map. -/
def bdryLin (G : FinGraph N M) : Cfg M →ₗ[ZMod 2] Cfg N where
  toFun := G.bdry
  map_add' := bdry_add G
  map_smul' := by
    intro c x
    funext i
    simp [bdry, Finset.mul_sum]

/-- `bdryLin` is `∂`. -/
lemma bdryLin_apply (G : FinGraph N M) (x : Cfg M) : bdryLin G x = G.bdry x := rfl

/-- `∂` of a single edge is the sum of its two endpoints. -/
lemma bdry_single (G : FinGraph N M) (e : Fin M) :
    G.bdry (Pi.single e 1) = Pi.single (G.src e) 1 + Pi.single (G.tgt e) 1 := by
  funext i
  have hne := G.src_ne_tgt e
  simp only [bdry, Pi.add_apply]
  rw [Finset.sum_pi_single']
  simp only [Finset.mem_filter, Finset.mem_univ, true_and, Pi.single_apply]
  by_cases h1 : G.src e = i
  · have h2 : ¬ G.tgt e = i := fun h2 => hne (h1.trans h2.symm)
    simp [h1, h2, Ne.symm h2]
  · by_cases h2 : G.tgt e = i
    · simp [h1, h2, Ne.symm h1]
    · simp [h1, h2, Ne.symm h1, Ne.symm h2]

/-- The pairing `z ↦ ∑ e, s e * z e` as a linear functional. -/
def pairLin (s : Cfg M) : Cfg M →ₗ[ZMod 2] ZMod 2 where
  toFun z := ∑ e, s e * z e
  map_add' x y := by simp only [Pi.add_apply, mul_add, Finset.sum_add_distrib]
  map_smul' a x := by
    simp only [Pi.smul_apply, smul_eq_mul, RingHom.id_apply, Finset.mul_sum]
    exact Finset.sum_congr rfl (fun e _ => by ring)

/-- `pairLin (s - s')` is the difference of the fluxes. -/
lemma pairLin_sub_apply (s s' z : Cfg M) :
    pairLin (s - s') z = fluxThrough s z - fluxThrough s' z := by
  simp only [pairLin, LinearMap.coe_mk, AddHom.coe_mk, fluxThrough, Pi.sub_apply, sub_mul,
    Finset.sum_sub_distrib]

/-- Composition of gauge transformations. -/
lemma gaugeAct_gaugeAct (G : FinGraph N M) (θ θ' : Cfg N) (s : Cfg M) :
    G.gaugeAct θ (G.gaugeAct θ' s) = G.gaugeAct (θ + θ') s := by
  funext e
  simp only [gaugeAct, Pi.add_apply]
  ring

/-- The trivial gauge transformation. -/
lemma gaugeAct_zero (G : FinGraph N M) (s : Cfg M) : G.gaugeAct 0 s = s := by
  funext e
  simp [gaugeAct]

/-- Two gauge orbits coincide iff the sign configurations are gauge equivalent. -/
lemma orbit_eq_iff (G : FinGraph N M) (s s' : Cfg M) :
    univ.image (fun θ : Cfg N => G.gaugeAct θ s) = univ.image (fun θ : Cfg N => G.gaugeAct θ s') ↔
      ∃ θ : Cfg N, s' = G.gaugeAct θ s := by
  constructor
  · intro h
    have hs : s' ∈ univ.image (fun θ : Cfg N => G.gaugeAct θ s') :=
      Finset.mem_image.2 ⟨0, Finset.mem_univ _, gaugeAct_zero G s'⟩
    rw [← h] at hs
    obtain ⟨θ, -, hθ⟩ := Finset.mem_image.1 hs
    exact ⟨θ, hθ.symm⟩
  · rintro ⟨θ, rfl⟩
    ext t
    simp only [Finset.mem_image, Finset.mem_univ, true_and]
    constructor
    · rintro ⟨θ', rfl⟩
      refine ⟨θ' + θ, ?_⟩
      rw [gaugeAct_gaugeAct, add_assoc, cfg_add_self, add_zero]
    · rintro ⟨θ', rfl⟩
      exact ⟨θ' + θ, (gaugeAct_gaugeAct G θ' θ s).symm⟩

/-- Two maps with the same fibres have images of the same cardinality. -/
lemma card_image_eq_of_iff {α β γ : Type*} [Fintype α] [DecidableEq β] [DecidableEq γ]
    (f : α → β) (g : α → γ) (h : ∀ a a', f a = f a' ↔ g a = g a') :
    (univ.image f).card = (univ.image g).card := by
  refine Finset.card_bij (fun b hb => g (Finset.mem_image.1 hb).choose) ?_ ?_ ?_
  · intro b hb
    exact Finset.mem_image_of_mem g (Finset.mem_univ _)
  · intro b₁ hb₁ b₂ hb₂ heq
    have e₁ := (Finset.mem_image.1 hb₁).choose_spec.2
    have e₂ := (Finset.mem_image.1 hb₂).choose_spec.2
    rw [← e₁, ← e₂]
    exact (h _ _).2 heq
  · intro c hc
    obtain ⟨a, -, rfl⟩ := Finset.mem_image.1 hc
    have hb : f a ∈ univ.image f := Finset.mem_image_of_mem f (Finset.mem_univ a)
    refine ⟨f a, hb, ?_⟩
    exact (h _ _).1 (Finset.mem_image.1 hb).choose_spec.2

end TwoDLatticesAux

open TwoDLatticesAux

/-- Two sign configurations are gauge equivalent iff they have the same fluxes through the cycles of a
cycle basis. -/
theorem gauge_iff_flux_basis (G : FinGraph N M) {ι : Type*} (c : ι → Cfg M) (hc : G.IsCycleBasis c)
    (s s' : Cfg M) :
    (∃ θ : Cfg N, s' = G.gaugeAct θ s) ↔ ∀ i, fluxThrough s (c i) = fluxThrough s' (c i) := by
  rw [← G.flux_eq_iff_gauge s s']
  constructor
  · intro h i
    exact h (c i) (hc.1 i)
  · intro h z hz
    have hle : Submodule.span (ZMod 2) (Set.range c) ≤ LinearMap.ker (pairLin (s - s')) := by
      rw [Submodule.span_le]
      rintro _ ⟨i, rfl⟩
      rw [SetLike.mem_coe, LinearMap.mem_ker, pairLin_sub_apply, h i, sub_self]
    have h2 := hle (hc.2.2 z hz)
    rw [LinearMap.mem_ker, pairLin_sub_apply] at h2
    exact sub_eq_zero.1 h2

/-- The fluxes through the cycles of a cycle basis are independent: every assignment occurs. -/
theorem exists_flux_basis (G : FinGraph N M) {ι : Type*} (c : ι → Cfg M) (hc : G.IsCycleBasis c)
    (φ : ι → ZMod 2) : ∃ s : Cfg M, ∀ i, fluxThrough s (c i) = φ i := by
  obtain ⟨-, hli, -⟩ := hc
  obtain ⟨g, hg⟩ := LinearMap.exists_extend ((Module.Basis.span hli).constr (ZMod 2) φ)
  refine ⟨fun e => g (fun j => if e = j then 1 else 0), fun i => ?_⟩
  have h1 : g (c i) = φ i := by
    have h2 : g (c i) =
        (g.comp (Submodule.span (ZMod 2) (Set.range c)).subtype) (Module.Basis.span hli i) := by
      simp [Module.Basis.span_apply]
    rw [h2, hg, Module.Basis.constr_basis]
  rw [← h1, LinearMap.pi_apply_eq_sum_univ g (c i)]
  simp only [fluxThrough, smul_eq_mul]
  exact Finset.sum_congr rfl (fun e _ => mul_comm _ _)

/-- A cycle basis of a connected graph has `β = M - N + 1` elements. -/
theorem card_cycleBasis (G : FinGraph N M) (hG : G.Connected) (hN : 0 < N) {ι : Type*} [Fintype ι]
    (c : ι → Cfg M) (hc : G.IsCycleBasis c) : Fintype.card ι = M + 1 - N := by
  classical
  have key : ∀ s s' : Cfg M,
      univ.image (fun θ : Cfg N => G.gaugeAct θ s) = univ.image (fun θ : Cfg N => G.gaugeAct θ s') ↔
        (fun i => fluxThrough s (c i)) = (fun i => fluxThrough s' (c i)) := by
    intro s s'
    rw [orbit_eq_iff, G.gauge_iff_flux_basis c hc, funext_iff]
  have h1 := card_image_eq_of_iff _ _ key
  have h2 : univ.image (fun s : Cfg M => fun i => fluxThrough s (c i)) = univ := by
    apply Finset.eq_univ_of_forall
    intro φ
    obtain ⟨s, hs⟩ := G.exists_flux_basis c hc φ
    exact Finset.mem_image.2 ⟨s, Finset.mem_univ _, funext hs⟩
  rw [h2, Finset.card_univ, Fintype.card_fun, ZMod.card, G.card_fluxClasses hG hN] at h1
  exact (Nat.pow_right_injective (le_refl 2) h1).symm

namespace TwoDLatticesAux

/-- The cycle space `ker ∂` of a connected graph has dimension `β = M - N + 1`. -/
lemma finrank_cycleSpace (G : FinGraph N M) (hG : G.Connected) (hN : 0 < N) :
    Module.finrank (ZMod 2) (LinearMap.ker (bdryLin G)) = M + 1 - N := by
  set K := LinearMap.ker (bdryLin G)
  let b := Module.finBasis (ZMod 2) K
  have hc : G.IsCycleBasis (fun i => (b i : Cfg M)) := by
    refine ⟨fun i => LinearMap.mem_ker.1 (b i).2, ?_, fun z hz => ?_⟩
    · exact b.linearIndependent.map' K.subtype (Submodule.ker_subtype K)
    · have h1 : (⟨z, LinearMap.mem_ker.2 hz⟩ : K) ∈ Submodule.span (ZMod 2) (Set.range b) :=
        b.mem_span _
      have h2 := Submodule.apply_mem_span_image_of_mem_span K.subtype h1
      rw [← Set.range_comp] at h2
      exact h2
  have := G.card_cycleBasis hG hN _ hc
  simpa using this

/-- An independent family of cycles of the right size `β = M - N + 1` is a cycle basis. -/
lemma isCycleBasis_of_card (G : FinGraph N M) (hG : G.Connected) (hN : 0 < N) {ι : Type*}
    [Fintype ι] (c : ι → Cfg M) (hcyc : ∀ i, G.bdry (c i) = 0)
    (hli : LinearIndependent (ZMod 2) c) (hcard : Fintype.card ι = M + 1 - N) :
    G.IsCycleBasis c := by
  refine ⟨hcyc, hli, fun z hz => ?_⟩
  have hle : Submodule.span (ZMod 2) (Set.range c) ≤ LinearMap.ker (bdryLin G) := by
    rw [Submodule.span_le]
    rintro _ ⟨i, rfl⟩
    exact LinearMap.mem_ker.2 (hcyc i)
  have heq := Submodule.eq_of_le_of_finrank_eq hle
    (by rw [finrank_span_eq_card hli, hcard, finrank_cycleSpace G hG hN])
  rw [heq]
  exact LinearMap.mem_ker.2 hz

end TwoDLatticesAux

/-! ### (c) The planar box -/

/-- The vertex `(x, y)` of the box. -/
def boxVert {L₁ L₂ : ℕ} (x : Fin L₁) (y : Fin L₂) : Fin (L₁ * L₂) := finProdFinEquiv (x, y)

/-- The number of edges of the box. -/
def boxM (L₁ L₂ : ℕ) : ℕ := (L₁ - 1) * L₂ + L₁ * (L₂ - 1)

/-- Edge labels of the box: horizontal edges `(x, y) — (x+1, y)` (`x < L₁ - 1`), then vertical edges
`(x, y) — (x, y+1)` (`y < L₂ - 1`). -/
def boxEdge (L₁ L₂ : ℕ) : (Fin (L₁ - 1) × Fin L₂) ⊕ (Fin L₁ × Fin (L₂ - 1)) ≃ Fin (boxM L₁ L₂) :=
  (Equiv.sumCongr finProdFinEquiv finProdFinEquiv).trans finSumFinEquiv

/-- Initial vertex of a box edge. -/
def boxSrc (L₁ L₂ : ℕ) : (Fin (L₁ - 1) × Fin L₂) ⊕ (Fin L₁ × Fin (L₂ - 1)) → Fin (L₁ * L₂)
  | Sum.inl e => boxVert ⟨e.1.val, by have := e.1.isLt; omega⟩ e.2
  | Sum.inr e => boxVert e.1 ⟨e.2.val, by have := e.2.isLt; omega⟩

/-- Final vertex of a box edge. -/
def boxTgt (L₁ L₂ : ℕ) : (Fin (L₁ - 1) × Fin L₂) ⊕ (Fin L₁ × Fin (L₂ - 1)) → Fin (L₁ * L₂)
  | Sum.inl e => boxVert ⟨e.1.val + 1, by have := e.1.isLt; omega⟩ e.2
  | Sum.inr e => boxVert e.1 ⟨e.2.val + 1, by have := e.2.isLt; omega⟩

theorem boxSrc_ne_boxTgt (L₁ L₂ : ℕ) (e : (Fin (L₁ - 1) × Fin L₂) ⊕ (Fin L₁ × Fin (L₂ - 1))) :
    boxSrc L₁ L₂ e ≠ boxTgt L₁ L₂ e := by
  rcases e with _e | _e <;> intro h <;>
    have h' := congrArg (fun v => (finProdFinEquiv.symm v)) h <;>
    simp [boxSrc, boxTgt, boxVert, Fin.ext_iff] at h'

/-- **The planar box** `{0,…,L₁-1} × {0,…,L₂-1}` with nearest-neighbour edges. -/
def box (L₁ L₂ : ℕ) : FinGraph (L₁ * L₂) (boxM L₁ L₂) where
  src e := boxSrc L₁ L₂ ((boxEdge L₁ L₂).symm e)
  tgt e := boxTgt L₁ L₂ ((boxEdge L₁ L₂).symm e)
  src_ne_tgt e := boxSrc_ne_boxTgt L₁ L₂ ((boxEdge L₁ L₂).symm e)

/-- The boundary `∂p` of the plaquette `p = (x, y)` (`x < L₁ - 1`, `y < L₂ - 1`): the edges
`(x,y)—(x+1,y)`, `(x,y+1)—(x+1,y+1)`, `(x,y)—(x,y+1)` and `(x+1,y)—(x+1,y+1)`. -/
def boxPlaq (L₁ L₂ : ℕ) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) : Cfg (boxM L₁ L₂) := fun e =>
  if e = boxEdge L₁ L₂ (Sum.inl (p.1, ⟨p.2.val, by have := p.2.isLt; omega⟩)) ∨
     e = boxEdge L₁ L₂ (Sum.inl (p.1, ⟨p.2.val + 1, by have := p.2.isLt; omega⟩)) ∨
     e = boxEdge L₁ L₂ (Sum.inr (⟨p.1.val, by have := p.1.isLt; omega⟩, p.2)) ∨
     e = boxEdge L₁ L₂ (Sum.inr (⟨p.1.val + 1, by have := p.1.isLt; omega⟩, p.2))
  then 1 else 0

namespace TwoDLatticesAux

/-- Adjacency is symmetric. -/
lemma adj_symm (G : FinGraph N M) {i j : Fin N} (h : G.Adj i j) : G.Adj j i := by
  obtain ⟨e, h⟩ := h
  exact ⟨e, h.elim (fun h => Or.inr ⟨h.2, h.1⟩) (fun h => Or.inl ⟨h.2, h.1⟩)⟩

/-- Connectedness by paths is symmetric. -/
lemma reflTransGen_adj_symm (G : FinGraph N M) {i j : Fin N}
    (h : Relation.ReflTransGen G.Adj i j) : Relation.ReflTransGen G.Adj j i := by
  induction h with
  | refl => exact Relation.ReflTransGen.refl
  | tail _ hbc ih => exact Relation.ReflTransGen.head (adj_symm G hbc) ih

/-- Initial vertex of a box edge. -/
lemma box_src (L₁ L₂ : ℕ) (e : (Fin (L₁ - 1) × Fin L₂) ⊕ (Fin L₁ × Fin (L₂ - 1))) :
    (box L₁ L₂).src (boxEdge L₁ L₂ e) = boxSrc L₁ L₂ e := by
  simp [box]

/-- Final vertex of a box edge. -/
lemma box_tgt (L₁ L₂ : ℕ) (e : (Fin (L₁ - 1) × Fin L₂) ⊕ (Fin L₁ × Fin (L₂ - 1))) :
    (box L₁ L₂).tgt (boxEdge L₁ L₂ e) = boxTgt L₁ L₂ e := by
  simp [box]

end TwoDLatticesAux

/-- The box is connected. -/
theorem box_connected (L₁ L₂ : ℕ) : (box L₁ L₂).Connected := by
  have hadj : ∀ e, (box L₁ L₂).Adj (boxSrc L₁ L₂ e) (boxTgt L₁ L₂ e) :=
    fun e => ⟨boxEdge L₁ L₂ e, Or.inl ⟨box_src L₁ L₂ e, box_tgt L₁ L₂ e⟩⟩
  -- moving left along a row
  have hx : ∀ (h₁ : 0 < L₁) (k : ℕ) (hk : k < L₁) (y : Fin L₂),
      Relation.ReflTransGen (box L₁ L₂).Adj (boxVert ⟨k, hk⟩ y) (boxVert ⟨0, h₁⟩ y) := by
    intro h₁ k
    induction k with
    | zero => intro hk y; exact Relation.ReflTransGen.refl
    | succ k ih =>
      intro hk y
      exact Relation.ReflTransGen.head (adj_symm _ (hadj (Sum.inl (⟨k, by omega⟩, y))))
        (ih (by omega) y)
  -- moving down along the column `x = 0`
  have hy : ∀ (h₁ : 0 < L₁) (h₂ : 0 < L₂) (k : ℕ) (hk : k < L₂),
      Relation.ReflTransGen (box L₁ L₂).Adj (boxVert ⟨0, h₁⟩ ⟨k, hk⟩)
        (boxVert ⟨0, h₁⟩ ⟨0, h₂⟩) := by
    intro h₁ h₂ k
    induction k with
    | zero => intro hk; exact Relation.ReflTransGen.refl
    | succ k ih =>
      intro hk
      exact Relation.ReflTransGen.head (adj_symm _ (hadj (Sum.inr (⟨0, h₁⟩, ⟨k, by omega⟩))))
        (ih (by omega))
  have hroot : ∀ (h₁ : 0 < L₁) (h₂ : 0 < L₂) (i : Fin (L₁ * L₂)),
      Relation.ReflTransGen (box L₁ L₂).Adj i (boxVert ⟨0, h₁⟩ ⟨0, h₂⟩) := by
    intro h₁ h₂ i
    have hi : i = boxVert (finProdFinEquiv.symm i).1 (finProdFinEquiv.symm i).2 :=
      (Equiv.apply_symm_apply finProdFinEquiv i).symm
    rw [hi]
    exact (hx h₁ _ _ _).trans (hy h₁ h₂ _ _)
  intro i j
  have h₁ : 0 < L₁ := Fin.pos (finProdFinEquiv.symm i).1
  have h₂ : 0 < L₂ := Fin.pos (finProdFinEquiv.symm i).2
  exact (hroot h₁ h₂ i).trans (reflTransGen_adj_symm _ (hroot h₁ h₂ j))

/-- **First Betti number of the box** (paper: Theorem `thm:2dmap`(c)): `β = (L₁-1)(L₂-1)`. -/
theorem box_beta (L₁ L₂ : ℕ) (h₁ : 0 < L₁) (h₂ : 0 < L₂) :
    boxM L₁ L₂ + 1 - L₁ * L₂ = (L₁ - 1) * (L₂ - 1) := by
  obtain ⟨a, rfl⟩ : ∃ a, L₁ = a + 1 := ⟨L₁ - 1, by omega⟩
  obtain ⟨b, rfl⟩ : ∃ b, L₂ = b + 1 := ⟨L₂ - 1, by omega⟩
  simp only [boxM, Nat.add_sub_cancel]
  have h : a * (b + 1) + (a + 1) * b + 1 = (a + 1) * (b + 1) + a * b := by ring
  rw [h, Nat.add_sub_cancel_left]

namespace TwoDLatticesAux

/-- The indicator of four distinct edges is the sum of the four unit vectors. -/
lemma ind4_eq {K : ℕ} (a b c d : Fin K) (hab : a ≠ b) (hac : a ≠ c) (had : a ≠ d) (hbc : b ≠ c)
    (hbd : b ≠ d) (hcd : c ≠ d) :
    (fun e => if e = a ∨ e = b ∨ e = c ∨ e = d then (1 : ZMod 2) else 0) =
      Pi.single a 1 + Pi.single b 1 + Pi.single c 1 + Pi.single d 1 := by
  funext e
  simp only [Pi.add_apply, Pi.single_apply]
  by_cases ha : e = a
  · subst ha
    simp [hab, hac, had]
  · by_cases hb : e = b
    · subst hb
      simp [ha, hbc, hbd]
    · by_cases hc : e = c
      · subst hc
        simp [ha, hb, hcd]
      · by_cases hd : e = d
        · subst hd
          simp [ha, hb, hc]
        · simp [ha, hb, hc, hd]

/-- `∂` of the boundary of a square `a — b`, `c — d`, `a — c`, `b — d` vanishes. -/
lemma bdry_square (G : FinGraph N M) {e₁ e₂ e₃ e₄ : Fin M} {a b c d : Fin N}
    (h₁ : G.src e₁ = a ∧ G.tgt e₁ = b) (h₂ : G.src e₂ = c ∧ G.tgt e₂ = d)
    (h₃ : G.src e₃ = a ∧ G.tgt e₃ = c) (h₄ : G.src e₄ = b ∧ G.tgt e₄ = d) :
    G.bdry (Pi.single e₁ 1 + Pi.single e₂ 1 + Pi.single e₃ 1 + Pi.single e₄ 1) = 0 := by
  simp only [bdry_add, bdry_single, h₁.1, h₁.2, h₂.1, h₂.2, h₃.1, h₃.2, h₄.1, h₄.2]
  have h : ∀ x y z w : Cfg N, x + y + (z + w) + (x + z) + (y + w) = 0 := by
    intro x y z w
    calc x + y + (z + w) + (x + z) + (y + w) = (x + x) + (y + y) + (z + z) + (w + w) := by abel
      _ = 0 := by simp only [cfg_add_self, add_zero]
  exact h _ _ _ _

/-- Evaluating a combination of unit vectors along an injective family. -/
lemma sum_smul_single_apply_self {α : Type*} [Fintype α] {K : ℕ} (f : α → Fin K)
    (hf : Function.Injective f) (c : α → ZMod 2) (q : α) :
    (∑ p, c p • (Pi.single (f p) 1 : Cfg K)) (f q) = c q := by
  classical
  rw [Finset.sum_apply, Finset.sum_eq_single q]
  · simp
  · intro p _ hpq
    simp [hf.ne hpq.symm]
  · simp

/-- A combination of unit vectors vanishes at `e` if all coefficients of `e` vanish. -/
lemma sum_smul_single_apply_eq_zero {α : Type*} [Fintype α] {K : ℕ} (f : α → Fin K)
    (c : α → ZMod 2) (e : Fin K) (h : ∀ p, f p = e → c p = 0) :
    (∑ p, c p • (Pi.single (f p) 1 : Cfg K)) e = 0 := by
  rw [Finset.sum_apply]
  apply Finset.sum_eq_zero
  intro p _
  by_cases hp : f p = e
  · simp [h p hp]
  · simp [Ne.symm hp]

/-- The bottom edge `(x,y)—(x+1,y)` of a box plaquette. -/
def boxBot (L₁ L₂ : ℕ) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) : Fin (boxM L₁ L₂) :=
  boxEdge L₁ L₂ (Sum.inl (p.1, ⟨p.2.val, by have := p.2.isLt; omega⟩))

/-- The top edge `(x,y+1)—(x+1,y+1)` of a box plaquette. -/
def boxTop (L₁ L₂ : ℕ) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) : Fin (boxM L₁ L₂) :=
  boxEdge L₁ L₂ (Sum.inl (p.1, ⟨p.2.val + 1, by have := p.2.isLt; omega⟩))

/-- The left edge `(x,y)—(x,y+1)` of a box plaquette. -/
def boxLeft (L₁ L₂ : ℕ) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) : Fin (boxM L₁ L₂) :=
  boxEdge L₁ L₂ (Sum.inr (⟨p.1.val, by have := p.1.isLt; omega⟩, p.2))

/-- The right edge `(x+1,y)—(x+1,y+1)` of a box plaquette. -/
def boxRight (L₁ L₂ : ℕ) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) : Fin (boxM L₁ L₂) :=
  boxEdge L₁ L₂ (Sum.inr (⟨p.1.val + 1, by have := p.1.isLt; omega⟩, p.2))

/-- A box plaquette is the sum of its four edges. -/
lemma boxPlaq_eq (L₁ L₂ : ℕ) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) :
    boxPlaq L₁ L₂ p = Pi.single (boxBot L₁ L₂ p) 1 + Pi.single (boxTop L₁ L₂ p) 1 +
      Pi.single (boxLeft L₁ L₂ p) 1 + Pi.single (boxRight L₁ L₂ p) 1 :=
  ind4_eq _ _ _ _ (by simp) (by simp) (by simp) (by simp) (by simp) (by simp)

/-- Box plaquettes are cycles. -/
lemma box_plaq_cycle (L₁ L₂ : ℕ) (p : Fin (L₁ - 1) × Fin (L₂ - 1)) :
    (box L₁ L₂).bdry (boxPlaq L₁ L₂ p) = 0 := by
  rw [boxPlaq_eq]
  exact bdry_square (box L₁ L₂) ⟨box_src L₁ L₂ _, box_tgt L₁ L₂ _⟩
    ⟨box_src L₁ L₂ _, box_tgt L₁ L₂ _⟩ ⟨box_src L₁ L₂ _, box_tgt L₁ L₂ _⟩
    ⟨box_src L₁ L₂ _, box_tgt L₁ L₂ _⟩

/-- Box plaquettes are linearly independent (peel off the rows from the bottom, using the bottom
edges). -/
lemma box_plaq_linIndep (L₁ L₂ : ℕ) : LinearIndependent (ZMod 2) (boxPlaq L₁ L₂) := by
  rw [Fintype.linearIndependent_iff]
  intro g hg
  have hsum : ∀ e, (∑ p, g p • (Pi.single (boxBot L₁ L₂ p) 1 : Cfg (boxM L₁ L₂))) e +
      (∑ p, g p • (Pi.single (boxTop L₁ L₂ p) 1 : Cfg (boxM L₁ L₂))) e +
      (∑ p, g p • (Pi.single (boxLeft L₁ L₂ p) 1 : Cfg (boxM L₁ L₂))) e +
      (∑ p, g p • (Pi.single (boxRight L₁ L₂ p) 1 : Cfg (boxM L₁ L₂))) e = 0 := by
    intro e
    have := congrFun hg e
    simp only [boxPlaq_eq, smul_add, Finset.sum_add_distrib, Pi.add_apply] at this
    exact this
  have hinj : Function.Injective (boxBot L₁ L₂) := by
    intro p q h
    simp only [boxBot, EmbeddingLike.apply_eq_iff_eq, Sum.inl.injEq, Prod.mk.injEq,
      Fin.mk.injEq] at h
    exact Prod.ext h.1 (Fin.ext h.2)
  have key : ∀ n, ∀ p : Fin (L₁ - 1) × Fin (L₂ - 1), p.2.val = n → g p = 0 := by
    intro n
    induction n using Nat.strong_induction_on with
    | _ n ih =>
      intro p hp
      have htop : (∑ p', g p' • (Pi.single (boxTop L₁ L₂ p') 1 : Cfg (boxM L₁ L₂)))
          (boxBot L₁ L₂ p) = 0 := by
        apply sum_smul_single_apply_eq_zero
        intro p' hp'
        simp only [boxTop, boxBot, EmbeddingLike.apply_eq_iff_eq, Sum.inl.injEq, Prod.mk.injEq,
          Fin.mk.injEq] at hp'
        exact ih p'.2.val (by omega) p' rfl
      have hleft : (∑ p', g p' • (Pi.single (boxLeft L₁ L₂ p') 1 : Cfg (boxM L₁ L₂)))
          (boxBot L₁ L₂ p) = 0 := by
        apply sum_smul_single_apply_eq_zero
        intro p' hp'
        simp [boxLeft, boxBot] at hp'
      have hright : (∑ p', g p' • (Pi.single (boxRight L₁ L₂ p') 1 : Cfg (boxM L₁ L₂)))
          (boxBot L₁ L₂ p) = 0 := by
        apply sum_smul_single_apply_eq_zero
        intro p' hp'
        simp [boxRight, boxBot] at hp'
      have h := hsum (boxBot L₁ L₂ p)
      rw [sum_smul_single_apply_self _ hinj g p, htop, hleft, hright] at h
      simpa using h
  exact fun p => key _ p rfl

end TwoDLatticesAux

/-- **The plaquettes form a cycle basis** (paper: Theorem `thm:2dmap`(c)). -/
theorem box_plaq_basis (L₁ L₂ : ℕ) : (box L₁ L₂).IsCycleBasis (boxPlaq L₁ L₂) := by
  rcases Nat.eq_zero_or_pos L₁ with h₁ | h₁
  · -- no vertices, no edges
    refine ⟨box_plaq_cycle L₁ L₂, box_plaq_linIndep L₁ L₂, fun z _ => ?_⟩
    have hM : boxM L₁ L₂ = 0 := by simp [boxM, h₁]
    have hz : z = 0 := funext fun e => absurd e.isLt (by omega)
    rw [hz]
    exact Submodule.zero_mem _
  rcases Nat.eq_zero_or_pos L₂ with h₂ | h₂
  · refine ⟨box_plaq_cycle L₁ L₂, box_plaq_linIndep L₁ L₂, fun z _ => ?_⟩
    have hM : boxM L₁ L₂ = 0 := by simp [boxM, h₂]
    have hz : z = 0 := funext fun e => absurd e.isLt (by omega)
    rw [hz]
    exact Submodule.zero_mem _
  exact isCycleBasis_of_card (box L₁ L₂) (box_connected L₁ L₂) (Nat.mul_pos h₁ h₂) _
    (box_plaq_cycle L₁ L₂) (box_plaq_linIndep L₁ L₂)
    (by rw [Fintype.card_prod, Fintype.card_fin, Fintype.card_fin, box_beta L₁ L₂ h₁ h₂])

/-- **All `2^β` flux labels occur** (paper: Theorem `thm:2dmap`(c)): the flux sectors of the box are
labelled by independent plaquette fluxes. -/
theorem box_card_fluxClasses (L₁ L₂ : ℕ) (h₁ : 0 < L₁) (h₂ : 0 < L₂) :
    (box L₁ L₂).fluxClasses.card = 2 ^ ((L₁ - 1) * (L₂ - 1)) := by
  rw [← box_beta L₁ L₂ h₁ h₂]
  exact (box L₁ L₂).card_fluxClasses (box_connected L₁ L₂) (Nat.mul_pos h₁ h₂)

/-! ### (d) The torus -/

/-- The number of edges of the torus. -/
def torusM (L₁ L₂ : ℕ) : ℕ := L₁ * L₂ + L₁ * L₂

/-- Edge labels of the torus: horizontal edges `(x, y) — (x+1, y)`, then vertical edges
`(x, y) — (x, y+1)` (coordinates mod `L₁`, `L₂`). -/
def torusEdge (L₁ L₂ : ℕ) : (Fin L₁ × Fin L₂) ⊕ (Fin L₁ × Fin L₂) ≃ Fin (torusM L₁ L₂) :=
  (Equiv.sumCongr finProdFinEquiv finProdFinEquiv).trans finSumFinEquiv

/-- `x + 1 mod L`. -/
def finSucc {L : ℕ} (x : Fin L) : Fin L := ⟨(x.val + 1) % L, Nat.mod_lt _ (by have := x.isLt; omega)⟩

/-- Initial vertex of a torus edge. -/
def torusSrc (L₁ L₂ : ℕ) : (Fin L₁ × Fin L₂) ⊕ (Fin L₁ × Fin L₂) → Fin (L₁ * L₂)
  | Sum.inl e => boxVert e.1 e.2
  | Sum.inr e => boxVert e.1 e.2

/-- Final vertex of a torus edge. -/
def torusTgt (L₁ L₂ : ℕ) : (Fin L₁ × Fin L₂) ⊕ (Fin L₁ × Fin L₂) → Fin (L₁ * L₂)
  | Sum.inl e => boxVert (finSucc e.1) e.2
  | Sum.inr e => boxVert e.1 (finSucc e.2)

theorem finSucc_ne {L : ℕ} (hL : 2 ≤ L) (x : Fin L) : finSucc x ≠ x := by
  intro h
  have h' := congrArg Fin.val h
  simp only [finSucc] at h'
  have hx := x.isLt
  by_cases hx1 : x.val + 1 < L
  · rw [Nat.mod_eq_of_lt hx1] at h'; omega
  · have : x.val + 1 = L := by omega
    rw [this, Nat.mod_self] at h'; omega

theorem torusSrc_ne_torusTgt (L₁ L₂ : ℕ) (h₁ : 2 ≤ L₁) (h₂ : 2 ≤ L₂)
    (e : (Fin L₁ × Fin L₂) ⊕ (Fin L₁ × Fin L₂)) : torusSrc L₁ L₂ e ≠ torusTgt L₁ L₂ e := by
  rcases e with e | e <;> intro h <;>
    have h' := congrArg (fun v => (finProdFinEquiv.symm v)) h
  · simp only [torusSrc, torusTgt, boxVert, Equiv.symm_apply_apply, Prod.mk.injEq] at h'
    exact finSucc_ne h₁ e.1 h'.1.symm
  · simp only [torusSrc, torusTgt, boxVert, Equiv.symm_apply_apply, Prod.mk.injEq] at h'
    exact finSucc_ne h₂ e.2 h'.2.symm

/-- **The torus** `(ℤ/L₁) × (ℤ/L₂)` with nearest-neighbour edges (`L₁, L₂ ≥ 3`). -/
def torus (L₁ L₂ : ℕ) (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) : FinGraph (L₁ * L₂) (torusM L₁ L₂) where
  src e := torusSrc L₁ L₂ ((torusEdge L₁ L₂).symm e)
  tgt e := torusTgt L₁ L₂ ((torusEdge L₁ L₂).symm e)
  src_ne_tgt e := torusSrc_ne_torusTgt L₁ L₂ (by omega) (by omega) ((torusEdge L₁ L₂).symm e)

/-- The boundary `∂p` of the plaquette `p = (x, y)`: the edges `(x,y)—(x+1,y)`, `(x,y+1)—(x+1,y+1)`,
`(x,y)—(x,y+1)`, `(x+1,y)—(x+1,y+1)`. -/
def torusPlaq (L₁ L₂ : ℕ) (p : Fin L₁ × Fin L₂) : Cfg (torusM L₁ L₂) := fun e =>
  if e = torusEdge L₁ L₂ (Sum.inl (p.1, p.2)) ∨ e = torusEdge L₁ L₂ (Sum.inl (p.1, finSucc p.2)) ∨
     e = torusEdge L₁ L₂ (Sum.inr (p.1, p.2)) ∨ e = torusEdge L₁ L₂ (Sum.inr (finSucc p.1, p.2))
  then 1 else 0

/-- The non-contractible loop `γ₁`: the horizontal edges of the row `y = 0`. -/
def torusGamma1 (L₁ L₂ : ℕ) (h₂ : 0 < L₂) : Cfg (torusM L₁ L₂) := fun e =>
  if ∃ x : Fin L₁, e = torusEdge L₁ L₂ (Sum.inl (x, ⟨0, h₂⟩)) then 1 else 0

/-- The non-contractible loop `γ₂`: the vertical edges of the column `x = 0`. -/
def torusGamma2 (L₁ L₂ : ℕ) (h₁ : 0 < L₁) : Cfg (torusM L₁ L₂) := fun e =>
  if ∃ y : Fin L₂, e = torusEdge L₁ L₂ (Sum.inr (⟨0, h₁⟩, y)) then 1 else 0

namespace TwoDLatticesAux

/-- Value of `finSucc` away from the wrap-around. -/
lemma finSucc_val_of_lt {L : ℕ} (x : Fin L) (h : x.val + 1 < L) : (finSucc x).val = x.val + 1 := by
  simp [finSucc, Nat.mod_eq_of_lt h]

/-- Value of `finSucc` at the wrap-around. -/
lemma finSucc_val_of_eq {L : ℕ} (x : Fin L) (h : x.val + 1 = L) : (finSucc x).val = 0 := by
  simp [finSucc, h]

/-- `finSucc` is injective. -/
lemma finSucc_injective {L : ℕ} : Function.Injective (finSucc : Fin L → Fin L) := by
  intro x y h
  have h' := congrArg Fin.val h
  have hx := x.isLt
  have hy := y.isLt
  apply Fin.ext
  by_cases hx1 : x.val + 1 < L <;> by_cases hy1 : y.val + 1 < L
  · rw [finSucc_val_of_lt x hx1, finSucc_val_of_lt y hy1] at h'; omega
  · rw [finSucc_val_of_lt x hx1, finSucc_val_of_eq y (by omega)] at h'; omega
  · rw [finSucc_val_of_eq x (by omega), finSucc_val_of_lt y hy1] at h'; omega
  · omega

/-- `finSucc` is bijective. -/
lemma finSucc_bijective {L : ℕ} : Function.Bijective (finSucc : Fin L → Fin L) :=
  Finite.injective_iff_bijective.1 finSucc_injective

/-- The shift `(x, y) ↦ (x, y + 1)` is bijective. -/
lemma shift2_bijective {L₁ L₂ : ℕ} :
    Function.Bijective (fun p : Fin L₁ × Fin L₂ => (p.1, finSucc p.2)) :=
  Function.Bijective.prodMap Function.bijective_id finSucc_bijective

/-- The shift `(x, y) ↦ (x + 1, y)` is bijective. -/
lemma shift1_bijective {L₁ L₂ : ℕ} :
    Function.Bijective (fun p : Fin L₁ × Fin L₂ => (finSucc p.1, p.2)) :=
  Function.Bijective.prodMap finSucc_bijective Function.bijective_id

/-- The horizontal torus edge `(x,y)—(x+1,y)`. -/
def torusH (L₁ L₂ : ℕ) (p : Fin L₁ × Fin L₂) : Fin (torusM L₁ L₂) := torusEdge L₁ L₂ (Sum.inl p)

/-- The vertical torus edge `(x,y)—(x,y+1)`. -/
def torusV (L₁ L₂ : ℕ) (p : Fin L₁ × Fin L₂) : Fin (torusM L₁ L₂) := torusEdge L₁ L₂ (Sum.inr p)

/-- Horizontal edges are labelled injectively. -/
lemma torusH_injective {L₁ L₂ : ℕ} : Function.Injective (torusH L₁ L₂) :=
  fun _ _ h => Sum.inl_injective ((torusEdge L₁ L₂).injective h)

/-- Vertical edges are labelled injectively. -/
lemma torusV_injective {L₁ L₂ : ℕ} : Function.Injective (torusV L₁ L₂) :=
  fun _ _ h => Sum.inr_injective ((torusEdge L₁ L₂).injective h)

/-- A vertical edge is not horizontal. -/
lemma torusV_ne_torusH {L₁ L₂ : ℕ} (p q : Fin L₁ × Fin L₂) : torusV L₁ L₂ p ≠ torusH L₁ L₂ q :=
  fun h => Sum.inr_ne_inl ((torusEdge L₁ L₂).injective h)

/-- A horizontal edge is not vertical. -/
lemma torusH_ne_torusV {L₁ L₂ : ℕ} (p q : Fin L₁ × Fin L₂) : torusH L₁ L₂ p ≠ torusV L₁ L₂ q :=
  fun h => Sum.inl_ne_inr ((torusEdge L₁ L₂).injective h)

/-- Initial vertex of a torus edge. -/
lemma torus_src {L₁ L₂ : ℕ} (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) (e : (Fin L₁ × Fin L₂) ⊕ (Fin L₁ × Fin L₂)) :
    (torus L₁ L₂ h₁ h₂).src (torusEdge L₁ L₂ e) = torusSrc L₁ L₂ e := by
  simp [torus]

/-- Final vertex of a torus edge. -/
lemma torus_tgt {L₁ L₂ : ℕ} (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) (e : (Fin L₁ × Fin L₂) ⊕ (Fin L₁ × Fin L₂)) :
    (torus L₁ L₂ h₁ h₂).tgt (torusEdge L₁ L₂ e) = torusTgt L₁ L₂ e := by
  simp [torus]

/-- A torus plaquette is the sum of its four edges. -/
lemma torusPlaq_eq {L₁ L₂ : ℕ} (h₁ : 2 ≤ L₁) (h₂ : 2 ≤ L₂) (p : Fin L₁ × Fin L₂) :
    torusPlaq L₁ L₂ p = Pi.single (torusH L₁ L₂ p) 1 + Pi.single (torusH L₁ L₂ (p.1, finSucc p.2)) 1 +
      Pi.single (torusV L₁ L₂ p) 1 + Pi.single (torusV L₁ L₂ (finSucc p.1, p.2)) 1 := by
  have hne1 : p.1 ≠ finSucc p.1 := fun h => finSucc_ne h₁ p.1 h.symm
  have hne2 : p.2 ≠ finSucc p.2 := fun h => finSucc_ne h₂ p.2 h.symm
  exact ind4_eq _ _ _ _ (by simp [Prod.ext_iff, hne2]) (by simp) (by simp) (by simp) (by simp)
    (by simp [Prod.ext_iff, hne1])

/-- Torus plaquettes are cycles. -/
lemma torus_plaq_cycle {L₁ L₂ : ℕ} (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) (p : Fin L₁ × Fin L₂) :
    (torus L₁ L₂ h₁ h₂).bdry (torusPlaq L₁ L₂ p) = 0 := by
  rw [torusPlaq_eq (by omega) (by omega)]
  exact bdry_square (torus L₁ L₂ h₁ h₂) ⟨torus_src h₁ h₂ _, torus_tgt h₁ h₂ _⟩
    ⟨torus_src h₁ h₂ _, torus_tgt h₁ h₂ _⟩ ⟨torus_src h₁ h₂ _, torus_tgt h₁ h₂ _⟩
    ⟨torus_src h₁ h₂ _, torus_tgt h₁ h₂ _⟩

/-- A sum of unit vectors along an injective family, evaluated on the family. -/
lemma sum_single_apply_self {α : Type*} [Fintype α] {K : ℕ} (f : α → Fin K)
    (hf : Function.Injective f) (q : α) : (∑ x, (Pi.single (f x) 1 : Cfg K)) (f q) = 1 := by
  have := sum_smul_single_apply_self f hf (fun _ => 1) q
  simp only [one_smul] at this
  exact this

/-- A sum of unit vectors along a family, evaluated off the family. -/
lemma sum_single_apply_eq_zero {α : Type*} [Fintype α] {K : ℕ} (f : α → Fin K) (e : Fin K)
    (he : ∀ x, f x ≠ e) : (∑ x, (Pi.single (f x) 1 : Cfg K)) e = 0 := by
  have := sum_smul_single_apply_eq_zero f (fun _ => 1) e (fun x hx => absurd hx (he x))
  simp only [one_smul] at this
  exact this

/-- `γ₁` is the sum of the horizontal edges of the row `y = 0`. -/
lemma torusGamma1_eq {L₁ L₂ : ℕ} (h : 0 < L₂) :
    torusGamma1 L₁ L₂ h = ∑ x : Fin L₁, Pi.single (torusH L₁ L₂ (x, ⟨0, h⟩)) 1 := by
  have hinj : Function.Injective (fun x : Fin L₁ => torusH L₁ L₂ (x, ⟨0, h⟩)) :=
    fun _ _ hxy => congrArg Prod.fst (torusH_injective hxy)
  funext e
  unfold torusGamma1
  split_ifs with he
  · obtain ⟨x, rfl⟩ := he
    exact (sum_single_apply_self (fun x : Fin L₁ => torusH L₁ L₂ (x, ⟨0, h⟩)) hinj x).symm
  · exact (sum_single_apply_eq_zero (fun x : Fin L₁ => torusH L₁ L₂ (x, ⟨0, h⟩)) e
      (fun x hx => he ⟨x, hx.symm⟩)).symm

/-- `γ₂` is the sum of the vertical edges of the column `x = 0`. -/
lemma torusGamma2_eq {L₁ L₂ : ℕ} (h : 0 < L₁) :
    torusGamma2 L₁ L₂ h = ∑ y : Fin L₂, Pi.single (torusV L₁ L₂ (⟨0, h⟩, y)) 1 := by
  have hinj : Function.Injective (fun y : Fin L₂ => torusV L₁ L₂ (⟨0, h⟩, y)) :=
    fun _ _ hxy => congrArg Prod.snd (torusV_injective hxy)
  funext e
  unfold torusGamma2
  split_ifs with he
  · obtain ⟨y, rfl⟩ := he
    exact (sum_single_apply_self (fun y : Fin L₂ => torusV L₁ L₂ (⟨0, h⟩, y)) hinj y).symm
  · exact (sum_single_apply_eq_zero (fun y : Fin L₂ => torusV L₁ L₂ (⟨0, h⟩, y)) e
      (fun y hy => he ⟨y, hy.symm⟩)).symm

/-- `γ₁` is a cycle. -/
lemma torus_gamma1_cycle {L₁ L₂ : ℕ} (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) (h : 0 < L₂) :
    (torus L₁ L₂ h₁ h₂).bdry (torusGamma1 L₁ L₂ h) = 0 := by
  rw [torusGamma1_eq, ← bdryLin_apply, map_sum]
  simp only [bdryLin_apply, bdry_single, torusH, torus_src, torus_tgt, torusSrc, torusTgt]
  rw [Finset.sum_add_distrib]
  have hs : ∑ x : Fin L₁, (Pi.single (boxVert (finSucc x) (⟨0, h⟩ : Fin L₂)) 1 : Cfg (L₁ * L₂)) =
      ∑ x : Fin L₁, Pi.single (boxVert x (⟨0, h⟩ : Fin L₂)) 1 :=
    finSucc_bijective.sum_comp (fun x => (Pi.single (boxVert x (⟨0, h⟩ : Fin L₂)) 1 : Cfg (L₁ * L₂)))
  rw [hs, cfg_add_self]

/-- `γ₂` is a cycle. -/
lemma torus_gamma2_cycle {L₁ L₂ : ℕ} (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) (h : 0 < L₁) :
    (torus L₁ L₂ h₁ h₂).bdry (torusGamma2 L₁ L₂ h) = 0 := by
  rw [torusGamma2_eq, ← bdryLin_apply, map_sum]
  simp only [bdryLin_apply, bdry_single, torusV, torus_src, torus_tgt, torusSrc, torusTgt]
  rw [Finset.sum_add_distrib]
  have hs : ∑ y : Fin L₂, (Pi.single (boxVert (⟨0, h⟩ : Fin L₁) (finSucc y)) 1 : Cfg (L₁ * L₂)) =
      ∑ y : Fin L₂, Pi.single (boxVert (⟨0, h⟩ : Fin L₁) y) 1 :=
    finSucc_bijective.sum_comp (fun y => (Pi.single (boxVert (⟨0, h⟩ : Fin L₁) y) 1 : Cfg (L₁ * L₂)))
  rw [hs, cfg_add_self]

end TwoDLatticesAux

/-- The torus is connected. -/
theorem torus_connected (L₁ L₂ : ℕ) (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) : (torus L₁ L₂ h₁ h₂).Connected := by
  have hadj : ∀ e, (torus L₁ L₂ h₁ h₂).Adj (torusSrc L₁ L₂ e) (torusTgt L₁ L₂ e) :=
    fun e => ⟨torusEdge L₁ L₂ e, Or.inl ⟨torus_src h₁ h₂ e, torus_tgt h₁ h₂ e⟩⟩
  have hsucc : ∀ {L : ℕ} (k : ℕ) (hk : k + 1 < L), finSucc (⟨k, by omega⟩ : Fin L) = ⟨k + 1, hk⟩ :=
    fun k hk => Fin.ext (finSucc_val_of_lt _ hk)
  -- moving left along a row
  have hx : ∀ (k : ℕ) (hk : k < L₁) (y : Fin L₂),
      Relation.ReflTransGen (torus L₁ L₂ h₁ h₂).Adj (boxVert ⟨k, hk⟩ y)
        (boxVert ⟨0, by omega⟩ y) := by
    intro k
    induction k with
    | zero => intro hk y; exact Relation.ReflTransGen.refl
    | succ k ih =>
      intro hk y
      have h := hadj (Sum.inl (⟨k, by omega⟩, y))
      simp only [torusSrc, torusTgt, hsucc k hk] at h
      exact Relation.ReflTransGen.head (adj_symm _ h) (ih (by omega) y)
  -- moving down along the column `x = 0`
  have hy : ∀ (k : ℕ) (hk : k < L₂),
      Relation.ReflTransGen (torus L₁ L₂ h₁ h₂).Adj (boxVert ⟨0, by omega⟩ ⟨k, hk⟩)
        (boxVert ⟨0, by omega⟩ ⟨0, by omega⟩) := by
    intro k
    induction k with
    | zero => intro hk; exact Relation.ReflTransGen.refl
    | succ k ih =>
      intro hk
      have h := hadj (Sum.inr (⟨0, by omega⟩, ⟨k, by omega⟩))
      simp only [torusSrc, torusTgt, hsucc k hk] at h
      exact Relation.ReflTransGen.head (adj_symm _ h) (ih (by omega))
  have hroot : ∀ i : Fin (L₁ * L₂), Relation.ReflTransGen (torus L₁ L₂ h₁ h₂).Adj i
      (boxVert ⟨0, by omega⟩ ⟨0, by omega⟩) := by
    intro i
    have hi : i = boxVert (finProdFinEquiv.symm i).1 (finProdFinEquiv.symm i).2 :=
      (Equiv.apply_symm_apply finProdFinEquiv i).symm
    rw [hi]
    exact (hx _ _ _).trans (hy _ _)
  intro i j
  exact (hroot i).trans (reflTransGen_adj_symm _ (hroot j))

/-- **First Betti number of the torus** (paper: Theorem `thm:2dmap`(d)): `β = L₁L₂ + 1`. -/
theorem torus_beta (L₁ L₂ : ℕ) (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) :
    torusM L₁ L₂ + 1 - L₁ * L₂ = L₁ * L₂ + 1 := by
  have _ := h₁
  have _ := h₂
  rw [torusM, add_assoc, Nat.add_sub_cancel_left]

/-- `∑_p ∂p = 0`: every edge lies on exactly two plaquettes. -/
theorem torus_sum_plaq (L₁ L₂ : ℕ) (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) :
    ∑ p : Fin L₁ × Fin L₂, torusPlaq L₁ L₂ p = 0 := by
  simp only [torusPlaq_eq (by omega : 2 ≤ L₁) (by omega : 2 ≤ L₂), Finset.sum_add_distrib]
  have hH : ∑ p : Fin L₁ × Fin L₂, (Pi.single (torusH L₁ L₂ (p.1, finSucc p.2)) 1 : Cfg (torusM L₁ L₂)) =
      ∑ p, Pi.single (torusH L₁ L₂ p) 1 :=
    shift2_bijective.sum_comp (fun p => (Pi.single (torusH L₁ L₂ p) 1 : Cfg (torusM L₁ L₂)))
  have hV : ∑ p : Fin L₁ × Fin L₂, (Pi.single (torusV L₁ L₂ (finSucc p.1, p.2)) 1 : Cfg (torusM L₁ L₂)) =
      ∑ p, Pi.single (torusV L₁ L₂ p) 1 :=
    shift1_bijective.sum_comp (fun p => (Pi.single (torusV L₁ L₂ p) 1 : Cfg (torusM L₁ L₂)))
  rw [hH, hV]
  simp only [cfg_add_self, zero_add]

namespace TwoDLatticesAux

open scoped Kronecker

/-- The translation `x ↦ x + z` of edge labels, as a matrix. -/
def transl (z : Cfg M) : Matrix (Cfg M) (Cfg M) ℝ :=
  fun x' x => if x' = x + z then 1 else 0

/-- Translations compose additively. -/
lemma transl_mul (z z' : Cfg M) : transl z * transl z' = transl (z + z') := by
  ext x' x
  simp only [transl, Matrix.mul_apply, mul_ite, mul_one, mul_zero]
  rw [Finset.sum_eq_single (x + z')]
  · simp only [ite_true, add_assoc, add_comm z' z]
  · intro y _ hy
    simp [hy]
  · simp

/-- The translation by `0` is the identity. -/
lemma transl_zero : transl (0 : Cfg M) = 1 := by
  ext x' x
  simp [transl, Matrix.one_apply]

/-- The Wilson loop is `1 ⊗ (translation by z)`. -/
lemma wilsonLoop_eq_transl (G : FinGraph N M) (z : Cfg M) :
    G.wilsonLoop z = (1 : Matrix (Cfg N) (Cfg N) ℝ) ⊗ₖ transl z := rfl

/-- `W_z W_{z'} = W_{z + z'}`. -/
lemma wilsonLoop_mul_wilsonLoop (G : FinGraph N M) (z z' : Cfg M) :
    G.wilsonLoop z * G.wilsonLoop z' = G.wilsonLoop (z + z') := by
  rw [wilsonLoop_eq_transl, wilsonLoop_eq_transl, wilsonLoop_eq_transl, ← Matrix.mul_kronecker_mul,
    Matrix.one_mul, transl_mul]

/-- `W_0 = 1`. -/
lemma wilsonLoop_zero (G : FinGraph N M) : G.wilsonLoop 0 = 1 := by
  rw [wilsonLoop_eq_transl, transl_zero, Matrix.one_kronecker_one]

/-- The product of the Wilson loops of a list is the Wilson loop of the sum. -/
lemma list_prod_wilsonLoop (G : FinGraph N M) {α : Type*} (f : α → Cfg M) (l : List α) :
    (l.map (fun a => G.wilsonLoop (f a))).prod = G.wilsonLoop (l.map f).sum := by
  induction l with
  | nil => simp [wilsonLoop_zero]
  | cons a l ih => rw [List.map_cons, List.prod_cons, ih, List.map_cons, List.sum_cons,
      wilsonLoop_mul_wilsonLoop]

end TwoDLatticesAux

/-- **`∏_p W_p = 1`** (paper: Theorem `thm:2dmap`(d)). -/
theorem torus_prod_wilson (L₁ L₂ : ℕ) (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) :
    ((univ : Finset (Fin L₁ × Fin L₂)).toList.map
      (fun p => (torus L₁ L₂ h₁ h₂).wilsonLoop (torusPlaq L₁ L₂ p))).prod = 1 := by
  rw [list_prod_wilsonLoop, Finset.sum_map_toList, torus_sum_plaq L₁ L₂ h₁ h₂, wilsonLoop_zero]

namespace TwoDLatticesAux

/-- A combination of torus plaquettes, evaluated on the horizontal edge `(x, y+1)—(x+1, y+1)`. -/
lemma torus_plaqSum_apply_H {L₁ L₂ : ℕ} (h₁ : 2 ≤ L₁) (h₂ : 2 ≤ L₂) (c : Fin L₁ × Fin L₂ → ZMod 2)
    (x : Fin L₁) (y : Fin L₂) :
    (∑ p, c p • torusPlaq L₁ L₂ p) (torusH L₁ L₂ (x, finSucc y)) = c (x, finSucc y) + c (x, y) := by
  simp only [torusPlaq_eq h₁ h₂, smul_add, Finset.sum_add_distrib, Pi.add_apply]
  have e1 : (∑ p, c p • (Pi.single (torusH L₁ L₂ p) 1 : Cfg (torusM L₁ L₂)))
      (torusH L₁ L₂ (x, finSucc y)) = c (x, finSucc y) :=
    sum_smul_single_apply_self _ torusH_injective c _
  have e2 : (∑ p, c p • (Pi.single (torusH L₁ L₂ (p.1, finSucc p.2)) 1 : Cfg (torusM L₁ L₂)))
      (torusH L₁ L₂ (x, finSucc y)) = c (x, y) :=
    sum_smul_single_apply_self (fun p : Fin L₁ × Fin L₂ => torusH L₁ L₂ (p.1, finSucc p.2))
      (torusH_injective.comp shift2_bijective.1) c (x, y)
  have e3 : (∑ p, c p • (Pi.single (torusV L₁ L₂ p) 1 : Cfg (torusM L₁ L₂)))
      (torusH L₁ L₂ (x, finSucc y)) = 0 :=
    sum_smul_single_apply_eq_zero _ c _ (fun p hp => absurd hp (torusV_ne_torusH _ _))
  have e4 : (∑ p, c p • (Pi.single (torusV L₁ L₂ (finSucc p.1, p.2)) 1 : Cfg (torusM L₁ L₂)))
      (torusH L₁ L₂ (x, finSucc y)) = 0 :=
    sum_smul_single_apply_eq_zero _ c _ (fun p hp => absurd hp (torusV_ne_torusH _ _))
  rw [e1, e2, e3, e4, add_zero, add_zero]

/-- A combination of torus plaquettes, evaluated on the vertical edge `(x+1, y)—(x+1, y+1)`. -/
lemma torus_plaqSum_apply_V {L₁ L₂ : ℕ} (h₁ : 2 ≤ L₁) (h₂ : 2 ≤ L₂) (c : Fin L₁ × Fin L₂ → ZMod 2)
    (x : Fin L₁) (y : Fin L₂) :
    (∑ p, c p • torusPlaq L₁ L₂ p) (torusV L₁ L₂ (finSucc x, y)) = c (finSucc x, y) + c (x, y) := by
  simp only [torusPlaq_eq h₁ h₂, smul_add, Finset.sum_add_distrib, Pi.add_apply]
  have e1 : (∑ p, c p • (Pi.single (torusH L₁ L₂ p) 1 : Cfg (torusM L₁ L₂)))
      (torusV L₁ L₂ (finSucc x, y)) = 0 :=
    sum_smul_single_apply_eq_zero _ c _ (fun p hp => absurd hp (torusH_ne_torusV _ _))
  have e2 : (∑ p, c p • (Pi.single (torusH L₁ L₂ (p.1, finSucc p.2)) 1 : Cfg (torusM L₁ L₂)))
      (torusV L₁ L₂ (finSucc x, y)) = 0 :=
    sum_smul_single_apply_eq_zero _ c _ (fun p hp => absurd hp (torusH_ne_torusV _ _))
  have e3 : (∑ p, c p • (Pi.single (torusV L₁ L₂ p) 1 : Cfg (torusM L₁ L₂)))
      (torusV L₁ L₂ (finSucc x, y)) = c (finSucc x, y) :=
    sum_smul_single_apply_self _ torusV_injective c _
  have e4 : (∑ p, c p • (Pi.single (torusV L₁ L₂ (finSucc p.1, p.2)) 1 : Cfg (torusM L₁ L₂)))
      (torusV L₁ L₂ (finSucc x, y)) = c (x, y) :=
    sum_smul_single_apply_self (fun p : Fin L₁ × Fin L₂ => torusV L₁ L₂ (finSucc p.1, p.2))
      (torusV_injective.comp shift1_bijective.1) c (x, y)
  rw [e1, e2, e3, e4]
  simp only [zero_add]

/-- `γ₁` on a horizontal edge. -/
lemma gamma1_H {L₁ L₂ : ℕ} (h : 0 < L₂) (x : Fin L₁) (y : Fin L₂) :
    torusGamma1 L₁ L₂ h (torusH L₁ L₂ (x, y)) = if y = ⟨0, h⟩ then 1 else 0 := by
  simp [torusGamma1, torusH]

/-- `γ₂` vanishes on horizontal edges. -/
lemma gamma2_H {L₁ L₂ : ℕ} (h : 0 < L₁) (x : Fin L₁) (y : Fin L₂) :
    torusGamma2 L₁ L₂ h (torusH L₁ L₂ (x, y)) = 0 := by
  simp [torusGamma2, torusH]

/-- `γ₁` vanishes on vertical edges. -/
lemma gamma1_V {L₁ L₂ : ℕ} (h : 0 < L₂) (x : Fin L₁) (y : Fin L₂) :
    torusGamma1 L₁ L₂ h (torusV L₁ L₂ (x, y)) = 0 := by
  simp [torusGamma1, torusV]

/-- `γ₂` on a vertical edge. -/
lemma gamma2_V {L₁ L₂ : ℕ} (h : 0 < L₁) (x : Fin L₁) (y : Fin L₂) :
    torusGamma2 L₁ L₂ h (torusV L₁ L₂ (x, y)) = if x = ⟨0, h⟩ then 1 else 0 := by
  simp [torusGamma2, torusV]

/-- A function on `ℤ/L` that is invariant under the shift `y ↦ y + 1` up to a defect `a` at the seam
is constant, and the defect vanishes. -/
lemma chain_const {L : ℕ} (hL : 0 < L) (f : Fin L → ZMod 2) (a : ZMod 2)
    (h : ∀ y : Fin L, f (finSucc y) + f y + (if finSucc y = ⟨0, hL⟩ then a else 0) = 0) :
    a = 0 ∧ ∀ y, f y = f ⟨0, hL⟩ := by
  have hstep : ∀ k (hk : k + 1 < L), f ⟨k + 1, hk⟩ = f ⟨k, by omega⟩ := by
    intro k hk
    have h1 := h ⟨k, by omega⟩
    have hs : finSucc (⟨k, by omega⟩ : Fin L) = ⟨k + 1, hk⟩ := Fin.ext (finSucc_val_of_lt _ hk)
    have hne : (⟨k + 1, hk⟩ : Fin L) ≠ ⟨0, hL⟩ := by simp
    rw [hs, ite_eq_right hne, add_zero] at h1
    exact z2_eq_of_add_eq_zero h1
  have hall : ∀ k (hk : k < L), f ⟨k, hk⟩ = f ⟨0, hL⟩ := by
    intro k
    induction k with
    | zero => intro hk; rfl
    | succ k ih => intro hk; rw [hstep k hk]; exact ih _
  refine ⟨?_, fun y => hall y.val y.isLt⟩
  have h1 := h ⟨L - 1, by omega⟩
  have hs : finSucc (⟨L - 1, by omega⟩ : Fin L) = ⟨0, hL⟩ :=
    Fin.ext (finSucc_val_of_eq _ (by simp only; omega))
  rw [hs, ite_eq_left rfl, hall (L - 1) (by omega), z2_add_self, zero_add] at h1
  exact h1

/-- The plaquettes `p ≠ p₀` and the loops `γ₁`, `γ₂` are linearly independent. -/
lemma torus_linIndep {L₁ L₂ : ℕ} (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) (p₀ : Fin L₁ × Fin L₂) (hL₁ : 0 < L₁)
    (hL₂ : 0 < L₂) :
    LinearIndependent (ZMod 2)
      (Sum.elim (fun p : {p : Fin L₁ × Fin L₂ // p ≠ p₀} => torusPlaq L₁ L₂ p.1)
        (fun k : Fin 2 => if k = 0 then torusGamma1 L₁ L₂ hL₂ else torusGamma2 L₁ L₂ hL₁)) := by
  rw [Fintype.linearIndependent_iff]
  intro g hg
  obtain ⟨c, hc0, hcp⟩ : ∃ c : Fin L₁ × Fin L₂ → ZMod 2,
      c p₀ = 0 ∧ ∀ p : {p : Fin L₁ × Fin L₂ // p ≠ p₀}, c p.1 = g (Sum.inl p) :=
    ⟨fun p => if hp : p = p₀ then 0 else g (Sum.inl ⟨p, hp⟩), dite_eq_left rfl,
      fun p => dite_eq_right p.2⟩
  have hsub : ∑ p, c p • torusPlaq L₁ L₂ p =
      ∑ p : {p : Fin L₁ × Fin L₂ // p ≠ p₀}, g (Sum.inl p) • torusPlaq L₁ L₂ p.1 := by
    rw [Fintype.sum_eq_add_sum_subtype_ne _ p₀, hc0, zero_smul, zero_add]
    exact Finset.sum_congr rfl (fun p _ => by rw [hcp p])
  have hF : ∑ p, c p • torusPlaq L₁ L₂ p + g (Sum.inr 0) • torusGamma1 L₁ L₂ hL₂ +
      g (Sum.inr 1) • torusGamma2 L₁ L₂ hL₁ = 0 := by
    rw [hsub, ← hg, Fintype.sum_sum_type, Fin.sum_univ_two, add_assoc]
    simp
  have hH : ∀ x y, c (x, finSucc y) + c (x, y) +
      (if finSucc y = ⟨0, hL₂⟩ then g (Sum.inr 0) else 0) = 0 := by
    intro x y
    have := congrFun hF (torusH L₁ L₂ (x, finSucc y))
    simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, Pi.zero_apply,
      torus_plaqSum_apply_H (by omega : 2 ≤ L₁) (by omega : 2 ≤ L₂), gamma1_H, gamma2_H, mul_zero,
      add_zero, mul_ite, mul_one] at this
    exact this
  have hV : ∀ x y, c (finSucc x, y) + c (x, y) +
      (if finSucc x = ⟨0, hL₁⟩ then g (Sum.inr 1) else 0) = 0 := by
    intro x y
    have := congrFun hF (torusV L₁ L₂ (finSucc x, y))
    simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, Pi.zero_apply,
      torus_plaqSum_apply_V (by omega : 2 ≤ L₁) (by omega : 2 ≤ L₂), gamma1_V, gamma2_V, mul_zero,
      add_zero, mul_ite, mul_one] at this
    exact this
  obtain ⟨ha, -⟩ := chain_const hL₂ (fun y => c (⟨0, hL₁⟩, y)) _ (hH ⟨0, hL₁⟩)
  have hcol : ∀ x y, c (x, y) = c (x, ⟨0, hL₂⟩) :=
    fun x => (chain_const hL₂ (fun y => c (x, y)) _ (hH x)).2
  obtain ⟨hb, hrow⟩ := chain_const hL₁ (fun x => c (x, ⟨0, hL₂⟩)) _ (fun x => hV x ⟨0, hL₂⟩)
  have hconst : ∀ p : Fin L₁ × Fin L₂, c p = c (⟨0, hL₁⟩, ⟨0, hL₂⟩) :=
    fun p => (hcol p.1 p.2).trans (hrow p.1)
  have hzero : ∀ p, c p = 0 := fun p => (hconst p).trans ((hconst p₀).symm.trans hc0)
  rintro (p | k)
  · rw [← hcp p]
    exact hzero p.1
  · fin_cases k
    · exact ha
    · exact hb

end TwoDLatticesAux

/-- **A cycle basis of the torus** (paper: Theorem `thm:2dmap`(d)): the plaquettes `p ≠ p₀` and the two
non-contractible loops `γ₁`, `γ₂`. -/
theorem torus_basis (L₁ L₂ : ℕ) (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) (p₀ : Fin L₁ × Fin L₂) :
    (torus L₁ L₂ h₁ h₂).IsCycleBasis
      (Sum.elim (fun p : {p : Fin L₁ × Fin L₂ // p ≠ p₀} => torusPlaq L₁ L₂ p.1)
        (fun k : Fin 2 => if k = 0 then torusGamma1 L₁ L₂ (by omega) else torusGamma2 L₁ L₂ (by omega))) := by
  refine isCycleBasis_of_card (torus L₁ L₂ h₁ h₂) (torus_connected L₁ L₂ h₁ h₂)
    (Nat.mul_pos (by omega) (by omega)) _ ?_ (torus_linIndep h₁ h₂ p₀ (by omega) (by omega)) ?_
  · rintro (p | k)
    · exact torus_plaq_cycle h₁ h₂ p.1
    · by_cases hk : k = 0
      · simp only [Sum.elim_inr, ite_eq_left hk]
        exact torus_gamma1_cycle h₁ h₂ _
      · simp only [Sum.elim_inr, ite_eq_right hk]
        exact torus_gamma2_cycle h₁ h₂ _
  · have hcard := Fintype.card_congr (Equiv.optionSubtypeNe p₀)
    rw [Fintype.card_option, Fintype.card_prod, Fintype.card_fin, Fintype.card_fin] at hcard
    rw [Fintype.card_sum, Fintype.card_fin, torus_beta L₁ L₂ h₁ h₂]
    omega

/-- **All `2^β` flux labels occur on the torus** (paper: Theorem `thm:2dmap`(d)). -/
theorem torus_card_fluxClasses (L₁ L₂ : ℕ) (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) :
    (torus L₁ L₂ h₁ h₂).fluxClasses.card = 2 ^ (L₁ * L₂ + 1) := by
  rw [← torus_beta L₁ L₂ h₁ h₂]
  exact (torus L₁ L₂ h₁ h₂).card_fluxClasses (torus_connected L₁ L₂ h₁ h₂)
    (Nat.mul_pos (by omega) (by omega))

end FinGraph

end DFL
