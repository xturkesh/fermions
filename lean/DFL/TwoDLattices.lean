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
lemma adj_symm (G : FinGraph N M) : Symmetric G.Adj := by
  rintro i j ⟨e, h⟩
  exact ⟨e, h.elim (fun h => Or.inr ⟨h.2, h.1⟩) (fun h => Or.inl ⟨h.2, h.1⟩)⟩

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
    have hi : i = boxVert (finProdFinEquiv.symm i).1 (finProdFinEquiv.symm i).2 := by
      simp [boxVert]
    rw [hi]
    exact (hx h₁ _ _ _).trans (hy h₁ h₂ _ _)
  intro i j
  have h₁ : 0 < L₁ := Fin.pos (finProdFinEquiv.symm i).1
  have h₂ : 0 < L₂ := Fin.pos (finProdFinEquiv.symm i).2
  exact (hroot h₁ h₂ i).trans (Relation.ReflTransGen.symmetric (adj_symm _) (hroot h₁ h₂ j))

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
  calc _ = (Pi.single a 1 + Pi.single a 1) + (Pi.single b 1 + Pi.single b 1) +
        (Pi.single c 1 + Pi.single c 1) + (Pi.single d 1 + Pi.single d 1) := by abel
    _ = 0 := by simp only [cfg_add_self, add_zero]

/-- Evaluating a combination of unit vectors along an injective family. -/
lemma sum_smul_single_apply_self {α : Type*} [Fintype α] {K : ℕ} (f : α → Fin K)
    (hf : Function.Injective f) (c : α → ZMod 2) (q : α) :
    (∑ p, c p • (Pi.single (f p) 1 : Cfg K)) (f q) = c q := by
  classical
  rw [Finset.sum_apply, Finset.sum_eq_single q]
  · simp
  · intro p _ hpq
    simp [Pi.single_apply, hf.ne hpq.symm]
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
  · simp [Pi.single_apply, Ne.symm hp]

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
  ind4_eq _ _ _ _ (by simp [boxBot, boxTop]) (by simp [boxBot, boxLeft])
    (by simp [boxBot, boxRight]) (by simp [boxTop, boxLeft]) (by simp [boxTop, boxRight])
    (by simp [boxLeft, boxRight])

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
      have h := hsum (boxBot L₁ L₂ p)
      rw [sum_smul_single_apply_self _ hinj g p, sum_smul_single_apply_eq_zero _ g _ ?_,
        sum_smul_single_apply_eq_zero _ g _ ?_, sum_smul_single_apply_eq_zero _ g _ ?_] at h
      · simpa using h
      · intro p' hp'
        simp only [boxTop, boxBot, EmbeddingLike.apply_eq_iff_eq, Sum.inl.injEq, Prod.mk.injEq,
          Fin.mk.injEq] at hp'
        exact ih p'.2.val (by omega) p' rfl
      · intro p' hp'
        simp [boxLeft, boxBot] at hp'
      · intro p' hp'
        simp [boxRight, boxBot] at hp'
  exact fun p => key _ p rfl

end TwoDLatticesAux

/-- **The plaquettes form a cycle basis** (paper: Theorem `thm:2dmap`(c)). -/
theorem box_plaq_basis (L₁ L₂ : ℕ) : (box L₁ L₂).IsCycleBasis (boxPlaq L₁ L₂) := by
  rcases Nat.eq_zero_or_pos L₁ with h₁ | h₁
  · -- no vertices, no edges
    refine ⟨box_plaq_cycle L₁ L₂, box_plaq_linIndep L₁ L₂, fun z _ => ?_⟩
    have hM : boxM L₁ L₂ = 0 := by simp [boxM, h₁]
    have hz : z = 0 := funext fun e => absurd e.isLt (by rw [hM]; exact Nat.not_lt_zero _)
    rw [hz]
    exact Submodule.zero_mem _
  rcases Nat.eq_zero_or_pos L₂ with h₂ | h₂
  · refine ⟨box_plaq_cycle L₁ L₂, box_plaq_linIndep L₁ L₂, fun z _ => ?_⟩
    have hM : boxM L₁ L₂ = 0 := by simp [boxM, h₂]
    have hz : z = 0 := funext fun e => absurd e.isLt (by rw [hM]; exact Nat.not_lt_zero _)
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

/-- The torus is connected. -/
theorem torus_connected (L₁ L₂ : ℕ) (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) : (torus L₁ L₂ h₁ h₂).Connected := by
  sorry

/-- **First Betti number of the torus** (paper: Theorem `thm:2dmap`(d)): `β = L₁L₂ + 1`. -/
theorem torus_beta (L₁ L₂ : ℕ) (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) :
    torusM L₁ L₂ + 1 - L₁ * L₂ = L₁ * L₂ + 1 := by
  sorry

/-- `∑_p ∂p = 0`: every edge lies on exactly two plaquettes. -/
theorem torus_sum_plaq (L₁ L₂ : ℕ) (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) :
    ∑ p : Fin L₁ × Fin L₂, torusPlaq L₁ L₂ p = 0 := by
  sorry

/-- **`∏_p W_p = 1`** (paper: Theorem `thm:2dmap`(d)). -/
theorem torus_prod_wilson (L₁ L₂ : ℕ) (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) :
    ((univ : Finset (Fin L₁ × Fin L₂)).toList.map
      (fun p => (torus L₁ L₂ h₁ h₂).wilsonLoop (torusPlaq L₁ L₂ p))).prod = 1 := by
  sorry

/-- **A cycle basis of the torus** (paper: Theorem `thm:2dmap`(d)): the plaquettes `p ≠ p₀` and the two
non-contractible loops `γ₁`, `γ₂`. -/
theorem torus_basis (L₁ L₂ : ℕ) (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) (p₀ : Fin L₁ × Fin L₂) :
    (torus L₁ L₂ h₁ h₂).IsCycleBasis
      (Sum.elim (fun p : {p : Fin L₁ × Fin L₂ // p ≠ p₀} => torusPlaq L₁ L₂ p.1)
        (fun k : Fin 2 => if k = 0 then torusGamma1 L₁ L₂ (by omega) else torusGamma2 L₁ L₂ (by omega))) := by
  sorry

/-- **All `2^β` flux labels occur on the torus** (paper: Theorem `thm:2dmap`(d)). -/
theorem torus_card_fluxClasses (L₁ L₂ : ℕ) (h₁ : 3 ≤ L₁) (h₂ : 3 ≤ L₂) :
    (torus L₁ L₂ h₁ h₂).fluxClasses.card = 2 ^ (L₁ * L₂ + 1) := by
  sorry

end FinGraph

end DFL
