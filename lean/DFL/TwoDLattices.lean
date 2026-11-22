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

/-- Two sign configurations are gauge equivalent iff they have the same fluxes through the cycles of a
cycle basis. -/
theorem gauge_iff_flux_basis (G : FinGraph N M) {ι : Type*} (c : ι → Cfg M) (hc : G.IsCycleBasis c)
    (s s' : Cfg M) :
    (∃ θ : Cfg N, s' = G.gaugeAct θ s) ↔ ∀ i, fluxThrough s (c i) = fluxThrough s' (c i) := by
  sorry

/-- The fluxes through the cycles of a cycle basis are independent: every assignment occurs. -/
theorem exists_flux_basis (G : FinGraph N M) {ι : Type*} (c : ι → Cfg M) (hc : G.IsCycleBasis c)
    (φ : ι → ZMod 2) : ∃ s : Cfg M, ∀ i, fluxThrough s (c i) = φ i := by
  sorry

/-- A cycle basis of a connected graph has `β = M - N + 1` elements. -/
theorem card_cycleBasis (G : FinGraph N M) (hG : G.Connected) (hN : 0 < N) {ι : Type*} [Fintype ι]
    (c : ι → Cfg M) (hc : G.IsCycleBasis c) : Fintype.card ι = M + 1 - N := by
  sorry

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

/-- The box is connected. -/
theorem box_connected (L₁ L₂ : ℕ) : (box L₁ L₂).Connected := by
  sorry

/-- **First Betti number of the box** (paper: Theorem `thm:2dmap`(c)): `β = (L₁-1)(L₂-1)`. -/
theorem box_beta (L₁ L₂ : ℕ) (h₁ : 0 < L₁) (h₂ : 0 < L₂) :
    boxM L₁ L₂ + 1 - L₁ * L₂ = (L₁ - 1) * (L₂ - 1) := by
  sorry

/-- **The plaquettes form a cycle basis** (paper: Theorem `thm:2dmap`(c)). -/
theorem box_plaq_basis (L₁ L₂ : ℕ) : (box L₁ L₂).IsCycleBasis (boxPlaq L₁ L₂) := by
  sorry

/-- **All `2^β` flux labels occur** (paper: Theorem `thm:2dmap`(c)): the flux sectors of the box are
labelled by independent plaquette fluxes. -/
theorem box_card_fluxClasses (L₁ L₂ : ℕ) (h₁ : 0 < L₁) (h₂ : 0 < L₂) :
    (box L₁ L₂).fluxClasses.card = 2 ^ ((L₁ - 1) * (L₂ - 1)) := by
  sorry

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
