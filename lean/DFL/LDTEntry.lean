import DFL.Resonance

/-!
# Large deviations and resonance elimination (paper: Proposition `prop:LDT`, Lemma `lem:entry`,
Lemma `lem:double`)

Probability is the uniform measure on the sign configurations `q ∈ Cfg N` (i.i.d. Bernoulli
charges); probabilities are written as cardinalities, `P(S) = #S / 2^N`.

* `LDT λ` is the literature input: a function `γ` (the Lyapunov exponent) bounded below on `I_λ` by
  some `γ_* > 0`, and uniform large deviations for `log ‖Φ_m(E) u‖ / m` around `γ(E)`, for all unit
  vectors `u`, `E ∈ I_λ` and `m ≥ 1`. For the Bernoulli–Anderson model this is
  Bucaj–Damanik–Fillman–Gerbuz–VandenBoom–Wang–Zhang, Prop. 3.6 (vectorwise uniform LDT), together with
  Furstenberg's theorem and the continuity of `γ` (their Thms. 2.1 and 2.6); see the paper,
  Propositions `prop:nocritical` and `prop:LDT`. Only the vector form of the LDT is used.
* `entry_ldt` (paper: Lemma `lem:entry`): matrix entries `⟨w, M e⟩` of products `M` of `m` transfer
  matrices over distinct sites are smaller than `e^{κ m}` with probability `≤ C₁ e^{-c₁ m}`.
* `single_box`, `two_box` (paper: Lemma `lem:double`): a fixed energy is resonant for a box with
  probability `≤ 2 C₁ e^{-c₁ n}`, and two disjoint boxes share a resonant energy with probability
  `≤ 4 (8n+6) C₁ e^{-c₁ n}`.
-/

open Matrix Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The product of the transfer matrices `T_E(λ q_z)` along the list of sites `zs` (first site
leftmost). -/
def tlist (lam E : ℝ) (q : Cfg N) (zs : List (Fin N)) : Matrix (Fin 2) (Fin 2) ℝ :=
  (zs.map (fun z => transferT E (lam * sgn (q z)))).prod

/-- **Uniform large deviations** (literature input; paper: Propositions `prop:nocritical`, `prop:LDT`,
Eq. `eq:LDTvec`). -/
def LDT (lam : ℝ) : Prop :=
  ∃ γ : ℝ → ℝ, ∃ γs : ℝ, 0 < γs ∧ (∀ E ∈ Ilam lam, γs ≤ γ E) ∧
    ∀ ε : ℝ, 0 < ε → ∃ C c : ℝ, 0 < c ∧ ∀ m : ℕ, 1 ≤ m → ∀ E ∈ Ilam lam, ∀ u : Fin 2 → ℝ,
      u 0 ^ 2 + u 1 ^ 2 = 1 →
        ((univ.filter (fun q : Cfg m =>
          ε ≤ |Real.log (enorm2 (tlist lam E q (List.finRange m) *ᵥ u)) / m - γ E|)).card : ℝ) ≤
          C * Real.exp (-c * m) * 2 ^ m

/-- The conclusion of Lemma `lem:entry` for given constants `κ, C₁, c₁, m₀`. -/
def EntryBound (lam κ C1 c1 : ℝ) (m0 : ℕ) : Prop :=
  ∀ (N : ℕ) (zs : List (Fin N)), zs.Nodup → m0 ≤ zs.length → ∀ E ∈ Ilam lam,
    ∀ e w : Fin 2 → ℝ, e 0 ^ 2 + e 1 ^ 2 = 1 → w 0 ^ 2 + w 1 ^ 2 = 1 →
      ((univ.filter (fun q : Cfg N =>
        |w ⬝ᵥ (tlist lam E q zs *ᵥ e)| ≤ Real.exp (κ * zs.length))).card : ℝ) ≤
        C1 * Real.exp (-c1 * zs.length) * 2 ^ N

/-! ### Auxiliary lemmas -/

namespace LDTEntryAux

/-! #### Counting configurations -/

/-- There are `2^N` configurations. -/
theorem card_cfg : Fintype.card (Cfg N) = 2 ^ N := by
  simp [Cfg, ZMod.card]

/-- Mixing two configurations: `q` on `S`, `q'` off `S`. -/
def mix (S : Finset (Fin N)) (q q' : Cfg N) : Cfg N := fun i => if i ∈ S then q i else q' i

/-- **Independence lemma**: if `g` depends only on the coordinates in `S` and each event `R a`
only on the coordinates outside `S`, with `#{R a} ≤ B` for all `a`, then `#{q | R (g q) q} ≤ B`. -/
theorem card_filter_indep {α : Type*} (S : Finset (Fin N)) (g : Cfg N → α)
    (R : α → Cfg N → Prop) [∀ a, DecidablePred (R a)] [DecidablePred (fun q => R (g q) q)]
    (B : ℝ) (hg : ∀ q q' : Cfg N, (∀ i ∈ S, q i = q' i) → g q = g q')
    (hR : ∀ a (q q' : Cfg N), (∀ i, i ∉ S → q i = q' i) → (R a q ↔ R a q'))
    (hB : ∀ a, ((univ.filter (R a)).card : ℝ) ≤ B) :
    ((univ.filter (fun q => R (g q) q)).card : ℝ) ≤ B := by
  classical
  set T := (univ : Finset (Cfg N × Cfg N)).filter (fun p => R (g p.1) p.1) with hT
  set T' := (univ : Finset (Cfg N × Cfg N)).filter (fun p => R (g p.1) p.2) with hT'
  have h1 : (T.card : ℝ) = ((univ.filter (fun q => R (g q) q)).card : ℝ) * 2 ^ N := by
    rw [hT, Finset.card_filter, Fintype.sum_prod_type, Finset.card_filter]
    push_cast
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl (fun q _ => ?_)
    split_ifs <;> simp
  have h2 : (T'.card : ℝ) = ∑ q : Cfg N, ((univ.filter (R (g q))).card : ℝ) := by
    rw [hT', Finset.card_filter, Fintype.sum_prod_type]
    push_cast
    refine Finset.sum_congr rfl (fun q _ => ?_)
    rw [Finset.card_filter]
    push_cast
    rfl
  have h3 : T.card ≤ T'.card := by
    refine Finset.card_le_card_of_injOn (fun p => (mix S p.1 p.2, mix S p.2 p.1)) ?_ ?_
    · intro p hp
      simp only [hT, hT', Finset.coe_filter, Finset.mem_univ, true_and] at hp ⊢
      have e1 : g (mix S p.1 p.2) = g p.1 := hg _ _ (fun i hi => by simp [mix, hi])
      show R (g (mix S p.1 p.2)) (mix S p.2 p.1)
      rw [e1]
      exact (hR _ _ _ (fun i hi => by simp [mix, hi])).2 hp
    · intro p _ p' _ h
      simp only [Prod.mk.injEq] at h
      obtain ⟨h1, h2⟩ := h
      ext i
      · by_cases hi : i ∈ S
        · have := congrFun h1 i; simpa [mix, hi] using this
        · have := congrFun h2 i; simpa [mix, hi] using this
      · by_cases hi : i ∈ S
        · have := congrFun h2 i; simpa [mix, hi] using this
        · have := congrFun h1 i; simpa [mix, hi] using this
  have h4 : (T'.card : ℝ) ≤ 2 ^ N * B := by
    rw [h2]
    calc ∑ q : Cfg N, ((univ.filter (R (g q))).card : ℝ) ≤ ∑ _q : Cfg N, B :=
          Finset.sum_le_sum (fun q _ => hB (g q))
      _ = 2 ^ N * B := by simp [Finset.card_univ]
  have h5 : ((univ.filter (fun q => R (g q) q)).card : ℝ) * 2 ^ N ≤ B * 2 ^ N := by
    rw [← h1]
    calc (T.card : ℝ) ≤ T'.card := by exact_mod_cast h3
      _ ≤ 2 ^ N * B := h4
      _ = B * 2 ^ N := by ring
  exact le_of_mul_le_mul_right h5 (by positivity)

/-- **Independence lemma, union-bound form**: if the finite set `Z q` depends only on the coordinates
in `S` and has at most `K` elements, and each event `R a` depends only on the coordinates outside `S`
with `#{R a} ≤ B`, then `#{q | ∃ z ∈ Z q, R z q} ≤ K B`. -/
theorem card_filter_indep_union {α : Type*} [DecidableEq α] (S : Finset (Fin N))
    (Z : Cfg N → Finset α) (R : α → Cfg N → Prop) [∀ a, DecidablePred (R a)] (K : ℕ) (B : ℝ)
    (hB0 : 0 ≤ B) (hZ : ∀ q q' : Cfg N, (∀ i ∈ S, q i = q' i) → Z q = Z q')
    (hK : ∀ q, (Z q).card ≤ K)
    (hR : ∀ a (q q' : Cfg N), (∀ i, i ∉ S → q i = q' i) → (R a q ↔ R a q'))
    (hB : ∀ a, ((univ.filter (R a)).card : ℝ) ≤ B)
    (P : Cfg N → Prop) [DecidablePred P] (hP : ∀ q, P q → ∃ z ∈ Z q, R z q) :
    ((univ.filter P).card : ℝ) ≤ K * B := by
  classical
  let R' : Finset α → Cfg N → Prop := fun Zs q => Zs.card ≤ K ∧ ∃ z ∈ Zs, R z q
  have h1 : univ.filter P ⊆ univ.filter (fun q => R' (Z q) q) := by
    intro q hq
    simp only [mem_filter, mem_univ, true_and] at hq ⊢
    exact ⟨hK q, hP q hq⟩
  refine (Nat.cast_le.2 (card_le_card h1)).trans ?_
  refine card_filter_indep S Z R' (K * B) hZ ?_ ?_
  · intro Zs q q' hqq'
    simp only [R']
    constructor
    · rintro ⟨hc, z, hz, hRz⟩; exact ⟨hc, z, hz, (hR z q q' hqq').1 hRz⟩
    · rintro ⟨hc, z, hz, hRz⟩; exact ⟨hc, z, hz, (hR z q q' hqq').2 hRz⟩
  · intro Zs
    by_cases hc : Zs.card ≤ K
    · have hsub : univ.filter (R' Zs) ⊆ Zs.biUnion (fun z => univ.filter (R z)) := by
        intro q hq
        simp only [mem_filter, mem_univ, true_and, R'] at hq
        obtain ⟨_, z, hz, hRz⟩ := hq
        exact mem_biUnion.2 ⟨z, hz, by simpa using hRz⟩
      calc ((univ.filter (R' Zs)).card : ℝ)
          ≤ (Zs.biUnion (fun z => univ.filter (R z))).card := by
            exact_mod_cast card_le_card hsub
        _ ≤ ∑ z ∈ Zs, ((univ.filter (R z)).card : ℝ) := by exact_mod_cast card_biUnion_le
        _ ≤ ∑ _z ∈ Zs, B := sum_le_sum (fun z _ => hB z)
        _ = Zs.card * B := by simp
        _ ≤ K * B := mul_le_mul_of_nonneg_right (by exact_mod_cast hc) hB0
    · have : univ.filter (R' Zs) = ∅ := by
        ext q; simp [R', hc]
      rw [this, card_empty, Nat.cast_zero]
      positivity

/-! #### Products over distinct sites -/

/-- Restriction of a configuration to the sites of a list. -/
def proj (zs : List (Fin N)) : Cfg N →+ Cfg zs.length where
  toFun q k := q (zs.get k)
  map_zero' := rfl
  map_add' _ _ := rfl

/-- The restriction to the sites of a list without repetitions is surjective. -/
theorem proj_surjective (zs : List (Fin N)) (hnd : zs.Nodup) :
    Function.Surjective (proj zs) := by
  have hinj : Function.Injective zs.get := List.nodup_iff_injective_get.1 hnd
  intro p
  obtain ⟨q, hq⟩ := hinj.surjective_comp_right p
  exact ⟨q, hq⟩

/-- The restriction to distinct sites maps the uniform law to the uniform law. -/
theorem card_filter_proj (zs : List (Fin N)) (hnd : zs.Nodup) (P : Cfg zs.length → Prop)
    [DecidablePred P] :
    (univ.filter (fun q : Cfg N => P (proj zs q))).card * 2 ^ zs.length =
      (univ.filter P).card * 2 ^ N := by
  classical
  set F := (univ.filter (fun q : Cfg N => proj zs q = 0)).card with hF
  have hfib : ∀ p : Cfg zs.length,
      (univ.filter (fun q : Cfg N => proj zs q = p)).card = F := fun p =>
    AddMonoidHom.card_fiber_eq_of_mem_range (proj zs) (proj_surjective zs hnd p)
      (proj_surjective zs hnd 0)
  have htot : 2 ^ N = 2 ^ zs.length * F := by
    have := Finset.card_eq_sum_card_fiberwise (f := proj zs) (s := univ) (t := univ)
      (fun _ _ => Finset.mem_coe.2 (Finset.mem_univ _))
    rw [Finset.card_univ, card_cfg] at this
    rw [this]
    simp [hfib]
  have hP : (univ.filter (fun q : Cfg N => P (proj zs q))).card = (univ.filter P).card * F := by
    rw [Finset.card_eq_sum_card_fiberwise (f := proj zs) (t := univ.filter P)
      (fun q hq => by simpa using hq)]
    rw [Finset.card_eq_sum_ones (univ.filter P), Finset.sum_mul, one_mul]
    refine Finset.sum_congr rfl (fun p hp => ?_)
    rw [← hfib p, Finset.filter_filter]
    congr 1
    refine Finset.filter_congr (fun q _ => ?_)
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hp
    constructor
    · exact fun h => h.2
    · intro h; exact ⟨h ▸ hp, h⟩
  rw [hP, htot]; ring

/-- `tlist` along `zs` is `tlist` along `0, …, m-1` of the restricted configuration. -/
theorem tlist_proj (lam E : ℝ) (q : Cfg N) (zs : List (Fin N)) :
    tlist lam E (proj zs q) (List.finRange zs.length) = tlist lam E q zs := by
  unfold tlist
  congr 1
  conv_rhs => rw [← List.map_get_finRange zs]
  rw [List.map_map]
  rfl

/-- `tlist` depends only on the configuration at the sites of the list. -/
theorem tlist_congr (lam E : ℝ) {q q' : Cfg N} (zs : List (Fin N))
    (h : ∀ z ∈ zs, q z = q' z) : tlist lam E q zs = tlist lam E q' zs := by
  unfold tlist
  congr 1
  exact List.map_congr_left (fun z hz => by rw [h z hz])

/-- The large deviation estimate, transferred to products over distinct sites of `Fin N`. -/
theorem ldt_transfer (lam : ℝ) (γ : ℝ → ℝ) (ε C c : ℝ)
    (hLDT : ∀ m : ℕ, 1 ≤ m → ∀ E ∈ Ilam lam, ∀ u : Fin 2 → ℝ, u 0 ^ 2 + u 1 ^ 2 = 1 →
      ((univ.filter (fun q : Cfg m =>
        ε ≤ |Real.log (enorm2 (tlist lam E q (List.finRange m) *ᵥ u)) / m - γ E|)).card : ℝ) ≤
        C * Real.exp (-c * m) * 2 ^ m)
    (zs : List (Fin N)) (hnd : zs.Nodup) (hm : 1 ≤ zs.length) (E : ℝ) (hE : E ∈ Ilam lam)
    (u : Fin 2 → ℝ) (hu : u 0 ^ 2 + u 1 ^ 2 = 1) :
    ((univ.filter (fun q : Cfg N =>
        ε ≤ |Real.log (enorm2 (tlist lam E q zs *ᵥ u)) / zs.length - γ E|)).card : ℝ) ≤
        C * Real.exp (-c * zs.length) * 2 ^ N := by
  have key := card_filter_proj zs hnd (fun p =>
    ε ≤ |Real.log (enorm2 (tlist lam E p (List.finRange zs.length) *ᵥ u)) / zs.length - γ E|)
  simp only [tlist_proj] at key
  have hb := hLDT zs.length hm E hE u hu
  have key' : ((univ.filter (fun q : Cfg N =>
        ε ≤ |Real.log (enorm2 (tlist lam E q zs *ᵥ u)) / zs.length - γ E|)).card : ℝ) *
        2 ^ zs.length = ((univ.filter (fun q : Cfg zs.length =>
        ε ≤ |Real.log (enorm2 (tlist lam E q (List.finRange zs.length) *ᵥ u)) / zs.length -
          γ E|)).card : ℝ) * 2 ^ N := by
    exact_mod_cast key
  have h2 : (0 : ℝ) < 2 ^ zs.length := by positivity
  refine le_of_mul_le_mul_right ?_ h2
  rw [key']
  calc _ ≤ C * Real.exp (-c * zs.length) * 2 ^ zs.length * 2 ^ N :=
        mul_le_mul_of_nonneg_right hb (by positivity)
    _ = _ := by ring

/-! #### Deterministic facts about transfer matrices -/

/-- `D = diag(1, -1)`. -/
def Dm : Matrix (Fin 2) (Fin 2) ℝ := !![1, 0; 0, -1]

/-- `D² = 1`. -/
theorem Dm_mul_Dm : Dm * Dm = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [Dm]

/-- `T_E(v)ᵀ = D T_E(v) D`. -/
theorem transferT_transpose' (E v : ℝ) : (transferT E v)ᵀ = Dm * transferT E v * Dm := by
  ext i j; fin_cases i <;> fin_cases j <;> simp [Dm, transferT, Matrix.mul_apply, Fin.sum_univ_two]

/-- Transfer matrices have determinant one. -/
theorem det_transferT' (E v : ℝ) : (transferT E v).det = 1 := by
  simp [transferT, Matrix.det_fin_two_of]

/-- `tlist` has determinant one. -/
theorem det_tlist (lam E : ℝ) (q : Cfg N) (zs : List (Fin N)) : (tlist lam E q zs).det = 1 := by
  induction zs with
  | nil => simp [tlist]
  | cons z zs ih =>
    simp only [tlist, List.map_cons, List.prod_cons, Matrix.det_mul, det_transferT'] at ih ⊢
    rw [ih, one_mul]

/-- `Mᵀ = D M̌ D` for `M = tlist` along `zs` and `M̌ = tlist` along `zs.reverse`. -/
theorem tlist_transpose (lam E : ℝ) (q : Cfg N) (zs : List (Fin N)) :
    (tlist lam E q zs)ᵀ = Dm * tlist lam E q zs.reverse * Dm := by
  induction zs with
  | nil => simp [tlist, Dm_mul_Dm]
  | cons z zs ih =>
    simp only [tlist, List.map_cons, List.prod_cons, List.reverse_cons, List.map_append,
      List.prod_append, List.map_nil, List.prod_nil, mul_one] at ih ⊢
    rw [Matrix.transpose_mul, ih, transferT_transpose']
    simp only [Matrix.mul_assoc]
    rw [← Matrix.mul_assoc Dm Dm, Dm_mul_Dm, Matrix.one_mul]

/-- `v ⬝ᵥ v = v₀² + v₁²`. -/
theorem dot_self_two (v : Fin 2 → ℝ) : v ⬝ᵥ v = v 0 ^ 2 + v 1 ^ 2 := by
  simp [dotProduct, Fin.sum_univ_two, sq]

/-- `‖v‖² = v₀² + v₁²`. -/
theorem enorm2_sq (v : Fin 2 → ℝ) : enorm2 v ^ 2 = v 0 ^ 2 + v 1 ^ 2 := by
  unfold enorm2; rw [Real.sq_sqrt (by positivity)]

/-- `D v = (v₀, -v₁)`. -/
theorem Dm_mulVec (v : Fin 2 → ℝ) : Dm *ᵥ v = ![v 0, -v 1] := by
  ext i; fin_cases i <;> simp [Dm, mulVec, dotProduct, Fin.sum_univ_two]

/-- Deterministic core of `entry_ldt`: if `‖Me‖, ‖Mᵀw‖ > X`, `‖M‖_F ≤ 2Y` and
`2KY + 1 < X² - K²`, then `|⟨w, M e⟩| > K` (from `sl2_entry_ineq`). -/
theorem entry_det (M : Matrix (Fin 2) (Fin 2) ℝ) (hM : M.det = 1) (e w : Fin 2 → ℝ)
    (he : e 0 ^ 2 + e 1 ^ 2 = 1) (hw : w 0 ^ 2 + w 1 ^ 2 = 1) (X K Y : ℝ) (hK : 0 ≤ K)
    (h1 : X ^ 2 < (M *ᵥ e) ⬝ᵥ (M *ᵥ e)) (h2 : X ^ 2 < (Mᵀ *ᵥ w) ⬝ᵥ (Mᵀ *ᵥ w))
    (h3 : ∑ i, ∑ j, M i j ^ 2 ≤ 4 * Y ^ 2) (h4 : 2 * K * Y + 1 < X ^ 2 - K ^ 2) (hY : 0 ≤ Y) :
    K < |w ⬝ᵥ (M *ᵥ e)| := by
  by_contra hcon
  push Not at hcon
  have hsl := sl2_entry_ineq M hM e w he hw
  set a := w ⬝ᵥ (M *ᵥ e) with ha
  have ha2 : a ^ 2 ≤ K ^ 2 := by
    rw [← sq_abs]; exact pow_le_pow_left₀ (abs_nonneg a) hcon 2
  have hsq : Real.sqrt (∑ i, ∑ j, M i j ^ 2) ≤ 2 * Y := by
    rw [Real.sqrt_le_left (by linarith)]; nlinarith
  have hP : 0 < X ^ 2 - K ^ 2 := by nlinarith
  have hA : X ^ 2 - K ^ 2 < (Mᵀ *ᵥ w) ⬝ᵥ (Mᵀ *ᵥ w) - a ^ 2 := by linarith
  have hB : X ^ 2 - K ^ 2 < (M *ᵥ e) ⬝ᵥ (M *ᵥ e) - a ^ 2 := by linarith
  have hlhs : (X ^ 2 - K ^ 2) ^ 2 <
      ((Mᵀ *ᵥ w) ⬝ᵥ (Mᵀ *ᵥ w) - a ^ 2) * ((M *ᵥ e) ⬝ᵥ (M *ᵥ e) - a ^ 2) := by
    rw [sq]; exact mul_lt_mul'' hA hB hP.le hP.le
  have hrhs : |a| * Real.sqrt (∑ i, ∑ j, M i j ^ 2) + 1 ≤ 2 * K * Y + 1 := by
    have := mul_le_mul hcon hsq (Real.sqrt_nonneg _) hK
    linarith
  have hrhs2 : (|a| * Real.sqrt (∑ i, ∑ j, M i j ^ 2) + 1) ^ 2 < (X ^ 2 - K ^ 2) ^ 2 := by
    have h0 : 0 ≤ |a| * Real.sqrt (∑ i, ∑ j, M i j ^ 2) + 1 := by positivity
    exact pow_lt_pow_left₀ (by linarith) h0 two_ne_zero
  linarith

/-- The elementary inequality behind the choice of `m₀` (`h = e^{γ_* m/8}`, `X = e^{(γ-γ_*/8)m}`). -/
theorem alg_ineq (X h : ℝ) (hh : 4 ≤ h) (hX : h ^ 7 ≤ X) :
    2 * h ^ 4 * (X * h ^ 2) + 1 < X ^ 2 - (h ^ 4) ^ 2 := by
  have h1 : 1 ≤ h := by linarith
  have hp6 : 0 < h ^ 6 := by positivity
  have h6 : 4 * h ^ 6 ≤ X := by
    calc 4 * h ^ 6 ≤ h * h ^ 6 := by nlinarith
      _ = h ^ 7 := by ring
      _ ≤ X := hX
  have hX0 : 0 < X := by nlinarith
  have hA : 4 * h ^ 6 * X ≤ X ^ 2 := by nlinarith
  have hh6 : h ≤ h ^ 6 := by
    calc h = h * 1 := by ring
      _ ≤ h * h ^ 5 := by
          apply mul_le_mul_of_nonneg_left _ (by linarith)
          exact one_le_pow₀ h1
      _ = h ^ 6 := by ring
  have hB : h ^ 8 ≤ h ^ 6 * X := by
    calc h ^ 8 = h * h ^ 7 := by ring
      _ ≤ h * X := by apply mul_le_mul_of_nonneg_left hX (by linarith)
      _ ≤ h ^ 6 * X := by apply mul_le_mul_of_nonneg_right hh6 hX0.le
  have hC : 1 < h ^ 6 * X := by
    have : 1 ≤ h ^ 6 := one_le_pow₀ h1
    nlinarith
  nlinarith

/-- Lower bound on `t` outside the large deviation event. -/
theorem lower_of_ldt {t m γ ε : ℝ} (hm : 0 < m) (hγε : ε < γ) (h : |Real.log t / m - γ| < ε)
    (ht : 0 ≤ t) : Real.exp ((γ - ε) * m) < t := by
  have h' := (abs_lt.1 h).1
  have hl : (γ - ε) * m < Real.log t := by
    have : γ - ε < Real.log t / m := by linarith
    rwa [lt_div_iff₀ hm] at this
  rcases ht.lt_or_eq with htp | ht0
  · rwa [Real.lt_log_iff_exp_lt htp] at hl
  · subst ht0
    rw [Real.log_zero] at hl
    have : 0 < (γ - ε) * m := mul_pos (by linarith) hm
    linarith

/-- Upper bound on `t` outside the large deviation event. -/
theorem upper_of_ldt {t m γ ε : ℝ} (hm : 0 < m) (h : |Real.log t / m - γ| < ε)
    (ht : 0 ≤ t) : t < Real.exp ((γ + ε) * m) := by
  have h' := (abs_lt.1 h).2
  have hl : Real.log t < (γ + ε) * m := by
    have : Real.log t / m < γ + ε := by linarith
    rwa [div_lt_iff₀ hm] at this
  rcases ht.lt_or_eq with htp | ht0
  · rwa [Real.log_lt_iff_lt_exp htp] at hl
  · subst ht0; exact Real.exp_pos _

/-! #### Boxes of the chain -/

/-- `potZ λ q` at `x` depends only on `q` at the site `x` (if any). -/
theorem potZ_congr (lam : ℝ) {q q' : Cfg N} {x : ℤ}
    (h : ∀ z : Fin N, (z : ℤ) = x → q z = q' z) : potZ lam q x = potZ lam q' x := by
  unfold potZ extZ
  split_ifs with hx
  · have hq := h ⟨x.toNat, by omega⟩ (by simp; omega)
    simp only [hq]
  · rfl

/-- The sites `s + m - 1, …, s` lying in `[0, N-1]`, later sites first. -/
def zsL (N : ℕ) (s : ℤ) : ℕ → List (Fin N)
  | 0 => []
  | (m + 1) =>
    (if h : 0 ≤ s + m ∧ s + m < N then [⟨(s + m).toNat, by omega⟩] else []) ++ zsL N s m

/-- The sites of `zsL N s m` lie in `[s, s+m-1]`. -/
theorem mem_zsL {s : ℤ} {m : ℕ} {z : Fin N} (hz : z ∈ zsL N s m) :
    s ≤ (z : ℤ) ∧ (z : ℤ) < s + m := by
  induction m with
  | zero => simp [zsL] at hz
  | succ m ih =>
    simp only [zsL, List.mem_append] at hz
    rcases hz with hz | hz
    · split_ifs at hz with h
      · simp only [List.mem_singleton] at hz
        subst hz
        simp only
        push_cast
        omega
      · simp at hz
    · have := ih hz
      push_cast
      omega

/-- `zsL N s m` has no repetitions. -/
theorem nodup_zsL (s : ℤ) (m : ℕ) : (zsL N s m).Nodup := by
  induction m with
  | zero => simp [zsL]
  | succ m ih =>
    simp only [zsL]
    rw [List.nodup_append]
    refine ⟨?_, ih, ?_⟩
    · split_ifs <;> simp
    · intro a ha b hb hab
      split_ifs at ha with h
      · simp only [List.mem_singleton] at ha
        subst ha
        subst hab
        have := mem_zsL hb
        simp only at this
        omega
      · simp at ha

/-- `zsL N s m` has length `m` if `[s, s+m-1] ⊆ [0, N-1]`. -/
theorem length_zsL {s : ℤ} {m : ℕ} (hs : 0 ≤ s) (hm : s + m ≤ N) : (zsL N s m).length = m := by
  induction m with
  | zero => simp [zsL]
  | succ m ih =>
    simp only [zsL, List.length_append]
    push_cast at hm
    rw [dite_eq_left ⟨by omega, by omega⟩, ih (by omega)]
    simp only [List.length_singleton]
    omega

/-- `tprodL` as a `tlist`. -/
theorem tprodL_eq_tlist (lam E : ℝ) (q : Cfg N) {s : ℤ} {m : ℕ} (hs : 0 ≤ s) (hm : s + m ≤ N) :
    tprodL (potZ lam q) E s m = tlist lam E q (zsL N s m) := by
  induction m with
  | zero => simp [tprodL, zsL, tlist]
  | succ m ih =>
    push_cast at hm
    have hm' : s + m < N := by omega
    rw [tprodL, ih (by omega), zsL, dite_eq_left ⟨by omega, hm'⟩]
    simp only [tlist, List.singleton_append, List.map_cons, List.prod_cons]
    congr 2
    unfold potZ extZ
    rw [dite_eq_left ⟨by omega, hm'⟩]

/-- `tprodR` as a `tlist`. -/
theorem tprodR_eq_tlist (lam E : ℝ) (q : Cfg N) {s : ℤ} {m : ℕ} (hs : 0 ≤ s) (hm : s + m ≤ N) :
    tprodR (potZ lam q) E s m = tlist lam E q (zsL N s m).reverse := by
  induction m with
  | zero => simp [tprodR, zsL, tlist]
  | succ m ih =>
    push_cast at hm
    have hm' : s + m < N := by omega
    rw [tprodR, ih (by omega), zsL, dite_eq_left ⟨by omega, hm'⟩]
    simp only [tlist, List.singleton_append, List.reverse_cons, List.map_append,
      List.prod_append, List.map_cons, List.map_nil, List.prod_cons, List.prod_nil, mul_one]
    congr 2
    unfold potZ extZ
    rw [dite_eq_left ⟨by omega, hm'⟩]

/-- `EntryBound` forces `C₁ ≥ 0`. -/
theorem C1_nonneg {lam κ C1 c1 : ℝ} {m0 : ℕ} (hEB : EntryBound lam κ C1 c1 m0) : 0 ≤ C1 := by
  have h0 : (0 : ℝ) ∈ Ilam lam := by
    unfold Ilam
    constructor <;> linarith [abs_nonneg lam]
  have h := hEB m0 (List.finRange m0) (List.nodup_finRange m0) (by simp) 0 h0 ![1, 0] ![1, 0]
    (by simp) (by simp)
  have h2 : 0 ≤ C1 * Real.exp (-c1 * (List.finRange m0).length) * 2 ^ m0 :=
    le_trans (Nat.cast_nonneg _) h
  have hpos : 0 < Real.exp (-c1 * (List.finRange m0).length) * 2 ^ m0 := by positivity
  rw [mul_assoc] at h2
  by_contra hneg
  push Not at hneg
  have := mul_neg_of_neg_of_pos hneg hpos
  linarith

/-- **One end of a box**: if `M = tlist` over `n ≥ m₀` distinct sites and `Y ≠ 0` depends only on
the coordinates in a set `S` disjoint from these sites, then `e^{-κn} |(M Y)₀| ≤ |Y₁|` has
probability at most `C₁ e^{-c₁ n}`. -/
theorem end_bound (lam κ C1 c1 : ℝ) (m0 : ℕ) (hEB : EntryBound lam κ C1 c1 m0) (E : ℝ)
    (hE : E ∈ Ilam lam) (n : ℕ) (zs : List (Fin N)) (hnd : zs.Nodup) (hlen : zs.length = n)
    (hn : m0 ≤ n) (S : Finset (Fin N)) (hS : ∀ z ∈ zs, z ∉ S) (Yf : Cfg N → (Fin 2 → ℝ))
    (hY : ∀ q q' : Cfg N, (∀ i ∈ S, q i = q' i) → Yf q = Yf q')
    (hY0 : ∀ q, Yf q 0 ≠ 0 ∨ Yf q 1 ≠ 0) (P : Cfg N → Prop) [DecidablePred P]
    (hP : ∀ q, P q → Real.exp (-(κ * n)) * |(tlist lam E q zs *ᵥ Yf q) 0| ≤ |Yf q 1|) :
    ((univ.filter P).card : ℝ) ≤ C1 * Real.exp (-c1 * n) * 2 ^ N := by
  classical
  have hC1 := C1_nonneg hEB
  have hB0 : 0 ≤ C1 * Real.exp (-c1 * n) * 2 ^ N := by positivity
  let g : Cfg N → (Fin 2 → ℝ) := fun q => (enorm2 (Yf q))⁻¹ • Yf q
  let R : (Fin 2 → ℝ) → Cfg N → Prop := fun y q =>
    y 0 ^ 2 + y 1 ^ 2 = 1 ∧ |![1, 0] ⬝ᵥ (tlist lam E q zs *ᵥ y)| ≤ Real.exp (κ * n)
  have hsub : univ.filter P ⊆ univ.filter (fun q => R (g q) q) := by
    intro q hq
    simp only [mem_filter, mem_univ, true_and, g, R] at hq ⊢
    have hpq := hP q hq
    set Y := Yf q with hYdef
    have hr2 := enorm2_sq Y
    have hr : 0 < enorm2 Y := by
      unfold enorm2
      apply Real.sqrt_pos.2
      rcases hY0 q with h | h
      · have := sq_pos_of_ne_zero h; positivity
      · have := sq_pos_of_ne_zero h; positivity
    refine ⟨?_, ?_⟩
    · have hne : enorm2 Y ≠ 0 := hr.ne'
      simp only [Pi.smul_apply, smul_eq_mul]
      rw [mul_pow, mul_pow, ← mul_add, ← hr2, inv_pow, inv_mul_cancel₀ (pow_ne_zero 2 hne)]
    · have hθ : Real.exp (-(κ * n)) * Real.exp (κ * n) = 1 := by
        rw [← Real.exp_add]; simp
      have hY1 : |Y 1| ≤ enorm2 Y := Real.abs_le_sqrt (by nlinarith [sq_nonneg (Y 0)])
      have h1 : ![1, 0] ⬝ᵥ (tlist lam E q zs *ᵥ ((enorm2 Y)⁻¹ • Y)) =
          (enorm2 Y)⁻¹ * (tlist lam E q zs *ᵥ Y) 0 := by
        rw [Matrix.mulVec_smul]; simp [dotProduct, Fin.sum_univ_two]
      rw [h1, abs_mul, abs_inv, abs_of_pos hr, inv_mul_le_iff₀ hr]
      calc |(tlist lam E q zs *ᵥ Y) 0|
          = Real.exp (κ * n) * (Real.exp (-(κ * n)) * |(tlist lam E q zs *ᵥ Y) 0|) := by
            rw [← mul_assoc, mul_comm (Real.exp (κ * n)), hθ, one_mul]
        _ ≤ Real.exp (κ * n) * enorm2 Y := by gcongr; linarith
        _ = enorm2 Y * Real.exp (κ * n) := by ring
  refine (Nat.cast_le.2 (card_le_card hsub)).trans ?_
  refine card_filter_indep S g R _ ?_ ?_ ?_
  · intro q q' h
    simp only [g]
    rw [hY q q' h]
  · intro y q q' h
    have : tlist lam E q zs = tlist lam E q' zs :=
      tlist_congr lam E zs (fun z hz => h z (hS z hz))
    simp only [R, this]
  · intro y
    by_cases hy : y 0 ^ 2 + y 1 ^ 2 = 1
    · have hEBy := hEB N zs hnd (hlen ▸ hn) E hE y ![1, 0] hy (by simp)
      rw [hlen] at hEBy
      refine le_trans ?_ hEBy
      exact_mod_cast card_le_card (fun q hq => by
        simp only [mem_filter, mem_univ, true_and, R] at hq ⊢
        exact hq.2)
    · have : univ.filter (R y) = ∅ := by
        ext q; simp [R, hy]
      rw [this, card_empty, Nat.cast_zero]
      exact hB0

end LDTEntryAux

open LDTEntryAux

/-- **Large deviations for entries in arbitrary directions** (paper: Lemma `lem:entry`). -/
theorem entry_ldt (lam : ℝ) (hL : LDT lam) :
    ∃ κ C1 c1 : ℝ, ∃ m0 : ℕ, 0 < κ ∧ 0 < c1 ∧ 0 ≤ C1 ∧ EntryBound lam κ C1 c1 m0 := by
  obtain ⟨γ, γs, hγs, hγ, hLDT⟩ := hL
  obtain ⟨C, c, hc, hC⟩ := hLDT (γs / 8) (by positivity)
  refine ⟨γs / 2, 4 * max C 0, c, ⌈24 / γs⌉₊ + 1, by positivity, hc, by positivity, ?_⟩
  intro N zs hnd hlen E hE e w he hw
  have hm1 : 1 ≤ zs.length := by omega
  have hmpos : (0 : ℝ) < zs.length := by exact_mod_cast hm1
  have hmr : 24 / γs < (zs.length : ℝ) := by
    have h1 := Nat.le_ceil (24 / γs)
    have h2 : ((⌈24 / γs⌉₊ + 1 : ℕ) : ℝ) ≤ zs.length := by exact_mod_cast hlen
    push_cast at h2
    linarith
  have hεm : 3 ≤ γs / 8 * zs.length := by
    rw [div_lt_iff₀ hγs] at hmr; nlinarith
  have hγE : γs ≤ γ E := hγ E hE
  set h := Real.exp (γs / 8 * zs.length) with hh
  have hh4 : 4 ≤ h := by have := Real.add_one_le_exp (γs / 8 * zs.length); linarith
  set X := Real.exp ((γ E - γs / 8) * zs.length) with hX
  set Y := Real.exp ((γ E + γs / 8) * zs.length) with hY
  set K := Real.exp (γs / 2 * zs.length) with hK
  have hYX : Y = X * h ^ 2 := by
    rw [hY, hX, hh, sq, ← Real.exp_add, ← Real.exp_add]; congr 1; ring
  have hKh : K = h ^ 4 := by
    rw [hK, hh, ← Real.exp_nat_mul]; congr 1; push_cast; ring
  have hXh : h ^ 7 ≤ X := by
    rw [hX, hh, ← Real.exp_nat_mul, Real.exp_le_exp]; push_cast; nlinarith
  have halg : 2 * K * Y + 1 < X ^ 2 - K ^ 2 := by
    rw [hYX, hKh]; exact alg_ineq X h hh4 hXh
  have b1 := ldt_transfer lam γ (γs / 8) C c hC zs hnd hm1 E hE e he
  have b2 := ldt_transfer lam γ (γs / 8) C c hC zs.reverse (List.nodup_reverse.2 hnd)
    (by simpa using hm1) E hE (Dm *ᵥ w) (by rw [Dm_mulVec]; simpa using hw)
  have b3 := ldt_transfer lam γ (γs / 8) C c hC zs hnd hm1 E hE ![1, 0] (by simp)
  have b4 := ldt_transfer lam γ (γs / 8) C c hC zs hnd hm1 E hE ![0, 1] (by simp)
  simp only [List.length_reverse] at b2
  set A1 := univ.filter (fun q : Cfg N =>
    γs / 8 ≤ |Real.log (enorm2 (tlist lam E q zs *ᵥ e)) / zs.length - γ E|) with hA1
  set A2 := univ.filter (fun q : Cfg N =>
    γs / 8 ≤ |Real.log (enorm2 (tlist lam E q zs.reverse *ᵥ (Dm *ᵥ w))) / zs.length - γ E|)
    with hA2
  set A3 := univ.filter (fun q : Cfg N =>
    γs / 8 ≤ |Real.log (enorm2 (tlist lam E q zs *ᵥ ![1, 0])) / zs.length - γ E|) with hA3
  set A4 := univ.filter (fun q : Cfg N =>
    γs / 8 ≤ |Real.log (enorm2 (tlist lam E q zs *ᵥ ![0, 1])) / zs.length - γ E|) with hA4
  have hsub : univ.filter (fun q : Cfg N => |w ⬝ᵥ (tlist lam E q zs *ᵥ e)| ≤ K) ⊆
      (A1 ∪ A2) ∪ (A3 ∪ A4) := by
    intro q hq
    simp only [hA1, hA2, hA3, hA4, mem_filter, mem_univ, true_and, mem_union] at hq ⊢
    by_contra hcon
    simp only [not_or, not_le] at hcon
    obtain ⟨⟨g1, g2⟩, g3, g4⟩ := hcon
    set M := tlist lam E q zs with hMdef
    have hε : γs / 8 < γ E := by linarith
    have l1 := lower_of_ldt hmpos hε g1 (Real.sqrt_nonneg _)
    have l2 := lower_of_ldt hmpos hε g2 (Real.sqrt_nonneg _)
    have u3 := upper_of_ldt hmpos g3 (Real.sqrt_nonneg _)
    have u4 := upper_of_ldt hmpos g4 (Real.sqrt_nonneg _)
    have hX0 : 0 < X := Real.exp_pos _
    have hY0 : 0 < Y := Real.exp_pos _
    have d1 : X ^ 2 < (M *ᵥ e) ⬝ᵥ (M *ᵥ e) := by
      rw [dot_self_two, ← enorm2_sq]; exact pow_lt_pow_left₀ l1 hX0.le two_ne_zero
    have d2 : X ^ 2 < (Mᵀ *ᵥ w) ⬝ᵥ (Mᵀ *ᵥ w) := by
      have hT : Mᵀ *ᵥ w = Dm *ᵥ (tlist lam E q zs.reverse *ᵥ (Dm *ᵥ w)) := by
        rw [hMdef, tlist_transpose, ← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec]
      rw [hT, dot_self_two, Dm_mulVec]
      simp only [Matrix.cons_val_zero, Matrix.cons_val_one, neg_sq]
      rw [← enorm2_sq]; exact pow_lt_pow_left₀ l2 hX0.le two_ne_zero
    have d3 : ∑ i, ∑ j, M i j ^ 2 ≤ 4 * Y ^ 2 := by
      have s3 : enorm2 (M *ᵥ ![1, 0]) ^ 2 < Y ^ 2 :=
        pow_lt_pow_left₀ u3 (Real.sqrt_nonneg _) two_ne_zero
      have s4 : enorm2 (M *ᵥ ![0, 1]) ^ 2 < Y ^ 2 :=
        pow_lt_pow_left₀ u4 (Real.sqrt_nonneg _) two_ne_zero
      rw [enorm2_sq] at s3 s4
      simp [Matrix.mulVec, dotProduct, Fin.sum_univ_two] at s3 s4
      simp only [Fin.sum_univ_two]
      nlinarith
    have := entry_det M (det_tlist lam E q zs) e w he hw X K Y (Real.exp_pos _).le d1 d2 d3 halg
      hY0.le
    linarith
  have hcard := (card_le_card hsub).trans ((card_union_le _ _).trans
    (add_le_add (card_union_le _ _) (card_union_le _ _)))
  have hCm : C * Real.exp (-c * zs.length) * 2 ^ N ≤
      max C 0 * Real.exp (-c * zs.length) * 2 ^ N := by
    gcongr; exact le_max_left _ _
  calc ((univ.filter (fun q : Cfg N => |w ⬝ᵥ (tlist lam E q zs *ᵥ e)| ≤ K)).card : ℝ)
      ≤ (A1.card + A2.card + (A3.card + A4.card) : ℕ) := by exact_mod_cast hcard
    _ = (A1.card : ℝ) + A2.card + (A3.card + A4.card) := by push_cast; ring
    _ ≤ _ := by linarith

open Classical in
/-- **Single-box estimate** (paper: Lemma `lem:double`, Eq. `eq:single`): for a box
`[a, a+L-1] ⊆ [0, N-1]` and `x = a + i` whose selected ends are at distance `n` from `x`, a fixed
`E ∈ I_λ` is resonant with probability at most `2 C₁ e^{-c₁ n}` (threshold `e^{-κn}`). -/
theorem single_box (lam κ C1 c1 : ℝ) (m0 : ℕ) (hEB : EntryBound lam κ C1 c1 m0) (a : ℤ)
    (L i n : ℕ) (endL endR : Bool) (E : ℝ) (ha : 0 ≤ a) (haL : a + L ≤ N) (hi : i < L)
    (hn : m0 ≤ n) (hR : endR = true → L - 1 - i = n) (hLe : endL = true → i = n)
    (hE : E ∈ Ilam lam) :
    ((univ.filter (fun q : Cfg N =>
      E ∈ resSet (potZ lam q) lam a L i endL endR (Real.exp (-(κ * n))))).card : ℝ) ≤
      2 * C1 * Real.exp (-c1 * n) * 2 ^ N := by
  have hC1 := C1_nonneg hEB
  set θ := Real.exp (-(κ * n)) with hθ
  have hsub : univ.filter (fun q : Cfg N => E ∈ resSet (potZ lam q) lam a L i endL endR θ) ⊆
      univ.filter (fun q : Cfg N => endR = true ∧
        θ * |dirL (potZ lam q) E a (L + 1)| ≤ |dirL (potZ lam q) E a (i + 1)|) ∪
      univ.filter (fun q : Cfg N => endL = true ∧
        θ * |dirL (potZ lam q) E a (L + 1)| ≤ |dirR (potZ lam q) E (a + L - 1) (L - i)|) := by
    intro q hq
    simp only [mem_filter, mem_univ, true_and, mem_union] at hq ⊢
    simp only [resSet, Set.mem_ofPred_eq] at hq
    exact hq.2
  have hRb : ((univ.filter (fun q : Cfg N => endR = true ∧
      θ * |dirL (potZ lam q) E a (L + 1)| ≤ |dirL (potZ lam q) E a (i + 1)|)).card : ℝ) ≤
      C1 * Real.exp (-c1 * n) * 2 ^ N := by
    by_cases hR' : endR = true
    · have hn' := hR hR'
      have hiL : i + 1 + n = L := by omega
      refine end_bound lam κ C1 c1 m0 hEB E hE n (zsL N (a + i + 1) n) (nodup_zsL _ _)
        (length_zsL (by omega) (by omega)) hn
        (univ.filter (fun z : Fin N => a ≤ (z : ℤ) ∧ (z : ℤ) ≤ a + i)) ?_
        (fun q => ![dirL (potZ lam q) E a (i + 2), dirL (potZ lam q) E a (i + 1)]) ?_ ?_ _ ?_
      · intro z hz hzS
        have := mem_zsL hz
        simp only [mem_filter, mem_univ, true_and] at hzS
        omega
      · intro q q' hqq'
        have h1 : dirL (potZ lam q) E a (i + 2) = dirL (potZ lam q') E a (i + 2) :=
          dirL_congr _ _ E a (i + 2) (fun j hj => potZ_congr lam (fun z hz =>
            hqq' z (by simp only [mem_filter, mem_univ, true_and]; omega)))
        have h2 : dirL (potZ lam q) E a (i + 1) = dirL (potZ lam q') E a (i + 1) :=
          dirL_congr _ _ E a (i + 1) (fun j hj => potZ_congr lam (fun z hz =>
            hqq' z (by simp only [mem_filter, mem_univ, true_and]; omega)))
        rw [h1, h2]
      · intro q
        have := dirL_pair_ne_zero (potZ lam q) E a (i + 1)
        simpa using this
      · intro q hq
        obtain ⟨_, hq⟩ := hq
        rw [box_right (potZ lam q) E a L i hi, hn',
          tprodL_eq_tlist lam E q (by omega) (by omega)] at hq
        simpa using hq
    · simp only [hR', Bool.false_eq_true, false_and, filter_false, card_empty, Nat.cast_zero]
      positivity
  have hLb : ((univ.filter (fun q : Cfg N => endL = true ∧
      θ * |dirL (potZ lam q) E a (L + 1)| ≤ |dirR (potZ lam q) E (a + L - 1) (L - i)|)).card : ℝ) ≤
      C1 * Real.exp (-c1 * n) * 2 ^ N := by
    by_cases hL' : endL = true
    · have hin := hLe hL'
      refine end_bound lam κ C1 c1 m0 hEB E hE n (zsL N a i).reverse
        (List.nodup_reverse.2 (nodup_zsL _ _))
        (by rw [List.length_reverse, length_zsL ha (by omega)]; exact hin) hn
        (univ.filter (fun z : Fin N => a + i ≤ (z : ℤ) ∧ (z : ℤ) ≤ a + L - 1)) ?_
        (fun q => ![dirR (potZ lam q) E (a + L - 1) (L - i + 1),
          dirR (potZ lam q) E (a + L - 1) (L - i)]) ?_ ?_ _ ?_
      · intro z hz hzS
        have := mem_zsL (List.mem_reverse.1 hz)
        simp only [mem_filter, mem_univ, true_and] at hzS
        omega
      · intro q q' hqq'
        have h1 : dirR (potZ lam q) E (a + L - 1) (L - i + 1) =
            dirR (potZ lam q') E (a + L - 1) (L - i + 1) :=
          dirR_congr _ _ E (a + L - 1) (L - i + 1) (fun j hj => potZ_congr lam (fun z hz =>
            hqq' z (by simp only [mem_filter, mem_univ, true_and]; omega)))
        have h2 : dirR (potZ lam q) E (a + L - 1) (L - i) =
            dirR (potZ lam q') E (a + L - 1) (L - i) :=
          dirR_congr _ _ E (a + L - 1) (L - i) (fun j hj => potZ_congr lam (fun z hz =>
            hqq' z (by simp only [mem_filter, mem_univ, true_and]; omega)))
        rw [h1, h2]
      · intro q
        have := dirR_pair_ne_zero (potZ lam q) E (a + L - 1) (L - i)
        simpa using this
      · intro q hq
        obtain ⟨_, hq⟩ := hq
        rw [box_left (potZ lam q) E a L i hi, tprodR_eq_tlist lam E q ha (by omega)] at hq
        simpa using hq
    · simp only [hL', Bool.false_eq_true, false_and, filter_false, card_empty, Nat.cast_zero]
      positivity
  calc ((univ.filter (fun q : Cfg N =>
        E ∈ resSet (potZ lam q) lam a L i endL endR θ)).card : ℝ)
      ≤ ((univ.filter (fun q : Cfg N => endR = true ∧
        θ * |dirL (potZ lam q) E a (L + 1)| ≤ |dirL (potZ lam q) E a (i + 1)|) ∪
      univ.filter (fun q : Cfg N => endL = true ∧
        θ * |dirL (potZ lam q) E a (L + 1)| ≤
          |dirR (potZ lam q) E (a + L - 1) (L - i)|)).card : ℝ) := by
        exact_mod_cast card_le_card hsub
    _ ≤ ((univ.filter (fun q : Cfg N => endR = true ∧
        θ * |dirL (potZ lam q) E a (L + 1)| ≤ |dirL (potZ lam q) E a (i + 1)|)).card : ℝ) +
        ((univ.filter (fun q : Cfg N => endL = true ∧
        θ * |dirL (potZ lam q) E a (L + 1)| ≤
          |dirR (potZ lam q) E (a + L - 1) (L - i)|)).card : ℝ) := by
        exact_mod_cast card_union_le _ _
    _ ≤ _ := by linarith

open Classical in
/-- **Two-box estimate** (paper: Lemma `lem:double`, Eq. `eq:double`): for two disjoint boxes of
length at most `2n + 1` as in `single_box`, the probability that they have a common resonant energy
is at most `4 (8n+6) C₁ e^{-c₁ n}`. -/
theorem two_box (lam κ C1 c1 : ℝ) (m0 : ℕ) (hEB : EntryBound lam κ C1 c1 m0) (hC1 : 0 ≤ C1)
    (a1 : ℤ) (L1 i1 : ℕ) (endL1 endR1 : Bool) (a2 : ℤ) (L2 i2 : ℕ) (endL2 endR2 : Bool) (n : ℕ)
    (ha1 : 0 ≤ a1) (h12 : a1 + L1 ≤ a2) (ha2 : a2 + L2 ≤ N) (hi1 : i1 < L1) (hi2 : i2 < L2)
    (hL1 : L1 ≤ 2 * n + 1) (hL2 : L2 ≤ 2 * n + 1) (hn : m0 ≤ n)
    (hR1 : endR1 = true → L1 - 1 - i1 = n) (hLe1 : endL1 = true → i1 = n)
    (hR2 : endR2 = true → L2 - 1 - i2 = n) (hLe2 : endL2 = true → i2 = n) :
    ((univ.filter (fun q : Cfg N =>
      (resSet (potZ lam q) lam a1 L1 i1 endL1 endR1 (Real.exp (-(κ * n))) ∩
        resSet (potZ lam q) lam a2 L2 i2 endL2 endR2 (Real.exp (-(κ * n)))).Nonempty)).card : ℝ) ≤
      4 * (8 * n + 6) * C1 * Real.exp (-c1 * n) * 2 ^ N := by
  set θ := Real.exp (-(κ * n)) with hθ
  have hθpos : 0 < θ := Real.exp_pos _
  set B := 2 * C1 * Real.exp (-c1 * n) * 2 ^ N with hB
  have hB0 : 0 ≤ B := by positivity
  have ha20 : 0 ≤ a2 := by omega
  have ha1N : a1 + L1 ≤ N := by omega
  have hs1 : ∀ z : ℝ, ((univ.filter (fun q : Cfg N =>
      z ∈ resSet (potZ lam q) lam a1 L1 i1 endL1 endR1 θ)).card : ℝ) ≤ B := by
    intro z
    by_cases hz : z ∈ Ilam lam
    · exact single_box lam κ C1 c1 m0 hEB a1 L1 i1 n endL1 endR1 z ha1 ha1N hi1 hn hR1 hLe1 hz
    · have : univ.filter (fun q : Cfg N =>
          z ∈ resSet (potZ lam q) lam a1 L1 i1 endL1 endR1 θ) = ∅ := by
        ext q
        simp only [mem_filter, mem_univ, true_and, Finset.notMem_empty, iff_false]
        intro h; exact hz h.1
      rw [this, card_empty, Nat.cast_zero]; exact hB0
  have hs2 : ∀ z : ℝ, ((univ.filter (fun q : Cfg N =>
      z ∈ resSet (potZ lam q) lam a2 L2 i2 endL2 endR2 θ)).card : ℝ) ≤ B := by
    intro z
    by_cases hz : z ∈ Ilam lam
    · exact single_box lam κ C1 c1 m0 hEB a2 L2 i2 n endL2 endR2 z ha20 ha2 hi2 hn hR2 hLe2 hz
    · have : univ.filter (fun q : Cfg N =>
          z ∈ resSet (potZ lam q) lam a2 L2 i2 endL2 endR2 θ) = ∅ := by
        ext q
        simp only [mem_filter, mem_univ, true_and, Finset.notMem_empty, iff_false]
        intro h; exact hz h.1
      rw [this, card_empty, Nat.cast_zero]; exact hB0
  have hsub : univ.filter (fun q : Cfg N =>
      (resSet (potZ lam q) lam a1 L1 i1 endL1 endR1 θ ∩
        resSet (potZ lam q) lam a2 L2 i2 endL2 endR2 θ).Nonempty) ⊆
      univ.filter (fun q : Cfg N => ∃ z ∈ resEnds (potZ lam q) lam a1 L1 i1 θ,
        z ∈ resSet (potZ lam q) lam a2 L2 i2 endL2 endR2 θ) ∪
      univ.filter (fun q : Cfg N => ∃ z ∈ resEnds (potZ lam q) lam a2 L2 i2 θ,
        z ∈ resSet (potZ lam q) lam a1 L1 i1 endL1 endR1 θ) := by
    intro q hq
    simp only [mem_filter, mem_univ, true_and, mem_union] at hq ⊢
    obtain ⟨E, hE1, hE2⟩ := hq
    obtain ⟨z1, hz1, hz1E, hI1⟩ :=
      exists_resEnd (potZ lam q) lam a1 L1 i1 endL1 endR1 θ hθpos hi1 E hE1
    obtain ⟨z2, hz2, hz2E, hI2⟩ :=
      exists_resEnd (potZ lam q) lam a2 L2 i2 endL2 endR2 θ hθpos hi2 E hE2
    rcases le_total z2 z1 with h | h
    · exact Or.inl ⟨z1, hz1, hI2 ⟨h, hz1E⟩⟩
    · exact Or.inr ⟨z2, hz2, hI1 ⟨h, hz2E⟩⟩
  have hp1 : ((univ.filter (fun q : Cfg N => ∃ z ∈ resEnds (potZ lam q) lam a1 L1 i1 θ,
      z ∈ resSet (potZ lam q) lam a2 L2 i2 endL2 endR2 θ)).card : ℝ) ≤
      ((8 * n + 6 : ℕ) : ℝ) * B := by
    refine card_filter_indep_union
      (univ.filter (fun z : Fin N => a1 ≤ (z : ℤ) ∧ (z : ℤ) < a1 + L1))
      (fun q => resEnds (potZ lam q) lam a1 L1 i1 θ)
      (fun z q => z ∈ resSet (potZ lam q) lam a2 L2 i2 endL2 endR2 θ) (8 * n + 6) B hB0
      ?_ ?_ ?_ hs2 _ (fun q hq => hq)
    · intro q q' h
      exact resEnds_congr _ _ lam a1 L1 i1 θ hi1 (fun j hj => potZ_congr lam (fun z hz =>
        h z (by simp only [mem_filter, mem_univ, true_and]; omega)))
    · intro q
      exact (card_resEnds _ lam a1 L1 i1 θ hi1).trans (by omega)
    · intro z q q' h
      rw [resSet_congr (potZ lam q) (potZ lam q') lam a2 L2 i2 endL2 endR2 θ hi2
        (fun j hj => potZ_congr lam (fun w hw =>
          h w (by simp only [mem_filter, mem_univ, true_and]; omega)))]
  have hp2 : ((univ.filter (fun q : Cfg N => ∃ z ∈ resEnds (potZ lam q) lam a2 L2 i2 θ,
      z ∈ resSet (potZ lam q) lam a1 L1 i1 endL1 endR1 θ)).card : ℝ) ≤
      ((8 * n + 6 : ℕ) : ℝ) * B := by
    refine card_filter_indep_union
      (univ.filter (fun z : Fin N => a2 ≤ (z : ℤ) ∧ (z : ℤ) < a2 + L2))
      (fun q => resEnds (potZ lam q) lam a2 L2 i2 θ)
      (fun z q => z ∈ resSet (potZ lam q) lam a1 L1 i1 endL1 endR1 θ) (8 * n + 6) B hB0
      ?_ ?_ ?_ hs1 _ (fun q hq => hq)
    · intro q q' h
      exact resEnds_congr _ _ lam a2 L2 i2 θ hi2 (fun j hj => potZ_congr lam (fun z hz =>
        h z (by simp only [mem_filter, mem_univ, true_and]; omega)))
    · intro q
      exact (card_resEnds _ lam a2 L2 i2 θ hi2).trans (by omega)
    · intro z q q' h
      rw [resSet_congr (potZ lam q) (potZ lam q') lam a1 L1 i1 endL1 endR1 θ hi1
        (fun j hj => potZ_congr lam (fun w hw =>
          h w (by simp only [mem_filter, mem_univ, true_and]; omega)))]
  have hcard := (card_le_card hsub).trans (card_union_le _ _)
  have hcast : ((8 * n + 6 : ℕ) : ℝ) = 8 * n + 6 := by push_cast; ring
  rw [hcast] at hp1 hp2
  calc ((univ.filter (fun q : Cfg N =>
      (resSet (potZ lam q) lam a1 L1 i1 endL1 endR1 θ ∩
        resSet (potZ lam q) lam a2 L2 i2 endL2 endR2 θ).Nonempty)).card : ℝ)
      ≤ ((univ.filter (fun q : Cfg N => ∃ z ∈ resEnds (potZ lam q) lam a1 L1 i1 θ,
          z ∈ resSet (potZ lam q) lam a2 L2 i2 endL2 endR2 θ)).card : ℝ) +
        ((univ.filter (fun q : Cfg N => ∃ z ∈ resEnds (potZ lam q) lam a2 L2 i2 θ,
          z ∈ resSet (potZ lam q) lam a1 L1 i1 endL1 endR1 θ)).card : ℝ) := by
        exact_mod_cast hcard
    _ ≤ (8 * n + 6) * B + (8 * n + 6) * B := add_le_add hp1 hp2
    _ = _ := by rw [hB]; ring

end DFL
