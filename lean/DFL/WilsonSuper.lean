import DFL.Wilson
import DFL.QLMSuper

/-!
# Wilson links: superpositions of background-charge sectors (paper: Remark `rem:wilson`(a))

A finite superposition of sectors is described by a finite family of background charges
`Q : ι → (Fin L → ℤ)` (distinct sectors: `Q` injective; the results hold for any family) and the
components `Ψ₀ a : Cfg L → ℂ` of the state in the sectors `Q a` (a Fock vector `|n⟩` of sector `a`
stands for `|n⟩ ⊗ |E(n)⟩`). The Hamiltonian is block diagonal, so the component in sector `a`
evolves with the sector Hamiltonian `schwinger … (Q a)` (`wStateS`), the energy is the sum of the
sector energies (`energyDensityS`), and every quantity diagonal in the occupations and fields is the
sum over the sectors.

* `wilson_field_bound_super`: Proposition `prop:wilson`(b) holds for superpositions;
* `wilson_numberEntropy_le_super`: `S_N(Ψ(t); [1,ℓ]) ≤ H_{Ψ(t)}(E_ℓ) + H_{Ψ₀}(UQ_{[1,ℓ]})`, by the summed
  Gauss law and subadditivity; the second term is independent of time (`bgEntropyS_wStateS`);
* `lawEntropy_le_massey`: Lemma `lem:intentropy` for integer-valued laws on a finite space;
* `wilson_numberEntropy_avg_super`: the time-dependent part `S_N - H_{Ψ₀}(UQ)` is bounded on average
  over the cuts, uniformly in `t` and `L`;
* `bgEntropyS_bare_vacuum`, `bgEntropyS_le_log_three`: for bare-vacuum initial states
  (`ρ_j Ψ₀ = 0`) the second term is the entropy of the initial field on the cut link, at most `ln 3`
  if the initial fields take values in `{-1, 0, 1}`.
-/

open Matrix Finset Real

noncomputable section

namespace DFL

variable {L : ℕ} {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The Shannon entropy of the law of an integer-valued function `f` on a finite configuration
space in the state `v` (weights `|v c|²`). -/
def lawEntropy {κ : Type*} [Fintype κ] (f : κ → ℤ) (v : κ → ℂ) : ℝ :=
  ∑ y ∈ univ.image f, Real.negMulLog (∑ c ∈ univ.filter (fun c => f c = y), Complex.normSq (v c))

/-- The state at time `t` of a superposition of the sectors `Q a`: each component evolves with its
sector Hamiltonian. -/
def wStateS (w m g θ : ℝ) (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ₀ : ι → Cfg L → ℂ) (t : ℝ) :
    ι → Cfg L → ℂ :=
  fun a => wState w m g θ e₀ (Q a) (Ψ₀ a) t

/-- The energy density `⟨Ψ₀, HΨ₀⟩ / L = ∑_a ⟨Ψ₀^a, H_{Q_a} Ψ₀^a⟩ / L` of a superposition. -/
def energyDensityS (w m g θ : ℝ) (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ₀ : ι → Cfg L → ℂ) : ℝ :=
  (∑ a, (expect (cplx (schwinger w m g θ e₀ (Q a))) (Ψ₀ a)).re) / L

/-- `⟨Ψ, (E_j + ϑ)² Ψ⟩` for a superposition of sectors. -/
def fieldSqS (θ : ℝ) (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ : ι → Cfg L → ℂ) (j : Fin L) : ℝ :=
  ∑ a, fieldSq θ e₀ (Q a) (Ψ a) j

/-- The law of `N_R` in a superposition of sectors. -/
def numDistS (R : Finset (Fin L)) (Ψ : ι → Cfg L → ℂ) (k : ℕ) : ℝ :=
  ∑ a, numDist R (Ψ a) k

/-- The background charge `UQ_{[1,ℓ]} = ∑_{k ≤ ℓ} Q_k` (paper indexing; Lean sites `< ℓ`). -/
def bgSum (Q : Fin L → ℤ) (ℓ : ℕ) : ℤ :=
  ∑ k ∈ univ.filter (fun k : Fin L => k.val < ℓ), Q k

/-- `H_Ψ(E_j)`: the Shannon entropy of the law of the field on the link `j` in a superposition. -/
def fieldEntropyS (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ : ι → Cfg L → ℂ) (j : Fin L) : ℝ :=
  lawEntropy (fun c : ι × Cfg L => wField e₀ (Q c.1) c.2 j) (fun c => Ψ c.1 c.2)

/-- `H_Ψ(UQ_{[1,ℓ]})`: the Shannon entropy of the law of the background charge of `[1, ℓ]`. -/
def bgEntropyS (Q : ι → Fin L → ℤ) (Ψ : ι → Cfg L → ℂ) (ℓ : ℕ) : ℝ :=
  lawEntropy (fun c : ι × Cfg L => bgSum (Q c.1) ℓ) (fun c => Ψ c.1 c.2)

/-! ### Auxiliary lemmas: laws of functions on a finite space -/

namespace WilsonSuperAux

section Laws

variable {κ β : Type*} [Fintype κ] [DecidableEq β]

/-- The mass `μ{F = y} = ∑_{F c = y} μ c` of a level set. -/
def mass (F : κ → β) (μ : κ → ℝ) (y : β) : ℝ := ∑ c ∈ univ.filter (fun c => F c = y), μ c

/-- The entropy `∑_{y ∈ s} η(μ{F = y})` of the law of `F`, summed over the values in `s`. -/
def ent (s : Finset β) (F : κ → β) (μ : κ → ℝ) : ℝ := ∑ y ∈ s, Real.negMulLog (mass F μ y)

omit [DecidableEq β] in
/-- `lawEntropy` is `ent` over the image. -/
theorem lawEntropy_eq_ent (f : κ → ℤ) (v : κ → ℂ) :
    lawEntropy f v = ent (univ.image f) f (fun c => Complex.normSq (v c)) := rfl

/-- The masses of a nonnegative weight are nonnegative. -/
theorem mass_nonneg (F : κ → β) {μ : κ → ℝ} (hμ : ∀ c, 0 ≤ μ c) (y : β) : 0 ≤ mass F μ y :=
  Finset.sum_nonneg (fun c _ => hμ c)

/-- The mass of a value taken at no point of the support vanishes. -/
theorem mass_eq_zero (F : κ → β) (μ : κ → ℝ) (y : β) (h : ∀ c, μ c ≠ 0 → F c ≠ y) :
    mass F μ y = 0 := by
  apply Finset.sum_eq_zero
  intro c hc
  by_contra hne
  exact h c hne (Finset.mem_filter.1 hc).2

/-- The masses only depend on `F` on the support of `μ`. -/
theorem mass_congr (F F' : κ → β) (μ : κ → ℝ) (h : ∀ c, μ c ≠ 0 → F c = F' c) (y : β) :
    mass F μ y = mass F' μ y := by
  unfold mass
  rw [Finset.sum_filter, Finset.sum_filter]
  refine Finset.sum_congr rfl (fun c _ => ?_)
  by_cases hc : μ c = 0
  · simp [hc]
  · rw [h c hc]

/-- Composing with an injective map relabels the masses. -/
theorem mass_comp_inj {β' : Type*} [DecidableEq β'] (F : κ → β) (σ : β → β')
    (hσ : Function.Injective σ) (μ : κ → ℝ) (y : β) :
    mass (fun c => σ (F c)) μ (σ y) = mass F μ y := by
  unfold mass
  congr 1
  ext c
  simp [hσ.eq_iff]

/-- Sums of `G (μ{F = y})` with `G 0 = 0` do not depend on the set of values summed over, as long
as it contains the values taken on the support. -/
theorem sum_mass_congr {G : ℝ → ℝ} (hG : G 0 = 0) (F : κ → β) (μ : κ → ℝ) {s t : Finset β}
    (hs : ∀ c, μ c ≠ 0 → F c ∈ s) (ht : ∀ c, μ c ≠ 0 → F c ∈ t) :
    ∑ y ∈ s, G (mass F μ y) = ∑ y ∈ t, G (mass F μ y) := by
  have key : ∀ u : Finset β, (∀ c, μ c ≠ 0 → F c ∈ u) → ∀ y, y ∉ u → G (mass F μ y) = 0 := by
    intro u hu y hy
    rw [mass_eq_zero F μ y (fun c hc hcy => hy (hcy ▸ hu c hc)), hG]
  have h1 : ∑ y ∈ s ∩ t, G (mass F μ y) = ∑ y ∈ s, G (mass F μ y) :=
    Finset.sum_subset Finset.inter_subset_left (fun y hys hyst =>
      key t ht y (fun hyt => hyst (Finset.mem_inter.2 ⟨hys, hyt⟩)))
  have h2 : ∑ y ∈ s ∩ t, G (mass F μ y) = ∑ y ∈ t, G (mass F μ y) :=
    Finset.sum_subset Finset.inter_subset_right (fun y hyt hyst =>
      key s hs y (fun hys => hyst (Finset.mem_inter.2 ⟨hys, hyt⟩)))
  rw [← h1, h2]

/-- Integrating a function of `F` against the law of `F`. -/
theorem sum_mass_mul (F : κ → β) (μ : κ → ℝ) {s : Finset β} (hs : ∀ c, F c ∈ s) (h : β → ℝ) :
    ∑ y ∈ s, mass F μ y * h y = ∑ c, μ c * h (F c) := by
  unfold mass
  simp_rw [Finset.sum_mul]
  rw [← Finset.sum_fiberwise_of_maps_to (s := univ) (t := s) (g := F) (fun c _ => hs c)]
  refine Finset.sum_congr rfl (fun y _ => Finset.sum_congr rfl (fun c hc => ?_))
  rw [(Finset.mem_filter.1 hc).2]

/-- The entropy is the expectation of `-log` of the mass of the observed value. -/
theorem ent_eq_sum (F : κ → β) (μ : κ → ℝ) {s : Finset β} (hs : ∀ c, F c ∈ s) :
    ent s F μ = ∑ c, μ c * -Real.log (mass F μ (F c)) := by
  unfold ent
  rw [← sum_mass_mul F μ hs (fun y => -Real.log (mass F μ y))]
  refine Finset.sum_congr rfl (fun y _ => ?_)
  rw [Real.negMulLog]
  ring

/-- Relabelling the values by an injective map does not change the entropy. -/
theorem ent_comp_inj {β' : Type*} [DecidableEq β'] (s : Finset β) (F : κ → β) (σ : β → β')
    (hσ : Function.Injective σ) (μ : κ → ℝ) :
    ent (s.image σ) (fun c => σ (F c)) μ = ent s F μ := by
  unfold ent
  rw [Finset.sum_image (fun x _ y _ hxy => hσ hxy)]
  refine Finset.sum_congr rfl (fun y _ => ?_)
  rw [mass_comp_inj F σ hσ μ y]

/-- **Grouping**: the entropy of a function of `H` is at most the entropy of `H`. -/
theorem ent_le_of_comp {β' : Type*} [DecidableEq β'] (N : κ → β') (H : κ → β) (φ : β → β')
    (hN : ∀ c, N c = φ (H c)) (μ : κ → ℝ) (hμ : ∀ c, 0 ≤ μ c) {K : Finset β'}
    (hK : ∀ c, N c ∈ K) {P : Finset β} (hP : ∀ c, H c ∈ P) : ent K N μ ≤ ent P H μ := by
  rw [ent_eq_sum N μ hK, ent_eq_sum H μ hP]
  refine Finset.sum_le_sum (fun c _ => ?_)
  rcases (hμ c).eq_or_lt with h | h
  · rw [← h, zero_mul, zero_mul]
  · have h1 : μ c ≤ mass H μ (H c) :=
      Finset.single_le_sum (f := μ) (fun c' _ => hμ c')
        (Finset.mem_filter.2 ⟨Finset.mem_univ c, rfl⟩)
    have h2 : mass H μ (H c) ≤ mass N μ (N c) := by
      apply Finset.sum_le_sum_of_subset_of_nonneg
      · intro c' hc'
        simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hc' ⊢
        rw [hN, hN, hc']
      · intro c' _ _
        exact hμ c'
    have h3 : Real.log (mass H μ (H c)) ≤ Real.log (mass N μ (N c)) :=
      Real.log_le_log (by linarith) h2
    nlinarith

end Laws

/-- **Subadditivity** for a joint law `m` on `X × Y` with marginals `p`, `q` of total mass one:
`∑ η(m) ≤ ∑ η(p) + ∑ η(q)`. -/
theorem sum_negMulLog_joint_le {β₁ β₂ : Type*} (X : Finset β₁) (Y : Finset β₂)
    (m : β₁ → β₂ → ℝ) (p : β₁ → ℝ) (q : β₂ → ℝ) (hm0 : ∀ x y, 0 ≤ m x y)
    (hmx : ∀ x ∈ X, ∑ y ∈ Y, m x y = p x) (hmy : ∀ y ∈ Y, ∑ x ∈ X, m x y = q y)
    (hp1 : ∑ x ∈ X, p x = 1) (hq1 : ∑ y ∈ Y, q y = 1) :
    ∑ x ∈ X, ∑ y ∈ Y, Real.negMulLog (m x y) ≤
      ∑ x ∈ X, Real.negMulLog (p x) + ∑ y ∈ Y, Real.negMulLog (q y) := by
  have hp0 : ∀ x ∈ X, 0 ≤ p x := fun x hx =>
    (hmx x hx) ▸ Finset.sum_nonneg (fun y _ => hm0 x y)
  have hq0 : ∀ y ∈ Y, 0 ≤ q y := fun y hy =>
    (hmy y hy) ▸ Finset.sum_nonneg (fun x _ => hm0 x y)
  -- termwise Gibbs inequality against `p x * q y`
  have hb : ∀ x ∈ X, ∀ y ∈ Y, Real.negMulLog (m x y) ≤
      p x * q y - m x y - m x y * Real.log (p x) - m x y * Real.log (q y) := by
    intro x hx y hy
    rcases (hm0 x y).eq_or_lt with h | h
    · rw [← h, Real.negMulLog_zero]
      simp only [zero_mul, sub_zero]
      exact mul_nonneg (hp0 x hx) (hq0 y hy)
    · have hmp : m x y ≤ p x := by
        rw [← hmx x hx]
        exact Finset.single_le_sum (fun y' _ => hm0 x y') hy
      have hmq : m x y ≤ q y := by
        rw [← hmy y hy]
        exact Finset.single_le_sum (fun x' _ => hm0 x' y) hx
      have hpx : 0 < p x := lt_of_lt_of_le h hmp
      have hqy : 0 < q y := lt_of_lt_of_le h hmq
      have := QLMSuperAux.negMulLog_le_gibbs (hm0 x y) (mul_pos hpx hqy)
      rw [Real.log_mul hpx.ne' hqy.ne'] at this
      linarith
  have e1 : ∑ x ∈ X, ∑ y ∈ Y, p x * q y = 1 := by
    simp_rw [← Finset.mul_sum, hq1, mul_one, hp1]
  have e2 : ∑ x ∈ X, ∑ y ∈ Y, m x y = 1 := by
    rw [Finset.sum_congr rfl hmx, hp1]
  have e3 : ∑ x ∈ X, ∑ y ∈ Y, m x y * Real.log (p x) = ∑ x ∈ X, p x * Real.log (p x) := by
    refine Finset.sum_congr rfl (fun x hx => ?_)
    rw [← Finset.sum_mul, hmx x hx]
  have e4 : ∑ x ∈ X, ∑ y ∈ Y, m x y * Real.log (q y) = ∑ y ∈ Y, q y * Real.log (q y) := by
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl (fun y hy => ?_)
    rw [← Finset.sum_mul, hmy y hy]
  have e5 : ∀ z : ℝ, Real.negMulLog z = -(z * Real.log z) := fun z => by
    rw [Real.negMulLog]
    ring
  calc ∑ x ∈ X, ∑ y ∈ Y, Real.negMulLog (m x y)
      ≤ ∑ x ∈ X, ∑ y ∈ Y,
          (p x * q y - m x y - m x y * Real.log (p x) - m x y * Real.log (q y)) :=
        Finset.sum_le_sum (fun x hx => Finset.sum_le_sum (fun y hy => hb x hx y hy))
    _ = ∑ x ∈ X, ∑ y ∈ Y, p x * q y - ∑ x ∈ X, ∑ y ∈ Y, m x y -
          ∑ x ∈ X, ∑ y ∈ Y, m x y * Real.log (p x) -
          ∑ x ∈ X, ∑ y ∈ Y, m x y * Real.log (q y) := by
        simp only [Finset.sum_sub_distrib]
    _ = ∑ x ∈ X, Real.negMulLog (p x) + ∑ y ∈ Y, Real.negMulLog (q y) := by
        rw [e1, e2, e3, e4]
        simp only [e5, Finset.sum_neg_distrib]
        ring

section Joint

variable {κ : Type*} [Fintype κ]

/-- **Subadditivity**: the joint entropy of `(F, G)` is at most `H(F) + H(G)` for a probability
weight `μ`. -/
theorem ent_prod_le {β₁ β₂ : Type*} [DecidableEq β₁] [DecidableEq β₂] (F : κ → β₁)
    (G : κ → β₂) (μ : κ → ℝ) (hμ : ∀ c, 0 ≤ μ c) (h1 : ∑ c, μ c = 1) :
    ent (univ.image F ×ˢ univ.image G) (fun c => (F c, G c)) μ ≤
      ent (univ.image F) F μ + ent (univ.image G) G μ := by
  have hFX : ∀ c, F c ∈ univ.image F := fun c => Finset.mem_image_of_mem F (Finset.mem_univ c)
  have hGY : ∀ c, G c ∈ univ.image G := fun c => Finset.mem_image_of_mem G (Finset.mem_univ c)
  unfold ent
  rw [Finset.sum_product]
  refine sum_negMulLog_joint_le _ _ (fun x y => mass (fun c => (F c, G c)) μ (x, y))
    (mass F μ) (mass G μ) (fun x y => mass_nonneg _ hμ _) ?_ ?_ ?_ ?_
  · intro x _
    unfold mass
    rw [← Finset.sum_fiberwise_of_maps_to (s := univ.filter (fun c => F c = x))
      (t := univ.image G) (g := G) (fun c _ => hGY c)]
    refine Finset.sum_congr rfl (fun y _ => ?_)
    rw [Finset.filter_filter]
    refine Finset.sum_congr (Finset.filter_congr (fun c _ => ?_)) (fun _ _ => rfl)
    rw [Prod.mk.injEq]
  · intro y _
    unfold mass
    rw [← Finset.sum_fiberwise_of_maps_to (s := univ.filter (fun c => G c = y))
      (t := univ.image F) (g := F) (fun c _ => hFX c)]
    refine Finset.sum_congr rfl (fun x _ => ?_)
    rw [Finset.filter_filter]
    refine Finset.sum_congr (Finset.filter_congr (fun c _ => ?_)) (fun _ _ => rfl)
    rw [Prod.mk.injEq]
    exact and_comm
  · have := sum_mass_mul F μ hFX (fun _ => 1)
    simp only [mul_one] at this
    rw [this, h1]
  · have := sum_mass_mul G μ hGY (fun _ => 1)
    simp only [mul_one] at this
    rw [this, h1]

/-- **Massey's bound** (`lem:intentropy`) for the law of an integer-valued function under a
probability weight. -/
theorem ent_le_massey (f : κ → ℤ) (μ : κ → ℝ) (hμ : ∀ c, 0 ≤ μ c) (h1 : ∑ c, μ c = 1)
    (θ : ℝ) :
    ent (univ.image f) f μ ≤ 1 / 2 * Real.log (2 * π * exp 1 *
      (∑ c, μ c * ((f c : ℝ) + θ) ^ 2 + 1 / 12)) := by
  have hne : (univ.image f).Nonempty := by
    by_contra h
    rw [Finset.not_nonempty_iff_eq_empty, Finset.image_eq_empty, Finset.univ_eq_empty_iff] at h
    rw [Finset.univ_eq_empty, Finset.sum_empty] at h1
    exact zero_ne_one h1
  obtain ⟨y₀, hy₀⟩ : ∃ y₀, y₀ = (univ.image f).min' hne := ⟨_, rfl⟩
  obtain ⟨M, hM⟩ : ∃ M : ℕ, M = ((univ.image f).max' hne - y₀).toNat := ⟨_, rfl⟩
  obtain ⟨σ, hσ⟩ : ∃ σ : ℕ → ℤ, ∀ k, σ k = y₀ + k := ⟨fun k => y₀ + k, fun _ => rfl⟩
  have hσinj : Function.Injective σ := fun a b h => by
    rw [hσ, hσ] at h
    omega
  have hmem : ∀ c, f c ∈ (range (M + 1)).image σ := by
    intro c
    have hc := Finset.mem_image_of_mem f (Finset.mem_univ c)
    have hlo : y₀ ≤ f c := hy₀ ▸ Finset.min'_le _ _ hc
    have hhi : f c ≤ (univ.image f).max' hne := Finset.le_max' _ _ hc
    refine Finset.mem_image.2 ⟨(f c - y₀).toNat, Finset.mem_range.2 (by omega), ?_⟩
    rw [hσ]
    omega
  have hp : IsProb M (fun k => mass f μ (σ k)) := by
    refine ⟨fun k => mass_nonneg f hμ _, ?_⟩
    rw [← Finset.sum_image (f := fun y => mass f μ y) (fun x _ y _ hxy => hσinj hxy)]
    have := sum_mass_mul f μ hmem (fun _ => 1)
    simp only [mul_one] at this
    rw [this, h1]
  have hshannon : shannon M (fun k => mass f μ (σ k)) = ent (univ.image f) f μ := by
    unfold shannon ent
    rw [← Finset.sum_image (f := fun y => Real.negMulLog (mass f μ y))
      (fun x _ y _ hxy => hσinj hxy)]
    exact sum_mass_congr Real.negMulLog_zero f μ (fun c _ => hmem c)
      (fun c _ => Finset.mem_image_of_mem f (Finset.mem_univ c))
  have hvar : pvar M (fun k => mass f μ (σ k)) ≤ ∑ c, μ c * ((f c : ℝ) + θ) ^ 2 := by
    refine (pvar_le_sum_sq hp (-((y₀ : ℝ) + θ))).trans (le_of_eq ?_)
    rw [← sum_mass_mul f μ hmem (fun y => ((y : ℝ) + θ) ^ 2),
      Finset.sum_image (fun x _ y _ hxy => hσinj hxy)]
    refine Finset.sum_congr rfl (fun k _ => ?_)
    rw [hσ]
    push_cast
    ring
  rw [← hshannon]
  exact shannon_le_of_pvar_le hp hvar

/-- A law supported on a set `T` of values has entropy at most `ln |T|`. -/
theorem ent_le_log_card (f : κ → ℤ) (μ : κ → ℝ) (hμ : ∀ c, 0 ≤ μ c) (h1 : ∑ c, μ c = 1)
    (T : Finset ℤ) (hT : ∀ c, μ c ≠ 0 → f c ∈ T) :
    ent (univ.image f) f μ ≤ Real.log T.card := by
  have himg : ∀ c, μ c ≠ 0 → f c ∈ univ.image f := fun c _ =>
    Finset.mem_image_of_mem f (Finset.mem_univ c)
  have hent : ent (univ.image f) f μ = ∑ y ∈ T, Real.negMulLog (mass f μ y) := by
    unfold ent
    exact sum_mass_congr Real.negMulLog_zero f μ himg hT
  have hsum : ∑ y ∈ T, mass f μ y = 1 := by
    have h := sum_mass_congr (G := id) rfl f μ hT himg
    simp only [id] at h
    rw [h]
    have := sum_mass_mul f μ (fun c => Finset.mem_image_of_mem f (Finset.mem_univ c))
      (fun _ => 1)
    simp only [mul_one] at this
    rw [this, h1]
  rw [hent]
  have := QLMSuperAux.sum_negMulLog_le T (mass f μ) (fun y _ => mass_nonneg f hμ y)
    (L := T.card) le_rfl
  rw [hsum, Real.negMulLog_one, one_mul, zero_add] at this
  exact this

/-- If `f = σ ∘ g` on the support of `v` with `σ` injective, the laws of `f` and `g` have the same
entropy. -/
theorem lawEntropy_congr_inj (f g : κ → ℤ) (σ : ℤ → ℤ) (hσ : Function.Injective σ) (v : κ → ℂ)
    (h : ∀ c, v c ≠ 0 → f c = σ (g c)) : lawEntropy f v = lawEntropy g v := by
  have hμ : ∀ c, Complex.normSq (v c) ≠ 0 → f c = σ (g c) := fun c hc =>
    h c (fun h0 => hc (by rw [h0, map_zero]))
  rw [lawEntropy_eq_ent, lawEntropy_eq_ent, ← ent_comp_inj (univ.image g) g σ hσ]
  unfold ent
  rw [sum_mass_congr Real.negMulLog_zero f _ (s := univ.image f) (t := (univ.image g).image σ)
    (fun c _ => Finset.mem_image_of_mem f (Finset.mem_univ c))
    (fun c hc => by
      rw [hμ c hc]
      exact Finset.mem_image_of_mem σ (Finset.mem_image_of_mem g (Finset.mem_univ c)))]
  refine Finset.sum_congr rfl (fun y _ => ?_)
  rw [mass_congr f (fun c => σ (g c)) _ hμ y]

end Joint

section Wilson

/-- Homogeneous field energy bound in a sector:
`(g²/2) ∑_j ⟨v, (E_j + ϑ)² v⟩ ≤ Re ⟨v, H v⟩ + L (w + |m|) ‖v‖²`. -/
theorem field_le_energy' (w m g θ : ℝ) (hw : 0 ≤ w) (e₀ : ℤ) (Q : Fin L → ℤ) (v : Cfg L → ℂ) :
    g ^ 2 / 2 * ∑ j, fieldSq θ e₀ Q v j ≤
      (expect (cplx (schwinger w m g θ e₀ Q)) v).re +
        L * (w + |m|) * ∑ n, Complex.normSq (v n) := by
  rw [WilsonAux.re_expect_schwinger]
  have hB := mul_le_mul_of_nonneg_left (WilsonAux.re_expect_hop_le v) hw
  have hM : -(|m| * L) * ∑ n, Complex.normSq (v n) ≤
      ∑ n, (m * ∑ j : Fin L, (-1 : ℝ) ^ (j.val + 1) * occ (n j)) * Complex.normSq (v n) := by
    rw [Finset.mul_sum]
    exact Finset.sum_le_sum (fun n _ =>
      mul_le_mul_of_nonneg_right (WilsonAux.mass_ge m n) (Complex.normSq_nonneg _))
  linarith

omit [Fintype ι] [DecidableEq ι] in
/-- Each component of the superposition keeps its norm. -/
theorem sum_normSq_wStateS (w m g θ : ℝ) (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ₀ : ι → Cfg L → ℂ)
    (t : ℝ) (a : ι) :
    ∑ n, Complex.normSq (wStateS w m g θ e₀ Q Ψ₀ t a n) = ∑ n, Complex.normSq (Ψ₀ a n) := by
  show ∑ n, Complex.normSq ((evo (schwinger w m g θ e₀ (Q a)) t *ᵥ Ψ₀ a) n) = _
  exact WilsonAux.sum_normSq_mulVec _
    (evo_unitary _ (WilsonAux.schwinger_isSymm w m g θ e₀ (Q a)) t) _

omit [DecidableEq ι] in
/-- The evolved superposition stays normalized. -/
theorem sum_normSq_wStateS_eq_one (w m g θ : ℝ) (e₀ : ℤ) (Q : ι → Fin L → ℤ)
    (Ψ₀ : ι → Cfg L → ℂ) (hΨ₀ : ∑ a, ∑ n, Complex.normSq (Ψ₀ a n) = 1) (t : ℝ) :
    ∑ a, ∑ n, Complex.normSq (wStateS w m g θ e₀ Q Ψ₀ t a n) = 1 := by
  rw [← hΨ₀]
  exact Finset.sum_congr rfl (fun a _ => sum_normSq_wStateS w m g θ e₀ Q Ψ₀ t a)

omit [DecidableEq ι] in
/-- `S_N(Ψ; [1,ℓ]) ≤ H_Ψ(E_ℓ) + H_Ψ(UQ_{[1,ℓ]})` for every normalized superposition `Ψ` (summed
Gauss law, grouping and subadditivity). -/
theorem shannon_numDistS_le (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ : ι → Cfg L → ℂ)
    (hΨ : ∑ a, ∑ n, Complex.normSq (Ψ a n) = 1) (ℓ : ℕ) (hℓ : 0 < ℓ) (hℓL : ℓ ≤ L) :
    shannon L (numDistS (segment ℓ) Ψ) ≤
      fieldEntropyS e₀ Q Ψ ⟨ℓ - 1, by omega⟩ + bgEntropyS Q Ψ ℓ := by
  have hμ0 : ∀ c : ι × Cfg L, 0 ≤ Complex.normSq (Ψ c.1 c.2) := fun c => Complex.normSq_nonneg _
  have hμ1 : ∑ c : ι × Cfg L, Complex.normSq (Ψ c.1 c.2) = 1 := by
    rw [Fintype.sum_prod_type]
    exact hΨ
  -- the law of `N_{[1,ℓ]}` as a law on `ι × Cfg L`
  have hS : shannon L (numDistS (segment ℓ) Ψ) =
      ent (range (L + 1)) (fun c : ι × Cfg L => WilsonAux.countN (segment ℓ) c.2)
        (fun c => Complex.normSq (Ψ c.1 c.2)) := by
    unfold shannon ent
    refine Finset.sum_congr rfl (fun k _ => ?_)
    congr 1
    unfold numDistS numDist mass
    rw [Finset.sum_filter, Fintype.sum_prod_type]
    refine Finset.sum_congr rfl (fun a _ => ?_)
    rw [Finset.sum_filter]
    refine Finset.sum_congr rfl (fun n _ => ?_)
    simp only [WilsonAux.count_eq_countN, Nat.cast_inj]
  -- summed Gauss law: `N_{[1,ℓ]}` is a function of `(E_ℓ, UQ_{[1,ℓ]})`
  have hN : ∀ c : ι × Cfg L, WilsonAux.countN (segment ℓ) c.2 =
      (wField e₀ (Q c.1) c.2 ⟨ℓ - 1, by omega⟩ - e₀ +
        ((univ.filter (fun k : Fin L => k.val < ℓ ∧ k.val % 2 = 0)).card : ℤ) -
          bgSum (Q c.1) ℓ).toNat := by
    intro c
    have h := wField_eq_count e₀ (Q c.1) c.2 ℓ hℓ hℓL
    have hc : wConst e₀ (Q c.1) ℓ = ((e₀ - ((univ.filter (fun k : Fin L =>
        k.val < ℓ ∧ k.val % 2 = 0)).card : ℤ) + bgSum (Q c.1) ℓ : ℤ) : ℝ) := by
      unfold wConst bgSum
      push_cast
      ring
    rw [hc, WilsonAux.count_eq_countN] at h
    have h' : wField e₀ (Q c.1) c.2 ⟨ℓ - 1, by omega⟩ =
        (WilsonAux.countN (segment ℓ) c.2 : ℤ) +
          (e₀ - ((univ.filter (fun k : Fin L => k.val < ℓ ∧ k.val % 2 = 0)).card : ℤ) +
            bgSum (Q c.1) ℓ) := by
      exact_mod_cast h
    omega
  rw [hS]
  have h1 := ent_le_of_comp (fun c : ι × Cfg L => WilsonAux.countN (segment ℓ) c.2)
    (fun c : ι × Cfg L => (wField e₀ (Q c.1) c.2 ⟨ℓ - 1, by omega⟩, bgSum (Q c.1) ℓ))
    (fun p : ℤ × ℤ => (p.1 - e₀ +
      ((univ.filter (fun k : Fin L => k.val < ℓ ∧ k.val % 2 = 0)).card : ℤ) - p.2).toNat)
    hN (fun c => Complex.normSq (Ψ c.1 c.2)) hμ0 (K := range (L + 1))
    (fun c => Finset.mem_range.2 (Nat.lt_succ_of_le (WilsonAux.countN_le _ _)))
    (P := univ.image (fun c : ι × Cfg L => wField e₀ (Q c.1) c.2 ⟨ℓ - 1, by omega⟩) ×ˢ
      univ.image (fun c : ι × Cfg L => bgSum (Q c.1) ℓ))
    (fun c => Finset.mem_product.2 ⟨Finset.mem_image_of_mem _ (Finset.mem_univ c),
      Finset.mem_image_of_mem _ (Finset.mem_univ c)⟩)
  have h2 := ent_prod_le (fun c : ι × Cfg L => wField e₀ (Q c.1) c.2 ⟨ℓ - 1, by omega⟩)
    (fun c : ι × Cfg L => bgSum (Q c.1) ℓ) (fun c => Complex.normSq (Ψ c.1 c.2)) hμ0 hμ1
  have h3 : fieldEntropyS e₀ Q Ψ ⟨ℓ - 1, by omega⟩ + bgEntropyS Q Ψ ℓ =
      ent (univ.image (fun c : ι × Cfg L => wField e₀ (Q c.1) c.2 ⟨ℓ - 1, by omega⟩))
        (fun c : ι × Cfg L => wField e₀ (Q c.1) c.2 ⟨ℓ - 1, by omega⟩)
        (fun c => Complex.normSq (Ψ c.1 c.2)) +
      ent (univ.image (fun c : ι × Cfg L => bgSum (Q c.1) ℓ))
        (fun c : ι × Cfg L => bgSum (Q c.1) ℓ) (fun c => Complex.normSq (Ψ c.1 c.2)) := rfl
  rw [h3]
  exact h1.trans h2

end Wilson

end WilsonSuperAux

/-! ### Integer-valued laws -/

/-- **Lemma `lem:intentropy` for integer-valued laws**: `H(X) ≤ ½ ln(2πe(E(X + ϑ)² + 1/12))`. -/
theorem lawEntropy_le_massey {κ : Type*} [Fintype κ] (f : κ → ℤ) (v : κ → ℂ)
    (hv : ∑ c, Complex.normSq (v c) = 1) (θ : ℝ) :
    lawEntropy f v ≤ 1 / 2 * Real.log (2 * π * exp 1 *
      (∑ c, Complex.normSq (v c) * ((f c : ℝ) + θ) ^ 2 + 1 / 12)) := by
  rw [WilsonSuperAux.lawEntropy_eq_ent]
  exact WilsonSuperAux.ent_le_massey f _ (fun c => Complex.normSq_nonneg _) hv θ

/-- A law supported on a set `T` of values has entropy at most `ln |T|`. -/
theorem lawEntropy_le_log_card {κ : Type*} [Fintype κ] (f : κ → ℤ) (v : κ → ℂ)
    (hv : ∑ c, Complex.normSq (v c) = 1) (T : Finset ℤ) (hT : ∀ c, v c ≠ 0 → f c ∈ T) :
    lawEntropy f v ≤ Real.log T.card := by
  rw [WilsonSuperAux.lawEntropy_eq_ent]
  refine WilsonSuperAux.ent_le_log_card f _ (fun c => Complex.normSq_nonneg _) hv T
    (fun c hc => hT c (fun h0 => hc ?_))
  simp [h0]

/-! ### Energy conservation for superpositions -/

set_option linter.unusedSectionVars false in
/-- **Proposition `prop:wilson`(b) for superpositions of sectors** (paper: Remark `rem:wilson`(a)):
`ε₀ + w + |m| ≥ 0` and `(1/L) ∑_j ⟨(E_j + ϑ)²⟩_t ≤ (2/g²)(ε₀ + w + |m|)`. -/
theorem wilson_field_bound_super (w m g θ : ℝ) (hw : 0 ≤ w) (hg : 0 < g) (e₀ : ℤ)
    (Q : ι → Fin L → ℤ) (hL : 0 < L) (Ψ₀ : ι → Cfg L → ℂ)
    (hΨ₀ : ∑ a, ∑ n, Complex.normSq (Ψ₀ a n) = 1) (t : ℝ) :
    0 ≤ energyDensityS w m g θ e₀ Q Ψ₀ + w + |m| ∧
      (1 / (L : ℝ)) * ∑ j, fieldSqS θ e₀ Q (wStateS w m g θ e₀ Q Ψ₀ t) j ≤
        2 / g ^ 2 * (energyDensityS w m g θ e₀ Q Ψ₀ + w + |m|) := by
  have hLpos : (0 : ℝ) < L := Nat.cast_pos.2 hL
  have hg2 : 0 < g ^ 2 := by positivity
  -- per sector, at time `t`: conservation of the energy and of the norm
  have hsec : ∀ a, g ^ 2 / 2 * ∑ j, fieldSq θ e₀ (Q a) (wStateS w m g θ e₀ Q Ψ₀ t a) j ≤
      (expect (cplx (schwinger w m g θ e₀ (Q a))) (Ψ₀ a)).re +
        L * (w + |m|) * ∑ n, Complex.normSq (Ψ₀ a n) := by
    intro a
    have h := WilsonSuperAux.field_le_energy' w m g θ hw e₀ (Q a) (wStateS w m g θ e₀ Q Ψ₀ t a)
    rwa [WilsonSuperAux.sum_normSq_wStateS w m g θ e₀ Q Ψ₀ t a,
      show wStateS w m g θ e₀ Q Ψ₀ t a = evo (schwinger w m g θ e₀ (Q a)) t *ᵥ Ψ₀ a from rfl,
      WilsonAux.expect_evo _ (WilsonAux.schwinger_isSymm w m g θ e₀ (Q a))] at h
  have hsum : g ^ 2 / 2 * ∑ j, fieldSqS θ e₀ Q (wStateS w m g θ e₀ Q Ψ₀ t) j ≤
      L * energyDensityS w m g θ e₀ Q Ψ₀ + L * (w + |m|) := by
    have h1 : ∑ j, fieldSqS θ e₀ Q (wStateS w m g θ e₀ Q Ψ₀ t) j =
        ∑ a, ∑ j, fieldSq θ e₀ (Q a) (wStateS w m g θ e₀ Q Ψ₀ t a) j := by
      unfold fieldSqS
      exact Finset.sum_comm
    have h2 : L * energyDensityS w m g θ e₀ Q Ψ₀ =
        ∑ a, (expect (cplx (schwinger w m g θ e₀ (Q a))) (Ψ₀ a)).re := by
      unfold energyDensityS
      field_simp
    rw [h1, h2, Finset.mul_sum]
    calc ∑ a, g ^ 2 / 2 * ∑ j, fieldSq θ e₀ (Q a) (wStateS w m g θ e₀ Q Ψ₀ t a) j
        ≤ ∑ a, ((expect (cplx (schwinger w m g θ e₀ (Q a))) (Ψ₀ a)).re +
          L * (w + |m|) * ∑ n, Complex.normSq (Ψ₀ a n)) := Finset.sum_le_sum (fun a _ => hsec a)
      _ = ∑ a, (expect (cplx (schwinger w m g θ e₀ (Q a))) (Ψ₀ a)).re + L * (w + |m|) := by
        rw [Finset.sum_add_distrib, ← Finset.mul_sum, hΨ₀, mul_one]
  have hS0 : 0 ≤ ∑ j, fieldSqS θ e₀ Q (wStateS w m g θ e₀ Q Ψ₀ t) j :=
    Finset.sum_nonneg (fun j _ =>
      Finset.sum_nonneg (fun a _ => WilsonAux.fieldSq_nonneg θ e₀ (Q a) _ j))
  set S := ∑ j, fieldSqS θ e₀ Q (wStateS w m g θ e₀ Q Ψ₀ t) j with hSdef
  set X := energyDensityS w m g θ e₀ Q Ψ₀ + w + |m| with hXdef
  have hF' : g ^ 2 / 2 * S ≤ L * X := by rw [hXdef]; linarith
  have hX : 0 ≤ X := by
    by_contra h
    have : L * X < 0 := mul_neg_of_pos_of_neg hLpos (not_le.1 h)
    nlinarith
  refine ⟨hX, ?_⟩
  rw [one_div_mul_eq_div, div_le_iff₀ hLpos, div_mul_eq_mul_div, div_mul_eq_mul_div,
    le_div_iff₀ hg2]
  nlinarith

/-! ### Number entropy of superpositions -/

set_option linter.unusedSectionVars false in
/-- The law of the background charge does not change in time. -/
theorem bgEntropyS_wStateS (w m g θ : ℝ) (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ₀ : ι → Cfg L → ℂ)
    (t : ℝ) (ℓ : ℕ) : bgEntropyS Q (wStateS w m g θ e₀ Q Ψ₀ t) ℓ = bgEntropyS Q Ψ₀ ℓ := by
  unfold bgEntropyS lawEntropy
  refine Finset.sum_congr rfl (fun y _ => ?_)
  congr 1
  rw [Finset.sum_filter, Finset.sum_filter, Fintype.sum_prod_type, Fintype.sum_prod_type]
  refine Finset.sum_congr rfl (fun a _ => ?_)
  by_cases ha : bgSum (Q a) ℓ = y
  · simp only [ha, ↓reduceIte]
    exact WilsonSuperAux.sum_normSq_wStateS w m g θ e₀ Q Ψ₀ t a
  · simp only [ha, ↓reduceIte, Finset.sum_const_zero]

/-- **Number entropy of superpositions** (paper: Remark `rem:wilson`(a)):
`S_N(Ψ(t); [1,ℓ]) ≤ H_{Ψ(t)}(E_ℓ) + H_{Ψ₀}(UQ_{[1,ℓ]})`. -/
theorem wilson_numberEntropy_le_super (w m g θ : ℝ) (e₀ : ℤ) (Q : ι → Fin L → ℤ)
    (Ψ₀ : ι → Cfg L → ℂ) (hΨ₀ : ∑ a, ∑ n, Complex.normSq (Ψ₀ a n) = 1) (t : ℝ) (ℓ : ℕ)
    (hℓ : 0 < ℓ) (hℓL : ℓ ≤ L) :
    shannon L (numDistS (segment ℓ) (wStateS w m g θ e₀ Q Ψ₀ t)) ≤
      fieldEntropyS e₀ Q (wStateS w m g θ e₀ Q Ψ₀ t) ⟨ℓ - 1, by omega⟩ + bgEntropyS Q Ψ₀ ℓ := by
  rw [← bgEntropyS_wStateS w m g θ e₀ Q Ψ₀ t ℓ]
  exact WilsonSuperAux.shannon_numDistS_le e₀ Q _
    (WilsonSuperAux.sum_normSq_wStateS_eq_one w m g θ e₀ Q Ψ₀ hΨ₀ t) ℓ hℓ hℓL

namespace WilsonSuperAux

omit [DecidableEq ι] in
/-- Summing the field profiles of the sectors over the cuts `ℓ ∈ [1, L]` is summing `fieldSqS` over
the links. -/
theorem sum_cutField_sectors (θ : ℝ) (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ : ι → Cfg L → ℂ) :
    ∑ ℓ ∈ Icc 1 L, ∑ a, WilsonAux.cutField θ e₀ (Q a) (Ψ a) ℓ = ∑ j, fieldSqS θ e₀ Q Ψ j := by
  unfold fieldSqS
  rw [Finset.sum_comm, Finset.sum_comm (s := (univ : Finset (Fin L)))]
  exact Finset.sum_congr rfl (fun a _ => WilsonAux.sum_cutField θ e₀ (Q a) (Ψ a))

/-- Per-cut bound for the time-dependent part:
`S_N(Ψ(t); [1,ℓ]) - H_{Ψ₀}(UQ_{[1,ℓ]}) ≤ ½ ln(2πe(⟨(E_ℓ + ϑ)²⟩_t + 1/12))`. -/
theorem numberEntropy_sub_le_cut (w m g θ : ℝ) (e₀ : ℤ) (Q : ι → Fin L → ℤ)
    (Ψ₀ : ι → Cfg L → ℂ) (hΨ₀ : ∑ a, ∑ n, Complex.normSq (Ψ₀ a n) = 1) (t : ℝ) (ℓ : ℕ)
    (hℓ : ℓ ∈ Icc 1 L) :
    shannon L (numDistS (segment ℓ) (wStateS w m g θ e₀ Q Ψ₀ t)) - bgEntropyS Q Ψ₀ ℓ ≤
      1 / 2 * Real.log (2 * π * exp 1 *
        (∑ a, WilsonAux.cutField θ e₀ (Q a) (wStateS w m g θ e₀ Q Ψ₀ t a) ℓ + 1 / 12)) := by
  rw [Finset.mem_Icc] at hℓ
  have h1 := wilson_numberEntropy_le_super w m g θ e₀ Q Ψ₀ hΨ₀ t ℓ (by omega) hℓ.2
  have h2 := lawEntropy_le_massey
    (fun c : ι × Cfg L => wField e₀ (Q c.1) c.2 ⟨ℓ - 1, by omega⟩)
    (fun c => wStateS w m g θ e₀ Q Ψ₀ t c.1 c.2)
    (by rw [Fintype.sum_prod_type]; exact sum_normSq_wStateS_eq_one w m g θ e₀ Q Ψ₀ hΨ₀ t) θ
  have h3 : ∑ c : ι × Cfg L, Complex.normSq (wStateS w m g θ e₀ Q Ψ₀ t c.1 c.2) *
      ((wField e₀ (Q c.1) c.2 ⟨ℓ - 1, by omega⟩ : ℝ) + θ) ^ 2 =
      ∑ a, WilsonAux.cutField θ e₀ (Q a) (wStateS w m g θ e₀ Q Ψ₀ t a) ℓ := by
    rw [Fintype.sum_prod_type]
    refine Finset.sum_congr rfl (fun a _ => ?_)
    unfold WilsonAux.cutField fieldSq
    rw [dite_eq_left ⟨by omega, hℓ.2⟩]
  rw [h3] at h2
  have h4 : fieldEntropyS e₀ Q (wStateS w m g θ e₀ Q Ψ₀ t) ⟨ℓ - 1, by omega⟩ =
      lawEntropy (fun c : ι × Cfg L => wField e₀ (Q c.1) c.2 ⟨ℓ - 1, by omega⟩)
        (fun c => wStateS w m g θ e₀ Q Ψ₀ t c.1 c.2) := rfl
  linarith

end WilsonSuperAux

/-- **The time-dependent part of `S_N` is bounded on average over the cuts** (paper: Remark
`rem:wilson`(a)): `(1/L) ∑_ℓ (S_N(Ψ(t); [1,ℓ]) - H_{Ψ₀}(UQ_{[1,ℓ]})) ≤ ½ ln(2πe((2/g²)(ε₀ + w + |m|) + 1/12))`,
uniformly in `t` and `L`. -/
theorem wilson_numberEntropy_avg_super (w m g θ : ℝ) (hw : 0 ≤ w) (hg : 0 < g) (e₀ : ℤ)
    (Q : ι → Fin L → ℤ) (hL : 0 < L) (Ψ₀ : ι → Cfg L → ℂ)
    (hΨ₀ : ∑ a, ∑ n, Complex.normSq (Ψ₀ a n) = 1) (t : ℝ) :
    (1 / (L : ℝ)) * ∑ ℓ ∈ Icc 1 L,
        (shannon L (numDistS (segment ℓ) (wStateS w m g θ e₀ Q Ψ₀ t)) - bgEntropyS Q Ψ₀ ℓ) ≤
      1 / 2 * Real.log (2 * π * exp 1 *
        (2 / g ^ 2 * (energyDensityS w m g θ e₀ Q Ψ₀ + w + |m|) + 1 / 12)) := by
  have hLpos : (0 : ℝ) < L := Nat.cast_pos.2 hL
  have hbound := (wilson_field_bound_super w m g θ hw hg e₀ Q hL Ψ₀ hΨ₀ t).2
  rw [← WilsonSuperAux.sum_cutField_sectors] at hbound
  set G : ℕ → ℝ := fun ℓ => ∑ a, WilsonAux.cutField θ e₀ (Q a) (wStateS w m g θ e₀ Q Ψ₀ t a) ℓ
    with hG
  have hG0 : ∀ ℓ, 0 ≤ G ℓ := fun ℓ =>
    Finset.sum_nonneg (fun a _ => WilsonAux.cutField_nonneg θ e₀ (Q a) _ ℓ)
  have hcard : ((Icc 1 L).card : ℝ) = L := by simp
  have hw1 : ∑ _ℓ ∈ Icc 1 L, (1 / (L : ℝ)) = 1 := by
    rw [Finset.sum_const, nsmul_eq_mul, hcard]
    field_simp
  have hJ : ∑ ℓ ∈ Icc 1 L, (1 / (L : ℝ)) • Real.log (2 * π * exp 1 * (G ℓ + 1 / 12)) ≤
      Real.log (∑ ℓ ∈ Icc 1 L, (1 / (L : ℝ)) • (2 * π * exp 1 * (G ℓ + 1 / 12))) :=
    strictConcaveOn_log_Ioi.concaveOn.le_map_sum (fun _ _ => by positivity) hw1
      (fun ℓ _ => by
        have := hG0 ℓ
        rw [Set.mem_Ioi]
        positivity)
  have havg : ∑ ℓ ∈ Icc 1 L, (1 / (L : ℝ)) • (2 * π * exp 1 * (G ℓ + 1 / 12)) =
      2 * π * exp 1 * ((1 / (L : ℝ)) * ∑ ℓ ∈ Icc 1 L, G ℓ + 1 / 12) := by
    have hL0 : (L : ℝ) ≠ 0 := hLpos.ne'
    simp only [smul_eq_mul]
    rw [← Finset.mul_sum, ← Finset.mul_sum, Finset.sum_add_distrib, Finset.sum_const,
      nsmul_eq_mul, hcard]
    field_simp
  have hA0 : 0 ≤ (1 / (L : ℝ)) * ∑ ℓ ∈ Icc 1 L, G ℓ :=
    mul_nonneg (by positivity) (Finset.sum_nonneg (fun ℓ _ => hG0 ℓ))
  calc (1 / (L : ℝ)) * ∑ ℓ ∈ Icc 1 L,
        (shannon L (numDistS (segment ℓ) (wStateS w m g θ e₀ Q Ψ₀ t)) - bgEntropyS Q Ψ₀ ℓ)
      ≤ (1 / (L : ℝ)) * ∑ ℓ ∈ Icc 1 L, 1 / 2 * Real.log (2 * π * exp 1 * (G ℓ + 1 / 12)) := by
        gcongr with ℓ hℓ
        exact WilsonSuperAux.numberEntropy_sub_le_cut w m g θ e₀ Q Ψ₀ hΨ₀ t ℓ hℓ
    _ = 1 / 2 * ∑ ℓ ∈ Icc 1 L, (1 / (L : ℝ)) • Real.log (2 * π * exp 1 * (G ℓ + 1 / 12)) := by
        simp only [smul_eq_mul]
        rw [Finset.mul_sum, Finset.mul_sum]
        exact Finset.sum_congr rfl (fun ℓ _ => by ring)
    _ ≤ 1 / 2 * Real.log (∑ ℓ ∈ Icc 1 L, (1 / (L : ℝ)) • (2 * π * exp 1 * (G ℓ + 1 / 12))) := by
        gcongr
    _ = 1 / 2 * Real.log (2 * π * exp 1 * ((1 / (L : ℝ)) * ∑ ℓ ∈ Icc 1 L, G ℓ + 1 / 12)) := by
        rw [havg]
    _ ≤ _ := WilsonAux.phi_mono hA0 hbound

/-! ### Bare-vacuum initial states -/

set_option linter.unusedSectionVars false in
/-- For bare-vacuum initial states (`ρ_j Ψ₀ = 0`), `UQ_{[1,ℓ]} = E_ℓ - e₀` on `Ψ₀`, so the background
term is the entropy of the initial field on the cut link. -/
theorem bgEntropyS_bare_vacuum (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ₀ : ι → Cfg L → ℂ)
    (hvac : ∀ a n, Ψ₀ a n ≠ 0 → ∀ j, stagRho n j = 0) (ℓ : ℕ) (hℓ : 0 < ℓ) (hℓL : ℓ ≤ L) :
    bgEntropyS Q Ψ₀ ℓ = fieldEntropyS e₀ Q Ψ₀ ⟨ℓ - 1, by omega⟩ := by
  have hset : univ.filter (fun k : Fin L => k ≤ ⟨ℓ - 1, by omega⟩) =
      univ.filter (fun k : Fin L => k.val < ℓ) := by
    ext k
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, Fin.le_def]
    omega
  refine (WilsonSuperAux.lawEntropy_congr_inj _ _ (fun y => e₀ + y)
    (fun a b h => by simpa using h) _ (fun c hc => ?_)).symm
  have h0 : ∀ k, stagRho c.2 k = 0 := hvac c.1 c.2 hc
  show e₀ + ∑ k ∈ univ.filter (fun k => k ≤ (⟨ℓ - 1, by omega⟩ : Fin L)),
      (stagRho c.2 k + Q c.1 k) = e₀ + bgSum (Q c.1) ℓ
  rw [hset]
  simp only [h0, zero_add, bgSum]

/-- For bare-vacuum initial states whose field on the cut link takes values in `{-1, 0, 1}`, the
background term is at most `ln 3`. -/
theorem bgEntropyS_le_log_three (e₀ : ℤ) (Q : ι → Fin L → ℤ) (Ψ₀ : ι → Cfg L → ℂ)
    (hΨ₀ : ∑ a, ∑ n, Complex.normSq (Ψ₀ a n) = 1)
    (hvac : ∀ a n, Ψ₀ a n ≠ 0 → ∀ j, stagRho n j = 0) (ℓ : ℕ) (hℓ : 0 < ℓ) (hℓL : ℓ ≤ L)
    (h3 : ∀ a n, Ψ₀ a n ≠ 0 → wField e₀ (Q a) n ⟨ℓ - 1, by omega⟩ ∈ ({-1, 0, 1} : Finset ℤ)) :
    bgEntropyS Q Ψ₀ ℓ ≤ Real.log 3 := by
  rw [bgEntropyS_bare_vacuum e₀ Q Ψ₀ hvac ℓ hℓ hℓL]
  have h := lawEntropy_le_log_card
    (fun c : ι × Cfg L => wField e₀ (Q c.1) c.2 ⟨ℓ - 1, by omega⟩) (fun c => Ψ₀ c.1 c.2)
    (by rw [Fintype.sum_prod_type]; exact hΨ₀) ({-1, 0, 1} : Finset ℤ)
    (fun c hc => h3 c.1 c.2 hc)
  have hcard : (({-1, 0, 1} : Finset ℤ).card : ℝ) = 3 := by
    rw [show ({-1, 0, 1} : Finset ℤ).card = 3 by decide]
    norm_num
  rw [hcard] at h
  exact h

end DFL
