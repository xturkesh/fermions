import DFL.Transfer

/-!
# The hypotheses of Furstenberg's theorem for the Bernoulli–Anderson model
(paper: Proposition `prop:nocritical`(i),(ii),(iv))

For `λ ≠ 0` and every energy `E`, the two transfer matrices `T_E(±λ)` (the support of the law of
`T_E(λ q₀)`) generate a non-compact, strongly irreducible subgroup of `SL(2, ℝ)`, and they do not
commute (no critical energies in the sense of Jitomirskaya–Schulz-Baldes–Stolz). With Furstenberg's
theorem (not formalized) this gives the positivity of the Lyapunov exponent, which is part of the
literature input `DFL.LDT`.

* `transferT_mul_inv`: `P := T_E(λ) T_E(-λ)⁻¹ = [[1, -2λ], [0, 1]]`.
* `transferT_commutator`: `T_E(v) T_E(v') - T_E(v') T_E(v) = (v - v') [[0, 1], [1, 0]]`.
* `unipotent_pow_unbounded`: the powers of `P` are unbounded (non-compactness, (i)).
* `no_finite_invariant_lines`: no finite non-empty set of lines in `ℝ²` is invariant under both
  `T_E(λ)` and `T_E(-λ)` (strong irreducibility, (ii); invariance under the generators is what has to
  be excluded).
-/

open Matrix

noncomputable section

namespace DFL

namespace FurstenbergAux

/-- Powers of the unipotent matrix: `[[1, -2λ], [0, 1]]^k = [[1, -2kλ], [0, 1]]`. -/
lemma unipotent_pow (lam : ℝ) (k : ℕ) :
    (!![1, -2 * lam; 0, 1] : Matrix (Fin 2) (Fin 2) ℝ) ^ k = !![1, -2 * k * lam; 0, 1] := by
  induction k with
  | zero =>
    ext i j
    fin_cases i <;> fin_cases j <;> simp
  | succ k ih =>
    rw [pow_succ, ih, Matrix.mul_fin_two]
    ext i j
    fin_cases i <;> fin_cases j <;> simp
    ring

/-- First coordinate of `T_E(v) u`. -/
lemma transferT_mulVec_zero (E v : ℝ) (u : Fin 2 → ℝ) :
    (transferT E v *ᵥ u) 0 = (E - v) * u 0 - u 1 := by
  simp [transferT, Matrix.mulVec, dotProduct, Fin.sum_univ_two]
  ring

/-- Second coordinate of `T_E(v) u`. -/
lemma transferT_mulVec_one (E v : ℝ) (u : Fin 2 → ℝ) : (transferT E v *ᵥ u) 1 = u 0 := by
  simp [transferT, Matrix.mulVec, dotProduct, Fin.sum_univ_two]

/-- The slope `u₀ / u₁` of the line spanned by `u` (`none` for the line `u₁ = 0`). -/
def lineSlope (u : Fin 2 → ℝ) : Option ℝ := if u 1 = 0 then none else some (u 0 / u 1)

/-- The action of `T_E(v)` on slopes: `s ↦ E - v - 1/s`, with `0 ↦ ∞ ↦ E - v`. -/
def mob (E v : ℝ) : Option ℝ → Option ℝ
  | none => some (E - v)
  | some s => if s = 0 then none else some (E - v - s⁻¹)

/-- `T_E(v)` acts on the slopes of non-zero vectors by `mob E v`. -/
lemma lineSlope_transferT (E v : ℝ) (u : Fin 2 → ℝ) (hu : u ≠ 0) :
    lineSlope (transferT E v *ᵥ u) = mob E v (lineSlope u) := by
  unfold lineSlope
  rw [transferT_mulVec_zero, transferT_mulVec_one]
  by_cases hb : u 1 = 0
  · have ha : u 0 ≠ 0 := by
      intro ha
      apply hu
      ext i
      fin_cases i <;> simp [ha, hb]
    simp only [hb, ha, ite_true, ite_false, mob, sub_zero, Option.some.injEq]
    field_simp
  · by_cases ha : u 0 = 0
    · simp [ha, hb, mob]
    · have hab : u 0 / u 1 ≠ 0 := div_ne_zero ha hb
      simp only [hb, ha, hab, ite_false, mob, inv_div, Option.some.injEq]
      field_simp

/-- Non-zero multiples have the same slope. -/
lemma lineSlope_smul (c : ℝ) (hc : c ≠ 0) (u : Fin 2 → ℝ) : lineSlope (c • u) = lineSlope u := by
  unfold lineSlope
  simp only [Pi.smul_apply, smul_eq_mul, mul_eq_zero, hc, false_or]
  split_ifs
  · rfl
  · rw [mul_div_mul_left _ _ hc]

/-- `mob E v` is injective (`T_E(v)` is invertible). -/
lemma mob_injective (E v : ℝ) : Function.Injective (mob E v) := by
  intro s t h
  rcases s with _ | s <;> rcases t with _ | t
  · rfl
  · by_cases ht : t = 0
    · simp [mob, ht] at h
    · simp only [mob, ht, ite_false, Option.some.injEq] at h
      exact absurd (inv_eq_zero.1 (by linarith)) ht
  · by_cases hs : s = 0 <;> simp [mob, hs] at h
  · by_cases hs : s = 0 <;> by_cases ht : t = 0 <;> simp [mob, hs, ht] at h ⊢
    exact h

/-- `T_E(λ) = P T_E(-λ)` on slopes, where `P` is the translation by `-2λ`. -/
lemma mob_eq_map (E lam : ℝ) (s : Option ℝ) :
    mob E lam s = (mob E (-lam) s).map (fun t => t - 2 * lam) := by
  rcases s with _ | s
  · simp [mob]
    ring
  · by_cases hs : s = 0
    · simp [mob, hs]
    · simp [mob, hs]
      ring

end FurstenbergAux

/-- `T_E(λ) T_E(-λ)⁻¹ = [[1, -2λ], [0, 1]]` (with the inverse written as the adjugate). -/
theorem transferT_mul_inv (E lam : ℝ) :
    transferT E lam * (transferT E (-lam)).adjugate = !![1, -2 * lam; 0, 1] := by
  rw [transferT, transferT, Matrix.adjugate_fin_two_of, Matrix.mul_fin_two]
  ext i j
  fin_cases i <;> fin_cases j <;> simp
  ring

/-- The adjugate is the inverse: `T_E(v) adj(T_E(v)) = 1`. -/
theorem transferT_mul_adjugate (E v : ℝ) : transferT E v * (transferT E v).adjugate = 1 := by
  rw [transferT, Matrix.adjugate_fin_two_of, Matrix.mul_fin_two]
  ext i j
  fin_cases i <;> fin_cases j <;> simp

/-- **No critical energies** (paper: Proposition `prop:nocritical`(iv)):
`T_E(v) T_E(v') - T_E(v') T_E(v) = (v - v') [[0, 1], [1, 0]]`; for `v ≠ v'` the transfer matrices do
not commute. -/
theorem transferT_commutator (E v v' : ℝ) :
    transferT E v * transferT E v' - transferT E v' * transferT E v = (v - v') • !![0, 1; 1, 0] := by
  rw [transferT, transferT, Matrix.mul_fin_two, Matrix.mul_fin_two]
  ext i j
  fin_cases i <;> fin_cases j <;> simp
  ring

/-- The transfer matrices `T_E(λ)`, `T_E(-λ)` do not commute for `λ ≠ 0`. -/
theorem transferT_not_commute (E lam : ℝ) (hlam : lam ≠ 0) :
    transferT E lam * transferT E (-lam) ≠ transferT E (-lam) * transferT E lam := by
  intro h
  have h1 := congrFun (congrFun (transferT_commutator E lam (-lam)) 0) 1
  rw [h, sub_self] at h1
  have h2 : (0 : ℝ) = 2 * lam := by simpa [two_mul] using h1
  exact hlam (by linarith)

/-- **Non-compactness** (paper: Proposition `prop:nocritical`(i)): the powers of
`P = [[1, -2λ], [0, 1]]` have unbounded entries. -/
theorem unipotent_pow_unbounded (lam : ℝ) (hlam : lam ≠ 0) (R : ℝ) :
    ∃ k : ℕ, R < |((!![1, -2 * lam; 0, 1] : Matrix (Fin 2) (Fin 2) ℝ) ^ k) 0 1| := by
  obtain ⟨k, hk⟩ := exists_nat_gt (R / (2 * |lam|))
  refine ⟨k, ?_⟩
  have h01 : ((!![1, -2 * lam; 0, 1] : Matrix (Fin 2) (Fin 2) ℝ) ^ k) 0 1 = -2 * k * lam := by
    rw [FurstenbergAux.unipotent_pow]
    simp
  rw [h01]
  have hl : 0 < 2 * |lam| := by positivity
  rw [div_lt_iff₀ hl] at hk
  calc R < k * (2 * |lam|) := hk
    _ = |-2 * k * lam| := by
      rw [abs_mul, abs_mul, abs_neg, abs_two, Nat.abs_cast]; ring

/-- **Strong irreducibility** (paper: Proposition `prop:nocritical`(ii)): a finite set of non-zero
vectors whose spanned lines are permuted by `T_E(λ)` and by `T_E(-λ)` is empty. Lines are encoded by
non-zero spanning vectors; "the image of the line of `u` is a line of the set" means that `T u` is a
non-zero multiple of some element. -/
theorem no_finite_invariant_lines (E lam : ℝ) (hlam : lam ≠ 0) (L : Finset (Fin 2 → ℝ))
    (hL0 : ∀ u ∈ L, u ≠ 0)
    (hinv : ∀ u ∈ L, ∀ v : ℝ, (v = lam ∨ v = -lam) →
      ∃ u' ∈ L, ∃ c : ℝ, c ≠ 0 ∧ transferT E v *ᵥ u = c • u') :
    L = ∅ := by
  classical
  -- the finite set of slopes of the lines of `L`
  set S : Finset (Option ℝ) := L.image FurstenbergAux.lineSlope with hSdef
  have hS : ∀ s ∈ S, ∀ v : ℝ, (v = lam ∨ v = -lam) → FurstenbergAux.mob E v s ∈ S := by
    intro s hs v hv
    obtain ⟨u, hu, rfl⟩ := Finset.mem_image.1 hs
    obtain ⟨u', hu', c, hc, h⟩ := hinv u hu v hv
    rw [← FurstenbergAux.lineSlope_transferT E v u (hL0 u hu), h, FurstenbergAux.lineSlope_smul c hc]
    exact Finset.mem_image_of_mem _ hu'
  -- `T_E(-λ)` permutes the finite set `S` (injective self-map of a finite set)
  have hsurj : S.image (FurstenbergAux.mob E (-lam)) = S := by
    apply Finset.eq_of_subset_of_card_le
    · intro s hs
      obtain ⟨s', hs', rfl⟩ := Finset.mem_image.1 hs
      exact hS s' hs' _ (Or.inr rfl)
    · rw [Finset.card_image_of_injective _ (FurstenbergAux.mob_injective E (-lam))]
  -- hence `S` is invariant under the translation by `-2λ` (the action of `P`)
  have hshift : ∀ t : ℝ, some t ∈ S → some (t - 2 * lam) ∈ S := by
    intro t ht
    rw [← hsurj] at ht
    obtain ⟨s', hs', hs'eq⟩ := Finset.mem_image.1 ht
    have h := hS s' hs' lam (Or.inl rfl)
    rw [FurstenbergAux.mob_eq_map, hs'eq] at h
    simpa using h
  -- no finite slopes: their translates `t - 2kλ` would be infinitely many
  have hfin : ∀ t : ℝ, some t ∉ S := by
    intro t ht
    have hk : ∀ k : ℕ, some (t - 2 * k * lam) ∈ (S : Set (Option ℝ)) := by
      intro k
      induction k with
      | zero => simpa using ht
      | succ k ih =>
        have h := hshift _ ih
        rw [Finset.mem_coe]
        convert h using 2
        push_cast
        ring
    have hinj : Function.Injective (fun k : ℕ => some (t - 2 * k * lam)) := by
      intro k k' h
      simp only [Option.some.injEq, sub_right_inj, mul_eq_mul_right_iff, hlam, or_false,
        mul_eq_mul_left_iff, two_ne_zero, Nat.cast_inj] at h
      exact h
    exact Set.infinite_of_injective_forall_mem hinj hk S.finite_toSet
  -- no vertical line either: `T_E(λ)` maps it to the slope `E - λ`
  have hnone : (none : Option ℝ) ∉ S := fun h => hfin _ (hS none h lam (Or.inl rfl))
  have hSe : S = ∅ := by
    rw [Finset.eq_empty_iff_forall_notMem]
    rintro (_ | t)
    · exact hnone
    · exact hfin t
  rwa [hSdef, Finset.image_eq_empty] at hSe

end DFL
