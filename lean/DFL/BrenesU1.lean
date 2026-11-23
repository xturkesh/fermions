import DFL.WilsonSuper

/-!
# The initial states of Brenes, Dalmonte, Heyl and Scardicchio (paper: Proposition `prop:brenes`)

Lattice Schwinger model with Wilson links (`schwinger`), initial state: bare vacuum of the staggered
fermions (`bareVac`, all `ρ_j = 0`) and every link in the equal superposition of `E ∈ {-1, 0, 1}`, with the
boundary field `E_0 = e₀` fixed. The state is the uniform superposition over the background field configurations
`ē ∈ {-1, 0, 1}^L` (`Fin L → Fin 3`, value `fieldOf`) of `|vac⟩ ⊗ |ē⟩`, and `|vac⟩ ⊗ |ē⟩` lies in the sector of
background charges `Q_j = ē_j - ē_{j-1}` (`brenesQ`): the background *fields* are i.i.d. uniform, the
background charges are not independent.

* `wField_brenes`: in the sector of `ē`, `E_j = ē_j + F_j`, `F_j = ∑_{k ≤ j} ρ_k`;
* `schwinger_brenes`: the sector Hamiltonian is the clean Schwinger model (`Q = 0`, `e₀ = 0`, `ϑ = 0`) plus the
  potential `g² A_k` on the charges, `A_k = ∑_{j ≥ k} (ē_j + ϑ)` (`rwPot`, a random walk with i.i.d. steps,
  `rwPot_succ`), plus a constant;
* `schwinger_meson`: a tightly bound meson on the bond `(j, j+1)` costs `2m + (g²/2)(1 + 2 s_j (ē_j + ϑ))`,
  an i.i.d. potential;
* `bgSum_brenes`, `bgEntropyS_brenes`: the background term of Remark `rem:wilson`(a) is exactly `ln 3` for every
  cut; `energyDensityS_brenes`: `ε₀ = (g²/2)(2/3 + ϑ²) - m ⌈L/2⌉/L`;
* `brenes_numberEntropy_avg`: the cut-averaged number entropy is at most
  `ln 3 + ½ ln(2πe((2/g²)(ε₀ + w + |m|) + 1/12))`, uniformly in time and `L`;
* `brenes_field_memory`: for all times, the electric energy density stays within `4(w + |m|)/g²` of its initial
  value `2/3 + ϑ²`.
-/

open Matrix Finset Real

noncomputable section

namespace DFL

variable {L : ℕ}

/-- The value `ē ∈ {-1, 0, 1}` of a background field label `a : Fin 3`. -/
def fieldOf (a : Fin 3) : ℤ := (a : ℤ) - 1

/-- The background charges of the field configuration `ē`: `Q_j = ē_j - ē_{j-1}`, with `ē_{-1} = e₀` (Lean
indexing; paper: `Q_j = ē_j - ē_{j-1}`, `ē_0 = e₀`). -/
def brenesQ (e₀ : ℤ) (eb : Fin L → Fin 3) (j : Fin L) : ℤ :=
  fieldOf (eb j) - (if hj : j.val = 0 then e₀ else fieldOf (eb ⟨j.val - 1, by omega⟩))

/-- The bare vacuum of the staggered fermions: occupied exactly on the paper's odd sites (Lean: even `j`). -/
def bareVac : Cfg L := fun j => if j.val % 2 = 0 then 1 else 0

/-- The initial state of Brenes et al.: its component in the sector of `ē` is `3^{-L/2} |vac⟩`. -/
def brenesInit (_eb : Fin L → Fin 3) (n : Cfg L) : ℂ :=
  if n = bareVac then (((Real.sqrt 3) ^ L : ℝ) : ℂ)⁻¹ else 0

/-- The random-walk potential `A_k = ∑_{j ≥ k} (ē_j + ϑ)` felt by the charges. -/
def rwPot (θ : ℝ) (eb : Fin L → Fin 3) (k : Fin L) : ℝ :=
  ∑ j ∈ univ.filter (fun j => k ≤ j), ((fieldOf (eb j) : ℝ) + θ)

namespace BrenesU1Aux

/-- `|brenesInit eb n|² = 3^{-L} [n = bareVac]`. -/
theorem normSq_brenesInit (eb : Fin L → Fin 3) (n : Cfg L) :
    Complex.normSq (brenesInit eb n) = if n = bareVac then ((3 : ℝ) ^ L)⁻¹ else 0 := by
  unfold brenesInit
  split_ifs with h
  · rw [map_inv₀, Complex.normSq_ofReal, ← mul_pow,
      Real.mul_self_sqrt (by norm_num : (0 : ℝ) ≤ 3)]
  · simp

/-- The number of background field configurations is `3^L`. -/
theorem card_eb : Fintype.card (Fin L → Fin 3) = 3 ^ L := by
  simp

/-- Each sector carries the weight `3^{-L}`. -/
theorem sum_normSq_brenesInit (eb : Fin L → Fin 3) :
    ∑ n, Complex.normSq (brenesInit eb n) = ((3 : ℝ) ^ L)⁻¹ := by
  simp only [normSq_brenesInit]
  simp

/-- The background fields in paper indexing: `ē_0 = e₀` and `ē_{k+1} = fieldOf (eb k)`. -/
def ebarN (e₀ : ℤ) (eb : Fin L → Fin 3) (k : ℕ) : ℤ :=
  if k = 0 then e₀ else if h : k - 1 < L then fieldOf (eb ⟨k - 1, h⟩) else 0

/-- `ē_{k+1} = fieldOf (eb k)`. -/
theorem ebarN_succ (e₀ : ℤ) (eb : Fin L → Fin 3) (k : Fin L) :
    ebarN e₀ eb (k.val + 1) = fieldOf (eb k) := by
  unfold ebarN
  rw [ite_eq_right (Nat.succ_ne_zero _), dite_eq_left (by simp)]
  simp

/-- `ē_ℓ = fieldOf (eb (ℓ - 1))` for `0 < ℓ ≤ L`. -/
theorem ebarN_of_pos (e₀ : ℤ) (eb : Fin L → Fin 3) (ℓ : ℕ) (hℓ : 0 < ℓ) (hℓL : ℓ ≤ L) :
    ebarN e₀ eb ℓ = fieldOf (eb ⟨ℓ - 1, by omega⟩) := by
  unfold ebarN
  rw [ite_eq_right (by omega), dite_eq_left (by omega)]

/-- `Q_k = ē_{k+1} - ē_k` (paper indexing for `ē`). -/
theorem brenesQ_eq (e₀ : ℤ) (eb : Fin L → Fin 3) (k : Fin L) :
    brenesQ e₀ eb k = ebarN e₀ eb (k.val + 1) - ebarN e₀ eb k.val := by
  rw [ebarN_succ]
  unfold brenesQ ebarN
  have hk := k.isLt
  split_ifs with h1 h2
  · rfl
  · rfl
  · omega

/-- Sums over the sites `k < ℓ` of `Fin L` are sums over `range ℓ`. -/
theorem sum_filter_lt {M : Type*} [AddCommMonoid M] (f : ℕ → M) (ℓ : ℕ) (hℓ : ℓ ≤ L) :
    ∑ k ∈ univ.filter (fun k : Fin L => k.val < ℓ), f k.val = ∑ i ∈ range ℓ, f i := by
  rw [Finset.sum_filter]
  rw [Fin.sum_univ_eq_sum_range (fun i => if i < ℓ then f i else 0) L, ← Finset.sum_filter]
  congr 1
  ext i
  simp only [Finset.mem_filter, Finset.mem_range]
  omega

/-- **Telescoping**: `∑_{k < ℓ} Q_k = ē_ℓ - e₀`. -/
theorem sum_brenesQ (e₀ : ℤ) (eb : Fin L → Fin 3) (ℓ : ℕ) (hℓ : ℓ ≤ L) :
    ∑ k ∈ univ.filter (fun k : Fin L => k.val < ℓ), brenesQ e₀ eb k = ebarN e₀ eb ℓ - e₀ := by
  simp only [brenesQ_eq]
  rw [sum_filter_lt (fun i => ebarN e₀ eb (i + 1) - ebarN e₀ eb i) ℓ hℓ, Finset.sum_range_sub]
  simp [ebarN]

end BrenesU1Aux

/-- The Brenes state is normalized. -/
theorem brenesInit_norm : ∑ eb, ∑ n, Complex.normSq (brenesInit (L := L) eb n) = 1 := by
  simp only [BrenesU1Aux.sum_normSq_brenesInit, Finset.sum_const, Finset.card_univ,
    BrenesU1Aux.card_eb, nsmul_eq_mul]
  push_cast
  exact mul_inv_cancel₀ (by positivity)

/-- The bare vacuum carries no charge. -/
theorem stagRho_bareVac (j : Fin L) : stagRho (bareVac (L := L)) j = 0 := by
  unfold stagRho bareVac
  by_cases h : j.val % 2 = 0 <;> simp [h]

/-- **Fields in the sector of `ē`**: `E_j = ē_j + ∑_{k ≤ j} ρ_k`. -/
theorem wField_brenes (e₀ : ℤ) (eb : Fin L → Fin 3) (n : Cfg L) (j : Fin L) :
    wField e₀ (brenesQ e₀ eb) n j =
      fieldOf (eb j) + ∑ k ∈ univ.filter (fun k => k ≤ j), stagRho n k := by
  have hset : univ.filter (fun k : Fin L => k ≤ j) =
      univ.filter (fun k : Fin L => k.val < j.val + 1) := by
    ext k
    simp only [mem_filter, mem_univ, true_and, Fin.le_def]
    omega
  unfold wField
  rw [Finset.sum_add_distrib, hset,
    BrenesU1Aux.sum_brenesQ e₀ eb (j.val + 1) (by have := j.isLt; omega),
    BrenesU1Aux.ebarN_succ]
  ring

/-- The steps of the random walk are the background fields: `A_k - A_{k+1} = ē_k + ϑ`. -/
theorem rwPot_succ (θ : ℝ) (eb : Fin L → Fin 3) (k : Fin L) (hk : k.val + 1 < L) :
    rwPot θ eb k - rwPot θ eb ⟨k.val + 1, hk⟩ = (fieldOf (eb k) : ℝ) + θ := by
  have hset : univ.filter (fun j : Fin L => k ≤ j) =
      insert k (univ.filter (fun j : Fin L => (⟨k.val + 1, hk⟩ : Fin L) ≤ j)) := by
    ext j
    simp only [mem_filter, mem_univ, true_and, mem_insert, Fin.le_def, Fin.ext_iff]
    omega
  have hnot : k ∉ univ.filter (fun j : Fin L => (⟨k.val + 1, hk⟩ : Fin L) ≤ j) := by
    simp [Fin.le_def]
  unfold rwPot
  rw [hset, Finset.sum_insert hnot]
  ring

namespace BrenesU1Aux

/-- In the clean sector (`e₀ = 0`, `Q = 0`) the field is the cumulative charge `F_j = ∑_{k ≤ j} ρ_k`. -/
theorem wField_clean (n : Cfg L) (j : Fin L) :
    (wField 0 (fun _ => 0) n j : ℝ) = ∑ k ∈ univ.filter (fun k => k ≤ j), (stagRho n k : ℝ) := by
  unfold wField
  push_cast
  simp

/-- In the sector of `ē`, `E_j = ē_j + F_j` with `F_j` the field of the clean sector. -/
theorem wField_brenes_clean (e₀ : ℤ) (eb : Fin L → Fin 3) (n : Cfg L) (j : Fin L) :
    (wField e₀ (brenesQ e₀ eb) n j : ℝ) = (fieldOf (eb j) : ℝ) + (wField 0 (fun _ => 0) n j : ℝ) := by
  rw [wField_brenes, wField_clean]
  push_cast
  rfl

/-- Exchanging the sums: `∑_j a_j F_j = ∑_k A_k ρ_k` with `A_k = ∑_{j ≥ k} a_j`. -/
theorem sum_mul_wField_clean (θ : ℝ) (eb : Fin L → Fin 3) (n : Cfg L) :
    ∑ j, ((fieldOf (eb j) : ℝ) + θ) * (wField 0 (fun _ => 0) n j : ℝ) =
      ∑ k, rwPot θ eb k * (stagRho n k : ℝ) := by
  simp only [wField_clean]
  unfold rwPot
  simp only [Finset.mul_sum, Finset.sum_mul, Finset.sum_filter]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun k _ => Finset.sum_congr rfl (fun j _ => ?_))
  split_ifs <;> simp

/-- The hopping has no diagonal entries. -/
theorem hop_diag (j : Fin L) (n : Cfg L) : hop j n n = 0 := by
  rw [hop_apply]
  split_ifs with hj h
  · exfalso
    obtain ⟨hne, heq⟩ := h
    apply hne
    have := congrFun heq j
    simpa [Equiv.swap_apply_left] using this
  · rfl
  · rfl

/-- The diagonal entries of the Schwinger Hamiltonian: mass plus field energy. -/
theorem schwinger_diag (w m g θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) (n : Cfg L) :
    schwinger w m g θ e₀ Q n n = m * ∑ j : Fin L, (-1 : ℝ) ^ (j.val + 1) * occ (n j) +
      g ^ 2 / 2 * ∑ j, ((wField e₀ Q n j : ℝ) + θ) ^ 2 := by
  unfold schwinger
  rw [Matrix.add_apply, Matrix.smul_apply, Matrix.sum_apply, diagonal_apply_eq]
  simp [hop_diag]

/-- `stagRho` at a site only depends on the occupation of that site. -/
theorem stagRho_congr (n n' : Cfg L) (k : Fin L) (h : n k = n' k) :
    stagRho n k = stagRho n' k := by
  unfold stagRho
  rw [h]

/-- `(-1)^n` by parity. -/
theorem neg_one_pow_eq_ite (n : ℕ) : (-1 : ℝ) ^ n = if n % 2 = 0 then 1 else -1 := by
  rcases Nat.even_or_odd n with h | h
  · rw [h.neg_one_pow, ite_eq_left (Nat.even_iff.mp h)]
  · rw [h.neg_one_pow, ite_eq_right (by rw [Nat.odd_iff.mp h]; norm_num)]

/-- The occupations of the bare vacuum. -/
theorem occ_bareVac (k : Fin L) : occ (bareVac k) = if k.val % 2 = 0 then 1 else 0 := by
  by_cases h : k.val % 2 = 0 <;> simp [occ, bareVac, h]

end BrenesU1Aux

/-- **Clean model plus a random-walk potential** (paper: Proposition `prop:brenes`(b)): the Hamiltonian of the
sector of `ē` is the Schwinger model without background charges, plus `g² ∑_k A_k ρ_k`, plus a constant. -/
theorem schwinger_brenes (w m g θ : ℝ) (e₀ : ℤ) (eb : Fin L → Fin 3) :
    schwinger w m g θ e₀ (brenesQ e₀ eb) =
      schwinger w m g 0 0 (fun _ => 0) +
        diagonal (fun n => g ^ 2 * ∑ k, rwPot θ eb k * (stagRho n k : ℝ) +
          g ^ 2 / 2 * ∑ j, ((fieldOf (eb j) : ℝ) + θ) ^ 2) := by
  unfold schwinger
  rw [add_assoc, diagonal_add]
  congr 2
  funext n
  simp only [BrenesU1Aux.wField_brenes_clean]
  have h2 : ∀ j, ((fieldOf (eb j) : ℝ) + (wField 0 (fun _ => 0) n j : ℝ) + θ) ^ 2 =
      ((wField 0 (fun _ => 0) n j : ℝ) + 0) ^ 2 +
        2 * (((fieldOf (eb j) : ℝ) + θ) * (wField 0 (fun _ => 0) n j : ℝ)) +
          ((fieldOf (eb j) : ℝ) + θ) ^ 2 := fun j => by ring
  simp only [h2, Finset.sum_add_distrib, ← Finset.mul_sum, BrenesU1Aux.sum_mul_wField_clean]
  ring

/-- **Meson energies** (paper: Proposition `prop:brenes`(b)): creating a meson on the bond `(j, j+1)` from the
bare vacuum costs `2m + (g²/2)(1 + 2 s_j (ē_j + ϑ))` with `s_j = (-1)^{j+1}` (Lean indexing). -/
theorem schwinger_meson (w m g θ : ℝ) (e₀ : ℤ) (eb : Fin L → Fin 3) (j : Fin L) (hj : j.val + 1 < L) :
    schwinger w m g θ e₀ (brenesQ e₀ eb) (bareVac ∘ Equiv.swap j ⟨j.val + 1, hj⟩)
        (bareVac ∘ Equiv.swap j ⟨j.val + 1, hj⟩) -
      schwinger w m g θ e₀ (brenesQ e₀ eb) bareVac bareVac =
      2 * m + g ^ 2 / 2 * (1 + 2 * (-1 : ℝ) ^ (j.val + 1) * ((fieldOf (eb j) : ℝ) + θ)) := by
  have hjj' : j ≠ ⟨j.val + 1, hj⟩ := by
    intro h
    have := congrArg Fin.val h
    simp at this
  have hvac : bareVac (L := L) j ≠ bareVac ⟨j.val + 1, hj⟩ := by
    unfold bareVac
    rcases Nat.mod_two_eq_zero_or_one j.val with h | h
    · simp [h, show (j.val + 1) % 2 = 1 by omega]
    · simp [h, show (j.val + 1) % 2 = 0 by omega]
  have hhop : hop j (bareVac ∘ Equiv.swap j ⟨j.val + 1, hj⟩) bareVac ≠ 0 := by
    rw [hop_apply, dite_eq_left hj, ite_eq_left ⟨hvac, rfl⟩]
    norm_num
  -- the mass term changes by `2`
  have hmass : ∑ i : Fin L, (-1 : ℝ) ^ (i.val + 1) *
        occ ((bareVac ∘ Equiv.swap j ⟨j.val + 1, hj⟩) i) -
      ∑ i : Fin L, (-1 : ℝ) ^ (i.val + 1) * occ (bareVac i) = 2 := by
    rw [← Finset.sum_sub_distrib, Finset.sum_eq_add j ⟨j.val + 1, hj⟩ hjj']
    · simp only [Function.comp_apply, Equiv.swap_apply_left, Equiv.swap_apply_right,
        BrenesU1Aux.occ_bareVac, BrenesU1Aux.neg_one_pow_eq_ite]
      rcases Nat.mod_two_eq_zero_or_one j.val with h | h
      · simp [h, show (j.val + 1) % 2 = 1 by omega, show (j.val + 1 + 1) % 2 = 0 by omega]
        norm_num
      · simp [h, show (j.val + 1) % 2 = 0 by omega, show (j.val + 1 + 1) % 2 = 1 by omega]
        norm_num
    · intro c _ hc
      simp [Equiv.swap_apply_of_ne_of_ne hc.1 hc.2]
    · intro h
      exact absurd (Finset.mem_univ _) h
    · intro h
      exact absurd (Finset.mem_univ _) h
  -- the field changes on the link `j` only, by `s_j`
  have hρ : (stagRho (bareVac ∘ Equiv.swap j ⟨j.val + 1, hj⟩) j : ℝ) =
      (-1 : ℝ) ^ (j.val + 1) := by
    unfold stagRho bareVac
    simp only [Function.comp_apply, Equiv.swap_apply_left]
    rw [BrenesU1Aux.neg_one_pow_eq_ite]
    rcases Nat.mod_two_eq_zero_or_one j.val with h | h
    · simp [h, show (j.val + 1) % 2 = 1 by omega]
    · simp [h, show (j.val + 1) % 2 = 0 by omega]
  have hsum : ∑ k ∈ univ.filter (fun k => k ≤ j),
      stagRho (bareVac ∘ Equiv.swap j ⟨j.val + 1, hj⟩) k =
        stagRho (bareVac ∘ Equiv.swap j ⟨j.val + 1, hj⟩) j := by
    rw [Finset.sum_eq_single j]
    · intro k hk hkj
      have hkj' : k ≠ ⟨j.val + 1, hj⟩ := by
        intro h
        have h1 := (Finset.mem_filter.1 hk).2
        rw [h, Fin.le_def] at h1
        simp at h1
      rw [BrenesU1Aux.stagRho_congr _ bareVac k
        (by simp [Equiv.swap_apply_of_ne_of_ne hkj hkj']), stagRho_bareVac]
    · intro h
      simp at h
  have hsq : (-1 : ℝ) ^ (j.val + 1) * (-1 : ℝ) ^ (j.val + 1) = 1 := by
    rw [← mul_pow]
    norm_num
  have hfield : ∑ i, ((wField e₀ (brenesQ e₀ eb) (bareVac ∘ Equiv.swap j ⟨j.val + 1, hj⟩) i : ℝ) +
        θ) ^ 2 - ∑ i, ((wField e₀ (brenesQ e₀ eb) bareVac i : ℝ) + θ) ^ 2 =
      1 + 2 * (-1 : ℝ) ^ (j.val + 1) * ((fieldOf (eb j) : ℝ) + θ) := by
    rw [← Finset.sum_sub_distrib, Finset.sum_eq_single j]
    · rw [wField_brenes, wField_brenes, hsum]
      simp only [stagRho_bareVac, Finset.sum_const_zero, add_zero]
      push_cast
      rw [hρ]
      linear_combination hsq
    · intro i _ hij
      rw [wField_hop e₀ _ j i _ bareVac hij hhop]
      ring
    · intro h
      exact absurd (Finset.mem_univ _) h
  rw [BrenesU1Aux.schwinger_diag, BrenesU1Aux.schwinger_diag]
  linear_combination m * hmass + g ^ 2 / 2 * hfield

/-- The cumulative background charge of `[1, ℓ]` is the background field on the cut link: `UQ_{[1,ℓ]} = ē_ℓ - e₀`. -/
theorem bgSum_brenes (e₀ : ℤ) (eb : Fin L → Fin 3) (ℓ : ℕ) (hℓ : 0 < ℓ) (hℓL : ℓ ≤ L) :
    bgSum (brenesQ e₀ eb) ℓ = fieldOf (eb ⟨ℓ - 1, by omega⟩) - e₀ := by
  unfold bgSum
  rw [BrenesU1Aux.sum_brenesQ e₀ eb ℓ hℓL, BrenesU1Aux.ebarN_of_pos e₀ eb ℓ hℓ hℓL]

namespace BrenesU1Aux

/-- The value map `fieldOf : Fin 3 → ℤ` is injective. -/
theorem fieldOf_injective : Function.Injective fieldOf := by
  intro a b h
  unfold fieldOf at h
  apply Fin.ext
  omega

/-- `3^L = 3 · 3^{L-1}` for `L ≥ 1`. -/
theorem three_pow_eq (hL : 0 < L) : (3 : ℝ) ^ L = 3 * 3 ^ (L - 1) := by
  rw [← pow_succ']
  congr 1
  omega

/-- Summing a function of one coordinate over the background configurations: every value of the
coordinate is taken by `3^{L-1}` configurations. -/
theorem sum_eb_coord (i : Fin L) (h : Fin 3 → ℝ) :
    ∑ eb : Fin L → Fin 3, h (eb i) = (3 : ℝ) ^ (L - 1) * ∑ a, h a := by
  rw [Fintype.sum_equiv (Equiv.funSplitAt i (Fin 3)) (fun eb => h (eb i)) (fun p => h p.1)
    (fun eb => rfl), Fintype.sum_prod_type]
  simp only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
  rw [← Finset.mul_sum]
  congr 1
  simp [Fintype.card_subtype_compl]

/-- `∑_{a} (fieldOf a + ϑ)² = 2 + 3ϑ²`. -/
theorem sum_fieldOf_sq (θ : ℝ) : ∑ a : Fin 3, ((fieldOf a : ℝ) + θ) ^ 2 = 2 + 3 * θ ^ 2 := by
  have h0 : fieldOf 0 = -1 := rfl
  have h1 : fieldOf 1 = 0 := rfl
  have h2 : fieldOf 2 = 1 := rfl
  rw [Fin.sum_univ_three, h0, h1, h2]
  push_cast
  ring

/-- `∑_{ē} (ē_j + ϑ)² = 3^{L-1} (2 + 3ϑ²)`. -/
theorem sum_eb_sq (θ : ℝ) (j : Fin L) :
    ∑ eb : Fin L → Fin 3, ((fieldOf (eb j) : ℝ) + θ) ^ 2 = (3 : ℝ) ^ (L - 1) * (2 + 3 * θ ^ 2) := by
  rw [sum_eb_coord j (fun a => ((fieldOf a : ℝ) + θ) ^ 2), sum_fieldOf_sq]

/-- The law of one background field is uniform: `∑_{ē_i = a} 3^{-L} = 1/3`. -/
theorem mass_coord (i : Fin L) (a : Fin 3) :
    ∑ c ∈ univ.filter (fun c : (Fin L → Fin 3) × Cfg L => fieldOf (c.1 i) = fieldOf a),
      Complex.normSq (brenesInit c.1 c.2) = 1 / 3 := by
  have hL : 0 < L := by have := i.isLt; omega
  rw [Finset.sum_filter, Fintype.sum_prod_type]
  have h1 : ∀ eb : Fin L → Fin 3,
      (∑ n, if fieldOf (eb i) = fieldOf a then Complex.normSq (brenesInit eb n) else 0) =
        if eb i = a then ((3 : ℝ) ^ L)⁻¹ else 0 := by
    intro eb
    by_cases hP : eb i = a
    · rw [ite_eq_left hP]
      simp only [hP, ite_true]
      exact sum_normSq_brenesInit eb
    · have hne : fieldOf (eb i) ≠ fieldOf a := fun h => hP (fieldOf_injective h)
      simp [hne, hP]
  simp only [h1]
  rw [sum_eb_coord i (fun a' => if a' = a then ((3 : ℝ) ^ L)⁻¹ else 0)]
  simp only [Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  rw [three_pow_eq hL]
  field_simp

/-- The law of one background field in the Brenes state has entropy `ln 3`. -/
theorem lawEntropy_coord (i : Fin L) :
    lawEntropy (fun c : (Fin L → Fin 3) × Cfg L => fieldOf (c.1 i))
      (fun c => brenesInit c.1 c.2) = Real.log 3 := by
  unfold lawEntropy
  have himg : univ.image (fun c : (Fin L → Fin 3) × Cfg L => fieldOf (c.1 i)) =
      univ.image fieldOf := by
    ext y
    simp only [Finset.mem_image, Finset.mem_univ, true_and]
    constructor
    · rintro ⟨c, rfl⟩
      exact ⟨c.1 i, rfl⟩
    · rintro ⟨a, rfl⟩
      exact ⟨(fun _ => a, 0), rfl⟩
  rw [himg, Finset.sum_image (fieldOf_injective.injOn)]
  simp only [mass_coord]
  rw [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  simp only [Real.negMulLog, one_div, Real.log_inv]
  push_cast
  ring

/-- The sum `∑_j (-1)^{j+1} [j even]` over `j < N` is `-⌈N/2⌉`. -/
theorem sum_range_mass (N : ℕ) :
    ∑ k ∈ range N, (-1 : ℝ) ^ (k + 1) * (if k % 2 = 0 then 1 else 0) =
      -(((N + 1) / 2 : ℕ) : ℝ) := by
  induction N with
  | zero => simp
  | succ N ih =>
    rw [Finset.sum_range_succ, ih, neg_one_pow_eq_ite]
    rcases Nat.mod_two_eq_zero_or_one N with h | h
    · rw [ite_eq_left h, ite_eq_right (by omega : ¬ (N + 1) % 2 = 0),
        show (N + 1 + 1) / 2 = (N + 1) / 2 + 1 by omega]
      push_cast
      ring
    · rw [ite_eq_right (by omega : ¬ N % 2 = 0), show (N + 1 + 1) / 2 = (N + 1) / 2 by omega]
      ring

/-- The staggered mass sum of the bare vacuum: `∑_j (-1)^{j+1} n_j = -⌈L/2⌉`. -/
theorem mass_bareVac :
    ∑ j : Fin L, (-1 : ℝ) ^ (j.val + 1) * occ (bareVac j) = -(((L + 1) / 2 : ℕ) : ℝ) := by
  simp only [occ_bareVac]
  rw [Fin.sum_univ_eq_sum_range (fun k => (-1 : ℝ) ^ (k + 1) * (if k % 2 = 0 then 1 else 0)) L]
  exact sum_range_mass L

/-- The hopping has no expectation in the bare vacuum. -/
theorem re_expect_hop_brenesInit (eb : Fin L → Fin 3) :
    (expect (cplx (∑ j : Fin L, hop j)) (brenesInit eb)).re = 0 := by
  rw [WilsonAux.re_expect_cplx]
  refine Finset.sum_eq_zero (fun a _ => Finset.sum_eq_zero (fun b _ => ?_))
  by_cases ha : a = bareVac
  · by_cases hb : b = bareVac
    · rw [ha, hb, Matrix.sum_apply]
      simp [hop_diag]
    · simp [brenesInit, hb]
  · simp [brenesInit, ha]

/-- The mass energy of the Brenes state in the sector of `ē`. -/
theorem mass_brenesInit (m : ℝ) (eb : Fin L → Fin 3) :
    ∑ n, (m * ∑ j : Fin L, (-1 : ℝ) ^ (j.val + 1) * occ (n j)) * Complex.normSq (brenesInit eb n) =
      ((3 : ℝ) ^ L)⁻¹ * (-(m * (((L + 1) / 2 : ℕ) : ℝ))) := by
  simp only [normSq_brenesInit, mul_ite, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  rw [mass_bareVac]
  ring

/-- The field energy of the Brenes state on the link `j` in the sector of `ē`. -/
theorem fieldSq_brenesInit (θ : ℝ) (e₀ : ℤ) (eb : Fin L → Fin 3) (j : Fin L) :
    fieldSq θ e₀ (brenesQ e₀ eb) (brenesInit eb) j =
      ((3 : ℝ) ^ L)⁻¹ * ((fieldOf (eb j) : ℝ) + θ) ^ 2 := by
  unfold fieldSq
  simp only [normSq_brenesInit, ite_mul, zero_mul, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
  rw [wField_brenes]
  simp [stagRho_bareVac]

/-- The initial field energy on every link: `⟨(E_j + ϑ)²⟩ = 2/3 + ϑ²`. -/
theorem fieldSqS_brenesInit (θ : ℝ) (e₀ : ℤ) (j : Fin L) :
    fieldSqS θ e₀ (brenesQ e₀) (brenesInit (L := L)) j = 2 / 3 + θ ^ 2 := by
  have hL : 0 < L := by have := j.isLt; omega
  unfold fieldSqS
  simp only [fieldSq_brenesInit, ← Finset.mul_sum, sum_eb_sq]
  rw [three_pow_eq hL]
  field_simp
  ring

/-- The initial field energy: `∑_j ⟨(E_j + ϑ)²⟩ = L (2/3 + ϑ²)`. -/
theorem sum_fieldSqS_brenesInit (θ : ℝ) (e₀ : ℤ) :
    ∑ j, fieldSqS θ e₀ (brenesQ e₀) (brenesInit (L := L)) j = L * (2 / 3 + θ ^ 2) := by
  simp only [fieldSqS_brenesInit, Finset.sum_const, Finset.card_univ, Fintype.card_fin,
    nsmul_eq_mul]

end BrenesU1Aux

/-- **The background term is exactly `ln 3`** for every cut (paper: Remark `rem:wilson`(a), Proposition
`prop:brenes`(c)). -/
theorem bgEntropyS_brenes (e₀ : ℤ) (ℓ : ℕ) (hℓ : 0 < ℓ) (hℓL : ℓ ≤ L) :
    bgEntropyS (brenesQ e₀) (brenesInit (L := L)) ℓ = Real.log 3 := by
  have hcong := WilsonSuperAux.lawEntropy_congr_inj
    (fun c : (Fin L → Fin 3) × Cfg L => bgSum (brenesQ e₀ c.1) ℓ)
    (fun c => fieldOf (c.1 ⟨ℓ - 1, by omega⟩)) (fun y => y - e₀)
    (fun a b h => by simpa using h) (fun c => brenesInit c.1 c.2)
    (fun c _ => bgSum_brenes e₀ c.1 ℓ hℓ hℓL)
  unfold bgEntropyS
  rw [hcong]
  exact BrenesU1Aux.lawEntropy_coord _

/-- The energy density of the Brenes state: `ε₀ = (g²/2)(2/3 + ϑ²) - m ⌈L/2⌉/L`. -/
theorem energyDensityS_brenes (w m g θ : ℝ) (e₀ : ℤ) (hL : 0 < L) :
    energyDensityS w m g θ e₀ (brenesQ e₀) (brenesInit (L := L)) =
      g ^ 2 / 2 * (2 / 3 + θ ^ 2) - m * (((L + 1) / 2 : ℕ) : ℝ) / L := by
  have hLpos : (0 : ℝ) < L := Nat.cast_pos.2 hL
  have hsec : ∀ eb : Fin L → Fin 3,
      (expect (cplx (schwinger w m g θ e₀ (brenesQ e₀ eb))) (brenesInit eb)).re =
        ((3 : ℝ) ^ L)⁻¹ * (-(m * (((L + 1) / 2 : ℕ) : ℝ))) +
          g ^ 2 / 2 * ∑ j, fieldSq θ e₀ (brenesQ e₀ eb) (brenesInit eb) j := by
    intro eb
    rw [WilsonAux.re_expect_schwinger, BrenesU1Aux.re_expect_hop_brenesInit,
      BrenesU1Aux.mass_brenesInit]
    ring
  have hF : ∑ eb : Fin L → Fin 3, ∑ j, fieldSq θ e₀ (brenesQ e₀ eb) (brenesInit eb) j =
      L * (2 / 3 + θ ^ 2) := by
    rw [Finset.sum_comm]
    exact BrenesU1Aux.sum_fieldSqS_brenesInit θ e₀
  unfold energyDensityS
  simp only [hsec]
  rw [Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ, BrenesU1Aux.card_eb,
    nsmul_eq_mul, ← Finset.mul_sum, hF]
  push_cast
  field_simp
  ring

/-- **Cut-averaged number entropy of the Brenes states** (paper: Proposition `prop:brenes`(c)): uniformly in
time and in `L`. -/
theorem brenes_numberEntropy_avg (w m g θ : ℝ) (hw : 0 ≤ w) (hg : 0 < g) (e₀ : ℤ) (hL : 0 < L)
    (t : ℝ) :
    (1 / (L : ℝ)) * ∑ ℓ ∈ Icc 1 L,
        shannon L (numDistS (segment ℓ) (wStateS w m g θ e₀ (brenesQ e₀) (brenesInit (L := L)) t)) ≤
      Real.log 3 + 1 / 2 * Real.log (2 * π * exp 1 *
        (2 / g ^ 2 * (g ^ 2 / 2 * (2 / 3 + θ ^ 2) - m * (((L + 1) / 2 : ℕ) : ℝ) / L + w + |m|) +
          1 / 12)) := by
  sorry

/-- **Electric-energy memory for all times** (paper: Proposition `prop:brenes`(d)): the electric energy density
stays within `4(w + |m|)/g²` of its initial value `2/3 + ϑ²`. -/
theorem brenes_field_memory (w m g θ : ℝ) (hw : 0 ≤ w) (hg : 0 < g) (e₀ : ℤ) (hL : 0 < L) (t : ℝ) :
    |(1 / (L : ℝ)) * ∑ j, fieldSqS θ e₀ (brenesQ e₀)
        (wStateS w m g θ e₀ (brenesQ e₀) (brenesInit (L := L)) t) j - (2 / 3 + θ ^ 2)| ≤
      4 * (w + |m|) / g ^ 2 := by
  sorry

end DFL
