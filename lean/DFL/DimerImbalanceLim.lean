import DFL.DimerImbalance

/-!
# The dimer imbalance vanishes in the thermodynamic limit (paper: Lemma `lem:fewapprox`,
Proposition `prop:dimerimbalance`(iii))

* `card_eig_near_le` (paper: Lemma `lem:fewapprox`): for every charge configuration `q` (no average), `h_q` has
  at most `2 (N/L + 1)` eigenvalues in `[E - ε, E + ε]` whenever `ε L (K + 1)^L < |J|`, with
  `K = (2|h| + |E|)/|J|`. A vector in the span of the corresponding eigenvectors that vanishes on the first
  two sites of every block of length `L` is `0`: on each block the eigenvalue recursion propagates the two
  zeros with errors controlled by `(K + 1)^L ‖(h_q - E) v‖ / |J|`.
* `dimerImbInf_le`: `Ī_N ≤ 8 L² (4|h|/|J| + 1)^{2L} / N + 4/L + 4/N` for every `L ≥ 1`, from the dimer identity
  of `prop:dimerimbalance`(i) (`dimer_stagWeight_identity`): eigenvectors away from the energies `±2h` have
  small staggered weight, and there are few eigenvalues near `±2h`.
* `dimerImbInf_tendsto_zero` (paper: `prop:dimerimbalance`(iii)): `Ī_N → 0` as `N → ∞` through even values.
  (Choosing `L ≍ log N` gives `Ī_N = O(1 / log N)`.)

This replaces the paper's earlier route through the density of states of the infinite-volume dimer model.
-/

open Matrix Finset Filter Topology

noncomputable section

namespace DFL

variable {N : ℕ}

namespace DimerImbalanceLimAux

open DimerImbalanceAux

/-- From two zeros, `|f(i+2)| ≤ |f i| + K|f(i+1)| + ρ(i+1)` gives `|f i| ≤ (K+1)^i ∑_{j≤i} ρ j`. -/
lemma rec_bound (f ρ : ℕ → ℝ) (K : ℝ) (hK : 0 ≤ K) (hρ : ∀ i, 0 ≤ ρ i)
    (h0 : f 0 = 0) (h1 : f 1 = 0)
    (hrec : ∀ i, |f (i + 2)| ≤ |f i| + K * |f (i + 1)| + ρ (i + 1)) (i : ℕ) :
    |f i| ≤ (K + 1) ^ i * ∑ j ∈ range (i + 1), ρ j ∧
      |f (i + 1)| ≤ (K + 1) ^ i * ∑ j ∈ range (i + 1), ρ j := by
  induction i with
  | zero =>
    simp only [h0, h1, abs_zero, pow_zero, one_mul, zero_add, sum_range_one]
    exact ⟨hρ 0, hρ 0⟩
  | succ i ih =>
    obtain ⟨ih1, ih2⟩ := ih
    have hT : 0 ≤ ∑ j ∈ range (i + 1), ρ j := sum_nonneg (fun j _ => hρ j)
    have hc : 1 ≤ (K + 1) ^ i := one_le_pow₀ (by linarith)
    have hρ1 := hρ (i + 1)
    have hr := hrec i
    rw [sum_range_succ, pow_succ]
    set T := ∑ j ∈ range (i + 1), ρ j
    set c := (K + 1) ^ i
    have hc0 : 0 ≤ c := by linarith
    constructor
    · nlinarith [mul_nonneg (mul_nonneg hc0 hK) hT, mul_nonneg (mul_nonneg hc0 hK) hρ1,
        mul_nonneg hc0 hρ1]
    · have e : |f (i + 1 + 1)| = |f (i + 2)| := rfl
      rw [e]
      calc |f (i + 2)| ≤ |f i| + K * |f (i + 1)| + ρ (i + 1) := hr
        _ ≤ c * T + K * (c * T) + ρ (i + 1) := by gcongr
        _ ≤ c * (K + 1) * (T + ρ (i + 1)) := by
          nlinarith [mul_nonneg hρ1 (show 0 ≤ c * K + c - 1 by nlinarith)]

/-- Block estimate: zeros at `s, s+1` and the recursion give `∑_{block} F² ≤ (L(K+1)^L)² ∑ ρ²`. -/
lemma block_bound (F ρ : ℕ → ℝ) (K : ℝ) (hK : 0 ≤ K) (hρ : ∀ i, 0 ≤ ρ i) (s L : ℕ)
    (h0 : F s = 0) (h1 : 1 < L → F (s + 1) = 0)
    (hrec : ∀ i, i + 2 < L →
      |F (s + i + 2)| ≤ |F (s + i)| + K * |F (s + i + 1)| + ρ (s + i + 1)) :
    ∑ i ∈ range L, F (s + i) ^ 2 ≤
      (L : ℝ) ^ 2 * ((K + 1) ^ L) ^ 2 * ∑ i ∈ range L, ρ (s + i) ^ 2 := by
  set f : ℕ → ℝ := fun i => if i < L then F (s + i) else 0 with hf
  set g : ℕ → ℝ := fun i => if i < L then ρ (s + i) else 0 with hg
  have hg0 : ∀ i, 0 ≤ g i := fun i => by
    simp only [hg]; split_ifs
    · exact hρ _
    · exact le_rfl
  have hf0 : f 0 = 0 := by
    simp only [hf]; split_ifs
    · simpa using h0
    · rfl
  have hf1 : f 1 = 0 := by
    simp only [hf]; split_ifs with h
    · exact h1 h
    · rfl
  have hfrec : ∀ i, |f (i + 2)| ≤ |f i| + K * |f (i + 1)| + g (i + 1) := by
    intro i
    have hnn : 0 ≤ |f i| + K * |f (i + 1)| + g (i + 1) := by
      have := hg0 (i + 1); positivity
    by_cases hi : i + 2 < L
    · have e2 : f (i + 2) = F (s + i + 2) := by simp only [hf, ite_eq_left hi]; ring_nf
      have e1 : f (i + 1) = F (s + i + 1) := by
        simp only [hf, ite_eq_left (show i + 1 < L by omega)]; ring_nf
      have e0 : f i = F (s + i) := by simp only [hf, ite_eq_left (show i < L by omega)]
      have eg : g (i + 1) = ρ (s + i + 1) := by
        simp only [hg, ite_eq_left (show i + 1 < L by omega)]; ring_nf
      rw [e2, e1, e0, eg]
      exact hrec i hi
    · have e2 : f (i + 2) = 0 := by simp only [hf, ite_eq_right hi]
      rw [e2, abs_zero]
      exact hnn
  set B := (K + 1) ^ L * ∑ j ∈ range L, ρ (s + j)
  have hK1 : 1 ≤ K + 1 := by linarith
  have hB : ∀ i ∈ range L, F (s + i) ^ 2 ≤ B ^ 2 := by
    intro i hi
    rw [mem_range] at hi
    have hb := (rec_bound f g K hK hg0 hf0 hf1 hfrec i).1
    have e0 : f i = F (s + i) := by simp only [hf, ite_eq_left hi]
    rw [e0] at hb
    have hsub : ∑ j ∈ range (i + 1), g j ≤ ∑ j ∈ range L, ρ (s + j) := by
      calc ∑ j ∈ range (i + 1), g j ≤ ∑ j ∈ range L, g j :=
            sum_le_sum_of_subset_of_nonneg (range_subset_range.mpr (by omega))
              (fun j _ _ => hg0 j)
        _ = ∑ j ∈ range L, ρ (s + j) :=
            sum_congr rfl fun j hj => by simp only [hg, ite_eq_left (mem_range.mp hj)]
    have hpow : (K + 1) ^ i ≤ (K + 1) ^ L := pow_le_pow_right₀ hK1 hi.le
    have hT0 : 0 ≤ ∑ j ∈ range (i + 1), g j := sum_nonneg fun j _ => hg0 j
    have hb2 : |F (s + i)| ≤ B :=
      hb.trans (mul_le_mul hpow hsub hT0 (by positivity))
    rw [← sq_abs]
    exact pow_le_pow_left₀ (abs_nonneg _) hb2 2
  calc ∑ i ∈ range L, F (s + i) ^ 2 ≤ ∑ _i ∈ range L, B ^ 2 := sum_le_sum hB
    _ = (L : ℝ) * (((K + 1) ^ L) ^ 2 * (∑ j ∈ range L, ρ (s + j)) ^ 2) := by
        rw [sum_const, card_range, nsmul_eq_mul]; ring
    _ ≤ (L : ℝ) * (((K + 1) ^ L) ^ 2 * ((L : ℝ) * ∑ j ∈ range L, ρ (s + j) ^ 2)) := by
        gcongr
        have := sq_sum_le_card_mul_sum_sq (s := range L) (f := fun j => ρ (s + j))
        simpa using this
    _ = (L : ℝ) ^ 2 * ((K + 1) ^ L) ^ 2 * ∑ i ∈ range L, ρ (s + i) ^ 2 := by ring

/-- `∑_{n < M L} f n = ∑_{b < M} ∑_{i < L} f (b L + i)`. -/
lemma sum_range_mul_eq (f : ℕ → ℝ) (M L : ℕ) :
    ∑ n ∈ range (M * L), f n = ∑ b ∈ range M, ∑ i ∈ range L, f (b * L + i) := by
  induction M with
  | zero => simp
  | succ M ih => rw [Nat.succ_mul, sum_range_add, ih, sum_range_succ]

/-- `(h_q v)(x) = -J (v(x-1) + v(x+1)) + 2h q_x v(x)` with zero boundary values. -/
lemma hq_mulVec_apply (J h : ℝ) (q : Cfg N) (v : Fin N → ℝ) (x : Fin N) :
    (hq J h q *ᵥ v) x =
      -J * (prev0 (ext0 v) x.val + ext0 v (x.val + 1)) + 2 * h * sgn (q x) * ext0 v x.val := by
  rw [DFL.hq, add_mulVec, smul_mulVec, smul_mulVec, Pi.add_apply, Pi.smul_apply, Pi.smul_apply,
    mulVec_diagonal, pathAdj_mulVec_apply, smul_eq_mul, smul_eq_mul, ext0_val]
  ring

/-- At most `2 (N/L + 1)` sites of `Fin N` are among the first two sites of an `L`-block. -/
lemma card_firstTwo_le (L : ℕ) :
    ((univ.filter (fun y : Fin N => y.val % L = 0 ∨ y.val % L = 1)).card : ℝ) ≤
      2 * ((N : ℝ) / L + 1) := by
  have h1 : (univ.filter (fun y : Fin N => y.val % L = 0 ∨ y.val % L = 1)).card ≤
      (range (N / L + 1) ×ˢ range 2).card := by
    apply card_le_card_of_injOn (fun y : Fin N => (y.val / L, y.val % L))
    · intro y hy
      simp only [coe_filter, mem_univ, true_and, Set.mem_ofPred_eq] at hy
      simp only [coe_product, coe_range, Set.mem_prod, Set.mem_Iio]
      constructor
      · have := Nat.div_le_div_right (c := L) y.isLt.le
        omega
      · omega
    · intro y _ z _ hyz
      simp only [Prod.mk.injEq] at hyz
      apply Fin.ext
      rw [← Nat.div_add_mod y.val L, ← Nat.div_add_mod z.val L, hyz.1, hyz.2]
  rw [card_product, card_range, card_range] at h1
  have h2 : ((N / L : ℕ) : ℝ) ≤ (N : ℝ) / L := Nat.cast_div_le
  calc ((univ.filter (fun y : Fin N => y.val % L = 0 ∨ y.val % L = 1)).card : ℝ)
      ≤ ((N / L + 1) * 2 : ℕ) := by exact_mod_cast h1
    _ = 2 * (((N / L : ℕ) : ℝ) + 1) := by push_cast; ring
    _ ≤ 2 * ((N : ℝ) / L + 1) := by linarith

/-- If `|P| < |S|`, a nonzero `c` supported on `S` has `∑_α c_α Φ_α` vanishing on `P`. -/
lemma exists_coeff (Φ : Fin N → Fin N → ℝ) (P S : Finset (Fin N)) (h : P.card < S.card) :
    ∃ c : Fin N → ℝ, c ≠ 0 ∧ (∀ y ∈ P, ∑ α, c α * Φ α y = 0) ∧ (∀ α ∉ S, c α = 0) := by
  let T : (Fin N → ℝ) →ₗ[ℝ] ((P → ℝ) × (↥Sᶜ → ℝ)) :=
    { toFun := fun c => (fun y => ∑ α, c α * Φ α y.val, fun α => c α.val)
      map_add' := fun c d => by
        ext y
        · simp [add_mul, sum_add_distrib]
        · simp
      map_smul' := fun a c => by
        ext y
        · simp [mul_sum, mul_assoc]
        · simp }
  have hrank : Module.finrank ℝ ((P → ℝ) × (↥Sᶜ → ℝ)) < Module.finrank ℝ (Fin N → ℝ) := by
    rw [Module.finrank_prod, Module.finrank_fintype_fun_eq_card,
      Module.finrank_fintype_fun_eq_card, Module.finrank_fin_fun, Fintype.card_coe,
      Fintype.card_coe, card_compl, Fintype.card_fin]
    have := card_le_univ S
    rw [Fintype.card_fin] at this
    omega
  obtain ⟨c, hc, hc0⟩ :=
    (Submodule.ne_bot_iff _).mp (LinearMap.ker_ne_bot_of_finrank_lt (f := T) hrank)
  refine ⟨c, hc0, fun y hy => ?_, fun α hα => ?_⟩
  · have := congrArg (fun p => p.1 ⟨y, hy⟩) (LinearMap.mem_ker.mp hc)
    simpa [T] using this
  · have := congrArg (fun p => p.2 ⟨α, mem_compl.mpr hα⟩) (LinearMap.mem_ker.mp hc)
    simpa [T] using this

/-- Parseval in the eigenbasis: `‖∑_α c_α φ_α‖² = ∑_α c_α²`. -/
lemma sum_sq_comb (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (c : Fin N → ℝ) :
    ∑ y, (∑ α, c α * hA.eigenvectorBasis α y) ^ 2 = ∑ α, c α ^ 2 := by
  have e : ∀ y, (∑ α, c α * hA.eigenvectorBasis α y) ^ 2 =
      ∑ α, ∑ β, c α * c β * (hA.eigenvectorBasis α y * hA.eigenvectorBasis β y) := by
    intro y
    rw [sq, sum_mul_sum]
    exact sum_congr rfl fun α _ => sum_congr rfl fun β _ => by ring
  simp_rw [e]
  rw [sum_comm]
  refine sum_congr rfl fun α _ => ?_
  rw [sum_comm]
  simp_rw [← mul_sum, ImbalanceAux.eig_orth A hA]
  simp [sq]

/-- `A (∑_α c_α φ_α) = ∑_α c_α E_α φ_α`. -/
lemma mulVec_comb (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (c : Fin N → ℝ)
    (y : Fin N) :
    (A *ᵥ fun x => ∑ α, c α * hA.eigenvectorBasis α x) y =
      ∑ α, c α * hA.eigenvalues α * hA.eigenvectorBasis α y := by
  simp only [mulVec, dotProduct, mul_sum]
  rw [sum_comm]
  refine sum_congr rfl fun α _ => ?_
  have e := ImbalanceAux.eig_eq A hA α y
  calc ∑ x, A y x * (c α * hA.eigenvectorBasis α x) =
        c α * ∑ x, A y x * hA.eigenvectorBasis α x := by
        rw [mul_sum]; exact sum_congr rfl fun x _ => by ring
    _ = _ := by rw [e]; ring

/-- The extension by zero at an index `n < N`. -/
lemma ext0_mk (φ : Fin N → ℝ) (n : ℕ) (hn : n < N) : ext0 φ n = φ ⟨n, hn⟩ :=
  ext0_val φ ⟨n, hn⟩

/-- With `r = (h_q - E) v`: `|v(x+1)| ≤ |v(x-1)| + K |v x| + |r x|/|J|`, `K = (2|h|+|E|)/|J|`. -/
lemma rec_step (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) (E K : ℝ) (hK : K = (2 * |h| + |E|) / |J|)
    (v r : Fin N → ℝ) (hr : ∀ y, r y = (hq J h q *ᵥ v) y - E * v y) (x : ℕ) (hx : x < N) :
    |ext0 v (x + 1)| ≤ |prev0 (ext0 v) x| + K * |ext0 v x| + |ext0 r x| / |J| := by
  have hJ0 : 0 < |J| := abs_pos.mpr hJ
  have e := hq_mulVec_apply J h q v ⟨x, hx⟩
  have hGx : ext0 r x = (hq J h q *ᵥ v) ⟨x, hx⟩ - E * v ⟨x, hx⟩ := by
    rw [ext0_mk r x hx, hr]
  have hFx : ext0 v x = v ⟨x, hx⟩ := ext0_mk v x hx
  simp only at e
  have key : J * ext0 v (x + 1) = -J * prev0 (ext0 v) x +
      (2 * h * sgn (q ⟨x, hx⟩) - E) * ext0 v x - ext0 r x := by
    linear_combination hGx + e + E * hFx
  have hs : |2 * h * sgn (q ⟨x, hx⟩) - E| ≤ 2 * |h| + |E| := by
    calc |2 * h * sgn (q ⟨x, hx⟩) - E| ≤ |2 * h * sgn (q ⟨x, hx⟩)| + |E| := abs_sub _ _
      _ = 2 * |h| + |E| := by rw [abs_mul, abs_mul, ImbalanceAux.abs_sgn, abs_two, mul_one]
  have h1 : |J| * |ext0 v (x + 1)| ≤
      |J| * |prev0 (ext0 v) x| + (2 * |h| + |E|) * |ext0 v x| + |ext0 r x| := by
    rw [← abs_mul, key]
    calc |-J * prev0 (ext0 v) x + (2 * h * sgn (q ⟨x, hx⟩) - E) * ext0 v x - ext0 r x|
        ≤ |-J * prev0 (ext0 v) x + (2 * h * sgn (q ⟨x, hx⟩) - E) * ext0 v x| + |ext0 r x| :=
          abs_sub _ _
      _ ≤ |-J * prev0 (ext0 v) x| + |(2 * h * sgn (q ⟨x, hx⟩) - E) * ext0 v x| + |ext0 r x| := by
          gcongr; exact abs_add_le _ _
      _ = |J| * |prev0 (ext0 v) x| + |2 * h * sgn (q ⟨x, hx⟩) - E| * |ext0 v x| + |ext0 r x| := by
          rw [abs_mul, abs_mul, abs_neg]
      _ ≤ |J| * |prev0 (ext0 v) x| + (2 * |h| + |E|) * |ext0 v x| + |ext0 r x| := by gcongr
  have hKJ : K * |J| = 2 * |h| + |E| := by rw [hK]; field_simp
  have hGJ : |J| * (|ext0 r x| / |J|) = |ext0 r x| := mul_div_cancel₀ _ hJ0.ne'
  refine le_of_mul_le_mul_left (h1.trans (le_of_eq ?_)) hJ0
  linear_combination (-|ext0 v x|) * hKJ - hGJ

/-- If `v` is zero at the first two sites of every block, `‖v‖² ≤ (L(K+1)^L/|J|)² ‖(h_q-E)v‖²`. -/
lemma sum_sq_le_of_vanish (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) (E : ℝ) (L : ℕ) (hL : 0 < L)
    (v : Fin N → ℝ) (hv : ∀ y : Fin N, (y.val % L = 0 ∨ y.val % L = 1) → v y = 0) :
    ∑ y, v y ^ 2 ≤ ((L : ℝ) * ((2 * |h| + |E|) / |J| + 1) ^ L / |J|) ^ 2 *
      ∑ y, ((hq J h q *ᵥ v) y - E * v y) ^ 2 := by
  have hJ0 : 0 < |J| := abs_pos.mpr hJ
  obtain ⟨K, hK⟩ : ∃ K : ℝ, K = (2 * |h| + |E|) / |J| := ⟨_, rfl⟩
  rw [← hK]
  have hK0 : 0 ≤ K := by rw [hK]; positivity
  obtain ⟨r, hr⟩ : ∃ r : Fin N → ℝ, r = fun y => (hq J h q *ᵥ v) y - E * v y := ⟨_, rfl⟩
  have hr' : ∀ y, r y = (hq J h q *ᵥ v) y - E * v y := fun y => by rw [hr]
  rw [show ∑ y, ((hq J h q *ᵥ v) y - E * v y) ^ 2 = ∑ y, r y ^ 2 by simp only [hr']]
  set ρ : ℕ → ℝ := fun x => |ext0 r x| / |J| with hρ
  have hρ0 : ∀ x, 0 ≤ ρ x := fun x => by positivity
  have hrec : ∀ x : ℕ, x < N → |ext0 v (x + 1)| ≤ |prev0 (ext0 v) x| + K * |ext0 v x| + ρ x :=
    fun x hx => rec_step J h hJ q E K hK v r hr' x hx
  -- the blocks `[bL, bL + L)`, `b < M`, cover `[0, N)`
  set M := N / L + 1 with hM
  have hNM : N ≤ M * L := by
    have := Nat.lt_div_mul_add (a := N) hL
    rw [hM, Nat.add_mul, one_mul]
    omega
  have hsum_ext : ∀ f : Fin N → ℝ,
      ∑ y, f y ^ 2 = ∑ b ∈ range M, ∑ i ∈ range L, ext0 f (b * L + i) ^ 2 := by
    intro f
    rw [← sum_range_mul_eq (fun n => ext0 f n ^ 2) M L]
    rw [show ∑ y, f y ^ 2 = ∑ y : Fin N, ext0 f y.val ^ 2 by simp only [ext0_val]]
    rw [Fin.sum_univ_eq_sum_range (fun n => ext0 f n ^ 2) N]
    apply sum_subset (range_subset_range.mpr hNM)
    intro n _ hn'
    rw [mem_range, not_lt] at hn'
    rw [ext0_of_le f hn']
    ring
  have hblock : ∀ b, ∑ i ∈ range L, ext0 v (b * L + i) ^ 2 ≤
      (L : ℝ) ^ 2 * ((K + 1) ^ L) ^ 2 * ∑ i ∈ range L, ρ (b * L + i) ^ 2 := by
    intro b
    apply block_bound (ext0 v) ρ K hK0 hρ0 (b * L) L
    · by_cases hb : b * L < N
      · rw [ext0_mk v _ hb]
        exact hv _ (Or.inl (by simp))
      · rw [ext0_of_le v (not_lt.mp hb)]
    · intro hL1
      by_cases hb : b * L + 1 < N
      · rw [ext0_mk v _ hb]
        exact hv _ (Or.inr (Nat.mul_add_mod_of_lt hL1))
      · rw [ext0_of_le v (not_lt.mp hb)]
    · intro i _
      by_cases hb : b * L + i + 1 < N
      · exact hrec (b * L + i + 1) hb
      · rw [ext0_of_le v (by omega : N ≤ b * L + i + 2), abs_zero]
        have := hρ0 (b * L + i + 1)
        positivity
  calc ∑ y, v y ^ 2 = ∑ b ∈ range M, ∑ i ∈ range L, ext0 v (b * L + i) ^ 2 := hsum_ext v
    _ ≤ ∑ b ∈ range M, (L : ℝ) ^ 2 * ((K + 1) ^ L) ^ 2 * ∑ i ∈ range L, ρ (b * L + i) ^ 2 :=
        sum_le_sum fun b _ => hblock b
    _ = (L : ℝ) ^ 2 * ((K + 1) ^ L) ^ 2 / J ^ 2 *
          ∑ b ∈ range M, ∑ i ∈ range L, ext0 r (b * L + i) ^ 2 := by
        rw [mul_sum]
        refine sum_congr rfl fun b _ => ?_
        rw [mul_sum, mul_sum]
        refine sum_congr rfl fun i _ => ?_
        simp only [hρ, div_pow, sq_abs]
        ring
    _ = ((L : ℝ) * (K + 1) ^ L / |J|) ^ 2 * ∑ y, r y ^ 2 := by
        rw [← hsum_ext r, div_pow, mul_pow, sq_abs]

end DimerImbalanceLimAux

/-- **Few eigenvalues near any energy** (paper: Lemma `lem:fewapprox`): if `ε L (K + 1)^L < |J|` with
`K = (2|h| + |E|)/|J|`, then `h_q` has at most `2 (N/L + 1)` eigenvalues in `[E - ε, E + ε]`. -/
theorem card_eig_near_le (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) (E ε : ℝ) (L : ℕ) (hL : 0 < L)
    (hε : ε * L * ((2 * |h| + |E|) / |J| + 1) ^ L < |J|) :
    ((univ.filter (fun α : Fin N => |(hq_isHermitian J h q).eigenvalues α - E| ≤ ε)).card : ℝ) ≤
      2 * ((N : ℝ) / L + 1) := by
  classical
  set hH := hq_isHermitian J h q with hHdef
  set S := univ.filter (fun α : Fin N => |hH.eigenvalues α - E| ≤ ε) with hS
  set P : Finset (Fin N) := univ.filter (fun y : Fin N => y.val % L = 0 ∨ y.val % L = 1) with hP
  have hPc : (P.card : ℝ) ≤ 2 * ((N : ℝ) / L + 1) := DimerImbalanceLimAux.card_firstTwo_le L
  suffices hSP : S.card ≤ P.card by
    calc (S.card : ℝ) ≤ P.card := by exact_mod_cast hSP
      _ ≤ _ := hPc
  by_contra hlt
  rw [not_le] at hlt
  obtain ⟨c, hc0, hcP, hcS⟩ :=
    DimerImbalanceLimAux.exists_coeff (fun α y => hH.eigenvectorBasis α y) P S hlt
  -- `ε ≥ 0` since `S` is nonempty
  have hε0 : 0 ≤ ε := by
    obtain ⟨α, hα⟩ : S.Nonempty := by
      rw [← card_pos]
      omega
    exact (abs_nonneg _).trans (mem_filter.mp hα).2
  -- the combination `v = ∑_{α ∈ S} c_α φ_α` vanishes on `P`
  set v : Fin N → ℝ := fun y => ∑ α, c α * hH.eigenvectorBasis α y with hvdef
  have hv : ∀ y : Fin N, (y.val % L = 0 ∨ y.val % L = 1) → v y = 0 :=
    fun y hy => hcP y (mem_filter.mpr ⟨mem_univ _, hy⟩)
  have hA := DimerImbalanceLimAux.sum_sq_le_of_vanish J h hJ q E L hL v hv
  have hnorm : ∑ y, v y ^ 2 = ∑ α, c α ^ 2 := DimerImbalanceLimAux.sum_sq_comb _ hH c
  -- `‖(h_q - E) v‖ ≤ ε ‖v‖`
  have hr : ∑ y, ((hq J h q *ᵥ v) y - E * v y) ^ 2 ≤ ε ^ 2 * ∑ α, c α ^ 2 := by
    have e : ∀ y, (hq J h q *ᵥ v) y - E * v y =
        ∑ α, (c α * (hH.eigenvalues α - E)) * hH.eigenvectorBasis α y := by
      intro y
      rw [DimerImbalanceLimAux.mulVec_comb _ hH c y]
      simp only [hvdef, mul_sum, ← sum_sub_distrib]
      exact sum_congr rfl fun α _ => by ring
    simp_rw [e]
    rw [DimerImbalanceLimAux.sum_sq_comb _ hH, mul_sum]
    refine sum_le_sum fun α _ => ?_
    by_cases hα : α ∈ S
    · have h1 : |hH.eigenvalues α - E| ≤ ε := (mem_filter.mp hα).2
      have h2 : (hH.eigenvalues α - E) ^ 2 ≤ ε ^ 2 := by
        rw [← sq_abs]; exact pow_le_pow_left₀ (abs_nonneg _) h1 2
      calc (c α * (hH.eigenvalues α - E)) ^ 2 = c α ^ 2 * (hH.eigenvalues α - E) ^ 2 := by ring
        _ ≤ c α ^ 2 * ε ^ 2 := by gcongr
        _ = ε ^ 2 * c α ^ 2 := by ring
    · rw [hcS α hα]; simp
  have hpos : 0 < ∑ α, c α ^ 2 := by
    obtain ⟨α, hα⟩ := Function.ne_iff.mp hc0
    have hα' : c α ≠ 0 := hα
    exact lt_of_lt_of_le (by positivity) (single_le_sum (fun β _ => sq_nonneg (c β)) (mem_univ α))
  have hJ0 : 0 < |J| := abs_pos.mpr hJ
  set D := (L : ℝ) * ((2 * |h| + |E|) / |J| + 1) ^ L / |J| with hD
  have hDε : D * ε < 1 := by
    rw [hD, div_mul_eq_mul_div, div_lt_one hJ0]
    linarith [hε]
  have hDε0 : 0 ≤ D * ε := by
    rw [hD]
    positivity
  have key : ∑ α, c α ^ 2 ≤ (D * ε) ^ 2 * ∑ α, c α ^ 2 := by
    calc ∑ α, c α ^ 2 = ∑ y, v y ^ 2 := hnorm.symm
      _ ≤ D ^ 2 * ∑ y, ((hq J h q *ᵥ v) y - E * v y) ^ 2 := hA
      _ ≤ D ^ 2 * (ε ^ 2 * ∑ α, c α ^ 2) := by gcongr
      _ = (D * ε) ^ 2 * ∑ α, c α ^ 2 := by ring
  have hlt1 : (D * ε) ^ 2 < 1 := by nlinarith
  nlinarith

namespace DimerImbalanceLimAux

open DimerImbalanceAux

/-- The staggered weight of a normalized eigenvector satisfies `|a_α| ≤ 1`. -/
lemma abs_stagWeight_le (J h : ℝ) (q : Cfg N) (α : Fin N) : |stagWeight J h q α| ≤ 1 := by
  unfold stagWeight
  calc |∑ x, (-1 : ℝ) ^ x.val * ((hq_isHermitian J h q).eigenvectorBasis α x) ^ 2|
      ≤ ∑ x, |(-1 : ℝ) ^ x.val * ((hq_isHermitian J h q).eigenvectorBasis α x) ^ 2| :=
        abs_sum_le_sum_abs _ _
    _ = ∑ x, ((hq_isHermitian J h q).eigenvectorBasis α x) ^ 2 := by
        refine sum_congr rfl fun x _ => ?_
        rw [abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul, abs_pow, sq_abs]
    _ = 1 := ImbalanceAux.sum_eig_sq' _ _ α

/-- Per-sector bound for a dimer `q`: `Ī_q ≤ 8 L² (4|h|/|J| + 1)^{2L} / N + 4/L + 4/N`. -/
lemma imbInf_le_of_dimer (J h : ℝ) (hJ : J ≠ 0) (hN : Even N) (hN0 : 0 < N) (L : ℕ)
    (hL : 0 < L) (q : Cfg N) (hq : IsDimer q) :
    imbInf J h q ≤ 8 * (L : ℝ) ^ 2 * (4 * |h| / |J| + 1) ^ (2 * L) / N + 4 / L + 4 / N := by
  set hH := hq_isHermitian J h q with hHdef
  have hJ0 : 0 < |J| := abs_pos.mpr hJ
  obtain ⟨K, hK⟩ : ∃ K : ℝ, K = 4 * |h| / |J| := ⟨_, rfl⟩
  rw [← hK]
  have hK0 : 0 ≤ K := by rw [hK]; positivity
  have hL0 : (0 : ℝ) < L := by exact_mod_cast hL
  have hN0' : (0 : ℝ) < N := by exact_mod_cast hN0
  have hX1 : 1 ≤ (K + 1) ^ L := one_le_pow₀ (by linarith)
  obtain ⟨ε, hε⟩ : ∃ ε : ℝ, ε = |J| / (2 * L * (K + 1) ^ L) := ⟨_, rfl⟩
  have hε0 : 0 < ε := by rw [hε]; positivity
  have hεL : ε * L * (K + 1) ^ L = |J| / 2 := by
    rw [hε]; field_simp
  set S₁ := univ.filter (fun α : Fin N => |hH.eigenvalues α - 2 * h| ≤ ε) with hS₁
  set S₂ := univ.filter (fun α : Fin N => |hH.eigenvalues α - -(2 * h)| ≤ ε) with hS₂
  have hc₁ : (S₁.card : ℝ) ≤ 2 * ((N : ℝ) / L + 1) := by
    apply card_eig_near_le J h hJ q (2 * h) ε L hL
    have e : (2 * |h| + |2 * h|) / |J| + 1 = K + 1 := by
      rw [hK, abs_mul, abs_two]; ring
    rw [e, hεL]; linarith
  have hc₂ : (S₂.card : ℝ) ≤ 2 * ((N : ℝ) / L + 1) := by
    apply card_eig_near_le J h hJ q (-(2 * h)) ε L hL
    have e : (2 * |h| + |-(2 * h)|) / |J| + 1 = K + 1 := by
      rw [hK, abs_neg, abs_mul, abs_two]; ring
    rw [e, hεL]; linarith
  -- pointwise bound on `a_α²`
  have hpt : ∀ α, stagWeight J h q α ^ 2 ≤
      (if α ∈ S₁ then 1 else 0) + (if α ∈ S₂ then 1 else 0) +
        J ^ 2 / ε ^ 2 * (hH.eigenvectorBasis α ⟨0, hN0⟩ ^ 2 +
          hH.eigenvectorBasis α ⟨N - 1, by omega⟩ ^ 2) := by
    intro α
    have ha1 := abs_stagWeight_le J h q α
    have hsq : stagWeight J h q α ^ 2 ≤ |stagWeight J h q α| := by
      rw [← sq_abs]; nlinarith [abs_nonneg (stagWeight J h q α)]
    have hnn : 0 ≤ J ^ 2 / ε ^ 2 * (hH.eigenvectorBasis α ⟨0, hN0⟩ ^ 2 +
        hH.eigenvectorBasis α ⟨N - 1, by omega⟩ ^ 2) := by positivity
    have hi1 : 0 ≤ (if α ∈ S₁ then (1 : ℝ) else 0) := by split_ifs <;> norm_num
    have hi2 : 0 ≤ (if α ∈ S₂ then (1 : ℝ) else 0) := by split_ifs <;> norm_num
    by_cases h1 : α ∈ S₁
    · rw [ite_eq_left h1]
      nlinarith
    by_cases h2 : α ∈ S₂
    · rw [ite_eq_left h2]
      nlinarith
    rw [ite_eq_right h1, ite_eq_right h2, zero_add, zero_add]
    have h1' : ε < |hH.eigenvalues α - 2 * h| := by
      simpa [hS₁] using h1
    have h2' : ε < |hH.eigenvalues α + 2 * h| := by
      simpa [hS₂, sub_neg_eq_add] using h2
    have hprod : ε ^ 2 < |(2 * h) ^ 2 - hH.eigenvalues α ^ 2| := by
      have e : |(2 * h) ^ 2 - hH.eigenvalues α ^ 2| =
          |hH.eigenvalues α - 2 * h| * |hH.eigenvalues α + 2 * h| := by
        rw [← abs_mul, ← abs_neg]; congr 1; ring
      rw [e, sq]
      exact mul_lt_mul'' h1' h2' hε0.le hε0.le
    have hid := dimer_stagWeight_identity J h hN hN0 q hq α
    have habs : |(2 * h) ^ 2 - hH.eigenvalues α ^ 2| * |stagWeight J h q α| ≤
        J ^ 2 * (hH.eigenvectorBasis α ⟨0, hN0⟩ ^ 2 +
          hH.eigenvectorBasis α ⟨N - 1, by omega⟩ ^ 2) := by
      rw [← abs_mul, hid, abs_mul, abs_of_nonneg (sq_nonneg J)]
      gcongr
      refine (abs_sub _ _).trans (le_of_eq ?_)
      rw [abs_of_nonneg (sq_nonneg _), abs_of_nonneg (sq_nonneg _)]
    have hfin : |stagWeight J h q α| ≤ J ^ 2 / ε ^ 2 * (hH.eigenvectorBasis α ⟨0, hN0⟩ ^ 2 +
        hH.eigenvectorBasis α ⟨N - 1, by omega⟩ ^ 2) := by
      rw [div_mul_eq_mul_div, le_div_iff₀ (by positivity)]
      nlinarith [abs_nonneg (stagWeight J h q α)]
    linarith
  have hsum : ∑ α, stagWeight J h q α ^ 2 ≤ S₁.card + S₂.card + J ^ 2 / ε ^ 2 * 2 := by
    calc ∑ α, stagWeight J h q α ^ 2 ≤ ∑ α, ((if α ∈ S₁ then 1 else 0) +
          (if α ∈ S₂ then 1 else 0) + J ^ 2 / ε ^ 2 * (hH.eigenvectorBasis α ⟨0, hN0⟩ ^ 2 +
            hH.eigenvectorBasis α ⟨N - 1, by omega⟩ ^ 2)) := sum_le_sum fun α _ => hpt α
      _ = S₁.card + S₂.card + J ^ 2 / ε ^ 2 * (∑ α, hH.eigenvectorBasis α ⟨0, hN0⟩ ^ 2 +
            ∑ α, hH.eigenvectorBasis α ⟨N - 1, by omega⟩ ^ 2) := by
        rw [sum_add_distrib, sum_add_distrib, ← mul_sum, sum_add_distrib]
        simp
      _ = S₁.card + S₂.card + J ^ 2 / ε ^ 2 * 2 := by
        rw [ImbalanceAux.sum_eig_sq, ImbalanceAux.sum_eig_sq]; ring
  have hJε : J ^ 2 / ε ^ 2 * 2 = 8 * (L : ℝ) ^ 2 * ((K + 1) ^ L) ^ 2 := by
    rw [hε, div_pow, mul_pow, mul_pow, sq_abs]
    field_simp
    ring
  rw [imbInf]
  calc (1 / N : ℝ) * ∑ α, stagWeight J h q α ^ 2
      ≤ (1 / N) * (S₁.card + S₂.card + J ^ 2 / ε ^ 2 * 2) := by gcongr
    _ ≤ (1 / N) * (2 * ((N : ℝ) / L + 1) + 2 * ((N : ℝ) / L + 1) +
          8 * (L : ℝ) ^ 2 * ((K + 1) ^ L) ^ 2) := by
        apply mul_le_mul_of_nonneg_left _ (by positivity)
        linarith
    _ = 8 * (L : ℝ) ^ 2 * (K + 1) ^ (2 * L) / N + 4 / L + 4 / N := by
        rw [pow_mul']
        field_simp
        ring

end DimerImbalanceLimAux

/-- **Quantitative vanishing of the dimer imbalance**: for every `L ≥ 1`,
`Ī_N ≤ 8 L² (4|h|/|J| + 1)^{2L} / N + 4/L + 4/N`. -/
theorem dimerImbInf_le (J h : ℝ) (hJ : J ≠ 0) (hN : Even N) (hN0 : 0 < N) (L : ℕ) (hL : 0 < L) :
    (1 / 2 ^ (N / 2) : ℝ) * ∑ q ∈ univ.filter (fun q : Cfg N => IsDimer q), imbInf J h q ≤
      8 * (L : ℝ) ^ 2 * (4 * |h| / |J| + 1) ^ (2 * L) / N + 4 / L + 4 / N := by
  have hcard : ((univ.filter (fun q : Cfg N => IsDimer q)).card : ℝ) = 2 ^ (N / 2) := by
    rw [card_isDimer hN]; push_cast; ring
  calc (1 / 2 ^ (N / 2) : ℝ) * ∑ q ∈ univ.filter (fun q : Cfg N => IsDimer q), imbInf J h q
      ≤ (1 / 2 ^ (N / 2) : ℝ) * ∑ _q ∈ univ.filter (fun q : Cfg N => IsDimer q),
          (8 * (L : ℝ) ^ 2 * (4 * |h| / |J| + 1) ^ (2 * L) / N + 4 / L + 4 / N) := by
        apply mul_le_mul_of_nonneg_left _ (by positivity)
        exact sum_le_sum fun q hq => DimerImbalanceLimAux.imbInf_le_of_dimer J h hJ hN hN0 L hL q
          (mem_filter.mp hq).2
    _ = 8 * (L : ℝ) ^ 2 * (4 * |h| / |J| + 1) ^ (2 * L) / N + 4 / L + 4 / N := by
        rw [sum_const, nsmul_eq_mul, hcard]
        field_simp

/-- **The imbalance of the dimer CDW state vanishes** (paper: Proposition `prop:dimerimbalance`(iii)):
`Ī_N → 0` as `N → ∞` through even values. -/
theorem dimerImbInf_tendsto_zero (J h : ℝ) (hJ : J ≠ 0) :
    Tendsto (fun M : ℕ => (1 / 2 ^ (2 * M / 2) : ℝ) *
      ∑ q ∈ univ.filter (fun q : Cfg (2 * M) => IsDimer q), imbInf J h q) atTop (𝓝 0) := by
  rw [Metric.tendsto_atTop]
  intro η hη
  obtain ⟨L, hL⟩ := exists_nat_gt (8 / η)
  have hL0' : (0 : ℝ) < L := lt_trans (by positivity) hL
  have hL0 : 0 < L := by exact_mod_cast hL0'
  obtain ⟨C, hC⟩ : ∃ C : ℝ, C = 8 * (L : ℝ) ^ 2 * (4 * |h| / |J| + 1) ^ (2 * L) + 4 := ⟨_, rfl⟩
  have hC0 : 0 < C := by rw [hC]; positivity
  obtain ⟨M₀, hM₀⟩ := exists_nat_gt (C / η)
  refine ⟨M₀ + 1, fun M hM => ?_⟩
  have hM0 : 0 < M := by omega
  have hM0' : (0 : ℝ) < M := by exact_mod_cast hM0
  have hMM : (M₀ : ℝ) + 1 ≤ M := by exact_mod_cast hM
  have hbound := dimerImbInf_le (N := 2 * M) J h hJ (even_two_mul M) (by omega) L hL0
  have hnn : 0 ≤ (1 / 2 ^ (2 * M / 2) : ℝ) *
      ∑ q ∈ univ.filter (fun q : Cfg (2 * M) => IsDimer q), imbInf J h q := by
    apply mul_nonneg (by positivity)
    refine sum_nonneg fun q _ => ?_
    unfold imbInf
    exact mul_nonneg (by positivity) (sum_nonneg fun _ _ => sq_nonneg _)
  rw [Real.dist_eq, sub_zero, abs_of_nonneg hnn]
  refine lt_of_le_of_lt hbound ?_
  push_cast
  have h1 : 4 / (L : ℝ) < η / 2 := by
    rw [div_lt_iff₀ hL0']
    rw [div_lt_iff₀ hη] at hL
    linarith
  have h2 : C / (2 * M) < η / 2 := by
    rw [div_lt_iff₀ (by positivity)]
    rw [div_lt_iff₀ hη] at hM₀
    nlinarith
  have e : 8 * (L : ℝ) ^ 2 * (4 * |h| / |J| + 1) ^ (2 * L) / (2 * M) + 4 / L + 4 / (2 * M) =
      C / (2 * M) + 4 / L := by
    rw [hC]; ring
  rw [e]
  linarith

end DFL
