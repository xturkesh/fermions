import DFL.Sector

/-!
# Disorder-free localization depends on the initial state (paper: Section `sec:dimer`)

A product gauge state with the even links polarized along `z` and the odd links in prescribed
`σˣ` eigenstates (labels alternating along the chain), together with density-wave matter, induces
the **random dimer** sector law: the charges come in equal pairs `q_{2k} = q_{2k+1}`, and the
`2^{N/2}` dimer configurations are equally likely. (The resulting one-body problem is the random
dimer model, which has critical energies `± 2h` for `0 < 2h < J`; the transport lower bounds of
Jitomirskaya–Schulz-Baldes–Stolz are infinite-volume statements and are not formalized here.)
-/

open Matrix Finset
open scoped Kronecker

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Prescribed label of the odd link `b = 2k+1`: `(k+1) mod 2`. -/
def dimerLabel (b : Fin N) : ZMod 2 := ((b.val / 2 + 1 : ℕ) : ZMod 2)

/-- The link configurations allowed by the dimer gauge state. -/
def dimerLinks (x : Cfg N) : Prop := ∀ b : Fin N, b.val % 2 = 1 → x b = dimerLabel b

instance (x : Cfg N) : Decidable (dimerLinks x) := by unfold dimerLinks; infer_instance

/-- The dimer gauge state in the `σˣ` basis (unnormalized amplitude `1` on allowed links). -/
def dimerGauge (x : Cfg N) : ℂ := if dimerLinks x then 1 else 0

/-- The full initial state: dimer gauge state ⊗ density-wave matter (unnormalized). -/
def dimerInit (N : ℕ) : Cfg N × Cfg N → ℂ :=
  fun p => (if p.1 = fockCfg (univ.filter (fun x : Fin N => Even x.val)) then 1 else 0) *
    dimerGauge p.2

/-- Dimer charge configurations: `q_{2k} = q_{2k+1}` for every pair inside the chain. -/
def IsDimer (q : Cfg N) : Prop :=
  ∀ k : Fin N, k.val % 2 = 0 → ∀ hk : k.val + 1 < N, q k = q ⟨k.val + 1, hk⟩

instance (q : Cfg N) : Decidable (IsDimer q) := by unfold IsDimer; infer_instance

/-! ### Auxiliary lemmas: configurations as sequences indexed by `ℕ` -/

/-- Extension of a configuration to a sequence indexed by `ℕ` (zero beyond the chain). -/
def dimerExt (x : Cfg N) (i : ℕ) : ZMod 2 := if h : i < N then x ⟨i, h⟩ else 0

/-- `dimerExt` agrees with the configuration on the chain. -/
theorem dimerExt_val (x : Cfg N) (b : Fin N) : dimerExt x b.val = x b := by
  simp [dimerExt, b.isLt]

/-- `dimerExt` agrees with the configuration on the chain (unbundled form). -/
theorem dimerExt_of_lt (x : Cfg N) (i : ℕ) (h : i < N) : dimerExt x i = x ⟨i, h⟩ := by
  simp [dimerExt, h]

/-- First link label: `x_{b_0} = q_0 + n_0`. -/
theorem dimerExt_linkOf_zero (n q : Cfg N) (h : 0 < N) :
    dimerExt (linkOf n q) 0 = dimerExt q 0 + dimerExt n 0 := by
  rw [dimerExt_of_lt _ _ h, dimerExt_of_lt _ _ h, dimerExt_of_lt _ _ h]
  unfold linkOf
  have hset : univ.filter (fun k : Fin N => k ≤ ⟨0, h⟩) = {⟨0, h⟩} := by
    ext k
    simp [Fin.le_def, Fin.ext_iff]
  rw [hset, sum_singleton]

/-- Recursion of the link labels: `x_{b_{i+1}} = x_{b_i} + q_{i+1} + n_{i+1}`. -/
theorem dimerExt_linkOf_succ (n q : Cfg N) (i : ℕ) (h : i + 1 < N) :
    dimerExt (linkOf n q) (i + 1) =
      dimerExt (linkOf n q) i + (dimerExt q (i + 1) + dimerExt n (i + 1)) := by
  rw [dimerExt_of_lt _ _ h, dimerExt_of_lt _ _ (by omega : i < N), dimerExt_of_lt _ _ h,
    dimerExt_of_lt _ _ h]
  unfold linkOf
  have hset : univ.filter (fun k : Fin N => k ≤ ⟨i + 1, h⟩) =
      insert ⟨i + 1, h⟩ (univ.filter (fun k : Fin N => k ≤ ⟨i, by omega⟩)) := by
    ext k
    simp only [mem_filter, mem_univ, true_and, mem_insert, Fin.le_def, Fin.ext_iff]
    omega
  rw [hset, sum_insert (by simp [Fin.le_def])]
  ring

/-- The density-wave configuration is occupied exactly on the even sites. -/
theorem dimerExt_cdw (i : ℕ) (h : i < N) :
    dimerExt (fockCfg (univ.filter (fun x : Fin N => Even x.val))) i = if Even i then 1 else 0 := by
  simp [dimerExt, h, fockCfg]

/-- Closed form of the link labels at the odd links for density-wave matter:
`x_{b_{2k+1}} = (k+1) + ∑_{i ≤ k} (q_{2i} + q_{2i+1})`. -/
theorem dimerExt_linkOf_cdw_odd (q : Cfg N) (k : ℕ) (h : 2 * k + 1 < N) :
    dimerExt (linkOf (fockCfg (univ.filter (fun x : Fin N => Even x.val))) q) (2 * k + 1) =
      ((k + 1 : ℕ) : ZMod 2) +
        ∑ i ∈ range (k + 1), (dimerExt q (2 * i) + dimerExt q (2 * i + 1)) := by
  induction k with
  | zero =>
    show dimerExt (linkOf _ q) (0 + 1) = _
    rw [dimerExt_linkOf_succ _ _ _ h, dimerExt_linkOf_zero _ _ (by omega),
      dimerExt_cdw _ (by omega), dimerExt_cdw _ h]
    simp
    ring
  | succ k ih =>
    have h' : 2 * k + 1 + 1 + 1 < N := by omega
    rw [show 2 * (k + 1) + 1 = (2 * k + 1 + 1) + 1 by ring,
      dimerExt_linkOf_succ _ _ (2 * k + 1 + 1) h', dimerExt_linkOf_succ _ _ (2 * k + 1) (by omega),
      ih (by omega), sum_range_succ _ (k + 1), dimerExt_cdw _ h', dimerExt_cdw _ (by omega)]
    have he1 : Even (2 * k + 1 + 1) := ⟨k + 1, by ring⟩
    have he2 : ¬ Even (2 * k + 1 + 1 + 1) := by
      rw [Nat.not_even_iff_odd]; exact ⟨k + 1, by ring⟩
    rw [ite_eq_left he1, ite_eq_right he2, show 2 * (k + 1) = 2 * k + 1 + 1 by ring]
    push_cast
    ring

/-- In `ZMod 2`, `a + b = 0` iff `a = b`. -/
theorem dimer_zmod2_add_eq_zero (a b : ZMod 2) : a + b = 0 ↔ a = b := by
  revert a b; decide

/-- Partial sums vanish up to a downward-closed bound iff all terms vanish. -/
theorem dimer_partial_sums (a : ℕ → ZMod 2) :
    (∀ k, 2 * k + 1 < N → ∑ i ∈ range (k + 1), a i = 0) ↔ (∀ k, 2 * k + 1 < N → a k = 0) := by
  constructor
  · intro hs k hk
    cases k with
    | zero => simpa using hs 0 hk
    | succ k =>
      have h1 := hs (k + 1) hk
      rw [sum_range_succ, hs k (by omega), zero_add] at h1
      exact h1
  · intro ha k hk
    refine sum_eq_zero fun i hi => ha i ?_
    rw [mem_range] at hi
    omega

/-- The dimer link constraints for density-wave matter, as conditions on partial sums. -/
theorem dimerLinks_linkOf_cdw_iff (q : Cfg N) :
    dimerLinks (linkOf (fockCfg (univ.filter (fun x : Fin N => Even x.val))) q) ↔
      ∀ k, 2 * k + 1 < N →
        ∑ i ∈ range (k + 1), (dimerExt q (2 * i) + dimerExt q (2 * i + 1)) = 0 := by
  constructor
  · intro hd k hk
    have h1 := hd ⟨2 * k + 1, hk⟩ (by show (2 * k + 1) % 2 = 1; omega)
    rw [← dimerExt_of_lt (linkOf (fockCfg (univ.filter (fun x : Fin N => Even x.val))) q)
      (2 * k + 1) hk, dimerExt_linkOf_cdw_odd q k hk] at h1
    have h2 : dimerLabel (⟨2 * k + 1, hk⟩ : Fin N) = ((k + 1 : ℕ) : ZMod 2) := by
      unfold dimerLabel
      congr 2
      show (2 * k + 1) / 2 = k
      omega
    rw [h2] at h1
    exact (add_eq_left).mp h1
  · intro hs b hb
    have hb' : b.val = 2 * (b.val / 2) + 1 := by omega
    have hlt : 2 * (b.val / 2) + 1 < N := by omega
    rw [← dimerExt_val (linkOf (fockCfg (univ.filter (fun x : Fin N => Even x.val))) q) b, hb',
      dimerExt_linkOf_cdw_odd q _ hlt, hs _ hlt, add_zero]
    unfold dimerLabel
    rfl

/-- Dimer charge configurations, as conditions on the pairs. -/
theorem isDimer_iff_pairs (q : Cfg N) :
    IsDimer q ↔ ∀ k, 2 * k + 1 < N → dimerExt q (2 * k) + dimerExt q (2 * k + 1) = 0 := by
  constructor
  · intro hd k hk
    rw [dimer_zmod2_add_eq_zero, dimerExt_of_lt _ _ (by omega), dimerExt_of_lt _ _ hk]
    exact hd ⟨2 * k, by omega⟩ (by simp) hk
  · intro hs k hk hk1
    have hk' : k.val = 2 * (k.val / 2) := by omega
    have h1 := (dimer_zmod2_add_eq_zero _ _).mp (hs (k.val / 2) (by omega))
    rw [← hk', dimerExt_val, dimerExt_of_lt _ _ hk1] at h1
    exact h1

/-- For density-wave matter, the dimer link constraints are exactly the dimer charge
constraints. -/
theorem dimerLinks_linkOf_cdw (q : Cfg N) :
    dimerLinks (linkOf (fockCfg (univ.filter (fun x : Fin N => Even x.val))) q) ↔ IsDimer q := by
  rw [dimerLinks_linkOf_cdw_iff, isDimer_iff_pairs]
  exact dimer_partial_sums (fun i => dimerExt q (2 * i) + dimerExt q (2 * i + 1))

/-- **The dimer sector law**: in the dual frame the dimer initial state is the uniform
superposition of density-wave matter over all dimer charge configurations. -/
theorem dimerInit_dual (hN : Even N) (n q : Cfg N) :
    dimerInit N (dualEquiv (n, q)) =
      if n = fockCfg (univ.filter (fun x : Fin N => Even x.val)) ∧ IsDimer q then 1 else 0 := by
  have _hN := hN  -- the sector law holds for every `N`; evenness is not needed here
  simp only [dimerInit, dualEquiv_apply]
  by_cases hn : n = fockCfg (univ.filter (fun x : Fin N => Even x.val))
  · subst hn
    simp only [ite_true, true_and, one_mul, dimerGauge]
    exact if_congr (dimerLinks_linkOf_cdw q) rfl rfl
  · simp [hn]

/-- The number of dimer charge configurations is `2^{N/2}`. -/
theorem card_isDimer (hN : Even N) :
    (univ.filter (fun q : Cfg N => IsDimer q)).card = 2 ^ (N / 2) := by
  obtain ⟨m, hm⟩ := hN
  have hcard : (univ : Finset (Fin (N / 2) → ZMod 2)).card = 2 ^ (N / 2) := by
    rw [card_univ, Fintype.card_fun, ZMod.card, Fintype.card_fin]
  rw [← hcard]
  refine card_nbij' (fun q k => q ⟨2 * k.val, by omega⟩)
    (fun r x => r ⟨x.val / 2, by omega⟩) ?_ ?_ ?_ ?_
  · intro q _
    exact mem_coe.mpr (mem_univ _)
  · intro r _
    rw [mem_coe, mem_filter]
    refine ⟨mem_univ _, ?_⟩
    intro k hk hk1
    show r ⟨k.val / 2, _⟩ = r ⟨(k.val + 1) / 2, _⟩
    congr 2
    omega
  · intro q hq
    rw [mem_coe, mem_filter] at hq
    have hq' : IsDimer q := hq.2
    funext x
    show q ⟨2 * (x.val / 2), _⟩ = q x
    by_cases hx : x.val % 2 = 0
    · congr 1
      exact Fin.ext (by show 2 * (x.val / 2) = x.val; omega)
    · have h1 := hq' ⟨2 * (x.val / 2), by omega⟩ (by show 2 * (x.val / 2) % 2 = 0; omega)
        (by show 2 * (x.val / 2) + 1 < N; omega)
      rw [h1]
      congr 1
      exact Fin.ext (by show 2 * (x.val / 2) + 1 = x.val; omega)
  · intro r _
    funext k
    show r ⟨2 * k.val / 2, _⟩ = r k
    congr 1
    exact Fin.ext (by show 2 * k.val / 2 = k.val; omega)

/-! ### Auxiliary lemmas: the dual-frame expectation formula -/

/-- Expectations are invariant under relabelling the basis by a bijection. -/
theorem dimer_expect_submatrix {κ ι : Type*} [Fintype κ] [Fintype ι] (M : Matrix ι ι ℂ)
    (e : κ ≃ ι) (v : ι → ℂ) : expect (M.submatrix e e) (v ∘ e) = expect M v := by
  unfold expect
  rw [submatrix_mulVec_equiv, Function.comp_assoc, Equiv.self_comp_symm, Function.comp_id]
  exact comp_equiv_dotProduct_comp_equiv (star v) (M *ᵥ v) e

/-- Matrix–vector product with a block-diagonal matrix, block by block. -/
theorem dimer_blockDiagonal_mulVec {m o : Type*} [Fintype m] [Fintype o] [DecidableEq o]
    (M : o → Matrix m m ℂ) (v : m × o → ℂ) (i : m) (k : o) :
    (blockDiagonal M *ᵥ v) (i, k) = (M k *ᵥ fun j => v (j, k)) i := by
  simp only [mulVec, dotProduct, blockDiagonal_apply', Fintype.sum_prod_type]
  refine Finset.sum_congr rfl fun j _ => ?_
  simp [ite_mul]

/-- Expectation of a block-diagonal observable is the sum of the blockwise expectations. -/
theorem dimer_expect_blockDiagonal {m o : Type*} [Fintype m] [Fintype o] [DecidableEq o]
    (M : o → Matrix m m ℂ) (v : m × o → ℂ) :
    expect (blockDiagonal M) v = ∑ k, expect (M k) (fun j => v (j, k)) := by
  unfold expect
  simp only [dotProduct, Fintype.sum_prod_type, dimer_blockDiagonal_mulVec, Pi.star_apply]
  exact Finset.sum_comm

/-- Dual-frame expectation formula: for every initial state `Ψ₀` and gauge-invariant `O`, the
expectation in the evolved LGT state is a sum over charge sectors of free-fermion expectations. -/
theorem dimer_expect_evo_dual (J h : ℝ) (Ψ₀ : Cfg N × Cfg N → ℂ)
    (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (H J h) t *ᵥ Ψ₀) =
      ∑ q : Cfg N, expect (sectorBlock O q)
        (evo (Hsec J h q) t *ᵥ fun n => Ψ₀ (dualEquiv (n, q))) := by
  have hU : blockDiagonal (fun q : Cfg N => evo (Hsec J h q) t) =
      (evo (H (N := N) J h) t).submatrix dualEquiv dualEquiv := by
    rw [evo_submatrix, duality, ← evo_blockDiagonal]
    rfl
  have hstate : (evo (H J h) t *ᵥ Ψ₀) ∘ dualEquiv =
      blockDiagonal (fun q : Cfg N => evo (Hsec J h q) t) *ᵥ (Ψ₀ ∘ dualEquiv) := by
    rw [hU, submatrix_mulVec_equiv, Function.comp_assoc, Equiv.self_comp_symm, Function.comp_id]
  rw [← dimer_expect_submatrix O dualEquiv, gaugeInv_submatrix O hO, hstate,
    dimer_expect_blockDiagonal]
  refine Finset.sum_congr rfl fun q _ => ?_
  congr 1
  ext n
  rw [dimer_blockDiagonal_mulVec]
  rfl

/-- **Sector averaging over the dimer law**: every gauge-invariant observable evolves as the
uniform average over dimer charge configurations of free-fermion expectations. -/
theorem sector_average_dimer (hN : Even N) (J h : ℝ)
    (O : Matrix (Cfg N × Cfg N) (Cfg N × Cfg N) ℂ) (hO : GaugeInv O) (t : ℝ) :
    expect O (evo (H J h) t *ᵥ dimerInit N) =
      ∑ q ∈ univ.filter (fun q : Cfg N => IsDimer q),
        expect (sectorBlock O q)
          (evo (Hsec J h q) t *ᵥ Pi.single (fockCfg (univ.filter (fun x : Fin N => Even x.val))) 1)
          := by
  rw [dimer_expect_evo_dual J h _ O hO t, Finset.sum_filter]
  refine Finset.sum_congr rfl fun q _ => ?_
  by_cases hq : IsDimer q
  · rw [ite_eq_left hq]
    congr 2
    ext n
    rw [dimerInit_dual hN, Pi.single_apply]
    simp [hq]
  · rw [ite_eq_right hq]
    have hvec : (fun n => dimerInit N (dualEquiv (n, q))) = 0 := by
      ext n
      rw [dimerInit_dual hN]
      simp [hq]
    rw [hvec, mulVec_zero]
    simp [expect]

end DFL
