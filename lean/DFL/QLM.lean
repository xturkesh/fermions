import DFL.Entropy
import DFL.Basic

/-!
# Number entropy in U(1) quantum link models: kinematic bound (paper: Proposition `prop:gausscap`)

In a U(1) quantum link model with spin-`S` links, in a fixed sector of background charges, the
summed Gauss law expresses the matter particle number of an interval through the electric fields
on its two boundary links. Hence, for *any* state of the sector (and therefore for any dynamics),
the particle-number distribution of an interval is supported on at most `4S + 1` values, and its
Shannon entropy is at most `ln(4S + 1)`.

Conventions: staggered fermions with charge `ρ_j = n_j - [j odd]`; electric fields
`E_b ∈ {-S, …, S}` on the links `b_j = (j, j+1)`, stored as `E_b + S ∈ Fin (2S+1)`; a fixed
boundary field `E_{-1}`; integer background charges `Q_j`; Gauss law
`E_{b_j} - E_{b_{j-1}} = ρ_j + Q_j`.
-/

open Finset

noncomputable section

namespace DFL

variable {N : ℕ}

/-- Electric field `E ∈ {-S, …, S}` encoded by `e ∈ Fin (2S+1)`. -/
def efield (S : ℕ) (e : Fin (2 * S + 1)) : ℤ := (e.val : ℤ) - S

/-- Electric field on the link to the left of site `j` (the boundary field `E0` for `j = 0`). -/
def prevField (S : ℕ) (E0 : ℤ) (e : Fin N → Fin (2 * S + 1)) (j : Fin N) : ℤ :=
  if hj : j.val = 0 then E0 else efield S (e ⟨j.val - 1, by have := j.isLt; omega⟩)

/-- Staggered matter charge `ρ_j = n_j - [j odd]`. -/
def stagCharge (n : Cfg N) (j : Fin N) : ℤ :=
  (if n j = 0 then 0 else 1) - (if j.val % 2 = 1 then 1 else 0)

/-- The Gauss law `E_{b_j} - E_{b_{j-1}} = ρ_j + Q_j` for all sites. -/
def GaussLaw (S : ℕ) (E0 : ℤ) (Q : Fin N → ℤ) (n : Cfg N) (e : Fin N → Fin (2 * S + 1)) :
    Prop :=
  ∀ j : Fin N, efield S (e j) - prevField S E0 e j = stagCharge n j + Q j

/-- Number of particles of `n` in the interval `[a, b]`, as an integer. -/
def countInt (a b : ℕ) (n : Cfg N) : ℤ :=
  ∑ j ∈ univ.filter (fun j : Fin N => a ≤ j.val ∧ j.val ≤ b), (if n j = 0 then 0 else 1)

/-- The one-site interval `[a, a]` is `{a}`. -/
private lemma filter_Icc_self (a : Fin N) :
    univ.filter (fun j : Fin N => a.val ≤ j.val ∧ j.val ≤ a.val) = {a} := by
  ext j
  simp only [mem_filter, mem_univ, true_and, mem_singleton, Fin.ext_iff]
  omega

/-- Peeling off the top site `B + 1` of the interval `[a, B + 1]`. -/
private lemma sum_filter_Icc_succ (a B : ℕ) (hB : B + 1 < N) (f : Fin N → ℤ)
    (haB : a ≤ B + 1) :
    ∑ j ∈ univ.filter (fun j : Fin N => a ≤ j.val ∧ j.val ≤ B + 1), f j =
      ∑ j ∈ univ.filter (fun j : Fin N => a ≤ j.val ∧ j.val ≤ B), f j +
        f ⟨B + 1, hB⟩ := by
  have hset : univ.filter (fun j : Fin N => a ≤ j.val ∧ j.val ≤ B + 1) =
      insert ⟨B + 1, hB⟩ (univ.filter (fun j : Fin N => a ≤ j.val ∧ j.val ≤ B)) := by
    ext j
    simp only [mem_filter, mem_univ, true_and, mem_insert, Fin.ext_iff]
    omega
  rw [hset, sum_insert (by simp), add_comm]

/-- The electric field lies in `[-S, S]`. -/
lemma efield_bounds (S : ℕ) (e : Fin (2 * S + 1)) :
    -(S : ℤ) ≤ efield S e ∧ efield S e ≤ S := by
  unfold efield
  have := e.isLt
  omega

/-- Away from the left edge, the field on the link to the left of `j` lies in `[-S, S]`. -/
lemma prevField_bounds (S : ℕ) (E0 : ℤ) (e : Fin N → Fin (2 * S + 1)) (j : Fin N)
    (hj : j.val ≠ 0) : -(S : ℤ) ≤ prevField S E0 e j ∧ prevField S E0 e j ≤ S := by
  unfold prevField
  rw [dite_eq_right hj]
  exact efield_bounds S _

/-- **Summed Gauss law**: the particle number of `[a, b]` is fixed by the boundary fields. -/
theorem countInt_gauss (S : ℕ) (E0 : ℤ) (Q : Fin N → ℤ) (n : Cfg N)
    (e : Fin N → Fin (2 * S + 1)) (hG : GaussLaw S E0 Q n e) (a b : Fin N) (hab : a ≤ b) :
    countInt a.val b.val n =
      efield S (e b) - prevField S E0 e a +
        ∑ j ∈ univ.filter (fun j : Fin N => a.val ≤ j.val ∧ j.val ≤ b.val),
          ((if j.val % 2 = 1 then 1 else 0) - Q j) := by
  have key : ∀ B : ℕ, a.val ≤ B → ∀ hB : B < N, countInt a.val B n =
      efield S (e ⟨B, hB⟩) - prevField S E0 e a +
        ∑ j ∈ univ.filter (fun j : Fin N => a.val ≤ j.val ∧ j.val ≤ B),
          ((if j.val % 2 = 1 then 1 else 0) - Q j) := by
    intro B haB
    induction B, haB using Nat.le_induction with
    | base =>
      intro hB
      have ha : (⟨a.val, hB⟩ : Fin N) = a := rfl
      rw [ha]
      unfold countInt
      rw [filter_Icc_self a, sum_singleton, sum_singleton]
      have h := hG a
      unfold stagCharge at h
      linarith
    | succ B haB ih =>
      intro hB
      have hB' : B < N := by omega
      have ih' := ih hB'
      unfold countInt at ih' ⊢
      rw [sum_filter_Icc_succ a.val B hB _ (by omega),
        sum_filter_Icc_succ a.val B hB _ (by omega), ih']
      have hprev : prevField S E0 e ⟨B + 1, hB⟩ = efield S (e ⟨B, hB'⟩) := by
        unfold prevField
        rw [dite_eq_right (by simp)]
        rfl
      have h := hG ⟨B + 1, hB⟩
      rw [hprev] at h
      unfold stagCharge at h
      linarith
  exact key b.val hab b.isLt

/-- Distribution of the particle number of `[a, b]` in a state `v` of the quantum link model. -/
def qlmNumDist (S : ℕ) (a b : ℕ) (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ) (k : ℕ) : ℝ :=
  ∑ c ∈ univ.filter (fun c : Cfg N × (Fin N → Fin (2 * S + 1)) => countInt a b c.1 = k),
    Complex.normSq (v c)

/-- The particle number of an interval is nonnegative. -/
lemma countInt_nonneg (a b : ℕ) (n : Cfg N) : 0 ≤ countInt a b n :=
  Finset.sum_nonneg fun j _ => by split_ifs <;> norm_num

/-- The particle number of an interval is at most `N`. -/
lemma countInt_le (a b : ℕ) (n : Cfg N) : countInt a b n ≤ N := by
  unfold countInt
  set A := univ.filter (fun j : Fin N => a ≤ j.val ∧ j.val ≤ b)
  calc ∑ j ∈ A, (if n j = 0 then (0 : ℤ) else 1)
      ≤ ∑ j ∈ A, (1 : ℤ) := Finset.sum_le_sum fun j _ => by split_ifs <;> norm_num
    _ = (A.card : ℤ) := by simp
    _ ≤ N := by
        have := Finset.card_le_univ A
        rw [Fintype.card_fin] at this
        exact_mod_cast this

/-- The particle-number distribution of a normalized state is a probability distribution. -/
lemma qlmNumDist_isProb (S a b : ℕ) (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ)
    (hnorm : ∑ c, Complex.normSq (v c) = 1) : IsProb N (qlmNumDist S a b v) where
  nonneg k := Finset.sum_nonneg fun c _ => Complex.normSq_nonneg _
  sum_one := by
    unfold qlmNumDist
    rw [← hnorm]
    have hfilt : ∀ k : ℕ,
        univ.filter (fun c : Cfg N × (Fin N → Fin (2 * S + 1)) => countInt a b c.1 = k) =
          univ.filter (fun c : Cfg N × (Fin N → Fin (2 * S + 1)) =>
            (countInt a b c.1).toNat = k) := by
      intro k
      apply Finset.filter_congr
      intro c _
      have := countInt_nonneg a b c.1
      omega
    simp_rw [hfilt]
    apply Finset.sum_fiberwise_of_maps_to
    intro c _
    rw [mem_range]
    have := countInt_le a b c.1
    omega

/-- If the particle number of every configuration in the support of `v` lies in the window
`[L, L + W]`, the number entropy is at most `ln(W + 1)`. -/
lemma shannon_qlmNumDist_le_of_window (S a b : ℕ)
    (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ)
    (hnorm : ∑ c, Complex.normSq (v c) = 1) (L : ℤ) (W : ℕ)
    (hwin : ∀ c, v c ≠ 0 → L ≤ countInt a b c.1 ∧ countInt a b c.1 ≤ L + W) :
    shannon N (qlmNumDist S a b v) ≤ Real.log ((W : ℝ) + 1) := by
  set s : Finset ℕ := (Finset.Icc L (L + W)).image Int.toNat
  have hs : s.Nonempty := ⟨L.toNat, Finset.mem_image.mpr ⟨L, by simp, rfl⟩⟩
  have hsupp : ∀ k ∈ range (N + 1), k ∉ s → qlmNumDist S a b v k = 0 := by
    intro k _ hks
    unfold qlmNumDist
    apply Finset.sum_eq_zero
    intro c hc
    rw [Finset.mem_filter] at hc
    by_contra hne
    have hv : v c ≠ 0 := fun h => hne (by rw [h]; simp)
    obtain ⟨h1, h2⟩ := hwin c hv
    rw [hc.2] at h1 h2
    exact hks (Finset.mem_image.mpr ⟨(k : ℤ), Finset.mem_Icc.mpr ⟨h1, h2⟩, by simp⟩)
  have hcard : s.card ≤ W + 1 := by
    calc s.card ≤ (Finset.Icc L (L + W)).card := Finset.card_image_le
      _ = W + 1 := by rw [Int.card_Icc]; omega
  calc shannon N (qlmNumDist S a b v) ≤ Real.log s.card :=
        shannon_le_log_card (qlmNumDist_isProb S a b v hnorm) s hs hsupp
    _ ≤ Real.log ((W : ℝ) + 1) := by
        apply Real.log_le_log (by exact_mod_cast hs.card_pos)
        exact_mod_cast hcard

/-- **Kinematic cap on the number entropy** (paper: Proposition `prop:gausscap`): for every
normalized state supported on a fixed background-charge sector, the number entropy of an interval
is at most `ln(4S + 1)`. -/
theorem qlm_numberEntropy_le (S : ℕ) (E0 : ℤ) (Q : Fin N → ℤ)
    (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ)
    (hnorm : ∑ c, Complex.normSq (v c) = 1)
    (hsupp : ∀ c, v c ≠ 0 → GaussLaw S E0 Q c.1 c.2) (a b : Fin N) (hab : a ≤ b) :
    shannon N (qlmNumDist S a.val b.val v) ≤ Real.log (4 * S + 1) := by
  set C : ℤ := ∑ j ∈ univ.filter (fun j : Fin N => a.val ≤ j.val ∧ j.val ≤ b.val),
    ((if j.val % 2 = 1 then 1 else 0) - Q j) with hC
  have hS : (0 : ℝ) ≤ S := Nat.cast_nonneg S
  by_cases ha : a.val = 0
  · -- left edge: the left boundary field is the fixed `E0`
    have h := shannon_qlmNumDist_le_of_window S a.val b.val v hnorm (C - E0 - S) (2 * S)
      (fun c hc => by
        have hG := countInt_gauss S E0 Q c.1 c.2 (hsupp c hc) a b hab
        have hp : prevField S E0 c.2 a = E0 := by
          unfold prevField
          rw [dite_eq_left ha]
        have he := efield_bounds S (c.2 b)
        rw [hG, hp, ← hC]
        push_cast
        constructor <;> linarith [he.1, he.2])
    refine h.trans (Real.log_le_log (by positivity) ?_)
    push_cast
    linarith
  · -- bulk: both boundary fields lie in `[-S, S]`
    have h := shannon_qlmNumDist_le_of_window S a.val b.val v hnorm (C - 2 * S) (4 * S)
      (fun c hc => by
        have hG := countInt_gauss S E0 Q c.1 c.2 (hsupp c hc) a b hab
        have hp := prevField_bounds S E0 c.2 a ha
        have he := efield_bounds S (c.2 b)
        rw [hG, ← hC]
        push_cast
        constructor <;> linarith [he.1, he.2, hp.1, hp.2])
    push_cast at h
    exact h

/-- At a fixed-flux edge (`a = 0`) the cap improves to `ln(2S + 1)`. -/
theorem qlm_numberEntropy_le_edge (S : ℕ) (E0 : ℤ) (Q : Fin N → ℤ)
    (v : Cfg N × (Fin N → Fin (2 * S + 1)) → ℂ)
    (hnorm : ∑ c, Complex.normSq (v c) = 1)
    (hsupp : ∀ c, v c ≠ 0 → GaussLaw S E0 Q c.1 c.2) (b : Fin N) :
    shannon N (qlmNumDist S 0 b.val v) ≤ Real.log (2 * S + 1) := by
  have hN : 0 < N := Nat.lt_of_le_of_lt (Nat.zero_le _) b.isLt
  let a : Fin N := ⟨0, hN⟩
  have hab : a ≤ b := Nat.zero_le _
  set C : ℤ := ∑ j ∈ univ.filter (fun j : Fin N => a.val ≤ j.val ∧ j.val ≤ b.val),
    ((if j.val % 2 = 1 then 1 else 0) - Q j) with hC
  have h := shannon_qlmNumDist_le_of_window S a.val b.val v hnorm (C - E0 - S) (2 * S)
    (fun c hc => by
      have hG := countInt_gauss S E0 Q c.1 c.2 (hsupp c hc) a b hab
      have hp : prevField S E0 c.2 a = E0 := by
        unfold prevField
        rw [dite_eq_left rfl]
      have he := efield_bounds S (c.2 b)
      rw [hG, hp, ← hC]
      push_cast
      constructor <;> linarith [he.1, he.2])
  push_cast at h
  exact h

end DFL
