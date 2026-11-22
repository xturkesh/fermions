import DFL.DimerLower
import DFL.Dimer
import DFL.OneBody

/-!
# Delocalization for the dimer sector law (paper: Theorem `thm:dimer`(a), and (d) for `p > 1`)

Let `0 < 2|h| < J` and `λ = 2|h|/J ∈ (0, 1)`. In every dimer sector `q` (`IsDimer q`: `q_{2k} = q_{2k+1}`) the
one-body operator is `h_q = J(-Δ + (2h/J) diag(q))`, a random-dimer operator with values `±λ`.

* (a) **Every sector, every volume**: for `p > 0` and `1 ≤ JT ≤ min(x, N-1-x)` (0-indexed sites),
  `∫_0^∞ (dt/T) e^{-t/T} ∑_y |y - x|^p |⟨δ_y, e^{-ith_q} δ_x⟩|² ≥ c_{λ,p} (JT)^{p-1}`
  (`dimer_sector_moment_lower`), by Proposition `prop:dimerlower` after the time change `s = Jt`.
* (d) for `p > 1`: **failure of uniform localization**, sector by sector: for every choice of dimer sectors
  `q^{(N)}`, `sup_{N, x, t} ∑_y |y - x|^p |⟨δ_y, e^{-ith_{q^{(N)}}} δ_x⟩|² = ∞` (`dimer_sup_moment_unbounded`);
  consequently there is no summable profile `F` with `∑_d (d+1)^p F(d) < ∞` such that
  `E_dim Q^{[1,N]}_q(x, y) ≤ F(|x - y|)` for all even `N` (`not_EC_dimer`), i.e. Hypothesis (EC) fails for the
  dimer law. (For `1/2 < p ≤ 1` the paper uses the transport bounds of Jitomirskaya–Schulz-Baldes–Stolz,
  which are not formalized.)
-/

open Matrix Finset MeasureTheory

noncomputable section

namespace DFL

variable {N : ℕ}

/-- The uniform average over the dimer charge configurations (the dimer sector law `μ_dim`). -/
def dimerAvg (f : Cfg N → ℝ) : ℝ :=
  (∑ q ∈ univ.filter (fun q : Cfg N => IsDimer q), f q) /
    ((univ.filter (fun q : Cfg N => IsDimer q)).card : ℝ)

namespace DimerLGTAux

/-- The relabelling `Fin N ≃ [0, N-1] ⊆ ℤ`, `y ↦ y`. -/
def finIcc (N : ℕ) : Fin N ≃ Set.Icc (0 : ℤ) ((N : ℤ) - 1) where
  toFun y := ⟨(y.val : ℤ), by have := y.isLt; constructor <;> omega⟩
  invFun z := ⟨z.1.toNat, by have := z.2.1; have := z.2.2; omega⟩
  left_inv y := by ext; simp
  right_inv z := by ext; have := z.2.1; simp; omega

/-- The relabelling is the identity on integers. -/
@[simp] lemma finIcc_apply (y : Fin N) :
    ((finIcc N y : Set.Icc (0 : ℤ) ((N : ℤ) - 1)) : ℤ) = y.val := rfl

/-- The potential of the sector `q`: `v(y) = (2h/J) sgn(q_y)` on `[0, N-1]`. -/
def potQ (J h : ℝ) (q : Cfg N) (y : ℤ) : ℝ := 2 * h / J * sgn (dimerExt q y.toNat)

/-- `h_q = J (-Δ + diag(v))` after relabelling `Fin N ≃ [0, N-1]`. -/
lemma hq_eq_submatrix (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) :
    hq J h q = (J • schrH 0 ((N : ℤ) - 1) (potQ J h q)).submatrix (finIcc N) (finIcc N) := by
  ext y z
  simp only [hq, schrH, pathAdj, potQ, Matrix.add_apply, Matrix.smul_apply, Matrix.submatrix_apply,
    diagonal_apply, smul_eq_mul, finIcc_apply]
  by_cases hyz : y = z
  · subst hyz
    simp [dimerExt_val]
    field_simp
  · have hyz' : (y.val : ℤ) ≠ z.val := by
      intro h; exact hyz (Fin.ext (by exact_mod_cast h))
    simp only [hyz, hyz', ite_false, mul_zero, add_zero]
    by_cases hadj : y.val + 1 = z.val ∨ z.val + 1 = y.val
    · have h2 : |(y.val : ℤ) - z.val| = 1 := by rw [abs_eq zero_le_one]; omega
      simp only [hadj, h2, ite_true]
      ring
    · have h2 : ¬ |(y.val : ℤ) - z.val| = 1 := by rw [abs_eq zero_le_one]; omega
      simp only [hadj, h2, ite_false]
      ring

/-- `e^{-it(J B)} = e^{-i(Jt)B}`. -/
lemma evo_smul {ι : Type*} [Fintype ι] [DecidableEq ι] (B : Matrix ι ι ℝ) (J t : ℝ) :
    evo (J • B) t = evo B (J * t) := by
  unfold evo
  congr 1
  ext i j
  simp only [cplx, Matrix.smul_apply, Matrix.map_apply, smul_eq_mul, Complex.ofReal_mul]
  ring

/-- The moment of `h_q` at time `t` is the moment of `-Δ + diag(v)` at time `Jt`. -/
lemma moment_eq (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) (x : Fin N) (p t : ℝ) :
    ∑ y : Fin N, |((y.val : ℝ) - x.val)| ^ p * ‖evo (hq J h q) t y x‖ ^ 2 =
      ∑ y : Set.Icc (0 : ℤ) ((N : ℤ) - 1), |((y.1 - (finIcc N x).1 : ℤ) : ℝ)| ^ p *
        ‖evo (schrH 0 ((N : ℤ) - 1) (potQ J h q)) (J * t) y (finIcc N x)‖ ^ 2 := by
  rw [hq_eq_submatrix J h hJ q, ← evo_submatrix, evo_smul, ← (finIcc N).sum_comp]
  refine Finset.sum_congr rfl fun y _ => ?_
  simp only [Matrix.submatrix_apply, finIcc_apply]
  push_cast
  rfl

/-- Time change `s = Jt` in the time average. -/
lemma integral_time_change (g : ℝ → ℝ) (J T : ℝ) (hJ : 0 < J) :
    ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) * g (J * t) =
      ∫ s in Set.Ioi (0 : ℝ), (1 / (J * T)) * Real.exp (-s / (J * T)) * g s := by
  have h := integral_comp_mul_left_Ioi'
    (fun s => (1 / (J * T)) * Real.exp (-s / (J * T)) * g s) 0 hJ
  rw [mul_zero] at h
  rw [← h, smul_eq_mul, ← integral_const_mul]
  congr 1
  funext t
  have e1 : -(J * t) / (J * T) = -t / T := by
    rw [neg_div, neg_div, mul_div_mul_left _ _ hJ.ne']
  have e2 : J * (1 / (J * T)) = 1 / T := by
    rw [one_div, mul_inv, ← mul_assoc, mul_inv_cancel₀ hJ.ne', one_mul, one_div]
  rw [e1, ← mul_assoc, ← mul_assoc, e2]

/-- The time-averaged moment of `h_q` is the moment `M_p(JT; x)` of `-Δ + diag(v)`. -/
lemma integral_moment_eq (J h : ℝ) (hJ : 0 < J) (q : Cfg N) (x : Fin N) (p T : ℝ) :
    ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) *
        ∑ y : Fin N, |((y.val : ℝ) - x.val)| ^ p * ‖evo (hq J h q) t y x‖ ^ 2 =
      dimerMom 0 ((N : ℤ) - 1) (potQ J h q) p (J * T) (finIcc N x) := by
  simp_rw [moment_eq J h hJ.ne' q x p]
  exact integral_time_change (fun s => ∑ y : Set.Icc (0 : ℤ) ((N : ℤ) - 1),
    |((y.1 - (finIcc N x).1 : ℤ) : ℝ)| ^ p *
      ‖evo (schrH 0 ((N : ℤ) - 1) (potQ J h q)) s y (finIcc N x)‖ ^ 2) J T hJ

/-- In a dimer sector, `v = (2h/J) sgn(q)` is a dimer configuration with values `±2|h|/J` and
alignment `0`. -/
lemma potQ_isDimerCfg (J h : ℝ) (q : Cfg N) (hq : IsDimer q) :
    IsDimerCfg (2 * |h| / J) 0 (potQ J h q) 0 ((N : ℤ) - 1 - 0 + 1).toNat := by
  have hL : (((N : ℤ) - 1 - 0 + 1).toNat : ℤ) = N := by omega
  refine ⟨fun y _ _ => ?_, fun k hk1 hk2 => ?_⟩
  · simp only [potQ, sgn]
    rcases abs_choice h with h' | h'
    · rw [h']
      split_ifs
      · left; ring
      · right; ring
    · rw [h']
      split_ifs
      · right; ring
      · left; ring
  · rw [hL] at hk2
    simp only [potQ]
    have e1 : (2 * k - 0).toNat = 2 * k.toNat := by omega
    have e2 : (2 * k - 0 + 1).toNat = 2 * k.toNat + 1 := by omega
    rw [e1, e2, dimerExt_of_lt _ _ (by omega), dimerExt_of_lt _ _ (by omega)]
    rw [hq ⟨2 * k.toNat, by omega⟩ (by simp) (by simp; omega)]

/-- `c_{λ,p} > 0`. -/
lemma dimerClp_pos (lam p : ℝ) (h0 : 0 < lam) (h1 : lam < 1) (hp : 0 < p) :
    0 < dimerClp lam p := by
  have hκ : 0 < dimerKappa lam := by linarith [DimerTransferAux.kappa_ge_one lam h0 h1]
  unfold dimerClp dimerC1
  positivity

/-- For `p > 1` and `c > 0`, `c s^{p-1}` exceeds any bound for some `s ≥ 1`. -/
lemma exists_large (p : ℝ) (hp : 1 < p) (c : ℝ) (hc : 0 < c) (C : ℝ) :
    ∃ s : ℝ, 1 ≤ s ∧ C < c * s ^ (p - 1) := by
  have h := (tendsto_rpow_atTop (by linarith : 0 < p - 1)).const_mul_atTop hc
  exact ((Filter.eventually_ge_atTop 1).and (h.eventually_gt_atTop C)).exists

/-- `∫_0^∞ (1/T) e^{-t/T} dt = 1`. -/
lemma integral_weight (T : ℝ) (hT : 0 < T) :
    ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) = 1 := by
  have h := integral_comp_mul_left_Ioi (fun s => Real.exp (-s)) 0 (inv_pos.mpr hT)
  simp only [mul_zero, integral_exp_neg_Ioi_zero, inv_inv, smul_eq_mul, mul_one] at h
  rw [integral_const_mul]
  have e : ∀ t : ℝ, -t / T = -(T⁻¹ * t) := fun t => by ring
  simp_rw [e, h]
  field_simp

/-- The weight `(1/T) e^{-t/T}` is integrable on `(0, ∞)`. -/
lemma integrableOn_weight (T : ℝ) (hT : 0 < T) :
    IntegrableOn (fun t => (1 / T) * Real.exp (-t / T)) (Set.Ioi 0) := by
  have h : IntegrableOn (fun t : ℝ => (1 / T) * Real.exp (-T⁻¹ * t)) (Set.Ioi 0) :=
    (exp_neg_integrableOn_Ioi 0 (inv_pos.mpr hT)).const_mul (1 / T)
  refine h.congr_fun (fun t _ => ?_) measurableSet_Ioi
  simp only
  congr 2
  ring

/-- A time average of a function with values in `[0, C]` is at most `C`. -/
lemma integral_weight_le (T C : ℝ) (hT : 0 < T) (D : ℝ → ℝ) (hD0 : ∀ t, 0 ≤ D t)
    (hDC : ∀ t, D t ≤ C) :
    ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) * D t ≤ C := by
  calc ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) * D t
      ≤ ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) * C := by
        refine integral_mono_of_nonneg (ae_of_all _ fun t => ?_)
          ((integrableOn_weight T hT).mul_const C) (ae_of_all _ fun t => ?_)
        · exact mul_nonneg (by positivity) (hD0 t)
        · exact mul_le_mul_of_nonneg_left (hDC t) (by positivity)
    _ = C := by rw [integral_mul_const, integral_weight T hT, one_mul]

/-- `|U_{yx}|² ≤ Q(y, x)`, since `|U_{yx}| ≤ 1` and `|U_{yx}| ≤ Q(y, x)`. -/
lemma norm_sq_evo_le_corrQ (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian) (hS : A.IsSymm)
    (t : ℝ) (y x : Fin N) : ‖evo A t y x‖ ^ 2 ≤ corrQ A hA y x := by
  have h1 := norm_evo_le_one A hS t y x
  have h2 := norm_evo_le_corrQ A hA t y x
  have h0 := norm_nonneg (evo A t y x)
  nlinarith

/-- `dimerAvg` is monotone. -/
lemma dimerAvg_mono {f g : Cfg N → ℝ} (hfg : ∀ q, IsDimer q → f q ≤ g q) :
    dimerAvg f ≤ dimerAvg g := by
  unfold dimerAvg
  gcongr with q hq
  exact hfg q (Finset.mem_filter.mp hq).2

/-- `dimerAvg` of a nonnegative function is nonnegative. -/
lemma dimerAvg_nonneg {f : Cfg N → ℝ} (hf : ∀ q, 0 ≤ f q) : 0 ≤ dimerAvg f := by
  unfold dimerAvg
  exact div_nonneg (Finset.sum_nonneg fun q _ => hf q) (Nat.cast_nonneg _)

/-- `dimerAvg` is linear. -/
lemma dimerAvg_sum_mul {ι : Type*} [Fintype ι] (a : ι → ℝ) (f : ι → Cfg N → ℝ) :
    dimerAvg (fun q => ∑ y, a y * f y q) = ∑ y, a y * dimerAvg (f y) := by
  unfold dimerAvg
  rw [Finset.sum_comm, Finset.sum_div]
  refine Finset.sum_congr rfl fun y _ => ?_
  rw [← Finset.mul_sum, mul_div_assoc]

/-- `∑_y |y - x|^p F(|y - x|) ≤ 2 ∑_d (d+1)^p F(d)`: at most two sites are at each distance. -/
lemma sum_dist_le (F : ℕ → ℝ) (hF0 : ∀ d, 0 ≤ F d) (p : ℝ) (hp : 0 ≤ p)
    (hFs : Summable (fun d : ℕ => ((d : ℝ) + 1) ^ p * F d)) (x : Fin N) :
    ∑ y : Fin N, |((y.val : ℝ) - x.val)| ^ p * F (Nat.dist y.val x.val) ≤
      2 * ∑' d : ℕ, ((d : ℝ) + 1) ^ p * F d := by
  have hg0 : ∀ d : ℕ, 0 ≤ ((d : ℝ) + 1) ^ p * F d := fun d =>
    mul_nonneg (by positivity) (hF0 d)
  have h1 : ∀ y : Fin N, |((y.val : ℝ) - x.val)| ^ p * F (Nat.dist y.val x.val) ≤
      (((Nat.dist y.val x.val : ℕ) : ℝ) + 1) ^ p * F (Nat.dist y.val x.val) := by
    intro y
    have e : |((y.val : ℝ) - x.val)| = ((Nat.dist y.val x.val : ℕ) : ℝ) := by
      rcases le_total x.val y.val with hxy | hxy
      · rw [Nat.dist_eq_sub_of_le_right hxy, Nat.cast_sub hxy,
          abs_of_nonneg (sub_nonneg.mpr (Nat.cast_le.mpr hxy))]
      · rw [Nat.dist_eq_sub_of_le hxy, Nat.cast_sub hxy,
          abs_of_nonpos (sub_nonpos.mpr (Nat.cast_le.mpr hxy)), neg_sub]
    rw [e]
    exact mul_le_mul_of_nonneg_right (Real.rpow_le_rpow (by positivity) (by linarith) hp) (hF0 _)
  have h2 : ∀ (s : Finset (Fin N)) (φ : Fin N → ℕ), Set.InjOn φ s →
      ∑ y ∈ s, (((φ y : ℕ) : ℝ) + 1) ^ p * F (φ y) ≤
        ∑' d : ℕ, ((d : ℝ) + 1) ^ p * F d := by
    intro s φ hφ
    rw [← Finset.sum_image (f := fun d : ℕ => ((d : ℝ) + 1) ^ p * F d) hφ]
    exact hFs.sum_le_tsum _ (fun d _ => hg0 d)
  have inj1 : Set.InjOn (fun y : Fin N => y.val - x.val)
      ↑(univ.filter (fun y : Fin N => x.val ≤ y.val)) := by
    intro a ha b hb hab
    simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq] at ha hb
    exact Fin.ext (by simp only at hab; omega)
  have inj2 : Set.InjOn (fun y : Fin N => x.val - y.val)
      ↑(univ.filter (fun y : Fin N => ¬ x.val ≤ y.val)) := by
    intro a ha b hb hab
    simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq] at ha hb
    exact Fin.ext (by simp only at hab; omega)
  calc ∑ y : Fin N, |((y.val : ℝ) - x.val)| ^ p * F (Nat.dist y.val x.val)
      ≤ ∑ y : Fin N,
          (((Nat.dist y.val x.val : ℕ) : ℝ) + 1) ^ p * F (Nat.dist y.val x.val) :=
        Finset.sum_le_sum fun y _ => h1 y
    _ = ∑ y ∈ univ.filter (fun y : Fin N => x.val ≤ y.val),
          (((y.val - x.val : ℕ) : ℝ) + 1) ^ p * F (y.val - x.val) +
        ∑ y ∈ univ.filter (fun y : Fin N => ¬ x.val ≤ y.val),
          (((x.val - y.val : ℕ) : ℝ) + 1) ^ p * F (x.val - y.val) := by
        rw [← Finset.sum_filter_add_sum_filter_not univ (fun y : Fin N => x.val ≤ y.val)]
        congr 1
        · refine Finset.sum_congr rfl fun y hy => ?_
          rw [Nat.dist_eq_sub_of_le_right (Finset.mem_filter.mp hy).2]
        · refine Finset.sum_congr rfl fun y hy => ?_
          rw [Nat.dist_eq_sub_of_le (not_le.mp (Finset.mem_filter.mp hy).2).le]
    _ ≤ ∑' d : ℕ, ((d : ℝ) + 1) ^ p * F d + ∑' d : ℕ, ((d : ℝ) + 1) ^ p * F d :=
        add_le_add (h2 _ _ inj1) (h2 _ _ inj2)
    _ = 2 * ∑' d : ℕ, ((d : ℝ) + 1) ^ p * F d := by ring

/-- If every dimer sector has time-averaged moment `≥ L > 0` while the dimer-averaged moment is
`≤ C` at all times, then `L ≤ C`. -/
lemma le_of_avg {n : ℕ} (hcard : 0 < ((univ.filter (fun q : Cfg n => IsDimer q)).card : ℝ))
    (T : ℝ) (hT : 0 < T) (D : Cfg n → ℝ → ℝ) (hD0 : ∀ q t, 0 ≤ D q t) (L C : ℝ) (hL : 0 < L)
    (hlow : ∀ q, IsDimer q →
      L ≤ ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) * D q t)
    (hup : ∀ t, dimerAvg (fun q => D q t) ≤ C) : L ≤ C := by
  -- each integrand is integrable: otherwise its integral would be `0`, contradicting `0 < L ≤ ∫`
  have hint : ∀ q ∈ univ.filter (fun q : Cfg n => IsDimer q),
      Integrable (fun t => (1 / T) * Real.exp (-t / T) * D q t) (volume.restrict (Set.Ioi 0)) := by
    intro q hq
    by_contra hni
    have h1 := hlow q (Finset.mem_filter.mp hq).2
    rw [integral_undef hni] at h1
    linarith
  have hswap : ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) * dimerAvg (fun q => D q t) =
      dimerAvg (fun q => ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) * D q t) := by
    unfold dimerAvg
    rw [← integral_finsetSum _ hint, ← integral_div]
    congr 1
    funext t
    rw [mul_div_assoc', Finset.mul_sum]
  have havg : L ≤
      dimerAvg (fun q => ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) * D q t) := by
    unfold dimerAvg
    rw [le_div_iff₀ hcard]
    calc L * ((univ.filter (fun q : Cfg n => IsDimer q)).card : ℝ)
        = ∑ q ∈ univ.filter (fun q : Cfg n => IsDimer q), L := by
          rw [Finset.sum_const, nsmul_eq_mul, mul_comm]
      _ ≤ _ := Finset.sum_le_sum fun q hq => hlow q (Finset.mem_filter.mp hq).2
  have hup' :
      ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) * dimerAvg (fun q => D q t) ≤ C :=
    integral_weight_le T C hT (fun t => dimerAvg (fun q => D q t))
      (fun t => dimerAvg_nonneg fun q => hD0 q t) hup
  linarith

end DimerLGTAux

/-- **Every dimer sector spreads at least diffusively, in every volume** (paper: Theorem `thm:dimer`(a)). -/
theorem dimer_sector_moment_lower (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J)
    (q : Cfg N) (hqd : IsDimer q) (x : Fin N) (p T : ℝ) (hp : 0 < p) (hT1 : 1 ≤ J * T)
    (hT2 : J * T ≤ x.val) (hT3 : J * T ≤ (N : ℝ) - 1 - x.val) :
    dimerClp (2 * |h| / J) p * (J * T) ^ (p - 1) ≤
      ∫ t in Set.Ioi (0 : ℝ), (1 / T) * Real.exp (-t / T) *
        ∑ y : Fin N, |((y.val : ℝ) - x.val)| ^ p * ‖evo (hq J h q) t y x‖ ^ 2 := by
  have hl0 : 0 < 2 * |h| / J := by positivity
  have hl1 : 2 * |h| / J < 1 := by rw [div_lt_one hJ]; exact hh1
  rw [DimerLGTAux.integral_moment_eq J h hJ q x p T]
  refine dimer_lower_bound (2 * |h| / J) hl0 hl1 p hp 0 0 ((N : ℤ) - 1) (DimerLGTAux.potQ J h q)
    (DimerLGTAux.potQ_isDimerCfg J h q hqd) (DimerLGTAux.finIcc N x) (J * T) hT1 ⟨?_, ?_⟩
  · simp only [DimerLGTAux.finIcc_apply]
    push_cast
    linarith
  · simp only [DimerLGTAux.finIcc_apply]
    push_cast
    linarith

/-- **Failure of uniform localization, sector by sector** (paper: Theorem `thm:dimer`(d), `p > 1`). -/
theorem dimer_sup_moment_unbounded (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J)
    (p : ℝ) (hp : 1 < p) (qs : ∀ N : ℕ, Cfg N) (hqs : ∀ N, IsDimer (qs N)) :
    ¬ BddAbove (Set.range (fun r : (Σ N : ℕ, Fin N × ℝ) =>
      ∑ y : Fin r.1, |((y.val : ℝ) - r.2.1.val)| ^ p * ‖evo (hq J h (qs r.1)) r.2.2 y r.2.1‖ ^ 2)) := by
  rintro ⟨C, hC⟩
  rw [mem_upperBounds] at hC
  have hl0 : 0 < 2 * |h| / J := by positivity
  have hl1 : 2 * |h| / J < 1 := by rw [div_lt_one hJ]; exact hh1
  have hc0 := DimerLGTAux.dimerClp_pos (2 * |h| / J) p hl0 hl1 (by linarith)
  obtain ⟨s, hs1, hsC⟩ := DimerLGTAux.exists_large p hp _ hc0 C
  -- the volume `n = 2M + 1` and the site `x = M`, with `M = ⌈s⌉`
  set M : ℕ := ⌈s⌉₊ with hM
  have hsM : s ≤ M := Nat.le_ceil s
  set x : Fin (2 * M + 1) := ⟨M, by omega⟩ with hx
  have hJT : J * (s / J) = s := by field_simp
  have hlow := dimer_sector_moment_lower J h hJ hh0 hh1 (qs (2 * M + 1)) (hqs _) x p (s / J)
    (by linarith) (by rw [hJT]; exact hs1) (by rw [hJT]; exact hsM)
    (by rw [hJT]; push_cast; linarith)
  have hup := DimerLGTAux.integral_weight_le (s / J) C (by positivity)
    (fun t => ∑ y : Fin (2 * M + 1), |((y.val : ℝ) - x.val)| ^ p *
      ‖evo (hq J h (qs (2 * M + 1))) t y x‖ ^ 2)
    (fun t => Finset.sum_nonneg fun y _ => by positivity)
    (fun t => hC _ ⟨⟨2 * M + 1, x, t⟩, rfl⟩)
  rw [hJT] at hlow
  linarith

/-- **Hypothesis (EC) fails for the dimer law** (paper: Theorem `thm:dimer`(d), `p > 1`): there is no profile
`F ≥ 0` with `∑_d (d+1)^p F(d) < ∞` such that `E_dim Q^{[1,N]}_q(x, y) ≤ F(|x - y|)` for all even `N`. -/
theorem not_EC_dimer (J h : ℝ) (hJ : 0 < J) (hh0 : 0 < |h|) (hh1 : 2 * |h| < J) (p : ℝ) (hp : 1 < p)
    (F : ℕ → ℝ) (hF0 : ∀ d, 0 ≤ F d) (hFs : Summable (fun d : ℕ => ((d : ℝ) + 1) ^ p * F d)) :
    ¬ ∀ (N : ℕ), Even N → ∀ x y : Fin N,
      dimerAvg (fun q : Cfg N => corrQ (hq J h q) (hq_isHermitian J h q) x y) ≤ F (Nat.dist x.val y.val) := by
  intro hEC
  have hl0 : 0 < 2 * |h| / J := by positivity
  have hl1 : 2 * |h| / J < 1 := by rw [div_lt_one hJ]; exact hh1
  have hc0 := DimerLGTAux.dimerClp_pos (2 * |h| / J) p hl0 hl1 (by linarith)
  obtain ⟨s, hs1, hsC⟩ := DimerLGTAux.exists_large p hp _ hc0
    (2 * ∑' d : ℕ, ((d : ℝ) + 1) ^ p * F d)
  -- the even volume `n = 2M + 2` and the site `x = M`, with `M = ⌈s⌉`
  set M : ℕ := ⌈s⌉₊ with hM
  have hsM : s ≤ M := Nat.le_ceil s
  have hn : Even (2 * M + 2) := ⟨M + 1, by ring⟩
  set x : Fin (2 * M + 2) := ⟨M, by omega⟩ with hx
  have hJT : J * (s / J) = s := by field_simp
  have hcard : 0 < ((univ.filter (fun q : Cfg (2 * M + 2) => IsDimer q)).card : ℝ) := by
    rw [card_isDimer hn]
    positivity
  refine absurd (DimerLGTAux.le_of_avg hcard (s / J) (by positivity)
    (fun q t => ∑ y : Fin (2 * M + 2), |((y.val : ℝ) - x.val)| ^ p *
      ‖evo (hq J h q) t y x‖ ^ 2)
    (fun q t => Finset.sum_nonneg fun y _ => by positivity)
    (dimerClp (2 * |h| / J) p * s ^ (p - 1)) (2 * ∑' d : ℕ, ((d : ℝ) + 1) ^ p * F d)
    (by positivity) (fun q hq => ?_) (fun t => ?_)) (not_le.mpr hsC)
  · have h1 := dimer_sector_moment_lower J h hJ hh0 hh1 q hq x p (s / J)
      (by linarith) (by rw [hJT]; exact hs1) (by rw [hJT]; exact hsM)
      (by rw [hJT]; push_cast; linarith)
    rwa [hJT] at h1
  · calc dimerAvg (fun q => ∑ y : Fin (2 * M + 2), |((y.val : ℝ) - x.val)| ^ p *
          ‖evo (hq J h q) t y x‖ ^ 2)
        ≤ dimerAvg (fun q => ∑ y : Fin (2 * M + 2), |((y.val : ℝ) - x.val)| ^ p *
            corrQ (hq J h q) (hq_isHermitian J h q) y x) := by
          refine DimerLGTAux.dimerAvg_mono fun q _ => Finset.sum_le_sum fun y _ => ?_
          exact mul_le_mul_of_nonneg_left
            (DimerLGTAux.norm_sq_evo_le_corrQ _ _ (hq_isSymm J h q) t y x) (by positivity)
      _ = ∑ y : Fin (2 * M + 2), |((y.val : ℝ) - x.val)| ^ p *
            dimerAvg (fun q => corrQ (hq J h q) (hq_isHermitian J h q) y x) :=
          DimerLGTAux.dimerAvg_sum_mul _ (fun y q => corrQ (hq J h q) (hq_isHermitian J h q) y x)
      _ ≤ ∑ y : Fin (2 * M + 2), |((y.val : ℝ) - x.val)| ^ p * F (Nat.dist y.val x.val) :=
          Finset.sum_le_sum fun y _ => mul_le_mul_of_nonneg_left (hEC _ hn y x) (by positivity)
      _ ≤ 2 * ∑' d : ℕ, ((d : ℝ) + 1) ^ p * F d :=
          DimerLGTAux.sum_dist_le F hF0 p (by linarith) hFs x

end DFL
