import DFL.ECProof
import DFL.Ring
import DFL.OpenImbalance

/-!
# Hypothesis (EC) on the open chain and on the ring (paper: Remark `rem:ringEC`)

* `ec_open_chain` (paper: Remark `rem:ringEC`(i)): the sector law of the open chain of Proposition
  `prop:variants`(i) is uniform on a parity hyperplane, so it has density at most `2` with respect to
  the i.i.d. law, and (EC) holds for it with the extra factor `2`.
* `ringEC_of_ldt` (paper: Remark `rem:ringEC`(ii)): the proof of Theorem `thm:EC` applies to the twisted
  ring operators `h_{q,w}` for pairs at ring distance `d < N/2`: the two boxes are disjoint arcs, both
  ends of each box are active, and an eigenvector of `h_{q,w}`, unwrapped to a function on `ℤ`
  (multiplied by `w` at every turn), solves the untwisted eigenvalue equation with the periodic
  potential, so the Poisson bound of `decay_of_not_res` applies; after a rotation of the ring the two
  arcs lie in `[0, N-1]` and the probabilistic estimate `two_box` applies.
  `ringEC_par_of_ldt`: the same for the ring sector law (uniform on a parity hyperplane, factor `2`).
-/

open Matrix Finset

noncomputable section

namespace DFL

/-- The ring distance `min(|x - y|, N - |x - y|)` on `ℤ/N`. -/
def ringDist {N : ℕ} (x y : Fin N) : ℕ := min (Nat.dist x.val y.val) (N - Nat.dist x.val y.val)

/-- The twisted ring operator `h_{q,w}` is symmetric. -/
theorem hqw_isHermitian {N : ℕ} [NeZero N] (J h : ℝ) (q : Cfg N) (ω : ZMod 2) :
    (hqw J h q ω).IsHermitian := by
  unfold Matrix.IsHermitian
  rw [Matrix.conjTranspose_eq_transpose_of_trivial]
  have hw : (wrapAdj (N := N)).IsSymm := Matrix.IsSymm.ext fun x y => by
    simp only [wrapAdj]
    congr 1
    exact propext ⟨fun h => h.elim (fun h => Or.inr ⟨h.2, h.1⟩) (fun h => Or.inl ⟨h.2, h.1⟩),
      fun h => h.elim (fun h => Or.inr ⟨h.2, h.1⟩) (fun h => Or.inl ⟨h.2, h.1⟩)⟩
  exact (hq_isSymm J h q).add (hw.smul _)

/-- **The open chain** (paper: Remark `rem:ringEC`(i)): under (EC), the parity-hyperplane average of
the eigenfunction correlator is at most `2 F(|x - y|)`. -/
theorem ec_open_chain (J h : ℝ) (F : ℕ → ℝ) (hEC : EC J h F) {N : ℕ} (hN : 0 < N) (π : ZMod 2)
    (x y : Fin N) :
    avgPar π (fun q : Cfg N => corrQ (hq J h q) (hq_isHermitian J h q) x y) ≤
      2 * F (Nat.dist x.val y.val) :=
  (OpenImbalanceAux.avgPar_le_two_avgQ hN π (fun q => corrQ_nonneg _ _ x y)).trans
    (mul_le_mul_of_nonneg_left (hEC.bound N x y) (by norm_num))

namespace RingECAux

/-! ### Sites of the ring as residues of integers -/

/-- The site `z mod N` of the ring. -/
def modFin (N : ℕ) [NeZero N] (z : ℤ) : Fin N :=
  ⟨(z % (N : ℤ)).toNat, by
    have hN : (0 : ℤ) < N := by exact_mod_cast NeZero.pos N
    have h1 := Int.emod_nonneg z hN.ne'
    have h2 := Int.emod_lt_of_pos z hN
    omega⟩

variable {N : ℕ} [NeZero N]

/-- The value of `modFin`. -/
lemma modFin_val (z : ℤ) : ((modFin N z : ℕ) : ℤ) = z % (N : ℤ) := by
  have hN : (0 : ℤ) < N := by exact_mod_cast NeZero.pos N
  have h1 := Int.emod_nonneg z hN.ne'
  simp only [modFin]
  omega

/-- `modFin` only depends on the residue. -/
lemma modFin_congr {z z' : ℤ} (h : z % (N : ℤ) = z' % (N : ℤ)) : modFin N z = modFin N z' := by
  apply Fin.ext
  have h1 := modFin_val (N := N) z
  have h2 := modFin_val (N := N) z'
  omega

/-- `modFin` at a site of the ring. -/
lemma modFin_natCast (i : Fin N) : modFin N (i.val : ℤ) = i := by
  apply Fin.ext
  have h := modFin_val (N := N) (i.val : ℤ)
  rw [Int.emod_eq_of_lt (by positivity) (by exact_mod_cast i.isLt)] at h
  exact_mod_cast h

/-- `modFin` is `N`-periodic. -/
lemma modFin_add_mul (z k : ℤ) : modFin N (z + N * k) = modFin N z :=
  modFin_congr (Int.add_mul_emod_self_left z N k)

/-- `modFin` of a shifted residue. -/
lemma modFin_val_add (z t : ℤ) : modFin N (((modFin N z : ℕ) : ℤ) + t) = modFin N (z + t) := by
  apply modFin_congr
  rw [modFin_val, Int.emod_add_emod]

/-- `modFin` in the window `[0, N-1]`. -/
lemma modFin_of_mem (z : ℤ) (h0 : 0 ≤ z) (h1 : z < N) :
    modFin N z = ⟨z.toNat, by omega⟩ := by
  apply Fin.ext
  have h := modFin_val (N := N) z
  rw [Int.emod_eq_of_lt h0 h1] at h
  simp only
  omega

/-! ### Unwrapping a function on the ring -/

/-- The unwrapping `φ̃(z) = w^{⌊z/N⌋} φ(z mod N)` of a function `φ` on the ring with twist `w`. -/
def unwrap (w : ℝ) (φ : Fin N → ℝ) (z : ℤ) : ℝ := w ^ (z / (N : ℤ)) * φ (modFin N z)

/-- Quasi-periodicity of the unwrapping: `φ̃(z + kN) = w^k φ̃(z)`. -/
lemma unwrap_add_mul (w : ℝ) (hw : w ≠ 0) (φ : Fin N → ℝ) (z k : ℤ) :
    unwrap w φ (z + N * k) = w ^ k * unwrap w φ z := by
  have hN : (N : ℤ) ≠ 0 := by exact_mod_cast NeZero.ne N
  simp only [unwrap, modFin_add_mul, Int.add_mul_ediv_left z k hN, zpow_add₀ hw]
  ring

/-- The unwrapping in the window `[0, N-1]`. -/
lemma unwrap_of_mem (w : ℝ) (φ : Fin N → ℝ) (z : ℤ) (h0 : 0 ≤ z) (h1 : z < N) :
    unwrap w φ z = φ (modFin N z) := by
  simp [unwrap, Int.ediv_eq_zero_of_lt h0 h1]

/-- The unwrapping in the window `[0, N-1]` is the extension by zero. -/
lemma unwrap_eq_extZ (w : ℝ) (φ : Fin N → ℝ) (z : ℤ) (h0 : 0 ≤ z) (h1 : z < N) :
    unwrap w φ z = extZ φ z := by
  rw [unwrap_of_mem w φ z h0 h1, modFin_of_mem z h0 h1, extZ, dite_eq_left ⟨h0, h1⟩]

/-- The unwrapping at a site of the ring. -/
lemma unwrap_natCast (w : ℝ) (φ : Fin N → ℝ) (i : Fin N) : unwrap w φ (i.val : ℤ) = φ i := by
  rw [unwrap_of_mem w φ _ (by positivity) (by exact_mod_cast i.isLt), modFin_natCast]

/-- For `|w| = 1`, the unwrapping has the modulus of `φ`. -/
lemma abs_unwrap (w : ℝ) (hw : |w| = 1) (φ : Fin N → ℝ) (z : ℤ) :
    |unwrap w φ z| = |φ (modFin N z)| := by
  rw [unwrap, abs_mul, abs_zpow, hw, _root_.one_zpow, one_mul]

/-- `sgn a ≠ 0`. -/
lemma sgn_ne_zero (a : ZMod 2) : sgn a ≠ 0 := by
  unfold sgn; split_ifs <;> norm_num

/-- `(sgn a)⁻¹ = sgn a`. -/
lemma sgn_inv (a : ZMod 2) : (sgn a)⁻¹ = sgn a := by
  unfold sgn; split_ifs <;> norm_num

/-- The unwrapping at `x + 1`: across the ring bond `(N - 1, 0)` it picks up the twist. -/
lemma unwrap_succ (w : ℝ) (hw : w ≠ 0) (φ : Fin N → ℝ) (x : Fin N) :
    unwrap w φ ((x.val : ℤ) + 1) =
      extZ φ ((x.val : ℤ) + 1) + (if x = lastSite then w * φ 0 else 0) := by
  by_cases hx : x.val + 1 < N
  · have hne : x ≠ lastSite := by
      intro h; rw [h, RingAux.lastSite_val] at hx; omega
    rw [ite_eq_right hne, add_zero, unwrap_eq_extZ w φ _ (by positivity) (by exact_mod_cast hx)]
  · have hx' : x = lastSite := by
      apply Fin.ext; rw [RingAux.lastSite_val]; have := x.isLt; omega
    rw [ite_eq_left hx', extZ, dite_eq_right (by omega), zero_add]
    have e : (x.val : ℤ) + 1 = 0 + N * 1 := by have := x.isLt; omega
    rw [e, unwrap_add_mul w hw, zpow_one, unwrap_of_mem w φ 0 le_rfl
      (by exact_mod_cast NeZero.pos N)]
    congr 2

/-- The unwrapping at `x - 1`: across the ring bond `(N - 1, 0)` it picks up the twist. -/
lemma unwrap_pred (w : ℝ) (hw : w ≠ 0) (φ : Fin N → ℝ) (x : Fin N) :
    unwrap w φ ((x.val : ℤ) - 1) =
      extZ φ ((x.val : ℤ) - 1) + (if x = 0 then w⁻¹ * φ lastSite else 0) := by
  by_cases hx : 1 ≤ x.val
  · have hne : x ≠ 0 := by
      intro h; rw [h, Fin.val_zero] at hx; omega
    rw [ite_eq_right hne, add_zero, unwrap_eq_extZ w φ _ (by omega)
      (by have := x.isLt; omega)]
  · have hx' : x = 0 := by
      apply Fin.ext; rw [Fin.val_zero]; omega
    rw [ite_eq_left hx', extZ, dite_eq_right (by omega), zero_add]
    have hN := NeZero.pos N
    have e : (x.val : ℤ) - 1 = ((N - 1 : ℕ) : ℤ) + N * (-1) := by omega
    rw [e, unwrap_add_mul w hw, _root_.zpow_neg_one, unwrap_of_mem w φ _ (by positivity) (by omega)]
    congr 2
    apply Fin.ext
    rw [modFin_of_mem _ (by positivity) (by omega), RingAux.lastSite_val]
    simp

/-! ### The eigenvalue equation of `h_{q,w}` -/

/-- The ring bond acting on a vector (`N ≥ 2`). -/
lemma wrapAdj_mulVec_apply (hN : 2 ≤ N) (φ : Fin N → ℝ) (x : Fin N) :
    (wrapAdj *ᵥ φ) x = (if x = lastSite then φ 0 else 0) + (if x = 0 then φ lastSite else 0) := by
  have hl : (lastSite : Fin N) ≠ 0 := RingAux.lastSite_ne_zero hN
  have hsplit : ∀ y, wrapAdj x y * φ y =
      (if y = 0 then (if x = lastSite then φ y else 0) else 0) +
        (if y = lastSite then (if x = 0 then φ y else 0) else 0) := by
    intro y
    simp only [wrapAdj]
    by_cases h1 : x = lastSite <;> by_cases h2 : y = 0 <;> by_cases h3 : x = 0 <;>
      by_cases h4 : y = lastSite <;> simp_all
  simp only [Matrix.mulVec, dotProduct, hsplit, Finset.sum_add_distrib, Finset.sum_ite_eq',
    Finset.mem_univ, ite_true]

/-- `h_{q,w}` acting on a vector, through the unwrapping (`N ≥ 2`). -/
lemma hqw_mulVec_apply (hN : 2 ≤ N) (J h : ℝ) (q : Cfg N) (ω : ZMod 2) (φ : Fin N → ℝ)
    (x : Fin N) :
    (hqw J h q ω *ᵥ φ) x =
      -J * (unwrap (sgn ω) φ ((x.val : ℤ) + 1) + unwrap (sgn ω) φ ((x.val : ℤ) - 1)) +
        2 * h * sgn (q x) * φ x := by
  rw [hqw, Matrix.add_mulVec, Matrix.smul_mulVec]
  simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul]
  rw [TransferAux.hq_mulVec_apply, wrapAdj_mulVec_apply hN, unwrap_succ _ (sgn_ne_zero ω),
    unwrap_pred _ (sgn_ne_zero ω), sgn_inv]
  split_ifs <;> ring

/-- An eigenvector of `h_{q,w}` satisfies the untwisted equation at the sites of the ring. -/
lemma unwrap_eq_site (hN : 2 ≤ N) (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) (ω : ZMod 2)
    (φ : Fin N → ℝ) (μ : ℝ) (hφ : hqw J h q ω *ᵥ φ = μ • φ) (x : Fin N) :
    unwrap (sgn ω) φ ((x.val : ℤ) + 1) + unwrap (sgn ω) φ ((x.val : ℤ) - 1) +
        -(2 * h / J) * sgn (q x) * unwrap (sgn ω) φ (x.val : ℤ) =
      (-μ / J) * unwrap (sgn ω) φ (x.val : ℤ) := by
  have hev := congrFun hφ x
  rw [hqw_mulVec_apply hN] at hev
  simp only [Pi.smul_apply, smul_eq_mul] at hev
  rw [unwrap_natCast]
  have hJJ : J * J⁻¹ = 1 := mul_inv_cancel₀ hJ
  linear_combination (-J⁻¹) * hev -
    (unwrap (sgn ω) φ ((x.val : ℤ) + 1) + unwrap (sgn ω) φ ((x.val : ℤ) - 1)) * hJJ

/-- **Unwrapping** (paper: Remark `rem:ringEC`(ii)): an eigenvector of `h_{q,w}`, unwrapped to `ℤ`,
solves `φ̃(z+1) + φ̃(z-1) + λ sgn(q_{z mod N}) φ̃(z) = E φ̃(z)` everywhere, `λ = -(2h/J)`,
`E = -E'/J`. -/
lemma unwrap_eq (hN : 2 ≤ N) (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) (ω : ZMod 2) (φ : Fin N → ℝ)
    (μ : ℝ) (hφ : hqw J h q ω *ᵥ φ = μ • φ) (z : ℤ) :
    unwrap (sgn ω) φ (z + 1) + unwrap (sgn ω) φ (z - 1) +
        -(2 * h / J) * sgn (q (modFin N z)) * unwrap (sgn ω) φ z =
      (-μ / J) * unwrap (sgn ω) φ z := by
  obtain ⟨x, k, rfl⟩ : ∃ (x : Fin N) (k : ℤ), z = (x.val : ℤ) + N * k :=
    ⟨modFin N z, z / N, by rw [modFin_val]; exact (Int.emod_add_mul_ediv z N).symm⟩
  have hw := sgn_ne_zero ω
  have e := unwrap_eq_site hN J h hJ q ω φ μ hφ x
  rw [modFin_add_mul, modFin_natCast,
    show (x.val : ℤ) + N * k + 1 = ((x.val : ℤ) + 1) + N * k by ring,
    show (x.val : ℤ) + N * k - 1 = ((x.val : ℤ) - 1) + N * k by ring,
    unwrap_add_mul _ hw, unwrap_add_mul _ hw, unwrap_add_mul _ hw]
  linear_combination (sgn ω ^ k) * e

/-! ### The deterministic step on the ring -/

/-- The rescaled eigenvalue `E = -E'/J` of `h_{q,w}` lies in `I_λ`, `λ = -(2h/J)` (`N ≥ 2`). -/
lemma eigval_mem_Ilam_ring (hN : 2 ≤ N) (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) (ω : ZMod 2)
    (α : Fin N) :
    -(hqw_isHermitian J h q ω).eigenvalues α / J ∈ Ilam (-(2 * h / J)) := by
  set hA := hqw_isHermitian J h q ω with hA_def
  set φ : Fin N → ℝ := fun i => hA.eigenvectorBasis α i with hφ_def
  obtain ⟨x, -, hxmax⟩ := Finset.exists_max_image Finset.univ (fun y => |φ y|)
    ⟨α, Finset.mem_univ _⟩
  have hnorm : ‖hA.eigenvectorBasis α‖ = 1 := hA.eigenvectorBasis.orthonormal.1 α
  have hpos : 0 < |φ x| := by
    refine lt_of_le_of_ne (abs_nonneg _) (fun hcon => ?_)
    have h0 : ∀ y, hA.eigenvectorBasis α y = 0 := fun y =>
      abs_nonpos_iff.mp ((hxmax y (Finset.mem_univ y)).trans hcon.symm.le)
    rw [EuclideanSpace.norm_eq] at hnorm
    simp [h0] at hnorm
  have hev : hqw J h q ω *ᵥ φ = hA.eigenvalues α • φ := hA.mulVec_eigenvectorBasis α
  have e := unwrap_eq_site hN J h hJ q ω φ _ hev x
  rw [unwrap_natCast] at e
  set E := -hA.eigenvalues α / J with hE
  set lam := -(2 * h / J) with hlam
  have hb : ∀ t : ℤ, |unwrap (sgn ω) φ t| ≤ |φ x| := fun t => by
    rw [abs_unwrap _ (TransferAux.abs_sgn ω)]
    exact hxmax _ (Finset.mem_univ _)
  have h1 : |E - lam * sgn (q x)| * |φ x| ≤ 2 * |φ x| := by
    rw [← abs_mul, show (E - lam * sgn (q x)) * φ x =
      unwrap (sgn ω) φ ((x.val : ℤ) + 1) + unwrap (sgn ω) φ ((x.val : ℤ) - 1) by
        linear_combination -e]
    calc |unwrap (sgn ω) φ ((x.val : ℤ) + 1) + unwrap (sgn ω) φ ((x.val : ℤ) - 1)|
        ≤ |unwrap (sgn ω) φ ((x.val : ℤ) + 1)| + |unwrap (sgn ω) φ ((x.val : ℤ) - 1)| :=
          abs_add_le _ _
      _ ≤ |φ x| + |φ x| := add_le_add (hb _) (hb _)
      _ = 2 * |φ x| := by ring
  have h2 : |E - lam * sgn (q x)| ≤ 2 := le_of_mul_le_mul_right h1 hpos
  have h3 : |lam * sgn (q x)| = |lam| := by rw [abs_mul, TransferAux.abs_sgn, mul_one]
  have h4 : |E| ≤ 2 + |lam| := by
    calc |E| = |(E - lam * sgn (q x)) + lam * sgn (q x)| := by ring_nf
      _ ≤ |E - lam * sgn (q x)| + |lam * sgn (q x)| := abs_add_le _ _
      _ ≤ 2 + |lam| := by rw [h3]; linarith
  exact ⟨by linarith [(abs_le.mp h4).1], (abs_le.mp h4).2⟩

/-- The configuration seen from the rotated ring: `(rotQ s q)_i = q_{(i + s) mod N}`. -/
def rotQ (s : ℤ) (q : Cfg N) : Cfg N := fun i => q (modFin N ((i.val : ℤ) + s))

/-- Rotating back. -/
lemma rotQ_neg_rotQ (s : ℤ) (q : Cfg N) : rotQ (-s) (rotQ s q) = q := by
  funext i
  simp only [rotQ]
  rw [modFin_val_add, show (i.val : ℤ) + -s + s = i.val by ring, modFin_natCast]

/-- The potential of the rotated configuration on `[0, N-1]` is the periodic potential. -/
lemma potZ_rotQ (lam : ℝ) (s : ℤ) (q : Cfg N) (z : ℤ) (h0 : 0 ≤ z) (h1 : z < N) :
    potZ lam (rotQ s q) z = lam * sgn (q (modFin N (s + z))) := by
  unfold potZ extZ
  rw [dite_eq_left ⟨h0, h1⟩]
  simp only [rotQ]
  rw [Int.toNat_of_nonneg h0, add_comm z s]

/-- **Decay across an arc** (paper: Remark `rem:ringEC`(ii) with Lemma `lem:badset`(a)): if the
rescaled eigenvalue is not resonant for the box `[a, a+2n]` of the configuration rotated by `s`,
then `|φ_α(x)| ≤ θ (|φ_α(u)| + |φ_α(v)|)` at the centre `x = s + a + n` of the arc, where
`u = s + a - 1` and `v = s + a + 2n + 1` are its outer neighbours (mod `N`). -/
lemma ring_box_decay (hN : 2 ≤ N) (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) (ω : ZMod 2) (α : Fin N)
    (s a : ℤ) (n : ℕ) (θ : ℝ) (ha : 0 ≤ a) (haN : a + ((2 * n + 1 : ℕ) : ℤ) ≤ N) (x : Fin N)
    (hx : modFin N (s + (a + n)) = x)
    (hres : -(hqw_isHermitian J h q ω).eigenvalues α / J ∉
      resSet (potZ (-(2 * h / J)) (rotQ s q)) (-(2 * h / J)) a (2 * n + 1) n true true θ) :
    |(hqw_isHermitian J h q ω).eigenvectorBasis α x| ≤
      θ * (|(hqw_isHermitian J h q ω).eigenvectorBasis α (modFin N (s + (a - 1)))| +
        |(hqw_isHermitian J h q ω).eigenvectorBasis α
          (modFin N (s + (a + ((2 * n + 1 : ℕ) : ℤ))))|) := by
  set hA := hqw_isHermitian J h q ω with hA_def
  set φ : Fin N → ℝ := fun i => hA.eigenvectorBasis α i with hφ_def
  have hev : hqw J h q ω *ᵥ φ = hA.eigenvalues α • φ := hA.mulVec_eigenvectorBasis α
  have hrec : ∀ z : ℤ, a ≤ z → z ≤ a + ((2 * n + 1 : ℕ) : ℤ) - 1 →
      unwrap (sgn ω) φ (s + (z + 1)) + unwrap (sgn ω) φ (s + (z - 1)) +
          potZ (-(2 * h / J)) (rotQ s q) z * unwrap (sgn ω) φ (s + z) =
        (-hA.eigenvalues α / J) * unwrap (sgn ω) φ (s + z) := by
    intro z hz1 hz2
    have e := unwrap_eq hN J h hJ q ω φ _ hev (s + z)
    rw [potZ_rotQ _ s q z (by omega) (by omega),
      show s + (z + 1) = s + z + 1 by ring, show s + (z - 1) = s + z - 1 by ring]
    exact e
  have key := decay_of_not_res (potZ (-(2 * h / J)) (rotQ s q)) (-(2 * h / J)) a (2 * n + 1) n
    true true θ (by omega) (Or.inl rfl) _ (eigval_mem_Ilam_ring hN J h hJ q ω α) hres
    (fun z => unwrap (sgn ω) φ (s + z)) hrec (fun hf => by cases hf) (fun hf => by cases hf)
  simp only [abs_unwrap _ (TransferAux.abs_sgn ω), hx] at key
  exact key

omit [NeZero N] in
/-- **Cauchy–Schwarz step** (paper: Eq. `eq:Qgood`): if every eigenvector decays at `x` (towards
`u₁, u₂`) or at `y` (towards `u₃, u₄`) with the factor `θ`, then `Q(x, y) ≤ 4θ`. -/
lemma corrQ_le_of_alt (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsHermitian)
    (x y u1 u2 u3 u4 : Fin N) (θ : ℝ) (hθ : 0 ≤ θ)
    (halt : ∀ α, |hA.eigenvectorBasis α x| ≤
        θ * (|hA.eigenvectorBasis α u1| + |hA.eigenvectorBasis α u2|) ∨
      |hA.eigenvectorBasis α y| ≤
        θ * (|hA.eigenvectorBasis α u3| + |hA.eigenvectorBasis α u4|)) :
    corrQ A hA x y ≤ 4 * θ := by
  set φ := hA.eigenvectorBasis with hφ_def
  have hpt : ∀ α, |φ α x| * |φ α y| ≤
      θ * (|φ α u1| * |φ α y| + |φ α u2| * |φ α y|) +
        θ * (|φ α x| * |φ α u3| + |φ α x| * |φ α u4|) := by
    intro α
    have hx0 := abs_nonneg (φ α x)
    have hy0 := abs_nonneg (φ α y)
    have h1 := abs_nonneg (φ α u1)
    have h2 := abs_nonneg (φ α u2)
    have h3 := abs_nonneg (φ α u3)
    have h4 := abs_nonneg (φ α u4)
    rcases halt α with hd | hd
    · have hd' := mul_le_mul_of_nonneg_right hd hy0
      have : 0 ≤ θ * (|φ α x| * |φ α u3| + |φ α x| * |φ α u4|) := by positivity
      have e : θ * (|φ α u1| + |φ α u2|) * |φ α y| =
          θ * (|φ α u1| * |φ α y| + |φ α u2| * |φ α y|) := by ring
      linarith
    · have hd' := mul_le_mul_of_nonneg_left hd hx0
      have : 0 ≤ θ * (|φ α u1| * |φ α y| + |φ α u2| * |φ α y|) := by positivity
      have e : |φ α x| * (θ * (|φ α u3| + |φ α u4|)) =
          θ * (|φ α x| * |φ α u3| + |φ α x| * |φ α u4|) := by ring
      linarith
  have s1 := ECProofAux.sum_abs_eigvec_mul_le_one A hA u1 y
  have s2 := ECProofAux.sum_abs_eigvec_mul_le_one A hA u2 y
  have s3 := ECProofAux.sum_abs_eigvec_mul_le_one A hA x u3
  have s4 := ECProofAux.sum_abs_eigvec_mul_le_one A hA x u4
  calc corrQ A hA x y
      = ∑ α, |φ α x| * |φ α y| := rfl
    _ ≤ ∑ α, (θ * (|φ α u1| * |φ α y| + |φ α u2| * |φ α y|) +
          θ * (|φ α x| * |φ α u3| + |φ α x| * |φ α u4|)) :=
        Finset.sum_le_sum fun α _ => hpt α
    _ = θ * ((∑ α, |φ α u1| * |φ α y|) + ∑ α, |φ α u2| * |φ α y|) +
          θ * ((∑ α, |φ α x| * |φ α u3|) + ∑ α, |φ α x| * |φ α u4|) := by
        rw [Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum, Finset.sum_add_distrib,
          Finset.sum_add_distrib]
    _ ≤ θ * (1 + 1) + θ * (1 + 1) := by gcongr
    _ = 4 * θ := by ring

/-- **Deterministic step on the ring** (paper: Remark `rem:ringEC`(ii), Eq. `eq:Qgood`): if no
energy is resonant for both arcs `[0, 2n]` and `[d, d+2n]` of the configuration rotated by `s`,
where `x = s + n` and `y = s + d + n` mod `N`, then `Q(x, y) ≤ 4θ`. -/
lemma corrQ_ring_le_of_good (hN : 2 ≤ N) (J h : ℝ) (hJ : J ≠ 0) (q : Cfg N) (ω : ZMod 2)
    (x y : Fin N) (s : ℤ) (n d : ℕ) (θ : ℝ) (hθ : 0 ≤ θ) (hx : modFin N (s + n) = x)
    (hy : modFin N (s + (d + n)) = y) (hdN : d + (2 * n + 1) ≤ N)
    (hgood : ¬ (resSet (potZ (-(2 * h / J)) (rotQ s q)) (-(2 * h / J)) 0 (2 * n + 1) n true true
        θ ∩ resSet (potZ (-(2 * h / J)) (rotQ s q)) (-(2 * h / J)) d (2 * n + 1) n true true
        θ).Nonempty) :
    corrQ (hqw J h q ω) (hqw_isHermitian J h q ω) x y ≤ 4 * θ := by
  refine corrQ_le_of_alt _ _ x y (modFin N (s + (0 - 1)))
    (modFin N (s + (0 + ((2 * n + 1 : ℕ) : ℤ)))) (modFin N (s + ((d : ℤ) - 1)))
    (modFin N (s + ((d : ℤ) + ((2 * n + 1 : ℕ) : ℤ)))) θ hθ fun α => ?_
  by_cases h1 : -(hqw_isHermitian J h q ω).eigenvalues α / J ∈
      resSet (potZ (-(2 * h / J)) (rotQ s q)) (-(2 * h / J)) 0 (2 * n + 1) n true true θ
  · right
    have h2 : -(hqw_isHermitian J h q ω).eigenvalues α / J ∉
        resSet (potZ (-(2 * h / J)) (rotQ s q)) (-(2 * h / J)) d (2 * n + 1) n true true θ :=
      fun h2 => hgood ⟨_, h1, h2⟩
    exact ring_box_decay hN J h hJ q ω α s d n θ (by positivity) (by omega) y hy h2
  · left
    exact ring_box_decay hN J h hJ q ω α s 0 n θ le_rfl (by omega) x (by rw [zero_add]; exact hx) h1

/-! ### The probabilistic step on the ring -/

/-- **Probabilistic step on the ring** (paper: Remark `rem:ringEC`(ii), Eq. `eq:ECn`): for
`x = s + n`, `y = s + d + n` mod `N` with `2n + 1 ≤ d` and `d + 2n + 1 ≤ N`,
`E_q Q_q(x, y) ≤ 4 e^{-κn} + 4(8n+6) C₁ e^{-c₁n}`. -/
lemma avgQ_ring_le_of_entry (hN : 2 ≤ N) (J h : ℝ) (hJ : J ≠ 0) (κ C1 c1 : ℝ) (m0 : ℕ)
    (hEB : EntryBound (-(2 * h / J)) κ C1 c1 m0) (hC1 : 0 ≤ C1) (ω : ZMod 2) (x y : Fin N)
    (s : ℤ) (n d : ℕ) (hx : modFin N (s + n) = x) (hy : modFin N (s + (d + n)) = y)
    (h2n : 2 * n + 1 ≤ d) (hdN : d + (2 * n + 1) ≤ N) (hn : m0 ≤ n) :
    avgQ (fun q => corrQ (hqw J h q ω) (hqw_isHermitian J h q ω) x y) ≤
      4 * Real.exp (-(κ * n)) + 4 * (8 * n + 6) * C1 * Real.exp (-c1 * n) := by
  have hB := two_box (N := N) (-(2 * h / J)) κ C1 c1 m0 hEB hC1 0 (2 * n + 1) n true true (d : ℤ)
    (2 * n + 1) n true true n le_rfl (by omega) (by omega) (by omega) (by omega) le_rfl le_rfl hn
    (fun _ => by omega) (fun _ => rfl) (fun _ => by omega) (fun _ => rfl)
  refine ECProofAux.avgQ_le_add _ (Finset.image (rotQ (-s)) _) _ _
    (fun q => ECProofAux.corrQ_le_one _ _ x y) (fun q hq => ?_) (by positivity)
    ((Nat.cast_le.mpr Finset.card_image_le).trans hB)
  refine corrQ_ring_le_of_good hN J h hJ q ω x y s n d _ (Real.exp_pos _).le hx hy hdN ?_
  intro hne
  apply hq
  rw [Finset.mem_image]
  exact ⟨rotQ s q, by simpa using hne, rotQ_neg_rotQ s q⟩

/-! ### Geometry of the ring -/

/-- Integer lifts at ring distance `d`: `y = x + d` or `x = y + d` modulo `N`. -/
lemma ringDist_cases (x y : Fin N) :
    modFin N ((x.val : ℤ) + ringDist x y) = y ∨ modFin N ((y.val : ℤ) + ringDist x y) = x := by
  have hx := x.isLt
  have hy := y.isLt
  have key : ((x.val : ℤ) + ringDist x y = y.val ∨ (x.val : ℤ) + ringDist x y = y.val + N * 1) ∨
      ((y.val : ℤ) + ringDist x y = x.val ∨ (y.val : ℤ) + ringDist x y = x.val + N * 1) := by
    unfold ringDist Nat.dist
    omega
  rcases key with (h | h) | (h | h)
  · left; rw [h, modFin_natCast]
  · left; rw [h, modFin_add_mul, modFin_natCast]
  · right; rw [h, modFin_natCast]
  · right; rw [h, modFin_add_mul, modFin_natCast]

end RingECAux

/-- **(EC) on the ring** (paper: Remark `rem:ringEC`(ii)): given the large-deviation estimate, the
sector-averaged eigenfunction correlator of the twisted ring operators `h_{q,w}` decays exponentially
in the ring distance, for pairs at ring distance `d < N/2`, uniformly in `N` and in the twist. -/
theorem ringEC_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ ∀ (N : ℕ) [NeZero N] (ω : ZMod 2) (x y : Fin N), 2 * ringDist x y < N →
      avgQ (fun q : Cfg N => corrQ (hqw J h q ω) (hqw_isHermitian J h q ω) x y) ≤
        C * Real.exp (-η * ringDist x y) := by
  obtain ⟨κ, C1, c1, m0, hκ, hc1, hC1, hEB⟩ := entry_ldt _ hL
  set c := min κ (c1 / 2) with hc_def
  have hc : 0 < c := lt_min hκ (by linarith)
  have hcκ : c ≤ κ := min_le_left _ _
  have hcc1 : c ≤ c1 / 2 := min_le_right _ _
  set C' := 4 + 4 * C1 * (6 + 16 / c1) with hC'_def
  have hC' : 0 ≤ C' := by positivity
  set η := c / 2 with hη_def
  have hη : 0 < η := by positivity
  refine ⟨max (C' * Real.exp c) (Real.exp (η * (2 * m0 + 1))), η, hη, ?_⟩
  intro N _ ω x y hdN
  set d := ringDist x y with hd_def
  rcases le_or_gt d (2 * m0) with hd | hd
  · -- short distances: `Q ≤ 1`
    refine (ECProofAux.avgQ_le_of_le _ 1 fun q => ECProofAux.corrQ_le_one _ _ x y).trans ?_
    refine le_trans ?_ (mul_le_mul_of_nonneg_right (le_max_right _ _) (Real.exp_pos _).le)
    rw [← Real.exp_add]
    apply Real.one_le_exp
    have : (d : ℝ) ≤ 2 * m0 := by exact_mod_cast hd
    nlinarith
  · -- `n = ⌊(d-1)/2⌋ ≥ m₀`: two disjoint arcs of length `2n + 1`
    obtain ⟨n, hn_def⟩ : ∃ n : ℕ, n = (d - 1) / 2 := ⟨_, rfl⟩
    have hn : m0 ≤ n := by omega
    have h2n : 2 * n + 1 ≤ d := by omega
    have hN2 : 2 ≤ N := by omega
    have hdN' : d + (2 * n + 1) ≤ N := by omega
    have hbound : avgQ (fun q : Cfg N => corrQ (hqw J h q ω) (hqw_isHermitian J h q ω) x y) ≤
        4 * Real.exp (-(κ * n)) + 4 * (8 * n + 6) * C1 * Real.exp (-c1 * n) := by
      rcases RingECAux.ringDist_cases x y with hxy | hxy
      · refine RingECAux.avgQ_ring_le_of_entry hN2 J h hJ κ C1 c1 m0 hEB hC1 ω x y
          ((x.val : ℤ) - n) n d ?_ ?_ h2n hdN' hn
        · rw [sub_add_cancel, RingECAux.modFin_natCast]
        · rw [show (x.val : ℤ) - n + (d + n) = x.val + d by ring]
          exact hxy
      · have hsym : avgQ (fun q : Cfg N => corrQ (hqw J h q ω) (hqw_isHermitian J h q ω) x y) =
            avgQ (fun q : Cfg N => corrQ (hqw J h q ω) (hqw_isHermitian J h q ω) y x) := by
          unfold avgQ
          congr 1
          exact Finset.sum_congr rfl fun q _ => ECProofAux.corrQ_symm _ _ _ _
        rw [hsym]
        refine RingECAux.avgQ_ring_le_of_entry hN2 J h hJ κ C1 c1 m0 hEB hC1 ω y x
          ((y.val : ℤ) - n) n d ?_ ?_ h2n hdN' hn
        · rw [sub_add_cancel, RingECAux.modFin_natCast]
        · rw [show (y.val : ℤ) - n + (d + n) = y.val + d by ring]
          exact hxy
    refine hbound.trans ?_
    have hn0 : (0 : ℝ) ≤ n := Nat.cast_nonneg n
    have hnR : (d : ℝ) ≤ 2 * (n : ℝ) + 2 := by
      have : d ≤ 2 * n + 2 := by omega
      exact_mod_cast this
    have e1 : Real.exp (-(κ * n)) ≤ Real.exp (-(c * n)) := by
      apply Real.exp_le_exp.mpr
      nlinarith
    have e2 := ECProofAux.poly_exp_le c1 c hc1 hcc1 n
    have e3 : Real.exp (-(c * n)) ≤ Real.exp c * Real.exp (-η * d) := by
      rw [← Real.exp_add]
      apply Real.exp_le_exp.mpr
      rw [hη_def]
      nlinarith
    calc 4 * Real.exp (-(κ * n)) + 4 * (8 * n + 6) * C1 * Real.exp (-c1 * n)
        ≤ 4 * Real.exp (-(c * n)) + 4 * C1 * ((6 + 16 / c1) * Real.exp (-(c * n))) := by
          have := mul_le_mul_of_nonneg_left e2 (by positivity : (0 : ℝ) ≤ 4 * C1)
          nlinarith
      _ = C' * Real.exp (-(c * n)) := by rw [hC'_def]; ring
      _ ≤ C' * (Real.exp c * Real.exp (-η * d)) := by gcongr
      _ = C' * Real.exp c * Real.exp (-η * d) := by ring
      _ ≤ max (C' * Real.exp c) (Real.exp (η * (2 * m0 + 1))) * Real.exp (-η * d) :=
          mul_le_mul_of_nonneg_right (le_max_left _ _) (Real.exp_pos _).le

/-- **(EC) for the ring sector law** (paper: Remark `rem:ringEC`(ii) and Proposition
`prop:variants`(ii)): the uniform law on a parity hyperplane costs a factor `2`. -/
theorem ringEC_par_of_ldt (J h : ℝ) (hJ : J ≠ 0) (hL : LDT (-(2 * h / J))) :
    ∃ C η : ℝ, 0 < η ∧ ∀ (N : ℕ) [NeZero N] (ω π : ZMod 2) (x y : Fin N), 2 * ringDist x y < N →
      avgPar π (fun q : Cfg N => corrQ (hqw J h q ω) (hqw_isHermitian J h q ω) x y) ≤
        C * Real.exp (-η * ringDist x y) := by
  obtain ⟨C, η, hη, hb⟩ := ringEC_of_ldt J h hJ hL
  refine ⟨2 * C, η, hη, fun N _ ω π x y hd => ?_⟩
  calc avgPar π (fun q : Cfg N => corrQ (hqw J h q ω) (hqw_isHermitian J h q ω) x y)
      ≤ 2 * avgQ (fun q : Cfg N => corrQ (hqw J h q ω) (hqw_isHermitian J h q ω) x y) :=
        OpenImbalanceAux.avgPar_le_two_avgQ (NeZero.pos N) π (fun q => corrQ_nonneg _ _ x y)
    _ ≤ 2 * (C * Real.exp (-η * ringDist x y)) :=
        mul_le_mul_of_nonneg_left (hb N ω x y hd) (by norm_num)
    _ = 2 * C * Real.exp (-η * ringDist x y) := by ring

end DFL
