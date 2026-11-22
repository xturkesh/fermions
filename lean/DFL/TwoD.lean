import DFL.Graph

/-!
# Two dimensions: no compactly supported eigenvectors; the checkerboard chiral gap
(paper: Proposition `prop:nocompact`; Proposition `prop:2dlaws`(c)(i),(ii), and the gauge covariance
`D_θ h_{q,s} D_θ = h_{q,θ·s}` below Eq. `eq:hqs`)

* `no_compact_eigenvector`: a nearest-neighbour operator on `ℤ²` with nonvanishing hoppings has no
  finitely supported eigenvector (look at the lexicographically extremal point of the support on its
  top antidiagonal). In particular no `h_{q,s}` on `ℓ²(ℤ²)` has one (`no_compact_eigenvector_hqs`).
* On a finite graph with link signs `s` (binary fluxes), `K_s = -J ∑_e s_e (|i⟩⟨j| + |j⟩⟨i|)` and
  `h_{q,s} = K_s + 2h diag(q)`. For the checkerboard state on a bipartite graph the charges are
  `q⁰ = -ε` (sublattice signs), `h_{q⁰,s} = K_s - 2h Γ_ε` with `K_s Γ_ε = -Γ_ε K_s`, hence
  `h_{q⁰,s}² = K_s² + 4h²`: every eigenvalue has modulus at least `2|h|`, and the eigenspace of
  `h²` at `4h²` (the band edges `±2|h|`) is `ker K_s`.
-/

open Matrix Finset

noncomputable section

namespace DFL

/-! ### No compactly supported eigenvectors on `ℤ²` -/

/-- Nearest neighbours in `ℤ²`. -/
def IsNN (x y : ℤ × ℤ) : Prop := (x.1 - y.1).natAbs + (x.2 - y.2).natAbs = 1

/-- **No compactly supported eigenvectors on `ℤ²`** (paper: Proposition `prop:nocompact`): if
`(Hu)(x) = ∑_{y ∼ x} a_{xy} u(y) + V(x) u(x)` with all nearest-neighbour `a_{xy} ≠ 0`, a finitely
supported solution of `Hu = Eu` vanishes. -/
theorem no_compact_eigenvector (a : ℤ × ℤ → ℤ × ℤ → ℂ) (V : ℤ × ℤ → ℂ) (E : ℂ)
    (ha : ∀ x y, IsNN x y → a x y ≠ 0) (u : ℤ × ℤ → ℂ) (hu : (Function.support u).Finite)
    (heig : ∀ x, a x (x + (1, 0)) * u (x + (1, 0)) + a x (x - (1, 0)) * u (x - (1, 0)) +
        a x (x + (0, 1)) * u (x + (0, 1)) + a x (x - (0, 1)) * u (x - (0, 1)) + V x * u x =
      E * u x) :
    u = 0 := by
  by_contra hne
  have memS : ∀ w, u w ≠ 0 → w ∈ hu.toFinset := fun w hw => hu.mem_toFinset.mpr hw
  have hSne : hu.toFinset.Nonempty := by
    obtain ⟨w, hw⟩ := Function.ne_iff.mp hne
    exact ⟨w, memS w hw⟩
  -- the top antidiagonal `z₁ + z₂ = m` of the support
  obtain ⟨z0, hz0S, hz0⟩ := hu.toFinset.exists_max_image (fun z : ℤ × ℤ => z.1 + z.2) hSne
  have hTne : (hu.toFinset.filter (fun z : ℤ × ℤ => z.1 + z.2 = z0.1 + z0.2)).Nonempty :=
    ⟨z0, Finset.mem_filter.mpr ⟨hz0S, rfl⟩⟩
  -- its point with the largest first coordinate
  obtain ⟨z, hzT, hz⟩ := (hu.toFinset.filter
    (fun z : ℤ × ℤ => z.1 + z.2 = z0.1 + z0.2)).exists_max_image (fun z : ℤ × ℤ => z.1) hTne
  obtain ⟨hzS, hzm⟩ := Finset.mem_filter.mp hzT
  have out1 : ∀ w : ℤ × ℤ, z0.1 + z0.2 < w.1 + w.2 → u w = 0 := by
    intro w hw
    by_contra hw0
    have := hz0 w (memS w hw0)
    omega
  have out2 : ∀ w : ℤ × ℤ, w.1 + w.2 = z0.1 + z0.2 → z.1 < w.1 → u w = 0 := by
    intro w hw1 hw2
    by_contra hw0
    have := hz w (Finset.mem_filter.mpr ⟨memS w hw0, hw1⟩)
    omega
  have h := heig (z + (1, 0))
  have e1 : u (z + (1, 0) + (1, 0)) = 0 := out1 _ (by simp; omega)
  have e2 : u (z + (1, 0) + (0, 1)) = 0 := out1 _ (by simp; omega)
  have e3 : u (z + (1, 0) - (0, 1)) = 0 := out2 _ (by simp; omega) (by simp)
  have e4 : u (z + (1, 0)) = 0 := out1 _ (by simp; omega)
  rw [e1, e2, e3, e4, add_sub_cancel_right] at h
  simp only [mul_zero, add_zero, zero_add] at h
  have hnn : IsNN (z + (1, 0)) z := by simp [IsNN]
  exact mul_ne_zero (ha _ _ hnn) (hu.mem_toFinset.mp hzS) h

/-- No `h_{q,s} = K_s + 2h diag(q)` on `ℓ²(ℤ²)` (signs `s_{xy} = ±1`, `J ≠ 0`, any potential `q`) has
a finitely supported eigenvector. -/
theorem no_compact_eigenvector_hqs (J h : ℝ) (hJ : J ≠ 0) (s : ℤ × ℤ → ℤ × ℤ → ℝ)
    (hs : ∀ x y, s x y = 1 ∨ s x y = -1) (q : ℤ × ℤ → ℝ) (E : ℂ) (u : ℤ × ℤ → ℂ)
    (hu : (Function.support u).Finite)
    (heig : ∀ x, ((-J : ℝ) : ℂ) * ((s x (x + (1, 0)) : ℂ) * u (x + (1, 0)) +
        (s x (x - (1, 0)) : ℂ) * u (x - (1, 0)) + (s x (x + (0, 1)) : ℂ) * u (x + (0, 1)) +
        (s x (x - (0, 1)) : ℂ) * u (x - (0, 1))) + ((2 * h * q x : ℝ) : ℂ) * u x = E * u x) :
    u = 0 := by
  refine no_compact_eigenvector (fun x y => ((-J : ℝ) : ℂ) * (s x y : ℂ))
    (fun x => ((2 * h * q x : ℝ) : ℂ)) E ?_ u hu ?_
  · intro x y _
    have hJ' : ((-J : ℝ) : ℂ) ≠ 0 := by exact_mod_cast neg_ne_zero.mpr hJ
    have hs' : (s x y : ℂ) ≠ 0 := by
      rcases hs x y with h1 | h1 <;> simp [h1]
    exact mul_ne_zero hJ' hs'
  · intro x
    linear_combination heig x

/-! ### Auxiliary lemmas: real symmetric matrices and their complexifications -/

namespace TwoDAux

/-- `(K + c•P)² = K² + c•(KP + PK) + c²•P²` for square real matrices. -/
lemma sq_add_smul {n : Type*} [Fintype n] [DecidableEq n] (K P : Matrix n n ℝ) (c : ℝ) :
    (K + c • P) * (K + c • P) = K * K + c • (K * P + P * K) + (c * c) • (P * P) := by
  simp only [add_mul, mul_add, smul_mul_assoc, mul_smul_comm, smul_smul, smul_add]
  abel

/-- A real symmetric matrix, viewed as a complex matrix, is self-adjoint. -/
lemma conjTranspose_map_ofReal {n : Type*} (A : Matrix n n ℝ) (hA : Aᵀ = A) :
    (A.map (fun r : ℝ => (r : ℂ)))ᴴ = A.map (fun r : ℝ => (r : ℂ)) := by
  ext i j
  have h := congrFun (congrFun hA j) i
  simp only [Matrix.transpose_apply] at h
  simp [Matrix.conjTranspose_apply, h]

/-- Complexification of real matrices is multiplicative. -/
lemma map_ofReal_mul {n : Type*} [Fintype n] (A B : Matrix n n ℝ) :
    (A * B).map (fun r : ℝ => (r : ℂ)) =
      A.map (fun r : ℝ => (r : ℂ)) * B.map (fun r : ℝ => (r : ℂ)) :=
  Matrix.map_mul (f := Complex.ofRealHom)

/-- Complexification of a real multiple of the identity. -/
lemma map_ofReal_smul_one {n : Type*} [DecidableEq n] (c : ℝ) :
    (c • (1 : Matrix n n ℝ)).map (fun r : ℝ => (r : ℂ)) = (c : ℂ) • (1 : Matrix n n ℂ) := by
  ext i j
  by_cases hij : i = j <;> simp [hij]

/-- For a self-adjoint `H`: `⟨Hv, Hv⟩ = ⟨v, H²v⟩`. -/
lemma star_mulVec_dotProduct {n : Type*} [Fintype n] (H : Matrix n n ℂ) (hH : Hᴴ = H)
    (v : n → ℂ) : star (H *ᵥ v) ⬝ᵥ (H *ᵥ v) = star v ⬝ᵥ ((H * H) *ᵥ v) := by
  rw [Matrix.star_mulVec, ← Matrix.dotProduct_mulVec, hH, Matrix.mulVec_mulVec]

/-- `⟨v, v⟩ = ∑ᵢ ‖vᵢ‖²`. -/
lemma star_dotProduct_self {n : Type*} [Fintype n] (v : n → ℂ) :
    star v ⬝ᵥ v = ((∑ i, ‖v i‖ ^ 2 : ℝ) : ℂ) := by
  simp only [dotProduct, Pi.star_apply, Complex.star_def, Complex.conj_mul']
  push_cast
  rfl

/-- For a real symmetric `K`: `K² v = 0 ↔ K v = 0`. -/
lemma sq_mulVec_eq_zero {n : Type*} [Fintype n] (K : Matrix n n ℝ) (hK : Kᵀ = K) (v : n → ℝ) :
    (K * K) *ᵥ v = 0 ↔ K *ᵥ v = 0 := by
  constructor
  · intro h0
    have : (K *ᵥ v) ⬝ᵥ (K *ᵥ v) = 0 := by
      calc (K *ᵥ v) ⬝ᵥ (K *ᵥ v) = (v ᵥ* Kᵀ) ⬝ᵥ (K *ᵥ v) := by rw [Matrix.vecMul_transpose]
        _ = v ⬝ᵥ (Kᵀ *ᵥ (K *ᵥ v)) := (Matrix.dotProduct_mulVec _ _ _).symm
        _ = 0 := by rw [hK, Matrix.mulVec_mulVec, h0, dotProduct_zero]
    exact dotProduct_self_eq_zero.mp this
  · intro h0
    rw [← Matrix.mulVec_mulVec, h0, Matrix.mulVec_zero]

end TwoDAux

/-! ### Binary fluxes on a finite graph; the checkerboard -/

namespace FinGraph

variable {N M : ℕ} (G : FinGraph N M)

/-- The hopping operator with link signs `s`: `K_s = -J ∑_e s_e (|s(e)⟩⟨t(e)| + |t(e)⟩⟨s(e)|)`. -/
def Ks (J : ℝ) (s : Fin M → ℝ) : Matrix (Fin N) (Fin N) ℝ :=
  fun i j => -J * ∑ e ∈ univ.filter (fun e => (G.src e = i ∧ G.tgt e = j) ∨
    (G.tgt e = i ∧ G.src e = j)), s e

/-- `h_{q,s} = K_s + 2h diag(q)`. -/
def hqs (J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) : Matrix (Fin N) (Fin N) ℝ :=
  G.Ks J s + (2 * h) • diagonal q

/-- Entries of `K_s`. -/
private lemma Ks_apply (J : ℝ) (s : Fin M → ℝ) (i j : Fin N) :
    G.Ks J s i j = -J * ∑ e ∈ univ.filter (fun e => (G.src e = i ∧ G.tgt e = j) ∨
      (G.tgt e = i ∧ G.src e = j)), s e := rfl

/-- `K_s` is symmetric. -/
private lemma Ks_transpose (J : ℝ) (s : Fin M → ℝ) : (G.Ks J s)ᵀ = G.Ks J s := by
  ext i j
  simp only [Matrix.transpose_apply, Ks]
  congr 2
  ext e
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]
  tauto

/-- `h_{q,s}` is symmetric. -/
private lemma hqs_transpose (J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) :
    (G.hqs J h q s)ᵀ = G.hqs J h q s := by
  rw [hqs, Matrix.transpose_add, Matrix.transpose_smul, Matrix.diagonal_transpose,
    Ks_transpose]

/-- With all signs `+1`, `h_{q,+}` is the Anderson–Bernoulli operator `h_q` of the trivial flux
sector. -/
theorem hqs_one (J h : ℝ) (q : Cfg N) :
    G.hqs J h (fun i => sgn (q i)) (fun _ => 1) = G.hq J h q := by
  have hK : G.Ks J (fun _ => 1) = (-J) • G.adj := by
    ext i j
    simp [Ks, adj]
  rw [hqs, hq, hK]

/-- **Gauge covariance**: `D_θ h_{q,s} D_θ = h_{q,θ·s}`, `(θ·s)_e = θ_{s(e)} s_e θ_{t(e)}`, for
`θ ∈ {±1}^V`. -/
theorem gauge_hqs (J h : ℝ) (q : Fin N → ℝ) (s : Fin M → ℝ) (θ : Fin N → ℝ)
    (hθ : ∀ i, θ i = 1 ∨ θ i = -1) :
    diagonal θ * G.hqs J h q s * diagonal θ =
      G.hqs J h q (fun e => θ (G.src e) * s e * θ (G.tgt e)) := by
  ext i j
  simp only [hqs, Matrix.diagonal_mul, Matrix.mul_diagonal, Matrix.add_apply, Matrix.smul_apply,
    smul_eq_mul, Matrix.diagonal_apply, Ks_apply]
  have hK : θ i * (-J * ∑ e ∈ univ.filter (fun e => (G.src e = i ∧ G.tgt e = j) ∨
      (G.tgt e = i ∧ G.src e = j)), s e) * θ j =
      -J * ∑ e ∈ univ.filter (fun e => (G.src e = i ∧ G.tgt e = j) ∨
      (G.tgt e = i ∧ G.src e = j)), θ (G.src e) * s e * θ (G.tgt e) := by
    simp only [Finset.mul_sum, Finset.sum_mul]
    refine Finset.sum_congr rfl fun e he => ?_
    rcases (Finset.mem_filter.mp he).2 with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · rw [h1, h2]; ring
    · rw [h1, h2]; ring
  have hD : θ i * (2 * h * (if i = j then q i else 0)) * θ j =
      2 * h * (if i = j then q i else 0) := by
    split_ifs with hij
    · subst hij
      rcases hθ i with h1 | h1 <;> rw [h1] <;> ring
    · ring
  linear_combination hK + hD

/-- On a bipartite graph (sublattice signs `ε`, every edge joins opposite signs), `K_s`
anticommutes with `Γ_ε = diag(ε)`. -/
theorem Ks_anticomm (J : ℝ) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε : ∀ e, ε (G.src e) = -ε (G.tgt e)) :
    G.Ks J s * diagonal ε + diagonal ε * G.Ks J s = 0 := by
  ext i j
  simp only [Matrix.add_apply, Matrix.mul_diagonal, Matrix.diagonal_mul, Matrix.zero_apply,
    Ks_apply]
  have h0 : (∑ e ∈ univ.filter (fun e => (G.src e = i ∧ G.tgt e = j) ∨
      (G.tgt e = i ∧ G.src e = j)), s e) * (ε i + ε j) = 0 := by
    rw [Finset.sum_mul]
    refine Finset.sum_eq_zero fun e he => ?_
    have := hε e
    rcases (Finset.mem_filter.mp he).2 with ⟨h1, h2⟩ | ⟨h1, h2⟩
    · rw [h1, h2] at this; rw [this]; ring
    · rw [h1, h2] at this; rw [this]; ring
  linear_combination (-J) * h0

/-- **Checkerboard identity** (paper: Proposition `prop:2dlaws`(c)(i)): for the checkerboard charges
`q⁰ = -ε`, `h_{q⁰,s}² = K_s² + 4h²`. -/
theorem checkerboard_sq (J h : ℝ) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε1 : ∀ i, ε i = 1 ∨ ε i = -1) (hε : ∀ e, ε (G.src e) = -ε (G.tgt e)) :
    G.hqs J h (fun i => -ε i) s * G.hqs J h (fun i => -ε i) s =
      G.Ks J s * G.Ks J s + (4 * h ^ 2) • (1 : Matrix (Fin N) (Fin N) ℝ) := by
  have hA := G.Ks_anticomm J s ε hε
  have hD : diagonal ε * diagonal ε = (1 : Matrix (Fin N) (Fin N) ℝ) := by
    rw [Matrix.diagonal_mul_diagonal, ← Matrix.diagonal_one]
    congr 1
    funext i
    rcases hε1 i with h1 | h1 <;> rw [h1] <;> norm_num
  have hP : diagonal (fun i => -ε i) = -diagonal ε := (Matrix.diagonal_neg ε).symm
  have hPP : (-diagonal ε) * (-diagonal ε) = (1 : Matrix (Fin N) (Fin N) ℝ) := by
    rw [neg_mul_neg, hD]
  have hKP : G.Ks J s * (-diagonal ε) + (-diagonal ε) * G.Ks J s = 0 := by
    rw [Matrix.mul_neg, Matrix.neg_mul, ← neg_add, hA, neg_zero]
  rw [hqs, hP, TwoDAux.sq_add_smul, hKP, hPP, smul_zero, add_zero,
    show 2 * h * (2 * h) = 4 * h ^ 2 by ring]

/-- **Chiral gap** (paper: Proposition `prop:2dlaws`(c)(ii)): every eigenvalue `E` of the checkerboard
operator `h_{q⁰,s}` (real signs `s`) has `|E| ≥ 2|h|`. -/
theorem checkerboard_gap (J h : ℝ) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε1 : ∀ i, ε i = 1 ∨ ε i = -1) (hε : ∀ e, ε (G.src e) = -ε (G.tgt e))
    (E : ℂ) (v : Fin N → ℂ) (hv : v ≠ 0)
    (heig : (G.hqs J h (fun i => -ε i) s).map (fun r => (r : ℂ)) *ᵥ v = E • v) :
    2 * |h| ≤ ‖E‖ := by
  have hsq := G.checkerboard_sq J h s ε hε1 hε
  have hAH := TwoDAux.conjTranspose_map_ofReal _ (hqs_transpose G J h (fun i => -ε i) s)
  have hKH := TwoDAux.conjTranspose_map_ofReal _ (Ks_transpose G J s)
  have h1 := TwoDAux.star_mulVec_dotProduct _ hAH v
  have h2 := TwoDAux.star_mulVec_dotProduct _ hKH v
  rw [heig, ← TwoDAux.map_ofReal_mul, hsq, Matrix.map_add _ (fun a b => Complex.ofReal_add a b),
    TwoDAux.map_ofReal_mul, TwoDAux.map_ofReal_smul_one, Matrix.add_mulVec, dotProduct_add, ← h2,
    Matrix.smul_mulVec, Matrix.one_mulVec, dotProduct_smul, star_smul, smul_dotProduct,
    dotProduct_smul, TwoDAux.star_dotProduct_self, TwoDAux.star_dotProduct_self, smul_eq_mul,
    smul_eq_mul, smul_eq_mul, ← mul_assoc, Complex.star_def, Complex.mul_conj'] at h1
  -- `‖E‖² ‖v‖² = ‖K v‖² + 4h² ‖v‖²`, with `‖v‖² > 0`
  have h3 : ‖E‖ ^ 2 * (∑ i, ‖v i‖ ^ 2) =
      (∑ i, ‖(((G.Ks J s).map fun r => (r : ℂ)) *ᵥ v) i‖ ^ 2) +
        4 * h ^ 2 * (∑ i, ‖v i‖ ^ 2) := by
    exact_mod_cast h1
  have hSpos : 0 < ∑ i, ‖v i‖ ^ 2 := by
    obtain ⟨i, hi⟩ := Function.ne_iff.mp hv
    exact Finset.sum_pos' (fun j _ => by positivity)
      ⟨i, Finset.mem_univ _, pow_pos (norm_pos_iff.mpr hi) 2⟩
  have hT : 0 ≤ ∑ i, ‖(((G.Ks J s).map fun r => (r : ℂ)) *ᵥ v) i‖ ^ 2 :=
    Finset.sum_nonneg (fun j _ => by positivity)
  have h4 : 4 * h ^ 2 ≤ ‖E‖ ^ 2 := by
    by_contra hc
    have := mul_lt_mul_of_pos_right (not_le.mp hc) hSpos
    linarith
  have h5 : (2 * |h|) ^ 2 ≤ ‖E‖ ^ 2 := by rw [mul_pow, sq_abs]; linarith
  nlinarith [norm_nonneg E, abs_nonneg h]

/-- **Band edges = chiral zero modes** (paper: Proposition `prop:2dlaws`(c)(ii)): the eigenspace of
`h_{q⁰,s}²` at `4h²` (the spectral subspace of `|h_{q⁰,s}|` at `2|h|`) is `ker K_s`. -/
theorem checkerboard_edge_iff (J h : ℝ) (s : Fin M → ℝ) (ε : Fin N → ℝ)
    (hε1 : ∀ i, ε i = 1 ∨ ε i = -1) (hε : ∀ e, ε (G.src e) = -ε (G.tgt e))
    (v : Fin N → ℝ) :
    G.Ks J s *ᵥ v = 0 ↔
      (G.hqs J h (fun i => -ε i) s * G.hqs J h (fun i => -ε i) s) *ᵥ v = (4 * h ^ 2) • v := by
  rw [G.checkerboard_sq J h s ε hε1 hε, Matrix.add_mulVec, Matrix.smul_mulVec, Matrix.one_mulVec,
    add_eq_right, TwoDAux.sq_mulVec_eq_zero _ (Ks_transpose G J s)]

end FinGraph

end DFL
