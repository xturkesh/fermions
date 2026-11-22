import DFL.CutLocal
import DFL.Main

/-!
# The number entropy is a lower bound to the entanglement entropy
(paper: Proposition `prop:entropies`(c))

If a density matrix commutes with a "charge" (here the particle number of the region), its von
Neumann entropy is at least the Shannon entropy of the charge distribution:
`S(ρ) = S_N + ∑_k p_k S(ρ^{(k)}) ≥ S_N`. We prove the inequality directly (each eigenvalue with
weight in the charge block `k` is at most the block weight `p_k`) and apply it to
* the entanglement of an evolved Fock state (Slater determinant) in every charge sector:
  `S(ρ^q_A(t)) ≥ H(p^q_t)`;
* the full matter + gauge entanglement of the LGT state started from `|↑…↑⟩_z ⊗ |n⟩`:
  `S_full(t) ≥ S_N(t)`.

Proof outline (auxiliary lemmas in `DFL.TaxonomyAux`):
* with `ρ = ∑_i λ_i v_i v_i†` and `w_{ik} = ‖P_k v_i‖²`, the block structure makes `P_k v_i` an
  eigenvector of `ρ`, so `λ_i w_{ik} = ⟨P_k v_i, ρ P_k v_i⟩ ≤ p_k w_{ik}` by Cauchy–Schwarz
  (`key_ineq`); hence `S(ρ) - H(p) = ∑_{i,k} λ_i w_{ik} (log p_k - log λ_i) ≥ 0` (`block_bound`;
  the normalization `Tr ρ = 1` is not needed);
* quadratic Hamiltonians and the LGT Hamiltonian conserve the (matter) particle number
  (`evo_eq_zero_of_conserve`), so the evolved states have a definite particle number and the
  reduced density matrices are block diagonal in the particle number of the region.
-/

open Matrix Finset
open scoped ComplexOrder

noncomputable section

namespace DFL

variable {N : ℕ}

namespace TaxonomyAux

/-! ### The block entropy bound -/

section Block

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- Spectral decomposition of a Hermitian matrix, entrywise: `ρ_{ab} = ∑_j U_{aj} λ_j conj U_{bj}`. -/
lemma spectral_apply {ρ : Matrix ι ι ℂ} (hH : ρ.IsHermitian) (a b : ι) :
    ρ a b = ∑ j, (hH.eigenvectorUnitary : Matrix ι ι ℂ) a j * (hH.eigenvalues j : ℂ) *
      star ((hH.eigenvectorUnitary : Matrix ι ι ℂ) b j) := by
  conv_lhs => rw [hH.spectral_theorem]
  rw [Unitary.conjStarAlgAut_apply, Matrix.mul_apply]
  simp only [Matrix.mul_diagonal, Matrix.star_apply, Function.comp_apply]
  rfl

/-- The columns of the eigenvector unitary are orthonormal. -/
lemma col_orth {ρ : Matrix ι ι ℂ} (hH : ρ.IsHermitian) (i j : ι) :
    ∑ a, star ((hH.eigenvectorUnitary : Matrix ι ι ℂ) a i) *
      (hH.eigenvectorUnitary : Matrix ι ι ℂ) a j = if i = j then 1 else 0 := by
  have h := congrFun (congrFun (Unitary.coe_star_mul_self hH.eigenvectorUnitary) i) j
  rw [Matrix.mul_apply] at h
  simpa [Matrix.star_apply, Matrix.one_apply] using h

/-- `conj z * z = |z|²`. -/
lemma star_mul_self_eq (z : ℂ) : star z * z = (Complex.normSq z : ℂ) := by
  rw [Complex.normSq_eq_conj_mul_self, Complex.star_def]

/-- Cauchy–Schwarz for complex sums over a finset. -/
lemma normSq_sum_le {κ : Type*} (s : Finset κ) (f g : κ → ℂ) :
    Complex.normSq (∑ b ∈ s, star (f b) * g b) ≤
      (∑ b ∈ s, Complex.normSq (f b)) * ∑ b ∈ s, Complex.normSq (g b) := by
  have h1 : ‖∑ b ∈ s, star (f b) * g b‖ ≤ ∑ b ∈ s, ‖f b‖ * ‖g b‖ := by
    refine (norm_sum_le _ _).trans (le_of_eq ?_)
    refine Finset.sum_congr rfl fun b _ => ?_
    rw [norm_mul, norm_star]
  have h2 : (∑ b ∈ s, ‖f b‖ * ‖g b‖) ^ 2 ≤
      (∑ b ∈ s, Complex.normSq (f b)) * ∑ b ∈ s, Complex.normSq (g b) := by
    refine Finset.sum_sq_le_sum_mul_sum_of_sq_le_mul s (fun b _ => Complex.normSq_nonneg _)
      (fun b _ => Complex.normSq_nonneg _) (fun b _ => ?_)
    rw [Complex.normSq_eq_norm_sq, Complex.normSq_eq_norm_sq, mul_pow]
  rw [Complex.normSq_eq_norm_sq]
  exact (pow_le_pow_left₀ (norm_nonneg _) h1 2).trans h2

/-- **Key inequality**: with `U` the eigenvector unitary of a block-diagonal positive semidefinite
`ρ` and `w_{ik}` the weight of the `i`-th eigenvector in the block `k`,
`λ_i w_{ik} = ⟨P_k v_i, ρ P_k v_i⟩ ≤ (∑_j λ_j w_{jk}) w_{ik}`. -/
lemma key_ineq (ρ : Matrix ι ι ℂ) (hρ : ρ.PosSemidef) (c : ι → ℕ)
    (hblock : ∀ a b, c a ≠ c b → ρ a b = 0) (i : ι) (k : ℕ) :
    hρ.isHermitian.eigenvalues i * ∑ a ∈ univ.filter (fun a => c a = k),
        Complex.normSq ((hρ.isHermitian.eigenvectorUnitary : Matrix ι ι ℂ) a i) ≤
      (∑ j, hρ.isHermitian.eigenvalues j * ∑ a ∈ univ.filter (fun a => c a = k),
        Complex.normSq ((hρ.isHermitian.eigenvectorUnitary : Matrix ι ι ℂ) a j)) *
      ∑ a ∈ univ.filter (fun a => c a = k),
        Complex.normSq ((hρ.isHermitian.eigenvectorUnitary : Matrix ι ι ℂ) a i) := by
  set hH := hρ.isHermitian with hHdef
  set U : Matrix ι ι ℂ := (hH.eigenvectorUnitary : Matrix ι ι ℂ) with hU
  set lam := hH.eigenvalues with hlam
  set F := univ.filter (fun a => c a = k) with hF
  have hlam0 : ∀ j, 0 ≤ lam j := fun j => hρ.eigenvalues_nonneg j
  have hspec : ∀ a b, ρ a b = ∑ j, U a j * (lam j : ℂ) * star (U b j) := spectral_apply hH
  have horth : ∀ i j, ∑ a, star (U a i) * U a j = if i = j then 1 else 0 := col_orth hH
  set S : ℂ := ∑ a ∈ F, ∑ b ∈ F, star (U a i) * ρ a b * U b i with hS
  -- (1) `P_k v_i` is an eigenvector of `ρ` with eigenvalue `λ_i`
  have hinner : ∀ a ∈ F, ∑ b ∈ F, ρ a b * U b i = (lam i : ℂ) * U a i := by
    intro a ha
    have hca : c a = k := (Finset.mem_filter.1 ha).2
    rw [Finset.sum_subset (Finset.subset_univ F) fun b _ hb => ?_]
    · simp only [hspec, Finset.sum_mul]
      rw [Finset.sum_comm]
      have : ∀ j, ∑ b, U a j * (lam j : ℂ) * star (U b j) * U b i =
          U a j * (lam j : ℂ) * (if j = i then 1 else 0) := by
        intro j
        rw [← horth j i, Finset.mul_sum]
        exact Finset.sum_congr rfl fun b _ => by ring
      simp only [this, mul_ite, mul_one, mul_zero, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
      ring
    · have hcb : c b ≠ k := fun h => hb (Finset.mem_filter.2 ⟨Finset.mem_univ _, h⟩)
      rw [hblock a b (by rw [hca]; exact fun h => hcb h.symm), zero_mul]
  have h1 : S = (lam i : ℂ) * ((∑ a ∈ F, Complex.normSq (U a i) : ℝ) : ℂ) := by
    rw [hS]
    have : ∀ a ∈ F, ∑ b ∈ F, star (U a i) * ρ a b * U b i =
        (lam i : ℂ) * (Complex.normSq (U a i) : ℂ) := by
      intro a ha
      simp only [mul_assoc, ← Finset.mul_sum]
      rw [hinner a ha, ← star_mul_self_eq]
      ring
    rw [Finset.sum_congr rfl this, ← Finset.mul_sum, Complex.ofReal_sum]
  -- (2) the same quantity through the spectral decomposition: `∑_j λ_j |⟨P_k v_j, P_k v_i⟩|²`
  set z : ι → ℂ := fun j => ∑ b ∈ F, star (U b j) * U b i with hz
  have h2 : S = ∑ j, (lam j : ℂ) * (Complex.normSq (z j) : ℂ) := by
    have hz2 : ∀ j, (lam j : ℂ) * (Complex.normSq (z j) : ℂ) =
        ∑ a ∈ F, ∑ b ∈ F, star (U a i) * (U a j * (lam j : ℂ) * star (U b j)) * U b i := by
      intro j
      rw [← star_mul_self_eq, hz, star_sum, Finset.sum_mul_sum, Finset.mul_sum]
      refine Finset.sum_congr rfl fun a _ => ?_
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl fun b _ => ?_
      rw [star_mul', star_star]
      ring
    simp only [hz2]
    rw [hS]
    simp only [hspec, Finset.mul_sum, Finset.sum_mul]
    calc _ = ∑ a ∈ F, ∑ j, ∑ b ∈ F,
          star (U a i) * (U a j * (lam j : ℂ) * star (U b j)) * U b i :=
          Finset.sum_congr rfl fun a _ => Finset.sum_comm
      _ = _ := Finset.sum_comm
  -- (3) Cauchy–Schwarz
  have h3 : ∀ j, Complex.normSq (z j) ≤
      (∑ a ∈ F, Complex.normSq (U a j)) * ∑ a ∈ F, Complex.normSq (U a i) :=
    fun j => normSq_sum_le F (fun b => U b j) (fun b => U b i)
  have hre1 : S.re = lam i * ∑ a ∈ F, Complex.normSq (U a i) := by
    rw [h1, ← Complex.ofReal_mul, Complex.ofReal_re]
  have hre2 : S.re = ∑ j, lam j * Complex.normSq (z j) := by
    rw [h2, Complex.re_sum]
    exact Finset.sum_congr rfl fun j _ => by rw [← Complex.ofReal_mul, Complex.ofReal_re]
  rw [← hre1, hre2, Finset.sum_mul]
  refine Finset.sum_le_sum fun j _ => ?_
  rw [mul_assoc]
  exact mul_le_mul_of_nonneg_left (h3 j) (hlam0 j)

/-- **Block entropy bound** (no normalization needed): for a positive semidefinite `ρ` which is
block diagonal with respect to `c`, `H(p) ≤ S(ρ)` with `p_k = ∑_{c a = k} ρ_{aa}`. -/
lemma block_bound (ρ : Matrix ι ι ℂ) (hρ : ρ.PosSemidef) (c : ι → ℕ) (M : ℕ)
    (hcM : ∀ a, c a ≤ M) (hblock : ∀ a b, c a ≠ c b → ρ a b = 0) :
    shannon M (fun k => (∑ a ∈ univ.filter (fun a => c a = k), ρ a a).re) ≤
      vonNeumann ρ hρ.isHermitian := by
  have hkey := key_ineq ρ hρ c hblock
  set hH := hρ.isHermitian with hHdef
  set U : Matrix ι ι ℂ := (hH.eigenvectorUnitary : Matrix ι ι ℂ) with hU
  set lam := hH.eigenvalues with hlam
  have hlam0 : ∀ j, 0 ≤ lam j := fun j => hρ.eigenvalues_nonneg j
  have hspec : ∀ a b, ρ a b = ∑ j, U a j * (lam j : ℂ) * star (U b j) := spectral_apply hH
  have horth : ∀ i j, ∑ a, star (U a i) * U a j = if i = j then 1 else 0 := col_orth hH
  -- the weight `w_{ik} = ‖P_k v_i‖²` of the `i`-th eigenvector in the block `k`
  set w : ι → ℕ → ℝ := fun i k =>
    ∑ a ∈ univ.filter (fun a => c a = k), Complex.normSq (U a i) with hw
  have hw0 : ∀ i k, 0 ≤ w i k := fun i k => Finset.sum_nonneg fun a _ => Complex.normSq_nonneg _
  have hw1 : ∀ i, ∑ k ∈ range (M + 1), w i k = 1 := by
    intro i
    simp only [hw]
    rw [Finset.sum_fiberwise_of_maps_to
      (fun a _ => Finset.mem_range.2 (Nat.lt_succ_of_le (hcM a)))]
    have := horth i i
    simp only [ite_true, star_mul_self_eq] at this
    exact_mod_cast this
  -- the block weights `p_k = ∑_i λ_i w_{ik}`
  have hp : ∀ k, (∑ a ∈ univ.filter (fun a => c a = k), ρ a a).re = ∑ j, lam j * w j k := by
    intro k
    simp only [hspec, Complex.re_sum, hw, Finset.mul_sum]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun a _ => ?_
    have : U a j * (lam j : ℂ) * star (U a j) = ((lam j * Complex.normSq (U a j) : ℝ) : ℂ) := by
      rw [Complex.ofReal_mul, ← star_mul_self_eq]; ring
    rw [this, Complex.ofReal_re]
  -- termwise comparison: `λ_i ≤ p_k` whenever `λ_i w_{ik} > 0`
  have hterm : ∀ i k, lam i * w i k * (-Real.log (∑ j, lam j * w j k)) ≤
      lam i * w i k * (-Real.log (lam i)) := by
    intro i k
    rcases (hlam0 i).eq_or_lt with h0 | hpos
    · rw [← h0]; simp
    rcases (hw0 i k).eq_or_lt with h0' | hwpos
    · rw [← h0']; simp
    have hle : lam i ≤ ∑ j, lam j * w j k := le_of_mul_le_mul_right (hkey i k) hwpos
    have hlog := Real.log_le_log hpos hle
    have hnn : 0 ≤ lam i * w i k := mul_nonneg hpos.le hwpos.le
    nlinarith
  unfold shannon vonNeumann
  simp only [hp]
  calc ∑ k ∈ range (M + 1), Real.negMulLog (∑ j, lam j * w j k)
      = ∑ k ∈ range (M + 1), ∑ i, lam i * w i k * (-Real.log (∑ j, lam j * w j k)) := by
        refine Finset.sum_congr rfl fun k _ => ?_
        rw [Real.negMulLog, ← Finset.sum_mul]
        ring
    _ ≤ ∑ k ∈ range (M + 1), ∑ i, lam i * w i k * (-Real.log (lam i)) :=
        Finset.sum_le_sum fun k _ => Finset.sum_le_sum fun i _ => hterm i k
    _ = ∑ i, Real.negMulLog (lam i) := by
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl fun i _ => ?_
        rw [Real.negMulLog]
        calc ∑ k ∈ range (M + 1), lam i * w i k * (-Real.log (lam i))
            = lam i * (∑ k ∈ range (M + 1), w i k) * (-Real.log (lam i)) := by
              rw [Finset.mul_sum, Finset.sum_mul]
          _ = -lam i * Real.log (lam i) := by rw [hw1 i]; ring

end Block

/-! ### Particle-number conservation -/

section Conservation

/-- If a real generator `A` conserves a real charge `f` (`A_{ab} = 0` unless `f a = f b`), so does
the propagator `e^{-itA}`. -/
lemma evo_eq_zero_of_conserve {κ : Type*} [Fintype κ] [DecidableEq κ] (A : Matrix κ κ ℝ)
    (f : κ → ℝ) (hA : ∀ a b, f a ≠ f b → A a b = 0) (t : ℝ) (a b : κ) (hab : f a ≠ f b) :
    evo A t a b = 0 := by
  set X : Matrix κ κ ℂ := (-(t : ℂ) * Complex.I) • cplx A with hX
  set d : κ → ℂ := fun m => (f m : ℂ) with hd
  have hcomm : Commute X (diagonal d) := by
    show X * diagonal d = diagonal d * X
    ext p q
    rw [mul_diagonal, diagonal_mul]
    by_cases hpq : f p = f q
    · simp only [hd, hpq, mul_comm]
    · simp [hX, cplx, hA p q hpq]
  have h2 := congrFun (congrFun hcomm.exp_left.eq a) b
  rw [mul_diagonal, diagonal_mul] at h2
  have h3 : NormedSpace.exp X a b * (d b - d a) = 0 := by rw [mul_sub, h2]; ring
  rcases mul_eq_zero.mp h3 with h4 | h4
  · exact h4
  · have h5 : (f a : ℂ) = f b := (sub_eq_zero.mp h4).symm
    exact absurd (by exact_mod_cast h5) hab

/-- Emptying an occupied mode lowers the particle number by one. -/
lemma count_univ_update (m : Cfg N) (y : Fin N) (hy : m y = 1) :
    count univ m = count univ (Function.update m y 0) + 1 := by
  unfold count
  rw [← Finset.add_sum_erase _ _ (mem_univ y), ← Finset.add_sum_erase _ _ (mem_univ y),
    Function.update_self, hy]
  have : ∑ k ∈ univ.erase y, occ (Function.update m y 0 k) = ∑ k ∈ univ.erase y, occ (m k) :=
    Finset.sum_congr rfl fun k hk => by rw [Function.update_of_ne (Finset.ne_of_mem_erase hk)]
  rw [this]
  have h1 : occ (1 : ZMod 2) = 1 := by simp [occ]
  have h0 : occ (0 : ZMod 2) = 0 := by simp [occ]
  rw [h1, h0]
  ring

/-- `c_x† c_y` conserves the particle number. -/
lemma cre_mul_ann_eq_zero_of_count (x y : Fin N) (m' m : Cfg N)
    (h : count univ m' ≠ count univ m) : (cre x * ann y) m' m = 0 := by
  simp only [cre, Matrix.mul_apply, Matrix.transpose_apply, ann]
  refine Finset.sum_eq_zero fun p _ => ?_
  by_cases h1 : m' x = 1 ∧ p = Function.update m' x 0
  · by_cases h2 : m y = 1 ∧ p = Function.update m y 0
    · exfalso
      apply h
      rw [count_univ_update m' x h1.1, count_univ_update m y h2.1, ← h1.2, ← h2.2]
    · simp [h2]
  · simp [h1]

/-- Quadratic Hamiltonians conserve the particle number. -/
lemma dGamma_eq_zero_of_count (A : Matrix (Fin N) (Fin N) ℝ) (m' m : Cfg N)
    (h : count univ m' ≠ count univ m) : dGamma A m' m = 0 := by
  simp only [dGamma, Matrix.sum_apply, Matrix.smul_apply,
    cre_mul_ann_eq_zero_of_count _ _ m' m h, smul_zero, Finset.sum_const_zero]

/-- The hopping terms conserve the particle number. -/
lemma hop_eq_zero_of_count (j : Fin N) (m' m : Cfg N) (h : count univ m' ≠ count univ m) :
    hop j m' m = 0 := by
  unfold hop
  split_ifs with hj
  · rw [Matrix.add_apply, cre_mul_ann_eq_zero_of_count _ _ m' m h,
      cre_mul_ann_eq_zero_of_count _ _ m' m h, add_zero]
  · rfl

/-- The LGT Hamiltonian conserves the matter particle number. -/
lemma H_eq_zero_of_count (J h : ℝ) (p' p : Cfg N × Cfg N)
    (hc : count univ p'.1 ≠ count univ p.1) : H J h p' p = 0 := by
  obtain ⟨n', x'⟩ := p'
  obtain ⟨n, x⟩ := p
  have hne : n' ≠ n := fun e => hc (by rw [e])
  simp only at hc
  simp only [H, Matrix.sub_apply, Matrix.smul_apply, Matrix.sum_apply, kronecker_apply,
    smul_eq_mul, hop_eq_zero_of_count _ n' n hc, zero_mul, Finset.sum_const_zero, mul_zero,
    Matrix.one_apply_ne hne, sub_zero]

/-- The evolved Fock state lives on configurations with the particle number of `n`. -/
lemma evo_fock_eq_zero (A : Matrix (Fin N) (Fin N) ℝ) (t : ℝ) (n m : Cfg N)
    (hm : count univ m ≠ count univ n) : (evo (dGamma A) t *ᵥ Pi.single n 1) m = 0 := by
  rw [CutLocalAux.mulVec_single_one_apply]
  exact evo_eq_zero_of_conserve _ (count univ) (fun a b hab => dGamma_eq_zero_of_count A a b hab)
    t m n hm

/-- The LGT state lives on matter configurations with the particle number of `n`. -/
lemma lgtState_eq_zero (J h : ℝ) (n : Cfg N) (t : ℝ) (p : Cfg N × Cfg N)
    (hp : count univ p.1 ≠ count univ n) : lgtState J h (Pi.single n 1) t p = 0 := by
  unfold lgtState
  simp only [Matrix.mulVec, dotProduct]
  refine Finset.sum_eq_zero fun p0 _ => ?_
  by_cases h0 : p0.1 = n
  · rw [evo_eq_zero_of_conserve (H J h) (fun q => count univ q.1)
      (fun a b hab => H_eq_zero_of_count J h a b hab) t p p0 (by rw [h0]; exact hp), zero_mul]
  · simp [initZ, h0]

end Conservation

/-! ### Particle numbers of split configurations -/

section Split

/-- Number of particles of a configuration, as a natural number. -/
def cnt {K : ℕ} (α : Cfg K) : ℕ := (univ.filter (fun i => α i ≠ 0)).card

/-- `cnt` is the particle number `count univ`. -/
lemma cnt_cast {K : ℕ} (α : Cfg K) : (cnt α : ℝ) = count univ α := by
  unfold cnt count occ
  rw [Finset.card_filter]
  push_cast
  refine Finset.sum_congr rfl fun x _ => ?_
  split_ifs <;> simp_all

/-- At most `K` particles fit in `K` modes. -/
lemma cnt_le {K : ℕ} (α : Cfg K) : cnt α ≤ K :=
  (Finset.card_filter_le _ _).trans (by simp)

/-- The particle number of a recombined configuration is the sum of the two parts. -/
lemma count_univ_splitCfg_symm (ℓ : ℕ) (hℓ : ℓ ≤ N) (α : Cfg ℓ) (β : Cfg (N - ℓ)) :
    count univ ((splitCfg ℓ hℓ).symm (α, β)) = count univ α + count univ β := by
  unfold count
  rw [AreaLawAux.sum_fin_split ℓ hℓ]
  congr 1
  · refine Finset.sum_congr rfl fun a _ => ?_
    rw [← AreaLawAux.splitCfg_fst ℓ hℓ ((splitCfg ℓ hℓ).symm (α, β)) a, Equiv.apply_symm_apply]
  · refine Finset.sum_congr rfl fun b _ => ?_
    rw [← AreaLawAux.splitCfg_snd ℓ hℓ ((splitCfg ℓ hℓ).symm (α, β)) b, Equiv.apply_symm_apply]

/-- The particle number in the left segment of a recombined configuration. -/
lemma count_segment_splitCfg_symm (ℓ : ℕ) (hℓ : ℓ ≤ N) (α : Cfg ℓ) (β : Cfg (N - ℓ)) :
    count (segment ℓ) ((splitCfg ℓ hℓ).symm (α, β)) = count univ α := by
  unfold count
  rw [AreaLawAux.sum_segment ℓ hℓ]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [← AreaLawAux.splitCfg_fst ℓ hℓ ((splitCfg ℓ hℓ).symm (α, β)) a, Equiv.apply_symm_apply]

/-- The matter part of a recombined LGT configuration. -/
lemma splitLGT_symm_fst (ℓ : ℕ) (hℓ : ℓ ≤ N) (a : Cfg ℓ × Cfg ℓ)
    (b : Cfg (N - ℓ) × Cfg (N - ℓ)) :
    ((splitLGT ℓ hℓ).symm (a, b)).1 = (splitCfg ℓ hℓ).symm (a.1, b.1) := rfl

/-- `Re (z conj z) = |z|²`. -/
lemma re_mul_star_self (z : ℂ) : (z * star z).re = Complex.normSq z := by
  rw [Complex.star_def, Complex.mul_conj, Complex.ofReal_re]

/-- Expectation of a complexified real diagonal matrix: `Re ⟨v, diag(d) v⟩ = ∑ d_i |v_i|²`. -/
lemma expect_diag_re {κ : Type*} [Fintype κ] [DecidableEq κ] (d : κ → ℝ) (v : κ → ℂ) :
    (expect (cplx (diagonal d)) v).re = ∑ i, d i * Complex.normSq (v i) := by
  rw [sec_cplx_diagonal]
  simp only [expect, dotProduct, mulVec_diagonal, Pi.star_apply, Complex.re_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [show star (v i) * ((d i : ℂ) * v i) = ((d i * Complex.normSq (v i) : ℝ) : ℂ) by
    rw [Complex.ofReal_mul, ← star_mul_self_eq]; ring, Complex.ofReal_re]

end Split

end TaxonomyAux

open TaxonomyAux

set_option linter.unusedVariables false in
/-- **Block entropy bound**: if the density matrix `ρ` is block diagonal with respect to a charge
`c : ι → ℕ` (with values `≤ M`), its von Neumann entropy is at least the Shannon entropy of the
distribution of the charge, `p_k = ∑_{a : c a = k} ρ_{aa}`. -/
theorem shannon_le_vonNeumann_of_block {ι : Type*} [Fintype ι] [DecidableEq ι]
    (ρ : Matrix ι ι ℂ) (hρ : ρ.PosSemidef) (htr : ρ.trace = 1) (c : ι → ℕ) (M : ℕ)
    (hcM : ∀ a, c a ≤ M) (hblock : ∀ a b, c a ≠ c b → ρ a b = 0) :
    shannon M (fun k => (∑ a ∈ univ.filter (fun a => c a = k), ρ a a).re) ≤
      vonNeumann ρ hρ.isHermitian :=
  -- the normalization `htr` is not needed
  block_bound ρ hρ c M hcM hblock

set_option linter.unusedVariables false in
/-- **Slater determinants** (every charge sector): the entanglement entropy of the left segment
`{0, …, ℓ-1}` is at least the number entropy of the segment. -/
theorem numDist_shannon_le_vonNeumann_slater (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm)
    (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    shannon N (numDist (segment ℓ) (evo (dGamma A) t *ᵥ Pi.single n 1)) ≤
      vonNeumann (rdm ℓ hℓ (evo (dGamma A) t *ᵥ Pi.single n 1)) (rdm_isHermitian ℓ hℓ _) := by
  -- the symmetry `hA` is not needed: `dΓ(A)` conserves the particle number for every `A`
  set v := evo (dGamma A) t *ᵥ Pi.single n 1 with hv
  have hsupp : ∀ m, v m ≠ 0 → count univ m = count univ n := fun m hm => by
    by_contra hc
    exact hm (evo_fock_eq_zero A t n m hc)
  -- `ρ_A = X X†` is positive semidefinite
  set X : Matrix (Cfg ℓ) (Cfg (N - ℓ)) ℂ := Matrix.of fun a b => v ((splitCfg ℓ hℓ).symm (a, b))
    with hX
  have hρX : rdm ℓ hℓ v = X * Xᴴ := by
    ext a a'
    simp [rdm, X, Matrix.mul_apply, Matrix.conjTranspose_apply]
  have hpsd : (rdm ℓ hℓ v).PosSemidef := hρX ▸ Matrix.posSemidef_self_mul_conjTranspose X
  -- `ρ_A` is block diagonal in the particle number of the segment
  have hblock : ∀ α α', cnt α ≠ cnt α' → rdm ℓ hℓ v α α' = 0 := by
    intro α α' hne
    refine Finset.sum_eq_zero fun β _ => ?_
    by_cases h1 : v ((splitCfg ℓ hℓ).symm (α, β)) = 0
    · rw [h1, zero_mul]
    by_cases h2 : v ((splitCfg ℓ hℓ).symm (α', β)) = 0
    · rw [h2, star_zero, mul_zero]
    exfalso
    apply hne
    have e1 := hsupp _ h1
    have e2 := hsupp _ h2
    rw [count_univ_splitCfg_symm] at e1 e2
    have : (cnt α : ℝ) = cnt α' := by rw [cnt_cast, cnt_cast]; linarith
    exact_mod_cast this
  have hmain := block_bound (rdm ℓ hℓ v) hpsd cnt N (fun α => (cnt_le α).trans hℓ) hblock
  -- the block weights are the number distribution of the segment
  have hdist : numDist (segment ℓ) v =
      fun k => (∑ α ∈ univ.filter (fun α => cnt α = k), rdm ℓ hℓ v α α).re := by
    funext k
    rw [numDist, Finset.sum_filter, ← (splitCfg ℓ hℓ).symm.sum_comp, Fintype.sum_prod_type,
      Complex.re_sum, Finset.sum_filter]
    refine Finset.sum_congr rfl fun α _ => ?_
    simp only [count_segment_splitCfg_symm, ← cnt_cast, Nat.cast_inj]
    split_ifs with hk
    · simp only [rdm, Complex.re_sum, re_mul_star_self]
    · simp
  rw [hdist]
  exact hmain

/-- **The LGT** (paper: Proposition `prop:entropies`(c) for `S_full`): for the state evolved from
`|↑…↑⟩_z ⊗ |n⟩`, the full matter + gauge entanglement entropy of the region
`A` = (sites `< ℓ`, links `< ℓ`) is at least the number entropy of the matter in `A`. -/
theorem lgt_numberEntropy_le_full (J h : ℝ) (n : Cfg N) (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ≤ N) :
    shannon N (lgtNumDist J h (segment ℓ) (Pi.single n 1) t) ≤
      vonNeumann (rdmLGT ℓ hℓ (lgtState J h (Pi.single n 1) t)) (rdmLGT_isHermitian ℓ hℓ _) := by
  set Ψ := lgtState J h (Pi.single n 1) t with hΨ
  have hsupp : ∀ p, Ψ p ≠ 0 → count univ p.1 = count univ n := fun p hp => by
    by_contra hc
    exact hp (lgtState_eq_zero J h n t p hc)
  -- `ρ_LGT = X X†` is positive semidefinite
  set X : Matrix (Cfg ℓ × Cfg ℓ) (Cfg (N - ℓ) × Cfg (N - ℓ)) ℂ :=
    Matrix.of fun a b => Ψ ((splitLGT ℓ hℓ).symm (a, b)) with hX
  have hρX : rdmLGT ℓ hℓ Ψ = X * Xᴴ := by
    ext a a'
    simp [rdmLGT, X, Matrix.mul_apply, Matrix.conjTranspose_apply]
  have hpsd : (rdmLGT ℓ hℓ Ψ).PosSemidef := hρX ▸ Matrix.posSemidef_self_mul_conjTranspose X
  -- `ρ_LGT` is block diagonal in the matter particle number of `A`
  have hblock : ∀ a a' : Cfg ℓ × Cfg ℓ, cnt a.1 ≠ cnt a'.1 → rdmLGT ℓ hℓ Ψ a a' = 0 := by
    intro a a' hne
    refine Finset.sum_eq_zero fun b _ => ?_
    by_cases h1 : Ψ ((splitLGT ℓ hℓ).symm (a, b)) = 0
    · rw [h1, zero_mul]
    by_cases h2 : Ψ ((splitLGT ℓ hℓ).symm (a', b)) = 0
    · rw [h2, star_zero, mul_zero]
    exfalso
    apply hne
    have e1 := hsupp _ h1
    have e2 := hsupp _ h2
    rw [splitLGT_symm_fst, count_univ_splitCfg_symm] at e1 e2
    have : (cnt a.1 : ℝ) = cnt a'.1 := by rw [cnt_cast, cnt_cast]; linarith
    exact_mod_cast this
  have hmain := block_bound (rdmLGT ℓ hℓ Ψ) hpsd (fun a => cnt a.1) N
    (fun a => (cnt_le a.1).trans hℓ) hblock
  -- the block weights are the distribution of the matter particle number of `A`
  have hdist : lgtNumDist J h (segment ℓ) (Pi.single n 1) t =
      fun k => (∑ a ∈ univ.filter (fun a : Cfg ℓ × Cfg ℓ => cnt a.1 = k),
        rdmLGT ℓ hℓ Ψ a a).re := by
    funext k
    simp only [lgtNumDist, matter, numProj]
    rw [← hΨ, sec_diagonal_kronecker_one, expect_diag_re, ← (splitLGT ℓ hℓ).symm.sum_comp,
      Fintype.sum_prod_type, Complex.re_sum, Finset.sum_filter]
    refine Finset.sum_congr rfl fun a _ => ?_
    simp only [splitLGT_symm_fst, count_segment_splitCfg_symm, ← cnt_cast, Nat.cast_inj]
    split_ifs with hk
    · simp only [one_mul, rdmLGT, Complex.re_sum, re_mul_star_self]
    · simp
  rw [hdist]
  exact hmain

end DFL
