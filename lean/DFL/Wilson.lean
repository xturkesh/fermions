import DFL.FreeDyn
import DFL.Entropy
import DFL.AreaLaw
import DFL.Duality

/-!
# Wilson links: cut-averaged number entropy from energy conservation
(paper: Proposition `prop:wilson`)

The lattice Schwinger model with Wilson links (paper: Eq. `eq:HU1`),
`H = -w ∑_j (f†_j U_j f_{j+1} + h.c.) + m ∑_j (-1)^j n_j + (g²/2) ∑_j (E_j + ϑ)²`,
in a fixed background-charge sector `Q`. The Gauss law `G_j = E_j - E_{j-1} - ρ_j = Q_j` (with
`E_0 = e₀` fixed) determines the electric fields from the occupations,
`E_j = e₀ + ∑_{k ≤ j} (ρ_k + Q_k)`, and the hopping `f†_j U_j f_{j+1}` raises exactly the field it
crosses (`wField_hop`). The sector is therefore the space of functions of the occupations, and `H`
acts on it as the `2^L × 2^L` real symmetric matrix `schwinger` (a Fock basis vector `|n⟩` stands for
`|n⟩ ⊗ |E(n)⟩`).

Indexing: Lean sites are `0, …, L-1`; Lean site `j` is paper site `j + 1`. So the paper's
`(-1)^j` is `(-1)^{j+1}` here, the staggered charge is `ρ_j = n_j - [j even]`, the paper's
interval `[1, ℓ]` is `segment ℓ`, and the paper's field `E_ℓ` (link `b_ℓ`) is `wField … ⟨ℓ - 1, _⟩`.
-/

open Matrix Finset Real

noncomputable section

namespace DFL

variable {L : ℕ}

/-- Staggered charge `ρ_j = n_j - [j even]` (0-indexed; paper: `n_j - ½(1 - (-1)^j)`, 1-indexed). -/
def stagRho (n : Cfg L) (j : Fin L) : ℤ :=
  (if n j = 1 then 1 else 0) - (if j.val % 2 = 0 then 1 else 0)

/-- The electric field on the link `b_j` in the sector `Q`: the Gauss law integrated from the left
edge, `E_j = e₀ + ∑_{k ≤ j} (ρ_k + Q_k)`. -/
def wField (e₀ : ℤ) (Q : Fin L → ℤ) (n : Cfg L) (j : Fin L) : ℤ :=
  e₀ + ∑ k ∈ univ.filter (fun k => k ≤ j), (stagRho n k + Q k)

/-- The lattice Schwinger Hamiltonian in the background-charge sector `Q`. -/
def schwinger (w m g θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) : Matrix (Cfg L) (Cfg L) ℝ :=
  (-w) • ∑ j, hop j +
    diagonal (fun n => m * ∑ j : Fin L, (-1 : ℝ) ^ (j.val + 1) * occ (n j) +
      g ^ 2 / 2 * ∑ j, ((wField e₀ Q n j : ℝ) + θ) ^ 2)

/-- The state at time `t`. -/
def wState (w m g θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) (ψ₀ : Cfg L → ℂ) (t : ℝ) : Cfg L → ℂ :=
  evo (schwinger w m g θ e₀ Q) t *ᵥ ψ₀

/-- The energy density `ε₀ = ⟨ψ₀, H ψ₀⟩ / L`. -/
def energyDensity (w m g θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) (ψ₀ : Cfg L → ℂ) : ℝ :=
  (expect (cplx (schwinger w m g θ e₀ Q)) ψ₀).re / L

/-- `⟨ψ, (E_j + ϑ)² ψ⟩`. -/
def fieldSq (θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) (ψ : Cfg L → ℂ) (j : Fin L) : ℝ :=
  ∑ n, Complex.normSq (ψ n) * ((wField e₀ Q n j : ℝ) + θ) ^ 2

/-- The constant in the summed Gauss law: `E_ℓ - N_{[1,ℓ]} = e₀ - #{odd paper sites in [1,ℓ]} + ∑_{k ≤ ℓ} Q_k`. -/
def wConst (e₀ : ℤ) (Q : Fin L → ℤ) (ℓ : ℕ) : ℝ :=
  (e₀ : ℝ) - ((univ.filter (fun k : Fin L => k.val < ℓ ∧ k.val % 2 = 0)).card : ℝ) +
    ∑ k ∈ univ.filter (fun k : Fin L => k.val < ℓ), (Q k : ℝ)

namespace WilsonAux

/-- Swapping two sites that both lie in `s` does not change the number of occupied sites of `s`. -/
theorem sum_occ_comp_swap (n : Cfg L) (a b : Fin L) (s : Finset (Fin L)) (ha : a ∈ s)
    (hb : b ∈ s) :
    ∑ i ∈ s, (if (n ∘ Equiv.swap a b) i = 1 then (1 : ℤ) else 0) =
      ∑ i ∈ s, (if n i = 1 then (1 : ℤ) else 0) := by
  refine Equiv.Perm.sum_comp (Equiv.swap a b) s (fun i => if n i = 1 then (1 : ℤ) else 0) ?_
  intro i hi
  by_contra hns
  have h1 : i ≠ a := by rintro rfl; exact hns ha
  have h2 : i ≠ b := by rintro rfl; exact hns hb
  exact hi (Equiv.swap_apply_of_ne_of_ne h1 h2)

/-- Swapping two sites outside `s` does not change the occupations in `s`. -/
theorem sum_occ_comp_swap_of_not_mem (n : Cfg L) (a b : Fin L) (s : Finset (Fin L))
    (ha : a ∉ s) (hb : b ∉ s) :
    ∑ i ∈ s, (if (n ∘ Equiv.swap a b) i = 1 then (1 : ℤ) else 0) =
      ∑ i ∈ s, (if n i = 1 then (1 : ℤ) else 0) := by
  refine Finset.sum_congr rfl (fun i hi => ?_)
  have h1 : i ≠ a := by rintro rfl; exact ha hi
  have h2 : i ≠ b := by rintro rfl; exact hb hi
  simp only [Function.comp_apply, Equiv.swap_apply_of_ne_of_ne h1 h2]

end WilsonAux

/-- The hopping across the bond `(j, j+1)` changes no field other than `E_j` (so `schwinger` is the
restriction of `H` to the sector). -/
theorem wField_hop (e₀ : ℤ) (Q : Fin L → ℤ) (j k : Fin L) (n' n : Cfg L) (hjk : k ≠ j)
    (hhop : hop j n' n ≠ 0) : wField e₀ Q n' k = wField e₀ Q n k := by
  rw [hop_apply] at hhop
  split_ifs at hhop with hj h
  · obtain ⟨-, rfl⟩ := h
    have key : ∑ i ∈ univ.filter (fun i => i ≤ k),
        (if (n ∘ Equiv.swap j ⟨j.val + 1, hj⟩) i = 1 then (1 : ℤ) else 0) =
          ∑ i ∈ univ.filter (fun i => i ≤ k), (if n i = 1 then (1 : ℤ) else 0) := by
      rcases lt_or_gt_of_ne hjk with hlt | hgt
      · refine WilsonAux.sum_occ_comp_swap_of_not_mem n _ _ _ ?_ ?_
        · simp only [mem_filter, mem_univ, true_and, not_le]
          exact hlt
        · simp only [mem_filter, mem_univ, true_and, not_le, Fin.lt_def]
          rw [Fin.lt_def] at hlt
          omega
      · refine WilsonAux.sum_occ_comp_swap n _ _ _ ?_ ?_
        · simp only [mem_filter, mem_univ, true_and]
          exact hgt.le
        · simp only [mem_filter, mem_univ, true_and, Fin.le_def]
          rw [Fin.lt_def] at hgt
          omega
    unfold wField stagRho
    simp only [Finset.sum_add_distrib, Finset.sum_sub_distrib, key]
  · exact absurd rfl hhop
  · exact absurd rfl hhop

/-- **Summed Gauss law** (paper: Proposition `prop:wilson`(a), via Eq. `eq:summedgauss`): in the
sector, `E_ℓ` is an injective affine function of `N_{[1,ℓ]}`, so their laws have the same Shannon
entropy. -/
theorem wField_eq_count (e₀ : ℤ) (Q : Fin L → ℤ) (n : Cfg L) (ℓ : ℕ) (hℓ : 0 < ℓ) (hℓL : ℓ ≤ L) :
    (wField e₀ Q n ⟨ℓ - 1, by omega⟩ : ℝ) = count (segment ℓ) n + wConst e₀ Q ℓ := by
  have hset : univ.filter (fun k : Fin L => k ≤ ⟨ℓ - 1, by omega⟩) =
      univ.filter (fun k : Fin L => k.val < ℓ) := by
    ext k
    simp only [mem_filter, mem_univ, true_and, Fin.le_def]
    omega
  have hocc : ∀ x : Fin L, ((if n x = 1 then 1 else 0 : ℤ) : ℝ) = occ (n x) := by
    intro x
    rcases FreeDynAux.zmod2_cases (n x) with h | h <;> simp [h, occ]
  have hodd : ∑ x ∈ univ.filter (fun k : Fin L => k.val < ℓ),
      ((if x.val % 2 = 0 then 1 else 0 : ℤ) : ℝ) =
        ((univ.filter (fun k : Fin L => k.val < ℓ ∧ k.val % 2 = 0)).card : ℝ) := by
    push_cast
    rw [Finset.sum_boole, Finset.filter_filter]
  unfold wField count segment wConst stagRho
  rw [hset]
  push_cast [Finset.sum_add_distrib, Finset.sum_sub_distrib]
  rw [← hodd]
  simp only [← hocc]
  push_cast
  ring

namespace WilsonAux

section Expect

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

omit [Fintype ι] [DecidableEq ι] in
/-- `cplx` is additive. -/
theorem cplx_add (A B : Matrix ι ι ℝ) : cplx (A + B) = cplx A + cplx B := by
  ext i j; simp [cplx]

omit [DecidableEq ι] in
/-- `expect` is additive in the observable. -/
theorem expect_add (M N : Matrix ι ι ℂ) (v : ι → ℂ) :
    expect (M + N) v = expect M v + expect N v := by
  unfold expect
  rw [Matrix.add_mulVec, dotProduct_add]

omit [DecidableEq ι] in
/-- `expect` is homogeneous in the observable. -/
theorem expect_smul_left (c : ℂ) (M : Matrix ι ι ℂ) (v : ι → ℂ) :
    expect (c • M) v = c * expect M v := by
  unfold expect
  rw [Matrix.smul_mulVec, dotProduct_smul, smul_eq_mul]

/-- Expectation of a complexified real diagonal matrix: `⟨v, diag(d) v⟩ = ∑ d_i |v_i|²`. -/
theorem expect_cplx_diagonal (d : ι → ℝ) (v : ι → ℂ) :
    expect (cplx (diagonal d)) v = ((∑ i, d i * Complex.normSq (v i) : ℝ) : ℂ) := by
  rw [expect, FreeDynAux.cplx_diagonal]
  simp only [dotProduct, mulVec_diagonal, Pi.star_apply, Complex.ofReal_sum, Complex.ofReal_mul]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [Complex.normSq_eq_conj_mul_self]
  simp only [Complex.star_def]
  ring

omit [DecidableEq ι] in
/-- Real part of the expectation of a complexified real matrix. -/
theorem re_expect_cplx (A : Matrix ι ι ℝ) (v : ι → ℂ) :
    (expect (cplx A) v).re = ∑ a, ∑ b, A a b * (star (v a) * v b).re := by
  unfold expect
  simp only [dotProduct, mulVec, Pi.star_apply, Finset.mul_sum, Complex.re_sum,
    FreeDynAux.cplx_apply]
  refine Finset.sum_congr rfl (fun a _ => Finset.sum_congr rfl (fun b _ => ?_))
  rw [mul_left_comm, Complex.re_ofReal_mul]

/-- `Re (conj x * y) ≤ (|x|² + |y|²) / 2`. -/
theorem re_star_mul_le (x y : ℂ) :
    (star x * y).re ≤ (Complex.normSq x + Complex.normSq y) / 2 := by
  simp only [Complex.star_def, Complex.mul_re, Complex.conj_re, Complex.conj_im,
    Complex.normSq_apply]
  nlinarith [sq_nonneg (x.re - y.re), sq_nonneg (x.im - y.im)]

omit [DecidableEq ι] in
/-- **Schur test**: a real matrix with nonnegative entries whose row and column sums are at most
`c` has `Re ⟨v, A v⟩ ≤ c ‖v‖²`. -/
theorem re_expect_le (A : Matrix ι ι ℝ) (c : ℝ) (h0 : ∀ a b, 0 ≤ A a b)
    (hrow : ∀ a, ∑ b, A a b ≤ c) (hcol : ∀ b, ∑ a, A a b ≤ c) (v : ι → ℂ) :
    (expect (cplx A) v).re ≤ c * ∑ a, Complex.normSq (v a) := by
  rw [re_expect_cplx]
  have h1 : ∑ a, ∑ b, A a b * Complex.normSq (v a) ≤ c * ∑ a, Complex.normSq (v a) := by
    rw [Finset.mul_sum]
    refine Finset.sum_le_sum (fun a _ => ?_)
    rw [← Finset.sum_mul]
    exact mul_le_mul_of_nonneg_right (hrow a) (Complex.normSq_nonneg _)
  have h2 : ∑ a, ∑ b, A a b * Complex.normSq (v b) ≤ c * ∑ b, Complex.normSq (v b) := by
    rw [Finset.sum_comm, Finset.mul_sum]
    refine Finset.sum_le_sum (fun b _ => ?_)
    rw [← Finset.sum_mul]
    exact mul_le_mul_of_nonneg_right (hcol b) (Complex.normSq_nonneg _)
  have h3 : ∑ a, ∑ b, A a b * (star (v a) * v b).re ≤
      ∑ a, ∑ b, (A a b * Complex.normSq (v a) + A a b * Complex.normSq (v b)) / 2 := by
    refine Finset.sum_le_sum (fun a _ => Finset.sum_le_sum (fun b _ => ?_))
    have h := re_star_mul_le (v a) (v b)
    have ha := h0 a b
    nlinarith
  have h4 : ∑ a, ∑ b, (A a b * Complex.normSq (v a) + A a b * Complex.normSq (v b)) / 2 =
      (∑ a, ∑ b, A a b * Complex.normSq (v a) + ∑ a, ∑ b, A a b * Complex.normSq (v b)) / 2 := by
    simp only [add_div, Finset.sum_div, Finset.sum_add_distrib]
  linarith

omit [DecidableEq ι] in
/-- `∑ |v_i|² = star v ⬝ᵥ v`. -/
theorem sum_normSq_eq_dotProduct (w : ι → ℂ) :
    ((∑ i, Complex.normSq (w i) : ℝ) : ℂ) = star w ⬝ᵥ w := by
  simp only [Complex.ofReal_sum, dotProduct, Pi.star_apply]
  refine Finset.sum_congr rfl (fun i _ => ?_)
  rw [Complex.normSq_eq_conj_mul_self]
  rfl

/-- An isometry preserves `∑ |v_i|²`. -/
theorem sum_normSq_mulVec (U : Matrix ι ι ℂ) (hU : Uᴴ * U = 1) (v : ι → ℂ) :
    ∑ i, Complex.normSq ((U *ᵥ v) i) = ∑ i, Complex.normSq (v i) := by
  have h1 : star (U *ᵥ v) ⬝ᵥ (U *ᵥ v) = star v ⬝ᵥ v := by
    rw [star_mulVec, dotProduct_mulVec, vecMul_vecMul, hU, vecMul_one]
  exact_mod_cast (sum_normSq_eq_dotProduct (U *ᵥ v)).trans
    (h1.trans (sum_normSq_eq_dotProduct v).symm)

/-- **Energy conservation**: the expectation of a real symmetric `A` is constant along `e^{-itA}`. -/
theorem expect_evo (A : Matrix ι ι ℝ) (hA : A.IsSymm) (t : ℝ) (v : ι → ℂ) :
    expect (cplx A) (evo A t *ᵥ v) = expect (cplx A) v := by
  rw [FreeDynAux.expect_mulVec]
  have hc : Commute (cplx A) (evo A t) := by
    unfold evo
    exact ((Commute.refl (cplx A)).smul_right _).exp_right
  rw [Matrix.mul_assoc, hc.eq, ← Matrix.mul_assoc, evo_unitary A hA t, Matrix.one_mul]

end Expect

/-- The hopping across a bond is a symmetric matrix. -/
theorem hop_isSymm (j : Fin L) : (hop j).IsSymm := by
  unfold hop
  split_ifs with hj
  · simp only [Matrix.IsSymm, Matrix.transpose_add, Matrix.transpose_mul, cre,
      Matrix.transpose_transpose]
    rw [add_comm]
  · exact Matrix.isSymm_zero

/-- The Schwinger Hamiltonian is a real symmetric matrix. -/
theorem schwinger_isSymm (w m g θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) :
    (schwinger w m g θ e₀ Q).IsSymm := by
  unfold schwinger
  refine Matrix.IsSymm.add (Matrix.IsSymm.smul ?_ _) (Matrix.isSymm_diagonal _)
  simp only [Matrix.IsSymm, Matrix.transpose_sum]
  exact Finset.sum_congr rfl (fun j _ => hop_isSymm j)

/-- The evolved state stays normalized. -/
theorem sum_normSq_wState (w m g θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) (ψ₀ : Cfg L → ℂ)
    (hψ₀ : ∑ n, Complex.normSq (ψ₀ n) = 1) (t : ℝ) :
    ∑ n, Complex.normSq (wState w m g θ e₀ Q ψ₀ t n) = 1 := by
  rw [wState, sum_normSq_mulVec _ (evo_unitary _ (schwinger_isSymm w m g θ e₀ Q) t), hψ₀]

/-- The entries of `hop j` are nonnegative. -/
theorem hop_nonneg (j : Fin L) (n' n : Cfg L) : 0 ≤ hop j n' n := by
  rw [hop_apply]
  split_ifs <;> norm_num

/-- Each row of `hop j` sums to at most `1`. -/
theorem hop_row_le (j : Fin L) (n' : Cfg L) : ∑ n, hop j n' n ≤ 1 := by
  by_cases hj : j.val + 1 < L
  · calc ∑ n, hop j n' n
        ≤ ∑ n, (if n = n' ∘ Equiv.swap j ⟨j.val + 1, hj⟩ then (1 : ℝ) else 0) := by
          refine Finset.sum_le_sum (fun n _ => ?_)
          rw [hop_apply, dite_eq_left hj]
          split_ifs with h1 h2 h2
          · exact le_rfl
          · exfalso
            apply h2
            rw [h1.2]
            funext x
            simp [Equiv.swap_apply_self]
          · norm_num
          · exact le_rfl
      _ = 1 := by simp
  · simp [hop_apply, hj]

/-- Each column of `hop j` sums to at most `1`. -/
theorem hop_col_le (j : Fin L) (n : Cfg L) : ∑ n', hop j n' n ≤ 1 := by
  by_cases hj : j.val + 1 < L
  · calc ∑ n', hop j n' n
        ≤ ∑ n', (if n' = n ∘ Equiv.swap j ⟨j.val + 1, hj⟩ then (1 : ℝ) else 0) := by
          refine Finset.sum_le_sum (fun n' _ => ?_)
          rw [hop_apply, dite_eq_left hj]
          split_ifs with h1 h2 h2
          · exact le_rfl
          · exact absurd h1.2 h2
          · norm_num
          · exact le_rfl
      _ = 1 := by simp
  · simp [hop_apply, hj]

/-- The hopping part satisfies `Re ⟨v, (∑_j hop j) v⟩ ≤ L ‖v‖²`. -/
theorem re_expect_hop_le (v : Cfg L → ℂ) :
    (expect (cplx (∑ j : Fin L, hop j)) v).re ≤ L * ∑ n, Complex.normSq (v n) := by
  refine re_expect_le _ _ ?_ ?_ ?_ v
  · intro a b
    rw [Matrix.sum_apply]
    exact Finset.sum_nonneg (fun j _ => hop_nonneg j a b)
  · intro a
    simp only [Matrix.sum_apply]
    rw [Finset.sum_comm]
    calc ∑ j : Fin L, ∑ b, hop j a b ≤ ∑ _j : Fin L, (1 : ℝ) :=
          Finset.sum_le_sum (fun j _ => hop_row_le j a)
      _ = L := by simp
  · intro b
    simp only [Matrix.sum_apply]
    rw [Finset.sum_comm]
    calc ∑ j : Fin L, ∑ a, hop j a b ≤ ∑ _j : Fin L, (1 : ℝ) :=
          Finset.sum_le_sum (fun j _ => hop_col_le j b)
      _ = L := by simp

/-- The staggered mass term is at least `-|m| L`. -/
theorem mass_ge (m : ℝ) (n : Cfg L) :
    -(|m| * L) ≤ m * ∑ j : Fin L, (-1 : ℝ) ^ (j.val + 1) * occ (n j) := by
  have hS : |∑ j : Fin L, (-1 : ℝ) ^ (j.val + 1) * occ (n j)| ≤ L := by
    calc |∑ j : Fin L, (-1 : ℝ) ^ (j.val + 1) * occ (n j)|
        ≤ ∑ j : Fin L, |(-1 : ℝ) ^ (j.val + 1) * occ (n j)| := Finset.abs_sum_le_sum_abs _ _
      _ ≤ ∑ _j : Fin L, (1 : ℝ) := by
          refine Finset.sum_le_sum (fun j _ => ?_)
          rw [abs_mul, abs_pow, abs_neg, abs_one, one_pow, one_mul]
          unfold occ
          split_ifs <;> norm_num
      _ = L := by simp
  have h1 : |m * ∑ j : Fin L, (-1 : ℝ) ^ (j.val + 1) * occ (n j)| ≤ |m| * L := by
    rw [abs_mul]
    exact mul_le_mul_of_nonneg_left hS (abs_nonneg m)
  linarith [neg_abs_le (m * ∑ j : Fin L, (-1 : ℝ) ^ (j.val + 1) * occ (n j))]

/-- `⟨(E_j + ϑ)²⟩ ≥ 0`. -/
theorem fieldSq_nonneg (θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) (v : Cfg L → ℂ) (j : Fin L) :
    0 ≤ fieldSq θ e₀ Q v j :=
  Finset.sum_nonneg (fun _ _ => mul_nonneg (Complex.normSq_nonneg _) (sq_nonneg _))

/-- The real part of the energy, split into hopping, mass and field parts. -/
theorem re_expect_schwinger (w m g θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) (v : Cfg L → ℂ) :
    (expect (cplx (schwinger w m g θ e₀ Q)) v).re =
      -w * (expect (cplx (∑ j : Fin L, hop j)) v).re +
        ∑ n, (m * ∑ j : Fin L, (-1 : ℝ) ^ (j.val + 1) * occ (n j)) * Complex.normSq (v n) +
          g ^ 2 / 2 * ∑ j, fieldSq θ e₀ Q v j := by
  unfold schwinger
  rw [cplx_add, FreeDynAux.cplx_smul, expect_add, expect_smul_left, expect_cplx_diagonal,
    Complex.add_re, Complex.re_ofReal_mul, Complex.ofReal_re]
  unfold fieldSq
  rw [add_assoc]
  congr 1
  simp only [add_mul, Finset.sum_add_distrib]
  congr 1
  calc ∑ n, (g ^ 2 / 2 * ∑ j, ((wField e₀ Q n j : ℝ) + θ) ^ 2) * Complex.normSq (v n)
      = ∑ n, ∑ j, g ^ 2 / 2 * (Complex.normSq (v n) * ((wField e₀ Q n j : ℝ) + θ) ^ 2) := by
        refine Finset.sum_congr rfl (fun n _ => ?_)
        rw [Finset.mul_sum, Finset.sum_mul]
        exact Finset.sum_congr rfl (fun j _ => by ring)
    _ = ∑ j, ∑ n, g ^ 2 / 2 * (Complex.normSq (v n) * ((wField e₀ Q n j : ℝ) + θ) ^ 2) :=
        Finset.sum_comm
    _ = g ^ 2 / 2 * ∑ j, ∑ n, Complex.normSq (v n) * ((wField e₀ Q n j : ℝ) + θ) ^ 2 := by
        rw [Finset.mul_sum]
        exact Finset.sum_congr rfl (fun j _ => by rw [Finset.mul_sum])

/-- **Field energy bound** for a normalized state: `(g²/2) ∑_j ⟨(E_j + ϑ)²⟩ ≤ ⟨H⟩ + L (w + |m|)`. -/
theorem field_le_energy (w m g θ : ℝ) (hw : 0 ≤ w) (e₀ : ℤ) (Q : Fin L → ℤ) (v : Cfg L → ℂ)
    (hv : ∑ n, Complex.normSq (v n) = 1) :
    g ^ 2 / 2 * ∑ j, fieldSq θ e₀ Q v j ≤
      (expect (cplx (schwinger w m g θ e₀ Q)) v).re + L * (w + |m|) := by
  rw [re_expect_schwinger]
  have hB := re_expect_hop_le v
  rw [hv, mul_one] at hB
  have hM : -(|m| * L) ≤
      ∑ n, (m * ∑ j : Fin L, (-1 : ℝ) ^ (j.val + 1) * occ (n j)) * Complex.normSq (v n) := by
    calc -(|m| * L) = ∑ n, (-(|m| * L)) * Complex.normSq (v n) := by
          rw [← Finset.mul_sum, hv, mul_one]
      _ ≤ _ := Finset.sum_le_sum (fun n _ =>
          mul_le_mul_of_nonneg_right (mass_ge m n) (Complex.normSq_nonneg _))
  nlinarith

end WilsonAux

/-- **Energy bound on the fields** (paper: Proposition `prop:wilson`(b)): `ε₀ + w + |m| ≥ 0` and
`(1/L) ∑_ℓ ⟨(E_ℓ + ϑ)²⟩_t ≤ (2/g²)(ε₀ + w + |m|)` for all `t`. -/
theorem wilson_field_bound (w m g θ : ℝ) (hw : 0 ≤ w) (hg : 0 < g) (e₀ : ℤ) (Q : Fin L → ℤ)
    (hL : 0 < L) (ψ₀ : Cfg L → ℂ) (hψ₀ : ∑ n, Complex.normSq (ψ₀ n) = 1) (t : ℝ) :
    0 ≤ energyDensity w m g θ e₀ Q ψ₀ + w + |m| ∧
      (1 / (L : ℝ)) * ∑ j, fieldSq θ e₀ Q (wState w m g θ e₀ Q ψ₀ t) j ≤
        2 / g ^ 2 * (energyDensity w m g θ e₀ Q ψ₀ + w + |m|) := by
  have hLpos : (0 : ℝ) < L := Nat.cast_pos.2 hL
  have hg2 : 0 < g ^ 2 := by positivity
  have hnorm := WilsonAux.sum_normSq_wState w m g θ e₀ Q ψ₀ hψ₀ t
  have hF := WilsonAux.field_le_energy w m g θ hw e₀ Q _ hnorm
  have hE : (expect (cplx (schwinger w m g θ e₀ Q)) (wState w m g θ e₀ Q ψ₀ t)).re =
      L * energyDensity w m g θ e₀ Q ψ₀ := by
    rw [wState, WilsonAux.expect_evo _ (WilsonAux.schwinger_isSymm w m g θ e₀ Q), energyDensity]
    field_simp
  rw [hE] at hF
  have hS0 : 0 ≤ ∑ j, fieldSq θ e₀ Q (wState w m g θ e₀ Q ψ₀ t) j :=
    Finset.sum_nonneg (fun j _ => WilsonAux.fieldSq_nonneg θ e₀ Q _ j)
  set S := ∑ j, fieldSq θ e₀ Q (wState w m g θ e₀ Q ψ₀ t) j with hSdef
  set X := energyDensity w m g θ e₀ Q ψ₀ + w + |m| with hXdef
  have hF' : g ^ 2 / 2 * S ≤ L * X := by rw [hXdef]; linarith
  have hX : 0 ≤ X := by
    by_contra h
    have : L * X < 0 := mul_neg_of_pos_of_neg hLpos (not_le.1 h)
    nlinarith
  refine ⟨hX, ?_⟩
  rw [one_div_mul_eq_div, div_le_iff₀ hLpos, div_mul_eq_mul_div, div_mul_eq_mul_div,
    le_div_iff₀ hg2]
  nlinarith

namespace WilsonAux

/-- The particle number of `R` in the configuration `n`, as a natural number. -/
def countN (R : Finset (Fin L)) (n : Cfg L) : ℕ := (R.filter (fun x => n x ≠ 0)).card

/-- `count` is the natural number `countN`. -/
theorem count_eq_countN (R : Finset (Fin L)) (n : Cfg L) : count R n = countN R n := by
  unfold count countN occ
  rw [Finset.card_filter]
  push_cast
  refine Finset.sum_congr rfl (fun x _ => ?_)
  split_ifs <;> simp_all

/-- The particle number is at most `L`. -/
theorem countN_le (R : Finset (Fin L)) (n : Cfg L) : countN R n ≤ L :=
  (Finset.card_filter_le _ _).trans ((Finset.card_le_univ R).trans (by simp))

/-- Summing over the level sets `{count R = k}`, `k ∈ {0, …, L}`, is summing over everything. -/
theorem sum_fiber_count (R : Finset (Fin L)) (φ : Cfg L → ℝ) :
    ∑ k ∈ range (L + 1), ∑ n ∈ univ.filter (fun n => count R n = k), φ n = ∑ n, φ n := by
  have h : ∀ k : ℕ, univ.filter (fun n : Cfg L => count R n = k) =
      univ.filter (fun n => countN R n = k) := by
    intro k; ext n; simp [count_eq_countN]
  simp_rw [h]
  exact Finset.sum_fiberwise_of_maps_to
    (fun n _ => mem_range.2 (Nat.lt_succ_of_le (countN_le R n))) φ

/-- The number distribution of a normalized vector is a probability distribution. -/
theorem numDist_isProb (R : Finset (Fin L)) (v : Cfg L → ℂ)
    (hv : ∑ n, Complex.normSq (v n) = 1) : IsProb L (numDist R v) where
  nonneg _ := Finset.sum_nonneg (fun _ _ => Complex.normSq_nonneg _)
  sum_one := by unfold numDist; rw [sum_fiber_count]; exact hv

/-- Second moment about `c` of the number distribution. -/
theorem sum_numDist_mul_sq (R : Finset (Fin L)) (v : Cfg L → ℂ) (c : ℝ) :
    ∑ k ∈ range (L + 1), numDist R v k * ((k : ℝ) - c) ^ 2 =
      ∑ n, Complex.normSq (v n) * (count R n - c) ^ 2 := by
  unfold numDist
  simp_rw [Finset.sum_mul]
  rw [← sum_fiber_count R]
  refine Finset.sum_congr rfl (fun k _ => Finset.sum_congr rfl (fun n hn => ?_))
  rw [(mem_filter.1 hn).2]

/-- **Per-cut bound**: `S_N([1,ℓ]) ≤ ½ log(2πe(⟨(E_ℓ + ϑ)²⟩ + 1/12))` (summed Gauss law + Massey). -/
theorem shannon_le_fieldSq (θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) (v : Cfg L → ℂ)
    (hv : ∑ n, Complex.normSq (v n) = 1) (ℓ : ℕ) (hℓ : 0 < ℓ) (hℓL : ℓ ≤ L) :
    shannon L (numDist (segment ℓ) v) ≤
      1 / 2 * Real.log (2 * π * exp 1 * (fieldSq θ e₀ Q v ⟨ℓ - 1, by omega⟩ + 1 / 12)) := by
  have hp := numDist_isProb (segment ℓ) v hv
  refine shannon_le_of_pvar_le hp
    ((pvar_le_sum_sq hp (-(wConst e₀ Q ℓ + θ))).trans (le_of_eq ?_))
  rw [sum_numDist_mul_sq]
  unfold fieldSq
  refine Finset.sum_congr rfl (fun n _ => ?_)
  rw [wField_eq_count e₀ Q n ℓ hℓ hℓL]
  ring

/-- The bound `x ↦ ½ log(2πe(x + 1/12))` is monotone on `[0, ∞)`. -/
theorem phi_mono {x y : ℝ} (hx : 0 ≤ x) (hxy : x ≤ y) :
    1 / 2 * Real.log (2 * π * exp 1 * (x + 1 / 12)) ≤
      1 / 2 * Real.log (2 * π * exp 1 * (y + 1 / 12)) := by
  have h1 : 0 < 2 * π * exp 1 * (x + 1 / 12) := by positivity
  gcongr

/-- The field profile along the cuts: `cutField θ e₀ Q v ℓ = ⟨(E_ℓ + ϑ)²⟩` for `ℓ ∈ [1, L]`. -/
def cutField (θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) (v : Cfg L → ℂ) (ℓ : ℕ) : ℝ :=
  if h : 0 < ℓ ∧ ℓ ≤ L then fieldSq θ e₀ Q v ⟨ℓ - 1, by omega⟩ else 0

/-- The field profile is nonnegative. -/
theorem cutField_nonneg (θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) (v : Cfg L → ℂ) (ℓ : ℕ) :
    0 ≤ cutField θ e₀ Q v ℓ := by
  unfold cutField
  split_ifs
  · exact fieldSq_nonneg θ e₀ Q v _
  · exact le_rfl

/-- Summing the field profile over the cuts `ℓ ∈ [1, L]` is summing over the links. -/
theorem sum_cutField (θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) (v : Cfg L → ℂ) :
    ∑ ℓ ∈ Icc 1 L, cutField θ e₀ Q v ℓ = ∑ j, fieldSq θ e₀ Q v j := by
  symm
  refine Finset.sum_bij (fun j _ => j.val + 1) ?_ ?_ ?_ ?_
  · intro j _
    simp only [mem_Icc]
    omega
  · intro a _ b _ h
    exact Fin.ext (by simpa using h)
  · intro ℓ hℓ
    simp only [mem_Icc] at hℓ
    exact ⟨⟨ℓ - 1, by omega⟩, mem_univ _, by simp only; omega⟩
  · intro j _
    unfold cutField
    rw [dite_eq_left ⟨by omega, by omega⟩]
    congr 1

/-- Per-cut bound in terms of the field profile. -/
theorem shannon_le_cutField (θ : ℝ) (e₀ : ℤ) (Q : Fin L → ℤ) (v : Cfg L → ℂ)
    (hv : ∑ n, Complex.normSq (v n) = 1) (ℓ : ℕ) (hℓ : ℓ ∈ Icc 1 L) :
    shannon L (numDist (segment ℓ) v) ≤
      1 / 2 * Real.log (2 * π * exp 1 * (cutField θ e₀ Q v ℓ + 1 / 12)) := by
  rw [mem_Icc] at hℓ
  unfold cutField
  rw [dite_eq_left ⟨by omega, hℓ.2⟩]
  exact shannon_le_fieldSq θ e₀ Q v hv ℓ (by omega) hℓ.2

/-- The summed field bound: `∑_j ⟨(E_j + ϑ)²⟩_t ≤ L (2/g²)(ε₀ + w + |m|)`. -/
theorem sum_fieldSq_le (w m g θ : ℝ) (hw : 0 ≤ w) (hg : 0 < g) (e₀ : ℤ) (Q : Fin L → ℤ)
    (hL : 0 < L) (ψ₀ : Cfg L → ℂ) (hψ₀ : ∑ n, Complex.normSq (ψ₀ n) = 1) (t : ℝ) :
    ∑ j, fieldSq θ e₀ Q (wState w m g θ e₀ Q ψ₀ t) j ≤
      L * (2 / g ^ 2 * (energyDensity w m g θ e₀ Q ψ₀ + w + |m|)) := by
  have hLpos : (0 : ℝ) < L := Nat.cast_pos.2 hL
  have h := (wilson_field_bound w m g θ hw hg e₀ Q hL ψ₀ hψ₀ t).2
  rw [one_div_mul_eq_div, div_le_iff₀ hLpos] at h
  linarith

end WilsonAux

/-- **Cut-averaged number entropy** (paper: Proposition `prop:wilson`(c)). -/
theorem wilson_numberEntropy_avg (w m g θ : ℝ) (hw : 0 ≤ w) (hg : 0 < g) (e₀ : ℤ)
    (Q : Fin L → ℤ) (hL : 0 < L) (ψ₀ : Cfg L → ℂ) (hψ₀ : ∑ n, Complex.normSq (ψ₀ n) = 1)
    (t : ℝ) :
    (1 / (L : ℝ)) * ∑ ℓ ∈ Icc 1 L, shannon L (numDist (segment ℓ) (wState w m g θ e₀ Q ψ₀ t)) ≤
      1 / 2 * Real.log (2 * π * exp 1 *
        (2 / g ^ 2 * (energyDensity w m g θ e₀ Q ψ₀ + w + |m|) + 1 / 12)) := by
  have hLpos : (0 : ℝ) < L := Nat.cast_pos.2 hL
  have hnorm := WilsonAux.sum_normSq_wState w m g θ e₀ Q ψ₀ hψ₀ t
  have hbound := (wilson_field_bound w m g θ hw hg e₀ Q hL ψ₀ hψ₀ t).2
  rw [← WilsonAux.sum_cutField] at hbound
  set ψ := wState w m g θ e₀ Q ψ₀ t with hψ
  set G := WilsonAux.cutField θ e₀ Q ψ with hG
  have hG0 : ∀ ℓ, 0 ≤ G ℓ := fun ℓ => WilsonAux.cutField_nonneg θ e₀ Q ψ ℓ
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
  calc (1 / (L : ℝ)) * ∑ ℓ ∈ Icc 1 L, shannon L (numDist (segment ℓ) ψ)
      ≤ (1 / (L : ℝ)) * ∑ ℓ ∈ Icc 1 L, 1 / 2 * Real.log (2 * π * exp 1 * (G ℓ + 1 / 12)) := by
        gcongr with ℓ hℓ
        exact WilsonAux.shannon_le_cutField θ e₀ Q ψ hnorm ℓ hℓ
    _ = 1 / 2 * ∑ ℓ ∈ Icc 1 L, (1 / (L : ℝ)) • Real.log (2 * π * exp 1 * (G ℓ + 1 / 12)) := by
        simp only [smul_eq_mul]
        rw [Finset.mul_sum, Finset.mul_sum]
        exact Finset.sum_congr rfl (fun ℓ _ => by ring)
    _ ≤ 1 / 2 * Real.log (∑ ℓ ∈ Icc 1 L, (1 / (L : ℝ)) • (2 * π * exp 1 * (G ℓ + 1 / 12))) := by
        gcongr
    _ = 1 / 2 * Real.log (2 * π * exp 1 * ((1 / (L : ℝ)) * ∑ ℓ ∈ Icc 1 L, G ℓ + 1 / 12)) := by
        rw [havg]
    _ ≤ _ := WilsonAux.phi_mono hA0 hbound

set_option linter.unusedVariables false in
/-- **Most cuts** (paper: Proposition `prop:wilson`(d)): for `δ ∈ (0, 1]`, the bound with `δ g²` in
place of `g²` fails for at most `δ L` of the cuts. -/
theorem wilson_numberEntropy_most (w m g θ : ℝ) (hw : 0 ≤ w) (hg : 0 < g) (e₀ : ℤ)
    (Q : Fin L → ℤ) (hL : 0 < L) (ψ₀ : Cfg L → ℂ) (hψ₀ : ∑ n, Complex.normSq (ψ₀ n) = 1)
    (t δ : ℝ) (hδ0 : 0 < δ) (hδ1 : δ ≤ 1) :
    (((Icc 1 L).filter (fun ℓ =>
        1 / 2 * Real.log (2 * π * exp 1 *
          (2 * (energyDensity w m g θ e₀ Q ψ₀ + w + |m|) / (δ * g ^ 2) + 1 / 12)) <
        shannon L (numDist (segment ℓ) (wState w m g θ e₀ Q ψ₀ t)))).card : ℝ) ≤ δ * L := by
  -- the Markov argument does not need `δ ≤ 1`
  have hLpos : (0 : ℝ) < L := Nat.cast_pos.2 hL
  have hg2 : 0 < g ^ 2 := by positivity
  have hnorm := WilsonAux.sum_normSq_wState w m g θ e₀ Q ψ₀ hψ₀ t
  have hX := (wilson_field_bound w m g θ hw hg e₀ Q hL ψ₀ hψ₀ t).1
  have hsum := WilsonAux.sum_fieldSq_le w m g θ hw hg e₀ Q hL ψ₀ hψ₀ t
  rw [← WilsonAux.sum_cutField] at hsum
  set ψ := wState w m g θ e₀ Q ψ₀ t with hψ
  set X := energyDensity w m g θ e₀ Q ψ₀ + w + |m| with hXdef
  set G := WilsonAux.cutField θ e₀ Q ψ with hG
  have hG0 : ∀ ℓ, 0 ≤ G ℓ := fun ℓ => WilsonAux.cutField_nonneg θ e₀ Q ψ ℓ
  set c := 2 * X / (δ * g ^ 2) with hc
  have hc0 : 0 ≤ c := by positivity
  have hsum' : ∑ ℓ ∈ Icc 1 L, G ℓ ≤ δ * L * c := by
    have : δ * L * c = L * (2 / g ^ 2 * X) := by
      rw [hc]
      field_simp
    rw [this]
    exact hsum
  set T := (Icc 1 L).filter (fun ℓ =>
    1 / 2 * Real.log (2 * π * exp 1 * (c + 1 / 12)) < shannon L (numDist (segment ℓ) ψ)) with hT
  have hTG : ∀ ℓ ∈ T, c < G ℓ := by
    intro ℓ hℓ
    rw [hT, Finset.mem_filter] at hℓ
    by_contra hle
    have h1 := WilsonAux.shannon_le_cutField θ e₀ Q ψ hnorm ℓ hℓ.1
    have h2 := WilsonAux.phi_mono (hG0 ℓ) (not_lt.1 hle)
    linarith [hℓ.2]
  by_cases hne : T.Nonempty
  · have hlt : ∑ _ℓ ∈ T, c < ∑ ℓ ∈ T, G ℓ := Finset.sum_lt_sum_of_nonempty hne hTG
    have hsub : ∑ ℓ ∈ T, G ℓ ≤ ∑ ℓ ∈ Icc 1 L, G ℓ :=
      Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _) (fun ℓ _ _ => hG0 ℓ)
    rw [Finset.sum_const, nsmul_eq_mul] at hlt
    have h1 : (T.card : ℝ) * c < δ * L * c := by linarith
    have hcpos : 0 < c := by
      rcases hc0.lt_or_eq with h | h
      · exact h
      · rw [← h, mul_zero, mul_zero] at h1
        exact absurd h1 (lt_irrefl 0)
    exact (lt_of_mul_lt_mul_right h1 hcpos.le).le
  · rw [Finset.not_nonempty_iff_eq_empty.1 hne, Finset.card_empty, Nat.cast_zero]
    positivity

/-- **Every cut** (paper: Proposition `prop:wilson`(e)). -/
theorem wilson_numberEntropy_each (w m g θ : ℝ) (hw : 0 ≤ w) (hg : 0 < g) (e₀ : ℤ)
    (Q : Fin L → ℤ) (hL : 0 < L) (ψ₀ : Cfg L → ℂ) (hψ₀ : ∑ n, Complex.normSq (ψ₀ n) = 1)
    (t : ℝ) (ℓ : ℕ) (hℓ : ℓ ∈ Icc 1 L) :
    shannon L (numDist (segment ℓ) (wState w m g θ e₀ Q ψ₀ t)) ≤
      1 / 2 * Real.log (2 * π * exp 1 *
        (2 * L * (energyDensity w m g θ e₀ Q ψ₀ + w + |m|) / g ^ 2 + 1 / 12)) := by
  have hnorm := WilsonAux.sum_normSq_wState w m g θ e₀ Q ψ₀ hψ₀ t
  have hsum := WilsonAux.sum_fieldSq_le w m g θ hw hg e₀ Q hL ψ₀ hψ₀ t
  rw [← WilsonAux.sum_cutField] at hsum
  refine (WilsonAux.shannon_le_cutField θ e₀ Q _ hnorm ℓ hℓ).trans
    (WilsonAux.phi_mono (WilsonAux.cutField_nonneg θ e₀ Q _ ℓ) ?_)
  calc WilsonAux.cutField θ e₀ Q (wState w m g θ e₀ Q ψ₀ t) ℓ
      ≤ ∑ ℓ' ∈ Icc 1 L, WilsonAux.cutField θ e₀ Q (wState w m g θ e₀ Q ψ₀ t) ℓ' :=
        Finset.single_le_sum (fun i _ => WilsonAux.cutField_nonneg θ e₀ Q _ i) hℓ
    _ ≤ L * (2 / g ^ 2 * (energyDensity w m g θ e₀ Q ψ₀ + w + |m|)) := hsum
    _ = 2 * L * (energyDensity w m g θ e₀ Q ψ₀ + w + |m|) / g ^ 2 := by ring

end DFL
