import DFL.SlaterBulk

/-!
# Quasi-free density matrices (paper: Lemma `lem:quasifree`, and `lem:slater`(a): `ρ_A = ρ[C_A]`)

For a matrix `C` on `ℓ²(Fin L)` let `Γ(C)` be its second quantization on the Fock space `𝓕_L`,
`Γ(C) c†_{x_1} ⋯ c†_{x_k} |0⟩ = c†(C δ_{x_1}) ⋯ c†(C δ_{x_k}) |0⟩` (`secondQ`), and (eq. `eq:main-rhoformula`)
`ρ[C] = ∑_{j=0}^{L} ((-1)^j / j!) ∑_{x_1, …, x_j} c_{x_1} ⋯ c_{x_j} Γ(C) c†_{x_j} ⋯ c†_{x_1}` (`qfState`).

* (a) For Hermitian `C` with orthonormal eigenvectors `v_i` and eigenvalues `ν_i`,
  `ρ[C] = ∏_i (ν_i c†(v_i)c(v_i) + (1 - ν_i) c(v_i)c†(v_i))` (`qfState_eq_prod`): `ρ[C]` does not depend on the
  eigenbasis; its matrix elements are polynomials in the entries of `C` (`continuous_qfState`), and `C ↦ ρ[C]` is
  uniformly continuous on `{0 ≤ C ≤ 1}` (`uniformContinuousOn_qfState`).
* (b) For `0 ≤ C ≤ 1`, `ρ[C]` is a density matrix (`qfState_posSemidef`, `trace_qfState`) with
  `S(ρ[C]) = ∑_i h(ν_i)` (`vonNeumann_qfState`); its generating function is
  `Tr(ρ[C] z^{N}) = det(1 + (z - 1) C)` (`trace_qfState_zN`).
* `lem:slater`(a): the reduced state of a Slater determinant on any region `R` is the quasi-free state with symbol
  `C_R` (`rdmF_eq_qfState`), after relabelling the modes of `R` increasingly by `Fin |R|`.
-/

open Matrix Finset
open scoped ComplexOrder

noncomputable section

namespace DFL

open AreaLawAux

variable {L : ℕ}

/-- The occupied sites of `n`. -/
def occSet (n : Cfg L) : Finset (Fin L) := univ.filter (fun x => n x = 1)

/-- **Second quantization** `Γ(C)`: `Γ(C)|n⟩ = c†(C δ_{x_1}) ⋯ c†(C δ_{x_k}) |0⟩` for the occupied sites
`x_1 < ⋯ < x_k` of `n`. -/
def secondQ (C : Matrix (Fin L) (Fin L) ℂ) : Matrix (Cfg L) (Cfg L) ℂ := fun n' n =>
  slaterAlt L (occSet n).card
    (fun j : Fin (occSet n).card => fun x => C x ((occSet n).orderEmbOfFin rfl j)) n'

/-- **The quasi-free density matrix** `ρ[C]` (paper: eq. `eq:main-rhoformula`). -/
def qfState (C : Matrix (Fin L) (Fin L) ℂ) : Matrix (Cfg L) (Cfg L) ℂ :=
  ∑ j ∈ range (L + 1), ((-1 : ℂ) ^ j / (j.factorial : ℂ)) •
    ∑ x : Fin j → Fin L,
      (List.ofFn fun i => cplx (ann (x i))).prod * secondQ C *
        (List.ofFn fun i => cplx (cre (x (Fin.rev i)))).prod

/-- The operator `z^{N}` on `𝓕_L`. -/
def zPowN (z : ℂ) : Matrix (Cfg L) (Cfg L) ℂ := diagonal (fun n => z ^ (occSet n).card)

namespace QuasiFreeAux

open SlaterExactAux SlaterGeneralAux

/-! ### Second quantization of Slater determinants -/

/-- `c†(δ_x) = c_x†`. -/
lemma crOp_single (x : Fin L) : crOp (Pi.single x (1 : ℂ)) = cplx (cre x) := by
  simp [crOp, Pi.single_apply]

/-- The occupied set of the Fock configuration of `S` is `S`. -/
lemma occSet_fockCfg (S : Finset (Fin L)) : occSet (fockCfg S) = S := by
  ext x; simp [occSet, FreeDynAux.fockCfg_eq_one]

/-- The image of `Fin.cons a s`. -/
lemma image_cons {α : Type*} [DecidableEq α] {m : ℕ} (a : α) (s : Fin m → α) :
    univ.image (Fin.cons a s : Fin (m + 1) → α) = insert a (univ.image s) := by
  ext x
  simp only [mem_image, mem_univ, true_and, mem_insert, Fin.exists_fin_succ, Fin.cons_zero,
    Fin.cons_succ]
  constructor
  · rintro (h | h)
    · exact Or.inl h.symm
    · exact Or.inr h
  · rintro (h | h)
    · exact Or.inl h.symm
    · exact Or.inr h

/-- A Slater determinant of increasing basis orbitals is a Fock basis vector. -/
lemma slaterAlt_single_of_strictMono : ∀ {m : ℕ} (u : Fin m → Fin L), StrictMono u →
    slaterAlt L m (fun j => Pi.single (u j) (1 : ℂ)) = fockVec (univ.image u)
  | 0, u, _ => by simp [slaterAlt_apply]
  | m + 1, u, hu => by
    have hu' : StrictMono (fun k : Fin m => u k.succ) := hu.comp (Fin.strictMono_succ)
    have himg : univ.image u = insert (u 0) (univ.image (fun k : Fin m => u k.succ)) := by
      have h := image_cons (u 0) (fun k : Fin m => u k.succ)
      have e : (Fin.cons (u 0) (fun k : Fin m => u k.succ) : Fin (m + 1) → Fin L) = u :=
        Fin.cons_self_tail u
      rwa [e] at h
    rw [slaterAlt_succ, slaterAlt_single_of_strictMono _ hu', crOp_single, himg,
      AreaLawAux.fockVec_insert]
    intro x hx
    obtain ⟨k, -, rfl⟩ := mem_image.mp hx
    exact hu (Fin.succ_pos k)

/-- The column of `Γ(C)` at `n`, written via any finset equal to the occupied set of `n`. -/
lemma secondQ_apply_of_eq (C : Matrix (Fin L) (Fin L) ℂ) (n' n : Cfg L) (T : Finset (Fin L))
    (hT : occSet n = T) :
    secondQ C n' n = slaterAlt L T.card (fun j x => C x (T.orderEmbOfFin rfl j)) n' := by
  subst hT; rfl

/-- `Γ(C) |S⟩ = c†(C δ_{x_1}) ⋯ c†(C δ_{x_k}) |0⟩`. -/
lemma secondQ_mulVec_fockVec (C : Matrix (Fin L) (Fin L) ℂ) (S : Finset (Fin L)) :
    secondQ C *ᵥ fockVec S = slaterAlt L S.card (fun j x => C x (S.orderEmbOfFin rfl j)) := by
  funext n'
  rw [fockVec, Matrix.mulVec_single_one]
  exact secondQ_apply_of_eq C n' _ S (occSet_fockCfg S)

/-- Changing the length argument of `orderEmbOfFin` along an equality. -/
lemma slaterAlt_orderEmb_cast (C : Matrix (Fin L) (Fin L) ℂ) (T : Finset (Fin L)) {m : ℕ}
    (h : T.card = m) :
    slaterAlt L T.card (fun j x => C x (T.orderEmbOfFin rfl j)) =
      slaterAlt L m (fun j x => C x (T.orderEmbOfFin h j)) := by
  subst h; rfl

/-- **Multiplicativity of the second quantization**:
`Γ(C) c†(w_1) ⋯ c†(w_m) |0⟩ = c†(C w_1) ⋯ c†(C w_m) |0⟩`. -/
lemma secondQ_mulVec_slaterAlt (C : Matrix (Fin L) (Fin L) ℂ) {m : ℕ} (w : Fin m → Fin L → ℂ) :
    secondQ C *ᵥ slaterAlt L m w = slaterAlt L m (fun j => C *ᵥ w j) := by
  have hFG : (secondQ C).mulVecLin.compAlternatingMap (slaterAlt L m) =
      (slaterAlt L m).compLinearMap C.mulVecLin := by
    refine Module.Basis.ext_alternating (Pi.basisFun ℂ (Fin L)) fun v hv => ?_
    set σ := Tuple.sort v
    have hmono : StrictMono (v ∘ σ) :=
      (Tuple.monotone_sort v).strictMono_of_injective (hv.comp σ.injective)
    rw [AlternatingMap.map_congr_perm _ _ σ,
      AlternatingMap.map_congr_perm ((slaterAlt L m).compLinearMap C.mulVecLin) _ σ]
    congr 1
    have hc : (univ.image (v ∘ σ)).card = m := by
      rw [card_image_of_injective _ hmono.injective, card_univ, Fintype.card_fin]
    have hs := Finset.orderEmbOfFin_unique hc (fun i => mem_image_of_mem (v ∘ σ) (mem_univ i))
      hmono
    have e1 : ((fun i => (Pi.basisFun ℂ (Fin L)) (v i)) ∘ σ) =
        fun j => Pi.single ((v ∘ σ) j) (1 : ℂ) := by
      funext j; simp
    rw [LinearMap.compAlternatingMap_apply, AlternatingMap.compLinearMap_apply, e1,
      slaterAlt_single_of_strictMono (v ∘ σ) hmono, Matrix.mulVecLin_apply,
      secondQ_mulVec_fockVec, slaterAlt_orderEmb_cast C _ hc, ← hs]
    congr 1
    funext j x
    simp
  exact congrArg (fun f => f w) hFG


/-! ### Slater determinants of an orthonormal family of eigenvectors -/

section Eig

variable (v : Fin L → Fin L → ℂ)

/-- The Slater determinant `c†(v_{s_0}) ⋯ c†(v_{s_{m-1}}) |0⟩` of orbitals of the family `v`. -/
def slat {m : ℕ} (s : Fin m → Fin L) : Cfg L → ℂ := slaterAlt L m (fun j => v (s j))

/-- `c†(v_a)` prepends the orbital `v_a`. -/
lemma crOp_mulVec_slat {m : ℕ} (s : Fin m → Fin L) (a : Fin L) :
    crOp (v a) *ᵥ slat v s = slat v (Fin.cons a s : Fin (m + 1) → Fin L) := by
  rw [slat, slat, slaterAlt_succ]
  simp only [Fin.cons_zero, Fin.cons_succ]

/-- A Slater determinant with a repeated orbital vanishes. -/
lemma slat_cons_eq_zero {m : ℕ} (s : Fin m → Fin L) (a : Fin L) (ha : a ∈ univ.image s) :
    slat v (Fin.cons a s : Fin (m + 1) → Fin L) = 0 := by
  obtain ⟨k, -, hk⟩ := mem_image.mp ha
  exact (slaterAlt L (m + 1)).map_eq_zero_of_eq _ (i := 0) (j := k.succ) (by simp [hk])
    (Fin.succ_ne_zero k).symm

/-- `Γ(C)` acts on a Slater determinant of eigenvectors by the product of the eigenvalues. -/
lemma secondQ_mulVec_slat (C : Matrix (Fin L) (Fin L) ℂ) (ν : Fin L → ℂ)
    (hCv : ∀ i, C *ᵥ v i = ν i • v i) {m : ℕ} (s : Fin m → Fin L) :
    secondQ C *ᵥ slat v s = (∏ k, ν (s k)) • slat v s := by
  rw [slat, secondQ_mulVec_slaterAlt]
  simp only [hCv]
  exact AlternatingMap.map_smul_univ _ _ _

variable {v}
variable (hv : ∀ i j, star (v i) ⬝ᵥ v j = if i = j then 1 else 0)
include hv

/-- `c(v_a)` kills a Slater determinant of the other orbitals. -/
lemma adj_mulVec_slat {m : ℕ} (s : Fin m → Fin L) (a : Fin L) (ha : a ∉ univ.image s) :
    (crOp (v a))ᴴ *ᵥ slat v s = 0 :=
  adj_slaterAlt _ _ fun k => by
    rw [hv, ite_eq_right_iff]
    intro h
    exact absurd (mem_image.mpr ⟨k, mem_univ k, h.symm⟩) ha

/-- CAR for a normalized orbital. -/
lemma car_self (a : Fin L) :
    (crOp (v a))ᴴ * crOp (v a) + crOp (v a) * (crOp (v a))ᴴ = 1 := by
  rw [crOp_CAR_gen, hv, ite_eq_left rfl, one_smul]

/-- `c(v_a) c†(v_a) |s⟩ = [a ∉ s] |s⟩`, written with the prepended orbital. -/
lemma adj_mulVec_slat_cons {m : ℕ} (s : Fin m → Fin L) (a : Fin L) :
    (crOp (v a))ᴴ *ᵥ slat v (Fin.cons a s : Fin (m + 1) → Fin L) =
      if a ∈ univ.image s then 0 else slat v s := by
  split_ifs with ha
  · rw [slat_cons_eq_zero v s a ha, mulVec_zero]
  · rw [← crOp_mulVec_slat, mulVec_mulVec]
    have h := car_self hv a
    rw [← eq_sub_iff_add_eq] at h
    rw [h, sub_mulVec, one_mulVec, ← mulVec_mulVec, adj_mulVec_slat hv s a ha, mulVec_zero,
      sub_zero]

/-- The number operator `n(v_a) = c†(v_a) c(v_a)` on a Slater determinant of the family. -/
lemma num_mulVec_slat {m : ℕ} (s : Fin m → Fin L) (a : Fin L) :
    (crOp (v a) * (crOp (v a))ᴴ) *ᵥ slat v s = if a ∈ univ.image s then slat v s else 0 := by
  split_ifs with ha
  · have h := car_self hv a
    rw [← eq_sub_iff_add_eq'] at h
    rw [h, sub_mulVec, one_mulVec, ← mulVec_mulVec, crOp_mulVec_slat,
      slat_cons_eq_zero v s a ha, mulVec_zero, sub_zero]
  · rw [← mulVec_mulVec, adj_mulVec_slat hv s a ha, mulVec_zero]

/-- The hole operator `1 - n(v_a) = c(v_a) c†(v_a)` on a Slater determinant of the family. -/
lemma hole_mulVec_slat {m : ℕ} (s : Fin m → Fin L) (a : Fin L) :
    ((crOp (v a))ᴴ * crOp (v a)) *ᵥ slat v s = if a ∈ univ.image s then 0 else slat v s := by
  rw [← mulVec_mulVec, crOp_mulVec_slat, adj_mulVec_slat_cons hv]

end Eig

/-! ### The operators `Ξ_j` -/

/-- `Ξ_j(Y) = ∑_{x_1, …, x_j} c_{x_1} ⋯ c_{x_j} Y c†_{x_j} ⋯ c†_{x_1}`, recursively. -/
def xiOp (Y : Matrix (Cfg L) (Cfg L) ℂ) : ℕ → Matrix (Cfg L) (Cfg L) ℂ
  | 0 => Y
  | j + 1 => ∑ x : Fin L, cplx (ann x) * xiOp Y j * cplx (cre x)

/-- The explicit sum of eq. `eq:main-rhoformula` is `Ξ_j`. -/
lemma sum_list_eq_xiOp (Y : Matrix (Cfg L) (Cfg L) ℂ) : ∀ j : ℕ,
    ∑ x : Fin j → Fin L, (List.ofFn fun i => cplx (ann (x i))).prod * Y *
      (List.ofFn fun i => cplx (cre (x (Fin.rev i)))).prod = xiOp Y j
  | 0 => by simp [xiOp]
  | j + 1 => by
    rw [← (Fin.consEquiv fun _ => Fin L).sum_comp, Fintype.sum_prod_type]
    simp only [Fin.consEquiv_apply, xiOp, ← sum_list_eq_xiOp Y j, Finset.mul_sum,
      Finset.sum_mul]
    refine Finset.sum_congr rfl fun a _ => Finset.sum_congr rfl fun x _ => ?_
    rw [List.ofFn_succ, List.prod_cons, List.ofFn_succ', List.prod_concat]
    simp only [Fin.cons_zero, Fin.cons_succ, Fin.rev_castSucc, Fin.rev_last]
    simp only [mul_assoc]

/-- **Completeness**: `∑_x c_x Z c†_x = ∑_i c(v_i) Z c†(v_i)` for an orthonormal basis `v`. -/
lemma sum_ann_cre_eq {v : Fin L → Fin L → ℂ}
    (hcomp : ∀ x y, ∑ i, star (v i x) * v i y = if x = y then 1 else 0)
    (Z : Matrix (Cfg L) (Cfg L) ℂ) :
    ∑ x : Fin L, cplx (ann x) * Z * cplx (cre x) = ∑ i, (crOp (v i))ᴴ * Z * crOp (v i) := by
  calc ∑ x : Fin L, cplx (ann x) * Z * cplx (cre x)
      = ∑ x : Fin L, ∑ y : Fin L, (if x = y then (1 : ℂ) else 0) •
          (cplx (ann x) * Z * cplx (cre y)) := by
        refine Finset.sum_congr rfl fun x _ => ?_
        simp [ite_smul, Finset.sum_ite_eq]
    _ = ∑ x : Fin L, ∑ y : Fin L, ∑ i, (star (v i x) * v i y) •
          (cplx (ann x) * Z * cplx (cre y)) := by
        simp only [← hcomp, Finset.sum_smul]
    _ = ∑ i, ∑ x : Fin L, ∑ y : Fin L, (star (v i x) * v i y) •
          (cplx (ann x) * Z * cplx (cre y)) := by
        rw [Finset.sum_congr rfl fun x _ => Finset.sum_comm]
        exact Finset.sum_comm
    _ = ∑ i, (crOp (v i))ᴴ * Z * crOp (v i) := by
        refine Finset.sum_congr rfl fun i _ => ?_
        rw [crOp_conjTranspose]
        simp only [crOp, Finset.sum_mul, Finset.mul_sum, smul_mul_assoc, mul_smul_comm,
          Finset.smul_sum, smul_smul]
        rw [Finset.sum_comm]
        refine Finset.sum_congr rfl fun y _ => Finset.sum_congr rfl fun x _ => ?_
        rw [mul_comm (star (v i x))]

/-! ### The inclusion–exclusion numbers -/

/-- `S_0(U) = 1`, `S_{j+1}(U) = ∑_{i ∉ U} ν_i S_j(U ∪ {i})`; `S_j(U) = j! e_j(ν|_{U^c})`. -/
def sFun (ν : Fin L → ℂ) : ℕ → Finset (Fin L) → ℂ
  | 0, _ => 1
  | j + 1, U => ∑ i ∈ univ \ U, ν i * sFun ν j (insert i U)

open Polynomial in
/-- `S_j(U)` is the `j`-th derivative at `0` of `∏_{i ∉ U} (1 + ν_i X)`. -/
lemma eval_iterate_derivative_prod (ν : Fin L → ℂ) (j : ℕ) : ∀ U : Finset (Fin L),
    (derivative^[j] (∏ i ∈ univ \ U, (C (ν i) * X + C 1))).eval 0 = sFun ν j U := by
  induction j with
  | zero => intro U; simp [sFun, eval_prod]
  | succ j ih =>
    intro U
    have hd : derivative (∏ i ∈ univ \ U, (C (ν i) * X + C 1)) =
        ∑ i ∈ univ \ U, C (ν i) * ∏ k ∈ univ \ insert i U, (C (ν k) * X + C 1) := by
      rw [derivative_prod_finset]
      refine Finset.sum_congr rfl fun i hi => ?_
      have he : (univ \ U).erase i = univ \ insert i U := by
        ext k
        simp only [mem_erase, mem_sdiff, mem_univ, true_and, mem_insert, not_or]
      rw [he, mul_comm]
      congr 1
      simp
    rw [Function.iterate_succ_apply, hd, iterate_derivative_sum, eval_finsetSum, sFun]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [iterate_derivative_C_mul, eval_mul, eval_C, ih]

open Polynomial in
/-- **Inclusion–exclusion**: `∑_{j ≤ L} (-1)^j / j! S_j(U) = ∏_{i ∉ U} (1 - ν_i)`. -/
lemma sum_sFun (ν : Fin L → ℂ) (U : Finset (Fin L)) :
    ∑ j ∈ range (L + 1), ((-1 : ℂ) ^ j / (j.factorial : ℂ)) * sFun ν j U =
      ∏ i ∈ univ \ U, (1 - ν i) := by
  set P : ℂ[X] := ∏ i ∈ univ \ U, (C (ν i) * X + C 1) with hP
  have hdeg : P.natDegree < L + 1 := by
    have h1 : P.natDegree ≤ ∑ i ∈ univ \ U, (C (ν i) * X + C 1).natDegree :=
      natDegree_prod_le _ _
    have h2 : ∑ i ∈ univ \ U, (C (ν i) * X + C 1).natDegree ≤ ∑ i ∈ univ \ U, 1 :=
      Finset.sum_le_sum fun i _ => natDegree_linear_le
    have h3 : ∑ i ∈ univ \ U, (1 : ℕ) ≤ L := by
      rw [Finset.sum_const, smul_eq_mul, mul_one]
      exact (card_le_univ _).trans (by simp)
    omega
  have hev : P.eval (-1) = ∏ i ∈ univ \ U, (1 - ν i) := by
    rw [hP, eval_prod]
    refine Finset.prod_congr rfl fun i _ => ?_
    simp only [eval_add, eval_mul, eval_C, eval_X]
    ring
  rw [← hev, eval_eq_sum_range' hdeg]
  refine Finset.sum_congr rfl fun j _ => ?_
  have hc := coeff_iterate_derivative (k := j) P 0
  rw [zero_add, Nat.descFactorial_self, coeff_zero_eq_eval_zero,
    eval_iterate_derivative_prod, nsmul_eq_mul] at hc
  have hf : (j.factorial : ℂ) ≠ 0 := by exact_mod_cast Nat.factorial_ne_zero j
  rw [hc]
  field_simp

/-- `∑_i [i ∉ U] f_i = ∑_{i ∈ U^c} f_i`. -/
lemma sum_ite_mem_zero {M : Type*} [AddCommMonoid M] (U : Finset (Fin L)) (f : Fin L → M) :
    ∑ i, (if i ∈ U then 0 else f i) = ∑ i ∈ univ \ U, f i := by
  rw [Finset.sdiff_eq_filter, Finset.sum_filter]
  exact Finset.sum_congr rfl fun i _ => by by_cases h : i ∈ U <;> simp [h]

/-! ### Action of `ρ[C]` on Slater determinants of eigenvectors -/

section Main

variable {v : Fin L → Fin L → ℂ}
  (hv : ∀ i j, star (v i) ⬝ᵥ v j = if i = j then 1 else 0)
  (hcomp : ∀ x y, ∑ i, star (v i x) * v i y = if x = y then 1 else 0)
  {C : Matrix (Fin L) (Fin L) ℂ} {ν : Fin L → ℂ} (hCv : ∀ i, C *ᵥ v i = ν i • v i)
include hv hcomp hCv

/-- `Ξ_j(Γ(C)) |s⟩ = ∏_k ν_{s_k} S_j(s) |s⟩` for a Slater determinant of eigenvectors. -/
lemma xiOp_mulVec_slat (j : ℕ) : ∀ {m : ℕ} (s : Fin m → Fin L),
    xiOp (secondQ C) j *ᵥ slat v s = ((∏ k, ν (s k)) * sFun ν j (univ.image s)) • slat v s := by
  induction j with
  | zero =>
    intro m s
    simp only [xiOp, sFun, mul_one]
    exact secondQ_mulVec_slat v C ν hCv s
  | succ j ih =>
    intro m s
    have h1 : ∀ i, ((crOp (v i))ᴴ * xiOp (secondQ C) j * crOp (v i)) *ᵥ slat v s =
        (if i ∈ univ.image s then 0 else
          ν i * ((∏ k, ν (s k)) * sFun ν j (insert i (univ.image s)))) • slat v s := by
      intro i
      rw [← mulVec_mulVec, ← mulVec_mulVec, crOp_mulVec_slat, ih, mulVec_smul,
        adj_mulVec_slat_cons hv, Fin.prod_univ_succ, image_cons]
      simp only [Fin.cons_zero, Fin.cons_succ]
      split_ifs <;> simp [mul_assoc]
    rw [xiOp, sum_ann_cre_eq hcomp, sum_mulVec]
    simp only [h1, ← Finset.sum_smul]
    congr 1
    rw [sum_ite_mem_zero, sFun, Finset.mul_sum]
    exact Finset.sum_congr rfl fun i _ => by ring

/-- **The quasi-free state on Slater determinants of eigenvectors** (paper: proof of
Lemma `lem:quasifree`(a)): `ρ[C] |s⟩ = ∏_k ν_{s_k} ∏_{i ∉ s} (1 - ν_i) |s⟩`. -/
lemma qfState_mulVec_slat {m : ℕ} (s : Fin m → Fin L) :
    qfState C *ᵥ slat v s =
      ((∏ k, ν (s k)) * ∏ i ∈ univ \ univ.image s, (1 - ν i)) • slat v s := by
  rw [qfState, sum_mulVec]
  simp only [smul_mulVec, sum_list_eq_xiOp, xiOp_mulVec_slat hv hcomp hCv, smul_smul,
    ← Finset.sum_smul]
  rw [← sum_sFun ν (univ.image s), Finset.mul_sum]
  congr 1
  exact Finset.sum_congr rfl fun j _ => by ring

end Main

/-! ### The Slater eigenbasis -/

/-- The Slater basis vector `|T⟩ = ∏_{i ∈ T} c†(v_i) |0⟩` (increasing order). -/
def bvec (v : Fin L → Fin L → ℂ) (T : Finset (Fin L)) : Cfg L → ℂ := slat v (T.orderEmbOfFin rfl)

/-- The increasing enumeration of `T` has image `T`. -/
lemma image_orderEmb (T : Finset (Fin L)) : univ.image (T.orderEmbOfFin rfl) = T :=
  image_orderEmbOfFin_univ T rfl

/-- `∏_k f(t_k) = ∏_{i ∈ T} f(i)` along the increasing enumeration. -/
lemma prod_orderEmb {M : Type*} [CommMonoid M] (T : Finset (Fin L)) (f : Fin L → M) :
    ∏ k, f ((T.orderEmbOfFin rfl) k) = ∏ i ∈ T, f i := by
  conv_rhs => rw [← image_orderEmb T]
  rw [Finset.prod_image fun x _ y _ h => (T.orderEmbOfFin rfl).injective h]

/-- The matrix of the Slater basis vectors (columns). -/
def bmat (v : Fin L → Fin L → ℂ) : Matrix (Cfg L) (Finset (Fin L)) ℂ :=
  Matrix.of fun n T => bvec v T n

/-- `#Cfg L = #Finset (Fin L) = 2^L`. -/
lemma card_cfg_eq : Fintype.card (Cfg L) = Fintype.card (Finset (Fin L)) := by
  simp [Fintype.card_finset, ZMod.card]

section Basis

variable {v : Fin L → Fin L → ℂ} (hv : ∀ i j, star (v i) ⬝ᵥ v j = if i = j then 1 else 0)
include hv

/-- **Orthonormality of the Slater basis**. -/
lemma bvec_orth (T T' : Finset (Fin L)) :
    star (bvec v T) ⬝ᵥ bvec v T' = if T = T' then 1 else 0 := by
  split_ifs with h
  · subst h
    rw [bvec, slat, gram_slaterAlt_diag]
    · exact Finset.prod_eq_one fun k _ => by rw [hv, ite_eq_left rfl]
    · intro k l hkl
      rw [hv, ite_eq_right_iff]
      intro h
      exact absurd ((T.orderEmbOfFin rfl).injective h) hkl
  · obtain ⟨i, hi⟩ : ∃ i, ¬ (i ∈ T ↔ i ∈ T') := by
      by_contra hc
      exact h (Finset.ext fun k => not_not.mp fun hk => hc ⟨k, hk⟩)
    by_cases hiT : i ∈ T
    · have hiT' : i ∉ T' := fun h' => hi ⟨fun _ => h', fun _ => hiT⟩
      have h1 := num_mulVec_slat hv (T.orderEmbOfFin rfl) i
      rw [image_orderEmb, ite_eq_left hiT] at h1
      rw [bvec, ← h1, ← mulVec_mulVec, dot_crOp_left, bvec,
        adj_mulVec_slat hv (T'.orderEmbOfFin rfl) i (by rwa [image_orderEmb]), dotProduct_zero]
    · have hiT' : i ∈ T' := by tauto
      have h1 := num_mulVec_slat hv (T'.orderEmbOfFin rfl) i
      rw [image_orderEmb, ite_eq_left hiT'] at h1
      rw [bvec, bvec, ← h1, ← mulVec_mulVec, dot_crOp_right,
        adj_mulVec_slat hv (T.orderEmbOfFin rfl) i (by rwa [image_orderEmb]), star_zero,
        zero_dotProduct]

/-- The Slater basis matrix is an isometry. -/
lemma bmat_ct_mul : (bmat v)ᴴ * bmat v = 1 := by
  refine Matrix.ext fun T T' => ?_
  rw [Matrix.mul_apply, Matrix.one_apply, ← bvec_orth hv T T']
  simp [bmat, dotProduct, conjTranspose_apply]

/-- The Slater basis matrix is unitary. -/
lemma bmat_mul_ct : bmat v * (bmat v)ᴴ = 1 :=
  (Matrix.mul_eq_one_comm_of_equiv (Fintype.equivOfCardEq card_cfg_eq.symm)).mp (bmat_ct_mul hv)

/-- An operator acting diagonally on the Slater basis is `W diag(d) W†`. -/
lemma eq_of_mulVec_bvec (X : Matrix (Cfg L) (Cfg L) ℂ) (d : Finset (Fin L) → ℂ)
    (h : ∀ T, X *ᵥ bvec v T = d T • bvec v T) : X = bmat v * diagonal d * (bmat v)ᴴ := by
  have hXW : X * bmat v = bmat v * diagonal d := by
    refine Matrix.ext fun n T => ?_
    rw [mul_diagonal]
    have := congrFun (h T) n
    simp only [mulVec, dotProduct, Pi.smul_apply, smul_eq_mul] at this
    rw [Matrix.mul_apply]
    simp only [bmat, of_apply]
    rw [this, mul_comm]
  calc X = X * (bmat v * (bmat v)ᴴ) := by rw [bmat_mul_ct hv, mul_one]
    _ = bmat v * diagonal d * (bmat v)ᴴ := by rw [← Matrix.mul_assoc, hXW]

end Basis

/-- The weights `w_T = ∏_{i ∈ T} ν_i ∏_{i ∉ T} (1 - ν_i)`. -/
def wT (ν : Fin L → ℂ) (T : Finset (Fin L)) : ℂ := ∏ i, if i ∈ T then ν i else 1 - ν i

/-- `w_T` as a product over `T` and its complement. -/
lemma wT_eq (ν : Fin L → ℂ) (T : Finset (Fin L)) :
    wT ν T = (∏ i ∈ T, ν i) * ∏ i ∈ univ \ T, (1 - ν i) := by
  rw [wT, ← Finset.prod_mul_prod_compl T, Finset.compl_eq_univ_sdiff]
  congr 1
  · exact Finset.prod_congr rfl fun k hk => ite_eq_left hk
  · exact Finset.prod_congr rfl fun k hk => ite_eq_right (Finset.mem_sdiff.mp hk).2

/-- **Diagonalization of `ρ[C]`** in the Slater eigenbasis of any orthonormal eigenbasis of `C`:
`ρ[C] = W diag(w) W†`. -/
theorem qfState_eq_conj {v : Fin L → Fin L → ℂ}
    (hv : ∀ i j, star (v i) ⬝ᵥ v j = if i = j then 1 else 0)
    (hcomp : ∀ x y, ∑ i, star (v i x) * v i y = if x = y then 1 else 0)
    {C : Matrix (Fin L) (Fin L) ℂ} {ν : Fin L → ℂ} (hCv : ∀ i, C *ᵥ v i = ν i • v i) :
    qfState C = bmat v * diagonal (wT ν) * (bmat v)ᴴ :=
  eq_of_mulVec_bvec hv _ _ fun T => by
    rw [bvec, qfState_mulVec_slat hv hcomp hCv, image_orderEmb, wT_eq, prod_orderEmb]

/-! ### Mathlib's eigenbasis -/

section Hermitian

variable {C : Matrix (Fin L) (Fin L) ℂ} (hC : C.IsHermitian)
include hC

/-- Orthonormality of the eigenvectors. -/
lemma eig_orth (i j : Fin L) :
    star (⇑(hC.eigenvectorBasis i) : Fin L → ℂ) ⬝ᵥ ⇑(hC.eigenvectorBasis j) =
      if i = j then 1 else 0 := by
  have h := orthonormal_iff_ite.mp hC.eigenvectorBasis.orthonormal i j
  rw [EuclideanSpace.inner_eq_star_dotProduct] at h
  rw [dotProduct_comm]
  exact h

/-- Completeness of the eigenvectors. -/
lemma eig_comp (x y : Fin L) :
    ∑ i, star ((⇑(hC.eigenvectorBasis i) : Fin L → ℂ) x) * ⇑(hC.eigenvectorBasis i) y =
      if x = y then 1 else 0 := by
  have h := congrFun (congrFun (Unitary.coe_mul_star_self hC.eigenvectorUnitary) y) x
  rw [Matrix.mul_apply, Matrix.one_apply] at h
  simp only [IsHermitian.eigenvectorUnitary_apply] at h
  calc ∑ i, star ((⇑(hC.eigenvectorBasis i) : Fin L → ℂ) x) * ⇑(hC.eigenvectorBasis i) y
      = ∑ i, (hC.eigenvectorBasis i).ofLp y * star ((hC.eigenvectorBasis i).ofLp x) :=
        Finset.sum_congr rfl fun i _ => mul_comm _ _
    _ = if y = x then 1 else 0 := h
    _ = if x = y then 1 else 0 := by
        by_cases hxy : x = y
        · subst hxy; simp
        · simp [hxy, Ne.symm hxy]

/-- The eigenvalue equation with complex scalars. -/
lemma eig_eq (i : Fin L) :
    C *ᵥ ⇑(hC.eigenvectorBasis i) = ((hC.eigenvalues i : ℝ) : ℂ) • ⇑(hC.eigenvectorBasis i) := by
  rw [hC.mulVec_eigenvectorBasis]
  funext x
  simp [Complex.real_smul]

/-- `ρ[C]` in Mathlib's eigenbasis. -/
lemma qfState_eq_conj_eig :
    qfState C = bmat (fun i => ⇑(hC.eigenvectorBasis i)) *
      diagonal (wT (fun i => (hC.eigenvalues i : ℂ))) * (bmat (fun i => ⇑(hC.eigenvectorBasis i)))ᴴ :=
  qfState_eq_conj (eig_orth hC) (eig_comp hC) (eig_eq hC)

/-- `det(1 + (z - 1) C) = ∏_i (ν_i z + 1 - ν_i)`. -/
lemma det_one_add_smul_eig (z : ℂ) :
    (1 + (z - 1) • C).det = ∏ i, ((hC.eigenvalues i : ℂ) * z + (1 - hC.eigenvalues i)) := by
  set U : Matrix (Fin L) (Fin L) ℂ := (hC.eigenvectorUnitary : Matrix (Fin L) (Fin L) ℂ)
  have hU1 : U * star U = 1 := Unitary.coe_mul_star_self _
  have hU2 : star U * U = 1 := Unitary.coe_star_mul_self _
  have hspec : C = U * diagonal (RCLike.ofReal ∘ hC.eigenvalues) * star U := by
    have := hC.spectral_theorem
    rwa [Unitary.conjStarAlgAut_apply] at this
  have hD : (1 : Matrix (Fin L) (Fin L) ℂ) + (z - 1) • diagonal (RCLike.ofReal ∘ hC.eigenvalues) =
      diagonal (fun i => (hC.eigenvalues i : ℂ) * z + (1 - hC.eigenvalues i)) := by
    ext i j
    by_cases h : i = j
    · subst h; simp; ring
    · simp [h]
  have h1 : 1 + (z - 1) • C =
      U * (1 + (z - 1) • diagonal (RCLike.ofReal ∘ hC.eigenvalues)) * star U := by
    rw [Matrix.mul_add, Matrix.add_mul, Matrix.mul_one, hU1, Matrix.mul_smul, Matrix.smul_mul,
      ← hspec]
  have h2 : (star U).det * U.det = 1 := by rw [← det_mul, hU2, det_one]
  rw [h1, det_mul, det_mul, hD, det_diagonal]
  calc U.det * (∏ i, ((hC.eigenvalues i : ℂ) * z + (1 - hC.eigenvalues i))) * (star U).det
      = (∏ i, ((hC.eigenvalues i : ℂ) * z + (1 - hC.eigenvalues i))) * ((star U).det * U.det) := by
        ring
    _ = _ := by rw [h2, mul_one]

end Hermitian

/-- The eigenvalues of `C` are at most `1` if `1 - C ≥ 0`. -/
lemma eig_le_one {C : Matrix (Fin L) (Fin L) ℂ} (hC0 : C.PosSemidef) (hC1 : (1 - C).PosSemidef)
    (i : Fin L) : hC0.isHermitian.eigenvalues i ≤ 1 := by
  have h := hC1.dotProduct_mulVec_nonneg (⇑(hC0.isHermitian.eigenvectorBasis i))
  rw [sub_mulVec, one_mulVec, eig_eq hC0.isHermitian, dotProduct_sub, dotProduct_smul,
    eig_orth hC0.isHermitian, ite_eq_left rfl, smul_eq_mul, mul_one] at h
  have h' : (0 : ℂ) ≤ ((1 - hC0.isHermitian.eigenvalues i : ℝ) : ℂ) := by
    push_cast; exact h
  exact sub_nonneg.mp (Complex.zero_le_real.mp h')

/-- Real weights. -/
lemma wT_real (ν : Fin L → ℝ) (T : Finset (Fin L)) :
    wT (fun i => (ν i : ℂ)) T = ((∏ i, if i ∈ T then ν i else 1 - ν i : ℝ) : ℂ) := by
  rw [wT, Complex.ofReal_prod]
  exact Finset.prod_congr rfl fun i _ => by split_ifs <;> simp

/-- The weights sum to one. -/
lemma sum_wT (ν : Fin L → ℂ) : ∑ T, wT ν T = 1 := by
  have h := Finset.prod_add ν (fun i => 1 - ν i) univ
  rw [powerset_univ] at h
  rw [Finset.sum_congr rfl fun T _ => wT_eq ν T, ← h]
  simp

/-- Generating function of the weights. -/
lemma sum_wT_pow (ν : Fin L → ℂ) (z : ℂ) :
    ∑ T, wT ν T * z ^ T.card = ∏ i, (ν i * z + (1 - ν i)) := by
  have h := Finset.prod_add (fun i => ν i * z) (fun i => 1 - ν i) univ
  rw [powerset_univ] at h
  rw [h]
  refine Finset.sum_congr rfl fun T _ => ?_
  rw [wT_eq, Finset.prod_mul_distrib, Finset.prod_const]
  ring

/-- The product of commuting-type factors acting diagonally. -/
lemma list_prod_mulVec {ι : Type*} [Fintype ι] [DecidableEq ι] {n : ℕ} (F : Fin n → Matrix ι ι ℂ)
    (c : Fin n → ℂ) (φ : ι → ℂ) (h : ∀ i, F i *ᵥ φ = c i • φ) :
    (List.ofFn F).prod *ᵥ φ = (∏ i, c i) • φ := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [List.ofFn_succ, List.prod_cons, ← mulVec_mulVec,
      ih (fun i => F i.succ) (fun i => c i.succ) (fun i => h i.succ), mulVec_smul, h 0, smul_smul,
      Fin.prod_univ_succ, mul_comm]

/-- The factor `ν_i n_i + (1 - ν_i)(1 - n_i)` on the Slater basis. -/
lemma factor_mulVec_bvec {v : Fin L → Fin L → ℂ}
    (hv : ∀ i j, star (v i) ⬝ᵥ v j = if i = j then 1 else 0) (ν : Fin L → ℝ) (i : Fin L)
    (T : Finset (Fin L)) :
    (((ν i : ℝ) : ℂ) • (crOp (v i) * (crOp (v i))ᴴ) +
        ((1 - ν i : ℝ) : ℂ) • ((crOp (v i))ᴴ * crOp (v i))) *ᵥ bvec v T =
      (if i ∈ T then (ν i : ℂ) else 1 - ν i) • bvec v T := by
  rw [add_mulVec, smul_mulVec, smul_mulVec, bvec, num_mulVec_slat hv, hole_mulVec_slat hv,
    image_orderEmb]
  split_ifs <;> simp

/-- `N = #occupied` as a real number. -/
lemma count_univ_eq (n : Cfg L) : count univ n = ((occSet n).card : ℝ) := by
  rw [count, occSet, card_filter, Nat.cast_sum]
  refine Finset.sum_congr rfl fun x _ => ?_
  rcases FreeDynAux.zmod2_cases (n x) with h | h <;> simp [occ, h]

/-- `z^N |T⟩ = z^{|T|} |T⟩`. -/
lemma zPowN_mulVec_bvec (v : Fin L → Fin L → ℂ) (z : ℂ) (T : Finset (Fin L)) :
    zPowN z *ᵥ bvec v T = (z ^ T.card) • bvec v T := by
  funext n
  rw [zPowN, mulVec_diagonal, Pi.smul_apply, smul_eq_mul]
  by_cases hn : (occSet n).card = T.card
  · rw [hn]
  · have hnum := isNumR_slaterAlt univ (fun j => v (T.orderEmbOfFin rfl j)) (fun _ => True)
      (fun _ _ x hx => absurd (mem_univ x) hx) (fun _ h => absurd trivial h)
    have h0 : bvec v T n = 0 := hnum n (by rw [count_univ_eq]; simpa using hn)
    rw [h0, mul_zero, mul_zero]

/-- `c†(·)` is continuous. -/
lemma continuous_crOp : Continuous (crOp : (Fin L → ℂ) → Matrix (Cfg L) (Cfg L) ℂ) := by
  unfold crOp
  exact continuous_finsetSum _ fun x _ => (continuous_apply x).smul continuous_const

/-- Slater determinants depend continuously on the orbitals. -/
lemma continuous_slaterAlt_comp {X : Type*} [TopologicalSpace X] {m : ℕ}
    (w : Fin m → X → Fin L → ℂ) (hw : ∀ j, Continuous (w j)) :
    Continuous fun a => slaterAlt L m (fun j => w j a) := by
  simp only [slaterAlt_apply]
  refine Continuous.matrix_mulVec ?_ continuous_const
  simp only [List.ofFn_eq_map]
  exact continuous_list_prod _ fun j _ => continuous_crOp.comp (hw j)

/-- The second quantization is continuous. -/
lemma continuous_secondQ :
    Continuous (secondQ : Matrix (Fin L) (Fin L) ℂ → Matrix (Cfg L) (Cfg L) ℂ) := by
  refine continuous_matrix fun n' n => ?_
  have h := continuous_slaterAlt_comp (m := (occSet n).card)
    (fun j (C : Matrix (Fin L) (Fin L) ℂ) x => C x ((occSet n).orderEmbOfFin rfl j))
    (fun j => continuous_pi fun x => Continuous.matrix_elem continuous_id x _)
  exact (continuous_apply n').comp h

/-- Entries of `U diag(d) U†` with `U` unitary and `|d_k| ≤ 1` are bounded by the dimension. -/
lemma norm_conj_diag_le {n : ℕ} (U : Matrix (Fin n) (Fin n) ℂ)
    (hU : U ∈ Matrix.unitaryGroup (Fin n) ℂ) (d : Fin n → ℂ) (hd : ∀ k, ‖d k‖ ≤ 1)
    (i j : Fin n) : ‖(U * diagonal d * star U) i j‖ ≤ n := by
  rw [Matrix.mul_apply]
  refine (norm_sum_le _ _).trans ?_
  have hb : ∀ k, ‖(U * diagonal d) i k * star U k j‖ ≤ 1 := by
    intro k
    rw [Matrix.mul_diagonal, Matrix.star_apply, norm_mul, norm_mul, norm_star]
    have h1 := entry_norm_bound_of_unitary hU i k
    have h2 := entry_norm_bound_of_unitary hU j k
    have h3 := hd k
    have h4 : 0 ≤ ‖U i k‖ := norm_nonneg _
    have h5 : 0 ≤ ‖d k‖ := norm_nonneg _
    calc ‖U i k‖ * ‖d k‖ * ‖U j k‖ ≤ 1 * 1 * 1 := by gcongr
      _ = 1 := by ring
  calc ∑ k, ‖(U * diagonal d) i k * star U k j‖ ≤ ∑ k : Fin n, (1 : ℝ) :=
        Finset.sum_le_sum fun k _ => hb k
    _ = n := by simp

/-- Entries of `0 ≤ C ≤ 1` are bounded by `L`. -/
lemma norm_entry_le {C : Matrix (Fin L) (Fin L) ℂ} (hC0 : C.PosSemidef)
    (hC1 : (1 - C).PosSemidef) (i j : Fin L) : ‖C i j‖ ≤ L := by
  have hspec := hC0.isHermitian.spectral_theorem
  rw [Unitary.conjStarAlgAut_apply] at hspec
  rw [hspec]
  refine norm_conj_diag_le _ hC0.isHermitian.eigenvectorUnitary.2 _ (fun k => ?_) i j
  simp only [Function.comp_apply]
  rw [RCLike.norm_ofReal, abs_le]
  constructor <;> linarith [hC0.eigenvalues_nonneg k, eig_le_one hC0 hC1 k]

end QuasiFreeAux

open QuasiFreeAux

/-- **Eigenmode form** (paper: Lemma `lem:quasifree`(a)): `ρ[C] = ∏_i (ν_i c†(v_i)c(v_i) + (1 - ν_i) c(v_i)c†(v_i))`
for Hermitian `C` (the factors commute); in particular `ρ[C]` does not depend on the choice of eigenbasis. -/
theorem qfState_eq_prod (C : Matrix (Fin L) (Fin L) ℂ) (hC : C.IsHermitian) :
    qfState C = (List.ofFn fun i : Fin L =>
      ((hC.eigenvalues i : ℝ) : ℂ) • (crOp (hC.eigenvectorBasis i) * (crOp (hC.eigenvectorBasis i))ᴴ) +
        ((1 - hC.eigenvalues i : ℝ) : ℂ) •
          ((crOp (hC.eigenvectorBasis i))ᴴ * crOp (hC.eigenvectorBasis i))).prod := by
  rw [qfState_eq_conj_eig hC]
  symm
  refine eq_of_mulVec_bvec (eig_orth hC) _ _ fun T => ?_
  rw [list_prod_mulVec _ (fun i => if i ∈ T then (hC.eigenvalues i : ℂ) else 1 - hC.eigenvalues i) _
    (fun i => factor_mulVec_bvec (eig_orth hC) hC.eigenvalues i T)]
  rfl

/-- The matrix elements of `ρ[C]` are polynomials in the entries of `C`: `C ↦ ρ[C]` is continuous. -/
theorem continuous_qfState : Continuous (qfState : Matrix (Fin L) (Fin L) ℂ → Matrix (Cfg L) (Cfg L) ℂ) := by
  unfold qfState
  refine continuous_finsetSum _ fun j _ =>
    Continuous.fun_const_smul (continuous_finsetSum _ fun x _ => ?_) _
  exact (continuous_const.matrix_mul continuous_secondQ).matrix_mul continuous_const

/-- **Uniform continuity** (paper: Lemma `lem:quasifree`(a)) on `{0 ≤ C ≤ 1}`. -/
theorem uniformContinuousOn_qfState :
    UniformContinuousOn (qfState : Matrix (Fin L) (Fin L) ℂ → Matrix (Cfg L) (Cfg L) ℂ)
      {C | C.PosSemidef ∧ (1 - C).PosSemidef} := by
  refine IsCompact.uniformContinuousOn_of_continuous ?_ continuous_qfState.continuousOn
  have hclosed : IsClosed {C : Matrix (Fin L) (Fin L) ℂ | C.PosSemidef ∧ (1 - C).PosSemidef} :=
    isClosed_setOfPred_posSemidef.inter
      (isClosed_setOfPred_posSemidef.preimage (continuous_const.sub continuous_id))
  have hK : IsCompact (Set.univ.pi fun _ : Fin L => Set.univ.pi fun _ : Fin L =>
      Metric.closedBall (0 : ℂ) L) :=
    isCompact_univ_pi fun _ => isCompact_univ_pi fun _ => isCompact_closedBall _ _
  refine hK.of_isClosed_subset hclosed ?_
  rintro C ⟨hC0, hC1⟩
  simp only [Set.mem_pi, Set.mem_univ, true_implies, Metric.mem_closedBall, dist_zero_right]
  intro i j
  exact norm_entry_le hC0 hC1 i j

/-- `ρ[C] ≥ 0` for `0 ≤ C ≤ 1` (paper: Lemma `lem:quasifree`(b)). -/
theorem qfState_posSemidef (C : Matrix (Fin L) (Fin L) ℂ) (hC0 : C.PosSemidef) (hC1 : (1 - C).PosSemidef) :
    (qfState C).PosSemidef := by
  rw [qfState_eq_conj_eig hC0.isHermitian]
  refine PosSemidef.mul_mul_conjTranspose_same ?_ _
  rw [posSemidef_diagonal_iff]
  intro T
  rw [wT_real]
  refine Complex.zero_le_real.mpr (Finset.prod_nonneg fun i _ => ?_)
  split_ifs
  · exact hC0.eigenvalues_nonneg i
  · linarith [eig_le_one hC0 hC1 i]

/-- `Tr ρ[C] = 1` for Hermitian `C` (paper: Lemma `lem:quasifree`(b)). -/
theorem trace_qfState (C : Matrix (Fin L) (Fin L) ℂ) (hC : C.IsHermitian) : (qfState C).trace = 1 := by
  rw [qfState_eq_conj_eig hC, trace_mul_comm, ← Matrix.mul_assoc, bmat_ct_mul (eig_orth hC),
    Matrix.one_mul, trace_diagonal]
  exact sum_wT _

/-- **Entropy of a quasi-free state** (paper: Lemma `lem:quasifree`(b)): `S(ρ[C]) = ∑_i h(ν_i)`. -/
theorem vonNeumann_qfState (C : Matrix (Fin L) (Fin L) ℂ) (hC0 : C.PosSemidef) (hC1 : (1 - C).PosSemidef) :
    vonNeumann (qfState C) (qfState_posSemidef C hC0 hC1).isHermitian =
      ∑ i, Real.binEntropy (hC0.isHermitian.eigenvalues i) := by
  set hC := hC0.isHermitian
  set W := bmat (fun i => ⇑(hC.eigenvectorBasis i))
  set d : Finset (Fin L) → ℝ := fun T => ∏ i, if i ∈ T then hC.eigenvalues i else 1 - hC.eigenvalues i
  have hw : wT (fun i => (hC.eigenvalues i : ℂ)) = fun T => (d T : ℂ) :=
    funext fun T => wT_real _ T
  have heq : qfState C = W * diagonal (fun T => (d T : ℂ)) * Wᴴ := by
    rw [qfState_eq_conj_eig hC, hw]
  have hdH : (diagonal (fun T => (d T : ℂ))).IsHermitian := by
    refine isHermitian_diagonal_of_self_adjoint _ (funext fun T => ?_)
    simp [Complex.conj_ofReal]
  have hWH : (W * diagonal (fun T => (d T : ℂ)) * Wᴴ).IsHermitian :=
    heq ▸ (qfState_posSemidef C hC0 hC1).isHermitian
  rw [CutLocalAux.vonNeumann_congr _ hWH heq,
    CutLocalAux.vonNeumann_conj_isometry W (bmat_ct_mul (eig_orth hC)) _ hdH hWH,
    vonNeumann_diagonal, ← SlaterExactAux.sum_negMulLog_prod_sel hC.eigenvalues,
    ← (SlaterExactAux.boolEquiv L).sum_comp]
  refine Finset.sum_congr rfl fun σ _ => ?_
  congr 1
  refine Finset.prod_congr rfl fun k _ => ?_
  simp [SlaterExactAux.boolEquiv]

/-- **Generating function** (paper: Lemma `lem:quasifree`(b)): `Tr(ρ[C] z^{N}) = det(1 + (z - 1) C)`. -/
theorem trace_qfState_zN (C : Matrix (Fin L) (Fin L) ℂ) (hC : C.IsHermitian) (z : ℂ) :
    (qfState C * zPowN z).trace = (1 + (z - 1) • C).det := by
  have hv := eig_orth hC
  have hWW := bmat_ct_mul hv
  have hZ := eq_of_mulVec_bvec hv (zPowN (L := L) z) (fun T : Finset (Fin L) => z ^ T.card)
    fun T => zPowN_mulVec_bvec _ z T
  rw [qfState_eq_conj_eig hC, hZ, det_one_add_smul_eig hC z, ← sum_wT_pow]
  simp only [Matrix.mul_assoc]
  rw [← Matrix.mul_assoc (bmat _)ᴴ (bmat _), hWW, Matrix.one_mul, trace_mul_comm,
    Matrix.mul_assoc, Matrix.mul_assoc, hWW, Matrix.mul_one, diagonal_mul_diagonal,
    trace_diagonal]

/-- **Reduced states of Slater determinants are quasi-free** (paper: Lemma `lem:slater`(a), `ρ_A = ρ[C_A]`): for any
region `R`, relabelled increasingly by `Fin |R|`. -/
theorem rdmF_eq_qfState {N : ℕ} (A : Matrix (Fin N) (Fin N) ℝ) (hA : A.IsSymm) (S : Finset (Fin N)) (t : ℝ)
    (R : Finset (Fin N)) :
    (rdmF R (evo (dGamma A) t *ᵥ fockVec S)).submatrix
        (fun n : Cfg R.card => fun a : R => n ((R.orderIsoOfFin rfl).symm a))
        (fun n : Cfg R.card => fun a : R => n ((R.orderIsoOfFin rfl).symm a)) =
      qfState ((corrR A S t R).submatrix (fun i => (R.orderIsoOfFin rfl) i) (fun i => (R.orderIsoOfFin rfl) i)) := by
  sorry

end DFL
