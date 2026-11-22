# Localization, slow dynamics and number entropy in lattice gauge theories: what can be proven (note, 2026-11-22)

## 0. Status, sources, conventions

**Network limits in this session.**
- arxiv.org did not resolve, and the shared WebSearch budget was already used up when I started.
- So **no reference was re-checked in this session**. Every literature statement below comes from the two colleague reports, with their tags ([read-abs], [snip], [memory]), or from my own memory, tagged [memory].
- Claims tagged [memory] must be checked before any citation. That applies in particular to the exact statements of Ge–Zhao, Macera–Sodin, Jitomirskaya–Schulz-Baldes–Stolz, L. Li and ADHH.

**Status labels used below.**
- **THEOREM (proved here):** complete proof given in this note.
- **THEOREM (conditional):** complete proof, given an explicitly stated hypothesis.
- **PROPOSITION / COROLLARY:** as above, smaller statements.
- **SKETCH:** gap identified explicitly.
- **CONJECTURE:** supported by numerics or heuristics only.

**Numerics.** All scripts are in `/tmp/claude-0/-home-user-fermions/e9e4eeb3-4dfa-5769-bb28-5161ded2b166/scratchpad/code/`, with a `.out` file next to each script (list in §8).

**Summary of results.**

| # | Statement | Status |
|---|---|---|
| T1 | Exact unitary duality: 1D Z2 LGT (SKKM type, two-link electric term) = free fermions coupled to static Ising charges. Explicit unitary; locality across a cut (one controlled-Z2 gate). | Proved. ED-verified to 1e-13 (spectra, N ≤ 6, three boundary conditions) and 1e-14 (dynamics, including string and Loschmidt observables) |
| C1 | z-polarized gauge state ⇒ exact i.i.d. Bernoulli(1/2) sector average for every gauge-invariant observable | Proved |
| T2a,b1 | Sector-averaged exponential dynamical localization; zero-velocity Lieb–Robinson bound for dressed fermions; bounded transport across any cut | Conditional on (DL_N). Unconditional in the thermodynamic limit, modulo Ge–Zhao [memory/secondary] |
| T3 | Strictly positive infinite-time CDW imbalance for **every** h ≠ 0, uniformly in N. New identity E_α a_α = 2h⟨φ_α, VΓφ_α⟩ | Unconditional inequality; uniformity in N conditional on (DL_N) |
| T4 | Uniform-in-time area law for the sector-averaged (conditional) entanglement, and for the full matter+gauge entanglement of definite-charge initial states | Conditional on (DL_N) |
| P5 | Structure of the three entropies (full, conditional, matter mixture); the duality shifts entanglement by ≤ ln 2 | Proved |
| CJ | S_full of the superposed state grows (∼ ln t); S_mix is volume law | Conjecture, with numerics and a sketch |
| T6 | Uniform-in-time bound on the number entropy of the full LGT state | Conditional on (DL_N) |
| P7 | Correlated labels: explicit product gauge state ⇒ random-dimer sector law. For h < J/2, DFL fails (superdiffusion) | Proved reduction (ED-verified). Delocalization via JSBS [memory]; numerics |
| T8 | 2D: exact mapping to free fermions with binary potential and binary π-flux | Proved; ED-verified (2×3 lattice) |
| P9 | 2D strong coupling: cluster confinement up to times ∼ h/J² | Proved |
| P10 | U(1) quantum link models: in a fixed charge sector, S_N ≤ ln(4S+1) (bulk interval), ≤ ln(2S+1) (interval at a fixed-flux edge), for any state | Proved; checked numerically |
| P11 | Wilson links: cut-averaged number entropy bounded uniformly in time by energy conservation | Proved |
| T12 | Prethermal conservation of electric energy or particle number (strong coupling or large mass) | Corollary of ADHH as recalled [memory] |
| T13 | Gauge protection by a quadratic penalty | Corollary of ADHH as recalled [memory] |
| — | One matter fermion in the *standard* (single-link) Z2 LGT = particle in a random-walk potential | Exact (ED-verified). Localization proof open |

---

## 1. L1: the 1D Z2 LGT with free matter

### 1.1 Model and boundary conditions

**Model M_N.**
- Sites j = 1..N carry spinless fermions f_j.
- Links b_j = (j, j+1), j = 1..N, carry spins 1/2. Link b_N is a **boundary link**: no hopping crosses it.
- Hilbert space: H_N = F_N ⊗ (C²)^{⊗N}.
- Hamiltonian:

      H = −J Σ_{j=1}^{N−1} σ^z_{b_j} (f†_j f_{j+1} + h.c.) − h Σ_{j=1}^{N} σ^x_{b_{j−1}} σ^x_{b_j},     with σ^x_{b_0} := 1.

- The two-link electric term is the form that makes the model free, as the colleagues' reconstruction indicates and as the ED below confirms. I could not compare it with SKKM's equations (arXiv:1701.04748 was unreachable).
- Gauss operators: G_j = σ^x_{b_{j−1}} (−1)^{n_j} σ^x_{b_j}. They commute with each other and with H, and G_j² = 1.
- With the boundary link, all 2^N joint eigenspaces P_q (q ∈ {±1}^N) are nonzero, and σ^x_{b_N} = (−1)^{N_f} Π_j G_j is conserved.

**Variants (all ED-verified, §1.7).**
- **"open":** no boundary link. Then Π_j G_j = (−1)^{N_f}, so q is constrained to a parity hyperplane.
- **"ring":** N links. Again Π G_j = (−1)^{N_f}, and in addition the Wilson loop W = Π_b σ^z_b is conserved.

### 1.2 THEOREM 1 (exact duality; proved)

**Dual space.** Let K_N = F_N ⊗ (C²)^{⊗N}, with c-fermions and τ-spins. Its basis is |n⟩ ⊗ |q⟩, where:
- τ^x_j |q⟩ = q_j |q⟩;
- |n⟩ = (c†_1)^{n_1} ⋯ (c†_N)^{n_N} |0⟩, and the same convention is used for the f-modes.

**The unitary.** Let |x⟩ be the σ^x eigenbasis, with |±⟩ = (|↑⟩ ± |↓⟩)/√2. Define U : K_N → H_N by

    U (|n⟩⊗|q⟩) = |n⟩ ⊗ |x(n,q)⟩,     x_{b_j}(n,q) := Π_{k≤j} q_k (−1)^{n_k}.

**Claims.** U is unitary, and:
1. U† G_j U = τ^x_j.
2. U† H U = Σ_q (dΓ(h_q) − h Σ_j q_j) ⊗ |q⟩⟨q|, where

       h_q = −J Δ_N + 2h diag(q)     on ℓ²({1..N}),

   and Δ_N is the open-chain adjacency matrix.
3. Dictionary:
   - U† f_j U = c_j τ^z_j;
   - U† σ^z_{b_j} U = τ^z_j τ^z_{j+1}, with τ^z_{N+1} := 1;
   - U† σ^x_{b_j} U = Π_{k≤j} τ^x_k (−1)^{n_k};
   - U† f†_j σ^z_{b_j} f_{j+1} U = c†_j c_{j+1};
   - U† n_j U = n_j.

   Equivalently, the **dressed fermions** c̃_j := U c_j U† = f_j Π_{k≥j} σ^z_{b_k} are gauge-invariant CAR fermions. Their Wilson line runs to the boundary link.
4. Locality across a cut: for 1 ≤ ℓ < N,

       U = C_ℓ (U_{≤ℓ} ⊗ U_{>ℓ}),     C_ℓ = P^+_{b_ℓ} ⊗ 1 + P^−_{b_ℓ} ⊗ Π_{k>ℓ} σ^z_{b_k},     P^± = (1 ± σ^x_{b_ℓ})/2.

   Here U_{≤ℓ} acts on (c_k, τ_k)_{k≤ℓ} → (f_k, σ_{b_k})_{k≤ℓ}, and U_{>ℓ} is defined in the same way on the right part, with x_{b_ℓ} := 1.

**Proof.**
- (n,q) ↦ (n, x(n,q)) is a bijection, because q_j = x_{b_{j−1}} x_{b_j} (−1)^{n_j}. So U maps one orthonormal basis onto another.
- Electric term and Gauss operators: σ^x_{b_{j−1}} σ^x_{b_j} acts on |x(n,q)⟩ as q_j (−1)^{n_j}. Hence U† G_j U = τ^x_j and U†(electric term)U = −h Σ q_j (1 − 2n_j).
- Hopping: let n_{j+1} = 1 and n_j = 0, and let n' = n − e_{j+1} + e_j. Then f†_j f_{j+1} |n⟩ = ε |n'⟩, with the same sign ε as for c†_j c_{j+1}. Moreover x(n',q) differs from x(n,q) only on b_j: for k > j, both parities have flipped. Since σ^z|±⟩ = |∓⟩, we get σ^z_{b_j} |x(n,q)⟩ = |x(n',q)⟩. This gives (2).
- (3): f_j changes n_j at fixed x. The compensating label change is q_j → −q_j, which is τ^z_j. The other entries follow the same way.
- (4): x_{b_k} for k ≤ ℓ depends only on (n,q)_{≤ℓ}. For k > ℓ, x_{b_k} = x_{b_ℓ} · y_k((n,q)_{>ℓ}). ∎

**Remark (correction to a claim in the colleague reports).** Take the "open" chain without a boundary link and left-anchored strings c_j = (Π_{k<j} σ^z_{b_k}) f_j. These c_j satisfy the CAR but **anticommute with G_1** (and so does every c_j); only their bilinears are gauge invariant. Right-anchored strings ending on the boundary link remove this defect.

**COROLLARY 1 (sector averaging; proved).** Let Ψ0 = |↑…↑⟩_z ⊗ ψ with ψ ∈ F_N arbitrary.
- Then U†Ψ0 = 2^{−N/2} Σ_q ψ ⊗ |q⟩. The charges are **exactly i.i.d. Bernoulli(1/2)**, and every sector carries the same matter state ψ.
- Hence for every gauge-invariant O, i.e. [O, G_j] = 0 for all j, so that U†OU = Σ_q O_q ⊗ |q⟩⟨q|:

      ⟨Ψ(t)|O|Ψ(t)⟩ = E_q ⟨ψ| e^{it dΓ(h_q)} O_q e^{−it dΓ(h_q)} |ψ⟩.

  The constants −hΣq drop out.
- Examples of O_q:
  - n_x and 1[N_A = k] are unchanged;
  - f†_x (Π_{x≤k<y} σ^z) f_y ↦ c†_x c_y;
  - σ^x_{b_j} ↦ (Π_{k≤j} q_k) (−1)^{N_{≤j}}, a Gaussian parity string.
- For |s⟩_z with some spins down: |↓⟩ = σ^x|↑⟩, so the sector matter state becomes Π_{j∈D} (−1)^{N_{≤j}} ψ, and the signs (Π_{k≤j} q_k) multiply the amplitudes. For a Fock state ψ this changes only signs.
- **Not sector-diagonal:** σ^z_{b_j} ↦ τ^z_j τ^z_{j+1}. Its expectation is a Loschmidt overlap, 2^{−N} Σ_q e^{iht(Σq − Σq')} ⟨ψ_{q'}(t)|ψ_q(t)⟩, where q' is q flipped at j and j+1. This was verified by ED, including the phase.

**Other boundary conditions (proved, ED-verified).**
- **"open":** the same statement holds, with q uniform on {Π q = (−1)^p} when ψ has fermion parity p. Its density relative to the product measure is ≤ 2.
- **"ring":** the sectors are (q, W). The one-body operator is the ring hopping with factor **W** on bond (N,1); there is no parity-dependent twist (checked sector by sector; the JW-type alternative fails with error 2.0). For |↑…↑⟩_z one has W = +1, so the boundary conditions are periodic.

### 1.3 The Anderson–Bernoulli input

**Hypothesis (DL_N).** For λ := 2h/J ≠ 0 there are C, η > 0 such that, for all N and all x, y ∈ {1..N},

    E_q sup_{t∈R} |⟨δ_x, e^{−ith_q} δ_y⟩| ≤ C e^{−η|x−y|},     q i.i.d. Bernoulli(1/2).

The stronger **eigenfunction-correlator form (EC_N)** is E Σ_α |φ_α(x) φ_α(y)| ≤ C e^{−η|x−y|}. It implies (DL_N), because for simple spectrum Σ_α |φ_α(x)φ_α(y)| = sup_{|f|≤1} |f(h)_{xy}| ≥ sup_t |U_t(x,y)|.

**What the literature gives** (status per the reports, not re-verified):
- **Carmona–Klein–Martinelli 1987** [snip]: almost-sure spectral localization on Z, every λ ≠ 0, all energies. Not dynamical, and not uniform in finite volume.
- **Bucaj et al. 2019; Jitomirskaya–Zhu 2019** [snip/memory]: spectral and dynamical localization on Z, via positivity plus large deviations. The exact dynamical form is unverified.
- **Ge–Zhao, J. Spectral Theory 10 (2020)** [secondary/memory]: exactly the form (DL) above, but **in infinite volume on Z**, for all bounded non-degenerate laws (Bernoulli included), all energies, all λ ≠ 0. This is the key reference for T2–T6 in the thermodynamic limit.
- **Macera–Sodin, AHP 23 (2022)** [snip]: eigenfunction-correlator bounds, with finite-interval and half-line operators and *arbitrary boundary conditions*. This is the most likely source of (EC_N) **uniform in N**. Unverified.
- **Energies and coupling.**
  - For i.i.d. Bernoulli there are **no critical energies**. A(E,v) A(E,v')^{−1} = [[1, v'−v],[0,1]] is parabolic, and A(E,v) does not fix e_1, so the group is non-compact and strongly irreducible (Zariski dense) for every E.
  - So **all h/J ≠ 0 and all energies** are covered.
  - The constants degenerate as λ → 0 (γ ∼ λ²). Macera 2024 [snip] shows that the annealed rate η is in general strictly below the Lyapunov exponent.
- **Periodic boundary conditions.** Not covered by the line theorems as stated, since the ring is a rank-2 perturbation of the open chain and dynamical localization is not stable under finite-rank perturbations in general. Folding the ring into a 2-channel strip with coupled ends may bring it under Macera–Sodin's quasi-1D framework (unverified).

**Numerical support for (DL_N)/(EC_N)** (`dfl_bounds.out`, 400 samples):
- E Q(x, x+d) is independent of N ∈ {100, 200, 400}, and nearly identical at the boundary site and at the centre.
- Fitted rates η ≈ 0.03, 0.085, 0.19, 0.25 for λ = 0.5, 1, 2, 4.

**How the results are stated.**
- In finite volume, uniformly in N, everything is conditional on (DL_N).
- In the **thermodynamic limit** at fixed t, finite-volume propagators converge entrywise, uniformly in q (Combes–Thomas plus finite speed). The sup_t bounds of T2(a), T2(b1), T4 and T6 then hold for the limit dynamics using only Ge–Zhao's infinite-volume bound. T3 is formulated in finite volume only.

### 1.4 THEOREM 2: localization and bounded transport (conditional on (DL_N))

Write U_q(t) = e^{−ith_q}. Let Ψ0 = |↑…↑⟩_z ⊗ ψ, and let E denote the i.i.d. sector average.

**(a) Propagator.** E sup_t |U_q(t)_{xy}| ≤ C e^{−η|x−y|}.

LGT content:
- The operator-valued propagator {c̃_x(t), c̃_y†} = Σ_q U_q(t)_{xy} P_q has sector-averaged modulus ≤ C e^{−η|x−y|}.
- For one matter fermion created at y, sup_t ⟨n_x(t)⟩ ≤ C e^{−η|x−y|}.
- **Zero-velocity Lieb–Robinson bound.** Let A be any even polynomial in c̃_y, y ∈ Y. These are local gauge-invariant matter observables, e.g. f†_x W_{xy} f_y; the strings cancel outside conv(Y). Then

      E sup_t ‖P_q [c̃_x(t), A] P_q‖ ≤ 2‖A‖ C Σ_{y∈Y} e^{−η|x−y|},

  because c̃_x(t) = Σ_z U_q(t)_{xz} c̃_z in sector q.

**(b1) Bounded transport.** Let ψ be a Fock state with occupied set S. For every cut m, every N and every t:

    E sup_t |⟨N_{≤m}(t)⟩_q − |S ∩ [1,m]|| ≤ E sup_t Σ_{x≤m<z or z≤m<x} |U_q(t)_{xz}|² ≤ 2C Σ_d d e^{−ηd}.

- The first inequality follows from ⟨n_x(t)⟩_q = Σ_{z∈S} |U_{xz}|² and Σ_x |U_{xz}|² = 1.
- Single-particle moments: E sup_t Σ_x |x−y|^p |U_q(t)_{xy}|² ≤ 2C Σ_d d^p e^{−ηd}.

*Proof.* Insert Corollary 1. Use |U|² ≤ |U| and sup Σ ≤ Σ sup. ∎

### 1.5 THEOREM 3: positive infinite-time imbalance for every h ≠ 0

**Setting.** N is even. The matter starts in the CDW (even sites occupied), and Γ := diag((−1)^x). Define

    I(t) := (2/N) Σ_x (−1)^x ⟨Ψ(t)|n_x|Ψ(t)⟩.

h_q is a Jacobi matrix with nonzero off-diagonal entries, so its spectrum is simple. Let (E_α, φ_α) be its real eigenpairs and a_α := ⟨φ_α, Γ φ_α⟩.

1. **(exact)** Ī := lim_{T→∞} T^{−1} ∫_0^T I(t) dt exists, and Ī = E_q [N^{−1} Σ_α a_α²].
2. **(identity)** E_α a_α = 2h ⟨φ_α, V_q Γ φ_α⟩, with V_q = diag(q). In particular Ī = 0 if h = 0, which is the chiral obstruction.
3. **(unconditional bound)** Let κ_N := sup over intervals B of E Σ_{x∈B, y∉B} R_q(x,y), where R_q(x,y) := Σ_α φ_α(x)² φ_α(y)² is the time-averaged transition probability. Then for every 2 ≤ ℓ ≤ N:

       Ī ≥ (⌊N/ℓ⌋ / N) · h² (√(ℓ/2) − 2κ_N)_+² / (ℓ (J+h)²),     and     Ī ≥ h² / (N (J+h)²).

4. **(uniform)** Under (DL_N), κ_N ≤ κ := 2C e^{−η} / (1 − e^{−η})². Hence

       inf_N Ī ≥ h² / (16 ℓ_* (J+h)²),     ℓ_* = max(2, ⌈32 κ²⌉).

**Proof.**
- **(1).** For Tr Γ = 0 we have I(t) = N^{−1} Tr(Γ U_q Γ U_q*). With simple spectrum its time average is N^{−1} Σ_α a_α². Then average over q, which is a finite sum.
- **(2).** Δ connects opposite sublattices, so ΓΔ = −ΔΓ. Hence hΓ + Γh = 4h V Γ. Take ⟨φ_α, · φ_α⟩ on both sides: the left side gives 2E_α a_α.
- **(3), block estimate.** Fix an interval B with |B| = ℓ, and let m_α(B) = ‖1_B φ_α‖². Write r_x = q_x (−1)^x, which is again i.i.d. symmetric ±1. By (2),

      X_B := Σ_α E_α a_α m_α(B) = 2h Σ_x r_x R_q(x,B).

  Cauchy–Schwarz with weights m_α, together with Σ_α E_α² m_α(B) = Tr(1_B h² 1_B) ≤ ℓ (2J+2h)², gives

      Σ_α a_α² m_α(B) ≥ X_B² / (ℓ (2J+2h)²).

- **(3), lower bound on |X_B|.** Pointwise, |X_B|/2h ≥ |Σ_{x∈B} r_x| − Σ_x |R(x,B) − 1_B(x)|. Moreover Σ_x |R(x,B) − 1_B(x)| = 2 Σ_{x∈B, y∉B} R(x,y), using Σ_y R(x,y) = 1 and the symmetry of R. Szarek's sharp Khintchine bound gives E|Σ_{x∈B} r_x| ≥ √(ℓ/2); Hölder gives √(ℓ/3) elementarily. Hence E X_B² ≥ (E|X_B|)² ≥ 4h² (√(ℓ/2) − 2κ_N)_+².
- **(3), summing.** Sum over ⌊N/ℓ⌋ disjoint blocks, using Σ_B m_α(B) ≤ 1.
- **(3), second bound.** Apply the same Cauchy–Schwarz globally: Σ a_α² ≥ 4h² (Σ_x r_x)² / Tr h², and E(Σ r_x)² = N.
- **(4).** E R(x,y) ≤ E sup_t |U_{xy}| ≤ C e^{−η|x−y|}, and an interval has at most 2d crossing pairs at distance d. Choose ℓ ≥ 32 κ², so that (√(ℓ/2) − 2κ)² ≥ ℓ/8. Use ⌊N/ℓ⌋ ≥ N/(2ℓ) for N ≥ 2ℓ; for N < 2ℓ use the global bound. ∎

**Remarks.**
- The proof uses only a time-averaged, eigenfunction-level localization input (E R(x,y) decays). It does not need sup_t.
- For the constrained "open" measure, blocks with ℓ < N are i.i.d. and E r_x r_y = 0 for N ≥ 3, so the same proof works.
- **Numerics** (`dfl_bounds.out`, `imbalance_bound.out`):
  - actual E Ī ≈ 0.015, 0.073, 0.166, 0.235 for λ = 0.5, 1, 2, 4, saturated in N;
  - the identity holds to 1.5e−14;
  - measured κ ≈ 16, 5.0, 2.5 for λ = 0.5, 1, 2, so the rigorous lower bound is positive but small (3e−7, 9e−6, 8e−5).
  - At λ = 0.2, N = 3200 is still not ≫ ξ.
- A uniform-in-time version, inf_t E I(t) > 0, would need explicit constants. It is open except heuristically at strong disorder.

### 1.6 Entanglement: which entropy (T4, P5, conjectures)

**Regions.**
- Region A = sites 1..ℓ plus links b_1..b_ℓ. Its dual is A' = (c_k, τ_k)_{k≤ℓ}.
- Fermionic reduced states are taken in Jordan–Wigner order, which is the canonical tensor factorization for a left segment.
- For an interior interval, every constant below doubles, and the duality defect becomes 2 ln 2 (two controlled gates).

**Definitions.**
- **S_G(t)** := E_q S(ρ^q_A(t)), where ρ^q_A is the reduced state on the c-modes of A of ψ_q(t) = Γ(U_q(t))ψ. This is the "free-fermion numerics" quantity: the entanglement **conditional on a measurement of all Gauss charges** (a gauge-invariant measurement).
- **S_full(t)** := S(Tr_{A^c} |Ψ(t)⟩⟨Ψ(t)|), the full matter + gauge entanglement in the original variables.
- **S_mix(t)** := S(E_q ρ^q_A(t)), the dressed-matter state with all charges traced out. S_mix − S_G = χ(t) ≥ 0 is the Holevo information that the matter in A holds about the charges, i.e. "classical mixing entropy".
- **S_f(t)** := entropy of the bare matter f_A with all links traced out (original variables). In a single sector it is the Shannon entropy of the occupation distribution of A, because the electric strings record n.
- **S_N(t)**: number entropy (§1.8).

**THEOREM 4 (area laws; conditional on (DL_N)).** Let ψ be any Fock state, and define

    M_η := Σ_{d≥1} d (d + 1 + 2/(e^η − 1)) e^{−ηd/2}.

1. E sup_t S(ρ^q_A(t)) ≤ 2 ln 2 · C M_η, uniformly in N, ℓ, t. In particular sup_t S_G(t) ≤ 2 ln 2 · C M_η.
2. **Definite-charge (quenched) initial states.** Let Ψ0 = |x⟩_σ ⊗ |n⟩, with σ^x-eigenstate links and Fock matter, so that q = q(n,x) is definite. Then S_full(t) = S(ρ^{q}_A(t)) exactly. For x i.i.d. uniform (equivalently q i.i.d.), E_x sup_t S_full(t) ≤ 2 ln 2 · C M_η.

**Proof of (1).**
- For a pure quasi-free state with correlation projection C, S(ρ_A) = Σ_i h(ν_i) over the eigenvalues ν_i of P_A C P_A.
- **Lemma:** h(ν) ≤ 2 ln 2 √(ν(1−ν)).
  - Put ν = (1−s)/2. Then h = ln 2 − Σ_k s^{2k}/(2k(2k−1)) and 2 ln 2 √(ν(1−ν)) = ln 2 − ln 2 Σ_k c_k s^{2k}, with c_k = |binom(1/2, k)|.
  - Let d_k := 1/(2k(2k−1)) − c_k ln 2. Then d_1 > 0, and d_k < 0 for k ≥ 2: the ratio e_k/c_k is decreasing and equals 2/3 < ln 2 at k = 2.
  - Σ_k d_k = 0, so Σ_k d_k s^{2k} ≥ d_1 (s² − s⁴) ≥ 0.
- Hence S ≤ 2 ln 2 Tr √(C_A(1−C_A)) = 2 ln 2 ‖P_A C P_B‖_1 ≤ 2 ln 2 Σ_{x∈A, y∉A} |C_xy|.
- For C(t) = U P_S U*: |C_xy| ≤ Σ_z min(|U_xz|, |U_yz|). Then

      E sup_t |C_xy| ≤ C Σ_z e^{−η max(|x−z|, |y−z|)} ≤ C (d + 1 + 2/(e^η − 1)) e^{−ηd/2},     d = |x−y|.

- A single cut has d pairs at distance d.

**Proof of (2).**
- U†Ψ(t) = ψ_q(t) ⊗ |q⟩, so the τ part is a product state.
- With fixed fermion parity π, the control bit of C_ℓ is x_{b_ℓ} = π (Π_{k≤ℓ} q_k) (−1)^{N_{>ℓ}}, an operator on the B side only. So C_ℓ is a local unitary on B, and the entropy is unchanged.
- ED check: E|S_full − S_G^{(q)}| = 0.000 at every t (`ent_small.out`). ∎

**PROPOSITION 5 (structure; proved).** For the superposed state Ψ0 = |↑…↑⟩_z ⊗ ψ with ψ Fock:
1. |S_full(t) − S_{A'}(Φ(t))| ≤ ln 2, where Φ(t) = U†Ψ(t) = 2^{−N/2} Σ_q ψ_q(t) ⊗ |q⟩.
   *Proof.* C_ℓ = P^+ ⊗ 1 + P^− ⊗ V leaves the P^± diagonal blocks of ρ_A unchanged. Pinching with 2 blocks raises the entropy by at most ln 2, applied in both directions.
2. S_G ≤ S_mix ≤ min(S_G + H(q), ℓ ln 2), by concavity, Holevo, and dimension.
3. N_A commutes with ρ_A, ρ^q_A and ρ_mix (fixed N_f). So S_full, S(ρ^q_A), S_mix ≥ S_N for the respective distributions. The distribution of N_A is the same in every picture.
4. The dynamics of every gauge-invariant observable is identical for the superposed state and for the ensemble of definite-charge states (Corollary 1). **Only entanglement-type quantities distinguish them.** This is the "superposition vs sampling" distinction of the Google experiment, made exact.

**Numerics** (`ent_small.out`; exact, N = 10, h = J = 1, ℓ = 5; maxima: S_full ≤ 10 ln 2 = 6.93, S_mix ≤ 5 ln 2 = 3.47):

| t | S_full | S_G | S_mix | S_f | S_N |
|---|---|---|---|---|---|
| 1 | 1.45 | 0.61 | 2.65 | 2.57 | 0.67 |
| 10 | 3.40 | 0.91 | 2.57 | 2.67 | 0.82 |
| 100 | 5.00 | 0.98 | 2.78 | 2.83 | 0.89 |
| 10^5 | 5.38 | 0.95 | 2.90 | 2.93 | 0.88 |

- At h = 2, S_full = 2.6, 3.0, 4.0, 4.5, 4.7, 5.0 at t = 10, 30, 100, 300, 10³, 10⁴. This is slow, log-like growth before finite-size saturation.
- At large N (`dfl_bounds.out`, up to N = 256, t ≤ 10⁸), S_G is bounded in t and N, e.g. ≈ 1.25 at λ = 2 and 3.2 at λ = 1.

**CONJECTURE A (log growth of S_full; SKETCH).** Under (DL_N): S_full(t) ≤ C_1 ln(2+t) + C_2, uniformly in N. Moreover S_full is unbounded in t as N → ∞.

*Rigorous ingredient.* Duhamel gives E sup_{s≤t} ‖(U_q(s) − U_{q'}(s)) δ_z‖ ≤ 4h t C Σ_{x∈D} e^{−η|x−z|}, where q, q' differ on a set D. So charges at distance R from the cut can influence modes near the cut only after times t ≳ e^{ηR}/h. The structure is that of an l-bit model with exactly conserved τ's.

*Gap.* This must be turned into an entanglement bound that is uniform in N. Entropy continuity costs ε ln(4^N). One needs a telescoping decoupling argument with summable errors (as in ANSS-type proofs), adapted to a coherent superposition of 2^N sectors. I have not done this.

**CONJECTURE B (volume law of S_mix and S_f).** The time-averaged S_mix ≥ c ℓ. Numerics: S_mix/(ℓ ln 2) ≈ 0.84 (h = 1) and 0.77 (h = 2) at N = 10.

Heuristic: at strong h each resonant pair contributes an O(1) mixture between the frozen and the oscillating configuration, and these contributions are independent across pairs.

**Consequence for "area law" claims.** Only S_G (equivalently, the quenched S_full) satisfies a uniform-in-time area law. S_mix and S_f do not, and S_full is not expected to. A "QDL-like" pattern (S_G area law, S_full and S_mix large) is automatic in this model.

### 1.7 ED verification of the mapping (`lgt1d_ed.out`, `ring_sectors.out`, `dyn_check.out`)

**Spectra.** Full LGT spectrum vs the union of mapped free-fermion sectors:
- "open with boundary link", N = 2..6 (dimensions up to 4096); "open", N = 2..6 (up to 2048); "ring", N = 3..6 (up to 4096): max |ΔE| ≤ 1.1e−13 in all cases (two independent runs).
- The union-of-spectra test cannot tell the two ring twists apart, because the union over W = ±1 is the same multiset. Hence the sector-resolved check below.

**Ring, sector by sector** (projecting onto joint eigenspaces of G_j and W, N = 3, 4, 5):
- twist = W: error ≤ 9e−15;
- twist = W(−1)^{N_f+1}: error 2.0, so this alternative is rejected.

**Dynamics** (N = 4, 5, 6; CDW and other Fock states; h ∈ {0.4, 0.7, 1.3}; t ≤ 10). Max deviation ≤ 9e−15 for:
- ⟨n_x⟩;
- the dressed current;
- ⟨σ^x_b⟩ (string observable, determinant formula);
- ⟨σ^z_b⟩ (Loschmidt overlap). This one required the sector-dependent phase e^{iht(Σq − Σq')}; omitting it gives O(1) errors.

**Entanglement.** Reconstructing Ψ = UΦ reproduces the ED values of S_full and S_f to 1e−15 (`validate_ent.py`, N = 6).

### 1.8 THEOREM 6: number entropy (conditional on (DL_N))

**Setting.** N_A = Σ_{x∈A} n_x is gauge invariant and conserved, and Ψ has fixed N_f. So ρ_A (matter + links) is block diagonal, S_full = S_N + Σ_k p_k S(ρ^{(k)}), and

    S_N(t) = H(p_t),     p_t(k) = E_q p^q_t(k).

**Statement.** For Fock ψ and A a half chain, uniformly in t and N,

    S_N(t) ≤ ½ ln(2πe (K_1 + K_2 + 1/12)),     K_1 = C M_η,     K_2 = (2√C Σ_{d≥1} d e^{−ηd/2})².

Also sup_t E_q S^q_N(t) ≤ ½ ln(2πe (K_1 + 1/12)).

**Proof.**
- Massey: H(X) ≤ ½ ln(2πe (Var X + 1/12)) for integer X. Proof: H(X) = h(X + U) for U uniform on [0,1), and the Gaussian maximizes differential entropy at fixed variance.
- Variance decomposition: Var_mix N_A = E_q Var_q + Var_q(⟨N_A⟩_q).
- First term: Var_q(N_A) = ‖P_A C P_B‖_2² = Σ_{x∈A, y∉A} |C_xy|² ≤ Σ |C_xy|, so E sup_t Var_q ≤ K_1.
- Second term: |⟨N_A⟩_q(t) − |S ∩ A|| ≤ Σ_{crossing pairs} s_{xz}, with s_{xz} = sup_t |U_xz| ∈ [0,1]. Since E s_i s_j ≤ √(E s_i E s_j), E(Σ s)² ≤ (Σ √(E s))² ≤ K_2. Use Var ≤ E(· − const)².
- The per-sector bound follows from Jensen. ∎

**Numerics** (`dfl_bounds.out`). At λ = 2, N = 256, t up to 10⁸: Var_mix ≈ 0.43, S_N(mixture) ≈ 0.99, Massey bound 1.09. The values are flat in t and identical for N = 64, 128, 256.

### 1.9 Subtleties

**Correlated vs i.i.d. labels.**
- **Bounded density.** If the sector law has density ≤ K relative to the product law, every "E sup" bound is multiplied by K (e.g. K = 2 for "open").
- **Biased σ^x-product gauge states** (amplitudes cos, sin) with Fock matter give q_j = (−1)^{n_j} x_{j−1} x_j with x i.i.d. biased. This is a 1-block factor of an i.i.d. sequence, not i.i.d.
  - Damanik–Gorodetski–Kleptsyn 2025 [snip] give spectral localization, and dynamical localization off a finite exceptional set.
  - The support is the full shift, so I expect no critical energies, but this is unverified.
- **PROPOSITION 7 (DFL can fail; proved reduction).**
  - Gauge state: odd links in |↑⟩_z; even links b_{2k} in the σ^x eigenstate x = (−1)^k. Matter: CDW on even sites.
  - This produces (q_{2k−1}, q_{2k}) = (ε_k, ε_k) with ε_k i.i.d. uniform, and the same matter state in every sector. ED, N = 6: exactly the 8 dimer sectors, each with weight 1/8.
  - So the sector-averaged dynamics is the **random dimer model** with energies ±2h.
  - For 0 < h < J/2 it has critical energies E_c = ±2h with γ(E_c) = 0. At E = ε_a the a-dimer transfer matrix is −1, and the b-dimer matrix is elliptic because |ε_a − ε_b| = 4h < 2J.
  - By Jitomirskaya–Schulz-Baldes–Stolz (CMP 233, 2003) [memory], time-averaged moments grow like T^{q−1/2−α}, i.e. superdiffusively. So **(DL) fails, and transport is unbounded**.
  - Numerics (N = 3000, single particle):
    - dimer law, h = 0.25: ⟨x²⟩ = 127, 5.3e3, 1.2e5 at t = 10, 10², 10³, i.e. ∼ t^{1.35–1.5};
    - i.i.d. law at the same h: ⟨x²⟩ saturates at ≈ 900;
    - dimer law, h = 0.75: ⟨x²⟩ saturates at ≈ 48.
  - **DFL is a property of the initial-state-induced sector law, not of the Hamiltonian.**

**Gauss-law projection.** In the strictly gauge-invariant subspace (G_j ≡ 1) there is a single sector and h_q = −JΔ + 2h, which is ballistic. DFL needs initial states that superpose sectors, or equivalently random static background charges.

**Boundary conditions.**
- The boundary link gives exact i.i.d. labels.
- "open" gives density 2.
- The ring needs (DL) for ring operators (unverified).
- The identity in T3 needs a bipartite lattice (ring with N even).

**Convention.** The potential is 2h q_j ∈ {±2h}; the disorder parameter is λ = 2h/J. Other papers write ±h depending on the normalization of the electric term.

---

## 2. L2: is any of this stated rigorously?

- Both colleague reports (abstract-level searches) found **no paper stating a rigorous dynamical-localization theorem for any LGT**, and no explicit statement of the corollary "SKKM + Ge–Zhao ⇒ averaged DFL".
- I could not search further (budget exhausted, arXiv unreachable).
- The duality itself is SKKM's mapping, presumably stated at the physics level of rigor. The specific items here appear new to me, but **this is unverified**:
  - the explicit unitary with the cut-locality (iv);
  - the boundary-link variant giving exact i.i.d. labels;
  - the imbalance theorem for all h;
  - the entropy taxonomy with the ln 2 duality bound;
  - the dimer counterexample;
  - the Gauss-law and energy caps on number entropy.

  The closest rigorous background is ANSS 2016, Hamza–Sims–Stolz 2012 and Pastur–Slavin/Elgart–Pastur–Shcherbina, all for absolutely continuous disorder.

---

## 3. L3: two dimensions

**THEOREM 8 (exact 2D mapping; proved, ED-verified).**
- Hamiltonian:

      H = −J Σ_{⟨ij⟩} σ^z_{ij} (f†_i f_j + h.c.) − h Σ_i Π_{j∼i} σ^x_{ij}     on a planar box.

- Conserved quantities: G_i = (−1)^{n_i} Π_{j∼i} σ^x_{ij} and W_p = Π_{∂p} σ^z.
- In each joint sector (q, W), H ≅ dΓ(h_{q,W}) restricted to the parity sector Π q = (−1)^{N_f}, minus hΣq. Here h_{q,W} = −J Σ s_ij (|i⟩⟨j| + h.c.) + 2h diag(q), with link signs s realizing the plaquette fluxes W.
- Proof as in T1, after gauge-fixing the σ^z.
- ED, 2×3 lattice (dimension 8192): max |ΔE| = 8.4e−14 (`lgt2d_ed.out`).

**Sector laws from the initial gauge state.**
- **z-polarized links:** W ≡ 1 and q i.i.d. This is the **2D Anderson–Bernoulli** model.
- **σ^x-polarized links with checkerboard CDW:** q = (−1)^{n_i} is deterministic, and the W_p are i.i.d. uniform. Then h = hop_W − 2hΓ with {hop_W, Γ} = 0, so **h² = hop_W² + 4h²**. The dynamics is a function of the binary random-π-flux hopping operator, and its chiral E = 0 (Gade) point is mapped to the band edges |E| = 2h. Proved identity.

**PROPOSITION 9 (strong-coupling cluster confinement up to t ∼ h/J²; proved).**
- Setup: on Z² (or any graph of degree 4), split h_q = h⁰ + h¹, where h¹ is the hopping across bonds with q_i ≠ q_j.
- For h > 2J, the spectra of h⁰ on "+" and "−" sites are separated by g = 4h − 8J.
- Solve [h⁰, S] = h¹ (Bhatia–Davis–McIntosh: ‖S‖ ≤ (π/2) 4J/g). Then e^S h e^{−S} = h⁰ + R with ‖R‖ ≤ C J²/g.
- Hence ‖e^{−ith_q} − e^{−ith⁰}‖ ≤ C(J/g + tJ²/g), uniformly in the volume.
- h⁰ is block diagonal over the maximal equal-charge clusters. So for t ≤ ε g/J², a particle started at x leaves its cluster C(x) with probability O((J/g + ε)²).
- At p = 1/2 < p_c^site(Z²) ≈ 0.5927, the cluster radius has exponential tails (Menshikov / Aizenman–Barsky [memory]).
- **Beyond t ∼ h/J² this argument cannot be pushed.** Second-order tunnelling connects equal-charge sites at distance 2, including diagonals. Equal-charge sites with NN+NNN adjacency percolate, since 1/2 > 0.4073 (the matching-lattice threshold), so the effective second-order network is a supercritical quantum-percolation problem with O(1) effective disorder.
- On the triangular lattice (p_c = 1/2) and on Z³ (p_c ≈ 0.31) even the first step fails.

**Available theorems and open problems.**
- **Ding–Smart 2020** [memory]: Anderson localization at the band edge of 2D Bernoulli. It is relevant only to the edge part of the spectral measure of an initial state; a Fock initial state has most of its weight elsewhere. A dynamical (in expectation) version from Bernoulli MSA is likely only polynomial (unverified).
- **L. Li** (arXiv:2002.11580) [read-abs/memory]: large-disorder localization outside small windows around finitely many energies (Dirichlet eigenvalues of finite clusters). If verified, this gives spectral localization of most of the spectrum at h ≫ J. The excluded windows carry O(1) spectral weight of a Fock state (e.g. monomer sites sit exactly at a cluster eigenvalue), so this is partial.
- **Binary random π-flux** (σ^x-polarized sector): **open**. Existing random-flux theorems need regular distributions [read-abs/memory].
- **Compactly supported eigenstates** on Z² do not exist for any potential or sign pattern. Proof: take the support's maximizer of x+y with the largest x; the neighbour to its right has a single neighbour in the support, and the eigenvalue equation fails there. Aharonov–Bohm caging needs other lattices, e.g. the rhombic or dice lattice, where π-flux sectors give flat bands. A sector-dependent strict localization theorem is then trivial.
- **Majorana analogue (Kitaev, class D):** a thermal metal is expected, i.e. delocalization [read-abs/N].

---

## 4. L4: 1D U(1) (Schwinger and quantum link) models

### 4.1 Gauss-law caps on number entropy

**Conventions.**
- Staggered fermions: ρ_j = n_j − (1 − (−1)^j)/2.
- Gauss law: E_j − E_{j−1} = ρ_j + Q_j, with integer background charges Q_j.
- E_0 is a fixed boundary flux.
- "Number" means the **conserved** fermion number N_A = Σ_{j∈A} n_j, equal to the U(1) charge up to a constant. The particle + antiparticle number is not conserved, and ρ_A is not block diagonal in it, so it admits no number entropy.

**PROPOSITION 10 (proved; checked in `qlm_gauss.out`).** For any state in a fixed Q-sector and any interval A = [a,b]:

    N_A = E_b − E_{a−1} + #{odd j ∈ A} − Σ_A Q_j.

*Proof:* sum Gauss's law over A.

Hence, for **spin-S quantum links** (|E| ≤ S):
- S_N(A) = H(E_b − E_{a−1}) ≤ **ln(4S+1)** for an interval in the bulk;
- S_N(A) ≤ **ln(2S+1)** if a = 1 (fixed-flux edge). This is the half-chain case.
- For a ring the half chain has two cuts, giving ln(4S+1).
- The task's ln((2S+1)²) is correct but not sharp.
- Superposition over Q-sectors: S_N ≤ ln(4S+1) + H(Σ_A Q). The second term is time independent.
- Numerics: random states in sectors with S = 1/2, 1, 3/2 give at most 2S+1 (edge) and 4S+1 (bulk) values. The bulk bound 5 = 4S+1 is attained for S = 1.

**PROPOSITION 11 (Wilson links; proved).**
- In a fixed sector with E_0 fixed, S_N([1,ℓ]) = H(E_ℓ) ≤ ½ ln(2πe (⟨E_ℓ²⟩ + 1/12)).
- Energy conservation, for H = −w Σ (f† U f + h.c.) + m Σ (−1)^j n_j + (g²/2) Σ E_j², gives (1/L) Σ_ℓ ⟨E_ℓ²⟩_t ≤ (2/g²)(⟨H⟩_0/L + 2w + |m|).
- Jensen then gives the **cut average** (1/L) Σ_ℓ S_N([1,ℓ])(t) ≤ ½ ln(2πe ((2/g²)(e_0 + 2w + |m|) + 1/12)), **uniformly in t and L**, with e_0 = ⟨H⟩_0/L. By Markov's inequality the same bound, up to constants, holds for all but a small fraction of cuts.
- An individual cut is only bounded by ½ ln O(L).

**Implications for diagnostics.**
- In QLMs, saturation of S_N is **kinematic**, whatever the dynamics. With gauge fields integrated out (e.g. Brenes et al.), S_N is the Shannon entropy of the electric field on the cut link.
- With Wilson links, unbounded growth of S_N on average over cuts is excluded by energy conservation, independently of localization.
- Only discrete gauge groups leave S_N essentially free: Z_n fixes N_A mod n. There, T6 is a genuine localization statement.
- I could not check what Brenes et al. or later papers actually computed under the name "number entropy".

### 4.2 THEOREM 12: prethermal emergent conservation (corollary of ADHH; ADHH form [memory])

**ADHH (CMP 354, 2017), as recalled.**
- Hypotheses: H = J N + Y, where N = Σ_x N_x with commuting, finite-range N_x and **integer spectrum**, and Y is local with finite local norm.
- Conclusion: a quasi-local unitary U, close to 1, with U*HU = JN + D + V, [D, N] = 0, and ‖V‖_loc ≤ C e^{−c n_*}, where n_* ∼ J/(‖Y‖ log³(J/‖Y‖)).
- Consequence: |⟨N(t)⟩ − ⟨N(0)⟩|/L ≤ C(‖Y‖/J + t e^{−c n_*}).

**(A) Strong coupling, U(1) QLM.**
- N = Σ_j (E_j + θ)² has integer spectrum up to a constant iff **2θ ∈ Z**. For spin-1/2 this is the number of links with E = +1/2.
- Y = hopping + mass is bounded. So the electric energy, i.e. string length, is conserved up to times exp(c g²/((w + |m|) log³)).
- D = (N-conserving part of Y) + O(‖Y‖²/g²). At θ = 1/2 this includes first-order resonant hoppings, the fragmented or PXP-type dynamics.

**(B) Large mass.**
- N = Σ_j (−1)^j n_j + const = particle + antiparticle number. Every hopping changes it by ±2, so D starts at order w²/m.
- Particle number is conserved up to exp(c m/((w + g²S²) log³)).
- In ADHH ‖Y‖ includes the electric term, so this needs m ≫ g²S² as well.

**Hypothesis check.**
- Integer spectrum and locality hold for truncated links.
- **Wilson (unbounded) links:** case (A) plausibly goes through, since only conjugation by e^{isN} (on-site, bounded) enters, but ADHH is formulated for bounded local terms, so this is unverified. Case (B) has an unbounded Y and is not covered.
- Incommensurate g² and m: only the dominant scale can be used.
- The result is prethermal (fragmentation), **not localization**.

### 4.3 Random-walk potentials

**Exact reduction (ED-verified, N = 3..6, error ≤ 2e−14).** One matter fermion in the **standard** Z2 LGT, with single-link term −h Σ σ^x_b and the boundary-link convention, in a uniform random sector, is equivalent to

    h = −JΔ + V,     V(x) = 2h Σ_{b≥x} Q_b,     Q_b = Π_{k≤b} q_k   i.i.d. ±1.

This is a random walk with steps ±2h, so |V| ∼ √|x|. Two or more fermions are interacting: in ED with N = 4 and 2 fermions the spectrum is not a free spectrum.

U(1) analogues:
- a tightly bound meson sees a random-walk potential;
- an isolated charge sees string tension plus an integrated walk, V ∼ |x|^{3/2}.

**What is provable.**
- Simon–Spencer (CMP 125, 1989) [memory]: lim sup |V| = ∞ on both sides implies σ_ac = ∅ almost surely.
- Pure point spectrum and dynamical localization are **open**:
  - the increments are independent but the values are not, so the Furstenberg / Le Page machinery is absent;
  - integer steps cause exact degeneracies between identical "wells" (returns of the walk), with no Wegner estimate.
- Suggested route for continuous increments: shifting one increment moves the potential rigidly on a half-line, which gives a one-variable Wegner estimate. Combine it with multiscale analysis, using |V − E| ≫ 1 off sparse returns, with local time ∼ √n. This seems within reach.
- For integer increments (the physical case) one would need a resonance-counting argument that exploits the super-exponential tunnelling suppression exp(−Σ ln |V|). This is hard.

**Numerics** (`rw_potential.out`): the eigenfunction correlator decays faster than in the i.i.d. Bernoulli model at large distance (h = 0.5, d = 64: 1.4e−4 vs 4e−3), consistent with super-exponential localization.

### 4.4 THEOREM 13: gauge protection (corollary of ADHH [memory])

**Setting.** H = H_0 + λH_1 + V Σ_j (G_j − g_j)², with:
- G_j integer-valued local generators (truncated U(1); for Z2 use (1 − G_j)/2);
- H_0 gauge invariant, H_1 arbitrary local.

**Statement.** Let N = Σ_j (G_j − g_j)². It is commuting and local, with integer spectrum, and N^{−1}(0) is exactly the target sector. ADHH, with J = V and Y = H_0 + λH_1, gives:
- for gauge-invariant initial states, the violation density obeys ⟨N(t)⟩/L ≤ C(‖Y‖/V + t e^{−c n_*});
- the rotated-frame dynamics on N = 0 is generated by the quasi-local gauge-invariant D|_{N=0} = P_0 (H_0 + λH_1) P_0 + O(‖Y‖²/V). This is the "renormalized gauge theory".

**Caveats.**
- V must dominate **all** local couplings, including H_0, not just λ.
- A **linear** penalty V Σ_j c_j G_j protects only Σ_j c_j G_j. Correlated violations that conserve it (e.g. G_i up, G_k down for c_i = c_k) are not suppressed exponentially. Exponential protection of each G_j additionally needs a combinatorial "compliance" property of (c_j) up to the relevant range, and unbounded c_j (Stark-type) needs an extension of ADHH.

---

## 5. L5: ranked feasible projects

| Rank | Project | Difficulty | Payoff |
|---|---|---|---|
| 1 | Locate or prove **(EC_N) uniform in N** for 1D Bernoulli, with boundary conditions, then publish T1–T6 + P7 as "Rigorous disorder-free localization in a Z2 LGT". Write the finite-volume version via Macera–Sodin, or redo the large-deviation proofs on intervals. Add the ring by folding. | Low–medium | High: first rigorous DFL theorem; also clarifies the entropies and the role of the initial state |
| 2 | **Conjecture A**: rigorous S_full(t) ≤ C ln t (exact l-bit structure with τ's; Duhamel sensitivity is already in hand) | Medium–high | High: first rigorous logarithmic entanglement-growth bound in a disorder-free system |
| 3 | Volume law for S_mix and the Holevo χ (Conjecture B), e.g. via Rényi-2 and approximate independence of blocks | Medium | Medium |
| 4 | Spin/Loschmidt observables ⟨σ^z_b(t)⟩: absence of orthogonality catastrophe for a two-site potential flip in the localized Bernoulli chain (Gebert–Küttler–Müller type) [memory] | Medium | Medium–high: rigorous memory of the gauge field |
| 5 | Write T12/T13 as stand-alone theorems with constants, including Wilson links in case (A) | Low–medium | Medium (mostly formalization) |
| 6 | Random-walk potential: continuous increments (Wegner by rigid shift + MSA), then integer increments | Medium / hard | Medium: the single-charge sector of the standard Z2 LGT and the U(1) meson |
| 7 | 2D: dynamical upgrade of L. Li at strong coupling beyond t ∼ h/J²; percolation-cluster MSA | Hard | High |
| 8 | Binary random π-flux localization; interacting DFL (= binary-disorder MBL); quantum percolation transition | Open / very hard | — |

---

## 6. Honest novelty assessment

- **T1/C1 (mapping).** Physically known (SKKM 2017). New here: the explicit unitary, the boundary link giving exact i.i.d. labels, the cut-locality via one controlled gate, the ring twist (= W, no parity twist) and the numerical verification. Modest novelty.
- **T2, T4(1), T6.** Elementary corollaries of (DL) plus standard quasi-free inequalities, ANSS-style. Probably unstated for LGTs [per reports]. Their value is in making DFL rigorous; there is **no new mathematics** beyond locating the finite-volume (DL_N).
- **T3 (imbalance for all h ≠ 0).** The identity E a = 2h⟨VΓ⟩ and the block Cauchy–Schwarz + Khintchine argument are, to my knowledge, new. I have not seen a proof of positive CDW imbalance at weak disorder; unverified. Small but genuine.
- **P5/T4(2) (entropy taxonomy, ln 2 duality bound, quenched equality).** The conceptual points are known in physics (QDL, Google "superposition vs sampling"). The precise statements are new as rigorous statements.
- **P7 (dimer counterexample).** Likely new as an LGT observation. The mathematics is JSBS's.
- **P9.** Elementary; probably not stated.
- **T8.** Physics-known; ED verification new.
- **P10/P11 (number-entropy caps).** Trivial consequences of Gauss's law and energy conservation. I could not check whether they are stated. They are relevant as a caveat for number-entropy diagnostics in LGTs.
- **T12/T13.** Direct ADHH applications, in the spirit of Halimeh et al.'s appendices [read-abs]. New only in formulation: the quadratic-penalty argument through N^{−1}(0), and the caveats on linear penalties.
- **Conjectures A/B and the random-walk problem.** Open.

---

## 7. What would change these conclusions

- If (EC_N) uniform in N turns out to be unavailable for Bernoulli, T2–T6 hold only in the thermodynamic limit, via Ge–Zhao, which itself still needs re-reading.
- If SKKM's electric term is the single-link one, their model is the interacting one of §4.3, and only the two-link model is covered by T1–T6.
- The JSBS statement is used in P7 from memory. The numerics support superdiffusion independently.

## 8. Files

All in `/tmp/claude-0/-home-user-fermions/e9e4eeb3-4dfa-5769-bb28-5161ded2b166/scratchpad/code/`:

| Script | Output | Content |
|---|---|---|
| `lgt1d_ed.py` | `lgt1d_ed.out` | Spectra, three boundary conditions; dynamics (with the Loschmidt-phase fix); standard model one-particle sector |
| `ring_sectors.py` | `ring_sectors.out` | Ring, sector by sector |
| (dynamics run) | `dyn_check.out` | Dynamics comparison |
| `ent_small.py`, `validate_ent.py` | `ent_small.out` | Entropies |
| `dfl_bounds.py` | `dfl_bounds.out` | Correlators, imbalance, number entropy, S_G |
| `imbalance_bound.py` | `imbalance_bound.out` | κ and the T3 lower bound |
| `dimer_check.py` | `dimer_check.out` | P7 |
| `lgt2d_ed.py` | `lgt2d_ed.out` | T8 |
| `qlm_gauss.py` | `qlm_gauss.out` | P10 |
| `rw_potential.py` | `rw_potential.out` | §4.3 |
