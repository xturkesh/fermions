# Mathematical toolkit for proving disorder-free localization after resolving superselection sectors: a map of the literature (status 2026-11-22)

## How to read this report
- **Tags.** [snip]: the bibliographic data and the content stated were confirmed in this session from web-search result snippets (abstract fragments). The paper itself was not opened. [memory]: from prior knowledge, not checked this session. [deriv]: a short argument written out here, which you should check.
- **No item is tagged [read].** Fetching pages failed: arxiv.org, ar5iv and link.springer.com gave DNS errors in this sandbox, and the Agent "web-fetch" subagent is not available here. The shared web-search budget ran out after about 25 searches, so most of items 2–5 rests on [memory].
- **arXiv ids.** An id is given only when it was confirmed [snip] or I am confident of it. Otherwise it reads "id?".
- **Claim types.** (P) proven theorem (assumptions stated), (N) numerical evidence, (H) heuristic, (E) experiment.

---

## 0. Executive summary
- **1D i.i.d. Bernoulli Anderson model: everything is in place in infinite volume.**
  - Spectral localization holds at all energies for every coupling λ≠0 (Carmona–Klein–Martinelli 1987).
  - Exponential dynamical localization in expectation, E sup_t |<x|e^{-itH}|y>| ≤ C e^{-c|x-y|}, holds at all energies (Ge–Zhao 2020). Macera–Sodin 2022 give the eigenfunction-correlator form, for quasi-1D models and random hopping.
  - **Gap:** none of the Bernoulli papers states the bound with constants uniform over finite boxes, which is what Abdul-Rahman–Nachtergaele–Sims–Stolz (ANSS) entanglement bounds need. Macera–Sodin's method, which allows arbitrary boundary conditions, is the natural source. Not verified.
- **Correlated sector measures** (Markov, block factor, random word):
  - Localization survives. Dynamical localization is proven only away from a finite set of critical energies.
  - At critical energies of dimer type, time-averaged moments provably grow at least like T^{q-1/2-α} (superdiffusive, Jitomirskaya–Schulz-Baldes–Stolz 2003). So no uniform dynamical-localization bound can hold there.
- **Random-hopping chain (your model).**
  - E=0 is a true critical point: the Lyapunov exponent vanishes there, γ(0)=0.
  - Rigorous dynamical localization exists only on energy windows that avoid 0.
  - I know of no rigorous result on transport, entanglement or number fluctuations at E=0.
  - Physics expectation: infinite randomness, polylogarithmic spreading, S(t) ~ ln ln t.
  - **Practical warning [deriv]:** in any bipartite chiral hopping model, the time-averaged imbalance of a density-wave state is exactly 0, localized or not (§2.4). Imbalance is not a usable probe there.
- **2D: no rigorous localization** for:
  - random π-flux;
  - supercritical quantum percolation;
  - the Bernoulli Anderson model at generic energies.

  What exists are band-edge results (Ding–Smart 2020 on Z²; Li–Zhang 2022 on Z³) and, unverified, a 2D large-disorder result by L. Li. The fractional-moment method (FMM) does not apply to Bernoulli disorder.
- **What dynamical localization buys for quasi-free fermions.** Dynamical localization in expectation, uniform in finite volume, gives uniform-in-time bounds on sector-wise entanglement and on number fluctuations (ANSS-type, elementary).
  - For a superposition over sectors, number fluctuations stay bounded.
  - The full entanglement entropy is not controlled by the sector-wise bounds. It is expected to grow like ln t, because sectors dephase against each other (l-bit-like charges).
- **Prethermalization.**
  - Abanin–De Roeck–Ho–Huveneers (ADHH) give emergent conservation of any integer-spectrum "N" up to times exponential in the large coupling. This covers the U(1) electric energy at strong coupling and the dipole moment in a strong Stark field.
  - I know of no localization theorem for random-walk potentials. Simon–Spencer gives absence of absolutely continuous spectrum for them.

---

## 1. 1D Anderson model with Bernoulli (binary) i.i.d. potential
Setup: (Hψ)(x) = ψ(x+1)+ψ(x−1)+λV_ω(x)ψ(x) on ℓ²(Z). V_ω(x) is i.i.d. with law μ; supp μ has ≥2 points and is compact.

**Furstenberg hypotheses at all energies (P, standard) [memory/deriv].**
- With A(E,v)=[[E−v,−1],[1,0]], one has A(E,v₁)A(E,v₂)^{-1} = [[1,v₂−v₁],[0,1]], a nontrivial parabolic matrix.
- So the group G_E is non-compact. Strong irreducibility follows because A(E,v₁) moves the only line that parabolic fixes.
- Hence γ(E)>0 for every E and every λ≠0. The i.i.d. Bernoulli model has **no critical energies**.
- Constants degenerate as λ→0, where γ ~ λ² (with the Kappus–Wegner anomaly at E=0 affecting only the weak-disorder expansion).

### 1.1 Spectral localization
- **(P) Carmona–Klein–Martinelli 1987** [snip for bibliographic data; content memory].
  - Statement: for any non-degenerate μ (compact support suffices; Bernoulli included), H_ω has almost surely pure point spectrum with exponentially decaying eigenfunctions.
  - Method: Furstenberg + Le Page (Hölder continuity of the integrated density of states, which replaces Wegner) + multiscale analysis.
- **Alternative proofs:**
  - Shubin–Vakilian–Wolff 1998 [memory].
  - Bucaj–Damanik–Fillman–Gerbuz–VandenBoom–Wang–Zhang 2019 and Jitomirskaya–Zhu 2019 (§1.2).
  - Gorodetski–Kleptsyn 2021: a geometric proof via the "parametric Furstenberg theorem" [snip].

### 1.2 Dynamical localization: which form is proven
- **(P) Polynomial moments, Damanik–Sims–Stolz, J. Funct. Anal. 2004 (random word models)** [memory; title snip].
  - Strong dynamical localization in expectation, E sup_t ‖|X|^{p/2} e^{-itH} χ_I(H) δ₀‖² < ∞ for all p.
  - Holds on compact I avoiding a finite critical set; that set is empty for i.i.d. Bernoulli.
  - Built on the Damanik–Stollmann principle "multiscale analysis implies strong dynamical localization" (GAFA 2001) [memory].
- **Continuum analogue:** Damanik–Sims–Stolz, Duke Math. J. 2002, arXiv:math-ph/0010016 [snip id/title]. Spectral and dynamical localization for continuum Bernoulli–Anderson models [memory for content].
- **(P) Bucaj et al., Trans. AMS 2019, arXiv:1706.06135** [snip].
  - Self-contained route: Furstenberg positivity + large deviations for the Lyapunov exponent.
  - Gives spectral and dynamical localization; the dynamical form is polynomial/almost sure [memory for the exact form].
- **(P) Jitomirskaya–Zhu, CMP 370 (2019) 311–324, arXiv:1803.10697** [snip].
  - Large-deviation proof without multiscale analysis.
  - Gives spectral localization and dynamical localization via a uniform large-deviation estimate [snip]. Exact form not verified.
- **(P) Ge–Zhao, J. Spectral Theory 10 (2020) 887–904, DOI 10.4171/JST/315, arXiv id? (submitted June 2018)** [snip for "exponential dynamical localization in expectation via positivity + uniform large deviations"; exact statement from memory].
  - Statement: for any non-degenerate compactly supported μ (Bernoulli included), every λ≠0 and all energies, there are C,c>0 with
    **E[ sup_t |<δ_x, e^{-itH} δ_y>| ] ≤ C e^{-c|x−y|} for all x,y ∈ Z.**
  - This is exactly the form you asked about. It is stated for the infinite-volume operator [memory].
- **(P) Macera–Sodin, Ann. Henri Poincaré 23 (2022) 4227–4247, arXiv:2110.00097** [snip].
  - Reproves the Klein–Lacroix–Speis 1990 strip theorem without multiscale analysis, and extends it to general quasi-1D models, random hopping allowed.
  - Bounds the **eigenfunction correlator**, E sup_{|f|≤1} |<δ_x, f(H) δ_y>| ≲ e^{-c|x−y|} [snip: "bounds the eigenfunction correlator ⇒ exponential dynamical localization and exponential decay of the Fermi projection"].
  - Also covers half-line operators with arbitrary deterministic self-adjoint boundary conditions [snip].
  - The exact hypotheses (Furstenberg-group conditions on an energy interval) are from memory.
- **(P) Rate is not the Lyapunov exponent.** Macera, J. Math. Phys. 65 (2024) 012103, arXiv:2307.07626 [snip].
  - The annealed decay rate c differs from γ in general. The first counterexample is Anderson–Bernoulli, and equality never holds for bounded potentials at large disorder.
  - So do not identify c with 1/ξ_Lyapunov in proofs.
- **(P) Non-stationary / correlated extensions** [snip]:
  - Gorodetski–Kleptsyn, "Non-stationary Anderson localization", Commun. AMS 5 (2025) 81–143, arXiv:2403.15970. Independent, non-identically distributed, uniformly non-degenerate bounded potentials, including non-stationary Anderson–Bernoulli. Spectral and dynamical localization.
  - Damanik–Gorodetski–Kleptsyn, "Localization for random Schrödinger operators defined by block factors", arXiv:2504.08153 (2025). Potentials given by a block code applied to an i.i.d. sequence; no regularity, monotonicity or continuity assumptions. Spectral localization, plus dynamical localization away from a finite exceptional set.
  - arXiv:2512.18520, "Localization of the 1D non-stationary Anderson model" (2025). Removes uniform boundedness (spectral localization) [snip; author unconfirmed].
  - Rangamani, "Exponential dynamical localization for random word models", AHP 23 (2022) 4171–4193, arXiv:1910.10214 [snip]. Exponential dynamical localization **in expectation** on any compact set avoiding a finite set of critical energies; covers random polymer and generalized Anderson models.
  - Boumaza, "Localization for random quasi-one-dimensional models", J. Math. Phys. 2023, arXiv:2305.05224 [snip]. Review plus a localization criterion through the Furstenberg group; all its Schrödinger-type results allow Bernoulli randomness.
  - Liu–Shi–Zhang, arXiv:2607.03472 (2026) [snip]. Anderson–Bernoulli on Z with long-range hopping: localization near the spectral edge for rational Laurent symbols, via quantitative unique continuation + multiscale analysis. Billed as the first long-range result with pure Bernoulli potential.
- **Ge–You–Zhou, "Exponential dynamical localization: criterion and applications", arXiv:1901.04258** [snip]. A general criterion for ergodic operators.

### 1.3 Are the constants uniform in finite volume?
- ANSS-type entanglement bounds (§4) need the following hypothesis:
  **(EL)** E sup_t |<δ_x, e^{-itH_Λ} δ_y>| ≤ C e^{-c|x−y|}, uniformly over intervals Λ ∋ x,y.
- **Ge–Zhao, Jitomirskaya–Zhu, Bucaj et al. and Gorodetski–Kleptsyn** state infinite-volume results [memory].
- **Macera–Sodin** work with half-line and finite-interval operators with arbitrary boundary conditions. I believe their correlator bound is uniform in the box [memory, verify in their Sect. 2–3].
- **For absolutely continuous μ**, the Kunz–Souillard method (CMP 78, 1980) and the fractional-moment method give (EL) directly [memory].
- **Assessment (H):** (EL) for Bernoulli is almost certainly true and provable by the large-deviation or Macera–Sodin route, but you should either locate the exact statement or write the finite-volume version. ANSS and related papers cite (EL) only for absolutely continuous distributions [memory].

### 1.4 Binary random hopping signs; correlated binary potentials
- **(P, elementary) Random signs are pure gauge in 1D.** t_j = ±t is removed by c_j → (∏_{k<j} s_k) c_j on an open chain. On a ring, a single flux ∏s_j = ±1 remains, which only switches periodic/antiperiodic boundary conditions. Dynamics is ballistic and there is no localization.
- **Random dimer / random polymer:**
  - (H/N) Dunlap–Wu–Phillips, PRL 65, 88 (1990): ~√N delocalized states near critical energies, ⟨x²⟩ ~ t^{3/2} [memory].
  - (P) De Bièvre–Germinet, arXiv:math-ph/9907006 [snip], J. Stat. Phys. 98 (2000) [memory for journal]. Dynamical localization on every energy interval that contains no critical energy.
  - (P) Jitomirskaya–Schulz-Baldes–Stolz, "Delocalization in random polymer models", CMP 233 (2003) 27–48, arXiv:math-ph/0405024 [snip]:
    - pure point spectrum with exponentially localized eigenstates almost surely;
    - at a generic critical energy, γ vanishes quadratically and the density of states is positive;
    - time-averaged moments satisfy ⟨|X|^q⟩_T ≳ T^{q−1/2−α} for all α>0, so almost surely overdiffusive for q=2.
  - Damanik–Sims–Stolz JFA 2004 and Rangamani 2022 give localization away from the critical set.
  - (P) "Eigenvalue statistics for random polymer models: localization and delocalization", arXiv:2508.21215 [snip title only].
  - I do not know a random-polymer paper by "Damanik–Lenz". The relevant one is Damanik–Sims–Stolz 2004.
- **Consequence for gauge models (P via JSBS):** if the charge-sector measure produced by your initial state has dimer-type correlations, uniform dynamical localization fails and transport is superdiffusive at the critical energies.

---

## 2. Random-hopping (off-diagonal disorder) chains
Model: h = Σ_i J_i(|i⟩⟨i+1|+h.c.), J_i > 0 i.i.d., E|ln J| < ∞.
- Your case J = u^δ gives −ln J ~ δ·Exp(1), so σ² := Var ln J = δ².
- This is the fixed-point family of Fisher's renormalization group with initial width Γ₀ = δ (the parametrization used by Iglói and coworkers) [memory].

### 2.1 Rigorous results
- **(P, elementary) Zero energy is critical.**
  - At E=0 the transfer matrices are anti-diagonal. Two-step products are diagonal, with entries ∏(−J_{2k−1}/J_{2k})^{±1}.
  - So G₀ preserves the two coordinate axes and is not strongly irreducible.
  - γ(0) = |E ln J_odd − E ln J_even| = 0 for i.i.d. hoppings.
  - The zero-energy solution satisfies ln|ψ(2n)| = a mean-zero random walk, so |ψ| ~ e^{±σ√n}. It is not in ℓ², so there is no atom at 0.
- **(P) Localization away from 0.**
  - Delyon–Simon–Souillard, "Localization for off-diagonal disorder and for continuous Schrödinger operators", CMP 109 (1987) 157–165 [memory]: pure point spectrum for absolutely continuous hopping laws, via Kotani theory + spectral averaging. E=0 has spectral measure zero.
  - For general (also singular, e.g. binary) hopping laws: Macera–Sodin 2022 (random hopping allowed) [snip] and Rangamani, "Singular-unbounded random Jacobi matrices", J. Math. Phys. 60 (2019), id? [memory].
  - These give spectral localization and exponential dynamical localization in expectation **on compact energy windows avoiding 0** (and possibly a finite exceptional set) [memory].
- **(P) No theorem on transport at E=0 is known to me.**
  - The Jitomirskaya–Schulz-Baldes–Stolz lower-bound technique needs bounded (elliptic/commuting) transfer matrices at the critical energy. Here they grow like e^{c√n}, so it does not apply.
  - Damanik–Tcheremchantsev upper bounds (JAMS 20, 2007) [memory] should give vanishing transport exponents, but this is not checked (H).
- **(P?) Dyson singularity.** N(E)−½ ∝ 1/ln²|E|, ρ(E) ∝ 1/(|E| |ln|E||³), γ(E) ∝ σ²/|ln|E||.
  - Physics derivations are exact (Dyson, Phys. Rev. 92, 1331 (1953); Theodorou–Cohen PRB 13, 4597 (1976); Eggarter–Riedinger PRB 18, 569 (1978); Balents–Fisher PRB 56, 12970 (1997)) [memory].
  - I could not verify a fully rigorous mathematical paper.
  - Krishna–Bhatt (arXiv:2004.00064, PRB 2020) [snip arXiv/claim; authors memory]: the Dyson form is universal for finite Var ln J and is modified for divergent hopping distributions. Your δ<∞ case is in the universal class.

### 2.2 Physics: infinite randomness and dynamics (H/N)
- Dasgupta–Ma PRB 22, 1305 (1980); Fisher PRB 50, 3799 (1994) [memory]:
  - the random XX chain (equivalent to random-hopping free fermions) flows to the infinite-randomness random-singlet fixed point;
  - ln(1/Ω) ~ σ√L.
- Ground state: Refael–Moore PRL 93, 260602 (2004): S ≃ (ln2/3) log₂ L, a logarithmic violation of the area law [memory].
- Quench entanglement:
  - Iglói–Szatmári–Lin PRB 85, 094417 (2012): S(t) ~ ln ln t at infinite-randomness criticality (random transverse-field Ising chain, free-fermion numerics + strong-disorder RG; the random XX chain maps to two Ising chains) [memory].
  - Vosk–Altman PRL 110, 067204 (2013): time-dependent strong-disorder RG; S ~ ln t needs Jz ≠ 0 [memory].
  - Zhao–Andraschko–Sirker PRB 93, 205146 (2016): global quenches in disordered XX/XXZ chains, including bond disorder [memory; their precise XX-bond result not verified].
  - Vasseur–Potter–Parameswaran PRL 114, 217201 (2015): RSRG-X, infinite randomness of excited states [memory].
- Number entropy:
  - Interacting MBL: Kiefer-Emmanouilidis–Unanyan–Fleischhauer–Sirker PRL 124, 243601 (2020) and PRB 103, 024203 (2021); counterpoints Luitz–Bar Lev PRB 102, 100202 (2020) and Ghosh–Žnidarič PRB 105, 144203 (2022) [memory].
  - **I know of no dedicated study of number entropy in free random-hopping chains** [memory].
- **My estimates (H)** from the renormalization-group scaling, for a site-localized initial state:
  - L(t) ~ ln²(t)/σ²;
  - ⟨x²⟩_t ~ ∫^{ln t} dℓ ℓ^{-3}·ℓ⁴ ~ ln² t, using dN ~ dℓ/ℓ³ with ℓ = ln(1/E), and L_E ~ ℓ²;
  - ⟨|x|⟩_t ~ ln ln t;
  - E sup_t |<x|e^{-ith}|y>| decays only polynomially, roughly 1/|x−y| up to logs. Summed over a cut, this gives ln L(t) ~ ln ln t, consistent with S(t) ~ ln ln t;
  - Var N_A(t) ~ ln ln t;
  - number entropy S_N ≲ ½ ln Var ~ ½ ln ln ln t.
  - At t=10¹²: ln ln t ≈ 3.3 and ln ln ln t ≈ 1.2. Unbounded growth is numerically very hard to distinguish from saturation.

### 2.3 Practical points for your Fortran code (P for the numerical-analysis fact; H for the scales)
- **Use high-relative-accuracy routines.**
  - Zero-diagonal tridiagonal h ≅ [[0,B],[Bᵀ,0]] with B bidiagonal, so eigenvalues = ±(singular values of B).
  - Late-time dynamics is dominated by |E| ~ e^{−σ√L}, and with J=u^δ the hoppings span many decades.
  - Standard dense or tridiagonal eigensolvers give only absolute accuracy ~ε_mach‖h‖.
  - Use bidiagonal SVD with high *relative* accuracy: Demmel–Kahan, SIAM J. Sci. Stat. Comput. 11 (1990) 873 [memory]; LAPACK DBDSQR, or the dqds algorithm DLASQ1.
  - Otherwise the phases e^{−iE t} at t ≲ 10¹² for the critical modes are wrong.
- **Finite-size saturation time (H).** t_L ~ exp(c σ√L). Reaching t=10¹² unsaturated needs σ√L ≳ 28, i.e. L ≳ 800/δ² up to an O(1) factor.

### 2.4 Imbalance identity for chiral models [deriv] (P, elementary)
- Assume a non-degenerate spectrum and a Fock initial state with occupations n_z = (1+(−1)^z)/2.
- The time-averaged imbalance is Ī = (1/2L) Σ_α a_α², with a_α = Σ_x (−1)^x |φ_α(x)|². (This uses Σ_α a_α = 0.)
- In a bipartite chiral model, every eigenvector with E≠0 has equal weight on the two sublattices, so a_α = 0 and **Ī = 0**.
- This holds for random-hopping chains, 2D random π-flux and quantum percolation without potential, in localized and delocalized regimes alike.
- So density-wave imbalance is useless there. Entanglement or number variance is needed.

---

## 3. Two dimensions
### 3.1 Random-sign hopping / random π-flux on Z²
- Real hoppings ±t, bipartite: symmetry class BDI. Independent uniform link signs give i.i.d. uniform plaquette fluxes [deriv].
- **(H/N) Physics:**
  - Gade–Wegner NPB 360, 213 (1991); Gade NPB 398, 499 (1993): critical, delocalized states at E=0 and a diverging density of states.
  - Corrected strong-disorder form ρ(E) ~ |E|^{-1} exp(−c|ln|E||^{2/3}): Motrunich–Damle–Huse PRB 65, 064206 (2002); Mudry–Ryu–Furusaki PRB 67, 064202 (2003).
  - States at E≠0 are localized (2D orthogonal class).
  - All [memory].
- **(P) Rigorous results** [memory]:
  - Lifshitz tails for 2D discrete random magnetic field: Nakamura, AHP 1 (2000) 823.
  - Localization near band edges: Klopp–Nakamura–Nakano–Nomura, AHP 4 (2003) 795. Their assumptions on the flux distribution very likely exclude binary Z2 fluxes (Wegner estimate); verify.
- **Open:** localization at E≠0 away from the edges, and any delocalization statement at E=0.
  - The only rigorous 2D delocalization mechanism I know is topological: Germinet–Klein–Schenker, Ann. Math. 166 (2007), random Landau Hamiltonians [memory]. It does not cover BDI.
- **Quasi-1D (ladders/strips) with π-flux:**
  - Brouwer–Mudry–Simons–Altland PRL 81, 862 (1998): chiral strips with an odd number of channels are delocalized at E=0 [memory].
  - So strip theorems (Klein–Lacroix–Speis JFA 94 (1990); Macera–Sodin; Boumaza) apply only away from E=0 for chiral ladders.

### 3.2 Quantum percolation
- (P) Chayes–Chayes–Franz–Sethna–Trugman, J. Phys. A 19, L1173 (1986): compactly supported eigenstates give a dense set of jumps in the integrated density of states [memory].
- (P) Veselić, "Spectral analysis of percolation Hamiltonians", Math. Ann. 331 (2005) 841: IDS, its jumps, and their location at energies of compactly supported eigenfunctions [memory].
- (P) Kirsch–Müller, Math. Z. 252 (2006) 899: Lifshitz tails for Laplacians on bond-percolation graphs [memory].
- (P) Müller–Stollmann, JFA 252 (2007) 233 (supercritical spectral asymptotics), and the review "Percolation Hamiltonians" (Birkhäuser, 2011) [memory].
- **(P, trivial) Subcritical percolation:** all clusters are finite almost surely. This gives pure point spectrum and E sup_t |U_{xy}| ≤ P(x↔y) ≤ e^{−c|x−y|}.
- **Supercritical in d=2:** no localization proof known to me; numerics conflict (N) [memory].

### 3.3 2D Anderson–Bernoulli and the fractional-moment method
- **(P) Bourgain–Kenig, Invent. Math. 161 (2005) 389** [memory]: continuum Bernoulli, localization near the bottom of the spectrum in any dimension.
- **(P) Ding–Smart, "Localization near the edge for the Anderson Bernoulli model on the two dimensional lattice", Invent. Math. 219 (2020) 467–506, arXiv:1809.09041** [memory].
  - Statement: for H = −Δ + λV on Z², V i.i.d. Bernoulli{0,1}, any λ>0, there is ε>0 such that H has almost surely Anderson localization in [0, ε].
  - Main tool: discrete unique continuation (Buhovsky–Logunov–Malinnikova–Sodin, Duke 2022, arXiv:1712.07902 [memory]).
- **(P) Li–Zhang, Duke Math. J. 171 (2022), arXiv:1906.04350** [memory]: the Z³ analogue (edge).
- **(P?) L. Li, "Anderson–Bernoulli localization at large disorder on the 2D lattice", CMP (2022), id?** [memory]. Exact energy range not verified; check before relying on it.
- **(P) Imbrie, Rev. Math. Phys. 33 (2021), arXiv:1705.01916** [memory]: localization and level-spacing bounds in any d for discrete distributions taking N values with N large, at large disorder. **Not binary.**
- **Recent titles, content not verified** [snip titles only]:
  - non-stationary models on Z²: arXiv:2405.10636;
  - non-stationary models in 3D: arXiv:2603.17810;
  - long-range hopping: arXiv:2412.17262, arXiv:2508.12714.
- **Fractional-moment method** (Aizenman–Molchanov CMP 157 (1993); Aizenman–Schenker–Friedrich–Hundertmark CMP 224 (2001)) [memory]:
  - the a priori bound E|G(x,y;z)|^s ≤ C needs τ-regular (e.g. Hölder or bounded-density) single-site laws;
  - with atoms, E|V−α|^{−s} = ∞ for α at an atom;
  - **so it does not apply to Bernoulli at any disorder strength.** In 1D one avoids it with Furstenberg + large deviations; for d≥2 no general replacement exists.
- **Strong binary disorder heuristic (H):**
  - on the square lattice with p=1/2, both "+" and "−" site clusters are subcritical, since p_c^site(Z²) ≈ 0.5927 > 1/2 [memory];
  - so for h ≫ J, eigenstates live on finite clusters, coupled by tunnelling ~(J/h)^{dist};
  - a cluster-expansion proof looks feasible, except for exact degeneracies between identical clusters, which is the Bernoulli Wegner obstruction again;
  - on the triangular lattice p=1/2 is critical, so this route fails.

---

## 4. Dynamical localization ⇒ bounds on entanglement and number fluctuations (quasi-free fermions)
### 4.1 Literature
- **(P) Pastur–Slavin, PRL 113, 150404 (2014), id?** [memory]. If E|P_{E_F}(x,y)| ≤ Ce^{−c|x−y|} (fractional-moment regime), the mean ground-state entanglement entropy obeys an area law: E S(Λ_L) ≤ C L^{d−1}.
- **(P) Elgart–Pastur–Shcherbina, J. Stat. Phys. 166 (2017) 1092, arXiv:1601.00294** [memory].
  - In 1D, E S_L stays bounded and S_L is not self-averaging (O(1) fluctuations).
  - For d ≥ 2, S/L^{d−1} self-averages.
  - Assumes exponential decay of the Fermi projection in expectation.
- **(P) Hamza–Sims–Stolz, CMP 315 (2012) 215, arXiv:1108.3811** [memory]. Zero-velocity Lieb–Robinson bound for the random XY chain; uses the fractional-moment method, so an absolutely continuous field distribution.
- **(P) Abdul-Rahman–Stolz, J. Math. Phys. 56 (2015) 121901, id?** [memory]. Uniform area law for all eigenstates of the disordered XY chain under eigencorrelator localization.
- **(P) Abdul-Rahman–Nachtergaele–Sims–Stolz, Lett. Math. Phys. 106 (2016) 649, arXiv:1510.00262** [memory]. Under (EL) uniform in finite volume, E sup_t S(ρ_A(t)) ≤ C uniformly in |A| and chain length, for product initial states.
- **(P?) Müller–Pastur–Schulte, "How much delocalisation is needed for an enhanced area law of the entanglement entropy?", CMP 376 (2020), id?** [memory]. Random dimer model: as I recall, delocalization at the critical energies does not produce a logarithmic enhancement. Verify.
  - Contrast: the random-hopping chain at E_F=0 does show a ln L enhancement (Refael–Moore, H/N).
- **Gebert–Müller, "Localization for random block operators", arXiv:1208.1947** [snip id/title; authors memory]. Fractional-moment method for Bogoliubov–de Gennes-type block operators (regular laws).

### 4.2 Elementary chain of inequalities usable directly [deriv] (P)
Setting: pure Gaussian state with correlation matrix C(t) = U C₀ U†, U = e^{−ith}, region A with complement B.
1. **Number variance:** Var N_A = Tr C_A(1−C_A) = Σ_{x∈A, y∈B} |C_{xy}|².
2. **Entropy upper bound:**
   - h(ν) ≤ 2ln2 √(ν(1−ν)) (checked numerically; tight at ν=½);
   - for a projection C, √(C_A(1−C_A)) has trace norm ‖P_A C P_B‖₁ ≤ Σ_{x∈A,y∈B} |C_{xy}|;
   - hence **S_A ≤ 2ln2 Σ_{x∈A,y∈B} |C_{xy}(t)|**.
3. **Entropy lower bound:** h(ν) ≥ 4ln2·ν(1−ν), so **S_A ≥ 4ln2·Var N_A**. There is no converse: bounded variance does not imply bounded entropy.
4. **Number entropy:** S_N ≤ ½ ln(2πe(Var N_A + 1/12)) (maximum entropy for integer-valued variables).
5. **Density-wave or Fock initial state:**
   - |C_{xy}(t)| ≤ Σ_z |U_{xz}||U_{yz}|, and E sup_t |U_{xz}||U_{yz}| ≤ min(E sup_t|U_{xz}|, E sup_t|U_{yz}|);
   - under (EL), E sup_t |C_{xy}(t)| ≤ C'(1+|x−y|) e^{−c|x−y|/2};
   - hence **E sup_t S_A(t) ≤ C''·#∂A and E sup_t Var N_A(t) ≤ C''·#∂A**, uniformly in time and volume.
   - Only (EL) is needed. No regularity of the distribution enters beyond (EL) itself.

### 4.3 Superpositions over charge sectors [deriv/H]
Write |Ψ(t)⟩ = Σ_q c_q |q⟩⊗|ψ_q(t)⟩, with local, conserved, commuting charges q and quasi-free sector dynamics h_q.
- **Number fluctuations (P).** N_A commutes with the q's, so its distribution is the mixture over sectors:
  - Var_Ψ N_A(t) = E_q Var_q + Var_q(⟨N_A⟩_q);
  - both terms are bounded uniformly in t under (EL), since |⟨N_A⟩_q(t) − ⟨N_A⟩_q(0)| ≤ Σ_{x∈A,y∈B} sup_t|U_{xy}|² terms;
  - so **number fluctuations and number entropy stay bounded**.
- **Mixed sector states (P).** For a mixture ρ = Σ p_q ρ_q, use Σ p S(ρ_q) ≤ S(Σ p ρ_q) ≤ Σ p S(ρ_q) + H(p). When the labels are locally readable, the left inequality is an equality plus H(p_A).
- **Full entanglement entropy of Ψ (H).** It is *not* bounded by the sector-wise bounds.
  - Coherences between sectors that differ at a charge ℓ sites from the cut decay once the R-reduced matter states differ. By Uhlmann, ‖Tr_R|ψ_q⟩⟨ψ_{q'}|‖₁ = F(ρ_R^q, ρ_R^{q'}).
  - Duhamel gives δU = −i∫₀^t U'(t−s)(h'−h)U(s)ds.
  - With (EL) and the symmetry of the Bernoulli(½) law under flipping q_j: E sup_{s≤t} |(δC_R)| ≲ h t e^{−cℓ}.
  - So charges within ℓ(t) ~ c^{−1} ln t of the cut dephase, giving **S(t) ≲ C ln(t|A|)** (plausible rigorous upper bound, sketch) and growth ~ ln t (MBL-like, the charges act as l-bits).
  - Consistent with Hart–De Tomasi–Castelnovo, "Logarithmic entanglement growth from disorder-free localization in the two-leg compass ladder", PRL 126, 227202 (2021) [memory].
- **Original vs. dual variables (P, elementary sketch).** For a half-chain cut, a Wilson-line/Jordan–Wigner map anchored at the left end acts across the cut as a single controlled-Z2 gate. Its operator Schmidt rank is 2, so it changes the entanglement entropy by at most 2ln2.

---

## 5. Prethermalization, Stark, and random-walk potentials
### 5.1 Exponentially long emergent conservation laws
- **(P) Abanin–De Roeck–Ho–Huveneers, "A rigorous theory of many-body prethermalization for periodically driven and closed quantum systems", CMP 354 (2017) 809–827, arXiv:1509.05386** [memory].
  - Setting: H = J·N + Y, where N = Σ_x N_x is a sum of commuting, finite-range terms with integer spectrum, Y is local with bounded local norm, and J ≫ 1.
  - Result: there is a quasi-local unitary U, close to 1, with U*HU = JN + D + V, [D,N]=0, and ‖V‖_local ≤ C e^{−cJ/(log J)³}.
  - Consequence: N is conserved up to times ~e^{cJ/log³J} (log corrections removed in some settings).
- **(P) Related theorems** [memory]:
  - Else–Bauer–Nayak, PRX 7, 011026 (2017), arXiv:1607.05277 (prethermal phases protected by emergent symmetry);
  - Abanin–De Roeck–Huveneers, PRL 115, 256803 (2015), arXiv:1507.01474;
  - Kuwahara–Mori–Saito, Ann. Phys. 367, 96 (2016), arXiv:1508.05797;
  - Mori–Kuwahara–Saito, PRL 116, 120401 (2016);
  - De Roeck–Verreet, arXiv:1911.01998 (weak driving);
  - Yin–Lucas, PRL 131, 050402 (2023), id? (local robustness of gapped systems).
- **Applications to gauge theories (H/N)** [memory]:
  - Halimeh–Hauke, PRL 125, 030503 (2020) ("reliability", energy-penalty protection);
  - Halimeh et al., PRX Quantum 2, 040311 (2021) (single-body protection);
  - Halimeh et al., PRX Quantum 3, 020345 (2022) (enhanced disorder-free localization through emergent local symmetries);
  - Lang et al., PRB 106, 174203 (2022) (Stark gauge protection).
- **Direct uses (P/H):**
  - **U(1) at strong coupling:** Σ_n E_n² has integer spectrum, so electric energy (string length) is conserved up to exp(c g²/w) times.
  - **Strong tilt F:** with N = Σ_x x n_x, the dipole moment is conserved up to exp(cF/J). The local terms are unbounded but ad_N acts boundedly on local operators; a routine adaptation, used e.g. by Khemani–Hermele–Nandkishore PRB 101, 174204 (2020) [memory].
  - These give long prethermal constrained dynamics (fragmentation), **not localization**.

### 5.2 Wannier–Stark
- **(P, classical) Single particle on the lattice** [memory]. H = hopping + F x:
  - spectrum FZ, pure point;
  - eigenvectors are Bessel functions J_{x−n}(2/F);
  - <x|e^{−itH}|y> = phase × J_{x−y}((4/F) sin(Ft/2));
  - so sup_t |U_{xy}| ≤ (2/F)^{|x−y|}/|x−y|!: deterministic, superexponential dynamical localization.
- **Stark many-body localization (N/E)** [memory]:
  - numerics: Schulz–Hooley–Moessner–Pollmann PRL 122, 040606 (2019); van Nieuwenburg–Baum–Refael PNAS 116, 9269 (2019);
  - experiments (E): Morong et al. Nature 599, 393 (2021); Scherg et al. Nat. Commun. 12, 4490 (2021); Guardado-Sanchez et al. PRX 10, 011042 (2020);
  - fragmentation: Sala et al. PRX 10, 011047 (2020).
- **No rigorous Stark-MBL theorem is known to me.** Imbrie's MBL proof (J. Stat. Phys. 163 (2016) 998) assumes i.i.d. random fields with a density plus a level-attraction assumption, and covers neither Bernoulli, Stark nor long-range settings [memory].

### 5.3 Random-walk potentials V(x) = λW_x, W_x = Σ_{j≤x} q_j
- **Where they come from in gauge theories [deriv]:**
  - **U(1) with E² energy, tightly bound meson** (dipole at x): energy g²(W_x + ½), so the potential is a random walk.
  - **U(1), isolated charge with partner fixed far away:** energy g²(x/2 + Σ_{n<x} W_n). That is string tension plus an integrated random walk, so the *force* is a random walk and V ~ |x|^{3/2}.
  - **1D Z2 theory with the standard single-link electric term −hΣσ^x_b:**
    - gauge fixing gives σ^x_b = σ^x_0 ∏_{k≤j} q_k (−1)^{n_k};
    - a pair of fermions at x<y costs 2h Σ_{j∈[x,y)} Q_j, with Q_j = ∏_{k≤j} q_k;
    - for i.i.d. uniform q, the Q_j are i.i.d. ±1, so the pair potential is a random walk;
    - for uniform q this reduces to linear confinement (Borla–Verresen–Grusdt–Moroz PRL 124, 120503 (2020) [memory]).
- **(P) Simon–Spencer, CMP 125 (1989) 113** [memory]. If limsup_{n→±∞} |V(n)| = ∞ on both sides, then σ_ac = ∅. This applies almost surely to λW_x and to integrated random walks: **no absolutely continuous spectrum**.
- **Pure point spectrum or dynamical localization: no theorem known to me (open).** Nearest rigorous results:
  - non-stationary Anderson localization (independent values, bounded: Gorodetski–Kleptsyn 2025; unbounded: arXiv:2512.18520 [snip]);
  - heavy-tailed i.i.d. potentials (arXiv:2505.09079 [snip title]);
  - Stark ladders.

  None covers values with independent increments.
- **Heuristic (H):** off resonances, |ψ(x)| ~ exp(−½|x| ln|x|), superexponential like Wannier–Stark, with sparse resonant returns of the walk to the level E/λ.
- **Proof routes (H):**
  - **(a) Increments with a bounded density.** Changing q_a shifts V *rigidly* on [a,b], so σ(H_{[a,b]}) moves rigidly. This gives an exact one-variable Wegner estimate P(dist(σ(H_Λ),E)<ε) ≤ 2‖ρ‖_∞|Λ|ε [deriv]. Combined with multiscale analysis, whose initial scale comes from |V−E| ≫ 1 on most sites, a proof looks within reach.
  - **(b) Integer increments (physical charges).** Exact degeneracies and no Wegner estimate: a Bernoulli-type difficulty without stationarity, so Furstenberg/Le Page are unavailable. Genuinely open.
  - Damanik–Tcheremchantsev-type upper bounds should still yield zero transport exponents (H).

---

## 6. Assessment: what plugs directly into a proof
### (i) 1D Z2 gauge theory with free fermionic matter (Smith–Knolle–Kovrizhin–Moessner, PRL 118, 266601 (2017), id?)
- **Which model is free.** [memory] SKKM's electric term is the two-link product hΣσ^x_{j−1,j}σ^x_{j,j+1}. Gauss's law, q_j = σ^xσ^x(−1)^{n_j}, then gives hΣ q_j(1−2n_j) [deriv].
  - So the sector Hamiltonian is **−JΔ + 2h q_j, an Anderson–Bernoulli chain**.
  - From a σ^z-polarized gauge state plus a Fock state for the fermions, q is i.i.d. Bernoulli(½) and the fermion initial state does not depend on the sector [deriv].
  - Check that this is SKKM's exact Hamiltonian. With the *standard* term −hΣσ^x the model is interacting (random-walk pair potentials, §5.3).
- **What plugs in directly:**
  - Carmona–Klein–Martinelli: spectral localization for all h/J>0 and all energies.
  - Ge–Zhao and Macera–Sodin: E sup_t |U_{xy}| ≤ Ce^{−c|x−y|}.
  - §4.2: bounded sector-wise entanglement and number fluctuations.
  - §4.3: bounded number fluctuations of the full superposition.
  - §2.4 does not apply (there is a potential), so imbalance persistence follows from the pure point spectrum.
- **To be supplied:**
  1. (EL) with constants uniform in finite volume for Bernoulli: likely routine via Macera–Sodin, but write or locate it.
  2. Correlated sector measures for other initial states: block-factor/random-word theorems give localization, but dimer-type correlations create critical energies with provable superdiffusion (JSBS).
  3. Non-local gauge-invariant observables (electric-field strings) need bounds on determinants/Pfaffians of localized correlation matrices.
  4. The full entanglement entropy is expected to grow ~ln t (sector dephasing). The best one can hope to prove is an upper bound ~ln t, not uniform boundedness.
- **Quasi-1D versions** (ladders, compass ladder): Klein–Lacroix–Speis / Macera–Sodin / Boumaza cover Bernoulli strips, except at chiral E=0 points (Brouwer et al.).

### (ii) 2D versions
- **Static gauge field** (no single-link σ^x): random π-flux fermions (class BDI). E=0 is expected to be critical, so no uniform dynamical localization; rigorous results only near the band edges (and probably not for binary flux).
- **With on-site binary potential from charges:** 2D Bernoulli Anderson model. Only band-edge theorems exist (Ding–Smart; Bourgain–Kenig in the continuum); the fractional-moment method is excluded; L. Li's large-disorder result needs checking.
- **Quantum percolation:** only the subcritical case is trivial.
- **Interacting 2D gauge models** (N): Karpov–Verdel–Huang–Schmitt–Heyl PRL 126, 130401 (2021); Chakraborty–Heyl–Karpov–Moessner PRB 106, L060308 (2022), a percolation picture [memory]. (E) Google Quantum AI 2024 preprint on disorder-free localization in a (2+1)D gauge theory, authors/id? [memory].
- **Verdict: open at the rigorous level.** The best target is strong binary disorder on the square lattice, via subcritical clusters (§3.3).

### (iii) 1D U(1) (Schwinger / quantum link)
- After Gauss's law, matter has a linear (confining) long-range interaction, so it is **not quasi-free**. None of §1–§4 applies beyond few-body sectors.
  - Numerics (N): Brenes–Dalmonte–Heyl–Scardicchio PRL 120, 030601 (2018), log entanglement growth with random background charges in exact diagonalization; Papaefstathiou–Smith–Knolle PRB 102, 165132 (2020); Giudici–Mazza–Scardicchio (2020) [memory, ids?].
  - Confinement-induced slow dynamics (H/N): Kormos et al. Nat. Phys. 13, 246 (2017); Lerose et al. PRB 102, 041118 (2020) [memory].
- **Rigorous tools available:**
  - ADHH: prethermal conservation of electric energy at strong coupling;
  - Wannier–Stark: two-body sector with uniform background, giving the relative coordinate exactly;
  - Simon–Spencer: no absolutely continuous spectrum for random-walk or integrated-random-walk single-particle reductions.
- **Missing:**
  - any localization theorem for random-walk potentials (open, §5.3; route (a) is likely provable with continuous increments, but physical integer charges give route (b));
  - any many-body localization result with long-range confining interactions.
- **Verdict: no plug-in theorem.** Rigorous statements are limited to prethermal (exponential-time) conservation laws and few-body spectral facts.

---

## 7. Reference list (status)
- Carmona, Klein, Martinelli, "Anderson localization for Bernoulli and other singular potentials", CMP 108, 41–66 (1987) [snip].
- Shubin, Vakilian, Wolff, GAFA 8, 932 (1998) [memory].
- Damanik, Sims, Stolz, "Localization for one dimensional, continuum, Bernoulli–Anderson models", Duke Math. J. 114 (2002), arXiv:math-ph/0010016 [snip id; journal memory].
- Damanik, Sims, Stolz, "Localization for discrete one-dimensional random word models", JFA 208 (2004) 423, id? [title snip; rest memory].
- Bucaj, Damanik, Fillman, Gerbuz, VandenBoom, Wang, Zhang, Trans. AMS 372 (2019), arXiv:1706.06135 [snip].
- Jitomirskaya, Zhu, CMP 370, 311–324 (2019), arXiv:1803.10697 [snip].
- Ge, Zhao, J. Spectral Theory 10, 887–904 (2020), DOI 10.4171/JST/315, id? [snip].
- Gorodetski, Kleptsyn, Adv. Math. 378, 107522 (2021), arXiv:1809.00416 [snip].
- Gorodetski, Kleptsyn, Commun. AMS 5, 81–143 (2025), arXiv:2403.15970 [snip].
- Damanik, Gorodetski, Kleptsyn, arXiv:2504.08153 (2025) [snip].
- Rangamani, AHP 23, 4171–4193 (2022), arXiv:1910.10214 [snip].
- Rangamani, "Singular-unbounded random Jacobi matrices", J. Math. Phys. 60 (2019), id? [memory].
- Macera, Sodin, AHP 23, 4227–4247 (2022), arXiv:2110.00097 [snip].
- Macera, J. Math. Phys. 65, 012103 (2024), arXiv:2307.07626 [snip].
- Boumaza, J. Math. Phys. (2023), arXiv:2305.05224 [snip].
- Liu, Shi, Zhang, arXiv:2607.03472 (2026) [snip].
- Ge, You, Zhou, arXiv:1901.04258 [snip].
- Jitomirskaya, Schulz-Baldes, Stolz, CMP 233, 27–48 (2003), arXiv:math-ph/0405024 [snip].
- De Bièvre, Germinet, arXiv:math-ph/9907006 (J. Stat. Phys. 98 (2000) [memory]) [snip].
- Delyon, Simon, Souillard, CMP 109, 157 (1987) [memory].
- Krishna, Bhatt, arXiv:2004.00064 (PRB 2020) [snip id; authors memory].
- Ding, Smart, Invent. Math. 219 (2020), arXiv:1809.09041 [memory].
- Li, Zhang, Duke Math. J. 171 (2022), arXiv:1906.04350 [memory].
- Imbrie, Rev. Math. Phys. 33 (2021), arXiv:1705.01916 [memory].
- Abdul-Rahman, Nachtergaele, Sims, Stolz, LMP 106 (2016), arXiv:1510.00262 [memory].
- Elgart, Pastur, Shcherbina, J. Stat. Phys. 166 (2017), arXiv:1601.00294 [memory].
- Hamza, Sims, Stolz, CMP 315 (2012), arXiv:1108.3811 [memory].
- Abanin, De Roeck, Ho, Huveneers, CMP 354 (2017), arXiv:1509.05386 [memory].
- Else, Bauer, Nayak, PRX 7, 011026 (2017), arXiv:1607.05277 [memory].
- Simon, Spencer, CMP 125, 113 (1989) [memory].
- Further physics references are as cited inline, all [memory]. Ids are omitted where not certain.

**Not checked because the search budget ran out:** the exact theorem forms in Ge–Zhao, Jitomirskaya–Zhu and Macera–Sodin (finite-volume uniformity); L. Li's 2D large-disorder statement; the SKKM Hamiltonian; Hart–De Tomasi–Castelnovo; Zhao–Andraschko–Sirker; Iglói–Szatmári–Lin; Klopp–Nakamura–Nakano–Nomura's assumptions; Müller–Pastur–Schulte. A follow-up run with arXiv fetching enabled should check these first.
