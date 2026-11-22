# What is proven about MBL and slow thermalization in interacting disordered or quasi-periodic lattice systems (as of Nov 2026)

**Limits of this research:**
- arxiv.org (and every other host I tried) could not be reached with the page-fetch tool in this session.
- The shared web-search budget ran out partway through.
- So I read no paper body and no complete abstract.

**Tags used below:**
- **[snippet]**: wording of an abstract or excerpt seen in this session through search-engine results. I did not read the full paper.
- **[memory]**: from my prior knowledge, not checked this session. Treat precise constants and journal details as unverified.
- **(P)** proven theorem; **(N)** numerical; **(H)** heuristic or phenomenological; **(E)** experiment.

---

## (a) Bottom line

1. **No unconditional proof of MBL in the thermodynamic limit.**
   - MBL here means infinite-time localization at positive energy density, in a short-range interacting disordered or quasi-periodic lattice Hamiltonian.
   - The only full-spectrum, all-times result is Imbrie (2016). It covers a 1D Ising-type chain with random fields and couplings and a weak transverse field. It is **conditional on the unproven LLA assumption**.
   - I found no proof of LLA, or of any substitute for it, in a lattice chain.
2. **What is proven without conditions** falls into seven groups:
   - **Quasi-free (XY/XX, Anderson fermions).** Essentially complete: zero-velocity Lieb–Robinson bounds, time-uniform entanglement and particle-transport bounds, and a uniform area law for all eigenstates. All of these hold for all times, uniformly in system size, in expectation.
   - **Interacting, zero energy density.** This covers the random XXZ chain in the Ising phase, in fixed energy windows at the bottom of the spectrum. Proven: localization, quasi-locality, slow propagation, a logarithmic light cone (2026), uniformly in system size and particle number (Elgart–Klein–Stolz; Beaud–Warzel; Elgart–Klein 2022–2026). The Holstein model at low energy (Mavi–Schenker) also belongs here.
   - **Interacting, positive energy density, finite times.** Quasi-localization holds up to t\* ~ exp{c(log W)^{2−ε}} (De Roeck–Huveneers–Meeus–Prośniak 2023).
   - **Interacting, positive energy density, rare events, infinite time.**
     - A finite chain of length L is fully MBL with probability ≥ e^{−cL}.
     - From this: the heat conductance is o(1/L) at any temperature, i.e. "absence of normal conduction" (De Roeck–Giacomin–Huveneers–Prośniak 2024/25).
     - This does **not** need LLA. It is the first unconditional transport statement at infinite time and positive temperature for an interacting disordered chain.
   - **Toy many-body models.**
     - Quantum Sun model: localization and Poisson statistics for small coupling (De Roeck–Hannani 2025).
     - Few-body Hamiltonians on expander or code geometries with a many-body mobility edge (Yin–Nandkishore–Lucas 2024).
   - **Prethermal theorems.** Exponentially long-lived emergent conservation laws when a large term has integer spectrum (Abanin–De Roeck–Ho–Huveneers; Else–Bauer–Nayak; Kuwahara–Mori–Saito). These are disorder-free and give finite times only.
   - **Interacting quasi-periodic fermions (Mastropietro).** Only ground-state (T = 0) statements: decay of correlations and vanishing Drude weight. Nothing about excited states.
3. **No rigorous proof of delocalization or thermalization** for any disordered interacting lattice Hamiltonian. All of the following are heuristic, with numerical or experimental support:
   - avalanche instability;
   - absence of MBL for d > 1, for long-range interactions, in the continuum, and with SU(2) symmetry;
   - absence of many-body mobility edges.
4. **Relevance to the two applications:**
   - U(1) number fluctuations and number entropy are bounded rigorously only in quasi-free models (for all times) and at zero energy density.
   - For interacting chains at positive density, the closest results are the quasi-localization time scale (DHMP) and the rare-region results (DGHP).
   - For disorder-free localization, only the quasi-free Z₂ case (reducible to the 1D Bernoulli Anderson model) and prethermal statements are within current techniques.

---

## (b) Per-result entries

### B1. Imbrie: conditional MBL for a spin chain (P, conditional)

**References**
- J.Z. Imbrie, *On many-body localization for quantum spin chains*, J. Stat. Phys. 163, 998 (2016), arXiv:1403.7837 [snippet: id and journal].
- J.Z. Imbrie, *Diagonalization and many-body localization for a disordered quantum spin chain*, PRL 117, 027201 (2016), arXiv:1605.03003 [snippet].

**Model** [memory]
- H = Σ_i h_i S^z_i + Σ_i J_i S^z_i S^z_{i+1} + Σ_i γ_i S^x_i on a finite interval.
- h_i and J_i are i.i.d. with bounded densities; γ_i = γΓ_i with Γ_i ∈ [−1, 1]; γ is small.
- There is **no U(1) symmetry**.

**Method** [snippet]
- A KAM-style sequence of quasi-local unitary rotations (multi-scale Jacobi) diagonalizes H.
- It connects every exact eigenstate to a basis product state, which gives a full set of LIOMs.

**Statement** [memory; check Thm 1.1–1.2 of the JSP paper]
- Assume LLA(ν, C) and take γ small enough.
- The diagonalizing rotation is quasi-local. Near any given site, the rotation is close to the identity with probability ≥ 1 − γ^κ; resonant spots have density ~ γ^κ.
- Consequences:
  - most eigenstates retain local S^z memory;
  - a many-body eigenfunction correlator decays exponentially;
  - hence absence of transport and local memory **for all times**.
- Bounds are uniform in chain length, hold in probability or expectation per site, and are **not almost-sure** statements about an infinite chain.

**The LLA assumption** [snippet + memory]
- In words (PRL): "with high probability the minimum level spacing should be no smaller than some exponential in the volume."
- My recollection of the formal version: for every box of n sites and every δ > 0, P(min_{α≠β}|E_α − E_β| < δ) ≤ δ^ν C^n.
- Imbrie states that "we lack the tools to adequately deal with such questions of level statistics."
- **Status: unproven.** I found no proof or partial proof for any lattice chain.
  - The nearest level-statistics theorem in a genuine many-body model is De Roeck–Hannani 2025 (B9, Quantum Sun model), which is not a short-range chain.
  - Elgart–Klein (B4) state that whether the XXZ spectrum is level-spaced with high probability is unknown [snippet].

**Scope**
- Only 1D with short-range (nearest-neighbour) interactions, and strong disorder in the sense of γ ≪ 1.
- It requires continuous distributions; binary disorder is not covered [memory].
- As far as I know, no published version exists for the U(1)-symmetric random-field XXZ chain [memory].

### B2. De Roeck–Imbrie 2017 and how avalanches interact with Imbrie's construction (H)

**References**
- W. De Roeck, J.Z. Imbrie, *Many-body localization: stability and instability*, Phil. Trans. R. Soc. A 375, 20160422 (2017), arXiv:1705.00756 [snippet].
- W. De Roeck, F. Huveneers, *Stability and instability towards delocalization in MBL systems*, PRB 95, 155129 (2017), arXiv:1608.01815 [snippet].
- Thiery, Huveneers, Müller, De Roeck, *Many-body delocalization as a quantum avalanche*, PRL 121, 140601 (2018), arXiv:1706.09338 [memory; id seen in search].

**Content**
- In 1D, the existence of MBL is "clarified" under the eigenvalue-statistics assumption, i.e. via Imbrie [snippet].
- In d > 1, "ergodicity is restored … although equilibration should be extremely slow, similar to glasses" [snippet]. This is heuristic.
- De Roeck–Huveneers: the "narrowing of the structure factor eventually brings the avalanche to a halt in 1d." Instability is predicted for d ≥ 2 and for interactions decaying sub-exponentially [snippet]. (H)

**Interaction with Imbrie's proof** (H; my synthesis)
- Imbrie isolates resonant blocks scale by scale and needs them to stay rare and small. LLA keeps the effective gaps inside a block from shrinking super-exponentially.
- Avalanche theory says a thermal block of size ℓ (level spacing ~2^{−ℓ}) thermalizes its neighbourhood if couplings decay more slowly than that spacing. In the l-bit language the condition is ξ > 1/ln 2, up to conventions.
- In Imbrie's regime (γ ≪ 1) the effective couplings decay much faster, so the avalanche stops. There is therefore no contradiction.
- Numerical estimates of the avalanche threshold for the random-field Heisenberg chain are W ≳ 20 (Morningstar et al., PRB 105, 174205 (2022), arXiv:2107.05642; N, H) [memory].
- So Imbrie's theorem says nothing about the W ~ 4–10 regime studied numerically. There is no conflict with the claims of slow drift (Šuntajs et al. PRE 2020; Sels–Polkovnikov 2021; both id?) [memory].

### B3. Quasi-free models (P, unconditional)

**Hamza, Sims, Stolz, *Dynamical localization in disordered quantum spin systems*, CMP 315, 215 (2012), arXiv:1108.3811**
- Abstract [snippet]: dynamical localization is defined as a zero-velocity LR bound. For short-range systems it implies exponential decay of ground-state correlations, up to an explicit correction. For random XY chains it reduces to dynamical localization of the effective one-particle Hamiltonian (via Jordan–Wigner). The isotropic XY chain in a random field is dynamically localized.
- Precise form [memory]: E sup_t ‖[τ_t(A), B]‖ ≤ C‖A‖‖B‖ Σ e^{−η|j−k|}.
- Holds for **all times**, uniformly in chain length, **in expectation**.
- For the isotropic chain it holds at any disorder strength, since 1D Anderson localization covers the whole spectrum. This needs regular, e.g. bounded-density, field distributions [memory].
- The anisotropic case needs localization of a block Jacobi operator: Chapman–Stolz, AHP 16 (2015), id? [memory].

**Abdul-Rahman, Nachtergaele, Sims, Stolz, *Entanglement dynamics of disordered quantum XY chains*, Lett. Math. Phys. 106, 649 (2016), arXiv:1510.00262 [snippet]**
- For a broad class of product initial states, bipartite entanglement stays **bounded for all times**. The assumption is dynamical localization of the effective one-particle operator [memory]; the bound is in expectation and uniform in system size.
- For the XX chain it also bounds **particle transport** [snippet].
- This is the rigorous statement behind the claim that "number entropy is bounded for free fermions" (S_N ≤ S_ent).

**Abdul-Rahman, Stolz, *A uniform area law for the entanglement of eigenstates in the disordered XY chain*, J. Math. Phys. 56, 121901 (2015), arXiv:1505.02117**
- Area law in expectation for **all eigenstates**, uniform in system size and, I believe, with the supremum over eigenstates inside the expectation [snippet; the sup placement is memory].
- Assumptions: almost-surely bounded coefficients and localization of the one-particle operator.
- Covers the isotropic chain in a random field, and the anisotropic chain in strong or strongly disordered fields [snippet].

**Pastur, Slavin, *Area law scaling for the entropy of disordered quasifree fermions*, PRL 113, 150404 (2014), arXiv:1408.2570**
- Free fermions in the Anderson model, Fermi sea (ground state).
- Area law in mean for d-dimensional cubes. In d = 1 the area law also holds for typical realizations, but the entropy is not self-averaging [snippet; journal from memory].

**Elgart, Pastur, Shcherbina, *Large block properties of the entanglement entropy of free disordered fermions*, J. Stat. Phys. 166, 1092 (2017), arXiv:1601.00294** [snippet; journal from memory]
- d = 1: almost-sure asymptotic form of the entanglement entropy.
- d ≥ 2 with i.i.d. potential: self-averaging.
- Area law follows from decay of the expected kernel of the Fermi projection.
- Related: Müller–Pastur–Schulte, arXiv:1903.05057, shows area-law violation at critical energies [snippet].

**Other quasi-free results** [memory]
- Bravyi–König, CMP 316 (2012), id?: zero-velocity LR for the disordered Kitaev/Majorana chain.
- Nachtergaele–Sims–Stolz, J. Stat. Phys. 149 (2012), id?: disordered harmonic lattices.
- Sims–Warzel, CMP 347 (2016), id?: decay of determinantal and Pfaffian correlations.
- Gebert–Lemm, J. Stat. Phys. 164 (2016), id?: polynomial LR in a decaying random field.
- Quasi-periodic XY chains: Damanik–Lemm–Lukic–Yessen, PRL 113, 127202 (2014), id?; Kachkovskiy, CMP 345 (2016), id?. Anomalous LR bounds and transport.

**Reviews**
- Abdul-Rahman, Nachtergaele, Sims, Stolz, Ann. Phys. (Berlin) 529, 1600280 (2017), arXiv:1610.01939 [title seen].
- Stolz, *Aspects of the mathematical theory of disordered quantum spin chains*, arXiv:1810.05047 [title seen].

### B4. Interacting models at zero energy density: random XXZ, droplets, clusters (P, unconditional)

**Common setting** [memory]
- H = Σ[(1/4 − S^zS^z) − Δ^{−1}(S^xS^x + S^yS^y)] + λΣ ω_i N_i.
- Δ > 1 (Ising phase, attractive hard-core bosons); ω_i i.i.d.
- Energy windows are **fixed while L → ∞**, i.e. zero energy density: states with a bounded number of clusters but **arbitrary particle number N**.

**Elgart, Klein, Stolz, *Many-body localization in the droplet spectrum of the random XXZ quantum spin chain*, J. Funct. Anal. 275 (2018), id?** [memory]
- Applies to the droplet band, below the two-cluster threshold.
- Exponential localization of the droplet spectral subspace; for observables or states restricted to it:
  - a zero-velocity LR bound, uniform in time;
  - dynamical exponential clustering.
- Holds for large Δ or large disorder, uniformly in L and N, in expectation.

**Related EKS papers** [snippet: abstracts mention "dynamical exponential clustering, non-spreading of information, zero-velocity LR bound"]
- Elgart, Klein, Stolz, *Droplet localization in the random XXZ model and its manifestations*, arXiv:1707.09665; J. Phys. A (2018) [journal from memory].
- Elgart, Klein, Stolz, *Manifestations of dynamical localization in the disordered XXZ spin chain*, arXiv:1708.00474. Journal: I recall CMP 361 (2018), while the task brief lists JFA. Verify.

**Beaud, Warzel** [memory; ids?]
- *Low-energy Fock-space localization for attractive hard-core particles in disorder*, AHP 18 (2017): Fock-space (configuration-space) localization of low-energy eigenstates, uniform in particle number.
- *Bounds on the entanglement entropy of droplet states in the XXZ spin chain*, JMP 59 (2018): area-law-type bounds, possibly with log corrections, for droplet-band states. I believe this also covers the non-random case.

**Abdul-Rahman, Fischbacher, Stolz, *Entanglement bounds in the XXZ quantum spin chain*, AHP 21 (2020), id?** [memory]
- Log-corrected area law for eigenstates in fixed energy intervals at the bottom of the spectrum (Ising phase). Improved bounds with disorder.

**Elgart, Klein, *Localization in the random XXZ quantum spin chain*, Forum Math. Sigma (2024), arXiv:2210.14873 (v3 Jan 2024)** [snippet]
- Localization "in any given energy interval at the bottom of the spectrum."
- The parameter region "includes weak interaction and strong disorder regimes, is independent of the size of the system and depends only on the energy interval."
- Method: the **fractional moment method applied directly to the many-body H**, with localization reformulated as quasi-locality of functions of H.
- v1 was phrased as localization in the "multi-cluster spectrum."

**Elgart, Klein, *Slow propagation of information on the random XXZ quantum spin chain*, CMP 405, 239 (2024), arXiv:2311.14188** [snippet]
- Quasi-locality in a fixed energy window implies slow propagation of information in that window.

**Elgart, Klein, *Localization phenomena in the random XXZ spin chain*, J. Funct. Anal. (2025/26; ScienceDirect S0022123625005026), arXiv:2502.01831 (v3 Nov 2025)** [snippet]
- **Infinite chain.** Spectral, eigenstate and weak dynamical localization in an arbitrary fixed window [0, E₀], uniform in N.
- One excerpt gives the condition as "λΔ² sufficiently large" (unverified).
- The authors state their methods do **not** give LIOM-type localization, and that it is unknown whether the spectrum is level-spaced with high probability. This is an LLA-type gap.

**Elgart, Klein, *Many-body localization for the random XXZ spin chain in fixed energy intervals*, arXiv:2602.01441 (Feb 2026, v3)** [snippet]
- **Logarithmic light cone** in any fixed energy interval at the bottom of the spectrum. The regime is set by the interval.

**Slow-propagation and LR results**
- Nachtergaele, Reschke, *Slow propagation in some disordered quantum spin chains*, J. Stat. Phys. (2021), arXiv:1906.10167 [snippet]:
  - (P) a zero-velocity LR bound implies a LIOM representation;
  - with sparse perturbations, transmission times diverge at least as a power of distance.
- Gebert, Moon, Nachtergaele, *A Lieb–Robinson bound for quantum spin chains with strong on-site impurities*, arXiv:2104.00968; Rev. Math. Phys. (2022?) [title seen; journal from memory].
- I could not identify a Gebert–Nachtergaele–Reschke–Sims paper on MBL. The paper I recall with those authors is *Lieb–Robinson bounds and strongly continuous dynamics for a class of many-body fermion systems in R^d*, AHP 21 (2020), id? [memory]. It is not about disorder.

**Mavi, Schenker, *Localization in the disordered Holstein model*, CMP 364 (2018)** [journal from memory], arXiv:1709.06621 [snippet]
- One electron coupled to oscillators, low energies, arbitrary dimension: localization bounds in Fock space, with no particle-number conservation for the phonons, at zero energy density.
- Follow-up: Mavi, *Log-dimensional bounds for the spectral measure of the disordered Holstein model*, arXiv:2009.03182 [snippet], at finite energies for small γ.

**Finite particle number** [memory; ids?]
- Chulaevsky–Suhov (2009); Aizenman–Warzel, CMP 290 (2009); Klein–Nguyen (2013–15).
- For a fixed number of particles n, localization is proven, but the constants degrade with n. These results say nothing in the thermodynamic limit.

### B5. Quasi-periodic interacting fermions (P; ground state only)

**Mastropietro**
- *Localization in the ground state of an interacting quasi-periodic fermionic chain*, CMP 342, 217 (2016), arXiv:1411.0201 [snippet]:
  - "zero temperature thermodynamical correlations are exponentially decaying … for certain densities and weak interactions."
  - Method: RG with fermionic cancellations controlling small divisors. Assumes Diophantine frequency and chemical potential [memory].
- *Localization of interacting fermions in the Aubry–André model*, PRL 115, 180401 (2015), id? [memory]: the same circle of results.
- *Localization in interacting fermionic chains with quasi-random disorder*, CMP 351 (2017) [journal from memory], arXiv:1604.08264 [snippet]: almost-Mathieu potential with short-range interaction; compared with the numerical MBL region in the (U/u, ε/u) plane. These are T = 0 statements.
- *Interacting spinning fermions with quasi-random disorder*, Ann. Phys. (Berlin) 529, 1600270 (2017), arXiv:1611.00748 [snippet]: with spin, an extra relevant quartic term appears.
- *Coupled identical localized fermionic chains with quasi-random disorder*, arXiv:1701.00543 [snippet]: one or two chains stay localized; more chains are generically different.
- *Vanishing of Drude weight in interacting fermions on Z^d with quasi-periodic disorder*, arXiv:2102.07822 (2021), journal? [snippet]: the T = 0 Drude weight vanishes.

**Status**
- All of these are ground-state or T = 0 statements in infinite volume. None concerns excited states, infinite temperature or long-time dynamics.
- My searches found no 2023–26 Mastropietro MBL papers, but I cannot confirm there are none.

**Two particles:** Bourgain–Kachkovskiy, GAFA 29 (2019), id?: localization for two interacting quasi-periodic particles [memory].

### B6. Quasi-localization and very slow thermalization (P, finite times or rare events)

**De Roeck, Huveneers, Meeus, Prośniak, *Rigorous and simple results on very slow thermalization, or quasi-localization, of the disordered quantum chain*, Physica A 631, 129245 (2023), arXiv:2308.08462**
- Abstract [snippet]: "For a strongly disordered quantum chain, we exhibit a full set of quasi-local quantities whose dynamics is negligible up to times of order exp{c(log W)^{2−ε}}."
- The proof is short; the paper grew out of 2021 Bruneck lectures.

**What the theorem gives** [memory; check the statement]
- Dressed quasi-local Z̃_i, obtained from about log W orders of Schrieffer–Wolff rotations, with E‖[H, Z̃_i]‖ ≲ exp{−c(log W)^{2−ε}}.
- Per-site errors, in expectation or probability, uniform in L.
- Hence |Z̃_i(t) − Z̃_i(0)| is negligible for t ≪ t\*.
- t\* = W^{c(log W)^{1−ε}} is **superpolynomial in W but not exponential**. The optimal perturbative order is ~log W because resonances at order n have probability ~(n!/W^n)-like.
- Not an infinite-time statement. No LLA is needed.
- I did not confirm the exact model class (random fields of strength W with bounded density plus O(1) nearest-neighbour couplings is my recollection), nor whether ε is arbitrary.

**De Roeck, Giacomin, Huveneers, Prośniak, *Absence of normal heat conduction in strongly disordered interacting quantum chains*, arXiv:2408.04338 (v2 substantially rewrote the second half), journal?** [snippet]
- "The conductance of a chain of length L vanishes faster than 1/L," at **any temperature**.
- As restated by Scardicchio–Sondhi: the time-averaged heat current is bounded by a decaying power of L, so lim L·E[J] = 0.
- **Key ingredient** (from Laumann's JCCM commentary): a chain of length L is fully MBL, i.e. diagonalized by a convergent iterated Schrieffer–Wolff flow, **with probability ≥ e^{−cL}**, for weak off-diagonal terms.
- Restricting to such exponentially rare non-resonant segments removes the need for LLA. These segments act as bottlenecks in a long chain.
- (P) for heat transport. The exact disorder and model class was not verified.

**Scardicchio, Sondhi, *On rare nonresonant regions and subdiffusive transport in an interacting disordered quantum chain*, arXiv:2609.13834 (Sep 2026; v2 Nov 2026)** [snippet] (N, H)
- Implements the DGHP flow numerically on finite systems and tracks the probability that a sample stays non-resonant scale by scale.
- If the series converges, the probability is exponential in length. If it does not, the probability is super-exponentially small and the rare-region mechanism fails.
- They suggest a possible sharp transition.
- They mention "a forthcoming proof by De Roeck, Huveneers and collaborators" giving "another route to the same result." I could not identify it.

**De Roeck, Huveneers, *Asymptotic quantum many-body localization from thermal disorder*, CMP 332, 1017 (2014), arXiv:1308.6263** [id seen; content from memory] (P)
- Translation-invariant chains with a small hopping parameter g, at high temperature.
- Relaxation and transport coefficients vanish faster than any power of g ("asymptotic localization").
- Related: De Roeck–Huveneers, CPAM (2015), arXiv:1305.5127 (oscillator chains); Bols–De Roeck, JMP (2018), arXiv:1612.04731 (Bose–Hubbard) [ids seen].
- A companion PRB (2014) argues that true MBL does not exist in translation-invariant systems (H).

**Huveneers, *Many-body localization: an introduction for mathematicians*, arXiv:2609.39462 (Les Houches 2026 notes)** [snippet]
- Covers resonances, avalanches and a rigorous treatment of localization in a resonance-free region. Probably the best current entry point.

### B7. Disordered Lieb–Robinson bounds without localization (P)

- Baldwin, Ehrenberg, Guo, Gorshkov, *Disordered Lieb–Robinson bounds in one dimension*, arXiv:2208.05509, PRX Quantum (2023)? [title seen; authors and journal from memory]. Light cones for random couplings; sub-ballistic when weak links are frequent.
- C.L. Baldwin, *Sub-ballistic operator growth in spin chains with heavy-tailed random fields*, PRB 111, 184204 (2025), arXiv:2409.17242 [snippet]. For power-law-distributed fields with tail exponent α:
  - ballistic growth is impossible for α < 1;
  - diffusive growth is impossible for α < 1/2;
  - any finite dynamical exponent is impossible for small α.
- These bounds hold for all states and all times. They give operator-spreading bounds, not localization.

### B8. Prethermalization theorems (P; disorder-free; finite but exponentially long times)

**Abanin, De Roeck, Ho, Huveneers, *A rigorous theory of many-body prethermalization for periodically driven and closed quantum systems*, CMP 354, 809 (2017), arXiv:1509.05386** [snippet + memory]
- Closed case: H = U·N + Y, where N is a sum of commuting local terms with **integer spectrum** and Y is local of strength J.
- A quasi-local unitary dresses N into N′, conserved up to τ\* ~ exp[c(U/J)/log³(U/J)] ("(almost) exponential in U/J" [snippet]).
- Driven case: an effective Hamiltonian exists and heating is suppressed up to exp[cω/(J log³)].
- Errors are controlled per unit volume.

**Abanin, De Roeck, Ho, Huveneers, *Effective Hamiltonians, prethermalization and slow energy absorption in periodically driven many-body systems*, PRB 95, 014112 (2017), arXiv:1510.03405** [title seen].

**Kuwahara–Mori–Saito, Ann. Phys. 367 (2016), id?; Mori–Kuwahara–Saito, PRL 116, 120401 (2016), id?** [memory]
- Floquet–Magnus bounds; energy absorption suppressed as exp(−cω/J).

**Else, Bauer, Nayak, *Prethermal phases of matter protected by time-translation symmetry*, PRX 7, 011026 (2017), arXiv:1607.05277** [snippet; theorem from memory]
- Generalizes ADHH: in a rotating frame there is an emergent Z_n symmetry up to exp(cω/log³ω), which gives prethermal time crystals and SPT phases.

**Else, Ho, Quinn, *Long-lived interacting phases of matter protected by multiple time-translation symmetries in quasiperiodically driven systems*, PRX 10 (2020), id?** [memory]
- Stretched-exponential lifetimes for quasi-periodic drives.

**Yin, Lucas, *Prethermalization and the local robustness of gapped systems*, PRL 131 (2023), id?** [memory, not verified]
- Prethermal bounds with **local** errors, independent of system size.

**Uses relevant to the user**
- **Stark.** In the interaction picture, a tilt F maps to a Floquet problem with frequency F. The theorems then give emergent dipole conservation up to ~exp(cF/J) (P via the known theorems) [memory].
  - Nandy, Herbrych, Lenarčič, Główkowski, Prelovšek, Mierzejewski, *Emergent dipole moment conservation and subdiffusion in tilted chains*, PRB 109, 115120 (2024), arXiv:2310.01862 [snippet]: bounds dipole-moment dynamics for arbitrary initial states, with a stronger bound at infinite temperature, "formally shown" for a broad class of models including higher dimensions.
  - **I found no rigorous proof of Stark MBL at infinite times.**
- **Gauge protection.** A protection term V Σ c_j G_j with commensurate (integer) c_j fits the ADHH form, which gives gauge violation suppressed until ~exp(cV/λ). Halimeh–Hauke and coworkers (PRL 125, 030503 (2020); PRX Quantum 2, 040311 (2021); ids?) build on this [memory]. Incommensurate "linear protection" is heuristic.
- **Confinement or strong coupling.** Electric energy g²ΣE², or a large longitudinal field ΣZ_i, has integer spectrum. ADHH then gives quasi-conservation of domain-wall or string number up to exp(c g²/J). (P for that quasi-conservation; the confinement-induced localization phenomenology itself is H.)

### B9. Toy many-body models (P)

**De Roeck, Hannani, *Many-body localization and Poisson statistics in the Quantum Sun model*, arXiv:2506.13511 (2025)** [snippet]
- Spins at distance n couple to a dot with strength α^n.
- For α ≪ 1/√2 they prove **localization and Poisson spectral statistics**. One secondary summary states the condition as α ≤ α₀ < 1/2, so check.
- The authors stress that the number of disorder variables grows only like log(Hilbert-space dimension), and propose the name "Fock-space localization."
- Avalanche theory puts the transition at α_c = 1/√2.
- This is the only rigorous many-body level-statistics result I found. It is relevant to LLA.

**Yin, Nandkishore, Lucas, *Eigenstate localization in a many-body quantum system*, PRL 133, 137101 (2024), arXiv:2405.12279** [snippet]
- Extensive few-body Hamiltonians built from good classical codes, with a **many-body mobility edge**: all eigenstates below a nonzero energy density are localized.
- The geometry is effectively infinite-dimensional (expander); there is no finite-d lattice. The proof is "related to … Anderson localization."

### B10. Negative and no-go results: status

- **No rigorous delocalization or ergodicity theorem** exists for disordered interacting lattice Hamiltonians (to my knowledge).
- Heuristic or numerical only (H or N):
  - **Avalanches and d > 1:** B2. Experiment: Léonard et al., Nat. Phys. 19 (2023), id? [memory] (E).
  - **Long-range interactions:** Burin; Yao et al. PRL 2014; Gutman et al. PRB 2016; ids? [memory] (H).
  - **Continuum:** Gornyi, Mirlin, Polyakov, Burin, *Absence of many-body localization in a continuum*, Ann. Phys. (Berlin) (2017), id? [memory] (H).
  - **SU(2):** Potter–Vasseur, PRB 94, 224206 (2016), id?; Protopopov–Ho–Abanin, PRB 96 (2017), id?; Protopopov et al. PRX (2020) [memory] (H, N).
    - The algebraic core is essentially rigorous and elementary: a complete set of commuting quasi-local LIOMs with non-degenerate product structure cannot be SU(2)-invariant.
    - That rules out only "full LIOM-MBL". It does not prove thermalization.
  - **No many-body mobility edges** in short-range lattice systems: De Roeck, Huveneers, Müller, Schiulaz, PRB 93 (2016), id? [memory] (H). Yin–Nandkishore–Lucas (B9) show this does not carry over to expander geometries.
  - **Translation-invariant systems:** no MBL expected (H); asymptotic localization proven (B6).
- **Perturbative check of weak-interaction avalanches:** Jiang, Surace, Motrunich, *Quasiconservation laws and suppressed transport in weakly interacting localized models*, PRB 112, 184201 (2025), arXiv:2507.03115 [snippet] (perturbative, N). At first order, transport capacity stays bounded for Anderson and Aubry–André models; avalanche instability would need higher orders or non-perturbative effects.

---

## (c) Technical obstacles, and intermediate statements within reach

### Obstacles

1. **Level statistics with few random variables (LLA).**
   - There are about 2^n levels but only n random variables.
   - Wegner-type density bounds are easy. Two-level (Minami-type) spacing bounds are not: for pairs of eigenstates inside a resonant block, the derivative of the energy difference with respect to local fields is not controlled.
   - Elgart–Klein hit the same wall, which blocks LIOM localization.
   - Progress so far: Poisson statistics for the Quantum Sun model (B9). Nothing for lattice chains.
2. **Resonances at all scales and non-perturbative resonant spots.**
   - Imbrie handles them only with LLA.
   - DHMP stop at the optimal perturbative order, which gives the finite t\*.
   - DGHP condition on exponentially rare fully non-resonant segments.
   - Scardicchio–Sondhi (N) question whether even that probability stays exponential at moderate W.
3. **Avalanches.** Controlling infinite-time dynamics needs a proof that thermal spots do not grow. That reduces again to level statistics or eigenstate-thermalization properties inside resonant blocks, which no one can prove.
4. **Positive energy density.** All unconditional infinite-time interacting results are at zero density (B4), on rare events (DGHP), or in toy geometries (B9).
5. **Discrete or correlated disorder.**
   - Imbrie and DHMP-type probability estimates use bounded densities: P(|h_i − h_j| < δ) ≤ Cδ.
   - Binary disorder (gauge sectors) and quasi-periodic potentials (no independence, Diophantine arithmetic) need different resonance control.
   - Mastropietro handles quasi-periodic potentials only at T = 0.
6. **U(1) symmetry** (S^z-conserving hopping). As far as I know, Imbrie's construction has not been published for random-field XXZ; it is a technical issue, not an obstruction of principle [memory]. EK and DHMP-type arguments do handle U(1) models or are compatible with them; this needs checking for DHMP.

### (i) Particle-number fluctuations and number entropy in U(1)-symmetric disordered chains

**Proven**
- **Quasi-free XX in a random field:** particle transport and entanglement are bounded uniformly in time and L, in expectation (ANSS 2016). Since S_N ≤ S_ent, the number entropy is bounded for all times.
- **Zero energy density** (EK 2024–26): quasi-locality, slow propagation and a log light cone for states in fixed windows. A bound on number fluctuations for such states is plausibly a short corollary (my assessment).

**Not proven:** any bound on Var N_A(t) or S_N(t) at positive energy density for t → ∞, in any interacting model. Numerical claims of unbounded S_N growth (Kiefer-Emmanouilidis et al., PRL 124 (2020), id?; Luitz–Bar Lev, PRB 102 (2020), id?; Ghosh–Žnidarič, PRB 105 (2022), id? [memory]) concern W ~ 4–10. That regime is far outside every theorem, so no rigorous statement conflicts with them.

**Within reach** (my assessment, H, flagged as such)
- **(a) Corollary of the DHMP construction for random-field XXZ.**
  - Take the U(1)-preserving dressing U, restricted near a cut. Then N_A − U N_A U† is an operator localized at the boundary, and [U†HU, N_A] involves only remainder terms that straddle the cut.
  - This should give sup_{|t|≤t\*} ‖N_A(t) − N_A(0)‖ ≤ B + t·ε(W). Here B is a random boundary constant with good tails, ε(W) ~ exp{−c(log W)^{2−ε}}, and the bound holds uniformly in L and for every initial state.
  - Then Var_{ψ(t)} N_A ≤ (B + tε + √Var_ψ N_A)².
  - With H(X) ≤ ½ log[2πe(σ² + 1/12)] for integer-valued X, this gives **S_N(t) bounded uniformly in L up to t\***.
  - Caveats: (1) check that DHMP's model class allows U(1)-symmetric hopping and that the dressing preserves N; (2) the constants are likely meaningful only at very large W.
- **(b) A spin- or charge-conductance version of DGHP:** σ_L = o(1/L) for XXZ. This gives subdiffusion bounds, not bounded S_N.
- **(c) A conditional theorem (P-level, short):** if U(1)-symmetric quasi-local LIOMs exist (Imbrie-type, with probabilistic tails), then sup_t ‖N_A(t) − N_A(0)‖ ≤ 2‖N_A − UN_AU†‖. That is O(1) in probability, so S_N is bounded for all times.
  - So any proof of unbounded S_N growth would also disprove the LIOM picture at that W.
- **(d) Bounds near rare non-resonant segments** (DGHP-type): transported charge stays bounded up to times ~e^{cℓ} across such a segment of length ℓ. This is not typical-cut information.

### (ii) Lattice gauge theories and disorder-free localization (DFL)

**Proven ingredients**
- Z₂ gauge field coupled to free fermions (Smith, Knolle, Kovrizhin, Moessner, PRL 118, 266601 (2017), id? [memory]): sector by sector this is the **Bernoulli Anderson model**.
  - In 1D it is localized at all energies with dynamical localization (Carmona–Klein–Martinelli, CMP 108 (1987); later proofs by Jitomirskaya–Zhu and by Bucaj et al.; ids? [memory]).
  - So for i.i.d. background-charge distributions (e.g. uniform superpositions), HSS/ANSS-type zero-velocity LR and number-fluctuation bounds should transfer. That is essentially available but not, to my knowledge, written down.
  - For correlated sector weights, positivity of the Lyapunov exponent for the induced potential is needed.
- In 2D, Bernoulli localization is proven only near band edges (Ding–Smart, Invent. Math. (2020), id? [memory]). Bulk 2D DFL is open.

**Interacting DFL** (Brenes, Dalmonte, Heyl, Scardicchio, PRL 120 (2018); Karpov et al., PRL 126 (2021); Halimeh et al.; ids? [memory]; N)
- **Nothing rigorous beyond prethermal statements.**
- The obstruction: with discrete background charges, resonances at fixed perturbative order have **finite probability**, so the DHMP-type smallness (resonance probability ≲ (J/W)^n) is lost.
- The accessible rigorous statements are ADHH-type:
  - when the gauge protection, string tension or electric energy is large with integer spectrum, gauge violation or domain-wall number is quasi-conserved up to exp(c·V/J), non-randomly;
  - with local errors à la Yin–Lucas (if verified).
- Adding a continuous random coupling would let DHMP-type quasi-localization apply. That is within reach but changes the model.

---

## (d) References (tag shows what was checked)

**Imbrie, LLA, avalanches**
- J.Z. Imbrie, J. Stat. Phys. 163, 998 (2016), arXiv:1403.7837 [snippet]
- J.Z. Imbrie, PRL 117, 027201 (2016), arXiv:1605.03003 [snippet]
- J.Z. Imbrie, *Multi-scale Jacobi method for Anderson localization*, CMP 341 (2016), id? [memory]
- Imbrie, Ros, Scardicchio, *Local integrals of motion in many-body localized systems*, Ann. Phys. (Berlin) 529 (2017), id? [memory]
- De Roeck, Imbrie, Phil. Trans. R. Soc. A 375, 20160422 (2017), arXiv:1705.00756 [snippet]
- De Roeck, Huveneers, PRB 95, 155129 (2017), arXiv:1608.01815 [snippet]
- Thiery, Huveneers, Müller, De Roeck, PRL 121, 140601 (2018), arXiv:1706.09338 [memory]
- Morningstar, Colmenarez, Khemani, Luitz, Huse, PRB 105, 174205 (2022), arXiv:2107.05642 [memory]
- Potirniche et al., *Exploration of the stability of MBL in d > 1*, PRB (2019), arXiv:1805.01475 [snippet]
- Sierant, Lewenstein, Scardicchio, Vidmar, Zakrzewski, *MBL in the age of classical computing*, Rep. Prog. Phys. 88 (2025), id? (I believe 2403.07111) [memory]

**Quasi-free**
- Hamza, Sims, Stolz, CMP 315, 215 (2012), arXiv:1108.3811 [snippet]
- Abdul-Rahman, Nachtergaele, Sims, Stolz, LMP 106, 649 (2016), arXiv:1510.00262 [snippet]
- Abdul-Rahman, Stolz, JMP 56, 121901 (2015), arXiv:1505.02117 [snippet]
- Pastur, Slavin, PRL 113, 150404 (2014), arXiv:1408.2570 [snippet]
- Elgart, Pastur, Shcherbina, JSP 166 (2017), arXiv:1601.00294 [snippet]
- Müller, Pastur, Schulte, arXiv:1903.05057 [snippet]
- Abdul-Rahman, Nachtergaele, Sims, Stolz (review), Ann. Phys. (Berlin) 529 (2017), arXiv:1610.01939 [title seen]
- Stolz, *Aspects of the mathematical theory of disordered quantum spin chains*, arXiv:1810.05047 [title seen]

**Interacting, zero energy density, XXZ, Holstein**
- Elgart, Klein, Stolz, JFA 275 (2018), id? [memory]
- Elgart, Klein, Stolz, arXiv:1707.09665 [snippet]
- Elgart, Klein, Stolz, arXiv:1708.00474; CMP 361 (2018)? [snippet / memory]
- Beaud, Warzel, AHP 18 (2017), id?; JMP 59 (2018), id? [memory]
- Abdul-Rahman, Fischbacher, Stolz, AHP 21 (2020), id? [memory]
- Elgart, Klein, Forum Math. Sigma (2024), arXiv:2210.14873 [snippet]
- Elgart, Klein, CMP 405, 239 (2024), arXiv:2311.14188 [snippet]
- Elgart, Klein, JFA (2025/26), arXiv:2502.01831 [snippet]
- Elgart, Klein, arXiv:2602.01441 (2026) [snippet]
- Nachtergaele, Reschke, JSP (2021), arXiv:1906.10167 [snippet]
- Gebert, Moon, Nachtergaele, arXiv:2104.00968 [title seen]
- Mavi, Schenker, CMP 364 (2018)?, arXiv:1709.06621 [snippet]
- Mavi, arXiv:2009.03182 [snippet]

**Quasi-periodic**
- Mastropietro, CMP 342, 217 (2016), arXiv:1411.0201 [snippet]
- Mastropietro, PRL 115, 180401 (2015), id? [memory]
- Mastropietro, arXiv:1604.08264 (CMP 351, 2017?) [snippet]
- Mastropietro, arXiv:1611.00748 [snippet]
- Mastropietro, arXiv:1701.00543 [snippet]
- Mastropietro, arXiv:2102.07822 [snippet]
- Bourgain, Kachkovskiy, GAFA 29 (2019), id? [memory]

**Quasi-localization, rare regions, toy models, LR**
- De Roeck, Huveneers, Meeus, Prośniak, Physica A 631, 129245 (2023), arXiv:2308.08462 [snippet]
- De Roeck, Giacomin, Huveneers, Prośniak, arXiv:2408.04338 [snippet]
- Scardicchio, Sondhi, arXiv:2609.13834 (2026) [snippet]
- Huveneers, arXiv:2609.39462 (2026) [snippet]
- De Roeck, Huveneers, CMP 332, 1017 (2014), arXiv:1308.6263 [memory]
- De Roeck, Huveneers, CPAM (2015), arXiv:1305.5127 [memory]
- Bols, De Roeck, JMP (2018), arXiv:1612.04731 [memory]
- De Roeck, Hannani, arXiv:2506.13511 (2025) [snippet]
- Yin, Nandkishore, Lucas, PRL 133, 137101 (2024), arXiv:2405.12279 [snippet]
- Baldwin, PRB 111, 184204 (2025), arXiv:2409.17242 [snippet]
- Baldwin, Ehrenberg, Guo, Gorshkov, arXiv:2208.05509 [memory]
- Jiang, Surace, Motrunich, PRB 112, 184201 (2025), arXiv:2507.03115 [snippet]

**Prethermalization and Stark**
- Abanin, De Roeck, Ho, Huveneers, CMP 354, 809 (2017), arXiv:1509.05386 [snippet]
- Abanin, De Roeck, Ho, Huveneers, PRB 95, 014112 (2017), arXiv:1510.03405 [title seen]
- Else, Bauer, Nayak, PRX 7, 011026 (2017), arXiv:1607.05277 [snippet]
- Kuwahara, Mori, Saito, Ann. Phys. 367 (2016), id? [memory]
- Mori, Kuwahara, Saito, PRL 116, 120401 (2016), id? [memory]
- Else, Ho, Quinn, PRX (2020), id? [memory]
- Yin, Lucas, PRL 131 (2023), id? [memory]
- Nandy et al., PRB 109, 115120 (2024), arXiv:2310.01862 [snippet]
- Halimeh, Hauke, PRL 125, 030503 (2020), id? [memory]

**Disorder-free localization and number entropy**
- Smith, Knolle, Kovrizhin, Moessner, PRL 118, 266601 (2017), id? [memory]
- Brenes et al., PRL 120 (2018), id? [memory]
- Ding, Smart, Invent. Math. (2020), id? [memory]
- Carmona, Klein, Martinelli, CMP 108 (1987) [memory]
- Kiefer-Emmanouilidis et al., PRL 124 (2020), id? [memory]
- Luitz, Bar Lev, PRB 102 (2020), id? [memory]
- Ghosh, Žnidarič, PRB 105 (2022), id? [memory]

**Not verified at all:** Gebert–Nachtergaele–Reschke–Sims (no MBL paper identified); Mastropietro 2023–26; the full statements of Yin–Lucas and Else–Ho–Quinn; any LLA progress (searched, nothing found).
