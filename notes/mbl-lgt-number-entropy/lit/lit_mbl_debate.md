# Is MBL a phase in the thermodynamic limit? Literature map (as of 2026-11-22)

**Verification status (read this first).** I could not read any source in this session. The shared web-search budget was already used up when I made my first query, and arxiv.org and scipost.org failed DNS resolution in this sandbox (they would have to be allowed in the session's network settings). So every claim below is tagged **[memory]** and none is [read]. Journal references and arXiv ids are from memory. Where I am not confident of an id I wrote "id?". Numbers (W_c estimates, system sizes, times) are also from memory and are approximate, so check them before citing. Coverage of 2025-26 is thin for the same reason.

**Conventions.** H = Σ_i [S^x_i S^x_{i+1} + S^y_i S^y_{i+1} + Δ S^z_i S^z_{i+1}] + Σ_i h_i S^z_i, with spin-1/2 operators, Δ = 1, and h_i uniform in [−W, W]. Papers that use Pauli matrices rescale W by a factor of 2. Evidence types: (P) proof, (N) numerics, (H) heuristic or phenomenology, (E) experiment.

---

## (a) Bottom line

As of 2026 the two camps agree on more than the 2019–21 exchanges suggested [memory].

1. **The ED "transition" is not the thermodynamic-limit transition.** ED places it at W ≈ 3.7 for L ≤ 22–26. Every finite-size estimator drifts to larger W as L grows, and in this window the crossover is controlled by rare many-body resonances.
2. **A broad slow regime exists.** Morningstar et al. call it the "MBL regime" and Long et al. call it the "prethermal MBL regime". In it, relaxation times grow at least exponentially with W, and ℓ-bit (LIOM) phenomenology describes every timescale reachable by numerics and experiment. No experiment can tell this regime apart from a true phase.
3. **Avalanches decide the asymptotic fate.** Both sides accept the avalanche mechanism (thermal inclusions spreading through the localized chain). Avalanche-based numerics push any 1D transition to W_c of order 10–20 (Morningstar et al. quote W_c ≳ 20), far beyond ED.
4. **Some settings are believed to have no asymptotic MBL.** This is mainly (H): d ≥ 2, slowly decaying power-law couplings, non-abelian symmetry, and the continuum with unbounded energies.

**What is still disputed** is the 1D short-range random chain at very strong disorder: does the relaxation time diverge at some finite W_c (a phase), or stay finite at every W (a crossover only)?
- **The majority view is that a phase exists.** It rests on three things: Imbrie's proof, which assumes "limited level attraction" (LLA, an unproven bound on how often many-body levels nearly coincide); all-orders perturbation theory; and the internal consistency of the avalanche theory and its Kosterlitz-Thouless (KT) renormalization-group picture.
- **The skeptics** include Šuntajs–Bonča–Prosen–Vidmar, Sels–Polkovnikov, Evers–Bera, and Prelovšek–Mierzejewski.
- **Nothing settles it.** There is no unconditional proof, and no numerics or experiment reaches the asymptotic regime at W ≳ 20.
- **The two strongest claims are not equally supported.** "No MBL at any W" rests on extrapolating W*(L) ∝ L from L ≤ 20–24. The defensible skeptic claim is narrower: there is no evidence of a phase for W ≲ 15–20.

The 2025 Rep. Prog. Phys. review was co-written by Sierant, Lewenstein, Scardicchio, Vidmar and Zakrzewski, so it spans both camps. From memory, it frames existence of the phase as an open question and adopts the regime-versus-phase distinction. That suggests convergence on how the question is posed [memory].

---

## (b) Who claims what

### 1. Standard picture (brief)

**Perturbative foundations**
- **Basko–Aleiner–Altshuler 2006** (H; perturbative to all orders): weakly interacting electrons in an Anderson insulator. There is a finite-T metal–insulator transition. Conductivity σ = 0 exactly for T < T_c ~ δ_ξ/(λ|ln λ|), where δ_ξ is the level spacing within a localization volume and λ the dimensionless interaction [memory].
- **Gornyi–Mirlin–Polyakov 2005** (H): the same conclusion, including in 1D [memory].

**Early numerics and the ℓ-bit picture**
- **Oganesyan–Huse 2007** (N; ED, L ≲ 16): at infinite temperature, the level statistics of a disordered interacting chain cross over from GOE to Poisson [memory].
- **Pal–Huse 2010** (N; ED, L ≤ 16; random-field Heisenberg chain): eigenstate-based diagnostics give W_c ≈ 3.5 ± 1 [memory].
- **Žnidarič–Prosen–Prelovšek 2008** (N; tDMRG): entanglement grows logarithmically in the random-field XXZ chain [memory].
- **Bardarson–Pollmann–Moore 2012** (N; TEBD): S(t) ~ ln t grows without bound and is absent at Δ = 0 [memory].
- **Serbyn–Papić–Abanin 2013 and Huse–Nandkishore–Oganesyan 2014** (H): ℓ-bits with H = Σ h̃_i τ^z_i + Σ J_ij τ^z_i τ^z_j + …, where J_ij ~ e^{−|i−j|/ξ}. Dephasing through these couplings explains the ln t growth. Ros–Müller–Scardicchio 2015 constructed the LIOMs order by order (H) [memory].

**Precision ED**
- **Luitz–Laflorencie–Alet 2015** (N; shift-invert ED, L ≤ 22): W_c ≈ 3.7 at infinite T, plus a many-body mobility edge.
- The fitted ν ≈ 1 violates the bound ν ≥ 2/d derived by Chandran–Laumann–Oganesyan (a Harris/CCFS-type bound). Khemani et al. (PRX 2017) also argued the sizes are not asymptotic. This was the first quantitative sign that the sizes are pre-asymptotic [memory].
- Later ED reached L = 24–26 (Pietracaprina et al. 2018, shift-invert; Sierant et al. 2020, POLFED, a polynomial-filtering ED method) [memory].

**Proof**
- **Imbrie 2016** (P, conditional):
  - Model: the 1D Ising chain Σ h_i σ^z_i + γ_i σ^x_i + J_i σ^z_i σ^z_{i+1}, with random h_i and J_i and weak transverse fields γ.
  - Result: almost surely there is a quasi-local unitary that diagonalizes H, so ℓ-bits exist with exponential tails.
  - Assumption: LLA, a probability bound on near-degeneracies of many-body levels in finite blocks. It is unproven.
  - Scope: this is a statement at strong disorder only, with no explicit W_c [memory].

**Reviews:** Nandkishore–Huse 2015; Alet–Laflorencie 2018; Abanin–Altman–Bloch–Serbyn 2019 (RMP) [memory].

**Experiments** (E):
- Cold atoms: Schreiber et al. 2015 (1D quasiperiodic fermions, imbalance); Choi et al. 2016 (2D bosons); Bordia et al. 2017 (2D quasiperiodic); Lukin et al. 2019 (number vs configurational entanglement); Rispoli et al. 2019 (critical regime); Rubio-Abadal et al. 2019 (coupling to a bath).
- Trapped ions: Smith et al. 2016 (L = 10).
- Superconducting qubits: Roushan et al. 2017; Guo et al. 2021.
- All of these reach at most ~10² to 10³ tunnelling times and a few tens of sites in 1D, so they cannot tell a phase from a regime [memory].

### 2. The challenge and the responses

**The original challenge**
- **Šuntajs–Bonča–Prosen–Vidmar, PRE 2020** (N; full ED, L ≲ 20; spectral form factor):
  - The Thouless time scales as t_Th ∝ L² e^{W/Ω}, i.e. diffusive with an exponentially large prefactor.
  - Setting t_Th equal to the Heisenberg time t_H ~ 2^L gives a crossover W* ∝ L, which drifts linearly.
  - Extrapolated, this means no MBL in the thermodynamic limit. The gap ratio ⟨r⟩ also drifts [memory].
- **Šuntajs et al., PRB 2020** (N): a BKT-like scaling ansatz for ⟨r⟩ fits the ED data, with W*(L) still drifting [memory].

**Responses defending a phase**
- **Abanin et al., Ann. Phys. 2021** (H, plus reanalysis):
  - The linear drift cannot continue. At strong disorder, relaxing over length L needs L-th-order processes, roughly t ~ (W/J)^{cL}, which is incompatible with L² e^{W/Ω}.
  - Drifts of this kind are expected generically from rare resonances; Anderson models on random graphs show analogous behaviour.
  - The data are compatible with MBL at larger W [memory].
- **Panda–Scardicchio–Schulz–Taylor–Žnidarič, EPL 2019/20** (H/N): ED-accessible sizes cannot locate the transition, and the drift is compatible with both scenarios [memory].
- **Sierant–Lewenstein–Zakrzewski, PRL 2020** (N; POLFED, L ≤ 24): the drift slows at larger L. Fitting W*(L) = W_∞ − a/L gives W_∞ ≈ 5.4 (number from memory; verify).
- **Sierant–Delande–Zakrzewski, PRL 2020** (N): Thouless-time analysis calibrated against Anderson models; also argues for a finite W_c [memory].

**Dynamical and spectral-function challenges**
- **Sels–Polkovnikov, PRE 2021, "dynamical obstruction"** (N; ED, L ≲ 20; plus H):
  - Tools: the norm of the adiabatic gauge potential (AGP; equivalently a fidelity susceptibility) and the low-frequency spectral functions of local operators.
  - Finding: the relaxation time grows roughly exponentially with W but stays finite at every W studied. The ED crossover is where the relaxation time reaches the Heisenberg time.
  - Conclusion: a transition, if it exists, lies far beyond the ED estimates [memory].
- **Sels, PRB 2022, bath-induced delocalization** (N/H): couples the chain to a small bath and applies a golden-rule analysis using exact spectral functions. The result is an avalanche-type instability up to W of order 20 (exact figure to check) [memory].
- **Sels–Polkovnikov, PRX 2023**: dilute impurities thermalize in 1D chains; this is related [memory].
- **Evers–Bera group** (N; Chebyshev and ED dynamics, L ≲ 20–24, t up to ~10³–10⁴):
  - Bera et al. PRL 2017: density propagator, transient subdiffusion, strong finite-size effects.
  - Weiner–Evers–Bera PRB 2019: imbalance keeps decaying slowly at W > W_c^ED, saturation values drift with L, for both random and quasiperiodic potentials.
  - Evers–Modak–Bera ("Internal clock of many-body delocalization", PRL ~2023, id?): imbalance curves collapse versus a logarithmic "internal clock", with no sign of freezing [memory].
- **Prelovšek, Mierzejewski and coworkers** (N):
  - Dynamical conductivity and density correlations relax slowly and subdiffusively, with no sharp transition (Mierzejewski–Herbrych–Prelovšek 2016; Prelovšek et al. Ann. Phys. 2017; Herbrych–Mierzejewski–Prelovšek 2022).
  - The Hubbard chain with SU(2) symmetry has no full MBL (Prelovšek–Barišić–Žnidarič 2016).
  - Vidmar–Krajewski–Bonča–Mierzejewski PRL 2021: spectral functions of local operators at infinite T show a broad ~1/ω regime.
  - Krajewski et al. PRL 2022, "Restoring ergodicity…": argue that relaxation times grow only about exponentially in W (details uncertain) [memory].
- **Number-entropy debate** (N):
  - Kiefer-Emmanouilidis–Unanyan–Fleischhauer–Sirker, PRL 2020 (ED, L ≲ 24, very long times): S_N ~ ln ln t, which they read as slow particle transport in the MBL regime.
  - Luitz–Bar Lev, PRB 2020: no slow particle transport; the growth is a finite-size effect.
  - Ghosh–Žnidarič, PRB 2022: the growth comes from rare resonances.
  - Unresolved, but mostly read as a resonance or finite-size effect [memory].

### 3. "MBL regime" versus "MBL phase"

- **Morningstar–Colmenarez–Khemani–Luitz–Huse, PRB 2022** (N plus H). Two analyses:
  - Avalanche analysis: couple one end of the chain to a Markovian (Lindblad) bath and measure how the slowest relaxation rate Γ(L) decays with L. Compare that decay with the many-body level-spacing scale and extrapolate to where avalanches stop. For the Heisenberg chain this gives W_c ≳ 20 (number to verify).
  - Resonance analysis: two drifting crossovers. W₁*(L) marks the onset of system-spanning resonances in rare samples; W₂*(L) marks where they are typical, which is the ED crossover.
  - Result: an "MBL regime" between the finite-size crossovers and W_c, in which finite systems look localized but the asymptotic fate is thermal. W_c is a lower bound, not a proof of existence [memory].
- **Crowley–Chandran, SciPost Phys. 2022** (H, compared with ED): a constructive theory of the numerically accessible crossover. Rare few-level resonances reproduce ED eigenstate statistics. The crossover drifts slowly (logarithmically in L in their theory; check this) and is not the asymptotic transition [memory].
- **Long–Crowley–Khemani–Chandran, PRL 2023** (H/N): a phenomenology of the prethermal MBL regime.
  - Dynamics proceeds through sparse many-body resonances of growing range.
  - It gives explicit crossover scales separating localized, resonant and thermalizing behaviour in the (W, L) and (W, t) planes.
  - The statistical Jacobi method of Long–Crowley–Chandran (SciPost 2023, id?) extends their reach in L [memory].
- **Ha–Morningstar–Huse, PRL 2023** (N): in a bath-coupled chain, the avalanche front advances through many-body resonances. This supports avalanches remaining effective far above W_c^ED [memory].
- **Sierant–Zakrzewski, PRB 2022, "Challenges to observation of MBL"** (N; tensor-network and ED dynamics): imbalance decay persists at W much larger than W_c^ED. Telling MBL apart from slow thermalization would need times far beyond experiments and numerics [memory].
- **Sierant–Lewenstein–Scardicchio–Zakrzewski, PRB 2023, Floquet stability** (id?; N): Floquet models reach larger effective sizes and show similar drifts [memory].
- **Estimates of W_c for the Heisenberg chain with h ∈ [−W, W]** (all [memory], verify):

  | Method | W_c |
  |---|---|
  | ED spectral, L ≤ 22 | ≈ 3.7 |
  | POLFED, L ≤ 24, 1/L extrapolation | ≈ 5.4 |
  | Tensor-network dynamics, L ≈ 50–200, t ≲ 10²–10³ (Doggen et al. 2018; Chanda et al. 2020): onset of apparently frozen dynamics | ≈ 4–5.5 |
  | Avalanches and bath coupling (Morningstar 2022; Sels 2022) | ≳ 18–20 |
  | Šuntajs et al. | none (W* ∝ L) |

- **Sierant–Lewenstein–Scardicchio–Vidmar–Zakrzewski, Rep. Prog. Phys. 2025, "MBL in the age of classical computing"**: a comprehensive review of ED (L ≤ 24–26), tensor networks, Floquet models and avalanches. From memory, it treats the MBL regime as established and the existence of the phase as open. I could not check its exact conclusions [memory].

### 4. Avalanche theory

**Core theory**
- **De Roeck–Huveneers, PRB 2017** (H). Setup and mechanism:
  - A rare ergodic grain of n₀ spins (probability ~e^{−c n₀}, so grains occur at finite density) acts as a bath.
  - It absorbs an ℓ-bit at distance ℓ when the coupling g e^{−ℓ/ζ} exceeds the grain's level spacing, roughly 2^{−(n₀+ℓ)}.
- Consequences:
  - The grain runs away if ζ > ζ_c = 1/ln 2. This uses the amplitude-decay convention with local dimension 2; other conventions give 2/ln 2.
  - MBL therefore requires ζ < ζ_c everywhere.
  - In d ≥ 2 the absorbed region grows faster than the decay, so avalanches always win and there is no asymptotic MBL [memory].
- **Luitz–Huveneers–De Roeck, PRL 2017** (N; ED, ~15 spins): an ℓ-bit chain coupled to a small random-matrix dot confirms the threshold ζ_c [memory].
- **Thiery–Huveneers–Müller–De Roeck, PRL 2018** (H/N; avalanche-based RG): avalanches drive the transition, and ζ stays finite (equal to ζ_c) at the critical point, unlike Anderson localization [memory].
- **Gopalakrishnan–Huse, PRB 2019** (H): the avalanche instability is a sharp transition in a non-standard thermodynamic limit (a finite MBL system coupled to a bath) [memory].
- **KT-type RG** (H/N on RG models):
  - Earlier RGs (Vosk–Huse–Altman 2015; Potter–Vasseur–Parameswaran 2015) gave ν ≈ 3.
  - Avalanche-informed RGs give KT-like scaling, with ξ ~ exp(c/√|W − W_c|): Goremykina–Vasseur–Serbyn 2019, Dumitrescu et al. 2019, Morningstar–Huse 2019, Morningstar–Huse–Imbrie 2020.
  - In this picture the MBL phase is a line of fixed points and the thermal fraction jumps at the transition.
  - ED data with KT fits (Šuntajs PRB 2020; Laflorencie–Lemarié–Macé PRR 2020, "chain breaking") are consistent, but the fitted W_c also drifts [memory].
- **Léonard et al., Nature Physics 2023** (E): ⁸⁷Rb in a 1D Bose–Hubbard chain imaged with a quantum-gas microscope. A clean, thermal region sits next to a disordered region. Bath-induced delocalization invades the localized region and slows down with W, which is the onset of an avalanche. The system is ~10 sites and runs to ~10² tunnelling times, so asymptotic stability is not accessible. Sels and Demler are co-authors (author list from memory) [memory].

**Extensions** (mostly H)
- **d > 1:**
  - Theory: the De Roeck–Huveneers argument says unstable.
  - Numerics see MBL-like behaviour on accessible scales: Wahl–Pal–Simon 2019 (2D tensor networks), Doggen et al. PRL 2020 (slow 2D delocalization), Théveniaut et al. 2020 (2D dimer model).
  - Experiments (Choi 2016; Bordia 2017) agree.
  - Common reading: a regime with astronomically long timescales [memory].
- **Continuum and mobility edges:** De Roeck–Huveneers–Müller–Schiulaz, PRB 2016. Thermal bubbles from the delocalized part of the spectrum destroy many-body mobility edges, so there is no MBL in the continuum with unbounded energy. Possible counterexample: "quantum sun" models with a mobility edge (Pawlik et al. 2024, id?) [memory].
- **Long-range couplings:**
  - Literature: Burin 2006; Yao et al. 2014; Gutman et al. 2016; Nandkishore–Sondhi PRX 2017; Tikhonov–Mirlin PRB 2018.
  - Mechanism: with 1/r^α couplings, resonances proliferate below thresholds (Burin: α < 2d for interactions).
  - Consequence: power-law ℓ-bit tails lose to the 2^ℓ growth of the bath, so MBL needs effectively exponential tails [memory].
- **Non-abelian symmetry:** Potter–Vasseur, PRB 2016. SU(2) multiplets are incompatible with an ℓ-bit structure, so the system either breaks the symmetry or thermalizes. At finite size there is a broad non-ergodic regime with slow thermalization (Protopopov–Ho–Abanin 2017; Protopopov–Panda–Abanin PRX 2020) [memory].
- **Quasiperiodic versus random:**
  - Khemani–Sheng–Huse PRL 2017 (N): two universality classes. Quasiperiodic potentials have no Griffiths rare regions, so the avalanche seed is absent, and ν ≈ 1 is allowed.
  - Agrawal–Gopalakrishnan–Vasseur, Nat. Commun. 2020 (RG). Agrawal–Vasseur–Gopalakrishnan, PRB 2022 (id?): quasiperiodic MBL may be stable even in d > 1.
  - Numerics: Žnidarič–Ljubotina PNAS 2018 (interacting Aubry–André); Doggen–Mirlin PRB 2019 (long chains); Weiner–Evers–Bera 2019; Aramthottil et al. 2021; Tu–Vu–Das Sarma 2023 (avalanches in quasiperiodic systems, id?).
  - Status: 1D quasiperiodic chains are the best candidate for a genuine phase, though finite-size drifts appear there too [memory].

### 5. Recent developments, 2023–2026 (coverage incomplete)

- **Avalanche toy models with clean transitions** (N), the "quantum sun" models:
  - Šuntajs–Vidmar PRL 2022, "Ergodicity breaking transition in zero dimensions": a model built on the avalanche mechanism has a sharp transition with an L-independent critical point.
  - Hopjan–Vidmar PRL 2023: scale-invariant survival probability at eigenstate transitions.
  - Šuntajs–Hopjan–De Roeck–Vidmar PRR 2024 (id?): link to ultrametric random matrices.
  - Pawlik–Sierant–Vidmar–Zakrzewski PRB 2024 (id?): a many-body mobility edge.
  - Significance: avalanche-type transitions exist in controlled models. This does not decide the XXZ chain [memory].
- **Thermal-inclusion numerics** (N): Peacock–Sels PRB 2023 (id?); Szołdra–Sierant–Lewenstein–Zakrzewski PRB 2024, "Catching thermal avalanches in the disordered XXZ model" (id?); Colmenarez–Luitz–De Roeck 2024 (inclusions; title and id uncertain). Inclusions spread for W well above W_c^ED; the quantitative thresholds differ between papers [memory].
- **Calibration on a known-localized system:** Šuntajs–Prosen–Vidmar PRB 2023, "Localization challenges quantum chaos in the finite two-dimensional Anderson model" (N). The 2D Anderson model, localized for every W in the thermodynamic limit, shows a drifting finite-size "transition". This weakens the drift argument in both directions [memory].
- **Fock-space and graph approaches** (H/N):
  - Background: Altshuler–Gefen–Kamenev–Levitov 1997 (Fock space as a tree).
  - Recent work: Roy–Logan (PRB 2020 and later; strong Fock-space correlations of the site energies); Tarzia PRB 2020; De Tomasi–Khaymovich–Pollmann–Warzel PRB 2021 (rare thermal bubbles seen in Fock space); Tikhonov–Mirlin Ann. Phys. 2021 (review linking random regular graphs to MBL); Scoquart–Gornyi–Mirlin PRB 2024 (id?); Vanoni–Altshuler–Kravtsov–Scardicchio PNAS 2024 (methodological RG on random regular graphs).
  - Mechanism: the Fock-space graph has connectivity ~L and strongly correlated energies, and those correlations are what keep the system localized.
  - Limitation: these theories reproduce the ED crossovers, but real-space rare-region avalanches are hard to capture in self-consistent Fock-space theories [memory].
- **Quantum processors** (E):
  - Experiments: Chiaro et al. PRR 2022 (Google; LIOM interactions measured); Mi et al. Nature 2022 (time-crystalline eigenstate order protected by MBL); Shtanko et al. 2023 (IBM, up to ~100 qubits, "uncovering local integrability", id?); Zhejiang/USTC 2D superconducting experiments (id?).
  - Reach: ~10–100 qubits and ≲10³ cycles, limited by noise.
  - None claims to settle the thermodynamic-limit question [memory].
- **Methods at scale:** flow equations (Thomson–Schiró 2018; Thomson–Eisert Nat. Phys. 2024, id?) truncate at quartic order around free fermions and reach L ~ 10². They miss resonances by construction [memory].
- **Rigorous** (P), titles and claims to verify:
  - De Roeck–Huveneers–Prośniak (2023/24, id?), "Rigorous and simple results on very slow thermalization, or quasi-localization, of the disordered quantum chain": unconditional bounds showing very slow dynamics at strong disorder, without LLA (quasi-localization, not localization).
  - Yin–Nandkishore–Lucas 2024 (id?), "Eigenstate localization in a many-body quantum system": a specific-model result; my memory of it is vague.
  - Mastropietro CMP 2017: localization of ground-state correlations for interacting quasiperiodic fermions [memory].
- **Claims that the question is settled:** I know of no credible one, in either direction, up to my training data [memory].

---

## (c) Strongest arguments for and against, and what each side concedes

**For a true 1D MBL phase (short-range, random):**
1. **(P, conditional) Imbrie 2016** establishes MBL at strong disorder, given LLA. This is the only theorem-level statement.
2. **(H) Resonance counting.** The probability of a resonance at range r is about (C J/W)^r times the number of channels, so ζ(W) ~ 1/ln(W/J) → 0 as W → ∞. Eventually ζ < ζ_c, which satisfies the avalanche criterion. The theory that describes how MBL fails also predicts where it survives.
3. **(H/N) Self-consistent critical theory.** The avalanche-based RG gives a KT-type transition with a stable MBL line of fixed points. Quantum-sun models show that avalanche-type transitions with clean scaling exist (N).
4. **(N) Strong-disorder numerics.** For W ≳ 20–30, no numerics show thermalization on accessible scales. This is weak evidence, since the accessible scales are short.
5. **(H) Quasiperiodic chains** remove the rare-region seed for avalanches.

**Against, i.e. only a regime:**
1. **(N) Drift.** Every finite-size estimator drifts monotonically with L, with no convergent crossing (Šuntajs et al.; also acknowledged by the other side).
2. **(N) No divergence seen.** At W much larger than W_c^ED, dynamics decay slowly with relaxation times ~e^{W/W₀}, and nothing indicates divergence at finite W (Sels–Polkovnikov; Vidmar–Mierzejewski; Evers–Bera).
3. **(N/H) Moving lower bound.** The lower bound on W_c has risen each time the analysis improved: 3.7 → 5.4 → ≳ 20.
4. **(P gap)** LLA is unproven. I know of no mechanism that violates it, but there is no proof even for nearly free blocks.
5. **(H) Pre-asymptotic data.** The resonance phenomenology (Crowley–Chandran; Long et al.) implies that all data at W ≲ 10 are pre-asymptotic. This cuts both ways: it also removes the direct numerical evidence for the phase.

**What each side concedes** [memory]:
- **Pro-phase side:** W_c^ED is not the transition; W_c, if it exists, is of order 10–20 or more; there is no direct numerical or experimental evidence for the phase; and MBL is fragile, absent in d ≥ 2, with long-range couplings, and with SU(2).
- **Skeptic side:** relaxation times are extremely long (≳ e^{W}) and LIOMs describe the accessible dynamics. Their methods cannot exclude a transition at W ≳ 20; Sels–Polkovnikov phrase their result as "if it exists, at much larger W". "No MBL at any W" is an extrapolation.
- **Common ground:** both sides agree on a large-but-finite-time (and finite-L) "MBL regime" and on avalanches and resonances as the mechanisms.

---

## (d) Open problems suited to a free-fermion, Gaussian-state, rigorous-minded theorist

1. **The resonance network at L ~ 10³–10⁴, built on Anderson orbitals (weak Δ).**
   - Setup: the random-field XX chain is, after Jordan–Wigner, the Anderson model. The Δ n_i n_{i+1} interaction has matrix elements V_αβγδ between exact orbitals.
   - Compute the statistics of resonant quartets and multi-step resonances on the Slater-determinant graph, feed them into an avalanche or percolation model (Crowley–Chandran, Long et al., or Thiery-type RG), and extract ζ(W, Δ).
   - Sharp prediction to test: the avalanche threshold depends on ζ, not on the prefactor Δ. So W_c(Δ → 0⁺) may stay finite and nonzero, a discontinuity at Δ = 0 tied to ξ_sp(W) ≈ O(1). Colbois–Alet–Laflorencie 2024 may have studied this ("Interaction-driven instabilities in the random-field XXZ chain", id?, unverified).
2. **The boundary-bath diagnostic beyond Gaussian order.**
   - Morningstar et al.'s diagnostic is the slowest rate of a boundary-Lindblad chain.
   - For quadratic Lindbladians it is exact via third quantization (Prosen 2008): Γ ∝ e^{−2L/ξ}, with no many-body enhancement.
   - Compute the O(Δ²) many-body channel analytically or semi-analytically at L ~ 10², to decide whether the decay length crosses ζ_c in a controlled corner.
3. **A rigorous number-entropy versus transport inequality.**
   - S_N is at most the entropy of the maximum-entropy integer distribution with the same variance, so S_N ~ ln ln t would require Var N_A(t) ~ ln t.
   - For Gaussian states, full counting statistics are sums of Bernoullis (Klich–Levitov), which gives closed forms.
   - Goal: reduce the Kiefer-Emmanouilidis versus Luitz–Bar Lev dispute to bounds on charge fluctuations across a cut. Prove those are bounded up to times exp(c L ln W) at strong disorder (Lieb–Robinson-type, building on Hamza–Sims–Stolz).
4. **Rigorously bound the drift.**
   - Prove finite-L, finite-time statements: for L ≤ L₀(W), transport or relaxation needs t ≥ (W/J)^{cL}.
   - This would show W*(L) grows at most logarithmically, refuting linear drift without settling the phase. It is the Abanin et al. argument made rigorous, and quasi-localization methods (De Roeck–Huveneers-type) are the natural starting point.
5. **LLA for nearly free blocks.**
   - Many-body levels are E = Σ ε_α n_α. Use Wegner and Minami estimates for Anderson eigenvalues, plus combinatorics of signed sums, to prove an LLA-type bound for the XX chain.
   - Then extend perturbatively to small Δ. That would make Imbrie-type MBL unconditional in a weak-interaction, strong-disorder corner (if the multiscale scheme can be adapted to that setting).
6. **Calibrate chaos diagnostics against exactly localized baselines.**
   - Compute AGP norms, spectral functions and Thouless-type scales for perturbed free systems with Gaussian methods at large L (e.g. 2D Anderson, cf. Šuntajs–Prosen–Vidmar 2023).
   - This quantifies which drifts a truly localized system produces, which is needed to interpret the Sels–Polkovnikov and Šuntajs diagnostics.
7. **A non-Gaussianity diagnostic for resonances.**
   - The one-particle density-matrix occupations (Bera–Schomerus–Heidrich-Meisner–Bardarson PRL 2015) and the distance of eigenstates or evolved states from the Gaussian manifold flag resonances.
   - Hartree–Fock and Gaussian time-dependent variational methods (Weidinger–Gopalakrishnan–Knap 2018) give large-L baselines with no resonances. Exact minus Gaussian isolates the resonant contribution, whose rare-sample statistics can be compared with Crowley–Chandran.
8. **Gaussian-closed "null models" for avalanches.**
   - Anderson chain plus dephasing (L_i = n_i) on a finite inclusion: the two-point function obeys closed equations at L ~ 10⁴ (cf. Taylor–Scardicchio 2021; Žnidarič–Scardicchio–Varma 2016).
   - This gives inclusion-induced relaxation without the 2^ℓ many-body enhancement, a baseline that separates genuinely many-body avalanche physics from single-particle leakage. Monitored Gaussian trajectories are a related variant.
9. **A minimal avalanche model, speculative (H).** An Anderson-localized free chain hybridized with a single small chaotic or interacting dot. The dot mediates all inelastic processes, and impurity-type tools (Bravyi–Gosset: low Gaussian rank) may make it tractable at large L, giving a candidate model where the avalanche criterion could be checked near-rigorously.
10. **Interacting quasiperiodic chains.** There is no rare-region seed. Use Gaussian-based perturbation theory and flow equations for the interacting Aubry–André model to compute ℓ-bit tails at large L, and try to extend Mastropietro-type constructions from ground states to excited states. This is plausibly the most tractable route to a genuine-phase theorem.

---

## (e) References

All entries are [memory]; journal details and ids are unverified this session.

**Standard picture and reviews**
- Basko, Aleiner, Altshuler, "Metal–insulator transition in a weakly interacting many-electron system with localized single-particle states", Ann. Phys. 321, 1126 (2006), cond-mat/0506617 [memory]
- Gornyi, Mirlin, Polyakov, "Interacting electrons in disordered wires: Anderson localization and low-T transport", PRL 95, 206603 (2005), cond-mat/0506411 [memory]
- Oganesyan, Huse, "Localization of interacting fermions at high temperature", PRB 75, 155111 (2007), cond-mat/0610854 [memory]
- Žnidarič, Prosen, Prelovšek, "Many-body localization in the Heisenberg XXZ magnet in a random field", PRB 77, 064426 (2008), 0706.2539 [memory]
- Pal, Huse, "Many-body localization phase transition", PRB 82, 174411 (2010), 1010.1992 [memory]
- Bardarson, Pollmann, Moore, "Unbounded growth of entanglement in models of many-body localization", PRL 109, 017202 (2012), 1202.5532 [memory]
- Serbyn, Papić, Abanin, "Universal slow growth of entanglement in interacting strongly disordered systems", PRL 110, 260601 (2013), 1304.4605 [memory]
- Serbyn, Papić, Abanin, "Local conservation laws and the structure of the many-body localized states", PRL 111, 127201 (2013), 1305.5554 [memory]
- Huse, Nandkishore, Oganesyan, "Phenomenology of fully many-body-localized systems", PRB 90, 174202 (2014), 1305.4915 [memory]
- Ros, Müller, Scardicchio, "Integrals of motion in the many-body localized phase", Nucl. Phys. B 891, 420 (2015), id? [memory]
- Luitz, Laflorencie, Alet, "Many-body localization edge in the random-field Heisenberg chain", PRB 91, 081103(R) (2015), 1411.0660 [memory]
- Chandran, Laumann, Oganesyan, "Finite size scaling bounds on many-body localized phase transitions" (2015), id? [memory]
- Khemani, Lim, Sheng, Huse, "Critical properties of the many-body localization transition", PRX 7, 021013 (2017), id? [memory]
- Pietracaprina, Macé, Luitz, Alet, "Shift-invert diagonalization of large many-body localizing spin chains", SciPost Phys. 5, 045 (2018), id? [memory]
- Imbrie, "On many-body localization for quantum spin chains", J. Stat. Phys. 163, 998 (2016), 1403.7837 [memory]
- Nandkishore, Huse, "Many-body localization and thermalization in quantum statistical mechanics", Annu. Rev. CMP 6, 15 (2015), 1404.0686 [memory]
- Alet, Laflorencie, "Many-body localization: an introduction and selected topics", C. R. Physique 19, 498 (2018), 1711.03145 [memory]
- Abanin, Altman, Bloch, Serbyn, "Colloquium: Many-body localization, thermalization, and entanglement", RMP 91, 021001 (2019), 1804.11065 [memory]

**Experiments**
- Schreiber et al., Science 349, 842 (2015), 1501.05661 [memory]
- Smith et al., Nat. Phys. 12, 907 (2016), 1508.07026 [memory]
- Choi et al., Science 352, 1547 (2016), 1604.04178 [memory]
- Bordia et al., PRX 7, 041047 (2017), id? [memory]
- Roushan et al., Science 358, 1175 (2017), id? [memory]
- Lukin et al., Science 364, 256 (2019), 1805.09819 [memory]
- Rispoli et al., Nature 573, 385 (2019), 1812.06959 [memory]
- Rubio-Abadal et al., PRX 9, 041014 (2019), id? [memory]
- Guo et al., Nat. Phys. 17, 234 (2021), id? [memory]
- Chiaro et al., PRR 4, 013148 (2022), id? [memory]
- Mi et al., Nature 601, 531 (2022), id? [memory]
- Léonard, Kim, Rispoli, Lukin, Schittko, Kwan, Demler, Sels, Greiner, "Probing the onset of quantum avalanches in a many-body localized system", Nat. Phys. 19, 481 (2023), 2012.15270 [memory]
- Shtanko et al., "Uncovering local integrability in quantum many-body dynamics" (IBM, 2023), id? [memory]

**The challenge and responses**
- Šuntajs, Bonča, Prosen, Vidmar, "Quantum chaos challenges many-body localization", PRE 102, 062144 (2020), 1905.06345 [memory]
- Šuntajs, Bonča, Prosen, Vidmar, "Ergodicity breaking transition in finite disordered spin chains", PRB 102, 064207 (2020), 2004.01719 [memory]
- Abanin, Bardarson, De Tomasi, Gopalakrishnan, Khemani, Parameswaran, Pollmann, Potter, Serbyn, Vasseur, "Distinguishing localization from chaos: challenges in finite-size systems", Ann. Phys. 427, 168415 (2021), 1911.04501 [memory]
- Panda, Scardicchio, Schulz, Taylor, Žnidarič, "Can we study the many-body localisation transition?", EPL 128, 67003 (2019; published 2020), 1911.07882 [memory]
- Sierant, Lewenstein, Zakrzewski, "Polynomially filtered exact diagonalization approach to many-body localization", PRL 125, 156601 (2020), 2005.09534 [memory]
- Sierant, Delande, Zakrzewski, "Thouless time analysis of Anderson and many-body localization transitions", PRL 124, 186601 (2020), id? [memory]
- Sels, Polkovnikov, "Dynamical obstruction to localization in a disordered spin chain", PRE 104, 054105 (2021), 2009.04501 [memory]
- Sels, "Bath-induced delocalization in interacting disordered spin chains", PRB 106, L020202 (2022), 2108.10796 [memory]
- Sels, Polkovnikov, "Thermalization of dilute impurities in one-dimensional spin chains", PRX 13, 011041 (2023), id? [memory]
- LeBlond, Sels, Polkovnikov, Rigol, "Universality in the onset of quantum chaos in many-body systems", PRB 104, L201117 (2021), id? [memory]
- Bera, De Tomasi, Weiner, Evers, "Density propagator for many-body localization: finite-size effects, transient subdiffusion, and exponential decay", PRL 118, 196801 (2017), id? [memory]
- Weiner, Evers, Bera, "Slow dynamics and strong finite-size effects in many-body localization with random and quasiperiodic potentials", PRB 100 (2019), id? [memory]
- Evers, Modak, Bera, "Internal clock of many-body delocalization", PRL (2023), id? (venue and year uncertain) [memory]
- Prelovšek, Barišić, Žnidarič, "Absence of full many-body localization in the disordered Hubbard chain", PRB 94, 241104(R) (2016), id? [memory]
- Mierzejewski, Herbrych, Prelovšek, "Universal dynamics of density correlations at the transition to the many-body localized state", PRB 94, 224207 (2016), id? [memory]
- Prelovšek, Mierzejewski, Barišić, Herbrych, "Density correlations and transport in models of many-body localization", Ann. Phys. (Berlin) 529, 1600362 (2017), id? [memory]
- Herbrych, Mierzejewski, Prelovšek, "Relaxation at different length-scales in models of many-body localization", PRB 105 (2022), id? [memory]
- Vidmar, Krajewski, Bonča, Mierzejewski, "Phenomenology of spectral functions in disordered spin chains at infinite temperature", PRL 127, 230603 (2021), id? [memory]
- Krajewski, Vidmar, Bonča, Mierzejewski, "Restoring ergodicity in a strongly disordered interacting chain", PRL 129, 260601 (2022), id? [memory]
- Kiefer-Emmanouilidis, Unanyan, Fleischhauer, Sirker, "Evidence for unbounded growth of the number entropy in many-body localized phases", PRL 124, 243601 (2020), id? [memory]
- Luitz, Bar Lev, "Absence of slow particle transport in the many-body localized phase", PRB 102, 100202(R) (2020), id? [memory]
- Ghosh, Žnidarič, "Resonance-induced growth of number entropy in strongly disordered systems", PRB 105, 144203 (2022), id? [memory]

**Regime versus phase**
- Morningstar, Colmenarez, Khemani, Luitz, Huse, "Avalanches and many-body resonances in many-body localized systems", PRB 105, 174205 (2022), 2107.05642 [memory]
- Crowley, Chandran, "A constructive theory of the numerically accessible many-body localized to thermal crossover", SciPost Phys. 12, 201 (2022), 2012.14393 [memory]
- Crowley, Chandran, "Avalanche induced coexisting localized and thermal regions in disordered chains", PRR 2, 033262 (2020), id? [memory]
- Long, Crowley, Khemani, Chandran, "Phenomenology of the prethermal many-body localized regime", PRL 131, 106301 (2023), 2207.05761 [memory]
- Long, Crowley, Chandran, "Beyond Fermi's golden rule with the statistical Jacobi approximation", SciPost Phys. (2023), id? [memory]
- Ha, Morningstar, Huse, "Many-body resonances in the avalanche instability of many-body localization", PRL 130, 250405 (2023), 2211.16465 [memory]
- Sierant, Zakrzewski, "Challenges to observation of many-body localization", PRB 105, 224203 (2022), id? [memory]
- Sierant, Lewenstein, Scardicchio, Zakrzewski, "Stability of many-body localization in Floquet systems", PRB 107 (2023), id? [memory]
- Sierant, Lewenstein, Scardicchio, Vidmar, Zakrzewski, "Many-body localization in the age of classical computing", Rep. Prog. Phys. 88, 026502 (2025), 2403.07111 [memory]
- Doggen, Schindler, Tikhonov, Mirlin, Neupert, Polyakov, Gornyi, "Many-body localization and delocalization in large quantum chains", PRB 98, 174202 (2018), id? [memory]
- Chanda, Sierant, Zakrzewski, "Time dynamics with matrix product states: many-body localization transition of large systems revisited", PRB 101, 035148 (2020), id? [memory]

**Avalanches, RG, extensions**
- De Roeck, Huveneers, "Stability and instability towards delocalization in many-body localization systems", PRB 95, 155129 (2017), 1608.01815 [memory]
- Luitz, Huveneers, De Roeck, "How a small quantum bath can thermalize long localized chains", PRL 119, 150602 (2017), 1705.10807 [memory]
- Thiery, Huveneers, Müller, De Roeck, "Many-body delocalization as a quantum avalanche", PRL 121, 140601 (2018), 1706.09338 [memory]
- Gopalakrishnan, Huse, "Instability of many-body localized systems as a phase transition in a nonstandard thermodynamic limit", PRB 99, 134305 (2019), id? [memory]
- Vosk, Huse, Altman, PRX 5, 031032 (2015), 1412.3117 [memory]
- Potter, Vasseur, Parameswaran, PRX 5, 031033 (2015), 1501.03501 [memory]
- Goremykina, Vasseur, Serbyn, PRL 122, 040601 (2019), id? [memory]
- Dumitrescu, Goremykina, Parameswaran, Serbyn, Vasseur, "Kosterlitz-Thouless scaling at many-body localization phase transitions", PRB 99, 094205 (2019), id? [memory]
- Morningstar, Huse, PRB 99, 224205 (2019), id? [memory]
- Morningstar, Huse, Imbrie, PRB 102, 125134 (2020), id? [memory]
- Laflorencie, Lemarié, Macé, "Chain breaking and Kosterlitz-Thouless scaling at the many-body localization transition…", PRR 2, 042033(R) (2020), id? [memory]
- De Roeck, Huveneers, Müller, Schiulaz, "Absence of many-body mobility edges", PRB 93, 014203 (2016), id? [memory]
- Burin, "Energy delocalization in strongly disordered systems induced by long-range many-body interaction" (2006), cond-mat/0611387 (id from memory) [memory]
- Yao et al., "Many-body localization in dipolar systems", PRL 113, 243002 (2014), id? [memory]
- Gutman et al., PRB 93, 245427 (2016), id? [memory]
- Nandkishore, Sondhi, "Many-body localization with long-range interactions", PRX 7, 041021 (2017), id? [memory]
- Tikhonov, Mirlin, PRB 97, 214205 (2018), id? [memory]
- Potter, Vasseur, "Symmetry constraints on many-body localization", PRB 94, 224206 (2016), id? [memory]
- Protopopov, Ho, Abanin, PRB 96, 041122(R) (2017), id? [memory]
- Protopopov, Panda, Abanin, PRX 10, 011025 (2020), id? [memory]
- Iyer, Oganesyan, Refael, Huse, "Many-body localization in a quasiperiodic system", PRB 87, 134202 (2013), 1212.4159 [memory]
- Khemani, Sheng, Huse, "Two universality classes for the many-body localization transition", PRL 119, 075702 (2017), id? [memory]
- Agrawal, Gopalakrishnan, Vasseur, "Universality and quantum criticality in quasiperiodic spin chains", Nat. Commun. 11, 2225 (2020), id? [memory]
- Agrawal, Vasseur, Gopalakrishnan, "Quasiperiodic many-body localization transition in dimension d > 1", PRB 106 (2022), id? [memory]
- Žnidarič, Ljubotina, "Interaction instability of localization in quasiperiodic systems", PNAS 115, 4595 (2018), id? [memory]
- Doggen, Mirlin, "Many-body delocalization dynamics in long Aubry-André quasiperiodic chains", PRB 100 (2019), id? [memory]
- Aramthottil, Chanda, Sierant, Zakrzewski, PRB 104 (2021), id? [memory]
- Tu, Vu, Das Sarma, "Avalanche stability transition in interacting quasiperiodic systems", PRB (2023), id? [memory]
- Wahl, Pal, Simon, Nat. Phys. 15, 164 (2019), id? [memory]
- Doggen et al., "Slow many-body delocalization beyond one dimension", PRL 125, 155701 (2020), id? [memory]
- Théveniaut, Lan, Meyer, Alet, PRR 2, 033154 (2020), id? [memory]

**Recent work, Fock space, rigorous results, methods**
- Šuntajs, Vidmar, "Ergodicity breaking transition in zero dimensions", PRL 129, 060602 (2022), id? [memory]
- Hopjan, Vidmar, "Scale-invariant survival probability at eigenstate transitions", PRL 131, 060404 (2023), id? [memory]
- Šuntajs, Hopjan, De Roeck, Vidmar, PRR (2024), id? [memory]
- Pawlik, Sierant, Vidmar, Zakrzewski, "Many-body mobility edge in quantum sun models", PRB 109 (2024), id? [memory]
- Šuntajs, Prosen, Vidmar, "Localization challenges quantum chaos in the finite two-dimensional Anderson model", PRB 107, 064205 (2023), id? [memory]
- Peacock, Sels, "Many-body delocalization from embedded thermal inclusion", PRB (2023), id? [memory]
- Szołdra, Sierant, Lewenstein, Zakrzewski, "Catching thermal avalanches in the disordered XXZ model", PRB (2024), id? [memory]
- Colbois, Alet, Laflorencie, "Interaction-driven instabilities in the random-field XXZ chain" (2024), id?; title and claim unverified [memory]
- Altshuler, Gefen, Kamenev, Levitov, PRL 78, 2803 (1997), id? [memory]
- Roy, Logan, "Fock-space correlations and the origins of many-body localization", PRB 101, 134202 (2020), id? [memory]
- Tarzia, "Many-body localization transition in Hilbert space", PRB 102, 014208 (2020), id? [memory]
- De Tomasi, Khaymovich, Pollmann, Warzel, "Rare thermal bubbles at the many-body localization transition from the Fock space point of view", PRB 104, 024202 (2021), id? [memory]
- Tikhonov, Mirlin, "From Anderson localization on random regular graphs to many-body localization", Ann. Phys. 435, 168525 (2021), id? [memory]
- Scoquart, Gornyi, Mirlin, "Role of Fock-space correlations in many-body localization", PRB 109 (2024), id? [memory]
- Vanoni, Altshuler, Kravtsov, Scardicchio, "Renormalization group analysis of the Anderson model on random regular graphs", PNAS 121 (2024), id? [memory]
- Macé, Alet, Laflorencie, "Multifractal scalings across the many-body localization transition", PRL 123, 180601 (2019), id? [memory]
- Hamza, Sims, Stolz, "Dynamical localization in disordered quantum spin systems", CMP 315, 215 (2012), id? [memory]
- Mastropietro, "Localization in interacting fermionic chains with quasi-random disorder", CMP 351, 283 (2017), id? [memory]
- De Roeck, Huveneers, Prośniak, "Rigorous and simple results on very slow thermalization, or quasi-localization, of the disordered quantum chain" (2023/24), id? [memory]
- Yin, Nandkishore, Lucas, "Eigenstate localization in a many-body quantum system" (2024), id?; claim unverified [memory]
- Imbrie, Ros, Scardicchio, "Local integrals of motion in many-body localized systems", Ann. Phys. (Berlin) 529, 1600278 (2017), id? [memory]
- Bera, Schomerus, Heidrich-Meisner, Bardarson, "Many-body localization characterized from a one-particle perspective", PRL 115, 046603 (2015), id? [memory]
- Weidinger, Gopalakrishnan, Knap, "Self-consistent Hartree-Fock approach to many-body localization", PRB 98, 224205 (2018), id? [memory]
- Thomson, Schiró, PRB 97, 060201(R) (2018), id?; Thomson, Eisert, Nat. Phys. (2024), id? [memory]
- Žnidarič, Scardicchio, Varma, PRL 117, 040601 (2016), id?; Taylor, Scardicchio, PRB 103, 184202 (2021), id? [memory]
- Prosen, "Third quantization", New J. Phys. 10, 043026 (2008), id?; Klich, Levitov, PRL 102, 100502 (2009), id?; Bravyi, Gosset, "Complexity of quantum impurity problems", CMP 356, 451 (2017), id? [memory]

**Re-verification needed.** To upgrade tags to [read], allow arxiv.org (and ideally journals.aps.org and scipost.org) in the session's network settings, or re-run with a fresh web-search budget. Priority checks: (i) the W_c figures in Morningstar et al. 2022 and Sels 2022; (ii) the 1/L-extrapolated value in Sierant et al. 2020; (iii) the drift law in Crowley–Chandran 2022; (iv) the conclusions of the Rep. Prog. Phys. 2025 review; (v) whether the 2023–24 rigorous items exist and what they prove; (vi) any 2025–26 work.
