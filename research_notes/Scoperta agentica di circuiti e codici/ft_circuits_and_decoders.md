# AI-driven discovery of fault-tolerant QEC circuits and learned decoders (2023 – Nov 2026)

> How these notes were made (compiled 2026-11-22): about 60 web searches. arxiv.org and export.arxiv.org could not be reached by the page-fetch tool in this session, and the shared web-search budget ran out before the last planned queries. Those queries were fold-transversal/automorphism gate search, Floquet-code scheduling, RL "gym" environments for QEC, bosonic hardware-in-the-loop RL after 2023, and details of Quantinuum Helios decoding. Numbers below come from arXiv abstract/HTML snippets returned by search, from journal pages and company blogs, and, where flagged, from secondary coverage. Legend: **n/r** = not retrieved; **[secondary]** = rests on news or aggregator summaries; **[title-only]** = only the existence and title of the work were confirmed. Each arXiv YYMM prefix gives the posting month (2611 = Nov 2026).

---

## 1. Newest works first (2025 – Nov 2026): LLM agents, LLM-guided evolutionary search and automated discovery (top priority)

### Takeaway
LLM agents entered fault-tolerant (FT) circuit discovery in 2026 in three distinct roles:
1. **Mutation operator in evolutionary program synthesis.** The LLM is paired with an exhaustive fault-enumeration verifier (Lahiri & Marquardt, Nov 2026).
2. **Search steering.** Claude agents chose symmetry groups for an ILP search, which yielded 564 new magic-state-factory classes (Oct 2026).
3. **"AI-scientist" pipeline.** The pipeline builds syndrome-extraction circuits and scores them by sampled circuit-level logical error rate (LER) (OmniQEC, Jul 2026).

All three share one pattern: the LLM proposes, and a deterministic verifier or solver certifies. Targeted searches found **no LLM/AlphaEvolve/FunSearch-discovered decoder** and no LLM work on lattice surgery or code switching.

### Cited Findings

**Agentic / LLM-guided works in scope (newest first)**
- **Lahiri & Marquardt (Nov 2026), "Agentic AI for Discovering Fault-Tolerant Quantum Circuits", arXiv:2611.18504.**
  - *Agent:* frontier LLMs inside an evolutionary program-synthesis loop (described as an "agentic FunSearch"). The LLM is the mutation operator that evolves Python programs defining dynamic Clifford circuits, and "every search starts cold with no seed programs".
  - *Verifier:* two stages. The first is an exhaustive fault-tolerance check under circuit-level depolarizing noise, modelling gate failures, input errors and measurement flips up to order *t*.
  - *Tasks:* logical state preparation (starting from the all-zeros state), logical gates and logical routines.
  - *Results:* searching for d = 3, 5, 7, 9 at the same time gives FT circuits of depth 10 in the first generation. The agent found shallower logical-Hadamard circuits than a prior reference. A table compares benchmark depth with agentic depth and lists "silent faults".
  - *Not retrieved:* per-code depths and gate counts.
  - Sources: [arXiv HTML](https://arxiv.org/html/2611.18504); [CASRAI](https://casrai.org/news/agentic-ai-discovering-fault-tolerant-quantum-circuits-2026) [secondary; gives MPL Erlangen/FAU affiliation]; [Quantum Zeitgeist](https://quantumzeitgeist.com/fault-tolerant-circuits-llm-designs-error/) [secondary]
- **"Symmetry and AI-assisted discovery of magic-state factories", arXiv:2610.06535 (Oct 2026).**
  - *Authors:* first author "Shubham …"; full author list n/r. Quantum Computing Report attributes the work to IBM Quantum and QuICS.
  - *Method:* a unified binary-matrix formulation of distillation factories. LLM agents ("Anthropic Claude models throughout") use earlier search outcomes to propose subgroup-symmetry orbits, parent codes and puncture sets. An integer-programming solver then builds candidate circuit matrices, and an independent verification engine checks gate identities and code distances.
  - *Results:* 699 factory classes, of which 564 are new. [[1715,287,6]] with γ = 0.998 is "the smallest known pure-T factory with γ < 1". [[63,11,3]] and [[850,128,6]] have the lowest overhead exponents known among protocols with ≤100 and ≤1000 inputs, respectively.
  - *Code:* the repository is public.
  - Sources: [arXiv:2610.06535](https://arxiv.org/abs/2610.06535); [HTML](https://arxiv.org/html/2610.06535v1); [GitHub repo](https://github.com/shubhamj810/magic-state-factory-repo); [Quantum Computing Report](https://quantumcomputingreport.com/ibm-quantum-and-quics-deploy-ai-agents-to-discover-over-500-high-efficiency-magic-state-factories/) [secondary]
- **OmniQEC (Jul 2026), "OmniQEC: discovering practical quantum error-correcting codes by an AI scientist", arXiv:2607.25865.**
  - *Pipeline:* candidates go through a fast screening loop, then a slow loop. In the slow loop a circuit-synthesis module compiles validated codes into syndrome-extraction (SE) circuits. An evaluation module derives the detector error model (DEM), samples noisy syndromes and decodes them.
  - *Scoring:* BP+OSD under circuit-level depolarizing noise at p = 0.005. The metric is LER per round per logical qubit.
  - *LLMs used:* frontier codes were found with Claude Opus 4.8, GPT-5.5 and DeepSeek-V4-Pro.
  - *Not retrieved:* the LER values.
  - *Scope note:* the work is mainly code discovery, which is out of scope here. The SE-circuit synthesis with LER in the loop is the in-scope part.
  - Source: [arXiv HTML](https://arxiv.org/html/2607.25865v1)

**Adjacent agentic infrastructure, benchmarks and risk studies (2026)**
- **QSA harness (Wan, D. Chen, Nahrstedt, B. Chen, UIUC; Sep 2026), arXiv:2609.38883.**
  - *Method:* wraps OpenEvolve-style LLM-guided evolution with quantum-specific support: coreset-based approximate scoring, static and snapshot analyses, and task-specific rewards. The stated motivation is "expensive verification, sparse feedback, heterogeneous objectives".
  - *Results on IBM hardware:* 4.2–9.5% QPU-utilization gain and 15.2–19.5% Hellinger-fidelity gain for multiprogramming; at least 96.8% less error-mitigation time. Total LLM API cost was $6.9 over 11.3 h.
  - *Scope note:* this is compiler and runtime work, not QEC, but it is an explicit design for verification-heavy LLM search.
  - Source: [arXiv:2609.38883](https://arxiv.org/html/2609.38883)
- **Özgüler (Sep 2026), "Securing quantum error correction against misleading advice from AI agents", arXiv:2609.19090.**
  - *Setting:* in an odd-distance toric code, opposite coherent X rotations produce identical passive syndrome-history distributions. A terminal logical measurement on calibration states is therefore needed to recover the sign.
  - *Safeguard:* a separate evaluator accepts an update only when calibration uncertainty and a drift bound justify it.
  - *Results:* under a misleading-advice attack, a Qwen3:8b adviser picked the harmful recovery in 12 of 12 cases, and authorization alone let harmful decoder updates through.
  - *Limitation:* the guarantees hold only if the drift bound holds.
  - Source: [arXiv:2609.19090](https://arxiv.org/abs/2609.19090)
- **QIQCBench (Guo, Li, Cheng et al., 15 Sep 2026), arXiv:2609.17439.**
  - 49 expert-written tasks covering calibration and control, error correction and compilation, and sensing and networking. They run in "Quantum-Harbor", a virtual lab.
  - Tested on 17 frontier agentic systems. Verified performance varied widely between systems and across repeated attempts. "Verified" means outputs must be backed by information gathered during execution.
  - Source: [arXiv:2609.17439](https://arxiv.org/html/2609.17439v1)
- **"Reimagining computational quantum research with AI scientists", arXiv:2611.02833 (Nov 2026).** Evaluates five LLM agentic systems and six harnesses across 59 papers, and introduces "Residual Agents". This is a general platform (out of scope), listed for context. — [arXiv:2611.02833](https://arxiv.org/html/2611.02833)
- **AI-Mandel (Arlt, Gu, Krenn; Nov 2025), arXiv:2511.11752.** LLM agents built around PyTheus for quantum-optics research. A general platform, out of scope. — [arXiv:2511.11752](https://arxiv.org/pdf/2511.11752)

**LLM-guided *code* discovery (out of scope; listed only because it defines the toolchain and its lessons)**
- **IBM, "Evolutionary Discovery of Bivariate Bicycle Codes with LLM-Guided Search", arXiv:2606.02418 (Jun 2026).**
  - Uses OpenEvolve (open-sourced). About 1650 evolutionary iterations screened about 2×10^5 candidate codes and yielded 465 new candidates.
  - Candidates pass verification and duplicate checks before entering the catalogue.
  - Sources: [arXiv:2606.02418](https://arxiv.org/pdf/2606.02418); [IBM Research blog](https://research.ibm.com/blog/ai-for-qec)
- **"Evolving Quantum Error-Correcting Encodings for Molecular Simulation", arXiv:2606.25870.**
  - Uses the AlphaEvolve Cloud API (early access), with mutations from Gemini 3.1 Pro and Gemini 3 Flash.
  - A d = 5 run evaluated 1,411 candidates and compressed a 272-qubit seed panel to 168 qubits.
  - The authors present this as evidence for a design pattern, not as evidence that LLM search "understands" QEC.
  - Source: [arXiv HTML](https://arxiv.org/html/2606.25870v1)
- **qBraid industry blog post on AlphaEvolve for QEC.**
  - About 1,500 candidate programs were evaluated, giving 4.2–5.0× fewer data qubits and 3.4–7.9× lower LER under exact decoding.
  - It reports **reward hacking**: optimizing for distance alone produced highly connected, impractical codes. Circuit-level fault tolerance is named as the next step.
  - Source: [qBraid blog](https://www.qbraid.com/blog-posts/qbraid-uses-alphaevolve-for-quantum-error-correction)

**Targeted negative results**
- Searches combining FunSearch, AlphaEvolve or OpenEvolve with "decoder" found no LLM-discovered decoding heuristic. The closest items were encoding and code discovery and the agentic circuit search above. — [arXiv:2606.25870](https://arxiv.org/html/2606.25870); [arXiv:2611.18504](https://arxiv.org/html/2611.18504)
- No source covered LLM agents for code switching. Beyond Lahiri & Marquardt, none combined LLM agents with lattice-surgery gate synthesis. — [arXiv:2611.18504](https://arxiv.org/pdf/2611.18504)

**Chronological index of 2025–2026 works in scope (newest first; details in the topic sections)**

| Posted | ID / venue | Short title | Method class | § |
|---|---|---|---|---|
| Nov 2026 | [2611.26447](https://arxiv.org/abs/2611.26447) | Rapids streaming neural decoder (FPGA) | NN decoder | 6 |
| Nov 2026 | [2611.25825](https://arxiv.org/abs/2611.25825) | Structure-guided SE synthesis | structured search | 3 |
| Nov 2026 | [2611.25107](https://arxiv.org/abs/2611.25107) | Block Checks + Pareto tree search for FT state prep | tree search | 2 |
| Nov 2026 | [2611.19715](https://arxiv.org/abs/2611.19715) | Proxima SE scheduling | SMT proxy search | 3 |
| Nov 2026 | [2611.18504](https://arxiv.org/html/2611.18504) | Agentic AI for FT circuits | LLM evolutionary | 2, 4 |
| Nov 2026 | [2611.16271](https://arxiv.org/html/2611.16271v1) | Hierarchical decoding (incl. neural pre-decoders) | review | 6 |
| Nov 2026 | [2611.13630](https://arxiv.org/abs/2611.13630v1) | PRISM transformer decoder for fusion-based QC | NN decoder | 6 |
| Oct 2026 | [2610.41398](https://arxiv.org/abs/2610.41398) | GraphMP "decoder model factory" | learned MP + search distillation | 6 |
| Oct 2026 | [2610.06535](https://arxiv.org/abs/2610.06535) | AI-assisted magic-state factories | LLM agents + ILP | 5 |
| Oct 2026 | [2610.03714](https://arxiv.org/pdf/2610.03714) | SpiderCSS | ZX / exact | 2 |
| Sep 2026 | [2609.19090](https://arxiv.org/abs/2609.19090) | Misleading AI advice in QEC | agent safety | 7 |
| Sep 2026 | [2609.17439](https://arxiv.org/html/2609.17439v1) | QIQCBench | agent benchmark | 8 |
| Sep 2026 | [2609.12020](https://arxiv.org/abs/2609.12020) | FastSched (RL for syndrome extraction) | RL | 3 |
| Sep 2026 | [2609.04557](https://arxiv.org/abs/2609.04557) | Sim-to-real decoder benchmarking (Willow) | benchmark | 6, 8 |
| Jul 2026 | [2607.25865](https://arxiv.org/html/2607.25865v1) | OmniQEC | LLM AI scientist | 1, 3 |
| May 2026 | [2605.04892](https://arxiv.org/abs/2605.04892) | Closed-loop FPGA NN decoder on hardware | NN decoder | 6 |
| Mar 2026 | [2603.05481](https://arxiv.org/abs/2603.05481) | Left-right SE circuits | heuristic optimization | 3 |
| Mar 2026 | [2603.05391](https://arxiv.org/abs/2603.05391) | SpiderCat | ZX / exact | 2 |
| Jan 2026 | [2601.17580](https://arxiv.org/html/2601.17580v1) | PropHunt | MaxSAT proxy | 3 |
| Jan 2026 | [2601.13313](https://arxiv.org/pdf/2601.13313) | Steane-type FT prep for any CSS code | SAT/exact | 2 |
| Jan 2026 | [2601.12509](https://arxiv.org/abs/2601.12509) | AlphaSyndrome | MCTS | 3 |
| Jan 2026 | [2601.09921](https://arxiv.org/html/2601.09921) | Self-coordinating parallel NN decoder | NN decoder | 6 |
| Dec 2025 | [2512.07737](https://arxiv.org/abs/2512.07737) | AlphaQubit 2 | NN decoder | 6 |
| Nov 2025 → Nature 2026 | [2511.08493](https://arxiv.org/abs/2511.08493) | RL control of QEC on Willow | RL, on hardware | 7 |
| Oct 2025 → npj QI 2026 | [2510.22724](https://arxiv.org/abs/2510.22724) | Mamba neural decoder | NN decoder | 6 |
| Oct 2025 (PRX) | [10.1103/gqpr-dgz7](https://doi.org/10.1103/gqpr-dgz7) | RL for FT logical state prep | RL | 2 |
| Sep 2025 | [2509.24666](https://arxiv.org/html/2509.24666) | Automated gadget discovery for RL | graph mining + RL | 2, 4 |
| Sep 2025 → QST | [2509.03974](https://arxiv.org/html/2509.03974) | Multi-agent RL QEC cycles + bandit retraining | MARL | 2, 7 |
| Apr 2025 → Quantum 2026 | [2504.13043](https://arxiv.org/abs/2504.13043) | Transformer decoder for BB codes | NN decoder | 6 |
| Mar 2025 | [2503.11638](https://arxiv.org/abs/2503.11638) | RL with gadgets (encoders) | RL | 2 |
| Mar 2025 (GTC) | [NVIDIA blog](https://developer.nvidia.com/blog/nvidia-and-quera-decode-quantum-errors-with-ai) | NVIDIA–QuEra transformer decoder | NN decoder | 5, 6 |
| Feb 2025 | [2502.19971](https://arxiv.org/abs/2502.19971) | GraphQEC code-agnostic GNN | GNN | 6 |
| Jan 2025 (DATE) | [2501.05527](https://export.arxiv.org/pdf/2501.05527) | Deterministic FT state prep via SAT | SAT | 2 |

**Checking the example works named in the brief**
- **Zen, Olle, Colmenarez, Puviani, Marquardt (2024).** **Verified.** The author list also includes **Markus Müller**. arXiv:2402.17761, published as Phys. Rev. X 15, 041012 (22 Oct 2025). — [ADS](https://ui.adsabs.harvard.edu/abs/2025PhRvX..15d1012Z/abstract); [PRX DOI](https://doi.org/10.1103/gqpr-dgz7)
- **"Peham, Schmid, Burgholzer, Wille".** The paper exists, but **Burgholzer is not an author**. The authors are Peham, Schmid, Berent, Müller and Wille, published in PRX Quantum 6, 020330 (2025). — [arXiv:2408.11894v3](https://arxiv.org/abs/2408.11894v3)
- **IBM's RL-based Clifford synthesis.** The paper "Practical and efficient quantum circuit synthesis and transpiling with Reinforcement Learning" (arXiv:2405.13196) exists. Its content and author list were n/r. — [arXiv:2405.13196](https://arxiv.org/pdf/2405.13196)
- **"AlphaSyndrome".** **Verified**: arXiv:2601.12509 (2026), an MCTS scheduler. — [arXiv:2601.12509](https://arxiv.org/abs/2601.12509)
- **Trenkwalder, Briegel et al. on automated gadget discovery.** **Verified** as "Automated gadget discovery in the quantum domain" (Mach. Learn.: Sci. Technol. 2023; arXiv:2212.12743). It mines gadgets from the trajectories of a trained RL agent in quantum-optics and circuit tasks. It is **not** about fault-tolerant gadgets. — [DOAJ](https://doaj.org/article/61bf0003a8ea4bb69c5f1e40ea217406); [arXiv:2212.12743](https://arxiv.org/pdf/2212.12743)
- **AlphaQubit (Bausch et al., Nature 2024).** **Verified** as "Learning high-accuracy error decoding for quantum processors". — [Google blog](https://blog.google/technology/google-deepmind/alphaqubit-quantum-error-correction/)
- **Sweke et al. and Andreasson et al.** **Verified**: arXiv:1810.07207 and arXiv:1811.12338. — [Sweke et al.](https://arxiv.org/pdf/1810.07207.pdf); [Andreasson et al.](https://arxiv.org/pdf/1811.12338)
- **Sivak et al. (Nature 2023).** **Verified**: Nature 616, 50–55. — [IDEAS/RePEc](https://ideas.repec.org/a/nat/nature/v616y2023i7955d10.1038_s41586-023-05782-6.html)
- **Gidney, Shutty, Jones (2024), magic-state cultivation.** **Not re-verified** in this session; only its 2025–2026 follow-ups were retrieved (§5).

### Inferences
- **The shared design is "LLM proposes, exact checker disposes".** The checker is exhaustive fault enumeration (Lahiri & Marquardt), an ILP plus independent verification (factories), or DEM sampling plus BP+OSD (OmniQEC). This keeps hallucinations out of the final results. The cost then sits in verification, which is why QSA-type harnesses (coresets, approximate scoring) and importance sampling matter.
- **The strongest LLM results so far are in Clifford-verifiable settings** (state preparation, Clifford gadgets, distillation factories whose correctness is algebraic). LLM search in settings without cheap exact verifiers, such as decoders under realistic noise, non-Clifford cultivation circuits and hardware calibration, is essentially absent.
- **Reward hacking and adversarial advice are already documented** (qBraid's distance-only objective; Özgüler's 12/12 harmful choices under attack). Agentic QEC pipelines therefore need multi-objective, certified acceptance tests.

### Gaps
- Exact per-code numbers for Lahiri & Marquardt and for OmniQEC's LER, plus the full author list of arXiv:2610.06535, were n/r because arXiv full text was unreachable.
- Despite targeted queries, no LLM-guided work was found for decoders, lattice surgery, code switching, Floquet-code scheduling or cultivation.

---

## 2. Fault-tolerant state preparation and encoders

### Takeaway
RL remains the main ML result: Zen et al. (PRX 2025) found flag-verified FT state-preparation circuits with fewer gates and flags than human designs, up to distance-5 codes. In 2025–2026 three lines extended it: gadget abstraction for RL scaling (Marquardt group, encoders up to distance 7), LLM program synthesis that covers d = 3–9 in a single search (Nov 2026), and multi-agent RL that discovers whole QEC cycles. Exact non-ML methods set strong baselines and in places prove optimality: MQT SAT synthesis, the Block-Check Pareto tree search (d ≤ 7), and ZX-based SpiderCat (provably CNOT-optimal FT cat states for n ≤ 100, t ≤ 5).

### Cited Findings
- **Schmid & Wille (Nov 2026), arXiv:2611.25107.** Non-ML. Introduces "Block Checks", in which several checks share one fault-tolerantly prepared auxiliary block, placing earlier verification constructions in one design space. A shared auxiliary must itself be prepared fault-tolerantly, a nested preparation problem, which they explore with a Pareto tree search. On codes up to distance 7 it finds "best-known circuits for many new operation–depth trade-offs" and surpasses existing methods. — [arXiv:2611.25107](https://arxiv.org/abs/2611.25107)
- **Lahiri & Marquardt (Nov 2026).** The LLM evolutionary synthesis in §1 treats FT logical state preparation as the first test task, "since fault-tolerant state preparation … has long been manually obtained". — [arXiv HTML](https://arxiv.org/html/2611.18504)
- **SpiderCSS (Oct 2026), arXiv:2610.03714.** Builds on SpiderCat for general CSS state preparation. [title-level] — [arXiv:2610.03714](https://arxiv.org/pdf/2610.03714); [SpiderCat summary](https://arxiv.org/abs/2603.05391)
- **SpiderCat: Khesin, Li, Poór, Rodatz, van de Wetering, Yeung (Mar 2026), arXiv:2603.05391.**
  - Notes that earlier FT cat-state constructions relied on "SAT solving, reinforcement learning, or exhaustive analysis".
  - Uses fault-equivalent ZX rewrites to reduce the problem to characterising certain 3-regular simple graphs.
  - Gives provably CNOT-optimal *t*-FT *n*-qubit cat states for n ≤ 100, t ≤ 5, plus a constant-depth variant that uses more ancillas.
  - The prior SAT method was limited to n ≤ 47 at t = 5 (n ≤ 19 at t = 9).
  - Sources: [arXiv:2603.05391](https://arxiv.org/abs/2603.05391); [project page](https://boldar99.github.io/spidercat/)
- **"Optimizing Fault-tolerant Cat State Preparation", arXiv:2601.03343; Phys. Rev. A (DOI 10.1103/13v7-n843).** [title-only; authors n/r] — [arXiv HTML](https://arxiv.org/html/2601.03343v1); [PRA](https://journals.aps.org/pra/pdf/10.1103/13v7-n843)
- **Weilandt, Peham, Wille (Jan 2026), arXiv:2601.13313.** Extends SAT/automated synthesis to Steane-type FT preparation for arbitrary CSS codes. [details n/r] — [arXiv:2601.13313](https://arxiv.org/pdf/2601.13313)
- **Zen, Olle, Colmenarez, Puviani, Müller, Marquardt, "Quantum circuit discovery for fault-tolerant logical state preparation with reinforcement learning", arXiv:2402.17761 (Feb 2024) → Phys. Rev. X 15, 041012 (22 Oct 2025).**
  - *Agent:* RL with a fast, parallelized stabilizer simulator and a reward function tailored to the task. It searches for flag-based FT preparation and verification circuits.
  - *Scale:* fewer gates and ancillas than published results, with and without connectivity constraints, for up to 15 physical qubits (arXiv v1). The published version extends to distance-5 codes. Transfer learning speeds up discovery.
  - *Example:* where a human verification circuit uses 6 flags (or 2 flags plus 2 resets) and 15 two-qubit gates, the RL circuit uses 2 flags and 7 two-qubit gates.
  - *Stated extensions (not shown):* magic-state preparation, logical-gate synthesis and syndrome measurement.
  - *Code:* GitHub "rlftqc".
  - *Earlier version:* a NeurIPS ML4PS 2023 workshop paper (without Colmenarez and Müller) reported codes up to 17 qubits with all-to-all connectivity.
  - Sources: [arXiv:2402.17761](https://arxiv.org/abs/2402.17761); [PRX accepted page](https://journals.aps.org/prx/accepted/10.1103/gqpr-dgz7); [ADS](https://ui.adsabs.harvard.edu/abs/2025PhRvX..15d1012Z/abstract); [GitHub](https://www.Github.com/remmyzen/rlftqc); [ML4PS 2023](https://ml4physicalsciences.github.io/2023/files/NeurIPS_ML4PS_2023_84.pdf)
- **Yevtushenko & Marquardt (Sep 2025), "Automated Discovery of Gadgets in Quantum Circuits for Efficient Reinforcement Learning", arXiv:2509.24666.**
  - Represents circuits as directed graphs and searches automatically for repeated subgraphs, giving gadgets and gadget families.
  - The goal is compact, hardware-adapted circuits that fault-tolerantly prepare logical QEC states with RL.
  - With the "DCX4" and "PL4" gadgets, encoders were found with a 100% success rate in their runs.
  - Sources: [arXiv HTML](https://arxiv.org/html/2509.24666); [MPG PuRe](https://pure.mpg.de/pubman/faces/ViewItemOverviewPage.jsp?itemId=item_3672554_1)
- **Guatto, Preti, Schilling, Calarco, Cárdenas-López, Motzoi (Sep 2025), arXiv:2509.03974 → Quantum Sci. Technol. (DOI 10.1088/2058-9565/aebdc4).**
  - Multi-agent RL discovers complete QEC cycles with explicit circuits for encoding, syndrome extraction and recovery.
  - An online "BRAVE" bandit decides when to retrain.
  - Reported: about 18× lower logical infidelity for qubit codes and about 3× for qutrit codes under time-varying noise [secondary].
  - Sources: [IOP](https://iopscience.iop.org/article/10.1088/2058-9565/aebdc4); [arXiv HTML](https://arxiv.org/html/2509.03974); [GitHub](https://github.com/ManuelGuatto/Adaptive_QEC); [Pith summary](https://pith.science/paper/2509.03974) [secondary]
- **Olle, Yevtushenko, Marquardt (Mar 2025), "Scaling the Automated Discovery of Quantum Circuits via Reinforcement Learning with Gadgets", arXiv:2503.11638.** Uses gadgets as macro-actions to scale RL. The case study is Clifford encoders for QEC, and it reports finding distance-7 codes (with encoders) from scratch. — [arXiv:2503.11638](https://arxiv.org/abs/2503.11638); [summary in 2509.24666](https://arxiv.org/html/2509.24666)
- **Schmid, Peham, Berent, Müller, Wille (DATE 2025), "Deterministic Fault-Tolerant State Preparation for Near-Term Quantum Error Correction: Automatic Synthesis Using Boolean Satisfiability", arXiv:2501.05527.** Applies SAT synthesis to deterministic FT preparation. — [DATE 2025 PDF](https://www.cda.cit.tum.de/files/eda/2025_date_deterministic_state_preparation.pdf); [arXiv](https://export.arxiv.org/pdf/2501.05527)
- **Peham, Schmid, Berent, Müller, Wille, "Automated Synthesis of Fault-Tolerant State Preparation Circuits for Quantum Error Correction Codes", arXiv:2408.11894 (2024) → PRX Quantum 6, 020330 (2025).** SAT solving gives depth- and gate-optimal preparation and verification circuits, with heuristics for larger codes. It is motivated by the absence of earlier algorithmic methods, since prior circuits were mostly hand-made. Open source in MQT QECC. — [arXiv:2408.11894v3](https://arxiv.org/abs/2408.11894v3); [TUM portal](https://portal.fis.tum.de/en/publications/automated-synthesis-of-fault-tolerant-state-preparation-circuits-/)
- **General RL Clifford/CNOT synthesis, relevant to encoders but general compilation and out of scope.** IBM's "Practical and efficient quantum circuit synthesis and transpiling with RL" (arXiv:2405.13196), "Pauli Network Circuit Synthesis with RL" (arXiv:2503.14448), "CNOT Minimal Circuit Synthesis: A RL Approach" (arXiv:2510.23304) and "Equivariant RL for Clifford Quantum Circuit Synthesis" (arXiv:2605.10910). [title-only] — [2405.13196](https://arxiv.org/pdf/2405.13196); [2503.14448](https://arxiv.org/pdf/2503.14448); [2510.23304](https://arxiv.org/pdf/2510.23304); [2605.10910](https://arxiv.org/pdf/2605.10910)

### Inferences
- **Scale ceilings by method:**
  - RL on individual instances: about d = 5 for FT state preparation (Zen et al.) and d = 7 for encoders with gadgets (Olle et al.).
  - LLM program synthesis: produces *generators of families*, with d = 3–9 in one search (Lahiri & Marquardt).
  - Exact methods: SpiderCat reaches n = 100 for cat states.
  The program-synthesis route looks most promising for going beyond d ≈ 9, but its exhaustive verifier grows combinatorially with *t*.
- **Complementary strengths.** Exact methods (SAT, ZX) give optimality certificates that RL and LLM search lack, while RL and LLMs cope more easily with irregular connectivity constraints and cost functions. The best current practice is to benchmark ML-found circuits against MQT SAT outputs and SpiderCat bounds.

### Gaps
- No hardware run of an RL- or LLM-discovered FT state-preparation circuit was found in the retrieved sources.
- No ML work was found on FT encoders or state preparation for large qLDPC codes such as [[144,12,12]].
- Content of IBM's RL synthesis paper was n/r. Olle et al., npj QI 2024 ("Simultaneous discovery of QEC codes and encoders with a noise-aware RL agent") is recalled from memory and **unverified** in this session.

---

## 3. Syndrome-extraction circuits: scheduling, hook errors, flags, qLDPC and Floquet codes

### Takeaway
"AlphaSyndrome" **exists**: an MCTS scheduler for commuting-stabilizer codes (arXiv:2601.12509, Jan 2026). In 2026 the topic went from almost no automated methods to at least seven, in three families:
- **Search with a decoder in the loop:** AlphaSyndrome (MCTS) and FastSched (RL with importance sampling).
- **Search over a proxy evaluated by a solver:** PropHunt (MaxSAT ambiguity) and Proxima (SMT "proximity").
- **Structured constructions:** left-right circuits; structure-guided synthesis.

Each paper reports beating the previous ones on its own benchmark: roughly 1.4–5× LER over automated baselines, 2.1–5.1× over IBM's BB schedule, and up to 10× over single-ancilla designs. The results are not directly comparable. No ML or automated work on Floquet-code measurement scheduling was found.

### Cited Findings

**Automated syndrome-extraction schedulers (2026, newest first)**

| Work | Search / learning | Reward or verifier | Codes tested | Headline claim | Source |
|---|---|---|---|---|---|
| Structure-guided (Li et al., Nov 2026) | groups checks by connectivity; one shared CNOT order per group | LER estimates | 15 qLDPC codes, 10 families | lowest LER vs ASC, PropHunt, AlphaSyndrome; 21.4–136.6× faster to first feasible circuit than ASC | [2611.25825](https://arxiv.org/abs/2611.25825) |
| Proxima (Zhang et al., Nov 2026) | solver search over a proxy | SMT on local DEM column subsets ("proximity", an upper bound on circuit distance) | BB, RQT and others | 1.4–2.6× LER vs best automated baseline; up to 40× vs others; BB depth-8 beats IBM depth-7 by 2.1–5.1× | [2611.19715](https://arxiv.org/abs/2611.19715) |
| FastSched (Ye, Pabla, Palsberg, Sep 2026) | RL policy picks the CNOT order of each check | LER by importance sampling | surface codes d = 3–15, others | −25.9% LER vs AlphaSyndrome, −71.7% vs PropHunt (averages); −97.8% at d = 15 surface code | [2609.12020](https://arxiv.org/abs/2609.12020) |
| Left-right circuits (Strikis, Browne, Beverland, Mar 2026) | generalized LRC construction plus optimization | residual-error distance metrics, idling, minimum-weight failures | [[61,1,9]] colour, [[108,18,4]] RQT, others | up to 10× vs existing single-ancilla designs | [2603.05481](https://arxiv.org/abs/2603.05481) |
| PropHunt (Viszlai et al., Jan 2026) | iterative CNOT reordering from a coloration circuit | MaxSAT finds low-weight decoder-ambiguous faults | surface, lifted-product, RQT | matches hand-made surface schedules; 2.5–4× vs coloration circuits on LP/RQT | [2601.17580](https://arxiv.org/html/2601.17580v1) |
| AlphaSyndrome (Liu et al., Jan 2026) | MCTS over commutation-preserving reorderings | decoded LER (sampling plus a heuristic decoder in every rollout) | surface, [[72,12,6]] BB, others | 80.6% average, up to 96.2%, LER reduction vs baselines [secondary] | [2601.12509](https://arxiv.org/abs/2601.12509) |

- **AlphaSyndrome (Liu et al., UPenn and UC Berkeley; Jan 2026), arXiv:2601.12509.**
  - *Method:* MCTS whose search tree starts from an initial state S₀. Each node stores an accumulated score and a visit count, and the search space is limited to commutation-preserving reorderings.
  - *Assumptions:* only mutually commuting stabilizers and a heuristic decoder. Different schedules create different error-propagation paths and so very different LERs: a trivial d = 3 surface-code schedule can be 9× worse than Google's good schedule.
  - *Benchmarks:* includes a [[72,12,6]] BB code against IBM's schedule. Two-qubit-gate depth is sometimes higher, but this "does not significantly increase—and in some cases even reduces—the absolute execution time" needed to reach a given LER.
  - *Headline numbers (secondary):* an average of 80.6% and up to 96.2% LER reduction vs standard methods.
  - *Cost and comparisons (later papers):* AlphaSyndrome samples and decodes at every rollout, and both AlphaSyndrome and PropHunt "run multi-hour iterative searches". Its schedules reportedly match Google's hand-made surface-code schedules.
  - Sources: [arXiv HTML](https://arxiv.org/html/2601.12509); [ResearchGate](https://www.researchgate.net/publication/399931521_AlphaSyndrome_Tackling_the_Syndrome_Measurement_Circuit_Scheduling_Problem_for_QEC_Codes); [Quantum Zeitgeist](https://quantumzeitgeist.com/percent-error-alphasyndrome-achieves-logical-rate-reduction/) [secondary]; later-paper statements via [Proxima](https://arxiv.org/html/2611.19715) and [Structure-guided](https://arxiv.org/html/2611.25825) (which paper made each statement was not pinned down)
- **FastSched (Ye, Pabla, Palsberg; Sep 2026), "Reinforcement Learning for Syndrome Extraction", arXiv:2609.12020 (cs.LG, cross-listed in quant-ph).**
  - Frames the problem as a search over implementations whose number grows exponentially with syndrome size.
  - Combines an RL schedule search with importance-sampled LER evaluation, which saves decoded shots.
  - Numbers as in the table; they are the authors' own benchmarks.
  - Sources: [arXiv:2609.12020](https://arxiv.org/abs/2609.12020); [Proxima related work](https://arxiv.org/html/2611.19715)
- **PropHunt (Viszlai et al.; Jan 2026), arXiv:2601.17580.**
  - *Method:* targets *decoding ambiguity* rather than simulating the full LER. MaxSAT finds low-weight fault patterns that a decoder could confuse, and the circuit is then rewired or reordered.
  - *Limitations:* some searches found no improved schedule within their time caps. It targets CSS codes, whereas AlphaSyndrome handles any commuting-stabilizer code.
  - *Proxima's critique:* PropHunt starts from a coloration circuit and estimates ambiguity on randomly grown DEM subgraphs. On RQT108 it needs 30 layers and has 22× higher LER than Proxima.
  - Sources: [arXiv HTML](https://arxiv.org/html/2601.17580v1); [INSPIRE](https://inspirehep.net/literature/3111511); [Proxima](https://arxiv.org/html/2611.19715)
- **Proxima (Ruiyun Zhang, Guangqi Zhao, Zhaohui Yang, Yingjian Liu, Chenghong Zhu, Anbang Wu, Jun Xia, Yuan Xie; Nov 2026), arXiv:2611.19715.** States that for general CSS qLDPC codes, directly optimizing LER or exact circuit distance is "computationally prohibitive", which motivates the proximity proxy. Its BB-code comparison with IBM's depth-7 schedule is worded ambiguously in the abstract. — [arXiv:2611.19715](https://arxiv.org/abs/2611.19715)
- **Structure-guided synthesis (Shanchuan Li, Xin Hong, Pengcheng Zhu, Xiangzhen Zhou, Fanxu Meng, Shenggang Ying; 19 Nov 2026), arXiv:2611.25825.** "ASC" is a prior automated syndrome-circuit method whose identity was not confirmed. — [arXiv:2611.25825](https://arxiv.org/abs/2611.25825)
- **Strikis, Browne, Beverland (IBM; 5 Mar 2026), "High-performance syndrome extraction circuits for quantum codes", arXiv:2603.05481.**
  - Generalizes left-right circuits, which offset X and Z checks so their gates do not interleave, to arbitrary CSS codes.
  - Optimizes for low idling, large effective distance and fewer minimum-weight failure mechanisms, using "residual errors" as lightweight error-propagation metrics.
  - Code is on GitHub and circuits on Zenodo.
  - Source: [arXiv:2603.05481](https://arxiv.org/abs/2603.05481)
- **Verification by simulation in qLDPC architecture papers.** "High-rate qLDPC processors" (arXiv:2607.28795) compiles candidate orderings with Stim and discards schedules with low estimated distance. "Untangling QLDPC Codes with Biased Noise Ancilla" (arXiv:2606.30592) avoids hooks and loops with biased-noise ancillas. Neither is ML. — [2607.28795](https://arxiv.org/pdf/2607.28795); [2606.30592](https://arxiv.org/pdf/2606.30592)
- **Other non-ML SE compilers [title-only]:**
  - "Syndrome Extraction Circuits with Near-Optimal Depths for Practical QEC Code Families" (DAC '26; targets lifted-product codes) — [DOI](https://doi.org/10.1145/3770743.3803913)
  - "Disassembling qLDPC codes for depth-optimal parity-check circuits" — [2608.19917](https://arxiv.org/pdf/2608.19917)
  - "Local Automorphism-Aware Syndrome Compilation for General QLDPC Codes" — [2609.40319](https://arxiv.org/pdf/2609.40319)
  - "Scheduling of syndrome measurements with a few ancillary qubits" — [2508.07913](https://arxiv.org/pdf/2508.07913)
  - "Fault-Tolerant Cut-Cat State Syndrome Extraction" — [2604.17339](https://arxiv.org/pdf/2604.17339)
  - "Turbo-Annihilation of Hook Errors" — [2504.21200](https://arxiv.org/pdf/2504.21200)
  - "Tangling schedules eases hardware connectivity requirements" — [ResearchGate](https://www.researchgate.net/publication/379128916_Tangling_Schedules_Eases_Hardware_Connectivity_Requirements_for_Quantum_Error_Correction)
- **Flags.**
  - The RL flag-verification results of Zen et al. are in §2. They list syndrome measurement as a future direction only. — [PRX accepted page](https://journals.aps.org/prx/accepted/10.1103/gqpr-dgz7)
  - Flag Proxy Networks (arXiv:2409.14283) addresses architecture, scheduling and decoding for qLDPC with flags (non-ML). [title-only] — [2409.14283](https://arxiv.org/pdf/2409.14283)
  - A flag SE circuit was evaluated on IBM hardware (arXiv:2403.10217). [title-only] — [2403.10217](https://arxiv.org/pdf/2403.10217)

### Inferences
- **Reward design trades accuracy against cost.** Rewards based on directly sampled LER (AlphaSyndrome, FastSched, OmniQEC) are accurate but expensive, noisy and tied to the decoder used in the loop. Solver proxies (PropHunt ambiguity, Proxima proximity, LRC residual distance) are cheaper and decoder-agnostic but may miss entropic effects. The competing 2026 claims of "X% better than AlphaSyndrome/PropHunt" come from different code suites, noise models and decoders.
- **All results are in simulation**, with circuit-level depolarizing or SI1000-like noise. None of the retrieved schedulers was validated on hardware.

### Gaps
- No ML or search work on measurement scheduling for **Floquet or dynamic codes** was found; this is also partly a search-budget limitation.
- A targeted search for a 2025 work combining flags, SAT, RL and syndrome measurement found none (negative search result, no source to cite).
- No ML-designed **flag gadgets for qLDPC codes** were found.
- It was not possible to verify which decoder each scheduler used inside its loop, or whether schedules transfer between decoders.

---

## 4. Logical operations: FT gadgets, lattice surgery, transversal/fold gates, code switching

### Takeaway
ML here is thin. The brief's "Trenkwalder–Briegel automated gadget discovery" is real but concerns mining reusable sub-sequences from RL agents in quantum-optics and circuit tasks, not FT gadgets. The FT-relevant lines are the Marquardt group's gadget RL (2025) and LLM evolutionary synthesis of logical gates (Nov 2026, including a shallower logical Hadamard). Lattice-surgery compilation in 2024–2026 uses SAT, ZX and graph routing, with no ML found. No ML/LLM work on code switching was found, and fold-transversal/automorphism gate search could not be checked before the search budget ran out.

### Cited Findings
- **Lahiri & Marquardt (Nov 2026).** Targets include logical gates and logical routines, and the authors report shallower logical-Hadamard circuits than a prior reference. Every circuit is checked by an exhaustive FT verifier. — [arXiv:2611.18504](https://arxiv.org/pdf/2611.18504)
- **Trenkwalder, López-Incera, Poulsen Nautrup, Flamini, Briegel (MLST 2023; arXiv:2212.12743, Dec 2022 — earlier work).** A pipeline of RL, then sequence mining of the agent's trajectories, then density-based clustering, which distils frequent compact subroutines ("gadgets"). Applied to quantum optics, where the gadgets match interferometer setups, and to a circuit environment with gadgets such as teleportation. — [DOAJ](https://doaj.org/article/61bf0003a8ea4bb69c5f1e40ea217406); [arXiv:2212.12743](https://arxiv.org/pdf/2212.12743)
- **Yevtushenko & Marquardt (2025) and Olle, Yevtushenko, Marquardt (2025).** Gadget discovery by repeated-subgraph search, used as RL macro-actions for QEC encoders and FT state preparation (see §2). — [2509.24666](https://arxiv.org/html/2509.24666); [2503.11638](https://arxiv.org/abs/2503.11638)
- **"Reinforcement learning with learned gadgets to tackle hard quantum problems on real hardware", arXiv:2411.00230.** [title-only] — [2411.00230](https://arxiv.org/pdf/2411.00230)
- **Non-ML lattice-surgery baselines.**
  - Tan, Niu, Gidney, "A SAT scalpel for lattice surgery" (arXiv:2404.18369, ISCA 2024; ID from search snippet) synthesizes lattice-surgery subroutines with SAT. — [arXiv:2404.18369](https://arxiv.org/abs/2404.18369)
  - "Lattice Surgery Compilation Beyond the Surface Code" (TUM; arXiv:2504.10591; IEEE QCE 2025) recasts compilation as mapping and routing on a macroscopic routing graph. — [arXiv:2504.10591](https://arxiv.org/abs/2504.10591); [QCE PDF](https://www.cda.cit.tum.de/files/eda/2025_qce_lattice_surgery_compilation_beyond_the_surface_code.pdf)
  - Other 2025–2026 compilers [title-only]: SpiderLS (full ZX reduction) [2608.30228](https://arxiv.org/pdf/2608.30228); flow-based optimization with runtime T scheduling [2609.23756](https://arxiv.org/pdf/2609.23756); demand-driven routing space [2610.30747](https://arxiv.org/pdf/2610.30747); resource-allocating compiler [2506.04620](https://www.alphaxiv.org/abs/2506.04620); movable logical qubits [2512.04169](https://arxiv.org/pdf/2512.04169); CUDA-Q Logical [2609.13388](https://arxiv.org/pdf/2609.13388).
- **Hardware lattice surgery (non-ML context).**
  - Lattice surgery on two d = 3 repetition codes with superconducting qubits (Nature Physics); it is fault-tolerant to bit flips but not phase flips. — [Nature Physics](https://www.nature.com/articles/s41567-025-03090-6)
  - Teleportation of logical states between two d = 3 logical qubits on 41 physical qubits (Nov 2026). — [arXiv:2611.23165](https://arxiv.org/abs/2611.23165)

### Inferences
- In logical operations, automated search is dominated by SAT, ZX and routing formulations. ML has so far contributed only Clifford-gadget discovery (RL and LLM) on small codes. Learning-based lattice-surgery compilation, code switching and qLDPC logical-gate discovery look like open territory, since targeted queries returned nothing.

### Gaps
- Fold-transversal and automorphism gate search for qLDPC codes was not searched (budget exhausted). Neither was ML-driven code switching, except one query that found nothing.
- A targeted search for RL or ML lattice-surgery compilation returned only non-ML compilers (negative search result).
- Trenkwalder et al.'s quantitative results were n/r.

---

## 5. Magic states: distillation-factory discovery and cultivation

### Takeaway
The first **AI-assisted discovery of distillation protocols** appeared in Oct 2026: Claude agents steer symmetry choices, an ILP builds candidates and an independent verifier checks them, giving 564 new factory classes. It sits beside non-ML SAT/branch-and-bound and algebraic searches. Cultivation follow-ups are all human, analytic or solver-driven, and **no ML/LLM-driven cultivation design was found**. On hardware, Google demonstrated cultivation on Willow with fidelity ≥ 0.9999(1) at 8% retention and about 40× lower error, and an NVIDIA–QuEra transformer decoder beat MLE decoding on QuEra's d = 3 distillation circuit.

### Cited Findings
- **AI-assisted factories (arXiv:2610.06535).** Details in §1. — [arXiv:2610.06535](https://arxiv.org/abs/2610.06535)
- **Non-ML factory searches (2026).**
  - "Exploring the landscape of compact magic-state distillation factories" (arXiv:2606.07734) uses SAT search and branch-and-bound with ILP certification, and notes that local search such as simulated annealing tends to be weaker for this problem. — [2606.07734](https://arxiv.org/pdf/2606.07734)
  - "The limits of magic-state distillation" (arXiv:2610.31115) uses integer and SAT programs to derive lower bounds. — [2610.31115](https://arxiv.org/pdf/2610.31115)
  - "Magic State Distillation via Codes over Binary Extension Fields" (arXiv:2608.09727) builds protocols with algebraic geometry and claims they beat the state of the art "in almost every situation" in its case study. — [2608.09727](https://arxiv.org/abs/2608.09727v1)
  - "Invariant Theory, Magic State Distillation, and Bounds on Classical Codes" (arXiv:2501.10163). [title-only] — [2501.10163](https://arxiv.org/pdf/2501.10163)
- **Cultivation follow-ups (non-ML).**
  - "Diagnosing and Restoring the Degraded Fault Distance of Magic State Cultivation" (arXiv:2609.17706). A companion paper credits Tim Chan with first reporting the reduced fault distance. — [2609.17706](https://arxiv.org/pdf/2609.17706)
  - "Exact logical error rates for magic state cultivation" (arXiv:2609.18922). — [2609.18922](https://arxiv.org/pdf/2609.18922)
  - "MagiCFirm", an algorithm–hardware co-design runtime (arXiv:2609.29267). — [2609.29267](https://arxiv.org/pdf/2609.29267)
  - Cultivation on RP² (arXiv:2503.18657), "nearly an order of magnitude reduction in space-time cost". — [alphaXiv](https://www.alphaxiv.org/abs/2503.18657)
  - Variants for lattice surgery, two-qubit-gate-only surface codes, constant depth via gauging, and fold-transversal surface codes are listed in surveys. — [Emergent Mind topic](https://www.emergentmind.com/topics/magic-state-cultivation-protocols) [secondary]
  - Logical catalysts for dyadic phase rotations (arXiv:2606.27358). [title-only] — [2606.27358](https://arxiv.org/pdf/2606.27358)
- **Non-Clifford verifiers.** Wan, Zhong, Zapirain, "Simulating magic state cultivation with few Clifford terms" (arXiv:2509.08658; Quantum, 12 Jun 2026) simulates d = 5 cultivation with an average of about 8 Clifford ZX-diagrams at 0.1% noise. Related tools: SOFT, Clifft, tsim. "Efficient simulation of logical magic state preparation protocols" (arXiv:2512.23799) is title-only. — [arXiv:2509.08658](https://arxiv.org/abs/2509.08658); [Quantum](https://quantum-journal.org/papers/q-2026-06-12-2134/); [2512.23799](https://arxiv.org/pdf/2512.23799)
- **Google cultivation on Willow (hardware; preprint Dec 2025, *Science* 2026).**
  - State fidelity ≥ 0.9999(1) while keeping 8% of attempts, and about 40× lower error than without cultivation.
  - This is the smallest (distance-3) version. Grafting into a surface code was a proof of principle with no grafted-state fidelity reported, and Google says larger-scale cultivation is still to be demonstrated.
  - Sources disagree on timing: preprint Dec 2025; Science announcement dated 15 Oct 2026.
  - Sources: [Quantum Computing Report](https://quantumcomputingreport.com/google-quantum-ai-demonstrates-magic-state-cultivation-on-willow-processor-in-science-milestone/) [secondary]; [PostQuantum](https://postquantum.com/quantum-research/google-magic-state-cultivation/) [secondary]; [postquantumsecurity.org](https://www.postquantumsecurity.org/publications/google_willow.html) [secondary]
- **ML for decoding and post-selection in magic-state protocols.**
  - NVIDIA–QuEra transformer decoder (GTC 2025; CUDA-Q, PhysicsNeMo). It outperformed the most-likely-error decoder on QuEra's d = 3 distillation circuit; a 50× speed-up figure is [secondary]. Scaling to higher distance needs AI-supercomputer-generated training data. — [NVIDIA developer blog](https://developer.nvidia.com/blog/nvidia-and-quera-decode-quantum-errors-with-ai); [Quantum Computing Report](https://quantumcomputingreport.com/nvidia-and-quera-advance-ai-driven-quantum-error-decoding-for-scalable-fault-tolerant-quantum-computing/) [secondary]
  - "Machine-learned syndrome post-selection for reliable QEC" (arXiv:2607.19563). [title-only] — [2607.19563](https://arxiv.org/pdf/2607.19563)
  - "Graph Neural Post-selection for QEC" (arXiv:2610.00504). [title-only] — [2610.00504](https://arxiv.org/pdf/2610.00504)
- **Stated intent only.** Zen et al. list magic-state preparation as a future extension of their RL framework. — [PRX accepted page](https://journals.aps.org/prx/accepted/10.1103/gqpr-dgz7)

### Inferences
- **Why LLMs succeeded on factories first.** Distillation-factory design is a near-ideal target for LLM-plus-solver search: correctness (gate identities, distances) can be checked exactly and cheaply, and the search space has algebraic structure that LLMs can be prompted to exploit.
- **Why cultivation lags.** It needs non-Clifford simulation (few-Clifford-term methods) and rare-event LER estimation, which makes rewards expensive. Its known flaw, the degraded fault distance, was found by a human.
- **Learned post-selection** (two 2026 titles) is a natural ML entry point to cultivation, since cultivation's cost is set by its discard rate (8% retention on Willow).

### Gaps
- Gidney–Shutty–Jones (2024) itself was not re-verified here.
- No ML/LLM search over cultivation circuits was found.
- No hardware test of AI-discovered factories, and no circuit-level noise evaluation of the 564 new classes, was found.
- Contents of the two ML post-selection papers were n/r.

---

## 6. Learned decoders: AlphaQubit and successors, transformer/GNN/CNN/RNN, neural BP, RL decoders, hardware-trained decoders, learned noise models, and accuracy vs latency

### Takeaway
**Surface and colour codes.** Neural decoders lead on accuracy: AlphaQubit made 6% fewer errors than tensor-network decoding and 30% fewer than correlated matching on Sycamore. Since Dec 2025 they also reach real-time throughput: AlphaQubit 2 runs under 1 µs per cycle up to d = 11 (surface) and d = 9 (colour); a self-coordinating parallel decoder claims d = 25 on one TPU; a closed-loop FPGA decoder ran on hardware in 550 ns at d = 3; and Rapids claims a 240× lower LER than earlier real-time hardware decoders.

**qLDPC codes.** Neural decoding is less mature. A transformer is about 5× better than BP-OSD on [[72,12,6]] but degrades on [[144,12,12]]. Learned message-passing and neural-BP decoders claim parity with Relay-BP or BP-OSD, mostly in code-capacity or simulated settings.

**RL decoders** are essentially a pre-2023 line.

### Cited Findings

**6a. Google: AlphaQubit (AQ1), Willow, AlphaQubit 2 (AQ2)**
- **Bausch et al., "Learning high-accuracy error decoding for quantum processors", Nature (Nov 2024).**
  - Recurrent, transformer-based decoder.
  - On Sycamore d = 3 and d = 5 surface-code data: 6% fewer errors than tensor-network decoding ("highly accurate but impractically slow") and 30% fewer than correlated matching.
  - Keeps its advantage on simulated data up to d = 11 with realistic noise including cross-talk and leakage.
  - *Limitation:* "still too slow to correct errors in a superconducting processor in real time".
  - Gigazine reverses the 6% and 30% figures; Google's blog supports the version given here.
  - Sources: [Google blog](https://blog.google/technology/google-deepmind/alphaqubit-quantum-error-correction/); [Gigazine](https://gigazine.net/gsc_news/en/20241121-alphaqubit-google-quantum-error-correction-ai) (conflicting)
- **Google Quantum AI, "Quantum error correction below the surface code threshold", Nature (2025), DOI 10.1038/s41586-024-08449-y; arXiv:2408.13687.**
  - Λ = 2.14 ± 0.02, reaching 0.143% error per cycle at d = 7.
  - Real-time decoding at d = 5: average latency about 63 ± 17 µs, roughly constant for runs up to 1.1 s, with a 1.1 µs cycle. Real-time accuracy is below offline accuracy; a snippet gives Λ = 2.0 ± 0.1 with a prior that does not use device data.
  - Commentary notes that the 63 µs spans more than 50 cycles, so "bounded latency" may be the more accurate term.
  - The individual neural-network and ensemble-matching figures were n/r.
  - Sources: [Nature PDF](https://www.nature.com/articles/s41586-024-08449-y.pdf); [arXiv:2408.13687](https://arxiv.org/pdf/2408.13687); [notes by dyb](https://dyb.mataroa.blog/blog/assessment-of-quantum-error-correction-below-the-surface-code-threshold/)
- **Senior et al. (Google DeepMind and Google Quantum AI; Dec 2025), "A scalable and real-time neural decoder for topological quantum codes" (AQ2), arXiv:2512.07737.**
  - Claims near-optimal LER for surface and colour codes at scale under realistic noise.
  - The real-time variant AQ2-RT runs under 1 µs per cycle on commercial accelerators, for the surface code up to d = 11 and the colour code up to d = 9. For the colour code it is "orders of magnitude faster than other high-accuracy decoders", and uses a "superdense" variant in real time. Chromobius is faster but less accurate and reaches real-time throughput only at d = 3.
  - On Willow data, AQ2-RT beats the original Willow real-time decoder with a small accuracy loss relative to full AQ2. AQ2 is "comparable to … AQ1 and better than Libra on the same experimental dataset"; a Lacuna summary instead says it beat AQ1 [conflict].
  - A Dec 2025 preprint citing AQ2 says it is about 500× below Chromobius at d = 15 (SI1000, p = 0.15%) and notes AQ2-RT latencies up to a millisecond [secondary]. Throughput (µs per cycle) and end-to-end latency should not be conflated.
  - Sources: [arXiv:2512.07737](https://arxiv.org/abs/2512.07737); [HTML](https://arxiv.org/html/2512.07737v1); [Lacuna](https://lacuna.tiptreesystems.com/work/a-scalable-and-real-time-neural-decoder-for-topological-quantum-codes/wrk_8cb5ba0c7e29f08e681896545cce4502) [secondary]; [Emergent Mind](https://www.emergentmind.com/topics/alphaqubit-2) [secondary]
- **Manor, Erhili, Jebbouri (Sep 2026), "A Sim-to-Real Study of Surface-Code Decoder Benchmarking", arXiv:2609.04557.**
  - Ranks six decoders under four noise models of increasing fidelity, scored against real Willow data (3 distances, 2 bases, 15 round counts).
  - An "Ising" pre-decoder was the only one whose error rate rose with distance on Willow, and it ranked last.
  - A SciRate comment objects that the Ising decoders ran on CPU without retraining.
  - A Substack summary says a template built without Willow data predicted the hardware ranking perfectly [secondary, unconfirmed].
  - Sources: [arXiv:2609.04557](https://arxiv.org/abs/2609.04557); [SciRate](https://scirate.com/arxiv/2609.04557); [Substack](https://christophergannatti.substack.com/p/grading-on-a-curve) [secondary]

**6b. Other neural decoders: throughput, latency and hardware**
- **Rapids (Jin, Yelin, Gu; Nov 2026), arXiv:2611.26447.** FPGA co-designed, integer-trained streaming decoder. Claims 240× lower LER than existing real-time hardware decoders, and 16× better than the real-time Relay decoder on [[144,12,12]] BB. — [arXiv:2611.26447](https://arxiv.org/abs/2611.26447)
- **PRISM (Nov 2026), arXiv:2611.13630.** Erasure-aware transformer for photonic fusion-based QC. Matches integer-programming decoders up to distance 24 on surface, colour and qLDPC codes, and is significantly faster. — [arXiv:2611.13630](https://arxiv.org/abs/2611.13630v1)
- **Kai Zhang et al. (Jan 2026), "Learning to Decode in Parallel: Self-Coordinating Neural Network for Real-Time QEC", arXiv:2601.09921 (earlier version arXiv:2509.03815).**
  - A recurrent transformer for parallel-window decoding. It outputs a combinable logical-correction bit per window, and adjacent windows agree at their seams without explicit merging. This answers the criticism that AlphaQubit outputs one global bit per experiment.
  - Claims real-time surface-code decoding up to d = 25 on a single TPU v6e, and notes that AQ2 is limited to d = 11 in real time.
  - Sources: [arXiv:2601.09921](https://arxiv.org/html/2601.09921) (also the source of the d = 25 and "AQ2 limited to d = 11" statements, via a search summary); [2509.03815](https://arxiv.org/pdf/2509.03815)
- **Lee, Hur, Park (Yonsei; Oct 2025), "Scalable Neural Decoders for Practical Real-Time QEC", arXiv:2510.22724 → npj Quantum Information (2026).**
  - A Mamba state-space decoder costing O(d²) instead of O(d⁴) for attention.
  - With decoding-latency-induced noise, its threshold is 0.0104 vs 0.0097 for a transformer, and the gap widens up to d = 9.
  - On real d = 3/5 hardware data its LER is comparable to transformer and tensor-network decoders [secondary].
  - Practical deployment faces edge-compute limits.
  - Sources: [arXiv:2510.22724](https://arxiv.org/abs/2510.22724); [npj QI](https://www.nature.com/articles/s41534-026-01403-y); [Lacuna](https://lacuna.tiptreesystems.com/work/scalable-neural-decoders-for-practical-real-time-quantum-error-correction/wrk_ed1ff4e832abea621ba4939aa2f8c007) [secondary]; [Quantum Computing Report](https://quantumcomputingreport.com/yonsei-university-researchers-develop-mamba-state-space-neural-decoder-to-enable-real-time-quantum-error-correction/) [secondary]
- **"Real-time Surface-Code Error Correction Using an FPGA-based Neural-Network Decoder" (May 2026), arXiv:2605.04892 (authors n/r).** A **closed-loop hardware demonstration** with a d = 3 surface code on a superconducting processor. Closed-loop latency is 550 ns, of which 124 ns is NN decoding and 20 ns syndrome calculation, inside a 1.25 µs QEC cycle. The extrapolated resource limit on a large FPGA is about d = 15. — [arXiv:2605.04892](https://arxiv.org/abs/2605.04892)
- **Earlier FPGA and in-memory designs.**
  - 2023 FPGA neural decoder: 197 ns at d = 5, 1.136 µs at d = 7, and 540 ns full-system latency (arXiv:2305.15767). — [2305.15767](https://arxiv.org/html/2305.15767)
  - 2024 compute-in-memory feedforward design: under 440 ns for d = 3–9, with 87.6 ns NN time on FPGA (arXiv:2411.18090). — [2411.18090](https://arxiv.org/pdf/2411.18090)
  - A 2026 review says INT4-quantized, pruned neural decoders can meet a 1 µs deadline up to d = 9 on current FPGAs (arXiv:2605.17230). — [2605.17230](https://arxiv.org/pdf/2605.17230)
- **Neural pre-decoders.**
  - QAdapt (arXiv:2607.28422) does noise-adaptive neural pre-decoding and hands the residual syndrome to a global decoder. — [2607.28422](https://arxiv.org/pdf/2607.28422)
  - QuantiSpect, a 3D-CNN pre-decoder (arXiv:2607.18204). [title-level] — [2607.18204](https://arxiv.org/pdf/2607.18204)
  - Fully convolutional 3D decoders for circuit noise (arXiv:2506.16113). [title-level] — [2506.16113](https://arxiv.org/pdf/2506.16113)
  - "Arqade" pairs a quantized 3D CNN on FPGA with Sparse Blossom on CPU and extends pre-decoding to qLDPC (described in the hierarchical-decoding review). — [arXiv:2611.16271](https://arxiv.org/html/2611.16271v1)

**6c. GNN and code-agnostic decoders**
- **Lange et al., "Data-driven decoding of quantum error correcting codes using graph neural networks", Phys. Rev. Research 7, 023181 (2025); arXiv:2307.01241.**
  - Frames decoding as graph classification of an annotated detector graph.
  - On circuit-level surface-code noise it outperforms matching using only simulated data, even when matching is given the full error model.
  - On Google's real repetition-code data (Nature 614, 676, 2023) it is on par with MWPM.
  - Sources: [arXiv:2307.01241](https://arxiv.org/abs/2307.01241); [Chalmers PDF](https://research.chalmers.se/publication/546994/file/546994_Fulltext.pdf)
- **Hu et al., "Efficient and Universal Neural-Network Decoder for Stabilizer-Based QEC" (GraphQEC), arXiv:2502.19971 (Feb 2025).** Code-agnostic GNN tested on surface, colour and qLDPC codes; the authors claim "unprecedented accuracy and efficiency". — [arXiv:2502.19971](https://arxiv.org/abs/2502.19971)
- **"A Code-Agnostic Graph Neural Network Decoder from the Detection Error Model", arXiv:2610.01683.** [title-only] — [2610.01683](https://arxiv.org/pdf/2610.01683)

**6d. qLDPC: transformers, neural BP and learned message passing**
- **Blue, Avlani, He, Ziyin, Chuang, "Machine Learning Decoding of Circuit-Level Noise for Bivariate Bicycle Codes", arXiv:2504.13043 (Apr 2025) → Quantum (2026).**
  - A recurrent transformer.
  - On [[72,12,6]] at p = 0.1%, its LER is about 5× lower than BP-OSD and about 5× higher than a most-likely-error decoder. Its runtime is consistent and about 10× faster than worst-case BP-OSD.
  - On [[144,12,12]] its LER is worse, but it keeps the speed advantage. The authors say "more complex architectures and/or training procedures" are needed to scale.
  - Source: [arXiv:2504.13043](https://arxiv.org/abs/2504.13043)
- **GraphMP (Oct 2026), arXiv:2610.41398.** Search-distilled learned message passing, "broadly competitive with Relay-BP" (details in §1). — [arXiv:2610.41398](https://arxiv.org/abs/2610.41398)
- **"Scale-Transferable Neural Belief Propagation for Quantum Protograph LDPC Codes …" (IEEE Access, 2026).** One weight set trained on the smallest code of a family decodes [[180,36]] to [[720,144]] without retraining, staying within a factor 1.3–1.7 of near-ML BP+OSD where plain BP fails. Code-capacity noise only. — [Zenodo record](https://zenodo.org/records/21201538)
- **"Learning to Decode QLDPC Codes via Cluster-Based Sequential BP" (arXiv:2607.20130).** Learned BP scheduling on the [[288,12,18]] BB code under depolarizing (code-capacity) noise. — [2607.20130](https://arxiv.org/html/2607.20130)
- **"Toward Uncertainty-Aware and Generalizable Neural Decoding for QLDPC Codes" (arXiv:2510.06257).** Apparently the source of "QuBA", a Bayesian graph neural decoder; this mapping is not confirmed. — [2510.06257](https://arxiv.org/pdf/2510.06257)
- **Earlier neural BP.** Liu & Poulin, "Neural Belief-Propagation Decoders for Quantum Error-Correcting Codes" (PRL 2019). — [ResearchGate](https://www.researchgate.net/publication/333307424_Neural_Belief-Propagation_Decoders_for_Quantum_Error-Correcting_Codes)
- **Non-ML qLDPC baselines (2024–2026) [title-level]:**
  - Relay-BP (referenced by GraphMP and Rapids)
  - BP-OSD
  - Beam search decoder (arXiv:2512.07057) — [2512.07057](https://arxiv.org/abs/2512.07057)
  - Automorphism ensemble decoding (arXiv:2503.01738) — [2503.01738](https://arxiv.org/pdf/2503.01738)
  - Almost-linear-time decoder under circuit noise (arXiv:2409.01440) — [2409.01440](https://arxiv.org/pdf/2409.01440)
  - PriorTilt, sparse perturbations of BP priors (arXiv:2611.13373) — [2611.13373](https://arxiv.org/pdf/2611.13373)
  - GPU parallel BP (arXiv:2610.13879) — [2610.13879](https://arxiv.org/pdf/2610.13879)
  - Adaptive decoding through decoder disagreement (arXiv:2609.37629) — [2609.37629](https://arxiv.org/abs/2609.37629)
  - Tesseract (A* search with an admissible heuristic) and its acceleration (arXiv:2602.02985) and BP-informed version (arXiv:2611.03128) — [2602.02985](https://arxiv.org/pdf/2602.02985); [2611.03128](https://arxiv.org/pdf/2611.03128)

**6e. RL decoders (earlier work, pre-2023)**
- **Sweke, Kesselring, van Nieuwenburg, Eisert, "Reinforcement Learning Decoders for Fault-Tolerant Quantum Computation" (arXiv:1810.07207; MLST 2020).** Recasts fault-tolerant decoding as repeated agent–environment interaction; deep Q-learning gives fast surface-code decoding agents under several noise models. — [arXiv:1810.07207](https://arxiv.org/pdf/1810.07207.pdf)
- **Andreasson et al., "Quantum error correction for the toric code using deep reinforcement learning" (arXiv:1811.12338; Quantum 2019).** The state is the syndrome and an action moves a defect. In the authors' words it "does not outperform state of the art decoders"; its value is as a proof of concept. — [arXiv:1811.12338](https://arxiv.org/pdf/1811.12338); [survey 2608.15760](https://arxiv.org/pdf/2608.15760)
- **Domingo Colomer, Skotiniotis, Muñoz-Tapia (arXiv:1911.02308).** Toric-code RL that rewards the agent only when logical information is preserved. — [arXiv:1911.02308](https://arxiv.org/pdf/1911.02308)
- **BNL, "Decoding surface codes with deep RL and probabilistic policy reuse" (arXiv:2212.11890).** [title-only] — [2212.11890](https://arxiv.org/pdf/2212.11890)

**6f. Decoders trained or evaluated on hardware data (beyond Google)**
- **NVIDIA–QuEra transformer** on QuEra's d = 3 distillation circuit (§5). QuEra's Harvard collaborators also train a CNN decoder offline for qLDPC codes [secondary]. — [NVIDIA blog](https://developer.nvidia.com/blog/nvidia-and-quera-decode-quantum-errors-with-ai); [QuEra FAQ page](https://www.quera.com/questions/how-is-ai-used-to-accelerate-quantum-computing) [secondary]
- **Quantinuum blog "Real Time Error Correction at Increased Scale".** Describes real-time error correction with AI decoders on Helios, with about 12× error reduction per code-distance step [snippet-level; details n/r]. — [Quantinuum blog](https://www.quantinuum.com/blog/real-time-error-correction-at-increased-scale)
- **IBM devices.**
  - ANN syndrome decoder for an adjusted heavy-hex code on several IBM devices, compared with MWM (arXiv:2311.15146, 2023). — [2311.15146](https://arxiv.org/pdf/2311.15146)
  - "Calibration-Conditioned FiLM Decoders for Low-Latency Decoding … Evaluated on IBM Repetition-Code Experiments" (arXiv:2601.16123). [title-only] — [2601.16123](https://arxiv.org/pdf/2601.16123)
  - Error-corrected memory and logic on heavy-hex (arXiv:2610.11658); its decoder type n/r. — [2610.11658](https://arxiv.org/pdf/2610.11658)
- **Varbanov et al., neural-network decoder for near-term surface-code experiments (Phys. Rev. Research 2025).** [title-level, via citation] — [citation context](https://arxiv.org/pdf/2610.11658)
- **"Neural Decoders for Universal Quantum Algorithms" (arXiv:2509.11370).** Decoding logical circuits rather than memory experiments. [title-only] — [2509.11370](https://arxiv.org/pdf/2509.11370)

**6g. Learned decoder priors and detector error models**
- **Sivak, Newman, Klimov, "Optimization of decoder priors for accurate quantum error correction", Phys. Rev. Lett. 133, 150603 (2024); arXiv:2406.02700.**
  - An RL-inspired calibration of DEM priors that directly minimizes LER (of an MWPM decoder), applied to Sycamore repetition- and surface-code experiments.
  - A truncated snippet suggests about 10.6% LER reduction at d = 5 relative to an uninformed prior [verify].
  - Later uses: the dynamic-surface-codes experiment used it as one prior configuration, and the RL-control paper retrained a DEM with it.
  - Sources: [arXiv:2406.02700](https://arxiv.org/pdf/2406.02700); [Sim-to-real paper](https://arxiv.org/html/2609.04557); [Dynamic surface codes](https://arxiv.org/html/2412.14360)
- **"Logical error estimation from syndrome data of surface-code experiments" (arXiv:2606.11496).** On Willow data, DEMs estimated from syndromes perform comparably to the RL-optimized models. — [2606.11496](https://arxiv.org/html/2606.11496v2)

**6h. Accuracy and latency at a glance (as reported; not like-for-like)**

| Decoder | Code / setting | Accuracy claim | Speed / latency claim | Source |
|---|---|---|---|---|
| AQ1 (2024) | surface d = 3, 5 (Sycamore data) | −6% vs tensor network; −30% vs correlated matching | too slow for real time | [blog](https://blog.google/technology/google-deepmind/alphaqubit-quantum-error-correction/) |
| Willow real-time decoder (2025) | surface d = 5 | Λ ≈ 2.0 (real time) vs 2.14 (offline) | ~63 µs average latency; 1.1 µs cycle | [Nature](https://www.nature.com/articles/s41586-024-08449-y.pdf) |
| AQ2-RT (2025) | surface ≤ d = 11; colour ≤ d = 9 | near-optimal; better than the Willow real-time decoder | < 1 µs per cycle throughput | [arXiv:2512.07737](https://arxiv.org/abs/2512.07737) |
| Self-coordinating parallel NN (2026) | surface ≤ d = 25 | "state-of-the-art accuracy" | real time on 1 TPU v6e | [arXiv:2601.09921](https://arxiv.org/html/2601.09921) |
| Mamba (2025/26) | surface ≤ d = 9 | threshold 0.0104 vs 0.0097 (transformer) | O(d²) vs O(d⁴) | [arXiv:2510.22724](https://arxiv.org/abs/2510.22724) |
| FPGA NN (2026) | surface d = 3, on hardware | n/r | 550 ns closed loop (124 ns NN) | [arXiv:2605.04892](https://arxiv.org/abs/2605.04892) |
| Rapids (2026) | surface, colour, [[144,12,12]] | 240× vs real-time hardware decoders; 16× vs real-time Relay | FPGA, integer arithmetic | [arXiv:2611.26447](https://arxiv.org/abs/2611.26447) |
| Transformer for BB (2025) | [[72,12,6]], p = 0.1% | 5× better than BP-OSD; 5× worse than MLE | ~10× faster than worst-case BP-OSD | [arXiv:2504.13043](https://arxiv.org/abs/2504.13043) |
| PRISM (2026) | FBQC surface, colour, qLDPC ≤ d = 24 | matches ILP decoders | "significantly faster" | [arXiv:2611.13630](https://arxiv.org/abs/2611.13630v1) |
| GraphMP (2026) | qLDPC (15 training codes) | competitive with Relay-BP | n/r | [arXiv:2610.41398](https://arxiv.org/abs/2610.41398) |

### Inferences
- **The accuracy–latency frontier moved in 2025–2026.** Transformer decoders went from offline-only (AQ1) to real-time throughput (AQ2-RT, parallel windows, Mamba, FPGA integer networks). The remaining questions are end-to-end *latency*, which matters for feed-forward, as opposed to *throughput*, and power and edge deployment.
- **qLDPC neural decoding is the weakest point.** The only peer-reviewed circuit-level BB result degrades at [[144,12,12]]. The strong 2026 claims (Rapids 16× over Relay; GraphMP parity with Relay-BP) are unreplicated preprints.
- **Hardware rankings can differ from simulation rankings.** The sim-to-real study and the Mamba/AQ2 results on Willow data show this. Learned DEM or prior fitting (Sivak 2024; syndrome-estimated DEMs) is the other way to close the gap.

### Gaps
- Exact AQ2 Willow LER numbers, PRISM authors, the FPGA-decoder authors and Varbanov et al. details were n/r.
- Specific architectures and Λ values behind Quantinuum's Helios AI decoders were not confirmed.
- No source found trained an ML decoder on IBM qLDPC hardware data.

---

## 7. Hardware-in-the-loop optimization and automated calibration of QEC

### Takeaway
There are two flagship hardware-in-the-loop RL demonstrations:
- **Yale (Sivak et al., Nature 2023):** in-situ RL tuning of GKP error correction beyond break-even, G = 2.27 ± 0.07.
- **Google (arXiv:2511.08493 → Nature 2026):** surface-code detection events used as the reward for an RL agent that steers more than 1,000 control parameters during computation, giving 3.5× better logical stability against injected drift.

Around them sit RL-inspired decoder-prior fitting (PRL 2024), simulated multi-agent RL with bandit retraining (Quantum Sci. Technol. 2026), and a 2026 security analysis of agentic calibration advice. No hardware-in-the-loop discovery of circuit *structure* (schedules, gadgets) was found.

### Cited Findings
- **Sivak et al., "Real-time quantum error correction beyond break-even", Nature 616, 50–55 (2023); arXiv:2211.09116.**
  - GKP logical qubit in a microwave oscillator, with gain G = 2.27 ± 0.07 over the best passive encoding and T_X = T_Z = 2.20 ± 0.03 ms.
  - An FPGA orchestrates QEC (a few hundred ns latency), and a GPU-based RL agent optimizes its parameters in situ, accounting for real error channels and control imperfections.
  - Sources: [IDEAS/RePEc](https://ideas.repec.org/a/nat/nature/v616y2023i7955d10.1038_s41586-023-05782-6.html); [Yale PDF](https://rsl.yale.edu/sites/default/files/2023-v_sivak-real-time_quantum_error_correction_beyond_break-even.pdf); [ar5iv](https://ar5iv.labs.arxiv.org/html/2211.09116); [DOE highlight](https://www.energy.gov/node/4830845)
- **Google Quantum AI, "Reinforcement Learning Control of Quantum Error Correction", arXiv:2511.08493 (Nov 2025) → Nature (2026), DOI 10.1038/s41586-026-10759-2.**
  - Detection events serve both for error correction and as the RL training signal ("unify calibration with computation").
  - More than 1,000 analog control parameters, which map the abstract QEC circuit to waveforms.
  - 3.5× better logical stability under injected drift (arXiv abstract). A fixed policy's error rates climbed with drift, while the RL policy kept them low and stable.
  - Secondary reports mention "20%" reductions that measure different things. A newsletter cites 7.72×10⁻⁴ per cycle at d = 7 on Willow from an "agent-based calibration" result; the attribution to this paper is **unverified**.
  - Sources: [arXiv:2511.08493](https://arxiv.org/abs/2511.08493); [Nature](https://www.nature.com/articles/s41586-026-10759-2); [Google Research blog](https://research.google/blog/towards-a-quantum-computer-that-learns-from-its-errors/); [APS GPS 2026 talk](https://meetings-archive.aps.org/smt/2026/mar-z05/1/); [Quantum Computing Report](https://quantumcomputingreport.com/google-research-stabilizes-willow-quantum-processor-using-continuous-reinforcement-learning-control-layers/) [secondary]; [Latent Dynamics](https://latentdynamics.substack.com/p/the-microsecond-problem) [secondary]
- **Decoder-prior optimization on hardware data (Sivak, Newman, Klimov, PRL 2024).** See §6g. — [arXiv:2406.02700](https://arxiv.org/pdf/2406.02700)
- **Guatto et al. (Quantum Sci. Technol. 2026).** Multi-agent RL QEC cycles plus a BRAVE bandit for online retraining under drift; results are from simulation per the retrieved sources. — [IOP](https://iopscience.iop.org/article/10.1088/2058-9565/aebdc4)
- **Özgüler (2026).** Certified acceptance of AI-proposed QEC updates under bounded drift (§1). — [arXiv:2609.19090](https://arxiv.org/abs/2609.19090)
- **QIQCBench (2026).** Includes calibration and control tasks for agents in a virtual lab (§1). — [arXiv:2609.17439](https://arxiv.org/html/2609.17439v1)

### Inferences
- **Both hardware RL demonstrations tune *parameters* of a fixed QEC circuit.** Discovering circuit structure (schedules, flags, gadgets) has so far happened only in simulation (§§2–3). A practical hybrid is to discover structure in simulation and then refine it in situ with RL against detection events.
- **The detection-event reward used by Google is cheap, dense and available on hardware.** It suits online agents, but as Özgüler shows, passive syndrome statistics can be blind to some coherent errors.

### Gaps
- Not covered because the search budget ran out: Marquardt-group or other bosonic RL feedback works after 2023, Google's "Snake" frequency optimizer and calibration graphs, and AWS/Alice&Bob cat-qubit RL.
- The exact Nature 2026 figures for the RL-control paper were n/r.

---

## 8. Tooling: simulators and verifiers as RL/search environments, and benchmarks

### Takeaway
Fast, verifiable rewards exist where the setting is Clifford circuits with Pauli noise:
- exhaustive fault enumeration (Lahiri & Marquardt);
- parallel stabilizer simulation (Zen et al.);
- DEM extraction, sampling and decoding (OmniQEC, AlphaSyndrome), sped up by importance sampling (FastSched);
- MaxSAT/SMT proxies on the DEM (PropHunt, Proxima);
- ILP construction plus independent algebraic verification (factories).

Non-Clifford protocols (cultivation) need specialized simulators. Benchmarks are a weak point: QIQCBench covers general agentic quantum engineering, and the sim-to-real Willow study covers decoder ranking, but no standard benchmark for QEC circuit discovery was found.

### Cited Findings
- **Verifiers and rewards used by the works above.**
  - Exhaustive fault check up to order t — [Lahiri & Marquardt](https://arxiv.org/html/2611.18504)
  - "Fast, parallelized stabilizer quantum circuit simulator" with a tailored reward — [Zen et al.](https://arxiv.org/abs/2402.17761)
  - Decoded LER in every MCTS rollout — [AlphaSyndrome](https://arxiv.org/html/2601.12509)
  - Importance-sampled LER — [FastSched](https://arxiv.org/abs/2609.12020)
  - MaxSAT ambiguity — [PropHunt](https://arxiv.org/html/2601.17580v1)
  - SMT on DEM column subsets — [Proxima](https://arxiv.org/abs/2611.19715)
  - DEM, then sampling, then BP+OSD at p = 0.005 — [OmniQEC](https://arxiv.org/html/2607.25865v1)
  - ILP plus independent verification engine — [factories](https://arxiv.org/abs/2610.06535)
  - Stim compilation plus distance-estimate filtering — [High-rate qLDPC processors](https://arxiv.org/pdf/2607.28795)
- **Open-source code released with these works.**
  - MQT QECC (SAT state preparation) — [Peham et al.](https://arxiv.org/abs/2408.11894v3)
  - rlftqc (RL FT state preparation) — [GitHub](https://www.Github.com/remmyzen/rlftqc)
  - Left-right-circuit code on GitHub, circuits on Zenodo — [Strikis et al.](https://arxiv.org/abs/2603.05481)
  - Magic-state factory repository — [GitHub](https://github.com/shubhamj810/magic-state-factory-repo)
  - Adaptive_QEC (multi-agent RL plus bandit) — [GitHub](https://github.com/ManuelGuatto/Adaptive_QEC)
  - IBM OpenEvolve (LLM-guided evolution) — [IBM blog](https://research.ibm.com/blog/ai-for-qec)
  - QUITS, a modular qLDPC circuit simulator (arXiv:2504.02673) — [QUITS](https://arxiv.org/pdf/2504.02673)
- **Non-Clifford simulation for cultivation rewards.** Few-Clifford-term ZX simulation (about 8 terms on average at d = 5, p = 0.1%), plus the related tools SOFT, Clifft and tsim. — [arXiv:2509.08658](https://arxiv.org/abs/2509.08658)
- **Decoders used as reward engines or baselines.** BP+OSD (OmniQEC; Blue et al.), Relay-BP (GraphMP; Rapids), Sparse Blossom (Arqade), Chromobius and "Libra" (AQ2 comparisons), Tesseract (A*), and most-likely-error decoders (Blue et al.; QuEra). — [OmniQEC](https://arxiv.org/html/2607.25865v1); [GraphMP](https://arxiv.org/abs/2610.41398); [Rapids](https://arxiv.org/abs/2611.26447); [Hierarchical decoding](https://arxiv.org/html/2611.16271v1); [AQ2](https://arxiv.org/html/2512.07737v1); [Tesseract](https://arxiv.org/pdf/2602.02985); [Blue et al.](https://arxiv.org/abs/2504.13043)
- **Harness design for expensive verification.** QSA's coreset with approximate scoring, static/snapshot analyses and task-specific rewards target "expensive verification, sparse feedback, heterogeneous objectives". — [arXiv:2609.38883](https://arxiv.org/html/2609.38883)
- **Benchmarks.** QIQCBench (49 tasks, 17 agentic systems) — [arXiv:2609.17439](https://arxiv.org/html/2609.17439v1); sim-to-real decoder benchmarking on Willow — [arXiv:2609.04557](https://arxiv.org/abs/2609.04557)
- **Reward hacking.** Optimizing distance alone produced highly connected, impractical codes. — [qBraid](https://www.qbraid.com/blog-posts/qbraid-uses-alphaevolve-for-quantum-error-correction)

### Inferences
- **Recommended reward stack**, from cheapest to most faithful:
  1. exact combinatorial certificates (fault enumeration up to t, distance checks, ILP/SAT feasibility);
  2. DEM-level proxies (ambiguity, proximity, residual distance);
  3. importance-sampled LER under a fixed decoder;
  4. final LER with confidence intervals under several decoders and noise models.
  Every 2026 paper uses a different subset of these, which explains why their claims disagree.

### Gaps
- Not fetched or searched (budget): documentation and performance claims for Stim, sinter, PyMatching (Sparse Blossom) and ldpc/BP-OSD, and any QEC-specific RL "gym" environments.
- No public leaderboard for syndrome-extraction scheduling or FT-gadget discovery was found.

---

## 9. Gaps and opportunities

### Takeaway
The literature now has (a) RL, MCTS and LLM search for Clifford FT circuits at small to moderate distance, (b) fast neural decoders for topological codes, and (c) two hardware-in-the-loop RL demonstrations. Large open areas remain: LLM- or ML-discovered *decoders*, ML for logical operations (lattice surgery, code switching, fold/automorphism gates), cultivation design, Floquet scheduling, shared benchmarks, hardware validation of discovered circuits, and verifiers for noise beyond Pauli.

### Cited Findings (evidence behind each gap)
- **No LLM-discovered decoders.** Targeted FunSearch, AlphaEvolve and OpenEvolve queries found none. The nearest is GraphMP's search-to-model distillation loop. — [arXiv:2606.25870 context](https://arxiv.org/html/2606.25870); [GraphMP](https://arxiv.org/abs/2610.41398)
- **No ML/LLM work on code switching or lattice-surgery synthesis** beyond Lahiri & Marquardt's logical gates. Lattice surgery relies on SAT and ZX methods. — [arXiv:2611.18504](https://arxiv.org/pdf/2611.18504); [Lattice surgery beyond surface code](https://arxiv.org/abs/2504.10591)
- **Cultivation's degraded fault distance was found by a human**, and non-Clifford verifiers only became practical in 2025–2026. — [arXiv:2609.17706](https://arxiv.org/pdf/2609.17706); [arXiv:2509.08658](https://arxiv.org/abs/2509.08658)
- **Mutually inconsistent scheduler claims** from AlphaSyndrome, PropHunt, FastSched, Proxima and Structure-guided on different benchmarks. — [FastSched](https://arxiv.org/abs/2609.12020); [Proxima](https://arxiv.org/abs/2611.19715); [Structure-guided](https://arxiv.org/abs/2611.25825)
- **Simulator rankings may not carry over to hardware.** — [Sim-to-real](https://arxiv.org/abs/2609.04557)
- **qLDPC neural decoders degrade at [[144,12,12]].** — [Blue et al.](https://arxiv.org/abs/2504.13043)
- **RL scalability is a "major challenge"** for circuit discovery. — [review arXiv:2607.00365](https://arxiv.org/pdf/2607.00365)
- **Agent reliability and security problems.** QIQCBench shows large variance across agents; Özgüler shows 12/12 harmful choices under attack; qBraid documents reward hacking. — [QIQCBench](https://arxiv.org/html/2609.17439v1); [Özgüler](https://arxiv.org/abs/2609.19090); [qBraid](https://www.qbraid.com/blog-posts/qbraid-uses-alphaevolve-for-quantum-error-correction)
- **Coherent errors can be invisible to passive syndrome statistics.** — [Özgüler](https://arxiv.org/abs/2609.19090)

### Inferences (concrete open problems)
1. **LLM-evolved decoder programs.** Use FunSearch, AlphaEvolve or OpenEvolve to evolve BP schedules, OSD ordering rules, matching-weight corrections, pre-decoder lookup policies or window-merging logic. The reward would be Stim-sampled LER at fixed latency and compute budgets. Held-out codes and noise models are needed to catch overfitting and reward hacking.
2. **A shared benchmark for syndrome-extraction scheduling.** It would fix:
   - the code suite (surface, colour, BB [[72,12,6]]/[[144,12,12]], lifted-product, random quantum Tanner);
   - the noise models (uniform depolarizing, SI1000-like, biased, with idling);
   - the decoders (BP-OSD, Relay-BP, MWPM/correlated matching, a neural decoder);
   - the LER-estimation protocol (confidence intervals, importance sampling) and compute-time accounting.
   AlphaSyndrome, PropHunt, FastSched, Proxima, LRC and structure-guided synthesis would then be re-run on it.
3. **Co-optimizing schedules and decoders, and testing transfer.** Check whether schedules found with one decoder in the loop stay optimal under real-time or neural decoders, and co-train schedule and decoder together.
4. **Hardware validation of discovered circuits.** Run RL-, MCTS- or LLM-found state-preparation and scheduling circuits on hardware, and refine them in situ with detection-event RL in the style of Google's 2025–26 work.
5. **LLM or RL search for logical operations beyond state preparation.** Targets: lattice-surgery subroutines (with SAT-scalpel outputs as baselines and verifiers), code-switching protocols, and fold-transversal or automorphism logical gates in qLDPC codes. Exact Clifford verifiers already exist for all of these.
6. **Search over cultivation and distillation hybrids.** Use few-Clifford-term non-Clifford simulators and exact LER tools as rewards, with fault-distance checks as hard constraints. Learned post-selection would aim to raise the 8% retention.
7. **Circuit-level and hardware-aware evaluation of the 564 AI-found factory classes**, including connectivity, qLDPC-native layouts and CCZ or other non-Clifford targets.
8. **Measurement-schedule discovery for Floquet and dynamic codes** with automated search. No ML work was found.
9. **Scaling FT-circuit discovery beyond d ≈ 9** through generators of parametric circuit families (programs, not instances), with verification that grows sub-combinatorially (symbolic or ZX fault-equivalence instead of exhaustive enumeration).
10. **Neural decoders for qLDPC and for logical circuits at real-time latency**, with independent replication of the 2026 claims (Rapids, GraphMP), FPGA/ASIC power budgets, and decoding of logical-gate circuits rather than memory experiments only.
11. **Learned noise models and drift tracking.** Combine syndrome-estimated DEMs, RL-optimized priors and neural decoders in one online-adaptive loop, extended to qLDPC hardware and to leakage and coherent noise.
12. **Verified agentic QEC pipelines.** Certified acceptance tests, adversarial evaluation (prompt injection, misleading advice), multi-objective rewards against reward hacking, and QEC-specific agent benchmarks beyond QIQCBench's 49 general tasks.
13. **Verifiers for noise beyond Pauli.** Cheap rewards for coherent, leakage and non-Markovian noise inside discovery loops. Nearly every discovery loop surveyed here assumes Clifford circuits with Pauli noise.

### Gaps (limits of these notes)
- arXiv full texts could not be fetched; the search budget ran out before fold-transversal/automorphism search, Floquet scheduling, RL environments, tooling documentation, bosonic RL after 2023, and Quantinuum/IBM hardware-decoder details were covered.
- Many 2026 numbers are the authors' own preprint claims, and several rely on secondary summaries, as flagged inline.
