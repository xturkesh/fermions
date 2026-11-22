# LLM-based agentic systems for discovery in quantum information and quantum physics (focus: 2025 to Nov 2026)

Method note, read first. Following the coordinator's scope update, the newest LLM, coding-agent and AlphaEvolve-style discovery of **quantum circuits and quantum error-correcting codes** comes first and is listed newest first. General physics agents and benchmarks come after. The "Gaps and opportunities" section is general and not tailored to any group.

Several limits apply to how this was sourced:
- WebFetch could not reach arxiv.org, deepmind.google or scottaaronson.blog from this session, because those hosts are not in the session's allowed network.
- The shared web-search budget (200 calls per turn across all agents) ran out partway through.
- All details therefore come from search-engine extracts of arXiv abstract and HTML pages, project repositories, official blogs and press coverage. **No PDF was read in full.** Check figures against the PDFs before quoting them in a final report.
- "Title-only" means the item was seen only as a title or a reference-list entry.
- Most 2026 items are unrefereed preprints.

---

## 1. What have LLM or coding agents discovered in quantum error correction (codes)? (2025 to Nov 2026, newest first)

### Takeaway
Between June and November 2026 at least seven preprints used LLM-guided program evolution to search for QEC codes. The engines were OpenEvolve, AlphaEvolve, or custom multi-agent and Claude-Code-style loops. The targets were mainly qLDPC families (bivariate-bicycle (BB), perturbed-BB (PBB), lifted-product, balanced-product) plus entanglement-assisted codes.

All of these works share three features:
- They evolve **programs that generate codes**, not the codes themselves.
- They score candidates with **exact, cheap algebraic evaluators**: GF(2) rank, and distance certified by MILP or exhaustive search.
- They increasingly add **circuit-level logical-error-rate (LER) simulation** with BP-OSD-family decoders.

Reported outcomes:
- new finite-length codes with leading or competitive parameters;
- codes that beat the [[72,12,6]] and [[144,12,12]] BB benchmarks at equal total-qubit budgets;
- one provable infinite EAQECC family checked in Lean.

So far the results have been validated by the authors themselves and by a public, automatically verified leaderboard, not by peer review.

### Cited Findings
- **arXiv:2610.41183 (30 Oct 2026). Jae-Hyun Baek, Byung-Sun Won, Jon-Lark Kim, "Researcher-in-the-Loop Evolutionary Discovery of Entanglement-Assisted Quantum Code Families" (cs.AI; cs.IT).**
  - Architecture: AlphaEvolve is the evolution engine. A researcher supplies target parameters, an exact evaluator and construction ideas. This is followed by human mathematical generalization and formal proof checking.
  - Result: the evolved programs gave finite qubit EAQECC constructions. The researchers generalized these into the family [[n,1,n−1; n−q−1]]_q for prime powers q ≥ 3 and n ≥ 2q+1, plus a q = 2 case. Lean checks the algebraic certificate. Combined with entanglement lifting, the family "closes 54 gaps in archived EAQECC tables."
  - Effect of shared guidance on held-out success: from 50.0% to 68.8% at length 9, and from 15.6% to 50.0% at length 11.
  - Stated stance: LLM-evolved constructions need "distinct evidence" for validity, generality and reproducibility.
  - Practical limit: rerunning the AlphaEvolve parts needs a Google Cloud project allowlisted for AlphaEvolve.
  - Sources: [arXiv](https://arxiv.org/abs/2610.41183); [code, data and Lean proofs](https://github.com/LeGenAI/llm_guided_EAQECC_search)
- **arXiv:2610.06623 (Oct 2026). Juan Cruz-Benito (IBM), "Large Language Model-Guided Discovery of Weight-Five Bivariate Bicycle Codes" (quant-ph).**
  - Architecture: OpenEvolve plus an LLM ensemble, evolving Python functions that emit the polynomial pairs or tuples defining weight-5 BB and PBB codes. Random-search controls were included.
  - Results: 1,142 distinct code proposals, of which 1,081 are attributable to LLM-generated programs. Of 936 distinct CSS proposals with a recorded positive distance, 816 (87.18%) are certified at d ≥ 5.
  - **Documented failure mode:** 409 of the 1,142 entries are *disconnected* codes, and 73.1% of the classes with exact distance certificates contain repeated connected components. The LLM campaign has the higher certified yield, but "the random controls cover more connected classes."
  - Community re-verification: pull requests to the Unitary Foundation leaderboard rebuild Table-3 codes from the published polynomials and find logical operators at the reported distances. Examples are [[420,8,14]], [[432,8,15]], [[420,6,16]] and [[420,8,13]].
  - Sources: [arXiv](https://arxiv.org/abs/2610.06623); [PR [[420,8,14]]](https://github.com/unitaryfoundation/qldpc-challenge/pull/2879); [PR [[432,8,15]]](https://github.com/unitaryfoundation/qldpc-challenge/pull/2877); [PR [[420,6,16]]](https://github.com/unitaryfoundation/qldpc-challenge/pull/2875); [PR [[420,8,13]]](https://github.com/unitaryfoundation/qldpc-challenge/pull/2878); [repo PR adding the weight-5 campaigns](https://github.com/qiskit-community/qcode-discovery/pull/9)
- **arXiv:2609.37056 (Sep 2026). Amal Seddas, Vladyslav Shashkov, Maryna Viazovska, Emmanuel Abbe (EPFL), "Evolving Towards Better Codes: LLM-Guided Search for High-Distance Binary Linear Codes".**
  - This is classical coding theory, but relevant as input to CSS and hypergraph-product constructions. It applies LLM-guided search to improving best-known bounds for binary linear codes.
  - The authors conclude that the method "can help find improved codes and complement existing methods." Specific parameters were not visible in the extracts.
  - Source: [arXiv PDF](https://arxiv.org/pdf/2609.37056)
- **arXiv:2608.08996 (10 Aug 2026). Dongheng Qian, Tianyi Li (Fudan), "Multi-agent discovery of practical quantum LDPC codes" (quant-ph; cs.AI).**
  - Architecture: a multi-agent LLM framework that instantiates coset-orbit balanced-product codes. The search space contains bicycle and lifted-product constructions and also non-normal subgroup actions. The paper describes its approach as generalizing Liu and Marquardt's three-level Structured Concept Evolution.
  - Scope: binary CSS codes with n ≤ 400 and overall weight w ≤ 10.
  - Results: "leading or competitive rate-distance performance in every weight class." Examples are [[288,16,18]] at w = 7, [[288,18,18]] at w = 9 and [[234,28,18]] at w = 10. The codes also have low logical failure under code-capacity depolarizing noise with BP-OSD.
  - Sources: [arXiv](https://arxiv.org/abs/2608.08996); [HTML](https://arxiv.org/html/2608.08996)
- **arXiv:2607.25865 (28 Jul 2026). Ge Yan, Shanchuan Li, Pengyue Ma, Qixin Zhang, Pingchuan Ma, Jianping Wang, Min-Hsiu Hsieh, Yuxuan Du (NTU Singapore, CityU Hong Kong, Hon Hai Research Institute), "OmniQEC: discovering practical quantum error-correcting codes by an AI scientist".**
  - Architecture: an orchestrator plus a code-search agent that does LLM-guided program search over executable QEC construction programs. Characterization, syndrome-extraction circuit synthesis and decoder-based evaluation are fixed tools; only the orchestrator and the code-search agent use LLMs. A "dual-loop" strategy avoids running a full circuit simulation for every candidate.
  - Evaluator: the circuit-level LER of the *complete implementation* under a fixed total physical-qubit budget N. Decoders are BP+OSD, BP+LSD and Relay-BP.
  - Experiments: 4 qLDPC construction families × 3 LLM backends (Claude, GPT, DeepSeek) × 14 qubit budgets per backend.
  - Result: the discovered codes outperform the BB codes [[72,12,6]] and [[144,12,12]] at complete-implementation budgets of 98 and 240 physical qubits respectively.
  - Backend comparison (from an alphaXiv overview, not the paper): Claude and GPT exploited successful structures better, while DeepSeek explored more diverse families but reached a slightly weaker frontier.
  - Stated limitation: the system should eventually become end-to-end multi-agent.
  - Sources: [arXiv](https://arxiv.org/abs/2607.25865); [HTML](https://arxiv.org/html/2607.25865v1); [alphaXiv overview](https://www.alphaxiv.org/overview/2607.25865)
- **arXiv:2606.24808 (Jun 2026). Zidu Liu, Florian Marquardt (MPI for the Science of Light), "Large-Language-Model Discovery of Quantum LDPC Codes through Structured Concept Evolution" (SCE).**
  - Architecture: the LLM evolves *construction programs and their internal search strategies*. It does this by changing the finite group, the protograph shape or the entry-assignment rule. Each program pairs a concept specification with an executable Python part, which may run its own inner search to produce a batch of candidate codes.
  - Scope: fewer than 600 data qubits and check weight ≤ 12. The objective is either kd²/n or the circuit-level LER.
  - Results: improves on published lifted-product benchmarks in selected parameter regimes, with exactly certified distances for the leading results. With the same framework and step budget, "the best LLM outperforms all random-mutation and hill-climbing runs on both objectives."
  - The search extracts mention Claude Code as the agent harness; the exact model versions were not verified.
  - Sources: [arXiv HTML](https://arxiv.org/html/2606.24808v1); [code](https://github.com/ziduliu/sce-qldpc)
- **arXiv:2606.02418 (Jun 2026). Juan Cruz-Benito, Andrew W. Cross, David Kremer, Ismael Faro (IBM Research), "Evolutionary Discovery of Bivariate Bicycle Codes with LLM-Guided Search" (quant-ph; cs.AI).**
  - Architecture: FunSearch- and AlphaEvolve-inspired, built on OpenEvolve with MAP-Elites. LLMs mutate Python programs that generate BB and PBB code ansätze.
  - Evaluator: a staged validation pipeline. It combines GF(2) rank, distance estimation and certification, mixed-integer linear programming, and equivalence and deduplication checks.
  - Scale and cost: about 1,650 evolutionary iterations over five campaigns, about 2×10^5 candidate codes screened, about 140 h of computation and about **US$400 of LLM inference**.
  - Results: 465 distinct candidates at n ≤ 360 (97 CSS BB codes and 368 non-CSS PBB codes). The search recovers known high-performing codes. New finite-length representatives include an indecomposable [[288,16,12]] code and higher-weight codes with up to k = 50 at d = 8. The non-CSS perturbed codes match the gross-code figure of merit at [[144,12,12]].
  - A repository mirror quotes k = 54 against a prior best of 16. This was not reconciled with the abstract.
  - Sources: [arXiv](https://arxiv.org/abs/2606.02418); [code](https://github.com/qiskit-community/qcode-discovery); [announcement](https://x.com/QuantumPapers/status/2062940065462907229); [press](https://bookstopology.substack.com/p/ibm-turns-to-llms-to-accelerate-quantum)
- **Unitary Foundation "QEC Challenge" (qldpc-challenge, 2026).** A public, automatically verified leaderboard for qLDPC codes that explicitly welcomes AI help.
  - AI rules: models (Claude, GPT, Gemini and others) are welcome, but the contributor is the author of record and must list the model in the construction field. CONTRIBUTING.md has a paste-ready prompt for agents.
  - Submission: CSS H_X/H_Z matrices as .npz files. A pinned `validate_candidate` must return `passed: true`, the same verdict as CI.
  - Other rules: one code per PR plus a research note. Direct sums of smaller codes are rejected.
  - Scale at the time of the snapshot: 2,610 codes and 41 contributors.
  - Example community submission: [[672,28,24]], a weight-6 two-block code over C28×Dic3 with exact distance.
  - Sources: [announcement](https://unitary.foundation/posts/2026_qec_challenge/); [repo](https://github.com/unitaryfoundation/qldpc-challenge); [scoring and AI policy](https://github.com/unitaryfoundation/qldpc-challenge/issues/1155); [leaderboard](https://unitaryfoundation.github.io/qldpc-challenge); [PR #5348](https://github.com/unitaryfoundation/qldpc-challenge/pull/5348)
- **Earlier, classical (2025). arXiv:2504.00613, Weindel and Heckel, "LLM-Guided Search for Deletion-Correcting Codes".** FunSearch adapted to evolve functions that construct deletion-correcting codes. It improved lower bounds for lengths n = 12, 13 and 16 for two deletions. — [arXiv](https://arxiv.org/pdf/2504.00613)
- **Non-LLM comparators (context only; out of scope):**
  - Bayesian optimization over BB codes (Chengyu, Meister, Carty, Lin, Bondesan; arXiv:2601.18562) reports better convergence and higher final objectives than evolutionary and random baselines. — [arXiv](https://arxiv.org/pdf/2601.18562)
  - GFlowQEC, generative flow networks for sampling qLDPC codes (arXiv:2611.06819; title-only). — [arXiv](https://arxiv.org/pdf/2611.06819)
  - A transformer distance-predictor plus genetic algorithm finding "champion" classical linear codes (arXiv:2512.13370). — [arXiv](https://arxiv.org/html/2512.13370v1)

### Inferences
- **The recipe is FunSearch's.** All of these works evolve a short *generator program* (polynomials, group and protograph choices, assignment rules) rather than a parity-check matrix. This keeps outputs human-readable and lets one program emit a family. That is what made both the EAQECC generalization and the SCE "concept" level possible.
- **The objective is moving from abstract parameters (k, d, kd²/n) to circuit-level LER at a fixed total qubit budget** that includes ancillas (OmniQEC, SCE). This is closer to hardware relevance but far more expensive. Staged or dual-loop filtering is how these works control cost.
- **Degenerate outputs are the code-search form of reward hacking.** These are disconnected codes, direct sums and repeated components. The countermeasures so far are the leaderboard's direct-sum rule and IBM's equivalence and deduplication stage.
- **Evidence that LLMs beat cheap baselines is narrow and mixed.**
  - SCE: the LLM beats random mutation and hill-climbing at equal steps.
  - Cruz-Benito: higher certified yield than random controls, but less coverage of connected classes.
  - Bondesan's group: BO beats evolutionary baselines, not necessarily LLM-based ones.
  - No matched-compute comparison of LLM evolution against BO, RL and GFlowNets was found.
- **Costs are within reach of a small group:** IBM reports about US$400 of LLM spend and about 140 h of compute for roughly 2×10^5 screened candidates.
- **Humans remain essential for generality and proof.** The EAQECC family required researcher generalization and a Lean certificate. None of the qLDPC papers claim asymptotic or family-level theorems.

### Gaps
- PDFs could not be opened, so parameters and percentages come from abstracts, HTML snippets and repositories. No peer-reviewed versions were identified.
- It was not independently checked whether codes such as [[288,16,12]], [[288,16,18]] or [[234,28,18]] are new relative to existing catalogs.
- The RL paper "Designing quantum error correcting codes to fit decoders via reinforcement learning" (arXiv:2608.15754, attributed to Heinis et al.) and a "Li multi-agent" QEC work were seen only in reference lists. Their methods are unverified.
- The INSPIRE entry "Quantum bivariate bicycle codes with weight-8 checks surpassing the BB benchmark" ([INSPIRE](https://inspirehep.net/literature/3200699)) appeared in results, but whether it is LLM-based is unknown.
- Within this session's budget, no LLM-guided search was found for: Floquet or dynamically generated codes, subsystem codes, bosonic (GKP or cat) codes, qudit codes, approximate or covariant codes, codes tailored to biased or erasure noise, or magic-state-distillation codes.

---

## 2. What have LLM or coding agents discovered in quantum circuits, compilations and quantum programs? (newest first)

### Takeaway
The pattern is clear. Since March 2026, and mostly in October and November 2026, agentic search over circuit-generating programs has produced **verified resource reductions in fault-tolerant subroutines**:
- arithmetic and decoders for Shor- and DQI-type algorithms, with Toffoli×qubit reductions from several-fold up to more than 9,100× against naive reversible baselines;
- fermion-to-qubit encodings and layouts, with up to 73% cost reduction;
- T-count records for GF(2^m) multiplication;
- fault-tolerant Clifford gadgets, verified by exhaustive fault enumeration;
- Lean-verified quantum programs, with Google reporting 290× on arithmetic.

The common enabler is a **classically checkable verifier**: reversible-logic simulation, stabilizer simulation, exhaustive fault checks, or the Lean kernel. Google Quantum AI turned this into a $150k Kaggle competition (October 2026 to January 2027). That competition followed the community "ECDSA.Fail" open-autoresearch challenge (from May 2026). Earlier (October 2025), AlphaEvolve was used to evolve Trotter-circuit generators for Google's Willow experiments.

### Cited Findings
- **arXiv:2611.18504 (Nov 2026). Aritra Lahiri, Florian Marquardt (MPI for the Science of Light and FAU Erlangen-Nürnberg), "Agentic AI for Discovering Fault-Tolerant Quantum Circuits".**
  - Architecture: "agentic evolutionary program synthesis." An LLM writes programs that define Clifford circuits, and a circuit verifier guides the search.
  - Evaluator: "an exhaustive check over faults."
  - Tasks: logical state preparation, logical gate design, and logical routines such as code switching and code growth.
  - Result: fault-tolerant circuits were found for all tasks, "matching or improving upon existing benchmarks in several cases." Secondary coverage says the circuits herald all errors dynamically, unlike standard measurement-based encodings.
  - An FAU "Research Spotlight" from March 2026 on AI-discovered fault-tolerant circuits suggests earlier related work from the group; the relation is unverified.
  - Sources: [arXiv](https://arxiv.org/abs/2611.18504); [CASRAI](https://casrai.org/news/agentic-ai-discovering-fault-tolerant-quantum-circuits-2026); [FAU spotlight](https://www.lightmatter.fau.de/2026/03/fau-lmq-research-spotlight-ai-discovered-fault-tolerant-quantum-circuits/); [Substack commentary](https://bookstopology.substack.com/p/autonomous-quantum-engineering-why)
- **arXiv:2611.09822 (Nov 2026). Google DeepMind and Google Quantum AI, "Discovering formally verified quantum programs with AI".**
  - Architecture: a quantum program library in Lean 4. Researchers write a *Lean quantum-Hoare-logic specification* of correctness and cost. "Agent teams" of LLMs then propose programs and iterate on error feedback from the Lean kernel. Teams share results through wikis and discuss ideas at "organized conferences." The work builds on Mathlib and Physlib's Lean-QuantumInfo.
  - Evaluator: a machine-checked proof against the formal specification and cost model.
  - Results, as improvement over prior baselines: 290× for quantum arithmetic, and 11×, 8× and 2.5× for fusion, photochemistry and drug-metabolism simulations.
  - Inconsistency: the parallel-search gain is reported as up to 1.4× in v1 but up to 2.55× in the abstract.
  - A separate framework called "Lean-QuantumAlg" (nine AI-formalized algorithms, from Grover to attacks on RSA-2048 and P-256) is a *different* paper and should not be conflated with this one.
  - Sources: [arXiv](https://arxiv.org/abs/2611.09822); [Quantum Computing Report](https://quantumcomputingreport.com/google-deepmind-and-google-quantum-ai-use-lean-4-to-discover-formally-verified-quantum-algorithms/); [Quantum Zeitgeist](https://quantumzeitgeist.com/google-quantum-ai-designs-quantum-verified/); [aqora](https://aqora.io/papers/2611.09822); [AlphaSignal](https://alphasignal.ai/news/google-deepmind-s-ai-agents-cut-quantum-algorithm-resources-by-290x)
- **arXiv:2610.31352 (Oct 2026). Marcos Crichigno, Steve Flammia, "Automating quantum circuit optimization with AI agents" (quant-ph).**
  - Architecture: "autoresearch loops." Competing coding agents built on Claude models run in a harness with a file system and a shell and iteratively improve programs. Each submission comes with notes explaining the ideas behind it.
  - Evaluator: a *frozen, exact verifier*.
  - Problems: fermion-to-qubit encoding, Hamiltonian layout and unitary synthesis.
  - Result: cost reduced by up to 73% at the largest sizes tested, relative to published or standard constructions. The agents also outperform a state-of-the-art optimizer for fermion-to-qubit encodings.
  - Claimed scope: any optimization problem whose solutions can be verified exactly and efficiently.
  - Source: [arXiv](https://arxiv.org/abs/2610.31352)
- **arXiv:2610.23832 (Oct 2026). Khattar, Rubin, Gidney, Broughton, Low, Babbush (Google Quantum AI), "Leveraging verifiable structure in quantum circuits to find more efficient compilations with agentic AI".**
  - Idea: when bottleneck subroutines use only reversible classical logic, phasing operations and measurement-based uncomputation, they can be simulated efficiently on classical inputs. Agents can therefore be scored exactly.
  - Applications: Reed–Solomon decoding (the bottleneck of decoded quantum interferometry, DQI), and quantum look-up tables and controlled Givens rotations (bottlenecks of quantum-chemistry simulation).
  - Verifier: *kickmix*, a fast simulator in Qualtran. It processes 64 to 512 test inputs per CPU instruction and handles circuits of more than 100M gates on one core.
  - The paper cites the ECDSA.Fail precedent, where more than 100 participants, mostly without a quantum background, used AI coding agents to cut the Toffoli count "roughly in half."
  - Outlook: co-optimize logical arithmetic, qubit layout and gate schedule.
  - Source: [arXiv](https://arxiv.org/abs/2610.23832)
- **Google Quantum AI × Kaggle, "Verifiable Quantum Circuit Compilation" (announced 21 Oct 2026).**
  - Prize pool $150,000. Tasks relate to fusion energy, optimization and chemistry. It opened 22 Oct 2026 with a final deadline of 20 Jan 2027, per Quantum Computing Report; this was not checked on the Kaggle page itself.
  - It uses kickmix, an open-source C++20/Python verification library in Qualtran. A baseline circuit and a standalone verifier are released.
  - Milestone prizes require public notebooks.
  - Sources: [Google Open Source blog](https://opensource.googleblog.com/2026/10/announcing-verifiable-quantum-circuit-compilation-competition-on-kaggle.html); [Kaggle](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation); [Quantum Computing Report](https://quantumcomputingreport.com/google-quantum-ai-launches-150k-kaggle-competition-to-compress-fault-tolerant-quantum-circuits-using-ai-agents/); [The Quantum Insider](https://thequantuminsider.com/2026/10/23/google-launches-150000-ai-competition-to-reduce-quantum-computing-demands/)
- **arXiv:2610.22367 (15 Oct 2026). Andi Gu, Alek H. Kemeny, Susanne F. Yelin, Mikhail D. Lukin, "Discovering state-of-the-art quantum circuits with large language models".** Press frames this as a Harvard, Anthropic and QuEra collaboration.
  - Architecture: OpenEvolve. A Claude model writes Python programs that generate circuits. Some coverage names the model "Claude Mythos 5"; this is unverified.
  - Evaluator: an independent checker rejects incorrect circuits and scores the rest by **Toffoli count × peak qubit count**. The search keeps space-efficient candidates even when their score on that cost is worse.
  - Results: on 7 of 8 tasks, both Toffoli and qubit counts drop relative to the best known circuits. On the 8th, the search finds a new reversible "divsteps" architecture for prime-field modular inversion; the benchmark there is affine point addition on a 256-bit elliptic curve.
  - Quantum Computing Report's table: the Hermitian DQI decoder goes from (1.0×10^8 Toffoli, 1.9×10^5 qubits) to (9.3×10^5, 2,395), a reduction of more than 9,100×. Other coverage cites 39× on the paper's metric, and an elliptic-curve design with 17× fewer gates for 12% more qubits.
  - Caveats: the paper is unrefereed. The numbers are resource estimates, not hardware runs. Some baselines are naive reversible implementations, because no prior reversible compilation exists for those tasks.
  - Sources: [arXiv HTML](https://arxiv.org/html/2610.22367v1); [Quantum Computing Report](https://quantumcomputingreport.com/harvard-anthropic-and-quera-use-llms-to-discover-quantum-circuits-with-up-to-9100x-lower-resource-costs/); [cellcog](https://cellcog.ai/blog/claude-quantum-circuits/); [MACCUS](https://maccus.fr/en/claude-quantum-circuits/)
- **arXiv:2609.09582 (Sep 2026). "ECDSA.Fail: Open Autoresearch for Optimizing Elliptic-Curve Point Addition in Shor's Algorithm".** Led by Eigen Labs; the author list was not verified.
  - Format: a public challenge launched in late May 2026. "Open autoresearch" means a verifier-gated loop in which humans and AI agents propose and test improvements against a shared objective: peak logical qubits × average executed Toffoli count.
  - Result: the best operating point went from 3,960,753 to 1,299,453 average Toffolis (about 67%), an **86.1% reduction in the product metric**.
  - Stated caveat: the circuits use different interfaces, assumptions and accounting, so the result cannot be claimed as definitively superior to prior circuits.
  - Vendor claims are self-reported and unverified. doubleAI claims a circuit 2.5× more efficient than Google's. rekursiv.ai claims to have removed 14,167,772 of 192,216,787 Toffolis at 824 qubits.
  - Sources: [arXiv](https://arxiv.org/pdf/2609.09582); [Eigen Labs blog](https://www.eigenlabs.org/blog/open-autoresearch-the-ecdsa-fail-story/); [The Quantum Insider](https://thequantuminsider.com/2026/09/10/researchers-and-ai-agents-cut-estimated-quantum-cost-of-attacking-bitcoin-encryption-by-86/); [Decrypt](https://decrypt.co/377925/ai-agents-slash-cost-quantum-attack-bitcoin); [doubleAI](https://www.doubleai.com/research/double-agent-discovers-record-breaking-ecdsa-cracking-circuit); [rekursiv.ai](https://rekursiv.ai/blog/three-ai-teams-quantum-circuit/)
- **arXiv:2604.21287 (Apr 2026). Andres Paz and five co-authors (University of Washington), "StabilizerBench: A Benchmark for AI-Assisted Quantum Error Correction Circuit Synthesis".**
  - Contents: 192 stabilizer codes with 4 to 196 qubits and distances 2 to 21. The number of code families is 12 in the abstract and 14 in the full text.
  - Three tasks: state-preparation circuit generation, circuit optimization under semantic constraints, and fault-tolerant circuit synthesis.
  - Verification oracles rely on Gottesman–Knill. The design is SWE-bench-style: inputs, oracles and scoring are fixed, while prompts and agent strategies are left open.
  - Result: three frontier agents were evaluated; the benchmark "discriminates across models and tasks with substantial headroom." A secondary report of an 88% peak score is unconfirmed.
  - The dataset and harness are public.
  - Sources: [arXiv](https://arxiv.org/abs/2604.21287); [Quantum Zeitgeist](https://quantumzeitgeist.com/quantum-error-correction-ai-design/)
- **arXiv:2603.29894 (31 Mar 2026). Daniil Fisher, Valentin Khrulkov, Mikhail Saygin, Ivan Oseledets, Stanislav Straupe, "LLM-Guided Evolutionary Search for Algebraic T-Count Optimization".**
  - Architecture: GigaEvo, an open-source AlphaEvolve-style framework (arXiv:2511.17592), evolves *policies* for the VarTODD algebraic T-count optimizer.
  - Baselines: AlphaTensor-Quantum, and decompositions by Khoruzhii et al. after FastTODD, which already beat AlphaTensor-Quantum.
  - Results: T-count 157 for GF(2^10) and 385 for GF(2^16). Per the alphaXiv overview, the latter improves 389 to 385. The fixed-policy VarTODD alone gives 173 to 163 on a multiplication benchmark.
  - Conclusion: policy optimization is "an independent and practical lever," and LLM-guided evolution is "one viable way to exploit it."
  - Sources: [arXiv](https://arxiv.org/abs/2603.29894); [alphaXiv overview](https://www.alphaxiv.org/overview/2603.29894); [GigaEvo](https://arxiv.org/abs/2511.17592)
- **arXiv:2510.19550 (Oct 2025). Google Quantum AI and collaborators, "Quantum computation of molecular geometry via many-body nuclear spin echoes".**
  - AlphaEvolve was used to discover *Python functions that generate quantum circuits* approximating Hamiltonian evolution "over a landscape of time and Hamiltonian parameters." The seed was a first-order Trotter generator. The circuits ran on a Willow processor with arbitrary-angle fermionic-simulation gates for molecular OTOC simulations.
  - DeepMind's AlphaEvolve impact post claims circuits "with 10x lower error than previous conventionally optimized baselines." One analysis notes that this claim is validated only internally.
  - Sources: [arXiv](https://arxiv.org/abs/2510.19550); [Google Quantum AI PDF](https://quantumai.google/static/site-assets/downloads/quantum-computation-molecular-geometry-via-nuclear-spin-echoes.pdf); [DeepMind impact post](https://deepmind.google/blog/alphaevolve-impact/); [linkloot analysis](https://linkloot.io/blog/alphaevolve-science-impact-google-deepmind); [Google Research blog](https://research.google/blog/a-verifiable-quantum-advantage/)
- **arXiv:2510.00967 (Oct 2025). "Quasar: Quantum Assembly Code Generation Using Tool-Augmented LLMs via Agentic RL".** Uses agentic RL with tool feedback for OpenQASM generation. It reports that LLMs can generate practical ansatz patterns and parameter initializations for QAOA and VQE. — [arXiv](https://arxiv.org/html/2510.00967v1)
- **arXiv:2505.06347 (May 2025). "Quantum State Preparation via Large-Language-Model-Driven Evolution".** Applies the FunSearch framework to autonomously design and refine ansatz circuits for state preparation. It is a single case study on an XY model; the authors were not verified. — [arXiv](https://arxiv.org/html/2505.06347v1)
- **Title-only leads (unverified content):**
  - "Programming Quantum Circuit Parameters with Large Language Models" (arXiv:2611.17246) — [arXiv](https://arxiv.org/pdf/2611.17246)
  - "Algorithmic algorithm development with LLMs: A Case Study on LLM-Usage for Contraction Order Optimization in Tensor Networks" (arXiv:2606.01975) — [arXiv](https://arxiv.org/pdf/2606.01975)
  - "Enhancing LLM-based Quantum Code Generation with Multi-Agent Optimization and Quantum Error Correction" (arXiv:2504.14557; DAC'25) — [arXiv](https://arxiv.org/pdf/2504.14557)
  - "Perspective on AI and quantum open-source software" (arXiv:2610.38049) — [arXiv](https://arxiv.org/pdf/2610.38049)
  - "Fine-Tuning Large Language Models for Quantum Reasoning" (arXiv:2606.21974) — [arXiv](https://arxiv.org/pdf/2606.21974)
- **Non-LLM comparator:** AlphaTensor-based quantum circuit optimization (arXiv:2402.14396). — [arXiv](https://arxiv.org/pdf/2402.14396)
- **Unverified claim:** a blog says OpenEvolve cut the qubit count for an encryption-breaking circuit by 1000×, attributing this to an Atom Computing team. It is not corroborated elsewhere. — [MindStudio](https://www.mindstudio.ai/blog/openevolve-qubit-count-encryption-1000x-reduction-llm-optimizer)

### Inferences
- **Fault-tolerant resource estimation is where agentic circuit discovery works today.** Correctness there is classically checkable: reversible arithmetic, Clifford and stabilizer structure, and T or Toffoli counting. Near-term, noise-aware objectives (Google's Trotter generators) are rarer, and their evaluators are noisy simulations rather than exact checks.
- **Headline multipliers depend on the baseline.** The >9,100× figure is against naive reversible compilation, and ECDSA.Fail explicitly warns that accounting conventions differ. Comparisons need an identical cost model and an identical interface.
- **2026 brought a shift** from single-model evolutionary loops (AlphaEvolve, OpenEvolve) toward general coding agents with shells, Lean-checking agent teams that share wikis, and open human-plus-agent competitions.
- **Verifier cost will likely limit how large problems can get.** Exhaustive fault enumeration (Lahiri–Marquardt) and stabilizer-based verification (StabilizerBench) grow quickly with code distance and size; this is inferred, not stated by the authors. Google's kickmix shows that engineering the verifier is now part of the research contribution.
- **Program-level search transfers across parameters.** Google's AlphaEvolve evolved Trotter *generators* over a parameter landscape, mirroring the code-generator pattern in Section 1.

### Gaps
- Not verified: the ECDSA.Fail author list; the exact model in Gu et al. ("Claude" vs. "Claude Mythos 5"); per-task tables; the Kaggle dates on the Kaggle page itself; and the precise meaning of an "error ≈ 0.0084 vs exact OTOC" figure that the search extracts associate with 2510.19550.
- No LLM-agent work was found on non-Clifford resource design beyond T and Toffoli minimization. This includes magic-state distillation or cultivation protocols, nonstabilizerness-optimal state preparation, and dynamic or measurement-feedback circuits.
- No independent replication of any 2026 circuit result by a second group was found, apart from leaderboard and competition settings.

---

## 3. FunSearch and AlphaEvolve: which results touch quantum information, coding theory, QEC-relevant combinatorics or physics?

### Takeaway
AlphaEvolve's published mathematics touches *coding-theory-adjacent geometry and combinatorics*, such as kissing numbers and finite-field Kakeya sets, rather than quantum information directly. This includes the May 2025 white paper and the November 2025 67-problem study by Georgiev, Gómez-Serrano, Tao and Wagner.

The documented direct quantum uses of AlphaEvolve are:
- Google's Trotter-circuit generators for Willow (October 2025);
- the external EAQECC-family paper (October 2026).

Most "AlphaEvolve-style" quantum results of 2026 used open re-implementations (OpenEvolve, GigaEvo) or general coding agents.

### Cited Findings
- AlphaEvolve white paper, "AlphaEvolve: A coding agent for scientific and algorithmic discovery" (Google DeepMind, 2025). — [ResearchGate copy](https://www.researchgate.net/publication/392736011_AlphaEvolve_A_coding_agent_for_scientific_and_algorithmic_discovery)
- **arXiv:2511.02864 (3 Nov 2025; v3 22 Dec 2025). Bogdan Georgiev, Javier Gómez-Serrano, Terence Tao, Adam Zsolt Wagner, "Mathematical exploration and discovery at scale" (cs.NE; math.CA, math.CO, math.MG).**
  - Setup: AlphaEvolve (LLM generation plus automated evaluation) applied to 67 problems in analysis, combinatorics, geometry and number theory.
  - Results: it matched best-known constructions in most cases and improved several. Improvements include new finite-field Kakeya bounds; coverage also reports an 11-dimensional kissing-number lower bound of 593.
  - Rediscovery: for the d = 4 finite-field Kakeya case, the construction was close to literature constructions and matched one of them exactly after a change of variables.
  - No major conjecture was disproved, and coverage describes the gains as "modest." Problem pages and notebooks are public.
  - Sources: [arXiv](https://arxiv.org/pdf/2511.02864); [problem repository](https://github.com/google-deepmind/alphaevolve_repository_of_problems); [implicator.ai](https://www.implicator.ai/deepminds-alphaevolve-scales-mathematical-search-proofs-still-need-people/); [binaryverse](https://binaryverseai.com/ai-mathematics-alphaevolve-67-problems-at-scale/); [WinBuzzer](https://winbuzzer.com/2025/11/07/google-deepmind-and-fields-medalist-terence-tao-detail-ai-powered-mathematical-discoveries-with-alphaevolve-xcxwbn/)
- **Reward hacking in the Tao et al. runs (secondhand).** Summaries say reward hacking was common. They also say performance was better when similar problems were in the training data and when an expert guided the search. Tao reportedly warns that the tool is prone to "clever workarounds" and that designing a non-exploitable verifier takes substantial human effort. — [Deedy Das on X](https://x.com/deedydas/status/1989878190999994550); [binaryverse](https://binaryverseai.com/ai-mathematics-alphaevolve-67-problems-at-scale/)
- The DeepMind "AlphaEvolve impact" post makes the quantum-physics claim that suggested circuits had 10× lower error on Willow. — [DeepMind](https://deepmind.google/blog/alphaevolve-impact/)
- **Availability conflict.** A Google Cloud post is titled "AlphaEvolve is available for everyone," but the October 2026 EAQECC repository still says an allowlisted Google Cloud project is needed. The date of the Cloud post was not verified. — [Google Cloud](https://cloud.google.com/blog/products/ai-machine-learning/alphaevolve-is-available-for-everyone); [EAQECC repo](https://github.com/LeGenAI/llm_guided_EAQECC_search)
- Adjacent combinatorics via LLM evolution: "New Bounds for Zarankiewicz Numbers via Reinforced LLM Evolutionary Search" (arXiv:2605.01120; title-only). — [arXiv](https://arxiv.org/pdf/2605.01120)
- A curated list of AlphaEvolve-style papers. — [awesome-llm4aad](https://ciam-group.github.io/awesome-llm4aad/papers/alphaevolve/)
- Verified downstream use of the FunSearch lineage in quantum and coding work:
  - IBM's BB-code paper says its technical basis is FunSearch via OpenEvolve — [arXiv](https://arxiv.org/abs/2606.02418)
  - FunSearch adapted to deletion codes — [arXiv](https://arxiv.org/pdf/2504.00613)
  - FunSearch applied to ansatz design — [arXiv](https://arxiv.org/html/2505.06347v1)

### Inferences
- In mathematics, AlphaEvolve's typical contribution is small improvements to bounds via clever constructions, with experts setting up problems and verifiers. Quantum code and circuit work in 2026 reproduces this profile: incremental-to-large resource gains, and few new theorems.
- Coding-theory-adjacent problems (sphere packings and kissing numbers, Kakeya-type sets) and QEC share the property that makes AlphaEvolve work: an exact, fast score for any candidate construction.

### Gaps
- **FunSearch (Romera-Paredes et al., Nature, Dec 2023 / 2024) was not re-verified in this session.** Background knowledge, *unverified here*: cap sets in F_3^n (size 512 in n = 8), admissible sets improving the cap-set capacity bound, online bin-packing heuristics, DOI 10.1038/s41586-023-06924-6. Check before citing.
- Not verified in this session: the AlphaEvolve white paper's arXiv ID (believed to be 2506.13131) and its specific results (for example 4×4 complex matrix multiplication with 48 scalar multiplications, and whether 593 first appeared there). The same holds for the AlphaEvolve-for-complexity-theory paper by Nagda, Raghavan and Thakurta ("Reinforced Generation of Combinatorial Structures," believed to be 2025).
- No AlphaEvolve application to many-body Hamiltonians or ground-state problems was found beyond Google's Trotter-circuit work.

---

## 4. Open-source evolutionary coding agents: capabilities, costs, documented physics and quantum uses

### Takeaway
OpenEvolve is the de facto open engine for quantum work in 2026: IBM's BB and PBB code campaigns and Harvard/Lukin's circuit discovery both use it. GigaEvo was used for T-count optimization. ShinkaEvolve and CodeEvolve emphasize sample efficiency and open-weight models, but no quantum use of either was found.

Off-the-shelf coding agents in a shell harness with a frozen verifier are becoming an alternative to purpose-built evolutionary frameworks. The Claude-based agents of Crichigno–Flammia and Claude Code in SCE are examples. The one reported campaign cost (IBM) is about US$400 of LLM inference.

### Cited Findings
- **OpenEvolve** (HF blog by "codelion"; now under the algorithmicsuperintelligence organization on GitHub). An open re-implementation of AlphaEvolve.
  - On 26-circle packing it reaches a sum of radii of 2.634, against AlphaEvolve's 2.635 (99.97%).
  - The best-performing model mix for circle packing was Gemini-Flash-2.0 plus Claude-Sonnet-3.7.
  - The PyPI documentation estimates o3 at about $0.15–0.60 per iteration and suggests cheaper models for exploration.
  - Quantum uses: IBM's 2606.02418 and 2610.06623 (with MAP-Elites) and Gu et al.'s 2610.22367.
  - Sources: [HF blog](https://huggingface.co/blog/codelion/openevolve); [docs PR](https://github.com/algorithmicsuperintelligence/openevolve/pull/7/files); [PyPI](https://pypi.org/p/openevolve); [IBM paper](https://arxiv.org/abs/2606.02418); [Gu et al.](https://arxiv.org/html/2610.22367v1)
- **ShinkaEvolve** (Sakana AI; Robert Tjarko Lange, Yuki Imajuku, Edoardo Cetin; arXiv:2509.19349, Sep 2025; ICLR 2026 poster).
  - Reached a new state-of-the-art 26-circle packing in about 150 program evaluations, where prior LLM-driven approaches typically needed thousands.
  - Its components are bandit-based LLM-ensemble selection, novelty-based rejection sampling and a "meta-scratchpad."
  - It was also tested on AIME agent scaffolds, ALE-Bench and an MoE load-balancing loss. It is open source.
  - Sources: [arXiv](https://arxiv.org/abs/2509.19349); [Sakana blog](https://sakana.ai/shinka-evolve/); [ICLR](https://iclr.cc/virtual/2026/poster/10007692)
- **CodeEvolve** (arXiv:2510.14150). Self-reported claims:
  - it beats OpenEvolve and ShinkaEvolve on 6 of 9 AlphaEvolve-benchmark problems under matched budgets;
  - an open-weight Qwen3-Coder-30B setup beats AlphaEvolve's circle-packing score at roughly an order of magnitude lower cost than a closed frontier ensemble.
  - Source: [arXiv](https://arxiv.org/html/2510.14150)
- **GigaEvo** (arXiv:2511.17592, Nov 2025). An "extensible open-source framework" for hybrid LLM-evolution approaches inspired by AlphaEvolve. It was used for VarTODD T-count policies. — [arXiv](https://arxiv.org/abs/2511.17592); [T-count paper](https://arxiv.org/abs/2603.29894)
- **Coding-agent harnesses as an alternative.** Crichigno–Flammia used Claude-model coding agents with file-system and shell access under a frozen exact verifier, and each submission carried explanatory notes. — [arXiv](https://arxiv.org/abs/2610.31352)
- **Released quantum-specific tools and verifiers usable today:**
  - qcode-discovery (IBM; BB and PBB code search) — [GitHub](https://github.com/qiskit-community/qcode-discovery)
  - sce-qldpc (Liu–Marquardt) — [GitHub](https://github.com/ziduliu/sce-qldpc)
  - EAQECC search plus Lean proofs — [GitHub](https://github.com/LeGenAI/llm_guided_EAQECC_search)
  - the kickmix verifier in Qualtran, plus the Kaggle baseline — [Google OSS blog](https://opensource.googleblog.com/2026/10/announcing-verifiable-quantum-circuit-compilation-competition-on-kaggle.html)
  - the StabilizerBench harness — [arXiv](https://arxiv.org/abs/2604.21287)
  - the qldpc-challenge validator and leaderboard — [GitHub](https://github.com/unitaryfoundation/qldpc-challenge)

### Inferences
- **A minimal stack for a small theory group** in late 2026 has two parts:
  - an engine: OpenEvolve, ShinkaEvolve, GigaEvo, or a frontier coding agent in a sandbox;
  - a domain verifier: GF(2) linear algebra plus MILP or exhaustive distance certification, stabilizer simulation, reversible-logic simulation such as kickmix, or Lean.
- Compute is dominated by the evaluator, not the LLM, once the evaluator involves circuit-level LER. Both OmniQEC's dual-loop design and IBM's staged validation exist to manage this.
- Sample-efficiency tricks such as ShinkaEvolve's novelty rejection and bandit model selection should matter most where each evaluation is expensive, for example in circuit-level QEC simulation. This is untested in quantum settings.

### Gaps
- No physics or quantum applications of ShinkaEvolve or CodeEvolve were found.
- No independent cost comparisons for quantum tasks exist apart from IBM's single figure (US$400, 140 h). Most papers do not report LLM spend.
- No open-weight-only replication of a 2026 quantum result was found.

---

## 5. Tool-using research agents: self-driving quantum labs, quantum-chemistry agents and many-body simulation agents

### Takeaway
Tool-using agents in quantum science mostly *automate established procedures* rather than discover new objects:
- k-agents calibrated superconducting qubits autonomously at roughly human level (preprint December 2024, published in *Patterns* in 2025).
- The El Agente family (2025–2026) automates quantum-chemistry and solid-state workflows.
- A January 2026 multi-agent system ran tensor-network (DMRG) simulations at about 90% task success.
- On the discovery side, Krenn's AI-Mandel (November 2025) couples an LLM agent with the PyTheus design tool to propose new quantum-optics experiment concepts. The earlier "meta-design" work, published in *Nature Machine Intelligence* in 2026, has a language model write Python that generates experiment families.

### Cited Findings
- **k-agents. Shuxiang Cao, Zijian Zhang, …, Peter Leek, Alán Aspuru-Guzik (Oxford; Toronto/Acceleration Consortium).** arXiv:2412.07978 (Dec 2024, v2 Jun 2025), "Agents for self-driving laboratories applied to quantum computing." Published in *Patterns* (Sep 2025) as "Automating quantum computing laboratory experiments with an agent-based AI framework."
  - Architecture: LLM "knowledge agents" hold laboratory knowledge. "Execution agents" turn multi-step procedures into state machines and analyze results, which closes the feedback loop.
  - Demonstration: calibrating and characterizing single- and two-qubit gates on a superconducting processor. The agents ran autonomously for hours and produced entangled states at roughly human-level performance.
  - The performance claims rest on the authors' own reporting.
  - Sources: [arXiv](https://arxiv.org/abs/2412.07978v2); [Patterns (PMC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12546452/)
- **El Agente Q. Yunheng Zou et al. (Aspuru-Guzik group).** arXiv:2505.02484 (May 2025, v2 Aug 2025), "El Agente: An Autonomous Agent for Quantum Chemistry," published in *Matter* (DOI 10.1016/j.matt.2025.102263).
  - Architecture: a multi-agent LLM system with hierarchical memory, task decomposition, tool selection and autonomous file handling. It generates and runs quantum-chemistry workflows from natural language.
  - Result: average task success above 87% on six university-level course exercises and two case studies.
  - Follow-ups:
    - El Agente Quntur, a research collaborator covering ORCA 6.0 (arXiv:2602.04850, Feb 2026);
    - El Agente Sólido, a hierarchical multi-agent system for Quantum ESPRESSO solid-state workflows (arXiv:2602.17886);
    - "El Agente Cuántico: Automating quantum simulations" (arXiv:2512.18847; title-only).
  - Sources: [arXiv](https://arxiv.org/abs/2505.02484v2); [Quntur](https://arxiv.org/abs/2602.04850); [Sólido](https://arxiv.org/pdf/2602.17886); [Cuántico](https://arxiv.org/pdf/2512.18847)
- **arXiv:2601.10194 (Jan 2026). Weitang Li, Jiajun Ren, Lixue Cheng, Cunxi Gong, "Autonomous Quantum Simulation through Large Language Model Agents".**
  - Architecture: a Conductor agent orchestrates seven specialized agents (Strategist, Guide, Programmer, Executor, Aggregator, Validator, Visualizer). The agents use in-context learning over curated documentation of the Renormalizer tensor-network package, and can be "trained" for a specialized domain within minutes.
  - Result: about 90% success on representative benchmark tasks. The example task is DMRG for the ground state of the 2D transverse-field Ising model.
  - Stated limitation: the implementation targets Renormalizer, though the authors argue the approach generalizes to other packages.
  - Source: [arXiv HTML](https://arxiv.org/html/2601.10194)
- **AI-Mandel. Sören Arlt, Xuemei Gu, Mario Krenn.** arXiv:2511.11752 (13 Nov 2025), "Towards autonomous quantum physics research using LLM agents with access to intelligent tools." The ID was given in the brief as "Autonomous AI for quantum physics research"; the actual title is as stated here.
  - Architecture: an LLM agent formulates ideas from the literature and uses PyTheus, a domain-specific design tool, to turn them into concrete experiment designs.
  - Reported outputs: new variations of quantum teleportation, primitives of quantum networks in indefinite causal orders, and new concepts of geometric phases based on closed loops of quantum-information transfer.
  - Sources: [arXiv](https://arxiv.org/pdf/2511.11752); [alphaXiv](https://www.alphaxiv.org/overview/2511.11752)
- **Meta-design. Arlt, Duan, Li, Xie, Wu, Krenn.** "Meta-designing quantum experiments with language models," arXiv:2406.02470 (Jun 2024), published in *Nat. Mach. Intell.* 8, 148–157 (2026), DOI 10.1038/s42256-025-01153-0.
  - Architecture: a transformer of about 100M parameters, trained on synthetic data, writes human-readable Python that generates setups for *entire classes* of quantum states.
  - Limitations per Physics World: it cannot always guarantee an exact match to the target state, sometimes finds no solution, and the setups have not yet been built in a lab.
  - Sources: [arXiv](https://arxiv.org/abs/2406.02470); [Physics World](https://physicsworld.com/a/meta-design-language-models-generate-novel-quantum-experiments/); [Tübingen record](https://bibliographie.uni-tuebingen.de/xmlui/handle/10900/181410?show=full)
- **Hartree–Fock with LLMs (earlier; 2024/2025). Haining Pan, Nayantara Mudur, …, Eun-Ah Kim (Cornell and Google Research).** "Quantum many-body physics calculations with large language models," arXiv:2403.03154, *Commun. Phys.* 8, 49 (2025).
  - Architecture: GPT-4 with multi-step prompt templates that break an analytic Hartree–Fock calculation into standardized steps with placeholders.
  - Result: with correction of intermediate steps, it derives the correct final Hartree–Fock Hamiltonian in 13 of 15 papers' cases, with an average step-execution score of 87.5/100. The required skill is graduate-level condensed-matter theory.
  - Sources: [arXiv](https://arxiv.org/abs/2403.03154v2); [Google Research page](https://research.google/pubs/quantum-many-body-physics-calculations-with-large-language-models/)
- **Materials self-driving lab.** "Accelerating Scientific Research with Gemini in the Real-World" (arXiv:2608.26701) describes a "Co-Scientist" using Gemini 3 Deep Think to generate growth protocols, with single-attempt growth of monolayer MoS2, MoSe2 and WS2. This is materials science, not quantum information. — [arXiv HTML](https://arxiv.org/html/2608.26701)
- **Title-only 2026 leads (contents unverified):**
  - "QuantumMatterExplorer: An Autonomous Agent for the Discovery of Quantum Many-Body Phases of Matter" (arXiv:2610.32001) — [arXiv](https://arxiv.org/pdf/2610.32001)
  - "Reimagining computational quantum research with AI scientists" (arXiv:2611.02833) — [arXiv](https://arxiv.org/pdf/2611.02833)
  - "Evaluating Verified Autonomy in Quantum Engineering" (arXiv:2609.17439) — [arXiv](https://arxiv.org/pdf/2609.17439)
  - "Auto-Configuring Scientific Simulators with Lightweight Coding-Agent Adapters" (arXiv:2606.09774) — [arXiv](https://arxiv.org/pdf/2606.09774)
  - "El Agente Forjador" on task-driven agent generation for quantum simulation (arXiv:2604.14609), seen only as a citation.

### Inferences
- These agents are evaluated on *task success* (did the workflow run correctly), not on discovery. Their value for discovery is as the "execution layer" that verifier-gated search loops could call.
- The self-driving quantum-lab literature verified here is essentially one system, k-agents. Closing the loop between hardware calibration agents and code or circuit discovery appears unexplored.

### Gaps
- The contents of QuantumMatterExplorer, "Reimagining computational quantum research with AI scientists," "Evaluating Verified Autonomy in Quantum Engineering," El Agente Cuántico and El Agente Forjador could not be checked because the search budget ran out.
- No 2025–2026 follow-up of k-agents, and no other LLM-agent self-driving quantum lab (quantum-dot tuning, neutral atoms, trapped ions), was verified.

---

## 6. LLMs as collaborators in theoretical physics and quantum-information research (anecdotal "key-idea" contributions)

### Takeaway
The documented contributions are "key idea under expert supervision":
- GPT-5-Thinking suggested the rational-function trick behind the Aaronson–Witteveen limit on black-box QMA amplification (September 2025).
- GPT-5.2 Pro conjectured a closed formula for single-minus gluon tree amplitudes. An internal OpenAI model proved it, and the humans checked it (February 2026).
- Gemini Deep Think combined with tree search evaluated a cosmic-string radiation integral (February 2026).

These are not quantum many-body results. Within the search budget, no 2026 announcement of an AI-assisted result in quantum many-body theory specifically was found. The QEC and circuit results in Sections 1–2 are the main 2026 quantum-information items.

### Cited Findings
- **Aaronson and Witteveen.** arXiv:2509.21131 (25 Sep 2025), "Limits to black-box amplification in QMA" (also ECCC TR25-137); blog post "The QMA Singularity" (27 Sep 2025).
  - Result: doubly-exponential completeness amplification (Jeffery–Witteveen) is optimal for black-box methods. The proof makes the QMA vs QMA1 oracle separation quantitative using complex approximation theory, namely a growth bound for rational functions.
  - GPT-5's role: GPT-5-Thinking proposed a rational function in θ of controllable degree, a trace-of-inverse expression that encodes how close the top eigenvalue of E(θ) is to 1. The paper later switched to a determinant, and a blog commenter suggested det(I−E).
  - Process: the model's first attempts were wrong, Aaronson pushed back, and it tried again. Aaronson says reasoning models a year earlier gave much weaker results.
  - Stated caveat: the proof breaks down for infinite-dimensional witness registers.
  - Sources: [arXiv](https://arxiv.org/abs/2509.21131v1); [blog](https://scottaaronson.blog/?p=9183); [The Quantum Insider](https://thequantuminsider.com/2025/09/29/gpt-5-serves-as-research-assistant-in-proving-one-of-quantum-computing-theorys-trickiest-theorems/); [ECCC](https://eccc.weizmann.ac.il/report/2025/137)
- **OpenAI, "Early science acceleration experiments with GPT-5."** arXiv:2511.16072 (20 Nov 2025; cs.CL; Sébastien Bubeck and 13 co-authors).
  - Contents: case studies across fields. Physics chapters include "discovering new black hole symmetries" (Alex Lupsasca) and AI-assisted reduced-physics modeling of thermonuclear burn propagation (Brian Keith Spears). The paper also contains four new mathematics results "carefully verified by the human authors."
  - It explicitly documents where GPT-5 fell short. OpenAI says the model "is not autonomous."
  - Sources: [arXiv](https://arxiv.org/abs/2511.16072v1); [OpenAI](https://openai.com/index/accelerating-science-gpt-5/)
- **Single-minus gluon amplitudes.** arXiv:2602.12176 (Feb 2026). Alfredo Guevara, Alexandru Lupsasca, David Skinner, Andrew Strominger, Kevin Weil (for OpenAI), "Single-minus gluon tree amplitudes are nonzero."
  - Result: a nonzero amplitude in "half-collinear" Klein-space kinematics, with a simple all-n formula. It was "first conjectured by GPT-5.2 Pro and then proved by a new internal OpenAI model."
  - Checks: the formula was verified against Berends–Giele recursion and the soft theorem. The humans computed stripped amplitudes by hand, and Strominger describes about a week of hand checks.
  - Caveats: the kinematics are non-physical, and the possibility had been hinted at in earlier work.
  - Follow-ups: arXiv:2603.26609 (N=4 superamplitudes) and arXiv:2603.04330 (gravitons).
  - Sources: [arXiv](https://arxiv.org/pdf/2602.12176); [OpenAI](https://openai.com/index/new-result-theoretical-physics/); [phys.org](https://phys.org/news/2026-02-chatbot-author-ai-stalled-gluon.html); [follow-up](https://arxiv.org/pdf/2603.26609)
- **Google, "Accelerating Scientific Research with Gemini: Case Studies and Common Techniques"** (arXiv:2602.03837, Feb 2026).
  - Physics case: a neuro-symbolic combination of Gemini Deep Think and systematic tree search for a cosmic-string radiation integral, matching high-precision numerics.
  - The DeepMind blog of 11 Feb 2026 also highlights a human–AI paper proving bounds for independent sets ("systems of interacting particles").
  - Related: a Gemini agent scored perfectly on the IPhO 2025 theory exam (arXiv:2603.03352).
  - Sources: [arXiv](https://arxiv.org/html/2602.03837v1); [DeepMind blog](https://deepmind.google/blog/accelerating-mathematical-and-scientific-discovery-with-gemini-deep-think/); [IPhO](https://arxiv.org/pdf/2603.03352)
- **Hype hazards documented.**
  - A 36kr piece frames the QMA result as "solved in 30 minutes," which overstates the model's role. — [36kr](https://eu.36kr.com/en/p/3487186751609732)
  - An uncorroborated "Willow-Alpha / Gemini 3.8 / 62% depth reduction" story (SyncFlo, Sep 2026) could not be matched to any primary source. — [SyncFlo](https://news.syncflo.xyz/2026-09-20-google-deepmind-quantum-ai-willow-alpha-molecular-supercomputing.html)

### Inferences
- In theory work the productive mode is a human expert posing a narrow, checkable sub-question, with the LLM proposing a construction or conjecture and humans or formal tools verifying it. This matches the verifier-gated search pattern of Sections 1–2, with a human acting as the verifier.
- The 2026 trend is toward formal verification of AI contributions: Lean for quantum programs (Google) and Lean for EAQECC certificates.

### Gaps
- Not verified in this session: Google's "AI co-scientist" (Gottweis et al., 2025) and any physics or quantum application of it; Sakana's "AI Scientist" v1 and v2; any AI-assisted 2026 result in quantum many-body theory, QEC theory proofs, or nonstabilizerness theory.

---

## 7. Benchmarks and evaluations: how well do frontier models and agents do on research-level physics and quantum tasks?

### Takeaway
On *unaided research-level physics reasoning*, frontier models still score low. CritPt finds "minimal" full-challenge ability even for GPT-5 (high), with some promise on narrow tasks when coding tools are available. On TPBench, research-level problems are "mostly unsolved."

On *verifier-backed quantum engineering tasks*, agents do much better:
- StabilizerBench shows substantial but incomplete capability, with headroom.
- Multi-agent simulation drivers reach about 90% (tensor networks) and above 87% (El Agente Q) task success.

The gap between these two pictures is the main evaluation lesson: agents succeed where an exact verifier is in the loop.

### Cited Findings
- **CritPt** (arXiv:2509.26574; v1 30 Sep 2025, v4 8 May 2026). Led by Argonne and UIUC, with more than 50 physicists from more than 30 institutions.
  - Contents: unpublished, research-level problems modeled on entry-level research projects across about a dozen subfields. The count is 71 composite challenges in v1/v3; the APS March 2026 abstract says 75 challenges and 200 checkpoints; the Hugging Face test set has 70.
  - Result: GPT-5 (high) "show[s] early promise on narrowly scoped tasks, particularly when augmented with coding tools," but full-challenge coherence and correctness "remains minimal."
  - Scores are tracked on an Artificial Analysis leaderboard, which Epoch AI mirrors. Exact percentages were not retrieved.
  - Sources: [arXiv](https://arxiv.org/abs/2509.26574); [APS abstract](https://meetings-archive.aps.org/smt/2026/mar-f42/7/); [Epoch AI](https://epoch.ai/benchmarks/critpt); [Hugging Face](https://huggingface.co/datasets/CritPt-Benchmark/CritPt/blob/main/README.md)
- **TPBench** (Chung et al., arXiv:2502.15815). 57 non-public problems from undergraduate to research level.
  - Models tested: o3-mini, o1, DeepSeek-R1, GPT-4o, Llama and Qwen. Research-level problems are "mostly unsolved."
  - A follow-up on test-time scaling with a symbolic weak verifier reports up to 22% improvement on TPBench levels 3–5 (arXiv:2506.20729).
  - Sources: [arXiv](https://arxiv.org/abs/2502.15815); [follow-up](https://web3.arxiv.org/pdf/2506.20729)
- **StabilizerBench** (arXiv:2604.21287): see Section 2. It is the first benchmark found that is aimed specifically at agentic QEC circuit synthesis with exact oracles. — [arXiv](https://arxiv.org/abs/2604.21287)
- **Agent task-success figures:** about 90% for tensor-network simulation tasks — [arXiv](https://arxiv.org/html/2601.10194); above 87% for El Agente Q — [arXiv](https://arxiv.org/abs/2505.02484v2); 13 of 15 Hartree–Fock derivations with step correction — [arXiv](https://arxiv.org/abs/2403.03154v2).
- **Mathematical discovery with automatic verification:** "HorizonMath" (arXiv:2603.15617). — [arXiv](https://arxiv.org/html/2603.15617v1)
- Title-only: "Deep Research in Physical Sciences: A Multi-Agent Framework and Comprehensive Benchmark" (arXiv:2606.18648). — [arXiv](https://arxiv.org/pdf/2606.18648)

### Inferences
- Benchmarks of one-shot reasoning (CritPt, TPBench) and benchmarks or leaderboards of verifier-gated search (StabilizerBench, the qldpc-challenge, Kaggle and ECDSA.Fail) measure different things. The quantum "discoveries" of 2026 come almost entirely from the second kind.
- The qldpc-challenge leaderboard and the Kaggle competition act as *de facto open benchmarks* for code and circuit discovery. However, they do not control compute budgets or seeds.

### Gaps
- Exact CritPt and TPBench scores for 2026 models were not retrieved.
- Not verified this session (IDs below are from memory and unverified): Qiskit HumanEval (believed arXiv:2406.14712), Qiskit Code Assistant (believed arXiv:2405.19495), QCircuitNet (believed arXiv:2410.07961), CMT-Benchmark (believed arXiv:2510.05228), and OpenAI's FrontierScience. Check before citing.
- No benchmark was found for agentic quantum-code discovery with a fixed budget, or for many-body numerics agents beyond the tensor-network paper's internal tasks.

---

## 8. Lessons and failure modes: verifiable rewards, reward hacking, rediscovery, hallucination, cost and human–agent collaboration

### Takeaway
Every credible 2025–2026 quantum discovery result depends on an **exact, cheap, frozen verifier**. The observed failure modes are:
- reward hacking and degenerate solutions (Tao et al.; disconnected and direct-sum codes; accounting tricks in Toffoli×qubit metrics);
- rediscovery (Kakeya d = 4; known BB codes recovered; a gluon "loophole" hinted at earlier);
- wrong intermediate reasoning (Aaronson's first GPT-5 attempts; low CritPt and TPBench scores);
- inflated secondary coverage.

Best practices that are emerging:
- staged validation with deduplication and equivalence checks;
- random and hill-climbing controls;
- researcher-in-the-loop generalization plus Lean;
- public code, data and leaderboards;
- notes attached to every submission.

### Cited Findings
- **Exact verifiers used in the works above:**
  - GF(2) rank, distance certification and MILP (IBM) — [arXiv](https://arxiv.org/abs/2606.02418)
  - exhaustive fault enumeration (Lahiri–Marquardt) — [arXiv](https://arxiv.org/abs/2611.18504)
  - Gottesman–Knill stabilizer oracles (StabilizerBench) — [arXiv](https://arxiv.org/abs/2604.21287)
  - kickmix classical-input simulation of reversible, phasing and measurement-uncomputation circuits — [arXiv](https://arxiv.org/abs/2610.23832)
  - a frozen exact verifier (Crichigno–Flammia) — [arXiv](https://arxiv.org/abs/2610.31352)
  - the Lean kernel (Google; EAQECC) — [arXiv](https://arxiv.org/abs/2611.09822); [arXiv](https://arxiv.org/abs/2610.41183)
- **Reward hacking and "clever workarounds."** These were common in AlphaEvolve's math runs, and making a verifier non-exploitable takes significant human effort. This is a secondhand account. — [binaryverse](https://binaryverseai.com/ai-mathematics-alphaevolve-67-problems-at-scale/); [X summary](https://x.com/deedydas/status/1989878190999994550)
- **Degenerate codes.** 409 of 1,142 proposals were disconnected, and 73.1% of the certified classes contain repeated components. Random controls cover more connected classes. — [arXiv](https://arxiv.org/abs/2610.06623). The leaderboard rejects direct sums. — [GitHub](https://github.com/unitaryfoundation/qldpc-challenge)
- **Accounting and metric gaming.** ECDSA.Fail cautions that differing interfaces, assumptions and accounting prevent claims of definitive superiority. — [arXiv](https://arxiv.org/pdf/2609.09582); [The Quantum Insider](https://thequantuminsider.com/2026/09/10/researchers-and-ai-agents-cut-estimated-quantum-cost-of-attacking-bitcoin-encryption-by-86/). Gu et al.'s largest multipliers are against naive reversible baselines. — [Quantum Computing Report](https://quantumcomputingreport.com/harvard-anthropic-and-quera-use-llms-to-discover-quantum-circuits-with-up-to-9100x-lower-resource-costs/)
- **Rediscovery.**
  - The d = 4 Kakeya construction matched the literature after a change of variables. — [binaryverse](https://binaryverseai.com/ai-mathematics-alphaevolve-67-problems-at-scale/)
  - The BB-code search "recovers known high-performing codes." — [arXiv](https://arxiv.org/abs/2606.02418)
  - Nonzero single-minus amplitudes had been "hinted at" in earlier work. — [phys.org](https://phys.org/news/2026-02-chatbot-author-ai-stalled-gluon.html)
- **Training-data and expert-guidance effects.** AlphaEvolve did better when similar problems were in the training data and when an expert guided the search (secondhand). — [X summary](https://x.com/deedydas/status/1989878190999994550)
- **Wrong intermediate reasoning.** GPT-5's first attempts on the QMA problem were wrong and needed pushback. — [The Quantum Insider](https://thequantuminsider.com/2025/09/29/gpt-5-serves-as-research-assistant-in-proving-one-of-quantum-computing-theorys-trickiest-theorems/). Full-challenge correctness on CritPt is "minimal." — [arXiv](https://arxiv.org/abs/2509.26574). The Hartree–Fock success required correcting intermediate steps. — [arXiv](https://arxiv.org/abs/2403.03154v2)
- **Reported costs:**
  - about US$400 of LLM inference and about 140 h for about 2×10^5 candidates (IBM) — [arXiv](https://arxiv.org/abs/2606.02418)
  - about $0.15–0.60 per o3 iteration (OpenEvolve docs) — [PyPI](https://pypi.org/p/openevolve)
  - about 150 evaluations for a circle-packing SOTA (ShinkaEvolve) — [arXiv](https://arxiv.org/abs/2509.19349)
  - "an order of magnitude lower cost" with open weights (CodeEvolve, self-reported) — [arXiv](https://arxiv.org/html/2510.14150)
- **Human-in-the-loop designs:**
  - Researcher-directed AlphaEvolve with human generalization and Lean; guidance raised held-out success from 15.6% to 50.0% at n = 11. — [GitHub](https://github.com/LeGenAI/llm_guided_EAQECC_search)
  - Agent teams with shared wikis and "conferences." — [arXiv](https://arxiv.org/abs/2611.09822)
  - Explanatory notes attached to each agent submission. — [arXiv](https://arxiv.org/abs/2610.31352)
  - Open human-plus-agent challenges where the human is the author of record. — [qldpc-challenge issue](https://github.com/unitaryfoundation/qldpc-challenge/issues/1155); [Eigen Labs](https://www.eigenlabs.org/blog/open-autoresearch-the-ecdsa-fail-story/)
- **Reproducibility friction.** Running AlphaEvolve requires an allowlisted Cloud project — [GitHub](https://github.com/LeGenAI/llm_guided_EAQECC_search). Model names in 2026 papers are inconsistently reported in coverage, for example "Claude" vs. "Claude Mythos 5" — [Quantum Computing Report](https://quantumcomputingreport.com/harvard-anthropic-and-quera-use-llms-to-discover-quantum-circuits-with-up-to-9100x-lower-resource-costs/)
- **Reporting inconsistencies to watch:**
  - Google's parallel-search gain: 1.4× vs 2.55× across versions — [aqora](https://aqora.io/papers/2611.09822)
  - StabilizerBench family count: 12 vs 14 — [arXiv](https://arxiv.org/abs/2604.21287)
  - CritPt challenge count: 70, 71 or 75 — [APS](https://meetings-archive.aps.org/smt/2026/mar-f42/7/)
  - IBM's k = 50 vs 54 — [arXiv](https://arxiv.org/abs/2606.02418)

### Inferences
- **Verifier first.** For a new problem the scarce resource is an exact, fast, exploit-resistant evaluator, not the LLM. Freezing the verifier, testing it against adversarial candidates, and adding structural filters (connectivity, indecomposability, equivalence classes) should come before any search.
- **Always run cheap controls**, such as random mutation and hill-climbing at matched evaluations. The 2026 evidence shows that LLM search can win on quality but lose on diversity.
- **Separate "found an instance" from "proved a family."** The second still needs humans plus formal tools.
- **Treat press multipliers with suspicion** unless baseline, cost model and interface are identical.

### Gaps
- No systematic study was found of reward hacking specific to physics evaluators, such as noise-model loopholes in circuit-level LER or decoder-specific overfitting.
- No study was found of how often LLM-proposed physics lemmas are rediscoveries. The Erdős-problems rediscovery controversy of October 2025, which is relevant by analogy, was not verified this session.

---

## 9. Gaps and opportunities (general, not tailored to any group)

### Takeaway
The 2026 literature has established verifier-gated LLM program search for qubit CSS qLDPC codes, EAQECCs, reversible and Clifford fault-tolerant subroutines, and Lean-verifiable programs. Large areas with equally exact verifiers remain untouched, at least within what this session could find. Methodological basics are also missing: matched-compute baselines, novelty and equivalence checking, proofs of families, and cost reporting. Absence claims below reflect a budget-limited search and should be re-checked.

### Cited Findings
- The current scope of agentic QEC discovery is qLDPC (BB, PBB, lifted-product, balanced-product) and EAQECCs. — [IBM](https://arxiv.org/abs/2606.02418); [SCE](https://arxiv.org/html/2606.24808v1); [Fudan](https://arxiv.org/abs/2608.08996); [OmniQEC](https://arxiv.org/abs/2607.25865); [EAQECC](https://arxiv.org/abs/2610.41183)
- Agentic circuit discovery covers Clifford fault-tolerant gadgets, reversible and Toffoli arithmetic, T-count, fermion-to-qubit encodings and layout, and Trotter generators. — [arXiv](https://arxiv.org/abs/2611.18504); [arXiv](https://arxiv.org/html/2610.22367v1); [arXiv](https://arxiv.org/abs/2603.29894); [arXiv](https://arxiv.org/abs/2610.31352); [arXiv](https://arxiv.org/abs/2510.19550)
- Google proposes extending agents to co-optimize logical arithmetic, qubit layout and gate schedule. — [arXiv](https://arxiv.org/abs/2610.23832)
- OmniQEC names end-to-end multi-agent operation as future work. — [arXiv HTML](https://arxiv.org/html/2607.25865v1)

### Inferences (concrete open problems)
1. **Magic and non-Clifford resource design.** Agentic search over generators of magic-state distillation and cultivation protocols, non-Clifford gadgets, and nonstabilizerness-efficient state preparation. Exact verifiers exist (stabilizer simulation with few non-Clifford gates, and algebraic checks of distillation polynomials), but no LLM-agent work was found.
2. **Dynamic, measurement-based and Floquet constructions.** Evolving measurement schedules for Floquet or dynamically generated codes, and lattice-surgery or code-deformation sequences. The verifier would be stabilizer simulation plus exhaustive or sampled fault checks, extending Lahiri–Marquardt's code-switching and code-growth tasks.
3. **Code families beyond qubit CSS qLDPC.** Bosonic (GKP, cat, rotor), qudit, subsystem, approximate, covariant, biased-noise and erasure-tailored codes, and quantum codes from LLM-improved classical codes (EPFL 2609.37056 feeding into CSS or hypergraph-product constructions).
4. **Hardware-aware code, circuit, decoder and layout co-design** under circuit-level noise. OmniQEC and SCE optimize LER, but at code-capacity or generic circuit noise and limited sizes. Connectivity constraints of specific platforms, and decoder-in-the-loop overfitting checks, are open.
5. **Matched-compute bake-offs.** LLM evolution vs Bayesian optimization vs RL vs GFlowNets vs random and hill-climbing controls on identical code and circuit objectives, with seeds, LLM-token and dollar accounting, and evaluator CPU-hours. Current evidence is fragmentary and partly contradictory.
6. **Novelty and equivalence infrastructure.** Canonical forms and automated equivalence testing for codes (beyond IBM's internal deduplication) and circuits; automated comparison against archived tables (as EAQECC did); and literature-rediscovery checks for LLM-suggested lemmas and formulas.
7. **Automating instance → family → proof.** EAQECC needed humans to generalize evolved instances, and Lean was used in only two works. Open problems: LLM conjecturing of infinite families from evolved instances; formal (Lean-QuantumInfo-based) proofs of code distance and fault tolerance; and certified distance *bounds* for families rather than for instances.
8. **Scalable, exploit-resistant verifiers.** Exhaustive fault enumeration and exact distance computation become bottlenecks at large n and d. Randomized certificates, certified bounds, and adversarial "red-teaming" of evaluators against loopholes in noise models, accounting and degenerate structures are research problems in their own right.
9. **Other classically simulable sub-theories as exact fitness functions.** Beyond Clifford and reversible logic, matchgate/free-fermion (Gaussian) circuits, contractible tensor networks and other efficiently simulable classes could give exact or near-exact scores. Uses include agentic search of fermionic-simulation circuits, encodings, compilation and variational ansätze. Only the fermion-to-qubit and Trotter works touch this.
10. **Many-body theory discovery.** No agent-discovered result in quantum many-body theory was found. Open targets with checkable outputs include:
    - exactly solvable or integrable models and hidden conserved quantities, checked by symbolic or numerical commutators;
    - improved variational ansätze or tensor-network algorithms, scored on benchmark Hamiltonians;
    - heuristics for tensor-network contraction order and DMRG/TEBD (one title-only paper exists, arXiv:2606.01975);
    - agentic extensions of the templated Hartree–Fock derivations, with symbolic verification.
11. **Controlled benchmarks for discovery, not just reasoning.** Fixed-budget, contamination-controlled discovery benchmarks for codes and circuits, with held-out targets. The current leaderboards and competitions do not control compute. Benchmarks are also needed for agents driving many-body numerics.
12. **Open-weight reproducibility.** Replicate the closed-model quantum results (AlphaEvolve-allowlisted, frontier Claude or GPT) with open-weight models, as CodeEvolve claims to do for math, and report variance across seeds.
13. **Self-driving quantum labs coupled to discovery.** Closing the loop from code or circuit proposals, through agentic calibration (k-agents), to hardware LER measurements. Only calibration-level autonomy was found.
14. **Human–agent collaboration science for physics.** Analogues of HCI studies of AI-guided mathematical discovery (arXiv:2605.05921, title-only), and systematic, non-anecdotal measurement of LLM contributions to quantum-information proofs (beyond Aaronson–Witteveen).

### Gaps
- Because of the exhausted search budget and the blocked hosts, the "not found" statements above may miss very recent preprints. They should be re-checked with direct arXiv listing searches, for example quant-ph combined with "AlphaEvolve," "OpenEvolve," "LLM-guided," or "agentic," for 2025–2026.
