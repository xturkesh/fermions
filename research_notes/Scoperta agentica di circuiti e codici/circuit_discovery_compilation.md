# AI-driven discovery, synthesis and optimization of quantum circuits and quantum experiments (outside QEC): an annotated catalogue, focused on 2025 to Nov 2026

*Compiled 2026-11-22. The coordinator's scope update is applied: no tailoring to any group; fermionic circuits are treated as one example among others; 2025 to Nov 2026 work is the top priority and is listed newest first, especially LLM agents, coding agents and AlphaEvolve/FunSearch/OpenEvolve-style search applied to circuits and experiments. QEC codes, decoders and fault-tolerant gadgets stay out of scope and appear only as one-line pointers.*

*Verification tags:*
- **[A]**: content confirmed from abstract-level or primary-source snippets.
- **[T]**: only the title and existence were confirmed.
- **[S]**: numbers come from secondary, news or vendor sources.
- **[U]**: unverified or conflicting.
- **[E]**: earlier (pre-2023) context.
- **n/r**: not retrieved.

*Entry format:* date (from the arXiv ID or venue), authors, title, ID/venue | agent or method | reward or verifier | scale | key result | limitations.

*Process limits:*
- arxiv.org failed DNS resolution in this container, so every arXiv detail comes from search-engine snippets of arXiv and mirror pages.
- The shared web-search budget (200 calls per turn) ran out before the final verification pass. The remaining unverified points are listed under each section's "Gaps".

## 0. (Priority, from the scope update) LLM agents, coding agents and LLM-guided evolutionary search applied to circuit discovery, synthesis, compilation and experiment design (2025 to Nov 2026, newest first)

### Takeaway
In 2026, LLM coding agents became a serious tool for circuit design. The strongest results all rely on exact, classical verifiers for restricted circuit classes:
- Gu, Kemeny, Yelin and Lukin with OpenEvolve: lower Toffoli×qubit cost than the best published design on 7 of 8 arithmetic tasks, up to 39×.
- Google DeepMind and Google Quantum AI: Lean-verified programs, up to 290× better for quantum arithmetic.
- Google Quantum AI: agentic compilation of "classically verifiable" kickmix circuits, used to seed a Kaggle competition.
- qBraid with AlphaEvolve: fermionic encodings of distance 5.

Almost all of these are unreplicated preprints or vendor reports as of Nov 2026.

### Cited Findings
- **[2026-11] Lahiri & Marquardt, "Agentic AI for Discovering Fault-Tolerant Quantum Circuits", arXiv:2611.18504.** *Out-of-scope pointer (FT circuits).*
  - Method: FunSearch-style evolutionary search. Frontier LLMs, with interactive tool execution, write Python programs that define stabilizer circuits, and an automated circuit verifier scores each one.
  - Result: multi-size searches found shallower circuits than a prior reference at some code distances. [A/S]
  - Sources: [arXiv](https://arxiv.org/html/2611.18504); [CASRAI summary](https://casrai.org/news/agentic-ai-discovering-fault-tolerant-quantum-circuits-2026)
- **[2026-11] Google DeepMind + Google Quantum AI, "Discovering formally verified quantum programs with AI", arXiv:2611.09822.**
  - Agent: LLM agents propose quantum programs against formal specifications in a Lean 4 library. The library builds on Mathlib and Physlib's Lean-QuantumInfo and uses quantum Hoare logic.
  - Verifier: machine-checked proofs of correctness and resource cost.
  - Results, each relative to prior baselines: parallel Grover search improved by up to 1.4×; quantum arithmetic improved by 290×; simulation subroutines for fusion, photochemistry and drug metabolism sped up by 11×, 8× and 2.5×.
  - Code is public at google-deepmind/formal-quantum-computation.
  - Limitations: preprint. Snippets disagree on which model was used ("Gemini 4" vs generic "LLM agents") [U].
  - Sources: [arXiv PDF](https://arxiv.org/pdf/2611.09822); [arXiv v1](https://arxiv.org/html/2611.09822v1); [Quantum Computing Report](https://quantumcomputingreport.com/google-deepmind-and-google-quantum-ai-use-lean-4-to-discover-formally-verified-quantum-algorithms/); [AlphaSignal](https://alphasignal.ai/news/google-deepmind-s-ai-agents-cut-quantum-algorithm-resources-by-290x)
- **[2026-10/11] Related formal-verification preprints [T]:**
  - "Lean-QuantumAlg: A Framework for AI Agents to Produce Quantum Algorithms that Humans Can Trust", arXiv:2610.30216. [arXiv](https://arxiv.org/pdf/2610.30216)
  - "The Future of AI-assisted Physics: Discovering, Verifying and Organizing Results via Formalization", arXiv:2610.14059. [arXiv](https://arxiv.org/pdf/2610.14059)
- **[2026-11] "Qubrio: High-Performance Quantum Compilation via Multi-Agent LLM Collaboration", arXiv:2611.02269.**
  - Agent: specialized LLM agents.
  - Verifier: the agents are grounded by deterministic verifiers and simulator feedback.
  - Target: neutral-atom hardware.
  - Stated finding: naive single-pass LLM compilation fails in practice.
  - Scale and numbers: n/r.
  - Source: [arXiv](https://arxiv.org/pdf/2611.02269)
- **[2026-11] "Reimagining computational quantum research with AI scientists", arXiv:2611.02833** (perspective) [T]. [arXiv](https://arxiv.org/pdf/2611.02833)
- **[2026-10] "Automating quantum circuit optimization with AI agents", arXiv:2610.31352.**
  - An autoresearch loop in which a frozen, exact verifier certifies every result.
  - Applied to fermion-to-qubit encoding, Hamiltonian layout and unitary synthesis.
  - Authors and numbers: n/r.
  - Source: [arXiv](https://arxiv.org/html/2610.31352v1)
- **[2026-10] Khattar, Rubin, Gidney, Broughton, Low, Babbush (Google Quantum AI), "Leveraging verifiable structure in quantum circuits to find more efficient compilations with agentic AI", arXiv:2610.23832.**
  - Method: AI agents compile bottleneck subroutines that are built only from reversible classical logic, phasing operations and measurement-based uncomputation. Such circuits can be simulated efficiently on classical inputs. Google calls them "kickmix" circuits, and they are supported by the open-source kickmix library in Qualtran.
  - Verifier: classical simulation.
  - Targets: Reed–Solomon decoding for decoded quantum interferometry, and quantum look-up tables with Givens rotations for chemistry simulation.
  - Output: constrained reference implementations that seed a Kaggle competition. A secondary report puts the prize at $150K [S].
  - Rationale: classical verification of arbitrary quantum circuits is intractable, so agentic results can only be trusted when candidates can be verified efficiently.
  - Sources: [arXiv](https://arxiv.org/abs/2610.23832); [Google Open Source blog](https://opensource.googleblog.com/2026/10/announcing-verifiable-quantum-circuit-compilation-competition-on-kaggle.html); [QCR](https://quantumcomputingreport.com/google-quantum-ai-launches-150k-kaggle-competition-to-compress-fault-tolerant-quantum-circuits-using-ai-agents/)
- **[2026-10-15] Gu, Kemeny, Yelin, Lukin, "Discovering state-of-the-art quantum circuits with large language models", arXiv:2610.22367.**
  - Agent: OpenEvolve (an open-source AlphaEvolve-style evolver). An Anthropic Claude model, named "Claude Mythos 5" in secondary coverage [S], proposes architectures and program rewrites.
  - Verifier: an independent evaluator checks outputs, unchanged inputs, clean scratch qubits and phase. It scores each candidate by Toffoli count × peak live qubits.
  - Scale: 8 reversible-arithmetic tasks drawn from chemistry and number theory.
  - Results: on 7 of 8 tasks, the discovered circuits reduce both Toffoli count and qubit count relative to the best published designs. The largest gain is 39× in Toffoli×qubits. On the 8th task, the search found a reversible modular-inversion architecture that makes the classical "divsteps" algorithm practical.
  - Limitations: the counts are fault-tolerant resource estimates, not hardware runs; the paper is not peer reviewed; secondary sources list different affiliations.
  - One blog headline says "seventeenfold", which conflicts with 39× (possibly a different task or metric) [U].
  - Sources: [arXiv](https://arxiv.org/html/2610.22367v1); [DAIR summary](https://academy.dair.ai/papers/discovering-state-of-the-art-quantum-circuits-with-large-language-models-2610.22367); [MACCUS](https://maccus.fr/en/claude-quantum-circuits/); conflicting figure: [Substack](https://brianroemmele.substack.com/p/the-gate-that-costs-forever-a-forgotten)
- **[2026-09] "Toward Quantum Software Automation: A Quantum-Aware Harness for LLM-Guided Evolution", arXiv:2609.38883.**
  - Argues that the base OpenEvolve prompt tells the LLM little about the quantum workload.
  - Adds static circuit analysis, a quantum-aware prompter, and cheap proxy scores in place of expensive fidelity checks.
  - Source: [arXiv](https://arxiv.org/html/2609.38883)
- **[2026-09] "A Protocol-Guided LLM Agent for Quantum Program Synthesis and Execution", arXiv:2609.20258.**
  - Tasks: VQE, QAOA and quantum kernels.
  - Protocol guidance raised evaluator completion from 80.0% to 99.2%.
  - Source: [arXiv](https://arxiv.org/pdf/2609.20258)
- **[2026-09] "Evaluating Verified Autonomy in Quantum Engineering", arXiv:2609.17439** [T]. [arXiv](https://arxiv.org/pdf/2609.17439)
- **[2026-07] "Hybrid LLM-Guided Search for Quantum Reservoir Architecture Design", arXiv:2607.19506.**
  - Tests whether LLM-guided controllers beat random, evolutionary or Bayesian search under the same evaluation budget, rather than treating the LLM as an oracle.
  - Outcome: n/r.
  - Source: [arXiv](https://arxiv.org/pdf/2607.19506)
- **[2026-07] Isogawa, Okabe, Phadetsuwannukun, Li, Cappellaro, "Agentic AI for Scientific Reasoning in Autonomous Quantum Sensing Experiments", arXiv:2607.25145.**
  - Experiment automation (adjacent to this topic).
  - An LLM agent autonomously selected a single NV center, calibrated its resonance, and measured T2* with Ramsey sequences.
  - Source: [arXiv](https://arxiv.org/abs/2607.25145)
- **[2026-06] Sakka, Mizukami, Mitarai, "An LLM System for Autonomous Variational Quantum Circuit Design", arXiv:2606.13380.**
  - A closed-loop LLM agent iterates circuit designs under explicit design constraints.
  - Evaluated on quantum-ML feature maps and VQE ansatz generation.
  - The authors say that realistic noisy hardware settings remain an open challenge.
  - Sources: [arXiv](https://arxiv.org/html/2606.13380); [abs](https://arxiv.org/abs/2606.13380v1)
- **[2026-06] qBraid with Google Cloud AlphaEvolve, "Evolving Quantum Error-Correcting Encodings for Molecular Simulation", arXiv:2606.25870.** *QEC-adjacent fermionic encoding.*
  - Agent: AlphaEvolve, seeded with a baseline Generalized Superfast Encoding (GSE) program. About 1,500 program variants were scored against an exact verifier.
  - Distance: exact code distance 5 on dense molecular Hamiltonians, where human engineering had reached 3.
  - Logical error and qubits: 3.4–7.9× lower logical error (code-capacity model, exact decoding) and 4.2–5.0× fewer data qubits than standard fault-tolerant compilation paths.
  - Generalization: distance 5 also held on BeH2 and H2O, which the search never saw.
  - Hardware: up to 54% lower logical error than direct physical execution, using Clifford noise reduction and in-circuit stabilizer checks. Most of the advantage disappeared when stabilizer readout was deferred to the end of the circuit.
  - Limitations: the results are the company's own, and circuit-level fault tolerance has not been tested.
  - Sources: [qBraid blog](https://www.qbraid.com/blog-posts/qbraid-uses-alphaevolve-for-quantum-error-correction); [arXiv](https://arxiv.org/pdf/2606.25870); [QCR](https://quantumcomputingreport.com/qbraid-integrates-nvidia-cuda-q-remote-targets-expands-gpu-fleet-and-deploys-google-cloud-alphaevolve-for-error-correction/)
- **[2026-05] Google DeepMind, AlphaEvolve one-year impact report.**
  - AlphaEvolve suggested circuits for complex molecular simulations on Google's Willow processor with about 10× lower error than earlier conventionally optimized baselines.
  - DeepMind says this contributed to first-of-a-kind experimental demonstrations.
  - Caveats: vendor claim, no absolute error rates given, and I found no technical paper [S].
  - Sources: [DeepMind blog](https://deepmind.google/blog/alphaevolve-impact/); [Digital Today](https://www.digitaltoday.co.kr/en/view/54172/google-deepmind-releases-one-year-results-for-alphaevolve-dna-variant-detection-errors-cut-30-percent)
- **[2026-03] Cain, Xu, King, Picard, Levine, Endres, Preskill, Huang, Bluvstein, "Shor's algorithm is possible with as few as 10,000 reconfigurable atomic qubits", arXiv:2603.28627.**
  - The paper estimates about 9,739 atomic qubits, combining qLDPC codes with neutral-atom reconfigurability.
  - Blogs claim OpenEvolve was used to improve the circuits, with "~1,000×" gains. They disagree on attribution (Atom Computing vs Caltech/Oratomic), and I could not confirm this from the paper [U].
  - Sources: [alphaXiv](https://www.alphaxiv.org/abs/2603.28627); [PostQuantum](https://postquantum.com/security-pqc/10000-qubits-shors/); [MindStudio](https://www.mindstudio.ai/blog/openevolve-qubit-count-encryption-1000x-reduction-llm-optimizer)
- **[2026-03] "Generative AI for Quantum Circuits and Quantum Code: A Technical Review and Taxonomy", arXiv:2603.16216** (review). It lists industrial Qiskit code LLMs post-trained with GRPO on quantum verifiable rewards (arXiv:2508.20907, 2025). [arXiv](https://arxiv.org/html/2603.16216)
- **[2026-02] "Quantum circuit generation via test-time learning with large language models", arXiv:2602.03466.**
  - No weight updates. The LLM gets textual feedback on how much a candidate circuit's entanglement improved or worsened, and restarts from the best circuit found so far.
  - Source: [arXiv](https://arxiv.org/html/2602.03466v1)
- **[2026, ID malformed] "Vista: Verifier-in-the-Loop Agentic Reinforcement Learning for Quantum Program Synthesis".**
  - GRPO over groups of candidate circuits.
  - Argues that collapsing all verifier outcomes into one scalar reward can destabilize learning.
  - Presents itself as an advance over the supervised fine-tuning of Agent-Q.
  - The alphaXiv ID "2605.vista" is not a valid arXiv ID [U].
  - Source: [alphaXiv](https://www.alphaxiv.org/abs/2605.vista)
- **[2025-12] "Anubuddhi: A Multi-Agent AI System for Designing and Simulating Quantum Optics Experiments", arXiv:2512.15736** [T]. [arXiv](https://arxiv.org/pdf/2512.15736)
- **[2025-11] "Towards autonomous quantum physics research using LLM agents with access to intelligent tools", arXiv:2511.11752** (cited as Krenn-group work) [T]. [arXiv](https://arxiv.org/pdf/2511.11752)
- **[2025-10] "QUASAR: Quantum Assembly Code Generation Using Tool-Augmented LLMs via Agentic RL", arXiv:2510.00967.**
  - Agent: an LLM post-trained with agentic RL to emit OpenQASM 3.0.
  - Verifier: an external quantum simulator, with a hierarchical reward.
  - Result: a 4B-parameter model reaches 99.31% validity at Pass@1 and 100% at Pass@10, beating GPT-4o, GPT-5 and DeepSeek-V3 on this metric.
  - Sources: [arXiv abs](https://arxiv.org/abs/2510.00967); [HTML](https://arxiv.org/html/2510.00967)
- **[2025-10] "Quantum Agents for Algorithmic Discovery", arXiv:2510.08159.** Proposes combining quantum computing and AI agents for algorithm discovery; I saw only the introduction [T]. [arXiv](https://arxiv.org/pdf/2510.08159)
- **[2025-06] "Advanced For-Loop for QML algorithm search", arXiv:2506.18260.** A FunSearch-inspired multi-agent system that searches program code for quantum-ML algorithms. [arXiv](https://arxiv.org/pdf/2506.18260)
- **[2025-06] Novikov et al., AlphaEvolve white paper, arXiv:2506.13131** (context only; the framework itself is covered by other researchers).
  - LLMs combined with evolutionary search; automatic evaluators feed results back into the next generation; a program database is reused across iterations.
  - Source: [arXiv](https://arxiv.org/html/2506.13131v1)
- **News report [S]:** a transformer-plus-RL framework that emits molecular ground-state circuits in one inference step, with an RL stage that rewards lower computed energies. I did not identify the primary paper. [The Quantum Insider, 2026-07-30](https://thequantuminsider.com/2026/07/30/researchers-ai-can-learn-to-build-quantum-circuits-for-drug-molecules-cutting-design-time-by-orders-of-magnitude/)
- **Out-of-scope pointers (QEC):**
  - IBM's LLM-guided evolutionary framework for QEC code discovery. A secondary source reports 465 bivariate-bicycle candidates. [NAND Research](https://nand-research.com/ibms-llm-guided-evolutionary-framework-for-quantum-error-correction-code-discovery/)
  - "Researcher-in-the-Loop Evolutionary Discovery of Entanglement-Assisted Quantum Code Families", arXiv:2610.41183. It used AlphaEvolve via Vertex AI, with a seed program, a task specification and an exact evaluator, and logged 401 completed candidate searches. [arXiv](https://arxiv.org/pdf/2610.41183)

### Inferences
- The unit of search has shifted. RL agents placing gates one at a time (2019–2025) are giving way to LLM agents that write programs which generate circuits, scored by exact evaluators. Both the Lean pipeline and OpenEvolve/AlphaEvolve work this way.
- The common enabler is choosing a circuit class that classical computers can verify exactly: reversible Toffoli arithmetic, kickmix circuits, stabilizer circuits, proof-checked programs, and exact code distance.
- The biggest headline gains (39×, 290×, 10×) are on logical-level resource counts for future fault-tolerant machines, or on vendor-measured error. None comes from an independent NISQ benchmark.
- Equal-compute comparisons against non-LLM search (simulated annealing, plain genetic programming, SAT/ILP) are almost absent. arXiv:2607.19506 is the clearest attempt, and its outcome was not retrieved.
- Evaluator design is becoming the main lever, as the quantum-aware harness paper and Vista's reward-collapse argument show. Google Quantum AI's Kaggle competition turns verifier design into public infrastructure.

### Gaps
- Model identity is unverified: "Claude Mythos 5" in Gu et al. and "Gemini 4" in arXiv:2611.09822 are reported only by secondary sources.
- The 39× vs "17×" discrepancy for Gu et al. is unresolved.
- Authors and numbers were not retrieved for arXiv:2610.31352 (AI-agent circuit optimization) and Qubrio.
- I found no technical write-up behind AlphaEvolve's "10× lower error" Willow claim.
- Seen only in reference lists, so unverified:
  - "Vibe Calibration", arXiv:2606.22376, described as autonomous bring-up of a 112-qubit superconducting processor by a language agent.
  - Neutral-atom agentic workflows (Jurczak et al.), arXiv:2607.25834.
  - An Aspuru-Guzik-group agent-based framework for quantum-lab automation (Patterns, 2025).
- Agent-Q (2025) was not verified beyond its description in Vista.
- Compute budgets (LLM calls, GPU hours) were not retrieved for most entries.

## 1. RL, AlphaZero and MCTS for circuit optimization, synthesis and compilation

### Takeaway
Learned synthesis and compilation have moved from proofs of concept (Fösel et al. 2021; Moro et al. 2021) to working tools:
- AlphaTensor-Quantum for T-count (Nature Machine Intelligence, 2025), with a 2026 reusability study extending it to general circuits.
- IBM's RL synthesis and routing, shipped in Qiskit's AI transpiler passes: near-optimal Clifford synthesis up to 11 qubits, linear functions up to 9, permutations up to 65, and routing up to 133 qubits.
- ZX-calculus RL agents that generalize from 5-qubit training circuits to 80 qubits.

The 2025–26 trend is model-based planning (AlphaZero/MCTS) and RL for non-Clifford count. Exact learned synthesis still stops at roughly 5–11 qubits, and most comparisons are against heuristics rather than optimal or superoptimizer baselines.

### Cited Findings
- **[2026-10] "AI-Assisted Discovery of CNOT Decompositions for Quantum Permutations on Constrained Topologies", arXiv:2610.37045.** Reports that AlphaCNOT gets a lower mean CNOT count than global transposition constructions on every topology tested, and found a 13-CNOT circuit for a full permutation. [arXiv](https://arxiv.org/html/2610.37045v1)
- **[2026-08] "AlphaClifford: Efficient Clifford Synthesis and Transpilation with Model-based RL", arXiv:2608.18946.** An author-profile blurb claims consistent reductions in total and two-qubit gate counts against state-of-the-art synthesis heuristics for unconstrained Clifford optimization, despite a less expressive gate set. The paper itself was not read [T/U]. [arXiv](https://arxiv.org/pdf/2608.18946)
- **[2026-04] "MonteQ: A Monte Carlo Tree Search Based Quantum Circuit Synthesis Framework", arXiv:2604.19029.** Plain MCTS over circuit-building decisions, apparently with no learned network [T]. [arXiv](https://arxiv.org/pdf/2604.19029)
- **[2026-04] "AlphaCNOT: Learning CNOT Minimization with Model-Based Planning", arXiv:2604.13812.**
  - Method: MCTS on top of RL, so the agent can look ahead over candidate CNOT sequences.
  - Benchmark: linear reversible synthesis against the PMH baseline. The v1 abstract reports up to 32% fewer CNOTs; the current abstract says up to 24.62%.
  - Scale: constrained topologies with up to 8 qubits; outperforms state-of-the-art RL methods.
  - Limitation: a third-party review notes the proposed Clifford extension is untested.
  - Sources: [arXiv](https://arxiv.org/abs/2604.13812); [v1](https://arxiv.org/pdf/2604.13812v1); [Pith review](https://pith.science/paper/2604.13812)
- **[2026-01] Zen, Nägele, Marquardt, "Reusability report: Optimizing T count in general quantum circuits with AlphaTensor-Quantum", Nature Machine Intelligence 8, 113–117 (2026), arXiv:2511.09951.**
  - Extends AlphaTensor-Quantum to general random circuits.
  - A single agent trained on 5- to 8-qubit circuits reduces T-count more than previous methods on a large fraction of circuits.
  - Sources: [arXiv](https://arxiv.org/abs/2511.09951); [FAU CRIS](https://cris.fau.de/publications/356862674)
- **[2025-12] F2, offline RL for compiling Hamiltonian simulation via free-fermionic subroutines, arXiv:2512.08023.** Covered in Section 7. [arXiv](https://arxiv.org/abs/2512.08023)
- **[2025-09] Kremer, Javadi-Abhari, Mukhopadhyay (IBM), "Optimizing the non-Clifford-count in unitary synthesis using Reinforcement Learning", arXiv:2509.21709.**
  - Target: T-count and CS-count of unitaries exactly implementable with Clifford+T or Clifford+CS.
  - Result: close-to-optimal decompositions for two-qubit unitaries with up to 100 T gates, which the authors describe as the largest instances achieved by any method to date.
  - Sources: [arXiv](https://arxiv.org/abs/2509.21709v1); [HTML](https://arxiv.org/html/2509.21709v1)
- **[2025-08] "Unitary Synthesis with AlphaZero via Dynamic Circuits", arXiv:2508.21217.**
  - An AlphaZero network plus MCTS builds circuits that match a target unitary, with dynamic-circuit elements (per the title).
  - The AlphaZero agents beat a plain-MCTS baseline.
  - Authors and scale: n/r.
  - Sources: [arXiv](https://arxiv.org/pdf/2508.21217); [HTML](https://arxiv.org/html/2508.21217v1)
- **[2025-07] IBM Qiskit AI transpiler passes (qiskit-ibm-transpiler).**
  - The passes are in beta. The AIRouting pass does both layout and routing.
  - Local-mode Pauli Network inference arrived in version 0.13.0 (2025-07-03).
  - The cloud Transpiler Service is no longer available as of July 2025, so local use is the main path.
  - IBM claims that combining AI passes with heuristic passes reduced two-qubit gate counts by 42% on average (vendor claim) [S].
  - Sources: [IBM docs](https://quantum.cloud.ibm.com/docs/en/guides/ai-transpiler-passes.md); [IBM blog](https://ibm.com/quantum/blog/ai-transpiler-passes); [release notes](https://quantum.cloud.ibm.com/docs/en/api/qiskit-ibm-transpiler/release-notes)
- **[2025-05] Riu et al., "Reinforcement Learning Based Quantum Circuit Optimization via ZX-Calculus", Quantum 9, 1758 (2025), arXiv:2312.11597.**
  - Agent: PPO with graph neural networks that choose ZX rewrite steps.
  - Training: small 5-qubit Clifford+T circuits with a few tens of gates.
  - Result: consistently improves on the state of the art for this circuit type, up to at least 80 qubits and 2,100 gates.
  - Sources: [Quantum](https://quantum-journal.org/papers/q-2025-05-28-1758/); [arXiv](https://arxiv.org/abs/2312.11597)
- **[2025-04] Mattick et al., "Optimizing Quantum Circuits via ZX Diagrams using Reinforcement Learning and Graph Neural Networks", arXiv:2504.03429.** An RL policy acts directly on ZX graphs and is combined with tree search over rewrite sequences to cut CNOT count. [arXiv](https://arxiv.org/abs/2504.03429)
- **[2025-03] Ruiz et al. (Google DeepMind, Quantinuum and others), "Quantum circuit optimization with AlphaTensor", Nature Machine Intelligence 7, 374–385 (2025), DOI 10.1038/s42256-025-01001-1, arXiv:2402.14396.**
  - Method: deep RL in the AlphaTensor style, exploiting the link between T-count and symmetric tensor rank.
  - Results: beats existing T-count optimizers on arithmetic benchmarks, even without gadgets; using gadgets cuts T-count substantially further.
  - It rediscovered a Karatsuba-like multiplication method for finite fields, and matched the best human designs for Shor-algorithm arithmetic and quantum-chemistry primitives.
  - Limitation: it does not yet handle T-depth.
  - Sources: [DeepMind](https://deepmind.google/research/publications/alphatensor-for-optimizing-quantum-computations/); [arXiv](https://arxiv.org/pdf/2402.14396); [The Quantum Insider](https://thequantuminsider.com/2025/03/21/quantum-optimization-gets-an-ai-boost-with-alphatensor-quantum/)
- **[2025-03] Dubal, Kremer et al. (IBM), RL synthesis of Pauli Network blocks, arXiv:2503.14448.** Packaged as a Qiskit transpiler pass, with an average 20% improvement in two-qubit gate count. [arXiv](https://export.arxiv.org/pdf/2503.14448)
- **[2025-02] "Quantum Circuit Design using a Progressive Widening Enhanced Monte Carlo Tree Search", arXiv:2502.03962** [T]. [arXiv](https://arxiv.org/abs/2502.03962)
- **[2024-08] "A Monte Carlo Tree Search approach to QAOA: finding a needle in the haystack", arXiv:2408.12648** [T]. [arXiv](https://arxiv.org/pdf/2408.12648)
- **[2024-05, v2 2025-02] Kremer et al. (IBM Quantum), "Practical and efficient quantum circuit synthesis and transpiling with Reinforcement Learning", arXiv:2405.13196.**
  - Scope: RL synthesis of Clifford, Linear Function and Permutation circuits.
  - Synthesis scale: near-optimal up to 11 (Clifford), 9 (linear function) and 65 (permutation) qubits.
  - Routing: significant cuts in two-qubit depth and count up to 133 qubits, compared with heuristics such as SABRE.
  - The authors say it improves on Qiskit's heuristics without a large runtime cost.
  - Sources: [arXiv](https://arxiv.org/abs/2405.13196); [v2](https://arxiv.org/html/2405.13196v2)
- **[2024-04] Rietsch et al. (Fraunhofer IIS), "Unitary Synthesis of Clifford+T Circuits with Reinforcement Learning", arXiv:2404.14865.**
  - Method: Gumbel AlphaZero, building circuits one gate at a time.
  - Scale: exact synthesis of a subset of Clifford+T unitaries up to 5 qubits, generated from random circuits with up to 60 gates.
  - Speed: about 30 s average inference on one GPU.
  - Comparison: beats QuantumCircuitOpt and MIN-T-SYNTH at higher qubit counts.
  - Sources: [arXiv](https://arxiv.org/pdf/2404.14865); [v4](https://arxiv.org/html/2404.14865v4)
- **[2023-07, OOPSLA 2024] Li et al., "Quarl: A Learning-Based Quantum Circuit Optimizer", arXiv:2307.10120.**
  - Method: RL chooses which Quartz transformations to apply, over the Nam gate set (CNOT, X, H, Rz). The action space is split into two parts, and the state is represented with graph neural networks.
  - Claim: significantly outperforms existing optimizers on almost all benchmarks (authors' claim).
  - Sources: [arXiv](https://arxiv.org/pdf/2307.10120); [OOPSLA version](https://www.wisdom.weizmann.ac.il/~padon/quarl-oopsla2024.pdf)
- **[2023-11] Reward shaping and related RL work:**
  - Moflic & Paler, exponential reward functions (NANOARCH '23). [arXiv](https://arxiv.org/pdf/2311.12509)
  - Moflic & Paler, letting cost temporarily "explode" before optimizing. [arXiv](https://arxiv.org/pdf/2311.12498)
  - "Optimizing ZX-Diagrams with Deep Reinforcement Learning", arXiv:2311.18588. [arXiv](https://arxiv.org/pdf/2311.18588)
  - Trenkwalder et al., RL compilation of product-formula Hamiltonian simulation, arXiv:2311.04285: about 12% fewer gates than the second-best method tested. [aggregator](https://awesomepapers.io/quantum-computing/papers/2311.04285)
- **[2024] Honest-assessment context:** a comprehensive review of circuit optimization notes that practical use of RL is held back by computational cost, scaling issues, and sensitivity to noise and hardware constraints. [arXiv:2408.08941](https://arxiv.org/pdf/2408.08941)
- **[E, 2021] Fösel, Niu, Marquardt, Li, "Quantum circuit optimization with deep reinforcement learning", arXiv:2103.07585.**
  - Attribution check: the fourth author is Li Li, not Sieberer.
  - Agent: a deep convolutional agent that learns generic strategies for any gate set and connectivity, with a free choice of optimization target.
  - Result: on 12-qubit random circuits, an average 27% depth reduction and 15% gate-count reduction.
  - Sources: [arXiv](https://www.arxiv.org/abs/2103.07585); [PDF](https://arxiv.org/pdf/2103.07585)
- **[E, 2021] Moro, Paris, Restelli, Prati, "Quantum compiling by deep reinforcement learning", Communications Physics 4, 178 (2021), arXiv:2105.15048.**
  - Deep RL approximates single-qubit unitaries as gate sequences, learning a general strategy in one pre-compilation step.
  - Examples: small rotations of π/128 and the Harrow–Recht–Chuang efficiently universal gate set.
  - Reduces execution time enough to make real-time compiling plausible.
  - Sources: [Nature](https://www.nature.com/articles/s42005-021-00684-3); [ar5iv](https://ar5iv.labs.arxiv.org/html/2105.15048)
- **[E, 2019] "Global optimization of quantum dynamics with AlphaZero deep exploration", arXiv:1907.05672.** AlphaZero-style search for control pulse sequences. [arXiv](https://arxiv.org/pdf/1907.05672)

### Inferences
- Synthesis has deterministic, fully known dynamics, which suits model-based planning. That likely explains the 2024–26 move from model-free RL to AlphaZero/MCTS variants (Rietsch, AlphaCNOT, AlphaClifford, dynamic-circuit AlphaZero).
- Scale depends on how structured the circuit class is:
  - Permutations (65 qubits) and routing (133) scale far beyond Clifford (11), linear functions (9) or Clifford+T exact synthesis (5).
  - Rewrite-based ZX agents generalize best (5 → 80 qubits), because each action preserves the circuit's semantics and needs no global verification.
- T-count and non-Clifford-count optimization (AlphaTensor-Quantum, Kremer 2025, Rietsch) amounts to minimizing magic (nonstabilizerness), the resource fault-tolerant machines pay for.
- Industrial deployment is limited to IBM's beta AI passes, and the 42% figure is a vendor number for AI passes combined with heuristics.
- Baselines are usually heuristics (SABRE, PMH, default Qiskit/Tket) rather than optimal SAT/ILP synthesis or strong superoptimizers. Claims of being "state of the art" should be read with that in mind.

### Gaps
- I could not find "Optimizing Quantum Circuits, Fast and Slow" (recalled as Xu, Molavi, Tannu, Albarghouthi, ASPLOS 2025, "GUOQ") by search. Its reported result that a simple rewrite-plus-resynthesis search beats learned optimizers is therefore unverified.
- AlphaTensor-Quantum's largest qubit counts and per-benchmark T-counts were not retrieved.
- Authors and scale were not retrieved for arXiv:2508.21217, AlphaClifford and MonteQ.
- The AlphaCNOT headline figure differs between versions (32% vs 24.62%).

## 2. Quantum architecture search (QAS) for variational algorithms: RL, MCTS, differentiable, evolutionary and LLM-based, and honest assessments

### Takeaway
RL, differentiable and evolutionary QAS can find compact ansätze for small chemistry and QML problems; the curriculum-RL method CRLQAS (ICLR 2024) is the best-documented example. The field's own benchmark papers show that gains over simple baselines are not established at scale:
- QAS-Bench (ICML 2023)
- BenchRL-QAS (2025)
- the surrogate-benchmark paper (2025)
- "Energy Accuracy Is Not Enough" (2026)
- an equal-budget comparison of LLM search against random search (2026)

LLM-driven QAS appeared in 2026 without a demonstrated advantage.

### Cited Findings
- **[2026-09] "Evaluation-efficient quantum architecture search with ZX-calculus-based topological reuse", arXiv:2609.29098** [T]. [arXiv](https://arxiv.org/pdf/2609.29098)
- **[2026-07] "DreamQAS: Learning a Decision-Useful World Model for VQE-Efficient Quantum Architecture Search", arXiv:2607.29491** (model-based QAS) [T]. [arXiv](https://arxiv.org/pdf/2607.29491)
- **[2026-07] "Hybrid LLM-Guided Search for Quantum Reservoir Architecture Design", arXiv:2607.19506.** Compares LLM controllers against random, evolutionary and Bayesian search at equal budget (see Section 0). [arXiv](https://arxiv.org/pdf/2607.19506)
- **[2026-07] "Energy Accuracy Is Not Enough: A Structure-Aware Benchmark and Evaluation Protocol for Quantum Architecture Search", arXiv:2607.04845** [T]. [arXiv](https://arxiv.org/pdf/2607.04845)
- **[2026-06] Sakka, Mizukami, Mitarai, LLM agent for ansatz and feature-map design, arXiv:2606.13380** (see Section 0). [arXiv](https://arxiv.org/html/2606.13380)
- **[2026-04] "Recent Advances in Quantum Architecture Search", arXiv:2604.26289** (review) [T]. [arXiv](https://arxiv.org/pdf/2604.26289)
- **[2026-01] "Noise-Aware Quantum Architecture Search Based on NSGA-II Algorithm", arXiv:2601.10965** [T]. [arXiv](https://arxiv.org/pdf/2601.10965)
- **[2026] Review of ansatz architectures (MDPI Processes).** Covers RL and evolutionary QAS (DDQN-guided exploration, AQEA-QAS with diversity-preserving evolution). [DOI](https://doi.org/10.3390/pr14132095)
- **[2025-11] Three QAS preprints [T]:**
  - "Distributed quantum architecture search using multi-agent reinforcement learning", arXiv:2511.22708. [arXiv](https://arxiv.org/pdf/2511.22708)
  - "QASER: Breaking the Depth vs. Accuracy Trade-Off for Quantum Architecture Search", arXiv:2511.16272. [arXiv](https://arxiv.org/pdf/2511.16272)
  - "Hybrid Action Reinforcement Learning for Quantum Architecture Search", arXiv:2511.04967. [arXiv](https://arxiv.org/pdf/2511.04967)
- **[2025-07] "BenchRL-QAS: Benchmarking Reinforcement Learning Algorithms for Quantum Architecture Search", arXiv:2507.12189 (AAAI Symposium 2025).** A unified framework for comparing RL algorithms on QAS for 2- to 8-qubit systems. Sources: [AAAI](https://ojs.aaai.org/index.php/AAAI-SS/article/view/36906); [arXiv](https://arxiv.org/pdf/2507.12189)
- **[2025-06] "Benchmarking Quantum Architecture Search with Surrogate Assistance", arXiv:2506.06762.** Findings n/r [T]. [arXiv](https://arxiv.org/pdf/2506.06762)
- **[2025-05] "Scalable Quantum Architecture Search via Landscape Analysis", arXiv:2505.05380.** Reports 1.41× better classical memory use and 25.23× shorter search time than the baselines DQAS and TF-QAS (TF-QAS uses random circuits). [arXiv](https://arxiv.org/html/2505.05380v1)
- **[2024-11, v3 2026-03] "Reinforcement learning with learned gadgets to tackle hard quantum problems on real hardware", arXiv:2411.00230.**
  - An RL agent builds a state-preparation circuit gate by gate, with the Hamiltonian's energy expectation as the cost.
  - Includes runs on real hardware, but ground-state tests use small systems such as a 3-qubit transverse-field Ising model.
  - Source: [arXiv](https://arxiv.org/html/2411.00230)
- **[2024-07] "Differentiable Quantum Architecture Search in Asynchronous Quantum Reinforcement Learning", arXiv:2407.18202** [T]. [arXiv](https://arxiv.org/pdf/2407.18202)
- **[2024-06] "Application of ZX-calculus to Quantum Architecture Search", arXiv:2406.01095** [T]. [arXiv](https://arxiv.org/pdf/2406.01095)
- **[2024-06] Martyniuk, Jung, Paschke, "Quantum Architecture Search: A Survey", IEEE QCE 2024, DOI 10.1109/QCE60285.2024.00198, arXiv:2406.06210.** Reviews automated design of parameterized quantum circuits. The authors are not Kundu or Dunjko. Sources: [arXiv](https://arxiv.org/abs/2406.06210v1); [Fraunhofer](https://publica.fraunhofer.de/handle/publica/509801)
- **[2024-02] Patel, Kundu, Ostaszewski, Bonet-Monroig, Dunjko, Danaci, "Curriculum reinforcement learning for quantum architecture search under hardware errors" (CRLQAS), ICLR 2024, arXiv:2402.03500.**
  - Agent: curriculum RL with a 3D architecture encoding, restrictions on the environment dynamics, and an episode-halting scheme that favors shorter circuits.
  - Result: outperforms existing QAS algorithms on quantum-chemistry tasks, in both noiseless and noisy simulation.
  - Motivation: the authors note that the effect of noise on architecture search is poorly understood.
  - Limitation: a secondary review says better simulation efficiency or hardware is needed for larger noisy problems.
  - Sources: [ICLR](https://iclr.cc/virtual/2024/poster/17697); [arXiv](https://arxiv.org/abs/2402.03500v1); [review](https://liner.com/review/curriculum-reinforcement-learning-for-quantum-architecture-search-under-hardware-errors)
- **[2024-03] "Distributed quantum architecture search", arXiv:2403.06214** [T]. [arXiv](https://arxiv.org/pdf/2403.06214)
- **[2023, ICML] Lu et al., "QAS-Bench: Rethinking Quantum Architecture Search and A Benchmark", PMLR v202.** A public benchmark of 900 random circuits and 400 random unitaries, evaluating six baseline algorithms including brute-force search. [PMLR](https://proceedings.mlr.press/v202/lu23f/lu23f.pdf)

### Inferences
- From 2023 onward, the field's self-assessment has moved from "new search algorithm X" to benchmarks and evaluation protocols (QAS-Bench, BenchRL-QAS, surrogate benchmarks, the structure-aware 2026 protocol). That suggests the community sees a reproducibility and baseline problem.
- Documented QAS gains are mostly fewer CNOTs and less depth for small molecules (around 8 qubits or fewer in the benchmarks retrieved). Hardware evidence is limited to toy sizes, such as the 3-qubit Ising model in arXiv:2411.00230.
- LLM-based QAS (2026) currently stands where RL-QAS stood around 2021: plausible, with no demonstrated advantage over equal-budget random or evolutionary search.

### Gaps
- Ostaszewski, Trenkwalder, Masarczyk, Scerri, Dunjko, "Reinforcement learning for optimization of variational quantum circuit architectures" (NeurIPS 2021, arXiv:2103.16089, from memory) was not re-verified, because the search budget ran out. Its LiH and H2 results are therefore not reported here.
- I found no paper explicitly showing that random search matches RL-QAS. QAS-Bench includes brute force, but its ranking results were not retrieved.
- No QAS result beyond about 12 qubits on chemistry with a fair baseline was confirmed.

## 3. Generative models for circuits: diffusion, transformers/GPT, GFlowNets, flow matching, and LLMs fine-tuned to emit circuits

### Takeaway
Generative circuit models have progressed through several families:
- discrete-gate diffusion (Nature Machine Intelligence, 2024);
- multimodal discrete-plus-continuous diffusion transformers (CirDiT, 257M parameters, 2025–26);
- GFlowNet samplers (FlowQ-Net 2025, QFlowNet 2026, FlowMeas 2026, Grinch 2026);
- LLMs post-trained with simulator-verified rewards (QUASAR 2025; Qiskit models trained with GRPO, 2025).

GPT-QE/GQE remains the main generative approach tied to ground states: it beats CCSD in the N2 bond-breaking regime and has run on 127-qubit IBM hardware. Scale and correctness guarantees remain the bottlenecks, and I found no flow-matching model that generates circuits.

### Cited Findings
- **[2026-10] "Mamba Quantum Eigensolver (MQE)", arXiv:2610.29779.** Notes that the transformer's quadratic attention cost is a bottleneck for large active spaces in GPT-QE-style models, and proposes a state-space-model alternative. [arXiv](https://arxiv.org/pdf/2610.29779)
- **[2026-09] "Bridge of Ψ's: Quantum Circuit Optimization with Schrödinger Bridges", arXiv:2609.25947** [T]. [arXiv](https://arxiv.org/pdf/2609.25947)
- **[2026-09] "Bottom-Up Design of Quantum Optical Experiments Using Discrete Generative Models" (Grinch), arXiv:2609.08073.** A GFlowNet-style model that samples optical circuits in proportion to a reward, building circuit graphs edge by edge from reward evaluations with no training dataset. Code: github.com/ishume94/Grinch. [arXiv](https://arxiv.org/html/2609.08073)
- **[2026-08] "Generative Learning for Quantum Measurement Design" (FlowMeas), arXiv:2608.11396.** A GFlowNet learns a distribution over ensembles of shallow, hardware-compatible measurement circuits, built in parallel under connectivity and depth masks. [arXiv](https://arxiv.org/html/2608.11396)
- **[2026-07] News report [S]:** a transformer-plus-RL model emits molecular ground-state circuits in one inference step, with RL rewarding lower computed energies. Primary paper not identified. [The Quantum Insider](https://thequantuminsider.com/2026/07/30/researchers-ai-can-learn-to-build-quantum-circuits-for-drug-molecules-cutting-design-time-by-orders-of-magnitude/)
- **[2026-06] LLM reasoning preprints [T]:**
  - "Fine-Tuning Large Language Models for Quantum Reasoning", arXiv:2606.21974. [arXiv](https://arxiv.org/pdf/2606.21974)
  - "Aligning Quantum Operators with Large Language Models", arXiv:2606.13811. [arXiv](https://arxiv.org/pdf/2606.13811)
- **[2026-05] "Physics Guided Generative Optimization for Trotter Suzuki Decomposition", arXiv:2605.13268.** A generative approach to Trotter–Suzuki decompositions that cites multimodal diffusion (D3PM + DDPM) as independent related work. [arXiv](https://arxiv.org/pdf/2605.13268)
- **[2026-03] "QFlowNet: Fast, Diverse, and Efficient Unitary Synthesis with Generative Flow Networks", arXiv:2603.03045.** A GFlowNet paired with a transformer decomposes unitaries into gate sequences, starting from an empty circuit and adding gates. [arXiv](https://arxiv.org/html/2603.03045)
- **[2026-03] Review: "Generative AI for Quantum Circuits and Quantum Code: A Technical Review and Taxonomy", arXiv:2603.16216.** [arXiv](https://arxiv.org/html/2603.16216)
- **[2026-02] LLM test-time learning with entanglement-based textual feedback, arXiv:2602.03466** (see Section 0). [arXiv](https://arxiv.org/html/2602.03466v1)
- **[2025-10] "FlowQ-Net: A Generative Framework for Automated Quantum Circuit Design", arXiv:2510.26688** (GFlowNet) [T]. [arXiv](https://arxiv.org/html/2510.26688v1)
- **[2025-10] QUASAR, arXiv:2510.00967.** RL post-training with a simulator verifier and hierarchical reward; a 4B model reaches 99.31% validity at Pass@1 (see Section 0). [arXiv](https://arxiv.org/abs/2510.00967)
- **[2025-09] "Discrete Flow-Based Generative Models for Measurement Optimization in Quantum Computing", arXiv:2509.15486** [T]. [arXiv](https://arxiv.org/pdf/2509.15486)
- **[2025-08] Qiskit code models post-trained with quantum verifiable rewards (GRPO), arXiv:2508.20907.** Described in IBM docs and in the 2026 review. Sources: [IBM docs](https://quantum.cloud.ibm.com/docs/en/guides/qiskit-code-assistant); [review](https://arxiv.org/html/2603.16216)
- **[2025-06, v2 2026-02, v3 2026-04] Fürrutter, Chandani, Hamamura, Briegel, Muñoz-Gil, "Synthesis of discrete-continuous quantum circuits with multimodal diffusion models", arXiv:2506.01666.**
  - Method: two coupled diffusion processes, one for discrete gate choice and one for continuous gate parameters, in a Circuit Diffusion Transformer (CirDiT).
  - The alphaXiv overview gives a size of 257M parameters.
  - Developed with NVIDIA's Quantum Algorithm Engineering team.
  - Motivation: earlier ML synthesis models were restricted to discrete gate sets.
  - Sources: [arXiv](https://arxiv.org/abs/2506.01666); [v3](https://arxiv.org/html/2506.01666v3); [alphaXiv](https://alphaxiv.org/overview/2506.01666v1); [genQC](https://florianfuerrutter.github.io/genQC/webpage/research.html)
- **[2025] Laird et al., "Generative Flow Networks with Parameterized Quantum Circuits".** Uses a quantum circuit as the GFlowNet policy: 6 layers and 151 parameters, against 9,741 for a classical GFlowNet. It matched the classical model on simple object generation but trained more slowly. This is not circuit design itself. [PDF](https://s2.smu.edu/~mitch/ftp_dir/pubs/dcas25d.pdf)
- **[2024, NeurIPS workshop] "GFlowNets for Hamiltonian decomposition in groups of compatible operators".** [NeurIPS](https://nips.cc/virtual/2024/100008)
- **[2024-10] "QCircuitNet: A Large-Scale Hierarchical Dataset for Quantum Algorithm Design", arXiv:2410.07961.** Renamed QCircuitBench in later versions; NeurIPS 2025 poster.
  - Format: Qiskit code converted to OpenQASM 3.0, with a problem description, circuit code, classical post-processing and a verification function per task.
  - Tasks range from oracle construction (Grover, Simon) to algorithm design up to the Generalized Simon's Problem.
  - Benchmarked LLMs include GPT-3.5-turbo and GPT-4.
  - The authors argue that lack of training data, not code-generation skill, is the main obstacle.
  - Sources: [v1](https://arxiv.org/html/2410.07961v1); [NeurIPS 2025](https://neurips.cc/virtual/2025/poster/121681); [v3](https://arxiv.org/html/2410.07961v3)
- **[2024-05] Dupuis et al. (IBM), "Qiskit Code Assistant: Training LLMs for generating Quantum Computing Code", arXiv:2405.19495.**
  - Base: an IBM Granite 8B code model (pre-trained on 1.6T tokens of code in 116 languages), fine-tuned on Qiskit code as granite-8b-qiskit.
  - The authors argue that training on code alone is not enough; the model also needs basic quantum-computing knowledge.
  - Evaluated on Qiskit HumanEval, an execution-based benchmark. IBM's claim to have "the best model for usable Qiskit code" is its own.
  - The current assistant runs granite-3.3-8b-qiskit, and Qwen2.5-Coder-14B-Qiskit and mistral-small-3.2-24b-qiskit are offered as GGUF models.
  - Sources: [arXiv](https://arxiv.org/abs/2405.19495v1); [IBM blog](https://www.ibm.com/quantum/blog/qiskit-code-assistant); [IBM docs](https://quantum.cloud.ibm.com/docs/en/guides/qiskit-code-assistant)
- **[2024-02] Apak, Bandic, Sarkar, Feld, "KetGPT – Dataset Augmentation of Quantum Circuits Using Transformers", ICCS 2024, DOI 10.1007/978-3-031-63778-0_17, arXiv:2402.13352.**
  - A transformer generates OpenQASM circuits that look like real algorithm circuits rather than random ones.
  - Validation: manual inspection with Qiskit execution, a transformer-based classifier, and structural analysis.
  - Data (per PennyLane): trained on 1,112 real algorithms from MQT Bench; the released set has 1,000 generated circuits with 2–117 qubits and 6–903 gates.
  - Motivation: existing datasets are too small and narrow, and random circuits do not reflect real algorithm structure.
  - Sources: [arXiv](https://arxiv.org/abs/2402.13352); [TU Delft](https://research.tudelft.nl/en/publications/ketgpt-dataset-augmentation-ofquantum-circuits-using-transformers/); [PennyLane](https://www.pennylane.ai/datasets/ketgpt)
- **[2024-01, v2 2025-09] Nakaji et al., "The generative quantum eigensolver (GQE) and its application for ground state search", arXiv:2401.09253.**
  - Agent: GPT-QE, a decoder-only, GPT-2-style transformer that samples sequences of operators from a fixed pool.
  - Reward: the energies of the generated trial states.
  - Results: beats CCSD in the strong bond-dissociation regime of N2. Run on IBM's 127-qubit ibm_kawasaki with error mitigation, it closely matched exact results.
  - The model can be pre-trained without quantum hardware. The authors claim the number of circuit runs per optimization step does not scale with parameter count or pool size.
  - Sources: [arXiv v2](https://arxiv.org/html/2401.09253v2); [abs](https://www.arxiv.org/abs/2401.09253); [HF](https://huggingface.co/papers/2401.09253)
- **[2023-11, Nature Machine Intelligence 2024] Fürrutter, Muñoz-Gil, Briegel, "Quantum circuit synthesis with diffusion models" (arXiv:2311.02041 from memory).** Its publication in Nature Machine Intelligence (2024) is confirmed through the 2025 follow-up's summaries; details are listed under Gaps. [genQC page](https://florianfuerrutter.github.io/genQC/webpage/research.html); [follow-up](https://arxiv.org/abs/2506.01666)
- **Flow matching:** the only flow-matching result found was Quantum MeanFlow, which is single-step generative sampling on NISQ hardware (MNIST, IBM devices), not circuit generation. [arXiv:2609.02186](https://arxiv.org/pdf/2609.02186)

### Inferences
- The field is converging on two templates:
  - Diffusion or GFlowNet models over circuit tokens or graphs, conditioned on a target such as a unitary or an entanglement property.
  - LLMs post-trained with RL against a simulator verifier, where validity is close to solved (QUASAR) but optimality and scale are not.
- GFlowNets are attractive when diversity matters, for example sampling many hardware-compatible measurement ensembles (FlowMeas) or many optical designs (Grinch). That fits discovery tasks better than single-optimum synthesis.
- GPT-QE's real-hardware run supports generative ansatz construction for ground states. The 2026 MQE paper shows that sequence length (active-space size) is the binding constraint.

### Gaps
- Not retrieved for the Nature Machine Intelligence 2024 diffusion paper: its qubit counts, gate set, entanglement-generation accuracy and unitary-compilation accuracy.
- The verifier, scale and results of Agent-Q and Vista were not confirmed.
- No flow-matching circuit generator was found; whether one exists remains open.
- Benchmark numbers for QCircuitBench and Qiskit HumanEval (pass rates) were not retrieved.

## 4. Automated design of quantum-optics experiments: MELVIN, projective simulation, THESEUS/PyTheus, meta-design with language models, and 2025–2026 results

### Takeaway
The Krenn-group line remains the reference: MELVIN (2016), projective-simulation RL (PNAS 2018), THESEUS (PRX 2021), PyTheus with 100 discoveries (Quantum 2023), and meta-design with language models (Nature Machine Intelligence, Feb 2026). Meta-design produced generalizing "meta-solutions" for 20 classes of states and two previously unknown construction rules.

2025–26 adds:
- generative GFlowNet-style designers (Grinch);
- differentiable structured-light discovery (esQueranto);
- high-probability heralded entanglement schemes;
- extraction of analytic NOON families from machine-discovered graphs;
- LLM multi-agent designers (Anubuddhi).

Most of these are preprints that I could confirm only at title or abstract level.

### Cited Findings
- **[2026-11] "From Machine-Discovered Graphs to Analytic NOON State Families", arXiv:2611.10685.** Turns discovered instances into analytic families [T]. [arXiv](https://arxiv.org/pdf/2611.10685)
- **[2026-10] "esQueranto: Differentiable Structured Quantum Light for Automated Scientific Discovery", arXiv:2610.00750** [T]. [arXiv](https://arxiv.org/pdf/2610.00750)
- **[2026-09] Grinch, arXiv:2609.08073.** Reward-proportional sampling of optical circuit graphs, built edge by edge with no dataset. [arXiv](https://arxiv.org/html/2609.08073)
- **[2026-05] "Automated experimental design for high-probability entanglement generation", arXiv:2605.02721** [T]. [arXiv](https://arxiv.org/pdf/2605.02721)
- **[2026-02] Arlt, Duan, Li, Xie, Wu, Krenn, "Meta-designing quantum experiments with language models", Nature Machine Intelligence 8, 148–157 (2026; issued 2026-02-19), DOI 10.1038/s42256-025-01153-0, arXiv:2406.02470 (v1 2024-06, revised 2025-07).**
  - Agent: a language model trained to write human-readable Python code that solves a whole class of problems in a single pass ("meta-design"). For example, one program produces the setups for a family of quantum states at every size.
  - Scale: meta-solutions found for a selection of 20 classes of quantum states.
  - Result: "We discovered two construction rules that we did not know of before" (Arlt).
  - Earlier versions appeared at ICML 2024 and NeurIPS 2024 workshops.
  - Phys.org (Feb 2026) presents it as a way to make AI's unorthodox solutions interpretable.
  - Sources: [Tübingen record](https://bibliographie.uni-tuebingen.de/xmlui/handle/10900/181410?show=full); [arXiv](https://arxiv.org/abs/2406.02470); [Physics World](https://physicsworld.com/a/meta-design-language-models-generate-novel-quantum-experiments/); [ICML 2024](https://icml.cc/virtual/2024/36793); [Phys.org](https://phys.org/news/2026-02-ai-easily-solutions-unusual-quantum.html)
- **[2025-12] Anubuddhi, a multi-agent AI system for designing and simulating quantum-optics experiments, arXiv:2512.15736** [T]. [arXiv](https://arxiv.org/pdf/2512.15736)
- **[2025-11] LLM agents with access to "intelligent tools" for autonomous quantum-physics research, arXiv:2511.11752** [T]. [arXiv](https://arxiv.org/pdf/2511.11752)
- **[2025-07] Quanta Magazine feature.** Krenn's student Sören Arlt applied PyTheus to entanglement swapping. The resulting layout looked nothing like Zeilinger's 1993 design, and a separate mathematical check confirmed it would still entangle particles with no shared past. [Quanta](https://www.quantamagazine.org/ai-comes-up-with-bizarre-physics-experiments-but-they-work-20250721/)
- **[2024-03] "Virtual Reality for Understanding Artificial-Intelligence-driven Scientific Discovery with an Application in Quantum Optics", arXiv:2403.00834** [T]. [arXiv](https://arxiv.org/pdf/2403.00834)
- **[2022-10, Quantum 2023] Ruiz-Gonzalez, Arlt, Petermann, Sayyad, Jaouni, Karimi, Tischler, Gu, Krenn, "Digital Discovery of 100 diverse Quantum Experiments with PyTheus", Quantum 7, 1204 (2023), DOI 10.22331/q-2023-12-12-1204, arXiv:2210.09980.**
  - An open-source, graph-based, highly efficient design framework.
  - Tasks: generating entangled and mixed states, maximizing entanglement, quantum communication, and quantum computation.
  - Result: 100 previously unknown or improved experimental implementations, which are interpretable.
  - Installed with pip as `pytheusQ`.
  - Sources: [arXiv](https://arxiv.org/abs/2210.09980v2); [MPG](https://pure.mpg.de/pubman/item/item_3450733_5); [GitHub](https://github.com/mahmoudifar/PyTheus)
- **[E, 2021] Krenn, Kottmann, Tischler, Aspuru-Guzik, "Conceptual Understanding through Efficient Automated Design of Quantum Optical Experiments" (THESEUS), PRX 11, 031044 (2021), arXiv:2005.06443.**
  - Uses a graph representation of experiments and is orders of magnitude faster than earlier design algorithms.
  - The interpretable solutions gave a more intuitive combinatorial picture of heralded photonic CNOT gates.
  - Sources: [DOAJ/PRX](https://doaj.org/article/a34eeb74fc0845a697308b7fa2729cdc); [arXiv](https://arxiv.org/abs/2005.06443v3)
- **[E, 2018] Melnikov, Poulsen Nautrup, Krenn, Dunjko, Tiersch, Zeilinger, Briegel, "Active learning machine learns to create new quantum experiments", PNAS 115(6), 1221–1226 (2018), arXiv:1706.00868.**
  - Agent: a projective-simulation RL agent designing photonic experiments that produce high-dimensional entangled multiphoton states.
  - Result: produced more distinct entangled states than the best earlier automated approaches.
  - It rediscovered techniques that were only then becoming standard in quantum optics.
  - Sources: [arXiv](https://arxiv.org/abs/1706.00868v3); [DeepAI](https://deepai.org/publication/active-learning-machine-learns-to-create-new-quantum-experiments)

### Inferences
- The progression runs from searching for instances (MELVIN, THESEUS, PyTheus) to searching for programs that generate families (meta-design), and in 2026 to extracting analytic families from discovered graphs (the NOON paper). This parallels the circuit-side move to LLMs writing generator programs (Section 0).
- The verifier in this line is exact: the graph and perfect-matching formalism computes the output state. That is why optical design was one of the first areas where AI produced human-interpretable discoveries.
- GFlowNet-style generative design (Grinch) adds a way to sample many diverse designs, where PyTheus mostly returned single optimized solutions.

### Gaps
- MELVIN (Krenn, Malik, Fickler, Lapkiewicz, Zeilinger, PRL 116, 090405 (2016)) was not re-verified in this session; the citation is from memory.
- No experimental realizations of designs from 2025–26 were confirmed. The Quanta article reports a mathematical check of the PyTheus entanglement-swapping design, not a lab test.
- Results for Grinch, esQueranto, arXiv:2605.02721, arXiv:2611.10685 and Anubuddhi were not retrieved.
- "Automated discovery of high-probability heralded schemes for path-entangled states" (arXiv:2607.25501) was seen only in a reference list.

## 5. RL and learning for quantum control and many-body state preparation

### Takeaway
RL for quantum control is well established in simulation:
- Bukov et al. (PRX 2018) found a glassy phase in protocol space.
- Metz & Bukov (Nature Machine Intelligence 2023) used MPS-based agents to prepare Ising-chain ground states, including critical ones.

On hardware, there are real-time, model-free demonstrations for a single qubit (Reuer et al. 2023, FPGA, sub-microsecond latency) and model-free bosonic control (Sivak et al., PRX 2022). In 2025–26, work moves to measurement-feedback RL under partial observability (Bose–Hubbard and GHZ targets, numerics only) and to small hardware runs. I found no RL-driven many-body state preparation on Rydberg-array hardware.

### Cited Findings
- **[2026-08] "Reinforcement-Learned Electric-Field Sensing with Asymmetrically Blockaded Rydberg Arrays", arXiv:2608.01832** [T]. [arXiv](https://arxiv.org/pdf/2608.01832)
- **[2026-07] Isogawa et al., agentic NV-center experiments, arXiv:2607.25145** (see Section 0). [arXiv](https://arxiv.org/abs/2607.25145)
- **[2026-06] "Experiment-compatible measurement–feedback quantum state preparation with reinforcement learning", arXiv:2606.13005.**
  - Adaptive measurement–feedback RL under partial observability, where the controller sees only measurement outcomes available in experiment.
  - Targets: Bose–Hubbard ground states and GHZ states.
  - Demonstrations are numerical.
  - Sources: [arXiv](https://arxiv.org/abs/2606.13005); [HTML](https://arxiv.org/html/2606.13005)
- **[2026-06] "Learning ground state observables from quantum computing experiments", arXiv:2606.15983** (not RL; hardware context).
  - Notes that earlier experiments were mostly limited to small systems or highly structured states, because many-body ground states are hard to prepare on quantum processors.
  - Uses approximate ground states of a 2D Heisenberg XXZ model with up to 115 qubits.
  - Source: [arXiv](https://arxiv.org/abs/2606.15983)
- **[2026-05] "Reinforcement Learning Assisted Quantum Simulation of …", arXiv:2605.18569.** Extends an RL-driven contracted quantum eigensolver from ground states to excited states and real-time dynamics. [arXiv](https://arxiv.org/pdf/2605.18569)
- **[2026-05] "Intelligent Optimal Control of Rydberg Gates with Incremental-Update Deep Reinforcement Learning", arXiv:2605.04628** [T]. [arXiv](https://arxiv.org/pdf/2605.04628)
- **[2026-01] "Reinforcement Learning for Quantum Technology", arXiv:2601.18953** (review). Describes a physics-informed PPO agent, with priors on the time scales of the dynamics, that prepares states and gates in a multi-level Λ system, a Rydberg atom and a transmon (simulation). [arXiv](https://arxiv.org/html/2601.18953v1)
- **[2025-07] "Preparing spin-squeezed states in Rydberg atom arrays via quantum optimal control", arXiv:2507.07875.** A non-RL optimal-control baseline for Rydberg arrays. [arXiv](https://arxiv.org/html/2507.07875v3)
- **[2025-03] "Robust Quantum Control using Reinforcement Learning from Demonstration", arXiv:2503.21085.** Benchmarks include preparing an even cat state |α⟩+|−α⟩ with α=2. [arXiv](https://arxiv.org/pdf/2503.21085)
- **[2024-11] Measurement-feedback and gadget-based RL:**
  - "Learning Feedback Mechanisms for Measurement-Based Variational Quantum State Preparation", arXiv:2411.19914 [T]. [arXiv](https://arxiv.org/pdf/2411.19914)
  - Gadget-RL on real hardware with small systems, arXiv:2411.00230 (see Section 2). [arXiv](https://arxiv.org/html/2411.00230)
- **[2024-09] Neural-network-driven cat-state preparation in a cavity coupled to a qubit (experiment), arXiv:2409.05557.** [arXiv](https://arxiv.org/html/2409.05557v1)
- **[2024-01] "Machine learning-enhanced optical tweezers for defect-free rearrangement", arXiv:2401.04893** [T]. [arXiv](https://arxiv.org/pdf/2401.04893)
- **[2023-11] Reuer, Landgraf, Fösel et al. (ETH Zurich and partners), "Realizing a deep reinforcement learning agent for real-time quantum feedback", Nature Communications 14, 7138 (2023), DOI 10.1038/s41467-023-42901-3.**
  - A small neural network on an FPGA acts within a sub-microsecond window to prepare (reset) a transmon in its ground state.
  - Training is model-free and uses only measurement data; the network adds about 48 ns of latency.
  - Episodes are sent to a PC, which computes rewards and updates the policy.
  - Sources: [DOI](https://doi.org/10.1038/s41467-023-42901-3); [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10628214/)
- **[2023-07] Metz & Bukov, "Self-correcting quantum many-body control using reinforcement learning with tensor networks", Nature Machine Intelligence 5, 780–791 (2023), arXiv:2201.11790.**
  - Agent: matrix product states both represent the quantum state and form part of the trainable RL architecture.
  - Task: preparing ground states of the quantum Ising chain, including critical states.
  - Results: finds universal controls, adapts protocols when the dynamics are perturbed, and was tested under experimentally relevant noise. It controls more interacting particles than standard neural-network methods.
  - Sources: [arXiv](https://arxiv.org/abs/2201.11790); [EPFL](https://infoscience.epfl.ch/record/304148); [Nature PDF](https://www.nature.com/articles/s42256-023-00687-5.pdf)
- **[E, 2022] Sivak et al., "Model-Free Quantum Control with Reinforcement Learning", PRX 12, 011059 (2022), arXiv:2104.14539.**
  - Agent: circuit-based, learning by trial and error with measurement outcomes as its only information.
  - Task: state preparation in an oscillator coupled to an ancilla qubit, rewarded only through ancilla measurements.
  - Shown in simulation with both open- and closed-loop (adaptive feedback) control.
  - Needs no averaging for tomography or fidelity estimation, and is much more sample-efficient than widely used model-free methods.
  - Sources: [ADS](https://ui.adsabs.harvard.edu/abs/10.1103/PhysRevX.12.011059/abstract); [arXiv](https://www.arxiv.org/pdf/2104.14539)
- **[E, 2018] Bukov, Day, Sels, Weinberg, Polkovnikov, Mehta, "Reinforcement Learning in Different Phases of Quantum Control", PRX 8, 031086 (2018), arXiv:1705.00565.**
  - RL performs comparably to optimal control at finding short, high-fidelity driving protocols in non-integrable many-body qubit systems.
  - It finds nearly optimal protocols even in the glassy phase, where optimal manipulation is exponentially hard.
  - The space of protocols shows a spin-glass-like transition as a function of protocol duration.
  - Sources: [arXiv](https://arxiv.org/pdf/1705.00565); [DOAJ](https://doaj.org/article/0177872bc18644bc8a7b543ac628db38)

### Inferences
- There is a scale gap between simulated many-body RL (MPS agents, Bose–Hubbard and GHZ feedback) and hardware RL (one transmon reset, a few-qubit Ising model). Making RL work with measurement outcomes only (Sivak; arXiv:2606.13005) is the necessary bridge, and it has not yet been demonstrated at many-body scale.
- On hardware, the strongest many-body state-preparation results in my sources (115-qubit XXZ approximate ground states) come from non-RL methods. RL has not shown a hardware advantage in many-body preparation.
- LLM agents are entering control and calibration through the lab-automation route (NV-center agents, and reported autonomous processor bring-up) rather than through pulse-level RL.

### Gaps
- Not retrieved: the largest system size in Metz & Bukov, and whether Sivak et al. (PRX 2022) included an experimental demonstration.
- No RL-for-Rydberg-array many-body state preparation on hardware was found. Whether one exists is open.
- "Vibe Calibration" (arXiv:2606.22376) and the Patterns 2025 agentic-lab paper were not verified (see Section 0).
- One autonomous multi-agent lab for spins in diamond appeared in search summaries, possibly arXiv:2610.40589 ("Collaborative AI agents for autonomous experimental research…"). The mapping of description to ID is uncertain [U].

## 6. AI discovery of new quantum algorithms or protocols (not just compiled circuits), 2024–2026

### Takeaway
There are credible examples of algorithmic content that originated with AI, but they are narrow:
- formally verified resource improvements to known subroutines (290× for arithmetic, 1.4× for parallel Grover);
- new reversible-arithmetic architectures, such as modular inversion via divsteps;
- new construction rules for families of states in quantum optics (meta-design);
- unconventional but verified entanglement-swapping layouts (PyTheus);
- a rediscovered Karatsuba-like multiplication (AlphaTensor-Quantum).

I found no AI-discovered quantum algorithm with a new asymptotic speedup. Google's 2025 headline algorithms (DQI, Quantum Echoes) were designed by humans.

### Cited Findings
- **[2026-11] DeepMind + Google Quantum AI, Lean-verified programs.**
  - Arithmetic improved 290×, parallel Grover up to 1.4×, and simulation subroutines for fusion, photochemistry and drug metabolism 11×, 8× and 2.5×.
  - All are improvements on known tasks, with machine-checked proofs.
  - Source: [arXiv](https://arxiv.org/html/2611.09822v1)
- **[2026-10] Gu, Kemeny, Yelin, Lukin.**
  - Found a reversible modular-inversion architecture that makes divsteps, a classical algorithm previously impractical for reversible integer arithmetic, viable.
  - On 7 of 8 tasks, cheaper than the best published designs.
  - Source: [arXiv](https://arxiv.org/html/2610.22367v1); [DAIR](https://academy.dair.ai/papers/discovering-state-of-the-art-quantum-circuits-with-large-language-models-2610.22367)
- **[2026-06] qBraid/AlphaEvolve.** A fermionic encoding with distance 5 where human designs had reached 3. This is a new encoding, not an algorithm, and QEC-adjacent. [qBraid](https://www.qbraid.com/blog-posts/qbraid-uses-alphaevolve-for-quantum-error-correction)
- **[2026-02] Meta-design.** Two previously unknown construction rules for families of quantum states. [Physics World](https://physicsworld.com/a/meta-design-language-models-generate-novel-quantum-experiments/)
- **[2025-07] PyTheus entanglement swapping.** An unrecognizable but verified alternative design for entanglement swapping. [Quanta](https://www.quantamagazine.org/ai-comes-up-with-bizarre-physics-experiments-but-they-work-20250721/)
- **[2025-03] AlphaTensor-Quantum.** Rediscovered a Karatsuba-like multiplication for finite fields, and matched the best human designs for Shor arithmetic and chemistry primitives. [DeepMind](https://deepmind.google/research/publications/alphatensor-for-optimizing-quantum-computations/)
- **[2025-10] "Quantum Agents for Algorithmic Discovery", arXiv:2510.08159.** A proposal paper; I saw no concrete discovered algorithm in the retrieved text. [arXiv](https://arxiv.org/pdf/2510.08159)
- **Human-designed counterexamples:**
  - Decoded Quantum Interferometry (Google Quantum AI and collaborators, Nature 2025) is a human-designed optimization algorithm. [Google Research](https://research.google/blog/a-new-quantum-toolkit-for-optimization/)
  - Quantum Echoes (an OTOC algorithm, Oct 2025) was presented as verifiable quantum advantage on Willow. [Google blog](https://blog.google/innovation-and-ai/technology/research/quantum-echoes-willow-verifiable-quantum-advantage/)

### Inferences
- Discovery is credible exactly where an exact verifier exists: Lean proofs, reversible simulation, optical graph semantics, exact code distance. The output is usually a better construction for a known task, not a new algorithmic paradigm.
- An emerging pattern is "AI search, then human or analytic generalization" (meta-design rules; NOON families from discovered graphs), which may be how AI-found constructions become citable algorithms.

### Gaps
- No independent replication was found for any 2026 claim.
- I found no AI-discovered protocol with a proven new complexity separation or asymptotic speedup. This is my conclusion from not finding one, not a sourced fact.

## 7. Fermionic circuits (one example among others): ML or search for fermionic Gaussian/matchgate compilation, Givens networks, swap networks, fermion-to-qubit mappings; plus entangling and pseudorandom circuits

### Takeaway
Fermionic compilation is one of the few NISQ-relevant areas with exact, efficient classical verification: matchgate circuits correspond to fermionic Gaussian unitaries, which correspond to O(2n). In 2025–26 the first learning and agentic methods appeared:
- offline RL over free-fermionic subroutines (F2);
- agentic compilation of Givens-rotation look-up tables (Google Quantum AI);
- AI-agent optimization of fermion-to-qubit encodings and Hamiltonian layout, with a frozen exact verifier;
- AlphaEvolve-evolved GSE encodings.

They sit alongside strong non-ML baselines: Fermihedral (SAT), Treespilation (annealing), HATT, and Majorana swap networks.

On entangling and pseudorandom circuits, AI work targets entanglement as a reward (LLM test-time learning, optical design). I found no ML work on discovering pseudorandom unitaries or t-designs.

### Cited Findings
- **[2026-10] AI-agent optimization with a frozen exact verifier, applied to fermion-to-qubit encoding and Hamiltonian layout, arXiv:2610.31352.** [arXiv](https://arxiv.org/html/2610.31352v1)
- **[2026-10] Khattar et al.** Agentic compilation targets include quantum look-up tables with Givens rotations for chemistry simulation, within the classically verifiable kickmix class. [arXiv](https://arxiv.org/abs/2610.23832)
- **[2026-07] COSMA, "Communication-aware Optimization of Fermionic Simulation Kernels for Modular Quantum Architectures", arXiv:2607.09381.** Describes Treespilation as a simulated-annealing search over PPTT mappings minimizing a user-defined cost. [arXiv](https://arxiv.org/pdf/2607.09381)
- **[2026-06] qBraid/AlphaEvolve evolved GSE fermionic encodings, arXiv:2606.25870** (see Section 0). [arXiv](https://arxiv.org/pdf/2606.25870)
- **[2026-05] "Fermion lattices can be simulated by same-size qubit lattices with O(1) interaction overhead", arXiv:2605.12600.** Notes that fermionic swap networks are a natural way to implement fermion-routing permutations. [arXiv](https://arxiv.org/pdf/2605.12600)
- **[2026-05] CO-MAP, RL for qubit allocation, arXiv:2605.13638** (adjacent: qubit mapping and SWAP reduction, not fermionic encodings). [arXiv](https://arxiv.org/html/2605.13638v1)
- **[2026-03] "Matchgate circuit representation of fermionic Gaussian states: optimal preparation, approximation, and classical simulation", arXiv:2603.05675** (non-ML baseline). Decomposes a state's covariance matrix directly into Givens rotations without diagonalizing it first, giving a shallow, depth-O(β) construction in the paper's notation. [arXiv](https://arxiv.org/html/2603.05675v1)
- **[2025-12] Decker, Watson and six co-authors, "F2: Offline Reinforcement Learning for Hamiltonian Simulation via Free-Fermionic Subroutine Compilation", arXiv:2512.08023.**
  - Environment: built on classically simulable free-fermionic subroutines, modeled as an MDP with discrete actions and deterministic dynamics.
  - Method: architecture and objective choices to stabilize long-horizon value learning, plus synthetic offline training trajectories.
  - Target: Trotter-based Hamiltonian-simulation circuits.
  - The authors argue that hand-engineered heuristics cannot learn input-dependent structure.
  - Benchmark numbers: n/r.
  - Sources: [arXiv](https://arxiv.org/abs/2512.08023); [HTML](https://arxiv.org/html/2512.08023)
- **[2025-09] "Improving fermionic variational quantum eigensolvers with Majorana swap networks", arXiv:2509.07855.**
  - Cyclic compilation localizes all of the up to O(M⁴) two-particle terms using O(M³) Majorana-swap transpositions.
  - About 50% less depth and 20% fewer two-qubit gates under all-to-all connectivity.
  - Source: [arXiv](https://arxiv.org/abs/2509.07855)
- **[2025] Fermionic learning theory and encodings [T]:**
  - "Efficiently learning fermionic unitaries with few non-Gaussian gates", arXiv:2504.15356. [arXiv](https://arxiv.org/pdf/2504.15356)
  - "Fermion-to-qubit encodings with arbitrary code distance", arXiv:2505.02916. [arXiv](https://arxiv.org/pdf/2505.02916)
- **[2025] "Unified framework for matchgate classical shadows", npj Quantum Information.** Random matchgate ensembles have an average gate count of n(n−1)/2 and a maximal depth of 2n layers of commuting Givens rotations. [npj QI](https://www.nature.com/articles/s41534-025-01015-y)
- **[2024-09] "Ternary Tree Fermion-to-Qubit Mapping with Hamiltonian Aware Optimization" (HATT), arXiv:2409.02010.** A polynomial-time, bottom-up ternary-tree construction that is aware of the Hamiltonian and reduces Pauli weight. It argues that Fermihedral's SAT approach has scalability problems. [arXiv](https://arxiv.org/html/2409.02010v1)
- **[2024-03] Liu, Che, Zhou, Shi, Li, "Fermihedral: On the Optimal Compilation for Fermion-to-Qubit Encoding", ASPLOS '24, DOI 10.1145/3620666.3651371, arXiv:2403.17794.**
  - Method: finding the optimal encoding is posed as a SAT problem using Pauli algebra.
  - Scaling: for larger cases it ignores algebraic independence and uses simulated annealing on Hamiltonian-independent encodings.
  - Result: beats asymptotically optimal encodings on Pauli weight, with better accuracy in noisy simulation and on an IonQ device.
  - A third-party "up to 60%" gate reduction figure is unverified.
  - Sources: [arXiv](https://arxiv.org/html/2403.17794v2); [GitHub](https://www.Github.com/acasta-yhliu/fermihedral); [Emergent Mind](https://www.emergentmind.com/papers/2403.17794)
- **[2024-03] Miller, Glos, Zimborás, "Treespilation: Architecture- and State-Optimised Fermion-to-Qubit Mappings", arXiv:2403.03992.** Cuts CNOT counts for ADAPT-VQE ground states by up to 74% with full connectivity, with similar reductions on IBM Eagle and Google Sycamore. [arXiv](https://arxiv.org/abs/2403.03992)
- **[2023-12] "Learning, Optimizing, and Simulating Fermions with Quantum Computers", arXiv:2312.10399.** A QR-style decomposition into nearest-neighbor Givens rotations: O(n²) rotations at depth O(n) (Reck-style), including rotations that do not conserve particle number. [arXiv](https://arxiv.org/pdf/2312.10399)
- **[E, 2022] Wan et al., "Matchgate Shadows for Fermionic Quantum Simulation", arXiv:2207.13723** (Comm. Math. Phys. 2023).
  - Under Jordan–Wigner, matchgate circuits correspond to fermionic Gaussian unitaries and, up to global phase, to O(2n).
  - Nearest-neighbor XX rotations plus Z rotations give Givens rotations that generate SO(2n).
  - Sources: [arXiv](https://arxiv.org/pdf/2207.13723); [CMP](https://link.springer.com/article/10.1007/s00220-023-04844-0)
- **[E, 2021] Hashim et al., "Optimized fermionic SWAP networks with equivalent circuit averaging for QAOA", arXiv:2111.04572.**
  - An n-qubit FSWAP network routes in n−1 steps with O(n) depth.
  - The FSWAP gate compiles to 3 native Sycamore gates.
  - Equivalent circuit averaging cut total variation distance by 60% on average for depth-1 QAOA on 4 transmons.
  - Source: [arXiv](https://arxiv.org/pdf/2111.04572)
- **[E, 2020] Jiang et al., "Optimal fermion-to-qubit mapping via ternary trees", Quantum 4, 276 (2020).** [Quantum](https://quantum-journal.org/papers/q-2020-06-04-276/)
- **Search result on RL for mappings:** a dedicated search found no paper applying RL directly to choosing fermion-to-qubit mappings. The nearest were F2, SAT, annealing and heuristic methods, and RL for generic qubit mapping. [F2](https://arxiv.org/abs/2512.08023); [CO-MAP](https://arxiv.org/html/2605.13638v1)
- **Ground states with learned or generated circuits:**
  - GPT-QE beats CCSD in N2 dissociation (Section 3). [arXiv](https://arxiv.org/html/2401.09253v2)
  - Treespilation lowers the cost of ADAPT-VQE ground states. [arXiv](https://arxiv.org/abs/2403.03992)
  - Metz & Bukov prepare Ising-chain (Jordan–Wigner free-fermion) critical ground states with MPS-RL. [arXiv](https://arxiv.org/abs/2201.11790)
- **Entangling circuits as targets:**
  - LLM test-time learning driven by entanglement feedback. [arXiv:2602.03466](https://arxiv.org/html/2602.03466v1)
  - Projective-simulation and PyTheus designs for high-dimensional and maximal entanglement. [PNAS preprint](https://arxiv.org/abs/1706.00868v3); [PyTheus](https://arxiv.org/abs/2210.09980v2)
- **Random-circuit data, and the bias it introduces:**
  - KetGPT argues that random circuits used as benchmarks do not reflect algorithm structure. [arXiv](https://arxiv.org/abs/2402.13352)
  - Fösel et al. trained on random 12-qubit circuits. [arXiv](https://arxiv.org/pdf/2103.07585)
  - The AlphaTensor-Quantum reusability study targets random 5–8-qubit circuits. [arXiv](https://arxiv.org/abs/2511.09951)

### Inferences
- Fermionic Gaussian compilation fits verifier-first AI search well: the target group is exactly known (SO(2n)), fidelity can be checked in polynomial time, and strong classical baselines exist (Givens QR depth O(n), SAT-optimal encodings, Treespilation). This explains why F2 and the 2026 agentic papers chose it.
- Co-optimizing fermion-to-qubit mapping, routing and Givens/swap-network layout is currently done by separate tools (Fermihedral or HATT for mappings, Treespilation for layout and mapping, FSWAP or Majorana networks for routing). Agentic and RL methods are only starting to combine them (arXiv:2610.31352).
- Pseudorandom circuits appear in this literature only as training or benchmark data, not as discovery targets.

### Gaps
- F2's quantitative gains and baselines were not retrieved.
- No ML or agentic work was found on discovering shallow pseudorandom unitaries, unitary t-designs, or magic-generating ensembles. This absence is a search finding with limited depth, because the search budget ran out before a dedicated query.

## Gaps and opportunities

### Takeaway
The main gaps in the literature:
- equal-budget benchmarking against strong classical baselines;
- exact verifiers for NISQ-relevant circuit classes beyond reversible and Clifford logic;
- generalization to larger sizes;
- learned many-body state preparation at hardware scale;
- AI search aimed at entanglement, magic or pseudorandomness;
- independent replication, with compute accounting, of the 2026 claims about LLM-driven discovery.

### Cited Findings
- Verification is the stated bottleneck for agentic compilation: classical verification of arbitrary quantum circuits is intractable, so agentic results are trusted only for verifiable classes. [Google Open Source blog](https://opensource.googleblog.com/2026/10/announcing-verifiable-quantum-circuit-compilation-competition-on-kaggle.html)
- AlphaTensor-Quantum does not yet handle T-depth. [The Quantum Insider](https://thequantuminsider.com/2025/03/21/quantum-optimization-gets-an-ai-boost-with-alphatensor-quantum/)
- Exact learned synthesis is still small:
  - Clifford+T exact synthesis up to 5 qubits. [arXiv:2404.14865](https://arxiv.org/pdf/2404.14865)
  - Clifford synthesis up to 11 qubits and linear functions up to 9. [arXiv:2405.13196](https://arxiv.org/abs/2405.13196)
  - CNOT minimization up to 8 qubits. [arXiv:2604.13812](https://arxiv.org/abs/2604.13812)
  - Non-Clifford-count RL on two qubits. [arXiv:2509.21709](https://arxiv.org/abs/2509.21709v1)
- LLM-designed variational circuits on noisy hardware remain an open challenge. [arXiv:2606.13380](https://arxiv.org/html/2606.13380)
- Noise-aware RL-QAS needs better simulation efficiency for larger problems. [Liner review of CRLQAS](https://liner.com/review/curriculum-reinforcement-learning-for-quantum-architecture-search-under-hardware-errors)
- Data scarcity is the main obstacle for LLM-based algorithm design. [QCircuitNet](https://arxiv.org/html/2410.07961v1)
- Random circuits used as data do not reflect algorithm structure. [KetGPT](https://arxiv.org/abs/2402.13352)
- A single scalar verifier reward can destabilize RL post-training. [Vista](https://www.alphaxiv.org/abs/2605.vista)
- The base OpenEvolve prompt carries too little quantum information. [arXiv:2609.38883](https://arxiv.org/html/2609.38883)
- Equal-budget comparisons with random, evolutionary and Bayesian search are rare. One explicit attempt is [arXiv:2607.19506](https://arxiv.org/pdf/2607.19506).
- Many-body ground-state preparation on hardware has mostly been limited to small or highly structured states. [arXiv:2606.15983](https://arxiv.org/abs/2606.15983)
- RL on hardware has been shown only at small scale. [arXiv:2411.00230](https://arxiv.org/html/2411.00230); [Reuer et al.](https://doi.org/10.1038/s41467-023-42901-3)
- Vendor or secondary-only figures:
  - AlphaEvolve's 10× lower error on Willow. [DeepMind](https://deepmind.google/blog/alphaevolve-impact/)
  - qBraid's results are company-reported, with circuit-level fault tolerance untested. [qBraid](https://www.qbraid.com/blog-posts/qbraid-uses-alphaevolve-for-quantum-error-correction)
  - The Gu et al. gain is reported as 39× in one place and "seventeenfold" in another. [DAIR](https://academy.dair.ai/papers/discovering-state-of-the-art-quantum-circuits-with-large-language-models-2610.22367); [Substack](https://brianroemmele.substack.com/p/the-gate-that-costs-forever-a-forgotten)
- Two examples of turning discovered instances into analytic families:
  - meta-design construction rules. [Physics World](https://physicsworld.com/a/meta-design-language-models-generate-novel-quantum-experiments/)
  - analytic NOON families. [arXiv:2611.10685](https://arxiv.org/pdf/2611.10685)
- Dynamic circuits are only beginning to enter learned synthesis. [arXiv:2508.21217](https://arxiv.org/pdf/2508.21217)

### Inferences
The following are concrete open problems, inferred from the findings above:

1. **Equal-compute benchmarking.** Build a public suite that runs LLM-evolution (OpenEvolve/AlphaEvolve), AlphaZero/MCTS, model-free RL, simulated annealing, genetic programming, SAT/ILP-optimal synthesis and rewrite-based superoptimizers under the same budget of evaluations and wall-clock time. Report Pareto fronts (gate count, depth, qubits, T-count/T-depth) and the cost in LLM tokens or GPU hours. The Kaggle kickmix competition and QAS benchmarks are partial starts.
2. **Verifiers for classes beyond reversible and Clifford circuits.** Agentic search needs exact, scalable checkers for classes that matter on NISQ devices. Candidates:
   - matchgate/fermionic Gaussian circuits (SO(2n) equality in polynomial time);
   - Clifford plus a few T gates (stabilizer-rank methods);
   - Gaussian bosonic optics;
   - circuits with bounded bond dimension (MPS);
   - Gaussian plus a few non-Gaussian gates, building on learning theory such as arXiv:2504.15356.

   Each would open a certified search space between "classically trivial" and "unverifiable".
3. **Programs that generalize across sizes.** Learned synthesis tops out at about 5–11 qubits. Meta-design shows that LLM-written generator programs can generalize to all N. Open problem: LLM or evolutionary search over parametric circuit families (arithmetic, Givens networks, state-preparation families), with proofs of correctness for all N, for example in Lean as in arXiv:2611.09822.
4. **Joint, hardware-aware objectives.** Optimize T-count, T-depth, CNOTs under connectivity, and qubit count together, including dynamic circuits (mid-circuit measurement and feedforward). No retrieved work optimizes all of these jointly. AlphaTensor-Quantum lacks T-depth, and AlphaCNOT handles connectivity only up to 8 qubits.
5. **Learned many-body state preparation at hardware scale.** Demonstrate RL that uses only measurement outcomes (as Sivak et al. and arXiv:2606.13005 do) for more than about 20 sites on hardware (Rydberg, superconducting, trapped ions). Compare it against adiabatic and optimal-control baselines (for example arXiv:2507.07875) and against non-RL preparation (arXiv:2606.15983).
6. **AI search aimed at entanglement, magic and pseudorandomness.** No retrieved work uses AI to discover shallow pseudorandom unitaries, t-designs, scrambling circuits, or circuits that maximize magic. Reward candidates include frame potentials, anticoncentration, OTOCs and stabilizer Rényi entropies, all computable for structured or classically simulable classes. T-count work (AlphaTensor-Quantum, Kremer 2025) covers only magic minimization.
7. **Unified search over fermionic simulation stacks.** Learn or evolve the fermion-to-qubit mapping, qubit layout and Givens/swap-network schedule together, with exact Gaussian verification. Existing tools handle one layer each (Fermihedral or HATT, Treespilation, FSWAP or Majorana networks), and RL for mapping selection is missing.
8. **Reward design for verifier-trained LLMs.** Study hierarchical and multi-signal rewards (QUASAR, Vista) and quantum-aware prompting (arXiv:2609.38883) systematically. Measure whether validity gains translate into better optimality or novelty.
9. **Replication and contamination audits.** Re-run 2026 LLM-discovery claims (Gu et al.; DeepMind/Google Quantum AI; qBraid; the AlphaEvolve Willow claim) independently. Check whether published optimal designs appear in training data, and report model identity and compute.
10. **From instances to theory.** Extend the "discovered graphs to analytic families" approach (meta-design, NOON families) to circuits found by RL or LLMs (AlphaTensor-Quantum decompositions, ZX-RL rewrites, divsteps arithmetic). Extracting closed-form families would turn search output into reusable constructions.
11. **Datasets with verified labels.** Build large circuit corpora that carry verified properties (unitary or state, entanglement spectrum, magic, fermionic structure, hardware cost), going beyond QCircuitBench and KetGPT, to train and evaluate generative models (diffusion, GFlowNet, LLM) without the bias of random-circuit data.
12. **Generative models with hard constraints.** Diffusion and GFlowNet circuit generators do not guarantee correctness. Combining them with projection onto group structure (for example SO(2n) for matchgates), or with verifier-guided sampling, is underexplored.

### Gaps
- Research-process limits:
  - arxiv.org could not be fetched (DNS-blocked), so arXiv content comes only from search snippets.
  - The shared web-search budget (200 calls per turn) ran out before final verification of several items: Ostaszewski et al. 2021, MELVIN 2016, Fürrutter et al. NMI 2024 details, the authors of arXiv:2610.31352, Qubrio and arXiv:2508.21217, the "Fast and Slow" (GUOQ) baseline, and searches specific to magic and pseudorandom circuits.
  - Further searching would need a follow-up session.
- Many 2026 entries are confirmed only by title or abstract [T]. They should be checked against the PDFs before quantitative claims are reused.
