# Newest works (Jan 2025 – Nov 2026) on agentic AI discovery of quantum circuits and quantum error-correcting codes

> **STATUS: SWEEP NOT PERFORMED. ZERO ITEMS VERIFIED. Read this before using the notes.**
> - **Search budget.** The shared WebSearch budget for this turn (200 calls, shared by every agent in the turn) was already used up when this reassigned task began. None of the queries ran, including "LLM agent quantum error correction code discovery 2026", "AlphaEvolve quantum circuit discovery", "LLM-guided evolutionary search quantum error-correcting codes qLDPC" and "agentic AI quantum error correction arXiv 2025". Each one was refused.
> - **Network.** Direct verification with WebFetch also failed. **arxiv.org** and **deepmind.google** do not resolve in this cloud session (`getaddrinfo ENOTFOUND`: the hosts are outside the session's allowed network domains). As the tool guidance says, I stopped trying other hosts after the second failure.
> - **What the leads are.** Nothing below is confirmed. The items listed under "Gaps" are leads I recalled from model memory (knowledge cutoff mid-2026), and each is marked **UNVERIFIED**. I give an arXiv ID only where I recalled it with fair confidence, and it may still be wrong. **The report writer should not cite any of these items as verified facts.**
> - **What finishing the task needs.** (1) A follow-up turn, which brings a fresh search budget, or a higher `CLAUDE_CODE_MAX_WEB_SEARCHES_PER_SESSION`. (2) arxiv.org allowed in the session's network settings, and also deepmind.google, research.google, nature.com, journals.aps.org and quantum-journal.org as needed.
> - **Scope change.** The previous topic (monitored many-body dynamics, free fermions, disorder) was dropped at the coordinator's request. No notes file for it had been saved (checked: `research_notes/` contained no files), so nothing needed deleting.

## 1. 2026 items (Nov 2026 back to Jan 2026)

### Takeaway
I did not sweep this period. I could find and verify no 2026 item, and I hold no specific 2026 item in memory that I can responsibly list. This is the period most likely to contain the "new works" the user mentions, so it should be the first priority of the follow-up.

### Cited Findings
- None. No search or fetch succeeded.

### Inferences
- None I can ground.

### Gaps
- **Queries still to run** (newest month first, extended mode):
  - "LLM agent quantum error correction code discovery 2026"
  - "AlphaEvolve quantum circuit"
  - "AlphaEvolve quantum error correction"
  - "agentic quantum error correction"
  - "AI scientist quantum code"
  - "LLM-guided search stabilizer codes"
  - "evolutionary LLM qLDPC"
  - "OpenEvolve quantum"
  - "ShinkaEvolve quantum"
  - "FunSearch quantum code"
  - "Claude Code quantum error correction"
  - "Codex agent quantum circuit"
  - "Gemini agent quantum code discovery"
  - "LLM agent syndrome extraction circuit"
  - "LLM magic state distillation protocol discovery"
  - "Monte Carlo tree search quantum code discovery 2026"
  - "reinforcement learning fault-tolerant gadget discovery 2026"
  - "LLM agent decoder discovery quantum"
- **2026 venues still to check:** ICLR 2026, ICML 2026, QIP 2026, TQC 2026, IEEE QCE 2026, and NeurIPS 2026 accepted papers if they were listed by Nov 2026.

## 2. 2025 items (Dec 2025 back to Jan 2025)

### Takeaway
I did not sweep this period. The leads I recall (all UNVERIFIED) suggest that 2025 work on AI-driven code and circuit discovery was mostly reinforcement-learning (RL) agents optimizing codes or circuits. LLM-guided evolutionary systems (AlphaEvolve, ShinkaEvolve and similar) were shown mainly on mathematics and algorithms. Whether any AlphaEvolve-, OpenEvolve- or ShinkaEvolve-style system was applied to QEC codes or quantum circuits in 2025 is an open item.

### Cited Findings
- None verified.

### Inferences
- I have no grounded inference. A tentative impression from unverified recall: before 2026, "agentic" discovery in QEC mostly meant RL/MCTS agents with exact verifiers, not LLM agents.

### Gaps: unverified leads, newest first
For each lead, the follow-up should confirm the ID and date, the claims, whether the method actually counts as agentic, and whether the results are new.

- **Nov 2025: Bogdan Georgiev, Javier Gómez-Serrano, Terence Tao, Adam Zsolt Wagner, "Mathematical exploration and discovery at scale".** arXiv ID recalled as 2511.02864. UNVERIFIED.
  - What it did: applied AlphaEvolve to a large set of mathematics problems. I recall no quantum-code or quantum-circuit problem in it.
  - Relevance: evidence about the method, namely that LLM-guided evolution with automated evaluators finds new constructions and also rediscovers known ones.
  - New or rediscovered: mixed, depending on the problem.
- **Nov 2025: OpenAI with academic co-authors, "Early science acceleration experiments with GPT-5".** UNVERIFIED; no ID given.
  - A collection of case studies across several sciences.
  - I do not know whether it includes any quantum-circuit or QEC case.
- **Sep 2025: Scott Aaronson and Freek Witteveen, "Limits to black-box amplification in QMA".** UNVERIFIED.
  - Publicly reported to have used GPT-5 to suggest a key technical step of the proof.
  - This is quantum complexity theory, not codes or circuits. I list it only as a boundary case of LLM-assisted quantum theory.
- **Sep 2025: Robert Tjarko Lange, Yuki Imajuku, Edoardo Cetin (Sakana AI), "ShinkaEvolve".** UNVERIFIED.
  - An open-source, sample-efficient system for LLM-driven program evolution.
  - I recall no quantum application, so it serves only as a method baseline.
- **Jun 2025 (blog May 2025): Novikov et al., Google DeepMind, "AlphaEvolve: A coding agent for scientific and algorithmic discovery".** arXiv ID recalled as 2506.13131. UNVERIFIED.
  - Architecture as I recall it: a set of Gemini models proposes code modifications, an evolutionary program database chooses which candidates to build on, and automated evaluators written by the user score each candidate.
  - I do not recall any quantum-circuit or QEC application in the white paper. Later Google DeepMind and Google Quantum AI posts should be checked for quantum uses.
- **Apr 2025: Linus Jern, Valter Uotila, Cong Yu, Bo Zhao (authors as recalled), "Agent-Q: Fine-tuning large language models for quantum circuit generation and optimization".** UNVERIFIED.
  - As recalled: a fine-tuned open LLM that writes parameterized, QAOA-type circuits in OpenQASM for optimization problems, checked by simulating the circuits.
  - It is a single fine-tuned model, not an agent loop. The target is circuit design, not QEC.
- **Mar 2025: IBM Quantum (Dubal, Kremer, Martiel, Villar, Wang, Cruz-Benito; author list as recalled), "Pauli network circuit synthesis with reinforcement learning".** UNVERIFIED.
  - As recalled: an RL agent that synthesizes Pauli-network circuits, continuing IBM's line of RL-based "AI transpiler" passes.
  - Verifier: exact equivalence checking of the synthesized circuit.
  - Method and target: agent-based RL for synthesis.
- **Feb 2025: Austin Yubo He and Zi-Wen Liu, "Discovering highly efficient low-weight quantum error-correcting codes with reinforcement learning".** UNVERIFIED.
  - As recalled: an RL agent modifies the structure of qLDPC codes to lower stabilizer (check) weight while keeping the code distance. It reports large savings in physical-qubit overhead compared with earlier weight-reduction methods. I recall no reliable figure.
  - Verifier: code parameters [[n, k, d]] and check weights, computed exactly or bounded.
  - Presented as new codes, not rediscoveries; this needs checking.
- **2025 journal version of a Feb 2024 preprint: Ruiz, Laakkonen, Bausch, Balog, …, van de Wetering, Fawzi, Meichanetzidis, Kohli, "Quantum circuit optimization with AlphaTensor" (AlphaTensor-Quantum).** Preprint ID recalled as 2402.14396; venue recalled as Nature Machine Intelligence (2025). UNVERIFIED.
  - As recalled: a deep-RL agent with tree search, in the style of AlphaZero/AlphaTensor, searches for low-rank decompositions of a circuit's signature tensor to cut the T-count, including the use of gadgets.
  - Reported results: T-counts matching or beating the best known on benchmark arithmetic circuits, and a Karatsuba-like construction for multiplication in finite fields (a known classical idea rediscovered in a new setting).
  - Verifier: exact, because a valid decomposition gives an equivalent circuit.
  - Only partly in the time window: it is a 2024 preprint.

## 3. Lab and company announcements (Google DeepMind, Google Quantum AI, IBM, Quantinuum, Microsoft, NVIDIA, OpenAI, Anthropic, AWS)

### Takeaway
I did not check any of these, because no blog or press page could be reached.

### Cited Findings
- None.

### Inferences
- None I can ground.

### Gaps
- For each lab, the blog and newsroom should be searched for 2025–2026 AI-discovered circuits, codes, encoders, syndrome-extraction schedules, fault-tolerant gadgets, magic-state protocols or decoders. Open questions:
  - Did Google DeepMind or Google Quantum AI apply AlphaEvolve to quantum circuits or codes?
  - Did IBM's RL-based transpiler passes gain LLM or agentic components?
  - Did Quantinuum, Microsoft, NVIDIA, AWS, OpenAI or Anthropic announce codes, gadgets or protocols discovered by an agent?
- I recall no specific Anthropic or AWS announcement of this kind. Not recalling one is not evidence that none exists.
- **Outside the method scope:** AlphaQubit (Google DeepMind and Google Quantum AI, Nature 2024, recalled) is a transformer decoder trained by supervised learning, not an agentic method. The follow-up should check whether any 2025–2026 follow-ups (UNVERIFIED) add RL or agentic components.

## 4. 2025–2026 conference papers (NeurIPS 2025, ICLR 2026, ICML 2026, QIP 2026, TQC 2026, IEEE QCE)

### Takeaway
I did not check these.

### Cited Findings
- None.

### Inferences
- None.

### Gaps
- **Venues to search:** NeurIPS 2025 (main track, Datasets & Benchmarks, and the AI4Science and ML4PS workshops), ICLR 2026, ICML 2026, QIP 2026, TQC 2026, IEEE QCE 2025 and 2026.
- **Terms to search accepted-paper lists for:** "quantum error correction", "stabilizer code", "qLDPC", "syndrome extraction", "magic state", "quantum circuit synthesis", "LLM agent", "evolutionary", "MCTS".
- **Benchmarks to check:** LLM-for-quantum benchmarks such as QCircuitNet and Qiskit HumanEval (recalled, UNVERIFIED), in case they later added agentic or discovery tracks.

## 5. Pre-2025 baselines for judging "new or rediscovered" (UNVERIFIED in this session)

### Takeaway
Each method family has pre-2025 baselines, and any 2025–2026 claim should be compared with them before it is called new.

### Cited Findings
- None verified.

### Inferences
- A 2025–2026 result from an RL agent should be compared at least with the Marquardt-group RL work and AlphaTensor-Quantum. A result from an LLM-evolutionary system should be compared with FunSearch and AlphaEvolve, and with non-LLM search (SAT/ILP, genetic algorithms, simulated annealing) at equal compute.

### Gaps: unverified recall
- **FunSearch** (Romera-Paredes et al., Nature, online Dec 2023): LLM-guided program search with automated evaluation, applied to cap sets and bin packing. It is the template for later LLM-evolutionary systems.
- **Olle, Zen, Puviani, Marquardt** (npj Quantum Information, 2024, recalled): an RL agent that discovers codes and their encoding circuits together under a noise model.
- **Zen, Olle, Colmenarez, Puviani, Müller, Marquardt** (preprint 2024, recalled; arXiv ID recalled as 2402.17761): RL discovery of fault-tolerant circuits for preparing logical states.
- **Kremer et al., IBM** (preprint 2024, recalled): RL for synthesizing Clifford, permutation and linear-function circuits, and for routing.
- **Cao et al., "k-agents"** (Dec 2024, recalled): LLM agents that automate calibration of superconducting-qubit experiments. This is not code or circuit discovery.

## 6. Gaps and opportunities (general)

### Takeaway
The most reliable route for agentic discovery in QEC and circuits is to choose targets that have exact, cheap verifiers and a built-in check for novelty. The weakest points of the field, judged without the sweep, are reward hacking, telling rediscovery apart from new results, and comparing against strong search methods that use no LLM.

### Cited Findings
- None. This section is opinion and was not checked against sources in this session.

### Inferences (opinion)
- **Choose problems by verifier first.** Prefer targets with exact, fast checks:
  - parameters of small stabilizer or subsystem codes, with exact distance found by exhaustive search, ILP or SAT;
  - circuit-level fault distance of syndrome-extraction schedules and fault-tolerant gadgets, found by enumerating faults in a stabilizer-circuit simulator;
  - exact equivalence for Clifford and Pauli-network synthesis, by comparing tableaux;
  - exact output-error polynomials and acceptance rates for small magic-state distillation or cultivation-style protocols.
- **Evaluate in two stages to resist reward hacking.** Search with a cheap proxy, such as code-capacity distance. Before claiming anything, run a sealed, held-out evaluation at circuit level with different noise models and decoders. For decoders, use held-out syndrome data and a latency budget.
- **Build in a novelty check.** Compare automatically against code tables and the Error Correction Zoo, and put results into a canonical form under qubit permutations and local Clifford equivalence. That separates rediscovered results from new ones; report both kinds.
- **Create open benchmarks with compute accounting.** Fix the evaluators and the budgets (LLM tokens, simulator time), and include strong baselines that use no LLM (SAT/ILP, genetic algorithms, simulated annealing, plain MCTS/RL). That would show whether LLM guidance actually adds value.
- **Search for families, not instances.** Agents that output programs or constructions generating families of codes or circuits produce results that scale to larger n and can be proved. Single instances do not.
- **Co-design for hardware.** Search codes, schedules, qubit layouts and native-gate constraints together (connectivity, measurement and reset latency, leakage). These are the areas where hand design is weakest.
- **Make results reproducible and attributable.** Release prompts, seeds, evaluator code and full search logs, and state clearly what the humans and the AI each contributed.

### Gaps
- Whether LLM guidance beats strong search without an LLM on QEC and circuit targets at equal compute is, as far as I know, not settled. I could not check this.
- The month-by-month recall sweep this task asked for is still entirely undone. See the status box at the top.
