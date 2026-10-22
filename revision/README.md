# Revisione QT10218 (PRX Quantum) — "Non-Local Magic Resources for Fermionic Gaussian States"

## File

| File | Contenuto |
|---|---|
| `prl/main.pdf` (`main.tex` + `endmatter.tex`) | Lettera in formato PRL: testo principale + End Matter (dopo le referenze) |
| `prl/supp.pdf` (`supp.tex`) | Supplemental Material = vecchie Appendici (S1 prova del Teorema 1, S2 verifica variazionale, S3 errori sperimentali, S4 generalizzazione α, S5 dettagli numerici SYK₂/circuiti) |
| `prl/response.pdf` (`response.tex`) | Lettera all'editor + lista delle modifiche + risposta punto per punto (stile Overleaf: referee in blu, risposte in nero, `[TO DO:]` in arancione) + tabella vecchia→nuova numerazione |
| `prl/build.sh` | Compila tutto (main prima di supp: i riferimenti SM→main usano `xr-hyper`) |
| `v5_long/main_v5.tex` | Versione lunga (formato originale) = versione sottomessa al giornale + tutte le modifiche di questo round |
| `v5_long/marked_up.pdf` | Diff (latexdiff) tra la versione **sottomessa al giornale** e `main_v5.tex`: aggiunte in rosso, cancellazioni barrate |
| `submitted/main_submitted.tex` | Sorgente della versione sottomessa al giornale (base del diff) |
| `v4_arxiv/` | Sorgenti arXiv 2604.27049v3 (riferimento) |

La base è la versione sottomessa al giornale (terminologia GNL / "Gaussian-restricted", paragrafo e nota su Fig. 2(a) con citazione Zenodo, grant DAAD, figura Odyssey ridisegnata, figure 1–2 con etichette GNL).

## Versione PRX Quantum (`prxq/`) — da usare per il resubmit/appeal

| File | Contenuto |
|---|---|
| `prxq/main_red.pdf` | Manoscritto con **tutte le modifiche in rosso** (cancellazioni omesse, come nella sottomissione precedente) |
| `prxq/main.pdf` (`main.tex`) | Manoscritto pulito (sorgente) = versione sottomessa + modifiche minime |
| `prxq/marked_up.pdf` | Diff con cancellazioni barrate (opzionale) |
| `prxq/response.pdf` (`response.tex`) | Lettera all'editor (appeal) + lista modifiche + risposte ai referee + Fig. R1 + tabella equazioni |
| `prxq/build.sh` | Compila tutto (main, marked_up, main_red, response) |
| `checks/small_N_check.py` | Fig. R1: campionamento Haar a N piccolo vs formula a N finito (peso corretto e stampato) |
| `email_editor.md` | Mail all'editor (richiesta di un terzo round) |

## Modifiche di contenuto di questo round

1. **Definizioni (Eq. 2–3) e Teorema 1 solo per i blocchi contigui** A={1,…,ℓ}, B={ℓ+1,…,N} (suggerimento del referee A): niente più T_π/π_AB né shift ciclico. Argomento aggiunto (End Matter "Local Gaussian orbit" + SM sotto Eq. S5): U_{O_A}⊗U_{O_B} ∝ U_{O_A⊕det(O_A)O_B}, quindi Eq. (3) = massimizzazione di ζ₂ sull'orbita (S5). Frase che dice che per bipartizioni non contigue l'identificazione fallisce (referee A).
2. **Peso di Jacobi**: x^{-1/2}(1-x)^{N-2ℓ} (typo), con P_n^{(-1/2,N-2ℓ)} e h_n corretti; limite termodinamico invariato.
3. **Appendice/SM** (referee A, punti 2–14): Eq. (124); spiegazione sotto (60); derivazione della (59); frase sul "transpose"; (80)–(81) semplificata; (148) rimossa; frasi riscritte; "invvariant"; **prova del Lemma 5 del referee**; "factorial normalization" rimosso; ‖·‖_HS; c_x^†; paragrafo di intuizione sulla prova.
4. Ringraziamenti ai referee e **dichiarazione sull'uso di AI (in rosso, da confermare)**.
5. Correzioni trovate nella revisione interna: espansione a tempi brevi nei circuiti (attorno a λ=1); unità di tempo delle formule del quench (evoluzione con H(μ)/2, autovalori di ih = ±2ε_k); soppressione nella fase ordinata solo per Ising; c_{N+1}=−c_1; didascalia Fig. 3 (ℓ=20,…,160); equazioni troppo larghe nel SM spezzate.

## Verifiche numeriche

- Peso di Jacobi: Monte Carlo O(2N) Haar (ℓ=1..4, N≤16): peso corretto compatibile, vecchio escluso; h_n a 1e-41; 2−log₂(2+√3) confermato.
- Teorema 1, blocco iniziale: U_{O_A}⊗U_{O_B} ∝ U_{O_A⊕det(O_A)O_B} (errore ~1e-15); bound mai violato e saturato; controesempio del referee confermato.
- Prova alternativa del Lemma 5: verificata fino a 1e-57.

## Da fare / controllare

- `[TO DO:]` arancioni in `response.pdf`: destinatario/giornale della lettera di accompagnamento; paragrafo sull'uso di AI; bipartizione contigua (metà catena) per Fig. 3(d); nessuna figura ha usato la formula a N finito col peso vecchio; metodo di campionamento in Fig. 1(a).
- Frase AI in rosso nei ringraziamenti di `main.tex` (e `main_v5.tex`).
- Unità di tempo del quench rispetto al codice di Fig. 3.
- **Figure**: `fig2_dynamics_4panel.pdf` (Fig. 3) ha ancora le etichette **FNL** sugli assi y (anche nella versione sottomessa): rigenerarla con GNL. In Fig. 2(a) e Fig. S2 la legenda "Eq.(6)" è coperta da un'etichetta TikZ "Eq.(4)": meglio rigenerare. In Fig. 2(c) la legenda dice β_GNL ln ℓ + b ma è il fit β_fit.
- Conteggio parole APS ~3340/3750; abstract 595 caratteri (≤600); End Matter ≈ 1.8 pagine (≤2).
