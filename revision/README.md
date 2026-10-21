# Revisione QT10218 — "Non-Local Magic Resources for Fermionic Gaussian States"

## Cosa c'è qui

| Cartella / file | Contenuto |
|---|---|
| `prl/main.tex` → `main.pdf` | Lettera PRL (testo principale + End Matter, `endmatter.tex`, dopo le referenze) |
| `prl/supp.tex` → `supp.pdf` | Supplemental Material = vecchie Appendici (S1 prova del Teorema 1, S2 verifica variazionale, S3 errori sperimentali, S4 generalizzazione α, S5 dettagli numerici SYK₂ / circuiti) |
| `prl/response.tex` → `response.pdf` | Risposta punto per punto ai referee A e B, con tabella vecchia→nuova numerazione |
| `prl/build.sh` | Compila tutto (main va compilato prima di supp: i riferimenti SM→main usano `xr-hyper`) |
| `v5_long/main_v5.tex` | Versione lunga (formato originale) con le sole correzioni di contenuto |
| `v5_long/diff_v4_v5.pdf` | latexdiff tra arXiv v3 (`main_v4.tex`) e `main_v5.tex`: mostra tutte le modifiche di contenuto |
| `v4_arxiv/` | Sorgenti arXiv 2604.27049v3 originali (riferimento) |

## Modifiche di contenuto (minime)

1. **Teorema 1** ristretto a bipartizioni in due blocchi contigui della catena di Jordan–Wigner, A={1,…,ℓ}, B={ℓ+1,…,N} (o viceversa). Aggiunto l'argomento per cui allora U_{O_A}⊗U_{O_B} = U_{O_A⊕det(O_A)O_B} a meno di una fase (orbita (50) ≡ definizione (5)), più l'osservazione che per bipartizioni non contigue l'identificazione fallisce (controesempio del referee A).
2. **Peso di Jacobi** (ex Eq. 11): x^{-1/2}(1-x)^{N-2ℓ} (era x^{-1/2}(1-x)^{N-2ℓ-1/2}, typo); corretti di conseguenza i parametri dei polinomi di Jacobi P_n^{(-1/2,N-2ℓ)} e h_n. Il limite termodinamico (incl. 2−log₂(2+√3)) non cambia.
3. **Appendice / SM** (referee A, punti 2–14): Eq. (124) corretta; spiegazione sotto (60); derivazione esplicita della (59) nel Lemma 1; frase sul "transpose" riscritta; prova (80)–(81) semplificata; (148) rimossa; frasi "LLM-ish" riscritte; typo "invvariant"; **prova del Lemma 5 sostituita con quella (più semplice) suggerita dal referee**; rimosso "with no factorial normalization"; notazione ‖·‖_HS uniforme; c_x^† = Σ x_j c_j^†; paragrafo di intuizione sulla logica della prova.
4. Ringraziamenti ai referee e **dichiarazione sull'uso di strumenti AI (in rosso, da confermare/modificare dagli autori)**.

## Verifiche numeriche fatte

- Peso di Jacobi: Monte Carlo su O(2N) Haar (ℓ=1..4, N≤16, 2·10⁵ campioni): peso corretto compatibile (|z|≤2.6), peso vecchio escluso (|z| fino a ~300); h_n verificato a 1e-41; limite termodinamico 0.1000313730 = 2−log₂(2+√3).
- Teorema 1 contiguo: identità operatoriale U_{O_A}⊗U_{O_B} ∝ U_{O_A⊕det(O_A)O_B} verificata (errore ~1e-15) per N=4,5; ζ₂ ≤ ∏(1−λ²+λ⁴) mai violato e saturato dall'ottimizzazione; controesempio del referee (A={1,3}) confermato.
- Prova alternativa del Lemma 5 e tutti i passaggi intermedi: verificati fino a 1e-57 (mpmath, 60 cifre).

## Da controllare da parte degli autori

- La dichiarazione sull'uso di AI (Acknowledgments, testo in rosso) e il paragrafo corrispondente nella risposta al referee A.
- I due `[...]` in rosso in `response.pdf` (bipartizione contigua nei dati di Fig. 3(d); metodo di campionamento in Fig. 1(a)).
- Conteggio parole stimato con le regole APS: ~3290/3750 (testo 2607 + 13 righe di equazioni + 3 figure). End Matter ≈ 1.75 pagine (≤ 2). Abstract: 575 caratteri (≤ 600).
