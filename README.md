# Spectrum of Choptuik's critical solution / Espectro da solução de Choptuik

[English](#english) · [Português](#português)

## English

A rigorous, computer-assisted count of the unstable modes of Choptuik's critical solution (the spherically
symmetric collapse of a massless scalar field), in the realization of Reiterer and Trubowitz (RT), *Choptuik's
critical spacetime exists*, Comm. Math. Phys. 368 (2019) 143–186.

### Result

**Under spherically symmetric linear perturbations that are smooth in time and analytic in the radial variable up
to and across the light cone, and modulo gauge, Choptuik's spacetime has exactly one unstable mode.** It is real and
algebraically simple, with exponent **s\* ∈ [0.73316949, 0.73318983]**; in logarithmic time its growth rate is
λ = 2πs\*/K ∈ [2.674040, 2.674115], so 1/λ ∈ [0.373955, 0.373966], consistent with the measured critical
exponent γ ≈ 0.374.

Per Floquet class, the linearized operator L has exactly four roots with Re s > 0, all simple:

| root | nature |
|---|---|
| μ ≈ 0.16831 (exact) | gauge mode: translation of the time of the singularity |
| s_B ≈ 0.04142 | violates the constraints (sector B) |
| s_K ≈ 0.40102 | violates the constraints |
| **s\* ≈ 0.73318** | **the physical unstable mode** |

On the imaginary axis the only root of L is the phase translation s = 0. The theorem is **spectral**: it says nothing
about the linearized evolution semigroup or the nonlinear dynamics, and it concerns spherically symmetric
perturbations only.

### The paper

**Preprint:** E. A. Gonçalves, *Exactly one unstable mode: a computer-assisted spectral count for Choptuik's
critical spacetime under spherically symmetric perturbations*, version 1 (October 2026), Zenodo,
DOI [10.5281/zenodo.23251930](https://doi.org/10.5281/zenodo.23251930). The LaTeX source is in [`paper/`](paper/)
(`pdflatex main && bibtex main && pdflatex main && pdflatex main`). The large data files are on Zenodo,
DOI [10.5281/zenodo.23222363](https://doi.org/10.5281/zenodo.23222363).

### How to verify

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
bash scripts/baixa_dados_rt.sh                      # RT's published data (arXiv:1203.3766v1), sha256-checked
.venv/bin/python scripts/verifica_T2.py              # level 1: re-checks the final inequalities in exact arithmetic
```

Level 1 re-checks, from the certificate files in `build/` and RT's published data (downloaded and sha256-checked by
the script above), the Perron bound at every point of the contour discs,
the covers and the finite counts, the windows, the four Newton–Kantorovich certificates, the witnesses, the 154
tiles, the deflation polynomial and the structural bounds. It takes about two hours on four cores (`--disco`
parallelizes the discs). The certificates rationalize floating-point quantities (Perron vectors, approximate
inverses), so on a machine with a different BLAS or CPU the recomputed rationals differ in the last bits; the verifier
then checks that each recomputed certificate is valid, is not worse than the bounds published in the paper, and agrees
with the stored one to 1e-9 relative. `--estrito` requires bit-for-bit identity, which holds on the machines where the
certificates were generated. Level 2, which recomputes the certificates from RT's data, and the map from each claim to
its files are described (in Portuguese) in [`docs/REPRODUCAO.md`](docs/REPRODUCAO.md) and
[`docs/REPRODUCAO_MAPA.md`](docs/REPRODUCAO_MAPA.md); the per-claim sha256 digests are in Appendix C of the paper.

### Repository layout

| directory | contents |
|---|---|
| [`paper/`](paper/) | the manuscript, and the dated account of the use of AI assistants ([`paper/AI_USE.md`](paper/AI_USE.md)) |
| [`scripts/`](scripts/) | code (comments in Portuguese); `verifica_T2.py` is the verification program |
| [`build/`](build/) | certificate files (JSON) and reproduction reports |
| [`docs/`](docs/) | working documents, in Portuguese; [`docs/PANORAMA.md`](docs/PANORAMA.md) is the overview |
| [`.codex-runs/`](.codex-runs/) | the independent AI-assisted checks: tasks, reports and artifacts |
| [`formal/`](formal/) | Lean code from an earlier approach; it does not cover the final proof |

### Use of AI assistants and independent checks

The work was done with extensive use of AI assistants (Claude, through Claude Code, and OpenAI's Codex), which
wrote most of the code and drafted the documents and the manuscript under the author's direction; the author is
responsible for the content. Each certificate and lemma was checked by an AI agent running a different model, with
a clean context and one verifiable artifact per claim, including a review of the whole manuscript. These checks are
not a substitute for peer review. See Section 7 of the paper.

### Credits and licences

The solution, the data `RefA` and `RefAplusB` and the existence code are by Reiterer and Trubowitz
(arXiv:1203.3766); they are not redistributed here (see `scripts/baixa_dados_rt.sh`). Licences: MIT for the code
([`LICENSE`](LICENSE)); CC BY 4.0 for documents, certificates and data ([`LICENSE-DOCS.md`](LICENSE-DOCS.md)). How to
cite: [`CITATION.cff`](CITATION.cff).

## Português

Contagem rigorosa, assistida por computador, dos modos instáveis da solução crítica de Choptuik (colapso
gravitacional esférico de um campo escalar sem massa), na realização de Reiterer–Trubowitz (RT).

### O resultado: T2

Na realização espectral de RT, **a solução crítica de Choptuik tem exatamente um modo instável físico**,
módulo gauge. Ele é real, simples e está em **s ∈ [0,73316949; 0,73318983]**.
- O enunciado completo, as hipóteses e a cadeia da prova estão em [`docs/T2_ENUNCIADO.md`](docs/T2_ENUNCIADO.md).
- A visão de conjunto está em [`docs/PANORAMA.md`](docs/PANORAMA.md).

| raiz em Re s > 0 | natureza | certificado |
|---|---|---|
| μ* ≈ 0,16831 | gauge: translação do instante T* da singularidade | NK, raio 6·10⁻⁵ |
| s_B ≈ 0,04142 | viola as constraints (setor B) | NK, raio 3,7·10⁻⁵ |
| s_K ≈ 0,40102 | viola as constraints | NK, raio 3,8·10⁻⁴ |
| **≈ 0,73318** | **o modo físico** | NK, raio 2,2·10⁻⁵ |

**Outras afirmações:**
- não há raízes em Re s >= 2,93;
- no eixo imaginário, a única raiz de L é a fase s = 0 (simples, de gauge);
- conferência externa: γ = 1/λ ∈ [0,373955; 0,373966], o expoente crítico conhecido.

**Escopo e modelagem:**
- perturbações esfericamente simétricas;
- regularidade C^∞ em τ e analítica em ξ;
- o fundo é a solução de RT (Comm. Math. Phys. 368, 2019).

O teorema é **espectral**: não diz nada sobre a evolução nem sobre a dinâmica não linear. As hipóteses
antes abertas (completude do gauge, reconstrução e cone móvel) são agora lemas provados.

**Revisão.** Cada peça da cadeia passou por checagem independente assistida por IA: agentes de outro
modelo, com contexto limpo, snapshot sha256 e artefato por alegação, inclusive uma revisão do manuscrito
inteiro, em [`.codex-runs/`](.codex-runs/). Nenhuma quebrou um resultado. Não substituem a revisão por pares.

**O artigo** (em inglês) foi publicado como preprint no Zenodo, DOI
[10.5281/zenodo.23251930](https://doi.org/10.5281/zenodo.23251930); o fonte LaTeX está em [`paper/`](paper/): `pdflatex main && bibtex main && pdflatex main && pdflatex main`.

### Como verificar

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
bash scripts/baixa_dados_rt.sh                      # dados publicados de RT (arXiv:1203.3766v1), com sha256
.venv/bin/python scripts/verifica_T2.py              # nível 1: refaz em racionais as verificações finais
```

O nível 1 refaz, a partir dos certificados gravados em `build/` e dos dados publicados de RT (baixados e conferidos
pelo script acima):
- o θ de cada ponto dos discos dos Rouchés;
- as coberturas e as contagens;
- o Perron dos quatro NK;
- os 154 tiles de S3a;
- os símbolos e os lemas.

Leva ~6 h num processo (os discos dominam; `--disco` permite paralelizar) e menos de 1 min para o resto. Em outra
máquina (outro BLAS ou CPU), os racionais refeitos mudam no último bit; o verificador confere então que cada
certificado refeito é válido, não é pior que as cotas publicadas no artigo e fica a 1e-9 relativo do gravado.
`--estrito` exige igualdade bit a bit.
Os dados pesados (A, T, U e os fatores F, 16 GiB) estão no Zenodo, DOI 10.5281/zenodo.23222363. O nível 2 (refazer os certificados a partir dos dados de RT) e o mapa de cada alegação para os seus
arquivos estão em [`docs/REPRODUCAO.md`](docs/REPRODUCAO.md) e [`docs/REPRODUCAO_MAPA.md`](docs/REPRODUCAO_MAPA.md).

### Por onde entrar

| se você quer… | leia |
|---|---|
| o conjunto: resultado, cronologia, o que falta e os caminhos | [`docs/PANORAMA.md`](docs/PANORAMA.md) |
| o teorema e a cadeia da prova | [`docs/T2_ENUNCIADO.md`](docs/T2_ENUNCIADO.md) |
| conferir os certificados | [`docs/REPRODUCAO.md`](docs/REPRODUCAO.md) |
| a formulação do problema espectral | [`docs/SPECTRAL_PROBLEM.md`](docs/SPECTRAL_PROBLEM.md) |
| o que não funcionou, e por quê | [`docs/NAO_FAZER.md`](docs/NAO_FAZER.md), [`docs/METODO.md`](docs/METODO.md) |
| as perturbações não esféricas (numérica) | [`docs/NAO_ESFERICO.md`](docs/NAO_ESFERICO.md) |
| todos os documentos | [`docs/INDEX.md`](docs/INDEX.md) |

**Perturbações não esféricas (fora de T2).** A numérica, não rigorosa, reproduz Martín-García & Gundlach
(1999). O modo polar ℓ = 2 é amortecido, mas está a 0,0048 do eixo, o que faz dele o candidato a um
próximo teorema.

### Ferramentas antigas (setembro)

O pipeline de setembro continua disponível: `make setup`, `make formal`, `make spectrum-small`,
`make verify-winding`, `make test`. A descrição dele e o histórico até 22/09 estão em
[`docs/historico/README_2026-09.md`](docs/historico/README_2026-09.md).

### Créditos

- A solução, os dados `RefA` e `RefAplusB` e o código de existência são de Reiterer & Trubowitz,
  arXiv:1203.3766, *Choptuik's critical spacetime exists*, Comm. Math. Phys. 368 (2019). Não são
  redistribuídos aqui; ver `scripts/baixa_dados_rt.sh`.
- O desenvolvimento foi feito com uso extensivo de assistentes de IA; a declaração completa, com os modelos,
  está na §7.6 do artigo e em [`paper/AI_USE.md`](paper/AI_USE.md).

**Licenças:** MIT para o código ([`LICENSE`](LICENSE)); CC BY 4.0 para documentos, certificados e dados
([`LICENSE-DOCS.md`](LICENSE-DOCS.md)). Como citar: [`CITATION.cff`](CITATION.cff).
