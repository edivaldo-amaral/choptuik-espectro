# Espectro da solução de Choptuik

Contagem rigorosa, assistida por computador, dos modos instáveis da solução crítica de Choptuik (colapso
gravitacional esférico de um campo escalar sem massa), na realização de Reiterer–Trubowitz (RT).

## O resultado: T2

Na realização espectral de RT, **a solução crítica de Choptuik tem exatamente um modo instável físico**,
módulo gauge. Ele é real, simples e está em **s ∈ [0,73316950; 0,73318982]**.
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
- conferência externa: γ = 1/λ ∈ [0,373939; 0,373950], o expoente crítico conhecido.

**Escopo e modelagem:**
- perturbações esfericamente simétricas;
- regularidade C^∞ em τ e analítica em ξ;
- o fundo é a solução de RT (Comm. Math. Phys. 368, 2019).

O teorema é **espectral**: não diz nada sobre a evolução nem sobre a dinâmica não linear. As hipóteses
antes abertas (completude do gauge, reconstrução e cone móvel) são agora lemas provados.

**Revisão.** Cada peça da cadeia passou por revisão independente: 10 revisões, por agentes de outro
modelo, com contexto limpo, snapshot sha256 e artefato por alegação, em [`.codex-runs/`](.codex-runs/).
Nenhuma quebrou um resultado.

## Como verificar

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
bash scripts/baixa_dados_rt.sh                      # dados publicados de RT (arXiv:1203.3766v1), com sha256
.venv/bin/python scripts/verifica_T2.py              # nível 1: refaz em racionais as verificações finais
```

O nível 1 refaz, a partir dos certificados gravados em `build/`:
- o θ de cada ponto dos discos dos Rouchés;
- as coberturas e as contagens;
- o Perron dos quatro NK;
- os 154 tiles de S3a;
- os símbolos e os lemas.

Leva ~6 h num processo (os discos dominam; `--disco` permite paralelizar) e menos de 1 min para o resto.
Os dados pesados (A, T, U e os fatores F, 16 GiB) estão no Zenodo, DOI 10.5281/zenodo.23222363 (aberto
junto com o artigo). O nível 2 (refazer os certificados a partir dos dados de RT) e o mapa de cada alegação para os seus
arquivos estão em [`docs/REPRODUCAO.md`](docs/REPRODUCAO.md) e [`docs/REPRODUCAO_MAPA.md`](docs/REPRODUCAO_MAPA.md).

## Por onde entrar

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

## Ferramentas antigas (setembro)

O pipeline de setembro continua disponível: `make setup`, `make formal`, `make spectrum-small`,
`make verify-winding`, `make test`. A descrição dele e o histórico até 22/09 estão em
[`docs/historico/README_2026-09.md`](docs/historico/README_2026-09.md).

## Créditos

- A solução, os dados `RefA` e `RefAplusB` e o código de existência são de Reiterer & Trubowitz,
  arXiv:1203.3766, *Choptuik's critical spacetime exists*, Comm. Math. Phys. 368 (2019). Não são
  redistribuídos aqui; ver `scripts/baixa_dados_rt.sh`.
- O desenvolvimento foi feito com assistência de modelos de IA. As revisões independentes foram feitas
  por agentes de outro modelo.

**Licenças:** MIT para o código ([`LICENSE`](LICENSE)); CC BY 4.0 para documentos, certificados e dados
([`LICENSE-DOCS.md`](LICENSE-DOCS.md)). Como citar: [`CITATION.cff`](CITATION.cff).
