# Revisão independente do manuscrito de T2 (08/10/2026)

Você é o revisor independente de um artigo de matemática assistida por computador, que vai ser submetido à
Communications in Mathematical Physics. O autor do texto é outro modelo. Seu trabalho é **tentar quebrar** o
artigo: achar afirmações falsas, números que os certificados não sustentam, provas com lacunas, hipóteses
escondidas e afirmações maiores do que o que se prova. Um "OK" sem artefato conferível (script mais saída, ou
derivação escrita passo a passo) não conta. Escreva o relatório em português do Brasil.

## O objeto

O manuscrito em `paper/`: `main.tex`, `macros.tex`, `refs.bib` e `sec/*.tex`. O PDF compilado está nesta pasta
(`main_revisado.pdf`, 33 páginas). O resultado é o Teorema Principal da §3: na realização de Reiterer–Trubowitz
(RT), a solução de Choptuik tem exatamente um modo instável físico sob perturbações esfericamente simétricas,
módulo gauge, mais a classificação das quatro raízes de L em Re s > 0.

Partes do trabalho por trás do artigo já passaram por checagens independentes. **Não as refaça do zero.** Elas
são: as contagens de Rouché, os NK em s_K, μ e s_B, os tiles de S3a, os lemas G, R e a volta. Foque em:
1. **fidelidade:** o artigo diz exatamente o que os certificados e as provas sustentam?
2. **números:** cada número certificado impresso confere com os arquivos e está arredondado para o lado seguro?
3. **as provas como estão escritas no artigo:** estão completas e corretas? Em especial as partes reescritas
   recentemente (listadas em R4–R6).
4. **consistência interna:** definições antes do uso, notação, referências cruzadas, enunciados iguais nos
   lugares em que se repetem.

## Regras

- **Pode ler:**
  - todo `paper/`, exceto `paper/ESTRUTURA.md`;
  - `scripts/*.py` e os certificados em `build/**` (JSON e TXT);
  - o artigo, o código e os dados de RT em `.cache/rt-1203.3766v1/` (o fonte TeX está no tarball
    `arXiv-1203.3766v1.tar.gz`, arquivo `choptuik.tex`);
  - o PDF de Gundlach (1997) nesta pasta;
  - como **alegações a conferir**, estes documentos técnicos: `docs/SPECTRAL_PROBLEM.md`, `docs/SECTOR_B.md`,
    `docs/RADIUS_RECOVERY.md`, `docs/SHARP_LINEARIZED_IDENTITY.md`, `docs/SHARP_DOMAIN.md`,
    `docs/GAUGE_RT_ACTION.md`, `docs/GAUGE_GLOBAL.md`, `docs/HREC_IDA.md`,
    `docs/GEOMETRIC_FLOQUET_RECONSTRUCTION.md`, `docs/GLOBAL_NULL_GAUGE_TRANSPORT.md`,
    `docs/S3B_ROTA.md`, `docs/S4_GAUGE.md`, `docs/S3_CONDICIONAMENTO_SHARP.md`;
  - a web, para conferir afirmações sobre a literatura.
- **NÃO pode ler** (são narrativa e conclusões do autor ou de checagens anteriores):
  - `docs/C1_REAVALIACAO.md`, `docs/T2_ENUNCIADO.md`, `docs/T2_HIPOTESES.md`, `docs/PANORAMA.md`,
    `docs/PROOF_OBLIGATIONS.md`, `docs/RETOMADA.md`, `docs/PRIORIDADES.md`, `docs/REPRODUCAO*.md`,
    `docs/AUDITORIA_*.md`;
  - nada em `.codex-runs/` além desta pasta, **exceto em R9**, e só para o que R9 pede;
  - a memória do autor e `~/.codex`.
- **Não altere nenhum arquivo existente.** Escreva só nesta pasta (`.codex-runs/2026-10-08-revisao-artigo/`).
- Ambiente: `.venv/bin/python` (numpy, scipy, python-flint) e `CHOPTUIK_EPS_FUNDO_L=4.7e-08`. No notebook (7 GB,
  4 núcleos) use `OPENBLAS_NUM_THREADS=2`. Se precisar de cálculo pesado, pode usar o **laboratorio2**
  (`ssh -o BatchMode=yes <usuario>@<host>`, projeto em `~/Choptuik`, escreva só em `~/revisao_artigo/`).
  Não use outras máquinas. Trabalhos longos: `setsid nohup ... &` e esperas com `until ...; do sleep 30; done`
  (os processos em segundo plano do harness morrem em ~30 min). Mate processos pelo PID.
- Snapshot: `.snapshots/2026-10-08-revisao-artigo/MANIFEST.sha256` (commit em `COMMIT.txt`). Confira-o no
  início e no fim (`LC_ALL=C sha256sum -c ... | grep -v ': OK$'`).
- Grave resultados parciais e vá atualizando `RELATORIO.md`. Orçamento: até ~4 h.

## Artefatos obrigatórios

- **R1 (auditoria dos números).** Extraia do manuscrito **todo** número que se apoia num certificado e compare
  com os arquivos. Entram:
  - enclausuramentos das raízes;
  - θ, Z2, Y, r e o erro de λ dos quatro NK;
  - cotas do símbolo e da parte livre;
  - θ máximos dos discos e quantidade de discos e pontos;
  - ρ das janelas;
  - ε_B, t e t·ε_B;
  - testemunhos e a distância da identificação de μ;
  - constantes do L-real e da recuperação de raio;
  - F1;
  - os intervalos de λ e 1/λ (refaça em Arb, com o K de RT).

  Entregue um script que lê os JSON e imprime uma tabela "impresso | exato | arredondado para o lado
  seguro?". Toda cota superior deve estar arredondada para cima e toda inferior para baixo. Aponte cada caso
  em que não está.
- **R2 (teorema contra prova).** Para cada afirmação do Teorema Principal e da versão informal da §1, indique
  onde ela é provada (lema, proposição ou certificado). Diga se a prova cobre exatamente a afirmação, com as
  mesmas hipóteses, ou se falta algo. Faça o mesmo para o resumo.
- **R3 (§2 contra RT).** Confira no fonte de RT:
  - as fórmulas usadas no artigo: D^σ, as variáveis ω, D_O, H_μ e o fato de ∂τ só aparecer em O_μ, a
    identidade ♯, e as cotas (eq:RTball) com 2⁻²⁷⁷ e 2⁻²⁵;
  - as seções citadas ("§2", "§3", "§3–4", "§4"), conferindo que cada citação aponta para a seção certa do
    arXiv v1.

  Refaça a prova do Lema 2.2 (L e C afins) e confira C1 = (0, −ξh2) (apêndice B) pelas fórmulas de RT.
- **R4 (§4, partes reescritas).** Escreva a prova completa ou ache a falha em:
  1. **Lema 4.5 (completude do gauge):** a média em SO(3), a conclusão módulo campos de Killing e os passos 1–6,
     inclusive a variante (b);
  2. **o passo R7 do Lema 4.9:** a fórmula e^{sτ}Σ_j (τ^j/j!) h_{k−1−j} e a dedução das constraints de jato.
     Teste numericamente com uma família afim pequena (matrizes) e um operador C afim;
  3. **a Proposição 4.12** (multiplicidade física), inclusive o cenário de cone fixo (Definição 4.3): a
     fixação de gauge fica tangente ao cone? A obstrução κ some? A correspondência vale com G_s = 0?
  4. **o Teorema 4.11** (correspondência): injetividade, sobrejetividade e a passagem ao quociente.
- **R5 (§5).** Confira, com prova ou contraexemplo:
  - o enunciado do Rouché de Gohberg–Sigal na forma de homotopia;
  - o Lema 5.4 (Brauer: o determinante 2×2 e a fórmula de multiplicidade);
  - o Lema 5.5 (bloco finito);
  - o Lema 5.6 (critério de Perron), inclusive o papel do vetor de Perron e da monotonicidade em t;
  - o Lema 5.7 (janelas, princípio do máximo);
  - o Lema 5.8 (contagem finita).

  Diga se o Teorema 5.1 (ii)–(iii) segue dos lemas e dos certificados tal como estão descritos no texto.
- **R6 (§6).** Confira:
  - o Lema 6.1 (bordejado ⇔ simples);
  - a forma do teorema de NK (a fórmula de r e o raio de unicidade);
  - a Proposição 6.3 (μ identificado);
  - a Proposição 6.5: a região certificada cobre de fato tudo o que o Corolário 6.6 usa? Atenção ao disco
    de raio 0,05 contra o disco de S3b, e aos semiplanos Im ≥ 0 e Im ≤ 0;
  - a Proposição 6.7 (realidade: em que espaço a conjugação é isometria, e o uso do Lema 2.6);
  - a montagem final da prova (§6.6).
- **R7 (§7 contra os arquivos).** Confira as afirmações factuais da §7 contra `build/reproducao/*`,
  `build/reproducao/digests.json` e os scripts:
  - os 175 pontos e os 174 discos;
  - o pior θ;
  - o nível 2 idêntico;
  - a tabela de certificados;
  - os digests do apêndice C (refaça com `scripts/digest_certificados.py` e confira).
- **R8 (literatura).** Confira:
  - a citação literal de RT na §1.2;
  - os valores de Gundlach (λ₁, γ, a convenção de sinal) no PDF desta pasta;
  - γ ≈ 0,37 de Choptuik;
  - os valores de Martín-García e Gundlach (1999) citados na §8.2 (κΔ ≈ −0,07 no ℓ = 2 polar e −2,30 no
    ℓ = 2 axial);
  - o resumo de Reiterer (2015, arXiv:1510.05310), e se a frase da §8.1 o descreve corretamente;
  - os metadados da bibliografia que puder.
- **R9 (§7.5 e §7.6 contra o repositório).** Só depois de R1–R8, confira as afirmações factuais sobre o
  processo contra o que existe no repositório:
  - as dez checagens e o que cada uma cobriu, lendo **só** as primeiras linhas (título e escopo) dos
    `RELATORIO.md` em `.codex-runs/*revisao*`;
  - os identificadores de modelo citados, contra `paper/AI_USE.md` e os trailers do `git log`.
- **R10 (síntese).** Uma lista única de achados, cada um com:
  - gravidade: **erro** (afirmação falsa ou não sustentada), **lacuna** (prova incompleta, mas reparável),
    **imprecisão** (enunciado ou número a ajustar) ou **cosmético**;
  - local (arquivo e linha, ou seção);
  - evidência (o artefato);
  - correção sugerida.

  Termine com um veredito: o artigo, tal como está, sustenta o Teorema Principal?

## Entrega

`RELATORIO.md` nesta pasta, com R1–R10 em seções, mais os scripts e saídas `rN_*`. No fim, a conferência do
manifesto do snapshot.
