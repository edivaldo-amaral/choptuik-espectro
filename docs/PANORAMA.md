# Panorama do projeto Choptuik (03/10/2026)

Documento para decidir o próximo passo. Junta numa só visão tudo o que foi produzido de 14/09 a 03/10:
o resultado, como se chegou a ele, o que está revisado, o que falta e os caminhos possíveis.
- A cronologia detalhada está em `RETOMADA.md`.
- O ledger de obrigações está em `PROOF_OBLIGATIONS.md`.
- A sequência viva está em `PRIORIDADES.md`.

## 1. O resultado: T2

**Em uma frase.** Na realização espectral de Reiterer–Trubowitz (RT), a solução crítica de Choptuik do
colapso esférico de um campo escalar sem massa tem **exatamente um modo instável físico**, módulo gauge.
Ele é real, simples e está em s ∈ [0,73316950; 0,73318982]. A prova é assistida por computador e todas
as peças estão revisadas de forma independente (`T2_ENUNCIADO.md`).

**O espectro em Re s > 0, por classe de Floquet do período completo:**

| raiz | setor | natureza | certificado |
|---|---|---|---|
| μ* ≈ 0,16831 | A | gauge: translação do instante T* da singularidade (Lema C) | NK, raio 6·10⁻⁵ |
| s_B ≈ 0,04142 | B | viola as constraints | NK, raio 3,7·10⁻⁵ |
| s_K ≈ 0,40102 | A | viola as constraints | NK, raio 3,8·10⁻⁴ |
| **≈ 0,73318** | A | **o modo físico** | NK, raio 2,2·10⁻⁵ |

**Outras afirmações:**
- Nada em Re s >= 2,93 (F3).
- No eixo imaginário, entre as raízes de L, só a fase s = 0, que é simples e de gauge.
- As contagens são N_A = 3 e N_B = 1, por Rouché de operadores com referência finita, deflação de
  Brauer e contagem por Schur.

**Conferência externa.** γ = 1/λ, com λ = 4πλ*/Δ, cai em [0,373939; 0,373950]. É o expoente crítico
medido por Choptuik e calculado por Gundlach (0,374). A relação γ = 1/λ é heurística, da teoria de escala,
e não faz parte do teorema.

**O que o enunciado pressupõe (modelagem):**
- simetria esférica;
- H-reg′: perturbações C^∞ em τ e analíticas em ξ;
- o fundo de RT é a solução de Choptuik. RT provam existência e unicidade local na bola deles; o artigo
  foi publicado em Comm. Math. Phys. 368 (2019).

**Natureza.** O teorema é espectral: não diz nada sobre a evolução (semigrupo) nem sobre a dinâmica não
linear.

**Saíram da lista de hipóteses em 03/10:**
- H-gauge, pelo Lema G;
- H-rec, pelo Lema R e pela volta;
- H-cone, pelo Lema C, com as duas contagens: 1 com o gauge padrão, 2 com cone fixo.

## 2. Como chegamos

| período | fase | produtos principais |
|---|---|---|
| 14–16/09 | **Interlúdio** (outro modelo) | formulação, setor B, reconstrução, gauge local, transporte nulo. Auditado em 16/09 (`AUDITORIA_*_16SET.md`) e depois revisado peça a peça |
| 17–22/09 | **Viabilidade e mapa negativo** | 9 alavancas de barateamento fechadas, entre 1,0× e 2,7×; 5 rotas de certificado de contagem negativas; F3 de 414 para 6,66 (`RELATORIO_SETEMBRO_2026.md`, `NAO_FAZER.md`, `METODO.md`) |
| 23–24/09 | **S3 e F3 pelo símbolo** | F3: R_L <= 2,93; R♯ <= 1,765 (Arb); S3a: K injetivo por 154 tiles de Perron |
| 24–26/09 | **S3b e S4 (NK)** | Rouché de K e NK bordejado de L em s_K, com testemunho de constraints; S4: μ* simples, identificada com o gauge |
| 26/09–01/10 | **Mudança de rota: Rouché de operadores** | C1/T1 sem majorante de exterior: faixa A com 96 discos, N_A = 3 |
| 02/10 | **Faixa B** e o **inventário de T2** | N_B = 1; lemas L-real, L-afim, L-iso e L-constr |
| 03/10 | **Hipóteses viram lemas** | NK em 0,7332; Lemas G, R e C; volta; eixo imaginário; H-reg′ |
| 03/10 | **Não esférico, Fases 0 e 1** | as equações de Gerlach–Sengupta no fundo de RT; o espectro reproduz Martín-García & Gundlach (`NAO_ESFERICO.md`) |

**A virada de 26/09.** A rota do majorante de exterior pedia HPC e ordens de grandeza. Trocá-la pelo
Rouché com referência finita fez tudo caber em quatro PCs comuns: notebook, PC do laboratório,
laboratorio2 e PC 3.

## 3. Inventário

- **179 commits** (17/09 a 03/10), **58 documentos**, **180 scripts** e **33 módulos de teste**.
- **Certificados** em `build/`, com sha256 e manifestos:
  - tiles de S3a;
  - discos dos Rouchés das faixas A (v2) e B;
  - janelas e caudas;
  - NK de S3b, S4, raiz B e 0,7332;
  - os lemas de T2.
- **Revisões independentes** em `.codex-runs/`, por agentes de outro modelo, com contexto limpo, snapshot
  sha256 e artefato por alegação:
  - S3a (duas rodadas mais uma revisão rasa registrada), S3b, S4;
  - faixa A (C1), faixa B;
  - T2 (lemas e cadeia);
  - Lema G, Lema R (H-rec ida), volta;
  - itens de modelagem (CRH).

  São 10 revisões. **Nenhuma quebrou um resultado.** Todas acharam lacunas, a maioria de rastreabilidade
  e de enunciado, e todas foram corrigidas.
- **Lean** (`formal/`, ~4.200 linhas): é das rotas anteriores a 17/09 (cauda, winding, certificados
  racionais), com o ledger de confiança em `KERNEL_TRUST.md`. **Não cobre a cadeia final de T2.**
- **Numérica não esférica** (`ns_axial.py`, `ns_polar*.py`): ℓ = 2 e 3, axial e polar.

## 4. O que é contribuição e o que é herança

- **Herança:**
  - a existência do fundo e a bola de RT, com RefA a ~80 dígitos;
  - a compacidade;
  - a formulação ω de RT;
  - as técnicas: Rouché e Gohberg–Sigal, NK, Perron, Loewner, campo de valores, Toeplitz + Hankel,
    elipse de Bernstein. Todas clássicas, como já registrava `ESTADO_DA_ARTE.md` em 22/09.
- **O que parece ser nosso:**
  - **a primeira contagem rigorosa dos modos instáveis da solução de Choptuik.** Até aqui só havia
    numérica, por exemplo Gundlach (1997) e Martín-García & Gundlach (1999);
  - o enclausuramento rigoroso do autovalor físico, |λ* − 0,73318| <= 10⁻⁵;
  - a classificação completa das raízes em Re s > 0: duas violam as constraints, uma é gauge, uma é
    física;
  - a ponte geométrica completa, com gauge, reconstrução e realização;
  - o mapa negativo, com as constantes.
- **Apresentação honesta:** um resultado novo com ferramentas conhecidas, aplicadas a um operador em que
  não tinham sido aplicadas, mais o mapa do que não funciona.

## 5. O que falta ou está fraco

| item | natureza | gravidade para T2 |
|---|---|---|
| NK em 0,7332 sem revisão própria | a mesma cadeia da faixa B, que foi revisada | baixa: T2 só usa a contagem |
| passos 1–3 e 12 só em papel | auditados e revisados, mas sem Lean | média para um árbitro exigente |
| pacote de reprodutibilidade | certificados espalhados em `build/` de várias máquinas, logs, manifestos; falta um verificador único "refaz tudo e confere" | **alta para publicação** |
| H-reg′ em ξ (C^∞ em vez de analítico) | problema de regularidade aberto | modelagem |
| neutros físicos em s = 0 (escala, translação de φ) | os Lemas G e R não cobrem Re s = 0 | só se o enunciado falar de neutros |
| evolução e não linear | fora do escopo espectral | é a distância entre T2 e a lei γ = 1/λ |
| não esférico | só numérico; ℓ = 2 polar a 0,0048 do eixo | é a distância entre T2 e "1 modo" em geral |

**Documentos desatualizados**, que dizem coisas anteriores a T2 e induzem a erro quem lê:
- `ESTADO_DA_ARTE.md` ("T2 segue aberta");
- `RELATORIO_SETEMBRO_2026.md`;
- `PROBLEMAS_ABERTOS.md`;
- `INDEX.md`;
- partes de `S4_GAUGE.md` e de `T2_HIPOTESES.md`, cujas hipóteses viraram lemas.

Precisam de cabeçalho de atualização ou de reescrita.

## 6. Caminhos possíveis

| caminho | o que é | custo | valor | depende de |
|---|---|---|---|---|
| **A. Artigo de T2** | escrever o resultado: enunciado, cadeia, certificados, lemas, revisões | 1–2 semanas de escrita | o resultado principal publicado | idealmente C antes |
| **B. T3 não esférico** | contagem rigorosa: nenhum modo instável com ℓ >= 1 | várias semanas; do tamanho do projeto esférico | resolve a questão "inconclusiva" que RT citam; ℓ = 2 polar a 0,0048 do eixo | formulação bem condicionada, vínculos, NK em ℓ = 2, cota em ℓ grande |
| **C. Consolidação e reprodutibilidade** | verificador único, pacote de certificados com sha, guia de reprodução, atualização dos documentos desatualizados | ~1 semana; pouco cálculo | condição para A ser aceito; protege o que já existe | — |
| **D. Lean da cadeia final** | formalizar os verificadores racionais (Perron, cobertura, contagem) e os lemas | longo | confiança máxima | C |
| **E. Matemática além do espectral** | H-reg′ em ξ, semigrupo, não linear | pesquisa aberta | liga T2 à física (γ) | — |
| **F. Limpeza** | apagar `~/choptuik_trabalho` do PC 3; arquivar resultados das máquinas | horas | — | — |

**Sugestão de sequência:**
1. C, com a limpeza F junto.
2. A.
3. B, como segundo projeto, com a Fase 1 já feita como estudo de viabilidade.

D e E ficam como linhas de longo prazo. A justificativa: T2 está completo e revisado, mas hoje a única
forma de alguém de fora conferir é refazer conosco, máquina a máquina. Um pacote de reprodutibilidade
transforma o resultado em algo publicável, e o artigo então se apoia nele.

## 7. Números de referência

- **Raízes:**
  - μ* = 0,168307078963…;
  - s_K ≈ 0,40102;
  - s_B ≈ 0,04142;
  - λ* ∈ [0,73316950; 0,73318982].
- **Contagens:** N_A = 3 (96 discos, θ <= 0,99842); N_B = 1 (78 discos efetivos).
- **Cotas:** F3: R_L <= 2,93; R♯ <= 1,765; K injetivo por 154 tiles (S3a), mais o Rouché de K no disco de
  S3b.
- **Não esférico:**
  - polar: ℓ = 2, s = −0,004764 ± 0,147345i (κΔ = −0,0599); ℓ = 3, κΔ = −1,655;
  - axial: ℓ = 2, −2,320; ℓ = 3, −3,292.
