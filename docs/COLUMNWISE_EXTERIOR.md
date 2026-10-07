# Majorante exterior por colunas d*(G): pré-registro

Autor: sub-agente MATEMÁTICO (Claude), rodada 2 de 16/09/2026.

Este documento foi escrito ANTES de qualquer cálculo certificado desta
tarefa. As seções "Enunciado", "Regra de decisão" e "Critérios de
aceitação" são o pré-registro e NÃO serão alteradas depois de ver o
resultado. Mudanças posteriores vão em seções datadas ao final.

**Transparência sobre o que já tinha sido rodado antes deste registro.**
Dois diagnósticos em float, rotulados, que motivaram a tarefa:

1. `independent_column_diagnostic.py` (rodada 1): normas reais de coluna
   ≈5× abaixo do majorante uniforme; majorante "por grau" ≈3× abaixo em
   n=145.
2. Nesta rodada, antes de escrever este arquivo: o supremo em m de
   (n+1)·‖Q e_{m,n}‖ para Q=O_μ^-1. Nos pesos fortes, J→≈67.0 e K→≈33.6,
   contra o bound analítico 81 atual. Nos pesos fracos, J→≈127 e K→≈63.5,
   contra 153.

Nenhum desses números é prova. A regra de decisão abaixo é a mesma
pré-registrada na §(c) de `AUDITORIA_MATEMATICA_16SET.md` (rodada 1),
antes de qualquer um deles.

## Enunciado

Seja G a caixa |m|<M_G, n<N_G e e uma coluna unitária FORA de G
(m≥M_G ou n≥N_G), em qualquer componente e qualquer parte real ou
imaginária.

**Setor A.** Realização H(s)=L_true(s)Q_ref, Q_ref=O_{μ_RefA}^-1. Norma ℓ¹
ponderada de RT, (κ1,κ2)=(65/64,5/4) × η=(916,4096,243,1233), com soma
por coluna. Encontrar um racional d*_A(G) tal que, para toda coluna e
fora de G e todo |s|≤5/4 (cobre a faixa A e o retângulo do setor B),

    ‖(H_true(s)-I) e‖ ≤ d*_A(G) ‖e‖,

com o erro de μ (ε_μ=43·2^-40+2^-277, RT §5 e (60)) e a bola
ε_ω=2^-25+2^-277 ((39b) e (60)) incluídos.

**Sharp.** Realização H♯(s)=K_true(s)Q♯_ref, pesos (129/128,9/8), sem
pesos de componentes, disco |s-733/1000|≤1/32 (portanto |s|≤4/5). Mesmo
enunciado, com d*_♯(G).

Consequência que será usada: d*(G) majora ‖H_TT-I‖ e ‖H_ST‖ por ser um
supremo de colunas em ℓ¹. Substitui o d(G) das pontes GUARDED e sharp.

## Método (Fase 1, semianalítica)

(H_true-I)e = (B_ref+s)Qe + (B_true-B_ref)Qe + (μ_true-μ_ref)D'Qe.

A coluna de Q é EXPLÍCITA por RT (22): c_n=1/D_n;
c_{n-k}=-2μn/(D_n D_{n-k}); c_j=-f_j c_{j+k}.

    ‖(B_ref+s)Qe‖ ≤ Σ_j (w_j/w_n)(|Re c_j| β(d,0,m,j)+|Im c_j| β(d,1,m,j)) + |s|‖Qe‖,

onde β(d,p,m,j)=‖B_ref e_{d,p,m,j}‖/‖e_{d,p,m,j}‖ são normas EXATAS de
coluna de B_ref (RefA tem suporte 41×101). Para j≥101 a norma não depende
de j; para m≥41, não depende de m. Logo há só finitos valores distintos,
todos calculáveis exatamente. As partes (B_true-B_ref) e D' entram pelos
bounds já auditados (A2.4 e A3.5), multiplicados pelo ‖Qe‖ da coluna.

O supremo sobre infinitas colunas é fechado por bounds analíticos em n e
m com constantes explícitas. Não é um máximo sobre colunas amostradas.

Fase 2 (normas de coluna exatas numa faixa, estilo RT (53)) só se a
Fase 1 não atingir N_G≤384.

## Regra de decisão (pré-registrada)

Seja N_G o menor corte radial certificado com d*<1, com M_G escolhido
livremente e reportado.

- Se **N_G≤384 no setor A E no sharp**, a rota de blocos continua. O
  próximo gargalo passa a ser uma parametriz certificada da faixa U.
- Se **N_G>384** em qualquer dos dois, a rota de caixas é declarada
  inviável neste hardware (7 GB de RAM, 4 núcleos) e o trabalho passa a
  uma parametriz exterior analítica que absorva B (características).

Se a Fase 1 falhar e a Fase 2 não couber no orçamento de crédito desta
rodada, o resultado é registrado como "N_G certificado = valor da Fase 1".
A regra se aplica a esse valor. Diagnósticos em float não contam.

## Critérios de aceitação (pré-registrados)

1. d* ≥ norma exata (Fraction) em ≥50 colunas exteriores amostradas,
   incluindo m=0, m ímpar, K e J (teste).
2. d* ≤ 2× o diagnóstico float nessas colunas (medida de folga; se
   falhar, é registrado e não invalida a prova).
3. Relatório com a menor caixa com d*<1, uniforme em s e com a bola RT:
   `build/spectrum/columnwise-exterior.json`.
4. Reexecução do Schur GUARDED e da ponte sharp com essa caixa, com as
   novas margens registradas.
5. Esta regra de decisão, aplicada sem alteração.

## Nota sobre a norma do testemunho de constraints (revisão cruzada, parte (i))

O testemunho da raiz 0,401 proposto pelo Físico mede C(s)h com PERDA DE
RAIO: entrada no domínio (graph norm) e saída sharp em pesos menores
(SHARP_DOMAIN). d*_A é na norma forte η. Se o testemunho for enclausurado
em pesos fracos, será preciso um d*_A nos pesos fracos, com β- e sem
misturar normas. Não troco de norma nesta tarefa; registro a necessidade.
Detalhes em `AUDITORIA_MATEMATICA_16SET.md` §(d).

## Estado

(atualizado a cada marco; ver seções datadas abaixo)

### Marco 1 (16/09, rodada 2): implementação da Fase 1, antes do resultado

- `scripts/columnwise_exterior.py`: betas EXATOS (inteiros em 2^-72, pesos por
  potências inteiras) nos representantes j∈{101,102} (janela L≤49 do termo
  R), para todo m≤40 e m=41/42; βF(j) para j≤100 em m=41/42. Mais as cotas
  analíticas de Q:
  - radial: (n+1)q ≤ (C+1/(4C))/μ_ref, com C_J=1+2√2·r/(1-r) e
    C_K=1+2√2·(N/(N-1))·r²/(1-r²);
  - Fourier: ‖c_n‖_C ≤ ((1+√2)/2)/b, |c_{n-1}| ≤ 1/b e
    |c_{n-2}| ≤ (2/√3)/b;
  - racionais usados para as constantes irracionais: 99/70 ≥ √2,
    1.207107 ≥ (1+√2)/2 e 1.155 ≥ 2/√3.
- Invariâncias usadas, conferidas exatamente antes da execução completa:
  - β(A,d=2,p=1,m=42,j=102) = β(m=60) = β(j=140), igualdade racional;
  - β(j=151, janela completa) - β(j=101, L≤49) = 3,04750e-10 ≤ cauda R
    3,04753e-10.
- Previsão em FLOAT, antes do certificado: N_G≈520 no setor A (dominado por
  J, d=2,3, com β_∞≈6,5), contra o limiar 384. Registro isto para que a
  leitura do resultado não seja ajustada depois.

### Marco 2 (17/09): resultado CERTIFICADO da Fase 1 e plano da Fase 2

Relatório: `build/spectrum/columnwise-exterior.json`. Blocos A e sharp
completos: 1005 e 431 betas exatos. Todas as desigualdades finais são
racionais exatas.

| sistema | N_G (Fase 1) | d_R(N_G) | d_R(384) | M_G (Fase 1) | d_F(M_G) | β_R por componente |
|---|---:|---:|---:|---:|---:|---|
| A (η, 65/64, 5/4, \|s\|≤5/4) | **586** | 0,99949 | 1,524 | 107 | 0,99820 | 8,074 / 5,772 / 6,756 / 6,513 |
| sharp (129/128, 9/8, \|s\|≤4/5) | **606** | 0,99903 | 1,575 | 107 | 0,99864 | 2,501 / 3,518 |

Outras cotas radiais da Fase 1, para orientar a Fase 2:

| sistema | d_R(650) | d_R(700) | d_R(800) |
|---|---:|---:|---:|
| A | 0,901 | 0,837 | 0,733 |
| sharp | 0,932 | 0,865 | 0,757 |

A previsão em float do Marco 1 (≈520) ficou abaixo do valor certificado
(586). A diferença vem das constantes analíticas √2 e (1+√2)/2 e do máximo
em m de β_R.

**Aplicação provisória da regra:** a Fase 1 dá N_G>384 nos dois sistemas.
Pelo pré-registro, a Fase 2 é executada antes da decisão.

**Plano da Fase 2 (método; a regra não muda).** Escolho M_G=107, o valor
certificado da Fase 1 e permitido pela regra ("M_G escolhido
livremente"). Assim a faixa de Fourier fica vazia. G=(107, 384). O
exterior de G divide-se em:

1. n≥N_far: cota radial da Fase 1. N_far=586 em A e 606 no sharp.
2. m≥107: cota de Fourier da Fase 1.
3. **faixa radial** 384≤n<N_far, m<107: TODAS as colunas uma a uma, com
   enclausuramento superior exato de ‖B_ref Qe‖+(s_max+e_B)‖Qe‖+e_μ.
   - Coluna de Q diádica, com raio em módulo.
   - Convolução exata por inteiros (substituição de Kronecker).
   - Cauda de profundidade com β_glob.

Tamanho da faixa: 75447 colunas (A) e 35631 (sharp). Para m≥41 a norma não
depende de p, pela simetria das imagens desdobradas; restam ≈52k e ≈25k
avaliações. Checkpoint por linha n em
`build/spectrum/columnwise-phase2/`. Limites da rodada: 2 núcleos e
`ulimit -v 3000000`.

d*(G) = max(máx. da faixa, d_R(N_far), d_F(107)) ≈ 0,9995. Isso é < 1,
mas com margem pequena. Se couber, a faixa será estendida até N_far=700 e
M=130 para melhorar a margem. Isso é registrado à parte e não muda a
decisão.

### Marco 3 (17/09): Fase 2 implementada e validada; custo medido

`scripts/columnwise_phase2.py` calcula, para cada coluna da faixa radial, um
majorante racional EXATO de ‖(H_true(s)-I)e‖, uniforme em |s|≤s_max.
- Coluna de Q diádica (2^-64), com raio em módulo propagado.
- Convolução EXATA por substituição de Kronecker, mais o termo μR.
- Resto de profundidade com β_glob.
- Imagens desdobradas em n (n-D≥101) e em m (m≥41). Para m≥41, p=1 é
  idêntico a p=0.
- Checkpoint por linha n em `build/spectrum/columnwise-phase2/`.

**Validação cruzada:** duas implementações independentes do mesmo
majorante, `columnwise_validate.upper` (dicionários, 4 imagens, sem atalhos
de simetria) e `columnwise_phase2.column_upper` (Kronecker, atalhos de
simetria). Em 5 colunas (A: m=0, m=6 p=1, m=41; sharp: m=4 p=1, m=50) as
duas diferem por ≤8·10^-16. A diferença é sempre para cima na versão
rápida, por causa da cota de raio mais grosseira, o que é esperado.

**Custo medido:** ≈1,2–2,8 s por coluna na versão rápida, contra 50–250 s
na direta. A faixa tem 52116 (A, n∈[384,586)) + 24642 (sharp,
n∈[384,606)) = 76758 avaliações. Com 2 núcleos isso dá **≈21–30 h** de
relógio. Não cabe nesta rodada. Execução em curso: as linhas de frente
n=384 e 385, nos dois sistemas.

### Marco 4 (17/09): núcleo exato em C, primeira linha certificada

Criei `scripts/columnwise_kernel.c`: C + GMP, `__int128` com cotas de
estouro verificadas em tempo de execução e OpenMP em 2 threads. O driver é
`columnwise_phase2.py --engine c`.
- A coluna de Q continua EXATA em Python, arredondada a 2^-40 com raio.
- O núcleo devolve o inteiro exato T2 da soma ponderada.
- O driver reconstrói a fração e soma raio, massa de Q, e_B, e_μ e resto.

**Validação:** motor C contra as duas implementações Python nas mesmas 6
colunas (A e sharp, m=0, m dobrado com p=1, m≥41). Os valores do C são
maiores por ≈10^-8, pelo arredondamento em 2^-40 contra 2^-64, e ficam
sempre acima, como deve ser.

**Custo:** 516 colunas em 20 s de relógio, ≈0,04 s por coluna. A faixa
inteira (76758 avaliações) cabe em ≈1 h.

**Linhas certificadas (setor A):**

| n | colunas | máx. U |
|---:|---:|---:|
| 384 | 221 | 0,4592 |
| 385 | 295 | 0,5049 |

Em n=385 o máximo está em d=0 (K), m=18, p=1.

Em execução: faixa completa A n∈[384,586) e sharp n∈[384,606). Checkpoint em
`build/spectrum/columnwise-phase2/{A-M107-D60,sharp-M107-D100}.jsonl`.

### Marco 5 (17/09): faixa A n∈[384,586) COMPLETA

Checkpoint `build/spectrum/columnwise-phase2/A-M107-D60.jsonl`: 202 linhas
contíguas e 52116 avaliações. TODAS com U<1.
- Máximo racional: **0,50493**, em n=385, d=0 (K), m=18, p=1.
- O máximo decresce com n: 0,351 em n=553.

**Observação de margem (antes do resultado sharp):** com N_far=586, d_R=0,99949
fica acima do limiar de Schur GUARDED (≈0,99944). Por isso a faixa A está
sendo estendida até N_far=606 (d_R=0,9666). A extensão não altera a
decisão: d*<1 já vale com 586.

### Marco 6 (17/09): d*(G) final e DECISÃO pela regra pré-registrada

Caixa G=(M_G, N_G)=(107, 384) nos dois sistemas. Relatórios:
- `build/spectrum/columnwise-exterior-phase2.json`;
- bloco `phase2` em `build/spectrum/columnwise-exterior.json`.

d*(G) é o máximo exato de três cotas, cada uma cobrindo uma parte do
exterior:
1. faixa n∈[384,606), m<107, coluna a coluna (Fase 2);
2. n≥606, qualquer m: cota radial da Fase 1;
3. m≥107, qualquer n: cota de Fourier da Fase 1.

Todas são uniformes em |s|≤s_max e incluem ε_μ e a bola RT ε_ω.

| sistema | colunas na faixa | máx. na faixa | d_R(606) | d_F(107) | **d*(G)** | DOFs de G |
|---|---:|---:|---:|---:|---:|---:|
| A | 57276 (n=384..605) | 0,50493 (n=385, K, m=18, p=1) | 0,96656 | 0,99820 | **0,99820** | 143424 |
| sharp | 24642 (n=384..605) | 0,46300 (n=385, d=1, m=104) | 0,99903 | 0,99864 | **0,99903** | 61632 |

A faixa A foi estendida até 605, como previsto no Marco 5. Todas as
desigualdades d*<1 foram verificadas em racionais.

**N_G certificado = 384 no setor A e no sharp.** Pela regra
pré-registrada, aplicada sem alteração: **a rota de blocos continua.**

**Critério 4, margens novas** (escalares exatos; c e κ como em
`AUDITORIA_MATEMATICA_16SET.md` A2/A4):

| subsistema | caixa | defeito | margem | DOFs restantes U | antes |
|---|---|---:|---:|---:|---|
| GUARDED S+T, κ=31,4266 (recertificado) | 107×384 | 0,998763 | 1,24e-3 | 143091 | 0,997676 em 176×944, U≈5,79e5 |
| GUARDED S+T, κ=36,2267 (relatório) | 107×384 | 0,998849 | 1,15e-3 | 143091 | |
| ponte sharp F+T, κ=11,4601 | 107×384 | 0,999033 | 9,7e-4 | 61038 | 0,915 em 160×960, U=228366 |

**Leitura honesta (não altera a decisão).**
- A faixa U encolheu 4× no setor A e 3,7× no sharp. O corte radial caiu de
  944/960 para 384, e na própria faixa as normas reais ficam ≤0,51.
- Quem limita agora são as cotas ANALÍTICAS da Fase 1 além da faixa:
  d_F(107) no A e d_R(606) no sharp. Por isso as margens do Schur são de
  ≈10^-3.
- M_G=107 domina a contagem de DOFs. O diagnóstico em float da rodada 1
  indica normas reais ≈0,44 já em m=48. A próxima redução natural é uma
  faixa de Fourier coluna a coluna, mesmo motor, com m∈[48,107).
- Mesmo assim, U com 1,4·10⁵ (A) e 6,1·10⁴ (sharp) DOFs continua grande
  demais para uma parametriz densa em 7 GB. A rota continua pela regra,
  mas o próximo gargalo, a parametriz da faixa U, é pesado.

### Marco 7 (17/09): critérios de aceitação

1. **Cumprido, exato.** Script `scripts/columnwise_criteria.py`, saída
   `build/spectrum/columnwise-criteria.json`.
   - 59 colunas exteriores a G=(107,384), cobrindo m=0, m ímpar, K e J,
     setor A e sharp.
   - Regiões: faixa, n=700/1001 além da faixa, e m≥107 com n pequeno
     (1,2,3,5).
   - Em todas vale U(e)≤d*(G), em racionais. U(e) é um majorante da norma
     de coluna, uniforme em |s| e com a bola RT.
   - Motores: núcleo C exato, e a implementação Python direta com 4 imagens
     para n pequeno.
   - `test_columnwise_exterior.py` reconfere tudo em racionais a partir do
     disco em ≈1 s: faixa contígua com U<1, d* = máximo das três cotas,
     Schur <1 e o critério 1.
2. **Falhou como medida de folga de d*** (diagnóstico float, setor A, 47
   colunas).
   - d*/float ≤2 em NENHUMA coluna; o máximo é 17,3, nas colunas de Fourier
     com n pequeno, onde a norma real é ≈0,06.
   - Em compensação, o majorante POR COLUNA é justo: U/float ≤1,42.
   - A folga está no supremo: as cotas analíticas d_F e d_R dominam d*,
     enquanto as colunas reais ficam ≤0,51. Conforme o pré-registro, isso
     não invalida a prova.
3. **Cumprido:** `build/spectrum/columnwise-exterior.json`, blocos `phase2`
   e `acceptance_criteria`.
4. **Cumprido:** margens no Marco 6, via script novo
   `scripts/columnwise_phase2_report.py`. `guarded_exterior.py` e
   `sharp_guarded_bridge.py` não foram editados.
5. **Cumprido:** regra aplicada sem alteração → a rota de blocos continua.

**Hipóteses da prova:**
- bola publicada de RT: ε_ω por (39b) e (60), μ por §5 e (60), μ∈[1/6,17/100];
- lemas analíticos de colunas deste documento;
- invariâncias exatas usadas: β independente de j≥101 e de m≥41 (cauda R
  explícita); imagens desdobradas para n-D≥101 e m≥41, com p=1≡p=0 nesse
  regime.

As invariâncias foram conferidas em racionais em casos-teste (Marco 1) e
por duas implementações independentes (Marcos 3 e 4). A demonstração é em
papel, neste documento. Nada foi formalizado em Lean.

### Marco 8 (17/09, pré-registro): baixar d*(G) de 0,998 para ≤0,7

Escrito ANTES de disparar qualquer cálculo desta etapa. G continua
**107×384**; só o exterior muda. A regra de decisão da rodada 2 não é
tocada; este marco tem meta própria.

**Meta.** d*(G)≤0,7 no setor A e no sharp, uniforme em |s|≤5/4 (A) e
|s|≤4/5 (sharp), com ε_μ e a bola RT, tudo em racionais.

**Cortes novos**, escolhidos pelas escalas EXATAS da Fase 1 (d_F é
exatamente Φ/b+e_μ, logo ∝1/M; d_R decresce ≈ como 1/N):

| sistema | d_R(900) | d_F(160) |
|---|---:|---:|
| A | 0,65116 | 0,66755 |
| sharp | 0,67304 | 0,66784 |

Portanto **N_far=900** e **M'=160**. Com N_far=880/M'=157 também passaria;
escolho 900 e 160 por margem e por serem redondos.

**Regiões de colunas a verificar** (motor C exato, uma coluna por vez).
Com o que já está certificado, o complemento de G fica coberto SEM buracos:

| região | n | m | estado |
|---|---|---|---|
| R0 | [384,606) | [0,107) | **já certificada** (rodada 2) |
| R1 | [384,606) | [107,160) | nova |
| R2 | [606,900) | [0,160) | nova |
| R3 | [0,384) | [107,160) | nova |
| cauda radial | ≥900 | qualquer | d_R(900), Fase 1 |
| cauda de Fourier | qualquer | ≥160 | d_F(160), Fase 1 |

Verificação de cobertura: uma coluna fora de G tem m≥107 ou n≥384. Se
m≥160 → cauda de Fourier. Se n≥900 → cauda radial. Senão m<160 e n<900, e
então: n<384 ⇒ m∈[107,160) ⇒ R3; 384≤n<606 ⇒ R0∪R1; 606≤n<900 ⇒ R2.

**Ajuste de código** (arquivos meus). `columnwise_phase2.py` ganha
`--m-lo` e escolhe automaticamente o modo do núcleo:
- n−D≥101: modo desdobrado, já validado (imagens τ=+ com peso 2, σ=+ só
  para m≥41, p=1≡p=0);
- n−D<101 (região R3 com n pequeno): **modo completo**, sem atalho de
  simetria, com as quatro imagens (σ=±, τ=±) e pesos κ2^|N|. É uma extensão
  nova do núcleo C e será validada contra a implementação Python direta
  (`columnwise_validate.upper`, 4 imagens) antes de valer como certificado.

Os checkpoints novos vão para
`build/spectrum/columnwise-phase2/{sistema}-m{m_lo}-M{M}-D{D}.jsonl`, uma
linha por n, com `m_lo` gravado.

**Montagem final** (feita depois do cálculo, por um script meu novo,
`columnwise_exterior_assemble.py`): lê TODOS os checkpoints, confere em
racionais que as regiões cobrem o complemento de G sem buracos (linha a
linha, unindo intervalos [m_lo,M)), toma
d*(G)=máx(máximo das regiões, d_R(900), d_F(160)) com d_R e d_F EXATOS da
Fase 1, e recalcula as margens de Schur do GUARDED e da ponte sharp.

**Critério de sucesso (pré-registrado):** d*(G)≤7/10 nos dois sistemas,
com cobertura completa verificada e todas as linhas com U<1. Se o máximo
das regiões novas exceder 0,7, registra-se o valor obtido e a meta falha —
o que não invalidaria nada da rodada 2, porque d*<1 já está certificado.
