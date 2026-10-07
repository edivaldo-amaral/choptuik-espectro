# Auditoria matemática de 16/09/2026: cauda, contração e Schur

Autor: sub-agente MATEMÁTICO (Claude). Escopo: os objetivos 1 e 3 do
`prompt.MD`, ou seja, a cadeia de cauda, contração e Schur escrita entre
14/09 16:40 e 16/09 13:27. Nenhum arquivo existente foi editado. Não usei
`PROOF_OBLIGATIONS.md` nem `RETOMADA.md` como evidência. Serviram só para
mapear as alegações.

Regra usada: um veredito CONFIRMADO exige recomputo independente, em
`Fraction` ou inteiros exatos, ou uma rederivação em papel escrita aqui.
Ponto flutuante aparece só como proposta (inversas candidatas), sempre
revalidada com resíduo exato, ou como DIAGNÓSTICO rotulado (§c).

Numeração de RT: compilei o fonte arXiv v1
(`.cache/rt-1203.3766v1/arXiv-1203.3766v1.tar.gz`) e li o `.aux`. Rótulos
usados: (22) inversas de O^A, O^B; (23) S Ξ Γ2 expandido; (33)
estimativa de Γ2; (34) Lipschitz 3; (35) Band/BandI; (36) K1; (37a,b) K2;
(38); (39a–f) L1…L6; (47) Γ2♯; (52a–d) estimativas exteriores; (53)
EXTERIOR_ESTIMATE; (60) R=2^-277. Seções: §3 séries de Fourier–Chebyshev
(e "The ♯ system"), §4 análise, §4.1 estimativas exteriores, §5 valores dos
parâmetros.

## 0. O que foi recomputado e o que foi só lido

**Recomputado de forma independente.** Os scripts novos não importam o
código auditado, salvo onde indicado como comparação.

| Script novo | O que recomputa | Tempo |
|---|---|---|
| `scripts/independent_free_inverse.py` | colunas EXATAS de Q=O_μ^-1 por substituição regressiva a partir de (22); ‖Q F‖ exata; rederivação da cauda livre com dominância sobre colunas exatas; defeitos de RADIUS_RECOVERY; versão só com constantes RT | ~60 s |
| `scripts/independent_component_beta.py` | majorante 4×4 de 2SΞΓ2(RefA,·) a partir de (23) do TeX (não de `GAMMA2_macro.c`); β_η; checagem de L1=10.5; linhas ‖THF‖ de GUARDED; majorante ♯ a partir de Γ2♯ (§3) | ~10 s |
| `scripts/independent_selected_matrix.py` | reconstrói O_μ+μSΓ1+2SΞΓ2(RefA,·) do TeX e compara TODAS as entradas de `rt-A-4x12.dat` e `rt-A-6x18.dat` | 10–20 s |
| `scripts/independent_sharp_matrix.py` | idem para `rt-sharp-8x24.dat` e `rt-sharp-12x36.dat` (exportador C novo) | 7–29 s |
| `scripts/independent_crown_tiles.py` | recertifica os 51 tiles da coroa com Q_ref exata, V própria, erros de arredondamento exatos e orçamento do fundo derivado de RT; checa a geometria | ~7 min, 2 núcleos |
| `scripts/independent_sharp_disk.py` | recertifica o disco sharp 8×24 (inversa à esquerda) e 12×36 (precondicionada à direita) com fundo | 11 s / ~3 min |
| `scripts/independent_schur_scalars.py` | Schur S+T (A2), ponte sharp F+T e graus de Neumann (A4), a partir dos valores acima | <1 s |
| `scripts/independent_column_diagnostic.py` | DIAGNÓSTICO float (não é prova): normas reais de coluna de (B+s)Q versus o majorante uniforme | ~2 min |
| `scripts/test_independent_audit.py` | 6 testes rápidos (<2 s) | |

**Só lido, sem recomputo:** tabelas §§4–5 de `FREQUENCY_BLOCKS.md`
(superadas); a versão SEM fundo de `crown-contour-6x18.json`;
`SHARP_LINEARIZED_IDENTITY.md`, `GAUGE_RT_ACTION.md` e
`CONSTRAINT_ROOT_QUOTIENT.md` (escopo do Físico, que são dependências de
A1.8 e A4.8); `explore_exterior_parametrix.py`; `check_free_preconditioner.py`
fora da função `input_tail`.

Lean: compilei `ResolventCompression.lean`, `AlgebraicTail.lean` e
`FrequencyBlocks.lean` com `lake env lean`. Todos passaram. Os axiomas são
só `propext, Classical.choice, Quot.sound`. Há zero `sorry`, `axiom`,
`admit` e `native_decide`.

## (a) Tabela de vereditos

Resultado geral: **nenhum item do meu escopo foi REFUTADO.** Todo número
que recomputei reproduz o valor alegado ou fica abaixo dele, com majorante
válido. Os problemas encontrados são de alcance e utilidade (SOBREVENDIDO),
não de aritmética.

### A1. Recuperação de raio e F1

| # | Alegação | Arquivo | Veredito | Evidência do recomputo |
|---|---|---|---|---|
| A1.1 | Q=J_{μ_true}^-1=diag(O^B,O^A,O^A,O^A)^-1 é limitada em Y+ (65/64,5/4) e Y- (129/128,9/8), coincide na interseção e preserva caixas | RADIUS_RECOVERY | CONFIRMADO | RT (22) e (37a): ‖Q‖≤2(1+r)/((1-r)μ), que dá 18/μ e 34/μ. Q preserva m e só baixa n, logo QG=GQG. A fórmula de coluna é a mesma nos dois espaços |
| A1.2 | defeitos 117774/512125 (forte) e 8046603/10242500 (fraco) em 1024×4096, \|s\|≤5/4 | radius-transfer.json | CONFIRMADO | recomputados exatamente (`declared_matches_report: true`). Cauda forte 81/4097, fraca 153/4097=9/241 |
| A1.3 | cauda de entrada ‖Q(I-G)‖≤(3/2)max(Fourier, radial); 15/M,81/(N+1) forte e 27/M,153/(N+1) fraco | STRUCTURED_EXTERIOR, `sharp_exterior.free_tail` | CONFIRMADO | rederivação em papel no §A1 abaixo. As colunas EXATAS amostradas (8 m × 24 n além do corte, μ∈{1/6,μ_RefA,17/100}, 3 caixas) ficam todas abaixo do bound (máx. 2.85 contra 4.26) |
| A1.4 | β+=‖B_true‖_η≤5191/500 | COMPONENT_EXTERIOR | CONFIRMADO | majorante 4×4 refeito a partir de (23) do TeX, sem o macro C: igualdade racional entrada a entrada com `component-exterior.json`. β_η=10.3810817…≤10.382. Checagem cruzada: o máximo dos 7 termos de (33) em RefA é 10.136≤L1=10.5 de RT §5 |
| A1.5 | β-≤(162/85)β+ | RADIUS_RECOVERY | CONFIRMADO | normas de campo decrescem com os pesos; ‖Ξ^-1_reg‖ passa de 40/9 a 144/17 por (35); 144/17÷40/9=162/85 |
| A1.6 | "existência forte + unicidade fraca" ⇒ kernels e cadeias de Jordan coincidem em \|s\|≤5/4 | RADIUS_RECOVERY | CONFIRMADO (papel) | argumento refeito no §A1. Usa só: Gx de suporte finito, H limitado em Y+, contração de T(H-I)T nos dois espaços e o fato de os operadores coincidirem em Y+⊂Y- |
| A1.7 | Q compacta; H(s)=I+(B+s)Q é Fredholm analítica de índice 0; domínio fechado | RADIUS_RECOVERY, ledger F1 | CONFIRMADO (papel), com a precisão de (b1) | ‖Q-QG‖=‖Q(I-G)‖≤t(G)→0, com QG de posto finito. B É limitado em Y por (33), (34) e (36), logo (B+s)Q é compacto. D:=Q(Y) com norma de gráfico é fechado porque J^-1=Q é limitado. H é invertível para s≥414 (ALGEBRAIC_TAIL) |
| A1.8 | o gerador gauge em s=μ pertence a D+ | RADIUS_RECOVERY, ledger S4 | CONFIRMADO CONDICIONALMENTE | a recuperação vale para QUALQUER h∈ker L(μ)⊂D-. Que h=(1+μ^-1∂τ+ξ∂ξ)ω* esteja em D- e em ker L(μ) vem de GAUGE_RT_ACTION, que não auditei |
| A1.9 | a matriz exportada é o operador selecionado de RT | rt-A-4x12.dat, rt-A-6x18.dat | CONFIRMADO | reconstrução pelas fórmulas do TeX com RefA COMPLETO: 0 divergências em 19044 e 110889 entradas. Confirma também que o corte de RefA em 2G não afeta o bloco G |

### A2. Acoplamento guardado e Schur 0.997676

| # | Alegação | Arquivo | Veredito | Evidência |
|---|---|---|---|---|
| A2.1 | ‖Q_ref F‖_η<5.942 em F=6×18 | crown report, GUARDED | CONFIRMADO | valor EXATO: ‖Q_ref F‖=1/μ_RefA=536870912/90359175≈5.94152073654944 (máximo na coluna (0,0)). O relatório dá 5.94152073659…≥exato |
| A2.2 | ‖Q_ref‖ global<5.942 como max(colunas finitas, cauda) | GUARDED, `audit_sharp_free.py` | CONFIRMADO | a norma ℓ¹ ponderada é o supremo por COLUNA, logo max e não soma. Cauda forte 81/19. No sharp (pesos menores, m par): máximo exato também 1/μ, cauda 153/37 |
| A2.3 | ‖THF‖ para G=6×18,…,176×944 (59.947 … 1.7909e-5) | guarded-exterior.json | CONFIRMADO | as 5 linhas recomputadas por filtragem própria do suporte de RefA e majorante próprio, com igualdade racional exata. O filtro m≥M'-M+1 ou n≥N'-N+1 está correto: SΞΓ2 não divide (RΞ=1, (19a)) e P não altera \|n\| |
| A2.4 | ‖H_TT-I‖,‖H_ST‖≤d=(β+5/4)t(G)+114ε_μ; ‖H_TS‖≤c | GUARDED | CONFIRMADO | T·S=0; Q, D', SΓ1 e s preservam F. A parcela de μ não decai e está corretamente fora da cauda. 114=(1+18)/μ_min é conservador frente a K2≈107 de (37b) |
| A2.5 | defeito de Schur d+κcd<0.997676 uniformemente no contorno | guarded-exterior.json | CONFIRMADO | recomputado com constantes próprias (ε_μ=43·2^-40+2^-277 de RT §5 e (60), ‖D'Q‖≤K2 de (37b)): 0.9976754 com o κ do relatório e 0.9975897 com o κ independente (A3). κ admissível até ≈166 |
| A2.6 | todas as normas na mesma norma (η × pesos RT, soma por coluna) | scripts do A2/A3 | CONFIRMADO | lidas todas as reduções: `weighted_bound`, `norm_real`, `norm_complex` e `weighted_norm` somam por coluna. A ordem de η (916,4096,243,1233) coincide com o índice d dos endereços. Arredondamentos saem por teto racional |
| A2.7 | orçamento ‖H_true-H_ref‖<1.7918e-5 | GUARDED | CONFIRMADO | 114ε_μ+‖Q_ref‖(g1_η ε_μ+6(4096/243)ε_ω), com ‖Q_ref‖ exato (A2.2) |
| A2.8 | "o subsistema S+T é uniformemente invertível" | ledger/RETOMADA | CONFIRMADO, utilidade baixa (ver b2) | vale no CONTORNO, não no interior. ‖R_T^-1‖≤431 e a margem em d é 0.0023 |

### A3. Coroa amortecida: 51 tiles

| # | Alegação | Arquivo | Veredito | Evidência |
|---|---|---|---|---|
| A3.1 | 51 tiles, defeito uniforme <0.503, sup_contorno ‖H_SS^-1‖<36.227 com fundo | crown-contour-6x18-background.json | CONFIRMADO | recertifiquei os 51 com Q_ref EXATA (sem o resíduo 108), V própria (proposta float, grade 2^-56), resíduo inteiro exato e erro do fundo de RT. Todos aceitos; defeito máx. 0.48850 e κ≤31.4266, ambos no tile 1/8+i/4. Minhas cotas ficam entre 0.87 e 1.00 das do relatório |
| A3.2 | a validação não usa float | `certify_crown_tiles.py` | CONFIRMADO (lido) | float só nas propostas V e Q_hat. Produtos inteiros por limbs de 15 bits com checagem de overflow int64; normas em `int`/`Fraction` |
| A3.3 | variação no tile controlada | idem | CONFIRMADO | H(s)=H(s_j)+(s-s_j)Q_SS é exato (pencil afim), logo o defeito é ≤a+r‖V Q_SS‖ |
| A3.4 | cobertura do retângulo 1/8≤Re s≤1, \|Im s\|≤1/4 | idem, `check_crown_cover.py` | CONFIRMADO | checagem própria: segmentos contínuos, cada um no disco do seu tile, três lados completos. A metade inferior vem por conjugação, válida porque A, Q_ref e B são reais |
| A3.5 | ‖H_true-H_ref‖<0.0003257 | reference-error-budget.json | CONFIRMADO | derivação própria por RT (37a,b), 43·2^-40 de §5 e (60): <0.0003225 |
| A3.6 | versão sem fundo, defeito <0.499 | crown-contour-6x18.json | NÃO VERIFICADO | não refeita. Irrelevante, porque a versão com fundo foi recertificada |

### A4. Cadeia sharp

| # | Alegação | Arquivo | Veredito | Evidência |
|---|---|---|---|---|
| A4.1 | as matrizes sharp exportadas são K(0)=O♯+μΞ^-1_reg Γ1♯+Γ2♯(RefA,·) | rt-sharp-8x24.dat, rt-sharp-12x36.dat | CONFIRMADO | 0 divergências em 63504 e 352836 entradas, contra reconstrução pelo TeX (§3, (26)) |
| A4.2 | disco \|s-0.733\|≤1/32 em 8×24 com fundo: defeito <0.263453, inversa <11.447 | sharp-disk-8x24-background.json | CONFIRMADO | recomputado: 0.2634524 e 11.44590 |
| A4.3 | "‖K_true,N-K_ref,N‖<1.040e-7" | SHARP_EXTERIOR_BOUNDS | CONFIRMADO COM NOTA | vale com ε_μ lido de RefAplusB.dat (3.8811e-11). Só com o bound publicado 43·2^-40 dá 1.04050e-7, acima de 1.040e-7. Ambas as leituras são legítimas, mas a nota deve dizer qual ε_μ usa |
| A4.4 | β♯≤4.62237 | sharp-exterior.json | CONFIRMADO | majorante 2×2 refeito do TeX: igual entrada a entrada ao do produtor, β♯=4.6223694 |
| A4.5 | inversa finita precondicionada 12×36 <11.461 | sharp-disk-12x36-preconditioned.json | CONFIRMADO | recomputado: defeito à direita 0.2658490 e ‖J_ref K^-1‖≤11.4601113 |
| A4.6 | ponte F+T: δ<0.915031 (160×960) e 0.989223 (148×840) | sharp-guarded-bridge-*.json | CONFIRMADO | recomputado: 0.9150305 e 0.9892222 |
| A4.7 | graus de Neumann 437 e 3785 para erro <1e-12 | sharp-shell-remainder-*.json | CONFIRMADO (aritmética) | refeitos: 437 e 3785. Os fatores ℓ≈289 e 312, r≈22.42 estão corretos (lido e recalculado) |
| A4.8 | "a exclusão sharp foi reduzida a um Schur finito de 228366 DOFs" | ledger, SHARP_INTERMEDIATE_SHELL | SOBREVENDIDO | é uma redução formal a um problema computacionalmente inviável: ver (b3) e §A4 |
| A4.9 | D_min=D_max para o principal sharp; bounds 8385 e 1440 | SHARP_DOMAIN | CONFIRMADO (papel) | argumento de ODE por modo refeito: (1+ξ)^α e ξ^α com Re α<-1 não são analíticos. Os valores 8385=129·130/2 e 1440 conferem |
| A4.10 | C(s):D_L(a,b)→D_K(a',b') e KC=ML nos domínios | SHARP_DOMAIN | NÃO VERIFICADO | depende de SHARP_LINEARIZED_IDENTITY (Físico) |
| A4.11 | cadeias: K0 c_j+c_{j-1}=0 com c_j=C0h_j+C1h_{j-1} | SHARP_LOCAL_EXCLUSION | CONFIRMADO (papel) | comparação de coeficientes em K(s)C(s)h(s)=M(s)L(s)h(s) |

### A5. Rota algébrica, blocos e auditoria analítica

| # | Alegação | Arquivo | Veredito | Evidência e fonte RT |
|---|---|---|---|---|
| A5.1 | ‖(O_μ+z)^-1‖≤18/(μ+Re z) | ALGEBRAIC_TAIL | CONFIRMADO | RT (22), com z incorporado: O+z=BandI(Band_{0,D+z}+Band_{1,-im/2+μn-z}). \|f_J\|,\|f_K\|≤1 verificados, com diferenças μ(2n+1)(μ+2a) e 4μ(n+1)a. Constante de (35), (37a), (38): 2(1+r)/(1-r)=18 |
| A5.2 | ‖B‖<23 | ALGEBRAIC_TAIL, AlgebraicTail.lean | CONFIRMADO | (17/100)K1 com K1=80/9 de (36) e §5; 2L1=21 de (39a) e §5; 6(L2+R) de (34), (39b) e (60). μ_true≤17/100 e ≥1/6 seguem de §5 e (60) |
| A5.3 | base 828 e sem espectro em Re s≥414 na realização | ALGEBRAIC_TAIL | CONFIRMADO | aritmética compilada em Lean. L(s)=(J+s)(I+Q(s)B) com 23·18/(414+1/6)=2484/2485 |
| A5.4 | piso de transferência 2592/11 e 2592/35; 865×2594 | ALGEBRAIC_TAIL | CONFIRMADO | Lean compilado; conferido à mão |
| A5.5 | det(I+ΠD Q0 ι)=det(A_N+sI)/det(J_N+z0) | ALGEBRAIC_TAIL | CONFIRMADO (papel) | Q0ι=ιQ0_N porque J preserva caixas. Sylvester dá Q0_N(A_N+s) |
| A5.6 | β<10.382 e TT<4.882 em 64×192; 117774/117875 em 175×942 | COMPONENT_EXTERIOR | CONFIRMADO | A1.4 e cauda A1.3: 11.632×81/193 |
| A5.7 | tabela de STRUCTURED_EXTERIOR (7857/292, 7857/772, 7857/7860) | STRUCTURED_EXTERIOR | CONFIRMADO | (23+5/4)×cauda |
| A5.8 | critério 2×2 ⇔ bc<(1-a)(1-d); Schur a+bc/(1-d); homotopia a+t²bc/(1-td) | FREQUENCY_BLOCKS §2, UNIFORM_COUPLED | CONFIRMADO (papel + Lean compilado) | refeito |
| A5.9 | tabelas §§4–5 de FREQUENCY_BLOCKS | FREQUENCY_BLOCKS | NÃO VERIFICADO | superadas pelos bounds de componentes. Custo: ~1 h |
| A5.10 | contraexemplo: comprimir a inversa ≠ inverter a compressão (2/3 contra 1/2) | ANALYTIC_AUDIT, ResolventCompression.lean | CONFIRMADO | Lean compilado. M_N-R_N=-ΠR0E_N rederivado |
| A5.11 | obstrução a P2 no modelo 1+∂τ em ℓ¹ | ANALYTIC_AUDIT | CONFIRMADO, diz exatamente o que afirma | ‖Rf_n‖_X≤√2 enquanto ‖Rf_n‖_σ=e^{σn}/√(1+n²)→∞. Refuta a INFERÊNCIA "fundo analítico ⇒ P2", não um teorema sobre Choptuik |
| A5.12 | os três módulos Lean verificam o que dizem | AlgebraicTail, FrequencyBlocks, ResolventCompression | CONFIRMADO | só aritmética escalar e matrizes 2×2. Nenhum teorema de operadores; os próprios cabeçalhos dizem isso |
| A5.13 | linhas C1 ("composição fechada e ancorada") e C2 ("só tailError é hipótese") do ledger | PROOF_OBLIGATIONS | SOBREVENDIDO | contradizem ANALYTIC_AUDIT (C1-G e falha de P2). Ver (b4) |

## Detalhes das rederivações

### A1: cauda livre e recuperação de raio

**Coluna exata.** Por (22), no bloco J: (O_μ c)_j=(im/2+μ(j+1))c_j+Σ_{i>j}2μ i c_i.
No bloco K a soma vale só para i≡j mod 2. Para Q e_n:

- c_n=1/D_n;
- c_{n-k}=-2μn/(D_n D_{n-k}), com k=1 para J e k=2 para K;
- c_j=(-f_j)…(-f_{n-k-1})·c_{n-k}, com \|f\|≤1.

Os pesos satisfazem w_j/w_n≤r^{n-j}, logo a coluna ponderada é
≤√2[1/\|D_n\|+\|c_{n-k}\| r^k/(1-r^k)], com √2≤3/2.

- **Faixa \|m\|≥M, b=m/2.** 2μn/(\|D_n\|\|D_{n-1}\|)≤2x/(b²+x²)≤1/b. Para K,
  2μn≤2(x+μ_max) com x=μ(n-1). Isso dá exatamente os dois termos de Fourier.
- **Faixa n≥N.** \|D_j\|≥μ(j+1) dá (1+2r/(1-r))/(μ(n+1)) para J e
  (1+2n r²/((n-1)(1-r²)))/(μ(n+1)) para K. O segundo é decrescente em n.

As duas faixas cobrem o exterior. `independent_free_inverse.py` verifica que
as colunas exatas amostradas ficam abaixo do bound.

**Lógica da recuperação.** Seja h∈ker L(s)∩D- e x=Jh∈Y-. Então Hx=0 e
T H T(Tx)=-T H G(Gx). Gx tem suporte finito, logo Gx∈Y+.
y:=-(H_TT,+)^-1 H_TG Gx está em Y+⊂Y- e resolve a mesma equação em Y-.
Pela unicidade em Y- (contração 0.786), Tx=y. Logo x∈Y+ e h=Qx∈D+.

Para cadeias L h_j=-h_{j-1}, o termo -T h_{j-1} é forte por indução e o
argumento se repete. Com as cadeias coincidindo, as multiplicidades
algébricas de raízes da família Fredholm coincidem.

**Fredholm.** B é limitado em Y. Q é limite em norma de QG. (B+s)Q é
compacto. H(s)=I+compacto é afim em s e invertível para s≥414, logo o
teorema analítico de Fredholm se aplica em ℂ.

Observação: para concluir Fredholm, basta t(G)→0. A caixa 1024×4096 e o
β por componentes só são necessários para a recuperação de raio.

**Robustez.** Trocando β_η pelo bound SÓ de RT (sem pesos de componentes)
μ_max K1+2L1+6(L2+R), os defeitos ficam 0.470 no forte e 0.938 no fraco.
A recuperação sobrevive sem confiar no majorante novo.

### A2 e A3: o que o Schur S+T controla

A eliminação de S dentro do bloco comprimido Z=S+T está correta. O
resultado é invertibilidade de H_ZZ NO CONTORNO. Por si só, isso não diz
nada sobre H inteiro.

O "problema restante" é o complemento de Schur em P+U, com
U=(176×944)-(6×18) contendo 578811 DOFs reais. Seus acoplamentos seriam
multiplicados por ‖H_ZZ^-1‖, e só a parte T já custa 431. Não é um
problema tratável.

### A4: a rota da faixa U é inviável

Para 160×960:

- **Dimensão.** E_U tem 228366 DOFs. A matriz densa em float64 ocupa ≈417
  GB. A máquina tem 7 GB; em racionais é impossível.
- **Construção de P_N.** P_N=Σ_{k≤437}(-D)^k age no exterior INFINITO do
  fundo exato. Só a parte RefA teria suporte finito, e ele cresce ~100
  índices radiais por potência. A bola ε não tem suporte finito e entra
  apenas como majorante.
- **Amplificação.** Os majorantes certificados de L_U R_T^-1 R_U somam
  ℓ·r/(1-δ)≈289·22.4/0.085≈7.6·10⁴. Nem a proximidade E_U≈H_UU está
  certificada. Qualquer parametriz de E_U teria de absorver essa
  amplificação.

Veredito: **beco sem saída** com os majorantes atuais. O grau de Neumann
não é o gargalo; o tamanho da caixa G é. G é imposto pelo majorante
uniforme d(G)=(β+\|s\|)t(G), que só fica <1 com N≳940 (§c).

## (b) Texto proposto para o ledger

**(b1) F1**, substituindo a célula de evidência:

> Demonstrada em papel na realização D:=J_μ^-1(Y), com Y o ℓ¹ ponderado
> de RT (κ1,κ2)=(65/64,5/4), complexificado. Q=J_μ^-1 é compacta
> (‖Q(I-G)‖≤t(G)→0), B é limitado por RT (33), (34) e (36), e
> H(s)=I+(B+s)Q é Fredholm de índice 0, invertível para s≥414. Hipóteses:
> 1/6≤μ≤17/100 e a bola RT (39b) e (60). RADIUS_RECOVERY acrescenta: em
> \|s\|≤5/4, kernels e cadeias coincidem nos pesos (65/64,5/4) e
> (129/128,9/8), recomputado de forma independente em 16/09. Não
> formalizado em Lean. Não identifica D com o domínio geométrico.

**(b2) GUARDED**, na frase do ledger:

> O bloco comprimido a Z=S+T (S=(6×18)-(4×12), T fora de 176×944) é
> invertível NO CONTORNO, com defeito de Schur <0.997676. Os valores foram
> recomputados de forma independente em 16/09: 0.99759 com κ<31.43. Isso
> não é uma afirmação sobre H. O complemento em P+U (U com ~5.8·10⁵ DOFs)
> é intratável, e ‖(H_ZZ)^-1‖ só é limitado por ~431·κ. Serve como
> verificação de consistência, não como passo de uma prova.

**(b3) SHARP_INTERMEDIATE_SHELL**, substituindo "A exclusão sharp foi
reduzida a um Schur finito…":

> F=12×36 mais T fora de 160×960 é invertível em todo o disco
> \|s-0.733\|≤1/32, condicionado à bola RT (δ<0.915031, recomputado). A
> exclusão sharp infinita EQUIVALE formalmente à invertibilidade de um
> Schur de 228366 DOFs cujos fatores contêm o exterior infinito do fundo
> exato. Na máquina disponível isso é inviável (matriz densa ≈417 GB;
> majorante de correção ≈7.6·10⁴). Os graus 437 e 3785 são orçamentos
> teóricos, não um plano computacional. A rota exige primeiro reduzir G
> (ver próximo passo).

**(b4) C1 e C2**, colunas de estado:

> C1: composição aritmética em Lean (TailBound, TailTileBridge) correta, mas
> as hipóteses P1–P3 NÃO são o caminho. P2 falha no modelo 1+∂τ
> (ANALYTIC_AUDIT), e a ponte compressão/inversa finita (C1-G) está aberta.
> Substituto em curso: cauda algébrica com pesos RT (ALGEBRAIC_TAIL e
> STRUCTURED_EXTERIOR), cujo majorante uniforme exige caixas ≳176×944.

> C2: além de `tailError`, são hipóteses externas ao kernel: C1-G, o erro
> do fundo, a cobertura do contorno e a identificação dos blocos com normas
> de operadores infinitos.

**(b5) A4.3**: acrescentar "com ε_μ=\|μ_RefAplusB-μ_RefA\|+2^-277 lido do
arquivo; com o bound publicado 43·2^-40 o valor é <1.041e-7".

## (c) Próximo passo de maior alavancagem

### Diagnóstico que motiva o passo

Isto é ponto flutuante e NÃO é prova. Usei
`independent_column_diagnostic.py`, validado contra cálculo exato em
`Fraction` em 4 colunas, com diferença relativa ≤1.5e-6.

Todas as rotas auditadas consomem o mesmo número: o majorante exterior
uniforme d(G)=(β+\|s\|)·t(G). Isso inclui a recuperação de raio, S+T, a
ponte sharp, as faixas U e ε_N do objetivo 1. Esse majorante só fica <1
com N≳940, o que impõe caixas de ~5.8·10⁵ DOFs.

Medindo as normas REAIS de coluna de (B+s)Q no setor A (pesos η, s=1 e
s=1/8+i/4):

| n | máx. real amostrado | majorante uniforme | razão |
|---:|---:|---:|---:|
| 36 | 4.69 | 24.9 | 5.3 |
| 144 | 1.25 | 6.36 | 5.0 |
| 288 | 0.63 | 3.19 | 5.0 |
| 576 | 0.32 | 1.60 | 5.0 |

As normas reais ficam abaixo de 1 a partir de n≈190. Na direção de Fourier
já dão 0.44 em m=48. Há ainda um majorante intermediário "por grau", sem
cancelamento entre graus, calculável de forma semianalítica. Ele fica ~3×
abaixo do uniforme (2.09 contra 6.31 em n=145), o que moveria o corte para
N≈300.

Uma caixa 24×192 tem 15648 DOFs; 176×944 tem 579144. São ~37× menos DOFs.

### Tarefa proposta

**"Majorante exterior por colunas d*(G) para T(H(s)-I), setor A e sharp"**
(análogo espectral de EXTERIOR_ESTIMATE, RT (52a) e (53)).

**Enunciado.** Construir um racional d*(G) tal que, para toda coluna e
fora de G (m≥M_G ou n≥N_G) e todo \|s\|≤5/4:

- ‖(H_true(s)-I)e‖_η≤d*(G)‖e‖_η, com erro de μ e bola ε_ω incluídos;
- idem para o sistema sharp nos pesos (129/128,9/8), no disco \|s-0.733\|≤1/32.

**Fase 1 (semianalítica, barata).**

1. Normas exatas de coluna de B por grau, ‖B e_{d,m,j}‖_η, calculadas em
   `Fraction` para j<J*, m<M*.
2. Um bound exato uniforme para j≥J*, m≥M*. Acima do suporte 41×101 de
   RefA não há dobra de índices, e a norma vira uma soma finita explícita
   que não depende de (m,j).
3. A estrutura exata das colunas de Q: c_n, c_{n-k} e produtos de -f.

Isso dá d_pd(G)=sup_e Σ_j\|c_j\|w_j/w_n·β_col(j)+\|s\|‖Qe‖+erros.

**Fase 2 (só se a Fase 1 não chegar a N_G≤384).** Normas de coluna exatas
na faixa, no estilo de RT, com aritmética diádica arredondada para fora.

**Arquivos:**
- `scripts/columnwise_exterior.py`
- `scripts/test_columnwise_exterior.py`
- `docs/COLUMNWISE_EXTERIOR.md`
- relatório `build/spectrum/columnwise-exterior.json`

**Critérios de aceitação:**

1. d* ≥ norma exata em `Fraction` em ≥50 colunas exteriores amostradas,
   incluindo m=0, m ímpar e K e J (teste).
2. d* ≤ 2× o diagnóstico float nessas colunas, para medir a folga.
3. Relatório da menor caixa com d*<1, uniforme em s e com a bola RT.
4. Reexecução do Schur de GUARDED e da ponte sharp com essa caixa,
   registrando as novas margens.
5. Regra de decisão registrada ANTES de rodar:
   - se N_G≤384 no setor A e no sharp, a rota de blocos continua, e o
     próximo gargalo passa a ser uma parametriz certificada de ~2–3·10⁴
     DOFs (ainda pesado, mas finito na máquina);
   - se N_G>384, a rota de caixas é declarada inviável neste hardware e o
     trabalho passa a uma parametriz exterior analítica que absorva B
     (características).

**Honestidade.** Mesmo no melhor caso, este passo não dá winding nem
exclusão. Ele mede se a transferência finita→infinita é computacionalmente
alcançável sem HPC. É a pergunta de que dependem todos os objetivos 1 e 3,
e hoje ela só tem resposta por majorantes 5× frouxos.

## Reprodução

```bash
cd scripts
../.venv/bin/python independent_free_inverse.py
../.venv/bin/python independent_component_beta.py
../.venv/bin/python independent_selected_matrix.py ../build/spectrum/rt-A-6x18.dat
../.venv/bin/python independent_sharp_matrix.py ../build/spectrum/rt-sharp-12x36.dat
OPENBLAS_NUM_THREADS=1 ../.venv/bin/python independent_crown_tiles.py --workers 2 --out crown.json
../.venv/bin/python independent_sharp_disk.py ../build/spectrum/rt-sharp-8x24.dat --mode left
../.venv/bin/python independent_sharp_disk.py ../build/spectrum/rt-sharp-12x36.dat --mode right --out sharp.json
../.venv/bin/python independent_schur_scalars.py --crown-json crown.json --sharp-json sharp.json
../.venv/bin/python independent_column_diagnostic.py        # DIAGNÓSTICO float
cd .. && .venv/bin/python -m unittest discover -s scripts -p 'test_independent_audit.py'
```

## (d) Revisão cruzada da auditoria física (rodada 2)

Li só as §§1, 2, 7 e 9 de `AUDITORIA_FISICA_16SET.md`. O item (iii), a raiz
da faixa B confirmada pelo orquestrador, não foi refeito.

### (i) Em que norma o testemunho de 0,401 tem de valer, e se a queda ameaça

Diagnóstico em float, NÃO é prova:
`scripts/independent_witness_norm_diagnostic.py`, com parsers próprios e
raiz real de L_A mais próxima de 0,401. Saída em
`scratchpad/matematico/witness_norm_0401.json`.

| caixa | ‖Ch‖/‖h‖ l2 | fraco/fraco | fraco/fraco na caixa fixa 4×12 | ‖h‖_fraco/‖h‖_2 | fração de ‖h‖_fraco na borda | ‖Ch‖_fraco/‖Ch‖_2 | Ch fraco / grafo forte (sem η) | Ch fraco / ‖Jh‖ forte·η |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 6×18 | 0,3195 | 0,1149 | 0,0881 | 33,4 | 22% | 12,0 | 0,0296 | 4,7e-5 |
| 8×24 | 0,3250 | 0,0742 | 0,0694 | 35,5 | 7,8% | 8,1 | 0,0168 | 2,4e-5 |
| 10×30 | 0,3243 | 0,0600 | 0,0683 | 36,2 | 2,6% | 6,7 | 0,0127 | 1,8e-5 |
| 12×36 | 0,3243 | 0,0559 | 0,0683 | 36,5 | 0,8% | 6,3 | 0,0115 | 1,6e-5 |

"Borda" é m≥M-2 ou n≥N-6. A coluna fraco/fraco reproduz a sequência do
Físico.

**Leitura.** A queda NÃO é C h→0.
- ‖h‖ ponderado converge: 33→36,5, com a fração na borda caindo 3× por
  passo.
- A razão restrita a uma caixa fixa (4×12) já estabilizou em 0,0683.
- O que cai é ‖Ch‖_fraco/‖Ch‖_2, de 12,0 para 6,3. É contaminação de
  modos altos de Ch pelo corte, desaparecendo com o refinamento.
- A extrapolação de Aitken da razão total dá ≈0,054>0.

Veredito: a queda não ameaça a EXISTÊNCIA do testemunho.

**Mas a norma decide a margem.** O testemunho só é compatível com a cadeia
que auditei se:
- h for medido na norma do domínio em que o autovetor for enclausurado;
- C(s)h for medido nos pesos sharp ESTRITAMENTE menores (perda de raio,
  SHARP_DOMAIN);
- c_op=‖C(s)‖ entre essas duas normas for majorado.

Na norma forte·η dos tiles, de GUARDED e de d*_A, com x=Jh, c_ref≈1,6e-5.
Nessa norma a cauda de x em 12×36 ainda carrega 6% da massa, então o
enclausuramento exigido (ε_v≲c_ref/c_op) está fora de alcance. Rotas
viáveis:
- (a) h enclausurado nos pesos fortes SEM η (norma equivalente; o fator
  4096/243 entra uma única vez na conversão de ε_v) e Ch nos pesos
  (129/128,9/8), com c_ref≈0,011;
- (b) h nos pesos fracos, com c_ref≈0,054, e Ch em pesos ainda menores,
  por exemplo (257/256,17/16). Isso exige tiles e d* do operador
  selecionado nos pesos fracos, com β-≤(162/85)β+, mais caro.

Em ambos os casos, o d*_A em norma η da tarefa (c) NÃO serve diretamente
ao testemunho. Isso fica registrado em `COLUMNWISE_EXTERIOR.md`.

### (ii) "multiplicidade 1 + lower>0 ⇒ sai da contagem física" é suficiente?

É suficiente para "a raiz contribui 0 para dim E_c", ou seja, para as
cadeias que satisfazem as constraints na realização RT. Com dim E=1 não
há cadeias: D h=C(s0)h≠0 dá rank(D|E)=1 e dim E_c=0 pela fórmula de
CONSTRAINT_ROOT_QUOTIENT. Faltam hipóteses a nomear:

- **H1.** A multiplicidade 1 é da família Fredholm INFINITA H(s)=I+(B+s)Q,
  na mesma realização em que se conta, certificada por winding 1 num disco
  com cauda e bola RT. Multiplicidade 1 da matriz de Galerkin não basta.
  Pelo meu A1, a realização forte e a fraca dão o mesmo espaço de raízes em
  |s|≤5/4, logo a escolha dos pesos é livre DENTRO dessa faixa.
- **H2.** O enclausuramento é de um vetor de ker L_true(s0) na raiz
  verdadeira s0 (desconhecida, dentro do disco): um par autovalor/vetor
  do operador infinito (Newton–Kantorovich ou Krawczyk com normalização
  e cauda), não o autovetor de Galerkin. Além disso ‖v‖>ε_v, para
  excluir h=0.
- **H3.** C_true(s):D→Y♯(a',b') é limitado com constante explícita c_op
  (perda de raio). SHARP_DOMAIN afirma isso a partir das fórmulas de
  SHARP_LINEARIZED_IDENTITY; nenhum de nós verificou a extensão aos
  domínios (A4.10 e B2.5 NÃO VERIFICADOS).
- **H4 (falta no texto do Físico).** Para dizer "contagem FÍSICA", é
  preciso a identificação entre a multiplicidade física e
  dim(E∩ker D) módulo gauge, com a escolha do cone. Isso é S4 e a
  equivalência geométrica global, ainda abertas. Sem H4, a conclusão
  correta é "sai da contagem das raízes que satisfazem as constraints na
  realização RT".

A subtração da linha de gauge só atua em s=μ e não interfere em 0,401.

### (iii) Consequência para as minhas contas

O candidato B em Re≈0,041 e o retângulo B até Re=1, Im=3/4 cabem em
|s|≤5/4, que é o raio usado em A1, A2 e na tarefa (c). Nenhum bound
auditado precisa ser refeito por isso.
