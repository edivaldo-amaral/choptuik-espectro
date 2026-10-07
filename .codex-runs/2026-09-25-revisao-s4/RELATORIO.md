# Revisão independente de S4: relatório (25–26/09/2026)

Revisor: agente independente. Objetivo: tentar quebrar o certificado. Todas as contas usaram
`.venv/bin/python` com `CHOPTUIK_EPS_FUNDO_L=4.7e-08`. As contas pesadas rodaram em `~/revisao_s4/`
no PC do laboratório (R3) e no laboratorio2 (R4-janela), sobre cópias de `scripts/`, de `RefA.dat` e dos
JSON/npy de `build/s4/rig`. Nenhum arquivo existente foi alterado: os sha256 de 31 originais antes e
depois são idênticos (`r7_sha256_antes.txt` = `r7_sha256_depois.txt`), e o `git status` não mostra
arquivos rastreados modificados. Não abri nenhum dos arquivos proibidos.

## Tabela de veredictos

| Item | Veredicto | Artefato | O que foi feito |
|---|---|---|---|
| R1 (B1: fundo, RefB, μ) | **OK** | `r1_refB.py`, `r1_saida.txt` | Enunciado da bola de RT localizado no artigo. Leitor próprio dos dois `.dat`, RefB exato em inteiros. Young com tabela própria transcrita do artigo. Conversão de Corr. MU_80 e DELTA_MU. Grep de todos os usos de `rt`/`EPS_FUNDO_L`/`eps_mu`. |
| R2 (B2: gerador de gauge) | **OK** (identidade e convenção); o domínio de g* entra em L2/L3 | `r2_identidade.py`, `r2_saida.txt`, `r2b_constraints.py`, `r2b_saida.txt` | Identidade L_w(μ)(Z+1)w = (Z+2)F_w testada nos operadores do código, com F_w e Z montados por mim e controles que falham. Mesmo teste para C_w(μ)(Z+1)w = (Z+1)c_w (dá D g* = 0). Resíduo de g_A na saída completa. |
| R3 (etapa A em μ) | **OK** | `r3_etapaA_remoto.py`, `r3_saida.txt`, `r3_resultado.json` | M̂ bordejado (14 017) montado por mim no laboratório. ‖M̂⁻¹‖, Y_Z e o acoplamento "a" calculados em ponto flutuante (ARPACK), sobre S = YYᵀ com 52 184 colunas de T. |
| R4 (etapas B e C em μ) | **OK** | `r4_json.py`, `r4_saida.txt`, `r4_janela_remoto.py`, `r4_janela_saida.txt`, `r4_janela_resultado.json`, `rig_copia/Q.json` | Q.json recalculado: idêntico bit a bit. λ0, rt, eps_mu, beta, pert e Qfar recalculados em μ_A. Normas verdadeiras de 4 blocos de E na janela 32..47 (no laboratorio2). Junção dos checkpoints e origem na mesma estimativa. |
| R5 (B4: identificação) | **OK com ressalva** (L2, L3, L4) | `r5_etapaE.py`, `r5_saida.txt` | Identidade de y_g − y0 e cotas derivadas. nz, cauda, h0, ‖∂τ RefA‖, e, ‖L_A(λ0)h0‖, δ e dist recalculados por conta própria. Checado ‖DT(x)‖ ≤ θ + 2 Z2 ‖x − x0‖ a partir das definições. Domínio de y_g. |
| R6 (B5: contagem) | **LACUNA** (L1, L5; hipóteses adicionais listadas) | `r6_alias.py`, `r6_saida.txt` e a análise abaixo | Lógica da contagem e os cinco documentos de gauge. Simetria de aliases de Floquet verificada no código. |
| R7 (reprodução) | **OK** | `r7_log.txt`, `rig_copia7/`, `r7_sha256_*.txt` | `nk_L3_C.py` e `nk_L3_E.py` rodados numa cópia: C3.json e E.json saíram idênticos aos originais (igualdade de JSON). |

## Resultados por item

### R1: fundo verdadeiro e RefB

- **A bola de RT** (choptuik.tex, eqs. `ref5ref5ref5ref5`, `corrinb`, `kdkkkskkskksksks`): a referência
  é (μ, ω)_ref = (RefA) + (RefB). Com 𝓡 = 2⁻²⁷⁷, o ponto fixo satisfaz
  |μ − (μ_A + μ_B)| ≤ 2⁻²⁷⁷ e ‖ω − (RefA + RefB)‖_SSPACE ≤ 2⁻²⁷⁷. A norma é
  ‖v‖_𝒱 = Σ_{m,n∈ℤ} κ1^{|m|} κ2^{|n|} (|Re| + |Im|), somada nas 4 componentes (é mais forte que "por
  componente"). O centro é portanto RefA + RefB, e o `choptuik.c` (`LoadAplusB`) confirma que
  `RefAplusB.dat` é a soma: refB = refAplusB − refA.
- **Leitor próprio.** Grupos aninhados explícitos, `num_m` grupos de `num_n` pares. RefA tem
  2⁻⁷⁰ e 41×101. RefA+RefB tem 2⁻⁶⁵⁰ e 450×1350. As paridades dos setores conferem: v1 só com n ímpar,
  v1 a v3 só com m par, v4 só com m ímpar. Com isso ficam confirmados o layout (m, n) e o
  índice interno n.
- **Resultados:**
  - ‖RefB‖_L = 2,350632·10⁻¹⁰ (autor 2,3506·10⁻¹⁰).
  - ‖(Z_A+1)RefB‖_L = 2,454216·10⁻⁸ (autor 2,4542·10⁻⁸). Minha fórmula de ξ∂ξ foi derivada e
    testada contra `chebder`, com erro 1e-13.
  - Young com a tabela transcrita do artigo (idêntica a `signed_operator_L.TABELA`):
    ‖B_L(RefB)‖ ≤ **4,670165·10⁻⁸** (autor 4,6748·10⁻⁸ = 4,6702·10⁻⁸ × 1,001).
  - Quase toda a norma de RefB está fora da caixa 41×101 de RefA (2,33·10⁻¹⁰ de 2,35·10⁻¹⁰). Os maiores
    termos de Young estão em d = 41 e 42.
- **Corr.** ‖B_L(Corr)‖ ≤ max_j W_j · 2⁻²⁷⁷ = 7,0·10⁻⁸³ < 10⁻⁶⁰, com W = (5, 17, 3, 12). Esses pesos
  são a soma de 2·fator·|c| da tabela. A conversão é válida porque os pesos de L com sinal,
  κ1^{m} ≤ κ1^{|m|}, e l1(a) = |a0| + 2Σκ2ⁿ|a_n| coincidem com a norma de RT quando se usa |z| ≤ ‖z‖_ℂ.
- **μ.** |MU_80 − μ_{A+B}| = 2,16·10⁻⁸¹, logo |μ* − MU_80| ≤ 2,2·10⁻⁸¹ ≤ 10⁻⁸⁰.
  |μ_A − μ*| ≤ 3,881119·10⁻¹¹ ≤ DELTA_MU = 3,881507·10⁻¹¹.
- **Onde o erro entra** (grep em `scripts/`). `rt` entra por `fundo_L.dB_L` e `eps_mu` por
  `eps_mu_L(λ0)`:
  - etapa A: dM, dr, aZ, epsZ;
  - etapa B: beta, desc, pert, resíduo das janelas;
  - etapa C: eg em todas as entradas, e rtY;
  - etapa E: rL e nB.

  Não encontrei constante de S3b esquecida. O que depende de s foi recalculado em μ_A e bate bit a bit
  (ver R4). O único ponto frágil é que o valor 4,7·10⁻⁸ só existe pela variável de ambiente (L5).

### R2: gerador de gauge

Montagem própria: F(w) = J(0)w + μSΓ1w + SΞΓ2(w,w) = J(0)w + (B_w + B_0)w/2, e Z = μ⁻¹∂τ + ξ∂ξ.

- (T0) Simetria de SΞΓ2: erro relativo 2,3·10⁻¹⁶.
- (T1) w polinomial aleatório, caixa exata, μ = 0,37:
  - com a convenção do autor: ‖L_w(μ)(Z+1)w − (Z+2)F_w‖/‖·‖ = **1,4·10⁻¹⁶**;
  - controle ∂τ → −i m/2: 1,5·10⁻²;
  - controle (Z+1)F no lado direito: 4,1·10⁻².

  O teste distingue as convenções: sinal e escala de ∂τ → +i m/2 são os mesmos de `J_livre`
  (σ = s + i m/2).
- (T2) RefA na saída completa, sem truncar:
  - F(RefA) = 6,3·10⁻¹¹: não há termo constante, e o F do código é o de RT;
  - ‖L_A(μ_A)g_A‖ = 6,9017·10⁻¹⁰ = ‖(Z_A+2)F(RefA)‖, com diferença de 1,2·10⁻¹⁵;
  - ‖L_A g_A‖/‖g_A‖ = 4,42·10⁻¹⁰ (autor 4,4·10⁻¹⁰).
- (T4) g_A de `gauge_vector_L` coincide com o meu (diferença 1,2·10⁻¹⁶).
- (T5) Constraints: ‖C_w(μ)(Z+1)w − (Z+1)c_w‖ relativo = 2,2·10⁻¹⁶, e os controles falham (10⁻²).
  c(RefA) = 1,6·10⁻¹¹. Logo D g* = C(μ)g* = (Z+1)c(w*) = 0, dado que a solução exata de RT satisfaz
  as constraints, o que é parte do teorema de RT.
- **Observação estrutural.** Nas componentes 1 a 3, a parte livre é J = ∂τ + μ((1+ξ)∂ξ + 1) (todas as
  linhas acima da diagonal). Na componente 0 é J = ∂τ + μ(ξ∂ξ + 1). Portanto μ(Z+1) ≠ J(0): a diferença
  é −μ∂ξ nas componentes 1 a 3. A identidade vale mesmo assim, como mostrado acima. Isso importa em L2.

### R3: etapa A em μ

Estimativas em ponto flutuante contra as cotas do autor:

| Grandeza | Minha estimativa | Cota do autor |
|---|---:|---:|
| ‖M̂⁻¹‖ (3 maiores valores singulares: 52,040; 20,31; 14,32) | 52,04043 | 52,14972 |
| Y_Z sem margens, ‖M̂⁻¹F_Z'(x0)‖ | 5,098·10⁻⁷ | Y_Z = 2,996·10⁻⁶ (as margens vêm de rt·tV) |
| a, parte de Loewner (colunas de T com \|m\| < 48 e banda completa) | 0,395262 | a_loewner = 0,404 (aZ = 0,40461) |

Minha estimativa de "a" coincide com a estimativa do autor em `nk3_mu.json` (0,395262). O valor
grande de ‖M̂⁻¹‖ reflete a raiz de fase s = 0, a 0,168 de λ0. Não se trata de uma subestimação.

### R4: etapas B e C em μ

- **Q.json** recalculado com `nk_L3_q.py --dir rig_copia` (λ0 lido de Z.json = μ_A): **idêntico** ao
  original.
- **Parâmetros em μ.** Z.json e C3.json têm lam0 = μ_A exatamente. Em Z.json, T3.json e E.json vale
  rt = 4,7·10⁻⁸. O eps_mu de Z e de T3 é igual a eps_mu_L(μ_A) = 2,82079·10⁻⁸; em 0,401 seria
  3,156·10⁻⁸. dB_L, beta e pert de T3 e Qfar de Z foram recalculados bit a bit em μ_A (em 0,401,
  Qfar seria 0,124892).
- **Normas verdadeiras na janela 32..47** (montagem própria, laboratorio2):

  | Bloco | Verdadeira | Cota de T3 |
  |---|---:|---:|
  | ‖V_J‖ | 1,74117 | 1,74134 |
  | T←Z' (banda 40) | 0,38047 | 0,39412 |
  | ← 48..63 | 0,25844 | 0,25919 |
  | ← 16..31 | 0,26343 | 0,26538 |

  Todas as cotas estão acima. Em S3b (λ0 = 0,401), ‖V‖ dessa janela era 1,6918, o que confirma que
  T3 foi calculado em μ.
- **Junção dos checkpoints.** Os checkpoints têm 11 e 15 janelas, com interseção vazia e união igual
  às 26 janelas. T3.json é a junção exata, com N recebendo ρ na diagonal. Nas 26 janelas, nV = est ×
  1,0001 × 1,002^k × (1 + 10⁻¹²) com k inteiro, e as outras entradas ficam ≥ à estimativa: as duas
  máquinas partiram da mesma estimativa `nk3_mu.json`.
- **Assimetria de Y** (6,1·10⁻⁶ em 32..47 contra 2,3·10⁻⁶ em −47..−32). É esperada: os pesos com sinal
  κ1^{2m} não são simétricos em m, e a razão κ1^{2m} ≈ 3 para m ≈ 36 explica o valor. Não é
  inconsistência entre máquinas.

### R5: identificação (etapa E)

- **Identidade** (derivada): y = (J_A + λ0)h; L_*(μ*)g_n = 0; J_* = J_A + δμJ1; B_* inclui μ*SΓ1.
  Isso dá y_g − y0 = −B_*(g_n − h0) − L_*(λ0)h0 + (λ0 − μ*)g_n − δμJ1(g_n − h0). Como J1(g_n − h0) =
  J1Q0(y_g − y0), a divisão por 1 − |δμ|‖J1Q0‖ está correta. Com g_A = nz h0 + t (t ⊥ Z) e ⟨h0, g*⟩ ∈
  [nz − ‖e‖, nz + ‖e‖], vale ‖g_n − h0‖ ≤ (‖t‖ + 2‖e‖)/(nz − ‖e‖). O e de ‖g* − g_A‖ decompõe-se
  corretamente em (Z_A+1)RefB + Δ(1/μ)∂τ(RefA + RefB) + (Z*+1)Corr. O fator 1000‖RefB‖ para ‖∂τRefB‖
  é válido porque m ≤ 449.
- **Recalculado por mim:**

  | Grandeza | Meu valor | E.json |
  |---|---:|---:|
  | nz | 2,6043546204 | idem |
  | cauda | 6,856173·10⁻⁷ | com as folgas do autor |
  | ‖h0w − P_Z g_A/nz‖ | 2,4·10⁻¹⁶ | — |
  | ‖∂τRefA‖_L | 0,373452 | — |
  | e | 2,50539·10⁻⁸ | 2,50562·10⁻⁸ (o autor multiplica ZB por 1,0001) |
  | ‖g_n − h0‖ | 2,82498·10⁻⁷ | — |
  | ‖L_A(λ0)h0‖ | **3,708878·10⁻⁶** | idem |
  | δ | 4,876413·10⁻⁶ | — |
  | dist | 1,815528·10⁻⁵ | — |
  | θ + 2 Z2 max(dist, r) | 0,902443 | — |

  A unicidade vale até ρ = 3,26·10⁻⁴, isto é, 18 vezes dist.
- **Normas.** O componente λ fica em Z' sem escala: a coluna bordejada é h0w e ⟨h0w, h0w⟩ = 1. Como w ≡ 1,
  ‖x_g − x0‖_v ≤ √(δ² + dλ²)/min v está correto, com min v = v_far = 0,2686.
- **Uniqueness.** F é quadrática, com DF(x) − DF(x0) = [(λ − λ0)Q0 dy + dλ Q0(y − y0); 0], e Q0 é
  diagonal em m e triangular superior em n. Refiz a cota por blocos: ‖A(DF(x) − DF(x0))‖_v ≤
  2v0 max(tV(nQZ v0 + qZT vT + qZf vf)/v0, maxwV(qTT vT + Qf vf)/vT, Qf)·‖x − x0‖_v, que é exatamente
  o Z2 de `nk_L3_C.py`. Então Z2 é a constante de Lipschitz inteira de A·DF. O teste NK Y + θr +
  Z2r² ≤ r e o critério θ + 2Z2ρ < 1 são conservadores (bastaria ½Z2r² e θ + Z2ρ). O critério vale
  para todo x.
- **Domínio.** Ver L2 e L3. Com g* no espaço forte e J g* = −(B_* + μ*)g* ∈ Y, y_g = (J_A + λ0)g_n ∈ Y.
  O argumento de domínio máximo = mínimo precisa valer no espaço com sinal de S4.

### R6: a contagem

**O que está provado** (computacionalmente ou por álgebra conferida):

- μ* é raiz algebricamente simples de L (realização RT com sinal, pesos (65/64, 5/4), setor A), com
  autoespaço span(g*), g* = (Z*+1)ω* ≠ 0;
- D g* = 0 (identidade de C verificada, mais c(w*) = 0 por RT);
- a álgebra do quociente: com G = span(g*) ⊂ ker(D|E) invariante, dim(E_constraints/G) =
  dim E − rank(D|E) − 1 é álgebra linear correta. Com dim E = 1 e D|E = 0, a contribuição é 0.
  Simplicidade e Dg* = 0 são **ambas** necessárias: sem simplicidade, uma cadeia sobre g* poderia
  contribuir.

**O que é hipótese:**

- (H-rec) e (H-gauge), como enunciadas;
- (H-cone), que é uma escolha e não um teorema;
- a completude da classificação Re s = μ(1 − n). Ela vale só para germes **analíticos** f(u) de
  reparametrizações nulas que preservam o centro. Germes não analíticos (u^α, logaritmos) e gauges
  fora do ansatz ficam de fora por hipótese (parte de H-gauge). As constantes δζ e δφ não têm perfil
  periódico com Re s > 0, segundo GEOMETRIC_FLOQUET_RECONSTRUCTION.

**Hipóteses escondidas ou faltantes:**

- **L1 (aliases de Floquet).** Verifiquei no código (`r6_saida.txt`) que J_m(s + i) = J_{m+2}(s) e que
  os blocos de B só dependem de d e da paridade de m_i. Logo σ(L) é invariante por s → s + i (a
  similaridade é limitada, com fator κ1⁴). "N_A = número de raízes em Re s > 0 com multiplicidade" é,
  literalmente, **infinito**: μ + ij, j ∈ ℤ, são todos raízes com autovetor de gauge. A contagem e a
  subtração "−1 (μ)" só fazem sentido numa faixa fundamental de altura 1 (p.ex. −1/2 < Im s ≤ 1/2).
  Além disso, o setor B (componentes com paridade de m trocada) equivale a L_A(s + i/2). Assim, contar
  L_A numa faixa de altura 1 cobre todos os expoentes físicos módulo i/2. Isso precisa estar enunciado.
- **(H-reg).** A classe física de perturbações é a analítica nesses raios. Alternativamente, é preciso
  que o espectro em Re s > 0 não dependa do espaço. RADIUS_RECOVERY cobre só a troca entre dois raios
  analíticos, e ainda noutro espaço (L3).
- **Mesma realização.** N_A de S3a/S3b precisa ter sido contado com multiplicidade algébrica no mesmo
  espaço, setor e realização "cone móvel" de S4. A escolha de cone precisa valer para todas as raízes,
  s_K inclusive.
- **Constraints nas outras raízes.** Para "N_físico = N_A − 2" é preciso rank(D|E_{s0}) = 0 em toda
  outra raiz da faixa: o espaço de raízes inteiro, não só o autovetor, satisfaz as constraints. Em s_K
  é preciso rank(D|E_{sK}) = dim E_{sK} (= 1 se s_K é simples com D h_K ≠ 0). Isso não está entre as
  três hipóteses e não foi (re)verificado aqui.

**O que faltaria para "N_físico = N_A − 2":** (a) N_A definido numa faixa fundamental e
certificado nela; (b) Σ_{s0 ≠ μ} rank(D|E_{s0}) = 1, com a contribuição só em s_K, e o espaço de raízes
de s_K certificado; (c) nenhum gauge além de μ na faixa, o que decorre de H-gauge; (d) H-rec,
incluindo cadeias e normalizações, e H-reg; (e) a escolha H-cone aplicada coerentemente.

### R7: reprodução

As duas saídas de `nk_L3_C.py --dir rig_copia7` (θ = 0,880172, Z2 = 183,7, Y = 6,591·10⁻⁶,
r = 6,064·10⁻⁵) e de `nk_L3_E.py` (dist = 1,816·10⁻⁵, 0,902443, IDENTIFICADO) são idênticas aos JSON
originais. Os sha256 dos originais ficaram inalterados.

## Lacunas e erros encontrados

Não achei **erro** numérico nem de lógica no certificado de simplicidade e identificação. Todos os
números reproduzem, e as normas verdadeiras que calculei estão abaixo das cotas.

| # | Gravidade | Lacuna | Conserto proposto |
|---|---|---|---|
| L1 | **Média** (enunciado da contagem) | A contagem N_A e a subtração de μ não fixam a faixa fundamental. σ(L) é invariante por s → s + i, e μ + ij são raízes de gauge. | Definir N_A numa faixa de altura 1, com convenção de fronteira. Conferir que a região de S3a/S3b é essa faixa e não conta aliases. Enunciar a equivalência do setor B com L_A(s + i/2). |
| L2 | Baixa–média (rigor escrito) | ‖(Z*+1)Corr‖_L ≤ 10⁻⁶⁰ é afirmado sem argumento. A bola de RT só controla Corr em ℓ¹ nos mesmos raios, e (Z+1) é derivada. Como μ(Z+1) = J(0) − μ∂ξ nas componentes 1 a 3, "O Corr ∈ SSPACE" sozinho não controla ξ∂ξ Corr. | Escrever a cota: O_{μ*}Corr = −[λa + A·Corr + SΩ⁺(Ref) + …] é limitado em ℓ¹ pelos bounds de RT (~10⁻⁸⁰). ∂τ e (D_ξ+1) saem por ‖O⁻¹‖, ‖O⁻¹DO‖ ≤ 𝒦2, e ∂ξ por divisão por (1+ξ) de uma função que se anula em ξ = −1 (a elipse de raio 5/4 dista 0,025 de −1: constante ~10²–10³). Qualquer constante < 10²⁰ basta. |
| L3 | Baixa | O fato de g* estar no domínio forte é citado de RADIUS_RECOVERY.md, escrito para outra realização (pesos com η, ‖B‖ ≤ 5191/500). Não é o espaço com sinal de S4 (‖B_L‖ ≤ 3,964). | Refazer o argumento no espaço de S4, ou usar o argumento de L2, que dá g* ∈ ℓ¹(κ) ⊂ Y diretamente. Mais J g* ∈ Y e D_min = D_max no espaço com sinal. |
| L4 | Cosmética | Em `nk_L3_E.py`, ‖B_*‖ ≤ NORMA_BL + rt omite δμ‖SΓ1‖ ≈ 3,4·10⁻¹⁰. A folga 4,7·10⁻⁸ − 4,6702·10⁻⁸ = 3,0·10⁻¹⁰ não cobre tudo. Efeito em δ: ~10⁻¹⁶. | Somar DELTA_MU·2·norma_divxi a nB (fiz isso em r5 e o resultado não muda). |
| L5 | Baixa (rastreabilidade) | `nk_L3_T.py` lê `--lam0` da linha de comando (padrão 0,401), não de Z.json. T3.json, Q.json e nk3_mu.json não registram λ0. Os 4,7·10⁻⁸ vêm só da variável de ambiente (padrão 40ε_ω). Um descuido geraria um certificado misto sem aviso. | Ler λ0 de Z.json e fazer assert; gravar lam0 e EPS_FUNDO_L em todo JSON. Tornar o erro do fundo um parâmetro explícito de S4. Conferi por vias indiretas (eps_mu, beta, Qfar, ‖V_J‖) que tudo está em μ. |
| L6 | Hipótese (conceitual) | Além de H-rec, H-gauge e H-cone: H-reg; mesma realização e setor para N_A; rank(D\|E) nas outras raízes e em s_K; faixa fundamental (L1). | Enunciar essas hipóteses e certificar rank(D\|E) nas outras raízes da faixa. |

## O que não consegui verificar

- A contagem N_A de S3a/S3b e a certificação de s_K (fora do escopo; não li S3B_ROTA).
- A validação de L e de C contra o exportador de RT (`rt-A-*.dat`). Usei os operadores do código como
  dados e verifiquei só sua consistência interna (identidades exatas, F(RefA) e c(RefA) pequenos).
- As etapas rigorosas de Loewner e Rump (Cholesky de 28 034): não as reexecutei. Comparei com
  estimativas verdadeiras em ponto flutuante, que ficaram abaixo das cotas. O valor NORMA_BL =
  3,964043 (símbolo em Arb, S3b) foi aceito.
- A prova de RT em si (bola 2⁻²⁷⁷) e as álgebras dos documentos de reconstrução geométrica (testes
  racionais do autor), além do argumento de RADIUS_RECOVERY no espaço original.
- A completude global do gauge (H-gauge) e a equivalência física (H-rec): são hipóteses.

## Pastas remotas

As cópias de trabalho ficaram em `~/revisao_s4/` nas duas máquinas. Contêm só cópias de scripts, de
`RefA.dat`, de h0w.npy, de Z.json e T3.json e as saídas. Nenhum processo meu continua rodando.
