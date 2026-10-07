# T2: o enunciado e as hipóteses (02/10/2026)

**T2.** O fundo crítico de Reiterer–Trubowitz (RT) tem exatamente um modo instável físico.

A parte computacional está fechada e revisada:
- faixa A: `C1_REAVALIACAO.md` §9–10;
- faixa B: §11 e §11.1.

Este documento separa o que é teorema (computacional ou em papel) do que é hipótese, e diz o que falta
para cada item.

## 1. O enunciado condicional

**Teorema (T2, condicional).** Seja L(s) o operador linearizado de RT em torno do fundo crítico, com
μ* > 0 o expoente do eco. Vale a modelagem H-reg (em 03/10, H-gauge virou o Lema G, H-rec o Lema R e H-cone o Lema C). Então, por classe de
Floquet do problema de período completo, o número de modos instáveis físicos (Re s > 0), contados com
multiplicidade algébrica e módulo gauge, é **1**. O modo é real, com s ≈ 0,7332 (realidade: `T2_ENUNCIADO.md`).

**Cadeia da prova:**
1. **Redução ao setor A** (`SECTOR_B.md` Lema 5.1): basta contar L_A na faixa
   Q̃ = {−1/4 < Im s <= 3/4} ∩ {Re s > 0}.
2. **Nada em Re s >= 2,93** (F3: campo de valores no toro, para todo Im s).
3. **Faixa A:** N_A = 3, com as raízes μ*, s_K ≈ 0,4010 e ≈ 0,7332 (Rouché de L, discos e contagem de
   Schur).
4. **Faixa B:** N_B = 1, com a raiz ≈ 0,0414 + 0,5i (mesmo maquinário).
5. **μ*:** raiz simples com autoespaço span(g*), g* = (Z* + 1)ω*, e Dg* = 0 (S4). É gauge pelo cone
   móvel e contribui 0.
6. **s_K:** raiz simples cujo autovetor viola as constraints (S3b: testemunho 0,1685), contribui 0.
7. **Raiz da faixa B:** simples, viola as constraints (testemunho 0,3547), contribui 0.
8. **0,7332:** K é injetivo ali (S3a), logo, por KC = ML, todo o espaço de raízes satisfaz as
   constraints. Contribui dim E − rank(D|E) = 1.
9. **Soma:** 1.

## 2. O que é teorema

| item | conteúdo | onde | natureza |
|---|---|---|---|
| E1 | existência do fundo DSS na bola de RT | artigo de RT | publicado (prova assistida) |
| F1 | H(s) = I + (B + s)Q família analítica de Fredholm de índice 0 | `RADIUS_RECOVERY.md`, ledger | papel, auditado em 16/09 |
| F3 | R_L <= 2,926 (campo de valores no toro) | `ALTERNATIVE_ROUTES.md` §5.26 | certificado em Arb |
| C1 | N_A = 3 | `C1_REAVALIACAO.md` §9–10 | certificado, revisado |
| C4 | N_B = 1, testemunho de violação | §11 | certificado, revisado |
| S3a | K injetivo em 0 <= Re s <= 1,765, \|Im s\| <= 1/4, fora do disco de S3b | `S3_CONDICIONAMENTO_SHARP.md` §10 | certificado, revisado |
| S3b | s_K raiz simples de K; testemunho em L | `S3B_ROTA.md` | certificado, revisado |
| S4 | μ* simples, autoespaço span(g*), Dg* = 0 | `S4_GAUGE.md` | certificado, revisado |
| KC = ML | identidade de propagação das constraints | `SHARP_LINEARIZED_IDENTITY.md` | papel, auditado |
| quociente | dim E − rank(D\|E) − [gauge] | `CONSTRAINT_ROOT_QUOTIENT.md` | papel (álgebra linear) |
| gauge local | classificação do gauge residual analítico; expoentes μ(1 − n) | `GAUGE_RT_ACTION.md`, `GAUGE_DOMAIN_CLASSIFICATION.md` | papel, local |
| redução B → A | M⁻¹L_B(s)M = L_A(s + i/2) e Lema 5.1 | `SECTOR_B.md` §2–5 | papel |

## 3. As hipóteses, uma a uma

### 3.1 Lacunas matemáticas: fechadas (02/10/2026)

**Lema L-real (as contagens no toro valem na realização de RT).**
- **Enunciado.** Seja |s| <= 3,1, o que inclui toda a região Ω_A ∪ Ω_B = [0, 3] × [−1/4, 3/4]. Os
  núcleos de L(s) e as cadeias de Jordan coincidem em duas realizações: no ℓ² do toro com sinal,
  ‖x‖² = Σ κ1^{2m} ω2(n) |x_mn|² com (κ1, κ2) = (65/64, 5/4), e na realização ℓ¹ de RT, Y+, com os
  mesmos raios. Logo os autovalores e as multiplicidades algébricas na região são os mesmos nas duas.
- **Prova** (recuperação de raio, como em `RADIUS_RECOVERY.md`):
  - Seja Q = J_μ*⁻¹, o mesmo operador nos dois espaços e independente de s. Escreva
    H(s) = I + (B + s)Q, G = {\|m\| < 1024, n < 4096} e T = I − G.
  - **Toro** (`scripts/l_real_toro.py`, `build/t2/l_real_toro.json`):
    ‖T(H − I)T‖ <= (‖B_L‖ + \|s\|)·‖J⁻¹T‖ <= (3,9641 + 3,1)·0,0525 = 0,371 < 1.
    - ‖B_L‖ <= 3,9641 vem do símbolo de F3, em Arb.
    - ‖J⁻¹T‖ <= 0,0525 sai das formas fechadas: a cauda n >= 4096 nos modos \|m\| < 400 (máx
      0,0153) e a cota do modo inteiro, (1 + 2c_w/(κ2 − 1))/(\|m\|/2), nos \|m\| >= 400.
    - μ entra com μ_RefA ± 10⁻¹⁰.
  - **RT** (`radius_transfer.py`, `build/t2/radius_transfer_s3_1.json`):
    ‖T(H − I)T‖_+ <= 0,2665 (forte), com β+ <= 5191/500 da bola publicada de RT.
  - **Inclusão.** Y+ ⊂ toro, com ‖x‖_toro <= ‖x‖_+/243, porque os pesos η de Y+ são >= 243 (revisão de T2, §2.3).
  - **Núcleos.** Seja H(s)x = 0, com x no toro. Gx é finito, logo está em Y+, e THGx está em Y+, porque
    J⁻¹ é triangular em n e preserva a caixa, e B é limitado em Y+.
    - Em Y+, THT é invertível (Neumann), e y = −(THT)₊⁻¹THGx está em Y+, logo também no toro.
    - No toro, THT também é invertível, e a solução de THT·z = −THGx é única. Logo Tx = y, x está em Y+
      e h = Qx está no domínio forte.
    - A recíproca é a inclusão.
  - **Cadeias.** Para L(s0)h_j = −h_{j−1}, por indução o lado direito é forte, e o mesmo argumento
    vale.
- **Consequências** (com a equivalência por Q0(s), que é invertível e analítica em Re s > −μ, e o
  lema de Brauer, N(L) = N(L̃) + 1 na faixa A):
  - os Rouchés no toro (N_A = 3 e N_B = 1) contam as raízes de L na realização de RT, com as mesmas
    multiplicidades;
  - os autovetores certificados no toro (S3b, S4 e o testemunho da faixa B) são os de RT;
  - a injetividade de K no toro (S3a) implica a injetividade em RT, pela inclusão. Mais precisamente
    (03/10): C(s)h ∈ D_K(129/128, 9/8) ⊂ toro♯(129/128, 9/8), onde estão os certificados de K
    (`SHARP_DOMAIN.md`, fechamento).
- **Condicional à bola publicada de RT** (E1), como todo o projeto.

**Lema L-afim (L(s) = L(0) + s·I).**
- O sistema de RT é O_μ ω + (não linear) = 0, com O_μ = ∂τ + μ D_O (`corr_gauge_bound.py`, a partir
  do TeX de RT).
- A perturbação e^{sτ}h, linearizada em torno de ω*, dá (∂τ + s + μ D_O + B_*)h: a única ∂τ está
  em O_μ, com coeficiente 1, e s entra com coeficiente identidade.
- Na base com sinal, ∂τ → im/2 e s entra só na diagonal de J: σ = s + im/2 em `J_livre`.
- Portanto L(s) = L(0) + s·I exatamente, e o mesmo vale para K(s) e C(s) = C0 + sC1
  (`signed_constraints.py`; C1 = derivada do resíduo ♯ em s).

**Lema L-iso e L-constr (a redução B → A na realização com sinal).**
- **Setores na base com sinal** (`signed_operator_L.Radial`):
  - setor A: componentes 1 a 3 nos m pares (a 1ª só com n ímpar) e componente 4 nos m ímpares;
  - setor B: as paridades trocadas.
- **M = multiplicação por e^{iτ/2}** é o deslocamento m → m + 1. Ele leva A em B bijetivamente:
  - no toro com sinal, ‖Mx‖ = κ1·‖x‖, uma isometria a menos de constante;
  - preserva o domínio, porque ∂τ M = M(∂τ + i/2);
  - comuta com a multiplicação por qualquer função de τ, e portanto com B e com Γ2(ω*, ·);
  - comuta com as operações em ξ.
- Logo M⁻¹L_B(s)M = L_A(s + i/2) e M⁻¹C_B(s)M = C_A(s + i/2). Isso é o Lema 2.3 de `SECTOR_B.md`, que
  só usa ∂τ ↦ ∂τ + i/2.
- **Na base:**
  - J_m(s) depende só de σ = s + im/2, então J_{m+1}(s) = J_m(s + i/2);
  - os blocos de B dependem só de d = m_out − m_in e da componente, que o deslocamento preserva
    (R1 da revisão da faixa B conferiu a i-periodicidade, m → m + 2, no código, com erro 10⁻¹⁶).
- **Consequências:**
  - por Gohberg–Sigal, a similaridade preserva as multiplicidades algébricas, e o Lema 5.1 de
    `SECTOR_B.md` vale nesta realização;
  - um autovetor do setor B satisfaz as constraints se e só se o transportado no setor A satisfaz;
  - logo o testemunho da faixa B mostra que o modo real do setor B viola as constraints.
- H-dof, a seleção de graus de liberdade no código de RT, não se aplica, porque os certificados usam a
  base própria.

### 3.2 Hipóteses físicas (ficam no enunciado)

- **(H-rec) Reconstrução: provada em 03/10.** Ida: Lema R (`HREC_IDA.md`, revisado). Volta:
  `GEOMETRIC_FLOQUET_RECONSTRUCTION.md`. Texto original do item:
  - Ida, a parte física: toda perturbação linear física (Einstein–campo escalar, esfericamente
    simétrica, na classe de H-reg) pode ser posta no ansatz de RT por um gauge admissível.
  - Volta: todo elemento de ker L(s) que satisfaz as constraints vem de uma perturbação física.
    `GEOMETRIC_FLOQUET_RECONSTRUCTION.md` faz a volta dentro do ansatz para Re s > 0.
  - Falta a ida, que é uma afirmação sobre o problema geométrico, não sobre L.
- ~~(H-gauge) Completude do gauge residual~~: **provada em 03/10** (Lema G, `GAUGE_GLOBAL.md`).
  - Para gauges C¹ que preservam o ansatz, os modos de gauge de Floquet com Re s > 0 são ℂ·g*.
  - A prova é global, sem analiticidade, com o fato certificado F1.
- **(H-cone) Cone móvel: não é mais hipótese desde 03/10** (Lema C, `GAUGE_GLOBAL.md` §6). Texto original: a instabilidade física é medida módulo a mudança do instante T* da
  singularidade. É a convenção da literatura de fenômenos críticos. Com cone fixo, μ* passaria a
  contar.
- **(H-reg) Regularidade:** a classe física de perturbações é a analítica nos raios de RT. Com
  L-real, as contagens não dependem de qual dos espaços analíticos (Y+, Y−, toro) se usa.

## 4. Ordem de trabalho

1. ~~L-real~~, ~~L-afim~~, ~~L-iso e L-constr~~: fechados (§3.1).
2. Enunciado final de T2 (`T2_ENUNCIADO.md`), com as hipóteses físicas H-rec, H-gauge, H-cone e H-reg
   explícitas e a cadeia de referências.
3. Revisão independente deste documento e do enunciado.
