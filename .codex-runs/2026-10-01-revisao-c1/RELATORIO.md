# Revisão independente do Rouché de L na faixa A (C1), 01/10/2026

Revisor: agente independente (outro modelo). Objetivo: tentar quebrar a alegação
**N(L; 0 < Re s < 3, |Im s| < 1/4) = 3** (com multiplicidade algébrica).
Estado revisado: commit `8c248b3` (o snapshot `.snapshots/2026-10-01-revisao-c1/MANIFEST.sha256` confere:
nenhum arquivo do manifesto mudou). Todos os artefatos estão nesta pasta (`rN_*.py`, `rN_saida.txt`,
`rN_resultado.json`). Contas remotas: laboratorio2 (`~/revisao_c1/`) e PC 3 (`~/revisao_c1/`).
O PC do laboratório não foi usado. Não li nenhum dos documentos proibidos.

Interrupção: o notebook e o laboratorio2 reiniciaram às ~09h20. Retomei com os mesmos scripts. O R2 do primeiro
ponto (4 valores de s) já tinha terminado no laboratorio2 e foi aproveitado.

## Resumo

Não achei erro que derrube a alegação. Achei **duas lacunas de derivação**. As duas são numericamente
inofensivas, e eu as quantifiquei com código próprio. Também achei uma lacuna de rastreabilidade, que fechei
recalculando, e alguns pontos que não pude verificar (lista no fim).

| R | Veredicto | Artefato | O que foi feito |
|---|---|---|---|
| R1 | **OK com ressalva** | este relatório (§R1) | Prova de que θ < 1 ⇒ H_t invertível. Gohberg–Sigal: H − I é compacto e analítico em Re s > −μ. A referência é uma só função analítica. Identidade da escala rederivada. Ressalva: a matriz de Perron não traz a deflação verdadeira fora das linhas de Z (lacuna L1). |
| R2 | **OK** | `r2_lu_densa.py`, `r2*_saida.txt`, `r2*_resultado.json`, `r2b.sh`, `r2d.sh` | LU densa no laboratorio2 em 27 valores de s, de 6 pontos p de 4 centros (0,25i, 0, −0,25i, 1,0), inclusive s = 0 e o canto 0,25i. Nenhuma violação: max ‖R‖/t(s) = 0,987 e max valor/E(s) = 0,99999. Taylor e os papéis de τ e JQ rederivados. |
| R3 | **OK** (rastreabilidade de 0 e 0,25i fechada por recálculo) | `r3_certF_proprio.py`, `r3_saida.txt`, `r3_resultado.json`, `r3b_*` | Desigualdade derivada. Certificado próprio (Φ por mínimos quadrados via QR, S montado do zero):
  - F do centro 0,25i (sha 4c7a…): ‖Φ‖ ≤ 1,014824 e ρ ≤ 9,2764·10⁻⁵, contra 1,014831 e 9,2774·10⁻⁵ do autor;
  - F do centro 0 (sha 75f8…): 1,014923 e 9,2854·10⁻⁵, contra 1,015157 e 9,2865·10⁻⁵.

  Os certF do autor nesses dois centros **não registram sha256_F**. |
| R4 | **OK** | `r4_sensibilidade.py`, `r4_saida.txt` | As três fórmulas rederivadas. Teste numérico na truncagem 19 modos × n < 40 com o fundo verdadeiro: as identidades valem com erro relativo ≤ 5·10⁻¹⁴ e as cotas valem em todos os s (pior razão verdadeiro/cota = 0,981). |
| R5 | **OK com ressalva** (lacunas L1 e L2) | `r5_perron_proprio.py`, `r5_saida.txt`, `r5_resultado.json`; `r5b_Qfar_s.py`, `r5b_saida.txt`; `r5c_Qfar_disco.py`, `r5c_saida.txt` | Remontagem própria em racionais de 7 discos (3 de 0,25i, 2 de 0, 2 de 1,0). Todas as entradas batem com as do autor (a única diferença é a cota de ‖·‖₂ do bloco T, a minha mais justa). Com os termos de deflação que faltam (L1), θ sobe ~1·10⁻⁴ e fica < 1 (pior 0,99890). A sensibilidade do far Qf_s não está justificada para as linhas de Z/T (L2). Com uma cota válida no disco, Qf diminui e θ cai (pior 0,9983 com L1 e L2 corrigidos). |
| R6 | **OK** | `r6_janelas_pontos.py`, `r6_saida.txt` | Analiticidade (inclusive em Re s = 0, já que Q0 é analítica em Re s > −μ) e partição de Ω conferidas. Princípio do máximo correto para funções analíticas matriciais. Cota recalculada em 2 amostras: 0,1488 e 0,1051, iguais às do autor. |
| R7 | **OK** | `r7_cobertura.py`, `r7_saida.txt` | Verificador próprio, racional, por varredura com teste exato |z − p|² ≤ r² e bisseção. As 4 arestas e os 4 cantos estão cobertos pelos 96 discos. Também conferi ok, θ racional < 1 e r t_p < 1 em cada disco. |
| R8 | **OK** | `r8_schur_residuo.py`, `r8_saida.txt`, `r8_resultado.json` | No laboratorio2 (A, T, U com os sha do manifesto): ‖AU − UT‖_F = 9,0·10⁻¹¹ em ponto flutuante, cota rigorosa ε_B ≤ 6,15·10⁻⁶ (10× abaixo dos 5,97·10⁻⁵ do autor), ‖I − U*U‖_F ≤ 1,3·10⁻⁸, T triangular exata e 2 valores −T_ii em Ω. Argumento de homotopia conferido. max t(s) = 560,3 reproduzido. |
| R9 | **OK com ressalva** | `r9_brauer.py`, `r9_saida.txt` | Identidade s(s − μ)D(s) = q(s) conferida exatamente em 200 pontos racionais. δ = (1, 1 + MU) confere. \|e_ij\| ≤ 1,4·10⁻⁷. min_{Re s ≥ 0} \|q\| ≥ 0,99999997. A indentação dá exatamente +1, e GS conta a multiplicidade algébrica. Ressalvas: ∂τCorr ≤ 10⁻⁶⁰ é uma constante embutida, e a localização do terceiro zero (0,733) não é certificada aqui. |
| R10 | **OK com ressalva** | `r10_consistencia.py`, `r10_consistencia_saida.txt`, `r10_reproducao_saida.txt`, `r10_disco_1_0_saida.txt`, `r10_copia/`, `r10_sha256_antes.txt`, `r10_sha256_depois.txt` | λ0 = centro nos 10 T3. Partição de janelas idêntica nos 10. Pontos do rig = plano. t_p do rig ≥ t_p dos JSONs `tp_*` das máquinas. Todos os discos com certF, e os sha dos certF que o registram conferem. disco 1,0, cobertura e contagem reproduzidos. sha de A igual no laboratorio2 e no João. Ressalva: o A do laboratório (7 dos 10 rigs) não pôde ser conferido. |

## Lacunas (gravidade e conserto)

**L1. Deflação verdadeira fora das linhas de Z' (gravidade baixa).**
- **O problema.** Em `rouche_L_disco.theta_disco`, a deflação verdadeira (x_k exatos, contra os x_k^A de A) só
  entra em Z'←Z', Z'←T, Z'←far (dEZ) e T←Z' (dET). Mas φ_k* Q0(s) **não** se anula nas colunas de T e do far,
  porque Q0 é triangular superior em n (Q_ZT, Q_Zfar ≠ 0). Além disso:
  - P_T x_k ≠ 0 (cauda de RefA nos modos 32..40 mais o erro);
  - P_far x_k = P_far(x_k − x_k^A) ≠ 0.
- **Termos que faltam.**
  - T←T: δ P_T x φ* Q_ZT(s), cota max(w nV)·dET·nQx;
  - T←far: idem;
  - far←Z', far←T e far←far: dEf ‖Q_ZZ(c)‖ e dEf nQx, com dEf = Σ δ‖φ‖‖x − x^A‖.

  Também falta o termo correspondente no bloco diagonal de janela (A6), da ordem de 10⁻⁷.
- **Efeito medido** (`r5_saida.txt`): Δθ ≈ 1,0–1,1·10⁻⁴ nos 7 discos. O pior caso vai a 0,998904 com a minha cota
  justa do bloco T, abaixo dos 0,999158 do autor.
- **Conserto:** somar esses termos em `theta_disco` (e em `rouche_L_janelas`), com dET e dEf já calculados.

**L2. Sensibilidade do far Qf_s = Qf(c)/(1 − d Qf(c)) sem justificativa (gravidade média, numericamente inofensiva).**
- **O problema.** Qf(c) = `cauda_Rinv` cota ‖Q0(c) P_far‖ com **todas** as linhas, porque Q0 é triangular
  superior. A fórmula de Neumann só vale para o bloco far←far:
  P_far Q0(s) P_far = Q_ff(c)(I + (s − c)Q_ff(c))⁻¹.
  As linhas de Z/T são outra coisa:
  P_N Q0(s) P_far = (I + (s − c)Q_NN(c))⁻¹ Q_Nf(c) (I + (s − c)Q_ff(c))⁻¹.
  O fator extra (I + (s − c)Q_NN(c))⁻¹ envolve ‖Q_NN‖, que é ≈ 5 nas linhas de Z.
- **Onde entra.** Qf_s entra em beta_s (far←far), Zf, Tf e desc_s. Uma cota válida porém grosseira, perto de 0,16,
  deixaria θ > 1, então a lacuna não é cosmética.
- **Efeito medido.**
  - `r5b_saida.txt`: o valor verdadeiro (truncagem n < 1200, cota inferior) fica em ≤ 0,92·Qf_s.
  - `r5c_saida.txt`: a forma fechada reescrita com majorantes no disco (|d_k| ≥ max(|d_k(p)| − r, μ(n+1)),
    |1 − v_k/d_k| ≤ min(1, (|d_k(p) − v_k| + r)/|d_k|)) dá uma cota **válida** e **menor** que a do autor
    (por exemplo, 0,131698 contra 0,13396). Com L1 incluído, θ fica entre 0,8508 e 0,9983 nos 7 discos remontados,
    sempre abaixo do θ do autor (tabela em §R5).
- **Conserto:** trocar Qf_s em `theta_disco` pela cota fechada no disco (função `Qfar_disco` de
  `r5c_Qfar_disco.py`) ou por outra derivação que cubra as linhas de Z/T.

**L3. Rastreabilidade do fator (gravidade baixa, fechada para 0,25i).**
- `certF_+0_0000+0_0000i.json` e `certF_+0_0000+0_2500i.json` não têm `sha256_F`: foram gerados por uma versão
  anterior do código.
- Para 0,25i, meu certF próprio no F com sha 4c7a… reproduz ‖Φ‖ e ρ. Esse F é o mesmo do laboratorio2, onde o rig
  rodou (manifesto remoto).
- Para 0, o certF próprio também reproduz ‖Φ‖ e ρ (§R3), no F com sha 75f8…, que é o do laboratorio2.
- **Conserto:** rodar de novo `certF` nesses dois centros com a versão atual (que grava o sha), ou anexar estes
  artefatos.

**L4. A do PC do laboratório não conferido (gravidade baixa a média, rastreabilidade).**
- Sete rigs (−0,25i, 0,55+0,25i, 1,0, 1,5, 2,0, 2,5, 3,0) rodaram no PC do laboratório, cujo A.npy eu não pude
  conferir (máquina proibida). Os JSON do rig não registram o sha de A.
- O schur.json do laboratório tem raízes diferentes (~10⁻¹³) das do laboratorio2/João. Logo T e U foram
  recalculados lá, o que é inofensivo, pois só servem de resolvedor. A pode ter sido regerado ou copiado.
- **Evidência de que A é o mesmo:**
  1. O João regerou A de forma independente (Ryzen) e obteve exatamente o mesmo sha que o laboratorio2 (Intel).
  2. `globais.json` do laboratório (‖B̂‖) é idêntico a 16 dígitos.
  3. Com o A do laboratorio2, o a0 verdadeiro no centro 1,0 (rig do laboratório) é 0,334932, contra a cota
     0,334935. No centro −0,25i é 0,406784, contra a cota 0,407202 (folga 20× maior que nos rigs do
     laboratorio2, mas no sentido certo). ‖R(p)‖ ≤ t_p nos dois.
- **Conserto:** gravar o sha256 de A (e de F) em cada `rig_*.json` e conferir o A do laboratório.

**L5. Arredondamentos não contados (cosmético).**
- Fase 1 do rig: o deslocamento `A[diag] += p` é arredondado antes de A A*, e eP não cobre isso
  (≈ u·max|A_ii + p| ≈ 10⁻¹⁴, efeito relativo ~10⁻¹³ em t_p).
- Em `rouche_L_rig`, o arredondamento de Q_ZZ^k F_b usa normas espectrais onde a cota elemento a elemento pede
  ‖|Q|‖.

Os dois itens ficam ordens de grandeza abaixo das folgas.

**L6. Constantes sem derivação no código (não verificado).**
- `errZ[0] = dtauB + 1e-60`: o 10⁻⁶⁰ cota ‖∂τ Corr‖, mas o código não deriva isso.
- dtauB = ‖∂τ RefB‖ = 1,4563·10⁻⁹ vem de `gauge_refB.py`.
- Não pude conferir a bola de Corr de RT, que está no documento proibido S4.

## R1. Lógica (A1)

**θ < 1 ⇒ H_t invertível.**
- **Norma.** ‖x‖_* = max(‖D_Z x_Z‖/v₀, ‖x_T‖_T/v₁, ‖x_f‖/v₂), com:
  - D_Z = J_Z(c)Q_ZZ(s), invertível;
  - ‖x_T‖_T = ‖(w_J‖x_J‖)_J‖₂, com 0 < w ≤ 1;
  - v racional > 0.
- **Inverso aproximado.** Seja V(s) = diag(Ĥ_ZZ(s)⁻¹, V_J(c_k), I). V é invertível: Ĥ_ZZ(s) = (A + s)Q_ZZ(s), com
  A + s invertível porque t_p é finito, e os V_J são inversos exatos de matrizes calculadas.
- **Bloco a bloco.** I − VH_t = (I − VA) + tV(A − H). Cada bloco tem norma ≤ M_ab para todo t ∈ [0,1], porque
  M_ab majora ‖(I − VA)_ab‖ + ‖V(H − A)_ab‖. Então ‖I − VH_t‖_* ≤ max_a (Mv)_a/v_a = θ.
- **Conclusão.** O autor calcula exatamente essa razão em racionais com um v qualquer, e não precisa do autovetor
  exato. Se θ < 1, VH_t é invertível por Neumann, e H_t = V⁻¹(VH_t) também é.
- **Escalas.**
  - Z←Z: D_Z Ĥ_ZZ⁻¹ ΔH D_Z⁻¹ = J_Z(c)R(s)ΔL Q_ZZ(c), porque Q_ZZ(s)J_Z(s)Q_ZZ(c) = Q_ZZ(c).
  - T←Z: V_J H_JZ D_Z⁻¹ = V_J (B + E)_JZ Q_ZZ(c), independente de s.
  - Z←T: D_Z Ĥ_ZZ⁻¹ H_ZT. Aqui H_ZT = (L̃_ZT − L̃_ZZ Q_ZZ J_ZT) Q_TT, pela inversa triangular por blocos
    (Q_ZT = −Q_ZZ J_ZT Q_TT, que é a identidade J_ZZ Q_ZT + J_ZT Q_TT = 0, pois Q_far,T = 0). Isso dá
    J_Z(c)[R(s)Ỹ − Q_ZZ(s)J_ZT]Q_TT(s).

  Tudo confere com o alegado. A deflação **verdadeira** nas linhas de T e do far nas colunas de T e do far
  não aparece: é a lacuna **L1**.

**Gohberg–Sigal.**
- **Analiticidade de Q0.** Q0(s) é limitada e analítica em Re s > −μ. Pela forma fechada,
  |1 − v_k/d_k| ≤ 1 ⟺ Re σ ≥ −μ.
- **Compacidade.** Q0 é compacta: ‖Q0 P_{n≥N ou |m|≥M}‖ → 0, já que cauda_Rinv ~ 1/(μN) e
  Rinv_uniforme ~ 1/|m|. B é limitado (‖B_L‖ ≤ 3,96). Logo H(s) − I = (B + E)Q0(s) é compacto e analítico perto
  de Ω̄ ⊂ {Re s ≥ 0}.
- **Referência.** A(s) − I tem posto finito: Z mais 26 janelas finitas.
- **Homotopia.** H_t = I + compacto, analítico em s, contínuo em t e invertível em Γ para todo t. Então
  N(H_t, Ω) é constante, pela continuidade do resíduo logarítmico.
- **Zeros de H e de L̃.** H = L̃Q0 tem os mesmos zeros, com as mesmas multiplicidades, que L̃: Q0(s) é um
  isomorfismo analítico X → D(J).

**Decomposição de N(A).**
- A(s) = diag(Ĥ_ZZ(s), Ĥ_J(s), I) é **uma só** função analítica: Ĥ_ZZ = (A + s)Q_ZZ(s) com o mesmo A (dado) e
  Ĥ_J = P_J H P_J.
- Os V_J(c_k) e as escalas D_Z(c_k) são só instrumentos da prova pontual de invertibilidade; trocar de centro ao
  longo do contorno não muda A(s).
- Isso exige a mesma partição de janelas em todos os centros, o que conferi em R10.
- Para operador bloco-diagonal com finitos blocos não triviais, N(A) = N(Ĥ_ZZ) + Σ_J N(Ĥ_J).
- N(Ĥ_ZZ) = #zeros de det(A + s) em Ω (com multiplicidade), porque det Q_ZZ(s) ≠ 0 em Re s > −μ.

## R2. Nível Z por ponto (A2)

**Derivação.** Com e = s − p, R(s) = R(p)(I + eR(p))⁻¹ = R(p) − eR(p)² + e²R(p)³ − e³R(p)³R(s) (resolventes
comutam), e o mesmo vale para Q_ZZ. Portanto:
‖J(c)[R(s)F_t − Q(s)F_b]‖ ≤ a0 + |e|b1 + |e|²b2 + |e|³(‖J(c)R(s)‖‖Y₃‖ + ‖J(c)Q(s)‖‖Q(p)³F_b‖).
- **τ.** J_Z(c)R(s) = I + (c − s − B̂)R(s), logo τ = 1 + (|c − s| + ‖B̂‖)t(s), com t(s) ≤ t_p/(1 − |e|t_p).
- **JQ.** J_Z(c)Q_ZZ(s) = I + (c − s)Q_ZZ(s), logo JQ = 1 + d nQZ/(1 − r nQZ).
- **Erros de Schur.** Y_k − X_k = R(p)(Y_{k−1} − X_{k−1} − Res_k).

Tudo confere com `rouche_L_rig.py`.

**Numérico** (LU densa e Lanczos no laboratorio2, A com sha 70791ce0):

| centro | p (θ do disco) | s testados | max ‖R(s)‖ / t(s) | max valor / E(s) |
|---|---|---|---|---|
| 0,25i | 0,1258i (0,99916, o pior) | p, p ± r/2, p ± 0,999r | 0,888 | 0,99996 |
| 0,25i | 0,25i (canto) | p, p − r/2, p − 0,999r, p + r/2, p + 0,999r (aresta superior) | 0,848 | 0,99995 |
| 0 | 0,0047i (contém s = 0) | p, p ± r/2, p ± 0,999r, **s = 0** | 0,909 | 0,99996 |
| 0 | 0,0882i (0,99915) | p, p ± r/2, p ± 0,999r | 0,891 | 0,99995 |
| 1,0 (rig do laboratório) | 0,8201+0,25i | p, p + r/2, p + 0,999r | 0,987 | 0,99999 |
| −0,25i (rig do laboratório) | −0,25i | p, p + r/2, p + 0,999r | 0,848 | 0,99897 |

Não houve nenhuma violação em 27 pontos. Em s = 0, ‖(A + 0)⁻¹‖ = 36,08 ≤ 48,76: a deflação funciona.
Observação: ‖J_Z(c)R(s)‖ verdadeiro ≈ 40–44, contra τ ≈ 1000–2100. τ é muito folgado e custa ≈ 0,14 na entrada
Z←T (termo (τ + JQ)ρ). Isso é melhoria possível, não erro.

## R3. Fator de posto baixo (A3)

**Desigualdade.** XS = XFΦ + X(S − FΦ), logo ‖XS‖ ≤ ‖XF‖‖Φ‖ + ‖X‖‖S − FΦ‖ para qualquer Φ. Com
X = [J_Z(c)R(s), −J_Z(c)Q_ZZ(s)]:
- XS_c = J_Z(c)[R(s)ỸQ_TT(c) − Q_ZZ(s)J_ZT Q_TT(c)];
- ‖X‖ ≤ τ + JQ.

O fator (I + (s − c)Q_TT(c))⁻¹ dá o 1/(1 − d qTT). Confere.

**Certificado próprio** (`r3_certF_proprio.py`, no PC 3). Montei S do zero: ctx.B como dado, Q0_m(c) por
substituição triangular com cota pelo resíduo e J pesado. Φ_m por mínimos quadrados (QR de F, sem Tikhonov).
‖Φ‖² por um G simetrizado, com autopares e cota de Ostrowski. ρ por Frobenius, com os erros de S e do produto.

| centro | F (sha) | ‖Φ‖ meu / autor (λ = 0) | ρ meu / autor | erro de S meu / autor |
|---|---|---|---|---|
| 0,25i | 4c7a7517… (= laboratorio2, onde o rig rodou) | 1,014824 / 1,014831 | 9,2764·10⁻⁵ / 9,2774·10⁻⁵ | 1,235·10⁻⁶ / 1,245·10⁻⁶ |
| 0 | 75f8bc16… (= laboratorio2) | 1,014923 / 1,015157 | 9,2854·10⁻⁵ / 9,2865·10⁻⁵ | 1,316·10⁻⁶ / 1,327·10⁻⁶ |

Também conferi:
- sha256 dos F locais = `sha256_F` dos 8 certF que o registram = manifesto do laboratorio2 (para 0, ±0,25i) e
  certF_extra (0,55+0,25i, 1,0…3,0);
- o `resto_F`/`qTT` gravados em cada rig = os do `F_*.json` correspondente, para os 10 centros.

## R4. Sensibilidade da cauda (A4)

**Derivação.** Para janela central, Q0(s)P_J = P_Z Q_ZJ(s) + P_J Q_JJ(s). Além disso:
- Q0(c)P_J Q_JJ(c) = Q0(c)²P_J − P_Z Q_ZZ(c)Q_ZJ(c);
- (I + (s − c)Q_JJ(c))⁻¹ = I − (s − c)Q_JJ(s).

Daí V_J P_J B(Q0(s) − Q0(c))P_J = −(s − c)[V_J P_J B P_Z Q_ZZ(c)Q_ZJ(s) + (V_J P_J B Q0(c)²P_J −
V_J P_J B P_Z Q_ZZ(c)Q_ZJ(c))(I − (s − c)Q_JJ(s))]. As normas dão a fórmula alegada, com fJ ≥ 1 + d‖Q_JJ(s)‖.
- **Janela não central:** Q0(s)P_J = P_J Q_JJ(s), o que dá d K_J fJ.
- **Fora da diagonal:** V_i P_i B Q0(s)P_j = N_ij(c) − (s − c)[TZ_i-termo Q_Zj(s) + V_i P_i B Q0(c)P_j Q_jj(s)],
  e N_ij(c) usa Q0 P_j com todas as linhas, o que conferi em `nk_L2_T.produto_cert`.
- **q_ZJ(s):** Q_ZJ(s) = −(I + (s − p)Q_ZZ(p))⁻¹[Q_ZZ(p)J_ZJ Q_JJ(c)](I − (s − c)Q_JJ(s)), donde
  q_ZJ(s) ≤ qz(p)/(1 − r nQZ) · fJ. Confere com o código.
- **TZ_J e K_J** cobrem o B verdadeiro, pois incluem pert·nQ.

**Numérico** (`r4_sensibilidade.py`). Truncagem pequena do operador verdadeiro: 19 modos, n < 40, Z = |m| < 4 com
n < 12, 6 janelas de 4 modos, c = 0,25i e d ∈ {0,02; 0,05; 0,1}.
- Erro relativo das identidades ≤ 5·10⁻¹⁴.
- Razão verdadeiro/cota ≤ 0,981 (janelas não centrais: a cota é quase exata, como esperado).
- Fora da diagonal: o incremento verdadeiro fica sempre abaixo do da cota.

## R5. Montagem por disco (A5)

**Remontagem.** `r5_perron_proprio.py` monta a matriz 3×3 em racionais a partir dos JSONs.
- Fiz a minha própria derivação das entradas (§R1) e a minha própria cota de ‖Ns‖₂: Collatz–Wielandt exato em
  Ns^T Ns.
- Só os q_ZJ certificados vêm das funções do autor. Conferi esses valores em ponto flutuante com código próprio:
  0,196631 contra 0,196700 e 0,196458 contra 0,196461.
- Comparei termo a termo com a docstring:
  - τ, E, EF (melhor par de certF);
  - (dB + dEZ)·nQx + emu_s;
  - 1/(1 − d qTT);
  - os pesos por wmin_SZ e wmin0;
  - eg = rt + emu_s em toda entrada;
  - dpert;
  - TZ + nV(dpert + dET)nQZ;
  - Zf com as descidas e nE·descida(NZ, Nc);
  - beta_s, desc_s, resto·nQb.

  Também conferi a perturbação do fundo (dB = arred40 + rt, resto = 0 na banda 40 de Z), μ (emu_c em Z←Z, emu_s
  nas demais) e dEZ/dET/nE.

| disco | θ autor | θ meu (mesmos termos) | θ + termos L1 | θ + L1 + Qf válido (L2) |
|---|---|---|---|---|
| 0,25i, p = 0,1258i | 0,999158 | 0,998797 | 0,998904 | 0,995666 |
| 0,25i, p = 0,25i | 0,998521 | 0,998318 | 0,998422 | 0,998313 |
| 0,25i, p = 0,1341+0,25i | 0,999073 | 0,998708 | 0,998815 | 0,992124 |
| 0, p = 0,0047i | 0,992247 | 0,992012 | 0,992115 | 0,991779 |
| 0, p = 0,0882i | 0,999154 | 0,998879* | 0,998984* | 0,996471 |
| 1,0, p = 0,8201+0,25i | 0,862385 | 0,862034 | 0,862131 | 0,858969 |
| 1,0, p = 0,8906−0,25i | 0,855099 | 0,854753 | 0,854849 | 0,850781 |

(*) Na primeira rodada, com a cota de ‖·‖₂ do bloco T em ponto flutuante, os valores foram 0,998879/0,998984. A
segunda rodada (Collatz–Wielandt) está em `r5_saida.txt`.

A última coluna (`r5c_saida.txt`) troca Qf_s pela cota válida no disco: 0,1197 a 0,1318, contra 0,134 do autor nos
discos de 0 e ±0,25i. Com L1 e L2 corrigidos, θ < 1 com folga ≥ 1,7·10⁻³ nos 7 discos.

Todas as entradas, menos T←T e far←T, coincidem com as do autor a ≤ 10⁻⁷ relativos. T←T e far←T diferem por
≤ 10⁻³ relativos: a cota rigorosa do autor (`cota_norma`) é mais folgada que a minha (Collatz–Wielandt, exata em
racionais). Termos procurados e ausentes: os de **L1** (deflação verdadeira fora de Z') e a justificativa de
**L2** (Qf_s). Não achei arredondamento faltando que importe (L5 é cosmético).

## R6. Janelas sem zeros em Ω (A6)

- **Partição.** Os 10 retângulos de `REGIOES` cobrem [0,3] × [−1/4,1/4]. As colunas são [0;0,3], [0,3;0,8],
  [0,8;1,2], [1,2;1,8], [1,8;2,2], [2,2;2,8], [2,8;3,0]; na faixa à esquerda, [0,15;0,25], [−0,15;0,15] e
  [−0,25;−0,15]; em [0,3;0,8], [0;0,25] e [−0,25;0].
- **Analiticidade.** F_k(s) = I − V_J(c_k)Ĥ_J(s) é analítica numa vizinhança de cada retângulo fechado, pois Q0 é
  analítica em Re s > −μ e o retângulo toca Re s = 0.
- **Princípio do máximo.** Vale para funções analíticas a valores matriciais: ‖F(s)‖ = sup_{|u|=|w|=1}|w*F(s)u| é
  supremo de módulos de funções analíticas, logo subarmônica. O argumento está certo.
- **Amostras do bordo.** A cota por amostra com Neumann em torno de q (Q_ZZ, Q_JJ, Q_ZJ(s) =
  (I + (s − q)Q_ZZ(q))⁻¹Q_ZJ(q)(I + (s − q)Q_JJ(q))⁻¹) está correta.
- **Recálculo próprio** (`r6_saida.txt`):
  - centro 0, s = 0,3+0,1i: 0,1488 (autor 0,1488);
  - centro 0,25i, s = 0,3+0,2i: 0,1051 (autor 0,1051);
  - extra, s = 0: 0,038.
- **Ressalva:** o termo de deflação de L1 também falta aqui. É da ordem de 10⁻⁷ e não muda nada.

## R7. Cobertura (A7)

`r7_cobertura.py` é independente do script do autor:
- varredura por aresta com teste racional exato |z − p|² ≤ r² e bisseção exata;
- cantos testados um a um;
- validação de ok, θ (Fraction) < 1, r > 0, r·t_p < 1 e p no contorno.

Resultado: 96 discos válidos, as 4 arestas cobertas (25, 17, 2 e 17 passos) e os 4 cantos dentro de 2, 3, 1 e 1
discos. CONTORNO COBERTO.

## R8. Contagem finita (A8)

- **Resíduo de Schur** (`r8_schur_residuo.py`, no laboratorio2). Calculei os sha de A, T e U no próprio processo:
  conferem com o manifesto remoto.
  - ‖AU − UT‖_F (ponto flutuante) = 9,05·10⁻¹¹; a cota global de arredondamento domina e dá
    ε_R ≤ 6,15·10⁻⁶;
  - ‖I − U*U‖_F ≤ 1,28·10⁻⁸;
  - ε_B ≤ 6,15·10⁻⁶;
  - T é exatamente triangular superior;
  - −T_ii em Ω: 0,401025 e 0,733180. O mais próximo do bordo de fora é −0,1683, a 0,168 do bordo.
- **Argumento.** B := UTU⁻¹ tem espectro exatamente {T_ii}, porque U é invertível (‖I − U*U‖ < 1). No contorno,
  ‖(A + s)⁻¹‖ ≤ t(s) e t(s)ε_B < 1. Então A + s + τ(B − A) é invertível no contorno para τ ∈ [0,1], e det(A + s)
  e det(B + s) têm o mesmo número de zeros em Ω. Correto.
- **max t(s).** 560,289 (reproduzido), com t·ε_B = 0,0334 (autor) ou 0,0034 (com o meu ε_B).
- **Atenção.** `rouche_L_contagem.rouche` pula discos com ok = False ao tomar o max t. Isso é inofensivo porque a
  cobertura só usa discos ok e todos os 96 são ok. Mas o script deveria falhar nesse caso, e não pular.
- **Outro caminho (não rigoroso, não feito).** Não rodei eigvals(A) nem o argumento de det por LU no contorno. A
  contagem pela diagonal de T, com o resíduo próprio, já é independente da do autor.

## R9. Brauer (A9)

- **D(s).** L(s)⁻¹x₀ = x₀/s e L(s)⁻¹x₁ = x₁/(s − μ*). Com P_ij = φ_i*x_j = δ_ij + e_ij, vem
  D(s) = det(I + [δ_j P_ij/(s − s_j)]) = q(s)/(s(s − μ*)). Conferido exatamente em 200 pontos racionais.
- **δ.** δ = (1; 1,1683070790022612) = (1, 1 + MU), com MU = 722873400/2³².
- **e_ij.** |e_ij| ≤ ‖φ_i‖‖x_j − x_j^A‖, com ‖φ‖ = (5,50; 0,789) e ‖x − x^A‖ ≤ (1,456·10⁻⁹ + 10⁻⁶⁰; 2,51·10⁻⁸):
  |e| ≤ (8,0·10⁻⁹, 1,4·10⁻⁷; 1,2·10⁻⁹, 2,0·10⁻⁸).
- **Zeros de q.** Para Re s ≥ 0, |q(s)| ≥ (1 − |e₀₀|)(1 − δ₁|e₁₁| − Δμ) − δ₀δ₁|e₀₁e₁₀| ≥ 0,99999997 > 0. q não tem
  zeros em Ω̄ (os zeros estão perto de −1).
- **Polos de D.** μ* ∈ [0,168307078963; 0,168307079041] ⊂ Ω, e 0 está no bordo.
- **Indentação.** Γ_ε exclui 0 por um semicírculo dentro de Ω. L̃ é invertível numa vizinhança de Γ
  (compacto), e L tem zeros isolados. Então M(L̃, Γ_ε) = M(L, Γ_ε) + [zeros − polos de D em Ω_ε] =
  M(L, Γ_ε) − 1, logo N(L, Ω) = N(L̃, Ω) + 1 = 2 + 1 = 3. É exatamente +1: o polo em 0 fica fora, e o de μ* é
  simples porque q(μ*) ≠ 0.
- **Multiplicidade.** GS entrega a soma das multiplicidades parciais, ou seja, a multiplicidade algébrica dos
  autovalores do pencil linear L(0) + s, isto é, de −L(0).
- **Ressalvas.**
  1. A alegação "μ ≈ 0,1683, 0,4010 e 0,7332" vai além do Rouché, que dá só a contagem 3. A identificação de
     μ* (S4) e de 0,401 (S3b) vem de outros certificados. O terceiro zero ≈ 0,733 só tem localização estimada:
     é o autovalor da matriz de ponto flutuante A, que é rigoroso para A mas não para L. Para L, o que se tem é
     "existe exatamente mais um zero em Ω, contado com multiplicidade".
  2. A constante 10⁻⁶⁰ (L6).

## R10. Consistência e reprodução

- **Consistência** (`r10_consistencia_saida.txt`): nenhuma inconsistência.
  - Os 10 T3 têm λ0 = centro, eps_fundo_L = 4,7·10⁻⁸ e rouche = true. Têm partição de janelas (nome, modos, SZ,
    folga_far) e dB_L idênticos.
  - Os pontos do rig são os do `plano_rig.txt`.
  - nBhat = 3,850313347957818 em todos.
  - t_p do rig ≥ t_p de `tp_*_rig.json` dos centros 0, ±0,25i e 0,55−0,25i (copiados do laboratorio2/João). A
    diferença é 0 ou +10⁻⁶ (o arredondamento para cima).
  - Todos os discos usam certF do mesmo centro, com pontos e t_p iguais aos do rig.
  - `defl_refA` e `globais.json` são idênticos nas 3 máquinas.
- **sha de A.** 70791ce0… no laboratorio2 (calculado em R2/R8) e no João (`sha256sum`). Laboratório: não
  verificado (L4).
- **Reprodução** (`r10_reproducao_saida.txt`):
  - `rouche_L_disco.py` no centro 1,0, numa cópia: os 10 raios idênticos; θ idêntico em 9 pontos e diferente
    em 7·10⁻¹⁶ no primeiro, por arredondamento de BLAS no `cota_norma`;
  - `rouche_L_cobertura.py`: saída idêntica à gravada;
  - `rouche_L_contagem.py rouche`: idêntico (560,289; 3,344·10⁻²).
- **sha256 dos originais.** Antes em `r10_sha256_antes.txt` e depois em `r10_sha256_depois.txt` (294 arquivos:
  `build/rouche_L_rig`, `build/rouche_L_est`, `scripts/*.py`, `E.json`). Os dois são **idênticos**, e o
  `MANIFEST.sha256` do snapshot confere sem falhas no fim da revisão.

## O que não consegui verificar

1. O A.npy (e o F usado) no **PC do laboratório**, onde rodaram 7 dos 10 rigs (L4). Só há evidência indireta.
2. A maquinaria já revisada em s real foi tomada como dada: `verified_inverse.cota_loewner`, `cota_norma`,
   `nk_L3_T` (Loewner por janela), `fundo_L.dB_L` e `NORMA_BL` (Arb), `closed_form` (descida, cauda_Rinv) e
   `eps_mu_L`. Só rederivei os usos novos.
3. A bola de Corr de RT e o 10⁻⁶⁰ de ∂τCorr (L6); o `E.json` de S4 foi tomado como dado.
4. Não reexecutei `nk_L3_T --rouche` (T3) nem `rouche_L_janelas.py` completo, que custam horas. Também não rodei
   o rig rigoroso (`rouche_L_rig ponto`) de nenhum ponto. Em vez disso, conferi os valores verdadeiros por LU
   densa (R2).
5. A contagem de autovalores de A por outro caminho (eigvals, ou argumento de det por LU no contorno). Usei o
   resíduo próprio do Schur.
6. Não remontei θ nos outros 89 discos (só em 7). Os de 0,25i e 0 têm folga de ~8·10⁻⁴, contra ~1·10⁻⁴ de L1,
   e L2 só diminui θ.
