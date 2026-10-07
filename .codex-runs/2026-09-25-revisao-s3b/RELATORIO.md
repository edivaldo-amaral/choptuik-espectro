# Revisão independente de S3b — relatório (25/09/2026)

Revisor: agente independente (Claude). Li só o permitido pela TAREFA: `scripts/*.py`,
`build/s3b/**`, o artigo e o código de RT e o PDF de Rump. Não abri `docs/S3B_ROTA.md`,
`RETOMADA.md`, `PROOF_OBLIGATIONS.md`, `AUDITORIA_*` nem relatórios anteriores. Também não
alterei nenhum arquivo existente: os sha256 de `build/s3b/rig/{Z,T3,C3,D}.json` e `h0w.npy`
continuam os mesmos antes e depois das rodadas (`r9_saida.txt`). Tudo o que escrevi está em
`.codex-runs/2026-09-25-revisao-s3b/`. No PC do laboratório trabalhei só em `~/revisao_s3b/`.

## Veredicto geral

**Não achei erro que derrube o certificado.** Reproduzi θ, Z2, Y e r bit a bit, a partir de
uma cópia dos dados. Também remontei os três níveis com código meu e cheguei aos mesmos números.
Medi a norma verdadeira de E = I − A DF(x0) e de ‖V_J‖ em 19 comparações, todas abaixo das cotas de
`T3.json`. Recontei o autovalor de K̂_Z com um montador independente e recalculei o testemunho.

Há cinco lacunas **formais** (L1–L5). São termos de E ou do resíduo que o código esquece ou
justifica mal. Quantifiquei todas: com as quatro primeiras corrigidas, θ não muda na sexta casa
e o NK continua fechando. A quinta é uma hipótese não provada no código, que conferi
numericamente.

O ponto que exige atenção é outro: **a margem é mínima**.

- Com Z2 e Y fixos, θ = 0,839826 ainda pode subir no máximo 7,5·10⁻⁴.
- O produto Y·Z2 pode crescer só 0,94%.
- Y é 99,93% o erro dos dados de RT (t_V · 40 ε_ω). O fator 40 pode subir no máximo até
  **40,38** antes de o NK falhar.

Por isso conferi a constante 40 a partir da definição da norma de RT (artigo, §"Weighted ℓ¹
norms"). A cota de Young dá ||δB|| ≤ max_j ||C_j||₂ · ||δω||_RT = **5,01 ε_ω**, oito vezes
abaixo de 40 (`r6_saida.txt`, parte c). Logo a escolha 40 é segura.

Uma entrada que não reverifiquei pesa no mesmo nível de sensibilidade: `NORMA_BL = 3,964043`,
a cota do símbolo em Arb. Um erro de ~0,2% nela bastaria para passar da margem de θ (ver
"O que não verifiquei").

## Tabela R1–R9

| R | Alegação | Veredicto | Artefato | O que fiz |
|---|---|---|---|---|
| R1 | A1 (Rouché de K) | **OK com ressalva** | `r1_rouche.py`, `r1_saida.txt` | Derivei N' (abaixo). Conferi em Arb, por método próprio (sobreposição dos arcos cobertos pelos discos), que os 16 discos cobrem a circunferência; a menor sobreposição é 3,9·10⁻⁷ rad. Conferi numericamente que Q0 P_Z cai em Z. Montei H(s) e A(s) em ponto flutuante, com código meu, numa caixa 40×200 no pior tile (j = 8, s = 0,351), e medi os blocos de A⁻¹(H − A) contra N' (resultado abaixo). Ressalva: a cobertura tem folga de só 3,9·10⁻⁷ rad (vem do fator 1 + 10⁻⁶ em rZ), que é rigorosa mas justa. |
| R2 | A2 (contagem finita) | **OK** | `r2_contagem.py`, `r2_saida.txt` | Montagem independente: `sharp_operator_builder`, na base real (cos, sen) com convolução 2D, e o truncamento (20, 40) feito por mim nos arrays Z². Há exatamente uma raiz no disco, em 0,4010247311. A condição do autovalor é 4,13 e a próxima raiz está a 0,616 do centro. O resultado é o mesmo com o fundo completo e com Z = 16×64 e 24×96. Conferi a lógica do Rouché escalar (abaixo): esquerda 3,9·10⁻⁴ ≪ direita 0,088. |
| R3 | A3 (lema do espaço simétrico) | **OK (prova abaixo)** | este relatório, seção R3 | Prova por ponto fixo na cauda de Fourier \|m\| ≥ M0 e unicidade nos dois espaços. Hipóteses explícitas: B com cota de Young nos dois pesos e Re λ ≥ 0. |
| R4 | A4 (estrutura do NK) | **OK com ressalva (lacunas L1–L5, inócuas)** | `r4_montagem.py`, `r4_saida.txt` | Derivei as condições do NK bordejado, a simplicidade algébrica, Z2 contra DF(x) − DF(x0) e todas as entradas dos três níveis. Remontei com código meu e reproduzi θ = 0,839826, Z2 = 116,805 e r = 6,1936·10⁻⁴. Achei L1–L4 e mostrei que, corrigidas, não mudam θ nem Y. Conferi L5 numericamente. |
| R5 | A4 (ataque numérico) | **OK** | `r5_blocos.py`, `r5_saida.txt` | Montagem própria (pesos, Q0, seleção de T e dos blocos), com B e J de `signed_operator_L`. Normas verdadeiras (potência com a LU de Hh), comparadas com `T3.json`: ver a tabela da seção R5. Todas ficaram abaixo das cotas. As mais justas: \|\|V\|\| da janela 32..47, 1,69163 contra 1,69180; e N[32..47 ← 16..31], 0,2508 contra 0,2525. |
| R6 | A5 (cotas verificadas) | **OK** | `r6_verified.py`, `r6_saida.txt` | (a) `verified_inverse` e `cota_norma` contra enclosures Arb de 400 bits, em 20 casos, incluindo cond 10⁴ a 10¹⁵, matriz triangular do tipo J_μ pesada, bordejada e com erro_X. Todas as cotas ficaram acima do valor verdadeiro. Nos casos mal condicionados (cond ≥ 10⁸) o método **falha sem dar cota falsa**. Conferi a constante de Rump contra o PDF (Teorema 2.3, cota (3.1), Corolário 2.4). (b) `fundo_L`: nos 81 d e duas paridades, a norma de cada bloco montado ficou abaixo da cota de Young (razão máxima 0,58). (c) O fator 40 ε_ω é justificado com folga: Young ≤ 5,01 ε_ω. |
| R7 | A6 (testemunho) | **OK** | `r7_testemunho.py`, `r7_saida.txt` | Com montagem própria, \|\|Π_F C(λ0) h0\|\| = 0,1717424 (D.json: c_ref ≥ 0,171739). c_J direto por SVD: janela 32..47, 5,73·10⁻⁹ (igual); −15..0, 2,4·10⁻⁸ (cota 3,5·10⁻⁷); 16..31, 1,2·10⁻¹⁵ (cota 1,1·10⁻¹¹). Os c_J são pequenos de verdade: o fundo tem coeficientes ≤ 6·10⁻⁸ no modo d = 21, e as colunas de cauda de \|m\| < 32 só descem de n ≥ 128. A seção finita de \|\|Π_F C1\|\| dá 0,96, abaixo do nC1 = 1,5 usado. A perda total (3,2·10⁻³) é muito menor que c_ref. |
| R8 | A7 (fechamento, espaços) | **OK com ressalva** | `r8_logica.py`, `r8_saida.txt`, seção R8 | D está no retângulo de S3a. O disco de λ* está no disco de Rouché. O centro float(0,401) difere de 0,401 por 2,3·10⁻¹⁷, e isso é coberto pelos próprios tiles. Y_L ⊂ X_K com norma ≤ 1, mas o espaço com sinal X_L não está contido em X_K, então o lema A3 é indispensável. Ressalva: a alegação vale para autovetores; cadeias de Jordan de outras raízes de L em D não ficam controladas (ver R8). |
| R9 | reprodução | **OK** | `r9_saida.txt`, `rig_copia/` | `nk_L3_C.py` e `nk_L3_D.py` rodados com `--dir` numa cópia: C3.json e D.json idênticos aos originais, campo a campo. θ = 0,8398257624513322, r = 6,193640·10⁻⁴, folga 0,1684947. Obs.: `build/s3b/rig/C.txt` está desatualizado (mostra θ 0,839802 e r 6,187·10⁻⁴); o JSON está certo. |

## Lacunas encontradas

| # | Gravidade | Onde | O que falta | Efeito medido | Conserto |
|---|---|---|---|---|---|
| L1 | Baixa (formal; seria real com outros pesos) | `nk_L3_C.py`, T←far; docstring de `nk_L2.py` ("o far só alcança as linhas altas de cada janela") | A frase é falsa nas janelas com \|m\| ≥ 160. As colunas do far de **Fourier** (\|m\| ≥ 200, **todo n**) alcançam, pela banda \|d\| ≤ 40, as linhas **baixas** (n < 240) dessas janelas. Esse acoplamento não é descida radial, então nem `desc` nem ‖V_J P_alto‖·β o cobrem. | Nenhum com os pesos 1 escolhidos: nessas janelas ‖V_J‖ ≤ 1,0887 < max ‖V_J P_alto‖ = 1,2808, e o máximo não muda. Com outros pesos w poderia pesar. | Em T←far, usar max(max_{não borda} w_J‖V_J P_alto‖, max_{\|m\|≥160} w_J‖V_J‖)·β + √Σ(w_J‖V_J‖)²·desc. Validade: as linhas altas das janelas internas e todas as linhas das janelas de borda são conjuntos disjuntos. |
| L2 | Desprezível | `nk_L3_T.py`, grupos de far←T | As linhas do far em 216 ≤ \|m\| < 240 recebem das janelas \|m\| ≥ 176 pelo resto da banda (17 ≤ \|d\| ≤ 40). Não estão em nenhum grupo, e o termo global de `nk_L3_C` só soma RT + μ, não o resto. | fT sobe ≤ 10⁻⁷ (as linhas entram em quadratura); θ inalterado (`r4_saida.txt`). | Incluir os grupos [216, 240) e (−240, −216] com entradas resto·‖Q0 P_J‖, ou ampliar os grupos [200, 216) para [200, 240). |
| L3 | Desprezível | `nk_L3_C.py`, Y | Falta a componente far do resíduo: P_far δB h0 (erro de RT e de μ não respeita a banda), ≤ 40 ε_ω·‖h0‖. Nas janelas não alcançadas pela banda (\|m\| ≥ 72) também falta o termo δB h0. | Y_far/v_f = 4,2·10⁻⁶ e a parte T sobe para 1,03·10⁻⁵, ambos ≪ Y_Z/v0 = 5,44·10⁻⁵. Y = max não muda. | Y = max(Y_Z/v0, (Y_T + max w‖V‖·rt)/v_T, rt/v_f). |
| L4 | Desprezível | `nk_L2_Z.py`, epsZ (Z'←far) | Falta a linha de bordejamento l Q0 P_far. Q0 desce de n ≥ 400 até as linhas de Z (n < 128). | ‖·‖ ≤ descida(128, 400) ≈ 10⁻²⁶. | Somar t_V·descida(NZ, Nc)·‖h0w‖ a epsZ. |
| L5 | Baixa (hipótese não provada, conferida) | `nk_L2_Z.py`, a_extra | Usa ‖Q0 P_{48≤\|m\|<72}‖ = ‖Q_{±48}‖ com a justificativa "‖Q0‖ decresce com \|m\|", que não é provada. | Medi os 48 modos: o máximo é 0,17184 (modo 48) e o valor usado é 0,17199. OK. | Tomar o máximo calculado sobre os 48 modos (é barato), com Q_cert, em vez da hipótese de monotonia. |

Nenhuma dessas lacunas muda o resultado: `r4_saida.txt` mostra o NK com L1–L4 corrigidas, e
L5 conferida, com θ = 0,839826, Y = 5,4398·10⁻⁵ e r = 6,1936·10⁻⁴.

## R1 — Rouché de K (A1): derivação

H(s) = K(s)Q0(c) = I + (B + s − c)Q0(c), com P = Q0(c) fixo e bijetor. Logo H e K têm os mesmos
zeros, com as mesmas multiplicidades. Seja A(s) = diag(Ĥ_ZZ(s), I, I, I). Então:

* **Linhas T1, T2 e far de A⁻¹(H − A):** nelas A = I, então as linhas são H_ik − δ_ik = −E_ik,
  com E = I − diag(V, I, I, I)H. As entradas são as mesmas de N, as do tile.
* **Linha Z:** seja e = N_ZZ ≥ ‖I − V Ĥ_ZZ(s)‖. O `rho` de `parte_Z` usa o erro dH ≥ dB‖Q‖ e
  vale para H_ZZ e para Ĥ_ZZ; o termo rZ‖V Q_ZZ‖ dá conta de s no tile. Então
  Ĥ_ZZ⁻¹ = (VĤ_ZZ)⁻¹V e ‖Ĥ_ZZ⁻¹X‖ ≤ ‖VX‖/(1 − e). Daí:
  * N'_Zk = N_Zk/(1 − e) para k ≠ Z;
  * N'_ZZ = ‖V‖(dB‖Q_ZZ‖ + ε_μ)/(1 − e), porque H_ZZ − Ĥ_ZZ = P_Z(B − B_t)Q0 P_Z mais o termo de
    μ, e dB majora resto + erro de RT (+ arredondamento, que é conservador).

  É exatamente o `Nr` de `circ()`.
* **Perron.** Com Perron em racionais, N'v ≤ θ'v com θ' ≤ 0,8126 < 1, e então
  ‖A⁻¹(H − A)‖ < 1 na norma max_i ‖x_i‖/v_i. O v varia de tile para tile, e isso é permitido:
  Gohberg–Sigal só precisa que A + t(H − A) = A(I + tA⁻¹(H − A)) seja invertível em cada s do
  contorno e para cada t ∈ [0, 1], e isso vale em qualquer norma equivalente.
* **Gohberg–Sigal.** H = I + compacto, porque Q0 é compacto; A = I + posto finito. Então H e A têm
  o mesmo número de zeros no disco.
* **Zeros de A = autovalores de K̂_Z.** K_livre é diagonal em m e triangular superior em n; logo
  Q0 P_Z = P_Z Q0 P_Z = J_Z(c)⁻¹, conferido numericamente: max|P_{não Z} Q0 P_Z| = 0. Daí
  Ĥ_ZZ(s) = (J_Z(s) + B_t,ZZ) J_Z(c)⁻¹ e det Ĥ_ZZ = det K̂_Z(s)/det J_Z(c).
* **Cobertura da circunferência.** Os discos dos tiles também provam K injetivo em si mesmos
  (θ_tile < 1). Com isso cobrem uma coroa em torno da circunferência, que absorve a diferença
  float(0,401) − 0,401 (`r8_saida.txt` iii).

**Medida numérica (tile 8, caixa 40×200, sem o far):** as normas dos blocos de A⁻¹(H − A) nos níveis (Z, T1, T2') ficaram
todas abaixo da N' do JSON, entrada a entrada. Z←Z: 3,7·10⁻⁶ ≤ 1,4·10⁻³; Z←T1: 0,393 ≤ 0,616;
T1←Z: 0,17808 ≤ 0,178125, praticamente exata; T1←T1: 0,452 ≤ 0,490; T1←T2': 0,223 ≤ 0,238;
T2'←T1: 0,037 ≤ 0,108. Com o v do tile, a razão de Perron numérica sem o far é 0,763, contra
θ' = 0,8126. A caixa F' = 40×200 é menor que 200×400, então as entradas com T2' são só
indicativas.

## R2 — Rouché escalar (A2)

t(s) = eᵀM(s)⁻¹e = det K̂_Z(s)/det M(s), pelo cofator da entrada (n+1, n+1). M(s) é invertível no
disco fechado, porque q = ε_c + ρβ = 0,257 < 1. Logo t e det K̂_Z têm os mesmos zeros, com as
mesmas ordens.

A expansão em série e o resto também estão certos: ‖E^j − E_puro^j‖ ≤ (ε_c + ρβ)^j − (ρβ)^j vale
para operadores que não comutam, e a cauda soma q^{K+1}/(1 − q). Como |g| ≥ ρ|τ1| − |τ0| na
circunferência, a desigualdade do código implica |t − g| < |g|. O zero de g está em
c + τ0/τ1, a 2,5·10⁻⁵ do centro.

O código só ganha folga ao incluir resto_banda e RT em dK, embora K̂_Z seja definido com B_t exato
(é conservador). A contagem independente (`r2_saida.txt`) confirma uma raiz simples.

## R3 — Prova do lema A3

**Enunciado.** Seja X_κ o ℓ² de (m, componente, n) com peso κ^{2m}ω2(n), e
Y = X_κ1 ∩ X_{1/κ1}, que tem peso κ1^{2|m|}. Seja L(λ) = O(λ) + B, com O(λ) = ⊕_m R_m(λ),
R_m = J_rad + λ + im/2 e B = Σ_d S_d B_d. Hipóteses:

- **(H1)** Σ_d κ1^{|d|}‖B_d‖ < ∞, com as normas radiais em ω2.
  * Vale para a parte de RefA, com banda \|d\| ≤ 40 e cota de Young de `fundo_L`/`fourier_tail_bound`.
  * Vale para a perturbação δB dos dados de RT, com ‖δB‖_Young ≤ 5,01‖δω‖_RT (R6c).
  * Vale para o termo de μ, que é diagonal em m.
  * Consequência: ‖B‖_{X_κ} ≤ Σ_d κ^{d}‖B_d‖ ≤ Σ_d κ1^{|d|}‖B_d‖ =: b para κ ∈ {κ1, 1/κ1}.
- **(H2)** Re λ ≥ 0. Então ‖R_m⁻¹‖ ≤ c/(\|m\|/2 − \|Im λ\|) para \|m\| > 2\|Im λ\|, com
  c = 1 + 2c_w/(κ2 − 1) (forma fechada; é o `Rinv_uniforme`).

**Afirmação.** Se h ∈ D(L) ⊂ X_κ1 e L(λ)h = f ∈ Y, então h ∈ Y e O h ∈ Y. Em particular,
autovetores (f = 0) e, por indução, cadeias de Jordan estão em Y. Nas cadeias,
L(λ)h_k = −h_{k−1}, porque L'(λ) = I.

**Prova.**

1. Escolha M0 > 2\|Im λ\| com q := b·sup_{\|m\|≥M0}‖R_m⁻¹‖ < 1. Seja P_> a projeção em
   \|m\| ≥ M0 e P_< = I − P_>.
2. Da equação vem P_>h = g − T P_>h, com T = R_>⁻¹P_>BP_> e g = R_>⁻¹P_>(f − BP_<h).
3. ‖T‖ ≤ q < 1 em P_>X_κ1 e também em P_>X_{1/κ1}.
4. P_<h tem um número finito de modos, e em cada modo a norma radial é a mesma nos dois espaços.
   Logo P_<h ∈ X_κ1 ∩ X_{1/κ1}, B P_<h ∈ Y por (H1), e portanto g ∈ Y.
5. A série de Neumann Σ(−T)^k g converge nos dois espaços. Cada coordenada é um funcional
   contínuo nos dois, e T age pela mesma fórmula coordenada a coordenada. Então os dois limites
   coincidem, e esse limite é a única solução em X_κ1. Logo P_>h ∈ Y.
6. Por fim, Oh = f − Bh ∈ Y. ∎

**Consequências.**

- ker_{X_κ1}L(λ) = ker_Y L(λ), e os sistemas de cadeias de Jordan coincidem.
- L(λ)O⁻¹ = I + BO⁻¹ é Fredholm de índice 0 nos dois espaços: O⁻¹ é compacto, porque
  ‖R_m⁻¹P_{n≥N}‖ → 0 e ‖R_m⁻¹‖ → 0.
- Logo as raízes e as multiplicidades algébricas (soma dos comprimentos das cadeias) coincidem.

O mesmo argumento vale para K, com κ1 = 129/128. **Limite de validade:** Re λ ≥ 0, que cobre D
(Re λ > 0,27). Para Re λ < 0 seria preciso outra cota de ‖R_m⁻¹‖.

## R4 — Derivações do NK (A4)

**Aplicação.** x = (y, λ), h = Q0 y, e

F(x) = ((I + BQ0 + (λ − λ0)Q0)y, ⟨h0w, Q0y⟩ − 1), com Q0y0 = h0w (y0 = (J + λ0)h0w, suporte em Z).

**Derivadas.**

- DF(x0) = [[I + BQ0, h0], [lQ0, 0]].
- DF(x) − DF(x0) = [(λ − λ0)Q0·dy + dλ·Q0(y − y0); 0], exato.
- F(x) − F(x0) − DF(x0)(x − x0) = ((λ − λ0)Q0(y − y0), 0).

**Condições.** Suponha ‖A b(u, w)‖ ≤ Z2‖u‖‖w‖ para o bilinear b. Então:

- ‖DT(x)‖ ≤ θ + Z2‖x − x0‖;
- ‖T(x) − x0‖ ≤ Y + θr + Z2r²/2.

As condições do código, Y + θr + Z2r² ≤ r e θ + 2Z2r < 1, são mais fortes que as necessárias
(r² em vez de r²/2 e 2Z2 em vez de Z2).

**A é invertível.** A = diag(Mh⁻¹, V_J, I), e Mh e Hh_J têm Loewner t²PP* − I ⪰ 0 verificado.
Logo ‖I − A DF(x*)‖ < 1 implica que DF(x*) é invertível.

**Simplicidade.** Seja G[dh, dλ] = (L(λ*)dh + dλh*, l(dh)), invertível.

- Se L(λ*)k = 0 e l(k) = 0, então k = 0; logo a multiplicidade geométrica é ≤ 1.
- Se existisse uma cadeia L(λ*)g = −h*, então G[g − l(g)h*, 1] = 0, contradição. Isso usa
  l(h*) = 1.
- Com índice 0, λ* é algebricamente simples.

**Z2.** Seja u = (λ − λ0)dy + dλ(y − y0). Pela estrutura de Q0 (diagonal em m, triangular
superior em n), P_TQ0P_Z = 0 e P_far Q0 P_{Z,T} = 0. Então:

- ‖(Ab)_{Z'}‖/v0 ≤ t_V(‖Q_ZZ‖‖u_Z‖ + ‖P_ZQ0P_T‖‖u_T‖₂ + Q_f‖u_f‖)/v0;
- ‖u_Z‖ ≤ 2v0²ρδ, ‖u_T‖_W ≤ 2v0v_Tρδ, ‖u_f‖ ≤ 2v0v_fρδ;
- ‖u_T‖₂ ≤ ‖u_T‖_W/min w;
- ‖P_ZQ0 u_T‖ ≤ max_{J∩Z≠∅}‖Q0P_J‖·‖u_T‖₂, porque os modos das janelas são disjuntos e Q0 é
  diagonal em m.

Isso dá exatamente a fórmula de `nk_L3_C.py`, e as linhas T e far saem do mesmo modo.

**Três níveis.**

- **T←T** = ‖WNW⁻¹‖₂. O vetor das normas por janela é dominado por N aplicado ao vetor ‖x_J‖.
- **T←Z'** = √Σ(w_J TZ_J)². A coluna λ de DF é (h0, 0), com P_T h0 = 0.
- **Z'←T** = a/min_{SZ} w, com SZ = janelas com min\|m\| < 72. Conferi que 64..79 ∈ SZ. Em a, a
  matriz Y = [P_Z B Q0 P_T; lQ0P_T] aplica B à coluna inteira de Q0, inclusive às linhas de Z, e a
  linha λ usa a descida de Q0 até Z. As colunas 48 ≤ \|m\| < 72 entram em quadratura via a_extra
  (L5).
- **T←far:** L1.
- **far←T:** L2.
- **far←Z'** = pert·‖Q_ZZ‖: só a perturbação global alcança o far a partir de Z (72 < 200 e
  228 < 400).
- **Z'←far** = epsZ: L4.
- **Perturbação global** (RT e μ, que não respeitam a banda): `nk_L3_C` soma-a de novo em todas as
  entradas, como cota de bloco inteiro (max w‖V‖·eg·‖Q0‖/min w). Isso é válido.
- **Linha de Z'** = [ρ_Z, Z'←T, epsZ]·v/v0 = 0,8398, igual às outras duas linhas (vetor de Perron
  equilibrado).

## R5 — Tabela das normas verdadeiras (PC do laboratório)

| Bloco | verdadeiro (potência) | cota T3 |
|---|---|---|
| ‖Hh⁻¹‖, 1..15 / −15..0 | 1,20389 / 1,20382 | 1,20728 / 1,20777 |
| ‖Hh⁻¹‖, 32..47 / −47..−32 | 1,69163 / 1,69163 | 1,69180 / 1,69180 |
| ‖Hh⁻¹‖, 64..79 | 1,33514 | 1,33530 |
| ‖Hh⁻¹P_alto‖, 64..79 (o máximo que entra em T←far) | 1,27552 | 1,28079 |
| ‖Hh⁻¹P_alto‖, −47..−32 | 1,17285 | 1,17900 |
| N[1..15 ← −15..0], N[1..15 ← 16..31] | 0,06722, 0,09412 | 0,07256, 0,09534 |
| N[−15..0 ← −31..−16], N[−15..0 ← 1..15] | 0,09855, 0,07182 | 0,09976, 0,07677 |
| N[32..47 ← 16..31], N[32..47 ← 48..63] | 0,25080, 0,24858 | 0,25252, 0,24927 |
| N[−47..−32 ← −63..−48] | 0,26197 | 0,26266 |
| N[64..79 ← 48..63] | 0,19532 | 0,19587 |
| N[32..47 ← 1..15], N[32..47 ← −15..0] (só \|d\| ≥ 17) | 1,1·10⁻⁵, 3,0·10⁻⁹ | 0,0065, 0,0075 |
| TZ: 32..47, −47..−32, −15..0 | 0,35968, 0,34355, 0,10858 | 0,36720, 0,35106, 0,11393 |

Dentro de uma janela de 16 modos, \|d\| ≤ 15; por isso ‖I − Hh⁻¹H_II‖ com banda 17..40 é
identicamente 0, e a diagonal do certificado (resto + RT + μ) é conservadora.

## R6 — Contabilidade de erros (resumo)

- **Loewner.** ‖H⁻¹X‖ ≤ t ⟺ XX* ⪯ t²HH*. Os erros entram na margem de Rump:
  - P: γ_k‖H‖_F², com γ complexo = γ_{4(k+1)};
  - t²P − S: 2u(t²‖P‖_F + ‖S‖_F);
  - simetrização: u‖G‖_F.
- **Realificação.** Conferi o sinal do bloco −Im.
- **Corolário 2.4.** É aplicado a A = Ar − extra·I, com diagonal real (não precisa ser
  representável). Isso é válido, e o 2u·max a_ii cobre o arredondamento de a_ii − c.
- **Arredondamento de t².** A cota sai para √fl(t²) em vez de t, uma diferença relativa de 10⁻¹⁶,
  coberta pelo FATOR_PESO.
- **`fundo_L`.** A estrutura M(d)[cout, (cin, par)] = 2·fator·l1(combo) majora os blocos montados
  (razão ≤ 0,58).
- **40 ε_ω.** A norma de RT é ‖δω‖ = Σ_j Σ_{m,n∈ℤ} κ1^{\|m\|}κ2^{\|n\|}‖v_mn‖_C, com ‖z‖_C ≥ \|z\|,
  somada nos 4 campos, e RT dão ‖δω‖ ≤ 2⁻²⁵ + 2⁻²⁷⁷. A cota de Young é subaditiva, e seu máximo
  na bola ℓ¹ é atingido num coeficiente só. Então ‖δB‖ ≤ max_j‖C_j‖₂·ε_ω = 5,01ε_ω.
- **ε_μ.** |RefBμ| ≤ 43·2⁻⁴⁰ em RT; `DELTA_MU` usa os 80 dígitos de μ. Não reconferi o valor de 80
  dígitos.

## R8 — Lógica final e espaços (A7)

1. **K em X_K com sinal (129/128, 9/8).**
   - S3a dá K injetivo em D \ {\|s − 0,401\| < 0,05}, e D está no retângulo.
   - Os tiles de A1 cobrem uma coroa em torno da circunferência, que inclui a diferença de centro.
   - Na circunferência, H é invertível.
   - Dentro do disco, A1 e A2 dão exatamente um zero de H, com multiplicidade 1.
   - Conclusão: K tem exatamente uma raiz s_K em D, algebricamente simples.
2. **L em X_L com sinal (65/64, 5/4).** O NK dá λ*, com |λ* − c| ≤ 3,05·10⁻⁴ (logo λ* ∈ D) e
   algebricamente simples.
3. **Passagem para o espaço de K.** X_L com sinal não está contido em X_K: em m → −∞ a razão de
   pesos é (130/129)^{2\|m\|}. Por isso o lema A3 é necessário. Com ele, h* ∈ Y_L ⊂ X_K (norma ≤ 1).
   - C tem banda 40 em m e é de primeira ordem em ξ. A folga radial (5/4 → 9/8) absorve a derivada,
     então C(λ*)h* ∈ X_K e K C = M L vale em h* ∈ D(L) ∩ Y_L.
   - A6 dá C(λ*)h* ≠ 0, e K(λ*)C(λ*)h* = M L(λ*)h* = 0.
   - Logo λ* é raiz de K em D, isto é, λ* = s_K.
4. **Unicidade.** Seja μ ∈ D raiz de L com um autovetor h tal que C(μ)h ≠ 0. Pelo lema A3,
   h ∈ Y_L; pela identidade, μ é raiz de K em D; logo μ = s_K.

**Ressalvas.**

- A alegação é sobre **autovetores**. Considere uma raiz μ ≠ s_K de L em D com autovetor que
  satisfaz as constraints, mas com um vetor generalizado g que as viola. Então
  K(μ)C(μ)g = −M h, que não precisa ser 0, e o argumento nada diz sobre μ. Se o uso a jusante
  contar multiplicidades algébricas de L em D, isso precisa ser tratado à parte.
- A pertinência de C(λ*)h* ao domínio de K e a validade da identidade no espaço certo foram
  conferidas por mim só qualitativamente, com a folga de pesos. Não li
  `SHARP_LINEARIZED_IDENTITY.md` a fundo.

## O que não verifiquei

- **`NORMA_BL` = 3,964043**, a cota do símbolo em Arb de `certify_symbol_numerical_range.py`.
  Entra em β = 0,4951, que pesa ~60% na linha do far e ~28% na linha T. Um erro de ≳0,2% passaria
  da folga de θ. Conferi só que o método do símbolo inclui a divisão por ξ (μSΓ1) e que a norma
  de Young (18,07) está acima. **É a principal entrada não reconferida deste certificado.**
- **Z'←T (a = 0,37889).** Não medi a norma verdadeira de Mh⁻¹[P_ZBQ0P_T; lQ0P_T]: exigiria uma
  potência implícita com a LU de Mh (14 017) e ~5·10⁴ colunas. Revisei a montagem de Y linha a
  linha, e os mesmos blocos (B, Q0, pesos) foram validados em R5 pelo lado T←Z'.
- **Formas fechadas `cauda_Rinv`, `descida` e `descida_divxi` para κ2 = 5/4.** São a regressão das
  versões de K auditadas em S3a. Só conferi Qfar por seção finita: 0,1167 ≤ 0,1249.
- **Validação de `signed_operator_L` e `signed_constraints` contra o exportador de RT.** Não rodei
  `--validar`, porque dependem de `build/spectrum/rt-A-*.dat`. Aceitei a validação do autor.
- **Dados de T1/T2 de `rouche_K_circ.json`** (partes T1, T2 e far de K). Estão cobertos pela
  revisão de S3a e por R1 só numericamente, numa caixa menor.
- **O valor de 80 dígitos de μ e a cota |μ − μ_80| ≤ 10⁻⁸⁰.**
