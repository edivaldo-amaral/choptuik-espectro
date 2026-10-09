# Revisão independente do NK de L em λ0 = 0,7331796596606924 (build/nk733)

Revisor: agente de contexto limpo (Claude), 08–09/10/2026. Relatório em construção; cada seção é gravada à
medida que o artefato correspondente fica pronto.

## 0. Snapshot (início)

- `COMMIT.txt`: d85f6490fa033870062368cdc04a0a317a9df9b8, igual ao `HEAD` do repositório no início.
- `sha256sum -c .snapshots/2026-10-09-revisao-nk733/MANIFEST.sha256`: 274/274 arquivos SUCESSO, status 0
  (saída completa em `n1_manifesto_snapshot_inicio.txt`).

## N1. Integridade

Artefatos: `n1_integridade.py` → `n1_integridade.out`; `n1_manifesto_nk733.txt`.

- `build/nk733/MANIFEST.sha256`: 18/18 SUCESSO (`n1_manifesto_nk733.txt`).
- Um só λ0: `Z.json`, `Q.json`, `C3.json` têm `0.7331796596606924`; `T3.json` tem `[0.7331796596606924, 0.0]`;
  `autovalor.txt` (Newton bordejado) termina em `LAMBDA 0.7331796596606924` (parte imaginária 2,3e-16 descartada).
- Caixas: `Z.json` Z = 32×128 (n = 14016), F = 200×400, banda 16, `eps_fundo_L` = 4,7e-8 em Z e em T3; `T3.json`
  `folga_far` 60, `rouche` false.
- Janelas: reconstruí as janelas de W = 16 **sem importar o código do autor** (central |m| < 16 dividida em
  −15..0 e 1..15; depois [16k, 16k+16) cortadas em 200, positivas e refletidas). As 26 de `T3.json` são exatamente
  essas, cobrem |m| < 200 sem sobreposição; `Nf` é 27×26 (25 grupos de F + 2 grupos do far de Fourier); `SZ` =
  janelas com min|m| < 72 (10 janelas), como deve.
- Constantes partilhadas: `dB_L` e `eps_mu` idênticos em Z e T3; `T3.json` traz `nQZ` e `Qfar` iguais aos de
  `Z.json` (o T3 leu este Z.json); `pert = dB_L.total + eps_mu`; `beta = (3,964043 + rt)·Qfar + eps_mu` bit a bit.
- `Q.json`: q_ZT, q_TT, q_Zf e `nQx_extra` batem com `q.txt`; `nQx_extra` é o máximo sobre 48 ≤ |m| < 72
  (= MZ + banda .. MZ + 40), como pede a correção L5 de `nk_L3_C.py`. `Z.json` tem `tQ` (de `tQ.txt`).
- Junção: `T3.json` é a junção **bit a bit** de `T3.parcial.json` (19 janelas: −15..0 a 192..199 e −31..−16 a
  −95..−80) e `T3.joao.json` (7 janelas: −111..−96 a −199..−192), com ρ_J na diagonal de N; os dois
  checkpoints são disjuntos e somam as 26. Os logs `etapaB.txt` e `etapaB_joao.txt` batem com o JSON nas 26
  janelas (‖V‖, ‖V P_alto‖, ρ, Y).
- `h0w.npy`: 14016 entradas complexas, ‖h0w‖ = 1 − 4,4e-16. Nenhum JSON grava um hash de h0w; a prova de que Z
  e T3 usaram **este** h0w é a reexecução de N4 (Y_Z e os Y das janelas refeitos com este arquivo batem com os
  gravados; ver N4).
- Logs × JSON: `etapaA.txt` (t_V, ρ_Z, Y_Z, a, ε_far), `tQ.txt` e `q.txt` batem com `Z.json`/`Q.json`.

Observação (sem gravidade): 7 das 26 janelas foram calculadas em outra máquina ("joao") e juntadas; os arquivos
não registram a máquina nem o ambiente (`CHOPTUIK_EPS_FUNDO_L`, λ0) dessa rodada. Por isso uma das janelas
refeitas em N4 é de lá (−191..−176).

**N1: OK.**

## N2. Montagem (racionais exatos, código próprio)

Artefatos: `n2_montagem.py` → `n2_montagem.out`, `n2_montagem.json`; reexecução do autor em cópia:
`n2_verifica/c3_rerun/` e `n2_verifica_T2.out` (`verifica_T2.py --so nk_0733` rodado numa cópia do projeto dentro
desta pasta, para não escrever em `build/reproducao/`).

Método do revisor (`n2_montagem.py`, sem importar `scripts/`): cada double dos JSON entra como o racional exato;
as normas espectrais ‖W N W⁻¹‖₂ (26×26) e ‖Nf W⁻¹‖₂ (27×26) saem de uma cota **diferente** da do autor
(Collatz–Wielandt exato em AᵀA ≥ 0: ‖A‖₂² ≤ max_i (AᵀA u)_i/u_i para u > 0; o autor usa Rump em ponto
flutuante); raízes quadradas com verificação s² ≥ x. Monto as 9 entradas com a derivação de N3, contando a
perturbação global (RT + μ) **uma** vez por bloco e sem os fatores de folga ARRED.

| quantidade | revisor (R) | gravado (C3.json) |
|---|---|---|
| ‖N‖₂ (T←T sem a perturbação global) | 0,37851409 (float 0,37851409) | 0,3787492 − 1,6e-7 |
| ‖Nf‖₂ | 0,09387085 (float 0,09387085) | ≈ 0,093964 |
| θ (Perron 3×3, w = 1) | 0,805500999 | 0,8057209580505 |
| Z2 | 10,325314 | 10,3232499 |
| Y | 4,187181e-6 (Z′ 8,96e-7; T 4,1872e-6; far 1,63e-7) | 4,1871805e-6 |
| discriminante (1−θ)² − 4 Z2 Y | 3,76569e-2 | 3,75714e-2 |
| r | 2,1552691e-5 | 2,1577165e-5 |
| \|λ* − λ0\| ≤ v0 r | 1,0153471e-5 | 1,0162255e-5 |

Conferências exatas:

- **A M gravada domina a do revisor entrada a entrada** (racionais): folga relativa entre 1e-9 e 9,9e-4; as
  maiores (6,2e-4 em T←T e 9,9e-4 em far←T) vêm da cota de Rump do autor para ‖N‖₂ e ‖Nf‖₂, mais frouxa que
  Collatz–Wielandt, e das somas duplicadas da perturbação global. Ambas a favor da segurança.
- θ gravado = max_i (M v)_i/v_i · (1 + 1e-9) **exatamente**, com a M e o v gravados.
- Com os números gravados: discriminante 3,757e-2 > 0; r gravado ≥ menor raiz (2,1577144e-5);
  Y + θr + Z2 r² ≤ r e θ + 2 Z2 r = 0,806166 < 1 (racionais); `erro_lambda` = v0·r exatamente; a unicidade vale
  até ρ < (1 − θ)/(2 Z2) = 9,41e-3.
- A montagem do revisor dá θ **menor** (0,805501) e r **menor**. Seu Z2 sai um pouco **maior** (10,3253 contra
  10,3232) só porque cada montagem usa o próprio vetor de Perron: v0 = 0,471100 (revisor) contra 0,470973
  (gravado), e Z2 = 2(tQ v0 + t_V q_ZT vT + …) cresce com v0; as duas são normas válidas; a do autor reescrita (fórmulas de `nk_L3_C.py`
  com as minhas normas) dá 0,805501316.
- Reexecução do código do autor numa cópia (`nk_L3_C.py --dir n2_verifica/c3_rerun`): o `C3.json` sai
  **idêntico bit a bit** ao gravado e a saída é idêntica a `etapaC.txt`. `verifica_T2.py --so nk_0733` (cópia do
  projeto em `n2_verifica/proj`): `OK`, `iguais: theta, r, erro_lambda = true`. Observação: `check_nk` compara só
  θ, r e erro_lambda; a comparação do `C3.json` inteiro (M, v, Z2, Y, w) foi feita aqui.

**N2: OK.** A montagem gravada é reproduzida e é dominada, entrada a entrada, por uma montagem independente.

## N3. Fórmulas (derivação do revisor)

Artefatos de apoio: `n3_triangular.py/.out` (estrutura de Q0), `n3_banda_n.py/.out` (largura radial de B),
`n3_simbolo_L.py/.out` e `n3_simbolo_L_arb.out` (a constante ‖B_L‖ = 3,964043).

**Montagem do sistema.** x = (y, λ), Q0 = (J_{μ_A} + λ0)⁻¹ (λ0 = o double, μ_A fixo), h = Q0 y,
F(x) = (L(λ) Q0 y, ℓ(Q0 y) − 1), L(λ) = J_μ + B + λ com μ e B **verdadeiros**, ℓ(h) = ⟨h0w, h⟩ (coordenadas
pesadas). y0 := Q0⁻¹ h0w, logo Q0 y0 = h0w exatamente. Como L(λ) = L(0) + λ:

- DF(x)[η, ν] = (L(λ) Q0 η + ν Q0 y, ℓ(Q0 η)), e L(λ0) Q0 = I + B Q0 + δμ J1 Q0 (δμ = μ − μ_A);
- DF(x) − DF(x0) = ((λ − λ0) Q0 η + ν Q0 (y − y0), 0);
- D²F[(η, ν), (η′, ν′)] = (ν′ Q0 η + ν Q0 η′, 0), constante; a linha de ℓ é linear e não tem 2ª derivada.

**Estrutura usada** (checada em `n3_triangular.out`, 10 modos, N = 128 e 400): J é diagonal em m e "triangular
superior em n" (J[i, j] = 0 se n_i > n_j), logo Q0 também; o bloco n < 128 de Q0(400) coincide com o inverso do
bloco truncado (erro relativo ≤ 1e-15). Daí P_T Q0 P_Z = 0, P_far Q0 P_{Z,T} = 0, e Q_ZZ é o Q0 exato restrito.
Em `n3_banda_n.out`: nos blocos montados de B, n_out − n_in ≤ 100 (o código supõe ≤ BAND_N + 1 = 101), o que
justifica truncar as linhas de B h0 em n < 229 e as do far em n < 501.

**Z2.** Na norma ‖x‖_v = max(‖x_Z′‖₂/v0, ‖x_T‖_W/vT, ‖x_f‖₂/vf), para ‖a‖_v, ‖b‖_v ≤ 1 vale |ν| ≤ v0, logo
‖A D²F[a, b]‖_v ≤ 2 v0 · sup_{‖η‖_v ≤ 1} ‖A (Q0 η, 0)‖_v, e por nível (A = diag(Mh⁻¹, V_J, I)):

- Z′: ‖Mh⁻¹ [P_Z Q0 η; 0]‖ ≤ tQ ‖η_Z‖ + t_V (q_ZT ‖η_T‖₂ + q_Zf ‖η_f‖) ≤ tQ v0 + t_V (q_ZT vT/w_min + q_Zf vf);
- T: ‖(V_J P_J Q0 η)_J‖_W ≤ max_J w_J‖V_J‖ · (q_TT ‖η_T‖₂ + ‖Q0 P_far‖ ‖η_f‖) (P_T Q0 P_Z = 0; as peças
  P_J Q0 η são ortogonais porque Q0 é diagonal em m);
- far: ‖P_f Q0 P_f η_f‖ ≤ Qfar vf.

Isso é exatamente `Z2 = 2 v0 max((tQ v0 + tV (qZT vT + qZf vf))/v0, maxwV (qTT vT + Qf vf)/vT, Qf)` de
`nk_L3_C.py` (com q_ZT, q_TT já divididos por w_min). Esse Z2 limita ‖A D²F‖ **inteiro**; a forma de NK do artigo
pede só ½‖A D²F‖ ≤ Z2, então há um fator 2 de folga, como o artigo declara. Correto.

**As 9 entradas de θ** (E = I − A DF(x0)). Conferi cada uma contra a derivação; todas cobrem o que devem:
Z′←Z′ = t_V‖M − Mh‖; Z′←T = a/min_{SZ} w (a por Loewner com banda completa |d| ≤ 40 nas colunas |m| < 48, mais
a_extra em quadratura para 48 ≤ |m| < 72, que só se ligam a Z por |d| > 16; o máximo de ‖Q0‖ nessa faixa vem de
`Q.json`); Z′←far = ε_Z + t_V q_Zf (linha de ℓ); T←Z′ = ‖(w_J TZ_J)‖₂; T←T = ‖W N W⁻¹‖₂ (diagonal ρ_J =
‖V_J‖·erro(Hh_J)); T←far = max_J w_J‖V_J P_alto‖ β (‖V_J‖ nas janelas com |m| ≥ 160, alcançadas pelo far de Fourier
em todas as linhas) + ‖(w_J‖V_J‖)‖₂ · desc; far←Z′ = pert·‖Q_ZZ‖; far←T = ‖Nf W⁻¹‖₂ + resto·‖Q0‖/w_min (linhas do
far em 216 ≤ |m| < 240, que nenhum grupo de Nf cobre); far←far = β. A perturbação **global** (erro de RT do fundo +
ε_μ), que não respeita a banda, é somada de novo em toda entrada com ‖Q0 P_T‖ ≤ max nQ (inclusive onde já estava,
o que é conservador). As descidas `closed_form.descida*` valem para Re σ ≥ 0; aqui σ = λ0 + i m/2 com
λ0 = 0,733 > 0.

**Constante ‖B_L‖ = 3,964043 (`nk_L.NORMA_BL`).** Entra em β, desc, ε_Z e nos erros. Refiz a folga do certificado
do símbolo em Arb (`n3_simbolo_L.out`): ‖B_L(RefA)‖ ≤ 1937/500 + 0,0900420665 = 3,9640420665 ≤ 3,964043 (folga
9,3e-7). O código soma o erro de RT à parte (β = (NORMA_BL + rt)·Qfar; nBL = NORMA_BL + dB_L.total), então está
coerente. Observação: com a correção genérica 40 ε_ω embutida, o símbolo daria 3,9640432586 > 3,964043; a constante
é, portanto, a do fundo **dos dados**, e só é válida porque o rt (4,7e-8) é somado separadamente. O comentário em
`nk_L.py` ("símbolo com folga em Arb (24/09)") não diz isso, e `l_real_toro.py` cita outra cota (3,9641333, com a
folga em float ×1,001 e 40 ε_ω). Ver achado I-1.

**Y.** F(x0) = (L(λ0) h0, ‖h0w‖² − 1). L(λ0) h0 = [L_calc(λ0) h0] + δμ (J1 + SΓ1) h0 + (B − B_RefA) h0 + (arredondamento
da montagem de B). A parte calculada vive em |m| < 72, n < 229 (Z e as 10 janelas de SZ). O termo de μ fica em Z
(J1 diagonal em m e triangular; SΓ1 em d = 0, só desce): ≤ Δμ (‖J1 h0‖ + 2 ‖divξ‖). O erro de RT não respeita
banda: ≤ rt ‖h0‖ = rt em **todos** os níveis. O código tem os três: Y_Z = ‖x_s‖ + t_V (resíduo da solução + dr), com
dr = produtos + U‖r‖ + (arred + rt) + termo de μ + 4nU; Y_J = Loewner(‖V_J x_J‖) + ‖V_J‖ e_J com e_J contendo
(arred + rt) e o termo de μ; `nk_L3_C.py` soma ainda max_J w_J‖V_J‖·rt em T (cobre as janelas fora de SZ) e
rt/vf no far. Y = max(Y_Z/v0, Y_T/vT, rt/vf). Correto (com dupla contagem conservadora de rt nas janelas de SZ).

**r e erro de λ.** r é a menor raiz de Z2 r² − (1 − θ) r + Y = 0; o código usa r em float ×(1 + 1e-6) e verifica
em racionais Y + θr + Z2 r² ≤ r e θ + 2 Z2 r < 1 (refeito em N2). Com ‖I − A DF(x)‖ ≤ θ + Z2 ‖x − x0‖_v < 1 na bola
e A injetivo (Mh e Hh_J invertíveis pelas próprias cotas de Loewner), sai zero único e DF(x*) invertível; pelo
Lema 6.1, λ* é raiz algebricamente simples. O λ é a última coordenada do bloco Z′, cuja norma é ‖·‖₂ sem peso, logo
|λ* − λ0| ≤ ‖(x* − x0)_{Z′}‖₂ ≤ v0 r. `erro_lambda = v0·r` (exato em N2). **Correto.**

**N3: OK**, com a imprecisão I-1 (documentação da constante ‖B_L‖).

## N5. Enunciado do artigo

Artefatos: `n5_artigo.py` → `n5_artigo.out`, `n5_artigo.json`. Li só a linha de s* da tabela da Proposição 6.2
(`paper/sec/roots.tex`, l. 59) e a linha de s* do Teorema Principal (`paper/sec/main_theorem.tex`, l. 44), e
comparei com os racionais de `C3.json`, tomando λ0 como o double exato
(0,73317965966069242877978…).

| | certificado (racional) | artigo | artigo ≥ certificado | folga relativa |
|---|---|---|---|---|
| θ | 0,80572095805 | 0,805721 | sim | 5,2e-8 |
| Z2 | 10,3232499 | 10,33 | sim | 6,5e-4 |
| Y | 4,1871805e-6 | 4,19·10⁻⁶ | sim | 6,7e-4 |
| r | 2,1577165e-5 | 2,16·10⁻⁵ | sim | 1,1e-3 |
| \|λ* − λ0\| | 1,0162255e-5 | 1,0163·10⁻⁵ | sim | 7,3e-5 |

Intervalo do Teorema Principal: λ0 ∓ 1,0163e-5 = [0,7331694966607; 0,7331898226607]; arredondado para fora em 8
casas dá exatamente [0,73316949; 0,73318983], que contém λ0 ∓ erro certificado com folgas de 7,4e-9 (embaixo) e
8,1e-9 (em cima). Tudo arredondado para o lado seguro.

Cosmético (C-2): a tabela escreve λ0 = 0,7331796596606924, que é o `repr` do double; o centro real é o double, que
difere desse decimal em 2,9e-17. Irrelevante para o intervalo (folgas de 7e-9), mas o texto diz que os valores são
"rounded outward from the exact rational values"; uma nota "λ0 é o double mais próximo" fecharia o ponto.

**N5: OK.**

## N4. Dados de entrada, por amostragem (laboratorio2)

Ambiente: cópia do snapshot em `~/revisao_nk733/proj` no laboratorio2 (`scripts/*.py`, `RefA.dat` e os dados de
`build/nk733`); `teste_local/sha_local.txt` × `teste_local/sha_lab2.txt`: 192 hashes idênticos. numpy 2.5.3 e
scipy 1.18.1 nas duas máquinas. `CHOPTUIK_EPS_FUNDO_L=4.7e-08`, `OPENBLAS_NUM_THREADS=4` no laboratorio2.
Cadeia: `n4_lab2.sh` (e `n4_lab2_tQ2.sh`); logs copiados em `lab2/`.

### N4.1 Janelas: reexecução de `nk_L3_T.py`

Escolhi −15..0 (central: tem acoplamento com Z′, resíduo Y e as duas vizinhas centrais), 32..47 (a mais pesada:
dá max_J ‖V_J‖ = 1,6315, que entra em Z2 e em T←far, e o maior Y_J = 3,84e-6, que domina Y) e −191..−176
(periférica, de borda, uma das 7 calculadas na outra máquina). Rodei `nk_L3_T.py --so-janelas` com um checkpoint
"falso" que contém as outras 23 janelas (`n4_fake.py`), de modo que só as três foram calculadas.

Resultado (`n4_compara_janelas.py` → `n4_compara_janelas.out`): **todas as quantidades das três janelas saem
idênticas bit a bit às de `T3.json`** — ‖V_J‖, ‖V_J P_alto‖, ρ_J, TZ_J, Y_J e todos os acoplamentos não nulos
(11 + 11 + 9 números). Como Y_J depende de h0w, isto confirma que o T3 foi feito com o `h0w.npy` gravado; e a
janela −191..−176 confirma que a rodada da outra máquina foi feita com os mesmos λ0, fundo e código.

### N4.2 Janelas: controles independentes em ponto flutuante

`n4_janelas_float.py` (laboratorio2) monta Hh_J com a mesma montagem e estima, por **ARPACK** (svds com a LU de
Hh_J — método diferente da potência e do Loewner do autor), ‖V_J‖, ‖V_J P_alto‖, ‖V_J X_JI‖, ‖V_J X_JZ‖ e
‖V_J x_res‖ (solução direta); depois aplica o critério de Loewner/Rump do autor em t = 0,999 × estimativa (controle
negativo). Saídas: `lab2/n4_janelas_float.json`; tabela em `n4_janelas_float_analise.py` → `n4_janelas_float_analise.out`.

- Todas as cotas rigorosas ficam **acima** das estimativas. A diferença rigoroso − float é, em cada acoplamento,
  ≈ ‖V_J‖·pert·‖Q0 P_c‖ (a perturbação global — resto da banda |d| > 16, erro de RT, μ — somada como cota), mais
  ≤ 0,2% do Loewner (passos geométricos de 1,002). Ex.: 32..47 ← 16..31: 0,234970 (float) → 0,236427 (rigoroso),
  diferença 1,458e-3 contra 1,434e-3 da perturbação.
- O controle negativo **falha** nas três janelas (t = 0,999·‖V_J‖_float não é aceito), como deve: o teste de
  definição positiva não é vazio.
- Folga não aproveitada (sem gravidade): na janela 32..47, Y_J = 3,84e-6 contra ‖V_J x‖ = 2,30e-6 em float. O
  Loewner do resíduo parte de t0 = ‖V_J‖·‖x‖ e aceita na primeira tentativa, então Y_J é, na prática, o produto de
  normas. É rigoroso; apertá-lo reduziria Y (que esta janela domina) e r em ~40%. Não afeta o veredito.

### N4.3 Nível Z: Mh, Y_Z, ρ_Z e a cota de Loewner de ‖Mh⁻¹‖

`n4_Z.py` (laboratorio2; log e JSON em `lab2/n4_Z.log`, `lab2/n4_Z.json`) remonta Mh = [[H_ZZ, h0w], [ℓ Q0, 0]]
com o `h0w.npy` gravado, refaz o resíduo e testa a cota.

- ‖M − Mh‖ (dM) = 9,506498120584314e-8, ‖Q_ZZ‖ = 1,1842320033535474, erro de Q_ZZ e ‖J1 h0‖: **idênticos bit a
  bit** aos de `Z.json`.
- Resíduo: ‖r_Z‖ = 1,5e-15 (h0w é autovetor da truncagem com precisão de máquina); ‖Mh⁻¹ r‖ ≈ 1,843e-15;
  dr = 4,7703e-8, dominado pelo erro de RT do fundo (4,700e-8), μ 6,2e-10, produtos 7,2e-11.
  **Y_Z = 4,22216912396724e-7 e ρ_Z = 8,414222057094064e-7, idênticos bit a bit a `Z.json`.**
- ‖Mh⁻¹‖ por ARPACK (3 maiores valores singulares): 8,850137 (7,4229; 6,6199). Cota gravada t_V = 8,85102:
  razão 1,0001.
- Loewner t²·Mh Mh* − I ≥ 0 no t gravado (t_V/(1 + 1e-12) = 8,851022): **verificado**. Controle negativo em
  t = 0,999 × 8,850137 = 8,841286: **rejeitado**, como deve.

### N4.4 Q.json e constantes que dependem de λ0 (notebook)

- `n4_Q.py` → `n4_Q.out`: normas restritas de Q0 pesado recalculadas por SVD (inversão e normas do revisor) em
  todos os modos |m| < 200: q_ZT 0,1516010 (gravado 0,1516012), q_TT 0,2738869 (0,2738872), max ‖Q0‖ em
  48 ≤ |m| < 72 0,1628312 (0,1629740), ‖Q_ZZ‖ 1,1830518 (1,1842320), max ‖Q0‖ em F 1,1830518 (1,1830518). Gravado
  ≥ float em todas, com razões 1,000001 a 1,001.
- `n4_constantes.py` → `n4_constantes.out`: `fundo_L.eps_mu_L(λ0)` = 3,6334382330162156e-8, idêntico ao gravado (e diferente do valor
  em 0,401 ou em μ, o que mostra que foi avaliado neste λ0); `dB_L(banda 16)` idêntico.
- `n4_Qfar.py` → `n4_Qfar.out`: ‖Q0 P_far‖ = max(radial 0,119916755485901; Fourier |m| ≥ 200 0,104992) = 0,119916755485901,
  idêntico a `Z.json`.

### N4.5 tQ ≥ ‖Mh⁻¹ [Q_ZZ; 0]‖ (o termo que domina Z2)

`n4_tQ.py` (laboratorio2; `lab2/n4_tQ.log`, `lab2/n4_tQ.json`): ARPACK dá ‖Mh⁻¹ [Q_ZZ; 0]‖ ≈ 8,102346; o Loewner
t²·Mh Mh* − X X* ≥ 0 (X = [Q_ZZ; 0]) no t gravado (tQ/(1 + 1e-12) − t_V·dQZ = 8,110448) é **verificado**; em
t = 0,999 × 8,102346 = 8,094243 é **rejeitado**. tQ = 8,11045 gravado é válido e justo (razão 1,001).

### N4.6 O acoplamento a (Z′←T): estimativa independente

O Loewner de a (`etapaA.txt`) passou **na primeira tentativa**, em t = 1,01 × chute = 0,37875, o que não diz
quão perto está do valor real. `n4_aZ.py` (laboratorio2; `lab2/n4_aZ.log`, `lab2/n4_aZ.json`) acumula
S = Y Y* com a mesma montagem (banda completa |d| ≤ 40) e calcula λ_max(Mh⁻¹ S Mh⁻*) por ARPACK (eigsh):

- colunas de T em |m| < 48 (o que o Loewner cobre): a ≈ 0,3463853;
- colunas em |m| < 72 (todas as que alcançam Z pela banda): a ≈ 0,3463853 (as de 48 ≤ |m| < 72 quase não pesam).

Contra a_loewner = 0,37875 e √(a_loewner² + a_extra²) = 0,3787646 usados em Z′←T: a cota é **válida** e ~9% folgada
(a favor da segurança). A estimativa coincide com o `DADOS-Z 0,346385` do `nkZ.txt` do autor (outro caminho de
código, `nk_L.acoplamentos_Z`).

**N4: OK.** Nada do que foi refeito diverge do gravado: as três janelas e Y_Z/ρ_Z/dM bit a bit; t_V, tQ e a
por Loewner (com controle negativo) e por ARPACK; Q.json, ε_μ, dB_L e Qfar.

Complemento de N4 (`n4_descidas.py` → `n4_descidas.out`): ε_Z = 0,002443870267187596, desc = 0,0008830602004921374 e
q_Zf = 3,2137e-25 refeitos das formas fechadas com os dados de `Z.json`: **idênticos bit a bit**.

## N6. Síntese

### Achados

| id | gravidade | local | evidência | correção sugerida |
|---|---|---|---|---|
| I-1 | imprecisão (documentação de uma constante de entrada; sem efeito no resultado) | `scripts/nk_L.py:35` (`NORMA_BL = 3.964043`), usada em β, desc, ε_Z e nBL; `scripts/l_real_toro.py:12,36` | `n3_simbolo_L.out`: com a folga refeita em Arb, ‖B_L(RefA)‖ ≤ 3,9640420665 ≤ 3,964043 (folga 9,3e-7); com a correção genérica 40 ε_ω embutida seria 3,9640432586 **>** 3,964043. A constante só é válida como cota do fundo **dos dados**, e o certificado está certo porque soma o erro de RT à parte (β = (NORMA_BL + rt)·Qfar). Mas o comentário diz só "símbolo com folga em Arb (24/09)", e `l_real_toro.py` chama de "cota certificada do símbolo" o valor 3,9641333 (que inclui a folga em float × 1,001 e 40 ε_ω) e registra que um valor abaixo dele foi apontado como problema (I-3). Um leitor que compare os dois vê uma constante "abaixo da cota certificada". | Documentar em `nk_L.py` a origem (`symbol_slack_arb.py` sobre `symbol-L-faixa*.json`, folga em Arb, **sem** a correção de RT) e que o erro de RT é somado à parte; harmonizar a citação em `l_real_toro.py`; acrescentar ao nível 1 um check `NORMA_BL ≥ t_norma + folga_Arb`. |
| C-1 | cosmético | `paper/sec/roots.tex:59` (coluna λ0) | `n5_artigo.out`: o centro do certificado é o double 0,73317965966069242877978…; a tabela imprime o `repr` 0,7331796596606924 (diferença 2,9e-17), enquanto o texto diz que tudo foi arredondado a partir dos racionais exatos. Folgas do intervalo: 7e-9. | Nota "λ0 é o double mais próximo do decimal mostrado". |
| C-2 | cosmético (rastreabilidade) | `build/nk733/T3.json`, `T3.joao.json`, `etapaB_joao.txt`, `C3.json` | 7 das 26 janelas vieram de outra máquina; nenhum arquivo registra máquina, ambiente (`CHOPTUIK_EPS_FUNDO_L`) ou o hash do `h0w.npy`/`Z.json` usados. Resolvido por N4.1 (a janela −191..−176 dessa rodada sai idêntica bit a bit) e pelos Y idênticos. | Gravar em T3.json/C3.json o sha256 de `h0w.npy`, `Z.json` e `Q.json` e o valor de `EPS_FUNDO_L`, e o nome da máquina nos checkpoints. |
| C-3 | cosmético (verificador) | `scripts/verifica_T2.py::check_nk` | compara só θ, r e erro_lambda (r já depende de Z2 e Y, então a cobertura é quase completa). Aqui comparei o `C3.json` inteiro: idêntico bit a bit (N2). | Comparar o `C3.json` inteiro. |

Observações sem gravidade (folgas não aproveitadas, todas a favor da segurança): o Loewner de a aceita
0,37875 contra a ≈ 0,34639 (N4.6); Y_J é, na prática, ‖V_J‖·‖x_J‖ (na janela 32..47, 3,84e-6 contra
‖V_J x_J‖ ≈ 2,30e-6, N4.2); as cotas de Rump do autor para ‖N‖₂ e ‖Nf‖₂ são 6e-4 e 1e-3 mais frouxas que
Collatz–Wielandt (N2); rt e o termo de μ entram em dobro nas janelas de SZ. Só apertar Y_J já reduziria r em algo como 40%
(estimativa grosseira, não calculada como cota).

Contagem: **0 erros, 0 lacunas, 1 imprecisão, 3 cosméticos.**

### O que não foi checado (escopo)

A maquinaria partilhada com S3b, S4 e a faixa B, que a tarefa manda não reauditar: o critério de Rump
(`verified_norms._margem_rump`) e o Loewner em si (só testei que ele aceita os valores gravados e rejeita valores
abaixo do verdadeiro); a montagem de B a partir de `RefA.dat` (`signed_operator_L.bloco_B`); as formas fechadas
(`closed_form`) em si; a cota 4,7e-8 do erro do fundo de RT (de S4); a cota do símbolo em Arb (só refiz a folga
com o script do autor). Das 26 janelas, refiz 3; a matriz Nf (far←T) e as outras 23 janelas entram como dados.

### Veredito

**O certificado sustenta a afirmação.** Com os dados gravados e refeitos: existe um único zero (y*, λ*) do sistema
bordejado na bola B(x0, r) da norma ‖·‖_v, com r = 2,1577e-5 ≤ 2,16·10⁻⁵ e |λ* − λ0| ≤ v0·r = 1,01623e-5 ≤
1,0163·10⁻⁵; θ + 2 Z2 r = 0,806 < 1, logo DF(x*) é invertível e λ* é raiz algebricamente simples de L; e
s* ∈ [0,73316949; 0,73318983], arredondado para o lado seguro. Uma montagem independente (racionais, outra cota de
norma, perturbação global contada uma vez) dá θ = 0,805501 e |λ* − λ0| ≤ 1,01535e-5, dentro do afirmado. Nenhum
dado de entrada refeito divergiu do gravado. A validade é condicional à maquinaria partilhada listada acima.

## Arquivos desta pasta

`n1_*` (integridade), `n2_montagem.*` e `n2_verifica/` (montagem; a cópia de `scripts/` e os links para `.venv`
e `.cache` usados para rodar `verifica_T2.py` e `nk_L3_C.py` isolados foram apagados depois do uso; ficam as saídas
`n2_verifica/c3_rerun/` e `n2_verifica/proj/build/reproducao/nivel1_nk_0733.json`), `n3_*` (fórmulas e
constantes), `n4_*` e `lab2/` (dados de entrada; scripts rodados no laboratorio2 e logs/JSON copiados de
`~/revisao_nk733/`; `n4_lab2_tQ.sh` foi a primeira versão da fila, encerrada pelo PID antes de rodar e substituída
por `n4_lab2_tQ2.sh`), `teste_local/` (teste de fumaça de `n4_janelas_float.py` na janela 192..199, no notebook, e
as listas de sha256 da cópia), `n5_artigo.*` (artigo), `n6_manifesto_snapshot_fim.txt`. No laboratorio2 escrevi
só em `~/revisao_nk733/` (2,6 MB no fim; os arquivos temporários P_Z.npy/P_tQ.npy foram apagados pelos scripts;
nenhum processo ficou rodando).

## Conferência do manifesto (fim)

- `build/nk733/MANIFEST.sha256`: 18/18 SUCESSO.
- `.snapshots/2026-10-09-revisao-nk733/MANIFEST.sha256`: **271/274 SUCESSO; 3 FALHOU**: `paper/main.tex`,
  `paper/sec/computation.tex`, `paper/sec/declarations.tex` (`n6_manifesto_snapshot_fim.txt`).
- Causa: o `HEAD` andou durante a revisão, por commits do autor, não desta revisão:
  `436434b` (20:23:59, README bilíngue), `199b251` (20:24:18, "versão de preprint"), `a6c372e` (20:32:17,
  agradecimentos). `git diff --stat d85f649 HEAD`: README.md, paper/main.tex (2 linhas), paper/sec/computation.tex
  (1 linha removida), paper/sec/declarations.tex (6 linhas). `git diff` da árvore de trabalho está vazio.
- Nenhum arquivo usado nesta revisão mudou: `scripts/*.py`, `build/nk733/**`, `paper/sec/roots.tex`,
  `paper/sec/main_theorem.tex` e `paper/sec/counting.tex` continuam com os hashes do snapshot. Os resultados
  acima valem para o snapshot d85f649 e, nesses arquivos, também para a6c372e.
