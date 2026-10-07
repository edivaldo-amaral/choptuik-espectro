# Relatório da revisão independente da faixa B (C4) — 02/10/2026

Revisor: agente independente (outro modelo). O objetivo foi tentar quebrar o certificado. Cada veredicto abaixo
tem script e saída nesta pasta.

- **Snapshot local:** `MANIFEST.sha256` conferido antes de qualquer conta, 413/413 OK (`r10_sha_antes_check.txt`), e
  depois de todas as contas, 413/413 OK (`r10_sha_depois_check.txt`).
- **Manifestos remotos, conferidos no fim:**
  - laboratorio2: 60/60 OK, `r10_sha_remoto_lab2.txt`;
  - João: 21/21 OK, `r10_sha_remoto_joao.txt`.
- **Máquinas:**
  - laboratorio2, para R2, R3 (a0), R7 (autovalores e diag de T);
  - PC 3, para R6 (conta densa 11200);
  - notebook, para o resto.
  - O PC do laboratório não foi usado.
  - Escrita remota só em `~/revisao_B/`, com uma exceção: um `r2.pid` caiu por engano em `~/` no laboratorio2 e foi
    movido para `~/revisao_B/` em seguida.
- **Ambiente:** `CHOPTUIK_EPS_FUNDO_L=4.7e-08`, `.venv/bin/python`.

## Tabela de veredictos

| R | Veredicto | Artefatos | O que foi feito |
|---|---|---|---|
| R1 (B1) | **OK com ressalva** | `r1_formulacao.py`, `r1_saida.txt`, `r1_resultado.json` | Ver nota R1 abaixo. |
| R2 (B2) | **OK** | `r2_resolvente.py` (laboratorio2), `r2_saida.txt`, `r2_resultado.json` | Ver nota R2 abaixo. |
| R3 (B3) | **OK com ressalva** | `r3_rig.py`, `r3_saida.txt`, `r3_resultado.json`; a0 em `r2_saida.txt` | Ver nota R3 abaixo. |
| R4 (B4) | **OK** | `r4_caudas.py`, `r4_saida.txt`, `r4_resultado.json`; `r4_nV.py`, `r4_nV_saida.txt`, `r4_nV_resultado.json` | Ver nota R4 abaixo. |
| R5 (B5) | **OK** | `r5_discos.py`, `r5_saida.txt`, `r5_resultado.json`, `r5_saida_normas_1_0.txt` | Ver nota R5 abaixo. |
| R6 (B6) | **OK** | `r6_janelas.py` (PC 3), `r6_saida.txt`, `r6_resultado.json` | Ver nota R6 abaixo. |
| R7 (B7) | **OK** | `r7_cobertura.py`, `r7_cobertura_saida.txt`, `r7_cobertura_resultado.json`; `r7_diagT.py`, `r7_diagT_saida.txt`; autovalores em `r2_saida.txt` | Ver nota R7 abaixo. |
| R8 (B8) | **OK com ressalva** | `r8_testemunho.py`, `r8_saida.txt`, `r8_resultado.json` | Ver nota R8 abaixo. |
| R9 (B9) | **OK com ressalva** (condicional) | este relatório, §R9 | Ver §R9. |
| R10 | **OK** | `r10_cobertura_contagem_saida.txt`, `r10_disco_0_25*.json/txt`, `r10_sha_antes_check.txt`, `r10_sha_depois_check.txt` | Ver nota R10 abaixo. |

### Notas da tabela

**R1 (B1).**
- **i-periodicidade no código:**
  - `max|L(s+i)[m,m'] − L(s)[m+2,m'+2]| = 1,1e-16`;
  - `bloco_B(m+2←m'+2) = bloco_B(m←m')` exatamente (as paridades por componente só dependem de m mod 2).
- **F3:** a parte hermitiana de W L(s) W⁻¹ tem os mesmos extremos para Im s = 0,1, 0,6 e −0,9, até 10 casas. Só
  Re s entra.
- **Fórmulas da montagem testadas para |Im| até 0,75 + r e além:**
  - `‖R⁻¹‖ ≤ c/max(μ,|Im σ|)`: razão máxima 0,35 em 90 casos;
  - `cauda_Rinv_disco`: razão numérico/cota máxima 0,90 em 10 discos;
  - `eps_mu_L` é monótona em |Re| e |Im|.
- **Brauer:** os polos de D(s) estão em 0 e 0,16831 e os zeros em −1. Todos são reais, logo ficam fora de Ω_B.
- **Mesma referência:** o mesmo sha de A nos 22 arquivos de discos e uma única partição de janelas nas 13 caudas.
- **Ressalva:** os números de F3 (2,926 etc.) não estão nos arquivos permitidos. Conferi só a independência de Im s.

**R2 (B2).**
- **Cota:** rederivei a cota de Loewner pré-condicionada (ver §Detalhes). Está correta.
- **‖(A+p)⁻¹‖ verdadeiro, com LU densa, Lanczos e cota inferior independente do solver:**
  - 0,502449i: verdadeiro 201,1476, t_p = 201,8549;
  - 0,47849i: verdadeiro 176,5567, t_p = 177,1690;
  - 0,464083i: verdadeiro 150,0446, t_p = 150,5558;
  - t_p/verdadeiro fica entre 1,0034 e 1,0035, ou seja, t_p está sempre acima.
- **dK:** não falta nenhum termo relevante (ver §Detalhes).

**R3 (B3).**
- **Arquivos:** os 12 rig da faixa B (8 em `z/`, 4 em `teste/`).
- **O que bate:**
  - sha de A = 70791ce0…;
  - sha de F = o do certF do centro;
  - versão 2;
  - pontos iguais aos do tp_*.json e aos dos discos;
  - t_p ≥ t_p da fase 1 do log, com diferença ≤ 4,9e-7, que é o arredondamento da impressão.
- **a0 refeito com LU densa e código próprio em 0,502449i:** 0,4252377. O rig dá 0,4253153, acima, como deve ser.
- **Ressalva:** a fase 2 sobrescreve o tp_*.json da fase 1 (lacuna L2).

**R4 (B4).**
- **As 8 janelas em comum** (80..95 até 192..199) entre T3.parcial (João) e T3.lab2 são iguais bit a bit.
- **O T3.json juntado:**
  - tem 26 janelas;
  - cada entrada vem de um checkpoint: 14 do João, 12 do laboratorio2;
  - é igual a `build/faixaB/cauda/T3_0_5.json`.
- **A cópia antiga de T3.lab2.json no João** é um subconjunto idêntico da do laboratorio2.
- **‖V_J‖ por SVD densa, com pesos, Q0 e montagem próprios:**
  - 192..199: 1,05272, contra 1,05409 do T3;
  - −199..−192: 1,05475, contra 1,05607.
- **md5 dos scripts:**
  - iguais nas três máquinas;
  - exceção: `nk_L2_Z.py` e `rouche_L_Z.py` do João, que são versões anteriores e só não têm `--so-tQ` e `certF`.

**R5 (B5).**
- **Os 3 discos, montados em racionais com código próprio:**

  | Disco | θ do revisor | θ do autor |
  |---|---|---|
  | 0,502449i, o pior: θ do autor 0,99907, r = 2,85e-3 | 0,998849 | 0,999074 |
  | 1 + 0,75i, cauda real, d = 0,873 | 0,980804 | 0,980574 |
  | 0,75i, canto | 0,998604 | 0,998836 |

- **Comparação com a matriz 3×3 do autor:**
  - Z'←Z', Z'←T, Z'←far, T←Z', T←far, far←Z' e far←far são iguais;
  - T←T e far←T diferem ≤ 1e-3, por causa de outra cota de norma (Collatz–Wielandt).
- **‖N_s‖₂ por SVD** = 0,61272. As duas cotas, a do autor e a minha, ficam acima.

**R6 (B6).**
- **Partição:** a de REGIOES_B cobre [0,3]×[1/4,3/4] (conferido em racionais).
- **Cota refeita com normas por modo próprias:**
  - 1,1+0,75i, região do centro 1, d = 0,857, janela −47..−32: 0,1517, igual à do autor;
  - 0,3+0,5625i, janela central −15..0: 0,1486, igual à do autor.
- **‖I − V_J(c)Ĥ_J(s)‖ direto** (LU 11200 e Lanczos), em s = 1,1, 0,8 e 1,2 + 0,75i: 0,1105, 0,1202 e 0,1112. Todos
  ficam abaixo de 0,1517.

**R7 (B7).**
- **Verificador próprio:** sem raiz quadrada, uma cadeia de pontos racionais em que cada par consecutivo está no
  mesmo disco, com os cantos incluídos. Resultado: CONTORNO COBERTO, com 93 discos nas cadeias; max t = 516,5 e
  t·ε_B ≤ 0,0308.
- **diag(T) do laboratorio2:**
  - é igual, bit a bit, ao `diagT_lab2.npy`;
  - os 14 valores de `contagem.json` aparecem bit a bit nela, logo ε_B e a diagonal vêm do mesmo Schur.
- **Autovalores, por shift-invert em dois deslocamentos:** exatamente 1 autovalor de A em −Ω_B, −0,0414194105−0,5i.

**R8 (B8).**
- **λ0 e λ*:**
  - λ0 difere do autovalor de A em 2,7e-11;
  - λ* ∈ Ω_B: Re ≥ 0,041383 e Im ∈ [0,499963, 0,500037].
- **Z2 com tQ:** rederivado, está correto (§Detalhes).
- **c_ref:**
  - 0,355192 com montagem e pesos próprios, contra 0,355189 do D.json;
  - estável na caixa F: 0,355190 (8×32) e 0,355192 (16×64).
- **Conferência independente com as matrizes de constraints exportadas pelo próprio RT** (10×30 e 12×36),
  ‖C(s)h‖/‖h‖:

  | Raiz | ‖C(s)h‖/‖h‖ |
  |---|---|
  | raiz da faixa B | 0,504 / 0,503 |
  | 0,401 | 0,324 |
  | 0,733, 0 e 0,168 | de 1e-4 a 8e-4 |

- **Ressalvas:**
  - c_Z e as perdas do D.json não foram refeitos;
  - c_ref usou os blocos não pesados `bloco_C` do autor.

**R10.**
- **cobertura e contagem reproduzidas exatamente:**
  - 174 discos, CONTORNO COBERTO;
  - max t 559,151, t·ε_B 3,337e-2, 1 zero em Ω_B.
- **`rouche_L_disco.py` para o centro 0,25i, numa cópia do `--dir`:**
  - com `--alvo` 0,99 (o padrão), os raios saem menores e θ ≤ 0,9893;
  - com `--alvo` 0,998, o 1º ponto é reproduzido exatamente: r = 1,811e-2, θ = 0,995933;
  - nos outros 5 pontos, os raios saem 0,2–0,6 % menores (o autor usou alvo maior; L3);
  - todos têm θ < 1.
- **sha256** antes e depois:
  - local 413/413;
  - laboratorio2 60/60;
  - João 21/21.

## Detalhes por alegação

### R1 — formulação (B1)

1. **Redução e periodicidade.**
   - No código, J_livre usa σ = s + i m/2.
   - A estrutura de componentes depende de m mod 2: m par carrega v1, v2, v3; m ímpar carrega v4.
   - O acoplamento só depende de d = m_o − m_i e da paridade.
   - Logo S⁻¹L(s)S = L(s+i) para o deslocamento S de m por 2. Os pesos κ1^{2m} mudam só por uma constante κ1⁴, e
     S é isomorfismo do espaço com peso. Verificado numericamente (1,1e-16).
   - A conjugação do setor B (SECTOR_B §2–§5) é uma conta de uma linha. Ela está correta como álgebra.
   - A redução ao setor A consome H-afim, H-iso, H-dof e H-constr (SECTOR_B §6). Elas não são verificáveis por
     conta; ver R9.
2. **F3.**
   - Re⟨W(J_μ + B + s)W⁻¹x, x⟩ = Re s‖x‖² + Re⟨(J_μ⁰ + B)x, x⟩, porque o termo i(Im s + m/2) é diagonal e
     imaginário puro no produto com peso.
   - Logo a cota F3 vale para todo Im s. Teste (2) de `r1_saida.txt`.
3. **Brauer.**
   - D(s) = Π(s − s_k + δ_k)/(s − s_k), com s_k ∈ {0, 0,16831} e δ_k ∈ {1, 1,16831}.
   - Polos e zeros são reais, logo não há correção em Ω_B: N(L̃, Ω_B) = N(L, Ω_B).
   - Os translados s_k + i (autovalores de L em Im s = 1) não são deflacionados. Estão fora de Ω_B, e o shift-invert
     do R7 os mostra em Im = 1.
4. **A mesma referência.**
   - A(s) = diag(Ĥ_ZZ(s), Ĥ_J(s), I) depende de:
     - A.npy, com o mesmo sha nos 22 arquivos de discos das duas faixas;
     - da partição das janelas: um único hash de `modos` nas 13 caudas;
     - de Z = 32×128 e F = 200×400.
   - Os discos da faixa A em Im = 1/4 valem para todo s do disco com Re s ≥ 0. Servem, portanto, no bordo de Ω_B.
   - A norma pode variar com s (pesos de Perron, escala D_Z): para Rouché basta a invertibilidade de A + t(H − A)
     ponto a ponto no contorno.
5. **Dependência da faixa A nas fórmulas.** Nenhuma encontrada:
   - `eps_mu_L` usa h = max(μ, |m|/2 − |Im s|) e o pior canto do disco;
   - `Qfar_disco` usa h = Mc/2 − |Im p| − r ≈ 99,1;
   - `cauda_Rinv_disco` só usa Re σ ≥ 0, porque |d_k|² − |d_k − v_k|² = 2μn(2μ + 2Re σ);
   - nQx = c/μ é uniforme.
   - Os testes numéricos (3) e (4) não falsificaram nenhuma delas.

### R2 — Loewner pré-condicionada (B2)

**Derivação.**
- J_Z(p) = J_Z(0) + pI vale exatamente: a similaridade por pesos preserva a identidade.
- Logo Q_e(A + p) = I + Q_e B̂_e = K_e e (A + p)⁻¹ = K_e⁻¹Q_e.
- **Loewner.** t²K_cK_c* ⪰ Q_cQ_c* dá ‖Q_c*y‖ ≤ t‖K_c*y‖. Como Q_c é invertível, K_c é invertível e
  ‖K_c⁻¹Q_c‖ ≤ t.
- **Inversa de Q_c.** Q_c⁻¹ = (I + J_Z(p)(Q_c − Q_e))⁻¹J_Z(p), logo ‖K_c⁻¹‖ ≤ t·nJ/(1 − nJ·dq) = kc.
- **Passagem para K_e.** Com K_e = K_c + ΔK:
  - ‖K_e⁻¹Q_c‖ ≤ t/(1 − kc·dK);
  - ‖K_e⁻¹(Q_e − Q_c)‖ ≤ kc·dq/(1 − kc·dK).
- A fórmula do código confere.

**Termos de dK conferidos:**
- fl(A − J0_c): u|X| na linha de blocos;
- J_Z(0) calculada contra a exata: 4u|J0|, que basta, pois são duas operações de arredondamento por entrada e μ(n+1),
  2μn são exatos;
- o produto Q_c·X por blocos: γ‖|Q||X|‖_F;
- a soma com I: u‖Y‖;
- ‖Q_c‖·erro(X);
- (Q_c − Q_e)B̂_e: dq·‖B̂‖, com ‖B̂‖ ≥ nBhat + 4u‖J0‖_F, que é até contado em dobro (conservador);
- erro de P = K Kᴴ: γ_n ‖|K|‖₁‖|K|‖_∞;
- a montagem de G e o critério de Rump em `pd_hermitiana`: revisados.

**Ponto fino.**
- O `Q_cert` cota ‖Q_c − fl(Jw)⁻¹‖, e não ‖Q_c − Jw_exata⁻¹‖.
- Além disso, os pesos do `Q_cert` (√ω) e os de `JZ_blocos` (κ1^m√ω) diferem por arredondamento.
- A diferença é ≤ 4u‖|Jw||Q|‖_F. Ela é absorvida pela folga de `err_mm`, que usa γ_{n+8}‖Jw‖_F‖Q‖_F em vez de
  γ_n‖|Jw||Q|‖_F; a folga é ≥ 36u‖Jw‖_F‖Q‖_F.
- Não é lacuna, mas vale um comentário no código.

**Ordem de grandeza.**
- kc·dK ≈ 5e-4 em t ≈ 200: kc ≈ 3,3e4 e dK ≈ 1,6e-8. A correção de t para t_p é de +0,05 %, e o assert < 1e-3 passa
  com folga de 2×.
- Os t_p verdadeiros confirmam que a cota final fica 0,35 % acima do valor verdadeiro.

### R3 — nível Z (B3)

- A fase 1 da Loewner pré-condicionada foi usada em todos os pontos com estimativa > 40 (log `rig_B.txt` do
  laboratorio2: "Loewner pre-condicionada").
- Os pontos com t_p < 40 dos rig `Bteste` (1,0, 1,5, 2,0 e 3,0 em Im = 3/4) vêm da forma direta (log
  `build/caudaB/teste_reuso.txt` do laboratorio2).
- O rig grava t_p corrigido pelo deslocamento arredondado. Ele trata os t_p reaproveitados como "não
  pré-condicionados", o que só aumenta t_p (conservador).
- a0 por LU: 0,4252377. O rig dá 0,4253153, com erros do Schur incluídos. Coerente.

### R4 — caudas (B4)

- Os checkpoints, a junção e os md5 conferem.
- Nas janelas recalculadas, o nV certificado fica 0,12–0,13 % acima do valor por SVD. É compatível com o passo 1,002
  da `cota_loewner`.

### R5 — discos (B5)

- **A escala do nível Z** (D_Z = J_Z(c)Q_ZZ(s)) foi rederivada:
  - D(J_Z(s)R(s)ΔB Q_ZZ(s))D⁻¹ = J_Z(c)R(s)ΔB Q_ZZ(c);
  - isso justifica ‖Q_ZZ(c)‖ e eps_mu(c) em Z'←Z';
  - e o TZ_J(c) exato em T←Z'.
- **Disco 1+0,75i** (cauda real de 1,0 na aresta Im = 3/4):
  - d = 0,873;
  - d·q_TT, d·nQ_JJ e r·t_p ficam abaixo de 0,9;
  - θ = 0,9808.

### R6 — janelas (B6)

- O princípio do máximo vale em cada retângulo: Q0(s) é analítica em Re s > −μ, e a deflação é constante em s.
- A cota no pior ponto é reproduzida.
- O valor direto de ‖F(s)‖ na aresta Im = 3/4 da região de centro real 1,0, com d até 0,78, é ≈ 0,11–0,12, contra
  a cota 0,1517.

### R7 — cobertura e contagem (B7)

- **Cobertura.** O verificador próprio e o do autor concordam.
- **Contagem:**
  - a contagem do autor usa todos os discos das duas faixas, inclusive Im = −1/4 (max t = 559, conservador);
  - só nos discos da cadeia de cobertura, max t = 516,5;
  - só nos da faixa B mais os de Im = 1/4: 538,7.
  - Em todos os casos t·ε_B < 0,034.
- **Os −T_ii mais próximos do bordo:**
  - −0,2453+0,5i, a 0,245 de Re = 0;
  - os translados em Im = 1, a 0,25 de Im = 3/4.

### R8 — NK e testemunho (B8)

**Segunda ordem.**
- DF(x) − DF(x0) = [(λ−λ0)Q0 dy + dλ Q0(y−y0); 0]. A linha l é linear.
- Com Q0 diagonal em m e triangular superior em n:
  - a parte de Z' é Mh⁻¹[P_Z Q0 v; 0], com
    ‖·‖ ≤ tQ‖v_Z‖ + t_V(q_ZT‖v_T‖ + q_Zf‖v_f‖);
  - a parte de T é ≤ maxwV(q_TT‖v_T‖ + Qf‖v_f‖);
  - a parte do far é ≤ Qf‖v_f‖.
- ‖v_nível‖ ≤ 2v0‖x−x0‖‖dx‖ × peso do nível. Isso dá exatamente o Z2 do `nk_L3_C`.

**tQ.**
- Loewner t²Mh Mhᴴ ⪰ [Q_ZZ;0][Q_ZZ;0]ᴴ, com a mesma Mh e o mesmo h0w das etapas.
- Mais t_V·dQZ.
- Correto.

**h0.**
- h0 é um dado. A qualidade da iteração inversa pesada só entra via Y (6,0e-6, certificado), então não precisa ser
  confiável em si.

**Identificação.**
- N(L, Ω_B) = 1 com multiplicidade.
- O NK dá um zero λ* ∈ Ω_B algebricamente simples: o Jacobiano bordejado é invertível, e o pencil é linear,
  L(λ) = L(0) + λ.
- Logo λ* é a única raiz em Ω_B. O espaço de raízes é E = span h*.
- **Constraints.** D|E = C0 − C1A restrito a E vale C(λ*)h* ≠ 0, logo dim E_constraints = 0
  (CONSTRAINT_ROOT_QUOTIENT).

## R9 — lógica física

**Contagem:**
- N(L, Ω_A) = 3, com raízes 0,401, 0,733 e μ;
- N(L, Ω_B) = 1, com raiz λ* ≈ 0,0414+0,5i;
- não há espectro em Re s ≥ 2,926 (F3, para todo Im s);
- não há zeros nas retas Im = −1/4, 1/4 e 3/4, nem em Re = 0, dentro das faixas (discos com θ < 1).

**Multiplicidades físicas,** por dim E − rank(D|E), menos o gauge:

| Raiz | Multiplicidade física | Por quê |
|---|---|---|
| 0,401 | 0 | testemunho, faixa A |
| 0,733 | 1 | K injetivo: KC = ML dá C(s)h = 0 |
| μ | 0 | gauge subtraído, S4 |
| λ* | 0 | testemunho e simplicidade |

**Total:** 1 por classe de Floquet. Q̃ = (−1/4, 3/4] é um domínio fundamental de iℤ para Σ_A, e Σ_B = Σ_A − i/2
(Lema 5.1).

**A conta está bem fundamentada como enunciado condicional.** O que é hipótese, e não está provado por este
certificado:

1. **Sector B:**
   - H-afim, H-iso, H-dof e H-constr (SECTOR_B §6).
   - Em particular, H-constr: o testemunho calculado com C_A em λ* precisa ser o de C_B para o modo do setor B. É
     uma conta de uma linha, mas depende da equivariância T do sistema ♯ de RT.
   - H-conv está confirmada pelo artigo.
2. **Identificação física (CONSTRAINT_ROOT_QUOTIENT, "Hipóteses ainda abertas"):**
   - toda perturbação admissível cabe no ansatz RT;
   - o gauge global permitido é o residual classificado;
   - cone móvel, para que a linha de gauge em μ seja subtraída.
3. **S4:** μ simples e igual ao gauge, com D g = 0.
4. **O espaço:** o espectro de L no espaço do toro com pesos (65/64, 5/4) é o espectro físico (recuperação de
   raio/analiticidade das autofunções).
5. **Dados:**
   - erro do fundo verdadeiro ε = 4,7e-8 (RefB, S4);
   - DELTA_MU para μ;
   - F3 (R_L ≤ 2,926), que não pude reconferir.
6. **Injetividade de K em 0,733** (faixa A).

**Sobre a injetividade de K em 1,765 < Re s < R.**
- A afirmação de que ela deixa de ser necessária é correta: com N_A = 3 e N_B = 1 exatos, não há outras raízes de L
  em 0 < Re s < 3.
- O espectro físico está contido no espectro de L com autovetor que satisfaz as constraints.

## Lacunas

Nenhuma lacuna compromete a conclusão N(L, Ω_B) = 1 nem o testemunho.

| # | Gravidade | Lacuna | Conserto proposto |
|---|---|---|---|
| L1 | baixa (documental) | Ver L1 abaixo. | Corrigir o texto para 78 discos com r > 0, mais 96 da faixa A, ou remover a entrada com r = 0. |
| L2 | baixa (rastreabilidade) | Ver L2 abaixo. | Ver L2 abaixo. |
| L3 | baixa (reprodutibilidade) | Ver L3 abaixo. | Gravar `alvo` e os argumentos no JSON do disco. |
| L4 | baixa (rastreabilidade) | Ver L4 abaixo. | Incluir o sha256 de T3_0_55p0_75.json, e de onde veio, no manifesto. |
| — | nota, não lacuna | Ver nota abaixo. | Anotar no código que o "+8" cobre esse termo, ou somar o termo explicitamente. |

**L1.** "79 discos, todos com θ < 1":
- os arquivos da faixa B têm 79 entradas, mas uma (disco_0_75, p = 0,291688+0,75i) tem r = 0;
- a cobertura usa 78 discos da faixa B mais 96 da faixa A, que dão os 174;
- a cobertura não depende da entrada com r = 0.

**L2.**
- A fase 2 de `rouche_L_rig ponto --tps` regrava o `tp_<c>_<rot>.json`. Os t_p ficam iguais, mas `precond` vira
  `false` em todos.
- O registro da fase 1 só sobrevive nos logs `rig_B.txt` (laboratorio2 e João) e `build/caudaB/teste_reuso.txt`
  (laboratorio2).
- Esses logs não estão no manifesto nem no repositório local.
- **Conserto proposto:**
  - a fase 2 não deve sobrescrever o tp da fase 1 (gravar `tp_..._fase2.json`, ou só ler);
  - arquivar os logs da fase 1 com sha no manifesto.

**L3.**
- O raio de cada disco depende de `--alvo`, que não é gravado no JSON.
- Com o padrão 0,99, `rouche_L_disco.py` dá raios menores que os do arquivo do autor. O autor usou ≈ 0,998.

**L4.**
- A cauda de 0,55+0,75i (`T3_0_55p0_75.json`) foi feita fora das duas máquinas permitidas.
- Não há manifesto remoto dela. Não pôde ser reconferida de ponta a ponta; o número de janelas e a partição foram
  conferidos.

**Nota (não lacuna).**
- `Q_cert` não soma explicitamente o erro de fl(Jw) contra a Jw exata, nem a diferença de pesos √ω contra κ1^m√ω.
- O `certF` soma esse erro.
- O termo, ≤ 4u‖|Jw||Q|‖, é coberto pela folga γ_{n+8} contra γ_n de `err_mm`.

## O que não foi verificado

- Os valores de F3: R_L ≤ 2,926, Re W(−B_L) ≤ 2,9492, parte livre ≤ −0,0232. O certificado não está nos arquivos
  permitidos. Conferi só a independência de Im s.
- Reimplementação de B e de C:
  - os blocos não pesados `signed_operator_L.bloco_B` e `signed_constraints.bloco_C` foram usados como dados;
  - a conferência de C contra o exportador de RT (R8) é qualitativa: outra base, outra norma.
- Das caudas: N_ij, TZ, K_J, nVhi e Nf não foram recalculados. Só ‖V_J‖ de 2 janelas.
- Os certF dos centros novos (0,5i, 0,75i, 0,55+0,75i): ‖Φ‖, ρ, q_TT.
- Do rig: b1, b2, ‖Y3‖ e ‖Q³F_b‖ não foram recalculados. Só a0 em 1 ponto.
- Do NK: as etapas A e B (t_V, ρ_Z, a, Y_Z, ε_Z, T3 do NK), tQ e as perdas c_Z e c_J do D.json não foram
  recalculadas. A montagem C3 e a fórmula de Z2 foram conferidas.
- As hipóteses físicas e de setor da §R9.
