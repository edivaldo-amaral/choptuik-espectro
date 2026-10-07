# Revisão independente de S4 (25/09/2026)

Você é o revisor independente de um certificado assistido por computador. O autor é outro
modelo. Seu trabalho é **tentar quebrar** o certificado, não confirmá-lo. Um "OK" sem artefato
conferível não conta. Escreva em português do Brasil.

## Regras

- **Pode ler:** `scripts/*.py`; `build/s4/**` (JSON, txt, npy); `build/s3b/rig/**` (o
  certificado de S3b, só para comparar); `docs/GAUGE_RT_ACTION.md`,
  `docs/GAUGE_DOMAIN_CLASSIFICATION.md`, `docs/CONSTRAINT_ROOT_QUOTIENT.md`,
  `docs/RADIUS_RECOVERY.md`, `docs/GEOMETRIC_FLOQUET_RECONSTRUCTION.md` (a base da parte
  conceitual; foram escritos por outro modelo em 14–16/09 e auditados em 16/09: trate-os como
  alegações a conferir, não como premissas); `docs/S3_CONDICIONAMENTO_SHARP.md` §10.1–10.3
  (definição do espaço com sinal); o artigo, o código e os dados de RT em
  `.cache/rt-1203.3766v1/` (inclusive `RefAplusB.dat` e `sourcecode/RefA.dat`); o PDF de Rump
  em `.codex-runs/2026-09-24-revisao-agente/`.
- **NÃO pode ler:** `docs/S4_GAUGE.md`, `docs/S3B_ROTA.md`, `docs/RETOMADA.md`,
  `docs/PROOF_OBLIGATIONS.md`, `docs/AUDITORIA_*.md`, nada em `.codex-runs/` além desta pasta e
  do PDF de Rump, nem a memória do autor. São a narrativa e as conclusões do autor (e de uma
  revisão anterior); as alegações estão enunciadas abaixo.
- **Não altere nenhum arquivo existente.** Escreva SÓ em `.codex-runs/2026-09-25-revisao-s4/`
  (scripts `rN_*.py`, saídas `rN_saida.txt`, e o `RELATORIO.md`). Os scripts que gravam em
  `--dir` devem receber uma CÓPIA dentro desta pasta.
- Toda conta de S4 usa `CHOPTUIK_EPS_FUNDO_L=4.7e-08` no ambiente (ver B1). Use
  `.venv/bin/python` (o Python do sistema não tem scipy) e `OPENBLAS_NUM_THREADS=2` no notebook.
- Máquinas: o notebook local tem 7 GB de RAM. Para contas pesadas (matrizes densas de 11.200 a
  14.017) use:
  - o PC do laboratório: `ssh <usuario>@<host>` (chave já instalada), 15 GB, 4 núcleos,
    `~/Choptuik/.venv/bin/python`;
  - o laboratorio2: `<usuario>@<host>`, 14 GB, 4 núcleos, por senha via
    `SSH_ASKPASS=<scratch> SSH_ASKPASS_REQUIRE=force DISPLAY=:0 ssh -o PubkeyAuthentication=no <usuario>@<host> ... < /dev/null`
    (não copie a senha para lugar nenhum).

  Nas duas, trabalhe SÓ numa pasta nova `~/revisao_s4/` (copie o que precisar; não altere
  `~/Choptuik`). Rode coisas longas com `setsid nohup ... &` e acompanhe; ao matar processos,
  use o PID (`pkill -f` casa com a própria linha de comando).
- Orçamento: até ~2 h de trabalho. Priorize o que pode estar errado.

## Contexto

L(s) é o operador selecionado de Reiterer–Trubowitz (RT) linearizado em torno da solução
crítica, na base de Fourier COM SINAL (`scripts/signed_operator_L.py`), com a norma do toro
`||x||² = Σ κ1^{2m} ω2(n) |x_mn|²`, (κ1, κ2) = (65/64, 5/4), ω2(n) = κ2^{2n} + κ2^{-2n}.
L(s) = J_μ + B + s, J_μ a parte livre (diagonal em m, triangular superior em n),
B = μ S Γ1 + 2 S Ξ Γ2(ω, ·) o fundo (linear em ω). Q0 = (J_μ + λ0)⁻¹. RT dão um fundo
aproximado RefA (μ_A = 722873400/2³², diádico) e um mais preciso RefA + RefB
(`RefAplusB.dat`); μ verdadeiro com 80 dígitos (`tile_perron.MU_80`).

Em s = μ há um modo de gauge (translação do instante da singularidade), com gerador
g* = (Z* + 1) ω*, Z* = μ*⁻¹ ∂τ + ξ ∂ξ. S4 quer: μ* é raiz algebricamente simples de L, com
autoespaço span(g*); então, com D g* = 0 (D = operador de constraints), a raiz μ não contribui
para a contagem física.

A maquinaria do NK de L (etapas A, B e C: `nk_L2_Z.py`, `nk_L3_T.py`, `nk_L3_C.py`,
`verified_inverse.py`, `fundo_L.py`) já foi revisada em s = 0,401 (S3b). **Foque no que muda em
μ e no que é novo:** o erro do fundo por RefB, o centro no vetor de gauge, λ0 = μ_A perto da
raiz de fase s = 0 (a 0,168), as normas restritas de Q0 em μ, a identificação (etapa E) e a
lógica da contagem.

## O que é alegado

**B1 (o fundo verdadeiro e RefB).** O fundo verdadeiro é w* = RefA + RefB + Corr, com RefB =
(RefA + RefB) − RefA calculado EXATAMENTE dos inteiros dos dois arquivos, e Corr dentro da bola
de RT (o autor diz ~2⁻²⁷⁷ por componente, nos pesos de RT). Então
‖B_L(w*) − B_L(RefA)‖ <= ‖B_L(RefB)‖ + ‖B_L(Corr)‖ <= 4,6748·10⁻⁸ + 10⁻⁶⁰
(`scripts/fundo_refB.py`, cota de Young, |d| <= 449). Esse valor (arredondado para 4,7·10⁻⁸)
substitui, via `CHOPTUIK_EPS_FUNDO_L`, o erro de RT genérico 40 ε_ω = 1,19·10⁻⁶ usado em S3b,
em todas as etapas. O μ verdadeiro: |μ* − MU_80| <= 10⁻⁸⁰, |μ_A − μ*| <= DELTA_MU = 3,88·10⁻¹¹,
coberto por `fundo_L.eps_mu_L`. Também: ‖RefB‖_L = 2,3506·10⁻¹⁰ e
‖(Z_A + 1) RefB‖_L = 2,4542·10⁻⁸ (`scripts/gauge_refB.py`).

**B2 (o gerador de gauge).** Para todo fundo w, L_w(μ)(Z + 1)w = (Z + 2)F_w (F_w = defeito de
w); logo L_*(μ*) g* = 0 com g* = (Z* + 1)ω* ≠ 0, no domínio forte; D g* = 0.
`scripts/gauge_vector_L.py` monta g_A = (Z_A + 1)ω_RefA na base de L (∂τ → i m/2 e ξ∂ξ nos
coeficientes simétricos). Conferências do autor: |cos(g_A, autovetor da truncagem perto de
μ)| = 0,99999992 (12×48); ‖L(μ) g_A‖/‖g_A‖ = 4,4·10⁻¹⁰ sem pesos; cauda de g_A fora de
Z = 32×128 igual a 2,6·10⁻⁷ relativa.

**B3 (NK de L em λ0 = μ_A).** Mesma formulação de S3b: x = (y, λ), h = Q0 y,
F(x) = (L(λ)h, ⟨h0w, h⟩ − 1), centro h0 = P_Z g_A normalizado (`build/s4/rig/h0w.npy`,
Z = 32×128), aproximação do inverso A = diag(M̂⁻¹, V_J, I) com 26 janelas densas na cauda
F \ Z (F = 200×400). Norma: max(‖x_Z'‖/v0, ‖x_T‖_W/vT, ‖x_f‖/vf).

- Etapa A (`nk_L2_Z.py --h0 build/s4/gauge_h0w_32x128.npy`, saída `build/s4/rig/Z.json`,
  log `build/s4/A.txt`): ‖M̂⁻¹‖ <= 52,14972, a <= 0,40461, ρ_Z <= 1,06·10⁻⁵, Y_Z <= 3,0·10⁻⁶,
  ε_far <= 0,0154.
- Normas restritas de Q0 (`nk_L3_q.py`, saída `build/s4/rig/Q.json`): q_ZT = 0,18412,
  q_TT = 0,32193, extra 0,17892. O arquivo não registra λ0.
- Etapa B (`nk_L3_T.py`): as janelas foram divididas entre duas máquinas, com checkpoints
  complementares `rig/T3.lab.json` (11 janelas) e `rig/T3.lab2.json` (15), juntados na passada
  final (`build/s4/B_final.txt`) → `rig/T3.json`. Logs: `build/s4/B_lab.txt`, `B_lab2.txt`.
- Etapa C (`nk_L3_C.py`, saída `rig/C3.json`, log `build/s4/C.txt`): θ = 0,880172 (Perron em
  racionais), Z2 = 183,7, Y = 6,59·10⁻⁶, r = 6,06·10⁻⁵ com Y + θr + Z2 r² <= r e
  θ + 2 Z2 r < 1. Logo: zero único x* em B(x0, r), DF(x*) invertível, λ* algebricamente simples.

**B4 (identificação; `nk_L3_E.py`, saída `rig/E.json`, log `build/s4/E.txt`).** O par de gauge
exato x_g = (Q0⁻¹ g_n, μ*), g_n = g*/⟨h0w, g*⟩, é zero de F. Pela equação de autovalor,
y_g − y0 = −B_*(g_n − h0) − L_*(λ0)h0 + (λ0 − μ*) g_n − δμ J1 (g_n − h0), então
‖y_g − y0‖ <= (‖B_*‖ ‖g_n − h0‖ + ‖L_*(λ0)h0‖ + |λ0 − μ*| ‖g_n‖)/(1 − |δμ| ‖J1 Q0‖) = 4,9·10⁻⁶,
com ‖g_n − h0‖ <= 2,83·10⁻⁷ e ‖L_A(λ0)h0‖ = 3,709·10⁻⁶ (arredondamento calculado 6,7·10⁻¹¹).
Com v e w do certificado, ‖x_g − x0‖_v <= 1,82·10⁻⁵ <= r. Critério de unicidade usado: como F
é quadrática, ‖DT(x)‖ <= θ + 2 Z2 ‖x − x0‖_v para todo x, e θ + 2 Z2 max(dist, r) = 0,902 < 1.
Conclusão: x_g = x*, então μ* é raiz algebricamente simples de L, com autoespaço span(g*).

**B5 (a contagem).** Com E = autoespaço generalizado de L em μ: dim E = 1, E = span(g*),
D g* = 0, e pela álgebra de `CONSTRAINT_ROOT_QUOTIENT.md` a contribuição física de μ é
dim E − rank(D|E) − 1 = 0. Na realização espectral (que contém g*, isto é, de cone móvel), com
N_A = número de raízes de L em Re s > 0 com multiplicidade, a contagem física é
N_A − 1 (s_K, raiz que viola as constraints, S3b) − 1 (μ, gauge). Para nenhuma outra raiz em
Re s > 0 ser gauge: a classificação de Floquet do gauge residual dá expoentes Re s = μ(1 − n),
n >= 0; só n = 0 é positivo; a fase (n = 1) fica em s = 0. Hipóteses físicas explícitas, não
computacionais:
- (H-rec) toda perturbação física linear se põe no ansatz de RT por um gauge admissível e,
  reciprocamente, todo elemento de ker L(s) que satisfaz as constraints vem de uma perturbação
  física;
- (H-gauge) o gauge residual que preserva o ansatz é o das reparametrizações nulas
  classificadas (provado localmente; a versão global é hipótese);
- (H-cone) a noção física de instabilidade quocienta a mudança do instante da singularidade.

## Artefatos obrigatórios (um por alegação; script + saída + veredicto)

- **R1 (B1).**
  - No artigo de RT, ache o enunciado da bola: centro, raio, norma e raios. Confira que o
    centro é RefA + RefB, e não RefA, e que o raio por componente é o que o autor usa.
  - Com um leitor PRÓPRIO dos dois `.dat` (layout dos índices m, n; `TwoExp`; `num_m` e `num_n`
    diferentes entre os arquivos), recalcule RefB exatamente.
  - Recalcule ‖B_L(RefB)‖ por um caminho independente: pode ser Young com as suas próprias
    tabelas, ou a norma de blocos montados de B(RefB) em alguns d, comparada com a cota.
  - Confira que ‖B_L(Corr)‖ fica de fato abaixo de 10⁻⁶⁰ nos pesos de L, incluindo a conversão
    de norma.
  - Confira DELTA_MU contra MU_80, e onde o erro do fundo entra em cada etapa: procure em
    `scripts/` todo uso de `EPS_FUNDO_L`, de `dB_L(...)['rt']` e de constantes de S3b que
    deveriam ter mudado.
- **R2 (B2).**
  - Confira (derive ou refaça) a identidade L_w(μ)(Z + 1)w = (Z + 2)F_w e a afirmação D g* = 0.
  - Com montagem própria, verifique ‖L(μ_A) g_A‖ e a convenção ∂τ → i m/2 (sinal e escala
    contra a série de Fourier de RT). Uma convenção errada daria outro vetor, que o NK poderia
    ainda assim "certificar" perto de μ.
- **R3 (B3, etapa A em μ).**
  - Numa máquina remota, monte M̂ (bordejado, 14.017) com código seu, ou reusando
    `signed_operator_L` para B e J, e calcule ‖M̂⁻¹‖ em ponto flutuante contra 52,14972.
  - Estime a norma verdadeira de pelo menos um termo de acoplamento (Z'←T, "a") e o resíduo Y_Z.
- **R4 (B3, etapas B e C).**
  - Recalcule `Q.json` com `nk_L3_q.py --dir <cópia>` e confira que bate, o que prova que foi
    feito em μ.
  - Confira em `T3.json`, `Z.json` e `C3.json` que o λ0, o `rt` registrado em `dB_L` e o `eps_mu`
    são os de μ.
  - Com montagem própria, calcule a norma verdadeira de 2 a 3 blocos de E = I − A DF(x0) em μ,
    por exemplo ‖V_J‖ e T←Z' da janela 32..47, e um acoplamento entre janelas vizinhas.
    Verifique que as cotas de `T3.json` estão acima.
  - Confira que a junção dos dois checkpoints não perde nem duplica janelas e que cada entrada
    veio da mesma estimativa (`build/s4/nk3_mu.json`).
- **R5 (B4).**
  - Derive a identidade e as cotas da etapa E.
  - Recalcule por conta própria nz = ‖P_Z g_A‖, a cauda, ‖g* − g_A‖ e ‖L_A(λ0)h0‖.
  - Confira a consistência de normas: pesos, normalização de h0w contra ⟨h0w, ·⟩, e o papel de
    v e w em `dist`.
  - Confira, a partir das definições de θ e Z2 em `nk_L3_C.py` (DF(x) − DF(x0) exato para F
    quadrática), se ‖DT(x)‖ <= θ + 2 Z2 ‖x − x0‖_v vale de fato para todo x.
  - Confira se y_g = Q0⁻¹ g_n está no espaço do NK (domínio).
- **R6 (B5).**
  - Examine a lógica da contagem e os cinco documentos de gauge. Separe o que está provado do
    que é hipótese, e diga se as três hipóteses bastam e se alguma está escondida. Exemplos a
    olhar: a subtração de μ exige D g* = 0 e simplicidade? O quociente dim E − rank(D|E) − 1 está
    certo? A classificação Re s = μ(1 − n) é completa?
  - Aponte o que faltaria para a conclusão "N_físico = N_A − 2".
- **R7 (reprodução).**
  - Rode `nk_L3_C.py` e `nk_L3_E.py` com `--dir` numa CÓPIA de `build/s4/rig` (com a variável
    de ambiente) e confira que reproduzem θ, r e dist.
  - Grave os sha256 dos originais antes e depois.

## Relatório

`.codex-runs/2026-09-25-revisao-s4/RELATORIO.md` deve trazer:
- uma tabela R1–R7 com veredicto (OK / OK com ressalva / LACUNA / ERRO), artefato e o que foi
  feito;
- a lista de lacunas, com gravidade e um conserto proposto;
- o que você não conseguiu verificar.
