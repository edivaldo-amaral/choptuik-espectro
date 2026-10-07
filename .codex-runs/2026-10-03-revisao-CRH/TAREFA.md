# Revisão independente leve: Lema C, eixo imaginário e H-reg′ (03/10/2026)

Você é o revisor independente de três complementos curtos a um resultado assistido por computador. O
autor é outro modelo. Seu trabalho é **tentar quebrá-los**. Um "OK" sem artefato (derivação escrita
passo a passo, ou script com saída) não conta. Escreva em português do Brasil.

## Os três objetos

1. **Lema C** (`docs/GAUGE_GLOBAL.md` §6). O modo de gauge e^{μτ}g* em s = μ é a derivada de uma família
   de soluções exatas ψ_ε*(ḡ, φ*), com ψ_ε a translação conjunta das coordenadas nulas. Essas soluções
   são isométricas ao fundo e estão na forma do ansatz de RT.

   Consequência alegada em `docs/T2_ENUNCIADO.md`, item 2 do teorema: com o gauge linearizado padrão
   (L_X, X campo C¹ em M̂) a contagem física é 1; com gauges que fixam o cone (f(0) = 0) é 2.
2. **Eixo imaginário** (`docs/C1_REAVALIACAO.md` §13 e item 3 de `T2_ENUNCIADO.md`). A alegação: a única
   raiz de L em Re s = 0 é s = 0, algebricamente simples. Ela usaria:
   - as coberturas das arestas Re s = 0 dos Rouchés das faixas A e B
     (`build/rouche_L_rig/cobertura_v2.txt`, `build/faixaB/cobertura.txt`);
   - a i-periodicidade de Σ_A;
   - a fórmula de multiplicidade da deflação de posto finito, com o quadrático q de §9.

   **Para este item você PODE ler `docs/C1_REAVALIACAO.md` §6, §9, §11 e §13**: é a fonte das
   alegações.
3. **H-reg′** (`docs/HREC_IDA.md`, R5 e §3). A recuperação de raio vale com peso de Fourier κ1' = 1, e por
   isso a classe física pode ser "C^∞ em τ, real-analítica em ξ".

## Regras

- **Pode ler:**
  - `docs/GAUGE_GLOBAL.md`, `docs/GAUGE_RT_ACTION.md`, `docs/GAUGE_DOMAIN_CLASSIFICATION.md`,
    `docs/HREC_IDA.md`, `docs/RADIUS_RECOVERY.md`, `docs/T2_ENUNCIADO.md`, `docs/T2_HIPOTESES.md`,
    `docs/SECTOR_B.md` e `docs/SPECTRAL_PROBLEM.md`;
  - `docs/C1_REAVALIACAO.md`, só as §6, §9, §11 e §13;
  - os scripts em `scripts/*.py` (por exemplo `rouche_L_cobertura.py`, `rouche_L_contagem.py`,
    `sharp_exterior.py` e `radius_recovery_geral.py`) e `build/**`;
  - o artigo de RT em `.cache/rt-1203.3766v1/arXiv-1203.3766v1.tar.gz`. Extraia para o scratchpad da sua
    sessão, não para esta pasta.
- **NÃO pode ler:**
  - as outras seções de `C1_REAVALIACAO.md`;
  - `docs/PROOF_OBLIGATIONS.md`, `docs/RETOMADA.md`, `docs/PRIORIDADES.md`;
  - `docs/AUDITORIA_*.md`;
  - outras pastas de `.codex-runs/`;
  - a memória do autor.
- **Não altere nenhum arquivo existente.** Escreva só em `.codex-runs/2026-10-03-revisao-CRH/`.
- **Use só o notebook,** com `.venv/bin/python` e `OPENBLAS_NUM_THREADS=2`, sem ssh. O `.venv` não tem
  sympy; se precisar, use `pip install --target <scratchpad>`.
- **Snapshot:** `.snapshots/2026-10-03-revisao-CRH/MANIFEST.sha256`, com o commit em `COMMIT.txt`. Confira
  no início e no fim.
- **Orçamento:** ~1 h.

## Artefatos

- **C1.** Derive que ψ_ε preserva du_± e o eixo, e que a métrica transladada está na forma de RT. Escreva
  ζ_ε e W_ε, e confira o fator (u_+ + u_-)².
  - Confira que d/dε ω_ε = δω_{X0} = e^{μτ}(Z + 1)ω*. Por exemplo, em sympy com funções genéricas, ou
    rodando `scripts/test_gauge_action.py`.
  - Diga se "contagem 1 com o gauge padrão, 2 com cone fixo" segue, e se X0 é de fato um campo C¹ em M̂.
- **C2.**
  - Confira nos arquivos de cobertura que as arestas Re s = 0 de A e de B estão cobertas e que há um
    disco centrado em 0, com o raio dele.
  - Confira a fórmula mult_L̃(s) = mult_L(s) + ord_s(q/(s(s − μ))). Sob quais hipóteses vale para
    perturbações de posto finito de famílias de Fredholm? Os vetores de deflação são os exatos?
  - Confira que q não se anula em Re s >= 0.
  - Confira que [−1/4, 3/4] em Im s, com a i-periodicidade, cobre o eixo inteiro do problema de período
    completo.
  - Procure furos: os cantos, o ponto s = 0 sobre o contorno e a convenção dos setores.
- **C3.** Confira que a cota de cauda free_tail não depende do peso de Fourier, que a norma de B com peso
  1 é <= a com peso κ1 > 1, e que um perfil C^∞(ℝ/4πℤ; A(E)) está em Y_{(1, κ2')} junto com J_μh.
  Confira também que as peças R1–R4 e R6–R7 do Lema R preservam C^∞ em τ.

## Relatório

`.codex-runs/2026-10-03-revisao-CRH/RELATORIO.md`:
- uma tabela C1–C3 com veredicto (OK / OK com ressalva / LACUNA / ERRO) e o artefato;
- as lacunas, com gravidade;
- o manifesto no início e no fim.
