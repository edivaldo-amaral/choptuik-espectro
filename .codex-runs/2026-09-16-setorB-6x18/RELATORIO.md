# Relatório: certificados de winding do setor B em 6×18

Escrito pelo orquestrador (Claude) em 17/09/2026. O Codex executou quase toda
a tarefa, mas sua cota de uso acabou antes do fim.

## Cronologia

1. **16/09, 19:29–20:07, Codex.** Calculou dois determinantes Bareiss exatos
   com uma extensão GMP em C (`bareiss_gmp.c`): `z=1/32` e `z=1/10`, cada um
   vinculado ao sha256 da matriz. O cálculo modular do polinômio chegou ao
   primo 500/877 e foi morto quando o sistema ficou sem memória.
2. **17/09, 01:26–02:42, Codex, sessão retomada.** Recalculou o polinômio
   com checkpoints (`checkpoints-execucao2/`) e conferiu a igualdade racional
   EXATA contra os dois Bareiss (`igualdades-exatas.json`). Iniciou o
   contorno grande e parou ao esgotar a cota da OpenAI (`codex-retomada.log`).
3. **17/09, 09:13–09:57, orquestrador.** Executou os mesmos passos de
   `produzir_certificados.py` via `produzir_claude.sh`. O pipeline não foi
   alterado, e os hashes da matriz e do cache foram conferidos antes.
   Resultado em `produzir-claude.log`.

## Proveniência

| item | valor |
|---|---|
| matriz | `build/spectrum/rt-A-6x18.dat`, sha256 `a3335d9550f9567db287c1ae105674a5914446f9d3772d3955e3628a747c5846` |
| polinômio | `polinomio-novo.txt`, sha256 `1c69fd380acbb94cc352261fc26891db99053e91716d88ece9358f08a3563d81`, grau 333, mônico |
| conferência | `det(A+zI)` por Bareiss = polinômio, igualdade racional exata em `z=1/32` e `z=1/10` |

## Certificados

| contorno | parâmetros | nós | classificam | após coarsen | winding | kernel |
|---|---|---:|---:|---:|---:|---|
| `contour-A-6x18-reduced-B-sigma1_64` | `[1/64,1]×[1/4,3/4]`, h=640, v=333 | 1946 | 1946 | 87 | 1 | aceito |
| `contour-A-6x18-B-root-box` | `[1/32,1/10]×[23/50,1/2]`, h=128, v=97 | 450 | 450 | 11 | 1 | aceito |

Ambos fecharam na primeira tentativa, sem refinamento. Diagnóstico em ponto
flutuante (não entra no certificado): a raiz está em `0.063467476+0.479488325i`.

## Reprodução

```bash
make verify-reduced-b-6x18
```

## Alcance

Estes são certificados sobre a MATRIZ FINITA 6×18. Eles não afirmam nada
sobre o operador infinito, sobre as constraints ou sobre a contagem física.
