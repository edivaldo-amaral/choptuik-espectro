# Tarefa para o Codex: certificado de winding do setor B em 6×18

Diretório: <repo>. NÃO é um repositório git; não rode git.
Idioma dos arquivos: português.

## Contexto (mínimo)

O projeto estuda o espectro do pencil `L_A(s) = L_A(0) + s I` da solução de
Choptuik. A contagem de raízes numa região usa o princípio do argumento: o
winding de `det(A + zI)` num contorno, com certificado de caixas racionais
conferido pelo KERNEL do Lean (`formal/ChoptuikFormal/WindingInterval.lean`).
Por redução, o setor B corresponde à faixa `1/4 <= Im s <= 3/4` de `L_A`.
Há um certificado para 4×12 (`make verify-reduced-b`, retângulo
`[1/8,1]×[1/4,3/4]`, winding 0). A auditoria de 16/09 mostrou que esse zero é
artefato da malha: em 6×18 (`build/spectrum/rt-A-6x18.dat`, dimensão 333) há
exatamente uma raiz na faixa, perto de `0.063467+0.479488i`, com
`Re s < 1/8`.

## O que fazer

Produza, com o pipeline EXISTENTE e SEM MODIFICÁ-LO
(`scripts/build_contour.py` e `scripts/make_winding_certificate.py`), dois
certificados de winding para a matriz 6×18, ambos aceitos pelo kernel:

1. **Retângulo grande** `∂([1/64, 1] × [1/4, 3/4])`: winding esperado **1**.
2. **Retângulo pequeno** com cantos racionais contendo a raiz
   `0.063467+0.479488i`, por exemplo `[1/32, 1/10] × [23/50, 1/2]`
   (ajuste se preciso, mantendo-o longe de outras raízes): winding esperado **1**.

Leia `Makefile` (alvo `verify-reduced-b`), `docs/SECTOR_B.md` e
`docs/WINDING_DATA.md` §2 para ver como o 4×12 foi gerado: parâmetros,
`--imag-center`, regras de número ímpar de subdivisões verticais e
`--coarsen`. Se uma aresta não classificar, REFINE o contorno; nunca afrouxe
critérios.

## Regras duras

- **Não edite, renomeie nem apague NENHUM arquivo existente.** Há um manifesto
  sha256 de todos os arquivos em `.snapshots/2026-09-16-pre-codex/`, e
  qualquer alteração será detectada e o trabalho descartado. Só crie arquivos
  novos, e SOMENTE estes:
  - `examples/contour-A-6x18-reduced-B-sigma1_64-nodes.txt` e `...-increments.txt`
  - `examples/contour-A-6x18-B-root-box-nodes.txt` e `...-increments.txt`
  - `scripts/verify_reduced_b_6x18.sh`: reproduz as duas verificações do jeito
    de `make verify-reduced-b` (verify_winding_certificate.py com o winding
    esperado, e depois o Lean emitido checado com `lake env lean` em `formal/`)
    e sai com código ≠ 0 se algo falhar
  - `.codex-runs/2026-09-16-setorB-6x18/RELATORIO.md`
  - arquivos temporários só dentro de `.codex-runs/2026-09-16-setorB-6x18/`
- **Cache de polinômio:** a leitura de `--cache` só confere dimensão e
  monicidade, então não prova que o polinômio veio desta matriz. NÃO use cache
  pré-existente. Se criar um nesta execução, grave junto o sha256 da matriz de
  origem e, antes de confiar nele, confira `det(A + zI)` por Bareiss exato
  (Fraction) em pelo menos 2 pontos racionais distintos contra o polinômio.
  Registre os pontos e a igualdade no relatório.
- **Recursos:** a máquina tem 4 núcleos e outro processo pesado já usa 2. Use
  `OMP_NUM_THREADS=2` e `nice -n 10`. O cálculo modular pode levar mais de uma
  hora: rode-o com `nohup` para um log em `.codex-runs/...` e acompanhe o log,
  em vez de esperar em foreground num comando com timeout.
- Nada de `sorry`, `axiom`, `admit` ou `native_decide` em Lean emitido. Não
  edite `formal/`.
- **Honestidade:** isto é um certificado sobre uma MATRIZ FINITA 6×18. Não
  escreva em lugar nenhum que prova algo sobre o operador infinito, sobre
  constraints ou sobre a contagem física. Se não fechar, diga que não fechou e
  por quê.

## RELATORIO.md (curto)

- Comandos exatos usados, com parâmetros.
- sha256 da matriz e checagem Bareiss×polinômio (pontos, valores iguais).
- Para cada contorno: número de nós e de arestas, arestas "none" (deve ser 0),
  winding, e saída final do `lake env lean`.
- Saída de `bash scripts/verify_reduced_b_6x18.sh`.
- Tempo total e o que ficou incompleto, se algo ficou.
