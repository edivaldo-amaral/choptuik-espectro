Continuação da tarefa em TAREFA.md (mesmas regras duras, mesma lista fechada de arquivos permitidos).

O que aconteceu: em 16/09 às ~20:07 o sistema ficou sem memória e seu processo foi morto, junto com o cálculo modular (parado no primo 500/877, sem checkpoint). Nenhum arquivo existente foi alterado; o manifesto confere.

O que continua válido: os dois determinantes Bareiss exatos em bareiss-ponto-1.json (z=1/32) e bareiss-ponto-2.json (z=1/10), vinculados ao sha256 da matriz. Não os recalcule.

O que fazer agora:
1. Recalcule o polinômio característico exato do zero. Como origem-polinomio.json já existe, renomeie SEU arquivo para origem-polinomio-execucao1.json ou grave com outro nome (só dentro de .codex-runs/2026-09-16-setorB-6x18/). Se puder, escreva uma variante do seu calcular_polinomio.py (arquivo novo, também dentro da pasta da execução) que salve os resíduos por primo em disco a cada 50 primos e retome deles, sem modificar scripts/build_contour.py. Mantenha o uso de memória baixo: a máquina tem 7,6 GB e outro processo pesado roda em paralelo, limitado a ~3 GB.
2. Confira o polinômio contra os dois determinantes Bareiss salvos (igualdade racional exata).
3. Produza os dois certificados (retângulo [1/64,1]×[1/4,3/4] com winding 1; caixa da raiz com winding 1), rode bash scripts/verify_reduced_b_6x18.sh e escreva RELATORIO.md.
Recursos: OMP_NUM_THREADS=2, nice -n 10, nohup para os cálculos longos.
