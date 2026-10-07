#!/usr/bin/env bash
# Tarefa longa: ||V|| na maior malha disponivel (12x36, dim 1422).
#
# A cadeia de malhas maiores foi abandonada: o fundo publicado RefA nao
# sustenta matriz espectral acima de 12x36 (todas abortam em Field_shrink).
# O que ainda rende e adensar a extrapolacao da rodada 6, que hoje se apoia em
# tres pontos (4x12, 6x18, 8x24), e obter o valor de borda na maior malha.
#
# Cada centro leva ~40-60 min por causa da inversa exata em racionais.
set -u
cd <repo>
LOG=build/spectrum/normas-$(date +%Y%m%d).log
exec >> "$LOG" 2>&1
echo "=== inicio $(date '+%F %T') ==="
echo "objetivo: 4o ponto da saturacao de ||V|| e valor de borda em 12x36"
for centro in "1" "1/8" "0,568:568/1000" "1/10"; do
  nome="${centro%%:*}"; valor="${centro##*:}"
  echo "[$(date +%T)] 12x36, centro $nome ..."
  S=$(date +%s)
  .venv/bin/python - "$valor" <<'PY'
import sys
from fractions import Fraction as Q
sys.path.insert(0, 'scripts')
from measure_parametrix_norm import norma_V
dim, V, eta = norma_V('12x36', Q(sys.argv[1]))
print(f"  dim={dim}  ||V||={float(V) if V else float('nan'):.2f}  defeito={eta:.2e}", flush=True)
PY
  echo "[$(date +%T)] concluido em $(( $(date +%s)-S ))s"
done
echo "=== fim $(date '+%F %T') ==="
