#!/usr/bin/env bash
# Cadeia longa: gera malhas cada vez maiores e verifica a rota 5 em cada uma.
#
# O limiar de contracao da rota 5 e n~51, e a zona de borda do truncamento
# contamina ~8 indices radiais, entao so uma malha com N bem acima de 60
# permite ler o limiar sem contaminacao. Esta cadeia sobe ate onde o fundo
# publicado permitir, e registra onde parou.
#
# Duas armadilhas ja pagas em 21/09, e blindadas aqui:
#  1. build_spectral.sh recompila o binario. Se um gerador anterior ainda o
#     estiver executando, a recompilacao corrompe a execucao e o processo
#     aborta em Field_shrink -- o que parece limite de truncamento e nao e.
#     Compila-se UMA vez, e matam-se geradores orfaos antes de comecar.
#  2. encanar a saida do gerador para grep troca o status de saida pelo do
#     grep, e uma falha vira "PRONTA". Aqui o status vem do proprio gerador.
set -u
cd <repo>
LOG=build/spectrum/noite-$(date +%Y%m%d).log
exec >> "$LOG" 2>&1
echo "=== inicio $(date '+%F %T') ==="

pkill -f choptuik-spectral 2>/dev/null && sleep 3
echo "[$(date +%T)] compilando o binario uma unica vez ..."
./scripts/build_spectral.sh > /dev/null 2>&1 || { echo "compilacao falhou"; exit 1; }
echo "[$(date +%T)] binario pronto"

for mn in "13 39" "14 42" "16 48" "18 54" "20 60" "24 72" "30 90"; do
  set -- $mn
  ALVO="build/spectrum/rt-A-${1}x${2}.dat"
  if [ -s "$ALVO" ]; then echo "[$(date +%T)] ${1}x${2} ja existe, pulando"; continue; fi
  rm -f "$ALVO"
  echo "[$(date +%T)] gerando ${1}x${2} ..."
  S=$(date +%s)
  OMP_NUM_THREADS=2 ./scripts/generate_spectrum.sh "$1" "$2" > /tmp/gen-passo.log 2>&1
  RC=$?
  E=$(date +%s)
  if [ "$RC" -eq 0 ] && [ -s "$ALVO" ]; then
    echo "[$(date +%T)] ${1}x${2} PRONTA em $((E-S))s  ($(du -m "$ALVO" | cut -f1) MB)"
    echo "[$(date +%T)] verificando ${1}x${2} ..."
    .venv/bin/python scripts/verify_apriori_regime.py "${1}x${2}" || echo "  verificacao falhou"
  else
    echo "[$(date +%T)] ${1}x${2} FALHOU (status $RC) apos $((E-S))s"
    grep -iE "assertion|error|erro" /tmp/gen-passo.log | head -3
    rm -f "$ALVO"
  fi
  sleep 5
done
echo "=== fim $(date '+%F %T') ==="
