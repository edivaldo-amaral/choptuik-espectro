#!/usr/bin/env bash
# Rodada 5a: regioes R1, R3, R2 do Marco 8 (docs/COLUMNWISE_EXTERIOR.md).
# Uso: columnwise_exterior_run.sh <A|sharp> <depth>
set -eu
sys="$1"; depth="$2"
cd "$(dirname "$0")"
ulimit -v 3000000
PY=../.venv/bin/python
nice -n 10 $PY columnwise_phase2.py "$sys" --n-lo 384 --n-hi 606 --M 160 --m-lo 107 --depth "$depth" --workers 1 --rows-per-call 4
nice -n 10 $PY columnwise_phase2.py "$sys" --n-lo 0   --n-hi 384 --M 160 --m-lo 107 --depth "$depth" --workers 1 --rows-per-call 4
nice -n 10 $PY columnwise_phase2.py "$sys" --n-lo 606 --n-hi 900 --M 160 --m-lo 0   --depth "$depth" --workers 1 --rows-per-call 4
echo "FIM $sys"
