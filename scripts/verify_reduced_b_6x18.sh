#!/usr/bin/env bash
# Confere os dois certificados da matriz finita 6×18 no Python e no kernel Lean.
set -euo pipefail
raiz=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$raiz"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 LEAN_NUM_THREADS=2
export PYTHONDONTWRITEBYTECODE=1 PYTHONINTMAXSTRDIGITS=0

temporarios=$(mktemp -d -t choptuik-setor-b-6x18-XXXXXX)
trap 'rm -rf "$temporarios"' EXIT
python="$raiz/.venv/bin/python"
lake=${LEAN_LAKE:-<home>/.elan/bin/lake}
for nome in contour-A-6x18-reduced-B-sigma1_64 contour-A-6x18-B-root-box; do
    nos="examples/$nome-nodes.txt"
    incrementos="examples/$nome-increments.txt"
    lean="$temporarios/$nome.lean"
    printf '\nVerificando %s (winding esperado: 1)\n' "$nome"
    nice -n 10 "$python" scripts/verify_winding_certificate.py "$incrementos" 1
    nice -n 10 "$python" scripts/make_winding_certificate.py "$nos" \
        --coarsen --standalone --quiet --lean "$lean" \
        --increments "$temporarios/$nome-increments.txt"
    cmp "$incrementos" "$temporarios/$nome-increments.txt"
    if rg -n '\b(sorry|axiom|admit|native_decide)\b' "$lean"; then
        printf 'ERRO: termo proibido no Lean emitido.\n' >&2
        exit 1
    fi
    (cd formal && nice -n 10 "$lake" env lean -j2 "$lean")
    printf 'OK: kernel Lean aceitou %s (código 0).\n' "$nome"
done
printf '\nOK: os dois certificados da matriz finita 6×18 foram aceitos.\n'
