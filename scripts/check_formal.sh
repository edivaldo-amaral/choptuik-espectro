#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
if rg -n --glob '*.lean' '\b(sorry|axiom|admit)\b' "$root_dir/formal"; then
  printf 'ERRO: placeholder ou axioma ad hoc encontrado.\n' >&2
  exit 1
fi

# `native_decide` nao e `sorry` nem `axiom`, mas delega a avaliacao ao
# compilador Lean fora do kernel e sela o resultado num axioma opaco gerado por
# declaracao. Nao falha o build: e contado e reportado, para que a base de
# confianca apareca no ledger em toda execucao. Ver docs/KERNEL_TRUST.md.
native_count=$(python3 "$root_dir/scripts/count_native_decide.py" "$root_dir/formal")
printf 'CONFIANCA: %s uso(s) reais de native_decide (nivel N, fora do kernel).\n' \
  "$native_count"

(cd "$root_dir/formal" && <home>/.elan/bin/lake build)
