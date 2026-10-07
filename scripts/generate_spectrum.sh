#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  printf 'uso: %s NUM_FOURIER NUM_CHEBYSHEV\n' "$0" >&2
  exit 2
fi

num_m=$1
num_n=$2
case "$num_m:$num_n" in
  *[!0-9:]*|0:*|*:0) printf 'os truncamentos devem ser inteiros positivos\n' >&2; exit 2 ;;
esac

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
binary=$("$root_dir/scripts/build_spectral.sh")
source_dir=$(dirname "$binary")
build_dir="$root_dir/build/spectrum"
output="$build_dir/rt-A-${num_m}x${num_n}.dat"
mkdir -p "$build_dir"

if [[ -e "$output" ]]; then
  printf 'a matriz já existe: %s\n' "$output" >&2
  exit 1
fi

threads=${OMP_NUM_THREADS:-1}
(
  cd "$source_dir"
  OMP_NUM_THREADS="$threads" ./choptuik-spectral SpectralMatrix \
    "0_0_${num_m}_${num_n}" "$output"
)

"$root_dir/.venv/bin/python" "$root_dir/scripts/analyze_spectrum.py" "$output"

