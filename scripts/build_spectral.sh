#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cache_dir=$("$root_dir/scripts/fetch_rt.sh")
source_dir="$cache_dir/sourcecode"
marker="$source_dir/.choptuik-spectral-v1"

if [[ ! -f "$marker" ]]; then
  # A saida do `patch` vai para stderr: a UNICA linha de stdout deste script e o
  # caminho do binario, consumido por `binary=$(build_spectral.sh)` em
  # generate_spectrum.sh. Sem este redirecionamento, o "patching file ..." entra
  # na substituicao de comando e a primeira execucao em checkout limpo falha
  # (a segunda funciona, porque o marcador ja existe e o patch nao roda).
  patch -d "$source_dir" -p1 < "$root_dir/patches/rt-spectral.patch" >&2
  touch "$marker"
fi

gcc -O3 -W -Wall -pedantic -std=c99 -fopenmp -fcommon \
  -o "$source_dir/choptuik-spectral" \
  "$source_dir/choptuik.c" "$source_dir/DyadicQ.c" \
  "$source_dir/Sector.c" "$source_dir/GI.c" "$source_dir/Field.c" \
  "$source_dir/MultiField.c" "$source_dir/Matrix.c" \
  "$source_dir/IndexDOF.c" "$source_dir/Field_approximate.c" \
  "$source_dir/MultiField_approximate.c" -lgmp

printf '%s\n' "$source_dir/choptuik-spectral"

