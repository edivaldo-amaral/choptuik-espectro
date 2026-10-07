#!/usr/bin/env bash
set -euo pipefail
root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source_dir="$root_dir/.cache/rt-1203.3766v1/sourcecode"
mkdir -p "$root_dir/build"
gcc -O3 -W -Wall -pedantic -std=c11 -fopenmp -fcommon -I "$source_dir" \
  "$root_dir/scripts/rt_sharp_matrix.c" \
  "$source_dir/DyadicQ.c" "$source_dir/Sector.c" "$source_dir/GI.c" \
  "$source_dir/Field.c" "$source_dir/MultiField.c" "$source_dir/Matrix.c" \
  "$source_dir/IndexDOF.c" "$source_dir/Field_approximate.c" \
  "$source_dir/MultiField_approximate.c" -lgmp -o "$root_dir/build/sharp-matrix"
printf '%s\n' "$root_dir/build/sharp-matrix"
