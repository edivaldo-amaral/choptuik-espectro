#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cache_dir=$("$root_dir/scripts/fetch_rt.sh")
source_dir="$cache_dir/sourcecode"
build_dir="$root_dir/build/rt-audit"
mkdir -p "$build_dir"

gcc -O3 -W -Wall -pedantic -std=c99 -fopenmp -fcommon \
  -o "$build_dir/choptuik" \
  "$source_dir/choptuik.c" "$source_dir/DyadicQ.c" \
  "$source_dir/Sector.c" "$source_dir/GI.c" "$source_dir/Field.c" \
  "$source_dir/MultiField.c" "$source_dir/Matrix.c" \
  "$source_dir/IndexDOF.c" "$source_dir/Field_approximate.c" \
  "$source_dir/MultiField_approximate.c" -lgmp

ln -sfn "$cache_dir/RefAplusB.dat" "$source_dir/RefAplusB.dat"
threads=${OMP_NUM_THREADS:-1}
for task in DedicatedGAMMA2 SharpIdentitiesA FailureToBeSolA LoadAplusB; do
  (
    cd "$source_dir"
    OMP_NUM_THREADS="$threads" "$build_dir/choptuik" "$task"
  ) | tee "$build_dir/$task.log"
done

grep -Fq 'operator norm of 2*SBold(Xi(GAMMA2(refA,.)))) <= 21*2^(0)' "$build_dir/DedicatedGAMMA2.log"
grep -Fq 'must be zero' "$build_dir/SharpIdentitiesA.log"
grep -Fq 'l1 norm of SBold(OmegaPlus(refA))) <= 2209*2^(-32)' "$build_dir/FailureToBeSolA.log"
grep -Fq '(abs value of mu) <= 43*2^(-40)' "$build_dir/LoadAplusB.log"

python3 "$root_dir/scripts/verify_rt_data.py" \
  "$cache_dir/INV_0_0_250_750.dat" \
  "$cache_dir/SHARP_INV_0_0_250_750.dat"

printf '%s\n' 'Reprodução parcial concluída.'
printf '%s\n' 'Limite explícito: os dois arquivos EXTERIOR_ESTIMATE não foram publicados no arXiv.'

