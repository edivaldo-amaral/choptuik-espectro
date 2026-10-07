#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cache_dir=${RT_CACHE_DIR:-"$root_dir/.cache/rt-1203.3766v1"}
archive="$cache_dir/arXiv-1203.3766v1.tar.gz"
source_dir="$cache_dir/sourcecode"
url=https://arxiv.org/e-print/1203.3766v1
expected_sha256=3001c85ab701b60262aa1de7e2f2511e6780d3d7f27cdbcfa9429994423e4295

mkdir -p "$cache_dir"
if [[ ! -f "$archive" ]]; then
  curl -L --fail --retry 5 -o "$archive.part" "$url"
  mv "$archive.part" "$archive"
fi

actual_sha256=$(sha256sum "$archive" | cut -d' ' -f1)
if [[ "$actual_sha256" != "$expected_sha256" ]]; then
  printf 'ERRO: SHA-256 inesperado para %s\nesperado: %s\nobtido:   %s\n' \
    "$archive" "$expected_sha256" "$actual_sha256" >&2
  exit 1
fi

if [[ ! -f "$source_dir/choptuik.c" ]]; then
  mkdir -p "$source_dir"
  tar -xzf "$archive" -C "$source_dir" --strip-components=2 anc/sourcecode
fi

for data_file in RefAplusB.dat INV_0_0_250_750.dat SHARP_INV_0_0_250_750.dat; do
  if [[ ! -f "$cache_dir/$data_file" ]]; then
    tar -xzf "$archive" -C "$cache_dir" --strip-components=1 "anc/$data_file"
  fi
done

printf '%s\n' "$cache_dir"

