#!/bin/bash
# Baixa os dados publicados de Reiterer & Trubowitz (arXiv:1203.3766v1, "Choptuik's critical spacetime exists",
# Comm. Math. Phys. 368 (2019)) e confere os sha256. O repositorio NAO redistribui esses arquivos.
# Uso: bash scripts/baixa_dados_rt.sh   (cria .cache/rt-1203.3766v1/)
set -euo pipefail
cd "$(dirname "$0")/.."
D=.cache/rt-1203.3766v1
TAR=$D/arXiv-1203.3766v1.tar.gz
mkdir -p "$D"
SHA_TAR=3001c85ab701b60262aa1de7e2f2511e6780d3d7f27cdbcfa9429994423e4295
if [ ! -f "$TAR" ]; then
  curl -sSL -o "$TAR" "https://arxiv.org/src/1203.3766v1"
fi
echo "$SHA_TAR  $TAR" | sha256sum -c -
TMP=$(mktemp -d); tar xzf "$TAR" -C "$TMP" anc
cp -r "$TMP"/anc/. "$D"/; rm -rf "$TMP"
for par in "sourcecode/RefA.dat 5776e6038f994e06" "RefAplusB.dat db50272899ebff4a"; do
  set -- $par
  h=$(sha256sum "$D/$1" | cut -c1-16)
  [ "$h" = "$2" ] && echo "ok  $1" || { echo "FALHA $1: $h != $2"; exit 1; }
done
echo "dados de RT prontos em $D"
