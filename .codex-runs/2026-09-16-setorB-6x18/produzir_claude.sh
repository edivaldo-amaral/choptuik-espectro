#!/usr/bin/env bash
# Retomada pelo orquestrador (Claude) em 17/09: mesmos passos de produzir_certificados.py,
# com logs novos. Pipeline inalterado; cache já conferido contra dois Bareiss exatos.
set -uo pipefail
cd <repo>
run=.codex-runs/2026-09-16-setorB-6x18
py=.venv/bin/python
export OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONINTMAXSTRDIGITS=0
mat=build/spectrum/rt-A-6x18.dat; cache=$run/polinomio-novo.txt
[ "$(sha256sum < $mat | cut -d' ' -f1)" = a3335d9550f9567db287c1ae105674a5914446f9d3772d3955e3628a747c5846 ] || { echo "ERRO: hash da matriz"; exit 1; }
[ "$(sha256sum < $cache | cut -d' ' -f1)" = 1c69fd380acbb94cc352261fc26891db99053e91716d88ece9358f08a3563d81 ] || { echo "ERRO: hash do cache"; exit 1; }
produzir() { # nome nomeLean h v parametros...
  local nome=$1 lname=$2 h=$3 v=$4; shift 4
  for t in 0 1 2 3 4 5 6; do
    local out=$run/$nome-claude-t$t-nodes.txt log=$run/$nome-claude-t$t.log
    echo "[$(date +%T)] $nome tentativa $t: h=$h v=$v"
    (ulimit -v 2500000; exec nice -n 10 $py scripts/build_contour.py $mat "$@" --horizontal $h --vertical $v --bits 64 \
      --name $lname --label 6x18 --cache $cache --out $out) > $log 2>&1
    local c=$?
    echo "[$(date +%T)] $nome tentativa $t: codigo $c"
    if [ $c -eq 0 ]; then
      grep -q '^target 1$' $out || { echo "ERRO: target != 1 em $out"; tail -5 $log; return 1; }
      [ -e examples/$nome-nodes.txt ] && { echo "ERRO: examples/$nome-nodes.txt ja existe"; return 1; }
      cp $out examples/$nome-nodes.txt
      nice -n 10 $py scripts/make_winding_certificate.py examples/$nome-nodes.txt --coarsen --standalone --quiet \
        --lean $run/$nome-claude.lean --increments examples/$nome-increments.txt > $run/$nome-claude-emissao.log 2>&1 \
        || { echo "ERRO: emissao"; tail -5 $run/$nome-claude-emissao.log; return 1; }
      echo "[$(date +%T)] $nome: certificado emitido"
      return 0
    fi
    [ $c -eq 2 ] || { echo "ERRO inesperado"; tail -15 $log; return 1; }
    h=$((h*2)); v=$((2*v+1))
  done
  echo "ERRO: $nome nao fechou em 7 refinamentos"; return 1
}
produzir contour-A-6x18-B-root-box setorB6x18Raiz 128 97 --sigma0 1/32 --radius 1/10 --height 1/50 --imag-center 12/25 || exit 1
produzir contour-A-6x18-reduced-B-sigma1_64 setorB6x18Grande 640 333 --sigma0 1/64 --radius 1 --height 1/4 --imag-center 1/2 || exit 1
echo "[$(date +%T)] verificacao final"
bash scripts/verify_reduced_b_6x18.sh > $run/verificacao-final-claude.log 2>&1; c=$?
tail -3 $run/verificacao-final-claude.log; echo "VERIFICACAO codigo $c"
