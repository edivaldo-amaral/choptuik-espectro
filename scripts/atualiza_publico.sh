#!/bin/bash
# Atualiza o clone publico (~/choptuik-publico, remoto edivaldo-amaral/choptuik-espectro) a partir do repositorio de
# trabalho: exporta a arvore limpa para um diretorio temporario e sincroniza (sem tocar no .git do clone), depois
# mostra o diff para revisao. O commit e o push ficam para o autor confirmar.
set -euo pipefail
cd "$(dirname "$0")/.."
CLONE=${1:-$HOME/choptuik-publico}
TMP=$(mktemp -d)
.venv/bin/python scripts/exporta_publico.py --destino "$TMP/arvore"
rsync -a --delete --exclude .git "$TMP/arvore/" "$CLONE/"
rm -rf "$TMP"
git -C "$CLONE" status --short | head -50
echo "Revise e, se estiver bom: git -C $CLONE add -A && git -C $CLONE commit -m '...' && git -C $CLONE push"
