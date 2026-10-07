#!/bin/bash
# Dispara o Codex destacado (setsid: sobrevive ao fim da sessão e à falta de memória).
cd <repo>
D=.codex-runs/2026-09-17-lab-benchmark
setsid nohup codex exec --skip-git-repo-check -C <repo> \
    -s workspace-write -o "$D/codex-saida.jsonl" - < "$D/TAREFA.md" \
    > "$D/codex.log" 2>&1 &
echo "Codex lancado. Acompanhe com: tail -f $D/codex.log"
