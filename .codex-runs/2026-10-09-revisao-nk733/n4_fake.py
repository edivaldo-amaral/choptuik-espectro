#!/usr/bin/env python3
"""Checkpoint 'falso' para nk_L3_T.py: todas as janelas de T3.json MENOS as escolhidas, no formato de feitos,
para que a reexecução (--juntar este arquivo, --so-janelas) calcule só as escolhidas."""
import json, sys
from pathlib import Path
T = json.loads(Path(sys.argv[1]).read_text()); escolhidas = sys.argv[3:]
feitos = {}
for i, nm in enumerate(T['nome']):
    if nm in escolhidas:
        continue
    lin = list(T['N'][i]); rho = lin[i]; lin[i] = 0.0
    feitos[nm] = dict(N=lin, TZ=T['TZ'][i], nV=T['nV'][i], nVhi=T['nVhi'][i], rho=rho, Y=T['Y'][i], K=T['K'][i], nQm=T['nQm'][i])
Path(sys.argv[2]).write_text(json.dumps(feitos))
print(len(feitos), 'janelas no checkpoint falso; a calcular:', escolhidas)
