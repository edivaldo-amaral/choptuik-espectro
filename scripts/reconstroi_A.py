#!/usr/bin/env python3
"""Nivel 2 (amostra): reconstroi a matriz de referencia A a partir dos dados de RT, com os parametros gravados no
schur.json (Z, F, deflacao, phi, vetores), e compara com a A.npy de referencia: sha256 do arquivo .npy e diferenca
maxima entrada a entrada. Confere tambem os dados da deflacao (defl_refA) gravados no schur.json.

Uso: .venv/bin/python scripts/reconstroi_A.py --schur <schur.json> --ref <A.npy> [--saida <A_nova.npy>]"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rouche_L_Z as rz


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--schur', type=Path, required=True)
    ap.add_argument('--ref', type=Path, required=True)
    ap.add_argument('--saida', type=Path, default=None)
    a = ap.parse_args()
    t0 = time.time()
    info = json.loads(a.schur.read_text())
    Zs = 'x'.join(map(str, info['Z'])); Fs = 'x'.join(map(str, info['F']))
    print(f"parametros: Z = {Zs}, F = {Fs}, deflacao = {info['deflacao']}, phi = {info['phi']}, vetores = {info['vetores']}", flush=True)
    A, z, extra = rz.monta_A(Zs, Fs, info['deflacao'], info['phi'], info['vetores'])
    print(f'A montada: n = {z["n"]} [{time.time() - t0:.0f}s]', flush=True)
    for k, v in extra.items():
        g = info['defl_refA'].get(k)
        print(f'  defl_refA.{k}: novo {v} | gravado {g} | {"igual" if v == g else "DIFERE"}')
    saida = a.saida or Path(tempfile.mkdtemp())/'A_nova.npy'
    np.save(saida, A)
    sh_novo, sh_ref = sha256(saida), sha256(a.ref)
    R = np.load(a.ref, mmap_mode='r')
    assert R.shape == A.shape, (R.shape, A.shape)
    dmax = 0.0; amax = 0.0
    for i in range(0, A.shape[0], 512):
        bl = np.asarray(R[i:i + 512])
        dmax = max(dmax, float(np.max(np.abs(A[i:i + 512] - bl))))
        amax = max(amax, float(np.max(np.abs(bl))))
    print(f'sha256 nova {sh_novo[:16]}  referencia {sh_ref[:16]}  -> {"IDENTICA" if sh_novo == sh_ref else "difere"}')
    print(f'max |A_nova - A_ref| = {dmax:.3e}  (max |A| = {amax:.3e}, relativo {dmax/amax:.2e}) [{time.time() - t0:.0f}s]')
    json.dump(dict(sha256_nova=sh_novo, sha256_ref=sh_ref, identica=sh_novo == sh_ref, dmax=dmax, amax=amax,
                   parametros=dict(Z=Zs, F=Fs, deflacao=info['deflacao'], phi=info['phi'], vetores=info['vetores'])),
              open(saida.with_suffix('.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
