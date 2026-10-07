#!/usr/bin/env python3
"""R7 (revisor, no laboratorio2): diagT_lab2.npy (copia local do autor) == diag(T.npy do laboratorio2)? e T e
triangular superior (amostra de blocos)? Grava ~/revisao_B/r7_diagT_saida.txt via stdout."""
import hashlib, numpy as np
from pathlib import Path
D = Path.home()/'Choptuik/build/rouche_L_v2'
T = np.load(D/'T.npy', mmap_mode='r')
d = np.array(np.diagonal(T))
loc = np.load(Path.home()/'revisao_B/diagT_lab2_local.npy')
print('n', len(d), 'diag(T) == diagT_lab2.npy (bit a bit):', np.array_equal(d, loc))
z = -d
em = [complex(x) for x in z if 0 < x.real < 3 and 0.25 < x.imag < 0.75]
print('-T_ii em Omega_B:', em)
perto = sorted([complex(x) for x in z], key=lambda x: min(abs(x.imag-0.25), abs(x.imag-0.75)) if 0 <= x.real <= 3 else 9)[:6]
print('mais perto das arestas horizontais (Re em [0,3]):', perto)
