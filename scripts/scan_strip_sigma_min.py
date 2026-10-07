#!/usr/bin/env python3
"""Varredura fina de sigma_min(H) na faixa 0<Re s<1/8, malha 6x18.

sigma_min(H) = 1/||V||. Onde sigma_min mergulha, ha espectro. A faixa esta
listada como DESCOBERTA em ALGEBRAIC_TAIL e em C4; o pico de ||V|| em
sigma0=1/32 sugeriu que nao esta vazia. Aqui se localiza.
"""
import sys, time
MALHA = sys.argv[1] if len(sys.argv)>1 else '6x18'
from fractions import Fraction as Q
sys.path.insert(0, 'scripts')
sys.path.insert(0, '<scratch>')
from sweep_V import norma_V

print(f"{'Re s':>8} {'||V||':>12} {'sigma_min':>12}")
print('-'*36)
for k in range(1, 13):
    s = Q(k, 100)
    t = time.time()
    dim, V, eta = norma_V(MALHA, s)
    if V is None:
        print(f"{float(s):>8.3f} {'DIVERGIU':>12} {'':>12}  ({time.time()-t:.0f}s)")
    else:
        print(f"{float(s):>8.3f} {float(V):>12.1f} {1/float(V):>12.3e}  ({time.time()-t:.0f}s)")
