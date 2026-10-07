#!/usr/bin/env python3
"""||V|| em funcao da borda sigma0 do contorno, na malha 6x18.

A explosao de ||V|| e governada pela distancia da borda a raiz em s=mu.
sigma0 e uma ESCOLHA: 1/8 nao e obrigatorio. Mede-se o premio de afasta-la.
"""
import sys, time
from fractions import Fraction as Q
sys.path.insert(0, 'scripts')
sys.path.insert(0, '<scratch>')
from sweep_V import norma_V

RAIZ = 0.144011          # raiz do 6x18; converge para mu=0,168311
MU = 0.168311
print(f"{'sigma0':<12} {'||V||':>12} {'d(raiz 6x18)':>14} {'d(mu)':>10} {'eta':>11}")
print('-'*64)
for nome, s in [('1/8', Q(1,8)), ('1/10', Q(1,10)), ('1/16', Q(1,16)),
                ('1/32', Q(1,32)), ('1/64', Q(1,64)), ('1/256', Q(1,256))]:
    t = time.time()
    dim, V, eta = norma_V('6x18', s)
    vs = f"{float(V):>12.1f}" if V is not None else f"{'DIVERGIU':>12}"
    print(f"{nome:<12} {vs} {abs(float(s)-RAIZ):>14.4f} {abs(float(s)-MU):>10.4f}"
          f" {eta:>11.2e}  ({time.time()-t:.0f}s)")
