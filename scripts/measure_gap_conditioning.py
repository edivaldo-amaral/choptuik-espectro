#!/usr/bin/env python3
"""||V|| no VAO entre a raiz de constraint (0,404) e a fisica (0,731), 6x18.

Se os modos nao fisicos forem removidos analiticamente ANTES da contagem, o
contorno so precisa cercar a raiz fisica, e a borda esquerda pode ficar no vao.
Mede-se quanto isso vale em condicionamento.
"""
import sys, time
from fractions import Fraction as Q
sys.path.insert(0, 'scripts')
sys.path.insert(0, '<scratch>')
from sweep_V import norma_V

R_S3, R_FIS = 0.404385, 0.731184      # raizes do 6x18
print(f"vao entre {R_S3} (S3) e {R_FIS} (fisica); meio em {(R_S3+R_FIS)/2:.3f}")
print(f"{'borda s':<10} {'||V||':>10} {'d(S3)':>8} {'d(fis)':>8} {'min':>8}")
print('-'*48)
for nome, s in [('0,45',Q(45,100)),('0,50',Q(1,2)),('0,55',Q(55,100)),
                ('0,568',Q(568,1000)),('0,60',Q(3,5)),('0,65',Q(65,100))]:
    t=time.time(); dim,V,eta = norma_V('6x18', s)
    d1,d2 = abs(float(s)-R_S3), abs(float(s)-R_FIS)
    vs = f"{float(V):>10.1f}" if V is not None else f"{'DIVERGIU':>10}"
    print(f"{nome:<10} {vs} {d1:>8.3f} {d2:>8.3f} {min(d1,d2):>8.3f}  ({time.time()-t:.0f}s)")
