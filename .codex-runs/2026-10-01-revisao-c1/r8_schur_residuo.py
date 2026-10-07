#!/usr/bin/env python3
"""R8 (revisao independente C1): codigo PROPRIO para
  eps_R >= ||A U - U T||_F,  eps_G >= ||I - U^* U||_F,  eps_B = eps_R/sqrt(1 - eps_G) >= ||A - U T U^-1||_2
e a contagem dos -T_ii em Omega = (0,3) x (-1/4,1/4), com a distancia ao bordo.
Diferencas para rouche_L_contagem.residuo: blocos de LINHAS de U em U^*U (outra ordem de soma), T lido por
memmap, cotas de arredondamento por Frobenius GLOBAIS (nao por bloco):
  |fl(XY) - XY| <= g_k |X||Y| (entrada a entrada, produto complexo de comprimento k, g_k = 4(k+1)u/(1-4(k+1)u)),
  || |X||Y| ||_F <= ||X||_F ||Y||_F;  a subtracao fl(P1 - P2) custa u |P1 - P2| <= u(|P1| + |P2|).
Tambem: T triangular superior exata (conferido), e os sha256 de A, T, U."""
import hashlib, json, math, sys, time
from pathlib import Path
import numpy as np

u = 2.0**-53
def g(k):
    kk = 4*(k + 1); return kk*u/(1 - kk*u)
t0 = time.time()
def log(*x):
    print(*x, f'[{time.time()-t0:.0f}s]', flush=True)

d = Path(sys.argv[1]); saida = Path(sys.argv[2])
sha = {}
for nome in ('A.npy', 'T.npy', 'U.npy'):
    h = hashlib.sha256()
    with open(d/nome, 'rb') as fh:
        for bl in iter(lambda: fh.read(1 << 24), b''):
            h.update(bl)
    sha[nome] = h.hexdigest(); log(nome, sha[nome])
A = np.load(d/'A.npy'); n = A.shape[0]
U = np.load(d/'U.npy')
T = np.load(d/'T.npy', mmap_mode='r')
fA = float(np.linalg.norm(A)); fU = float(np.linalg.norm(U))
b = 1024
fT2 = 0.0; tri = True; diag = np.empty(n, complex)
R2 = 0.0; G2 = 0.0; P1f2 = 0.0; P2f2 = 0.0; Gf2 = 0.0
for j0 in range(0, n, b):
    j1 = min(n, j0 + b)
    Tb = np.array(T[:, j0:j1])                     # colunas j0..j1 de T (todas as linhas)
    # abaixo da diagonal: linha i > coluna j  <=>  i > j0 + jj
    for jj in range(j1 - j0):
        if np.any(Tb[j0 + jj + 1:, jj] != 0):
            tri = False
        diag[j0 + jj] = Tb[j0 + jj, jj]
    fT2 += float(np.linalg.norm(Tb))**2
    P1 = A @ U[:, j0:j1]
    P2 = U @ Tb                                    # inclui os zeros abaixo da diagonal (exatos)
    P1f2 += float(np.linalg.norm(P1))**2; P2f2 += float(np.linalg.norm(P2))**2
    P1 -= P2
    R2 += float(np.linalg.norm(P1))**2
    del P1, P2, Tb
    log(f'  AU - UT: colunas {j0}..{j1}')
# I - U^*U por blocos de colunas, mas somando sobre blocos de linhas (ordem diferente da do autor)
for j0 in range(0, n, b):
    j1 = min(n, j0 + b)
    Gb = np.zeros((n, j1 - j0), complex)
    for i0 in range(0, n, 2048):
        i1 = min(n, i0 + 2048)
        Gb += U[i0:i1].conj().T @ U[i0:i1, j0:j1]
    Gb[j0:j1] -= np.eye(j1 - j0)
    Gf2 += float(np.linalg.norm(Gb))**2
    del Gb
    log(f'  I - U*U: colunas {j0}..{j1}')
fT = math.sqrt(fT2)
nsom = n//2048 + 1
# cotas globais de arredondamento (Frobenius)
arR = g(n)*fA*fU + g(n)*fU*fT + u*(math.sqrt(P1f2) + math.sqrt(P2f2))
arG = (g(2048) + nsom*u)*fU*fU*(1 + 1e-6) + u*math.sqrt(Gf2)
eR = math.sqrt(R2)*(1 + 1e-9) + arR
eG = math.sqrt(Gf2)*(1 + 1e-9) + arG
assert eG < 0.5
eB = eR/math.sqrt(1 - eG)*(1 + 1e-9)
z = -diag
dentro = [complex(x) for x in z if 0 < x.real < 3 and -0.25 < x.imag < 0.25]
def dist(x):
    dx = max(-x.real, 0, x.real - 3); dy = max(-0.25 - x.imag, 0, x.imag - 0.25)
    if dx or dy: return math.hypot(dx, dy)
    return min(x.real, 3 - x.real, x.imag + 0.25, 0.25 - x.imag)
perto = sorted(z, key=lambda x: dist(complex(x)))[:8]
out = dict(sha256=sha, n=n, T_triangular=tri, fA=fA, fU=fU, fT=fT, R_fl=math.sqrt(R2), G_fl=math.sqrt(Gf2),
           arred_R=arR, arred_G=arG, eps_R=eR, eps_G=eG, eps_B=eB,
           em_Omega=[[x.real, x.imag] for x in sorted(dentro, key=lambda x: x.real)],
           mais_perto_do_bordo=[[complex(x).real, complex(x).imag, dist(complex(x))] for x in perto])
json.dump(out, open(saida, 'w'), indent=1)
log(json.dumps(out, indent=1))
