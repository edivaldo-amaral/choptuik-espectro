#!/usr/bin/env python3
# DIAGNOSTICO EM PONTO FLUTUANTE -- NAO E PROVA. Ver S3_CONDICIONAMENTO_SHARP.md 10.2.
"""Cauda de Fourier, modo a modo: composicao radial 1D contra produto de normas."""
import sys; sys.path.insert(0, 'scripts')
import numpy as np
from pathlib import Path
from component_exterior import read_fields
MU = 722873400/2**32; K2 = 9/8; NR = 300
F = read_fields(Path('.cache/rt-1203.3766v1/sourcecode/RefA.dat'))
def radial(j, d):
    """coeficientes simetricos (n>=0) do modo de Fourier d do campo j (complexo)."""
    g = np.zeros(101, complex)
    for (m, n, part), v in F[j].items():
        if m == d: g[n] += float(v)*(1j if part else 1)
    return g
def Pn(g):  # reflexao: n -> (-1)^n
    return g*(-1.0)**np.arange(len(g))
def mult(g, nin, nout):
    """matriz de multiplicacao por g (simetrico) de entradas n in nin para saidas nout."""
    G = lambda k: g[abs(k)] if abs(k) < len(g) else 0.0
    M = np.zeros((len(nout), len(nin)), complex)
    for c, n1 in enumerate(nin):
        for r, n in enumerate(nout):
            M[r, c] = G(n - n1) + (G(n + n1) if n1 > 0 else 0.0)
    return M
def livre(ns, passo, sigma):
    J = np.diag(MU*(ns + 1.0)).astype(complex) + sigma*np.eye(len(ns))
    for b in range(len(ns)):
        J[:b, b] += 2*MU*ns[b]
    return J
n1 = np.arange(1, NR, 2); n2 = np.arange(NR)                       # c1 impar, c2 todos
v = [radial(j, 0) for j in range(3)]; P = [Pn(x) for x in v]
v1, v2, v3 = v; P1, P2, P3 = P
F1 = (v2+P2-P3-v3, v1-v2+P2); F2 = (-v1+P3-v3, -v1+v2+3*P2); F3 = (v2+P2+P3+v3, -v1+v2+3*P2)
ev = (n2 % 2 == 0); od = ~ev
def bloco(Fk, nout):   # Gamma: (1/2)F1 x1 + (1/2)F2 E x2 + (1/2)F3 O x2
    A = 0.5*mult(Fk[0], n1, nout)
    Bm = 0.5*mult(Fk[1], n2, nout)*ev[None, :] + 0.5*mult(Fk[2], n2, nout)*od[None, :]
    return A, Bm
B11, B12 = bloco((F1[0], F2[0], F3[0]), n1)
B21, B22 = bloco((F1[1], F2[1], F3[1]), n2)
Rx = np.zeros((len(n2), len(n1)))
for c, n in enumerate(n1):
    k = (n - 1)//2
    for jj in range(k + 1):
        Rx[2*jj, c] = 2*(-1)**(k - jj)
B0 = np.block([[B11, B12], [B21 + MU*Rx, B22]])
om = lambda n: np.where(n == 0, 1.0, K2**(2.0*n) + K2**(-2.0*n))
dv = np.sqrt(np.concatenate([om(n1), om(n2)]))
W = lambda X: (dv[:, None]*X)/dv[None, :]
nb = np.linalg.norm(W(B0), 2)
print(f'||B0 (d=0, com mu R_x)|| radial = {nb:.4f}')
for s in (0.0, 0.25):
    for m in (10, 20, 30, 40, 60):
        sig = s + 1j*m/2
        Q = np.block([[np.linalg.inv(livre(n1, 2, sig)), np.zeros((len(n1), len(n2)))],
                      [np.zeros((len(n2), len(n1))), np.linalg.inv(livre(n2, 1, sig))]])
        nq = np.linalg.norm(W(Q), 2); nc = np.linalg.norm(W(B0 @ Q), 2)
        print(f's={s} m={m:3d}: ||R_m^-1|| = {nq:.4f}   ||B0 R_m^-1|| = {nc:.4f}   produto = {nb*nq:.4f}   ganho = {nb*nq/nc:.2f}x')

print('\n--- d != 0: B_d (modo de Fourier d do fundo, sem R_x) composto com R_m^-1 ---')
def Bd(d):
    vv = [radial(j, abs(d)) for j in range(3)]
    if d < 0: vv = [x.conj() for x in vv]
    PP = [Pn(x) for x in vv]; a1, a2, a3 = vv; q1, q2, q3 = PP
    G1 = (a2+q2-q3-a3, a1-a2+q2); G2 = (-a1+q3-a3, -a1+a2+3*q2); G3 = (a2+q2+q3+a3, -a1+a2+3*q2)
    b11, b12 = bloco((G1[0], G2[0], G3[0]), n1); b21, b22 = bloco((G1[1], G2[1], G3[1]), n2)
    return np.block([[b11, b12], [b21, b22]])
Qc = {}
def Qm(m, s):
    if (m, s) not in Qc:
        sig = s + 1j*m/2
        Qc[(m, s)] = np.block([[np.linalg.inv(livre(n1, 2, sig)), np.zeros((len(n1), len(n2)))],
                               [np.zeros((len(n2), len(n1))), np.linalg.inv(livre(n2, 1, sig))]])
    return Qc[(m, s)]
for M in (20, 30):
    s = 0.0
    ms = [M, M+10, M+30]
    tot_comp = tot_prod = 0.0
    linhas = []
    for d in range(-40, 41, 2):
        B = Bd(d) if d else B0
        nbd = np.linalg.norm(W(B), 2)
        if nbd < 1e-9: continue
        comp = max(np.linalg.norm(W(B @ Qm(m, s)), 2) for m in ms)
        prod = nbd*max(np.linalg.norm(W(Qm(m, s)), 2) for m in ms)
        tot_comp += comp; tot_prod += prod
        if abs(d) <= 6: linhas.append(f'd={d:+d}: {comp:.3f} (prod {prod:.3f})')
    print(f'M={M}, s=0: ' + ' | '.join(linhas))
    print(f'   soma sobre d: composicoes {tot_comp:.4f}   produtos {tot_prod:.4f}')
