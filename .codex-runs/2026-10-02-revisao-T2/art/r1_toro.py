#!/usr/bin/env python3
"""R1 (revisao T2): checagem INDEPENDENTE do lado do toro do L-real.

1. ||W J_m^-1 P_{n>=N} W^-1|| (norma pesada do toro, m fixo) calculada pelo operador EXATO
   (recursao triangular O(n), sem formar matriz) + Lanczos (svds) numa truncagem grande,
   comparada com closed_form.cauda_Rinv (cota de Schur usada em l_real_toro.py).
2. ||W J_m^-1 W^-1|| (modo inteiro) para |m| >= 400, contra (1 + 2 c_w/(k2 - 1))/(|m|/2).
3. Estrutura de J: confere que signed_operator_L.J_livre e exatamente a forma usada aqui
   (diagonal mu(n+1) + s + i m/2; 2 mu n_c em TODAS as linhas j < c da lista da componente).

O peso kappa1^m cancela dentro de um modo porque J e diagonal em m (bloco).
Truncagem: J e triangular superior em n, logo o bloco lider de J^-1 e exato; a norma truncada
cresce monotonamente com a truncagem e converge para a norma infinita (cota inferior).
"""
import math
import sys
import json
import numpy as np
from scipy.sparse.linalg import LinearOperator, svds

sys.path.insert(0, 'scripts')
import closed_form as cf
import signed_operator_L as sl

MU = sl.MU
K2 = 5/4


def razoes(ns, k2=K2):
    """w_{ns[k]}/w_{ns[k+1]}, w = sqrt(omega2(n)), sem overflow."""
    def logw(n):
        n = np.asarray(n, float)
        return np.where(n == 0, 0.0, n*math.log(k2) + 0.5*np.log1p(k2**(-4.0*n)))
    lw = logw(ns)
    return np.exp(lw[:-1] - lw[1:])


class Modo:
    def __init__(self, ns, sigma, mu=MU):
        self.ns = np.asarray(ns)
        self.d = mu*(self.ns + 1.0) + sigma
        self.v = 2*mu*self.ns.astype(float)
        self.r = razoes(self.ns)          # r[k] = w_k / w_{k+1}
        self.n = len(self.ns)

    def solve(self, b):
        """y~ = W J^-1 W^-1 b~ (recursao para tras)."""
        n = self.n; d = self.d; v = self.v; r = self.r
        y = np.zeros(n, complex)
        St = 0.0 + 0.0j                   # S~_j = w_j sum_{c>j} v_c y_c
        y[n-1] = b[n-1]/d[n-1]
        for j in range(n - 2, -1, -1):
            St = r[j]*(St + v[j+1]*y[j+1])
            y[j] = (b[j] - St)/d[j]
        return y

    def solve_adj(self, b):
        """y~ = W^-1 J^-* W b~ (recursao para frente): J^* = conj(D) + U^T, (U^T z)_c = v_c sum_{j<c} z_j."""
        n = self.n; dc = np.conj(self.d); v = self.v; r = self.r
        y = np.zeros(n, complex)
        Tt = 0.0 + 0.0j                   # T~_c = (1/w_c) sum_{j<c} y_j (y nao pesado)
        y[0] = b[0]/dc[0]
        for c in range(1, n):
            Tt = r[c-1]*(Tt + y[c-1])
            y[c] = (b[c] - v[c]*Tt)/dc[c]
        return y


def norma(modo: Modo, N0=None, k=1):
    """sigma_max de W J^-1 W^-1 restrito as colunas com n >= N0 (None: todas)."""
    n = modo.n
    mask = np.ones(n, bool) if N0 is None else (modo.ns >= N0)

    def mv(x):
        x = np.asarray(x).ravel().astype(complex)
        x = np.where(mask, x, 0)
        return modo.solve(x)

    def rmv(x):
        x = np.asarray(x).ravel().astype(complex)
        return np.where(mask, modo.solve_adj(x), 0)

    A = LinearOperator((n, n), matvec=mv, rmatvec=rmv, dtype=complex)
    s = svds(A, k=k, return_singular_vectors=False, tol=1e-10, maxiter=5000)
    return float(max(s))


def checa_adjunto(modo: Modo, rng):
    x = rng.standard_normal(modo.n) + 1j*rng.standard_normal(modo.n)
    y = rng.standard_normal(modo.n) + 1j*rng.standard_normal(modo.n)
    a = np.vdot(y, modo.solve(x)); b = np.vdot(modo.solve_adj(y), x)
    return abs(a - b)/abs(a)


def checa_Jlivre(rng):
    """A recursao inverte a J_livre do codigo (pesada) numa truncagem pequena."""
    pior = 0.0
    for m in (0, 1, 2, -3):
        r = sl.Radial(m, 60)
        J = sl.J_livre(r, 0.3 + 0.2j)
        for c, ns in r.comps:
            o0, o1 = r.off[c]
            Jc = J[o0:o1, o0:o1]
            w = np.sqrt(cf.omega(ns, K2))
            mo = Modo(ns, 0.3 + 0.2j + 1j*m/2)
            b = rng.standard_normal(len(ns)) + 0j
            y = mo.solve(b)
            x = y/w                           # x = J^-1 W^-1 b
            pior = max(pior, np.abs(Jc @ x - b/w).max()/np.abs(b/w).max())
    return pior


def main():
    rng = np.random.default_rng(1)
    out = dict()
    out['J_livre_recursao_erro_rel'] = checa_Jlivre(rng)
    print('recursao inverte J_livre: erro relativo', out['J_livre_recursao_erro_rel'], flush=True)
    N = 4096
    cw = math.sqrt(1 + K2**-4)
    linhas = []
    casos = [(0, 1), (0, 2), (1, 1), (2, 1), (2, 2), (-7, 1), (40, 1), (-150, 2), (399, 1)]
    for m, passo in casos:
        for extra in (3000, 12000):
            Ntot = N + extra
            ns = np.arange(Ntot) if passo == 1 else np.arange(1, Ntot, 2)
            mo = Modo(ns, 1j*m/2)
            if extra == 3000:
                eadj = checa_adjunto(mo, rng)
            v = norma(mo, N)
            cota = cf.cauda_Rinv(1j*m/2, N, passo, MU, K2, N + 2500)
            linhas.append(dict(m=m, passo=passo, Ntot=Ntot, norma_exata_truncada=v, cota_cauda_Rinv=cota,
                               razao=v/cota, erro_adjunto=eadj))
            print(f'm={m:5d} passo={passo} Ntot={Ntot}: ||J^-1 P_(n>=N)|| = {v:.6f}   cota cauda_Rinv = {cota:.6f}'
                  f'   razao {v/cota:.3f}   (adj {eadj:.1e})', flush=True)
    out['cauda'] = linhas
    inteiro = []
    for m in (400, -400, 401, 700, 1023):
        for passo in ((1, 2) if m % 2 == 0 else (1,)):
            Ntot = 6000
            ns = np.arange(Ntot) if passo == 1 else np.arange(1, Ntot, 2)
            v = norma(Modo(ns, 1j*m/2))
            unif = (1 + 2*cw/(K2 - 1))/(abs(m)/2)
            inteiro.append(dict(m=m, passo=passo, Ntot=Ntot, norma_modo_inteiro=v, cota_uniforme=unif, razao=v/unif))
            print(f'm={m:5d} passo={passo}: ||J_m^-1|| (modo inteiro, trunc {Ntot}) = {v:.6f}   cota uniforme = {unif:.6f}'
                  f'   razao {v/unif:.3f}', flush=True)
    out['modo_inteiro'] = inteiro
    json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else 'r1_toro.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
