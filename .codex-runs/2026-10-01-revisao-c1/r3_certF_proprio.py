#!/usr/bin/env python3
"""R3 (revisao independente C1): certificado PROPRIO do fator de posto baixo F de um centro c.

S = [[(B + J)_ZT Q0_TT(c)], [J_ZT Q0_TT(c)]] (2n x colunas de T acopladas a Z), B = blocos ctx.B (dados).
Para QUALQUER Phi:  ||X S|| <= ||X F|| ||Phi|| + ||X|| ||S - F Phi||   (desigualdade triangular +
submultiplicatividade: X S = X F Phi + X (S - F Phi)).
Aqui Phi_m = minimos quadrados (QR de F, sem Tikhonov), e opcionalmente Tikhonov por QR aumentado.
Codigo proprio para: Q0_m(c) (inversao triangular + cota pelo residuo), montagem de S, Phi, rho e ||Phi||
(Ostrowski + residuo de autovetores). Usa do projeto SO: signed_operator_L (Radial, J_livre), closed_form.omega,
nk_L.Contexto.B (os blocos de fundo, que sao DADOS do certificado) e nk_L.pesos.
Arredondamento (Frobenius, entrada a entrada): |fl(XY) - XY| <= g_k |X||Y|, g_k = 4(k+1)u/(1 - 4(k+1)u) (complexo)."""
import argparse, json, math, sys, time, hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import solve_triangular, qr

sys.path.insert(0, str(Path.home()/'choptuik_trabalho/scripts'))
sys.path.insert(0, str(Path.home()/'Choptuik/scripts'))
import signed_operator_L as sl
import closed_form as cf
import nk_L as nk

u = 2.0**-53
def g(k):
    kk = 4*(k + 1); return kk*u/(1 - kk*u)
MZ, NZ, Mc, Nc, K2 = 32, 128, 200, 400, 5/4
t0 = time.time()
def log(*x):
    print(*x, f'[{time.time()-t0:.0f}s]', flush=True)
fro = lambda X: float(np.linalg.norm(X))

def Q0_cert(m, c):
    """Q^ ~ (J_mu + c)^-1 pesado (sqrt omega) do modo m, n < Nc; dq >= ||Q0 - Q^||_2."""
    r = sl.Radial(m, Nc)
    w = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))
    Jw = (w[:, None]*sl.J_livre(r, c))/w[None, :]
    k = len(Jw)
    Qh = solve_triangular(Jw, np.eye(k, dtype=complex), check_finite=False)
    aJ, aQ = np.abs(Jw), np.abs(Qh)
    E = np.eye(k) - Jw @ Qh
    # erro de E: produto g_k |J||Q| + subtracao u|E| + pesos de Jw (4u |J| |Q|)
    eE = (g(k) + 4*u)*fro(aJ @ aQ)*(1 + 1e-6) + u*fro(E)
    e = fro(E)*(1 + 1e-9) + eE
    assert e < 0.5, e
    nQ = fro(Qh)*(1 + 1e-9)
    return r, Qh, nQ*e/(1 - e), e

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--F', type=Path, required=True)
    ap.add_argument('--c', type=complex, required=True)
    ap.add_argument('--tik', type=float, default=0.0, help='lambda de Tikhonov relativo a ||F||_2 (0 = minimos quadrados)')
    ap.add_argument('--saida', type=Path, required=True)
    a = ap.parse_args()
    c = a.c
    with open(a.F, 'rb') as fh:
        sha = hashlib.sha256(fh.read()).hexdigest()
    log('sha256 F', sha)
    F = np.load(a.F); n2, r = F.shape; n = n2//2
    ms = list(range(-(MZ - 1), MZ)); rsZ = [sl.Radial(m, NZ) for m in ms]
    offs = np.cumsum([0] + [x.dim for x in rsZ]); assert offs[-1] == n
    idx = {m: i for i, m in enumerate(ms)}
    fF = fro(F); fFt = fro(F[:n]); fFb = fro(F[n:])
    # QR de F (aumentado com lambda I para Tikhonov)
    sv = np.linalg.svd(F, compute_uv=False); s1 = float(sv[0])
    lam = a.tik*s1
    log(f'F {F.shape}, sigma {s1:.4e} .. {sv[-1]:.4e}, lambda {lam:.3e}')
    Fa = np.vstack([F, lam*np.eye(r)]) if lam > 0 else F
    Qf, Rf = qr(Fa, mode='economic', check_finite=False); del Fa
    QfH = Qf[:n2].conj().T; del Qf                 # r x 2n (so as linhas de F: o bloco lambda I multiplica 0)
    ctx = nk.Contexto(c, Nc)
    G = np.zeros((r, r), complex); eG = 0.0; rho2 = 0.0; erho2 = 0.0; eS2 = 0.0; qTT = 0.0; ncol = 0
    for m in range(-(min(MZ + nk.BAND_M, Mc) - 1), min(MZ + nk.BAND_M, Mc)):
        rr, Qh, dq, e = Q0_cert(m, c)
        nn = np.concatenate([ns for _, ns in rr.comps])
        cT = np.where(nn >= NZ)[0] if abs(m) < MZ else np.arange(rr.dim)
        QT = Qh[np.ix_(cT, cT)]
        # ||Q0_TT|| <= ||Q^_TT||_F + dq  (Frobenius: grosseiro mas valido; o certF usa espectral)
        qTT = max(qTT, float(np.linalg.norm(QT, 2))*(1 + 1e-9) + dq)   # espectral em ponto flutuante + dq (estimativa)
        Y = np.zeros((n, len(cT)), complex)
        for mo in ms:
            if abs(mo - m) <= nk.BAND_M:
                i = idx[mo]
                Y[offs[i]:offs[i+1]] = ctx.B(mo, m, NZ, Nc)[:, cT]
        Jm = None
        if m in idx:
            i = idx[m]
            w = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in rr.comps]))
            Jfull = (w[:, None]*sl.J_livre(rr, 0.0))/w[None, :]
            lz = np.where(nn < NZ)[0]
            Jm = Jfull[np.ix_(lz, cT)]
            Y[offs[i]:offs[i+1]] += Jm
        aQT = np.abs(QT)
        top = Y @ QT
        eS = fro(Y)*dq + (g(len(cT)) + 6*u)*fro(np.abs(Y) @ aQT)*(1 + 1e-6)
        S = np.zeros((n2, len(cT)), complex); S[:n] = top; del top, Y
        if Jm is not None:
            bot = Jm @ QT
            eb = fro(Jm)*dq + (g(len(cT)) + 5*u)*fro(np.abs(Jm) @ aQT)*(1 + 1e-6)
            eS = math.hypot(eS, eb)
            S[n + offs[i]:n + offs[i+1]] = bot; del bot
        eS2 += eS**2
        # Phi_m = R^-1 Q^* S_m (minimos quadrados / Tikhonov)
        Phi = solve_triangular(Rf, QfH @ S, check_finite=False)
        fP = fro(Phi)
        Rm = F @ Phi; Rm -= S
        rho2 += fro(Rm)**2
        erho2 += (g(r)*fF*fP + u*(fro(S) + fro(Rm)))**2
        G += Phi @ Phi.conj().T; eG += g(len(cT) + 300)*fP**2
        ncol += len(cT)
        del S, Rm, Phi
        ctx.Q.cache_clear()
        if m % 8 == 0:
            log(f'  modo {m}: rho parcial {math.sqrt(rho2):.4e}')
    ctx._B.cache_clear()
    rho = math.sqrt(rho2)*(1 + 1e-9) + math.sqrt(erho2) + math.sqrt(eS2)
    # ||Phi||^2 = lambda_max(G) (G exata = sum Phi_m Phi_m^*; |G^ - G| <= eG em Frobenius)
    G = (G + G.conj().T)/2                       # Hermitiana: ||G_sim - G_exata|| <= ||G^ - G_exata||
    lamG, V = np.linalg.eigh(G)
    Res = G @ V - V*lamG
    eRes = g(r)*fro(np.abs(G) @ np.abs(V)) + u*fro(Res) + u*fro(V*lamG)
    nRes = fro(Res) + eRes
    D = V.conj().T @ V - np.eye(r); dD = fro(D) + g(r)*fro(V)**2 + u*fro(D)
    assert dD < 0.5
    # Ostrowski: lambda_max(G^) <= lambda_max(V^* G^ V)/sigma_min(V)^2; V^* G^ V = Lambda D' + V^* Res, ...
    # mais simples: lambda_max(V^*G^V) <= max lambda^ (1 + dD) + ||V|| ||Res||  (V^* G^ V = V^*V Lambda + V^* Res)
    nV = math.sqrt(1 + dD)
    lmax = float(np.abs(lamG).max())*(1 + dD) + nV*nRes
    lmaxG = lmax/(1 - dD) + eG
    phi = math.sqrt(lmaxG)*(1 + 1e-12)
    out = dict(c=[c.real, c.imag], sha256_F=sha, tik_rel=a.tik, lam=lam, phi=phi, rho=rho, rho_fl=math.sqrt(rho2),
               erro_rho=math.sqrt(erho2), erro_S=math.sqrt(eS2), qTT_est=qTT, colunas=ncol, posto=r,
               eG=eG, dD=dD, nRes=nRes, lam_max_fl=float(lamG.max()))
    json.dump(out, open(a.saida, 'w'), indent=1)
    log(json.dumps(out))

if __name__ == '__main__':
    main()
