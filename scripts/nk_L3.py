#!/usr/bin/env python3
"""S3b: certificado NK bordejado de L com BLOCO-JACOBI EM TODAS AS JANELAS (inteiras, n < Nc) e a
cauda como UM nivel com norma l2 ponderada. ESTIMATIVA em ponto flutuante. S3B_ROTA.md, 10-11.

Por que: a norma verdadeira da cauda com bloco-Jacobi nas 8 janelas fortes e 0,39, mas a cota por
normas de blocos da 0,74 (o operador e muito nao normal: norma 0,64, raio espectral 0,15); e o far
(beta = 0,495) acoplado a cada nivel, na norma do maximo, amplifica ~sqrt(niveis). Com todas as
janelas densas, a diagonal de blocos some e a agregacao fica numa cadeia de diagonal nula.

Niveis: Z' (bordejado; dados de nk_L.py), janelas J da cauda (inteiras em F, V_J = (I + A_JJ)^-1;
a central |m| < c0 dividida em duas), far. Norma: max(||x_Z'||/v_0, ||x_T||_W/v_T, ||x_far||/v_f),
||x_T||_W^2 = sum w_J^2 ||x_J||^2. Entradas:
  T<-T  = ||W N W^-1||_2, N_IJ = ||V_I A_IJ|| (I != J), N_JJ = rho_J (~0);
  T<-Z' = sqrt(sum w_J^2 ||V_J P_J B Q0 P_Z||^2);   Z'<-T = a/min w (a conjunta, nk_L.py);
  T<-far = max_J w_J ||V_J P_(J, n >= Nc - BAND_N - g)|| beta + sqrt(sum w_J^2 ||V_J||^2) desc(Nc - BAND_N - g)
           (o far so alcanca as linhas altas de cada janela; nas baixas, so por descida);
  far<-T = ||Nf W^-1||_2 (Nf: grupos de linhas do far x janelas);  far<-far = beta.
As entradas densas saem por potencia com a LU de Hh_J (sem formar V_J A_IJ)."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.linalg import lu_factor, lu_solve

sys.path.insert(0, str(Path(__file__).resolve().parent))
import closed_form as cf
import nk_L as nk
import nk_L2
from nk_L2 import BANDA, mapa, janelas
from tile_certificate import norma2

MU, K2 = nk.MU, nk.K2


def potencia_lu(lu, X, iters=60):
    """||H^-1 X|| por potencia (H = LU), sem formar H^-1 X."""
    rng = np.random.default_rng(0)
    v = rng.standard_normal(X.shape[1]) + 1j*rng.standard_normal(X.shape[1]); v /= np.linalg.norm(v)
    for _ in range(iters):
        u = lu_solve(lu, X @ v)
        v = X.conj().T @ lu_solve(lu, u, trans=2)
        nv = np.linalg.norm(v)
        if nv == 0:
            return 0.0
        v /= nv
    return float(np.linalg.norm(lu_solve(lu, X @ v)))


def janelas_L3(MZ, Mc, W):
    js = janelas(MZ, Mc, W)
    c = js[0]                                   # central |m| < c0: divide em duas
    meio = len(c)//2
    return [c[:meio + 1], c[meio + 1:]] + js[1:]


class Cert3(nk_L2.Certificado):
    def __init__(self, a):
        self.a = a
        self.MZ, self.NZ = map(int, a.Z.split('x')); self.Mc, self.Nc = map(int, a.F.split('x'))
        self.ctx = nk.Contexto(a.lam0, self.Nc)
        self.No = self.Nc + nk.BAND_N + 1
        self.t0 = time.time(); self._mapas = {}
        ctx = self.ctx
        self.niveis = []
        for ms in janelas_L3(self.MZ, self.Mc, a.W):
            sel = {m: nk.colsT(ctx, m, self.MZ, self.NZ, self.Nc) for m in ms}
            self.niveis.append(dict(nome=f'{ms[0]}..{ms[-1]}', tipo='D', modos=ms, sel=sel,
                                    dim=int(sum(len(x) for x in sel.values()))))
        self.log(f'{len(self.niveis)} janelas densas: ' + ', '.join(f"{lv['nome']} ({lv['dim']})" for lv in self.niveis))

    def H(self, lv):
        X = self.produto(self.lin_nivel(lv, self.Nc), lv['sel'], self.Nc)
        X[np.diag_indices(len(X))] += 1.0
        return X

    def rodar(self):
        a, ctx, MZ, NZ, Mc, Nc, No = self.a, self.ctx, self.MZ, self.NZ, self.Mc, self.Nc, self.No
        aZ, nVb, rhob, YZ, rT = map(float, a.dados_Z.split(','))
        niv = self.niveis; k = len(niv)
        N = np.zeros((k, k)); TZ = np.zeros(k); nV = np.zeros(k); nVhi = np.zeros(k)
        colZ = {m: nk.nZ_local(ctx, m, NZ, Nc) for m in nk.sl.modos(MZ)}
        g = a.folga_far
        lim = Nc - nk.BAND_N - g                      # linhas n >= lim: alcancadas pelo far
        ck = Path(str(a.saida) + '.parcial') if a.saida else None
        feitos = json.loads(ck.read_text()) if ck is not None and ck.exists() else {}
        for i, r in enumerate(niv):
            if r['nome'] in feitos:                          # checkpoint (reinicio da maquina)
                f = feitos[r['nome']]
                N[i] = f['N']; TZ[i] = f['TZ']; nV[i] = f['nV']; nVhi[i] = f['nVhi']
                self.log(f"  {r['nome']}: do checkpoint")
                continue
            Hr = self.H(r)
            lu = lu_factor(Hr, overwrite_a=True, check_finite=False); del Hr
            n = r['dim']
            # ||V_r|| e ||V_r P_alto|| (colunas das linhas altas da janela)
            nn = np.concatenate([np.concatenate([ns for _, ns in ctx.rad(m, Nc).comps])[r['sel'][m]]
                                 for m in r['modos'] if len(r['sel'][m])])
            I = np.eye(n)
            nV[i] = potencia_lu(lu, I)
            alto = np.where(nn >= lim)[0]
            nVhi[i] = potencia_lu(lu, I[:, alto]) if len(alto) else 0.0
            del I
            for j, c in enumerate(niv):
                if j == i or not self.perto(r['modos'], c['modos']):
                    continue
                X = self.produto(self.lin_nivel(r, Nc), c['sel'], Nc)
                N[i, j] = potencia_lu(lu, X); del X
            if self.perto(r['modos'], list(colZ)):
                X = self.produto(self.lin_nivel(r, Nc), colZ, Nc)
                TZ[i] = potencia_lu(lu, X); del X
            del lu
            self.log(f"  {r['nome']} ({n}): ||V|| {nV[i]:.3f}, ||V P_alto|| {nVhi[i]:.3f}; "
                     + ', '.join(f"<- {niv[j]['nome']} {N[i, j]:.3f}" for j in range(k) if N[i, j] > 5e-4)
                     + (f'; <- Z\' {TZ[i]:.3f}' if TZ[i] else ''))
            ctx.Q.cache_clear(); ctx._B.cache_clear()
            if ck is not None:
                feitos[r['nome']] = dict(N=N[i].tolist(), TZ=float(TZ[i]), nV=float(nV[i]), nVhi=float(nVhi[i]))
                ck.write_text(json.dumps(feitos))
        # far: linhas por grupo x janelas (sem V)
        grupos = janelas(MZ, Mc, a.W) + [list(range(Mc, Mc + BANDA)), list(range(-(Mc + BANDA - 1), -Mc + 1))]
        Nf = np.zeros((len(grupos), k))
        for j, c in enumerate(niv):
            for gi, w in enumerate(grupos):
                Xg = []
                for mo in w:
                    if min(abs(mo - mi) for mi in c['modos']) > BANDA:
                        continue
                    nn = np.concatenate([ns for _, ns in ctx.rad(mo, No).comps])
                    lf = np.where(nn >= Nc)[0] if abs(mo) < Mc else np.arange(len(nn))
                    Xg.append(self.produto({mo: lf}, c['sel'], No))
                if Xg:
                    Nf[gi, j] = norma2(np.concatenate(Xg))
            ctx.Q.cache_clear(); ctx._B.cache_clear()
        self.log('far <- janelas: ' + ', '.join(f"{c['nome']} {np.linalg.norm(Nf[:, j]):.3f}" for j, c in enumerate(niv)))
        # far
        rad = max(max(cf.cauda_Rinv(a.lam0 + 1j*m/2, Nc, p, MU, K2, 1500) for p in ((2, 1, 1) if m % 2 == 0 else (1,)))
                  for m in range(-(Mc - 1), Mc))
        Qf = max(rad, cf.Rinv_uniforme(Mc, MU, a.lam0, K2)); beta = nk.NORMA_BL*Qf
        desc = nk.NORMA_BL*(cf.descida(lim + nk.BAND_N + 1, Nc, MU, K2) + 6*cf.descida_divxi(lim, Nc, K2))
        epsZ = nVb*nk.NORMA_BL*(cf.descida(NZ + nk.BAND_N + 1, Nc, MU, K2) + 6*cf.descida_divxi(NZ, Nc, K2))
        SZ = [i for i, lv in enumerate(niv) if min(abs(m) for m in lv['modos']) < MZ + nk.BAND_M]
        # pesos: 1 e sqrt(u/z) (vetores de Perron de N)
        ws = {'1': np.ones(k)}
        lam, Vr = np.linalg.eig(N); z = np.abs(Vr[:, np.argmax(np.abs(lam))])
        lam, Vl = np.linalg.eig(N.T); u = np.abs(Vl[:, np.argmax(np.abs(lam))])
        w = np.sqrt(np.maximum(u, 1e-9)/np.maximum(z, 1e-9)); ws['sqrt(u/z)'] = w/w.max()
        res = {}
        for nomew, w in ws.items():
            nT = float(np.linalg.norm((w[:, None]*N)/w[None, :], 2))
            TZw = math.sqrt(float(np.sum((w*TZ)**2)))
            ZT = aZ/min(w[j] for j in SZ)
            Tf = max(w*nVhi)*beta + math.sqrt(float(np.sum((w*nV)**2)))*desc
            fT = float(np.linalg.norm(Nf/w[None, :], 2))
            M = np.array([[rhob, ZT, epsZ], [TZw, nT, Tf], [0.0, fT, beta]])
            th = float(max(abs(np.linalg.eigvals(M))))
            lamM, vec = np.linalg.eig(M); v = np.abs(vec[:, np.argmax(np.abs(lamM))]); v /= v.max()
            Y = max(YZ/v[0], rT*max(nV)*math.sqrt(len(SZ))/v[1])
            self.log(f'pesos {nomew}: T<-T {nT:.4f} (raio {max(abs(np.linalg.eigvals(N))):.4f}), T<-Z\' {TZw:.3f}, '
                     f'Z\'<-T {ZT:.3f}, T<-far {Tf:.3f}, far<-T {fT:.3f}, beta {beta:.3f}:  THETA0 = {th:.4f}; '
                     f'raio NK ~ {Y/(1 - th) if th < 1 else float("inf"):.1e}')
            res[nomew] = dict(M=M.tolist(), theta0=th, nT=nT, w=w.tolist())
        out = dict(nome=[lv['nome'] for lv in niv], N=N.tolist(), TZ=TZ.tolist(), nV=nV.tolist(), nVhi=nVhi.tolist(),
                   Nf=Nf.tolist(), beta=beta, desc=desc, epsZ=epsZ, aZ=aZ, SZ=SZ, res=res)
        if a.saida:
            a.saida.write_text(json.dumps(out, indent=1) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--lam0', type=complex, default=0.4010247330)
    ap.add_argument('--Z', default='32x128'); ap.add_argument('--F', default='200x400')
    ap.add_argument('--W', type=int, default=16)
    ap.add_argument('--folga-far', type=int, default=60, help='g: o far alcanca as linhas n >= Nc - BAND_N - g')
    ap.add_argument('--dados-Z', required=True, help='a,||V_b||,rho,Y_Z,r_T de nk_L.py')
    ap.add_argument('--saida', type=Path)
    a = ap.parse_args()
    Cert3(a).rodar()


if __name__ == '__main__':
    main()
