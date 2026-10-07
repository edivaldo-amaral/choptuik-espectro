#!/usr/bin/env python3
"""S3b: certificado NK BORDEJADO de L com BLOCO-JACOBI nas janelas fortes da cauda.
ESTIMATIVA em ponto flutuante. Ver docs/S3B_ROTA.md, secao 9.

x = (y, lambda), h = Q0 y, F(x) = (L(lambda) h, l(h) - 1); DF(x0) = [[H, h0], [l Q0, 0]],
H = I + B Q0 (norma do toro de L, (65/64, 5/4), m com sinal). Niveis:
  Z'   = Z + {lambda}: inverso bordejado denso V_b (dados de uma rodada de nk_L.py: --dados-Z);
  D_J  = janela forte J (|m| em [--densas]) restrita a n < N1: V_J = (I + A_JJ)^-1 denso;
  R_J  = o resto de F = Mc x Nc dentro da janela J: identidade;
  far  = complemento de F: ||B_L|| (simbolo) vezes ||Q0 P_far|| (forma fechada).
A = diag(V_b, V_J, I, I), E = I - A DF. Janelas de W modos alinhadas com a borda de Z:
|m| < c0 (central) e [c0 + kW, c0 + (k+1)W) dos dois lados (c0 = MZ mod W, ou W).

O acoplamento A_JJ da cauda tem norma ~0,64 e raio espectral ~0,15 (muito nao normal): com V_J
a diagonal de Perron zera nas janelas fortes e sobram acoplamentos ~0,25. A linha de Z' usa a
norma CONJUNTA a = ||V_b [P_Z B Q0 P_T; l Q0 P_T]|| com agregacao l2: sua entrada e
a sqrt(sum_{j perto de Z} v_j^2). Banda exata |d| <= BANDA; o resto do fundo (|d| > BANDA) entra
em todas as entradas como resto ||Q0 P_c|| (vezes ||V_r|| nas linhas densas)."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import closed_form as cf
import nk_L as nk
from tile_certificate import norma2

BANDA = 16
MU, K2 = nk.MU, nk.K2


def janelas(MZ, Mc, W):
    c0 = MZ % W or W
    pos = [list(range(k, min(k + W, Mc))) for k in range(c0, Mc, W)]
    return [list(range(-(c0 - 1), c0))] + pos + [[-m for m in reversed(j)] for j in pos]


def mapa(ctx, m, Na, Nb):
    """Posicoes, no layout radial (m, Nb), dos indices do layout (m, Na) (Na <= Nb)."""
    ra, rb = ctx.rad(m, Na), ctx.rad(m, Nb)
    pos = {}; o = 0
    for c, ns in rb.comps:
        for k, n in enumerate(ns):
            pos[(c, int(n))] = o + k
        o += len(ns)
    return np.array([pos[(c, int(n))] for c, ns in ra.comps for n in ns], int)


def potencia(mv, rmv, n, iters):
    if n == 0:
        return 0.0
    return nk._potencia(mv, rmv, n, iters=iters)


class Certificado:
    def __init__(self, a):
        self.a = a
        self.MZ, self.NZ = map(int, a.Z.split('x')); self.Mc, self.Nc = map(int, a.F.split('x'))
        self.N1 = a.N1
        self.ctx = nk.Contexto(a.lam0, self.Nc)
        self.No = self.Nc + nk.BAND_N + 1
        self.t0 = time.time()
        self._mapas = {}
        self.dirV = Path(getattr(a, 'dirV', 'build/s3b/tmpV')); self.dirV.mkdir(parents=True, exist_ok=True)
        d_lo, d_hi = map(int, a.densas.split(','))
        ctx, MZ, NZ, Nc = self.ctx, self.MZ, self.NZ, self.Nc
        self.niveis = []
        for J, ms in enumerate(janelas(MZ, self.Mc, a.W)):
            forte = min(abs(m) for m in ms) >= d_lo and max(abs(m) for m in ms) < d_hi
            selD, selR = {}, {}
            for m in ms:
                t = nk.colsT(ctx, m, MZ, NZ, Nc)
                n = np.concatenate([ns for _, ns in ctx.rad(m, Nc).comps])[t]
                selD[m] = t[n < self.N1] if forte else t[:0]
                selR[m] = t[n >= self.N1] if forte else t
            nome = f'{ms[0]}..{ms[-1]}'
            if forte:
                self.niveis.append(dict(nome='D ' + nome, tipo='D', modos=ms, sel=selD))
            if sum(len(x) for x in selR.values()):
                self.niveis.append(dict(nome='R ' + nome, tipo='R', modos=ms, sel=selR))
        for lv in self.niveis:
            lv['dim'] = int(sum(len(x) for x in lv['sel'].values()))
        self.log(f'{len(self.niveis)} niveis de cauda: '
                 + ', '.join(f"{lv['nome']} ({lv['dim']})" for lv in self.niveis))

    def log(self, s):
        print(f'{s}   [{time.time()-self.t0:.0f}s]', flush=True)

    def perto(self, ms1, ms2, b=BANDA):
        return min(abs(x - y) for x in ms1 for y in ms2) <= b

    def produto(self, lin, col, No):
        """Bloco denso P_lin B Q0 P_col (banda |d| <= BANDA); lin = {m: linhas no layout (m, No)}.
        So as linhas pedidas (memoria: a janela inteira de linhas chegaria a ~6 GB)."""
        ctx, Nc = self.ctx, self.Nc
        mc = [m for m in col if len(col[m])]; ml = [m for m in lin if len(lin[m])]
        oc = np.cumsum([0] + [len(col[m]) for m in mc]); ol = np.cumsum([0] + [len(lin[m]) for m in ml])
        X = np.zeros((int(ol[-1]), int(oc[-1])), complex)
        for q, mi in enumerate(mc):
            Qi = ctx.Q(mi, Nc)[:, col[mi]]
            for p, mo in enumerate(ml):
                if abs(mo - mi) <= BANDA:
                    X[ol[p]:ol[p+1], oc[q]:oc[q+1]] = ctx.B(mo, mi, No, Nc)[lin[mo]] @ Qi
        return X

    def lin_nivel(self, lv, No):
        return {m: self.mapas(m, No)[lv['sel'][m]] for m in lv['modos'] if len(lv['sel'][m])}

    def mapas(self, m, No):
        k = (m, No)
        if k not in self._mapas:
            self._mapas[k] = mapa(self.ctx, m, self.Nc, No)
        return self._mapas[k]

    def densos(self):
        """V_J = (I + A_JJ)^-1 nas janelas fortes; rho_J = ||I - V_J H_JJ||."""
        ctx, Nc = self.ctx, self.Nc
        for lv in self.niveis:
            if lv['tipo'] != 'D':
                continue
            ms = [m for m in lv['modos'] if len(lv['sel'][m])]
            o = np.cumsum([0] + [len(lv['sel'][m]) for m in ms]); n = int(o[-1])
            H = np.eye(n, dtype=complex)
            for q, mi in enumerate(ms):
                Qi = ctx.Q(mi, Nc)[:, lv['sel'][mi]]
                for p, mo in enumerate(ms):
                    if abs(mo - mi) <= BANDA:
                        H[o[p]:o[p+1], o[q]:o[q+1]] += ctx.B(mo, mi, Nc, Nc)[lv['sel'][mo]] @ Qi
            nA = norma2(H - np.eye(n))
            V = np.linalg.inv(H)
            X = V @ H; X[np.diag_indices(n)] -= 1
            lv['nV'] = norma2(V); lv['rho'] = norma2(X); lv['nA'] = nA
            lv['arqV'] = self.dirV/f"V_{lv['nome'].replace(' ', '_')}.npy"
            np.save(lv['arqV'], V)
            del H, X, V
            self.log(f"  {lv['nome']} ({n}): ||A_JJ|| {nA:.4f}, ||V_J|| {lv['nV']:.4f}, rho {lv['rho']:.1e}")

    def resto_banda(self):
        """sum_{BANDA < |d| <= BAND_M} max_paridade ||B_d|| (pesado); o de |d| > BAND_M ja e
        desprezivel (~1e-9)."""
        tot = 0.0
        for d in range(-nk.BAND_M, nk.BAND_M + 1):
            if abs(d) > BANDA:
                tot += max(norma2(self.ctx._B(d, par, self.No, self.Nc)) for par in (0, 1))
        self.ctx._B.cache_clear()
        return tot

    def q0(self, lv):
        return max(norma2(self.ctx.Q(m, self.Nc)[:, lv['sel'][m]]) for m in lv['modos'] if len(lv['sel'][m]))

    def rodar(self):
        a, ctx, MZ, NZ, Mc, Nc, No = self.a, self.ctx, self.MZ, self.NZ, self.Mc, self.Nc, self.No
        aZ, nVb, rhob, YZ, rT = map(float, a.dados_Z.split(','))
        self.log(f'nivel Z (rodada anterior): a {aZ}, ||V_b|| {nVb}, rho {rhob}, Y_Z {YZ}')
        self.densos()
        resto = self.resto_banda()
        self.log(f'resto da banda (|d| > {BANDA}): {resto:.2e}')
        niv = self.niveis; k = len(niv); iZ, iF = 0, k + 1
        N = np.zeros((k + 2, k + 2))
        nome = ['Z\''] + [lv['nome'] for lv in niv] + ['far']
        # far
        rad = 0.0
        for m in range(-(Mc - 1), Mc):
            sig = a.lam0 + 1j*m/2
            rad = max(rad, max(cf.cauda_Rinv(sig, Nc, p, MU, K2, 1500) for p in ((2, 1, 1) if m % 2 == 0 else (1,))))
        Qf = max(rad, cf.Rinv_uniforme(Mc, MU, a.lam0, K2))
        beta = nk.NORMA_BL*Qf
        desc = lambda J: nk.NORMA_BL*(cf.descida(J + nk.BAND_N + 1, Nc, MU, K2) + 6*cf.descida_divxi(J, Nc, K2))
        N[iZ, iF] = nVb*desc(NZ)
        N[iF, iF] = beta
        for i, lv in enumerate(niv, 1):
            N[i, iF] = lv['nV']*desc(self.N1) if lv['tipo'] == 'D' else beta
        self.log(f'far: ||Q0 P_far|| {Qf:.4f}, beta {beta:.4f}; Z\' <- far {N[iZ, iF]:.1e}')
        def norma_linha(i, X):
            lv = niv[i - 1]
            if lv['tipo'] != 'D':
                return norma2(X)
            V = np.load(lv['arqV']); v = norma2(V @ X); del V
            return v
        # colunas de Z': so linhas de cauda (o far nao alcanca Z com o fundo truncado)
        colZ = {m: nk.nZ_local(ctx, m, NZ, Nc) for m in nk.sl.modos(MZ)}
        q0Z = max(norma2(ctx.Q(m, Nc)[:, colZ[m]]) for m in colZ)
        for i, lv in enumerate(niv, 1):
            if self.perto(lv['modos'], list(colZ)):
                X = self.produto(self.lin_nivel(lv, Nc), colZ, Nc)
                N[i, iZ] = norma_linha(i, X); del X
        for i, lv in enumerate(niv, 1):
            if self.perto(lv['modos'], list(colZ), nk.BAND_M):
                N[i, iZ] += resto*q0Z*(lv['nV'] if lv['tipo'] == 'D' else 1.0)
        self.log('T <- Z\': ' + ', '.join(f'{nome[i]} {N[i, iZ]:.3f}' for i in range(1, k + 1) if N[i, iZ] > 5e-4))
        # grupos de linhas do far: as janelas de F e o far de Fourier dos dois lados
        grupos = janelas(MZ, Mc, a.W) + [list(range(Mc, Mc + BANDA)), list(range(-(Mc + BANDA - 1), -Mc + 1))]
        Nf = np.zeros((len(grupos), k + 2))
        # colunas de cauda
        for j, c in enumerate(niv, 1):
            q0 = self.q0(c)
            for i, r in enumerate(niv, 1):
                if (i == j and r['tipo'] == 'D') or not self.perto(r['modos'], c['modos']):
                    continue
                X = self.produto(self.lin_nivel(r, No), c['sel'], No)
                N[i, j] = norma_linha(i, X); del X
            # linhas de far (n >= Nc nos modos de F, todas nos de fora), por grupo de linhas
            # (janelas de F e far de Fourier): Nf[g, j] para a agregacao l2; N[far, j] conjunta
            Xall = []
            for g, w in enumerate(grupos):
                Xg = []
                for mo in w:
                    if min(abs(mo - mi) for mi in c['modos']) > BANDA:
                        continue
                    nn = np.concatenate([ns for _, ns in ctx.rad(mo, No).comps])
                    lf = np.where(nn >= Nc)[0] if abs(mo) < Mc else np.arange(len(nn))
                    Xg.append(self.produto({mo: lf}, c['sel'], No))
                if Xg:
                    Xg = np.concatenate(Xg); Nf[g, j] = norma2(Xg); Xall.append(Xg)
            if Xall:
                Xall = np.concatenate(Xall); N[iF, j] = norma2(Xall)
            del Xall
            for i, r in enumerate(niv, 1):
                if self.perto(r['modos'], c['modos'], nk.BAND_M):
                    N[i, j] += resto*q0*(r['nV'] if r['tipo'] == 'D' else 1.0)
            N[iF, j] += resto*q0
            self.log(f"  coluna {c['nome']} (||Q0 P|| {q0:.3f}): " + ', '.join(
                f'{nome[i]} {N[i, j]:.3f}' for i in range(1, k + 2) if N[i, j] > 5e-4))
            ctx.Q.cache_clear(); ctx._B.cache_clear()
        for i, lv in enumerate(niv, 1):
            if lv['tipo'] == 'D':
                N[i, i] += lv['rho']          # ja contem o resto da banda (V_J inverte so |d| <= BANDA)
        N[iZ, iZ] = rhob
        # linha de Z' (l2): niveis alcancados por a (colunas com |m| < MZ + BAND_M)
        SZ = [i for i, lv in enumerate(niv, 1) if min(abs(m) for m in lv['modos']) < MZ + nk.BAND_M]
        # residuo por nivel, CONSERVADOR: a cota conjunta r_T em cada nivel que o residuo alcanca
        # (|m| < MZ + BAND_M, n < NZ + BAND_N + 1), vezes ||V_J|| nos densos
        NZb = NZ + nk.BAND_N + 1
        Y = []
        for lv in niv:
            alcanca = any(abs(m) < MZ + nk.BAND_M and len(lv['sel'][m]) and
                          np.concatenate([ns for _, ns in ctx.rad(m, Nc).comps])[lv['sel'][m]].min() < NZb
                          for m in lv['modos'])
            Y.append(rT*(lv['nV'] if lv['tipo'] == 'D' else 1.0) if alcanca else 0.0)
        self.Nf = Nf
        return self.perron(N, aZ, SZ, Y, YZ, nome)

    def tres_niveis(self, N, aZ, SZ, nome):
        k = len(self.niveis); T = list(range(1, k + 1))
        D = [lv['tipo'] == 'D' for lv in self.niveis]
        NTT = N[np.ix_(T, T)]
        ws = {'1': np.ones(k)}
        lam, Vr = np.linalg.eig(NTT); z = np.abs(Vr[:, np.argmax(np.abs(lam))])
        lam, Vl = np.linalg.eig(NTT.T); u = np.abs(Vl[:, np.argmax(np.abs(lam))])
        w = np.sqrt(np.maximum(u, 1e-9)/np.maximum(z, 1e-9)); ws['sqrt(u/z)'] = w/w.max()
        melhor = None
        for nomew, w in ws.items():
            th, M, nT = _tres(N, self.Nf, aZ, SZ, w, D)
            self.log(f'3 niveis (T em l2, pesos {nomew}): T<-T {nT:.4f}, theta0 {th:.4f}')
            if melhor is None or th < melhor[0]:
                melhor = (th, M, nomew)
        np.set_printoptions(precision=4, suppress=True)
        print(f'  M3 ({melhor[2]}) =\n{melhor[1]}')
        return dict(theta0_3niveis=melhor[0], M3=melhor[1].tolist(), raio_TT=float(max(abs(np.linalg.eigvals(NTT)))))

    def perron(self, N, aZ, SZ, Y, YZ, nome):
        k = len(self.niveis); iF = k + 1
        def Fv(v):
            w = N @ v
            w[0] = N[0, 0]*v[0] + aZ*math.sqrt(sum(v[j]**2 for j in SZ)) + N[0, iF]*v[iF]
            return w
        v = np.ones(k + 2)
        for _ in range(3000):
            w = Fv(v) + 1e-12
            v = w/w.max()
        th = float(max(Fv(v)/v))
        Nl = N.copy(); Nl[0, SZ] = aZ
        thl = float(max(abs(np.linalg.eigvals(Nl))))
        np.set_printoptions(precision=3, suppress=True, linewidth=200)
        self.log(f'theta0 (linha de Z\' com l2) = {th:.4f};  com a em cada nivel (linear) {thl:.4f}')
        print('  v = ' + ', '.join(f'{nome[i]} {v[i]:.3f}' for i in range(k + 2)))
        out = dict(nome=nome, N=N.tolist(), aZ=aZ, SZ=SZ, v=v.tolist(), theta0=th, theta0_linear=thl,
                   Nf=self.Nf.tolist())
        out.update(self.tres_niveis(N, aZ, SZ, nome))
        if Y is not None:
            Yv = [YZ] + Y + [0.0]
            Ymax = max(Yv[i]/v[i] for i in range(k + 2))
            r = Ymax/(1 - th) if th < 1 else float('inf')
            self.log(f'residuo por nivel: ' + ', '.join(f'{nome[i]} {Yv[i]:.1e}' for i in range(k + 2) if Yv[i]))
            self.log(f'Y = {Ymax:.2e};  raio NK ~ Y/(1 - theta0) = {r:.2e}')
            out.update(Y=Yv, Ymax=Ymax, raio=r)
        return out


def _tres(N, Nf, aZ, SZ, w, D):
    """theta0 de 3 niveis (Z', T, far) com T na norma l2 ponderada ||x||_W^2 = sum w_J^2 ||x_J||^2."""
    k = N.shape[0] - 2; T = list(range(1, k + 1)); iF = k + 1
    NTT = N[np.ix_(T, T)]
    nT = float(np.linalg.norm((w[:, None]*NTT)/w[None, :], 2))
    TZ = math.sqrt(sum((w[t]*N[1 + t, 0])**2 for t in range(k)))
    ZT = aZ/min(w[j - 1] for j in SZ)
    wR = max([w[t] for t in range(k) if not D[t]] or [0.0])
    Tf = math.sqrt((wR*N[iF, iF])**2 + sum((w[t]*N[1 + t, iF])**2 for t in range(k) if D[t]))
    fT = float(np.linalg.norm(Nf[:, 1:k + 1]/w[None, :], 2))
    M = np.array([[N[0, 0], ZT, N[0, iF]], [TZ, nT, Tf], [N[iF, 0], fT, N[iF, iF]]])
    return float(max(abs(np.linalg.eigvals(M)))), M, nT


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--lam0', type=float, default=0.4010247330)
    ap.add_argument('--Z', default='32x128'); ap.add_argument('--F', default='200x400')
    ap.add_argument('--W', type=int, default=16)
    ap.add_argument('--densas', default='16,80', help='|m| das janelas fortes: lo,hi')
    ap.add_argument('--N1', type=int, default=200, help='as janelas fortes sao densas em n < N1')
    ap.add_argument('--iters', type=int, default=40)
    ap.add_argument('--dados-Z', required=True, help='a,||V_b||,rho,Y_Z,r_T de uma rodada de nk_L.py')
    ap.add_argument('--saida', type=Path)
    ap.add_argument('--dirV', default='build/s3b/tmpV', help='onde gravar os V_J densos')
    a = ap.parse_args()
    cert = Certificado(a)
    out = cert.rodar()
    if a.saida:
        out.update(vars(a) | dict(saida=str(a.saida)))
        a.saida.write_text(json.dumps(out, indent=1, default=str) + '\n')


if __name__ == '__main__':
    main()
