#!/usr/bin/env python3
"""S3b, L, modo rigoroso, ETAPA B: os niveis de cauda (D_J densos, R_J, far) do certificado NK.

Mesmos niveis de nk_L2.py (janelas alinhadas com a borda de Z; densas as fortes em n < N1). O
inverso aproximado nas janelas densas e o EXATO da matriz calculada Hh_J = I + fl(A_JJ); nas
demais, a identidade. Cada entrada de Perron e uma cota VERIFICADA:

  linha densa r:  ||Hh_r^-1 X|| por Loewner (verified_inverse), mais ||Hh_r^-1|| erro(X);
  linha R, far:   cota_norma (Rump no Gram do lado menor), mais erro(X);
  erro(X)       = produto (gamma ||B||_F ||Q||_F) + ||B|| erro(Q0) + (dB_L + eps_mu) ||Q0 P_c||,
                  com dB_L = arredondamento da montagem + resto da banda (|d| > BANDA) + erro de RT;
  diagonal densa: ||Hh_J^-1|| (erro de Hh_J);
  far:            coluna far <- (descidas + perturbacao) nas linhas densas; beta nas R e no far.
Residuo por nivel com h0w (etapa A), banda completa (|d| <= 40), mais a perturbacao global.
Grava T.json para a etapa C (montagem, Perron em racionais, NK)."""
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
import fundo_L
import nk_L as nk
import nk_L2
import signed_operator_L as sl
from nk_L2_Z import Q_cert, err_mm, FATOR_PESO
from tile_perron import DELTA_MU
from verified_inverse import cota_inv, cota_inv_X
from verified_norms import cota_norma, gamma, U

MU, K2 = sl.MU, sl.K2


class Cauda(nk_L2.Certificado):
    def __init__(self, a):
        super().__init__(a)
        self.fb = fundo_L.dB_L(self.ctx.g, nk_L2.BANDA)
        self.emu = fundo_L.eps_mu_L(complex(a.lam0))
        self.nB = nk.NORMA_BL + self.fb['total']          # ||B|| (simbolo em Arb) + perturbacao
        self.pert = self.fb['total'] + self.emu              # perturbacao global (por ||Q0 P_c||)
        self._Q = {}

    def Qc(self, m):
        if m not in self._Q:
            self._Q = {k: v for k, v in self._Q.items() if abs(k - m) <= 48}   # memoria
            self._Q[m] = Q_cert(m, self.Nc, self.a.lam0)
        return self._Q[m]

    def produto_cert(self, lin, col, No):
        """X = fl(P_lin B_BANDA Q0 P_col), erro(X) (sem a perturbacao global) e ||Q0 P_col||."""
        ctx, Nc = self.ctx, self.Nc
        mc = [m for m in col if len(col[m])]; ml = [m for m in lin if len(lin[m])]
        oc = np.cumsum([0] + [len(col[m]) for m in mc]); ol = np.cumsum([0] + [len(lin[m]) for m in ml])
        X = np.zeros((int(ol[-1]), int(oc[-1])), complex)
        e2 = 0.0; dQ = 0.0; nQ = 0.0
        for q, mi in enumerate(mc):
            Q, dq, nq = self.Qc(mi)
            Qi = Q[:, col[mi]]
            dQ = max(dQ, dq); nQ = max(nQ, nq + dq)
            for p, mo in enumerate(ml):
                if abs(mo - mi) <= nk_L2.BANDA:
                    Bb = ctx.B(mo, mi, No, Nc)[lin[mo]]
                    X[ol[p]:ol[p+1], oc[q]:oc[q+1]] = Bb @ Qi
                    e2 += err_mm(Bb, Qi)**2
        erro = math.sqrt(e2) + self.nB*dQ + self.fb['arred']*nQ
        return X, erro, nQ

    def densos_cert(self):
        ctx, Nc = self.ctx, self.Nc
        for lv in self.niveis:
            if lv['tipo'] != 'D':
                continue
            X, erro, nQ = self.produto_cert(self.lin_nivel(lv, Nc), lv['sel'], Nc)
            H = X; H[np.diag_indices(len(H))] += 1.0
            lv['arqH'] = self.a.dir/f"H_{lv['nome'].replace(' ', '_')}.npy"
            np.save(lv['arqH'], H)
            lv['nV'] = cota_inv(H)*FATOR_PESO
            lv['dH'] = erro + U*float(np.linalg.norm(H)) + self.pert*nQ
            del H, X
            self.log(f"  {lv['nome']} ({lv['dim']}): ||Hh^-1|| <= {lv['nV']:.4f}, erro de H {lv['dH']:.1e}")

    def norma_libera(self, caixa, erro):
        """cota_norma com a matriz numa caixa (lista): o Gram e formado e a matriz sai da memoria
        antes do Rump (central: 14,5 mil)."""
        from verified_norms import cota_gram, _potencia
        X = caixa.pop()
        if X.size == 0:
            return erro
        Y = X if X.shape[1] <= X.shape[0] else X.conj().T
        del X
        G = Y.conj().T @ Y
        eG = gamma(Y.shape[0])*float(np.sum(np.abs(Y)**2))
        ch = _potencia(Y)
        del Y
        return cota_gram(G, eG, chute=ch) + erro

    def entrada(self, r, caixa, erro, nQ):
        """Cota verificada de ||(Hh_r^-1) P_r (B Q0) P_c|| para o operador VERDADEIRO."""
        erro_tot = erro + self.pert*nQ
        if r['tipo'] == 'D':
            X = caixa.pop()
            H = np.load(r['arqH'])
            v = cota_inv_X(H, X, erro_X=erro_tot, nHinv=r['nV'])
            del H, X
            return v*FATOR_PESO
        return self.norma_libera(caixa, erro_tot)*FATOR_PESO

    def rodar(self):
        a, ctx, MZ, NZ, Mc, Nc, No = self.a, self.ctx, self.MZ, self.NZ, self.Mc, self.Nc, self.No
        z = json.loads((a.dir/'Z.json').read_text())
        self.log(f"etapa A: ||Mh^-1|| {z['tV']:.4f}, a {z['aZ']:.4f}; dB_L {self.fb['total']:.2e}")
        self.densos_cert()
        niv = self.niveis; k = len(niv); iZ, iF = 0, k + 1
        N = np.zeros((k + 2, k + 2))
        nome = ['Z\''] + [lv['nome'] for lv in niv] + ['far']
        Qf = z['Qfar']
        beta = (nk.NORMA_BL + self.fb['rt'])*Qf + self.emu     # o simbolo ja e o da banda completa
        desc = lambda J: self.nB*(cf.descida(J + nk.BAND_N + 1, Nc, MU, K2) + 6*cf.descida_divxi(J, Nc, K2))
        N[iZ, iF] = z['epsZ']; N[iF, iF] = beta
        N[iF, iZ] = self.pert*z['nQZ']              # so a perturbacao global liga Z ao far
        for i, lv in enumerate(niv, 1):
            N[i, iF] = lv['nV']*(desc(self.N1) + self.pert*Qf) if lv['tipo'] == 'D' else beta
        # colunas de Z'
        colZ = {m: nk.nZ_local(ctx, m, NZ, Nc) for m in nk.sl.modos(MZ)}
        for i, lv in enumerate(niv, 1):
            if self.perto(lv['modos'], list(colZ), nk.BAND_M):
                if self.perto(lv['modos'], list(colZ)):
                    X, erro, nQ = self.produto_cert(self.lin_nivel(lv, Nc), colZ, Nc)
                    cx = [X]; del X
                    N[i, iZ] = self.entrada(lv, cx, erro, nQ)
                else:
                    nQ = max(self.Qc(m)[2] for m in colZ)
                    N[i, iZ] = self.pert*nQ*(lv['nV'] if lv['tipo'] == 'D' else 1.0)
        self.log('T <- Z\': ' + ', '.join(f'{nome[i]} {N[i, iZ]:.4f}' for i in range(1, k + 1) if N[i, iZ] > 5e-4))
        # colunas de cauda
        for j, c in enumerate(niv, 1):
            for i, r in enumerate(niv, 1):
                if i == j and r['tipo'] == 'D':
                    N[i, i] = r['nV']*r['dH']*FATOR_PESO
                    continue
                if self.perto(r['modos'], c['modos']):
                    X, erro, nQ = self.produto_cert(self.lin_nivel(r, No), c['sel'], No)
                    cx = [X]; del X
                    N[i, j] = self.entrada(r, cx, erro, nQ)
                elif self.perto(r['modos'], c['modos'], nk.BAND_M):
                    nQ = max(self.Qc(m)[2] for m in c['modos'] if len(c['sel'][m]))
                    N[i, j] = self.pert*nQ*(r['nV'] if r['tipo'] == 'D' else 1.0)
            Xf = []; ef2 = 0.0; nQ = 0.0
            for mo in range(-(Mc + nk_L2.BANDA - 1), Mc + nk_L2.BANDA):
                if min(abs(mo - mi) for mi in c['modos']) > nk_L2.BANDA:
                    continue
                nn = np.concatenate([ns for _, ns in ctx.rad(mo, No).comps])
                lf = np.where(nn >= Nc)[0] if abs(mo) < Mc else np.arange(len(nn))
                X, erro, nq = self.produto_cert({mo: lf}, c['sel'], No)
                Xf.append(X); ef2 += erro**2; nQ = max(nQ, nq)
            c['nQ'] = nQ                                    # ||Q0 P_c|| (para o termo de 2a ordem do NK)
            Xf = np.concatenate(Xf) if len(Xf) > 1 else Xf[0]
            N[iF, j] = self.norma_libera([Xf], math.sqrt(ef2) + self.pert*nQ)*FATOR_PESO
            del Xf
            self.log(f"  coluna {c['nome']}: " + ', '.join(f'{nome[i]} {N[i, j]:.4f}' for i in range(1, k + 2) if N[i, j] > 5e-4))
            ctx.Q.cache_clear(); ctx._B.cache_clear()
        Y = self.residuo_cert(z)
        out = dict(nome=nome, N=N.tolist(), Y=Y, SZ=[i for i, lv in enumerate(niv, 1)
                                                        if min(abs(m) for m in lv['modos']) < MZ + nk.BAND_M],
                   beta=beta, pert=self.pert, dB_L=self.fb, eps_mu=self.emu,
                   niveis=[dict(nome=lv['nome'], tipo=lv['tipo'], dim=lv['dim'], nV=lv.get('nV'), modos=lv['modos'],
                                nQ=lv['nQ'])
                           for lv in niv])
        (a.dir/'T.json').write_text(json.dumps(out, indent=1) + '\n')
        self.log('etapa B gravada (T.json)')

    def residuo_cert(self, z):
        """||(Hh_r^-1) P_r L(lambda0) h0|| por nivel: B h0 nas linhas de cauda (banda completa),
        com erro de produto, mais a perturbacao global (arredondamento, RT, mu) em cada nivel."""
        ctx, MZ, NZ, Nc = self.ctx, self.MZ, self.NZ, self.Nc
        h0w = np.load(self.a.dir/'h0w.npy')
        ms = sl.modos(MZ); rs = [ctx.rad(m, NZ) for m in ms]
        offs = np.cumsum([0] + [r.dim for r in rs])
        NZb = NZ + nk.BAND_N + 1
        assert NZb <= Nc, 'o residuo de h0 (n < NZ + 101) precisa caber em F'
        res = {}; er = {}
        for mo in range(-(MZ + nk.BAND_M - 1), MZ + nk.BAND_M):
            acc = np.zeros(ctx.rad(mo, NZb).dim, complex); e2 = 0.0
            for j, mi in enumerate(ms):
                if abs(mo - mi) <= nk.BAND_M:
                    Bb = ctx.B(mo, mi, NZb, NZ); h = h0w[offs[j]:offs[j+1]]
                    acc += Bb @ h; e2 += err_mm(Bb, h[:, None])**2
            full = np.zeros(ctx.rad(mo, Nc).dim, complex); full[nk_L2.mapa(ctx, mo, NZb, Nc)] = acc
            res[mo] = full; er[mo] = math.sqrt(e2) + U*float(np.linalg.norm(acc))
        ctx._B.cache_clear()
        fb40 = fundo_L.dB_L(ctx.g, nk.BAND_M)
        pert = (fb40['arred'] + fb40['rt'])*1.0001 + DELTA_MU*(z['J1h'] + 2*cf.norma_divxi(K2))
        Y = []
        for lv in self.niveis:
            partes = [res[m][lv['sel'][m]] if m in res else np.zeros(len(lv['sel'][m]), complex)
                      for m in lv['modos'] if len(lv['sel'][m])]
            x = np.concatenate(partes)
            e = math.sqrt(sum(er[m]**2 for m in lv['modos'] if m in res)) + pert
            if lv['tipo'] == 'D':
                H = np.load(lv['arqH'])
                Y.append(cota_inv_X(H, x[:, None], erro_X=e, nHinv=lv['nV'])*FATOR_PESO)
                del H
            else:
                Y.append((float(np.linalg.norm(x))*(1 + 1e-6) + e)*FATOR_PESO)
        return Y


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--lam0', type=float, default=0.4010247330)
    ap.add_argument('--Z', default='32x128'); ap.add_argument('--F', default='200x400')
    ap.add_argument('--W', type=int, default=16)
    ap.add_argument('--densas', default='16,80')
    ap.add_argument('--N1', type=int, default=200)
    ap.add_argument('--iters', type=int, default=40)
    ap.add_argument('--dir', type=Path, default=Path('build/s3b/rig'))
    a = ap.parse_args()
    Cauda(a).rodar()


if __name__ == '__main__':
    main()
