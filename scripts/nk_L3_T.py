#!/usr/bin/env python3
"""S3b, L, modo rigoroso, ETAPA B da arquitetura nk_L3: cauda com bloco-Jacobi em TODAS as janelas.

Por janela de linha r (V_r = Hh_r^-1, o inverso EXATO da matriz calculada Hh_r = I + fl(A_rr)):
  P_r = Hh_r Hh_r^* uma vez; cada cota por Loewner (verified_inverse.cota_loewner):
    ||V_r||              (S = I),
    ||V_r P_alto||       (S = projecao nas linhas n >= Nc - BAND_N - g),
    ||V_r A_rc||         (S = X X^*, X = fl(P_r B Q0 P_c), c vizinha) + ||V_r|| (erro(X) + pert ||Q0 P_c||),
    ||V_r P_r B Q0 P_Z|| (idem, colunas de Z),
    rho_r = ||V_r|| erro(Hh_r)   (diagonal de E);
  far: Nf[g, c] = cota_norma das linhas do far do grupo g <- janela c, mais a perturbacao.
Residuo por janela: ||V_r P_r L(lambda0) h0|| (h0w da etapa A), Loewner com S = x x^*.
Grava T3.json para a montagem (nk_L3_C.py)."""
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
import nk_L3
from nk_L2 import BANDA
from nk_L2_T import Cauda
from nk_L2_Z import FATOR_PESO
from tile_perron import DELTA_MU
from verified_inverse import cota_loewner, gram_linhas
from verified_norms import U, gamma
from nk_L2_Z import err_mm

MU, K2 = nk.MU, nk.K2


class Cauda3(Cauda):
    """Niveis de nk_L3 (janelas inteiras, central dividida); maquinaria de cotas de nk_L2_T."""
    def __init__(self, a):
        base = nk_L3.Cert3(a)
        self.__dict__.update(base.__dict__)
        self.fb = fundo_L.dB_L(self.ctx.g, BANDA)
        self.emu = fundo_L.eps_mu_L(complex(a.lam0))
        self.nB = nk.NORMA_BL + self.fb['total']
        self.pert = self.fb['total'] + self.emu
        self._Q = {}

    def loewner(self, P, eP, caixa, erro, nV, chute):
        """||H^-1 X_exata|| <= t + ||H^-1|| erro, t por Loewner com S = X X^* (X numa caixa: sai da
        memoria antes da verificacao)."""
        X = caixa.pop()
        S, eS = gram_linhas(X); del X
        t = cota_loewner(P, eP, S, eS, max(chute, 1e-12)*(1 + 1e-4))
        return (t + nV*erro)*FATOR_PESO

    def rodar(self):
        a, ctx, MZ, NZ, Mc, Nc, No = self.a, self.ctx, self.MZ, self.NZ, self.Mc, self.Nc, self.No
        if getattr(a, 'rouche', False):
            z = self.dados_rouche()                                # sem etapa A: so ||Q_ZZ|| e Qfar em lambda0
        else:
            z = json.loads((a.dir/'Z.json').read_text())
        est = json.loads(Path(a.estimativa).read_text())          # chutes (ponto flutuante, nk_L3)
        niv = self.niveis; k = len(niv)
        N = np.zeros((k, k)); TZ = np.zeros(k); nV = np.zeros(k); nVhi = np.zeros(k); rho = np.zeros(k)
        K = np.zeros(k); nQm = np.zeros(k)
        colZ = {m: nk.nZ_local(ctx, m, NZ, Nc) for m in nk.sl.modos(MZ)}
        lim = Nc - nk.BAND_N - a.folga_far
        if getattr(a, 'rouche', False):
            Yr = [None]*k                                        # Rouche: nao ha residuo
        else:
            h0w = np.load(a.dir/'h0w.npy')
            Yr = self.residuo_janelas(z, h0w)
        Y = np.zeros(k)
        ck = Path(a.ckpt) if a.ckpt else a.dir/'T3.parcial.json'
        feitos = {}
        for arq in [ck] + [Path(x) for x in (a.juntar or [])]:
            if arq.exists():
                feitos.update(json.loads(arq.read_text()))
        ordem = list(range(k))[::-1] if a.reverso else list(range(k))
        for i in ordem:
            r = niv[i]
            if r['nome'] in feitos:                          # checkpoint (reinicio da maquina)
                f = feitos[r['nome']]
                N[i] = f['N']; TZ[i] = f['TZ']; nV[i] = f['nV']; nVhi[i] = f['nVhi']; rho[i] = f['rho']; Y[i] = f['Y']
                K[i] = f.get('K', 0.0); nQm[i] = f.get('nQm', 0.0)
                self.log(f"  {r['nome']}: do checkpoint")
                continue
            X, erro, nQ = self.produto_cert(self.lin_nivel(r, Nc), r['sel'], Nc)
            H = X; H[np.diag_indices(len(H))] += 1.0; del X
            dH = erro + U*float(np.linalg.norm(H)) + self.pert*nQ
            P, eP = gram_linhas(H); del H
            n = r['dim']
            nV[i] = cota_loewner(P, eP, None, 0.0, est['nV'][i]*(1 + 1e-4))*FATOR_PESO
            nn = np.concatenate([np.concatenate([ns for _, ns in ctx.rad(m, Nc).comps])[r['sel'][m]]
                                 for m in r['modos'] if len(r['sel'][m])])
            proj = np.zeros((n, n)); alto = np.where(nn >= lim)[0]; proj[alto, alto] = 1.0
            nVhi[i] = cota_loewner(P, eP, proj, 0.0, max(est['nVhi'][i], 1e-6)*(1 + 1e-4))*FATOR_PESO if len(alto) else 0.0
            del proj
            if getattr(a, 'rouche', False):                  # sensibilidade a s: K_J = ||V_J P_J B Q0(c)^2 P_J||
                X2, erro2, nQm[i] = self.produto_cert_Q2(self.lin_nivel(r, Nc), r['sel'], Nc)
                ch = self.chute_VX(P, X2)
                K[i] = self.loewner(P, eP, [X2], erro2, nV[i], ch)
            rho[i] = nV[i]*dH*FATOR_PESO
            for j, c in enumerate(niv):
                if j == i:
                    continue
                if self.perto(r['modos'], c['modos']):
                    X, erro, nQ = self.produto_cert(self.lin_nivel(r, Nc), c['sel'], Nc)
                    cx = [X]; del X
                    N[i, j] = self.loewner(P, eP, cx, erro + self.pert*nQ, nV[i], est['N'][i][j])
                elif self.perto(r['modos'], c['modos'], nk.BAND_M):
                    nQ = max(self.Qc(m)[2] for m in c['modos'] if len(c['sel'][m]))
                    N[i, j] = nV[i]*self.pert*nQ*FATOR_PESO
            if self.perto(r['modos'], list(colZ), nk.BAND_M):
                if self.perto(r['modos'], list(colZ)):
                    X, erro, nQ = self.produto_cert(self.lin_nivel(r, Nc), colZ, Nc)
                    cx = [X]; del X
                    TZ[i] = self.loewner(P, eP, cx, erro + self.pert*nQ, nV[i], est['TZ'][i])
                else:
                    TZ[i] = nV[i]*self.pert*z['nQZ']*FATOR_PESO
            if Yr[i] is not None:
                x, e = Yr[i]
                Y[i] = self.loewner(P, eP, [x[:, None]], e, nV[i], max(float(np.linalg.norm(x))*nV[i], 1e-14))
            del P
            self.log(f"  {r['nome']}: ||V|| <= {nV[i]:.4f}, ||V P_alto|| <= {nVhi[i]:.4f}, rho {rho[i]:.1e}; "
                     + ', '.join(f"<- {niv[j]['nome']} {N[i, j]:.4f}" for j in range(k) if N[i, j] > 5e-4)
                     + (f"; <- Z' {TZ[i]:.4f}" if TZ[i] else '') + f'; Y {Y[i]:.1e}')
            ctx.Q.cache_clear(); ctx._B.cache_clear()
            feitos[r['nome']] = dict(N=N[i].tolist(), TZ=float(TZ[i]), nV=float(nV[i]), nVhi=float(nVhi[i]),
                                     rho=float(rho[i]), Y=float(Y[i]), K=float(K[i]), nQm=float(nQm[i]))
            ck.write_text(json.dumps(feitos))
        if a.so_janelas:
            self.log('so janelas: checkpoint gravado; a montagem fica para a passada de juncao')
            return
        np.fill_diagonal(N, rho)
        grupos = nk_L2.janelas(MZ, Mc, a.W) + [list(range(Mc, Mc + BANDA)), list(range(-(Mc + BANDA - 1), -Mc + 1))]
        Nf = np.zeros((len(grupos), k))
        for j, c in enumerate(niv):
            for gi, w in enumerate(grupos):
                Xg = []; e2 = 0.0; nQ = 0.0
                for mo in w:
                    if min(abs(mo - mi) for mi in c['modos']) > BANDA:
                        continue
                    nn = np.concatenate([ns for _, ns in ctx.rad(mo, No).comps])
                    lf = np.where(nn >= Nc)[0] if abs(mo) < Mc else np.arange(len(nn))
                    X, erro, nq = self.produto_cert({mo: lf}, c['sel'], No)
                    Xg.append(X); e2 += erro**2; nQ = max(nQ, nq)
                if Xg:
                    Xg = np.concatenate(Xg) if len(Xg) > 1 else Xg[0]
                    Nf[gi, j] = self.norma_libera([Xg], math.sqrt(e2) + self.pert*nQ)*FATOR_PESO
                elif min(min(abs(mo - mi) for mi in c['modos']) for mo in w) <= nk.BAND_M:
                    Nf[gi, j] = self.pert*max(self.Qc(m)[2] for m in c['modos'] if len(c['sel'][m]))
            ctx.Q.cache_clear(); ctx._B.cache_clear()
        self.log('far <- janelas: ' + ', '.join(f"{c['nome']} {np.linalg.norm(Nf[:, j]):.4f}" for j, c in enumerate(niv)))
        Qf = z['Qfar']
        beta = (nk.NORMA_BL + self.fb['rt'])*Qf + self.emu
        desc = self.nB*(cf.descida(lim + nk.BAND_N + 1, Nc, MU, K2) + 6*cf.descida_divxi(lim, Nc, K2)) + self.pert*Qf
        out = dict(nome=[lv['nome'] for lv in niv], modos=[lv['modos'] for lv in niv], N=N.tolist(),
                   TZ=TZ.tolist(), nV=nV.tolist(), nVhi=nVhi.tolist(), Nf=Nf.tolist(), Y=Y.tolist(),
                   beta=beta, desc=desc, pert=self.pert, dB_L=self.fb, eps_mu=self.emu,
                   nQ=[max(self.Qc(m)[2] + self.Qc(m)[1] for m in lv['modos'] if len(lv['sel'][m])) for lv in niv],
                   SZ=[i for i, lv in enumerate(niv) if min(abs(m) for m in lv['modos']) < MZ + nk.BAND_M],
                   folga_far=a.folga_far, lam0=[complex(a.lam0).real, complex(a.lam0).imag],
                   eps_fundo_L=fundo_L.EPS_FUNDO_L, rouche=bool(getattr(a, 'rouche', False)), nQZ=z['nQZ'], Qfar=z['Qfar'],
                   K=K.tolist(), nQm=nQm.tolist())
        (a.dir/'T3.json').write_text(json.dumps(out, indent=1) + '\n')
        self.log('etapa B (nk_L3) gravada: T3.json')

    def produto_cert_Q2(self, lin, col, No):
        """X = fl(P_lin B_BANDA Q0^2 P_col) (Q0^2 no bloco inteiro de cada modo), erro(X) e max ||Q0|| nos
        modos de col: a constante K_J = ||V_J P_J B Q0(c)^2 P_J|| da sensibilidade da cauda a s (Rouche)."""
        ctx, Nc = self.ctx, self.Nc
        mc = [m for m in col if len(col[m])]; ml = [m for m in lin if len(lin[m])]
        oc = np.cumsum([0] + [len(col[m]) for m in mc]); ol = np.cumsum([0] + [len(lin[m]) for m in ml])
        X = np.zeros((int(ol[-1]), int(oc[-1])), complex)
        e2 = 0.0; dQ2 = 0.0; nQ2 = 0.0; nQm = 0.0
        for q, mi in enumerate(mc):
            Q, dq, nq = self.Qc(mi)
            Q2 = Q @ Q[:, col[mi]]
            dQ2 = max(dQ2, (nq + dq)**2 - nq**2 + gamma(Q.shape[0])*nq*nq)
            nQ2 = max(nQ2, (nq + dq)**2); nQm = max(nQm, nq + dq)
            for p_, mo in enumerate(ml):
                if abs(mo - mi) <= nk_L2.BANDA:
                    Bb = ctx.B(mo, mi, No, Nc)[lin[mo]]
                    X[ol[p_]:ol[p_+1], oc[q]:oc[q+1]] = Bb @ Q2
                    e2 += err_mm(Bb, Q2)**2
        erro = math.sqrt(e2) + self.nB*dQ2 + self.fb['arred']*nQ2 + self.pert*nQ2/max(nQm, 1e-300)*nQm
        return X, erro, nQm

    @staticmethod
    def chute_VX(P, X, it=30):
        """Estimativa de ||H^-1 X|| com P = H H^* (Cholesky): potencia em X^* P^-1 X."""
        from scipy.linalg import cho_factor, cho_solve
        cf_ = cho_factor(P, lower=True, check_finite=False)
        v = np.random.default_rng(3).standard_normal(X.shape[1]) + 0j; v /= np.linalg.norm(v); lam = 0.0
        for _ in range(it):
            w = X.conj().T @ cho_solve(cf_, X @ v, check_finite=False); lam = float(np.linalg.norm(w)); v = w/lam
        return math.sqrt(lam)

    def dados_rouche(self):
        """Para o Rouche no contorno (C1_REAVALIACAO.md): ||Q_ZZ(lambda0)|| certificado (Q_cert) e
        ||Q0(lambda0) P_far|| (formas fechadas, como em nk_L2_Z)."""
        from nk_L2_Z import Q_cert
        a, MZ, NZ, Mc, Nc = self.a, self.MZ, self.NZ, self.Mc, self.Nc
        nQZ = 0.0
        for m in nk.sl.modos(MZ):
            _, dq, nq = Q_cert(m, NZ, a.lam0); nQZ = max(nQZ, nq + dq)
        Qf = max(max(cf.cauda_Rinv(a.lam0 + 1j*m/2, Nc, p, MU, K2, 1500) for p in ((2, 1, 1) if m % 2 == 0 else (1,)))
                 for m in range(-(Mc - 1), Mc))
        Qf = max(Qf, cf.Rinv_uniforme(Mc, MU, a.lam0, K2))
        self.log(f'Rouche em lambda0 = {a.lam0}: ||Q_ZZ|| <= {nQZ:.4f}, ||Q0 P_far|| <= {Qf:.4f}')
        return dict(nQZ=nQZ, Qfar=Qf)

    def residuo_janelas(self, z, h0w):
        """Por janela: (x, erro) com x = P_r B h0 (linhas de cauda; banda completa), ou None."""
        ctx, MZ, NZ, Nc = self.ctx, self.MZ, self.NZ, self.Nc
        ms = nk.sl.modos(MZ); rs = [ctx.rad(m, NZ) for m in ms]
        offs = np.cumsum([0] + [r.dim for r in rs])
        NZb = NZ + nk.BAND_N + 1
        assert NZb <= Nc
        res = {}; er = {}
        for mo in range(-(MZ + nk.BAND_M - 1), MZ + nk.BAND_M):
            acc = np.zeros(ctx.rad(mo, NZb).dim, complex); e2 = 0.0
            for j, mi in enumerate(ms):
                if abs(mo - mi) <= nk.BAND_M:
                    Bb = ctx.B(mo, mi, NZb, NZ); h = h0w[offs[j]:offs[j+1]]
                    acc += Bb @ h; e2 += nk_L2_err(Bb, h)
            full = np.zeros(ctx.rad(mo, Nc).dim, complex); full[nk_L2.mapa(ctx, mo, NZb, Nc)] = acc
            res[mo] = full; er[mo] = math.sqrt(e2) + U*float(np.linalg.norm(acc))
        ctx._B.cache_clear()
        fb40 = fundo_L.dB_L(ctx.g, nk.BAND_M)
        pert = (fb40['arred'] + fb40['rt'])*1.0001 + DELTA_MU*(z['J1h'] + 2*cf.norma_divxi(K2))
        out = []
        for lv in self.niveis:
            if not any(m in res for m in lv['modos']):
                out.append(None); continue
            x = np.concatenate([res[m][lv['sel'][m]] if m in res else np.zeros(len(lv['sel'][m]), complex)
                                for m in lv['modos'] if len(lv['sel'][m])])
            e = math.sqrt(sum(er[m]**2 for m in lv['modos'] if m in res)) + pert
            out.append((x, e))
        return out


def nk_L2_err(B, h):
    return (gamma(B.shape[1] + 8)*float(np.linalg.norm(B))*float(np.linalg.norm(h)))**2


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--lam0', type=complex, default=None, help='padrao: o lam0 de dir/Z.json (se dado, tem de bater)')
    ap.add_argument('--rouche', action='store_true', help='Rouche no contorno: lambda0 complexo, sem etapa A nem residuo')
    ap.add_argument('--Z', default='32x128'); ap.add_argument('--F', default='200x400')
    ap.add_argument('--W', type=int, default=16)
    ap.add_argument('--folga-far', type=int, default=60)
    ap.add_argument('--estimativa', required=True, help='JSON de nk_L3.py (chutes)')
    ap.add_argument('--dir', type=Path, default=Path('build/s3b/rig'))
    ap.add_argument('--ckpt', default=None, help='arquivo de checkpoint (padrao: dir/T3.parcial.json)')
    ap.add_argument('--juntar', nargs='*', help='outros checkpoints a carregar (passada de juncao)')
    ap.add_argument('--reverso', action='store_true', help='janelas em ordem reversa (segunda maquina)')
    ap.add_argument('--so-janelas', action='store_true', help='so as janelas (sem far nem T3.json)')
    a = ap.parse_args()
    if a.rouche:
        assert a.lam0 is not None, '--rouche exige --lam0'
        a.dir.mkdir(parents=True, exist_ok=True)
    else:
        lz = json.loads((a.dir/'Z.json').read_text())['lam0']      # rastreabilidade: um so lambda0 por certificado
        lz = complex(*lz) if isinstance(lz, list) else lz            # faixa B: lambda0 complexo
        if a.lam0 is None:
            a.lam0 = lz
        assert a.lam0 == lz, f'--lam0 {a.lam0} difere do lam0 de Z.json ({lz})'
    Cauda3(a).rodar()


if __name__ == '__main__':
    main()
