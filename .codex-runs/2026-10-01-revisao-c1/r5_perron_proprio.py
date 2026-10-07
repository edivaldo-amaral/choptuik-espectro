#!/usr/bin/env python3
"""R5 (revisao independente C1): remontagem PROPRIA, em racionais, da matriz de Perron 3x3 (Z', T, far)
de alguns discos, a partir dos JSONs (rig, cauda T3, certF, schur.json, E.json), e comparacao com a
matriz gravada pelo autor (disco_*.json, campo det.M).

Derivacao (ver RELATORIO.md, R1/R5). Norma: Z' com escala D_Z = J_Z(c) Q_ZZ(s), T com
||x_T||_T = ||(w_J ||x_J||)_J||_2 (w <= 1), far simples. Entradas (B calculado = dados; dB, emu, dE como
perturbacoes):
  ZZ  = tau ((dB + dEZ) ||Q_ZZ(c)|| + emu_c)                  [J_Z(c) R(s) dL Q_ZZ(c)]
  ZT  = [E ||Phi|| + (tau + JQ) rho + tau((dB + dEZ) nQx + emu_s)]/(1 - d qTT)/wmin_SZ + tau eg nQmax/wmin0
  Zf  = tau (||B|| descidas + (dB + dEZ) Qf_s + emu_s + nE desc(NZ, Nc))
  TZ  = ||(w_i [TZ_i + nV_i (dpert + dET) nQZ])||_2 + max(w nV) eg nQZ
  TT  = ||(w_i Ns_ij / w_j)||_2 + max(w nV) eg nQmax/wmin0      (Ns: sensibilidade, minha derivacao = a do autor)
  Tf, fT, ff, fZ: como em nk_L3_C (maquinaria ja revisada)
TERMOS QUE O AUTOR NAO INCLUI (deflacao VERDADEIRA fora das linhas de Z' <- Z'/T <- Z'):
  phi* Q0(s) nao se anula nas colunas de T e do far (Q0 e triangular superior em n: Q_ZT, Q_Zfar != 0), e
  P_T x_k != 0 (cauda de RefA nos modos 32..40 + erro), P_far x_k = P_far (x_k - x_k^A):
  T <- T  : + max(w nV) dET nQx / wmin0        (dET = sum delta ||phi|| (cauda + erro))
  T <- far: + max(w nV) dET nQx
  far <- Z: + dEf ||Q_ZZ(c)||                  (dEf = sum delta ||phi|| erro)
  far <- T: + dEf nQx / wmin0 ;  far <- far: + dEf nQx
(nQx = c_w/mu >= ||Q0(s)|| para Re s >= 0, a mesma cota uniforme usada pelo autor.)"""
import json, math, sys
from fractions import Fraction as Fr
from pathlib import Path
import numpy as np

sys.path.insert(0, 'scripts')
import closed_form as cf
import fundo_L
import nk_L as nk
import signed_operator_L as sl
from rouche_L_disco import dados_centrais, qz_ponto   # q_ZJ certificados (Q_cert); conferidos em ponto flutuante abaixo
from verified_norms import U as UNIT, cota_norma

MU, K2 = nk.MU, nk.K2
ARRED = 1 + Fr(1, 10**9)
def Q(x):
    return Fr(float(x))*ARRED
def raiz(x):
    s = Fr(math.sqrt(float(x)))*(1 + Fr(1, 10**12))
    while s*s < x: s *= 1 + Fr(1, 10**9)
    return s
def n2(M):   # cota da norma espectral de matriz pequena >= 0: sqrt(||M^T M||_1) (racional)
    Mq = [[Q(x) for x in l] for l in M]
    k = len(Mq)
    G = [[sum(Mq[a][i]*Mq[a][j] for a in range(k)) for j in range(k)] for i in range(k)]
    return raiz(max(sum(G[i][j] for i in range(k)) for j in range(k)))

ROT = {0.25j: ('0_25', '+0_0000+0_2500i', 'lab2'), 0j: ('0', '+0_0000+0_0000i', 'lab2'), 1+0j: ('1_0', '+1_0000+0_0000i', 'lab')}
base = Path('build/rouche_L_rig')
info = json.load(open(base/'z/lab2/schur.json')); MZ, NZ = info['Z']; Mc, Nc = info['F']
dfr = info['defl_refA']
eE = json.load(open('build/s4/rig/E.json'))['e']
dtauB = 1.4563e-09
fb40 = fundo_L.dB_L(sl.campos(), nk.BAND_M)
errZ = [dtauB + 1e-60 + UNIT*dfr['norma_x'][0]*10, eE + UNIT*dfr['norma_x'][1]*10]
dEZ = sum(dk*nphi*e for dk, nphi, e in zip(dfr['delta'], dfr['norma_phi'], errZ))
dET = sum(dk*nphi*(cx + e) for dk, nphi, cx, e in zip(dfr['delta'], dfr['norma_phi'], dfr['cauda_x'], errZ))
nE = sum(dk*nphi*nx for dk, nphi, nx in zip(dfr['delta'], dfr['norma_phi'], dfr['norma_x']))
dEf = dEZ
nQx = (1 + 2*math.sqrt(1 + K2**-4)/(K2 - 1))/MU*1.0001
print(f'dEZ {dEZ:.3e}  dET {dET:.3e}  nE {nE:.3f}  nQx {nQx:.3f}  rt {fb40["rt"]:.2e} arred40 {fb40["arred"]:.2e}')

def norma_CW(Ns):
    """Cota racional de ||Ns||_2 para Ns >= 0 (lista de listas de Fraction): Collatz-Wielandt em Ns^T Ns."""
    a, b = len(Ns), len(Ns[0])
    G = [[sum(Ns[t][i]*Ns[t][j] for t in range(a)) for j in range(b)] for i in range(b)]
    Gf = np.array([[float(x) for x in l] for l in G])
    lam, V = np.linalg.eigh(Gf); v = np.abs(V[:, -1]); v = np.maximum(v/v.max(), 1e-6)
    vq = [Fr(float(x)) for x in v]
    lmax = max(sum(G[i][j]*vq[j] for j in range(b))/vq[i] for i in range(b))
    return raiz(lmax)


def perron(M):
    Mf = np.array([[float(x) for x in l] for l in M])
    lam, vec = np.linalg.eig(Mf); v = np.abs(vec[:, np.argmax(np.abs(lam))]); v = np.maximum(v/v.max(), 1e-9)
    vq = [Fr(float(x)) for x in v]
    return max(sum(M[i][j]*vq[j] for j in range(3))/vq[i] for i in range(3))

def disco(c, ips, saida):
    rot, lab, maq = ROT[c]
    t3 = json.load(open(base/f'cauda/T3_{rot}.json')); assert complex(*t3['lam0']) == c
    zr = json.load(open(base/f'z/{maq}/rig_{lab}_rig.json')); assert complex(*zr['c']) == c
    cF = json.load(open(base/f'F/certF_{lab}.json')); assert complex(*cF['c']) == c
    dd = json.load(open(base/f'discos/disco_{rot}.json'))
    cent = dados_centrais(t3, c, MZ, NZ, Nc)
    k = len(t3['nV'])
    for ip in ips:
        zp = zr['pontos'][ip]; dp = dd['pontos'][ip]
        p = complex(*zp['p']); assert complex(*dp['p']) == p
        r = Q(dp['r']); dq_ = Q(abs(p - c)) + r
        qz = qz_ponto(cent, p, NZ)
        tp = Q(zp['t_p']); nQZp = Q(zp['nQZ'])
        assert r*tp < 1 and r*nQZp < 1
        ts = tp/(1 - r*tp); tau = 1 + (dq_ + Q(zp['nBhat']))*ts
        JQ = 1 + dq_*nQZp/(1 - r*nQZp)
        E = Q(zp['a0']) + r*Q(zp['b1']) + r*r*Q(zp['b2']) + r**3*(tau*Q(zp['nY3']) + JQ*Q(zp['nQ3F']))
        qTT = Q(cF['qTT']); assert dq_*qTT < 1
        EF = min(E*Q(pr['phi']) + (tau + JQ)*Q(pr['rho']) for pr in cF['pares'])
        emu_s = Q(fundo_L.eps_mu_L(complex(abs(p.real) + float(r), abs(p.imag) + float(r))))
        emu_c = Q(t3['eps_mu']); rt = Q(fb40['rt']); dB = Q(fb40['arred']) + rt
        nQZc = Q(t3['nQZ']); Qfar = Q(t3['Qfar']); Qf_s = Qfar/(1 - dq_*Qfar)
        if (c, ip) in QF_OVERRIDE:                                     # r5c: cota VALIDA de ||Q0(s) P_far|| no disco
            print(f'   Qf_s do autor {float(Qf_s):.6f} -> substituido por {QF_OVERRIDE[(c, ip)]:.6f}')
            Qf_s = Q(QF_OVERRIDE[(c, ip)])
        eg = rt + emu_s
        # T: sensibilidade (minha derivacao; ver relatorio)
        N = t3['N']; TZ = t3['TZ']; nV = t3['nV']; nVhi = t3['nVhi']; K = t3['K']; nQm = t3['nQm']; nQl = t3['nQ']
        nQJJ = [Q(cent[i]['nQJJ']) if i in cent else Q(nQm[i]) for i in range(k)]
        assert all(dq_*x < 1 for x in nQJJ)
        fJ = [1/(1 - dq_*x) for x in nQJJ]
        qZs = {i: Q(qz[i])/(1 - r*nQZp)*fJ[i] for i in cent}            # (1 + d nQ fJ) = fJ
        dpert = max(Fr(0), emu_s - emu_c); nQmax = max(Q(x) for x in nQl)
        Ns = [[Fr(0)]*k for _ in range(k)]
        for i in range(k):
            if i in cent:
                Ns[i][i] = Q(N[i][i]) + dq_*(Q(TZ[i])*qZs[i] + (Q(K[i]) + Q(TZ[i])*Q(cent[i]['qZc']))*fJ[i])
            else:
                Ns[i][i] = Q(N[i][i]) + dq_*Q(K[i])*fJ[i]
            Ns[i][i] += Q(nV[i])*dpert*nQmax
            for j in range(k):
                if j != i and N[i][j] != 0:
                    Ns[i][j] = Q(N[i][j]) + dq_*((Q(TZ[i])*qZs[j] if j in cent else 0) + Q(N[i][j])*nQJJ[j]*fJ[j]) + Q(nV[i])*dpert*nQmax
        w = [Fr(1)]*k                                                  # pesos 1 (os que o autor escolheu nestes discos)
        wmin0 = min(w); wminSZ = min(w[j] for j in t3['SZ'])
        maxwV = max(w[i]*Q(nV[i]) for i in range(k))
        # ||Ns||_2 rigorosa: Ns >= 0 entrada a entrada, G = Ns^T Ns >= 0 exata em racionais e, por
        # Collatz-Wielandt, lambda_max(G) <= max_i (G v)_i / v_i para qualquer v > 0 (v = Perron de G em fl)
        nT = norma_CW(Ns)
        nT_fl = float(np.linalg.norm(np.array([[float(x) for x in l] for l in Ns]), 2)); nT_cru = nT
        nT += maxwV*eg*nQmax/wmin0
        TZs = [Q(TZ[i]) + Q(nV[i])*(dpert + Q(dET))*nQZc for i in range(k)]
        TZw = raiz(sum((w[i]*TZs[i])**2 for i in range(k))) + maxwV*eg*nQZc
        ZT = (EF + tau*((dB + Q(dEZ))*Q(nQx) + emu_s))/(1 - dq_*qTT)/wminSZ + tau*eg*nQmax/wmin0
        ZZ = tau*((dB + Q(dEZ))*nQZc + emu_c)
        Zf = tau*(Q(nk.NORMA_BL)*Q(cf.descida(NZ + nk.BAND_N + 1, Nc, MU, K2) + 6*cf.descida_divxi(NZ, Nc, K2))
                  + (dB + Q(dEZ))*Qf_s + emu_s + Q(nE)*Q(cf.descida(NZ, Nc, MU, K2)))
        # far e T<-far (maquinaria de nk_L3_C, reproduzida)
        fb_tot = Q(t3['dB_L']['total']); resto = Q(t3['dB_L']['resto'])
        beta_s = (Q(nk.NORMA_BL) + rt)*Qf_s + emu_s
        lim = Nc - nk.BAND_N - t3.get('folga_far', 60)
        desc_s = (Q(nk.NORMA_BL) + fb_tot)*Q(cf.descida(lim + nk.BAND_N + 1, Nc, MU, K2) + 6*cf.descida_divxi(lim, Nc, K2)) + (fb_tot + emu_s)*Qf_s
        borda = [max(abs(m) for m in t3['modos'][i]) >= Mc - nk.BAND_M for i in range(k)]
        nQb = max(Q(nQl[i]) for i in range(k) if max(abs(m) for m in t3['modos'][i]) >= Mc - 24)
        Tf = max(w[i]*Q(nV[i] if borda[i] else nVhi[i]) for i in range(k))*beta_s \
            + raiz(sum((w[i]*Q(nV[i]))**2 for i in range(k)))*desc_s + maxwV*eg*Qf_s
        Nf = np.array(t3['Nf'])*np.array([float(x) for x in fJ])[None, :]
        fT = norma_CW([[Q(t3['Nf'][a][b])*fJ[b] for b in range(k)] for a in range(len(t3['Nf']))]) + resto*nQb/wmin0 + eg*nQmax/wmin0
        fZ = (fb_tot + emu_s)*nQZc
        M_aut = [[ZZ, ZT, Zf], [TZw, nT, Tf], [fZ, fT, beta_s]]
        th_aut = perron(M_aut)
        # termos de deflacao que faltam
        M_ext = [[ZZ, ZT, Zf],
                 [TZw, nT + maxwV*Q(dET)*Q(nQx)/wmin0, Tf + maxwV*Q(dET)*Q(nQx)],
                 [fZ + Q(dEf)*nQZc, fT + Q(dEf)*Q(nQx)/wmin0, beta_s + Q(dEf)*Q(nQx)]]
        th_ext = perron(M_ext)
        Mreg = dp['det']['M']
        dif = max(abs(float(M_aut[i][j]) - Mreg[i][j])/max(Mreg[i][j], 1e-300) for i in range(3) for j in range(3))
        lin = (f'c={c} p={p:.6f} r={float(r):.4e}: theta_autor(JSON) {dp["theta_f"]:.6f} | minha remontagem {float(th_aut):.6f}'
               f' (max dif. relativa das entradas {dif:.1e}) | com termos de deflacao faltantes {float(th_ext):.6f}'
               f' (delta {float(th_ext - th_aut):.1e})')
        print(lin, flush=True)
        print('   M (minha)   = ' + str([[f'{float(x):.6g}' for x in l] for l in M_aut]))
        print('   M (autor)   = ' + str([[f'{x:.6g}' for x in l] for l in Mreg]))
        print(f'   tau {float(tau):.2f} E {float(E):.6f} EF {float(EF):.6f} nT: fl {nT_fl:.6f} CW racional {float(nT_cru):.6f}')
        saida.append(dict(c=[c.real, c.imag], p=[p.real, p.imag], r=float(r), theta_autor=dp['theta_f'], theta_meu=str(th_aut),
                          theta_meu_f=float(th_aut), theta_com_deflacao=str(th_ext), theta_com_deflacao_f=float(th_ext),
                          M_meu=[[float(x) for x in l] for l in M_aut], M_autor=Mreg, dif_rel=dif))
        # conferencia em ponto flutuante PROPRIA de qz_J(p) = max_m ||Q_ZZ,m(p) J_ZJ,m Q_JJ,m(c)|| (janelas centrais)
        if ip == ips[0]:
            for i in list(cent)[:2]:
                v = 0.0
                for m in [m for m in t3['modos'][i] if abs(m) < MZ]:
                    rr = sl.Radial(m, Nc); wv = np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in rr.comps]))
                    nn = np.concatenate([ns for _, ns in rr.comps]); iz = np.where(nn < NZ)[0]; iT = np.where(nn >= NZ)[0]
                    Jp = (wv[:, None]*sl.J_livre(rr, p))/wv[None, :]; Jc = (wv[:, None]*sl.J_livre(rr, c))/wv[None, :]
                    Qp = np.linalg.inv(Jp[np.ix_(iz, iz)]); Qc = np.linalg.inv(Jc[np.ix_(iT, iT)])
                    v = max(v, float(np.linalg.norm(Qp @ Jc[np.ix_(iz, iT)] @ Qc, 2)))
                print(f'   qz janela {i}: autor (cert) {qz[i]:.6f}, meu (fl) {v:.6f}')

QF_OVERRIDE = {}
if __name__ == '__main__':
  saida = []
  for c, ips in ((0.25j, [0, 8, 15]), (0j, [6, 10]), (1+0j, [0, 3])):
      disco(c, ips, saida)
      json.dump(saida, open('.codex-runs/2026-10-01-revisao-c1/r5_resultado.json', 'w'), indent=1)
