#!/usr/bin/env python3
"""R5c (revisao independente C1): cota VALIDA de sup_{|s - p| <= r, Re s >= 0} ||Q0(s) P_far|| pela forma fechada
(closed_form.cauda_Rinv reescrita com majorantes no disco) e reavaliacao de theta com ela no lugar de
Qf_s = Qf(c)/(1 - d Qf(c)) (cuja derivacao nao cobre as linhas de Z/T de Q0(s) P_far).
Forma fechada (closed_form): (J + sigma)^-1_{jc} = -(1/d_j)(v_c/d_c) prod_{j<k<c}(1 - v_k/d_k), d_k = mu(n_k + 1) + sigma.
No disco sigma = sigma_p + e, |e| <= r, Re sigma >= Re(s) >= 0:
  |d_k| >= max(|d_k(p)| - r, mu(n_k + 1))     (Re sigma >= 0 da a segunda),
  |1 - v_k/d_k| = |d_k - v_k|/|d_k| <= min(1, (|d_k(p) - v_k| + r)/|d_k|_inf)   (o 1 vale para Re sigma >= 0).
Resto (cauda n > Nbig, Rinv_uniforme para |m| >= Mc) como no original (uniformes em Re sigma >= 0; h = |m|/2 - |Im s|)."""
import math, sys, json
import numpy as np
sys.path.insert(0, 'scripts')
import closed_form as cf

MU = 722873400/2**32; K2 = 5/4; Nc = 400; Mc = 200

def cauda_disco(sig_p, r, N, passo, mu=MU, k2=K2, Nbig=1500):
    ns = np.arange(Nbig + 1) if passo == 1 else np.arange(1, 2*Nbig + 2, 2)
    dkp = mu*(ns + 1.0) + sig_p
    vk = 2*mu*ns
    dinf = np.maximum(np.abs(dkp) - r, mu*(ns + 1.0))
    rho = np.minimum(1.0, (np.abs(dkp - vk) + r)/dinf)*(1 + 1e-12)
    rho[0] = 1.0
    lw = np.where(ns == 0, 0.0, ns*math.log(k2) + 0.5*np.log1p(k2**(-4.0*ns)))
    S = np.cumsum(np.log(np.maximum(rho, 1e-300)))
    Cmax = 0.0; rows_acc = np.zeros(len(ns))
    for c in np.where(ns >= N)[0]:
        j = np.arange(c + 1)
        prod = np.exp(np.where(j < c, S[c-1] - S[j], 0.0)) if c > 0 else np.ones(1)
        ent = (1/dinf[j])*(vk[c]/dinf[c])*prod
        ent[c] = 1/dinf[c]
        wt = ent*np.exp(lw[j] - lw[c])
        Cmax = max(Cmax, float(wt.sum())); rows_acc[:c+1] += wt
    cw = math.sqrt(1 + k2**-4); nb = ns[-1]
    Ctail = (1/(mu*(nb + 1)) + 2*cw*(1/((nb/2 + 1)*(k2 - 1)) + k2**(-nb/2)*(1 + math.log(nb + 1)))/mu)
    Cmax = max(Cmax, Ctail)
    rows_acc += (2/(mu*(ns + 1)))*cw*k2**(ns - nb)/(k2 - 1)
    Rmax = max(float(rows_acc.max()), 2*cw*k2/(mu*(nb + 1)*(k2 - 1)))
    return math.sqrt(Cmax*Rmax)*1.0001

def Qfar_disco(p, r):
    q = 0.0
    for m in range(-(Mc - 1), Mc):
        for passo in ((2, 1) if m % 2 == 0 else (1,)):
            q = max(q, cauda_disco(p + 1j*m/2, r, Nc, passo))
    cw = math.sqrt(1 + K2**-4)
    q = max(q, (1 + 2*cw/(K2 - 1))/(Mc/2 - abs(p.imag) - r)*1.0001)
    return q

if __name__ == '__main__':
    import r5_perron_proprio as r5
    casos = [(0.25j, 'build/rouche_L_rig/discos/disco_0_25.json', [0, 8, 15]), (0j, 'build/rouche_L_rig/discos/disco_0.json', [6, 10]),
             (1+0j, 'build/rouche_L_rig/discos/disco_1_0.json', [0, 3])]
    for c, arq, ips in casos:
        pts = json.load(open(arq))['pontos']
        for ip in ips:
            p = complex(*pts[ip]['p']); r = pts[ip]['r']
            r5.QF_OVERRIDE[(c, ip)] = Qfar_disco(p, r)
            print(f'c={c} p={p:.4f} r={r:.3e}: sup_disco ||Q0(s) P_far|| <= {r5.QF_OVERRIDE[(c, ip)]:.6f}', flush=True)
    saida = []
    for c, arq, ips in casos:
        r5.disco(c, ips, saida)
    json.dump(saida, open('.codex-runs/2026-10-01-revisao-c1/r5c_resultado.json', 'w'), indent=1)
