#!/usr/bin/env python3
"""R4 (revisao independente C1): teste numerico das tres formulas de sensibilidade da cauda (A4) numa
truncagem PEQUENA do operador verdadeiro (fundo de RefA, pesos do toro), em ponto flutuante:
  janela central J:     ||V_J P_J B (Q0(s) - Q0(c)) P_J|| <= d [TZ_J q_ZJ(s) + (K_J + TZ_J q_ZJ(c)) fJ]
  janela nao central J: ||V_J P_J B (Q0(s) - Q0(c)) P_J|| <= d K_J fJ
  fora da diagonal:     ||V_i P_i B Q0(s) P_j|| <= N_ij(c) + d [TZ_i q_Zj(s) + N_ij(c) nQ_jj fJ_j]
com TZ_J = ||V_J P_J B P_Z Q_ZZ(c)||, K_J = ||V_J P_J B Q0(c)^2 P_J||, q_ZJ = ||Q_ZJ||, fJ = 1/(1 - d ||Q_JJ(c)||),
V_J = (I + P_J B Q0(c) P_J)^-1. Tambem confere a IDENTIDADE (diferenca entre os dois lados, em norma).
Truncagem: modos |m| < Mc = 10, n < Nc = 40; Z = |m| < 4, n < 12; janelas de 4 modos (T = F \\ Z)."""
import sys
import numpy as np
sys.path.insert(0, 'scripts')
import nk_L as nk
import signed_operator_L as sl

MZ, NZ, Mc, Nc = 4, 12, 10, 40
c = 0.25j
ctx = nk.Contexto(c, Nc)
ms = list(range(-(Mc - 1), Mc))
rads = {m: sl.Radial(m, Nc) for m in ms}
off = np.cumsum([0] + [rads[m].dim for m in ms]); pos = {m: (off[i], off[i+1]) for i, m in enumerate(ms)}
n = int(off[-1])
B = np.zeros((n, n), complex)
for mo in ms:
    for mi in ms:
        if abs(mo - mi) <= nk.BAND_M:
            B[pos[mo][0]:pos[mo][1], pos[mi][0]:pos[mi][1]] = ctx.B(mo, mi, Nc, Nc)
def Q0(s):
    Q = np.zeros((n, n), complex)
    for m in ms:
        r = rads[m]; w = nk.pesos(r)
        J = (w[:, None]*sl.J_livre(r, s))/w[None, :]
        Q[pos[m][0]:pos[m][1], pos[m][0]:pos[m][1]] = np.linalg.inv(J)
    return Q
# indices de Z e das janelas
nloc = {m: np.concatenate([ns for _, ns in rads[m].comps]) for m in ms}
iZ = np.concatenate([pos[m][0] + np.where(nloc[m] < NZ)[0] for m in ms if abs(m) < MZ])
def idx_jan(modos):
    out = []
    for m in modos:
        loc = np.where(nloc[m] >= NZ)[0] if abs(m) < MZ else np.arange(rads[m].dim)
        out.append(pos[m][0] + loc)
    return np.concatenate(out)
janelas = {'-3..0 (central)': range(-3, 1), '1..3 (central)': range(1, 4), '4..7': range(4, 8), '-7..-4': range(-7, -3),
           '8..9': range(8, 10), '-9..-8': range(-9, -7)}
J = {k: idx_jan(v) for k, v in janelas.items()}
central = {k: min(abs(m) for m in v) < MZ for k, v in janelas.items()}
nrm = lambda X: float(np.linalg.norm(X, 2))
Qc = Q0(c)
V = {k: np.linalg.inv(np.eye(len(ix)) + B[np.ix_(ix, ix)] @ Qc[np.ix_(ix, ix)]) for k, ix in J.items()}
pior = 0.0
for d_ in (0.02, 0.05, 0.1):
    for ang in np.linspace(0, 2*np.pi, 7)[:-1]:
        s = c + d_*np.exp(1j*ang)
        if s.real < 0:
            continue
        Qs = Q0(s); d = abs(s - c)
        for k, ix in J.items():
            Vj = V[k]
            verd = nrm(Vj @ B[np.ix_(ix, np.arange(n))] @ (Qs - Qc)[:, ix])
            K = nrm(Vj @ B[np.ix_(ix, np.arange(n))] @ (Qc @ Qc)[:, ix])
            nQJ = nrm(Qc[np.ix_(ix, ix)]); fJ = 1/(1 - d*nQJ)
            if central[k]:
                TZ = nrm(Vj @ B[np.ix_(ix, iZ)] @ Qc[np.ix_(iZ, iZ)])
                qs = nrm(Qs[np.ix_(iZ, ix)]); qc = nrm(Qc[np.ix_(iZ, ix)])
                cota = d*(TZ*qs + (K + TZ*qc)*fJ)
                # identidade (sem normas): lado direito
                QJJs = Qs[np.ix_(ix, ix)]
                rhs = -(s - c)*(Vj @ B[np.ix_(ix, iZ)] @ Qc[np.ix_(iZ, iZ)] @ Qs[np.ix_(iZ, ix)]
                                + (Vj @ B[np.ix_(ix, np.arange(n))] @ (Qc @ Qc)[:, ix] - Vj @ B[np.ix_(ix, iZ)] @ Qc[np.ix_(iZ, iZ)] @ Qc[np.ix_(iZ, ix)])
                                @ (np.eye(len(ix)) - (s - c)*QJJs))
            else:
                cota = d*K*fJ
                rhs = -(s - c)*(Vj @ B[np.ix_(ix, np.arange(n))] @ (Qc @ Qc)[:, ix]) @ (np.eye(len(ix)) - (s - c)*Qs[np.ix_(ix, ix)])
            lhs = Vj @ B[np.ix_(ix, np.arange(n))] @ (Qs - Qc)[:, ix]
            erro_id = nrm(lhs - rhs)/max(nrm(lhs), 1e-300)
            pior = max(pior, verd/cota)
            print(f'd={d:.3f} s={s:.4f} janela {k:16s}: verdadeiro {verd:.3e} <= cota {cota:.3e} (razao {verd/cota:.3f}); erro rel. da identidade {erro_id:.1e}')
        # fora da diagonal: i = '4..7' <- j = '1..3 (central)'  e  i = '-3..0' <- j = '1..3'
        for ki, kj in (('4..7', '1..3 (central)'), ('-3..0 (central)', '1..3 (central)'), ('1..3 (central)', '4..7')):
            ixi, ixj = J[ki], J[kj]; Vi = V[ki]
            Nc_ = nrm(Vi @ B[np.ix_(ixi, np.arange(n))] @ Qc[:, ixj]); Ns_ = nrm(Vi @ B[np.ix_(ixi, np.arange(n))] @ Qs[:, ixj])
            TZi = nrm(Vi @ B[np.ix_(ixi, iZ)] @ Qc[np.ix_(iZ, iZ)])
            qZj = nrm(Qs[np.ix_(iZ, ixj)]) if central[kj] else 0.0
            nQjj = nrm(Qc[np.ix_(ixj, ixj)]); fj = 1/(1 - d*nQjj)
            cota = Nc_ + d*(TZi*qZj + Nc_*nQjj*fj)
            print(f'   N[{ki} <- {kj}]: N(s) {Ns_:.4e} <= cota {cota:.4e}  (N(c) {Nc_:.4e}; razao do incremento {(Ns_ - Nc_)/(cota - Nc_):.3f})')
            pior = max(pior, (Ns_ - Nc_)/(cota - Nc_) if cota > Nc_ else 0)
print(f'pior razao verdadeiro/cota (diagonal) e incremento/cota (fora): {pior:.3f}  ->  {"OK" if pior <= 1 else "VIOLACAO"}')
