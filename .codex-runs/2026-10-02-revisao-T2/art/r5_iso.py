#!/usr/bin/env python3
"""R5 (revisao T2): L-iso e L-constr NOS BLOCOS DO CODIGO.

(a) Os blocos de B (signed_operator_L.bloco_B) dependem so de d = m_out - m_in e das componentes
    (para m_in de mesma paridade): comparacao exata entre varios m.
(b) Setor B construido com o MESMO montador (componentes 0,1,2 nos m impares, 3 nos pares):
    L_B(s) nos modos (m + 1) coincide com L_A(s + i/2) nos modos m, bloco a bloco (J e B).
    E C_B(s)(m + 1) = C_A(s + i/2)(m) para o operador de constraints (signed_constraints.bloco_C).
(c) Norma de M (deslocamento m -> m + 1) no toro com sinal: ||M x|| = kappa1 ||x|| exatamente.
(d) i-periodicidade: L_A(s + i) nos modos m + 2 = L_A(s) nos modos m.
"""
import sys
import json
import numpy as np

sys.path.insert(0, 'scripts')
import signed_operator_L as sl
import signed_constraints as sc
import closed_form as cf

K1, K2 = 65/64, 5/4


class RadialB(sl.Radial):
    """Setor B: componentes 0 (n impar), 1, 2 nos m IMPARES; componente 3 nos m PARES."""
    def __init__(self, m, N):
        self.m = m
        if m % 2 != 0:
            self.comps = [(0, np.arange(1, N, 2)), (1, np.arange(N)), (2, np.arange(N))]
        else:
            self.comps = [(3, np.arange(N))]
        self.off = {}
        o = 0
        for c, ns in self.comps:
            self.off[c] = (o, o + len(ns)); o += len(ns)
        self.dim = o


class SaidaB(sc.Saida):
    def __init__(self, m, N):
        super().__init__(m, N)


def main():
    g = sl.campos()
    N = 40
    out = {}
    # (a)
    pior_a = 0.0
    for d in (-40, -17, -3, 0, 1, 2, 9, 40):
        for par in (0, 1):
            ref = sl.bloco_B(g, sl.Radial(par + d, N), sl.Radial(par, N))
            for m in (par - 6, par + 4, par + 12, par - 30):
                X = sl.bloco_B(g, sl.Radial(m + d, N), sl.Radial(m, N))
                pior_a = max(pior_a, float(np.abs(X - ref).max()))
    out['a_blocos_so_dependem_de_d'] = pior_a
    print('(a) max |B(m+d,m) - B(par+d,par)| =', pior_a, flush=True)
    # (b) L_B(s) em (m+1) contra L_A(s + i/2) em m
    pior_bJ = pior_bB = pior_bC = 0.0
    for s in (0.3, 0.0414 + 0.0j, 0.7 - 0.2j):
        for m in (-5, -2, 0, 1, 3, 8):
            for d in (0, 1, -1, 2, 7, -12):
                rA_i, rA_o = sl.Radial(m, N), sl.Radial(m + d, N)
                rB_i, rB_o = RadialB(m + 1, N), RadialB(m + 1 + d, N)
                assert [c for c, _ in rA_i.comps] == [c for c, _ in rB_i.comps]
                if d == 0:
                    JA = sl.J_livre(rA_i, s + 0.5j); JB = sl.J_livre(rB_i, s)
                    pior_bJ = max(pior_bJ, float(np.abs(JA - JB).max()))
                BA = sl.bloco_B(g, rA_o, rA_i); BB = sl.bloco_B(g, rB_o, rB_i)
                pior_bB = max(pior_bB, float(np.abs(BA - BB).max()))
    # C: a saida de K no setor A vive nos m pares (comp. 0 e 1); no setor B, nos m impares.
    for s in (0.3, 0.0414 + 0.0j, 0.7 - 0.2j):
        for m in (-4, 0, 2, 6):           # m par: saida A
            for d in (0, 1, -1, 2, -3, 5):
                mi = m - d
                CA = sc.bloco_C(g, sc.Saida(m, N), sl.Radial(mi, N), s + 0.5j)
                CB = sc.bloco_C(g, SaidaB(m + 1, N), RadialB(mi + 1, N), s)
                pior_bC = max(pior_bC, float(np.abs(CA - CB).max()))
    out['b_J_B_menos_J_A'] = pior_bJ; out['b_B_B_menos_B_A'] = pior_bB; out['b_C_B_menos_C_A'] = pior_bC
    print('(b) J: max |J_B(s)(m+1) - J_A(s+i/2)(m)| =', pior_bJ)
    print('    B: max |B_B(m+1+d, m+1) - B_A(m+d, m)| =', pior_bB)
    print('    C: max |C_B(s)(m+1 <- mi+1) - C_A(s+i/2)(m <- mi)| =', pior_bC, flush=True)
    # (c) norma de M
    rng = np.random.default_rng(3)
    ms = list(range(-9, 10))
    x = {m: rng.standard_normal(sl.Radial(m, N).dim) + 1j*rng.standard_normal(sl.Radial(m, N).dim) for m in ms}

    def norma(xd, setor):
        t = 0.0
        for m, v in xd.items():
            r = (sl.Radial if setor == 'A' else RadialB)(m, N)
            w = K1**m*np.sqrt(np.concatenate([cf.omega(ns, K2) for _, ns in r.comps]))
            t += float(np.sum(np.abs(w*v)**2))
        return t**0.5
    Mx = {m + 1: v for m, v in x.items()}
    razao = norma(Mx, 'B')/norma(x, 'A')
    out['c_norma_M'] = razao
    print('(c) ||Mx||/||x|| =', razao, ' kappa1 =', K1, flush=True)
    # (d) i-periodicidade
    pior_d = 0.0
    for s in (0.3, 0.5 + 0.1j):
        for m in (-4, 0, 3):
            r0, r2 = sl.Radial(m, N), sl.Radial(m + 2, N)
            pior_d = max(pior_d, float(np.abs(sl.J_livre(r0, s + 1j) - sl.J_livre(r2, s)).max()))
    out['d_i_periodicidade_J'] = pior_d
    print('(d) max |J_A(s+i)(m) - J_A(s)(m+2)| =', pior_d)
    json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else 'r5_iso.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
