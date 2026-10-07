#!/usr/bin/env python3
"""R1 (revisao independente de S3b): Rouche de operadores para K (A1).

 cobertura: a circunferencia |s - c| = rho (c = float(0,401), rho = 0,05) esta na uniao dos 16
   discos |s - c - dV_j| <= rZ do JSON. Metodo PROPRIO (diferente do do codigo, que testa as
   pontas de arcos fixos): cada disco cobre o arco de angulo [phi_j - a_j, phi_j + a_j],
   cos a_j = (rho^2 + |dV_j|^2 - rZ^2)/(2 rho |dV_j|); em Arb (200 bits) confere-se que arcos
   consecutivos se sobrepoem (volta completa). Tambem |dV_j| + rZ <= R (racionais).

 numerico (--numerico): num ponto s da circunferencia, monta EM PONTO FLUTUANTE, com codigo
   proprio, H(s) = I + (B + s - c) Q0(c) e A(s) = diag(Hc_ZZ(s), I, I) numa caixa moderada
   (Z = 20x80, G+ = 40x160, F' = 40x200 no lugar de 200x400; o far nao entra), e mede as normas
   dos blocos de A^-1 (H - A) nos niveis (Z, T1, T2'), comparando com a N' do JSON (entrada a
   entrada) e a razao de Perron max_i (N_num v)_i/v_i com o v do tile.
"""
import json
import math
import sys
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path('scripts').resolve()))


def cobertura(js):
    import flint
    flint.ctx.prec = 200
    rho = flint.arb(js['rho']); rZ = flint.arb(js['rZ'])
    pi = flint.arb.pi()
    arcos = []
    for t in js['tiles']:
        x, y = (flint.arb(v) for v in t['dV'])
        r2 = x*x + y*y; r = r2.sqrt()
        ca = (rho*rho + r2 - rZ*rZ)/(2*rho*r)
        assert ca < 1 and ca > -1
        a = ca.acos(); phi = flint.arb.atan2(y, x)
        arcos.append((phi, a))
        # tile dentro do disco R dos niveis T (racionais exatos)
        assert (Fr(js['R']) - Fr(js['rZ']))**2 >= Fr(t['dV'][0])**2 + Fr(t['dV'][1])**2
    n = len(arcos)
    ok = True; folga_min = None
    ordem = sorted(range(n), key=lambda j: float(arcos[j][0].mid()) % (2*math.pi))
    for k in range(n):
        j, j2 = ordem[k], ordem[(k + 1) % n]
        (p1, a1), (p2, a2) = arcos[j], arcos[j2]
        dp = p2 - p1                         # reduzido a (0, 2 pi): o passo angular entre vizinhos
        while float(dp.mid()) <= 0:
            dp = dp + 2*pi
        while float(dp.mid()) >= float((2*pi).mid()):
            dp = dp - 2*pi
        folga = (a1 + a2) - dp               # > 0: arcos consecutivos se sobrepoem
        ok = ok and bool(folga > 0)
        folga_min = folga if folga_min is None or float(folga.lower()) < float(folga_min.lower()) else folga_min
    print(f'cobertura (Arb): {"OK" if ok else "FALHA"}; menor sobreposicao angular entre discos vizinhos: '
          f'{float(folga_min.lower()):.3e} rad (semi-arco coberto ~ {float(arcos[0][1].mid()):.5f}, passo {2*math.pi/n:.5f})')
    return ok


def numerico(js, j=0):
    import scipy.linalg as sla
    import fourier_tail_bound as ft
    import signed_operator as so
    import tile_perron as tp
    c = js['c']; t = js['tiles'][j]
    s = c + complex(*t['dV'])
    K1, K2 = ft.K1, ft.K2
    MZ, NZ = 20, 80; MP, NP = 40, 160; MF, NF = 40, 200
    ms = list(range(-(MF - 2), MF - 1, 2))
    n1 = np.arange(1, NF, 2); n2 = np.arange(NF)
    nn = np.concatenate([n1, n2]); k = len(nn)
    w_r = np.sqrt(np.concatenate([ft.omega(n1), ft.omega(n2)]))
    g = ft.modos_campos()
    gt = tp.truncar(g, 20, 40)
    dim = len(ms)*k
    print(f'  caixa F\' = {MF}x{NF}: dim {dim}; ponto s = {s:.6f} (tile {j})', flush=True)
    # Q0(c) pesado por modo e B pesado (banda completa e truncada)
    Q = {}
    for m in ms:
        J = so.Jrad(m, c, NF)
        Q[m] = (w_r[:, None]*np.linalg.solve(J, np.eye(k)))/w_r[None, :]
    def Bbloco(gg, d):
        X = ft.B_rad(gg, d, n1, n2, n1, n2, ft.MU)
        return K1**d*(w_r[:, None]*X)/w_r[None, :]
    niv = np.zeros(dim, int)                  # 0 = Z, 1 = T1, 2 = T2'
    for ii, m in enumerate(ms):
        loc = np.where(nn < NZ, 0, np.where(nn < NP, 1, 2)) if abs(m) < MZ else \
            (np.where(nn < NP, 1, 2) if abs(m) < MP else np.full(k, 2))
        niv[ii*k:(ii+1)*k] = loc
    iZ = np.where(niv == 0)[0]
    # H = I + (B + s - c) Q0(c), por blocos (Q0 bloco-diagonal: sem matriz densa de Q0)
    H = np.zeros((dim, dim), complex)
    cache = {}
    for jj, mi in enumerate(ms):
        for ii, mo in enumerate(ms):
            d = mo - mi
            if abs(d) <= 40:
                if d not in cache:
                    cache[d] = Bbloco(g, d)
                H[ii*k:(ii+1)*k, jj*k:(jj+1)*k] = cache[d] @ Q[mi]
        H[jj*k:(jj+1)*k, jj*k:(jj+1)*k] += (s - c)*Q[mi]
    del cache
    H[np.diag_indices(dim)] += 1.0
    # Hc_ZZ = I + P_Z (B_t + s - c) Q0(c) P_Z (fundo truncado), so nos modos de Z
    msZ = [m for m in ms if abs(m) < MZ]
    locZ = np.where(nn < NZ)[0]; kz = len(locZ)
    HcZZ = np.zeros((len(msZ)*kz, len(msZ)*kz), complex)
    cache = {}
    for jj, mi in enumerate(msZ):
        QZ = Q[mi][np.ix_(locZ, locZ)]
        for ii, mo in enumerate(msZ):
            d = mo - mi
            if d not in cache:
                cache[d] = Bbloco(gt, d)[np.ix_(locZ, locZ)]
            HcZZ[ii*kz:(ii+1)*kz, jj*kz:(jj+1)*kz] = cache[d] @ QZ
        HcZZ[jj*kz:(jj+1)*kz, jj*kz:(jj+1)*kz] += (s - c)*QZ
    HcZZ[np.diag_indices(len(HcZZ))] += 1.0
    # ordem de iZ (modo a modo, n < NZ) = ordem de HcZZ
    assert len(iZ) == len(HcZZ)
    # sanidade: Q0 P_Z cai em Z
    fora = max(np.abs(Q[m][np.ix_(np.where(nn >= NZ)[0], locZ)]).max() for m in msZ)
    print(f'  max |P_(nao Z) Q0 P_Z| = {fora:.1e} (tem de ser 0: K_livre triangular superior)', flush=True)
    E = H; del H
    E[np.ix_(iZ, iZ)] -= HcZZ
    for lv in (1, 2):
        idx = np.where(niv == lv)[0]
        E[idx, idx] -= 1.0
    E[iZ] = np.linalg.solve(HcZZ, E[iZ])
    Nn = np.zeros((3, 3))
    for a in range(3):
        ia = np.where(niv == a)[0]
        for b in range(3):
            ib = np.where(niv == b)[0]
            Nn[a, b] = np.linalg.norm(E[np.ix_(ia, ib)], 2)
    Nr = np.array([[float(Fr(x)) for x in l] for l in t['N_rouche']])
    v = np.array(t['v'])
    np.set_printoptions(precision=5, suppress=False)
    print('  normas numericas dos blocos de A^-1 (H - A) (Z, T1, T2\'):\n', Nn)
    print('  cotas N\' do JSON (Z, T1, T2, far):\n', Nr)
    acima = all(Nn[a, b] <= Nr[a, b] for a in range(3) for b in range(3))
    print(f'  todas as normas numericas <= N\' correspondente: {acima}')
    th = max((Nn @ v[:3])/v[:3])
    print(f'  razao de Perron numerica (sem o far) com o v do tile: {th:.4f};  theta\' do JSON: {t["theta_rouche"]:.4f}')
    print(f'  ||A^-1 (H - A)|| na norma max_i ||x_i||/v_i (so Z, T1, T2\'), direto: '
          f'{max(sum(Nn[a, b]*v[b] for b in range(3))/v[a] for a in range(3)):.4f} (cota superior da norma induzida)')


if __name__ == '__main__':
    js = json.loads(Path('build/s3b/rouche_K_circ.json').read_text() if Path('build/s3b/rouche_K_circ.json').exists()
                    else Path('rouche_K_circ.json').read_text())
    cobertura(js)
    if '--numerico' in sys.argv:
        numerico(js, int(sys.argv[sys.argv.index('--numerico') + 1]) if len(sys.argv) > sys.argv.index('--numerico') + 1 else 0)
