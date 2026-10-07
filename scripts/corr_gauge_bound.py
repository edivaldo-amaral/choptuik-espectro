#!/usr/bin/env python3
"""S4, lema de L2/L3 (revisao de S4): cota de ||(Z* + 1) Corr||_L, com Corr = omega* - (RefA + RefB)
a correcao da bola de RT, Z* = mu*^-1 d_tau + xi d_xi. Tudo em racionais, a partir das constantes
PUBLICADAS por RT (arXiv 1203.3766v1, choptuik.tex; rotulos entre colchetes) e da norma de RefA.dat.

Espaco e norma de RT: ||v||_V = sum_{m,n in Z} k1^|m| k2^|n| ||v_mn||_C, ||z||_C = |Re z| + |Im z|,
||(v1..v4)||_SSPACE = soma; bola [kdkkkskkskksksks]: kappa_* |Corr_mu| + ||Corr_omega|| <= R = 2^-277,
kappa_* = 1.

(1) Equacao da correcao. Por [fdrrihifuhr78], C = S Omega+(mu*, omega*) - S Omega+(Ref) com
    S Omega+(mu*, omega*) = 0 e O_mu = d_tau + mu D_O afim em mu:
      O_{mu*} Corr = G := - S Omega+(Ref) - Corr_mu (D_O + S Gamma1) Ref - mu* S Gamma1 Corr
                          - 2 S Xi Gamma2(Ref, Corr) - S Xi Gamma2(Corr, Corr),
    ||G|| <= L3 + R (||D_O Ref|| + K1 ||Ref||) + mu* K1 R + 2 (L1 + 3 L2) R + 3 R^2
    ([dkjfhk1-3], [fdhkjhjdh4z894z89rhiu], [dkhfdkdhk3iuzi3i]; Ref = RefA + RefB, ||RefB|| <= L2).
    ||D_O v|| <= 9 N ||v|| para v com suporte em n < N (D_OA = BandI_{1,-1}(Band_{0,n+1} + Band_{1,n}),
    D_OB = BandI_{2,-1}(Band_{0,n+1} + Band_{2,n+1}) [dkhddjjjjjff], cotas de banda [fdkhdfkjhk464huih]).
(2) d_tau O_mu^-1 = BandI[k, f] Band[0, (im/2)/(im/2 + mu(n+1))] (1 + Band_{k,-1}) (k = 1 em O_A, 2 em O_B;
    d_tau comuta com as bandas): ||d_tau O^-1|| <= 2 (1 - k2^-1)^-1 (1 + k2^-1) = 18 =: c_tau.
(3) Componente B (1): D_OB = xi d_xi + 1, logo mu*(Z* + 1) = O_B e (Z* + 1) Corr_1 = G_1/mu*.
    Componentes A (2..4): D_OA = (1 + xi) d_xi + 1. Com h = (1 + xi) d_xi Corr_j = (G_j - d_tau Corr_j)/mu*
    - Corr_j: ||h|| <= (1 + c_tau) ||G_j||/mu* + ||Corr_j||. Como (1 + xi) d_xi = (2 BandI_{1,-1} - 1) Band_{0,n}
    [3765444f4jgjgh] = (1 + S)(1 - S)^-1 n, S = shift n -> n + 1: ||n Corr_j|| <= ||1 - S|| ||(1 + S)^-1|| ||h||
    <= (1 + k2^-1)(1 - k2^-1)^-1 ||h|| = 9 ||h||. E xi d_xi = D_OB - 1: ||xi d_xi v|| <= c_B (||n v|| + ||v||)
    + ||v||, c_B = (1 - k2^-2)^-1 (1 + k2^-2). Entao
      ||(Z* + 1) Corr_j|| <= c_tau ||G_j||/mu* + ||xi d_xi Corr_j|| + ||Corr_j||.
(3b) Fase (Rouche de L, C1_REAVALIACAO.md §9): Corr = O_{mu*}^-1 G, logo d_tau Corr = (d_tau O^-1) G e
    ||d_tau Corr||_SSPACE <= c_tau ||G|| (o mesmo (2), com G em todas as componentes).
(4) Norma de L: ||x||_L <= sqrt(2) ||x||_SSPACE (pesos com sinal k1^m <= k1^|m|, sqrt(omega2(n)) <= sqrt(2) k2^n,
    |z| <= ||z||_C; as coordenadas de L sao os coeficientes simetricos de RT).
(5) Dominio (L3): g* = g_A + e com ||e||_L finito, logo g* em Y. A identidade L_*(mu*) g* = 0 vale coeficiente a
    coeficiente, entao (J + lambda0) g_n = -B_* g_n + (lambda0 - mu*) g_n esta em Y; J + lambda0 e injetivo no
    dominio maximo de Y (as solucoes de (O + lambda0)u = 0 sao (1 + xi)^-(1 + (lambda0 + im/2)/mu) em O_A e
    xi^-(...) em O_B, singulares dentro da elipse), logo Q0 y_g = g_n e y_g = Q0^-1 g_n esta em Y."""
from __future__ import annotations

import json
from fractions import Fraction as Fr
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gauge_refB import ler

K1_, K2_ = Fr(65, 64), Fr(5, 4)
R = Fr(1, 2**277)                        # [kdkkkskkskksksks], kappa_* = 1
L1 = Fr(21, 2)                           # ||S Xi Gamma2(RefA, .)||   [dkjfhk1]
L2 = Fr(1, 2**25)                        # ||RefB||                   [dkjfhk2]
L3 = Fr(1, 2**294)                       # ||S Omega+(RefA + RefB)||  [dkjfhk3]
K1 = Fr(80, 9)                           # ||S Gamma1||               [fdhkjhjdh4z894z89rhiu]
MU_A = Fr(722873400, 2**32)
MU_B = Fr(43, 2**40)                     # |RefB_mu| <= 43 2^-40
SQRT2 = Fr(14143, 10000)                 # >= sqrt(2)


def norma_RT(path):
    """||RefA||_SSPACE exata (racional) e o suporte em n (N = max n + 1)."""
    campos, _ = ler(path)
    tot = Fr(0); N = 0
    for h, pares in campos:
        esc = Fr(2)**h['TwoExp']
        for i in range(h['num_m']):
            m = h['off_m'] + i
            wm = (2 - (m == 0))*K1_**m
            for k in range(h['num_n']):
                n = h['off_n'] + k
                re, im = pares[i*h['num_n'] + k]
                a = abs(int(re)) + abs(int(im))
                if a:
                    tot += wm*(2 - (n == 0))*K2_**n*a*esc
                    N = max(N, n + 1)
    return tot, N


def main():
    nA, NA = norma_RT('.cache/rt-1203.3766v1/sourcecode/RefA.dat')
    hAB = ler('.cache/rt-1203.3766v1/RefAplusB.dat')[0][0][0]
    NAB = hAB['off_n'] + hAB['num_n']
    nRef = nA + L2
    DORef = 9*NA*nA + 9*NAB*L2                                  # ||D_O RefA|| + ||D_O RefB||
    mu_max = MU_A + MU_B + R; mu_min = MU_A - MU_B - R
    G = L3 + R*(DORef + K1*nRef) + mu_max*K1*R + 2*(L1 + 3*L2)*R + 3*R*R
    c_tau = 2*(1 - 1/K2_)**-1*(1 + 1/K2_)
    c_B = (1 - K2_**-2)**-1*(1 + K2_**-2)
    c_n = (1 + 1/K2_)*(1 - 1/K2_)**-1
    # componentes A: h, ||n Corr||, ||xi d_xi Corr||; componente B: G_1/mu*. Soma com ||G_j|| <= ||G||, ||Corr_j|| <= R
    h = (1 + c_tau)*G/mu_min + R
    xdx = c_B*(c_n*h + R) + R
    ZA = c_tau*G/mu_min + xdx + R
    Zc = G/mu_min + 3*ZA                                       # soma sobre as 4 componentes (cota grosseira)
    ZL = SQRT2*Zc
    print(f'||RefA||_SSPACE = {float(nA):.6f} (suporte n < {NA}); RefAplusB: n < {NAB}')
    print(f'||D_O Ref|| <= {float(DORef):.4e};  ||G|| <= {float(G):.3e}  (L3 = 2^-294 = {float(L3):.2e}, R = 2^-277 = {float(R):.2e})')
    print(f'c_tau = {float(c_tau):.4g}, c_n = {float(c_n):.4g}, c_B = {float(c_B):.4g}; mu* >= {float(mu_min):.12f}')
    print(f'||(Z* + 1) Corr||_SSPACE <= {float(Zc):.3e};  ||(Z* + 1) Corr||_L <= {float(ZL):.3e}   (usado: 1e-60)')
    ok = ZL <= Fr(1, 10**60)
    print('OK: a cota 1e-60 de nk_L3_E.py vale' if ok else 'FALHA: 1e-60 nao cobre')
    TL = SQRT2*c_tau*G                                         # ||d_tau Corr||_L (fase; rouche_L_disco usa 1e-60)
    ok_t = TL <= Fr(1, 10**60)
    print(f'||d_tau Corr||_L <= sqrt(2) c_tau ||G|| = {float(TL):.3e}  (usado em rouche_L_disco: 1e-60) '
          + ('OK' if ok_t else 'FALHA'))
    ok = ok and ok_t
    Path('build/s4/corr_gauge.json').write_text(json.dumps(dict(
        norma_RefA_RT=float(nA), N_RefA=NA, N_RefAB=NAB, G=float(G), c_tau=float(c_tau), c_n=float(c_n),
        c_B=float(c_B), ZCorr_SSPACE=float(Zc), ZCorr_L=float(ZL), ZCorr_L_frac=str(ZL),
        dtauCorr_L=float(TL), dtauCorr_L_frac=str(TL), cobre_1e60=bool(ok)), indent=1) + '\n')


if __name__ == '__main__':
    main()
