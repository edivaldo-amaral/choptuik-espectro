#!/usr/bin/env python3
"""R8 (revisao independente de S3b): conferencias geometricas e de espacos do fechamento (A7),
em racionais exatos.
 (i)   D = {|s - 0,401| < 1/8} esta no retangulo de S3a {0 <= Re s <= 1,765, |Im s| <= 1/4};
 (ii)  o disco de lambda* (|lambda* - lambda0| <= v0 r, C3.json) esta dentro do disco de Rouche e de D;
 (iii) o centro do Rouche e float(0,401) != 0,401: a diferenca e coberta pelos proprios tiles
       (cada tile prova K injetivo no disco inteiro dele; uniao contem a coroa
       |rho - |s - c|| <= rZ - sagita);
 (iv)  pesos: Y_L = X_(65/64) inter X_(64/65) (radial 5/4) esta contido em X_K = X_(129/128) (radial
       9/8) com norma <= 1; o espaco com sinal X_L NAO esta contido em X_K (m -> -infinito)."""
import json
from fractions import Fraction as F
from pathlib import Path


def main():
    c = F(401, 1000); R8 = F(1, 8)
    # (i)
    ok1 = (c - R8 >= 0) and (c + R8 <= F(1765, 1000)) and (R8 <= F(1, 4))
    print(f'(i)   D contido no retangulo de S3a: {ok1}  (Re em [{float(c-R8)}, {float(c+R8)}], |Im| < {float(R8)})')
    # (ii)
    c3 = json.loads(Path('build/s3b/rig/C3.json').read_text())
    lam0 = F(c3['lambda0']); e = F(c3['erro_lambda'])
    cf = F(0.401)                              # centro do Rouche em ponto flutuante
    ok2 = abs(lam0 - cf) + e < F(5, 100) and abs(lam0 - c) + e < R8
    print(f'(ii)  |lambda* - c| <= {float(abs(lam0 - cf) + e):.3e} < 0,05 e < 1/8: {ok2}')
    # (iii)
    js = json.loads(Path('build/s3b/rouche_K_circ.json').read_text())
    import flint
    flint.ctx.prec = 200
    rho, rZ = flint.arb(js['rho']), flint.arb(js['rZ']); n = js['n']
    dif = abs(cf - c)
    ok3 = True
    for t in (flint.arb(10)**-9, -flint.arb(10)**-9):
        for tile in js['tiles']:
            x0, y0 = (flint.arb(v) for v in tile['dV'])
            phi = flint.arb.atan2(y0, x0)
            for sgn in (1, -1):
                ang = phi + sgn*flint.arb.pi()/n
                x = (rho + t)*ang.cos() - x0; y = (rho + t)*ang.sin() - y0
                ok3 = ok3 and bool(x*x + y*y < rZ*rZ)
    print(f'(iii) |float(0,401) - 0,401| = {float(dif):.2e}; os 16 discos cobrem as circunferencias |s - c| = rho +- 1e-9'
          f' (pontas dos arcos, Arb; quase-convexidade como no codigo): {ok3}')
    # (iv)
    k1L, k2L, k1K, k2K = F(65, 64), F(5, 4), F(129, 128), F(9, 8)
    # Y_L: peso k1L^(2|m|); X_K: peso k1K^(2m). k1K^(2m) <= k1L^(2|m|) para todo m inteiro: m >= 0 porque
    # k1K <= k1L; m < 0 porque k1K^(2m) <= 1 <= k1L^(2|m|). Radial: omega(9/8) <= omega(5/4).
    ok4 = k1K <= k1L and k2K <= k2L
    print(f'(iv)  Y_L -> X_K com norma <= 1: {ok4};  X_L (com sinal) -> X_K: razao de pesos em m = -M vale '
          f'(130/129)^(2M) -> infinito (ex.: M = 200: {float((F(130, 129))**400):.2f}), logo X_L nao esta em X_K: '
          'o lema A3 e necessario para levar o autovetor de L ao espaco de K.')


if __name__ == '__main__':
    main()
