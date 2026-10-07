#!/usr/bin/env python3
"""R4 -- verified_norms.cota_gram / cota_norma contra a verdade RIGOROSA (python-flint, Arb).

(mpmath nao esta instalado no .venv; python-flint 0.9 esta, e da bolas rigorosas, que e mais forte
que 50 digitos: prec = 200 bits.)

Verdade: para G hermitiana (entradas = floats exatos), com G ~ U Lambda U^* (eigh em float),
em Arb: F = U^* G U - Lambda, D = U^* U - I. Para x = U y:
   x^* G x <= (lam_max(Lambda) + ||F||_F) |y|^2,   x^* x >= (1 - ||D||_F) |y|^2,
logo lam_max(G) <= (lam_max(Lambda) + ||F||_F)/(1 - ||D||_F)   (se ||D||_F < 1).
Cota inferior: quociente de Rayleigh do autovetor dominante, em Arb.

Testes:
  A. cota_gram com Gram NAO hermitiano (um triangulo perturbado ~1e-15), quase singular
     (lam_min/lam_max ~ 1e-10), n = 10..400;
  B. topo aglomerado (metade dos autovalores iguais a lam_max): o caso dificil para Rump;
  C. falso positivo: chute tal que o primeiro t tenha t^2 = lam_max (1 - eps); se passar de
     primeira com t^2 < cota inferior rigorosa de lam_max, e ERRO;
  D. cota_norma(X) para X complexo ate 400 x 300 (valores singulares aglomerados/quase singular);
  E. escalas extremas (1e-280, 1e+140) para o termo de underflow.
"""
from __future__ import annotations

import math
from pathlib import Path
import sys
import time

import numpy as np
import flint

RAIZ = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RAIZ/'scripts'))
import verified_norms as vn

flint.ctx.prec = 200
arb, acb, acb_mat = flint.arb, flint.acb, flint.acb_mat


def para_acb(A):
    A = np.asarray(A, complex)
    return acb_mat([[acb(complex(x)) for x in linha] for linha in A])


def fro_sup(M):
    """cota superior (float) de ||M||_F para acb_mat M."""
    s = arb(0)
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            e = M[i, j]
            s += e.real*e.real + e.imag*e.imag
    s = s.sqrt()
    return float(s.mid() + s.rad())*(1 + 1e-15) + 1e-300


def lam_max_rigoroso(G_ex, Ga=None):
    """(inf, sup) rigorosos de lam_max(G_ex); G_ex hermitiana exata (np) ou Ga (acb_mat) dada."""
    Gf = np.asarray(G_ex, complex) if G_ex is not None else None
    if Ga is None:
        Ga = para_acb(Gf)
    if Gf is None:
        Gf = np.array([[complex(float(Ga[i, j].real.mid()), float(Ga[i, j].imag.mid()))
                        for j in range(Ga.ncols())] for i in range(Ga.nrows())])
    Gf = (Gf + Gf.conj().T)/2
    w, U = np.linalg.eigh(Gf)
    Ua = para_acb(U)
    UH = Ua.conjugate().transpose()
    n = len(w)
    L = acb_mat(n, n)
    for i in range(n):
        L[i, i] = acb(float(w[i]))
    F = UH*Ga*Ua - L
    I = acb_mat(n, n)
    for i in range(n):
        I[i, i] = acb(1)
    D = UH*Ua - I
    nF, nD = fro_sup(F), fro_sup(D)
    assert nD < 0.5
    sup = (float(w[-1]) + nF)/(1 - nD)*(1 + 1e-14)
    # inferior: Rayleigh do autovetor dominante, em Arb
    u = acb_mat([[acb(complex(x))] for x in U[:, -1]])
    uH = u.conjugate().transpose()
    num = (uH*Ga*u)[0, 0].real; den = (uH*u)[0, 0].real
    q = num/den
    inf = float(q.mid() - q.rad())*(1 - 1e-15)
    return inf, sup


def hermitiana_exata(A):
    """Hermitiana exata em float a partir do triangulo superior de A."""
    U = np.triu(A, 1)
    H = U + U.conj().T
    H[np.diag_indices(len(A))] = np.real(np.diag(A))
    return H


def espectro_G(n, lams, rng):
    Q, _ = np.linalg.qr(rng.standard_normal((n, n)) + 1j*rng.standard_normal((n, n)))
    return hermitiana_exata((Q*lams) @ Q.conj().T)


def perturba_triangulo(G, rel, rng):
    """G + perturbacao so no triangulo superior estrito (Gram nao hermitiano), relativa ~ rel."""
    P = np.triu(rng.standard_normal(G.shape) + 1j*rng.standard_normal(G.shape), 1)*rel*np.abs(G).max()
    Gc = G + P
    D = Gc - G                         # diferenca exata? cota rigorosa via Arb abaixo
    return Gc


def erro_rigoroso(Gc, G):
    """cota superior de ||Gc - G||_2 <= ||Gc - G||_F calculada em Arb."""
    return fro_sup(para_acb(Gc) - para_acb(G))


def teste_A_B(rng, log):
    log('== A/B: cota_gram com Gram nao hermitiano; quase singular / topo aglomerado')
    pior = 0.0; linhas = []
    for n in (10, 50, 120, 250, 400):
        for tipo in ('quase singular', 'topo aglomerado', 'topo 1e-13'):
            if tipo == 'quase singular':
                lams = np.logspace(0, -10, n)
            elif tipo == 'topo aglomerado':
                lams = np.concatenate([np.ones(n//2), np.logspace(-1, -10, n - n//2)])
            else:
                lams = np.concatenate([1 + 1e-13*rng.standard_normal(n//3)**2, np.logspace(-3, -10, n - n//3)])
            G = espectro_G(n, lams, rng)
            Gc = perturba_triangulo(G, 1e-15, rng)
            e = erro_rigoroso(Gc, G)
            t = vn.cota_gram(Gc, e)
            inf, sup = lam_max_rigoroso(G)
            ok = t*t >= sup
            pior = max(pior, math.sqrt(inf)/t)
            linhas.append((n, tipo, t, math.sqrt(inf), math.sqrt(sup), ok))
            log(f'   n={n:3d} {tipo:16s}: cota {t:.12f}  verdade em [{math.sqrt(inf):.12f}, {math.sqrt(sup):.12f}]'
                f'  cota/verdade-1 = {t/math.sqrt(sup) - 1:.2e}  {"OK" if ok else "*** VIOLACAO ***"}')
            assert ok
    return pior


def teste_C(rng, log):
    log('\n== C: falso positivo? primeiro t com t^2 = lam_max (1 - eps); passar de primeira seria ERRO')
    falsos = 0; n_testes = 0
    for n in (20, 100, 300):
        for tipo in ('espalhado', 'topo aglomerado'):
            lams = np.logspace(0, -6, n) if tipo == 'espalhado' else np.concatenate([np.ones(n//2), np.logspace(-1, -6, n - n//2)])
            G = espectro_G(n, lams, rng)
            inf, sup = lam_max_rigoroso(G)
            for eps in (1e-3, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 0.0):
                t0 = math.sqrt(inf*(1 - eps))
                chute = t0/(1 + 1e-6)
                t = vn.cota_gram(G, 0.0, chute=chute)
                t_ini = max(chute*(1 + 1e-6), 1e-150)
                passou_de_primeira = (t == t_ini)
                n_testes += 1
                if passou_de_primeira and t*t < inf:
                    falsos += 1
                    log(f'   *** FALSO POSITIVO: n={n} {tipo} eps={eps}: t^2 = {t*t!r} < inf = {inf!r}')
            log(f'   n={n:3d} {tipo:16s}: eps em 1e-3..0 -> nenhum aceite abaixo da verdade'
                if falsos == 0 else f'   n={n} {tipo}: {falsos} falsos')
    log(f'   {n_testes} tentativas, falsos positivos: {falsos}')
    assert falsos == 0
    return n_testes


def teste_D(rng, log):
    log('\n== D: cota_norma(X), X complexo (verdade: lam_max(X^* X) com X^* X em Arb)')
    pior = 0.0
    for (m, n, tipo) in ((60, 40, 'aglomerado'), (200, 120, 'quase singular'), (400, 300, 'aglomerado'),
                         (300, 400, 'aglomerado'), (150, 150, 'posto 1 + ruido')):
        k = min(m, n)
        if tipo == 'aglomerado':
            sv = np.concatenate([np.ones(k//2), np.logspace(-1, -5, k - k//2)])
        elif tipo == 'quase singular':
            sv = np.logspace(0, -5, k)
        else:
            sv = np.concatenate([[1.0], 1e-8*np.ones(k - 1)])
        Uq, _ = np.linalg.qr(rng.standard_normal((m, k)) + 1j*rng.standard_normal((m, k)))
        Vq, _ = np.linalg.qr(rng.standard_normal((n, k)) + 1j*rng.standard_normal((n, k)))
        X = (Uq*sv) @ Vq.conj().T
        t = vn.cota_norma(X)
        Xa = para_acb(X if n <= m else X.conj().T)
        Ga = Xa.conjugate().transpose()*Xa
        inf, sup = lam_max_rigoroso(None, Ga=Ga)
        ok = t*t >= sup
        pior = max(pior, math.sqrt(inf)/t)
        log(f'   {m}x{n} {tipo:15s}: cota {t:.12f}  verdade em [{math.sqrt(inf):.12f}, {math.sqrt(sup):.12f}]'
            f'  cota/verdade-1 = {t/math.sqrt(sup) - 1:.2e}  {"OK" if ok else "*** VIOLACAO ***"}')
        assert ok
    return pior


def teste_E(rng, log):
    log('\n== E: escalas extremas (underflow / grande)')
    for esc in (1e-280, 1e-150, 1e140):
        G = espectro_G(40, np.logspace(0, -10, 40), rng)*esc
        G = hermitiana_exata(G)
        t = vn.cota_gram(G, 0.0)
        inf, sup = lam_max_rigoroso(G)
        ok = t*t >= sup if sup > 0 else True
        log(f'   escala {esc:.0e}: cota^2/sup = {t*t/sup:.6f}  {"OK" if ok else "*** VIOLACAO ***"}')
        assert ok


if __name__ == '__main__':
    t0 = time.time()
    rng = np.random.default_rng(2026)
    pA = teste_A_B(rng, print)
    nC = teste_C(rng, print)
    pD = teste_D(rng, print)
    teste_E(rng, print)
    print(f'\nRESUMO R4: nenhuma cota abaixo da verdade rigorosa; max sqrt(inf)/cota = {max(pA, pD):.10f};'
          f' {nC} tentativas de falso positivo, nenhuma  [{time.time()-t0:.0f}s]')
