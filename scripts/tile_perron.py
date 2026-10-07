#!/usr/bin/env python3
"""Tile de S3a pelo criterio de Perron, na norma do toro com m COM SINAL.

Niveis (caixas Z c G+ c F; K_livre preserva caixas -- diagonal em m, triangular
superior em n -- logo e bloco-triangular superior na decomposicao [Z, T1, T2, far]):
  Z    = {|m| < MZ, n < NZ}: inverso aproximado V = H_ZZ^-1 (denso);
  T1   = G+ \\ Z: casca EXATA -- blocos finitos calculados inteiros (sem soma sobre modos);
  T2   = F \\ G+: composicoes radiais 1D modo a modo; T2 -> T2 por Schur sobre d;
         T2 -> Z e T2 -> T1 conjuntos (Gram acumulado);
  far  = complemento de F: ||B|| (simbolo, 1,65) vezes ||Q0 P_far||.

Um centro por nivel. Com P = Q0(sZ) P_Z + Q0(s1) P_T1 + Q0(s2) P_(T2+far) -- bloco-
triangular com diagonais invertiveis, logo uma bijecao sobre o dominio --

    H(s) = K(s) P = I + sum_l (B + s - s_l) Q0(s_l) P_l,

e um tile |s - sZ| <= rZ exige |s - s1| <= r1 e |s - s2| <= r2 nele. So as colunas de Z
sao sensiveis a s (||V Q0_ZZ|| ~ 18 perto do eixo imaginario); T1 e T2 saem uma vez
por tile medio/grande (cache em disco).

Com A = diag(V, I, I, I), E = I - A H; se a matriz N de cotas dos blocos de E tem raio
espectral < 1, existe v > 0 com N v < v, E contrai na norma max_i ||P_i x||/v_i e H e
injetivo, logo K(s) tambem.

MODOS. Sem --rigoroso: estimativa (iteracao de potencia), para escolher as caixas.
Com --rigoroso: cotas SUPERIORES verificadas (verified_norms) e contabilidade de erros:
  * pesos calculados: a norma com eles difere da do toro por <= 1 + 10u (FATOR_PESO);
  * fundo truncado na banda (BAND_M, BAND_N) para fora; o resto (Young), o erro de
    arredondamento das entradas (gamma_8 Young_abs) e o erro de RT entram em dB;
    mu R_x e exato (mu e diadico);
  * R_m^-1: pelo residuo E = I - J R com pesos, ||R_exato - R|| <= ||R|| ||E||/(1-||E||);
  * produtos: ||fl(AB) - AB|| <= gamma_k ||A||_F ||B||_F; Grams acumulados idem;
  * Perron: v exibido e N v <= theta v verificado em racionais exatos.

Registro (s0 = 0, F = 200x400): sem casca, Z = 12x36 da 2,22, 20x80 da 1,27, 30x120
da 0,950; com casca 40x160 e Z = 20x80, 0,866. Um nivel por modo, ou faixas de modos,
piora: as linhas de Z e de far viram somas sobre modos.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import math
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fourier_tail_bound as ft
from signed_operator import radial, Rinv, Jrad, modos, MU, K1, K2
from tile_certificate import norma2
from verified_norms import cota_norma, cota_gram, gamma, U

NORMA_B = 1.64967          # ||B#|| na norma do toro, simbolo certificado (1024x2048)
BAND_M, BAND_N = 20, 40    # banda do fundo truncado para fora (suporte completo: 40, 100)
FATOR_PESO = 1 + 1e-12     # pesos calculados com erro relativo <= 4u
EPS_FUNDO = 40*2.0**-25    # erro do fundo de RT na norma do toro (precedente: simbolo)
RIG = False


class Modo:
    """Indices radiais de um modo truncado em n < N."""
    def __init__(self, N):
        self.n1, self.n2 = radial(N)
        self.n = np.concatenate([self.n1, self.n2])
        self.w = np.sqrt(np.concatenate([ft.omega(self.n1), ft.omega(self.n2)]))

    def faixa(self, lo, hi):
        return np.where((self.n >= lo) & (self.n < hi))[0]


def Bw(g, d, ent: Modo, sai: Modo, mu=MU):
    """Bloco radial pesado kappa1^d W B_d W^-1 de ent para sai."""
    X = ft.B_rad(g, d, ent.n1, ent.n2, sai.n1, sai.n2, mu)
    return K1**d*(sai.w[:, None]*X)/ent.w[None, :]


def Rw(m, s0, md: Modo, N):
    return (md.w[:, None]*Rinv(m, s0, N))/md.w[None, :]


def lmax(G):
    return math.sqrt(max(float(np.linalg.eigvalsh(G)[-1]), 0.0)) if G.size else 0.0


# ---------- normas: estimativa ou cota verificada
def fro(X):
    return float(np.linalg.norm(X)) if X.size else 0.0


def err_mm(A, B):
    """||fl(A B) - A B||_2 (inclui o arredondamento dos pesos em A, relativo <= 8u)."""
    return gamma(A.shape[1] + 8)*fro(A)*fro(B) if RIG else 0.0


def N2(X, erro=0.0):
    return cota_norma(X, erro)*FATOR_PESO if RIG else norma2(X)


def L2(G, erro=0.0):
    return cota_gram(G, erro)*FATOR_PESO if RIG else lmax(G)


def err_gram(blocos_fro2, dim_interna):
    return gamma(dim_interna)*blocos_fro2 if RIG else 0.0


def Rw_cert(m, s0, md: Modo, N):
    """R_m^-1 pesado e cota de ||R_exato - R_calculado|| (0 na estimativa)."""
    R = Rw(m, s0, md, N)
    if not RIG:
        return R, 0.0
    Jw = (md.w[:, None]*Jrad(m, s0, N))/md.w[None, :]
    E = np.eye(len(R)) - Jw @ R
    e = cota_norma(E, err_mm(Jw, R) + U*fro(E))
    assert e < 0.5, (m, s0, e)
    return R, cota_norma(R)*e/(1 - e)


def young_abs(g, d):
    """Young do fundo com campos em valor absoluto (majora o erro de arredondamento)."""
    v = [np.abs(ft.modo(g, j, d)) for j in range(3)]
    n = np.arange(len(v[0]))
    l1 = lambda a: float(np.sum(a*np.where(n == 0, 1.0, 2*K2**n)))
    a1, a2, a3 = (l1(x) for x in v)
    M = np.array([[0.5*(2*a2 + 2*a3), 0.5*(a1 + 2*a3), 0.5*(2*a2 + 2*a3)],
                  [0.5*(a1 + 2*a2), 0.5*(a1 + 4*a2), 0.5*(a1 + 4*a2)]])
    return float(np.linalg.norm(M, 2))*1.0001


def truncar(g, Dm, Dn):
    """Fundo truncado: modos de Fourier |d| <= Dm e coeficientes radiais n <= Dn. A banda
    radial so limita a subida (Toeplitz); Hankel e mu R_x so descem."""
    out = []
    for campo in g:
        h = {}
        for m, a in campo.items():
            if m <= Dm:
                b = a.copy(); b[Dn + 1:] = 0
                h[m] = b
        out.append(h)
    return out


def resto_banda(g, Dm, Dn):
    """Cota de Young de ||B - B_truncado|| na norma do toro (estrutura de norma_Bd_young)."""
    n = np.arange(101)
    wn = np.where(n == 0, 1.0, 2*K2**n)
    def young_parte(d, nmin):
        F1, F2, F3 = ft.combos(g, d)
        def l1(a):
            b = np.abs(a).copy(); b[:nmin] = 0
            return float(np.sum(b*wn))
        M = np.array([[0.5*l1(F1[0]), 0.5*l1(F2[0]), 0.5*l1(F3[0])],
                      [0.5*l1(F1[1]), 0.5*l1(F2[1]), 0.5*l1(F3[1])]])
        return float(np.linalg.norm(M, 2))
    r = sum(K1**abs(d)*young_parte(d, 0) for d in range(-40, 41, 2) if abs(d) > Dm)
    r += sum(K1**abs(d)*young_parte(d, Dn + 1) for d in range(-Dm, Dm + 1, 2))
    return r*1.0001


def descida(J, Nc, mu=MU):
    """Cota de ||P_{n < J} R_m^-1 P_{n >= Nc}|| (pesada, uniforme em m e em Re s >= 0) pela
    forma fechada: |entrada| <= (2/(mu(n_j + 1))) c_w kappa2^(n_j - n_c). Teste de Schur:
    colunas <= (2 c_w/mu) H_J kappa2^(J - Nc), linhas <= (2 c_w/mu) kappa2^(J - Nc) kappa2/(kappa2 - 1)."""
    cw = math.sqrt(1 + K2**-4)
    HJ = sum(1/(j + 1) for j in range(J))
    base = (2*cw/mu)*K2**(J - Nc)
    return math.sqrt(base*HJ*base*K2/(K2 - 1))*1.0001


def descida_rx(J, Nc):
    """Cota de ||P_{n < J} mu R_x Q0 P_{n >= Nc}|| (pesada, uniforme em m e em Re s >= 0).

    mu R_x (divisao por xi, d = 0) leva a coluna c1 de n impar a TODAS as linhas c2 pares
    2j < n: desce de qualquer n (a banda truncada so limita a subida). Entradas pesadas:
    |mu R_x (2j <- n)| <= 2 mu c_w kappa2^(2j - n); |Q0 (n <- c)| <= (2 c_w/(mu(n+1))) kappa2^(n - c).
    Composicao (2j <- c): <= 4 c_w^2 (1 + ln(c+1)) kappa2^(2j - c). Schur com f(c) = (1 + ln(c+1))
    kappa2^-c decrescente, razao <= q = (1 + 1/((Nc+1)(1 + ln(Nc+1))))/kappa2 < 1:
    linhas <= 4 c_w^2 kappa2^J f(Nc)/(1 - q), colunas <= 4 c_w^2 f(Nc) kappa2^J/(1 - kappa2^-2)."""
    cw = math.sqrt(1 + K2**-4)
    f = (1 + math.log(Nc + 1))*K2**(-Nc)
    q = (1 + 1/((Nc + 1)*(1 + math.log(Nc + 1))))/K2
    assert q < 1
    base = 4*cw*cw*f*K2**J
    return math.sqrt(base/(1 - q)*base/(1 - K2**-2))*1.0001


# mu verdadeiro de RT: |mu - MU_80| <= 1e-80 (RT, eq. dnkdfkjhfkdhfkdhdfdk2), com 80 digitos
MU_80 = Fraction('0.16830707896344996951013497904285742072100199080892966476395293134873313662587505')
DELTA_MU = float(abs(Fraction(722873400, 2**32) - MU_80) + Fraction(1, 10**80))*1.0001   # ~3,88e-11


def eps_mu(s):
    """Cota de ||delta_mu (J1 + R_x) Q0(s)|| (pesada), com K(mu) = K(mu_RefA) + delta_mu (J1 + R_x),
    J1 = dK_livre/dmu (diag(n+1) + 2n acima da diagonal), |delta_mu| <= DELTA_MU. Modo a modo,
    R_m = MU J1 + sigma_m (sigma_m = s + i m/2), logo J1 R_m^-1 = (I - sigma_m R_m^-1)/MU.
    Forma fechada com |d_k| >= max(MU, |m|/2 - |Im s|) (Re s >= 0): ||R_m^-1|| <= c/max(...),
    c = 1 + 2 c_w/(kappa2 - 1); ||R_x|| <= 2 c_w/(kappa2 (1 - kappa2^-2)) (Schur)."""
    cw = math.sqrt(1 + K2**-4)
    c = 1 + 2*cw/(K2 - 1)
    nrx = 2*cw/(K2*(1 - K2**-2))
    pior = 0.0
    for m in list(range(0, 4002, 2)):
        h = max(MU, m/2 - abs(s.imag))
        nR = c/h
        sig = abs(s.real) + m/2 + abs(s.imag)
        pior = max(pior, (1 + sig*nR)/MU + nrx*nR)
    # |m| > 4000: (1 + sig nR)/MU decresce para (1 + c)/MU e nrx nR para 0 (ja dominado acima)
    pior = max(pior, (1 + c*(1 + 1e-3))/MU)
    return DELTA_MU*pior*1.0001


def verifica_perron(N, theta):
    """Exibe v > 0 com N v <= theta v, conferido em racionais exatos. Devolve (ok, v, max razao)."""
    M = np.array(N, float)
    lam, vec = np.linalg.eig(M + 1e-9)
    v = np.abs(vec[:, np.argmax(np.abs(lam))]); v = v/v.max()
    v = np.maximum(v, 1e-6)
    Nq = [[Fraction(x) for x in linha] for linha in M.tolist()]
    vq = [Fraction(x) for x in v.tolist()]
    razoes = [sum(Nq[i][j]*vq[j] for j in range(len(vq)))/vq[i] for i in range(len(vq))]
    rmax = max(razoes)
    return rmax <= Fraction(theta), v, float(rmax)


class Ctx:
    """Caixas, fundo truncado, perturbacao global dB e indices de Z e T1 dentro de G+."""
    def __init__(self, Z, Gp, F, D):
        (MZ, NZ), (MP, NP), (Mc, Nc) = self.Z, self.Gp, self.F = Z, Gp, F
        assert MZ <= MP <= Mc and NZ <= NP <= Nc
        assert Mc >= MP + BAND_M and Nc >= NP + BAND_N + 1, 'F precisa conter G+ mais a banda do fundo'
        self.D = D
        self.g0 = ft.modos_campos()
        self.g = truncar(self.g0, BAND_M, BAND_N)
        self.young = {d: ft.norma_Bd_young(self.g, d, MU) for d in range(-BAND_M, BAND_M + 1, 2)}
        self.resto_banda = resto_banda(self.g0, BAND_M, BAND_N)
        self.dB = ((sum(K1**abs(d)*gamma(8)*young_abs(self.g, d) for d in self.young) if RIG else 0.0)
                   + EPS_FUNDO + self.resto_banda)
        self.nB = NORMA_B + self.dB
        self.chave = (tuple(Z), tuple(Gp), tuple(F), D, BAND_M, BAND_N, RIG)
        self.mP = Modo(NP); self.msP = modos(MP); self.k = len(self.mP.n)
        self.iP = {m: i for i, m in enumerate(self.msP)}
        self.mZ = Modo(NZ); self.msZ = modos(MZ); self.kZ = len(self.mZ.n)
        zloc = self.mP.faixa(0, NZ)
        self.iZ = np.concatenate([self.iP[m]*self.k + zloc for m in self.msZ])
        em_Z = np.zeros(len(self.msP)*self.k, bool); em_Z[self.iZ] = True
        self.i1 = np.where(~em_Z)[0]
        # colunas de T1 por modo (indices locais do modo em G+)
        self.cols1 = {m: (self.mP.faixa(NZ, NP) if abs(m) < MZ else np.arange(self.k)) for m in self.msP}
        self.t0 = time.time()

    def tempo(self):
        return f'[{time.time()-self.t0:.0f}s]'


def linhas_fora(ctx, r_, lo_n):
    """Linhas do modo r_ fora de G+ alcancaveis pela banda, na base Modo(NP + BAND_N + 1)."""
    (MP, NP) = ctx.Gp
    m2 = Modo(NP + BAND_N + 1)
    return m2, (m2.faixa(NP, NP + BAND_N + 1) if abs(r_) < MP else m2.faixa(lo_n, NP + BAND_N + 1))


def parte_Z(ctx, sZ, log=print, dV=0.0):
    """Colunas de Z com o precondicionador Q0(sZ): V, rho, ||V Q_ZZ||, b1 (T1 <- Z), b2 (T2 <- Z).
    V inverte H_ZZ(sZ + dV) = I + (B + dV) Q_ZZ(sZ) (dV = 0 em S3a; no Rouche de S3b o
    precondicionador fica no centro do disco e V anda pela circunferencia). b1, b2 e ||Q_ZZ||
    nao dependem de dV: Q0 P_Z cai em Z, logo P_T (s - sZ) Q0 P_Z = 0."""
    (MZ, NZ), (MP, NP) = ctx.Z, ctx.Gp
    g, mZ, msZ, kZ = ctx.g, ctx.mZ, ctx.msZ, ctx.kZ
    QC = [Rw_cert(m, sZ, mZ, NZ) for m in msZ]
    Q = [x[0] for x in QC]; dQ = max(x[1] for x in QC)
    nQ = max(N2(x) for x in Q) + dQ
    q = {'nQZ': nQ}
    # H_ZZ
    H = np.eye(len(msZ)*kZ, dtype=complex)
    Bc = {}; e2 = 0.0
    for j, mi in enumerate(msZ):
        for i, mo in enumerate(msZ):
            d = mo - mi
            if abs(d) <= BAND_M:
                if d not in Bc:
                    Bc[d] = Bw(g, d, mZ, mZ)
                H[i*kZ:(i+1)*kZ, j*kZ:(j+1)*kZ] += Bc[d] @ Q[j]
                e2 += err_mm(Bc[d], Q[j])**2
    dH = math.sqrt(e2) + (U*fro(H) if RIG else 0.0) + ctx.dB*nQ + ctx.nB*dQ
    if dV != 0:
        # dV Q_ZZ nos blocos diagonais: produto complexo (gamma(4) relativo), mais uma soma
        for i in range(len(msZ)):
            H[i*kZ:(i+1)*kZ, i*kZ:(i+1)*kZ] += dV*Q[i]
        dH += abs(dV)*(dQ + (gamma(4)*math.sqrt(sum(fro(x)**2 for x in Q)) if RIG else 0.0)) \
            + (U*fro(H) if RIG else 0.0)
    V = np.linalg.inv(H)
    q['V'] = nV = N2(V)
    X = np.eye(len(H)) - V @ H
    q['rho'] = N2(X, err_mm(V, H) + (U*fro(X) if RIG else 0.0) + nV*dH)
    QZ = np.zeros_like(H)
    for i in range(len(msZ)):
        QZ[i*kZ:(i+1)*kZ, i*kZ:(i+1)*kZ] = Q[i]
    q['VQ_ZZ'] = N2(V @ QZ, err_mm(V, QZ) + nV*dQ)
    del QZ, X
    # b1: linhas de T1 (dentro de G+) a partir de Z; b2: linhas fora de G+
    mP = ctx.mP
    G1 = np.zeros((len(H), len(H)), complex); f1 = 0.0; l1 = 0; e1 = 0.0
    G2 = np.zeros((len(H), len(H)), complex); f2 = 0.0; l2 = 0; e2b = 0.0
    BcP = {}; Bc2 = {}
    for mo in modos(MZ + BAND_M):
        # dentro de G+ (T1)
        linhas = mP.faixa(NZ, NP) if abs(mo) < MZ else mP.faixa(0, NP)
        linhas = linhas[mP.n[linhas] < NZ + BAND_N + 1]
        if abs(mo) < MP and len(linhas):
            X = np.zeros((len(linhas), len(H)), complex)
            for j, mi in enumerate(msZ):
                d = mo - mi
                if abs(d) <= BAND_M:
                    if d not in BcP:
                        BcP[d] = Bw(g, d, mZ, mP)
                    A = BcP[d][linhas]
                    X[:, j*kZ:(j+1)*kZ] = A @ Q[j]
                    e1 += err_mm(A, Q[j])**2
            G1 += X.conj().T @ X; f1 += fro(X)**2; l1 += len(linhas)
        # fora de G+ (T2)
        m2, lf = linhas_fora(ctx, mo, 0)
        lf = lf[m2.n[lf] < NZ + BAND_N + 1]
        if len(lf):
            X = np.zeros((len(lf), len(H)), complex)
            for j, mi in enumerate(msZ):
                d = mo - mi
                if abs(d) <= BAND_M:
                    if d not in Bc2:
                        Bc2[d] = Bw(g, d, mZ, m2)
                    A = Bc2[d][lf]
                    X[:, j*kZ:(j+1)*kZ] = A @ Q[j]
                    e2b += err_mm(A, Q[j])**2
            G2 += X.conj().T @ X; f2 += fro(X)**2; l2 += len(lf)
    extra = ctx.dB*nQ + ctx.nB*dQ
    q['b1'] = L2(G1, err_gram(f1, max(l1, 1))) + math.sqrt(e1) + extra
    q['b2'] = L2(G2, err_gram(f2, max(l2, 1))) + math.sqrt(e2b) + extra
    log(f'  Z ({len(H)}) em sZ = {sZ}' + (f' + {dV:.6f}' if dV != 0 else '') + f': ||V|| {nV:.3f}, rho {q["rho"]:.1e}, ||V Q_ZZ|| {q["VQ_ZZ"]:.2f},'
        f' b1 {q["b1"]:.4f}, b2 {q["b2"]:.4f}, erro de H_ZZ {dH:.1e}  {ctx.tempo()}', flush=True)
    return q, V


def parte_T1(ctx, s1, log=print):
    """Colunas de T1 no centro s1: alfa11, ||Q_11||, c21 (T2 <- T1) e as matrizes H_Z,T1 e
    Q_Z,T1 (esta so nas colunas radiais dos modos de Z), que V multiplica em cada tile."""
    (MZ, NZ), (MP, NP) = ctx.Z, ctx.Gp
    g, mP, msP, k, iP = ctx.g, ctx.mP, ctx.msP, ctx.k, ctx.iP
    QC = {m: Rw_cert(m, s1, mP, NP) for m in msP}
    dQ = max(x[1] for x in QC.values())
    nQ = max(N2(QC[m][0][:, ctx.cols1[m]]) for m in msP) + dQ
    c0 = np.cumsum([0] + [len(ctx.cols1[m]) for m in msP])
    n1 = int(c0[-1])
    assert n1 == len(ctx.i1)
    # H_(G+),T1 = P_T1 + P_(G+) B Q0 P_T1
    HG = np.zeros((len(msP)*k, n1), complex)
    Bc = {}; e2 = 0.0
    for j, mi in enumerate(msP):
        RT = QC[mi][0][:, ctx.cols1[mi]]
        for i, mo in enumerate(msP):
            d = mo - mi
            if abs(d) <= BAND_M:
                if d not in Bc:
                    Bc[d] = Bw(g, d, mP, mP)
                HG[i*k:(i+1)*k, c0[j]:c0[j+1]] += Bc[d] @ RT
                e2 += err_mm(Bc[d], RT)**2
        # identidade de T1: linha global iP[mi]*k + cols1 local
        HG[iP[mi]*k + ctx.cols1[mi], np.arange(c0[j], c0[j+1])] += 1.0
    del Bc
    dH = math.sqrt(e2) + (U*fro(HG) if RIG else 0.0) + ctx.dB*nQ + ctx.nB*dQ
    q = {'nQ1': nQ, 'dH1': dH, 'dQ1': dQ}
    HZ1 = HG[ctx.iZ].copy()
    X = HG[ctx.i1]                           # copia; H inteiro sai da memoria antes da norma
    del HG
    X[np.arange(n1), np.arange(n1)] -= 1.0
    q['alfa11'] = N2(X, dH + (U*fro(X) if RIG else 0.0))
    del X
    # Q0 P_T1: bloco T1 x T1 e bloco Z x (colunas radiais dos modos de Z)
    QG = np.zeros((len(msP)*k, n1), complex)
    for j, mi in enumerate(msP):
        QG[iP[mi]*k:(iP[mi]+1)*k, c0[j]:c0[j+1]] = QC[mi][0][:, ctx.cols1[mi]]
    q['Q_11'] = N2(QG[ctx.i1], dQ)
    colsZ = np.concatenate([np.arange(c0[j], c0[j+1]) for j, mi in enumerate(msP) if abs(mi) < MZ])
    QZ1 = QG[ctx.iZ][:, colsZ].copy()
    del QG
    # c21: linhas fora de G+
    G1 = np.zeros((n1, n1), complex); f1 = 0.0; lt = 0; e3 = 0.0
    Bc = {}
    for r_ in modos(MP + BAND_M):
        m2, lf = linhas_fora(ctx, r_, 0)
        if not len(lf):
            continue
        X = np.zeros((len(lf), n1), complex)
        for j, mi in enumerate(msP):
            d = r_ - mi
            if abs(d) <= BAND_M:
                if d not in Bc:
                    Bc[d] = Bw(g, d, mP, m2)
                A = Bc[d][lf]
                RT = QC[mi][0][:, ctx.cols1[mi]]
                X[:, c0[j]:c0[j+1]] = A @ RT
                e3 += err_mm(A, RT)**2
        G1 += X.conj().T @ X; f1 += fro(X)**2; lt += len(lf)
    q['c21'] = L2(G1, err_gram(f1, lt)) + math.sqrt(e3) + ctx.dB*nQ + ctx.nB*dQ
    del G1
    log(f'  T1 ({n1}) em s1 = {s1}: alfa11 {q["alfa11"]:.4f}, ||Q_11|| {q["Q_11"]:.3f},'
        f' c21 {q["c21"]:.4f}, erro de H_.,T1 {dH:.1e}  {ctx.tempo()}', flush=True)
    return q, HZ1, QZ1


def parte_T2(ctx, s2, log=print):
    """Colunas de T2 e far no centro s2. Devolve escalares e os Grams S_Z, S_QZ (linhas de Z),
    que o sanduiche com V de cada tile usa."""
    (MZ, NZ), (MP, NP), (Mc, Nc) = ctx.Z, ctx.Gp, ctx.F
    g, msP, iP, k, mP, iZ, i1, young, D = ctx.g, ctx.msP, ctx.iP, ctx.k, ctx.mP, ctx.iZ, ctx.i1, ctx.young, ctx.D
    q = {}
    mc = Modo(Nc); No = Nc + BAND_N + 1; mo_ = Modo(No)
    SZ = np.zeros((len(iZ), len(iZ)), complex)
    S1 = np.zeros((len(i1), len(i1)), complex)
    SQZ = np.zeros((len(iZ), len(iZ)), complex)
    SQ1 = np.zeros((len(i1), len(i1)), complex)
    acc = dict(eY=0.0, fYZ=0.0, fY1=0.0, fyZ=0.0, fy1=0.0, cols=0, dQ2=0.0)
    Bc = {}; BcP = {}
    sup22 = {}; supf = {}; arg22 = {}
    QT2 = 0.0; dQ2max = 0.0
    for m in modos(Mc):
        R, dq = Rw_cert(m, s2, mc, Nc)
        cols = mc.faixa(NP, Nc) if abs(m) < MP else mc.faixa(0, Nc)
        RT = R[:, cols]
        QT2 = max(QT2, N2(RT, dq)); dQ2max = max(dQ2max, dq)
        acc['cols'] += len(cols)
        if abs(m) < MP:                    # ds Q0: a descida radial cai em G+
            y = np.zeros((len(msP)*k, len(cols)), complex)
            y[iP[m]*k:(iP[m]+1)*k] = RT[mc.faixa(0, NP)]
            yZ = y[iZ]; y1 = y[i1]
            SQZ += yZ @ yZ.conj().T; SQ1 += y1 @ y1.conj().T
            acc['fyZ'] += fro(yZ)**2; acc['fy1'] += fro(y1)**2; acc['dQ2'] += dq*dq
        if abs(m) < MP + BAND_M:           # B Q0 de volta para G+ (Z e T1), conjunto
            Y = np.zeros((len(msP)*k, len(cols)), complex)
            for mo in msP:
                d = mo - m
                if abs(d) <= BAND_M:
                    if d not in BcP:
                        BcP[d] = Bw(g, d, mc, mP)
                    Y[iP[mo]*k:(iP[mo]+1)*k] = BcP[d] @ RT
                    acc['eY'] += err_mm(BcP[d], RT)**2
            YZ = Y[iZ]; Y1 = Y[i1]
            SZ += YZ @ YZ.conj().T; S1 += Y1 @ Y1.conj().T
            acc['fYZ'] += fro(YZ)**2; acc['fY1'] += fro(Y1)**2
        for d in range(-D, D + 1, 2):
            rr = m + d
            if d not in Bc:
                Bc[d] = Bw(g, d, mc, mo_)
            X = Bc[d] @ RT
            # ||B_d|| <= ||B|| (a extracao do modo d e media de conjugacoes unitarias por rotacao
            # em theta); norma_Bd_young nao serve para d = 0: o termo de mu R_x nao e cota provada
            erro = err_mm(Bc[d], RT) + ctx.nB*dq
            if abs(rr) >= Mc:
                l2 = np.array([], int); lf = np.arange(len(mo_.n))
            elif abs(rr) < MP:
                l2 = mo_.faixa(NP, Nc); lf = mo_.faixa(Nc, No)
            else:
                l2 = mo_.faixa(0, Nc); lf = mo_.faixa(Nc, No)
            if len(l2):
                v = N2(X[l2], erro)
                if v > sup22.get(d, 0.0):
                    sup22[d] = v; arg22[d] = m
            if len(lf):
                supf[d] = max(supf.get(d, 0.0), N2(X[lf], erro))
        ft._RCACHE.clear()
    del Bc, BcP
    resto = sum(K1**abs(d)*young[d] for d in young if abs(d) > D)*QT2
    extra = math.sqrt(acc['eY']) + ctx.dB*QT2 + ctx.nB*dQ2max
    q.update(cols=acc['cols'], fYZ=acc['fYZ'], fyZ=acc['fyZ'], extra=extra, dQ2=math.sqrt(acc['dQ2']))
    q['c12'] = L2(S1, err_gram(acc['fY1'], acc['cols'])) + extra
    q['Q_12'] = L2(SQ1, err_gram(acc['fy1'], acc['cols'])) + math.sqrt(acc['dQ2'])
    q['alfa22'] = sum(sup22.values()) + resto + ctx.dB*QT2
    q['alfa_f2'] = sum(supf.values()) + resto + ctx.dB*QT2
    q['Q_22'] = QT2
    # far
    rad = max(max(ft.cauda_Rinv(s2 + 1j*m/2, Nc, 1, MU), ft.cauda_Rinv(s2 + 1j*m/2, Nc, 2, MU))
              for m in modos(Mc))
    fou = ft.Rinv_uniforme(Mc, MU, s2)*1.0001      # formula fechada avaliada em ponto flutuante
    q['Qfar'] = max(rad, fou)
    q['beta'] = ctx.nB*q['Qfar']
    log(f'  T2 em s2 = {s2}: T1 <- T2 {q["c12"]:.4f}; alfa22 {q["alfa22"]:.4f} (|d|>{D}: {resto:.1e});'
        f' far <- T2 {q["alfa_f2"]:.4f}; ||Q0 P_T2|| {QT2:.3f}; ||Q0 P_far|| <= max({rad:.4f}, {fou:.4f});'
        f' beta {q["beta"]:.4f}  {ctx.tempo()}')
    log('     alfa22 por d: ' + ' '.join(f'{d:+d}:{v:.3f}@{arg22[d]}' for d, v in sorted(sup22.items())))
    return q, SZ, SQZ


def com_cache(caminho, s, ctx, calcula, nomes):
    """Le a parte de um nivel do cache (.npz) se for do mesmo centro e caixas; senao calcula."""
    if caminho is not None and Path(caminho).exists():
        dados = np.load(caminho, allow_pickle=True)
        q = dados['q'].item()
        if q['_s'] == complex(s) and q['_chave'] == ctx.chave:
            return (q,) + tuple(dados[n] for n in nomes)
        # cache de outro centro ou de outras caixas: recalcula e sobrescreve
    out = calcula()
    q = out[0]; q['_s'] = complex(s); q['_chave'] = ctx.chave
    if caminho is not None:
        np.savez(caminho, q=np.array(q, dtype=object), **{n: a for n, a in zip(nomes, out[1:])})
    return out


def montar(ctx, qZ, V, q1, HZ1, QZ1, q2, SZ, SQZ, rZ, r1, r2, log=print):
    """Matriz de Perron N0 + rZ NZ + r1 N1 + r2 N2 nos niveis (Z, T1, T2, far)."""
    (MZ, NZ), (MP, NP), (Mc, Nc) = ctx.Z, ctx.Gp, ctx.F
    nV = qZ['V']; dB, nB = ctx.dB, ctx.nB
    a1 = N2(V @ HZ1, err_mm(V, HZ1) + nV*q1['dH1'])
    VQ_Z1 = N2(V @ QZ1, err_mm(V, QZ1) + nV*q1['dQ1'])
    def sanduiche(S, fS2, extra_op):
        W = V @ S
        G = W @ V.conj().T
        eG = (err_mm(V, S)*fro(V) + err_mm(W, V.conj().T) + nV*nV*err_gram(fS2, q2['cols'])) if RIG else 0.0
        return L2(G, eG) + nV*extra_op
    a2 = sanduiche(SZ, q2['fYZ'], q2['extra'])
    VQ_Z2 = sanduiche(SQZ, q2['fyZ'], q2['dQ2'])
    dZ, dZ0 = descida(NZ + BAND_N + 1, Nc), descida(NZ, Nc)
    d1, d10 = descida(NP + BAND_N + 1, Nc), descida(NP, Nc)
    Qf = q2['Qfar']
    # far -> G+: o fundo truncado so SOBE BAND_N em n, mas mu R_x (d = 0) DESCE de qualquer n
    # (revisao independente, 24/09); entao P_G+ B Q0 P_far = [banda: descida(J + BAND_N + 1)]
    # + [mu R_x: descida_rx(J)], mais a perturbacao global dB (resto, erro de RT, arredondamento)
    rxZ, rx1 = descida_rx(NZ, Nc), descida_rx(NP, Nc)
    N0 = np.array([
        [qZ['rho'], a1, a2, nV*(nB*dZ + rxZ + dB*Qf)],
        [qZ['b1'], q1['alfa11'], q2['c12'], nB*d1 + rx1 + dB*Qf],
        [qZ['b2'], q1['c21'], q2['alfa22'], q2['beta']],
        [dB*qZ['nQZ'], dB*q1['nQ1'], q2['alfa_f2'], q2['beta']]])
    NZm = np.zeros((4, 4)); NZm[0, 0] = qZ['VQ_ZZ']
    N1m = np.zeros((4, 4)); N1m[0, 1] = VQ_Z1; N1m[1, 1] = q1['Q_11']
    N2m = np.zeros((4, 4)); N2m[:, 2] = [VQ_Z2, q2['Q_12'], q2['Q_22'], 0.0]
    N2m[:, 3] = [nV*dZ0, d10, Qf, Qf]
    # cada entrada e uma soma/produto de poucas dezenas de cotas positivas em ponto flutuante,
    # que pode arredondar para baixo por ~1e-14 relativo: inflar cobre isso
    ARRED = 1 + 1e-9
    N0, NZm, N1m, N2m = N0*ARRED, NZm*ARRED, N1m*ARRED, N2m*ARRED
    return N0, NZm, N1m, N2m, dict(a1=a1, VQ_Z1=VQ_Z1, a2=a2, VQ_Z2=VQ_Z2)


def tile(sZ, rZ, Z, Gp, F, D=12, s1=None, r1=None, s2=None, r2=None, cache1=None, cache2=None,
         log=print, ctx=None, alvo=0.95):
    s1 = sZ if s1 is None else s1; r1 = rZ if r1 is None else r1
    s2 = sZ if s2 is None else s2; r2 = rZ if r2 is None else r2
    assert abs(sZ - s1) + rZ <= r1 + 1e-15 and abs(sZ - s2) + rZ <= r2 + 1e-15, \
        'o tile precisa caber nos discos dos centros de T1 e T2'
    if RIG:
        assert min(sZ.real, s1.real, s2.real) >= 0, 'as cotas da forma fechada pedem Re s >= 0'
    if ctx is None:
        ctx = Ctx(Z, Gp, F, D)
    ctx.t0 = time.time()
    log(f'  fundo: resto da banda {ctx.resto_banda:.1e}, dB {ctx.dB:.1e}')
    qZ, V = parte_Z(ctx, sZ, log)
    q1, HZ1, QZ1 = com_cache(cache1, s1, ctx, lambda: parte_T1(ctx, s1, log), ('HZ1', 'QZ1'))
    q2, SZ, SQZ = com_cache(cache2, s2, ctx, lambda: parte_T2(ctx, s2, log), ('SZ', 'SQZ'))
    N0, NZm, N1m, N2m, qa = montar(ctx, qZ, V, q1, HZ1, QZ1, q2, SZ, SQZ, rZ, r1, r2, log)
    # mu verdadeiro de RT: |mu - mu_RefA| <= DELTA_MU ~ 3,9e-11; a coluna do nivel l ganha eps_mu(s_l),
    # vezes ||V|| na linha de Z (revisao independente, 24/09)
    for j, c in enumerate((sZ, s1, s2, s2)):
        e = eps_mu(c)*(1 + 1e-9)
        N0[0, j] += qZ['V']*e
        N0[1:, j] += e
    raio = lambda M: float(max(abs(np.linalg.eigvals(M))))
    um_centro = (s1 == sZ and s2 == sZ and r1 == rZ and r2 == rZ)
    Nde = (lambda x: N0 + x*(NZm + N1m + N2m)) if um_centro else (lambda x: N0 + x*NZm + r1*N1m + r2*N2m)
    Nfinal = Nde(rZ)
    q = dict(rigoroso=RIG, banda=(BAND_M, BAND_N), dB=ctx.dB, Z=qZ, T1=q1, T2=q2, montagem=qa,
             N0=N0.tolist(), NZ=NZm.tolist(), N1=N1m.tolist(), N2=N2m.tolist())
    q['perron0'] = raio(N0); q['perron'] = raio(Nfinal)
    lo = 0.0
    hi = 1.0 if um_centro else max(0.0, min(r1 - abs(sZ - s1), r2 - abs(sZ - s2)))
    if raio(Nde(0.0)) < alvo:
        for _ in range(40):
            mid = (lo + hi)/2
            lo, hi = (mid, hi) if raio(Nde(mid)) < alvo else (lo, mid)
    q['r_adm'] = lo
    np.set_printoptions(precision=4, suppress=True)
    log(f'  Z <- T1 {qa["a1"]:.4f}, Z <- T2 {qa["a2"]:.4f}, ||V Q_Z1|| {qa["VQ_Z1"]:.3f}, ||V Q_Z2|| {qa["VQ_Z2"]:.3f}')
    log(f'  N0 (Z, T1, T2, far) =\n{N0}')
    log(f'  raio de Perron: r = 0: {q["perron0"]:.4f};  tile: {q["perron"]:.4f};'
        f'  rZ admissivel (raio <= {alvo}): {q["r_adm"]:.4f}   {ctx.tempo()}')
    if RIG:
        ok, v, th = verifica_perron(Nfinal, 0.999)
        q['perron_v'] = v.tolist(); q['perron_theta'] = th; q['certificado'] = bool(ok)
        log(f'  VERIFICACAO: v = {np.round(v, 4)}, max (N v)_i / v_i = {th:.6f}  ->  '
            + ('K(s) INJETIVO no tile' if ok else 'NAO certificado'))
    for parte in ('T1', 'T2'):
        q[parte] = {k_: (str(v_) if isinstance(v_, (complex, tuple)) else v_) for k_, v_ in q[parte].items()}
    return q


def main():
    global RIG, BAND_M, BAND_N
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--s0', type=complex, default=0.0, help='centro das colunas de Z (o tile)')
    ap.add_argument('--r', type=float, default=0.0)
    ap.add_argument('--s1', type=complex, default=None, help='centro das colunas de T1')
    ap.add_argument('--r1', type=float, default=None)
    ap.add_argument('--s2', type=complex, default=None, help='centro das colunas de T2 e far')
    ap.add_argument('--r2', type=float, default=None)
    ap.add_argument('--cache1', type=Path, default=None, help='.npz da parte de T1 (le ou grava)')
    ap.add_argument('--cache2', type=Path, default=None, help='.npz da parte de T2 (le ou grava)')
    ap.add_argument('--Z', default='20x80')
    ap.add_argument('--Gp', default=None, help='casca exata G+ (padrao: G+ = Z)')
    ap.add_argument('--F', default='200x400', help='regiao finita Mc x Nc')
    ap.add_argument('--D', type=int, default=12)
    ap.add_argument('--rigoroso', action='store_true')
    ap.add_argument('--banda', default='20,40', help='fundo truncado para fora: |d| <= Dm, |dn| <= Dn')
    ap.add_argument('--saida', type=Path, help='JSON com as cotas')
    a = ap.parse_args()
    RIG = a.rigoroso
    BAND_M, BAND_N = map(int, a.banda.split(','))
    cx = lambda t: tuple(map(int, t.split('x')))
    Z = cx(a.Z); Gp = cx(a.Gp) if a.Gp else Z; F = cx(a.F)
    print(f'tile |s - {a.s0}| <= {a.r}; T1 em {a.s1} (r1 {a.r1}); T2 em {a.s2} (r2 {a.r2});'
          f' Z = {Z}, G+ = {Gp}, F = {F}, banda {BAND_M},{BAND_N}, rigoroso = {RIG}', flush=True)
    q = tile(a.s0, a.r, Z, Gp, F, D=a.D, s1=a.s1, r1=a.r1, s2=a.s2, r2=a.r2,
             cache1=a.cache1, cache2=a.cache2)
    if a.saida:
        import json
        q.update(s0=[a.s0.real, a.s0.imag], r=a.r, Zc=Z, Gp=Gp, F=F, D=a.D)
        a.saida.write_text(json.dumps(q, indent=1, default=str) + '\n')


if __name__ == '__main__':
    main()
