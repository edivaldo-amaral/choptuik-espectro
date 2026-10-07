#!/usr/bin/env python3
"""F1 recalculado com codigo proprio (revisao independente do Lema G, 03/10/2026).

Nao reutiliza scripts/gauge_global_centro.py.  Le RefA.dat e RefAplusB.dat diretamente,
no formato documentado em RT (secao "Data structures"): Field com TwoExp, off_m, off_n,
num_m, num_n e data[u][v] = Re + i Im (inteiros), v_{off_m+u, off_n+v} = 2^TwoExp data[u][v].

Calcula, em racionais exatos (Fraction / inteiros):
  c_m(xi0) = sum_{n in Z} v_{mn} e^{i n theta0},  theta0 = arccos xi0,
para o campo de indice 3 (omega4), m = 1, em xi0 = 0 (theta = pi/2, e^{in theta} = i^n)
e, como bonus, em xi0 = -1 (theta = pi; usado para omega4^+ no cone xi = 1).
Usa v_{m,-n} = v_{mn} (definicao do espaco V de RT).

Tambem:
  - confere paridades (omega4 so com m impar; omega1 so com m par e n impar; omega2,3 so m par);
  - confere Re v_{1,0}(omega4) = 0 (espaco Gauged de RT);
  - calcula c_1 de RefAplusB (exato) e |c_1(RefAplusB) - c_1(RefA)|;
  - calcula a norma SSPACE exata de RefB = RefAplusB - RefA (todas as 4 componentes), para
    conferir L2 = 2^-25 da bola de RT;
  - confere por avaliacao numerica (mpmath) de omega4(tau, 0) e FFT que c_1 e o coeficiente
    de Fourier e^{i tau/2} de omega4(., 0) (convencao da base e^{i m tau/2 + i n theta}).
"""
import json
import re
import sys
import time
from fractions import Fraction

import mpmath as mp

RAIZ = "<repo>/.cache/rt-1203.3766v1/"
ARQ_A = RAIZ + "sourcecode/RefA.dat"
ARQ_AB = RAIZ + "RefAplusB.dat"
PAR = re.compile(r"\((-?\d+),(-?\d+)\)")


def ler_multifield(caminho, componentes=None, linhas_m=None):
    """Retorna (mu como Fraction, lista de campos).

    Cada campo: dict com TwoExp, off_m, off_n, num_m, num_n, e 'dados' = lista de linhas
    (u = indice m - off_m) de listas de (Re, Im) inteiros.  Se componentes for dado, so
    estas componentes sao decodificadas.  Se linhas_m for dado, so estas linhas m."""
    with open(caminho) as f:
        linhas = f.read().split("\n")
    i = 0
    assert linhas[i].strip() == "mu_MultiField"; i += 1
    assert linhas[i].strip() == "DyadicQ"; i += 1
    te = int(linhas[i].split()[1]); i += 1
    co = int(linhas[i].split()[1]); i += 1
    mu = Fraction(co) * Fraction(2) ** te
    assert linhas[i].strip() == "MultiField"; i += 1
    nc = int(linhas[i].split()[1]); i += 1
    campos = []
    for c in range(nc):
        assert linhas[i].strip() == "Field"; i += 1
        cab = {}
        for chave in ("TwoExp", "off_m", "off_n", "num_m", "num_n"):
            k, v = linhas[i].split()
            assert k == chave, (k, chave)
            cab[chave] = int(v)
            i += 1
        texto = linhas[i]; i += 1
        if componentes is not None and c not in componentes:
            campos.append(cab)
            continue
        assert texto.startswith("(((") and texto.endswith(")))")
        # o texto e "(" + linha0 + "," + linha1 + ... + ")", linha = "(" + pares + ")",
        # par = "(Re,Im)", na ordem u (= m - off_m) e depois v (= n - off_n) (Field_to_file de RT).
        # Confere a estrutura contando as linhas: num_m - 1 separadores ")),((".
        assert texto.count(")),((") == cab["num_m"] - 1
        todos = PAR.findall(texto)
        assert len(todos) == cab["num_m"] * cab["num_n"], (len(todos), cab["num_m"], cab["num_n"])
        N = cab["num_n"]
        dados = []
        for u in range(cab["num_m"]):
            m = cab["off_m"] + u
            if linhas_m is not None and m not in linhas_m:
                dados.append(None)
                continue
            dados.append([(int(a), int(b)) for a, b in todos[u * N:(u + 1) * N]])
        del todos
        cab["dados"] = dados
        campos.append(cab)
    return mu, campos


def coef(campo, m, n):
    """v_{mn} para m, n >= 0 como (Re, Im) Fraction."""
    u = m - campo["off_m"]
    v = n - campo["off_n"]
    if u < 0 or v < 0 or u >= campo["num_m"] or v >= campo["num_n"]:
        return (Fraction(0), Fraction(0))
    a, b = campo["dados"][u][v]
    esc = Fraction(2) ** campo["TwoExp"]
    return (a * esc, b * esc)


def c_m_xi0(campo, m, fase):
    """sum_{n in Z} v_{mn} e^{i n theta0}, com v_{m,-n} = v_{mn}.

    fase(n) = cos(n theta0) racional (theta0 = pi/2: cos(n pi/2); theta0 = pi: (-1)^n).
    sum_{n in Z} v_{mn} e^{i n th} = v_{m0} + sum_{n>=1} v_{mn} (e^{i n th} + e^{-i n th})
                                  = v_{m0} + 2 sum_{n>=1} v_{mn} cos(n th)."""
    re_, im_ = coef(campo, m, 0)
    for n in range(1, campo["off_n"] + campo["num_n"]):
        a, b = coef(campo, m, n)
        f = 2 * fase(n)
        re_ += a * f
        im_ += b * f
    return re_, im_


def cos_npi2(n):
    return [1, 0, -1, 0][n % 4]


def cos_npi(n):
    return 1 if n % 2 == 0 else -1


def paridades(campos):
    res = {}
    for c, campo in enumerate(campos):
        m_par_nz = m_impar_nz = n_par_nz = n_impar_nz = 0
        for u, linha in enumerate(campo["dados"]):
            m = campo["off_m"] + u
            for v, (a, b) in enumerate(linha):
                n = campo["off_n"] + v
                if a or b:
                    if m % 2 == 0:
                        m_par_nz += 1
                    else:
                        m_impar_nz += 1
                    if n % 2 == 0:
                        n_par_nz += 1
                    else:
                        n_impar_nz += 1
        # v_{0n} real (definicao de V para m = 0)
        im_m0 = sum(1 for (a, b) in campo["dados"][0] if b != 0) if campo["off_m"] == 0 else None
        res[c] = dict(nz_m_par=m_par_nz, nz_m_impar=m_impar_nz, nz_n_par=n_par_nz,
                      nz_n_impar=n_impar_nz, im_nao_nulo_em_m0=im_m0)
    return res


def omega4_eixo_numerico(campo, taus, dps=40):
    """omega4(tau, 0) = sum_{m,n in Z} v_{mn} e^{i m tau/2} i^n, com v_{-m,-n} = conj(v_{mn}),
    v_{m,-n} = v_{mn}.  Soma completa sobre m em Z (m >= 0 armazenado; m < 0 por conjugacao)."""
    mp.mp.dps = dps
    # c_m para m >= 0
    cms = {}
    for u in range(campo["num_m"]):
        m = campo["off_m"] + u
        re_, im_ = c_m_xi0(campo, m, cos_npi2)
        cms[m] = mp.mpc(mp.mpf(re_.numerator) / re_.denominator, mp.mpf(im_.numerator) / im_.denominator)
    vals = []
    for t in taus:
        s = mp.mpc(0)
        for m, c in cms.items():
            if m == 0:
                s += c
            else:
                # termo m e termo -m: c_{-m} = sum_n v_{-m,n} i^n = sum_n conj(v_{m,-n}) i^n
                #   = sum_n conj(v_{m,n}) i^n  (v_{m,-n}=v_{mn}) = conj(sum_n v_{mn} (-i)^n)
                #   = conj(sum_n v_{mn} i^{-n}) = conj(sum_n v_{m,-n} i^{n}) = conj(c_m)
                s += c * mp.expj(m * t / 2) + mp.conj(c) * mp.expj(-m * t / 2)
        vals.append(s)
    return cms, vals


def norma_V_exata(dif_campo_int, two_exp, num_m, num_n):
    """||v||_V = sum_{m,n>=0} (2-d_m0)(2-d_n0) (65/64)^m (5/4)^n (|Re|+|Im|) * 2^two_exp,
    em inteiros exatos.  dif_campo_int[u][v] = (Re, Im) inteiros (m = u, n = v)."""
    # S_m = sum_n c_n 5^n 4^(N-1-n) (|a|+|b|)  -> dividido por 4^(N-1)
    N = num_n
    pot5 = [5 ** n for n in range(N)]
    pot4 = [4 ** (N - 1 - n) for n in range(N)]
    M = num_m
    total = 0
    for m in range(M):
        linha = dif_campo_int[m]
        Sm = 0
        for n in range(N):
            a, b = linha[n]
            w = abs(a) + abs(b)
            if w:
                Sm += (1 if n == 0 else 2) * pot5[n] * pot4[n] * w
        total += (1 if m == 0 else 2) * 65 ** m * 64 ** (M - 1 - m) * Sm
    den = 4 ** (N - 1) * 64 ** (M - 1)
    return Fraction(total, den) * Fraction(2) ** two_exp


def main():
    t0 = time.time()
    out = {}
    muA, A = ler_multifield(ARQ_A)
    out["mu_RefA"] = str(muA)
    out["mu_RefA_float"] = float(muA)
    out["mu_RefA_esperado_RT"] = "722873400*2^-32"
    assert muA == Fraction(722873400, 2 ** 32)
    out["cabecalhos_RefA"] = [{k: c[k] for k in ("TwoExp", "off_m", "off_n", "num_m", "num_n")} for c in A]
    par = paridades(A)
    out["paridades_RefA"] = par
    # omega4 = componente 3: so m impar (TauAntiPer); omega1: m par, n impar; omega2,3: m par
    assert par[3]["nz_m_par"] == 0 and par[3]["nz_m_impar"] > 0
    assert par[0]["nz_m_impar"] == 0 and par[0]["nz_n_par"] == 0
    assert par[1]["nz_m_impar"] == 0 and par[2]["nz_m_impar"] == 0
    w4 = A[3]
    re10, im10 = coef(w4, 1, 0)
    out["RefA4_10"] = [str(re10), str(im10)]
    assert re10 == 0, "Re v_{1,0}(omega4) deveria ser 0 (Gauged)"
    c1r, c1i = c_m_xi0(w4, 1, cos_npi2)
    out["c1_RefA_xi0"] = {"re": str(c1r), "im": str(c1i), "re_float": float(c1r), "im_float": float(c1i),
                          "abs_float": float(abs(complex(float(c1r), float(c1i))))}
    # bonus: xi = -1 (theta = pi): omega4(tau,-1); omega4^+(tau, 1) = -omega4(tau, -1)
    c1r_m1, c1i_m1 = c_m_xi0(w4, 1, cos_npi)
    out["c1_RefA_xi_menos1"] = {"re_float": float(c1r_m1), "im_float": float(c1i_m1)}
    # c_m em xi = 0 para m = 1, 3, 5 (diagnostico)
    out["c_m_RefA_xi0_float"] = {}
    for m in (1, 3, 5, 7):
        a, b = c_m_xi0(w4, m, cos_npi2)
        out["c_m_RefA_xi0_float"][m] = [float(a), float(b)]

    # Conferencia numerica da convencao: omega4(tau,0) por soma direta e coeficiente de Fourier m=1
    K = 64
    taus = [4 * mp.pi * k / K for k in range(K)]
    cms, vals = omega4_eixo_numerico(w4, taus)
    # coeficiente de e^{i tau/2}: (1/K) sum_k omega(tau_k) e^{-i tau_k/2}
    coefF = sum(v * mp.expj(-t / 2) for v, t in zip(vals, taus)) / K
    out["conferencia_fft"] = {
        "c1_por_fft": [float(coefF.real), float(coefF.imag)],
        "max_abs_Im_omega4": float(max(abs(v.imag) for v in vals)),
        "antiperiodicidade_max": float(max(abs(vals[k] + vals[(k + K // 2) % K]) for k in range(K))),
        "omega4_eixo_amostras": [float(vals[k].real) for k in range(0, K, 8)],
    }
    # zeros de omega4(., 0) em [0, 4pi): mudancas de sinal na amostragem fina
    K2 = 2048
    taus2 = [4 * mp.pi * k / K2 for k in range(K2)]
    _, vals2 = omega4_eixo_numerico(w4, taus2, dps=20)
    sinais = [1 if v.real > 0 else -1 for v in vals2]
    trocas = sum(1 for k in range(K2) if sinais[k] != sinais[(k + 1) % K2])
    out["omega4_eixo_trocas_de_sinal_em_[0,4pi)"] = trocas
    out["omega4_eixo_min_abs_amostra"] = float(min(abs(v.real) for v in vals2))
    print("RefA ok", time.time() - t0, file=sys.stderr)

    # RefAplusB: so a componente 3, linha m = 1, para c_1 exato
    muAB, AB = ler_multifield(ARQ_AB, componentes={3}, linhas_m={1})
    out["mu_RefAplusB_float"] = float(muAB)
    out["mu_RefB_float"] = float(muAB - muA)
    w4AB = AB[3]
    c1r_ab, c1i_ab = c_m_xi0(w4AB, 1, cos_npi2)
    d_re = c1r_ab - c1r
    d_im = c1i_ab - c1i
    out["c1_RefAplusB_xi0"] = {"re_float": float(c1r_ab), "im_float": float(c1i_ab)}
    out["c1_RefAplusB_menos_RefA"] = {"re_float": float(d_re), "im_float": float(d_im),
                                      "cota_l1": float(abs(d_re) + abs(d_im)),
                                      "cota_l1_menor_que_2^-25": bool(abs(d_re) + abs(d_im) < Fraction(1, 2 ** 25))}
    # cota inferior rigorosa de |c_1(omega*)| usando so a bola R = 2^-277 em torno de RefAplusB
    # |c1(w*)| >= |Im c1(RefAplusB)| - 2^-277
    lb_AB = abs(c1i_ab) - Fraction(1, 2 ** 277)
    lb_A = abs(c1i) - Fraction(1, 2 ** 25) - Fraction(1, 2 ** 277)
    out["cota_inferior_|c1(omega*)|_via_RefAplusB"] = float(lb_AB)
    out["cota_inferior_|c1(omega*)|_via_RefA_e_L2"] = float(lb_A)
    assert lb_AB > 0 and lb_A > 0
    print("c1 RefAplusB ok", time.time() - t0, file=sys.stderr)

    with open("<repo>/.codex-runs/2026-10-03-revisao-G/f1_proprio.json", "w") as f:
        json.dump(out, f, indent=1)

    # Norma SSPACE exata de RefB (todas as componentes), para conferir L2 = 2^-25
    if "--refB" in sys.argv:
        normas = []
        for c in range(4):
            _, ABc = ler_multifield(ARQ_AB, componentes={c})
            fab = ABc[c]
            fa = A[c]
            assert fab["off_m"] == 0 and fab["off_n"] == 0
            te = fab["TwoExp"]
            sh = fa["TwoExp"] - te  # RefA em unidades 2^te
            assert sh >= 0
            dif = []
            for u in range(fab["num_m"]):
                linha = []
                for v in range(fab["num_n"]):
                    a, b = fab["dados"][u][v]
                    if u < fa["num_m"] and v < fa["num_n"]:
                        a2, b2 = fa["dados"][u][v]
                        a -= a2 << sh
                        b -= b2 << sh
                    linha.append((a, b))
                dif.append(linha)
            nv = norma_V_exata(dif, te, fab["num_m"], fab["num_n"])
            normas.append(nv)
            print("componente", c, float(nv), time.time() - t0, file=sys.stderr)
            del ABc, fab, dif
        tot = sum(normas)
        out["norma_SSPACE_RefB"] = {"por_componente": [float(x) for x in normas], "total": float(tot),
                                    "total_vezes_2^25": float(tot * 2 ** 25),
                                    "<=2^-25": bool(tot <= Fraction(1, 2 ** 25))}
        with open("<repo>/.codex-runs/2026-10-03-revisao-G/f1_proprio.json", "w") as f:
            json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
    print("tempo", time.time() - t0, file=sys.stderr)


if __name__ == "__main__":
    main()
