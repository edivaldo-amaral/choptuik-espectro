#!/usr/bin/env python3
"""Recomputo INDEPENDENTE (auditoria de 16/09) dos majorantes do fundo.

Nao importa modulos do repositorio nem le GAMMA2_macro.c. As formulas vem
diretamente do TeX de RT (arXiv:1203.3766v1):
  * eq. (23): S Xi Gamma2(v,w) em convolucoes com campos coeficiente;
  * definicao de Gamma1 e S (secao 3), dando S Gamma1 h = (0, R(h1+(1-P)h2), -R h1, R(1-P)h4);
  * Gamma2^sharp e Gamma1^sharp (secao 3, "The sharp system");
  * norma (secao 4): sum (2-d_m0)(2-d_n0) kappa1^m kappa2^n (|Re|+|Im|);
  * (34) Lipschitz 3; (39b) L2=2^-25; (60) R=2^-277; secao 5: L1=10.5.

Saidas: majorante 4x4 do bloco A (pesos fortes), beta na norma eta,
checagem cruzada com L1=10.5 de RT, majorante sharp (pesos menores),
versoes filtradas por suporte (acoplamento F->T), e comparacao EXATA com os
relatorios auditados.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
REFA = ROOT/'.cache/rt-1203.3766v1/sourcecode/RefA.dat'
REFA_SHA = '5776e6038f994e06b4b82fc638d797154f630aeee37b251058512d3cb4cd450b'
EPS_OMEGA = Q(1, 2**25)+Q(1, 2**277)
ETA = (Q(916), Q(4096), Q(243), Q(1233))


def parse_refa(path: Path = REFA):
    """Quatro campos: dict (m,n) -> (re,im) inteiros, escala 2^TwoExp comum."""
    text = path.read_text()
    parts = text.split('\nField\n')[1:]
    if len(parts) != 4:
        raise ValueError('esperados 4 campos')
    fields, exps = [], set()
    for part in parts:
        head = dict(re.findall(r'^(TwoExp|off_m|off_n|num_m|num_n)\s+(-?\d+)$', part, re.M))
        head = {k: int(v) for k, v in head.items()}
        body = part[part.index('('):]
        pairs = re.findall(r'\((-?\d+),(-?\d+)\)', body)
        if len(pairs) != head['num_m']*head['num_n'] or head['off_m'] or head['off_n']:
            raise ValueError('formato inesperado')
        exps.add(head['TwoExp'])
        field = {}
        for idx, (a, b) in enumerate(pairs):
            m, n = divmod(idx, head['num_n'])
            if int(a) or int(b):
                field[m, n] = (int(a), int(b))
        fields.append(field)
    if len(exps) != 1:
        raise ValueError('expoentes distintos')
    # Validacao estrutural da leitura (m externo, n interno):
    # campo1 in TP cap XiOdd, campos 2,3 in TP, campo 4 in TAP.
    assert all(m % 2 == 0 and n % 2 == 1 for m, n in fields[0])
    assert all(m % 2 == 0 for m, n in fields[1]) and all(m % 2 == 0 for m, n in fields[2])
    assert all(m % 2 == 1 for m, n in fields[3])
    assert all(b == 0 for f in fields for (m, n), (a, b) in f.items() if m == 0)
    return fields, Q(2)**exps.pop()


def lin(*terms):
    """Combinacao linear de campos: terms = (coef, field, apply_P)."""
    out = {}
    for coef, field, parity in terms:
        for (m, n), (a, b) in field.items():
            s = coef*((-1)**n if parity else 1)
            x = out.get((m, n), (Q(0), Q(0)))
            out[m, n] = (x[0]+s*a, x[1]+s*b)
    return out


def norm(field, scale, k1, k2, keep=lambda m, n: True):
    return scale*sum(((2-(m == 0))*(2-(n == 0))*k1**m*k2**n*(abs(a)+abs(b))
                      for (m, n), (a, b) in field.items() if keep(m, n)), Q(0))


def selected_coefficients(v):
    """Tabela de (23): para cada entrada (j, paridade) o vetor de 4 campos.

    Retorna coeficientes de S Xi Gamma2(v, .) aplicados a
      w1 (fator 1/4), w_j^even=(w_j+Pw_j)/2 e w_j^odd (fator 1/2).
    """
    v1, v2, v3, v4 = v
    Z = {}
    P = True
    tab = {
        ('w1', None): (Q(1, 4), [lin((1, v2, P), (1, v2, 0), (-1, v3, P), (-1, v3, 0)),
                                 lin((1, v1, 0), (-1, v2, 0), (1, v2, P)),
                                 lin((-1, v1, 0), (1, v2, 0), (-1, v2, P)), Z]),
        ('w2', 'e'): (Q(1, 2), [lin((1, v1, 0), (1, v3, P), (-1, v3, 0)),
                                lin((2, v2, 0), (2, v2, P)),
                                lin((-1, v2, 0), (-1, v2, P)),
                                lin((1, v4, P), (1, v4, 0))]),
        ('w2', 'o'): (Q(1, 2), [Z, lin((-1, v1, 0), (1, v2, P), (-1, v2, 0)),
                                lin((1, v1, 0)), lin((1, v4, P), (-1, v4, 0))]),
        ('w3', 'e'): (Q(1, 2), [lin((-1, v1, 0)), Z, Z, Z]),
        ('w3', 'o'): (Q(1, 2), [lin((-1, v2, 0), (-1, v2, P)), Z, Z, Z]),
        ('w4', 'e'): (Q(1, 2), [lin((1, v4, 0), (-1, v4, P)), Z,
                                lin((1, v4, P), (1, v4, 0)), lin((1, v2, P), (1, v2, 0))]),
        ('w4', 'o'): (Q(1, 2), [lin((1, v4, 0), (1, v4, P)), Z,
                                lin((1, v4, P), (-1, v4, 0)), lin((1, v2, P), (-1, v2, 0))]),
    }
    return tab


def sharp_coefficients(v):
    """Gamma2^sharp(v,.) (secao 3): entradas w1, w2#^even, w2#^odd, fator 1/2."""
    v1, v2, v3, v4 = v
    P = True
    return {
        ('w1', None): (Q(1, 2), [lin((1, v2, 0), (1, v2, P), (-1, v3, P), (-1, v3, 0)),
                                 lin((1, v1, 0), (-1, v2, 0), (1, v2, P))]),
        ('w2', 'e'): (Q(1, 2), [lin((-1, v1, 0), (1, v3, P), (-1, v3, 0)),
                                lin((-1, v1, 0), (1, v2, 0), (3, v2, P))]),
        ('w2', 'o'): (Q(1, 2), [lin((1, v2, 0), (1, v2, P), (1, v3, P), (1, v3, 0)),
                                lin((-1, v1, 0), (1, v2, 0), (3, v2, P))]),
    }


def filtered(fields, dm, dn):
    return [{k: x for k, x in f.items() if k[0] >= dm or k[1] >= dn} for f in fields]


def selected_majorant(fields, scale, k1=Q(65, 64), k2=Q(5, 4), factor=2):
    """A[i][j] com ||(factor*S Xi Gamma2(v,h))_i|| <= sum_j A_ij ||h_j||."""
    tab = selected_coefficients(fields)
    A = [[Q(0)]*4 for _ in range(4)]
    seven = []
    for (inp, par), (c, vec) in tab.items():
        j = int(inp[1])-1
        norms = [factor*c*norm(f, scale, k1, k2) for f in vec]
        seven.append(sum(norms)/factor)  # norma da multicampo de (33), fator 1
        for i in range(4):
            A[i][j] = max(A[i][j], norms[i])
    return A, max(seven)


def gamma1_majorant(mu_max, k2):
    R = 2/k2/(1-1/k2**2)       # ||Xi^-1_reg|| (RT, secao 4, via (35))
    d = mu_max*R
    return [[Q(0)]*4, [d, 2*d, Q(0), Q(0)], [d, Q(0), Q(0), Q(0)], [Q(0), Q(0), Q(0), 2*d]]


def eta_norm(A, eta=ETA):
    return max(sum(eta[i]*A[i][j] for i in range(len(A)))/eta[j] for j in range(len(A)))


def sharp_majorant(fields, scale, mu_max, k1=Q(129, 128), k2=Q(9, 8)):
    tab = sharp_coefficients(fields)
    A = [[Q(0)]*2 for _ in range(2)]
    for (inp, par), (c, vec) in tab.items():
        j = int(inp[1])-1
        for i in range(2):
            A[i][j] = max(A[i][j], c*norm(vec[i], scale, k1, k2))
    A[1][0] += mu_max*(2/k2/(1-1/k2**2))
    return A


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    raw = REFA.read_bytes()
    if hashlib.sha256(raw).hexdigest() != REFA_SHA:
        raise SystemExit('RefA diferente do auditado')
    fields, scale = parse_refa()
    A, seven = selected_majorant(fields, scale)
    lin_part = gamma1_majorant(Q(17, 100), Q(5, 4))
    total = [[a+b for a, b in zip(x, y)] for x, y in zip(A, lin_part)]
    ratio = max(ETA)/min(ETA)
    beta = eta_norm(total)+6*EPS_OMEGA*ratio
    unit = eta_norm(total, (1, 1, 1, 1))
    out = dict(refa_sha256=REFA_SHA,
               rt_L1_crosscheck_max_of_33=str(seven), rt_L1_crosscheck_decimal=float(seven),
               rt_L1_crosscheck_passes=seven <= Q(21, 2),
               nonlinear_majorant=[[str(x) for x in row] for row in A],
               beta_eta=str(beta), beta_eta_decimal=float(beta),
               beta_eta_le_5191_500=beta <= Q(5191, 500),
               unit_weight_bound=float(unit))
    # Comparacao exata com o produtor auditado.
    rep = json.loads((ROOT/'build/spectrum/component-exterior.json').read_text())
    out['nonlinear_equals_component_exterior_json'] = all(
        Q(rep['nonlinear_majorant'][i][j]) == A[i][j] for i in range(4) for j in range(4))
    out['beta_equals_component_exterior_json'] = Q(rep['weighted_beta']) == beta
    # Acoplamento F=6x18 -> T=fora de G, so parte bilinear de RefA filtrada.
    guarded = json.loads((ROOT/'build/spectrum/guarded-exterior.json').read_text())
    qF = Q(guarded['free_bound'])
    rows = []
    for (M, N), row in zip(((6, 18), (12, 36), (24, 72), (64, 192), (176, 944)), guarded['rows']):
        Ag, _ = selected_majorant(filtered(fields, M-6+1, N-18+1), scale)
        bound = (eta_norm(Ag)+6*EPS_OMEGA*ratio)*qF
        rows.append(dict(outer=[M, N], T_from_F=float(bound),
                         equals_report=Q(row['T_from_G']) == bound))
    out['guarded_rows'] = rows
    # Sharp: pesos menores, sem reponderacao, fator 1, + 3 eps.
    S = sharp_majorant(fields, scale, Q(17, 100))
    beta_s = max(S[0][j]+S[1][j] for j in range(2))+3*EPS_OMEGA
    out['sharp_majorant'] = [[float(x) for x in row] for row in S]
    out['sharp_beta'] = str(beta_s)
    out['sharp_beta_decimal'] = float(beta_s)
    out['sharp_beta_le_4_62237'] = beta_s <= Q(462237, 100000)
    Sg = sharp_majorant(filtered(fields, 160-12+1, 960-36+1), scale, Q(0))
    out['sharp_guard_160x960_reference_part'] = str(max(Sg[0][j]+Sg[1][j] for j in range(2)))
    text = json.dumps(out, indent=2)
    if args.out:
        args.out.write_text(text+'\n')
    print(text)


if __name__ == '__main__':
    main()
