#!/usr/bin/env python3
"""Verifica a rota 5 no regime que nenhuma malha do projeto alcancou.

O limiar de contracao da recursao u = -(J_mu+s)^-1 B u foi estimado em
n ~ 51 (radial) e |m| ~ 18 (Fourier), a partir de malhas que vao so ate
n=35 e m=12. Tudo em ALTERNATIVE_ROUTES 5.3 e extrapolacao.

Este script responde tres perguntas numa malha que cruze os limiares:
  1. C(n), a maior soma de coluna ponderada de B por indice radial, satura
     em ~8,8 ou cresce?
  2. o fator de contracao C(n)/|D(n)| cruza 1 onde previsto e FICA abaixo?
  3. o decaimento observado das autofuncoes entra no regime previsto?

DIAGNOSTICO. Nao certifica nada.
"""
from __future__ import annotations
import sys, json
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from weighted_norms import enderecos, MU
from check_free_preconditioner import read_addresses, free_matrix

ETA = [229/1024, 1.0, 243/4096, 1233/4096]


def le_matriz_rapido(path: Path):
    """Leitor em ponto flutuante: Fraction e inviavel acima de ~10^6 entradas.

    Formato: linhas `i j numerador expoente`, valor = numerador * 2^expoente.
    O numerador e inteiro grande; int->float arredonda para o double mais
    proximo, o que basta para um DIAGNOSTICO (precisamos de 3 digitos).
    """
    dim = None
    with path.open(encoding='ascii') as fh:
        for linha in fh:
            if linha.startswith('dimension'):
                dim = int(linha.split()[1])
            elif linha.strip() == 'BEGIN_ENTRIES':
                break
        if dim is None:
            raise ValueError('dimensao ausente')
        A = np.zeros((dim, dim))
        for linha in fh:
            campos = linha.split()
            if len(campos) != 4:
                continue
            i, j, num, exp = campos
            A[int(i), int(j)] = float(int(num)) * 2.0 ** int(exp)
    return dim, A


def main():
    malha = sys.argv[1] if len(sys.argv) > 1 else '12x36'
    p = Path(f'build/spectrum/rt-A-{malha}.dat')
    import time as _t
    _t0 = _t.time()
    dim, A = le_matriz_rapido(p)
    print(f'matriz lida: {dim}x{dim} em {_t.time()-_t0:.0f}s', flush=True)
    addr_raw, addr = read_addresses(p, dim), enderecos(p, dim)
    J = np.array([[float(x) for x in row] for row in free_matrix(addr_raw, MU)])
    B = np.abs(A - J)
    w = np.array([(2-(m == 0))*(2-(n == 0))*(65/64)**m*(5/4)**n*ETA[c]
                  for (c, _, m, n) in addr])
    Bw = B*w[:, None]/w[None, :]
    n_i = np.array([a[3] for a in addr]); m_i = np.array([a[2] for a in addr])
    col = Bw.sum(0)
    mu = float(MU)

    print(f"=== malha {malha}, dim {dim}, n ate {n_i.max()}, m ate {m_i.max()} ===")
    print("\n1+2. C(n) e o fator de contracao radial")
    print(f"{'n':>4} {'C(n)':>9} {'mu(n+1)':>9} {'contracao':>10}")
    print('-'*36)
    radial = {}
    cruzou = None
    for n in range(int(n_i.max())+1):
        mask = n_i == n
        if not mask.any(): continue
        c = float(col[mask].max()); k = c/(mu*(n+1))
        radial[n] = dict(C=c, contracao=k)
        if n % max(1, int(n_i.max())//14) == 0 or (cruzou is None and k < 1):
            print(f"{n:>4} {c:>9.3f} {mu*(n+1):>9.3f} {k:>10.3f}"
                  + ("   <- cruza 1" if cruzou is None and k < 1 else ""))
        if cruzou is None and k < 1: cruzou = n
    print(f"\n  cruzamento observado: n = {cruzou}")
    print(f"  C maximo sobre toda a malha: {max(v['C'] for v in radial.values()):.3f}")
    altos = [v['C'] for n, v in radial.items() if n > 0.6*n_i.max()]
    print(f"  C medio no terco superior de n: {np.mean(altos):.3f}"
          f"  (satura?  desvio {np.std(altos):.3f})")

    print("\n  fator de contracao APOS o cruzamento (fica abaixo de 1?)")
    depois = {n: v['contracao'] for n, v in radial.items() if cruzou and n >= cruzou}
    if depois:
        print(f"    maximo apos o cruzamento: {max(depois.values()):.3f}")
        print(f"    valor no maior n: {depois[max(depois)]:.3f}")

    print("\n3. decaimento das autofuncoes")
    vals, vecs = np.linalg.eig(A)
    alvos = {'gauge': -0.1683, 'S3': -0.4010, 'fisica': -0.7332}
    saida_modos = {}
    for nome, alvo in alvos.items():
        k = int(np.argmin(np.abs(vals - alvo)))
        v = np.abs(vecs[:, k]); v = v/v.max()
        perfil = {}
        for n in range(int(n_i.max())+1):
            mask = n_i == n
            if mask.any(): perfil[n] = float(v[mask].max())
        ns = sorted(perfil)
        razoes = [perfil[b]/perfil[a] for a, b in zip(ns, ns[1:]) if perfil[a] > 0]
        meio = razoes[len(razoes)//3:2*len(razoes)//3]
        fim = razoes[int(0.75*len(razoes)):]
        print(f"  {nome:<8} autovalor {vals[k].real:>9.5f}"
              f"  razao media (meio) {np.mean(meio):.3f}"
              f"  (ultimo quarto) {np.mean(fim):.3f}")
        saida_modos[nome] = dict(autovalor=float(vals[k].real),
                                 razao_meio=float(np.mean(meio)),
                                 razao_fim=float(np.mean(fim)),
                                 perfil={str(n): perfil[n] for n in ns})

    rep = dict(malha=malha, dim=dim, n_max=int(n_i.max()), m_max=int(m_i.max()),
               cruzamento_radial=cruzou,
               C_max=max(v['C'] for v in radial.values()),
               C_terco_superior=float(np.mean(altos)),
               radial={str(n): v for n, v in radial.items()}, modos=saida_modos)
    saida = Path(f'build/spectrum/apriori-regime-{malha}.json')
    saida.write_text(json.dumps(rep, indent=2)+'\n')
    print(f"\nescrito: {saida}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
