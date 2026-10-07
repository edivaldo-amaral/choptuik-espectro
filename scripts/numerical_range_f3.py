#!/usr/bin/env python3
"""F3 pelo campo de valores: nao ha espectro a direita de R = max Re W(L).

Para T num espaco de Hilbert, spec(T) subset fecho(W(T)), com
W(T) = {<Tx,x>/<x,x>}. A borda direita de W da diretamente um R com
"nao ha espectro para Re s > R" -- e e uma estimativa de ENERGIA, que e
exatamente a evidencia que a linha F3 do ledger exige.

A realizacao e l^2 com os pesos de RT. A ponte para o l^1 do projeto e
gratuita NA DIRECAO QUE IMPORTA: para pesos positivos,
sum (w_j |x_j|)^2 <= (sum w_j |x_j|)^2, logo ||x||_{2,w} <= ||x||_{1,w} e
l^1(w) subset l^2(w). Toda autofuncao em l^1 e autofuncao em l^2, entao
excluir espectro em l^2 exclui em l^1.

DIAGNOSTICO. Nada aqui e certificado; os blocos sao finitos.
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_apriori_regime import le_matriz_rapido, ETA
from weighted_norms import enderecos, MU
from check_free_preconditioner import read_addresses, free_matrix

K2 = 1.25


def R_livre_1d(N: int = 400) -> float:
    """R_J pelo modelo 1-D da componente 3, onde o minimo de fato vive.

    Diagonal mu(n+1); acoplamento 2*mu*n para TODO n'<n (passo 1 fora da
    componente 0); pesos l^2 sqrt((2-[n=0]) kappa2^n). Estabiliza em N=20.
    """
    mu = float(MU)
    w = np.array([(2 if n else 1)*K2**n for n in range(N)])
    S = np.sqrt(w)
    J = np.zeros((N, N))
    for n in range(N):
        J[n, n] = mu*(n+1)
        for lo in range(n-1, -1, -1):
            J[lo, n] = 2*mu*n
    Jw = J*S[:, None]/S[None, :]
    return float(-np.linalg.eigvalsh((Jw+Jw.T)/2)[0])


def R_malha(malha: str) -> dict:
    p = Path(f'build/spectrum/rt-A-{malha}.dat')
    dim, A = le_matriz_rapido(p)
    addr_raw, addr = read_addresses(p, dim), enderecos(p, dim)
    J = np.array([[float(x) for x in row] for row in free_matrix(addr_raw, MU)])
    w = np.array([(2-(m == 0))*(2-(n == 0))*(65/64)**m*K2**n*ETA[c]
                  for (c, _, m, n) in addr])
    S = np.sqrt(w)
    T = S[:, None]/S[None, :]
    herm = lambda X: (X+X.conj().T)/2
    return dict(malha=malha, dim=dim,
                R=float(np.linalg.eigvalsh(herm(-A*T))[-1]),
                R_J=float(np.linalg.eigvalsh(herm(-J*T))[-1]),
                B2=float(np.linalg.norm((A-J)*T, 2)))


def main():
    malhas = sys.argv[1:] or ['4x12', '6x18', '8x24', '12x36']
    print(f"R_J pelo modelo 1-D da componente 3: {R_livre_1d():.9f}")
    print(f"\n{'malha':>7} {'dim':>6} {'R':>9} {'R_J':>11} {'||B||_2':>9}")
    print('-'*46)
    saida = []
    for m in malhas:
        r = R_malha(m)
        saida.append(r)
        print(f"{m:>7} {r['dim']:>6} {r['R']:>9.4f} {r['R_J']:>11.7f} {r['B2']:>9.4f}")
    print("\nR certificado hoje por Neumann grosseiro: 414")
    Path('build/spectrum/numerical-range-f3.json').write_text(
        json.dumps(dict(R_J_modelo_1d=R_livre_1d(), malhas=saida,
                        R_atual_por_neumann=414,
                        nota='diagnostico; l^2 com pesos RT; l^1(w) subset l^2(w)'),
                   indent=2)+'\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
