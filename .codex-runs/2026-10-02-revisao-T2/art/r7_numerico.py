#!/usr/bin/env python3
"""R7 (revisao T2): consequencia testavel de K C = M L nos blocos do codigo (truncagem).

M(s) nao esta implementado na base com sinal, entao K(s)C(s) - M(s)L(s) nao pode ser montado.
O que se testa e a consequencia: para h em ker L(s0) (truncado), K(s0) C(s0) h = M(s0) L(s0) h ~ 0
(a menos de truncagem e do defeito G(h, e(RefA)) do fundo aproximado). Para raizes que VIOLAM as
constraints (C h != 0), isso obriga C(s0) h a estar no nucleo (aproximado) de K(s0): testamos o
cosseno entre C(s0)h e o vetor singular inferior de K(s0), e ||K C h||/(||K|| ||C h||), comparando com
um vetor h aleatorio (controle negativo).

Operadores: signed_operator_L.operador (L), signed_constraints.operador_C (C), signed_operator.operador (K);
mesmos layouts de saida de C e de entrada de K (m par, comp 0 com n impar, comp 1 com todo n).
"""
import sys
import json
import numpy as np

sys.path.insert(0, 'scripts')
import signed_operator_L as sl
import signed_constraints as sc
import signed_operator as so


def main():
    out = []
    for M, N in ((12, 36), (16, 48)):
        L0, rs, offs = sl.operador(M, N, 0.0)
        lam, V = np.linalg.eig(L0)
        raizes = -lam                                    # L(s) = L0 + s singular em s = -lam
        rng = np.random.default_rng(5)
        for alvo in (0.40102, 0.73318, 0.0414 + 0.5j, 0.16831):
            k = int(np.argmin(np.abs(raizes - alvo)))
            s0 = raizes[k]; h = V[:, k]
            res_L = np.linalg.norm((L0 + s0*np.eye(len(L0))) @ h)/np.linalg.norm(h)
            C, sos, ris, oo, oi = sc.operador_C(M, N, s0)
            K, _ = so.operador(M, N, s0)
            assert K.shape[0] == C.shape[0]
            c = C @ h
            Kc = K @ c
            nK = np.linalg.norm(K, 2)
            U, S, Vh = np.linalg.svd(K)
            vmin = Vh[-1].conj()
            cos = abs(np.vdot(vmin, c))/np.linalg.norm(c)
            hr = rng.standard_normal(len(h)) + 1j*rng.standard_normal(len(h))
            cr = C @ hr
            cos_r = abs(np.vdot(vmin, cr))/np.linalg.norm(cr)
            linha = dict(M=M, N=N, alvo=str(alvo), s0=str(s0), residuo_L=float(res_L),
                         Ch_sobre_h=float(np.linalg.norm(c)/np.linalg.norm(h)),
                         KCh_rel=float(np.linalg.norm(Kc)/(nK*np.linalg.norm(c))),
                         sigma_min_K=float(S[-1]), sigma_min_K_rel=float(S[-1]/S[0]),
                         cos_Ch_nucleoK=float(cos), controle_aleatorio_KCh_rel=float(np.linalg.norm(K @ cr)/(nK*np.linalg.norm(cr))),
                         controle_aleatorio_cos=float(cos_r))
            out.append(linha)
            print(f'{M}x{N} s0={s0:.6f}: ||L h||/||h|| {res_L:.1e}; ||C h||/||h|| {linha["Ch_sobre_h"]:.3e}; '
                  f'||K C h||/(||K|| ||C h||) {linha["KCh_rel"]:.2e} (aleatorio {linha["controle_aleatorio_KCh_rel"]:.2e}); '
                  f'sigma_min(K)/||K|| {linha["sigma_min_K_rel"]:.2e}; cos(Ch, nucleo K) {cos:.6f} (aleatorio {cos_r:.2e})',
                  flush=True)
    json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else 'r7_numerico.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
