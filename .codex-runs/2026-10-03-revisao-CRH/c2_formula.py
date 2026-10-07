"""C2 (revisao CRH): teste em dimensao finita da contabilidade da deflacao de posto 2 usada em
C1_REAVALIACAO.md par. 9 e 13, e de quando ela FALHA.

L(s) = L0 + s I (como em L-afim), x0, x1 autovetores EXATOS (L(s)x0 = s x0, L(s)x1 = (s - mu)x1),
phi_i com phi_i^* x_j = delta_ij + e_ij, E = d0 x0 phi0^* + d1 x1 phi1^*, Lt = L + E.
  (a) det Lt(s) / det L(s) = q(s)/(s(s - mu)), q(s) = (s + d0(1 + e00))(s - mu + d1(1 + e11)) - d0 d1 e01 e10;
  (b) se 0 e' raiz SIMPLES de L, Lt(0) e' invertivel;
  (c) se 0 tem bloco de Jordan 2 (x0 autovetor) ou e' semissimples dupla, Lt(0) e' SINGULAR:
      ou seja, 'Lt(0) invertivel' certifica mult_L(0) = 1;
  (d) se o vetor deflacionado NAO e' exato (x0 + eta), o determinante deixa de ser racional com polo
      simples em 0: a razao det Lt/det L muda O(|eta|/|s|) perto de 0 - por isso o certificado precisa
      dos vetores exatos (o codigo os usa via dEZ/dET).
"""
import numpy as np

rng = np.random.default_rng(20261003)
n, mu = 9, 0.16831
d0, d1 = 1.0, 1.0 + mu


def montar(jordan=None):
    """L0 = S diag(...) S^-1 com raizes de L(s) em 0, mu e outras; opcionalmente Jordan em 0."""
    lam = np.array([0.0, mu, 0.7, 1.3 + 0.2j, 1.3 - 0.2j, 2.1, -0.5, 3.0, 0.4j])  # raizes de L(s) em s = lam (L0 tem autovalores -lam)
    D = np.diag(lam).astype(complex)
    D = -D                                  # L0 x = -lam x  => L(s)x = (s - lam)x
    if jordan == 'jordan':
        D[6, 6] = 0.0                       # segunda raiz em 0 ...
        D[0, 6] = 1.0                       # ... num bloco de Jordan com o autovetor e0
    elif jordan == 'dupla':
        D[6, 6] = 0.0
    S = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
    L0 = S @ D @ np.linalg.inv(S)
    x0, x1 = S[:, 0], S[:, 1]               # autovetores exatos: raiz 0 e raiz mu
    return L0, x0, x1, S


def deflacao(L0, x0, x1, S, e=1e-7, eta=0.0):
    Sinv = np.linalg.inv(S)
    phi0 = Sinv[0] + e*(rng.normal(size=n) + 1j*rng.normal(size=n))
    phi1 = Sinv[1] + e*(rng.normal(size=n) + 1j*rng.normal(size=n))
    xx0 = x0 + eta*(rng.normal(size=n) + 1j*rng.normal(size=n))
    E = d0*np.outer(xx0, phi0) + d1*np.outer(x1, phi1)
    eij = np.array([[phi0 @ x0 - 1, phi0 @ x1], [phi1 @ x0, phi1 @ x1 - 1]])
    return E, eij


I = np.eye(n)
L0, x0, x1, S = montar()
E, eij = deflacao(L0, x0, x1, S)
q = lambda s: (s + d0*(1 + eij[0, 0]))*(s - mu + d1*(1 + eij[1, 1])) - d0*d1*eij[0, 1]*eij[1, 0]
err = 0.0
for s in [0.3 + 0.1j, 0.05j, 1e-4 + 0.2j, 2.0 - 0.7j, 0.6j]:
    r1 = np.linalg.det(L0 + s*I + E)/np.linalg.det(L0 + s*I)
    r2 = q(s)/(s*(s - mu))
    err = max(err, abs(r1 - r2)/abs(r2))
print(f'(a) max erro relativo det(Lt)/det(L) vs q/(s(s-mu)): {err:.2e}')
print(f'    raizes de Lt (autovalores de -(L0+E)): {np.sort_complex(-np.linalg.eigvals(L0 + E)).round(6)}')
print(f'    raizes de q: {np.sort_complex(np.roots([1, d0*(1+eij[0,0]) - mu + d1*(1+eij[1,1]), q(0)])).round(6)}')
sv = np.linalg.svd(L0 + E, compute_uv=False)
print(f'(b) 0 simples: sigma_min(Lt(0)) = {sv[-1]:.3e} (invertivel)')
for caso in ('jordan', 'dupla'):
    L0c, x0c, x1c, Sc = montar(caso)
    Ec, _ = deflacao(L0c, x0c, x1c, Sc)
    svc = np.linalg.svd(L0c + Ec, compute_uv=False)
    svL = np.linalg.svd(L0c, compute_uv=False)
    print(f'(c) 0 com multiplicidade 2 ({caso}): sigma_min(Lt(0)) = {svc[-1]:.3e}  (sigma_min(L(0)) = {svL[-1]:.1e}): '
          f'{"SINGULAR, como previsto" if svc[-1] < 1e-8 else "INVERTIVEL (contradiz)"}')
for eta in (1e-3, 1e-6):
    Ee, eij2 = deflacao(L0, x0, x1, S, eta=eta)
    q2 = lambda s: (s + d0*(1 + eij2[0, 0]))*(s - mu + d1*(1 + eij2[1, 1])) - d0*d1*eij2[0, 1]*eij2[1, 0]
    devs = []
    for s in (1e-2j, 1e-3j, 1e-4j):
        r1 = np.linalg.det(L0 + s*I + Ee)/np.linalg.det(L0 + s*I)
        devs.append(abs(r1*s*(s - mu) - q2(s)))
    print(f'(d) vetor de fase inexato (|eta| ~ {eta:.0e}): |det Lt/det L * s(s-mu) - q(s)| em s = 1e-2i, 1e-3i, 1e-4i: '
          + ', '.join(f'{x:.1e}' for x in devs) + '  (nao some quando o vetor e inexato; some com eta = 0)')
