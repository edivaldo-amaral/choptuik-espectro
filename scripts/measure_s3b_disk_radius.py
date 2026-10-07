#!/usr/bin/env python3
"""DIAGNOSTICO EM PONTO FLUTUANTE -- NAO E PROVA. Escolha do raio do disco D de S3b em
torno da raiz 0,401 (setor A).

(1) raizes de L (operador selecionado, exportador de RT) perto de 0,401, em varias malhas;
(2) max ||V(s)|| = ||J (op + s)^-1||_w (norma que dimensiona a rota certificada de L) na
    circunferencia |s - 0,401| = r, metade superior (conjugacao), para L e para K.
"""
import sys
import numpy as np
sys.path.insert(0, 'scripts')
from measure_sharp_vs_selected_V import selecionado, sharp, norma_V


def main():
    for malha in ('8x24', '10x30', '11x33', '12x36'):
        A, _, _ = selecionado(malha)
        s = -np.linalg.eigvals(A)
        d = np.abs(s - 0.401)
        viz = [f'{s[i].real:+.4f}{s[i].imag:+.4f}i' for i in np.argsort(d)[:4]]
        print(f'raizes de L {malha}: ' + ', '.join(viz))
    angs = np.linspace(0, np.pi, 49)
    raios = (0.02, 0.05, 0.08, 0.10, 0.11, 0.125, 0.135, 0.15)
    for nome, dados in (('L 11x33', selecionado('11x33')), ('L 12x36', selecionado('12x36')), ('K 12x36', sharp('12x36'))):
        A, J, w = dados
        print(f'{nome}: ' + '  '.join(
            f'r {r:.3f}: {max(norma_V(A, J, w, 0.401 + r*np.exp(1j*t)) for t in angs):.1f}' for r in raios))


if __name__ == '__main__':
    main()
