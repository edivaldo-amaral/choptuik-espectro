#!/usr/bin/env python3
"""R6 (revisao S4): a simetria de aliases de Floquet do operador L na base com sinal.
J(s) no modo m depende de s + i m/2; os blocos de B dependem so de d = mo - mi e da paridade de mi.
Logo o deslocamento m -> m + 2 (preserva o setor: m par <-> componentes 0,1,2; m impar <-> 3)
conjuga L(s + i) a L(s): o espectro e invariante por s -> s + i, e a similaridade e limitada nos pesos
kappa1^(2m) (fator constante kappa1^4). 'Numero de raizes em Re s > 0' so e finito numa faixa
fundamental de altura 1 em Im s."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
import signed_operator_L as sl
g = sl.campos(); mu = sl.MU; N = 30
out = []
for m in (-5, -2, 0, 3, 7):
    dJ = np.abs(sl.J_livre(sl.Radial(m, N), mu + 1j) - sl.J_livre(sl.Radial(m + 2, N), mu)).max()
    dB = max(np.abs(sl.bloco_B(g, sl.Radial(m + d, N), sl.Radial(m, N)) - sl.bloco_B(g, sl.Radial(m + 2 + d, N), sl.Radial(m + 2, N))).max() for d in (-3, 0, 2, 5))
    out.append(f'm = {m:3d}: max|J_m(mu + i) - J_(m+2)(mu)| = {dJ:.1e};  max_d |B(m+d <- m) - B(m+2+d <- m+2)| = {dB:.1e}')
out.append('=> sigma(L) invariante por s -> s + i; mu + i j (j inteiro) sao raizes de L com o mesmo autovetor deslocado.')
print('\n'.join(out)); Path(__file__).with_name('r6_saida.txt').write_text('\n'.join(out) + '\n')
