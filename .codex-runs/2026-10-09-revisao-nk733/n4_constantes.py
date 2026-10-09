#!/usr/bin/env python3
"""N4: eps_mu_L e dB_L(banda 16) refeitos para este lambda0 (rodar da raiz do projeto, CHOPTUIK_EPS_FUNDO_L=4.7e-08)."""
import sys, json; sys.path.insert(0, 'scripts')
import fundo_L, signed_operator_L as sl
z = json.load(open('build/nk733/Z.json'))
for s in (0.7331796596606924, 0.401024733, 0.16830707):
    print(s, repr(fundo_L.eps_mu_L(complex(s))))
print('Z.json eps_mu', repr(z['eps_mu']))
g = sl.campos(); print('dB_L banda 16', fundo_L.dB_L(g, 16)); print('Z.json', z['dB_L'])
