"""Testes rapidos (<5 s) dos ingredientes analiticos de columnwise_exterior.py.

Os betas exatos (lentos) e os criterios 1-2 ficam em columnwise_exterior.py e
columnwise_validate.py; aqui so as cotas de Q contra colunas EXATAS em Fraction.
"""
from fractions import Fraction as Q
import unittest

import columnwise_exterior as CE
from columnwise_validate import q_column, P_BITS
from independent_free_inverse import exact_column, column_norm


def s_hat(kind, r, N):
    if kind == 'J':
        C = 1+2*CE.SQ2*r/(1-r)
    else:
        C = 1+2*CE.SQ2*Q(N, N-1)*r*r/(1-r*r)
    return (C+1/(4*C))/CE.MU


class ColumnwiseTests(unittest.TestCase):
    def test_radial_q_bound_dominates_exact_columns(self):
        for k2 in (Q(5, 4), Q(9, 8)):
            r = 1/k2
            for kind in ('J', 'K'):
                N = 9
                for m in (0, 1, 2, 3, 7, 30):
                    for n in (9, 11, 17):
                        if kind == 'K' and n % 2 == 0:
                            continue
                        exact = column_norm(exact_column(kind, m, n, CE.MU), k2)
                        self.assertLessEqual((n+1)*exact, s_hat(kind, r, N), (kind, m, n))

    def test_fourier_q_bound_dominates_exact_columns(self):
        for k2 in (Q(5, 4), Q(9, 8)):
            r = 1/k2
            for kind, k, nu in (('J', 1, Q(1)), ('K', 2, CE.NU2)):
                for m in (42, 43, 60):
                    for n in (0, 1, 3, 5, 9, 20):
                        if kind == 'K' and n % 2 == 0:
                            continue
                        exact = column_norm(exact_column(kind, m, n, CE.MU), k2)
                        bound = (CE.G0+CE.SQ2*nu*r**k/(1-r**k))/Q(m, 2)
                        self.assertLessEqual(exact, bound, (kind, m, n))

    def test_constants_are_upper_bounds(self):
        self.assertGreater(CE.SQ2**2, 2)
        self.assertGreater((2*CE.G0-1)**2, 2)
        self.assertGreater(CE.NU2**2*3, 4)

    def test_enclosure_contains_exact_column(self):
        for sysname, d, m, n in (('A', 1, 3, 12), ('A', 0, 2, 11), ('sharp', 0, 5, 9)):
            kind = CE.SYSTEMS[sysname]['kinds'][d]
            cols, _ = q_column(sysname, d, m, n, depth=50)
            exact = exact_column(kind, m, n, CE.MU)
            for j, ((xr, xi), rad) in cols.items():
                err = abs(exact[j].re-Q(xr, 1 << P_BITS))+abs(exact[j].im-Q(xi, 1 << P_BITS))
                self.assertLessEqual(err, 2*rad+Q(2, 1 << P_BITS), (sysname, j))


if __name__ == '__main__':
    unittest.main()


class ColumnwiseCertificateFiles(unittest.TestCase):
    """Reconfere, em racionais, os artefatos da Fase 2 e dos criterios 1-2 (leitura de disco, rapido)."""

    @classmethod
    def setUpClass(cls):
        import json
        root = CE.ROOT/'build/spectrum'
        cls.rep = json.loads((root/'columnwise-exterior.json').read_text())
        crit = root/'columnwise-criteria.json'
        cls.crit = json.loads(crit.read_text()) if crit.exists() else None
        cls.ck = {}
        for s, name in (('A', 'A-M107-D60.jsonl'), ('sharp', 'sharp-M107-D100.jsonl')):
            path = root/'columnwise-phase2'/name
            cls.ck[s] = [json.loads(l) for l in path.read_text().splitlines() if l.strip()] if path.exists() else None

    def test_band_contiguous_and_below_one(self):
        for s in ('A', 'sharp'):
            rows = self.ck[s]
            if rows is None:
                self.skipTest('checkpoint ausente')
            ns = sorted(r['n'] for r in rows)
            self.assertEqual(ns, list(range(384, 606)))
            self.assertTrue(all(Q(r['max_U']) < 1 for r in rows))
            self.assertEqual(Q(self.rep['phase2'][s]['band_max_U']), max(Q(r['max_U']) for r in rows))

    def test_d_star_is_max_of_three_bounds_and_rule(self):
        p2 = self.rep['phase2']
        for s in ('A', 'sharp'):
            e = p2[s]
            dstar = max(Q(e['band_max_U']), Q(e['d_R_far']), Q(e['d_F']))
            self.assertEqual(dstar, Q(e['d_star']))
            self.assertLess(dstar, 1)
            self.assertEqual(e['N_G_certified'], 384)
        self.assertLess(Q(p2['guarded_schur_kappa_crown_relatorio']['schur']), 1)
        self.assertLess(Q(p2['sharp_bridge']['delta']), 1)

    def test_criterion1_exact_on_sampled_exterior_columns(self):
        if self.crit is None:
            self.skipTest('rode scripts/columnwise_criteria.py')
        res = self.crit['results']
        self.assertGreaterEqual(len(res), 50)
        for r in res:
            self.assertLessEqual(Q(r['U']), Q(self.rep['phase2'][r['system']]['d_star']))
            self.assertTrue(r['m'] >= 107 or r['n'] >= 384)
        self.assertTrue(any(r['m'] == 0 for r in res) and any(r['m'] % 2 for r in res))
        self.assertTrue(any(r['d'] == 0 for r in res) and any(r['d'] > 0 for r in res))
