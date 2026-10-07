"""Testes rapidos (<5 s) dos recomputos independentes da auditoria de 16/09.

Os recomputos pesados (51 tiles, matrizes 6x18/12x36, beta de RefA) ficam nos
scripts independent_*.py; aqui so as partes baratas e as identidades-chave.
"""
from fractions import Fraction as Q
import unittest

from independent_free_inverse import (MU_REFA, exact_column, column_norm, tail_bound,
                                      sampled_exterior_max, radius_recovery, rt_K1, rt_K2)
from independent_crown_tiles import free_matrix, exact_inverse_columns
from independent_schur_scalars import crown_schur, sharp_bridge, KAPPA_CROWN_REPORT


def box_addresses(M, N):
    out = []
    for m in range(M):
        for n in range(N):
            for d in range(4):
                if m % 2 != (1 if d == 3 else 0) or (d == 0 and n % 2 == 0):
                    continue
                for part in ((0,) if m == 0 else (0, 1)):
                    out.append((d, part, m, n))
    return out


class IndependentAuditTests(unittest.TestCase):
    def test_free_inverse_is_exact_inverse_on_small_box(self):
        addresses = box_addresses(3, 5)
        J = free_matrix(addresses, MU_REFA)
        Qm = exact_inverse_columns(addresses, MU_REFA)
        n = len(addresses)
        for i in range(n):
            for j in range(n):
                self.assertEqual(sum((J[i][k]*Qm[k][j] for k in range(n)), Q(0)), 1 if i == j else 0)

    def test_column_zero_is_inverse_mu(self):
        self.assertEqual(column_norm(exact_column('J', 0, 0, MU_REFA), Q(5, 4)), 1/MU_REFA)

    def test_tail_values(self):
        self.assertEqual(tail_bound(176, 944, Q(5, 4)), Q(3, 35))
        self.assertEqual(tail_bound(6, 18, Q(5, 4)), Q(81, 19))
        self.assertEqual(tail_bound(1024, 4096, Q(9, 8)), Q(153, 4097))
        self.assertLessEqual(sampled_exterior_max(4, 12, MU_REFA, Q(5, 4), 2, 4),
                             tail_bound(4, 12, Q(5, 4)))

    def test_radius_recovery_declared_and_rt_only(self):
        result = radius_recovery()
        self.assertTrue(result['declared_matches_report'])
        self.assertTrue(result['rt_only_both_contract'])

    def test_rt_constants(self):
        self.assertEqual(rt_K1(Q(5, 4)), Q(80, 9))
        self.assertEqual(rt_K2(Q(1), Q(5, 4)), 18)
        self.assertEqual(rt_K2(Q(1), Q(9, 8)), 34)

    def test_schur_scalars(self):
        self.assertTrue(crown_schur(KAPPA_CROWN_REPORT)['schur_lt_0_997676'])
        kappa = Q(1146011126718314, 10**14)  # ~11.4601 (relatorio sharp precondicionado)
        self.assertEqual(sharp_bridge(kappa, (160, 960))['degree_1e12'], 437)


if __name__ == '__main__':
    unittest.main()
