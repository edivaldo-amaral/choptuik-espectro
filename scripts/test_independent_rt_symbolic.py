"""Regressao da verificacao simbolica independente (auditoria fisica 16/09).

Alem de exigir que as identidades valham, exige que MUTACOES plausiveis
falhem: isso mostra que o verificador nao e vacuo (ao contrario de testes que
conferem uma formula contra ela mesma).
"""
import unittest

import independent_rt_symbolic as m


def fails(fn, degrees):
    for deg, order in degrees:
        m.CAP[0] = deg
        if not m.all_zero(fn(order)):
            return True
    return False


class UniversalIdentities(unittest.TestCase):
    def test_all_identities_hold_universally(self):
        for name, fn, degrees in m.CHECKS:
            with self.subTest(name=name):
                self.assertFalse(fails(fn, degrees), name)

    def test_macro_tables_agree_with_tex(self):
        self.assertEqual(m.check_macro_tables(), 0)


class Mutations(unittest.TestCase):
    def test_background_correction_is_needed(self):
        def mutated(order):
            src = m.SymbolSource()
            h, w = src.multifield(order), src.multifield(order)
            F, C = m.select(m.lin_omega_fc(w, h, m.S))
            return [a - b for a, b in zip(m.K_doc(w, C, m.S), m.M_doc(w, F, m.S))]
        self.assertTrue(fails(mutated, [(2, 1)]))

    def test_factor_one_gamma2_in_linearization_fails(self):
        def mutated(order):
            src = m.SymbolSource()
            w = src.multifield(order)
            h = [m.muZ1(p) for p in w]
            lin = m.omega_fc(h, m.MU, lin_only=True)
            lhs = [a + m.X * b for a, b in zip(lin, m.gamma2_tex(w, h))]
            return [a - m.muZ1(b) for a, b in zip(lhs, m.omega_fc(w))]
        self.assertTrue(fails(mutated, [(2, 1)]))

    def test_wrong_exponent_for_translation_fails(self):
        def mutated(order):
            src = m.SymbolSource()
            w = src.multifield(order)
            h = [m.muZ1(p) for p in w]
            lhs = m.lin_omega_fc(w, h, m.S)  # s generico em vez de s=mu
            return [a - m.muZ1(b) for a, b in zip(lhs, m.omega_fc(w))]
        self.assertTrue(fails(mutated, [(1, 2)]))

    def test_sharp_coefficient_mutation_fails(self):
        original = m.sharp_G_tex

        def mutated_G(w, Om):
            out = original(w, Om)
            return [out[0], out[1] + m.half(w[1] * m.half(Om[1] + Om[1].P()))]
        try:
            m.sharp_G_tex = mutated_G
            self.assertTrue(fails(m.check_rt_sharp_identity, [(2, 1)]))
        finally:
            m.sharp_G_tex = original

    def test_sign_of_M_source_matters(self):
        def mutated(order):
            src = m.SymbolSource()
            h, w = src.multifield(order), src.multifield(order)
            F, C = m.select(m.lin_omega_fc(w, h, m.S))
            e = [-p.P() for p in m.omega_fc(w)]
            corr = m.sharp_G_tex(h, e)
            return [a + b + c for a, b, c in zip(m.K_doc(w, C, m.S), m.M_doc(w, F, m.S), corr)]
        self.assertTrue(fails(mutated, [(1, 2)]))


if __name__ == '__main__':
    unittest.main()
