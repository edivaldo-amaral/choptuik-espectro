from fractions import Fraction as Q
import unittest

from component_exterior import (field_expression, field_norm, component_majorant,
                               gamma1_majorant, weighted_bound, propose_weights)
from sharp_propagation import SOURCE, X, TAU, gamma1, odd


class ComponentExteriorTests(unittest.TestCase):
    def test_combine_signs_before_absolute_values(self):
        fields = [{(0, 1, 0): Q(1)}, {(0, 1, 0): Q(1)}, {}, {}]
        self.assertEqual(field_norm(field_expression('A(t,v1,1); A(t,v2,-1);', fields)), 0)
        self.assertEqual(field_expression('A(t,v1,1); P(t);', fields), {(0, 1, 0): -1})

    def test_parity_projection_and_multiplicity(self):
        field = {(0, 0, 0): Q(2), (0, 1, 0): Q(3)}
        self.assertEqual(field_norm(field, 0), 2)
        self.assertEqual(field_norm(field, 1), Q(15, 2))

    @unittest.skipUnless((SOURCE/'GAMMA2_macro.c').exists(), 'fonte RT nao disponivel')
    def test_zero_reference(self):
        self.assertEqual(component_majorant([{}, {}, {}, {}]), [[0]*4 for _ in range(4)])

    def test_gamma1_first_component_vanishes(self):
        h = [X*(1+TAU), 1+X, X*X, TAU+X]
        selected = [odd(gamma1(h)[0].div_x()), gamma1(h)[1].div_x(),
                    gamma1(h)[3].div_x(), gamma1(h)[4].div_x()]
        self.assertEqual(selected[0], 0)
        self.assertEqual(selected[1], h[0].div_x()+(h[1]-h[1].parity()).div_x())
        self.assertEqual(selected[2], -h[0].div_x())
        self.assertEqual(weighted_bound(gamma1_majorant(), [Q(1)]*4), Q(68, 45))

    def test_weights_validated_by_columns(self):
        matrix = [[Q(0), Q(4)], [Q(1, 4), Q(0)]]
        weights = propose_weights(matrix)
        self.assertLess(weighted_bound(matrix, weights), Q(11, 10))
        self.assertEqual(weighted_bound(matrix, [Q(1)]*2), 4)
        with self.assertRaises(ValueError):
            weighted_bound(matrix, [Q(0), Q(1)])


if __name__ == '__main__':
    unittest.main()
