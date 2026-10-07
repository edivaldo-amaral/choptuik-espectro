from fractions import Fraction as Q
import unittest

from sharp_shell_remainder import (schur_factors, remainder, required_degree,
                                   close_shell)
from build_tile import exact_inverse


class ShellRemainderTests(unittest.TestCase):
    def test_three_block_schur_against_full_exact_elimination(self):
        # Ordem F,U,T; sinais diferentes exercitam ambos os acoplamentos.
        ff, fu, ft = Q(2), Q(1,5), Q(-1,7)
        uf, uu, ut = Q(-1,3), Q(3), Q(1,4)
        tf, tu, tt = Q(1,6), Q(-1,8), Q(4,5)
        inverse_z = exact_inverse([[ff,ft],[tf,tt]])
        exact = uu-sum(a*inverse_z[i][j]*b for i,a in enumerate((uf,ut))
                       for j,b in enumerate((fu,tu)))
        a = 1/ff
        rt = tt-tf*a*ft
        left, right = ut-uf*a*ft, tu-tf*a*fu
        lb, rb = schur_factors(abs(a),abs(ft),abs(tf),abs(uf),abs(fu),abs(ut),abs(tu))
        for degree in (0,1,5,12):
            inverse_approx = sum((1-rt)**j for j in range(degree+1))
            approx = uu-uf*a*fu-left*inverse_approx*right
            self.assertLessEqual(abs(exact-approx), remainder(abs(rt-1),lb,rb,degree))

    def test_budget_includes_both_couplings_and_candidate_norm(self):
        degree = required_degree(Q(9,10), 20, 30, Q(1,10**8))
        self.assertLessEqual(remainder(Q(9,10),20,30,degree), Q(1,10**8))
        self.assertGreater(remainder(Q(9,10),20,30,degree-1), Q(1,10**8))
        self.assertFalse(close_shell(100,Q(1,2),Q(1,100))['accepted'])
        self.assertTrue(close_shell(100,Q(1,2),Q(1,1000))['accepted'])
        self.assertFalse(close_shell(100,Q(1,2),Q(1,200))['accepted'])

    def test_invalid_or_noncontractive_inputs(self):
        for delta in (Q(-1,10), Q(1), Q(2)):
            with self.assertRaises(ValueError):
                required_degree(delta,1,1,Q(1,100))
        with self.assertRaises(ValueError):
            required_degree(Q(1,2),1,1,0)
        self.assertEqual(required_degree(0,2,3,Q(1,100)),0)


if __name__ == '__main__':
    unittest.main()
