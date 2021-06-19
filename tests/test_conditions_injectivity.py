from sage.all import *
from bijectivity_exponential_maps.conditions_injectivity import cond_inj_intersection, cond_inj_minors, geq, leq, geq_leq
import unittest

class Tests(unittest.TestCase):
    def test_conditions_injtivity_equivalence_random(self):
        m = 3
        n = 4

        for i in range(50):
            A = random_matrix(ZZ,m,n)
            B = random_matrix(ZZ,m,n)
            if A.rank() != B.rank():
                print('different rank')
            else:
                self.assertEqual(cond_inj_intersection(A, B), cond_inj_minors(A, B))

    def test_cond_inj_intersection_empty(self):
        A = identity_matrix(3)
        B = A # kernel of B is empty
        self.assertTrue(cond_inj_intersection(A, B))

    def test_geq(self):
        l = [0, 5, 1]
        self.assertTrue(geq(l))
        l = [0, 0]
        self.assertTrue(geq(l))
        l = [0, -5]
        self.assertFalse(geq(l))
        var('x')
        l = [x, x**2 + 1, -1, 5]
        self.assertFalse(geq(l))

    def test_geq_leq_symbolic(self):
        var('x')
        l = [x, x**2 + 1, -1, 5]
        self.assertFalse(geq_leq(l))

if __name__ == '__main__':
    unittest.main()
