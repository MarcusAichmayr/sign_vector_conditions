from sage.all import *
from bijectivity_exponential_maps.conditions_bijectivity_robust import cond_closure_sign_vectors, cond_closure_minors
import unittest

class Tests(unittest.TestCase):
    def test_conditions_bijectivity_robust_equivalence_random(self):
        m = 3
        n = 4

        for i in range(50):
            A = random_matrix(ZZ, m, n)
            B = random_matrix(ZZ, m, n)
            if A.rank() != 3 or B.rank() != 3:
                print('not full rank')
            else:
                self.assertEqual(cond_closure_sign_vectors(A, B), cond_closure_minors(A, B))

if __name__ == '__main__':
    unittest.main()