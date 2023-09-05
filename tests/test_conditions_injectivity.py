from sage.all import *
from sign_vector_conditions.conditions_injectivity import condition_inj_intersection, condition_inj_minors, geq, leq, geq_leq
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
                self.assertEqual(condition_inj_intersection(A, B), condition_inj_minors(A, B))

if __name__ == '__main__':
    unittest.main()
