from sage.all import *
from bijectivity_exponential_maps.functions import *
import unittest

# Todo: more meaningful tests
class Tests(unittest.TestCase):
    def test_f_pol(self):
        W = matrix([[1,0,-1],[0,1,-1]])
        W
        Wt = matrix([[1,0,-1],[0,1,0]])
        Wt
        c = vector([1,2,4])
        fc = f_pol(W, Wt, c)
        self.assertEqual(fc(1,2), vector([-3, 0]))

        fc = f_pol(W, Wt) # without c
        self.assertEqual(fc(1,2), vector([0, 1]))

    # Todo: Tests for f_exp
    def test_f_pol(self):
        W = matrix([[1,0,-1],[0,1,-1]])
        W
        Wt = matrix([[1,0,-1],[0,1,0]])
        Wt
        c = vector([1,2,4])
        Fc = f_exp(W, Wt, c)
        Fc(1,2)

        Fc = f_exp(W, Wt) # without c
        Fc(1,2)

if __name__ == '__main__':
    unittest.main()
