from sage.all import *
from bijectivity_exponential_maps.conditions_bijectivity import cond_faces, nondegenerate, nondeg_cond1, nondeg_cond2
import unittest

class Tests(unittest.TestCase):
    def test_cond_faces(self):
        W = matrix([[1,0,-1,0],[0,1,0,-1]]).right_kernel_matrix()
        Wt = matrix([[1,0,-1,1],[0,1,-1,0]]).right_kernel_matrix()
        self.assertTrue(cond_faces(W, Wt))
    
    def test_nondeg_cond1(self):
        W = matrix([[-4,2,-7,1],[-9,-1,-1,-1],[-1,0,-1,1]]).right_kernel_matrix()
        Wt = matrix([[-5,-1,2,2],[1,0,2,21],[-2,0,0,2]]).right_kernel_matrix()

        self.assertFalse(nondeg_cond1(W, Wt))

        W = matrix([[-4,2,-7,1],[-9,-1,-1,-1],[-1,0,-1,1]]).right_kernel_matrix()
        Wt = matrix([[-5,1,-2,2],[1,0,-2,21],[-2,0,0,2]]).right_kernel_matrix()
        self.assertTrue(nondeg_cond1(W, Wt))

        W = matrix(2,4,[1,2,0,0,0,0,5,1]).right_kernel_matrix()
        A = matrix([[1,1,2,2],[1,0,1,-1]]).right_kernel_matrix()
        Wt = A.right_kernel_matrix()

        self.assertFalse(nondeg_cond1(W, Wt))

        W = matrix(3,5,[1,2,0,0,0,0,0,5,1,0,0,0,0,0,1]).right_kernel_matrix()
        A = matrix([[1,1,-2,-2,2],[1,0,1,-1,2]]).right_kernel_matrix()
        Wt = A.right_kernel_matrix()

        self.assertFalse(nondeg_cond1(W, Wt))

        A = matrix([[1, 0, 0, 1], [0, 1, 1, -1]]).right_kernel_matrix()
        B = matrix([[1,1,0,0],[0,0,1,1]])

        self.assertTrue(nondeg_cond1(B, A))

# Todo: test nondeg_cond2

if __name__ == '__main__':
    unittest.main()
