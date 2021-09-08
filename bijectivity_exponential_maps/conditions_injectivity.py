#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from .utility import normalize
from sign_vectors.oriented_matroids import covectors_from_matrix

from sage.rings.real_mpfr import RR # used for casting

def cond_inj_intersection(W, Wt):
    r"""
    Returns whether the intersection of the oriented matroids corresponding to
    ``W`` and ``right_kernel(Wt)`` consists of the zero sign vector only.
    
    INPUT:
    
    - ``W`` -- a matrix with ``n`` columns
    
    - ``Wt`` -- a matrix with ``n`` columns
    
    OUTPUT:
    a boolean
    """
    if W.ncols() != Wt.ncols():
        raise ValueError('Matrices have different number of columns.')
    
    cvW = covectors_from_matrix(W)
    try:
        cvWt = covectors_from_matrix(Wt, kernel=True)
    except ValueError: # kernel matrix might be empty, no elementary vectors
        return True
    cvWn = normalize(cvW)
    cvWtn = normalize(cvWt)
    for X in cvWtn[1:]: # first element is zero vector
        if X in cvWn[1:]: # first element is zero vector
            return False
    return True

def max_minors_prod(A, B):
    r"""Multiplies the maximal minors of two matrices component-wise."""
    A1 = A.matrix_from_rows(A.pivot_rows())
    B1 = B.matrix_from_rows(B.pivot_rows())
    if A1.dimensions() != B1.dimensions():
        raise ValueError('Matrices must have same dimensions.')
    r = A1.nrows() # should be equal to B1.nrows()
    mA = A1.minors(r)
    mB = B1.minors(r)
    
    return [mA[i]*mB[i] for i in range(len(mA))]

def compare_all(v, rel):
    r"""Checks whether all entries satisfy the relation."""
    l = []
    for a in v:
        if rel(a) == False:
#            print('system cannot be satisfied')
            return False
        # here a >= 0 is either True or there are indeterminants
        elif rel(a) != True:
            l.append(rel(a))
    if len(l) == 0:
##        print('all entries are non-negative')
        return True
    else:
#        print('all conditions of this list should hold:')
        return l


def geq(v):
    r"""Checks whether all entries are non-negative."""
    def rel(a):
        try:
            return RR(a) >= 0 # cast to a real number
        except TypeError:
            return a >= 0
    return compare_all(v, rel)

def leq(v):
    r"""Checks whether all entries are non-positive."""
    def rel(a):
        try:
            return RR(a) <= 0 # cast to a real number
        except TypeError:
            return a <= 0
    return compare_all(v, rel)

def geq_leq(v):
    r"""
    Returns true if either each element of ``v`` is greater than or equal to zero or less than or equal to zero.
    Supports symbolic expressions.
    
    INPUT:
    
    - ``v`` -- a list of numbers or symbolic expressions.
    
    OUTPUT:
    Depending on the input, the output can have several appearances:
    
    - a boolean
        
        - if true, no symbolic expressions have occured and either each entry of ``v`` is greater or equal zero
          or each entry is less or equal zero.
          
        - if false, there are entries with opposing signs or every entry is zero.
    
    - a list of inequalities
    
        - if all inequalities are satisfied,
          then either each element of ``v`` is greater or equal zero or less or equal zero.
    
    - a list of two lists of inequalities
    
        - if the inequalities of exactly one of these lists are satisfied,
          then either each element of ``v`` is greater or equal zero or less or equal zero.
    """
    ge = geq(v)
    le = leq(v)
    
    if ge == True and le == True: # all entries are zero
        return False
    elif ge == False and le == False: # mixed signs
        return False
    elif ge == False:
        if le == True:
            return True
        else:
#            print('all of these conditions have to be satisfied:')
            return le
    elif le == False:
        if ge == True:
            return True
        else:
#            print('all of these conditions have to be satisfied 2:')
            return ge
    else:
#        print('either all of the first or all of the second conditions have to be satisfied:')
        return [ge, le]

# Corollary 4 condition 2
def cond_inj_minors(W, Wt):
    r"""
    Returns whether the products of the corresponding maximal minors of the
    matrices ``W`` and ``Wt`` have the same sign.
    
    INPUT:
    
    - ``W`` -- a matrix
    
    - ``Wt`` -- a matrix with the same dimensions as ``W``
    
    OUTPUT:
    A boolean or a symbolic expression if variables occur.
    """
    m = max_minors_prod(W, Wt)
#    print(m)
    return geq_leq(m)

