r"""
In this module, we check robust bijectivity for exponential maps.

EXAMPLES::

    sage: from bijectivity_exponential_maps import *
    sage: W = matrix([[1,0,-1,0],[0,1,0,0]])
    sage: W
    [ 1  0 -1  0]
    [ 0  1  0  0]
    sage: Wt = matrix([[1,0,-1,0],[0,1,-1,1]])
    sage: Wt
    [ 1  0 -1  0]
    [ 0  1 -1  1]
    sage: from sign_vectors.oriented_matroids import *
    sage: topes_from_matrix(W, kernel=True)
    [(+0+-), (-0--), (-0-+), (+0++)]
    sage: topes_from_matrix(Wt, kernel=True)
    [(++++), (+-++), (---+), (+++-), (-+--), (----)]
    sage: cond_closure_sign_vectors(W, Wt)
    True
    sage: W.minors(2)
    [1, 0, 0, 1, 0, 0]
    sage: Wt.minors(2)
    [1, -1, 1, 1, 0, -1]
    sage: cond_closure_minors(W, Wt)
    True
    sage: var('a,b,c')
    (a, b, c)
    sage: W = matrix([[1,0,-1],[0,c,-1]])
    sage: W
    [ 1  0 -1]
    [ 0  c -1]
    sage: Wt = matrix([[1,0,a],[0,1,b]])
    sage: Wt
    [1 0 a]
    [0 1 b]
    sage: W.minors(2)
    [c, -1, c]
    sage: Wt.minors(2)
    [1, b, -a]
    sage: cond_closure_minors(W, Wt)
    [[-b > 0, [c == 0, 'or', c > 0], [c == 0, 'or', -a*c > 0]], 'or', [-b < 0, [c == 0, 'or', c < 0], [c == 0, 'or', -a*c < 0]]]
    sage: cond_closure_minors(W, Wt(b=0))
    False
    sage: cond_closure_minors(W, Wt(b=2))
    [[c == 0, 'or', c < 0], [c == 0, 'or', -a*c < 0]]
    sage: cond_closure_minors(W(c=1), Wt(b=2))
    False
    sage: cond_closure_minors(W(c=0), Wt(b=2))
    True
    sage: cond_closure_minors(W(c=-1), Wt(b=2))
    [a < 0]
"""

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

from sign_vectors.oriented_matroids import topes_from_matrix

from sage.functions.generalized import sign

from sage.rings.real_mpfr import RR # used for casting
from sage.rings.integer_ring import ZZ

# Cor33cond2
def cond_closure_sign_vectors(W, Wt):
    r"""
    Returns whether the oriented matroid corresponding to ``W`` is a subset of
    the closure of the oriented matroid corresponding to ``Wt``.

    INPUT:
    
    - ``W`` -- a matrix with ``n`` columns
    
    - ``Wt`` -- a matrix with ``n`` columns
    
    OUTPUT:
    a boolean
    """
    tW  = topes_from_matrix(W, kernel=True)
    tWt = topes_from_matrix(Wt, kernel=True)
    tWn = normalize(tW) # 0++0
#    tWtn = normalize(tWt) # +--+ Do not normalize second list of Topes.
#    print('tWn:', tWn)
#    print('tWtn:', tWtn)
    for X in tWn:
        val = True
        for Y in tWt:
            if X <= Y:
                val = False
                break
        if val:
            return False
    return True

# Cor33cond3
def cond_closure_minors(W, Wt):
    r"""
    Returns whether for each non-zero maximal minor ``m`` of ``W``, we have
    either ``m mt > 0`` for each corresponding maximal minor ``mt`` of ``Wt`` or
    ``m mt < 0`` for each corresponding maximal minor ``mt`` of ``Wt``.
    
    INPUT:
    
    - ``W`` -- a matrix
    
    - ``Wt`` -- a matrix with the same dimensions as ``W``
    
    OUTPUT:
    A boolean or a symbolic expression if variables occur.
    """
    if W.dimensions() != Wt.dimensions():
        raise ValueError('Matrices must have same dimensions.')
    d  = W.nrows()
    m  = W.minors(d)
    mt = Wt.minors(d)
#    print('minors of W: ', m)
#    print('minors of Wt:', mt)
    
    def eq(a):
        r"""
        Compares an expression with ``0``.
        
        INPUT:
        
        - ``a`` -- a real number or symbolic expression
        
        OUTPUT:
        
        - ``True``  -- The number ``a`` is equal to zero.
        
        - ``False`` -- The number ``a`` is a non-zero real number.
        
        - ``None``  -- The number ``a`` is a symbolic expression.
        """
        try:
            if RR(a) == 0:
                return True # equal to zero
            else:
                return False # non-zero real number
        except TypeError:
            return None # symbolic expression

    s = 0 # will be set to the sign (1 or -1) of m[i] mt[i]
    L = [] # will be set to a list of symbolic entries. Depending on s, each entry is ``> 0`` or ``< 0``.
    L_eq = [] # will be set to a list of equalities.
    # L_left, L_right will be the ``or`` pairs
    L_left = [] # will be set to a list of equalities.
    L_right = [] # will be set list of symbolic entries. Depending on s, each entry is ``> 0`` or ``< 0``.
    
    for i in range(len(m)):
        if eq(m[i]) == None: # mi symbolic
            if eq(mt[i]) == True: # mti zero, hence, product zero and cannot be > or < 0
                L_eq.append(m[i] == 0)
            else:
                L_left.append(m[i] == 0)
                L_right.append(m[i]*mt[i])
        #if eq(m[i]) != None:
        elif eq(m[i]) == False: # mi non-zero, not symbolic
            if eq(mt[i]) == None: # mti symbolic
                L.append(m[i]*mt[i]) # needs > or < later
            elif eq(mt[i]) == True: # mti zero
                return False
            elif s == 0: # mti non-zero, not symbolic
                s = ZZ(sign(m[i]*mt[i]))
            else:
                if (s*m[i]*mt[i]) <= 0: # cast
                    return False
    
    # construction of the output list
    if L == [] and L_eq == [] and L_left == []: # hence, L_right is also []
        return True
    if L == [] and L_left == []:
        return L_eq
    if s == 0:
        if L != []:
            E_eq_left = []
            E_eq_right = []
            for l in L:
                E_eq_left.append(l > 0)
                E_eq_right.append(l < 0)
                
        if L_left != []:
            E_or_left = []
            E_or_right = []
            for i in range(len(L_left)):
                E_or_left.append([L_left[i], 'or', L_right[i] > 0])
                E_or_right.append([L_left[i], 'or', L_right[i] < 0])
        
        if L == []:
            if len(E_or_left) == 1: # hence, also len(E_or_right) == 1
                out = [flatten(E_or_left), 'or', flatten(E_or_right)]
            else:
                out = [E_or_left, 'or', E_or_right]
        elif L_left == []:
            out = [E_eq_left, 'or', E_eq_right]
        else:
            out = [E_eq_left + E_or_left, 'or', E_eq_right + E_or_right]
    else:
        if s > 0:
            def rel(a):
                return a > 0
        else:
            def rel(a):
                return a < 0
        
        if L != []:
            E_eq = []
            for l in L:
                E_eq.append(rel(l))
        if L_left != []:
            E_or = []
            for i in range(len(L_left)):
                E_or.append([L_left[i], 'or', rel(L_right[i])])
        if L == []:
            if len(E_or) == 1:
                out = flatten(E_or)
            else:
                out = E_or
        elif L_left ==[]:
            out = E_eq
        else:
            out = E_eq + E_or
    
    if L_eq == []:
        return out
    else:
        return L_eq + [out]
