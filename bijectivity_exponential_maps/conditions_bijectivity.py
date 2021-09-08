#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from .utility import normalize, pos_cocircuits_from_matrix, pos_covectors_from_matrix
from elementary_vectors import elementary_vectors, exists_vector
from sign_vectors.oriented_matroids import cocircuits_from_matrix

from sage.modules.free_module_element import vector, zero_vector
from sage.matrix.constructor import matrix
from sage.rings.infinity import Infinity
from sage.misc.flatten import flatten
from sage.geometry.polyhedron.constructor import Polyhedron

def cond_faces(W, Wt):
    r"""
    Returns whether every positive sign vector ``X`` corresponding to the rows of
    ``Wt`` has a positive sign vector ``Y`` corresponding to the rows of ``W``
    such that ``Y <= X``.
    
    INPUT:
    
    - ``W`` -- a matrix with ``n`` columns
    
    - ``Wt`` -- a matrix with ``n`` columns
    
    OUTPUT:
    a boolean
    """
    ccWp  = pos_cocircuits_from_matrix(W,  kernel=False)
    ccWtp = pos_cocircuits_from_matrix(Wt, kernel=False)

    for X in ccWtp:
        val = True
        for Y in ccWp:
            if Y <= X:
                val = False
                break
        if val:
            return False
    return True


def cond_nondegenerate(W, Wt, certificate=False):
    r"""
    Let ``S``, ``St`` be the sub spaces corresponding to the matrices ``W``,
    ``Wt``. Returns whether the pair ``(S, St)`` is non-degenerate.
    
    INPUT:
    
    - ``W`` -- a matrix with ``n`` columns
    
    - ``Wt`` -- a matrix with ``n`` columns
    
    - ``certificate`` -- a boolean (default: ``False``)

    OUTPUT:
    a boolean
    
    - If ``certificate`` is true:
    
      - If the result is true, returns a vector as a certificate.
    """
    return nondegenerate(W, Wt, certificate=certificate)


def nondegenerate(W, Wt, certificate=False):
    r"""
    Let ``S``, ``St`` be the sub spaces corresponding to the matrices ``W``,
    ``Wt``. Returns whether the pair ``(S, St)`` is non-degenerate.
    
    INPUT:
    
    - ``W`` -- a matrix with ``n`` columns
    
    - ``Wt`` -- a matrix with ``n`` columns
    
    - ``certificate`` -- a boolean (default: ``False``)

    OUTPUT:
    a boolean

    - If ``certificate`` is true:
    
      - If the result is true, returns a vector as a certificate.
    """
    c2 = nondeg_cond2(W, Wt)
    if c2:
        return True
    
    c1 = nondeg_cond1(W, Wt, certificate=certificate)
    return c1


def nondeg_cond1(W, Wt, certificate=False):
    r"""
    Returns whether the first condition of ``cond_nondegenerate`` is satisfied.
    
    INPUT:
    
    - ``W`` -- a matrix with ``n`` columns
    
    - ``Wt`` -- a matrix with ``n`` columns
    
    - ``certificate`` -- a boolean (default: ``False``)
    
    OUTPUT:
    a boolean

    - If ``certificate`` is true:
    
      - If the result is true, returns a vector as a certificate.
    """
    
    if W.ncols() != Wt.ncols():
        raise ValueError('Matrices have different number of columns.')
    # eventuell auf disjoint support achten. Sind +0 und 0+ dabei, dann müssen wir ++ nicht betrachten.
    P = pos_covectors_from_matrix(W, kernel=True)
    
    if P == []:
        return True
    
    n = Wt.ncols() # length of vectors
    degenerate = False # might change in recursion
    
    L  = [-Infinity for i in range(n)]
    R = [0   for i in range(n)]
    inf = [Infinity  for i in range(n)] # does not change
    
#    l = [True for i in range(n)]
    proof = []
    
    M = Wt.right_kernel_matrix()
    
    def rec(P, M, I, L, R):
        r"""
        recursive function
        
        INPUT:
        
        - ``P`` -- a list of positive sign vectors
        
        - ``M`` -- a matrix
        
        - ``I`` -- a list of indices
        
        - ``L`` -- a list of values ``-Infinity`` and ``1``

        - ``R`` -- a list of values ``0`` and ``Infinity``
        """
        nonlocal degenerate # use nonlocal to use the same ``degenerate``
        nonlocal proof
        
        while P:
            X = P.pop()
            if set(flatten(I)).issubset(X.zero_support()): # P[0].positive_support() subseteq J, insbesondere muss X ein + auf J haben
                L_ = L[:]
                R_ = R[:]
                for i in X.support():
                    L_[i]  = 1
#                     L_[i]  = 0 # hier sollte auch 1 gehen und abgeschlossene Intervallhälfte.
#                     l_[i]  = False # open interval
                    R_[i] = Infinity
                
                M_ = equal_components(M, X.support()) # make to echelon form?
                A = M_.right_kernel_matrix()
                
                if A: # A is not empty matrix
                    evs = elementary_vectors(A, kernel=True)
                
                    if exists_vector(evs, L_, R_):
                        degenerate = True

                        I = I + [X.support()]
                        if certificate:
                            proof = [I, find_vector(Wt, I)]
                        return
                    elif exists_vector(evs, L_, inf):
                        rec(P[:], M_, I + [X.support()], L_, R_)
                    else:
                        if certificate:
                            v = exists_vector(A, L_, inf, kernel=False, certificate=True)[1]
                            proof.append([v, I + [X.support()]])
                else:
                    if certificate:
                        v = exists_vector(A, L_, inf, kernel=False, certificate=True)[1]
                        proof.append([v, I + [X.support()]])

            if degenerate:
                return
        return
    
    rec(P, M, [], L, R)
    
    if certificate:
        return [not degenerate, proof]
    return not degenerate


def equal_components(M, I):
    r"""
    The kernel matrix has equal components on ``I``.
    
    INPUT:
    
    - ``M`` -- a matrix
    
    - ``I`` -- a list of indices
    """
    
    n = M.ncols()
    
    if len(I) >= 2:
        i = I[0]
        v = zero_vector(n)
        v[i] = 1
        for j in I[1:]:
            w = v[:]
            w[j] = -1
            M = matrix(list(M) + [w])
    return M


def find_vector(M, I):
    r"""
    Returns a vector
    
    INPUT:
    
    - ``M`` -- a matrix
    
    - ``I`` -- a list of lists of indices
    
    OUTPUT:
    
    Returns a vector $v$ in the row space of ``M`` such that:
    
    - For each $I_k$ in ``I``, the entries $v_i$ with $i \in I_k$ are equal and positive.
    
    - The remaining entries are zero or negative.
    """
    n = M.ncols()
    ieqs = []
    eqns = []
    A = M.right_kernel_matrix()
    # lie in rowspace of M
    for v in A.rows():
        eqns.append([0] + list(v))
        # subspace with equal components

    # components need to be equal that is opposite sign for kernel
    for Ik in I:
        if len(Ik) >= 2:
            i = Ik[0]
            v = zero_vector(n+1)
            v[i+1] = 1
            for j in Ik[1:]:
                w = v[:]
                w[j+1] = -1
                eqns.append(w)
    
    J = flatten(I)
    for i in range(n):
        v = zero_vector(n+1)
        if i in J:
            v[0] = -1 # component - 1 >= 0 # instead of > 0
            v[i+1] = 1
            ieqs.append(v)
        else:
            v[i+1] = -1 # component <= 0
            ieqs.append(v)
    
#     print('ieqs:', ieqs)
#     print('eqns:', eqns)
    P = Polyhedron(ieqs=ieqs, eqns=eqns)
    return P.representative_point()


def nondeg_cond2(W, Wt):
    ccWt = normalize(cocircuits_from_matrix(Wt))
    ccWp = pos_cocircuits_from_matrix(W)
    for X in ccWt:
        val = False
        for Y in ccWp:
            if set(X.zero_support()).issubset(Y.zero_support()):
                val = True
                break
        if val == False:
            return False
    return True
