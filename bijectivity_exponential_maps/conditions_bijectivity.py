r"""
In this module, we check bijectivity for exponential maps by verifying conditions from [MHR19]_.

EXAMPLES::

    sage: from bijectivity_exponential_maps import *

Let us consider the following matrices::

    sage: W = matrix([[1,0,1,0],[0,1,0,1]])
    sage: W
    [1 0 1 0]
    [0 1 0 1]
    sage: Wt = matrix([[1,0,0,-1],[0,1,1,1]])
    sage: Wt
    [ 1  0  0 -1]
    [ 0  1  1  1]

    sage: var('x1, x2, c1, c2, c3, c4')
    (x1, x2, c1, c2, c3, c4)
    sage: c = [c1, c2, c3, c4]

Therefore, we obtain the following exponential map::

    sage: Fc = f_exp(W, Wt, c)
    sage: Fc(x1, x2)
    (c1*e^x1 + c3*e^x2, c4*e^(-x1 + x2) + c2*e^x2)

Hence, we obtain the oriented matroids::

    sage: from sign_vectors.oriented_matroids import *
    sage: covectors_from_matrix(W, kernel=True, algorithm='fe', separate=True)
    [[(0000)], [(0+0-), (-0+0), (+0-0), (0-0+)], [(-++-), (++--), (+--+), (--++)]]
    sage: covectors_from_matrix(Wt, algorithm='fe', separate=True)
    [[(0000)], [(+00-), (+++0), (-00+), (---0), (0---), (0+++)], [(+++-), (---+), (----), (+---), (++++), (-+++)]]

We can check injectivity by using the function :func:`~cond_inj_intersection`::

    sage: cond_inj_intersection(W, Wt)
    True

Therefore, the corresponding exponential map is injective for all vectors ``c > 0``.
To check surjectivity, we need to verify two conditions.
First, we check the face condition.
For this purpose, we compute the cocircuits of the oriented matroids
corresponding to the matrices::

    sage: cc1 = cocircuits_from_matrix(W, kernel=False)
    sage: cc1
    [(+0+0), (0+0+), (-0-0), (0-0-)]
    sage: cc2 = cocircuits_from_matrix(Wt, kernel=False)
    sage: cc2
    [(+++0), (-00+), (0+++), (---0), (+00-), (0---)]

Here, we are only interested in the positive cocircuits::

    sage: cc1p = [X for X in cc1 if X > 0]
    sage: cc1p
    [(+0+0), (0+0+)]
    sage: cc2p = [X for X in cc2 if X > 0]
    sage: cc2p
    [(+++0), (0+++)]

Since every sign vector in ``cc2p`` has a smaller element in ``cc1p``,
the face condition is satisfied.
There is also a function in the package that can be used directly
to check whether this condition is fulfilled::

    sage: cond_faces(W, Wt)
    True

We need to check a third condition to verify surjectivity.
For this purpose, we consider again the oriented matroid determined by ``W``::

    sage: covectors_from_matrix(W, kernel=True)
    [(0000), (-0+0), (0-0+), (+0-0), (0+0-), (-++-), (++--), (+--+), (--++)]

Since there is no positive covector,
the exponential map is surjective.
The package offers a function to check this condition condition::

    sage: cond_nondegenerate(W, Wt)
    True

Hence, the exponential map is bijective.

Let us consider another example.
We swap the two matrices from before::

    sage: W = matrix([[1,0,0,-1],[0,1,1,1]])
    sage: W
    [ 1  0  0 -1]
    [ 0  1  1  1]
    sage: Wt = matrix([[1,0,1,0],[0,1,0,1]])
    sage: Wt
    [1 0 1 0]
    [0 1 0 1]

Because of symmetry, the corresponding exponential map is injective::

    sage: cond_inj_intersection(W, Wt)
    True

Now, we attempt to check the face condition::

    sage: cc1 = cocircuits_from_matrix(W, kernel=False)
    sage: cc1
    [(+++0), (-00+), (0+++), (---0), (+00-), (0---)]
    sage: cc2 = cocircuits_from_matrix(Wt, kernel=False)
    sage: cc2
    [(+0+0), (0+0+), (-0-0), (0-0-)]

Again, we are only interested in the positive cocircuits::

    sage: cc1p = [X for X in cc1 if X > 0]
    sage: cc1p
    [(+++0), (0+++)]
    sage: cc2p = [X for X in cc2 if X > 0]
    sage: cc2p
    [(+0+0), (0+0+)]

Therefore, condition does not hold.
We also apply the corresponding function from the package::

    sage: cond_faces(W, Wt)
    False

Consequently, this map is not bijective.

Now, we consider Example 20 from [MHR19]_.
Here, we have a parameter ``wt > 0``.
Depending on this parameter, the resulting exponential map will be bijective::

    sage: var('wt')
    wt
    sage: W = matrix(3, 6, [0,0,1,1,-1,0,1,-1,0,0,0,-1,0,0,1,-1,0,0])
    sage: W
    [ 0  0  1  1 -1  0]
    [ 1 -1  0  0  0 -1]
    [ 0  0  1 -1  0  0]
    sage: Wt = matrix(3, 6, [1,1,0,0,-1,wt,1,-1,0,0,0,0,0,0,1,-1,0,0])
    sage: Wt
    [ 1  1  0  0 -1 wt]
    [ 1 -1  0  0  0  0]
    [ 0  0  1 -1  0  0]

The first two conditions depend on the sign vectors of the corresponding oriented matroids.
Consequently, the choice of the positive parameter ``wt`` does not affect the result.
In order to compute the sign vectors, we set ``wt`` to ``1``::

    sage: cond_inj_intersection(W, Wt(wt=1))
    True

Hence, the map is injective.
Also the face condition is satisfied::

    sage: cond_faces(W, Wt(wt=1))
    True

For specific values of ``wt``, the pair of subspaces
determined by kernels of the matrices is non-degenerate.
This is the case for :math:`wt \in (0, 1) \cup (1, 2)`::

    sage: cond_nondegenerate(W, Wt(wt=1/2))
    True
    sage: cond_nondegenerate(W, Wt(wt=3/2))
    True

On the other hand, this condition does not hold if
:math:`wt \in {1} \cup [2, \infty)`.
In this case, the exponential map is injective but not surjective::

    sage: cond_nondegenerate(W, Wt(wt=1))
    False
    sage: cond_nondegenerate(W, Wt(wt=2))
    False
    sage: cond_nondegenerate(W, Wt(wt=3))
    False

We consider some final example::

    sage: from bijectivity_exponential_maps.conditions_bijectivity import nondeg_cond1
    sage: W = matrix([[1,1,0,0],[0,0,1,0]]).right_kernel_matrix()
    sage: Wt = matrix([[1,0,2,0],[0,1,0,-1]])

Now, we check whether the non-degeneracy condition is satisfied::

    sage: nondeg_cond1(W, Wt, certificate=True)
    [False, [[[2], [0, 1]], (2, 2, 4, -2)]]

From the output, we see that the condition is violated.
The vector ``(2, 2, 4, −2)`` lies in the row space
of ``Wt`` and corresponds to the positive sign vectors
``(00+0)`` and ``(++00)``
that are represented by the sets ``{2}`` and ``{0, 1}``.
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

from .utility import normalize, pos_cocircuits_from_matrix, pos_covectors_from_matrix
from elementary_vectors import elementary_vectors, exists_vector
from sign_vectors.oriented_matroids import cocircuits_from_matrix

from sage.modules.free_module_element import zero_vector
from sage.matrix.constructor import matrix
from sage.rings.infinity import Infinity
from sage.misc.flatten import flatten
from sage.geometry.polyhedron.constructor import Polyhedron


def cond_faces(W, Wt):
    r"""
    Check a condition on positive sign vectors of two matrices.

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    OUTPUT:
    Returns whether every positive sign vector ``X`` corresponding to the rows of
    ``Wt`` has a positive sign vector ``Y`` corresponding to the rows of ``W``
    such that ``Y <= X``.

    Returns a boolean.

    EXAMPLES::

        sage: from bijectivity_exponential_maps.conditions_bijectivity import cond_faces
        sage: W = matrix([[1,0,-1,0],[0,1,0,-1]]).right_kernel_matrix()
        sage: W
        [1 0 1 0]
        [0 1 0 1]
        sage: Wt = matrix([[1,0,-1,1],[0,1,-1,0]]).right_kernel_matrix()
        sage: Wt
        [ 1  0  0 -1]
        [ 0  1  1  1]
        sage: cond_faces(W, Wt)
        True
    """
    ccWp = pos_cocircuits_from_matrix(W,  kernel=False)
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
    Check whether the pair of the given matrices is non-degenerate.

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    - ``certificate`` -- a boolean (default: ``False``)

    OUTPUT:
    Let ``S``, ``St`` be the sub spaces corresponding to the matrices ``W``,
    ``Wt``. Returns whether the pair ``(S, St)`` is non-degenerate.

    Returns a boolean.

    - If ``certificate`` is true:

      - If the result is true, returns a vector as a certificate.

    .. SEEALSO::

        :func:`~nondeg_cond1`
        :func:`~nondeg_cond2`
    """
    return nondegenerate(W, Wt, certificate=certificate)


def nondegenerate(W, Wt, certificate=False):
    r"""
    Check whether the pair of the given matrices is non-degenerate.

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    - ``certificate`` -- a boolean (default: ``False``)

    OUTPUT:

    Let ``S``, ``St`` be the sub spaces corresponding to the matrices ``W``,
    ``Wt``. Returns whether the pair ``(S, St)`` is non-degenerate.

    Returns a boolean.

    - If ``certificate`` is true:

      - If the result is true, returns a vector as a certificate.

    .. SEEALSO::

        :func:`~cond_nondegenerate`
        :func:`~nondeg_cond1`
        :func:`~nondeg_cond2`
    """
    c2 = nondeg_cond2(W, Wt)
    if c2:
        return True

    c1 = nondeg_cond1(W, Wt, certificate=certificate)
    return c1


def nondeg_cond1(W, Wt, certificate=False):
    r"""
    Return whether the first condition of ``cond_nondegenerate`` is satisfied.

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    - ``certificate`` -- a boolean (default: ``False``)

    OUTPUT:
    a boolean

    - If ``certificate`` is true:

      - If the result is true, returns a vector as a certificate.

    .. SEEALSO::

        :func:`~cond_nondegenerate`
        :func:`~nondeg_cond2`

    EXAMPLES::

        sage: from bijectivity_exponential_maps.conditions_bijectivity import nondeg_cond1
        sage: W = matrix([[-4,2,-7,1],[-9,-1,-1,-1],[-1,0,-1,1]]).right_kernel_matrix()
        sage: W
        [ 10 -54 -23 -13]
        sage: Wt = matrix([[-5,-1,2,2],[1,0,2,21],[-2,0,0,2]]).right_kernel_matrix()
        sage: Wt
        [  1 -25 -11   1]
        sage: nondeg_cond1(W, Wt)
        False

        sage: W = matrix([[-4,2,-7,1],[-9,-1,-1,-1],[-1,0,-1,1]]).right_kernel_matrix()
        sage: W
        [ 10 -54 -23 -13]
        sage: Wt = matrix([[-5,1,-2,2],[1,0,-2,21],[-2,0,0,2]]).right_kernel_matrix()
        sage: Wt
        [ 1 25 11  1]
        sage: nondeg_cond1(W, Wt)
        True

        sage: W = matrix(2,4,[1,2,0,0,0,0,5,1]).right_kernel_matrix()
        sage: W
        [ 2 -1  0  0]
        [ 0  0  1 -5]
        sage: A = matrix([[1,1,2,2],[1,0,1,-1]]).right_kernel_matrix()
        sage: A
        [ 1  1 -1  0]
        [ 0  4 -1 -1]
        sage: Wt = A.right_kernel_matrix()
        sage: Wt
        [ 1  0  1 -1]
        [ 0  1  1  3]
        sage: nondeg_cond1(W, Wt)
        False

        sage: W = matrix(3,5,[1,2,0,0,0,0,0,5,1,0,0,0,0,0,1]).right_kernel_matrix()
        sage: W
        [ 2 -1  0  0  0]
        [ 0  0  1 -5  0]
        sage: A = matrix([[1,1,-2,-2,2],[1,0,1,-1,2]]).right_kernel_matrix()
        sage: A
        [ 1  1  0  1  0]
        [ 0  2  0  2  1]
        [ 0  0  1 -3 -2]
        sage: Wt = A.right_kernel_matrix()
        sage: Wt
        [ 1  0  1 -1  2]
        [ 0  1 -3 -1  0]
        sage: nondeg_cond1(W, Wt)
        False

        sage: A = matrix([[1, 0, 0, 1], [0, 1, 1, -1]]).right_kernel_matrix()
        sage: A
        [ 1  0 -1 -1]
        [ 0  1 -1  0]
        sage: B = matrix([[1,1,0,0],[0,0,1,1]])
        sage: B
        [1 1 0 0]
        [0 0 1 1]
        sage: nondeg_cond1(B, A)
        True
    """
    if W.ncols() != Wt.ncols():
        raise ValueError('Matrices have different number of columns.')
    # eventuell auf disjoint support achten. Sind +0 und 0+ dabei, dann müssen wir ++ nicht betrachten.
    P = pos_covectors_from_matrix(W, kernel=True)

    if P == []:
        return True

    n = Wt.ncols()  # length of vectors
    degenerate = False  # might change in recursion

    L = [-Infinity for i in range(n)]
    R = [0 for i in range(n)]
    inf = [Infinity for i in range(n)]  # does not change

#    l = [True for i in range(n)]
    proof = []

    M = Wt.right_kernel_matrix()

    def rec(P, M, I, L, R):
        r"""
        Recursive function.

        INPUT:

        - ``P`` -- a list of positive sign vectors

        - ``M`` -- a matrix

        - ``I`` -- a list of indices

        - ``L`` -- a list of values ``-Infinity`` and ``1``

        - ``R`` -- a list of values ``0`` and ``Infinity``
        """
        nonlocal degenerate  # use nonlocal to use the same ``degenerate``
        nonlocal proof

        while P:
            X = P.pop()
            if set(flatten(I)).issubset(X.zero_support()):  # P[0].positive_support() subseteq J, insbesondere muss X ein + auf J haben
                L_ = L[:]
                R_ = R[:]
                for i in X.support():
                    L_[i] = 1
                    R_[i] = Infinity

                M_ = equal_components(M, X.support())  # make to echelon form?
                A = M_.right_kernel_matrix()

                if A:  # A is not empty matrix
                    evs = elementary_vectors(A, kernel=True)

                    if exists_vector(evs, L_, R_):
                        degenerate = True

                        I += [X.support()]
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
    Return a vector.

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
            v[0] = -1  # component - 1 >= 0 # instead of > 0
            v[i+1] = 1
            ieqs.append(v)
        else:
            v[i+1] = -1  # component <= 0
            ieqs.append(v)

#     print('ieqs:', ieqs)
#     print('eqns:', eqns)
    P = Polyhedron(ieqs=ieqs, eqns=eqns)
    return P.representative_point()


def nondeg_cond2(W, Wt):
    r"""
    Check a condition an matrices.

    .. SEEALSO::

        :func:`~cond_nondegenerate`
        :func:`~nondeg_cond1`
    """
    ccWt = normalize(cocircuits_from_matrix(Wt))
    ccWp = pos_cocircuits_from_matrix(W)
    for X in ccWt:
        val = False
        for Y in ccWp:
            if set(X.zero_support()).issubset(Y.zero_support()):
                val = True
                break
        if not val:
            return False
    return True
