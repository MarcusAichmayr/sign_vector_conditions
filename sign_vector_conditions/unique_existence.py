r"""
Existence and uniqueness of equilibria

EXAMPLES::

    sage: from sign_vector_conditions import *

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
    [{(0000)}, {(0-0+), (0+0-), (+0-0), (-0+0)}, {(--++), (-++-), (++--), (+--+)}]
    sage: covectors_from_matrix(Wt, kernel=False, algorithm='fe', separate=True)
    [{(0000)},
     {(+++0), (0+++), (+00-), (-00+), (0---), (---0)},
     {(-+++), (----), (---+), (+---), (++++), (+++-)}]


We can check injectivity by using the function :func:`~condition_uniqueness_signvectors`::

    sage: condition_uniqueness_signvectors(W, Wt)
    True

Therefore, the corresponding exponential map is injective for all vectors ``c > 0``.
To check surjectivity, we need to verify two conditions.
First, we check the face condition.
For this purpose, we compute the cocircuits of the oriented matroids
corresponding to the matrices::

    sage: cc1 = cocircuits_from_matrix(W, kernel=False)
    sage: cc1
    {(+0+0), (0-0-), (-0-0), (0+0+)}
    sage: cc2 = cocircuits_from_matrix(Wt, kernel=False)
    sage: cc2
    {(+++0), (0+++), (+00-), (-00+), (0---), (---0)}

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

    sage: condition_faces(W, Wt)
    True

We need to check a third condition to verify surjectivity.
For this purpose, we consider again the oriented matroid determined by ``W``::

    sage: covectors_from_matrix(W, kernel=True)
    {(0000), (+--+), (++--), (-0+0), (0-0+), (--++), (0+0-), (-++-), (+0-0)}

Since there is no positive covector,
the exponential map is surjective.
The package offers a function to check this condition condition::

    sage: condition_nondegenerate(W, Wt)
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

    sage: condition_uniqueness_signvectors(W, Wt)
    True

Now, we attempt to check the face condition::

    sage: cc1 = cocircuits_from_matrix(W, kernel=False)
    sage: cc1
    {(+++0), (0+++), (+00-), (-00+), (0---), (---0)}
    sage: cc2 = cocircuits_from_matrix(Wt, kernel=False)
    sage: cc2
    {(+0+0), (0-0-), (-0-0), (0+0+)}

Again, we are only interested in the positive cocircuits::

    sage: cc1p = [X for X in cc1 if X > 0]
    sage: cc1p
    [(+++0), (0+++)]
    sage: cc2p = [X for X in cc2 if X > 0]
    sage: cc2p
    [(+0+0), (0+0+)]

Therefore, condition does not hold.
We also apply the corresponding function from the package::

    sage: condition_faces(W, Wt)
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

    sage: condition_uniqueness_signvectors(W, Wt(wt=1))
    True

Hence, the map is injective.
Also the face condition is satisfied::

    sage: condition_faces(W, Wt(wt=1))
    True

For specific values of ``wt``, the pair of subspaces
determined by kernels of the matrices is non-degenerate.
This is the case for :math:`wt \in (0, 1) \cup (1, 2)`::

    sage: condition_nondegenerate(W, Wt(wt=1/2))
    True
    sage: condition_nondegenerate(W, Wt(wt=3/2))
    True

On the other hand, this condition does not hold if
:math:`wt \in {1} \cup [2, \infty)`.
In this case, the exponential map is injective but not surjective::

    sage: condition_nondegenerate(W, Wt(wt=1))
    False
    sage: condition_nondegenerate(W, Wt(wt=2))
    False
    sage: condition_nondegenerate(W, Wt(wt=3))
    False

We consider some final example::

    sage: from sign_vector_conditions.unique_existence import nondeg_cond1
    sage: W = matrix([[1,1,0,0],[0,0,1,0]]).right_kernel_matrix()
    sage: Wt = matrix([[1,0,2,0],[0,1,0,-1]])

Now, we check whether the non-degeneracy condition is satisfied::

    sage: nondeg_cond1(W, Wt, certificate=True)
    [False, [[[2], [0, 1]], (1, 1, 2, -1)]]

From the output, we see that the condition is violated.
The vector ``(1, 1, 2, -1)`` lies in the row space
of ``Wt`` and corresponds to the positive sign vectors
``(00+0)`` and ``(++00)``
that are represented by the sets ``{2}`` and ``{0, 1}``.
"""

#############################################################################
#  Copyright (C) 2023                                                       #
#                Marcus Aichmayr (aichmayr@mathematik.uni-kassel.de)        #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from .utility import normalize, pos_cocircuits_from_matrix, pos_covectors_from_matrix
from elementary_vectors import elementary_vectors
from elementary_vectors import setup_intervals, exists_vector, exists_orthogonal_vector, construct_vector
from sign_vectors.oriented_matroids import cocircuits_from_matrix

from sage.modules.free_module_element import zero_vector
from sage.matrix.constructor import matrix
from sage.rings.infinity import Infinity
from sage.misc.flatten import flatten


def condition_faces(W, Wt):
    r"""
    Condition on positive sign vectors for existence and uniqueness of equilibria

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    OUTPUT:
    Returns whether every positive sign vector ``X`` corresponding to the rows of
    ``Wt`` has a positive sign vector ``Y`` corresponding to the rows of ``W``
    such that ``Y <= X``.

    Returns a boolean.

    EXAMPLES::

        sage: from sign_vector_conditions.unique_existence import condition_faces
        sage: W = matrix([[1,0,-1,0],[0,1,0,-1]]).right_kernel_matrix()
        sage: W
        [1 0 1 0]
        [0 1 0 1]
        sage: Wt = matrix([[1,0,-1,1],[0,1,-1,0]]).right_kernel_matrix()
        sage: Wt
        [ 1  0  0 -1]
        [ 0  1  1  1]
        sage: condition_faces(W, Wt)
        True
    """
    positive_cocircuits = pos_cocircuits_from_matrix(W,  kernel=False)

    for cocircuit1 in pos_cocircuits_from_matrix(Wt, kernel=False):
        value = True
        for cocircuit2 in positive_cocircuits:
            if cocircuit2 <= cocircuit1:
                value = False
                break
        if value:
            return False
    return True


def condition_nondegenerate(W, Wt, certificate=False):
    r"""
    Non-degeneracy condition for existence and uniqueness of equilibria

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    - ``certificate`` -- a boolean (default: ``False``)

    OUTPUT:
    Let ``S``, ``St`` be the subspaces corresponding to the matrices ``W``,
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
    Check whether the pair of given matrices is non-degenerate.

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    - ``certificate`` -- a boolean (default: ``False``)

    OUTPUT:
    Let ``S``, ``St`` be the subspaces corresponding to the matrices ``W``,
    ``Wt``. Returns whether the pair ``(S, St)`` is non-degenerate.

    Returns a boolean.

    - If ``certificate`` is true:

      - If the result is true, returns a vector as a certificate.

    .. SEEALSO::

        :func:`~condition_nondegenerate`
        :func:`~nondeg_cond1`
        :func:`~nondeg_cond2`
    """
    if nondeg_cond2(W, Wt):
        return True

    return nondeg_cond1(W, Wt, certificate=certificate)


def nondeg_cond1(W, Wt, certificate=False):
    r"""
    Return whether the first condition of ``condition_nondegenerate`` is satisfied.

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    - ``certificate`` -- a boolean (default: ``False``)

    OUTPUT:
    a boolean

    - If ``certificate`` is true:

      - If the result is true, returns a vector as a certificate.

    .. SEEALSO::

        :func:`~condition_nondegenerate`
        :func:`~nondeg_cond2`

    EXAMPLES::

        sage: from sign_vector_conditions.unique_existence import nondeg_cond1
        sage: W = matrix([[-4, 2, -7, 1], [-9, -1, -1, -1], [-1, 0, -1, 1]]).right_kernel_matrix()
        sage: W
        [ 10 -54 -23 -13]
        sage: Wt = matrix([[-5, -1, 2, 2], [1, 0, 2, 21], [-2, 0, 0, 2]]).right_kernel_matrix()
        sage: Wt
        [  1 -25 -11   1]
        sage: nondeg_cond1(W, Wt)
        False
        sage: W = matrix([[-4, 2, -7, 1], [-9, -1, -1, -1], [-1, 0, -1, 1]]).right_kernel_matrix()
        sage: W
        [ 10 -54 -23 -13]
        sage: Wt = matrix([[-5, 1, -2, 2], [1, 0, -2, 21], [-2, 0, 0, 2]]).right_kernel_matrix()
        sage: Wt
        [ 1 25 11  1]
        sage: nondeg_cond1(W, Wt)
        True
        sage: W = matrix(2, 4, [1, 2, 0, 0, 0, 0, 5, 1]).right_kernel_matrix()
        sage: W
        [ 2 -1  0  0]
        [ 0  0  1 -5]
        sage: A = matrix([[1, 1, 2, 2], [1, 0, 1, -1]]).right_kernel_matrix()
        sage: A
        [ 1  1 -1  0]
        [ 0  4 -1 -1]
        sage: Wt = A.right_kernel_matrix()
        sage: Wt
        [ 1  0  1 -1]
        [ 0  1  1  3]
        sage: nondeg_cond1(W, Wt)
        False
        sage: W = matrix(3, 5, [1, 2, 0, 0, 0, 0, 0, 5, 1, 0, 0, 0, 0, 0, 1]).right_kernel_matrix()
        sage: W
        [ 2 -1  0  0  0]
        [ 0  0  1 -5  0]
        sage: A = matrix([[1, 1, -2, -2, 2], [1, 0, 1, -1, 2]]).right_kernel_matrix()
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
        sage: B = matrix([[1, 1, 0, 0], [0, 0, 1, 1]])
        sage: B
        [1 1 0 0]
        [0 0 1 1]
        sage: nondeg_cond1(B, A)
        True
    """
    if W.ncols() != Wt.ncols():
        raise ValueError('Matrices have different number of columns.')
    # TODO: consider disjoint support: If we have "+0" and "0+", then we do not need to consider "++".
    positive_covectors = pos_covectors_from_matrix(W, kernel=True)

    if not positive_covectors:
        return True

    length = Wt.ncols()
    degenerate = False

    lower_bounds = [-Infinity for i in range(length)]
    upper_bounds = [0 for i in range(length)]
    inf = [Infinity for i in range(length)]

    proof = []

    kernel_matrix = Wt.right_kernel_matrix()

    def rec(positive_covectors, kernel_matrix, indices, lower_bounds, upper_bounds):
        r"""
        Recursive function.

        INPUT:

        - ``positive_covectors`` -- a list of positive sign vectors

        - ``kernel_matrix`` -- a matrix

        - ``indices`` -- a list of indices

        - ``lower_bounds`` -- a list of values ``-Infinity`` and ``1``

        - ``upper_bounds`` -- a list of values ``0`` and ``Infinity``
        """
        nonlocal degenerate
        nonlocal proof

        while positive_covectors:
            covector = positive_covectors.pop()
            if set(flatten(indices)).issubset(covector.zero_support()):
                # covector must have a "+" on zero support
                lower_bounds_new = lower_bounds[:]
                upper_bounds_new = upper_bounds[:]
                for i in covector.support():
                    lower_bounds_new[i] = 1
                    upper_bounds_new[i] = Infinity

                matrix_equal_components = matrix_with_equal_components(kernel_matrix, covector.support())
                evs = elementary_vectors(matrix_equal_components.right_kernel_matrix())
                intervals = setup_intervals(lower_bounds_new, upper_bounds_new)

                if exists_vector(evs, intervals):
                    degenerate = True
                    indices += [covector.support()]
                    if certificate:
                        proof = [indices, construct_vector(Wt, intervals)]
                    return
                intervals = setup_intervals(lower_bounds_new, inf)
                if exists_vector(evs, intervals):
                    rec(positive_covectors[:], matrix_equal_components, indices + [covector.support()], lower_bounds_new, upper_bounds_new)
                else:
                    if certificate:
                        for element in evs:
                            if exists_orthogonal_vector(element, intervals):
                                proof.append([element, indices + [covector.support()]])
                                break

            if degenerate:
                return
        return

    rec(positive_covectors, kernel_matrix, [], lower_bounds, upper_bounds)

    if certificate:
        return [not degenerate, proof]
    return not degenerate


def matrix_with_equal_components(M, indices):
    r"""
    Insert additional rows to a matrix, such that given components of the kernel matrix are equal.

    INPUT:

    - ``M`` -- a matrix

    - ``indices`` -- a list of indices

    OUTPUT:
    a matrix with additional rows.

    EXAMPLES::

        sage: M = matrix([[1, 0, 1, 0], [0, 0 ,-1, 2]])
        sage: M
        [ 1  0  1  0]
        [ 0  0 -1  2]
        sage: from sign_vector_conditions.unique_existence import matrix_with_equal_components
        sage: matrix_with_equal_components(M, [0, 1])
        [ 1  0  1  0]
        [ 0  0 -1  2]
        [ 1 -1  0  0]
        sage: _.right_kernel_matrix()
        [ 2  2 -2 -1]
        sage: matrix_with_equal_components(M, [0, 1, 2])
        [ 1  0  1  0]
        [ 0  0 -1  2]
        [ 1 -1  0  0]
        [ 1  0 -1  0]
        sage: _.right_kernel_matrix()
        []
    """
    length = M.ncols()

    if len(indices) >= 2:
        i = indices[0]
        element = zero_vector(length)
        element[i] = 1
        for j in indices[1:]:
            row = element[:]
            row[j] = -1
            M = matrix(list(M) + [row])
    return M


def nondeg_cond2(W, Wt):
    r"""
    Check a condition an matrices.

    .. SEEALSO::

        :func:`~condition_nondegenerate`
        :func:`~nondeg_cond1`

    EXAMPLES::

        sage: from sign_vector_conditions.unique_existence import nondeg_cond2
        sage: W = matrix([[-4,2,-7,1],[-9,-1,-1,-1],[-1,0,-1,1]]).right_kernel_matrix()
        sage: W
        [ 10 -54 -23 -13]
        sage: Wt = matrix([[-5,-1,2,2],[1,0,2,21],[-2,0,0,2]]).right_kernel_matrix()
        sage: Wt
        [  1 -25 -11   1]
        sage: nondeg_cond2(W, Wt)
        False
        sage: W = matrix([[-4,2,-7,1],[-9,-1,-1,-1],[-1,0,-1,1]]).right_kernel_matrix()
        sage: W
        [ 10 -54 -23 -13]
        sage: Wt = matrix([[-5,1,-2,2],[1,0,-2,21],[-2,0,0,2]]).right_kernel_matrix()
        sage: Wt
        [ 1 25 11  1]
        sage: nondeg_cond2(W, Wt)
        False
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
        sage: nondeg_cond2(W, Wt)
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
        sage: nondeg_cond2(W, Wt)
        False
        sage: A = matrix([[1, 0, 0, 1], [0, 1, 1, -1]]).right_kernel_matrix()
        sage: A
        [ 1  0 -1 -1]
        [ 0  1 -1  0]
        sage: B = matrix([[1,1,0,0],[0,0,1,1]])
        sage: B
        [1 1 0 0]
        [0 0 1 1]
        sage: nondeg_cond2(B, A)
        False
    """
    positive_cocircuits = pos_cocircuits_from_matrix(W, kernel=False)
    for cocircuit1 in normalize(cocircuits_from_matrix(Wt, kernel=False)):
        value = False
        for cocircuit2 in positive_cocircuits:
            if set(cocircuit1.zero_support()).issubset(cocircuit2.zero_support()):
                value = True
                break
        if not value:
            return False
    return True
