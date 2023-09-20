r"""
Existence and uniqueness of equilibria

EXAMPLES::

    sage: from sign_vector_conditions import *

Let us consider the following matrices::

    sage: W = matrix([[1, 0, 1, 0], [0, 1, 0, 1]])
    sage: W
    [1 0 1 0]
    [0 1 0 1]
    sage: Wt = matrix([[1, 0, 0, -1], [0, 1, 1, 1]])
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

    sage: condition_subspaces_nondegenerate(W, Wt)
    True

Hence, the exponential map is bijective.

Let us consider another example.
We swap the two matrices from before::

    sage: W = matrix([[1, 0, 0, -1], [0, 1, 1, 1]])
    sage: W
    [ 1  0  0 -1]
    [ 0  1  1  1]
    sage: Wt = matrix([[1, 0, 1, 0], [0, 1, 0, 1]])
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
    sage: W = matrix(3, 6, [0, 0, 1, 1, -1, 0, 1, -1, 0, 0, 0, -1, 0, 0, 1, -1, 0, 0])
    sage: W
    [ 0  0  1  1 -1  0]
    [ 1 -1  0  0  0 -1]
    [ 0  0  1 -1  0  0]
    sage: Wt = matrix(3, 6, [1, 1, 0, 0, -1, wt, 1, -1, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0])
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

    sage: condition_subspaces_nondegenerate(W, Wt(wt=1/2))
    True
    sage: condition_subspaces_nondegenerate(W, Wt(wt=3/2))
    True

On the other hand, this condition does not hold if
:math:`wt \in {1} \cup [2, \infty)`::

    sage: condition_subspaces_nondegenerate(W, Wt(wt=1))
    False

To certify the result, we call::

    sage: condition_subspaces_degenerate(W, Wt(wt=1), certify=True)
    [True, [[[0, 1, 5]], (1, 1, 0, 0, -1, 1)]]

Hence, the positive support of the vector ``v = (1, 1, 0, 0, -1, 1)`` of ``Wt``
can be covered by a sign vector ``(++000+)`` corresponding to ``ker(W)``.
Further, ``v`` does not satisfy the support condition.

    sage: condition_subspaces_nondegenerate(W, Wt(wt=2))
    False
    sage: condition_subspaces_nondegenerate(W, Wt(wt=3))
    False
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

from copy import copy

from sage.matrix.constructor import matrix
from sage.misc.flatten import flatten
from sage.rings.infinity import Infinity

from elementary_vectors import elementary_vectors
from elementary_vectors import setup_intervals, exists_vector, exists_orthogonal_vector, construct_vector

from .utility import positive_cocircuits_from_matrix, positive_covectors_from_matrix, equal_entries_lists, degenerate_condition_support


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
        sage: W = matrix([[1, 0, -1, 0], [0, 1, 0, -1]]).right_kernel_matrix()
        sage: W
        [1 0 1 0]
        [0 1 0 1]
        sage: Wt = matrix([[1, 0, -1, 1], [0, 1, -1, 0]]).right_kernel_matrix()
        sage: Wt
        [ 1  0  0 -1]
        [ 0  1  1  1]
        sage: condition_faces(W, Wt)
        True
    """
    positive_cocircuits = positive_cocircuits_from_matrix(W, kernel=False)

    for cocircuit1 in positive_cocircuits_from_matrix(Wt, kernel=False):
        value = True
        for cocircuit2 in positive_cocircuits:
            if cocircuit2 <= cocircuit1:
                value = False
                break
        if value:
            return False
    return True


def condition_subspaces_nondegenerate(W, Wt):
    r"""
    Return whether the subspaces given by to matrices are non-degenerate

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    OUTPUT:
    a boolean
    
    .. SEEALSO::
    
        :func:`~condition_subspaces_degenerate`
    """
    return not condition_subspaces_degenerate(W, Wt)


def condition_subspaces_degenerate(W, Wt, certify=False):
    r"""
    Return whether the subspaces given by to matrices are degenerate

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    - ``certify`` -- a boolean (default: ``False``)

    OUTPUT:
    a boolean

    If ``certify`` is true, the result will be certified.

    EXAMPLES::

        sage: from sign_vector_conditions.unique_existence import *
        sage: W = matrix([[-4, 2, -7, 1], [-9, -1, -1, -1], [-1, 0, -1, 1]]).right_kernel_matrix()
        sage: W
        [ 10 -54 -23 -13]
        sage: Wt = matrix([[-5, -1, 2, 2], [1, 0, 2, 21], [-2, 0, 0, 2]]).right_kernel_matrix()
        sage: Wt
        [  1 -25 -11   1]
        sage: condition_subspaces_degenerate(W, Wt)
        True
        sage: W = matrix([[-4, 2, -7, 1], [-9, -1, -1, -1], [-1, 0, -1, 1]]).right_kernel_matrix()
        sage: W
        [ 10 -54 -23 -13]
        sage: Wt = matrix([[-5, 1, -2, 2], [1, 0, -2, 21], [-2, 0, 0, 2]]).right_kernel_matrix()
        sage: Wt
        [ 1 25 11  1]
        sage: condition_subspaces_degenerate(W, Wt)
        False
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
        sage: condition_subspaces_degenerate(W, Wt)
        True
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
        sage: condition_subspaces_degenerate(W, Wt)
        True
        sage: A = matrix([[1, 0, 0, 1], [0, 1, 1, -1]]).right_kernel_matrix()
        sage: A
        [ 1  0 -1 -1]
        [ 0  1 -1  0]
        sage: B = matrix([[1, 1, 0, 0], [0, 0, 1, 1]])
        sage: B
        [1 1 0 0]
        [0 0 1 1]
        sage: condition_subspaces_degenerate(B, A)
        False
    """
    if W.ncols() != Wt.ncols():
        raise ValueError('Matrices have different number of columns.')
    # TODO: consider disjoint support: If we have "+0" and "0+", then we do not need to consider "++".
    positive_covectors = positive_covectors_from_matrix(W, kernel=True)

    if not positive_covectors:
        if certify:
            return [False, "no positive covectors"]
        return False

    length = Wt.ncols()
    degenerate = False

    lower_bounds = [-Infinity] * length
    upper_bounds = [0] * length
    inf = [Infinity] * length

    certificate = []

    kernel_matrix = Wt.right_kernel_matrix()
    
    covectors_support_condition = positive_cocircuits_from_matrix(W, kernel=False)

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
        nonlocal certificate

        while positive_covectors:
            covector = positive_covectors.pop()
            if set(flatten(indices)).issubset(covector.zero_support()):
                lower_bounds_new = copy(lower_bounds)
                upper_bounds_new = copy(upper_bounds)
                for i in covector.support():
                    lower_bounds_new[i] = 1
                    upper_bounds_new[i] = Infinity

                matrix_equal_components = matrix(kernel_matrix.rows() + equal_entries_lists(length, covector.support()))
                evs = elementary_vectors(matrix_equal_components.right_kernel_matrix())
                intervals = setup_intervals(lower_bounds_new, upper_bounds_new)

                if exists_vector(evs, intervals):
                    if degenerate_condition_support(matrix_equal_components, intervals, covectors_support_condition):
                        degenerate = True
                        indices += [covector.support()]
                        if certify:
                            certificate = [indices, construct_vector(Wt, intervals)]
                        return
                intervals = setup_intervals(lower_bounds_new, inf)
                if exists_vector(evs, intervals):
                    rec(copy(positive_covectors),
                            matrix_equal_components,
                            indices + [covector.support()],
                            lower_bounds_new,
                            upper_bounds_new
                    )
                elif certify:
                    for element in evs:
                        if exists_orthogonal_vector(element, intervals):
                            certificate.append([element, indices + [covector.support()]])
                            break

            if degenerate:
                return
        return

    rec(positive_covectors, kernel_matrix, [], lower_bounds, upper_bounds)

    if certify:
        return [degenerate, certificate]
    return degenerate
