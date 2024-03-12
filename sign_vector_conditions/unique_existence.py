r"""
Existence and uniqueness of equilibria
======================================

Let us consider the following matrices to describe a chemical reaction network::

    sage: W = matrix([[1, 0, 1, 0], [0, 1, 0, 1]])
    sage: W
    [1 0 1 0]
    [0 1 0 1]
    sage: Wt = matrix([[1, 0, 0, -1], [0, 1, 1, 1]])
    sage: Wt
    [ 1  0  0 -1]
    [ 0  1  1  1]

To check whether a unique equilibrium exists, we apply :func:`~condition_uniqueness_minors`::

    sage: from sign_vector_conditions import *
    sage: condition_uniqueness_minors(W, Wt)
    True

This means that the chemical reaction network has at most one equilibrium.
Next, we verify whether an equilibrium exists.
First, we check the face condition.
For this purpose, we compute the cocircuits of the oriented matroids
corresponding to the matrices::

    sage: from sign_vectors.oriented_matroids import *
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

Since there are no nonnegative covectors, the chemical reaction network has at least one equilibrium.
The package offers a function to check this condition condition::

    sage: condition_nondegenerate(W, Wt)
    True

Hence, the chemical reaction network has a unique equilibrium.

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

Because of symmetry, there is at most one equilibrium::

    sage: condition_uniqueness_sign_vectors(W, Wt)
    True

Now, we check the face condition::

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

Therefore, the condition does not hold.
We also apply the corresponding function from the package::

    sage: condition_faces(W, Wt)
    False

Consequently, there exists no unique equilibrium.

Now, we consider Example 20 from [MHR19]_.
Here, we have a parameter ``a > 0``.
Depending on this parameter, the chemical reaction network has a unique equilibrium::

    sage: var('a')
    a
    sage: W = matrix(3, 6, [0, 0, 1, 1, -1, 0, 1, -1, 0, 0, 0, -1, 0, 0, 1, -1, 0, 0])
    sage: W
    [ 0  0  1  1 -1  0]
    [ 1 -1  0  0  0 -1]
    [ 0  0  1 -1  0  0]
    sage: Wt = matrix(3, 6, [1, 1, 0, 0, -1, a, 1, -1, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0])
    sage: Wt
    [ 1  1  0  0 -1  a]
    [ 1 -1  0  0  0  0]
    [ 0  0  1 -1  0  0]

The first two conditions depend on the sign vectors of the corresponding oriented matroids.
Consequently, the choice of the positive parameter ``a`` does not affect the result::

    sage: assume(a > 0)
    sage: condition_uniqueness_sign_vectors(W, Wt)
    True

Hence, there exists at most one equilibrium.
Also the face condition is satisfied::

    sage: condition_faces(W, Wt)
    True

For specific values of ``a``, the pair of subspaces
determined by kernels of the matrices is nondegenerate.
This is the case for :math:`a \in (0, 1) \cup (1, 2)`::

    sage: condition_nondegenerate(W, Wt(a=1/2))
    True
    sage: condition_nondegenerate(W, Wt(a=3/2))
    True

On the other hand, this condition does not hold if
:math:`a \in {1} \cup [2, \infty)`::

    sage: condition_nondegenerate(W, Wt(a=1))
    False

To certify the result, we call::

    sage: condition_degenerate(W, Wt(a=1), certify=True)
    (True, (1, 1, 0, 0, -1, 1))

Hence, the positive support of the vector ``v = (1, 1, 0, 0, -1, 1)`` of ``Wt``
can be covered by a sign vector ``(++000+)`` corresponding to ``ker(W)``.
Further, ``v`` does not satisfy the support condition.

    sage: condition_nondegenerate(W, Wt(a=2))
    False
    sage: condition_nondegenerate(W, Wt(a=3))
    False
"""

#############################################################################
#  Copyright (C) 2024                                                       #
#                Marcus Aichmayr (aichmayr@mathematik.uni-kassel.de)        #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from copy import copy

from sage.matrix.constructor import matrix
from sage.rings.infinity import Infinity

from elementary_vectors import elementary_vectors
from vectors_in_intervals import intervals_from_bounds, exists_vector, vector_from_sign_vector, sign_vectors_in_intervals, intervals_from_sign_vector

from .utility import non_negative_cocircuits_from_matrix, equal_entries_lists


def condition_faces(W, Wt):
    r"""
    Condition on positive sign vectors for existence and uniqueness of equilibria

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    OUTPUT:
    Return whether every positive sign vector ``X`` corresponding to the rows of
    ``Wt`` has a positive sign vector ``Y`` corresponding to the rows of ``W``
    such that ``Y <= X``.

    Return a boolean.

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
    non_negative_cocircuits = non_negative_cocircuits_from_matrix(W, kernel=False)

    for cocircuit1 in non_negative_cocircuits_from_matrix(Wt, kernel=False):
        if not any(cocircuit2 <= cocircuit1 for cocircuit2 in non_negative_cocircuits):
            return False
    return True


def condition_nondegenerate(W, Wt):
    r"""
    Return whether a pair of subspaces given by matrices is nondegenerate.

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    OUTPUT:
    a boolean
    
    .. SEEALSO::
    
        :func:`~condition_degenerate`
    """
    return not condition_degenerate(W, Wt)


def condition_degenerate(W, Wt, certify=False):
    r"""
    Return whether a pair of subspaces given by matrices is degenerate.

    This condition is about whether all positive equal components of a vector in ``Wt``
    can be covered by covectors corresponding to the kernel of ``W``.

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    - ``certify`` -- a boolean (default: ``False``)

    OUTPUT:
    a boolean

    If ``certify`` is true, a list is returned to certify the result.
    (see the examples)

    EXAMPLES::

        sage: from sign_vector_conditions.unique_existence import *

    Next, we certify our results. In the first examples, the subspaces are trivially nondegenerate
    since there are no nonnegative covectors in the kernel of ``W``::

        sage: W = matrix([[1, 1, 0, 0], [0, 0, 1, 1]])
        sage: Wt = matrix([[1, 1, 0, -1], [0, 0, 1, 0]])
        sage: condition_degenerate(W, Wt, certify=True)
        (False, 'no nonnegative covectors')

    Here, we have a pair of degenerate subspaces::

        sage: W = matrix([[1, -1, 0], [0, 0, 1]])
        sage: Wt = matrix([[1, 0, 0], [0, 1, 0]])
        sage: condition_degenerate(W, Wt, certify=True)
        (True, (1, 1, 0))
    
    The resulting vector lies in the row space of ``Wt``.
    The nonnegative covector ``(++0)`` in the kernel of ``W`` covers the first two equal components.

    In the following, we have another example for nondegenerate subspaces::

        sage: W = matrix([[1, 1, 0, -1, 0], [0, 0, 1, -1, -1]])
        sage: Wt = matrix([[1, 1, 0, -1, 0], [0, 0, 1, 1, 1]])
        sage: condition_degenerate(W, Wt, certify=True)
        (False, ([[[0, 2, 3]], [[1, 2, 3]]], [[[2, 4]]], []))

    The certificate tells us that there is no vector in the row space of ``Wt``
    with positive support on the components ``0, 2, 3`` and ``1, 2, 3``.
    Positive equal components can partially be covered by a covector ``(00+0+)``
    which corresponds to ``[[2, 4]]``.
    However, it is impossible to fully cover the positive support.
    
    In the next example, there exists a partial cover::

        sage: W = matrix([[1, 1, 0, 0], [0, 0, 1, -1]])
        sage: Wt = matrix([[1, 1, 0, -1], [0, 0, 1, 0]])
        sage: condition_degenerate(W, Wt, certify=True)
        (False, ([], [[[2, 3]]], [[[[2, 3]], [(--++)]]]))

    In fact, a vector in ``Wt`` with equal positive components on ``[2, 3]``
    corresponding to ``(--++)`` can be fully covered by covectors.
    However, this vector would not satisfy the support condition.
    """
    if W.ncols() != Wt.ncols():
        raise ValueError('Matrices have different number of columns.')
    non_negative_cocircuits = non_negative_cocircuits_from_matrix(W, kernel=True)

    if not non_negative_cocircuits:
        if certify:
            return False, "no nonnegative covectors"
        return False

    non_negative_cocircuits = sorted(non_negative_cocircuits, key=lambda covector: len(covector.support()))
    length = Wt.ncols()
    is_degenerate = False

    lower_bounds = [-Infinity] * length
    upper_bounds = [0] * length
    upper_bounds_inf = [Infinity] * length

    kernel_matrix = Wt.right_kernel_matrix()
    covectors_support_condition = non_negative_cocircuits_from_matrix(W, kernel=False)

    if certify:
        certificate = []
        certificates_zero_equal_components = []
        certificates_partial_cover = []
        certificate_support_condition = []

    def recursive_degenerate(non_negative_cocircuits, kernel_matrix, indices, lower_bounds, upper_bounds):
        r"""
        Recursive function.

        INPUT:

        - ``non_negative_cocircuits`` -- a list of positive sign vectors

        - ``kernel_matrix`` -- a matrix

        - ``indices`` -- a list of indices

        - ``lower_bounds`` -- a list of values ``-Infinity`` and ``1``

        - ``upper_bounds`` -- a list of values ``0`` and ``Infinity``
        """
        nonlocal is_degenerate
        nonlocal certificate

        while non_negative_cocircuits:
            covector = non_negative_cocircuits.pop()
            lower_bounds_new = copy(lower_bounds)
            upper_bounds_new = copy(upper_bounds)
            for i in covector.support():
                lower_bounds_new[i] = 1
                upper_bounds_new[i] = Infinity

            intervals = intervals_from_bounds(lower_bounds_new, upper_bounds_new)
            indices_new = indices + [covector.support()]
            kernel_matrix_new = matrix(kernel_matrix.rows() + equal_entries_lists(length, covector.support())).echelon_form()
            evs_kernel = elementary_vectors(kernel_matrix_new.right_kernel_matrix())

            if exists_vector(evs_kernel, intervals):
                if certify:
                    covectors_certificate_support_condition = []
                for sign_pattern in sign_vectors_in_intervals(intervals):
                    if not exists_vector(evs_kernel, intervals_from_sign_vector(sign_pattern)):
                        continue
                    if not any(set(cocircuit.support()).issubset(sign_pattern.support()) for cocircuit in covectors_support_condition):
                        is_degenerate = True
                        if certify:
                            certificate = vector_from_sign_vector(sign_pattern, elementary_vectors(kernel_matrix_new))
                        return
                    if certify:
                        covectors_certificate_support_condition.append(sign_pattern)
                if certify:
                    certificate_support_condition.append([indices_new, covectors_certificate_support_condition])

            if exists_vector(evs_kernel, intervals_from_bounds(lower_bounds_new, upper_bounds_inf)):
                if certify:
                    certificates_partial_cover.append(indices_new)
                recursive_degenerate(copy(non_negative_cocircuits), kernel_matrix_new, indices_new, lower_bounds_new, upper_bounds_new)
            elif certify:
                certificates_zero_equal_components.append(indices_new)

            if is_degenerate:
                return
        return

    recursive_degenerate(non_negative_cocircuits, kernel_matrix, [], lower_bounds, upper_bounds)

    if certify:
        if is_degenerate:
            return is_degenerate, certificate
        return is_degenerate, (certificates_zero_equal_components, certificates_partial_cover, certificate_support_condition)
    return is_degenerate
