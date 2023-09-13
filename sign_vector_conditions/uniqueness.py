r"""
Uniqueness of equilibria

EXAMPLES::

    sage: from sign_vector_conditions import *

We define some matrices::

    sage: W = matrix([[1,0,-1],[0,1,-1]])
    sage: W
    [ 1  0 -1]
    [ 0  1 -1]
    sage: Wt = matrix([[1,0,-1],[0,1,0]])
    sage: Wt
    [ 1  0 -1]
    [ 0  1  0]

Therefore, we obtain the following exponential map::

    sage: var('x1, x2, c1, c2, c3')
    (x1, x2, c1, c2, c3)
    sage: c = [c1, c2, c3]
    sage: Fc = f_exp(W, Wt, c)
    sage: Fc(x1, x2)
    (-c3*e^(-x1) + c1*e^x1, -c3*e^(-x1) + c2*e^x2)

We want to check whether this map is injective for each vector ``c > 0``.
For this purpose, we compute the corresponding oriented matroids::

    sage: from sign_vectors.oriented_matroids import *
    sage: cvW = covectors_from_matrix(W, kernel=False, algorithm='fe')
    sage: cvW
    {(000),
     (+-+),
     (-+0),
     (-+-),
     (0-+),
     (-++),
     (+--),
     (0+-),
     (+-0),
     (+0-),
     (-0+),
     (--+),
     (++-)}
    sage: cvWt = covectors_from_matrix(Wt, kernel=True, algorithm='fe')
    sage: cvWt
    {(000), (+0+), (-0-)}

The intersection of these oriented matroids consists only of the zero sign vector.
We can compute the intersection directly by applying the built in method intersection::

    sage: set(cvW).intersection(cvWt)
    {(000)}

Therefore, the corresponding exponential map is injective.
This means that the chemical reaction network has at most one solution.
We can also check this condition in the following way::

    sage: condition_uniqueness_signvectors(W, Wt)
    True

There is another way to check injectivity for exponential maps
that involves the computation of maximal minors of the corresponding matrices::

    sage: m1 = W.minors(2)
    sage: m1
    [1, -1, 1]
    sage: m2 = Wt.minors(2)
    sage: m2
    [1, 0, 1]

We multiply those minors component-wise::

    sage: [m1[i]*m2[i] for i in range(len(m1))]
    [1, 0, 1]

Since all arguments are greater or equal zero, the map is injective.
We can also check this condition by applying the following function
from this package::

    sage: condition_uniqueness_minors(W, Wt)
    True

Now, we consider another example::

    sage: W = matrix([[1,0,-1],[0,1,-1]])
    sage: W
    [ 1  0 -1]
    [ 0  1 -1]
    sage: Wt = matrix([[1,0,-1],[0,1,1]])
    sage: Wt
    [ 1  0 -1]
    [ 0  1  1]

Next, we compute the corresponding oriented matroids::

    sage: covectors_from_matrix(W, kernel=False, algorithm='fe', separate=True)
    [{(000)},
     {(-+0), (0-+), (0+-), (+-0), (+0-), (-0+)},
     {(+-+), (-+-), (--+), (-++), (++-), (+--)}]
    sage: covectors_from_matrix(Wt, kernel=True, algorithm='fe', separate=True)
    [{(000)}, {(+-+), (-+-)}]

Now, we check the condition from before::

    sage: condition_uniqueness_signvectors(W, Wt)
    False

Therefore, the corresponding exponential map is not injective.
Furthermore, we obtain the following minors::

    sage: m1 = W.minors(2)
    sage: m1
    [1, -1, 1]
    sage: m2 = Wt.minors(2)
    sage: m2
    [1, 1, 1]
    sage: [m1[i]*m2[i] for i in range(len(m1))]
    [1, -1, 1]

There are positive and negative elements in the resulting list.
Hence, this condition also states that the map is not injective::

    sage: condition_uniqueness_minors(W, Wt)
    False

Finally, we consider an example with variables::

    sage: var('a,b')
    (a, b)
    sage: W = matrix([[1,0,-1],[0,1,-1]])
    sage: W
    [ 1  0 -1]
    [ 0  1 -1]
    sage: Wt = matrix([[1,0,a],[0,1,b]])
    sage: Wt
    [1 0 a]
    [0 1 b]

The matrix ``Wt`` contains variables :math:`a, b \in \mathbb{R}`.
Consequently, we cannot compute the corresponding oriented matroids.
On the other hand, we can still compute the minors of ``W`` and ``Wt``, that is::

    sage: m1 = W.minors(2)
    sage: m1
    [1, -1, 1]
    sage: m2 = Wt.minors(2)
    sage: m2
    [1, b, -a]
    sage: [m1[i] * m2[i] for i in range(len(m1))]
    [1, -b, -a]

Therefore, there is at most one equilibrium if and only if
:math:`a \leq 0` and :math:`b \leq 0`.
The function :func:`~condition_uniqueness_minors` also works for matrices with symbolic entries.
In this case, it returns a system of inequalities::

    sage: condition_uniqueness_minors(W, Wt)
    {-b >= 0, -a >= 0}
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

from .utility import normalize
from sign_vectors.oriented_matroids import covectors_from_matrix
from sage.rings.real_mpfr import RR


def condition_uniqueness_signvectors(W, Wt):
    r"""
    Use sign vectors to verify that the chemical reaction network has at most one equilibrium.

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    OUTPUT:
    Returns whether the intersection of the oriented matroids corresponding to
    ``W`` and ``right_kernel(Wt)`` consists of the zero sign vector only.

    EXAMPLES::

        sage: from sign_vector_conditions import *
        sage: W = matrix([[1,0,-1],[0,1,-1]])
        sage: W
        [ 1  0 -1]
        [ 0  1 -1]
        sage: Wt = matrix([[1,0,-1],[0,1,0]])
        sage: Wt
        [ 1  0 -1]
        [ 0  1  0]
        sage: condition_uniqueness_signvectors(W, Wt)
        True
        sage: W = matrix([[1,0,-1],[0,1,-1]])
        sage: W
        [ 1  0 -1]
        [ 0  1 -1]
        sage: Wt = matrix([[1,0,-1],[0,1,1]])
        sage: Wt
        [ 1  0 -1]
        [ 0  1  1]
        sage: condition_uniqueness_signvectors(W, Wt)
        False

    TESTS::

        sage: from sign_vector_conditions.uniqueness import condition_uniqueness_signvectors
        sage: A = identity_matrix(3)
        sage: B = A # kernel of B is empty
        sage: condition_uniqueness_signvectors(A, B)
        True
    """
    if W.ncols() != Wt.ncols():
        raise ValueError('Matrices have different number of columns.')

    return len(
        normalize(covectors_from_matrix(W, kernel=False)).intersection(
            normalize(covectors_from_matrix(Wt, kernel=True))
        )
    ) == 1


def max_minors_prod(A, B):
    r"""
    Multiply the maximal minors of two matrices component-wise.

    EXAMPLES::

        sage: A = matrix([[1, 0, 1], [0, 1, 2]])
        sage: B = matrix([[1, 0, -1], [0, 1, 3]])
        sage: A.minors(2)
        [1, 2, -1]
        sage: B.minors(2)
        [1, 3, 1]
        sage: from sign_vector_conditions.uniqueness import max_minors_prod
        sage: max_minors_prod(A, B)
        [1, 6, -1]

    TESTS::

        sage: A = matrix([[1, 0, 1], [0, 1, 2], [0, 1, 2]])
        sage: B = matrix([[1, 0, -1], [0, 1, 3]])
        sage: max_minors_prod(A, B)
        [1, 6, -1]
        sage: A = matrix([[1, 0, 1], [0, 1, 2], [0, 0, 1]])
        sage: max_minors_prod(A, B)
        Traceback (most recent call last):
        ...
        ValueError: Matrices must have same rank and number of columns.
    """
    A = A.matrix_from_rows(A.pivot_rows())
    B = B.matrix_from_rows(B.pivot_rows())
    if A.dimensions() != B.dimensions():
        raise ValueError('Matrices must have same rank and number of columns.')

    return [m1 * m2 for m1, m2 in zip(A.minors(A.nrows()), B.minors(A.nrows()))]


def compare_all(iterable, relation):
    r"""
    Check whether all entries satisfy the relation.

    This is an auxiliary function used by :func:`entries_non_negative` and :func:`entries_non_positive`.
    """
    output = set()
    for value in iterable:
        if relation(value) is False:
            return False
        if relation(value) is not True:  # if variables occur, we will return the expression
            output.add(relation(value))
    if len(output) == 0:
        return True
    return output


def entries_non_negative(iterable):
    r"""
    Check whether all entries are non-negative.

    INPUT:

    - ``iterable`` -- an iterable

    OUTPUT:
    a boolean or symbolic expression

    EXAMPLES::

        sage: from sign_vector_conditions.uniqueness import entries_non_negative
        sage: entries_non_negative([0, 5, 1])
        True
        sage: entries_non_negative([0, 0])
        True
        sage: entries_non_negative([0, -5])
        False
        sage: entries_non_negative([x, x^2 + 1, -1, 5])
        False
        sage: entries_non_negative([x, x^2 + 1]) # random
        {x^2 + 1 >= 0, x >= 0}
    """
    def relation(value):
        try:
            return RR(value) >= 0
        except TypeError:
            return value >= 0
    return compare_all(iterable, relation)


def entries_non_positive(iterable):
    r"""
    Check whether all entries are non-positive.

    INPUT:

    - ``iterable`` -- an iterable

    OUTPUT:
    a boolean or symbolic expression

    EXAMPLES::

        sage: from sign_vector_conditions.uniqueness import entries_non_positive
        sage: entries_non_positive([0, 5, 1])
        False
        sage: entries_non_positive([0, 0])
        True
        sage: entries_non_positive([0, -5])
        True
        sage: entries_non_positive([x, x^2 + 1, -1, 5])
        False
        sage: entries_non_positive([x, x^2 + 1])
        {x^2 + 1 <= 0, x <= 0}
    """
    def relation(value):
        try:
            return RR(value) <= 0
        except TypeError:
            return value <= 0
    return compare_all(iterable, relation)


def entries_non_negative_or_non_positive(iterable):
    r"""
    Return whether each component of a given vector is non-negative or non-positive.

    INPUT:

    - ``iterable`` -- an iterable of numbers or symbolic expressions.

    OUTPUT:
    Returns true if either each element of ``iterable`` is greater than or equal to zero or less than or equal to zero.
    Supports symbolic expressions.

    Depending on the input, the output can have several appearances:

    - a boolean

        - if true, no symbolic expressions have occurred and either each entry of ``iterable`` is greater or equal zero
          or each entry is less or equal zero.

        - if false, there are entries with opposing signs or every entry is zero.

    - a set of inequalities

        - if all inequalities are satisfied,
          then either each element of ``iterable`` is greater or equal zero or less or equal zero.

    - a list of two sets of inequalities

        - if the inequalities of exactly one of these sets are satisfied,
          then either each element of ``iterable`` is greater or equal zero or less or equal zero.

    EXAMPLES::

        sage: from sign_vector_conditions.uniqueness import entries_non_negative_or_non_positive
        sage: entries_non_negative_or_non_positive([0, 5, 1])
        True
        sage: entries_non_negative_or_non_positive([0, 0])
        False
        sage: entries_non_negative_or_non_positive([0, -5])
        True
        sage: entries_non_negative_or_non_positive([x, x^2 + 1, -1, 5])
        False
        sage: entries_non_negative_or_non_positive([x, x^2 + 1]) # random
        [{x^2 + 1 >= 0, x >= 0}, {x^2 + 1 <= 0, x <= 0}]
    """
    entries_nn = entries_non_negative(iterable)
    entries_np = entries_non_positive(iterable)

    if entries_nn is True and entries_np is True:  # all entries are zero
        return False
    if entries_nn is False and entries_np is False:  # mixed signs
        return False
    if entries_nn is False:
        return True if entries_np is True else entries_np
    if entries_np is False:
        return True if entries_nn is True else entries_nn
    return [entries_nn, entries_np]


def condition_uniqueness_minors(W, Wt):
    r"""
    Use maximal minors to verify that the chemical reaction network has at most one equilibrium.

    INPUT:

    - ``W`` -- a matrix

    - ``Wt`` -- a matrix with the same dimensions as ``W``

    OUTPUT:
    Returns whether the products of the corresponding maximal minors of the
    matrices ``W`` and ``Wt`` have the same sign.

    Returns a boolean or a symbolic expression if variables occur.

    EXAMPLES::

        sage: from sign_vector_conditions import *
        sage: W = matrix([[1, 0, -1], [0, 1, -1]])
        sage: W
        [ 1  0 -1]
        [ 0  1 -1]
        sage: Wt = matrix([[1, 0, -1], [0, 1, 0]])
        sage: Wt
        [ 1  0 -1]
        [ 0  1  0]
        sage: condition_uniqueness_minors(W, Wt)
        True
        sage: W = matrix([[1, 0, -1], [0, 1, -1]])
        sage: W
        [ 1  0 -1]
        [ 0  1 -1]
        sage: Wt = matrix([[1, 0, -1], [0, 1, 1]])
        sage: Wt
        [ 1  0 -1]
        [ 0  1  1]
        sage: condition_uniqueness_minors(W, Wt)
        False
        sage: var('a, b')
        (a, b)
        sage: W = matrix([[1, 0, -1], [0, 1, -1]])
        sage: W
        [ 1  0 -1]
        [ 0  1 -1]
        sage: Wt = matrix([[1, 0, a], [0, 1, b]])
        sage: Wt
        [1 0 a]
        [0 1 b]
        sage: condition_uniqueness_minors(W, Wt)
        {-b >= 0, -a >= 0}
    """
    return entries_non_negative_or_non_positive(max_minors_prod(W, Wt))
