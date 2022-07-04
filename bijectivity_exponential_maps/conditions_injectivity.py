r"""
In this module, we check injectivity for exponential maps by verifying conditions from [MHR19]_.

EXAMPLES::

    sage: from bijectivity_exponential_maps import *

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
    sage: cvW = covectors_from_matrix(W, algorithm='fe')
    sage: cvW
    [(000),
     (0-+),
     (+-0),
     (-0+),
     (-+0),
     (+0-),
     (0+-),
     (+-+),
     (-++),
     (+--),
     (--+),
     (-+-),
     (++-)]
    sage: cvWt = covectors_from_matrix(Wt, kernel=True, algorithm='fe')
    sage: cvWt
    [(000), (+0+), (-0-)]

The intersection of these oriented matroids consists only of the zero sign vector.
We can compute the intersection directly by applying the built in method intersection::

    sage: set(cvW).intersection(cvWt)
    {(000)}

Therefore, the corresponding exponential map is injective.
We can also check this condition in the following way::

    sage: cond_inj_intersection(W, Wt)
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

    sage: cond_inj_minors(W, Wt)
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

    sage: covectors_from_matrix(W, algorithm='fe', separate=True)
    [[(000)],
     [(0-+), (+-0), (-0+), (-+0), (+0-), (0+-)],
     [(+-+), (-++), (+--), (--+), (-+-), (++-)]]
    sage: covectors_from_matrix(Wt, kernel=True, algorithm='fe', separate=True)
    [[(000)], [(+-+), (-+-)]]

Now, we check the condition from before::

    sage: cond_inj_intersection(W, Wt)
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

    sage: cond_inj_minors(W, Wt)
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
    sage: [m1[i]*m2[i] for i in range(len(m1))]
    [1, -b, -a]

Therefore, the corresponding exponential map is injective if and only if
:math:`a \leq 0` and :math:`b \leq 0`.
The function :func:`~cond_inj_minors` also works for matrices with symbolic entries.
In this case, it returns a system of inequalities::

    sage: cond_inj_minors(W, Wt)
    [-b >= 0, -a >= 0]
"""

#############################################################################
#  Copyright (C) 2022                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from .utility import normalize
from sign_vectors.oriented_matroids import covectors_from_matrix
from sage.rings.real_mpfr import RR  # used for casting


def cond_inj_intersection(W, Wt):
    r"""
    Return whether the intersection of two oriented matroids consists of the zero sign vector only.

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    OUTPUT:
    Returns whether the intersection of the oriented matroids corresponding to
    ``W`` and ``right_kernel(Wt)`` consists of the zero sign vector only.

    Returns a boolean.

    EXAMPLES::

        sage: from bijectivity_exponential_maps import *
        sage: W = matrix([[1,0,-1],[0,1,-1]])
        sage: W
        [ 1  0 -1]
        [ 0  1 -1]
        sage: Wt = matrix([[1,0,-1],[0,1,0]])
        sage: Wt
        [ 1  0 -1]
        [ 0  1  0]
        sage: cond_inj_intersection(W, Wt)
        True
        sage: W = matrix([[1,0,-1],[0,1,-1]])
        sage: W
        [ 1  0 -1]
        [ 0  1 -1]
        sage: Wt = matrix([[1,0,-1],[0,1,1]])
        sage: Wt
        [ 1  0 -1]
        [ 0  1  1]
        sage: cond_inj_intersection(W, Wt)
        False

    TESTS::

        sage: from bijectivity_exponential_maps.conditions_injectivity import cond_inj_intersection
        sage: A = identity_matrix(3)
        sage: B = A # kernel of B is empty
        sage: cond_inj_intersection(A, B)
        True
    """
    if W.ncols() != Wt.ncols():
        raise ValueError('Matrices have different number of columns.')

    cvW = covectors_from_matrix(W)
    try:
        cvWt = covectors_from_matrix(Wt, kernel=True)
    except ValueError:  # kernel matrix might be empty, no elementary vectors
        return True
    cvWn = normalize(cvW)
    cvWtn = normalize(cvWt)
    for X in cvWtn[1:]:  # first element is zero vector
        if X in cvWn[1:]:  # first element is zero vector
            return False
    return True


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
        sage: from bijectivity_exponential_maps.conditions_injectivity import max_minors_prod
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
    A1 = A.matrix_from_rows(A.pivot_rows())
    B1 = B.matrix_from_rows(B.pivot_rows())
    if A1.dimensions() != B1.dimensions():
        raise ValueError('Matrices must have same rank and number of columns.')
    r = A1.nrows()
    mA = A1.minors(r)
    mB = B1.minors(r)

    return [mA[i]*mB[i] for i in range(len(mA))]


def compare_all(v, rel):
    r"""
    Check whether all entries satisfy the relation.

    This is an auxiliary function used by :func:`geq` and :func:`leq`.
    """
    l = []
    for a in v:
        if rel(a) is False:
            return False
        elif rel(a) is not True:  # if variables occur, we will return the expression
            l.append(rel(a))
    if len(l) == 0:
        return True
    else:
        return l


def geq(v):
    r"""
    Check whether all entries are non-negative.

    INPUT:

    - ``v`` -- an iterable

    OUTPUT:
    a boolean or symbolic expression

    EXAMPLES::

        sage: from bijectivity_exponential_maps.conditions_injectivity import geq
        sage: geq([0, 5, 1])
        True
        sage: geq([0, 0])
        True
        sage: geq([0, -5])
        False
        sage: geq([x, x^2 + 1, -1, 5])
        False
        sage: geq([x, x^2 + 1])
        [x >= 0, x^2 + 1 >= 0]
    """
    def rel(a):
        try:
            return RR(a) >= 0  # cast to a real number
        except TypeError:
            return a >= 0
    return compare_all(v, rel)


def leq(v):
    r"""
    Check whether all entries are non-positive.

    INPUT:

    - ``v`` -- an iterable

    OUTPUT:
    a boolean or symbolic expression

    EXAMPLES::

        sage: from bijectivity_exponential_maps.conditions_injectivity import leq
        sage: leq([0, 5, 1])
        False
        sage: leq([0, 0])
        True
        sage: leq([0, -5])
        True
        sage: leq([x, x^2 + 1, -1, 5])
        False
        sage: leq([x, x^2 + 1])
        [x <= 0, x^2 + 1 <= 0]
    """
    def rel(a):
        try:
            return RR(a) <= 0  # cast to a real number
        except TypeError:
            return a <= 0
    return compare_all(v, rel)


def geq_leq(v):
    r"""
    Return whether each component of a given vector is non-negative or non-positive.

    INPUT:

    - ``v`` -- a list of numbers or symbolic expressions.

    OUTPUT:
    Returns true if either each element of ``v`` is greater than or equal to zero or less than or equal to zero.
    Supports symbolic expressions.

    Depending on the input, the output can have several appearances:

    - a boolean

        - if true, no symbolic expressions have occurred and either each entry of ``v`` is greater or equal zero
          or each entry is less or equal zero.

        - if false, there are entries with opposing signs or every entry is zero.

    - a list of inequalities

        - if all inequalities are satisfied,
          then either each element of ``v`` is greater or equal zero or less or equal zero.

    - a list of two lists of inequalities

        - if the inequalities of exactly one of these lists are satisfied,
          then either each element of ``v`` is greater or equal zero or less or equal zero.

    EXAMPLES::

        sage: from bijectivity_exponential_maps.conditions_injectivity import geq_leq
        sage: geq_leq([0, 5, 1])
        True
        sage: geq_leq([0, 0])
        False
        sage: geq_leq([0, -5])
        True
        sage: geq_leq([x, x^2 + 1, -1, 5])
        False
        sage: geq_leq([x, x^2 + 1])
        [[x >= 0, x^2 + 1 >= 0], [x <= 0, x^2 + 1 <= 0]]
    """
    ge = geq(v)
    le = leq(v)

    if ge is True and le is True:  # all entries are zero
        return False
    elif ge is False and le is False:  # mixed signs
        return False
    elif ge is False:
        if le is True:
            return True
        else:
            return le
    elif le is False:
        if ge is True:
            return True
        else:
            return ge
    else:
        return [ge, le]


# Corollary 4 condition 2
def cond_inj_minors(W, Wt):
    r"""
    Return whether the products of maximal minors have the same sign.

    INPUT:

    - ``W`` -- a matrix

    - ``Wt`` -- a matrix with the same dimensions as ``W``

    OUTPUT:
    Returns whether the products of the corresponding maximal minors of the
    matrices ``W`` and ``Wt`` have the same sign.

    Returns a boolean or a symbolic expression if variables occur.

    EXAMPLES::

        sage: from bijectivity_exponential_maps import *
        sage: W = matrix([[1,0,-1],[0,1,-1]])
        sage: W
        [ 1  0 -1]
        [ 0  1 -1]
        sage: Wt = matrix([[1,0,-1],[0,1,0]])
        sage: Wt
        [ 1  0 -1]
        [ 0  1  0]
        sage: cond_inj_minors(W, Wt)
        True
        sage: W = matrix([[1,0,-1],[0,1,-1]])
        sage: W
        [ 1  0 -1]
        [ 0  1 -1]
        sage: Wt = matrix([[1,0,-1],[0,1,1]])
        sage: Wt
        [ 1  0 -1]
        [ 0  1  1]
        sage: cond_inj_minors(W, Wt)
        False
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
        sage: cond_inj_minors(W, Wt)
        [-b >= 0, -a >= 0]
    """
    m = max_minors_prod(W, Wt)
    return geq_leq(m)
