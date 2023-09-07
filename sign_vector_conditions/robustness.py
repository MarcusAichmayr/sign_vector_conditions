r"""
Robustness of existence and uniqueness of equilibria

EXAMPLES::

    sage: from sign_vector_conditions import *

Let us consider the following matrices::

    sage: W = matrix([[1,0,-1,0],[0,1,0,0]])
    sage: W
    [ 1  0 -1  0]
    [ 0  1  0  0]
    sage: Wt = matrix([[1,0,-1,0],[0,1,-1,1]])
    sage: Wt
    [ 1  0 -1  0]
    [ 0  1 -1  1]

Therefore, we consider the following exponential map::

    sage: var('x1, x2, c1, c2, c3, c4')
    (x1, x2, c1, c2, c3, c4)
    sage: c = [c1, c2, c3, c4]
    sage: Fc = f_exp(W, Wt, c)
    sage: Fc(x1, x2)
    (c1*e^x1 - c3*e^(-x1 - x2), c2*e^x2)

To check, whether this map is a diffeomorphism for all ``c > 0``
and all small perturbations of ``Wt``,
we consider the topes of the corresponding oriented matroids::

    sage: from sign_vectors.oriented_matroids import *
    sage: topes_from_matrix(W, kernel=True)
    {(+0+-), (+0++), (-0--), (-0-+)}
    sage: topes_from_matrix(Wt, kernel=True)
    {(---+), (-+--), (++++), (----), (+-++), (+++-)}

One can see that for every tope ``X`` of the oriented matroid corresponding to ``W`` there is a
tope ``Y`` corresponding to ``Wt`` such that ``X`` conforms to ``Y``.
Therefore, the exponential map is a diffeomorphism for all ``c > 0``
and all small perturbations of ``Wt``.
The package offers a function that checks this condition directly::

    sage: condition_closure_sign_vectors(W, Wt)
    True

There is an equivalent condition.
To verify it, we compute the maximal minors of the two matrices::

    sage: W.minors(2)
    [1, 0, 0, 1, 0, 0]
    sage: Wt.minors(2)
    [1, -1, 1, 1, 0, -1]

From the output, we see whenever a minor of ``W`` is non-zero,
the corresponding minor of ``Wt`` has the same sign.
Hence, this condition is fulfilled.
This condition can also be checked directly with the package::

    sage: condition_closure_minors(W, Wt)
    True

Now, we consider matrices with variables::

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

We cannot check the first condition since there are variables in ``W`` and ``Wt``.
Therefore, we want to obtain equations on the variables ``a``, ``b``, ``c``
such that this condition is satisfied.
First, we compute the minors of the matrices::

    sage: W.minors(2)
    [c, -1, c]
    sage: Wt.minors(2)
    [1, b, -a]

The function from the package supports symbolic matrices as input.
In this case, we obtain the following equations on the variables::

    sage: condition_closure_minors(W, Wt) # random
    [{-b > 0, c == 0},
     {-b < 0, c == 0},
     {-b > 0, c > 0, -a*c > 0},
     {-b < 0, c < 0, -a*c < 0}]

Thus, there are four possibilities to set the variables:
From the first two sets of conditions, we see that the closure condition is satisfied
if ``c`` is zero and ``b`` is non-zero.
The closure condition is also satisfied if ``a`` and ``b`` are negative and ``c`` is positive
or if ``a`` and ``b`` are positive and ``c`` is negative.
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
from sign_vectors.oriented_matroids import topes_from_matrix
from sage.rings.real_mpfr import RR
from sage.misc.flatten import flatten


def condition_closure_sign_vectors(W, Wt):
    r"""
    Closure condition on sign vectors for robustness

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    OUTPUT:
    Return whether the oriented matroid corresponding to ``W`` is a subset of the closure
    of the oriented matroid corresponding to ``Wt``.
    """
    tW = topes_from_matrix(W, kernel=True)
    tWt = topes_from_matrix(Wt, kernel=True)
    tWn = normalize(tW)
    # Do not normalize second list of topes: 0++0, +--+
    for X in tWn:
        val = True
        for Y in tWt:
            if X <= Y:
                val = False
                break
        if val:
            return False
    return True


def condition_closure_minors(W, Wt):
    r"""
    Closure condition on maximal minors for robustness

    INPUT:

    - ``W`` -- a matrix

    - ``Wt`` -- a matrix with the same dimensions as ``W``

    OUTPUT:
    Returns either a boolean or sets of conditions on variables occurring in the input.
    """
    if W.dimensions() != Wt.dimensions():
        raise ValueError('Matrices must have same dimensions.')

    return conditions_on_products(W.minors(W.nrows()), Wt.minors(W.nrows()))


def conditions_on_products(list1, list2):
    r"""
    Return whether all products of components are positive (or negative) if first element is non-zero
    
    INPUT:
    Two lists of the same length.
    
    OUTPUT:
    Returns either a boolean or sets of conditions on variables occurring in the input.
    If the conditions of one of these sets are satisfied,
    then for all non-zero elements of the first list,
    the product with the corresponding element of the second list is positive.
    (Or all products are negative.)
    
    .. SEEALSO::

        :func:`~condition_closure_minors`
    
    TESTS::
    
        sage: from sign_vector_conditions.robustness import conditions_on_products
        sage: var('a,b,c')
        (a, b, c)
        sage: conditions_on_products([0, a], [1, 1])
        [{a == 0}, {a > 0}, {a < 0}]
        sage: len(_)
        3
        sage: conditions_on_products([c, -1, c], [1, b, -a]) # random
        [{-b > 0, c == 0},
         {-b < 0, c == 0},
         {-b > 0, c > 0, -a*c > 0},
         {-b < 0, c < 0, -a*c < 0}]
        sage: len(_)
        4
        sage: conditions_on_products([c, -1, a], [1, b, -c]) # random
        [{-b > 0, a == 0, c == 0},
         {-b < 0, a == 0, c == 0},
         {-b > 0, a == 0, c > 0},
         {-b < 0, a == 0, c < 0},
         {-b > 0, a != 0, c > 0, -a*c > 0},
         {-b < 0, a != 0, c < 0, -a*c < 0}]
        sage: len(_[4])
        4
        sage: conditions_on_products([-1, -1], [1, 1])
        True
        sage: conditions_on_products([-1, 1], [1, 1])
        False
        sage: conditions_on_products([0, 1], [1, 1])
        True
        sage: conditions_on_products([1], [0])
        False
    """
    def rec(list1, list2, zero_expressions, non_zero_expressions):
        r"""Recursive call"""
        pairs = [
            (elem1, elem2) for elem1, elem2 in zip(list1, list2)
            if not elem1.is_zero() and not elem1 in zero_expressions
        ]

        for elem1, _ in pairs:
            if is_symbolic(elem1) and not elem1 in non_zero_expressions:
                return [
                    rec(list1, list2, zero_expressions.union([elem1]), non_zero_expressions),
                    rec(list1, list2, zero_expressions, non_zero_expressions.union([elem1]))
                ]

        products = set(
            substitute_and_simplify(elem1 * elem2, [value == 0 for value in zero_expressions])
        for elem1, elem2 in pairs)

        equalities = set(value == 0 for value in zero_expressions)
        non_equalities = set(value != 0 for value in non_zero_expressions if not value in products)

        positive_inequalities = set(value > 0 for value in products)
        negative_inequalities = set(value < 0 for value in products)

        if True in positive_inequalities:
            positive_inequalities.remove(True)
        if True in negative_inequalities:
            negative_inequalities.remove(True)

        return [
            positive_inequalities.union(equalities).union(non_equalities),
            negative_inequalities.union(equalities).union(non_equalities)
        ]

    output = flatten(rec(list1, list2, set(), set()))
    for conditions in output.copy():
        if False in conditions:
            output.remove(conditions)
    if not output: # e.g. [1, -1], [1, 1]
        return False
    output = remove_duplicates(output)
    if output == [set()]: # e.g. [1], [1] or [0], [1]
        return True
    return output


def is_symbolic(value):
    r"""Return whether this element is symbolic"""
    try:
        return value.is_symbol()
    except AttributeError:
        return False


def substitute_and_simplify(expression, *args):
    r"""Substitute arguments and make result a real number if possible"""
    value = expression.substitute(*args)
    try:
        return RR(value)
    except TypeError:
        return value


def remove_duplicates(iterable):
    r"""Remove duplicates from a list of iterables"""
    seen = set()
    result = []
    for item in iterable:
        marker = frozenset(item) # only works if item is an iterable
        if marker in seen:
            continue
        seen.add(marker)
        result.append(item)
    return result
