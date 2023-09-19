r"""This module offers utility functions that are used by the other functions of this package."""

#############################################################################
#  Copyright (C) 2023                                                       #
#                Marcus Aichmayr (aichmayr@mathematik.uni-kassel.de)        #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.functions.generalized import sign
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import zero_vector
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR

from sign_vectors import zero_sign_vector
from sign_vectors.oriented_matroids import cocircuits_from_matrix


def positive_cocircuits_from_matrix(M, kernel=True):
    r"""
    Compute the positive cocircuits determined by a given matrix.

    INPUT:

    - ``A`` -- a matrix with real arguments.

    - ``kernel`` -- a boolean (default: ``True``)

    OUTPUT:

    - If ``kernel`` is true, returns a set of positive cocircuits determined by the kernel of the matrix ``A``.

    - If ``kernel`` is false, returns a set of positive cocircuits determined by the rows of the matrix ``A``.

    EXAMPLES::

        sage: M = matrix([[2, -1, -1, 0]])
        sage: M
        [ 2 -1 -1  0]
        sage: from sign_vectors.oriented_matroids import cocircuits_from_matrix
        sage: cocircuits_from_matrix(M)
        {(--00), (000-), (0-+0), (+0+0), (++00), (-0-0), (000+), (0+-0)}
        sage: from sign_vector_conditions.utility import positive_cocircuits_from_matrix
        sage: positive_cocircuits_from_matrix(M)
        {(+0+0), (++00), (000+)}
    """
    return set(X for X in cocircuits_from_matrix(M, kernel=kernel) if X > 0)


def positive_covectors_from_matrix(M, kernel=True):
    r"""
    Use a list of cocircuits to compute all covectors of the corresponding oriented matroid.

    INPUT:

    - ``M`` -- a matrix

    - ``kernel`` -- a boolean (default: ``True``)

    OUTPUT:

    - a list of all covectors of the oriented matroid.

      - If ``kernel`` is true, returns a list of cocircuits determined by the kernel of the matrix ``A``.

      - If ``kernel`` is false, returns a list of cocircuits determined by the rows of the matrix ``A``.

    EXAMPLES::

        sage: M = matrix([[2, -1, -1, 0]])
        sage: M
        [ 2 -1 -1  0]
        sage: from sign_vectors.oriented_matroids import covectors_from_matrix
        sage: covectors_from_matrix(M)
        {(0000),
         (0-+0),
         (0-+-),
         (++++),
         (--0+),
         (000-),
         (----),
         (++0-),
         (+-++),
         ...
         (-+-0),
         (--0-)}
        sage: from sign_vector_conditions.utility import positive_covectors_from_matrix
        sage: positive_covectors_from_matrix(M)
        [(+0++), (+++0), (++++), (+0+0), (++00), (++0+), (000+)]
    """
    cocircuits = [X for X in cocircuits_from_matrix(M, kernel=kernel) if not X < 0]

    if not cocircuits:
        raise ValueError('List of cocircuits is empty.')
    output = set()
    new_elements = {zero_sign_vector(M.ncols())}
    while new_elements:
        Y = new_elements.pop()
        for X in cocircuits:
            if X >= 0 and not X <= Y: # otherwise Z = X.compose(Y) = Y in F
                Z = X.compose(Y)
                if Z not in output and Z >= 0:
                    output.add(Z)
                    new_elements.add(Z)
    return list(output)


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
        sage: from sign_vector_conditions.utility import max_minors_prod
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

        sage: from sign_vector_conditions.utility import entries_non_negative
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

        sage: from sign_vector_conditions.utility import entries_non_positive
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

        sage: from sign_vector_conditions.utility import entries_non_negative_or_non_positive
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


def condition_on_products(list1, list2):
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
    
        sage: from sign_vector_conditions.utility import condition_on_products
        sage: var('a,b,c')
        (a, b, c)
        sage: condition_on_products([0, a], [1, 1])
        [{a == 0}, {a > 0}, {a < 0}]
        sage: len(_)
        3
        sage: condition_on_products([c, -1, c], [1, b, -a]) # random
        [{-b > 0, c == 0},
         {-b < 0, c == 0},
         {-b > 0, c > 0, -a*c > 0},
         {-b < 0, c < 0, -a*c < 0}]
        sage: len(_)
        4
        sage: condition_on_products([c, -1, a], [1, b, -c]) # random
        [{-b > 0, a == 0, c == 0},
         {-b < 0, a == 0, c == 0},
         {-b > 0, a == 0, c > 0},
         {-b < 0, a == 0, c < 0},
         {-b > 0, a != 0, c > 0, -a*c > 0},
         {-b < 0, a != 0, c < 0, -a*c < 0}]
        sage: len(_[4])
        4
        sage: condition_on_products([-1, -1], [1, 1])
        True
        sage: condition_on_products([-1, 1], [1, 1])
        False
        sage: condition_on_products([0, 1], [1, 1])
        True
        sage: condition_on_products([1], [0])
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
                yield from rec(list1, list2, zero_expressions.union([elem1]), non_zero_expressions)
                yield from rec(list1, list2, zero_expressions, non_zero_expressions.union([elem1]))

        products = set(
            sign_or_symbolic((elem1 * elem2).substitute([value == 0 for value in zero_expressions]))
        for elem1, elem2 in pairs)

        equalities = set(value == 0 for value in zero_expressions)
        non_equalities = set(value != 0 for value in non_zero_expressions if not value in products)

        positive_inequalities = set(value > 0 for value in products)
        negative_inequalities = set(value < 0 for value in products)

        if True in positive_inequalities:
            positive_inequalities.remove(True)
        if True in negative_inequalities:
            negative_inequalities.remove(True)

        yield positive_inequalities.union(equalities).union(non_equalities)
        yield negative_inequalities.union(equalities).union(non_equalities)

    output = list(rec(list1, list2, set(), set()))
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


def sign_or_symbolic(expression):
    r"""Return the sign if defined"""
    try:
        return ZZ(sign(expression))
    except TypeError:
        return expression


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


def equal_entries_lists(length, indices):
    r"""
    Return a list of lists such that the corresponding kernel matrix has equal entries.
    
    INPUT:
    
    - ``length`` -- an integer
    
    - ``indices`` -- a list of integers
    
    OUTPUT:
    a list of lists
    
    EXAMPLES::
    
        sage: from sign_vector_conditions.utility import equal_entries_lists
        sage: equal_entries_lists(5, [1, 2, 3])
        [[0, 1, -1, 0, 0], [0, 1, 0, -1, 0]]
        sage: equal_entries_lists(3, [0])
        []
        sage: equal_entries_lists(3, [0, 1])
        [[1, -1, 0]]
    """
    if len(indices) < 2:
        return []

    one_position = indices[0]
    return [
        [1 if i == one_position else
            (-1 if i == minus_one_position else 0)
            for i in range(length)
        ] for minus_one_position in indices[1:]
    ]
