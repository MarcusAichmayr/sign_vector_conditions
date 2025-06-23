r"""Utility functions"""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.functions.generalized import sign
from sage.rings.integer_ring import ZZ

from elementary_vectors.utility import is_symbolic
from sign_vectors import sign_vector, zero_sign_vector
from sign_vectors.oriented_matroids import cocircuits_from_matrix


def non_negative_cocircuits_from_matrix(M, dual: bool = True) -> set:
    r"""
    Compute nonnegative cocircuits.

    INPUT:

    - ``M`` -- a matrix with real arguments.

    - ``dual`` -- a boolean (default: ``True``)

    OUTPUT:

    Return a set of nonnegative cocircuits determined by the kernel of ``M``. (default)
    If ``dual`` is false, considers the row space of ``M``.

    EXAMPLES::

        sage: M = matrix([[2, -1, -1, 0]])
        sage: M
        [ 2 -1 -1  0]
        sage: from sign_vectors.oriented_matroids import cocircuits_from_matrix
        sage: cocircuits_from_matrix(M)
        {(0+-0), (--00), (000+), (++00), (0-+0), (+0+0), (-0-0), (000-)}
        sage: from sign_vector_conditions.utility import non_negative_cocircuits_from_matrix
        sage: non_negative_cocircuits_from_matrix(M)
        {(+0+0), (000+), (++00)}
    """
    return set(X for X in cocircuits_from_matrix(M, dual=dual) if X > 0)


def non_negative_covectors_from_matrix(M, dual: bool = True) -> set:
    r"""
    Compute all nonnegative covectors.

    INPUT:

    - ``M`` -- a matrix

    - ``dual`` -- a boolean (default: ``True``)

    OUTPUT:

    Return a set of nonnegative covectors determined by the kernel of ``M``. (default)
    If ``dual`` is false, considers the row space of ``M``.

    EXAMPLES::

        sage: M = matrix([[2, -1, -1, 0]])
        sage: M
        [ 2 -1 -1  0]
        sage: from sign_vectors.oriented_matroids import covectors_from_matrix
        sage: covectors_from_matrix(M)
        {(0000),
         (++-0),
         (--+0),
         (000+),
         (--0+),
         (+-+-),
         (-0-0),
         (+-++),
         (000-),
         (-+-0),
         (--0-),
         (0-+0),
         (++00),
         (--++),
         (+0+0),
         (++--),
         (--00),
         (--+-),
         (-0-+),
         (++-+),
         (---0),
         (0+-0),
         (-+-+),
         (-0--),
         (+++0),
         (-+--),
         (0-++),
         (+-+0),
         (0-+-),
         (++0-),
         (++0+),
         (0+-+),
         (---+),
         (+0+-),
         (+0++),
         (++++),
         (----),
         (0+--),
         (+++-)}
        sage: from sign_vector_conditions.utility import non_negative_covectors_from_matrix
        sage: non_negative_covectors_from_matrix(M)
        {(++0+), (+0+0), (+++0), (000+), (+0++), (++00), (++++)}
    """
    cocircuits = [
        cocircuit
        for cocircuit in cocircuits_from_matrix(M, dual=dual)
        if not cocircuit < 0
    ]

    if not cocircuits:
        raise ValueError("List of cocircuits is empty.")
    output = set()
    new_elements = {zero_sign_vector(M.ncols())}
    while new_elements:
        covector1 = new_elements.pop()
        for covector2 in cocircuits:
            if not covector2 >= 0:
                continue
            if covector2 <= covector1:
                continue
            composition = covector2.compose(covector1)
            if composition not in output and composition >= 0:
                output.add(composition)
                new_elements.add(composition)
    return output


def closure_minors_utility(pairs, positive_only: bool = False, negative_only: bool = False) -> list:
    r"""
    Return whether all products of components are positive (or negative) if first element is nonzero.

    INPUT:

    - ``pairs`` -- an iterable of pairs consisting of a minor and a product

    - ``positive_only`` -- a boolean, considers only positive products if true

    - ``negative_only`` -- a boolean, considers only negative products if true

    OUTPUT:
    Returns either a boolean or sets of conditions on variables occurring in the input.
    If the conditions of one of these sets are satisfied,
    then for all nonzero elements of the first list,
    the product with the corresponding element of the second list is positive.
    (Or all products are negative.)

    TESTS::

        sage: from sign_vector_conditions.utility import closure_minors_utility
        sage: var('a, b, c')
        (a, b, c)
        sage: closure_minors_utility(zip([0, a], [0, a]), positive_only=True)
        [{a == 0}, {a > 0}]
        sage: len(_) # for testing
        2
        sage: closure_minors_utility(zip([c, -1, c], [c, -b, -a * c])) # random
        [{-b > 0, c == 0},
         {-b < 0, c == 0},
         {-b > 0, c > 0, -a*c > 0},
         {-b < 0, c < 0, -a*c < 0}]
        sage: len(_) # for testing
        4
        sage: closure_minors_utility(zip([c, -1, a], [c, -b, -a * c])) # random
        [{-b > 0, a == 0, c == 0},
         {-b < 0, a == 0, c == 0},
         {-b > 0, a == 0, c > 0},
         {-b < 0, a == 0, c < 0},
         {-b > 0, a != 0, c > 0, -a*c > 0},
         {-b < 0, a != 0, c < 0, -a*c < 0},
         {-a*c > 0, c > 0, -b > 0},
         {-a*c < 0, c < 0, -b < 0}]]
        sage: len(_) # for testing
        8
        sage: closure_minors_utility(zip([-1, -1], [-1, -1]))
        True
        sage: closure_minors_utility(zip([-1, 1], [-1, 1]))
        False
        sage: closure_minors_utility(zip([0, 1], [0, 1]))
        True
        sage: closure_minors_utility([(1, 0)])
        False
    """

    def recursive(pairs, zero_expressions, non_zero_expressions):
        r"""Recursive call"""
        pairs = [
            (minor, product)
            for minor, product in pairs
            if not minor in zero_expressions and not minor.is_zero()
        ]
        for minor, _ in pairs:
            if is_symbolic(minor) and not minor in non_zero_expressions:
                yield from recursive(
                    pairs, zero_expressions.union([minor]), non_zero_expressions
                )
                yield from recursive(
                    pairs, zero_expressions, non_zero_expressions.union([minor])
                )

        products = set(
            sign_or_symbolic(
                product.substitute([value == 0 for value in zero_expressions])
            )
            for _, product in pairs
        )
        equalities = set(value == 0 for value in zero_expressions)
        non_equalities = set(
            value != 0 for value in non_zero_expressions if not value in products
        )

        if not negative_only:
            positive_inequalities = set(value > 0 for value in products)
            if True in positive_inequalities:
                positive_inequalities.remove(True)
            yield positive_inequalities.union(equalities).union(non_equalities)

        if not positive_only:
            negative_inequalities = set(value < 0 for value in products)
            if True in negative_inequalities:
                negative_inequalities.remove(True)
            yield negative_inequalities.union(equalities).union(non_equalities)

    output = list(recursive(pairs, set(), set()))
    for conditions in output.copy():
        if False in conditions:
            output.remove(conditions)
    if not output:  # e.g. [1, -1], [1, 1]
        return False
    output = remove_duplicates(output)
    if output == [set()]:  # e.g. [1], [1] or [0], [1]
        return True
    return output


def sign_or_symbolic(expression):
    r"""Return the sign of an expression if defined."""
    if is_symbolic(expression):
        return expression
    return ZZ(sign(expression))


def remove_duplicates(iterable):
    r"""Remove duplicates from a list of iterables."""
    seen = set()
    result = []
    for item in iterable:
        marker = frozenset(item)  # only works if item is an iterable
        if marker in seen:
            continue
        seen.add(marker)
        result.append(item)
    return result


def equal_entries_lists(length: int, indices: list[int]) -> list[list[int]]:
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
    return [[
        1 if i == one_position else (-1 if i == minus_one_position else 0)
        for i in range(length)
    ] for minus_one_position in indices[1:]]


def non_negative_vectors(vectors) -> list:
    r"""
    Return nonnegative vectors.

    INPUT:

    - ``vectors`` -- an iterable of vectors

    OUTPUT:

    Return all vectors of ``vectors`` that are
    - non_negative in each component; or
    - negative in each component. Those will be multiplied by ``-1``; or
    - containing variables such that no opposing signs occur.

    EXAMPLES::

        sage: from sign_vector_conditions.utility import non_negative_vectors
        sage: l = [vector([1, 1, 0, -1]), vector([0, 0, 0, 0]), vector([1, 0, 0, 1])]
        sage: l
        [(1, 1, 0, -1), (0, 0, 0, 0), (1, 0, 0, 1)]
        sage: non_negative_vectors(l)
        [(0, 0, 0, 0), (1, 0, 0, 1)]
        sage: from elementary_vectors import elementary_vectors
        sage: from elementary_vectors.reductions import *
        sage: var('a')
        a
        sage: A = matrix([[a, 0, 0, 0, 1], [0, 1, 0, 0, 1]])
        sage: evs = reduce_vectors(elementary_vectors(A), cancel_factors=True)
        sage: evs
        [(0, 0, 1, 0, 0), (0, 0, 0, 1, 0), (-1, -a, 0, 0, a)]
        sage: non_negative_vectors(evs)
        ...
        UserWarning: Cannot determine sign of symbolic expression, using ``0`` instead.
        [(0, 0, 1, 0, 0), (0, 0, 0, 1, 0), (1, a, 0, 0, -a)]
        sage: assume(a > 0)
        sage: non_negative_vectors(evs)
        [(0, 0, 1, 0, 0), (0, 0, 0, 1, 0)]

    TESTS::

        sage: l = [vector([x, 0, 0])]
        sage: non_negative_vectors(l)
        [(x, 0, 0)]
    """
    result = []
    for element in vectors:
        if sign_vector(element) >= 0:
            result.append(element)
        elif sign_vector(element) < 0:
            result.append(-element)
    return result
