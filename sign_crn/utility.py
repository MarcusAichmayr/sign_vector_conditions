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

from collections.abc import Generator
from typing import Iterator

from sage.functions.generalized import sign
from sage.matrix.constructor import Matrix
from sage.misc.mrange import cartesian_product_iterator
from sage.modules.free_module_element import vector, zero_vector
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import minus_infinity, Infinity

from elementary_vectors import CircuitEnumerator
from elementary_vectors.utility import is_constant
from certlin import Interval, Intervals
from sign_vectors import sign_vector, zero_sign_vector, SignVector, OrientedMatroid


def non_negative_cocircuits_from_matrix(matrix: Matrix) -> set[SignVector]:
    r"""
    Compute nonnegative cocircuits.

    OUTPUT:

    Return a set of nonnegative cocircuits determined by the kernel of ``matrix``.

    EXAMPLES::

        sage: M = matrix([[1, 0, 2, 0], [0, 1, -1, 0], [0, 0, 0, 1]])
        sage: from sign_vectors import *
        sage: from sign_crn.utility import non_negative_cocircuits_from_matrix
        sage: OrientedMatroid(M).cocircuits()
        {(0+-0), (--00), (0-+0), (000+), (++00), (+0+0), (-0-0), (000-)}
        sage: non_negative_cocircuits_from_matrix(M)
        {(+0+0), (000+), (++00)}
    """
    return set(X for X in OrientedMatroid(matrix).cocircuits() if X > 0)


def non_negative_circuits_from_matrix(matrix: Matrix) -> set[SignVector]:
    r"""
    Compute all nonnegative circuits.

    OUTPUT:

    Return a set of nonnegative circuits determined by the kernel of ``matrix``.

    EXAMPLES::

        sage: M = matrix([[2, -1, -1, 0]])
        sage: from sign_vectors import *
        sage: from sign_crn.utility import non_negative_circuits_from_matrix
        sage: OrientedMatroid(M).circuits()
        {(0+-0), (--00), (000+), (++00), (0-+0), (+0+0), (-0-0), (000-)}
        sage: non_negative_circuits_from_matrix(M)
        {(+0+0), (000+), (++00)}
    """
    return set(X for X in OrientedMatroid(matrix).circuits() if X > 0)


def non_negative_covectors_from_cocircuits(cocircuits: set[SignVector], length: int) -> set[SignVector]:
    r"""Compute all nonnegative covectors from a set of cocircuits."""
    if not cocircuits:
        raise ValueError("List of cocircuits is empty.")
    output = {zero_sign_vector(length)}
    new_elements = {zero_sign_vector(length)}
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


def non_negative_covectors_from_matrix(matrix: Matrix) -> set[SignVector]:
    r"""
    Compute all nonnegative covectors.

    OUTPUT:

    Return a set of nonnegative covectors determined by the kernel of ``matrix``.

    EXAMPLES::

        sage: M = matrix([[1, 0, 2, 0], [0, 1, -1, 0], [0, 0, 0, 1]])
        sage: from sign_vectors.oriented_matroids import OrientedMatroid
        sage: from sign_crn.utility import non_negative_covectors_from_matrix
        sage: OrientedMatroid(M).covectors()
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
         (---+),
         (+0+-),
         (0+-+),
         (+0++),
         (++++),
         (----),
         (0+--),
         (+++-)}
        sage: non_negative_covectors_from_matrix(M)
        {(0000), (++00), (++0+), (+0+0), (+++0), (000+), (+0++), (++++)}
    """
    return non_negative_covectors_from_cocircuits(OrientedMatroid(matrix).cocircuits(), matrix.ncols())


def non_negative_vectors_from_matrix(matrix: Matrix) -> set[SignVector]:
    r"""
    Compute all nonnegative covectors from a matrix.

    OUTPUT:

    Return a set of nonnegative covectors determined by the kernel of ``matrix``.

    EXAMPLES::

        sage: M = matrix([[2, -1, -1, 0]])
        sage: from sign_vectors.oriented_matroids import OrientedMatroid
        sage: from sign_crn.utility import non_negative_vectors_from_matrix
        sage: OrientedMatroid(M).vectors()
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
        sage: non_negative_vectors_from_matrix(M)
        {(0000), (++00), (++0+), (+0+0), (+++0), (000+), (+0++), (++++)}
    """
    return non_negative_covectors_from_cocircuits(OrientedMatroid(matrix).circuits(), matrix.ncols())


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

        sage: from sign_crn.utility import closure_minors_utility
        sage: var("a, b, c")
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
            if not is_constant(minor) and not minor in non_zero_expressions:
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
    if is_constant(expression):
        return ZZ(sign(expression))
    return expression


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

    EXAMPLES::

        sage: from sign_crn.utility import equal_entries_lists
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


def vector_from_sign_vector(data, sv: SignVector) -> vector:
    r"""
    Find a vector in the row space of a matrix that has given signs.

    INPUT:

    - ``data`` -- either a real matrix with ``n`` columns or a list of
                elementary vectors of length ``n``
    - ``sv`` -- a sign vector of length ``n``

    OUTPUT:
    Return a conformal sum of elementary vectors that lies in the given subspace.

    If ``data`` is a matrix, the elementary vectors in the kernel of this matrix are used for the result.
    If ``data`` is a list of elementary vectors, those are used.

    .. NOTE::

        A ``ValueError`` is raised if no solution exists.

    EXAMPLES::

        sage: from sign_crn.utility import vector_from_sign_vector
        sage: from elementary_vectors import *
        sage: from sign_vectors import *
        sage: M = matrix([[1, 0, 2, 0], [0, 1, 1, 0], [0, 0, 0, 1]])
        sage: vector_from_sign_vector(M, zero_sign_vector(4))
        (0, 0, 0, 0)
        sage: vector_from_sign_vector(M, sign_vector("+-+0"))
        (2, -2, 2, 0)
        sage: vector_from_sign_vector(M, sign_vector("+0+0"))
        (1, 0, 2, 0)
        sage: vector_from_sign_vector(M, sign_vector("+-0+"))
        (1, -2, 0, 1)
        sage: vector_from_sign_vector(cocircuits(M), sign_vector("+-0+"))
        (1, -2, 0, 1)
        sage: vector_from_sign_vector(M, sign_vector("+0-0"))
        Traceback (most recent call last):
        ...
        ValueError: Cannot find vector corresponding to given sign vector.
        sage: vector_from_sign_vector([], zero_sign_vector(4))
        (0, 0, 0, 0)
    """
    if isinstance(data, list):
        evs = data
        try:
            result = data[0].parent().zero_vector()
        except IndexError:
            result = zero_vector(sv.length())
    elif isinstance(data, Generator):
        evs = data
        result = zero_vector(sv.length())
    else:
        evs_object = CircuitEnumerator(data)
        # evs_object.set_combinations_dual(Combinations(upper.support(), evs_object.length - evs_object.rank + 1))
        evs = evs_object.cocircuit_generator()
        result = zero_vector(data.base_ring(), sv.length())

    if sign_vector(result) == sv:
        return result
    for v in evs:
        for w in [v, -v]:
            if sign_vector(w) <= sv:
                result += w
                if sign_vector(result) == sv:
                    return result
                break

    raise ValueError("Cannot find vector corresponding to given sign vector.")


def intervals_to_sign_vectors(intervals: Intervals) -> Iterator[SignVector]:
    r"""
    Generate all sign vectors that correspond to a vector with components in given intervals.

    INPUT:

    - ``intervals`` -- an `Intervals` object

    EXAMPLES::

        sage: from certlin import *
        sage: from sign_crn.utility import intervals_to_sign_vectors
        sage: intervals = Intervals.from_bounds([-1, 1], [0, 1])
        sage: list(intervals_to_sign_vectors(intervals))
        [(0+), (-+)]
        sage: intervals = Intervals.from_bounds([-1, -2], [0, 1])
        sage: list(intervals_to_sign_vectors(intervals))
        [(00), (0+), (0-), (-0), (-+), (--)]
        sage: intervals = Intervals.from_bounds([-1, -1, 0], [0, 5, 0])
        sage: list(intervals_to_sign_vectors(intervals))
        [(000), (0+0), (0-0), (-00), (-+0), (--0)]
        sage: intervals = Intervals.from_bounds([-1, -1, -1], [0, 1, 0], False, False)
        sage: list(intervals_to_sign_vectors(intervals))
        [(-0-), (-+-), (---)]

    TESTS::

        sage: intervals = Intervals.from_bounds([-1, 0], [1, 0], False, False)
        sage: list(intervals_to_sign_vectors(intervals))
        []
        sage: intervals = Intervals.from_bounds([], [])
        sage: list(intervals_to_sign_vectors(intervals))
        []
    """
    list_of_signs = []
    if intervals.is_empty():
        def empty():
            yield from ()
        return empty()
    for interval in intervals:
        available_signs = []
        if 0 in interval:
            available_signs.append(0)
        if interval.supremum() > 0:
            available_signs.append(1)
        if interval.infimum() < 0:
            available_signs.append(-1)
        list_of_signs.append(available_signs)

    return (
        sign_vector(signs) for signs in cartesian_product_iterator(list_of_signs)
    )


def sign_vector_to_intervals(sv: SignVector) -> Intervals:
    r"""
    Return intervals that correspond to a sign vector.

    EXAMPLES::

        sage: from certlin import *
        sage: from sign_vectors import *
        sage: from sign_crn.utility import sign_vector_to_intervals
        sage: sv = sign_vector("+0-")
        sage: sign_vector_to_intervals(sv)
        (0, +oo) x {0} x (-oo, 0)
    """
    return Intervals([
        Interval(
            0 if element > 0 else (minus_infinity if element < 0 else 0),
            Infinity if element > 0 else (0 if element < 0 else 0),
            element == 0,
            element == 0
        )
        for element in sv
    ])
