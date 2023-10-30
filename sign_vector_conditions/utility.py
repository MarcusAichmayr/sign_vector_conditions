r"""Utility functions."""

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
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR

from sign_vectors import sign_vector, zero_sign_vector
from sign_vectors.oriented_matroids import cocircuits_from_matrix


def non_negative_cocircuits_from_matrix(M, kernel=True):
    r"""
    Compute non-negative cocircuits.

    INPUT:

    - ``M`` -- a matrix with real arguments.

    - ``kernel`` -- a boolean (default: ``True``)

    OUTPUT:

    Return a set of non-negative cocircuits determined by the kernel of ``M``. (default)
    If ``kernel`` is false, considers the row space of ``M``.

    EXAMPLES::

        sage: M = matrix([[2, -1, -1, 0]])
        sage: M
        [ 2 -1 -1  0]
        sage: from sign_vectors.oriented_matroids import cocircuits_from_matrix
        sage: cocircuits_from_matrix(M)
        {(--00), (000-), (0-+0), (+0+0), (++00), (-0-0), (000+), (0+-0)}
        sage: from sign_vector_conditions.utility import non_negative_cocircuits_from_matrix
        sage: non_negative_cocircuits_from_matrix(M)
        {(+0+0), (++00), (000+)}
    """
    return set(X for X in cocircuits_from_matrix(M, kernel=kernel) if X > 0)


def non_negative_covectors_from_matrix(M, kernel=True):
    r"""
    Compute all non-negative covectors.

    INPUT:

    - ``M`` -- a matrix

    - ``kernel`` -- a boolean (default: ``True``)

    OUTPUT:

    Return a set of non-negative covectors determined by the kernel of ``M``. (default)
    If ``kernel`` is false, considers the row space of ``M``.

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
        sage: from sign_vector_conditions.utility import non_negative_covectors_from_matrix
        sage: non_negative_covectors_from_matrix(M)
        {(+0+0), (+0++), (++00), (+++0), (000+), (++0+), (++++)}
    """
    cocircuits = [cocircuit for cocircuit in cocircuits_from_matrix(M, kernel=kernel) if not cocircuit < 0]

    if not cocircuits:
        raise ValueError('List of cocircuits is empty.')
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


def condition_closure_minors_utility(pairs, positive_only=False, negative_only=False):
    r"""
    Return whether all products of components are positive (or negative) if first element is non-zero.

    INPUT:

    - ``pairs`` -- an iterable of pairs consisting of a minor and a product

    - ``positive_only`` -- a boolean, considers only positive products if true

    - ``negative_only`` -- a boolean, considers only negative products if true

    OUTPUT:
    Returns either a boolean or sets of conditions on variables occurring in the input.
    If the conditions of one of these sets are satisfied,
    then for all non-zero elements of the first list,
    the product with the corresponding element of the second list is positive.
    (Or all products are negative.)

    TESTS::

        sage: from sign_vector_conditions.utility import condition_closure_minors_utility
        sage: var('a, b, c')
        (a, b, c)
        sage: condition_closure_minors_utility(zip([0, a], [0, a]), positive_only=True)
        [{a == 0}, {a > 0}]
        sage: len(_) # for testing
        2
        sage: condition_closure_minors_utility(zip([c, -1, c], [c, -b, -a * c])) # random
        [{-b > 0, c == 0},
         {-b < 0, c == 0},
         {-b > 0, c > 0, -a*c > 0},
         {-b < 0, c < 0, -a*c < 0}]
        sage: len(_) # for testing
        4
        sage: condition_closure_minors_utility(zip([c, -1, a], [c, -b, -a * c])) # random
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
        sage: condition_closure_minors_utility(zip([-1, -1], [-1, -1]))
        True
        sage: condition_closure_minors_utility(zip([-1, 1], [-1, 1]))
        False
        sage: condition_closure_minors_utility(zip([0, 1], [0, 1]))
        True
        sage: condition_closure_minors_utility([(1, 0)])
        False
    """
    def rec(pairs, zero_expressions, non_zero_expressions):
        r"""Recursive call"""
        pairs = [
            (minor, product) for minor, product in pairs
            if not minor in zero_expressions and not minor.is_zero()
        ]
        for minor, _ in pairs:
            if is_symbolic(minor) and not minor in non_zero_expressions:
                yield from rec(pairs, zero_expressions.union([minor]), non_zero_expressions)
                yield from rec(pairs, zero_expressions, non_zero_expressions.union([minor]))

        products = set(
            sign_or_symbolic(product.substitute([value == 0 for value in zero_expressions]))
            for _, product in pairs
        )
        equalities = set(value == 0 for value in zero_expressions)
        non_equalities = set(value != 0 for value in non_zero_expressions if not value in products)

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

    output = list(rec(pairs, set(), set()))
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
    r"""
    Return whether this element is a symbolic expression.

    EXAMPLES::

        sage: from sign_vector_conditions.utility import is_symbolic
        sage: is_symbolic(5)
        False
        sage: var('a, b')
        (a, b)
        sage: is_symbolic(a)
        True
        sage: is_symbolic(-a)
        True
        sage: is_symbolic(b^2 - a)
        True
        sage: is_symbolic(SR(5))
        False
    """
    try:
        return bool(value.variables())
    except AttributeError:
        return False


def sign_or_symbolic(expression):
    r"""Return the sign of an expression if defined."""
    try:
        return ZZ(sign(expression))
    except TypeError:
        return expression


def remove_duplicates(iterable):
    r"""Remove duplicates from a list of iterables."""
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


def non_negative_vectors(vectors):
    r"""
    Return non-negative vectors.

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
        UserWarning: Cannot determine sign of symbolic expression, using 0 for sign vector instead.
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


# TODO: should assumptions be an optional argument?
# def positive_elementary_vectors(data, dim=None, kernel=True, return_minors=False, ring=None, **kwargs):
#     r"""
#     Compute positive elementary vectors.

#     INPUT:

#     - ``data`` -- a matrix or a list of maximal minors

#     - ``dim`` -- Not needed if ``data`` is a matrix.
#                  Otherwise, this is a list of dimensions of the matrix
#                  corresponding to the list of maximal minors ``data``.

#     - ``kernel`` -- a boolean (default: ``True``)

#     - ``return_minors`` -- a boolean (default: ``False``)

#     - ``ring`` -- Parent of the entries of the elementary vectors.
#                   By default, determine this from ``data``.

#     .. NOTE::

#         Keyword arguments may be specified to apply certain reductions to the output.
#         By default, all those reductions (like canceling factors) are applied.
#         Possible keyword arguments are the same as in the function
#         :func:`elementary_vectors.reductions.reduce_vectors`.

#     OUTPUT:

#     The output is a list of pairs.
#     Each pair consists of a list of assumptions and the corresponding positive elementary vectors of ``data``.

#     - If ``kernel`` is true, returns a list of elementary vectors lying in
#       the kernel of ``data``. (default)
#       Otherwise, returns a list of elementary vectors lying in
#       the row space of ``data``.
#       This argument is ignored if ``data`` is not a matrix.

#     - If ``return_minors`` is true, a list is returned where the first
#       element is the list of elementary vectors and the second element is
#       the list of maximal minors used to compute the elementary vectors.

#     EXAMPLES::

#         sage: from elementary_vectors import positive_elementary_vectors
#         sage: A = matrix([[1, -1, 0]])
#         sage: positive_elementary_vectors(A)
#         [[[], [(1, 1, 0), (0, 0, 1)]]]
#         sage: positive_elementary_vectors(A, return_minors=True)
#         [[[], [(1, 1, 0), (0, 0, 1)], [1, -1, 0]]]
#         sage: M = matrix([[0, 0, 1, -1, 0], [2, 0, 0, 0, 2], [1, 1, 1, 1, 1]])
#         sage: positive_elementary_vectors(M)
#         [[[], []]]

#         TODO: Do more examples also symbolic examples
#     """
#     kwargs["cancel_factors"] = False

#     def rec(i, l, eq):
#         if i < len(m):
#             if not sign_determined(m[i]):
#                 a = SR(m[i])
#                 try:
#                     expr = a > 0
#                     assume(expr)
#                     rec(i + 1, l + [expr], eq)
#                     forget(expr)
#                 except ValueError:
#                     pass
#                 try:
#                     expr = a < 0
#                     assume(expr)
#                     rec(i + 1, l + [expr], eq)
#                     forget(expr)
#                 except ValueError:
#                     pass
#                 try:
#                     expr = a == 0
#                     assume(expr)
#                     rec(i + 1, l + [expr], eq + [expr, -a == 0])  # a.substitute(-a == 0) results in a instead of 0.
#                     forget(expr)
#                 except ValueError:
#                     pass
#             else:
#                 rec(i+1, l, eq)
#         else:
#             m_eq = reduce_vector_using_equalities(vector(m), eq=eq)
#             if m_eq == 0:  # minors are all zero
#                 warnings.warn('For ' + str(l) + ' all maximal minors are zero. There might be missing positive elementary vectors.')
# #                 if isinstance(data, list):
# #                     result.append([l, 'maximal minors are zero, could not compute elementary vectors for this branch'])
# #                     return
# #                 else:
# #                     print('TODO')
# #                     # substitute in data
# #                     # compute positive_elementary_vectors from resulting matrix, eventuell auch l übergeben.
# #                     # append each pair of result to result, also append assumptions l to first arguments
# #                     pass
#             # Do not overwrite ``evs`` here! It might be used in another branch.

#             L = reduce_vectors(evs, eq=eq, **kwargs)  # TODO: this might result in an empty list. Instead, compute elementary vectors again if data is matrix
#             if return_minors:
#                 result.append([l, non_negative_vectors(L), list(m_eq)])
#             else:
#                 result.append([l, non_negative_vectors(L)])
#             return

    # def sign_determined(expression):
    #     r"""
    #     Check whether the sign of a number or symbolic expression ``expression`` is uniquely determined.

    #     EXAMPLES::

    #         sage: from elementary_vectors.utility import sign_determined

    #     Integers have always a unique sign::

    #         sage: sign_determined(2)
    #         True
    #         sage: sign_determined(-5)
    #         True

    #     Now, we consider a variable::

    #         sage: var('a')
    #         a
    #         sage: sign_determined(a)
    #         False
    #         sage: assume(a >= 0)
    #         sage: sign_determined(a)
    #         False
    #         sage: assume(a != 0)
    #         sage: sign_determined(a)
    #         True
    #         sage: sign_determined(a - 1)
    #         False
    #     """
    #     return bool(SR(expression) > 0 or SR(expression) < 0 or SR(expression) == 0)

#     evs, m = elementary_vectors(data, dim=dim, kernel=kernel, return_minors=True, ring=ring, **kwargs)
#     if not evs:
#         return []
#     result = []
#     rec(0, [], [])
#     return result
