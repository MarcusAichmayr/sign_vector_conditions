r"""
Robustness of existence and uniqueness of equilibria
====================================================

Let us consider the following matrices::

    sage: W = matrix([[1, 0, -1, 0], [0, 1, 0, 0]])
    sage: W
    [ 1  0 -1  0]
    [ 0  1  0  0]
    sage: Wt = matrix([[1, 0, -1, 0], [0, 1, -1, 1]])
    sage: Wt
    [ 1  0 -1  0]
    [ 0  1 -1  1]

To check, whether the corresponding chemical reaction network
has a unique equilibrium for all rate constants and all small perturbations of ``Wt``,
we consider the topes of the corresponding oriented matroids::

    sage: from sign_vectors.oriented_matroids import *
    sage: topes_from_matrix(W, dual=True)
    {(+0+-), (+0++), (-0--), (-0-+)}
    sage: topes_from_matrix(Wt, dual=True)
    {(---+), (-+--), (++++), (----), (+-++), (+++-)}

One can see that for every tope ``X`` of the oriented matroid corresponding to ``W`` there is a
tope ``Y`` corresponding to ``Wt`` such that ``X`` conforms to ``Y``.
Therefore, the exponential map is a diffeomorphism for all ``c > 0``
and all small perturbations of ``Wt``.
The package offers a function that checks this condition directly::

    sage: from sign_vector_conditions import *
    sage: condition_closure_sign_vectors(W, Wt)
    True

There is an equivalent condition.
To verify it, we compute the maximal minors of the two matrices::

    sage: W.minors(2)
    [1, 0, 0, 1, 0, 0]
    sage: Wt.minors(2)
    [1, -1, 1, 1, 0, -1]

From the output, we see whenever a minor of ``W`` is nonzero,
the corresponding minor of ``Wt`` has the same sign.
Hence, this condition is fulfilled.
This condition can also be checked directly with the package::

    sage: condition_closure_minors(W, Wt)
    True

Now, we consider matrices with variables::

    sage: var('a, b, c')
    (a, b, c)
    sage: W = matrix([[1, 0, -1], [0, c, -1]])
    sage: W
    [ 1  0 -1]
    [ 0  c -1]
    sage: Wt = matrix([[1, 0, a], [0, 1, b]])
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
if ``c`` is zero and ``b`` is nonzero.
The closure condition is also satisfied if ``a`` and ``b`` are negative and ``c`` is positive
or if ``a`` and ``b`` are positive and ``c`` is negative.

We can also apply the built-in function ``solve_ineq`` to the resulting sets of inequalities.
For instance, the last set can be equivalently written as::

    sage: solve_ineq(list(condition_closure_minors(W, Wt)[3])) # random
    [[c < 0, 0 < b, a < 0]]
"""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.combinat.combination import Combinations

from sign_vectors.oriented_matroids import topes_from_matrix
from .utility import closure_minors_utility
from elementary_vectors.utility import is_symbolic


def condition_closure_sign_vectors(
    stoichiometric_kernel_matrix, kinetic_order_kernel_matrix
) -> bool:
    r"""
    Closure condition for robustness using sign vectors.

    INPUT:

    - ``stoichiometric_kernel_matrix`` -- a matrix with ``n`` columns

    - ``kinetic_order_kernel_matrix`` -- a matrix with ``n`` columns

    OUTPUT:
    Return whether the closure condition for robustness regarding small perturbations is satisfied.

    .. NOTE::

        This implementation is inefficient and should not be used for large examples.
        Instead, use :func:`~condition_closure_minors`.
    """
    topes = topes_from_matrix(kinetic_order_kernel_matrix, dual=True)
    for covector1 in topes_from_matrix(stoichiometric_kernel_matrix, dual=True):
        if not any(covector1 <= covector2 for covector2 in topes):
            return False
    return True


def condition_closure_minors(stoichiometric_kernel_matrix, kinetic_order_kernel_matrix):
    r"""
    Closure condition for robustness using maximal maximal minors.

    INPUT:

    - ``stoichiometric_kernel_matrix`` -- a matrix

    - ``kinetic_order_kernel_matrix`` -- a matrix with the same dimensions as ``W``

    OUTPUT:
    Return whether the closure condition for robustness regarding small perturbations is satisfied.
    If the result depends on variables, a list of sets is returned.
    The condition holds if the inequalities in (at least) one of these sets are satisfied.

    .. NOTE::

        The matrices need to have the same rank and number of columns.
        Otherwise, a ``ValueError`` is raised.
    """
    stoichiometric_kernel_matrix = stoichiometric_kernel_matrix.matrix_from_rows(
        stoichiometric_kernel_matrix.pivot_rows()
    )
    kinetic_order_kernel_matrix = kinetic_order_kernel_matrix.matrix_from_rows(
        kinetic_order_kernel_matrix.pivot_rows()
    )
    if (
        stoichiometric_kernel_matrix.dimensions()
        != kinetic_order_kernel_matrix.dimensions()
    ):
        raise ValueError("Matrices must have same rank and number of columns.")

    positive_found = False
    negative_found = False
    symbolic_pairs = set()
    for indices in Combinations(
        stoichiometric_kernel_matrix.ncols(), stoichiometric_kernel_matrix.nrows()
    ):
        minor1 = stoichiometric_kernel_matrix.matrix_from_columns(indices).det()
        if not minor1:
            continue
        minor2 = kinetic_order_kernel_matrix.matrix_from_columns(indices).det()
        if not minor2:
            return False
        product = minor1 * minor2
        if is_symbolic(product):
            symbolic_pairs.add((minor1, product))
            continue
        if product > 0:
            positive_found = True
        elif product < 0:
            negative_found = True
        if positive_found and negative_found:
            return False

    return closure_minors_utility(symbolic_pairs, positive_found, negative_found)
