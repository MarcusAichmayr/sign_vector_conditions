r"""
Robustness of existence and uniqueness of equilibria
====================================================

Let us consider the following matrices::

    sage: S = matrix([[1, 0, 1, 0], [0, 0, 0, 1]])
    sage: S
    [1 0 1 0]
    [0 0 0 1]
    sage: St = matrix([[1, 0, 1, 1], [0, 1, 0, -1]])
    sage: St
    [ 1  0  1  1]
    [ 0  1  0 -1]


To check, whether the corresponding chemical reaction network
has a unique equilibrium for all rate constants and all small perturbations of ``St``,
we consider the topes of the corresponding oriented matroids::

    sage: from sign_vectors.oriented_matroids import *
    sage: OrientedMatroid(S).topes()
    {(-0-+), (-0--), (+0++), (+0+-)}
    sage: OrientedMatroid(St).topes()
    {(---+), (-+--), (++++), (----), (+-++), (+++-)}

One can see that for every tope ``X`` of the oriented matroid corresponding to ``S`` there is a
tope ``Y`` corresponding to ``St`` such that ``X`` conforms to ``Y``.
Therefore, the exponential map is a diffeomorphism for all ``c > 0``
and all small perturbations of ``St``.
The package offers a function that checks this condition directly::

    sage: from sign_vector_conditions import *
    sage: condition_closure_sign_vectors(S, St)
    True

There is an equivalent condition.
To verify it, we compute the maximal minors of the two matrices::

    sage: S.minors(2)
    [0, 0, 1, 0, 0, 1]
    sage: St.minors(2)
    [1, 0, -1, -1, -1, -1]

From the output, we see whenever a minor of ``S`` is nonzero,
the corresponding minor of ``St`` has the same sign.
Hence, this condition is fulfilled.
This condition can also be checked directly with the package::

    sage: condition_closure_minors(S, St)
    True

Now, we consider matrices with variables::

    sage: var('a, b, c')
    (a, b, c)
    sage: S = matrix([[c, 1, c]])
    sage: S
    [c 1 c]
    sage: St = matrix([[a, b, -1]])
    sage: St
    [ a  b -1]

We cannot check the first condition since there are variables in ``S`` and ``St``.
Therefore, we want to obtain equations on the variables ``a``, ``b``, ``c``
such that this condition is satisfied.
First, we compute the minors of the matrices::

    sage: S.minors(1)
    [c, 1, c]
    sage: St.minors(1)
    [a, b, -1]

The function from the package supports symbolic matrices as input.
In this case, we obtain the following equations on the variables::

    sage: condition_closure_minors(S, St) # random
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

    sage: solve_ineq(list(condition_closure_minors(S, St)[3])) # random
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
from sage.matrix.constructor import Matrix

from sign_vectors.oriented_matroids import OrientedMatroid
from .utility import closure_minors_utility
from elementary_vectors.utility import is_symbolic


def condition_closure_sign_vectors(stoichiometric_matrix: Matrix, kinetic_order_matrix: Matrix) -> bool:
    r"""
    Closure condition for robustness using sign vectors.

    OUTPUT:
    Return whether the closure condition for robustness regarding small perturbations is satisfied.

    .. NOTE::

        This implementation is inefficient and should not be used for large examples.
        Instead, use :func:`~condition_closure_minors`.
    """
    topes = OrientedMatroid(kinetic_order_matrix).topes()
    for covector1 in OrientedMatroid(stoichiometric_matrix).topes():
        if not any(covector1 <= covector2 for covector2 in topes):
            return False
    return True


def condition_closure_minors(stoichiometric_matrix: Matrix, kinetic_order_matrix: Matrix):
    r"""
    Closure condition for robustness using maximal maximal minors.

    OUTPUT:
    Return whether the closure condition for robustness regarding small perturbations is satisfied.
    If the result depends on variables, a list of sets is returned.
    The condition holds if the inequalities in (at least) one of these sets are satisfied.

    .. NOTE::

        The matrices need to have maximal rank and the same dimensions.
        Otherwise, a ``ValueError`` is raised.
    """
    positive_found = False
    negative_found = False
    symbolic_pairs = set()
    for indices in Combinations(stoichiometric_matrix.ncols(), stoichiometric_matrix.nrows()):
        minor1 = stoichiometric_matrix.matrix_from_columns(indices).det()
        if not minor1:
            continue
        minor2 = kinetic_order_matrix.matrix_from_columns(indices).det()
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
