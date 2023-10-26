r"""
Robustness of existence and uniqueness of equilibria.

EXAMPLES::

    sage: from sign_vector_conditions import *

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

from sign_vectors.oriented_matroids import topes_from_matrix
from .utility import condition_on_products


def condition_closure_sign_vectors(W, Wt):
    r"""
    Closure condition for robustness using sign vectors.

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    OUTPUT:
    Return whether the closure condition for robustness regarding small perturbations is satisfied.

    .. NOTE::

        This implementation is inefficient and should not be used for large examples.
        Instead, use :func:`~condition_closure_minors`.
    """
    topes = topes_from_matrix(Wt, kernel=True)
    for covector1 in topes_from_matrix(W, kernel=True):
        if not any(covector1 <= covector2 for covector2 in topes):
            return False
    return True


def condition_closure_minors(W, Wt):
    r"""
    Closure condition for robustness using maximal maximal minors.

    INPUT:

    - ``W`` -- a matrix

    - ``Wt`` -- a matrix with the same dimensions as ``W``

    OUTPUT:
    Return whether the closure condition for robustness regarding small perturbations is satisfied.
    If the result depends on variables, a list of sets is returned.
    The condition holds if the inequalities in (at least) one of these sets are satisfied.

    .. NOTE::

        The matrices need to have the same rank and number of columns.
        Otherwise, a ``ValueError`` is raised.
    """
    W = W.matrix_from_rows(W.pivot_rows())
    Wt = Wt.matrix_from_rows(Wt.pivot_rows())
    if W.dimensions() != Wt.dimensions():
        raise ValueError('Matrices must have same rank and number of columns.')

    return condition_on_products(W.minors(W.nrows()), Wt.minors(W.nrows()))
