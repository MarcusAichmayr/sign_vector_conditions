r"""
Uniqueness of equilibria
========================

We define some matrices::

    sage: W = matrix([[1, 0, -1], [0, 1, -1]])
    sage: W
    [ 1  0 -1]
    [ 0  1 -1]
    sage: Wt = matrix([[1, 0, -1], [0, 1, 0]])
    sage: Wt
    [ 1  0 -1]
    [ 0  1  0]

We want to check whether the corresponding chemical reaction network
has at most one equilibrium for all rate constants.
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

Therefore, there is at most one equilibrium.
We can also check this condition in the following way::

    sage: from sign_vector_conditions import *
    sage: condition_uniqueness_sign_vectors(W, Wt)
    True

There is another way to check this condition
that involves the computation of maximal minors of the corresponding matrices::

    sage: m1 = W.minors(2)
    sage: m1
    [1, -1, 1]
    sage: m2 = Wt.minors(2)
    sage: m2
    [1, 0, 1]

We multiply those minors component-wise::

    sage: [m1[i] * m2[i] for i in range(len(m1))]
    [1, 0, 1]

Since all arguments are greater or equal zero, there is at most one equilibrium.
We can also check this condition by applying the following function from this package::

    sage: condition_uniqueness_minors(W, Wt)
    True

Now, we consider another example::

    sage: W = matrix([[1, 0, -1], [0, 1, -1]])
    sage: W
    [ 1  0 -1]
    [ 0  1 -1]
    sage: Wt = matrix([[1, 0, -1], [0, 1, 1]])
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

    sage: condition_uniqueness_sign_vectors(W, Wt)
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
Hence, this condition also states that there is no unique equilibrium::

    sage: condition_uniqueness_minors(W, Wt)
    False

Finally, we consider an example with variables::

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

Therefore, there is at most one equilibrium if and only if :math:`a, b \leq 0`.
The function :func:`~condition_uniqueness_minors` also works for matrices with symbolic entries.
In this case, it returns a system of inequalities::

    sage: condition_uniqueness_minors(W, Wt)
    [{-a >= 0, -b >= 0}]
"""

#############################################################################
#  Copyright (C) 2024                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.combinat.combination import Combinations

from sign_vectors.oriented_matroids import covectors_from_matrix
from elementary_vectors.utility import is_symbolic


def condition_uniqueness_sign_vectors(
    stoichiometric_kernel_matrix, kinetic_order_kernel_matrix
) -> bool:
    r"""
    Uniqueness condition for existence of an equilibrium using sign vectors.

    INPUT:

    - ``stoichiometric_kernel_matrix`` -- a matrix with ``n`` columns

    - ``kinetic_order_kernel_matrix`` -- a matrix with ``n`` columns

    OUTPUT:
    Return whether there exists at most one equilibrium.

    .. NOTE::

        This implementation is inefficient and should not be used for large examples.
        Instead, use :func:`~condition_uniqueness_minors`.

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
        sage: condition_uniqueness_sign_vectors(W, Wt)
        True
        sage: W = matrix([[1, 0, -1], [0, 1, -1]])
        sage: W
        [ 1  0 -1]
        [ 0  1 -1]
        sage: Wt = matrix([[1, 0, -1], [0, 1, 1]])
        sage: Wt
        [ 1  0 -1]
        [ 0  1  1]
        sage: condition_uniqueness_sign_vectors(W, Wt)
        False

    TESTS::

        sage: from sign_vector_conditions.uniqueness import condition_uniqueness_sign_vectors
        sage: A = identity_matrix(3)
        sage: B = A # kernel of B is empty
        sage: condition_uniqueness_sign_vectors(A, B)
        True
    """
    if stoichiometric_kernel_matrix.ncols() != kinetic_order_kernel_matrix.ncols():
        raise ValueError("Matrices have different number of columns.")

    return (
        len(
            covectors_from_matrix(
                stoichiometric_kernel_matrix, kernel=False
            ).intersection(
                covectors_from_matrix(kinetic_order_kernel_matrix, kernel=True)
            )
        )
        == 1
    )


def condition_uniqueness_minors(
    stoichiometric_kernel_matrix, kinetic_order_kernel_matrix
):
    r"""
    Uniqueness condition for existence of an equilibrium using maximal minors.

    INPUT:

    - ``stoichiometric_kernel_matrix`` -- a matrix

    - ``kinetic_order_kernel_matrix`` -- a matrix

    OUTPUT:
    Return whether there exists at most one equilibrium.
    If the result depends on variables, a list of sets is returned.
    The condition holds if the inequalities in exactly one of these sets are satisfied.

    .. NOTE::

        The matrices need to have the same rank and number of columns.
        Otherwise, a ``ValueError`` is raised.

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
        [{-a >= 0, -b >= 0}]
        sage: W = matrix([
        ....:     [a, 0, 1, 0],
        ....:     [0, 1, -1, 0],
        ....:     [0, 0, 0, 1]
        ....: ])
        sage: Wt = matrix([
        ....:     [1, 0, 0, -1],
        ....:     [0, b, 1, 1],
        ....:     [0, 0, a, 1]
        ....: ])
        sage: condition_uniqueness_minors(W, Wt) # random
        [{(a - 1)*a >= 0, a*b >= 0}, {(a - 1)*a <= 0, a*b <= 0}]
        sage: len(_), len(_[0]) # for testing
        (2, 2)

    We can also apply the built-in function ``solve_ineq`` to the resulting sets of inequalities.
    For instance, the first set can be equivalently written as::

        sage: solve_ineq(list(condition_uniqueness_minors(W, Wt)[0])) # random
        [[b == 0, a == 0],
        [a == 0],
        [b == 0, a == 1],
        [a == 1, 0 < b],
        [b == 0, 1 < a],
        [0 < b, 1 < a],
        [b == 0, a < 0],
        [b < 0, a < 0]]
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

    positive_product_found = False
    negative_product_found = False
    symbolic_expressions = set()

    for indices in Combinations(
        stoichiometric_kernel_matrix.ncols(), stoichiometric_kernel_matrix.nrows()
    ):
        minor1 = stoichiometric_kernel_matrix.matrix_from_columns(indices).det()
        if not minor1:
            continue
        product = (
            minor1 * kinetic_order_kernel_matrix.matrix_from_columns(indices).det()
        )
        if is_symbolic(product):
            symbolic_expressions.add(product)
        elif product > 0:
            positive_product_found = True
        elif product < 0:
            negative_product_found = True
        if positive_product_found and negative_product_found:
            return False
    if positive_product_found:
        if symbolic_expressions:
            return [set(expression >= 0 for expression in symbolic_expressions)]
        return True
    if negative_product_found:
        if symbolic_expressions:
            return [set(expression <= 0 for expression in symbolic_expressions)]
        return True
    if symbolic_expressions:
        return [
            set(expression >= 0 for expression in symbolic_expressions),
            set(expression <= 0 for expression in symbolic_expressions),
        ]
    return False
