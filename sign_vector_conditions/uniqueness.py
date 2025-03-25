r"""
Uniqueness of equilibria
========================

We define some matrices::

    sage: S = matrix([[1, 1, 1]])
    sage: S
    [1 1 1]
    sage: St = matrix([[1, 0, 1]])
    sage: St
    [1 0 1]

We want to check whether the corresponding chemical reaction network
has at most one equilibrium for all rate constants.
For this purpose, we compute the corresponding oriented matroids::

    sage: from sign_vectors.oriented_matroids import *
    sage: cvS = covectors_from_matrix(S, dual=True, algorithm='fe')
    sage: cvS
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
    sage: cvSt = covectors_from_matrix(St, dual=False, algorithm='fe')
    sage: cvSt
    {(000), (+0+), (-0-)}

The intersection of these oriented matroids consists only of the zero sign vector.
We can compute the intersection directly by applying the built in method intersection::

    sage: set(cvS).intersection(cvSt)
    {(000)}

Therefore, there is at most one equilibrium.
We can also check this condition in the following way::

    sage: from sign_vector_conditions import *
    sage: condition_uniqueness_sign_vectors(S, St)
    True

There is another way to check this condition
that involves the computation of maximal minors of the corresponding matrices::

    sage: m1 = S.minors(1)
    sage: m1
    [1, 1, 1]
    sage: m2 = St.minors(1)
    sage: m2
    [1, 0, 1]

We multiply those minors component-wise::

    sage: [m1[i] * m2[i] for i in range(len(m1))]
    [1, 0, 1]

Since all arguments are greater or equal zero, there is at most one equilibrium.
We can also check this condition by applying the following function from this package::

    sage: condition_uniqueness_minors(S, St)
    True

Now, we consider another example::

    sage: S = matrix([[1, 1, 1]])
    sage: S
    [1 1 1]
    sage: St = matrix([[1, 0, -1], [0, 1, 1]])
    sage: St = matrix([[1, -1, 1]])
    sage: St
    [ 1 -1  1]

Next, we compute the corresponding oriented matroids::

    sage: covectors_from_matrix(S, dual=True, algorithm='fe', separate=True)
    [{(000)},
     {(-+0), (0-+), (0+-), (+-0), (+0-), (-0+)},
     {(+-+), (-+-), (--+), (-++), (++-), (+--)}]
    sage: covectors_from_matrix(St, dual=False, algorithm='fe', separate=True)
    [{(000)}, {(+-+), (-+-)}]

Now, we check the condition from before::

    sage: condition_uniqueness_sign_vectors(S, St)
    False

Therefore, the corresponding exponential map is not injective.
Furthermore, we obtain the following minors::

    sage: m1 = S.minors(1)
    sage: m1
    [1, 1, 1]
    sage: m2 = St.minors(1)
    sage: m2
    [1, -1, 1]
    sage: [m1[i] * m2[i] for i in range(len(m1))]
    [1, -1, 1]

There are positive and negative elements in the resulting list.
Hence, this condition also states that there is no unique equilibrium::

    sage: condition_uniqueness_minors(S, St)
    False

Finally, we consider an example with variables::

    sage: var('a, b')
    (a, b)
    sage: S = matrix([[1, 1, 1]])
    sage: S
    [1 1 1]
    sage: St = matrix([[a, b, -1]])
    sage: St
    [ a  b -1]

The matrix ``St`` contains variables :math:`a, b \in \mathbb{R}`.
Consequently, we cannot compute the corresponding oriented matroids.
On the other hand, we can still compute the minors of ``S`` and ``St``, that is::

    sage: m1 = S.minors(1)
    sage: m1
    [1, 1, 1]
    sage: m2 = St.minors(1)
    sage: m2
    [a, b, -1]
    sage: [m1[i] * m2[i] for i in range(len(m1))]
    [a, b, -1]

Therefore, there is at most one equilibrium if and only if :math:`a, b \leq 0`.
The function :func:`~condition_uniqueness_minors` also works for matrices with symbolic entries.
In this case, it returns a system of inequalities::

    sage: condition_uniqueness_minors(S, St) # random order
    [{-a >= 0, -b >= 0}]
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

from sign_vectors.oriented_matroids import covectors_from_matrix
from elementary_vectors.utility import is_symbolic


def condition_uniqueness_sign_vectors(stoichiometric_matrix, kinetic_order_matrix) -> bool:
    r"""
    Uniqueness condition for existence of an equilibrium using sign vectors.

    INPUT:

    - ``stoichiometric_matrix`` -- a matrix

    - ``kinetic_order_matrix`` -- a matrix

    OUTPUT:
    Return whether there exists at most one equilibrium.

    .. NOTE::

        This implementation is inefficient and should not be used for large examples.
        Instead, use :func:`~condition_uniqueness_minors`.

    EXAMPLES::

        sage: from sign_vector_conditions import *
        sage: S = matrix([[1, 1, 1]])
        sage: S
        [1 1 1]
        sage: St = matrix([[1, 0, 1]])
        sage: St
        [1 0 1]
        sage: condition_uniqueness_sign_vectors(S, St)
        True
        sage: S = matrix([[1, 0, -1], [0, 1, -1]])
        sage: S
        [ 1  0 -1]
        [ 0  1 -1]
        sage: St = matrix([[1, 0, -1], [0, 1, 1]])
        sage: St
        [ 1  0 -1]
        [ 0  1  1]
        sage: condition_uniqueness_sign_vectors(S, St)
        False

    TESTS::

        sage: from sign_vector_conditions.uniqueness import condition_uniqueness_sign_vectors
        sage: A = identity_matrix(3)
        sage: B = A # kernel of B is empty
        sage: condition_uniqueness_sign_vectors(A, B)
        True
    """
    if stoichiometric_matrix.ncols() != kinetic_order_matrix.ncols():
        raise ValueError("Matrices have different number of columns.")

    covectors = covectors_from_matrix(stoichiometric_matrix, dual=True)

    counter = 0
    for covector in covectors_from_matrix(kinetic_order_matrix, dual=False):
        if covector in covectors:
            counter += 1
            if counter > 1:
                return False
    return True


def condition_uniqueness_minors(stoichiometric_matrix, kinetic_order_matrix):
    r"""
    Uniqueness condition for existence of an equilibrium using maximal minors.

    INPUT:

    - ``stoichiometric_matrix`` -- a matrix

    - ``kinetic_order_matrix`` -- a matrix

    OUTPUT:
    Return whether there exists at most one equilibrium.
    If the result depends on variables, a list of sets is returned.
    The condition holds if the inequalities in exactly one of these sets are satisfied.

    .. NOTE::

        The matrices need to have the same rank and number of columns.
        Otherwise, a ``ValueError`` is raised.

    EXAMPLES::

        sage: from sign_vector_conditions import *
        sage: S = matrix([[1, 0, -1], [0, 1, -1]])
        sage: S
        [ 1  0 -1]
        [ 0  1 -1]
        sage: St = matrix([[1, 0, -1], [0, 1, 0]])
        sage: St
        [ 1  0 -1]
        [ 0  1  0]
        sage: condition_uniqueness_minors(S, St)
        True
        sage: S = matrix([[1, 0, -1], [0, 1, -1]])
        sage: S
        [ 1  0 -1]
        [ 0  1 -1]
        sage: St = matrix([[1, 0, -1], [0, 1, 1]])
        sage: St
        [ 1  0 -1]
        [ 0  1  1]
        sage: condition_uniqueness_minors(S, St)
        False
        sage: var('a, b')
        (a, b)
        sage: S = matrix([[1, 0, -1], [0, 1, -1]])
        sage: S
        [ 1  0 -1]
        [ 0  1 -1]
        sage: St = matrix([[1, 0, a], [0, 1, b]])
        sage: St
        [1 0 a]
        [0 1 b]
        sage: condition_uniqueness_minors(S, St) # random order
        [{-a >= 0, -b >= 0}]
        sage: conditions = condition_uniqueness_minors(S, St)[0]
        sage: conditions # random order
        sage: (-a >= 0) in conditions and (-b >= 0) in conditions # 
        True
        sage: S = matrix([[a, 0, 1, 0], [0, 1, -1, 0], [0, 0, 0, 1]])
        sage: St = matrix([[1, 0, 0, -1], [0, b, 1, 1], [0, 0, a, 1]])
        sage: condition_uniqueness_minors(S, St) # random
        [{(a - 1)*a >= 0, a*b >= 0}, {(a - 1)*a <= 0, a*b <= 0}]
        sage: len(_), len(_[0]) # for testing
        (2, 2)

    We can also apply the built-in function ``solve_ineq`` to the resulting sets of inequalities.
    For instance, the first set can be equivalently written as::

        sage: solve_ineq(list(condition_uniqueness_minors(S, St)[0])) # random
        [[b == 0, a == 0],
        [a == 0],
        [b == 0, a == 1],
        [a == 1, 0 < b],
        [b == 0, 1 < a],
        [0 < b, 1 < a],
        [b == 0, a < 0],
        [b < 0, a < 0]]
    """
    stoichiometric_matrix = stoichiometric_matrix.matrix_from_rows(stoichiometric_matrix.pivot_rows())
    kinetic_order_matrix = kinetic_order_matrix.matrix_from_rows(kinetic_order_matrix.pivot_rows())
    if stoichiometric_matrix.dimensions() != kinetic_order_matrix.dimensions():
        raise ValueError("Matrices must have same rank and number of columns.")

    positive_product_found = False
    negative_product_found = False
    symbolic_expressions = set()

    for indices in Combinations(stoichiometric_matrix.ncols(), stoichiometric_matrix.nrows()):
        minor1 = stoichiometric_matrix.matrix_from_columns(indices).det()
        if not minor1:
            continue
        product = (
            minor1 * kinetic_order_matrix.matrix_from_columns(indices).det()
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
