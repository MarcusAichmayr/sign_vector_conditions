r"""
Conditions for CBE
==================

Uniqueness of CBE
-----------------

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
    sage: cvS = OrientedMatroid(S).vectors()
    sage: cvS
    {(000),
     (0+-),
     (+-0),
     (-+0),
     (++-),
     (-0+),
     (--+),
     (-++),
     (+0-),
     (0-+),
     (+--),
     (-+-),
     (+-+)}
    sage: cvSt = OrientedMatroid(St).covectors()
    sage: cvSt
    {(000), (+0+), (-0-)}

The intersection of these oriented matroids consists only of the zero sign vector.
We can compute the intersection directly by applying the built in method intersection::

    sage: set(cvS).intersection(cvSt)
    {(000)}

Therefore, there is at most one equilibrium.
We can also check this condition in the following way::

    sage: from sign_crn import *
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

    sage: OrientedMatroid(S).dual().faces()
    [{(000)},
     {(0+-), (+-0), (-+0), (+0-), (0-+), (-0+)},
     {(++-), (--+), (-++), (+--), (-+-), (+-+)}]
    sage: OrientedMatroid(St).faces()
    [{(000)}, {(-+-), (+-+)}]

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

Existence and uniqueness of CBE
-------------------------------

Let us consider the following matrices to describe a chemical reaction network::

    sage: S = matrix([[1, 0, -1, 0], [0, 1, 0, -1]])
    sage: S
    [ 1  0 -1  0]
    [ 0  1  0 -1]
    sage: St = matrix([[1, 0, -1, 1], [0, 1, -1, 0]])
    sage: St
    [ 1  0 -1  1]
    [ 0  1 -1  0]


To check whether a unique equilibrium exists, we apply :func:`~condition_uniqueness_minors`::

    sage: from sign_crn import *
    sage: condition_uniqueness_minors(S, St)
    True

This means that the chemical reaction network has at most one equilibrium.
Next, we verify whether an equilibrium exists.
First, we check the face condition.
For this purpose, we compute the cocircuits of the oriented matroids
corresponding to the matrices::

    sage: from sign_vectors.oriented_matroids import *
    sage: cc1 = OrientedMatroid(S).circuits()
    sage: cc1
    {(0+0+), (+0+0), (0-0-), (-0-0)}
    sage: cc2 = OrientedMatroid(St).circuits()
    sage: cc2
    {(---0), (-00+), (+00-), (+++0), (0+++), (0---)}

Here, we are only interested in the positive cocircuits::

    sage: cc1p = [X for X in cc1 if X > 0]
    sage: cc1p
    [(0+0+), (+0+0)]
    sage: cc2p = [X for X in cc2 if X > 0]
    sage: cc2p
    [(+++0), (0+++)]

Since every sign vector in ``cc2p`` has a smaller element in ``cc1p``,
the face condition is satisfied.
There is also a function in the package that can be used directly
to check whether this condition is fulfilled::

    sage: condition_faces(S, St)
    True

We need to check a third condition to verify surjectivity.
For this purpose, we consider again the oriented matroid determined by ``S``::

    sage: OrientedMatroid(S).covectors()
    {(0000), (-++-), (+0-0), (0-0+), (+--+), (-0+0), (--++), (0+0-), (++--)}

Since there are no nonnegative covectors, the chemical reaction network has at least one equilibrium.
The package offers a function to check this condition condition::

    sage: condition_nondegenerate(S, St)
    True

Hence, the chemical reaction network has a unique equilibrium.

Let us consider another example.
We swap the two matrices from before::

    sage: S, St = St, S

Because of symmetry, there is at most one equilibrium::

    sage: condition_uniqueness_sign_vectors(S, St)
    True

Now, we check the face condition::

    sage: cc1 = OrientedMatroid(S).circuits()
    sage: cc1
    {(---0), (-00+), (+00-), (+++0), (0+++), (0---)}
    sage: cc2 = OrientedMatroid(St).circuits()
    sage: cc2
    {(0+0+), (+0+0), (0-0-), (-0-0)}

Again, we are only interested in the positive cocircuits::

    sage: cc1p = [X for X in cc1 if X > 0]
    sage: cc1p
    [(+++0), (0+++)]
    sage: cc2p = [X for X in cc2 if X > 0]
    sage: cc2p
    [(0+0+), (+0+0)]

Therefore, the condition does not hold.
We also apply the corresponding function from the package::

    sage: condition_faces(S, St)
    False

Consequently, there exists no unique equilibrium.

Now, we consider a network given by two matrices involving a parameter.
(see Example 20 of [MHR19]_)
Depending on this parameter, the network has a unique positive CBE::

    sage: var("a")
    a
    sage: assume(a > 0)
    sage: S = matrix([[1, 0, 0, 0, 0, 1], [0, 1, 0, 0, 0, -1], [0, 0, 1, 1, 2, 0]])
    sage: S
    [ 1  0  0  0  0  1]
    [ 0  1  0  0  0 -1]
    [ 0  0  1  1  2  0]
    sage: St = matrix([[-1, -1, 0, 0, -2, 0], [0, 0, 1, 1, 0, 0], [0, 0, 0, 0, a, 1]])
    sage: St
    [-1 -1  0  0 -2  0]
    [ 0  0  1  1  0  0]
    [ 0  0  0  0  a  1]

The first two conditions depend on the sign vectors corresponding
to the rows of these matrices which are independent of the specific value for :math:`a`::

    sage: condition_uniqueness_sign_vectors(S, St)
    True

Hence, there exists at most one equilibrium.
Also the face condition is satisfied::

    sage: condition_faces(S, St)
    True

For specific values of :math:`a`, the pair of subspaces
determined by kernels of the matrices is nondegenerate.
This is exactly the case for :math:`a \in (0, 1) \cup (1, 2)`.
We demonstrate this for specific values::

    sage: condition_nondegenerate(S, St(a=1/2))
    True
    sage: condition_nondegenerate(S, St(a=3/2))
    True

On the other hand, this condition does not hold if
:math:`a \in {1} \cup [2, \infty)`::

    sage: condition_nondegenerate(S, St(a=1))
    False

To certify the result, we call::

    sage: condition_degenerate(S, St(a=1), certify=True)
    (True, (1, 1, 0, 0, -1, 1))

Hence, the positive support of the vector ``v = (1, 1, 0, 0, -1, 1)`` of ``St``
can be covered by a sign vector ``(++000+)`` corresponding to ``ker(S)``.
Further, ``v`` does not satisfy the support condition.

    sage: condition_nondegenerate(S, St(a=2))
    False
    sage: condition_nondegenerate(S, St(a=3))
    False

Robustness of existence and uniqueness of CBE
---------------------------------------------

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

    sage: from sign_crn import *
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

from copy import copy

from sage.combinat.combination import Combinations
from sage.matrix.constructor import Matrix
from sage.rings.infinity import Infinity

from sign_vectors import SignVector, OrientedMatroid
from elementary_vectors.utility import is_constant
from certlin import Intervals, LinearInequalitySystem
from .utility import (
    closure_minors_utility,
    intervals_to_sign_vectors,
    sign_vector_to_intervals,
    non_negative_circuits_from_matrix,
    non_negative_cocircuits_from_matrix,
    equal_entries_lists,
    vector_from_sign_vector
)


def condition_uniqueness_sign_vectors(stoichiometric_matrix: Matrix, kinetic_order_matrix: Matrix) -> bool:
    r"""
    Uniqueness condition for existence of an equilibrium using sign vectors.

    OUTPUT:
    Return whether there exists at most one equilibrium.

    .. NOTE::

        This implementation is inefficient and should not be used for large examples.
        Instead, use :func:`~condition_uniqueness_minors`.

    EXAMPLES::

        sage: from sign_crn import *
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

        sage: from sign_crn.uniqueness import condition_uniqueness_sign_vectors
        sage: A = identity_matrix(3)
        sage: B = A # kernel of B is empty
        sage: condition_uniqueness_sign_vectors(A, B)
        True
    """
    covectors = OrientedMatroid(stoichiometric_matrix).covectors()
    counter = 0
    for covector in OrientedMatroid(kinetic_order_matrix).vectors():
        if covector in covectors:
            counter += 1
            if counter > 1:
                return False
    return True


def condition_uniqueness_minors(stoichiometric_matrix: Matrix, kinetic_order_matrix: Matrix):
    r"""
    Uniqueness condition for existence of an equilibrium using maximal minors.

    OUTPUT:
    Return whether there exists at most one equilibrium.
    If the result depends on variables, a list of sets is returned.
    The condition holds if the inequalities in exactly one of these sets are satisfied.

    .. NOTE::

        The matrices need to have maximal rank and the same dimensions.
        Otherwise, a ``ValueError`` is raised.

    EXAMPLES::

        sage: from sign_crn import *
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
        if not is_constant(product):
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


def condition_faces(stoichiometric_matrix: Matrix, kinetic_order_matrix: Matrix) -> bool:
    r"""
    Condition on positive sign vectors for existence and uniqueness of equilibria

    OUTPUT:
    TODO
    Return whether every positive sign vector ``X`` corresponding to the rows of
    ``St`` has a positive sign vector ``Y`` corresponding to the rows of ``S``
    such that ``Y <= X``.

    Return a boolean.

    EXAMPLES::

        sage: from sign_crn.unique_existence import condition_faces
        sage: S = matrix([[1, 0, -1, 0], [0, 1, 0, -1]])
        sage: S
        [ 1  0 -1  0]
        [ 0  1  0 -1]
        sage: St = matrix([[1, 0, -1, 1], [0, 1, -1, 0]])
        sage: St
        [ 1  0 -1  1]
        [ 0  1 -1  0]
        sage: condition_faces(S, St)
        True
    """
    non_negative_cocircuits = non_negative_circuits_from_matrix(stoichiometric_matrix)

    for cocircuit1 in non_negative_circuits_from_matrix(kinetic_order_matrix):
        if not any(cocircuit2 <= cocircuit1 for cocircuit2 in non_negative_cocircuits):
            return False
    return True


def condition_nondegenerate(stoichiometric_matrix: Matrix, kinetic_order_matrix: Matrix) -> bool:
    r"""
    Return whether a pair of subspaces given by matrices is nondegenerate.

    OUTPUT:
    a boolean

    .. SEEALSO::

        :func:`~condition_degenerate`
    """
    return not condition_degenerate(stoichiometric_matrix, kinetic_order_matrix)


def condition_degenerate(stoichiometric_matrix: Matrix, kinetic_order_matrix: Matrix, certify: bool = False) -> bool:
    r"""
    Return whether a pair of subspaces given by matrices is degenerate.

    This condition is about whether all positive equal components of a vector in ``St``
    can be covered by covectors corresponding to the kernel of ``S``.

    OUTPUT:
    a boolean

    If ``certify`` is true, a list is returned to certify the result.
    (see the examples)

    EXAMPLES::

        sage: from sign_crn.unique_existence import *

    Next, we certify our results. In the first examples, the subspaces are trivially nondegenerate
    since there are no nonnegative covectors in the kernel of ``S``::

        sage: S = matrix([[1, 0, -1, 0], [0, 1, 0, -1]])
        sage: St = matrix([[1, 0, 0, 1], [0, 1, 0, 1]])
        sage: condition_degenerate(S, St, certify=True)
        (False, 'no nonnegative covectors')

    Here, we have a pair of degenerate subspaces::

        sage: S = matrix([[1, 1, 0]])
        sage: St = matrix([[0, 0, 1]])
        sage: condition_degenerate(S, St, certify=True)
        (True, (1, 1, 0))

    The resulting vector lies in the row space of ``St``.
    The nonnegative covector ``(++0)`` in the kernel of ``S`` covers the first two equal components.

    In the following, we have another example for nondegenerate subspaces::

        sage: S = matrix([[1, 0, 0, 1, -1], [0, 1, 0, 1, -1], [0, 0, 1, 0, 1]])
        sage: S
        [ 1  0  0  1 -1]
        [ 0  1  0  1 -1]
        [ 0  0  1  0  1]
        sage: St = matrix([[1, 0, 0, 1, -1], [0, 1, 0, 1, -1], [0, 0, 1, 0, -1]])
        sage: St
        [ 1  0  0  1 -1]
        [ 0  1  0  1 -1]
        [ 0  0  1  0 -1]
        sage: condition_degenerate(S, St, certify=True)
        (False, ([[[1, 2, 3]], [[0, 2, 3]]], [[[2, 4]]], []))

    The certificate tells us that there is no vector in the row space of ``St``
    with positive support on the components ``0, 2, 3`` and ``1, 2, 3``.
    Positive equal components can partially be covered by a covector ``(00+0+)``
    which corresponds to ``[[2, 4]]``.
    However, it is impossible to fully cover the positive support.

    In the next example, there exists a partial cover::

        sage: S = matrix([[1, -1, 0, 0], [0, 0, 1, 1]])
        sage: St = matrix([[1, 0, 0, 1], [0, 1, 0, 1]])
        sage: condition_degenerate(S, St, certify=True)
        (False, ([], [[[2, 3]]], [[[[2, 3]], [(--++)]]]))

    In fact, a vector in ``St`` with equal positive components on ``[2, 3]``
    corresponding to ``(--++)`` can be fully covered by covectors.
    However, this vector would not satisfy the support condition.
    """
    if stoichiometric_matrix.ncols() != kinetic_order_matrix.ncols():
        raise ValueError("Matrices have different number of columns.")
    non_negative_cocircuits = non_negative_cocircuits_from_matrix(stoichiometric_matrix)

    if not non_negative_cocircuits:
        if certify:
            return False, "no nonnegative covectors"
        return False

    non_negative_cocircuits = sorted(non_negative_cocircuits, key=lambda covector: len(covector.support()))
    length = kinetic_order_matrix.ncols()
    degenerate = False

    lower_bounds = [-Infinity] * length
    upper_bounds = [0] * length
    upper_bounds_inf = [Infinity] * length

    kernel_matrix = kinetic_order_matrix
    covectors_support_condition = non_negative_circuits_from_matrix(stoichiometric_matrix)

    if certify:
        certificate = []
        certificates_zero_equal_components = []
        certificates_partial_cover = []
        certificate_support_condition = []

    def recursive_degenerate(
        non_negative_cocircuits: set[SignVector],
        matrix_old: Matrix,
        indices: list[int],
        lower_bounds: list[int],
        upper_bounds: list[int]
    ):
        r"""
        Recursive function.

        INPUT:

        - ``non_negative_cocircuits`` -- a list of positive sign vectors
        - ``lower_bounds`` -- a list of values ``-Infinity`` and ``1``
        - ``upper_bounds`` -- a list of values ``0`` and ``Infinity``
        """
        nonlocal degenerate
        nonlocal certificate

        while non_negative_cocircuits:
            cocircuit = non_negative_cocircuits.pop()
            lower_bounds_new = copy(lower_bounds)
            upper_bounds_new = copy(upper_bounds)
            for i in cocircuit.support():
                lower_bounds_new[i] = 1
                upper_bounds_new[i] = Infinity

            intervals = Intervals.from_bounds(lower_bounds_new, upper_bounds_new)
            indices_new = indices + [cocircuit.support()]
            matrix_new = Matrix(
                matrix_old.rows() + equal_entries_lists(length, cocircuit.support())
            ).echelon_form()
            # TODO don't use kernel matrix? consider evs in row space
            system = LinearInequalitySystem(matrix_new.right_kernel_matrix().T, intervals)

            if system.certify()[0]:
                if certify:
                    covectors_certificate_support_condition = []
                for sv in intervals_to_sign_vectors(intervals):
                    if not system.with_intervals(sign_vector_to_intervals(sv)).certify()[0]:
                        continue
                    if not any(
                        set(cocircuit.support()).issubset(sv.support())
                        for cocircuit in covectors_support_condition
                    ):
                        degenerate = True
                        if certify:
                            certificate = vector_from_sign_vector(
                                system._evs_generator(kernel=False),
                                sv
                            )
                        return
                    if certify:
                        covectors_certificate_support_condition.append(sv)
                if certify:
                    certificate_support_condition.append(
                        [indices_new, covectors_certificate_support_condition]
                    )

            if system.with_intervals(Intervals.from_bounds(lower_bounds_new, upper_bounds_inf)).certify()[0]:
                if certify:
                    certificates_partial_cover.append(indices_new)
                recursive_degenerate(
                    copy(non_negative_cocircuits),
                    matrix_new,
                    indices_new,
                    lower_bounds_new,
                    upper_bounds_new,
                )
            elif certify:
                certificates_zero_equal_components.append(indices_new)

            if degenerate:
                return
        return

    recursive_degenerate(
        non_negative_cocircuits, kernel_matrix, [], lower_bounds, upper_bounds
    )

    if certify:
        if degenerate:
            return degenerate, certificate
        return degenerate, (
            certificates_zero_equal_components,
            certificates_partial_cover,
            certificate_support_condition,
        )
    return degenerate


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
        if not is_constant(product):
            symbolic_pairs.add((minor1, product))
            continue
        if product > 0:
            positive_found = True
        elif product < 0:
            negative_found = True
        if positive_found and negative_found:
            return False

    return closure_minors_utility(symbolic_pairs, positive_found, negative_found)
