r"""
Conditions for CBE
==================

We consider polynomial and exponential maps given by two matrices as in [MHR19]_
and study conditions for injectivity, bijectivity and robustness.

These matrices also describe a chemical reaction network.
In this context, the conditions are equivalent to conditions for positive complex-balanced equilibria (CBE).

Uniqueness
----------

We define some matrices::

    sage: from sign_crn import *
    sage: S = matrix([[1, 1, 1]])
    sage: St = matrix([[1, 0, 1]])

We want to check whether the corresponding map is injective.
For this purpose, we consider the corresponding oriented matroids::

    sage: from sign_vectors.oriented_matroids import *
    sage: OrientedMatroid(S).vectors()
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
    sage: OrientedMatroid(St).covectors()
    {(000), (+0+), (-0-)}

Since the intersection of theses sets of sign vectors contains only the zero sign vector,
the corresponding map is injective.
The package offers a function to check this condition directly::

    sage: condition_uniqueness_sign_vectors(S, St)
    True

Instead of computing sign vectors, we can instead use maximal minors to check this condition::

    sage: condition_uniqueness_minors(S, St)
    True

Now, we consider another example::

    sage: S = matrix([[1, 1, 1]])
    sage: St = matrix([[1, -1, 1]])

The condition is violated::

    sage: condition_uniqueness_sign_vectors(S, St)
    False

In fact, the intersection of the sets of sign vectors is nontrivial::

    sage: OrientedMatroid(S).dual().faces()
    [{(000)},
     {(0+-), (+-0), (-+0), (+0-), (0-+), (-0+)},
     {(++-), (--+), (-++), (+--), (-+-), (+-+)}]
    sage: OrientedMatroid(St).faces()
    [{(000)}, {(-+-), (+-+)}]

Therefore, the corresponding exponential map is not injective.

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
However, we can still compare the products of the maximal minors of ``S`` and ``St``, that is::

    sage: m1 = S.minors(1)
    sage: m1
    [1, 1, 1]
    sage: m2 = St.minors(1)
    sage: m2
    [a, b, -1]
    sage: [m1[i] * m2[i] for i in range(len(m1))]
    [a, b, -1]

The corresponding map is injective, 
if all these products are smaller than or equal to zero if :math:`a, b \leq 0`.
The function :func:`~condition_uniqueness_minors` also works for matrices with symbolic entries.
In this case, it returns a system of inequalities::

    sage: condition_uniqueness_minors(S, St) # random order
    [{-a >= 0, -b >= 0}]

Existence and uniqueness
------------------------

We consider the following matrices to describe an exponential map::

    sage: S = matrix([[1, 0, -1, 0], [0, 1, 0, -1]])
    sage: S
    [ 1  0 -1  0]
    [ 0  1  0 -1]
    sage: St = matrix([[1, 0, -1, 1], [0, 1, -1, 0]])
    sage: St
    [ 1  0 -1  1]
    [ 0  1 -1  0]

The corresponding map is injective::

    sage: condition_uniqueness_minors(S, St)
    True

To examine bijectivity, we first check the face condition.
For this purpose, we compute the circuits of the corresponding oriented matroids::

    sage: circuits = OrientedMatroid(S).circuits()
    sage: circuits
    {(0+0+), (+0+0), (0-0-), (-0-0)}
    sage: circuits_t = OrientedMatroid(St).circuits()
    sage: circuits_t
    {(---0), (-00+), (+00-), (+++0), (0+++), (0---)}

Here, we are only interested in the positive elements::

    sage: [X for X in circuits if X > 0]
    [(0+0+), (+0+0)]
    sage: [X for X in circuits_t if X > 0]
    [(+++0), (0+++)]

The fact condition is satisfied since every positive circuit of ``St`` has a smaller element in ``S``.
The package also offers a function to check this condition directly::

    sage: condition_faces(S, St)
    True

We need to check a third condition to verify bijectivity.
For this purpose, we consider again the oriented matroid determined by ``S``::

    sage: OrientedMatroid(S).covectors()
    {(0000), (-++-), (+0-0), (0-0+), (+--+), (-0+0), (--++), (0+0-), (++--)}

Since there are no nonnegative covectors, the chemical reaction network has at least one equilibrium.
The package offers a function to check this condition condition::

    sage: condition_nondegenerate(S, St)
    True

Hence, the exponential map is bijective.

Let us consider another example.
We swap the two matrices from before::

    sage: S, St = St, S

Because of symmetry, the map is injective::

    sage: condition_uniqueness_sign_vectors(S, St)
    True

The face condition is violated::

    sage: condition_faces(S, St)
    False

Consequently, the map is not bijective.

Now, we consider a map given by two matrices involving a parameter.
(see Example 20 of [MHR19]_)
Depending on this parameter, the map is bijective::

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

Hence, the map is injective.
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

Robustness of existence and uniqueness
--------------------------------------

We consider the following matrices::

    sage: S = matrix([[1, 0, 1, 0], [0, 0, 0, 1]])
    sage: S
    [1 0 1 0]
    [0 0 0 1]
    sage: St = matrix([[1, 0, 1, 1], [0, 1, 0, -1]])
    sage: St
    [ 1  0  1  1]
    [ 0  1  0 -1]

To check, whether the corresponding map is bijective for small perturbations of ``St``,
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

    sage: condition_closure_sign_vectors(S, St)
    True

There is an equivalent condition.
To verify it, we compute the maximal minors of the two matrices::

    sage: S.minors(2)
    [0, 0, 1, 0, 0, 1]
    sage: St.minors(2)
    [1, 0, -1, -1, -1, -1]

From the output, we see whenever a minor of ``S`` is nonzero,
the corresponding minor of ``St`` has the opposite sign.
Hence, this condition is fulfilled::

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

    TESTS::

        sage: from sign_crn import *
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
