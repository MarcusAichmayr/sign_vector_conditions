r"""
Conditions for polynomial and exponential maps
==============================================

We consider polynomial and exponential maps given by two matrices as in [MHR19]_
and study conditions for injectivity, bijectivity and robustness.

These matrices also describe a chemical reaction network.
In this context, the conditions are equivalent to conditions for positive complex-balanced equilibria (CBE).

Uniqueness
----------

We define matrices that describe an exponential map::

    sage: from sign_crn import *
    sage: S = matrix([[1, 1, 1]])
    sage: St = matrix([[1, 0, 1]])

The package uses maximal minors to study injectivity::

    sage: uniqueness_condition(S, St)
    True

Therefore, the map is injective.
Instead, we can consider the oriented matroids determined by these matrices::

    sage: from sign_crn.conditions import uniqueness_condition_sign_vectors
    sage: uniqueness_condition_sign_vectors(S, St)
    True

Now, we consider another example::

    sage: S = matrix([[1, 1, 1]])
    sage: St = matrix([[1, -1, 1]])

The condition is violated::

    sage: uniqueness_condition(S, St)
    False

Therefore, the corresponding exponential map is not injective.

Finally, we consider an example with parameters :math:`a, b \in \mathbb{R}`::

    sage: var("a, b")
    (a, b)
    sage: S = matrix([[1, 1, 1]])
    sage: S
    [1 1 1]
    sage: St = matrix([[a, b, -1]])
    sage: St
    [ a  b -1]

Here, the function returns a system of inequalities::

    sage: uniqueness_condition(S, St) # random order
    [{-a >= 0, -b >= 0}]

Hence, the map is injective if :math:`a, b \leq 0`.

Note that we cannot apply :func:`~uniqueness_condition_sign_vectors` because of the parameters.

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

    sage: uniqueness_condition(S, St)
    True

To examine bijectivity, we first check the face condition::

    sage: face_condition(S, St)
    True

Finally, we apply the degeneracy condition::

    sage: degeneracy_condition(S, St)
    False

Hence, the exponential map is bijective.

Let us consider another example.
We swap the two matrices from before::

    sage: S, St = St, S

Because of symmetry, the map is injective::

    sage: uniqueness_condition_sign_vectors(S, St)
    True

The face condition is violated::

    sage: face_condition(S, St)
    False

Consequently, the map is not bijective.

Now, we consider a map involving a parameter.
(see Example 20 of [MHR19]_)::

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

    sage: uniqueness_condition_sign_vectors(S, St)
    True

Hence, the map is injective.
Also the face condition is satisfied::

    sage: face_condition(S, St)
    True

For specific values of :math:`a`, the pair of subspaces
determined by kernels of the matrices is nondegenerate.
This is exactly the case for :math:`a \in (0, 1) \cup (1, 2)`.
We demonstrate this for specific values::

    sage: degeneracy_condition(S, St(a=1/2))
    False
    sage: degeneracy_condition(S, St(a=3/2))
    False

On the other hand, this condition does not hold if
:math:`a \in {1} \cup [2, \infty)`::

    sage: degeneracy_condition(S, St(a=1))
    True

To certify the result, we call::

    sage: degeneracy_condition(S, St(a=1), certify=True)
    (True, (1, 1, 0, 0, -1, 1))

Hence, the positive support of the vector ``v = (1, 1, 0, 0, -1, 1)`` of ``St``
can be covered by a sign vector ``(++000+)`` corresponding to ``ker(S)``.
Further, ``v`` does not satisfy the support condition::

    sage: degeneracy_condition(S, St(a=2))
    True
    sage: degeneracy_condition(S, St(a=3))
    True

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

To study robustness of the corresponding map,
we consider again a condition involving maximal minors::

    sage: closure_condition(S, St)
    True

Hence, the map is bijective for small perturbations of ``St``.
There is also an equivalent condition using sign vectors::

    sage: from sign_crn.conditions import closure_condition_sign_vectors
    sage: closure_condition_sign_vectors(S, St)
    True

Now, we consider an example involving parameters::

    sage: var("a, b, c")
    (a, b, c)
    sage: S = matrix([[c, 1, c]])
    sage: S
    [c 1 c]
    sage: St = matrix([[a, b, -1]])
    sage: St
    [ a  b -1]

We obtain the following condition on the variables::

    sage: closure_condition(S, St) # random
    [{-b > 0, c == 0},
     {-b < 0, c == 0},
     {-b > 0, c > 0, -a*c > 0},
     {-b < 0, c < 0, -a*c < 0}]

Thus, there are four possibilities to set the variables:
From the first two sets of conditions, we see that the closure condition is satisfied
if :math:`c` is zero and :math:`b` is nonzero.
The closure condition is also satisfied if :math:`a` and :math:`b` are negative and :math:`c` is positive
or if :math:`a` and :math:`b` are positive and :math:`c` is negative.
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


def uniqueness_condition_sign_vectors(stoichiometric_matrix: Matrix, kinetic_order_matrix: Matrix) -> bool:
    r"""
    Uniqueness condition for existence of an equilibrium using sign vectors.

    OUTPUT:
    Return whether there exists at most one equilibrium.

    .. SEEALSO::

        :func:`uniqueness_condition`

    .. NOTE::

        This implementation is inefficient and should not be used for large examples.
        Instead, use :func:`~uniqueness_condition`.

    EXAMPLES::

        sage: from sign_crn.conditions import uniqueness_condition_sign_vectors
        sage: S = matrix([[1, 1, 1]])
        sage: S
        [1 1 1]
        sage: St = matrix([[1, 0, 1]])
        sage: St
        [1 0 1]
        sage: uniqueness_condition_sign_vectors(S, St)
        True
        sage: S = matrix([[1, 0, -1], [0, 1, -1]])
        sage: S
        [ 1  0 -1]
        [ 0  1 -1]
        sage: St = matrix([[1, 0, -1], [0, 1, 1]])
        sage: St
        [ 1  0 -1]
        [ 0  1  1]
        sage: uniqueness_condition_sign_vectors(S, St)
        False

    TESTS::

        sage: from sign_crn.conditions import uniqueness_condition_sign_vectors
        sage: A = identity_matrix(3)
        sage: B = A # kernel of B is empty
        sage: uniqueness_condition_sign_vectors(A, B)
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


def uniqueness_condition(stoichiometric_matrix: Matrix, kinetic_order_matrix: Matrix):
    r"""
    Uniqueness condition for existence of an equilibrium using maximal minors.

    OUTPUT:
    Return whether there exists at most one equilibrium.
    If the result depends on variables, a list of sets is returned.
    The condition holds if the inequalities in exactly one of these sets are satisfied.

    .. SEEALSO::

        :func:`uniqueness_condition_sign_vectors`

    .. NOTE::

        The matrices need to have maximal rank and the same dimensions.
        Otherwise, a ``ValueError`` is raised.

    TESTS::

        sage: from sign_crn import *
        sage: var("a, b")
        (a, b)
        sage: S = matrix([[1, 0, -1], [0, 1, -1]])
        sage: S
        [ 1  0 -1]
        [ 0  1 -1]
        sage: St = matrix([[1, 0, a], [0, 1, b]])
        sage: St
        [1 0 a]
        [0 1 b]
        sage: uniqueness_condition(S, St) # random order
        [{-a >= 0, -b >= 0}]
        sage: conditions = uniqueness_condition(S, St)[0]
        sage: conditions # random order
        sage: (-a >= 0) in conditions and (-b >= 0) in conditions
        True
        sage: S = matrix([[a, 0, 1, 0], [0, 1, -1, 0], [0, 0, 0, 1]])
        sage: St = matrix([[1, 0, 0, -1], [0, b, 1, 1], [0, 0, a, 1]])
        sage: uniqueness_condition(S, St) # random order
        [{(a - 1)*a >= 0, a*b >= 0}, {(a - 1)*a <= 0, a*b <= 0}]
        sage: len(_), len(_[0])
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


def face_condition(stoichiometric_matrix: Matrix, kinetic_order_matrix: Matrix) -> bool:
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


def degeneracy_condition(stoichiometric_matrix: Matrix, kinetic_order_matrix: Matrix, certify: bool = False) -> bool:
    r"""
    Return whether a pair of subspaces given by matrices is degenerate.

    This condition is about whether all positive equal components of a vector in ``St``
    can be covered by covectors corresponding to the kernel of ``S``.

    OUTPUT:
    a boolean

    If ``certify`` is true, a list is returned to certify the result.
    (see the examples)

    EXAMPLES:

    We consider the following matrices::

        sage: from sign_crn import *
        sage: S = matrix([[1, 0, -1, 0], [0, 1, 0, -1]])
        sage: St = matrix([[1, 0, 0, 1], [0, 1, 0, 1]])
        sage: degeneracy_condition(S, St)
        False

    Next, we certify the result.
    The corresponding subspaces are trivially nondegenerate
    since there are no nonnegative covectors in the kernel of ``S``::

        sage: degeneracy_condition(S, St, certify=True)
        (False, 'no nonnegative covectors')

    Now, we consider an example of degenerate subspaces::

        sage: S = matrix([[1, 1, 0]])
        sage: St = matrix([[0, 0, 1]])
        sage: degeneracy_condition(S, St, certify=True)
        (True, (1, 1, 0))

    The resulting vector lies in the row space of ``St``.
    The nonnegative covector ``(++0)`` in the kernel of ``S`` covers the first two equal components.

    We have another example for nondegenerate subspaces::

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
        sage: degeneracy_condition(S, St, certify=True)
        (False, ([[[1, 2, 3]], [[0, 2, 3]]], [[[2, 4]]], []))

    The certificate tells us that there is no vector in the row space of ``St``
    with positive support on the components ``0, 2, 3`` and ``1, 2, 3``.
    Positive equal components can partially be covered by a covector ``(00+0+)``
    which corresponds to ``[[2, 4]]``.
    However, it is impossible to fully cover the positive support.

    In the next example, there exists a partial cover::

        sage: S = matrix([[1, -1, 0, 0], [0, 0, 1, 1]])
        sage: St = matrix([[1, 0, 0, 1], [0, 1, 0, 1]])
        sage: degeneracy_condition(S, St, certify=True)
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


def closure_condition_sign_vectors(stoichiometric_matrix: Matrix, kinetic_order_matrix: Matrix) -> bool:
    r"""
    Closure condition for robustness using sign vectors.

    OUTPUT:
    Return whether the closure condition for robustness regarding small perturbations is satisfied.

    .. SEEALSO::

        :func:`closure_condition`

    .. NOTE::

        This implementation is inefficient and should not be used for large examples.
        Instead, use :func:`~closure_condition`.
    """
    topes = OrientedMatroid(kinetic_order_matrix).topes()
    for covector1 in OrientedMatroid(stoichiometric_matrix).topes():
        if not any(covector1 <= covector2 for covector2 in topes):
            return False
    return True


def closure_condition(stoichiometric_matrix: Matrix, kinetic_order_matrix: Matrix):
    r"""
    Closure condition for robustness of bijectivity using maximal minors.

    OUTPUT:
    Return whether the closure condition for robustness regarding small perturbations is satisfied.
    If the result depends on variables, a list of sets is returned.
    The condition holds if the inequalities in (at least) one of these sets are satisfied.

    .. SEEALSO::

        :func:`closure_condition_sign_vectors`

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
