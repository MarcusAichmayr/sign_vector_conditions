r"""
Robustness of existence and uniqueness of equilibria

EXAMPLES::

    sage: from sign_vector_conditions import *

Let us consider the following matrices::

    sage: W = matrix([[1,0,-1,0],[0,1,0,0]])
    sage: W
    [ 1  0 -1  0]
    [ 0  1  0  0]
    sage: Wt = matrix([[1,0,-1,0],[0,1,-1,1]])
    sage: Wt
    [ 1  0 -1  0]
    [ 0  1 -1  1]

Therefore, we consider the following exponential map::

    sage: var('x1, x2, c1, c2, c3, c4')
    (x1, x2, c1, c2, c3, c4)
    sage: c = [c1, c2, c3, c4]
    sage: Fc = f_exp(W, Wt, c)
    sage: Fc(x1, x2)
    (c1*e^x1 - c3*e^(-x1 - x2), c2*e^x2)

To check, whether this map is a diffeomorphism for all ``c > 0``
and all small perturbations of ``Wt``,
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

    sage: var('a,b,c')
    (a, b, c)
    sage: W = matrix([[1,0,-1],[0,c,-1]])
    sage: W
    [ 1  0 -1]
    [ 0  c -1]
    sage: Wt = matrix([[1,0,a],[0,1,b]])
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

    sage: condition_closure_minors(W, Wt)
    [[-b > 0, [c == 0, 'or', c > 0], [c == 0, 'or', -a*c > 0]], 'or', [-b < 0, [c == 0, 'or', c < 0], [c == 0, 'or', -a*c < 0]]]

From this system of equations, we can induce that ``b``  must be non-zero.
Indeed, if ``b = 0``, we obtain::

    sage: condition_closure_minors(W, Wt(b=0))
    False

Therefore, assume for instance ``b = 2``::

    sage: condition_closure_minors(W, Wt(b=2))
    [[c == 0, 'or', c < 0], [c == 0, 'or', -a*c < 0]]

From the output, we see that either ``c = 0`` or ``c < 0``.
Otherwise, we obtain::

    sage: condition_closure_minors(W(c=1), Wt(b=2))
    False

If ``c = 0``, then the result will be ``True``::

    sage: condition_closure_minors(W(c=0), Wt(b=2))
    True

For ``c < 0``, we are left with a condition on ``a``::

    sage: condition_closure_minors(W(c=-1), Wt(b=2))
    [a < 0]

Analogously, we can assume ``b < 0`` and we obtain similar conditions on the other variables.
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

from .utility import normalize
from sign_vectors.oriented_matroids import topes_from_matrix
from sage.functions.generalized import sign
from sage.rings.real_mpfr import RR  # used for casting
from sage.rings.integer_ring import ZZ


def condition_closure_sign_vectors(W, Wt):
    r"""
    Closure condition on sign vectors for robustness

    INPUT:

    - ``W`` -- a matrix with ``n`` columns

    - ``Wt`` -- a matrix with ``n`` columns

    OUTPUT:
    Return whether the oriented matroid corresponding to ``W`` is a subset of the closure
    of the oriented matroid corresponding to ``Wt``.
    """
    tW = topes_from_matrix(W, kernel=True)
    tWt = topes_from_matrix(Wt, kernel=True)
    tWn = normalize(tW)
    # Do not normalize second list of topes: 0++0, +--+
    for X in tWn:
        val = True
        for Y in tWt:
            if X <= Y:
                val = False
                break
        if val:
            return False
    return True


def condition_closure_minors(W, Wt):
    r"""
    Closure condition on maximal minors for robustness

    INPUT:

    - ``W`` -- a matrix

    - ``Wt`` -- a matrix with the same dimensions as ``W``

    OUTPUT:
    Returns whether for each non-zero maximal minor ``m`` of ``W``, we have
    either ``m mt > 0`` for each corresponding maximal minor ``mt`` of ``Wt``
    or ``m mt < 0`` for each corresponding maximal minor ``mt`` of ``Wt``.

    Returns a boolean or a symbolic expression if variables occur.
    """
    if W.dimensions() != Wt.dimensions():
        raise ValueError('Matrices must have same dimensions.')
    m = W.minors(W.nrows())
    mt = Wt.minors(W.nrows())

    def equals_zero(value):
        r"""
        Compare an expression with ``0``.

        INPUT:

        - ``a`` -- a real number or symbolic expression

        OUTPUT:

        - ``True``  -- The number ``a`` is equal to zero.

        - ``False`` -- The number ``a`` is a non-zero real number.

        - ``None``  -- The number ``a`` is a symbolic expression.
        """
        try:
            if RR(value) == 0:
                return True  # equal to zero
            return False  # non-zero real number
        except TypeError:
            return None  # symbolic expression

    sign_product_minors = 0
    symbolic_entries = []  # Depending on s, each entry is ``> 0`` or ``< 0``.
    equalities = []
    # equalities_left, equalities_right will be the ``or`` pairs
    equalities_left = []
    equalities_right = []

    for m_i, mt_i in zip(m, mt):
        if equals_zero(m_i) is None:  # mi symbolic
            if equals_zero(mt_i) is True:  # mti zero, hence, product zero and cannot be > or < 0
                equalities.append(m_i == 0)
            else:
                equalities_left.append(m_i == 0)
                equalities_right.append(m_i * mt_i)
        elif equals_zero(m_i) is False:  # mi non-zero, not symbolic
            if equals_zero(mt_i) is None:  # mti symbolic
                symbolic_entries.append(m_i * mt_i)  # needs > or < later
            elif equals_zero(mt_i) is True:  # mti zero
                return False
            elif sign_product_minors == 0:  # mti non-zero, not symbolic
                sign_product_minors = ZZ(sign(m_i * mt_i))
            else:
                if (sign_product_minors * m_i * mt_i) <= 0:  # cast
                    return False

    # construction of the output list
    if not symbolic_entries and not equalities and not equalities_left:  # hence, equalities_right is also []
        return True
    if not symbolic_entries and not equalities_left:
        return equalities
    if sign_product_minors == 0:
        if symbolic_entries:
            output_equalities_left = []
            output_equalities_right = []
            for l in symbolic_entries:
                output_equalities_left.append(l > 0)
                output_equalities_right.append(l < 0)

        if equalities_left:
            output_disjunction_left = []
            output_disjunction_right = []
            for equalities_left_i, equalities_right_i in zip(equalities_left, equalities_right):
                output_disjunction_left.append([equalities_left_i, 'or', equalities_right_i > 0])
                output_disjunction_right.append([equalities_left_i, 'or', equalities_right_i < 0])

        if not symbolic_entries:
            if len(output_disjunction_left) == 1:  # hence, also len(output_disjunction_right) == 1
                out = [flatten(output_disjunction_left), 'or', flatten(output_disjunction_right)]
            else:
                out = [output_disjunction_left, 'or', output_disjunction_right]
        elif not equalities_left:
            out = [output_equalities_left, 'or', output_equalities_right]
        else:
            out = [output_equalities_left + output_disjunction_left, 'or', output_equalities_right + output_disjunction_right]
    else:
        if sign_product_minors > 0:
            def rel(value):
                return value > 0
        else:
            def rel(value):
                return value < 0

        if symbolic_entries:
            output_equalities = []
            for l in symbolic_entries:
                output_equalities.append(rel(l))
        if equalities_left:
            output_disjunctions = []
            for equalities_left_i, equalities_right_i in zip(equalities_left, equalities_right):
                output_disjunctions.append([equalities_left_i, 'or', rel(equalities_right_i)])
        if not symbolic_entries:
            if len(output_disjunctions) == 1:
                out = flatten(output_disjunctions)
            else:
                out = output_disjunctions
        elif not equalities_left:
            out = output_equalities
        else:
            out = output_equalities + output_disjunctions

    if not equalities:
        return out
    return equalities + [out]
