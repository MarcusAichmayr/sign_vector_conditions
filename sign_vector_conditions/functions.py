r"""
In this module, we compute polynomial and exponential maps corresponding to given matrices.

EXAMPLES::

    sage: from sign_vector_conditions.functions import *

We define some matrices and a vector c::

    sage: var('x1, x2')
    (x1, x2)
    sage: W = matrix([[1,0,-1],[0,1,-1]])
    sage: W
    [ 1  0 -1]
    [ 0  1 -1]
    sage: Wt = matrix([[1,0,-1],[0,1,0]])
    sage: Wt
    [ 1  0 -1]
    [ 0  1  0]
    sage: c = vector([1,2,4])
    sage: c
    (1, 2, 4)

Next, we compute the polynomial map corresponding to ``W``, ``Wt`` and ``c``::

    sage: fc = f_pol(W, Wt, c)
    sage: fc(x1, x2)
    (x1 - 4/x1, 2*x2 - 4/x1)
    sage: fc(1,2)
    (-3, 0)

We can also omit the argument ``c``.
In this case, ``f_pol`` uses the vector that is one at every component::

    sage: fc = f_pol(W, Wt)
    sage: fc(x1, x2)
    (x1 - 1/x1, x2 - 1/x1)
    sage: fc(1,2)
    (0, 1)

Similarly, we can compute the corresponding exponential map::

    sage: Fc = f_exp(W, Wt, c)
    sage: Fc(x1, x2)
    (-4*e^(-x1) + e^x1, -4*e^(-x1) + 2*e^x2)
    sage: Fc(1,2)
    (e - 4*e^(-1), 2*e^2 - 4*e^(-1))
    sage: Fc = f_exp(W, Wt)
    sage: Fc(x1, x2)
    (-e^(-x1) + e^x1, -e^(-x1) + e^x2)
    sage: Fc(1,2)
    (e - e^(-1), e^2 - e^(-1))
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

from sage.misc.misc_c import prod

from sage.functions.log import exp
from sage.modules.free_module_element import vector
from sage.matrix.special import ones_matrix


def f_exp_pol(W, Wt, c=None, mode="exp"):
    r"""
    Return the exponential or polynomial map determined by the matrices ``W`` and ``Wt``.

    INPUT:

    - ``W`` -- a matrix

    - ``Wt`` -- a matrix

    - ``c`` -- a vector (optional)

    - ``mode`` -- a string. Either 'exp' (default) or 'pol'.

    OUTPUT:

    - If ``c`` is omitted, the result take the vector consisting of ones.
    """
    if W.dimensions() != Wt.dimensions():
        raise ValueError('Matrices must have same dimensions.')

    if c is None:
        c = vector(ones_matrix(1, W.ncols()))
    elif len(c) != W.ncols():
        raise ValueError('Number of columns and dimension of ``c`` do not match.')

    (d, n) = W.dimensions()

    if mode == "exp":
        def function(*x):
            return sum(
                    c_i * exp(Wt_i.dot_product(vector(x))) * W_i
                    for c_i, W_i, Wt_i in zip(c, W.columns(), Wt.columns())
            )
    elif mode == "pol":
        def function(*x):
            return vector(sum(W[i, j] * c[j] * prod([x[k]**Wt[k, j] for k in range(d)])
                    for j in range(n)) for i in range(d)
            )
    else:
        raise ValueError("Mode '" + mode + "' is not available. Try 'exp' or 'pol'.")
    return function


def f_exp(W, Wt, c=None):
    r"""
    Return the exponential map determined by the matrices ``W`` and ``Wt``.

    INPUT:

    - ``W`` -- a matrix

    - ``Wt`` -- a matrix

    - ``c`` -- a vector (optional)

    OUTPUT:

    - If ``c`` is omitted, the result take the vector consisting of ones.

    EXAMPLES::

        sage: from sign_vector_conditions.functions import f_exp
        sage: W = matrix([[1,0,-1],[0,1,-1]])
        sage: W
        [ 1  0 -1]
        [ 0  1 -1]
        sage: Wt = matrix([[1,0,-1],[0,1,0]])
        sage: Wt
        [ 1  0 -1]
        [ 0  1  0]
        sage: c = vector([1,2,4])
        sage: c
        (1, 2, 4)
        sage: Fc = f_exp(W, Wt, c)
        sage: Fc(1,2)
        (e - 4*e^(-1), 2*e^2 - 4*e^(-1))
        sage: var('x,y')
        (x, y)
        sage: Fc(x,y)
        (-4*e^(-x) + e^x, -4*e^(-x) + 2*e^y)

    We can also omit the argument ``c``::

        sage: Fc = f_exp(W, Wt)
        sage: Fc(1,2)
        (e - e^(-1), e^2 - e^(-1))
        sage: Fc(x,y)
        (-e^(-x) + e^x, -e^(-x) + e^y)
    """
    return f_exp_pol(W, Wt, c, mode="exp")


def f_pol(W, Wt, c=None):
    r"""
    Return the polynomial map determined by the matrices ``W`` and ``Wt``.

    INPUT:

    - ``W`` -- a matrix

    - ``Wt`` -- a matrix

    - ``c`` -- a vector (optional)

    OUTPUT:

    - If ``c`` is omitted, the result take the vector consisting of ones.

    EXAMPLES::

        sage: from sign_vector_conditions.functions import f_pol
        sage: W = matrix([[1,0,-1],[0,1,-1]])
        sage: W
        [ 1  0 -1]
        [ 0  1 -1]
        sage: Wt = matrix([[1,0,-1],[0,1,0]])
        sage: Wt
        [ 1  0 -1]
        [ 0  1  0]
        sage: c = vector([1,2,4])
        sage: fc = f_pol(W, Wt, c)
        sage: fc(1,2)
        (-3, 0)
        sage: var('x,y')
        (x, y)
        sage: fc(x,y)
        (x - 4/x, 2*y - 4/x)

    We can also omit the argument ``c``::

        sage: fc = f_pol(W, Wt)
        sage: fc(1,2)
        (0, 1)
        sage: fc(x,y)
        (x - 1/x, y - 1/x)
    """
    return f_exp_pol(W, Wt, c, mode="pol")