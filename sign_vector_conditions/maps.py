r"""
Polynomial and exponential functions corresponding to chemical reaction networks.

Using matrices, we can represent polynomial and exponential maps.
To check whether such a map is injective or surjective,
we can use maximal minors of the matrices or sign vectors of the corresponding oriented matroids.
Check out the other modules in this package.
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

        sage: from sign_vector_conditions.maps import f_exp
        sage: W = matrix([[1, 0, -1], [0, 1, -1]])
        sage: W
        [ 1  0 -1]
        [ 0  1 -1]
        sage: Wt = matrix([[1, 0, -1], [0, 1, 0]])
        sage: Wt
        [ 1  0 -1]
        [ 0  1  0]
        sage: c = vector([1, 2, 4])
        sage: c
        (1, 2, 4)

    Next, we compute the corresponding exponential map corresponding::
        
        sage: f = f_exp(W, Wt, c)
        sage: f(1, 2)
        (e - 4*e^(-1), 2*e^2 - 4*e^(-1))
        sage: var('x, y')
        (x, y)
        sage: f(x, y)
        (-4*e^(-x) + e^x, -4*e^(-x) + 2*e^y)

    We can also omit the argument ``c``.
    In this case, ``f_exp`` uses the vector that is one at every component::

        sage: f = f_exp(W, Wt)
        sage: f(1,2)
        (e - e^(-1), e^2 - e^(-1))
        sage: f(x, y)
        (-e^(-x) + e^x, -e^(-x) + e^y)
    """
    if W.dimensions() != Wt.dimensions():
        raise ValueError('Matrices must have same dimensions.')

    if c is None:
        c = vector(ones_matrix(1, W.ncols()))
    elif len(c) != W.ncols():
        raise ValueError('Number of columns and dimension of ``c`` do not match.')

    def f(*x):
        return sum(
                c_i * exp(Wt_i.dot_product(vector(x))) * W_i
                for c_i, W_i, Wt_i in zip(c, W.columns(), Wt.columns())
        )
    return f


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

        sage: from sign_vector_conditions.maps import f_pol
        sage: W = matrix([[1, 0, -1], [0, 1, -1]])
        sage: W
        [ 1  0 -1]
        [ 0  1 -1]
        sage: Wt = matrix([[1, 0, -1], [0, 1, 0]])
        sage: Wt
        [ 1  0 -1]
        [ 0  1  0]
        sage: c = vector([1, 2, 4])

    Next, we compute the corresponding polynomial map::

        sage: f = f_pol(W, Wt, c)
        sage: f(1, 2)
        (-3, 0)
        sage: var('x, y')
        (x, y)
        sage: f(x, y)
        (x - 4/x, 2*y - 4/x)

    We can also omit the argument ``c``.
    In this case, ``f_pol`` uses the vector that is one at every component::

        sage: f = f_pol(W, Wt)
        sage: f(1, 2)
        (0, 1)
        sage: f(x, y)
        (x - 1/x, y - 1/x)
    """
    if W.dimensions() != Wt.dimensions():
        raise ValueError('Matrices must have same dimensions.')

    if c is None:
        c = vector(ones_matrix(1, W.ncols()))
    elif len(c) != W.ncols():
        raise ValueError('Number of columns and dimension of ``c`` do not match.')

    (d, n) = W.dimensions()

    def f(*x):
        return vector(sum(W[i, j] * c[j] * prod([x[k]**Wt[k, j] for k in range(d)])
                for j in range(n)) for i in range(d)
        )
    return f
