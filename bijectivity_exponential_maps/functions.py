#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
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

def f_pol(W, Wt, c=None):
    r"""
    Returns the polynomial map determined by the matrices ``W`` and ``Wt``.
    
    INPUT:
    
    - ``W`` -- a matrix
    
    - ``Wt`` -- a matrix
    
    - ``c`` -- a vector (optional)
    
    OUTPUT:
    
    - If ``c`` is omitted, the result take the vector consisting of ones.
    """
    if W.dimensions() != Wt.dimensions():
        raise ValueError('Matrices must have same dimensions.')

    if c == None:
        c = vector(ones_matrix(1, W.ncols()))
    elif len(c) != W.ncols():
        raise ValueError('Number of columns and dimension of ``c`` do not match.')
    
    (d, n) = W.dimensions()
    
    def f(x):
        return vector([sum([W[i,j] * c[j] * prod([x[k]**Wt[k,j] for k in range(d)]) for j in range(n)]) for i in range(d)])
    return f


def f_exp(W, Wt, c=None):
    r"""
    Returns the exponential map determined by the matrices ``W`` and ``Wt``.
    
    INPUT:
    
    - ``W`` -- a matrix
    
    - ``Wt`` -- a matrix
    
    - ``c`` -- a vector (optional)
    
    OUTPUT:
    
    - If ``c`` is omitted, the result take the vector consisting of ones.
    """
    if W.dimensions() != Wt.dimensions():
        raise ValueError('Matrices must have same dimensions.')

    if c == None:
        c = vector(ones_matrix(1, W.ncols()))
    elif len(c) != W.ncols():
        raise ValueError('Number of columns and dimension of ``c`` do not match.')
    
    (d, n) = W.dimensions()
    
    def f(x):
        return sum([c[i] * exp(Wt.column(i).dot_product(vector(x))) * W.column(i) for i in range(n)])
    return f
