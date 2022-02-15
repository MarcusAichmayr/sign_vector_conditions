r"""This module offers utility functions that are used by the other functions of this package."""

#############################################################################
#  Copyright (C) 2022                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from elementary_vectors import elementary_vectors
from sign_vectors import sign_vector, zero_sign_vector


def normalize(L):
    r"""
    Compute a normalized list of sign vectors.

    .. NOTE::

        A sign vector is normalized if it is the zero sign vector or the first non-zero entry is positive.
    """
    L_new = []
    for X in L:
        val = True
        for Xi in X:
            if Xi > 0:
                break
            elif Xi < 0:
                val = False
                break
        if val:
            L_new.append(X)
    return L_new


def pos_cocircuits_from_matrix(A, kernel=False):
    r"""
    Compute a list of positive cocircuits determined by a given matrix.

    INPUT:

    - ``A`` -- a matrix with real arguments.

    - ``kernel`` -- a boolean (default: False)

    OUTPUT:

    - If ``kernel`` is false, returns a list of positive cocircuits determined by the rows of the matrix ``A``.

    - If ``kernel`` is true, returns a list of positive cocircuits determined by the kernel of the matrix ``A``.
    """
    L = elementary_vectors(A, kernel=kernel)
    return [sign_vector(v) for v in L if sign_vector(v) > 0] + [sign_vector(-v) for v in L if sign_vector(-v) > 0]


def pos_covectors_from_matrix(A, kernel=False):
    r"""
    Use a list of cocircuits to compute all covectors of the corresponding oriented matroid.

    INPUT:

    - ``L`` -- a list of cocircuits of an oriented matroid.

    - ``kernel`` -- a boolean (default: False)

    OUTPUT:

    - a list of all covectors of the oriented matroid.

      - If ``kernel`` is false, returns a list of cocircuits determined by the rows of the matrix ``A``.

      - If ``kernel`` is true, returns a list of cocircuits determined by the kernel of the matrix ``A``.
    """
    ev = elementary_vectors(A, kernel=kernel)
    L = [sign_vector(v) for v in ev if not sign_vector(v) < 0] + [sign_vector(-v) for v in ev if not sign_vector(-v) < 0]

    if not L:
        raise ValueError('List of cocircuits is empty.')
    n = L[0].length()
    F = []
    F_new = [zero_sign_vector(n)]
    while F_new != []:
        Y = F_new.pop()
        for X in L:
            if X >= 0:
                if not X <= Y:  # otherwise Z = X.compose(Y) = Y in F
                    Z = X.compose(Y)
                    if Z not in F:
                        if Z >= 0:
                            F.append(Z)
                            F_new.append(Z)
    return F
