r"""This module offers utility functions that are used by the other functions of this package."""

#############################################################################
#  Copyright (C) 2023                                                       #
#                Marcus Aichmayr (aichmayr@mathematik.uni-kassel.de)        #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from elementary_vectors import elementary_vectors
from sign_vectors import sign_vector, zero_sign_vector
from sign_vectors.oriented_matroids import cocircuits_from_matrix


def normalize(L):
    r"""
    Compute a normalized set of sign vectors.

    INPUT:

    - ``L`` -- an iterable of sign vectors

    OUTPUT:
    a set of sign vectors consisting only of those sign vectors of ``L``
    that are normalized

    .. NOTE::

        A sign vector is normalized if it is the zero sign vector
        or the first non-zero entry is positive.

    EXAMPLES::

        sage: from sign_vectors import sign_vector
        sage: L = [
        ....:     sign_vector('++0'), sign_vector('--0'), sign_vector('-0+'),
        ....:     sign_vector('+0-'), sign_vector('0++'), sign_vector('0--')
        ....: ]
        sage: L
        [(++0), (--0), (-0+), (+0-), (0++), (0--)]
        sage: from sign_vector_conditions.utility import normalize
        sage: normalize(L)
        {(0++), (++0), (+0-)}
        sage: L = [sign_vector('000'), sign_vector('-++')]
        sage: normalize(L)
        {(000)}
    """
    def is_normalized(X):
        for Xi in X:
            if Xi > 0:
                return True
            if Xi < 0:
                return False
        return True

    return set(X for X in L if is_normalized(X))


def pos_cocircuits_from_matrix(M, kernel=True):
    r"""
    Compute the positive cocircuits determined by a given matrix.

    INPUT:

    - ``A`` -- a matrix with real arguments.

    - ``kernel`` -- a boolean (default: ``True``)

    OUTPUT:

    - If ``kernel`` is true, returns a set of positive cocircuits determined by the kernel of the matrix ``A``.

    - If ``kernel`` is false, returns a set of positive cocircuits determined by the rows of the matrix ``A``.

    EXAMPLES::

        sage: M = matrix([[2, -1, -1, 0]])
        sage: M
        [ 2 -1 -1  0]
        sage: from sign_vectors.oriented_matroids import cocircuits_from_matrix
        sage: cocircuits_from_matrix(M)
        {(--00), (000-), (0-+0), (+0+0), (++00), (-0-0), (000+), (0+-0)}
        sage: from sign_vector_conditions.utility import pos_cocircuits_from_matrix
        sage: pos_cocircuits_from_matrix(M)
        {(+0+0), (++00), (000+)}
    """
    return set(X for X in cocircuits_from_matrix(M, kernel=kernel) if X > 0)


def pos_covectors_from_matrix(M, kernel=True):
    r"""
    Use a list of cocircuits to compute all covectors of the corresponding oriented matroid.

    INPUT:

    - ``M`` -- a matrix

    - ``kernel`` -- a boolean (default: ``True``)

    OUTPUT:

    - a list of all covectors of the oriented matroid.

      - If ``kernel`` is true, returns a list of cocircuits determined by the kernel of the matrix ``A``.

      - If ``kernel`` is false, returns a list of cocircuits determined by the rows of the matrix ``A``.

    EXAMPLES::

        sage: M = matrix([[2, -1, -1, 0]])
        sage: M
        [ 2 -1 -1  0]
        sage: from sign_vectors.oriented_matroids import covectors_from_matrix
        sage: covectors_from_matrix(M)
        {(0000),
         (0-+0),
         (0-+-),
         (++++),
         (--0+),
         (000-),
         (----),
         (++0-),
         (+-++),
         ...
         (-+-0),
         (--0-)}
        sage: from sign_vector_conditions.utility import pos_covectors_from_matrix
        sage: pos_covectors_from_matrix(M)
        [(+0++), (+++0), (++++), (+0+0), (++00), (++0+), (000+)]
    """
    cocircuits = [X for X in cocircuits_from_matrix(M, kernel=kernel) if not X < 0]

    if not cocircuits:
        raise ValueError('List of cocircuits is empty.')
    output = set()
    new_elements = {zero_sign_vector(M.ncols())}
    while new_elements:
        Y = new_elements.pop()
        for X in cocircuits:
            if X >= 0 and not X <= Y: # otherwise Z = X.compose(Y) = Y in F
                Z = X.compose(Y)
                if Z not in output and Z >= 0:
                    output.add(Z)
                    new_elements.add(Z)
    return list(output)