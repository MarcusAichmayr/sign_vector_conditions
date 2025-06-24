r"""General equilibria"""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.modules.free_module_element import vector

from sign_vectors import sign_vector

from .utility import non_negative_covectors_from_matrix


def condition_properness(A, B):
    A_bar = A.insert_row(0, vector([1] * A.ncols()))
    B_bar = B.insert_row(0, vector([1] * A.ncols()))
    covectors = non_negative_covectors_from_matrix(B_bar, dual=False)
    for face in non_negative_covectors_from_matrix(A_bar, dual=True):
        if face == 0:
            continue
        if not face >= 0:
            continue
        if sign_vector(1 if Fe == 0 else 0 for Fe in face) in covectors:
            return False
    return True
