r"""Class for setting up chemical reaction networks with mass-action kinetics."""

#############################################################################
#  Copyright (C) 2024                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from sage.structure.sage_object import SageObject
from sage.matrix.constructor import matrix

from elementary_vectors import kernel_matrix_using_elementary_vectors

from .uniqueness import condition_uniqueness_minors
from .unique_existence import condition_faces, condition_nondegenerate
from .robustness import condition_closure_minors


class GMAKSystem(SageObject):
    r"""
    Class for chemical reaction networks with generalized mass-action kinetics.

    A generalized mass-action system is represented by a (weighted) directed graph
    and stoichiometric and kinetic-order labels of the vertices.

    EXAMPLES::

        sage: G = DiGraph({1: [2], 2: [1, 3], 3: [1], 4: [5], 5: [4]})
        sage: Y = matrix(5, 5, [1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
        sage: Y
        [1 1 0 0 0]
        [0 0 1 0 0]
        [0 0 0 1 0]
        [1 0 0 0 0]
        [0 0 0 0 1]
        sage: var('a, b, c')
        (a, b, c)
        sage: Yt = matrix(5, 5, [a, b, 0, 0, 0, 0, 0, 1, 0, 0, c, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1])
        sage: Yt
        [a b 0 0 0]
        [0 0 1 0 0]
        [c 0 0 1 0]
        [1 0 0 0 0]
        [0 0 0 0 1]
        sage: from sign_vector_conditions import *
        sage: crn = GMAKSystem(G, Y, Yt)
        sage: crn.incidence_matrix()
        [-1  1  0  1  0  0]
        [ 1 -1 -1  0  0  0]
        [ 0  0  1 -1  0  0]
        [ 0  0  0  0 -1  1]
        [ 0  0  0  0  1 -1]
        sage: crn.source_matrix()
        [1 0 0 0 0 0]
        [0 1 1 0 0 0]
        [0 0 0 1 0 0]
        [0 0 0 0 1 0]
        [0 0 0 0 0 1]
        sage: crn.stoichiometric_matrix
        [-1 -1  1  0  0]
        [ 0  0 -1  1  0]
        [-1  0  0  0  1]
        sage: crn.kinetic_order_matrix
        [-a -b  1  0  0]
        [ c  0 -1  1  0]
        [-1  0  0  0  1]
        sage: crn.stoichiometric_matrix_kernel
        [1 0 1 1 1]
        [0 1 1 1 0]
        sage: crn.kinetic_order_matrix_kernel
        [    1     0     a a - c     1]
        [    0     1     b     b     0]
        sage: crn.has_robust_CBE()
        [{a > 0, a - c > 0, b > 0}]
        sage: crn.has_at_most_1_CBE()
        [{a >= 0, a - c >= 0, b >= 0}]
    """

    def __init__(self, graph, stoichiometric_labels, kinetic_order_labels):
        self.graph = graph
        self.stoichiometric_labels = stoichiometric_labels
        self.kinetic_order_labels = kinetic_order_labels
        self.stoichiometric_matrix = self._stoichiometric_matrix()
        self.kinetic_order_matrix = self._kinetic_order_matrix()
        self.stoichiometric_matrix_kernel = self._stoichiometric_matrix_kernel()
        self.kinetic_order_matrix_kernel = self._kinetic_order_matrix_kernel()

    # def _repr_(self):
    #     return graph

    def incidence_matrix(self):
        return self.graph.incidence_matrix()

    def source_matrix(self):
        return matrix(
            (1 if value == -1 else 0 for value in row)
            for row in self.graph.incidence_matrix()
        )

    def _stoichiometric_matrix(self):
        M = self.incidence_matrix().T * self.stoichiometric_labels
        return M.matrix_from_rows(M.pivot_rows())

    def _kinetic_order_matrix(self):
        M = self.incidence_matrix().T * self.kinetic_order_labels
        return M.matrix_from_rows(M.pivot_rows())

    def _stoichiometric_matrix_kernel(self):
        return kernel_matrix_using_elementary_vectors(self.stoichiometric_matrix)

    def _kinetic_order_matrix_kernel(self):
        return kernel_matrix_using_elementary_vectors(self.kinetic_order_matrix)

    def has_robust_CBE(self):
        r"""Check whether there is a unique positive CBE with regards to small perturbations."""
        return condition_closure_minors(
            self.stoichiometric_matrix_kernel, self.kinetic_order_matrix_kernel
        )

    def has_at_most_1_CBE(self):
        r"""Check whether there is at most one positive CBE."""
        return condition_uniqueness_minors(
            self.stoichiometric_matrix_kernel, self.kinetic_order_matrix_kernel
        )

    def condition_faces(self):
        r"""Check whether the system satisfies the face condition for existence of a unique positive CBE."""
        return condition_faces(
            self.stoichiometric_matrix_kernel, self.kinetic_order_matrix_kernel
        )

    def are_subspaces_nondegenerate(self):
        r"""Check whether the system satisfies the nondegenerate condition for existence of a unique positive CBE."""
        return condition_nondegenerate(
            self.stoichiometric_matrix_kernel, self.kinetic_order_matrix_kernel
        )
