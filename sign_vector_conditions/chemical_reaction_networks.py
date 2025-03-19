r"""Class for setting up chemical reaction networks with mass-action kinetics."""

#############################################################################
#  Copyright (C) 2025                                                       #
#          Marcus S. Aichmayr (aichmayr@mathematik.uni-kassel.de)           #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from __future__ import annotations

from copy import copy
from sage.graphs.digraph import DiGraph
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

    EXAMPLES:

    We define a chemical reaction network with generalized mass-action kinetics involving 5 complexes and 2 connected components::

        sage: from sign_vector_conditions import *
        sage: edge_destinations = [[1], [0, 2], [0], [4], [3]]
        sage: Y = matrix([[1, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [1, 0, 0, 0, 0], [0, 0, 0, 0, 1]])
        sage: Y
        [1 1 0 0 0]
        [0 0 1 0 0]
        [0 0 0 1 0]
        [1 0 0 0 0]
        [0 0 0 0 1]
        sage: var('a, b, c')
        (a, b, c)
        sage: Yt = matrix([[a, b, 0, 0, 0], [0, 0, 1, 0, 0], [c, 0, 0, 1, 0], [1, 0, 0, 0, 0], [0, 0, 0, 0, 1]])
        sage: Yt
        [a b 0 0 0]
        [0 0 1 0 0]
        [c 0 0 1 0]
        [1 0 0 0 0]
        [0 0 0 0 1]
        sage: crn = GMAKSystem(edge_destinations, Y, Yt)

    We compute the incidence and source matrices of the directed graph::

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

    We describe the stoichiometric and kinetic-order subspaces using matrices::

        sage: crn.stoichiometric_matrix
        [-1 -1  1  0  0]
        [ 0  0 -1  1  0]
        [-1  0  0  0  1]
        sage: crn.kinetic_order_matrix
        [-a -b  1  0  0]
        [ c  0 -1  1  0]
        [-1  0  0  0  1]
        sage: crn.stoichiometric_kernel_matrix
        [1 0 1 1 1]
        [0 1 1 1 0]
        sage: crn.kinetic_order_kernel_matrix
        [    1     0     a a - c     1]
        [    0     1     b     b     0]

    We check some conditions for our system::

        sage: crn.are_deficiencies_zero()
        True
        sage: crn.is_weakly_reversible()
        True
        sage: crn(a=2, b=1, c=1).has_robust_CBE()
        True
        sage: crn.has_robust_CBE() # random order
        [{a > 0, a - c > 0, b > 0}]
        sage: crn.has_at_most_1_CBE() # random order
        [{a >= 0, a - c >= 0, b >= 0}]
    """
    def __init__(self, edge_destinations: list[list[int]], stoichiometric_labels, kinetic_order_labels, set_matrices=True) -> None:
        r"""
        Initialize a chemical reaction network with generalized mass-action kinetics.

        INPUT:

        - ``edge_destinations`` -- a list of lists of integers representing the directed graph
        - ``stoichiometric_labels`` -- a matrix of stoichiometric labels
        - ``kinetic_order_labels`` -- a matrix of kinetic-order labels
        - ``set_matrices`` -- whether to compute the stoichiometric and kinetic-order matrices

        The vertices of the graph are ``0, ..., n - 1``.
        The ``i``th list in ``edge_destinations`` contains the destinations of the vertex ``i``.
        The ``i``th row of ``stoichiometric_labels`` and ``kinetic_order_labels`` contain the labels of the vertex ``i``.
        """
        self.edge_destinations = edge_destinations
        self.graph = DiGraph(dict(enumerate(edge_destinations)))
        self.stoichiometric_labels = stoichiometric_labels
        self.kinetic_order_labels = kinetic_order_labels
        if set_matrices:
            self.stoichiometric_matrix = self._stoichiometric_matrix()
            self.kinetic_order_matrix = self._kinetic_order_matrix()
            self.stoichiometric_kernel_matrix = self._stoichiometric_kernel_matrix()
            self.kinetic_order_kernel_matrix = self._kinetic_order_kernel_matrix()
        else:
            self.stoichiometric_matrix = None
            self.kinetic_order_matrix = None
            self.stoichiometric_kernel_matrix = None
            self.kinetic_order_kernel_matrix = None

    def _repr_(self) -> str:
        return f"System of GMAK with {self.stoichiometric_labels.nrows()} reactions and {self.stoichiometric_labels.ncols()} species"

    def __copy__(self) -> GMAKSystem:
        new = GMAKSystem(self.edge_destinations, self.stoichiometric_labels, self.kinetic_order_labels, set_matrices=False)
        new.stoichiometric_matrix = self.stoichiometric_matrix
        new.kinetic_order_matrix = self.kinetic_order_matrix
        new.stoichiometric_kernel_matrix = self.stoichiometric_kernel_matrix
        new.kinetic_order_kernel_matrix = self.kinetic_order_kernel_matrix
        return new

    def __call__(self, **kwargs) -> GMAKSystem:
        new = copy(self)
        for attribute in [
            "stoichiometric_labels",
            "kinetic_order_labels",
            "stoichiometric_matrix",
            "kinetic_order_matrix",
            "stoichiometric_kernel_matrix",
            "kinetic_order_kernel_matrix",
        ]:
            try:
                setattr(new, attribute, getattr(self, attribute)(**kwargs))
            except TypeError:
                pass
        return new

    def incidence_matrix(self):
        r"""Return the incidence matrix of the graph."""
        return self.graph.incidence_matrix()

    def source_matrix(self):
        r"""Return the source matrix of the graph."""
        return matrix((1 if value == -1 else 0 for value in row) for row in self.graph.incidence_matrix())

    def number_of_species(self) -> int:
        r"""Return the number of species."""
        return self.stoichiometric_matrix.ncols()

    def deficiency_stoichiometric(self):
        r"""Return the stoichiometric deficiency."""
        return self.graph.num_verts() - self.graph.connected_components_number() - self.stoichiometric_matrix.rank()

    def deficiency_kinetic_order(self):
        r"""Return the kinetic-order deficiency."""
        return self.graph.num_verts() - self.graph.connected_components_number() - self.kinetic_order_matrix.rank()

    def _stoichiometric_matrix(self):
        M = self.incidence_matrix().T * self.stoichiometric_labels
        return M.matrix_from_rows(M.pivot_rows())

    def _kinetic_order_matrix(self):
        M = self.incidence_matrix().T * self.kinetic_order_labels
        return M.matrix_from_rows(M.pivot_rows())

    def _stoichiometric_kernel_matrix(self):
        return kernel_matrix_using_elementary_vectors(self.stoichiometric_matrix)

    def _kinetic_order_kernel_matrix(self):
        return kernel_matrix_using_elementary_vectors(self.kinetic_order_matrix)

    def are_deficiencies_zero(self) -> bool:
        r"""Return whether both deficiencies are zero."""
        return self.deficiency_stoichiometric() == self.deficiency_kinetic_order() == 0

    def is_weakly_reversible(self) -> bool:
        r"""Return whether each component of the system is strongly connected."""
        return all(g.is_strongly_connected() for g in self.graph.connected_components_subgraphs())

    def has_robust_CBE(self):
        r"""Check whether there is a unique positive CBE with regards to small perturbations."""
        return condition_closure_minors(self.stoichiometric_kernel_matrix, self.kinetic_order_kernel_matrix)

    def has_at_most_1_CBE(self):
        r"""Check whether there is at most one positive CBE."""
        return condition_uniqueness_minors(self.stoichiometric_kernel_matrix, self.kinetic_order_kernel_matrix)

    def condition_faces(self) -> bool:
        r"""Check whether the system satisfies the face condition for existence of a unique positive CBE."""
        return condition_faces(self.stoichiometric_kernel_matrix, self.kinetic_order_kernel_matrix)

    def are_subspaces_nondegenerate(self) -> bool:
        r"""Check whether the system satisfies the nondegenerate condition for existence of a unique positive CBE."""
        return condition_nondegenerate(self.stoichiometric_kernel_matrix, self.kinetic_order_kernel_matrix)
