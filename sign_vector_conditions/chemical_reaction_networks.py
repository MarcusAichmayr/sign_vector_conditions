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

        sage: crn.matrix_stoichiometric
        [-1 -1  1  0  0]
        [ 0  0 -1  1  0]
        [-1  0  0  0  1]
        sage: crn.matrix_kinetic_order
        [-a -b  1  0  0]
        [ c  0 -1  1  0]
        [-1  0  0  0  1]
        sage: crn.kernel_matrix_stoichiometric
        [1 0 1 1 1]
        [0 1 1 1 0]
        sage: crn.kernel_matrix_kinetic_order
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
    def __init__(self, edge_destinations: list[list[int]], matrix_of_complexes_stoichiometric, matrix_of_complexes_kinetic_order, set_matrices=True) -> None:
        r"""
        Initialize a chemical reaction network with generalized mass-action kinetics.

        INPUT:

        - ``edge_destinations`` -- a list of lists of integers representing the directed graph
        - ``matrix_of_complexes_stoichiometric`` -- a matrix of stoichiometric labels
        - ``matrix_of_complexes_kinetic_order`` -- a matrix of kinetic-order labels
        - ``set_matrices`` -- whether to compute the stoichiometric and kinetic-order matrices

        The vertices of the graph are ``0, ..., n - 1``.
        The ``i``th list in ``edge_destinations`` contains the destinations of the vertex ``i``.
        The ``i``th row of ``matrix_of_complexes_stoichiometric`` and ``matrix_of_complexes_kinetic_order`` contain the labels of the vertex ``i``.
        """
        self.edge_destinations = edge_destinations
        self.graph = DiGraph(dict(enumerate(edge_destinations)))
        self.matrix_of_complexes_stoichiometric = matrix_of_complexes_stoichiometric
        self.matrix_of_complexes_kinetic_order = matrix_of_complexes_kinetic_order
        if set_matrices:
            self.set_matrix_stoichiometric()
            self.set_matrix_kinetic_order()
            self.set_kernel_matrix_stoichiometric()
            self.set_kernel_matrix_kinetic_order()
        else:
            self.matrix_stoichiometric = None
            self.matrix_kinetic_order = None
            self.kernel_matrix_stoichiometric = None
            self.kernel_matrix_kinetic_order = None

    def _repr_(self) -> str:
        return f"System of GMAK with {self.matrix_of_complexes_stoichiometric.nrows()} reactions and {self.matrix_of_complexes_stoichiometric.ncols()} species"

    def __copy__(self) -> GMAKSystem:
        new = GMAKSystem(self.edge_destinations, self.matrix_of_complexes_stoichiometric, self.matrix_of_complexes_kinetic_order, set_matrices=False)
        new.matrix_stoichiometric = self.matrix_stoichiometric
        new.matrix_kinetic_order = self.matrix_kinetic_order
        new.kernel_matrix_stoichiometric = self.kernel_matrix_stoichiometric
        new.kernel_matrix_kinetic_order = self.kernel_matrix_kinetic_order
        return new

    def __call__(self, **kwargs) -> GMAKSystem:
        new = copy(self)
        for attribute in [
            "matrix_of_complexes_stoichiometric",
            "matrix_of_complexes_kinetic_order",
            "matrix_stoichiometric",
            "matrix_kinetic_order",
            "kernel_matrix_stoichiometric",
            "kernel_matrix_kinetic_order",
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
        return matrix((1 if value == -1 else 0 for value in row) for row in self.incidence_matrix())

    def number_of_species(self) -> int:
        r"""Return the number of species."""
        return self.matrix_stoichiometric.ncols()

    def deficiency_stoichiometric(self):
        r"""Return the stoichiometric deficiency."""
        return self.graph.num_verts() - self.graph.connected_components_number() - self.matrix_stoichiometric.rank()

    def deficiency_kinetic_order(self):
        r"""Return the kinetic-order deficiency."""
        return self.graph.num_verts() - self.graph.connected_components_number() - self.matrix_kinetic_order.rank()

    def set_matrix_stoichiometric(self):
        M = self.incidence_matrix().T * self.matrix_of_complexes_stoichiometric
        self.matrix_stoichiometric = M.matrix_from_rows(M.pivot_rows())

    def set_matrix_kinetic_order(self):
        M = self.incidence_matrix().T * self.matrix_of_complexes_kinetic_order
        self.matrix_kinetic_order = M.matrix_from_rows(M.pivot_rows())

    def set_kernel_matrix_stoichiometric(self):
        self.kernel_matrix_stoichiometric = kernel_matrix_using_elementary_vectors(self.matrix_stoichiometric)

    def set_kernel_matrix_kinetic_order(self):
        self.kernel_matrix_kinetic_order = kernel_matrix_using_elementary_vectors(self.matrix_kinetic_order)

    def are_deficiencies_zero(self) -> bool:
        r"""Return whether both deficiencies are zero."""
        return self.deficiency_stoichiometric() == self.deficiency_kinetic_order() == 0

    def is_weakly_reversible(self) -> bool:
        r"""Return whether each component of the system is strongly connected."""
        return all(g.is_strongly_connected() for g in self.graph.connected_components_subgraphs())

    def has_robust_CBE(self):
        r"""Check whether there is a unique positive CBE with regards to small perturbations."""
        return condition_closure_minors(self.kernel_matrix_stoichiometric, self.kernel_matrix_kinetic_order)

    def has_at_most_1_CBE(self):
        r"""Check whether there is at most one positive CBE."""
        return condition_uniqueness_minors(self.kernel_matrix_stoichiometric, self.kernel_matrix_kinetic_order)

    def condition_faces(self) -> bool:
        r"""Check whether the system satisfies the face condition for existence of a unique positive CBE."""
        return condition_faces(self.kernel_matrix_stoichiometric, self.kernel_matrix_kinetic_order)

    def are_subspaces_nondegenerate(self) -> bool:
        r"""Check whether the system satisfies the nondegenerate condition for existence of a unique positive CBE."""
        return condition_nondegenerate(self.kernel_matrix_stoichiometric, self.kernel_matrix_kinetic_order)
