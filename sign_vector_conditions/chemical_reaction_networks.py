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
        sage: var('a, b, c')
        (a, b, c)
        sage: species = var('A, B, C, D, E')
        sage: crn = GMAKSystem(species)
        sage: crn.add_complex(0, A + B, a * A + b * B)
        sage: crn.add_complex(1, C)
        sage: crn.add_complex(2, D, c * A + D)
        sage: crn.add_complex(3, A)
        sage: crn.add_complex(4, E)
        sage: crn.add_reaction(0, 1)
        sage: crn.add_reaction(1, 0, var("h"))
        sage: crn.add_reaction(1, 2)
        sage: crn.add_reaction(2, 0)
        sage: crn.add_reaction(3, 4)
        sage: crn.add_reaction(4, 3)
        sage: crn
        Reaction network with 5 complexes and 6 reactions.

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

        sage: crn.set_matrices()
        sage: crn.matrix_of_complexes_stoichiometric
        [1 1 0 0 0]
        [0 0 1 0 0]
        [0 0 0 1 0]
        [1 0 0 0 0]
        [0 0 0 0 1]
        sage: crn.matrix_of_complexes_kinetic_order
        [a b 0 0 0]
        [0 0 1 0 0]
        [c 0 0 1 0]
        [1 0 0 0 0]
        [0 0 0 0 1]
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
    def __init__(self, species: list) -> None:
        r"""
        A (chemical) reaction network with (generalized) mass-action kinetics.

        INPUT:

        - ``species`` -- a list of species.
        """
        self.species = species
        self.complexes = {}
        self.complexes_kinetic_order = {}
        self.graph = DiGraph()

        self.matrix_of_complexes_stoichiometric = None
        self.matrix_of_complexes_kinetic_order = None
        self.matrix_stoichiometric = None
        self.matrix_kinetic_order = None
        self.kernel_matrix_stoichiometric = None
        self.kernel_matrix_kinetic_order = None

    def _repr_(self) -> str:
        return f"Reaction network with {self.graph.num_verts()} complexes and {self.graph.num_edges()} reactions."

    def __copy__(self) -> GMAKSystem:
        new = GMAKSystem(self.species)
        new.complexes = self.complexes.copy()
        new.complexes_kinetic_order = self.complexes_kinetic_order.copy()
        new.graph = self.graph.copy()
        new.matrix_stoichiometric = self.matrix_stoichiometric
        new.matrix_kinetic_order = self.matrix_kinetic_order
        new.kernel_matrix_stoichiometric = self.kernel_matrix_stoichiometric
        new.kernel_matrix_kinetic_order = self.kernel_matrix_kinetic_order
        return new

    def __call__(self, **kwargs) -> GMAKSystem:
        new = copy(self)
        new.complexes = {i: complex(**kwargs) for i, complex in self.complexes.items()}
        new.complexes_kinetic_order = {i: complex(**kwargs) for i, complex in self.complexes_kinetic_order.items()}
        for attribute in [
            "matrix_of_complexes_stoichiometric",
            "matrix_of_complexes_kinetic_order",
            "matrix_stoichiometric",
            "matrix_kinetic_order",
            "kernel_matrix_stoichiometric",
            "kernel_matrix_kinetic_order",
        ]:
            attr = getattr(self, attribute)
            setattr(new, attribute, attr(**kwargs) if callable(attr) else attr)
        return new

    def add_complexes(self, complexes: list[tuple]) -> None:
        for complex in complexes:
            self.add_complex(*complex)

    def add_complex(self, i: int, complex, complex_kinetic_order=None) -> None:
        self.complexes[i] = complex
        self.complexes_kinetic_order[i] = complex if complex_kinetic_order is None else complex_kinetic_order
        self.graph.add_vertex(i)

    def remove_complex(self, i: int) -> None:
        self.complexes.pop(i)
        self.complexes_kinetic_order.pop(i)
        self.graph.delete_vertex(i)

    def add_reactions(self, reactions: list[tuple]) -> None:
        for reaction in reactions:
            self.add_reaction(*reaction)

    def add_reaction(self, start: int, end: int, label: str = None) -> None:
        for vertex in (start, end):
            if vertex not in self.complexes:
                self.add_complex(vertex, 0)
        if label is None:
            label = f"$k_{{{start}{end}}}$"
            # label = var(f"k_{start}{end}")
        self.graph.add_edge(start, end, label)

    def remove_reaction(self, start: int, end: int) -> None:
        self.graph.delete_edge(start, end)

    def reactions(self) -> list:
        return self.graph.edges()

    def add_species(self, species) -> None:
        # TODO do we need this
        self.species.append(species)

    def set_matrices(self) -> None:
        self.set_matrices_of_complexes()
        self.set_matrix_stoichiometric()
        self.set_matrix_kinetic_order()
        self.set_kernel_matrix_stoichiometric()
        self.set_kernel_matrix_kinetic_order()
        # TODO (re)computes matrices if changes have been detected?

    def set_matrices_of_complexes(self) -> None:
        self.matrix_of_complexes_stoichiometric = self._matrix_from_complexes(self.complexes)
        self.matrix_of_complexes_kinetic_order = self._matrix_from_complexes(self.complexes_kinetic_order)

    def _matrix_from_complexes(self, complexes: list):
        return matrix(
            [0 if complex == 0 else complex.coefficient(s) for s in self.species]
            for _, complex in sorted(complexes.items())
        )

    def set_matrix_stoichiometric(self):
        product = self.incidence_matrix().T * self.matrix_of_complexes_stoichiometric
        self.matrix_stoichiometric = product.matrix_from_rows(product.pivot_rows())

    def set_matrix_kinetic_order(self):
        product = self.incidence_matrix().T * self.matrix_of_complexes_kinetic_order
        self.matrix_kinetic_order = product.matrix_from_rows(product.pivot_rows())

    def set_kernel_matrix_stoichiometric(self):
        self.kernel_matrix_stoichiometric = kernel_matrix_using_elementary_vectors(self.matrix_stoichiometric)

    def set_kernel_matrix_kinetic_order(self):
        self.kernel_matrix_kinetic_order = kernel_matrix_using_elementary_vectors(self.matrix_kinetic_order)

    def incidence_matrix(self):
        r"""Return the incidence matrix of the graph."""
        return self.graph.incidence_matrix()

    def source_matrix(self):
        r"""Return the source matrix of the graph."""
        return matrix((1 if value == -1 else 0 for value in row) for row in self.incidence_matrix())

    def plot(self, kinetic_order: bool = True):
        return self.graph.plot(
            vertex_labels={i: self.vertex_label(i, kinetic_order=kinetic_order) for i in self.graph.vertices()},
            edge_labels=True,
            # edge_labels_background="transparent",
            vertex_colors="white",
            vertex_size=5000,
        )

    def vertex_label(self, i: int, kinetic_order: bool = False) -> str:
        if not kinetic_order or self.complexes[i] == self.complexes_kinetic_order[i]:
            return f"{self.complexes[i]}"
        return f"{self.complexes[i]}\n({self.complexes_kinetic_order[i]})"

    def edge_labels(self):
        # TODO remove?
        return self.graph.edge_labels()

    def deficiency_stoichiometric(self):
        r"""Return the stoichiometric deficiency."""
        return self.graph.num_verts() - self.graph.connected_components_number() - self.matrix_stoichiometric.rank()

    def deficiency_kinetic_order(self):
        r"""Return the kinetic-order deficiency."""
        return self.graph.num_verts() - self.graph.connected_components_number() - self.matrix_kinetic_order.rank()

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
