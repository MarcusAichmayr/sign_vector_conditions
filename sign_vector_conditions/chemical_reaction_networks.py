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
from sage.misc.latex import latex
from sage.calculus.var import var

from elementary_vectors import kernel_matrix_using_elementary_vectors

from .uniqueness import condition_uniqueness_minors
from .unique_existence import condition_faces, condition_nondegenerate
from .robustness import condition_closure_minors


class ReactionNetwork(SageObject):
    r"""
    Class for chemical reaction networks with generalized mass-action kinetics.

    A generalized mass-action system is represented by a (weighted) directed graph
    and stoichiometric and kinetic-order labels of the vertices.

    EXAMPLES:

    We define a chemical reaction network with generalized mass-action kinetics involving 5 complexes and 2 connected components::

        sage: from sign_vector_conditions import *
        sage: var('a, b')
        (a, b)
        sage: species = var('A, B, C')
        sage: rn = ReactionNetwork(species)
        sage: rn.add_complex(0, A + B, a * A + b * B)
        sage: rn.add_complex(1, C)
        sage: rn.add_reactions([(0, 1), (1, 0)])
        sage: rn
        Reaction network with 2 complexes and 2 reactions.
        sage: rn.plot()
        Graphics object consisting of 8 graphics primitives

    We describe the stoichiometric and kinetic-order subspaces using matrices::

        sage: rn.matrix_of_complexes_stoichiometric
        [1 0]
        [1 0]
        [0 1]
        sage: rn.matrix_of_complexes_kinetic_order
        [a 0]
        [b 0]
        [0 1]
        sage: rn.matrix_stoichiometric
        [-1  1]
        [-1  1]
        [ 1 -1]
        sage: rn.matrix_kinetic_order
        [-a  a]
        [-b  b]
        [ 1 -1]
        sage: rn.matrix_stoichiometric_as_kernel
        [1 0 1]
        [0 1 1]
        sage: rn.matrix_kinetic_order_as_kernel
        [1 0 a]
        [0 1 b]

    We check some conditions for our system::

        sage: rn.are_both_deficiencies_zero()
        True
        sage: rn.is_weakly_reversible()
        True
        sage: rn(a=2, b=1).has_robust_cbe()
        True
        sage: rn.has_robust_cbe()
        [{a > 0, b > 0}]
        sage: rn.has_at_most_one_cbe()
        [{a >= 0, b >= 0}]

    We extend our network by adding further complexes and reactions::

        sage: var('c')
        c
        sage: var('D, E')
        (D, E)
        sage: rn.add_species(D, E)
        sage: rn.add_complexes([(2, D, c * A + D), (3, A), (4, E)])
        sage: rn.add_reactions([(1, 2), (3, 4), (4, 3)])
        sage: rn.reactions
        [(0, 1), (1, 0), (1, 2), (3, 4), (4, 3)]
        sage: rn.rate_constants()
        (k_0_1, k_1_0, k_1_2, k_3_4, k_4_3)
        sage: rn
        Reaction network with 5 complexes and 5 reactions.
        sage: rn.plot()
        Graphics object consisting of 20 graphics primitives

    To make this system weakly reversible, we add another reaction.
    We specify a name for it::

        sage: rn.is_weakly_reversible()
        False
        sage: rn.add_reaction(2, 0)
        sage: rn
        Reaction network with 5 complexes and 6 reactions.
        sage: rn.is_weakly_reversible()
        True

    We compute the incidence and source matrices of the directed graph::

        sage: rn.incidence_matrix
        [-1  1  0  1  0  0]
        [ 1 -1 -1  0  0  0]
        [ 0  0  1 -1  0  0]
        [ 0  0  0  0 -1  1]
        [ 0  0  0  0  1 -1]
        sage: rn.source_matrix
        [1 0 0 0 0 0]
        [0 1 1 0 0 0]
        [0 0 0 1 0 0]
        [0 0 0 0 1 0]
        [0 0 0 0 0 1]

    We describe the stoichiometric and kinetic-order subspaces using matrices::

        sage: rn.matrix_of_complexes_stoichiometric
        [1 0 0 1 0]
        [1 0 0 0 0]
        [0 1 0 0 0]
        [0 0 1 0 0]
        [0 0 0 0 1]
        sage: rn.matrix_of_complexes_kinetic_order
        [a 0 c 1 0]
        [b 0 0 0 0]
        [0 1 0 0 0]
        [0 0 1 0 0]
        [0 0 0 0 1]
        sage: rn.matrix_stoichiometric
        [-1  1  0  1 -1  1]
        [-1  1  0  1  0  0]
        [ 1 -1 -1  0  0  0]
        [ 0  0  1 -1  0  0]
        [ 0  0  0  0  1 -1]
        sage: rn.matrix_kinetic_order
        [   -a     a     c a - c    -1     1]
        [   -b     b     0     b     0     0]
        [    1    -1    -1     0     0     0]
        [    0     0     1    -1     0     0]
        [    0     0     0     0     1    -1]
        sage: rn.matrix_stoichiometric_reduced
        [-1  0 -1]
        [-1  0  0]
        [ 1 -1  0]
        [ 0  1  0]
        [ 0  0  1]
        sage: rn.matrix_kinetic_order_reduced
        [-a  c -1]
        [-b  0  0]
        [ 1 -1  0]
        [ 0  1  0]
        [ 0  0  1]
        sage: rn.matrix_stoichiometric_as_kernel
        [1 0 1 1 1]
        [0 1 1 1 0]
        sage: rn.matrix_kinetic_order_as_kernel
        [    1     0     a a - c     1]
        [    0     1     b     b     0]

    We check some conditions for our system::

        sage: rn.are_both_deficiencies_zero()
        True
        sage: rn.is_weakly_reversible()
        True
        sage: rn(a=2, b=1, c=1).has_robust_cbe()
        True
        sage: rn.has_robust_cbe() # random order
        [{a > 0, a - c > 0, b > 0}]
        sage: rn.has_at_most_one_cbe() # random order
        [{a >= 0, a - c >= 0, b >= 0}]
        sage: rn.has_exactly_one_cbe()
        Traceback (most recent call last):
        ...
        ValueError: Method does not support variables!
        sage: rn(a=2, b=1, c=1).has_exactly_one_cbe()
        True

    We remove one component and a reaction of our system::

        sage: rn.remove_complex(3)
        sage: rn.remove_complex(4)
        sage: rn.remove_reaction(1, 0)
        sage: rn.remove_species(E)
        sage: rn
        Reaction network with 3 complexes and 3 reactions.
        sage: rn.plot()
        Graphics object consisting of 10 graphics primitives
        sage: rn.is_weakly_reversible()
        True
        sage: rn.has_at_most_one_cbe() # random order
        [{a >= 0, a - c >= 0, b >= 0}]

    If a species is not specified, it will be ignored::

        sage: rn.add_complex(4, E)
        sage: rn.add_reaction(1, 4)
        sage: rn
        Reaction network with 4 complexes and 4 reactions.
        sage: rn.species
        [A, B, C, D]
        sage: rn.plot()
        Graphics object consisting of 13 graphics primitives

    To overwrite complexes we add them as usual::

        sage: rn.add_species(E)
        sage: rn.add_complex(4, E)
        sage: rn
        Reaction network with 4 complexes and 4 reactions.
        sage: rn.species
        [A, B, C, D, E]
        sage: rn.plot()
        Graphics object consisting of 13 graphics primitives

    We can set the rate constant variable to a different name::

        sage: rn.set_rate_constant_variable("t")
        sage: rn.rate_constants()
        (t_0_1, t_1_2, t_1_4, t_2_0)
        sage: rn.plot()
        Graphics object consisting of 13 graphics primitives
    """
    def __init__(self, species: list) -> None:
        r"""
        A (chemical) reaction network with (generalized) mass-action kinetics.

        INPUT:

        - ``species`` -- a list of species.
        """
        self.species = species if isinstance(species, list) else list(species)
        self.complexes = {}
        self.complexes_kinetic_order = {}
        self._rate_constant_variable = "k"
        self.graph = DiGraph()

        self._matrix_of_complexes_stoichiometric = None
        self._matrix_of_complexes_kinetic_order = None
        self._matrix_stoichiometric = None
        self._matrix_kinetic_order = None
        self._matrix_stoichiometric_reduced = None
        self._matrix_kinetic_order_reduced = None
        self._kernel_matrix_stoichiometric = None # TODO remove member
        self._kernel_matrix_kinetic_order = None # TODO remove member

        self._update_needed = True

    def _repr_(self) -> str:
        return f"Reaction network with {self.graph.num_verts()} complexes and {self.graph.num_edges()} reactions."

    def __copy__(self) -> ReactionNetwork:
        new = ReactionNetwork(self.species)
        new.complexes = self.complexes.copy()
        new.complexes_kinetic_order = self.complexes_kinetic_order.copy()
        new.graph = self.graph.copy()
        new._matrix_stoichiometric = self._matrix_stoichiometric
        new._matrix_kinetic_order = self._matrix_kinetic_order
        new._kernel_matrix_stoichiometric = self._kernel_matrix_stoichiometric
        new._kernel_matrix_kinetic_order = self._kernel_matrix_kinetic_order

        new._update_needed = False
        return new

    def __call__(self, **kwargs) -> ReactionNetwork:
        new = copy(self)
        new.complexes = {i: complex(**kwargs) for i, complex in self.complexes.items()}
        new.complexes_kinetic_order = {i: complex(**kwargs) for i, complex in self.complexes_kinetic_order.items()}
        for attribute in [
            "_matrix_of_complexes_stoichiometric",
            "_matrix_of_complexes_kinetic_order",
            "_matrix_stoichiometric",
            "_matrix_kinetic_order",
            "_matrix_stoichiometric_reduced",
            "_matrix_kinetic_order_reduced",
            "_kernel_matrix_stoichiometric",
            "_kernel_matrix_kinetic_order",
        ]:
            attr = getattr(self, attribute)
            setattr(new, attribute, attr(**kwargs) if callable(attr) else attr)
        return new

    def add_complexes(self, complexes: list[tuple]) -> None:
        r"""Add complexes to system."""
        for complex in complexes:
            self.add_complex(*complex)

    def add_complex(self, i: int, complex, complex_kinetic_order=None) -> None:
        r"""Add complex to system."""
        self.complexes[i] = complex
        self.complexes_kinetic_order[i] = complex if complex_kinetic_order is None else complex_kinetic_order
        self.graph.add_vertex(i)
        self._update_needed = True

    def remove_complex(self, i: int) -> None:
        r"""Remove complex from system."""
        self.complexes.pop(i)
        self.complexes_kinetic_order.pop(i)
        self.graph.delete_vertex(i)
        self._update_needed = True

    def add_reactions(self, reactions: list[tuple]) -> None:
        r"""Add reactions to system."""
        for reaction in reactions:
            self.add_reaction(*reaction)

    def add_reaction(self, start: int, end: int) -> None:
        r"""Add reaction to system."""
        for vertex in (start, end):
            if vertex not in self.complexes:
                self.add_complex(vertex, 0)
        self.graph.add_edge(start, end, label=f"${self._rate_constant_variable}_{{{start}, {end}}}$")
        self._update_needed = True

    def remove_reaction(self, start: int, end: int) -> None:
        r"""Remove reaction from system."""
        self.graph.delete_edge(start, end)
        self._update_needed = True

    @property
    def reactions(self) -> list:
        r"""Return reactions."""
        return [(start, end) for start, end, _ in self.graph.edges()]

    def set_rate_constant_variable(self, variable: str) -> None:
        r"""Set rate constant variable."""
        self._rate_constant_variable = variable
        for edge in self.graph.edges():
            self.graph.set_edge_label(edge[0], edge[1], f"${variable}_{{{edge[0]}, {edge[1]}}}$")

    def rate_constants(self) -> tuple:
        r"""Return rate constants."""
        return tuple(var(f"{self._rate_constant_variable}_{start}_{end}") for start, end in self.reactions)

    def add_species(self, *species) -> None:
        r"""Add one or more species."""
        for s in species:
            self.species.append(s)
        self._update_needed = True

    def remove_species(self, *species) -> None:
        r"""Remove one or more species."""
        for s in species:
            self.species.remove(s)
        self._update_needed = True

    def _update_matrices(self) -> None:
        r"""Set stoichiometric and kinetic-order matrices."""
        if not self._update_needed:
            return
        self._matrix_of_complexes_stoichiometric = self._matrix_from_complexes(self.complexes)
        self._matrix_of_complexes_kinetic_order = self._matrix_from_complexes(self.complexes_kinetic_order)

        self._matrix_stoichiometric = self.incidence_matrix.T * self._matrix_of_complexes_stoichiometric
        self._matrix_kinetic_order = self.incidence_matrix.T * self._matrix_of_complexes_kinetic_order

        self._matrix_stoichiometric_reduced = self._matrix_stoichiometric.matrix_from_rows(self._matrix_stoichiometric.pivot_rows())
        self._matrix_kinetic_order_reduced = self._matrix_kinetic_order.matrix_from_rows(self._matrix_kinetic_order.pivot_rows())

        try:
            self._kernel_matrix_stoichiometric = kernel_matrix_using_elementary_vectors(self._matrix_stoichiometric_reduced)
            self._kernel_matrix_kinetic_order = kernel_matrix_using_elementary_vectors(self._matrix_kinetic_order_reduced)
        except ValueError:
            print("TODO what should we do, when no zero minor?")
        self._update_needed = False

    def _matrix_from_complexes(self, complexes: list):
        return matrix(
            [0 if complex == 0 else complex.coefficient(s) for s in self.species]
            for _, complex in sorted(complexes.items())
        )

    def _get_matrix(self, matrix_name):
        self._update_matrices()
        return getattr(self, matrix_name)

    @property
    def matrix_of_complexes_stoichiometric(self):
        r"""Return the matrix that decodes the stoichiometric complexes of the reaction network."""
        return self._get_matrix('_matrix_of_complexes_stoichiometric').T

    @property
    def matrix_of_complexes_kinetic_order(self):
        r"""Return the matrix that decodes the kinetic-order complexes of the reaction network."""
        return self._get_matrix('_matrix_of_complexes_kinetic_order').T

    @property
    def matrix_stoichiometric(self):
        r"""Return the stoichiometric matrix."""
        return self._get_matrix('_matrix_stoichiometric').T

    @property
    def matrix_kinetic_order(self):
        r"""Return the kinetic-order matrix."""
        return self._get_matrix('_matrix_kinetic_order').T

    @property
    def matrix_stoichiometric_reduced(self):
        r"""Return the reduced stoichiometric matrix."""
        return self._get_matrix('_matrix_stoichiometric_reduced').T

    @property
    def matrix_kinetic_order_reduced(self):
        r"""Return the reduced kinetic-order matrix."""
        return self._get_matrix('_matrix_kinetic_order_reduced').T

    @property
    def matrix_stoichiometric_as_kernel(self):
        r"""Return the kernel matrix of the stoichiometric matrix."""
        return self._get_matrix('_kernel_matrix_stoichiometric')

    @property
    def matrix_kinetic_order_as_kernel(self):
        r"""Return the kernel matrix of the kinetic-order matrix."""
        return self._get_matrix('_kernel_matrix_kinetic_order')

    @property
    def incidence_matrix(self):
        r"""Return the incidence matrix of the graph."""
        return self.graph.incidence_matrix()

    @property
    def source_matrix(self):
        r"""Return the source matrix of the graph."""
        return matrix((1 if value == -1 else 0 for value in row) for row in self.incidence_matrix)

    def plot(
            self,
            kinetic_order: bool = True,
            layout="circular",
            edge_labels=True,
            vertex_colors="white",
            vertex_size=5000,
            **kwargs
        ):
        r"""Plot the reaction network."""
        return self.graph.plot(
            vertex_labels={i: self._vertex_label(i, kinetic_order=kinetic_order) for i in self.graph.vertices()},
            layout=layout,
            edge_labels=edge_labels,
            # edge_labels_background="transparent",
            vertex_colors=vertex_colors,
            vertex_size=vertex_size,
            **kwargs
        )

    def _vertex_label(self, i: int, kinetic_order: bool = False) -> str:
        if not kinetic_order or self.complexes[i] == self.complexes_kinetic_order[i]:
            return f"${latex(self.complexes[i])}$"
        return f"${latex(self.complexes[i])}$\n$({latex(self.complexes_kinetic_order[i])})$"

    def deficiency_stoichiometric(self):
        r"""Return the stoichiometric deficiency."""
        self._update_matrices()
        return self.graph.num_verts() - self.graph.connected_components_number() - self._matrix_stoichiometric.rank()

    def deficiency_kinetic_order(self):
        r"""Return the kinetic-order deficiency."""
        self._update_matrices()
        return self.graph.num_verts() - self.graph.connected_components_number() - self._matrix_kinetic_order.rank()

    def are_both_deficiencies_zero(self) -> bool:
        r"""Return whether both deficiencies are zero."""
        return self.deficiency_stoichiometric() == self.deficiency_kinetic_order() == 0

    def is_weakly_reversible(self) -> bool:
        r"""Return whether each component of the system is strongly connected."""
        return all(g.is_strongly_connected() for g in self.graph.connected_components_subgraphs())

    def _check_network_conditions(self):
        r"""Perform common network checks for uniqueness and existence of CBE."""
        self._update_matrices()
        if self.deficiency_stoichiometric() != 0:
            raise ValueError(f"Stoichiometric deficiency should be zero and not {self.deficiency_stoichiometric()}!")
        if self.deficiency_kinetic_order() != 0:
            raise ValueError(f"Kinetic-order deficiency should be zero and not {self.deficiency_kinetic_order()}!")
        if not self.is_weakly_reversible():
            raise ValueError("Network is not weakly reversible!")

    def has_robust_cbe(self):
        r"""Check whether there is a unique positive CBE with regards to small perturbations."""
        self._check_network_conditions()
        return condition_closure_minors(self._kernel_matrix_stoichiometric, self._kernel_matrix_kinetic_order)

    def has_at_most_one_cbe(self):
        r"""Check whether there is at most one positive CBE."""
        self._check_network_conditions()
        return condition_uniqueness_minors(self._kernel_matrix_stoichiometric, self._kernel_matrix_kinetic_order)

    def _condition_faces(self) -> bool:
        r"""Check whether the system satisfies the face condition for existence of a unique positive CBE."""
        self._check_network_conditions()
        return condition_faces(self._kernel_matrix_stoichiometric, self._kernel_matrix_kinetic_order)

    def _are_subspaces_nondegenerate(self) -> bool:
        r"""Check whether the system satisfies the nondegenerate condition for existence of a unique positive CBE."""
        self._check_network_conditions()
        return condition_nondegenerate(self._kernel_matrix_stoichiometric, self._kernel_matrix_kinetic_order)

    def has_exactly_one_cbe(self) -> bool:
        r"""Check whether there is exactly one positive CBE."""
        self._check_network_conditions()
        at_most_one = self.has_at_most_one_cbe()
        if at_most_one not in [True, False]:
            raise ValueError("Method does not support variables!")
        return (
            at_most_one
            and self._condition_faces()
            and self._are_subspaces_nondegenerate()
        )
