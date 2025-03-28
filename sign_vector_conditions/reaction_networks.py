r"""
Module for setting up reaction networks.

This module provides tools for defining species, complexes, and reaction networks
with (generalized) mass-action kinetics. It includes functionality for analyzing
reaction networks, such as checking weak reversibility, computing deficiencies,
and verifying conditions for unique positive complex-balanced equilibria (CBE).

Key Classes and Functions:
- :class:`ReactionNetwork`: Represents a chemical reaction network.
- :func:`species`: Utility function for defining species.

For detailed examples and usage, see :class:`ReactionNetwork`.

EXAMPLES:

We define species for a reaction network::

    sage: from sign_vector_conditions import *
    sage: A, B, C = species("A, B, C")
    sage: rn = ReactionNetwork()
    sage: rn.add_complex(0, A + B)
    sage: rn.add_complex(1, C)
    sage: rn.add_reactions([(0, 1), (1, 0)])
    sage: rn.plot()
    Graphics object consisting of 6 graphics primitives

We can analyze the reaction network::

    sage: rn.is_weakly_reversible()
    True
    sage: rn.deficiency_stoichiometric
    0
"""

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

import inspect

from copy import copy
from typing import NamedTuple, Tuple, List, Dict, Union
from sage.calculus.var import var
from sage.graphs.digraph import DiGraph
from sage.structure.sage_object import SageObject
from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix
from sage.modules.free_module_element import vector
from sage.misc.latex import latex
from sage.misc.misc_c import prod

from elementary_vectors import kernel_matrix_using_elementary_vectors

from .uniqueness import condition_uniqueness_minors
from .unique_existence import condition_faces, condition_nondegenerate
from .robustness import condition_closure_minors


def species(names: str):
    r"""
    Define species from a string of names.

    See :class:`Complex` for operations and more details.

    EXAMPLES::

        sage: from sign_vector_conditions import *
        sage: species("A")
        A
        sage: A
        A
        sage: species("A, B, C")
        (A, B, C)
        sage: A
        A
        sage: B
        B
        sage: C
        C
        sage: species("A,B,C")
        (A, B, C)
        sage: species("A B C")
        (A, B, C)
    """
    names = names.strip()
    if "," in names:
        names_list = [s.strip() for s in names.split(",") if s.strip()]
    elif " " in names:
        names_list = [s.strip() for s in names.split()]
    else:
        names_list = [names] if names else []

    caller_globals = inspect.currentframe().f_back.f_globals

    def define_species_globally(name: str) -> Complex:
        complex = Complex({Species(name): 1})
        caller_globals[name] = complex
        return complex

    if len(names_list) == 1:
        return define_species_globally(names_list[0])
    return tuple(define_species_globally(name) for name in names_list)


class Species(NamedTuple):
    r"""
    Auxiliary class for species.

    To compute with species, use the :func:`species` function.
    """
    name: str
    def __str__(self) -> str:
        return self.name

    def _latex_(self) -> str:
        return self.name


class Complex(SageObject):
    r"""
    A complex involving species.

    EXAMPLES:

    First, we define some species::

        sage: from sign_vector_conditions import *
        sage: species("A, B, C")
        (A, B, C)

    Usual operations like addition and multiplication are supported::

        sage: A + B
        A + B
        sage: 2 * A + 3 * B
        2*A + 3*B

    Symbolic expressions are also supported::

        sage: var("a")
        a
        sage: 2 * a * A
        2*a*A
        sage: (2 + a) * A
        (a + 2)*A
        sage: a * (A + B)
        a*A + a*B

    TESTS::

        sage: A * B
        Traceback (most recent call last):
        ...
        TypeError: Cannot multiply species by species.
        sage: B + A
        A + B
        sage: A * 2
        2*A
        sage: 0 * A
        0
        sage: 1 * A
        A
        sage: A + 0
        A
        sage: A + 1
        Traceback (most recent call last):
        ...
        TypeError: Cannot add <class 'sage.rings.integer.Integer'> to species.
        sage: -A
        -A
        sage: (-a - 1) * A
        (-a - 1)*A
        sage: A - B
        A - B
        sage: A - 2 * B
        A - 2*B
        sage: (2 * A + 3 * B)[A]
        2
        sage: A in 2 * A + B
        True
        sage: A in 3 * B
        False
        sage: s = a * A - B
        sage: s
        a*A - B
        sage: s(a=0)
        -B
        sage: s(a=1)
        A - B
        sage: (a * A)(a=0)
        0
        sage: A(a=1)
        A

    We test the latex representation::

        sage: species("A, B")
        (A, B)
        sage: 2 * A + 3 * B
        2*A + 3*B
        sage: (2 * A + 3 * B)._latex_()
        '2 \\, A + 3 \\, B'
        sage: (A - B)._latex_()
        'A - B'
        sage: (A + B)._latex_()
        'A + B'
        sage: (2 * A - 3 * B)._latex_()
        '2 \\, A - 3 \\, B'
        sage: (A)._latex_()
        'A'
        sage: (0 * A)._latex_()
        '0'
    """
    def __init__(self, species_dict: Dict[Species, Union[int, float, var]]) -> None:
        self.species_dict = {}
        for key, value in species_dict.items():
            if not isinstance(key, Species):
                raise TypeError(f"Key {key} is not a Species")
            if value == 0:
                continue
            self.species_dict[key] = value

    def __str__(self) -> str:
        return self._repr_()

    def __copy__(self) -> Complex:
        return Complex(self.species_dict)

    def __call__(self, **kwargs) -> Complex:
        return Complex({
            key: (value(**kwargs) if callable(value) else value)
            for key, value in self.species_dict.items()
        })

    def __contains__(self, species: Union[Species, Complex]) -> bool:
        if isinstance(species, Species):
            return species in self.species_dict
        if isinstance(species, Complex):
            return species._to_species() in self.species_dict
        return False

    def __getitem__(self, species: Union[Species, Complex]) -> Union[int, float]:
        if isinstance(species, Species):
            return self.species_dict.get(species, 0)
        if isinstance(species, Complex):
            return self.species_dict.get(species._to_species(), 0)
        raise TypeError(f"Cannot get {type(species)} from species.")

    def __add__(self, other) -> Complex:
        if other == 0:
            return copy(self)
        if not isinstance(other, Complex):
            raise TypeError(f"Cannot add {type(other)} to species.")
        species_dict = self.species_dict.copy()
        for key, value in other.species_dict.items():
            if key in species_dict:
                species_dict[key] += value
            else:
                species_dict[key] = value
        return Complex(species_dict)

    def __sub__(self, other) -> Complex:
        return self + (-other)

    def __mul__(self, other) -> Complex:
        if isinstance(other, Complex):
            raise TypeError("Cannot multiply species by species.")
        species_dict = {key: value * other for key, value in self.species_dict.items()}
        return Complex(species_dict)

    def __rmul__(self, other) -> Complex:
        return self.__mul__(other)

    def __neg__(self) -> Complex:
        return -1 * self

    def __pos__(self) -> Complex:
        return copy(self)

    def _repr_(self) -> str:
        return self._format_repr_(self._repr_coefficient)

    def _latex_(self) -> str:
        return self._format_repr_(self._latex_coefficient)

    def _format_repr_(self, coefficient_function) -> str:
        if not self.species_dict:
            return "0"
        terms = []
        for key, _ in sorted(self.species_dict.items()):
            summand = coefficient_function(key)
            if not terms:
                terms.append(summand)
            elif str(summand)[0] == "-":
                terms.append(f"- {summand[1:]}")
            else:
                terms.append(f"+ {summand}")
        return " ".join(terms)

    def _repr_coefficient(self, key: Species) -> str:
        return self._format_coefficient(key, str)

    def _latex_coefficient(self, key: Species) -> str:
        return self._format_coefficient(key, latex)

    def _format_coefficient(self, key: Species, formatter) -> str:
        value = self.species_dict[key]
        formatted_key = formatter(key)
        formatted_value = formatter(value)

        if value == 1:
            return formatted_key
        if value == -1:
            return f"-{formatted_key}"
        if "+" in str(value) or " - " in str(value):
            return f"({formatted_value})*{formatted_key}" if formatter == str else rf"({formatted_value}) \, {formatted_key}"
        return f"{formatted_value}*{formatted_key}" if formatter == str else rf"{formatted_value} \, {formatted_key}"

    def _to_species(self) -> Species:
        if len(self.species_dict) != 1:
            raise ValueError("Complex must contain exactly one species.")
        return next(iter(self.species_dict.keys()))

    def involved_species(self) -> set[Species]:
        r"""Return the species involved in the complex."""
        return set(self.species_dict.keys())

    @staticmethod
    def from_species(species: Species) -> Complex:
        r"""Return a complex from a species."""
        return Complex({species: 1})


class ReactionNetwork(SageObject):
    r"""
    A reaction network with (generalized) mass-action kinetics.

    This class represents a reaction network, where complexes are connected
    by directed reactions. It supports generalized mass-action kinetics, allowing
    for symbolic rate constants and kinetic orders.

    The `ReactionNetwork` class provides tools for:
    - Adding and removing complexes and reactions.
    - Computing stoichiometric and kinetic-order matrices.
    - Analyzing network properties, such as weak reversibility and deficiencies.
    - Checking conditions for unique positive complex-balanced equilibria (CBE).
    - Visualizing the reaction network as a directed graph.

    Key Attributes:
    - `graph`: The directed graph representing the reaction network.
    - `complexes_stoichiometric`: A dictionary mapping complex indices to stoichiometric complexes.
    - `complexes_kinetic_order`: A dictionary mapping complex indices to kinetic-order complexes.
    - `species`: A tuple of all species involved in the network.

    Key Methods:
    - `add_complex`, `remove_complex`: Add or remove complexes from the network.
    - `add_reaction`, `remove_reaction`: Add or remove reactions between complexes.
    - `plot`: Visualize the reaction network as a directed graph.

    EXAMPLES:

    We define a reaction network with two complexes involving variables in the kinetic orders::

        sage: from sign_vector_conditions import *
        sage: var("a, b")
        (a, b)
        sage: species("A, B, C")
        (A, B, C)
        sage: rn = ReactionNetwork()
        sage: rn.add_complex(0, A + B, a * A + b * B)
        sage: rn.add_complex(1, C)
        sage: rn.add_reactions([(0, 1), (1, 0)])
        sage: rn
        Reaction network with 2 complexes and 2 reactions.
        sage: rn.plot()
        Graphics object consisting of 6 graphics primitives

    We describe the reaction network using matrices::

        sage: rn.matrix_of_complexes_stoichiometric
        [1 0]
        [1 0]
        [0 1]
        sage: rn.matrix_of_complexes_kinetic_order
        [a 0]
        [b 0]
        [0 1]

    The stoichiometric and kinetic-order matrices are given by::

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

        sage: var("c")
        c
        sage: species("D, E")
        (D, E)
        sage: rn.add_complexes([(2, D, c * A + D), (3, A), (4, E)])
        sage: rn.add_reactions([(1, 2), (3, 4), (4, 3)])
        sage: rn
        Reaction network with 5 complexes and 5 reactions.
        sage: rn.plot()
        Graphics object consisting of 15 graphics primitives

    The network involves the following species::

        sage: rn.species
        (A, B, C, D, E)

    To make this system weakly reversible, we add another reaction::

        sage: rn.is_weakly_reversible()
        False
        sage: rn.add_reaction(2, 0)
        sage: rn.is_weakly_reversible()
        True

    Now, our network consists of 6 reactions::

        sage: rn.reactions
        [(0, 1), (1, 0), (1, 2), (2, 0), (3, 4), (4, 3)]

    The corresponding rate constants are::

        sage: rn.rate_constants()
        (k_0_1, k_1_0, k_1_2, k_2_0, k_3_4, k_4_3)

    We compute the incidence and source matrices of the directed graph::

        sage: rn.incidence_matrix()
        [-1  1  0  1  0  0]
        [ 1 -1 -1  0  0  0]
        [ 0  0  1 -1  0  0]
        [ 0  0  0  0 -1  1]
        [ 0  0  0  0  1 -1]
        sage: rn.source_matrix()
        [1 0 0 0 0 0]
        [0 1 1 0 0 0]
        [0 0 0 1 0 0]
        [0 0 0 0 1 0]
        [0 0 0 0 0 1]

    The Laplacian matrix involving the rate constants is given by::

        sage: rn.laplacian_matrix()
        [        -k_0_1          k_1_0          k_2_0              0              0]
        [         k_0_1 -k_1_0 - k_1_2              0              0              0]
        [             0          k_1_2         -k_2_0              0              0]
        [             0              0              0         -k_3_4          k_4_3]
        [             0              0              0          k_3_4         -k_4_3]
        sage: rn.differential_equation()
        (-k_0_1*x_0^a*x_2^c*x_3 + k_1_0*x_0^b + k_2_0*x_1 - k_3_4*x_2 + k_4_3*x_4,
         -k_0_1*x_0^a*x_2^c*x_3 + k_1_0*x_0^b + k_2_0*x_1,
         k_0_1*x_0^a*x_2^c*x_3 - (k_1_0 + k_1_2)*x_0^b,
         k_1_2*x_0^b - k_2_0*x_1,
         k_3_4*x_2 - k_4_3*x_4)

    The network is described by the following matrices::

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
        sage: rn.matrix_stoichiometric_as_kernel
        [1 0 1 1 1]
        [0 1 1 1 0]
        sage: rn.matrix_kinetic_order_as_kernel
        [    1     0     a a - c     1]
        [    0     1     b     b     0]

    We check some conditions for our system::

        sage: rn.deficiency_stoichiometric
        0
        sage: rn.deficiency_kinetic_order
        0
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
        sage: rn
        Reaction network with 3 complexes and 3 reactions.

    We can change the names of the rate constants::

        sage: rn.set_rate_constant_variable(var("tau"))
        sage: rn.rate_constants()
        (tau_0_1, tau_1_2, tau_2_0)
        sage: rn.plot(edge_labels=True)
        Graphics object consisting of 10 graphics primitives

    ::

        sage: A, B, C = species("H_2, O_2, H_2O")
        sage: var('a')
        a
        sage: rn = ReactionNetwork()
        sage: rn.add_complex(0, 2 * A + B, 2 * a * A + a * B)
        sage: rn.add_complex(1, 2 * C)
        sage: rn.species
        (H_2, H_2O, O_2)
        sage: rn.add_reactions([(0, 1), (1, 0)])
        sage: rn.plot()
        Graphics object consisting of 6 graphics primitives

    ::

        sage: fox, rabbit = species("Fox, Rabbit")
        sage: rn = ReactionNetwork()
        sage: rn.add_complex(0, rabbit + fox)
        sage: rn.add_complex(1, 2 * fox)
        sage: rn.add_reactions([(0, 1), (1, 0)])
        sage: rn.plot()
        Graphics object consisting of 6 graphics primitives
    """
    def __init__(self) -> None:
        r"""
        A (chemical) reaction network with (generalized) mass-action kinetics.

        INPUT:

        - ``species`` -- a list of species.
        """
        self._update_needed: bool = True
        self._rate_constant_variable: var = var("k")

        self.graph: DiGraph = DiGraph()
        self.complexes_stoichiometric: Dict[int, Complex] = {}
        self.complexes_kinetic_order: Dict[int, Complex] = {}

        self._species: List[Species] = []

        self._matrix_of_complexes_stoichiometric = None
        self._matrix_of_complexes_kinetic_order = None
        self._matrix_stoichiometric = None
        self._matrix_kinetic_order = None
        self._matrix_stoichiometric_reduced = None
        self._matrix_kinetic_order_reduced = None

        self._deficiency_stoichiometric: int = 0
        self._deficiency_kinetic_order: int = 0

    def _repr_(self) -> str:
        return f"Reaction network with {self.graph.num_verts()} complexes and {self.graph.num_edges()} reactions."

    def __copy__(self) -> ReactionNetwork:
        new = ReactionNetwork()
        for attribute in vars(self):
            setattr(new, attribute, copy(getattr(self, attribute)))
        return new

    def __call__(self, **kwargs) -> ReactionNetwork:
        new = copy(self)
        new.complexes_stoichiometric = {i: complex(**kwargs) for i, complex in self.complexes_stoichiometric.items()}
        new.complexes_kinetic_order = {i: complex(**kwargs) for i, complex in self.complexes_kinetic_order.items()}
        new._update_needed = True
        return new

    def add_complexes(self, complexes: List[Tuple[int, Complex, Union[Complex, None]]]) -> None:
        r"""Add complexes to system."""
        for element in complexes:
            self.add_complex(*element)

    def add_complex(self, i: int, complex_stoichiometric: Complex, complex_kinetic_order: Union[Complex, None] = None) -> None:
        r"""Add complex to system."""
        self.complexes_stoichiometric[i] = complex_stoichiometric
        self.complexes_kinetic_order[i] = complex_stoichiometric if complex_kinetic_order is None else complex_kinetic_order
        self.graph.add_vertex(i)
        self._update_needed = True

    def remove_complex(self, i: int) -> None:
        r"""Remove complex from system."""
        self.complexes_stoichiometric.pop(i)
        self.complexes_kinetic_order.pop(i)
        self.graph.delete_vertex(i)
        self._update_needed = True

    def add_reactions(self, reactions: List[Tuple[int, int]]) -> None:
        r"""Add reactions to system."""
        for reaction in reactions:
            self.add_reaction(*reaction)

    def add_reaction(self, start: int, end: int) -> None:
        r"""Add reaction to system."""
        for vertex in (start, end):
            if vertex not in self.complexes_stoichiometric:
                self.add_complex(vertex, 0)
        self.graph.add_edge(start, end)
        self._update_needed = True

    def remove_reaction(self, start: int, end: int) -> None:
        r"""Remove reaction from system."""
        self.graph.delete_edge(start, end)
        self._update_needed = True

    @property
    def reactions(self) -> List[Tuple[int, int]]:
        r"""Return reactions."""
        return [(start, end) for start, end, _ in self.graph.edges()]

    @property
    def species(self) -> Tuple[Complex, ...]:
        r"""Return the species of the reaction network as a tuple of complexes."""
        self._update()
        return tuple(Complex.from_species(s) for s in self._species)

    @property
    def matrix_of_complexes_stoichiometric(self) -> matrix:
        r"""Return the matrix that decodes the stoichiometric complexes of the reaction network."""
        return self._get("_matrix_of_complexes_stoichiometric").T

    @property
    def matrix_of_complexes_kinetic_order(self) -> matrix:
        r"""Return the matrix that decodes the kinetic-order complexes of the reaction network."""
        return self._get("_matrix_of_complexes_kinetic_order").T

    @property
    def matrix_stoichiometric(self) -> matrix:
        r"""Return the stoichiometric matrix where the columns correspond to the reactions."""
        return self._get("_matrix_stoichiometric").T

    @property
    def matrix_kinetic_order(self) -> matrix:
        r"""Return the kinetic-order matrix where the columns correspond to the reactions."""
        return self._get("_matrix_kinetic_order").T

    @property
    def matrix_stoichiometric_as_kernel(self) -> matrix:
        r"""Return the kernel matrix of the stoichiometric matrix."""
        self._update()
        return kernel_matrix_using_elementary_vectors(self._matrix_stoichiometric_reduced)

    @property
    def matrix_kinetic_order_as_kernel(self) -> matrix:
        r"""Return the kernel matrix of the kinetic-order matrix."""
        self._update()
        return kernel_matrix_using_elementary_vectors(self._matrix_kinetic_order_reduced)

    @property
    def deficiency_stoichiometric(self) -> int:
        r"""Return the stoichiometric deficiency."""
        return self._get("_deficiency_stoichiometric")

    @property
    def deficiency_kinetic_order(self) -> int:
        r"""Return the kinetic-order deficiency."""
        return self._get("_deficiency_kinetic_order")

    def incidence_matrix(self, **kwargs) -> matrix:
        r"""Return the incidence matrix of the graph."""
        return self.graph.incidence_matrix(**kwargs)

    def source_matrix(self, **kwargs) -> matrix:
        r"""Return the source matrix of the graph."""
        return matrix((1 if value == -1 else 0 for value in row) for row in self.incidence_matrix(**kwargs))

    def laplacian_matrix(self) -> matrix:
        r"""Return the Laplacian matrix of the graph."""
        return self.incidence_matrix() * diagonal_matrix(self.rate_constants()) * self.source_matrix().T

    def differential_equation(self) -> vector:
        r"""Return the differential equation of the system."""
        self._update()
        x = vector(var(f"x_{i}") for i in range(len(self.complexes_stoichiometric)))
        return (
            self._matrix_of_complexes_stoichiometric.T * self.laplacian_matrix() * vector(
                prod(xi ** yi for xi, yi in zip(x, y))
                for y in self._matrix_of_complexes_kinetic_order.columns()
            )
        )

    def rate_constants(self) -> Tuple[var, ...]:
        r"""Return rate constants."""
        return tuple(self._rate_constant(*edge) for edge in self.reactions)

    def _rate_constant(self, start: int, end: int) -> var:
        return var(
            f"{self._rate_constant_variable}_{start}_{end}",
            latex_name=f"{latex(self._rate_constant_variable)}_{{{start}, {end}}}"
        )

    def set_rate_constant_variable(self, variable: var) -> None:
        r"""
        Set rate constant variable.
        This method allows you to set a custom variable for the rate constants.

        EXAMPLES::

            sage: from sign_vector_conditions import *
            sage: species("A, B, C")
            (A, B, C)
            sage: rn = ReactionNetwork()
            sage: rn.add_complexes([(0, A + B), (1, C)])
            sage: rn.add_reactions([(0, 1), (1, 0)])
            sage: rn.rate_constants()
            (k_0_1, k_1_0)
            sage: rn.plot(edge_labels=True)
            Graphics object consisting of 8 graphics primitives

        You can also use a variable with a LaTeX name::

            sage: rn.set_rate_constant_variable(var("tau"))
            sage: rn.rate_constants()
            (tau_0_1, tau_1_0)
            sage: rn.plot(edge_labels=True)
            Graphics object consisting of 8 graphics primitives
            sage: var("k", latex_name=r"\kappa")
            k
            sage: rn.set_rate_constant_variable(k)
            sage: rn.rate_constants()
            (k_0_1, k_1_0)
            sage: rn.plot(edge_labels=True)
            Graphics object consisting of 8 graphics primitives
        """
        self._rate_constant_variable = variable

    def plot(
            self,
            kinetic_order: bool = True,
            layout: str = "spring",
            edge_labels: bool = False,
            vertex_colors: Union[str, List[str]] = "white",
            vertex_size: int = 6000,
            **kwargs
        ) -> None:
        r"""
        Plot the reaction network.

        This method visualizes the reaction network as a directed graph.
        The vertices represent complexes, and the edges represent reactions.
        The layout, labels, and other visual properties can be customized.

        INPUT:

        - ``kinetic_order`` (bool, default: True):
          If True, displays both stoichiometric and kinetic-order complexes
          for each vertex. If False, only stoichiometric complexes are shown.

        - ``layout`` (str, default: "spring"):
          Specifies the layout of the graph. Common options include
          "circular" and "spring".

        - ``edge_labels`` (bool, default: False):
          If True, displays the rate constants as labels on the edges.

        - ``vertex_colors`` (str or list, default: "white"):
          Specifies the color of the vertices. Can be a single color or a
          list of colors corresponding to each vertex.

        - ``vertex_size`` (int, default: 6000):
          Specifies the size of the vertices in the plot.

        - ``**kwargs``:
          Additional keyword arguments passed to the underlying graph plotting
          function.

        OUTPUT:

        - A graphical representation of the reaction network.

        EXAMPLES::

            sage: from sign_vector_conditions import *
            sage: species("A, B, C")
            (A, B, C)
            sage: rn = ReactionNetwork()
            sage: rn.add_complex(0, A + B)
            sage: rn.add_complex(1, C)
            sage: rn.add_reactions([(0, 1), (1, 0)])
            sage: rn.plot()
            Graphics object consisting of 6 graphics primitives
        """
        if edge_labels:
            self._update_edge_labels()
        return self.graph.plot(
            vertex_labels={i: self._vertex_label(i, kinetic_order=kinetic_order) for i in self.graph.vertices()},
            layout=layout,
            edge_labels=edge_labels,
            # edge_labels_background="transparent",
            vertex_colors=vertex_colors,
            vertex_size=vertex_size,
            **kwargs
        )

    def _update_edge_labels(self) -> None:
        for edge in self.reactions:
            # plot does not use latex representation for edge labels
            # f-string would mess up braces
            self.graph.set_edge_label(*edge, "$" + latex(self._rate_constant(*edge)) + "$")

    def _vertex_label(self, i: int, kinetic_order: bool = False) -> str:
        if not kinetic_order or self.complexes_stoichiometric[i] == self.complexes_kinetic_order[i]:
            return f"${latex(self.complexes_stoichiometric[i])}$"
        return f"${latex(self.complexes_stoichiometric[i])}$\n${latex(self.complexes_kinetic_order[i])}$"

    def _update(self) -> None:
        if not self._update_needed:
            return
        self._update_species()
        self._update_matrices()
        self._compute_deficiencies()
        self._update_needed = False

    def _update_species(self) -> None:
        self._species = sorted(
            set(
                species
                for complex in self.complexes_stoichiometric.values()
                for species in complex.involved_species()
            ).union(
                species
                for complex in self.complexes_kinetic_order.values()
                for species in complex.involved_species()
            )
        )

    def _update_matrices(self) -> None:
        self._matrix_of_complexes_stoichiometric = self._matrix_from_complexes(self.complexes_stoichiometric)
        self._matrix_of_complexes_kinetic_order = self._matrix_from_complexes(self.complexes_kinetic_order)
        self._matrix_stoichiometric = self.incidence_matrix().T * self._matrix_of_complexes_stoichiometric
        self._matrix_kinetic_order = self.incidence_matrix().T * self._matrix_of_complexes_kinetic_order
        self._matrix_stoichiometric_reduced = self._matrix_stoichiometric.matrix_from_rows(self._matrix_stoichiometric.pivot_rows())
        self._matrix_kinetic_order_reduced = self._matrix_kinetic_order.matrix_from_rows(self._matrix_kinetic_order.pivot_rows())

    def _compute_deficiencies(self) -> None:
        self._deficiency_stoichiometric = len(self.complexes_stoichiometric) - self.graph.connected_components_number() - self._matrix_stoichiometric_reduced.nrows()
        self._deficiency_kinetic_order = len(self.complexes_stoichiometric) - self.graph.connected_components_number() - self._matrix_kinetic_order_reduced.nrows()

    def _get(self, element: str):
        self._update()
        return getattr(self, element)

    def _matrix_from_complexes(self, complexes: Dict[int, Complex]) -> matrix:
        return matrix([complex[s] for s in self._species] for _, complex in sorted(complexes.items()))

    def are_both_deficiencies_zero(self) -> bool:
        r"""Return whether both deficiencies are zero."""
        self._update()
        return self._deficiency_stoichiometric == self._deficiency_kinetic_order == 0

    def is_weakly_reversible(self) -> bool:
        r"""Return whether each component of the system is strongly connected."""
        return all(g.is_strongly_connected() for g in self.graph.connected_components_subgraphs())

    def _check_network_conditions(self) -> None:
        self._update()
        if self._deficiency_stoichiometric != 0:
            raise ValueError(
                f"Stoichiometric deficiency should be zero, but got {self._deficiency_stoichiometric}. "
                "Ensure the network satisfies the deficiency-zero condition."
            )
        if self._deficiency_kinetic_order != 0:
            raise ValueError(
                f"Kinetic-order deficiency should be zero, but got {self._deficiency_kinetic_order}. "
                "Ensure the network satisfies the deficiency-zero condition."
            )
        if not self.is_weakly_reversible():
            raise ValueError("The network is not weakly reversible. Ensure all components are strongly connected.")

    def has_robust_cbe(self) -> bool:
        r"""
        Check whether there is a unique positive CBE in every stoichiometric class,
        for all rate constants and for all small perturbations of the kinetic orders.
        """
        self._check_network_conditions()
        return condition_closure_minors(self._matrix_stoichiometric_reduced, self._matrix_kinetic_order_reduced)

    def has_at_most_one_cbe(self) -> bool:
        r"""
        Check whether there is at most one positive CBE in every stoichiometric class,
        for all rate constants.
        """
        return condition_uniqueness_minors(self._matrix_stoichiometric_reduced, self._matrix_kinetic_order_reduced)

    def _condition_faces(self) -> bool:
        r"""Check whether the system satisfies the face condition for existence of a unique positive CBE."""
        self._check_network_conditions()
        return condition_faces(self._matrix_stoichiometric_reduced, self._matrix_kinetic_order_reduced)

    def _are_subspaces_nondegenerate(self) -> bool:
        r"""Check whether the system satisfies the nondegenerate condition for existence of a unique positive CBE."""
        self._check_network_conditions()
        return condition_nondegenerate(self._matrix_stoichiometric_reduced, self._matrix_kinetic_order_reduced)

    def has_exactly_one_cbe(self) -> bool:
        r"""
        Check whether there is a unique positive CBE in every stoichiometric class
        and for all rate constants.

        .. NOTE::

            This method does not support symbolic expressions (variables) in the complexes.
        """
        self._check_network_conditions()
        at_most_one = self.has_at_most_one_cbe()
        if at_most_one not in [True, False]:
            raise ValueError("Method does not support variables!")
        return at_most_one and self._condition_faces() and self._are_subspaces_nondegenerate()
