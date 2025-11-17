r"""
Define species, complexes, and reaction networks.

Reaction networks
=================

This module provides tools for defining species, complexes, and reaction networks
with (generalized) mass-action kinetics.
It includes functionality for analyzing reaction networks, such as checking weak reversibility,
computing deficiencies, and verifying conditions for positive complex-balanced equilibria (CBE)
based on [MHR19]_ and [AMR24]_.

EXAMPLES:

First, we define species for a reaction network::

    sage: from sign_crn import *
    sage: A, B, C = species("A, B, C")

Now, we create a reaction network and add complexes and reactions::

    sage: rn = ReactionNetwork()
    sage: rn.add_complex(0, A + B)
    sage: rn.add_complex(1, C)
    sage: rn.add_reactions([(0, 1), (1, 0)])
    sage: rn.plot()
    Graphics object consisting of 6 graphics primitives

We analyze the reaction network::

    sage: rn.is_weakly_reversible()
    True
    sage: rn.deficiency_stoichiometric
    0
    sage: rn.has_exactly_one_cbe()
    True
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
from sage.matrix.constructor import Matrix
from sage.matrix.special import diagonal_matrix
from sage.modules.free_module_element import vector
from sage.misc.latex import latex
from sage.misc.misc_c import prod

from elementary_vectors import circuit_kernel_matrix
from sign_vectors import sign_vector

from .uniqueness import condition_uniqueness_minors
from .unique_existence import condition_faces, condition_nondegenerate
from .robustness import condition_closure_minors
from .utility import non_negative_covectors_from_matrix, non_negative_vectors_from_matrix


def species(names: str) -> Union[Complex, Tuple[Complex, ...]]:
    r"""
    Define species from a string of names.

    The string can contain species names separated by commas or spaces.
    The function defines each species globally similar as :func:`var`.

    See :class:`Complex` for operations and more details.

    EXAMPLES::

        sage: from sign_crn import *
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
        sage: A, B = species("H_2, O_2")
        sage: A
        H_2
        sage: B
        O_2
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
        complex = Complex({_Species(name): 1})
        caller_globals[name] = complex
        return complex

    if len(names_list) == 1:
        return define_species_globally(names_list[0])
    return tuple(define_species_globally(name) for name in names_list)


class _Species(NamedTuple):
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

    This class represents a linear combination of species with coefficients,
    supporting various operations such as addition, subtraction, and scalar multiplication.
    It also supports symbolic expressions and substitution of values for variables.

    EXAMPLES:

    First, we define some species::

        sage: from sign_crn import *
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

    Similar to polynomials, we can substitute values for variables in a complex::

        sage: complex = a * A - B
        sage: complex
        a*A - B
        sage: complex(a=0)
        -B
        sage: complex(a=1)
        A - B
        sage: (a * A)(a=0)
        0
        sage: A(a=1)
        A

    TESTS:

    Operations with invalid types raise appropriate errors::

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
        sage: (2 * A + 3 * B).get_coefficient(A)
        2

    LaTeX representation is supported for better visualization::

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
    def __init__(self, species_dict: Dict[_Species, Union[int, float, var]]) -> None:
        self.species_dict = {}
        for key, value in species_dict.items():
            if not isinstance(key, _Species):
                raise TypeError(f"Key {key} is not a _Species")
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

    def get_coefficient(self, species: Union[_Species, Complex]) -> Union[int, float]:
        r"""
        Return the coefficient of the species in the complex.

        If a species is not present, return 0.
        """
        if isinstance(species, _Species):
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

    def __eq__(self, other) -> bool:
        r"""
        Check equality of two complexes.

        TESTS::

            sage: from sign_crn import *
            sage: species("A")
            A
            sage: var("a")
            a
            sage: complex1 = A
            sage: complex2 = a * A
            sage: complex1 == complex2(a=1)
            True
        """
        if isinstance(other, Complex):
            return self.species_dict == other.species_dict
        return False

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

    def _repr_coefficient(self, key: _Species) -> str:
        return self._format_coefficient(key, str)

    def _latex_coefficient(self, key: _Species) -> str:
        return self._format_coefficient(key, latex)

    def _format_coefficient(self, key: _Species, formatter) -> str:
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

    def _to_species(self) -> _Species:
        if len(self.species_dict) != 1:
            raise ValueError("Complex must contain exactly one species.")
        return next(iter(self.species_dict.keys()))

    def involved_species(self) -> set[_Species]:
        r"""Return the species involved in the complex."""
        return set(self.species_dict.keys())

    @staticmethod
    def from_species(species: _Species) -> Complex:
        r"""Return a complex from a species."""
        return Complex({species: 1})


class ReactionNetwork(SageObject):
    r"""
    A reaction network with (generalized) mass-action kinetics.

    This class represents a reaction network, where complexes are connected
    by directed reactions. It supports generalized mass-action kinetics, allowing
    for symbolic rate constants and kinetic orders.

    The ``ReactionNetwork`` class provides tools for:

    - Adding and removing complexes and reactions.
    - Computing stoichiometric and kinetic-order matrices.
    - Analyzing network properties, such as weak reversibility and deficiencies.
    - Checking conditions for unique positive complex-balanced equilibria (CBE).
    - Visualizing the reaction network as a directed graph.

    EXAMPLES:

    We define a reaction network with two complexes involving variables in the kinetic orders::

        sage: from sign_crn import *
        sage: var("a, b")
        (a, b)
        sage: species("A, B, C")
        (A, B, C)
        sage: rn = ReactionNetwork()
        sage: rn.add_complex(0, A + B, a * A + b * B)
        sage: rn.add_complex(1, C)
        sage: rn.add_reactions([(0, 1), (1, 0)])
        sage: rn
        Reaction network with 2 complexes, 2 reactions and 3 species.
        sage: rn.complexes_stoichiometric
        {0: A + B, 1: C}
        sage: rn.complexes_kinetic_order
        {0: a*A + b*B, 1: C}
        sage: rn.reactions
        [(0, 1), (1, 0)]
        sage: rn.species
        (A, B, C)
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

        sage: rn.stoichiometric_matrix
        [-1  1]
        [-1  1]
        [ 1 -1]
        sage: rn.kinetic_order_matrix
        [-a  a]
        [-b  b]
        [ 1 -1]
        sage: rn.stoichiometric_matrix_as_kernel
        [1 0 1]
        [0 1 1]
        sage: rn.kinetic_order_matrix_as_kernel
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
        Reaction network with 5 complexes, 5 reactions and 5 species.
        sage: rn.plot()
        Graphics object consisting of 15 graphics primitives

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
        sage: rn.ode_rhs()
        (-k_0_1*x_A^a*x_B^b + k_2_0*x_A^c*x_D - k_3_4*x_A + k_1_0*x_C + k_4_3*x_E,
         -k_0_1*x_A^a*x_B^b + k_2_0*x_A^c*x_D + k_1_0*x_C,
         k_0_1*x_A^a*x_B^b - (k_1_0 + k_1_2)*x_C,
         -k_2_0*x_A^c*x_D + k_1_2*x_C,
         k_3_4*x_A - k_4_3*x_E)

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
        sage: rn.stoichiometric_matrix
        [-1  1  0  1 -1  1]
        [-1  1  0  1  0  0]
        [ 1 -1 -1  0  0  0]
        [ 0  0  1 -1  0  0]
        [ 0  0  0  0  1 -1]
        sage: rn.kinetic_order_matrix
        [   -a     a     c a - c    -1     1]
        [   -b     b     0     b     0     0]
        [    1    -1    -1     0     0     0]
        [    0     0     1    -1     0     0]
        [    0     0     0     0     1    -1]
        sage: rn.stoichiometric_matrix_as_kernel
        [1 0 1 1 1]
        [0 1 1 1 0]
        sage: rn.kinetic_order_matrix_as_kernel
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
        sage: rn(a=2, b=1, c=1).has_at_most_one_cbe()
        True
        sage: rn.has_at_most_one_cbe() # random order
        [{a >= 0, a - c >= 0, b >= 0}]
        sage: rn.has_exactly_one_cbe()
        Traceback (most recent call last):
        ...
        ValueError: Method does not support variables.
        sage: rn(a=2, b=1, c=1).has_exactly_one_cbe()
        True

    We remove one component and a reaction of our system::

        sage: rn.remove_complex(3)
        sage: rn.remove_complex(4)
        sage: rn.remove_reaction(1, 0)
        sage: rn
        Reaction network with 3 complexes, 3 reactions and 4 species.

    Here is an example involving molecules::

        sage: A, B, C = species("H_2, O_2, H_2O")
        sage: var("a")
        a
        sage: rn = ReactionNetwork()
        sage: rn.add_complex(0, 2 * A + B, 2 * a * A + a * B)
        sage: rn.add_complex(1, 2 * C)
        sage: rn.species
        (H_2, H_2O, O_2)
        sage: rn.add_reactions([(0, 1), (1, 0)])
        sage: rn.plot()
        Graphics object consisting of 6 graphics primitives

    We can also define an ecosystem involving animals as species::

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

        self._species: List[_Species] = []

        self._matrix_of_complexes_stoichiometric = None
        self._matrix_of_complexes_kinetic_order = None
        self._stoichiometric_matrix = None
        self._kinetic_order_matrix = None
        self._stoichiometric_matrix_reduced = None
        self._kinetic_order_matrix_reduced = None

        self._deficiency_stoichiometric: int = 0
        self._deficiency_kinetic_order: int = 0

    def _repr_(self) -> str:
        return f"Reaction network with {self.graph.num_verts()} complexes, {self.graph.num_edges()} reactions and {len(self.species)} species."

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
    def matrix_of_complexes_stoichiometric(self) -> Matrix:
        r"""
        Return the matrix that decodes the stoichiometric complexes of the reaction network.

        Each column stands for a complex, and each row stands for a species.
        """
        return self._get("_matrix_of_complexes_stoichiometric").T

    @property
    def matrix_of_complexes_kinetic_order(self) -> Matrix:
        r"""
        Return the matrix that decodes the kinetic-order complexes of the reaction network.

        Each column stands for a complex, and each row stands for a species.
        """
        return self._get("_matrix_of_complexes_kinetic_order").T

    @property
    def stoichiometric_matrix(self) -> Matrix:
        r"""
        Return the stoichiometric matrix where the columns correspond to the reactions.

        Each columns stands for a reaction, and each row stands for a species.
        """
        return self._get("_stoichiometric_matrix").T

    @property
    def kinetic_order_matrix(self) -> Matrix:
        r"""
        Return the kinetic-order matrix where the columns correspond to the reactions.

        Each columns stands for a reaction, and each row stands for a species.
        """
        return self._get("_kinetic_order_matrix").T

    @property
    def stoichiometric_matrix_as_kernel(self) -> Matrix:
        r"""Return the kernel matrix of the stoichiometric matrix."""
        self._update()
        return circuit_kernel_matrix(self._stoichiometric_matrix_reduced)

    @property
    def kinetic_order_matrix_as_kernel(self) -> Matrix:
        r"""Return the kernel matrix of the kinetic-order matrix."""
        self._update()
        return circuit_kernel_matrix(self._kinetic_order_matrix_reduced)

    @property
    def deficiency_stoichiometric(self) -> int:
        r"""Return the stoichiometric deficiency."""
        return self._get("_deficiency_stoichiometric")

    @property
    def deficiency_kinetic_order(self) -> int:
        r"""Return the kinetic-order deficiency."""
        return self._get("_deficiency_kinetic_order")

    def incidence_matrix(self, **kwargs) -> Matrix:
        r"""Return the incidence matrix of the graph."""
        return self.graph.incidence_matrix(**kwargs)

    def source_matrix(self, **kwargs) -> Matrix:
        r"""Return the source matrix of the graph."""
        return Matrix((1 if value == -1 else 0 for value in row) for row in self.incidence_matrix(**kwargs))

    def laplacian_matrix(self) -> Matrix:
        r"""Return the Laplacian matrix of the graph."""
        return self.incidence_matrix() * diagonal_matrix(self.rate_constants()) * self.source_matrix().T

    def ode_rhs(self) -> vector:
        r"""Return the right hand side of the ordinary differential equation of this system."""
        self._update()
        x = vector(var(f"x_{s}", latex_name=f"x_{{{s}}}") for s in self.species)
        return (
            self._matrix_of_complexes_stoichiometric.T * self.laplacian_matrix() * vector(
                prod(xi ** yi for xi, yi in zip(x, y))
                for y in self._matrix_of_complexes_kinetic_order.rows()
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

            sage: from sign_crn import *
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

        self._stoichiometric_matrix = self.incidence_matrix().T * self._matrix_of_complexes_stoichiometric
        self._kinetic_order_matrix = self.incidence_matrix().T * self._matrix_of_complexes_kinetic_order

        self._stoichiometric_matrix_reduced = self._stoichiometric_matrix.matrix_from_rows(self._stoichiometric_matrix.pivot_rows())
        self._kinetic_order_matrix_reduced = self._kinetic_order_matrix.matrix_from_rows(self._kinetic_order_matrix.pivot_rows())

    def _matrix_from_complexes(self, complexes: Dict[int, Complex]) -> Matrix:
        return Matrix([complexes[v].get_coefficient(s) for s in self._species] for v in self.graph.vertices())

    def _compute_deficiencies(self) -> None:
        connected_components_number = self.graph.connected_components_number()
        self._deficiency_stoichiometric = len(self.complexes_stoichiometric) - connected_components_number - self._stoichiometric_matrix_reduced.nrows()
        self._deficiency_kinetic_order = len(self.complexes_stoichiometric) - connected_components_number - self._kinetic_order_matrix_reduced.nrows()

    def _get(self, element: str):
        self._update()
        return getattr(self, element)

    def plot(
            self,
            kinetic_orders: bool = True,
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

            sage: from sign_crn import *
            sage: species("A, B, C")
            (A, B, C)
            sage: rn = ReactionNetwork()
            sage: rn.add_complex(0, A + B)
            sage: rn.add_complex(1, C)
            sage: rn.add_reactions([(0, 1), (1, 0)])
            sage: rn.plot()
            Graphics object consisting of 6 graphics primitives

        We can customize plotting::

            sage: rn.plot(edge_labels=True)
            Graphics object consisting of 8 graphics primitives
            sage: rn.plot(kinetic_orders=False)
            Graphics object consisting of 6 graphics primitives
            sage: rn.plot(vertex_size=3000)
            Graphics object consisting of 6 graphics primitives
            sage: rn.plot(layout="circular")
            Graphics object consisting of 6 graphics primitives
        """
        if edge_labels:
            self._update_edge_labels()
        return self.graph.plot(
            vertex_labels={i: self._vertex_label(i, show_kinetic_order=kinetic_orders) for i in self.graph.vertices()},
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

    def _vertex_label(self, i: int, show_kinetic_order: bool = False) -> str:
        if not show_kinetic_order or self.complexes_stoichiometric[i] == self.complexes_kinetic_order[i]:
            return f"${latex(self.complexes_stoichiometric[i])}$"
        return f"${latex(self.complexes_stoichiometric[i])}$\n$({latex(self.complexes_kinetic_order[i])})$"

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
        return condition_closure_minors(self._stoichiometric_matrix_reduced, self._kinetic_order_matrix_reduced)

    def has_at_most_one_cbe(self) -> bool:
        r"""
        Check whether there is at most one positive CBE in every stoichiometric class
        and for all rate constants.
        """
        self._update()
        return condition_uniqueness_minors(self._stoichiometric_matrix_reduced, self._kinetic_order_matrix_reduced)

    def _condition_faces(self) -> bool:
        r"""Check whether the system satisfies the face condition for existence of a unique positive CBE."""
        self._check_network_conditions()
        return condition_faces(self._stoichiometric_matrix_reduced, self._kinetic_order_matrix_reduced)

    def _are_subspaces_nondegenerate(self) -> bool:
        r"""Check whether the system satisfies the nondegenerate condition for existence of a unique positive CBE."""
        self._check_network_conditions()
        return condition_nondegenerate(self._stoichiometric_matrix_reduced, self._kinetic_order_matrix_reduced)

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
            raise ValueError("Method does not support variables.")
        return at_most_one and self._condition_faces() and self._are_subspaces_nondegenerate()

    def has_exactly_one_equilibrium(self) -> bool:
        r"""
        Check whether there is a unique positive complex-balanced equilibrium.

        .. NOTE::

            This method does not support symbolic expressions (variables) in the complexes.

        EXAMPLES::

            sage: from sign_crn import *
            sage: species("X, Y")
            (X, Y)
            sage: rn = ReactionNetwork()
            sage: rn.add_complexes([(1, 2 * X + Y), (2, 2 * Y), (3, 3 * X + Y), (4, 4 * X), (5, 3 * X + 2 * Y), (6, 3 * X + 3 * Y), (7, 2 * X + 4 * Y)])
            sage: rn.add_reactions([(1, 2), (3, 4), (5, 6), (6, 7)])
            sage: rn.has_exactly_one_equilibrium()
            True
            sage: rn.remove_reaction(6, 7)
            sage: rn.has_exactly_one_equilibrium()
            True
            sage: rn.remove_reaction(1, 2)
            sage: rn.has_exactly_one_equilibrium()
            False

        ::

            sage: rn = ReactionNetwork()
            sage: rn.add_complexes([(0, 5 * X + Y), (1, 6 * X), (2, 4 * X + 3 * Y), (3, 5 * X + 6 * Y), (4, 2 * X + 3 * Y), (5, 4 * Y), (6, 4 * X + 2 * Y), (8, 3 * X + 3 * Y), (9, 2 * X + 5 * Y)])
            sage: rn.add_reactions([(0, 1), (2, 3), (4, 5), (6, 2), (8, 9)])
            sage: rn.has_exactly_one_equilibrium()
            True
            sage: rn.remove_reaction(0, 1)
            sage: rn.has_exactly_one_equilibrium()
            False
        """
        self._update()

        def matrix_with_one_row(A):
            return A.insert_row(0, vector([1] * A.ncols()))

        A = self._stoichiometric_matrix.T
        B = self._matrix_of_complexes_kinetic_order.T * self.source_matrix()
        A_bar = matrix_with_one_row(A)
        B_bar = matrix_with_one_row(B)

        first_condition = condition_uniqueness_minors(A_bar, B_bar)
        if first_condition not in [True, False]:
            raise ValueError("Method does not support variables.")
        if not first_condition:
            return False

        covectors = non_negative_covectors_from_matrix(B_bar)
        all_positive_found = False
        for face in non_negative_vectors_from_matrix(A):
            if face == 0:
                continue
            mirrored_face = sign_vector(1 if fe == 0 else 0 for fe in face)
            if mirrored_face == 0:
                all_positive_found = True
                continue
            if mirrored_face in covectors:
                return False
        return all_positive_found
