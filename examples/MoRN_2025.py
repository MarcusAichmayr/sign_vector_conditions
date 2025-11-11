r"""
MoRN seminar software showcase 2025.

============================================
Sign Vector Conditions for Reaction Networks
============================================

A simple example
================

We define a reaction network::

    sage: from sign_crn import *
    sage: species("H_2, O_2, H_2O")
    (H_2, O_2, H_2O)
    sage: rn = ReactionNetwork()
    sage: rn.add_complex(0, 2 * H_2 + O_2)
    sage: rn.add_complex(1, 2 * H_2O)
    sage: rn.add_reactions([(0, 1), (1, 0)])
    sage: rn
    Reaction network with 2 complexes, 2 reactions and 3 species.
    sage: rn.plot()
    Graphics object consisting of 6 graphics primitives

Mass-Action Kinetics (MAK)
==========================

Kinetic orders (exponents) are stoichiometric coefficients.

:math:`1 A + 1 B \to C` with reaction rate

.. MATH::

    k (x_A)^1 (x_B)^1

We define a reaction network with MAK::

    sage: species("A, B, C, D, E")
    (A, B, C, D, E)
    sage: rn = ReactionNetwork()
    sage: rn.add_complexes([(0, A + B), (1, C), (2, D), (3, A), (4, E)])
    sage: rn.add_reactions([(0, 1), (1, 0), (1, 2), (2, 0), (3, 4), (4, 3)])
    sage: rn
    Reaction network with 5 complexes, 6 reactions and 5 species.
    sage: rn.plot()
    Graphics object consisting of 16 graphics primitives

The exponents are the stoichiometric coefficients of the reaction network::

    sage: rn.ode_rhs()
    (-k_0_1*x_A*x_B - k_3_4*x_A + k_1_0*x_C + k_2_0*x_D + k_4_3*x_E,
     -k_0_1*x_A*x_B + k_1_0*x_C + k_2_0*x_D,
     k_0_1*x_A*x_B - (k_1_0 + k_1_2)*x_C,
     k_1_2*x_C - k_2_0*x_D,
     k_3_4*x_A - k_4_3*x_E)

Generalized Mass-Action Kinetics (GMAK)
=======================================

Kinetic orders are arbitrary real numbers.
(different from stoichiometric coefficients)

:math:`1 A + 1 B \to C` with reaction rate

.. MATH::

    k (x_A)^a (x_B)^b

denoted as :math:`(a A + b B)`

We define a reaction network with GMAK.
First, we define symbolic exponents::

    sage: var("a, b, c")
    (a, b, c)

We update our reaction network with kinetic-order complexes::

    sage: rn.add_complexes([(0, A + B, a * A + b * B), (2, D, c * A + D)])
    sage: rn.plot()
    Graphics object consisting of 16 graphics primitives
    sage: rn.ode_rhs()
    (-k_0_1*x_A^a*x_B^b + k_2_0*x_A^c*x_D - k_3_4*x_A + k_1_0*x_C + k_4_3*x_E,
     -k_0_1*x_A^a*x_B^b + k_2_0*x_A^c*x_D + k_1_0*x_C,
     k_0_1*x_A^a*x_B^b - (k_1_0 + k_1_2)*x_C,
     -k_2_0*x_A^c*x_D + k_1_2*x_C,
     k_3_4*x_A - k_4_3*x_E)

Matrices
========

The network is defined by the following matrices::

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
    sage: rn.incidence_matrix()
    [-1  1  0  1  0  0]
    [ 1 -1 -1  0  0  0]
    [ 0  0  1 -1  0  0]
    [ 0  0  0  0 -1  1]
    [ 0  0  0  0  1 -1]
    sage: rn.laplacian_matrix()
    [        -k_0_1          k_1_0          k_2_0              0              0]
    [         k_0_1 -k_1_0 - k_1_2              0              0              0]
    [             0          k_1_2         -k_2_0              0              0]
    [             0              0              0         -k_3_4          k_4_3]
    [             0              0              0          k_3_4         -k_4_3]

Network properties
==================

We show some properties of the network::

    sage: rn.is_weakly_reversible()
    True
    sage: rn.deficiency_stoichiometric
    0
    sage: rn.deficiency_kinetic_order
    0

Sign vector conditions
======================

Uniqueness and existence of complex-balanced equilibria (CBE)
-------------------------------------------------------------

See [MHR19]_.

Existence
^^^^^^^^^

First, we instantiate the symbolic exponents::

    sage: rn_instantiated = rn(a=2, b=1, c=1)
    sage: rn_instantiated.plot()
    Graphics object consisting of 16 graphics primitives

There is at most one CBE::

    sage: rn_instantiated.has_at_most_one_cbe()
    True

We can also apply this method for reaction networks involving symbolic exponents.
In that case, we obtain conditions on these::

    sage: rn.has_at_most_one_cbe() # random order
    [{a >= 0, a - c >= 0, b >= 0}]

Robustness
^^^^^^^^^^

There is at most one CBE for all small perturbations of kinetic orders::

    sage: rn.has_robust_cbe() # random order
    [{a > 0, a - c > 0, b > 0}]

Exactly one CBE (Deficiency Zero Theorem for GMAK)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To check this condition, we need sign vector conditions.
There is also an involved nondegeneracy condition
which we check using a recursive algorithm
that certifies (un)solvability of linear inequality systems with elementary vectors.
See [AMR24]_ for details.

There is exactly one CBE::

    sage: rn_instantiated.has_exactly_one_cbe()
    True

Elementary and sign vectors
---------------------------

We define a matrix::

    sage: from elementary_vectors import *
    sage: from sign_vectors import *
    sage: from sign_vectors.oriented_matroids import *
    sage: P = matrix([[1, 2, 0, 0], [0, 1, 2, 3]])
    sage: P
    [1 2 0 0]
    [0 1 2 3]

We compute the elements with minimal-support in the kernel of `P`::

    sage: circuits(P)
    [(4, -2, 1, 0), (6, -3, 0, 1), (0, 0, -3, 2)]

Next, we compute the sign vectors with minimal support::

    sage: om = OrientedMatroid(P)
    sage: om.circuits()
    {(00-+), (+-0+), (-+0-), (00+-), (-+-0), (+-+0)}

The sign vectors of the corresponding oriented matroid are::

    sage: om.vectors()
    {(0000),
     (00-+),
     (-+-+),
     (+-0+),
     (+--+),
     (-+0-),
     (-+--),
     (00+-),
     (+-+0),
     (+-+-),
     (+-++),
     (-+-0),
     (-++-)}
    sage: plot_sign_vectors(om.vectors()) # random
"""
