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


class GMAKSystem(SageObject):
    r"""
    Class for chemical reaction networks with generalized mass-action kinetics.

    A generalized mass-action system is represented by a (weighted) directed graph
    and stoichiometric and kinetic-order labels of the vertices.

    EXAMPLES::

        sage: from sign_vector_conditions.chemical_reaction_networks import *
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
        sage: crn.stoichiometric_matrix()
        [-1 -1  1  0  0]
        [ 1  1 -1  0  0]
        [ 0  0 -1  1  0]
        [ 1  1  0 -1  0]
        [-1  0  0  0  1]
        [ 1  0  0  0 -1]
        sage: crn.kinetic_order_matrix()
        [   -a    -b     1     0     0]
        [    a     b    -1     0     0]
        [    c     0    -1     1     0]
        [a - c     b     0    -1     0]
        [   -1     0     0     0     1]
        [    1     0     0     0    -1]
    """

    def __init__(self, graph, stoichiometric_labels, kinetic_order_labels):
        self.graph = graph
        self.stoichiometric_labels = stoichiometric_labels
        self.kinetic_order_labels = kinetic_order_labels

    # def _repr_(self):
    #     return graph

    def incidence_matrix(self):
        return self.graph.incidence_matrix()

    def source_matrix(self):
        return matrix(
            (1 if value == -1 else 0 for value in row)
            for row in self.graph.incidence_matrix()
        )

    def stoichiometric_matrix(self):
        return self.incidence_matrix().T * self.stoichiometric_labels

    def kinetic_order_matrix(self):
        return self.incidence_matrix().T * self.kinetic_order_labels
