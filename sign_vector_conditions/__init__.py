r"""
Conditions from [MHR19]_ for chemical reaction networks.

.. [MHR19] Müller, S.; Hofbauer, J., and Regensburger, G.:
   "On the bijectivity of families of exponential/generalized polynomial maps".
   In: SIAM Journal on Applied Algebra and Geometry 3.3 (2019),
   pp. 412--438. doi: 10.1137/18M1178153.

.. autosummary::
    :toctree: generated

    reaction_networks
    uniqueness
    unique_existence
    robustness
    utility
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

from __future__ import absolute_import

from .reaction_networks import ReactionNetwork, species
from .uniqueness import condition_uniqueness_sign_vectors, condition_uniqueness_minors
from .unique_existence import condition_faces, condition_nondegenerate, condition_degenerate
from .robustness import condition_closure_sign_vectors, condition_closure_minors
