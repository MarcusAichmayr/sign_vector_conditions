r"""
In this package, we check several conditons from [MHR19]_ for exponential maps.

.. [MHR19] Müller, S.; Hofbauer, J., and Regensburger, G.:
   „On the bijectivity of families of exponential/generalized polynomial maps“.
   In: SIAM Journal on Applied Algebra and Geometry 3.3 (2019),
   pp. 412–438. doi: 10.1137/18M1178153.

.. autosummary::
    :toctree: generated

    bijectivity_exponential_maps.functions
    bijectivity_exponential_maps.conditions_injectivity
    bijectivity_exponential_maps.conditions_bijectivity
    bijectivity_exponential_maps.conditions_bijectivity_robust
    bijectivity_exponential_maps.utility
"""

#############################################################################
#  Copyright (C) 2023                                                       #
#                Marcus Aichmayr (aichmayr@mathematik.uni-kassel.de)        #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 3, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from __future__ import absolute_import

from .conditions_injectivity import cond_inj_intersection, cond_inj_minors
from .conditions_bijectivity import cond_faces, cond_nondegenerate, nondegenerate
from .conditions_bijectivity_robust import cond_closure_sign_vectors, cond_closure_minors

from .functions import f_pol, f_exp