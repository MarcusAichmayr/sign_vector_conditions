r"""
In this package, we check several conditons from [MHR19]_ for exponential maps.

.. [MHR19] Müller, S.; Hofbauer, J., and Regensburger, G.:
   „On the bijectivity of families of exponential/generalized polynomial maps“.
   In: SIAM Journal on Applied Algebra and Geometry 3.3 (2019),
   pp. 412–438. doi: 10.1137/18M1178153.

.. autosummary::
    :toctree: generated

    sign_vector_conditions.functions
    sign_vector_conditions.conditions_injectivity
    sign_vector_conditions.conditions_bijectivity
    sign_vector_conditions.conditions_bijectivity_robust
    sign_vector_conditions.utility
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

from .conditions_injectivity import condition_inj_intersection, condition_inj_minors
from .conditions_bijectivity import condition_faces, condition_nondegenerate, nondegenerate
from .conditions_bijectivity_robust import condition_closure_sign_vectors, condition_closure_minors

from .functions import f_pol, f_exp