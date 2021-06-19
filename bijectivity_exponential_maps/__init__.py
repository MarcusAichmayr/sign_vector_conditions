#############################################################################
#  Copyright (C) 2021                                                       #
#                Marcus Aichmayr (aichmayr.marcus@gmail.com)                #
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

