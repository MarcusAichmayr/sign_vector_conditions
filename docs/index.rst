Sign vector conditions for (chemical) reaction networks
=======================================================

This package provides functionality for analyzing reaction networks,
such as checking weak reversibility, computing deficiencies, and verifying conditions
for unique positive complex-balanced equilibria (CBE) based on [MHR19]_.
See [AMR24]_ for more details.

To install this package, visit the `repository on GitHub <https://github.com/MarcusAichmayr/sign_crn>`_.

.. rubric:: Reaction networks

.. autosummary::
    :toctree: generated

    sign_crn.reaction_networks

.. rubric:: Conferences and presentations

.. autosummary::
    :toctree: generated

    examples.MoRN_2025
    examples.ICMS_2024

.. rubric:: Sign vector conditions

.. autosummary::
    :toctree: generated

    sign_crn.uniqueness
    sign_crn.unique_existence
    sign_crn.robustness
    sign_crn.utility

.. rubric:: References

.. [AMR24] Marcus S. Aichmayr, Stefan Müller, and Georg Regensburger.
    "A SageMath package for elementary and sign vectors with applications to chemical reaction networks".
    In: Mathematical software—ICMS 2024 (2024),
    pp. 155--164. doi: `10.1007/978-3-031-64529-7_17 <https://doi.org/10.1007/978-3-031-64529-7_17>`_.

.. [MHR19] Müller, S.; Hofbauer, J., and Regensburger, G.:
    "On the bijectivity of families of exponential/generalized polynomial maps".
    In: SIAM Journal on Applied Algebra and Geometry 3.3 (2019),
    pp. 412--438. doi: `10.1137/18M1178153 <https://doi.org/10.1137/18M1178153>`_.
