Sign Vector Conditions for Reaction Networks
============================================

To install this package, visit the `repository on GitHub <https://github.com/MarcusAichmayr/sign_vector_conditions>`_.

.. rubric:: Reaction networks
.. autosummary::
    :toctree: generated

    sign_vector_conditions.reaction_networks

.. rubric:: Sign vector conditions

The conditions for existence and uniqueness of equilibria for (chemical) reaction networks
are based on [MHR19]_.
See [AMR24]_ for a more details.

.. autosummary::
    :toctree: generated

    sign_vector_conditions.uniqueness
    sign_vector_conditions.unique_existence
    sign_vector_conditions.robustness

.. rubric:: Conferences and presentations

.. autosummary::
    :toctree: generated

    examples.ICMS_2024
    examples.MoRN_2025

.. rubric:: References

.. [AMR24] Marcus S. Aichmayr, Stefan Müller, and Georg Regensburger.
    "A SageMath package for elementary and sign vectors with applications to chemical reaction networks".
    In: Mathematical software—ICMS 2024 (2024),
    pp. 155--164. doi: 10.1007/978-3-031-64529-7_17.

.. [MHR19] Müller, S.; Hofbauer, J., and Regensburger, G.:
    "On the bijectivity of families of exponential/generalized polynomial maps".
    In: SIAM Journal on Applied Algebra and Geometry 3.3 (2019),
    pp. 412--438. doi: 10.1137/18M1178153.
