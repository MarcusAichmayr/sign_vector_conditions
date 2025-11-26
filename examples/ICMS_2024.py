r"""
Up-to-date examples of [AMR24]_.

A SageMath Package for Elementary and Sign Vectors with Applications to Chemical Reaction Networks
--------------------------------------------------------------------------------------------------

Here are the up-to-date examples appearing in [AMR24]_ for `ICMS 2024 <https://icms-conference.org/2024/>`_.
The paper is also available at `ARXIV <https://arxiv.org/abs/2407.12660>`_.

Elementary vectors
~~~~~~~~~~~~~~~~~~

Functions dealing with elementary vectors (circuits of a subspace given by a matrix)
are implemented in the package `elementary_vectors <https://github.com/MarcusAichmayr/elementary_vectors>`_.

We compute elementary vectors (circuits), using maximal minors::

    sage: from elementary_vectors import *
    sage: M = matrix([[1, 1, 2, 0], [0, 0, 1, 2]])
    sage: M
    [1 1 2 0]
    [0 0 1 2]
    sage: M.minors(2)
    [0, 1, 2, 1, 2, 4]
    sage: circuits(M)
    [(1, -1, 0, 0), (4, 0, -2, 1), (0, 4, -2, 1)]

Solvability of linear inequality systems
****************************************

Our package `certlin <https://github.com/MarcusAichmayr/certlin>`_
provides tools for the solvability of linear inequality systems.
We state linear inequality systems as intersection of a vector space and a Cartesian product of intervals.
To represent these objects, we use a matrix and a list of intervals::

    sage: from certlin import *
    sage: M = matrix([[1, 0], [0, 1], [1, 1], [0, 1]])
    sage: M
    [1 0]
    [0 1]
    [1 1]
    [0 1]
    sage: lower = [2, 5, 0, -oo]
    sage: upper = [5, oo, 8, 5]
    sage: lower_closed = [True, True, False, False]
    sage: upper_closed = [False, False, False, True]
    sage: I = Intervals.from_bounds(lower, upper, lower_closed, upper_closed)
    sage: I
    [2, 5) x [5, +oo) x (0, 8) x (-oo, 5]
    sage: sys = LinearInequalitySystem(M, I)
    sage: sys
    [1 0]  x in [2, 5)
    [0 1]  x in [5, +oo)
    [1 1]  x in (0, 8)
    [0 1]  x in (-oo, 5]

The package can certify whether the system has a solution or not::

    sage: sys.certify()
    (True, (5/2, 5))
    sage: sys.find_solution()
    (5/2, 5)

Therefore, the system has a solution.

Sign vectors and oriented matroids
**********************************

The package `sign_vectors <https://github.com/MarcusAichmayr/sign_vectors>`_
provides functions for sign vectors and oriented matroids.
We consider an oriented matroid given by a matrix and compute the cocircuits and covectors::

    sage: from sign_vectors import *
    sage: M = matrix([[1, 3, -2, 1], [0, 4, -2, 1]])
    sage: M
    [ 1  3 -2  1]
    [ 0  4 -2  1]
    sage: om = OrientedMatroid(M)
    sage: om.cocircuits()
    {(0-+-), (+-00), (-+00), (-0+-), (0+-+), (+0-+)}
    sage: om.covectors()
    {(0000),
     (0-+-),
     (-+-+),
     (+-00),
     (-+00),
     (-0+-),
     (0+-+),
     (+--+),
     (+-+-),
     (--+-),
     (-++-),
     (+0-+),
     (++-+)}

Chemical reaction networks
~~~~~~~~~~~~~~~~~~~~~~~~~~

See :func:`sign_crn.reaction_networks.ReactionNetwork`
for a user-friendly class to define (chemical) reaction networks.

Several sign vector conditions for such networks are implemented
in the package `sign_crn <https://github.com/MarcusAichmayr/sign_crn>`_.

Robustness
**********

Given is a chemical reaction network involving five complexes.
To examine robustness of CBE, we compute the covectors corresponding to the resulting subspaces::

    sage: from sign_vectors import *
    sage: S = matrix([[-1, -1, 1, 0, 0], [0, 0, -1, 1, 0], [-1, 0, 0, 0, 1]])
    sage: S
    [-1 -1  1  0  0]
    [ 0  0 -1  1  0]
    [-1  0  0  0  1]
    sage: om = OrientedMatroid(S)
    sage: om.covectors()
    {(00000),
     (--+-+),
     (0++-+),
     (++-00),
     (--+00),
     (++-++),
     (--0+0),
     (0-+0-),
     (--++-),
     (--+++),
     (++--+),
     (0-0+-),
     (++--0),
     (++---),
     (0-+--),
     (--+--),
     (-000+),
     (-0-++),
     (+0-+-),
     (0+-0+),
     (-+-++),
     (++-+0),
     (++-+-),
     (-+--+),
     (--++0),
     (00+-0),
     (--+-0),
     (0-++-),
     (+000-),
     (+-+0-),
     (++0--),
     (++0-+),
     (+0+--),
     (00-+0),
     (---++),
     (0+-++),
     (---+-),
     (+--+-),
     (+-++-),
     (-+0-+),
     (0+--+),
     (--+0+),
     (-0+-+),
     (-++-+),
     (++0-0),
     (+++--),
     (+++-+),
     (++-0+),
     (--+0-),
     (---+0),
     (++-0-),
     (--0+-),
     (--0++),
     (+-0+-),
     (+-+--),
     (0--+-),
     (+++-0),
     (-+-0+),
     (0+0-+)}
    sage: var("a, b, c")
    (a, b, c)
    sage: St = matrix([[-a, -b, 1, 0, 0], [c, 0, -1, 1, 0], [-1, 0, 0, 0, 1]])
    sage: St
    [-a -b  1  0  0]
    [ c  0 -1  1  0]
    [-1  0  0  0  1]
    sage: OrientedMatroid(St(a=2, b=1, c=1)).covectors()
    {(00000),
     (--+-+),
     (0++-+),
     (++-00),
     (--+00),
     (-0+-0),
     (00+--),
     (0++--),
     (++-++),
     (--0+0),
     (0-+0-),
     (--++-),
     (--+++),
     (0--+0),
     (++---),
     (0-0+-),
     (--+--),
     (0+-0+),
     (++--+),
     (-000+),
     (-0-++),
     (+0-+-),
     (+0-++),
     (-+-++),
     (++-+0),
     (++-+-),
     (0++-0),
     (-+--+),
     (--++0),
     (0+0-+),
     (++--0),
     (--+-0),
     (0-++-),
     (+-+0-),
     (+0-+0),
     (0-+--),
     (++0--),
     (++0-+),
     (00-++),
     (+0+--),
     (---++),
     (---+-),
     (+--+-),
     (0+-++),
     (+-++-),
     (+--++),
     (-+0-+),
     (0+--+),
     (--+0+),
     (-0+-+),
     (-++-+),
     (-++--),
     (++0-0),
     (+++-0),
     (+++--),
     (+++-+),
     (++-0+),
     (--+0-),
     (---+0),
     (+--+0),
     (+000-),
     (++-0-),
     (-0+--),
     (0--++),
     (--0+-),
     (--0++),
     (+-0+-),
     (+-+--),
     (0--+-),
     (-+-0+),
     (-++-0)}

For :math:`a = 2`, :math:`b = 1` and :math:`c = 1`, the covectors of :math:`S` are included in the closure of the covectors of :math:`\widetilde{S}`.
To consider the general case, we compute the maximal minors of :math:`S` and :math:`\widetilde{S}`::

    sage: S  = matrix([[1, 0, 1, 1, 1], [0, 1, 1, 1, 0]])
    sage: S
    [1 0 1 1 1]
    [0 1 1 1 0]
    sage: var("a, b, c")
    (a, b, c)
    sage: St = matrix([[1, 0, a, a - c, 1], [0, 1, b, b, 0]])
    sage: St
    [    1     0     a a - c     1]
    [    0     1     b     b     0]
    sage: from sign_crn import *
    sage: condition_closure_minors(S, St) # random order
    [{a > 0, b > 0, a - c > 0}]

Hence, the network has a unique positive CBE if and only if :math:`a, b > 0` and :math:`a > c`.

Uniqueness
**********

We can also use the maximal minors to study uniqueness of CBE::

    sage: condition_uniqueness_minors(S, St) # random order
    [{a >= 0, b >= 0, a - c >= 0}]

Hence, positive CBE are unique if and only if :math:`a, b \geq 0` and :math:`a \geq c`.

Unique existence of CBE
***********************

Now, we consider Example 20 from [MHR19]_.
Here, we have a parameter :math:`a > 0`.
Depending on this parameter, the network has a unique positive CBE::

    sage: var("a")
    a
    sage: assume(a > 0)
    sage: S = matrix([[1, 0, 0, 0, 0, 1], [0, 1, 0, 0, 0, -1], [0, 0, 1, 1, 2, 0]])
    sage: S
    [ 1  0  0  0  0  1]
    [ 0  1  0  0  0 -1]
    [ 0  0  1  1  2  0]
    sage: St = matrix([[-1, -1, 0, 0, -2, 0], [0, 0, 1, 1, 0, 0], [0, 0, 0, 0, a, 1]])
    sage: St
    [-1 -1  0  0 -2  0]
    [ 0  0  1  1  0  0]
    [ 0  0  0  0  a  1]

The first two conditions depend on the sign vectors corresponding
to the rows of these matrices which are independent of the specific value for :math:`a`::

    sage: condition_uniqueness_sign_vectors(S, St)
    True

Hence, there exists at most one equilibrium.
Also the face condition is satisfied::

    sage: condition_faces(S, St)
    True

For specific values of :math:`a`, the pair of subspaces
determined by kernels of the matrices is nondegenerate.
This is exactly the case for :math:`a \in (0, 1) \cup (1, 2)`.
We demonstrate this for specific values::

    sage: condition_nondegenerate(S, St(a=1/2))
    True
    sage: condition_nondegenerate(S, St(a=3/2))
    True
    sage: condition_nondegenerate(S, St(a=1))
    False
    sage: condition_nondegenerate(S, St(a=2))
    False
    sage: condition_nondegenerate(S, St(a=3))
    False
"""
