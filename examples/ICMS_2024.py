r"""
Examples for `ICMS 2024 <https://icms-conference.org/2024/>`_.

A SageMath Package for Elementary and Sign Vectors with Applications to Chemical Reaction Networks
--------------------------------------------------------------------------------------------------

Here are the up-to-date examples appearing
in TODO (link to extended abstract).

Elementary vectors
~~~~~~~~~~~~~~~~~~
::

    sage: from elementary_vectors import *
    sage: M = matrix([[1, 1, 2, 0], [0, 0, 1, 2]])
    sage: M
    [1 1 2 0]
    [0 0 1 2]
    sage: M.minors(2)
    [0, 1, 2, 1, 2, 4]
    sage: elementary_vectors(M)
    [(1, -1, 0, 0), (4, 0, -2, 1), (0, 4, -2, 1)]

Solvability of linear inequality systems
****************************************
::

    sage: from vectors_in_intervals import *
    sage: M = matrix([[1, 0, 1, 0], [0, 1, 1, 1]])
    sage: M
    [1 0 1 0]
    [0 1 1 1]
    sage: I = intervals_from_bounds([2, 5, 0, -oo], [5, oo, 8, 5], [True, True, False, False], [False, False, False, True])
    sage: I
    [[2, 5), [5, +oo), (0, 8), (-oo, 5]]
    sage: exists_vector(M, I)
    True

Therefore, the system has a solution.

Sign vectors and oriented matroids
**********************************

We consider an oriented matroid given by a matrix and compute the cocircuits and covectors::

    sage: from sign_vectors.oriented_matroids import *
    sage: M = matrix([[1, 1, 2, 0], [0, 0, 1, 2]])
    sage: M
    [1 1 2 0]
    [0 0 1 2]
    sage: cocircuits_from_matrix(M)
    {(-+00), (-0+-), (+0-+), (+-00), (0-+-), (0+-+)}
    sage: covectors_from_matrix(M)
    {(0000),
     (-+00),
     (--+-),
     (+-00),
     (+--+),
     (-0+-),
     (+0-+),
     (0-+-),
     (0+-+),
     (++-+),
     (-+-+),
     (+-+-),
     (-++-)}

For further examples on elementary vectors, solvability of linear inequality systems, sign vectors and oriented matroids, see `<https://marcusaichmayr.github.io/elementary_vectors/>`_.

Chemical reaction networks
~~~~~~~~~~~~~~~~~~~~~~~~~~

Robustness
**********

Given is a chemical reaction network involving five complexes.
To examine robustness of CBE, we compute the covectors corresponding to the resulting subspaces::

    sage: from sign_vectors.oriented_matroids import *
    sage: S = matrix([[-1, -1, 1, 0, 0], [0, 0, -1, 1, 0], [-1, 0, 0, 0, 1]])
    sage: S
    [-1 -1  1  0  0]
    [ 0  0 -1  1  0]
    [-1  0  0  0  1]
    sage: covectors_from_matrix(S)
    {(00000),
     (-+00-),
     (+-00+),
     (0---0),
     (-+++-),
     (-0---),
     (+0+++),
     (0+++0),
     (+---+),
     (+++++),
     (-----),
     (-+---),
     (+-+++)}
    sage: var('a, b, c')
    (a, b, c)
    sage: St = matrix([[-a, -b, 1, 0, 0], [c, 0, -1, 1, 0], [-1, 0, 0, 0, 1]])
    sage: St
    [-a -b  1  0  0]
    [ c  0 -1  1  0]
    [-1  0  0  0  1]
    sage: covectors_from_matrix(St(a=2, b=1, c=1))
    {(00000),
     (+-+0+),
     (0---0),
     (+-0-+),
     (-+0+-),
     (0+++0),
     (-+++-),
     (-0---),
     (+0+++),
     (+-+-+),
     (-+-0-),
     (+---+),
     (-+-+-),
     (+++++),
     (-----),
     (-+---),
     (+-+++)}

For :math:`a = 2`, :math:`b = 1` and :math:`c = 1`, the covectors of :math:`S` are included in the closure of the covectors of :math:`\widetilde{S}`.
To consider the general case, we compute the maximal minors of :math:`S` and :math:`\widetilde{S}`::

    sage: W  = matrix([[1, 0, 1, 1, 1], [0, 1, 1, 1, 0]])
    sage: W
    [1 0 1 1 1]
    [0 1 1 1 0]
    sage: var('a, b, c')
    (a, b, c)
    sage: Wt = matrix([[1, 0, a, a - c, 1], [0, 1, b, b, 0]])
    sage: Wt
    [    1     0     a a - c     1]
    [    0     1     b     b     0]
    sage: from sign_vector_conditions import *
    sage: condition_closure_minors(W, Wt)
    [{a - c > 0, b > 0, a > 0}]

Hence, the network has a unique positive CBE if and only if :math:`a, b > 0` and :math:`a > c`.

Uniqueness
**********
::

    sage: condition_uniqueness_minors(W, Wt)
    [{a - c >= 0, a >= 0, b >= 0}]

Hence, positive CBE are unique if and only if :math:`a, b \geq 0` and :math:`a \geq c`.

Unique existence of CBE
***********************

Now, we consider Example 20 from [MHR19]_.
Here, we have a parameter :math:`a > 0`.
Depending on this parameter, the network has a unique positive CBE::

    sage: var('a')
    a
    sage: W = matrix(3, 6, [0, 0, 1, 1, -1, 0, 1, -1, 0, 0, 0, -1, 0, 0, 1, -1, 0, 0])
    sage: W
    [ 0  0  1  1 -1  0]
    [ 1 -1  0  0  0 -1]
    [ 0  0  1 -1  0  0]
    sage: Wt = matrix(3, 6, [1, 1, 0, 0, -1, a, 1, -1, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0])
    sage: Wt
    [ 1  1  0  0 -1  a]
    [ 1 -1  0  0  0  0]
    [ 0  0  1 -1  0  0]

The first two conditions depend on the sign vectors corresponding
to the rows of these matrices which are independent of the specific value for :math:`a`::

    sage: assume(a > 0)
    sage: condition_uniqueness_sign_vectors(W, Wt)
    True

Hence, there exists at most one equilibrium.
Also the face condition is satisfied::

    sage: condition_faces(W, Wt)
    True

For specific values of ``a``, the pair of subspaces
determined by kernels of the matrices is nondegenerate.
This is exactly the case for :math:`a \in (0, 1) \cup (1, 2)`.
We demonstrate this for specific values::

    sage: condition_nondegenerate(W, Wt(a=1/2))
    True
    sage: condition_nondegenerate(W, Wt(a=3/2))
    True
    sage: condition_nondegenerate(W, Wt(a=1))
    False
    sage: condition_nondegenerate(W, Wt(a=2))
    False
    sage: condition_nondegenerate(W, Wt(a=3))
    False
"""
