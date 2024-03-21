r"""
Examples for ICMS 2024.

A SageMath Package for Elementary and Sign Vectors with Applications to Chemical Reaction Networks
--------------------------------------------------------------------------------------------------

Here are the up-to-date examples appearing
in TODO (link to extended abstract)
for `ICMS 2024 <https://icms-conference.org/2024/>`_.

Elementary vectors
~~~~~~~~~~~~~~~~~~

Functions dealing with elementary vectors, solvability of linear inequality systems and oriented matroids
are implemented in the package `elementary_vectors <https://github.com/MarcusAichmayr/elementary_vectors>`_.

We compute elementary vectors, using maximal minors::

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

We state linear inequality systems as intersection of a vector space and a Cartesian product of intervals.
To represent these objects, we use a matrix and a list of intervals::

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

Here, we give further details to the chemical reaction network appearing in the extended abstract.
The incidence and source matrix are given by::

    sage: I_E = matrix(5, 6, [-1, 1, 0, 1, 0, 0, 1, -1, -1, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 1, -1])
    sage: I_E
    [-1  1  0  1  0  0]
    [ 1 -1 -1  0  0  0]
    [ 0  0  1 -1  0  0]
    [ 0  0  0  0 -1  1]
    [ 0  0  0  0  1 -1]
    sage: I_S = matrix(5, 6, [1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1])
    sage: I_S
    [1 0 0 0 0 0]
    [0 1 1 0 0 0]
    [0 0 0 1 0 0]
    [0 0 0 0 1 0]
    [0 0 0 0 0 1]

By introducing reaction rates, we obtain the Laplacian matrix::

    sage: var('k12, k21, k23, k31, k45, k54')
    (k12, k21, k23, k31, k45, k54)
    sage: k = [k12, k21, k23, k31, k45, k54]
    sage: A_k = I_E * diagonal_matrix(k) * I_S.T
    sage: A_k
    [      -k12        k21        k31          0          0]
    [       k12 -k21 - k23          0          0          0]
    [         0        k23       -k31          0          0]
    [         0          0          0       -k45        k54]
    [         0          0          0        k45       -k54]

We have the following matrices of stoichiometric and kinetic-order complexes::

    sage: Y = matrix(5, 5, [1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]).T
    sage: Y
    [1 0 0 1 0]
    [1 0 0 0 0]
    [0 1 0 0 0]
    [0 0 1 0 0]
    [0 0 0 0 1]
    sage: var('a,b,c')
    (a, b, c)
    sage: Yt = matrix(5, 5, [a, b, 0, 0, 0, 0, 0, 1, 0, 0, c, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]).T
    sage: Yt
    [a 0 c 1 0]
    [b 0 0 0 0]
    [0 1 0 0 0]
    [0 0 1 0 0]
    [0 0 0 0 1]

The associated ODE system for the concentrations :math:`x` is given by::

    sage: var('x1, x2, x3, x4, x5')
    (x1, x2, x3, x4, x5)
    sage: x = vector([x1, x2, x3, x4, x5])
    sage: x_Yt = vector(prod(xi^yi for xi, yi in zip(x, y)) for y in Yt.columns())
    sage: x_Yt
    (x1^a*x2^b, x3, x1^c*x4, x1, x5)
    sage: Y * A_k * x_Yt
    (-k12*x1^a*x2^b + k31*x1^c*x4 - k45*x1 + k21*x3 + k54*x5, -k12*x1^a*x2^b + k31*x1^c*x4 + k21*x3, k12*x1^a*x2^b - (k21 + k23)*x3, -k31*x1^c*x4 + k23*x3, k45*x1 - k54*x5)

To study CBE, we consider the column spaces of the following matrices::

    sage: Y * I_E
    [-1  1  0  1 -1  1]
    [-1  1  0  1  0  0]
    [ 1 -1 -1  0  0  0]
    [ 0  0  1 -1  0  0]
    [ 0  0  0  0  1 -1]
    sage: Yt * I_E
    [   -a     a     c a - c    -1     1]
    [   -b     b     0     b     0     0]
    [    1    -1    -1     0     0     0]
    [    0     0     1    -1     0     0]
    [    0     0     0     0     1    -1]

We write the column spaces as matrices::

    sage: S = matrix([[-1, -1, 1, 0, 0], [0, 0, -1, 1, 0], [-1, 0, 0, 0, 1]])
    sage: S
    [-1 -1  1  0  0]
    [ 0  0 -1  1  0]
    [-1  0  0  0  1]
    sage: St = matrix([[-a, -b, 1, 0, 0], [c, 0, -1, 1, 0], [-1, 0, 0, 0, 1]])
    sage: St
    [-a -b  1  0  0]
    [ c  0 -1  1  0]
    [-1  0  0  0  1]

By computing the sign vectors of these matrices, we can investigate existence and uniqueness of CBE.
Several sign vector conditions for chemical reaction networks are implemented
in the package `sign_vector_conditions <https://github.com/MarcusAichmayr/elementary_vectors>`_.

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
    sage: covectors_from_matrix(S, kernel=False)
    {(00000),
     (00+-0),
     (0-0+-),
     (+000-),
     (+-++-),
     (0--+-),
     (+0+--),
     (+++--),
     (0-+--),
     (---+-),
     (--+0-),
     (--+--),
     (++0-0),
     (+++-+),
     (++-0-),
     (++---),
     (+-+--),
     (--0+0),
     (--0++),
     (+0-+-),
     (0+-0+),
     (+++-0),
     (+--+-),
     (---+0),
     (+-+0-),
     (++--0),
     (0+0-+),
     (0+--+),
     (++--+),
     (--+-0),
     (0-++-),
     (--++-),
     (-000+),
     (++0--),
     (--0+-),
     (++-0+),
     (0++-+),
     (-0+-+),
     (0+-++),
     (-+--+),
     (-++-+),
     (++-00),
     (00-+0),
     (++-+0),
     (++-+-),
     (-+-++),
     (--+00),
     (--++0),
     (+-0+-),
     (--+-+),
     (++0-+),
     (--+++),
     (-0-++),
     (-+-0+),
     (0-+0-),
     (++-++),
     (--+0+),
     (---++),
     (-+0-+)}
    sage: var('a, b, c')
    (a, b, c)
    sage: St = matrix([[-a, -b, 1, 0, 0], [c, 0, -1, 1, 0], [-1, 0, 0, 0, 1]])
    sage: St
    [-a -b  1  0  0]
    [ c  0 -1  1  0]
    [-1  0  0  0  1]
    sage: covectors_from_matrix(St(a=2, b=1, c=1), kernel=False)
    {(00000),
     (+-++-),
     (0--+-),
     (+0+--),
     (+++--),
     (0-0+-),
     (0-+--),
     (---+-),
     (00-++),
     (+000-),
     (+0-+0),
     (--+--),
     (0++-0),
     (+++-0),
     (+++-+),
     (0++--),
     (++---),
     (0--+0),
     (+-+--),
     (--0+0),
     (--0++),
     (-++--),
     (++0-0),
     (++-0-),
     (+0-+-),
     (0+-0+),
     (--+0-),
     (+--+-),
     (---+0),
     (+-+0-),
     (++--+),
     (--+-0),
     (0-++-),
     (++--0),
     (00+--),
     (--++-),
     (-++-0),
     (-000+),
     (++0--),
     (--0+-),
     (++-0+),
     (+0-++),
     (0++-+),
     (0+--+),
     (+--+0),
     (0+0-+),
     (0+-++),
     (-0+-+),
     (-++-+),
     (-+--+),
     (++-+-),
     (-+-++),
     (+--++),
     (++-00),
     (--+00),
     (--++0),
     (+-0+-),
     (--+0+),
     (--+-+),
     (--+++),
     (++0-+),
     (-0-++),
     (-+-0+),
     (0-+0-),
     (++-++),
     (0--++),
     (---++),
     (++-+0),
     (-0+-0),
     (-+0-+),
     (-0+--)}

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
    sage: condition_closure_minors(W, Wt) # random order
    [{a > 0, b > 0, a - c > 0}]

Hence, the network has a unique positive CBE if and only if :math:`a, b > 0` and :math:`a > c`.

Uniqueness
**********

We can also use the maximal minors to study uniqueness of CBE::

    sage: condition_uniqueness_minors(W, Wt) # random order
    [{a >= 0, b >= 0, a - c >= 0}]

Hence, positive CBE are unique if and only if :math:`a, b \geq 0` and :math:`a \geq c`.

Unique existence of CBE
***********************

Now, we consider Example 20 from [MHR19]_.
Here, we have a parameter :math:`a > 0`.
Depending on this parameter, the network has a unique positive CBE::

    sage: var('a')
    a
    sage: assume(a > 0)
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
