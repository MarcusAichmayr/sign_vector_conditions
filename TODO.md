# TODO

* [x] cite MHR19 (once)

## functions (f_pol, f_exp)

* [x] add examples and explanations (like in thesis) to top of file
* Lists are no longer accepted as input.
  * [ ] adapt files in CoCalc
  * [ ] adapt in Master's thesis
* These functions should return something like a polynomial.
  e.g. `f_exp(W, Wt)` could return `e^x + 2 e` instead of `<function f_exp_pol.<locals>.f at 0x7f31b88f9940>`

## utility

* [x] add examples

## injectivity

* [x] add examples and explanations (like in thesis) to top of file
* [ ] add examples
    * (optional) `cond_inj_intersection`
    * [ ] `max_minors_prod`
    * [ ] `compare_all`
    * [ ] `leq`
    * [ ] `geq_leq` add more examples
    * [ ] `cond_inj_minors`
* [ ] remove random tests
* [x] support latex statements in documentation for inline math
  (`a, b \in \mathcal{R}`, `\leq`, `\geq`)

## bijectivity

* [x] add examples and explanations (like in thesis) to top of file
* [ ] add examples
    * (optional) `cond_nondegenerate`, (but not for `nondegenerate`)
    * [ ] `nondeg_cond1` add explanations
    * [ ] `equal_components`
    * [ ] `find_vector`
    * [ ] `nondeg_cond2`

## robust bijectivity

* [x] add examples and explanations (like in thesis) to top of file
* [ ] add examples
    * (optional) `cond_closure_sign_vectors`
    * (optional) `cond_closure_minors`
* [ ] remove random tests
* [ ] add further results, perturbations in coefficients and exponents

