# MC++ Toolkit for Construction, Manipulation and Evaluation of Factorable Functions #

MC++ provides a collection of classes to support the construction and manipulation of factorable functions and their evaluation using a range of arithmetics. It is primarily written in C++ to promote execution speed, while the module pyMC includes Python binders through [pybind11](https://pybind11.readthedocs.io/en/stable/).

Expression trees in MC++ can comprise any finite combination of unary and binary operations from a default library. MC++ also features a mechanism to extend expression trees with external operations, currently including affine and polynomial subexpressions as well as multi-layer perceptrons (MLP). Through the library [CRONOS](https://github.com/omega-icl/cronos), MC++ can also be extended with systems of algebraic and differential equations. Expression trees generated with MC++ can be used by the library [CANON](https://github.com/omega-icl/canon) for local and global numerical optimization; and the library [MAGNUS](https://github.com/omega-icl/magnus) for the development and analysis of mathematical models.

The present version 4.0 of MC++ has capability for:

* Construction and differentiation of expression trees (both forward and reverse accumulation modes)
* Recursive decomposition of factorable expressions into linear/polynomial subexpressions and transcendental operations
* Elimination of variables in factorable expressions
* Decomposition of sparse polynomial expressions into quadratic forms
* Detection of reduction constraints in sparse polynomial expressions

The main bounding components in version 4.0 of MC++ include:

* Interval arithmetic
* Eigenvalue arithmetic
* Ellipsoidal arithmetic
* McCormick relaxations
* Taylor and Chebyshev models
* Polyhedral relaxations
* Superposition relaxations

A range of python notebooks are provided in `src/pymc` to illustrate these capabilities.

---
### Setting up MC++ ###

Refer to [INSTALL.md](./INSTALL.md) for instructions.

### Contacts ###

* Repo owner: [Benoit C. Chachuat](https://profiles.imperial.ac.uk/b.chachuat)
* OMEGA Research Group

---
### References ###

* Bompadre, A., A. Mitsos, [Convergence rate of McCormick relaxations](http://dx.doi.org/10.1007/s10898-011-9685-2), _Journal of Global Optimization_ **52**(1), 1-28, 2012
* Bompadre, A., A. Mitsos, B. Chachuat, [Convergence analysis of Taylor models and McCormick-Taylor models](http://dx.doi.org/10.1007/s10898-012-9998-9), _Journal of Global Optimization_, **57**(1), 75-114, 2013
* Bongartz, D., J. Najman, S. Sass, A. Mitsos, [MAiNGO – McCormick-based Algorithm for mixed-integer Nonlinear Global Optimization](http://permalink.avt.rwth-aachen.de/?id=729717), _Technical Report_, Process Systems Engineering (AVT.SVT), RWTH Aachen University, 2018
* Chachuat, B, B. Houska, R. Paulen, N. Peric, J. Rajyaguru, M.E. Villanueva, [Set-theoretic approaches in analysis, estimation and control of nonlinear systems](http://dx.doi.org/10.1016/j.ifacol.2015.09.097), _IFAC-PapersOnLine_, **48**(8), 981-995, 2015
* Karia, T., C.S. Adjiman, B. Chachuat, [Assessment of a two-step approach for global optimization of mixed-integer polynomial programs using quadratic reformulation](https://doi.org/10.1016/j.compchemeng.2022.107909), _Computers & Chemical Engineering_, **165**, 107909, 2022
* Mitsos, A., B. Chachuat, P.I. Barton, [McCormick-based relaxations of algorithms](http://dx.doi.org/10.1137/080717341), _SIAM Journal on Optimization_, **20**(2):573-601, 2009
* Rajyaguru, J., Villanueva M.E., Houska B., Chachuat B., [Chebyshev model arithmetic for factorable functions](http://dx.doi.org/10.1007/s10898-016-0474-9), _Journal of Global Optimization_, *68*, 413-438, 2017
* Tsoukalas, A., A. Mitsos, [Multi-variate McCormick relaxations](https://doi.org/10.1007/s10898-014-0176-0), _Journal of Global Optimization_, **59**(2), 633-662, 2014
* Villanueva, M.E., [Set-Theoretic Methods for Analysis, Estimation and Control of Nonlinear Systems](https://doi.org/10.25560/32528), PhD Thesis, Department of Chemical Engineering, Imperial College London, 2016
* Villanueva, M.E., J. Rajyaguru, B. Houska, B. Chachuat, [Ellipsoidal arithmetic for multivariate systems](https://doi.org/10.1016/B978-0-444-63578-5.50123-7), _Computer Aided Chemical Engineering_, **37**, 767-772, 2015
* Wechsung, A., P.I. Barton, [Global optimization of bounded factorable functions with discontinuities](http://dx.doi.org/10.1007/s10898-013-0060-3), _Journal of Global Optimization_, **58**(1), 1-30, 2014

