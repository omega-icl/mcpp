# MC++: Toolkit for Construction, Manipulation and Evaluation of Factorable Functions

MC++ provides a collection of classes to support the construction and manipulation of factorable functions and
their evaluation using a range of arithmetics.  It is primarily written in C++ to promote execution speed, while the
library PyMC++ provides Python binders through [pybind11](https://pybind11.readthedocs.io/en/stable/) (the module
is imported as `pymcpp`).

Expression trees in MC++ can comprise any finite combination of unary and binary operations from a default
library and are described by directed acyclic graphs (DAGs).  MC++ also features a mechanism to extend DAGs with
*external operations*, currently including affine and polynomial subexpressions, multi-layer perceptrons (MLP) and
nested DAGs.  Through the library [CRONOS](https://github.com/omega-icl/cronos), MC++ can also be extended with
systems of algebraic and differential equations.  Expression trees generated with MC++ are used by the library
[CANON](https://github.com/omega-icl/canon) for local and global numerical optimization, and by the library
[MAGNUS](https://github.com/omega-icl/magnus) for the development and analysis of mathematical models.

## Capabilities

Version 5 of MC++ supports:

* **Expression trees**: construction of factorable functions as DAGs, with common-subexpression elimination;
  differentiation in both forward and reverse accumulation modes; evaluation in any of the arithmetics below;
  vectorized and multi-threaded evaluation over many points.
* **Decomposition**: recursive decomposition of factorable expressions into linear/polynomial subexpressions and
  transcendental operations, the basis for tailored relaxation and reformulation strategies.
* **Nested expression trees**: a DAG as an operation of another DAG, to enable tailored bounding strategies.
* **External operations**: user-defined operations with their own evaluation and differentiation rules -- the
  mechanism CRONOS uses to embed differential-equation solvers.

The bounding arithmetics, each a header that can be used on its own or through the DAGs:

| arithmetic | purpose |
|---|---|
| Interval arithmetic | rigorous enclosures (Boost, PROFIL/BIAS or FILIB++ as backend) |
| Eigenvalue arithmetic | spectral bounds of Hessians |
| Ellipsoidal arithmetic | ellipsoidal enclosures of multivariate systems |
| McCormick relaxations | convex/concave relaxations and their subgradients |
| Taylor and Chebyshev models | polynomial models with rigorous remainder bounds |
| Polyhedral relaxations | linear relaxations for use in LP/MILP solvers |
| Superposition relaxations | relaxation via separable estimators |

A range of Python scripts and notebooks in `notebook/` illustrate these capabilities.

## Setting up MC++

PyMC++ installs from PyPI, with pre-built wheels for Python 3.10 to 3.14 on Linux, macOS and Windows:

```
pip install pymcpp
```

To use the C++ headers, a different interval backend, or the HSL and Torch options, build from source with CMake;
refer to [INSTALL.md](INSTALL.md) for the requirements, options and instructions.

## Layout

```
src/mc/              the MC++ headers (header-only library)
src/3rdparty/        bundled third-party headers (FADBAD++, used only with MC__USE_FADBAD)
src/pymcpp/          the Python binders (PyMC++, module `pymcpp`)
test/                C++ test and example programs
notebook/            Python scripts and Jupyter notebooks
```

## Contacts

* Repo owner: [Benoit C. Chachuat](https://profiles.imperial.ac.uk/b.chachuat)
* OMEGA Research Group, Imperial College London

## References

Methods implemented in MC++:

* Mitsos, A., B. Chachuat, P.I. Barton, [McCormick-based relaxations of algorithms](http://dx.doi.org/10.1137/080717341), *SIAM Journal on Optimization*, **20**(2):573-601, 2009
* Bompadre, A., A. Mitsos, [Convergence rate of McCormick relaxations](http://dx.doi.org/10.1007/s10898-011-9685-2), *Journal of Global Optimization* **52**(1), 1-28, 2012
* Bompadre, A., A. Mitsos, B. Chachuat, [Convergence analysis of Taylor models and McCormick-Taylor models](http://dx.doi.org/10.1007/s10898-012-9998-9), *Journal of Global Optimization*, **57**(1), 75-114, 2013
* Tsoukalas, A., A. Mitsos, [Multi-variate McCormick relaxations](https://doi.org/10.1007/s10898-014-0176-0), *Journal of Global Optimization*, **59**(2), 633-662, 2014
* Wechsung, A., P.I. Barton, [Global optimization of bounded factorable functions with discontinuities](http://dx.doi.org/10.1007/s10898-013-0060-3), *Journal of Global Optimization*, **58**(1), 1-30, 2014
* Villanueva, M.E., J. Rajyaguru, B. Houska, B. Chachuat, [Ellipsoidal arithmetic for multivariate systems](https://doi.org/10.1016/B978-0-444-63578-5.50123-7), *Computer Aided Chemical Engineering*, **37**, 767-772, 2015
* Chachuat, B, B. Houska, R. Paulen, N. Peric, J. Rajyaguru, M.E. Villanueva, [Set-theoretic approaches in analysis, estimation and control of nonlinear systems](http://dx.doi.org/10.1016/j.ifacol.2015.09.097), *IFAC-PapersOnLine*, **48**(8), 981-995, 2015
* Villanueva, M.E., [Set-Theoretic Methods for Analysis, Estimation and Control of Nonlinear Systems](https://doi.org/10.25560/32528), PhD Thesis, Department of Chemical Engineering, Imperial College London, 2016
* Rajyaguru, J., Villanueva M.E., Houska B., Chachuat B., [Chebyshev model arithmetic for factorable functions](http://dx.doi.org/10.1007/s10898-016-0474-9), *Journal of Global Optimization*, **68**, 413-438, 2017
* Karia, T., C.S. Adjiman, B. Chachuat, [Assessment of a two-step approach for global optimization of mixed-integer polynomial programs using quadratic reformulation](https://doi.org/10.1016/j.compchemeng.2022.107909), *Computers & Chemical Engineering*, **165**, 107909, 2022
* Zha, Y., M.E. Villanueva, B. Houska, B. Chachuat, [Relaxation via separable estimators: Arithmetic and implementation](https://doi.org/10.1007/s10898-026-01648-z), *Journal of Global Optimization*, in press, 2026

Software building on MC++:

* Bongartz, D., J. Najman, S. Sass, A. Mitsos, [MAiNGO -- McCormick-based Algorithm for mixed-integer Nonlinear Global Optimization](http://permalink.avt.rwth-aachen.de/?id=729717), *Technical Report*, Process Systems Engineering (AVT.SVT), RWTH Aachen University, 2018
* [Pyomo](https://pyomo.readthedocs.io/en/latest/explanation/solvers/mcpp.html), whose MC++ interface bounds factorable functions with McCormick relaxations (built on an earlier MC++ release)

## License

MC++ is published under the Eclipse Public License.  The bundled FADBAD++ (`MC__USE_FADBAD`, off by default) is
a separate work distributed for non-commercial use only: commercial use requires a license from its authors.
