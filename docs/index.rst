PyMC++
======

PyMC++ (``pymcpp``) is the Python interface to `MC++
<https://github.com/omega-icl/mcpp>`_, a C++ library for constructing,
manipulating and bounding factorable functions. It provides:

- **Expression DAGs** (``FFGraph``, ``FFVar``): build computational graphs of
  factorable functions, evaluate them and differentiate them symbolically
  (forward and backward automatic differentiation).
- **Set arithmetics** for rigorous bounding: intervals (``Interval``),
  McCormick relaxations (``McCormick``), dense and sparse Taylor/Chebyshev
  models (``TModel``, ``CModel``, ``SCModel``), spectral bounds (``Specbnd``),
  interval superposition models and polyhedral relaxations (``PolImg``).
- **Machine-learning models in DAGs**: multilayer perceptrons imported from
  PyTorch (``MLP``, ``FFMLP``), for evaluation, differentiation and bounding
  (requires a torch-enabled build, see below).

Why PyMC++
----------

Deterministic global optimization, reachability analysis of dynamic systems and
rigorous model validation all rest on the same primitive: given a factorable
function and a host set for its variables, produce a guaranteed enclosure of its
range, or a convex relaxation a solver can minimize. These are neither point
evaluations nor derivatives but set-valued computations, propagated operation by
operation through the expression, each step carrying a correctness guarantee.

Automatic-differentiation tools and modelling languages build expression graphs
and differentiate them efficiently, but none of them propagates set-valued
information. The libraries that do are specialized to one arithmetic apiece, and
the global optimization solvers that construct relaxations keep them inside
their own branch-and-bound loop, where a user cannot call them on an arbitrary
expression. PyMC++ exposes several such arithmetics over one shared expression
graph, composable and extensible, from Python.

It is intended for researchers and practitioners in global and robust
optimization, optimal control and estimation, reachability and verification, and
anyone who needs guaranteed bounds rather than sampled ones -- including on
models that embed trained neural networks. MC++ has underpinned this kind of
work for over a decade: `MAiNGO
<https://git.rwth-aachen.de/avt-svt/public/maingo>`_ is built on it and `Pyomo
<https://pyomo.readthedocs.io/en/latest/explanation/solvers/mcpp.html>`_ wraps
it.

Installation
------------

.. code-block:: bash

   pip install pymcpp

Prebuilt wheels cover Linux (x86_64) and macOS (arm64) for Python 3.10+.

To use the PyTorch features, pymcpp must be compiled locally against your
installed torch (binary wheels are torch-free, since libtorch is only
ABI-compatible with the exact torch version it was built against):

.. code-block:: bash

   pip install torch
   pip install scikit-build-core pybind11 pybind11-stubgen ninja cmake
   pip install pymcpp --no-binary pymcpp --no-build-isolation -C cmake.define.ENABLE_TORCH=ON

A local build requires a C++17 compiler plus Boost, Armadillo and a
BLAS/LAPACK on the system.

Contents
--------

.. toctree::
   :maxdepth: 1

   quickstart
   tutorials
   api
