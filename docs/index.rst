PyMC
====

PyMC (``pymcpp``) is the Python interface to `MC++
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
