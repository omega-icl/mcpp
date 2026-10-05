# Installing MC++

## What MC++ is

MC++ is a toolkit for the construction, manipulation and evaluation of factorable functions.  Expression trees
are described by directed acyclic graphs (DAGs) and can comprise any finite combination of unary and binary
operations from a default library; the DAGs can be extended with external operations -- affine and polynomial
subexpressions, multi-layer perceptrons, nested DAGs, and through [CRONOS](https://github.com/omega-icl/cronos)
systems of algebraic and differential equations.  Version 5 provides the construction and differentiation of
expression trees (forward and reverse accumulation), their recursive decomposition into linear/polynomial
subexpressions and transcendental operations, nested expression trees, and a range of bounding arithmetics:
intervals, eigenvalues, ellipsoids, McCormick relaxations, Taylor and Chebyshev models, polyhedral and
superposition relaxations.

MC++ is a header-only C++ library in the `mc::` namespace.  Its Python binders form the library PyMC++, built
with [pybind11](https://pybind11.readthedocs.io/) and imported as the module `pymcpp`.  Expression trees built with MC++ are used by
[CANON](https://github.com/omega-icl/canon) for local and global optimisation and by
[MAGNUS](https://github.com/omega-icl/magnus) for model development and analysis.

## Option A: PyMC++ from PyPI

```
pip install pymcpp
```

installs pre-built wheels for Python 3.10 to 3.14 on Linux x86-64 (glibc >= 2.27), macOS 14+ (Apple silicon)
and Windows x86-64, with no further dependency.  The notebooks in `notebook/` run as they are.  Build from
source when you need the C++ headers, a different interval backend, the HSL or Torch options, or a `pymcpp`
built against the same headers and pybind11 release as another module that exchanges objects with it (CRONOS's
`cronos`, for instance).

## Option B: building from source with CMake

### Requirements

| dependency | version | role |
|---|---|---|
| C++ compiler | C++20 | GCC, Clang or MSVC (the PyPI wheels are built on Linux, macOS and Windows) |
| CMake | >= 3.18 | |
| Boost | | the default interval arithmetic (`MC_INTERVAL_LIBRARY=BOOST`) |
| BLAS / LAPACK | | |
| Armadillo | | dense linear algebra |
| PROFIL/BIAS or FILIB++ | optional | alternative interval backends (`PROFIL_HOME`, `FILIB_HOME`) |
| HSL MC13 / MC21 / MC33 | optional | structural analysis routines (`ENABLE_HSL`, needs a Fortran runtime) |
| PyTorch (libtorch, CPU) | optional | the Torch interface (`ENABLE_TORCH`) |
| Python | >= 3.10, development headers | for PyMC++ (`pymcpp`) |
| pybind11 | **>= 3.0.3, not 3.1.0** (3.0.4 recommended) | the binders |
| pybind11-stubgen | optional | generates the `pymcpp.pyi` type stub |

pybind11 is taken from the vendored `extern/pybind11` (a git submodule: `git submodule update --init`) if
present, otherwise from an installed pybind11 package.  Release 3.1.0 has a `keep_alive` regression on rejected
overloads that breaks the binders; 3.0.3 is the first release with the virtual-base detection and the
string-ownership fix the binders rely on.

### Configure

```
git clone --recurse-submodules https://github.com/omega-icl/mcpp.git
cmake -S mcpp -B build
```

Everything is found without help when installed in standard places.  The main options:

| option | default | |
|---|---|---|
| `MC_INTERVAL_LIBRARY` | `BOOST` | interval backend: `BOOST`, `PROFIL` (set `PROFIL_HOME`), `FILIB` (set `FILIB_HOME`), `NONVERIFIED` |
| `ENABLE_HSL` | OFF | link the HSL routines MC13, MC21, MC33 (found by `find_library`, else assumed on the linker path) |
| `ENABLE_TORCH` | OFF | the Torch interface; `TORCH_PYTHON_PREFIX` points at the torch package or its CMake directory, else torch is imported from the build's Python |
| `MC__USE_FADBAD` | OFF | FADBAD++ forward AD (`src/3rdparty/fadbad++`) in place of the default implementation -- FADBAD++ is restricted to non-commercial use, see Licence |
| `ENABLE_EXAMPLES` | OFF | build every program under `test/` (the FADBAD ones need `MC__USE_FADBAD`) |
| `PYMCPP_STUBS` | ON | generate and install `pymcpp.pyi` (skipped, with a warning, if `pybind11-stubgen` is not importable) |
| `CUSTOM_PYTHON_PATH` | | a Python executable or a virtual-environment root, to build against a specific Python |
| `CMAKE_BUILD_TYPE` | `Release` | `Release` is `-O2`; `Debug` is `-O0 -g` |
| `PYMCPP_INSTALL_DIR`, `MCPP_NOTEBOOK_INSTALL_DIR` | `lib`, `notebook` | installation subdirectories |

### Build, install, test

```
cmake --build build -j
cmake --install build --prefix <prefix>
```

installs the headers in `<prefix>/include`, the `pymcpp` module and its `.pyi` stub in `<prefix>/lib`, and the
notebooks and scripts in `<prefix>/notebook`.  Then, with `<prefix>/lib` on the Python path:

```python
import pymcpp
```

With `-DENABLE_EXAMPLES=ON`, every program under `test/` is built as an executable in the build directory.

For a wheel build with scikit-build-core, pass `-DPYMCPP_INSTALL_DIR=.` and install only the `python_modules`
component.

### Using the headers in another project

MC++ v5 installs its headers but exports no CMake package yet.  A project using them adds `<prefix>/include`
(or `<MC++>/src/mc` from the source tree, plus `src/3rdparty/fadbad++` under `MC__USE_FADBAD`) to its include
path, and defines the macros its MC++ build used: `MC__USE_THREAD`, `MC__USE_TADIFF`, `MC__USE_ARMADILLO`, the
interval backend (`MC__USE_BOOST`, `MC__USE_PROFIL` or `MC__USE_FILIB`), and `MC__USE_HSL`, `MC__USE_TORCH` or
`MC__USE_FADBAD` when enabled.  CRONOS does this through its `MCPP_ROOT` variable.

### Modules that exchange objects with `pymcpp`

A module built on `pymcpp`'s types (CRONOS's `cronos`, for one) receives `FFVar`, `FFGraph`, ... objects from it.
pybind11 shares types between modules only when both were built against the same internals, so build the two
modules with **the same pybind11 release** and **the same MC++ headers**, and rebuild both after any header
change.  A mismatch typically shows as `free(): invalid pointer` or a segmentation fault the first time the two
modules exchange an object.

## Troubleshooting

| symptom | cause |
|---|---|
| `extern/pybind11` is empty | the submodule is not checked out: `git submodule update --init`, or clone with `--recurse-submodules` |
| `pybind11-stubgen` warning at configure time | `pip install pybind11-stubgen` into the build's Python, or `-DPYMCPP_STUBS=OFF` |
| `ENABLE_TORCH=ON, but torch could not be imported` | install `torch` into the build's Python, or set `TORCH_PYTHON_PREFIX` |
| `PROFIL_HOME not set` / `FILIB_HOME not set` | the chosen interval backend's prefix must be given (cache or environment variable) |
| crash when another module exchanges objects with `pymcpp` | the two modules were built from different MC++ headers or pybind11 releases: rebuild both from clean |
| the Python docstrings or stub look stale after a header change | rebuild `pymcpp` (the stub is regenerated after each build) |

## Licence

MC++ is published under the Eclipse Public License.

**FADBAD++** (`MC__USE_FADBAD=ON`, off by default) is a separate work by Ole Stauning and Claus Bendtsen,
distributed free of charge for non-commercial use only: commercial use requires a licence from its authors.  MC++
is built and distributed without it, with its own forward-AD implementation; anyone enabling FADBAD++ takes on
that restriction for the resulting binaries.
