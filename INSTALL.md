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
| `MC_ARMA_WRAPPER` | OFF | link Armadillo's runtime wrapper library instead of using Armadillo header-only with BLAS/LAPACK linked directly (the wrapper also pulls in whatever Armadillo was built with, e.g. ARPACK and MPI) |
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

Layouts depend on the **configuration** as well as on the version: `MC__USE_THREADLOCAL` changes `FFBase` (and so
`FFGraph`), `MC__USE_THREAD` changes `Vect`, `MC__USE_ARMADILLO` the Chebyshev and Taylor model classes, and the
interval library the template operations.  Build a module that exchanges objects with `pymcpp` with the definitions of
MC++'s `mc` CMake target (`MC__USE_THREAD MC__USE_TADIFF MC__USE_ARMADILLO`, the Boost interval library, not
`MC__USE_THREADLOCAL`), as `pymcpp`'s wheels are.

## Releasing

A module built on `pymcpp` (CRONOS's `cronos` requires `pymcpp~=5.0.4`: any 5.0.x from 5.0.4) is compiled against
MC++'s headers and then uses objects created by a separately compiled `pymcpp`.  That holds only if the two agree on
those classes' binary layout, hence the rule for **patch releases** (x.y.z to x.y.z+1):

* **allowed**: bug fixes inside functions, new free functions, new classes, documentation;
* **a new minor version** (x.y+1.0) for anything else that a dependent module can see: a member added, removed or
  reordered in a class shared with other modules (`FFBase`, `FFGraph`, `FFVar`, `FFNum`, `FFOp`, `FFSubgraph`, the
  operation classes `FFPartial`, `FFEval`, `FFIntegral`, `FFCustom`, `FFDAGEXT`, `FFLin`, `FFMLP`, `FFVect`, ...), a
  virtual function added or reordered in one, a change of the configuration macros above, or **another pybind11
  release** (pinned in `pyproject.toml`).

**The check.**  `tools/abi/abi_fingerprint.cpp` prints the layout of those classes (`sizeof`, `alignof`, whether
polymorphic), the configuration and the pybind11 version, built with `pymcpp`'s own definitions:

```bash
cmake --build build --target abi_fingerprint
python tools/abi/check_abi.py build/abi_fingerprint          # compares with tools/abi/fingerprint-X.Y-<platform>.txt
```

It passes when the fingerprint matches the one recorded for the minor version, and fails otherwise, showing what
changed.  The CI job `abi` runs it on every push, and publishing a release requires it.  At a new minor version,
record its fingerprint (`check_abi.py build/abi_fingerprint --record`) and commit the file.  Only Linux x86-64 is
recorded: a source change shows on any platform.  The fingerprint cannot see everything -- a virtual function
reordered without a size change, a member's type changed to one of the same size -- so the rule still applies where
the check is silent.

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

MC++ is published under the Eclipse Public License 2.0, with the GNU General Public License, version 2 or later,
as a Secondary License (see [LICENSE](LICENSE)).

**FADBAD++** (`MC__USE_FADBAD=ON`, off by default) is a separate work by Ole Stauning and Claus Bendtsen,
distributed free of charge for non-commercial use only: commercial use requires a licence from its authors.  MC++
is built and distributed without it, with its own forward-AD implementation; anyone enabling FADBAD++ takes on
that restriction for the resulting binaries.
