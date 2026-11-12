# Changes

## 5.1.1 (2026-11-12)

A patch release, by the rule INSTALL.md's "Releasing" states: no class shared with dependent modules changed
layout and pybind11 is still 3.0.4, so the 5.1 ABI fingerprint holds (`check_abi.py`: "MC++ 5.1.1 keeps the 5.1
fingerprint") and a module built against 5.1.0 -- CRONOS's `cronos` -- needs no rebuild.

### Wheels for Python 3.15; macOS 15 as the floor

Wheels are built for CPython 3.10 to 3.15 (cibuildwheel 4.3).  Under 3.15 the interpreter aborted at exit
once `pymcpp` had been imported: the `EllipsoidError` and `EllImgError` types were held in function-local
`static py::exception` objects, destroyed after the interpreter had finalized.  They are held as released
handles now, as pybind11's `docs/advanced/exceptions.rst` prescribes.

GitHub retired its macOS 14 runner image on 2026-11-02, so the macOS wheels are built on `macos-15` and
`MACOSX_DEPLOYMENT_TARGET` follows to 15.0: the Homebrew Boost and Armadillo bottles target the runner's OS,
and delocate refuses a wheel that claims less.  **The macOS wheels need macOS 15 or later.**  On macOS 14,
`pip install pymcpp` picks the 5.1.1 sdist over the 5.1.0 wheel and builds from source, which needs Boost and
Armadillo; `pip install --prefer-binary pymcpp`, or `pymcpp==5.1.0`, keeps the wheel.

`wheels.yml` filtered its push and pull-request triggers on `branches: [version-5.0]`, which does not exist,
so only tags and manual runs ever built wheels.  It names `version-5` now.

### Variable-length arrays removed

`CModel::_interpolation`, `TModel::_scale_bernstein` and both variants of `PolImg::_pwmccormick_cuts`
declared arrays sized at run time: a GCC/Clang extension, not C++, which `cl.exe` rejects -- the reason the
Windows wheels are pinned to clang-cl.  In `_pwmccormick_cuts` they held `PolVar<T>` (168 bytes with
`Interval`), as many as the product of the two variables' breakpoint counts, on the stack, where an overflow
is silent.  They are `std::vector`s now, and the library and its tests compile clean under `-Werror=vla`.
The external operation in `test/DAG/DAG_external.cpp` and `DAG_doxygen.cpp`, which the Doxygen pages
reproduce, did the same in `deriv()` and is fixed likewise.  Results are unchanged: the tests that reach
these paths (`CM-1D`, `CM-2D`, `TM-1D`, `TM-2D`, `SCM-Bernstein`, `LINREL_1d`, `LINREL_2d`, `PWLREL_1d`,
`DAG_external`) print byte-identical output.

### Compiler warnings

`override` on the overriding members of `FFPartial`, `FFIntegral` and `FFEval` (`ocbase.hpp`); the unused
`fact` in `fabs(CVar)` (`cmodel.hpp`) removed; the misleading indentation in `BASE_OC::eval_lagrange`
corrected.  No change in behaviour: the `RESCALE` scaling already ran once per row, after the inner loop, as
the indentation now shows.

### The release gate, CI and documentation

`check` no longer runs the `notebook/*.py` examples: they print rather than assert.  `check.yml` runs the gate
in CI on every push and pull request, followed by the JOSS paper's notebook, whose last cell asserts every
figure the paper quotes.  `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md` and `CITATION.cff` are added, the JOSS paper
sources are under `paper/`, and the Sphinx documentation is titled PyMC++.

## 5.1.0 (2026-10-26)

A minor version, by the rule INSTALL.md's "Releasing" states: a class shared with dependent modules changed
layout.  `FFSubgraph` gained two members, so it -- and `FFDAGEXT` and `FFMLP`, which hold one by value --
are not binary-compatible with 5.0.x, and a module built on `pymcpp` (CRONOS's `cronos`) must be rebuilt
against 5.1.  The ABI fingerprint for 5.1 is recorded in `tools/abi/`.

### The evaluation paths restore the value slots they found

Every evaluation of a subgraph -- `eval`, `SDFAD`, `SDBAD`, `TAD`, `compose`, `substitute` -- binds each
operation output's value slot to a working entry for the sweep and used to NULL the slots on exit.  A slot
can legitimately hold something before the sweep: a solver embedded as an external operation whose setup
evaluates the model's subgraph inside the enclosing evaluation, with the model's own variables as the
operation's inputs.  The nested sweep nulled slots the outer sweep was still reading through: a segfault.
`FFSubgraph::save_val()` records the slots as a sweep finds them and `unset_wk()` puts them back; all
seven paths use them, on the exception path too.

### `FFSubgraph` from Python

`pymcpp.FFSubgraph` is bound: `FFBase.subgraph()` returns one, it can be held, inspected (`len`, indexing
and iteration over its `FFOp`s, `l_op`, `len_tap`, `len_wrk`) and passed to `FFGraph.eval`, `reval` and
`veval`, so an expression evaluated repeatedly is traversed once.  The subgraph keeps its DAG alive.
`notebook/dag_incidence.py` reads a DAG's incidence pattern off one.

### The release gate

`cmake --build build --target check` runs, in order and stopping at the first failure: the version
written in three places agreeing (`tools/ci/check_version.py` -- `doc/MC++.dox` had stayed at 5.0.3
through 5.0.4), the ABI fingerprint against the recorded one, the Python checks (`src/pymcpp/check_*.py`,
starting with `check_ffgraph.py`: construction, evaluation, the subgraph as an object, the slot guard
from Python, forward against backward derivatives), the Python examples, and Doxygen with no warning
allowed.  Fourteen documentation defects that the last step found are fixed.

### Sparse differentiation: a correctness fix and the cost brought down to the nonzeros

`SBAD` truncated the adjoint of a dependent to its seed, skipping what the node's consumers hand back.
Right only for a dependent nothing else reads: a dependent consumed by another dependent's expression had
a wrong row (`{s, t}` with `t = 2s + x_0 x_1` gave `(2, 1)` for `(6, 9)`), and a dependent that is itself an
independent lost its column from every other row (`{s, x_1}`).  Both silent.  The seed now starts the
accumulation.  `check_ffgraph.py` section E gates both shapes against `SFAD`.

`SDFAD` and `SDBAD` swept every operation of the subgraph once per independent, or once per dependent:
`n x |ops|`, quadratic where the Jacobian's nonzero count is linear.  Each column of `SDFAD` now sweeps
the independent's forward cone (through `FFVar::opuse`), each row of `SDBAD` the dependent's backward
cone, in topological order within the subgraph.  On a tridiagonal system with `nnz = 3n`, measured by
`test/DAG/DAG_sparsead_scaling`: `SFAD` from `n^2.16` to `n^1.25` (298 ms to 5.8 ms at `n = 800`), `SBAD`
from `n^2.49` to `n^1.19` (1024 ms to 5.4 ms); a dense row from 450 ms to 2.0 ms in `SBAD`.  The symbolic
entries and their order are unchanged (`DAG_differentiate` byte-identical).  What remains above linear
at large `n` is the cost of inserting each entry into a DAG that grows with it, which the driver now
separates out as microseconds per entry.

### Known

* On a dense row formed as a left-deep chain (`s += sqr(x_i)`), `SFAD` stays quadratic: the forward cone
  of `x_i` holds every later addition for one entry.  That is the forward/reverse asymmetry, and `SBAD`
  is the routine for that row.

## 5.0.4 (2026-10-08)

Type comparisons that hold across extension modules (`mc::type_cast`); the patch-release rule, the ABI
fingerprint and its CI check.

## 5.0.3 (2026-10-07)

License and third-party notices.
