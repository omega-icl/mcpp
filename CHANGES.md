# Changes

## 5.2.0 (2026-11-18)

A minor version, by the rule INSTALL.md's "Releasing" states: classes shared with dependent modules changed.
`FFOp` has four more virtual functions (the tangent and adjoint hooks below), and `FFMLP` and `FFDAGEXT`
changed layout (1632 to 1152 bytes).  A module built on `pymcpp` (CRONOS's `cronos`) must be rebuilt against
5.2 and require `pymcpp~=5.2.0`.  The ABI fingerprint for 5.2 is recorded in `tools/abi/`.

### Tangent and adjoint hooks for external operations

An external operation could only give its partial derivatives, entry by entry, and the differentiation
formed the products with them: an operation whose natural derivative is a product -- a linear solve, a
sensitivity system, a pass through a network -- had to build its whole Jacobian.  `FFOp` gains two
optional hooks: `tangent()`, which inserts the directional derivative of the results along the tangents
of the operands, and `adjoint()`, which inserts the contribution of the operation to the adjoints of its
operands.  `has_tangent(nSweep)` and `has_adjoint(nSweep)` say whether to use them, given the number of
sweeps that will go through the operation, so that an operation can keep its Jacobian where that is the
cheaper.  The forward engines (`FAD`, `SFAD`, `DFAD`) and the backward ones (`BAD`, `SBAD`, `DBAD`) call
them; an operation without them is differentiated as before.  The Doxygen section on external operations
has a complete example, which is also `test/DAG/DAG_doxygen_hooks.cpp`.

`BAD` and its variants refused any subgraph holding a linear operator without an adjoint (`ocbase.hpp`).
Such an operator is now differentiated with respect to the independents in forward mode, inside the
backward sweep: in `DAG_linear_operator` the line "BAD d E/du (evolution, thru OpP)" goes from an
exception to the result `FAD` gives.

### Multilayer perceptrons

* `MLP::Options::TANGENT` and `ADJOINT` (off by default) have the differentiation of an `FFMLP` insert
  `FFTanMLP` and `FFAdjMLP` -- one pass through the network per sweep -- in place of the Jacobian
  `FFGradMLP`, where that is estimated to be the cheaper from the sizes of the network and the number of
  sweeps.  `MLP::tang` and `MLP::adj` are the passes.
* `MLP::grad` in double computes the Jacobian by passes through the layers -- forward, or one evaluation
  with the rows carried backward where the outputs are no more than the inputs -- instead of FADBAD++'s
  forward mode: 88 to 10 us for 2 inputs and 4 layers of 40, 6.9 to 0.14 ms for 200 inputs and 5 outputs,
  1.1 s to 8 ms for 500 inputs and layers of 1000.  `grad` in interval and other arithmetics is unchanged.
* Torch (`read_data`, `EVALTORCH`): a read that failed left the new module attached to the previous
  network; activations were matched by substring (`LeakyReLU` read as ReLU, `Tanhshrink` as Tanh); two
  `Linear` layers in a row had their activations shifted; a layer without bias was refused; a bias read
  before its weight wrote out of bounds.  All fixed, and a module is now checked against the network read
  from it at two points.  The parameters no longer accumulate gradients, the Jacobian takes one backward
  pass where the module takes batches, and modules are evaluated in double precision.
* `MLP` has an assignment operator (the implicit one was deleted).

### External operations: one home per object and per DAG

* An `FFMLP` or `FFDAGEXT` operation held its own copy of its network or expression tree, taken at the
  insertion, and was identified by the address of the object passed: a network changed and inserted again
  returned the operation of the former network, and **two different expression trees of the same shape on
  the same operands were taken for one operation** (`x0*x1` and `x0+x1` both evaluated to `x0*x1`).
  Silent.  An object now has one home per DAG -- the object, copied once (policy `COPY`) or the one given
  (`SHALLOW`, `TRANSFER`), and the options of the operations on it -- shared by all its operations,
  derivative operations included.  An object is found again by a stamp of its content (`MLP::uid()`,
  `DAGEXT::uid()`), so a changed object is another operation.
* The options can be changed after the insertion, and every operation of the home follows: the
  relaxation options through the handle's `options`, the network's or tree's own through `MLP()` or
  `DAGEXT()`, which return the object of the DAG by reference (`pMLP()` and `pDAGEXT()` are deprecated).
  The Jacobian of an expression tree is differentiated once per tree instead of at every call.
* Python: `FFMLP.options` and `FFDAGEXT.options` can be changed after the call; `FFMLP.MLP`,
  `FFDAGEXT.DAGEXT` and `attached` are new.  `options` is a reference member now: a binding of a class
  derived from `FFEXTERN` by `def_readwrite("options", ...)` must become a property.

### Compiler warnings

`-Wmaybe-uninitialized` from the default constructor of `mc::B<T>` (`bdiff.hpp`) wherever an `MLP` was
differentiated in interval or McCormick arithmetic, and from `inter()` of two `SCVar` (`scmodel.hpp`);
`-Walloc-size-larger-than=` from GCC 15 on the allocation in the backward sweep through linear operators
(`ffunc.hpp`).

### Tests

`test/DAG/DAG_tangent`, `DAG_adjoint`, `DAG_doxygen_hooks`, `DAG_dagext_home`; `test/ML/test_MLP_tangent`,
`test_MLP_adjoint`, `test_MLP_grad`, `test_MLP_copy`, `test_MLP_torch` (with `test_MLP_torch.py` and the
modules it generates), `test_MLP_ident`, `test_MLP_home`; `src/pymcpp/check_ffextern.py` in the release
gate.  `test_MLP`, `test_MLPREG`, `DAG_external`, `DAG_unittests` and `FFDAGEXT_2d` print the same output
as 5.1.1.

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
