# Changes

## 5.2.0 (2026-11-18)

A minor version, by the rule INSTALL.md's "Releasing" states: classes shared with dependent modules changed.
`FFOp` has four more virtual functions (the tangent and adjoint hooks below), and `FFMLP` and `FFDAGEXT`
changed layout (1632 to 1168 bytes, now derived from `FFHOME`); `FFMLP` is also another type, an alias of
`FFMAP` (below).  A module built on `pymcpp` (CRONOS's
`cronos`) must be rebuilt against 5.2 and require `pymcpp~=5.2.0`.  Code using `MLPREG` follows its new
interface (below).  The ABI fingerprint for 5.2 is recorded in `tools/abi/`.

### HSL removed

`mchsl.hpp` (`mc::FFHSL`, wrapping the HSL routines MC13, MC21 and MC33) is removed, with the CMake option
`ENABLE_HSL`, the macro `MC__USE_HSL` and the `FLAG_HSL`/`LIB_HSL` settings of `makeoptions.mk`: nothing in MC++
used it.  Its successor is CRONOS's `ffblock.hpp` (`mc::FFBLOCK`), which uses HSL under `CRONOS__WITH_HSL` and
SuiteSparse BTF otherwise.

### Constraint propagation through products with Boost intervals

The backward propagation of `FFGraph::reval` tightened each factor of `z = x*y` with `z/y` (and the divisor
of `z = x/y` with `x/z`).  When the result and the other operand both contain zero, the tightened operand is
unconstrained, since `y` may be zero.  `mc::Interval` throws on a divisor containing zero, and the operation
was skipped; Boost returns `[0,0]/[0,1] = [0,0]`, and the factor was forced to zero.  With Boost intervals,
constraint propagation could then report a feasible model as infeasible, e.g. with complementarity
constraints `mu*g = 0` and `mu` in `[0,1]` (CANON's first-order optimality cuts).  The tightening is now
skipped in that case, whatever the interval arithmetic; results with `mc::Interval` are unchanged.

### Hyperbolic tangent and cosine with PROFIL intervals

`mcprofil.hpp` passed `tanh` and `cosh` on to PROFIL, which evaluates them from `exp` in interval arithmetic,
`tanh` as `(e^x-e^-x)/(e^x+e^-x)`.  On an interval this overestimates without limit (`tanh([-2,2])` came out as
`[-26.8,26.8]`, `tanh([-400,400])` as `[-inf,inf]`) and gives NaN once `exp` overflows; `cosh` was up to twice
too wide and could go below one.  A network with `tanh` activations then had infinite or NaN bounds, and
a model built on it could be reported infeasible.  Both are now bounded from their values at points, still in
PROFIL's rounded arithmetic: the end points for `tanh`, kept within `[-1,1]`, and the points of least and
greatest magnitude for `cosh`.

### Constraint propagation through a linear sum in linear time

The backward propagation through `FFLin` summed the terms of the other operands afresh for each operand: a
number of interval operations quadratic in the number of operands.  The constraint propagation of an `MLP`
runs through such sums, about 127,000 interval operations a pass for three hidden layers of 40 neurons.  The
terms are now summed once, backwards then forwards, which gives the same bounds in a linear number of
operations.  In CANON's `ann_peaks` example, for the same search, constraint propagation takes 10 s instead
of 68 s with Boost intervals, and 1.9 s instead of 7.5 s with `mc::Interval`.

### Constraint propagation: a bound that is not a number is no information

`FFGraph::reval` intersected the ranges it propagates with whatever an operation returned.  A result with a
NaN bound -- an operation that fails in the arithmetic: infinity over infinity after an overflow, zero times
an unbounded range, `tanh` of a wide range with PROFIL before the fix above -- made the propagation report
infeasibility with Boost and PROFIL intervals, and put NaN bounds on the variables with `mc::Interval`.  Such
a result is now left unbounded in the forward pass and skipped in the backward pass, as when the operation
throws.  `FFOp::_inter` is the intersection to use in the `reval` of an external operation, as `FFLin` does.
Boost also returns such bounds for an operation undefined on a whole range, e.g. the square root of a
negative range: this alone no longer makes the propagation infeasible, as with `mc::Interval` and PROFIL,
where the operation throws.  Test: `test/DAG/DAG_reval_nan`.

### Cuts for squares and products of auxiliary variables in a polyhedral image

`PolImg::append_cuts_SQR` and `PolImg::append_cuts_TIMES` add the cuts of a square or a product between
variables of an image, outside the evaluation of a DAG.  They found their operands through the DAG variable,
so an operand that is an auxiliary variable of the image (`PolVar(PolImg*, range)`, e.g. a linear combination
of variables defined by a cut) failed an assertion, and was a null pointer dereference under `NDEBUG`.  Each
operand may now be a DAG variable, an auxiliary variable or a constant (`PolVar(double)`); a variable that
does not belong to the image throws `PolImg::Exceptions::ENVMIS`.  An auxiliary variable now holds the
subdivision of its own range: without it the semilinear cuts failed on such a variable with
`BREAKPOINT_TYPE` `CONT` or `BIN` (the default).  The cuts for DAG variables are unchanged, and no class
changes layout.  `test/RELAX/LINREL_aux.cpp` checks the new cuts on a grid of points.  CANON's cuts from the
eigenvalue decomposition of quadratic forms (`MINLPBND::Options::DCQUADCUTS`) need this.

### Functions that do not depend on the variables in a polyhedral image or a sparse lifting

`PolImg::generate_cuts` reached the subtree of each dependent through its DAG variable, which a constant
dependent -- the zero objective of a feasibility problem, say -- does not have: a null pointer dereference.
A dependent without DAG variable, a constant or an auxiliary variable of the image, is now skipped, as there
are no cuts to generate for it.  `SLiftEnv::insert_dag` asserted on the zero polynomial, which it now
inserts as the constant 0.  No class changes layout.  `test/RELAX/LINREL_constant.cpp` checks both (it
crashed before the change).  CANON's models without objective (`MINLPBND`, `MINLPREF`) rely on this.

### Multilayer perceptrons in the reverse mode of AD

`FFMLP` and `FFMLPREG` could not be evaluated in `BADType<double>` (`mc::B`, or `fadbad::B` under
`MC__USE_FADBAD`): a numeric reverse-mode gradient of a DAG holding either threw, so that CANON's `NLPSLV`
with `GRADMETH = BAD` stopped before its first iteration.  Each output of the network is now one node of the
tape, whose partial derivatives are the Jacobian of the network computed by `MLP::grad` or `MLPREG::grad` in
double precision, rather than one node per neuron.  `mc::B::_nary` makes a node with n operands, and
`mc::bad_node` does so for either reverse-mode type.  On the training loss of a network with 2 inputs, 3
layers of 40 and one output over 100 samples (3441 weights), the gradient takes 7 ms, against 53 ms for the
symbolic gradient `FFGraph::BAD` (345,000 operations, 1.8 s to build); the two agree to 2e-17.  No class
changes layout.  `test/ML/test_MLP_reverse.cpp` compares the reverse mode with the forward mode and the
symbolic gradient, also under `MC__USE_FADBAD`.

### Piecewise-linear univariates: `PWLU` on the real line, `FFPWLU`

* `mc::PWLU` also holds a continuous piecewise-linear function on the whole real line, extended beyond
  its end breakpoints by its end segments: a constructor from breakpoints, values and end slopes,
  `relu` (from a sum of ReLU terms), `val(x, true)`, `slope`, `restrict` (exact range on an interval),
  `hull` (convex or concave envelope), `bkpts` and `compose`.
* `mc::FFPWLU` (`ffpwlu.hpp`) is the external operation for such functions of one operand, with one
  result per function.  It is exact where the arithmetic allows: the range in interval arithmetic, the
  envelopes in McCormick arithmetic, and the summand of a variable in superposition models (reduced to
  `SupModel::Options::MAX_SUBDIV` segments); other operands and Chebyshev models go through the
  abs-normal form.  The static `FFPWLU::value` and `slope` evaluate a single function.
* In a polyhedral image, `FFPWLU::Options::ENCODING` selects the convex hull of each function by cuts
  (`HULL`), or a subdivision of the operand at the breakpoints with continuous, binary or SOS2
  variables (`CONT`, `BIN`, `SOS2`); `DEFAULT` follows `PolImg::Options::BREAKPOINT_TYPE`.  The
  functions of an operation share one subdivision, unless `SEPARATE` is set.  A function with a
  breakpoint that the image leaves out, within its tolerances of another one, is enclosed between two
  piecewise-linear functions on the subdivision, so that the relaxation remains valid.
* `PolImg::add_subdivision` and `add_subdivision_cut` build a subdivision once and express several
  functions on it; `add_semilinear_cuts` calls them, with the same cuts as before.
* Tests: `test/DAG/DAG_pwlu.cpp` (35 checks), `test/RELAX/LINREL_pwlu.cpp` (104 relaxations).

### Kolmogorov-Arnold networks: `KAN`

`mc::KAN` (`ffkan.hpp`) mirrors `mc::MLP` for networks whose layers sum edge functions of one variable:
`set_data`, `append_data` (a layer given as its edges) and `read_data`; evaluation in the arithmetics
of MC++, each edge as in `FFPWLU`; `grad`, `tang` and `adj`; `reval`; and an expression tree with one
`FFPWLU` per input of a layer and one `FFLin` per output, whose encoding in a polyhedral image is set
by `KAN::Options::ENCODING` and `SEPARATE`.  The edges are piecewise linear (`mc::PWLU`); `KAN::Layer`
is where other types -- splines, polynomials -- would be added.  On a network with widths 2-5-5-5-1
and 520 segments: 0.4 us for a value, 0.8 us for the Jacobian, 12 us for interval bounds.  Test:
`test/ML/test_KAN.cpp` (60 checks).

`test/ML/sim_NN.cpp` and `notebook/sim_NN.ipynb` simulate surrogates of the Ackley function of two
variables -- three Kolmogorov-Arnold networks and three dense networks, `test/ML/ackley2d_*.pt` -- as
`FFKAN` and `FFMLP` operations, against the function itself.

### Networks as external operations: `FFMAP`, with `FFMLP` and `FFKAN`

`mc::FFMAP<T, MAP>` (`ffmap.hpp`), with `FFGradMAP`, `FFTanMAP` and `FFAdjMAP`, is the external operation
on a network held by an object of class `MAP`: the former classes `FFMLP`, `FFGradMLP`, `FFTanMLP` and
`FFAdjMLP`, with the network as a template argument.  `mc::FFMLP<T>` is the alias `FFMAP<T, MLP<T>>`
and `mc::FFKAN<T>` is `FFMAP<T, KAN<T>>`, likewise for the derivative operations, so that a
Kolmogorov-Arnold network is inserted in a DAG, differentiated, relaxed and propagated backward as a
multilayer perceptron is.  The network of the DAG is returned by `map()`; `MLP()` remains for an
`mc::MLP`.  A forward declaration `template <typename T> class FFMLP;` no longer matches the alias.
With `Options::RELAX` set to `AUX`, the relaxation of an `FFKAN` goes through the encodings of its
edges; it is by McCormick relaxations by default.  Test: `test/ML/test_FFKAN.cpp` (14 checks).

### Network files: `NETIO`

`mc::NETIO` (`netio.hpp`) holds the readers of network files, with no other library: ONNX files and
TorchScript Modules as before, and now the dictionaries saved by `torch.save` (`model_type`, `config`,
`state_dict`, `normalization`), with a reader of the zip archive and of the pickle format of its own
(`read_pickle`, `read_dict`).  A `KANmodel` is read by `KAN::read_data` and an `ANNmodel` by
`MLP::read_data`, which tells a dictionary from a TorchScript Module by its content.  The min-max
scalings of the file are folded into the first and last layers, so that the network maps raw inputs to
raw outputs.  The activation of a dense network is the one named in `config`, else the one matching
the bounds held in the file, else the one given to `read_dict`.  Where the type, the sizes or the
scalings of the file are also wanted, the file is read once with `NETIO::read_dict` and the dictionary
given to `KAN::set_data` or `MLP::set_data`; `read_data` does both.  Test: `test/ML/test_KAN.cpp`, on
files made by `test/ML/test_KAN_dict.py`.

### Python: `KAN` and `FFKAN`

`pymcpp.KAN` and `pymcpp.FFKAN` mirror `MLP` and `FFMLP`: a network is given by its layers of `PWLU`
edges or read from a dictionary of `torch.save`, and inserted in a DAG by `FFKAN()(X, NN)`.  The
network that the DAG holds is `FFKAN.map` (also `FFKAN.KAN`); `FFMLP.map` is added beside `FFMLP.MLP`.
`PWLU` gains the constructor from breakpoints, values and end slopes, `relu`, `val(x, extend)`, `slope`
and `bkpts`.  `pymcpp.NETIO.read_dict` returns the dictionary of a file (`type`, `widths`, scalings),
which `KAN.set_data` and `MLP.set_data` accept.  `src/pymcpp/check_ffmap.py`, in the `check` target, runs the same checks on `FFMLP` and
`FFKAN` against their networks evaluated in plain Python.

### The page of `pymcpp` on PyPI

`README_PYPI.md` is the page of `pymcpp` on PyPI, written for Python users: installation, what the
module does class by class, an example, the notebooks.  `README.md` remains the page of the
repository, C++ first.

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
* `MLP::Options::RELU2ABS` and `SIG2EXP` changed after the expression tree of the network was built -- after
  a polyhedral relaxation or a constraint propagation, say -- were ignored by these, which kept the former
  forms of the ReLU and sigmoid activations while the other evaluations used the new ones.  The tree is built
  again when they change, in `MLP` and `MLPREG`; all the options of an `MLP` in a DAG can now be changed after
  the insertion (`Op.MLP().options`).
* **ONNX as a second front end.** `MLP::read_data` reads a file whose name ends with `.onnx` with a protobuf
  decoder of its own (`mc::NETIO`, `netio.hpp`) and no other library, so that the pymcpp wheels, built
  without Torch, read trained networks too.  The graph must be a chain from its input to its output: `Gemm`
  (`transA` = 0, any `alpha`, `beta`, `transB`), or `MatMul` and `Add`, for each dense layer, with the weights
  as initializers or `Constant` nodes in float, double or float16; `BatchNormalization` after a dense layer,
  folded into it; `Relu`, `Tanh`, `Sigmoid`; and `Identity`, `Dropout`, `Flatten`, `Reshape`, `Squeeze`,
  `Unsqueeze`, `Cast`.  Anything else -- another node, a skip connection, a second input, weights stored
  outside the file -- is refused with the reason, and the MLP is left as it was.  This reads the exports of
  PyTorch (`torch.onnx.export`, Gemm or MatMul and Add), scikit-learn (skl2onnx) and Keras (tf2onnx), whose
  values agree with ONNX Runtime to its precision.  A network read from ONNX is evaluated by MC++
  (`EVALTORCH` has no effect).  The reader of TorchScript Modules moved to `NETIO` as well, unchanged.

### Networks whose weights are variables: `MLPREG` and `FFMLPREG` mirror `MLP` and `FFMLP`

* **Wrong derivatives, silently.**  The Jacobian operation of `FFMLPREG` wrote its entries row-wise and
  `FFMLPREG::deriv` read them column-wise: with more than one output, `FAD`, `BAD` and their variants
  returned the derivatives of other outputs in other operands (off by 2.5 on a 2-output test network).  And
  `FFMLPREG` had no bounds or relaxations at all: an interval or McCormick evaluation threw.
* `FFMLPREG` derives from `FFEXTERN`, as `FFMLP`: interval bounds, McCormick, polyhedral (AUX), spectral,
  Chebyshev and superposition relaxations by the options of the handle, constraint propagation, one home per
  `MLPREG` and per DAG, `MLPREG()` by reference, options changed after the insertion.  Its differentiation
  inserts the Jacobian `FFGradMLPREG`, or with `MLPREG::Options::TANGENT` and `ADJOINT` (off by default, as
  for `MLP`) the directional derivative `FFTanMLPREG` and the adjoint `FFAdjMLPREG`, one pass through the
  network per sweep -- always the cheaper here, the Jacobian having `nout*(nin+nwei)` entries.
* `MLPREG::grad`, `tang` and `adj` work by passes through the layers: the Jacobian in the inputs and the
  weights costs one backward pass per output, in place of `FADType` with `nin+nwei` directions (a cost in
  the square of the number of weights).  `AUTODIFF` selects the mode in other arithmetics; its default is
  `B` for `MLPREG`.
* `MLPREG::read_data` reads the structure and the weights of a trained network from ONNX or TorchScript, as
  `MLP::read_data` does; `weights()` holds the weights read, in the ordering of the operands, and
  `data(w)` gives the layers of `MLP::set_data` for any weights -- a trained `MLPREG` becomes an `MLP`.
* Torch (`MLPREG::Options::EVALTORCH`, for a Module whose layers have biases): the weights of the operands
  are written into a copy of the Module per thread, and the Jacobian, tangent and adjoint come from Torch's
  backward passes; they agree with MC++ to 1e-15.  Writing the weights in costs as much as the pass, and MC++
  is the faster: 63 against 514 us for the Jacobian of 12,705 weights and 5 outputs, 6.7 against 15 ms for
  a million weights.
* **Interface changed** to that of `MLP`: an `MLPREG` is evaluated on its operands, the inputs then the
  weights -- `eval(valin, valout)` in place of `eval(y, x, w)`; `set_structure(nin, {{nneu, activ}, ...})`
  in place of the pair, and `append_layer`; `MLPREG::ACTIV_TYPE` in place of
  `MLPREG::Options::ACTIVTYPE` (same values); `nvar()`, `uid()`, `layers()`; `FFGRADMLPREG` is
  `FFGradMLPREG`.  The call operators of `FFMLPREG` keep their arguments, with a copy policy and vector
  forms added.
* Python: `MLPREG` and `FFMLPREG` are bound as `MLP` and `FFMLP` are (`FFMLPREG.MLPREG` and `options` after
  the insertion, `MLPREG.data` for `MLP.set_data`), and checked by `check_ffextern.py`.  The notebook
  `notebook/tutorial_mlp.ipynb` samples MATLAB's `peaks` at Sobol points, reads the network from ONNX as an
  `MLPREG`, trains its 3,441 weights with PyTorch's Adam on the loss and gradient of the MC++ DAG (`bdiff` with
  `ADJOINT`, `veval` over mini-batches), then evaluates, differentiates, bounds and relaxes the trained `MLP`
  in a DAG.

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
* The homes are a base class of their own, `mc::FFHOME<O, OPT>` (`ffhome.hpp`), from which `FFEXTERN`
  derives, so that the external operations of the libraries built on MC++ that hold an object get the
  same homes -- CANON's `FFNLPSLV` holds its NLP solver this way.  `ffhome_traits<O>::copy` copies an
  object that its copy constructor cannot (a solver that needs a DAG of its own for its copy, say), and
  `Home::aux` keeps what the operations of a class share on the object (a cache of solutions, the tree of
  a Jacobian).  `FFEXTERN<T, O>::Options` is `FFEXTERN_OPTIONS<T, O>`, now at namespace scope, and
  `ffextern_uid()` is `ffhome_uid()`.

### DAG evaluation

`FFGraph::veval` evaluated each scenario with the `FFGraph::eval` that builds, at every call, a hash map from
the variables to their values, with one allocation per variable: on a subgraph with thousands of variables -- a
network whose weights are variables -- this cost more than the evaluation itself.  Each thread of `veval` now
builds that map once and keeps it for all its scenarios, through the overload of `eval` that takes it
(`ndxVar`), which `veval` used only when `MC__VEVAL_USE_NDXVAR` was defined; the macro is gone.  For one
`FFMLPREG` operation on 3,443 variables, 300 to 130 us per scenario for the loss and its gradient in the
weights.  For repeated calls of `eval` with the same variables, pass `ndxVar` likewise.

### Compiler warnings

`-Wmaybe-uninitialized` from the default constructor of `mc::B<T>` (`bdiff.hpp`) wherever an `MLP` was
differentiated in interval or McCormick arithmetic, and from `inter()` of two `SCVar` (`scmodel.hpp`);
`-Walloc-size-larger-than=` from GCC 15 on the allocation in the backward sweep through linear operators
(`ffunc.hpp`); `-Wmaybe-uninitialized` from the constructors of a variable or constant in a DAG (`ffunc.hpp`),
which copied the object before two of its members were set.

### Tests

`test/DAG/DAG_tangent`, `DAG_adjoint`, `DAG_doxygen_hooks`, `DAG_dagext_home`; `test/ML/test_MLP_tangent`,
`test_MLP_adjoint`, `test_MLP_grad`, `test_MLP_copy`, `test_MLP_torch` (with `test_MLP_torch.py` and the
modules it generates), `test_MLP_ident`, `test_MLP_home`, `test_MLP_onnx` (with `test_MLP_onnx.py` and the ONNX
models it generates), `test_MLPREG` (rewritten), `test_MLPREG_torch`; `src/pymcpp/check_ffextern.py` in the release
gate.  `test_MLP`, `DAG_external`, `DAG_unittests` and `FFDAGEXT_2d` print the same output
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
