# Changes

## Unreleased

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

### Piecewise-linear univariates on the whole real line; shared subdivisions in polyhedral images

`mc::PWLU` also holds a continuous piecewise-linear function on the whole real line, such as an edge
of a Kolmogorov-Arnold network.  A constructor from breakpoints, values and end slopes, and
`PWLU::relu` from a sum of ReLU terms, build one whose first and last segments carry the end slopes;
`relu` gives each segment the sum of the slopes of the terms active on it, rather than a difference of
values, so that breakpoints close to one another cost no accuracy.
`val(x, true)`, `slope(x)` and `restrict(xL, xU)` extend it linearly beyond its domain; `restrict`
returns it on any interval, with its exact range there by `l()` and `u()`; `hull(under)` gives its
convex or concave envelope, `bkpts()` its breakpoints where the slope changes, and `compose(g)` the
exact composite of a piecewise-linear `g` with the estimator, with the crossings of the breakpoints
of `g` added.  The estimators of
the superposition models are unaffected: `val(x)` still throws outside the domain, and no class
changes layout.

`PolImg::add_semilinear_cuts` is split into `add_subdivision`, which subdivides a variable at given
breakpoints with continuous, binary or SOS2 auxiliary variables, and `add_subdivision_cut`, which
adds the cut of one function from its values at the knots kept, so that several functions can share
one subdivision.  Both forms of `add_semilinear_cuts` now call them.  The cuts are unchanged on the
univariate relaxations with breakpoints and on the piecewise-linear superposition relaxation of an
`FFMLP` (identical output for every encoding), except that a subdivision left on a variable by an
earlier call is now discarded, rather than reused, when none of the new breakpoints is kept.

### Piecewise-linear univariates as external operations

`ffpwlu.hpp` adds `mc::FFPWLU`, an external operation for one or several continuous
piecewise-linear functions of the same operand, each an `mc::PWLU` extended linearly beyond its
domain (see above).  In interval arithmetic the ranges are exact (`PWLU::restrict`), and the
backward propagation contracts the operand to the hull of the preimage of each range.  The McCormick
relaxations are the convex and concave envelopes on the range of the operand (`PWLU::hull`).
Derivatives are right derivatives at a breakpoint; symbolic ones are sums of `fstep`.  In a
polyhedral image each function is encoded exactly on the range of its operand
(`FFPWLU::Options::ENCODING`): by the convex hull of its graph, with cuts only, or on a subdivision
of the operand at its breakpoints (`PolImg::add_subdivision`), with continuous, binary or SOS2
auxiliary variables, as set by `PolImg::Options::BREAKPOINT_TYPE` by default (`DEFAULT`).  The
functions of an operation share one subdivision, at the union of their breakpoints, unless
`FFPWLU::Options::SEPARATE` is set: one set of auxiliary variables for all of them, and with `CONT`
the convex hull of their graphs jointly, where `HULL` gives the convex hull of each.  The image
leaves out of a subdivision the breakpoints within its tolerances (`BREAKPOINT_ATOL`,
`BREAKPOINT_RTOL`) of another one or of a bound; a function with a breakpoint left out is not its
interpolant on the subdivision, and is enclosed between two piecewise-linear functions there, so
that the relaxation remains valid (and by its convex hull where no breakpoint is kept).  In a
superposition model the
result is exact where the operand is a single-dependency estimator whose under- and overestimators
coincide, a variable of the model say: a `PWLU` summand is the composite of the function with the
operand's (`PWLU::compose`), a `PWCU` summand holds the exact range of the function on each cell.  As
after any composition, the summands are then reduced to `SupModel::Options::MAX_SUBDIV` segments,
which leaves them exact as long as the composite has no more segments (`SupVar::mod()` is new, to
read the options).  Other operands, and Chebyshev models, go through the abs-normal form
`y_0 + s_0 (x - x_0) + sum_k ds_k (z_k + |z_k|)/2`, with `z_k = x - x_k`.  Each evaluation is also there for a single function, as
the static `FFPWLU::value` (`slope` for the derivative), which `mc::KAN` uses for its edges.  `test/DAG/DAG_pwlu.cpp` checks every
evaluation and the new `PWLU` methods (35 checks); `test/RELAX/LINREL_pwlu.cpp` checks the polyhedral
encodings: every cut holds on the graph of the functions, the HULL cuts give the envelopes, and, with
Gurobi, the relaxations at a fixed operand give the function (BIN, SOS2) or its envelopes (HULL, CONT).

### Kolmogorov-Arnold networks

`ffkan.hpp` adds `mc::KAN`, a Kolmogorov-Arnold network: layers whose outputs are the sums over their
inputs of edge functions of one variable.  It mirrors `mc::MLP`: `append_data` (a layer given as its
edges), `set_data`, `read_data`; evaluation in the arithmetics of MC++, each edge with
`FFPWLU::value` -- exact range of each edge in interval arithmetic, envelopes in McCormick
arithmetic, exact summands of a variable in superposition models; Jacobian (`grad`, accumulated
forward or backward through the layers), directional derivative (`tang`) and adjoint (`adj`);
expression tree (`dag`), with one `FFPWLU` operation per input of a layer and one `FFLin` per output,
through which the polyhedral image and the backward propagation (`reval`) go.  The encoding of the
edges in a polyhedral relaxation is set by `KAN::Options::ENCODING` and `SEPARATE`.  The edges of a
layer are of one type; the type available is piecewise linear (`mc::PWLU`), and `KAN::Layer` is where
others -- splines, polynomials -- would be added.  The external operation on a `mc::KAN` in a DAG
is `mc::FFKAN` (see the entry on `mc::FFMAP` below).  On a network with widths 2-5-5-5-1 and
520 segments: 0.4 us for a value, 0.8 us for the Jacobian, 12 us for interval bounds.
`test/ML/test_KAN.cpp` checks the above (57 checks with the entry below).

### Networks saved as dictionaries by torch.save, read without Torch

`mc::NETIO::read_dict` reads a network saved by `torch.save` as a dictionary -- `model_type`,
`config`, `state_dict`, `normalization` -- with a reader of the zip archive and of the pickle format
of its own (`NETIO::read_pickle`, which returns the whole content as a tree of `NETIO::PValue`), as
`read_onnx` does for ONNX files: a Kolmogorov-Arnold network with edges `c + v u + sum_r a_r max(w_r
u + b_r, 0)` (`KANmodel`), read by `KAN::read_data`, or a dense network (`ANNmodel`), read by
`MLP::read_data`, which tells such a file from a TorchScript Module by its content
(`NETIO::is_torch_dict`).  The min-max scalings of the file are folded into the first and last
layers, so that the network read maps raw inputs to raw outputs; the layers in the scaled coordinates
of the file and the scalings are those returned by `read_dict`, for `set_data`.  The arguments of
`MLP::read_data` are unchanged.  The activation of the hidden layers of a dense network is the one
named in `config` (`activation`), else the one that matches the bounds held in the file, else the
one given to `read_dict`.  Files that are not such dictionaries, hold tensors of other shapes, or
are compressed or truncated are refused with the reason.  `test/ML/test_KAN_dict.py` generates the
files of the test and their values and bounds.

### FFMAP: networks as external operations, FFMLP and FFKAN

`ffmap.hpp` holds `mc::FFMAP<T, MAP>`, the external operation on a network held by an object of
class `MAP`, and `FFGradMAP`, `FFTanMAP`, `FFAdjMAP` for its Jacobian, directional derivative and
adjoint: the classes `FFMLP`, `FFGradMLP`, `FFTanMLP` and `FFAdjMLP` of `ffmlp.hpp`, which used
nothing of `mc::MLP` beyond `nin`, `nout`, `nwei`, `uid`, `eval`, `reval`, `grad`, `tang`, `adj` and
the options `TANGENT` and `ADJOINT`, with the network as a template argument.  `mc::FFMLP<T>` is now
the alias `FFMAP<T, MLP<T>>` (likewise `FFGradMLP`, `FFTanMLP`, `FFAdjMLP`), and `mc::FFKAN<T>` is
`FFMAP<T, KAN<T>>` (`FFGradKAN`, `FFTanKAN`, `FFAdjKAN`): a Kolmogorov-Arnold network is inserted
in a DAG, differentiated, relaxed and propagated backward as a multilayer perceptron is, with the
same policies for the network that the DAG holds.  Code on `FFMLP` is unchanged, with two
provisos: the network of the DAG is returned by `map()` -- `MLP()` remains for an `mc::MLP`, and
`pMLP()`, deprecated -- and a forward declaration `template <typename T> class FFMLP;` no longer
matches the alias.  The name of the operations is given by the static `tag()` of the network
("MLP", "KAN").  The polyhedral relaxation of an `FFKAN` is by McCormick relaxations by default, as
for any external operation; with `Options::RELAX` set to `AUX`, it goes through the expression tree
of the network, with the encodings of its edges.  `test/ML/test_FFKAN.cpp` checks the operation
against the network (14 checks).

`mc::MLPIO` (`mlpio.hpp`), which now reads Kolmogorov-Arnold networks as well, is renamed
`mc::NETIO` (`netio.hpp`).

## 5.2.0 (2026-11-18)

A minor version, by the rule INSTALL.md's "Releasing" states: classes shared with dependent modules changed.
`FFOp` has four more virtual functions (the tangent and adjoint hooks below), and `FFMLP` and `FFDAGEXT`
changed layout (1632 to 1168 bytes, now derived from `FFHOME`).  A module built on `pymcpp` (CRONOS's
`cronos`) must be rebuilt against 5.2 and require `pymcpp~=5.2.0`.  Code using `MLPREG` follows its new
interface (below).  The ABI fingerprint for 5.2 is recorded in `tools/abi/`.

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
  decoder of its own (`mc::MLPIO`, `mlpio.hpp`) and no other library, so that the pymcpp wheels, built
  without Torch, read trained networks too.  The graph must be a chain from its input to its output: `Gemm`
  (`transA` = 0, any `alpha`, `beta`, `transB`), or `MatMul` and `Add`, for each dense layer, with the weights
  as initializers or `Constant` nodes in float, double or float16; `BatchNormalization` after a dense layer,
  folded into it; `Relu`, `Tanh`, `Sigmoid`; and `Identity`, `Dropout`, `Flatten`, `Reshape`, `Squeeze`,
  `Unsqueeze`, `Cast`.  Anything else -- another node, a skip connection, a second input, weights stored
  outside the file -- is refused with the reason, and the MLP is left as it was.  This reads the exports of
  PyTorch (`torch.onnx.export`, Gemm or MatMul and Add), scikit-learn (skl2onnx) and Keras (tf2onnx), whose
  values agree with ONNX Runtime to its precision.  A network read from ONNX is evaluated by MC++
  (`EVALTORCH` has no effect).  The reader of TorchScript Modules moved to `MLPIO` as well, unchanged.

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
(`ffunc.hpp`).

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
