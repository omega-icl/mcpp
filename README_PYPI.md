<p align="center">
  <img src="https://raw.githubusercontent.com/omega-icl/mcpp/version-5/doc/mcpp-logo-256.png" alt="MC++ logo" width="160">
</p>

# pymcpp: factorable functions, their derivatives and their bounds in Python

`pymcpp` (PyMC++) is the Python distribution of [MC++](https://github.com/omega-icl/mcpp).  It builds factorable
functions as expression graphs, differentiates them, and evaluates them in a range of arithmetics: real numbers,
but also intervals, convex relaxations, polynomial models and polyhedral relaxations, which enclose the range of
a function over a set instead of sampling it.

These set-valued computations are what deterministic global optimisation, reachability analysis and rigorous
model validation rest on.  `pymcpp` exposes several of them over one shared expression graph, composable and
extensible, including on models that embed trained neural networks.  The computations are compiled in C++.

## Installation

```
pip install pymcpp
```

Wheels are built for CPython 3.10 to 3.15 on Linux (x86-64, glibc 2.27 or later), macOS 15 or later (Apple
silicon) and Windows (x86-64), with no further dependency.

The wheels are built without PyTorch.  Networks are read from ONNX files and from dictionaries saved by
`torch.save` with no other library; reading TorchScript modules, or evaluating a network by PyTorch, needs a
build from source with the Torch option.  See
[INSTALL.md](https://github.com/omega-icl/mcpp/blob/version-5/INSTALL.md) for that, for the other build options,
and for the C++ headers.

## What it does

* **Expression graphs** (`FFGraph`, `FFVar`): factorable functions as directed acyclic graphs, with
  common-subexpression elimination; differentiation in forward and reverse modes, built on the graph (`fdiff`,
  `bdiff`); evaluation in any of the arithmetics below, at one point or vectorised over many (`eval`, `veval`).
* **Interval bounds** (`Interval`) and **spectral bounds** of Hessians (`Specbnd`).
* **Convex and concave relaxations** with their subgradients (`McCormick`).
* **Polynomial models** with rigorous remainder bounds: Taylor and Chebyshev (`TModel`, `CModel`), and their
  sparse variants (`SCModel`, `SICModel`).
* **Superposition relaxations**, by sums of piecewise-constant or piecewise-linear functions of one variable
  (`PWCSModel`, `PWLSModel`).
* **Ellipsoidal enclosures** (`Ellipsoid`, `EllImg`).
* **Polyhedral relaxations** (`PolImg`): linear cuts, with auxiliary and binary variables where asked for, to
  hand to an LP or MILP solver.
* **Structure and reformulation**: dependence and sparsity (`FFDep`), sparse polynomials (`SPoly`), lifting into
  polynomial and univariate subexpressions (`SLift`).
* **External operations**, with their own evaluation, differentiation and relaxation rules: affine expressions
  (`FFLin`), a graph nested in another (`FFDAGEXT`), operations written in Python (`FFCustom`).
* **Neural networks as operations of a graph**: multilayer perceptrons (`MLP`, `FFMLP`), read from ONNX or
  PyTorch files; networks whose weights are variables, to train on the graph (`MLPREG`, `FFMLPREG`); and
  Kolmogorov-Arnold networks with piecewise-linear edges (`KAN`, `FFKAN`).  Each is evaluated, differentiated,
  bounded and relaxed like any other operation.

Every class and option carries its documentation, and a `.pyi` stub serves editors.

## An example

One function, on one graph, evaluated in five ways:

```python
import pymcpp as mc

# f(x, y) = exp(x y) - 2 x^2 + 3 on a graph
G = mc.FFGraph()
x, y = mc.FFVar(G, "x"), mc.FFVar(G, "y")
f = mc.exp(x * y) - 2 * x**2 + 3

# Value and gradient at a point; the gradient is built on the graph by reverse differentiation
_, _, df = G.bdiff([f], [x, y])
print(G.eval([f], [x, y], [1., 0.5]))        # [2.6487...]
print(G.eval(df, [x, y], [1., 0.5]))         # [-3.1756..., 1.6487...]

# Guaranteed bounds on the range over the box [0,1] x [0,2]
box = [mc.Interval(0., 1.), mc.Interval(0., 2.)]
print(G.eval([f], [x, y], box))              # [2 : 10.389...]

# Convex and concave relaxations at the centre of the box, with their subgradients
mx = [mc.McCormick(b, 0.5 * (b.l + b.u)).sub(2, i) for i, b in enumerate(box)]
r = G.eval([f], [x, y], mx)[0]
print(r.cv, r.cc)                            # 3.0 6.6945...

# Polyhedral relaxation, as linear cuts for an LP or MILP solver
img = mc.PolImg()
p = [mc.PolVar(img, v, b) for v, b in zip([x, y], box)]
img.generate_cuts(G.eval([f], [x, y], p))
print(len(img.cuts), "cuts")                 # 19 cuts
```

## Tutorials

Jupyter notebooks in the repository's [notebook](https://github.com/omega-icl/mcpp/tree/version-5/notebook)
folder, one per module:

* [ffunc](https://github.com/omega-icl/mcpp/blob/version-5/notebook/ffunc.ipynb) -- expression graphs:
  construction, differentiation, evaluation;
* [ffdep](https://github.com/omega-icl/mcpp/blob/version-5/notebook/ffdep.ipynb) -- dependence and sparsity;
* [interval](https://github.com/omega-icl/mcpp/blob/version-5/notebook/interval.ipynb),
  [mccormick](https://github.com/omega-icl/mcpp/blob/version-5/notebook/mccormick.ipynb),
  [specbnd](https://github.com/omega-icl/mcpp/blob/version-5/notebook/specbnd.ipynb) -- interval bounds,
  McCormick relaxations, spectral bounds;
* [cmodel](https://github.com/omega-icl/mcpp/blob/version-5/notebook/cmodel.ipynb),
  [scmodel](https://github.com/omega-icl/mcpp/blob/version-5/notebook/scmodel.ipynb) -- Chebyshev models, dense
  and sparse;
* [supmodel](https://github.com/omega-icl/mcpp/blob/version-5/notebook/supmodel.ipynb) -- superposition
  relaxations;
* [polimage](https://github.com/omega-icl/mcpp/blob/version-5/notebook/polimage.ipynb),
  [ellimage](https://github.com/omega-icl/mcpp/blob/version-5/notebook/ellimage.ipynb) -- polyhedral
  relaxations, ellipsoidal enclosures;
* [ffmlp](https://github.com/omega-icl/mcpp/blob/version-5/notebook/ffmlp.ipynb),
  [tutorial_mlp](https://github.com/omega-icl/mcpp/blob/version-5/notebook/tutorial_mlp.ipynb),
  [sim_NN](https://github.com/omega-icl/mcpp/blob/version-5/notebook/sim_NN.ipynb) -- neural networks in a
  graph: a multilayer perceptron as an operation, training one on the graph, and simulating Kolmogorov-Arnold
  and dense surrogates of a test function.

[joss_paper](https://github.com/omega-icl/mcpp/blob/version-5/notebook/joss_paper.ipynb) gathers the listings of
the paper describing the library, and
[runs in Colab](https://colab.research.google.com/github/omega-icl/mcpp/blob/version-5/notebook/joss_paper.ipynb)
with no installation.

## Related packages

[CRONOS](https://github.com/omega-icl/cronos) builds on `pymcpp` graphs to simulate ODE, DAE and PDE models and
compute their sensitivities, and embeds the solvers as operations of a graph.
[CANON](https://github.com/omega-icl/canon) uses the graphs and their relaxations for local and global
optimisation.

## Documentation and source

The documentation of the Python module is at [pymcpp.readthedocs.io](https://pymcpp.readthedocs.io).  MC++ is
also a header-only C++ library; the source, its documentation, the release notes
([CHANGES.md](https://github.com/omega-icl/mcpp/blob/version-5/CHANGES.md)) and the issue tracker are on
[GitHub](https://github.com/omega-icl/mcpp).

## Citing pymcpp

If MC++ or PyMC++ supports your research, please cite it.  A paper describing the library is under review at the
[Journal of Open Source Software](https://joss.theoj.org/); until it appears, cite the software release:

```bibtex
@software{mcpp,
  author  = {Chachuat, Beno{\^i}t and Ruth, Kai and Sandrin, Marco},
  title   = {{MC++}: a toolkit for the construction, manipulation and bounding of factorable functions},
  version = {5.2.0},
  year    = {2026},
  url     = {https://github.com/omega-icl/mcpp}
}
```

[`CITATION.cff`](https://github.com/omega-icl/mcpp/blob/version-5/CITATION.cff) carries the same metadata.  The
methods MC++ implements have papers of their own, listed under
[References](https://github.com/omega-icl/mcpp/blob/version-5/README.md#references) in the repository's README --
please cite the one your results rest on alongside the software.

## Licence

`pymcpp` is distributed under the Eclipse Public License 2.0, with the GNU General Public License, version 2 or
later, as a Secondary License.  The third-party components are listed in THIRD_PARTY_NOTICES.txt.  MC++ is
developed by Benoit Chachuat and the OMEGA Research Group at Imperial College London.
