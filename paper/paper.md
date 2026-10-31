---
title: 'MC++ and pyMC++: construction, manipulation and bounding of factorable functions in C++ and Python'
tags:
  - C++
  - Python
  - factorable functions
  - expression graphs
  - automatic differentiation
  - interval arithmetic
  - convex relaxation
  - global optimization
authors:
  - name: Benoit Chachuat
    orcid: 0000-0003-4780-9686
    corresponding: true
    affiliation: 1
  - name: Kai Ruth
    orcid: 0009-0002-4371-3484
    affiliation: "1, 2"
  - name: Marco Sandrin
    orcid: 0009-0007-0157-4152
    affiliation: 1
affiliations:
  - name: Sargent Centre for Process Systems Engineering, Department of Chemical Engineering, Imperial College London, United Kingdom
    index: 1
  - name: Department of Chemistry and Applied Biosciences, ETH Zurich, Switzerland
    index: 2
date: 28 October 2026
bibliography: paper.bib
---

# Summary

Some questions about a mathematical model cannot be settled by merely
evaluating it, however many times: the range of values an output can take over
a whole set of inputs, whether a constraint holds throughout a prescribed
domain, or provable bounds on the global optimum of a nonconvex problem. They
call for enclosures and convex relaxations, obtained by propagating sets rather
than numbers through the model's expressions.

MC++ is a C++ toolkit that does this for factorable functions, the class
obtained by composing arithmetic operations with standard intrinsic functions,
which covers most algebraic models in science and engineering. A function is
represented once, as a directed acyclic graph (DAG) of its operations, and that
graph is then evaluated in whichever arithmetic a task needs: real numbers,
derivatives, interval enclosures, convex and concave relaxations, polynomial
models with remainder bounds, and others. The companion library pyMC++
(module `pymcpp`) brings the toolkit to Python as prebuilt wheels, one
`pip install pymcpp` away.

# Statement of need

Deterministic global optimization, reachability analysis of dynamic systems and
rigorous model validation all rest on the same primitive: given a factorable
function and a host set for its variables, produce a guaranteed enclosure of its
range, or a convex relaxation a solver can minimize. These are set-valued
computations, propagated operation by operation through the expression with a
correctness guarantee at each step [@chachuat2015].

The audience is researchers and practitioners in global and robust
optimization, optimal control and estimation, and reachability and
verification, together with the developers of the solvers and modeling tools
that serve them. They have had to choose between implementing an arithmetic
themselves and working inside a solver that keeps its own private. What has been
missing is a library offering several such arithmetics over one shared
expression graph, composable and extensible, and reachable from Python, where
much of the modeling work that now needs rigorous bounds is done.

# State of the field

Two mature software ecosystems sit either side of this gap.

On one side are expression-graph and automatic-differentiation tools: CasADi
[@casadi], CppAD, ADOL-C [@adolc], and the algebraic modeling languages Pyomo
[@pyomo; @pyomobook], JuMP [@jump], AMPL and GAMS, together with the tracing
autodiff of PyTorch and JAX. All build expression graphs and differentiate them
efficiently: they answer what $f(x)$ and $\nabla f(x)$ are *here*. What values
$f$ can take over a set, or what convex function underestimates it there, is at
most the business of an extension; none offers a breadth of independently usable
and interchangeable set-valued arithmetics.

On the other side are rigorous-bounding libraries, each specialized to one
arithmetic and most of them ageing. Intervals are served by filib++ [@filib],
PROFIL/BIAS [@profilbias], Boost.Interval [@boostinterval] and IBEX [@ibex], the
last pairing them with contractor programming. Taylor models are served by COSY
Infinity, which is not openly licensed, by TaylorModels.jl [@taylormodels_jl;
@taylorseries_jl], and inside reachability tools such as Flow\* [@flowstar] and
CORA [@cora], where the arithmetic serves reachable sets of dynamic systems
rather than arbitrary expressions. Chebyshev-model arithmetic has, to our
knowledge, no open alternative.

Deterministic global optimization solvers --- BARON [@baron], ANTIGONE
[@antigone], SCIP [@scip], Couenne [@couenne], Alpine [@alpine], Gurobi
[@gurobi], EAGO [@eago] and MAiNGO [@maingo] --- do relax factorable functions,
but inside a branch-and-bound loop, and several are closed source.

The nearest reusable alternative, the Julia stack around EAGO.jl, McCormick.jl
and IntervalArithmetic.jl, offers fewer arithmetics as separate packages rather
than over one DAG, and no Python interface; there, the only prior route has been
Pyomo, which exposes an earlier MC++ release through a separately compiled
`contrib.mcpp` extension.

Contributing these arithmetics to an existing package was not a realistic
alternative to building a library. Several of the graph-building tools admit
user-defined value types, so an enclosure could in principle be propagated
through them; but none was designed to host a collection of independently
usable, interchangeable and composable bounding arithmetics, with one serving
as the coefficient type of another --- which is why Pyomo delegates its
bounding to MC++ rather than to its own expression system. The solvers do have
the machinery, but specialized to their own search and not addressable from
outside.

# Capability

Functions are defined on a DAG environment, with participating variables
connected by operations; expressions are built by writing ordinary arithmetic on
them. Identical subexpressions are detected and shared as they are created. For
instance, a DAG for the two factorable functions $f(x,y) = e^{xy} - 2x^2y + 3$
and $g(x,y) = \sqrt{x^2y + 4} + e^{xy}$ shares the subexpressions $xy$, $e^{xy}$
and $x^2$, and comprises 16 operations rather than 21. A *subgraph* --- the
operation list for a chosen set of dependents, in topological order, what an
automatic-differentiation tool would call a tape or a codelist --- is a
first-class object that can be held, inspected and re-used, so repeated
evaluation reuses one precomputed sequence instead of re-extracting it from the
graph.

A single evaluation routine sweeps a subgraph and returns its dependents; a
vectorized form handles many points at once over several threads, and a reverse
form contracts the variable ranges against bounds on the dependents by interval
propagation.

What that routine computes is decided by the type of the values it is given, not
by the function called. Floats evaluate $f$; intervals enclose its range over
the host set; and the relaxation types return convex and concave estimators
[@mccormick1976] with their subgradients, polynomial models with remainder
bounds, separable
superposition estimators, or the linear cuts of a polyhedral relaxation.
Spectral bounds on Hessians and ellipsoidal enclosures of multivariate systems
propagate the same way.

Derivatives are generated symbolically, in forward or reverse accumulation, as a
sparse triplet of *DAG nodes* that share subexpressions with the function they
differentiate. Being expressions in the same graph, the entries can be
differentiated again for Hessians, and evaluated in any of the arithmetics
above. Underneath lie MC++'s own operator-overloading differentiation types
templated on the evaluation type.

A DAG can also be rewritten, through recursive decomposition of the dependents
into sparse polynomial or rational subexpressions and transcendental ones, with
an auxiliary variable introduced for each lifted term --- the reformulation step
that precedes tailored relaxation or quadratic reformulation in global
optimization. The sparse polynomial arithmetic, in a monomial or a Chebyshev
basis, is also exposed directly.

# User interface

The listings below build the DAG of $f$ and $g$, evaluate and bound $f$,
differentiate $f$, and decompose both dependents. A DAG is an `FFGraph`, its
variables are `FFVar` objects connected by `FFOp` operations, and `subgraph`
returns the `FFSubgraph` that `eval`, `veval` and `reval` take; derivatives come
from `fdiff` or `bdiff`. Each listing was executed against `pymcpp` 5.1.0; the
values in the comments are what it returned.

The listings are Python. The C++ interface exposes the same classes and calls in
the `mc` namespace, with each arithmetic's template arguments written out where
Python has a single bound type, so it is in C++ that the composition described
in the next section is chosen.

```python
import pymcpp as mc

DAG = mc.FFGraph()
X, Y = DAG.add_var("X"), DAG.add_var("Y")
F = mc.exp(X*Y) - 2*X**2*Y + 3
G = mc.sqrt(X**2*Y + 4) + mc.exp(X*Y)
F.set("F"); G.set("G")

SGF = DAG.subgraph([F])                               # 11 operations
DAG.eval(SGF, [F], [X, Y], [0.5, 1.0])                # [4.148721270700128]
DAG.veval(SGF, [F], [X, Y], [[0.5, 1.], [1., 2.]])    # one entry per point
```

Bounding uses the very same call: only the type of the values changes. All of
the following describe $f$ over the host set $x \in [-1,1]$, $y \in [0,2]$.

```python
IX, IY = mc.Interval(-1., 1.), mc.Interval(0., 2.)
DAG.eval(SGF, [F], [X, Y], [IX, IY])          # [-0.8647, 10.3891]

xMC = mc.McCormick(IX, 0.5)                   # range [-1,1], point x = 0.5
yMC = mc.McCormick(IY, 1.0)                   # range  [0,2], point y = 1.0
xMC.sub(2, 0); yMC.sub(2, 1)                  # seed subgradients
[FMC] = DAG.eval(SGF, [F], [X, Y], [xMC, yMC])
FMC.cv, FMC.cc, FMC.cvsub                     # 2.0, 8.5756, [2.0, -1.0]

CM = mc.SCModel(3)                            # sparse Chebyshev model, order 3
[FCM] = DAG.eval(SGF, [F], [X, Y],
                 [mc.SCVar(CM, 0, IX), mc.SCVar(CM, 1, IY)])
FCM.B, FCM.R                                  # [-1.5230, 7.7258], [-0.6674, 0.6674]

SUP = mc.PWLSModel(2)                         # superposition relaxation,
[FSU] = DAG.eval(SGF, [F], [X, Y],            #   8 segments per variable
                 [mc.PWLSVar(SUP, 0, IX, 8), mc.PWLSVar(SUP, 1, IY, 8)])
FSU.l(), FSU.u()                              # -0.8647, 9.3891
FSU.uval({0: 0.5, 1: 1.0}), FSU.oval({0: 0.5, 1: 1.0})    # 2.9699, 4.8069

IMG = mc.PolImg()                             # polyhedral relaxation
XP, YP = mc.PolVar(IMG, X, IX, True), mc.PolVar(IMG, Y, IY, True)
[FP] = DAG.eval([F], [X, Y], [XP, YP])
IMG.generate_cuts([FP])                       # 23 cuts over 9 lifted variables
```

Where a McCormick relaxation carries the convex and concave relaxations of $f$
on the box evaluated at a given point, a Taylor or Chebyshev model and a
superposition relaxation describe their enclosure over the whole host set and
can also be evaluated pointwise. The subgradients from a McCormick relaxation
additionally supply the cut a bundle method or an LP relaxation needs, and the
cuts generated by a polyhedral relaxation are ready for an LP or MILP bounding
model.

Differentiation and decomposition are called on the DAG in the same way:

```python
rows, cols, dF = DAG.bdiff([F], [X, Y])   # sparse Jacobian as DAG nodes
[e.str() for e in dF]
# ['Y * EXP( X * Y ) + ( X * 2 ) * ( - Y ) * 2',
#  'X * EXP( X * Y ) + ( - SQR( X ) * 2 )']
DAG.eval(dF, [X, Y], [0.5, 1.0])          # [-0.35127872929987, 0.32436063535006]

env = mc.SLift(DAG)
env.process([F, G])
print(env)
# 4 AUXILIARY VARIABLES:
#   V2 <- X * Y      V3 <- EXP( X * Y )
#   V4 <- Y * SQR( X ) + 4            V5 <- SQRT( Y * SQR( X ) + 4 )
# 2 AUXILIARY POLYNOMIAL CONSTRAINTS:
#   0 = V2 - X * Y          0 = V4 - ( Y * SQR( X ) + 4 )
# 2 AUXILIARY TRANSCENDENTAL CONSTRAINTS:
#   0 = V3 - EXP( V2 )      0 = SQR( V5 ) - V4
# 2 DEPENDENT EXPRESSIONS:
#   V3 + 3 - ( Y * SQR( X ) ) * 2           V3 + V5
```

# Software design

MC++ is header-only, and each arithmetic is usable on its own, through the DAGs,
or as the coefficient type of another. The hinge is a single traits class,
`mc::Op<T>`: any type specializing it can be propagated through a DAG.
Composing arithmetics is supported in C++ as a template choice. Examples include
Taylor and Chebyshev models propagating their remainder in McCormick arithmetic
rather than in intervals [@sahlodin2011; @bompadre2013], and superposition
relaxations parameterizing their univariate summands as piecewise-constant or
continuous piecewise-linear functions [@zha2026]. Each header (`interval.hpp`,
`mccormick.hpp`, `scmodel.hpp`, `polimage.hpp`, `fdiff.hpp`, ...) can equally
be included alone and used without any DAG, and the interval type beneath them
all is interchangeable --- MC++'s own, Boost, PROFIL/BIAS or FILIB++. In
contrast, pyMC++, built with pybind11 [@pybind11], binds a fixed selection of
these instantiations --- the relaxations and polynomial models over the default
interval backend, and both superposition variants --- while the full template
space is open to C++ users. When supported by the selected interval backend,
directed rounding provides verified interval enclosures. However, the rest of
the enclosures remain subject to rounding errors.

The DAG is itself extensible, through *external operations*: user-defined nodes
carrying their own evaluation and differentiation rules in each arithmetic.
One such operation encapsulates external code that the DAG treats as opaque,
such as a Python function or a proprietary numerical routine, provided that
suitable evaluation rules are supplied for the required arithmetics and,
optionally, a differentiation rule. The same mechanism supports affine and
polynomial subexpressions, multilayer perceptrons imported from PyTorch, and
nested DAGs, and it is how the companion library CRONOS
(<https://github.com/omega-icl/cronos>, `pip install cronos-mcpp`) embeds
solvers for systems of algebraic, differential and partial differential
equations as single DAG nodes. CANON (<https://github.com/omega-icl/canon>,
local and global optimization) and MAGNUS
(<https://github.com/omega-icl/magnus>, modeling and analysis of uncertain
systems) build on the same graphs.

MC++ is released under the Eclipse Public License 2.0, with the GNU GPL v2 or
later as a secondary license. Documentation, Python scripts and Jupyter
notebooks covering each arithmetic accompany the source, among them a companion
notebook that reproduces every listing in this paper and runs in the browser
without installation.[^colab]

[^colab]: <https://colab.research.google.com/github/omega-icl/mcpp/blob/version-5/notebook/joss_paper.ipynb>


# Research impact statement

MC++ is the implementation of record for a line of set arithmetics first
published in peer-reviewed journals: McCormick-based relaxations of algorithms
[@mitsos2009], multivariate McCormick relaxations [@tsoukalas2014],
McCormick--Taylor models [@sahlodin2011; @bompadre2013], Chebyshev model
arithmetic [@rajyaguru2017], ellipsoidal arithmetic [@villanueva2015] and, most
recently, relaxation via separable estimators [@zha2026]. Each was developed and
validated in MC++ --- in its predecessor libMC in the earliest case --- so the
library is where these methods are available to others rather than only
described.

Its use extends well beyond its authors. MAiNGO [@maingo] is built on MC++ and
vendors it as a dependency; Pyomo [@pyomo] wraps it to bound factorable
functions with McCormick relaxations; and the companion libraries CRONOS, CANON
and MAGNUS build on its DAGs. Applications include surrogate-based global
optimization, with trained neural networks [@schweidtmann2019] and Gaussian
processes [@schweidtmann2021] embedded through MeLOn
(<https://git.rwth-aachen.de/avt-svt/public/MeLOn>), which builds on MAiNGO,
and model-based design of experiments, whose computations run entirely on MC++
DAGs [@sandrin2025]. Despite more than a decade of such use, MC++ has never had
a software paper, and users have cited method papers in its place --- the gap
this paper closes.

# AI usage disclosure

Generative AI was used in both the software and this manuscript. Anthropic's
Claude models --- Claude Opus 4.8, Claude Opus 5, Claude Opus 5.5 and Claude
Fable 5.1 --- assisted in debugging and improving parts of the C++ code,
including the Python binders, and in drafting sections of this paper; OpenAI's
ChatGPT was used to review drafts of the manuscript.

The authors reviewed, edited and validated all AI-assisted output, and made the
core design decisions: the architecture of the library, the arithmetics it
implements and the interfaces it exposes are the authors' own. AI-assisted code
was accepted only after passing the repository's `check` target, which gates
version consistency, the recorded ABI fingerprint, the Python check scripts in
`src/pymcpp/` and a warning-free documentation build, with forward and reverse
derivatives cross-checked against each other. Every code listing in this paper
was executed against the released `pymcpp` wheel, and the companion notebook
asserts each quoted figure.

# Acknowledgements

This work was supported by UKRI Engineering and Physical Sciences Research
Council grant EP/W003317/1 (ADOPT). We thank the contributors listed in the
`AUTHORS` file, in particular Dominik Bongartz, Boris Houska, Alexander Mitsos,
Jaromil Najman, Nikola Peric, Jai Rajyaguru, Ali M. Sahlodin, Mario E.
Villanueva and Yanlin Zha, for the arithmetic implementations. We are grateful
to Terrence Crombie, Tanuj Karia, Radoslav Paulen, Yukun Wang, Hongzhen Zhang
and Lifeng Zhang for exercising the library and reporting the defects they
found.

# References
