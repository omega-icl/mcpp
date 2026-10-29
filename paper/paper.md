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

A *factorable* function is one obtained by finite recursive composition of
binary sums and products with a library of univariate intrinsics. Nearly every
algebraic model written in science and engineering is factorable, and it is the
class on which rigorous bounding theory is built: interval arithmetic, convex
and concave relaxations, polynomial models with remainder bounds, and polyhedral
relaxations all propagate operation by operation through a factored expression
[@chachuat2015].

MC++ is a C++ toolkit that represents factorable functions as directed acyclic
graphs (DAGs) and propagates through them not only numbers and derivatives but
*sets*. The same expression graph can be evaluated in real arithmetic,
differentiated symbolically in forward or reverse mode, decomposed into sparse
polynomial and transcendental subexpressions, enclosed in intervals, relaxed by
McCormick convex/concave estimators with subgradients [@mccormick1976;
@mitsos2009; @tsoukalas2014], approximated by Taylor or Chebyshev models with
rigorous remainders [@bompadre2013; @rajyaguru2017], bounded by ellipsoids
[@villanueva2015] or Hessian spectra, and relaxed into polyhedra or by separable
estimators [@zha2026]. The companion library pyMC++ (module `pymcpp`) exposes
all of this to Python through pybind11 [@pybind11], as prebuilt wheels for
CPython 3.10--3.14 on Linux, macOS and Windows: the whole toolkit is `pip
install pymcpp` away.

# Statement of need

Deterministic global optimization, reachability analysis of dynamic systems and
rigorous model validation all rest on the same primitive: given a factorable
function and a host set for its variables, produce a guaranteed enclosure of its
range, or a convex relaxation a solver can minimize. These are neither point
evaluations nor derivatives. They are set-valued computations, propagated
operation by operation through the expression, each step carrying a correctness
guarantee.

Researchers needing them have had to choose between implementing an arithmetic
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
efficiently; none propagates set-valued information. They answer what $f(x)$ and
$\nabla f(x)$ are *here*, never what values $f$ can take over a set, or what
convex function underestimates it there.

On the other side are rigorous-bounding libraries, each specialized to one
arithmetic and most of them ageing. Intervals are served by filib++ [@filib],
PROFIL/BIAS [@profilbias], Boost.Interval [@boostinterval] and IBEX [@ibex], the
last pairing them with contractor programming. Taylor models are served by COSY
Infinity, which is not openly licensed, by TaylorModels.jl [@taylormodels_jl;
@taylorseries_jl], and inside reachability tools such as Flow\* [@flowstar] and
CORA [@cora]; in those the set arithmetic serves reachable sets of dynamic
systems rather than arbitrary factorable expressions. Chebyshev-model arithmetic
has essentially no open alternative.

Deterministic global optimization solvers --- BARON [@baron], ANTIGONE
[@antigone], SCIP [@scip], Couenne [@couenne], Alpine [@alpine], Gurobi
[@gurobi], EAGO [@eago] and MAiNGO [@maingo] --- do relax factorable functions,
but inside a branch-and-bound loop, and several are closed source.

The nearest reusable multi-arithmetic alternative is the Julia stack around
EAGO.jl, McCormick.jl and IntervalArithmetic.jl --- separate packages rather
than one DAG shared by every arithmetic, with no Chebyshev, ellipsoidal,
spectral or polyhedral relaxation, and no Python interface. For Python users the
only prior route has been Pyomo's `contrib.mcpp` module, which wraps an earlier
MC++ release through a shim the user compiles by hand.

Contributing these arithmetics to an existing package was not a realistic
alternative to building a library. The differentiation tools have no set-valued
type to extend: propagating an enclosure changes what a node computes rather
than adding an intrinsic. The solvers do have the machinery, but specialized to
their own search and not addressable from outside.

# Capability

Functions are defined on a DAG environment, an `FFGraph` object, with
participating variables of type `FFVar` connected by operations of type `FFOp`;
expressions are built by writing ordinary arithmetic on them. Identical
subexpressions are detected and shared as they are created. For instance, a DAG
for the two factorable functions $f(x,y) = e^{xy} - 2x^2y + 3$ and
$g(x,y) = \sqrt{x^2y + 4} + e^{xy}$ shares the subexpressions $xy$, $e^{xy}$ and
$x^2$, and comprises 16 operations rather than 21. A *subgraph* (`FFSubgraph`)
--- the `FFOp` list for a chosen set of dependents, in topological order, what
an automatic-differentiation tool would call a tape or a codelist --- is a
first-class object that can be held, inspected and re-used, so a function
evaluated many times is traversed once.

`eval` sweeps a subgraph once and returns the values of its dependents; `veval`
does the same for many points at once, over up to `DAG.options.MAXTHREAD`
threads; and `reval` sweeps in both directions, contracting the variable ranges
against prescribed bounds on the dependents by reverse interval propagation.

What `eval` computes is decided by the type of the values it is given, not by
the function called. Floats evaluate $f$; intervals enclose its range over the
host set; and the relaxation types return convex and concave estimators with
their subgradients, polynomial models with remainder bounds, separable
superposition estimators, or the linear cuts of a polyhedral relaxation.
Spectral bounds on Hessians and ellipsoidal enclosures of multivariate systems
propagate the same way.

Derivatives are generated symbolically, in forward (`fdiff`) or reverse
(`bdiff`) accumulation, as a sparse triplet of *new DAG nodes*. Being
expressions in the same graph, the entries can be differentiated again for
Hessians, and evaluated in any of the arithmetics above.

A DAG can also be rewritten. `SLift` decomposes the dependents recursively into
sparse polynomial or rational subexpressions and transcendental ones,
introducing an auxiliary variable for each lifted term --- the reformulation
step that precedes tailored relaxation or quadratic reformulation in global
optimization. `FFPoly` and `FFMon` expose the sparse polynomial arithmetic
directly, with conversion between monomial and Chebyshev bases and factorization
in a chosen variable.

# User interface

The listings below build the DAG of $f$ and $g$, evaluate and bound $f$,
differentiate $f$, and decompose both dependents. Each was executed against
`pymcpp` 5.1.0; the values in the comments are what it returned.

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
the following describe $f$ over the host set $x \in [-1,1]$, $y \in [0,2]$. A
`McCormick` object carries the convex and concave relaxations of $f$ on the box
evaluated at a given point --- here $(0.5, 1.0)$, the second constructor
argument.

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

Each brackets $f(0.5,1) = 4.1487$, at a different cost: the order-3 Chebyshev
model narrows the range estimate from a width of 11.25 to 9.25, and the
superposition relaxation brackets the point within $[2.97, 4.81]$ where
McCormick gives $[2.00, 8.58]$. The subgradient `cvsub` supplies the cut a
bundle method or an LP relaxation needs, and the `PolCut` objects are ready for
an LP or MILP model.

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
`mc::Op<T>`, declaring what an arithmetic must provide: an interval hull, a
midpoint, a diameter, and the intrinsics. Any type specializing `mc::Op` can be
propagated through a DAG. Composing arithmetics is supported in C++ as a
template choice. Examples include McCormick relaxations and Taylor or Chebyshev
models templated on a chosen interval backend; Taylor and Chebyshev models
propagating their remainder in McCormick arithmetic rather than in intervals
[@sahlodin2011; @bompadre2013]; superposition relaxations parameterizing their
univariate summands as piecewise-constant or continuous piecewise-linear
functions [@zha2026]. Each header (`interval.hpp`, `mccormick.hpp`,
`scmodel.hpp`, `polimage.hpp`, ...) can equally be included alone and used
without any DAG. In contrast, pyMC++ binds a fixed selection of these
instantiations --- McCormick relaxations and the polynomial models over the
default interval backend, sparse Chebyshev models with interval coefficients,
and both superposition variants --- while the full template space is open to C++
users.

The DAG is itself extensible, through *external operations*: user-defined nodes
carrying their own evaluation and differentiation rules in each arithmetic. The
mechanism supports affine and polynomial subexpressions, multilayer perceptrons
imported from PyTorch, and nested DAGs, and it is how the companion library
CRONOS (<https://github.com/omega-icl/cronos>, `pip install cronos-mcpp`) embeds
solvers for systems of algebraic, differential and partial differential
equations as single DAG nodes. CANON (<https://github.com/omega-icl/canon>), for
local and global search in nonlinear optimization, and MAGNUS
(<https://github.com/omega-icl/magnus>), for modeling, analysis and design of
uncertain systems, build on the same graphs. At build time the interval backend
is interchangeable --- Boost by default, PROFIL/BIAS or FILIB++ on request ---
and the PyTorch interface is optional.

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
and MAGNUS build on its DAGs. Despite more than a decade of such use, MC++ has
never been described in a software paper, and users have cited method papers in
its place --- a gap this paper closes. The PyPI wheels remove the last barrier
for Python users.

# AI usage disclosure

Generative AI was used in both the software and this manuscript. Anthropic's
Claude models --- Claude Opus 4.8, Claude Opus 5, Claude Opus 5.5 and Claude
Fable 5.1 --- assisted in debugging and improving parts of the C++ code,
including the Python binders, and in drafting sections of this paper.

The authors reviewed, edited and validated all AI-assisted output, and made the
core design decisions: the architecture of the library, the arithmetics it
implements and the interfaces it exposes are the authors' own. AI-assisted code
was accepted only after passing the repository's `check` target, which gates
version consistency across source, documentation and build metadata, the
recorded ABI fingerprint, the Python check scripts in `src/pymcpp/`, and a
warning-free documentation build; derivative routines are additionally
cross-checked forward against reverse accumulation. Every code listing in this
paper was executed against the released `pymcpp` wheel, and the companion
notebook asserts each quoted figure.

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
