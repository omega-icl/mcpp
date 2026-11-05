---
title: 'MC^++^ and PyMC^++^: construction, manipulation and bounding of factorable functions'
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
evaluating it, however many times: what values an output can take over a whole
set of inputs, whether a constraint holds throughout a prescribed domain, or
within what bounds the global optimum of a nonconvex problem provably lies. Such
questions call for enclosures and convex relaxations, obtained by propagating
sets rather than numbers through the model's expressions.

MC^++^ is a C++ toolkit that does this for factorable functions, the class
obtained by composing arithmetic operations with standard intrinsic functions,
which covers most algebraic models in science and engineering. A function is
represented once, as a directed acyclic graph (DAG) of its operations, and that
graph is then evaluated in whichever arithmetic a task needs: real numbers,
derivatives, interval enclosures, convex and concave relaxations, polynomial
models with remainder bounds, and others. PyMC^++^ (module `pymcpp`) brings the
toolkit to Python as prebuilt wheels on PyPI.

# Statement of need

Deterministic global optimization, reachability analysis of dynamic systems and
rigorous model validation all rest on the same primitive: given a factorable
function and a host set for its variables, produce a guaranteed enclosure of its
range, or a convex relaxation a solver can minimize [@chachuat2015]. Researchers
in these fields, and the developers of the solvers and modeling tools that serve
them, have so far had to choose between implementing an arithmetic from scratch
and using one locked inside a solver. What has been missing is a library that
offers several such arithmetics, composable and extensible, over one shared
expression graph, and is reachable from Python, where much of the modeling work
that needs rigorous bounds is done.

# State of the field

Expression-graph and automatic-differentiation tools --- CasADi [@casadi], CppAD
[@cppad], ADOL-C [@adolc], the modeling languages Pyomo [@pyomo; @pyomobook],
JuMP [@jump], AMPL [@ampl] and GAMS [@gams], and the tracing autodiff of PyTorch
[@pytorch] and JAX [@jax] --- answer what $f(x)$ and $\nabla f(x)$ are *here*.
What values $f$ can take over a set, or what convex function underestimates it
there, is at best handled by an extension.

Rigorous-bounding libraries are each specialized to one arithmetic: intervals in
filib++ [@filib], PROFIL/BIAS [@profilbias], Boost.Interval [@boostinterval] and
IBEX [@ibex]; Taylor models in COSY Infinity [@cosy], which is not openly
licensed, TaylorModels.jl [@taylormodels_jl; @taylorseries_jl], and
reachability tools such as Flow\* [@flowstar] and CORA [@cora]. Chebyshev model
arithmetic has, to our knowledge, no open alternative. Global optimization
solvers --- BARON [@baron], ANTIGONE [@antigone], SCIP [@scip], Couenne
[@couenne], Alpine [@alpine], Gurobi [@gurobi], EAGO [@eago] and MAiNGO
[@maingo] --- relax factorable functions, but inside their own search and not
addressable from outside. The nearest reusable alternative, the Julia stack
around EAGO.jl, McCormick.jl [@eago] and IntervalArithmetic.jl
[@intervalarithmetic_jl], offers fewer arithmetics as separate packages and no
Python interface. For Python users, the only prior route has been Pyomo's
separately compiled `contrib.mcpp` extension to an earlier MC^++^ release.

Contributing these arithmetics to an existing package was not realistic.
Several graph-building tools admit user-defined value types, but none was
designed to host a collection of interchangeable and composable arithmetics,
with one serving as the bounding type of another --- which is why Pyomo
delegates its bounding to MC^++^.

# Software design

Expressions are built by writing ordinary arithmetic on the variables of a DAG,
and identical subexpressions are shared as they are created. For instance, the
DAG of $f(x,y) = e^{xy} - 2x^2y + 3$ and $g(x,y) = \sqrt{x^2y + 4} + e^{xy}$
shares $xy$, $e^{xy}$ and $x^2$. A *subgraph*, the operation list for chosen
dependents in topological order, is a first-class object, so repeated evaluation
reuses one precomputed sequence. A single routine sweeps it, with a vectorized,
multithreaded form for many points and a reverse form that contracts variable
ranges against bounds on the dependents.

What that routine computes is decided by the type of the values it is given.
Floats evaluate the function; intervals enclose its range over the host set
[@moore2009]; and the relaxation types return convex and concave estimators
[@mccormick1976] with their subgradients [@mitsos2009], polynomial models with
remainder bounds [@neumaier2003; @rajyaguru2017], or separable superposition
estimators [@zha2026]. Spectral bounds on Hessians [@monnigmann2011], the linear
cuts of a polyhedral relaxation [@smith1999; @baron], and ellipsoidal enclosures
[@villanueva2015] propagate the same way. Derivatives are generated
symbolically, in forward or reverse accumulation [@griewank2008], as sparse DAG
nodes that can be differentiated again or evaluated in any of these arithmetics.
A DAG can also be decomposed into polynomial and transcendental subexpressions,
the reformulation step that precedes relaxation in global optimization.

MC^++^ is header-only. The hinge is a single traits class, `mc::Op<T>`: any type
specializing it can be propagated through a DAG, used on its own, or serve as
the bounding type of another --- for instance, Taylor and Chebyshev models
propagating their remainder in McCormick arithmetic [@sahlodin2011;
@bompadre2013], or superposition relaxations parameterizing their univariate
summands as piecewise-constant or piecewise-linear functions [@zha2026]. The
interval type beneath them all is interchangeable --- MC^++^'s own, Boost,
PROFIL/BIAS or filib++ --- and directed rounding, where the backend supports it,
makes the interval enclosures verified; the other enclosures remain subject to
rounding errors. PyMC^++^, built with pybind11 [@pybind11], binds a fixed
selection of these instantiations, while the full template space is open to C++
users.

The DAG is extensible through *external operations*, user-defined nodes
carrying their own evaluation and differentiation rules. One such operation
encapsulates opaque external code, such as a Python function or a proprietary
routine; others support affine and polynomial subexpressions, multilayer
perceptrons imported from PyTorch, and nested DAGs. The companion library CRONOS
(<https://github.com/omega-icl/cronos>) uses them to embed solvers for
algebraic, differential and partial differential equations as single DAG nodes,
while CANON (<https://github.com/omega-icl/canon>, optimization) and MAGNUS
(<https://github.com/omega-icl/magnus>, uncertainty analysis) build on the same
graphs.

MC^++^ is released under the Eclipse Public License 2.0, with the GNU GPL v2 or
later as a secondary license.

# Example

The Python listing below builds the DAG of $f(x,y) = e^{xy} - 2x^2y + 3$; then,
the same `eval` call returns a value or an enclosure, depending only on the type
of its inputs. \autoref{fig:relax3d} plots the McCormick, Chebyshev and
superposition enclosures over the host set $[-1,1]\times[0,2]$.

![Enclosures of $f$ on $[-1,1]\times[0,2]$: (a) McCormick relaxations, (b) third-order sparse Chebyshev model, (c) piecewise-linear superposition estimators with eight segments per variable. The exact function is the black wireframe; under- and overestimators are in blue and orange.\label{fig:relax3d}](figures/relaxations_3d.png)

```python
import pymcpp as mc

DAG = mc.FFGraph()
X, Y = DAG.add_var("X"), DAG.add_var("Y")
F = mc.exp(X*Y) - 2*X**2*Y + 3
SGF = DAG.subgraph([F])

DAG.eval(SGF, [F], [X, Y], [0.5, 1.0])     # [4.1487]

IX, IY = mc.Interval(-1., 1.), mc.Interval(0., 2.)
DAG.eval(SGF, [F], [X, Y], [IX, IY])       # [-0.8647, 10.3891]

xMC, yMC = mc.McCormick(IX, 0.5), mc.McCormick(IY, 1.0)
[FMC] = DAG.eval(SGF, [F], [X, Y], [xMC, yMC])
FMC.cv, FMC.cc                             # 2.0, 8.5756 at (0.5, 1)

CM = mc.SCModel(3)                         # Chebyshev model, order 3
[FCM] = DAG.eval(SGF, [F], [X, Y], [mc.SCVar(CM, 0, IX), mc.SCVar(CM, 1, IY)])
FCM.B                                      # [-1.5230, 7.7258]

SUP = mc.PWLSModel(2)                      # superposition, 8 segments/variable
[FSU] = DAG.eval(SGF, [F], [X, Y],
                 [mc.PWLSVar(SUP, 0, IX, 8), mc.PWLSVar(SUP, 1, IY, 8)])
FSU.l(), FSU.u()                           # -0.8647, 9.3891
```

A companion notebook reproduces these results, with differentiation,
decomposition and other arithmetics; it runs in the browser without
installation.[^colab]

[^colab]: <https://colab.research.google.com/github/omega-icl/mcpp/blob/version-5/notebook/joss_paper.ipynb>

# Research impact statement

MC^++^ is the implementation of record for a line of set arithmetics published
in peer-reviewed journals: McCormick-based relaxations of algorithms
[@mitsos2009], multivariate McCormick relaxations [@tsoukalas2014],
McCormick--Taylor models [@sahlodin2011; @bompadre2013], Chebyshev model
arithmetic [@rajyaguru2017], ellipsoidal arithmetic [@villanueva2015] and
relaxation via separable estimators [@zha2026], each developed and validated in
MC^++^ or its predecessor libMC.

Its use extends beyond its authors. MAiNGO [@maingo] is built on MC^++^ and
vendors it as a dependency, and Pyomo [@pyomo] wraps it for McCormick
relaxations. Applications include surrogate-based global optimization with
embedded neural networks [@schweidtmann2019] and Gaussian processes
[@schweidtmann2021], through MeLOn
(<https://git.rwth-aachen.de/avt-svt/public/MeLOn>), and model-based design of
experiments running entirely on MC^++^ DAGs [@sandrin2025]. Despite more than a
decade of such use, MC^++^ has never had a software paper, and users have cited
method papers in its place.

# AI usage disclosure

Anthropic's Claude models (Claude Opus 4.8, Opus 5, Opus 5.5 and Fable 5.1)
assisted in debugging and improving parts of the C++ code and Python binders,
and in drafting sections of this paper; OpenAI's ChatGPT was used to review
drafts. The authors reviewed, edited and validated all AI-assisted output and
made the core design decisions. AI-assisted code was accepted only after passing
the repository's `check` target, which gates version consistency, the ABI
fingerprint, the Python check scripts and a warning-free documentation build.
Every listing in this paper was executed against the released `pymcpp` wheel,
and the companion notebook reproduces each quoted value.

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
