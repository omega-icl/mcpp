#!/usr/bin/env python3
"""check_ffextern.py -- external operations on an object (FFMLP on an MLP, FFDAGEXT on a DAGEXT) from
Python: the home of the object in a DAG, and options changed after the insertion.

Run against the freshly built module (the `check` CMake target does), it exits non-zero when any check
fails.  Every check names what it tests and prints what it got, so a failure is a diagnosis, not a flag.

What is checked: the object passed to an operation is copied once per DAG, and that copy -- `FFMLP.mlp`,
`FFDAGEXT.dagext` -- is what every operation on it evaluates; `options` of the operation object are the
options of the operations it inserted, so that a change made after the insertion changes the next
relaxation; the options of the copy change the next evaluation or differentiation; an object changed and
passed again, or another object of the same shape, is another operation; and what Python holds of a DAG
(options, the copy of an object) stays valid after the DAG and the operation object are gone.
"""
import gc
import math
import sys

import pymcpp as mc

n_pass = n_fail = 0


def check(what, got, want, tol=0.):
    global n_pass, n_fail
    ok = (abs(got - want) <= tol) if isinstance(want, (int, float)) and not isinstance(want, bool) else (got == want)
    n_pass += ok; n_fail += (not ok)
    print("  %-66s got=%-18s want=%-18s %s" % (what, repr(got) if not isinstance(got, float) else "%.12g" % got,
                                                 repr(want) if not isinstance(want, float) else "%.12g" % want,
                                                 "PASS" if ok else "FAIL"))


def section(title):
    print("\n--- " + title)


def relax(G, F, X):
    """Polyhedral relaxation of F over a box: (number of cuts, number of auxiliary variables)."""
    img = mc.PolImg()
    box = [mc.Interval(-1., 1.), mc.Interval(-0.5, 2.)]
    PX = [mc.PolVar(img, x, b) for x, b in zip(X, box)]
    img.generate_cuts(G.eval(F, X, PX))
    return (len(img.cuts), len(img.aux))


def operations(G, F):
    """Names of the external operations that F depends on."""
    return " ".join(sorted({str(op).split("[")[0].strip() for op in G.subgraph(F).l_op if "[" in str(op)}))


def jacobian(G, F, X, x):
    """Dense Jacobian of F at x by forward differentiation, by rows."""
    rows, cols, dF = G.fdiff(F, X)
    J = [[0.] * len(X) for _ in F]
    for r, c, v in zip(rows, cols, G.eval(dF, X, x)):
        J[r][c] = v
    return J


def network(scale=1.):
    NN = mc.MLP()
    NN.set_data([([[0.1, 1., -1.], [0.2, 0.5, 1.], [-0.3, -1., 0.7]], mc.MLP.TANH),
                 ([[0.1 * scale, 1., -2., 0.5], [0.3, -1., scale, 2.]], mc.MLP.LINEAR)])
    return NN


def value(NN, x):
    """The network evaluated on its own, through a DAG that holds nothing else."""
    G = mc.FFGraph()
    X = G.add_vars(2)
    return G.eval(mc.FFMLP()(X, NN), X, x)


x0 = [0.3, -0.7]

# ------------------------------------------------------------------ A. FFMLP: the copy of the network in the DAG
section("A. FFMLP: one copy of the network per DAG")
NN = network()
G = mc.FFGraph()
X, Z = G.add_vars(2), G.add_vars(2)
Op = mc.FFMLP()
check("before any call: attached", Op.attached, False)
check("before any call: mlp is None", Op.mlp is None, True)
F = Op(X, NN)
net = Op.mlp
check("after a call: attached", Op.attached, True)
check("mlp is a copy, not the network passed", net is not None and net is not NN, True)
check("mlp is the same object at each access", Op.mlp is net, True)
W = Op(Z, NN)
check("the same copy for the network on other variables", Op.mlp is net, True)
Op2 = mc.FFMLP()
F2 = Op2(X, NN)
check("another FFMLP object: the same operation", F2[0].id == F[0].id, True)
check("another FFMLP object: the same copy", Op2.mlp is net, True)
check("values in the DAG are those of the network", G.eval(F, X, x0) == value(NN, x0), True)

section("A'. FFMLP: a network changed and passed again is another operation")
y_old = G.eval(F, X, x0)
NN.set_data([([[0.5, -1., 1.], [0.2, 0.5, 1.], [-0.3, -1., 0.7]], mc.MLP.TANH),
             ([[0.1, 1., -2., 0.5], [0.3, -1., 1., 2.]], mc.MLP.LINEAR)])
Fnew = Op(X, NN)
check("another DAG variable", Fnew[0].id != F[0].id, True)
check("another copy", Op.mlp is not net, True)
check("values of the changed network", G.eval(Fnew, X, x0) == value(NN, x0), True)
check("values of the first network, unchanged", G.eval(F, X, x0) == y_old, True)

# ------------------------------------------------------------------ B. FFMLP: relaxation options after the insertion
section("B. FFMLP: relaxation options changed after the insertion")
NN = network()
G = mc.FFGraph()
X, Z = G.add_vars(2), G.add_vars(2)
Op = mc.FFMLP()
Op.options.RELAX = {Op.Options.INT}
F, W = Op(X, NN), Op(Z, NN)
check("options is the same object at each access", Op.options is Op.options, True)
r_int = relax(G, F, X)
check("RELAX = {INT}, set before: no auxiliary variable", r_int[1], 0)
Op.options.RELAX = {Op.Options.AUX}
r_aux = relax(G, F, X)
check("RELAX = {AUX}, set after: auxiliary variables", r_aux[1] > 0, True)
check("the same for the network on other variables", relax(G, W, Z) == r_aux, True)
GB = mc.FFGraph()
XB = GB.add_vars(2)
OpB = mc.FFMLP()
OpB.options.RELAX = {OpB.Options.AUX}
check("the same relaxation as with {AUX} set before", relax(GB, OpB(XB, NN), XB) == r_aux, True)
Op.options.RELAX = [Op.Options.INT]
check("and back to {INT}", relax(G, F, X) == r_int, True)
new = mc.FFMLP.Options()
new.RELAX = {mc.FFMLP.Options.AUX}
new.SCMODORD = 6
Op.options = new
check("options assigned as a whole: relaxation", relax(G, F, X) == r_aux, True)
check("options assigned as a whole: SCMODORD", Op.options.SCMODORD, 6)
check("the object assigned from is not the options", Op.options is not new, True)
check("another FFMLP object keeps its own options", OpB.options.SCMODORD != 6, True)

# ------------------------------------------------------------------ C. FFMLP: options of the network in the DAG
section("C. FFMLP: options of the network, changed in the DAG")
NN = network()
G = mc.FFGraph()
X, Z = G.add_vars(2), G.add_vars(2)
Op = mc.FFMLP()
F, W = Op(X, NN), Op(Z, NN)
_, _, dF = G.fdiff(F, X)
check("fdiff by default: the Jacobian operation", operations(G, dF), "GradMLP")
NN.options.TANGENT = True
_, _, dW = G.fdiff(W, Z)
check("option set on the network passed: no effect on the DAG", operations(G, dW), "GradMLP")
NN.options.TANGENT = False
Op.mlp.options.TANGENT = True
_, _, tF = G.fdiff(F, X)
_, _, tW = G.fdiff(W, Z)
check("option set on mlp: tangent operations at the next fdiff", operations(G, tF), "TanMLP")
check("for the network on other variables too", operations(G, tW), "TanMLP")
d, t = G.eval(dF, X, x0), G.eval(tF, X, x0)
check("same derivative values either way", max(abs(a - b) for a, b in zip(d, t)), 0., 1e-12)

# ------------------------------------------------------------------ D. FFMLP: a copy of the DAG, and lifetimes
section("D. FFMLP: each DAG has its own copy and its own options")
NN = network()
G = mc.FFGraph()
X = G.add_vars(2)
Op = mc.FFMLP()
Op.options.SCMODORD = 7
F = Op(X, NN)
opt1, net1 = Op.options, Op.mlp
G2 = mc.FFGraph()
X2 = G2.add_vars(2)
F2 = G2.insert(G, F)
opt1.SCMODORD = 9
F2again = Op(X2, NN)                      # the network is found in the copy; Op is now attached there
check("the network is found in the copy of the DAG", F2again[0].id == F2[0].id, True)
check("the copy has its own network", Op.mlp is not net1, True)
check("the copy has its own options", Op.options is not opt1, True)
check("with the values at the time of the copy", Op.options.SCMODORD, 7)
check("values in the copy", G2.eval(F2, X2, x0) == value(NN, x0), True)

section("D'. what Python holds stays valid after the DAG and the operation object are gone")
opt2, net2 = Op.options, Op.mlp
del Op, F, F2, F2again, X, X2, G, G2
gc.collect()
check("options of the first DAG", opt1.SCMODORD, 9)
check("options of the copy", opt2.SCMODORD, 7)
check("network of the first DAG", (net1.nin, net1.nout), (2, 2))
check("network of the copy", (net2.nin, net2.nout), (2, 2))

# ------------------------------------------------------------------ E. FFDAGEXT
section("E. FFDAGEXT: expression trees of the same shape are distinct operations")
G = mc.FFGraph()
X, Z = G.add_vars(2), G.add_vars(2)
G1, G2 = mc.FFGraph(), mc.FFGraph()
U, V = G1.add_vars(2), G2.add_vars(2)
D1, D2 = mc.DAGEXT(G1, U, [U[0] * U[1]]), mc.DAGEXT(G2, V, [V[0] + V[1]])
Op = mc.FFDAGEXT()
Y1, Y2 = Op(X, D1), Op(X, D2)
check("x0*x1 at (3,5)", G.eval(Y1, X, [3., 5.])[0], 15.)
check("x0+x1 at (3,5)", G.eval(Y2, X, [3., 5.])[0], 8.)
check("their derivatives", jacobian(G, Y1 + Y2, X, [3., 5.]), [[5., 3.], [1., 1.]])
D1.set(G1, U, [U[0] - U[1]])
Y3 = Op(X, D1)
check("the same object set to x0-x1: another operation", G.eval(Y1 + Y3, X, [3., 5.]), [15., -2.])

section("E'. FFDAGEXT: the copy in the DAG and the options, after the insertion")
D = mc.DAGEXT(G1, U, [mc.tanh(U[0] * U[1]) + mc.sqr(U[0])])
D.options.CPMAX = 3
Op = mc.FFDAGEXT()
check("before any call: dagext is None", Op.dagext is None, True)
Op.options.RELAX = {Op.Options.INT}
Y, W = Op(X, D), Op(Z, D)
tree = Op.dagext
check("dagext is a copy, the same for other variables", tree is not None and tree is not D and Op.dagext is tree, True)
check("the copy has the options of the object passed", tree.options.CPMAX, 3)
r_int = relax(G, Y, X)
check("RELAX = {INT}, set before: no auxiliary variable", r_int[1], 0)
Op.options.RELAX = {Op.Options.AUX}
r_aux = relax(G, Y, X)
check("RELAX = {AUX}, set after: auxiliary variables", r_aux[1] > 0, True)
check("the same for the tree on other variables", relax(G, W, Z) == r_aux, True)
Op.options.RELAX = {Op.Options.INT}
check("and back to {INT}", relax(G, Y, X) == r_int, True)
_, _, dY = G.fdiff(Y, X)
Op.dagext.options.AUTODIFF = mc.DAGEXT.Options.B
_, _, bY = G.bdiff(Y, X)
d, b = G.eval(dY, X, x0), G.eval(bY, X, x0)
ref = [x0[1] * (1. - math.tanh(x0[0] * x0[1]) ** 2) + 2. * x0[0], x0[0] * (1. - math.tanh(x0[0] * x0[1]) ** 2)]
check("derivatives, forward then backward in the tree", max(abs(u - v) for u, v in zip(d + b, ref + ref)), 0., 1e-12)
opt, tree = Op.options, Op.dagext
del Op, Y, W, dY, bY, X, Z, G
gc.collect()
check("options and tree still valid once the DAG is gone", (opt.SCMODORD, tree.nin, tree.nout), (3, 2, 1))

print("\n%d passed, %d failed" % (n_pass, n_fail))
sys.exit(1 if n_fail else 0)
