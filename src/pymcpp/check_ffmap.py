#!/usr/bin/env python3
"""check_ffmap.py -- networks as external operations of a DAG from Python: mc::FFMAP through its two
specializations, FFMLP on an MLP and FFKAN on a KAN.

Run against the freshly built module (the `check` CMake target does), it exits non-zero when any check
fails.  Every check names what it tests and prints what it got, so a failure is a diagnosis, not a flag.

The same checks are run on both operations, each against its network evaluated in plain Python: the
values in a DAG, alone and inside an expression; the network that the DAG holds -- a copy, returned by
`map` and by `MLP` or `KAN` -- and what becomes of it when the network is changed and passed again; the
interval bounds; the Jacobian, forward and backward, and the tangent and adjoint operations that the
options of the network ask for; the polyhedral relaxations selected by `options.RELAX`; and what Python
holds of a DAG after the DAG is gone.  Then what is specific to a KAN: its layers given as PWLU edges, the
PWLU methods that define an edge on the real line, the encodings of the edges, and the data that are
refused.  Last, `NETIO.read_dict`: a file read once, then given to a KAN or an MLP by `set_data`.
"""
import gc
import math
import os
import sys

import pymcpp as mc

n_pass = n_fail = 0


def check(what, got, want, tol=0.):
    global n_pass, n_fail
    ok = (abs(got - want) <= tol) if isinstance(want, (int, float)) and not isinstance(want, bool) else (got == want)
    n_pass += ok; n_fail += (not ok)
    print("  %-70s got=%-18s want=%-18s %s" % (what, repr(got) if not isinstance(got, float) else "%.12g" % got,
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


def jacobian(G, F, X, x, backward=False):
    """Dense Jacobian of F at x by forward or backward differentiation, by rows, and the operations in it."""
    rows, cols, dF = (G.bdiff if backward else G.fdiff)(F, X)
    J = [[0.] * len(X) for _ in F]
    for r, c, v in zip(rows, cols, G.eval(dF, X, x)):
        J[r][c] = v
    return J, operations(G, dF)


def dist(A, B):
    """Largest difference between two lists, or lists of lists, of numbers."""
    if A and isinstance(A[0], (list, tuple)):
        return max(dist(a, b) for a, b in zip(A, B))
    return max(abs(a - b) for a, b in zip(A, B))


# ------------------------------------------------------------------ the two networks, and their values in Python
MLP_W = ([[0.1, 1., -1.], [0.2, 0.5, 1.], [-0.3, -1., 0.7]], [[0.1, 1., -2., 0.5], [0.3, -1., 1., 2.]])


def mlp_network(scale=1.):
    NN = mc.MLP()
    NN.set_data([(MLP_W[0], mc.MLP.TANH),
                 ([[0.1 * scale, 1., -2., 0.5], [0.3, -1., scale, 2.]], mc.MLP.LINEAR)])
    return NN


def mlp_reference(x, scale=1.):
    """Values and Jacobian of mlp_network(scale) at x."""
    W2 = [[0.1 * scale, 1., -2., 0.5], [0.3, -1., scale, 2.]]
    h = [math.tanh(r[0] + r[1] * x[0] + r[2] * x[1]) for r in MLP_W[0]]
    y = [r[0] + sum(r[1 + k] * h[k] for k in range(3)) for r in W2]
    J = [[sum(r[1 + k] * (1. - h[k] ** 2) * MLP_W[0][k][1 + i] for k in range(3)) for i in range(2)] for r in W2]
    return y, J


# edges as (c, v, a, w, b): c + v u + sum_r a[r] max(w[r] u + b[r], 0); layer[j][i] from input i to output j
KAN_E = [[[(0.1, 0.5, [1., -0.8], [1., -1.5], [0.2, 0.3]), (-0.2, 0.3, [0.6], [2.], [0.5])],
          [(0., -0.4, [0.7, 0.4], [-1., 1.], [0.1, -0.6]), (0.3, 0.2, [-0.5], [1.], [-0.25])],
          [(0.2, 0.1, [0.9], [-2.], [0.4]), (0., 1., [-1.2, 0.3], [1., 1.], [0.15, -0.45])]],
         [[(0.05, 1., [0.5], [1.], [-0.1]), (0., -0.6, [0.8], [-1.], [0.2]), (0.1, 0.4, [-0.3], [1.5], [0.3])],
          [(-0.1, 0.2, [1.1], [-1.], [0.35]), (0.2, 0.7, [-0.4], [2.], [-0.5]), (0., -1., [0.6], [1.], [0.05])]]]


def kan_layers(scale=1.):
    return [[[mc.PWLU.relu(c * (scale if l else 1.), v, a, w, b) for (c, v, a, w, b) in row] for row in layer]
            for l, layer in enumerate(KAN_E)]


def kan_network(scale=1.):
    NN = mc.KAN()
    NN.set_data(kan_layers(scale))
    return NN


def kan_reference(x, scale=1.):
    """Values and Jacobian of kan_network(scale) at x (not at a breakpoint)."""
    def edge(e, u, s):
        c, v, a, w, b = e
        return (c * s + v * u + sum(ar * max(wr * u + br, 0.) for ar, wr, br in zip(a, w, b)),
                v + sum(ar * wr for ar, wr, br in zip(a, w, b) if wr * u + br > 0.))
    u, J = list(x), [[1. if i == k else 0. for k in range(len(x))] for i in range(len(x))]
    for l, layer in enumerate(KAN_E):
        s = scale if l else 1.
        ev = [[edge(e, u[i], s) for i, e in enumerate(row)] for row in layer]
        J = [[sum(ev[j][i][1] * J[i][k] for i in range(len(u))) for k in range(len(x))] for j in range(len(layer))]
        u = [sum(e[0] for e in row) for row in ev]
    return u, J


SPECS = [dict(name="FFMLP", tag="MLP", Op=mc.FFMLP, network=mlp_network, reference=mlp_reference),
         dict(name="FFKAN", tag="KAN", Op=mc.FFKAN, network=kan_network, reference=kan_reference)]
x0 = [0.3, -0.7]
points = [[0.3, -0.7], [-0.6, 0.45], [0.85, 0.1], [-0.15, -0.9]]

# ------------------------------------------------------------------ A-F. the same checks on both operations
for S in SPECS:
    name, tag = S["name"], S["tag"]

    section("A. %s: values in a DAG" % name)
    NN = S["network"]()
    G = mc.FFGraph()
    X, Z = G.add_vars(2), G.add_vars(2)
    Op = S["Op"]()
    check("before any call: attached", Op.attached, False)
    check("before any call: map is None", Op.map is None, True)
    F = Op(X, NN)
    check("number of results", len(F), 2)
    check("name of the operation", operations(G, F), tag)
    check("values: largest error over the points", max(dist(G.eval(F, X, x), S["reference"](x)[0]) for x in points),
          0., 1e-13)
    E = [F[0] * X[1] + mc.exp(F[1])]
    y0 = S["reference"](x0)[0]
    check("inside an expression", G.eval(E, X, x0)[0], y0[0] * x0[1] + math.exp(y0[1]), 1e-13)
    check("one result by its index", G.eval([Op(1, X, NN)], X, x0)[0], y0[1], 1e-13)

    section("B. %s: the network of the DAG" % name)
    net = Op.map
    check("after a call: attached", Op.attached, True)
    check("map is a copy, not the network passed", net is not None and net is not NN, True)
    check("map is the same object at each access", Op.map is net, True)
    check("%s is map" % tag, getattr(Op, tag) is net, True)
    W = Op(Z, NN)
    check("the same copy for the network on other variables", Op.map is net, True)
    Op2 = S["Op"]()
    F2 = Op2(X, NN)
    check("another %s object: the same operation" % name, F2[0].id == F[0].id, True)
    check("another %s object: the same copy" % name, Op2.map is net, True)
    NN2 = S["network"](scale=2.)
    Fnew = Op(X, NN2)
    check("another network: another DAG variable", Fnew[0].id != F[0].id, True)
    check("another network: another copy", Op.map is not net, True)
    check("values of the other network", dist(G.eval(Fnew, X, x0), S["reference"](x0, 2.)[0]), 0., 1e-13)
    check("values of the first network, unchanged", dist(G.eval(F, X, x0), y0), 0., 1e-13)

    section("C. %s: interval bounds" % name)
    NN = S["network"]()
    G = mc.FFGraph()
    X, Z = G.add_vars(2), G.add_vars(2)
    Op = S["Op"]()
    F, W = Op(X, NN), Op(Z, NN)
    B = G.eval(F, X, [mc.Interval(-1., 1.), mc.Interval(-1., 1.)])
    inside = all(B[k].l - 1e-12 <= S["reference"](x)[0][k] <= B[k].u + 1e-12 for x in points for k in range(2))
    check("bounds on [-1,1]^2 enclose the values at the points", inside, True)
    check("bounds are finite", all(math.isfinite(b.l) and math.isfinite(b.u) for b in B), True)

    section("D. %s: derivatives" % name)
    Jref = S["reference"](x0)[1]
    J, ops = jacobian(G, F, X, x0)
    check("forward: the Jacobian operation", ops, "Grad" + tag)
    check("forward: error", dist(J, Jref), 0., 1e-12)
    J, ops = jacobian(G, F, X, x0, backward=True)
    check("backward: the Jacobian operation", ops, "Grad" + tag)
    check("backward: error", dist(J, Jref), 0., 1e-12)
    NN.options.TANGENT = True
    check("TANGENT set on the network passed: no effect on the copy", jacobian(G, W, Z, x0)[1], "Grad" + tag)
    NN.options.TANGENT = False
    Op.map.options.TANGENT = True
    J, ops = jacobian(G, F, X, x0)
    check("TANGENT set on map: tangent operations at the next fdiff", ops, "Tan" + tag)
    check("tangent operations: error", dist(J, Jref), 0., 1e-12)
    Op.map.options.TANGENT = False
    Op.map.options.ADJOINT = True
    J, ops = jacobian(G, F, X, x0, backward=True)
    check("ADJOINT set on map: adjoint operations at the next bdiff", ops, "Adj" + tag)
    check("adjoint operations: error", dist(J, Jref), 0., 1e-12)
    Op.map.options.ADJOINT = False

    section("E. %s: relaxation options changed after the insertion" % name)
    check("options is the same object at each access", Op.options is Op.options, True)
    r_mc = relax(G, F, X)
    check("RELAX by default (MC): cuts, no auxiliary variable", r_mc[0] > 0 and r_mc[1] == 0, True)
    Op.options.RELAX = {Op.Options.INT}
    r_int = relax(G, F, X)
    check("RELAX = {INT}, set after: no auxiliary variable", r_int[1], 0)
    Op.options.RELAX = {Op.Options.AUX}
    r_aux = relax(G, F, X)
    check("RELAX = {AUX}, set after: auxiliary variables", r_aux[1] > 0, True)
    check("the same for the network on other variables", relax(G, W, Z) == r_aux, True)

    section("F. %s: what Python holds after the DAG is gone" % name)
    net, opt = Op.map, Op.options
    del G, X, Z, F, W, E, Fnew, F2, B, Op, Op2
    gc.collect()
    check("the network of the DAG: its sizes", (net.nin, net.nout), (2, 2))
    check("the options of the operation", opt.AUX in opt.RELAX and len(opt.RELAX) == 1, True)
    G = mc.FFGraph()
    X = G.add_vars(2)
    check("the network inserted in another DAG", dist(G.eval(S["Op"]()(X, net), X, x0), y0), 0., 1e-13)

# ------------------------------------------------------------------ G. KAN: edges, encodings, data refused
section("G. PWLU: an edge on the real line")
f = mc.PWLU([-1., 0., 1.], [1., 0., 1.], -1., 1.)
check("breakpoints: where the slope changes", f.bkpts, [0.])
check("value inside", f.val(0.5), 0.5, 1e-15)
check("value beyond the range, extended", f.val(7., True), 7., 1e-14)
check("slope beyond the range", f.slope(-7.), -1., 1e-15)
check("slope at a breakpoint: right derivative", f.slope(0.), 1., 1e-15)
g = mc.PWLU.relu(0.1, 0.5, [1., -0.8], [1., -1.5], [0.2, 0.3])
check("relu: breakpoints at -b/w", dist(g.bkpts, [-0.2, 0.2]), 0., 1e-15)
ref = lambda u: 0.1 + 0.5 * u + max(u + 0.2, 0.) - 0.8 * max(-1.5 * u + 0.3, 0.)
check("relu: values, extended", max(abs(g.val(u, True) - ref(u)) for u in (-5., -0.2, 0., 0.1, 0.2, 3.)), 0., 1e-14)

section("G'. KAN: layers and options")
NN = kan_network()
check("sizes", (NN.nin, NN.nout, NN.nlay), (2, 2, 2))
check("number of segments", NN.nwei > 0, True)
check("data: layers, outputs and inputs", [(len(L), len(L[0])) for L in NN.data], [(3, 2), (2, 3)])
check("data: an edge", dist(NN.data[0][0][0].bkpts, [-0.2, 0.2]), 0., 1e-15)
NB = mc.KAN()
for L in kan_layers():
    NB.append_data(L)
G = mc.FFGraph()
X = G.add_vars(2)
check("append_data layer by layer: the same network", G.eval(mc.FFKAN()(X, NB), X, x0) == G.eval(mc.FFKAN()(X, NN), X, x0),
      True)
try:
    NB.append_data([[mc.PWLU.relu(0., 1., [1.], [1.], [0.])] * 3])
    check("append_data of a layer with another number of inputs raises", False, True)
except ValueError:
    check("append_data of a layer with another number of inputs raises, the KAN unchanged", NB.nlay, 2)
try:
    NB.read_data("missing.pt")
    check("read_data of a missing file raises", False, True)
except RuntimeError:
    check("read_data of a missing file raises, the KAN unchanged", NB.nlay, 2)
check("ENCODING by default", NN.options.ENCODING, int(NN.Options.DEFAULT))
check("dag of the network: inputs and outputs", (len(NN.varin), len(NN.varout)), (2, 2))

section("G''. FFKAN: encodings of the edges")
Op = mc.FFKAN()
F = Op(X, NN)
Op.options.RELAX = {Op.Options.AUX}
res = {}
for enc in ("HULL", "CONT", "BIN", "SOS2"):
    Op.map.options.ENCODING = getattr(NN.Options, enc)
    res[enc] = relax(G, F, X)
    print("  %-5s %d cuts, %d auxiliary variables" % ((enc,) + res[enc]))
check("ENCODING is an option of the network of the DAG", Op.map.options.ENCODING, int(NN.Options.SOS2))
check("HULL: the fewest auxiliary variables", all(res["HULL"][1] < res[e][1] for e in ("CONT", "BIN", "SOS2")), True)
check("BIN: more auxiliary variables than CONT", res["BIN"][1] > res["CONT"][1], True)
Op.map.options.ENCODING = NN.Options.BIN
Op.map.options.SEPARATE = True
r_sep = relax(G, F, X)
check("SEPARATE: more auxiliary variables than on one subdivision", r_sep[1] > res["BIN"][1], True)

# ------------------------------------------------------------------ H. NETIO: a file read once
section("H. NETIO: a dictionary read once, then given to the network")
data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "test", "ML")
for fname, Net, Op, kind in (("dict_kan_scaled.pt", mc.KAN, mc.FFKAN, "KANmodel"), ("dict_ann_tanh.pt", mc.MLP, mc.FFMLP, "ANNmodel")):
    D = mc.NETIO.read_dict(os.path.join(data, fname))
    check("%s: type" % fname, D.type, kind)
    NF, ND = Net(), Net()
    NF.read_data(os.path.join(data, fname))
    ND.set_data(D)
    check("%s: sizes, from widths" % fname, (ND.nin, ND.nout), (D.widths[0], D.widths[-1]))
    check("%s: scaling of the inputs, one pair of bounds per input or none" % fname,
          len(D.xlower) in (0, ND.nin) and len(D.xupper) == len(D.xlower), True)
    G = mc.FFGraph()
    X = G.add_vars(ND.nin)
    x = [(D.xlower[i] + 0.3 * (D.xupper[i] - D.xlower[i])) if D.xlower else 0.3 for i in range(ND.nin)]
    check("%s: set_data(D) is the network of read_data" % fname, G.eval(Op()(X, ND), X, x) == G.eval(Op()(X, NF), X, x), True)
DK, DA = mc.NETIO.read_dict(os.path.join(data, "dict_kan_scaled.pt")), mc.NETIO.read_dict(os.path.join(data, "dict_ann_tanh.pt"))
NK = mc.KAN()
NK.set_data(DK)
try:
    NK.set_data(DA)
    check("KAN.set_data of a dense network raises", False, True)
except ValueError:
    check("KAN.set_data of a dense network raises, the KAN unchanged", NK.nin, DK.widths[0])
try:
    mc.MLP().set_data(DK)
    check("MLP.set_data of a Kolmogorov-Arnold network raises", False, True)
except ValueError:
    check("MLP.set_data of a Kolmogorov-Arnold network raises", True, True)
try:
    mc.NETIO.read_dict(os.path.join(data, "missing.pt"))
    check("read_dict of a missing file raises", False, True)
except RuntimeError:
    check("read_dict of a missing file raises", True, True)

print("\n%d passed, %d failed" % (n_pass, n_fail))
sys.exit(1 if n_fail else 0)
