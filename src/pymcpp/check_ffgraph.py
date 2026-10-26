#!/usr/bin/env python3
"""check_ffgraph.py -- FFGraph and FFSubgraph from Python: construction, evaluation, differentiation, and
the subgraph as a first-class object.

Run against the freshly built module (the `check` CMake target does), it exits non-zero when any check
fails.  Every check names what it tests and prints what it got, so a failure is a diagnosis, not a flag.

The checks on FFSubgraph are the gate of its binding (2026-10-26): a subgraph formed from Python can be
held, inspected and handed back to eval(), and an evaluation through a held subgraph gives the same
numbers as one that forms its own.  The check on the slot guard is the one the binding makes possible
from Python: an evaluation nested inside another's subgraph -- a subgraph held, then evaluated twice over
without clearing -- must leave the first subgraph's numbers intact.
"""
import math
import sys

import pymcpp as mc

n_pass = n_fail = 0


def check(what, got, want, tol=0.):
    global n_pass, n_fail
    ok = (abs(got - want) <= tol) if isinstance(want, (int, float)) and not isinstance(want, bool) else (got == want)
    n_pass += ok; n_fail += (not ok)
    print("  %-56s got=%-22s want=%-22s %s" % (what, repr(got) if not isinstance(got, float) else "%.12g" % got,
                                                 repr(want) if not isinstance(want, float) else "%.12g" % want,
                                                 "PASS" if ok else "FAIL"))


def section(title):
    print("\n--- " + title)


# ------------------------------------------------------------------ A. construction and evaluation
section("A. a DAG of three variables, evaluated at a point")
G = mc.FFGraph()
x, y, z = G.add_var("x"), G.add_var("y"), G.add_var("z")
f1 = x * y + mc.sqr(z)                   # 1*2 + 9 = 11
f2 = mc.exp(x) - y * z                   # e - 6
f3 = x + y + z                           # 6
vals = G.eval([f1, f2, f3], [x, y, z], [1., 2., 3.])
check("eval f1 = x*y + z^2 at (1,2,3)", vals[0], 11., 1e-12)
check("eval f2 = exp(x) - y*z", vals[1], math.exp(1.) - 6., 1e-12)
check("eval f3 = x + y + z", vals[2], 6., 1e-12)

# ------------------------------------------------------------------ B. the subgraph as an object
section("B. FFSubgraph: formed, held, inspected, and passed back to eval")
sg = G.subgraph([f1, f2, f3])
check("subgraph is an FFSubgraph", isinstance(sg, mc.FFSubgraph), True)
check("len(subgraph) == number of operations in l_op", len(sg), len(sg.l_op))
check("len(subgraph) > 0", len(sg) > 0, True)
check("subgraph[0] is an FFOp", isinstance(sg[0], mc.FFOp), True)
check("subgraph[-1] is the last of l_op", sg[-1].name(), sg.l_op[-1].name())
check("len_tap >= number of operations", sg.len_tap >= len(sg), True)
n_ops_first = len(sg)
vals_sg = G.eval(sg, [f1, f2, f3], [x, y, z], [1., 2., 3.])
check("eval through the held subgraph == eval without one (f1)", vals_sg[0], vals[0], 0.)
check("eval through the held subgraph == eval without one (f2)", vals_sg[1], vals[1], 0.)
check("eval through the held subgraph == eval without one (f3)", vals_sg[2], vals[2], 0.)
vals_sg2 = G.eval(sg, [f1, f2, f3], [x, y, z], [0.5, -1., 2.])
check("the held subgraph evaluates again at a new point (f3)", vals_sg2[2], 1.5, 1e-12)
check("the subgraph is unchanged by evaluation", len(sg), n_ops_first)
sub1 = G.subgraph([f3])
check("a dependent's own subgraph is smaller than the three's", len(sub1) < len(sg), True)
check("x + y + z: 2 additions and 3 variables, 5 operations", len(sub1), 5)
try:
    _ = sg[len(sg)]
    check("subgraph[len] raises IndexError", False, True)
except IndexError:
    check("subgraph[len] raises IndexError", True, True)
sg_copy = mc.FFSubgraph(sg)
check("the copy has the same operation count", len(sg_copy), len(sg))
sg_copy.clear()
check("clear() empties the copy", len(sg_copy), 0)
check("... and leaves the original alone", len(sg), n_ops_first)

# ------------------------------------------------------------------ C. the slot guard, from Python
section("C. the value slots survive a nested evaluation (FFSubgraph::save_val)")
# A subgraph of f1 held; the same graph evaluated on a DIFFERENT subgraph that shares operations with it
# (f1 and f3 share x, y, z); then the held subgraph evaluated again.  Before the guard, the second
# evaluation nulled the slots the first had bound, and the third read through null.
sg1 = G.subgraph([f1])
a = G.eval(sg1, [f1], [x, y, z], [1., 2., 3.])[0]
_ = G.eval([f3, f1], [x, y, z], [4., 5., 6.])
b = G.eval(sg1, [f1], [x, y, z], [1., 2., 3.])[0]
check("f1 through its held subgraph, before and after another evaluation", b, a, 0.)
check("... and it is the right number", b, 11., 1e-12)

# ------------------------------------------------------------------ D. differentiation
section("D. forward and backward derivatives agree with each other and with the hand values")
dfdx = G.fdiff([f1, f2], [x, y, z])              # (rows, cols, entries) sparse triplets
dfdx_b = G.bdiff([f1, f2], [x, y, z])
J_f = {}
for i, j, e in zip(*dfdx): J_f[(i, j)] = e
J_b = {}
for i, j, e in zip(*dfdx_b): J_b[(i, j)] = e
check("fdiff and bdiff have the same sparsity", sorted(J_f), sorted(J_b))
want = {(0, 0): 2., (0, 1): 1., (0, 2): 6., (1, 0): math.exp(1.), (1, 1): -3., (1, 2): -2.}
check("fdiff finds the 6 nonzeros of a 2x3 Jacobian", len(J_f), 6)
ent = sorted(J_f)
got_f = G.eval([J_f[k] for k in ent], [x, y, z], [1., 2., 3.])
got_b = G.eval([J_b[k] for k in ent], [x, y, z], [1., 2., 3.])
for k, gf, gb in zip(ent, got_f, got_b):
    check("d f%d / d %s  forward" % (k[0] + 1, "xyz"[k[1]]), gf, want[k], 1e-12)
    check("d f%d / d %s  backward" % (k[0] + 1, "xyz"[k[1]]), gb, want[k], 1e-12)

# ------------------------------------------------------------------ E. a dependent read by another dependent
section("E. backward derivatives when a dependent is consumed by another, or is itself an independent")
# SBAD truncated a dependent's adjoint to its seed, skipping what its consumers hand back (fixed
# 2026-10-26): the row of t = 2s + x*y lost the 2s path when s was also a dependent, and the row of s
# lost column y when y was also a dependent.  Both shapes are common: an output defined from another.
s = mc.sqr(x) + mc.sqr(y)
t2 = 2. * s + x * y
for deps, label in (([s, t2], "{s, t}"), ([t2, s], "{t, s}"), ([s, y], "{s, y}")):
    Jf = {}; Jb = {}
    for i, j, e in zip(*G.fdiff(deps, [x, y])): Jf[(i, j)] = e
    for i, j, e in zip(*G.bdiff(deps, [x, y])): Jb[(i, j)] = e
    check("%s: bdiff has fdiff's sparsity" % label, sorted(Jb), sorted(Jf))
    ks = sorted(Jf)
    vf = G.eval([Jf[k] for k in ks], [x, y], [1., 2.])
    vb = G.eval([Jb.get(k, 0. * x) for k in ks], [x, y], [1., 2.])
    for k, a, b in zip(ks, vf, vb):
        check("%s: entry %s backward == forward" % (label, k), b, a, 1e-12)

# ------------------------------------------------------------------
print("\n==========  check_ffgraph: %d passed, %d failed -- %s  ==========\n"
      % (n_pass, n_fail, "FAILURES" if n_fail else "ALL PASS"))
sys.exit(1 if n_fail else 0)
