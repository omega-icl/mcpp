#!/usr/bin/env python3
"""Sparse operand/operation incidence matrix of an MC++ DAG, via FFSubgraph.

Builds a DAG for a given expression, extracts its FFSubgraph, and walks the
subgraph's operation list to assemble the incidence pattern as a list of
index pairs ``(i, j)`` meaning "operand i is an input of operation j".

    rows    i = DAG nodes (FFVar), numbered in order of definition
    columns j = operations, numbered in evaluation order (= FFSubgraph order)

Requires the FFSubgraph binder additions (l_op / tape / types / __iter__).

Usage
-----
    python3 dag_incidence.py
    python3 dag_incidence.py --vars x,y --expr "x*exp(x+sqr(y)) - sqr(y)" \
                            --expr "sqrt(x+1) + sin(y)"
    python3 dag_incidence.py --example signomial --with-outputs --no-map
"""

import argparse
import sys

import pymcpp as mc

# ---------------------------------------------------------------------------
# The expression. Edit DEFAULT_* or pass --vars / --expr.
# ---------------------------------------------------------------------------
DEFAULT_VARS = "x,y"
DEFAULT_EXPRS = ["x*exp(x + sqr(y)) - sqr(y)", "sqrt(x + 1) + sin(y)"]

EXAMPLES = {
    # sqr(y) is shared between the two dependents -> that row has fan-out 2
    "default": (DEFAULT_VARS, DEFAULT_EXPRS),
    "chain": ("x", ["exp(sin(sqr(x)))"]),
    "signomial": ("x,y,z", ["pow(x, 3) * pow(y, 2) / z", "x*y + y*z + z*x"]),
    "shared": ("x,y", ["sqr(x + y)", "sqr(x + y) + x", "log(sqr(x + y))"]),
}

# Functions made available to --expr, all from pymcpp.
_FUNCS = (
    "inv sqr sqrt exp log xlog cos sin tan acos asin atan cosh sinh tanh "
    "fabs relu fstep bstep erf erfc pow cheb max min inter"
).split()


def build_dag(varnames, exprs):
    """Return (DAG, variables, dependents) for the given expression strings."""
    DAG = mc.FFGraph()
    variables = [mc.FFVar(DAG, n) for n in varnames]
    env = {n: v for n, v in zip(varnames, variables)}
    env.update({f: getattr(mc, f) for f in _FUNCS if hasattr(mc, f)})
    deps = []
    for e in exprs:
        value = eval(e, {"__builtins__": {}}, env)  # noqa: S307 - user's own expression
        if not isinstance(value, mc.FFVar):
            raise TypeError(f"expression {e!r} did not evaluate to an FFVar")
        deps.append(value)
    return DAG, variables, deps


# ---------------------------------------------------------------------------
# Incidence assembly
# ---------------------------------------------------------------------------
def incidence(sg, with_outputs=False):
    """Assemble the incidence pattern of a subgraph from sg.tape().

    tape() is the wrapper-free view: one (op_type, varin_ids, varout_ids)
    tuple per operation, where each id is the (FFVar.TYPE, index) pair that
    uniquely identifies a node -- the index alone restarts per TYPE, so
    x is (VAR, 0) and the first auxiliary is (AUX, 0).

    Returns a dict with:
        pairs      sorted [(i, j)] : operand i is an INPUT of operation j
        out_pairs  sorted [(i, j)] : operand i is an OUTPUT of operation j
        row_id     [node_id] indexed by row
        def_col    [j or None] indexed by row: operation defining that node
        op_types   [FFOp.TYPE] indexed by column
        repeats    [(i, j, multiplicity)] for operands appearing >1x in one op
        dangling   [(node_id, j)] inputs used before being defined
    """
    tape = sg.tape()

    row_of = {}          # node id -> row index
    row_id = []          # row index -> node id
    def_col = []         # row index -> defining operation column (or None)
    pairs, out_pairs, repeats, dangling = set(), set(), [], []
    op_types = []

    def row(node_id):
        i = row_of.get(node_id)
        if i is None:
            i = row_of[node_id] = len(row_id)
            row_id.append(node_id)
            def_col.append(None)
        return i

    for j, (op_type, varin, varout) in enumerate(tape):
        op_types.append(op_type)

        # Outputs first, so rows are numbered in order of definition.
        for vid in varout:
            i = row(vid)
            if def_col[i] is None:
                def_col[i] = j
            if with_outputs:
                out_pairs.add((i, j))

        # Inputs. In evaluation order every input must already be defined.
        seen = {}
        for vid in varin:
            if vid not in row_of:
                dangling.append((vid, j))
            i = row(vid)
            seen[i] = seen.get(i, 0) + 1
            pairs.add((i, j))
        # FFOp::varin may hold a repeated operand (FFOp::dupfirst exists for
        # that case); the incidence pattern is boolean, so collapse it here
        # and report the multiplicity rather than emitting a duplicate pair.
        for i, k in seen.items():
            if k > 1:
                repeats.append((i, j, k))

    return {
        "pairs": sorted(pairs),
        "out_pairs": sorted(out_pairs),
        "row_id": row_id,
        "def_col": def_col,
        "op_types": op_types,
        "repeats": repeats,
        "dangling": dangling,
    }


def incidence_via_ops(sg):
    """Same input-incidence pattern, taken through sg.l_op / FFOp / FFVar.id.

    Slower -- it materialises one FFOp wrapper per operation and one FFVar
    wrapper per operand -- and used here only to cross-check incidence().
    """
    row_of, pairs = {}, set()

    def row(node_id):
        return row_of.setdefault(node_id, len(row_of))

    for j, op in enumerate(sg.l_op):
        for v in op.varout:
            row(v.id)
        for v in op.varin:
            pairs.add((row(v.id), j))
    return sorted(pairs)


def node_labels(sg):
    """Row labels (FFVar names) keyed by node id, via one pass over l_op."""
    label = {}
    for op in sg:
        for v in list(op.varout) + list(op.varin):
            label.setdefault(v.id, str(v))
    return label


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def fmt_id(node_id):
    return f"{str(node_id[0]).split('.')[-1]}{node_id[1]}"


def fmt_type(op_type):
    return str(op_type).split(".")[-1]


def print_map(inc, labels, max_rows=60, max_cols=60):
    nrow, ncol = len(inc["row_id"]), len(inc["op_types"])
    if nrow > max_rows or ncol > max_cols:
        print(f"  (pattern map suppressed: {nrow}x{ncol} exceeds "
              f"{max_rows}x{max_cols})")
        return
    cols = set(inc["pairs"])
    outs = set(inc["out_pairs"])
    width = max((len(labels.get(v, fmt_id(v))) for v in inc["row_id"]), default=4)
    print("      " + " ".join(f"{j%10}" for j in range(ncol)))
    print("      " + "-" * (2 * ncol - 1))
    for i, vid in enumerate(inc["row_id"]):
        name = labels.get(vid, fmt_id(vid))
        cells = []
        for j in range(ncol):
            if (i, j) in cols:
                cells.append("X")
            elif (i, j) in outs:
                cells.append("o")
            #elif inc["def_col"][i] == j:
            #    cells.append("d")
            else:
                cells.append(".")
        print(f"  {name:>{width}} |" + " ".join(cells))
    #print("      " + "-" * (2 * ncol - 1))
    #print("      " + " ".join(fmt_type(t)[0] for t in inc["op_types"]))
    legend = "X = operand of op"#   d = defined by op"
    print(f"\n  {legend}" + ("   o = output edge" if outs else ""))


def verify(sg, inc, via_ops):
    """Checks on the assembled pattern. Returns the number of failures."""
    fails = 0

    def check(ok, what):
        nonlocal fails
        print(f"  [{' ok ' if ok else 'FAIL'}] {what}")
        if not ok:
            fails += 1

    check(inc["pairs"] == via_ops,
          "tape() and l_op routes give the same pattern")

    # Fan-in per column must match the operation's own varin count.
    want = [len({v.id for v in op.varin}) for op in sg.l_op]
    got = [0] * len(want)
    for _, j in inc["pairs"]:
        got[j] += 1
    check(want == got, "column counts match each operation's distinct fan-in")

    check(not inc["dangling"],
          "every operand is defined before it is used (evaluation order is "
          "topological)")

    # Stronger form of the same invariant, on the assembled pattern.
    ok = all(inc["def_col"][i] is not None and inc["def_col"][i] < j
             for i, j in inc["pairs"])
    check(ok, "def_col[i] < j for every incidence pair")

    # Rows and columns are dense and consistent.
    nrow, ncol = len(inc["row_id"]), len(inc["op_types"])
    check(ncol == len(sg), "one column per operation in the subgraph")
    check(len(set(inc["row_id"])) == nrow, "node ids are unique per row")
    check(all(0 <= i < nrow and 0 <= j < ncol for i, j in inc["pairs"]),
          "all pairs are in range")

    # Rows defined by a VAR operation should be exactly the subgraph's own
    # list of independents.
    indep = {str(v) for v in sg.v_indep}
    labels = node_labels(sg)
    from_pattern = {labels[inc["row_id"][i]]
                    for i in range(nrow)
                    if inc["def_col"][i] is not None
                    and fmt_type(inc["op_types"][inc["def_col"][i]]) == "VAR"}
    check(from_pattern == indep,
          f"VAR-defined rows match sg.v_indep ({sorted(indep)})")

    return fails


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vars", help="comma-separated variable names")
    p.add_argument("--expr", action="append", help="dependent expression "
                   "(repeat for several); uses pymcpp functions")
    p.add_argument("--example", choices=sorted(EXAMPLES), default="default",
                   help="built-in expression set (default: default)")
    p.add_argument("--with-outputs", action="store_true",
                   help="also list (operand, operation) output edges")
    p.add_argument("--no-map", action="store_true", help="skip the pattern map")
    args = p.parse_args()

    if args.expr:
        varnames = (args.vars or DEFAULT_VARS).split(",")
        exprs = args.expr
    else:
        v, exprs = EXAMPLES[args.example]
        varnames = (args.vars or v).split(",")
    varnames = [n.strip() for n in varnames if n.strip()]

    if not hasattr(mc.FFSubgraph, "tape"):
        sys.exit("This pymcpp lacks the FFSubgraph binder additions "
                 "(l_op/tape/types). Apply the binder patch and rebuild.")

    # --- DAG ---------------------------------------------------------------
    DAG, variables, deps = build_dag(varnames, exprs)
    print("variables  :", ", ".join(varnames))
    for k, e in enumerate(exprs):
        print(f"dependent {k}: {e}")

    # --- subgraph ----------------------------------------------------------
    sg = DAG.subgraph(deps)
    print("\nsubgraph   :", repr(sg))

    # --- incidence ---------------------------------------------------------
    inc = incidence(sg, with_outputs=args.with_outputs)
    labels = node_labels(sg)
    nrow, ncol = len(inc["row_id"]), len(inc["op_types"])

    print(f"\nincidence  : {nrow} operands x {ncol} operations, "
          f"{len(inc['pairs'])} nonzeros "
          f"(density {len(inc['pairs'])/(nrow*ncol):.3f})")

    print("\nrows (operand i):")
    for i, vid in enumerate(inc["row_id"]):
        d = inc["def_col"][i]
        origin = "independent" if d is not None and \
            fmt_type(inc["op_types"][d]) == "VAR" else f"defined by op {d}"
        print(f"  {i:3d}  {labels.get(vid, fmt_id(vid)):>6s}  "
              f"[{fmt_id(vid):>7s}]  {origin}")

    print("\ncolumns (operation j):")
    for j, t in enumerate(inc["op_types"]):
        ins = [i for i, jj in inc["pairs"] if jj == j]
        print(f"  {j:3d}  {fmt_type(t):<8s} operands={ins}")

    print(f"\nincidence pairs (i, j) -- operand i participates in operation j:")
    print(" ", inc["pairs"])

    if args.with_outputs:
        print("\noutput edges (i, j) -- operand i is produced by operation j:")
        print(" ", inc["out_pairs"])

    if inc["repeats"]:
        print("\nrepeated operands (i, j, multiplicity), collapsed in the "
              "pattern:")
        print(" ", inc["repeats"])

    if not args.no_map:
        print("\npattern:")
        print_map(inc, labels)

    # --- verification ------------------------------------------------------
    print("\nverification:")
    fails = verify(sg, inc, incidence_via_ops(sg))
    print(f"\n{'OK' if not fails else str(fails) + ' CHECK(S) FAILED'}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
