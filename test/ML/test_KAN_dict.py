# Generates the dictionaries of torch.save read by test_KAN.cpp without Torch, and the values and
# bounds of the networks they hold:
#   python3 test_KAN_dict.py
# needs torch and numpy. Networks that mc::MLPIO::read_dict reads (dict_*.pt, values and bounds in
# dict_values.txt), and files that it must refuse (dict_bad_*.pt).
# A dictionary holds model_type, config, state_dict, normalization and bounds: a Kolmogorov-Arnold
# network with edges c + v u + sum_r a_r relu(w_r u + b_r) ("KANmodel"), or a dense network
# ("ANNmodel"). The values are computed here in double precision, and the bounds of a KAN are the
# exact ranges of its edges on the ranges of their inputs, summed over the edges entering a node.
import zipfile
import numpy as np
import torch

rng = np.random.default_rng(7)
values = []


def kan_layers(state, widths):
    return [{"c": state[f"layers.{l}.edge_bias"].double().numpy(),
             "v": state[f"layers.{l}.linear_weight"][..., 0].double().numpy(),
             "a": state[f"layers.{l}.relu_output_weight"].double().numpy(),
             "w": state[f"layers.{l}.relu_weight"].double().numpy(),
             "b": state[f"layers.{l}.relu_bias"].double().numpy()} for l in range(len(widths) - 1)]


def kan_forward(layers, x):
    for L in layers:
        pre = L["w"] * x[:, None, :, None] + L["b"]
        x = (L["c"] + L["v"] * x[:, None, :] + (np.maximum(pre, 0.0) * L["a"]).sum(-1)).sum(-1)
    return x


def edge_range(L, j, i, lo, hi):
    w, b, a = L["w"][j, i], L["b"][j, i], L["a"][j, i]
    valid = (w != 0) & (a != 0)
    knots = -b[valid] / w[valid]
    knots = knots[(knots > lo) & (knots < hi)]
    x = np.unique(np.concatenate(([lo, hi], knots)))
    y = L["c"][j, i] + L["v"][j, i] * x + (np.maximum(np.outer(x, w) + b, 0.0) * a).sum(1)
    return y.min(), y.max()


def kan_bounds(layers, box):
    """Node bounds of each layer from the exact ranges of the edges."""
    out, node = [], np.array(box, dtype=float)
    for L in layers:
        nout, nin = L["c"].shape
        edge = np.array([[edge_range(L, j, i, node[i, 0], node[i, 1]) for i in range(nin)] for j in range(nout)])
        node = edge.sum(1)
        out.append({"edge": torch.tensor(edge), "node": torch.tensor(node)})
    return out


def scale(norm, side, n):
    if norm[side] is None:
        return np.zeros(n), np.ones(n)
    return np.array(norm[side]["lower"], dtype=float), np.array(norm[side]["upper"], dtype=float)


def record(name, widths, norm, forward, bound=None, npts=4):
    """Lines 'value name nin nout x.. z.. y..' (raw inputs, raw outputs, outputs in the scaled
    coordinates) and 'bound name nout lower.. upper..' (outputs in the scaled coordinates on the
    whole input box)."""
    nin, nout = widths[0], widths[-1]
    xl, xu = scale(norm, "x", nin)
    zl, zu = scale(norm, "z", nout)
    for k in range(npts):
        u = 0.5 + 0.45 * np.cos(np.arange(nin) + 1.3 * k + 0.4)     # scaled inputs in [0.05,0.95]
        y = forward(u[None, :])[0]
        values.append(" ".join(["value", name, str(nin), str(nout)]
                               + [repr(float(v)) for v in xl + (xu - xl) * u]
                               + [repr(float(v)) for v in zl + (zu - zl) * y]
                               + [repr(float(v)) for v in y]))
    if bound is not None:
        values.append(" ".join(["bound", name, str(nout)] + [repr(float(v)) for v in bound[:, 0]]
                               + [repr(float(v)) for v in bound[:, 1]]))


def make_kan(name, widths, intervals, norm, zero_terms=False):
    state = {}
    for l in range(len(widths) - 1):
        no, ni, R = widths[l + 1], widths[l], intervals[l] - 1
        w = rng.uniform(-2, 2, (no, ni, R))
        b = rng.uniform(-1.5, 1.5, (no, ni, R))
        a = rng.uniform(-1, 1, (no, ni, R))
        if zero_terms:              # terms without breakpoint, and two terms with the same one
            a[0, 0, 0] = 0.0
            w[0, 1, 0] = 0.0
            w[1, 0, 1], b[1, 0, 1] = 2.0 * w[1, 0, 0], 2.0 * b[1, 0, 0]
        state[f"layers.{l}.relu_weight"] = torch.tensor(w, dtype=torch.float32)
        state[f"layers.{l}.relu_bias"] = torch.tensor(b, dtype=torch.float32)
        state[f"layers.{l}.relu_output_weight"] = torch.tensor(a, dtype=torch.float32)
        state[f"layers.{l}.linear_weight"] = torch.tensor(rng.uniform(-1, 1, (no, ni, 1)), dtype=torch.float32)
        state[f"layers.{l}.edge_bias"] = torch.tensor(rng.uniform(-0.5, 0.5, (no, ni)), dtype=torch.float32)
    layers = kan_layers(state, widths)
    box = [[0.0, 1.0]] * widths[0]
    lb = kan_bounds(layers, box)
    payload = {"model_type": "KANmodel",
               "config": {"widths": widths, "intervals": intervals},
               "state_dict": state,
               "normalization": norm,
               "bounds": {"CPWL": {"input_box": box, "constraints": None,
                                   "bounds": {"input": torch.tensor(box), "layers": lb}}}}
    torch.save(payload, name)
    record(name, widths, norm, lambda u: kan_forward(layers, u), lb[-1]["node"].numpy())
    return payload


def make_ann(name, widths, activ, norm, named=True, with_bounds=True):
    g = {"relu": lambda z: np.maximum(z, 0), "tanh": np.tanh, "sigmoid": lambda z: 1 / (1 + np.exp(-z))}[activ]
    state = {}
    for l in range(len(widths) - 1):
        state[f"layers.{l}.weight"] = torch.tensor(rng.uniform(-1, 1, (widths[l + 1], widths[l])), dtype=torch.float32)
        state[f"layers.{l}.bias"] = torch.tensor(rng.uniform(-0.5, 0.5, widths[l + 1]), dtype=torch.float32)

    def forward(u):
        for l in range(len(widths) - 1):
            u = u @ state[f"layers.{l}.weight"].double().numpy().T + state[f"layers.{l}.bias"].double().numpy()
            if l < len(widths) - 2:
                u = g(u)
        return u

    # interval bounds of the pre-activations and activations of each layer
    lb, lo, hi = [], np.zeros(widths[0]), np.ones(widths[0])
    for l in range(len(widths) - 1):
        W, b = state[f"layers.{l}.weight"].double().numpy(), state[f"layers.{l}.bias"].double().numpy()
        pl = b + np.minimum(W * lo, W * hi).sum(1)
        pu = b + np.maximum(W * lo, W * hi).sum(1)
        lo, hi = (g(pl), g(pu)) if l < len(widths) - 2 else (pl, pu)
        lb.append({"pre": torch.tensor(np.stack([pl, pu], 1)), "post": torch.tensor(np.stack([lo, hi], 1))})
    config = {"widths": widths}
    if named:
        config["activation"] = activ
    payload = {"model_type": "ANNmodel", "config": config, "state_dict": state, "normalization": norm}
    if with_bounds:
        payload["bounds"] = {"IA": {"input_box": [[0.0, 1.0]] * widths[0], "constraints": None,
                                    "bounds": {"layers": lb}}}
    torch.save(payload, name)
    record(name, widths, norm, forward, lb[-1]["post"].numpy())
    return payload


# Kolmogorov-Arnold networks: both sides scaled; inputs not scaled, two outputs, terms without
# breakpoint; no normalization entry at all, a single layer
xz = {"method": "min_max", "x": {"columns": ["x1", "x2"], "lower": [-1.0, 0.5], "upper": [2.0, 3.0]},
      "z": {"columns": ["z1"], "lower": [1.0], "upper": [4.0]}}
make_kan("dict_kan_scaled.pt", [2, 3, 1], [4, 3], xz)
make_kan("dict_kan_two.pt", [3, 2, 2], [3, 4],
         {"method": "min_max", "x": None, "z": {"columns": ["z1", "z2"], "lower": [0.0, -1.0], "upper": [5.0, 1.0]}},
         zero_terms=True)
p = make_kan("dict_kan_single.pt", [2, 2], [5], {"method": "min_max", "x": None, "z": None})

# Dense networks: activation named in config; recognised from the bounds; neither (read with the
# activation given)
make_ann("dict_ann_tanh.pt", [3, 4, 2], "tanh", {"method": "min_max", "x": None, "z": None})
make_ann("dict_ann_relu.pt", [2, 5, 5, 1], "relu", xz, named=False)
make_ann("dict_ann_sigmoid.pt", [2, 3, 1], "sigmoid", xz, named=False)
make_ann("dict_ann_plain.pt", [2, 3, 1], "relu", xz, named=False, with_bounds=False)

# Files that must be refused: unknown type of model; tensor of the wrong shape; missing tensor;
# scaling with upper not above lower; a tensor alone, not a dictionary; an unknown activation;
# truncated file; compressed archive
q = dict(p); q["model_type"] = "RBFmodel"; torch.save(q, "dict_bad_type.pt")
q = dict(p); q["state_dict"] = dict(p["state_dict"])
q["state_dict"]["layers.0.relu_bias"] = q["state_dict"]["layers.0.relu_bias"][:, :, :2].clone()
torch.save(q, "dict_bad_shape.pt")
q = dict(p); q["state_dict"] = {k: v for k, v in p["state_dict"].items() if "edge_bias" not in k}
torch.save(q, "dict_bad_missing.pt")
q = dict(p); q["normalization"] = {"method": "min_max", "x": {"lower": [0.0, 1.0], "upper": [1.0, 1.0]}, "z": None}
torch.save(q, "dict_bad_scaling.pt")
torch.save(torch.zeros(3), "dict_bad_tensor.pt")
q = torch.load("dict_ann_tanh.pt", weights_only=True); q["config"]["activation"] = "softplus"
torch.save(q, "dict_bad_activ.pt")
raw = open("dict_kan_single.pt", "rb").read()
open("dict_bad_truncated.pt", "wb").write(raw[:len(raw) // 2])
with zipfile.ZipFile("dict_kan_single.pt") as zin, zipfile.ZipFile("dict_bad_deflated.pt", "w", zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        zout.writestr(item.filename, zin.read(item.filename))

open("dict_values.txt", "w").write("\n".join(values) + "\n")
print(len(values), "lines in dict_values.txt")
