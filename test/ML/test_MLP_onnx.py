# Generates the ONNX models read by test_MLP_onnx.cpp and test_MLPREG.cpp, and the values of the
# networks they hold at sample points, computed by ONNX Runtime:
#   python3 test_MLP_onnx.py
# needs torch, onnx, onnxruntime, scikit-learn and skl2onnx. Models that mc::MLPIO::read_onnx reads
# (onnx_*.onnx, values in onnx_values.txt), and models that it must refuse (onnx_bad_*.onnx).
import warnings
import numpy as np
import torch, torch.nn as nn
import onnx
from onnx import helper, numpy_helper, TensorProto
import onnxruntime as ort

warnings.filterwarnings("ignore")
torch.manual_seed(3)
rng = np.random.default_rng(3)
values = []


class Net(nn.Module):
    def __init__(self, seq):
        super().__init__()
        self.net = seq

    def forward(self, x):
        return self.net(x)


def record(name, nin, npts=3):
    """Values of the model at npts points, by ONNX Runtime, in the input shape of the model."""
    sess = ort.InferenceSession(name, providers=["CPUExecutionProvider"])
    inp = sess.get_inputs()[0]
    dtype = np.float64 if inp.type == "tensor(double)" else np.float32
    for k in range(npts):
        x = (0.8 * np.cos(np.arange(nin) + 1.3 * k + 0.4)).astype(dtype)
        shape = [d if isinstance(d, int) else 1 for d in inp.shape]
        y = sess.run(None, {inp.name: x.reshape(shape)})[0].astype(np.float64).ravel()
        values.append(name + " " + str(nin) + " " + str(len(y)) + " "
                      + " ".join(repr(float(v)) for v in x) + " "
                      + " ".join(repr(float(v)) for v in y))


def torch_export(name, seq, nin, batch=True, double=False):
    m = Net(seq).eval()
    if double:
        m = m.double()
    x = torch.zeros((1, nin) if batch else (nin,), dtype=torch.float64 if double else torch.float32)
    torch.onnx.export(m, (x,), name, input_names=["x"], output_names=["y"], dynamo=False)
    record(name, nin)


# PyTorch, TorchScript exporter: Gemm (batch of inputs), MatMul + Add (a single input)
torch_export("onnx_torch_tanh.onnx", nn.Sequential(nn.Linear(3, 8), nn.Tanh(), nn.Linear(8, 8), nn.Tanh(), nn.Linear(8, 2)), 3)
torch_export("onnx_torch_relu.onnx", nn.Sequential(nn.Linear(4, 6), nn.ReLU(), nn.Dropout(0.2), nn.Linear(6, 1), nn.Sigmoid()), 4)
torch_export("onnx_torch_vector.onnx", nn.Sequential(nn.Linear(3, 5), nn.Sigmoid(), nn.Linear(5, 2)), 3, batch=False)
torch_export("onnx_torch_double.onnx", nn.Sequential(nn.Flatten(), nn.Linear(2, 7), nn.Tanh(), nn.Linear(7, 1)), 2, double=True)
torch_export("onnx_torch_nobias.onnx", nn.Sequential(nn.Linear(3, 4, bias=False), nn.ReLU(), nn.Linear(4, 2)), 3)
bn = nn.Sequential(nn.Linear(3, 5), nn.BatchNorm1d(5), nn.Tanh(), nn.Linear(5, 1))
with torch.no_grad():
    bn[1].running_mean.copy_(torch.randn(5) * 0.3)
    bn[1].running_var.copy_(torch.rand(5) + 0.5)
    bn[1].weight.copy_(torch.randn(5))
    bn[1].bias.copy_(torch.randn(5) * 0.2)
# without constant folding, so that the BatchNormalization node stays in the graph
m = Net(bn).eval()
torch.onnx.export(m, (torch.zeros(1, 3),), "onnx_torch_batchnorm.onnx", input_names=["x"], output_names=["y"],
                  do_constant_folding=False, dynamo=False)
record("onnx_torch_batchnorm.onnx", 3)

# scikit-learn MLPRegressor, by skl2onnx: MatMul, Add, Relu / Tanh / Sigmoid
try:
    from sklearn.neural_network import MLPRegressor
    from skl2onnx import to_onnx
    X = rng.uniform(-1, 1, (200, 3))
    Y = np.column_stack([np.sin(X[:, 0]) + X[:, 1] * X[:, 2], np.cos(X[:, 1])])
    for act in ["relu", "tanh", "logistic"]:
        reg = MLPRegressor(hidden_layer_sizes=(6, 5), activation=act, max_iter=200, random_state=0).fit(X, Y)
        name = "onnx_sklearn_" + act + ".onnx"
        onx = to_onnx(reg, X[:1].astype(np.float32), target_opset=17)
        onnx.save(onx, name)
        record(name, 3)
except ImportError as e:
    print("scikit-learn models not generated:", e)


# Graphs written by hand
def save_graph(name, nodes, inits, nin, nout, dtype=TensorProto.FLOAT, inputs=None, external=False):
    g = helper.make_graph(nodes, "g", inputs or [helper.make_tensor_value_info("x", dtype, [1, nin])],
                          [helper.make_tensor_value_info("y", dtype, [1, nout])], inits)
    model = helper.make_model(g, opset_imports=[helper.make_opsetid("", 17)])
    model.ir_version = 8
    if external:
        onnx.save_model(model, name, save_as_external_data=True, all_tensors_to_one_file=True,
                        location=name + ".data", size_threshold=0)
    else:
        onnx.save(model, name)


def T(name, arr, dtype=np.float32):
    return numpy_helper.from_array(np.asarray(arr, dtype=dtype), name)


W1, b1 = rng.normal(size=(3, 4)), rng.normal(size=4)   # [nin, nneu]: MatMul, or Gemm with transB = 0
W2, b2 = rng.normal(size=(2, 4)), rng.normal(size=2)   # [nout, nprv]: Gemm with transB = 1
# Gemm with transB = 0, alpha and beta, the bias as a Constant node, a Cast and an Identity
save_graph("onnx_hand_gemm.onnx",
           [helper.make_node("Constant", [], ["c1"], value=T("c1v", b1)),
            helper.make_node("Gemm", ["x", "W1", "c1"], ["h1"], alpha=0.5, beta=2.0, transB=0),
            helper.make_node("Cast", ["h1"], ["h1c"], to=TensorProto.FLOAT),
            helper.make_node("Tanh", ["h1c"], ["a1"]),
            helper.make_node("Identity", ["a1"], ["a1i"]),
            helper.make_node("Gemm", ["a1i", "W2", "b2"], ["y"], transB=1)],
           [T("W1", W1), T("W2", W2), T("b2", b2)], 3, 2)
record("onnx_hand_gemm.onnx", 3)
# MatMul followed by Add with the bias first, in double precision
save_graph("onnx_hand_matmul.onnx",
           [helper.make_node("MatMul", ["x", "W1"], ["m1"]),
            helper.make_node("Add", ["b1", "m1"], ["h1"]),
            helper.make_node("Sigmoid", ["h1"], ["a1"]),
            helper.make_node("MatMul", ["a1", "W2t"], ["m2"]),
            helper.make_node("Add", ["m2", "b2"], ["y"])],
           [T("W1", W1, np.float64), T("b1", b1, np.float64), T("W2t", W2.T, np.float64), T("b2", b2, np.float64)],
           3, 2, dtype=TensorProto.DOUBLE)
record("onnx_hand_matmul.onnx", 3)
# float16 weights
save_graph("onnx_hand_float16.onnx",
           [helper.make_node("Cast", ["x"], ["x16"], to=TensorProto.FLOAT16),
            helper.make_node("Gemm", ["x16", "W1", "b1"], ["h1"]),
            helper.make_node("Relu", ["h1"], ["a1"]),
            helper.make_node("Gemm", ["a1", "W2", "b2"], ["y16"], transB=1),
            helper.make_node("Cast", ["y16"], ["y"], to=TensorProto.FLOAT)],
           [T("W1", W1, np.float16), T("b1", b1, np.float16), T("W2", W2, np.float16), T("b2", b2, np.float16)], 3, 2)
record("onnx_hand_float16.onnx", 3)

# Models to refuse
save_graph("onnx_bad_leakyrelu.onnx",
           [helper.make_node("Gemm", ["x", "W1", "b1"], ["h1"]), helper.make_node("LeakyRelu", ["h1"], ["a1"]),
            helper.make_node("Gemm", ["a1", "W2", "b2"], ["y"], transB=1)],
           [T("W1", W1), T("b1", b1), T("W2", W2), T("b2", b2)], 3, 2)
save_graph("onnx_bad_skip.onnx",
           [helper.make_node("Gemm", ["x", "W1", "b1"], ["h1"]), helper.make_node("Tanh", ["h1"], ["a1"]),
            helper.make_node("Add", ["a1", "h1"], ["s1"]),
            helper.make_node("Gemm", ["s1", "W2", "b2"], ["y"], transB=1)],
           [T("W1", W1), T("b1", b1), T("W2", W2), T("b2", b2)], 3, 2)
save_graph("onnx_bad_twoinputs.onnx",
           [helper.make_node("Gemm", ["x", "W1", "b1"], ["h1"]), helper.make_node("Add", ["h1", "z"], ["y"])],
           [T("W1", W1), T("b1", b1)], 3, 4,
           inputs=[helper.make_tensor_value_info("x", TensorProto.FLOAT, [1, 3]),
                   helper.make_tensor_value_info("z", TensorProto.FLOAT, [1, 4])])
save_graph("onnx_bad_size.onnx",
           [helper.make_node("Gemm", ["x", "W1", "b1"], ["h1"]),
            helper.make_node("Gemm", ["h1", "W2", "b2"], ["y"], transB=0)],
           [T("W1", W1), T("b1", b1), T("W2", W2), T("b2", b2)], 3, 2)
save_graph("onnx_bad_external.onnx",
           [helper.make_node("Gemm", ["x", "W1", "b1"], ["y"])],
           [T("W1", W1), T("b1", b1)], 3, 4, external=True)
with open("onnx_bad_truncated.onnx", "wb") as f:
    f.write(open("onnx_hand_gemm.onnx", "rb").read()[:150])

with open("onnx_values.txt", "w") as f:
    f.write("\n".join(values) + "\n")
print(len(values), "values written to onnx_values.txt")
