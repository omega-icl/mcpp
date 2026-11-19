#!/usr/bin/env python
# coding: utf-8

# # Training a neural network in an MC++ DAG
# 
# This tutorial trains a multilayer perceptron (MLP) whose value and gradient are computed by MC++, and then
# uses the trained network as an ordinary MC++ operation. It goes through five steps:
# 
# 1. sample a ground-truth function -- MATLAB's `peaks` -- at Sobol points;
# 2. define the network with ONNX, read it as an `MLPREG` -- a network whose **weights are variables** -- and
#    embed it in a DAG with an `FFMLPREG` operation;
# 3. train the weights with PyTorch's Adam optimizer, the loss and its gradient in the weights being computed by
#    evaluating and differentiating the MC++ DAG;
# 4. define an `MLP` -- a network whose weights are **data** -- from the trained weights;
# 5. embed the `MLP` in a DAG with an `FFMLP` operation and evaluate it: values, derivatives, bounds and
#    relaxations.
# 
# `MLPREG` and `FFMLPREG` mirror `MLP` and `FFMLP`: the operands of an `MLPREG` are the network inputs followed by
# its weights, so that the same network can be evaluated, differentiated and relaxed in its weights. PyTorch is
# used here for the definition of the network and for the optimizer only.
# 
# **Requirements.** pymcpp 5.2 or later, PyTorch, SciPy and Matplotlib; the `onnx` package is used only if the
# installed PyTorch cannot export ONNX files itself.

# In[ ]:


import time
import numpy as np
import torch
import matplotlib.pyplot as plt
from scipy.stats import qmc
import pymcpp as mc

print("pymcpp", mc.__version__, "| torch", torch.__version__)


# ## 1. Ground truth and training data
# 
# MATLAB's `peaks` function on $[-3,3]^2$:
# $$
# f(x,y) = 3(1-x)^2 e^{-x^2-(y+1)^2} - 10\left(\frac{x}{5}-x^3-y^5\right) e^{-x^2-y^2} - \frac{1}{3} e^{-(x+1)^2-y^2}
# $$
# is sampled at $2^8 = 256$ points of a scrambled Sobol sequence, which fill the domain more evenly than random
# points. The network is trained on scaled data: inputs in $[-1,1]^2$ and outputs of zero mean and unit
# variance. The scaling is undone in the DAG of the trained network (step 5).

# In[ ]:


def peaks(x, y):
    return (3 * (1 - x)**2 * np.exp(-x**2 - (y + 1)**2)
            - 10 * (x / 5 - x**3 - y**5) * np.exp(-x**2 - y**2)
            - np.exp(-(x + 1)**2 - y**2) / 3)

lo, up = -3., 3.
sobol = qmc.Sobol(d=2, scramble=True, seed=1)
P = qmc.scale(sobol.random_base2(m=8), [lo, lo], [up, up])   # 256 points in [-3,3]^2
F = peaks(P[:, 0], P[:, 1])

xscale = up                       # inputs scaled to [-1,1]
fmean, fstd = F.mean(), F.std()   # outputs scaled to zero mean, unit variance
Xs, Ts = P / xscale, (F - fmean) / fstd
print(f"{len(P)} samples, f in [{F.min():.3f}, {F.max():.3f}]")

g = np.linspace(lo, up, 121)
GX, GY = np.meshgrid(g, g)
FG = peaks(GX, GY)
fig, ax = plt.subplots(figsize=(5, 4.2))
cs = ax.contourf(GX, GY, FG, levels=30, cmap="viridis")
ax.plot(P[:, 0], P[:, 1], "w.", ms=3)
fig.colorbar(cs, ax=ax)
ax.set_title("peaks and the Sobol samples"); ax.set_xlabel("x"); ax.set_ylabel("y")
plt.show()


# ## 2. The network: from ONNX to an `MLPREG` in a DAG
# 
# The architecture is defined in PyTorch -- 2 inputs, 3 hidden layers of 40 neurons with tanh activation (ReLU
# works as well), 1 linear output -- and exported to an ONNX file. `MLPREG.read_data` reads its structure, and its
# initial weights, which PyTorch drew at random, into `MLPREG.weights`. The ONNX file could as well come from Keras,
# scikit-learn or any other framework that exports ONNX.

# In[ ]:


ACTIV = "tanh"   # or "relu"
act = {"tanh": torch.nn.Tanh, "relu": torch.nn.ReLU}[ACTIV]

torch.manual_seed(0)
net = torch.nn.Sequential(torch.nn.Linear(2, 40), act(),
                          torch.nn.Linear(40, 40), act(),
                          torch.nn.Linear(40, 40), act(),
                          torch.nn.Linear(40, 1))

def export_onnx(net, path, nin):
    # ONNX file of a Sequential of Linear layers and activations
    try:
        torch.onnx.export(net, (torch.zeros(1, nin),), path, input_names=["x"], output_names=["y"],
                          dynamo=False)
    except Exception:
        # PyTorch without its TorchScript exporter: the graph written with the onnx package, one Gemm node
        # per Linear layer and one node per activation
        import onnx
        from onnx import helper, numpy_helper, TensorProto
        nodes, inits, cur = [], [], "x"
        for k, m in enumerate(net):
            out = f"h{k}"
            if isinstance(m, torch.nn.Linear):
                inits += [numpy_helper.from_array(m.weight.detach().numpy(), f"W{k}"),
                          numpy_helper.from_array(m.bias.detach().numpy(), f"b{k}")]
                nodes.append(helper.make_node("Gemm", [cur, f"W{k}", f"b{k}"], [out], transB=1))
            else:
                nodes.append(helper.make_node({"Tanh": "Tanh", "ReLU": "Relu"}[type(m).__name__], [cur], [out]))
            cur = out
        nodes[-1].output[0] = "y"
        graph = helper.make_graph(nodes, "mlp", [helper.make_tensor_value_info("x", TensorProto.FLOAT, [1, nin])],
                                  [helper.make_tensor_value_info("y", TensorProto.FLOAT, [1, 1])], inits)
        onnx.save(helper.make_model(graph, opset_imports=[helper.make_opsetid("", 17)]), path)

export_onnx(net, "peaks_mlp.onnx", 2)

REG = mc.MLPREG()
REG.read_data("peaks_mlp.onnx")
print("inputs:", REG.nin, "| layers (neurons, activation):", REG.layers, "| outputs:", REG.nout)
print("weights and biases:", REG.nwei, "| operands (inputs, then weights):", REG.nvar)


# The weights of an `MLPREG` are ordered layer after layer and, in each layer, neuron after neuron: the bias of
# the neuron followed by its weights, one per neuron of the previous layer. These are the rows `[bias, w_1, ...,
# w_n]` that `MLP.set_data` takes, end to end. The weights read from the ONNX file are those of the PyTorch
# network in this ordering:

# In[ ]:


def flatten(net):
    # weights of a Sequential in the ordering of MLPREG
    return np.concatenate([np.hstack([m.bias.detach().numpy()[:, None], m.weight.detach().numpy()]).ravel()
                           for m in net if isinstance(m, torch.nn.Linear)])

print("largest difference:", np.max(np.abs(flatten(net) - np.array(REG.weights))))


# The DAG has variables for the scaled inputs $X$, the weights $W$ and the scaled target $T$. The `FFMLPREG`
# operation gives the network output $Y(X,W)$, and the loss at one data point is the squared error
# $\ell = (Y-T)^2$. Its gradient in the weights is a symbolic derivative in the same DAG (`bdiff`). With
# `MLPREG.options.ADJOINT`, the backward differentiation inserts the *adjoint* of the network -- one backward pass
# through the layers, an `AdjMLPREG` operation -- instead of its whole Jacobian, which has as many entries as there
# are weights: always the cheaper here.

# In[ ]:


REG.options.ADJOINT = True

DAG = mc.FFGraph()
X = DAG.add_vars(REG.nin, "X")    # scaled inputs
W = DAG.add_vars(REG.nwei, "W")   # weights and biases
T = DAG.add_var("T")              # scaled target

OpREG = mc.FFMLPREG()
Y = OpREG(X, W, REG)[0]           # network output
Loss = mc.sqr(Y - T)              # squared error at one data point

_, cols, dLoss = DAG.bdiff([Loss], W)     # gradient in the weights: sparse triplets
assert list(cols) == list(range(REG.nwei))
deps = [Loss] + list(dLoss)
sgLoss = DAG.subgraph(deps)       # the operations to evaluate, extracted once

print("external operations:", sorted({mc.FFOp.name(op).split("[")[0] for op in sgLoss.l_op if op.type >= int(op.EXTERN)}))
print("operations in the subgraph:", len(sgLoss.l_op))


# The data points are *scenarios* of the variables $X$ and $T$, and the weights are shared by all of them:
# `FFGraph.veval` evaluates the loss and its gradient at a whole batch of points in one call. As a check, the mean
# loss and gradient over a batch agree with those of PyTorch's own automatic differentiation, in single precision:

# In[ ]:


data = np.column_stack([Xs, Ts]).tolist()      # one scenario per data point: x, y, target

def loss_and_gradient(batch, w):
    # mean loss and gradient over the batch, by the MC++ DAG
    res = np.asarray(DAG.veval(sgLoss, deps, X + [T], batch, W, w))
    return res[:, 0].mean(), res[:, 1:].mean(axis=0)

L0, G0 = loss_and_gradient(data[:32], list(REG.weights))

net.zero_grad()
xt, tt = torch.tensor(Xs[:32], dtype=torch.float32), torch.tensor(Ts[:32], dtype=torch.float32)
lt = ((net(xt)[:, 0] - tt)**2).mean()
lt.backward()
Gt = np.concatenate([np.hstack([m.bias.grad.numpy()[:, None], m.weight.grad.numpy()]).ravel()
                     for m in net if isinstance(m, torch.nn.Linear)])
print(f"loss: MC++ {L0:.8f}, PyTorch {lt.item():.8f}")
print(f"gradient: largest relative difference {np.max(np.abs(G0 - Gt)) / np.max(np.abs(Gt)):.1e}")


# ## 3. Training with Adam, the gradient from the MC++ DAG
# 
# The weights are a PyTorch tensor, so that `torch.optim.Adam` updates them. At each step, the loss and its
# gradient over a mini-batch of 32 points come from the MC++ DAG; the gradient is handed to the optimizer through
# `w.grad`. The learning rate decays from $10^{-2}$ to $10^{-4}$ over 250 epochs.

# In[ ]:


w = torch.tensor(REG.weights, dtype=torch.float64, requires_grad=True)
nepoch, nbatch = 250, 32
opt = torch.optim.Adam([w], lr=1e-2)
sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=nepoch * len(data) // nbatch, eta_min=1e-4)

rng = np.random.default_rng(0)
history = []
start = time.time()
for epoch in range(nepoch):
    perm, eloss = rng.permutation(len(data)), 0.
    for b in range(0, len(data), nbatch):
        batch = [data[i] for i in perm[b:b + nbatch]]
        L, G = loss_and_gradient(batch, w.detach().tolist())
        opt.zero_grad()
        w.grad = torch.from_numpy(G)
        opt.step()
        sched.step()
        eloss += L * len(batch)
    history.append(eloss / len(data))
    if epoch % 50 == 0 or epoch == nepoch - 1:
        print(f"epoch {epoch:4d}   mean squared error (scaled) {history[-1]:.3e}   {time.time() - start:6.1f} s")

fig, ax = plt.subplots(figsize=(5, 3))
ax.semilogy(history)
ax.set_xlabel("epoch"); ax.set_ylabel("mean squared error (scaled)")
plt.show()


# ## 4. An `MLP` from the trained weights
# 
# `MLPREG.data(w)` gives the layers of the network for the weights `w`, in the form that `MLP.set_data` takes. The
# `MLP` holds them as data. At any input, it evaluates as the `MLPREG` with the trained weights:

# In[ ]:


wtrained = w.detach().tolist()
NN = mc.MLP()
NN.set_data(REG.data(wtrained))
print("MLP:", NN.nin, "inputs,", NN.nhid, "hidden layers,", NN.nout, "output")

xs = [0.25, -0.5]
D1 = mc.FFGraph(); X1 = D1.add_vars(2)
y_mlp = D1.eval(mc.FFMLP()(X1, NN), X1, xs)[0]
y_reg = REG.eval(xs + wtrained)[0]
print(f"at (x,y) = {xs}: MLP {y_mlp:.15f}, MLPREG {y_reg:.15f}")


# ## 5. The trained `MLP` in a DAG
# 
# The `FFMLP` operation embeds the network in a DAG, with the scaling of the inputs and the output written as
# ordinary DAG expressions: $\hat f(x,y) = \bar f + \sigma\, \mathrm{MLP}(x/3, y/3)$. The DAG is evaluated on a
# 121 x 121 grid in one `veval` call, and compared with `peaks`:

# In[ ]:


DAG2 = mc.FFGraph()
Z = DAG2.add_vars(2, "Z")                    # x, y in [-3,3]^2
OpNN = mc.FFMLP()
Fhat = fmean + fstd * OpNN([Z[0] / xscale, Z[1] / xscale], NN)[0]

grid = np.column_stack([GX.ravel(), GY.ravel()]).tolist()
FH = np.asarray(DAG2.veval([Fhat], Z, grid))[:, 0].reshape(GX.shape)
err = FH - FG
print(f"on the grid: RMS error {np.sqrt(np.mean(err**2)):.3f}, largest error {np.max(np.abs(err)):.3f},"
      f" for f in [{FG.min():.2f}, {FG.max():.2f}]")

fig, axs = plt.subplots(1, 3, figsize=(14, 4))
emax = np.max(np.abs(err))
for ax, Z_, title in zip(axs, [FG, FH, err], ["peaks", "trained MLP in the MC++ DAG", "error"]):
    if title == "error":
        cs = ax.contourf(GX, GY, Z_, levels=np.linspace(-emax, emax, 31), cmap="RdBu_r")
    else:
        cs = ax.contourf(GX, GY, Z_, levels=30, cmap="viridis")
    fig.colorbar(cs, ax=ax); ax.set_title(title); ax.set_xlabel("x"); ax.set_ylabel("y")
axs[2].plot(P[:, 0], P[:, 1], "k.", ms=2)
plt.show()


# The network is now an operation like any other in the DAG. It can be differentiated symbolically -- here the
# gradient of $\hat f$, against central finite differences -- and bounded or relaxed over a box: interval bounds and
# McCormick relaxations, at the centre of a box and over the whole domain, against the range of $\hat f$ sampled on
# the grid.

# In[ ]:


# gradient of the trained model, by forward symbolic differentiation of the DAG
_, _, dF = DAG2.fdiff([Fhat], Z)
z0, h = [0.5, -1.0], 1e-6
grad = DAG2.eval(list(dF), Z, z0)
fd = [(DAG2.eval([Fhat], Z, [z0[0] + h, z0[1]])[0] - DAG2.eval([Fhat], Z, [z0[0] - h, z0[1]])[0]) / (2 * h),
      (DAG2.eval([Fhat], Z, [z0[0], z0[1] + h])[0] - DAG2.eval([Fhat], Z, [z0[0], z0[1] - h])[0]) / (2 * h)]
print("gradient at", z0, ":", np.round(grad, 6), " finite differences:", np.round(fd, 6))

# bounds and relaxations over boxes
for box in ([(lo, up), (lo, up)], [(0.5, 1.0), (1.5, 2.0)]):
    IZ = [mc.Interval(a, b) for a, b in box]
    MZ = [mc.McCormick(IZ[i], 0.5 * sum(box[i])).sub(2, i) for i in range(2)]
    IF = DAG2.eval([Fhat], Z, IZ)[0]
    MF = DAG2.eval([Fhat], Z, MZ)[0]
    inbox = (GX >= box[0][0]) & (GX <= box[0][1]) & (GY >= box[1][0]) & (GY <= box[1][1])
    print(f"box {box}: sampled range [{FH[inbox].min():.3f}, {FH[inbox].max():.3f}]")
    print(f"   interval bounds [{IF.l:.3f}, {IF.u:.3f}];"
          f" McCormick relaxations at the centre: convex {MF.cv:.3f} <= {DAG2.eval([Fhat], Z, [0.5 * sum(b) for b in box])[0]:.3f} <= concave {MF.cc:.3f}")

