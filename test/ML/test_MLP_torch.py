# Generates the TorchScript models read by test_MLP_torch.cpp:
#   python3 test_MLP_torch.py
# Models that mc::MLP::read_data reads, and models that it must refuse.
import torch, torch.nn as nn

torch.manual_seed(1)


class Net(nn.Module):
    """The form read_data expects: a Sequential in the attribute 'net'."""

    def __init__(self, seq):
        super().__init__()
        self.net = seq

    def forward(self, x):
        return self.net(x)


class Scaled(Net):
    """forward does more than go through 'net'."""

    def forward(self, x):
        return 2.0 * self.net(x)


class NoBatch(Net):
    """forward takes one input only, not a batch."""

    def forward(self, x):
        return self.net(x.view(-1))


class NoNet(nn.Module):
    """No attribute 'net'."""

    def __init__(self, seq):
        super().__init__()
        self.seq = seq

    def forward(self, x):
        return self.seq(x)


class MyLinear(nn.Module):
    """A layer of another class, which registers its bias before its weight."""

    def __init__(self, nin, nout):
        super().__init__()
        self.bias = nn.Parameter(torch.randn(nout) * 0.1)
        self.weight = nn.Parameter(torch.randn(nout, nin) * 0.5)

    def forward(self, x):
        return torch.nn.functional.linear(x, self.weight, self.bias)


def save(name, module, double=False, train=False):
    if double:
        module = module.double()
    module.train() if train else module.eval()
    torch.jit.script(module).save("torch_" + name + ".pt")


def mlp(activ, nin=2, nhid=5, nout=2):
    return nn.Sequential(nn.Linear(nin, nhid), activ, nn.Linear(nhid, nout))


# read
save("tanh", Net(nn.Sequential(nn.Linear(2, 5), nn.Tanh(), nn.Linear(5, 4), nn.Tanh(), nn.Linear(4, 2))))
save("tanh_double", Net(nn.Sequential(nn.Linear(2, 5), nn.Tanh(), nn.Linear(5, 4), nn.Sigmoid(), nn.Linear(4, 2))), double=True)
save("two_linear", Net(nn.Sequential(nn.Linear(2, 5), nn.Linear(5, 4), nn.Tanh(), nn.Linear(4, 2))))
save("out_activ", Net(nn.Sequential(nn.Linear(2, 5), nn.ReLU(), nn.Linear(5, 3), nn.Sigmoid())))
save("nobias", Net(nn.Sequential(nn.Linear(3, 5, bias=False), nn.Tanh(), nn.Linear(5, 2))))
save("dropout", Net(nn.Sequential(nn.Linear(2, 6), nn.Tanh(), nn.Dropout(0.5), nn.Linear(6, 2))), train=True)
save("nobatch", NoBatch(nn.Sequential(nn.Linear(3, 5), nn.ReLU(), nn.Linear(5, 4), nn.Tanh(), nn.Linear(4, 2))))
# refused
save("leakyrelu", Net(mlp(nn.LeakyReLU(0.3))))
save("tanhshrink", Net(mlp(nn.Tanhshrink())))
save("softplus", Net(mlp(nn.Softplus())))
save("scaled", Scaled(mlp(nn.Tanh())))
save("nonet", NoNet(mlp(nn.Tanh())))
save("mylinear", Net(nn.Sequential(MyLinear(2, 5), nn.Tanh(), MyLinear(5, 2))))
