import torch
from torch import nn
from nflows.transforms.base import Transform


# ─────────────────────────────────────────────
# UMNN pieces
# ─────────────────────────────────────────────
class _ANet1D(nn.Module):
    def __init__(self, eta_dim, hidden=128):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(eta_dim, hidden), nn.LeakyReLU(0.01),
            nn.Linear(hidden, hidden), nn.LeakyReLU(0.01),
            nn.Linear(hidden, 1)
        )
        for m in self.net:
            if isinstance(m, nn.Linear):
                nn.init.normal_(m.weight, 0, 0.01)
                nn.init.constant_(m.bias, 0)

    def forward(self, eta):          
        return self.net(eta).squeeze(-1)

class _HNet1D(nn.Module):
    def __init__(self, eta_dim, hidden=128):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(1 + eta_dim, hidden), nn.LeakyReLU(0.01),
            nn.Linear(hidden, hidden),     nn.LeakyReLU(0.01),
            nn.Linear(hidden, hidden),     nn.LeakyReLU(0.01),
            nn.Linear(hidden, 1)
        )

    def forward(self, t, eta):       
        B, N = t.shape
        inp = torch.cat([t.unsqueeze(-1),
                         eta.unsqueeze(1).expand(B, N, -1)], dim=-1)
        return torch.exp(self.net(inp).squeeze(-1)) 


class UMNN1DTransform(Transform):
    def __init__(self, eta_dim, t0=0.0, steps=500, hidden=128):
        super().__init__()
        self.t0 = t0
        self.register_buffer('u', torch.linspace(0., 1., steps))
        self.a = _ANet1D(eta_dim, hidden)
        self.h = _HNet1D(eta_dim, hidden)

    def forward(self, inputs, context=None):
        z = inputs
        assert z.dim() == 2 and z.size(1) == 1, "UMNN1DTransform expects [B,1] inputs"
        if context is None:
            raise ValueError("UMNN1DTransform requires context η")

        B = z.size(0)
        z0 = z[:, 0]                   
        t0 = self.t0
        steps = self.u.numel()
        du = 1.0 / (steps - 1)

        tg = t0 + (z0 - t0).unsqueeze(1) * self.u.unsqueeze(0) 
        hv = self.h(tg, context)                               

        trap = 0.5 * (hv[:, 1:] + hv[:, :-1]) * du             
        integral = trap.sum(1) * (z0 - t0)                   
        a = self.a(context)                                     

        theta = a + integral                                   
        hz = self.h(z0.unsqueeze(1), context).squeeze(-1)       
        logabsdet = torch.log(hz.clamp(min=1e-12))              

        return theta.unsqueeze(1), logabsdet

    def _call(self, inputs, context=None):
        return self.forward(inputs, context)

    def inverse(self, inputs, context=None):
        raise NotImplementedError("UMNN1DTransform.inverse is not implemented.")
