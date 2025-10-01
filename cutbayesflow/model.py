import torch
from torch import nn
from nflows.flows import Flow
from nflows.distributions.normal import StandardNormal
from nflows.transforms.base import CompositeTransform
from nflows.transforms.permutations import ReversePermutation, RandomPermutation
from nflows.transforms.autoregressive import MaskedPiecewiseRationalQuadraticAutoregressiveTransform
from nflows.transforms.coupling import PiecewiseRationalQuadraticCouplingTransform, AffineCouplingTransform
from nflows.nn.nets import ResidualNet


def generate_cyclic_masks(dim, num_layers):
    masks = []
    for i in range(num_layers):
        mask = torch.zeros(dim, dtype=torch.bool)
        mask[i % dim] = True
        masks.append(mask)
    return masks


class CutBayesFlow(nn.Module):
    def __init__(
        self,
        prior_log_prob_fn,
        likelihood_log_prob_fn,
        eta_dim,
        theta_dim,
        num_layers=8,
        hidden_features=64,
        num_bins=16,
        tail_bound=10.0,
        use_residual_blocks=True,
        num_blocks=2,
        tails="linear",
        map_to_simplex=False,
        flow_type="NSF-AR",
        use_alr=False,
        activation=nn.ReLU()
    ):
        super().__init__()
        self.p_eta = eta_dim
        self.p_theta = theta_dim
        self.map_to_simplex = map_to_simplex
        self.use_alr = use_alr
        self.flow_type = flow_type.upper()
        self.prior_log_prob_fn = prior_log_prob_fn
        self.likelihood_log_prob_fn = likelihood_log_prob_fn
        self.activation = activation

        if (map_to_simplex or use_alr) and theta_dim < 2:
            raise ValueError("Simplex or ALR transform requires theta_dim ≥ 2.")

        # Flow dimension is K-1 when mapping to simplex or using ALR; else K
        self.flow_dim = theta_dim - 1 if (map_to_simplex or use_alr) else theta_dim
        self.flow = self._build_flow(
            num_layers,
            hidden_features,
            num_bins,
            tail_bound,
            use_residual_blocks,
            num_blocks,
            tails
        )

    def _build_flow(
        self,
        num_layers,
        hidden_features,
        num_bins,
        tail_bound,
        use_residual_blocks,
        num_blocks,
        tails
    ):
        transforms = []
        # Use cyclic masks for NSF and RealNVP coupling
        if self.flow_type in ("NSF-C", "REALNVP"):
            masks = generate_cyclic_masks(self.flow_dim, num_layers)

        for i in range(num_layers):
            if self.flow_type == "NSF-AR":
                transforms.append(ReversePermutation(features=self.flow_dim))
                transforms.append(
                    MaskedPiecewiseRationalQuadraticAutoregressiveTransform(
                        features=self.flow_dim,
                        hidden_features=hidden_features,
                        num_bins=num_bins,
                        tails=tails,
                        tail_bound=tail_bound,
                        context_features=self.p_eta,
                        use_residual_blocks=use_residual_blocks,
                        num_blocks=num_blocks,
                        min_bin_width=1e-8,
                        min_bin_height=1e-8,
                        min_derivative=1e-8,
                    )
                )

            elif self.flow_type == "NSF-C":
                mask = masks[i]
                def make_net(in_f, out_f):
                    return ResidualNet(
                        in_features=in_f,
                        out_features=out_f,
                        hidden_features=hidden_features,
                        context_features=self.p_eta,
                        num_blocks=num_blocks,
                        activation=self.activation,
                    )
                transforms.append(
                    PiecewiseRationalQuadraticCouplingTransform(
                        mask=mask,
                        transform_net_create_fn=make_net,
                        num_bins=num_bins,
                        tails=tails,
                        tail_bound=tail_bound,
                    )
                )
                transforms.extend([
                    ReversePermutation(features=self.flow_dim),
                    RandomPermutation(features=self.flow_dim),
                ])

            elif self.flow_type == "REALNVP":
                mask = masks[i]
                def make_net(in_f, out_f):
                    return ResidualNet(
                        in_features=in_f,
                        out_features=out_f,
                        hidden_features=hidden_features,
                        context_features=self.p_eta,
                        num_blocks=num_blocks,
                        activation=self.activation,
                    )
                transforms.append(
                    AffineCouplingTransform(
                        mask=mask,
                        transform_net_create_fn=make_net
                    )
                )
                transforms.extend([
                    ReversePermutation(features=self.flow_dim),
                    RandomPermutation(features=self.flow_dim),
                ])

            else:
                raise NotImplementedError(f"Unsupported flow_type '{self.flow_type}'")

        return Flow(CompositeTransform(transforms), StandardNormal([self.flow_dim]))

    def alr_inverse(self, z):
        exp_z = torch.exp(z)
        denom = 1.0 + exp_z.sum(dim=-1, keepdim=True)
        return torch.cat([exp_z / denom, 1.0 / denom], dim=-1)

    def alr_log_det_jacobian(self, z):
        eps = 1e-10
        theta = self.alr_inverse(z).clamp(min=eps, max=1 - eps)
        return torch.log(theta).sum(dim=-1)

    def stick_breaking_transform(self, x):
        """More numerically stable stick-breaking with log-space computation"""
        # Compute v in log-space for numerical stability
        log_v = -torch.nn.functional.softplus(-x)  # log(sigmoid(x))
        log_1_v = -torch.nn.functional.softplus(x)  # log(1-sigmoid(x))
        
        batch, d = x.shape
        device = x.device
        
        # Compute stick-breaking weights in log-space
        # log(prefixes) = cumsum(log(1-v_{1..k-1}))
        log_prefixes = torch.cat([
            torch.zeros(batch, 1, device=device),
            torch.cumsum(log_1_v[:, :-1], dim=1)
        ], dim=1)
        
        # theta_k = v_k * prod_{j<k} (1-v_j) for k <= d
        log_theta = log_v + log_prefixes
        
        # Last component is prod_{j=1..d} (1-v_j)
        log_last = torch.sum(log_1_v, dim=1, keepdim=True)
        
        # Combine all components
        log_theta = torch.cat([log_theta, log_last], dim=1)
        theta = torch.exp(log_theta)
        
        # Compute log det Jacobian
        weights = torch.arange(d, 0, -1, device=device, dtype=torch.float)
        log_det = torch.sum(log_v, dim=1) + torch.sum(weights * log_1_v, dim=1)
        
        return theta, log_det

    def sample_q_theta_given_eta(self, eta):
        N = eta.shape[0]
        z = torch.randn(N, self.flow_dim, device=eta.device)
        x, _ = self.flow._transform(z, context=eta)

        if self.use_alr:
            return self.alr_inverse(x)
        if self.map_to_simplex:
            theta, _ = self.stick_breaking_transform(x)
            return theta
        return x

    def compute_loss(self, eta_batch, data_D2):
        N = eta_batch.shape[0]
        z = torch.randn(N, self.flow_dim, device=eta_batch.device)
        x, flow_log_det = self.flow._transform(z, context=eta_batch)

        if self.use_alr:
            theta = self.alr_inverse(x)
            sb_log_det = self.alr_log_det_jacobian(x)
            total_log_det = flow_log_det + sb_log_det
        elif self.map_to_simplex:
            theta, sb_log_det = self.stick_breaking_transform(x)
            total_log_det = flow_log_det + sb_log_det
        else:
            theta = x
            total_log_det = flow_log_det

        log_qz = self.flow._distribution.log_prob(z)
        log_q = log_qz - total_log_det
        log_p = self.prior_log_prob_fn(theta)
        log_lik = self.likelihood_log_prob_fn(theta, eta_batch, data_D2)

        return (log_q - (log_p + log_lik)).mean()
