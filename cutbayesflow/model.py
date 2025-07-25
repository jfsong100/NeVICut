import torch
from torch import nn
from nflows.flows import Flow
from nflows.distributions.normal import StandardNormal
from nflows.transforms.base import CompositeTransform
from nflows.transforms.permutations import ReversePermutation
from nflows.transforms.autoregressive import MaskedPiecewiseRationalQuadraticAutoregressiveTransform


class CutBayesFlow(nn.Module):
    def __init__(self, prior_log_prob_fn, likelihood_log_prob_fn, eta_dim, theta_dim,
                 num_layers=8, hidden_features=32, num_bins=16, tail_bound=10.0,
                 use_residual_blocks=True, num_blocks=2, tails="linear",
                 map_to_simplex=False):
        super().__init__()
        self.p_eta = eta_dim
        self.p_theta = theta_dim
        self.map_to_simplex = map_to_simplex
        self.prior_log_prob_fn = prior_log_prob_fn
        self.likelihood_log_prob_fn = likelihood_log_prob_fn
        self.flow_dim = theta_dim - 1 if map_to_simplex else theta_dim

        self.flow = self._build_flow(
            self.flow_dim, num_layers, hidden_features, num_bins, tail_bound,
            use_residual_blocks, num_blocks, tails
        )

    def _build_flow(self, flow_dim, num_layers, hidden_features, num_bins, tail_bound,
                    use_residual_blocks, num_blocks, tails):
        transforms = []
        for _ in range(num_layers):
            transforms.append(ReversePermutation(features=flow_dim))
            transforms.append(
                MaskedPiecewiseRationalQuadraticAutoregressiveTransform(
                    features=flow_dim,
                    hidden_features=hidden_features,
                    num_bins=num_bins,
                    tails=tails,
                    tail_bound=tail_bound,
                    context_features=self.p_eta,
                    use_residual_blocks=use_residual_blocks,
                    num_blocks=num_blocks,
                )
            )
        return Flow(CompositeTransform(transforms), StandardNormal([flow_dim]))

    def stick_breaking_transform(self, x):
        """
        Convert [N, K-1] unconstrained x to [N, K] simplex via stick-breaking.
        """
        eps = 1e-8
        v = torch.sigmoid(x).clamp(eps, 1 - eps)
        
        # Use log-space for numerical stability
        log_one_minus_v = torch.log(1 - v + eps)
        log_cumprod = torch.cumsum(log_one_minus_v, dim=1)
        
        # Compute simplex coordinates
        theta = torch.zeros(x.shape[0], self.p_theta, device=x.device, dtype=x.dtype)
        theta[:, 0] = v[:, 0]
        for k in range(1, self.p_theta - 1):
            theta[:, k] = v[:, k] * torch.exp(log_cumprod[:, k - 1])
        theta[:, -1] = torch.exp(log_cumprod[:, -1])

        # Correct log determinant calculation
        log_det = torch.sum(torch.log(v + eps) + log_one_minus_v, dim=1)  # sigmoid derivatives
        
        # Add stick-breaking dependency terms: Σ_{k=1}^{K-1} (K-1-k) * log(1-v_k)
        weights = torch.arange(self.p_theta - 2, -1, -1, device=x.device, dtype=x.dtype)  # [K-2, K-3, ..., 0]
        log_det += torch.sum(weights.unsqueeze(0) * log_one_minus_v, dim=1)
        
        return theta, log_det

    def sample_q_theta_given_eta(self, eta):
        N = eta.shape[0]
        z = torch.randn((N, self.flow_dim), device=eta.device)
        x, _ = self.flow._transform(z, context=eta)  # Fixed: removed .forward
        if self.map_to_simplex:
            theta, _ = self.stick_breaking_transform(x)
        else:
            theta = x
        return theta

    def compute_loss(self, eta_batch, data_D2):
        if eta_batch.dim() != 2 or eta_batch.shape[1] != self.p_eta:
            raise ValueError(f"eta_batch shape must be [N, p_eta], got {eta_batch.shape}")

        N = eta_batch.shape[0]
        z = torch.randn((N, self.flow_dim), device=eta_batch.device)
        x, flow_log_det = self.flow._transform(z, context=eta_batch)  # Fixed: removed .forward

        if self.map_to_simplex:
            theta, stick_log_det = self.stick_breaking_transform(x)
            total_log_det = flow_log_det + stick_log_det
        else:
            theta = x
            total_log_det = flow_log_det

        log_q_z = self.flow._distribution.log_prob(z).sum(dim=-1)
        log_q_theta = log_q_z - total_log_det

        log_p_theta = self.prior_log_prob_fn(theta)
        log_p_likelihood = self.likelihood_log_prob_fn(theta, eta_batch, data_D2)

        kl = log_q_theta - (log_p_theta + log_p_likelihood)
        return kl.mean()