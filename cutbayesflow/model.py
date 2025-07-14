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
        """
        prior_log_prob_fn: callable, log p(θ)
        likelihood_log_prob_fn: callable, log p(y|θ,η)
        map_to_simplex: bool, whether to apply softmax to θ before prior/likelihood
        """
        super().__init__()
        self.p_eta = eta_dim
        self.p_theta = theta_dim
        self.prior_log_prob_fn = prior_log_prob_fn
        self.likelihood_log_prob_fn = likelihood_log_prob_fn
        self.num_layers = num_layers
        self.hidden_features = hidden_features
        self.num_bins = num_bins
        self.tail_bound = tail_bound
        self.use_residual_blocks = use_residual_blocks
        self.num_blocks = num_blocks
        self.tails = tails
        self.map_to_simplex = map_to_simplex  # NEW

        self.flow = self._build_flow()

    def _build_flow(self):
        transforms = []
        for _ in range(self.num_layers):
            transforms.append(ReversePermutation(features=self.p_theta))
            transforms.append(
                MaskedPiecewiseRationalQuadraticAutoregressiveTransform(
                    features=self.p_theta,
                    hidden_features=self.hidden_features,
                    num_bins=self.num_bins,
                    tails=self.tails,
                    tail_bound=self.tail_bound,
                    context_features=self.p_eta,
                    use_residual_blocks=self.use_residual_blocks,
                    num_blocks=self.num_blocks,
                )
            )
        return Flow(CompositeTransform(transforms), StandardNormal([self.p_theta]))

    def sample_q_theta_given_eta(self, eta):
        """Sample θ ~ q(θ|η) for each η"""
        N = eta.shape[0]
        z = torch.randn((N, self.p_theta))
        theta, _ = self.flow._transform(z, context=eta)
        if self.map_to_simplex:
            theta = torch.nn.functional.softmax(theta, dim=1)
        return theta  # shape: [N, p_theta]

    def compute_loss(self, eta_batch, data_D2):
        """
        Approximate E_{p(η)} [KL(q(θ|η)||p(y|θ,η)p(θ))]
        """
        if eta_batch.dim() != 2 or eta_batch.shape[1] != self.p_eta:
            raise ValueError(f"eta_batch shape must be [N, p_eta], got {eta_batch.shape}")

        N = eta_batch.shape[0]
        z = torch.randn((N, self.p_theta))
        theta, log_det = self.flow._transform(z, context=eta_batch)
        if self.map_to_simplex:
            theta = torch.nn.functional.softmax(theta, dim=1)
        log_q = self.flow._distribution.log_prob(z).sum(dim=-1) - log_det
        log_p_theta = self.prior_log_prob_fn(theta)
        log_p_likelihood = self.likelihood_log_prob_fn(theta, eta_batch, data_D2)
        kl = log_q - (log_p_theta + log_p_likelihood)
        return kl.mean()
