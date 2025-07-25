import torch

def get_marginal_samples(model, eta, n_iterations=1000):
    """
    Generate repeated samples from marginal q(theta) ≈ E_{p(eta)}[q(theta|eta)]
    Args:
        model: CutBayesFlow instance
        eta: Tensor of eta samples [N_MC, p_eta]
        n_iterations: Number of repeated samples
    Returns:
        Tensor of shape [n_iterations, N_MC, p_theta]
    """
    theta_samples_list = []
    with torch.no_grad():
        for _ in range(n_iterations):
            theta = model.sample_q_theta_given_eta(eta)
            theta_samples_list.append(theta)
    
    return torch.stack(theta_samples_list, dim=0)  # [n_iterations, N_MC, p_theta]