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
    model.eval()

    with torch.no_grad():
        for i in range(n_iterations):
            theta = model.sample_q_theta_given_eta(eta)
            theta_samples_list.append(theta)
            if (i + 1) % 100 == 0:
                print(f"Generated {i + 1}/{n_iterations} samples")

    return torch.stack(theta_samples_list, dim=0)  # [n_iterations, N_MC, p_theta]
