import torch
import numpy as np
import torch.nn as nn
import torch.optim as optim

def train_cut_bayes_flow(model, eta_samples, data_D2, epochs=3000, batch_size=None, lr=1e-3,
                         patience = 100, log_interval=100, seed = None, verbose=True):
    """
    Train CutBayesFlow model.

    Args:
        model: CutBayesFlow instance
        eta_samples: Tensor of η samples [num_samples, eta_dim]
        data_D2: Data used for likelihood evaluation
        epochs: Number of training epochs
        batch_size: Batch size (default full-batch if None)
        lr: Learning rate
        log_interval: Steps between progress logs
        verbose: Print progress if True
    Returns:
        List of loss values over training
    """

    if seed is not None:
        torch.manual_seed(seed)
        np.random.seed(seed)
    
    # Default to full batch
    if batch_size is None:
        batch_size = len(eta_samples)
        if verbose:
            print(f"⚡ Using full-batch mode with batch_size={batch_size}")

    optimizer = optim.Adam(model.parameters(), lr=lr)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=patience)
    loss_history = []

    for step in range(epochs):
        if batch_size == len(eta_samples):
            eta_batch = eta_samples  # Full batch
        else:
            idx = torch.randint(0, len(eta_samples), (batch_size,))
            eta_batch = eta_samples[idx]

        # Compute loss
        loss = model.compute_loss(eta_batch, data_D2)

        # Backpropagation
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()
        scheduler.step(loss)

        loss_history.append(loss.item())

        # Logging
        if verbose and (step % log_interval == 0 or step == epochs - 1):
            print(f"Step {step:5d} | Loss: {loss.item():.6f}")

    return loss_history
