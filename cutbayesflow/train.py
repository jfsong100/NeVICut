import torch
import numpy as np
import torch.nn as nn
import torch.optim as optim

def train_cut_bayes_flow(
    model,
    eta_samples,
    data_D2,
    epochs=3000,
    batch_size=None,
    lr=1e-3,
    patience=100,
    log_interval=100,
    seed=None,
    verbose=True,
    early_stopping=False,
    es_patience=None
):
    """
    Train CutBayesFlow model with optional early stopping.
    Args:
        model: CutBayesFlow instance
        eta_samples: Tensor of eta samples [num_samples, eta_dim]
        data_D2: Data used for likelihood evaluation
        epochs: Number of training epochs
        batch_size: Batch size (default full-batch if None)
        lr: Learning rate
        patience: Patience for LR scheduler
        log_interval: Steps between progress logs
        seed: Random seed
        verbose: Print progress if True
        early_stopping: Whether to stop early when loss plateaus
        es_patience: Number of no‐improve steps before stopping (defaults to 2*patience)
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
        print(f"Using batch_size={batch_size} (full-batch: {batch_size == len(eta_samples)})")

    optimizer = optim.Adam(model.parameters(), lr=lr)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min',
        factor=0.5, patience=patience,
        min_lr=1e-6, verbose=verbose
    )

    # Early‐stopping setup
    if early_stopping:
        best_loss = float('inf')
        no_improve = 0
        if es_patience is None:
            es_patience = patience * 2

    loss_history = []
    for step in range(epochs):
        # --- batch sampling ---
        if batch_size == len(eta_samples):
            eta_batch = eta_samples
        else:
            idx = torch.randint(0, len(eta_samples), (batch_size,))
            eta_batch = eta_samples[idx]

        # --- forward & loss ---
        loss = model.compute_loss(eta_batch, data_D2)
        if torch.isnan(loss):
            print(f"Warning: NaN loss at step {step}, stopping training")
            break

        # --- backward & step ---
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
        optimizer.step()

        # --- scheduler ---
        scheduler.step(loss)

        loss_val = loss.item()
        loss_history.append(loss_val)

        # --- early stopping check ---
        if early_stopping:
            if loss_val < best_loss - 1e-8:
                best_loss = loss_val
                no_improve = 0
            else:
                no_improve += 1
                if no_improve >= es_patience:
                    if verbose:
                        print(f"Early stopping at step {step}: no improvement for {no_improve} steps.")
                    break

        # --- logging ---

        if verbose and (step % log_interval == 0 or step == epochs - 1):
            current_lr = optimizer.param_groups[0]['lr']
            print(f"Step {step:5d} | Loss: {loss_val:.6f} | LR: {current_lr:.2e}")

    return loss_history
