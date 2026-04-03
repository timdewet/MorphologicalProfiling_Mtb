#!/usr/bin/env python3
"""
Phase 2: Supervised Contrastive Learning + Optimal Transport Matching
=====================================================================
Trains a ResNet-18 backbone with SupCon loss on single-cell phase-contrast
crops, extracts 512-d embeddings, computes pairwise Wasserstein distances
between conditions, and runs permutation tests for significance.

Usage:
    python supcon_ot_pipeline.py                # full pipeline
    python supcon_ot_pipeline.py --train-only   # just train
    python supcon_ot_pipeline.py --embed-only   # extract embeddings (needs checkpoint)
    python supcon_ot_pipeline.py --ot-only      # OT + permutation (needs embeddings)

Requirements:
    pip install torch torchvision h5py albumentations
    pip install pytorch-metric-learning POT umap-learn
    pip install matplotlib seaborn pandas scikit-learn scipy
"""

import argparse
import time
from collections import Counter
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.spatial.distance import squareform
from scipy.stats import false_discovery_control
from torch.utils.data import DataLoader, Dataset, WeightedRandomSampler

from config_supcon import CONFIG


# ==============================================================================
# DATASET
# ==============================================================================

class CropHDF5Dataset(Dataset):
    """
    Lazy-loading dataset backed by an HDF5 file of cell crops.

    Each item returns (crop_tensor, label_index) where crop_tensor is
    (2, 96, 96) float32 and label_index is an integer.

    Crops are loaded into memory at init for fast random access (the HDF5
    gzip compression makes lazy random reads extremely slow).
    """

    def __init__(self, h5_path, indices=None, label_to_idx=None,
                 transform=None):
        """
        Args:
            h5_path:      path to the HDF5 file from extract_crops.py
            indices:       subset of row indices to use (None = all)
            label_to_idx:  dict mapping condition_label → int (built if None)
            transform:     albumentations transform (or None)
        """
        self.transform = transform

        # Load everything into memory
        with h5py.File(str(h5_path), "r") as f:
            all_labels = f["condition_labels"][:].astype(str)
            all_is_control = f["is_control"][:]

            if indices is not None:
                self.indices = np.array(indices)
            else:
                self.indices = np.arange(len(all_labels))

            # Load only the crops we need
            print(f"  Loading {len(self.indices):,} crops into memory ...")
            self.crops = f["crops"][self.indices]  # (N, 2, 96, 96)

        self.labels = all_labels[self.indices]
        self.is_control = all_is_control[self.indices]

        # Build or reuse label mapping
        if label_to_idx is not None:
            self.label_to_idx = label_to_idx
        else:
            unique = sorted(set(self.labels))
            self.label_to_idx = {lbl: i for i, lbl in enumerate(unique)}

        self.targets = np.array([self.label_to_idx.get(l, -1)
                                 for l in self.labels])

    def __len__(self):
        return len(self.indices)

    def __getitem__(self, idx):
        crop = self.crops[idx].copy()  # (2, 96, 96) float32

        if self.transform is not None:
            # Two-view SupCon: generate two differently-augmented views
            img = crop.transpose(1, 2, 0)  # (96, 96, 2)
            view1 = self.transform(image=img)["image"].transpose(2, 0, 1)
            view2 = self.transform(image=img)["image"].transpose(2, 0, 1)
            crop = np.stack([view1, view2], axis=0)  # (2, 2, 96, 96)
            return torch.from_numpy(crop), self.targets[idx]

        return torch.from_numpy(crop), self.targets[idx]


# ── Data splitting & augmentation ────────────────────────────────────────────

def build_fov_split(h5_path, val_fraction, seed):
    """
    FOV-level train/val split to prevent data leakage.

    Returns (train_indices, val_indices) as numpy int arrays.
    """
    rng = np.random.default_rng(seed)

    with h5py.File(str(h5_path), "r") as f:
        fov_ids = f["fov_ids"][:]
        is_control = f["is_control"][:]

    unique_fovs = np.unique(fov_ids)
    rng.shuffle(unique_fovs)
    n_val = max(1, int(len(unique_fovs) * val_fraction))
    val_fovs = set(unique_fovs[:n_val])

    all_idx = np.arange(len(fov_ids))
    val_mask = np.array([fid in val_fovs for fid in fov_ids])

    train_indices = all_idx[~val_mask]
    val_indices = all_idx[val_mask]

    print(f"FOV split: {len(unique_fovs)} FOVs → "
          f"{len(unique_fovs) - n_val} train, {n_val} val")
    print(f"  Train cells: {len(train_indices):,}  |  Val cells: {len(val_indices):,}")

    return train_indices, val_indices


def build_augmentation():
    """Build albumentations transforms for SupCon training."""
    import albumentations as A

    train_transform = A.Compose([
        A.HorizontalFlip(p=0.5),
        A.VerticalFlip(p=0.5),
        A.RandomRotate90(p=0.5),
        A.Affine(
            shift_limit=0.1, scale_limit=0.1, rotate_limit=30,
            border_mode=0, p=0.5,
        ),
        A.ElasticTransform(alpha=30, sigma=4, p=0.3),
        A.RandomBrightnessContrast(
            brightness_limit=0.15, contrast_limit=0.15, p=0.4,
        ),
        A.GaussNoise(std_range=(0.0, 0.08), p=0.3),
        A.GaussianBlur(blur_limit=(3, 5), p=0.2),
    ])

    val_transform = None  # no augmentation for validation

    return train_transform, val_transform


def build_sampler(dataset, max_per_class, exclude_controls=True):
    """
    Weighted random sampler with class balancing.

    Controls are excluded from training if exclude_controls=True.
    """
    targets = dataset.targets
    is_control = dataset.is_control

    # Compute weights
    weights = np.zeros(len(targets), dtype=np.float64)
    counts = Counter()

    for i, (t, ctrl) in enumerate(zip(targets, is_control)):
        if exclude_controls and ctrl:
            weights[i] = 0.0
        else:
            counts[t] += 1

    # Inverse frequency, capped
    for i, (t, ctrl) in enumerate(zip(targets, is_control)):
        if exclude_controls and ctrl:
            continue
        weights[i] = 1.0 / max(counts[t], 1)

    # Effective epoch size: sum of min(count, max_per_class) per class
    epoch_size = sum(min(c, max_per_class) for c in counts.values())

    sampler = WeightedRandomSampler(
        weights=weights,
        num_samples=epoch_size,
        replacement=True,
    )
    return sampler


# ==============================================================================
# MODEL
# ==============================================================================

class SupConResNet(nn.Module):
    """
    ResNet-18 encoder + MLP projection head for Supervised Contrastive Learning.

    The encoder outputs 512-d features (L2-normalized).
    The projector maps to 128-d projections (L2-normalized, used only for loss).
    """

    def __init__(self, in_channels=2, embedding_dim=512, projection_dim=128,
                 pretrained=True):
        super().__init__()
        from torchvision.models import resnet18, ResNet18_Weights

        weights = ResNet18_Weights.DEFAULT if pretrained else None
        backbone = resnet18(weights=weights)

        # Replace first conv: 3ch → in_channels
        old_conv = backbone.conv1
        new_conv = nn.Conv2d(in_channels, 64, kernel_size=7, stride=2,
                             padding=3, bias=False)
        with torch.no_grad():
            # Average RGB weights for phase channel
            new_conv.weight[:, 0] = old_conv.weight.mean(dim=1)
            if in_channels >= 2:
                # Green channel weights for mask channel
                new_conv.weight[:, 1] = old_conv.weight[:, 1]
        backbone.conv1 = new_conv

        # Remove final FC, keep avgpool → (batch, 512, 1, 1)
        self.encoder = nn.Sequential(*list(backbone.children())[:-1])

        # MLP projection head
        self.projector = nn.Sequential(
            nn.Linear(embedding_dim, embedding_dim),
            nn.ReLU(inplace=True),
            nn.Linear(embedding_dim, projection_dim),
        )

        self.embedding_dim = embedding_dim
        self.projection_dim = projection_dim

    def encode(self, x):
        """Extract L2-normalised 512-d embeddings (for inference)."""
        h = self.encoder(x)
        h = h.squeeze(-1).squeeze(-1)  # (batch, 512)
        return F.normalize(h, dim=1)

    def forward(self, x):
        """Returns (features_512d, projections_128d), both L2-normalised."""
        h = self.encoder(x)
        h = h.squeeze(-1).squeeze(-1)      # (batch, 512)
        features = F.normalize(h, dim=1)
        projections = self.projector(h)
        projections = F.normalize(projections, dim=1)
        return features, projections


# ==============================================================================
# TRAINING
# ==============================================================================

def get_device():
    """Auto-select best available device."""
    if torch.cuda.is_available():
        return torch.device("cuda")
    elif hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def train_supcon(config):
    """Train SupCon model and save best checkpoint."""
    from pytorch_metric_learning.losses import SupConLoss

    device = get_device()
    print(f"Device: {device}")

    h5_path = config["crop_output_h5"]
    seed = config["random_seed"]

    # Split
    train_idx, val_idx = build_fov_split(h5_path, config["val_fraction"], seed)

    # Build label mapping from non-control conditions only
    with h5py.File(str(h5_path), "r") as f:
        all_labels = f["condition_labels"][:].astype(str)
        all_is_control = f["is_control"][:]

    # Label mapping: only non-control conditions
    if config["exclude_controls_from_training"]:
        non_ctrl_labels = set(all_labels[~all_is_control])
    else:
        non_ctrl_labels = set(all_labels)
    label_to_idx = {lbl: i for i, lbl in enumerate(sorted(non_ctrl_labels))}
    print(f"Training classes: {len(label_to_idx)}")

    # Augmentation
    train_tf, val_tf = build_augmentation()

    # Datasets
    train_ds = CropHDF5Dataset(h5_path, indices=train_idx,
                               label_to_idx=label_to_idx, transform=train_tf)
    val_ds = CropHDF5Dataset(h5_path, indices=val_idx,
                             label_to_idx=label_to_idx, transform=val_tf)

    # Sampler
    train_sampler = build_sampler(
        train_ds,
        max_per_class=config["max_samples_per_class"],
        exclude_controls=config["exclude_controls_from_training"],
    )

    # DataLoaders
    train_loader = DataLoader(
        train_ds, batch_size=config["batch_size"], sampler=train_sampler,
        num_workers=0, drop_last=True,
    )

    # For validation: filter out controls
    if config["exclude_controls_from_training"]:
        val_valid = np.array([not c for c in val_ds.is_control])
        val_valid_idx = np.where(val_valid & (val_ds.targets >= 0))[0]
    else:
        val_valid_idx = np.arange(len(val_ds))

    val_loader = DataLoader(
        val_ds, batch_size=config["batch_size"], shuffle=False,
        num_workers=0,
    )

    # Model
    model = SupConResNet(
        in_channels=config["supcon_input_channels"],
        embedding_dim=config["embedding_dim"],
        projection_dim=config["projection_dim"],
        pretrained=config["pretrained"],
    ).to(device)

    # Freeze early ResNet layers to prevent overfitting on small dataset.
    # conv1, bn1, relu, maxpool = children 0-3; layer1 = child 4; layer2 = child 5
    freeze_layers = config.get("freeze_layers", 0)
    if freeze_layers > 0:
        children = list(model.encoder.children())
        # Freeze first N children (conv1/bn1/relu/maxpool count as separate)
        for child in children[:freeze_layers]:
            for param in child.parameters():
                param.requires_grad = False
        n_frozen = sum(1 for p in model.parameters() if not p.requires_grad)
        n_total = sum(1 for p in model.parameters())
        print(f"  Frozen {n_frozen}/{n_total} parameter groups "
              f"(first {freeze_layers} encoder children)")

    # Loss, optimizer, scheduler
    criterion = SupConLoss(temperature=config["temperature"])
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=config["lr"],
        weight_decay=config["weight_decay"],
    )
    warmup_epochs = config.get("warmup_epochs", 0)
    warmup_scheduler = torch.optim.lr_scheduler.LinearLR(
        optimizer, start_factor=0.01, end_factor=1.0,
        total_iters=warmup_epochs,
    ) if warmup_epochs > 0 else None
    cosine_scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=config["epochs"] - warmup_epochs,
    )
    scheduler = torch.optim.lr_scheduler.SequentialLR(
        optimizer,
        schedulers=[warmup_scheduler, cosine_scheduler],
        milestones=[warmup_epochs],
    ) if warmup_scheduler else cosine_scheduler

    # Checkpoint dir
    ckpt_dir = Path(config["checkpoint_dir"])
    ckpt_dir.mkdir(parents=True, exist_ok=True)

    # Training loop
    best_val_loss = float("inf")
    patience_counter = 0
    history = {"train_loss": [], "val_loss": [], "lr": []}

    print(f"\nStarting training for up to {config['epochs']} epochs ...")
    print(f"  Batch size: {config['batch_size']}")
    print(f"  Train batches/epoch: {len(train_loader)}")

    for epoch in range(config["epochs"]):
        t0 = time.time()

        # ── Train ────────────────────────────────────────────────────────
        model.train()
        train_losses = []

        for crops, targets in train_loader:
            # Skip batches with controls (target = -1)
            valid = targets >= 0
            if valid.sum() < 2:
                continue
            targets = targets[valid].to(device)

            # Two-view SupCon: crops is (B, 2, C, H, W)
            # Unpack views and concatenate along batch dim
            crops = crops[valid]
            view1 = crops[:, 0].contiguous().to(device)  # (B, C, H, W)
            view2 = crops[:, 1].contiguous().to(device)  # (B, C, H, W)
            all_views = torch.cat([view1, view2], dim=0)  # (2B, C, H, W)
            all_targets = torch.cat([targets, targets], dim=0)  # (2B,)

            _, projections = model(all_views)
            loss = criterion(projections, all_targets)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            train_losses.append(loss.item())

        avg_train = np.mean(train_losses) if train_losses else float("nan")

        # ── Validate ─────────────────────────────────────────────────────
        model.eval()
        val_losses = []

        with torch.no_grad():
            for crops, targets in val_loader:
                valid = targets >= 0
                if valid.sum() < 2:
                    continue
                crops = crops[valid].to(device)
                targets = targets[valid].to(device)

                _, projections = model(crops)
                loss = criterion(projections, targets)
                val_losses.append(loss.item())

        avg_val = np.mean(val_losses) if val_losses else float("nan")
        lr = scheduler.get_last_lr()[0]
        scheduler.step()

        history["train_loss"].append(avg_train)
        history["val_loss"].append(avg_val)
        history["lr"].append(lr)

        elapsed = time.time() - t0
        print(f"  Epoch {epoch + 1:3d}/{config['epochs']}  "
              f"train={avg_train:.4f}  val={avg_val:.4f}  "
              f"lr={lr:.2e}  ({elapsed:.1f}s)")

        # ── Checkpoint ───────────────────────────────────────────────────
        if avg_val < best_val_loss:
            best_val_loss = avg_val
            patience_counter = 0
            torch.save({
                "epoch": epoch + 1,
                "model_state_dict": model.state_dict(),
                "optimizer_state_dict": optimizer.state_dict(),
                "val_loss": avg_val,
                "label_to_idx": label_to_idx,
                "config": {k: v for k, v in config.items()
                           if isinstance(v, (int, float, str, bool, list))},
            }, ckpt_dir / "best_supcon.pth")
            print(f"    ✓ New best (val={avg_val:.4f})")
        else:
            patience_counter += 1
            if patience_counter >= config["patience"]:
                print(f"  Early stopping at epoch {epoch + 1} "
                      f"(patience={config['patience']})")
                break

    # ── Save training curves ─────────────────────────────────────────────
    fig_dir = Path(config["fig_dir"])
    fig_dir.mkdir(parents=True, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    epochs_x = range(1, len(history["train_loss"]) + 1)
    ax1.plot(epochs_x, history["train_loss"], label="Train")
    ax1.plot(epochs_x, history["val_loss"], label="Val")
    ax1.set_xlabel("Epoch")
    ax1.set_ylabel("SupCon Loss")
    ax1.legend()
    ax1.set_title("Training Curves")
    sns.despine(ax=ax1)

    ax2.plot(epochs_x, history["lr"])
    ax2.set_xlabel("Epoch")
    ax2.set_ylabel("Learning Rate")
    ax2.set_title("LR Schedule")
    sns.despine(ax=ax2)

    fig.tight_layout()
    fig.savefig(fig_dir / "supcon_training_curves.png",
                dpi=config["fig_dpi"], bbox_inches="tight")
    plt.close()

    print(f"\nTraining complete. Best val loss: {best_val_loss:.4f}")
    print(f"Checkpoint: {ckpt_dir / 'best_supcon.pth'}")

    return model, label_to_idx


# ==============================================================================
# EMBEDDING EXTRACTION
# ==============================================================================

def extract_embeddings(config):
    """
    Load best checkpoint and extract 512-d embeddings for ALL cells.

    Saves embeddings + metadata to HDF5.
    """
    device = get_device()
    ckpt_path = Path(config["checkpoint_dir"]) / "best_supcon.pth"

    print(f"Loading checkpoint: {ckpt_path}")
    ckpt = torch.load(str(ckpt_path), map_location=device, weights_only=False)
    label_to_idx = ckpt["label_to_idx"]

    model = SupConResNet(
        in_channels=config["supcon_input_channels"],
        embedding_dim=config["embedding_dim"],
        projection_dim=config["projection_dim"],
        pretrained=False,
    ).to(device)
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()

    # Load ALL cells (including controls)
    h5_path = config["crop_output_h5"]
    dataset = CropHDF5Dataset(h5_path, label_to_idx=label_to_idx, transform=None)

    loader = DataLoader(
        dataset, batch_size=config["batch_size"], shuffle=False,
        num_workers=0,
    )

    all_embeddings = []
    print(f"Extracting embeddings for {len(dataset):,} cells ...")

    with torch.no_grad():
        for i, (crops, _) in enumerate(loader):
            crops = crops.to(device)
            features = model.encode(crops)  # (batch, 512)
            all_embeddings.append(features.cpu().numpy())
            if (i + 1) % 100 == 0:
                print(f"  Batch {i + 1}/{len(loader)}")

    embeddings = np.concatenate(all_embeddings, axis=0)  # (N, 512)
    print(f"Embeddings shape: {embeddings.shape}")

    # Save
    emb_dir = Path(config["embedding_dir"])
    emb_dir.mkdir(parents=True, exist_ok=True)
    emb_path = emb_dir / "cell_embeddings.h5"

    str_dt = h5py.string_dtype()
    with h5py.File(str(h5_path), "r") as src, \
         h5py.File(str(emb_path), "w") as dst:
        dst.create_dataset("embeddings", data=embeddings,
                           chunks=(min(1000, len(embeddings)), 512),
                           compression="gzip")
        # Copy metadata from crop H5
        for key in ["condition_labels", "reporters", "experiment_types",
                    "knockdowns", "tiff_files", "fov_ids",
                    "is_control", "is_drug", "areas_px"]:
            if key in src:
                dst.create_dataset(key, data=src[key][:])

        dst.attrs["embedding_dim"] = config["embedding_dim"]
        dst.attrs["checkpoint"] = str(ckpt_path)

    print(f"Saved to {emb_path}")
    return embeddings


# ==============================================================================
# OT DISTANCE COMPUTATION
# ==============================================================================

def compute_ot_distances(config):
    """
    Compute pairwise Sinkhorn distances between non-control conditions
    using cell-level embeddings.
    """
    import ot

    emb_path = Path(config["embedding_dir"]) / "cell_embeddings.h5"
    print(f"Loading embeddings from {emb_path}")

    with h5py.File(str(emb_path), "r") as f:
        embeddings = f["embeddings"][:]
        labels = f["condition_labels"][:].astype(str)
        is_control = f["is_control"][:]
        is_drug = f["is_drug"][:]

    # Optional PCA dimensionality reduction
    pca_dims = config.get("pca_dims")
    if pca_dims and pca_dims < embeddings.shape[1]:
        from sklearn.decomposition import PCA
        print(f"Applying PCA: {embeddings.shape[1]} -> {pca_dims} dimensions")
        pca = PCA(n_components=pca_dims, random_state=config["random_seed"])
        embeddings = pca.fit_transform(embeddings)
        var_explained = pca.explained_variance_ratio_.sum()
        print(f"  Variance explained: {var_explained:.1%}")

    # Group by condition (non-control only)
    non_ctrl_mask = ~is_control
    conditions = sorted(set(labels[non_ctrl_mask]))

    rng = np.random.default_rng(config["random_seed"])
    n_sub = config["subsample_n"]
    n = len(conditions)

    # Pre-extract and subsample
    arrays = {}
    for cond in conditions:
        mask = (labels == cond) & non_ctrl_mask
        X = embeddings[mask]
        if len(X) == 0:
            raise ValueError(f"No cells for condition '{cond}'")
        if len(X) >= n_sub:
            idx = rng.choice(len(X), n_sub, replace=False)
        else:
            idx = rng.choice(len(X), n_sub, replace=True)
        arrays[cond] = X[idx].astype(np.float64)

    # Pairwise Sinkhorn distances
    dist_matrix = np.zeros((n, n))
    n_pairs = n * (n - 1) // 2
    print(f"Computing {n_pairs} pairwise Sinkhorn distances "
          f"({n} conditions, {n_sub} cells each) ...")

    pair_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            X_a, X_b = arrays[conditions[i]], arrays[conditions[j]]
            M = ot.dist(X_a, X_b, metric="sqeuclidean")
            M_max = M.max()
            if M_max > 0:
                M /= M_max
            a = np.ones(len(X_a)) / len(X_a)
            b = np.ones(len(X_b)) / len(X_b)
            w = ot.sinkhorn2(
                a, b, M,
                reg=config["sinkhorn_reg"],
                numItermax=config["sinkhorn_max_iter"],
                warn=False,
            )
            dist_matrix[i, j] = dist_matrix[j, i] = float(w)
            pair_count += 1
            if pair_count % 50 == 0 or pair_count == n_pairs:
                print(f"  {pair_count}/{n_pairs} pairs done")

    dist_df = pd.DataFrame(dist_matrix, index=conditions, columns=conditions)

    # Identify drug vs gene conditions
    drug_conds = set()
    for cond in conditions:
        mask = labels == cond
        if is_drug[mask].any():
            drug_conds.add(cond)

    dist_df.to_csv(config["distance_csv"])
    print(f"Distance matrix saved to {config['distance_csv']}")

    return dist_df, sorted(drug_conds)


# ==============================================================================
# DRUG–GENE RANKING
# ==============================================================================

def rank_matches(dist_df, drug_conditions, gene_conditions, config):
    """Rank gene knockdowns by ascending distance for each drug,
    with pathway annotation if gene_pathway_map.csv is available."""
    gene_to_pathway = load_pathway_map(config) or {}

    rows = []
    for drug in sorted(drug_conditions):
        distances = dist_df.loc[drug, gene_conditions].sort_values()
        for rank, (gene, dist) in enumerate(distances.items(), 1):
            rows.append({
                "drug": drug,
                "rank": rank,
                "gene": gene,
                "pathway": gene_to_pathway.get(gene, ""),
                "distance": dist,
            })
    return pd.DataFrame(rows)


# ==============================================================================
# PERMUTATION TEST
# ==============================================================================

def permutation_test(config, observed_dist_df, drug_conditions, gene_conditions):
    """
    Permutation test: shuffle condition labels at cell level, recompute
    OT distances, test whether observed drug–gene distances are smaller
    than expected by chance.
    """
    import ot

    emb_path = Path(config["embedding_dir"]) / "cell_embeddings.h5"
    with h5py.File(str(emb_path), "r") as f:
        embeddings = f["embeddings"][:]
        labels = f["condition_labels"][:].astype(str)
        is_control = f["is_control"][:]

    # Apply same PCA reduction as in compute_ot_distances
    pca_dims = config.get("pca_dims")
    if pca_dims and pca_dims < embeddings.shape[1]:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=pca_dims, random_state=config["random_seed"])
        embeddings = pca.fit_transform(embeddings)

    non_ctrl_mask = ~is_control
    non_ctrl_emb = embeddings[non_ctrl_mask]
    non_ctrl_labels = labels[non_ctrl_mask]

    rng = np.random.default_rng(config["random_seed"])
    top_k = config["permutation_top_k"]
    n_perm = config["n_permutations"]
    n_sub = config["permutation_subsample_n"]

    all_conditions = sorted(set(drug_conditions) | set(gene_conditions))

    # Pairs to test
    pairs_to_test = []
    for drug in drug_conditions:
        top_genes = observed_dist_df.loc[drug, gene_conditions].nsmallest(top_k).index
        for gene in top_genes:
            pairs_to_test.append((drug, gene))

    observed = {(d, g): observed_dist_df.loc[d, g] for d, g in pairs_to_test}
    counts = {pair: 0 for pair in pairs_to_test}

    needed = set()
    for d, g in pairs_to_test:
        needed.add(d)
        needed.add(g)

    print(f"Running {n_perm} permutations for {len(pairs_to_test)} "
          f"drug–gene pairs ...")

    for perm_i in range(n_perm):
        # Shuffle condition labels
        shuffled_labels = non_ctrl_labels.copy()
        rng.shuffle(shuffled_labels)

        # Subsample per condition
        arrays = {}
        for cond in needed:
            mask = shuffled_labels == cond
            X = non_ctrl_emb[mask]
            if len(X) == 0:
                continue
            if len(X) >= n_sub:
                idx = rng.choice(len(X), n_sub, replace=False)
            else:
                idx = rng.choice(len(X), n_sub, replace=True)
            arrays[cond] = X[idx].astype(np.float64)

        for drug, gene in pairs_to_test:
            if drug not in arrays or gene not in arrays:
                continue
            X_a, X_b = arrays[drug], arrays[gene]
            M = ot.dist(X_a, X_b, metric="sqeuclidean")
            M_max = M.max()
            if M_max > 0:
                M /= M_max
            a = np.ones(n_sub) / n_sub
            b = np.ones(n_sub) / n_sub
            w = float(ot.sinkhorn2(a, b, M, reg=config["sinkhorn_reg"],
                                   numItermax=config["sinkhorn_max_iter"],
                                   warn=False))
            if w <= observed[(drug, gene)]:
                counts[(drug, gene)] += 1

        if (perm_i + 1) % 50 == 0 or perm_i == n_perm - 1:
            print(f"  permutation {perm_i + 1}/{n_perm}")

    # p-values + FDR
    p_values = {pair: (counts[pair] + 1) / (n_perm + 1) for pair in pairs_to_test}
    pv_df = pd.DataFrame([
        {"drug": d, "gene": g, "p_value": p_values[(d, g)]}
        for d, g in pairs_to_test
    ])
    pv_df["p_adj"] = false_discovery_control(pv_df["p_value"].values, method="bh")
    pv_df["significant"] = pv_df["p_adj"] < config["alpha"]

    n_sig = pv_df["significant"].sum()
    print(f"Permutation test: {n_sig}/{len(pv_df)} pairs significant "
          f"(FDR < {config['alpha']})")

    return pv_df


# ==============================================================================
# PATHWAY-LEVEL ANALYSIS
# ==============================================================================

def load_pathway_map(config):
    """Load gene→pathway mapping from CSV. Returns dict {gene: pathway}."""
    csv_path = config.get("gene_pathway_csv")
    if not csv_path or not Path(csv_path).exists():
        print("  No gene_pathway_map.csv found — skipping pathway analysis")
        return None
    df = pd.read_csv(csv_path)
    return dict(zip(df["Gene"], df["Pathway"]))


def pathway_level_analysis(config, dist_df, drug_conditions, gene_conditions):
    """
    Compute pathway-level OT distances by averaging gene-level distances
    within each pathway, then run a separate permutation test.

    Returns a DataFrame with drug–pathway ranked matches + p-values.
    """
    import ot as pot

    gene_to_pathway = load_pathway_map(config)
    if gene_to_pathway is None:
        return None, None

    # Build pathway → gene list (only for genes in our data)
    pathway_genes = {}
    for gene in gene_conditions:
        pw = gene_to_pathway.get(gene)
        if pw:
            pathway_genes.setdefault(pw, []).append(gene)

    pathways = sorted(pathway_genes.keys())
    print(f"\nPathway-level analysis: {len(pathways)} pathways")
    for pw in pathways:
        print(f"  {pw}: {pathway_genes[pw]}")

    # ── Observed pathway distances (mean of gene distances) ──────────────
    pw_rows = []
    for drug in sorted(drug_conditions):
        for pw in pathways:
            genes = pathway_genes[pw]
            mean_dist = dist_df.loc[drug, genes].mean()
            pw_rows.append({
                "drug": drug,
                "pathway": pw,
                "distance": mean_dist,
                "genes": ", ".join(genes),
                "n_genes": len(genes),
            })

    pw_df = pd.DataFrame(pw_rows)

    # Add rank per drug
    pw_df["rank"] = pw_df.groupby("drug")["distance"].rank(method="min").astype(int)
    pw_df = pw_df.sort_values(["drug", "rank"])

    # ── Permutation test at pathway level ────────────────────────────────
    n_perm = config["n_permutations"]
    if n_perm <= 0:
        return pw_df, gene_to_pathway

    emb_path = Path(config["embedding_dir"]) / "cell_embeddings.h5"
    with h5py.File(str(emb_path), "r") as f:
        embeddings = f["embeddings"][:]
        labels = f["condition_labels"][:].astype(str)
        is_control = f["is_control"][:]

    # Apply same PCA reduction
    pca_dims = config.get("pca_dims")
    if pca_dims and pca_dims < embeddings.shape[1]:
        from sklearn.decomposition import PCA
        pca = PCA(n_components=pca_dims, random_state=config["random_seed"])
        embeddings = pca.fit_transform(embeddings)

    non_ctrl_mask = ~is_control
    non_ctrl_emb = embeddings[non_ctrl_mask]
    non_ctrl_labels = labels[non_ctrl_mask]

    rng = np.random.default_rng(config["random_seed"])
    n_sub = config["permutation_subsample_n"]

    # Only test rank-1 pathway per drug
    top_k_pw = config.get("permutation_top_k", 1)
    pairs_to_test = []
    for drug in sorted(drug_conditions):
        drug_pw = pw_df[pw_df["drug"] == drug].nsmallest(top_k_pw, "distance")
        for _, row in drug_pw.iterrows():
            pairs_to_test.append((drug, row["pathway"]))

    observed = {}
    for drug, pw in pairs_to_test:
        mask = (pw_df["drug"] == drug) & (pw_df["pathway"] == pw)
        observed[(drug, pw)] = pw_df.loc[mask, "distance"].values[0]

    counts = {pair: 0 for pair in pairs_to_test}

    # All conditions needed for permutation
    all_needed = set()
    for d, _ in pairs_to_test:
        all_needed.add(d)
    for pw in set(p for _, p in pairs_to_test):
        for g in pathway_genes[pw]:
            all_needed.add(g)

    print(f"Running {n_perm} pathway-level permutations for "
          f"{len(pairs_to_test)} drug–pathway pairs ...")

    for perm_i in range(n_perm):
        shuffled_labels = non_ctrl_labels.copy()
        rng.shuffle(shuffled_labels)

        # Compute perm distances for needed conditions
        arrays = {}
        for cond in all_needed:
            mask = shuffled_labels == cond
            X = non_ctrl_emb[mask]
            if len(X) == 0:
                continue
            if len(X) >= n_sub:
                idx = rng.choice(len(X), n_sub, replace=False)
            else:
                idx = rng.choice(len(X), n_sub, replace=True)
            arrays[cond] = X[idx].astype(np.float64)

        # Compute perm OT distances for drug–gene pairs, then average by pathway
        for drug, pw in pairs_to_test:
            if drug not in arrays:
                continue
            gene_dists = []
            for gene in pathway_genes[pw]:
                if gene not in arrays:
                    continue
                X_a, X_b = arrays[drug], arrays[gene]
                M = pot.dist(X_a, X_b, metric="sqeuclidean")
                M_max = M.max()
                if M_max > 0:
                    M /= M_max
                a = np.ones(n_sub) / n_sub
                b = np.ones(n_sub) / n_sub
                w = float(pot.sinkhorn2(a, b, M, reg=config["sinkhorn_reg"],
                                        numItermax=config["sinkhorn_max_iter"],
                                        warn=False))
                gene_dists.append(w)

            if gene_dists:
                perm_pw_dist = np.mean(gene_dists)
                if perm_pw_dist <= observed[(drug, pw)]:
                    counts[(drug, pw)] += 1

        if (perm_i + 1) % 500 == 0 or perm_i == n_perm - 1:
            print(f"  permutation {perm_i + 1}/{n_perm}")

    # p-values + FDR
    p_values = {pair: (counts[pair] + 1) / (n_perm + 1) for pair in pairs_to_test}
    for drug, pw in pairs_to_test:
        mask = (pw_df["drug"] == drug) & (pw_df["pathway"] == pw)
        pw_df.loc[mask, "p_value"] = p_values[(drug, pw)]

    tested_mask = pw_df["p_value"].notna()
    if tested_mask.any():
        pw_df.loc[tested_mask, "p_adj"] = false_discovery_control(
            pw_df.loc[tested_mask, "p_value"].values, method="bh"
        )
        pw_df["significant"] = pw_df["p_adj"] < config["alpha"]
        n_sig = pw_df["significant"].sum()
        print(f"Pathway permutation test: {n_sig}/{tested_mask.sum()} "
              f"pairs significant (FDR < {config['alpha']})")

    return pw_df, gene_to_pathway


# ==============================================================================
# VISUALISATION
# ==============================================================================

def generate_plots(dist_df, match_table, drug_conditions, gene_conditions,
                   config):
    """Generate heatmap, UMAPs, and ranked-match charts."""
    import umap

    fig_dir = Path(config["fig_dir"])
    fig_dir.mkdir(parents=True, exist_ok=True)
    dpi = config["fig_dpi"]

    palette = {"gene": "#4C72B0", "drug": "#DD8452"}

    # ── 1. Distance heatmap with clustering ──────────────────────────────
    print("Generating distance heatmap ...")
    cond_type = pd.Series("gene", index=dist_df.index)
    cond_type[cond_type.index.isin(drug_conditions)] = "drug"
    row_colors = cond_type.map(palette).rename("type")

    g = sns.clustermap(
        dist_df, method="ward", cmap="viridis_r",
        figsize=(12, 10),
        row_colors=row_colors, col_colors=row_colors,
        linewidths=0, xticklabels=True, yticklabels=True,
        dendrogram_ratio=(0.12, 0.12),
        cbar_pos=(0.02, 0.82, 0.03, 0.15),
    )
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.tick_params(labelsize=8)
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=palette["gene"], label="Gene KD"),
                       Patch(facecolor=palette["drug"], label="Drug")]
    g.ax_heatmap.legend(handles=legend_elements, loc="upper left",
                        bbox_to_anchor=(1.02, 1), frameon=False, fontsize=9)
    g.savefig(fig_dir / "supcon_distance_heatmap.png", dpi=dpi,
              bbox_inches="tight")
    plt.close()

    # ── 2. Cell-level UMAP ───────────────────────────────────────────────
    print("Generating cell-level UMAP ...")
    emb_path = Path(config["embedding_dir"]) / "cell_embeddings.h5"

    with h5py.File(str(emb_path), "r") as f:
        embeddings = f["embeddings"][:]
        labels = f["condition_labels"][:].astype(str)
        is_control = f["is_control"][:]

    # Sample ~5000 non-control cells for visualisation
    non_ctrl_idx = np.where(~is_control)[0]
    rng = np.random.default_rng(config["random_seed"])
    n_sample = min(5000, len(non_ctrl_idx))
    sample_idx = rng.choice(non_ctrl_idx, n_sample, replace=False)

    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1,
                        random_state=config["random_seed"])
    umap_emb = reducer.fit_transform(embeddings[sample_idx])
    sample_labels = labels[sample_idx]

    # Color by condition
    unique_labels = sorted(set(sample_labels))
    cmap = plt.cm.get_cmap("tab20", len(unique_labels))
    label_colors = {lbl: cmap(i) for i, lbl in enumerate(unique_labels)}

    fig, ax = plt.subplots(figsize=(12, 10))
    for lbl in unique_labels:
        mask = sample_labels == lbl
        ax.scatter(umap_emb[mask, 0], umap_emb[mask, 1],
                   c=[label_colors[lbl]], s=8, alpha=0.5, label=lbl)
    ax.legend(fontsize=6, ncol=3, markerscale=3, frameon=False,
              loc="upper left", bbox_to_anchor=(1.02, 1))
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("Cell-level embeddings (SupCon)")
    sns.despine(ax=ax)
    fig.savefig(fig_dir / "supcon_umap_cells.png", dpi=dpi, bbox_inches="tight")
    plt.close()

    # ── 3. Condition-level UMAP from distance matrix ─────────────────────
    print("Generating condition UMAP ...")
    n_conditions = len(dist_df)
    n_neighbors = min(5, n_conditions - 1)

    reducer = umap.UMAP(
        metric="precomputed", n_neighbors=n_neighbors, min_dist=0.3,
        random_state=config["random_seed"],
    )
    cond_umap = reducer.fit_transform(dist_df.values)

    fig, ax = plt.subplots(figsize=(10, 8))
    for ctype, color in palette.items():
        mask = cond_type == ctype
        ax.scatter(cond_umap[mask, 0], cond_umap[mask, 1],
                   c=color, s=80,
                   label="Gene KD" if ctype == "gene" else "Drug",
                   edgecolors="white", linewidth=0.5, zorder=3)
    for i, label in enumerate(dist_df.index):
        ax.annotate(label, (cond_umap[i, 0], cond_umap[i, 1]),
                    fontsize=7, ha="center", va="bottom",
                    xytext=(0, 5), textcoords="offset points")
    ax.legend(frameon=False)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title("Condition-level UMAP (SupCon + Wasserstein)")
    sns.despine(ax=ax)
    fig.savefig(fig_dir / "supcon_umap_conditions.png", dpi=dpi,
                bbox_inches="tight")
    plt.close()

    # ── 4. Ranked match bar charts ───────────────────────────────────────
    print("Generating ranked match chart ...")
    top_n = min(5, len(gene_conditions))
    top_matches = match_table[match_table["rank"] <= top_n].copy()
    n_drugs = top_matches["drug"].nunique()

    fig, axes = plt.subplots(1, n_drugs, figsize=(4 * n_drugs, 5), sharey=False)
    if n_drugs == 1:
        axes = [axes]

    for ax, (drug, grp) in zip(axes, top_matches.groupby("drug")):
        grp = grp.sort_values("distance")
        if "significant" in grp.columns:
            colors = ["#55A868" if sig else "#CCCCCC"
                      for sig in grp["significant"]]
        else:
            colors = ["#DD8452"] * len(grp)
        ax.barh(grp["gene"], grp["distance"], color=colors, edgecolor="white")
        ax.set_xlabel("Wasserstein distance")
        ax.set_title(drug, fontweight="bold")
        ax.invert_yaxis()
        sns.despine(ax=ax)

    fig.suptitle("Top gene matches per drug (SupCon + OT)",
                 fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(fig_dir / "supcon_ranked_matches.png", dpi=dpi,
                bbox_inches="tight")
    plt.close()

    print(f"Figures saved to {fig_dir}/")


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="SupCon training + OT morphological profiling")
    parser.add_argument("--train-only", action="store_true",
                        help="Only train the model")
    parser.add_argument("--embed-only", action="store_true",
                        help="Only extract embeddings (needs checkpoint)")
    parser.add_argument("--ot-only", action="store_true",
                        help="Only run OT + permutation (needs embeddings)")
    args = parser.parse_args()

    print("=" * 60)
    print("SupCon + OT Morphological Profiling Pipeline")
    print("=" * 60)

    run_all = not (args.train_only or args.embed_only or args.ot_only)

    # ── Train ────────────────────────────────────────────────────────────
    if run_all or args.train_only:
        print("\n── Phase 2A: SupCon Training ──")
        train_supcon(CONFIG)

    # ── Extract embeddings ───────────────────────────────────────────────
    if run_all or args.embed_only:
        print("\n── Phase 2B: Embedding Extraction ──")
        extract_embeddings(CONFIG)

    # ── OT distances + permutation + plots ───────────────────────────────
    if run_all or args.ot_only:
        print("\n── Phase 2C: OT Distance Computation ──")
        dist_df, drug_conds = compute_ot_distances(CONFIG)

        all_conditions = dist_df.index.tolist()
        gene_conds = sorted(set(all_conditions) - set(drug_conds))

        print(f"\nDrug conditions ({len(drug_conds)}): {drug_conds}")
        print(f"Gene conditions ({len(gene_conds)}): {gene_conds}")

        # Rank (with pathway annotation)
        match_table = rank_matches(dist_df, drug_conds, gene_conds, CONFIG)

        # Permutation test (gene-level)
        if CONFIG["n_permutations"] > 0:
            print("\n── Phase 2D: Permutation Test (gene-level) ──")
            pv_df = permutation_test(CONFIG, dist_df, drug_conds, gene_conds)
            match_table = match_table.merge(pv_df, on=["drug", "gene"],
                                            how="left")

        # Save gene-level results
        match_table.to_csv(CONFIG["output_csv"], index=False)
        print(f"\nRanked matches saved to {CONFIG['output_csv']}")
        print("\nTop match per drug:")
        for _, row in match_table[match_table["rank"] == 1].iterrows():
            sig = ""
            if "p_adj" in row and pd.notna(row["p_adj"]):
                sig = f"  (p_adj={row['p_adj']:.4f})"
            pw = f"  [{row['pathway']}]" if row.get("pathway") else ""
            print(f"  {row['drug']:15s} -> {row['gene']:15s}  "
                  f"dist={row['distance']:.4f}{sig}{pw}")

        # Pathway-level analysis
        print("\n── Phase 2D+: Pathway-level Analysis ──")
        pw_df, _ = pathway_level_analysis(
            CONFIG, dist_df, drug_conds, gene_conds)
        if pw_df is not None:
            pw_df.to_csv(CONFIG["pathway_output_csv"], index=False)
            print(f"\nPathway matches saved to {CONFIG['pathway_output_csv']}")
            print("\nTop pathway per drug:")
            for _, row in pw_df[pw_df["rank"] == 1].iterrows():
                sig = ""
                if "p_adj" in row and pd.notna(row["p_adj"]):
                    sig = f"  (p_adj={row['p_adj']:.4f})"
                print(f"  {row['drug']:15s} -> {row['pathway']:30s}  "
                      f"dist={row['distance']:.4f}{sig}  "
                      f"[{row['genes']}]")

        # Plots
        print("\n── Phase 2E: Visualisation ──")
        generate_plots(dist_df, match_table, drug_conds, gene_conds, CONFIG)

    print("\n" + "=" * 60)
    print("Pipeline complete.")
    print("=" * 60)


if __name__ == "__main__":
    main()
