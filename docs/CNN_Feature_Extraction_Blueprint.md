# CNN Feature Extraction Blueprint for Mycobacterial Morphological Profiling

## Goal
Train a self-supervised model on individual cell images to learn a rich feature
representation. Use these features to build a reference atlas of known
perturbations (knockdowns/drugs), then predict which known perturbation a new
sample most resembles by projecting it into the same feature space.

## Dataset
- **Organism**: M. smegmatis (test dataset), M. tuberculosis (target)
- **Scale**: ~138GB images, ~260 strains
- **Available per cell**: full field-of-view images, segmentation masks (filtered
  through a classifier), coordinate list
- **Channels**: brightfield/phase + fluorescence reporter (ch1), possibly ch2/ch3

---

## Architecture: Pretrained Encoder + Autoencoder

```
Input: 2-channel cell thumbnail (brightfield + fluorescence)
       e.g. 128×128×2

  ┌─────────────────────────────────────────────┐
  │  ENCODER (ResNet-18, pretrained on ImageNet) │
  │  - Replace first conv layer: 3ch → 2ch      │
  │  - Freeze early layers (first 2-3 blocks)   │
  │  - Fine-tune later layers                    │
  │  - Output: 512-d feature vector              │
  └──────────────────┬──────────────────────────┘
                     │
              512-d bottleneck  ← THIS IS YOUR FEATURE VECTOR
                     │
  ┌──────────────────┴──────────────────────────┐
  │  DECODER (transposed convolutions)           │
  │  - Upsample back to 128×128×2               │
  │  - Reconstruct the input image              │
  └─────────────────────────────────────────────┘

Loss: MSE between input and reconstruction
      (optionally + perceptual loss for sharper reconstructions)
```

### Why this architecture
- **Autoencoder** forces the 512-d bottleneck to contain enough information to
  reconstruct the cell image — so it must capture morphology and intensity
- **Pretrained encoder** gives you good low-level features (edges, textures)
  for free, even though ImageNet has no bacteria
- **Fine-tuning later layers** adapts the high-level features to your domain
- **Simple to implement**, easy to debug, interpretable (you can look at
  reconstructions to see what the model has learned)

---

## Pipeline

### Phase 1: Preprocessing

```
Full field-of-view images + coordinate list + masks
    │
    ▼
For each detected cell:
    1. Get centroid from coordinate list
    2. Crop fixed-size patch (128×128) from brightfield + fluorescence channels
    3. Optionally: multiply by dilated mask to zero out neighbouring cells
    4. Normalise per-channel (zero mean, unit variance per image or per dataset)
    5. Save as .npy or HDF5 (faster than loading individual PNGs)
    │
    ▼
Cell thumbnail dataset (~millions of 128×128×2 images)
with metadata: strain, experiment, field-of-view, cell_id
```

**Decisions:**
- **Crop size**: should comfortably contain the largest cells with some margin.
  Check your typical cell length in pixels. 64×64 may suffice if cells are
  small; 128×128 is safer. Larger = more context but slower training.
- **Mask multiplication**: pros — isolates the target cell, removes neighbours.
  Cons — destroys edge information. Recommendation: try both, compare
  reconstruction quality.
- **Storage format**: with 138GB of raw images and potentially millions of
  crops, use HDF5 with chunked storage. A single HDF5 file with all crops
  indexed by strain/cell_id is much faster than individual files.

### Phase 2: Model Training

```python
import torch
import torch.nn as nn
from torchvision.models import resnet18, ResNet18_Weights

class CellAutoencoder(nn.Module):
    def __init__(self, feature_dim=512):
        super().__init__()

        # Encoder: pretrained ResNet-18 with modified input
        encoder = resnet18(weights=ResNet18_Weights.DEFAULT)

        # Replace first conv: 3 channels → 2 channels
        # Initialise by averaging the RGB weights into 2 channels
        old_conv = encoder.conv1
        new_conv = nn.Conv2d(2, 64, kernel_size=7, stride=2, padding=3, bias=False)
        with torch.no_grad():
            # Use mean of RGB weights for brightfield, green for fluorescence
            new_conv.weight[:, 0] = old_conv.weight.mean(dim=1)
            new_conv.weight[:, 1] = old_conv.weight[:, 1]
        encoder.conv1 = new_conv

        # Remove the final FC layer — we want the 512-d feature vector
        self.encoder = nn.Sequential(*list(encoder.children())[:-1])
        # Output: (batch, 512, 1, 1) after adaptive avg pool

        # Decoder: upsample back to 128×128×2
        self.decoder = nn.Sequential(
            nn.ConvTranspose2d(512, 256, 4, 2, 1),  # 1→2 (or 4→8 etc.)
            nn.BatchNorm2d(256),
            nn.ReLU(),
            nn.ConvTranspose2d(256, 128, 4, 2, 1),  # →4 (or →16)
            nn.BatchNorm2d(128),
            nn.ReLU(),
            nn.ConvTranspose2d(128, 64, 4, 2, 1),   # →8 (or →32)
            nn.BatchNorm2d(64),
            nn.ReLU(),
            nn.ConvTranspose2d(64, 32, 4, 2, 1),    # →16 (or →64)
            nn.BatchNorm2d(32),
            nn.ReLU(),
            nn.ConvTranspose2d(32, 2, 4, 2, 1),     # →32 (or →128)
            nn.Sigmoid()  # if input is normalised to [0,1]
        )
        # NOTE: exact upsampling steps depend on the spatial size after encoder
        # ResNet-18 with 128×128 input → 4×4 after the last block
        # So you need: 4→8→16→32→64→128 = 5 upsample steps

    def encode(self, x):
        h = self.encoder(x)
        return h.squeeze(-1).squeeze(-1)  # (batch, 512)

    def forward(self, x):
        h = self.encoder(x)             # (batch, 512, 1, 1)
        # Reshape for decoder — need spatial dimensions
        h_spatial = h.expand(-1, -1, 4, 4)  # broadcast to 4×4
        reconstruction = self.decoder(h_spatial)
        features = h.squeeze(-1).squeeze(-1)
        return reconstruction, features
```

**Training details:**
- **Optimiser**: AdamW, lr=1e-4, weight_decay=1e-5
- **Schedule**: cosine annealing over 50-100 epochs
- **Batch size**: 256-512 (adjust to GPU memory)
- **Freeze strategy**: freeze encoder blocks 1-2 for first 10 epochs, then
  unfreeze all. This lets the decoder learn to work with pretrained features
  before fine-tuning the encoder.
- **Data augmentation**: random rotation (bacteria have no canonical
  orientation), random flip, small brightness/contrast jitter. Do NOT use
  aggressive colour augmentation — your fluorescence intensity is meaningful.
- **Validation**: hold out ~10% of cells (entire fields of view, not random
  cells, to avoid data leakage from neighbouring cells in the same image)

**Training time estimate**: with ~1M cell crops at 128×128, ResNet-18 encoder,
batch size 512 on a single modern GPU (A100/4090): roughly 2-4 hours for
50 epochs. The 138GB raw data will produce a much larger crop dataset, but
training is on the crops, not the raw images.

### Phase 3: Feature Extraction

```python
model.eval()
all_features = []
all_metadata = []

with torch.no_grad():
    for batch_images, batch_meta in dataloader:
        features = model.encode(batch_images.cuda())  # (batch, 512)
        all_features.append(features.cpu().numpy())
        all_metadata.extend(batch_meta)

# Shape: (n_cells, 512)
feature_matrix = np.vstack(all_features)
```

### Phase 4: Sample-Level Summarisation

```python
import pandas as pd

df = pd.DataFrame(feature_matrix, columns=[f"feat_{i}" for i in range(512)])
df["strain"] = [m["strain"] for m in all_metadata]

# Mean feature vector per strain (like your current S-score means)
strain_profiles = df.groupby("strain").mean()

# Optional: also compute per-strain standard deviation or other statistics
strain_std = df.groupby("strain").std()
```

This gives you a (260 strains × 512 features) matrix — analogous to your
current S-score matrix, but with learned features.

### Phase 5: Downstream Analysis (slots into existing R pipeline)

Export `strain_profiles` as a CSV. Then in R:

```r
# Load CNN features instead of computing S-scores
cnn_features <- read_csv("cnn_strain_profiles.csv")
feature_mat <- as.matrix(cnn_features %>% select(-strain))
rownames(feature_mat) <- cnn_features$strain

# Everything from here is identical to MorphologicalProfiling_PCA_Pipeline.R:
# PCA → Ward's D2 → silhouette/gap → visualise → bootstrap
pca <- prcomp(feature_mat, center = TRUE, scale. = TRUE)
# ... etc
```

### Phase 6: Predicting New Samples

```python
# Extract features for cells from a new strain/drug treatment
new_features = model.encode(new_cell_images)        # (n_cells, 512)
new_profile = new_features.mean(dim=0, keepdim=True)  # (1, 512)

# Find nearest known strains
from sklearn.metrics.pairwise import cosine_similarity
similarities = cosine_similarity(new_profile, strain_profiles)
top_matches = similarities.argsort()[0][::-1][:5]
```

Or in R, project into the existing PCA space and find nearest cluster
centroid — exactly as Mode A currently works.

---

## Quality Checks

1. **Reconstruction quality**: visually inspect input vs reconstruction for
   random cells. If reconstructions are blurry but capture overall shape and
   fluorescence intensity, the bottleneck is working. If they're garbage, the
   model isn't learning.

2. **Feature space structure**: run UMAP on the 512-d features (display only)
   coloured by strain. If known-similar strains cluster together without
   supervision, the features are biologically meaningful.

3. **Retrieval test**: for each strain, find its 5 nearest neighbours in
   feature space. Do they make biological sense? Do target-gene-related
   knockdowns group together?

4. **Comparison to hand-crafted features**: run your existing PCA pipeline on
   both the CNN features and the SHAPE S-scores. Do they agree? Where do they
   disagree? Disagreements are the interesting cases — either the CNN is
   picking up something SHAPE misses, or it's capturing noise.

5. **Held-out strain prediction**: leave out 10% of strains entirely during
   training. After training, extract features for the held-out strains and
   check whether they map to biologically sensible neighbours. This tests
   generalisation.

---

## Practical Considerations

### Compute
- **Training**: single GPU, 2-4 hours. Can be done on a university HPC node
  or cloud (Google Colab Pro would work for prototyping)
- **Feature extraction**: ~30 min for 1M cells on GPU
- **Preprocessing (cropping)**: CPU-bound, parallelisable. Budget 1-2 hours
  for 138GB of images with multiprocessing

### Storage
- Raw images: 138GB (existing)
- Cell crops (HDF5): estimate ~20-50GB depending on crop size and number of
  cells. With 260 strains and potentially 5-10K cells per strain, that is
  ~1-2.5M crops at 128×128×2 × float32 = ~125-300GB uncompressed.
  With HDF5 LZF compression: ~30-80GB.
- Feature vectors: trivial (~5MB for 1M cells × 512 features × float32,
  or ~500KB for 260 strain-level profiles)

### Python Packages
- PyTorch + torchvision (model, training)
- h5py (HDF5 crop storage)
- scikit-image or PIL (image loading, cropping)
- numpy, pandas (data handling)
- scikit-learn (optional: for nearest neighbours, UMAP via umap-learn)

### Scaling to M. tuberculosis
The model trained on M. smegmatis will need to be evaluated on M. tb before
assuming transfer. Smeg is faster-growing and morphologically different.
Options:
- **Direct transfer**: extract features using the smeg-trained model on tb
  images. Test whether the features still separate known tb knockdowns.
- **Fine-tune**: take the smeg-trained model and fine-tune on tb data
  (fewer epochs, lower learning rate). This is the most likely path.
- **Train from scratch on tb**: only if the above fail. Requires sufficient
  tb imaging data.

---

## Summary: What to Do First

1. Pick one strain set (or a subset of ~20-30 strains) from the 138GB smeg
   dataset for prototyping
2. Write the crop extraction script (coordinates + masks → HDF5)
3. Train the autoencoder (~2-4 hours on GPU)
4. Extract features, compute strain profiles, export to CSV
5. Run through your existing PCA + clustering R pipeline
6. Compare to SHAPE-based results — does the CNN find the same clusters?
   Different ones? More structure?
7. If it works: scale to all 260 strains, then test transfer to M. tb
