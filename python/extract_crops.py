#!/usr/bin/env python3
"""
Phase 1: Crop Extraction (M. smegmatis)
========================================
Walks all segmented TIFFs, extracts 128×128 single-cell crops
(phase + ParB fluorescence + mask), and stores them in a single HDF5
file with condition metadata.

The segmented TIFFs already have quality-filtered masks (QC was applied
during the segmentation pipeline), so no additional filtering is needed.

Usage:
    python extract_crops.py
    python extract_crops.py --dry-run         # just print stats, no H5

Requirements:
    pip install tifffile numpy scikit-image h5py
"""

import argparse
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
from skimage.measure import label as sk_label
from skimage.transform import resize

from config_supcon import CONFIG

# Add ImagingPipeline to path so we can import load_hyperstack
_IMAGING_DIR = Path(CONFIG["imaging_pipeline_dir"])
if str(_IMAGING_DIR) not in sys.path:
    sys.path.insert(0, str(_IMAGING_DIR))

from label_cells import load_hyperstack

# ── Constants ────────────────────────────────────────────────────────────────

CROP_SIZE = CONFIG["crop_size"]     # 128
CROP_PAD = CONFIG["crop_pad"]       # 10
N_CHANNELS = CONFIG["supcon_input_channels"]  # 3


# ── Condition parsing ────────────────────────────────────────────────────────

# After stripping "Labelled__" and "Drugs__" prefixes:
#   Mutant:  MSMEG_XXXX_RY       (e.g. MSMEG_1234_R1)
#   Drug:    ZZZ_AXX_RY          (e.g. BDQ_025X_R1, INH_1X_R2)
#   Control: Plasmid_RY          (e.g. Plasmid_R1)
_MUTANT_RE = re.compile(r"^(MSMEG_\d+)_R(\d+)$")
_DRUG_RE = re.compile(r"^([A-Z][A-Z0-9]{1,4})_(\d+X)_R(\d+)$")
_CONTROL_RE = re.compile(r"^(Plasmid)_R(\d+)$")


def _strip_prefix(stem):
    """Strip 'Labelled__' and 'Drugs__' prefixes from filename stem."""
    if stem.startswith("Labelled__Drugs__"):
        return stem[len("Labelled__Drugs__"):]
    if stem.startswith("Labelled__"):
        return stem[len("Labelled__"):]
    return stem


def parse_condition_from_filename(stem):
    """
    Parse condition metadata from a TIFF filename stem.

    Raw filenames have prefixes (Labelled__, Labelled__Drugs__) that are
    stripped before parsing. The meaningful part follows these patterns:
        Mutants:  MSMEG_XXXX_RY   → gene accession + replica
        Drugs:    ZZZ_AX_RY       → drug code + concentration (xMIC) + replica
        Controls: Plasmid_RY      → empty vector control + replica

    Returns dict with keys: condition_label, condition_type, gene, drug,
        concentration, replica, is_control, is_drug.
    """
    name = _strip_prefix(stem)

    # Apply name corrections
    for corr in CONFIG["name_corrections"]:
        if stem == corr.get("stem"):
            name = _strip_prefix(corr.get("corrected_stem", stem))
            break

    meta = {
        "gene": "",
        "drug": "",
        "concentration": "",
        "replica": "",
        "condition_type": "",
        "condition_label": "",
        "is_control": False,
        "is_drug": False,
    }

    # Try control pattern (Plasmid_RY)
    m = _CONTROL_RE.match(name)
    if m:
        meta["replica"] = m.group(2)
        meta["condition_type"] = "control"
        meta["condition_label"] = "Plasmid"
        meta["is_control"] = True
        return meta

    # Try mutant pattern (MSMEG_XXXX_RY)
    m = _MUTANT_RE.match(name)
    if m:
        gene = m.group(1)
        meta["gene"] = gene
        meta["replica"] = m.group(2)
        meta["condition_type"] = "mutant"
        meta["condition_label"] = gene
        return meta

    # Try drug pattern (ZZZ_AX_RY)
    m = _DRUG_RE.match(name)
    if m:
        drug = m.group(1)
        concentration = m.group(2)
        meta["drug"] = drug
        meta["concentration"] = concentration
        meta["replica"] = m.group(3)
        meta["condition_type"] = "drug"
        meta["condition_label"] = f"{drug}_{concentration}"
        meta["is_drug"] = True
        return meta

    raise ValueError(
        f"Cannot parse filename '{stem}' (stripped: '{name}') — "
        f"expected MSMEG_XXXX_RY, ZZZ_AX_RY, or Plasmid_RY"
    )


# ── TIFF discovery ───────────────────────────────────────────────────────────

def discover_tiff_files():
    """
    Find all segmented TIFFs across configured directories.

    Deduplication: if the same stem exists in multiple dirs, keep the
    version from the first directory listed (segmentedData preferred).
    """
    seen_stems = {}  # stem → (path, meta)
    all_files = []

    for tiff_dir in CONFIG["tiff_dirs"]:
        tiff_dir = Path(tiff_dir)
        if not tiff_dir.exists():
            print(f"WARNING: directory not found: {tiff_dir}")
            continue

        for tiff_path in sorted(tiff_dir.glob("*.tif")):
            stem = tiff_path.stem
            if stem in seen_stems:
                continue  # keep first-seen version
            try:
                meta = parse_condition_from_filename(stem)
            except ValueError as e:
                print(f"  SKIP: {e}")
                continue
            seen_stems[stem] = (tiff_path, meta)
            all_files.append((tiff_path, meta))

    print(f"Discovered {len(all_files)} TIFFs across "
          f"{len(CONFIG['tiff_dirs'])} directories")
    return all_files


# ── Crop extraction (128×128, 3-channel) ──────────────────────────────────

def extract_cell_crop(image_channels, labeled_mask, cell_label,
                      phase_ch, fluor_ch, pad=CROP_PAD,
                      mask_background=True, dilate_px=3):
    """
    Extract a 128×128 crop of a single cell from multi-channel image data.

    Args:
        image_channels: (C, Y, X) numpy array (all channels including mask)
        labeled_mask:   (Y, X) integer array (0=bg, N=cell label)
        cell_label:     which cell to extract
        phase_ch:       index of the phase contrast channel
        fluor_ch:       index of the ParB fluorescence channel
        pad:            pixels of padding around bounding box
        mask_background: if True, replace background pixels with local mean
                         (keeps the phase halo around the cell edge via dilation)
        dilate_px:      dilation radius in pixels for the background mask

    Returns:
        crop: (3, 128, 128) float32 — [phase, parb, mask], normalised to [0, 1]
        area_px: int — true pixel area of cell in original mask
        None, None if cell_label not found
    """
    from scipy.ndimage import binary_dilation

    h, w = labeled_mask.shape
    cell_pixels = labeled_mask == cell_label

    if not cell_pixels.any():
        return None, None

    ys, xs = np.where(cell_pixels)
    y_min = max(ys.min() - pad, 0)
    y_max = min(ys.max() + pad + 1, h)
    x_min = max(xs.min() - pad, 0)
    x_max = min(xs.max() + pad + 1, w)

    # Phase channel crop
    phase_crop = image_channels[phase_ch, y_min:y_max,
                                x_min:x_max].astype(np.float32)

    # ParB fluorescence channel crop
    fluor_crop = image_channels[fluor_ch, y_min:y_max,
                                x_min:x_max].astype(np.float32)

    # Binary mask for this cell only
    mask_crop = cell_pixels[y_min:y_max, x_min:x_max]

    # Background masking: replace non-cell pixels with local mean,
    # keeping a dilated halo to preserve the phase-contrast edge
    if mask_background:
        dilated = binary_dilation(mask_crop, iterations=dilate_px)
        bg_phase = np.median(phase_crop[~dilated]) if (~dilated).any() else 0.0
        bg_fluor = np.median(fluor_crop[~dilated]) if (~dilated).any() else 0.0
        phase_crop = np.where(dilated, phase_crop, bg_phase)
        fluor_crop = np.where(dilated, fluor_crop, bg_fluor)

    mask_crop = mask_crop.astype(np.float32)

    # Stack: (3, h_crop, w_crop)
    combined = np.stack([phase_crop, fluor_crop, mask_crop], axis=0)

    # Pad to square
    _, ch, cw = combined.shape
    side = max(ch, cw)
    padded = np.zeros((3, side, side), dtype=np.float32)
    y_off = (side - ch) // 2
    x_off = (side - cw) // 2
    padded[:, y_off:y_off + ch, x_off:x_off + cw] = combined

    # Resize to CROP_SIZE × CROP_SIZE
    resized = np.zeros((3, CROP_SIZE, CROP_SIZE), dtype=np.float32)
    for c in range(3):
        resized[c] = resize(padded[c], (CROP_SIZE, CROP_SIZE),
                            order=1, preserve_range=True, anti_aliasing=True)

    # Normalise phase channel to [0, 1]
    pmin, pmax = resized[0].min(), resized[0].max()
    if pmax > pmin:
        resized[0] = (resized[0] - pmin) / (pmax - pmin)
    else:
        resized[0] = 0.0

    # Normalise fluorescence channel to [0, 1]
    fmin, fmax = resized[1].min(), resized[1].max()
    if fmax > fmin:
        resized[1] = (resized[1] - fmin) / (fmax - fmin)
    else:
        resized[1] = 0.0

    # Re-threshold mask after resize
    resized[2] = (resized[2] > 0.5).astype(np.float32)

    area_px = int(cell_pixels.sum())
    return resized, area_px


# ── Process a single TIFF ───────────────────────────────────────────────────

def _process_fov(args):
    """Process a single FOV — designed for use with multiprocessing.Pool."""
    fov_data, fov_idx, mask_ch, phase_ch, fluor_ch, condition_meta, tiff_name = args

    # Convert binary mask → labeled mask via connected components
    binary_mask = fov_data[mask_ch] > 0
    labeled_mask = sk_label(binary_mask, connectivity=1)
    n_cells = int(labeled_mask.max())

    if n_cells == 0:
        return [], [], n_cells

    crops = []
    metadata = []
    cell_labels = set(np.unique(labeled_mask)) - {0}

    for cell_label in sorted(cell_labels):
        crop, area_px = extract_cell_crop(
            fov_data, labeled_mask, cell_label,
            phase_ch=phase_ch, fluor_ch=fluor_ch,
        )
        if crop is None:
            continue

        crops.append(crop)
        metadata.append({
            "condition_label": condition_meta["condition_label"],
            "condition_type": condition_meta["condition_type"],
            "gene": condition_meta["gene"],
            "drug": condition_meta["drug"],
            "concentration": condition_meta["concentration"],
            "replica": condition_meta["replica"],
            "tiff_file": tiff_name,
            "fov_index": fov_idx,
            "cell_label": int(cell_label),
            "area_px": area_px,
            "is_control": condition_meta["is_control"],
            "is_drug": condition_meta["is_drug"],
        })

    return crops, metadata, n_cells


def process_single_tiff(tiff_path, condition_meta, n_workers=None):
    """
    Load a segmented TIFF and extract 128×128 crops for all cells.

    The masks are already quality-filtered from the segmentation pipeline,
    so no additional QC is applied here. FOVs are processed in parallel.

    Returns:
        crops:    list of (3, 128, 128) float32 arrays
        metadata: list of dicts with per-cell info
        stats:    dict with processing statistics
    """
    from multiprocessing import Pool, cpu_count

    data, tiff_meta = load_hyperstack(tiff_path)
    n_fov, n_ch = data.shape[0], data.shape[1]

    # Mask is always the last channel
    mask_ch = n_ch - 1

    # Identify phase and fluorescence channels among non-mask channels.
    # Phase contrast has a bright background (high median); fluorescence is
    # mostly dark with sparse bright spots (low median).
    candidate_chs = [c for c in range(n_ch) if c != mask_ch]
    if len(candidate_chs) == 2:
        medians = {c: np.median(data[:, c]) for c in candidate_chs}
        phase_ch = max(medians, key=medians.get)
        fluor_ch = min(medians, key=medians.get)
    elif len(candidate_chs) == 1:
        phase_ch = candidate_chs[0]
        fluor_ch = candidate_chs[0]  # fallback: duplicate phase as fluor
        print("  WARNING: only 1 non-mask channel found, using it for both phase and fluorescence")
    else:
        raise ValueError(f"Expected 2 non-mask channels, got {len(candidate_chs)}")

    print(f"    Phase channel: {phase_ch}, Fluorescence channel: {fluor_ch}  "
          f"(of {n_ch} channels)")

    if n_workers is None:
        n_workers = min(cpu_count(), n_fov)

    # Build args for each FOV
    fov_args = [
        (data[fov_idx], fov_idx, mask_ch, phase_ch, fluor_ch,
         condition_meta, tiff_path.name)
        for fov_idx in range(n_fov)
    ]

    crops = []
    metadata = []
    total_cells = 0

    if n_workers > 1 and n_fov > 1:
        with Pool(n_workers) as pool:
            for fov_crops, fov_meta, n_cells in pool.imap_unordered(
                    _process_fov, fov_args):
                crops.extend(fov_crops)
                metadata.extend(fov_meta)
                total_cells += n_cells
    else:
        for args in fov_args:
            fov_crops, fov_meta, n_cells = _process_fov(args)
            crops.extend(fov_crops)
            metadata.extend(fov_meta)
            total_cells += n_cells

    del data  # free TIFF memory immediately

    stats = {"n_fov": n_fov, "total_cells": total_cells}
    return crops, metadata, stats


# ── HDF5 storage ────────────────────────────────────────────────────────────

def _write_metadata_to_h5(h5_file, all_metadata, n):
    """
    Write metadata datasets and attributes to an already-open HDF5 file.
    The 'crops' dataset must already exist and be populated.
    """
    print(f"\nWriting metadata for {n:,} cells ...")

    # Build globally unique FOV IDs
    fov_key_to_id = {}
    fov_ids = np.empty(n, dtype=np.int32)
    for i, meta in enumerate(all_metadata):
        key = (meta["tiff_file"], meta["fov_index"])
        if key not in fov_key_to_id:
            fov_key_to_id[key] = len(fov_key_to_id)
        fov_ids[i] = fov_key_to_id[key]

    str_dt = h5py.string_dtype()

    h5_file.create_dataset("condition_labels",
                           data=[m["condition_label"] for m in all_metadata],
                           dtype=str_dt)
    h5_file.create_dataset("condition_types",
                           data=[m["condition_type"] for m in all_metadata],
                           dtype=str_dt)
    h5_file.create_dataset("genes",
                           data=[m["gene"] for m in all_metadata],
                           dtype=str_dt)
    h5_file.create_dataset("drugs",
                           data=[m["drug"] for m in all_metadata],
                           dtype=str_dt)
    h5_file.create_dataset("concentrations",
                           data=[m["concentration"] for m in all_metadata],
                           dtype=str_dt)
    h5_file.create_dataset("replicas",
                           data=[m["replica"] for m in all_metadata],
                           dtype=str_dt)
    h5_file.create_dataset("tiff_files",
                           data=[m["tiff_file"] for m in all_metadata],
                           dtype=str_dt)
    h5_file.create_dataset("fov_ids", data=fov_ids)
    h5_file.create_dataset("is_control",
                           data=[m["is_control"] for m in all_metadata],
                           dtype=bool)
    h5_file.create_dataset("is_drug",
                           data=[m["is_drug"] for m in all_metadata],
                           dtype=bool)
    h5_file.create_dataset("areas_px",
                           data=[m["area_px"] for m in all_metadata],
                           dtype=np.int32)

    # Attributes
    h5_file.attrs["crop_size"] = CROP_SIZE
    h5_file.attrs["n_channels"] = N_CHANNELS
    h5_file.attrs["channel_names"] = ["phase", "parb", "mask"]
    h5_file.attrs["total_cells"] = n
    h5_file.attrs["n_fovs"] = len(fov_key_to_id)


# ── Summary ─────────────────────────────────────────────────────────────────

def print_summary(all_metadata):
    """Print per-condition cell counts."""
    counts = defaultdict(int)
    drug_flags = {}
    ctrl_flags = {}
    for m in all_metadata:
        cl = m["condition_label"]
        counts[cl] += 1
        drug_flags[cl] = m["is_drug"]
        ctrl_flags[cl] = m["is_control"]

    print("\n" + "=" * 60)
    print(f"{'Condition':<20} {'Type':<10} {'Cells':>8}")
    print("-" * 60)

    for cl in sorted(counts, key=lambda x: (-counts[x])):
        ctype = "control" if ctrl_flags[cl] else ("drug" if drug_flags[cl] else "mutant")
        print(f"{cl:<20} {ctype:<10} {counts[cl]:>8,}")

    print("-" * 60)
    n_ctrl = sum(c for cl, c in counts.items() if ctrl_flags[cl])
    n_drug = sum(c for cl, c in counts.items() if drug_flags[cl])
    n_mutant = sum(c for cl, c in counts.items()
                   if not ctrl_flags[cl] and not drug_flags[cl])
    print(f"{'TOTAL':<20} {'':10} {sum(counts.values()):>8,}")
    print(f"  Controls: {n_ctrl:,}  |  Drugs: {n_drug:,}  |  Mutants: {n_mutant:,}")
    print(f"  Conditions: {len(counts)} "
          f"({sum(1 for v in drug_flags.values() if v)} drugs, "
          f"{sum(1 for cl in counts if not ctrl_flags[cl] and not drug_flags[cl])} mutants, "
          f"{sum(1 for v in ctrl_flags.values() if v)} controls)")

    # Flag conditions with few cells
    min_cells = CONFIG["min_cells_per_condition"]
    low = {cl: c for cl, c in counts.items() if c < min_cells}
    if low:
        print(f"\n  WARNING: {len(low)} condition(s) below {min_cells} cells:")
        for cl, c in sorted(low.items(), key=lambda x: x[1]):
            print(f"    {cl}: {c}")
    print("=" * 60)


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Extract single-cell crops from segmented TIFFs")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print stats only, don't write HDF5")
    parser.add_argument("--tiff-dir", type=str, action="append", default=None,
                        help="Override tiff_dirs from config (can repeat)")
    parser.add_argument("--tiff-file", type=str, action="append", default=None,
                        help="Process specific TIFF file(s) instead of dirs")
    parser.add_argument("--workers", type=int, default=None,
                        help="Parallel workers per TIFF (default: cpu_count)")
    args = parser.parse_args()

    # Determine which TIFFs to process
    if args.tiff_file:
        # Specific files provided
        tiff_files = []
        for fpath in args.tiff_file:
            p = Path(fpath)
            if not p.exists():
                print(f"WARNING: file not found: {p}")
                continue
            try:
                meta = parse_condition_from_filename(p.stem)
            except ValueError as e:
                print(f"  SKIP: {e}")
                continue
            tiff_files.append((p, meta))
        print(f"Processing {len(tiff_files)} specified TIFF file(s)")
    elif args.tiff_dir:
        # Override directories
        orig = CONFIG["tiff_dirs"]
        CONFIG["tiff_dirs"] = args.tiff_dir
        tiff_files = discover_tiff_files()
        CONFIG["tiff_dirs"] = orig  # restore
    else:
        tiff_files = discover_tiff_files()

    if not tiff_files:
        print("No TIFFs found. Check tiff_dirs in config_supcon.py.")
        sys.exit(1)

    all_metadata = []
    total_cells = 0

    output_path = Path(CONFIG["crop_output_h5"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    h5_file = None
    crops_ds = None
    n_written = 0

    if not args.dry_run:
        h5_file = h5py.File(str(output_path), "w")
        # Resizable dataset — start empty, grow as crops arrive
        crops_ds = h5_file.create_dataset(
            "crops",
            shape=(0, N_CHANNELS, CROP_SIZE, CROP_SIZE),
            maxshape=(None, N_CHANNELS, CROP_SIZE, CROP_SIZE),
            dtype=np.float32,
            chunks=(256, N_CHANNELS, CROP_SIZE, CROP_SIZE),
            compression="gzip", compression_opts=4,
        )

    t0 = time.time()

    try:
        for i, (tiff_path, condition_meta) in enumerate(tiff_files):
            print(f"\n[{i + 1}/{len(tiff_files)}] {tiff_path.name}")
            crops, metadata, stats = process_single_tiff(tiff_path, condition_meta,
                                                           n_workers=args.workers)
            total_cells += stats["total_cells"]
            print(f"  FOVs: {stats['n_fov']}  |  Cells: {stats['total_cells']:,}")

            all_metadata.extend(metadata)

            # Stream crops to HDF5 immediately, then free memory
            if crops_ds is not None and crops:
                batch = np.stack(crops, axis=0)
                new_size = n_written + len(crops)
                crops_ds.resize(new_size, axis=0)
                crops_ds[n_written:new_size] = batch
                n_written = new_size
                del batch
            del crops  # free crop memory

        elapsed = time.time() - t0
        print(f"\n{'=' * 60}")
        print(f"Processed {len(tiff_files)} TIFFs in {elapsed:.0f}s")
        print(f"Total cells: {total_cells:,}")

        print_summary(all_metadata)

        if h5_file is not None and n_written > 0:
            _write_metadata_to_h5(h5_file, all_metadata, n_written)
            print(f"\nDone. {n_written:,} crops saved to {output_path}")
        elif args.dry_run:
            print("\n(Dry run — no HDF5 written)")
    finally:
        if h5_file is not None:
            h5_file.close()


if __name__ == "__main__":
    main()
