"""
Mapping stemness- and EMT-associated gene expression in HER2+ breast cancer
at single-cell resolution.

This script reproduces every processed dataset, figure and table in the
project from the raw GEO files (GSE176078, Wu et al. 2021, Nature Genetics).

What it does, in order:
  1. Screens all five HER2+ patients in GSE176078, using metadata only, and
     records which ones have annotated malignant epithelial cells.
  2. For each analysed patient (CID3921, CID4066, CID45171), one at a time:
       - builds a cells x genes AnnData object from the raw count matrix
       - calculates QC metrics and applies exploratory QC thresholds
       - saves the raw-count object (before any transformation)
       - normalises to 10,000 counts per cell and log1p-transforms ONCE,
         keeping the raw counts in layers["counts"]
       - saves the normalised object
       - keeps malignant epithelial cells in the two states shared by all
         three patients (Cancer Her2 SC and Cancer Cycling)
       - summarises the 12 candidate genes per population and saves
         four tables and a dot plot
  3. Combines the three patients into one figure on a shared colour scale,
     and one long-format summary table.
  4. Saves sequencing-depth summaries, which are needed to interpret
     detection rates.

How to run (Terminal, not the Python interpreter):
    conda activate brast-ml
    python /Volumes/Expansion/BRaST-ML/breast-cancer-stemness-ml/scripts/run_analysis.py

Expected input layout (each GEO archive unpacked with `tar -xzf`):
    data/raw/<PATIENT>/count_matrix_sparse.mtx
    data/raw/<PATIENT>/count_matrix_genes.tsv
    data/raw/<PATIENT>/count_matrix_barcodes.tsv
    data/raw/<PATIENT>/metadata.csv

Scientific limitations (also in README.md):
  - Expression of individual markers does not establish cancer stem cell
    identity, and EMT cannot be classified from a few genes.
  - scRNA-seq dropout means zero detection is not proof of absence, and
    detection rates depend on sequencing depth.
  - Cells are not independent replicates; the patient is the unit of
    replication, and n = 3 patients allows description, not inference.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # draw figures to files, without opening windows

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

# Paths are worked out from this file's location, so the script behaves the
# same whatever folder Terminal is in when it is run.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
FIGURES_DIR = PROJECT_ROOT / "results" / "figures"
TABLES_DIR = PROJECT_ROOT / "results" / "tables"

# All five clinically HER2+ tumours in GSE176078.
SCREENED_PATIENTS = ["CID3586", "CID3838", "CID3921", "CID4066", "CID45171"]

# HER2+ tumours with annotated malignant epithelial cells (from screening).
ANALYSED_PATIENTS = ["CID3921", "CID4066", "CID45171"]

# Malignant epithelial states present in all three analysed patients.
# Smaller labels (Cancer LumA/LumB/Basal SC, at most 9 cells per patient) are
# excluded because there are too few cells to summarise.
POPULATIONS = ["Cancer Her2 SC", "Cancer Cycling"]

STEMNESS_GENES = ["SOX2", "POU5F1", "NANOG", "ALDH1A1", "PROM1", "CD44"]
EMT_GENES = ["VIM", "CDH1", "CDH2", "ZEB1", "SNAI1", "FN1"]
CANDIDATE_GENES = STEMNESS_GENES + EMT_GENES

# Exploratory QC thresholds: remove likely empty droplets / low-quality cells
# (too few genes) and damaged or dying cells (high mitochondrial fraction).
MIN_GENES = 200
MAX_PCT_MT = 20

# Library-size normalisation target (counts per cell) before log1p.
TARGET_SUM = 10_000

FIGURE_DPI = 300


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def check_output_folders():
    """Create output folders if needed, and confirm that the raw data folder exists."""
    if not RAW_DIR.is_dir():
        raise FileNotFoundError(f"Raw data folder not found: {RAW_DIR}")
    for folder in (PROCESSED_DIR, FIGURES_DIR, TABLES_DIR):
        folder.mkdir(parents=True, exist_ok=True)


def save_current_figure(filename):
    """Save the active matplotlib figure into results/figures and close it."""
    path = FIGURES_DIR / filename
    plt.savefig(path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close("all")
    print(f"    saved figure: {path.relative_to(PROJECT_ROOT)}")


def save_table(df, filename, index=True):
    """Save a DataFrame as CSV into results/tables."""
    path = TABLES_DIR / filename
    df.to_csv(path, index=index)
    print(f"    saved table:  {path.relative_to(PROJECT_ROOT)}")


def screen_her2_patients():
    """Summarise subtype and malignant cell counts for every HER2+ patient.

    Only metadata.csv is read here, so this step is fast and uses little memory.
    Patients whose folders are missing are reported and skipped.
    """
    print("\n[1] Screening HER2+ patients (metadata only)")
    rows = []
    for patient in SCREENED_PATIENTS:
        meta_path = RAW_DIR / patient / "metadata.csv"
        if not meta_path.exists():
            print(f"    {patient}: metadata.csv not found, skipped")
            continue
        meta = pd.read_csv(meta_path, index_col=0)
        malignant = meta[meta["celltype_major"] == "Cancer Epithelial"]
        minor_counts = malignant["celltype_minor"].value_counts()
        other_labels = minor_counts.drop(POPULATIONS, errors="ignore")
        rows.append({
            "patient": patient,
            "clinical_subtype": ";".join(sorted(meta["subtype"].astype(str).unique())),
            "total_cells": len(meta),
            "cancer_epithelial_cells": len(malignant),
            "cancer_her2_sc_cells": int(minor_counts.get("Cancer Her2 SC", 0)),
            "cancer_cycling_cells": int(minor_counts.get("Cancer Cycling", 0)),
            "other_cancer_labels": "; ".join(f"{k} ({v})" for k, v in other_labels.items()),
            "included_in_analysis": patient in ANALYSED_PATIENTS,
        })
    screening = pd.DataFrame(rows)
    print(screening.to_string(index=False))
    save_table(screening, "HER2_patient_screening.csv", index=False)


def load_raw_counts(patient):
    """Build a cells x genes AnnData object of raw counts for one patient.

    The GEO matrix is stored genes x cells, so it is transposed.
    """
    folder = RAW_DIR / patient
    matrix = mmread(folder / "count_matrix_sparse.mtx").tocsr()  # genes x cells
    genes = pd.read_csv(folder / "count_matrix_genes.tsv", header=None)[0].astype(str)
    barcodes = pd.read_csv(folder / "count_matrix_barcodes.tsv", header=None)[0].astype(str)
    metadata = pd.read_csv(folder / "metadata.csv", index_col=0)

    # Safety checks: the matrix, gene list, barcode list and metadata must agree.
    assert matrix.shape == (len(genes), len(barcodes)), "matrix shape does not match genes/barcodes"
    assert (barcodes.values == metadata.index.values).all(), "barcode order differs from metadata"

    adata = sc.AnnData(
        X=matrix.T.tocsr(),                      # cells x genes
        obs=metadata.copy(),
        var=pd.DataFrame(index=genes.values),
    )
    adata.obs_names = barcodes.values
    adata.var.index.name = None
    adata.var_names_make_unique()
    print(f"    loaded {adata.n_obs} cells x {adata.n_vars} genes (max raw count {adata.X.max():.0f})")
    return adata


def run_qc(adata, patient):
    """Calculate QC metrics, save a QC violin plot and apply thresholds."""
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True)

    # Our recalculated totals should equal the authors' nCount_RNA.
    assert np.allclose(adata.obs["total_counts"], adata.obs["nCount_RNA"]), "QC totals differ from metadata"

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4, multi_panel=True, show=False,
    )
    save_current_figure(f"{patient}_QC_violin.png")

    keep = (adata.obs["n_genes_by_counts"] >= MIN_GENES) & (adata.obs["pct_counts_mt"] < MAX_PCT_MT)
    n_removed = int((~keep).sum())
    print(f"    {int(adata.var['mt'].sum())} mitochondrial genes; "
          f"{int(keep.sum())}/{adata.n_obs} cells pass QC ({n_removed} removed)")
    return adata[keep].copy()


def normalise(adata):
    """Library-size normalise and log1p-transform once, keeping raw counts."""
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    sc.pp.log1p(adata)
    print(f"    normalised: max log-expression {adata.X.max():.3f}")
    return adata


def malignant_subset(adata):
    """Keep malignant epithelial cells in the two shared populations."""
    is_malignant = adata.obs["celltype_major"] == "Cancer Epithelial"
    in_populations = adata.obs["celltype_minor"].isin(POPULATIONS)
    cancer = adata[is_malignant & in_populations].copy()
    # Plain-text copy of the labels, so unused categories never appear.
    cancer.obs["population"] = cancer.obs["celltype_minor"].astype(str)
    counts = cancer.obs["population"].value_counts().to_dict()
    print(f"    malignant subset: {cancer.n_obs} cells {counts}")
    return cancer


def summarise_genes(cancer, patient):
    """Per-population summaries of the candidate genes.

    detection %      = % of cells with expression > 0 (how widespread)
    mean             = mean log-normalised expression over ALL cells, zeros
                       included (average signal)
    mean in positive = mean over expressing cells only (how strong, where
                       detected). Interpret only alongside n positive.
    """
    expr = pd.DataFrame(
        cancer[:, CANDIDATE_GENES].X.toarray(),
        index=cancer.obs_names,
        columns=CANDIDATE_GENES,
    )
    groups = pd.Series(cancer.obs["population"].values, index=expr.index, name="celltype_minor")

    detected = expr > 0
    detection = detected.groupby(groups).mean() * 100
    means = expr.groupby(groups).mean()
    n_positive = detected.groupby(groups).sum()
    positive_means = expr.where(detected).groupby(groups).mean()

    # Internal consistency check: mean = mean in positive cells x fraction positive.
    assert np.allclose((positive_means.fillna(0) * detection / 100).values, means.values)

    save_table(detection, f"{patient}_gene_detection_percentages.csv")
    save_table(means, f"{patient}_mean_gene_expression.csv")
    save_table(positive_means, f"{patient}_mean_expression_positive_cells.csv")
    save_table(n_positive, f"{patient}_n_positive_cells.csv")

    # Long-format rows for the combined table.
    n_cells = groups.value_counts()
    long_rows = []
    for population in detection.index:
        for gene in CANDIDATE_GENES:
            long_rows.append({
                "patient": patient,
                "population": population,
                "gene": gene,
                "gene_set": "Stemness" if gene in STEMNESS_GENES else "EMT",
                "n_cells": int(n_cells[population]),
                "n_positive": int(n_positive.loc[population, gene]),
                "pct_detected": detection.loc[population, gene],
                "mean_expression": means.loc[population, gene],
                "mean_in_positive_cells": positive_means.loc[population, gene],
            })
    return pd.DataFrame(long_rows)


def depth_summary(cancer, patient):
    """Median sequencing depth per population (needed to interpret detection rates)."""
    depth = (
        cancer.obs.groupby("population")[["nCount_RNA", "nFeature_RNA"]]
        .median()
        .rename(columns={"nCount_RNA": "median_counts_per_cell", "nFeature_RNA": "median_genes_per_cell"})
        .reset_index()
    )
    depth.insert(0, "patient", patient)
    depth.insert(2, "n_cells", depth["population"].map(cancer.obs["population"].value_counts()).values)
    return depth


def process_patient(patient):
    """Run the full per-patient workflow and return a small 12-gene subset."""
    print(f"\n[2] Processing {patient}")
    adata = load_raw_counts(patient)
    adata = run_qc(adata, patient)

    # Save raw counts BEFORE any transformation.
    raw_path = PROCESSED_DIR / f"{patient}.h5ad"
    adata.write_h5ad(raw_path)
    print(f"    saved raw object: {raw_path.relative_to(PROJECT_ROOT)}")

    adata = normalise(adata)
    norm_path = PROCESSED_DIR / f"{patient}_normalised.h5ad"
    adata.write_h5ad(norm_path)
    print(f"    saved normalised object: {norm_path.relative_to(PROJECT_ROOT)}")

    cancer = malignant_subset(adata)
    del adata  # free memory before the next steps

    summary = summarise_genes(cancer, patient)
    depth = depth_summary(cancer, patient)

    sc.pl.dotplot(cancer, var_names=CANDIDATE_GENES, groupby="population", standard_scale=None, show=False)
    save_current_figure(f"{patient}_Cancer_Epithelial_Stemness_EMT_DotPlot.png")

    if patient == "CID3921":
        # Exploratory distribution plot for the four most informative genes.
        sc.pl.violin(cancer, keys=["CD44", "PROM1", "VIM", "CDH1"], groupby="population", rotation=90, show=False)
        save_current_figure("CID3921_Cancer_Epithelial_Violin_CD44_PROM1_VIM_CDH1.png")

    # Keep only the 12 candidate genes for the combined figure (small in memory).
    return cancer[:, CANDIDATE_GENES].copy(), summary, depth


def combined_outputs(subsets, summaries, depths):
    """Three-patient figure on one shared colour scale, plus combined tables."""
    print("\n[3] Combining patients")
    combined = sc.concat(subsets)
    combined.obs["group"] = (
        combined.obs["orig.ident"].astype(str) + " | " + combined.obs["population"].astype(str)
    )
    print(f"    combined: {combined.n_obs} cells x {combined.n_vars} genes")

    sc.pl.dotplot(
        combined,
        var_names={"Stemness": STEMNESS_GENES, "EMT": EMT_GENES},
        groupby="group",
        standard_scale=None,
        show=False,
    )
    save_current_figure("AllPatients_HER2_Cancer_Epithelial_Stemness_EMT_DotPlot.png")

    save_table(pd.concat(summaries, ignore_index=True), "AllPatients_gene_expression_summary.csv", index=False)
    save_table(pd.concat(depths, ignore_index=True), "AllPatients_sequencing_depth.csv", index=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    sc.settings.verbosity = 1  # warnings only
    print(f"Project root: {PROJECT_ROOT}")
    check_output_folders()

    screen_her2_patients()

    subsets, summaries, depths = [], [], []
    for patient in ANALYSED_PATIENTS:  # one patient at a time to limit memory use
        subset, summary, depth = process_patient(patient)
        subsets.append(subset)
        summaries.append(summary)
        depths.append(depth)

    combined_outputs(subsets, summaries, depths)
    print("\nDone. All outputs are in data/processed/, results/figures/ and results/tables/.")


if __name__ == "__main__":
    main()
