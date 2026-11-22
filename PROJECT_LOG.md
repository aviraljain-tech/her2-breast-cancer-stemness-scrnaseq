# PROJECT LOG — Mapping Stemness-Associated Gene Expression in Breast Cancer at Single-Cell Resolution

**Author:** Aviral Jain (GitHub: `aviraljain-tech`)
**Project root:** `/Volumes/Expansion/BRaST-ML/breast-cancer-stemness-ml`
**Log started:** 2026-11-22
**Last updated:** 2026-11-22 (reproducible script + README)

How to read this log:
- **[DONE]** = completed and verified (file exists on disk, or output was checked).
- **[DONE – reported]** = completed in an earlier session; recorded from notes, not re-checked since.
- **[TO CHECK]** = something that needs verifying before relying on it.
- **[PLANNED]** = not yet started.

---

## 1. Objective and research questions

**Objective:** Describe the expression of selected stemness-associated and epithelial-to-mesenchymal transition (EMT)-associated genes in annotated malignant epithelial cells from public breast cancer scRNA-seq data.

**Research questions**
1. Which candidate stemness- and EMT-associated genes are detectable in malignant epithelial cells?
2. Do detection rates and mean expression differ between annotated malignant epithelial subpopulations within a sample?
3. (Planned) Are any within-sample patterns consistent across patients?

**Biological motivation.** Breast tumours are heterogeneous. Cancer stem cell-like properties (self-renewal, tumour initiation) are linked to progression, therapy resistance and recurrence. EMT programmes overlap with, but are not the same as, stemness programmes.

**Scientific limitations (apply to all results below)**
- Expression of individual markers does not establish cancer stem cell identity.
- EMT is continuous and multi-gene; it cannot be classified from one or two genes.
- scRNA-seq has substantial dropout: zero detection ≠ biological absence.
- Cells from one patient are not independent biological replicates.
- Within-patient comparisons are exploratory observations, not generalisable conclusions.

---

## 2. Computing environment

| Item | Value |
|---|---|
| Machine | MacBook Air M2, 8 GB RAM, macOS 26.5 |
| Conda | Anaconda at `/opt/anaconda3`; project env `brast-ml` (base env Python 3.13.5, not modified) |
| Python | 3.12.16 |
| Scanpy | 1.12.5 |
| AnnData | 0.13.4 |
| NumPy | 2.5.5 |
| scikit-learn | 1.9.1 |
| Storage | External SSD (exFAT). macOS creates `._*` sidecar files on exFAT; these are excluded by `.gitignore`. |

Git: initialised on branch `main`, one initial commit (2026-11-21). Not yet pushed to GitHub. `.gitignore` excludes `data/`, `*.h5ad`, `*.h5`, `*.mtx`, compressed matrices, Python caches, environments, `.DS_Store` and `._*`.

---

## 3. Dataset

**Source:** GEO accession **GSE176078** — Wu et al. 2021, *A single-cell and spatially resolved atlas of human breast cancers* (Nature Genetics).

Files per sample: `count_matrix_sparse.mtx` (genes × cells), `count_matrix_genes.tsv`, `count_matrix_barcodes.tsv`, `metadata.csv`.

Metadata columns: `orig.ident`, `nCount_RNA`, `nFeature_RNA`, `percent.mito`, `subtype`, `celltype_subset`, `celltype_minor`, `celltype_major`.

### On disk [DONE – verified 2026-11-22]
```
data/raw/CID3921.tar.gz            (17.4 MB)
data/raw/CID3921/                  (4 files above)
data/raw/CID3946.tar.gz            (1.6 MB)
data/raw/CID3946/                  (4 files above)
```

### CID3946 — pilot, excluded [DONE – reported]
774 cells × 29,733 genes. Populations: PVL 248, CAFs 167, Myeloid 157, Endothelial 110, T-cells 92.
**Decision:** excluded from the malignant epithelial analysis because it contains no annotated epithelial cells.

### CID3921 — main sample [DONE – reported]
HER2-positive. 3,024 cells × 29,733 genes.

| Major population | Cells |
|---|---:|
| T-cells | 1,473 |
| Cancer Epithelial | 441 |
| Myeloid | 385 |
| Endothelial | 210 |
| Plasmablasts | 175 |
| B-cells | 162 |
| CAFs | 106 |
| PVL | 72 |

| Cancer epithelial subpopulation (`celltype_minor`) | Cells |
|---|---:|
| Cancer Her2 SC | 377 |
| Cancer Cycling | 64 |
| **Total** | **441** |

Metadata barcode order was checked and matches the count-matrix barcode order exactly.

### Additional HER2+ patients: screening (2026-11-22) [DONE – verified from metadata.csv]
**Decision:** add only clinically HER2+ patients. Cancer-cell states in this dataset are labelled by subtype-like programmes (e.g. "Cancer Her2 SC"), so mixing clinical subtypes would confound the comparison. Wu et al. report 5 HER2+ tumours.

Downloaded from GEO (saved to `data/raw/` and unpacked):
- `GSM5354513_CID3586.tar.gz`
- `GSM5354514_CID3838.tar.gz`
- `GSM5354521_CID4066.tar.gz`
- `GSM5354535_CID45171.tar.gz`

The `subtype` column in each `metadata.csv` confirmed HER2+. Together with CID3921, these are all 5 HER2+ patients.

| Patient | Cells | Cancer Epithelial | Cancer Her2 SC | Cancer Cycling | Other cancer labels | Use? |
|---|---:|---:|---:|---:|---|---|
| CID3921 | 3,024 | 441 | 377 | 64 | – | Yes (done) |
| CID3586 | 6,178 | 0 | 0 | 0 | 698 *Normal* Epithelial only | **No**: no annotated cancer epithelial cells |
| CID3838 | 2,353 | 0 | 0 | 0 | no epithelial cells | **No**: no cancer epithelial cells |
| CID4066 | 5,309 | 521 | 393 | 112 | LumA SC 8, LumB SC 6, Basal SC 2 | Yes |
| CID45171 | 2,447 | 813 | 567 | 236 | LumB SC 9, Basal SC 1 | Yes |

**Decisions:**
- The multi-patient analysis uses **3 patients** (CID3921, CID4066, CID45171).
- Comparisons are restricted to the two states present in all three: **Cancer Her2 SC** and **Cancer Cycling**.
- The small LumA/LumB/Basal SC groups (≤ 9 cells each) are excluded from comparisons: too few cells.

**Observation:** in all five samples, the minimum detected genes per cell is ≥ 201 and the maximum mitochondrial % is < 20. This is *consistent with* the authors having applied similar filters before deposition, but it does not prove it.

CID4066 checks:
- 29,733 genes with no duplicate gene names.
- Barcode order in `count_matrix_barcodes.tsv` matches `metadata.csv`.

---

## 4. Preprocessing

### 4.1 AnnData construction [DONE – reported]
The matrix was transposed from genes × cells to cells × genes:
```python
adata = sc.AnnData(X=matrix.T, obs=metadata.copy(), var=pd.DataFrame(index=genes.astype(str)))
adata.obs_names = barcodes.astype(str)
adata.var.index.name = None
adata.var_names_make_unique()
```
Saved as `data/processed/CID3921.h5ad` (intended to hold **raw counts**).

### 4.2 Quality control [DONE – reported]
```python
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True)
```
- 13 mitochondrial genes identified.
- Recalculated metrics matched the original metadata exactly.

| Metric | Mean | Median | Min | Max |
|---|---:|---:|---:|---:|
| Total counts | 8796.40 | 3301 | 263 | 117276 |
| Detected genes | 1822.38 | 1156.5 | 202 | 7977 |
| Mitochondrial % | 5.92 | 5.23 | 0 | 19.82 |

**Thresholds chosen:** ≥ 200 detected genes; < 20% mitochondrial counts.
**Why:** common exploratory thresholds for removing empty droplets/low-quality cells (low gene count) and damaged or dying cells (high mitochondrial fraction).
**Outcome:** all 3,024 cells passed; none removed. This shows only that all cells passed *our* thresholds; it does not prove the deposited data were pre-filtered.

### 4.3 Normalisation [DONE – reported]
Run once:
```python
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=10000)
sc.pp.log1p(adata)
```
- `adata.X` = library-size-normalised (10,000 counts per cell), log1p-transformed.
- `adata.layers["counts"]` = original raw counts.
- Saved as `data/processed/CID3921_normalised.h5ad`.
- **Rule:** never normalise or log-transform this object again. To restart, load the raw object.

### 4.4 Processed files on disk [DONE – verified 2026-11-22]
| File | Contents | Size | Last modified |
|---|---|---:|---|
| `data/processed/CID3921.h5ad` | **Raw counts** in `X`; no `counts` layer; obs/var include QC metric columns | 69,432,086 B | 2026-11-22 15:06 (rebuilt) |
| `data/processed/CID3921_normalised.h5ad` | log1p-normalised `X` + `layers["counts"]` | 113,581,882 B | 2026-11-22 14:34 |
| `data/processed/CID3921_cancer_epithelial.h5ad` | 441-cell cancer subset (contents not yet re-checked) | 41,965,806 B | 2026-11-22 14:34 |

#### Incident and fix: raw file overwritten (2026-11-22) [RESOLVED]
- **Problem:** `CID3921.h5ad` had been overwritten with the normalised object. Evidence: it had the same byte size as the normalised file, and on reload `X.max()` = 9.000352 (log scale) with a `counts` layer present.
- **Verification that raw counts survived:** `adata.layers["counts"]` has max 14133, and its per-cell row sums equal the original `nCount_RNA` (`np.allclose` → `True`).
- **Side issue:** in AnnData 0.13.4, `adata.layers[None]` is an alias for `adata.X` (`raw.X is raw.layers[None]` → `True`). It is not a duplicate layer.
  - Running `del adata.layers[None]` by mistake deleted `X` from `adata` and `cancer` in memory.
  - `X` was restored from `CID3921_normalised.h5ad` after checking that cell order, gene order and the `counts` layer were identical (`nnz` of differences = 0). This restored saved values; it was not a second normalisation.
  - The restored `cancer` values match `cancer_expression` (`np.allclose` → `True`).
  - **Rule: never run `del adata.layers[None]`.**
- **Fix:**
  - The raw object was rebuilt as a copy of `adata`, with `X` = `layers["counts"]`, the `counts` layer removed, and `uns["log1p"]` removed.
  - It was saved to `CID3921.h5ad` and reloaded to verify: shape (3024, 29733), `X.max()` = 14133, layers `[None]` (X only), and row sums equal `nCount_RNA` (`True`).
  - The rebuilt file also contains the QC columns added by `calculate_qc_metrics`. This does not affect the raw counts.
---

## 5. Candidate genes [DONE – reported]
```python
stemness_genes = ["SOX2", "POU5F1", "NANOG", "ALDH1A1", "PROM1", "CD44"]
emt_genes      = ["VIM", "CDH1", "CDH2", "ZEB1", "SNAI1", "FN1"]
candidate_genes = stemness_genes + emt_genes
```
All 12 present in `adata.var_names` (`present_genes` = all 12).

Note: CDH1 (E-cadherin) is an *epithelial* marker. It is included in the EMT panel because its loss is a hallmark of EMT, not because its expression indicates EMT.

---

## 6. Analyses and results

### 6.1 Candidate genes across all 3,024 cells (preliminary) [DONE – reported]
Includes immune and stromal cells, so the signal may come from non-malignant cells. Treated as preliminary only.

| Gene | Cells expressing (>0) | Mean log-norm expression |
|---|---:|---:|
| VIM | 2440 | 2.557 |
| CD44 | 1975 | 1.147 |
| FN1 | 580 | 0.416 |
| CDH1 | 351 | 0.096 |
| SNAI1 | 320 | 0.098 |
| ZEB1 | 315 | 0.134 |
| PROM1 | 290 | 0.063 |
| ALDH1A1 | 151 | 0.038 |
| CDH2 | 17 | 0.003 |
| POU5F1 | 13 | 0.002 |
| SOX2 | 5 | 0.000 |
| NANOG | 1 | 0.000 |

### 6.2 Malignant epithelial subset (n = 441) [DONE – verified via saved tables]
```python
cancer = adata[adata.obs["celltype_major"] == "Cancer Epithelial"].copy()
```
Metrics per `celltype_minor`:
- **Detection %** = percentage of cells with expression > 0 (how widespread).
- **Mean expression** = mean of log-normalised values over *all* cells in the group, zeros included (average signal).

| Gene | Set | Detection % Cycling (n=64) | Detection % Her2 SC (n=377) | Mean Cycling | Mean Her2 SC |
|---|---|---:|---:|---:|---:|
| SOX2 | Stemness | 1.56 | 0.80 | 0.005 | 0.002 |
| POU5F1 | Stemness | 1.56 | 2.12 | 0.003 | 0.005 |
| NANOG | Stemness | 0.00 | 0.00 | 0.000 | 0.000 |
| ALDH1A1 | Stemness | 6.25 | 8.75 | 0.035 | 0.050 |
| PROM1 | Stemness | 59.38 | 66.58 | 0.399 | 0.433 |
| CD44 | Stemness | 60.94 | 69.23 | 0.457 | 0.532 |
| VIM | EMT | 15.62 | 28.38 | 0.076 | 0.178 |
| CDH1 | EMT (epithelial) | 65.62 | 79.58 | 0.545 | 0.651 |
| CDH2 | EMT | 1.56 | 2.12 | 0.004 | 0.007 |
| ZEB1 | EMT | 0.00 | 0.53 | 0.000 | 0.002 |
| SNAI1 | EMT | 23.44 | 20.69 | 0.141 | 0.074 |
| FN1 | EMT | 14.06 | 14.32 | 0.062 | 0.136 |

Small-number caution: Cancer Cycling has only 64 cells, so one cell = 1.56 percentage points. Differences of a few percent in rare genes (SOX2, POU5F1, CDH2, ZEB1) mean one or two cells.

### 6.3 Dot plot of 12 genes × 2 populations [DONE – verified 2026-11-22]
```python
sc.pl.dotplot(cancer, var_names=present_genes, groupby="celltype_minor", standard_scale=None)
```
Dot size = % cells expressing; colour = mean expression (unscaled, 0–~0.65). Saved file was opened and checked: it shows both populations and all 12 genes, and matches the table in §6.2.

### 6.4 Violin plots: CD44, PROM1, VIM, CDH1 [DONE – verified 2026-11-22]
```python
plt.close("all")
sc.pl.violin(cancer, keys=["CD44", "PROM1", "VIM", "CDH1"], groupby="celltype_minor", rotation=90, show=False)
plt.savefig(".../results/figures/CID3921_Cancer_Epithelial_Violin_CD44_PROM1_VIM_CDH1.png", dpi=300, bbox_inches="tight")
plt.close("all")
```
- Values are log-normalised `X`; dots = individual cells; Scanpy defaults otherwise.
- **Why these genes:** CD44 and PROM1 are the most-detected stemness candidates; CDH1 and VIM are the epithelial and mesenchymal markers with the clearest between-population difference.
- **Observations (visual, not tested):**
  - CD44, PROM1 and CDH1 each show a mass at 0 (dropout or non-expression) plus a second mode around 0.5–0.8 among expressing cells.
  - The height of that expressing-cell mode looks similar in both populations. This suggests the dot-plot differences come mainly from detection frequency rather than expression level in positive cells.
  - VIM is zero in most cells of both populations. A small number of Her2 SC cells reach values up to ~3.
- **Cautions:**
  - Violin widths are scaled to the same maximum by default, so width does not reflect cell number.
  - The Cancer Cycling outline (n = 64) is less reliable.

### 6.5 Mean expression among positive cells only [DONE – verified 2026-11-22]
```python
expr = cancer_expression[present_genes]
pos = expr.where(expr > 0)            # zeros -> NaN, skipped by mean()
pos["celltype_minor"] = cancer_expression["celltype_minor"].values
pos_means = pos.groupby("celltype_minor")[present_genes].mean()
pos_n = pos.groupby("celltype_minor")[present_genes].count()
```
**Checks:**
- Positive-only mean × detection fraction = overall mean (`np.allclose` → `True`).
- Positive-cell counts ÷ group size × 100 = detection % (`True`).
- The second check was done on the saved CSVs.

| Gene | Mean in positive cells: Cycling (n pos / 64) | Mean in positive cells: Her2 SC (n pos / 377) |
|---|---:|---:|
| SOX2 | 0.307 (1) | 0.257 (3) |
| POU5F1 | 0.210 (1) | 0.218 (8) |
| NANOG | – (0) | – (0) |
| ALDH1A1 | 0.554 (4) | 0.572 (33) |
| PROM1 | 0.671 (38) | 0.651 (251) |
| CD44 | 0.750 (39) | 0.769 (261) |
| VIM | 0.489 (10) | 0.626 (107) |
| CDH1 | 0.830 (42) | 0.818 (300) |
| CDH2 | 0.254 (1) | 0.333 (8) |
| ZEB1 | – (0) | 0.358 (2) |
| SNAI1 | 0.603 (15) | 0.358 (78) |
| FN1 | 0.440 (9) | 0.946 (54) |

**Interpretation (descriptive, one patient, no statistical test):**
- **CD44, PROM1, CDH1:** positive-cell means are almost identical between the populations (differences ≤ 0.02). This confirms the §6.4 impression: their between-population differences come from detection frequency, not expression level in expressing cells.
- **FN1:** detection is similar (~14% in both), but positive Her2 SC cells average about twice the level (0.946 vs 0.440). The higher overall FN1 mean in Her2 SC comes from expression *level*. However, it rests on only 9 positive Cycling cells.
- **SNAI1:** positive Cycling cells have a higher level (0.603 vs 0.358; 15 positive cells). This explains most of SNAI1's higher overall mean in Cycling.
- **VIM:** a modest level difference (0.489 vs 0.626), plus the detection difference. Only 10 positive Cycling cells.
- **Genes with ≤ 8 positive cells** (SOX2, POU5F1, CDH2, ZEB1): positive-only means are not interpretable.

### 6.6 Second patient: CID4066 (HER2+) [DONE – verified 2026-11-22]
**Processing:** identical to CID3921.
- Built from `data/raw/CID4066/`. Gene order is identical to CID3921 (`True`).
- 13 MT genes. Recalculated `total_counts` equal `nCount_RNA` (`True`).
- QC (≥200 genes, <20% mito): 5,309 of 5,309 cells passed; none removed.
- Raw object saved **before** normalisation: `data/processed/CID4066.h5ad` (100.5 MB; `X.max()` = 5456).
- Normalised once (10,000 counts per cell + log1p; raw counts in `layers["counts"]`), `X.max()` = 8.509, saved as `data/processed/CID4066_normalised.h5ad` (165.1 MB).
- Subset `cancer4066`: Cancer Her2 SC (393) + Cancer Cycling (112) = 505 cells. LumA/LumB/Basal SC (16 cells) excluded.
- `celltype_minor` was converted with `.astype(str)` before grouping. Saving converts it to a categorical type, which would otherwise create empty rows for the excluded labels.

**Tables:**
- `results/tables/CID4066_gene_detection_percentages.csv`
- `results/tables/CID4066_mean_gene_expression.csv`
- `results/tables/CID4066_mean_expression_positive_cells.csv`
- `results/tables/CID4066_n_positive_cells.csv`

Checks on the saved CSVs: positive-mean × detection = mean (`True`); n positive ÷ group size = detection (`True`).

**Detection % — CID3921 vs CID4066**
| Gene | 3921 Cycling | 4066 Cycling | 3921 Her2 SC | 4066 Her2 SC |
|---|---:|---:|---:|---:|
| SOX2 | 1.6 | 3.6 | 0.8 | 3.6 |
| POU5F1 | 1.6 | 0.0 | 2.1 | 0.0 |
| NANOG | 0.0 | 0.9 | 0.0 | 0.3 |
| ALDH1A1 | 6.3 | 13.4 | 8.8 | 15.5 |
| PROM1 | 59.4 | 46.4 | 66.6 | 44.3 |
| CD44 | 60.9 | 36.6 | 69.2 | 33.8 |
| VIM | 15.6 | 13.4 | 28.4 | 12.7 |
| CDH1 | 65.6 | 60.7 | 79.6 | 64.1 |
| CDH2 | 1.6 | 0.9 | 2.1 | 1.0 |
| ZEB1 | 0.0 | 0.0 | 0.5 | 1.5 |
| SNAI1 | 23.4 | 12.5 | 20.7 | 11.7 |
| FN1 | 14.1 | 8.0 | 14.3 | 6.6 |

**Figure:** `results/figures/CID4066_Cancer_Epithelial_Stemness_EMT_DotPlot.png` [DONE – verified 2026-11-22]
- Same settings as CID3921, but grouped by a string copy of the labels (`cancer4066.obs["population"]`) so unused categories cannot appear.
- The image shows 2 rows and 12 genes, consistent with the tables.
- Its colour scale runs to ~0.5 (CID3921: ~0.65), and its size legend is 10–70% (CID3921: 20–80%). Per-patient dot plots are therefore **not** directly comparable by colour or size. A combined plot on a shared scale is planned.

**Mean in positive cells — CID4066** (n positive: Cycling / Her2 SC)
| Gene | Cycling | Her2 SC |
|---|---:|---:|
| ALDH1A1 | 0.496 (15) | 0.797 (61) |
| PROM1 | 0.508 (52) | 0.710 (174) |
| CD44 | 0.608 (41) | 0.624 (133) |
| VIM | 1.278 (15) | 1.112 (50) |
| CDH1 | 0.742 (68) | 0.792 (252) |
| SNAI1 | 0.326 (14) | 0.484 (46) |
| FN1 | 0.681 (9) | 0.828 (26) |
| SOX2 | 0.669 (4) | 0.553 (14) |

Genes with ≤ 6 positive cells per group (NANOG, CDH2, ZEB1) are not interpretable; POU5F1 had 0 positive cells.

**Sequencing depth** (median per cell; from `metadata.csv`)
| Patient | Cycling: counts / genes | Her2 SC: counts / genes |
|---|---|---|
| CID3921 | 25,800 / 4,898 | 28,827 / 5,052 |
| CID4066 | 20,737 / 4,198 | 16,687 / 3,690 |
| CID45171 | 10,658 / 2,726 | 13,723 / 3,022 |

**Observations (descriptive, two patients, not tested):**
- **Replicated:** the same overall pattern appears in both patients. CDH1, PROM1 and CD44 are the most widely detected candidates. SOX2, POU5F1, NANOG, CDH2 and ZEB1 are rare or absent.
- **Not replicated:** in CID3921, Her2 SC had higher detection than Cycling for CD44, PROM1 and VIM. In CID4066 those differences vanish or slightly reverse (CD44 33.8 vs 36.6; PROM1 44.3 vs 46.4; VIM 12.7 vs 13.4).
- **Likely technical contributor:** detection % depends on sequencing depth.
  - In CID3921, Her2 SC cells have slightly more genes per cell than Cycling cells.
  - In CID4066 it is the reverse (3,690 vs 4,198).
  - The direction of the detection differences follows the direction of the depth differences. This suggests within-patient detection differences may partly reflect depth, not biology.
  - The same applies between patients: CID4066 has lower depth and generally lower detection of CD44 and PROM1.
- **Implication:**
  - Compare *patterns* across patients, not raw detection percentages.
  - Treat small within-patient detection differences as unreliable.
  - Mean normalised expression adjusts for library size but is still affected by dropout.

### 6.7 Third patient: CID45171 (HER2+) [DONE – verified 2026-11-22]
**Processing:** identical to CID4066.
- 13 MT genes. `total_counts` equal `nCount_RNA` (`True`).
- 2,447 of 2,447 cells passed QC.
- Raw saved first: `data/processed/CID45171.h5ad` (50.6 MB; counts max 8601).
- Normalised once (`X.max()` = 8.516): `data/processed/CID45171_normalised.h5ad` (82.3 MB).
- Subset `cancer45171`: Cancer Her2 SC (567) + Cancer Cycling (236) = 803 cells. LumB/Basal SC (10 cells) excluded.
- Gene/barcode order check: the result was not pasted into the conversation.

**Tables:**
- `results/tables/CID45171_gene_detection_percentages.csv`
- `results/tables/CID45171_mean_gene_expression.csv`
- `results/tables/CID45171_mean_expression_positive_cells.csv`
- `results/tables/CID45171_n_positive_cells.csv`

Checks on the saved CSVs: positive-mean × detection = mean (`True`); n positive ÷ group size = detection (`True`).

**Figure:** `results/figures/CID45171_Cancer_Epithelial_Stemness_EMT_DotPlot.png` [DONE – verified]. 2 rows, 12 genes; colour scale to ~0.58; size legend 10–60%.

### 6.8 Three-patient comparison [DONE – descriptive, 2026-11-22]
**Detection % (Cycling / Her2 SC)**
| Gene | CID3921 | CID4066 | CID45171 |
|---|---|---|---|
| SOX2 | 1.6 / 0.8 | 3.6 / 3.6 | 1.3 / 1.2 |
| POU5F1 | 1.6 / 2.1 | 0.0 / 0.0 | 0.4 / 0.0 |
| NANOG | 0.0 / 0.0 | 0.9 / 0.3 | 0.0 / 0.0 |
| ALDH1A1 | 6.2 / 8.8 | 13.4 / 15.5 | 9.3 / 4.2 |
| PROM1 | 59.4 / 66.6 | 46.4 / 44.3 | **8.1 / 6.9** |
| CD44 | 60.9 / 69.2 | 36.6 / 33.8 | 58.9 / 54.7 |
| VIM | 15.6 / 28.4 | 13.4 / 12.7 | 11.9 / 13.2 |
| CDH1 | 65.6 / 79.6 | 60.7 / 64.1 | 43.6 / 52.9 |
| CDH2 | 1.6 / 2.1 | 0.9 / 1.0 | 0.4 / 0.7 |
| ZEB1 | 0.0 / 0.5 | 0.0 / 1.5 | 0.8 / 0.5 |
| SNAI1 | 23.4 / 20.7 | 12.5 / 11.7 | 6.8 / 6.3 |
| FN1 | 14.1 / 14.3 | 8.0 / 6.6 | 0.8 / 1.2 |

**Her2 SC minus Cycling, detection percentage points** (genes with meaningful detection)
| Gene | CID3921 | CID4066 | CID45171 | Same sign in all 3? |
|---|---:|---:|---:|---|
| CDH1 | +14.0 | +3.4 | +9.3 | Yes (Her2 SC higher) |
| CD44 | +8.3 | −2.8 | −4.2 | No |
| PROM1 | +7.2 | −2.2 | −1.2 | No |
| VIM | +12.8 | −0.7 | +1.4 | No (only CID3921 is large) |
| ALDH1A1 | +2.5 | +2.1 | −5.1 | No |
| SNAI1 | −2.7 | −0.8 | −0.4 | Yes, but all differences are small |

**Findings (descriptive; n = 3 patients; no formal test):**
1. **Consistent across patients:**
   - CD44 and CDH1 are widely detected in malignant epithelial cells of all three HER2+ tumours.
   - The pluripotency factors SOX2, POU5F1 and NANOG are rare or undetected in all three, as are CDH2 and ZEB1.
2. **PROM1 varies strongly between patients:** 44–67% in CID3921 and CID4066, but only 7–8% in CID45171. Sequencing depth alone is an unlikely explanation, because CD44 is still detected in ~55–59% of CID45171 cells at the same depth. This suggests real inter-patient heterogeneity in PROM1 (to be stated cautiously).
3. **The Her2 SC vs Cycling differences seen in CID3921 for CD44, PROM1 and VIM do not replicate.**
4. **CDH1 is the only well-detected gene with higher detection in Her2 SC than Cycling in all three patients.**
   - This holds even in CID4066, where Cycling cells were sequenced *more* deeply, so depth does not explain it there.
   - Plausible biology: cycling cells might express epithelial adhesion genes slightly less. This is a hypothesis only.
   - Statistical caution: with 3 patients, the chance of the same direction in all 3 by luck is 2 × (1/2)³ = 1/4. This is **not** statistical evidence.
5. **Detection of SNAI1 and FN1 falls with sequencing depth across patients** (CID3921 > CID4066 > CID45171), consistent with depth-dependent dropout.

### 6.9 Combined three-patient figure and table [DONE – verified 2026-11-22]
```python
c1 = cancer[:, present_genes].copy(); c2 = cancer4066[:, present_genes].copy(); c3 = cancer45171[:, present_genes].copy()
combined = sc.concat([c1, c2, c3])                      # 1,749 cells x 12 genes
combined.obs["group"] = combined.obs["orig.ident"].astype(str) + " | " + combined.obs["celltype_minor"].astype(str)
sc.pl.dotplot(combined, var_names={"Stemness": stemness_genes, "EMT": emt_genes}, groupby="group", standard_scale=None, show=False)
```
Each patient was normalised identically (10,000 counts per cell + log1p), so the values share one scale.

**Figure:** `results/figures/AllPatients_HER2_Cancer_Epithelial_Stemness_EMT_DotPlot.png`
- 6 rows (3 patients × 2 populations), 12 genes with Stemness/EMT brackets.
- One shared colour scale (0 to ~0.65) and one shared size legend (20–80%).
- Checked visually.

**Table:** `results/tables/AllPatients_gene_expression_summary.csv`
- Long format, 72 rows × 8 columns: `patient, population, gene, n_cells, n_positive, pct_detected, mean_expression, mean_in_positive_cells`.
- `mean_in_positive_cells` = `mean_expression` / (`pct_detected` / 100); NaN where no cells are positive.
- Every value was cross-checked against the 12 per-patient CSVs: **all match**.

**What the shared-scale figure shows:**
- CDH1 and CD44 are prominent in every row.
- PROM1 almost disappears in CID45171.
- Pluripotency factors are near-absent throughout.
- CID45171 CD44 dots are darker than the other patients' despite similar or lower detection, because positive-cell CD44 levels are higher (~0.95–1.0 vs ~0.6–0.77).
- The Her2 SC vs Cycling rows within each patient look broadly similar, apart from CDH1 and the CID3921-specific differences.

### 6.10 Reproducible pipeline: `scripts/run_analysis.py` [DONE – verified in a test environment 2026-11-22]
A single script reproduces every result above from the raw GEO files.
- It screens all 5 HER2+ patients.
- For each of the 3 analysed patients, one at a time, it runs: load → QC (with a QC violin plot) → save raw `.h5ad` → normalise once → save normalised `.h5ad` → malignant subset → 4 tables → dot plot. It also draws the CID3921 violin plot.
- It then makes the combined dot plot and the combined tables.
- Paths are derived from the script's own location (`Path(__file__).resolve().parents[1]`).
- Built-in safety checks (`assert`):
  - barcode order matches the metadata;
  - recalculated totals equal `nCount_RNA`;
  - mean = positive mean × detection fraction.
- **Verification:** run end-to-end on a copy of the raw data in a separate test environment (Scanpy 1.13.0, AnnData 0.13.4, Python 3.13). All 12 per-patient CSVs and all 72 rows of the combined table **match the interactively produced tables exactly**. Peak memory was ~0.94 GB.
- **[DONE 2026-11-22]** Run end-to-end on the Mac in `brast-ml` (Python 3.12, Scanpy 1.12.5, AnnData 0.13.4, NumPy 2.5.5, pandas 3.1.1, SciPy 1.18.1, matplotlib 3.11.2).
  - It completed with `Done.`, and the per-step outputs matched the interactive results (cell counts, QC, max values, subset sizes).
  - The combined summary and depth tables written on the Mac are identical to those from the test environment.
  - The only warning was a SciPy `spmatrix` DeprecationWarning, which is harmless.

New outputs from the script:
- `results/tables/HER2_patient_screening.csv`: the 5-patient screening table (§3).
- `results/tables/AllPatients_sequencing_depth.csv`: median counts and genes per patient × population (§6.6).
- `results/figures/CID3921_QC_violin.png`, `CID4066_QC_violin.png`, `CID45171_QC_violin.png`: n genes, total counts and % MT per patient. These replace the earlier manually saved QC plots.
- `AllPatients_gene_expression_summary.csv` gains a `gene_set` column.

Not reproduced by the script:
- The preliminary all-cell CID3921 analysis (§6.1), which was superseded by the malignant-cell analysis.
- `data/processed/CID3921_cancer_epithelial.h5ad` (an intermediate file that is not needed).

Also created:
- `README.md`: project description, findings, limitations, reproduction steps.
- `requirements.txt`: all six packages pinned to the versions above.

**Correction to §6.6 / §6.8:** the "detection differences follow sequencing depth" observation holds for CID3921 and CID4066. It does **not** clearly hold for CID45171, where Her2 SC cells are deeper but CD44 and PROM1 detection is slightly lower. The README states it for the two patients only.

---

## 7. Interpretation (exploratory, single patient)

**Stemness-associated genes**
- CD44 and PROM1 are the most widely detected stemness candidates in both populations, slightly higher in Cancer Her2 SC (CD44 69.2% vs 60.9%; PROM1 66.6% vs 59.4%).
- ALDH1A1 is detected in a minority of cells (6–9%). SOX2, POU5F1 and NANOG are essentially undetected. This is consistent with the general observation that these pluripotency factors are rarely detected in adult tumour scRNA-seq, although dropout prevents firm conclusions.
- CD44 and PROM1 are therefore the most informative candidates for further exploration in this sample. Their expression does not establish functional stemness.

**EMT-associated genes**
- CDH1 (epithelial) is the strongest signal in both populations, as expected for cells annotated as epithelial.
- VIM is detected in more Her2 SC cells (28.4%) than Cycling cells (15.6%).
- SNAI1 is slightly higher in Cancer Cycling (detection 23.4% vs 20.7%; mean 0.141 vs 0.074).
- FN1 detection is similar (~14%), but mean expression is higher in Her2 SC.
- CDH2 and ZEB1 are barely detected.
- CDH1 and VIM are detected in the same populations, which illustrates why EMT should not be treated as binary. Whether individual cells co-express both genes has **not** been tested yet.

**Overall:** differences between the two populations are modest. With one patient and no statistical testing, these are descriptive observations only.

---

## 8. Figures and tables

### Figures
| File | Status |
|---|---|
| `results/figures/CID3921_Cancer_Epithelial_Stemness_EMT_DotPlot.png` | [DONE – verified] 2026-11-22 |
| `results/figures/CID3921_QC_violin.png` | [DONE – verified, file exists] |
| `results/figures/CID4066_Cancer_Epithelial_Stemness_EMT_DotPlot.png` | [DONE – verified] 2026-11-22 |
| `results/figures/CID45171_Cancer_Epithelial_Stemness_EMT_DotPlot.png` | [DONE – verified] 2026-11-22 |
| `results/figures/AllPatients_HER2_Cancer_Epithelial_Stemness_EMT_DotPlot.png` | [DONE – verified] 2026-11-22 — **main figure** |
| `results/figures/CID3921_Cancer_Epithelial_Violin_CD44_PROM1_VIM_CDH1.png` | [DONE – verified] 2026-11-22 (8181×1472 px) |
| `CID3921_QC_violin.png` (project root) | [TO CHECK] Different file size from the copy in `results/figures/`, so it is a different version. Decide which to keep. |
| `CID3921_QC_counts_vs_genes.png` (project root) | [TO CHECK] Move to `results/figures/` |
| `CID3921_QC_counts_vs_mito.png` (project root) | [TO CHECK] Move to `results/figures/` |

### Tables
| File | Contents | Status |
|---|---|---|
| `results/tables/CID3921_gene_detection_percentages.csv` | Detection % by `celltype_minor`, 12 genes | [DONE – verified, values match §6.2] |
| `results/tables/CID3921_mean_gene_expression.csv` | Mean log-norm expression by `celltype_minor`, 12 genes | [DONE – verified, values match §6.2] |
| `results/tables/CID3921_mean_expression_positive_cells.csv` | Mean among expressing cells only | [DONE – verified] 2026-11-22 |
| `results/tables/CID3921_n_positive_cells.csv` | Number of expressing cells per gene/population | [DONE – verified] 2026-11-22 |
| `results/tables/CID4066_*.csv`, `results/tables/CID45171_*.csv` | Same four tables per patient | [DONE – verified] 2026-11-22 |
| `results/tables/AllPatients_gene_expression_summary.csv` | All patients, long format (72 rows) | [DONE – verified against per-patient tables] 2026-11-22 |

---

## 9. Open issues

1. ~~Raw file overwritten~~ **[RESOLVED 2026-11-22]**, see §4.4.
2. **[TO CHECK]** Contents of `data/processed/CID3921_cancer_epithelial.h5ad`: confirm 441 cells and that `X` is log-normalised with a `counts` layer.
3. **[TO CHECK]** Tidy up the QC figures in the project root (§8).
4. `README.md` currently contains only the title `# BRaST-ML`.
---

## 10. Decision record

| Date | Decision | Reason |
|---|---|---|
| 2026-11-21/22 | Use CID3921 as the main sample | HER2+, 441 annotated cancer epithelial cells |
| 2026-11-21/22 | Exclude CID3946 | No annotated epithelial cells |
| 2026-11-22 | QC: ≥200 genes, <20% mito; no cells removed | Standard exploratory thresholds; all cells passed |
| 2026-11-22 | Normalise to 10,000 counts per cell + log1p; keep raw counts in `layers["counts"]` | Standard Scanpy workflow; keeps raw counts for later methods |
| 2026-11-22 | Use the author-provided `celltype_major`/`celltype_minor` annotations rather than re-clustering | Keeps the project focused and comparable to the original study |
| 2026-11-22 | Report detection % and mean expression (including zeros) | Two complementary, easy-to-explain summaries |
| 2026-11-22 | No ML, trajectory inference or CytoTRACE2 at this stage | Keep the analysis defensible and explainable |
| 2026-11-22 | Treat detection-% differences as depth-sensitive; record median counts/genes per group | The CID3921 Her2 SC > Cycling differences reversed in CID4066 in line with the depth difference |
| 2026-11-22 | Report positive-cell means together with n positive cells | Means from very few cells are unreliable and must be visible as such |
| 2026-11-22 | Rebuild raw `CID3921.h5ad` from `layers["counts"]` rather than re-reading the `.mtx` | Counts were verified identical to the original (row sums equal `nCount_RNA`); faster and less memory |

---

## 11. Next steps [PLANNED]
1. ~~Reproducible script~~ **[DONE, §6.10]**: verified on the Mac; versions pinned.
2. Move the three old QC plots from the project root into an untracked folder.
3. Review the README wording.
4. Publish to GitHub (no data files): add the remote, push, and add a short repository description.
5. Optional future work (listed in the README): multi-gene signature scores, depth-aware or pseudobulk comparisons, other subtypes, spatial data.
