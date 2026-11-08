# Alsema Dataset A preprocessing for LaCIS-DA

## Overview

Dataset A uses the publicly available processed count files from [GEO GSE146639](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146639), described by [Alsema et al. (2020)](https://doi.org/10.3389/fnmol.2020.00134).

**Stage 1** preprocesses the publicly available GEO count data and published donor metadata to construct a Seurat object. It does not use the original RDS or previously reconstructed expression matrices as inputs.

**Stage 2** prepares the preprocessed dataset for the existing LaCIS-DA simulation code. It reorders genes using the original ordered list of 2,000 simulation genes and publicly available gene annotations, removes cells with missing sex information, and assigns simulation cell-group identities using a separately provided annotation file derived from exact matching with the original LaCIS-DA dataset. All retained genes and existing expression values are preserved.

These scripts do not reconstruct counts from FASTQ reads or run the LaCIS-DA simulations.

## Directory structure

```text
Preprocess/
├── README.md
├── preprocess_Alsema_DatasetA.R
├── prepare_Alsema_simulation_fullgene.R
└── metadata/
    ├── simulation_genes_2000.txt
    ├── simulation_gene_mapping.csv
    └── simulation_cell_groups_exact.csv
```

- `preprocess_Alsema_DatasetA.R`: Preprocesses the public GEO data and constructs a Seurat object.
- `prepare_Alsema_simulation_fullgene.R`: Reorders genes, removes cells with missing sex information, and assigns simulation cell-group identities.
- `simulation_genes_2000.txt`: Contains the original 2,000 simulation gene symbols in their original order, one per line.
- `simulation_gene_mapping.csv`: Contains mappings between the original gene symbols and Ensembl IDs, including unresolved mappings.
- `simulation_cell_groups_exact.csv`: Contains original cell-group annotations for 12,516 cells matched by barcode, sex, and exact counts across 1,991 mapped genes.

Downloaded data and generated outputs are not included in this directory.

## Requirements

The scripts were tested with the following R environment:

| Component | Tested version | Used in |
|---|---|---|
| R | 4.5.1 (macOS, Apple silicon) | Both stages |
| Matrix | 1.7-4 | Both stages |
| SeuratObject | 5.3.0 | Both stages |
| digest | 0.6.39 | Both stages |
| readxl | 1.6.0 | Stage 1 |

These are the tested versions; other versions may also work. Stage 2 requires SeuratObject version 5 or later. The full Seurat package is not required to run the two preprocessing scripts, but is used by the downstream LaCIS-DA simulation code.

Required packages must be installed before running the scripts.

**Computational resources:** Stage 1 processes 162 SmartSeq2 count files and two 10x libraries. Sufficient memory and disk space are needed for the input data and generated outputs.

## Public data sources

Download the following files from their official sources:

1. **GEO count data:** [GSE146639_RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146639/suppl/GSE146639_RAW.tar). This archive contains processed count data, not raw FASTQ files.
2. **Donor metadata:** [Supplementary materials from Alsema et al. (2020)](https://www.frontiersin.org/journals/molecular-neuroscience/articles/10.3389/fnmol.2020.00134/full#supplementary-material). Download the original Excel workbook (`Table 1.xlsx`) containing Table S2, which provides donor information for the single-cell data.

Place both files in a separate input directory, such as `../preprocess_data/public_inputs/` when running from `Preprocess/`.

Stage 1 verifies file integrity using the following SHA-256 hashes:

| Filename | SHA-256 |
|---|---|
| `GSE146639_RAW.tar` | `b5ee658b48554fa4ec38ed0282c6a9ffe361f943804b34240369acd3e452e251` |
| `Table 1.xlsx` | `943f7682fc04cd05466b66656e841d91c0d35e6a3ae4e75240998d85b69cdee8` |

Use the original downloaded files without modifying or resaving them. If a checksum does not match, verify that the correct file version was downloaded.

The scripts use local files and do not download data automatically.

## Stage 1 — Public-source preprocessing

The preprocessing script performs the following steps:

1. **Load count data:** Read 162 SmartSeq2 count files and two 10x libraries from the GEO archive.
2. **SmartSeq2 quality control:** For each donor, retain cells with positive total counts whose log-transformed library sizes fall within the median ±3 scaled median absolute deviations (MAD, scaling constant 1.4826). Additionally, require no more than 3,000 detected genes and at most 10% mitochondrial counts, calculated using 13 mitochondrial protein-coding genes.
3. **10x quality control:** Retain cells with at least 200 detected genes.
4. **Merge count matrices:** Combine SmartSeq2 and 10x counts by Ensembl gene ID, using the union of genes and filling absent entries with zeros. Genes are ordered by Ensembl ID.
5. **Gene filtering:** Retain genes detected in at least three cells across the combined dataset.
6. **Metadata:** Construct cell metadata using GEO filenames and the published Table S2. The 73 cells from donor `2018-077`, whose sex information is unavailable in Table S2, are retained with missing sex values.
7. **Normalization:** Apply LogNormalize with a scale factor of 10,000 using all retained genes.
8. **Save results:** Construct and save a Seurat object containing the count matrix, normalized expression data, and cell metadata.

### Preprocessing choices

The preprocessing follows the published methods where the necessary details are available. Some implementation choices, including the log-transformed library sizes, MAD scaling convention, mitochondrial gene list, and application of gene filtering to the combined dataset, are specified in this reconstruction.

For 10x data, the minimum threshold of 200 detected genes was informed by comparisons with the original dataset. This procedure does not reproduce all QC steps described in the publication, including its UMI-based, mitochondrial, and multiplet filtering procedures.

### Run Stage 1

From the repository root, run:

```bash
cd Preprocess

# Create a directory for the downloaded public data.
mkdir -p ../preprocess_data/public_inputs

# Place GSE146639_RAW.tar and "Table 1.xlsx" in the directory above.

# Run preprocessing.
Rscript --vanilla preprocess_Alsema_DatasetA.R \
  ../preprocess_data/public_inputs \
  ../preprocess_data/trackA
```

The script creates a new run directory under `../preprocess_data/trackA/`. The output directory contains `TrackA_full_gene_Seurat.rds`, gene and cell metadata, QC summaries, and validation records.

A successful run is indicated by `COMPLETE` in `STATUS.txt`.

### Expected results

| Quantity | Value |
|---|---:|
| Initial genes | 45,590 |
| Retained genes | 28,689 |
| Retained cells | 16,793 |
| SmartSeq2 cells | 11,312 |
| 10x cells | 5,481 |
| Total retained counts | 197,362,249 |
| Cells with missing sex | 73 |

## Stage 2 — Simulation preparation

The original LaCIS-DA simulation code selects the first 2,000 genes from the input Seurat object. Since the GEO-derived dataset uses Ensembl IDs rather than the original gene symbols, this script reorders the genes to make the preprocessed dataset compatible with the existing simulation code.

Of the original 2,000 simulation genes, 1,991 could be uniquely mapped to Ensembl IDs using publicly available gene annotations. The remaining nine could not be mapped unambiguously.

The script prepares the simulation input as follows:

1. **Original simulation genes (1,991):** Place the successfully mapped genes first, preserving their relative order in `simulation_genes_2000.txt`.
2. **Additional genes (9):** Fill the remaining nine positions among the first 2,000 genes with the first nine unselected genes in the preprocessed dataset's original order. These genes are not intended as biological replacements for the unresolved genes.
3. **Remaining genes (26,689):** Append all other genes, preserving their original order.
4. **Cell selection:** Remove the 73 cells with missing sex information.
5. **Cell-group assignment:** Assign the original simulation identities (`healthy_Microglia`, `ad1_Microglia`, and `ad2_Microglia`) to cells with exact-matched annotations. Assign `unmatched` to all remaining cells.

The script retains all 28,689 genes and preserves the existing count and normalized expression values for the remaining 16,720 cells. No additional normalization or quality control is performed.

### Cell-group annotations

The original LaCIS-DA simulation code uses three cell-group identities for downsampling: `healthy_Microglia`, `ad1_Microglia`, and `ad2_Microglia`.

These identities are not directly available from the public GEO metadata. To preserve compatibility with the existing simulation code, cell-group annotations were transferred from the original LaCIS-DA Seurat object to the reconstructed dataset through exact cell matching.

Matching used cell barcodes, sex information, and count profiles across 1,991 verified gene mappings. Exact matches were required to have unique matching keys in both datasets, and count values were checked directly after matching.

This procedure identified 12,516 exact-matched cells. The resulting annotations are provided in `metadata/simulation_cell_groups_exact.csv`.

| Cell identity | Number of cells |
|---|---:|
| `healthy_Microglia` | 1,568 |
| `ad1_Microglia` | 3,245 |
| `ad2_Microglia` | 7,703 |
| `unmatched` | 4,204 |
| **Total** | **16,720** |

The 4,204 unmatched cells remain in the Stage 2 Seurat object but are not assigned to any of the three original simulation groups.

The annotation file was derived using the original LaCIS-DA dataset. It is supplied as a separate input so that Stage 2 can be run without requiring the original RDS file.

### Run Stage 2

After completing Stage 1, identify the generated run directory containing `TrackA_full_gene_Seurat.rds`.

From `Preprocess/`, run:

```sh
# Replace the path below with the run directory reported by Stage 1.
STAGE1_RUN="../preprocess_data/trackA/run_YYYYMMDD_HHMMSS_suffix"

Rscript --vanilla prepare_Alsema_simulation_fullgene.R \
  "$STAGE1_RUN/TrackA_full_gene_Seurat.rds" \
  metadata/simulation_genes_2000.txt \
  metadata/simulation_gene_mapping.csv \
  metadata/simulation_cell_groups_exact.csv \
  ../preprocess_data/simulation_prepared
```

The Stage 2 output directory must not already exist.

### Expected results

The main output is `TrackA_simulation_ordered_2000_Seurat.rds`, containing:

| Quantity | Value |
|---|---:|
| Retained genes | 28,689 |
| Retained cells | 16,720 |
| `healthy_Microglia` cells | 1,568 |
| `ad1_Microglia` cells | 3,245 |
| `ad2_Microglia` cells | 7,703 |
| `unmatched` cells | 4,204 |
| Genes in the initial simulation selection | 2,000 |
| Genes passing the full-cell detection diagnostic | 1,940 |

The output directory also contains a gene-ordering manifest, validation summaries, and execution information.

A successful run is indicated by `COMPLETE` in `STATUS.txt`.

### Gene mapping

The gene mapping was constructed using the gene annotations provided with the GEO 10x data and publicly available HGNC approved gene symbols.

Only unambiguous, one-to-one mappings to Ensembl IDs present in the preprocessed dataset were accepted. The resulting mappings are provided in `metadata/simulation_gene_mapping.csv`. No additional annotation downloads are required to run Stage 2.

Stage 2 verifies the following gene metadata files against fixed SHA-256 checksums:

| File | SHA-256 |
|---|---|
| `simulation_genes_2000.txt` | `6ad96328dbf3de4a01c6ef1af178579256dc14f26f710f8fbe5b39e402f9a972` |
| `simulation_gene_mapping.csv` | `d4d5984361528c0494fbaba629e80521b0bc58a6bffc99d6e2424018df625a74` |

These two files should be used without modification. The SHA-256 checksum of `simulation_cell_groups_exact.csv` is recorded in the Stage 2 output, while its cell IDs, group labels, and group counts are validated separately.

## Relationship to the original LaCIS-DA simulation

The original simulation code selects the first 2,000 genes from the input Seurat object and subsequently filters genes detected in more than 33 cells when 30 principal components are used.

Stage 2 prepares a Seurat object compatible with this gene-selection procedure and assigns the three cell-group identities expected by the existing simulation code.

With `down = 500`, the simulation code samples 500 cells from each of the three original groups. Cells labeled `unmatched` are not included in this downsampling.

With `down = NULL`, the simulation code uses all 16,720 cells, including those labeled `unmatched`.

The generated RDS can be used by updating the input file path in the existing simulation code. The initial data loading, group-based downsampling, gene selection, and RNA layer access have been tested for compatibility.

The resulting dataset is not identical to the original LaCIS-DA RDS because the count preprocessing, retained cell population, and group-annotation coverage differ.

## Reproducibility and limitations

Stage 1 uses publicly available GEO count data and published supplementary donor metadata. It does not require the original LaCIS-DA RDS.

Because not all preprocessing details from the original study can be recovered, the resulting dataset is not an exact reconstruction of the original analysis dataset. The preprocessing choices and gene-mapping procedures are described above.

Stage 2 uses the Stage 1 output, the supplied simulation gene list and mapping file, and the cell-group annotation file. The cell-group annotations were derived from the original LaCIS-DA dataset through exact count-profile matching and therefore are not independently reconstructed from GEO metadata alone.

The annotation file allows Stage 2 to run without the original RDS, but the original dataset was used to establish the annotations. Cells without exact-matched annotations are retained as `unmatched`.

The reconstructed dataset and its simulation results may differ from those obtained using the original LaCIS-DA analysis dataset.

## Citation and license

If you use these scripts or the resulting data, please cite:

- **Alsema et al. (2020):** Alsema AM, Jiang Q, Kracht L, et al. *Profiling Microglia From Alzheimer's Disease Donors and Non-demented Elderly in Acute Human Postmortem Cortical Tissue.* Frontiers in Molecular Neuroscience, 13:134. [https://doi.org/10.3389/fnmol.2020.00134](https://doi.org/10.3389/fnmol.2020.00134)
- **Data source:** [GEO GSE146639](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146639).
- **Gene annotations:** [HGNC](https://www.genenames.org/) (RRID:SCR_002827).

Please also cite the LaCIS-DA publication as indicated in the main repository.

The preprocessing scripts are covered by the license specified in the parent repository's LICENSE file. The downloaded GEO files and generated Seurat objects are not redistributed here.