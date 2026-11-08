#!/usr/bin/env Rscript

# ============================================================================
# Alsema Dataset A preprocessing for LaCIS-DA
# ============================================================================
# Dataset: GEO GSE146639; Alsema et al. (2020)
# DOI: https://doi.org/10.3389/fnmol.2020.00134
#
# This script preprocesses publicly available GEO count data and donor
# metadata to construct a Seurat object.
#
# The original RDS and previously reconstructed expression matrices are
# not used as inputs.
#
# Required input files:
#   1. GSE146639_RAW.tar
#      Public processed count data (not raw FASTQ reads).
#      https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146639/suppl/GSE146639_RAW.tar
#
#   2. Table 1.xlsx
#      Published supplementary workbook containing worksheet "Table S2".
#      https://www.frontiersin.org/journals/molecular-neuroscience/articles/10.3389/fnmol.2020.00134/full#supplementary-material
#
# Usage:
#   Rscript --vanilla preprocess_Alsema_DatasetA.R [input_dir] [output_root]
#
# Output:
#   A Seurat RDS containing RNA counts, normalized expression, and metadata,
#   along with QC summaries and preprocessing information.
#
# Required R packages: Matrix, SeuratObject, readxl, digest.

# ============================================================================
# Step 0 -- Configuration
# ============================================================================

# Default paths assume execution from the Preprocess/ directory.
# Command-line arguments can override these defaults.
input_dir <- "../preprocess_data/public_inputs"
output_root <- "../preprocess_data/trackA"

arguments <- commandArgs(trailingOnly = TRUE)

if (length(arguments) > 2L) {
  stop(
    "Usage: Rscript preprocess_Alsema_DatasetA.R [input_dir] [output_root]"
  )
}

if (length(arguments) >= 1L) {
  input_dir <- arguments[1L]
}

if (length(arguments) >= 2L) {
  output_root <- arguments[2L]
}

# Fixed preprocessing parameters used in the validated reconstruction.
# Some choices were informed by comparisons with the original dataset.
# Changing these parameters may produce different preprocessing results.
mad_multiplier <- 3
mad_constant <- 1.4826
smart_max_detected <- 3000L
smart_max_mt_percent <- 10
tenx_min_detected <- 200L
global_min_cells <- 3L
normalization_scale <- 10000

for (package in c("Matrix", "SeuratObject", "readxl", "digest")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Required package is not installed: ", package, ". Install it separately.")
  }
}
options(Seurat.object.assay.version = "v5")

# The 13 mitochondrial protein-coding genes used for QC.
# Gene IDs and symbols were obtained from public HGNC annotations
# and cross-checked with the GEO 10x gene annotations.
#
# The original publication specifies a mitochondrial QC threshold,
# but the exact gene panel used by the authors is not established.
#
# HGNC annotation source:
# https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
# SHA-256 of the annotation snapshot:
# 12363916be36f8f6f91879b0faade4eeaac973722de7136ec61329b77f2c2f91
mt13 <- data.frame(
  ensembl_id = c("ENSG00000198695", "ENSG00000198712", "ENSG00000198727",
                 "ENSG00000198763", "ENSG00000198786", "ENSG00000198804",
                 "ENSG00000198840", "ENSG00000198886", "ENSG00000198888",
                 "ENSG00000198899", "ENSG00000198938", "ENSG00000212907",
                 "ENSG00000228253"),
  symbol = c("MT-ND6", "MT-CO2", "MT-CYB", "MT-ND2", "MT-ND5", "MT-CO1",
             "MT-ND3", "MT-ND4", "MT-ND1", "MT-ATP6", "MT-CO3", "MT-ND4L", "MT-ATP8")
)

# Small helpers used below: validity checks and lossless sparse row alignment.
check <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}
check_ids <- function(ids, label) {
  check(!anyNA(ids) && all(nzchar(ids)) && !anyDuplicated(ids),
        paste(label, "must be nonmissing, nonempty and unique"))
}
check_counts <- function(m) {
  check(inherits(m, "dgCMatrix"), "Counts must be a sparse dgCMatrix")
  check(all(is.finite(m@x)) && all(m@x >= 0) && all(m@x == round(m@x)),
        "Counts must be finite, nonnegative integers")
  check_ids(rownames(m), "Gene IDs")
  check_ids(colnames(m), "Source cell IDs")
}
merge_on_ensembl <- function(matrices) {
  # Exact Ensembl IDs only: no alias mapping, symbol-based deletion or aggregation.
  genes <- sort(unique(unlist(lapply(matrices, rownames), use.names = FALSE)), method = "radix")
  check_ids(unlist(lapply(matrices, colnames), use.names = FALSE), "Merged cell IDs")
  aligned <- lapply(matrices, function(m) {
    positions <- match(rownames(m), genes)
    check(!anyNA(positions), "A source gene is missing from the union")
    entries <- Matrix::summary(m)
    Matrix::sparseMatrix(i = positions[entries$i], j = entries$j, x = entries$x,
                         dims = c(length(genes), ncol(m)),
                         dimnames = list(genes, colnames(m)))
  })
  do.call(cbind, aligned)
}
write_csv <- function(x, filename) write.csv(x, filename, row.names = FALSE, na = "NA")
sha256 <- function(filename) digest::digest(file = filename, algo = "sha256", serialize = FALSE)

# A fresh run directory prevents replacement of prior outputs. This function only
# groups one sequential execution so that failures can leave a clear status file.
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
run_dir <- tempfile(paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_"), tmpdir = output_root)
check(dir.create(run_dir), "Cannot create a fresh output directory")
writeLines("RUNNING", file.path(run_dir, "STATUS.txt"))

preprocess <- function() {
  # ===========================================================================
  # Step 1 -- Public data loading
  # ===========================================================================
  archive <- file.path(input_dir, "GSE146639_RAW.tar")
  workbook <- file.path(input_dir, "Table 1.xlsx")
  check(file.exists(archive) && file.exists(workbook), "Supply both documented public input files")
  expected_hashes <- c(
    GSE146639_RAW.tar = "b5ee658b48554fa4ec38ed0282c6a9ffe361f943804b34240369acd3e452e251",
    `Table 1.xlsx` = "943f7682fc04cd05466b66656e841d91c0d35e6a3ae4e75240998d85b69cdee8"
  )
  actual_hashes <- vapply(c(archive, workbook), sha256, character(1))
  check(identical(unname(actual_hashes), unname(expected_hashes)),
        "Public input bytes differ from the validated sources; investigate, do not replace expected hashes silently")
  write_csv(
    data.frame(
      filename = c("GSE146639_RAW.tar", "Table 1.xlsx"),
      sha256 = unname(actual_hashes)
    ),
    file.path(run_dir, "public_source_files.csv")
  )

  # Select the SmartSeq2 count files and the two 10x libraries from the
  # public GEO archive. Unrelated files in the archive are not used.
  # All 162 publicly available SmartSeq2 count files are included.
  members <- utils::untar(archive, list = TRUE)
  smart_files <- members[grepl("^GSM[0-9]+_.*_completeCounts\\.txt\\.gz$", members)]
  smart_gsm <- sub("_.*$", "", smart_files)
  ordering <- order(smart_gsm)
  smart_files <- smart_files[ordering]
  smart_gsm <- smart_gsm[ordering]
  tenx_files <- c(
    "GSM4403285_10X_Donor2019-010_barcodes.tsv.gz",
    "GSM4403285_10X_Donor2019-010_genes.tsv.gz",
    "GSM4403285_10X_Donor2019-010_matrix.mtx.gz",
    "GSM4403286_10X_Donor2018-135_barcodes.tsv.gz",
    "GSM4403286_10X_Donor2018-135_genes.tsv.gz",
    "GSM4403286_10X_Donor2018-135_matrix.mtx.gz"
  )
  selected_files <- c(smart_files, tenx_files)
  check(length(smart_files) == 162L && !anyDuplicated(smart_gsm), "Unexpected Smart file set")
  check(all(selected_files %in% members) && all(selected_files == basename(selected_files)),
        "Archive members are missing or are not flat filenames")
  extracted <- file.path(run_dir, "extracted_public_inputs")
  check(dir.create(extracted), "Cannot create extraction directory")
  utils::untar(archive, files = selected_files, exdir = extracted)
  write_csv(data.frame(filename = selected_files,
                       role = c(rep("SmartSeq2", length(smart_files)), rep("10x", 6)),
                       sha256 = vapply(file.path(extracted, selected_files), sha256, character(1))),
            file.path(run_dir, "input_manifest.csv"))

  # Smart files are tab-separated gene-by-barcode tables. Preserve original column
  # labels (check.names=FALSE), numeric values, and within-file barcode ordering.
  smart_matrices <- smart_metadata <- vector("list", length(smart_files))
  for (i in seq_along(smart_files)) {
    if (i == 1L || i %% 20L == 0L) cat("Reading Smart pool", i, "of", length(smart_files), "\n")
    d <- read.delim(gzfile(file.path(extracted, smart_files[i])), header = TRUE,
                    sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
    check(names(d)[1] == "gene" && ncol(d) >= 2L, "Invalid Smart table header")
    genes <- as.character(d[[1]])
    barcodes <- names(d)[-1]
    check_ids(genes, "Smart genes")
    check_ids(barcodes, "Smart barcodes")
    dense <- as.matrix(d[, -1, drop = FALSE])
    storage.mode(dense) <- "double"
    m <- as(Matrix::Matrix(dense, sparse = TRUE), "dgCMatrix")
    dimnames(m) <- list(genes, paste(smart_gsm[i], barcodes, sep = "|"))
    check_counts(m)
    descriptor <- sub("_completeCounts\\.txt\\.gz$", "", sub("^GSM[0-9]+_", "", smart_files[i]))
    check(grepl("_(CTR|AD)_", descriptor), "Unexpected Smart donor/condition filename")
    smart_matrices[[i]] <- m
    smart_metadata[[i]] <- data.frame(
      source_cell_id = colnames(m), technology = "SmartSeq2", gsm = smart_gsm[i],
      barcode = barcodes, source_donor = sub("_(CTR|AD)_.*$", "", descriptor),
      source_condition = sub("^.*_(CTR|AD)_.*$", "\\1", descriptor),
      sample_descriptor = descriptor, selection_rule = "FULL_QC_FIRST")
  }
  smart_counts <- merge_on_ensembl(smart_matrices)
  smart_md <- do.call(rbind, smart_metadata)
  check(identical(colnames(smart_counts), smart_md$source_cell_id), "Smart metadata ordering differs")
  smart_feature_ids <- rownames(smart_counts)
  stages <- data.frame(stage = "SmartSeq2_all", genes = nrow(smart_counts),
                       cells = ncol(smart_counts), total_counts = sum(smart_counts@x))
  rm(smart_matrices, smart_metadata, d, dense, m)

  # ===========================================================================
  # Step 2 -- SmartSeq2 QC
  # ===========================================================================
  # Calculate QC metrics using all SmartSeq2 genes before gene filtering.
  # Unlike 10x cells, SmartSeq2 cells are not subject to a minimum
  # detected-gene threshold of 200.
  total <- as.numeric(Matrix::colSums(smart_counts))
  detected <- as.integer(Matrix::colSums(smart_counts > 0))
  check(all(mt13$ensembl_id %in% rownames(smart_counts)), "Missing MT13 Smart features")
  mt_percent <- 100 * as.numeric(Matrix::colSums(smart_counts[mt13$ensembl_id, , drop = FALSE])) / total
  log_total <- rep(NA_real_, length(total))
  log_total[total > 0] <- log(total[total > 0])
  # Apply donor-specific median +/- 3 MAD filtering to log library sizes.
  # MAD is calculated using R's scaling constant of 1.4826:
  # MAD = 1.4826 * median(abs(x - median(x))).
  #
  # The log transformation and MAD scaling convention are implementation
  # choices, as the exact conventions used in the original study are unclear.
  # Cells on the threshold boundaries are retained.
  thresholds <- do.call(rbind, lapply(sort(unique(smart_md$source_donor)), function(donor) {
    v <- log_total[smart_md$source_donor == donor & !is.na(log_total)]
    check(length(v) > 2L, "Too few positive cells in donor")
    center <- median(v)
    spread <- mad(v, center = center, constant = mad_constant)
    check(is.finite(spread) && spread > 0, "Donor MAD is zero or nonfinite")
    data.frame(source_donor = donor, median_log_total = center, MAD_log_total = spread,
               lower = center - mad_multiplier * spread, upper = center + mad_multiplier * spread)
  }))
  donor_index <- match(smart_md$source_donor, thresholds$source_donor)
  pass_mad <- !is.na(log_total) & log_total >= thresholds$lower[donor_index] &
    log_total <= thresholds$upper[donor_index]
  pass_detected <- detected <= smart_max_detected
  pass_mt <- !is.na(mt_percent) & mt_percent <= smart_max_mt_percent
  keep_smart <- pass_mad & pass_detected & pass_mt
  smart_qc <- data.frame(source_cell_id = smart_md$source_cell_id, total_counts = total,
                         detected_genes = detected, pct_MT13 = mt_percent,
                         pass_MAD = pass_mad, pass_detected = pass_detected,
                         pass_MT13 = pass_mt, retained = keep_smart)
  write_csv(smart_qc, file.path(run_dir, "SmartSeq2_cell_QC.csv"))
  write_csv(thresholds, file.path(run_dir, "SmartSeq2_donor_thresholds.csv"))
  smart_counts <- smart_counts[, keep_smart, drop = FALSE]
  smart_md <- smart_md[keep_smart, , drop = FALSE]
  stages <- rbind(stages, data.frame(stage = "SmartSeq2_QC", genes = nrow(smart_counts),
                                     cells = ncol(smart_counts), total_counts = sum(smart_counts@x)))

  # ===========================================================================
  # Step 3 -- 10x QC (loading each public triplet immediately before its QC)
  # ===========================================================================
  # Retain 10x cells with at least 200 detected genes, calculated using
  # all features in each source 10x count matrix.
  #
  # This threshold was informed by comparisons with the original dataset.
  # The procedure does not implement all QC steps described in the publication,
  # including its UMI-based, mitochondrial, and multiplet filtering procedures.
  tenx_matrices <- tenx_metadata <- tenx_qc <- list()
  tenx_features <- NULL
  for (gsm in c("GSM4403285", "GSM4403286")) {
    files <- tenx_files[startsWith(tenx_files, paste0(gsm, "_"))]
    barcode_file <- files[grepl("_barcodes.tsv.gz$", files)]
    gene_file <- files[grepl("_genes.tsv.gz$", files)]
    matrix_file <- files[grepl("_matrix.mtx.gz$", files)]
    cat("Reading 10x library", gsm, "\n")
    con <- gzfile(file.path(extracted, barcode_file))
    barcodes <- readLines(con, warn = FALSE)
    close(con)
    ft <- read.delim(gzfile(file.path(extracted, gene_file)), header = FALSE,
                     quote = "", comment.char = "", stringsAsFactors = FALSE)
    check(ncol(ft) >= 2L, "Invalid 10x genes table")
    con <- gzfile(file.path(extracted, matrix_file))
    m <- as(as(Matrix::readMM(con), "generalMatrix"), "CsparseMatrix")
    try(close(con), silent = TRUE)
    check(identical(dim(m), c(nrow(ft), length(barcodes))), "10x dimensions differ from labels")
    dimnames(m) <- list(as.character(ft[[1]]), paste(gsm, barcodes, sep = "|"))
    check_counts(m)
    features <- data.frame(GEO_Ensembl = as.character(ft[[1]]), gene_symbol = ft[[2]],
                            feature_type = if (ncol(ft) > 2L) ft[[3]] else "Gene Expression")
    if (is.null(tenx_features)) {
      tenx_features <- features
    } else {
      k <- match(tenx_features$GEO_Ensembl, features$GEO_Ensembl)
      check(!anyNA(k) && nrow(tenx_features) == nrow(features), "10x feature universes differ")
      check(identical(tenx_features$gene_symbol, features$gene_symbol[k]) &&
              identical(tenx_features$feature_type, features$feature_type[k]), "10x annotations differ")
      m <- m[k, , drop = FALSE]
    }
    check(identical(as.character(tenx_features$gene_symbol[match(mt13$ensembl_id, tenx_features$GEO_Ensembl)]),
                    mt13$symbol), "MT13 differs from original 10x annotation")
    total <- as.numeric(Matrix::colSums(m))
    detected <- as.integer(Matrix::colSums(m > 0))
    keep <- detected >= tenx_min_detected
    donor <- sub("^.*_Donor([^_]+)_barcodes.*$", "\\1", barcode_file)
    check(donor != barcode_file, "10x donor is absent from filename")
    tenx_qc[[gsm]] <- data.frame(source_cell_id = colnames(m), total_counts = total,
                                 detected_genes = detected, retained = keep)
    tenx_matrices[[gsm]] <- m[, keep, drop = FALSE]
    tenx_metadata[[gsm]] <- data.frame(
      source_cell_id = colnames(m)[keep], technology = "10x", gsm = gsm, barcode = barcodes[keep],
      source_donor = donor, source_condition = NA_character_,
      sample_descriptor = paste0("10X_Donor", donor), selection_rule = "S1_MIN200")
  }
  tenx_counts <- do.call(cbind, tenx_matrices)
  tenx_md <- do.call(rbind, tenx_metadata)
  con <- gzfile(file.path(run_dir, "10x_all_barcode_QC.csv.gz"), "wt")
  write.csv(do.call(rbind, tenx_qc), con, row.names = FALSE, na = "NA")
  close(con)
  stages <- rbind(stages, data.frame(stage = "10x_S1", genes = nrow(tenx_counts),
                                     cells = ncol(tenx_counts), total_counts = sum(tenx_counts@x)))
  rm(tenx_matrices, tenx_metadata, tenx_qc, m, ft)

  # ===========================================================================
  # Step 4 -- Feature alignment and merging
  # ===========================================================================
  # Merge SmartSeq2 and 10x count matrices using exact Ensembl gene IDs.
  # Use the union of genes, adding zeros where genes are absent from a
  # source dataset. These zeros do not necessarily indicate zero expression.
  #
  # Genes are sorted by Ensembl ID, while cells retain their source order.
  # No batch correction or integration is performed.
  raw <- merge_on_ensembl(list(smart_counts, tenx_counts))
  check_counts(raw)
  md <- rbind(smart_md, tenx_md)
  k <- match(colnames(raw), md$source_cell_id)
  check(!anyNA(k) && !anyDuplicated(md$source_cell_id), "Missing or duplicate metadata")
  md <- md[k, , drop = FALSE]
  stages <- rbind(stages, data.frame(stage = "joint_raw", genes = nrow(raw),
                                     cells = ncol(raw), total_counts = sum(raw@x)))

  # ===========================================================================
  # Step 5 -- Gene filtering: GLOBAL_MIN3
  # ===========================================================================
  # Retain genes detected in at least three cells across all QC-retained
  # SmartSeq2 and 10x cells, using positive counts to define detection.
  #
  # This filter is applied to the complete merged gene set, before any
  # simulation-specific gene selection.
  gene_detection <- as.integer(Matrix::rowSums(raw > 0))
  keep_gene <- gene_detection >= global_min_cells
  gene_manifest <- data.frame(
    GEO_Ensembl = rownames(raw), n_detected_cells = gene_detection, pass_GLOBAL_MIN3 = keep_gene,
    in_SmartSeq2_feature_union = rownames(raw) %in% smart_feature_ids,
    in_10x_feature_list = rownames(raw) %in% tenx_features$GEO_Ensembl,
    GEO_10x_symbol = tenx_features$gene_symbol[match(rownames(raw), tenx_features$GEO_Ensembl)])
  counts <- raw[keep_gene, , drop = FALSE]
  stages <- rbind(stages, data.frame(stage = "GLOBAL_MIN3", genes = nrow(counts),
                                     cells = ncol(counts), total_counts = sum(counts@x)))
  rm(raw, smart_counts, tenx_counts)
  invisible(gc())

  # ===========================================================================
  # Step 6 -- Public metadata
  # ===========================================================================
  # Match donor IDs from GEO filenames to the published Table S2 metadata.
  # Sex and clinical group are obtained from Table S2, while the condition
  # recorded in GEO filenames is preserved separately.
  #
  # Donor 2018-077 is absent from Table S2. Its 73 cells are retained
  # with missing sex information.
  donor_table <- as.data.frame(readxl::read_excel(workbook, sheet = "Table S2",
                                                 range = "A2:I16", col_types = "text"))
  check(all(c("DonorID", "Gender", "Group", "Region") %in% names(donor_table)), "Table S2 headers differ")
  donor_table$DonorID <- trimws(donor_table$DonorID)
  check_ids(donor_table$DonorID, "Table S2 donor IDs")
  gender <- toupper(trimws(donor_table$Gender))
  check(all(gender %in% c("F", "M")) && all(trimws(donor_table$Region) == "LPS"), "Unexpected Table S2 values")
  k <- match(md$source_donor, donor_table$DonorID)
  # Construct unique Seurat cell IDs from the original GEO source IDs.
  # Check for duplicate IDs after replacing the separators.
  md$seurat_cell_id <- paste0("SRC-", gsub("_", "-", gsub("|", "-", md$source_cell_id, fixed = TRUE), fixed = TRUE))
  check_ids(md$seurat_cell_id, "Seurat cell IDs")
  md$sex <- ifelse(gender[k] == "F", "female", "male")
  md$sex_numeric <- ifelse(is.na(md$sex), NA_integer_, as.integer(md$sex == "female"))
  md$sex_provenance <- ifelse(is.na(k), "NOT_IN_TABLE_S2", "TABLE_S2_VERIFIED")
  md$published_group <- donor_table$Group[k]
  md$Table_S2_row <- k + 2L
  rownames(md) <- md$seurat_cell_id
  colnames(counts) <- md$seurat_cell_id
  check(identical(colnames(counts), rownames(md)), "Metadata does not align with count columns")

  # ===========================================================================
  # Step 7 -- Normalization
  # ===========================================================================
  # Apply LogNormalize with a scale factor of 10,000:
  # normalized[g,j] = log(1 + 10000 * counts[g,j] / library_size[j]).
  #
  # Library sizes are calculated using all genes retained after gene filtering,
  # before any simulation-specific gene selection.
  # The transformation preserves the sparse matrix structure.
  library_size <- as.numeric(Matrix::colSums(counts))
  check(all(is.finite(library_size) & library_size > 0), "Nonpositive retained library size")
  normalized <- counts
  normalized@x <- log1p(normalized@x * rep(normalization_scale / library_size,
                                          times = diff(normalized@p)))

  # ===========================================================================
  # Step 8 -- Seurat object construction
  # ===========================================================================
  
  # Create a Seurat object with the retained counts and cell metadata.
  # Disable additional gene and cell filtering during object construction.
  object <- SeuratObject::CreateSeuratObject(
    counts = counts,
    meta.data = md,
    project = "Alsema_GEO_TrackA",
    min.cells = 0,
    min.features = 0
  )
  
  # Store the normalized expression values in the RNA data layer.
  SeuratObject::LayerData(
    object,
    assay = "RNA",
    layer = "data"
  ) <- normalized
  
  # No variable feature selection, PCA, clustering, or simulation is performed.

  # ===========================================================================
  # Step 9 -- Validation and saving
  # ===========================================================================
  
  # Verify the final count matrix.
  check_counts(counts)
  
  # Summarize the preprocessing results.
  run_summary <- data.frame(
    initial_genes = nrow(gene_manifest),
    retained_genes = nrow(counts),
    cells = ncol(counts),
    SmartSeq2 = sum(md$technology == "SmartSeq2"),
    tenx = sum(md$technology == "10x"),
    raw_total_counts = stages$total_counts[stages$stage == "joint_raw"],
    retained_total_counts = sum(counts@x),
    female = sum(md$sex == "female", na.rm = TRUE),
    male = sum(md$sex == "male", na.rm = TRUE),
    missing_sex = sum(is.na(md$sex)),
    scale_factor = normalization_scale
  )
  
  # Validate expected results and data consistency.
  # Expected values are used only for validation, not to modify the data.
  checks <- c(
    initial_union = nrow(gene_manifest) == 45590L,
    dimensions = identical(dim(counts), c(28689L, 16793L)),
    total_counts = sum(counts@x) == 197362249,
    
    technology_counts =
      run_summary$SmartSeq2 == 11312L &&
      run_summary$tenx == 5481L,
    
    sex_counts =
      run_summary$female == 10861L &&
      run_summary$male == 5859L &&
      run_summary$missing_sex == 73L,
    
    missing_sex_identity =
      identical(is.na(md$sex), md$source_donor == "2018-077"),
    
    metadata_order =
      identical(rownames(object@meta.data), colnames(counts)),
    
    gene_manifest_order =
      identical(rownames(counts), gene_manifest$GEO_Ensembl[keep_gene]),
    
    finite_data = all(is.finite(normalized@x)),
    
    layer_IDs =
      identical(dimnames(counts), dimnames(normalized))
  )
  
  # Save the validation results and stop if any check fails.
  write_csv(
    data.frame(check = names(checks), passed = unname(checks)),
    file.path(run_dir, "basic_checks.csv")
  )
  
  write_csv(run_summary, file.path(run_dir, "summary.csv"))
  
  check(
    all(checks),
    "Preprocessing validation failed; inspect basic_checks.csv"
  )
  
  # Save the preprocessed Seurat object.
  saveRDS(
    object,
    file.path(run_dir, "TrackA_full_gene_Seurat.rds"),
    compress = "gzip"
  )
  
  # Save metadata and preprocessing summaries.
  write_csv(md, file.path(run_dir, "cell_metadata.csv"))
  write_csv(gene_manifest, file.path(run_dir, "gene_manifest.csv"))
  write_csv(stages, file.path(run_dir, "stage_summary.csv"))
  
  write_csv(
    data.frame(
      source_cell_id = md$source_cell_id,
      library_size = library_size
    ),
    file.path(run_dir, "normalization_denominators.csv")
  )
  
  # Record the R environment and completion status.
  writeLines(
    capture.output(sessionInfo()),
    file.path(run_dir, "session_info.txt")
  )
  
  writeLines("COMPLETE", file.path(run_dir, "STATUS.txt"))
  
  cat("\nCompleted Dataset A preprocessing:\n")
  print(run_summary, row.names = FALSE)
  cat("Output directory:", normalizePath(run_dir), "\n")
}

tryCatch(preprocess(), error = function(e) {
  writeLines(c("FAILED", conditionMessage(e)), file.path(run_dir, "STATUS.txt"))
  stop(conditionMessage(e), call. = FALSE)
})
