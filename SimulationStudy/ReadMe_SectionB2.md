# ReadMe_Section B2

- This document describes how to use the following simulation scripts:
  - `01_mother_code_asdp.R`
  - `02_mother_code_LR.R`
  - `03_04_mother_code_multiLR.R`
  - `05_mother_code_decomp.R`
  - `06_07_mother_code_multidecomp.R`
  - `08_mother_code_BH.R`
  - `09_mother_code_naive.R`

- This document also describes the following summary scripts:
  - `99_SummaryAll_0_Simulation B2.R`
  - `99_SummaryAll_0_Simulation B2_appendix.R`

- Related simulation scripts use the leading numbers `1`, `5`, and `7`, as illustrated by the following examples. These variants are not covered by the instructions below:
  - **1**: `11_mother_code_asdp_PoisNegbin.R`
    - Simulation code for Section B5, using ASDP knockoffs and synthetic data generated from a negative binomial distribution.
  - **5**: `53_54_mother_code_multiLR_debias.R`
    - Simulation code for Section E5, using de-biased imputation and multi-LR knockoffs.
  - **7**: `75_mother_code_decomp_directKnockoffs.R`
    - Simulation code for Section B3, targeting pairwise exchangeability for the observed components. For computational efficiency, missing components are conditionally sampled before generating Gaussian knockoffs, and only the knockoff values at the observed indices are retained.

## 1. Simulation scripts

- The simulation scripts implement the following procedures:
  - `01_mother_code_asdp.R`: **ASDP knockoffs**
  - `02_mother_code_LR.R`: **LR knockoffs**
  - `03_04_mother_code_multiLR.R`: **Multi-LR knockoffs** (with e-BH applied in the summary scripts)
  - `05_mother_code_decomp.R`: **Decomp knockoffs**
  - `06_07_mother_code_multidecomp.R`: **Multi-decomp knockoffs** (with e-BH applied in the summary scripts)
  - `08_mother_code_BH.R`: BH procedure without knockoffs
  - `09_mother_code_naive.R`: **Standard Gaussian knockoffs**

- Required dataset:
  - `combined_superior_parietal_lobe.rds`

- Simulation settings to specify:
  - Dataset path: Specify the file path directly in `superior_parietal_lobe = readRDS(file = 'path/to/combined_superior_parietal_lobe.rds')`.

  - `down`: The number of cells sampled from each group (`healthy_Microglia`, `ad1_Microglia`, and `ad2_Microglia`).
    - When `down` is not `NULL`, the total number of sampled cells is `3 * down`.
    - Set to `NULL` to use all cells.

  - `gene.index`: A vector specifying the genes to be used (*p*).
    - For Section B2, `gene.index = 1:2000`.

  - `sign_strength`: The signal strength used to create synthetic signals.
    - In Section B2, this setting is referred to as **`sgn`**.

  - `n_signal`: The number of signals in the simulation.

  - `PC`: The number of latent factors.
    - For Section B2, `PC = 30`.

  - `llam`: A scalar that controls the lambda value for sc-softImpute.
    - It can take any nonnegative value. The final lambda is computed by multiplying this value by the largest singular value of the data matrix, i.e., `llam*softImpute::lambda0(data.matrix)`.
    - This setting is not used for imputation in the BH or standard Gaussian knockoff scripts.

  - `m_kos`: The number of knockoff copies to generate for the multiple-knockoff procedure.
    - This setting applies only to the multi-LR and multi-decomp scripts.

  - `nSim`: The number of simulation replicates.
    - Set `nSim = 20` to reproduce the results in Section B2.

  - `testUse`: The method used for statistical testing. One of `"LCD"`, `"LR"`, `"MAST"`, or `"wilcox_limma"`.
    - `"LCD"` cannot be used with the BH procedure.
    - `"wilcox_limma"` does not use covariate adjustment through the `latent.vars` argument in `FindMarkers()`.

  - `max_iter_imp`: The maximum number of iterations allowed for sc-softImpute.
    - This setting applies to scripts that perform imputation.

  - `seed_num`: The random seed used for reproducibility.

- The simulation scripts follow different procedures depending on the method:

  - `01_mother_code_asdp.R`, `02_mother_code_LR.R`, `03_04_mother_code_multiLR.R`, `05_mother_code_decomp.R`, and `06_07_mother_code_multidecomp.R`:
    1. Imputation and knockoff construction
       - sc-softImpute is applied to the preprocessed dataset to estimate the latent factor structure and impute missing components.
       - Knockoffs are constructed using the estimated model and imputation results, according to the specified knockoff construction method.
       - These steps are performed once before the simulation repetitions.

    2. Simulation and variable selection
       - In each simulation replicate, synthetic signals and binary response labels are generated.
       - Test statistics are computed according to the method specified in `testUse` and stored in `W_imp` and `W_imp_b`.
       - For single-knockoff methods (ASDP, LR, and Decomp):
         - For `"LR"`, `"MAST"`, and `"wilcox_limma"`:
           - `W_imp` contains the differences in negative log unadjusted p-values between the original variables and their knockoffs.
           - `W_imp_b` contains the differences in negative log adjusted p-values between the original variables and their knockoffs.
         - For `"LCD"`, the statistics are constructed using differences in absolute Lasso coefficient values between the original variables and their knockoffs.
       - For multiple-knockoff methods (Multi-LR and Multi-decomp):
         - `W_imp` and `W_imp_b` are matrices containing feature importance statistics for the original variables and each knockoff copy.
         - For `"LR"`, `"MAST"`, and `"wilcox_limma"`, `W_imp` uses negative log unadjusted p-values, while `W_imp_b` uses negative log adjusted p-values.
         - For `"LCD"`, the matrices contain absolute Lasso coefficient values.
       - The knockoff procedure is applied to `W_imp` and `W_imp_b` separately, and the code prints the empirical FDP and power for the resulting selections.

  - `09_mother_code_naive.R` (Standard Gaussian knockoffs):
    1. Knockoff construction
      - Missing entries are treated as zeros, and standard second-order Gaussian knockoffs are generated across all feature coordinates without imputation. The knockoff values at unobserved indices are then set to zero.
    2. Simulation and variable selection
       - Synthetic signals and binary response labels are generated in each replicate.
       - Test statistics, knockoff-based selections, empirical FDP, and power are computed as described above for single-knockoff methods.

  - `08_mother_code_BH.R` (BH procedure):
    1. Data preprocessing
       - The dataset is preprocessed, including cell downsampling, gene filtering, and preparation of covariates.
       - No imputation or knockoff construction is performed.
    2. Simulation and variable selection
       - Synthetic signals and binary response labels are generated in each replicate.
       - P-values are computed using the method specified in `testUse`.
       - `W_imp` contains the unadjusted p-values, and `W_imp_b` contains the Bonferroni-adjusted p-values returned by `FindMarkers()`.
       - BH correction is applied to `W_imp` to select significant variables at a target FDR of 0.1.
       - For internal checks, the code also applies BH correction to `W_imp_b`; the summary scripts use `W_imp`, not `W_imp_b`, for the BH results.
       - The code prints the empirical FDP and power for both selection results.

- For all seven simulation scripts, the results of the simulation replicates are stored in `FinalResult` and saved as an `.rda` file in the current working directory (`getwd()`).

## 2. Summary scripts

- `99_SummaryAll_0_Simulation B2.R`
  - This script generates Figure S9 in Section B2 and saves `SimB2_FigureS9.jpeg` and `SimB2_FigureS9.csv`.
  - Set `target_FDR_set = c(0.1, 0.01)` to reproduce Figure S9.
  - The figure summary uses the knockoff-based results and excludes BH and standard Gaussian knockoffs.

- `99_SummaryAll_0_Simulation B2_appendix.R`
  - This script generates the results presented in Table S6 of Section B2 of the Appendix and saves `SimB2_TableS6.csv`.
  - Set `target_FDR_set = c(0.05)` to reproduce Table S6.
  - This summary includes BH and standard Gaussian knockoff results in addition to the other methods.

- Set the current working directory to the folder containing the relevant `.rda` simulation results.
  - Both summary scripts search for `.rda` files recursively, including those in subdirectories.
  - Ensure that the folder and its subdirectories contain only the simulation results intended for the corresponding figure or table.
  - Use `setwd()` to set the current working directory.
  - Use `getwd()` to verify the current working directory.

- Both summary scripts use `W_imp_b` to summarize knockoff results; the Table S6 summary uses `W_imp` for the BH procedure. The scripts compute empirical FDP and power at the target FDR levels specified by `target_FDR_set`.

- `target_FDR_set`: A numeric vector specifying the target FDR levels.
  - For Figure S9, use `target_FDR_set = c(0.1, 0.01)`.
  - For Table S6, use `target_FDR_set = c(0.05)`.