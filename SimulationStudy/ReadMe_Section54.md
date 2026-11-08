# ReadMe_Section 54

- This document describes how to use the following simulation scripts:
  - `35_mother_code_decomp_confounding.R`
  - `38_mother_code_BH_confounding.R`
  - `85_mother_code_decomp_confounding_debiased.R`

- This document also describes the following summary script:
  - `99_SummaryAll_3_Simulation 54_confounding.R`

- `85_mother_code_decomp_confounding_debiased.R` performs the same simulation as `35_mother_code_decomp_confounding.R`, but uses the de-biased imputation algorithm.

## 1. Simulation scripts

- The simulation scripts implement the following procedures:
  - `35_mother_code_decomp_confounding.R`: **Decomp knockoffs**
  - `38_mother_code_BH_confounding.R`: BH procedure without knockoffs
  - `85_mother_code_decomp_confounding_debiased.R`: **Decomp knockoffs** with de-biased imputation

- Required datasets:
  - `HuVascAD_seuratcluster3.rds`

- Simulation settings to specify:
  - Dataset path: Specify the file path directly in `HuVascAD = readRDS(file = "path/to/file.rds")`.

  - `down`: The number of cells sampled from each group (Control and AD).
    - When `down` is not `NULL`, the total number of sampled cells is `2 * down`.
    - Set to `NULL` to use all cells.

  - `gene.index`: A vector specifying the genes to be used (*p*).
    - For Section 5.4, `gene.index` is `1:3000`.

  - `sign_strength`: The signal strength used to create synthetic signals.
    - In Section 5.4, this parameter is referred to as **`sgn`**.

  - `n_signal`: The number of signals in the simulation.

  - `PC`: The number of latent factors.
    - For Section 5.4, `PC` is `13`.

  - `llam`: A scalar that controls the lambda value for sc-softImpute.
    - It can take any nonnegative value. The final lambda is computed by multiplying this value by the largest singular value of the data matrix, i.e., `llam*softImpute::lambda0(data.matrix)`.

  - `nSim`: The number of simulation replicates.
    - Set `nSim = 100` for decomp knockoffs and `nSim = 20` for the BH procedure to reproduce the results in Section 5.4.

  - `testUse`: The method used for statistical testing. One of `"LCD"`, `"LR"`, `"MAST"`, or `"wilcox_limma"`.
    - `"LCD"` cannot be used with the BH procedure.
    - `"wilcox_limma"` does not use covariate adjustment through the `latent.vars` argument in `FindMarkers()`.

  - `max_iter_imp`: The maximum number of iterations allowed for sc-softImpute.

  - The following settings control which covariates are used for imputation and the computation of test statistics. Imputation-related settings apply only to the decomp-knockoff simulations:
    - `batch_option`: One of `"with"`, `"without"`, or `"permute"`.
      - Set `batch_option = "with"` to use the batch variable in the applicable imputation and test-statistic computation steps.
      - Set `batch_option = "without"` to exclude the batch variable.
      - Set `batch_option = "permute"` to permute the batch variable.
    - `CDR_impute_option`: One of `"with"` or `"without"`.
      - Set `CDR_impute_option = "with"` to use the CDR variable in the imputation step.
      - Set `CDR_impute_option = "without"` to exclude the CDR variable in the imputation step.
    - `CDR_p_compute_option`: One of `"with"` or `"without"`.
      - Set `CDR_p_compute_option = "with"` to use the CDR variable in the computation step.
      - Set `CDR_p_compute_option = "without"` to exclude the CDR variable in the computation step.

  - `seed_num`: The random seed used for reproducibility.

- The simulation scripts follow different procedures depending on the method.

  - `35_mother_code_decomp_confounding.R` and `85_mother_code_decomp_confounding_debiased.R` (Decomp knockoffs):
    1. Imputation and knockoff construction
       - sc-softImpute is applied to the preprocessed dataset to estimate the latent factor structure and impute missing components.
       - The covariance matrix is estimated using the imputation results, and Gaussian knockoffs are constructed.
       - These steps are performed once before the simulation repetitions.

    2. Simulation and variable selection
       - In each simulation replicate, synthetic signals and binary response labels are generated.
       - Test statistics are computed according to the method specified in `testUse` and stored in `W_imp` and `W_imp_b`.
         - For `"LR"`, `"MAST"`, and `"wilcox_limma"`:
           - `W_imp` contains the differences in negative log unadjusted p-values between the original variables and their knockoffs.
           - `W_imp_b` contains the differences in negative log adjusted p-values between the original variables and their knockoffs.
         - For `"LCD"`, the statistics are constructed using differences in absolute Lasso coefficient values between the original variables and their knockoffs.
       - The knockoff procedure is applied using `W_imp` and `W_imp_b` separately to select significant variables.
       - The code prints the empirical FDP and power for each selection result.

  - `38_mother_code_BH_confounding.R` (BH procedure):
    1. Data preprocessing
       - The dataset is preprocessed, including cell downsampling, gene filtering, and preparation of covariates.
       - No imputation or knockoff construction is performed.

    2. Simulation and variable selection
       - In each simulation replicate, synthetic signals and binary response labels are generated.
       - P-values are computed using the method specified in `testUse`.
       - `W_imp` contains the unadjusted p-values, and `W_imp_b` contains the Bonferroni-adjusted p-values returned by `FindMarkers()`.
       - BH correction is applied to `W_imp`, and significant variables are selected using a target FDR of 0.1.
       - For internal checks, BH correction is also applied to `W_imp_b`, but these results are not used in the paper.
       - The code prints the empirical FDP and power for both selection results.

- For all three simulation scripts, the results of all simulation replicates are stored in `FinalResult` and saved as an `.rda` file in the current working directory (`getwd()`).

## 2. Summary script

- `99_SummaryAll_3_Simulation 54_confounding.R`
  - This script generates the results presented in Table 2 of Section 5.4 and Table S5 of the Appendix.
  - This script consists of two parts:
    - Part for Table 2 (**Decomp knockoffs**)
    - Part for Table S5 (BH procedure)

- Save the simulation results in separate folders for each simulation script and test statistic.
  - The result folders shown below are separate from the folders containing the generated `.R` and `.sh` files. Organize the resulting `.rda` files into these folders before running the summary script.
  - For example:
    - `35_mother_code_decomp_confounding_LCD_rda`
    - `85_mother_code_decomp_confounding_debiased_LCD_rda`
    - `85_mother_code_decomp_confounding_debiased_LR_rda`
  - This organization prevents results from different simulation methods or test statistics from being combined during summarization.
  - Keep results with different signal strengths or sample sizes in separate folders as well.

- The current working directory must be set to the folder containing the relevant `.rda` files.
  - For the Table 2 and Table S5 parts, the working directory must be set separately.
    - For the first part, set the working directory to the folder containing the `.rda` files generated by the Decomp knockoffs.
    - For the second part, set it to the folder containing the `.rda` files generated by the BH procedure.
  - Use `setwd()` to set the current working directory.
  - Use `getwd()` to verify the current working directory.

- The script uses `W_imp_b` (or `W_imp` for the BH procedure) from the `.rda` files to compute the FDP and power for the target FDR value specified by `target_FDR_set`.

- `target_FDR_set`: A scalar specifying the target FDR value.
  - For Tables 2 and S5, `target_FDR_set = 0.1`.

- For each method and setting, the average FDP, the standard deviation of the FDP, the average power, and the standard deviation of the power are stored in four separate `.csv` files:
  - Before running the summary script, update `final_summarize_name` according to the test statistic (`"LCD"`, `"LR"`, `"MAST"`, or `"WRT"`) being summarized.
  - For example, The Table 2 block writes four CSV files: 
    - `Section54_Table2_LR_DownNULL_Avg_FDP.csv`: Average FDP
    - `Section54_Table2_LR_DownNULL_SD_FDP.csv`: SD of FDP
    - `Section54_Table2_LR_DownNULL_Avg_Power.csv`: Average power
    - `Section54_Table2_LR_DownNULL_SD_Power.csv`: SD of power
  - The Table S5 block writes four CSV files with the prefix `Section54_TableS5_BH_DownNULL`. Each prefix is also followed by `_Avg_FDP.csv`, `_SD_FDP.csv`, `_Avg_Power.csv`, or `_SD_Power.csv`.