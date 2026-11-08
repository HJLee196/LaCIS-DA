# ReadMe_Section 54 - Write All

- This document describes how to use the following R scripts:
  - `99_WriteAll_3_Simulation 54_confounding.R`
  - `99_WriteAll_8_Simulation 54_confounding_debiased.R`

- These scripts generate simulation scripts with different simulation settings used in the paper.
  - The configurable simulation settings are as follows:
    - `testUse`
    - `llam`
    - `sign_strength`
    - `batch_option`
    - `CDR_impute_option`
    - `CDR_p_compute_option`
    - `nSim`
    - `down`
    - `PC`
    - `seed_num`
    - `file_path`
    - `max_iter_imp`
  - For additional descriptions of these settings, please refer to `ReadMe_Section54.md`.

- The settings `nSim`, `down`, `PC`, `seed_num`, `file_path`, and `max_iter_imp` can be modified by directly assigning values to the corresponding variables.
  - For example, setting `nSim = 20` generates simulation scripts with `nSim = 20`.
  - To reproduce the results in Section 5.4, use `nSim = 100` for decomp knockoffs and `nSim = 20` for the BH procedure.
  - `file_path` specifies the full path to the `.rds` dataset, including the file name. This path is inserted into the `readRDS()` call in each generated script.
- The following settings can be specified via corresponding `_set` variables:
  - `testUse` -> `testUse_set`
  - `llam` -> `llam_set`
  - `sign_strength` -> `sign_strength_set`
  - `batch_option` -> `batch_option_set`
  - `CDR_impute_option` -> `CDR_impute_option_set`
  - `CDR_p_compute_option` -> `CDR_p_compute_option_set`
- These variables accept vectors, and simulation scripts are generated for all combinations of `testUse_set`, `llam_set`, `sign_strength_set`, `batch_option_set`, `CDR_impute_option_set`, and `CDR_p_compute_option_set`, subject to the restrictions described below.
- For example:
  ```r
  testUse_set = c("LCD", "LR")
  llam_set = c("0.1", "0.2")
  sign_strength_set = "2"
  batch_option_set = c("with", "without")
  CDR_impute_option_set = c("with")
  CDR_p_compute_option_set = c("with")
  ```
  This setting generates eight simulation scripts corresponding to:
  - (LCD, 0.1, 2, with, with, with)
  - (LCD, 0.1, 2, without, with, with)
  - (LCD, 0.2, 2, with, with, with)
  - (LCD, 0.2, 2, without, with, with)
  - (LR, 0.1, 2, with, with, with)
  - (LR, 0.1, 2, without, with, with)
  - (LR, 0.2, 2, with, with, with)
  - (LR, 0.2, 2, without, with, with)

- The current working directory must be set to the directory containing the `.R` and `.sh` files.
  - The `.R` files are simulation scripts, and the `.sh` files are the corresponding SLURM batch scripts.
  - Each simulation `.R` file must have a corresponding `.sh` file with the same base file name.
    - For example: `35_mother_code_decomp_confounding.R` and `35_mother_code_decomp_confounding.sh`
  - Use `setwd()` to set the current working directory.
  - Use `getwd()` to verify the current working directory.

- The user must define `mother_file_name` as the name of the simulation script to be used as a template.
  - For example:
    ```r
    mother_file_name = "85_mother_code_decomp_confounding_debiased.R"
    ```

- After specifying all applicable simulation settings and `mother_file_name`, run the code.
  - The script generates copies of the simulation code with different simulation settings, along with the corresponding SLURM batch scripts.
  - The generated files are saved in folders named based on `mother_file_name` and `testUse`.
    - For example, when `mother_file_name = "85_mother_code_decomp_confounding_debiased.R"` and `testUse = "LR"`, the generated files are saved in `85_mother_code_decomp_confounding_debiased_LR`.
  - The file names include information about `testUse`, `llam`, `down`, `sign_strength`, `batch_option`, `CDR_impute_option`, and `CDR_p_compute_option`.
    - For example, when `testUse = "LR"`, `llam = "0.2"`, `down = NULL`, `sign_strength = "2"`, `batch_option = "permute"`, `CDR_impute_option = "without"`, and `CDR_p_compute_option = "with"`, the generated files are:
      - `85_child_code_decomp_confounding_debiased_downNULL_testUseLR_10llam2_10sgnStrength20_batchpermute_impCDRwithout_compCDRwith.R`
      - `85_child_code_decomp_confounding_debiased_downNULL_testUseLR_10llam2_10sgnStrength20_batchpermute_impCDRwithout_compCDRwith.sh`

- The code excludes certain combinations of simulation settings that are not applicable or not considered in the paper.
  - For the BH procedure:
    - `"LCD"` is removed from `testUse_set` because it cannot be applied to the BH procedure.
    - `llam` is fixed to `0`, regardless of the values specified in `llam_set`, because the BH procedure does not involve imputation.
    - Since the BH procedure does not involve imputation, `CDR_impute_option` does not affect the results. Therefore, `CDR_impute_option_set` is restricted to `"with"` to avoid generating redundant scripts.
  - For non-BH simulations:
    - Combinations with `CDR_impute_option = "with"` and `CDR_p_compute_option = "without"` are excluded because they are not considered in the paper.

- The file `call_all.sh` can be used to submit all generated simulation jobs to the cluster.