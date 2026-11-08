# ReadMe_Section B2 - Write All

- This document describes how to use the following WriteAll scripts:
  - `99_WriteAll_0_Simulation B2.R`
  - `99_WriteAll_1_Simulation B2_Negbinomial.R`
  - `99_WriteAll_5_Simulation B2_Debiased.R`
  - `99_WriteAll_7_Simulation B2_Observed.R`

- These scripts share the same structure and file-generation procedure, differing only in the mother scripts used and the default values of `llam_set`.
  - The instructions below apply to all four scripts, with the appropriate `mother_file_name` and simulation settings.

- These WriteAll scripts generate simulation scripts with different simulation settings used in the paper.
  - The configurable simulation settings are as follows:
    - `testUse`
    - `llam`
    - `sign_strength`
    - `nSim`
    - `down`
    - `PC`
    - `seed_num`
    - `file_path`
    - `max_iter_imp`
  - For detailed descriptions of these settings, please refer to `ReadMe_SectionB2.md`.

- The settings `nSim`, `down`, `PC`, `seed_num`, `file_path`, and `max_iter_imp` can be modified by directly assigning values to the corresponding variables.
  - For example, setting `nSim = 20` generates simulation scripts with `nSim = 20`, as used for Section B2.
  - `file_path` specifies the full path to the `.rds` dataset, including the file name. This path is inserted into the `readRDS()` call in each generated script.
  - When using a relative `file_path`, ensure that the dataset is accessible from the working directory used to execute the generated scripts.
- The following settings can be specified via corresponding `_set` variables:
  - `testUse` -> `testUse_set`
  - `llam` -> `llam_set`
  - `sign_strength` -> `sign_strength_set`
- These variables accept vectors, and simulation scripts are generated for all combinations of `testUse_set`, `llam_set`, and `sign_strength_set`, subject to the restrictions described below.
- For example, when using a non-BH and non-naive script:
  ```r
  testUse_set = c("LCD", "LR")
  llam_set = c("0.1", "0.2")
  sign_strength_set = "2"
  ```
  This setting generates four simulation scripts corresponding to:
  - (LCD, 0.1, 2)
  - (LCD, 0.2, 2)
  - (LR, 0.1, 2)
  - (LR, 0.2, 2)

- The current working directory must be set to the directory containing the mother `.R` and `.sh` files.
  - The `.R` files are simulation scripts, and the `.sh` files are the corresponding SLURM batch scripts.
  - Each mother `.R` script must have a corresponding `.sh` file with the same base file name.
    - For example: `01_mother_code_asdp.R` and `01_mother_code_asdp.sh`
  - Use `setwd()` to set the current working directory.
  - Use `getwd()` to verify the current working directory.

- The user must define `mother_file_name` as the name of the simulation script to be used as a template.
  - For example:
    ```r
    mother_file_name = "01_mother_code_asdp.R"
    ```

- After specifying all applicable simulation settings and `mother_file_name`, run the code.
  - The script generates copies of the simulation code with different simulation settings, along with the corresponding SLURM batch scripts.
  - The generated files are saved in a folder named based on `mother_file_name`.
    - For example, `mother_file_name = "01_mother_code_asdp.R"` creates the folder `01_mother_code_asdp`.
  - The file names include information about `testUse`, `llam`, and `sign_strength`.
    - For example, when `testUse = "LCD"`, `llam = "0.1"`, and `sign_strength = "2"`, the generated files are:
      - `01_child_code_asdp_testUseLCD_10llam1_10sgnStrength20.R`
      - `01_child_code_asdp_testUseLCD_10llam1_10sgnStrength20.sh`

- The code excludes certain combinations of simulation settings that are not applicable.
  - For the BH procedure in `99_WriteAll_0_Simulation B2.R`:
    - `"LCD"` is removed from `testUse_set` because it cannot be applied to the BH procedure.
    - `llam` is fixed to `0` because the BH procedure does not involve imputation.
  - For standard Gaussian knockoffs in `99_WriteAll_0_Simulation B2.R`:
    - `llam` is fixed to `0` because this method does not involve imputation.

- The file `call_all.sh` can be used to submit all generated simulation jobs to the cluster.