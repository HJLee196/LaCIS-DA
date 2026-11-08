# ReadMe_Section 52 - Real Data Analysis Code - Write All

- This document describes how to use the following R scripts:
  - `99_WriteAll_4_Real 52_realData.R`
  - `99_WriteAll_6_Real 52_realData_without.R`
  - `99_WriteAll_48_Real 52_realData.R`

- These scripts generate real-data analysis scripts with different analysis settings used in the paper.
  - The configurable analysis settings are as follows:
    - `testUse`
    - `llam`
    - `nSim`
    - `down`
    - `PC`
    - `cell_name`
    - `file_path`
    - `max_iter_imp`
    - `part1`
    - `part2`
  - For detailed descriptions of these settings, please refer to `ReadMe_Section52_RealDataAnalysisCode.md`.
  - The applicable settings depend on the script. In particular, `part1` and `part2` apply only to scripts whose filenames contain `AllPart12`.
  - `99_WriteAll_48_Real 52_realData.R` does not vary `testUse` across generated scripts because `testUse` is fixed to `"BH"` in `48_mother_code_BH_realdata.R`. The BH script performs analyses using `"LR"`, `"MAST"`, and `"wilcox_limma"` within a single script. The BH procedure does not use `"LCD"`.

- The settings `nSim`, `down`, `file_path`, and `max_iter_imp` can be modified by directly assigning values to the corresponding variables.
  - For example, setting `nSim = 1` generates real-data analysis scripts with `nSim = 1`, as used in the paper.
- The following settings can be specified via corresponding `_set` variables:
  - `testUse` -> `testUse_set`
  - `llam` -> `llam_set`
  - `PC` -> `PC_set`
  - `cell_name` -> `cell_name_set`
- Each element of `PC_set` is used as the value of `PC` for the corresponding element of `cell_name_set`, based on their indices.
- For non-BH scripts, these variables accept vectors, and real-data analysis scripts are generated for **all combinations** of `testUse_set`, `llam_set`, and the paired values of `cell_name_set` and `PC_set`.
- The vectors `cell_name_set` and `PC_set` must have the same length.
- For example:
  ```r
  testUse_set = c("LCD", "LR")
  llam_set = c("0.1")
  cell_name_set = c("astrocyte", "oligo")
  PC_set = c("43", "64")
  ```
  This setting generates four real-data analysis scripts corresponding to:
  - (LCD, 0.1, astrocyte_43)
  - (LCD, 0.1, oligo_64)
  - (LR, 0.1, astrocyte_43)
  - (LR, 0.1, oligo_64)

- The current working directory must be set to the directory containing the `.R` and `.sh` files.
  - The `.R` files are real-data analysis code files, and the `.sh` files are the corresponding SLURM batch scripts.
  - Each real-data analysis `.R` file must have a corresponding `.sh` file with the same base file name.
    - For example: `43_44_mother_code_multiLR_realdata-AllPart12.R` and `43_44_mother_code_multiLR_realdata-AllPart12.sh`
  - Use `setwd()` to set the current working directory.
  - Use `getwd()` to verify the current working directory.

- The user must define `mother_file_name` as the name of the real-data analysis script to be used as a template.
  - For example:
    ```r
    mother_file_name = "43_44_mother_code_multiLR_realdata-AllPart12.R"
    ```

- After specifying all applicable analysis settings and `mother_file_name`, run the code.
  - The script generates copies of real-data analysis code files with different analysis settings, along with the corresponding SLURM batch scripts.
  - The generated files are saved in a folder named based on `mother_file_name`.
  - The file names include information about `cell_name`, `down`, `testUse`, `llam`, and `PC`.
    - For example, when `cell_name = "astrocyte"`, `down = 1000`, `testUse = "MAST"`, `llam = "0.1"`, and `PC = 43`, the generated files are:
      - `43_44_child_code_multiLR_realdata-AllPart12_cellastrocyte_down1000_testUseMAST_10llam1_PC43.R`
      - `43_44_child_code_multiLR_realdata-AllPart12_cellastrocyte_down1000_testUseMAST_10llam1_PC43.sh`

- The file `call_all.sh` can be used to submit all generated real-data analysis jobs to the cluster.