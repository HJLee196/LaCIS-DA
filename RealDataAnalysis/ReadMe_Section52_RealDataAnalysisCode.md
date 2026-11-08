# ReadMe_Section 52 - Real Data Analysis Code

- This document describes how to use the following real-data analysis scripts:
  - `43_44_mother_code_multiLR_realdata-AllPart12.R`
  - `46_47_mother_code_multidecomp_realdata-AllPart12.R`
  - `48_mother_code_BH_realdata.R`
  - `66_67_mother_code_multidecomp_realdata-AllPart12.R`

- This document also describes the following summary script:
  - `99_SummaryAll_4_Real 52_realdata.R`

## 1. Real-data analysis scripts

- `43_44_mother_code_multiLR_realdata-AllPart12.R`
  - Code for reproducing the multi-LR/e-LR results in the real-data analysis in Section 5.2.
- `46_47_mother_code_multidecomp_realdata-AllPart12.R`
  - Code for reproducing the multi-decomp/e-decomp results in the real-data analysis in Section 5.2.
- `48_mother_code_BH_realdata.R`
  - Code for reproducing the BH results in the real-data analysis in Section A4.
  - This script performs BH analyses using `"LR"`, `"MAST"`, and `"wilcox_limma"` within a single script. The BH procedure does not use `"LCD"`.
- `66_67_mother_code_multidecomp_realdata-AllPart12.R`
  - Code similar to `46_47_mother_code_multidecomp_realdata-AllPart12.R`, but using the covariance matrix without conditioning on the covariates. It generates the multi-decomp/e-decomp results in Section A3.

- Required datasets:
  - `HuVascAD_astrocyte.rds`
  - `HuVascAD_microglia.rds`
  - `HuVascAD_neuron.rds`
  - `HuVascAD_oligo.rds`
  - `HuVascAD_pericyte.rds`

- Analysis settings to specify:
  - The applicable settings depend on the script. In particular, `part1` and `part2` apply only to scripts whose filenames contain `AllPart12`.

  - `file_path`: The directory where the dataset is stored, excluding the dataset file name (e.g., `path/to/your/rds_folder`, not `path/to/your/rds_folder/rds_file.rds`).

  - `cell_name`: A name that specifies the dataset. One of `"astrocyte"`, `"microglia"`, `"neuron"`, `"oligo"`, and `"pericyte"`.

  - `down`: The number of cells sampled from each group (Control and AD).
    - When `down` is not `NULL`, the total number of sampled cells is `2 * down`.
    - Set to `NULL` to use all cells.

  - `max_gene`: The number of genes to be used (*p*).
    - Use `max_gene = rownames(HuVascAD) %>% length()` to use all genes.

  - `PC`: The number of latent factors.
    - For Section 5, the values of `PC` used are as follows:
      - Astrocyte: 43
      - Microglia: 11
      - Neuron: 13
      - Oligo: 64
      - Pericyte: 76

  - `llam`: A scalar that controls the lambda value for sc-softImpute.
    - It can take any nonnegative value. The final lambda is computed by multiplying this value by the largest singular value of the data matrix, i.e., `llam*softImpute::lambda0(data.matrix)`.

  - `nSim`: The number of repetitions of the statistical testing and variable selection procedures.
    - For scripts whose filenames contain `AllPart12`, the repetitions are performed in Part 2 using the saved knockoffs.
    - For the BH script, the differential analysis and BH correction are repeated without repeating preprocessing or latent factor estimation.
    - Set to `1` to reproduce the results reported in the paper.

  - `m_kos`: The number of knockoffs to generate for the multiple-knockoff procedure.

  - `testUse`: The analysis method specified in the script.
    - For scripts whose filenames contain `AllPart12`, it specifies the method used to compute test statistics: `"LCD"`, `"LR"`, `"MAST"`, or `"wilcox_limma"`.
    - For `48_mother_code_BH_realdata.R`, `testUse` is fixed to `"BH"`. This script performs BH analyses using `"LR"`, `"MAST"`, and `"wilcox_limma"`. The BH procedure does not use `"LCD"`.
    - `"wilcox_limma"` does not use covariate adjustment through the `latent.vars` argument in `FindMarkers()`.

  - `max_iter_imp`: The maximum number of iterations allowed for sc-softImpute.

  - `part1`: If `part1 = TRUE`, the code performs preprocessing and knockoff construction.

  - `part2`: If `part2 = TRUE`, the code computes test statistics and identifies significant variables.
    - If both `part1` and `part2` are `TRUE`, both parts are executed.
    - If only `part1` is `TRUE`, only Part 1 is executed.
    - If only `part2` is `TRUE`, only Part 2 is executed.
    - If both are `FALSE`, neither Part 1 nor Part 2 is executed.
    - The initialization stage is executed regardless of the values of `part1` and `part2`.
    - If knockoffs have already been generated and saved, only Part 2 needs to be executed.
    - These settings apply only to scripts whose filenames contain `AllPart12`.

- The scripts `43_44_mother_code_multiLR_realdata-AllPart12.R`, `46_47_mother_code_multidecomp_realdata-AllPart12.R`, and `66_67_mother_code_multidecomp_realdata-AllPart12.R` share the following execution structure:

  - Initialization:
    - The analysis settings are specified.
    - This stage is executed regardless of the values of `part1` and `part2`.

  - Part 1: Preprocessing and knockoff construction
    - This part is executed when `part1 = TRUE`.
    - sc-softImpute is applied to a given dataset, and knockoffs are constructed using the imputed data matrix.
    - The processed expression matrix, generated knockoffs, and relevant analysis information are saved in `First_Step_List_Knockoffs` in the directory specified by `file_path`.
    - The saved results can be reused in Part 2 without repeating preprocessing and knockoff construction.

  - Part 2: Computation of test statistics and selection of significant variables
    - This part is executed when `part2 = TRUE`.
    - First, the list saved from Part 1 is loaded into the R session.
    - Test statistics are computed according to the method specified in `testUse` and stored in `W_imp` and `W_imp_b`.
      - For `"LR"`, `"MAST"`, and `"wilcox_limma"`:
        - `W_imp` contains `-log(p-values)`.
        - `W_imp_b` contains `-log(adjusted p-values)`.
      - For `"LCD"`:
        - `W_imp` contains the test statistics computed with the `lambda.min` selected by `glmnet::cv.glmnet`.
        - `W_imp_b` contains the test statistics computed with the `lambda.1se` selected by `glmnet::cv.glmnet`.
    - The code also prints the results of the multiple-knockoff and e-BH procedures based on the target FDR specified by `target_FDR_temp` (default: 0.1).
      - `"Selected Genes (Not Selected Setting)"` represents the genes selected by multiple knockoffs using `W_imp`.
      - `"Selected Genes (Selected Setting)"` represents the genes selected by multiple knockoffs using `W_imp_b`.
      - `"Selected Genes (ebh)"` represents the genes selected by the e-BH procedure using `W_imp_b`.
    - The results, including `W_imp`, `W_imp_b`, and the selected variables, are stored in `FinalResult` and saved as an `.rda` file in the current working directory (`getwd()`).

  - Computational efficiency:
    - This structure improves computational efficiency when running on a high-performance computing platform.
    - For each knockoff construction method, Part 1 needs to be executed only once. The saved results can then be reused by running Part 2 separately for `"LCD"`, `"LR"`, `"MAST"`, and `"wilcox_limma"`, without repeating preprocessing and knockoff construction.
    - When running only Part 2, the results from Part 1 must already be saved and accessible.
    - The relevant analysis settings, including `cell_name`, `down`, `max_gene`, `PC`, `llam`, `seed_num`, `max_iter_imp`, and `m_kos`, should be consistent with those used in Part 1.
    - The same `file_path` must be used to load the saved results.

## 2. Summary script

- `99_SummaryAll_4_Real 52_realdata.R`
  - This script generates Table 1 in Section 5.
  - The current working directory must be set to the directory where the `.rda` files are stored.
    - These `.rda` files contain the results generated in the final part of:
      - `43_44_mother_code_multiLR_realdata-AllPart12.R`
      - `46_47_mother_code_multidecomp_realdata-AllPart12.R`
    - Use `setwd()` to set the current working directory.
    - Use `getwd()` to verify the current working directory.
