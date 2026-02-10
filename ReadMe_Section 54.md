# ReadMe_Section 54
#Knockoffs

- This document describes how to use the following simulation files:
  - `35_mother_code_decomp_confounding.R`
  - `38_mother_code_BH_confounding.R`

- It also describes the following summary files:
  - `99_SummaryAll_3_Simulation 54_confounding.R`

1. Simulation files
   - Each simulation file runs simulations using one of the following knockoff methods:
     - `35_mother_code_decomp_confounding.R`: *Decomp knockoffs*
     - `38_mother_code_BH_confounding.R`: No knockoffs, BH procedure

   - Required datasets:
     - `HuVascAD_seuratcluster3.rds` 

   - Hyperparameters the user needs to specify:
     - `file_path`: For this code, specify the file path directly in `HuVascAD = readRDS(file = "path/to/file.rds")`.
     - `down`: The number of cells to be used (*n*).
       - Set to `NULL` to use all cells.
     - `gene.index`: A vector specifying the genes to be used (*p*).
       - For Section 5.4, `gene.index` is `1:3000`
     - `sign_strength`: The signal strength used to create synthetic signals.  
       - In Section 5.4, this parameter is referred to as **`sgn`**.
     - `n_signal`: The number of signals in the simulation.
     - `PC`: The number of latent factors.
       - For Section 5.4, `PC` is `13`
     - `llam`: A scalar that controls the lambda value for sc-softImpute.
       - It can take any nonnegative value. The final lambda is computed by multiplying this value by the largest singular value of the data matrix. i.e., `llam*softImpute::lambda0(data.matrix)` is used.
     - `nSim`: The number of simulation replicates.
       - For Section 5.4, `nSim = 100` is used for *decomp knockoffs*, and `nSim = 20` is used for the BH procedure.
     - `testUse`: The method used to compute test statistics.  One of `"LCD"`, `"LR"`, `"MAST"`, or `"wilcox_limma"`. `"LCD"` cannot be used with the BH procedure. `"wilcox_limma"` does not use latent variables. 
     - `max_iter_imp`: The maximum number of iterations allowed for sc-softImpute.
     - The following hyperparameters control which variables are used for the imputation step and the computation of test statistics:
       - `batch_option`: One of `"with"`, `"without"`, or `"permute"`.
         - Set `batch_option = "with"` to use the batch variable in both the imputation and computation steps.
         - Set `batch_option = "without"` to exclude the batch variable.
         - Set `batch_option = "permute"` to permute the batch variable.
       - `"CDR_impute_option"`: One of `"with"` or `"without"`.
         - Set `CDR_impute_option = "with"` to use the CDR variable in the imputation step.
         - Set `CDR_impute_option = "without"` to exclude the CDR variable in the imputation step.
       - `"CDR_p_compute_option"`: One of `"with"` or `"without"`.
         - Set `CDR_p_commpute_option = "with"` to use the CDR variable in the computation step.
         - Set `CDR_p_commpute_option = "without"` to exclude the CDR variable in the computation step.

   - Each simulation file consists of two main parts:
     1. Imputation and knockoff construction
        - In this part, sc-softImpute is applied to a given dataset and knockoffs are constructed using the imputed data matrix.
        - For the BH procedure, there are no imputation or knockoff construction steps.
       
     2. Computation of test statistics and applying the knockoff procedure
        - Test statistics are computed according to the method specified in `testUse` and stored in `W_imp` and `W_imp_b`.
          - For `"LR"`, `"MAST"`, and `"wilcox_limma"`:
            - `W_imp` contains `-log(p-values)`
            - `W_imp_b` contains `-log(adjusted p-values)`
          - For `LCD`:
            - `W_imp` contains the test statistics computed with the `lambda.min` selected by `glmnet::cv.glmnet`.
            - `W_imp_b` contains the test statistics computed with the `lambda.1se` selected by `glmnet::cv.glmnet`.
        
        - The code also prints out the results of the knockoffs and the e-BH procedures (only for multiple knockoffs) based on the target FDR equal to 0.1. 
          - `"Selected Genes (Not Selected Setting)"` represents the genes selected by multiple knockoffs using `W_imp`.
          - `"Selected Genes (Selected Setting)"` represents the genes selected by multiple knockoffs using `W_imp_b`.
          - `"Selected Genes (ebh)"` represents the genes selected by the e-BH procedure using `W_imp_b`. (only for multiple knockoffs)
          
        - The list containing `W_imp` and `W_imp_b` is saved in the current directory (`getwd()`).

2. Summary files
   - `99_SummaryAll_3_Simulation 54_confounding.R`
     - This code reproduces Table 2 in Section 5.4 and Table S4 in Appendix.
     - This code consists of two parts:
       - Part for Table 2 (*decomp knockoffs*)
       - Part for Table S4 (BH procedure)

   - The current working directory must be set to the directory where the `.rda` files are stored.
     - For the Table 2 and Table S4 parts, the working directory must be set differently.
       - For the first part, set the working directory to the folder containing the `.rda` files generated by the *decomp knockoffs*.
       - For the second part, set it to the folder containing the `.rda` files generated by the BH procedure.
     - Use `setwd()` to set the current working directory.
     - Use `getwd()` to check the current working directory.

   - The script uses `W_imp_b` (or `W_imp` for the BH procedure) from the `.rda` files to compute the FDP and power for the target FDR value specified by `target_FDR_set`.

   - `target_FDR_set`: A scalar specifying the target FDR value.
     - For Tables 2 and S4, `target_FDR_set = 0.1`.

   - For each method and setting, the average FDP, the standard deviation of the FDP, the average power, and the standard deviation of the power are stored in four separate `.csv` files:
     - `Section54_Table2_LR_DownNULL_Avg_FDP.csv`: Average FDP
     - `Section54_Table2_LR_DownNULL_SD_FDP.csv`: SD of FDP
     - `Section54_Table2_LR_DownNULL_Avg_Power.csv`: Average power
     - `Section54_Table2_LR_DownNULL_SD_Power.csv`: SD of power
