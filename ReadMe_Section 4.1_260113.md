# Section 4.1
#Knockoffs

- This document describes how to use the following simulation files:
  - `01_mother_code_asdp.R`
  - `02_mother_code_LR.R`
  - `03_04_mother_code_multiLR.R`
  - `05_mother_code_decomp.R`
  - `06_07_mother_code_multidecomp.R`
  - `08_mother_code_BH.R`
  - `09_mother_code_naive.R`

- It also describes the following summary files:
  - `99_SummaryAll_0_Simulation4.1.R`
  - `99_SummaryAll_0_Simulation4.1_appendix.R`

1. Simulation files
   - Each simulation file runs simulations using one of the following knockoff methods:
     - `01_mother_code_asdp.R`: *ASDP knockoffs*
     - `02_mother_code_LR.R`: *LR knockoffs*
     - `03_04_mother_code_multiLR.R`: *Multi-LR knockoffs*, e-BH procedure
     - `05_mother_code_decomp.R`: *Decomp knockoffs*
     - `06_07_mother_code_multidecomp.R`: *Multi-decomp knockoffs*, e-BH procedure
     - `08_mother_code_BH.R`: No knockoffs, BH procedure
     - `09_mother_code_naive.R`: *Standard Gaussian knockoffs*
     
   - Required datasets:
     - `combined_superior_parietal_lobe.rds`
     
   - Hyperparameters the user needs to specify:
     - `file_path`: The directory where the dataset is stored.
     - `down`: The number of cells to be used (*n*).
       - Set to `NULL` to use all cells.
     - `gene.index`: A vector specifying the genes to be used (*p*).
       - For Section 4.1, `gene.index` is `1:2000`
     - `sign_strength`: The signal strength used to create synthetic signals.  
       - In Section 4.1, this parameter is referred to as **`sgn`**.
     - `n_signal`: The number of signals in the simulation.
     - `PC`: The number of latent factors.
       - For Section 4.1, `PC` is `30`
     - `llam`: A scalar that controls the lambda value for sc-softImpute.
       - It can take any nonnegative value. The final lambda is computed by multiplying this value by the largest singular value of the data matrix. i.e., `llam*softImpute::lambda0(data.matrix)` is used.
     - `m_kos`: The number of knockoffs to generate for the multiple-knockoff procedure. (only for the multi-LR and multi-decomp codes).
     - `nSim`: The number of simulation replicates.
       - For Section 4.1, `nSim = 20`.
     - `testUse`: The method used to compute test statistics.  One of `"LCD"`, `"LR"`, `"MAST"`, or `"wilcox_limma"`. `"LCD"` cannot be used with the BH procedure. `"wilcox_limma"` does not use latent variables.
     - `max_iter_imp`: The maximum number of iterations allowed for sc-softImpute.
     
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
   - `99_SummaryAll_0_Simulation4.1.R`
     - This code reproduces Figure 1 in Section 4.1.
   - `99_SummaryAll_0_Simulation4.1_appendix.R`
     - This code reproduces Table S4 in Section B.3 of the Appendix.

   - The current working directory must be set to the directory where the `.rda` files are stored.
     - These `.rda` files contain the results generated in the final step of the simulation files described above.
     - Use `setwd()` to set the current working directory.
     - Use `getwd()` to check the current working directory.

   - Both scripts use `W_imp_b` (or `W_imp` for the BH procedure) from the `.rda` files to compute the FDP and power for the target FDR values specified in the `target_FDR_set` vector.

   - `target_FDR_set`: A vector containing target FDR values.
     - For Figure 1, `target_FDR_set = c(0.1, 0.01)` is used.
     - For Table S4, `target_FDR_set = 0.05` is used.