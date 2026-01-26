# ReadMe_Section 5 - Real Data Analysis Code
#Knockoffs

- This document describes how to use the following simulation files:
  - `43_44_mother_code_multiLR_realdata-AllPart12.R`
  - `46_47_mother_code_multidecomp_realdata-AllPart12.R`
  - `48_mother_code_BH_realdata.R` 

- It also describes the following summary files:
  - `99_SummaryAll_4_Simulation 5_realdata.R`.

1. `43_44_mother_code_multiLR_realdata-AllPart12.R`,`46_47_mother_code_multidecomp_realdata-AllPart12.R`, and `48_mother_code_BH_realdata.R`
   - `43_44_mother_code_multiLR_realdata-AllPart12.R`
     - Code for reproducing the multi-LR/e-LR results in the real data analysis in Section 5.
   - `46_47_mother_code_multidecomp_realdata-AllPart12.R`
     - Code for reproducing the multi-decomp/e-decomp results in the real data analysis in Section 5.
   - `48_mother_code_BH_realdata.R`
     - Code for reproducing the BH results in the real data analysis in Section 5.
   - `66_67_mother_code_multidecomp_realdata-AllPart12.R`
     - Code similar to `46_47_mother_code_multidecomp_realdata-AllPart12.R`, but using the covariance matrix without conditioning on the covariates.
     
   - Required datasets:
     - `HuVascAD_astrocyte.rds`
     - `HuVascAD_microglia.rds`
     - `HuVascAD_neuron.rds`
     - `HuVascAD_oligo.rds`
     - `HuVascAD_pericyte.rds`
     
   - Hyperparameters the user needs to specify:
     - `file_path`: The directory where the dataset is stored, excluding the dataset file name (e.g., `path/to/your/rds_folder`, not `path/to/your/rds_folder/rds_file.rds`).
     - `cell_name`: A name that specifies the dataset. One of `"astrocyte"`, `"microglia"`, `"neuron"`, `"oligo"`, and `"pericyte"`.
     - `down`: The number of cells to be used (*n*).
       - Set to `NULL` to use all cells.
     - `max_gene`: The number of genes to be used (*p*).
       - Use `max_gene = rownames(HuVascAD) %>% length()` to use all genes.
     - `PC`: The number of latent factors.
       - For Section 5, PCs used are as follows:
         - Astrocyte:  43
         - Microglia: 11
         - Neuron: 13
         - Oligo: 64
         - Pericyte: 76
     - `llam`: A scalar that controls the lambda value for sc-softImpute.
       - It can take any nonnegative value. The final lambda is computed by multiplying this value by the largest singular value of the data matrix. i.e., `llam*softImpute::lambda0(data.matrix)` is used.
     - `m_kos`: The number of knockoffs to generate for the multiple-knockoff procedure.
     - `testUse`: The method used to compute test statistics.  One of `"LCD"`, `"LR"`, `"MAST"`, or `"wilcox_limma"`. `"LCD"` cannot be used with the BH procedure. `"wilcox_limma"` does not use latent variables.
     - `max_iter_imp`: The maximum number of iterations allowed for sc-softImpute.
     
   - The scripts `43_44_mother_code_multiLR_realdata-AllPart12.R` and `46_47_mother_code_multidecomp_realdata-AllPart12.R` consist of two parts:
     1. Imputation and knockoff construction
        - In this part, sc-softImpute is applied to a given dataset and knockoffs are constructed using the imputed data matrix.
        - The list containing the resulting data matrix is saved in the directory specified by`file_path`.
        
     2. Computation of test statistics and applying the knockoff procedure
        - First, the list saved from the first part is loaded into the R session.
        
        - Test statistics are computed according to the method specified in `testUse` and stored in `W_imp` and `W_imp_b`.
          - For `"LR"`, `"MAST"`, and `"wilcox_limma"`:
            - `W_imp` contains `-log(p-values)`
            - `W_imp_b` contains `-log(adjusted p-values)`
          - For `LCD`:
            - `W_imp` contains the test statistics computed with the `lambda.min` selected by `glmnet::cv.glmnet`.
            - `W_imp_b` contains the test statistics computed with the `lambda.1se` selected by `glmnet::cv.glmnet`.
            
        - The code also prints out the results of the multiple knockoffs and e-BH procedures based on the target FDR specified in `target_FDR_temp` (default: 0.1). 
          - `"Selected Genes (Not Selected Setting)"` represents the genes selected by multiple knockoffs using `W_imp`.
          - `"Selected Genes (Selected Setting)"` represents the genes selected by multiple knockoffs using `W_imp_b`.
          - `"Selected Genes (ebh)"` represents the genes selected by the e-BH procedure using `W_imp_b`.
          
        - The list containing `W_imp` and `W_imp_b` is saved in the current directory (`getwd()`).
        
     - Each script is split into two parts for computational efficiency when running on a high-performance computing platform. For each knockoff construction method, knockoffs need to be created only once to compute`"LCD"`, `"LR"`, `"MAST"`, and `"wilcox_limma"`.

2. `99_SummaryAll_4_Simulation 5_realdata.R`
   - This code produces Table 3 in Section 5.
   - The current working directory must be set to the directory where the `.rda` files are stored.
     - These `.rda` files are the results generated in the final part of 
       `43_44_mother_code_multiLR_realdata-AllPart12.R` and 
       `46_47_mother_code_multidecomp_realdata-AllPart12.R`.
     - Use `setwd()` to set the current working directory.
     - Use `getwd()` to check the current working directory.   