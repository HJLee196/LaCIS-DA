rm(list=ls())

###
# 1. Preparation
###

###
# 1-1. Load Necessary R Packages
###

library(tidyverse)
library(Seurat)
library(knockoff)
library(MAST)
library(Matrix)
library(softImpute)
library(glmnet)
library(mvnfast) # create.gaussian

####
# 1-2. Read the dataset
#### 

HuVascAD = readRDS(file = 'path/to/HuVascAD_seuratcluster3.rds')
HuVascAD = subset(HuVascAD, cells = which(HuVascAD$region %in% "Hippocampus"))
table_batch_sample = table(HuVascAD$sample, HuVascAD$batch)

####
# 1-3. Set hyperparameters
####

seed_num = 2017
set.seed(seed_num)
down = NULL
gene.index = 1:3000

# Prespecify how many variables are fixed during the imputation
batch_option = "with" # "with" "without" "permute"
CDR_impute_option = "with"
CDR_p_compute_option = "with"

impute_cov_names = c("CDR", "age", "age2", "gender", "batch")
covariate_names = c("CDR", "age", "age2", "gender", "batch")

if (batch_option == "permute") {
  batch_permuted = sample(HuVascAD$batch)
  names(batch_permuted) = names(HuVascAD$batch)
  HuVascAD$batch = batch_permuted
  
} else if (batch_option == "without") {
  impute_cov_names = impute_cov_names[impute_cov_names != "batch"]
  covariate_names = covariate_names[covariate_names != "batch"]
}

if (CDR_impute_option == "without") {
  impute_cov_names = impute_cov_names[impute_cov_names != "CDR"]
}

if (CDR_p_compute_option == "without") {
  covariate_names = covariate_names[covariate_names != "CDR"] 
}

cat("batch option:", batch_option, "\n")
cat("CDR impute option:", CDR_impute_option, "\n")
print(impute_cov_names)
cat("CDR p compute option:", CDR_p_compute_option, "\n")
print(covariate_names)

k_num = 3 + length(unique(HuVascAD$gender)) - 1 + length(unique(HuVascAD$batch)) - 1 # 7

####
# 1-4. Downsample the data
####

if(is.null(down)){
  test_data <- HuVascAD # use full dataset
} else {
  # Manually downsample the data:
  healthy_M <- names(HuVascAD$treat)[HuVascAD$treat %in% "Control"] 
  ad_M <- names(HuVascAD$treat)[HuVascAD$treat %in% "AD"] # length(ad_M) + length(healthy_M) 12193
  
  down_ind <- c(sample(healthy_M,down),
                sample(ad_M,down))
  
  test_data <- subset(HuVascAD, 
                      cells = which(names(HuVascAD$treat) %in% down_ind))
}
rm(HuVascAD) # rm dataset to save memory

# Take a subset of genes
feature.names <- rownames(test_data)[gene.index]
xp_data.sub <- subset(x = test_data, features = feature.names)

####
# 1-5. Create the X variables
####

# Calculate CDR ####
xp_data.count <- GetAssayData(object = xp_data.sub, layer = "count") %>% as(Class = "dgCMatrix")
# rotate to match downstream analyses.
xp_data.matrix.exp <- t(as.matrix(xp_data.count) > 0) %>% as(Class = "dgCMatrix") # indicator matrix for expressed genes, 
xp_data.matrix.exp.count <- colSums(xp_data.matrix.exp)

CDR = rowMeans(xp_data.matrix.exp) # expressed genes in each CELL

# xp_data <- GetAssayData(object = xp_data.sub, layer = "data") 
xp_data.matrix <- t(GetAssayData(object = xp_data.sub, layer = "data")) %>% as(Class = "dgCMatrix")# Normalized data matrix

print(dim(xp_data.matrix))
print(dim(xp_data.matrix.exp))

# Include the covariates
## 1. gender
xp_data.sub$gender = xp_data.sub$gender %>% as.factor()
xp_data.sub$gender %>% table() %>% barplot()
## 2. CDR
xp_data.sub$CDR = CDR
## 3. Sample ID
xp_data.sub$sample = xp_data.sub$sample %>% as.factor()
xp_data.sub$sample %>% table() %>% barplot()
## 4. batch label
xp_data.sub$batch = xp_data.sub$batch %>% as.factor()
xp_data.sub$batch %>% table() %>% barplot()
## 5. age
xp_data.sub$age %>% table() %>% barplot()
## 6. age squared 
xp_data.sub$age2 = xp_data.sub$age^2
xp_data.sub$age2 %>% table() %>% barplot() 

# xp_data.sub$sample[xp_data.sub$batch == 4]

# Sample, Batch, Gender are categorical variables
# Sample and Batch have more than 2 categories.
samp = model.matrix(y~x-1,data.frame(x = xp_data.sub$sample, y = 1))
samp = samp[,-1] # To set the first group to the reference group
colnames(samp) = 
  colnames(samp) %>%
  str_replace(pattern = "^x", replacement = "sample_")

batch = model.matrix(y~x-1,data.frame(x = xp_data.sub$batch, y = 1))
batch = batch[,-1] # To set the first group to the reference group
colnames(batch) = 
  colnames(batch) %>%
  str_replace(pattern = "^x", replacement = "batch_")

# Generate the covariate matrix
# Standardize the age variable
covariate_mat = 
  cbind(xp_data.sub$CDR, 
        xp_data.sub$age,
        as.numeric(xp_data.sub$gender),
        samp, 
        batch)
colnames(covariate_mat)[1:4] = c("CDR", "age", "age2", "gender")

covariate_imp_pattern = impute_cov_names %>% str_c(collapse = "|") %>% str_replace(pattern = "(^.*$)", replacement = "\\^(\\1)")

selected_covariate_imp_ind = 
  colnames(covariate_mat) %>%
  str_starts(pattern = covariate_imp_pattern)
selected_covariate_imp = colnames(covariate_mat)[selected_covariate_imp_ind]

covariate_p_pattern = covariate_names %>% str_c(collapse = "|") %>% str_replace(pattern = "(^.*$)", replacement = "\\^(\\1)")

selected_covariate_p_ind = 
  colnames(covariate_mat) %>%
  str_starts(pattern = covariate_p_pattern)
selected_covariate_p = colnames(covariate_mat)[selected_covariate_p_ind]

covariate_imp_mat = covariate_mat[,selected_covariate_imp]
covariate_p_mat = covariate_mat[,selected_covariate_p]

# Filter the genes that are expressed in more than 10%. 
xp_data.matrix.exp.count <- colSums(xp_data.matrix.exp)
threshold_for_genes = 0.01
feature_ind = (colMeans(xp_data.matrix.exp) >= threshold_for_genes)
filtered_feature_names = colnames(xp_data.matrix.exp)[feature_ind]

xp_data.sub <- subset(x = test_data, features = filtered_feature_names)

xp_data.matrix = GetAssayData(object = xp_data.sub, layer = "data") %>% t() %>% as(Class = "dgCMatrix")
xp_data.matrix.count = GetAssayData(object = xp_data.sub, layer = "count") %>% as(Class = "dgCMatrix")
xp_data.matrix.exp = t(xp_data.matrix.count != 0) %>% as(Class = "dgCMatrix")
xp_data.matrix.exp.count <- colSums(xp_data.matrix.exp)

print("Min Count: ")
xp_data.matrix.exp.count %>% min() %>% print()
print("0.01 X nrow(xp_data.matrix)")
(0.01*(xp_data.matrix %>% nrow())) %>% print()

# Center the data matrix and the covariate matrix

xp_data.matrix.fix <- (colSums(xp_data.matrix, na.rm = T))/xp_data.matrix.exp.count
xp_data.matrix <- sweep(xp_data.matrix,2,xp_data.matrix.fix,FUN = "-")*(xp_data.matrix.exp)
xp_data.matrix %>% apply(2, sum) %>% abs() %>% max() %>% print()

# Regress Out CDR and sex from the data matrix
X_0 = covariate_p_mat
X_0 <- sweep(X_0,2,colMeans(X_0),FUN = "-")
B_0 = solve(t(X_0)%*%X_0)%*%t(X_0)%*%xp_data.matrix
xp_data.matrix_regressed_new = (xp_data.matrix - X_0 %*% B_0) %>% as.matrix()
xp_data.matrix_regressed_new %>% apply(2, mean) %>% abs() %>% max()

xp_data.matrix_regressed_new[,1][xp_data.matrix_regressed_new[,1] != 0 ] %>% hist()
xp_data.matrix_regressed_new[,1][xp_data.matrix_regressed_new[,1] != 0 ] %>% sd()
xp_data.matrix_regressed_new %>% 
  # apply(MARGIN = 2, 
  #       FUN = function(xp_row){
  #         xp_row[xp_row != 0] %>% sd()
  #       }) %>% 
  colSds() %>% # this is different from the apply function above. For simplicity and memory issue, I use this instead
  hist()

data_name = 
  paste("Section53_xp_data_down", 
        ifelse(is.null(down), yes = "NULL", no = down),
        "_regressed_new_p",
        round(threshold_for_genes*100), 
        ".csv", sep = "")

print(dim(xp_data.matrix_regressed_new))

rda_name =
  data_name %>%
  str_replace(pattern = ".csv$", replacement = ".rda")

save(xp_data.matrix_regressed_new, file = rda_name)

