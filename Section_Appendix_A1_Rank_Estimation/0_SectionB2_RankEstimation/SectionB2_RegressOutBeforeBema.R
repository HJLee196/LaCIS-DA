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

superior_parietal_lobe = readRDS(file = 'path/to/combined_superior_parietal_lobe.rds')
superior_parietal_lobe = UpdateSeuratObject(superior_parietal_lobe) 
superior_parietal_lobe$sex <- ifelse(superior_parietal_lobe$sex=='female',1,0)

####
# 1-3. Set hyperparameters
####

seed_num = 2017
set.seed(seed_num)
down = NULL
gene.index = 1:2000

####
# 1-4. Downsample the data
####

if(is.null(down)){
  test_data <- superior_parietal_lobe # use full dataset
} else {
  # Manually downsample the data:
  healthy_M <- names(superior_parietal_lobe@active.ident)[superior_parietal_lobe@active.ident %in% "healthy_Microglia"] 
  ad1_M <- names(superior_parietal_lobe@active.ident)[superior_parietal_lobe@active.ident %in% "ad1_Microglia"] 
  ad2_M <- names(superior_parietal_lobe@active.ident)[superior_parietal_lobe@active.ident %in% "ad2_Microglia"]

  down_ind <- c(sample(healthy_M, size = down),
                sample(ad1_M, size = down),
                sample(ad2_M, size = down))
  # which(names(superior_parietal_lobe@active.ident) %in% down_ind)
  
  test_data <- subset(superior_parietal_lobe, 
                      cells = which(names(superior_parietal_lobe@active.ident) %in% down_ind))
}

rm(superior_parietal_lobe)


# Take a subset of genes
feature.names <- rownames(test_data)[gene.index]

xp_data.sub <- subset(x = test_data, features = feature.names)

# Generate CDR and sex from the data matrix
xp_data.matrix = GetAssayData(object = xp_data.sub, layer = "data") %>% t() %>% as(Class = "dgCMatrix")
xp_data.matrix.count = GetAssayData(object = xp_data.sub, layer = "count")
xp_data.matrix.exp = t(xp_data.matrix.count != 0) %>% as(Class = "dgCMatrix")
CDR = xp_data.matrix.exp %>% rowMeans()
CDR_sex_matrix = cbind(CDR, xp_data.sub$sex)

# Filter the genes that are expressed in more than 10%. 
xp_data.matrix.exp.count <- colSums(xp_data.matrix.exp)
threshold_for_genes = 0.01
feature_ind = (colMeans(xp_data.matrix.exp) >= threshold_for_genes)
filtered_feature_names = colnames(xp_data.matrix.exp)[feature_ind]

xp_data.sub <- subset(x = test_data, features = filtered_feature_names)

xp_data.matrix = GetAssayData(object = xp_data.sub, layer = "data") %>% t() %>% as(Class = "dgCMatrix")
xp_data.matrix.count = GetAssayData(object = xp_data.sub, layer = "count")
xp_data.matrix.exp = t(xp_data.matrix.count != 0) %>% as(Class = "dgCMatrix")
xp_data.matrix.exp.count <- colSums(xp_data.matrix.exp)

print("Min Count: ")
xp_data.matrix.exp.count %>% min() %>% print()
print("0.01 X nrow(xp_data.matrix)")
(0.01*(xp_data.matrix %>% nrow())) %>% print()

# Center the data matrix and the covariate matrix

xp_data.matrix.fix <- (colSums(xp_data.matrix, na.rm = T))/xp_data.matrix.exp.count
xp_data.matrix <- sweep(xp_data.matrix,2,xp_data.matrix.fix,FUN = "-")*(xp_data.matrix.exp)
xp_data.matrix %>% colSums() %>% abs() %>% max()

# Regress Out CDR and sex from the data matrix
X_0 = CDR_sex_matrix
X_0 <- sweep(X_0,2,colMeans(X_0),FUN = "-")
B_0 = solve(t(X_0)%*%X_0)%*%t(X_0)%*%xp_data.matrix
xp_data.matrix_regressed_new = (xp_data.matrix - X_0 %*% B_0) %>% as.matrix()
xp_data.matrix_regressed_new %>% colMeans() %>% abs() %>% max()

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
  paste("Section52_xp_data_down", 
        ifelse(is.null(down), yes = "NULL", no = down),
        "_regressed_new_p",
        round(threshold_for_genes*100), 
        ".csv", sep = "")

print(dim(xp_data.matrix_regressed_new))

rda_name =
  data_name %>%
  str_replace(pattern = ".csv$", replacement = ".rda")

save(xp_data.matrix_regressed_new, file = rda_name)
