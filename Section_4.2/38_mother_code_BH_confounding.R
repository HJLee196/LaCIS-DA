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

lm_impute = function(Yl_imp,Xl,Bl,Yl,Yl_exp,PC,q_ncol,lambda){
  
  ## Fit linear model 
  n <- dim(Yl_imp)[1]
  p <- dim(Yl_imp)[2]
  Q_index = 1:q_ncol
  A_index = 1:PC +q_ncol
  
  # Fix Bl
  B_q = Bl[Q_index,]
  B_A = Bl[A_index,]
  Q = Xl[,Q_index]
  
  Y_hat = Xl%*%Bl
  Yl_imp = (Yl - Y_hat)*Yl_exp + Y_hat - Q%*%(B_q)
  I = diag(1,PC,PC)
  Al = solve((B_A)%*%t(B_A) + lambda*I)%*%(B_A)%*%t(Yl_imp)
  Al = t(Al)
  
  # Fix Al
  Xl[,A_index] = Al
  
  Y_hat = Xl%*%Bl
  Yl_imp = (Yl - Y_hat)*Yl_exp + Y_hat
  I = diag(1,PC+q_ncol,PC+q_ncol)
  Bl = solve(t(Xl)%*%Xl + lambda*I)%*%t(Xl)%*%Yl_imp
  
  Y_hat = Xl%*%Bl
  Yl_imp = (Yl - Y_hat)*Yl_exp + Y_hat
  
  return(list(Yl_imp,Xl,Bl))
}



####
# 1-2. Read the dataset
#### 

HuVascAD = readRDS(file = 'path/to/HuVascAD_seuratcluster3.rds')
# HuVascAD = readRDS(file = './HuVascAD_seuratcluster3.rds')
HuVascAD = subset(HuVascAD, cells = which(HuVascAD$region %in% "Hippocampus"))
table_batch_sample = table(HuVascAD$sample, HuVascAD$batch)

####
# 1-3. Set hyperparameters
####

seed_num = 2017
set.seed(seed_num)
down = 3000
gene.index = 1:3000
sign_strength = 3
n_signal = 50
PC = 13
llam = 0.1 
nSim = 1

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

k_num = 3 + length(unique(HuVascAD$gender)) - 1 + length(unique(HuVascAD$batch)) - 1

testUse = "LR" # c("wilcox_limma", "LR", "MAST", "LCD")
LCD_crit_use = "lambda.1se"
max_iter_imp = 100

hyperparam_used_list = 
  list(seed_num = seed_num,
       down = down,
       max_gene.index = max(gene.index),
       sign_strength = sign_strength,
       n_signal = n_signal,
       PC = PC,
       llam = llam,
       nSim = nSim,
       k_num = k_num,
       max_iter_imp = max_iter_imp,
       testUse = testUse,
       covariate_names = covariate_names,
       impute_cov_names = impute_cov_names
  )

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
xp_data <- GetAssayData(object = xp_data.sub, layer = "count")

xp_data.matrix.exp <- t(as.matrix(xp_data) > 0) # indicator matrix for expressed genes, 
xp_data.matrix.exp.count <- apply(xp_data.matrix.exp,2,sum)

xp_data.matrix <- (as.matrix(xp_data)) # do not transpose here. Keep the genes as rows.

# exclude empty variables
feature.exp.ind <- which(xp_data.matrix.exp.count > (PC+k_num+1)) # PC + num of X variables + intercept
feature.names.new <- feature.names[feature.exp.ind]  # length(feature.names.new)

xp_data.matrix.exp <- xp_data.matrix.exp[,feature.exp.ind]
xp_data.matrix.exp.count <- xp_data.matrix.exp.count[feature.exp.ind]
CDR = rowMeans(xp_data.matrix.exp) # expressed genes in each CELL

# refresh subset
xp_data.sub <- subset(x = test_data, features = feature.names.new)
xp_data <- GetAssayData(object = xp_data.sub, layer = "data") # Normalized data matrix

xp_data.matrix <- t(as.matrix(xp_data))

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
## 6. age squared (dont use as written in Section 5.3)
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

####
# 2. Run simulations
####

if (testUse == "wilcox_limma"){
  covariate_names = NULL
}

# hyperparam_used_list$imp_total_seconds = imp_total_seconds
# hyperparam_used_list$KO_total_seconds = KO_total_seconds

xp_data.matrix = xp_data.matrix %>% as(Class = "dgCMatrix")
# xp_data.knockoff.sav = xp_data.knockoff.sav %>% as(Class = "dgCMatrix")

FinalResult = list()

start_time_nsim_set = NULL
end_time_nsim_set = NULL

for(nsim in 1:nSim){
  # Computing Time
  start_time_nsim = Sys.time()
  start_time_nsim_set = c(start_time_nsim_set, start_time_nsim)
  
  n <- dim(xp_data.matrix)[1]
  p <- dim(xp_data.matrix)[2]

  set.seed(seed = seed_num + nsim)
  # generate coefficients ####
  
  # Scaling Method 1: Simple
  
  xp_data.matrix_sd = xp_data.matrix %>% apply(2, sd)
  xp_data.matrix_scale = 
    xp_data.matrix %>%
    sweep(MARGIN = 2, 
          STATS = xp_data.matrix_sd, 
          FUN = "/")
  
  norm_coef <- rep(0, p)
  tmp_coef <- rep(0,n_signal)
  
  for(i in 1:n_signal){
    new_coef <- 0
    # truncate the coefficients which are too small
    while(isTRUE((new_coef)^2 < (sign_strength^2)*2*log(p)/n)){
      new_coef <- rnorm(1,0,sign_strength*sqrt(2*log(p)/n))
    }
    
    tmp_coef[i] <- new_coef
  }
  
  ######################
  # Filtering Scheme 1 #
  ######################
  # xp_data.matrix_expressed_prop = xp_data.matrix_denominator/n
  # filtered_idx = ((xp_data.matrix_expressed_prop > 0.1) & (xp_data.matrix_expressed_prop < 0.3)) %>% which()
  filtered_idx = 1:p
  
  sig_coef <- sample(filtered_idx, n_signal)
  norm_coef[sig_coef] <- tmp_coef
  
  
  # introduce the logit link
  label_p = 
    (1/(1 + exp(- (xp_data.matrix_scale %*% norm_coef + 2*sign_strength*sqrt(2*log(p)/n)*as.numeric(xp_data.sub$batch)))
    )) %>% 
    as.vector()
  rm(xp_data.matrix_scale)
  
  # Generate labels
  label_binom <- rbinom(n = length(label_p), size = 1, prob = label_p)
  label_ad <- ifelse((label_binom == 1), yes = "ad", no = "healthy")
  
  xp_data.sub <- SetIdent(object = xp_data.sub, value = label_ad)
  
  # Calculate importance statistics
  
  ## LRT, MAST, WRT ##
  if (testUse != "LCD")
  {
    # use observed values to find p-values
    # shouldn't matter which slot I change as long as I use it. # layer = slot Seurat v5 use layer instead but FindMarkers use slot for the both.
    xp_data.sub <- SetAssayData(object = xp_data.sub, layer = "data", new.data = t(xp_data.matrix))
    
    result.sub <- FindMarkers(xp_data.sub, ident.1 = 'healthy', ident.2 = 'ad', slot = "data",
                              min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use = testUse, latent.vars = covariate_names) 
  
    # P-values
    W_imp = result.sub$p_val[match(feature.names.new, rownames(result.sub))]
    
    W_imp_b = result.sub$p_val_adj[match(feature.names.new, rownames(result.sub))]
    
    cat("####", testUse, "####", "\n")
    print(covariate_names)
  } else if (testUse == "LCD") {
    # No LCD for the BH procedure
  }
  
  # calculate fdr and power (for testing purpose)
  
  cat("# Iter: ", nsim, "\n")
  
  if (sum(is.na(W_imp)) != 0 | sum(is.na(W_imp_b)) != 0) {
    # BC knockoffs have negative values
    cat("# na in W_imp: ", sum(is.na(W_imp)), "\n")
    cat("# na in W_imp_b: ", sum(is.na(W_imp_b)), "\n")
    cat("Set it zero", "\n")
    W_imp[is.na(W_imp)] = 0
    W_imp_b[is.na(W_imp_b)] = 0
  }
  
  if (testUse != "LCD") 
  {
    BH_p = p.adjust(p = W_imp, method = "BH")
    selected = (1:length(BH_p))[BH_p < 0.1]
    
    print((length(selected) - sum(selected %in% sig_coef))/length(selected)) # FDR
    print(sum(selected %in% sig_coef)/length(sig_coef)) # Power
  } else if (testUse == "LCD") {
    W_imp = NULL
    selected = NULL
  }
  
  BH_p_b = p.adjust(p = W_imp_b, method = "BH")
  selected_b = (1:length(BH_p_b))[BH_p_b < 0.1]
  
  print((length(selected_b) - sum(selected_b %in% sig_coef))/length(selected_b)) # FDR
  
  print(sum(selected_b %in% sig_coef)/length(sig_coef)) # Power
  
  # Computation Time
  end_time_nsim = Sys.time()
  end_time_nsim_set = c(end_time_nsim_set, end_time_nsim)
  
  nsim_total_seconds = as.numeric(difftime(end_time_nsim , start_time_nsim, units = "secs"))
  nsim_minutes <- floor(nsim_total_seconds / 60)
  nsim_seconds <- round(nsim_total_seconds %% 60)
  
  print(sprintf("Spent Computation Time (each sim): %02d:%02d", nsim_minutes, nsim_seconds))
  
  FinalResult[[nsim]] <- (list(sig_coef = sig_coef, 
                               W_imp = W_imp, 
                               selected = selected, 
                               W_imp_b = W_imp_b, 
                               selected_b = selected_b, 
                               feature.names.new = feature.names.new,
                               hyperparam_used_list = hyperparam_used_list,
                               nsim_total_seconds = nsim_total_seconds))
}

start_time_allsim = min(start_time_nsim_set)
end_time_allsim = max(end_time_nsim_set)

allsim_total_seconds = as.numeric(difftime(end_time_allsim , start_time_allsim, units = "secs"))
allsim_hours <- floor(allsim_total_seconds / 3600)
allsim_minutes <- floor((allsim_total_seconds %% 3600) / 60)
allsim_seconds <- round(allsim_total_seconds %% 60)
print(sprintf("Spent Computation Time (all sims): %02d:%02d:%02d", allsim_hours, allsim_minutes, allsim_seconds))

####
# 3. Save the results
####

data_name = "HuVascAD"
method_name = "BH"

if ((LCD_crit_use == "lambda.1se")&(testUse == "LCD"))
{
  testUseFull = str_c(testUse, "1")
} else {
  testUseFull = testUse
}

# batch_option # "with" "without" "permute"
# CDR_impute_option # "with" "without"
# CDR_p_compute_option # "with" "without"

file_name = 
  str_c(data_name,
        "_method", method_name,
        "_testUse", testUseFull,
        "_down", ifelse(is.null(down), yes = "NULL", no = down),
        "_gene", max(gene.index),
        "_sgn", sign_strength,
        "_PC", PC,
        "_10llam", 10*llam,
        "_batch", batch_option,
        "_impCDR", CDR_impute_option,
        "_compCDR", CDR_p_compute_option,
        "_nSim", nSim, 
        "_seed", seed_num,
        "_maxIterImp", max_iter_imp,
        ".rda")

save(FinalResult,
     file = file_name)

