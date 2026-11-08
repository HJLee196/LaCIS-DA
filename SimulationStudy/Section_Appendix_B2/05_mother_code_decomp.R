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

# superior_parietal_lobe = readRDS(file = 'path/to/combined_superior_parietal_lobe.rds')
superior_parietal_lobe = readRDS(file = './combined_superior_parietal_lobe.rds')
superior_parietal_lobe = UpdateSeuratObject(superior_parietal_lobe) 
superior_parietal_lobe$sex <- ifelse(superior_parietal_lobe$sex=='female',1,0)

####
# 1-3. Set hyperparameters
####

seed_num = 2017
set.seed(seed_num)
down = 500
gene.index = 1:2000
sign_strength = 3
n_signal = 50
PC = 30
llam = 0.1 # max_singular X llam (0.1 in the original one)
nSim = 1
testUse = "LCD" # c("wilcox_limma", "LR", "MAST", "LCD")
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
       max_iter_imp = max_iter_imp,
       testUse = testUse)

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
  
  test_data <- subset(superior_parietal_lobe, 
                      cells = which(names(superior_parietal_lobe@active.ident) %in% down_ind))
}

rm(superior_parietal_lobe)

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
feature.exp.ind <- which(xp_data.matrix.exp.count > (PC+2+1)) # PC + num of X variables + intercept
feature.names.new <- feature.names[feature.exp.ind]  # length(feature.names.new)

xp_data.matrix.exp <- xp_data.matrix.exp[,feature.exp.ind]
xp_data.matrix.exp.count <- xp_data.matrix.exp.count[feature.exp.ind]
CDR = rowMeans(xp_data.matrix.exp) # expressed genes in each CELL

# refresh subset
xp_data.sub <- subset(x = test_data, features = feature.names.new)
xp_data <- GetAssayData(object = xp_data.sub, layer = "data") # Normalized data matrix

xp_data.matrix <- t(as.matrix(xp_data))

xp_data.sub$CDR = CDR

covariate_mat = cbind((xp_data.sub$sex),CDR)
covariate_names = c("sex", "CDR")
colnames(covariate_mat) = covariate_names

####
# 2. Impute missing components
####

# Center target matrix, but only center the expressed parts.
xp_data.matrix.fix <- (colSums(xp_data.matrix, na.rm = T))/xp_data.matrix.exp.count
xp_data.matrix <- sweep(xp_data.matrix,2,xp_data.matrix.fix,FUN = "-")*(xp_data.matrix.exp)

xp_data.matrix.imp = xp_data.matrix

tmp = RunPCA(t(xp_data.matrix.imp), npcs = PC, verbose = F) # slot = "data", ???
X_0 = tmp@cell.embeddings
rm(tmp)

X_0 <- cbind(covariate_mat,X_0)
X_0 <- sweep(X_0,2,colMeans(X_0),FUN = "-")
B_0 = solve(t(X_0)%*%X_0)%*%t(X_0)%*%xp_data.matrix.imp

start_time_imp = Sys.time()

lambda = lambda0(xp_data.matrix.imp)*llam 
print(paste0("lambda: ",round(lambda/llam,2)))

xp_data.matrix.imp_l <- lm_impute(Yl_imp = xp_data.matrix.imp,
                                  Yl = xp_data.matrix,
                                  Yl_exp = xp_data.matrix.exp,
                                  Bl = B_0,
                                  Xl = X_0,
                                  PC = PC,
                                  q_ncol = 2,
                                  lambda = lambda)

xp_data.matrix.imp_l <- lm_impute(Yl_imp = xp_data.matrix.imp_l[[1]],
                                  Yl = xp_data.matrix,
                                  Yl_exp = xp_data.matrix.exp,
                                  Bl = xp_data.matrix.imp_l[[3]],
                                  Xl = xp_data.matrix.imp_l[[2]],
                                  PC = PC,
                                  q_ncol = 2,
                                  lambda = lambda)


stop_crit_compare <- sum((xp_data.matrix.imp_l[[1]] - xp_data.matrix.imp)^2)

stop_crit <- TRUE
i <- 1
while (stop_crit) {
  xp_data.matrix.imp_l <- lm_impute(Yl_imp = xp_data.matrix.imp_l[[1]],
                                    Yl = xp_data.matrix,
                                    Yl_exp = xp_data.matrix.exp,
                                    Bl = xp_data.matrix.imp_l[[3]],
                                    Xl = xp_data.matrix.imp_l[[2]],
                                    PC = PC,
                                    q_ncol = 2,
                                    lambda = lambda)
  
  loop_compare <- sum((xp_data.matrix.imp_l[[1]] - xp_data.matrix.imp_l[[2]] %*% (xp_data.matrix.imp_l[[3]]))^2*
                        xp_data.matrix.exp)
  cat("Iter - sc-softImpute:", i, "\n")
  print(loop_compare)
  stop_crit <- (loop_compare > stop_crit_compare/10) && (i < max_iter_imp)
  
  i <- i + 1
}

xp_data.matrix.imp = xp_data.matrix.imp_l[[1]]

# Add noise to the imputed values
err <- colSums((xp_data.matrix.imp - xp_data.matrix.imp_l[[2]] %*% (xp_data.matrix.imp_l[[3]]))^2)
err <- err/(xp_data.matrix.exp.count-PC-2-1)

n <- dim(xp_data.matrix.imp)[1]
for(i in 1:length(err)){
  xp_data.matrix.imp[,i] <- xp_data.matrix.imp[,i] + rnorm(n,0,sqrt(err[i]))*(1-xp_data.matrix.exp[,i])
}

end_time_imp = Sys.time()
imp_total_seconds = as.numeric(difftime(end_time_imp , start_time_imp, units = "secs"))
imp_hours <- floor(imp_total_seconds / 3600)
imp_minutes <- floor((imp_total_seconds %% 3600) / 60)
imp_seconds <- round(imp_total_seconds %% 60)
print(sprintf("Spent Computation Time (Imputation): %02d:%02d:%02d", imp_hours, imp_minutes, imp_seconds))

####
# 3. Generate knockoffs
####

####
# 3-1. Compute the covariance matrix
####

# Assemble X
X <- xp_data.matrix.imp_l[[2]]
coln.X <- colnames(X)
coln.X[1] = "sex"
colnames(X) = coln.X

Q_index = 1:(ncol(X)-PC)
A_index = (max(Q_index)+1):(max(Q_index)+PC)

# X_hat: fixed and given
X_mat = X[,Q_index]
A_hat = X[,A_index]

B0_hat = t(xp_data.matrix.imp_l[[3]])[,Q_index]
B1_hat = t(xp_data.matrix.imp_l[[3]])[,A_index]

Gamma_hat = (solve(t(X_mat) %*% X_mat) %*% t(X_mat) %*% A_hat) %>% t()
U_hat = A_hat - (X_mat %*% t(Gamma_hat))
SigmaA_hat = (t(U_hat) %*% U_hat)/nrow(U_hat)

W_hat = xp_data.matrix.imp - (X_mat %*% t(B0_hat)) - (X_mat %*% t(Gamma_hat) %*% t(B1_hat))
SigmaG_hat = diag(err) + (B1_hat %*% SigmaA_hat %*% t(B1_hat))

W <- diag(err)

## exlcude NAs
excl_id <- which(is.na(diag(SigmaG_hat)))
if(length(excl_id) > 0){
  xp_data.matrix <- xp_data.matrix[,-excl_id]
  SigmaG_hat <- SigmaG_hat[-excl_id,-excl_id]
  feature.names.new <- feature.names.new[-excl_id]
  xp_data.sub <- subset(x = test_data, features = feature.names.new)
  xp_data.matrix.fix <- xp_data.matrix.fix[-excl_id]
  xp_data.matrix.exp <- xp_data.matrix.exp[,-excl_id]
  # W <- diag(W)[-excl_id]
}

rm(test_data)

####
# 3-2. Sample knockoffs
####

W_hat_imp = W_hat

iter_time_start_KO = Sys.time()

xp_data.knockoff = 
  create.gaussian(W_hat_imp, 
                  mu = rep(0, times = ncol(W_hat_imp)), 
                  Sigma = SigmaG_hat,
                  method = "asdp",
                  diag_s = 1.95*W)

iter_time_end_KO = Sys.time()

KO_total_seconds = as.numeric(difftime(iter_time_end_KO, iter_time_start_KO, units = "secs"))
KO_hours <- floor(KO_total_seconds / 3600)
KO_minutes <- floor((KO_total_seconds %% 3600) / 60)
KO_seconds <- round(KO_total_seconds %% 60)
print(sprintf("Spent Computation Time (KO): %02d:%02d:%02d", KO_hours, KO_minutes, KO_seconds))

# scale back the knockoff

xp_data.matrix <- sweep(xp_data.matrix,2,xp_data.matrix.fix,FUN = "+")*(xp_data.matrix.exp)

xp_data.knockoff.sav0 = xp_data.knockoff + (X_mat %*% t(B0_hat)) + (X_mat %*% t(Gamma_hat) %*% t(B1_hat))
xp_data.knockoff.sav <- t(t(xp_data.knockoff.sav0) + xp_data.matrix.fix)

# Use only the expressed part
xp_data.knockoff.sav <- xp_data.knockoff.sav*(xp_data.matrix.exp) + 0*(1 - xp_data.matrix.exp) 

rm(xp_data.matrix.imp_l)
rm(xp_data.knockoff)
rm(xp_data.knockoff.sav0)
rm(W)
rm(SigmaG_hat)
rm(Q_index, A_index)
rm(X_mat, A_hat)
rm(B0_hat, B1_hat)
rm(Gamma_hat)
rm(U_hat)
rm(SigmaA_hat)
rm(W_hat)

####
# 4. Run simulations
####

if (testUse == "wilcox_limma"){
  covariate_names = NULL
}

hyperparam_used_list$imp_total_seconds = imp_total_seconds
hyperparam_used_list$KO_total_seconds = KO_total_seconds

xp_data.matrix = xp_data.matrix %>% as(Class = "dgCMatrix")
xp_data.knockoff.sav = xp_data.knockoff.sav %>% as(Class = "dgCMatrix")

FinalResult = list()

start_time_nsim_set = NULL
end_time_nsim_set = NULL

for(nsim in 1:nSim){
  # Computation Time
  start_time_nsim = Sys.time()
  start_time_nsim_set = c(start_time_nsim_set, start_time_nsim)
  
  n <- dim(xp_data.matrix)[1]
  p <- dim(xp_data.matrix)[2]
  cat("n: ", n, "\n")
  cat("p: ", p, "\n")
  
  set.seed(seed = seed_num + nsim)
  # Generate coefficients
  
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
  
  # xp_data.matrix_expressed_prop = xp_data.matrix_denominator/n
  # filtered_idx = ((xp_data.matrix_expressed_prop > 0.1) & (xp_data.matrix_expressed_prop < 0.3)) %>% which()
  filtered_idx = 1:p
  
  sig_coef <- sample(filtered_idx, n_signal)
  norm_coef[sig_coef] <- tmp_coef
  
  
  # Introduce the logit link
  label_p <- 1/(1 + exp(- xp_data.matrix_scale %*% norm_coef)) %>% as.vector()
  rm(xp_data.matrix_scale)
  
  # Generate labels
  label_binom <- rbinom(n = length(label_p), size = 1, prob = label_p)
  label_ad <- ifelse((label_binom == 1), yes = "ad", no = "healthy")
  
  xp_data.sub <- SetIdent(object = xp_data.sub, value = label_ad)
  
  # Calculate importance statistics
  
  ## LRT, MAST, WRT, LCD ##
  if (testUse != "LCD")
  {
    # use observed values to find p-values
    # shouldn't matter which slot I change as long as I use it. # layer = slot Seurat v5 use layer instead but FindMarkers use slot for the both.
    xp_data.sub <- SetAssayData(object = xp_data.sub, layer = "data", new.data = t(xp_data.matrix))
    
    result.sub <- FindMarkers(xp_data.sub, ident.1 = 'healthy', ident.2 = 'ad', slot = "data",
                              min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use = testUse, latent.vars = covariate_names) 
    
    xp_data.sub <- SetAssayData(object = xp_data.sub, layer = "data", new.data = t(xp_data.knockoff.sav))
    
    # Find p-values for knockoffs
    result.sub.k <- FindMarkers(xp_data.sub, ident.1 = 'healthy', ident.2 = 'ad', slot = "data",
                                min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use = testUse, latent.vars = covariate_names) 
    
    # Importance statistics
    W_imp = -log(result.sub$p_val[match(feature.names.new, rownames(result.sub))]) + 
      log(result.sub.k$p_val[match(feature.names.new, rownames(result.sub.k))]) # 0<= -log(p-value)
    W_imp_b = -log(result.sub$p_val_adj[match(feature.names.new, rownames(result.sub))]) + 
      log(result.sub.k$p_val_adj[match(feature.names.new, rownames(result.sub.k))]) # 0<= -log(p-value)
    
    cat("####", testUse, "####", "\n")
    print(covariate_names)
  } else if (testUse == "LCD") {
    ## LCD ##
    set.seed(seed = seed_num + nsim) # for reproducibility of cv.glmnet
    
    xp_data.matrix_sd = 
      rbind(xp_data.matrix, xp_data.knockoff.sav) %>%
      apply(2, sd)
    
    xp_data.matrix_scale =
      xp_data.matrix %>%
      sweep(2, xp_data.matrix_sd, "/")
    
    xp_data.knockoff.sav_scale = 
      xp_data.knockoff.sav %>%
      sweep(2, xp_data.matrix_sd, "/")
    
    colnames(xp_data.knockoff.sav_scale) = str_c(colnames(xp_data.matrix_scale), "_1")
    
    covariate_mat_center = covariate_mat %>% scale()
    cat("####", testUse, "####", "\n")
    colnames(covariate_mat) %>% print()
    
    xp_data.full = cbind(covariate_mat_center,
                         xp_data.matrix_scale,
                         xp_data.knockoff.sav_scale) %>%
      as(Class = "dgCMatrix")
    
    cv_lasso_fit = 
      cv.glmnet(y = label_binom, x = xp_data.full, 
                nfolds = 10,
                alpha = 1, # alpha = 1 means lasso
                family = "binomial",
                standardize = TRUE) 
    
    #############################
    # 1. Other than LCD_crit_use
    
    LCD_crit_not_use = setdiff(x = c("lambda.1se", "lambda.min"), y = LCD_crit_use)
    cv_lasso_fit_best_lambda = cv_lasso_fit[[LCD_crit_not_use]] #lambda.1se #lambda.min
    cv_lasso_fit_coef = coef(cv_lasso_fit, s = cv_lasso_fit_best_lambda) %>% t()
    cv_lasso_fit_coef = cv_lasso_fit_coef[,colnames(cv_lasso_fit_coef) != "(Intercept)", drop = F]
    
    
    cv_lasso_fit_coef_orig = cv_lasso_fit_coef[, match(colnames(xp_data.matrix_scale), colnames(cv_lasso_fit_coef))]
    # max(abs(cv_lasso_fit_coef_orig - cv_lasso_fit_coef_orig2))
    cv_lasso_fit_coef_ko = cv_lasso_fit_coef[,match(colnames(xp_data.knockoff.sav_scale), colnames(cv_lasso_fit_coef))]
    # max(abs(cv_lasso_fit_coef_ko - cv_lasso_fit_coef_ko2))
    
    print("Nonzero original")
    print(LCD_crit_not_use)
    (cv_lasso_fit_coef_orig != 0) %>% sum() %>% print()
    
    W_imp = abs(cv_lasso_fit_coef_orig) - abs(cv_lasso_fit_coef_ko)
    
    ##############################
    # 2. LCD_crit_use  
    
    cv_lasso_fit_best_lambda_b = cv_lasso_fit[[LCD_crit_use]] #lambda.1se #lambda.min
    cv_lasso_fit_coef_b = coef(cv_lasso_fit, s = cv_lasso_fit_best_lambda_b) %>% t()
    cv_lasso_fit_coef_b = cv_lasso_fit_coef_b[,colnames(cv_lasso_fit_coef_b) != "(Intercept)", drop = F]
    
    cv_lasso_fit_coef_orig_b = cv_lasso_fit_coef_b[,match(colnames(xp_data.matrix_scale), colnames(cv_lasso_fit_coef_b))]
    # cv_lasso_fit_coef_orig_b2 = cv_lasso_fit_coef_b[,colnames(cv_lasso_fit_coef_b) %in% colnames(xp_data.matrix_scale)]
    # max(abs(cv_lasso_fit_coef_orig_b - cv_lasso_fit_coef_orig_b2))
    
    cv_lasso_fit_coef_ko_b = cv_lasso_fit_coef_b[,match(colnames(xp_data.knockoff.sav_scale), colnames(cv_lasso_fit_coef_b))]
    # cv_lasso_fit_coef_ko_b2 = cv_lasso_fit_coef_b[,colnames(cv_lasso_fit_coef_b) %in% colnames(xp_data.knockoff.sav_scale)]
    # max(abs(cv_lasso_fit_coef_ko_b - cv_lasso_fit_coef_ko_b2))
    
    print("Nonzero original")
    print(LCD_crit_use)
    (cv_lasso_fit_coef_orig_b != 0) %>% sum() %>% print()
    
    W_imp_b = abs(cv_lasso_fit_coef_orig_b) - abs(cv_lasso_fit_coef_ko_b)
    
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
  
  threshold <- knockoff.threshold(W_imp)
  
  selected <- which(W_imp >= threshold)
  
  print((length(selected) - sum(selected %in% sig_coef))/length(selected)) # FDR
  
  print(sum(selected %in% sig_coef)/length(sig_coef)) # Power
  
  threshold <- knockoff.threshold(W_imp_b)
  
  selected_b <- which(W_imp_b >= threshold)
  
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
# 5. Save the results
####

data_name = "SupParLob"
method_name = "decomp"

if ((LCD_crit_use == "lambda.1se")&(testUse == "LCD"))
{
  testUseFull = str_c(testUse, "1")
} else {
  testUseFull = testUse
}

file_name = 
  str_c(data_name,
        "_method", method_name,
        "_testUse", testUseFull,
        "_down", ifelse(is.null(down), yes = "NULL", no = down),
        "_gene", max(gene.index),
        "_sgn", sign_strength,
        "_PC", PC,
        "_10llam", 10*llam,
        "_nSim", nSim, 
        "_seed", seed_num,
        "_maxIterImp", max_iter_imp,
        ".rda")

save(FinalResult,
     file = file_name)
