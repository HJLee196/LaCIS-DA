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

ebh <- function(E, alpha){
  
  p <- length(E)
  E_ord <- order(E, decreasing = TRUE)
  E <- sort(E, decreasing = TRUE)
  comp <- E >= (p / alpha / (1:p))
  id <- suppressWarnings(max(which(comp>0)))
  if(id > 0){
    rej <- E_ord[1:id]
  }else{
    rej <- NULL
  }
  return(list(rej = rej))
}

####
# 1-2. Read the dataset
#### 

file_path = "path/to/your/rds_folder"
# file_path = "./"
cell_name = 'microglia'
data_name = 
  paste0("HuVascAD_", 
         cell_name,
         ".rds")
full_file_path = paste(file_path, data_name, sep = "/")
HuVascAD = readRDS(file = full_file_path)
options(future.globals.maxSize= +Inf)
max_gene = rownames(HuVascAD) %>% length()

####
# 1-3. Set hyperparameters
####

seed_num = 2017
set.seed(seed_num)
down = 200
gene.index = 1:max_gene
n_signal = 50
PC = 11
llam = 0.1 
nSim = 1
m_kos = 10 # the number of the sets of knockoffs
testUse = "LR" # c("wilcox_limma", "LR", "MAST", "LCD")
LCD_crit_use = "lambda.1se"
max_iter_imp = 100

hyperparam_used_list = 
  list(seed_num = seed_num,
       down = down,
       max_gene.index = max(gene.index),
       # sign_strength = sign_strength,
       n_signal = n_signal,
       PC = PC,
       llam = llam,
       nSim = nSim,
       max_iter_imp = max_iter_imp,
       testUse = testUse)

data_name = "HuVascAD"
method_name = "multidecomp"

if ((LCD_crit_use == "lambda.1se")&(testUse == "LCD"))
{
  testUseFull = str_c(testUse, "1")
} else {
  testUseFull = testUse
}

file_name =
  str_c(data_name,
        "_cell", cell_name,
        "_method", method_name,
        "_testUse", testUseFull,
        "_down", ifelse(is.null(down), yes = "NULL", no = down),
        "_gene", max(gene.index),
        # "_sgn", sign_strength,
        "_PC", PC,
        "_10llam", 10*llam,
        # "_nSim", nSim,
        "_seed", seed_num,
        "_maxIterImp", max_iter_imp,
        ".rda")

knockoffs_name =
  str_c("Knockoffs_",
        data_name,
        "_cell", cell_name,
        "_method", method_name,
        # "_testUse", testUseFull,
        "_down", ifelse(is.null(down), yes = "NULL", no = down),
        "_gene", max(gene.index),
        # "_sgn", sign_strength,
        "_PC", PC,
        "_10llam", 10*llam,
        # "_nSim", nSim,
        "_seed", seed_num,
        "_maxIterImp", max_iter_imp,
        ".rda")

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
# 3: CDR, age, age2
k_num = 3 + length(unique(HuVascAD$gender)) - 1 + length(unique(HuVascAD$batch)) - 1

####
# 1-4. Downsample the data
####

if(is.null(down)){
  test_data <- HuVascAD # use full dataset
} else {

  healthy_M <- names(HuVascAD$treat)[HuVascAD$treat %in% "Control"]
  ad_M <- names(HuVascAD$treat)[HuVascAD$treat %in% "AD"]
  
  set.seed(seed_num)
  down_ind <- c(sample(healthy_M,down),
                sample(ad_M,down))
  which(names(HuVascAD$treat) %in% down_ind)

  test_data <- subset(HuVascAD,
                      cells = which(names(HuVascAD$treat) %in% down_ind))
}
rm(HuVascAD) # rm dataset to save memory

# Take a subset of genes
feature.names <- rownames(test_data)[gene.index]

DefaultAssay(test_data) <- "RNA"
xp_data.sub <- subset(x = test_data, features = feature.names)
DefaultAssay(xp_data.sub) <- "RNA"

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

xp_data.sub$CDR = CDR

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
        xp_data.sub$age2,
        as.numeric(xp_data.sub$gender),
        samp,
        batch)
colnames(covariate_mat)[1:4] = c("CDR", "age", "age2", "gender")


covariate_pattern = impute_cov_names %>% str_c(collapse = "|") %>% str_replace(pattern = "(^.*$)", replacement = "\\^(\\1)")

selected_covariate_ind =
  colnames(covariate_mat) %>%
  str_starts(pattern = covariate_pattern)
selected_covariate = colnames(covariate_mat)[selected_covariate_ind]

covariate_mat = covariate_mat[,selected_covariate]

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

X_0 <- cbind(covariate_mat, X_0)
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
                                  q_ncol = k_num,
                                  lambda = lambda)

xp_data.matrix.imp_l <- lm_impute(Yl_imp = xp_data.matrix.imp_l[[1]],
                                  Yl = xp_data.matrix,
                                  Yl_exp = xp_data.matrix.exp,
                                  Bl = xp_data.matrix.imp_l[[3]],
                                  Xl = xp_data.matrix.imp_l[[2]],
                                  PC = PC,
                                  q_ncol = k_num,
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
                                    q_ncol = k_num,
                                    lambda = lambda)

  loop_compare <- sum((xp_data.matrix.imp_l[[1]] - xp_data.matrix.imp_l[[2]] %*% (xp_data.matrix.imp_l[[3]]))^2*
                        xp_data.matrix.exp)
  cat("Iter - sc-softImpute:", i, "\n")
  print(loop_compare)
  stop_crit <- (loop_compare > stop_crit_compare/10) && (i < max_iter_imp)

  i <- i + 1
}

xp_data.matrix.imp = xp_data.matrix.imp_l[[1]]

# Add noise to the imputed values # for LR it is not necessary
err <- colSums((xp_data.matrix.imp - xp_data.matrix.imp_l[[2]] %*% (xp_data.matrix.imp_l[[3]]))^2)
err <- err/(xp_data.matrix.exp.count-PC-k_num-1)
print("Summary err")
summary(err)

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
# 3. Generate multiple knockoffs
####

xp_data.knockoff.sav_list = list()

# Assemble X
X <- xp_data.matrix.imp_l[[2]]

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

W_hat_imp = W_hat

iter_time_start_KO = Sys.time()

# Use He et al.'s method to generate multiple Gaussian knockoffs
# The details are in the appendix
p <- dim(xp_data.matrix.imp)[2]

S_multi_decomp = diag(((m_kos+1)/m_kos-0.05)*(err))
SigmaG_hat_inv = solve(SigmaG_hat)

C_Cov = (2*S_multi_decomp) - (S_multi_decomp %*% SigmaG_hat_inv %*% S_multi_decomp)
V1_Cov = C_Cov - ((m_kos - 1)/m_kos)*S_multi_decomp
rm(C_Cov)
# V1 = MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = V1_Cov)
V1_cho = chol(V1_Cov)
V1 = matrix(rnorm(n*p), nrow =n, ncol =p)
V1 = V1 %*% V1_cho
rm(V1_Cov)
rm(V1_cho)
V2_array = array(dim = c(n, p, m_kos))
V2_sd = sqrt(diag(S_multi_decomp))
for (m_ko in 1:m_kos)
{
  # V2_array[,,m_ko] = MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = S_multi_decomp)
  V2_tmp = matrix(rnorm(p*n), nrow =p, ncol = n)
  V2_array[,,m_ko] <- t(V2_tmp * V2_sd)
}
rm(V2_sd)

V2bar = apply(X = V2_array, MARGIN = c(1,2), FUN = mean)
V2_array = sweep(V2_array, c(1, 2), V2bar, "-")
V2_array = sweep(V2_array, c(1, 2), V1, "+")
# V2_array[,1948,9] %>% head(6)
rm(V2bar)
rm(V1)

start_time_ko = Sys.time()
for (m_ko in 1:m_kos)
{
  xp_data.knockoff = (W_hat_imp %*% t(diag(1,p,p) - (S_multi_decomp %*% SigmaG_hat_inv)) ) + V2_array[,,m_ko]

  colnames(xp_data.knockoff) <- colnames(xp_data.matrix)
  rownames(xp_data.knockoff) <- rownames(xp_data.matrix)

  ## scale back the knockoff

  xp_data.knockoff.sav0 = xp_data.knockoff + (X_mat %*% t(B0_hat)) + (X_mat %*% t(Gamma_hat) %*% t(B1_hat))
  xp_data.knockoff.sav <- t(t(xp_data.knockoff.sav0) + xp_data.matrix.fix)
  xp_data.knockoff.sav <- xp_data.knockoff.sav*(xp_data.matrix.exp) + 0*(1 - xp_data.matrix.exp)

  xp_data.knockoff.sav_list[[m_ko]] = xp_data.knockoff.sav
  print(m_ko)
}

iter_time_end_KO = Sys.time()

KO_total_seconds = as.numeric(difftime(iter_time_end_KO, iter_time_start_KO, units = "secs"))
KO_hours <- floor(KO_total_seconds / 3600)
KO_minutes <- floor((KO_total_seconds %% 3600) / 60)
KO_seconds <- round(KO_total_seconds %% 60)
print(sprintf("Spent Computation Time (KO): %02d:%02d:%02d", KO_hours, KO_minutes, KO_seconds)) 

## scale back the original
xp_data.matrix <- sweep(xp_data.matrix,2,xp_data.matrix.fix,FUN = "+")*(xp_data.matrix.exp)

rm(xp_data.matrix.imp_l)
rm(xp_data.knockoff)
rm(xp_data.knockoff.sav)
rm(xp_data.knockoff.sav0)
rm(V2_array)
rm(SigmaG_hat)
rm(Q_index, A_index)
rm(X_mat, A_hat)
rm(B0_hat, B1_hat)
rm(Gamma_hat)
rm(U_hat)
rm(SigmaA_hat)
rm(W_hat)

rm(xp_data)
rm(xp_data.matrix.imp)
rm(xp_data.matrix.exp)
rm(B_0)
rm(X_0)

hyperparam_used_list$imp_total_seconds = imp_total_seconds
hyperparam_used_list$KO_total_seconds = KO_total_seconds

First_Step_List_Knockoffs =
  list(
    hyperparam_used_list = hyperparam_used_list,
    xp_data.knockoff.sav_list = xp_data.knockoff.sav_list,
    xp_data.matrix = xp_data.matrix,
    xp_data.sub = xp_data.sub,
    file_name = file_name,
    knockoffs_name = knockoffs_name,
    feature.names.new = feature.names.new,
    covariate_mat = covariate_mat,
    covariate_names = covariate_names
  )

save(First_Step_List_Knockoffs,
     file = paste(file_path, knockoffs_name, sep = "/"))

####
# 4. Run simulations
####

if ("HuVascAD" %in% ls()){
  rm(HuVascAD)
}
if ("test_data" %in% ls()){
  rm(test_data)
}

# Second Step Initiallize

load(file = paste(file_path, knockoffs_name, sep = "/"))

hyperparam_used_list = First_Step_List_Knockoffs$hyperparam_used_list
xp_data.knockoff.sav_list = First_Step_List_Knockoffs$xp_data.knockoff.sav_list %>% lapply(FUN = as, Class = "dgCMatrix")
xp_data.matrix = First_Step_List_Knockoffs$xp_data.matrix %>% as(Class = "dgCMatrix")
xp_data.sub = First_Step_List_Knockoffs$xp_data.sub

knockoffs_name = First_Step_List_Knockoffs$knockoffs_name
feature.names.new = First_Step_List_Knockoffs$feature.names.new
covariate_mat = First_Step_List_Knockoffs$covariate_mat

covariate_names = First_Step_List_Knockoffs$covariate_names
if (testUse == "wilcox_limma"){
  covariate_names = NULL
}

xp_data.matrix.exp = (xp_data.matrix != 0)
n = xp_data.matrix.exp %>% nrow()

rm(First_Step_List_Knockoffs)

FinalResult = list()

start_time_nsim_set = NULL
end_time_nsim_set = NULL
target_FDR_temp = 0.1

for(nsim in 1:nSim){
  # Computing Time
  start_time_nsim = Sys.time()
  start_time_nsim_set = c(start_time_nsim_set, start_time_nsim)
  
  n <- dim(xp_data.matrix)[1]
  p <- dim(xp_data.matrix)[2]
  cat("n: ", n, "\n")
  cat("p: ", p, "\n")
  
  set.seed(seed = seed_num + nsim)
  
  # Set labels to be observed ones
  Idents(xp_data.sub) <- 'treat'
  ## For LCD
  label_binom = 
    ifelse(xp_data.sub$treat == "Control", 
           yes = 0, no = 1)
  
  W_imp <- matrix(NA,nrow = m_kos + 1, ncol = length(feature.names.new))
  W_imp_b <- matrix(NA,nrow = m_kos + 1, ncol = length(feature.names.new))
  
  # Calculate importance statistic for each knockoffs
  
  ## LRT, MAST, WRT, "LCD" ##
  if (testUse != "LCD")
  {
    # use observed values to find p-values
    xp_data.sub <- SetAssayData(object = xp_data.sub, layer = "data", new.data = t(xp_data.matrix))
    
    result.sub <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
                              min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use = testUse, latent.vars = covariate_names)
    
    W_imp[1,] = -log(result.sub$p_val[match(feature.names.new, rownames(result.sub))])
    W_imp_b[1,] = -log(result.sub$p_val_adj[match(feature.names.new, rownames(result.sub))])
    
    for (m_ko in 1:m_kos)
    {
      print("m_knockoffs")
      print(m_ko)
      xp_data.knockoff.sav = xp_data.knockoff.sav_list[[m_ko]]
      
      xp_data.sub <- SetAssayData(object = xp_data.sub, layer = "data", new.data = t(xp_data.knockoff.sav))
      
      # Find p-values for knockoffs
      result.sub.k <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
                                  min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use = testUse, latent.vars = covariate_names)
      
      # Importance statistics
      W_imp[1+m_ko,] = -log(result.sub.k$p_val[match(feature.names.new, rownames(result.sub.k))])
      W_imp_b[1+m_ko,] = -log(result.sub.k$p_val_adj[match(feature.names.new, rownames(result.sub.k))])
    }
    
    cat("####", testUse, "####", "\n")
    print(covariate_names)
  } else if (testUse == "LCD") {
    ## LCD ##
    set.seed(seed = seed_num + nsim) # for reproducibility of cv.glmnet
    
    xp_data.matrix_sq_sum = (xp_data.matrix^2) %>% colSums()
    xp_data.knockoff.sav_sq_sum = 
      xp_data.knockoff.sav_list %>% 
      lapply(FUN = function(x){
        apply(X = x^2, MARGIN = 2, FUN = sum)
      }) %>%
      do.call(what = rbind) %>%
      apply(2, sum)
    
    xp_data.matrix_sd = ((xp_data.matrix_sq_sum + xp_data.knockoff.sav_sq_sum)/(n*(m_kos + 1) - 1)) %>% sqrt()
  
    xp_data.full = 
      xp_data.matrix %>%
      sweep(MARGIN = 2, STATS = xp_data.matrix_sd, FUN = "/")
   
    for (m_ko in 1:m_kos)
    {
      xp_data.knockoff.sav_scale = 
        xp_data.knockoff.sav_list[[m_ko]] %>%
        sweep(MARGIN = 2, STATS = xp_data.matrix_sd, FUN = "/")
      # xp_data.knockoff.sav_scale %>% apply(2, sd) %>% head()
      colnames(xp_data.knockoff.sav_scale) = str_c(colnames(xp_data.matrix), "_", m_ko)
      xp_data.full = cbind(xp_data.full, xp_data.knockoff.sav_scale)
    }
    
    covariate_mat_center = covariate_mat %>% scale()
    cat("####", testUse, "####", "\n")
    colnames(covariate_mat) %>% print()
    
    xp_data.full = 
      cbind(covariate_mat_center, xp_data.full)
    
    cv_lasso_fit =
      cv.glmnet(y = label_binom,
                x = xp_data.full,
                nfolds = 5,
                alpha = 1, # alpha = 1 means lasso
                family = "binomial")
   
    #############################
    # 1. Other than LCD_crit_use
    
    LCD_crit_not_use = setdiff(x = c("lambda.1se", "lambda.min"), y = LCD_crit_use)
    cv_lasso_fit_best_lambda = cv_lasso_fit[[LCD_crit_not_use]] #lambda.1se #lambda.min
    cv_lasso_fit_coef = coef(cv_lasso_fit, s = cv_lasso_fit_best_lambda) %>% t()
    cv_lasso_fit_coef = cv_lasso_fit_coef[,colnames(cv_lasso_fit_coef) != "(Intercept)", drop = F]
    
    cv_lasso_fit_coef_orig = cv_lasso_fit_coef[,colnames(cv_lasso_fit_coef) %in% colnames(xp_data.matrix)]
    
    print("Nonzero original")
    print(LCD_crit_not_use)
    (cv_lasso_fit_coef_orig != 0) %>% sum() %>% print()
    
    W_imp[1, ] = cv_lasso_fit_coef_orig %>% abs()
    for (m_ko in 1:m_kos)
    {
      colnames_mth_ko = str_c(colnames(xp_data.matrix), "_", m_ko)
      
      cv_lasso_fit_coef_mth_ko = cv_lasso_fit_coef[,colnames(cv_lasso_fit_coef) %in% colnames_mth_ko]
      
      W_imp[1 + m_ko, ] = cv_lasso_fit_coef_mth_ko %>% abs()
    }
    
    ##############################
    # 2. LCD_crit_use  
    cv_lasso_fit_best_lambda_b = cv_lasso_fit[[LCD_crit_use]] #lambda.1se #lambda.min
    cv_lasso_fit_coef_b = coef(cv_lasso_fit, s = cv_lasso_fit_best_lambda_b) %>% t()
    cv_lasso_fit_coef_b = cv_lasso_fit_coef_b[,colnames(cv_lasso_fit_coef_b) != "(Intercept)", drop = F]

    cv_lasso_fit_coef_orig_b = cv_lasso_fit_coef_b[,colnames(cv_lasso_fit_coef_b) %in% colnames(xp_data.matrix)]
    print("Nonzero original_b")
    print(LCD_crit_use)
    (cv_lasso_fit_coef_orig_b != 0) %>% sum() %>% print()

    W_imp_b[1, ] = cv_lasso_fit_coef_orig_b %>% abs()
    for (m_ko in 1:m_kos)
    {
      colnames_mth_ko = str_c(colnames(xp_data.matrix), "_", m_ko)
      cv_lasso_fit_coef_mth_ko_b = cv_lasso_fit_coef_b[,colnames(cv_lasso_fit_coef_b) %in% colnames_mth_ko]
      W_imp_b[1 + m_ko, ] = cv_lasso_fit_coef_mth_ko_b %>% abs()
    }
    
    # W_imp_b = W_imp
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
  
  W_max_ind = apply(W_imp,2,which.max)
  W_max = sapply(1:length(W_max_ind), function(x) W_imp[W_max_ind[x],x])
  W_med = sapply(1:length(W_max_ind), function(x) median(W_imp[-W_max_ind[x],x]))
  
  tau_stab <- (W_max - W_med)
  kap = (W_max_ind==1 & tau_stab>0)
  
  ts = sort(c(0, abs(tau_stab)))
  ratio = sapply(ts, function(t) (1/m_kos + 1/m_kos*sum((tau_stab >= t)*(1-kap)))/max(1, 
                                                                                      sum((tau_stab >= t)*(kap))))
  ok = which(ratio <= target_FDR_temp)
  threshold <- ifelse(length(ok) > 0, ts[ok[1]], Inf)
  
  selected <- which(tau_stab*kap >= threshold)

  print("Selected Genes (Not Selected Setting):")
  print(feature.names.new[selected])
  cat("Number of selected genes:", length(selected), "\n")
  
  # bonferroni 
  W_max_ind = apply(W_imp_b,2,which.max)
  W_max = sapply(1:length(W_max_ind), function(x) W_imp_b[W_max_ind[x],x])
  W_med = sapply(1:length(W_max_ind), function(x) median(W_imp_b[-W_max_ind[x],x]))
  
  tau_stab <- (W_max - W_med)
  kap = (W_max_ind==1 & tau_stab>0)
  
  ts = sort(c(0, abs(tau_stab)))
  ratio = sapply(ts, function(t) (1/m_kos + 1/m_kos*sum((tau_stab >= t)*(1-kap)))/max(1, 
                                                                                      sum((tau_stab >= t)*(kap))))
  ok = which(ratio <= target_FDR_temp)
  threshold <- ifelse(length(ok) > 0, ts[ok[1]], Inf)
  
  selected_b <- which(tau_stab*kap >= threshold)
  
  print("Selected Genes (Selected Setting):")
  print(feature.names.new[selected_b])
  cat("Number of selected genes:", length(selected_b), "\n")
  
  # 2. ebh # (1. Multi 2. ebh) 
  
  E_temp = matrix(NA, m_kos, dim(W_imp_b)[2])
  for(m_ko in 1:m_kos){
    
    E_temp[m_ko,] = W_imp_b[1,] - W_imp_b[(m_ko+1),]
    tau = knockoff.threshold(E_temp[m_ko,],
                             fdr = target_FDR_temp/2)
    E_temp[m_ko,] = (E_temp[m_ko,] >= tau) / (1 + sum(E_temp[m_ko,] <= -tau))
  }
  
  E_temp <- (dim(W_imp_b)[2]) * colMeans(E_temp)
  
  selected_b_ebh = ebh(E_temp, alpha = target_FDR_temp)$rej
  print("Selected Genes (ebh):")
  print(feature.names.new[selected_b_ebh])
  cat("Number of selected genes:", length(selected_b_ebh), "\n")
  
  # Computation Time
  end_time_nsim = Sys.time()
  end_time_nsim_set = c(end_time_nsim_set, end_time_nsim)
  
  nsim_total_seconds = as.numeric(difftime(end_time_nsim , start_time_nsim, units = "secs"))
  nsim_minutes <- floor(nsim_total_seconds / 60)
  nsim_seconds <- round(nsim_total_seconds %% 60)
  
  print(sprintf("Spent Computation Time (each sim): %02d:%02d", nsim_minutes, nsim_seconds))
  
  FinalResult[[nsim]] <- (list(
    W_imp = W_imp, 
    selected = selected, 
    W_imp_b = W_imp_b, 
    selected_b = selected_b, 
    selected_b_ebh = selected_b_ebh,
    feature.names.new = feature.names.new,
    hyperparam_used_list = hyperparam_used_list,
    nsim_total_seconds = nsim_total_seconds,
    target_FDR = target_FDR_temp))
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

save(FinalResult,
     file = file_name)

