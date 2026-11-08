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
down = NULL
gene.index = 1:2000
sign_strength = 3
n_signal = 50
PC = 30
llam = 0.1 # max_singular X llam (0.1 in the original one)
nSim = 1
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
# rotate to match downstream analyses.
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

tmp = RunPCA(t(xp_data.matrix.imp), slot = "data", npcs = PC, verbose = F) # slot = "data", ???
X_0 = tmp@cell.embeddings
rm(tmp)

X_0 <- cbind((xp_data.sub$sex),CDR,X_0)
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
  xp_data.matrix.imp[,i] <-  xp_data.matrix.imp[,i] + rnorm(n,0,sqrt(err[i]))*(1-xp_data.matrix.exp[,i])
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
# 3-1. Generate knockoffs with imputation
####

####
# 3-1-1. Compute the covariance matrix
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

## scale back the knockoff

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
# 3-2. Generate knockoffs without imputation
####

# Center target matrix, but only center the expressed parts.
xp_data.matrix.fix <- (colSums(xp_data.matrix, na.rm = T))/xp_data.matrix.exp.count
xp_data.matrix <- sweep(xp_data.matrix,2,xp_data.matrix.fix,FUN = "-")*(xp_data.matrix.exp)

xp_data.matrix.imp = xp_data.matrix

tmp = RunPCA(t(xp_data.matrix.imp), slot = "data",npcs = PC, verbose = F) # slot = "data", ???
X_0 = tmp@cell.embeddings
rm(tmp)

X_0 <- cbind((xp_data.sub$sex),CDR,X_0)
X_0 <- sweep(X_0,2,colMeans(X_0),FUN = "-")
B_0 = solve(t(X_0)%*%X_0)%*%t(X_0)%*%xp_data.matrix.imp

start_time_imp = Sys.time()

print(paste0("lambda: ",round(lambda/llam,2)))

xp_data.matrix.imp_l <- lm_impute(Yl_imp = xp_data.matrix.imp,
                                  Yl = xp_data.matrix,
                                  Yl_exp = xp_data.matrix.exp,
                                  Bl = B_0,
                                  Xl = X_0,
                                  PC = PC,
                                  q_ncol = 2,
                                  lambda = lambda)

xp_data.matrix.imp = xp_data.matrix.imp_l[[1]]

# Add noise to the imputed values
err <- colSums((xp_data.matrix.imp - xp_data.matrix.imp_l[[2]] %*% (xp_data.matrix.imp_l[[3]]))^2)
err <- err/(xp_data.matrix.exp.count-PC-2-1)

n <- dim(xp_data.matrix.imp)[1]
for(i in 1:length(err)){
  xp_data.matrix.imp[,i] <-  xp_data.matrix.imp[,i] + rnorm(n,0,sqrt(err[i]))*(1-xp_data.matrix.exp[,i])
}

end_time_imp = Sys.time()
imp_total_seconds = as.numeric(difftime(end_time_imp , start_time_imp, units = "secs"))
imp_hours <- floor(imp_total_seconds / 3600)
imp_minutes <- floor((imp_total_seconds %% 3600) / 60)
imp_seconds <- round(imp_total_seconds %% 60)
print(sprintf("Spent Computation Time (Imputation): %02d:%02d:%02d", imp_hours, imp_minutes, imp_seconds))

####
# 3-2-1. Compute the covariance matrix
####

# Assemble X
X <- xp_data.matrix.imp_l[[2]]
coln.X <- colnames(X)
coln.X[1] = "sex"
colnames(X) = coln.X

colnames(X)

sigma.X <- cov(X)
B <- t(xp_data.matrix.imp_l[[3]])
W <- diag(err)
sigma.Y <- W + (B)%*%sigma.X%*%t(B)

## exlcude NAs
excl_id <- which(is.na(diag(sigma.Y)))
if(length(excl_id) > 0){
  xp_data.matrix <- xp_data.matrix[,-excl_id]
  sigma.Y <- sigma.Y[-excl_id,-excl_id]
  feature.names.new <- feature.names.new[-excl_id]
  xp_data.sub <- subset(x = test_data, features = feature.names.new)
  xp_data.matrix.fix <- xp_data.matrix.fix[-excl_id]
  xp_data.matrix.exp <- xp_data.matrix.exp[,-excl_id]
  # W <- diag(W)[-excl_id]
}

iter_time_start_KO = Sys.time()

xp_data.knockoff.noimpuate = 
  create.gaussian(xp_data.matrix, 
                  mu = rep(0, times = ncol(xp_data.matrix)), 
                  Sigma = sigma.Y,
                  method = "asdp",
                  diag_s = 1.95*W)

iter_time_end_KO = Sys.time()

KO_total_seconds = as.numeric(difftime(iter_time_end_KO, iter_time_start_KO, units = "secs"))
KO_hours <- floor(KO_total_seconds / 3600)
KO_minutes <- floor((KO_total_seconds %% 3600) / 60)
KO_seconds <- round(KO_total_seconds %% 60)
print(sprintf("Spent Computation Time (KO): %02d:%02d:%02d", KO_hours, KO_minutes, KO_seconds)) 

## scale back the knockoff

xp_data.matrix <- sweep(xp_data.matrix,2,xp_data.matrix.fix,FUN = "+")*(xp_data.matrix.exp)

xp_data.knockoff.sav.noimpute <- t(t(xp_data.knockoff.noimpuate) + xp_data.matrix.fix)
# Use only the expressed part
xp_data.knockoff.sav.noimpute <- xp_data.knockoff.sav.noimpute*(xp_data.matrix.exp) + 0*(1 - xp_data.matrix.exp) # why the later part? just for the clarity?

rm(xp_data.matrix.imp_l)
rm(W)

##############################################
## 4. Generate the plot for B.5 in Appendix ##
##############################################

orig_cov <- cov(xp_data.matrix)
ko_cov <- cov(xp_data.knockoff.sav)
ko_cov_noimpute <- cov(xp_data.knockoff.sav.noimpute)

p = dim(orig_cov)[1]

cov_df <- data.frame(orig_var = c(diag(orig_cov), diag(orig_cov)), 
                     ko_var = c(diag(ko_cov), diag(ko_cov_noimpute)), 
                     impute = c(rep(FALSE, p), rep(TRUE, p)))

p_A <- ggplot(cov_df, aes(x = orig_var, y = ko_var, colour = impute)) +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "grey80",
    linetype = "dashed",
    linewidth = 0.8
  )+
  geom_point(size = 1, alpha= 0.6) +
  theme_bw() +
  labs(
    x = "Covariance of original variables",
    y = "Covariance of knockoff variables "
  ) +
  scale_color_brewer(
    name = "Imputation method",
    labels = c("With imputation", "No imputation"),
    palette = "Set1"
  )+
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12),
    plot.title = element_text(size = 16, hjust = 0.5)
  )

p_A 

# Generate random indices of covariance matrix
M <- p * (p - 1) / 2
ind <- sample.int(M, size = 5000, replace = FALSE)

k_to_ij <- function(ind, p) {
  i <- ceiling((2*p - 1 - sqrt((2*p - 1)^2 - 8*ind)) / 2)
  prev <- (i - 1) * (2*p - i) / 2
  j <- ind - prev + i
  cbind(i, j)
}

idx_sample <- k_to_ij(ind, p)

orig_off_cov <- rep(0, 5000)
for(i in 1:5000) {
  orig_off_cov[i] <- orig_cov[ idx_sample[i,1], idx_sample[i, 2]]
}

ko_off_cov <- rep(0, 5000)
for(i in 1:5000) {
  ko_off_cov[i] <- ko_cov[ idx_sample[i,1], idx_sample[i, 2]]
}

ko_off_cov_noimpute <- rep(0, 5000)
for(i in 1:5000) {
  ko_off_cov_noimpute[i] <- ko_cov_noimpute[ idx_sample[i,1], idx_sample[i, 2]]
}

cov_df <- data.frame(orig_off_cov = c(orig_off_cov, orig_off_cov),
                     ko_off_cov = c(ko_off_cov, ko_off_cov_noimpute),
                     impute = c(rep(FALSE, 5000), rep(TRUE, 5000)))

p_B<- ggplot(cov_df, aes(x = orig_off_cov, y = ko_off_cov, colour = impute)) +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "grey80",
    linetype = "dashed",
    linewidth = 0.8
  )+
  geom_point(size = 1, alpha= 0.6) +
  theme_bw() +
  xlim(-0.1, 0.2)+
  ylim(-0.1, 0.2)+
  labs(
    x = "Covariance of original variables",
    y = "Covariance of knockoff variables "
  ) +
  scale_color_brewer(
    name = "Imputation method",
    labels = c("With imputation", "No imputation"),
    palette = "Set1"
  )+
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12),
    plot.title = element_text(size = 16, hjust = 0.5)
  )

p_B

library(patchwork)

(p_A | p_B ) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 12),)

ggsave(
  "scatter_plot_imputation.pdf",
  (p_A | p_B ) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "right",
          legend.title = element_text(size = 12),
          legend.text  = element_text(size = 12),),
  width = 10,
  height = 4.6
)
