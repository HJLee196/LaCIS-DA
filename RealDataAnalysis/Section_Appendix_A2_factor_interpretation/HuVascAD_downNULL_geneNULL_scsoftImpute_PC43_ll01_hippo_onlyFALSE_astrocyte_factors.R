# lmimpute Imputation
# Use lm to generate coefficients, and subtract the mean from the dataset.
# Fit lm with intercept.
# Include factor variables (just sex in this case).
# Only use expressed genes when calculating the variance in W.

# Include factor variables when calculating covariance.

# 'data' when generating labels.

# Repeating the knockoff experiment in order to check the FDR

rm(list=ls())
# down=500
# gene.index = 1:2000
# method="random"
# sign_strength=3
# n_signal = 50
# expresc2 = T
# PC = 10
# llam = 0.1
# nSim = 1
# W_IT = 5
#################################################################################
# Loading Other necessary R Libraries:
#################################################################################
#if (!requireNamespace("tidyverse", quietly=TRUE)) install.packages("tidyverse", repos="https://cloud.r-project.org")
#if (!requireNamespace("knockoff", quietly=TRUE)) install.packages("knockoff")
#if (!requireNamespace("softImpute", quietly=TRUE)) install.packages("softImpute")

suppressPackageStartupMessages({
  library(Seurat)
  library(knockoff)
  library(MAST)
  library(Matrix)
  library(softImpute)
  library(tidyverse)
})

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

#setwd("/Users/linxiliu/Dropbox/Research/Knockoff_scRNAseq/Data/seperate_data/")

seed_num = 17
set.seed(seed_num)
down = NULL
gene.index = NULL
m_kos = 10 # the number of the sets of knockoffs
# method="random" # Should be removed
llam = 0.1 # max_eigen X llam (0.1 in the original one)
nSim = 1
celltype = "astrocyte"
PC = 43
#"wilcox"
#"LR"
#"MAST"
#"LCD"
sign_strength=3
n_signal = 50
hyperparam_used_list = 
  list(seed_num = seed_num,
       down = down,
       max_gene.index = max(gene.index),
       sign_strength = sign_strength,
       n_signal = n_signal,
       PC = PC,
       llam = llam,
       nSim = nSim,
       cellType = celltype)

#################################################################################
#  Load External Data
#################################################################################
# This seems to be the only data set I can upload. The others are too large and I don't have enough memory.
# WARNING: The wrapper function 'Seurat_repeat' is written for this specific dataset. Changes will be required if 
# using a different dataset.
#HuVascAD <- readRDS(file = "./HuVascAD_neuron.rds")
datafilename <- paste0("HuVascAD_", celltype, ".rds")
HuVascAD <- readRDS(file = datafilename)
options(future.globals.maxSize= +Inf)
# down: Downsample size, should be less than the smallest group, or NULL to use the full sample.
# gene.index: Subsample gene wise.
# method: Generate labels in a random or deterministic way.
# sign_strength: signal strength used as cutoff for label generation
# n_signal: Number of genes considered as true coefficients.
# PC: Number of principal components that will be used for imputation and covariance estimation
# nSim: Number of times the signals will be regenerated.
# 
# Data that we work on will be "data"

if(is.null(down)){
  test_data <- HuVascAD # use full dataset
} else {
  # Manually downsample the data:
  # healthy_M <- names(HuVascAD@active.ident)[HuVascAD@active.ident %in% "no"] 
  # ad_M <- names(HuVascAD@active.ident)[HuVascAD@active.ident %in% "yes"] 
  # 
  # # the function subset() might have a set.seed function somewhere, 
  # # returning non-random results after the first run.
  # set.seed(NULL)
  # down_ind <- c(sample(healthy_M,down),
  #               sample(ad_M,down))
  # which(names(HuVascAD@active.ident) %in% down_ind)
  # 
  # test_data <- subset(HuVascAD, 
  #                     cells = which(names(HuVascAD@active.ident) %in% down_ind))
  
  healthy_M <- names(HuVascAD$treat)[HuVascAD$treat %in% "Control"] 
  ad_M <- names(HuVascAD$treat)[HuVascAD$treat %in% "AD"] 
  
  # the function subset() might have a set.seed function somewhere, 
  # returning non-random results after the first run.
  set.seed(NULL)
  down_ind <- c(sample(healthy_M,down),
                sample(ad_M,down))
  which(names(HuVascAD$treat) %in% down_ind)
  
  test_data <- subset(HuVascAD, 
                      cells = which(names(HuVascAD$treat) %in% down_ind))
}
rm(HuVascAD) # rm dataset to save memory

# table(test_data$treat, test_data$sample) # this is the indicator?

# Take a subset of genes
if(length(gene.index) > 0){
  feature.names <- rownames(test_data)[gene.index]
} else {
  feature.names <- rownames(test_data)
}

xp_data.sub <- subset(x = test_data, features = feature.names)

DefaultAssay(xp_data.sub) <- "RNA"
# DefaultAssay(xp_data.sub) <- "SCT"

#### Calculate CDR ####
xp_data <- GetAssayData(object = xp_data.sub, slot = "count")
# rotate to match downstream analyses.
xp_data.matrix.exp <- t(as.matrix(xp_data) > 0) # indicator matrix for expressed genes, 
xp_data.matrix.exp.count <- apply(xp_data.matrix.exp,2,sum)

# xp_data.matrix <- (as.matrix(xp_data)) # do not transpose here. Keep the genes as rows.


# Prespecify how many variables are fixed during the imputation
k_num = 7

# exclude empty variables
# feature.exp.ind <- which(apply(xp_data.matrix,1,max) > (apply(xp_data.matrix,1,min)))
#n = nrow(xp_data.matrix.exp)
#feature.exp.ind <- which(xp_data.matrix.exp.count > n*0.02) 
feature.exp.ind <- which(xp_data.matrix.exp.count > (PC+k_num+1)) 
feature.names.new <- feature.names[feature.exp.ind]  # length(feature.names.new)

xp_data.matrix.exp <- xp_data.matrix.exp[,feature.exp.ind]
xp_data.matrix.exp.count <- xp_data.matrix.exp.count[feature.exp.ind]
CDR = rowMeans(xp_data.matrix.exp) # expressed genes in each CELL


# refresh subset
xp_data.sub <- subset(x = test_data, features = feature.names.new)
# rm(test_data) # release memory space

xp_data <- GetAssayData(object = xp_data.sub, slot = "data") # Normalized data matrix


xp_data.matrix <- (as.matrix(xp_data))
# xp_data.matrix <- ScaleData(xp_data.matrix)
xp_data.matrix <- t(xp_data.matrix)

print(dim(xp_data.matrix))
print(dim(xp_data.matrix.exp))

#### Generate knockoff ####

# Center target matrix, but only center the expressed parts.
xp_data.matrix.fix <- (colSums(xp_data.matrix, na.rm = T))/xp_data.matrix.exp.count
xp_data.matrix <- sweep(xp_data.matrix,2,xp_data.matrix.fix,FUN = "-")*(xp_data.matrix.exp)

# Impute missing values

# excl_id <- which(apply(xp_data.matrix.exp, 2, sum) < (PC+2))
# if(length(excl_id) > 0){
#   xp_data.matrix <- xp_data.matrix[,-excl_id]
#   feature.names.new <- feature.names.new[-excl_id]
#   xp_data.sub <- subset(x = test_data, features = feature.names.new)
#   xp_data.matrix.fix <- xp_data.matrix.fix[-excl_id]
#   xp_data.matrix.exp <- xp_data.matrix.exp[,-excl_id]
# }
# 
# excl_id <- which(apply(xp_data.matrix.exp, 1, sum) < (PC+2))
# if(length(excl_id) > 0){
#   xp_data.matrix <- xp_data.matrix[-excl_id,]
#   xp_data.matrix.exp <- xp_data.matrix.exp[-excl_id,]
# }
start.time = Sys.time()
xp_data.matrix.imp = xp_data.matrix

tmp = RunPCA(t(xp_data.matrix.imp), slot = "data", npcs = PC)

gender = as.numeric(xp_data.sub$gender)
age = xp_data.sub$age
age2 = xp_data.sub$age^2
#filler
samp = model.matrix(y~x-1,data.frame(x = xp_data.sub$sample, y = 1))
samp = samp[,-1]
#filler
batch = model.matrix(y~x-1,data.frame(x = as.factor(xp_data.sub$batch), y = 1))
batch = batch[,-1]

xp_data.sub$age2 = xp_data.sub$age^2
xp_data.sub$CDR = CDR

X_0 = tmp@cell.embeddings
# X_0 <- cbind(gender,age,age2,samp,CDR,
#              X_0) # Error in solve.default(t(X_0) %*% X_0)
# X_0 <- cbind(gender,age,samp,CDR,
#              X_0) # Error in solve.default(t(X_0) %*% X_0)
# X_0 <- cbind(gender,samp,CDR,
#              X_0)
# X_0 <- cbind(gender,samp,
#              X_0)
# X_0 <- cbind(gender,CDR,
#              X_0)
X_0 <- cbind(gender,age, age2, batch,CDR,
             X_0)
# X_0 <- cbind(age,samp,CDR,
#              X_0) # Error in solve.default(t(X_0) %*% X_0)
# X_0 <- cbind(samp,CDR,
#              X_0)
dim(X_0)
X_0 <- sweep(X_0,2,colMeans(X_0),FUN = "-")
B_0 = solve(t(X_0)%*%X_0)%*%t(X_0)%*%xp_data.matrix.imp
# lambda0(xp_data.matrix.imp) #221.6794 gene3000
#xp_data.matrix  - X_0 %*% t(B_0)
lambda =  lambda0(xp_data.matrix.imp)*llam
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

# stop_crit_compare <- (sum((xp_data.matrix.imp_l[[1]] - xp_data.matrix.imp_l[[2]] %*% t(xp_data.matrix.imp_l[[3]]))^2*
#                     xp_data.matrix.exp))
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
  print(loop_compare)
  stop_crit <- (loop_compare > stop_crit_compare/10) && (i < 100)
  
  # loop_compare <- sum((xp_data.matrix.imp_l[[1]] - xp_data.matrix.imp)^2)
  # print(loop_compare)
  # stop_crit <- (loop_compare > (stop_crit_compare/10)) && (i < 10)
  i <- i + 1
}

xp_data.matrix.imp = xp_data.matrix.imp_l[[1]]

# Add noise to the imputed values
err <- colSums((xp_data.matrix.imp - xp_data.matrix.imp_l[[2]] %*% (xp_data.matrix.imp_l[[3]]))^2)
err <- err/(xp_data.matrix.exp.count-PC-k_num-1)
summary(err)

n <- dim(xp_data.matrix.imp)[1]
for(i in 1:length(err)){
  xp_data.matrix.imp[,i] <-  xp_data.matrix.imp[,i] + rnorm(n,0,sqrt(err[i]))*(1-xp_data.matrix.exp[,i])
}

end.time = Sys.time()

#################################################################################
#  Save the result to a folder on the cluster:
#################################################################################

FinalResult = list(xp_data.matrix.imp_l[[2]], t(xp_data.matrix.imp_l[[3]]), err)
names(FinalResult)<- c("A", "B", "D")

outputname = paste0("HuVascAD_downNULL_geneNULL_scsoftImpute_PC", PC, "_hippo_onlyFALSE_", celltype, "_factors.rda")
save(FinalResult,file = outputname)

